"""Check the six paper products and any available detailed working tables."""

import json
from pathlib import Path
import sys

import pandas as pd
import numpy as np
from scipy.special import logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from settings import (FIGURE_ROOT, IDENTIFIER_DTYPES, ONE_ROOT, OUTPUT_ROOT,
                      FINAL_CANDIDATES_PATH, ONE_ARCSECOND_SCORES_PATH,
                      IMAGING_SCORES_PATH, IMAGING_CUTFLOW_PATH,
                      THRESHOLD_CALIBRATION_PATH, PM_SIGNIFICANCE_MAX)
from score_selection import build_secondary_pm_flags, has_two_distinct_gaia_sources

calibration = json.loads(THRESHOLD_CALIBRATION_PATH.read_text())
imaging = pd.read_csv(IMAGING_SCORES_PATH, dtype=IDENTIFIER_DTYPES, low_memory=False)
imaging_cutflow = pd.read_csv(IMAGING_CUTFLOW_PATH)
joint = pd.read_csv(FINAL_CANDIDATES_PATH, dtype=IDENTIFIER_DTYPES)
one_arcsecond = pd.read_csv(ONE_ARCSECOND_SCORES_PATH, dtype=IDENTIFIER_DTYPES)
expected_products = {
    'final_candidates.csv', 'one_arcsecond_scores.csv', 'imaging_scores.csv',
    'one_arcsecond_cutflow.csv', 'imaging_cutflow.csv', 'threshold_calibration.json',
}
product_entries = [path for path in OUTPUT_ROOT.iterdir() if path.name != '.DS_Store']
pm_flags_consistent = True
for systems, primary_id, secondary_id, significance in [
    (one_arcsecond, 'primary_component_gaia_source_id', 'secondary_component_gaia_source_id', 'secondary_pm_significance'),
    (imaging, 'primary_gaia_source_id', 'secondary_gaia_source_id', 'secondary_pm_sig')]:
    two_gaia = has_two_distinct_gaia_sources(systems[primary_id], systems[secondary_id])
    expected = build_secondary_pm_flags(systems[significance],
                                         two_gaia, PM_SIGNIFICANCE_MAX)
    pm_flags_consistent &= expected.has_two_distinct_gaia_sources.equals(systems.has_two_distinct_gaia_sources)
    pm_flags_consistent &= expected.secondary_pm_pass.equals(systems.secondary_pm_pass)
    pm_flags_consistent &= expected.secondary_pm_missing.equals(systems.secondary_pm_missing)
checks = {
    'parent_systems': len(one_arcsecond) == one_arcsecond.qso_group_id.nunique() == 712,
    'scoring_routes': one_arcsecond.scoring_route.value_counts().to_dict() == {
        'resolved_original': 324, 'resolved_supplemental': 52,
        'unresolved': 321, 'currently_unscoreable': 15},
    'threshold': calibration['model'] == 'either' and calibration['adopted_log_r_threshold'] == -6.81,
    'calibration_settings': (calibration['fit_systems'] == 109
                             and calibration['pm_significance_upper_limit'] == 70),
    'only_six_products': ({path.name for path in product_entries} == expected_products
                          and all(path.is_file() for path in product_entries)),
    'imaging_score_and_pm': (imaging.joint_score_pass.sum(), imaging.joint_secondary_pm_pass.sum(),
                            imaging.joint_pre_morphology_pass.sum()) == (406, 340, 340),
    'dcr_once_one_gaia_only': (imaging.dcr_cut_applied.equals(imaging.candidate_tier.eq('one_gaia'))
                              and imaging.loc[imaging.candidate_tier.eq('one_gaia'), 'pair_dcr_delta_deg'].gt(20).all()
                              and imaging_cutflow['Selection stage'].str.contains('DCR').sum() == 1
                              and 'joint_dcr_pass' not in imaging),
    'two_gaia_without_dcr_requirement': int((imaging.joint_final_pass
                                             & imaging.candidate_tier.eq('two_gaia')
                                             & imaging.pair_dcr_delta_deg.le(20)).sum()) == 28,
    'one_gaia_selection_unchanged': int((imaging.joint_final_pass
                                        & imaging.candidate_tier.eq('one_gaia')).sum()) == 233,
    'secondary_pm_flags': bool(pm_flags_consistent),
    'selected_secondary_pm': bool(joint.secondary_pm_pass.all()),
    'missing_secondary_pm_retained': bool(joint.secondary_pm_missing.any()),
    'redshift_error_metadata': (joint.adopted_redshift_error_used.notna().all()
                                and joint.loc[joint.redshift_catalogue.str.contains('quaia', case=False),
                                              'adopted_redshift_error_used'].gt(0).all()),
    'joint_counts': (len(joint), joint.candidate_id.nunique(), joint.selected_by_resolved.sum(),
                     joint.selected_by_unresolved.sum(), joint.selected_by_imaging.sum(),
                     joint.selected_by_both_searches.sum()) == (443, 443, 100, 44, 334, 35),
}
skipped_checks = []
fit_path = ONE_ROOT / 'threshold_fit_candidates.csv'
if fit_path.is_file():
    fit = pd.read_csv(fit_path)
    checks['calibration_sample_before_secondary_pm'] = bool(
        len(fit) == 109 and fit.secondary_pm_significance.lt(70).all()
        and fit.secondary_pm_significance.ge(3).any())
else:
    skipped_checks.append('Detailed calibration sample: requires the working fit table.')
blend_path = ONE_ROOT / 'blend_scores_all_nonresolved_by_model_and_band_set.csv'
if blend_path.is_file():
    blends = pd.read_csv(blend_path)
    uncertain_blends = blends.loc[blends.redshift_error.gt(0)]
    grids_match = True
    for row in uncertain_blends.itertuples():
        grid = pd.DataFrame(json.loads(row.redshift_grid_json))
        finite = np.isfinite(grid.log_like_qq) & np.isfinite(grid.log_like_qs)
        expected = [logsumexp(grid.loc[finite, field] + np.log(grid.loc[finite, 'weight']))
                    for field in ('log_like_qq', 'log_like_qs')]
        grids_match &= (len(grid) == 9 and np.isclose(grid.weight.sum(), 1.0)
                        and np.allclose(expected, [row.log_like_qq, row.log_like_qs], atol=1e-10, rtol=0))
    checks['blend_grid_likelihoods'] = bool(grids_match) and len(uncertain_blends) == 516
else:
    skipped_checks.append('Detailed blend likelihood grids: requires working per-model scores.')

figures = [
    'one_arcsecond/resolved_pm_calibration.pdf',
    'one_arcsecond/resolved_colour_redshift.pdf',
    'imaging/colour_redshift.pdf',
    'imaging/pre_morphology_separation.pdf',
    'joint/galactic_latitude.pdf',
]
for name in figures:
    path = FIGURE_ROOT / name
    checks[name] = path.is_file() and path.stat().st_size > 1000
summary = {'passed': all(checks.values()), 'checks': {name: bool(value) for name, value in checks.items()}}
summary['skipped_checks'] = skipped_checks
print(json.dumps(summary, indent=2))
raise SystemExit(0 if summary['passed'] else 1)
