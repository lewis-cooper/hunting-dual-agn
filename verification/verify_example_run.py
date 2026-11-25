"""Check the six paper products and any available detailed working tables."""

import json
import hashlib
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
from settings import (IMAGING_SEARCH_ROOT, LITERATURE_ROOT, LITERATURE_WORK_ROOT,
                      ONE_ARCSECOND_PARENT_PATH, ONE_ARCSECOND_SEARCH_ROOT,
                      ONE_ARCSECOND_CUTFLOW_PATH)
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
# Invariants apply to the automatic release, including fresh native inputs.
# Historical visually vetted counts are deliberately not acceptance criteria.
checks = {
    'unique_parent_entries': one_arcsecond.qso_group_id.is_unique,
    'known_scoring_routes': one_arcsecond.scoring_route.isin([
        'resolved_original', 'resolved_supplemental', 'unresolved', 'currently_unscoreable']).all(),
    'threshold_calibrated': calibration['model'] == 'either' and np.isfinite(calibration['adopted_log_r_threshold']),
    'calibration_before_pm_selection': calibration['pm_significance_upper_limit'] == 70,
    'only_six_products': ({path.name for path in product_entries} == expected_products
                          and all(path.is_file() for path in product_entries)),
    'dcr_once_one_gaia_only': (imaging.dcr_cut_applied.equals(imaging.candidate_tier.eq('one_gaia'))
                              and imaging.loc[imaging.candidate_tier.eq('one_gaia'), 'pair_dcr_delta_deg'].gt(20).all()
                              and imaging_cutflow['Selection stage'].str.contains('DCR').sum() == 1
                              and 'joint_dcr_pass' not in imaging),
    'secondary_pm_flags': bool(pm_flags_consistent),
    'selected_secondary_pm': bool(joint.secondary_pm_pass.all()),
    'redshift_error_metadata': (joint.adopted_redshift_error_used.notna().all()
                                and joint.loc[joint.redshift_catalogue.str.contains('quaia', case=False),
                                              'adopted_redshift_error_used'].gt(0).all()),
    'unique_final_candidates': joint.candidate_id.is_unique,
    'confirmed_contaminants_removed': not joint.literature_contaminant_removed.any(),
    'visual_inspection_outside_workflow': not any(
        'visual_inspection' in column for frame in (one_arcsecond, imaging, joint) for column in frame.columns),
    'automatic_imaging_selection': imaging.joint_final_pass.equals(imaging.joint_morphology_pass),
}
skipped_checks = []
one_summary_path = ONE_ARCSECOND_SEARCH_ROOT/'search_summary.json'
if one_summary_path.is_file() and ONE_ARCSECOND_PARENT_PATH.is_file():
    discovery = json.loads(one_summary_path.read_text())
    parent = pd.read_csv(ONE_ARCSECOND_PARENT_PATH, dtype={'qso_group_id': 'string', 'matched_gaia_source_id': 'string'})
    counts = discovery['counts']
    checks['one_arcsecond_candidate_discovery'] = (
        discovery['network_access'] is False
        and counts['parent_entries'] == len(parent)
        and parent.gaia_edr3_cnt1.gt(1).all()
        and set(parent.qso_group_id) == set(one_arcsecond.qso_group_id))
    checks['one_arcsecond_discovery_output_hash'] = hashlib.sha256(ONE_ARCSECOND_PARENT_PATH.read_bytes()).hexdigest() == discovery['output_sha256']
else:
    skipped_checks.append('Detailed one-arcsecond discovery: requires its generated parent and search summary.')
search_summary_path = IMAGING_SEARCH_ROOT / 'search_summary.json'
if search_summary_path.is_file():
    search_summary = json.loads(search_summary_path.read_text())
    search_counts = search_summary['counts']
    checks['imaging_candidate_search_counts'] = (
        search_summary['stage'] == 'all' and search_summary['network_access'] is False
        and search_counts['primary_has_gaia'] == search_counts['one_gaia'] + search_counts['two_gaia']
        and search_counts['initial_psf_pairs'] >= search_counts['validated_psf_pairs'] >= search_counts['primary_has_gaia'])
    checks['imaging_search_generated_inputs'] = all(
        (IMAGING_SEARCH_ROOT/name).is_file()
        and hashlib.sha256((IMAGING_SEARCH_ROOT/name).read_bytes()).hexdigest() == expected
        for name,expected in search_summary['output_sha256'].items())
else:
    skipped_checks.append('Detailed imaging candidate search: requires the search working tables.')
literature_audit_path = LITERATURE_WORK_ROOT / 'candidate_audit.csv'
literature_summary_path = LITERATURE_WORK_ROOT / 'summary.json'
if literature_audit_path.is_file() and literature_summary_path.is_file():
    literature_audit = pd.read_csv(literature_audit_path, dtype=IDENTIFIER_DTYPES)
    literature_summary = json.loads(literature_summary_path.read_text())
    removed = literature_audit.literature_contaminant_removed
    checks['literature_exclusion_audit'] = bool(
        len(literature_audit) - int(removed.sum()) == len(joint)
        and set(literature_audit.loc[~removed, 'candidate_id']) == set(joint.candidate_id)
        and literature_summary['counts']['final_candidates'] == len(joint))
    checks['literature_snapshot_hashes'] = all(
        hashlib.sha256((LITERATURE_ROOT / name).read_bytes()).hexdigest() == expected
        for name, expected in literature_summary['reference_sha256'].items())
else:
    skipped_checks.append('Detailed literature associations: requires the literature working audit.')
one_cutflow = pd.read_csv(ONE_ARCSECOND_CUTFLOW_PATH)
literature_stage = 'Confirmed literature contaminant removal'
one_final_rows = one_cutflow.loc[one_cutflow['Selection stage'].eq(literature_stage)]
imaging_final_row = imaging_cutflow.loc[imaging_cutflow['Selection stage'].eq(literature_stage)]
checks['literature_cutflow_counts'] = (
    one_final_rows.set_index('Route')['Either model'].to_dict() == {
        'Resolved': int(one_arcsecond.literature_final_pass.where(one_arcsecond.selected_by_resolved, False).sum()),
        'Unresolved': int(one_arcsecond.literature_final_pass.where(one_arcsecond.selected_by_unresolved, False).sum())}
    and len(imaging_final_row) == 1
    and int(imaging_final_row.Total.iloc[0]) == int(imaging.literature_final_pass.sum()))
fit_path = ONE_ROOT / 'threshold_fit_candidates.csv'
if fit_path.is_file():
    fit = pd.read_csv(fit_path)
    checks['calibration_sample_before_secondary_pm'] = bool(
        len(fit) == calibration['fit_systems'] and fit.secondary_pm_significance.lt(70).all()
        )
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
    checks['blend_grid_likelihoods'] = bool(grids_match)
else:
    skipped_checks.append('Detailed blend likelihood grids: requires working per-model scores.')

figures = [
    'one_arcsecond/resolved_pm_calibration.pdf',
    'one_arcsecond/resolved_colour_redshift.pdf',
    'imaging/colour_redshift.pdf',
    'imaging/pre_morphology_separation.pdf',
    'joint/galactic_latitude.pdf',
    'joint/separation_redshift.pdf',
]
for name in figures:
    path = FIGURE_ROOT / name
    checks[name] = path.is_file() and path.stat().st_size > 1000
summary = {'passed': all(checks.values()), 'checks': {name: bool(value) for name, value in checks.items()}}
summary['skipped_checks'] = skipped_checks
print(json.dumps(summary, indent=2))
raise SystemExit(0 if summary['passed'] else 1)
