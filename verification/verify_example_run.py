"""Check products, scientific flags and working audits for the requested pipelines."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.special import logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dual_agn.settings import (FIGURE_ROOT, IDENTIFIER_DTYPES, ONE_ROOT, OUTPUT_ROOT,
    FINAL_CANDIDATES_PATH, ONE_ARCSECOND_SCORES_PATH, IMAGING_SCORES_PATH,
    IMAGING_CUTFLOW_PATH, THRESHOLD_CALIBRATION_PATH, PM_SIGNIFICANCE_MAX,
    IMAGING_SEARCH_ROOT, LITERATURE_ROOT, LITERATURE_WORK_ROOT,
    ONE_ARCSECOND_PARENT_PATH, ONE_ARCSECOND_SEARCH_ROOT, ONE_ARCSECOND_CUTFLOW_PATH)
from dual_agn.score_selection import build_secondary_pm_flags, has_two_distinct_gaia_sources
from dual_agn.workflow import PIPELINE_CHOICES, product_names, selected_pipelines


def checksum(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_outputs(pipeline='both'):
    selected = selected_pipelines(pipeline)
    calibration = json.loads(THRESHOLD_CALIBRATION_PATH.read_text())
    joint = pd.read_csv(FINAL_CANDIDATES_PATH, dtype=IDENTIFIER_DTYPES)
    products = [path for path in OUTPUT_ROOT.iterdir() if path.name != '.DS_Store']
    checks = {
        'expected_products': ({path.name for path in products} == product_names(pipeline)
                              and all(path.is_file() for path in products)),
        'threshold_calibrated': calibration['model'] == 'either' and np.isfinite(calibration['adopted_log_r_threshold']),
        'calibration_before_pm_selection': calibration['pm_significance_upper_limit'] == 70,
        'unique_final_candidates': joint.candidate_id.is_unique,
        'selected_secondary_pm': bool(joint.secondary_pm_pass.all()),
        'confirmed_contaminants_removed': not joint.literature_contaminant_removed.any(),
        'redshift_error_metadata': (joint.adopted_redshift_error_used.notna().all()
            and joint.loc[joint.redshift_catalogue.str.contains('quaia', case=False, na=False),
                          'adopted_redshift_error_used'].gt(0).all()),
    }
    if pipeline == 'imaging':
        checks['requested_pipeline_only'] = bool(joint.selected_by_imaging.all()
            and not (joint.selected_by_resolved | joint.selected_by_unresolved).any())
    elif pipeline == 'gaia_neighbour':
        checks['requested_pipeline_only'] = bool(not joint.selected_by_imaging.any()
            and (joint.selected_by_resolved | joint.selected_by_unresolved).all())
    frames = [joint]
    pm_inputs = []
    figures = ['joint/galactic_latitude.pdf', 'joint/separation_redshift.pdf']
    literature_stage = 'Confirmed literature contaminant removal'
    if 'gaia_neighbour' in selected:
        one = pd.read_csv(ONE_ARCSECOND_SCORES_PATH, dtype=IDENTIFIER_DTYPES, low_memory=False)
        one_cutflow = pd.read_csv(ONE_ARCSECOND_CUTFLOW_PATH)
        frames.append(one)
        pm_inputs.append((one, 'primary_component_gaia_source_id', 'secondary_component_gaia_source_id', 'secondary_pm_significance'))
        checks['unique_parent_entries'] = one.qso_group_id.is_unique
        checks['known_scoring_routes'] = one.scoring_route.isin([
            'resolved_original', 'resolved_supplemental', 'unresolved', 'currently_unscoreable']).all()
        final_rows = one_cutflow.loc[one_cutflow['Selection stage'].eq(literature_stage)]
        checks['gaia_literature_cutflow_counts'] = final_rows.set_index('Route')['Either model'].to_dict() == {
            'Resolved': int(one.literature_final_pass.where(one.selected_by_resolved, False).sum()),
            'Unresolved': int(one.literature_final_pass.where(one.selected_by_unresolved, False).sum())}
        figures += ['one_arcsecond/resolved_pm_calibration.pdf', 'one_arcsecond/resolved_colour_redshift.pdf']
    if 'imaging' in selected:
        imaging = pd.read_csv(IMAGING_SCORES_PATH, dtype=IDENTIFIER_DTYPES, low_memory=False)
        cutflow = pd.read_csv(IMAGING_CUTFLOW_PATH)
        frames.append(imaging)
        pm_inputs.append((imaging, 'primary_gaia_source_id', 'secondary_gaia_source_id', 'secondary_pm_sig'))
        checks['dcr_once_one_gaia_only'] = (imaging.dcr_cut_applied.equals(imaging.candidate_tier.eq('one_gaia'))
            and imaging.loc[imaging.candidate_tier.eq('one_gaia'), 'pair_dcr_delta_deg'].gt(20).all()
            and cutflow['Selection stage'].str.contains('DCR').sum() == 1 and 'joint_dcr_pass' not in imaging)
        checks['automatic_imaging_selection'] = imaging.joint_final_pass.equals(imaging.joint_morphology_pass)
        final_row = cutflow.loc[cutflow['Selection stage'].eq(literature_stage)]
        checks['imaging_literature_cutflow_counts'] = len(final_row) == 1 and int(final_row.Total.iloc[0]) == int(imaging.literature_final_pass.sum())
        search = json.loads((IMAGING_SEARCH_ROOT/'search_summary.json').read_text())
        counts = search['counts']
        checks['imaging_candidate_search_counts'] = (search['stage'] == 'all' and search['network_access'] is False
            and counts['primary_has_gaia'] == counts['one_gaia'] + counts['two_gaia']
            and counts['initial_psf_pairs'] >= counts['validated_psf_pairs'] >= counts['primary_has_gaia'])
        checks['imaging_search_generated_inputs'] = all(
            (IMAGING_SEARCH_ROOT/name).is_file() and checksum(IMAGING_SEARCH_ROOT/name) == expected
            for name, expected in search['output_sha256'].items())
        figures += ['imaging/colour_redshift.pdf', 'imaging/pre_morphology_separation.pdf']
    for systems, primary, secondary, significance in pm_inputs:
        expected = build_secondary_pm_flags(systems[significance],
            has_two_distinct_gaia_sources(systems[primary], systems[secondary]), PM_SIGNIFICANCE_MAX)
        checks['secondary_pm_flags_' + primary] = all(expected[column].equals(systems[column])
            for column in ('has_two_distinct_gaia_sources', 'secondary_pm_pass', 'secondary_pm_missing'))
    checks['visual_inspection_outside_workflow'] = not any(
        'visual_inspection' in column for frame in frames for column in frame.columns)

    # Gaia discovery supplies calibration even in an imaging-only run.
    discovery = json.loads((ONE_ARCSECOND_SEARCH_ROOT/'search_summary.json').read_text())
    parent = pd.read_csv(ONE_ARCSECOND_PARENT_PATH, dtype={'qso_group_id':'string'})
    scoring = pd.read_csv(ONE_ROOT/'scoring_inputs.csv', usecols=['qso_group_id'], dtype={'qso_group_id':'string'})
    checks['gaia_discovery_for_calibration'] = (discovery['network_access'] is False
        and discovery['counts']['parent_entries'] == len(parent) and parent.gaia_edr3_cnt1.gt(1).all()
        and set(parent.qso_group_id) == set(scoring.qso_group_id)
        and checksum(ONE_ARCSECOND_PARENT_PATH) == discovery['output_sha256'])
    fit = pd.read_csv(ONE_ROOT/'threshold_fit_candidates.csv')
    checks['calibration_sample_before_secondary_pm'] = (len(fit) == calibration['fit_systems']
        and fit.secondary_pm_significance.lt(70).all())

    audit = pd.read_csv(LITERATURE_WORK_ROOT/'candidate_audit.csv', dtype=IDENTIFIER_DTYPES)
    literature = json.loads((LITERATURE_WORK_ROOT/'summary.json').read_text())
    removed = audit.literature_contaminant_removed
    checks['literature_exclusion_audit'] = (len(audit) - int(removed.sum()) == len(joint)
        and set(audit.loc[~removed, 'candidate_id']) == set(joint.candidate_id)
        and literature['counts']['final_candidates'] == len(joint) and literature.get('pipeline', 'both') == pipeline)
    checks['literature_snapshot_hashes'] = all(checksum(LITERATURE_ROOT/name) == expected
        for name, expected in literature['reference_sha256'].items())
    if 'gaia_neighbour' in selected:
        blends = pd.read_csv(ONE_ROOT/'blend_scores_all_nonresolved_by_model_and_band_set.csv')
        grids_match = True
        for row in blends.loc[blends.redshift_error.gt(0)].itertuples():
            grid = pd.DataFrame(json.loads(row.redshift_grid_json))
            finite = np.isfinite(grid.log_like_qq) & np.isfinite(grid.log_like_qs)
            expected = [logsumexp(grid.loc[finite, field] + np.log(grid.loc[finite, 'weight']))
                        for field in ('log_like_qq', 'log_like_qs')]
            grids_match &= (len(grid) == 9 and np.isclose(grid.weight.sum(), 1.0)
                and np.allclose(expected, [row.log_like_qq, row.log_like_qs], atol=1e-10, rtol=0))
        checks['blend_grid_likelihoods'] = bool(grids_match)
    for name in figures:
        path = FIGURE_ROOT/name
        checks[name] = path.is_file() and path.stat().st_size > 1000
    return {'pipeline': pipeline, 'passed': bool(all(checks.values())),
            'checks': {name: bool(value) for name, value in checks.items()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pipeline', choices=PIPELINE_CHOICES, default='both')
    args = parser.parse_args()
    summary = verify_outputs(args.pipeline)
    print(json.dumps(summary, indent=2))
    raise SystemExit(0 if summary['passed'] else 1)


if __name__ == '__main__':
    main()
