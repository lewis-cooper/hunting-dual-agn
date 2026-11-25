import argparse
from dual_agn.redshift_inputs import redshift_errors
from pathlib import Path
import json
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18
from astropy.table import Table
from dual_agn.settings import *
from dual_agn.workflow import PIPELINE_CHOICES, selected_pipelines


def gaia_candidate_records(scores):
    one_arcsecond_selected = scores.loc[scores.selected_by_resolved | scores.selected_by_unresolved].copy()
    one_arcsecond_points = pd.DataFrame({
        'source_candidate_id': 'one_arcsecond:' + one_arcsecond_selected.qso_group_id,
        'pipeline': 'one_arcsecond',
        'selection_route': np.where(one_arcsecond_selected.selected_by_resolved, 'resolved', 'unresolved'),
        'system_name': one_arcsecond_selected.primary_qso_name,
        'input_redshift': one_arcsecond_selected.primary_qso_z,
        'original_redshift': one_arcsecond_selected.primary_qso_z_original,
        'input_redshift_release': one_arcsecond_selected.redshift_catalogue_release,
        'input_redshift_catalogue': one_arcsecond_selected.primary_qso_catalog,
        'redshift_error_used': redshift_errors(one_arcsecond_selected, 'primary_qso_zerr', 'primary_qso_catalog', 'primary_qso_redshift_kind'),
        'secondary_pm_significance': one_arcsecond_selected.secondary_pm_significance,
        'secondary_pm_missing': one_arcsecond_selected.secondary_pm_missing,
        'secondary_pm_pass': one_arcsecond_selected.secondary_pm_pass,
        'qso_group_id': one_arcsecond_selected.qso_group_id,
        'resolved_log_r_magdep': one_arcsecond_selected.resolved_log_r_per_unit_z_magdep,
        'resolved_log_r_magindep': one_arcsecond_selected.resolved_log_r_per_unit_z_magindep,
        'resolved_log_r_either': one_arcsecond_selected.resolved_log_r_per_unit_z_either,
        'resolved_pass_magdep': one_arcsecond_selected.selected_by_resolved_magdep,
        'resolved_pass_magindep': one_arcsecond_selected.selected_by_resolved_magindep,
        'blend_log_bf_qq_qs_magdep': one_arcsecond_selected.blend_log_bf_qq_qs_magdep,
        'blend_log_bf_qq_qs_magindep': one_arcsecond_selected.blend_log_bf_qq_qs_magindep
    })
    for role in ('primary', 'secondary'):
        for coordinate in ('ra', 'dec'):
            one_arcsecond_points[f'{role}_{coordinate}_deg'] = one_arcsecond_selected[f'{role}_component_{coordinate}']
        one_arcsecond_points[f'{role}_gaia_source_id'] = one_arcsecond_selected[f'{role}_component_gaia_source_id']
    return one_arcsecond_points


def imaging_candidate_records(scores):
    imaging_selected = scores.loc[scores.joint_final_pass].copy()
    imaging_points = pd.DataFrame({
        'source_candidate_id': 'imagery:' + imaging_selected.qso_parent_index.astype(str),
        'pipeline': 'imaging',
        'selection_route': imaging_selected.candidate_tier,
        'system_name': imaging_selected.sdss_name,
        'input_redshift': imaging_selected.adopted_redshift,
        'input_redshift_catalogue': imaging_selected.adopted_redshift_catalogue,
        'redshift_error_used': imaging_selected.adopted_redshift_error,
        'secondary_pm_significance': imaging_selected.secondary_pm_significance,
        'secondary_pm_missing': imaging_selected.secondary_pm_missing,
        'secondary_pm_pass': imaging_selected.secondary_pm_pass,
        'qso_parent_index': imaging_selected.qso_parent_index,
        'n_psf_components_within_aperture': imaging_selected.n_psf_components_within_aperture,
        'imaging_primary_log_r_magdep': imaging_selected.primary_log_r_per_unit_z_magdep,
        'imaging_secondary_log_r_magdep': imaging_selected.secondary_log_r_per_unit_z_magdep,
        'imaging_primary_log_r_magindep': imaging_selected.primary_log_r_per_unit_z_magindep,
        'imaging_secondary_log_r_magindep': imaging_selected.secondary_log_r_per_unit_z_magindep,
        'imaging_minimum_component_log_r_either': imaging_selected.minimum_component_log_r_either,
        'imaging_pass_magdep': imaging_selected.passes_pair_score_magdep,
        'imaging_pass_magindep': imaging_selected.passes_pair_score_magindep,
        'imaging_pair_dcr_angle_deg': imaging_selected.pair_dcr_delta_deg,
        'imaging_dcr_cut_applied': imaging_selected.dcr_cut_applied,
        'imaging_r_perp_sigma_ratio': imaging_selected.r_perp_sigma_ratio,
        'imaging_morphology_source': imaging_selected.morphology_measurement_source
    })
    for role in ('primary', 'secondary'):
        for coordinate in ('ra', 'dec'):
            imaging_points[f'{role}_{coordinate}_deg'] = imaging_selected[f'{role}_{coordinate}']
        imaging_points[f'{role}_gaia_source_id'] = imaging_selected[f'{role}_gaia_source_id']
    return imaging_points


def combine_candidates(pipeline='both'):
    selected = selected_pipelines(pipeline)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    records = []
    if 'gaia_neighbour' in selected:
        scores = pd.read_csv(ONE_ARCSECOND_SCORES_PATH, dtype=IDENTIFIER_DTYPES, float_precision='round_trip')
        records.append(gaia_candidate_records(scores))
    if 'imaging' in selected:
        scores = pd.read_csv(IMAGING_SCORES_PATH, dtype=IDENTIFIER_DTYPES, float_precision='round_trip')
        records.append(imaging_candidate_records(scores))
    pair_match_tolerance_arcsec = PAIR_MATCH_TOLERANCE_ARCSEC
    source_records = pd.concat(records, ignore_index=True)
    source_records['qso_parent_index'] = source_records.get('qso_parent_index', pd.Series(pd.NA, index=source_records.index, dtype='Int64')).astype('Int64')
    coordinate_columns = [f'{role}_{coordinate}_deg' for role in ('primary', 'secondary') for coordinate in ('ra', 'dec')]
    print('Selected input records before cross-pipeline deduplication:', len(source_records))
    primary_sky = SkyCoord(
        source_records.primary_ra_deg.to_numpy(dtype=float) * u.deg,
        source_records.primary_dec_deg.to_numpy(dtype=float) * u.deg
    )
    secondary_sky = SkyCoord(
        source_records.secondary_ra_deg.to_numpy(dtype=float) * u.deg,
        source_records.secondary_dec_deg.to_numpy(dtype=float) * u.deg
    )
    aligned_offsets = np.maximum(
        primary_sky[:, None].separation(primary_sky[None, :]).arcsec,
        secondary_sky[:, None].separation(secondary_sky[None, :]).arcsec
    )
    swapped_offsets = np.maximum(
        primary_sky[:, None].separation(secondary_sky[None, :]).arcsec,
        secondary_sky[:, None].separation(primary_sky[None, :]).arcsec
    )
    pair_offsets = np.minimum(aligned_offsets, swapped_offsets)
    first_indices, second_indices = np.where(np.triu(pair_offsets <= pair_match_tolerance_arcsec, k=1))
    pair_matches = pd.DataFrame({
        'one_arcsecond_id': source_records.source_candidate_id.iloc[first_indices].to_numpy(),
        'imaging_id': source_records.source_candidate_id.iloc[second_indices].to_numpy(),
        'maximum_component_offset_arcsec': pair_offsets[first_indices, second_indices],
        'components_swapped': swapped_offsets[first_indices, second_indices] < aligned_offsets[first_indices, second_indices]
    })
    if (not pair_matches.one_arcsecond_id.str.startswith('one_arcsecond:').all()
            or not pair_matches.imaging_id.str.startswith('imagery:').all()):
        raise ValueError('Within-search duplicate coordinates must be deduplicated upstream')
    source_records['candidate_id'] = 'joint:' + source_records.source_candidate_id
    shared_id_map = pair_matches.set_index('imaging_id').one_arcsecond_id
    shared_imaging = source_records.source_candidate_id.isin(shared_id_map.index)
    if shared_imaging.any():
        source_records.loc[shared_imaging, 'candidate_id'] = 'joint:' + source_records.loc[shared_imaging, 'source_candidate_id'].map(shared_id_map)
    pair_matches['candidate_id'] = 'joint:' + pair_matches.one_arcsecond_id.astype('string')
    joint_candidates = source_records.loc[
        ~shared_imaging,
        [
            'candidate_id',
            'system_name',
            'input_redshift',
            'input_redshift_catalogue',
            'primary_gaia_source_id',
            'secondary_gaia_source_id',
            *coordinate_columns
        ]
    ].rename(columns={'input_redshift': 'redshift', 'input_redshift_catalogue': 'redshift_catalogue'}).copy()
    one_records = source_records.loc[source_records.pipeline.eq('one_arcsecond')].set_index('candidate_id')
    image_records = source_records.loc[source_records.pipeline.eq('imaging')].set_index('candidate_id')
    joint_candidates['n_psf_components_within_aperture'] = joint_candidates.candidate_id.map(
        image_records.n_psf_components_within_aperture
        if 'n_psf_components_within_aperture' in image_records else pd.Series(dtype='Int64')).astype('Int64')
    for column in (
        'source_candidate_id',
        'qso_group_id',
        'input_redshift',
        'original_redshift',
        'input_redshift_release',
        'input_redshift_catalogue',
        'redshift_error_used',
        'secondary_pm_significance',
        'secondary_pm_missing',
        'secondary_pm_pass',
        'selection_route',
        'resolved_log_r_magdep',
        'resolved_log_r_magindep',
        'resolved_log_r_either',
        'resolved_pass_magdep',
        'resolved_pass_magindep',
        'blend_log_bf_qq_qs_magdep',
        'blend_log_bf_qq_qs_magindep'
    ):
        joint_candidates['one_arcsecond_' + column] = joint_candidates.candidate_id.map(one_records[column] if column in one_records else pd.Series(dtype=float))
    for column in (
        'source_candidate_id',
        'qso_parent_index',
        'input_redshift',
        'input_redshift_catalogue',
        'redshift_error_used',
        'secondary_pm_significance',
        'secondary_pm_missing',
        'secondary_pm_pass',
        'selection_route',
        'imaging_primary_log_r_magdep',
        'imaging_secondary_log_r_magdep',
        'imaging_primary_log_r_magindep',
        'imaging_secondary_log_r_magindep',
        'imaging_minimum_component_log_r_either',
        'imaging_pass_magdep',
        'imaging_pass_magindep',
        'imaging_pair_dcr_angle_deg',
        'imaging_dcr_cut_applied',
        'imaging_r_perp_sigma_ratio',
        'imaging_morphology_source'
    ):
        output_column = column if column.startswith('imaging_') else 'imaging_' + column
        joint_candidates[output_column] = joint_candidates.candidate_id.map(image_records[column] if column in image_records else pd.Series(dtype=float))
    joint_candidates['selected_by_resolved'] = joint_candidates.one_arcsecond_selection_route.eq('resolved')
    joint_candidates['selected_by_unresolved'] = joint_candidates.one_arcsecond_selection_route.eq('unresolved')
    joint_candidates['selected_by_imaging'] = joint_candidates.imaging_source_candidate_id.notna()
    joint_candidates['selected_by_both_searches'] = joint_candidates.one_arcsecond_source_candidate_id.notna() & joint_candidates.selected_by_imaging
    joint_candidates['pair_coordinate_source'] = np.where(
        joint_candidates.one_arcsecond_source_candidate_id.notna(),
        'Gaia DR3',
        'DECaLS Tractor'
    )
    joint_candidates['cross_pipeline_component_offset_arcsec'] = joint_candidates.candidate_id.map(pair_matches.set_index('candidate_id').maximum_component_offset_arcsec)
    for role in ('primary', 'secondary'):
        column = f'{role}_gaia_source_id'
        valid_source_id = joint_candidates[column].astype('string').str.fullmatch('[1-9][0-9]*', na=False)
        joint_candidates[column] = joint_candidates[column].where(valid_source_id)
    same_gaia_source = joint_candidates.primary_gaia_source_id.eq(joint_candidates.secondary_gaia_source_id)
    joint_candidates.loc[same_gaia_source.fillna(False), 'secondary_gaia_source_id'] = pd.NA
    joint_candidates['has_two_distinct_gaia_sources'] = joint_candidates.primary_gaia_source_id.notna() & joint_candidates.secondary_gaia_source_id.notna()
    for field in ('secondary_pm_significance', 'secondary_pm_missing', 'secondary_pm_pass'):
        joint_candidates[field] = joint_candidates['one_arcsecond_' + field].where(
            joint_candidates.one_arcsecond_source_candidate_id.notna(),
            joint_candidates['imaging_' + field])
    joint_candidates['secondary_pm_missing'] = joint_candidates.secondary_pm_missing.fillna(False).astype(bool)
    joint_candidates['secondary_pm_pass'] = joint_candidates.secondary_pm_pass.fillna(False).astype(bool)
    assert joint_candidates.secondary_pm_pass.all(), 'Selected systems must pass the secondary-PM policy'
    joint_candidates['plot_group'] = np.select(
        [
            joint_candidates.selected_by_both_searches,
            joint_candidates.selected_by_resolved,
            joint_candidates.selected_by_unresolved
        ],
        ['Both searches', 'One-arcsecond resolved', 'One-arcsecond unresolved'],
        default='Imaging only'
    )
    print(f'Cross-pipeline duplicates: {len(pair_matches)}; unique pairs: {len(joint_candidates)}')
    shared = joint_candidates.selected_by_both_searches
    redshift_values_differ = shared & ~np.isclose(
        joint_candidates.one_arcsecond_input_redshift,
        joint_candidates.imaging_input_redshift,
        rtol=0,
        atol=1e-06
    )
    both_photometric = joint_candidates.one_arcsecond_redshift_error_used.gt(0) & joint_candidates.imaging_redshift_error_used.gt(0)
    joint_candidates['redshift_values_differ'] = redshift_values_differ
    joint_candidates['photometric_redshift_discrepancy'] = redshift_values_differ & both_photometric
    joint_candidates['redshift_difference_imaging_minus_one_arcsecond'] = joint_candidates.imaging_input_redshift - joint_candidates.one_arcsecond_input_redshift
    photometric_discrepancies = joint_candidates.photometric_redshift_discrepancy
    joint_candidates['score_redshift_differs_from_adopted'] = False
    joint_candidates['adopted_redshift_origin'] = np.where(
        joint_candidates.one_arcsecond_source_candidate_id.notna(),
        'one_arcsecond_corrected_input',
        'imaging_master_input'
    )
    joint_candidates['adopted_redshift_error_used'] = np.where(
        joint_candidates.one_arcsecond_source_candidate_id.notna(),
        joint_candidates.one_arcsecond_redshift_error_used,
        joint_candidates.imaging_redshift_error_used)
    joint_candidates['redshift_treatment'] = np.where(
        joint_candidates.adopted_redshift_error_used.gt(0),
        'Gaussian_photometric_redshift_average', 'fixed_spectroscopic_redshift')
    primary_coordinates = SkyCoord(
        joint_candidates.primary_ra_deg.to_numpy(dtype=float) * u.deg,
        joint_candidates.primary_dec_deg.to_numpy(dtype=float) * u.deg
    )
    secondary_coordinates = SkyCoord(
        joint_candidates.secondary_ra_deg.to_numpy(dtype=float) * u.deg,
        joint_candidates.secondary_dec_deg.to_numpy(dtype=float) * u.deg
    )
    joint_candidates['separation_arcsec'] = primary_coordinates.separation(secondary_coordinates).arcsec
    joint_candidates['projected_separation_kpc'] = joint_candidates.separation_arcsec.to_numpy() * Planck18.kpc_proper_per_arcmin(joint_candidates.redshift.to_numpy()).value / 60
    joint_candidates = joint_candidates.sort_values(['primary_ra_deg', 'primary_dec_deg'], kind='stable').reset_index(drop=True)
    redshift_audit = joint_candidates.loc[
        joint_candidates.selected_by_both_searches,
        [
            'candidate_id',
            'one_arcsecond_qso_group_id',
            'imaging_qso_parent_index',
            'one_arcsecond_input_redshift',
            'imaging_input_redshift',
            'redshift',
            'one_arcsecond_input_redshift_catalogue',
            'imaging_input_redshift_catalogue',
            'redshift_difference_imaging_minus_one_arcsecond',
            'redshift_values_differ',
            'photometric_redshift_discrepancy',
            'adopted_redshift_origin'
        ]
    ].copy()
    joint_candidates['plot_group'] = np.select(
        [
            joint_candidates.selected_by_both_searches,
            joint_candidates.selected_by_resolved,
            joint_candidates.selected_by_unresolved,
            joint_candidates.imaging_selection_route.eq('one_gaia')
        ],
        [
            'Both searches',
            'One-arcsecond resolved',
            'One-arcsecond unresolved',
            'Imaging: one Gaia'
        ],
        default='Imaging: two Gaia'
    )
    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    joint_candidates.to_csv(PRE_LITERATURE_CANDIDATES_PATH, index=False)
    joint_candidates.to_csv(FINAL_CANDIDATES_PATH, index=False)
    print(
        'Output:',
        FINAL_CANDIDATES_PATH,
        '| selected detections:',
        len(source_records),
        '| shared:',
        len(pair_matches),
        '| unique:',
        len(joint_candidates)
    )

    return joint_candidates


def main():
    parser = argparse.ArgumentParser(description='Combine candidates from the requested pipelines.')
    parser.add_argument('--pipeline', choices=PIPELINE_CHOICES, default='both')
    args = parser.parse_args()
    combine_candidates(args.pipeline)


if __name__ == '__main__':
    main()
