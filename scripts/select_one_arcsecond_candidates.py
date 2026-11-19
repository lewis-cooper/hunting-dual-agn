from pathlib import Path
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.table import Table
from settings import *
from photometry import read_fits
from one_arcsecond import combine_current_evidence
from score_selection import (build_secondary_pm_flags, eligible_model_scores,
                             has_two_distinct_gaia_sources, read_either_model_threshold)
log_r_per_unit_z_threshold = read_either_model_threshold(THRESHOLD_CALIBRATION_PATH)
systems = pd.read_csv(
    ONE_ROOT / 'scoring_input_audit_712.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)
resolved_all = pd.read_csv(
    ONE_ROOT / 'ranked_resolved_extended.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)
blend_ranked = pd.read_csv(
    ONE_ROOT / 'ranked_unresolved.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)
evidence = combine_current_evidence(systems, resolved_all, blend_ranked)
duplicate_columns = evidence.columns.intersection(systems.columns).drop('qso_group_id')
one_arcsecond = evidence.merge(
    systems.drop(columns=duplicate_columns),
    on='qso_group_id',
    how='left',
    validate='one_to_one'
)
one_arcsecond = one_arcsecond.copy()
visual_rejections = pd.read_csv(VISUAL_REJECTIONS_PATH, dtype='string').fillna('')
rejected = visual_rejections.loc[visual_rejections.one_arcsecond_qso_group_id.ne('')].set_index('one_arcsecond_qso_group_id')
visual_inspection_key = one_arcsecond.get('historical_qso_group_id', one_arcsecond.qso_group_id)
one_arcsecond['visual_inspection_rejected'] = visual_inspection_key.astype('string').isin(rejected.index)
one_arcsecond['visual_inspection_rejection_reason'] = visual_inspection_key.map(rejected.reason).fillna('')
two_gaia = has_two_distinct_gaia_sources(one_arcsecond.primary_component_gaia_source_id,
                                        one_arcsecond.secondary_component_gaia_source_id)
pm_flags = build_secondary_pm_flags(one_arcsecond.secondary_pm_significance,
                                    two_gaia, PM_SIGNIFICANCE_MAX)
one_arcsecond[pm_flags.columns] = pm_flags
resolved_route = one_arcsecond.scoring_route.isin(['resolved_original', 'resolved_supplemental'])
unresolved_route = one_arcsecond.scoring_route.eq('unresolved')
unresolved_finite_both = unresolved_route & np.isfinite(one_arcsecond.blend_log_bf_qq_qs_magdep) & np.isfinite(one_arcsecond.blend_log_bf_qq_qs_magindep)
qq_either = unresolved_finite_both & (one_arcsecond.blend_log_bf_qq_qs_magdep.gt(0) | one_arcsecond.blend_log_bf_qq_qs_magindep.gt(0))
qq_both = unresolved_finite_both & (one_arcsecond.blend_log_bf_qq_qs_magdep.gt(0) & one_arcsecond.blend_log_bf_qq_qs_magindep.gt(0))
qq_either &= one_arcsecond.secondary_pm_pass & ~one_arcsecond.visual_inspection_rejected
qq_both &= one_arcsecond.secondary_pm_pass & ~one_arcsecond.visual_inspection_rejected
one_arcsecond['selected_by_unresolved'] = qq_either
one_arcsecond['qq_both_models'] = qq_both

def model_counts(route, stage, magdep, magindep):
    return {
        'Route': route,
        'Selection stage': stage,
        'Magnitude dependent': int(magdep.sum()),
        'Magnitude independent': int(magindep.sum()),
        'Either model': int((magdep | magindep).sum())
    }
quality_metadata = pd.DataFrame(
    {
        'reference_snr': one_arcsecond.secondary_decals_flux_r * np.sqrt(one_arcsecond.secondary_decals_flux_ivar_r),
        'abs_b_deg': np.abs(SkyCoord(
            one_arcsecond.secondary_component_ra.to_numpy(float) * u.deg,
            one_arcsecond.secondary_component_dec.to_numpy(float) * u.deg
        ).galactic.b.deg),
        'secondary_decals_type': one_arcsecond.secondary_decals_type
    },
    index=one_arcsecond.index
)
resolved_cumulative = resolved_route.copy()
one_arcsecond_count_rows = [model_counts('Resolved', 'Assigned to route', resolved_route, resolved_route)]
for label, selection in [
    (
        'Companion reference-band S/N >= 10',
        np.isfinite(quality_metadata.reference_snr) & quality_metadata.reference_snr.ge(10)
    ),
    (
        '|b| >= 25 degrees',
        np.isfinite(quality_metadata.abs_b_deg) & quality_metadata.abs_b_deg.ge(25)
    ),
    ('Companion PSF morphology', quality_metadata.secondary_decals_type.eq('PSF'))
]:
    resolved_cumulative &= selection
    one_arcsecond_count_rows.append(model_counts('Resolved', label, resolved_cumulative, resolved_cumulative))
resolved_magdep = resolved_cumulative & one_arcsecond.resolved_eligible_magdep.eq(True)
resolved_magindep = resolved_cumulative & one_arcsecond.resolved_eligible_magindep.eq(True)
one_arcsecond_count_rows.append(model_counts(
    'Resolved',
    'Remaining model eligibility/support requirements',
    resolved_magdep,
    resolved_magindep
))
resolved_scores = eligible_model_scores(one_arcsecond, 'resolved_eligible', 'resolved_log_r_per_unit_z')
one_arcsecond['resolved_log_r_per_unit_z_either'] = resolved_scores.max(axis=1).where(resolved_route)
resolved_magdep &= resolved_scores.magdep.notna()
resolved_magindep &= resolved_scores.magindep.notna()
one_arcsecond_count_rows.append(model_counts('Resolved', 'Finite ln R', resolved_magdep, resolved_magindep))
resolved_magdep &= resolved_scores.magdep.ge(log_r_per_unit_z_threshold)
resolved_magindep &= resolved_scores.magindep.ge(log_r_per_unit_z_threshold)
one_arcsecond_count_rows.append(model_counts(
    'Resolved',
    f'ln R >= {log_r_per_unit_z_threshold:.2f}',
    resolved_magdep,
    resolved_magindep
))
resolved_magdep &= one_arcsecond.secondary_pm_pass
resolved_magindep &= one_arcsecond.secondary_pm_pass
one_arcsecond_count_rows.append(model_counts(
    'Resolved', 'Secondary Gaia PM significance < 3 (missing retained)',
    resolved_magdep, resolved_magindep))
resolved_magdep &= ~one_arcsecond.visual_inspection_rejected
resolved_magindep &= ~one_arcsecond.visual_inspection_rejected
one_arcsecond_count_rows.append(model_counts('Resolved', 'Visual inspection', resolved_magdep, resolved_magindep))
one_arcsecond['selected_by_resolved_magdep'] = resolved_magdep
one_arcsecond['selected_by_resolved_magindep'] = resolved_magindep
resolved_selected = resolved_magdep | resolved_magindep
one_arcsecond['selected_by_resolved'] = resolved_selected
one_arcsecond_selected = one_arcsecond.loc[resolved_selected | qq_either].copy()
print(f'Resolved: {resolved_selected.sum()}; unresolved: {qq_either.sum()}; one-arcsecond total: {len(one_arcsecond_selected)}.')
blend_magdep = unresolved_route.copy()
blend_magindep = unresolved_route.copy()
one_arcsecond_count_rows.append(model_counts('Unresolved', 'Assigned to route', blend_magdep, blend_magindep))
latitude_pass = np.isfinite(one_arcsecond.blend_galactic_abs_b_deg) & one_arcsecond.blend_galactic_abs_b_deg.ge(25)
blend_magdep &= latitude_pass
blend_magindep &= latitude_pass
one_arcsecond_count_rows.append(model_counts('Unresolved', '|b| >= 25 degrees', blend_magdep, blend_magindep))
for label, failed_status in [
    ('Usable reference-band photometry', 'no_reference_band'),
    ('Redshift within model support', 'redshift_outside_support'),
    ('Reference-band S/N >= 20', 'too_faint_for_two_components')
]:
    blend_magdep &= one_arcsecond.blend_status_magdep.ne(failed_status)
    blend_magindep &= one_arcsecond.blend_status_magindep.ne(failed_status)
    one_arcsecond_count_rows.append(model_counts('Unresolved', label, blend_magdep, blend_magindep))
blend_magdep &= unresolved_finite_both
blend_magindep &= unresolved_finite_both
one_arcsecond_count_rows.append(model_counts('Unresolved', 'Finite QQ-QS log Bayes factor', blend_magdep, blend_magindep))
blend_magdep &= one_arcsecond.blend_log_bf_qq_qs_magdep.gt(0)
blend_magindep &= one_arcsecond.blend_log_bf_qq_qs_magindep.gt(0)
one_arcsecond_count_rows.append(model_counts('Unresolved', 'QQ favoured', blend_magdep, blend_magindep))
blend_magdep &= one_arcsecond.secondary_pm_pass
blend_magindep &= one_arcsecond.secondary_pm_pass
one_arcsecond_count_rows.append(model_counts(
    'Unresolved', 'Secondary Gaia PM significance < 3 (missing retained)',
    blend_magdep, blend_magindep))
blend_magdep &= ~one_arcsecond.visual_inspection_rejected
blend_magindep &= ~one_arcsecond.visual_inspection_rejected
one_arcsecond_count_rows.append(model_counts('Unresolved', 'Visual inspection', blend_magdep, blend_magindep))
assert (blend_magdep | blend_magindep).equals(qq_either), 'Unresolved cutflow and adopted selection differ'
one_arcsecond_cutflow = pd.DataFrame(one_arcsecond_count_rows)
print('Unscoreable parent systems:', (~(resolved_route | unresolved_route)).sum())
print(
    'Adopted one-arcsecond selection, either model in each route:',
    len(one_arcsecond_selected)
)
one_arcsecond.to_csv(ONE_ARCSECOND_SCORES_PATH, index=False)
one_arcsecond_cutflow.to_csv(ONE_ARCSECOND_CUTFLOW_PATH, index=False)
