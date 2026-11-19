from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from settings import IDENTIFIER_DTYPES, MATCH_HALF_WIDTH_KMS, CLOSE_SEPARATION_ARCSEC
from photometry import RESOLVED_POLICY
from redshift_uncertainty import score_components_with_redshift_errors
MODEL_POINTERS = {'magdep': '20261005_magdep_qlf', 'magindep': '20261005_magindep_qlf'}
QUAIA_RELEASE = '1.0.0'
QUAIA_DOI = '10.5281/zenodo.10403370'
QUAIA_G20_MD5 = '72531bc67bde1b08a69d5aeae03fb26e'

def attach_tractor_diagnostics(
    systems: pd.DataFrame,
    component_metadata: pd.DataFrame,
    tractor_sources: pd.DataFrame
) -> pd.DataFrame:
    """Join morphology and contamination by exact Tractor identity, preserving rows."""
    identity_columns = ['release', 'brickid', 'objid']
    morphology = component_metadata[['decals_release', 'decals_brickid', 'decals_objid', 'decals_type']].rename(columns=lambda column: column.removeprefix('decals_'))
    for column in identity_columns:
        morphology[column] = pd.to_numeric(morphology[column], errors='coerce').astype('Int64')
    morphology['type'] = morphology['type'].map(lambda value: value.decode().strip() if isinstance(value, bytes) else value)
    morphology = morphology.drop_duplicates()
    morphology = morphology.loc[morphology.release.isin([9010, 9011, 9012])]
    contamination = tractor_sources[identity_columns + ['fracflux_r']].drop_duplicates()
    for column in identity_columns:
        contamination[column] = pd.to_numeric(contamination[column], errors='raise').astype('Int64')
    diagnostics = morphology.merge(contamination, on=identity_columns, how='outer', validate='one_to_one')
    enriched = systems.copy()
    for component in ('primary', 'secondary'):
        for column in identity_columns:
            key = f'{component}_decals_{column}'
            enriched[key] = pd.to_numeric(enriched[key], errors='coerce').astype('Int64')
        renamed = diagnostics.rename(columns={column: f'{component}_decals_{column}' for column in diagnostics})
        enriched = enriched.merge(
            renamed,
            on=[f'{component}_decals_{column}' for column in identity_columns],
            how='left',
            validate='many_to_one',
            sort=False
        )
    return enriched


def score_resolved_candidates(candidates, models):
    """Score native companion fluxes with both external model bundles."""
    from qso_pcolor import Photometry, RedshiftMatch
    blend_policy = RESOLVED_POLICY
    minimum_separation_arcsec = CLOSE_SEPARATION_ARCSEC
    match_half_width_kms = MATCH_HALF_WIDTH_KMS
    score_rows = []
    for model_name, model in models.items():
        for system_name, group in candidates.groupby('photometric_system', sort=False):
            bands = tuple((f'{system_name}:{band}' for band in 'grz'))
            flux = group[[f'secondary_decals_flux_{band}' for band in 'grz']].to_numpy(float)
            inverse_variance = group[[f'secondary_decals_flux_ivar_{band}' for band in 'grz']].to_numpy(float)
            valid_measurement = np.isfinite(flux) & np.isfinite(inverse_variance)
            valid_measurement &= inverse_variance > 0
            variance = np.divide(1.0, inverse_variance, out=np.full_like(flux, np.inf), where=valid_measurement)
            photometry = Photometry(np.where(valid_measurement, flux, np.nan), variance, bands)
            redshift_errors = group.primary_qso_zerr.fillna(0.0).to_numpy(float)
            redshift_errors = np.where(group.primary_qso_catalog.eq('QUAIA_G20_V1'), redshift_errors, 0.0)
            scores, decision = score_components_with_redshift_errors(
                model, photometry, redshift_errors,
                ra_deg=group.secondary_component_ra.to_numpy(float),
                dec_deg=group.secondary_component_dec.to_numpy(float),
                z_primary=group.primary_qso_z.to_numpy(float),
                morphology=group.secondary_decals_type.fillna('UNKNOWN').to_numpy(str),
                separation_arcsec=group.component_separation_arcsec.to_numpy(float),
                fracflux=group.secondary_decals_fracflux_r.to_numpy(float),
                candidate_id=group.secondary_component_gaia_source_id.to_numpy(str),
                primary_id=group.primary_component_gaia_source_id.to_numpy(str),
                match=RedshiftMatch(half_width_kms=match_half_width_kms),
                blend_policy=blend_policy,
                ood_flag_sigma=4.0
            )
            quasar_probability = model.quasar_probability(scores)
            for index, (group_id, score) in enumerate(zip(group.qso_group_id, scores)):
                row = {} if score is None else score.as_row()
                eligible = bool(decision['eligible'][index])
                reason = str(decision['reason'][index])
                row.update(
                    qso_group_id=group_id,
                    primary_redshift_error_used=float(redshift_errors[index]),
                    model=model_name,
                    eligible=eligible,
                    decision_reason=reason,
                    status=reason or row.get('status', 'not_scored'),
                    qso_support_percentile=float(decision['qso_support']['percentile'][index]),
                    p_quasar=float(quasar_probability[index]) if eligible else np.nan,
                    scoring_model_id=row.pop('photometric_system', '')
                )
                for field in ('log_r_per_unit_z', 'p_sameq'):
                    row[field] = row.get(field, np.nan) if eligible else np.nan
                row.setdefault('quality_flags', '')
                separation = group.component_separation_arcsec.iloc[index]
                if np.isfinite(separation) and (0 <= separation < minimum_separation_arcsec):
                    flags = [flag for flag in row['quality_flags'].split(',') if flag]
                    row['quality_flags'] = ','.join(dict.fromkeys([*flags, 'blended', 'close_separation']))
                for field in ('ref_mag', 'qso_ood_sigma', 'dz_match_eff'):
                    row.setdefault(field, np.nan)
                score_rows.append(row)
    scored = pd.DataFrame(score_rows)
    return scored


def join_resolved_model_scores(candidates: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """Expose both variants; unsuffixed plotting fields explicitly use magdep."""
    fields = [
        'eligible',
        'decision_reason',
        'status',
        'log_r_per_unit_z',
        'p_sameq',
        'p_quasar',
        'quality_flags',
        'qso_support_percentile',
        'ref_mag',
        'qso_ood_sigma',
        'dz_match_eff'
    ]
    retained_columns = [column for column in candidates if column not in fields and (not column.startswith('rank_'))]
    ranked = candidates[retained_columns].copy()
    for model_name in MODEL_POINTERS:
        model_scores = scores.loc[scores.model.eq(model_name), ['qso_group_id', *fields]]
        model_scores = model_scores.rename(columns={field: f'{field}_{model_name}' for field in fields})
        ranked = ranked.merge(model_scores, on='qso_group_id', how='left', validate='one_to_one')
    for field in fields:
        ranked[field] = ranked[f'{field}_magdep']
    ranked['rank_log_r'] = ranked.log_r_per_unit_z.rank(ascending=False, method='min')
    ranked['rank_log_r_magindep'] = ranked.log_r_per_unit_z_magindep.rank(ascending=False, method='min')
    return ranked.sort_values('rank_log_r', na_position='last')


def combine_current_evidence(systems, resolved_all, blend_ranked):
    """Keep all systems; use a valid blend when no eligible resolved result exists."""
    columns = [
        'qso_group_id',
        'primary_qso_catalog',
        'primary_qso_z',
        'decals_resolution_class',
        'scoring_route',
        'route_reason'
    ]
    combined = systems[columns].copy()
    resolved_lookup = resolved_all.set_index('qso_group_id')
    blend_lookup = blend_ranked.set_index('qso_group_id')
    for field in resolved_all:
        if field in columns or field == 'photometric_system':
            continue
        combined[f'resolved_{field}'] = combined.qso_group_id.map(resolved_lookup[field])
    for field in blend_ranked:
        if field in columns:
            continue
        combined[f'blend_{field}'] = combined.qso_group_id.map(blend_lookup[field])
    combined['redshift_scope'] = np.where(
        combined.primary_qso_catalog.isin(['QUAIA_G20_V1', 'QUAIA']),
        'Quaia_Gaussian_redshift_average_exploratory',
        'spectroscopic'
    )
    resolved_any = np.isfinite(combined.resolved_log_r_per_unit_z_magdep)
    resolved_any |= np.isfinite(combined.resolved_log_r_per_unit_z_magindep)
    blend_any = np.isfinite(combined.blend_log_bf_qq_qs_magdep)
    blend_any |= np.isfinite(combined.blend_log_bf_qq_qs_magindep)
    combined['evidence_route'] = np.select([resolved_any, blend_any], ['resolved', 'blend'], default='unscored')
    combined['evidence_parameter'] = combined.evidence_route.map({'resolved': 'log_r_per_unit_z', 'blend': 'log_bf_qq_qs', 'unscored': 'none'})
    for model_name in MODEL_POINTERS:
        combined[f'evidence_log_score_{model_name}'] = np.select(
            [combined.evidence_route.eq('resolved'), combined.evidence_route.eq('blend')],
            [
                combined[f'resolved_log_r_per_unit_z_{model_name}'],
                combined[f'blend_log_bf_qq_qs_{model_name}']
            ],
            default=np.nan
        )
        combined[f'rank_within_route_and_redshift_scope_{model_name}'] = combined.groupby(['evidence_route', 'redshift_scope'])[f'evidence_log_score_{model_name}'].rank(ascending=False, method='min')
    combined['has_numeric_evidence'] = combined.evidence_route.ne('unscored')
    combined['evidence_log_score'] = combined.evidence_log_score_magdep
    combined['rank_within_route_and_redshift_scope'] = combined.rank_within_route_and_redshift_scope_magdep
    reason = combined.resolved_decision_reason_magdep.fillna('')
    blend_reason = combined.blend_status_magdep.fillna('')
    combined['not_ranked_reason'] = reason.where(reason.ne(''), blend_reason)
    combined.loc[
        combined.not_ranked_reason.eq('') & ~combined.has_numeric_evidence,
        'not_ranked_reason'
    ] = combined.route_reason
    combined.loc[combined.has_numeric_evidence, 'not_ranked_reason'] = ''
    return combined

def build_redshift_lookup(systems: pd.DataFrame, original_parent: pd.DataFrame, quaia_path: Path) -> pd.DataFrame:
    """Keep the parent sample fixed and replace Quaia z/zerr by exact Gaia ID."""
    lookup = systems[[
        'qso_group_id',
        'primary_qso_catalog',
        'primary_qso_z',
        'primary_component_gaia_source_id'
    ]].copy()
    parent = original_parent.set_index('qso_group_id')
    lookup['primary_qso_z_original'] = lookup.primary_qso_z
    lookup['primary_qso_zerr_original'] = lookup.qso_group_id.map(parent.primary_qso_zerr)
    lookup['primary_qso_zerr'] = lookup.primary_qso_zerr_original
    lookup['redshift_catalogue_release'] = lookup.primary_qso_catalog.map({
        'SDSS_DR16Q': 'SDSS DR16Q, original observation',
        'DESI_DR1_QSO': 'DESI DR1, original observation',
        'QUAIA_G20_V1': QUAIA_RELEASE
    })
    lookup['redshift_source_doi'] = ''
    lookup['redshift_source_file_md5'] = ''
    quaia_mask = lookup.primary_qso_catalog.eq('QUAIA_G20_V1')
    source_ids = lookup.loc[quaia_mask, 'primary_component_gaia_source_id']
    with fits.open(quaia_path, memmap=True) as hdul:
        table = hdul[1].data
        selected = np.isin(table['source_id'], source_ids.astype('int64'))
        reference = pd.DataFrame({
            'source_id': table['source_id'][selected].astype(str),
            'redshift': table['redshift_quaia'][selected].astype(float),
            'redshift_error': table['redshift_quaia_err'][selected].astype(float)
        }).set_index('source_id')
    lookup.loc[quaia_mask, 'primary_qso_z'] = source_ids.map(reference.redshift)
    lookup.loc[quaia_mask, 'primary_qso_zerr'] = source_ids.map(reference.redshift_error)
    quaia_errors = lookup.loc[quaia_mask, 'primary_qso_zerr']
    if not np.isfinite(quaia_errors).all() or quaia_errors.lt(0).any():
        raise ValueError('Every adopted Quaia redshift requires a published finite error')
    lookup.loc[quaia_mask, 'redshift_source_doi'] = QUAIA_DOI
    lookup.loc[quaia_mask, 'redshift_source_file_md5'] = QUAIA_G20_MD5
    lookup['redshift_changed'] = ~np.isclose(lookup.primary_qso_z_original, lookup.primary_qso_z, rtol=0, atol=1e-12)
    lookup['redshift_uncertainty_changed'] = ~np.isclose(
        lookup.primary_qso_zerr_original,
        lookup.primary_qso_zerr,
        rtol=0,
        atol=1e-12,
        equal_nan=True
    )
    return lookup


def apply_redshift_lookup(frame: pd.DataFrame, lookup: pd.DataFrame) -> pd.DataFrame:
    """Update every scoring view by system ID, preserving row order and inputs."""
    reference = lookup.set_index('qso_group_id')
    expected_original = frame.qso_group_id.map(reference.primary_qso_z_original)
    expected_adopted = frame.qso_group_id.map(reference.primary_qso_z)
    consistent = np.isclose(frame.primary_qso_z, expected_original, rtol=0, atol=1e-12) | np.isclose(frame.primary_qso_z, expected_adopted, rtol=0, atol=1e-12)
    corrected = frame.copy()
    columns = [name for name in lookup.columns if name not in ('qso_group_id', 'primary_qso_catalog', 'primary_component_gaia_source_id')]
    for name in columns:
        corrected[name] = corrected.qso_group_id.map(reference[name])
    return corrected
