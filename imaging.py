import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
MODEL_POINTERS = {'magdep': '20261005_magdep_qlf', 'magindep': '20261005_magindep_qlf'}
SCORE_FIELDS = ['eligible', 'status', 'decision_reason', 'log_r_per_unit_z', 'p_sameq',
                'p_quasar', 'quality_flags', 'qso_support_percentile']


def attach_quaia_redshift_errors(candidates, reference):
    """Join the published uncertainty by the source supplying the adopted z."""
    candidates = candidates.copy()
    lookup = reference.copy()
    lookup['source_id'] = lookup.source_id.astype('string')
    if not lookup.source_id.is_unique:
        raise ValueError('Quaia reference source_id must be unique')
    lookup = lookup.set_index('source_id')
    uses_quaia = candidates.adopted_redshift_catalogue.str.contains('quaia', case=False, na=False)
    source_ids = candidates.source_id.astype('string')
    reference_z = source_ids.map(lookup.redshift_quaia)
    errors = source_ids.map(lookup.redshift_quaia_err)
    invalid = uses_quaia & (~np.isfinite(errors) | errors.lt(0)
                           | ~np.isclose(candidates.adopted_redshift, reference_z, rtol=0, atol=1e-12))
    if invalid.any():
        raise ValueError('Missing or inconsistent Quaia redshift errors for parent indices: '
                         + str(candidates.loc[invalid, 'qso_parent_index'].tolist()))
    candidates['adopted_redshift_error'] = errors.where(uses_quaia, 0.0)
    candidates['redshift_error_source_id'] = source_ids.where(uses_quaia)
    return candidates

def prepare_component_photometry(candidates, components, tractor):
    """Map fixed roles by coordinates, then recover errors by exact Tractor identity."""
    selected = components.loc[components.qso_parent_index.isin(candidates.qso_parent_index)].copy()
    counts = selected.groupby('qso_parent_index').size()
    metadata = [
        'qso_parent_index',
        'physical_system_key',
        'candidate_tier',
        'primary_ra',
        'primary_dec',
        'secondary_ra',
        'secondary_dec',
        'adopted_redshift',
        'adopted_redshift_catalogue',
        'adopted_redshift_error',
        'redshift_error_source_id'
    ]
    selected = selected.merge(candidates[metadata], on='qso_parent_index', validate='many_to_one')
    component_coordinates = SkyCoord(selected.component_ra.to_numpy() * u.deg, selected.component_dec.to_numpy() * u.deg)
    role_matches = {}
    for role in ('primary', 'secondary'):
        coordinates = SkyCoord(
            selected[f'{role}_ra'].to_numpy() * u.deg,
            selected[f'{role}_dec'].to_numpy() * u.deg
        )
        role_matches[role] = component_coordinates.separation(coordinates).arcsec < 0.0001
    selected = selected.rename(columns={'component_role': 'original_r_ranked_role'})
    selected['component_role'] = np.where(role_matches['primary'], 'primary', 'secondary')
    selected['candidate_id'] = selected.physical_system_key + ':' + selected.component_role
    selected['redshift_scope'] = np.where(
        selected.adopted_redshift_catalogue.str.lower().str.contains('quaia'),
        'Quaia_Gaussian_redshift_average_exploratory',
        'spectroscopic'
    )
    identity = ['release', 'brickid', 'objid']
    fields = ['type', 'ra', 'dec', 'maskbits', 'fracflux_r']
    fields += [f'{field}_{band}' for band in 'grz' for field in ('flux', 'flux_ivar', 'nobs')]
    local = tractor[identity + fields].copy()
    for column in identity:
        local[column] = pd.to_numeric(local[column], errors='raise').astype('int64')
        selected[f'tractor_{column}'] = pd.to_numeric(selected[f'tractor_{column}'], errors='raise').astype('int64')
    local = local.drop_duplicates()
    original_fluxes = selected[[f'tractor_flux_{band}' for band in 'grz']].copy()
    selected = selected.drop(columns=original_fluxes.columns)
    local = local.rename(columns={name: f'tractor_{name}' for name in local})
    selected = selected.merge(
        local,
        on=[f'tractor_{name}' for name in identity],
        how='left',
        validate='many_to_one',
        indicator=True
    )
    selected = selected.drop(columns='_merge')
    selected['photometric_system'] = selected.tractor_release.map({9010: 'decals_dr9_south', 9012: 'decals_dr9_south', 9011: 'decals_dr9_north'})
    for band in 'grz':
        flux, ivar = (selected[f'tractor_flux_{band}'], selected[f'tractor_flux_ivar_{band}'])
        usable = np.isfinite(flux) & np.isfinite(ivar) & ivar.gt(0)
        usable &= selected[f'tractor_nobs_{band}'].gt(0) & selected.tractor_maskbits.eq(0)
        selected[f'usable_{band}'] = usable
        selected[f'flux_{band}'] = flux.where(usable, np.nan)
        selected[f'variance_{band}'] = np.divide(
            1.0,
            ivar.to_numpy(float),
            out=np.full(len(selected), np.inf),
            where=usable.to_numpy()
        )
    primary = SkyCoord(selected.primary_ra.to_numpy() * u.deg, selected.primary_dec.to_numpy() * u.deg)
    secondary = SkyCoord(selected.secondary_ra.to_numpy() * u.deg, selected.secondary_dec.to_numpy() * u.deg)
    selected['pair_separation_arcsec'] = primary.separation(secondary).arcsec
    selected['galactic_abs_b_deg'] = np.abs(component_coordinates.galactic.b.deg)
    selected['r_snr'] = selected.flux_r / np.sqrt(selected.variance_r)
    selected['usable_band_count'] = selected[[f'usable_{band}' for band in 'grz']].sum(axis=1)
    return selected


def join_system_scores(candidates, scores):
    systems = candidates.copy()
    for model in MODEL_POINTERS:
        for role in ('primary', 'secondary'):
            subset = scores.loc[
                scores.model.eq(model) & scores.component_role.eq(role),
                ['qso_parent_index', *SCORE_FIELDS]
            ]
            subset = subset.rename(columns={field: f'{role}_{field}_{model}' for field in SCORE_FIELDS})
            systems = systems.merge(subset, on='qso_parent_index', how='left', validate='one_to_one')
        eligible = systems[f'primary_eligible_{model}'].fillna(False).astype(bool)
        eligible &= systems[f'secondary_eligible_{model}'].fillna(False).astype(bool)
        systems[f'both_components_eligible_{model}'] = eligible
        systems[f'minimum_component_log_r_{model}'] = np.minimum(
            systems[f'primary_log_r_per_unit_z_{model}'],
            systems[f'secondary_log_r_per_unit_z_{model}']
        ).where(eligible)
        systems[f'rank_{model}'] = systems[f'minimum_component_log_r_{model}'].rank(ascending=False, method='min')
    return systems


