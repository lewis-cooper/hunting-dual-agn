import numpy as np
import pandas as pd
from dual_agn.redshift_inputs import photometric_redshifts, redshift_errors
from astropy import units as u
from astropy.coordinates import SkyCoord
MODEL_POINTERS = {'magdep': '20261005_magdep_qlf', 'magindep': '20261005_magindep_qlf'}
SCORE_FIELDS = ['eligible', 'status', 'decision_reason', 'log_r_per_unit_z', 'p_sameq',
                'p_quasar', 'quality_flags', 'qso_support_percentile']
ROLE_MATCH_TOLERANCE_ARCSEC = 0.0001


def validate_component_counts(candidates, components):
    """Require two prepared rows for every parent, without dropping missing systems."""
    if candidates.qso_parent_index.isna().any() or not candidates.qso_parent_index.is_unique:
        raise ValueError('Imaging parent indices must be unique and non-missing')
    if components.qso_parent_index.isna().any():
        raise ValueError('Imaging component parent indices must be non-missing')
    counts = components.groupby('qso_parent_index').size().reindex(
        candidates.qso_parent_index, fill_value=0)
    invalid = counts.ne(2)
    unexpected = set(components.qso_parent_index) - set(candidates.qso_parent_index)
    if invalid.any() or unexpected:
        raise ValueError('Exactly two imaging components required for parent indices: '
                         + str(counts.index[invalid].tolist() + sorted(unexpected)))


def validate_component_redshifts(candidates, components):
    """Compare adopted system metadata, not independently measured component redshifts."""
    parents = candidates.set_index('qso_parent_index')
    for field in ('adopted_redshift', 'adopted_redshift_error'):
        values = pd.to_numeric(parents[field], errors='coerce')
        invalid_parent = ~np.isfinite(values) | values.lt(0)
        if invalid_parent.any():
            raise ValueError(f'Invalid {field} for parent indices: '
                             + str(parents.index[invalid_parent].tolist()))
        # Raw discovery components use z; prepared components use adopted_redshift.
        component_field = 'z' if field == 'adopted_redshift' and field not in components else field
        if component_field not in components:
            continue
        expected = components.qso_parent_index.map(parents[field])
        actual = pd.to_numeric(components[component_field], errors='coerce')
        inconsistent = ~np.isfinite(actual) | ~np.isclose(actual, expected, rtol=0, atol=1e-12)
        if inconsistent.any():
            raise ValueError(f'Inconsistent {field} for imaging parent indices: '
                             + str(components.loc[inconsistent, 'qso_parent_index'].unique().tolist()))
    for field in ('adopted_redshift_catalogue', 'redshift_kind'):
        if field not in candidates or field not in components:
            continue
        expected = components.qso_parent_index.map(parents[field]).astype('string').str.lower()
        actual = components[field].astype('string').str.lower()
        inconsistent = actual.isna() | expected.isna() | actual.ne(expected)
        if inconsistent.any():
            raise ValueError(f'Inconsistent {field} for imaging parent indices: '
                             + str(components.loc[inconsistent, 'qso_parent_index'].unique().tolist()))


def validate_component_photometry(candidates, components):
    """Check identities, fixed coordinate roles and adopted redshifts before scoring."""
    validate_component_counts(candidates, components)
    validate_component_redshifts(candidates, components)
    roles = components.groupby('qso_parent_index').component_role.agg(list)
    invalid_roles = roles.map(lambda values: set(values) != {'primary', 'secondary'})
    if invalid_roles.any():
        raise ValueError('One unique primary and secondary required for parent indices: '
                         + str(roles.index[invalid_roles].tolist()))
    if components.candidate_id.isna().any() or not components.candidate_id.is_unique:
        raise ValueError('Imaging candidate_id must be unique and non-missing')
    identity = ['tractor_release', 'tractor_brickid', 'tractor_objid']
    invalid_identity = components[identity].isna().any(axis=1)
    duplicate_identity = components.duplicated(['qso_parent_index', *identity], keep=False)
    if (invalid_identity | duplicate_identity).any():
        raise ValueError('Distinct, non-missing Tractor identities required for parent indices: '
                         + str(components.loc[invalid_identity | duplicate_identity,
                                               'qso_parent_index'].unique().tolist()))
    parents = candidates.set_index('qso_parent_index')
    measured = SkyCoord(components.component_ra.to_numpy(float)*u.deg,
                        components.component_dec.to_numpy(float)*u.deg)
    matches = {}
    for role in ('primary', 'secondary'):
        expected = SkyCoord(
            components.qso_parent_index.map(parents[f'{role}_ra']).to_numpy(float)*u.deg,
            components.qso_parent_index.map(parents[f'{role}_dec']).to_numpy(float)*u.deg)
        matches[role] = measured.separation(expected).arcsec < ROLE_MATCH_TOLERANCE_ARCSEC
    unique_match = matches['primary'] ^ matches['secondary']
    correct_role = np.where(components.component_role.eq('primary'), matches['primary'], matches['secondary'])
    invalid_position = ~unique_match | ~correct_role
    if invalid_position.any():
        raise ValueError('Components must uniquely match their adopted coordinate roles; parent indices: '
                         + str(components.loc[invalid_position, 'qso_parent_index'].unique().tolist()))


def imaging_pair_quality(components):
    """Derive both-component quality cuts from measurements, not cached flags."""
    counts = components.groupby('qso_parent_index').size()
    if counts.ne(2).any():
        raise ValueError('Pair quality requires exactly two component rows per parent')
    snr = pd.to_numeric(components.r_snr, errors='coerce')
    latitude = pd.to_numeric(components.galactic_abs_b_deg, errors='coerce')
    flags = pd.DataFrame({
        'reference_snr_pass': np.isfinite(snr) & snr.ge(10),
        'latitude_pass': np.isfinite(latitude) & latitude.ge(25),
        'psf_pass': components.tractor_type.eq('PSF').fillna(False),
    }, index=components.index)
    return flags.groupby(components.qso_parent_index).all()


def attach_quaia_redshift_errors(candidates, reference):
    """Join the published uncertainty by the source supplying the adopted z."""
    candidates = candidates.copy()
    if 'redshift_kind' in candidates:
        candidates['adopted_redshift_error'] = redshift_errors(candidates, 'adopted_redshift_error', 'adopted_redshift_catalogue', 'redshift_kind')
        candidates['redshift_error_source_id'] = candidates.source_id.astype('string').where(candidates.redshift_kind.eq('photometric'))
        return candidates
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
    validate_component_counts(candidates, selected)
    validate_component_redshifts(candidates, selected)
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
    if 'redshift_kind' in candidates:
        metadata.append('redshift_kind')
    selected = selected.drop(columns=[field for field in metadata
                                      if field != 'qso_parent_index' and field in selected])
    selected = selected.merge(candidates[metadata], on='qso_parent_index', validate='many_to_one')
    component_coordinates = SkyCoord(selected.component_ra.to_numpy() * u.deg, selected.component_dec.to_numpy() * u.deg)
    role_matches = {}
    for role in ('primary', 'secondary'):
        coordinates = SkyCoord(
            selected[f'{role}_ra'].to_numpy() * u.deg,
            selected[f'{role}_dec'].to_numpy() * u.deg
        )
        role_matches[role] = component_coordinates.separation(coordinates).arcsec < ROLE_MATCH_TOLERANCE_ARCSEC
    unique_role_match = role_matches['primary'] ^ role_matches['secondary']
    if not unique_role_match.all():
        raise ValueError('Components do not uniquely match adopted primary/secondary coordinates; parent indices: '
                         + str(selected.loc[~unique_role_match, 'qso_parent_index'].unique().tolist()))
    selected = selected.rename(columns={'component_role': 'original_r_ranked_role'})
    selected['component_role'] = np.where(role_matches['primary'], 'primary', 'secondary')
    selected['candidate_id'] = selected.physical_system_key + ':' + selected.component_role
    selected['redshift_scope'] = np.where(
        photometric_redshifts(selected, 'adopted_redshift_catalogue', 'redshift_kind'),
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
    if selected._merge.ne('both').any():
        raise ValueError('Missing exact Tractor measurements for imaging parent indices: '
                         + str(selected.loc[selected._merge.ne('both'), 'qso_parent_index'].unique().tolist()))
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
    validate_component_photometry(candidates, selected)
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


