"""Build Gaia component, Tractor resolution and survey photometry from native rows."""

from collections import Counter
import json
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import search_around_sky

from dual_agn.imaging_search.pairs import pm_sig
from dual_agn.measurement_inputs import (sky, read_local_sources, nearest_sources, require_coverage)

TRACTOR_ID = ['release', 'brickid', 'objid']
TRACTOR_FIELDS = [*TRACTOR_ID, 'type', 'ra', 'dec', 'maskbits', 'fracflux_r',
                  *[f'{quantity}_{band}' for quantity in ('flux', 'flux_ivar', 'nobs') for band in 'grz']]


def desi_target_ids(parent):
    """Retain all native target IDs when a master row combines observations."""
    identifiers = {}
    for row in parent.loc[parent.primary_qso_catalog.eq('DESI_DR1_QSO')].itertuples():
        if hasattr(row, 'qso_catalogue_object_ids_json'):
            values = json.loads(row.qso_catalogue_object_ids_json)
        else:
            values = [str(row.primary_qso_objid)]
        identifiers[row.qso_group_id] = list(dict.fromkeys(str(value) for value in values))
    return identifiers


def build_desi_fallback(parent, photometry):
    """Use exact target photometry; conflicting alternatives stay unavailable and flagged."""
    fields = ['release', *[f'{quantity}_{band}' for quantity in ('flux','flux_ivar') for band in 'grz']]
    lookup = photometry.set_index('targetid')
    targets = desi_target_ids(parent)
    rows = []
    for system in parent.itertuples():
        identifiers = targets.get(system.qso_group_id, [])
        available = lookup.loc[lookup.index.isin(identifiers), fields].drop_duplicates()
        status = ('matched_exact_target' if len(available)==1 else
                  'conflicting_native_targets' if len(available)>1 else 'no_target_photometry')
        values = available.iloc[0] if len(available)==1 else pd.Series(np.nan,index=fields)
        rows.append({'qso_group_id':system.qso_group_id,
            **{'desi_'+field:values[field] for field in fields},
            'desi_targetids_json':json.dumps(identifiers), 'desi_target_photometry_status':status})
    return pd.DataFrame(rows)


def component_quality(component):
    if component['decals_assignment_status'] != 'assigned_unique_nearest_decals_object':
        return False, component['decals_assignment_status']
    if component['decals_type'] != 'PSF':
        return False, 'matched_decals_type_not_psf'
    flux = np.array([component[f'decals_flux_{band}'] for band in 'grz'], float)
    ivar = np.array([component[f'decals_flux_ivar_{band}'] for band in 'grz'], float)
    nobs = np.array([component[f'decals_nobs_{band}'] for band in 'grz'], float)
    for passed, reason in [
        (np.isfinite(flux).all() and (flux > 0).all(), 'nonfinite_or_nonpositive_grz_flux'),
        (np.isfinite(ivar).all() and (ivar > 0).all(), 'nonpositive_grz_inverse_variance'),
        ((nobs > 0).all(), 'no_positive_grz_observations'),
        (np.isfinite(ivar).all() and (flux * np.sqrt(np.maximum(ivar, 0)) >= 3).all(), 'below_minimum_grz_snr')]:
        if not passed:
            return False, reason
    return True, 'passes_resolved_decals_quality'


def build_components(parent, gaia, tractor, component_radius_arcsec=1., tractor_radius_arcsec=.1):
    """Keep anchor first, other Gaia components by distance; shared Tractor withheld."""
    required_gaia = {'source_id', 'ra', 'dec', 'pmra', 'pmdec', 'pmra_error', 'pmdec_error', 'pmra_pmdec_corr', 'phot_g_mean_mag'}
    if not required_gaia.issubset(gaia) or not set(TRACTOR_FIELDS).issubset(tractor):
        raise ValueError('Native Gaia/Tractor inputs are missing required measurement columns')
    gaia = gaia.copy()
    gaia['source_id'] = gaia.source_id.astype('string')
    if not gaia.source_id.str.fullmatch('[1-9][0-9]*').all() or not gaia.source_id.is_unique:
        raise ValueError('Gaia source_id must contain unique exact integer identifiers')
    gaia_lookup = gaia.set_index('source_id')
    anchor = gaia_lookup.reindex(parent.matched_gaia_source_id)
    if anchor.ra.isna().any():
        raise ValueError('Native Gaia data do not contain every discovered anchor')
    target_index, source_index, separation, _ = search_around_sky(sky(anchor), sky(gaia), component_radius_arcsec * u.arcsec)
    grouped = pd.DataFrame({'target': target_index, 'source': source_index, 'separation': separation.arcsec})
    rows = []
    tractor_sky = sky(tractor)
    for target, members in grouped.groupby('target', sort=True):
        system = parent.iloc[target]
        members = members.copy()
        members['source_id'] = gaia.iloc[members.source].source_id.to_numpy()
        members['not_anchor'] = members.source_id.ne(system.matched_gaia_source_id)
        members = members.sort_values(['not_anchor', 'separation', 'source_id'])
        sources = gaia.iloc[members.source].reset_index(drop=True)
        component_index, tractor_index, distance, _ = search_around_sky(sky(sources), tractor_sky, tractor_radius_arcsec * u.arcsec)
        candidates = pd.DataFrame({'component': component_index, 'tractor': tractor_index, 'separation': distance.arcsec})
        for name in TRACTOR_ID:
            candidates[name] = tractor.iloc[tractor_index][name].to_numpy()
        closest = candidates.sort_values(['component', 'separation', *TRACTOR_ID]).drop_duplicates('component').set_index('component')
        nearest_identities = {index: tuple(match[name] for name in TRACTOR_ID) for index, match in closest.iterrows()}
        counts = Counter(nearest_identities.values())
        for index, source in sources.iterrows():
            identity = nearest_identities.get(index)
            assigned = identity is not None and counts[identity] == 1
            record = dict(qso_group_id=system.qso_group_id, component_index=index + 1,
                          component_role='primary' if index == 0 else ('secondary' if index == 1 else f'companion_{index}'),
                          component_gaia_source_id=source.source_id, component_ra=source.ra, component_dec=source.dec,
                          component_anchor_sep_arcsec=float(members.separation.iloc[index]),
                          gaia_phot_g_mean_mag=source.get('phot_g_mean_mag', np.nan),
                          expected_gaia_component_count=int(system.gaia_edr3_cnt1), actual_component_count=len(sources),
                          component_count_matches_expected=len(sources) == int(system.gaia_edr3_cnt1),
                          decals_assignment_status=('assigned_unique_nearest_decals_object' if assigned else
                              'shared_unresolved_nearest_decals_object' if identity is not None else 'no_decals_object_within_radius'))
            for name in ('pmra', 'pmdec', 'pmra_error', 'pmdec_error', 'pmra_pmdec_corr'):
                record[f'component_gaia_{name}'] = source.get(name, np.nan)
            for name in TRACTOR_FIELDS:
                record['decals_' + name] = tractor.iloc[int(closest.loc[index, 'tractor'])][name] if assigned else ('' if name == 'type' else np.nan)
            quality, reason = component_quality(record)
            record['decals_component_photometry_quality'] = quality
            record['decals_component_quality_status'] = reason
            rows.append(record)
    return pd.DataFrame(rows)


def build_system_measurements(parent, components):
    """Derive wide measurements and structural classes without colour-locus cuts."""
    systems = parent.copy()
    roles = components.loc[components.component_role.isin(['primary', 'secondary'])]
    for role in ('primary', 'secondary'):
        selected = roles.loc[roles.component_role.eq(role)].set_index('qso_group_id')
        if len(selected) != len(parent):
            raise ValueError(f'Every discovered entry needs a {role} Gaia component; check input releases/coverage')
        for name in selected:
            if name in ('component_role', 'component_index'):
                continue
            systems[f'{role}_{name}'] = systems.qso_group_id.map(selected[name])
    secondary = systems[[f'secondary_component_gaia_{name}' for name in
                         ('pmra', 'pmdec', 'pmra_error', 'pmdec_error', 'pmra_pmdec_corr')]].to_numpy(float)
    systems['secondary_pm_significance'] = pm_sig(*secondary.T)
    primary_quality = systems.primary_decals_component_photometry_quality
    secondary_quality = systems.secondary_decals_component_photometry_quality
    systems['eligible_both_components'] = primary_quality & secondary_quality
    systems['eligibility_reason'] = np.select([~primary_quality, ~secondary_quality],
        ['primary_quality_failed', 'secondary_quality_failed'], default='eligible')
    assigned = 'assigned_unique_nearest_decals_object'
    primary_assigned = systems.primary_decals_assignment_status.eq(assigned)
    secondary_assigned = systems.secondary_decals_assignment_status.eq(assigned)
    exactly_two = systems.gaia_edr3_cnt1.eq(2) & systems.primary_actual_component_count.eq(2)
    systems['decals_resolution_class'] = np.select([
        exactly_two & primary_assigned & secondary_assigned,
        exactly_two & primary_assigned & ~secondary_assigned],
        ['decals_resolved', 'decals_unresolved_primary_only'], default='unclassified')
    return systems


def build_survey_photometry(parent, config, survey, specification, bands):
    """Match native survey positions, then compute the original per-band flags."""
    targets = pd.DataFrame({'source_id': parent.matched_gaia_source_id,
                            'ra': parent.primary_qso_ra, 'dec': parent.primary_qso_dec})
    require_coverage(targets, config, specification, 1.)
    required = ['source_id', *[name for pair in bands.values() for name in pair]]
    identity = specification.get('record_identity', ['source_id'])
    if not isinstance(identity, list) or not identity or not all(isinstance(name, str) for name in identity):
        raise ValueError(f'{survey}: record_identity must be a nonempty list of canonical native-row columns')
    sources = read_local_sources(targets, config, specification, 1., identity, required)
    if not set(required).issubset(sources):
        raise ValueError(f'{survey} native column map needs {required}')
    matched = nearest_sources(targets, sources, 1., identity)
    result = pd.DataFrame({'qso_group_id': parent.qso_group_id})
    present = matched.source_id.notna()
    result[f'{survey}_source_id'] = matched.source_id.to_numpy()
    result[f'{survey}_source_ra_deg'] = matched.ra.to_numpy()
    result[f'{survey}_source_dec_deg'] = matched.dec.to_numpy()
    result[f'{survey}_match_sep_arcsec'] = matched.match_sep_arcsec.to_numpy()
    result[f'{survey}_source_matched'] = present.to_numpy()
    for band, (magnitude_column, error_column) in bands.items():
        magnitude = pd.to_numeric(matched[magnitude_column], errors='coerce')
        error = pd.to_numeric(matched[error_column], errors='coerce')
        result[magnitude_column] = magnitude.to_numpy()
        result[error_column] = error.to_numpy()
        quality = present & magnitude.gt(0) & magnitude.lt(40) & error.ge(0) & error.le(.3)
        result[f'{survey}_{band}_band_quality'] = quality.to_numpy()
    return result
