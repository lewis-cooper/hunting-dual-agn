from pathlib import Path
import json
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.table import Table
from settings import *
from photometry import read_fits
from one_arcsecond import build_redshift_lookup, apply_redshift_lookup, attach_tractor_diagnostics
ONE_ROOT.mkdir(parents=True, exist_ok=True)
(ONE_ROOT / 'inputs').mkdir(exist_ok=True)
systems = pd.read_csv(INPUT_ROOT / 'eligibility_audit_712.csv', dtype=IDENTIFIER_DTYPES)
candidates = pd.read_csv(
    INPUT_ROOT / 'secondary_measurements_324.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)
original_parent = pd.read_csv(INPUT_ROOT / 'original_parent_712.csv', dtype=IDENTIFIER_DTYPES)
lookup = build_redshift_lookup(systems, original_parent, INPUT_ROOT / 'quaia_712_reference.fits')
systems = apply_redshift_lookup(systems, lookup)
candidates = apply_redshift_lookup(candidates, lookup)
lookup.to_csv(ONE_ROOT / 'inputs/adopted_redshifts_712.csv', index=False)
systems.to_csv(ONE_ROOT / 'inputs/eligibility_audit_712.csv', index=False)
candidates.to_csv(ONE_ROOT / 'inputs/secondary_measurements_324.csv', index=False)
# Reproduce the notebook's corrected-input CSV boundary before any calculations.
# The external support sampler hashes exact floating-point inputs.
systems = pd.read_csv(ONE_ROOT / 'inputs/eligibility_audit_712.csv', dtype=IDENTIFIER_DTYPES)
candidates = pd.read_csv(ONE_ROOT / 'inputs/secondary_measurements_324.csv', dtype=IDENTIFIER_DTYPES)
resolution_path = INPUT_ROOT / 'one_arcsecond_resolution.csv'
tractor_path = INPUT_ROOT / 'one_arcsecond_tractor.csv'
resolution = pd.read_csv(resolution_path, dtype={'qso_group_id': str})
systems = systems.merge(
    resolution[['qso_group_id', 'decals_resolution_class']],
    on='qso_group_id',
    validate='one_to_one'
)
desi_target_photometry = pd.read_csv(INPUT_ROOT / 'desi_target_photometry.csv', dtype=IDENTIFIER_DTYPES)
systems = systems.merge(desi_target_photometry, on='qso_group_id', how='left', validate='one_to_one')
primary_only_mask = systems.decals_resolution_class.eq('decals_unresolved_primary_only')
nonresolved_class_mask = ~systems.decals_resolution_class.eq('decals_resolved')
resolved_mask = systems.qso_group_id.isin(candidates.qso_group_id)
measured_decals_r = {}
for role in ('primary', 'secondary'):
    flux_r = pd.to_numeric(systems[f'{role}_decals_flux_r'], errors='coerce')
    ivar_r = pd.to_numeric(systems[f'{role}_decals_flux_ivar_r'], errors='coerce')
    measured_decals_r[role] = systems[f'{role}_decals_release'].isin((9010, 9011, 9012)) & np.isfinite(flux_r) & np.isfinite(ivar_r) & ivar_r.gt(0)
component_r_available = measured_decals_r['primary'] | measured_decals_r['secondary']
desi_r_flux = pd.to_numeric(systems.desi_flux_r, errors='coerce')
desi_r_ivar = pd.to_numeric(systems.desi_flux_ivar_r, errors='coerce')
systems['desi_target_r_available'] = systems.primary_qso_catalog.eq('DESI_DR1_QSO') & systems.desi_release.isin((9010, 9011, 9012)) & np.isfinite(desi_r_flux) & np.isfinite(desi_r_ivar) & desi_r_ivar.gt(0)
desi_fallback_mask = nonresolved_class_mask & ~component_r_available & systems.desi_target_r_available
legacy_r_available_mask = component_r_available | systems.desi_target_r_available
blend_input_mask = nonresolved_class_mask
two_source_mask = systems.decals_resolution_class.eq('decals_resolved') | resolved_mask
secondary_flux = systems[[f'secondary_decals_flux_{band}' for band in 'grz']].to_numpy(float)
secondary_ivar = systems[[f'secondary_decals_flux_ivar_{band}' for band in 'grz']].to_numpy(float)
three_measured_secondary_bands = (np.isfinite(secondary_flux) & np.isfinite(secondary_ivar) & (secondary_ivar > 0)).all(axis=1)
supplemental_mask = two_source_mask & ~resolved_mask & three_measured_secondary_bands
source_columns = [
    'release',
    'brickid',
    'objid',
    'ra',
    'dec',
    'maskbits',
    'nobs_g',
    'nobs_r',
    'nobs_z',
    'type',
    'fracflux_r'
]
sources = pd.read_csv(tractor_path, usecols=source_columns)
sources = sources.drop_duplicates(['release', 'brickid', 'objid'])
primary_sky = SkyCoord(
    systems.primary_component_ra.to_numpy(float) * u.deg,
    systems.primary_component_dec.to_numpy(float) * u.deg
)
secondary_sky = SkyCoord(
    systems.secondary_component_ra.to_numpy(float) * u.deg,
    systems.secondary_component_dec.to_numpy(float) * u.deg
)
source_sky = SkyCoord(sources.ra.to_numpy(float) * u.deg, sources.dec.to_numpy(float) * u.deg)
primary_index, primary_source_index, _, _ = search_around_sky(primary_sky, source_sky, 0.5 * u.arcsec)
secondary_index, secondary_source_index, _, _ = search_around_sky(secondary_sky, source_sky, 0.5 * u.arcsec)
near_primary = pd.DataFrame({'system_index': primary_index, 'source_index': primary_source_index})
near_secondary = pd.DataFrame({'system_index': secondary_index, 'source_index': secondary_source_index})
near_either = pd.concat([near_primary, near_secondary]).drop_duplicates()
near_both = near_primary.merge(near_secondary, on=['system_index', 'source_index'])
systems['n_sources_near_either'] = systems.index.to_series().map(near_either.groupby('system_index').size()).fillna(0).astype(int)
systems['n_sources_near_both'] = systems.index.to_series().map(near_both.groupby('system_index').size()).fillna(0).astype(int)
systems['one_local_source_near_both'] = systems.n_sources_near_either.eq(1) & systems.n_sources_near_both.eq(1)
from astropy.table import Table
component_metadata_path = INPUT_ROOT / 'one_arcsecond_components.fits'
component_metadata = Table.read(component_metadata_path)[['decals_release', 'decals_brickid', 'decals_objid', 'decals_type']].to_pandas()
component_metadata['decals_type'] = component_metadata.decals_type.map(lambda value: value.decode().strip() if isinstance(value, bytes) else str(value).strip())
systems = attach_tractor_diagnostics(systems, component_metadata, sources)
systems['component_separation_arcsec'] = primary_sky.separation(secondary_sky).arcsec
new_columns = [
    'qso_group_id',
    'primary_decals_type',
    'secondary_decals_type',
    'primary_decals_fracflux_r',
    'secondary_decals_fracflux_r',
    'component_separation_arcsec'
]
candidates = candidates.merge(systems[new_columns], on='qso_group_id', validate='one_to_one')
print('Measured component separations:', systems.component_separation_arcsec.agg(['min', 'max']).to_dict())
print('Secondary morphology missing:', int(systems.secondary_decals_type.isna().sum()))
unresolved_mask = blend_input_mask & ~resolved_mask
systems['scoring_route'] = np.select(
    [resolved_mask, supplemental_mask, unresolved_mask],
    ['resolved_original', 'resolved_supplemental', 'unresolved'],
    default='currently_unscoreable'
)
systems['route_reason'] = np.select(
    [
        resolved_mask,
        supplemental_mask,
        two_source_mask,
        desi_fallback_mask,
        unresolved_mask & legacy_r_available_mask,
        unresolved_mask
    ],
    [
        'two_components_passed_original_quality_rules',
        'three_measured_secondary_grz_bands_outside_original_audit',
        'two_sources_but_component_quality_failed',
        'desi_target_Legacy_DR9_photometry_used_as_system_flux',
        'primary_or_secondary_DR9_assignment_used_as_system_flux',
        'no_legacy_r_try_survey_fallback'
    ],
    default='two_sources_but_component_quality_failed'
)
systems.to_csv(ONE_ROOT / 'scoring_input_audit_712.csv', index=False)
candidates.to_csv(ONE_ROOT / 'original_resolved_measurements.csv', index=False)
print('Input:', INPUT_ROOT, 'Output:', ONE_ROOT)
print(systems.scoring_route.value_counts().to_string())
