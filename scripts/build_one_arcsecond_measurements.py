"""Generate component associations, resolution and photometry from native catalogues."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from imaging_search.catalogues import read_chunks, write_fits, file_hash
from measurement_inputs import (catalogue_files, read_config, local_path, read_local_sources,
                                require_coverage)
from one_arcsecond_measurements import (build_components, build_system_measurements,
                                       build_survey_photometry, TRACTOR_FIELDS, desi_target_ids,
                                       build_desi_fallback)
from settings import (ROOT, ONE_ARCSECOND_PARENT_PATH, ONE_INPUT_ROOT, survey_products,
                      IDENTIFIER_DTYPES)


def build_measurements(config, parent, output, parent_path=None):
    output.mkdir(parents=True, exist_ok=True)
    quasar_positions = pd.DataFrame({'source_id': parent.matched_gaia_source_id,
                                     'ra': parent.primary_qso_ra, 'dec': parent.primary_qso_dec})
    # Anchors can be offset by up to 1 arcsec from the quasar coordinates.
    gaia = read_local_sources(quasar_positions, config, config['gaia'], 2., ['source_id'],
        ['pmra', 'pmdec', 'pmra_error', 'pmdec_error', 'pmra_pmdec_corr', 'phot_g_mean_mag'])
    anchors = gaia.set_index('source_id').reindex(parent.matched_gaia_source_id).rename_axis('source_id').reset_index()
    if anchors.ra.isna().any():
        raise ValueError('Native Gaia data lack discovered anchor IDs; check catalogue release and coverage')
    for name, radius in [('gaia', 1.), ('tractor', 1.5)]:
        require_coverage(anchors, config, config[name], radius)
    tractor = read_local_sources(anchors, config, config['tractor'], 1.5, ['release', 'brickid', 'objid'], TRACTOR_FIELDS)
    # No preference for an alternative Tractor object if the nearest is shared.
    components = build_components(parent, gaia, tractor)
    systems = build_system_measurements(parent, components)
    components.to_csv(output / 'components.csv', index=False)
    write_fits(components, output / 'one_arcsecond_components.fits')
    systems.to_csv(output / 'component_quality_audit.csv', index=False)
    systems.loc[systems.eligible_both_components].to_csv(output / 'quality_resolved_measurements.csv', index=False)
    systems[['qso_group_id', 'decals_resolution_class']].to_csv(output / 'one_arcsecond_resolution.csv', index=False)
    tractor.to_csv(output / 'one_arcsecond_tractor.csv', index=False)
    # The adopted source ID, rather than a positional substitute, selects Quaia z/errors.
    quaia_rows = []
    legacy_quaia = parent.primary_qso_catalog.eq('QUAIA_G20_V1')
    if 'primary_qso_redshift_kind' in parent:
        legacy_quaia &= False
    quaia_ids = set(parent.loc[legacy_quaia, 'matched_gaia_source_id'])
    if quaia_ids:
        quaia_path = local_path(config, config['quaia_catalogue'])
        for chunk in read_chunks(quaia_path, ['source_id', 'redshift_quaia', 'redshift_quaia_err']):
            quaia_rows.append(chunk.loc[chunk.source_id.isin(quaia_ids)])
    if not quaia_rows:
        quaia_rows = [pd.DataFrame({'source_id':pd.Series(dtype='string'),
            'redshift_quaia':pd.Series(dtype=float), 'redshift_quaia_err':pd.Series(dtype=float)})]
    write_fits(pd.concat(quaia_rows, ignore_index=True), output / 'quaia_reference.fits')
    desi_ids = {value for values in desi_target_ids(parent).values() for value in values}
    desi_rows = []
    fields = ['targetid', 'release', *[f'{quantity}_{band}' for quantity in ('flux', 'flux_ivar') for band in 'grz']]
    if desi_ids:
        desi_spec = config['desi']
        rename = {native: name for name, native in desi_spec.get('columns', {}).items()}
        native_fields = [desi_spec.get('columns', {}).get(name, name) for name in fields]
        for path in catalogue_files(config, desi_spec):
            for chunk in read_chunks(path, native_fields):
                chunk = chunk.rename(columns=rename)
                if not set(fields).issubset(chunk):
                    raise ValueError('DESI native catalogue requires '+', '.join(fields))
                chunk['targetid'] = chunk.targetid.astype('string')
                desi_rows.append(chunk.loc[chunk.targetid.isin(desi_ids), fields])
    if not desi_rows:
        desi_rows = [pd.DataFrame(columns=fields)]
    desi = pd.concat(desi_rows, ignore_index=True).drop_duplicates()
    if not desi.targetid.is_unique:
        raise ValueError('DESI target IDs have conflicting native photometry')
    fallback = build_desi_fallback(parent, desi)
    fallback.to_csv(output / 'desi_target_photometry.csv', index=False)
    missing_surveys = set(survey_products).difference(config.get('surveys', {}))
    if missing_surveys:
        raise ValueError('Supply native photometry/coverage for surveys: '+', '.join(sorted(missing_surveys)))
    survey_counts = {}
    for name, specification in survey_products.items():
        if config['surveys'][name].get('enabled') is False:
            photometry = parent[['qso_group_id']].copy()
            photometry[f'{name}_source_matched'] = False
            photometry[specification['source_flag']] = False
            for band, (magnitude, error) in specification['bands'].items():
                photometry[magnitude] = np.nan
                photometry[error] = np.nan
                photometry[f'{name}_{band}_band_quality'] = False
            if name == 'sdss':
                photometry['sdss_clean'] = False
                photometry['sdss_mode'] = np.nan
        else:
            photometry = build_survey_photometry(parent, config, name, config['surveys'][name], specification['bands'])
        photometry.to_csv(output / f'{name}_measurements.csv', index=False)
        survey_counts[name] = int(photometry[f'{name}_source_matched'].sum())
    summary = dict(parent_entries=len(parent), gaia_components=len(components),
                   component_count_mismatches=int((~systems.primary_component_count_matches_expected).sum()),
                   original_quality_resolved=int(systems.eligible_both_components.sum()),
                   resolution_counts=systems.decals_resolution_class.value_counts().to_dict(),
                   survey_matches=survey_counts, input_parent_sha256=file_hash(parent_path) if parent_path else None,
                   desi_target_photometry_status=fallback.desi_target_photometry_status.value_counts().to_dict(),
                   output_root=str(output), visual_inspection_applied=False)
    (output / 'measurement_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT / 'measurement_inputs.local.json')
    parser.add_argument('--parent', type=Path, default=ONE_ARCSECOND_PARENT_PATH)
    parser.add_argument('--output-root', type=Path, default=ONE_INPUT_ROOT)
    args = parser.parse_args()
    if not args.config.is_file():
        parser.error('Copy measurement_inputs.example.json and supply native catalogue paths/coverage')
    parent = pd.read_csv(args.parent, dtype={**IDENTIFIER_DTYPES, 'matched_gaia_source_id': 'string', 'primary_qso_objid': 'string'})
    build_measurements(read_config(args.config), parent, args.output_root, args.parent)


if __name__ == '__main__':
    main()
