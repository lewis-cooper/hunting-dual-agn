"""Export raw WSDB records, then hand local configurations to the shared workflow.

The master catalogue, model assets and FITS images remain user inputs. Sources
omitted from the WSDB configuration use the supplied local input configuration.
"""

import argparse
import json
from pathlib import Path

import pandas as pd

from imaging_search.catalogues import read_frame, read_master, resolve_config as imaging_config
from imaging_search.pipeline import discover
from measurement_inputs import local_path, read_config, sky
from one_arcsecond_search import resolve_config as one_config, run_search as search_one
from one_arcsecond_measurements import TRACTOR_FIELDS, desi_target_ids
from scripts.prepare_quasar_catalogue import adapt_catalogue, read_quasars
from settings import ROOT, WORK_ROOT, survey_products
from wsdb_inputs import WsdbCatalogues, WsdbInputError, source_tables
from catalogue_downloads import export_public_sources


def absolute_paths(config):
    """Keep relative local inputs valid when writing generated configs elsewhere."""
    path_keys = {'path', 'paths', 'coverage', 'path_template', 'catalogue', 'gaia_catalogue',
                 'gaia_coverage', 'master_catalogue', 'tractor_catalogues', 'tractor_coverage',
                 'ccd_catalogue', 'ccd_coverage', 'quaia_catalogue'}

    def convert(value):
        if not isinstance(value, dict):
            return value
        result = {}
        for key, item in value.items():
            if key == '_directory':
                continue
            if key in path_keys:
                result[key] = [str(local_path(config, path)) for path in item] if isinstance(item, list) else str(local_path(config, item))
            else:
                result[key] = convert(item)
        return result

    return convert(config)


def write_config(path, config):
    path.write_text(json.dumps(config, indent=2, default=str) + '\n')
    return path


def anchor_coverage(anchors, quasars, radius_arcsec=2.):
    """Record the fully enclosed anchor cones within successful quasar cones."""
    _, offsets, _ = sky(anchors).match_to_catalog_sky(sky(quasars))
    available = radius_arcsec - offsets.arcsec
    if (available < 1. - 1e-6).any():
        raise ValueError('Some discovered anchors lack a complete 1-arcsec Gaia neighbourhood')
    return pd.DataFrame({'source_id': anchors.source_id, 'ra_deg': anchors.ra,
        'dec_deg': anchors.dec, 'query_radius_arcsec': available,
        'coverage_basis': 'enclosed_by_successful_quasar_cone'})


def validate_sources(sources):
    allowed = {'gaia', 'tractor', 'ccd', 'desi', *survey_products}
    if not sources or set(sources).difference(allowed):
        raise ValueError('Configure WSDB sources using gaia, tractor, ccd, desi or survey names')
    required = {
        'gaia': {'source_id', 'ra', 'dec', 'phot_g_mean_mag', 'pmra', 'pmdec', 'pmra_error',
                 'pmdec_error', 'pmra_pmdec_corr', 'parallax_over_error', 'ruwe'},
        'tractor': set(TRACTOR_FIELDS),
        'ccd': {'ra', 'dec', 'filter', 'ccd_cuts', 'camera', 'width', 'height', 'crpix1', 'crpix2',
                'crval1', 'crval2', 'cd1_1', 'cd1_2', 'cd2_1', 'cd2_2', 'mjd_obs', 'exptime', 'sig1'},
        'desi': {'targetid', 'release', *[f'{quantity}_{band}' for quantity in ('flux','flux_ivar') for band in 'grz']},
    }
    for name, specification in survey_products.items():
        required[name] = {'source_id', 'ra', 'dec', *[field for fields in specification['bands'].values() for field in fields]}
    required['sdss'].update(['sdss_clean', 'sdss_mode'])
    for name, specification in sources.items():
        for table in source_tables(specification):
            missing = required[name].difference(table['columns'])
            if missing:
                raise ValueError(f'WSDB {name} mapping lacks canonical fields: '+', '.join(sorted(missing)))


def prepare_inputs(config, output, resume=False, client=None, stage='all'):
    allowed = {'_directory', 'quasar_catalogue_config', 'measurement_config', 'imaging_search_config',
               'quasar_catalogue', 'morphology_images', 'imaging_parameters', 'sources', 'public_sources'}
    if set(config).difference(allowed):
        raise ValueError('WSDB config accepts only input config paths and sources; keep connection settings outside JSON')
    sources = config.get('sources', {})
    public = config.get('public_sources', {})
    if (set(public).difference({'sdss','ccd'}) or set(public).intersection(sources)
            or public.get('sdss', 'sdss_dr17') != 'sdss_dr17'
            or public.get('ccd', 'legacy_dr9_ccds') != 'legacy_dr9_ccds'):
        raise ValueError('public_sources supports sdss=sdss_dr17 and ccd=legacy_dr9_ccds, without duplicating WSDB sources')
    if stage == 'one-arcsecond-search':
        sources = {name: specification for name, specification in sources.items() if name == 'gaia'}
    validate_sources(sources)
    inline = 'quasar_catalogue' in config
    if inline:
        if any(name in config for name in ('quasar_catalogue_config','measurement_config','imaging_search_config')):
            raise ValueError('Use the inline WSDB configuration or separate local configs, not both')
        quasar_config = dict(config['quasar_catalogue'], _directory=config['_directory'])
        required = {'gaia'}
        if stage != 'one-arcsecond-search':
            required.add('tractor')
        if stage == 'all':
            required.update({'desi',*survey_products})
        missing = required.difference(sources).difference(public)
        if missing:
            raise ValueError('Configure raw inputs for every required source: '+', '.join(sorted(missing)))
        if stage != 'one-arcsecond-search' and 'ccd' not in sources and 'ccd' not in public:
            raise ValueError('Configure ccd in sources or public_sources for DCR metadata')
    else:
        quasar_config = read_config(local_path(config, config['quasar_catalogue_config']))
    quasars = read_quasars(quasar_config)
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    raw = output / 'raw_catalogues'
    raw.mkdir(exist_ok=True)
    adapted = output / 'quasar_catalogue'
    database = client if client is not None else WsdbCatalogues(resume=resume)
    with database:
        database.check_sources(sources)
        if 'gaia' in sources:
            path, coverage = database.export_cones(sources['gaia'],
                quasars[['object_id', 'ra', 'dec']], raw/'gaia_discovery.csv', 2.)
            quasar_config.update(gaia_catalogue=str(path), gaia_coverage=str(coverage))
        adapt_catalogue(quasar_config, adapted, quasars=quasars)
        # Reuse the ordinary discovery implementation to define photometry requests.
        search_one(one_config(adapted/'one_arcsecond_search.json'), output/'request_positions/one_arcsecond')
        parent = pd.read_csv(output/'request_positions/one_arcsecond/discovered_parent.csv',
            dtype={'matched_gaia_source_id': 'string', 'primary_qso_objid': 'string'})
        write_config(output/'quasar_catalogue.json', absolute_paths(quasar_config))
        write_config(output/'one_arcsecond_search.json', json.loads((adapted/'one_arcsecond_search.json').read_text()))
        if stage == 'one-arcsecond-search':
            summary = {'stage': stage, 'input_quasars': len(quasars),
                'one_arcsecond_parent_entries': len(parent), 'queried_sources': sorted(sources)}
            write_config(output/'input_summary.json', summary)
            print(json.dumps(summary, indent=2), flush=True)
            return summary
        if inline:
            images = {'path_template': str(output/'rband_images/qso_{qso_parent_index:06d}.fits'),
                      **config.get('morphology_images', {})}
            measurements = absolute_paths(dict(morphology_images=images, surveys={}, _directory=config['_directory']))
            imaging = dict(config.get('imaging_parameters', {}))
            if set(imaging).difference({'aperture_arcsec','tractor_radius_arcsec','gaia_radius_arcsec','colour_redshift_half_width'}):
                raise ValueError('imaging_parameters accepts the search radii and colour reference width only')
        else:
            measurements = absolute_paths(read_config(local_path(config, config['measurement_config']))) if stage == 'all' else {}
            imaging = absolute_paths(read_config(local_path(config, config['imaging_search_config'])))
        master = read_master(adapted/'quasar_master.csv')
        gaia = read_frame(local_path(quasar_config, quasar_config['gaia_catalogue'])).drop_duplicates()
        if not gaia.source_id.is_unique:
            raise ValueError('Conflicting Gaia records in raw discovery exports')
        anchor_ids = pd.concat([master.source_id, parent.matched_gaia_source_id]).drop_duplicates()
        anchors = gaia.set_index('source_id').reindex(anchor_ids).rename_axis('source_id').reset_index()[['source_id', 'ra', 'dec']]
        if anchors.ra.isna().any():
            raise ValueError('Gaia extract is missing a discovered anchor')
        if 'gaia' in sources:
            coverage_path = raw/'gaia_anchor_coverage.csv'
            anchor_coverage(anchors, quasars).to_csv(coverage_path, index=False)
            measurements['gaia'] = {'path': str(path), 'coverage': str(coverage_path)}
        imaging['master_catalogue'] = str(adapted/'quasar_master.csv')
        if 'tractor' in sources:
            path, coverage = database.export_cones(sources['tractor'], anchors, raw/'tractor.csv', 10.)
            measurements['tractor'] = {'path': str(path), 'coverage': str(coverage)}
            imaging.update(tractor_catalogues=[str(path)], tractor_coverage=str(coverage))
        imaging_path = write_config(output/'imaging_search.json', imaging)
        request_root = output/'request_positions/imaging'
        request_root.mkdir(parents=True, exist_ok=True)
        discover(imaging_config(imaging_path), request_root)
        if 'gaia' in sources:
            components = pd.read_csv(request_root/'gaia_query_positions.csv')
            targets = components[['component_id', 'component_ra', 'component_dec']].rename(
                columns={'component_ra': 'ra', 'component_dec': 'dec'})
            path, coverage = database.export_cones(sources['gaia'], targets, raw/'gaia_imaging.csv', imaging.get('gaia_radius_arcsec', .7))
            imaging.update(gaia_catalogue=str(path), gaia_coverage=str(coverage))
        ccd_inputs = None
        if 'ccd' in sources or 'ccd' in public:
            targets = pd.read_csv(request_root/'ccd_query_positions.csv').rename(columns={'ra_deg': 'ra', 'dec_deg': 'dec'})
            if 'ccd' in sources:
                ccd_inputs = database.export_cones(sources['ccd'], targets[['qso_parent_index', 'ra', 'dec']], raw/'g_ccds.csv', .65*3600.)
            else:
                ccd_inputs = export_public_sources(public['ccd'], targets[['qso_parent_index', 'ra', 'dec']], raw/'g_ccds.csv', .65*3600., resume)
        if ccd_inputs is not None:
            path, coverage = ccd_inputs
            ccd_coverage = pd.read_csv(coverage)
            ccd_coverage['query_radius_deg'] = ccd_coverage.pop('query_radius_arcsec') / 3600.
            ccd_coverage.to_csv(coverage.with_name('g_ccds_degree_coverage.csv'), index=False)
            imaging.update(ccd_catalogue=str(path), ccd_coverage=str(coverage.with_name('g_ccds_degree_coverage.csv')))
        # Radius 2 around an anchor encloses the 1-arcsec photometry cone around
        # its quasar, including when several catalogue positions share an anchor.
        photometry_targets = anchors.loc[anchors.source_id.isin(parent.matched_gaia_source_id)]
        for name in survey_products:
            if stage == 'all' and name in sources:
                path, coverage = database.export_cones(sources[name], photometry_targets,
                    raw/f'{name}.csv', 2.)
                measurements.setdefault('surveys', {})[name] = {'path': str(path), 'coverage': str(coverage)}
            elif stage == 'all' and name in public:
                path, coverage = export_public_sources(public[name], photometry_targets, raw/f'{name}.csv', 2., resume)
                measurements.setdefault('surveys', {})[name] = {'path': str(path), 'coverage': str(coverage)}
        if stage == 'all' and 'desi' in sources:
            identifiers = [value for values in desi_target_ids(parent).values() for value in values]
            path, _ = database.export_ids(sources['desi'], identifiers, raw/'desi.csv')
            measurements['desi'] = {'path': str(path)}
    write_config(imaging_path, imaging)
    write_config(output/'measurement_inputs.json', measurements)
    write_config(output/'quasar_catalogue.json', absolute_paths(quasar_config))
    # The neighbour-count file here was calculated by adapt_catalogue in this run.
    one = json.loads((adapted/'one_arcsecond_search.json').read_text())
    write_config(output/'one_arcsecond_search.json', one)
    summary = {'input_quasars': len(quasars), 'imaging_anchors': len(master),
        'one_arcsecond_parent_entries': len(parent), 'queried_sources': sorted(sources),
        'public_sources': public,
        'local_replay_configs': ['quasar_catalogue.json', 'one_arcsecond_search.json', 'imaging_search.json', 'measurement_inputs.json'],
        'selection_implementation': 'shared local-file discovery and vetting', 'visual_inspection_applied': False}
    write_config(output/'input_summary.json', summary)
    print(json.dumps(summary, indent=2), flush=True)
    print('Local inputs ready:', output, flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT/'wsdb_inputs.local.json')
    parser.add_argument('--output-root', type=Path, default=WORK_ROOT/'wsdb_inputs')
    parser.add_argument('--resume', action='store_true', help='Reuse hash-verified raw query exports')
    parser.add_argument('--stage', choices=['all', 'one-arcsecond-search', 'imaging-search'], default='all', help='Export only inputs needed by the requested search or full analysis')
    parser.add_argument('--check-schema', action='store_true', help='Check access and native column names without downloading source rows')
    args = parser.parse_args()
    try:
        config = read_config(args.config)
        if args.check_schema:
            validate_sources(config['sources'])
            with WsdbCatalogues() as database:
                database.check_sources(config['sources'])
            print('Configured WSDB tables and columns are accessible; no catalogue rows exported.')
        else:
            prepare_inputs(config, args.output_root, args.resume, stage=args.stage)
    except WsdbInputError as error:
        parser.exit(2, str(error) + '\n')


if __name__ == '__main__':
    main()
