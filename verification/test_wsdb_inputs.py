"""WSDB/local equivalence and failure handling, without live database queries."""

from contextlib import nullcontext
import importlib.util
import json
from pathlib import Path
import socket
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import search_around_sky

from dual_agn.imaging_search.catalogues import read_frame, resolve_config as imaging_config
from dual_agn.imaging_search.pipeline import run_search as search_imaging
from dual_agn.measurement_inputs import read_config, sky
from dual_agn.one_arcsecond_search import resolve_config as one_config, run_search as search_one
from scripts.gaia_neighbour.build_measurements import build_measurements
from scripts.shared.acquire_inputs import absolute_paths, prepare_inputs, validate_sources
from scripts.shared.prepare_quasar_catalogue import adapt_catalogue
from run_analysis import commands
from dual_agn.wsdb_inputs import WsdbCatalogues, WsdbInputError, source_tables
from dual_agn.settings import survey_products


def source_spec(name, frame):
    return dict(table=f'example.{name}', release='test release',
                columns={name: name for name in frame})


class FixtureCatalogues:
    """Stand in for a database using complete native catalogues and real cones."""
    def __init__(self, frames):
        self.frames = frames
        self.requests = []

    def __enter__(self):
        return self

    def __exit__(self, *exception):
        pass

    def check_sources(self, sources):
        for specification in sources.values():
            source_tables(specification)

    def export_cones(self, specification, targets, path, radius_arcsec):
        frame = self.frames[specification['table'].split('.')[-1]]
        _, source_index, _, _ = search_around_sky(sky(targets), sky(frame), radius_arcsec*u.arcsec)
        selected = frame.iloc[np.unique(source_index)].copy()
        path.parent.mkdir(parents=True, exist_ok=True)
        selected.to_csv(path, index=False)
        coverage = targets.drop(columns=['ra', 'dec']).copy()
        coverage['ra_deg'] = targets.ra
        coverage['dec_deg'] = targets.dec
        coverage['query_radius_arcsec'] = radius_arcsec
        coverage_path = path.with_name(path.stem+'_coverage.csv')
        coverage.to_csv(coverage_path, index=False)
        self.requests.append((path.name, len(targets), radius_arcsec))
        return path, coverage_path

    def export_ids(self, specification, identifiers, path, identity='targetid'):
        frame = self.frames[specification['table'].split('.')[-1]]
        frame.loc[frame[identity].isin(identifiers)].to_csv(path, index=False)
        return path, None


class WsdbWorkflowTests(unittest.TestCase):
    def test_wsdb_exports_feed_the_same_offline_discovery_and_measurements(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            quasars = pd.DataFrame(dict(id=['pair', 'isolated', 'empty'], ra=[10., 20., 30.],
                dec=[0., 0., 0.], z=[1., 1.2, 2.]))
            quasars.to_csv(root/'quasars.csv', index=False)
            gaia = pd.DataFrame(dict(source_id=['9007199254740993', '9007199254740994', '9007199254740995'],
                ra=[10., 10.+.72/3600., 20.], dec=[0., 0., 0.], phot_g_mean_mag=[18., 19., 18.],
                pmra=[0., np.nan, 0.], pmdec=[0., 0., 0.], pmra_error=[1., 1., 1.],
                pmdec_error=[1., 1., 1.], pmra_pmdec_corr=[0., 0., 0.],
                parallax_over_error=[0., 0., 0.], ruwe=[1., 1., 1.]))
            tractor = pd.DataFrame(dict(release=[9010]*3, brickid=[1, 1, 2], objid=[1, 2, 1],
                brick_primary=True, type='PSF', ra=gaia.ra, dec=gaia.dec, maskbits=0, fracflux_r=.1))
            for band in 'grz':
                tractor[f'flux_{band}'] = [100., 50., 100.]
                tractor[f'flux_ivar_{band}'] = 100.
                tractor[f'nobs_{band}'] = 3
            ccd = pd.DataFrame([dict(ra=10., dec=0., filter='g', ccd_cuts=0, camera='decam',
                width=2048, height=2048, crpix1=1024., crpix2=1024., crval1=10., crval2=0.,
                cd1_1=-.262/3600, cd1_2=0., cd2_1=0., cd2_2=.262/3600,
                mjd_obs=60000., exptime=100., sig1=1.)])
            ps1 = gaia[['source_id', 'ra', 'dec']].copy()
            for band in 'grizy':
                ps1[f'ps1_{band}kronmag'] = 19.
                ps1[f'ps1_{band}kronmagerr'] = .1
            frames = dict(gaia=gaia, tractor=tractor, ccd=ccd, ps1=ps1)
            for name, specification in survey_products.items():
                frame = gaia[['source_id', 'ra', 'dec']].copy()
                for magnitude, error in specification['bands'].values():
                    frame[magnitude], frame[error] = 19., .1
                frames[name] = frame
            frames['desi'] = pd.DataFrame(columns=['targetid', 'release',
                *[f'{quantity}_{band}' for quantity in ('flux', 'flux_ivar') for band in 'grz']])
            config = dict(_directory=str(root), quasar_catalogue=dict(
                catalogue='quasars.csv', catalogue_label='test', redshift_kind='spectroscopic',
                columns=dict(object_id='id', ra='ra', dec='dec', redshift='z')),
                sources={name: source_spec(name, frame) for name, frame in frames.items()})
            database = FixtureCatalogues(frames)
            with patch.object(socket.socket, 'connect', side_effect=AssertionError('Live query forbidden')):
                result = prepare_inputs(config, root/'generated', client=database)
                self.assertEqual(result['input_quasars'], 3)
                self.assertEqual(result['one_arcsecond_parent_entries'], 1)
                self.assertEqual(result['imaging_anchors'], 2)
                generated = root/'generated'
                # This is the same entry point used by a user with local data.
                offline = search_one(one_config(generated/'one_arcsecond_search.json'), root/'offline_one')
                self.assertEqual(offline['counts']['parent_entries'], 1)
                parent = pd.read_csv(root/'offline_one/discovered_parent.csv',
                    dtype={'matched_gaia_source_id':'string', 'primary_qso_objid':'string'})
                build_measurements(read_config(generated/'measurement_inputs.json'), parent, root/'measurements')
                systems = pd.read_csv(root/'measurements/component_quality_audit.csv')
                self.assertEqual(systems.decals_resolution_class.tolist(), ['decals_resolved'])
                self.assertTrue(systems.secondary_pm_significance.isna().all())
                surveys = pd.read_csv(root/'measurements/ps1_measurements.csv')
                self.assertTrue(surveys.ps1_source_matched.all())
                imaging_result = search_imaging(imaging_config(generated/'imaging_search.json'), root/'offline_imaging')
                self.assertEqual(imaging_result['counts']['two_gaia'], 1)
                self.assertEqual(imaging_result['counts']['systems_with_dcr'], 1)
                raw_gaia = read_frame(generated/'raw_catalogues/gaia_discovery.csv')
                pd.testing.assert_frame_equal(raw_gaia, gaia.astype({'source_id':'string'}), check_exact=False)
                pd.testing.assert_frame_equal(read_frame(generated/'raw_catalogues/tractor.csv'), tractor, check_exact=False)
                # Recompute discovery from the original master plus raw exports,
                # without reading the earlier run's neighbour-count table.
                adapt_catalogue(read_config(generated/'quasar_catalogue.json'), root/'fresh_local_master',
                    generated/'imaging_search.json')
                pd.testing.assert_frame_equal(pd.read_csv(root/'fresh_local_master/gaia_neighbour_counts.csv'),
                    pd.read_csv(generated/'quasar_catalogue/gaia_neighbour_counts.csv'))
                fresh = search_one(one_config(root/'fresh_local_master/one_arcsecond_search.json'), root/'fresh_local_one')
                self.assertEqual(fresh['counts'], offline['counts'])
            # Discovery queried the empty and isolated fields too, rather than
            # starting from the known pair. Imaging Gaia is queried at components.
            self.assertEqual(database.requests[0], ('gaia_discovery.csv', 3, 2.))
            coverage = pd.read_csv(generated/'raw_catalogues/gaia_discovery_coverage.csv')
            self.assertEqual(coverage.object_id.tolist(), ['pair', 'isolated', 'empty'])

    def test_runner_exports_before_local_builders_and_preserves_offline_mode(self):
        sequence = commands('vetting', wsdb_config=Path('wsdb.json'))
        self.assertIn('--models-only', sequence[0])
        self.assertIn('scripts.shared.acquire_inputs', sequence[1])
        modules = [command[command.index('-m')+1] for command in sequence if '-m' in command]
        self.assertLess(modules.index('scripts.gaia_neighbour.search_candidates'), modules.index('scripts.gaia_neighbour.score_resolved'))
        self.assertFalse(any('scripts.shared.prepare_quasar_catalogue' in command for command in sequence))
        offline = commands('vetting', data_source='local')
        self.assertFalse(any('scripts.shared.acquire_inputs' in command for command in offline))
        with self.assertRaisesRegex(ValueError, 'do not combine'):
            commands('vetting', wsdb_config=Path('w.json'), measurement_config=Path('m.json'))
        self.assertIn('--resume', commands('vetting', wsdb_config=Path('w.json'), resume=True)[1])

    def test_default_uses_wsdb_and_local_choice_is_explicit(self):
        sequence = commands('vetting')
        self.assertIn('--models-only', sequence[0])
        self.assertIn('scripts.shared.acquire_inputs', sequence[1])
        self.assertEqual(Path(sequence[1][-1]).name, 'wsdb_inputs.local.json')
        for stage in ('one-arcsecond-search', 'imaging-search'):
            sequence = commands(stage)
            self.assertEqual(sequence[0][-2:], ['--stage', stage])
            self.assertEqual(len(sequence), 2)
        with self.assertRaisesRegex(ValueError, 'cannot be used'):
            commands('vetting', data_source='local', wsdb_config=Path('w.json'))
        controls = commands('chance-alignments', chance_config=Path('controls.json'))
        self.assertEqual(len(controls), 1)
        self.assertIn('run_chance_alignments', controls[0])

    def test_one_arcsecond_only_does_not_query_imaging_or_load_its_inputs(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            pd.DataFrame(dict(id=['pair'], ra=[10.], dec=[0.], z=[1.])).to_csv(root/'quasars.csv', index=False)
            quasar_config = dict(catalogue='quasars.csv', catalogue_label='test', redshift_kind='spectroscopic',
                columns=dict(object_id='id', ra='ra', dec='dec', redshift='z'))
            gaia = pd.DataFrame(dict(source_id=['9007199254740993','9007199254740994'],
                ra=[10.,10.+.7/3600], dec=[0.,0.], phot_g_mean_mag=[18.,19.],
                pmra=[0.,0.], pmdec=[0.,0.], pmra_error=[1.,1.], pmdec_error=[1.,1.],
                pmra_pmdec_corr=[0.,0.], parallax_over_error=[0.,0.], ruwe=[1.,1.]))
            config = dict(_directory=str(root), quasar_catalogue=quasar_config,
                sources={'gaia':source_spec('gaia',gaia)})
            database = FixtureCatalogues({'gaia':gaia})
            result = prepare_inputs(config, root/'generated', client=database, stage='one-arcsecond-search')
            self.assertEqual(result['one_arcsecond_parent_entries'],1)
            self.assertEqual(database.requests,[('gaia_discovery.csv',1,2.)])
            self.assertTrue((root/'generated/one_arcsecond_search.json').is_file())

    def test_relocated_configs_keep_original_relative_native_paths(self):
        config = dict(_directory='/example', surveys={'ps1': {'path':'ps1.csv', 'columns': {'source_id':'objid'}}},
            morphology_images={'path_template':'images/{qso_parent_index}.fits'})
        actual = absolute_paths(config)
        self.assertEqual(actual['surveys']['ps1']['path'], '/example/ps1.csv')
        self.assertEqual(actual['surveys']['ps1']['columns'], {'source_id':'objid'})
        self.assertEqual(actual['morphology_images']['path_template'], '/example/images/{qso_parent_index}.fits')

    def test_mapping_rejects_sql_expressions_and_credentials(self):
        for specification in [dict(table='x.y; DROP TABLE x', release='one', columns={'ra':'ra'}),
                dict(table='x.y', release='one', columns={'ra':'coalesce(ra,0)'}),
                dict(table='x.y', release='one', columns={'ra':'ra'}, password='private')]:
            with self.subTest(specification=specification), self.assertRaises(ValueError):
                source_tables(specification)

    def test_missing_measurement_fields_fail_before_any_connection(self):
        with self.assertRaisesRegex(ValueError, 'pmra'):
            validate_sources({'gaia': source_spec('gaia', pd.DataFrame(columns=['source_id','ra','dec']))})

    def test_connection_errors_do_not_expose_credentials(self):
        def failing_connect(**options):
            raise OSError('private_user super_secret_password')
        with self.assertRaises(WsdbInputError) as caught:
            WsdbCatalogues(connect=failing_connect).connection()
        self.assertNotIn('private_user', str(caught.exception))
        self.assertNotIn('super_secret_password', str(caught.exception))


class FakeConnection:
    def __init__(self, csv_bytes=None, fail=False):
        self.csv_bytes = csv_bytes
        self.fail = fail
        self.queries = []

    def transaction(self):
        return nullcontext()

    def cursor(self):
        return nullcontext(self)

    def execute(self, query):
        self.queries.append(str(query))

    def copy(self, query):
        class Copied:
            def write_row(self, row):
                pass
        if isinstance(query, str):
            return nullcontext(Copied())
        self.queries.append(query.as_string())
        if self.fail:
            raise OSError('private_user private_password')
        return nullcontext(iter([self.csv_bytes]))

    def close(self):
        pass


@unittest.skipUnless(importlib.util.find_spec('psycopg'), 'Optional WSDB driver is not installed')
class WsdbExportTests(unittest.TestCase):
    def test_exact_identity_export_and_empty_identity_request(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'desi.csv'
            specification = source_spec('desi', pd.DataFrame(columns=['targetid','flux_g']))
            connection = FakeConnection(b'targetid,flux_g\n9007199254740993,100\n')
            database = WsdbCatalogues(connect=lambda **options: connection)
            database.export_ids(specification, ['9007199254740993'], path)
            self.assertIn('s."targetid" = t.identifier::bigint', connection.queries[-1])
            self.assertEqual(pd.read_csv(path, dtype='string').targetid.item(), '9007199254740993')
            calls = len(connection.queries)
            database.export_ids(specification, [], path)
            self.assertEqual(calls, len(connection.queries))
            self.assertEqual(pd.read_csv(path).columns.tolist(), ['targetid','flux_g'])

    def test_multi_table_export_quotes_identifiers_and_unions_native_rows(self):
        with TemporaryDirectory() as directory:
            north = source_spec('north', pd.DataFrame(columns=['ra','dec']))
            south = source_spec('south', pd.DataFrame(columns=['ra','dec']))
            connection = FakeConnection(b'ra,dec\n0,0\n')
            database = WsdbCatalogues(connect=lambda **options: connection)
            database.export_cones({'tables':[north,south]}, pd.DataFrame(dict(ra=[0.], dec=[0.])), Path(directory)/'tractor.csv', 10.)
            query = connection.queries[-1]
            self.assertIn('"example"."north"', query)
            self.assertIn('"example"."south"', query)
            self.assertIn(' UNION ', query)

    def test_empty_successful_query_records_coverage_and_ids_remain_exact(self):
        targets = pd.DataFrame(dict(object_id=['empty'], ra=[359.99999], dec=[0.]))
        with TemporaryDirectory() as directory:
            path = Path(directory)/'gaia.csv'
            specification = source_spec('gaia', pd.DataFrame(columns=['source_id', 'ra', 'dec']))
            connection = FakeConnection(b'source_id,ra,dec\n')
            with WsdbCatalogues(connect=lambda **options: connection) as database:
                _, coverage = database.export_cones(specification, targets, path, 2.)
            self.assertTrue(path.is_file())
            self.assertEqual(pd.read_csv(coverage).object_id.tolist(), ['empty'])
            connection.csv_bytes = b'source_id,ra,dec\n9007199254740993,0,0\n'
            with WsdbCatalogues(connect=lambda **options: connection) as database:
                database.export_cones(specification, targets, path, 2.)
            self.assertEqual(read_frame(path).source_id.tolist(), ['9007199254740993'])
            self.assertIn('q3c_join', connection.queries[-1])
            self.assertNotIn('LIMIT', connection.queries[-1])

    def test_failed_query_has_no_coverage_or_partial_output(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            database = WsdbCatalogues(connect=lambda **options: FakeConnection(fail=True))
            with self.assertRaises(WsdbInputError) as caught:
                database.export_cones(source_spec('gaia', pd.DataFrame(columns=['ra','dec'])),
                    pd.DataFrame(dict(ra=[0.], dec=[0.])), root/'gaia.csv', 2.)
            self.assertNotIn('private_password', str(caught.exception))
            self.assertFalse(list(root.iterdir()))

    def test_resume_checks_targets_configuration_and_output_hash(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            specification = source_spec('gaia', pd.DataFrame(columns=['ra','dec']))
            targets = pd.DataFrame(dict(ra=[0.], dec=[0.]))
            connection = FakeConnection(b'ra,dec\n0,0\n')
            database = WsdbCatalogues(resume=True, connect=lambda **options: connection)
            database.export_cones(specification, targets, root/'gaia.csv', 2.)
            previous = len(connection.queries)
            database.export_cones(specification, targets, root/'gaia.csv', 2.)
            self.assertEqual(previous, len(connection.queries))
            (root/'gaia.csv').write_text('changed')
            database.export_cones(specification, targets, root/'gaia.csv', 2.)
            self.assertGreater(len(connection.queries), previous)
            previous = len(connection.queries)
            database.export_cones(specification, targets.assign(ra=1.), root/'gaia.csv', 2.)
            self.assertGreater(len(connection.queries), previous)


if __name__ == '__main__':
    unittest.main()
