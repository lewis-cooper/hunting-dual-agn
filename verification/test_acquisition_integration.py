"""Exercise raw acquisition, native builders and image downloads without a network."""

from io import BytesIO
import json
import importlib.util
from pathlib import Path
import socket
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd
from astropy.io import fits

from catalogue_downloads import export_public_sources, download_legacy_r_cutout, public_query, CCD_FIELDS, SDSS_FIELDS
from imaging_search.catalogues import resolve_config as imaging_config
from imaging_search.pipeline import run_search as search_imaging
from measurement_inputs import read_config
from one_arcsecond_measurements import build_desi_fallback
from scripts.build_one_arcsecond_measurements import build_measurements
from scripts.prepare_quasar_catalogue import read_quasars
from scripts.query_wsdb_inputs import prepare_inputs
from scripts.measure_imaging_morphology import measure_candidates
from settings import survey_products
from verification import test_native_measurements as native_tests
from verification.test_wsdb_inputs import FixtureCatalogues, FakeConnection, source_spec
from wsdb_inputs import WsdbCatalogues


class FixtureWithIds(FixtureCatalogues):
    def export_ids(self, specification, identifiers, path, identity='targetid'):
        frame = self.frames[specification['table'].split('.')[-1]]
        frame.loc[frame[identity].isin(identifiers)].to_csv(path,index=False)
        return path, None


class AcquisitionIntegrationTests(unittest.TestCase):
    def test_one_inline_config_builds_both_searches_and_native_measurements(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            _, gaia, tractor = native_tests.native_fixture()
            gaia['parallax_over_error'], gaia['ruwe'] = 0., 1.
            pd.DataFrame(dict(object_id=['master:one'], ra=[30.], dec=[0.], z=[1.],
                native_id=['["400","401"]'])).to_csv(root/'quasars.csv',index=False)
            frames = dict(gaia=gaia, tractor=tractor)
            for name, specification in survey_products.items():
                frame = gaia[['source_id','ra','dec']].copy()
                for magnitude, error in specification['bands'].values():
                    frame[magnitude], frame[error] = 19., .1
                if name == 'sdss':
                    frame['sdss_clean'], frame['sdss_mode'] = True, 1
                frames[name] = frame
            frames['desi'] = pd.DataFrame({'targetid':['400','401'], 'release':[9010,9010],
                **{f'{quantity}_{band}':[100.,100.] if quantity=='flux' else [1.,1.]
                   for quantity in ('flux','flux_ivar') for band in 'grz'}})
            frames['ccd'] = pd.DataFrame([dict(ra=30.,dec=0.,filter='g',ccd_cuts=0,camera='decam',
                width=2048,height=2048,crpix1=1024.,crpix2=1024.,crval1=30.,crval2=0.,
                cd1_1=-.262/3600,cd1_2=0.,cd2_1=0.,cd2_2=.262/3600,mjd_obs=60000.,exptime=100.,sig1=1.)])
            database = FixtureWithIds(frames)
            config = dict(_directory=str(root),quasar_catalogue=dict(catalogue='quasars.csv',
                catalogue_label='DESI_DR1_QSO',redshift_kind='spectroscopic',
                columns=dict(object_id='object_id',ra='ra',dec='dec',redshift='z',native_id='native_id'),
                catalogue_object_id_columns={'DESI_DR1_QSO':'native_id'},catalogue_object_id_format='json_array'),
                sources={name:source_spec(name,frame) for name,frame in frames.items() if name not in ('ccd','sdss')},
                public_sources={'ccd':'legacy_dr9_ccds','sdss':'sdss_dr17'},morphology_images={'download_missing':True})

            def public_export(service, targets, path, radius, resume):
                name = 'ccd' if service == 'legacy_dr9_ccds' else 'sdss'
                return database.export_cones(source_spec(name,frames[name]),targets,path,radius)

            with (patch.object(socket.socket,'connect',side_effect=AssertionError('Network forbidden')),
                  patch('scripts.query_wsdb_inputs.export_public_sources',side_effect=public_export)):
                summary = prepare_inputs(config,root/'generated',client=database)
                parent = pd.read_csv(root/'generated/request_positions/one_arcsecond/discovered_parent.csv',
                                     dtype={'matched_gaia_source_id':'string','primary_qso_objid':'string'})
                self.assertEqual(summary['one_arcsecond_parent_entries'],1)
                self.assertEqual(json.loads(parent.qso_catalogue_object_ids_json.item()),['400','401'])
                build_measurements(read_config(root/'generated/measurement_inputs.json'),parent,root/'measurements')
                fallback = pd.read_csv(root/'measurements/desi_target_photometry.csv')
                self.assertEqual(fallback.desi_target_photometry_status.item(),'matched_exact_target')
                self.assertEqual(fallback.desi_flux_r.item(),100.)
                result = search_imaging(imaging_config(root/'generated/imaging_search.json'),root/'imaging')
                self.assertEqual(result['counts']['two_gaia'],1)
                self.assertEqual(result['counts']['systems_with_dcr'],1)

    def test_conflicting_desi_targets_are_flagged_instead_of_selected_arbitrarily(self):
        parent = pd.DataFrame(dict(qso_group_id=['pair'],primary_qso_catalog=['DESI_DR1_QSO'],
            primary_qso_objid=['master:one'],qso_catalogue_object_ids_json=['["400","401"]']))
        photometry = pd.DataFrame({'targetid':['400','401'],'release':[9010,9010],
            **{f'{quantity}_{band}':[100.,200.] if quantity=='flux' else [1.,1.]
               for quantity in ('flux','flux_ivar') for band in 'grz'}})
        result = build_desi_fallback(parent,photometry)
        self.assertEqual(result.desi_target_photometry_status.item(),'conflicting_native_targets')
        self.assertTrue(np.isnan(result.desi_flux_r.item()))

    def test_mixed_master_redshift_types_and_errors_are_preserved(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            pd.DataFrame(dict(id=['spec','photo'],ra=[10.,20.],dec=[0.,0.],z=[1.,2.],
                catalogue=['sdss','quaia'],error=['[]','[0.3]'])).to_csv(root/'master.csv',index=False)
            config=dict(_directory=str(root),catalogue='master.csv',
                columns=dict(object_id='id',ra='ra',dec='dec',redshift='z',catalogue_label='catalogue',redshift_error='error'),
                redshift_kind_by_catalogue={'sdss':'spectroscopic','quaia':'photometric'},
                redshift_error_format='json_singleton',catalogue_labels={'sdss':'SDSS_DR16Q','quaia':'QUAIA_G20_V1'})
            quasars=read_quasars(config)
            self.assertEqual(quasars.redshift_error.tolist(),[0.,.3])
            self.assertEqual(quasars.redshift_kind.tolist(),['spectroscopic','photometric'])

    def test_missing_cutout_download_is_measured_and_then_reused(self):
        image,header,ra,dec = native_tests.MorphologyTests.image_fixture()
        payload = BytesIO()
        fits.PrimaryHDU(image,header).writeto(payload)
        response = Mock(content=payload.getvalue())
        candidates = pd.DataFrame(dict(qso_parent_index=[0],primary_ra=[ra[0]],primary_dec=[dec[0]],
                                      secondary_ra=[ra[1]],secondary_dec=[dec[1]],gaia_ra=[ra[0]],gaia_dec=[dec[0]]))
        with TemporaryDirectory() as directory:
            config=dict(_directory=directory,morphology_images=dict(path_template='{qso_parent_index}.fits',download_missing=True))
            with patch('catalogue_downloads.requests.get',return_value=response) as request:
                measured=measure_candidates(candidates,config)
                reused=measure_candidates(candidates,config)
            self.assertEqual(request.call_count,1)
            self.assertGreater(measured.r_perp_sigma_ratio.item(),1.05)
            self.assertEqual(measured.r_perp_sigma_ratio.item(),reused.r_perp_sigma_ratio.item())
            self.assertEqual(request.call_args.kwargs['params']['bands'],'r')
            self.assertEqual(float(request.call_args.kwargs['params']['pixscale']),.1)
            self.assertEqual(int(request.call_args.kwargs['params']['size']),100)
            self.assertEqual(request.call_args.kwargs['params']['ra'],f'{ra[0]:.7f}')

    def test_download_rejects_wrong_field_before_writing_image(self):
        image,header,ra,dec=native_tests.MorphologyTests.image_fixture()
        payload=BytesIO();fits.PrimaryHDU(image,header).writeto(payload)
        row=pd.Series(dict(primary_ra=100.,primary_dec=0.,secondary_ra=100.0001,secondary_dec=0.))
        with TemporaryDirectory() as directory:
            path=Path(directory)/'image.fits'
            session=Mock();session.get.return_value=Mock(content=payload.getvalue())
            with self.assertRaisesRegex(ValueError,'contain both'):
                download_legacy_r_cutout(row,path,session=session)
            self.assertFalse(path.exists())

    def test_public_empty_results_have_coverage_but_errors_do_not(self):
        targets=pd.DataFrame(dict(source_id=['9007199254740993'],ra=[10.],dec=[0.]))
        _,parameters,_=public_query('sdss_dr17',targets,2.)
        self.assertTrue(parameters['cmd'].lstrip().startswith('SELECT'))
        with TemporaryDirectory() as directory:
            root=Path(directory)
            session=Mock();session.get.return_value=Mock(text=','.join(SDSS_FIELDS)+'\n')
            path,coverage=export_public_sources('sdss_dr17',targets,root/'sdss.csv',2.,session=session)
            self.assertTrue(pd.read_csv(path).empty)
            self.assertEqual(pd.read_csv(coverage,dtype={'source_id':'string'}).source_id.item(),'9007199254740993')
            session.get.return_value=Mock(text='service error')
            with self.assertRaises(ValueError):
                export_public_sources('sdss_dr17',targets,root/'failed.csv',2.,session=session)
            self.assertFalse((root/'failed_coverage.csv').exists())
            session.get.return_value=Mock(text=','.join(['qso_parent_index','target_ra','target_dec',*CCD_FIELDS])+'\n')
            export_public_sources('legacy_dr9_ccds',targets.rename(columns={'source_id':'qso_parent_index'}).assign(qso_parent_index=0),root/'ccds.csv',2340.,session=session)
            self.assertEqual(session.get.call_args.kwargs['headers']['X-DL-AuthToken'],'anonymous.0.0.anon_access')

    @unittest.skipUnless(importlib.util.find_spec('psycopg'), 'Optional WSDB driver is not installed')
    def test_wsdb_batches_have_one_header_and_resume_verified_raw_chunks(self):
        with TemporaryDirectory() as directory:
            path=Path(directory)/'gaia.csv'
            connection=FakeConnection(b'source_id,ra,dec\n9007199254740993,10,0\n')
            database=WsdbCatalogues(resume=True,batch_size=1,connect=lambda **options:connection)
            specification=source_spec('gaia',pd.DataFrame(columns=['source_id','ra','dec']))
            targets=pd.DataFrame(dict(object_id=['one','two'],ra=[10.,20.],dec=[0.,0.]))
            database.export_cones(specification,targets,path,2.)
            self.assertEqual(path.read_text().count('source_id,ra,dec'),1)
            before=len(connection.queries)
            path.unlink()
            database.export_cones(specification,targets,path,2.)
            self.assertEqual(before,len(connection.queries))


if __name__=='__main__':
    unittest.main()
