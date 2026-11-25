"""Native association and morphology tests, without survey downloads or model fits."""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

from dual_agn.imaging_morphology import measure_image
from dual_agn.measurement_inputs import read_config, read_local_sources, require_coverage
from dual_agn.one_arcsecond_measurements import (build_components, build_system_measurements,
                                       build_survey_photometry)
from scripts.run_analysis import commands
from scripts.gaia_neighbour.build_measurements import build_measurements
from scripts.imaging.measure_morphology import measure_candidates
from dual_agn.settings import survey_products
from dual_agn.imaging_search.catalogues import write_fits
from scripts.shared.prepare_quasar_catalogue import adapt_catalogue
from dual_agn.one_arcsecond_search import resolve_config as discovery_config, run_search
from dual_agn.redshift_inputs import redshift_errors


def native_fixture(shared=False):
    parent = pd.DataFrame(dict(qso_group_id=['test'], matched_gaia_source_id=['100000000000000001'],
        gaia_edr3_cnt1=[2], primary_qso_ra=[30.], primary_qso_dec=[0.]))
    separation = .15 if shared else .7
    gaia = pd.DataFrame(dict(source_id=['100000000000000001', '100000000000000002'],
        ra=[30., 30. + separation / 3600], dec=[0., 0.],
        pmra=[0., 2.], pmdec=[0., 0.], pmra_error=[1., 1.], pmdec_error=[1., 1.],
        pmra_pmdec_corr=[0., .5], phot_g_mean_mag=[18., 19.]))
    tractor = pd.DataFrame(dict(release=[9010, 9010], brickid=[100, 100], objid=[1, 2],
        type=['PSF', 'PSF'], ra=gaia.ra, dec=gaia.dec, maskbits=[0, 0], fracflux_r=[.2, .2]))
    if shared:
        tractor = tractor.iloc[:1].copy()
        tractor['ra'] = 30. + .075 / 3600
    for band in 'grz':
        tractor[f'flux_{band}'] = 100.
        tractor[f'flux_ivar_{band}'] = 1.
        tractor[f'nobs_{band}'] = 2
    return parent, gaia, tractor


class NativeMeasurementTests(unittest.TestCase):
    def test_calibration_only_build_needs_no_supplementary_catalogues(self):
        parent,gaia,tractor=native_fixture()
        parent=parent.assign(primary_qso_catalog='test', primary_qso_objid='test',
            primary_qso_z=1., primary_qso_redshift_kind='spectroscopic')
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            gaia.to_csv(root/'gaia.csv',index=False)
            tractor.to_csv(root/'tractor.csv',index=False)
            pd.DataFrame(dict(source_id=parent.matched_gaia_source_id,ra_deg=[30.],
                dec_deg=[0.],query_radius_arcsec=[2.])).to_csv(root/'coverage.csv',index=False)
            config=dict(_directory=directory,gaia=dict(path='gaia.csv',coverage='coverage.csv'),
                tractor=dict(path='tractor.csv',coverage='coverage.csv'))
            build_measurements(config,parent,root/'measurements',calibration_only=True)
            self.assertTrue(pd.read_csv(root/'measurements/component_quality_audit.csv').eligible_both_components.all())
            for survey in survey_products:
                photometry=pd.read_csv(root/'measurements'/f'{survey}_measurements.csv')
                self.assertFalse(photometry[f'{survey}_source_matched'].any())

    def test_distinct_assignments_and_covariant_pm(self):
        parent, gaia, tractor = native_fixture()
        components = build_components(parent, gaia, tractor)
        systems = build_system_measurements(parent, components)
        self.assertEqual(systems.decals_resolution_class.iloc[0], 'decals_resolved')
        self.assertTrue(systems.eligible_both_components.iloc[0])
        self.assertAlmostEqual(systems.secondary_pm_significance.iloc[0], 2 / np.sqrt(.75))
        self.assertEqual(systems.secondary_component_gaia_source_id.iloc[0], '100000000000000002')

    def test_shared_source_is_withheld_from_both_components(self):
        parent, gaia, tractor = native_fixture(shared=True)
        components = build_components(parent, gaia, tractor)
        systems = build_system_measurements(parent, components)
        self.assertTrue(components.decals_assignment_status.eq('shared_unresolved_nearest_decals_object').all())
        self.assertTrue(components.decals_flux_r.isna().all())
        self.assertFalse(systems.eligible_both_components.iloc[0])
        self.assertNotEqual(systems.decals_resolution_class.iloc[0], 'decals_resolved')

    def test_two_detections_are_structurally_resolved_even_when_flux_fails(self):
        parent, gaia, tractor = native_fixture()
        tractor.loc[1, 'flux_g'] = -1.
        systems = build_system_measurements(parent, build_components(parent, gaia, tractor))
        self.assertEqual(systems.decals_resolution_class.iloc[0], 'decals_resolved')
        self.assertFalse(systems.eligible_both_components.iloc[0])

    def test_missing_pm_is_preserved_for_downstream_flagging(self):
        parent, gaia, tractor = native_fixture()
        gaia.loc[1, 'pmra_error'] = np.nan
        systems = build_system_measurements(parent, build_components(parent, gaia, tractor))
        self.assertTrue(np.isnan(systems.secondary_pm_significance.iloc[0]))

    def test_missing_coverage_is_an_error(self):
        targets = pd.DataFrame(dict(source_id=['100'], ra=[30.], dec=[0.]))
        with tempfile.TemporaryDirectory() as directory:
            coverage = Path(directory) / 'coverage.csv'
            pd.DataFrame(dict(source_id=['101'], ra_deg=[30.], dec_deg=[0.], query_radius_arcsec=[2.])).to_csv(coverage, index=False)
            with self.assertRaisesRegex(ValueError, 'lack complete'):
                require_coverage(targets, {'_directory':directory}, {'coverage':'coverage.csv'}, 1.)

    def test_native_survey_match_and_band_flags_are_derived(self):
        parent, _, _ = native_fixture()
        with tempfile.TemporaryDirectory() as directory:
            pd.DataFrame(dict(native_id=['a', 'b'], ra=[30. + .2/3600, 30. + .8/3600], dec=[0., 0.],
                              rmag=[19., 18.], rerr=[.5, .1])).to_csv(Path(directory)/'survey.csv', index=False)
            specification = {'path':'survey.csv', 'full_catalogue':True,
                             'columns':{'source_id':'native_id', 'ps1_rkronmag':'rmag', 'ps1_rkronmagerr':'rerr'}}
            result = build_survey_photometry(parent, {'_directory':directory}, 'ps1', specification,
                                            {'r':('ps1_rkronmag', 'ps1_rkronmagerr')})
            self.assertEqual(result.ps1_source_id.iloc[0], 'a')
            self.assertTrue(result.ps1_source_matched.iloc[0])
            self.assertFalse(result.ps1_r_band_quality.iloc[0])

    def test_runner_generates_measurements_before_validation_and_scoring(self):
        sequence = commands('vetting', measurement_config=Path('native.json'), data_source='local')
        modules = [command[command.index('-m')+1] for command in sequence if '-m' in command]
        self.assertLess(modules.index('scripts.gaia_neighbour.build_measurements'), modules.index('scripts.gaia_neighbour.prepare_inputs'))
        self.assertLess(modules.index('scripts.imaging.measure_morphology'), modules.index('scripts.imaging.select_candidates'))
        for command in sequence:
            if any(name in command for name in ('scripts.gaia_neighbour.build_measurements', 'scripts.imaging.measure_morphology')):
                self.assertEqual(Path(command[command.index('--config')+1]).name, 'native.json')

    def test_native_files_build_all_required_intermediates(self):
        parent, gaia, tractor = native_fixture()
        parent['primary_qso_catalog'] = 'DESI_DR1_QSO'
        parent['primary_qso_objid'] = '12345'
        parent['primary_qso_z'] = 1.
        parent['primary_qso_zerr'] = 0.
        parent['primary_qso_name'] = 'test'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gaia.to_csv(root/'gaia.csv', index=False)
            tractor.to_csv(root/'tractor.csv', index=False)
            pd.DataFrame({'targetid':['12345'], 'release':[9010],
                **{f'{quantity}_{band}':[100. if quantity == 'flux' else 1.] for quantity in ('flux','flux_ivar') for band in 'grz'}}).to_csv(root/'desi.csv', index=False)
            write_fits(pd.DataFrame({'source_id':['100000000000000001'], 'redshift_quaia':[1.], 'redshift_quaia_err':[.1]}), root/'quaia.fits')
            pd.DataFrame({'source_id':['100000000000000001'], 'ra_deg':[30.], 'dec_deg':[0.], 'query_radius_arcsec':[2.]}).to_csv(root/'coverage.csv',index=False)
            config = {'_directory':directory, 'gaia':{'path':'gaia.csv','coverage':'coverage.csv'},
                      'tractor':{'path':'tractor.csv','coverage':'coverage.csv'},
                      'desi':{'path':'desi.csv'}, 'quaia_catalogue':'quaia.fits', 'surveys':{}}
            for survey, specification in survey_products.items():
                source = {'source_id':['native'], 'ra':[30.], 'dec':[0.]}
                for magnitude, error in specification['bands'].values():
                    source[magnitude], source[error] = [19.], [.1]
                pd.DataFrame(source).to_csv(root/f'{survey}.csv', index=False)
                config['surveys'][survey] = {'path':f'{survey}.csv','full_catalogue':True}
            generated = root/'work/one_arcsecond_inputs'
            build_measurements(config, parent, generated)
            self.assertEqual(pd.read_csv(generated/'component_quality_audit.csv').eligible_both_components.tolist(), [True])
            self.assertEqual(pd.read_csv(generated/'one_arcsecond_resolution.csv').decals_resolution_class.tolist(), ['decals_resolved'])
            self.assertEqual(pd.read_csv(generated/'desi_target_photometry.csv').desi_flux_r.tolist(), [100.])
            for survey in survey_products:
                self.assertTrue((generated/f'{survey}_measurements.csv').is_file())
            parent_path = root/'work/one_arcsecond_search/discovered_parent.csv'
            parent_path.parent.mkdir(parents=True)
            parent.to_csv(parent_path,index=False)
            environment = dict(os.environ, PAPER1_WORK_ROOT=str(root/'work'),
                               PAPER1_INPUT_ROOT=str(root/'empty_external_inputs'))
            completed = subprocess.run([sys.executable,'-B','-m','scripts.gaia_neighbour.prepare_inputs'],
                cwd=Path(__file__).resolve().parents[1], env=environment, capture_output=True,text=True)
            self.assertEqual(completed.returncode,0,completed.stdout+'\n'+completed.stderr)
            prepared = pd.read_csv(root/'work/one_arcsecond/scoring_inputs.csv')
            self.assertEqual(prepared.scoring_route.tolist(),['resolved_original'])


class MorphologyTests(unittest.TestCase):
    @staticmethod
    def image_fixture():
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [51., 51.]
        wcs.wcs.crval = [30., 0.]
        wcs.wcs.cdelt = [-.1/3600, .1/3600]
        wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        y, x = np.indices((101, 101), dtype=float)
        image = np.exp(-((x-45)**2 + (y-50)**2)/18.) + np.exp(-((x-55)**2 + (y-50)**2)/18.)
        ra, dec = wcs.pixel_to_world_values([45., 55.], [50., 50.])
        return image, wcs.to_header(), ra, dec

    def test_binary_is_elongated_and_background_invariant(self):
        image, header, ra, dec = self.image_fixture()
        result = measure_image(image, header, ra[0], dec[0], ra[1], dec[1])
        offset = measure_image(image + 50., header, ra[0], dec[0], ra[1], dec[1])
        self.assertGreater(result['r_perp_sigma_ratio'], 1.05)
        self.assertAlmostEqual(result['r_perp_sigma_ratio'], offset['r_perp_sigma_ratio'], places=10)

    def test_missing_wcs_or_wrong_field_fails(self):
        image, header, ra, dec = self.image_fixture()
        with self.assertRaises(ValueError):
            measure_image(image, {}, ra[0], dec[0], ra[1], dec[1])
        with self.assertRaisesRegex(ValueError, 'inside'):
            measure_image(image, header, ra[0], dec[0], ra[1] + 1, dec[1])

    def test_missing_images_stop_the_stage_without_dropping_objects(self):
        candidates = pd.DataFrame({'qso_parent_index':[1], 'primary_ra':[30.], 'primary_dec':[0.],
                                   'secondary_ra':[30.0001], 'secondary_dec':[0.]})
        with tempfile.TemporaryDirectory() as directory:
            config = {'_directory':directory, 'morphology_images':{
                'path_template':'{qso_parent_index}.fits'}}
            with self.assertRaisesRegex(FileNotFoundError, 'No objects were silently excluded'):
                measure_candidates(candidates, config)

    def test_supplied_fits_image_produces_morphology_and_provenance(self):
        image, header, ra, dec = self.image_fixture()
        candidates = pd.DataFrame({'qso_parent_index':[1], 'primary_ra':[ra[0]], 'primary_dec':[dec[0]],
                                   'secondary_ra':[ra[1]], 'secondary_dec':[dec[1]]})
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'1.fits'
            fits.writeto(path,image,header)
            config = {'_directory':directory, 'morphology_images':{'path_template':'{qso_parent_index}.fits'}}
            measured = measure_candidates(candidates,config)
            self.assertEqual(measured.qso_parent_index.tolist(),[1])
            self.assertEqual(measured.morphology_measurement_status.tolist(),['measured'])
            self.assertGreater(measured.r_perp_sigma_ratio.iloc[0],1.05)
            self.assertEqual(measured.morphology_image_path.iloc[0],str(path))


class OtherQuasarCatalogueTests(unittest.TestCase):
    def test_new_catalogue_maps_columns_and_discovers_without_fixed_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            quasars = pd.DataFrame({'ID':['0001','0002','0003'], 'RA':[30.,40.,50.],
                                    'DEC':[0.,0.,0.], 'Z':[1.,2.,3.], 'ZERR':[.1,.2,.3]})
            quasars.to_csv(root/'other.csv',index=False)
            gaia = pd.DataFrame({'source_id':['100000000000000001','100000000000000002',
                                            '200000000000000001','200000000000000002'],
                'ra':[30.,30.+.7/3600,40.,40.+.7/3600], 'dec':[0.]*4, 'phot_g_mean_mag':[18.]*4})
            gaia.to_csv(root/'gaia.csv',index=False)
            pd.DataFrame({'object_id':quasars.ID, 'ra_deg':quasars.RA,'dec_deg':quasars.DEC,
                          'query_radius_arcsec':[2.]*3}).to_csv(root/'coverage.csv',index=False)
            (root/'tractor.csv').write_text('unused\n')
            imaging = root/'imaging.json'
            imaging.write_text(json.dumps({'master_catalogue':'other.csv',
                'tractor_catalogues':['tractor.csv'],'tractor_coverage':'coverage.csv'}))
            config = {'_directory':directory,'catalogue':'other.csv','catalogue_label':'Future_QSO_Release',
                      'redshift_kind':'photometric','columns':{'object_id':'ID','ra':'RA','dec':'DEC',
                          'redshift':'Z','redshift_error':'ZERR'},
                      'gaia_catalogue':'gaia.csv','gaia_coverage':'coverage.csv'}
            summary = adapt_catalogue(config, root/'generated', imaging)
            self.assertEqual(summary['input_rows'],3)
            self.assertEqual(summary['unmatched_quasar_rows'],1)
            self.assertEqual(summary['imaging_gaia_anchors'],2)
            discovered = run_search(discovery_config(root/'generated/one_arcsecond_search.json'), root/'discovered')
            self.assertEqual(discovered['counts']['parent_entries'],2)
            parent = pd.read_csv(root/'discovered/discovered_parent.csv',dtype={'primary_qso_objid':'string'})
            self.assertEqual(parent.primary_qso_objid.tolist(), ['0001','0002'])
            self.assertTrue(parent.primary_qso_redshift_kind.eq('photometric').all())
            self.assertEqual(redshift_errors(parent,'primary_qso_zerr','primary_qso_catalog','primary_qso_redshift_kind').tolist(), [.1,.2])

    def test_explicit_redshift_kind_overrides_catalogue_name(self):
        frame = pd.DataFrame({'catalogue':['Future_QSOs','Quaia_spectroscopic_followup'],
                              'kind':['photometric','spectroscopic'],'error':[.2,.01]})
        self.assertEqual(redshift_errors(frame,'error','catalogue','kind').tolist(), [.2,0.])
