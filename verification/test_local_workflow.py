"""Check local discovery and explicit failure on missing observational inputs."""

import json
from pathlib import Path
import socket
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from astropy.io import fits

from run_analysis import commands
from scripts.shared.check_local_inputs import check_inputs
from scripts.shared.prepare_quasar_catalogue import adapt_catalogue
from scripts.imaging.measure_morphology import measure_candidates
from dual_agn.one_arcsecond_search import resolve_config, run_search
from verification import test_native_measurements


class LocalWorkflowTests(unittest.TestCase):
    def test_default_plan_constructs_neighbour_counts_without_acquisition(self):
        sequence=commands('vetting')
        modules=[command[command.index('-m')+1] for command in sequence if '-m' in command]
        self.assertIn('scripts.shared.prepare_quasar_catalogue',modules)
        self.assertIn('scripts.gaia_neighbour.search_candidates',modules)
        self.assertLess(modules.index('scripts.shared.prepare_quasar_catalogue'),
                        modules.index('scripts.gaia_neighbour.search_candidates'))
        self.assertFalse(any('acquire' in name or 'wsdb' in name for name in modules))

    def test_raw_catalogues_determine_multiplicity_and_keep_isolated_inputs_in_audit(self):
        with TemporaryDirectory() as directory:
            root=Path(directory)
            pd.DataFrame(dict(id=['pair','isolated'],ra=[30.,40.],dec=[0.,0.],
                redshift=[1.,1.2])).to_csv(root/'quasars.csv',index=False)
            # Two companions are ordinary Gaia sources, absent from the quasar list.
            pd.DataFrame(dict(source_id=['100','101','102','200'],ra=[30.,30.+.3/3600,30.+.6/3600,40.],
                dec=[0.,0.,0.,0.],phot_g_mean_mag=[18.,19.,20.,18.])).to_csv(root/'gaia.csv',index=False)
            pd.DataFrame(dict(object_id=['pair','isolated'],ra_deg=[30.,40.],dec_deg=[0.,0.],
                query_radius_arcsec=[2.,2.])).to_csv(root/'coverage.csv',index=False)
            config=dict(_directory=directory,catalogue='quasars.csv',catalogue_label='test',
                redshift_kind='spectroscopic',gaia_catalogue='gaia.csv',gaia_coverage='coverage.csv',
                columns=dict(object_id='id',ra='ra',dec='dec',redshift='redshift'))
            with patch.object(socket.socket,'connect',side_effect=AssertionError('Network forbidden')):
                adapted=adapt_catalogue(config,root/'adapted')
                result=run_search(resolve_config(root/'adapted/one_arcsecond_search.json'),root/'search')
            parent=pd.read_csv(root/'search/discovered_parent.csv')
            self.assertEqual(adapted['input_rows'],2)
            self.assertEqual(result['counts']['parent_entries'],1)
            self.assertEqual(parent.gaia_edr3_cnt1.tolist(),[3])
            counts=pd.read_csv(root/'adapted/gaia_neighbour_counts.csv')
            self.assertEqual(counts.loc[counts.gaia_source_id.eq(200),'gaia_edr3_cnt1'].item(),1)

    def test_preflight_reports_missing_data_without_search_or_network(self):
        with TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'quasars.json').write_text(json.dumps(dict(catalogue='missing_quasars.csv',
                gaia_catalogue='missing_gaia.fits',gaia_coverage='missing_coverage.csv')))
            (root/'measurements.json').write_text(json.dumps(dict(gaia=dict(path='missing_gaia.fits',
                full_catalogue=True),tractor=dict(path='missing_tractor.fits',full_catalogue=True))))
            with patch.object(socket.socket,'connect',side_effect=AssertionError('Network forbidden')):
                with self.assertRaisesRegex(FileNotFoundError,'missing_gaia.fits'):
                    check_inputs(root/'quasars.json',root/'measurements.json',root/'unused.json','gaia_neighbour')

    def test_imaging_discovery_does_not_require_gaia_calibration_inputs(self):
        sequence=commands('discovery',pipeline='imaging')
        self.assertEqual(len(sequence),1)
        self.assertIn('scripts.imaging.search_candidates',sequence[0])

    def test_supplied_imaging_master_and_its_parent_ids_are_preserved(self):
        with TemporaryDirectory() as directory:
            root=Path(directory)
            pd.DataFrame(dict(id=['qso'],ra=[30.],dec=[0.],z=[1.])).to_csv(root/'quasars.csv',index=False)
            pd.DataFrame(dict(source_id=['100'],ra=[30.],dec=[0.],phot_g_mean_mag=[18.])).to_csv(root/'gaia.csv',index=False)
            pd.DataFrame(dict(object_id=['qso'],ra_deg=[30.],dec_deg=[0.],query_radius_arcsec=[2.])).to_csv(root/'coverage.csv',index=False)
            pd.DataFrame(dict(qso_parent_index=[99],source_id=['100'],gaia_ra=[30.],gaia_dec=[0.],
                adopted_redshift=[1.],adopted_redshift_catalogue=['test'])).to_csv(root/'imaging_master.csv',index=False)
            pd.DataFrame(columns=['release','brickid','objid']).to_csv(root/'tractor.csv',index=False)
            (root/'imaging.json').write_text(json.dumps(dict(master_catalogue='imaging_master.csv',
                tractor_catalogues=['tractor.csv'],tractor_coverage='coverage.csv')))
            config=dict(_directory=directory,catalogue='quasars.csv',catalogue_label='test',
                redshift_kind='spectroscopic',gaia_catalogue='gaia.csv',gaia_coverage='coverage.csv',
                columns=dict(object_id='id',ra='ra',dec='dec',redshift='z'))
            adapt_catalogue(config,root/'adapted',root/'imaging.json')
            generated=json.loads((root/'adapted/imaging_search.json').read_text())
            self.assertEqual(Path(generated['master_catalogue']),(root/'imaging_master.csv').resolve())
            self.assertEqual(pd.read_csv(generated['master_catalogue']).qso_parent_index.tolist(),[99])

    def test_morphology_requires_supplied_fits_and_has_no_download_option(self):
        candidates=pd.DataFrame(dict(qso_parent_index=[0],primary_ra=[30.],primary_dec=[0.],
            secondary_ra=[30.0002],secondary_dec=[0.]))
        with TemporaryDirectory() as directory:
            config=dict(_directory=directory,morphology_images=dict(path_template='{qso_parent_index}.fits'))
            with self.assertRaisesRegex(FileNotFoundError,'Supply the required local FITS'):
                measure_candidates(candidates,config)
            config['morphology_images']['download_missing']=True
            with self.assertRaisesRegex(ValueError,'Unsupported morphology_images settings'):
                measure_candidates(candidates,config)

    def test_wrong_field_image_stops_instead_of_becoming_a_morphology_rejection(self):
        image,header,_,_=test_native_measurements.MorphologyTests.image_fixture()
        candidates=pd.DataFrame(dict(qso_parent_index=[0],primary_ra=[60.],primary_dec=[0.],
            secondary_ra=[60.0002],secondary_dec=[0.]))
        with TemporaryDirectory() as directory:
            root=Path(directory)
            fits.PrimaryHDU(image,header).writeto(root/'0.fits')
            config=dict(_directory=directory,morphology_images=dict(path_template='{qso_parent_index}.fits'))
            with self.assertRaisesRegex(ValueError,'does not contain both components'):
                measure_candidates(candidates,config)
