"""Check local pair discovery, survey resolution, coverage and Gaia/DCR geometry."""

import json
from pathlib import Path
import socket
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from astropy.table import Table
from astropy.wcs import WCS
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dual_agn.imaging_search.catalogues import (find_tractor_pairs, read_master, require_tractor_coverage,
                                      resolve_config, resolve_tractor)
from dual_agn.imaging_search.dcr import axis_delta, calculate_dcr, footprint_contains
from dual_agn.imaging_search.pairs import (build_candidate_summary, build_component_table,
                                 local_color_stats, match_gaia_components, pm_sig,
                                 require_gaia_coverage)
from dual_agn.imaging_search.pipeline import GENERATED_INPUTS, run_search
from scripts.run_analysis import commands


class ImagingSearchTests(unittest.TestCase):
    def master(self):
        return pd.DataFrame([dict(qso_parent_index=0,source_id='9007199254740993',
            gaia_ra=10.,gaia_dec=0.,q_ra=10.,q_dec=0.,z=1.,adopted_redshift=1.,
            adopted_redshift_catalogue='quaia',sdss_name='example')])

    def tractor(self):
        rows=[]
        for index,(ra,flux) in enumerate([(10.,100.),(10.0002,50.)]):
            rows.append(dict(release=9010,brickid=1,objid=index,brick_primary=True,
                type='PSF',ra=ra,dec=0.,flux_g=flux*.7,flux_r=flux,flux_z=flux*1.2,
                flux_ivar_g=100.,flux_ivar_r=100.,flux_ivar_z=100.,
                nobs_g=3,nobs_r=3,nobs_z=3,maskbits=0,fracflux_r=.1))
        return pd.DataFrame(rows)

    def coverage(self):
        return pd.DataFrame([dict(source_id='9007199254740993',ra_deg=10.,dec_deg=0.,query_radius_arcsec=10.)])

    def ccd(self):
        return pd.DataFrame([dict(qso_parent_index=0,target_ra=10.,target_dec=0.,
            filter='g',ccd_cuts=0,camera='decam',width=2048,height=2048,
            crpix1=1024.,crpix2=1024.,crval1=10.,crval2=0.,
            cd1_1=-.262/3600,cd1_2=0.,cd2_1=0.,cd2_2=.262/3600,
            mjd_obs=60000.,exptime=100.,sig1=1.)])

    def components(self):
        with TemporaryDirectory() as directory:
            path=Path(directory)/'tractor.csv'
            self.tractor().to_csv(path,index=False)
            channels,local,_,_=find_tractor_pairs(self.master(),[path],10.,1.5)
        return build_component_table(Table.from_pandas(channels),Table.from_pandas(local),1.5,0)

    def gaia(self):
        return pd.DataFrame([dict(source_id='9007199254740993',ra=10.0001,dec=0.,
            phot_g_mean_mag=18.,pmra=0.,pmdec=0.,pmra_error=1.,pmdec_error=1.,
            pmra_pmdec_corr=0.,parallax_over_error=0.,ruwe=1.)])

    def fixture_config(self, root):
        self.master().to_csv(root/'master.csv',index=False)
        self.tractor().to_csv(root/'tractor.csv',index=False)
        self.coverage().to_csv(root/'tractor_coverage.csv',index=False)
        self.gaia().to_csv(root/'gaia.csv',index=False)
        pd.DataFrame({'ra_deg':[10.,10.0002],'dec_deg':[0.,0.],
            'query_radius_arcsec':[.7,.7]}).to_csv(root/'gaia_coverage.csv',index=False)
        self.ccd().to_csv(root/'ccds.csv',index=False)
        pd.DataFrame([dict(qso_parent_index=0,ra_deg=10.,dec_deg=0.,
            query_radius_deg=.65)]).to_csv(root/'ccd_coverage.csv',index=False)
        pd.DataFrame([dict(source_id='9007199254740993',redshift_quaia=1.,
            redshift_quaia_err=.1)]).to_csv(root/'quaia.csv',index=False)
        config=dict(master_catalogue='master.csv',tractor_catalogues=['tractor.csv'],
            tractor_coverage='tractor_coverage.csv',gaia_catalogue='gaia.csv',
            gaia_coverage='gaia_coverage.csv',ccd_catalogue='ccds.csv',
            ccd_coverage='ccd_coverage.csv',quaia_catalogue='quaia.csv')
        path=root/'config.json'
        path.write_text(json.dumps(config))
        return resolve_config(path)

    def test_two_sources_are_found_from_unassociated_catalogue_rows(self):
        with TemporaryDirectory() as directory:
            path=Path(directory)/'tractor.csv'
            self.tractor().to_csv(path,index=False)
            channel,local,reference,summary=find_tractor_pairs(self.master(),[path],10.,1.5)
        self.assertEqual(len(channel),1)
        self.assertEqual(channel.type_pair.iloc[0],'PSF-PSF')
        self.assertAlmostEqual(channel.component_sep_arcsec.iloc[0],.72,places=7)
        self.assertEqual(len(local),2)
        self.assertEqual(summary['initial_psf_pairs'],1)

    def test_parent_cache_rekeys_by_exact_gaia_id(self):
        with TemporaryDirectory() as directory:
            path=Path(directory)/'tractor.csv'
            raw=self.tractor().assign(qso_parent_index=999,parent_gaia_source_id='9007199254740993')
            raw.to_csv(path,index=False)
            channel,_,_,_=find_tractor_pairs(self.master(),[path],10.,1.5)
        self.assertEqual(channel.qso_parent_index.iloc[0],0)

    def test_master_retains_original_redshift_parsing_and_exact_identifiers(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'master.csv'
            self.master().assign(adopted_redshift=.43030370980814325).to_csv(path,index=False)
            expected = pd.read_csv(path).adopted_redshift.iloc[0]
            exact = pd.read_csv(path,float_precision='round_trip').adopted_redshift.iloc[0]
            actual = read_master(path)
        self.assertNotEqual(expected,exact)
        self.assertEqual(actual.adopted_redshift.iloc[0],expected)
        self.assertEqual(actual.source_id.iloc[0],'9007199254740993')

    def test_missing_query_coverage_is_not_a_non_detection(self):
        coverage=self.coverage().assign(query_radius_arcsec=1.)
        with self.assertRaisesRegex(ValueError,'lack valid, complete'):
            require_tractor_coverage(self.master(),coverage,10.)

    def test_wrong_query_position_fails_coverage_validation(self):
        with self.assertRaisesRegex(ValueError,'lack valid, complete'):
            require_tractor_coverage(self.master(),self.coverage().assign(ra_deg=11.),10.)

    def test_north_south_overlap_uses_resolved_dr9_area(self):
        rows=pd.concat([self.tractor().iloc[:1]]*2,ignore_index=True)
        rows['ra']=180.;rows['dec']=50.;rows['release']=[9010,9011]
        resolved=resolve_tractor(rows)
        self.assertEqual(resolved.release.tolist(),[9011])

    def test_same_gaia_detection_is_not_forced_into_a_two_gaia_pair(self):
        components=self.components()
        matches=match_gaia_components(components,self.gaia(),.7)
        summary=build_candidate_summary(components,matches)
        self.assertTrue(summary['primary_has_gaia'][0])
        self.assertTrue(summary['secondary_has_gaia'][0])
        self.assertFalse(summary['secondary_has_distinct_gaia'][0])

    def test_gaia_non_detection_needs_complete_component_coverage(self):
        components=self.components()
        coverage=pd.DataFrame([dict(ra_deg=11.,dec_deg=0.,query_radius_arcsec=.7)])
        with self.assertRaisesRegex(ValueError,'lack complete Gaia coverage'):
            require_gaia_coverage(components,coverage,.7)

    def test_covariance_is_used_in_pm_significance(self):
        actual=pm_sig([1.],[1.],[1.],[1.],[.5])
        self.assertAlmostEqual(actual[0],np.sqrt(4/3),places=12)

    def test_colour_reference_excludes_the_query_system(self):
        ids=np.arange(25)
        colours=np.linspace(-1,1,25)
        median,_,_,_,count=local_color_stats(np.array([0]),np.array([1.]),ids,np.ones(25),colours,.1)
        self.assertEqual(count[0],24)
        self.assertAlmostEqual(median[0],np.median(colours[1:]))

    def test_vectorized_ccd_footprint_matches_astropy_wcs(self):
        ccds=pd.concat([self.ccd()]*4,ignore_index=True)
        ras=np.array([10.,10.03,9.97,11.])
        decs=np.array([0.,.02,-.02,0.])
        actual=footprint_contains(ccds,ras,decs)
        header={'CTYPE1':'RA---TAN','CTYPE2':'DEC--TAN'}
        for key in ['crpix1','crpix2','crval1','crval2','cd1_1','cd1_2','cd2_1','cd2_2']:
            header[key.upper()]=ccds[key].iloc[0]
        x,y=WCS(header).world_to_pixel_values(ras,decs)
        np.testing.assert_array_equal(actual,(x>=0)&(x<2048)&(y>=0)&(y<2048))

    def test_dcr_average_uses_an_axis_not_a_linear_angle(self):
        target=Table(dict(qso_parent_index=[0],primary_ra=[10.],primary_dec=[0.],secondary_ra=[10.0002],secondary_dec=[0.]))
        ccds=pd.concat([self.ccd()]*2,ignore_index=True)
        coverage=pd.DataFrame([dict(qso_parent_index=0,ra_deg=10.,dec_deg=0.,query_radius_deg=.65)])
        with patch('dual_agn.imaging_search.dcr.parallactic_pa',return_value=np.array([179.,1.])):
            result=calculate_dcr(target,ccds,coverage)
        self.assertLess(axis_delta(result.dcr_pa_deg.iloc[0],0.),1e-10)

    def test_end_to_end_search_is_offline_and_writes_vetting_inputs(self):
        with TemporaryDirectory() as directory:
            root=Path(directory)
            config=self.fixture_config(root)
            with patch.object(socket.socket,'connect',side_effect=AssertionError('Network forbidden')):
                result=run_search(config,root/'work')
            self.assertFalse(result['network_access'])
            self.assertEqual(result['counts']['validated_psf_pairs'],1)
            self.assertTrue(all((root/'work/inputs'/name).is_file() for name in GENERATED_INPUTS))

    def test_acquisition_discovery_is_reused_and_matches_a_fresh_search(self):
        with TemporaryDirectory() as directory:
            root=Path(directory)
            config=self.fixture_config(root)
            with patch('dual_agn.imaging_search.pipeline.find_tractor_pairs', wraps=find_tractor_pairs) as finder:
                run_search(config,root/'work',stage='discover')
                reused=run_search(config,root/'work',resume=True)
                self.assertEqual(finder.call_count,1)
            fresh=run_search(config,root/'fresh')
            self.assertEqual(reused['counts'],fresh['counts'])
            self.assertEqual(reused['output_sha256'],fresh['output_sha256'])

    def test_changed_gaia_rebuilds_associations_without_repeating_tractor_search(self):
        with TemporaryDirectory() as directory:
            root=Path(directory)
            config=self.fixture_config(root)
            with patch('dual_agn.imaging_search.pipeline.find_tractor_pairs', wraps=find_tractor_pairs) as finder:
                initial=run_search(config,root/'work')
                self.gaia().assign(pmra=2.).to_csv(root/'gaia.csv',index=False)
                updated=run_search(config,root/'work',resume=True)
                self.assertEqual(finder.call_count,1)
            self.assertNotEqual(initial['fingerprint'],updated['fingerprint'])
            self.assertNotEqual(initial['output_sha256']['inputs/imaging_labelled.fits'],
                                updated['output_sha256']['inputs/imaging_labelled.fits'])

    def test_changed_tractor_or_damaged_discovery_outputs_repeat_discovery(self):
        with TemporaryDirectory() as directory:
            root=Path(directory)
            config=self.fixture_config(root)
            with patch('dual_agn.imaging_search.pipeline.find_tractor_pairs', wraps=find_tractor_pairs) as finder:
                run_search(config,root/'work')
                self.tractor().assign(flux_r=[120.,50.]).to_csv(root/'tractor.csv',index=False)
                run_search(config,root/'work',resume=True)
                self.assertEqual(finder.call_count,2)
                (root/'work/inputs/imaging_channel.fits').write_bytes(b'incomplete output')
                run_search(config,root/'work',resume=True)
                self.assertEqual(finder.call_count,3)
            self.assertEqual(len(Table.read(root/'work/inputs/imaging_channel.fits')),1)

    def test_full_runner_searches_before_prepared_input_validation(self):
        sequence=commands('vetting',imaging_search_config=Path('search.json'),data_source='local')
        search=next(i for i,command in enumerate(sequence) if 'scripts.imaging.search_candidates' in command)
        prepared=next(i for i,command in enumerate(sequence)
                      if 'verification/verify_input_contracts.py' in command and '--external-only' not in command)
        self.assertLess(search,prepared)


if __name__=='__main__':
    unittest.main()
