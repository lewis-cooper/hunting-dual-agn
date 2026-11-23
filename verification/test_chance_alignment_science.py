"""Scientific regression tests carried over from the control experiment."""
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import sys
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'chance_alignments/engine'))
from native_photometry import additional_photometry
from build_parent_control_trials import resolve_photometric_area
from prepare_matched_parent_controls import allocate_quotas
from score_controls import match_catalogues
from summarize_repaired_controls import estimate_repaired_rate
from build_repaired_blended_controls import augment_background_fluxes
from cache_blend_priors import CachedMagnitudePrior
from summarize_parent_densities import density_interval
from score_parent_blended_controls import selected_trials

class ScientificRegressionTests(unittest.TestCase):
    def test_quotas_preserve_population_and_sample_totals(self):
        inputs=pd.DataFrame(dict(route=['resolved','blended'],primary_cell=['a','b']))
        quotas=allocate_quotas({'a':100,'b':10,'c':2},inputs,30,minimum=2)
        self.assertEqual(quotas.sum(),30)
        self.assertTrue((quotas<=pd.Series({'a':100,'b':10,'c':2})).all())
        self.assertEqual(quotas['c'],2)

    def test_legacy_area_boundary_respects_galactic_hemisphere(self):
        sources=pd.DataFrame(dict(ra=[180.,180.,0.,0.],dec=[40.]*4,
                                 release=[9010,9011,9010,9011]))
        selected=resolve_photometric_area(sources)
        self.assertEqual(selected.index.tolist(),[1,2])

    def test_same_gaia_in_different_fields_is_not_ambiguous(self):
        tractor=pd.DataFrame(dict(target_id=['a','b'],release=[9010]*2,brickid=[1]*2,
            objid=[1]*2,ra=[180.]*2,dec=[0.]*2,flux_r=[1.]*2))
        gaia=pd.DataFrame(dict(target_id=['a','b'],source_id=['123']*2,ra=[180.]*2,dec=[0.]*2,
            phot_g_mean_mag=[20.]*2,pmra=[0.]*2,pmdec=[0.]*2,pmra_error=[1.]*2,
            pmdec_error=[1.]*2,pmra_pmdec_corr=[0.]*2))
        result=match_catalogues(tractor,gaia)
        self.assertTrue(result.has_gaia.all())
        self.assertFalse(result.gaia_match_ambiguous.any())

    def test_low_snr_component_photometry_is_kept_for_mock_sum(self):
        row={'primaryDetection':1}
        for band in 'grizy':
            row.update({f'{band}PSFMag':22.,f'{band}PSFMagErr':1.,
                        f'{band}KronMag':22.,f'{band}KronMagErr':1.})
        with TemporaryDirectory() as directory,patch('native_photometry.fetch_ps1',return_value=pd.DataFrame([row])):
            strict,_=additional_photometry(180.,0.,Path(directory),surveys=['ps1'])
            raw,_=additional_photometry(180.,0.,Path(directory),component_quality=False,surveys=['ps1'])
        self.assertEqual(len(strict),0)
        self.assertEqual(len(raw),5)
        self.assertFalse(raw[0]['component_band_quality'])
        self.assertGreater(raw[0]['quality_variance'],0)

    def test_ps1_negative_flux_is_measured_noise_not_a_missing_magnitude(self):
        row={'primaryDetection':1}
        for band in 'grizy':
            row.update({f'{band}PSFMag':-999.,f'{band}PSFMagErr':-999.,
                f'{band}KronMag':-999.,f'{band}KronMagErr':-999.,
                f'{band}KronFlux':-1e-6,f'{band}KronFluxErr':2e-6,
                f'{band}PSFFlux':-1e-6,f'{band}PSFFluxErr':2e-6})
        with TemporaryDirectory() as directory,patch('native_photometry.fetch_ps1',return_value=pd.DataFrame([row])):
            raw,_=additional_photometry(180.,0.,Path(directory),component_quality=False,surveys=['ps1'])
        self.assertEqual(len(raw),5)
        self.assertAlmostEqual(raw[0]['flux'],-1e-6/(3631e-9))
        self.assertGreater(raw[0]['variance'],0)

    def test_expanding_parent_quotas_keeps_previous_fields_available(self):
        inputs=pd.DataFrame(dict(route=['resolved','blended'],primary_cell=['a','b']))
        previous=allocate_quotas({'a':100,'b':10,'c':2},inputs,30,minimum=2)
        expanded=allocate_quotas({'a':100,'b':10,'c':2},inputs,50,minimum=2,previous=previous)
        self.assertEqual(expanded.sum(),50)
        self.assertTrue((expanded>=previous).all())

    def test_repairing_one_source_many_times_does_not_improve_precision(self):
        controls=pd.DataFrame(dict(input_pair_id=['input']*2,mock_weight=[.5,.5],pair_factor=[1.,1.],
            diagnostic_only=[False]*2,individual_cuts_pass=[True,False],
            background_field_id=['field_a','field_b'],background_source_id=['a','b'],
            selection_mode=['enumerated']*2))
        first=estimate_repaired_rate(controls,1,500,np.random.default_rng(7))
        repeated=pd.concat([controls]*100,ignore_index=True)
        repeated['mock_weight']/=100
        second=estimate_repaired_rate(repeated,1,500,np.random.default_rng(7))
        self.assertAlmostEqual(first['f_pass_known'],second['f_pass_known'])
        self.assertAlmostEqual(first['f_pass_lower_95'],second['f_pass_lower_95'])
        self.assertAlmostEqual(first['f_pass_upper_95'],second['f_pass_upper_95'])
        self.assertEqual(second['N_background_sources'],2)
        self.assertAlmostEqual(second['effective_fields'],2)

    def test_repaired_unknown_outcomes_remain_in_bounds(self):
        controls=pd.DataFrame(dict(input_pair_id=['input']*2,mock_weight=[.5,.5],pair_factor=[.8,.8],
            diagnostic_only=[False,True],individual_cuts_pass=[True,False],
            background_field_id=['field_a','field_b'],background_source_id=['a','b'],
            selection_mode=['sampled']*2))
        result=estimate_repaired_rate(controls,1,100,np.random.default_rng(7))
        self.assertAlmostEqual(result['f_pass_population_lower'],.4)
        self.assertAlmostEqual(result['f_pass_population_upper'],.8)
        self.assertAlmostEqual(result['known_weight_fraction'],.5)
        controls.loc[1,'individual_cuts_pass']=True
        diagnostic=estimate_repaired_rate(controls,1,100,np.random.default_rng(7))
        self.assertAlmostEqual(diagnostic['f_pass_population_lower'],.4)
        self.assertAlmostEqual(diagnostic['f_pass_population_upper'],.8)

    def test_missing_native_flux_repair_rejects_shared_measurements(self):
        background=pd.DataFrame(dict(background_source_id=['a','b'],photometric_system=[None,None],
                                     flux_r=[np.nan,np.nan]))
        native=pd.DataFrame(dict(target_id=['a','b'],unique_gaia_association=[True,False],
            unique_native_association=[True,True],photometric_system=['decals_dr9_south']*2,flux_r=[1.,2.]))
        with TemporaryDirectory() as directory:
            path=Path(directory)/'native.parquet';native.to_parquet(path)
            repaired=augment_background_fluxes(background,path)
        self.assertEqual(repaired.loc[0,'flux_r'],1.)
        self.assertTrue(pd.isna(repaired.loc[1,'flux_r']))

    def test_population_prior_cache_keeps_exact_values_and_distinct_coordinates(self):
        class Prior:
            mag_edges=np.array([0.,30.])
            calls=0
            def __call__(self,z,magnitude):
                self.calls+=1
                return np.asarray(z)*np.exp(-magnitude)
        prior=Prior();cached=CachedMagnitudePrior(prior,'quasar')
        first=cached(np.array([1.1]),20.)
        np.testing.assert_array_equal(first,np.array([1.1])*np.exp(-20.))
        np.testing.assert_array_equal(first,cached(np.array([1.1]),20.))
        self.assertEqual(prior.calls,1)
        self.assertNotEqual(cached(np.array([1.2]),20.)[0],first[0])

    def test_zero_density_events_do_not_imply_zero_upper_limit(self):
        fields=pd.DataFrame(dict(sampling_weight=[10.,10.],sampling_cell=['a','a']))
        point,lower,upper=density_interval(fields,[0,0],1/3500,100,np.random.default_rng(3))
        self.assertEqual(point,0.)
        self.assertEqual(lower,0.)
        self.assertGreater(upper,0.)

    def test_fixed_mock_subset_preserves_one_unit_of_weight_per_real_input(self):
        frame=pd.DataFrame(dict(trial_id=[f'blended:3-10:input_a:{i}' for i in range(100)]+['blended:3-10:input_b:unmatched'],
            input_pair_id=['a']*100+['b'],annulus=['3-10']*101,mock_weight=[.01]*100+[1.]))
        with TemporaryDirectory() as directory:
            root=Path(directory);frame.to_parquet(root/'blended_control_trials.parquet')
            subset=selected_trials(SimpleNamespace(output=root,draws_per_input=20))
        self.assertEqual(len(subset),21)
        np.testing.assert_allclose(subset.groupby('input_pair_id').mock_weight.sum(),1.)

    def test_unequal_draw_counts_do_not_upweight_legacy_inputs(self):
        records=[]
        for identity,policy in [('a','legacy'),('b','no_legacy')]:
            records.extend(dict(trial_id=f'blended:3-10:{identity}:{i}',input_pair_id=identity,
                annulus='3-10',mock_weight=.01,photometry_class=policy) for i in range(100))
        with TemporaryDirectory() as directory:
            root=Path(directory);pd.DataFrame(records).to_parquet(root/'blended_control_trials.parquet')
            subset=selected_trials(SimpleNamespace(output=root,draws_per_input=20,legacy_draws_per_input=100))
        self.assertEqual(subset.groupby('input_pair_id').size().to_dict(),{'a':100,'b':20})
        np.testing.assert_allclose(subset.groupby('input_pair_id').mock_weight.sum(),1.)

    def test_batched_prior_cache_preserves_interpolation_and_support_boundaries(self):
        from scipy.interpolate import RegularGridInterpolator
        class GridQSOPrior:
            mag_edges=np.array([16.,24.]);mag_centres=np.array([17.,18.5,20.,23.])
            z_edges=np.array([.15,4.35]);z_centres=np.array([.25,1.,2.,4.])
            sigma=np.array([[2.,3.,1.,8.],[5.,1.,7.,4.],[9.,2.,5.,3.],[1.,8.,2.,6.]])
            def __init__(self):
                self.interpolator=RegularGridInterpolator((self.z_centres,self.mag_centres),self.sigma)
            def __call__(self,redshift,magnitude):
                redshift=np.atleast_1d(redshift)
                positions=np.column_stack([np.clip(redshift,self.z_centres[0],self.z_centres[-1]),
                    np.full(len(redshift),np.clip(magnitude,self.mag_centres[0],self.mag_centres[-1]))])
                valid=(redshift>=self.z_edges[0])&(redshift<=self.z_edges[-1])
                valid &= self.mag_edges[0]<=magnitude<=self.mag_edges[-1]
                return np.where(valid,self.interpolator(positions),0.)
        prior=GridQSOPrior();cached=CachedMagnitudePrior(prior,'quasar')
        magnitudes=[*np.arange(16.,24.+1e-9,.01)[::113],15.9,18.12345,24.,24.1]
        for redshift in [.149,.15,.5,1.,2.7,4.35,4.351]:
            for magnitude in magnitudes:
                np.testing.assert_allclose(cached(np.array([redshift]),magnitude),
                    prior(np.array([redshift]),magnitude),rtol=0,atol=1e-12)

