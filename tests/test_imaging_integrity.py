"""Reject malformed imaging pairs and enforce quality cuts despite good scores."""

from pathlib import Path
import json
import runpy
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from dual_agn import settings
from dual_agn.imaging import prepare_component_photometry, validate_component_photometry


def imaging_fixture():
    candidates = pd.DataFrame([dict(qso_parent_index=1, physical_system_key='pair:1',
        candidate_tier='one_gaia', primary_ra=180., primary_dec=30.,
        secondary_ra=180.0002, secondary_dec=30., adopted_redshift=1.1,
        adopted_redshift_error=.1, adopted_redshift_catalogue='quaia',
        redshift_error_source_id='9007199254740993', redshift_kind='photometric')])
    components, tractor = [], []
    for index, role in enumerate(('primary', 'secondary')):
        ra = candidates[f'{role}_ra'].iloc[0]
        components.append(dict(qso_parent_index=1, component_role=role,
            component_ra=ra, component_dec=30., tractor_release=9010,
            tractor_brickid=1, tractor_objid=index, z=1.1,
            tractor_flux_g=100., tractor_flux_r=100., tractor_flux_z=100.))
        row = dict(release=9010, brickid=1, objid=index, type='PSF',
                   ra=ra, dec=30., maskbits=0, fracflux_r=.1)
        for band in 'grz':
            row.update({f'flux_{band}':100., f'flux_ivar_{band}':1., f'nobs_{band}':2})
        tractor.append(row)
    return candidates, pd.DataFrame(components), pd.DataFrame(tractor)


class ImagingIntegrityTests(unittest.TestCase):
    def test_valid_measurements_preserve_roles_and_flux(self):
        candidates, components, tractor = imaging_fixture()
        prepared = prepare_component_photometry(candidates, components, tractor)
        self.assertEqual(prepared.component_role.tolist(), ['primary', 'secondary'])
        self.assertEqual(prepared.flux_r.tolist(), [100., 100.])
        validate_component_photometry(candidates, prepared)

    def test_missing_or_extra_component_is_an_error(self):
        candidates, components, tractor = imaging_fixture()
        for rows in (components.iloc[:1], components.iloc[:0],
                     pd.concat([components, components.iloc[:1]], ignore_index=True)):
            with self.subTest(count=len(rows)), self.assertRaisesRegex(ValueError, 'Exactly two'):
                prepare_component_photometry(candidates, rows, tractor)

    def test_unmatched_or_ambiguous_coordinates_are_an_error(self):
        for ambiguous in (False, True):
            candidates, components, tractor = imaging_fixture()
            if ambiguous:
                candidates['secondary_ra'] = candidates.primary_ra
            else:
                components.loc[1, 'component_ra'] += .01
            with self.subTest(ambiguous=ambiguous), self.assertRaisesRegex(ValueError, 'uniquely match'):
                prepare_component_photometry(candidates, components, tractor)

    def test_two_components_at_primary_position_are_an_error(self):
        candidates, components, tractor = imaging_fixture()
        components.loc[1, 'component_ra'] = components.loc[0, 'component_ra']
        with self.assertRaisesRegex(ValueError, 'unique primary and secondary'):
            prepare_component_photometry(candidates, components, tractor)

    def test_repeated_tractor_identity_is_an_error(self):
        candidates, components, tractor = imaging_fixture()
        components.loc[1, 'tractor_objid'] = components.loc[0, 'tractor_objid']
        with self.assertRaisesRegex(ValueError, 'Tractor identities'):
            prepare_component_photometry(candidates, components, tractor)

    def test_missing_exact_tractor_measurement_is_an_error(self):
        candidates, components, tractor = imaging_fixture()
        with self.assertRaisesRegex(ValueError, 'Missing exact Tractor'):
            prepare_component_photometry(candidates, components, tractor.iloc[:1])

    def test_contradictory_adopted_redshift_metadata_is_an_error(self):
        for field, value in [('z', 2.), ('adopted_redshift', 2.),
                             ('adopted_redshift_error', .3), ('redshift_kind', 'spectroscopic')]:
            candidates, components, tractor = imaging_fixture()
            if field not in components:
                components[field] = candidates[field].iloc[0]
            components.loc[1, field] = value
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, 'Inconsistent'):
                prepare_component_photometry(candidates, components, tractor)

    def test_prepared_component_roles_cannot_be_swapped(self):
        candidates, components, tractor = imaging_fixture()
        prepared = prepare_component_photometry(candidates, components, tractor)
        prepared['component_role'] = ['secondary', 'primary']
        with self.assertRaisesRegex(ValueError, 'coordinate roles'):
            validate_component_photometry(candidates, prepared)

    def run_selection(self, component_changes):
        candidates, components, tractor = imaging_fixture()
        prepared = prepare_component_photometry(candidates, components, tractor)
        for index, column, value in component_changes:
            prepared.loc[index, column] = value
        # Deliberately stale cached flags must not override actual bad measurements.
        for column in ('reference_snr_pass', 'latitude_pass', 'psf_pass'):
            prepared[column] = True
        scores = candidates.assign(primary_gaia_source_id='9007199254740993',
            secondary_gaia_source_id='-1', secondary_pm_sig=np.nan,
            both_components_eligible_magdep=True, both_components_eligible_magindep=True)
        for model in ('magdep', 'magindep'):
            scores[f'primary_log_r_per_unit_z_{model}'] = -1.
            scores[f'secondary_log_r_per_unit_z_{model}'] = -1.
        with TemporaryDirectory() as directory:
            root = Path(directory)
            prepared.to_csv(root/'component_photometry.csv', index=False)
            scores.to_csv(root/'system_scores.csv', index=False)
            pd.DataFrame([dict(qso_parent_index=1, image_axis_sigma_arcsec=1.2,
                image_perp_sigma_arcsec=1., r_perp_sigma_ratio=1.2,
                morphology_measurement_source='fixture')]).to_csv(root/'morphology.csv', index=False)
            (root/'calibration.json').write_text(json.dumps(dict(model='either', adopted_log_r_threshold=-6.81)))
            (root/'upstream_counts.json').write_text(json.dumps(dict(initial_psf_pairs=1,
                validated_pairs=1, primary_gaia=1)))
            pd.DataFrame([dict(stage=stage, tier=tier, systems=int(tier=='one_gaia'))
                for stage in ('Gaia tier','Primary PM','Delta(g-r)','Tier cut','Deduplicated pre-colour')
                for tier in ('one_gaia','two_gaia')]).to_csv(root/'upstream_cut_flow.csv', index=False)
            with patch.multiple(settings, IMAGING_ROOT=root, IMAGING_MORPHOLOGY_PATH=root/'morphology.csv',
                    THRESHOLD_CALIBRATION_PATH=root/'calibration.json',
                    IMAGING_SCORES_PATH=root/'scores.csv', IMAGING_CUTFLOW_PATH=root/'cutflow.csv'):
                runpy.run_module('scripts.imaging.select_candidates', run_name='__main__')
            return pd.read_csv(root/'scores.csv'), pd.read_csv(root/'cutflow.csv')

    def test_good_scores_do_not_override_either_components_quality_failure(self):
        for index in (0, 1):
            for column, value in [('r_snr',9.99), ('r_snr',np.nan),
                                  ('galactic_abs_b_deg',24.99), ('galactic_abs_b_deg',np.inf),
                                  ('tractor_type','EXP')]:
                with self.subTest(component=index, column=column, value=value):
                    scores, cutflow = self.run_selection([(index,column,value)])
                    self.assertFalse(scores.joint_quality_pass.any())
                    self.assertFalse(scores.joint_final_pass.any())
                    self.assertEqual(cutflow.iloc[-1].Total, 0)

    def test_quality_boundaries_and_good_scores_remain_selected(self):
        scores, cutflow = self.run_selection([(0,'r_snr',10.), (1,'galactic_abs_b_deg',25.)])
        self.assertTrue(scores.joint_quality_pass.all())
        self.assertTrue(scores.joint_final_pass.all())
        self.assertEqual(cutflow.iloc[-1].Total, 1)


if __name__ == '__main__':
    unittest.main()
