"""Check the historical discovery order, exact identities and workflow connections."""

import json
from pathlib import Path
import socket
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dual_agn.one_arcsecond_search import (deduplicate_quasars, match_local_gaia, resolve_config,
    run_search, select_neighbour_matches, validate_measurement_identities)
from run_analysis import commands


class OneArcsecondDiscoveryTests(unittest.TestCase):
    def matches(self):
        return pd.DataFrame([
            dict(qso_catalog='SDSS_DR16Q', qso_uid='SDSS:one', qso_objid='1',
                 qso_name='first', qso_ra=10., qso_dec=0., qso_z=1., qso_zerr=.01,
                 gaia_source_id='9007199254740993', gaia_sep_arcsec=.1),
            dict(qso_catalog='SDSS_DR16Q', qso_uid='SDSS:one', qso_objid='1',
                 qso_name='first', qso_ra=10., qso_dec=0., qso_z=1., qso_zerr=.01,
                 gaia_source_id='9007199254740994', gaia_sep_arcsec=.8),
            dict(qso_catalog='DESI_DR1_QSO', qso_uid='DESI:one', qso_objid='2',
                 qso_name='second', qso_ra=10.00001, qso_dec=0., qso_z=1.1, qso_zerr=.02,
                 gaia_source_id='9007199254740994', gaia_sep_arcsec=.2),
            dict(qso_catalog='QUAIA_G20_V1', qso_uid='QUAIA:two', qso_objid='3',
                 qso_name='third', qso_ra=20., qso_dec=1., qso_z=2., qso_zerr=.2,
                 gaia_source_id='9007199254740995', gaia_sep_arcsec=0.),
        ])

    def counts(self):
        return pd.DataFrame(dict(gaia_source_id=['9007199254740993', '9007199254740994',
                                                '9007199254740995'], gaia_edr3_cnt1=[1, 2, 3]))

    def parent(self):
        _, nearest = select_neighbour_matches(self.matches(), self.counts(), 1.)
        return deduplicate_quasars(nearest, 1.)

    def test_filter_precedes_nearest_qualifying_match(self):
        selected, nearest = select_neighbour_matches(self.matches(), self.counts(), 1.)
        self.assertEqual(len(selected), 3)
        self.assertEqual(nearest.loc[nearest.qso_uid.eq('SDSS:one'), 'gaia_source_id'].item(), '9007199254740994')

    def test_triples_are_included_and_no_final_pm_or_latitude_cut_is_added(self):
        matches = self.matches().assign(secondary_pm_significance=1000., galactic_b=0.)
        _, nearest = select_neighbour_matches(matches, self.counts(), 1.)
        parent = deduplicate_quasars(nearest, 1.)
        self.assertEqual(parent.gaia_edr3_cnt1.tolist(), [2, 3])

    def test_nonprimary_desi_qso_records_are_not_excluded(self):
        matches = self.matches().loc[lambda frame: frame.qso_catalog.eq('DESI_DR1_QSO')]
        matches = matches.assign(qso_zcat_primary=False, qso_main_primary=False,
                                 qso_zwarn=0, qso_spectype='QSO')
        _, nearest = select_neighbour_matches(matches, self.counts(), 1.)
        parent = deduplicate_quasars(nearest, 1.)
        self.assertEqual(parent.primary_qso_objid.tolist(), ['2'])
        self.assertFalse(parent.qso_zcat_primary.any())

    def test_grouping_uses_survey_positions_even_when_gaia_positions_would_merge(self):
        matches = self.matches().iloc[[0, 2]].copy()
        matches['qso_ra'] = [10., 10. + 1.15/3600.]
        matches['gaia_ra'] = [10., 10. + .99/3600.]
        counts = self.counts().assign(gaia_edr3_cnt1=2)
        _, nearest = select_neighbour_matches(matches, counts, 1.)
        self.assertEqual(len(deduplicate_quasars(nearest, 1.)), 2)

    def test_cross_catalogue_duplicates_use_historical_primary_priority(self):
        parent = self.parent()
        self.assertEqual(len(parent), 2)
        self.assertEqual(parent.primary_qso_catalog.iloc[0], 'SDSS_DR16Q')
        self.assertEqual(parent.n_qso_catalog_rows.iloc[0], 2)
        self.assertEqual(parent.primary_qso_z.iloc[0], 1.)
        self.assertEqual(parent.qso_group_id.tolist(), ['QSO_CNT1_GT1_000001', 'QSO_CNT1_GT1_000002'])

    def test_ids_remain_exact_above_float_precision(self):
        parent = self.parent()
        self.assertEqual(parent.matched_gaia_source_id.tolist(), ['9007199254740994', '9007199254740995'])
        counts = self.counts().assign(gaia_source_id=['9007199254740993.0', '9007199254740994', '9007199254740995'])
        with self.assertRaisesRegex(ValueError, 'not floats'):
            select_neighbour_matches(self.matches(), counts, 1.)

    def test_missing_counts_are_not_treated_as_no_companion(self):
        with self.assertRaisesRegex(ValueError, 'neighbour-count coverage'):
            select_neighbour_matches(self.matches(), self.counts().iloc[1:], 1.)

    def test_prebuilt_parent_is_not_a_discovery_input(self):
        with self.assertRaisesRegex(ValueError, 'previously selected parent table'):
            select_neighbour_matches(self.matches().assign(qso_group_id='already_selected'), self.counts(), 1.)

    def test_duplicate_and_invalid_counts_fail(self):
        with self.assertRaisesRegex(ValueError, 'unique Gaia'):
            select_neighbour_matches(self.matches(), pd.concat([self.counts()]*2), 1.)
        for invalid in [np.nan, 0, 1.5]:
            counts = self.counts().astype({'gaia_edr3_cnt1': float})
            counts.loc[0, 'gaia_edr3_cnt1'] = invalid
            with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError, 'positive integer'):
                select_neighbour_matches(self.matches(), counts, 1.)

    def test_old_count_columns_cannot_override_supplied_counts(self):
        matches = self.matches().assign(gaia_edr3_cnt1=99)
        _, nearest = select_neighbour_matches(matches, self.counts(), 1.)
        self.assertEqual(set(nearest.gaia_edr3_cnt1), {2, 3})

    def test_spherical_duplicate_matching_wraps_ra(self):
        _, nearest = select_neighbour_matches(self.matches(), self.counts(), 1.)
        nearest.loc[nearest.qso_catalog.eq('SDSS_DR16Q'), 'qso_ra'] = 359.99995
        nearest.loc[nearest.qso_catalog.eq('DESI_DR1_QSO'), 'qso_ra'] = .00005
        parent = deduplicate_quasars(nearest, 1.)
        self.assertEqual(len(parent), 2)

    def test_local_match_keeps_all_gaia_matches_not_only_the_nearest(self):
        quasars = self.matches().iloc[:1].drop(columns=['gaia_source_id', 'gaia_sep_arcsec'])
        gaia = pd.DataFrame(dict(source_id=['9007199254740993', '9007199254740994'],
                                  ra=[10.00002, 10.0002], dec=[0., 0.]))
        coverage = pd.DataFrame(dict(qso_uid=['SDSS:one'], ra_deg=[10.], dec_deg=[0.], query_radius_arcsec=[1.]))
        matches = match_local_gaia(quasars, gaia, coverage, 1.)
        self.assertEqual(len(matches), 2)
        _, nearest = select_neighbour_matches(matches, self.counts(), 1.)
        self.assertEqual(nearest.gaia_source_id.item(), '9007199254740994')

    def test_missing_or_wrong_local_coverage_fails(self):
        quasars = self.matches().iloc[:1].drop(columns=['gaia_source_id', 'gaia_sep_arcsec'])
        gaia = pd.DataFrame(dict(source_id=['9007199254740993'], ra=[10.], dec=[0.]))
        for radius, ra in [(.5, 10.), (1., 11.)]:
            coverage = pd.DataFrame(dict(qso_uid=['SDSS:one'], ra_deg=[ra], dec_deg=[0.], query_radius_arcsec=[radius]))
            with self.subTest(radius=radius, ra=ra), self.assertRaisesRegex(ValueError, 'complete local Gaia coverage'):
                match_local_gaia(quasars, gaia, coverage, 1.)

    def test_quaia_uses_its_published_source_id(self):
        quasars = self.matches().iloc[3:].drop(columns=['gaia_source_id', 'gaia_sep_arcsec'])
        quasars['qso_parent_gaia_source_id'] = '9007199254740995'
        gaia = pd.DataFrame(dict(source_id=['9007199254740993'], ra=[20.], dec=[1.]))
        coverage = pd.DataFrame(columns=['qso_uid', 'ra_deg', 'dec_deg', 'query_radius_arcsec'])
        matches = match_local_gaia(quasars, gaia, coverage, 1.)
        self.assertEqual(matches.gaia_source_id.item(), '9007199254740995')

    def measurements(self):
        return self.parent()[['qso_group_id', 'primary_qso_catalog', 'primary_qso_ra',
                              'primary_qso_dec', 'matched_gaia_source_id']].rename(
                                  columns={'matched_gaia_source_id': 'primary_component_gaia_source_id'})

    def test_shifted_group_ids_cannot_attach_another_systems_measurements(self):
        measurements = self.measurements()
        validate_measurement_identities(self.parent(), measurements)
        measurements['primary_component_gaia_source_id'] = measurements.primary_component_gaia_source_id.iloc[::-1].to_numpy()
        with self.assertRaisesRegex(ValueError, 'identities disagree'):
            validate_measurement_identities(self.parent(), measurements)

    def test_changed_discovery_membership_requires_rebuilt_measurements(self):
        with self.assertRaisesRegex(ValueError, 'membership differs'):
            validate_measurement_identities(self.parent(), self.measurements().iloc[:1])

    def test_offline_search_and_resume_validate_output_hash(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            self.matches().to_csv(root/'matches.csv', index=False)
            self.counts().to_csv(root/'counts.csv', index=False)
            path = root/'config.json'
            path.write_text(json.dumps(dict(input_mode='survey_gaia_matches',
                quasar_gaia_matches=['matches.csv'], gaia_neighbour_counts='counts.csv')))
            config = resolve_config(path)
            with patch.object(socket.socket, 'connect', side_effect=AssertionError('Network forbidden')):
                summary = run_search(config, root/'work')
                resumed = run_search(config, root/'work', resume=True)
                self.assertEqual(summary, resumed)
                (root/'work/discovered_parent.csv').write_text('changed\n')
                rebuilt = run_search(config, root/'work', resume=True)
            self.assertEqual(summary['output_sha256'], rebuilt['output_sha256'])
            self.assertFalse(summary['network_access'])
            self.assertEqual(summary['counts']['parent_entries'], 2)

    def test_full_local_catalogue_search_includes_empty_covered_fields(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            quasars = self.matches().iloc[[0, 2, 3]].drop(columns=['gaia_source_id', 'gaia_sep_arcsec']).copy()
            no_match = quasars.iloc[:1].assign(qso_uid='SDSS:empty', qso_ra=30.)
            quasars = pd.concat([quasars, no_match], ignore_index=True)
            quasars['qso_parent_gaia_source_id'] = pd.Series([pd.NA]*len(quasars), dtype='string')
            quasars.loc[quasars.qso_catalog.eq('QUAIA_G20_V1'), 'qso_parent_gaia_source_id'] = '9007199254740995'
            quasars.to_csv(root/'quasars.csv', index=False)
            pd.DataFrame(dict(source_id=['9007199254740993', '9007199254740994'],
                              ra=[10.00002, 10.0002], dec=[0., 0.])).to_csv(root/'gaia.csv', index=False)
            coverage = quasars.loc[~quasars.qso_catalog.eq('QUAIA_G20_V1'), ['qso_uid', 'qso_ra', 'qso_dec']]
            coverage.rename(columns={'qso_ra': 'ra_deg', 'qso_dec': 'dec_deg'}).assign(
                query_radius_arcsec=1.).to_csv(root/'coverage.csv', index=False)
            self.counts().to_csv(root/'counts.csv', index=False)
            path = root/'config.json'
            path.write_text(json.dumps(dict(input_mode='local_catalogues',
                quasar_catalogues=['quasars.csv'], gaia_catalogue='gaia.csv',
                gaia_coverage='coverage.csv', gaia_neighbour_counts='counts.csv')))
            with patch.object(socket.socket, 'connect', side_effect=AssertionError('Network forbidden')):
                summary = run_search(resolve_config(path), root/'work')
            self.assertEqual(summary['counts']['quasar_catalogue_rows'], 4)
            self.assertEqual(summary['counts']['gaia_match_rows'], 5)
            self.assertEqual(summary['counts']['parent_entries'], 2)

    def test_full_workflow_discovers_both_samples_before_preparation(self):
        sequence = commands('vetting', quasar_catalogue_config=Path('quasars.json'))
        modules = [command[command.index('-m')+1] for command in sequence if '-m' in command]
        self.assertLess(modules.index('scripts.gaia_neighbour.search_candidates'),
                        modules.index('scripts.gaia_neighbour.prepare_inputs'))
        search=next(command for command in sequence if 'scripts.gaia_neighbour.search_candidates' in command)
        self.assertEqual(Path(search[search.index('--config')+1]).name,'one_arcsecond_search.json')
        discovery=commands('discovery',pipeline='gaia_neighbour')
        self.assertEqual(len(discovery),2)
        self.assertIn('scripts.shared.prepare_quasar_catalogue',discovery[0])
        self.assertIn('scripts.gaia_neighbour.search_candidates',discovery[1])


if __name__ == '__main__':
    unittest.main()
