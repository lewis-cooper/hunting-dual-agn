"""Check identity-based tiers without changing candidate selection."""

import unittest

import pandas as pd

from dual_agn.candidate_tiers import assign_candidate_tiers


def example_candidates():
    return pd.DataFrame([
        dict(candidate_id='shared', primary_gaia_source_id='9007199254740993',
             secondary_gaia_source_id='9007199254740995', has_two_distinct_gaia_sources=True,
             selected_by_resolved=True, selected_by_unresolved=False,
             selected_by_imaging=True, imaging_selection_route='two_gaia',
             selected_by_both_searches=True, redshift=1., separation_arcsec=.5,
             projected_separation_kpc=4., literature_confirmed_dual=True,
             secondary_pm_missing=True),
        dict(candidate_id='additional', primary_gaia_source_id='9007199254740997',
             secondary_gaia_source_id='', has_two_distinct_gaia_sources=False,
             selected_by_resolved=False, selected_by_unresolved=False,
             selected_by_imaging=True, imaging_selection_route='one_gaia',
             selected_by_both_searches=False, redshift=2., separation_arcsec=1.,
             projected_separation_kpc=8., literature_confirmed_dual=False,
             secondary_pm_missing=False),
    ])


class CandidateTierTests(unittest.TestCase):
    def test_tiers_use_exact_ids_without_requiring_proper_motion(self):
        candidates = example_candidates()
        before = candidates.copy(deep=True)
        result = assign_candidate_tiers(candidates)
        self.assertEqual(result.tier.tolist(), ['main', 'additional'])
        pd.testing.assert_frame_equal(candidates, before)
        pd.testing.assert_frame_equal(result.drop(columns='tier'), before)

    def test_csv_flags_and_swapped_id_order_give_the_same_tiers(self):
        candidates = example_candidates().astype(str)
        candidates.loc[0, ['primary_gaia_source_id', 'secondary_gaia_source_id']] = ['9007199254740995', '9007199254740993']
        self.assertEqual(assign_candidate_tiers(candidates).tier.tolist(), ['main', 'additional'])

    def test_inexact_identifier_and_inconsistent_multiplicity_are_rejected(self):
        for column, value in [('primary_gaia_source_id', '9007199254740992.0'),
                              ('has_two_distinct_gaia_sources', False)]:
            with self.subTest(column=column):
                candidates = example_candidates()
                candidates.loc[0, column] = value
                with self.assertRaises(ValueError):
                    assign_candidate_tiers(candidates)



if __name__ == '__main__':
    unittest.main()
