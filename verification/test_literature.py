"""Exercise pair identity, uncertain classifications and exclusion safeguards."""

from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from literature import annotate_candidates, match_references, read_references
from scripts.run_analysis import VETTING_STAGES
from settings import LITERATURE_ROOT


class LiteratureTests(unittest.TestCase):
    def candidates(self):
        return pd.DataFrame([dict(candidate_id='pair1', system_name='target',
            primary_ra_deg=10., primary_dec_deg=0., secondary_ra_deg=10.0002,
            secondary_dec_deg=0., separation_arcsec=.72, redshift=2.,
            primary_gaia_source_id='9007199254740993',
            secondary_gaia_source_id='9007199254740995')])

    def reference(self, **changes):
        values = dict(reference_id='paper:pair1', system_id='published-pair',
            system_name='published target', ra_deg=10., dec_deg=0.,
            secondary_ra_deg=10.0002, secondary_dec_deg=0.,
            redshift=.5, separation_arcsec=.72, position_angle_deg=np.nan,
            primary_gaia_source_id='', secondary_gaia_source_id='',
            classification='confirmed_lens', author_year='Example (2026)',
            publication_url='https://example.org/paper', source_table='Table 1',
            reference_catalogue='known_lenses_contaminants', match_rule='component_pair',
            position_tolerance_arcsec=.1, separation_tolerance_arcsec=.02)
        values.update(changes)
        return pd.DataFrame([values])

    def audit(self, reference):
        candidates = self.candidates()
        associations = match_references(candidates, reference)
        return associations, annotate_candidates(candidates, associations)

    def test_swapped_published_components_identify_the_same_pair(self):
        associations, audit = self.audit(self.reference(ra_deg=10.0002, secondary_ra_deg=10.))
        self.assertTrue(associations.match_verified.all())
        self.assertTrue(audit.literature_contaminant_removed.all())

    def test_nearby_primary_does_not_verify_an_unrelated_companion(self):
        associations, audit = self.audit(self.reference(secondary_ra_deg=10.0004))
        self.assertEqual(len(associations), 1)
        self.assertFalse(associations.match_verified.any())
        self.assertFalse(audit.literature_contaminant_removed.any())
        self.assertTrue(audit.literature_possible_association.all())

    def test_one_shared_gaia_id_is_insufficient(self):
        associations, audit = self.audit(self.reference(secondary_ra_deg=10.0004,
            primary_gaia_source_id='9007199254740993', secondary_gaia_source_id='9007199254740997'))
        self.assertFalse(associations.match_verified.any())
        self.assertFalse(audit.literature_contaminant_removed.any())

    def test_exact_gaia_pair_preserves_large_integer_identity(self):
        associations, _ = self.audit(self.reference(ra_deg=10.00008, secondary_ra_deg=10.00028,
            primary_gaia_source_id='9007199254740995', secondary_gaia_source_id='9007199254740993'))
        self.assertEqual(associations.verification_method.iloc[0], 'exact_distinct_gaia_pair')

    def test_redshift_disagreement_does_not_veto_a_positional_pair(self):
        associations, audit = self.audit(self.reference(redshift=3.9))
        self.assertTrue(associations.match_verified.all())
        self.assertTrue(audit.literature_contaminant_removed.all())

    def test_a_confirmed_classification_without_pair_evidence_is_only_a_flag(self):
        _, audit = self.audit(self.reference(secondary_ra_deg=np.nan, secondary_dec_deg=np.nan,
                                            match_rule='position_only'))
        self.assertFalse(audit.literature_contaminant_removed.any())
        self.assertTrue(audit.literature_possible_association.all())

    def test_likely_lenses_and_other_uncertain_classes_are_retained(self):
        for classification in ['likely_lens', 'lens_candidate', 'possible_quasar_star',
                               'disputed_dual_lens', 'unclassified_quasar_pair']:
            with self.subTest(classification=classification):
                _, audit = self.audit(self.reference(classification=classification))
                self.assertFalse(audit.literature_contaminant_removed.any())
                self.assertTrue(audit.literature_uncertain_classification.all())

    def test_conflicting_confirmed_classes_are_retained_for_review(self):
        references = pd.concat([self.reference(), self.reference(reference_id='dual:pair1',
            classification='confirmed_dual_agn', reference_catalogue='known_binaries')], ignore_index=True)
        _, audit = self.audit(references)
        self.assertFalse(audit.literature_contaminant_removed.any())
        self.assertFalse(audit.literature_confirmed_dual.any())
        self.assertTrue(audit.literature_uncertain_classification.all())

    def test_published_geometry_requires_separation_agreement(self):
        _, audit = self.audit(self.reference(secondary_ra_deg=np.nan, secondary_dec_deg=np.nan,
            match_rule='published_pair_geometry', separation_arcsec=.9))
        self.assertFalse(audit.literature_contaminant_removed.any())

    def test_missing_reference_positions_remain_unassessed(self):
        associations, audit = self.audit(self.reference(ra_deg=np.nan, dec_deg=np.nan,
            secondary_ra_deg=np.nan, secondary_dec_deg=np.nan, match_rule='position_only'))
        self.assertTrue(associations.empty)
        self.assertFalse(audit.literature_contaminant_removed.any())

    def test_reviewed_system_geometry_accepts_a_published_system_position(self):
        _, audit = self.audit(self.reference(ra_deg=10.0001,
            secondary_ra_deg=np.nan, secondary_dec_deg=np.nan,
            match_rule='reviewed_system_geometry'))
        self.assertTrue(audit.literature_contaminant_removed.all())

    def test_bundled_reference_scope_and_qualified_lens_policy(self):
        references = read_references(LITERATURE_ROOT/'known_binaries.csv',
                                     LITERATURE_ROOT/'known_lenses_contaminants.csv')
        binaries = references.loc[references.reference_catalogue.eq('known_binaries')]
        self.assertEqual(len(binaries), 55)
        self.assertEqual(set(binaries.loc[binaries.ra_deg.isna(), 'system_name']),
                         {'J0400+2754', 'J2332+2243'})
        j0823 = references.loc[references.system_id.eq('J0823+2418')]
        self.assertFalse(j0823.empty)
        self.assertTrue(j0823.classification.eq('likely_lens').all())

    def test_literature_stage_precedes_final_candidate_plots(self):
        self.assertLess(VETTING_STAGES.index('combine_candidate_catalogues'),
                        VETTING_STAGES.index('crossmatch_literature'))
        self.assertLess(VETTING_STAGES.index('crossmatch_literature'),
                        VETTING_STAGES.index('plot_candidate_separation_redshift'))

    def test_dual_or_lensed_wording_is_not_a_lens_confirmation(self):
        references = read_references(LITERATURE_ROOT/'known_binaries.csv',
                                     LITERATURE_ROOT/'known_lenses_contaminants.csv')
        qualified = references.loc[references.reference_id.str.startswith('Chen2025:')
                                   & references.notes.str.contains('dual/lensed quasar', regex=False)]
        self.assertFalse(qualified.empty)
        self.assertFalse(qualified.classification.eq('confirmed_lens').any())


if __name__ == '__main__':
    unittest.main()
