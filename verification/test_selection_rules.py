"""Check DCR and PM boundaries and the source of adopted Quaia uncertainties."""

from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dual_agn.imaging import attach_quaia_redshift_errors
from dual_agn.score_selection import (build_imaging_dcr_flags, build_secondary_pm_flags,
                             has_two_distinct_gaia_sources)


class ImagingDCRSelectionTests(unittest.TestCase):
    def test_one_gaia_requires_finite_angle_strictly_above_twenty(self):
        angle = pd.Series([np.nan, 0.0, 20.0, 20.01, np.inf])
        flags = build_imaging_dcr_flags(angle, pd.Series(True, index=angle.index), 20.0)
        self.assertEqual(flags.dcr_pass.tolist(), [False, False, False, True, False])
        self.assertTrue(flags.dcr_cut_applied.all())

    def test_two_gaia_is_exempt_even_when_the_angle_is_missing(self):
        angle = pd.Series([np.nan, 0.0, 20.0, 20.01, np.inf])
        flags = build_imaging_dcr_flags(angle, pd.Series(False, index=angle.index), 20.0)
        self.assertTrue(flags.dcr_pass.all())
        self.assertFalse(flags.dcr_cut_applied.any())


class SecondaryPMSelectionTests(unittest.TestCase):
    def test_measured_threshold_and_missing_measurements(self):
        significance = pd.Series([2.99, 3.0, 10.0, np.nan, np.inf, -np.inf])
        flags = build_secondary_pm_flags(significance, pd.Series(True, index=significance.index), 3.0)
        self.assertEqual(flags.secondary_pm_pass.tolist(), [True, False, False, True, True, True])
        self.assertEqual(flags.secondary_pm_missing.tolist(), [False, False, False, True, True, True])

    def test_one_gaia_systems_have_no_secondary_pm_cut(self):
        significance = pd.Series([10.0, np.nan])
        flags = build_secondary_pm_flags(significance, pd.Series([False, False]), 3.0)
        self.assertTrue(flags.secondary_pm_pass.all())
        self.assertFalse(flags.secondary_pm_missing.any())

    def test_gaia_ids_are_exact_and_distinct(self):
        primary = pd.Series(['1000993221264706432'] * 4)
        secondary = pd.Series(['1000993221264706433', '1000993221264706432', '0', None])
        self.assertEqual(has_two_distinct_gaia_sources(primary, secondary).tolist(),
                         [True, False, False, False])

    def test_missing_and_placeholder_ids_work_with_both_string_backends(self):
        for storage in ['python', 'pyarrow']:
            with self.subTest(storage=storage):
                dtype = pd.StringDtype(storage=storage)
                primary = pd.Series(['1000993221264706432'] * 8, dtype=dtype)
                secondary = pd.Series(['1000993221264706433', '1000993221264706432',
                                       '0', None, '', 'nan', '1.0', '-2'], dtype=dtype)
                self.assertEqual(has_two_distinct_gaia_sources(primary, secondary).tolist(),
                                 [True, False, False, False, False, False, False, False])


class QuaiaUncertaintyJoinTests(unittest.TestCase):
    def setUp(self):
        self.candidates = pd.DataFrame({
            'qso_parent_index': [1, 2], 'source_id': ['1000993221264706432', '20'],
            'primary_gaia_source_id': ['1000993221264706433', '20'],
            'adopted_redshift_catalogue': ['quaia', 'desi'], 'adopted_redshift': [1.5, 2.0]})
        self.reference = pd.DataFrame({'source_id': ['1000993221264706432'],
                                       'redshift_quaia': [1.5], 'redshift_quaia_err': [0.1]})

    def test_use_the_record_supplying_redshift_after_relabelling(self):
        enriched = attach_quaia_redshift_errors(self.candidates, self.reference)
        self.assertEqual(enriched.adopted_redshift_error.tolist(), [0.1, 0.0])
        self.assertEqual(enriched.redshift_error_source_id.iloc[0], '1000993221264706432')

    def test_missing_quaia_reference_is_an_error(self):
        with self.assertRaises(ValueError):
            attach_quaia_redshift_errors(self.candidates, self.reference.iloc[0:0])

    def test_inconsistent_adopted_redshift_is_an_error(self):
        self.reference.loc[0, 'redshift_quaia'] = 1.6
        with self.assertRaises(ValueError):
            attach_quaia_redshift_errors(self.candidates, self.reference)


if __name__ == '__main__':
    unittest.main()
