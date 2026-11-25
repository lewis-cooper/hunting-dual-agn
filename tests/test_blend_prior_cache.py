"""Check that cached prior interpolation preserves the pinned model's values."""

from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dual_agn.blend_priors import CachedMagnitudePrior


class BlendPriorCacheTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            from qso_pcolor.priors import GridQSOPrior
        except ImportError as error:
            raise unittest.SkipTest('Install the pinned qso_p_color dependency to check prior interpolation') from error
        cls.prior = GridQSOPrior(
            z_centres=np.array([.5, 1.5, 2.5]),
            mag_centres=np.array([17.5, 18.5, 20.]),
            sigma=np.array([[1., 3., 7.], [2., 8., 4.], [9., 5., 6.]]),
            z_edges=np.array([0., 1., 2., 3.]),
            mag_edges=np.array([17., 18., 19., 21.]))

    def test_cached_scalar_grid_and_off_grid_values_equal_original(self):
        cached = CachedMagnitudePrior(self.prior, 'quasar', mag_step=.1)
        grid = np.arange(17., 21.+1e-9, .1)
        magnitudes = np.r_[grid, 16.99, 21.01, 18.123456789]
        redshifts = np.array([-.01, 0., .25, .5, 1.25, 2.75, 3., 3.01])
        for magnitude in magnitudes:
            expected = self.prior(redshifts, float(magnitude))
            actual = np.array([cached([z], float(magnitude))[0] for z in redshifts])
            np.testing.assert_array_equal(actual, expected)

    def test_vector_input_keeps_original_evaluation(self):
        cached = CachedMagnitudePrior(self.prior, 'quasar', mag_step=.1)
        redshifts = np.array([0., .75, 1.9, 3.])
        np.testing.assert_array_equal(cached(redshifts, 18.25), self.prior(redshifts, 18.25))


if __name__ == '__main__':
    unittest.main()
