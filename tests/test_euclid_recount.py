"""Distinguish footprint membership, actual image pixels and visual review."""

from pathlib import Path
import tempfile
import unittest

import healpy as hp
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from recount_euclid import in_moc, measure_vis_pixels, read_moc_intervals, summarize_audit
from test_publication_catalogue import sample_candidates


class EuclidRecountTests(unittest.TestCase):
    def test_nuniq_mixed_orders_cover_the_expected_sky_cells(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'coverage.fits'
            table = fits.BinTableHDU.from_columns([
                fits.Column(name='PIXEL', format='K', array=np.array([4, 36], dtype=np.int64))])
            table.header.update(ORDERING='NUNIQ', COORDSYS='C', MOCORDER=1)
            fits.HDUList([fits.PrimaryHDU(), table]).writeto(path)
            inside_coarse = hp.pix2ang(1, 0, nest=True, lonlat=True)
            inside_fine = hp.pix2ang(2, 20, nest=True, lonlat=True)
            outside = hp.pix2ang(1, 11, nest=True, lonlat=True)
            coordinates = np.array([inside_coarse, inside_fine, outside])
            self.assertEqual(in_moc(read_moc_intervals(path), coordinates[:, 0], coordinates[:, 1]).tolist(),
                             [True, True, False])

    def test_partial_vis_patch_is_usable_but_zero_image_is_not(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'image.npz'
            primary = SkyCoord(10 * u.deg, 20 * u.deg)
            secondary = primary.directional_offset_by(90 * u.deg, .5 * u.arcsec)
            record = dict(mid_ra=10., mid_dec=20., grid=dict(
                pixel_arcsec=.1, north_up=True, east_left=True))
            image = np.ones((60, 60))
            wcs = WCS(naxis=2)
            wcs.wcs.crpix = [30.5, 30.5]
            wcs.wcs.cdelt = [-.1 / 3600, .1 / 3600]
            wcs.wcs.crval = [10., 20.]
            wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
            x, y = wcs.world_to_pixel_values(secondary.ra.deg, secondary.dec.deg)
            image[int(round(float(y))), int(round(float(x)))] = np.nan
            np.savez(path, VIS=image)
            measured = measure_vis_pixels(record, path, primary, secondary)
            self.assertTrue(measured['usable_vis'])
            self.assertGreater(measured['secondary_patch_valid_pixels'], 0)
            np.savez(path, VIS=np.zeros((60, 60)))
            self.assertFalse(measure_vis_pixels(record, path, primary, secondary)['usable_vis'])

    def test_routes_overlap_but_tiers_count_unique_systems(self):
        candidates = sample_candidates()
        audit = pd.DataFrame([
            dict(candidate_id='system:0', covered_both=True, usable_vis=True, vis_class='resolved_pair'),
            dict(candidate_id='system:1', covered_both=True, usable_vis=True, vis_class='potential_blend'),
            dict(candidate_id='system:2', covered_both=True, usable_vis=False, vis_class='no_usable_vis'),
        ])
        counts = summarize_audit(audit, candidates)
        self.assertEqual(counts['NeuclidVIS'], 2)
        self.assertEqual(counts['NeuclidShared'], 1)
        self.assertEqual(counts['NeuclidMain'] + counts['NeuclidAdditional'], counts['NeuclidVIS'])
        self.assertEqual(counts['NeuclidResolved'] + counts['NeuclidBlended'] + counts['NeuclidImaging']
                         - counts['NeuclidShared'], counts['NeuclidVIS'])
        audit.loc[1, 'vis_class'] = 'unreviewed'
        with self.assertRaisesRegex(ValueError, 'visual classifications'):
            summarize_audit(audit, candidates)


if __name__ == '__main__':
    unittest.main()
