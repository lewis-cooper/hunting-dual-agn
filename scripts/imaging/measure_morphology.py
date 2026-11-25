"""Compute r-band morphology from user-supplied local FITS cutouts."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

from dual_agn.imaging_morphology import measure_image
from dual_agn.measurement_inputs import read_config, local_path
from dual_agn.imaging_search.catalogues import file_hash
from dual_agn.score_selection import (imaging_model_scores, read_either_model_threshold,
                             has_two_distinct_gaia_sources, build_secondary_pm_flags)
from dual_agn.settings import (ROOT, IMAGING_ROOT, IMAGING_MORPHOLOGY_PATH, IDENTIFIER_DTYPES,
                      THRESHOLD_CALIBRATION_PATH, PM_SIGNIFICANCE_MAX)


def image_path(config, row):
    specification = config['morphology_images']
    values = row.to_dict()
    values['qso_parent_index'] = int(row.qso_parent_index)
    return local_path(config, specification['path_template'].format(**values))


def measure_candidates(candidates, config):
    specification = config['morphology_images']
    unsupported = set(specification).difference({'path_template', 'r_band_index'})
    if unsupported:
        raise ValueError('Unsupported morphology_images settings: ' + ', '.join(sorted(unsupported)))
    rows, unavailable = [], []
    for _, system in candidates.iterrows():
        path = image_path(config, system)
        if not path.is_file():
            unavailable.append(str(path))
            continue
        with fits.open(path, memmap=False) as hdus:
            image = np.asarray(hdus[0].data, float)
            if image.ndim == 3:
                band_index = 0 if image.shape[0] == 1 else specification.get('r_band_index')
                if band_index is None:
                    if image.shape[0] != 1:
                        raise ValueError('Set r_band_index for a local multiband image cube')
                    band_index = 0
                image = image[int(band_index)]
            if image.ndim != 2:
                raise ValueError(f'Expected a 2D image or 3D cube: {path}')
            x,y=WCS(hdus[0].header).celestial.world_to_pixel_values(
                [system.primary_ra,system.secondary_ra],[system.primary_dec,system.secondary_dec])
            inside=(np.isfinite(x).all() and np.isfinite(y).all() and (x>=0).all()
                    and (y>=0).all() and (x<image.shape[1]).all() and (y<image.shape[0]).all())
            if not inside:
                raise ValueError(f'Supplied image does not contain both components: {path}')
            metrics = measure_image(image, hdus[0].header, system.primary_ra, system.primary_dec,
                                    system.secondary_ra, system.secondary_dec)
        rows.append(dict(qso_parent_index=int(system.qso_parent_index), **metrics,
                         morphology_measurement_source='computed_rband_profile_moments',
                         morphology_image_path=str(path), morphology_image_sha256=file_hash(path),
                         morphology_measurement_status='measured' if np.isfinite(metrics['r_perp_sigma_ratio']) else 'no_finite_profile_width'))
    if unavailable:
        raise FileNotFoundError(f'{len(unavailable)} selected fields lack r-band images. Supply the required local FITS images. No objects were silently excluded. Examples:\n' + '\n'.join(unavailable[:5]))
    columns = ['qso_parent_index', 'image_axis_sigma_arcsec', 'image_perp_sigma_arcsec',
               'r_perp_sigma_ratio', 'morphology_measurement_source']
    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=columns)


def morphology_candidates(systems, threshold):
    """Return the same score/PM-qualified fields used by the morphology stage."""
    two_gaia = has_two_distinct_gaia_sources(systems.primary_gaia_source_id, systems.secondary_gaia_source_id)
    pm = build_secondary_pm_flags(systems.secondary_pm_sig, two_gaia, PM_SIGNIFICANCE_MAX)
    requires_morphology = imaging_model_scores(systems).max(axis=1).ge(threshold) & pm.secondary_pm_pass
    return systems.loc[requires_morphology].copy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT / 'measurement_inputs.local.json')
    args = parser.parse_args()
    config = read_config(args.config)
    systems = pd.read_csv(IMAGING_ROOT / 'system_scores.csv', dtype=IDENTIFIER_DTYPES, float_precision='round_trip')
    threshold = read_either_model_threshold(THRESHOLD_CALIBRATION_PATH)
    measurements = measure_candidates(morphology_candidates(systems, threshold), config)
    IMAGING_MORPHOLOGY_PATH.parent.mkdir(parents=True, exist_ok=True)
    measurements.to_csv(IMAGING_MORPHOLOGY_PATH, index=False)
    print(f'Morphology from r-band images: {len(measurements)} systems; output={IMAGING_MORPHOLOGY_PATH}; finite ratios={np.isfinite(measurements.r_perp_sigma_ratio).sum()}; threshold={threshold:.2f}')


if __name__ == '__main__':
    main()
