from pathlib import Path
import json
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18
from settings import *
from photometry import configure_dustmaps, saved_true
from qso_pcolor import Photometry
from qso_pcolor.extinction import deredden
from qso_pcolor.unified import UnifiedPSFModel
from qso_pcolor.blend import BlendModel
import colour_projection as colour_tools
from plot_style import plt, mpl, save_figure
configure_dustmaps()
log_r_per_unit_z_threshold = json.loads(THRESHOLD_CALIBRATION_PATH.read_text())['adopted_log_r_threshold']

scoring_inputs = pd.read_csv(
    ONE_ROOT / 'scoring_inputs.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)
candidate_ids = scoring_inputs.qso_group_id.to_numpy()
measurements = []

def add_band(ids, survey, source, band, snr, usable):
    usable = np.asarray(usable, dtype=bool)
    measurements.append(pd.DataFrame({
        'qso_group_id': np.asarray(ids)[usable],
        'survey': survey,
        'source': source,
        'band': band,
        'snr': np.asarray(snr, dtype=float)[usable]
    }))
# Source cell ec0b9b51
for source in ('primary', 'secondary'):
    release = scoring_inputs[f'{source}_decals_release']
    for band in 'grz':
        flux = pd.to_numeric(scoring_inputs[f'{source}_decals_flux_{band}'], errors='coerce')
        ivar = pd.to_numeric(scoring_inputs[f'{source}_decals_flux_ivar_{band}'], errors='coerce')
        usable = release.isin((9010, 9011, 9012)) & np.isfinite(flux) & np.isfinite(ivar) & ivar.gt(0)
        add_band(candidate_ids, 'DECaLS', source, band, flux * np.sqrt(ivar.clip(lower=0)), usable)


# Source cell 61e6096d
for survey, spec in survey_products.items():
    phot = pd.read_csv(spec['path'], dtype={'qso_group_id': str}, low_memory=False)
    phot = phot.set_index('qso_group_id').reindex(candidate_ids)
    source_ok = saved_true(phot[spec['source_flag']])
    if survey == 'sdss':
        source_ok &= saved_true(phot.sdss_clean)
        source_ok &= pd.to_numeric(phot.sdss_mode, errors='coerce').eq(1)
    for band, (mag_column, error_column) in spec['bands'].items():
        mag = pd.to_numeric(phot[mag_column], errors='coerce')
        mag_error = pd.to_numeric(phot[error_column], errors='coerce')
        usable = source_ok & mag.between(-50, 90, inclusive='neither') & mag_error.gt(0) & mag_error.lt(90)
        if survey in ('ps1', 'skymapper', 'nsc', 'vhs'):
            usable &= saved_true(phot[f'{survey}_{band}_band_quality'])
        snr = 2.5 / np.log(10) / mag_error
        add_band(candidate_ids, survey.upper(), 'survey_match', band, snr, usable)


# Source cell 7850b63a
band_snr_df = pd.concat(measurements, ignore_index=True)
band_snr_df = scoring_inputs[['qso_group_id']].merge(band_snr_df, on='qso_group_id', how='left', validate='one_to_many').sort_values(['qso_group_id', 'survey', 'source', 'band'], na_position='last')
print(f'{band_snr_df.qso_group_id.nunique()} candidates; {band_snr_df.snr.notna().sum()} usable band measurements')

band_snr_df.to_csv(ONE_ROOT / 'photometry_band_inventory.csv', index=False)
print('Output:', ONE_ROOT / 'photometry_band_inventory.csv')
