from pathlib import Path
import json
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from settings import *
from photometry import read_fits
from imaging import attach_quaia_redshift_errors, prepare_component_photometry
from photometry import RESOLVED_POLICY
candidate_systems = pd.read_csv(
    IMAGING_ROOT / 'precolour_systems.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)
quaia_reference = read_fits(INPUT_ROOT / 'quaia_imaging_reference.fits')
candidate_systems = attach_quaia_redshift_errors(candidate_systems, quaia_reference)
candidate_systems.to_csv(IMAGING_ROOT / 'precolour_systems.csv', index=False)
component_rows = read_fits(INPUT_ROOT / 'imaging_components.fits')
tractor_sources = read_fits(INPUT_ROOT / 'imaging_tractor.fits')
component_photometry = prepare_component_photometry(candidate_systems, component_rows, tractor_sources)
component_photometry['reference_snr_pass'] = component_photometry.r_snr.ge(10)
component_photometry['latitude_pass'] = component_photometry.galactic_abs_b_deg.ge(25)
component_photometry['psf_pass'] = component_photometry.tractor_type.eq('PSF')
component_photometry['policy_pass'] = ~RESOLVED_POLICY.violations(
    component_photometry.pair_separation_arcsec.to_numpy(float),
    component_photometry.tractor_fracflux_r.to_numpy(float),
    len(component_photometry)
)
component_photometry['close_separation'] = np.isfinite(component_photometry.pair_separation_arcsec) & component_photometry.pair_separation_arcsec.ge(0) & component_photometry.pair_separation_arcsec.lt(CLOSE_SEPARATION_ARCSEC)
component_photometry.to_csv(IMAGING_ROOT / 'component_photometry.csv', index=False)
print('Components:', len(component_photometry), '| Output:', IMAGING_ROOT / 'component_photometry.csv')
