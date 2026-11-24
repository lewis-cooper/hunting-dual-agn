"""Offline DCR axes from local, quality-selected g-band CCD metadata."""

from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.time import Time
from astropy.utils import iers
import numpy as np
import pandas as pd


OBSERVATORIES = {'decam':(-30.169,-70.806),
                 '90prime':(31.9583,-111.5986),
                 'mosaic':(31.9583,-111.5986)}
iers.conf.auto_download = False


def axis_delta(first_deg, second_deg):
    return np.abs((np.asarray(first_deg)-np.asarray(second_deg)+90.) % 180. - 90.)


def pair_pa(primary_ra,primary_dec,secondary_ra,secondary_dec):
    east = (np.asarray(secondary_ra)-np.asarray(primary_ra))*np.cos(np.deg2rad(primary_dec))
    north = np.asarray(secondary_dec)-np.asarray(primary_dec)
    return np.rad2deg(np.arctan2(east,north)) % 180.


def footprint_contains(ccds, ra_deg, dec_deg):
    """Vectorized inverse of the published TAN WCS, using zero-based pixels."""
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    ra0 = np.deg2rad(ccds.crval1.to_numpy())
    dec0 = np.deg2rad(ccds.crval2.to_numpy())
    delta_ra = ra-ra0
    denominator = np.sin(dec)*np.sin(dec0)+np.cos(dec)*np.cos(dec0)*np.cos(delta_ra)
    with np.errstate(divide='ignore',invalid='ignore'):
        east_deg = np.rad2deg(np.cos(dec)*np.sin(delta_ra)/denominator)
        north_deg = np.rad2deg((np.sin(dec)*np.cos(dec0)-np.cos(dec)*np.sin(dec0)*np.cos(delta_ra))/denominator)
        determinant = ccds.cd1_1.to_numpy()*ccds.cd2_2.to_numpy()-ccds.cd1_2.to_numpy()*ccds.cd2_1.to_numpy()
        x = (ccds.cd2_2.to_numpy()*east_deg-ccds.cd1_2.to_numpy()*north_deg)/determinant+ccds.crpix1.to_numpy()-1.
        y = (-ccds.cd2_1.to_numpy()*east_deg+ccds.cd1_1.to_numpy()*north_deg)/determinant+ccds.crpix2.to_numpy()-1.
    return ((denominator > 0) & np.isfinite(x) & np.isfinite(y)
            & (x >= 0) & (x < ccds.width.to_numpy()) & (y >= 0) & (y < ccds.height.to_numpy()))


def parallactic_pa(ra_deg, dec_deg, mjd, camera):
    latitude, longitude = OBSERVATORIES[camera]
    sidereal_deg = Time(np.asarray(mjd),format='mjd',scale='utc').sidereal_time('apparent',longitude=longitude*u.deg).deg
    hour_angle = np.deg2rad((sidereal_deg-np.asarray(ra_deg)+180.) % 360.-180.)
    dec = np.deg2rad(dec_deg)
    return np.rad2deg(np.arctan2(np.sin(hour_angle),
        np.tan(np.deg2rad(latitude))*np.cos(dec)-np.sin(dec)*np.cos(hour_angle))) % 180.


def calculate_dcr(summary, ccds, coverage):
    targets = summary.to_pandas().set_index('qso_parent_index',drop=False)
    required_coverage = {'qso_parent_index','ra_deg','dec_deg','query_radius_deg'}
    if not required_coverage.issubset(coverage):
        raise ValueError('CCD coverage needs qso_parent_index, ra_deg, dec_deg, query_radius_deg, including successful empty queries.')
    if not coverage.qso_parent_index.is_unique:
        raise ValueError('CCD coverage must have one row per parent')
    covered = coverage.set_index('qso_parent_index').reindex(targets.index)
    valid = np.isfinite(covered[['ra_deg','dec_deg']].to_numpy(float)).all(axis=1)
    offsets = np.full(len(targets),np.inf)
    if valid.any():
        offsets[valid] = SkyCoord(covered.ra_deg.to_numpy()[valid]*u.deg,covered.dec_deg.to_numpy()[valid]*u.deg).separation(
            SkyCoord(targets.primary_ra.to_numpy()[valid]*u.deg,targets.primary_dec.to_numpy()[valid]*u.deg)).arcsec
    covered_pass = (offsets <= 1e-6) & covered.query_radius_deg.ge(.65).to_numpy()
    if not covered_pass.all():
        raise ValueError(f'{(~covered_pass).sum()} pairs lack complete, coordinate-consistent 0.65-degree CCD query coverage. No download was attempted.')
    if 'qso_parent_index' not in ccds:
        if not {'ra','dec'}.issubset(ccds):
            raise ValueError('Unassociated CCD catalogues need centre ra and dec in degrees')
        target_sky = SkyCoord(targets.primary_ra.to_numpy()*u.deg,targets.primary_dec.to_numpy()*u.deg)
        ccd_sky = SkyCoord(ccds.ra.to_numpy()*u.deg,ccds.dec.to_numpy()*u.deg)
        target_indices, ccd_indices, _, _ = search_around_sky(target_sky,ccd_sky,.65*u.deg)
        ccds = ccds.iloc[ccd_indices].copy()
        ccds['qso_parent_index'] = targets.index.to_numpy()[target_indices]
    ccds = ccds.loc[ccds.qso_parent_index.isin(targets.index)].copy()
    if {'target_ra','target_dec'}.issubset(ccds):
        expected = targets.loc[ccds.qso_parent_index]
        offsets = SkyCoord(ccds.target_ra.to_numpy()*u.deg,ccds.target_dec.to_numpy()*u.deg).separation(
            SkyCoord(expected.primary_ra.to_numpy()*u.deg,expected.primary_dec.to_numpy()*u.deg)).arcsec
        if (offsets > 1e-6).any():
            raise ValueError('Cached CCD target coordinates do not identify the discovered r-ranked primary')
    quality = ccds['filter'].astype(str).str.strip().eq('g') & ccds.ccd_cuts.eq(0)
    ccds['camera'] = ccds.camera.astype(str).str.strip().str.lower()
    quality &= ccds.camera.isin(OBSERVATORIES)
    ccds = ccds.loc[quality].copy()
    row_targets = targets.loc[ccds.qso_parent_index]
    inside = footprint_contains(ccds,row_targets.primary_ra.to_numpy(),row_targets.primary_dec.to_numpy())
    ccds = ccds.loc[inside].reset_index(drop=True)
    row_targets = targets.loc[ccds.qso_parent_index]
    dcr_pa = np.full(len(ccds),np.nan)
    for camera in OBSERVATORIES:
        selected = ccds.camera.eq(camera).to_numpy()
        if selected.any():
            dcr_pa[selected] = parallactic_pa(row_targets.primary_ra.to_numpy()[selected],
                row_targets.primary_dec.to_numpy()[selected],ccds.mjd_obs.to_numpy()[selected],camera)
    exptime = ccds.exptime.to_numpy(float)
    exptime = np.where(np.isfinite(exptime),exptime,1.)
    sig1 = ccds.sig1.to_numpy(float)
    measured_noise = np.isfinite(sig1) & (sig1 > 0)
    weights = exptime.copy()
    weights[measured_noise] /= sig1[measured_noise]**2
    valid = np.isfinite(dcr_pa) & np.isfinite(weights) & (weights > 0)
    moments = pd.DataFrame({'qso_parent_index':ccds.qso_parent_index,
        'cosine':np.where(valid,weights*np.cos(2*np.deg2rad(dcr_pa)),0.),
        'sine':np.where(valid,weights*np.sin(2*np.deg2rad(dcr_pa)),0.),
        'weight':np.where(valid,weights,0.)})
    moments = moments.groupby('qso_parent_index')[['cosine','sine','weight']].sum()
    mean = .5*np.rad2deg(np.arctan2(moments.sine,moments.cosine)) % 180.
    mean = mean.where(moments.weight.gt(0))
    result = pd.DataFrame({'qso_parent_index':targets.index})
    result['dcr_pa_deg'] = result.qso_parent_index.map(mean)
    pa = pair_pa(targets.primary_ra,targets.primary_dec,targets.secondary_ra,targets.secondary_dec)
    result['pair_dcr_delta_deg'] = axis_delta(pa,result.dcr_pa_deg.to_numpy())
    result['n_g_ccd'] = result.qso_parent_index.map(ccds.groupby('qso_parent_index').size()).fillna(0).astype(int)
    return result
