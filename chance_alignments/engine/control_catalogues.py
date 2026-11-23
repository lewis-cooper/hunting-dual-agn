"""Public catalogue queries and deliberately explicit control-sample preparation."""

from __future__ import annotations

import hashlib
import io
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

LATITUDE_EDGES = np.array([25., 40., 55., 90.000001])
ROUTES = ('resolved', 'blended', 'imaging_two_gaia', 'imaging_one_gaia')
QUERY_URL = 'https://datalab.noirlab.edu/query/query'
# This is the service's public anonymous token, not a user's credential.
PUBLIC_HEADERS = {'X-DL-AuthToken': 'anonymous.0.0.anon_access',
                  'X-DL-TimeoutRequest': '55'}
TRACTOR_COLUMNS = ('release,brickid,objid,ra,dec,type,flux_g,flux_r,flux_z,'
                   'flux_ivar_g,flux_ivar_r,flux_ivar_z,nobs_g,nobs_r,nobs_z,'
                   'maskbits,fracflux_r')
GAIA_COLUMNS = ('source_id,ra,dec,phot_g_mean_mag,pmra,pmra_error,pmdec,'
                'pmdec_error,pmra_pmdec_corr,parallax,parallax_error')


def read_csv(path):
    return pd.read_csv(path, dtype={'source_id': 'string', 'qso_group_id': 'string',
                                   'primary_gaia_source_id': 'string',
                                   'gaia_source_id': 'string', 'primary_source_id':'string',
                                   'secondary_gaia_source_id': 'string',
                                   'primary_component_gaia_source_id':'string',
                                   'secondary_component_gaia_source_id':'string'},
                       float_precision='round_trip', low_memory=False)


def read_fits(path, columns=None):
    table = Table.read(path, memmap=True)
    frame = (table if columns is None else table[columns]).to_pandas()
    for column in frame.select_dtypes('object'):
        frame[column] = frame[column].map(
            lambda value: value.decode().strip() if isinstance(value, bytes) else value)
    return frame


def latitude_bin(ra, dec):
    b = np.abs(SkyCoord(np.asarray(ra)*u.deg, np.asarray(dec)*u.deg).galactic.b.deg)
    return b, pd.cut(b, LATITUDE_EDGES, right=False, labels=['25-40', '40-55', '55-90'])


def magnitude(flux):
    flux = np.asarray(flux, float)
    return 22.5 - 2.5*np.log10(np.where(flux > 0, flux, np.nan))


def validate_cone_coverage(targets,outdir,prefix=''):
    coverage=read_csv(outdir/(prefix+'cone_coverage.csv')).set_index('target_id')
    if not targets.target_id.is_unique or not coverage.index.is_unique:
        raise ValueError('Control target identifiers must be unique')
    coverage=coverage.reindex(targets.target_id)
    valid=np.isclose(coverage.ra.to_numpy(),targets.ra.to_numpy(),rtol=0,atol=1e-10)
    valid &= np.isclose(coverage.dec.to_numpy(),targets.dec.to_numpy(),rtol=0,atol=1e-10)
    valid &= coverage.tractor_query_complete.eq(True).to_numpy() & coverage.gaia_query_complete.eq(True).to_numpy()
    if not valid.all():raise ValueError('Cones are missing or use different target coordinates; rerun fetch')


def pm_significance(frame):
    """The production covariance-aware sqrt(mu.T C^-1 mu) definition."""
    a = frame.pmra.to_numpy(float) / frame.pmra_error.to_numpy(float)
    b = frame.pmdec.to_numpy(float) / frame.pmdec_error.to_numpy(float)
    rho = frame.pmra_pmdec_corr.to_numpy(float)
    valid = np.isfinite(a+b+rho) & (np.abs(rho) < 1)
    valid &= frame.pmra_error.gt(0).to_numpy() & frame.pmdec_error.gt(0).to_numpy()
    value = np.full(len(frame), np.nan)
    value[valid] = np.sqrt(np.maximum(0, (a[valid]**2+b[valid]**2
                           -2*rho[valid]*a[valid]*b[valid]) / (1-rho[valid]**2)))
    return value






def public_query(sql, cache):
    """Cache successful CSV responses by query hash; errors never become empty cones."""
    cache.mkdir(parents=True,exist_ok=True)
    stem = hashlib.sha256(sql.encode()).hexdigest()
    path = cache/(stem+'.csv')
    if path.exists():
        return read_csv(path)
    (cache/(stem+'.sql')).write_text(sql+'\n')
    for attempt in range(5):
        try:
            response = requests.get(QUERY_URL,params={'sql':sql,'ofmt':'csv','async':'False'},
                                    headers=PUBLIC_HEADERS,timeout=(15,60))
            if response.status_code == 400:
                if 'too many clients' in response.text or 'timeout' in response.text.lower():
                    raise RuntimeError(response.text[:500])
                raise ValueError('Catalogue query rejected: '+response.text[:800])
            response.raise_for_status()
            if response.text.lstrip().startswith(('Error','<')):
                raise RuntimeError(response.text[:500])
            frame = pd.read_csv(io.StringIO(response.text),dtype={'source_id':'string'},float_precision='round_trip')
            path.write_text(response.text)
            return frame
        except (requests.RequestException,RuntimeError):
            if attempt == 4: raise
            time.sleep(3*(attempt+1))


def cone_sql(targets, survey, radius=60., point_sources_only=False, gaia_quality=False):
    columns, table = (TRACTOR_COLUMNS,'ls_dr9.tractor') if survey=='tractor' else (GAIA_COLUMNS,'gaia_dr3.gaia_source')
    if survey=='gaia' and gaia_quality:
        columns+=',ruwe,astrometric_params_solved'
    selected = ','.join('s.'+c for c in columns.split(','))
    primary = ' AND s.brick_primary=1' if survey=='tractor' else ''
    # Constant-centre cones keep the remote Q3C spatial index usable.
    queries=[]
    for row in targets.itertuples():
        point_filter=(f" AND (s.type='PSF' OR q3c_radial_query(s.ra,s.dec,{row.ra:.10f},"
            f'{row.dec:.10f},{3/3600:.12f}))') if survey=='tractor' and point_sources_only else ''
        queries.append(f"SELECT '{row.target_id}' AS target_id,{selected},"
            f'q3c_dist({row.ra:.10f},{row.dec:.10f},s.ra,s.dec)*3600 AS separation_arcsec '
            f'FROM {table} s WHERE q3c_radial_query(s.ra,s.dec,{row.ra:.10f},'
            f'{row.dec:.10f},{radius/3600:.12f}){primary}{point_filter}')
    return ' UNION ALL '.join(queries)




