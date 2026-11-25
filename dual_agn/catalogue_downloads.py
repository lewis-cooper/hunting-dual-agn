"""Public raw inputs used where the original workflow did not query WSDB."""

import hashlib
from io import StringIO, BytesIO
import json
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

from dual_agn.imaging_search.catalogues import file_hash
from dual_agn.measurement_inputs import sky


CCD_FIELDS = ['ra','dec','filter','ccd_cuts','camera','width','height','crpix1','crpix2',
              'crval1','crval2','cd1_1','cd1_2','cd2_1','cd2_2','mjd_obs','exptime','sig1']


def download_legacy_r_cutout(system, path, size_pixels=100, pixel_scale_arcsec=.1, session=None):
    """Download a DR9 r image atomically; reject invalid or wrong-field FITS data."""
    if int(size_pixels) != size_pixels or size_pixels < 21 or not np.isfinite(pixel_scale_arcsec) or pixel_scale_arcsec <= 0:
        raise ValueError('Cutouts require integer size_pixels >= 21 and positive pixel_scale_arcsec')
    primary = SkyCoord(system.primary_ra*u.deg, system.primary_dec*u.deg)
    secondary = SkyCoord(system.secondary_ra*u.deg, system.secondary_dec*u.deg)
    centre = primary.directional_offset_by(primary.position_angle(secondary), primary.separation(secondary)/2)
    for ra_field, dec_field in [('cutout_ra','cutout_dec'),('gaia_ra','gaia_dec')]:
        coordinates = pd.to_numeric([system.get(ra_field,np.nan), system.get(dec_field,np.nan)],errors='coerce')
        if np.isfinite(coordinates).all():
            centre = SkyCoord(coordinates[0]*u.deg,coordinates[1]*u.deg)
            break
    url = 'https://www.legacysurvey.org/viewer/fits-cutout'
    # Preserve the centre/coordinate rounding and image sampling of the original
    # download_decals_fits_cutouts_for_targets.py acquisition script.
    parameters = dict(ra=f'{centre.ra.deg:.7f}', dec=f'{centre.dec.deg:.7f}', layer='ls-dr9', bands='r',
                      pixscale=f'{pixel_scale_arcsec:.6f}', size=str(int(size_pixels)))
    response = (session or requests).get(url, params=parameters, timeout=(20,180))
    response.raise_for_status()
    with fits.open(BytesIO(response.content), memmap=False) as hdus:
        image = np.asarray(hdus[0].data)
        if image.ndim == 3 and image.shape[0] == 1:
            image = image[0]
        if image.ndim != 2 or not np.isfinite(image).any():
            raise ValueError('Legacy cutout response is not a usable r-band image')
        wcs = WCS(hdus[0].header).celestial
        x, y = wcs.world_to_pixel_values([primary.ra.deg, secondary.ra.deg], [primary.dec.deg, secondary.dec.deg])
        if not (np.isfinite(x).all() and np.isfinite(y).all() and (x>=0).all()
                and (y>=0).all() and (x<image.shape[1]).all() and (y<image.shape[0]).all()):
            raise ValueError('Downloaded r-band image does not contain both component coordinates')
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.partial')
    try:
        temporary.write_bytes(response.content)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    path.with_suffix('.download.json').write_text(json.dumps({'url':url,'parameters':parameters,
        'sha256':file_hash(path)},indent=2)+'\n')
    print('Downloaded DR9 r-band cutout:', path, flush=True)
    return path


def public_query(service, targets, radius_arcsec):
    if service == 'legacy_dr9_ccds':
        values = ','.join(f'({int(row.qso_parent_index)},{float(row.ra):.15g},{float(row.dec):.15g})'
                          for row in targets.itertuples())
        fields = ','.join('c.'+field for field in CCD_FIELDS)
        query = f'''WITH targets(qso_parent_index,ra,dec) AS (VALUES {values})
            SELECT t.qso_parent_index, t.ra AS target_ra, t.dec AS target_dec, {fields}
            FROM targets AS t JOIN ls_dr9.ccds_annotated AS c
            ON q3c_radial_query(c.ra,c.dec,t.ra,t.dec,{radius_arcsec/3600.:.15g})
            WHERE c.filter='g' AND c.ccd_cuts=0'''
        return ('https://datalab.noirlab.edu/query/query',
                {'sql': query, 'ofmt': 'csv', 'async': 'False'},
                ['qso_parent_index','target_ra','target_dec',*CCD_FIELDS])
    raise ValueError('The supported public catalogue input is legacy_dr9_ccds')


def export_public_sources(service, targets, path, radius_arcsec, resume=False, session=None):
    """Save raw sources and successful cones; reuse only hash-verified responses."""
    sky(targets)
    if not np.isfinite(radius_arcsec) or radius_arcsec <= 0:
        raise ValueError('Public query radius must be finite and positive')
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    coverage_path = path.with_name(path.stem+'_coverage.csv')
    summary_path = path.with_suffix('.query.json')
    fingerprint = {'service': service, 'radius_arcsec': radius_arcsec,
        'targets_sha256': hashlib.sha256(targets.to_csv(index=False).encode()).hexdigest(),
        'exporter_sha256': file_hash(Path(__file__))}
    if resume and path.is_file() and coverage_path.is_file() and summary_path.is_file():
        previous = json.loads(summary_path.read_text())
        if (previous.get('fingerprint') == fingerprint and previous.get('sha256') == file_hash(path)
                and previous.get('coverage_sha256') == file_hash(coverage_path)):
            return path, coverage_path
    temporary = path.with_suffix('.partial')
    chunk_root = path.parent/(path.stem+'_query_chunks')
    chunk_root.mkdir(exist_ok=True)
    session = session or requests.Session()
    total_rows = 0
    wrote_header = False
    try:
        for start in range(0, len(targets), 30):
            chunk = targets.iloc[start:start+30]
            url, parameters, required = public_query(service, chunk, radius_arcsec)
            query_hash = hashlib.sha256(json.dumps([url,parameters],sort_keys=True).encode()).hexdigest()
            cache = chunk_root/(query_hash+'.csv')
            cache_hash = cache.with_suffix('.sha256')
            if resume and cache.is_file() and cache_hash.is_file() and file_hash(cache) == cache_hash.read_text().strip():
                payload = cache.read_text()
            else:
                # The service's documented anonymous token is public, as in
                # the existing control-catalogue downloader; no user login.
                headers = {'X-DL-AuthToken':'anonymous.0.0.anon_access',
                           'X-DL-TimeoutRequest':'120'} if service == 'legacy_dr9_ccds' else {}
                response = session.get(url, params=parameters, headers=headers, timeout=(20,180))
                response.raise_for_status()
                payload = response.text
            frame = pd.read_csv(StringIO(payload), comment='#', dtype={'source_id':'string'})
            if not set(required).issubset(frame):
                raise ValueError(f'{service} response lacks required native fields; response was not recorded as coverage')
            cache.write_text(payload)
            cache_hash.write_text(file_hash(cache)+'\n')
            frame.to_csv(temporary, index=False, mode='a' if wrote_header else 'w', header=not wrote_header)
            wrote_header = True
            total_rows += len(frame)
            print(f'{service}: {min(start+30,len(targets)):,}/{len(targets):,} query fields', flush=True)
        if not wrote_header:
            columns = ['qso_parent_index','target_ra','target_dec',*CCD_FIELDS]
            pd.DataFrame(columns=columns).to_csv(temporary,index=False)
        temporary.replace(path)
        coverage = targets.drop(columns=['ra','dec']).copy()
        coverage['ra_deg'], coverage['dec_deg'] = targets.ra.to_numpy(), targets.dec.to_numpy()
        coverage['query_radius_arcsec'] = radius_arcsec
        coverage.to_csv(coverage_path,index=False)
        summary_path.write_text(json.dumps({'fingerprint':fingerprint, 'source_rows':total_rows,
            'sha256':file_hash(path),'coverage_sha256':file_hash(coverage_path)},indent=2)+'\n')
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return path, coverage_path
