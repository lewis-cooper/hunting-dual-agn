"""Local catalogue input, coverage validation and initial Tractor pair finding."""

from pathlib import Path
import glob
import hashlib
import json

from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.io import fits
from astropy.table import Table
import numpy as np
import pandas as pd

from .pairs import channel_row, flux_colour


MASTER_COLUMNS = ['qso_parent_index', 'source_id', 'gaia_ra', 'gaia_dec',
                  'adopted_redshift', 'adopted_redshift_catalogue', 'sdss_name',
                  'z', 'q_ra', 'q_dec', 'catalogue_memberships',
                  'in_sdss', 'in_desi', 'in_quaia', 'phot_g_mean_mag']
TRACTOR_COLUMNS = ['release', 'brickid', 'objid', 'brick_primary', 'type', 'ra', 'dec',
                   'flux_g', 'flux_r', 'flux_z', 'flux_ivar_g', 'flux_ivar_r', 'flux_ivar_z',
                   'nobs_g', 'nobs_r', 'nobs_z', 'maskbits', 'fracflux_r',
                   'parent_gaia_source_id', 'qso_parent_index']
PATH_KEYS = ['master_catalogue', 'tractor_coverage', 'gaia_catalogue', 'gaia_coverage',
             'ccd_catalogue', 'ccd_coverage', 'quaia_catalogue']


def resolve_config(path):
    supplied = json.loads(path.read_text())
    config = dict(aperture_arcsec=1.5, tractor_radius_arcsec=10.,
                  gaia_radius_arcsec=.7, colour_redshift_half_width=.1)
    config.update(supplied)
    required = ['master_catalogue', 'tractor_catalogues', 'tractor_coverage']
    if any(not config.get(name) for name in required):
        raise ValueError('Set master_catalogue, tractor_catalogues and tractor_coverage in the imaging-search config.')
    for name in PATH_KEYS:
        if config.get(name):
            candidate = Path(config[name]).expanduser()
            config[name] = (path.resolve().parent / candidate).resolve() if not candidate.is_absolute() else candidate.resolve()
    patterns = config['tractor_catalogues']
    if not isinstance(patterns, list) or not patterns:
        raise ValueError('tractor_catalogues must be a nonempty list of local files or glob patterns.')
    files = []
    for pattern in patterns:
        candidate = Path(pattern).expanduser()
        candidate = path.resolve().parent / candidate if not candidate.is_absolute() else candidate
        matched = sorted(glob.glob(str(candidate)))
        if not matched:
            raise FileNotFoundError(f'No local Tractor files match {candidate}; supply survey data separately. No download was attempted.')
        files.extend(Path(name).resolve() for name in matched)
    config['tractor_catalogues'] = list(dict.fromkeys(files))
    for name in ['aperture_arcsec', 'tractor_radius_arcsec', 'gaia_radius_arcsec', 'colour_redshift_half_width']:
        if not np.isfinite(config[name]) or config[name] <= 0:
            raise ValueError(f'{name} must be finite and positive')
    if config['tractor_radius_arcsec'] < config['aperture_arcsec']:
        raise ValueError('Tractor coverage radius must include the pair-search aperture.')
    return config


def input_paths(config, discovery_only=False):
    names = ['master_catalogue', 'tractor_coverage']
    if not discovery_only:
        names += ['gaia_catalogue', 'gaia_coverage', 'ccd_catalogue', 'ccd_coverage', 'quaia_catalogue']
    missing_keys = [name for name in names if not config.get(name)]
    if missing_keys:
        raise ValueError('Missing local input settings: '+', '.join(missing_keys)+'. Use --stage discover to find pairs first without Gaia/CCD inputs.')
    paths = [config[name] for name in names] + config['tractor_catalogues']
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError('Missing local imaging-search data. No download was attempted:\n'+'\n'.join(missing))
    return list(dict.fromkeys(paths))


def file_hash(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read_chunks(path, columns=None, chunksize=100_000, float_precision='round_trip'):
    """Read only local CSV/FITS data, without opening survey services."""
    if path.suffix.lower() == '.csv':
        names = pd.read_csv(path, nrows=0).columns
        selected = [name for name in columns if name in names] if columns else None
        identifier_types = {name: 'string' for name in names if name.endswith('source_id')}
        yield from pd.read_csv(path, usecols=selected, dtype=identifier_types,
                               chunksize=chunksize, float_precision=float_precision)
    elif path.suffix.lower() in {'.fits', '.fit', '.fz'}:
        with fits.open(path, memmap=True) as hdus:
            data = hdus[1].data
            selected = [name for name in columns if name in data.names] if columns else data.names
            for start in range(0, len(data), chunksize):
                values = {}
                for name in selected:
                    column = np.asarray(data[name][start:start+chunksize])
                    if column.dtype.kind in 'iufc' and not column.dtype.isnative:
                        column = column.astype(column.dtype.newbyteorder('='))
                    if column.dtype.kind == 'S':
                        column = np.char.decode(column).astype(str)
                    values[name] = column
                frame = pd.DataFrame(values)
                for name in frame.columns:
                    if name.endswith('source_id'):
                        frame[name] = frame[name].astype('string')
                yield frame
    else:
        raise ValueError(f'Expected a local CSV or FITS catalogue: {path}')


def read_frame(path, columns=None, float_precision='round_trip'):
    frames = list(read_chunks(path, columns,float_precision=float_precision))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=columns)


def decoded_frame(table):
    frame = table.to_pandas()
    for name in frame.select_dtypes(include=['object','string']).columns:
        frame[name] = frame[name].map(lambda value: value.decode().strip() if isinstance(value,bytes) else value)
    return frame


def write_fits(frame, path):
    table = Table.from_pandas(frame)
    for name in table.colnames:
        if table[name].dtype.kind == 'O':
            values = frame[name].fillna('').astype(str).to_numpy(dtype=str)
            table[name] = values
    table.write(path, overwrite=True)


def read_master(path):
    # The original search used pandas' standard CSV parser. Its last-bit z values
    # must be retained because the scorer's support sampler hashes exact inputs.
    master = read_frame(path, MASTER_COLUMNS,float_precision='high')
    required = ['qso_parent_index', 'source_id', 'gaia_ra', 'gaia_dec',
                'adopted_redshift', 'adopted_redshift_catalogue']
    if not set(required).issubset(master.columns):
        raise ValueError('Master catalogue needs columns: '+', '.join(required))
    for name in ['qso_parent_index', 'source_id']:
        if master[name].isna().any() or not master[name].is_unique:
            raise ValueError(f'Master {name} must be non-missing and unique')
    parent_index = pd.to_numeric(master.qso_parent_index, errors='raise')
    if not (np.isfinite(parent_index) & parent_index.ge(0) & parent_index.eq(np.floor(parent_index))).all():
        raise ValueError('Master qso_parent_index must contain nonnegative integers')
    master['qso_parent_index'] = parent_index.astype('int64')
    if not master.source_id.str.fullmatch('[1-9][0-9]*').all():
        raise ValueError('Master source_id values must be exact positive Gaia identifiers')
    if not np.isfinite(master[['gaia_ra', 'gaia_dec']].to_numpy(float)).all():
        raise ValueError('Master Gaia coordinates must be finite degrees')
    master['z'] = master.adopted_redshift
    if 'sdss_name' not in master:
        master['sdss_name'] = master.source_id
    master['sdss_name'] = master.sdss_name.fillna('')
    for coordinate in ['ra', 'dec']:
        if 'q_'+coordinate not in master:
            master['q_'+coordinate] = master['gaia_'+coordinate]
    return master


def require_tractor_coverage(master, coverage, radius_arcsec):
    required = {'source_id', 'ra_deg', 'dec_deg', 'query_radius_arcsec'}
    if not required.issubset(coverage):
        raise ValueError('Tractor coverage needs source_id, ra_deg, dec_deg and query_radius_arcsec, including successful empty cones.')
    if not coverage.source_id.is_unique:
        raise ValueError('Tractor coverage source_id must be unique')
    covered = coverage.set_index('source_id').reindex(master.source_id)
    available = covered.query_radius_arcsec.ge(radius_arcsec).to_numpy()
    coordinates = np.isfinite(covered[['ra_deg','dec_deg']].to_numpy(float)).all(axis=1)
    offsets = np.full(len(master), np.inf)
    if coordinates.any():
        offsets[coordinates] = SkyCoord(covered.ra_deg.to_numpy()[coordinates]*u.deg,
                                       covered.dec_deg.to_numpy()[coordinates]*u.deg).separation(
            SkyCoord(master.gaia_ra.to_numpy()[coordinates]*u.deg,
                     master.gaia_dec.to_numpy()[coordinates]*u.deg)).arcsec
    passed = available & (offsets <= 1e-6)
    if not passed.all():
        raise ValueError(f'{(~passed).sum()} master fields lack valid, complete {radius_arcsec:g}-arcsec Tractor coverage. Missing coverage is not a non-detection.')


def resolve_tractor(frame):
    required = {'release','brickid','objid','ra','dec','type', 'flux_g','flux_r','flux_z',
                'flux_ivar_g','flux_ivar_r','flux_ivar_z','nobs_g','nobs_r','nobs_z',
                'maskbits','fracflux_r'}
    if not required.issubset(frame):
        raise ValueError('Tractor catalogue missing '+', '.join(sorted(required.difference(frame.columns))))
    finite = np.isfinite(frame[['ra','dec']].to_numpy(float)).all(axis=1)
    frame = frame.loc[finite].copy()
    galactic_b = SkyCoord(frame.ra.to_numpy()*u.deg, frame.dec.to_numpy()*u.deg).galactic.b.deg
    north = frame.dec.ge(32.375).to_numpy() & (galactic_b > 0)
    resolved = (frame.release.eq(9011).to_numpy() & north) | (frame.release.isin([9010,9012]).to_numpy() & ~north)
    if 'brick_primary' in frame:
        resolved &= frame.brick_primary.astype(str).str.lower().isin(['true','t','1','1.0']).to_numpy()
    frame = frame.loc[resolved].copy()
    frame['brick_primary'] = True
    frame['type'] = frame.type.astype(str).str.strip()
    return frame


def find_tractor_pairs(master, paths, radius_arcsec, aperture_arcsec):
    """Search every covered master field; cache associations are re-keyed by Gaia ID."""
    parent_sky = SkyCoord(master.gaia_ra.to_numpy()*u.deg, master.gaia_dec.to_numpy()*u.deg)
    parent_lookup = pd.Index(master.source_id)
    nearest_separation = np.full(len(master), np.inf)
    nearest_g = np.full(len(master), np.nan)
    nearest_r = np.full(len(master), np.nan)
    aperture_sources = []
    n_input = n_resolved = 0
    for file_index, path in enumerate(paths, start=1):
        for raw in read_chunks(path, TRACTOR_COLUMNS):
            n_input += len(raw)
            local = resolve_tractor(raw)
            n_resolved += len(local)
            if local.empty:
                continue
            source_sky = SkyCoord(local.ra.to_numpy()*u.deg, local.dec.to_numpy()*u.deg)
            if 'parent_gaia_source_id' in local:
                parent_indices = parent_lookup.get_indexer(local.parent_gaia_source_id)
                found = parent_indices >= 0
                local = local.loc[found].reset_index(drop=True)
                parent_indices = parent_indices[found]
                separation = source_sky[found].separation(parent_sky[parent_indices]).arcsec
            elif 'qso_parent_index' in local:
                raise ValueError('Parent-keyed Tractor caches need parent_gaia_source_id for safe re-keying. Supply unassociated survey rows or exact Gaia parent identifiers.')
            else:
                parent_indices, source_indices, separation, _ = search_around_sky(parent_sky, source_sky, radius_arcsec*u.arcsec)
                local = local.iloc[source_indices].reset_index(drop=True)
                separation = separation.arcsec
            within_radius = np.isfinite(separation) & (separation <= radius_arcsec)
            local = local.loc[within_radius].reset_index(drop=True)
            parent_indices = parent_indices[within_radius]
            separation = separation[within_radius]
            if local.empty:
                continue
            local['qso_parent_index'] = master.qso_parent_index.to_numpy()[parent_indices]
            local['target_sep_arcsec'] = separation
            order = np.lexsort((separation, parent_indices))
            first = np.r_[True, np.diff(parent_indices[order]) != 0]
            nearest_rows = order[first]
            improving = separation[nearest_rows] < nearest_separation[parent_indices[nearest_rows]]
            nearest_rows = nearest_rows[improving]
            nearest_separation[parent_indices[nearest_rows]] = separation[nearest_rows]
            nearest_g[parent_indices[nearest_rows]] = local.flux_g.to_numpy()[nearest_rows]
            nearest_r[parent_indices[nearest_rows]] = local.flux_r.to_numpy()[nearest_rows]
            aperture_sources.append(local.loc[separation <= aperture_arcsec].copy())
        if file_index % 20 == 0 or file_index == len(paths):
            print(f'Tractor files: {file_index}/{len(paths)}; raw rows: {n_input:,}', flush=True)
    if not aperture_sources:
        raise ValueError('No Tractor sources within the search aperture in the supplied covered fields')
    local = pd.concat(aperture_sources, ignore_index=True)
    identity = ['qso_parent_index','release','brickid','objid']
    duplicate = local.duplicated(identity)
    n_duplicates = int(duplicate.sum())
    local = local.loc[~duplicate].copy()
    finite_flux = np.isfinite(local.flux_r)
    counts = local.loc[finite_flux].groupby('qso_parent_index').size()
    multi_ids = counts.loc[counts.ge(2)].index
    local = local.loc[local.qso_parent_index.isin(multi_ids)].copy()
    parent_by_index = master.set_index('qso_parent_index',drop=False)
    channels = []
    for parent_id, group in local.groupby('qso_parent_index',sort=True):
        row, _ = channel_row(parent_by_index.loc[parent_id].to_dict(), group, aperture_arcsec)
        if row is not None:
            channels.append(row)
    channel = pd.DataFrame(channels)
    reference_colour = np.full(len(master), np.nan)
    usable = np.isfinite(nearest_g) & np.isfinite(nearest_r) & (nearest_g > 0) & (nearest_r > 0)
    reference_colour[usable] = -2.5*np.log10(nearest_g[usable]/nearest_r[usable])
    references = pd.DataFrame({'qso_parent_index':master.qso_parent_index,
        'z':master.z,'nearest_g_minus_r':reference_colour})
    references = references.loc[np.isfinite(references.z) & np.isfinite(references.nearest_g_minus_r)].copy()
    summary = dict(master_rows=len(master), raw_tractor_rows=n_input, resolved_tractor_rows=n_resolved,
                   duplicate_aperture_rows_removed=n_duplicates,
                   covered_fields_with_sources=int(np.isfinite(nearest_separation).sum()),
                   nearest_colour_references=len(references), initial_pair_systems=len(channel),
                   initial_psf_pairs=int(channel.type_pair.eq('PSF-PSF').sum()))
    return channel, local.reset_index(drop=True), references, summary
