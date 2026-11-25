"""Read native catalogue extracts and make reproducible positional associations.

Column maps rename native fields; they do not contain analysis decisions.
Coverage is explicit so an unqueried field cannot become a non-detection.
"""

from pathlib import Path
import glob
import json

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky

from dual_agn.imaging_search.catalogues import read_chunks


def read_config(path):
    config = json.loads(Path(path).read_text())
    config['_directory'] = str(Path(path).resolve().parent)
    return config


def local_path(config, value):
    path = Path(value).expanduser()
    return path if path.is_absolute() else Path(config['_directory']) / path


def catalogue_files(config, specification):
    patterns = specification.get('paths', [specification.get('path')])
    if not patterns or any(pattern is None for pattern in patterns):
        raise ValueError('Every native catalogue needs path or paths')
    paths = []
    for pattern in patterns:
        matches = sorted(glob.glob(str(local_path(config, pattern))))
        if not matches:
            raise FileNotFoundError(f'No native catalogue matches {pattern}; no download attempted')
        paths.extend(Path(match) for match in matches)
    return list(dict.fromkeys(paths))


def sky(frame, ra='ra', dec='dec'):
    coordinates = frame[[ra, dec]].to_numpy(float)
    if not (np.isfinite(coordinates).all() and (np.abs(coordinates[:, 1]) <= 90).all()):
        raise ValueError(f'Invalid coordinates in {ra}, {dec}')
    return SkyCoord(coordinates[:, 0] * u.deg, coordinates[:, 1] * u.deg)


def require_coverage(targets, config, specification, radius_arcsec):
    """Require a successful cone, including empty cones, for every exact anchor."""
    if specification.get('full_catalogue') is True:
        return
    if not specification.get('coverage'):
        raise ValueError('Supply coverage CSV or explicitly set full_catalogue=true for a complete survey catalogue')
    coverage = pd.read_csv(local_path(config, specification['coverage']), dtype={'source_id': 'string'})
    required = {'source_id', 'ra_deg', 'dec_deg', 'query_radius_arcsec'}
    if not required.issubset(coverage) or not coverage.source_id.is_unique:
        raise ValueError('Coverage needs unique source_id, ra_deg, dec_deg, query_radius_arcsec')
    covered = coverage.set_index('source_id').reindex(targets.source_id).reset_index()
    coordinates = covered[['ra_deg', 'dec_deg']].to_numpy(float)
    finite = np.isfinite(coordinates).all(axis=1)
    separation = np.full(len(targets), np.inf)
    if finite.any():
        separation[finite] = sky(covered.loc[finite], 'ra_deg', 'dec_deg').separation(sky(targets.loc[finite])).arcsec
    complete = (covered.query_radius_arcsec.to_numpy(float) - separation) >= radius_arcsec - 1e-6
    if not complete.all():
        raise ValueError(f'{(~complete).sum()} fields lack complete {radius_arcsec:g}-arcsec coverage; missing coverage is not a non-detection')


def read_local_sources(targets, config, specification, radius_arcsec, identity, fields=None):
    """Stream native files, keeping only sources near the supplied coordinates."""
    target_sky = sky(targets)
    retained = []
    column_map = specification.get('columns', {})
    # The map is canonical -> native; unmatched native fields remain available.
    rename = {native: canonical for canonical, native in column_map.items()}
    selected = None if fields is None else [column_map.get(name, name) for name in dict.fromkeys(['ra', 'dec', *identity, *fields])]
    for path in catalogue_files(config, specification):
        for chunk in read_chunks(path, selected):
            chunk = chunk.rename(columns=rename)
            if not {'ra', 'dec', *identity}.issubset(chunk):
                raise ValueError(f'{path} needs ra, dec and identity {identity}, or a columns mapping')
            finite_position = np.isfinite(chunk[['ra', 'dec']].to_numpy(float)).all(axis=1)
            if not finite_position.all():
                raise ValueError(f'{path} contains nonfinite native source positions')
            if chunk.empty:
                continue
            _, source_index, _, _ = search_around_sky(target_sky, sky(chunk), radius_arcsec * u.arcsec)
            retained.append(chunk.iloc[np.unique(source_index)].copy())
    if not retained:
        # Retain schema even when the queried region has no sources.
        first = next(read_chunks(catalogue_files(config, specification)[0], selected)).rename(columns=rename)
        return first.iloc[:0].copy()
    sources = pd.concat(retained, ignore_index=True).drop_duplicates()
    if sources.duplicated(identity).any():
        raise ValueError(f'Conflicting native records for identity {identity}')
    return sources.reset_index(drop=True)


def nearest_sources(targets, sources, radius_arcsec, identity):
    """Return nearest native source per target, breaking ties by exact identity."""
    target_index, source_index, separation, _ = search_around_sky(sky(targets), sky(sources), radius_arcsec * u.arcsec)
    matches = sources.iloc[source_index].copy().reset_index(drop=True)
    matches['_target_index'] = target_index
    matches['match_sep_arcsec'] = separation.arcsec
    matches = matches.sort_values(['_target_index', 'match_sep_arcsec', *identity], kind='stable')
    return matches.drop_duplicates('_target_index').set_index('_target_index').reindex(range(len(targets)))
