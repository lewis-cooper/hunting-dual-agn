#!/usr/bin/env python3
"""Recover native Tractor fluxes for Gaia controls omitted by the PSF-only cones."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from control_catalogues import cone_sql, public_query
from query_parent_control_cones import normalize_columns
from score_controls import HEMISPHERE


def query(args):
    root = args.source
    library = pd.read_parquet(args.library)
    backgrounds = library.loc[library.route.eq('blended')].copy()
    eligible = ~np.isfinite(backgrounds.pm_sig) | backgrounds.pm_sig.lt(3)
    missing = backgrounds.photometric_system.isna()
    identified=backgrounds.background_source_id.fillna('').ne('')
    if (identified&~np.isfinite(backgrounds.ra+backgrounds.dec)).any():
        raise ValueError('A real background source lacks finite coordinates')
    targets = backgrounds.loc[eligible & missing & identified].drop_duplicates('background_source_id')
    targets = targets.assign(target_id=targets.background_source_id)
    print(f'Gaia controls needing an all-morphology Tractor check: {len(targets)}', flush=True)
    jobs = [targets.iloc[start:start + args.batch_size] for start in range(0, len(targets), args.batch_size)]

    def fetch(group):
        sql = cone_sql(group, 'tractor', radius=.5)
        return normalize_columns(public_query(sql, root / 'query_cache'), 'tractor')

    pieces = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for number, frame in enumerate(pool.map(fetch, jobs), 1):
            pieces.append(frame)
            if number % 25 == 0 or number == len(jobs):
                print(f'Native-flux batches {number}/{len(jobs)}', flush=True)
    result = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()
    if len(result):
        # A source must have exactly one native measurement and no competing Gaia
        # association. Shared Tractor measurements cannot be added independently.
        from build_parent_control_trials import resolve_photometric_area
        result = resolve_photometric_area(result).reset_index(drop=True)
        gaia = pd.read_parquet(root / 'gaia_neighbours.parquet').drop_duplicates('source_id')
        source_coordinates = SkyCoord(gaia.ra.to_numpy(float)*u.deg, gaia.dec.to_numpy(float)*u.deg)
        coordinates = SkyCoord(result.ra.to_numpy(float)*u.deg, result.dec.to_numpy(float)*u.deg)
        source_index, measurement_index, _, _ = coordinates.search_around_sky(source_coordinates, .5*u.arcsec)
        association_counts = np.bincount(measurement_index, minlength=len(result))
        nearest, _, _ = coordinates.match_to_catalog_sky(source_coordinates)
        requested = result.target_id.str.removeprefix('gaia_').to_numpy()
        result['unique_gaia_association'] = (association_counts == 1) & (
            gaia.iloc[nearest].source_id.astype(str).to_numpy() == requested)
        result['unique_native_association'] = ~result.duplicated('target_id', keep=False)
        result['photometric_system'] = result.release.map(HEMISPHERE)
    result.to_parquet(args.output, index=False, compression='zstd')
    args.output.with_suffix('.json').write_text(json.dumps(dict(
        requested_sources=len(targets), returned_measurements=len(result),
        radius_arcsec=.5, morphology='all Tractor types; brick_primary; DR9 native area',
        independent='one native measurement and exactly one Gaia source within 0.5 arcsec',
        purpose='repair query omission, without filling truly missing fluxes with zero'), indent=2)+'\n')
    print(f'Wrote {len(result)} native measurements to {args.output}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--batch-size', type=int, default=12)
    query(parser.parse_args())


if __name__ == '__main__':
    main()
