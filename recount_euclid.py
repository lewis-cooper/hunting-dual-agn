"""Recount Euclid footprint, image availability and reviewed multiplicity.

All inputs are local. Footprint membership is recalculated from the supplied
MOCs. Image usability is measured from cached VIS pixels around both component
positions. Multiplicity is a supplied visual classification, not an automated
claim of two AGN. No images or catalogues are downloaded.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from dual_agn.candidate_tiers import boolean_flag


def read_moc_intervals(path):
    """Reuse the sky-coverage plot's NUNIQ-to-NESTED interval calculation.

    NUNIQ encoding: https://www.ivoa.net/documents/MOC/20220727/REC-moc-2.0-20220727.pdf
    """
    with fits.open(path) as hdus:
        table = hdus[1]
        if table.header['ORDERING'] != 'NUNIQ' or table.header['COORDSYS'] != 'C':
            raise ValueError(f'Expected a celestial NUNIQ MOC: {path}')
        max_order = int(table.header['MOCORDER'])
        unique_pixels = np.asarray(table.data.field(0), dtype=np.int64)
    if len(unique_pixels) == 0 or (unique_pixels < 4).any():
        raise ValueError(f'Empty or invalid MOC: {path}')
    order = (np.floor(np.log2(unique_pixels)).astype(int) - 2) // 2
    if (order > max_order).any():
        raise ValueError('MOCORDER is smaller than an encoded cell order')
    pixel = unique_pixels - (4 << (2 * order))
    factor = 1 << (2 * (max_order - order))
    sort_index = np.argsort(pixel * factor)
    starts = (pixel * factor)[sort_index]
    stops = ((pixel + 1) * factor)[sort_index]
    if (stops[:-1] > starts[1:]).any():
        raise ValueError(f'Overlapping MOC cells: {path}')
    return max_order, starts, stops


def in_moc(moc, ra_deg, dec_deg):
    from healpy import ang2pix

    order, starts, stops = moc
    pixels = ang2pix(1 << order, ra_deg, dec_deg, nest=True, lonlat=True)
    indices = np.searchsorted(starts, pixels, side='right') - 1
    return (indices >= 0) & (pixels < stops[indices])


def measure_vis_pixels(record, array_path, primary, secondary):
    """A partial image is usable if both positions have data in a 3x3 patch."""
    with np.load(array_path) as arrays:
        image = np.array(arrays['VIS'], dtype=float)
    if image.ndim != 2 or image.shape[0] != image.shape[1]:
        raise ValueError('Expected a square aligned VIS image')
    grid = record['grid']
    if not grid['north_up'] or not grid['east_left']:
        raise ValueError('Unexpected cached image orientation')
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [(image.shape[1] + 1) / 2, (image.shape[0] + 1) / 2]
    scale_deg = float(grid['pixel_arcsec']) / 3600
    wcs.wcs.cdelt = [-scale_deg, scale_deg]
    wcs.wcs.crval = [record['mid_ra'], record['mid_dec']]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    x, y = wcs.world_to_pixel_values([primary.ra.deg, secondary.ra.deg],
                                    [primary.dec.deg, secondary.dec.deg])
    valid = np.isfinite(image) & (image != 0)
    patches = []
    for xx, yy in zip(x, y):
        ix, iy = int(round(xx)), int(round(yy))
        inside = 0 <= ix < image.shape[1] and 0 <= iy < image.shape[0]
        patch = valid[max(0, iy - 1):iy + 2, max(0, ix - 1):ix + 2] if inside else np.empty(0)
        patches.append(int(patch.sum()))
    return {'vis_valid_pixels': int(valid.sum()), 'primary_patch_valid_pixels': patches[0],
            'secondary_patch_valid_pixels': patches[1],
            'usable_vis': bool(all(count > 0 for count in patches)),
            'vis_image_sha256': hashlib.sha256(array_path.read_bytes()).hexdigest()}


def summarize_audit(audit, candidates):
    joined = candidates.merge(audit, on='candidate_id', validate='one_to_one')
    if len(joined) != len(candidates):
        raise ValueError('Incomplete Euclid audit')
    usable = joined.usable_vis.astype(bool)
    resolved = boolean_flag(joined, 'selected_by_resolved')
    unresolved = boolean_flag(joined, 'selected_by_unresolved')
    imaging = boolean_flag(joined, 'selected_by_imaging')
    shared = imaging & (resolved | unresolved)
    if (resolved & unresolved).any():
        raise ValueError('A system cannot be both resolved and unresolved in the Gaia route')
    counts = {
        'NeuclidCovered': int(joined.covered_both.sum()),
        'NeuclidVIS': int(usable.sum()),
        'NeuclidTwo': int((usable & joined.vis_class.eq('resolved_pair')).sum()),
        'NeuclidBlend': int((usable & joined.vis_class.eq('potential_blend')).sum()),
        'NeuclidResolved': int((usable & resolved).sum()),
        'NeuclidBlended': int((usable & unresolved).sum()),
        'NeuclidImaging': int((usable & imaging).sum()),
        'NeuclidShared': int((usable & shared).sum()),
        'NeuclidMain': int((usable & joined.tier.eq('main')).sum()),
        'NeuclidAdditional': int((usable & joined.tier.eq('additional')).sum()),
    }
    if counts['NeuclidTwo'] + counts['NeuclidBlend'] != counts['NeuclidVIS']:
        raise ValueError('Usable VIS systems need explicit, exhaustive visual classifications')
    if counts['NeuclidMain'] + counts['NeuclidAdditional'] != counts['NeuclidVIS']:
        raise ValueError('Euclid tier counts do not reconcile')
    if (counts['NeuclidResolved'] + counts['NeuclidBlended'] + counts['NeuclidImaging']
            - counts['NeuclidShared']) != counts['NeuclidVIS']:
        raise ValueError('Euclid route counts do not reconcile after deduplication')
    return counts


def recount(catalogue_path, wide_moc, deep_moc, grid_records_path, classifications_path,
            output_dir, position_tolerance_arcsec=0.1):
    candidates = pd.read_csv(catalogue_path, dtype=str, keep_default_na=False)
    classifications = pd.read_csv(classifications_path, dtype=str, keep_default_na=False)
    if not candidates.candidate_id.is_unique or not classifications.candidate_id.is_unique:
        raise ValueError('Candidate and visual-review identifiers must be unique')
    review = classifications.set_index('candidate_id')
    records = json.loads(grid_records_path.read_text())
    by_id = {record['candidate_id']: record for record in records}
    if len(by_id) != len(records):
        raise ValueError('Duplicate cached image identities')
    mocs = [read_moc_intervals(path) for path in [wide_moc, deep_moc]]
    covered = {}
    coordinates = {}
    for role in ['primary', 'secondary']:
        ra = pd.to_numeric(candidates[f'{role}_ra_deg'], errors='raise').to_numpy()
        dec = pd.to_numeric(candidates[f'{role}_dec_deg'], errors='raise').to_numpy()
        if not (np.isfinite(ra).all() and np.isfinite(dec).all()):
            raise ValueError(f'Nonfinite {role} coordinates')
        coordinates[role] = SkyCoord(ra * u.deg, dec * u.deg)
        covered[role] = in_moc(mocs[0], ra, dec) | in_moc(mocs[1], ra, dec)
    audit_rows = []
    for index, row in enumerate(candidates.itertuples(index=False)):
        record = dict(candidate_id=row.candidate_id, covered_primary=bool(covered['primary'][index]),
                      covered_secondary=bool(covered['secondary'][index]),
                      covered_both=bool(covered['primary'][index] & covered['secondary'][index]),
                      usable_vis=False, vis_class='outside_footprint',
                      vis_valid_pixels=0, primary_patch_valid_pixels=0,
                      secondary_patch_valid_pixels=0, vis_image_sha256='', grid_row='',
                      maximum_position_offset_arcsec='', classification_basis='')
        if record['covered_both']:
            if row.candidate_id not in by_id:
                raise ValueError(f'Covered system has no checked image cache: {row.candidate_id}')
            cached = by_id[row.candidate_id]
            offsets = []
            for role in ['primary', 'secondary']:
                previous = SkyCoord(float(cached[f'{role}_ra_deg']) * u.deg,
                                    float(cached[f'{role}_dec_deg']) * u.deg)
                offsets.append(float(previous.separation(coordinates[role][index]).arcsec))
            if max(offsets) > position_tolerance_arcsec:
                raise ValueError(f'Stale cached component positions: {row.candidate_id}')
            path = Path(cached['array_path'])
            if not path.exists():
                path = grid_records_path.parent / path.name
            if not path.exists():
                raise FileNotFoundError(path)
            record.update(measure_vis_pixels(cached, path, coordinates['primary'][index],
                                             coordinates['secondary'][index]))
            record.update(grid_row=cached['row_number'], maximum_position_offset_arcsec=max(offsets),
                          vis_class='no_usable_vis')
            if record['usable_vis']:
                if row.candidate_id not in review.index:
                    raise ValueError(f'Usable VIS system lacks a visual review: {row.candidate_id}')
                decision = review.loc[row.candidate_id]
                if decision.vis_class not in ['resolved_pair', 'potential_blend']:
                    raise ValueError(f'Invalid visual classification: {row.candidate_id}')
                record.update(vis_class=decision.vis_class, classification_basis=decision.review_basis)
        audit_rows.append(record)
    audit = pd.DataFrame(audit_rows)
    counts = summarize_audit(audit, candidates)
    output_dir.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output_dir / 'euclid_coverage.csv', index=False)
    pd.DataFrame([{'macro': name, 'count': value} for name, value in counts.items()]).to_csv(
        output_dir / 'euclid_counts.csv', index=False)
    (output_dir / 'euclid_macros.tex').write_text('\n'.join(
        rf'\newcommand{{\{name}}}{{\prov{{{value}}}\xspace}}' for name, value in counts.items()) + '\n')
    summary = dict(candidate_count=len(candidates), counts=counts,
                   catalogue_sha256=hashlib.sha256(catalogue_path.read_bytes()).hexdigest(),
                   footprint_sha256={name: hashlib.sha256(path.read_bytes()).hexdigest()
                                     for name, path in [('wide', wide_moc), ('deep', deep_moc)]},
                   visual_review_sha256=hashlib.sha256(classifications_path.read_bytes()).hexdigest(),
                   grid_records_sha256=hashlib.sha256(grid_records_path.read_bytes()).hexdigest(),
                   usability='Both component positions have at least one finite, nonzero pixel in a 3x3 patch of the aligned VIS cutout.',
                   classification='Supplied visual review; source multiplicity is not confirmation of two AGN.')
    (output_dir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(counts, indent=2))
    print('Euclid audit:', output_dir.resolve())
    return audit, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', type=Path, required=True)
    parser.add_argument('--wide-moc', type=Path, required=True)
    parser.add_argument('--deep-moc', type=Path, required=True)
    parser.add_argument('--grid-records', type=Path, required=True)
    parser.add_argument('--classifications', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--position-tolerance-arcsec', type=float, default=0.1)
    args = parser.parse_args()
    recount(args.catalogue, args.wide_moc, args.deep_moc, args.grid_records,
            args.classifications, args.output_dir, args.position_tolerance_arcsec)


if __name__ == '__main__':
    main()
