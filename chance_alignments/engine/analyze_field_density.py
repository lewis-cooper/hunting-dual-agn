"""Compare detected Gaia companions with outer densities in the same fields.

The geometric prediction assumes complete close-pair detection (C=1). The
observed/predicted ratio is neither completeness nor a contamination fraction.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil

import numpy as np
import pandas as pd


def prepare_fields(source):
    fields = pd.read_csv(source/'master_targets.csv', dtype={'target_id': 'string', 'source_id': 'string'},
                         float_precision='round_trip')
    coverage = pd.read_csv(source/'cone_coverage.csv', dtype={'target_id': 'string'},
                           float_precision='round_trip').set_index('target_id')
    if not len(fields) or fields.target_id.isna().any() or fields.source_id.isna().any():
        raise ValueError('Supply nonempty fields with exact target and central Gaia identities')
    if not fields.target_id.is_unique or not coverage.index.is_unique:
        raise ValueError('Field and coverage identities must be unique')
    covered = coverage.reindex(fields.target_id)
    complete = covered.gaia_query_complete.eq(True).to_numpy(copy=True)
    for coordinate in ('ra', 'dec'):
        complete &= np.isfinite(fields[coordinate]).to_numpy()
        complete &= np.isclose(fields[coordinate], covered[coordinate], rtol=0, atol=1e-10)
    if not complete.all():
        raise ValueError('Every field needs a successful Gaia cone at its exact coordinates')
    query = json.loads((source/'query_configuration.json').read_text())
    digest = hashlib.sha256((source/'master_targets.csv').read_bytes()).hexdigest()
    if query['radius_arcsec'] < 60 or query['parent_targets_sha256'] != digest:
        raise ValueError('Cones must cover 60 arcsec and the same sampled parent fields')
    weights = fields.sampling_weight.to_numpy(float)
    if not (np.isfinite(weights) & (weights > 0)).all() or fields.sampling_cell.isna().any():
        raise ValueError('Sampling cells and finite positive inclusion weights are required')
    neighbours = pd.read_parquet(source/'gaia_neighbours.parquet')
    if pd.api.types.is_float_dtype(neighbours.source_id):
        raise ValueError('Gaia source IDs must be exact integers or strings, never floating point')
    for name in ('source_id', 'target_id'):
        neighbours[name] = neighbours[name].astype('string')
        if neighbours[name].isna().any():
            raise ValueError(f'Missing neighbour {name}')
    if neighbours.duplicated(['target_id', 'source_id']).any():
        raise ValueError('Duplicate Gaia identities within a field')
    if not neighbours.target_id.isin(fields.target_id).all():
        raise ValueError('Unknown field identity in Gaia cones')
    central_ids = neighbours.target_id.map(fields.set_index('target_id').source_id)
    central = neighbours.source_id.eq(central_ids)
    central_count = neighbours.loc[central].groupby('target_id').size().reindex(fields.target_id, fill_value=0)
    if not central_count.eq(1).all():
        raise ValueError('Every field needs exactly one central Gaia detection for identity exclusion')
    separation = neighbours.separation_arcsec.to_numpy(float)
    if not (np.isfinite(separation) & (separation >= 0) & (separation <= 60.000001)).all():
        raise ValueError('Neighbour distances must be finite and within the 60-arcsec cones')
    return fields.reset_index(drop=True), neighbours.loc[~central].copy()


def export_replay_inputs(source, destination):
    """Export raw close/outer detections, including exact central IDs and coverage."""
    source, destination = source.resolve(), destination.resolve()
    if source == destination or source.is_relative_to(destination) or destination.is_relative_to(source):
        raise ValueError('Use a separate replay input directory outside the source cones')
    prepare_fields(source)
    destination.mkdir(parents=True, exist_ok=True)
    for name in ('master_targets.csv', 'cone_coverage.csv'):
        shutil.copy2(source/name, destination/name)
    columns = ['target_id', 'source_id', 'separation_arcsec']
    pd.read_parquet(source/'gaia_neighbours.parquet', columns=columns).to_parquet(
        destination/'gaia_neighbours.parquet', index=False, compression='zstd')
    metadata = dict(radius_arcsec=60,
        parent_targets_sha256=hashlib.sha256((destination/'master_targets.csv').read_bytes()).hexdigest(),
        central_sources_included=True, purpose='Matched-field density replay; no selection or scoring')
    (destination/'query_configuration.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print('Density replay inputs:', destination, flush=True)


def field_counts(fields, neighbours, mask):
    return neighbours.loc[mask].groupby('target_id').size().reindex(fields.target_id, fill_value=0).to_numpy(float)


def finite_interval(values):
    finite = values[np.isfinite(values)]
    return np.quantile(finite, [.025, .975]) if len(finite) else (np.nan, np.nan)


def compare_counts(fields, neighbours, radii, bootstrap, seed):
    radii = np.asarray(radii, float)
    if not len(radii) or not (np.isfinite(radii) & (radii > 0) & (radii <= 10)).all():
        raise ValueError('Close radii must be positive and at most 10 arcsec, outside the outer annuli')
    if bootstrap < 100:
        raise ValueError('Use at least 100 bootstrap draws')
    weights = fields.sampling_weight.to_numpy(float)
    field_table = fields[['target_id', 'source_id', 'sampling_weight', 'sampling_cell']].copy()
    close_counts = {}
    for radius in radii:
        close_counts[radius] = field_counts(fields, neighbours, neighbours.separation_arcsec.le(radius))
        field_table[f'neighbours_within_{radius:g}_arcsec'] = close_counts[radius].astype(int)
    specifications = []
    for lower, upper in [(10., 60.), (10., 30.), (30., 60.)]:
        outer = neighbours.separation_arcsec.gt(lower) & neighbours.separation_arcsec.le(upper)
        outer_counts = field_counts(fields, neighbours, outer)
        field_table[f'neighbours_{lower:g}_{upper:g}_arcsec'] = outer_counts.astype(int)
        for radius in radii:
            geometry = radius**2/(upper**2-lower**2)
            specifications.append((lower, upper, radius, close_counts[radius], outer_counts, geometry))
    observed = np.column_stack([row[3] for row in specifications])
    expected = np.column_stack([row[4]*row[5] for row in specifications])
    weighted_observed, weighted_expected = weights @ observed, weights @ expected
    cells, names = pd.factorize(fields.sampling_cell)
    population = np.bincount(cells, weights=weights, minlength=len(names))
    rng = np.random.default_rng(seed)
    observed_draws = np.empty((bootstrap, len(specifications)))
    expected_draws = np.empty_like(observed_draws)
    for draw in range(bootstrap):
        multiplier = rng.exponential(size=len(fields))
        mass = np.bincount(cells, weights=multiplier*weights, minlength=len(names))
        multiplier *= population[cells]/mass[cells]
        observed_draws[draw] = (weights*multiplier) @ observed
        expected_draws[draw] = (weights*multiplier) @ expected
    ratio_draws = np.divide(observed_draws, expected_draws,
                           out=np.full_like(observed_draws, np.nan), where=expected_draws > 0)
    summaries = []
    for column, (lower, upper, radius, close, outer, geometry) in enumerate(specifications):
        ratio = weighted_observed[column]/weighted_expected[column] if weighted_expected[column] > 0 else np.nan
        ratio_low, ratio_high = finite_interval(ratio_draws[:, column])
        observed_low, observed_high = finite_interval(observed_draws[:, column])
        expected_low, expected_high = finite_interval(expected_draws[:, column])
        summaries.append(dict(outer_lower_arcsec=lower, outer_upper_arcsec=upper,
            close_radius_arcsec=radius, n_fields=len(fields), parent_population=weights.sum(),
            observed_sample_companions=int(close.sum()), observed_sample_fields_with_companions=int((close > 0).sum()),
            outer_sample_neighbours=int(outer.sum()), expected_sample_companions_C1=float(outer.sum()*geometry),
            weighted_observed_companions=weighted_observed[column],
            weighted_observed_lower_95=observed_low, weighted_observed_upper_95=observed_high,
            weighted_expected_companions_C1=weighted_expected[column],
            weighted_expected_lower_95=expected_low, weighted_expected_upper_95=expected_high,
            observed_to_C1_ratio=ratio, ratio_lower_95=ratio_low, ratio_upper_95=ratio_high,
            C1_within_ratio_interval=bool(ratio_low <= 1 <= ratio_high)))
    return pd.DataFrame(summaries), field_table


def write_report(summary, output):
    main = summary.loc[summary.outer_lower_arcsec.eq(10) & summary.outer_upper_arcsec.eq(60)]
    lines = ['## Matched close-neighbour and 10–60 arcsec density comparison', '',
        f'The same {int(main.n_fields.iloc[0]):,} Gaia-centred quasar fields supply the '
        'close counts and outer densities. Central Gaia sources are excluded by exact identity. '
        'Fields with zero neighbours remain included. No photometric, PM or final-candidate cut is applied.', '',
        'The geometric prediction is N_outer × theta²/(60² − 10²). It assumes complete '
        'close-pair detection (C=1). Sample counts and inclusion-weighted population estimates '
        'are different quantities and are reported separately.', '',
        '| Radius [arcsec] | Detected in sampled fields | Predicted in sampled fields (C=1) | '
        'Weighted detected | Weighted predicted (C=1) | Weighted observed/predicted (95% interval) |',
        '| ---: | ---: | ---: | ---: | ---: | --- |']
    for row in main.itertuples():
        lines.append(f'| {row.close_radius_arcsec:g} | {row.observed_sample_companions} | '
            f'{row.expected_sample_companions_C1:.2f} | {row.weighted_observed_companions:.1f} | '
            f'{row.weighted_expected_companions_C1:.1f} | {row.observed_to_C1_ratio:.3f} '
            f'({row.ratio_lower_95:.3f}–{row.ratio_upper_95:.3f}) |')
    lines += ['', 'Outer-annulus check for a 1-arcsec aperture:', '']
    for row in summary.loc[summary.close_radius_arcsec.eq(1)].itertuples():
        lines.append(f'- {row.outer_lower_arcsec:g}–{row.outer_upper_arcsec:g} arcsec: '
                     f'{row.expected_sample_companions_C1:.2f} predicted in the sampled fields (C=1).')
    lines += ['', 'The observed/predicted ratio is not a measured unrelated fraction or detection '
        'completeness. N_detected = C × N_random(C=1) + N_associated_detected. Close-pair '
        'detection losses and associated companions are not independently constrained. '
        'Even agreement with the complete-detection prediction would not establish that '
        'most input pairs are chance alignments. Do not multiply these densities by pooled '
        'control pass rates to report actual contaminant counts.', '',
        'Intervals use a stratified field multiplier bootstrap preserving each sampling cell’s '
        'population mass. They exclude detection-completeness and population-transfer uncertainty; '
        'empty annuli cannot establish a zero upper bound on true density. Counts are incidences '
        'around parent anchors, which can include repeated physical systems.', '']
    (output/'FIELD_DENSITY.md').write_text('\n'.join(lines))


def plot_comparison(summary, output):
    os.environ.setdefault('MPLCONFIGDIR', str(output/'matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    main = summary.loc[summary.outer_lower_arcsec.eq(10) & summary.outer_upper_arcsec.eq(60)]
    finite = np.isfinite(main[['observed_to_C1_ratio', 'ratio_lower_95', 'ratio_upper_95']]).all(axis=1)
    rows = main.loc[finite]
    figure, axis = plt.subplots(figsize=(5.5, 4.5), layout='constrained')
    ratio = rows.observed_to_C1_ratio.to_numpy()
    axis.errorbar(rows.close_radius_arcsec, ratio,
        yerr=[np.maximum(0, ratio-rows.ratio_lower_95), np.maximum(0, rows.ratio_upper_95-ratio)],
        marker='o', capsize=3, color='#0072B2')
    axis.axhline(1, color='0.4', linestyle='--', linewidth=1)
    axis.set(xlabel='Close-neighbour radius [arcsec]',
        ylabel='Observed / geometric random prediction (C = 1)',
        title=f'Matched comparison in {int(main.n_fields.iloc[0]):,} quasar fields')
    axis.set_ylim(bottom=0)
    figure.savefig(output/'observed_vs_field_density.pdf')
    figure.savefig(output/'observed_vs_field_density.png', dpi=180)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--export-inputs', type=Path, help='Separate directory of replayable raw density inputs')
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=20261121)
    args = parser.parse_args()
    if args.source.resolve() == args.output.resolve():
        parser.error('Use an output directory separate from the source cones')
    fields, neighbours = prepare_fields(args.source)
    print(f'Read {len(fields):,} fields and {len(neighbours):,} non-central Gaia detections', flush=True)
    summary, counts = compare_counts(fields, neighbours, [.5, .75, 1., 1.5, 2., 3., 5., 10.], args.bootstrap, args.seed)
    args.output.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output/'matched_density_comparison.csv', index=False)
    counts.to_csv(args.output/'matched_density_field_counts.csv', index=False)
    write_report(summary, args.output)
    plot_comparison(summary, args.output)
    metadata = dict(seed=args.seed, bootstrap=args.bootstrap, n_fields=len(fields),
        inputs_sha256={name: hashlib.sha256((args.source/name).read_bytes()).hexdigest() for name in
                      ['master_targets.csv', 'cone_coverage.csv', 'gaia_neighbours.parquet', 'query_configuration.json']},
        interpretation='C=1 geometric comparison; no actual contamination or completeness inferred')
    (args.output/'matched_density_configuration.json').write_text(json.dumps(metadata, indent=2)+'\n')
    if args.export_inputs:
        export_replay_inputs(args.source, args.export_inputs)
    print(summary.loc[summary.close_radius_arcsec.eq(1), ['outer_lower_arcsec', 'outer_upper_arcsec',
        'observed_sample_companions', 'expected_sample_companions_C1', 'observed_to_C1_ratio']].to_string(index=False))
    print('Density outputs:', args.output, flush=True)


if __name__ == '__main__':
    main()
