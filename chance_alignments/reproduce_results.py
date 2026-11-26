#!/usr/bin/env python3
"""Reproduce method pass rates from saved controls, without models or network access."""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

ENGINE = Path(__file__).resolve().parent/'engine'


def validate_controls(root):
    inputs = pd.read_csv(root/'input_pairs.csv')
    columns = ['trial_id', 'input_pair_id', 'route', 'annulus', 'mock_weight']
    frames = [pd.read_parquet(root/name, columns=columns) for name in
              ['component_control_scores.parquet', 'blend_control_scores.parquet']]
    controls = pd.concat(frames, ignore_index=True)
    if not controls.trial_id.is_unique:
        raise ValueError('Duplicate trial IDs would count a control more than once.')
    if not (np.isfinite(controls.mock_weight) & controls.mock_weight.gt(0)).all():
        raise ValueError('Control weights must be finite and positive.')
    np.testing.assert_allclose(
        controls.groupby(['input_pair_id', 'annulus']).mock_weight.sum(), 1., atol=1e-12)
    for route, group in inputs.loc[inputs.latitude_bin.notna()].groupby('route'):
        for annulus in ['3-10', '10-30']:
            observed = controls.loc[controls.route.eq(route) & controls.annulus.eq(annulus)]
            if set(observed.input_pair_id) != set(group.pair_id):
                raise ValueError(f'{route}, {annulus}: inputs are missing or unexpected.')
    print(f'Validated {len(controls):,} mock rows; source identities and unit per-input '
          'weights are retained. Mock rows are not independent observations.', flush=True)


def compare_reference(output, reference):
    for name in ['repaired_summary_by_route.csv', 'repaired_summary_by_latitude.csv',
                 'parent_densities_by_latitude.csv', 'parent_densities_by_magnitude.csv',
                 'control_magnitude_matching.csv']:
        expected = pd.read_csv(reference/name)
        actual = pd.read_csv(output/name)
        pd.testing.assert_frame_equal(actual, expected, check_dtype=False,
                                      check_exact=False, rtol=1e-10, atol=1e-12)
    print('All five numerical result tables reproduce the reference values.', flush=True)


def reproduce(data, output, bootstrap=1000, seed=20261121, verify_reference=False, density_source=None):
    data = data.resolve()
    output = output.resolve()
    if output == data or data.is_relative_to(output) or output.is_relative_to(data):
        raise ValueError('Use a separate output folder, outside the supplied control data.')
    density_source = Path(density_source).resolve() if density_source else data/'field_density'
    required_density = ['master_targets.csv', 'cone_coverage.csv',
                        'gaia_neighbours.parquet', 'query_configuration.json']
    missing_density = [str(density_source/name) for name in required_density if not (density_source/name).is_file()]
    if missing_density:
        raise FileNotFoundError('Matched-field replay needs complete close/outer Gaia cones with exact '
            'central IDs and coverage. Older control archives lack these inputs. Supply '
            '--density-source /path/to/fields_60000 or add the exported field_density_inputs '
            'directory to control_data as field_density. Missing files:\n'+'\n'.join(missing_density))
    output.mkdir(parents=True, exist_ok=True)
    for name in ['input_pairs.csv', 'pair_cut_factors.csv',
                 'component_control_scores.parquet', 'blend_control_scores.parquet']:
        source = data/'replay'/name
        if not source.is_file():
            raise FileNotFoundError(f'Unpack the optional control-data ZIP first: {source}')
        shutil.copy2(source, output/name)
    validate_controls(output)
    commands = [
        ('analyze_field_density.py', ['--source', density_source, '--output', output,
                                     '--bootstrap', bootstrap, '--seed', seed,
                                     '--export-inputs', output/'field_density_inputs']),
        ('summarize_repaired_controls.py', ['--output', output, '--bootstrap', bootstrap, '--seed', seed]),
        ('summarize_parent_densities.py', ['--source', data/'density', '--output', output,
                                         '--bootstrap', bootstrap, '--seed', seed]),
        ('plot_control_matching.py', ['--output', output,
                                     '--blended-trials', output/'blend_control_scores.parquet']),
        ('write_precision_report.py', ['--output', output]),
    ]
    for script, arguments in commands:
        print('Running', script, flush=True)
        subprocess.run([sys.executable, '-B', str(ENGINE/script),
                        *map(str, arguments)], check=True)
    if verify_reference:
        if bootstrap != 1000 or seed != 20261121:
            raise ValueError('Reference intervals use bootstrap=1000 and seed=20261121.')
        compare_reference(output, data/'reference')
    print('Report:', output/'RESULTS.md', flush=True)
    print('These are method false-positive acceptance rates, not selected-sample contamination.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True, help='Unpacked control_data folder')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=20261121)
    parser.add_argument('--verify-reference', action='store_true')
    parser.add_argument('--density-source', type=Path, help='Complete Gaia cones for older replay archives')
    args = parser.parse_args()
    if args.bootstrap < 100:
        parser.error('Use at least 100 replicates; the reference calculation uses 1000.')
    reproduce(args.data, args.output, args.bootstrap, args.seed, args.verify_reference, args.density_source)
