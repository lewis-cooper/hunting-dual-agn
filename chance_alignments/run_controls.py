#!/usr/bin/env python3
"""Rebuild the general-method experiment using explicit inputs and production scorers."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

ENGINE = Path(__file__).resolve().parent/'engine'


def load_config(path):
    config = json.loads(path.read_text())
    for name in ['data_directory', 'output_directory', 'production_repo',
                 'model_directory', 'dustmaps_directory', 'master_catalogue', 'quaia_catalogue']:
        if config.get(name):
            value = Path(config[name]).expanduser()
            config[name] = (path.resolve().parent/value).resolve() if not value.is_absolute() else value
    config.setdefault('workers', 4)
    config.setdefault('bootstrap', 1000)
    if config['workers'] < 1 or config['bootstrap'] < 100:
        raise ValueError('Require workers >= 1 and bootstrap >= 100.')
    return config


def merge_components(previous, expanded):
    original = pd.read_parquet(previous/'component_control_scores.parquet')
    extension = pd.read_parquet(expanded/'component_control_scores.parquet')
    extension = extension.loc[extension.route.isin(['resolved', 'imaging_two_gaia'])]
    retained = original.loc[original.route.eq('imaging_one_gaia')]
    scores = pd.concat([extension, retained], ignore_index=True)
    if not scores.trial_id.is_unique:
        raise ValueError('Component batches contain duplicate trials.')
    np.testing.assert_allclose(scores.groupby(['input_pair_id', 'annulus']).mock_weight.sum(), 1.)
    scores.to_parquet(expanded/'component_control_scores.parquet', index=False, compression='zstd')


def run(config, stage, resume=False, dry_run=False, resample=False):
    data = config['data_directory']
    root = config['output_directory']
    repo = config['production_repo']
    if root == repo or root == data or root.is_relative_to(data) or data.is_relative_to(root):
        raise ValueError('Use a separate control output folder, distinct from source data and the repo root.')
    if not dry_run:
        root.mkdir(parents=True, exist_ok=True)
        version = json.loads((data/'reference/production_version.json').read_text())
        for name in ['dual_agn/one_arcsecond.py', 'dual_agn/photometry.py', 'dual_agn/redshift_uncertainty.py', 'dual_agn/score_selection.py']:
            actual = hashlib.sha256((repo/name).read_bytes()).hexdigest()
            if actual != version[name]:
                raise ValueError(f'Production scorer differs from the measured experiment: {name}')
        calibration = hashlib.sha256((data/'reference/threshold_calibration.json').read_bytes()).hexdigest()
        if calibration != version['calibration_sha256']:
            raise ValueError('Supply the calibrated threshold product used for this experiment.')
        fingerprint = json.dumps({key: str(value) if isinstance(value, Path) else value
                                  for key, value in config.items()}, sort_keys=True, indent=2)
        settings = root/'run_settings.json'
        if settings.exists() and settings.read_text() != fingerprint:
            raise ValueError('Configuration changed. Use a new output directory.')
        settings.write_text(fingerprint)

    def execute(label, script, arguments, expected=None):
        command = [sys.executable, '-B', str(ENGINE/script), *map(str, arguments)]
        if dry_run:
            print(label+':', ' '.join(command), flush=True)
            return
        if resume and expected and Path(expected).exists():
            print('Using completed step:', label, flush=True)
            return
        print('Running:', label, flush=True)
        with (root/(label+'.log')).open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)

    def copy(source, destination):
        if not dry_run:
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)

    source30, source60 = root/'fields_30000', root/'fields_60000'
    paired30, paired60 = root/'paired_30000', root/'paired_60000'
    coverage, legacy = root/'blend_coverage_extension', root/'blend_legacy_extension'
    common = ['--repo', repo, '--models', config['model_directory'],
              '--dustmaps', config['dustmaps_directory']]
    component_common = [*common, '--calibration', data/'reference/threshold_calibration.json']
    phases = ['prepare', 'fetch', 'components', 'blends', 'report'] if stage == 'all' else [stage]
    for phase in phases:
        if phase == 'prepare':
            for number, destination in [(30000, source30), (60000, source60)]:
                for name in ['input_pairs.csv', 'pair_cut_factors.csv']:
                    copy(data/'prepared'/name, destination/name)
                if resample:
                    if not config.get('master_catalogue') or not config.get('quaia_catalogue'):
                        raise ValueError('Resampling requires master_catalogue and quaia_catalogue paths.')
                    previous = data/'prepared' if number == 30000 else source30
                    execute(f'sample_{number}', 'prepare_matched_parent_controls.py',
                            ['--previous', previous, '--master', config['master_catalogue'],
                             '--quaia', config['quaia_catalogue'], '--output', destination,
                             '--n-parents', number], destination/'master_targets.csv')
                else:
                    copy(data/'prepared'/f'parent_fields_{number}.csv', destination/'master_targets.csv')
                    copy(data/'prepared'/f'parent_sampling_cells_{number}.csv',
                         destination/'parent_sampling_cells.csv')
        elif phase == 'fetch':
            for number, destination in [(30000, source30), (60000, source60)]:
                reuse = ['--reuse-cones', source30] if number == 60000 else []
                execute(f'query_{number}', 'query_parent_control_cones.py',
                        ['--output', destination, '--workers', config['workers'], *reuse],
                        destination/'query_configuration.json')
                execute(f'build_fields_{number}', 'build_parent_control_trials.py',
                        ['--output', destination], destination/'parent_control_trials.parquet')
        elif phase == 'components':
            for number, source, destination, routes in [
                (30000, source30, paired30, []),
                (60000, source60, paired60, ['--routes', 'resolved', 'imaging_two_gaia'])]:
                execute(f'pair_components_{number}', 'build_repaired_component_controls.py',
                        ['--source', source, '--output', destination, *component_common, *routes],
                        destination/'component_control_trials.parquet')
                execute(f'score_components_{number}', 'score_repaired_component_controls.py',
                        ['--output', destination, *component_common, '--workers', config['workers']],
                        destination/'component_control_scores.parquet')
            if not dry_run:
                merge_components(paired30, paired60)
        elif phase == 'blends':
            for number, source, paired in [(30000, source30, paired30), (60000, source60, paired60)]:
                execute(f'native_flux_{number}', 'query_missing_legacy_photometry.py',
                        ['--source', source, '--library', paired/'background_source_library.parquet',
                         '--output', source/'supplemental_legacy_photometry.parquet',
                         '--workers', config['workers']], source/'supplemental_legacy_photometry.parquet')
            batches = [
                ('main', source30, paired30, 100, 20261122, None),
                ('coverage', source60, coverage, 20, 20261122, coverage/'input_ids.csv'),
                ('legacy', source60, legacy, 100, 20261124, legacy/'input_ids.csv')]
            for label, source, destination, draws, seed, input_ids in batches:
                if label != 'main':
                    copy(paired60/'input_pairs.csv', destination/'input_pairs.csv')
                    copy(paired60/'background_source_library.parquet', destination/'background_source_library.parquet')
                    if not dry_run:
                        if label == 'coverage':
                            audit = pd.read_csv(paired30/'blended_matching_audit.csv')
                            ids = audit.loc[audit.annulus.eq('3-10') &
                                ((audit.n_primary_donors == 0) | (audit.n_background_sources == 0)), 'input_pair_id']
                        else:
                            inputs = pd.read_csv(destination/'input_pairs.csv')
                            ids = inputs.loc[inputs.route.eq('blended') & inputs.latitude_bin.notna()
                                             & inputs.photometry_class.eq('legacy'), 'pair_id']
                        ids.rename('input_pair_id').to_csv(input_ids, index=False)
                selection = ['--input-ids', input_ids] if input_ids else []
                execute(f'pair_blends_{label}', 'build_repaired_blended_controls.py',
                        ['--source', source, '--output', destination, '--repo', repo,
                         '--dustmaps', config['dustmaps_directory'], '--draws-per-input', draws,
                         '--seed', seed, *selection], destination/'blended_control_trials.parquet')
                if not dry_run:
                    sys.path.insert(0, str(ENGINE))
                    from types import SimpleNamespace
                    from score_parent_blended_controls import selected_trials
                    subset = selected_trials(SimpleNamespace(output=destination,
                        draws_per_input=20 if label == 'main' else draws,
                        legacy_draws_per_input=100 if label == 'main' else 0))
                    subset.to_parquet(destination/'photometry_trials.parquet', index=False)
                execute(f'query_fluxes_{label}', 'query_parent_blend_photometry.py',
                        ['--output', destination, '--trials', destination/'photometry_trials.parquet',
                         '--workers', config['workers']], destination/'independent_blend_photometry.csv')
                allocation = ['--draws-per-input', 20, '--legacy-draws-per-input', 100] if label == 'main' else []
                execute(f'score_blends_{label}', 'score_parent_blended_controls.py',
                        ['--output', destination, *common, '--workers', config['workers'], *allocation],
                        destination/'blend_control_scores.parquet')
            execute('assemble_batches', 'assemble_precision_results.py',
                    ['--output', paired60, '--blended', paired30, '--extension', coverage,
                     '--legacy-batch', legacy, '--blended-score-file', 'blend_control_scores.parquet',
                     '--draws-per-input', 20, '--legacy-draws-per-input', 100,
                     '--bootstrap', config['bootstrap']])
        elif phase == 'report':
            execute('summarize_rates', 'summarize_repaired_controls.py',
                    ['--output', paired60, '--bootstrap', config['bootstrap']])
            execute('summarize_density', 'summarize_parent_densities.py',
                    ['--source', source60, '--output', paired60, '--bootstrap', config['bootstrap']])
            execute('plot_matching', 'plot_control_matching.py',
                    ['--output', paired60, '--blended-trials', paired60/'blend_control_scores.parquet'])
            execute('write_report', 'write_precision_report.py', ['--output', paired60])
    print('Output directory:', root, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--stage', choices=['all', 'prepare', 'fetch', 'components', 'blends', 'report'], default='all')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--dry-run', action='store_true', help='Show commands without querying or scoring')
    parser.add_argument('--resample', action='store_true', help='Resample the original master rather than reuse frozen field identities')
    args = parser.parse_args()
    run(load_config(args.config), args.stage, args.resume, args.dry_run, args.resample)
