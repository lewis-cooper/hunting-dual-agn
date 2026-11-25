"""Run chance-alignment controls from the repository's vetting outputs."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

from settings import DUSTMAP_ROOT, INPUT_ROOT, IMAGING_INPUT_ROOT, MODEL_ROOT, OUTPUT_ROOT, ROOT, WORK_ROOT


def resolve_config(path):
    supplied = json.loads(path.read_text())
    resolved = dict(input_directory=INPUT_ROOT.resolve(), imaging_input_directory=IMAGING_INPUT_ROOT.resolve(),
        products_directory=OUTPUT_ROOT.resolve(),
        vetting_work_directory=WORK_ROOT.resolve(), model_directory=MODEL_ROOT.resolve(),
        dustmaps_directory=DUSTMAP_ROOT.resolve(), workspace_directory=ROOT/'chance_alignment_work',
        workers=4, bootstrap=1000)
    resolved.update(supplied)
    required = ['master_catalogue', 'quaia_catalogue']
    if any(not resolved.get(name) for name in required):
        raise ValueError('Set master_catalogue and quaia_catalogue in the chance-alignment config.')
    path_keys = [name for name in resolved if name.endswith(('_directory', '_catalogue'))]
    for name in path_keys:
        value = Path(resolved[name]).expanduser()
        resolved[name] = (path.resolve().parent/value).resolve() if not value.is_absolute() else value.resolve()
    if resolved['workers'] < 1 or resolved['bootstrap'] < 100:
        raise ValueError('Require workers >= 1 and bootstrap >= 100.')
    workspace = resolved['workspace_directory']
    protected = [ROOT, resolved['input_directory'], resolved['imaging_input_directory'], resolved['products_directory'],
                 resolved['vetting_work_directory'], resolved['model_directory'],
                 resolved['dustmaps_directory']]
    if workspace == ROOT or any(workspace == parent or workspace.is_relative_to(parent)
                                for parent in protected[1:]):
        raise ValueError('Use a distinct chance-alignment workspace, separate from vetting data and products.')
    if any(parent.is_relative_to(workspace) for parent in protected):
        raise ValueError('The chance-alignment workspace must not contain its input directories.')
    return resolved


def source_paths(config):
    products = config['products_directory']
    inputs = config['input_directory']
    imaging_inputs = config['imaging_input_directory']
    return [products/'one_arcsecond_scores.csv', products/'imaging_scores.csv',
        products/'threshold_calibration.json', config['vetting_work_directory']/'one_arcsecond_inputs/one_arcsecond_components.fits',
        imaging_inputs/'imaging_components.fits', imaging_inputs/'imaging_tractor.fits',
        config['vetting_work_directory']/'imaging/labelled_systems_audit.csv',
        config['master_catalogue'], config['quaia_catalogue']]


def file_hash(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def prepare_runtime(config):
    paths = source_paths(config)
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError('Run QSO vetting and supply the control parent/Quaia catalogues first.\n'
                                'Missing inputs:\n'+'\n'.join(missing))
    workspace = config['workspace_directory']
    inputs = workspace/'control_inputs'
    prepared = inputs/'prepared'
    reference = inputs/'reference'
    provenance = workspace/'input_sources.json'
    current = {str(path): file_hash(path) for path in paths}
    if provenance.exists():
        if json.loads(provenance.read_text()) != current:
            raise ValueError('Vetting inputs or control catalogues changed. Use a new workspace_directory.')
    else:
        prepared.mkdir(parents=True, exist_ok=True)
        reference.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, '-B', '-m', 'chance_alignments.prepare_inputs',
            '--products', str(config['products_directory']), '--input-root', str(config['input_directory']),
            '--work-root', str(config['vetting_work_directory']), '--quaia', str(config['quaia_catalogue']),
            '--imaging-input-root', str(config['imaging_input_directory']),
            '--output', str(prepared)]
        subprocess.run(command, cwd=ROOT, check=True)
        modules = ['one_arcsecond.py', 'photometry.py', 'redshift_uncertainty.py', 'redshift_inputs.py', 'score_selection.py']
        version = {name: file_hash(ROOT/name) for name in modules}
        calibration = config['products_directory']/'threshold_calibration.json'
        version['calibration_sha256'] = file_hash(calibration)
        shutil.copy2(calibration, reference/'threshold_calibration.json')
        (reference/'production_version.json').write_text(json.dumps(version, indent=2)+'\n')
        provenance.write_text(json.dumps(current, indent=2)+'\n')
    runtime = dict(data_directory=str(inputs), output_directory=str(workspace/'run'),
        production_repo=str(ROOT), model_directory=str(config['model_directory']),
        dustmaps_directory=str(config['dustmaps_directory']),
        master_catalogue=str(config['master_catalogue']), quaia_catalogue=str(config['quaia_catalogue']),
        workers=config['workers'], bootstrap=config['bootstrap'])
    path = workspace/'runtime_config.json'
    text = json.dumps(runtime, indent=2)+'\n'
    if path.exists() and path.read_text() != text:
        raise ValueError('Run settings changed. Use a new workspace_directory.')
    path.write_text(text)
    return path


def run_full(config, stage='all', resume=False, dry_run=False):
    workspace = config['workspace_directory']
    if dry_run:
        print('Prepare controls from:', config['products_directory'], 'and', config['vetting_work_directory'])
        runtime_path = workspace/'runtime_config.json'
    elif stage in ['all', 'prepare']:
        runtime_path = prepare_runtime(config)
    else:
        runtime_path = workspace/'runtime_config.json'
        if not runtime_path.is_file():
            raise FileNotFoundError('Run --stage prepare first: '+str(runtime_path))
    command = [sys.executable, '-B', '-m', 'chance_alignments.run_controls',
               '--config', str(runtime_path), '--stage', stage, '--resample']
    if resume:
        command.append('--resume')
    if dry_run:
        print('Then:', ' '.join(command))
        print('Stages: prepare -> fetch -> components -> blends -> report')
        return
    subprocess.run(command, cwd=ROOT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, help='Config for a fresh control construction and rescore')
    parser.add_argument('--stage', choices=['all', 'prepare', 'fetch', 'components', 'blends', 'report'], default='all')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--replay-data', type=Path, help='Optional previously supplied control_data directory')
    parser.add_argument('--replay-output', type=Path, default=ROOT/'chance_alignment_work/replay')
    parser.add_argument('--verify-reference', action='store_true')
    args = parser.parse_args()
    if args.replay_data:
        if args.config or args.stage != 'all' or args.dry_run:
            parser.error('Use replay options separately from a full-run config/stage.')
        from chance_alignments.reproduce_results import reproduce
        reproduce(args.replay_data, args.replay_output, verify_reference=args.verify_reference)
        return
    if not args.config:
        parser.error('Supply --config for a fresh run, or --replay-data for an offline replay.')
    run_full(resolve_config(args.config), args.stage, args.resume, args.dry_run)


if __name__ == '__main__':
    main()
