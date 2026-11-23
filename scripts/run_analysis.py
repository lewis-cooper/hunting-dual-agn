"""Run QSO vetting and, optionally, the unchanged chance-alignment calculations."""

import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
VETTING_STAGES = [
    'prepare_one_arcsecond_inputs', 'score_one_arcsecond_resolved',
    'calibrate_resolved_threshold', 'score_one_arcsecond_unresolved',
    'select_one_arcsecond_candidates', 'select_imaging_precolour',
    'prepare_imaging_components', 'score_imaging_components', 'select_imaging_candidates',
    'combine_candidate_catalogues', 'prepare_photometry_inventory',
    'plot_one_arcsecond_results', 'plot_one_arcsecond_unresolved',
    'plot_imaging_results', 'plot_joint_results', 'plot_candidate_separation_redshift',
]


def commands(stage, chance_config=None, resume=False):
    result = []
    if stage in ['all', 'vetting']:
        result += [[sys.executable, '-B', '-m', 'unittest', 'discover', '-s', 'verification', '-p', 'test_*.py'],
                   [sys.executable, '-B', 'verification/verify_redshift_uncertainty.py'],
                   [sys.executable, '-B', 'verification/verify_input_contracts.py']]
        result += [[sys.executable, '-u', '-m', 'scripts.'+name] for name in VETTING_STAGES]
        result.append([sys.executable, '-B', 'verification/verify_example_run.py'])
    if stage in ['all', 'chance-alignments']:
        if chance_config is None:
            raise ValueError('Supply --chance-config when including chance alignments.')
        command = [sys.executable, '-u', '-m', 'scripts.run_chance_alignments',
                   '--config', str(chance_config.resolve())]
        if resume:
            command.append('--resume')
        result.append(command)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['all', 'vetting', 'chance-alignments'], default='all')
    parser.add_argument('--chance-config', type=Path)
    parser.add_argument('--resume', action='store_true', help='Reuse completed chance-alignment checkpoints')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    try:
        sequence = commands(args.stage, args.chance_config, args.resume)
    except ValueError as error:
        parser.error(str(error))
    for command in sequence:
        print('Running:', ' '.join(command), flush=True)
        if not args.dry_run:
            subprocess.run(command, cwd=ROOT, check=True)


if __name__ == '__main__':
    main()
