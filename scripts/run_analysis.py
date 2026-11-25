"""Run QSO vetting and, optionally, the unchanged chance-alignment calculations."""

import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from settings import WORK_ROOT, FINAL_CANDIDATES_PATH
VETTING_STAGES = [
    'search_one_arcsecond_candidates', 'search_imaging_candidates',
    'build_one_arcsecond_measurements', 'prepare_one_arcsecond_inputs', 'score_one_arcsecond_resolved',
    'calibrate_resolved_threshold', 'score_one_arcsecond_unresolved',
    'select_one_arcsecond_candidates', 'select_imaging_precolour',
    'prepare_imaging_components', 'score_imaging_components', 'measure_imaging_morphology', 'select_imaging_candidates',
    'combine_candidate_catalogues', 'crossmatch_literature', 'prepare_photometry_inventory',
    'plot_one_arcsecond_results', 'plot_one_arcsecond_unresolved',
    'plot_imaging_results', 'plot_joint_results', 'plot_candidate_separation_redshift',
]


def commands(stage, chance_config=None, resume=False, imaging_search_config=None, one_arcsecond_search_config=None, measurement_config=None, quasar_catalogue_config=None, wsdb_config=None, data_source='wsdb'):
    result = []
    if data_source not in ('wsdb', 'local'):
        raise ValueError('Choose data_source wsdb or local')
    if data_source == 'local' and wsdb_config is not None:
        raise ValueError('--wsdb-config cannot be used with --data-source local')
    if data_source == 'wsdb' and stage != 'chance-alignments':
        if any(value is not None for value in (imaging_search_config, one_arcsecond_search_config, measurement_config, quasar_catalogue_config)):
            raise ValueError('WSDB is the default. Set paths inside its config, or choose --data-source local; do not combine WSDB with local input CLI options')
        wsdb_config = (wsdb_config or ROOT/'wsdb_inputs.local.json').resolve()
        generated = (WORK_ROOT/'wsdb_inputs').resolve()
        query_command = [sys.executable, '-u', '-m', 'scripts.query_wsdb_inputs', '--config', str(wsdb_config.resolve())]
        if stage in ('one-arcsecond-search', 'imaging-search'):
            query_command += ['--stage', stage]
        if resume:
            query_command.append('--resume')
        if stage in ('all','vetting'):
            result.append([sys.executable,'-B','verification/verify_input_contracts.py','--models-only'])
        result.append(query_command)
        imaging_search_config = generated/'imaging_search.json'
        one_arcsecond_search_config = generated/'one_arcsecond_search.json'
        measurement_config = generated/'measurement_inputs.json'
    measurement_config = (measurement_config or ROOT/'measurement_inputs.local.json').resolve()
    if quasar_catalogue_config is not None:
        adapted_root = (WORK_ROOT/'quasar_catalogue').resolve()
        result.append([sys.executable, '-u', '-m', 'scripts.prepare_quasar_catalogue',
            '--config', str(quasar_catalogue_config.resolve()), '--imaging-search-config',
            str((imaging_search_config or ROOT/'imaging_search.local.json').resolve())])
        one_arcsecond_search_config = adapted_root/'one_arcsecond_search.json'
        imaging_search_config = adapted_root/'imaging_search.json'
    search_command = [sys.executable, '-u', '-m', 'scripts.search_imaging_candidates',
        '--config', str((imaging_search_config or ROOT/'imaging_search.local.json').resolve())]
    one_search_command = [sys.executable, '-u', '-m', 'scripts.search_one_arcsecond_candidates',
        '--config', str((one_arcsecond_search_config or ROOT/'one_arcsecond_search.local.json').resolve())]
    if resume:
        search_command.append('--resume')
        one_search_command.append('--resume')
    if stage == 'one-arcsecond-search':
        return result + [one_search_command]
    if stage == 'imaging-search':
        return result + [search_command]
    if stage in ['all', 'vetting']:
        result += [[sys.executable, '-B', '-m', 'unittest', 'discover', '-s', 'verification', '-p', 'test_*.py'],
                   [sys.executable, '-B', 'verification/verify_redshift_uncertainty.py'],
                   [sys.executable, '-B', 'verification/verify_input_contracts.py', '--external-only', '--measurement-config', str(measurement_config)]]
        result += [one_search_command, search_command]
        for name in VETTING_STAGES[2:]:
            command = [sys.executable, '-u', '-m', 'scripts.'+name]
            if name in ('build_one_arcsecond_measurements', 'measure_imaging_morphology'):
                command += ['--config', str(measurement_config)]
            result.append(command)
            if name == 'build_one_arcsecond_measurements':
                result.append([sys.executable, '-B', 'verification/verify_input_contracts.py', '--measurement-config', str(measurement_config)])
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
    parser.add_argument('--stage', choices=['all', 'vetting', 'one-arcsecond-search', 'imaging-search', 'chance-alignments'], default='all')
    parser.add_argument('--chance-config', type=Path)
    parser.add_argument('--imaging-search-config', type=Path, help='Local imaging-search catalogue configuration')
    parser.add_argument('--one-arcsecond-search-config', type=Path, help='Local Gaia one-arcsecond discovery configuration')
    parser.add_argument('--quasar-catalogue-config', type=Path, help='Optional column mapping and Gaia data for another quasar catalogue')
    parser.add_argument('--measurement-config', type=Path, help='Native photometry and morphology-image configuration')
    parser.add_argument('--data-source', choices=['wsdb', 'local'], default='wsdb', help='Candidate data source (default: wsdb). Choose local to run without database access.')
    parser.add_argument('--wsdb-config', type=Path, help='WSDB export configuration (default: wsdb_inputs.local.json)')
    parser.add_argument('--resume', action='store_true', help='Reuse verified discovery and chance-alignment checkpoints')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    try:
        sequence = commands(args.stage, args.chance_config, args.resume, args.imaging_search_config, args.one_arcsecond_search_config, args.measurement_config, args.quasar_catalogue_config, args.wsdb_config, args.data_source)
    except ValueError as error:
        parser.error(str(error))
    for command in sequence:
        print('Running:', ' '.join(command), flush=True)
        if not args.dry_run:
            subprocess.run(command, cwd=ROOT, check=True)
    if not args.dry_run and args.stage in ('all', 'vetting'):
        import csv
        with FINAL_CANDIDATES_PATH.open(newline='') as stream:
            count = sum(1 for _ in csv.DictReader(stream))
        print(f'Final candidate count: {count:,}\nCatalogue: {FINAL_CANDIDATES_PATH}', flush=True)


if __name__ == '__main__':
    main()
