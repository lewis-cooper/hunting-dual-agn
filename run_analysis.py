"""Run either candidate pipeline, or both, followed by optional control calculations."""

import argparse
import csv
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from dual_agn.settings import WORK_ROOT, OUTPUT_ROOT, FINAL_CANDIDATES_PATH
from dual_agn.workflow import PIPELINE_CHOICES, product_names, selected_pipelines

CALIBRATION_STAGES = [
    'gaia_neighbour.build_measurements', 'gaia_neighbour.prepare_inputs',
    'gaia_neighbour.score_resolved', 'gaia_neighbour.calibrate_threshold',
]
GAIA_SELECTION_STAGES = ['gaia_neighbour.score_unresolved', 'gaia_neighbour.select_candidates']
IMAGING_SELECTION_STAGES = [
    'imaging.select_precolour', 'imaging.prepare_components', 'imaging.score_components',
    'imaging.measure_morphology', 'imaging.select_candidates',
]
VETTING_STAGES = [
    'gaia_neighbour.search_candidates', 'imaging.search_candidates',
    *CALIBRATION_STAGES, *GAIA_SELECTION_STAGES, *IMAGING_SELECTION_STAGES,
    'shared.combine_candidates', 'shared.crossmatch_literature',
    'gaia_neighbour.plot_resolved', 'gaia_neighbour.plot_unresolved',
    'imaging.plot_results', 'shared.plot_galactic_latitude', 'shared.plot_separation_redshift',
]


def commands(stage, chance_config=None, resume=False, imaging_search_config=None,
             one_arcsecond_search_config=None, measurement_config=None,
             quasar_catalogue_config=None, wsdb_config=None, data_source='wsdb', pipeline='both'):
    selected = selected_pipelines(pipeline)
    if stage in ('all', 'chance-alignments'):
        if chance_config is None:
            raise ValueError('Supply --chance-config when including chance alignments.')
        if pipeline != 'both':
            raise ValueError('The independent chance-alignment experiment compares both pipelines; use --pipeline both, or --stage vetting for a single candidate pipeline.')
    if data_source not in ('wsdb', 'local'):
        raise ValueError('Choose data_source wsdb or local')
    if data_source == 'local' and wsdb_config is not None:
        raise ValueError('--wsdb-config cannot be used with --data-source local')
    result = []
    if data_source == 'wsdb' and stage != 'chance-alignments':
        if any(value is not None for value in (imaging_search_config, one_arcsecond_search_config,
                                               measurement_config, quasar_catalogue_config)):
            raise ValueError('WSDB is the default. Set paths inside its config, or choose --data-source local; do not combine WSDB with local input CLI options')
        wsdb_config = (wsdb_config or ROOT/'wsdb_inputs.local.json').resolve()
        generated = (WORK_ROOT/'wsdb_inputs').resolve()
        query = [sys.executable, '-u', '-m', 'scripts.shared.acquire_inputs', '--config', str(wsdb_config)]
        if stage in ('one-arcsecond-search', 'imaging-search'):
            query += ['--stage', stage]
        elif pipeline != 'both':
            query += ['--pipeline', pipeline]
        if resume:
            query.append('--resume')
        if stage in ('all', 'vetting'):
            result.append([sys.executable, '-B', 'verification/verify_input_contracts.py', '--models-only'])
        result.append(query)
        imaging_search_config = generated/'imaging_search.json'
        one_arcsecond_search_config = generated/'one_arcsecond_search.json'
        measurement_config = generated/'measurement_inputs.json'
    measurement_config = (measurement_config or ROOT/'measurement_inputs.local.json').resolve()
    if quasar_catalogue_config is not None:
        adapted = (WORK_ROOT/'quasar_catalogue').resolve()
        result.append([sys.executable, '-u', '-m', 'scripts.shared.prepare_quasar_catalogue',
            '--config', str(quasar_catalogue_config.resolve()), '--pipeline', pipeline,
            '--imaging-search-config', str((imaging_search_config or ROOT/'imaging_search.local.json').resolve())])
        one_arcsecond_search_config = adapted/'one_arcsecond_search.json'
        imaging_search_config = adapted/'imaging_search.json'
    gaia_search = [sys.executable, '-u', '-m', 'scripts.gaia_neighbour.search_candidates',
        '--config', str((one_arcsecond_search_config or ROOT/'one_arcsecond_search.local.json').resolve())]
    imaging_search = [sys.executable, '-u', '-m', 'scripts.imaging.search_candidates',
        '--config', str((imaging_search_config or ROOT/'imaging_search.local.json').resolve())]
    if resume or data_source == 'wsdb':
        gaia_search.append('--resume')
        imaging_search.append('--resume')
    if stage == 'one-arcsecond-search':
        return result + [gaia_search]
    if stage == 'imaging-search':
        return result + [imaging_search]
    if stage in ('all', 'vetting'):
        if pipeline == 'imaging':
            print('Imaging selection: Gaia-neighbour resolved scoring is used for threshold calibration only.', flush=True)
        contract = [sys.executable, '-B', 'verification/verify_input_contracts.py',
                    '--measurement-config', str(measurement_config), '--pipeline', pipeline]
        result += [[sys.executable, '-B', '-m', 'unittest', 'discover', '-s', 'verification', '-p', 'test_*.py'],
                   [sys.executable, '-B', 'verification/verify_redshift_uncertainty.py'],
                   contract + ['--external-only'], gaia_search]
        if 'imaging' in selected:
            result.append(imaging_search)
        stages = list(CALIBRATION_STAGES)
        if 'gaia_neighbour' in selected:
            stages += GAIA_SELECTION_STAGES
        if 'imaging' in selected:
            stages += IMAGING_SELECTION_STAGES
        stages += ['shared.combine_candidates', 'shared.crossmatch_literature']
        if 'gaia_neighbour' in selected:
            stages += ['gaia_neighbour.plot_resolved', 'gaia_neighbour.plot_unresolved']
        if 'imaging' in selected:
            stages += ['imaging.plot_results']
        stages += ['shared.plot_galactic_latitude', 'shared.plot_separation_redshift']
        for name in stages:
            command = [sys.executable, '-u', '-m', 'scripts.'+name]
            if name in ('gaia_neighbour.build_measurements', 'imaging.measure_morphology'):
                command += ['--config', str(measurement_config)]
            if name == 'gaia_neighbour.build_measurements' and pipeline == 'imaging':
                command.append('--calibration-only')
            if name in ('shared.combine_candidates', 'shared.crossmatch_literature'):
                command += ['--pipeline', pipeline]
            result.append(command)
            if name == 'gaia_neighbour.build_measurements':
                result.append(contract)
        result.append([sys.executable, '-B', 'verification/verify_example_run.py', '--pipeline', pipeline])
    if stage in ('all', 'chance-alignments'):
        command = [sys.executable, '-u', '-m', 'run_chance_alignments', '--config', str(chance_config.resolve())]
        if resume:
            command.append('--resume')
        result.append(command)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pipeline', choices=PIPELINE_CHOICES, default='both', help='Candidate pipeline (default: both). Imaging calibrates on Gaia neighbours automatically.')
    parser.add_argument('--stage', choices=['all', 'vetting', 'one-arcsecond-search', 'imaging-search', 'chance-alignments'], help='Default: vetting, or all when --chance-config is supplied')
    parser.add_argument('--chance-config', type=Path)
    parser.add_argument('--imaging-search-config', type=Path)
    parser.add_argument('--one-arcsecond-search-config', type=Path)
    parser.add_argument('--quasar-catalogue-config', type=Path)
    parser.add_argument('--measurement-config', type=Path)
    parser.add_argument('--data-source', choices=['wsdb', 'local'], default='wsdb')
    parser.add_argument('--wsdb-config', type=Path)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    args.stage = args.stage or ('all' if args.chance_config else 'vetting')
    try:
        sequence = commands(args.stage, args.chance_config, args.resume, args.imaging_search_config,
            args.one_arcsecond_search_config, args.measurement_config, args.quasar_catalogue_config,
            args.wsdb_config, args.data_source, args.pipeline)
        if not args.dry_run and args.stage in ('all', 'vetting'):
            inactive = product_names('both') - product_names(args.pipeline)
            if any((OUTPUT_ROOT/name).exists() for name in inactive):
                raise ValueError('This output folder contains products from another pipeline. Set PAPER1_OUTPUT_ROOT, PAPER1_WORK_ROOT and PAPER1_FIGURE_ROOT to a separate run location.')
    except ValueError as error:
        parser.error(str(error))
    for command in sequence:
        print('Running:', ' '.join(command), flush=True)
        if not args.dry_run:
            subprocess.run(command, cwd=ROOT, check=True)
    if not args.dry_run and args.stage in ('all', 'vetting'):
        with FINAL_CANDIDATES_PATH.open(newline='') as stream:
            count = sum(1 for _ in csv.DictReader(stream))
        print(f'Pipeline: {args.pipeline}\nFinal candidate count: {count:,}\nCatalogue: {FINAL_CANDIDATES_PATH}', flush=True)


if __name__ == '__main__':
    main()
