"""Discover and vet candidate pairs from user-supplied local catalogues."""

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
from dual_agn.provenance import verify_imported_scorer
from dual_agn.progress import command_label, report_progress

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
             measurement_config=None, quasar_catalogue_config=None, pipeline='both',
             one_arcsecond_search_config=None):
    selected = selected_pipelines(pipeline)
    if stage in ('all', 'chance-alignments'):
        if chance_config is None:
            raise ValueError('Supply --chance-config when including chance alignments.')
        if pipeline != 'both':
            raise ValueError('The independent chance-alignment experiment compares both pipelines; use --pipeline both, or --stage vetting for a single candidate pipeline.')
    if stage == 'chance-alignments':
        command = [sys.executable, '-u', '-m', 'run_chance_alignments', '--config', str(chance_config.resolve())]
        return [command + (['--resume'] if resume else [])]
    if stage == 'discovery' and pipeline == 'imaging':
        command = [sys.executable, '-u', '-m', 'scripts.imaging.search_candidates',
            '--config', str((imaging_search_config or ROOT/'imaging_search.local.json').resolve())]
        return [command + (['--resume'] if resume else [])]
    result = []
    if stage in ('all', 'vetting'):
        result.append([sys.executable, '-B', 'tests/verify_input_contracts.py', '--models-only'])
    quasar_catalogue_config = (quasar_catalogue_config or ROOT/'quasar_catalogue.local.json').resolve()
    measurement_config = (measurement_config or ROOT/'measurement_inputs.local.json').resolve()
    adapted = (WORK_ROOT/'quasar_catalogue').resolve()
    if one_arcsecond_search_config is None:
        result.append([sys.executable, '-u', '-m', 'scripts.shared.prepare_quasar_catalogue',
            '--config', str(quasar_catalogue_config.resolve()), '--pipeline', pipeline,
            '--imaging-search-config', str((imaging_search_config or ROOT/'imaging_search.local.json').resolve())])
        one_arcsecond_search_config = adapted/'one_arcsecond_search.json'
        imaging_search_config = adapted/'imaging_search.json'
    gaia_search = [sys.executable, '-u', '-m', 'scripts.gaia_neighbour.search_candidates',
        '--config', str((one_arcsecond_search_config or ROOT/'one_arcsecond_search.local.json').resolve())]
    imaging_search = [sys.executable, '-u', '-m', 'scripts.imaging.search_candidates',
        '--config', str((imaging_search_config or ROOT/'imaging_search.local.json').resolve())]
    if resume:
        gaia_search.append('--resume')
        imaging_search.append('--resume')
    if stage == 'discovery':
        return result + ([gaia_search] if pipeline == 'gaia_neighbour' else [imaging_search] if pipeline == 'imaging' else [gaia_search, imaging_search])
    if stage in ('all', 'vetting'):
        if pipeline == 'imaging':
            print('Imaging selection: Gaia-neighbour resolved scoring is used for threshold calibration only.', flush=True)
        contract = [sys.executable, '-B', 'tests/verify_input_contracts.py',
                    '--measurement-config', str(measurement_config), '--pipeline', pipeline]
        result += [[sys.executable, '-B', '-m', 'unittest', 'discover', '-s', 'tests', '-p', 'test_*.py'],
                   [sys.executable, '-B', 'tests/verify_redshift_uncertainty.py'],
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
        result.append([sys.executable, '-B', 'tests/verify_example_run.py', '--pipeline', pipeline])
    if stage in ('all', 'chance-alignments'):
        command = [sys.executable, '-u', '-m', 'run_chance_alignments', '--config', str(chance_config.resolve())]
        if resume:
            command.append('--resume')
        result.append(command)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pipeline', choices=PIPELINE_CHOICES, default='both', help='Candidate pipeline (default: both). Imaging calibrates on Gaia neighbours automatically.')
    parser.add_argument('--stage', choices=['all', 'vetting', 'discovery', 'chance-alignments'], help='Default: vetting, or all when --chance-config is supplied')
    parser.add_argument('--chance-config', type=Path)
    parser.add_argument('--imaging-search-config', type=Path)
    parser.add_argument('--quasar-catalogue-config', type=Path)
    parser.add_argument('--one-arcsecond-search-config', type=Path,
                        help='Supplied quasar–Gaia associations and neighbour counts; skips raw-Gaia catalogue adaptation')
    parser.add_argument('--measurement-config', type=Path)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--check-inputs', action='store_true', help='Check supplied local files without searching or scoring')
    args = parser.parse_args()
    args.stage = args.stage or ('all' if args.chance_config else 'vetting')
    try:
        if args.check_inputs:
            from scripts.shared.check_local_inputs import check_inputs
            check_inputs(args.quasar_catalogue_config or ROOT/'quasar_catalogue.local.json',
                args.measurement_config or ROOT/'measurement_inputs.local.json',
                args.imaging_search_config or ROOT/'imaging_search.local.json', args.pipeline, args.stage,
                one_arcsecond_search_config=args.one_arcsecond_search_config)
            return
        sequence = commands(args.stage, chance_config=args.chance_config, resume=args.resume,
            imaging_search_config=args.imaging_search_config, measurement_config=args.measurement_config,
            quasar_catalogue_config=args.quasar_catalogue_config, pipeline=args.pipeline,
            one_arcsecond_search_config=args.one_arcsecond_search_config)
        if not args.dry_run and args.stage in ('all', 'vetting'):
            inactive = product_names('both') - product_names(args.pipeline)
            if any((OUTPUT_ROOT/name).exists() for name in inactive):
                raise ValueError('This output folder contains products from another pipeline. Set PAPER1_OUTPUT_ROOT, PAPER1_WORK_ROOT and PAPER1_FIGURE_ROOT to a separate run location.')
    except (ValueError, OSError) as error:
        parser.error(str(error))
    if not args.dry_run and args.stage in ('all', 'vetting', 'chance-alignments'):
        verify_imported_scorer()
    for number, command in enumerate(sequence, 1):
        print('Running:', ' '.join(command), flush=True)
        if not args.dry_run:
            with report_progress(f'Analysis {number}/{len(sequence)}: {command_label(command)}'):
                subprocess.run(command, cwd=ROOT, check=True)
    if not args.dry_run and args.stage in ('all', 'vetting'):
        with FINAL_CANDIDATES_PATH.open(newline='') as stream:
            count = sum(1 for _ in csv.DictReader(stream))
        print(f'Pipeline: {args.pipeline}\nFinal candidate count: {count:,}\nCatalogue: {FINAL_CANDIDATES_PATH}', flush=True)


if __name__ == '__main__':
    main()
