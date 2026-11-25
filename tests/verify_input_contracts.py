"""Check separately supplied analysis inputs and pinned model assets."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dual_agn.settings import (DUSTMAP_ROOT, IMAGING_MORPHOLOGY_PATH, INPUT_ROOT,
                      IMAGING_INPUT_ROOT, IMAGING_SEARCH_ROOT,
                      ONE_ARCSECOND_PARENT_PATH, ONE_ARCSECOND_SEARCH_ROOT,
                      MODEL_BUNDLES,
                      ONE_INPUT_ROOT, ROOT, survey_products)
from dual_agn.workflow import PIPELINE_CHOICES
from dual_agn.provenance import verify_imported_scorer


def sha256(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def verify_models():
    verify_imported_scorer()
    expected = {
        'magdep': 'e6ff672f2cfae305a5c5e792197f9740a2695cc1dd77290ad1e9376488cb59c6',
        'magindep': '3b92fe80e1b7054da02830a16e013b6f7a06c16a6f2087d0360d668ae8fab6bb',
    }
    for name, path in MODEL_BUNDLES.items():
        if not (path/'manifest.json').is_file() or sha256(path/'manifest.json') != expected[name]:
            raise SystemExit(f'Missing or unexpected model manifest at {path}; follow README.md before acquiring survey data')
    for hemisphere in ('ngp','sgp'):
        path = DUSTMAP_ROOT/f'sfd/SFD_dust_4096_{hemisphere}.fits'
        if not path.is_file():
            raise SystemExit(f'Missing extinction map: {path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--external-only',action='store_true',help='Check external inputs before generating either discovery stage')
    parser.add_argument('--models-only',action='store_true',help='Check model/dust prerequisites without survey inputs or network access')
    parser.add_argument('--measurement-config', type=Path, default=ROOT/'measurement_inputs.local.json')
    parser.add_argument('--pipeline', choices=PIPELINE_CHOICES, default='both')
    args = parser.parse_args()
    if args.models_only:
        verify_models()
        print('Verified pinned model assets/revision and extinction maps; no survey query performed.')
        return
    from dual_agn.measurement_inputs import read_config, catalogue_files, local_path
    config = read_config(args.measurement_config)
    if args.pipeline != 'gaia_neighbour' and not config.get('morphology_images', {}).get('path_template'):
        raise ValueError('Set morphology_images.path_template for the r-band morphology measurement')
    native = [local_path(config, config['quaia_catalogue'])] if config.get('quaia_catalogue') else []
    surveys = [] if args.pipeline == 'imaging' else list(config['surveys'].values())
    for specification in [config['gaia'], config['tractor'], *([config['desi']] if config.get('desi') else []), *surveys]:
        if specification.get('enabled') is False:
            continue
        native.extend(catalogue_files(config, specification))
        if specification.get('coverage'):
            native.append(local_path(config, specification['coverage']))
        elif specification is not config.get('desi') and specification.get('full_catalogue') is not True:
            raise ValueError('Every positional native input needs explicit coverage or full_catalogue=true')
    catalogue_names = [
        'component_quality_audit.csv',
        'quality_resolved_measurements.csv', 'one_arcsecond_resolution.csv',
        'one_arcsecond_components.fits', 'one_arcsecond_tractor.csv',
        'desi_target_photometry.csv', 'quaia_reference.fits',
    ]
    required_paths = list(dict.fromkeys(native))

    generated = []
    if not args.external_only:
        generated = [] if args.pipeline == 'gaia_neighbour' else [IMAGING_INPUT_ROOT/name for name in ['imaging_channel.fits',
            'imaging_components.fits','imaging_labelled.fits','imaging_tractor.fits',
            'quaia_imaging_reference.fits']]
        required_paths.extend(generated)
        required_paths.extend(ONE_INPUT_ROOT / name for name in catalogue_names)
        required_paths.extend(specification['path'] for specification in survey_products.values())
        required_paths.append(ONE_ARCSECOND_PARENT_PATH)
    missing = [str(path) for path in required_paths if not path.is_file()]
    if missing:
        raise SystemExit(
            'Supply the external inputs and run the offline discovery stages as described '
            'in DATA.md. Set PAPER1_INPUT_ROOT and the optional path overrides.\n'
            'Missing input files:\n' + '\n'.join(missing))
    if not args.external_only:
        summary_path = ONE_ARCSECOND_SEARCH_ROOT/'search_summary.json'
        if not summary_path.is_file():
            raise SystemExit('Run the one-arcsecond discovery first; its input/output audit is missing.')
        summary = json.loads(summary_path.read_text())
        if summary.get('output_sha256') != sha256(ONE_ARCSECOND_PARENT_PATH):
            raise SystemExit('Generated one-arcsecond parent changed; rerun discovery before vetting.')
    if generated and IMAGING_INPUT_ROOT == IMAGING_SEARCH_ROOT/'inputs':
        summary_path = IMAGING_SEARCH_ROOT/'search_summary.json'
        if not summary_path.is_file():
            raise SystemExit('Run the imaging search first; its completed local-data audit is missing.')
        summary = json.loads(summary_path.read_text())
        if summary.get('stage') != 'all' or any(
                summary.get('output_sha256',{}).get('inputs/'+path.name) != sha256(path)
                for path in generated):
            raise SystemExit('Imaging search is incomplete or generated files changed. Rerun the search before vetting.')
    verify_models()
    print(f'Found {len(set(native))} configured native input files; '
          + ('external prerequisites checked' if args.external_only else 'generated discovery and component tables checked')
          + '; verified pinned model manifests/revision and extinction maps. Morphology images are checked for selected systems during measurement.')


if __name__ == '__main__':
    main()
