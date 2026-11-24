"""Check separately supplied analysis inputs and pinned model assets."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from settings import (DUSTMAP_ROOT, IMAGING_MORPHOLOGY_PATH, INPUT_ROOT,
                      IMAGING_INPUT_ROOT, IMAGING_SEARCH_ROOT,
                      MODEL_BUNDLES, MODEL_REVISION, PACKAGE_ROOT,
                      VISUAL_REJECTIONS_PATH, survey_products)


def sha256(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--external-only',action='store_true',help='Check external inputs before generating the imaging-search tables')
    args = parser.parse_args()
    catalogue_names = [
        'original_parent_712.csv', 'eligibility_audit_712.csv',
        'secondary_measurements_324.csv', 'one_arcsecond_resolution.csv',
        'one_arcsecond_components.fits', 'one_arcsecond_tractor.csv',
        'desi_target_photometry.csv', 'quaia_712_reference.fits',
    ]
    required_paths = [INPUT_ROOT / name for name in catalogue_names]
    required_paths.extend(specification['path'] for specification in survey_products.values())
    required_paths.extend([IMAGING_MORPHOLOGY_PATH, VISUAL_REJECTIONS_PATH])
    generated = []
    if not args.external_only:
        generated = [IMAGING_INPUT_ROOT/name for name in ['imaging_channel.fits',
            'imaging_components.fits','imaging_labelled.fits','imaging_tractor.fits',
            'quaia_imaging_reference.fits']]
        required_paths.extend(generated)
    missing = [str(path) for path in required_paths if not path.is_file()]
    if missing:
        raise SystemExit(
            'Supply the external inputs and run the offline imaging search as described '
            'in data/README.md. Set PAPER1_INPUT_ROOT and the optional path overrides.\n'
            'Missing input files:\n' + '\n'.join(missing))
    if generated and IMAGING_INPUT_ROOT == IMAGING_SEARCH_ROOT/'inputs':
        summary_path = IMAGING_SEARCH_ROOT/'search_summary.json'
        if not summary_path.is_file():
            raise SystemExit('Run the imaging search first; its completed local-data audit is missing.')
        summary = json.loads(summary_path.read_text())
        if summary.get('stage') != 'all' or any(
                summary.get('output_sha256',{}).get('inputs/'+path.name) != sha256(path)
                for path in generated):
            raise SystemExit('Imaging search is incomplete or generated files changed. Rerun the search before vetting.')
    expected_manifests = {
        'magdep': 'e6ff672f2cfae305a5c5e792197f9740a2695cc1dd77290ad1e9376488cb59c6',
        'magindep': '3b92fe80e1b7054da02830a16e013b6f7a06c16a6f2087d0360d668ae8fab6bb',
    }
    for name, path in MODEL_BUNDLES.items():
        if sha256(path / 'manifest.json') != expected_manifests[name]:
            raise SystemExit(f'Unexpected model manifest at {path}')
    source_revision = subprocess.check_output(
        ['git', '-C', str(PACKAGE_ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    if source_revision != MODEL_REVISION:
        raise SystemExit(f'Expected model code at {MODEL_REVISION}, found {source_revision}')
    for hemisphere in ('ngp', 'sgp'):
        path = DUSTMAP_ROOT / f'sfd/SFD_dust_4096_{hemisphere}.fits'
        if not path.is_file():
            raise SystemExit(f'Missing extinction map: {path}')
    print(f'Found {len(required_paths)-len(generated)} external input files and {len(generated)} generated imaging-search tables; '
          'verified pinned model manifests/revision and extinction-map availability.')


if __name__ == '__main__':
    main()
