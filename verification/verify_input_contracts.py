"""Check separately supplied analysis inputs and pinned model assets."""

import hashlib
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from settings import (DUSTMAP_ROOT, IMAGING_MORPHOLOGY_PATH, INPUT_ROOT,
                      MODEL_BUNDLES, MODEL_REVISION, PACKAGE_ROOT,
                      VISUAL_REJECTIONS_PATH, survey_products)


def sha256(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    catalogue_names = [
        'original_parent_712.csv', 'eligibility_audit_712.csv',
        'secondary_measurements_324.csv', 'one_arcsecond_resolution.csv',
        'one_arcsecond_components.fits', 'one_arcsecond_tractor.csv',
        'desi_target_photometry.csv', 'quaia_712_reference.fits',
        'imaging_channel.fits', 'imaging_components.fits',
        'imaging_labelled.fits', 'imaging_tractor.fits',
        'quaia_imaging_reference.fits',
    ]
    required_paths = [INPUT_ROOT / name for name in catalogue_names]
    required_paths.extend(specification['path'] for specification in survey_products.values())
    required_paths.extend([IMAGING_MORPHOLOGY_PATH, VISUAL_REJECTIONS_PATH])
    missing = [str(path) for path in required_paths if not path.is_file()]
    if missing:
        raise SystemExit(
            'Input data are supplied separately. Set PAPER1_INPUT_ROOT and the optional '
            'measurement/visual-exclusion path overrides described in data/README.md.\n'
            'Missing input files:\n' + '\n'.join(missing))
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
    print(f'Found {len(required_paths)} separately supplied input files; '
          'verified pinned model manifests/revision and extinction-map availability.')


if __name__ == '__main__':
    main()
