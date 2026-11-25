"""Verify the scorer that Python actually imported before running the paper analysis."""

import importlib
from pathlib import Path
import subprocess

from dual_agn.settings import MODEL_REVISION, PACKAGE_ROOT


def verify_imported_scorer(package_root=PACKAGE_ROOT, expected_revision=MODEL_REVISION):
    """Require the intended import, the full pinned Git revision and a clean checkout."""
    package_root = Path(package_root).expanduser().resolve()
    module = importlib.import_module('qso_pcolor')
    origin = getattr(module, '__file__', None)
    if origin is None:
        raise RuntimeError('The imported qso_pcolor package has no identifiable source file')
    imported_file = Path(origin).resolve()
    expected_file = package_root/'src/qso_pcolor/__init__.py'
    if imported_file != expected_file.resolve():
        raise RuntimeError(f'Python imported qso_pcolor from {imported_file}, expected {expected_file}. Check the active environment and PYTHONPATH.')
    try:
        checkout = Path(subprocess.check_output(
            ['git', '-C', str(imported_file.parent), 'rev-parse', '--show-toplevel'],
            text=True, stderr=subprocess.PIPE).strip()).resolve()
        revision = subprocess.check_output(
            ['git', '-C', str(checkout), 'rev-parse', 'HEAD'],
            text=True, stderr=subprocess.PIPE).strip()
        changes = subprocess.check_output(
            ['git', '-C', str(checkout), 'status', '--porcelain', '--untracked-files=normal'],
            text=True, stderr=subprocess.PIPE)
    except (OSError, subprocess.CalledProcessError):
        raise RuntimeError('Cannot verify the imported scorer checkout; install the pinned Git checkout described in README.md') from None
    if checkout != package_root:
        raise RuntimeError('The imported scorer is not rooted in the configured qso_pcolor Git checkout')
    if revision != expected_revision:
        raise RuntimeError(f'This analysis requires qso_pcolor {expected_revision}; Python imported revision {revision}. Do not use the later scoring rules for this paper.')
    if changes.strip():
        raise RuntimeError('The imported qso_pcolor checkout has local changes. Use a clean checkout of the pinned revision; the workflow will not alter it automatically.')
    record = dict(package='qso_pcolor', revision=revision, clean_checkout=True,
                  import_relative_path=imported_file.relative_to(checkout).as_posix())
    print(f'Verified imported scorer: {imported_file}; revision {revision}; clean checkout', flush=True)
    return record
