"""Build blend audit metadata and retain completed scoring batches."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


BLEND_METADATA_COLUMNS = (
    'qso_group_id', 'primary_qso_catalog', 'primary_qso_z', 'redshift_scope',
    'targetid', 'desi_release', 'desi_brickid', 'desi_objid', 'scoring_route',
    'decals_system_assignment', 'both_decals_r_assignments', 'blend_r_snr',
    'fallback_r_band_input', 'fallback_r_snr', 'galactic_abs_b_deg',
    'n_decals_bands', 'n_extra_bands', 'extra_surveys_used',
    'one_local_source_near_both', 'selected_decals_component_photometry_quality',
)
OPTIONAL_DESI_METADATA = {'targetid', 'desi_brickid', 'desi_objid'}


def build_blend_metadata(candidates):
    """Native preparation may omit old DESI labels; scientific inputs remain required."""
    missing = set(BLEND_METADATA_COLUMNS) - OPTIONAL_DESI_METADATA - set(candidates.columns)
    if missing:
        raise ValueError(f'Missing required blend output fields before scoring: {sorted(missing)}')
    columns = list(BLEND_METADATA_COLUMNS)
    columns += [field for field in ('desi_targetids_json', 'desi_target_photometry_status')
                if field in candidates]
    return candidates.reindex(columns=columns).copy()


def scoring_fingerprint(paths, revision):
    """Invalidate checkpoints if any supplied measurement, model or scoring code changes."""
    digest = hashlib.sha256(revision.encode())
    for path in sorted(set(map(Path, paths))):
        digest.update(str(path.resolve()).encode())
        with path.open('rb') as stream:
            digest.update(hashlib.file_digest(stream, 'sha256').digest())
    return digest.hexdigest()


def _json_scalar(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f'Unsupported score value: {type(value).__name__}')


def checkpoint_blend_scores(path, model_name, group_ids, compute):
    """Reuse only a complete batch with the same ordered identities and model."""
    path = Path(path)
    group_ids = [str(value) for value in group_ids]
    if path.exists():
        saved = json.loads(path.read_text())
        if saved['model'] != model_name or saved['group_ids'] != group_ids:
            raise ValueError(f'Blend checkpoint identities do not match: {path}')
        rows = saved['rows']
        if ([str(row['qso_group_id']) for row in rows] != group_ids
                or any(row['model'] != model_name for row in rows)):
            raise ValueError(f'Incomplete or inconsistent blend checkpoint: {path}')
        for row in rows:
            for field in ('bands_used', 'alpha_range', 'notes'):
                if field in row:
                    row[field] = tuple(row[field])
        print(f'Reusing completed {model_name} batch: {len(rows)} systems', flush=True)
        return rows
    rows = compute()
    if ([str(row['qso_group_id']) for row in rows] != group_ids
            or any(row['model'] != model_name for row in rows)):
        raise ValueError('Computed blend batch does not match the requested identities')
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.partial')
    temporary.write_text(json.dumps(dict(model=model_name, group_ids=group_ids, rows=rows),
                                    default=_json_scalar))
    temporary.replace(path)
    return rows
