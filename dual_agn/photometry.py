from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from qso_pcolor import BlendPolicy
from dual_agn.settings import CLOSE_SEPARATION_ARCSEC

class ResolvedMeasurementPolicy(BlendPolicy):
    """Flag close fitted pairs, ignore fracflux, refuse invalid separations."""
    def violations(self, separation_arcsec, fracflux, n):
        separation = np.broadcast_to(np.asarray(separation_arcsec, float), (n,))
        return ~np.isfinite(separation) | (separation < 0)

    def describe(self):
        return {**super().describe(), 'separation_action': 'flag', 'fracflux_action': 'ignore'}

RESOLVED_POLICY = ResolvedMeasurementPolicy(
    min_separation_arcsec=CLOSE_SEPARATION_ARCSEC, max_fracflux=0.2, action='exclude')

def read_fits(path, columns=None):
    table = Table.read(path, memmap=True)
    if columns is not None:
        table = table[columns]
    frame = table.to_pandas()
    for column in frame.select_dtypes(include='object'):
        frame[column] = frame[column].map(lambda value: value.decode().strip() if isinstance(value, bytes) else value)
    return frame


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def configure_dustmaps():
    """Keep the per-run dustmap configuration with intermediate working files."""
    from dustmaps.config import config
    from dual_agn.settings import DUSTMAP_ROOT, WORK_ROOT
    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    configuration_file = WORK_ROOT / 'dustmaps_configuration.json'
    configuration_file.write_text(json.dumps({'data_dir': str(DUSTMAP_ROOT)}) + '\n')
    config.fname = str(configuration_file)
    config.load()


def model_bundle_records():
    from dual_agn.settings import MODEL_BUNDLES
    return {name: {'bundle': str(path), 'manifest_sha256': file_sha256(path / 'manifest.json')}
            for name, path in MODEL_BUNDLES.items()}

def saved_true(values):
    return values.astype('string').str.strip().str.lower().isin(('true', 't', '1', '1.0', 'yes'))


