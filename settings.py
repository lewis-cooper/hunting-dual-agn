"""Editable paths and the adopted scientific settings; no original-project imports."""
from pathlib import Path
import os

ROOT = Path(__file__).resolve().parent
INPUT_ROOT = Path(os.environ.get('PAPER1_INPUT_ROOT', ROOT / 'data/inputs')).expanduser()
OUTPUT_ROOT = Path(os.environ.get('PAPER1_OUTPUT_ROOT', ROOT / 'products')).expanduser()
FIGURE_ROOT = Path(os.environ.get('PAPER1_FIGURE_ROOT', ROOT / 'figures')).expanduser()
WORK_ROOT = Path(os.environ.get(
    'PAPER1_WORK_ROOT', OUTPUT_ROOT.parent / f'{OUTPUT_ROOT.name}_work')).expanduser()
FINAL_CANDIDATES_PATH = OUTPUT_ROOT / 'final_candidates.csv'
PRE_LITERATURE_CANDIDATES_PATH = WORK_ROOT / 'joint_candidates_before_literature.csv'
LITERATURE_ROOT = ROOT / 'data/literature'
LITERATURE_WORK_ROOT = WORK_ROOT / 'literature'
ONE_ARCSECOND_SCORES_PATH = OUTPUT_ROOT / 'one_arcsecond_scores.csv'
IMAGING_SCORES_PATH = OUTPUT_ROOT / 'imaging_scores.csv'
ONE_ARCSECOND_CUTFLOW_PATH = OUTPUT_ROOT / 'one_arcsecond_cutflow.csv'
IMAGING_CUTFLOW_PATH = OUTPUT_ROOT / 'imaging_cutflow.csv'
THRESHOLD_CALIBRATION_PATH = OUTPUT_ROOT / 'threshold_calibration.json'
IMAGING_MORPHOLOGY_PATH = Path(os.environ.get(
    'PAPER1_IMAGING_MORPHOLOGY_PATH', INPUT_ROOT / 'imaging_morphology_measurements.csv')).expanduser()
VISUAL_REJECTIONS_PATH = Path(os.environ.get(
    'PAPER1_VISUAL_REJECTIONS_PATH', INPUT_ROOT / 'visual_inspection_rejections.csv')).expanduser()
PACKAGE_ROOT = Path(os.environ.get('QSO_PCOLOR_PACKAGE_ROOT', ROOT / '../qso_p_color')).expanduser()
DUSTMAP_ROOT = Path(os.environ.get('QSO_PCOLOR_DUSTMAPS_DIR', INPUT_ROOT / 'dustmaps')).expanduser()
MODEL_ROOT = Path(os.environ.get('PAPER1_MODEL_ROOT', PACKAGE_ROOT / 'models/multisurvey_psf')).expanduser()
MODEL_BUNDLES = {name: MODEL_ROOT / f'20261005_{name}_qlf' for name in ('magdep', 'magindep')}
MODEL_REVISION = 'dbc392623d27e22744d5224f9f1a3fd52f6230a8'
ONE_ROOT = WORK_ROOT / 'one_arcsecond'
IMAGING_ROOT = WORK_ROOT / 'imaging'
MATCH_HALF_WIDTH_KMS = 2000.0
CLOSE_SEPARATION_ARCSEC = 3.0
PM_SIGNIFICANCE_MAX = 3.0
DELTA_G_MINUS_R_MIN = -1.0
DELTA_G_MINUS_R_MAX = 1.0
DCR_OFFSET_MIN_DEG = 20.0
R_PERP_MIN = 1.05
PAIR_MATCH_TOLERANCE_ARCSEC = 0.1
IDENTIFIER_DTYPES = {name: 'string' for name in (
 'qso_group_id', 'primary_component_gaia_source_id', 'secondary_component_gaia_source_id',
 'source_id', 'original_gaia_source_id', 'primary_gaia_source_id', 'secondary_gaia_source_id',
 'redshift_error_source_id', 'targetid')}

survey_products = {
    'sdss': {
        'path': INPUT_ROOT / 'sdss_measurements.csv',
        'bands': {band: (f'{band}mag', f'{band}mag_err') for band in 'ugriz'},
        'source_flag': 'passes_sdss_imaging_match'
    },
    'allwise': {
        'path': INPUT_ROOT / 'allwise_measurements.csv',
        'bands': {band: (f'allwise_{band}mpro', f'allwise_{band}sigmpro') for band in ('w1', 'w2')},
        'source_flag': 'wise_w1_minus_w2_photometry_quality'
    },
    'ps1': {
        'path': INPUT_ROOT / 'ps1_measurements.csv',
        'bands': {band: (f'ps1_{band}kronmag', f'ps1_{band}kronmagerr') for band in 'grizy'},
        'source_flag': 'ps1_source_matched'
    },
    'skymapper': {
        'path': INPUT_ROOT / 'skymapper_measurements.csv',
        'bands': {band: (f'skymapper_{band}_petro', f'skymapper_e_{band}_petro') for band in 'uvgriz'},
        'source_flag': 'skymapper_source_matched'
    },
    'nsc': {
        'path': INPUT_ROOT / 'nsc_measurements.csv',
        'bands': {band: (f'nsc_{band}mag', f'nsc_{band}err') for band in ('u', 'g', 'r', 'i', 'z', 'y', 'vr')},
        'source_flag': 'nsc_source_matched'
    },
    'vhs': {
        'path': INPUT_ROOT / 'vhs_measurements.csv',
        'bands': {band: (f'vhs_{band}apermag3', f'vhs_{band}apermag3err') for band in ('y', 'j', 'h', 'ks')},
        'source_flag': 'vhs_source_matched'
    }
}
