"""Shared photometric, proper-motion and DCR selection rules."""

import json

import numpy as np
import pandas as pd


MODELS = ('magdep', 'magindep')


def has_two_distinct_gaia_sources(primary_source_id, secondary_source_id):
    """Compare exact Gaia identifiers, excluding absent or placeholder IDs."""
    primary = primary_source_id.astype('string')
    secondary = secondary_source_id.astype('string')
    primary_valid = primary.str.fullmatch('[1-9][0-9]*', na=False)
    secondary_valid = secondary.str.fullmatch('[1-9][0-9]*', na=False)
    # Use one nullable-boolean backend before combining string comparisons.
    different = primary.ne(secondary).fillna(False).astype('boolean')
    return primary_valid.astype('boolean') & secondary_valid.astype('boolean') & different


def build_secondary_pm_flags(pm_significance, two_gaia, threshold):
    """Reject measured high-PM companions; retain and flag missing measurements."""
    significance = pd.to_numeric(pm_significance, errors='coerce')
    available = np.isfinite(significance)
    applies = two_gaia.fillna(False).astype(bool)
    return pd.DataFrame({
        'has_two_distinct_gaia_sources': applies,
        'secondary_pm_missing': applies & ~available,
        'secondary_pm_pass': ~applies | ~available | significance.lt(threshold),
    }, index=pm_significance.index)


def build_imaging_dcr_flags(pair_dcr_delta_deg, one_gaia, threshold):
    """Apply the DCR angle requirement to one-Gaia imaging pairs only."""
    angle = pd.to_numeric(pair_dcr_delta_deg, errors='coerce')
    applies = one_gaia.fillna(False).astype(bool)
    measured_angle_pass = np.isfinite(angle) & angle.gt(threshold)
    return pd.DataFrame({
        'dcr_cut_applied': applies,
        'dcr_pass': ~applies | measured_angle_pass,
    }, index=pair_dcr_delta_deg.index)


def read_either_model_threshold(path):
    """Load a threshold fitted to the same statistic used for selection."""
    calibration = json.loads(path.read_text())
    if calibration.get('model') != 'either':
        raise ValueError(f'{path} must be recalibrated with scripts.calibrate_resolved_threshold')
    return calibration['adopted_log_r_threshold']


def eligible_model_scores(rows, eligibility_prefix, score_prefix):
    """Return one score per model, with ineligible or non-finite scores masked."""
    scores = {
        model: rows[f'{score_prefix}_{model}'].where(
            rows[f'{eligibility_prefix}_{model}'].eq(True)
        )
        for model in MODELS
    }
    scores = pd.DataFrame(scores, index=rows.index)
    return scores.where(np.isfinite(scores))


def imaging_model_scores(rows):
    """Use the lower component score within each model; require both scores."""
    primary = eligible_model_scores(rows, 'both_components_eligible', 'primary_log_r_per_unit_z')
    secondary = eligible_model_scores(rows, 'both_components_eligible', 'secondary_log_r_per_unit_z')
    return np.minimum(primary, secondary)
