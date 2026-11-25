from pathlib import Path
import json
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from dual_agn.settings import *
from dual_agn.photometry import read_fits
from dual_agn.imaging import imaging_pair_quality, validate_component_photometry
from dual_agn.score_selection import (MODELS, build_secondary_pm_flags, has_two_distinct_gaia_sources,
                             imaging_model_scores, read_either_model_threshold)
log_r_per_unit_z_threshold = read_either_model_threshold(THRESHOLD_CALIBRATION_PATH)
r_perp_min = R_PERP_MIN
imaging_scores = pd.read_csv(
    IMAGING_ROOT / 'system_scores.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)
imaging_components = pd.read_csv(
    IMAGING_ROOT / 'component_photometry.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)
validate_component_photometry(imaging_scores, imaging_components)
component_quality = imaging_pair_quality(imaging_components).reindex(imaging_scores.qso_parent_index)
component_quality.index = imaging_scores.index
quality_columns = component_quality.add_prefix('both_components_')
imaging_scores = pd.concat([
    imaging_scores.drop(columns=quality_columns.columns, errors='ignore'), quality_columns
], axis=1).copy()
imaging_quality_pass = component_quality.all(axis=1)
imaging_scores['joint_quality_pass'] = imaging_quality_pass
imaging_upstream = pd.read_csv(IMAGING_ROOT / 'upstream_cut_flow.csv', float_precision='round_trip')
morphology = pd.read_csv(IMAGING_MORPHOLOGY_PATH, float_precision='round_trip')
imaging_scores = imaging_scores.drop(
    columns=['image_axis_sigma_arcsec', 'image_perp_sigma_arcsec', 'r_perp_sigma_ratio'],
    errors='ignore'
).merge(
    morphology[[
        'qso_parent_index',
        'image_axis_sigma_arcsec',
        'image_perp_sigma_arcsec',
        'r_perp_sigma_ratio',
        'morphology_measurement_source'
    ]],
    on='qso_parent_index',
    how='left',
    validate='one_to_one'
)
imaging_scores = imaging_scores.copy()
two_gaia = has_two_distinct_gaia_sources(imaging_scores.primary_gaia_source_id,
                                        imaging_scores.secondary_gaia_source_id)
pm_flags = build_secondary_pm_flags(imaging_scores.secondary_pm_sig, two_gaia,
                                    PM_SIGNIFICANCE_MAX)
imaging_scores[pm_flags.columns] = pm_flags
imaging_scores['secondary_pm_significance'] = pd.to_numeric(imaging_scores.secondary_pm_sig, errors='coerce')
model_scores = imaging_model_scores(imaging_scores)
imaging_scores['minimum_component_log_r_either'] = model_scores.max(axis=1)
for model in MODELS:
    imaging_scores[f'passes_pair_score_{model}'] = model_scores[model].ge(log_r_per_unit_z_threshold)
imaging_score_pass = imaging_quality_pass & imaging_scores.minimum_component_log_r_either.ge(log_r_per_unit_z_threshold)
imaging_pm_pass = imaging_score_pass & imaging_scores.secondary_pm_pass
pre_morphology_pass = imaging_pm_pass
morphology_available = np.isfinite(imaging_scores.r_perp_sigma_ratio)
imaging_morphology_pass = pre_morphology_pass & morphology_available & imaging_scores.r_perp_sigma_ratio.ge(r_perp_min)
imaging_final_pass = imaging_morphology_pass.copy()
imaging_scores['joint_morphology_pass'] = imaging_morphology_pass
imaging_scores['joint_score_pass'] = imaging_score_pass
imaging_scores['joint_secondary_pm_pass'] = imaging_pm_pass
imaging_scores['joint_pre_morphology_pass'] = pre_morphology_pass
imaging_scores['joint_final_pass'] = imaging_final_pass
imaging_selected = imaging_scores.loc[imaging_final_pass].copy()
print(f'Imaging: score {imaging_score_pass.sum()} -> secondary PM {imaging_pm_pass.sum()} -> morphology {imaging_morphology_pass.sum()}; automatically selected {len(imaging_selected)}')

def imaging_counts(stage, selection):
    selected = imaging_scores.loc[selection]
    return {
        'Selection stage': stage,
        'One-Gaia': int(selected.candidate_tier.eq('one_gaia').sum()),
        'Two-Gaia': int(selected.candidate_tier.eq('two_gaia').sum()),
        'Total': len(selected)
    }
upstream_counts = json.loads((IMAGING_ROOT / 'upstream_counts.json').read_text())
n_initial_psf_pairs = upstream_counts['initial_psf_pairs']
n_validated_pairs = upstream_counts['validated_pairs']
n_primary_gaia = upstream_counts['primary_gaia']
imaging_count_rows = [
    {
        'Selection stage': 'Initial two-component PSF systems',
        'One-Gaia': pd.NA,
        'Two-Gaia': pd.NA,
        'Total': n_initial_psf_pairs
    },
    {
        'Selection stage': 'Two validated components; fixed primary assignment',
        'One-Gaia': pd.NA,
        'Two-Gaia': pd.NA,
        'Total': n_validated_pairs
    },
    {
        'Selection stage': 'Fixed primary has Gaia',
        'One-Gaia': pd.NA,
        'Two-Gaia': pd.NA,
        'Total': n_primary_gaia
    }
]
upstream_labels = {
    'Gaia tier': 'Gaia-tier definition',
    'Primary PM': 'Primary PM significance < 3',
    'Delta(g-r)': '|Delta(g-r)| <= 1',
    'Tier cut': 'One-Gaia pair-DCR angle > 20 degrees',
    'Deduplicated pre-colour': 'Physical-system deduplication'
}
for stage, label in upstream_labels.items():
    counts = imaging_upstream.loc[imaging_upstream.stage.eq(stage)].set_index('tier').systems
    imaging_count_rows.append({
        'Selection stage': label,
        'One-Gaia': int(counts.one_gaia),
        'Two-Gaia': int(counts.two_gaia),
        'Total': int(counts.sum())
    })
cumulative_imaging = pd.Series(True, index=imaging_scores.index)
for column, label in [
    ('reference_snr_pass', 'Both components: reference-band S/N >= 10'),
    ('latitude_pass', 'Both components: |b| >= 25 degrees'),
    ('psf_pass', 'Both components: PSF morphology')
]:
    cumulative_imaging &= imaging_scores[f'both_components_{column}']
    imaging_count_rows.append(imaging_counts(label, cumulative_imaging))
cumulative_imaging &= imaging_scores[[f'both_components_eligible_{model}' for model in MODELS]].eq(True).any(axis=1)
imaging_count_rows.append(imaging_counts('Remaining eligibility requirements met by both components in either model', cumulative_imaging))
cumulative_imaging &= model_scores.notna().any(axis=1)
imaging_count_rows.append(imaging_counts('Both component scores finite in either model', cumulative_imaging))
cumulative_imaging &= imaging_score_pass
imaging_count_rows.append(imaging_counts(f'Both component ln R >= {log_r_per_unit_z_threshold:.2f} in either model', cumulative_imaging))
cumulative_imaging &= imaging_scores.secondary_pm_pass
imaging_count_rows.append(imaging_counts('Two-Gaia secondary PM significance < 3 (missing retained)', cumulative_imaging))
imaging_count_rows.append(imaging_counts('Morphology measured (coverage)', cumulative_imaging & morphology_available))
cumulative_imaging &= imaging_morphology_pass
imaging_count_rows.append(imaging_counts('Measured R_perp >= 1.05', cumulative_imaging))
imaging_cutflow = pd.DataFrame(imaging_count_rows)
assert cumulative_imaging.equals(imaging_final_pass), 'Imaging cutflow and adopted selection differ'
imaging_cutflow[['One-Gaia', 'Two-Gaia', 'Total']] = imaging_cutflow[['One-Gaia', 'Two-Gaia', 'Total']].astype('Int64')
imaging_scores['selection_status'] = np.select(
    [
        ~imaging_quality_pass,
        ~imaging_score_pass,
        imaging_score_pass & ~imaging_pm_pass,
        pre_morphology_pass & ~morphology_available,
        pre_morphology_pass & morphology_available & ~imaging_morphology_pass
    ],
    [
        'component_quality_rejected',
        'score_rejected_or_ineligible',
        'secondary_pm_rejected',
        'morphology_unmeasurable',
        'morphology_rejected'
    ],
    default='selected'
)
imaging_scores['morphology_available'] = morphology_available
imaging_scores.to_csv(IMAGING_SCORES_PATH, index=False)
imaging_cutflow.to_csv(IMAGING_CUTFLOW_PATH, index=False)
