from pathlib import Path
import json
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from settings import *
from photometry import read_fits
from score_selection import build_imaging_dcr_flags

IMAGING_ROOT.mkdir(parents=True, exist_ok=True)
channel_systems = read_fits(INPUT_ROOT / 'imaging_channel.fits')
imaging_candidates_start = channel_systems.loc[channel_systems.type_pair.eq('PSF-PSF')].copy()
component_rows = read_fits(INPUT_ROOT / 'imaging_components.fits')
labelled_systems = read_fits(INPUT_ROOT / 'imaging_labelled.fits')
# Source cell c41e267e
component_count_by_system = component_rows.groupby('qso_parent_index').size()
validated_component_ids = set(component_count_by_system.loc[component_count_by_system.eq(2)].index)
passes_component_validation = imaging_candidates_start['qso_parent_index'].isin(validated_component_ids)
validated_channel1a_systems = imaging_candidates_start.loc[passes_component_validation].copy()
removed_during_component_validation = imaging_candidates_start.loc[~passes_component_validation].copy()
labelled_ids = set(labelled_systems['qso_parent_index'])
initial_validation_flow = pd.DataFrame({
    'Stage': ['Channel 1a PSF–PSF', 'Two validated components and relabelled'],
    'Systems remaining': [len(imaging_candidates_start), len(validated_channel1a_systems)],
    'Removed at stage': [0, len(removed_during_component_validation)]
})


# Source cell be950603
primary_has_gaia = labelled_systems['primary_has_gaia'].fillna(False).astype(bool)
secondary_has_distinct_gaia = labelled_systems['secondary_has_distinct_gaia'].fillna(False).astype(bool)
two_gaia_tier = primary_has_gaia & secondary_has_distinct_gaia
one_gaia_tier = primary_has_gaia & ~secondary_has_distinct_gaia
no_primary_gaia = ~primary_has_gaia
tier_counts = pd.DataFrame({
    'Group': [
        'Validated/relabelled',
        'Primary has Gaia',
        'One-Gaia tier',
        'Two-Gaia tier',
        'No Gaia match for QSO-colour primary'
    ],
    'Systems': [
        len(labelled_systems),
        int(primary_has_gaia.sum()),
        int(one_gaia_tier.sum()),
        int(two_gaia_tier.sum()),
        int(no_primary_gaia.sum())
    ]
})


# Source cell 04bd346c
primary_pm_significance = pd.to_numeric(labelled_systems['primary_pm_sig'], errors='coerce')
delta_g_minus_r = pd.to_numeric(labelled_systems['delta_g_minus_r'], errors='coerce')
primary_pm_available = np.isfinite(primary_pm_significance)
primary_pm_pass = primary_pm_available & (primary_pm_significance < PM_SIGNIFICANCE_MAX)
one_gaia_after_primary_pm = one_gaia_tier & primary_pm_pass
two_gaia_after_primary_pm = two_gaia_tier & primary_pm_pass


# Source cell c724a0e6
delta_g_minus_r_pass = np.isfinite(delta_g_minus_r) & (delta_g_minus_r >= DELTA_G_MINUS_R_MIN) & (delta_g_minus_r <= DELTA_G_MINUS_R_MAX)
one_gaia_after_delta_gr = one_gaia_after_primary_pm & delta_g_minus_r_pass
two_gaia_after_delta_gr = two_gaia_after_primary_pm & delta_g_minus_r_pass


# Source cell c30fbe05
dcr_flags = build_imaging_dcr_flags(labelled_systems.pair_dcr_delta_deg,
                                   one_gaia_tier, DCR_OFFSET_MIN_DEG)
labelled_systems[dcr_flags.columns] = dcr_flags
one_gaia_before_deduplication = one_gaia_after_delta_gr & dcr_flags.dcr_pass
two_gaia_before_deduplication = two_gaia_after_delta_gr.copy()
# Secondary PM is applied after scoring, with missing values retained and flagged.


# Source cell 3e95a868
def coordinate_pair_key(row):
    """Build the unordered rounded-coordinate key used by the pipeline."""
    component_coordinates = sorted([
        f'{row.primary_ra:.6f},{row.primary_dec:.6f}',
        f'{row.secondary_ra:.6f},{row.secondary_dec:.6f}'
    ])
    return 'coord:' + '|'.join(component_coordinates)

def physical_system_key(row, candidate_tier):
    if candidate_tier == 'two_gaia':
        source_ids = sorted([int(row.primary_gaia_source_id), int(row.secondary_gaia_source_id)])
        if source_ids[0] > 0 and source_ids[0] != source_ids[1]:
            return f'gaia:{source_ids[0]}|{source_ids[1]}'
    return coordinate_pair_key(row)

def deduplicate_candidates(selected_rows, candidate_tier):
    """Keep the preferred parent row for each physical binary system."""
    candidates = selected_rows.copy()
    candidates['candidate_tier'] = candidate_tier
    candidates['physical_system_key'] = [physical_system_key(row, candidate_tier) for row in candidates.itertuples(index=False)]
    primary_ra = pd.to_numeric(candidates['primary_ra'], errors='coerce').to_numpy()
    primary_dec = pd.to_numeric(candidates['primary_dec'], errors='coerce').to_numpy()
    quasar_ra = pd.to_numeric(candidates['q_ra'], errors='coerce').to_numpy()
    quasar_dec = pd.to_numeric(candidates['q_dec'], errors='coerce').to_numpy()
    delta_ra = (primary_ra - quasar_ra + 180.0) % 360.0 - 180.0
    mean_dec_rad = np.deg2rad(0.5 * (primary_dec + quasar_dec))
    candidates['primary_qso_parent_sep_arcsec'] = 3600.0 * np.hypot(delta_ra * np.cos(mean_dec_rad), primary_dec - quasar_dec)
    candidates['abs_delta_g_minus_r'] = np.abs(pd.to_numeric(candidates['delta_g_minus_r'], errors='coerce'))
    preference_columns = [
        'physical_system_key',
        'primary_qso_parent_sep_arcsec',
        'abs_delta_g_minus_r',
        'qso_color_primary_distance_sigma',
        'qso_parent_index'
    ]
    candidates = candidates.sort_values(preference_columns, kind='mergesort')
    return candidates.drop_duplicates('physical_system_key', keep='first').copy()
one_gaia_candidates = deduplicate_candidates(labelled_systems.loc[one_gaia_before_deduplication], 'one_gaia')
two_gaia_candidates = deduplicate_candidates(labelled_systems.loc[two_gaia_before_deduplication], 'two_gaia')
candidate_systems = pd.concat([one_gaia_candidates, two_gaia_candidates], ignore_index=True)
upstream_cut_flow = pd.DataFrame([
    {'tier': tier, 'stage': stage, 'systems': int(count)} for tier,
    stages in {
        'one_gaia': [
            ('Gaia tier', one_gaia_tier.sum()),
            ('Primary PM', one_gaia_after_primary_pm.sum()),
            ('Delta(g-r)', one_gaia_after_delta_gr.sum()),
            ('Tier cut', one_gaia_before_deduplication.sum()),
            ('Deduplicated pre-colour', len(one_gaia_candidates))
        ],
        'two_gaia': [
            ('Gaia tier', two_gaia_tier.sum()),
            ('Primary PM', two_gaia_after_primary_pm.sum()),
            ('Delta(g-r)', two_gaia_after_delta_gr.sum()),
            ('Tier cut', two_gaia_before_deduplication.sum()),
            ('Deduplicated pre-colour', len(two_gaia_candidates))
        ]
    }.items() for stage,
    count in stages
])

candidate_systems.to_csv(IMAGING_ROOT / 'precolour_systems.csv', index=False)
labelled_systems.to_csv(IMAGING_ROOT / 'labelled_systems_audit.csv', index=False)
upstream_cut_flow.to_csv(IMAGING_ROOT / 'upstream_cut_flow.csv', index=False)
(IMAGING_ROOT / 'upstream_counts.json').write_text(json.dumps({
    'initial_psf_pairs': len(imaging_candidates_start), 'validated_pairs': len(labelled_systems),
    'primary_gaia': int(primary_has_gaia.sum()), 'precolour_systems': len(candidate_systems)}, indent=2) + '\n')
print('Input:', INPUT_ROOT, '| Output:', IMAGING_ROOT, '| pre-colour:', len(candidate_systems))
