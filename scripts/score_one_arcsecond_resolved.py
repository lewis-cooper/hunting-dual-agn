from pathlib import Path
import json
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.table import Table
from settings import *
from photometry import read_fits
from qso_pcolor.unified import UnifiedPSFModel
from photometry import configure_dustmaps
from one_arcsecond import score_resolved_candidates, join_resolved_model_scores
configure_dustmaps()
models = {name: UnifiedPSFModel.load(path) for name, path in MODEL_BUNDLES.items()}
systems = pd.read_csv(
    ONE_ROOT / 'scoring_inputs.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)
candidates = pd.read_csv(
    ONE_ROOT / 'original_resolved_measurements.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)
supplemental_mask = systems.scoring_route.eq('resolved_supplemental')
resolved_scores_original = score_resolved_candidates(candidates, models)
ranked = join_resolved_model_scores(candidates, resolved_scores_original)
additional_candidates = systems.loc[supplemental_mask].copy()
additional_candidates['photometric_system'] = np.select(
    [
        additional_candidates.secondary_decals_release.isin([9010, 9012]),
        additional_candidates.secondary_decals_release.eq(9011)
    ],
    ['decals_dr9_south', 'decals_dr9_north'],
    default='unknown'
)
additional_sky = SkyCoord(
    additional_candidates.secondary_component_ra.to_numpy(float) * u.deg,
    additional_candidates.secondary_component_dec.to_numpy(float) * u.deg
)
additional_candidates['secondary_galactic_l_deg'] = additional_sky.galactic.l.deg
additional_candidates['secondary_galactic_b_deg'] = additional_sky.galactic.b.deg
resolved_scores_additional = score_resolved_candidates(additional_candidates, models)
additional_ranked = join_resolved_model_scores(additional_candidates, resolved_scores_additional)
resolved_scores = pd.concat(
    [
        resolved_scores_original.assign(source_subset='original_audit'),
        resolved_scores_additional.assign(source_subset='additional_measured')
    ],
    ignore_index=True
)
resolved_columns = [
    'qso_group_id',
    'primary_qso_catalog',
    'primary_qso_z',
    'photometric_system',
    'log_r_per_unit_z',
    'p_sameq',
    'status',
    'quality_flags',
    'eligible',
    'decision_reason',
    'log_r_per_unit_z_magdep',
    'log_r_per_unit_z_magindep',
    'p_sameq_magdep',
    'p_sameq_magindep',
    'p_quasar_magdep',
    'p_quasar_magindep',
    'eligible_magdep',
    'eligible_magindep',
    'decision_reason_magdep',
    'decision_reason_magindep',
    'status_magdep',
    'status_magindep'
]
resolved_all = pd.concat(
    [
        ranked[resolved_columns].assign(source_subset='original_audit'),
        additional_ranked[resolved_columns].assign(source_subset='additional_measured')
    ],
    ignore_index=True
)
resolved_all['rank_all_resolved'] = resolved_all.log_r_per_unit_z.rank(ascending=False, method='min')
print(
    'Eligible resolved pairs:',
    {name: int(resolved_all[f'eligible_{name}'].sum()) for name in MODEL_BUNDLES}
)
resolved_scores.to_csv(ONE_ROOT / 'resolved_scores_by_model.csv', index=False)
ranked.to_csv(ONE_ROOT / 'ranked_resolved.csv', index=False)
additional_ranked.to_csv(ONE_ROOT / 'ranked_resolved_supplemental.csv', index=False)
resolved_all.to_csv(ONE_ROOT / 'ranked_resolved_extended.csv', index=False)
print('Output:', ONE_ROOT, '| original/supplemental:', len(ranked), len(additional_ranked))
