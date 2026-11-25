import numpy as np
import pandas as pd
from dual_agn.settings import *
from qso_pcolor import Photometry, RedshiftMatch
from qso_pcolor.unified import UnifiedPSFModel
from dual_agn.photometry import configure_dustmaps, RESOLVED_POLICY
from dual_agn.imaging import join_system_scores, validate_component_photometry
from dual_agn.redshift_uncertainty import score_components_with_redshift_errors
configure_dustmaps()
components = pd.read_csv(
    IMAGING_ROOT / 'component_photometry.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)
candidate_systems = pd.read_csv(
    IMAGING_ROOT / 'precolour_systems.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)
validate_component_photometry(candidate_systems, components)
models = {name: UnifiedPSFModel.load(path) for name, path in MODEL_BUNDLES.items()}
score_rows = []
for model_name, model in models.items():
    for native_system, group in components.groupby('photometric_system', sort=True):
        group = group.sort_values('candidate_id')
        for start in range(0, len(group), 64):
            batch = group.iloc[start:start + 64]
            bands = tuple((f'{native_system}:{band}' for band in 'grz'))
            photometry = Photometry(
                batch[[f'flux_{b}' for b in 'grz']].to_numpy(float),
                batch[[f'variance_{b}' for b in 'grz']].to_numpy(float),
                bands
            )
            redshift_errors = batch.adopted_redshift_error.fillna(0.0).to_numpy(float)
            scores, decision = score_components_with_redshift_errors(
                model, photometry, redshift_errors,
                ra_deg=batch.component_ra.to_numpy(float),
                dec_deg=batch.component_dec.to_numpy(float),
                z_primary=batch.adopted_redshift.to_numpy(float),
                morphology=batch.tractor_type.to_numpy(str),
                separation_arcsec=batch.pair_separation_arcsec.to_numpy(float),
                fracflux=batch.tractor_fracflux_r.to_numpy(float),
                candidate_id=batch.candidate_id.to_numpy(str),
                primary_id=batch.physical_system_key.to_numpy(str),
                match=RedshiftMatch(half_width_kms=MATCH_HALF_WIDTH_KMS),
                blend_policy=RESOLVED_POLICY,
                ood_flag_sigma=4.0
            )
            probabilities = model.quasar_probability(scores)
            for index, (candidate, score) in enumerate(zip(batch.itertuples(), scores)):
                row = {} if score is None else score.as_row()
                eligible = bool(decision['eligible'][index])
                reason = str(decision['reason'][index])
                row.update(
                    candidate_id=candidate.candidate_id,
                    primary_redshift_error_used=float(redshift_errors[index]),
                    qso_parent_index=candidate.qso_parent_index,
                    component_role=candidate.component_role,
                    model=model_name,
                    eligible=eligible,
                    decision_reason=reason,
                    status=reason or row.get('status', 'not_scored'),
                    p_quasar=float(probabilities[index]) if eligible else np.nan,
                    qso_support_percentile=float(decision['qso_support']['percentile'][index])
                )
                for field in ('log_r_per_unit_z', 'p_sameq'):
                    row[field] = row.get(field, np.nan) if eligible else np.nan
                row.setdefault('quality_flags', '')
                separation = candidate.pair_separation_arcsec
                if np.isfinite(separation) and 0 <= separation < CLOSE_SEPARATION_ARCSEC:
                    flags = [flag for flag in row['quality_flags'].split(',') if flag]
                    row['quality_flags'] = ','.join(dict.fromkeys([*flags, 'blended', 'close_separation']))
                score_rows.append(row)
            print(f'{model_name} {native_system}: {min(start + 64, len(group))}/{len(group)}', flush=True)
component_scores = pd.DataFrame(score_rows)
system_scores = join_system_scores(candidate_systems, component_scores)
component_scores.to_csv(IMAGING_ROOT / 'component_scores_by_model.csv', index=False)
system_scores.to_csv(IMAGING_ROOT / 'system_scores.csv', index=False)
print('Output:', IMAGING_ROOT, '| systems:', len(system_scores))
