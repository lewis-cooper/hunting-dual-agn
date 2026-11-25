from dual_agn.redshift_inputs import photometric_redshifts, redshift_errors
import json
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from dual_agn.settings import *
from dataclasses import asdict
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from multiprocessing import get_context
import os
from qso_pcolor import Photometry
from qso_pcolor.blend import BlendModel
from dual_agn.photometry import configure_dustmaps, saved_true
from dual_agn.redshift_uncertainty import blend_score_with_redshift_error
from dual_agn.blend_priors import cache_population_priors


def load_worker_dustmaps():
    from dustmaps.config import config

    config.fname = str(WORK_ROOT / 'dustmaps_configuration.json')
    config.load()


@lru_cache(maxsize=2)
def cached_blend_model(model_path):
    return cache_population_priors(BlendModel.load(model_path))


def score_batch(blend_model, model_name, photometry, coordinates, group_ids):
    rows = []
    for index, (ra, dec, redshift, error) in enumerate(coordinates):
        score, audit, grid = blend_score_with_redshift_error(
            blend_model, photometry.subset([index]), ra, dec, redshift, error)
        rows.append({'qso_group_id': group_ids[index], 'model': model_name,
                     'band_set': 'multisurvey', **asdict(score), **audit,
                     'redshift_grid_json': json.dumps(grid)})
    return rows


def score_worker_batch(batch):
    model_path, *arguments = batch
    return score_batch(cached_blend_model(model_path), *arguments)


def main():
    configure_dustmaps()
    systems = pd.read_csv(
        ONE_ROOT / 'scoring_inputs.csv',
        dtype=IDENTIFIER_DTYPES,
        float_precision='round_trip'
    )
    blend_input_mask = ~systems.decals_resolution_class.eq('decals_resolved')
    measured_decals_r = {}
    for role in ('primary', 'secondary'):
        flux_r = systems[f'{role}_decals_flux_r']
        ivar_r = systems[f'{role}_decals_flux_ivar_r']
        measured_decals_r[role] = systems[f'{role}_decals_release'].isin((9010, 9011, 9012)) & np.isfinite(flux_r) & np.isfinite(ivar_r) & ivar_r.gt(0)
    model_bundles = {name: {'bundle': path} for name, path in MODEL_BUNDLES.items()}
    blend_candidates = systems.loc[blend_input_mask].copy().reset_index(drop=True)
    component_primary_r = measured_decals_r['primary'].loc[blend_input_mask].to_numpy()
    component_secondary_r = measured_decals_r['secondary'].loc[blend_input_mask].to_numpy()
    blend_candidates['decals_system_assignment'] = np.select(
        [
            component_primary_r,
            component_secondary_r,
            blend_candidates.desi_target_r_available.to_numpy()
        ],
        ['primary', 'secondary', 'desi_target'],
        default='none'
    )
    blend_candidates['both_decals_r_assignments'] = component_primary_r & component_secondary_r
    selected_primary = blend_candidates.decals_system_assignment.eq('primary')
    selected_secondary = blend_candidates.decals_system_assignment.eq('secondary')
    selected_desi = blend_candidates.decals_system_assignment.eq('desi_target')
    blend_candidates['selected_decals_component_photometry_quality'] = np.select(
        [selected_primary, selected_secondary],
        [
            blend_candidates.primary_decals_component_photometry_quality,
            blend_candidates.secondary_decals_component_photometry_quality
        ],
        default=np.nan
    )
    blend_candidates['has_spectroscopic_z'] = ~photometric_redshifts(blend_candidates, 'primary_qso_catalog', 'primary_qso_redshift_kind')
    blend_candidates['redshift_scope'] = np.where(
        blend_candidates.has_spectroscopic_z,
        'spectroscopic',
        'Quaia_Gaussian_redshift_average_exploratory'
    )
    for band in 'grz':
        blend_candidates[f'system_flux_{band}'] = np.select(
            [selected_primary, selected_secondary, selected_desi],
            [
                blend_candidates[f'primary_decals_flux_{band}'],
                blend_candidates[f'secondary_decals_flux_{band}'],
                blend_candidates[f'desi_flux_{band}']
            ],
            default=np.nan
        )
        blend_candidates[f'system_flux_ivar_{band}'] = np.select(
            [selected_primary, selected_secondary, selected_desi],
            [
                blend_candidates[f'primary_decals_flux_ivar_{band}'],
                blend_candidates[f'secondary_decals_flux_ivar_{band}'],
                blend_candidates[f'desi_flux_ivar_{band}']
            ],
            default=np.nan
        )
    blend_candidates['blend_r_snr'] = blend_candidates.system_flux_r * np.sqrt(blend_candidates.system_flux_ivar_r)
    selected_release = np.select(
        [selected_primary, selected_secondary, selected_desi],
        [
            blend_candidates.primary_decals_release,
            blend_candidates.secondary_decals_release,
            blend_candidates.desi_release
        ],
        default=np.nan
    )
    blend_candidates['photometric_system'] = np.select(
        [np.isin(selected_release, [9010, 9012]), selected_release == 9011],
        ['decals_dr9_south', 'decals_dr9_north'],
        default='none'
    )
    blend_models = {name: cached_blend_model(record['bundle']) for name, record in model_bundles.items()}
    band_labels = blend_models['magdep'].bands
    system_count = len(blend_candidates)
    decals_flux = np.full((system_count, len(band_labels)), np.nan)
    decals_variance = np.full_like(decals_flux, np.inf)
    for hemisphere in ('south', 'north'):
        matching_hemisphere = blend_candidates.photometric_system.eq(f'decals_dr9_{hemisphere}').to_numpy()
        for band in 'grz':
            label = f'decals_dr9_{hemisphere}:{band}'
            band_index = band_labels.index(label)
            measured_flux = pd.to_numeric(blend_candidates[f'system_flux_{band}'], errors='coerce').to_numpy(float)
            inverse_variance = pd.to_numeric(blend_candidates[f'system_flux_ivar_{band}'], errors='coerce').to_numpy(float)
            measured = matching_hemisphere & np.isfinite(measured_flux) & np.isfinite(inverse_variance) & (inverse_variance > 0)
            decals_flux[measured, band_index] = measured_flux[measured]
            decals_variance[measured, band_index] = 1 / inverse_variance[measured]
    multisurvey_flux = decals_flux.copy()
    multisurvey_variance = decals_variance.copy()
    blend_total_surveys = ('ps1', 'nsc', 'skymapper')
    survey_coverage = {}
    for survey, specification in survey_products.items():
        columns = ['qso_group_id', specification['source_flag']]
        columns.extend((name for pair in specification['bands'].values() for name in pair))
        columns.extend((f'{survey}_{band}_band_quality' for band in specification['bands']))
        frame = pd.read_csv(
            specification['path'],
            usecols=columns,
            dtype={'qso_group_id': str},
            low_memory=False
        )
        frame = frame.set_index('qso_group_id').reindex(blend_candidates.qso_group_id)
        source_usable = saved_true(frame[specification['source_flag']])
        added_bands = np.zeros(system_count, dtype=int)
        for band, (magnitude_column, error_column) in specification['bands'].items():
            label = f'{survey}:{band}'
            magnitude = pd.to_numeric(frame[magnitude_column], errors='coerce')
            magnitude_error = pd.to_numeric(frame[error_column], errors='coerce')
            usable = source_usable & magnitude.between(-50, 90, inclusive='neither') & magnitude_error.gt(0) & magnitude_error.lt(90)
            usable &= saved_true(frame[f'{survey}_{band}_band_quality'])
            usable_rows = usable.to_numpy(dtype=bool)
            native_flux = np.full(system_count, np.nan)
            native_error = np.full(system_count, np.nan)
            native_flux[usable_rows] = 10 ** ((22.5 - magnitude.to_numpy(float)[usable_rows]) / 2.5)
            native_error[usable_rows] = np.log(10) / 2.5 * native_flux[usable_rows] * magnitude_error.to_numpy(float)[usable_rows]
            band_index = band_labels.index(label)
            multisurvey_flux[usable_rows, band_index] = native_flux[usable_rows]
            multisurvey_variance[usable_rows, band_index] = native_error[usable_rows] ** 2
            added_bands += usable_rows
        survey_coverage[survey] = int(np.count_nonzero(added_bands))
        blend_candidates[f'n_{survey}_bands'] = added_bands
    blend_candidates['n_extra_bands'] = sum((blend_candidates[f'n_{survey}_bands'] for survey in survey_products))
    blend_candidates['extra_surveys_used'] = [','.join((survey for survey in blend_total_surveys if blend_candidates.at[row_index, f'n_{survey}_bands'] > 0)) for row_index in range(system_count)]
    blend_candidates['n_decals_bands'] = np.isfinite(decals_flux).sum(axis=1)
    fallback_labels = ('ps1:r', 'nsc:r', 'skymapper:r')
    fallback_indices = [band_labels.index(label) for label in fallback_labels]
    fallback_flux = multisurvey_flux[:, fallback_indices]
    fallback_variance = multisurvey_variance[:, fallback_indices]
    fallback_observed = np.isfinite(fallback_flux) & np.isfinite(fallback_variance)
    fallback_observed &= fallback_variance > 0
    first_fallback = fallback_observed.argmax(axis=1)
    rows = np.arange(system_count)
    has_fallback = fallback_observed.any(axis=1)
    missing_legacy_r = blend_candidates.decals_system_assignment.eq('none').to_numpy()
    blend_candidates['fallback_r_band_input'] = np.where(missing_legacy_r & has_fallback, np.asarray(fallback_labels)[first_fallback], '')
    blend_candidates['fallback_r_snr'] = np.where(
        missing_legacy_r & has_fallback,
        fallback_flux[rows, first_fallback] / np.sqrt(fallback_variance[rows, first_fallback]),
        np.nan
    )
    blend_sky = SkyCoord(
        blend_candidates.primary_component_ra.to_numpy(float) * u.deg,
        blend_candidates.primary_component_dec.to_numpy(float) * u.deg
    )
    blend_candidates['galactic_abs_b_deg'] = np.abs(blend_sky.galactic.b.deg)
    multisurvey_photometry = Photometry(multisurvey_flux, multisurvey_variance, band_labels)
    print('Unresolved-class systems passed to BlendModel:', system_count)
    print('With Legacy r:', int((~missing_legacy_r).sum()))
    print(
        'Without Legacy r but with fallback r S/N >= 20:',
        int((missing_legacy_r & blend_candidates.fallback_r_snr.ge(20)).sum())
    )
    print('At Galactic |b| < 25 degrees:', int(blend_candidates.galactic_abs_b_deg.lt(25).sum()))
    print(
        'Additional total-flux survey coverage:',
        {survey: survey_coverage[survey] for survey in blend_total_surveys}
    )
    worker_count = int(os.environ.get('PAPER1_BLEND_WORKERS', '1'))
    if worker_count < 1:
        raise ValueError('PAPER1_BLEND_WORKERS must be positive')
    blend_candidates['redshift_error_used'] = redshift_errors(blend_candidates, 'primary_qso_zerr', 'primary_qso_catalog', 'primary_qso_redshift_kind')
    coordinates = blend_candidates[['primary_component_ra', 'primary_component_dec',
                                    'primary_qso_z', 'redshift_error_used']].to_numpy(float)
    blend_score_rows = []
    if worker_count == 1:
        for model_name, blend_model in blend_models.items():
            print(f'Scoring {model_name}, multisurvey: {len(blend_candidates)} blends', flush=True)
            blend_score_rows.extend(score_batch(blend_model, model_name, multisurvey_photometry,
                                                coordinates, blend_candidates.qso_group_id.tolist()))
    else:
        indices = np.arange(len(blend_candidates))
        batch_count = min(len(indices), worker_count * 3)
        batches = []
        for model_name, model_path in MODEL_BUNDLES.items():
            for batch in np.array_split(indices, batch_count):
                batches.append((model_path, model_name, multisurvey_photometry.subset(batch),
                                coordinates[batch], blend_candidates.qso_group_id.iloc[batch].tolist()))
        print(f'Scoring {len(batches)} independent blend batches with {worker_count} workers.', flush=True)
        with ProcessPoolExecutor(max_workers=worker_count, mp_context=get_context('spawn'),
                                 initializer=load_worker_dustmaps) as pool:
            for batch_index, rows in enumerate(pool.map(score_worker_batch, batches), start=1):
                blend_score_rows.extend(rows)
                if batch_index % batch_count == 0:
                    print(f'Completed {batch_index}/{len(batches)} blend batches.', flush=True)
    blend_scores = pd.DataFrame(blend_score_rows)
    blend_ranked = blend_candidates[[
        'qso_group_id',
        'primary_qso_catalog',
        'primary_qso_z',
        'redshift_scope',
        'targetid',
        'desi_release',
        'desi_brickid',
        'desi_objid',
        'scoring_route',
        'decals_system_assignment',
        'both_decals_r_assignments',
        'blend_r_snr',
        'fallback_r_band_input',
        'fallback_r_snr',
        'galactic_abs_b_deg',
        'n_decals_bands',
        'n_extra_bands',
        'extra_surveys_used',
        'one_local_source_near_both',
        'selected_decals_component_photometry_quality'
    ]].copy()
    for model_name in blend_models:
        scores = blend_scores.loc[blend_scores.model.eq(model_name)].set_index('qso_group_id')
        for field in (
            'status',
            'reference_band',
            'bands_used',
            'notes',
            'log_bf_qq_qs',
            'alpha_mean_qq',
            'alpha_mean_qs',
            'redshift_error',
            'finite_nodes',
            'clipped_nodes',
            'retained_weight'
        ):
            blend_ranked[f'{field}_{model_name}'] = blend_ranked.qso_group_id.map(scores[field])
        blend_ranked[f'p_qq_odds_1_{model_name}'] = blend_ranked.qso_group_id.map(scores['p_qq'])
    blend_ranked['fallback_band_used'] = blend_ranked.notes_magdep.map(lambda notes: next(
        (note.split(':', 1)[1] for note in notes if note.startswith('legacy_r_marginalised:')),
        ''
    ) if isinstance(notes, (tuple, list)) else '')
    blend_ranked['rank_blend'] = blend_ranked.groupby('redshift_scope').log_bf_qq_qs_magdep.rank(ascending=False, method='min')
    blend_ranked = blend_ranked.sort_values(['redshift_scope', 'rank_blend', 'qso_group_id'], na_position='last')
    blend_scores.to_csv(ONE_ROOT / 'blend_scores_all_nonresolved_by_model_and_band_set.csv', index=False)
    blend_ranked.to_csv(ONE_ROOT / 'ranked_unresolved.csv', index=False)
    blend_candidates.to_csv(ONE_ROOT / 'unresolved_photometry_audit.csv', index=False)
    print(
        'Output:',
        ONE_ROOT,
        '| diagnostic rows:',
        len(blend_ranked),
        '| finite multisurvey:',
        {name: int(np.isfinite(blend_ranked[f'log_bf_qq_qs_{name}']).sum()) for name in MODEL_BUNDLES}
    )


if __name__ == '__main__':
    main()
