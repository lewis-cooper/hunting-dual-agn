"""Reuse production likelihoods and flags, with explicit control-photometry provenance."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from control_catalogues import read_csv, magnitude, pm_significance,validate_cone_coverage

HEMISPHERE = {9010:'decals_dr9_south',9011:'decals_dr9_north',9012:'decals_dr9_south'}


def match_catalogues(tractor, gaia):
    """One record per Tractor source; preserve unmatched and ambiguous associations."""
    rows=tractor.drop_duplicates(['target_id','release','brickid','objid']).copy().reset_index(drop=True)
    rows['control_id']=(rows.release.astype(str)+':'+rows.brickid.astype(str)+':'+rows.objid.astype(str))
    rows['gaia_source_id']=pd.Series(pd.NA,index=rows.index,dtype='string')
    rows['gaia_g']=np.nan;rows['pm_sig']=np.nan
    sources=gaia.drop_duplicates('source_id').reset_index(drop=True)
    if len(rows) and len(sources):
        source_sky=SkyCoord(sources.ra.to_numpy(float)*u.deg,sources.dec.to_numpy(float)*u.deg)
        index,separation,_=SkyCoord(rows.ra.to_numpy(float)*u.deg,rows.dec.to_numpy(float)*u.deg).match_to_catalog_sky(source_sky)
        matched=separation.arcsec<=.5
        associated=sources.iloc[index[matched]]
        rows.loc[matched,'gaia_source_id']=associated.source_id.astype('string').to_numpy()
        rows.loc[matched,'gaia_g']=associated.phot_g_mean_mag.to_numpy()
        rows.loc[matched,'pm_sig']=pm_significance(associated)
    rows['gaia_match_ambiguous']=rows.gaia_source_id.notna() & rows.duplicated(['target_id','gaia_source_id'],keep=False)
    rows['has_gaia']=rows.gaia_source_id.notna() & ~rows.gaia_match_ambiguous
    rows['photometric_system']=rows.release.map(HEMISPHERE)
    rows['r_mag']=magnitude(rows.flux_r)
    return rows


def component_scores(frame, models, threshold, role='secondary',skip_known_failures=False):
    """Native grz policy, per-object Quaia error, support test and both model variants."""
    from dual_agn.one_arcsecond import score_resolved_candidates
    if not np.isfinite(frame.redshift_error).all() or frame.redshift_error.lt(0).any():
        raise ValueError('Controls need finite, non-negative redshift uncertainties')
    model_input=frame.copy()
    model_input['qso_group_id']=frame.trial_id
    model_input['primary_qso_z']=frame.redshift
    model_input['primary_qso_zerr']=frame.redshift_error
    model_input['primary_qso_redshift_kind']=np.where(frame.redshift_error.gt(0),'photometric','spectroscopic')
    model_input['primary_qso_catalog']=np.where(
        frame.redshift_catalogue.str.contains('quaia',case=False,na=False),'QUAIA_G20_V1','SPECTROSCOPIC')
    model_input['primary_component_gaia_source_id']=frame.target_id
    model_input['secondary_component_gaia_source_id']=frame.control_id
    model_input['secondary_component_ra']=frame.ra
    model_input['secondary_component_dec']=frame.dec
    model_input['secondary_decals_type']=frame.type
    model_input['component_separation_arcsec']=frame.separation_arcsec
    model_input['secondary_decals_fracflux_r']=frame.fracflux_r
    for band in 'grz':
        model_input[f'secondary_decals_flux_{band}']=frame[f'flux_{band}']
        ivar=frame[f'flux_ivar_{band}']
        if role=='imaging':
            usable=frame.maskbits.eq(0) & frame[f'nobs_{band}'].gt(0)
            ivar=ivar.where(usable,0.)
        model_input[f'secondary_decals_flux_ivar_{band}']=ivar
    scored=frame.copy()
    reference_snr=frame.flux_r*np.sqrt(frame.flux_ivar_r.where(frame.flux_ivar_r.gt(0)))
    latitude=np.abs(SkyCoord(frame.ra.to_numpy(float)*u.deg,
                             frame.dec.to_numpy(float)*u.deg).galactic.b.deg)
    scored['reference_snr']=reference_snr
    scored['reference_snr_pass']=np.isfinite(reference_snr)&reference_snr.ge(10)
    scored['psf_pass']=frame.type.str.strip().eq('PSF')
    scored['latitude_pass']=np.isfinite(latitude)&(latitude>=25)
    usable=model_input.photometric_system.notna()
    if skip_known_failures:
        usable &= scored.reference_snr_pass&scored.psf_pass&scored.latitude_pass
        if 'secondary_pm_pass' in frame:
            usable &= frame.secondary_pm_pass.eq(True)
        if 'preselection_pass' in frame:
            usable &= frame.preselection_pass.eq(True)
    result=score_resolved_candidates(model_input.loc[usable],models) if usable.any() else pd.DataFrame()
    for model in models:
        subset=result.loc[result.model.eq(model)].set_index('qso_group_id') if len(result) else pd.DataFrame()
        for field in ['log_r_per_unit_z','eligible','status','quality_flags']:
            scored[f'{field}_{model}']=scored.trial_id.map(subset[field]) if len(subset) else np.nan
            if field in ['status','quality_flags']:
                scored[f'{field}_{model}']=scored[f'{field}_{model}'].astype('string')
        if skip_known_failures:
            scored.loc[~usable,f'status_{model}']='not_scored_known_preselection_failure'
        scored[f'score_pass_{model}']=scored[f'eligible_{model}'].eq(True)
        scored[f'score_pass_{model}'] &= np.isfinite(scored[f'log_r_per_unit_z_{model}'])
        scored[f'score_pass_{model}'] &= scored[f'log_r_per_unit_z_{model}'].ge(threshold)
        scored[f'score_pass_{model}'] &= scored.reference_snr_pass & scored.psf_pass & scored.latitude_pass
    return scored


def cached_scores(frame,models,threshold,role,root,chunk_size=2500):
    path=root/f'{role}_model_score_cache.parquet'
    cached=pd.read_parquet(path) if path.exists() else pd.DataFrame()
    known=set(cached.trial_id) if len(cached) else set()
    missing=frame.loc[~frame.trial_id.isin(known)]
    outputs=[cached] if len(cached) else []
    for start in range(0,len(missing),chunk_size):
        batch=missing.iloc[start:start+chunk_size]
        scored=component_scores(batch,models,threshold,role=role,skip_known_failures=True)
        fields=['trial_id']+[column for column in scored if column not in batch]
        outputs.append(scored[fields])
        pd.concat(outputs,ignore_index=True).to_parquet(path,index=False,compression='zstd')
        print(f'{role} component models: {min(start+chunk_size,len(missing))}/{len(missing)} new rows',flush=True)
    cache=pd.concat(outputs,ignore_index=True) if outputs else pd.DataFrame()
    if not len(frame):
        return frame
    return frame.merge(cache,on='trial_id',how='left',validate='one_to_one')


def configure_production(args):
    """Redirect every production runtime write before importing the scorer."""
    os.environ['PAPER1_OUTPUT_ROOT']=str(args.output.resolve()/'scorer_runtime')
    os.environ['PAPER1_WORK_ROOT']=str(args.output.resolve()/'scorer_runtime')
    os.environ['PAPER1_MODEL_ROOT']=str(args.models.resolve())
    os.environ['QSO_PCOLOR_DUSTMAPS_DIR']=str(args.dustmaps.resolve())
    os.environ['MPLCONFIGDIR']=str(args.output.resolve()/'scorer_runtime/matplotlib')
    package_source=args.models.resolve().parents[1]/'src'
    if (package_source/'qso_pcolor').is_dir():
        sys.path.insert(0,str(package_source))
    sys.path.insert(0,str(args.repo.resolve()))
    runtime = args.output.resolve()/'scorer_runtime'
    runtime.mkdir(parents=True, exist_ok=True)
    # Workers must not truncate the same configuration file while another reads it.
    with NamedTemporaryFile(mode='w', dir=runtime, prefix='dustmaps_', suffix='.json',
                            delete=False) as stream:
        json.dump({'data_dir': str(args.dustmaps.resolve())}, stream)
        configuration_file = stream.name
    os.environ['DUSTMAPS_CONFIG_FNAME'] = configuration_file
    from dustmaps.config import config
    config.fname = configuration_file
    config.load()


def verify_production_scorer(args):
    from dual_agn.provenance import verify_imported_scorer
    return verify_imported_scorer(package_root=args.models.resolve().parents[1])


def score_blends(trials,args,models):
    from qso_pcolor import Photometry
    from qso_pcolor.blend import BlendModel
    from dual_agn.redshift_uncertainty import blend_score_with_redshift_error
    if not args.blend_photometry and not args.blend_grz_diagnostic:
        return pd.DataFrame()
    supplied=read_csv(args.blend_photometry) if args.blend_photometry else None
    if supplied is not None:
        required={'trial_id','component','band','flux','variance','independent'}
        if not required.issubset(supplied):raise ValueError(f'Blend photometry requires {sorted(required)}')
        if not supplied.independent.eq(True).all():
            raise ValueError('Supplied blend components must be certified independent measurements')
        if supplied.duplicated(['trial_id','component','band']).any():raise ValueError('Duplicate blend bands')
    path=args.output/'blend_control_scores.csv'
    previous=read_csv(path) if getattr(args,'resume',False) and path.exists() else pd.DataFrame()
    rows=previous.to_dict('records')
    completed=set(previous.trial_id) if len(previous) else set()
    blend_models={}
    def save_rows():
        temporary=path.with_suffix('.tmp')
        pd.DataFrame(rows).to_csv(temporary,index=False)
        temporary.replace(path)
    for number,trial in enumerate(trials.itertuples(),1):
        if trial.trial_id in completed:
            continue
        if args.max_blend_controls and number>args.max_blend_controls:break
        row={'trial_id':trial.trial_id,'diagnostic_only':supplied is None,
             'runtime_capped':bool(args.max_blend_controls)}
        if getattr(trial,'selection_mode','')=='unmatched':
            row.update(diagnostic_only=True,photometric_pass=False,missing_sum_bands='unmatched_control_population',
                       photometry_scope='no matched independent control')
            for name in models:
                row[f'log_bf_qq_qs_{name}']=np.nan
                row[f'status_{name}']='unknown_unmatched_control_population'
            rows.append(row)
            continue
        if getattr(args,'skip_known_pm_failures',False) and not trial.secondary_pm_pass:
            row.update(diagnostic_only=False,photometric_pass=False,missing_sum_bands='',
                       photometry_scope='known secondary-PM rejection; flux sum unnecessary')
            for name in models:
                row[f'log_bf_qq_qs_{name}']=np.nan
                row[f'status_{name}']='not_scored_secondary_pm_failure'
            rows.append(row)
            continue
        if not np.isfinite(trial.redshift_error) or trial.redshift_error<0:
            raise ValueError('Blended controls need a measured finite, non-negative redshift uncertainty')
        for name in models:
            if name not in blend_models:
                blend_models[name]=BlendModel.load(args.models/f'20261005_{name}_qlf')
                if getattr(args,'cache_population_priors',False):
                    from dual_agn.blend_priors import cache_population_priors
                    blend_models[name]=cache_population_priors(blend_models[name])
            model=blend_models[name]
            flux=np.full((1,len(model.bands)),np.nan); variance=np.full_like(flux,np.inf)
            if supplied is not None:
                measurements=supplied.loc[supplied.trial_id.eq(trial.trial_id)]
                primary=measurements.loc[measurements.component.eq('primary')].set_index('band')
                secondary=measurements.loc[measurements.component.eq('secondary')].set_index('band')
                # The frozen production inputs query NSC/SkyMapper after a missing PS1 match.
                uses_ps1=any(band.startswith('ps1:') for band in primary.index)
                permitted=lambda band: band.startswith('decals_dr9_') or (
                    band.startswith('ps1:') if uses_ps1 else band.startswith(('nsc:','skymapper:')))
                if getattr(args,'production_legacy_first',False):
                    if trial.photometry_class=='legacy':
                        permitted=lambda band:band.startswith('decals_dr9_')
                    else:
                        permitted=lambda band:(band.startswith('ps1:') if uses_ps1 else band.startswith(('nsc:','skymapper:')))
                if getattr(args,'fixed_survey_policy',False):
                    prefixes=('decals_dr9_',) if trial.band_policy=='legacy' else tuple(
                        survey+':' for survey in str(trial.band_policy).split(','))
                    permitted=lambda band:band.startswith(prefixes)
                intended=[band for band in primary.index if permitted(band)]
                row['missing_sum_bands']=','.join(band for band in intended if band not in secondary.index)
                row['diagnostic_only']=bool(row['missing_sum_bands'])
                common=[band for band in primary.index.intersection(secondary.index) if permitted(band)]
                if not any(band.endswith(':r') for band in common):
                    row['diagnostic_only']=True
                    row['missing_sum_bands']+=',unmeasured_reference_band'
                for band in primary.index.intersection(secondary.index):
                    if not permitted(band):continue
                    if band not in model.bands:raise ValueError(f'Unknown native band {band}')
                    if band.startswith(('ps1:','nsc:','skymapper:')) and getattr(args,'summed_band_quality',False):
                        quality_flux=primary.at[band,'quality_flux']+secondary.at[band,'quality_flux']
                        quality_variance=primary.at[band,'quality_variance']+secondary.at[band,'quality_variance']
                        quality_mag=22.5-2.5*np.log10(quality_flux) if quality_flux>0 else np.nan
                        quality_error=2.5/np.log(10)*np.sqrt(quality_variance)/quality_flux if quality_flux>0 else np.inf
                        if not (np.isfinite(quality_mag+quality_error) and 0<quality_mag<40 and quality_error<=.3):
                            continue
                    index=model.bands.index(band)
                    flux[0,index]=primary.at[band,'flux']+secondary.at[band,'flux']
                    variance[0,index]=primary.at[band,'variance']+secondary.at[band,'variance']
                row['photometry_scope']='independent components; production Legacy-first/fallback survey policy; common measured bands'
            else:
                for band in 'grz':
                    if trial.primary_release not in HEMISPHERE or pd.isna(trial.photometric_system):continue
                    if HEMISPHERE[trial.primary_release]!=trial.photometric_system:continue
                    first_ivar=getattr(trial,f'primary_ivar_{band}')
                    second_ivar=getattr(trial,f'flux_ivar_{band}')
                    if not (np.isfinite(first_ivar+second_ivar) and first_ivar>0 and second_ivar>0):continue
                    index=model.bands.index(f'{trial.photometric_system}:{band}')
                    flux[0,index]=getattr(trial,f'primary_flux_{band}')+getattr(trial,f'flux_{band}')
                    variance[0,index]=1/first_ivar+1/second_ivar
                row['photometry_scope']='grz_only; primary may include existing companion; diagnostic only'
            if row['diagnostic_only'] and getattr(args,'skip_incomplete_photometry',False):
                row[f'log_bf_qq_qs_{name}']=np.nan
                row[f'status_{name}']='unknown_incomplete_independent_photometry'
                row[f'bands_used_{name}']=''
                continue
            score,audit,grid=blend_score_with_redshift_error(model,
                Photometry(flux,variance,model.bands),trial.primary_ra,trial.primary_dec,
                trial.redshift,trial.redshift_error)
            row[f'log_bf_qq_qs_{name}']=score.log_bf_qq_qs
            row[f'status_{name}']=score.status
            row[f'bands_used_{name}']=','.join(score.bands_used)
            row[f'redshift_grid_{name}']=json.dumps(grid)
        values=np.array([row[f'log_bf_qq_qs_{name}'] for name in models])
        row['photometric_pass']=bool(np.isfinite(values).all() and (values>0).any())
        rows.append(row)
        save_rows()
        if number%getattr(args,'progress_every',1)==0:
            print(f'Blend control {number}/{len(trials)}; exact components={supplied is not None}',flush=True)
    save_rows()
    (args.output/'blend_scoring_configuration.json').write_text(json.dumps({
        'models':str(args.models.resolve()),'repo':str(args.repo.resolve()),
        'input_photometry':str(args.blend_photometry.resolve()) if args.blend_photometry else None,
        'input_photometry_sha256':hashlib.sha256(args.blend_photometry.read_bytes()).hexdigest() if args.blend_photometry else None,
        'scorer_sha256':hashlib.sha256((args.repo/'dual_agn/redshift_uncertainty.py').read_bytes()).hexdigest(),
        'n_input_controls':len(trials),'n_scored_controls':len(rows),'runtime_cap':args.max_blend_controls,
        'selection':'finite Bayes factors in both variants; QQ favoured in either; secondary PM applied in summary',
        'photometry_policy':'Legacy plus PS1 when present for primary, otherwise NSC/SkyMapper; summed independent native fluxes and variances',
        'redshift_policy':'production per-object uncertainty dispatcher',
        'population_prior_memoization':getattr(args,'cache_population_priors',False)},indent=2)+'\n')
    return pd.DataFrame(rows)


