#!/usr/bin/env python3
"""Apply the unchanged production component models to scalable parent controls."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from score_controls import component_scores,configure_production
from control_catalogues import magnitude


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


def score(args):
    configure_production(args)
    from qso_pcolor.unified import UnifiedPSFModel
    from score_selection import read_either_model_threshold
    root=args.output
    threshold=read_either_model_threshold(
        getattr(args,'calibration',None) or args.repo/'products/threshold_calibration.json')
    models={name:UnifiedPSFModel.load(args.models/f'20261005_{name}_qlf') for name in ['magdep','magindep']}
    trials=pd.read_parquet(root/'parent_control_trials.parquet')
    primary=pd.read_parquet(root/'parent_primary_photometry.parquet')
    primary=primary.loc[primary.primary_unique&primary.primary_type.eq('PSF')].copy()
    primary['trial_id']='primary:'+primary.target_id
    primary['control_id']=primary.primary_control_id
    primary['type']=primary.primary_type
    primary['photometric_system']=primary.primary_photometric_system
    primary['fracflux_r']=np.nan
    primary['separation_arcsec']=primary.target_id.map(trials.groupby('target_id').separation_arcsec.min())
    primary=primary.loc[np.isfinite(primary.separation_arcsec)].copy()
    primary['maskbits']=primary.primary_maskbits
    for band in 'grz':
        primary[f'flux_{band}']=primary[f'primary_flux_{band}']
        primary[f'flux_ivar_{band}']=primary[f'primary_flux_ivar_{band}']
        primary[f'nobs_{band}']=primary[f'primary_nobs_{band}']
    primary['primary_pm_pass']=np.isfinite(primary.primary_pm_sig)&primary.primary_pm_sig.lt(3)
    primary['secondary_pm_pass']=primary.primary_pm_pass
    cache_root=root/'score_cache';cache_root.mkdir(exist_ok=True)
    (cache_root/'primary').mkdir(exist_ok=True)
    primary=cached_scores(primary,models,threshold,'imaging',cache_root/'primary')
    primary.to_parquet(root/'primary_component_scores.parquet',index=False,compression='zstd')
    primary=primary.set_index('target_id')
    resolved=trials.loc[trials.route.eq('resolved')].copy()
    resolved=cached_scores(resolved,models,threshold,'secondary',cache_root)
    imaging=trials.loc[trials.route.str.startswith('imaging')].copy()
    imaging['model_key']='imaging:'+imaging.target_id+':'+imaging.control_id
    unique=imaging.drop_duplicates('model_key').copy()
    unique['trial_id']=unique.model_key
    unique=unique.drop(columns=['secondary_pm_pass'])
    primary_any=primary[['score_pass_magdep','score_pass_magindep']].any(axis=1)
    unique['preselection_pass']=unique.target_id.map(primary_any).eq(True)
    unique=cached_scores(unique,models,threshold,'imaging',cache_root)
    model_fields=[column for column in unique if column.startswith(('log_r_per_unit_z_','eligible_',
        'status_','quality_flags_','score_pass_'))]
    model_fields+=['reference_snr','reference_snr_pass','psf_pass','latitude_pass']
    imaging=imaging.merge(unique[['model_key',*model_fields]],on='model_key',how='left',validate='many_to_one')
    for name in models:
        imaging[f'primary_pass_{name}']=imaging.target_id.map(primary[f'score_pass_{name}']).eq(True)
        imaging[f'score_pass_{name}'] &= imaging[f'primary_pass_{name}']
    primary_gr=magnitude(imaging.primary_flux_g)-magnitude(imaging.primary_flux_r)
    secondary_gr=magnitude(imaging.flux_g)-magnitude(imaging.flux_r)
    imaging['delta_gr_pass']=np.isfinite(primary_gr+secondary_gr)&(np.abs(primary_gr-secondary_gr)<=1)
    imaging['primary_pm_pass']=np.isfinite(imaging.primary_pm_sig)&imaging.primary_pm_sig.lt(3)
    scored=pd.concat([resolved,imaging],ignore_index=True)
    scored['photometric_pass']=scored[['score_pass_magdep','score_pass_magindep']].any(axis=1)
    scored['individual_cuts_pass']=scored.photometric_pass&scored.secondary_pm_pass
    image_mask=scored.route.str.startswith('imaging')
    scored.loc[image_mask,'individual_cuts_pass'] &= scored.loc[image_mask,'delta_gr_pass'].eq(True)
    scored.loc[image_mask,'individual_cuts_pass'] &= scored.loc[image_mask,'primary_pm_pass'].eq(True)
    scored.to_parquet(root/'component_control_scores.parquet',index=False,compression='zstd')
    configuration=dict(threshold=threshold,repo=str(args.repo.resolve()),models=str(args.models.resolve()),
        primary_scorer_sha256=hashlib.sha256((args.repo/'one_arcsecond.py').read_bytes()).hexdigest(),
        redshift_policy_sha256=hashlib.sha256((args.repo/'redshift_uncertainty.py').read_bytes()).hexdigest(),
        rule='either model; same variant for both imaging components; secondary PM<3 when finite; missing retained',
        known_failure_shortcut='fixed calibrated threshold; deterministic PM/PSF/SNR/primary failures retained with failure status')
    (root/'component_scoring_configuration.json').write_text(json.dumps(configuration,indent=2)+'\n')
    print(scored.groupby(['route','annulus']).individual_cuts_pass.agg(['size','sum']).to_string(),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--models',type=Path,required=True)
    parser.add_argument('--dustmaps',type=Path,required=True)
    parser.add_argument('--calibration',type=Path)
    score(parser.parse_args())


if __name__=='__main__':
    main()
