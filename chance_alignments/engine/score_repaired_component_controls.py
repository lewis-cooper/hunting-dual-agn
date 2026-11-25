#!/usr/bin/env python3
"""Score re-paired sources at the actual input redshifts and fixed primary fluxes."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from control_catalogues import magnitude
from score_controls import configure_production
from score_controls import cached_scores


def score_worker(payload):
    arguments, input_path, cache_root, threshold, role = payload
    args = SimpleNamespace(**arguments)
    configure_production(args)
    from qso_pcolor.unified import UnifiedPSFModel
    models = {name: UnifiedPSFModel.load(args.models/f'20261005_{name}_qlf')
              for name in ['magdep','magindep']}
    frame = pd.read_parquet(input_path)
    cached_scores(frame, models, threshold, role, cache_root)
    return cache_root/f'{role}_model_score_cache.parquet'


def parallel_scores(frame, models, threshold, role, cache_root, args):
    path = cache_root/f'{role}_model_score_cache.parquet'
    previous = pd.read_parquet(path) if path.exists() else pd.DataFrame()
    known = set(previous.trial_id) if len(previous) else set()
    missing = frame.loc[~frame.trial_id.isin(known)].copy()
    if not len(missing):
        return frame.merge(previous, on='trial_id', validate='one_to_one')
    partition = pd.util.hash_pandas_object(missing.trial_id, index=False).to_numpy()%args.workers
    jobs = []
    for worker in range(args.workers):
        subset = missing.loc[partition == worker]
        if not len(subset):
            continue
        worker_root = cache_root/f'{role}_workers/worker_{worker}'
        worker_root.mkdir(parents=True, exist_ok=True)
        input_path = worker_root/'input.parquet'
        subset.to_parquet(input_path, index=False, compression='zstd')
        jobs.append((dict(vars(args)), input_path, worker_root, threshold, role))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        paths = list(pool.map(score_worker, jobs))
    complete = pd.concat([previous, *[pd.read_parquet(part) for part in paths]], ignore_index=True)
    if not complete.trial_id.is_unique:
        raise ValueError('Duplicate cached re-paired scores')
    temporary = path.with_suffix('.tmp')
    complete.to_parquet(temporary, index=False, compression='zstd')
    temporary.replace(path)
    return frame.merge(complete, on='trial_id', how='left', validate='one_to_one')


def score(args):
    os.environ['OPENBLAS_NUM_THREADS']='1'
    os.environ['OMP_NUM_THREADS']='1'
    configure_production(args)
    from qso_pcolor.unified import UnifiedPSFModel
    from dual_agn.score_selection import read_either_model_threshold
    root=args.output
    models={name:UnifiedPSFModel.load(args.models/f'20261005_{name}_qlf') for name in ['magdep','magindep']}
    threshold=read_either_model_threshold(
        getattr(args,'calibration',None) or args.repo/'products/threshold_calibration.json')
    trials=pd.read_parquet(root/'component_control_trials.parquet')
    if args.routes:
        trials=trials.loc[trials.route.isin(args.routes)].copy()
    primary=pd.read_parquet(root/'fixed_primary_scores.parquet').set_index('pair_id')
    outputs=[]
    for role,routes in [('secondary',['resolved']),('imaging',['imaging_one_gaia','imaging_two_gaia'])]:
        frame=trials.loc[trials.route.isin(routes)].copy()
        if not len(frame):
            continue
        unique=frame.drop_duplicates('model_key').copy()
        unique['trial_id']=unique.model_key
        cache_root=root/'score_cache/repaired';cache_root.mkdir(parents=True,exist_ok=True)
        unique=parallel_scores(unique,models,threshold,role,cache_root,args) if args.workers>1 else cached_scores(unique,models,threshold,role,cache_root)
        fields=[column for column in unique if column.startswith(('log_r_per_unit_z_','eligible_',
            'status_','quality_flags_','score_pass_'))]
        fields+=['reference_snr','reference_snr_pass','psf_pass','latitude_pass']
        frame=frame.merge(unique[['model_key',*fields]],on='model_key',how='left',validate='many_to_one')
        if role=='imaging':
            for name in models:
                frame[f'primary_pass_{name}']=frame.input_pair_id.map(primary[f'score_pass_{name}']).eq(True)
                frame[f'score_pass_{name}'] &= frame[f'primary_pass_{name}']
            primary_gr=magnitude(frame.primary_flux_g)-magnitude(frame.primary_flux_r)
            secondary_gr=magnitude(frame.flux_g)-magnitude(frame.flux_r)
            frame['delta_gr_pass']=np.isfinite(primary_gr+secondary_gr)&(np.abs(primary_gr-secondary_gr)<=1)
        frame['photometric_pass']=frame[['score_pass_magdep','score_pass_magindep']].any(axis=1)
        frame['individual_cuts_pass']=frame.photometric_pass&frame.secondary_pm_pass&~frame.known_primary_rejection
        if role=='imaging':
            frame['individual_cuts_pass'] &= frame.delta_gr_pass
        outputs.append(frame)
    scored=pd.concat(outputs,ignore_index=True)
    destination=root/args.score_file
    scored.to_parquet(destination,index=False,compression='zstd')
    destination.with_suffix('.json').write_text(json.dumps(dict(
        threshold=threshold,repo=str(args.repo.resolve()),models=str(args.models.resolve()),
        scorer_sha256=hashlib.sha256((args.repo/'dual_agn/one_arcsecond.py').read_bytes()).hexdigest(),
        redshift_policy_sha256=hashlib.sha256((args.repo/'dual_agn/redshift_uncertainty.py').read_bytes()).hexdigest(),
        primary='actual real input primaries; no substitution of donor redshifts or primary SEDs',
        weighting='per-input empirical source matching; reused sources retain their bootstrap identities'),indent=2)+'\n')
    print(scored.groupby(['route','annulus']).agg(mocks=('trial_id','size'),
        passing=('individual_cuts_pass','sum'),unknown=('diagnostic_only','sum')).to_string(),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--models',type=Path,required=True)
    parser.add_argument('--dustmaps',type=Path,required=True)
    parser.add_argument('--calibration',type=Path)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--routes',nargs='+')
    parser.add_argument('--score-file',default='component_control_scores.parquet')
    score(parser.parse_args())


if __name__=='__main__':
    main()
