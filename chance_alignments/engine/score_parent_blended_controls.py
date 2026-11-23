#!/usr/bin/env python3
"""Score independent flux sums in resumable workers with unchanged production models."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from score_controls import configure_production,score_blends
from control_catalogues import read_csv


def score_worker(payload):
    trials,arguments=payload
    args=SimpleNamespace(**arguments)
    args.output.mkdir(parents=True,exist_ok=True)
    configure_production(args)
    return score_blends(trials,args,{'magdep':None,'magindep':None})


def selected_trials(args):
    root=args.output
    trials=(pd.read_parquet(root/'blended_control_trials.parquet') if (root/'blended_control_trials.parquet').exists()
            else read_csv(root/'isolated_blend_trials.csv'))
    if args.draws_per_input:
        indices=pd.to_numeric(trials.trial_id.str.extract(r':(\d+)(?::unmatched)?$',expand=False),errors='coerce')
        limits=pd.Series(args.draws_per_input,index=trials.index)
        if getattr(args,'legacy_draws_per_input',0):
            limits.loc[trials.photometry_class.eq('legacy')]=args.legacy_draws_per_input
        trials=trials.loc[indices.isna()|indices.lt(limits)].copy()
        totals=trials.groupby(['input_pair_id','annulus']).mock_weight.transform('sum')
        trials['mock_weight']/=totals
    return trials.sample(frac=1,random_state=20261121)


def write_scores(scores,trials,args):
    root=args.output
    scores=scores.loc[scores.trial_id.isin(trials.trial_id)].copy()
    if len(scores)!=len(trials) or not scores.trial_id.is_unique:
        raise ValueError('Not every independent blended trial was scored or classified')
    metadata=trials.drop(columns=trials.columns.intersection(scores.columns).difference(['trial_id']))
    scores=scores.merge(metadata,on='trial_id',validate='one_to_one')
    scores['individual_cuts_pass']=scores.photometric_pass.eq(True)&scores.secondary_pm_pass.eq(True)
    destination=root/args.score_file
    temporary=destination.with_suffix('.tmp')
    scores.to_parquet(temporary,index=False,compression='zstd')
    temporary.replace(destination)
    destination.with_suffix('.json').write_text(json.dumps(dict(
        n_trials=len(trials),n_scored_or_known_rejected=len(scores),workers=args.workers,
        fallback_draws_per_input=getattr(args,'draws_per_input',0),
        legacy_draws_per_input=getattr(args,'legacy_draws_per_input',0) or getattr(args,'draws_per_input',0),
        n_unknown_photometry=int(scores.diagnostic_only.eq(True).sum()),
        policy='production QQ/QS; finite in both, positive in either; secondary PM applied',
        quality='native total fluxes; PSF quality on summed flux; independent variances added',
        shortcut='PM failures are known rejections; no claim of computed Bayes factors for these rows'),indent=2)+'\n')
    print(scores.groupby('annulus').agg(trials=('trial_id','size'),unknown=('diagnostic_only','sum'),
        passing=('individual_cuts_pass','sum')).to_string(),flush=True)


def score(args):
    root=args.output
    os.environ['OPENBLAS_NUM_THREADS']='1'
    os.environ['OMP_NUM_THREADS']='1'
    trials=selected_trials(args)
    assignment=pd.util.hash_pandas_object(trials.trial_id,index=False).to_numpy()%args.workers
    jobs=[]
    for worker in range(args.workers):
        if not (assignment==worker).any():
            continue
        arguments=dict(vars(args),output=root/'blend_workers'/f'worker_{worker}',
            blend_photometry=root/'independent_blend_photometry.csv',blend_grz_diagnostic=False,
            max_blend_controls=0,resume=True,skip_known_pm_failures=True,summed_band_quality=True)
        arguments['production_legacy_first']=True
        arguments['skip_incomplete_photometry']=True
        arguments['fixed_survey_policy']='band_policy' in trials
        arguments['cache_population_priors']=True
        arguments['progress_every']=25
        jobs.append((trials.loc[assignment==worker].copy(),arguments))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        outputs=list(pool.map(score_worker,jobs))
    scores=pd.concat(outputs,ignore_index=True)
    write_scores(scores,trials,args)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--models',type=Path,required=True)
    parser.add_argument('--dustmaps',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--draws-per-input',type=int,default=0)
    parser.add_argument('--legacy-draws-per-input',type=int,default=0)
    parser.add_argument('--score-file',default='blend_control_scores.parquet')
    score(parser.parse_args())


if __name__=='__main__':
    main()
