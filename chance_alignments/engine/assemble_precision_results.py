#!/usr/bin/env python3
"""Combine completed route-specific batches, preserving per-input mock weights."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pandas as pd


def assemble(args):
    main_path=args.blended/args.blended_score_file
    extension_path=args.extension/'blend_control_scores.parquet'
    legacy_path=args.legacy_batch/'blend_control_scores.parquet' if args.legacy_batch else None
    ready=lambda:main_path.exists() and extension_path.exists() and (legacy_path is None or legacy_path.exists())
    waiting=not ready()
    if waiting:
        print('Waiting for the two independent blended-score batches to complete.',flush=True)
    while not ready():
        time.sleep(15)
    main=pd.read_parquet(main_path)
    extension=pd.read_parquet(extension_path)
    main['control_source_batch']=str(args.blended)
    legacy_ids=set()
    if legacy_path:
        legacy=pd.read_parquet(legacy_path)
        if not legacy.photometry_class.eq('legacy').all():
            raise ValueError('The Legacy precision extension contains another survey policy')
        legacy_ids=set(legacy.input_pair_id)
        legacy['control_source_batch']=str(args.legacy_batch)
        main=pd.concat([main.loc[~main.input_pair_id.isin(legacy_ids)],legacy],ignore_index=True)
        extension=extension.loc[~extension.input_pair_id.isin(legacy_ids)].copy()
    identities=set(extension.input_pair_id)
    retained=main.loc[~main.input_pair_id.isin(identities)].copy()
    extension['control_source_batch']=str(args.extension)
    scores=pd.concat([retained,extension],ignore_index=True)
    if not scores.trial_id.is_unique:
        raise ValueError('A blended mock has been duplicated across batches')
    np.testing.assert_allclose(scores.groupby(['input_pair_id','annulus']).mock_weight.sum(),1.)
    pairs=pd.read_csv(args.output/'input_pairs.csv')
    expected=set(pairs.loc[pairs.route.eq('blended')&pairs.latitude_bin.notna(),'pair_id'])
    for annulus in ['3-10','10-30']:
        if set(scores.loc[scores.annulus.eq(annulus),'input_pair_id'])!=expected:
            raise ValueError('The assembled blended controls do not represent every eligible input')
    scores.to_parquet(args.output/'blend_control_scores.parquet',index=False,compression='zstd')
    old_audit=pd.read_csv(args.blended/'blended_matching_audit.csv')
    if legacy_path:
        legacy_audit=pd.read_csv(args.legacy_batch/'blended_matching_audit.csv')
        old_audit=pd.concat([old_audit.loc[~old_audit.input_pair_id.isin(legacy_ids)],legacy_audit],ignore_index=True)
    new_audit=pd.read_csv(args.extension/'blended_matching_audit.csv')
    new_audit=new_audit.loc[new_audit.input_pair_id.isin(identities)]
    audit=pd.concat([old_audit.loc[~old_audit.input_pair_id.isin(identities)],new_audit],ignore_index=True)
    audit.to_csv(args.output/'blended_matching_audit.csv',index=False)
    (args.output/'control_batches.json').write_text(json.dumps(dict(
        component_extension=str(args.output.resolve()),retained_one_gaia_batch=str(args.blended.resolve()),
        blended_main_score_file=str(main_path.resolve()),blended_coverage_extension=str(extension_path.resolve()),
        legacy_precision_batch=str(legacy_path.resolve()) if legacy_path else None,
        coverage_allocation='only inputs lacking main-annulus matches were checked in the larger field library; no allocation by score outcome',
        mock_draws_per_input=dict(legacy=args.legacy_draws_per_input or args.draws_per_input,
                                 fallback=args.draws_per_input,coverage_extension=20),
        interpretation='all rates remain conditional on the documented matching and photometry-transfer assumptions'),indent=2)+'\n')
    scripts=Path(__file__).resolve().parent
    commands=[('summarize_repaired_controls.py',['--output',str(args.output),
        '--bootstrap',str(args.bootstrap)]),
        ('write_precision_report.py',['--output',str(args.output)])]
    for script,arguments in commands:
        subprocess.run([sys.executable,'-B',str(scripts/script),*arguments],check=True)
    print('Assembled rates, missing-outcome bounds and diagnostic figures:',args.output,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['output','blended','extension']:
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--blended-score-file',default='blend_control_scores_20_draws.parquet')
    parser.add_argument('--draws-per-input',type=int,default=20)
    parser.add_argument('--legacy-draws-per-input',type=int,default=0)
    parser.add_argument('--legacy-batch',type=Path)
    parser.add_argument('--bootstrap',type=int,default=1000)
    assemble(parser.parse_args())
