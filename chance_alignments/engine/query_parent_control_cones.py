#!/usr/bin/env python3
"""Fetch public control cones with per-field checkpoints and bounded memory."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from control_catalogues import cone_sql,public_query,read_csv


def field_path(root,survey,target_id):
    return root/'field_cones'/survey/(target_id+'.parquet')


def normalize_columns(frame,survey):
    frame=frame.copy()
    strings={'target_id','source_id','type'}
    integers={'release','brickid','objid','nobs_g','nobs_r','nobs_z','maskbits','astrometric_params_solved'}
    for column in frame:
        if column in strings:
            frame[column]=frame[column].astype('string')
        elif column in integers:
            frame[column]=pd.to_numeric(frame[column],errors='raise').astype('Int64')
        else:
            frame[column]=pd.to_numeric(frame[column],errors='raise').astype(float)
    return frame


def combine_fields(root,survey,targets):
    destination=root/(survey+'_neighbours.parquet')
    writer=None;blocks=[];n_rows=0
    try:
        for number,target in enumerate(targets.itertuples(),1):
            blocks.append(pd.read_parquet(field_path(root,survey,target.target_id)))
            if len(blocks)<128 and number<len(targets):
                continue
            frame=pd.concat(blocks,ignore_index=True)
            table=pa.Table.from_pandas(frame,preserve_index=False)
            if writer is None:
                writer=pq.ParquetWriter(destination,table.schema,compression='zstd')
            writer.write_table(table)
            n_rows+=len(frame);blocks=[]
    finally:
        if writer is not None:
            writer.close()
    return n_rows


def query_controls(args):
    root=args.output
    targets=read_csv(root/'master_targets.csv')
    if not targets.target_id.is_unique:
        raise ValueError('Target IDs must be unique')
    if args.reuse_cones:
        previous=read_csv(args.reuse_cones/'master_targets.csv').set_index('target_id')
        shared=targets.loc[targets.target_id.isin(previous.index)].set_index('target_id')
        for coordinate in ['ra','dec']:
            if not (shared[coordinate]-previous.loc[shared.index,coordinate]).abs().le(1e-10).all():
                raise ValueError('Reused cone identity has different sky coordinates')
        for survey in ['tractor','gaia']:
            (root/'field_cones'/survey).mkdir(parents=True,exist_ok=True)
            for identity in shared.index:
                source=field_path(args.reuse_cones,survey,identity)
                destination=field_path(root,survey,identity)
                if source.exists() and not destination.exists():
                    os.link(source,destination)
    jobs=[]
    for survey in ['tractor','gaia']:
        (root/'field_cones'/survey).mkdir(parents=True,exist_ok=True)
        missing=targets.loc[~targets.target_id.map(lambda target:field_path(root,survey,target).exists())]
        jobs.extend((survey,missing.iloc[start:start+args.batch_size])
                    for start in range(0,len(missing),args.batch_size))
    print('Parent fields:',len(targets),'| Missing query batches:',len(jobs),'| Workers:',args.workers,flush=True)
    def fetch(job):
        survey,group=job
        sql=cone_sql(group,survey,point_sources_only=True,gaia_quality=True)
        frame=normalize_columns(public_query(sql,root/'query_cache'),survey)
        if not set(frame.target_id).issubset(set(group.target_id)):
            raise ValueError('Query returned an unexpected target identity')
        for target in group.itertuples():
            field=frame.loc[frame.target_id.eq(target.target_id)]
            destination=field_path(root,survey,target.target_id)
            temporary=destination.with_suffix('.tmp')
            field.to_parquet(temporary,index=False,compression='zstd')
            temporary.replace(destination)
        return survey,len(frame)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for number,(survey,n_rows) in enumerate(pool.map(fetch,jobs),1):
            if number%25==0 or number==len(jobs):
                print(f'Completed query batches {number}/{len(jobs)}; latest {survey}: {n_rows} rows',flush=True)
    totals={survey:combine_fields(root,survey,targets) for survey in ['tractor','gaia']}
    targets[['target_id','ra','dec']].assign(tractor_query_complete=True,gaia_query_complete=True).to_csv(
        root/'cone_coverage.csv',index=False)
    (root/'query_configuration.json').write_text(json.dumps(dict(
        n_parent_fields=len(targets),neighbour_rows=totals,radius_arcsec=60,
        tractor_pool='brick_primary; all detections within 3 arcsec, PSF detections outside',
        gaia_pool='all Gaia DR3 detections; includes RUWE and solution type',
        checkpoint='one file per field and survey, including successful empty cones',
        parent_targets_sha256=hashlib.sha256((root/'master_targets.csv').read_bytes()).hexdigest()),indent=2)+'\n')
    print('Completed public cones:',totals,'| Output:',root,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--batch-size',type=int,default=8)
    parser.add_argument('--reuse-cones',type=Path)
    query_controls(parser.parse_args())


if __name__=='__main__':
    main()
