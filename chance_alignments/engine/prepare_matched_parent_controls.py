#!/usr/bin/env python3
"""Draw expandable parent controls with known inclusion weights and stable IDs."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from control_catalogues import latitude_bin, read_csv, read_fits

PARENT_COLUMNS = ['qso_parent_index','source_id','gaia_ra','gaia_dec',
                  'phot_g_mean_mag','adopted_redshift','adopted_redshift_catalogue']


def catalogue_family(values):
    names=values.fillna('').astype(str).str.upper()
    return pd.Series(np.select([names.str.contains('QUAIA'),names.str.contains('DESI'),
                               names.str.contains('SDSS')],['QUAIA','DESI','SDSS'],default='OTHER'),index=values.index)


def primary_cells(frame):
    """Match source family, broad sky region, primary G and redshift explicitly."""
    coordinates=SkyCoord(frame.ra.to_numpy(float)*u.deg,frame.dec.to_numpy(float)*u.deg)
    longitude_quadrant=np.floor(coordinates.galactic.l.deg/90).astype(int)
    cells=pd.Series(pd.NA,index=frame.index,dtype='string')
    finite=np.isfinite(frame.primary_gaia_g)&np.isfinite(frame.redshift)&frame.latitude_bin.notna()
    cells.loc[finite]=(frame.loc[finite,'latitude_bin'].astype(str)+'|'
        +catalogue_family(frame.loc[finite,'redshift_catalogue'])+'|'
        +pd.Series(longitude_quadrant,index=frame.index).loc[finite].astype(str)+'|'
        +np.floor(frame.loc[finite,'primary_gaia_g']).astype(int).astype(str)+'|'
        +np.floor(frame.loc[finite,'redshift']/.5).astype(int).astype(str))
    return cells


def annotate_input_photometry(inputs,one_arcsecond_scores):
    """Blended mocks must represent the actual Legacy versus fallback access pattern."""
    scores=read_csv(one_arcsecond_scores).set_index('qso_group_id')
    inputs=inputs.copy()
    inputs['photometry_class']='not_applicable'
    mask=inputs.route.eq('blended')
    identifiers=inputs.loc[mask,'pair_id'].str.removeprefix('blended:')
    assignment=identifiers.map(scores.blend_decals_system_assignment)
    inputs.loc[mask,'photometry_class']=np.where(assignment.eq('none'),'no_legacy','legacy')
    if assignment.isna().any():
        raise ValueError('Missing native-photometry assignment for a blended input')
    return inputs


def parent_chunks(path):
    for chunk in pd.read_csv(path,usecols=PARENT_COLUMNS,dtype={'source_id':'string'},chunksize=100000):
        chunk=chunk.rename(columns={'gaia_ra':'ra','gaia_dec':'dec',
            'phot_g_mean_mag':'primary_gaia_g','adopted_redshift':'redshift',
            'adopted_redshift_catalogue':'redshift_catalogue'})
        finite=np.isfinite(chunk.ra)&np.isfinite(chunk.dec)
        chunk=chunk.loc[finite].copy()
        chunk['abs_b_deg'],chunk['latitude_bin']=latitude_bin(chunk.ra,chunk.dec)
        chunk=chunk.loc[chunk.latitude_bin.notna()].copy()
        chunk['primary_cell']=primary_cells(chunk)
        # Missing-primary cells stay in the density population, never masquerading as scored controls.
        chunk['sampling_cell']=chunk.primary_cell.fillna(chunk.latitude_bin.astype(str)+'|missing')
        yield chunk


def allocate_quotas(population,inputs,total,minimum=2,previous=None):
    """Mix a representative density sample with equal emphasis on each scoring route."""
    population=pd.Series(population,dtype='int64').sort_index()
    minimum_quota=np.minimum(population.to_numpy(),minimum)
    if previous is not None:
        minimum_quota=np.maximum(minimum_quota,previous.reindex(population.index,fill_value=0).to_numpy())
    target=pd.Series(0.,index=population.index)
    for _,group in inputs.loc[inputs.primary_cell.notna()].groupby('route'):
        histogram=group.primary_cell.value_counts(normalize=True)
        target=target.add(histogram.reindex(target.index,fill_value=0),fill_value=0)
    target=target/target.sum()
    probability=.5*population/population.sum()+.5*target
    if total<int(minimum_quota.sum()):
        raise ValueError(f'Sample must contain at least {minimum_quota.sum()} fields for the requested cell floor')
    total=min(total,int(population.sum()))
    quotas=pd.Series(minimum_quota,index=population.index,dtype='int64')
    remaining=total-int(quotas.sum())
    while remaining:
        available=population-quotas
        weights=probability.where(available.gt(0),0.)
        expected=remaining*weights/weights.sum()
        increment=np.minimum(np.floor(expected).astype(int),available)
        if increment.sum()==0:
            eligible=expected.loc[available.gt(0)].sort_values(ascending=False)
            increment.loc[eligible.index[:remaining]]=1
        quotas+=increment
        remaining=total-int(quotas.sum())
    return quotas


def prepare(args):
    args.output.mkdir(parents=True,exist_ok=True)
    inputs=read_csv(args.previous/'input_pairs.csv')
    inputs['primary_cell']=primary_cells(inputs)
    if args.one_arcsecond_scores:
        inputs=annotate_input_photometry(inputs,args.one_arcsecond_scores)
    inputs.to_csv(args.output/'input_pairs.csv',index=False)
    shutil.copy2(args.previous/'pair_cut_factors.csv',args.output/'pair_cut_factors.csv')
    population=Counter()
    for chunk in parent_chunks(args.master):
        population.update(chunk.sampling_cell)
    previous_path=args.output/'parent_sampling_cells.csv'
    if not previous_path.exists() and (args.previous/'parent_sampling_cells.csv').exists():
        previous_path=args.previous/'parent_sampling_cells.csv'
    previous=(pd.read_csv(previous_path,index_col='sampling_cell').sample_count if previous_path.exists() else None)
    quotas=allocate_quotas(population,inputs,args.n_parents,args.minimum_per_cell,previous)
    rng=np.random.default_rng(args.seed)
    reservoirs={}
    for chunk in parent_chunks(args.master):
        chunk['priority']=rng.random(len(chunk))
        for cell,group in chunk.groupby('sampling_cell',sort=False):
            retained=pd.concat([reservoirs.get(cell,pd.DataFrame()),group],ignore_index=True)
            reservoirs[cell]=retained.nsmallest(int(quotas[cell]),'priority')
    targets=pd.concat(reservoirs.values(),ignore_index=True).sort_values('qso_parent_index')
    sample_counts=targets.sampling_cell.value_counts()
    targets['sampling_weight']=targets.sampling_cell.map(population)/targets.sampling_cell.map(sample_counts)
    targets['target_id']='parent_'+targets.qso_parent_index.astype(str)
    targets['target_kind']='master';targets['route']='master'
    if not targets.target_id.is_unique:
        raise ValueError('Parent identities are not unique')
    reference=read_fits(args.quaia,['source_id','redshift_quaia','redshift_quaia_err']).set_index('source_id')
    uses_quaia=catalogue_family(targets.redshift_catalogue).eq('QUAIA')
    identifiers=pd.to_numeric(targets.source_id,errors='raise').astype('int64')
    targets['redshift_error']=identifiers.map(reference.redshift_quaia_err).where(uses_quaia,0.)
    reference_z=identifiers.map(reference.redshift_quaia)
    inconsistent=uses_quaia & (~np.isfinite(reference_z)
        | ~np.isclose(targets.redshift,reference_z,rtol=0,atol=1e-12)
        | ~np.isfinite(targets.redshift_error) | targets.redshift_error.lt(0))
    if inconsistent.any():
        raise ValueError(f'{inconsistent.sum()} Quaia controls have missing/inconsistent published redshifts/errors')
    targets.to_csv(args.output/'master_targets.csv',index=False)
    pd.DataFrame({'parent_count':pd.Series(population),'sample_count':quotas}).rename_axis('sampling_cell').to_csv(
        args.output/'parent_sampling_cells.csv')
    metadata=dict(master=str(args.master.resolve()),previous=str(args.previous.resolve()),
        quaia=str(args.quaia.resolve()),n_sampled_parents=len(targets),
        n_parent_population=sum(population.values()),seed=args.seed,
        parent_cells='latitude bin, redshift catalogue family, Galactic longitude quadrant, 1-mag primary G, 0.5 redshift',
        allocation='50% parent-distribution + 50% equal mixture of four input-route distributions; cell floor',
        weights='parent population in cell / number sampled in cell',
        reuse='Stable IDs; deterministic priorities make larger samples nested',
        scope='Independent controls; production selection and reviewed candidate list unchanged')
    (args.output/'parent_sampling.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('Sampled parents:',len(targets),'| Parent population:',sum(population.values()),
          '| Sampling cells:',len(quotas),'| Output:',args.output,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous',type=Path,required=True)
    parser.add_argument('--master',type=Path,required=True)
    parser.add_argument('--quaia',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--n-parents',type=int,default=30000)
    parser.add_argument('--minimum-per-cell',type=int,default=2)
    parser.add_argument('--seed',type=int,default=20261121)
    parser.add_argument('--one-arcsecond-scores',type=Path)
    prepare(parser.parse_args())


if __name__=='__main__':
    main()
