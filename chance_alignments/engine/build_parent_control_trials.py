#!/usr/bin/env python3
"""Build main and check-annulus trials around independently sampled quasars."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from control_catalogues import read_csv,pm_significance
from score_controls import match_catalogues,HEMISPHERE


def resolve_photometric_area(tractor):
    """Use the same DR9 North/South boundary as the research input builder."""
    galactic_latitude=SkyCoord(tractor.ra.to_numpy(float)*u.deg,tractor.dec.to_numpy(float)*u.deg).galactic.b.deg
    north=tractor.dec.ge(32.375).to_numpy()&(galactic_latitude>0)
    keep=(tractor.release.eq(9011).to_numpy()&north
          |tractor.release.isin([9010,9012]).to_numpy()&~north)
    return tractor.loc[keep].copy()


def assemble_trials(root,check_limit=10000,seed=20261121):
    targets=read_csv(root/'master_targets.csv')
    tractor=resolve_photometric_area(pd.read_parquet(root/'tractor_neighbours.parquet'))
    gaia=pd.read_parquet(root/'gaia_neighbours.parquet')
    matched=match_catalogues(tractor,gaia)
    matched['type']=matched.type.str.strip()
    matched.to_parquet(root/'matched_neighbours.parquet',index=False,compression='zstd')
    central=matched.loc[matched.separation_arcsec.le(.5)].copy()
    central_count=central.groupby('target_id').size()
    central=central.sort_values('separation_arcsec').drop_duplicates('target_id').set_index('target_id')
    primary=targets.copy()
    for field in ['control_id','type','release','photometric_system','ra','dec',
                  'flux_g','flux_r','flux_z','flux_ivar_g','flux_ivar_r','flux_ivar_z',
                  'nobs_g','nobs_r','nobs_z','maskbits']:
        primary['primary_'+field]=primary.target_id.map(central[field])
    primary['primary_unique']=primary.target_id.map(central_count).eq(1)
    central_gaia=gaia.merge(targets[['target_id','source_id']],on=['target_id','source_id'],validate='many_to_one')
    if central_gaia.target_id.duplicated().any():
        raise ValueError('Duplicate central Gaia detections')
    central_gaia['pm_sig']=pm_significance(central_gaia)
    primary['primary_pm_sig']=primary.target_id.map(central_gaia.set_index('target_id').pm_sig)
    inner_tractor=matched.loc[matched.separation_arcsec.le(3)].copy()
    inner_tractor=inner_tractor.merge(primary[['target_id','primary_control_id']],on='target_id',validate='many_to_one')
    extra_tractor=inner_tractor.loc[inner_tractor.control_id.ne(inner_tractor.primary_control_id)].groupby('target_id').size()
    inner_gaia=gaia.loc[gaia.separation_arcsec.le(3)].merge(targets[['target_id','source_id']],
        on='target_id',suffixes=('','_primary'),validate='many_to_one')
    extra_gaia=inner_gaia.loc[inner_gaia.source_id.ne(inner_gaia.source_id_primary)].groupby('target_id').size()
    primary['primary_isolated']=primary.primary_unique|primary.target_id.map(central_count).isna()
    primary['primary_isolated'] &= primary.target_id.map(extra_tractor).fillna(0).eq(0)
    primary['primary_isolated'] &= primary.target_id.map(extra_gaia).fillna(0).eq(0)
    primary.to_parquet(root/'parent_primary_photometry.parquet',index=False,compression='zstd')
    metadata=primary.drop(columns=['ra','dec']).rename(columns={'primary_ra':'primary_tractor_ra',
        'primary_dec':'primary_tractor_dec','source_id':'primary_gaia_source_id'})
    metadata=metadata.merge(targets[['target_id','ra','dec']].rename(columns={'ra':'primary_ra','dec':'primary_dec'}),
                            on='target_id',validate='one_to_one')
    # Gaia routes retain Gaia detections without Tractor measurements as explicit controls.
    associations=matched.loc[matched.has_gaia].copy()
    associations=associations.drop_duplicates(['target_id','gaia_source_id'])
    gaia_pool=gaia.loc[gaia.separation_arcsec.gt(3)&gaia.separation_arcsec.le(30)].copy()
    gaia_pool=gaia_pool.rename(columns={'ra':'gaia_ra','dec':'gaia_dec'})
    measurements=associations.drop(columns=['separation_arcsec','pm_sig','gaia_g'],errors='ignore')
    gaia_pool=gaia_pool.merge(measurements,left_on=['target_id','source_id'],
        right_on=['target_id','gaia_source_id'],how='left',validate='one_to_one')
    gaia_pool['ra']=gaia_pool.ra.fillna(gaia_pool.gaia_ra)
    gaia_pool['dec']=gaia_pool.dec.fillna(gaia_pool.gaia_dec)
    gaia_pool['control_id']='gaia_'+gaia_pool.source_id.astype(str)
    gaia_pool['gaia_source_id']=gaia_pool.source_id
    gaia_pool['gaia_g']=gaia_pool.phot_g_mean_mag
    gaia_pool['pm_sig']=pm_significance(gaia_pool)
    gaia_pool['has_gaia']=True
    gaia_pool['type']=gaia_pool.type.fillna('UNKNOWN')
    records=[]
    for route in ['resolved','blended','imaging_two_gaia','imaging_one_gaia']:
        if route in ['resolved','blended']:
            pool=gaia_pool.copy()
        else:
            pool=matched.loc[matched.separation_arcsec.gt(3)&matched.separation_arcsec.le(30)
                &matched.type.eq('PSF')].copy()
            if route=='imaging_two_gaia':
                pool=pool.loc[pool.has_gaia].copy()
        pool=pool.drop(columns=pool.columns.intersection(metadata.columns).difference(['target_id']))
        frame=pool.merge(metadata,on='target_id',validate='many_to_one')
        if route.startswith('imaging'):
            # These inputs have already entered the PSF-pair route; condition controls the same way.
            frame=frame.loc[frame.primary_unique&frame.primary_type.eq('PSF')].copy()
        elif route=='blended':
            frame=frame.loc[frame.primary_isolated].copy()
        frame['route']=route
        frame['photometry_class']='not_applicable'
        if route=='blended':
            has_legacy=np.isfinite(frame.primary_flux_r)&np.isfinite(frame.primary_flux_ivar_r)&frame.primary_flux_ivar_r.gt(0)
            frame['photometry_class']=np.where(has_legacy,'legacy','no_legacy')
        frame['companion_mag']=frame.gaia_g if route in ['resolved','blended'] else frame.r_mag
        frame['annulus']=np.where(frame.separation_arcsec.le(10),'3-10','10-30')
        frame['secondary_pm_missing']=(~np.isfinite(frame.pm_sig) if route!='imaging_one_gaia' else False)
        frame['secondary_pm_pass']=True if route=='imaging_one_gaia' else (~np.isfinite(frame.pm_sig)|frame.pm_sig.lt(3))
        frame['control_population']='matched_parent_donors'
        frame['trial_id']=route+':'+frame.annulus+':'+frame.target_id+':'+frame.control_id
        frame['n_controls_in_host_annulus']=frame.groupby(['target_id','annulus']).trial_id.transform('count')
        frame['diagnostic_only']=False
        frame['control_subsampling_weight']=1.
        check=frame.loc[frame.annulus.eq('10-30')]
        if check_limit and len(check)>check_limit:
            chosen=check.sample(check_limit,random_state=seed)
            frame=pd.concat([frame.loc[frame.annulus.eq('3-10')],chosen],ignore_index=True)
            frame.loc[frame.annulus.eq('10-30'),'control_subsampling_weight']=len(check)/check_limit
        records.append(frame)
    trials=pd.concat(records,ignore_index=True)
    if not trials.trial_id.is_unique:
        raise ValueError('Control trial identities are not unique')
    trials.to_parquet(root/'parent_control_trials.parquet',index=False,compression='zstd')
    counts=trials.groupby(['route','annulus']).size()
    (root/'trial_preparation.json').write_text(json.dumps(dict(n_parent_fields=len(targets),
        n_matched_tractor_rows=len(matched),n_isolated_primaries=int(primary.primary_isolated.sum()),
        check_limit_per_route=check_limit,seed=seed,
        primary_conditioning='PSF unique parent for imaging; isolated Gaia parent for blends, with a unique Tractor primary or no Tractor source, and no extra primary morphology cut',
        annuli={'main':[3,10],'check':[10,30]},
        trial_counts={f'{route}:{annulus}':int(value) for (route,annulus),value in counts.items()}),indent=2)+'\n')
    print(counts.to_string(),flush=True)
    return trials,primary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--check-limit',type=int,default=10000)
    parser.add_argument('--seed',type=int,default=20261121)
    args=parser.parse_args()
    assemble_trials(args.output,args.check_limit,args.seed)


if __name__=='__main__':
    main()
