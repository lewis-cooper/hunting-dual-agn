#!/usr/bin/env python3
"""Pair measured field sources with fixed real input primaries, preserving provenance."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from control_catalogues import read_csv,pm_significance,magnitude
from score_controls import HEMISPHERE,configure_production
from score_controls import cached_scores


def source_library(source_root,dustmaps):
    from dustmaps.sfd import SFDQuery
    matched=pd.read_parquet(source_root/'matched_neighbours.parquet')
    gaia=pd.read_parquet(source_root/'gaia_neighbours.parquet')
    primary=pd.read_parquet(source_root/'parent_primary_photometry.parquet').set_index('target_id')
    associations=matched.loc[matched.has_gaia].drop_duplicates(['target_id','gaia_source_id'])
    measurements=associations.drop(columns=['separation_arcsec','pm_sig','gaia_g'])
    gaia_pool=gaia.loc[gaia.separation_arcsec.gt(3)&gaia.separation_arcsec.le(30)].copy()
    gaia_pool=gaia_pool.rename(columns={'ra':'gaia_ra','dec':'gaia_dec'})
    gaia_pool=gaia_pool.merge(measurements,left_on=['target_id','source_id'],
        right_on=['target_id','gaia_source_id'],how='left',validate='one_to_one')
    gaia_pool['ra']=gaia_pool.ra.fillna(gaia_pool.gaia_ra)
    gaia_pool['dec']=gaia_pool.dec.fillna(gaia_pool.gaia_dec)
    gaia_pool['control_id']='gaia_'+gaia_pool.source_id.astype(str)
    gaia_pool['gaia_source_id']=gaia_pool.source_id
    gaia_pool['gaia_g']=gaia_pool.phot_g_mean_mag
    gaia_pool['pm_sig']=pm_significance(gaia_pool)
    gaia_pool['has_gaia']=True;gaia_pool['type']=gaia_pool.type.fillna('UNKNOWN')
    records=[]
    for route in ['resolved','blended','imaging_two_gaia','imaging_one_gaia']:
        pool=gaia_pool.copy() if route in ['resolved','blended'] else matched.loc[
            matched.separation_arcsec.gt(3)&matched.separation_arcsec.le(30)&matched.type.eq('PSF')].copy()
        if route in ['resolved','imaging_two_gaia','imaging_one_gaia']:
            pool=pool.loc[pool.target_id.map(primary.primary_unique).eq(True)]
        if route.startswith('imaging'):
            pool=pool.loc[pool.target_id.map(primary.primary_type).eq('PSF')]
        if route=='imaging_two_gaia':
            pool=pool.loc[pool.has_gaia]
        pool['route']=route
        pool['annulus']=np.where(pool.separation_arcsec.le(10),'3-10','10-30')
        pool['companion_mag']=pool.gaia_g if route in ['resolved','blended'] else pool.r_mag
        pool['background_source_id']=pool.control_id
        pool['background_field_id']=pool.target_id
        counts=pool.groupby(['target_id','annulus']).control_id.transform('count')
        pool['source_weight']=pool.target_id.map(primary.sampling_weight)/counts
        # Combine repeated observations of the same physical source; retain a single bootstrap identity.
        exposure=pool.groupby(['annulus','background_source_id']).source_weight.sum()
        pool=pool.drop_duplicates(['annulus','background_source_id']).copy()
        pool['source_weight']=[exposure.loc[(annulus,identity)] for annulus,identity
                               in zip(pool.annulus,pool.background_source_id)]
        records.append(pool)
    library=pd.concat(records,ignore_index=True)
    coordinates=SkyCoord(library.ra.to_numpy(float)*u.deg,library.dec.to_numpy(float)*u.deg)
    library['galactic_b_deg']=coordinates.galactic.b.deg
    north=library.dec.ge(32.375)&library.galactic_b_deg.gt(0)
    library['native_region']=np.where(north,'north','south')
    library['ebv']=SFDQuery(map_dir=str(dustmaps/'sfd'))(coordinates)
    library['mag_cell']=np.floor(library.companion_mag*2)/2
    return library


def real_primary_scores(inputs,args):
    configure_production(args)
    from qso_pcolor.unified import UnifiedPSFModel
    from dual_agn.score_selection import read_either_model_threshold
    models={name:UnifiedPSFModel.load(args.models/f'20261005_{name}_qlf') for name in ['magdep','magindep']}
    threshold=read_either_model_threshold(
        getattr(args,'calibration',None) or args.repo/'products/threshold_calibration.json')
    primaries=inputs.loc[inputs.route.str.startswith('imaging')&inputs.latitude_bin.notna()].copy()
    primaries['trial_id']=primaries.pair_id+':fixed_primary'
    primaries['target_id']=primaries.pair_id
    primaries['control_id']=primaries.pair_id+':primary'
    primaries['type']=primaries.primary_type
    primaries['photometric_system']=primaries.primary_release.map(HEMISPHERE)
    primaries['fracflux_r']=np.nan
    primaries['separation_arcsec']=primaries.pair_separation_arcsec
    primaries['maskbits']=0
    for band in 'grz':
        primaries[f'flux_{band}']=primaries[f'primary_flux_{band}']
        primaries[f'flux_ivar_{band}']=primaries[f'primary_ivar_{band}']
        primaries[f'nobs_{band}']=1
    primaries['secondary_pm_pass']=np.isfinite(primaries.primary_pm_sig)&primaries.primary_pm_sig.lt(3)
    root=args.output/'score_cache/fixed_primary';root.mkdir(parents=True,exist_ok=True)
    scored=cached_scores(primaries,models,threshold,'imaging',root)
    scored.to_parquet(args.output/'fixed_primary_scores.parquet',index=False,compression='zstd')
    return scored


def select_sources(template,library,max_distance_deg,latitude_tolerance_deg,ebv_tolerance):
    if not np.isfinite(template.companion_mag):
        return library.iloc[:0]
    cell=np.floor(template.companion_mag*2)/2
    source=library.loc[library.mag_cell.eq(cell)&library.native_region.eq(template.native_region)
        &(library.galactic_b_deg*template.galactic_b_deg).gt(0)
        &(library.galactic_b_deg.abs()-abs(template.galactic_b_deg)).abs().le(latitude_tolerance_deg)
        &(library.ebv-template.ebv).abs().le(ebv_tolerance)].copy()
    if not len(source):
        return source
    distance=SkyCoord(template.ra*u.deg,template.dec*u.deg).separation(
        SkyCoord(source.ra.to_numpy(float)*u.deg,source.dec.to_numpy(float)*u.deg)).deg
    source['sky_match_distance_deg']=distance
    return source.loc[(distance>1/60)&(distance<=max_distance_deg)]


def prepare(args):
    root=args.output;root.mkdir(parents=True,exist_ok=True)
    inputs=read_csv(args.source/'input_pairs.csv')
    inputs.to_csv(root/'input_pairs.csv',index=False)
    shutil.copy2(args.source/'pair_cut_factors.csv',root/'pair_cut_factors.csv')
    primary=real_primary_scores(inputs,args).set_index('pair_id')
    from dustmaps.sfd import SFDQuery
    primary_any=primary[['score_pass_magdep','score_pass_magindep']].any(axis=1)
    inputs=inputs.loc[inputs.latitude_bin.notna()].copy()
    coordinates=SkyCoord(inputs.ra.to_numpy(float)*u.deg,inputs.dec.to_numpy(float)*u.deg)
    inputs['galactic_b_deg']=coordinates.galactic.b.deg
    inputs['native_region']=np.where(inputs.dec.ge(32.375)&inputs.galactic_b_deg.gt(0),'north','south')
    inputs['ebv']=SFDQuery(map_dir=str(args.dustmaps/'sfd'))(coordinates)
    library=source_library(args.source,args.dustmaps)
    library.to_parquet(root/'background_source_library.parquet',index=False,compression='zstd')
    rng=np.random.default_rng(args.seed)
    records=[];matching=[]
    for route in (args.routes or ['resolved','imaging_two_gaia','imaging_one_gaia']):
        real=inputs.loc[inputs.route.eq(route)]
        for annulus in ['3-10','10-30']:
            pool=library.loc[library.route.eq(route)&library.annulus.eq(annulus)]
            for number,template in enumerate(real.itertuples(),1):
                known_rejection=route.startswith('imaging') and not bool(primary_any.get(template.pair_id,False))
                sources=pool.iloc[:0] if known_rejection else select_sources(template,pool,
                    args.max_sky_distance,args.latitude_tolerance,args.ebv_tolerance)
                n_available=len(sources)
                if n_available and n_available<=args.draws_per_input:
                    chosen=sources.copy()
                    chosen['mock_weight']=chosen.source_weight/chosen.source_weight.sum()
                    chosen['draw_multiplicity']=1
                    mode='enumerated'
                elif n_available:
                    probabilities=sources.source_weight.to_numpy(float)
                    draw=rng.choice(n_available,args.draws_per_input,replace=True,p=probabilities/probabilities.sum())
                    indices,counts=np.unique(draw,return_counts=True)
                    chosen=sources.iloc[indices].copy()
                    chosen['mock_weight']=counts/args.draws_per_input
                    chosen['draw_multiplicity']=counts
                    mode='sampled'
                else:
                    chosen=pd.DataFrame([dict(control_id='',background_source_id='',background_field_id='fixed',
                        ra=template.ra,dec=template.dec,type='UNKNOWN',photometric_system=np.nan,
                        flux_g=np.nan,flux_r=np.nan,flux_z=np.nan,flux_ivar_g=0.,flux_ivar_r=0.,flux_ivar_z=0.,
                        nobs_g=0,nobs_r=0,nobs_z=0,maskbits=0,fracflux_r=np.nan,pm_sig=np.nan,
                        separation_arcsec=np.nan,mock_weight=1.,draw_multiplicity=0)])
                    mode='known_primary_rejection' if known_rejection else 'unmatched'
                chosen['input_pair_id']=template.pair_id
                chosen['target_id']=template.pair_id
                chosen['route']=route;chosen['annulus']=annulus
                chosen['primary_ra']=template.ra;chosen['primary_dec']=template.dec
                chosen['redshift']=template.redshift;chosen['redshift_error']=template.redshift_error
                chosen['redshift_catalogue']=template.redshift_catalogue
                chosen['latitude_bin']=template.latitude_bin
                chosen['primary_gaia_g']=template.primary_gaia_g
                chosen['companion_mag_target']=template.companion_mag
                chosen['known_primary_rejection']=known_rejection
                chosen['diagnostic_only']=mode=='unmatched'
                chosen['selection_mode']=mode
                chosen['trial_id']=route+':'+annulus+':'+template.pair_id+':'+chosen.control_id
                chosen['model_key']=route+':'+template.pair_id+':'+chosen.control_id
                chosen['secondary_pm_pass']=True if route=='imaging_one_gaia' else (~np.isfinite(chosen.pm_sig)|chosen.pm_sig.lt(3))
                chosen['secondary_pm_missing']=False if route=='imaging_one_gaia' else ~np.isfinite(chosen.pm_sig)
                chosen['preselection_pass']=not known_rejection and mode!='unmatched'
                chosen['primary_pm_sig']=template.primary_pm_sig
                for band in 'grz':
                    chosen[f'primary_flux_{band}']=getattr(template,f'primary_flux_{band}')
                    chosen[f'primary_ivar_{band}']=getattr(template,f'primary_ivar_{band}')
                records.append(chosen)
                matching.append(dict(input_pair_id=template.pair_id,route=route,annulus=annulus,
                    n_available_sources=n_available,n_unique_drawn=len(chosen) if n_available else 0,
                    known_primary_rejection=known_rejection,selection_mode=mode))
                if number%1000==0:
                    print(f'Re-pairing {route} {annulus}: {number}/{len(real)} inputs',flush=True)
    trials=pd.concat(records,ignore_index=True)
    if not trials.trial_id.is_unique:
        raise ValueError('Re-paired trial identities are not unique')
    np.testing.assert_allclose(trials.groupby(['input_pair_id','annulus']).mock_weight.sum(),1.)
    trials.to_parquet(root/'component_control_trials.parquet',index=False,compression='zstd')
    pd.DataFrame(matching).to_csv(root/'component_matching_audit.csv',index=False)
    (root/'repair_configuration.json').write_text(json.dumps(dict(source=str(args.source.resolve()),
        seed=args.seed,draws_per_input=args.draws_per_input,max_sky_distance_deg=args.max_sky_distance,
        latitude_tolerance_deg=args.latitude_tolerance,ebv_tolerance=args.ebv_tolerance,
        matching='same 0.5-mag companion bin, Galactic hemisphere, native Legacy region; nearby sky and similar reddening',
        primary='actual input redshift/errors, photometry and PM held fixed for component routes',
        variance='reused background sources retain source/field IDs; combinations are not independent data',
        known_failures='fixed primary rejection retained without requiring a matching background source'),indent=2)+'\n')
    print(trials.groupby(['route','annulus']).size().to_string(),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--models',type=Path,required=True)
    parser.add_argument('--dustmaps',type=Path,required=True)
    parser.add_argument('--calibration',type=Path)
    parser.add_argument('--draws-per-input',type=int,default=100)
    parser.add_argument('--max-sky-distance',type=float,default=40.)
    parser.add_argument('--latitude-tolerance',type=float,default=10.)
    parser.add_argument('--ebv-tolerance',type=float,default=.02)
    parser.add_argument('--seed',type=int,default=20261121)
    parser.add_argument('--routes',nargs='+')
    prepare(parser.parse_args())


if __name__=='__main__':
    main()
