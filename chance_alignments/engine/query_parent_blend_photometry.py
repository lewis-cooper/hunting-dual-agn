#!/usr/bin/env python3
"""Fetch independent native fluxes; assess band quality after adding the components."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from native_photometry import additional_photometry


def query_photometry(args):
    root=args.output
    trials=pd.read_parquet(args.trials or root/'parent_control_trials.parquet')
    trials=trials.loc[trials.route.eq('blended')].copy()
    trials.to_csv(root/'isolated_blend_trials.csv',index=False)
    required=trials.loc[trials.secondary_pm_pass & ~trials.diagnostic_only.eq(True)].copy() if 'diagnostic_only' in trials else trials.loc[trials.secondary_pm_pass].copy()
    sources={};measurements=[];pair_sources=[];source_surveys={}
    for trial in required.itertuples():
        for component in ['primary','secondary']:
            source_id=trial.primary_gaia_source_id if component=='primary' else trial.gaia_source_id
            ra=(trial.primary_tractor_ra if np.isfinite(trial.primary_tractor_ra) else trial.primary_ra) if component=='primary' else trial.ra
            dec=(trial.primary_tractor_dec if np.isfinite(trial.primary_tractor_dec) else trial.primary_dec) if component=='primary' else trial.dec
            source_key=str(source_id)
            if trial.photometry_class=='no_legacy':
                sources[source_key]=(float(ra),float(dec))
                policy=getattr(trial,'band_policy','ps1,nsc,skymapper')
                source_surveys.setdefault(source_key,set()).update(str(policy).split(','))
                pair_sources.append(dict(trial_id=trial.trial_id,component=component,source_key=source_key))
                continue
            system=trial.primary_photometric_system if component=='primary' else trial.photometric_system
            if not isinstance(system,str):
                continue
            for band in 'grz':
                flux=getattr(trial,f'primary_flux_{band}' if component=='primary' else f'flux_{band}')
                ivar=getattr(trial,f'primary_flux_ivar_{band}' if component=='primary' else f'flux_ivar_{band}')
                if np.isfinite(flux+ivar) and ivar>0:
                    measurements.append(dict(trial_id=trial.trial_id,component=component,band=f'{system}:{band}',
                        flux=flux,variance=1/ivar,quality_flux=np.nan,quality_variance=np.nan,
                        independent=True,measurement='Tractor component'))
    cache=root/'independent_photometry_cache';cache.mkdir(exist_ok=True)
    def fetch(item):
        key,(ra,dec)=item
        digest=hashlib.sha256(f'{key}:{ra:.10f}:{dec:.10f}:raw_native_flux_v3'.encode()).hexdigest()
        existing=[cache/(digest+'.json')]
        if args.cache_source:
            existing.append(args.cache_source/'independent_photometry_cache'/(digest+'.json'))
        for path in existing:
            if not path.exists():
                continue
            try:
                return key,json.loads(path.read_text())
            except json.JSONDecodeError:
                pass
        surveys=sorted(source_surveys[key])
        digest=hashlib.sha256(f'{key}:{ra:.10f}:{dec:.10f}:raw_native_flux_v4:{surveys}'.encode()).hexdigest()
        path=cache/(digest+'.json')
        if path.exists():
            return key,json.loads(path.read_text())
        records,audit=additional_photometry(ra,dec,root,component_quality=False,surveys=surveys)
        result=dict(records=records,audit=audit,ra=ra,dec=dec,source_key=key,
                    quality_policy='retain finite native component measurements; apply quality to mock sum')
        temporary=path.with_suffix('.tmp')
        temporary.write_text(json.dumps(result)+'\n')
        temporary.replace(path)
        return key,result
    source_photometry={}
    print('Blend trials:',len(trials),'| PM-eligible flux sums:',len(required),
          '| Unique fallback-photometry sources:',len(sources),flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for number,(key,result) in enumerate(pool.map(fetch,sources.items()),1):
            source_photometry[key]=result
            if number%25==0 or number==len(sources):
                print(f'Independent multi-survey sources: {number}/{len(sources)}',flush=True)
    audit=[]
    for pair in pair_sources:
        result=source_photometry[pair['source_key']]
        measurements.extend(dict(trial_id=pair['trial_id'],component=pair['component'],independent=True,**record)
                            for record in result['records'])
        audit.extend(dict(trial_id=pair['trial_id'],component=pair['component'],**record) for record in result['audit'])
    columns=['trial_id','component','band','flux','variance','quality_flux','quality_variance','independent','measurement']
    photometry=pd.DataFrame(measurements).reindex(columns=columns)
    if photometry.duplicated(['trial_id','component','band']).any():
        raise ValueError('Duplicate independent component measurements')
    photometry.to_csv(root/'independent_blend_photometry.csv',index=False)
    pd.DataFrame(audit).to_csv(root/'isolated_blend_catalogue_matches.csv',index=False)
    (root/'blend_photometry_configuration.json').write_text(json.dumps(dict(
        n_trials=len(trials),n_known_pm_rejections=int((~trials.secondary_pm_pass).sum()),
        n_trials_requiring_photometry=len(required),n_unique_queried_sources=len(sources),
        n_measurements=len(photometry),
        band_quality='PS1/NSC/SkyMapper PSF quality assessed after summing fluxes and independent variances',
        survey_policy='Legacy grz only where Legacy r is assigned; otherwise PS1-first then NSC/SkyMapper, matching the frozen production inputs',
        missing_measurements='never filled with zero; unresolved sums flagged unknown when required native bands are missing',
        secondary_pm='known PM failures retained as rejections without unnecessary flux queries'),indent=2)+'\n')
    print('Blend trials:',len(trials),'| Need flux sums:',len(required),'| Unique sources:',len(sources),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--trials',type=Path)
    parser.add_argument('--cache-source',type=Path)
    query_photometry(parser.parse_args())


if __name__=='__main__':
    main()
