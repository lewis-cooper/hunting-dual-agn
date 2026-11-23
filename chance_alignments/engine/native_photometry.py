#!/usr/bin/env python3
"""Build independent multi-survey flux sums around isolated master quasars."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import requests

from control_catalogues import read_csv,read_fits,public_query,pm_significance,validate_cone_coverage
from score_controls import match_catalogues,HEMISPHERE


def flux_measurement(mag,error):
    if not (np.isfinite(mag+error) and -50<mag<90 and 0<error<90):return None
    flux=10**((22.5-mag)/2.5)
    return flux,(np.log(10)/2.5*flux*error)**2


def fetch_ps1(ra,dec,cache,include_fluxes=False):
    # Stack Kron totals are the same measurement family used in production blends.
    columns=['objID','raMean','decMean','primaryDetection']
    columns += [f'{b}{kind}' for b in 'grizy' for kind in ['PSFMag','PSFMagErr','KronMag','KronMagErr']]
    if include_fluxes:
        columns += [f'{b}{kind}' for b in 'grizy' for kind in ['PSFFlux','PSFFluxErr','KronFlux','KronFluxErr']]
    params={'ra':ra,'dec':dec,'radius':1/3600,'columns':'['+','.join(columns)+']','pagesize':100}
    key=hashlib.sha256(json.dumps(params,sort_keys=True).encode()).hexdigest()
    path=cache/(key+'.csv')
    cache.mkdir(parents=True,exist_ok=True)
    if path.exists():return pd.read_csv(path)
    (cache/(key+'.json')).write_text(json.dumps(params,indent=2)+'\n')
    url='https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/stack.csv'
    for attempt in range(4):
        try:
            response=requests.get(url,params=params,timeout=(15,60))
            response.raise_for_status()
            frame=pd.read_csv(io.StringIO(response.text)) if response.text.strip() else pd.DataFrame(columns=columns)
            if len(frame)>=100:raise ValueError('PS1 response may be truncated')
            frame.to_csv(path,index=False)
            return frame
        except requests.RequestException:
            if attempt==3:raise
            time.sleep(3*(attempt+1))


def additional_photometry(ra,dec,outdir,component_quality=True,surveys=None):
    records=[];audit=[]
    for survey in (surveys or ['ps1','nsc','skymapper']):
        if survey=='ps1':
            frame=fetch_ps1(ra,dec,outdir/'ps1_query_cache',include_fluxes=not component_quality)
            if len(frame):frame=frame.loc[frame.primaryDetection.eq(1)]
            bands='grizy'
        else:
            if survey=='nsc':
                table='nsc_dr2.object';racol='ra';deccol='dec';bands=['u','g','r','i','z','y','vr']
                columns=['id','ra','dec']+[f'{b}{kind}' for b in bands for kind in ['mag','err']]
            else:
                table='skymapper_dr4.master';racol='raj2000';deccol='dej2000';bands='uvgriz'
                columns=['object_id',racol,deccol]+[c for b in bands for c in [f'{b}_psf',f'e_{b}_psf',f'{b}_petro',f'e_{b}_petro']]
            sql=(f"SELECT {','.join(columns)} FROM {table} WHERE q3c_radial_query("
                 f'{racol},{deccol},{ra:.10f},{dec:.10f},{1/3600:.12f}) LIMIT 3')
            frame=public_query(sql,outdir/'query_cache')
        audit.append(dict(survey=survey,ra=ra,dec=dec,n_matches=len(frame),used=len(frame)==1))
        # Ambiguous matches do not silently acquire another object's flux.
        if len(frame)!=1:continue
        row=frame.iloc[0]
        for band in bands:
            if survey=='ps1':
                quality_mag,quality_error=row[f'{band}PSFMag'],row[f'{band}PSFMagErr']
                mag,error=row[f'{band}KronMag'],row[f'{band}KronMagErr']
            elif survey=='nsc':
                mag,error=row[f'{band}mag'],row[f'{band}err']
                quality_mag,quality_error=mag,error
            else:
                quality_mag,quality_error=row[f'{band}_psf'],row[f'e_{band}_psf']
                mag,error=row[f'{band}_petro'],row[f'e_{band}_petro']
            # Reconstruct the saved per-band eligibility flags, then use total flux.
            quality=np.isfinite(quality_mag+quality_error) and 0<quality_mag<40 and 0<=quality_error<=.3
            if survey=='ps1' and not component_quality and f'{band}KronFlux' in row.index:
                # MAST StackObjectView fluxes are Jy; keep negative noise measurements, excluding -999 sentinels.
                raw=np.array([row[f'{band}{name}'] for name in ['KronFlux','KronFluxErr','PSFFlux','PSFFluxErr']],float)
                measured=np.isfinite(raw).all() and (raw!=-999).all() and raw[1]>0 and raw[3]>0
                if measured:
                    converted=raw/(3631e-9)
                    records.append(dict(band=f'ps1:{band}',flux=converted[0],variance=converted[1]**2,
                        quality_flux=converted[2],quality_variance=converted[3]**2,
                        component_band_quality=bool(quality),measurement='PS1 stack Kron flux; Jy / (3631e-9)'))
                continue
            quality_measured=np.isfinite(quality_mag+quality_error) and 0<quality_mag<40 and 0<=quality_error<90
            values=flux_measurement(mag,error) if (quality if component_quality else quality_measured) else None
            if values:records.append(dict(band=f'{survey}:{band}',flux=values[0],variance=values[1],
                                           quality_flux=10**((22.5-quality_mag)/2.5) if quality_measured else np.nan,
                                           quality_variance=(np.log(10)/2.5*10**((22.5-quality_mag)/2.5)*quality_error)**2 if quality_measured else np.nan,
                                           component_band_quality=bool(quality),
                                           measurement='Kron' if survey=='ps1' else 'Petrosian' if survey=='skymapper' else 'NSC native'))
    return records,audit






