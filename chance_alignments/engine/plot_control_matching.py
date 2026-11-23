#!/usr/bin/env python3
"""Check companion-magnitude matching before applying the control vetting cuts."""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

from control_catalogues import ROUTES,read_csv


def plot(args):
    os.environ.setdefault('MPLCONFIGDIR',str(args.output/'scorer_runtime/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    inputs=read_csv(args.output/'input_pairs.csv').set_index('pair_id')
    components=pd.read_parquet(args.output/'component_control_scores.parquet')
    blends=pd.read_parquet(args.blended_trials)
    trials=pd.concat([components,blends],ignore_index=True)
    trials=trials.loc[trials.annulus.eq('3-10')&trials.background_source_id.fillna('').ne('')].copy()
    trials['target_companion_mag']=trials.input_pair_id.map(inputs.companion_mag)
    figure,axes=plt.subplots(2,2,figsize=(10,8),layout='constrained')
    results=[]
    for axis,route in zip(axes.flat,ROUTES):
        rows=trials.loc[trials.route.eq(route)].copy()
        valid=np.isfinite(rows.companion_mag)&np.isfinite(rows.target_companion_mag)&rows.mock_weight.gt(0)
        rows=rows.loc[valid]
        if not len(rows):
            continue
        difference=rows.companion_mag-rows.target_companion_mag
        bins=np.arange(np.floor(min(rows.companion_mag.min(),rows.target_companion_mag.min())*2)/2,
                       np.ceil(max(rows.companion_mag.max(),rows.target_companion_mag.max())*2)/2+.25,.25)
        axis.hist(rows.target_companion_mag,bins=bins,weights=rows.mock_weight,density=True,
                  histtype='step',linewidth=1.8,color='#333333',label='Real companion targets')
        axis.hist(rows.companion_mag,bins=bins,weights=rows.mock_weight,density=True,
                  histtype='step',linewidth=1.6,color='#0072B2',label='Matched field sources')
        axis.set_title(route.replace('_',' ').replace('gaia','Gaia'))
        axis.set_xlabel('Gaia G [mag]' if route in ['resolved','blended'] else 'Native Legacy r [mag]')
        axis.set_ylabel('Weighted probability density');axis.grid(alpha=.2);axis.legend(fontsize=8)
        results.append(dict(route=route,n_represented_inputs=rows.input_pair_id.nunique(),
            n_distinct_background_sources=rows.background_source_id.nunique(),
            weighted_mean_magnitude_residual=float(np.average(difference,weights=rows.mock_weight)),
            maximum_absolute_magnitude_residual=float(difference.abs().max()),
            matching='same 0.5-mag bin; this diagnostic precedes all vetting cuts'))
    for extension in ['pdf','png']:
        figure.savefig(args.output/f'control_magnitude_matching.{extension}',dpi=180)
    plt.close(figure)
    pd.DataFrame(results).to_csv(args.output/'control_magnitude_matching.csv',index=False)
    print(pd.DataFrame(results).to_string(index=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--blended-trials',type=Path,required=True)
    plot(parser.parse_args())
