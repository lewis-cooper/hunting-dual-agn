#!/usr/bin/env python3
"""Measure weighted 10–60 arcsec densities, keeping close-pair completeness explicit."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import gamma

from control_catalogues import ROUTES, read_csv


def density_interval(fields, counts, geometry, draws, rng):
    weights = fields.sampling_weight.to_numpy(float)
    counts = np.asarray(counts, float)
    values = weights*counts*geometry
    point = values.sum()
    codes, cells = pd.factorize(fields.sampling_cell)
    population_mass = np.bincount(codes, weights=weights, minlength=len(cells))
    replicates = []
    for _ in range(draws):
        multiplier = rng.exponential(size=len(fields))
        normalization = np.bincount(codes, weights=multiplier*weights, minlength=len(cells))
        multiplier *= population_mass[codes]/normalization[codes]
        replicates.append(np.sum(multiplier*values))
    lower, upper = np.quantile(replicates, [.025,.975])
    weighted_events = np.sum(weights*counts)
    if weighted_events:
        effective_events = weighted_events**2/np.sum(weights**2*counts)
        scale = weighted_events/effective_events*geometry
        guard = gamma.ppf([.025,.975], effective_events+.5)*scale
        lower = min(lower, guard[0]); upper = max(upper, guard[1])
    else:
        lower = 0.
        upper = gamma.ppf(.975,.5)*weights.mean()*geometry
    return point,float(lower),float(upper)


def summarize(args):
    root = args.source
    fields = read_csv(root/'master_targets.csv').set_index('target_id')
    pairs = read_csv(root/'input_pairs.csv')
    gaia = pd.read_parquet(root/'gaia_neighbours.parquet')
    gaia = gaia.loc[gaia.separation_arcsec.gt(10) & gaia.separation_arcsec.le(60)].copy()
    gaia['magnitude'] = gaia.phot_g_mean_mag
    tractor = pd.read_parquet(root/'matched_neighbours.parquet')
    tractor = tractor.loc[tractor.separation_arcsec.gt(10) & tractor.separation_arcsec.le(60)
                          & tractor.type.str.strip().eq('PSF')].copy()
    tractor['magnitude'] = tractor.r_mag
    rng = np.random.default_rng(args.seed)
    result = []; density_cells = []
    for route in ROUTES:
        pool = tractor if route.startswith('imaging') else gaia
        if route == 'imaging_two_gaia':
            pool = pool.loc[pool.has_gaia]
        radius = 1.5 if route.startswith('imaging') else 1.
        geometry = radius**2/(60**2-10**2)
        route_magnitudes=pairs.loc[pairs.route.eq(route)&pairs.latitude_bin.notna(),'companion_mag']
        route_magnitudes=route_magnitudes[np.isfinite(route_magnitudes)]
        low=np.floor(route_magnitudes.min()*2)/2
        high=np.floor(route_magnitudes.max()*2)/2+.5
        for latitude in ['25-40','40-55','55-90']:
            real = pairs.loc[pairs.route.eq(route) & pairs.latitude_bin.eq(latitude)]
            magnitudes = real.companion_mag[np.isfinite(real.companion_mag)]
            if not len(magnitudes):
                continue
            local_fields = fields.loc[fields.latitude_bin.eq(latitude)]
            nearby = pool.loc[pool.target_id.isin(local_fields.index)
                              & pool.magnitude.ge(low) & pool.magnitude.lt(high)]
            counts = nearby.groupby('target_id').size().reindex(local_fields.index, fill_value=0)
            point,lower,upper = density_interval(local_fields,counts,geometry,args.bootstrap,rng)
            result.append(dict(route=route,latitude_bin=latitude,N_input=len(real),
                n_parent_fields=len(local_fields),parent_population_weight=local_fields.sampling_weight.sum(),
                n_outer_sources=len(nearby),magnitude_lower=low,magnitude_upper=high,
                radius_arcsec=radius,expected_chance_pairs_C1=point,C1_lower_95=lower,C1_upper_95=upper))
            for edge in np.arange(low,high,.5):
                cell = nearby.loc[nearby.magnitude.ge(edge) & nearby.magnitude.lt(edge+.5)]
                expected = cell.target_id.map(local_fields.sampling_weight).sum()*geometry
                density_cells.append(dict(route=route,latitude_bin=latitude,magnitude_lower=edge,
                    magnitude_upper=edge+.5,n_outer_sources=len(cell),expected_chance_pairs_C1=expected))
    args.output.mkdir(parents=True,exist_ok=True)
    frame = pd.DataFrame(result)
    frame.to_csv(args.output/'parent_densities_by_latitude.csv',index=False)
    frame.groupby('route')[['N_input','n_parent_fields','n_outer_sources','expected_chance_pairs_C1']].sum().to_csv(
        args.output/'parent_densities_by_route.csv')
    pd.DataFrame(density_cells).to_csv(args.output/'parent_densities_by_magnitude.csv',index=False)
    (args.output/'parent_density_configuration.json').write_text(json.dumps(dict(
        annulus_arcsec=[10,60],n_parent_fields=len(fields),seed=args.seed,
        normalization='known parent sampling weights; Gaia-linked SDSS+DESI+Quaia(G<20) master at |b|>=25 degrees; no global G<20 cut on SDSS/DESI',
        uncertainty='stratified field bootstrap and finite-event Poisson guard; not completeness uncertainty',
        inference='C=1 geometric scenarios only; route detection completeness and full search eligibility unvalidated',
        prohibition='do not multiply pooled densities and pooled pass rates or report these as measured fake counts'),indent=2)+'\n')
    print(frame.groupby('route')[['n_outer_sources','expected_chance_pairs_C1']].sum().to_string(),flush=True)
    plot_densities(frame,args.output)


def plot_densities(frame,root):
    os.environ.setdefault('MPLCONFIGDIR',str(root/'scorer_runtime/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure,axes=plt.subplots(2,2,figsize=(10,8),layout='constrained')
    for axis,route in zip(axes.flat,ROUTES):
        rows=frame.loc[frame.route.eq(route)]
        scale=1000/rows.parent_population_weight.to_numpy(float)
        value=rows.expected_chance_pairs_C1.to_numpy(float)*scale
        errors=[np.maximum(0,value-rows.C1_lower_95.to_numpy(float)*scale),
                np.maximum(0,rows.C1_upper_95.to_numpy(float)*scale-value)]
        axis.errorbar(np.arange(3),value,yerr=errors,fmt='o-',capsize=3,color='#0072B2')
        axis.set_xticks([0,1,2],['25–40','40–55','55–90'])
        axis.set_title(route.replace('_',' ').replace('gaia','Gaia')+f' (radius {rows.radius_arcsec.iloc[0]:g} arcsec)')
        axis.set_xlabel(r'$|b|$ [degrees]')
        axis.set_ylabel('Geometric neighbours per 1000 quasars (C=1)')
        axis.set_ylim(bottom=0);axis.grid(alpha=.2)
    for extension in ['pdf','png']:
        figure.savefig(root/f'local_density_scenarios.{extension}',dpi=180)
    plt.close(figure)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--bootstrap',type=int,default=1000)
    parser.add_argument('--seed',type=int,default=20261121)
    summarize(parser.parse_args())


if __name__=='__main__':
    main()
