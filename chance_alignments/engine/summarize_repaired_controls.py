#!/usr/bin/env python3
"""Summarize re-paired controls without counting reused objects as independent."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

from control_catalogues import ROUTES, read_csv


def estimate_repaired_rate(sample, n_input, draws, rng):
    if not len(sample) or not n_input:
        return dict(f_pass_known=np.nan, f_pass_population_lower=0., f_pass_population_upper=1.,
                    f_pass_lower_95=np.nan, f_pass_upper_95=np.nan, relative_halfwidth_95=np.inf,
                    N_mocks=0, N_background_sources=0, N_background_fields=0, N_primary_fields=0,
                    known_weight_fraction=0., matching_coverage=0., effective_fields=0.)
    sample = sample.copy()
    weight = sample.mock_weight.to_numpy(float)
    factor = sample.pair_factor.to_numpy(float)
    unknown = ~sample.diagnostic_only.eq(False).to_numpy()
    passed = sample.individual_cuts_pass.eq(True).to_numpy(float)*~unknown
    known = ~unknown
    mass = weight.sum()
    known_mass = np.sum(weight*known)
    probability = np.sum(weight*passed*factor)/known_mass if known_mass else np.nan
    lower_bound = np.sum(weight*passed*factor)/mass
    upper_bound = lower_bound + np.sum(weight*unknown*factor)/mass
    input_codes, input_ids = pd.factorize(sample.input_pair_id)
    background_codes, background_fields = pd.factorize(sample.background_field_id.fillna('unmatched'))
    primary_codes = None
    if 'primary_field_id' in sample and sample.primary_field_id.notna().any():
        primary_codes, primary_fields = pd.factorize(sample.primary_field_id.fillna('fixed'))
        shared_fields=pd.Index(background_fields).union(pd.Index(primary_fields),sort=False)
        background_shared=shared_fields.get_indexer(background_fields)[background_codes]
        primary_shared=shared_fields.get_indexer(primary_fields)[primary_codes]
    else:
        primary_fields = []
    stochastic = sample.background_source_id.fillna('').ne('').to_numpy()
    # Weight carried by one queried sky field is one information cluster, even
    # when its sources appear beside many input quasars or primary donors.
    field_mass = np.bincount(background_codes, weights=weight*known*stochastic,
                             minlength=len(background_fields))
    effective = field_mass.sum()**2/np.sum(field_mass**2) if np.sum(field_mass**2) else 0.
    if primary_codes is not None:
        primary_mass = np.bincount(primary_codes, weights=weight*known*stochastic,
                                  minlength=len(primary_fields))
        primary_effective = primary_mass.sum()**2/np.sum(primary_mass**2) if np.sum(primary_mass**2) else 0.
        effective = min(effective, primary_effective)
    bootstrap = []
    sampled=sample.selection_mode.eq('sampled').to_numpy()
    multiplicity=sample.draw_multiplicity.fillna(1).to_numpy(float) if 'draw_multiplicity' in sample else np.ones(len(sample))
    multiplicity=np.where(sampled,np.maximum(1,multiplicity),1.)
    for _ in range(draws):
        if primary_codes is not None:
            field_multipliers=rng.exponential(size=len(shared_fields))
            multipliers=field_multipliers[background_shared]*field_multipliers[primary_shared]
        else:
            multipliers = rng.exponential(size=len(background_fields))[background_codes]
        multipliers = np.where(stochastic, multipliers, 1.)
        resampled = weight*multipliers
        if sampled.any():
            # Finite mock draws add Monte Carlo error on top of field sampling.
            jitter=rng.poisson(multiplicity)/multiplicity
            resampled*=np.where(sampled,jitter,1.)
        input_mass = np.bincount(input_codes, weights=resampled, minlength=len(input_ids))
        zero=input_mass[input_codes]<=0
        if zero.any():
            resampled[zero]=weight[zero]
            input_mass=np.bincount(input_codes,weights=resampled,minlength=len(input_ids))
        # Keep the distribution of real input quasars fixed in every replicate.
        resampled /= input_mass[input_codes]
        denominator = np.sum(resampled*known)
        if denominator:
            bootstrap.append(np.sum(resampled*passed*factor)/denominator)
    lower, upper = np.quantile(bootstrap, [.025, .975]) if bootstrap else [np.nan, np.nan]
    stochastic_mass = np.sum(weight*known*stochastic)
    raw_probability = np.sum(weight*passed)/stochastic_mass if stochastic_mass else np.nan
    effective_factor = probability/raw_probability if raw_probability > 0 else np.max(factor)
    # A finite number of independent fields cannot establish a zero upper rate.
    # The guard is conservative; mock combinations never increase its sample size.
    if effective > 0 and np.isfinite(raw_probability):
        guard = beta.ppf([.025, .975], raw_probability*effective+.5,
                         (1-raw_probability)*effective+.5)*effective_factor
        lower = min(lower, guard[0]); upper = max(upper, guard[1])
    relative = max(probability-lower, upper-probability)/probability if probability > 0 else np.inf
    unmatched = sample.selection_mode.eq('unmatched').to_numpy()
    missing_pm_fraction=np.nan
    if 'secondary_pm_missing' in sample and 'route' in sample and not sample.route.eq('imaging_one_gaia').all():
        pass_weight=np.sum(weight*passed*factor)
        missing_pm_fraction=(np.sum(weight*passed*factor*sample.secondary_pm_missing.eq(True).to_numpy())/pass_weight
                             if pass_weight else np.nan)
    return dict(f_pass_known=probability, f_pass_population_lower=lower_bound,
        f_pass_population_upper=upper_bound, f_pass_lower_95=lower, f_pass_upper_95=upper,
        relative_halfwidth_95=relative, N_mocks=len(sample),
        N_background_sources=sample.loc[stochastic, 'background_source_id'].nunique(),
        N_background_fields=sample.loc[stochastic, 'background_field_id'].nunique(),
        N_primary_fields=sample.loc[stochastic,'primary_field_id'].nunique() if primary_codes is not None else 0,
        known_weight_fraction=known_mass/mass,
        matching_coverage=1-np.sum(weight*unmatched)/mass, effective_fields=effective,
        passing_weight_with_missing_secondary_pm_fraction=missing_pm_fraction,
        N_passing_mocks=int(passed.sum()), unknown_weight_fraction=np.sum(weight*unknown)/mass)


def summarize(args):
    root = args.output
    pairs = read_csv(root/'input_pairs.csv')
    paths = [root/args.component_scores, root/args.blend_scores]
    scores = pd.concat([pd.read_parquet(path) for path in paths if path.exists()], ignore_index=True)
    factors = read_csv(root/'pair_cut_factors.csv')
    factors['pair_factor'] = factors.dcr_factor*factors.morphology_factor
    scores = scores.merge(factors[['route','latitude_bin','pair_factor']], on=['route','latitude_bin'], validate='many_to_one')
    rng = np.random.default_rng(args.seed)
    rows = []
    for route in ROUTES:
        for annulus in ['3-10', '10-30']:
            for latitude in ['all', '25-40', '40-55', '55-90']:
                real = pairs.loc[pairs.route.eq(route) & pairs.latitude_bin.notna()]
                group = scores.loc[scores.route.eq(route) & scores.annulus.eq(annulus)]
                if latitude != 'all':
                    real = real.loc[real.latitude_bin.eq(latitude)]
                    group = group.loc[group.latitude_bin.eq(latitude)]
                row = estimate_repaired_rate(group, len(real), args.bootstrap, rng)
                row.update(route=route, annulus=annulus, latitude_bin=latitude,
                           N_input=len(real))
                row['precision_met'] = bool(row['relative_halfwidth_95'] <= args.relative_precision)
                rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(root/'repaired_summary_by_latitude.csv', index=False)
    pooled = summary.loc[summary.latitude_bin.eq('all')]
    pooled.to_csv(root/'repaired_summary_by_route.csv', index=False)
    (root/'repaired_precision_status.json').write_text(json.dumps(dict(
        relative_halfwidth_95_target=args.relative_precision,
        pending_score_files=[path.name for path in paths if not path.exists()],
        uncertainty='shared field-cluster multipliers across donor/background roles; finite mock-draw error and finite-field guard',
        missing_fluxes='explicit lower/upper outcome bounds, not assumed negative',
        interpretation='P(pass | unrelated matched control, known outcome); not selected-sample contamination',
        all_main_routes_statistical_target_met=bool(pooled.loc[pooled.annulus.eq('3-10'),'precision_met'].all())), indent=2)+'\n')
    plot(summary, root)
    print(pooled[['route','annulus','N_input','N_background_sources','N_background_fields',
                  'f_pass_known','f_pass_lower_95','f_pass_upper_95','relative_halfwidth_95',
                  'matching_coverage','known_weight_fraction','precision_met']].to_string(index=False), flush=True)


def plot(summary, root):
    os.environ.setdefault('MPLCONFIGDIR', str(root/'scorer_runtime/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 2, figsize=(10, 8), layout='constrained')
    for axis, route in zip(axes.flat, ROUTES):
        for annulus, offset, colour in [('3-10',-.06,'#0072B2'),('10-30',.06,'#D55E00')]:
            group = summary.loc[summary.route.eq(route) & summary.annulus.eq(annulus)
                & summary.latitude_bin.ne('all')]
            value = group.f_pass_known.to_numpy(float)
            errors = [np.maximum(0,value-group.f_pass_lower_95), np.maximum(0,group.f_pass_upper_95-value)]
            axis.errorbar(np.arange(len(group))+offset, value, yerr=errors, fmt='o',
                          color=colour, capsize=3, label=f'{annulus} arcsec')
        axis.set_title(route.replace('_',' ').replace('gaia','Gaia'))
        axis.set_xticks([0,1,2], ['25–40','40–55','55–90'])
        axis.set_xlabel(r'$|b|$ [degrees]'); axis.set_ylabel('Conditional control pass fraction')
        axis.set_ylim(bottom=0); axis.grid(alpha=.2); axis.legend(fontsize=9)
    for extension in ['pdf','png']:
        figure.savefig(root/f'repaired_control_pass_rates.{extension}', dpi=180)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--relative-precision', type=float, default=.2)
    parser.add_argument('--seed', type=int, default=20261121)
    parser.add_argument('--component-scores',default='component_control_scores.parquet')
    parser.add_argument('--blend-scores',default='blend_control_scores.parquet')
    summarize(parser.parse_args())


if __name__ == '__main__':
    main()
