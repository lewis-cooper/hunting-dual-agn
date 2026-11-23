#!/usr/bin/env python3
"""Rerun six component and six blended examples with the actual production models."""

import argparse
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

if __package__:
    from .run_controls import ENGINE, load_config
else:
    from run_controls import ENGINE, load_config


def verify(config, output):
    sys.path.insert(0, str(ENGINE))
    from score_controls import component_scores, configure_production, score_blends
    output.mkdir(parents=True, exist_ok=True)
    examples = config['data_directory']/'scoring_examples'
    args = SimpleNamespace(output=output, repo=config['production_repo'],
        models=config['model_directory'], dustmaps=config['dustmaps_directory'],
        blend_photometry=examples/'independent_blend_photometry.csv',
        blend_grz_diagnostic=False, max_blend_controls=0, resume=False,
        skip_known_pm_failures=True, summed_band_quality=True, production_legacy_first=True,
        skip_incomplete_photometry=True, fixed_survey_policy=True, cache_population_priors=True)
    configure_production(args)
    from qso_pcolor.unified import UnifiedPSFModel
    from score_selection import read_either_model_threshold
    models = {name: UnifiedPSFModel.load(args.models/f'20261005_{name}_qlf')
              for name in ['magdep', 'magindep']}
    threshold = read_either_model_threshold(config['data_directory']/'reference/threshold_calibration.json')
    expected_components = pd.read_parquet(examples/'component_examples.parquet')
    scored_components = component_scores(expected_components, models, threshold, role='imaging')
    differences = {}
    for model in models:
        column = f'log_r_per_unit_z_{model}'
        np.testing.assert_allclose(scored_components[column], expected_components[column], rtol=0, atol=1e-10)
        differences[column] = float(np.max(np.abs(scored_components[column]-expected_components[column])))
    trials = pd.read_parquet(examples/'blended_trials.parquet')
    expected_blends = pd.read_csv(examples/'expected_blended_scores.csv').set_index('trial_id')
    scored_blends = score_blends(trials, args, {'magdep': None, 'magindep': None}).set_index('trial_id')
    expected_blends = expected_blends.loc[scored_blends.index]
    for model in models:
        column = f'log_bf_qq_qs_{model}'
        np.testing.assert_allclose(scored_blends[column], expected_blends[column], rtol=0, atol=1e-10)
        if not scored_blends[f'status_{model}'].eq(expected_blends[f'status_{model}']).all():
            raise AssertionError('The production score status changed.')
        differences[column] = float(np.max(np.abs(scored_blends[column]-expected_blends[column])))
    result = dict(component_examples=len(scored_components), blended_examples=len(scored_blends),
                  threshold=threshold, maximum_absolute_score_differences=differences)
    (output/'model_comparison.txt').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    verify(load_config(args.config), args.output.resolve())
