"""Check portable paths and control bookkeeping before a scientific replay."""

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from chance_alignments.reproduce_results import validate_controls
from chance_alignments.run_controls import load_config, merge_components


class WorkflowTests(unittest.TestCase):
    def write_controls(self, root, duplicate=False, weight=1.):
        pd.DataFrame(dict(pair_id=['resolved:a', 'blended:b'], route=['resolved', 'blended'],
                          latitude_bin=['25-40']*2)).to_csv(root/'input_pairs.csv', index=False)
        for route, identity, filename in [
            ('resolved', 'resolved:a', 'component_control_scores.parquet'),
            ('blended', 'blended:b', 'blend_control_scores.parquet')]:
            rows = pd.DataFrame(dict(trial_id=[identity+':'+annulus for annulus in ['3-10', '10-30']],
                input_pair_id=[identity]*2, route=[route]*2,
                annulus=['3-10', '10-30'], mock_weight=[weight]*2))
            if duplicate and route == 'resolved':
                rows = pd.concat([rows, rows.iloc[:1]], ignore_index=True)
            rows.to_parquet(root/filename, index=False)

    def test_duplicate_control_cannot_enter_replay(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_controls(root, duplicate=True)
            with self.assertRaisesRegex(ValueError, 'Duplicate'):
                validate_controls(root)

    def test_missing_input_cohort_cannot_enter_replay(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_controls(root)
            rows = pd.read_parquet(root/'blend_control_scores.parquet').iloc[:1]
            rows.to_parquet(root/'blend_control_scores.parquet', index=False)
            with self.assertRaisesRegex(ValueError, 'missing or unexpected'):
                validate_controls(root)

    def test_nonunit_input_weight_cannot_enter_replay(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_controls(root, weight=.5)
            with self.assertRaises(AssertionError):
                validate_controls(root)

    def test_config_paths_are_relative_to_config_not_current_directory(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path = root/'config.json'
            path.write_text(json.dumps(dict(data_directory='../controls', output_directory='output',
                production_repo='../production', model_directory='../models', dustmaps_directory='../dustmaps')))
            config = load_config(path)
            self.assertEqual(config['data_directory'], (root/'../controls').resolve())
            self.assertEqual(config['output_directory'], (root/'output').resolve())

    def test_route_extension_replaces_only_requested_routes(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            old, expanded = root/'old', root/'expanded'
            old.mkdir(); expanded.mkdir()
            def rows(batch):
                return pd.DataFrame([dict(route=route, trial_id=batch+':'+route,
                    input_pair_id=route+':input', annulus='3-10', mock_weight=1.)
                    for route in ['resolved', 'imaging_two_gaia', 'imaging_one_gaia']])
            rows('old').to_parquet(old/'component_control_scores.parquet', index=False)
            rows('expanded').to_parquet(expanded/'component_control_scores.parquet', index=False)
            merge_components(old, expanded)
            actual = pd.read_parquet(expanded/'component_control_scores.parquet').set_index('route')
            self.assertEqual(actual.at['imaging_one_gaia', 'trial_id'], 'old:imaging_one_gaia')
            self.assertEqual(actual.at['resolved', 'trial_id'], 'expanded:resolved')
            self.assertEqual(actual.at['imaging_two_gaia', 'trial_id'], 'expanded:imaging_two_gaia')
