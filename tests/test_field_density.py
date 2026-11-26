"""Same-field density counts, exact identities and portable replay checks."""

import hashlib
import io
import json
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from chance_alignments.engine.analyze_field_density import (
    compare_counts, export_replay_inputs, prepare_fields, write_report,
)
from chance_alignments.reproduce_results import reproduce
from chance_alignments.run_controls import run


class MatchedDensityTests(unittest.TestCase):
    def write_cones(self, root):
        identities = ['9007199254740993', '9007199254740995', '9007199254740997']
        fields = pd.DataFrame(dict(target_id=['a', 'b', 'c'], source_id=identities,
            ra=[1., 2., 3.], dec=[30.]*3, sampling_weight=[2., 8., 5.], sampling_cell=['cell']*3))
        fields.to_csv(root/'master_targets.csv', index=False)
        coverage = fields[['target_id', 'ra', 'dec']].copy()
        coverage['gaia_query_complete'] = True
        coverage.to_csv(root/'cone_coverage.csv', index=False)
        neighbours = pd.DataFrame(dict(target_id=['a', 'b', 'c', 'a', 'a', 'b', 'b'],
            source_id=identities+['9100000000000001', '9100000000000003', '9100000000000005', '9100000000000007'],
            separation_arcsec=[.02, .02, .02, .5, 20., 35., 40.]))
        neighbours.to_parquet(root/'gaia_neighbours.parquet', index=False)
        self.write_metadata(root)

    def write_metadata(self, root, radius=60):
        (root/'query_configuration.json').write_text(json.dumps(dict(radius_arcsec=radius,
            parent_targets_sha256=hashlib.sha256((root/'master_targets.csv').read_bytes()).hexdigest())))

    def test_central_identity_exclusion_and_zero_neighbour_fields(self):
        with TemporaryDirectory() as directory:
            root = Path(directory); self.write_cones(root)
            fields, neighbours = prepare_fields(root)
            self.assertEqual(fields.source_id.iloc[0], '9007199254740993')
            self.assertEqual(len(neighbours), 4)
            summary, counts = compare_counts(fields, neighbours, [1., 1.5], 100, 7)
            row = summary.iloc[0]
            self.assertEqual(row.n_fields, 3)
            self.assertEqual(row.observed_sample_companions, 1)
            self.assertEqual(row.outer_sample_neighbours, 3)
            self.assertAlmostEqual(row.expected_sample_companions_C1, 3/3500)
            self.assertAlmostEqual(row.weighted_observed_companions, 2)
            self.assertAlmostEqual(row.weighted_expected_companions_C1, 18/3500)
            self.assertEqual(counts.loc[counts.target_id.eq('c'), 'neighbours_10_60_arcsec'].iloc[0], 0)
            # Sample and weighted ratios differ for unequal inclusion weights.
            self.assertNotAlmostEqual(row.observed_to_C1_ratio, 3500/3)

    def test_area_scaling_and_split_annuli(self):
        with TemporaryDirectory() as directory:
            root = Path(directory); self.write_cones(root)
            fields, neighbours = prepare_fields(root)
            summary, _ = compare_counts(fields, neighbours, [1., 1.5], 100, 7)
            full = summary.loc[summary.outer_upper_arcsec.eq(60) & summary.outer_lower_arcsec.eq(10)]
            self.assertAlmostEqual(full.iloc[1].expected_sample_companions_C1 /
                                   full.iloc[0].expected_sample_companions_C1, 1.5**2)
            self.assertEqual(set(zip(summary.outer_lower_arcsec, summary.outer_upper_arcsec)),
                             {(10., 60.), (10., 30.), (30., 60.)})

    def test_exported_raw_inputs_reproduce_counts_and_bootstrap(self):
        with TemporaryDirectory() as directory:
            source, replay = Path(directory)/'source', Path(directory)/'replay'
            source.mkdir(); self.write_cones(source)
            export_replay_inputs(source, replay)
            original = compare_counts(*prepare_fields(source), [1., 1.5], 100, 7)
            rebuilt = compare_counts(*prepare_fields(replay), [1., 1.5], 100, 7)
            pd.testing.assert_frame_equal(original[0], rebuilt[0])
            pd.testing.assert_frame_equal(original[1], rebuilt[1])

    def test_failed_or_shifted_cone_is_rejected(self):
        with TemporaryDirectory() as directory:
            root = Path(directory); self.write_cones(root)
            coverage = pd.read_csv(root/'cone_coverage.csv')
            for change in ('failure', 'coordinates'):
                invalid = coverage.copy()
                if change == 'failure':
                    invalid.loc[0, 'gaia_query_complete'] = False
                else:
                    invalid.loc[0, 'ra'] += .001
                invalid.to_csv(root/'cone_coverage.csv', index=False)
                with self.assertRaisesRegex(ValueError, 'successful Gaia cone'):
                    prepare_fields(root)

    def test_metadata_rejects_short_or_different_field_cones(self):
        with TemporaryDirectory() as directory:
            root = Path(directory); self.write_cones(root)
            self.write_metadata(root, radius=30)
            with self.assertRaisesRegex(ValueError, '60 arcsec'):
                prepare_fields(root)
            self.write_metadata(root)
            with (root/'master_targets.csv').open('a') as stream:
                stream.write('\n')
            with self.assertRaisesRegex(ValueError, 'same sampled parent'):
                prepare_fields(root)

    def test_missing_duplicate_or_unknown_gaia_identity_is_rejected(self):
        with TemporaryDirectory() as directory:
            root = Path(directory); self.write_cones(root)
            original = pd.read_parquet(root/'gaia_neighbours.parquet')
            invalid = [original.iloc[1:], pd.concat([original, original.iloc[:1]], ignore_index=True),
                       original.assign(target_id=original.target_id.replace('a', 'unknown'))]
            for neighbours in invalid:
                neighbours.to_parquet(root/'gaia_neighbours.parquet', index=False)
                with self.assertRaises(ValueError):
                    prepare_fields(root)

    def test_floating_gaia_ids_and_invalid_weights_are_rejected(self):
        with TemporaryDirectory() as directory:
            root = Path(directory); self.write_cones(root)
            neighbours = pd.read_parquet(root/'gaia_neighbours.parquet')
            neighbours.source_id = neighbours.source_id.astype(float)
            neighbours.to_parquet(root/'gaia_neighbours.parquet', index=False)
            with self.assertRaisesRegex(ValueError, 'never floating point'):
                prepare_fields(root)
            self.write_cones(root)
            fields = pd.read_csv(root/'master_targets.csv', dtype={'source_id': 'string'})
            fields.loc[0, 'sampling_weight'] = 0
            fields.to_csv(root/'master_targets.csv', index=False); self.write_metadata(root)
            with self.assertRaisesRegex(ValueError, 'positive inclusion weights'):
                prepare_fields(root)

    def test_empty_outer_annulus_does_not_claim_a_measured_ratio(self):
        with TemporaryDirectory() as directory:
            root = Path(directory); self.write_cones(root)
            fields, neighbours = prepare_fields(root)
            summary, _ = compare_counts(fields, neighbours.loc[neighbours.separation_arcsec.lt(10)], [1.], 100, 7)
            self.assertTrue(summary.observed_to_C1_ratio.isna().all())
            self.assertTrue(summary.ratio_upper_95.isna().all())

    def test_density_report_separates_raw_and_weighted_quantities(self):
        with TemporaryDirectory() as directory:
            root = Path(directory); self.write_cones(root)
            summary, _ = compare_counts(*prepare_fields(root), [1., 1.5], 100, 7)
            write_report(summary, root)
            report = (root/'FIELD_DENSITY.md').read_text()
            self.assertIn('same 3 Gaia-centred', report)
            self.assertIn('Predicted in sampled fields', report)
            self.assertIn('Weighted predicted', report)
            self.assertIn('not a measured unrelated fraction', report)

    def test_report_stage_runs_comparison_and_exports_replay_inputs(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            config = {key: root/key for key in ['data_directory', 'output_directory',
                'production_repo', 'model_directory', 'dustmaps_directory']}
            config.update(workers=1, bootstrap=100)
            output = io.StringIO()
            with redirect_stdout(output):
                run(config, 'report', dry_run=True)
            plan = output.getvalue()
            self.assertIn('analyze_field_density.py', plan)
            self.assertIn('--export-inputs', plan)
            self.assertLess(plan.index('analyze_field_density.py'), plan.index('write_precision_report.py'))

    def test_old_replay_archive_requires_explicit_density_inputs(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(FileNotFoundError, '--density-source'):
                reproduce(root/'old_data', root/'output')

    def test_replay_runs_matched_comparison_with_supplied_cones(self):
        with TemporaryDirectory() as directory:
            root = Path(directory); data = root/'data'; source = root/'cones'
            (data/'replay').mkdir(parents=True); source.mkdir(); self.write_cones(source)
            for name in ['input_pairs.csv', 'pair_cut_factors.csv',
                         'component_control_scores.parquet', 'blend_control_scores.parquet']:
                (data/'replay'/name).touch()
            with patch('chance_alignments.reproduce_results.validate_controls'), \
                 patch('chance_alignments.reproduce_results.subprocess.run') as execute:
                reproduce(data, root/'output', density_source=source)
            commands = [call.args[0] for call in execute.call_args_list]
            self.assertTrue(commands[0][2].endswith('analyze_field_density.py'))
            self.assertIn(str(source.resolve()), commands[0])
            self.assertIn('--export-inputs', commands[0])


if __name__ == '__main__':
    unittest.main()
