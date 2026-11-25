"""Regress native blend metadata compatibility and exact checkpoint recovery."""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock

import numpy as np
import pandas as pd

from dual_agn.blend_outputs import (BLEND_METADATA_COLUMNS, OPTIONAL_DESI_METADATA,
    build_blend_metadata, checkpoint_blend_scores, scoring_fingerprint)


class BlendOutputTests(unittest.TestCase):
    def test_native_metadata_accepts_missing_historical_desi_labels(self):
        candidates = pd.DataFrame({field: [1] for field in BLEND_METADATA_COLUMNS
                                   if field not in OPTIONAL_DESI_METADATA})
        candidates['desi_targetids_json'] = '["39627392414648236"]'
        metadata = build_blend_metadata(candidates)
        self.assertTrue(metadata[list(OPTIONAL_DESI_METADATA)].isna().all().all())
        self.assertEqual(metadata.desi_targetids_json.iloc[0], '["39627392414648236"]')
        self.assertEqual(metadata.qso_group_id.tolist(), candidates.qso_group_id.tolist())

    def test_missing_scientific_fields_raise_before_scoring(self):
        with self.assertRaisesRegex(ValueError, 'desi_release'):
            build_blend_metadata(pd.DataFrame({'qso_group_id': ['a']}))

    def test_historical_metadata_is_unchanged(self):
        candidates = pd.DataFrame({field: [1] for field in BLEND_METADATA_COLUMNS})
        pd.testing.assert_frame_equal(candidates, build_blend_metadata(candidates))

    def test_completed_batch_is_reused_without_computation_or_precision_loss(self):
        rows = [dict(qso_group_id='a', model='magdep',
                     log_bf_qq_qs=-1.2345678901234567, finite_nodes=np.int64(9),
                     notes=('legacy_r_marginalised:ps1:r',), redshift_error=np.nan)]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'batch.json'
            checkpoint_blend_scores(path, 'magdep', ['a'], lambda: rows)
            compute = Mock(side_effect=AssertionError('Must not recompute'))
            recovered = checkpoint_blend_scores(path, 'magdep', ['a'], compute)
            compute.assert_not_called()
            self.assertEqual(recovered[0]['log_bf_qq_qs'], rows[0]['log_bf_qq_qs'])
            self.assertEqual(recovered[0]['notes'], rows[0]['notes'])
            self.assertTrue(np.isnan(recovered[0]['redshift_error']))

    def test_changed_identity_cannot_reuse_checkpoint(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'batch.json'
            checkpoint_blend_scores(path, 'magdep', ['a'],
                lambda: [dict(qso_group_id='a', model='magdep')])
            with self.assertRaisesRegex(ValueError, 'identities'):
                checkpoint_blend_scores(path, 'magdep', ['b'], lambda: [])

    def test_failed_batch_does_not_write_a_completed_checkpoint(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'batch.json'
            compute = Mock(side_effect=RuntimeError('Interrupted'))
            with self.assertRaises(RuntimeError):
                checkpoint_blend_scores(path, 'magdep', ['a'], compute)
            self.assertFalse(path.exists())

    def test_changed_inputs_or_scorer_invalidate_checkpoints(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'measurements.csv'
            path.write_text('flux\n1\n')
            before = scoring_fingerprint([path], 'revision_a')
            self.assertNotEqual(before, scoring_fingerprint([path], 'revision_b'))
            path.write_text('flux\n2\n')
            self.assertNotEqual(before, scoring_fingerprint([path], 'revision_a'))


if __name__ == '__main__':
    unittest.main()
