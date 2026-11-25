"""Prevent a shadow import, a newer scorer or local scorer edits from being used."""

from pathlib import Path
import json
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from dual_agn import provenance


class ScorerProvenanceTests(unittest.TestCase):
    def test_scoring_records_verified_provenance_in_working_directory(self):
        from dual_agn import settings
        from dual_agn.photometry import configure_dustmaps
        record = dict(package='qso_pcolor', revision=provenance.MODEL_REVISION,
                      clean_checkout=True, import_relative_path='src/qso_pcolor/__init__.py')
        with TemporaryDirectory() as directory, \
             patch.object(provenance, 'verify_imported_scorer', return_value=record), \
             patch('dustmaps.config.config'), \
             patch.multiple(settings, WORK_ROOT=Path(directory), DUSTMAP_ROOT=Path(directory)/'maps'):
            configure_dustmaps()
            self.assertEqual(json.loads((Path(directory)/'scorer_provenance.json').read_text()), record)

    def test_blend_workers_verify_their_actual_import_before_loading_dustmaps(self):
        from scripts.gaia_neighbour.score_unresolved import load_worker_dustmaps
        with patch('dual_agn.provenance.verify_imported_scorer',
                   side_effect=RuntimeError('wrong worker import')) as verify, \
             patch('dustmaps.config.config') as dustmaps:
            with self.assertRaisesRegex(RuntimeError, 'wrong worker import'):
                load_worker_dustmaps()
            verify.assert_called_once_with()
            dustmaps.load.assert_not_called()

    def check(self, root, revision=None, changes='', origin=None, checkout=None):
        root = Path(root).resolve()
        module = SimpleNamespace(__file__=str(origin or root/'src/qso_pcolor/__init__.py'))
        answers = [str(checkout or root)+'\n', (revision or provenance.MODEL_REVISION)+'\n', changes]
        with patch.object(provenance.importlib, 'import_module', return_value=module), \
             patch.object(provenance.subprocess, 'check_output', side_effect=answers) as git:
            result = provenance.verify_imported_scorer(root)
        self.assertEqual(git.call_args_list[-1].args[0][-3:],
                         ['status', '--porcelain', '--untracked-files=normal'])
        return result

    def test_clean_pinned_import_records_the_actual_source(self):
        with TemporaryDirectory() as root:
            result = self.check(root)
        self.assertEqual(result['revision'], provenance.MODEL_REVISION)
        self.assertTrue(result['clean_checkout'])
        self.assertEqual(result['import_relative_path'], 'src/qso_pcolor/__init__.py')

    def test_newer_revision_is_refused(self):
        with TemporaryDirectory() as root, self.assertRaisesRegex(RuntimeError, 'requires qso_pcolor'):
            self.check(root, revision='a'*40)

    def test_tracked_and_untracked_changes_are_refused(self):
        for changes in [' M src/qso_pcolor/unified.py\n', '?? custom_scorer.py\n']:
            with TemporaryDirectory() as root, self.subTest(changes=changes), \
                 self.assertRaisesRegex(RuntimeError, 'local changes'):
                self.check(root, changes=changes)

    def test_another_checkout_cannot_be_hidden_by_a_correct_configured_checkout(self):
        with TemporaryDirectory() as root, self.assertRaisesRegex(RuntimeError, 'Python imported'):
            self.check(root, origin=Path(root)/'other/qso_pcolor/__init__.py')

    def test_configured_directory_must_be_the_actual_git_root(self):
        with TemporaryDirectory() as root, self.assertRaisesRegex(RuntimeError, 'not rooted'):
            self.check(root, checkout=Path(root).parent)

    def test_namespace_or_unidentifiable_import_is_refused(self):
        with patch.object(provenance.importlib, 'import_module', return_value=SimpleNamespace()), \
             self.assertRaisesRegex(RuntimeError, 'no identifiable'):
            provenance.verify_imported_scorer('/tmp/scorer')
