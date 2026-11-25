"""Check that public root scripts work when launched outside the repository."""

from pathlib import Path
import os
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]


class RootEntrypointTests(unittest.TestCase):
    def test_analysis_script_resolves_its_repository_from_another_directory(self):
        environment = {**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'}
        with TemporaryDirectory() as directory:
            for pipeline in ('both', 'gaia_neighbour', 'imaging'):
                with self.subTest(pipeline=pipeline):
                    completed = subprocess.run([sys.executable, '-B', str(ROOT/'run_analysis.py'),
                        '--pipeline', pipeline, '--dry-run'], cwd=directory, env=environment,
                        capture_output=True, text=True, check=True)
                    self.assertIn(str(ROOT/'quasar_catalogue.local.json'), completed.stdout)
                    self.assertIn('scripts.shared.combine_candidates', completed.stdout)

    def test_control_runner_can_be_invoked_from_another_directory(self):
        with TemporaryDirectory() as directory:
            completed = subprocess.run([sys.executable, '-B', str(ROOT/'run_chance_alignments.py'),
                '--help'], cwd=directory, capture_output=True, text=True, check=True)
            self.assertIn('--replay-data', completed.stdout)


if __name__ == '__main__':
    unittest.main()
