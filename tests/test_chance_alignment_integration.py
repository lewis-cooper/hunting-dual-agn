"""Check sequential execution and chance-alignment input connections."""

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_analysis import VETTING_STAGES, commands
from run_chance_alignments import resolve_config, source_paths


class ChanceAlignmentIntegrationTests(unittest.TestCase):
    def configuration(self, directory, **changes):
        values = dict(master_catalogue='master.csv', quaia_catalogue='quaia.fits',
                      workspace_directory='controls_work')
        values.update(changes)
        path = directory/'config.json'
        path.write_text(json.dumps(values))
        return path

    def test_sequential_runner_places_controls_after_vetting_verification(self):
        sequence = commands('all', Path('controls.json'), resume=True)
        modules = [command[command.index('-m')+1] for command in sequence if '-m' in command]
        self.assertEqual([name for name in modules if name.startswith('scripts.')],
                         ['scripts.shared.prepare_quasar_catalogue']+['scripts.'+name for name in VETTING_STAGES])
        self.assertIn('tests/verify_example_run.py', sequence[-2])
        self.assertEqual(sequence[-2][-2:], ['--pipeline', 'both'])
        self.assertEqual(modules[-1], 'run_chance_alignments')
        self.assertIn('--resume', sequence[-1])

    def test_vetting_only_needs_no_control_catalogue_config(self):
        sequence = commands('vetting')
        self.assertTrue(sequence)
        self.assertFalse(any('run_chance_alignments' in command for command in sequence))

    def test_all_stages_require_an_explicit_control_config(self):
        with self.assertRaisesRegex(ValueError, 'chance-config'):
            commands('all')

    def test_relative_input_paths_follow_the_config_location(self):
        with TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            config = resolve_config(self.configuration(root))
            self.assertEqual(config['master_catalogue'], root/'master.csv')
            self.assertEqual(config['quaia_catalogue'], root/'quaia.fits')

    def test_controls_cannot_overwrite_the_vetting_products(self):
        with TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            path = self.configuration(root, products_directory='products',
                                      workspace_directory='products/controls')
            with self.assertRaisesRegex(ValueError, 'distinct chance-alignment workspace'):
                resolve_config(path)

    def test_control_preparation_uses_inputs_and_scores_not_selected_catalogues(self):
        with TemporaryDirectory() as directory:
            config = resolve_config(self.configuration(Path(directory)))
            names = {path.name for path in source_paths(config)}
            self.assertIn('one_arcsecond_scores.csv', names)
            self.assertIn('imaging_scores.csv', names)
            self.assertIn('labelled_systems_audit.csv', names)
            self.assertNotIn('final_candidates.csv', names)
