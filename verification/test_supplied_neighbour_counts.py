"""Supplied association/count exports must enter discovery without raw-Gaia adaptation."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from run_analysis import commands
from scripts.shared.check_local_inputs import check_inputs


class SuppliedNeighbourCountsTests(unittest.TestCase):
    def test_all_scoring_modes_use_supplied_counts_without_loading_a_raw_master(self):
        for pipeline in ('both', 'gaia_neighbour', 'imaging'):
            sequence = commands('vetting', pipeline=pipeline,
                one_arcsecond_search_config=Path('counts_input.json'),
                imaging_search_config=Path('imaging_input.json'))
            self.assertFalse(any('scripts.shared.prepare_quasar_catalogue' in row for row in sequence))
            discovery = next(row for row in sequence if 'scripts.gaia_neighbour.search_candidates' in row)
            self.assertEqual(Path(discovery[-1]).name, 'counts_input.json')
            if pipeline != 'gaia_neighbour':
                imaging = next(row for row in sequence if 'scripts.imaging.search_candidates' in row)
                self.assertEqual(Path(imaging[-1]).name, 'imaging_input.json')

    def test_discovery_preflight_accepts_exports_without_quasar_gaia_placeholder(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'matches.csv').write_text('qso_uid\n')
            (root/'counts.csv').write_text('gaia_source_id,gaia_edr3_cnt1\n')
            config = root/'search.json'
            config.write_text(json.dumps(dict(input_mode='saved_qualifying_matches',
                quasar_gaia_matches=['matches.csv'], gaia_neighbour_counts='counts.csv')))
            checked = check_inputs(root/'absent_quasars.json', None, None,
                'gaia_neighbour', 'discovery', one_arcsecond_search_config=config)
            self.assertEqual(checked, {(root/'matches.csv').resolve(), (root/'counts.csv').resolve()})

    def test_imaging_discovery_remains_independent_of_gaia_calibration_inputs(self):
        sequence = commands('discovery', pipeline='imaging',
            one_arcsecond_search_config=Path('absent_counts.json'),
            imaging_search_config=Path('imaging.json'))
        self.assertEqual(len(sequence), 1)
        self.assertIn('scripts.imaging.search_candidates', sequence[0])


if __name__ == '__main__':
    unittest.main()
