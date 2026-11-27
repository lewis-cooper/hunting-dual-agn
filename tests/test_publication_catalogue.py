"""Check identity, measurement and selection preservation in publication exports."""

from pathlib import Path
import json
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from export_catalogue import (COLUMN_DESCRIPTIONS, build_publication_catalogue,
                              coordinate_names, curate_candidates, export_catalogue)


def sample_candidates():
    rows = []
    for number, tier in enumerate(['main', 'additional', 'main']):
        row = {name: '' for name in COLUMN_DESCRIPTIONS}
        row.update(candidate_id=f'system:{number}', system_name=f'native:{number}',
                   primary_gaia_source_id=str(9007199254740993 + 2 * number),
                   secondary_gaia_source_id='' if tier == 'additional' else str(9007199254741009 + number),
                   primary_ra_deg=str(15 + number), primary_dec_deg='-0.0005',
                   secondary_ra_deg=str(15.0001 + number), secondary_dec_deg='-0.0005',
                   separation_arcsec='0.36000000001', redshift='1.23456789012345',
                   projected_separation_kpc='3.14159265358979',
                   selected_by_resolved=str(number == 0),
                   selected_by_unresolved=str(number == 2),
                   selected_by_imaging=str(number != 2),
                   selected_by_both_searches=str(number == 0),
                   imaging_selection_route='two_gaia' if number == 0 else 'one_gaia' if number == 1 else '',
                   literature_confirmed_dual=str(number == 0), literature_confirmed_lens='False',
                   literature_uncertain_classification='False', literature_possible_association='False',
                   tier=tier)
        rows.append(row)
    return pd.DataFrame(rows)


class PublicationCatalogueTests(unittest.TestCase):
    def test_exclusion_checks_gaia_identity_not_just_run_ordinal(self):
        source = sample_candidates()
        decision = pd.DataFrame([dict(candidate_id='system:2',
            primary_gaia_source_id=source.iloc[2].primary_gaia_source_id,
            secondary_gaia_source_id=source.iloc[2].secondary_gaia_source_id,
            classification='imaging_artefact', reason='Corrupted data')])
        retained, rejected = curate_candidates(source, decision)
        self.assertEqual(set(retained.candidate_id), {'system:0', 'system:1'})
        self.assertEqual(rejected.candidate_id.tolist(), ['system:2'])
        decision.loc[0, 'secondary_gaia_source_id'] = '9007199254741999'
        with self.assertRaisesRegex(ValueError, 'Gaia identity mismatch'):
            curate_candidates(source, decision)

    def test_renaming_preserves_measurements_ids_and_shared_routes(self):
        source = sample_candidates()
        original = source.copy(deep=True)
        retained, _ = curate_candidates(source, pd.DataFrame())
        publication = build_publication_catalogue(retained)
        pd.testing.assert_frame_equal(source, original)
        indexed = publication.set_index('candidate_id')
        for row in source.itertuples(index=False):
            for field in ['primary_gaia_source_id', 'secondary_gaia_source_id', 'redshift',
                          'separation_arcsec', 'projected_separation_kpc', 'primary_ra_deg']:
                self.assertEqual(indexed.loc[row.candidate_id, field], getattr(row, field))
            self.assertEqual(indexed.loc[row.candidate_id, 'original_system_name'], row.system_name)
        self.assertEqual(indexed.loc['system:0', 'selection_route'], 'gaia_resolved;imaging_two_gaia')

    def test_coordinate_names_handle_ra_wrap_and_negative_near_zero_dec(self):
        names = coordinate_names([0., 359.99999999], [-0.0005, 12.])
        self.assertEqual(names[0], 'J000000.00-000001.8')
        self.assertEqual(names[1], 'J000000.00+120000.0')

    def test_stored_tier_mismatch_and_duplicate_display_names_are_errors(self):
        source = sample_candidates()
        source.loc[0, 'tier'] = 'additional'
        with self.assertRaisesRegex(ValueError, 'Stored tier'):
            curate_candidates(source, pd.DataFrame())
        source = sample_candidates()
        source.loc[1, ['primary_ra_deg', 'primary_dec_deg']] = source.loc[0, ['primary_ra_deg', 'primary_dec_deg']]
        with self.assertRaisesRegex(ValueError, 'names collide'):
            build_publication_catalogue(source)

    def test_saved_tier_files_reconcile_and_automatic_input_is_unchanged(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = sample_candidates()
            input_path = root / 'automatic.csv'
            source.to_csv(input_path, index=False)
            before = input_path.read_bytes()
            publication, summary = export_catalogue(input_path, None, root / 'catalogue', root / 'audit')
            main = pd.read_csv(root / 'catalogue/main_tier.csv', dtype=str, keep_default_na=False)
            additional = pd.read_csv(root / 'catalogue/additional_tier.csv', dtype=str, keep_default_na=False)
            self.assertEqual(input_path.read_bytes(), before)
            self.assertFalse(set(main.candidate_id) & set(additional.candidate_id))
            self.assertEqual(set(main.candidate_id) | set(additional.candidate_id), set(publication.candidate_id))
            self.assertEqual(main.iloc[0].primary_gaia_source_id, source.iloc[0].primary_gaia_source_id)
            self.assertEqual(summary['publication_candidates'], len(source))

    def test_optional_euclid_counts_do_not_replace_catalogue_tier_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / 'automatic.csv'
            sample_candidates().to_csv(input_path, index=False)
            config_path = root / 'euclid.json'
            config_path.write_text(json.dumps({name: name for name in
                ['wide_moc', 'deep_moc', 'grid_records', 'classifications']}))
            with patch('recount_euclid.recount', return_value=(pd.DataFrame(), {'NeuclidVIS': 1})):
                _, summary = export_catalogue(input_path, None, root / 'catalogue',
                                              root / 'audit', euclid_config=config_path)
            self.assertEqual(summary['tiers'], {'main': 2, 'additional': 1})
            self.assertEqual(summary['euclid'], {'NeuclidVIS': 1})
            self.assertIn('2 in the main tier and 1 in the additional tier',
                          (root / 'catalogue/README.md').read_text())


if __name__ == '__main__':
    unittest.main()
