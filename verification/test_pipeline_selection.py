"""Exercise pipeline routing, isolated outputs, shared pairs and literature removal."""

from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from dual_agn.workflow import product_names
from run_analysis import commands
from scripts.shared import combine_candidates as catalogue
from scripts.shared import crossmatch_literature as literature
from verification import test_literature


def gaia_scores():
    rows = []
    for index, ra in enumerate((10., 20.)):
        rows.append(dict(qso_group_id=f'g{index}', primary_qso_name=f'gaia{index}',
            primary_qso_z=1., primary_qso_z_original=1., primary_qso_zerr=0.,
            primary_qso_redshift_kind='spectroscopic', primary_qso_catalog='SDSS_DR16Q',
            redshift_catalogue_release='SDSS DR16Q', selected_by_resolved=True,
            selected_by_unresolved=False, selected_by_resolved_magdep=True,
            selected_by_resolved_magindep=True, secondary_pm_significance=0.,
            secondary_pm_missing=False, secondary_pm_pass=True,
            resolved_log_r_per_unit_z_magdep=-1., resolved_log_r_per_unit_z_magindep=-1.,
            resolved_log_r_per_unit_z_either=-1., blend_log_bf_qq_qs_magdep=np.nan,
            blend_log_bf_qq_qs_magindep=np.nan, primary_component_ra=ra,
            primary_component_dec=0., secondary_component_ra=ra+.0002,
            secondary_component_dec=0., primary_component_gaia_source_id=str(100+2*index),
            secondary_component_gaia_source_id=str(101+2*index)))
    return pd.DataFrame(rows)


def imaging_scores():
    rows = []
    for index, ra in enumerate((10., 30.)):
        rows.append(dict(qso_parent_index=index, joint_final_pass=True,
            candidate_tier='two_gaia', sdss_name=f'image{index}', adopted_redshift=1.,
            adopted_redshift_catalogue='SDSS_DR16Q', adopted_redshift_error=0.,
            secondary_pm_significance=0., secondary_pm_missing=False, secondary_pm_pass=True,
            primary_log_r_per_unit_z_magdep=-1., secondary_log_r_per_unit_z_magdep=-1.,
            primary_log_r_per_unit_z_magindep=-1., secondary_log_r_per_unit_z_magindep=-1.,
            minimum_component_log_r_either=-1., passes_pair_score_magdep=True,
            passes_pair_score_magindep=True, pair_dcr_delta_deg=10., dcr_cut_applied=False,
            r_perp_sigma_ratio=1.2, morphology_measurement_source='fixture',
            primary_ra=ra, primary_dec=0., secondary_ra=ra+.0002, secondary_dec=0.,
            primary_gaia_source_id=str(100+4*index), secondary_gaia_source_id=str(101+4*index)))
    return pd.DataFrame(rows)


class PipelineSelectionTests(unittest.TestCase):
    def modules(self, pipeline):
        return [command[command.index('-m')+1] for command in commands('vetting', pipeline=pipeline) if '-m' in command]

    def test_gaia_only_omits_imaging_search_and_selection(self):
        modules = self.modules('gaia_neighbour')
        self.assertIn('scripts.gaia_neighbour.score_unresolved', modules)
        self.assertIn('scripts.gaia_neighbour.select_candidates', modules)
        self.assertFalse(any(name.startswith('scripts.imaging.') for name in modules))

    def test_imaging_only_calibrates_without_selecting_gaia_candidates(self):
        modules = self.modules('imaging')
        self.assertIn('scripts.gaia_neighbour.score_resolved', modules)
        self.assertIn('scripts.gaia_neighbour.calibrate_threshold', modules)
        self.assertIn('scripts.imaging.select_candidates', modules)
        self.assertNotIn('scripts.gaia_neighbour.score_unresolved', modules)
        self.assertNotIn('scripts.gaia_neighbour.select_candidates', modules)
        self.assertLess(modules.index('scripts.gaia_neighbour.calibrate_threshold'),
                        modules.index('scripts.imaging.score_components'))
        preparation = next(command for command in commands('vetting', pipeline='imaging')
                           if 'scripts.gaia_neighbour.build_measurements' in command)
        self.assertIn('--calibration-only', preparation)

    def test_pipeline_is_explicit_in_combination_literature_and_verification(self):
        for pipeline in ('both', 'gaia_neighbour', 'imaging'):
            for command in commands('vetting', pipeline=pipeline):
                if any(name in command for name in ('scripts.shared.combine_candidates',
                    'scripts.shared.crossmatch_literature', 'verification/verify_example_run.py')):
                    self.assertEqual(command[-2:], ['--pipeline', pipeline])
            self.assertEqual(len(product_names(pipeline)), 6 if pipeline == 'both' else 4)

    def test_control_comparison_requires_both_candidate_pipelines(self):
        with self.assertRaisesRegex(ValueError, 'compares both pipelines'):
            commands('all', chance_config=Path('controls.json'), pipeline='imaging')

    def run_output_stages(self, root, pipeline):
        products, work = root/'products', root/'work'
        products.mkdir(); work.mkdir()
        paths = dict(ONE_ARCSECOND_SCORES_PATH=products/'one_arcsecond_scores.csv',
            IMAGING_SCORES_PATH=products/'imaging_scores.csv',
            ONE_ARCSECOND_CUTFLOW_PATH=products/'one_arcsecond_cutflow.csv',
            IMAGING_CUTFLOW_PATH=products/'imaging_cutflow.csv',
            FINAL_CANDIDATES_PATH=products/'final_candidates.csv',
            PRE_LITERATURE_CANDIDATES_PATH=work/'before_literature.csv')
        for selected in ('gaia_neighbour', 'imaging'):
            active = pipeline in ('both', selected)
            score_path = paths['ONE_ARCSECOND_SCORES_PATH' if selected == 'gaia_neighbour' else 'IMAGING_SCORES_PATH']
            if active:
                (gaia_scores() if selected == 'gaia_neighbour' else imaging_scores()).to_csv(score_path,index=False)
                cutflow_path = paths['ONE_ARCSECOND_CUTFLOW_PATH' if selected == 'gaia_neighbour' else 'IMAGING_CUTFLOW_PATH']
                columns = ['Route', 'Selection stage', 'Magnitude dependent', 'Magnitude independent', 'Either model'] if selected == 'gaia_neighbour' else ['Selection stage', 'One-Gaia', 'Two-Gaia', 'Total']
                pd.DataFrame(columns=columns).to_csv(cutflow_path,index=False)
            else:
                # An old unrelated output must never become an implicit input.
                score_path.write_text('invalid inactive-pipeline file\n')
        reference = test_literature.LiteratureTests().reference()
        reference['secondary_redshift'] = np.nan
        reference.to_csv(root/'contaminants.csv',index=False)
        reference.iloc[:0].to_csv(root/'binaries.csv',index=False)
        with ExitStack() as stack:
            stack.enter_context(patch.multiple(catalogue, OUTPUT_ROOT=products, WORK_ROOT=work,
                **{name:value for name,value in paths.items() if hasattr(catalogue,name)}))
            stack.enter_context(patch.multiple(literature, LITERATURE_WORK_ROOT=work/'literature', **paths))
            combined = catalogue.combine_candidates(pipeline)
            stack.enter_context(patch('sys.argv', ['crossmatch', '--pipeline', pipeline,
                '--binaries', str(root/'binaries.csv'), '--contaminants', str(root/'contaminants.csv')]))
            literature.main()
        final = pd.read_csv(paths['FINAL_CANDIDATES_PATH'])
        self.assertEqual(len(combined), 3 if pipeline == 'both' else 2)
        self.assertEqual(len(final), 2 if pipeline == 'both' else 1)
        self.assertFalse(final.literature_contaminant_removed.any())
        if pipeline == 'gaia_neighbour':
            self.assertFalse(final.selected_by_imaging.any())
            self.assertTrue(final.selected_by_resolved.all())
        elif pipeline == 'imaging':
            self.assertFalse((final.selected_by_resolved | final.selected_by_unresolved).any())
            self.assertTrue(final.selected_by_imaging.all())
        else:
            self.assertEqual(int(combined.selected_by_both_searches.sum()),1)

    def test_each_mode_assembles_only_requested_candidates_and_removes_confirmed_lenses(self):
        for pipeline in ('both', 'gaia_neighbour', 'imaging'):
            with self.subTest(pipeline=pipeline), TemporaryDirectory() as directory:
                self.run_output_stages(Path(directory), pipeline)
