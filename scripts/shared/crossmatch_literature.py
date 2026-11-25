"""Crossmatch literature, retaining known lenses and removing non-lens contaminants."""

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from dual_agn.literature import (LITERATURE_FLAG_COLUMNS, LITERATURE_TEXT_COLUMNS,
                        annotate_candidates, match_references, read_references)
from dual_agn.settings import (FINAL_CANDIDATES_PATH, IDENTIFIER_DTYPES, IMAGING_CUTFLOW_PATH,
                      IMAGING_SCORES_PATH, LITERATURE_ROOT, LITERATURE_WORK_ROOT,
                      ONE_ARCSECOND_CUTFLOW_PATH, ONE_ARCSECOND_SCORES_PATH,
                      PRE_LITERATURE_CANDIDATES_PATH)
from dual_agn.workflow import PIPELINE_CHOICES, selected_pipelines


EXCLUSION_STAGE = 'Confirmed non-lens literature contaminant removal'


def add_score_flags(scores, candidate_audit, score_key, candidate_key):
    source = candidate_audit.loc[candidate_audit[candidate_key].notna()].copy()
    if score_key == 'qso_parent_index':
        source[candidate_key] = pd.to_numeric(source[candidate_key]).astype('Int64').astype('string')
        keys = pd.to_numeric(scores[score_key]).astype('Int64').astype('string')
    else:
        keys = scores[score_key].astype('string')
        source[candidate_key] = source[candidate_key].astype('string')
    if not source[candidate_key].is_unique:
        raise ValueError(f'A combined candidate maps to duplicate {candidate_key} values')
    source = source.set_index(candidate_key)
    result = scores.copy()
    for column in LITERATURE_FLAG_COLUMNS:
        result[column] = keys.map(source[column]).fillna(False).astype(bool)
    for column in LITERATURE_TEXT_COLUMNS:
        result[column] = keys.map(source[column]).fillna('')
    return result


def update_score_audits(candidate_audit, pipeline='both'):
    selected = selected_pipelines(pipeline)
    one = imaging = one_cutflow = imaging_cutflow = None
    if 'gaia_neighbour' in selected:
        one = pd.read_csv(ONE_ARCSECOND_SCORES_PATH, dtype=IDENTIFIER_DTYPES, float_precision='round_trip')
        one = add_score_flags(one, candidate_audit, 'qso_group_id', 'one_arcsecond_qso_group_id')
        one['literature_final_pass'] = (one.selected_by_resolved | one.selected_by_unresolved) & ~one.literature_contaminant_removed
        one_cutflow = pd.read_csv(ONE_ARCSECOND_CUTFLOW_PATH)
        one_cutflow = one_cutflow.loc[~one_cutflow['Selection stage'].isin([EXCLUSION_STAGE, 'Confirmed literature contaminant removal'])].copy()
        one_rows = []
        for route in ['Resolved', 'Unresolved']:
            if route == 'Resolved':
                magdep = one.selected_by_resolved_magdep & ~one.literature_contaminant_removed
                magindep = one.selected_by_resolved_magindep & ~one.literature_contaminant_removed
            else:
                magdep = one.selected_by_unresolved & one.blend_log_bf_qq_qs_magdep.gt(0) & ~one.literature_contaminant_removed
                magindep = one.selected_by_unresolved & one.blend_log_bf_qq_qs_magindep.gt(0) & ~one.literature_contaminant_removed
            one_rows.append({'Route': route, 'Selection stage': EXCLUSION_STAGE,
                             'Magnitude dependent': int(magdep.sum()),
                             'Magnitude independent': int(magindep.sum()),
                             'Either model': int((magdep | magindep).sum())})
        one_cutflow = pd.concat([one_cutflow, pd.DataFrame(one_rows)], ignore_index=True)
    if 'imaging' in selected:
        imaging = pd.read_csv(IMAGING_SCORES_PATH, dtype=IDENTIFIER_DTYPES, float_precision='round_trip')
        imaging = add_score_flags(imaging, candidate_audit, 'qso_parent_index', 'imaging_qso_parent_index')
        imaging['literature_final_pass'] = imaging.joint_final_pass & ~imaging.literature_contaminant_removed
        imaging_cutflow = pd.read_csv(IMAGING_CUTFLOW_PATH)
        imaging_cutflow = imaging_cutflow.loc[~imaging_cutflow['Selection stage'].isin([EXCLUSION_STAGE, 'Confirmed literature contaminant removal'])].copy()
        row = {'Selection stage': EXCLUSION_STAGE,
               'One-Gaia': int((imaging.literature_final_pass & imaging.candidate_tier.eq('one_gaia')).sum()),
               'Two-Gaia': int((imaging.literature_final_pass & imaging.candidate_tier.eq('two_gaia')).sum()),
               'Total': int(imaging.literature_final_pass.sum())}
        imaging_cutflow = pd.concat([imaging_cutflow, pd.DataFrame([row])], ignore_index=True)
    return one, imaging, one_cutflow, imaging_cutflow


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidates', type=Path, default=PRE_LITERATURE_CANDIDATES_PATH,
                        help='Unfiltered combined working catalogue from combine_candidate_catalogues')
    parser.add_argument('--binaries', type=Path, default=LITERATURE_ROOT/'known_binaries.csv')
    parser.add_argument('--contaminants', type=Path, default=LITERATURE_ROOT/'known_lenses_contaminants.csv')
    parser.add_argument('--search-radius-arcsec', type=float, default=1.0)
    parser.add_argument('--pipeline', choices=PIPELINE_CHOICES, default='both')
    args = parser.parse_args()
    if not args.candidates.is_file():
        parser.error(f'Missing {args.candidates}; run python -m scripts.shared.combine_candidates first.')
    candidates = pd.read_csv(args.candidates, dtype=IDENTIFIER_DTYPES, float_precision='round_trip')
    references = read_references(args.binaries, args.contaminants)
    associations = match_references(candidates, references, args.search_radius_arcsec)
    candidate_audit = annotate_candidates(candidates, associations)
    removed = candidate_audit.literature_contaminant_removed
    final_candidates = candidate_audit.loc[~removed].copy()
    associations['exclusion_applied'] = associations.candidate_id.isin(candidate_audit.loc[removed, 'candidate_id']) & associations.remove_candidate
    one, imaging, one_cutflow, imaging_cutflow = update_score_audits(candidate_audit, args.pipeline)
    counts = {'before_literature': len(candidate_audit), 'confirmed_contaminants_removed': int(removed.sum()),
              'final_candidates': len(final_candidates),
              'resolved': int(final_candidates.selected_by_resolved.sum()),
              'unresolved': int(final_candidates.selected_by_unresolved.sum()),
              'imaging': int(final_candidates.selected_by_imaging.sum()),
              'shared': int(final_candidates.selected_by_both_searches.sum()),
              'confirmed_dual_matches': int(final_candidates.literature_confirmed_dual.sum()),
              'confirmed_lenses_retained': int(final_candidates.literature_confirmed_lens.sum()),
              'uncertain_classification_flags': int(final_candidates.literature_uncertain_classification.sum())}
    if counts['resolved'] + counts['unresolved'] + counts['imaging'] - counts['shared'] != counts['final_candidates']:
        raise ValueError('Cleaned route counts do not reconcile with the unique catalogue')
    gaia_count = int(one.literature_final_pass.sum()) if one is not None else 0
    imaging_count = int(imaging.literature_final_pass.sum()) if imaging is not None else 0
    if counts['resolved'] + counts['unresolved'] != gaia_count or counts['imaging'] != imaging_count:
        raise ValueError('Score-audit flags differ from the cleaned combined catalogue')
    summary = {
        'pipeline': args.pipeline, 'counts': counts, 'search_radius_arcsec': args.search_radius_arcsec,
        'redshift_used_for_matching': False,
        'binary_reference_records': int(references.reference_catalogue.eq('known_binaries').sum()),
        'contaminant_reference_records': int(references.reference_catalogue.eq('known_lenses_contaminants').sum()),
        'references_without_positions': references.loc[references.ra_deg.isna(), ['reference_id', 'system_name']].to_dict('records'),
        'reference_sha256': {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in [args.binaries,args.contaminants]},
        'input_sha256': hashlib.sha256(args.candidates.read_bytes()).hexdigest(),
        'removal_policy': 'Confirmed lenses retained and flagged; only verified confirmed non-lens contaminants removed; uncertain or conflicting cases retained.',
    }
    LITERATURE_WORK_ROOT.mkdir(parents=True, exist_ok=True)
    associations.to_csv(LITERATURE_WORK_ROOT/'associations.csv', index=False)
    candidate_audit.to_csv(LITERATURE_WORK_ROOT/'candidate_audit.csv', index=False)
    (LITERATURE_WORK_ROOT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    if one is not None:
        one.to_csv(ONE_ARCSECOND_SCORES_PATH, index=False)
        one_cutflow.to_csv(ONE_ARCSECOND_CUTFLOW_PATH, index=False)
    if imaging is not None:
        imaging.to_csv(IMAGING_SCORES_PATH, index=False)
        imaging_cutflow.to_csv(IMAGING_CUTFLOW_PATH, index=False)
    final_candidates.to_csv(FINAL_CANDIDATES_PATH, index=False)
    print('Input:', args.candidates.resolve())
    print('References:', args.binaries.resolve(), args.contaminants.resolve())
    print(json.dumps(summary, indent=2))
    print('Final catalogue:', FINAL_CANDIDATES_PATH)
    print('Match/exclusion audit:', LITERATURE_WORK_ROOT)


if __name__ == '__main__':
    main()
