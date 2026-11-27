"""Export a curated publication catalogue, its two tiers and paper previews.

This is a post-selection export. It does not score sources, alter the automatic
catalogue, or change the selection rules. Exclusions are explicit local inputs.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from dual_agn.candidate_tiers import assign_candidate_tiers, boolean_flag


ROOT = Path(__file__).resolve().parent
COLUMN_DESCRIPTIONS = {
    'candidate_id': ('', 'Stable identifier linking the system to the automatic score audits.'),
    'system_name': ('', 'J-coordinate name from the adopted primary ICRS position; RA seconds rounded to two decimals and Dec seconds to one.'),
    'original_system_name': ('', 'Unmodified name or native catalogue identifier in the automatic catalogue.'),
    'tier': ('', 'main: two distinct Gaia detections; additional: one-Gaia imaging selection.'),
    'selection_route': ('', 'All selecting routes, separated by semicolons: gaia_resolved, gaia_unresolved, imaging_one_gaia or imaging_two_gaia.'),
    'primary_ra_deg': ('deg', 'Adopted primary right ascension, ICRS.'),
    'primary_dec_deg': ('deg', 'Adopted primary declination, ICRS.'),
    'secondary_ra_deg': ('deg', 'Adopted secondary right ascension, ICRS.'),
    'secondary_dec_deg': ('deg', 'Adopted secondary declination, ICRS.'),
    'separation_arcsec': ('arcsec', 'Angular component separation from the automatic catalogue.'),
    'redshift': ('', 'Adopted system redshift; a shared component redshift is an assumption, not confirmation.'),
    'redshift_catalogue': ('', 'Catalogue supplying the adopted redshift.'),
    'adopted_redshift_error_used': ('', 'Additional redshift uncertainty used in scoring; zero for spectroscopic redshifts.'),
    'redshift_treatment': ('', 'Fixed spectroscopic redshift or marginalisation over a photometric redshift uncertainty.'),
    'projected_separation_kpc': ('kpc', 'Projected separation assuming both components share the adopted redshift; Planck18 cosmology.'),
    'primary_gaia_source_id': ('', 'Exact Gaia DR3 identifier of the primary, stored as text.'),
    'secondary_gaia_source_id': ('', 'Exact Gaia DR3 identifier of the secondary; blank when absent.'),
    'selected_by_resolved': ('', 'True when selected by the resolved Gaia-neighbour route.'),
    'selected_by_unresolved': ('', 'True when selected by the unresolved Gaia-neighbour route.'),
    'selected_by_imaging': ('', 'True when selected by the imaging search.'),
    'selected_by_both_searches': ('', 'True when selected by both Gaia-neighbour and imaging searches.'),
    'imaging_selection_route': ('', 'one_gaia or two_gaia for imaging selections; blank otherwise.'),
    'secondary_pm_significance': ('', 'Secondary Gaia proper-motion significance, when measured.'),
    'secondary_pm_missing': ('', 'True when a two-Gaia system has no usable secondary proper-motion significance.'),
    'n_psf_components_within_aperture': ('', 'Number of Tractor PSF sources inside the imaging search aperture, when available.'),
    'literature_class': ('', 'dual_agn, lens, uncertain or unclassified; literature identifications are not new confirmations.'),
    'literature_reference_ids': ('', 'Identifiers of the matching records in the literature comparison catalogues.'),
    'literature_author_year': ('', 'Author/year labels for literature classifications.'),
    'literature_publication_urls': ('', 'Source-paper links for literature classifications.'),
}
ROUTE_LABELS = {
    'gaia_resolved': 'R', 'gaia_unresolved': 'U',
    'imaging_one_gaia': 'I1', 'imaging_two_gaia': 'I2',
}


def read_csv(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def coordinate_names(ra_deg, dec_deg):
    """Use explicit, uniform display precision while retaining full coordinates."""
    coordinates = SkyCoord(ra=np.asarray(ra_deg, float) * u.deg,
                           dec=np.asarray(dec_deg, float) * u.deg)
    ra = coordinates.ra.to_string(unit=u.hourangle, sep='', precision=2, pad=True)
    dec = coordinates.dec.to_string(unit=u.deg, sep='', precision=1,
                                    pad=True, alwayssign=True)
    return ['J' + ('000000.00' if r == '240000.00' else r) + d
            for r, d in zip(np.atleast_1d(ra), np.atleast_1d(dec))]


def curate_candidates(candidates, exclusions):
    """Exclude only explicitly identified rows; require their Gaia identities too."""
    if candidates.candidate_id.eq('').any() or not candidates.candidate_id.is_unique:
        raise ValueError('Candidate IDs must be nonempty and unique')
    if not exclusions.empty:
        required = ['candidate_id', 'primary_gaia_source_id', 'secondary_gaia_source_id',
                    'classification', 'reason']
        if not set(required).issubset(exclusions):
            raise ValueError(f'Exclusions require {required}')
        if not exclusions.candidate_id.is_unique:
            raise ValueError('Duplicate exclusion identifiers')
        unknown = set(exclusions.candidate_id) - set(candidates.candidate_id)
        if unknown:
            raise ValueError(f'Exclusions absent from this catalogue: {sorted(unknown)}')
        indexed = candidates.set_index('candidate_id')
        for row in exclusions.itertuples(index=False):
            actual = indexed.loc[row.candidate_id]
            observed_ids = {actual.primary_gaia_source_id, actual.secondary_gaia_source_id}
            expected_ids = {row.primary_gaia_source_id, row.secondary_gaia_source_id}
            if '' in expected_ids or observed_ids != expected_ids:
                raise ValueError(f'Exclusion Gaia identity mismatch: {row.candidate_id}')
    excluded_ids = set(exclusions.candidate_id) if not exclusions.empty else set()
    rejected = candidates.candidate_id.isin(excluded_ids)
    retained = candidates.loc[~rejected].copy()
    assigned = assign_candidate_tiers(retained)
    if 'tier' in retained and not assigned.tier.eq(retained.tier).all():
        raise ValueError('Stored tier disagrees with Gaia component identities')
    rejected_rows = candidates.loc[rejected].copy()
    if not exclusions.empty:
        rejected_rows = rejected_rows.merge(
            exclusions[['candidate_id', 'classification', 'reason']],
            on='candidate_id', validate='one_to_one')
    return assigned, rejected_rows


def build_publication_catalogue(candidates):
    result = candidates.copy()
    numeric_columns = ['primary_ra_deg', 'primary_dec_deg', 'secondary_ra_deg',
                       'secondary_dec_deg', 'redshift', 'separation_arcsec',
                       'projected_separation_kpc']
    numeric = {name: pd.to_numeric(result[name], errors='raise') for name in numeric_columns}
    for name, values in numeric.items():
        if not np.isfinite(values).all():
            raise ValueError(f'Nonfinite publication measurement: {name}')
    for role in ['primary', 'secondary']:
        if not numeric[f'{role}_ra_deg'].between(0, 360, inclusive='left').all():
            raise ValueError(f'Invalid {role} RA')
        if not numeric[f'{role}_dec_deg'].between(-90, 90).all():
            raise ValueError(f'Invalid {role} Dec')
    if not numeric['separation_arcsec'].gt(0).all():
        raise ValueError('Angular separations must be positive')
    result['original_system_name'] = result.system_name
    result['system_name'] = coordinate_names(numeric['primary_ra_deg'], numeric['primary_dec_deg'])
    if not result.system_name.is_unique:
        raise ValueError('Coordinate names collide at the requested precision; inspect the IDs')
    flags = {name: boolean_flag(result, name) for name in [
        'selected_by_resolved', 'selected_by_unresolved', 'selected_by_imaging',
        'literature_confirmed_dual', 'literature_confirmed_lens',
        'literature_uncertain_classification', 'literature_possible_association']}
    if (flags['literature_confirmed_dual'] & flags['literature_confirmed_lens']).any():
        raise ValueError('Conflicting confirmed dual/lens classifications')
    routes = []
    for index, row in result.iterrows():
        selected = []
        if flags['selected_by_resolved'].loc[index]:
            selected.append('gaia_resolved')
        if flags['selected_by_unresolved'].loc[index]:
            selected.append('gaia_unresolved')
        if flags['selected_by_imaging'].loc[index]:
            if row.imaging_selection_route not in ['one_gaia', 'two_gaia']:
                raise ValueError(f'Unknown imaging route: {row.candidate_id}')
            selected.append('imaging_' + row.imaging_selection_route)
        if not selected:
            raise ValueError(f'No selecting route: {row.candidate_id}')
        routes.append(';'.join(selected))
    result['selection_route'] = routes
    result['literature_class'] = np.select(
        [flags['literature_confirmed_dual'], flags['literature_confirmed_lens'],
         flags['literature_uncertain_classification'] | flags['literature_possible_association']],
        ['dual_agn', 'lens', 'uncertain'], default='unclassified')
    result = result.sort_values(['primary_ra_deg', 'primary_dec_deg'],
                                key=pd.to_numeric, kind='stable')
    return result[list(COLUMN_DESCRIPTIONS)].reset_index(drop=True)


def latex_escape(value):
    replacements = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%',
                    '_': r'\_', '#': r'\#', '$': r'\$',
                    '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}',
                    '^': r'\textasciicircum{}'}
    return ''.join(replacements.get(character, character) for character in str(value))


def write_previews(catalogue, output_dir, preview_rows):
    output_dir.mkdir(parents=True, exist_ok=True)
    for tier, count_macro in [('main', r'\Nmain'), ('additional', r'\NimagingOne')]:
        selected = catalogue.loc[catalogue.tier.eq(tier)].head(preview_rows)
        lines = [r'\begin{table*}', r'\centering', r'\footnotesize',
                 r'\setlength{\tabcolsep}{4pt}',
                 r'\caption{Example entries from the ' + tier + ' tier. The complete machine-readable table contains '
                 + count_macro + r' systems. Column definitions are given in Table~\ref{tab:catalogue_columns}. '
                 r'Route codes are R (resolved Gaia neighbour), U (unresolved Gaia neighbour), '
                 r'I1 (one-Gaia imaging) and I2 (two-Gaia imaging); a plus sign marks selection by both searches.}',
                 r'\label{tab:catalogue_' + tier + r'}',
                 r'\begin{tabular}{@{}lclrrrrrl@{}}', r'\toprule',
                 r'Name & Tier & Route & RA & Dec & $\Delta\theta$ & $z$ & $r_p$ & Literature class \\',
                 r' & & & [deg] & [deg] & [$^{\prime\prime}$] & & [kpc] & \\', r'\midrule']
        for row in selected.itertuples(index=False):
            route = '+'.join(ROUTE_LABELS[name] for name in row.selection_route.split(';'))
            values = [latex_escape(row.system_name), tier, route,
                      f'{float(row.primary_ra_deg):.6f}', f'{float(row.primary_dec_deg):+.6f}',
                      f'{float(row.separation_arcsec):.3f}', f'{float(row.redshift):.3f}',
                      f'{float(row.projected_separation_kpc):.2f}',
                      {'dual_agn': 'Dual AGN', 'lens': 'Lens', 'uncertain': 'Uncertain',
                       'unclassified': 'Unclassified'}[row.literature_class]]
            lines.append(' & '.join(values) + r' \\')
        lines += [r'\bottomrule', r'\end{tabular}', r'\end{table*}']
        (output_dir / f'catalogue_{tier}_preview.tex').write_text('\n'.join(lines) + '\n')
    lines = [r'\begin{table*}', r'\centering', r'\footnotesize',
             r'\caption{Columns in the machine-readable candidate catalogues. Blank fields indicate unavailable measurements or inapplicable routes. '
             r'Names use the primary position; the numeric coordinates retain their original precision. Gaia identifiers must be read as text.}',
             r'\label{tab:catalogue_columns}',
             r'\begin{tabularx}{\linewidth}{@{}llX@{}}', r'\toprule',
             r'Column & Unit & Description \\', r'\midrule']
    for column, (unit, description) in COLUMN_DESCRIPTIONS.items():
        lines.append(latex_escape(column) + ' & ' + latex_escape(unit or '---') + ' & '
                     + latex_escape(description) + r' \\')
    lines += [r'\bottomrule', r'\end{tabularx}', r'\end{table*}']
    (output_dir / 'catalogue_columns.tex').write_text('\n'.join(lines) + '\n')


def export_catalogue(candidates_path, exclusions_path, output_dir, audit_dir,
                     latex_output_dir=None, preview_rows=10, euclid_config=None):
    if preview_rows < 1:
        raise ValueError('Preview row count must be positive')
    if candidates_path.resolve() == (output_dir / 'final_candidates.csv').resolve():
        raise ValueError('Publication export must not overwrite its automatic input')
    candidates = read_csv(candidates_path)
    exclusions = read_csv(exclusions_path) if exclusions_path else pd.DataFrame()
    retained, rejected = curate_candidates(candidates, exclusions)
    catalogue = build_publication_catalogue(retained)
    counts = catalogue.tier.value_counts().to_dict()
    output_dir.mkdir(parents=True, exist_ok=True)
    catalogue.to_csv(output_dir / 'final_candidates.csv', index=False)
    for tier in ['main', 'additional']:
        catalogue.loc[catalogue.tier.eq(tier)].to_csv(output_dir / f'{tier}_tier.csv', index=False)
    pd.DataFrame([{'column': column, 'unit': unit, 'description': description}
                  for column, (unit, description) in COLUMN_DESCRIPTIONS.items()]).to_csv(
                      output_dir / 'columns.csv', index=False)
    if latex_output_dir:
        write_previews(catalogue, latex_output_dir, preview_rows)
    audit_dir.mkdir(parents=True, exist_ok=True)
    rejected.to_csv(audit_dir / 'excluded_candidates.csv', index=False)
    summary = dict(input_candidates=len(candidates), excluded_candidates=len(rejected),
                   publication_candidates=len(catalogue), tiers=counts,
                   input_sha256=hashlib.sha256(candidates_path.read_bytes()).hexdigest(),
                   exclusion_sha256=hashlib.sha256(exclusions_path.read_bytes()).hexdigest() if exclusions_path else None,
                   excluded_ids=rejected.candidate_id.tolist(),
                   catalogue_sha256=hashlib.sha256((output_dir / 'final_candidates.csv').read_bytes()).hexdigest())
    if euclid_config:
        from recount_euclid import recount

        config = json.loads(euclid_config.read_text())
        required = {'wide_moc', 'deep_moc', 'grid_records', 'classifications'}
        if set(config) != required:
            raise ValueError(f'Euclid configuration requires exactly {sorted(required)}')
        paths = {}
        for name, value in config.items():
            path = Path(value).expanduser()
            paths[name] = path if path.is_absolute() else euclid_config.parent / path
        audit, euclid_counts = recount(output_dir / 'final_candidates.csv',
            paths['wide_moc'], paths['deep_moc'], paths['grid_records'],
            paths['classifications'], audit_dir / 'euclid')
        audit.to_csv(output_dir / 'euclid_coverage.csv', index=False)
        pd.DataFrame([{'quantity': name, 'count': value} for name, value in euclid_counts.items()]).to_csv(
            output_dir / 'euclid_counts.csv', index=False)
        summary['euclid'] = euclid_counts
    (audit_dir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (output_dir / 'README.md').write_text(
        '# Candidate catalogue\n\n'
        f'The publication catalogue contains {len(catalogue)} unique systems: '
        f'{counts.get("main", 0)} in the main tier and {counts.get("additional", 0)} in the additional tier.\n\n'
        '- `final_candidates.csv`: complete curated catalogue.\n'
        '- `main_tier.csv` and `additional_tier.csv`: disjoint parts of that same catalogue.\n'
        '- `columns.csv`: field definitions and units.\n\n'
        'The main tier requires two distinct Gaia detections. The additional tier contains one-Gaia imaging selections. '
        'The original identifiers are retained in `original_system_name` and `candidate_id`; coordinate names are display labels. '
        'Read Gaia source IDs as text, not floating-point numbers.\n\n'
        f'This release applies {len(rejected)} explicit post-selection exclusion(s) to the automatic catalogue. '
        'The decisions and reasons are supplied to `export_catalogue.py` through its exclusions file. '
        'Automatic scores, selection rules and cutflows are unchanged. '
        'In this repository, `products/final_candidates.csv` is the automatic, pre-visual-inspection result; '
        'use the files in this folder for the publication sample.\n')
    if euclid_config:
        with (output_dir / 'README.md').open('a') as stream:
            stream.write('\n`euclid_coverage.csv` records footprint membership, VIS pixel availability '
                         'and the supplied visual classification for each publication candidate. '
                         '`euclid_counts.csv` gives the resulting counts. Footprint membership '
                         'does not guarantee usable imaging, and a resolved image pair does not '
                         'confirm two AGN.\n')
    print(f'{len(candidates)} automatic - {len(rejected)} exclusions = {len(catalogue)} publication candidates')
    print('Tier counts:', counts)
    print('Publication catalogues:', output_dir.resolve())
    return catalogue, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidates', type=Path, default=ROOT / 'products/final_candidates.csv')
    parser.add_argument('--exclusions', type=Path, help='Explicit post-selection decisions with exact Gaia component IDs')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'catalogue')
    parser.add_argument('--audit-dir', type=Path, default=ROOT / 'products_work/publication')
    parser.add_argument('--latex-output-dir', type=Path, help='Optional paper-only table fragments; no manuscript is edited')
    parser.add_argument('--preview-rows', type=int, default=10)
    parser.add_argument('--euclid-config', type=Path,
                        help='Optional local MOC, image-cache and visual-review paths for a Euclid recount')
    args = parser.parse_args()
    export_catalogue(args.candidates, args.exclusions, args.output_dir, args.audit_dir,
                     args.latex_output_dir, args.preview_rows, args.euclid_config)


if __name__ == '__main__':
    main()
