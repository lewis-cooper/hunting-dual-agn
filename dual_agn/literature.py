"""Match published systems while keeping finder associations distinct from pair identities."""

from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import pandas as pd


CONFIRMED_CONTAMINANTS = {
    'confirmed_lens', 'confirmed_quasar_star', 'confirmed_projected_quasar_pair',
    'confirmed_quasar_galaxy', 'confirmed_stellar_pair',
}
UNCERTAIN_CLASSIFICATIONS = {
    'likely_lens', 'lens_candidate', 'disputed_dual_lens',
    'possible_quasar_star', 'unclassified_quasar_pair',
}
REFERENCE_NUMERIC_COLUMNS = [
    'ra_deg', 'dec_deg', 'secondary_ra_deg', 'secondary_dec_deg', 'redshift',
    'secondary_redshift', 'separation_arcsec', 'position_angle_deg',
    'position_tolerance_arcsec', 'separation_tolerance_arcsec',
]
LITERATURE_FLAG_COLUMNS = [
    'literature_assessed', 'literature_confirmed_dual', 'literature_possible_association',
    'literature_uncertain_classification', 'literature_contaminant_removed',
]
LITERATURE_TEXT_COLUMNS = [
    'literature_reference_ids', 'literature_system_names', 'literature_classifications',
    'literature_author_year', 'literature_publication_urls', 'literature_exclusion_reason',
]
MATCH_COLUMNS = [
    'candidate_id', 'reference_id', 'system_id', 'system_name', 'classification',
    'reference_catalogue', 'author_year', 'publication_url', 'source_table',
    'match_offset_arcsec', 'maximum_pair_offset_arcsec', 'separation_difference_arcsec',
    'axis_difference_deg', 'verification_method', 'match_verified',
    'confirmed_contaminant', 'uncertain_classification', 'remove_candidate',
]


def read_references(binaries_path, contaminants_path):
    frames = []
    for path, catalogue in [(binaries_path, 'known_binaries'),
                            (contaminants_path, 'known_lenses_contaminants')]:
        frame = pd.read_csv(path, dtype='string').fillna('')
        required = {'reference_id', 'system_id', 'system_name', 'classification',
                    'primary_gaia_source_id', 'secondary_gaia_source_id',
                    'author_year', 'publication_url', 'source_table', 'match_rule',
                    *REFERENCE_NUMERIC_COLUMNS}
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f'{path}: missing reference columns {sorted(missing)}')
        for column in REFERENCE_NUMERIC_COLUMNS:
            numeric = pd.to_numeric(frame[column], errors='coerce')
            if (frame[column].ne('') & ~np.isfinite(numeric)).any():
                raise ValueError(f'{path}: invalid numeric values in {column}')
            frame[column] = numeric.astype(float)
        if frame.reference_id.eq('').any() or not frame.reference_id.is_unique:
            raise ValueError(f'{path}: reference_id must be nonempty and unique')
        for column in ['primary_gaia_source_id', 'secondary_gaia_source_id']:
            if (frame[column].ne('') & ~frame[column].str.fullmatch('[1-9][0-9]*')).any():
                raise ValueError(f'{path}: Gaia IDs must be exact integer strings in {column}')
        allowed_classes = {'confirmed_dual_agn'} if catalogue == 'known_binaries' else CONFIRMED_CONTAMINANTS | UNCERTAIN_CLASSIFICATIONS
        if not frame.classification.isin(allowed_classes).all():
            raise ValueError(f'{path}: unknown or inappropriate classification')
        if frame[['author_year', 'publication_url', 'source_table']].eq('').any().any():
            raise ValueError(f'{path}: every reference needs publication provenance')
        for ra_column, dec_column in [('ra_deg', 'dec_deg'),
                                      ('secondary_ra_deg', 'secondary_dec_deg')]:
            ra, dec = frame[ra_column], frame[dec_column]
            if (ra.notna() != dec.notna()).any():
                raise ValueError(f'{path}: incomplete coordinate pair')
            if (ra.notna() & ~(ra.ge(0) & ra.lt(360) & dec.between(-90, 90))).any():
                raise ValueError(f'{path}: coordinates outside degree ranges')
        if not frame.match_rule.isin(['position_only', 'component_pair', 'published_pair_geometry',
                                      'reviewed_system_geometry']).all():
            raise ValueError(f'{path}: unknown match_rule')
        if not frame[['position_tolerance_arcsec', 'separation_tolerance_arcsec']].gt(0).all().all():
            raise ValueError(f'{path}: matching tolerances must be finite and positive')
        has_secondary = frame.secondary_ra_deg.notna()
        if (frame.match_rule.eq('component_pair') & ~has_secondary).any():
            raise ValueError(f'{path}: component_pair rule requires both positions')
        if catalogue == 'known_binaries':
            eligible = frame.redshift.between(0, 4) & frame.separation_arcsec.gt(0) & frame.separation_arcsec.le(1.5)
            if not eligible.all():
                print(f'{path}: excluded {(~eligible).sum()} binary references outside 0 <= z <= 4 and 0 < separation <= 1.5 arcsec.')
            frame = frame.loc[eligible].copy()
        frame['reference_catalogue'] = catalogue
        frames.append(frame)
    references = pd.concat(frames, ignore_index=True)
    if not references.reference_id.is_unique:
        raise ValueError('Reference IDs overlap between the two catalogues')
    return references


def match_references(candidates, references, search_radius_arcsec=1.0):
    if not np.isfinite(search_radius_arcsec) or search_radius_arcsec <= 0:
        raise ValueError('The finder radius must be finite and positive')
    if candidates.candidate_id.isna().any() or not candidates.candidate_id.is_unique:
        raise ValueError('Candidate IDs must be non-missing and unique')
    coordinate_columns = ['primary_ra_deg', 'primary_dec_deg', 'secondary_ra_deg', 'secondary_dec_deg']
    if not np.isfinite(candidates[coordinate_columns].to_numpy(float)).all():
        raise ValueError('Candidate components must have finite degree coordinates')
    primary = SkyCoord(candidates.primary_ra_deg.to_numpy(float)*u.deg,
                       candidates.primary_dec_deg.to_numpy(float)*u.deg)
    secondary = SkyCoord(candidates.secondary_ra_deg.to_numpy(float)*u.deg,
                         candidates.secondary_dec_deg.to_numpy(float)*u.deg)
    candidate_pa_deg = primary.position_angle(secondary).deg
    midpoint = SkyCoord(primary.cartesian + secondary.cartesian)
    matches = []
    for reference in references.itertuples(index=False):
        if not np.isfinite([reference.ra_deg, reference.dec_deg]).all():
            continue
        first = SkyCoord(reference.ra_deg*u.deg, reference.dec_deg*u.deg)
        first_primary = first.separation(primary).arcsec
        first_secondary = first.separation(secondary).arcsec
        midpoint_offset = first.separation(midpoint).arcsec
        finder_offset = np.minimum(first_primary, first_secondary)
        pair_offset = np.full(len(candidates), np.nan)
        if np.isfinite([reference.secondary_ra_deg, reference.secondary_dec_deg]).all():
            second = SkyCoord(reference.secondary_ra_deg*u.deg, reference.secondary_dec_deg*u.deg)
            second_primary = second.separation(primary).arcsec
            second_secondary = second.separation(secondary).arcsec
            finder_offset = np.minimum(finder_offset, np.minimum(second_primary, second_secondary))
            pair_offset = np.minimum(np.maximum(first_primary, second_secondary),
                                     np.maximum(first_secondary, second_primary))
        for index in np.flatnonzero(finder_offset <= search_radius_arcsec):
            candidate = candidates.iloc[index]
            candidate_ids = [str(candidate.primary_gaia_source_id), str(candidate.secondary_gaia_source_id)]
            published_ids = [str(reference.primary_gaia_source_id), str(reference.secondary_gaia_source_id)]
            exact_gaia_pair = (all(value.isdigit() and int(value) > 0 for value in published_ids)
                               and len(set(published_ids)) == 2
                               and set(candidate_ids) == set(published_ids))
            separation_difference = abs(candidate.separation_arcsec - reference.separation_arcsec)
            axis_difference = (abs((candidate_pa_deg[index] - reference.position_angle_deg + 90) % 180 - 90)
                               if np.isfinite(reference.position_angle_deg) else np.nan)
            method = 'possible_positional_association'
            if exact_gaia_pair:
                method = 'exact_distinct_gaia_pair'
            elif np.isfinite(pair_offset[index]) and pair_offset[index] <= reference.position_tolerance_arcsec:
                method = 'published_component_pair'
            elif reference.match_rule in {'published_pair_geometry', 'reviewed_system_geometry'}:
                system_offset = min(first_primary[index], first_secondary[index])
                if reference.match_rule == 'reviewed_system_geometry':
                    system_offset = min(system_offset, midpoint_offset[index])
                primary_agrees = system_offset <= reference.position_tolerance_arcsec
                separation_agrees = np.isfinite(separation_difference) and separation_difference <= reference.separation_tolerance_arcsec
                angle_agrees = not np.isfinite(axis_difference) or axis_difference <= 5.0
                if primary_agrees and separation_agrees and angle_agrees:
                    method = reference.match_rule
            verified = method != 'possible_positional_association'
            confirmed_contaminant = reference.classification in CONFIRMED_CONTAMINANTS
            uncertain = reference.classification in UNCERTAIN_CLASSIFICATIONS
            matches.append({
                'candidate_id': candidate.candidate_id, 'reference_id': reference.reference_id,
                'system_id': reference.system_id, 'system_name': reference.system_name,
                'classification': reference.classification, 'reference_catalogue': reference.reference_catalogue,
                'author_year': reference.author_year, 'publication_url': reference.publication_url,
                'source_table': reference.source_table, 'match_offset_arcsec': finder_offset[index],
                'maximum_pair_offset_arcsec': pair_offset[index],
                'separation_difference_arcsec': separation_difference, 'axis_difference_deg': axis_difference,
                'verification_method': method, 'match_verified': verified,
                'confirmed_contaminant': confirmed_contaminant,
                'uncertain_classification': uncertain, 'remove_candidate': verified and confirmed_contaminant,
            })
    associations = pd.DataFrame(matches, columns=MATCH_COLUMNS)
    return associations.sort_values(['candidate_id', 'reference_id'], kind='stable').reset_index(drop=True)


def annotate_candidates(candidates, associations):
    result = candidates.drop(columns=[*LITERATURE_FLAG_COLUMNS, *LITERATURE_TEXT_COLUMNS], errors='ignore').copy()
    result['literature_assessed'] = True
    for column in LITERATURE_FLAG_COLUMNS[1:]:
        result[column] = False
    for column in LITERATURE_TEXT_COLUMNS:
        result[column] = ''
    for candidate_id, group in associations.groupby('candidate_id', sort=False):
        selected = result.candidate_id.eq(candidate_id)
        verified = group.loc[group.match_verified]
        confirmed_binary = verified.classification.eq('confirmed_dual_agn').any()
        excluded = verified.loc[verified.confirmed_contaminant]
        # Conflicting classifications require review; they never trigger automatic removal.
        disputed = verified.classification.eq('disputed_dual_lens').any()
        conflict = confirmed_binary and not excluded.empty
        remove = not excluded.empty and not disputed and not conflict
        result.loc[selected, 'literature_confirmed_dual'] = confirmed_binary and not conflict and not disputed
        result.loc[selected, 'literature_possible_association'] = (~group.match_verified).any()
        result.loc[selected, 'literature_uncertain_classification'] = group.uncertain_classification.any() or conflict or disputed
        result.loc[selected, 'literature_contaminant_removed'] = remove
        for destination, source in [
            ('literature_reference_ids', 'reference_id'), ('literature_system_names', 'system_name'),
            ('literature_classifications', 'classification'), ('literature_author_year', 'author_year'),
            ('literature_publication_urls', 'publication_url'),
        ]:
            result.loc[selected, destination] = '; '.join(dict.fromkeys(group[source].astype(str)))
        if remove:
            result.loc[selected, 'literature_exclusion_reason'] = '; '.join(sorted(set(excluded.classification)))
    return result
