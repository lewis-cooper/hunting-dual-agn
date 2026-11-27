"""Assign catalogue priority tiers from exact Gaia component identities."""

import numpy as np
import pandas as pd


def boolean_flag(frame, column):
    """Read boolean and CSV-string flags without treating 'False' as true."""
    values = frame[column].astype('string')
    if not values.isin(['True', 'False']).all():
        raise ValueError(f'{column} must contain explicit True/False flags')
    return values.eq('True')


def assign_candidate_tiers(candidates):
    """Return a copy with main/additional labels; do not change selection."""
    result = candidates.copy()
    source_ids = {}
    for role in ('primary', 'secondary'):
        column = f'{role}_gaia_source_id'
        values = result[column].astype('string').fillna('')
        valid = values.str.fullmatch('[1-9][0-9]*', na=False)
        if (values.ne('') & ~valid).any():
            raise ValueError(f'{column} contains invalid or inexact Gaia identifiers')
        source_ids[role] = values
    primary_present = source_ids['primary'].ne('')
    secondary_present = source_ids['secondary'].ne('')
    distinct_sources = source_ids['primary'].ne(source_ids['secondary'])
    main = primary_present & secondary_present & distinct_sources
    if 'has_two_distinct_gaia_sources' in result:
        if not boolean_flag(result, 'has_two_distinct_gaia_sources').eq(main).all():
            raise ValueError('Gaia multiplicity flag disagrees with component identifiers')
    additional = (~main & primary_present
                  & boolean_flag(result, 'selected_by_imaging')
                  & result.imaging_selection_route.eq('one_gaia'))
    if not (main | additional).all():
        raise ValueError('Every candidate must have two Gaia sources or a one-Gaia imaging selection')
    for column in ('selected_by_resolved', 'selected_by_unresolved'):
        if column in result and (boolean_flag(result, column) & ~main).any():
            raise ValueError(f'{column} requires two distinct Gaia sources')
    result['tier'] = np.where(main, 'main', 'additional')
    return result
