"""Select redshift uncertainty treatment from explicit catalogue metadata."""

import numpy as np
import pandas as pd


def photometric_redshifts(frame, catalogue_column, kind_column):
    if kind_column in frame:
        kind = frame[kind_column].astype('string').str.lower()
        if not kind.isin(['spectroscopic', 'photometric']).all():
            raise ValueError(f'{kind_column} must be spectroscopic or photometric for every row')
        return kind.eq('photometric')
    # Historical tables predate the explicit type column.
    return frame[catalogue_column].astype('string').str.contains('quaia', case=False, na=False)


def redshift_errors(frame, error_column, catalogue_column, kind_column):
    photometric = photometric_redshifts(frame, catalogue_column, kind_column)
    errors = pd.to_numeric(frame[error_column], errors='coerce')
    if (photometric & (~np.isfinite(errors) | errors.lt(0))).any():
        raise ValueError('Photometric redshifts require finite nonnegative published uncertainties')
    return errors.where(photometric, 0.)
