"""Adapt a user quasar catalogue to both searches using configurable column maps.

This optional stage is for new catalogues. The historical discovery exports
remain available through the original configs, preserving their input definition.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import search_around_sky

from imaging_search.catalogues import read_frame, file_hash
from imaging_search.catalogues import resolve_config as imaging_config
from measurement_inputs import read_config, local_path, sky, nearest_sources
from settings import WORK_ROOT, ROOT


def singleton_json(value):
    values = json.loads(value)
    if not isinstance(values, list) or len(values) != 1:
        raise ValueError('A mapped JSON field requires exactly one value; ambiguous or missing values must be resolved in the input catalogue')
    return values[0]


def read_quasars(config):
    """Validate a mapped master before either local or WSDB data access."""
    columns = config['columns']
    required = {'object_id', 'ra', 'dec', 'redshift'}
    if not required.issubset(columns):
        raise ValueError('Map object_id, ra, dec and redshift in the quasar config')
    path = local_path(config, config['catalogue'])
    if path.suffix.lower() == '.csv':
        quasars = pd.read_csv(path, usecols=list(columns.values()),
            dtype={columns['object_id']:'string'}, float_precision='high')
    else:
        quasars = read_frame(path, list(columns.values()), float_precision='high')
    quasars = pd.DataFrame({name: quasars[native] for name, native in columns.items()})
    if not required.issubset(quasars):
        raise ValueError('The quasar input is missing mapped columns')
    quasars['object_id'] = quasars.object_id.astype('string')
    if not quasars.object_id.notna().all() or not quasars.object_id.is_unique:
        raise ValueError('Quasar object_id must be nonmissing and unique; choose a native row identifier')
    sky(quasars)
    if not np.isfinite(quasars.redshift).all():
        raise ValueError('Every input quasar requires a finite redshift')
    labels = quasars.get('catalogue_label', pd.Series(config.get('catalogue_label',''), index=quasars.index)).astype('string').str.strip()
    if not (labels.notna() & labels.ne('')).all():
        raise ValueError('Set a nonempty catalogue_label or map a catalogue_label column')
    if 'redshift_kind' in quasars:
        kinds = quasars.redshift_kind.astype('string').str.lower()
    elif 'redshift_kind_by_catalogue' in config:
        kinds = labels.map(config['redshift_kind_by_catalogue']).astype('string')
    else:
        kinds = pd.Series(config['redshift_kind'], index=quasars.index, dtype='string')
    if not kinds.isin(['spectroscopic', 'photometric']).all():
        raise ValueError('Redshift kind must be spectroscopic or photometric, as a config value or mapped column')
    if kinds.eq('photometric').any() and 'redshift_error' not in quasars:
        raise ValueError('Map redshift_error for photometric redshifts')
    quasars['redshift_error'] = quasars.get('redshift_error', 0.)
    if config.get('redshift_error_format', 'number') == 'json_singleton':
        photometric = kinds.eq('photometric')
        errors = pd.Series(0., index=quasars.index)
        errors.loc[photometric] = quasars.loc[photometric, 'redshift_error'].map(singleton_json).astype(float)
        quasars['redshift_error'] = errors
    valid_error = np.isfinite(quasars.redshift_error) & quasars.redshift_error.ge(0)
    if (kinds.eq('photometric') & ~valid_error).any():
        raise ValueError('Photometric redshift errors must be finite and nonnegative')
    quasars.loc[kinds.eq('spectroscopic'), 'redshift_error'] = 0.
    quasars['catalogue_object_id'] = quasars.object_id
    quasars['catalogue_object_ids_json'] = quasars.object_id.map(lambda value: json.dumps([str(value)]))
    for label, field in config.get('catalogue_object_id_columns', {}).items():
        selected = labels.eq(label)
        if field not in quasars:
            raise ValueError(f'Map the native object ID field {field} in columns')
        values = quasars.loc[selected, field]
        if config.get('catalogue_object_id_format') == 'json_array':
            identifiers = values.map(json.loads)
            if not identifiers.map(lambda value: isinstance(value,list) and len(value)>0).all():
                raise ValueError(f'{field} must contain a nonempty JSON array of native IDs')
            identifiers = identifiers.map(lambda values: list(dict.fromkeys(str(value) for value in values)))
            quasars.loc[selected, 'catalogue_object_ids_json'] = identifiers.map(json.dumps)
            single = identifiers.map(len).eq(1)
            quasars.loc[identifiers.index[single], 'catalogue_object_id'] = identifiers.loc[single].map(lambda values:values[0]).astype('string')
            continue
        if config.get('catalogue_object_id_format') == 'json_singleton':
            values = values.map(singleton_json)
        quasars.loc[selected, 'catalogue_object_id'] = values.astype('string')
        quasars.loc[selected, 'catalogue_object_ids_json'] = values.map(lambda value:json.dumps([str(value)]))
    if 'catalogue_labels' in config:
        labels = labels.map(config['catalogue_labels']).astype('string')
        if labels.isna().any():
            raise ValueError('catalogue_labels must map every native catalogue label')
    quasars['catalogue_label'] = labels
    quasars['redshift_kind'] = kinds
    quasars['qso_uid'] = labels + ':' + quasars.object_id
    return quasars


def adapt_catalogue(config, output, imaging_config_path=None, quasars=None):
    output.mkdir(parents=True, exist_ok=True)
    if quasars is None:
        quasars = read_quasars(config)
    labels, kinds = quasars.catalogue_label, quasars.redshift_kind
    gaia_path = local_path(config, config['gaia_catalogue'])
    gaia = read_frame(gaia_path, ['source_id','ra','dec','phot_g_mean_mag'])
    gaia['source_id'] = gaia.source_id.astype('string')
    gaia = gaia.drop_duplicates()
    if not gaia.source_id.is_unique or not gaia.source_id.str.fullmatch('[1-9][0-9]*').all():
        raise ValueError('Gaia requires unique exact source IDs')
    coverage_path = local_path(config, config['gaia_coverage'])
    coverage = pd.read_csv(coverage_path, dtype={'object_id':'string'})
    required_coverage = {'object_id','ra_deg','dec_deg','query_radius_arcsec'}
    if not required_coverage.issubset(coverage) or not coverage.object_id.is_unique:
        raise ValueError('Gaia coverage needs unique object_id, ra_deg, dec_deg, query_radius_arcsec')
    covered = coverage.set_index('object_id').reindex(quasars.object_id).reset_index()
    valid = np.isfinite(covered[['ra_deg','dec_deg','query_radius_arcsec']].to_numpy(float)).all(axis=1)
    offsets = np.full(len(quasars), np.inf)
    if valid.any():
        offsets[valid] = sky(covered.loc[valid], 'ra_deg', 'dec_deg').separation(sky(quasars.loc[valid])).arcsec
    # Two arcsec encloses every potential 1-arcsec anchor's 1-arcsec neighbours.
    complete = valid & covered.query_radius_arcsec.ge(2.) & (offsets < 1e-5)
    if not complete.all():
        raise ValueError(f'{(~complete).sum()} quasars lack complete 2-arcsec Gaia data')
    nearest = nearest_sources(quasars, gaia, 1., ['source_id'])
    normalized = pd.DataFrame(dict(qso_catalog=labels, qso_uid=quasars.qso_uid,
        qso_objid=quasars.catalogue_object_id, qso_name=quasars.get('name', quasars.object_id),
        qso_ra=quasars.ra, qso_dec=quasars.dec, qso_z=quasars.redshift,
        qso_zerr=quasars.redshift_error, qso_redshift_kind=kinds,
        qso_catalogue_object_ids_json=quasars.catalogue_object_ids_json))
    normalized.to_csv(output/'quasars.csv', index=False)
    generated_coverage = covered.drop(columns='object_id').copy()
    generated_coverage['qso_uid'] = quasars.qso_uid.to_numpy()
    generated_coverage.to_csv(output/'gaia_coverage.csv', index=False)
    matched = nearest.source_id.notna()
    # Imaging is anchored to the nearest Gaia source; keep this association
    # separate from the one-arcsecond all-match-then-cnt1 selection.
    master = pd.DataFrame(dict(source_id=nearest.loc[matched, 'source_id'].to_numpy(),
        gaia_ra=nearest.loc[matched, 'ra'].to_numpy(), gaia_dec=nearest.loc[matched, 'dec'].to_numpy(),
        phot_g_mean_mag=nearest.loc[matched, 'phot_g_mean_mag'].to_numpy(),
        adopted_redshift=quasars.loc[matched, 'redshift'].to_numpy(),
        adopted_redshift_error=quasars.loc[matched, 'redshift_error'].to_numpy(),
        adopted_redshift_catalogue=labels.loc[matched].to_numpy(), redshift_kind=kinds.loc[matched].to_numpy(),
        q_ra=quasars.loc[matched, 'ra'].to_numpy(), q_dec=quasars.loc[matched, 'dec'].to_numpy(),
        sdss_name=normalized.loc[matched, 'qso_name'].to_numpy(),
        input_object_id=quasars.loc[matched, 'object_id'].to_numpy(),
        anchor_match_sep_arcsec=nearest.loc[matched, 'match_sep_arcsec'].to_numpy()))
    master = master.sort_values(['anchor_match_sep_arcsec','input_object_id'], kind='stable').drop_duplicates('source_id')
    master = master.sort_values('source_id').reset_index(drop=True)
    master.insert(0, 'qso_parent_index', np.arange(len(master), dtype=np.int64))
    master.to_csv(output/'quasar_master.csv', index=False)
    # Counts cover ALL matched Gaia anchors, not only nearest anchors or candidates.
    _, matched_gaia, _, _ = search_around_sky(sky(quasars), sky(gaia), 1.*u.arcsec)
    anchors = gaia.iloc[np.unique(matched_gaia)].reset_index(drop=True)
    anchor_index, _, _, _ = search_around_sky(sky(anchors), sky(gaia), 1.*u.arcsec)
    counts = np.bincount(anchor_index, minlength=len(anchors))
    pd.DataFrame({'gaia_source_id':anchors.source_id, 'gaia_edr3_cnt1':counts}).to_csv(output/'gaia_neighbour_counts.csv', index=False)
    one_config = dict(input_mode='local_catalogues', quasar_catalogues=[str(output/'quasars.csv')],
        gaia_catalogue=str(gaia_path), gaia_coverage=str(output/'gaia_coverage.csv'),
        gaia_neighbour_counts=str(output/'gaia_neighbour_counts.csv'), gaia_radius_arcsec=1., dedupe_radius_arcsec=1.)
    (output/'one_arcsecond_search.json').write_text(json.dumps(one_config,indent=2)+'\n')
    if imaging_config_path is not None:
        imaging = imaging_config(imaging_config_path)
        imaging['master_catalogue'] = output/'quasar_master.csv'
        (output/'imaging_search.json').write_text(json.dumps(imaging, indent=2, default=str)+'\n')
    summary = dict(input_rows=len(quasars), matched_quasar_rows=int(matched.sum()),
                   unmatched_quasar_rows=int((~matched).sum()), imaging_gaia_anchors=len(master),
                   native_catalogue_sha256=file_hash(local_path(config, config['catalogue'])), catalogue_counts=labels.value_counts().to_dict(), redshift_kind_counts=kinds.value_counts().to_dict())
    (output/'catalogue_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))
    print('Generated catalogue/configs:', output)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--imaging-search-config', type=Path, default=ROOT/'imaging_search.local.json')
    parser.add_argument('--output-root', type=Path, default=WORK_ROOT/'quasar_catalogue')
    args = parser.parse_args()
    adapt_catalogue(read_config(args.config), args.output_root.resolve(), args.imaging_search_config)


if __name__ == '__main__':
    main()
