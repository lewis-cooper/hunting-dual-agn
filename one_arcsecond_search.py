"""Reproduce the historical Gaia neighbour-count discovery from local catalogues."""

import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky


CATALOGUE_PRIORITY = {'SDSS_DR16Q': 0, 'DESI_DR1_QSO': 1, 'QUAIA_G20_V1': 2}
INPUT_MODES = ('local_catalogues', 'survey_gaia_matches', 'saved_qualifying_matches')
QUASAR_COLUMNS = ['qso_catalog', 'qso_uid', 'qso_objid', 'qso_name',
                  'qso_ra', 'qso_dec', 'qso_z', 'qso_zerr']
IDENTIFIER_COLUMNS = ['qso_uid', 'qso_objid', 'gaia_source_id', 'source_id',
                      'qso_parent_gaia_source_id', 'qso_group_id']
PARENT_FILENAME = 'discovered_parent.csv'


def file_hash(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read_catalogue(path):
    # Match the original discovery's standard pandas floating-point parsing.
    if Path(path).suffix.lower() in ('.fits', '.fit', '.fz'):
        from imaging_search.catalogues import read_frame
        return read_frame(Path(path))
    return pd.read_csv(path, dtype={name: 'string' for name in IDENTIFIER_COLUMNS},
                       low_memory=False)


def resolve_config(path):
    path = Path(path).resolve()
    config = json.loads(path.read_text())
    mode = config.get('input_mode', 'survey_gaia_matches')
    if mode not in INPUT_MODES:
        raise ValueError(f'input_mode must be one of {INPUT_MODES}')
    config['input_mode'] = mode
    file_keys = ['gaia_neighbour_counts']
    list_key = 'quasar_gaia_matches'
    if mode == 'local_catalogues':
        list_key = 'quasar_catalogues'
        file_keys += ['gaia_catalogue', 'gaia_coverage']
    values = config.get(list_key)
    if not isinstance(values, list) or not values:
        raise ValueError(f'Supply a nonempty {list_key} list')
    for name in file_keys:
        if not config.get(name):
            raise ValueError(f'Missing {name} in {path}')
        config[name] = (path.parent / Path(config[name]).expanduser()).resolve()
    config[list_key] = [(path.parent / Path(value).expanduser()).resolve() for value in values]
    for name in ('gaia_radius_arcsec', 'dedupe_radius_arcsec'):
        config[name] = float(config.get(name, 1.0))
        if not np.isfinite(config[name]) or config[name] <= 0:
            raise ValueError(f'{name} must be finite and positive')
    return config


def input_paths(config):
    paths = [config['gaia_neighbour_counts']]
    if config['input_mode'] == 'local_catalogues':
        paths += config['quasar_catalogues'] + [config['gaia_catalogue'], config['gaia_coverage']]
    else:
        paths += config['quasar_gaia_matches']
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise ValueError('Missing local discovery inputs:\n' + '\n'.join(missing))
    return paths


def require_columns(frame, columns, label):
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError(f'{label} lacks columns: {sorted(missing)}')


def require_positions(frame, ra_name, dec_name, label):
    ra_deg = pd.to_numeric(frame[ra_name], errors='coerce')
    dec_deg = pd.to_numeric(frame[dec_name], errors='coerce')
    valid = np.isfinite(ra_deg) & np.isfinite(dec_deg)
    valid &= ra_deg.ge(0) & ra_deg.lt(360) & dec_deg.abs().le(90)
    if not valid.all():
        raise ValueError(f'{label} has {int((~valid).sum())} invalid ICRS positions')


def require_gaia_ids(values, label):
    valid = values.notna() & values.str.fullmatch(r'[1-9][0-9]*', na=False)
    if not valid.all():
        raise ValueError(f'{label} requires exact, positive integer Gaia IDs (not floats)')
    values.astype('int64')


def read_quasars(paths):
    quasars = pd.concat([read_catalogue(path) for path in paths], ignore_index=True)
    require_columns(quasars, QUASAR_COLUMNS, 'Quasar catalogue')
    require_positions(quasars, 'qso_ra', 'qso_dec', 'Quasar catalogue')
    if quasars.qso_uid.isna().any() or quasars.qso_catalog.isna().any() or quasars.qso_catalog.astype(str).str.strip().eq('').any():
        raise ValueError('Every quasar needs a catalogue UID and a nonempty catalogue label')
    return quasars


def match_local_gaia(quasars, gaia, coverage, radius_arcsec):
    """All spatial matches for SDSS/DESI; exact published Gaia IDs for Quaia."""
    require_columns(gaia, ['source_id', 'ra', 'dec'], 'Gaia catalogue')
    require_gaia_ids(gaia.source_id, 'Gaia catalogue')
    require_positions(gaia, 'ra', 'dec', 'Gaia catalogue')
    if not gaia.source_id.is_unique:
        raise ValueError('Gaia catalogue source_id must be unique')
    positional = quasars.loc[~quasars.qso_catalog.eq('QUAIA_G20_V1')].reset_index(drop=True)
    require_columns(coverage, ['qso_uid', 'ra_deg', 'dec_deg', 'query_radius_arcsec'], 'Gaia coverage')
    if not coverage.qso_uid.is_unique:
        raise ValueError('Coverage must have one complete cone per positional quasar UID')
    checked = positional.merge(coverage, on='qso_uid', how='left', validate='many_to_one')
    finite = np.isfinite(checked[['ra_deg', 'dec_deg', 'query_radius_arcsec']].to_numpy(float)).all(axis=1)
    agrees = np.isclose(checked.qso_ra.to_numpy(float), checked.ra_deg.to_numpy(float), atol=1e-10, rtol=0)
    agrees &= np.isclose(checked.qso_dec.to_numpy(float), checked.dec_deg.to_numpy(float), atol=1e-10, rtol=0)
    complete = finite & agrees & checked.query_radius_arcsec.ge(radius_arcsec)
    if not complete.all():
        raise ValueError(f'{int((~complete).sum())} quasars lack complete local Gaia coverage')
    matches = []
    if len(positional) and len(gaia):
        qso_sky = SkyCoord(positional.qso_ra.to_numpy(float)*u.deg,
                           positional.qso_dec.to_numpy(float)*u.deg)
        gaia_sky = SkyCoord(gaia.ra.to_numpy(float)*u.deg, gaia.dec.to_numpy(float)*u.deg)
        qso_index, gaia_index, separations, _ = search_around_sky(qso_sky, gaia_sky, radius_arcsec*u.arcsec)
        frame = positional.iloc[qso_index].reset_index(drop=True)
        frame['gaia_source_id'] = gaia.source_id.iloc[gaia_index].to_numpy()
        frame['gaia_sep_arcsec'] = separations.arcsec
        matches.append(frame)
    quaia = quasars.loc[quasars.qso_catalog.eq('QUAIA_G20_V1')].copy()
    if len(quaia):
        require_columns(quaia, ['qso_parent_gaia_source_id'], 'Quaia catalogue')
        require_gaia_ids(quaia.qso_parent_gaia_source_id, 'Quaia catalogue')
        quaia['gaia_source_id'] = quaia.qso_parent_gaia_source_id
        quaia['gaia_sep_arcsec'] = 0.0
        matches.append(quaia)
    if not matches:
        return quasars.iloc[:0].assign(gaia_source_id=pd.Series(dtype='string'), gaia_sep_arcsec=np.nan)
    return pd.concat(matches, ignore_index=True)


def select_neighbour_matches(matches, counts, radius_arcsec):
    if 'qso_group_id' in matches or 'primary_qso_catalog' in matches:
        raise ValueError('Supply pre-grouping quasar/Gaia associations, not a previously selected parent table')
    require_columns(matches, QUASAR_COLUMNS + ['gaia_source_id', 'gaia_sep_arcsec'], 'Quasar/Gaia matches')
    require_positions(matches, 'qso_ra', 'qso_dec', 'Quasar/Gaia matches')
    require_gaia_ids(matches.gaia_source_id, 'Quasar/Gaia matches')
    if matches.qso_uid.isna().any() or matches.qso_catalog.isna().any() or matches.qso_catalog.astype(str).str.strip().eq('').any():
        raise ValueError('Missing quasar UID or empty catalogue label')
    require_columns(counts, ['gaia_source_id', 'gaia_edr3_cnt1'], 'EDR3 neighbour counts')
    require_gaia_ids(counts.gaia_source_id, 'EDR3 neighbour counts')
    if not counts.gaia_source_id.is_unique:
        raise ValueError('EDR3 neighbour counts must have unique Gaia source IDs')
    missing = ~matches.gaia_source_id.isin(counts.gaia_source_id)
    if missing.any():
        raise ValueError(f'{matches.loc[missing, "gaia_source_id"].nunique()} matched Gaia IDs lack neighbour-count coverage')
    # Counts, not duplicated columns in an association export, define selection.
    matches = matches.drop(columns=[name for name in counts if name != 'gaia_source_id' and name in matches])
    combined = matches.merge(counts, on='gaia_source_id', how='left', validate='many_to_one', sort=False)
    multiplicity = pd.to_numeric(combined.gaia_edr3_cnt1, errors='coerce')
    valid_count = np.isfinite(multiplicity) & multiplicity.ge(1) & multiplicity.mod(1).eq(0)
    if not valid_count.all():
        raise ValueError('Matched EDR3 cnt1 values must be finite positive integer counts including the central source')
    separation_arcsec = pd.to_numeric(combined.gaia_sep_arcsec, errors='coerce')
    valid_separation = np.isfinite(separation_arcsec) & separation_arcsec.ge(0)
    if not valid_separation.all():
        raise ValueError('Quasar/Gaia match separations must be finite and nonnegative')
    keep = multiplicity.gt(1) & separation_arcsec.le(radius_arcsec)
    combined = combined.loc[keep].copy()
    combined['_gaia_numeric_id'] = combined.gaia_source_id.astype('int64')
    combined = combined.sort_values(['qso_catalog', 'qso_uid', 'gaia_sep_arcsec', '_gaia_numeric_id'])
    combined = combined.drop(columns='_gaia_numeric_id')
    # Historical ordering: filter cnt1 first, then choose the nearest QUALIFYING match.
    nearest = combined.drop_duplicates('qso_uid', keep='first').reset_index(drop=True)
    return combined, nearest


def angular_sep_arcsec(ra1, dec1, ra2, dec2):
    ra1_rad, dec1_rad = math.radians(float(ra1)), math.radians(float(dec1))
    ra2_rad, dec2_rad = math.radians(float(ra2)), math.radians(float(dec2))
    value = (math.sin((dec2_rad-dec1_rad)/2)**2
             + math.cos(dec1_rad)*math.cos(dec2_rad)*math.sin((ra2_rad-ra1_rad)/2)**2)
    return math.degrees(2*math.asin(min(1.0, math.sqrt(max(0.0, value)))))*3600


def deduplicate_quasars(nearest, radius_arcsec):
    """Historical position groups and primary priorities; retain original group IDs."""
    nearest = nearest.reset_index(drop=True)
    catalogue_priority = dict(CATALOGUE_PRIORITY)
    for name in sorted(set(nearest.qso_catalog).difference(catalogue_priority)):
        catalogue_priority[name] = len(catalogue_priority)
    if nearest.empty:
        raise ValueError('No quasars have a qualifying Gaia companion in the supplied discovery data')
    parents = list(range(len(nearest)))

    def find(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    sky = SkyCoord(nearest.qso_ra.to_numpy(float)*u.deg, nearest.qso_dec.to_numpy(float)*u.deg)
    first, second, _, _ = search_around_sky(sky, sky, radius_arcsec*u.arcsec)
    for left, right in zip(first, second):
        parents[find(int(left))] = find(int(right))
    groups = defaultdict(list)
    for index in range(len(nearest)):
        groups[find(index)].append(index)
    ordered = sorted(groups.values(), key=lambda indices: (
        nearest.loc[indices, 'qso_ra'].min(), nearest.loc[indices, 'qso_dec'].min()))
    output = []
    for number, indices in enumerate(ordered, 1):
        group = nearest.loc[indices].copy()
        group['_priority'] = group.qso_catalog.map(catalogue_priority)
        group['_gaia_numeric_id'] = group.gaia_source_id.astype('int64')
        group = group.sort_values(['_priority', 'gaia_sep_arcsec', 'qso_uid', '_gaia_numeric_id'])
        primary = group.iloc[0]
        secondary = None if len(group) == 1 else group.iloc[1]
        catalogues = sorted(group.qso_catalog.unique(), key=catalogue_priority.get)
        row = dict(qso_group_id=f'QSO_CNT1_GT1_{number:06d}',
                   n_qso_catalog_rows=len(group), n_qso_catalogs=len(catalogues),
                   has_sdss='SDSS_DR16Q' in catalogues, has_desi='DESI_DR1_QSO' in catalogues,
                   has_quaia='QUAIA_G20_V1' in catalogues)
        parent_fields = ['catalog', 'uid', 'objid', 'name', 'ra', 'dec', 'z', 'zerr']
        if 'qso_redshift_kind' in group:
            parent_fields.append('redshift_kind')
        for field in parent_fields:
            row[f'primary_qso_{field}'] = primary[f'qso_{field}']
        for field in ('catalog', 'uid', 'objid', 'name', 'z', 'zerr'):
            row[f'secondary_qso_{field}'] = (
                ('NONE' if field not in ('z', 'zerr') else -999.0) if secondary is None
                else secondary[f'qso_{field}'])
        row['secondary_qso_sep_arcsec'] = -999.0 if secondary is None else angular_sep_arcsec(
            primary.qso_ra, primary.qso_dec, secondary.qso_ra, secondary.qso_dec)
        for field, name in [('catalog', 'catalogs'), ('uid', 'uids'), ('objid', 'objids'), ('z', 'redshifts')]:
            row[f'all_qso_{name}'] = ';'.join(group[f'qso_{field}'].fillna(-999.0).astype(str))
        row['all_gaia_source_ids'] = ';'.join(group.gaia_source_id)
        row['n_distinct_gaia_source_ids'] = group.gaia_source_id.nunique()
        row.update(primary.drop(labels=['_priority', '_gaia_numeric_id']).to_dict())
        output.append(row)
    parent = pd.DataFrame(output).rename(columns={'gaia_source_id': 'matched_gaia_source_id'})
    for name in parent:
        if parent[name].isna().any():
            parent[name] = parent[name].fillna(-999.0 if pd.api.types.is_numeric_dtype(parent[name]) else 'NONE')
    return parent


def validate_measurement_identities(parent, measurements):
    """Do not let shifted group numbers silently attach the wrong photometry."""
    require_columns(measurements, ['qso_group_id', 'primary_qso_catalog',
        'primary_qso_ra', 'primary_qso_dec', 'primary_component_gaia_source_id'], 'Associated measurements')
    if not measurements.qso_group_id.is_unique or not parent.qso_group_id.is_unique:
        raise ValueError('Discovery and measurement group IDs must be unique')
    if set(parent.qso_group_id) != set(measurements.qso_group_id):
        raise ValueError('Discovered membership differs from supplied measurement associations; rebuild those associations')
    matched = parent.set_index('qso_group_id').loc[measurements.qso_group_id]
    same_catalogue = matched.primary_qso_catalog.to_numpy() == measurements.primary_qso_catalog.to_numpy()
    same_id = matched.matched_gaia_source_id.to_numpy(str) == measurements.primary_component_gaia_source_id.to_numpy(str)
    same_ra = np.isclose(matched.primary_qso_ra, measurements.primary_qso_ra, atol=1e-10, rtol=0)
    same_dec = np.isclose(matched.primary_qso_dec, measurements.primary_qso_dec, atol=1e-10, rtol=0)
    consistent = same_catalogue & same_id & same_ra & same_dec
    if not consistent.all():
        raise ValueError(f'{int((~consistent).sum())} discovered identities disagree with measurement associations; rebuild them')


def run_search(config, output_root, resume=False):
    paths = input_paths(config)
    output_root = Path(output_root).resolve()
    if any(path.is_relative_to(output_root) for path in paths):
        raise ValueError('Discovery outputs must be separate from source inputs')
    fingerprint = dict(inputs={str(path): file_hash(path) for path in paths},
        code=file_hash(__file__), numpy=np.__version__, pandas=pd.__version__,
        input_mode=config['input_mode'], gaia_radius_arcsec=config['gaia_radius_arcsec'],
        dedupe_radius_arcsec=config['dedupe_radius_arcsec'])
    summary_path = output_root/'search_summary.json'
    parent_path = output_root/PARENT_FILENAME
    if resume and summary_path.is_file() and parent_path.is_file():
        summary = json.loads(summary_path.read_text())
        if (summary.get('fingerprint') == fingerprint
                and summary.get('output_sha256') == file_hash(parent_path)):
            print('Reused verified one-arcsecond discovery:', parent_path, flush=True)
            print(json.dumps(summary['counts'], indent=2), flush=True)
            return summary
    if config['input_mode'] == 'local_catalogues':
        quasars = read_quasars(config['quasar_catalogues'])
        matches = match_local_gaia(quasars, read_catalogue(config['gaia_catalogue']),
                                   read_catalogue(config['gaia_coverage']), config['gaia_radius_arcsec'])
        quasar_rows = len(quasars)
    else:
        matches = pd.concat([read_catalogue(path) for path in config['quasar_gaia_matches']], ignore_index=True)
        quasar_rows = None
    selected, nearest = select_neighbour_matches(matches, read_catalogue(config['gaia_neighbour_counts']),
                                                config['gaia_radius_arcsec'])
    parent = deduplicate_quasars(nearest, config['dedupe_radius_arcsec'])
    counts = dict(quasar_catalogue_rows=quasar_rows, gaia_match_rows=len(matches),
        qualifying_match_rows=len(selected), nearest_per_quasar_rows=len(nearest), parent_entries=len(parent),
        cnt1_distribution={str(key): int(value) for key, value in parent.gaia_edr3_cnt1.value_counts().sort_index().items()},
        primary_catalogue_distribution={key: int(value) for key, value in parent.primary_qso_catalog.value_counts().items()})
    output_root.mkdir(parents=True, exist_ok=True)
    parent.to_csv(parent_path, index=False)
    summary = dict(counts=counts, fingerprint=fingerprint, output_sha256=file_hash(parent_path),
        network_access=False, input_mode=config['input_mode'],
        note='Historical catalogue-position deduplication; parent entries are not guaranteed unique physical Gaia pairs.')
    summary_path.write_text(json.dumps(summary, indent=2)+'\n')
    print('One-arcsecond discovery inputs:', ', '.join(map(str, paths)), flush=True)
    print('Generated parent:', parent_path, flush=True)
    print(json.dumps(summary, indent=2), flush=True)
    return summary
