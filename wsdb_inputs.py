"""Optional WSDB exports of raw survey rows for the ordinary local-file builders.

Only temporary target tables are written on the server. Connection details never
enter configuration, query summaries or output files. psycopg is loaded lazily.
"""

import hashlib
import json
import os
from pathlib import Path
import re
import shutil

import numpy as np
import pandas as pd

from imaging_search.catalogues import file_hash
from measurement_inputs import sky


class WsdbInputError(RuntimeError):
    """A database failure whose message is safe to include in a public run log."""


def source_tables(specification):
    if set(specification).difference({'tables', 'table', 'release', 'columns'}):
        raise ValueError('WSDB sources accept table, release and column mappings only')
    tables = specification.get('tables', [specification])
    if not tables:
        raise ValueError('Supply at least one WSDB source table')
    fields = None
    for table in tables:
        if set(table).difference({'table', 'release', 'columns'}):
            raise ValueError('WSDB table settings cannot contain credentials or SQL expressions')
        if not table.get('release'):
            raise ValueError('Record the survey release for every WSDB table')
        names = table.get('columns', {})
        if not names:
            raise ValueError('Map canonical fields to native WSDB column names')
        identifiers = [*table.get('table', '').split('.'), *names, *names.values()]
        if not all(re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*', value) for value in identifiers):
            raise ValueError('WSDB table/column mappings must be identifiers, not SQL expressions')
        if fields is not None and list(names) != fields:
            raise ValueError('Multiple WSDB tables must map the same ordered canonical fields')
        fields = list(names)
    return tables


class WsdbCatalogues:
    def __init__(self, resume=False, connect=None, batch_size=2000):
        if int(batch_size) != batch_size or batch_size < 1:
            raise ValueError('WSDB batch_size must be a positive integer')
        self.resume = resume
        self.batch_size = int(batch_size)
        self._connect = connect
        self._connection = None

    def __enter__(self):
        return self

    def __exit__(self, *exception):
        if self._connection is not None:
            self._connection.close()

    def connection(self):
        if self._connection is None:
            try:
                if self._connect is None:
                    import psycopg
                    self._connect = psycopg.connect
                # libpq also supports PGSERVICE and an external password file.
                options = {key: os.environ[variable] for key, variable in {
                    'host': 'WSDB_HOST', 'port': 'WSDB_PORT', 'dbname': 'WSDB_DATABASE',
                    'user': 'WSDB_USER', 'password': 'WSDB_PASSWORD',
                }.items() if os.environ.get(variable)}
                if not os.environ.get('PGSERVICE'):
                    options.setdefault('host', os.environ.get('PGHOST', 'wsdb.ast.cam.ac.uk'))
                    options.setdefault('dbname', os.environ.get('WSDB_NAME', os.environ.get('PGDATABASE', 'wsdb')))
                self._connection = self._connect(connect_timeout=15,
                    application_name='hunting_dual_agn_inputs', **options)
            except ImportError:
                raise WsdbInputError('Install requirements-wsdb.txt to use WSDB; local-file runs do not need it') from None
            except Exception as error:
                raise WsdbInputError(f'WSDB connection failed ({type(error).__name__}). Check your private connection settings and network access.') from None
        return self._connection

    def check_sources(self, sources):
        """Check configured tables/columns before starting any catalogue exports."""
        try:
            connection = self.connection()
            from psycopg import sql
            with connection.transaction(), connection.cursor() as cursor:
                for specification in sources.values():
                    for table in source_tables(specification):
                        fields = sql.SQL(', ').join(map(sql.Identifier, table['columns'].values()))
                        cursor.execute(sql.SQL('SELECT {} FROM {} LIMIT 0').format(
                            fields, sql.Identifier(*table['table'].split('.'))))
        except WsdbInputError:
            raise
        except Exception as error:
            raise WsdbInputError(f'WSDB schema check failed ({type(error).__name__}); check table access and configured column names.') from None

    def export_cones(self, specification, targets, path, radius_arcsec):
        """Fetch every source in each cone; an empty successful cone has coverage."""
        sky(targets)
        if not np.isfinite(radius_arcsec) or radius_arcsec <= 0:
            raise ValueError('Query radius must be finite and positive in arcsec')
        for table in source_tables(specification):
            if not {'ra', 'dec'}.issubset(table['columns']):
                raise ValueError('Spatial WSDB inputs need canonical ra and dec columns')
        coverage = targets.drop(columns=['ra', 'dec']).copy()
        coverage['ra_deg'] = targets.ra.to_numpy()
        coverage['dec_deg'] = targets.dec.to_numpy()
        coverage['query_radius_arcsec'] = radius_arcsec
        return self._export(specification, targets, Path(path), coverage,
                            radius_arcsec=radius_arcsec)

    def export_ids(self, specification, identifiers, path, identity='targetid'):
        targets = pd.DataFrame({'identifier': pd.Series(identifiers, dtype='string')}).drop_duplicates()
        if targets.identifier.isna().any():
            raise ValueError('WSDB identity queries require nonmissing identifiers')
        if identity == 'targetid' and not targets.identifier.str.fullmatch('[1-9][0-9]*').all():
            raise ValueError('DESI target IDs must be exact positive integers')
        for table in source_tables(specification):
            if identity not in table['columns']:
                raise ValueError(f'WSDB identity query requires canonical {identity}')
        return self._export(specification, targets, Path(path), None, identity=identity)

    def _export(self, specification, targets, path, coverage, radius_arcsec=None, identity=None):
        tables = source_tables(specification)
        target_bytes = targets.to_csv(index=False).encode()
        fingerprint = {'source': specification, 'target_sha256': hashlib.sha256(target_bytes).hexdigest(),
                       'radius_arcsec': radius_arcsec, 'identity': identity,
                       'exporter_sha256': file_hash(Path(__file__))}
        summary_path = path.with_suffix('.query.json')
        coverage_path = path.with_name(path.stem + '_coverage.csv')
        if self.resume and summary_path.is_file() and path.is_file():
            previous = json.loads(summary_path.read_text())
            intact = previous.get('fingerprint') == fingerprint and previous.get('sha256') == file_hash(path)
            if coverage is not None:
                intact &= coverage_path.is_file() and previous.get('coverage_sha256') == file_hash(coverage_path)
            if intact:
                print('Reusing verified raw export:', path.name, flush=True)
                return path, coverage_path if coverage is not None else None
        path.parent.mkdir(parents=True, exist_ok=True)
        if len(targets) > self.batch_size:
            return self._export_batches(specification, targets, path, coverage, fingerprint,
                                        radius_arcsec, identity)
        temporary = path.with_suffix('.partial')
        columns = list(tables[0]['columns'])
        rows = 0
        try:
            if targets.empty:
                pd.DataFrame(columns=columns).to_csv(temporary, index=False)
            else:
                connection = self.connection()
                from psycopg import sql
                with connection.transaction(), connection.cursor() as cursor:
                    if identity is None:
                        cursor.execute('CREATE TEMP TABLE qso_input_targets (ra double precision, dec double precision) ON COMMIT DROP')
                        with cursor.copy('COPY qso_input_targets (ra, dec) FROM STDIN') as copied:
                            for row in targets[['ra', 'dec']].itertuples(index=False, name=None):
                                copied.write_row(row)
                    else:
                        cursor.execute('CREATE TEMP TABLE qso_input_targets (identifier text) ON COMMIT DROP')
                        with cursor.copy('COPY qso_input_targets (identifier) FROM STDIN') as copied:
                            for value in targets.identifier:
                                copied.write_row((value,))
                    cursor.execute('ANALYZE qso_input_targets')
                    selects = []
                    for table in tables:
                        selected = sql.SQL(', ').join(sql.SQL('s.{} AS {}').format(
                            sql.Identifier(native), sql.Identifier(canonical))
                            for canonical, native in table['columns'].items())
                        if identity is None:
                            condition = sql.SQL('q3c_join(t.ra, t.dec, s.{}, s.{}, {})').format(
                                sql.Identifier(table['columns']['ra']), sql.Identifier(table['columns']['dec']),
                                sql.Literal(radius_arcsec / 3600.))
                        else:
                            # Cast the small target list, leaving the indexed native
                            # DESI column unchanged, as in the original ID query.
                            expression = 's.{} = t.identifier::bigint' if identity == 'targetid' else 's.{}::text = t.identifier'
                            condition = sql.SQL(expression).format(sql.Identifier(table['columns'][identity]))
                        selects.append(sql.SQL('SELECT {} FROM {} AS s JOIN qso_input_targets AS t ON {}').format(
                            selected, sql.Identifier(*table['table'].split('.')), condition))
                    # One native record may occur in several overlapping query cones.
                    query = sql.SQL('COPY ({}) TO STDOUT WITH (FORMAT CSV, HEADER)').format(sql.SQL(' UNION ').join(selects))
                    if len(selects) == 1:
                        query = sql.SQL('COPY (SELECT DISTINCT * FROM ({}) AS matched) TO STDOUT WITH (FORMAT CSV, HEADER)').format(selects[0])
                    with cursor.copy(query) as copied, temporary.open('wb') as stream:
                        for block in copied:
                            stream.write(block)
            rows = sum(len(chunk) for chunk in pd.read_csv(temporary, chunksize=100_000, dtype='string'))
            temporary.replace(path)
            if coverage is not None:
                coverage.to_csv(coverage_path, index=False)
            summary = {'fingerprint': fingerprint, 'queried_fields': len(targets), 'source_rows': rows,
                       'sha256': file_hash(path)}
            if coverage is not None:
                summary['coverage_sha256'] = file_hash(coverage_path)
            summary_path.write_text(json.dumps(summary, indent=2) + '\n')
        except WsdbInputError:
            temporary.unlink(missing_ok=True)
            raise
        except Exception as error:
            temporary.unlink(missing_ok=True)
            raise WsdbInputError(f'Raw export {path.name} failed ({type(error).__name__}); no successful coverage was recorded. Check the configured table and column mappings.') from None
        print(f'Raw export {path.name}: {len(targets):,} fields; {rows:,} source rows', flush=True)
        return path, coverage_path if coverage is not None else None

    def _export_batches(self, specification, targets, path, coverage, fingerprint, radius_arcsec, identity):
        """Checkpoint the same cone query in bounded batches, as in the original scripts."""
        chunks = path.parent/(path.stem+'_query_chunks')
        chunks.mkdir(exist_ok=True)
        temporary = path.with_suffix('.partial')
        source_rows = 0
        try:
            with temporary.open('wb') as destination:
                for start in range(0, len(targets), self.batch_size):
                    stop = min(start+self.batch_size, len(targets))
                    chunk_path = chunks/f'{start:09d}_{stop:09d}.csv'
                    chunk_coverage = coverage.iloc[start:stop].reset_index(drop=True) if coverage is not None else None
                    self._export(specification, targets.iloc[start:stop].reset_index(drop=True), chunk_path,
                                 chunk_coverage, radius_arcsec, identity)
                    with chunk_path.open('rb') as source:
                        if start:
                            source.readline()
                        shutil.copyfileobj(source, destination)
                    summary = json.loads(chunk_path.with_suffix('.query.json').read_text())
                    source_rows += summary['source_rows']
                    print(f'{path.name}: {stop:,}/{len(targets):,} query fields', flush=True)
            temporary.replace(path)
            coverage_path = path.with_name(path.stem+'_coverage.csv')
            summary = {'fingerprint':fingerprint,'queried_fields':len(targets),
                       'source_rows':source_rows,'sha256':file_hash(path),'batch_size':self.batch_size,
                       'note':'Overlapping query batches can repeat native rows; the shared builders deduplicate exact native records.'}
            if coverage is not None:
                coverage.to_csv(coverage_path,index=False)
                summary['coverage_sha256'] = file_hash(coverage_path)
            path.with_suffix('.query.json').write_text(json.dumps(summary,indent=2)+'\n')
        finally:
            temporary.unlink(missing_ok=True)
        return path, coverage_path if coverage is not None else None
