"""Check user-supplied file locations without running discovery or any network call."""

from pathlib import Path
import glob

from dual_agn.measurement_inputs import read_config, local_path
from dual_agn.workflow import selected_pipelines


def check_inputs(quasar_path, measurement_path, imaging_path, pipeline='both', stage='vetting',
                 one_arcsecond_search_config=None):
    selected = selected_pipelines(pipeline)
    missing = []
    checked = set()

    def read(path):
        path = Path(path)
        if not path.is_file():
            missing.append(str(path))
            return None
        return read_config(path)

    def require(config, value):
        path = local_path(config, value)
        matches = [Path(name) for name in glob.glob(str(path))]
        if not matches or not all(path.is_file() for path in matches):
            missing.append(str(path))
        else:
            checked.update(matches)

    imaging_discovery_only = stage == 'discovery' and pipeline == 'imaging'
    if one_arcsecond_search_config is not None and not imaging_discovery_only:
        from dual_agn.one_arcsecond_search import resolve_config, input_paths
        checked.update(input_paths(resolve_config(one_arcsecond_search_config)))
    quasar = None if imaging_discovery_only or one_arcsecond_search_config is not None else read(quasar_path)
    if quasar:
        for key in ('catalogue', 'gaia_catalogue', 'gaia_coverage'):
            if key not in quasar:
                raise ValueError(f'{quasar_path}: set {key}')
            require(quasar, quasar[key])
    measurements = None if stage == 'discovery' else read(measurement_path)
    if measurements:
        specifications = [measurements['gaia'], measurements['tractor']]
        if measurements.get('desi'):
            specifications.append(measurements['desi'])
        if 'gaia_neighbour' in selected:
            specifications.extend(measurements.get('surveys', {}).values())
        for specification in specifications:
            if specification.get('enabled') is False:
                continue
            for value in specification.get('paths', [specification.get('path')]):
                if not value:
                    raise ValueError('Every local source needs a path or paths setting')
                require(measurements, value)
            if specification.get('coverage'):
                require(measurements, specification['coverage'])
        if measurements.get('quaia_catalogue'):
            require(measurements, measurements['quaia_catalogue'])
        if 'imaging' in selected and not measurements.get('morphology_images', {}).get('path_template'):
            raise ValueError('Set morphology_images.path_template for supplied FITS images')
    if 'imaging' in selected:
        imaging = read(imaging_path)
        if imaging:
            require(imaging, imaging['master_catalogue'])
            for key in ('tractor_coverage', 'gaia_catalogue', 'gaia_coverage', 'ccd_catalogue', 'ccd_coverage'):
                if key not in imaging:
                    raise ValueError(f'{imaging_path}: set {key}')
                require(imaging, imaging[key])
            for value in imaging.get('tractor_catalogues', []):
                require(imaging, value)
            if imaging.get('quaia_catalogue'):
                require(imaging, imaging['quaia_catalogue'])
    if missing:
        raise FileNotFoundError('Supply these local inputs before running discovery:\n' + '\n'.join(sorted(set(missing))))
    print(f'Found {len(checked)} local catalogue/coverage files for {pipeline}. No discovery, scoring or download was run.')
    if 'imaging' in selected:
        print('Individual FITS cutouts are checked after scoring identifies the systems requiring morphology.')
    return checked
