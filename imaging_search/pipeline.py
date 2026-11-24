"""Run local imaging discovery and build the inputs consumed by QSO vetting."""

import json
from pathlib import Path

import astropy
from astropy.table import Table
import numpy as np
import pandas as pd

from .catalogues import (decoded_frame, file_hash, find_tractor_pairs, input_paths, read_chunks, read_frame,
                         read_master, require_tractor_coverage, write_fits)
from .dcr import calculate_dcr
from .pairs import (build_candidate_summary, build_component_table, match_gaia_components,
                    relabel_by_qso_colour, require_gaia_coverage)


GENERATED_INPUTS = ['imaging_channel.fits','imaging_components.fits','imaging_labelled.fits',
                    'imaging_tractor.fits','quaia_imaging_reference.fits']


def discover(config, root):
    master = read_master(config['master_catalogue'])
    coverage = read_frame(config['tractor_coverage'])
    require_tractor_coverage(master,coverage,config['tractor_radius_arcsec'])
    channel, tractor, reference, summary = find_tractor_pairs(master,
        config['tractor_catalogues'],config['tractor_radius_arcsec'],config['aperture_arcsec'])
    if channel.empty:
        raise ValueError('No initial Tractor pairs were found in the supplied fields')
    # Retain the master's original provenance without loading its full-width table at once.
    parent_ids = set(channel.qso_parent_index)
    metadata = pd.concat([chunk.loc[chunk.qso_parent_index.isin(parent_ids)]
                          for chunk in read_chunks(config['master_catalogue'],float_precision='high')],ignore_index=True)
    metadata = metadata.set_index('qso_parent_index')
    for name in metadata.columns:
        if name not in channel:
            channel[name] = metadata.loc[channel.qso_parent_index,name].to_numpy()
    components = build_component_table(Table.from_pandas(channel),Table.from_pandas(tractor),
                                       config['aperture_arcsec'],max_targets=0)
    if not len(components):
        raise ValueError('No PSF pairs have two finite, positive r-band components')
    separation = np.asarray(components['component_sep_arcsec'],dtype=float)
    if not (np.isfinite(separation) & (separation > 0)).all():
        raise ValueError('Repeated or zero-separation component positions in initial pairs')
    inputs = root/'inputs'
    inputs.mkdir(parents=True,exist_ok=True)
    write_fits(channel,inputs/'imaging_channel.fits')
    write_fits(tractor,inputs/'imaging_tractor.fits')
    write_fits(components.to_pandas(),inputs/'imaging_components.fits')
    write_fits(reference,root/'nearest_colour_reference.fits')
    queries = components.to_pandas()[['component_id','qso_parent_index','tractor_release',
        'tractor_brickid','tractor_objid','component_ra','component_dec']].copy()
    queries['query_radius_arcsec'] = config['gaia_radius_arcsec']
    queries.to_csv(root/'gaia_query_positions.csv',index=False)
    primary = queries.loc[queries.component_id.str.endswith('_primary')].copy()
    primary = primary.rename(columns={'component_ra':'ra_deg','component_dec':'dec_deg'})
    primary[['qso_parent_index','ra_deg','dec_deg']].assign(query_radius_deg=.65).to_csv(root/'ccd_query_positions.csv',index=False)
    summary['validated_psf_pairs'] = len(components)//2
    summary['component_rows'] = len(components)
    print(json.dumps(summary,indent=2),flush=True)
    return summary


def associate(config, root, summary):
    inputs = root/'inputs'
    channel = Table.read(inputs/'imaging_channel.fits',memmap=False)
    components = Table.read(inputs/'imaging_components.fits',memmap=False)
    gaia_coverage = read_frame(config['gaia_coverage'])
    require_gaia_coverage(components,gaia_coverage,config['gaia_radius_arcsec'])
    gaia = read_frame(config['gaia_catalogue'],['source_id','ra','dec','phot_g_mean_mag',
        'pmra','pmdec','pmra_error','pmdec_error','pmra_pmdec_corr','parallax_over_error','ruwe'])
    matches = match_gaia_components(components,gaia,config['gaia_radius_arcsec'])
    component_summary = build_candidate_summary(components,matches)
    dcr = calculate_dcr(component_summary,read_frame(config['ccd_catalogue']),read_frame(config['ccd_coverage']))
    reference = read_frame(root/'nearest_colour_reference.fits')
    labelled, label_summary = relabel_by_qso_colour(component_summary,channel,
        reference=reference,z_tol=config['colour_redshift_half_width'])
    dcr_lookup = dcr.set_index('qso_parent_index')
    labelled['pair_dcr_delta_deg'] = pd.Series(np.asarray(labelled['qso_parent_index'])).map(dcr_lookup.pair_dcr_delta_deg).to_numpy()
    channel_lookup = decoded_frame(channel).set_index('qso_parent_index')
    parent_indices = np.asarray(labelled['qso_parent_index'])
    for name in channel.colnames:
        if name not in labelled.colnames:
            labelled[name] = channel_lookup.loc[parent_indices,name].to_numpy()
    labelled_frame = decoded_frame(labelled)
    labelled_frame['source_id'] = labelled_frame.source_id.astype('string')
    quaia = read_frame(config['quaia_catalogue'],['source_id','redshift_quaia','redshift_quaia_err'])
    quaia['source_id'] = quaia.source_id.astype('string')
    if not quaia.source_id.is_unique:
        raise ValueError('Quaia source_id must be unique')
    quaia = quaia.loc[quaia.source_id.isin(labelled_frame.source_id)].copy()
    lookup = quaia.set_index('source_id')
    uses_quaia = labelled_frame.adopted_redshift_catalogue.str.contains('quaia',case=False,na=False)
    values = labelled_frame.source_id.map(lookup.redshift_quaia)
    errors = labelled_frame.source_id.map(lookup.redshift_quaia_err)
    consistent = np.isfinite(errors) & errors.ge(0) & np.isclose(labelled_frame.adopted_redshift,values,rtol=0,atol=1e-12)
    if (uses_quaia & ~consistent).any():
        raise ValueError('Published Quaia errors/redshifts are missing or inconsistent for some discovered systems')
    write_fits(labelled_frame,inputs/'imaging_labelled.fits')
    write_fits(quaia,inputs/'quaia_imaging_reference.fits')
    dcr.to_csv(root/'dcr_axes.csv',index=False)
    summary.update(primary_has_gaia=int(labelled_frame.primary_has_gaia.sum()),
        one_gaia=int((labelled_frame.primary_has_gaia & ~labelled_frame.secondary_has_distinct_gaia).sum()),
        two_gaia=int((labelled_frame.primary_has_gaia & labelled_frame.secondary_has_distinct_gaia).sum()),
        systems_with_dcr=int(np.isfinite(labelled_frame.pair_dcr_delta_deg).sum()),
        fixed_primary_assignment=label_summary)
    return summary


def run_search(config, root, stage='all', resume=False):
    discovery_only = stage == 'discover'
    paths = input_paths(config,discovery_only)
    root = root.resolve()
    if any(path.is_relative_to(root) or root.is_relative_to(path) for path in paths):
        raise ValueError('Keep search outputs separate from their source catalogues')
    fingerprint = {'inputs':{str(path):file_hash(path) for path in paths},
        'code':{path.name:file_hash(path) for path in Path(__file__).parent.glob('*.py')},
        'numerical_versions':{'numpy':np.__version__,'pandas':pd.__version__,'astropy':astropy.__version__},
        'parameters':{name:config[name] for name in ['aperture_arcsec','tractor_radius_arcsec',
                                                    'gaia_radius_arcsec','colour_redshift_half_width']}}
    summary_path = root/'search_summary.json'
    if resume and summary_path.is_file():
        previous = json.loads(summary_path.read_text())
        if previous.get('stage') == stage and previous.get('fingerprint') == fingerprint:
            outputs = previous.get('output_sha256',{})
            if outputs and all((root/name).is_file() and file_hash(root/name)==expected for name,expected in outputs.items()):
                print('Reused verified local imaging search:',root,flush=True)
                print(json.dumps(previous['counts'],indent=2),flush=True)
                return previous
    root.mkdir(parents=True,exist_ok=True)
    counts = discover(config,root)
    if not discovery_only:
        counts = associate(config,root,counts)
    output_paths = list((root/'inputs').glob('*.fits'))+[root/'nearest_colour_reference.fits',
        root/'gaia_query_positions.csv',root/'ccd_query_positions.csv']
    if not discovery_only:
        output_paths.append(root/'dcr_axes.csv')
    summary = dict(stage=stage,counts=counts,fingerprint=fingerprint,
                   output_sha256={str(path.relative_to(root)):file_hash(path) for path in output_paths},
                   network_access=False,
                   note='Gaia and DCR use complete local coverage; PM, colour and DCR selection cuts are applied by the following vetting stages.')
    summary_path.write_text(json.dumps(summary,indent=2)+'\n')
    print('Search inputs:',config['master_catalogue'])
    print('Generated vetting inputs:',root/'inputs')
    print(json.dumps(counts,indent=2),flush=True)
    return summary
