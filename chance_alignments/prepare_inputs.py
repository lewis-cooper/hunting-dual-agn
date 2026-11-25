"""Prepare chance-alignment control inputs from the vetting tables."""

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from redshift_inputs import redshift_errors

ENGINE = Path(__file__).resolve().parent/'engine'
sys.path.insert(0, str(ENGINE))
from control_catalogues import ROUTES, latitude_bin, read_csv, read_fits


def prepare_control_inputs(products, input_root, work_root, quaia_path, outdir, imaging_input_root=None):
    """Use input pairs before colour/PM filtering, never the final sample as N_input."""
    one = read_csv(products/'one_arcsecond_scores.csv')
    components = read_fits(work_root/'one_arcsecond_inputs/one_arcsecond_components.fits',
                           ['component_gaia_source_id', 'gaia_phot_g_mean_mag'])
    components['component_gaia_source_id'] = components.component_gaia_source_id.astype('string')
    gaia_mag = components.drop_duplicates('component_gaia_source_id').set_index(
        'component_gaia_source_id').gaia_phot_g_mean_mag
    frames = []
    for route, mask in [('resolved', one.scoring_route.str.startswith('resolved')),
                        ('blended', one.scoring_route.eq('unresolved'))]:
        group = one.loc[mask]
        pairs = pd.DataFrame({
            'pair_id': route+':'+group.qso_group_id.astype(str), 'route': route,
            'ra': group.primary_component_ra, 'dec': group.primary_component_dec,
            'redshift': group.primary_qso_z, 'redshift_error': group.primary_qso_zerr,
            'redshift_catalogue': group.primary_qso_catalog,
            'primary_gaia_g': group.primary_component_gaia_source_id.astype('string').map(gaia_mag),
            'companion_mag': group.secondary_component_gaia_source_id.astype('string').map(gaia_mag),
            'magnitude_band': 'Gaia_G', 'primary_pm_sig': np.nan,
            'pair_separation_arcsec': group.component_separation_arcsec,
            # A measured unresolved system is NOT certified to be primary-only flux.
            'primary_flux_independent': route == 'resolved',
        })
        for band in 'grz':
            pairs[f'primary_flux_{band}'] = group[f'primary_decals_flux_{band}']
            pairs[f'primary_ivar_{band}'] = group[f'primary_decals_flux_ivar_{band}']
        pairs['primary_release'] = group.primary_decals_release
        pairs['primary_type'] = group.primary_decals_type
        pairs['redshift_error'] = redshift_errors(group, 'primary_qso_zerr', 'primary_qso_catalog', 'primary_qso_redshift_kind')
        pairs['photometry_class'] = 'not_applicable'
        pairs['band_policy'] = ''
        if route == 'blended':
            pairs['photometry_class'] = np.where(group.blend_decals_system_assignment.eq('none'),
                                                 'no_legacy', 'legacy')
            pairs['band_policy'] = group.blend_extra_surveys_used.fillna('')
            pairs.loc[pairs.photometry_class.eq('legacy'), 'band_policy'] = 'legacy'
        frames.append(pairs)
    imaging = read_csv(work_root/'imaging/labelled_systems_audit.csv')
    evidence = read_csv(products/'imaging_scores.csv')
    imaging_input_root = imaging_input_root or input_root
    components = read_fits(imaging_input_root/'imaging_components.fits')
    tractor = read_fits(imaging_input_root/'imaging_tractor.fits')
    primary = components.merge(tractor.drop_duplicates(['release','brickid','objid']),
        left_on=['tractor_release','tractor_brickid','tractor_objid'],
        right_on=['release','brickid','objid'], validate='many_to_one', suffixes=('','_source'))
    # Fixed primary roles can differ from the original r-band ranking.
    primary = primary.merge(imaging[['qso_parent_index','primary_ra','primary_dec']],
                            on='qso_parent_index', validate='many_to_one')
    offset = SkyCoord(primary.component_ra.to_numpy()*u.deg, primary.component_dec.to_numpy()*u.deg).separation(
        SkyCoord(primary.primary_ra.to_numpy()*u.deg, primary.primary_dec.to_numpy()*u.deg)).arcsec
    primary = primary.loc[offset < .0001].set_index('qso_parent_index')
    if not primary.index.is_unique:
        raise ValueError('Fixed primary assignment is not unique by imaging parent identifier')
    for two_gaia, route in [(True,'imaging_two_gaia'),(False,'imaging_one_gaia')]:
        mask = imaging.primary_has_gaia.eq(True) & imaging.secondary_has_distinct_gaia.eq(two_gaia)
        group = imaging.loc[mask].copy()
        pairs = pd.DataFrame({
            'pair_id': route+':'+group.qso_parent_index.astype(str), 'route': route,
            'ra': group.primary_ra, 'dec': group.primary_dec,
            'redshift': group.adopted_redshift, 'redshift_catalogue': group.adopted_redshift_catalogue,
            'redshift_error': 0., 'primary_gaia_g': group.primary_gaia_g,
            'companion_mag': group.secondary_r_mag, 'magnitude_band': 'Legacy_r',
            'primary_pm_sig': group.primary_pm_sig,
            'pair_separation_arcsec': group.component_sep_arcsec,
            'primary_flux_independent': True,
        })
        if 'redshift_kind' in group:
            pairs['redshift_error'] = redshift_errors(group, 'adopted_redshift_error', 'adopted_redshift_catalogue', 'redshift_kind')
        else:
            quaia = read_fits(quaia_path or imaging_input_root/'quaia_imaging_reference.fits',
                             ['source_id','redshift_quaia','redshift_quaia_err']).set_index('source_id')
            quaia.index = quaia.index.astype('string')
            errors = group.source_id.astype('string').map(quaia.redshift_quaia_err)
            uses_quaia = group.adopted_redshift_catalogue.str.contains('quaia', case=False)
            reference_z=group.source_id.astype('string').map(quaia.redshift_quaia)
            inconsistent=uses_quaia & np.isfinite(reference_z) & ~np.isclose(group.adopted_redshift,reference_z,rtol=0,atol=1e-12)
            if inconsistent.any():raise ValueError('Quaia reference redshifts differ from the adopted primary redshifts')
            pairs['redshift_error'] = errors.where(uses_quaia, 0.)
        for band in 'grz':
            pairs[f'primary_flux_{band}'] = group.qso_parent_index.map(primary[f'flux_{band}'])
            ivar = group.qso_parent_index.map(primary[f'flux_ivar_{band}'])
            usable = group.qso_parent_index.map(primary.maskbits).eq(0)
            usable &= group.qso_parent_index.map(primary[f'nobs_{band}']).gt(0)
            pairs[f'primary_ivar_{band}'] = ivar.where(usable, 0.)
        pairs['primary_release'] = group.qso_parent_index.map(primary.release)
        pairs['primary_type'] = group.qso_parent_index.map(primary.type)
        frames.append(pairs)
    pairs = pd.concat(frames, ignore_index=True)
    pairs['photometry_class'] = pairs.photometry_class.fillna('not_applicable')
    pairs['band_policy'] = pairs.band_policy.fillna('')
    pairs['abs_b_deg'], pairs['latitude_bin'] = latitude_bin(pairs.ra, pairs.dec)
    # Outside-training-latitude inputs remain in the audit, but are not scored.
    pairs.to_csv(outdir/'input_pairs.csv', index=False)
    # Empirical pair-only factors, separate from photometric and PM selection.
    factors = []
    evidence=evidence.copy()
    evidence['abs_b_deg'], evidence['latitude_bin'] = latitude_bin(evidence.primary_ra, evidence.primary_dec)
    imaging['abs_b_deg'], imaging['latitude_bin'] = latitude_bin(imaging.primary_ra, imaging.primary_dec)
    for route in ROUTES:
        for label in ['25-40','40-55','55-90']:
            dcr_num = dcr_den = morphology_num = morphology_den = 0
            dcr = morphology = 1.
            if route.startswith('imaging'):
                tier = 'two_gaia' if route.endswith('two_gaia') else 'one_gaia'
                rows = imaging.loc[imaging.latitude_bin.eq(label) & imaging.primary_has_gaia.eq(True)]
                rows = rows.loc[rows.secondary_has_distinct_gaia.eq(tier=='two_gaia')]
                before_dcr = np.isfinite(rows.primary_pm_sig) & rows.primary_pm_sig.lt(3)
                before_dcr &= np.isfinite(rows.delta_g_minus_r) & rows.delta_g_minus_r.abs().le(1)
                if tier == 'one_gaia':
                    dcr_den = int(before_dcr.sum())
                    dcr_num = int((before_dcr & rows.dcr_pass.eq(True)).sum())
                    dcr = dcr_num/dcr_den if dcr_den else np.nan
                rows = evidence.loc[evidence.candidate_tier.eq(tier) & evidence.latitude_bin.eq(label)]
                morphology_den = int(rows.joint_pre_morphology_pass.eq(True).sum())
                morphology_num = int(rows.joint_morphology_pass.eq(True).sum())
                morphology = morphology_num/morphology_den if morphology_den else np.nan
            factors.append(dict(route=route,latitude_bin=label,dcr_factor=dcr,
                morphology_factor=morphology,dcr_numerator=dcr_num,dcr_denominator=dcr_den,
                morphology_numerator=morphology_num,morphology_denominator=morphology_den))
    pd.DataFrame(factors).to_csv(outdir/'pair_cut_factors.csv', index=False)
    return pairs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['products', 'input-root', 'work-root', 'quaia', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--imaging-input-root',type=Path,
                        help='Generated imaging-search inputs; defaults to input-root for standalone compatibility')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    pairs = prepare_control_inputs(args.products, args.input_root, args.work_root, args.quaia, args.output,
                                  imaging_input_root=args.imaging_input_root)
    print('Control matching inputs:', pairs.groupby('route').size().to_dict())
    print('Input products:', args.products, '| Prepared controls:', args.output)
