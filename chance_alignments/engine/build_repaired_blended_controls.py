#!/usr/bin/env python3
"""Match isolated quasars and independent field sources to each blended input."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from build_repaired_component_controls import select_sources
from control_catalogues import read_csv
from prepare_matched_parent_controls import catalogue_family


def augment_background_fluxes(library, path):
    library = library.copy()
    if not path.exists():
        raise FileNotFoundError(f'Complete the all-morphology native-flux query first: {path}')
    native = pd.read_parquet(path)
    if not len(native):
        return library
    native = native.loc[native.unique_gaia_association & native.unique_native_association].set_index('target_id')
    for column in ['photometric_system', 'release', 'type', 'maskbits', 'fracflux_r',
                   *[f'{prefix}_{band}' for prefix in ['flux', 'flux_ivar', 'nobs'] for band in 'grz']]:
        if column not in native:
            continue
        replacement = library.background_source_id.map(native[column])
        use = library.photometric_system.isna() if column == 'photometric_system' else library.background_source_id.isin(native.index)
        library.loc[use, column] = replacement.loc[use]
    return library


def prepare(args):
    os.environ.setdefault('MPLCONFIGDIR',str(args.output/'scorer_runtime/matplotlib'))
    from dustmaps.sfd import SFDQuery
    root = args.output
    root.mkdir(parents=True, exist_ok=True)
    inputs = read_csv(root/'input_pairs.csv')
    real = inputs.loc[inputs.route.eq('blended') & inputs.latitude_bin.notna()].copy()
    if args.input_ids:
        requested=read_csv(args.input_ids).input_pair_id
        real=real.loc[real.pair_id.isin(requested)].copy()
    if 'band_policy' not in real:
        production = read_csv(args.repo/'products/one_arcsecond_scores.csv').set_index('qso_group_id')
        real['band_policy'] = real.pair_id.str.removeprefix('blended:').map(production.blend_extra_surveys_used).fillna('')
    real.loc[real.photometry_class.eq('legacy'), 'band_policy'] = 'legacy'
    for frame in [real]:
        sky = SkyCoord(frame.ra.to_numpy(float)*u.deg, frame.dec.to_numpy(float)*u.deg)
        frame['galactic_b_deg'] = sky.galactic.b.deg
        frame['native_region'] = np.where(frame.dec.ge(32.375) & frame.galactic_b_deg.gt(0), 'north', 'south')
        frame['ebv'] = SFDQuery(map_dir=str(args.dustmaps/'sfd'))(sky)
        frame['source_family'] = catalogue_family(frame.redshift_catalogue)

    primaries = pd.read_parquet(args.source/'parent_primary_photometry.parquet')
    primaries = primaries.loc[primaries.primary_isolated].copy()
    sky = SkyCoord(primaries.ra.to_numpy(float)*u.deg, primaries.dec.to_numpy(float)*u.deg)
    primaries['galactic_b_deg'] = sky.galactic.b.deg
    primaries['native_region'] = np.where(primaries.dec.ge(32.375) & primaries.galactic_b_deg.gt(0), 'north', 'south')
    primaries['ebv'] = SFDQuery(map_dir=str(args.dustmaps/'sfd'))(sky)
    primaries['source_family'] = catalogue_family(primaries.redshift_catalogue)
    primaries.to_parquet(root/'isolated_primary_library.parquet', index=False, compression='zstd')
    library = pd.read_parquet(root/'background_source_library.parquet')
    library = augment_background_fluxes(library.loc[library.route.eq('blended')],
                                        args.supplemental_photometry or args.source/'supplemental_legacy_photometry.parquet')
    library.to_parquet(root/'blended_background_library.parquet', index=False, compression='zstd')
    rng = np.random.default_rng(args.seed)
    trials = []; audit = []
    for template in real.itertuples():
        eligible = primaries.native_region.eq(template.native_region)
        eligible &= primaries.source_family.eq(template.source_family)
        eligible &= (primaries.primary_gaia_g-template.primary_gaia_g).abs().le(args.primary_magnitude_tolerance)
        eligible &= (primaries.redshift-template.redshift).abs().le(args.redshift_tolerance)
        eligible &= (primaries.galactic_b_deg.abs()-abs(template.galactic_b_deg)).abs().le(args.latitude_tolerance)
        eligible &= (primaries.ebv-template.ebv).abs().le(args.ebv_tolerance)
        if template.photometry_class == 'legacy':
            eligible &= primaries.primary_unique & primaries.primary_photometric_system.notna()
            eligible &= primaries.primary_flux_ivar_r.gt(0) & np.isfinite(primaries.primary_flux_r)
        donors = primaries.loc[eligible].copy()
        if len(donors):
            distance = SkyCoord(template.ra*u.deg, template.dec*u.deg).separation(
                SkyCoord(donors.ra.to_numpy(float)*u.deg, donors.dec.to_numpy(float)*u.deg)).deg
            donors = donors.loc[(distance > 1/60) & (distance <= args.max_sky_distance)]
        for annulus in ['3-10', '10-30']:
            sources = select_sources(template, library.loc[library.annulus.eq(annulus)],
                args.max_sky_distance, args.latitude_tolerance, args.ebv_tolerance)
            audit.append(dict(input_pair_id=template.pair_id, annulus=annulus,
                              n_primary_donors=len(donors), n_background_sources=len(sources)))
            if not len(donors) or not len(sources):
                trials.append(dict(trial_id=f'blended:{annulus}:{template.pair_id}:unmatched',
                    input_pair_id=template.pair_id, route='blended', annulus=annulus,
                    latitude_bin=template.latitude_bin, mock_weight=1., diagnostic_only=True,
                    selection_mode='unmatched', known_primary_rejection=False,
                    primary_field_id='unmatched', background_field_id='unmatched',
                    background_source_id='', primary_source_id='', secondary_pm_pass=True,
                    photometry_class=template.photometry_class, band_policy=template.band_policy))
                continue
            background_probability = sources.source_weight.to_numpy(float)
            primary_probability = donors.sampling_weight.to_numpy(float)
            source_draws = rng.choice(len(sources), args.draws_per_input, replace=True,
                                     p=background_probability/background_probability.sum())
            primary_draws = rng.choice(len(donors), args.draws_per_input, replace=True,
                                      p=primary_probability/primary_probability.sum())
            for number, (source_index, donor_index) in enumerate(zip(source_draws, primary_draws)):
                source = sources.iloc[source_index]
                donor = donors.iloc[donor_index]
                # Exclude a donor's own field and accidental physical overlap.
                if source.background_field_id == donor.target_id or SkyCoord(
                    source.ra*u.deg, source.dec*u.deg).separation(SkyCoord(donor.ra*u.deg, donor.dec*u.deg)).arcsec <= 60:
                    candidate = donors.loc[donors.target_id.ne(source.background_field_id)]
                    if len(candidate):
                        distances = SkyCoord(source.ra*u.deg,source.dec*u.deg).separation(
                            SkyCoord(candidate.ra.to_numpy(float)*u.deg,candidate.dec.to_numpy(float)*u.deg)).arcsec
                        candidate = candidate.loc[distances>60]
                    if not len(candidate):
                        trials.append(dict(trial_id=f'blended:{annulus}:{template.pair_id}:{number}:unmatched',
                            input_pair_id=template.pair_id,route='blended',annulus=annulus,
                            latitude_bin=template.latitude_bin,mock_weight=1/args.draws_per_input,
                            diagnostic_only=True,selection_mode='unmatched',known_primary_rejection=False,
                            primary_field_id='unmatched',background_field_id='unmatched',background_source_id='',
                            primary_source_id='',secondary_pm_pass=True,photometry_class=template.photometry_class,
                            band_policy=template.band_policy))
                        continue
                    donor = candidate.sample(n=1, weights='sampling_weight', random_state=int(rng.integers(2**31))).iloc[0]
                trial = source.to_dict()
                trial.update(trial_id=f'blended:{annulus}:{template.pair_id}:{number}',
                    input_pair_id=template.pair_id, target_id=template.pair_id, route='blended',
                    primary_field_id=donor.target_id, primary_source_id=str(donor.source_id),
                    primary_gaia_source_id=str(donor.source_id), primary_ra=donor.ra, primary_dec=donor.dec,
                    primary_tractor_ra=donor.primary_ra, primary_tractor_dec=donor.primary_dec,
                    primary_gaia_g=donor.primary_gaia_g, primary_photometric_system=donor.primary_photometric_system,
                    primary_release=donor.primary_release, redshift=donor.redshift,
                    redshift_error=donor.redshift_error, redshift_catalogue=donor.redshift_catalogue,
                    redshift_match_residual=donor.redshift-template.redshift,
                    primary_g_match_residual=donor.primary_gaia_g-template.primary_gaia_g,
                    companion_g_match_residual=source.companion_mag-template.companion_mag,
                    latitude_bin=template.latitude_bin, mock_weight=1/args.draws_per_input,
                    secondary_pm_pass=not np.isfinite(source.pm_sig) or source.pm_sig < 3,
                    secondary_pm_missing=not np.isfinite(source.pm_sig), primary_flux_independent=True,
                    known_primary_rejection=False, diagnostic_only=False, selection_mode='sampled',
                    photometry_class=template.photometry_class, band_policy=template.band_policy)
                for band in 'grz':
                    trial[f'primary_flux_{band}'] = donor[f'primary_flux_{band}']
                    trial[f'primary_flux_ivar_{band}'] = donor[f'primary_flux_ivar_{band}']
                trials.append(trial)
    frame = pd.DataFrame(trials)
    np.testing.assert_allclose(frame.groupby(['input_pair_id', 'annulus']).mock_weight.sum(), 1.)
    frame.to_parquet(root/'blended_control_trials.parquet', index=False, compression='zstd')
    frame.to_csv(root/'isolated_blend_trials.csv', index=False)
    pd.DataFrame(audit).to_csv(root/'blended_matching_audit.csv', index=False)
    (root/'blend_repair_configuration.json').write_text(json.dumps(dict(
        draws_per_input=args.draws_per_input, seed=args.seed,
        primary_magnitude_tolerance_mag=args.primary_magnitude_tolerance,
        redshift_tolerance=args.redshift_tolerance, latitude_tolerance_deg=args.latitude_tolerance,
        max_sky_distance_deg=args.max_sky_distance, ebv_tolerance=args.ebv_tolerance,
        primary='isolated measured quasar; same source family; own redshift and uncertainty',
        photometry='mask donor surveys to reproduce actual input Legacy/fallback access pattern',
        independence='sources and quasar fields retain bootstrap IDs; repeated combinations are not independent'), indent=2)+'\n')
    print(frame.groupby('annulus').agg(mocks=('trial_id','size'), eligible_pm=('secondary_pm_pass','sum'),
                                     unmatched=('diagnostic_only','sum')).to_string(), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--dustmaps', type=Path, required=True)
    parser.add_argument('--draws-per-input', type=int, default=100)
    parser.add_argument('--primary-magnitude-tolerance', type=float, default=.5)
    parser.add_argument('--redshift-tolerance', type=float, default=.25)
    parser.add_argument('--max-sky-distance', type=float, default=40.)
    parser.add_argument('--latitude-tolerance', type=float, default=10.)
    parser.add_argument('--ebv-tolerance', type=float, default=.02)
    parser.add_argument('--seed', type=int, default=20261122)
    parser.add_argument('--input-ids',type=Path)
    parser.add_argument('--supplemental-photometry',type=Path)
    prepare(parser.parse_args())


if __name__ == '__main__':
    main()
