#!/usr/bin/env python3
"""Report method false-positive pass rates without a selected-candidate comparison."""

import argparse
from pathlib import Path

import pandas as pd


def write_report(root):
    summary = pd.read_csv(root/'repaired_summary_by_route.csv')
    main = summary.loc[summary.annulus.eq('3-10')]
    lines = ['# Random-control acceptance by the selection method', '',
        'The measured quantity is P(pass | unrelated matched control, known outcome). '
        'It describes how often unrelated sources pass this method for the supplied input '
        'brightness, redshift and sky distribution. It is not the fraction of selected '
        'systems that are unrelated. No final-candidate catalogue is read.', '',
        '| Route | Annulus [arcsec] | Pre-vetting inputs | Distinct background sources | '
        'Pass rate (95% interval) | Relative 95% uncertainty | Known outcome weight |',
        '| --- | --- | ---: | ---: | --- | ---: | ---: |']
    for row in summary.itertuples():
        lines.append(f'| {row.route} | {row.annulus} | {row.N_input} | '
            f'{row.N_background_sources} | {100*row.f_pass_known:.2f}% '
            f'({100*row.f_pass_lower_95:.2f}–{100*row.f_pass_upper_95:.2f}%) | '
            f'{100*row.relative_halfwidth_95:.1f}% | {row.known_weight_fraction:.1%} |')
    lines += ['', '## Missing-control sensitivity ranges', '',
        'The ranges below allow unknown outcomes to take any outcome permitted by the '
        'empirical pair-only pass factor. They hold the measured outcomes fixed and are '
        'not 95% statistical intervals. Missing measurements are never treated as rejections.', '',
        '| Route | Pass-rate range allowing unknown outcomes | Matching coverage |',
        '| --- | ---: | ---: |']
    for row in main.itertuples():
        lines.append(f'| {row.route} | {row.f_pass_population_lower:.2%}–'
                     f'{row.f_pass_population_upper:.2%} | {row.matching_coverage:.1%} |')
    failed = main.loc[~main.precision_met.eq(True)]
    lines += ['', f'{int(main.precision_met.eq(True).sum())} of {len(main)} routes meet '
        'the strict 20% relative-precision target for the main 3–10 arcsec annulus. '
        + ('Remaining: '+', '.join(f'{row.route}: {100*row.relative_halfwidth_95:.1f}%'
            for row in failed.itertuples())+'.' if len(failed) else 'All pass.'),
        'Relative uncertainty is the larger distance from the estimate to either 95% '
        'endpoint, divided by the estimate. Individual latitude-bin intervals are wider.', '',
        '## Assumptions and interpretation', '',
        '- Wide sources represent unrelated neighbours; some physical wide pairs may remain.',
        '- Sources are matched in brightness, native survey region, Galactic hemisphere, '
        'latitude, sky region and foreground reddening. This is a transfer assumption.',
        '- Component controls keep the actual pre-vetting input primary and its redshift/error. '
        'Blended controls use isolated, independently measured quasar donors, matched in '
        'brightness, redshift and catalogue family, retaining their own redshift/errors.',
        '- Component backgrounds retain their coordinates for dust and Galactic scoring; '
        'blends use donor coordinates. The sky matching is not exact flux relocation.',
        '- Native component fluxes and variances are added separately. Missing fluxes are '
        'not zero-filled. The production Legacy/fallback survey policy is retained; common '
        'measured band masks can still differ from those of the original inputs.',
        '- The production either-model rules and calibrated threshold are reused. Imaging '
        'requires both components to pass under the same model variant. Secondary PM must '
        'be below 3 when finite; missing PM is retained and flagged. Imaging primary PM '
        'must be finite and below 3. DCR applies only to one-Gaia imaging.',
        '- Empirical DCR and morphology pass factors are taken from the real cut flow. '
        'Their transfer uncertainty is not included in the statistical interval.',
        '- Reused background sources and donor quasars keep shared field identities in '
        'the uncertainty calculation. Mock combinations are not independent observations.',
        '- The 10–30 arcsec annulus is a separate check and is not pooled with 3–10 arcsec.',
        '- These rates depend on the tested input population and photometric coverage; '
        'they are not universal constants for every survey or magnitude range.', '',
        '## What this measures', '',
        'A false-positive acceptance probability is P(pass | unrelated). A contamination '
        'fraction is P(unrelated | pass), which additionally depends on how common unrelated '
        'and associated pairs are before selection and on acceptance of associated pairs. '
        'This study reports the former. It does not estimate a selected-sample contaminant '
        'count or purity. Lenses and physical duals are not distinguished by this test.', '',
        'The 10–60 arcsec density products are separate neighbourhood diagnostics under '
        'C=1 geometric scaling. They are not combined with these pass rates to infer a '
        'candidate count. The stellar close/outer ratio is confounded by physical stellar '
        'binaries and is not used as a measured detection probability.', '']
    density_report = root/'FIELD_DENSITY.md'
    if density_report.is_file():
        lines += density_report.read_text().splitlines() + ['']
    else:
        lines += ['## Matched field-density check', '',
                  'Not yet evaluated. Run the report stage with complete close/outer Gaia '
                  'cones before drawing conclusions about close-neighbour field density.', '']
    lines += ['## Products', '',
        '- `repaired_summary_by_route.csv`: pass rates, intervals, coverage and missing-outcome bounds.',
        '- `repaired_summary_by_latitude.csv`: the same quantities split by Galactic latitude.',
        '- `repaired_control_pass_rates.pdf` / `.png`: main and check annuli by latitude.',
        '- `control_magnitude_matching.csv` / `.pdf` / `.png`: magnitude-matching diagnostics.',
        '- `parent_densities_by_latitude.csv` / `local_density_scenarios.pdf`: neighbourhood diagnostics.',
        '- `matched_density_comparison.csv`: same-field close counts and geometric predictions.',
        '- `matched_density_field_counts.csv`: per-field counts, exact IDs and inclusion weights.',
        '- `observed_vs_field_density.pdf` / `.png` and `FIELD_DENSITY.md`: density comparison.',
        '- `field_density_inputs/`: raw close/outer Gaia data and coverage for offline replay.',
        '- Individual control score Parquets retain source identities, scores and flags.', '']
    (root/'RESULTS.md').write_text('\n'.join(lines))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    write_report(parser.parse_args().output)
