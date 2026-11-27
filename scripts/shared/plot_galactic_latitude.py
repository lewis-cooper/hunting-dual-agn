"""Plot route and tier distributions per geometric Galactic-latitude area."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from dual_agn.candidate_tiers import assign_candidate_tiers, boolean_flag
from dual_agn.plot_style import plt
from dual_agn.settings import FIGURE_ROOT, IDENTIFIER_DTYPES, FINAL_CANDIDATES_PATH, WORK_ROOT


def plot_latitude(candidates, figure_root, summary_dir, axis_break=True):
    candidates = assign_candidate_tiers(candidates)
    if candidates.candidate_id.isna().any() or not candidates.candidate_id.is_unique:
        raise ValueError('Candidate identities must be unique and nonmissing')
    candidate_sky = SkyCoord(
        pd.to_numeric(candidates.primary_ra_deg, errors='raise').to_numpy() * u.deg,
        pd.to_numeric(candidates.primary_dec_deg, errors='raise').to_numpy() * u.deg)
    abs_latitude_deg = np.abs(candidate_sky.galactic.b.deg)
    if not np.isfinite(abs_latitude_deg).all():
        raise ValueError('Nonfinite candidate latitude; no rows were discarded')
    latitude_bin_edges_deg = np.arange(0, 95, 5)
    latitude_bin_area_deg2 = (4 * np.pi * (180 / np.pi) ** 2
                             * np.diff(np.sin(np.deg2rad(latitude_bin_edges_deg))))
    latitude_area_summary = pd.DataFrame({
        'abs_b_min_deg': latitude_bin_edges_deg[:-1],
        'abs_b_max_deg': latitude_bin_edges_deg[1:], 'geometric_area_deg2': latitude_bin_area_deg2})
    imaging = boolean_flag(candidates, 'selected_by_imaging')
    groups = [
        ('resolved', boolean_flag(candidates, 'selected_by_resolved'),
         'Main: Gaia neighbour (resolved)', 'tab:blue', '-'),
        ('blended', boolean_flag(candidates, 'selected_by_unresolved'),
         'Main: Gaia neighbour (unresolved)', 'tab:green', '--'),
        ('imaging_two_gaia', imaging & candidates.imaging_selection_route.eq('two_gaia'),
         'Main: Imaging (two-Gaia)', 'tab:red', '-.'),
        ('imaging_one_gaia', candidates.tier.eq('additional'),
         'Additional: Imaging (one-Gaia)', 'tab:red', ':'),
    ]
    if axis_break:
        break_start_deg, break_end_deg = 1., 24.
        hidden = (abs_latitude_deg > break_start_deg) & (abs_latitude_deg < break_end_deg)
        if hidden.any():
            raise ValueError('The requested axis break would hide candidates; use --no-axis-break')
        figure, axes = plt.subplots(1, 2, figsize=(8, 5), sharey=True,
                                   width_ratios=[1, 66])
        # Equal visible spans beside the gap centre the break between the 0 and 25 ticks.
        figure.subplots_adjust(left=.12, right=.98, bottom=.15, top=.97, wspace=.025)
        axes[0].set_xlim(0, break_start_deg)
        axes[0].set_xticks([0])
        axes[1].set_xlim(break_end_deg, 90)
        axes[1].set_xticks([25, 30, 40, 50, 60, 70, 80, 90])
        axes[0].spines['right'].set_visible(False)
        axes[1].spines['left'].set_visible(False)
        for axis, edge in zip(axes, (1, 0)):
            axis.plot([edge, edge], [0, 1], transform=axis.transAxes,
                      marker=[(-1, -.5), (1, .5)], markersize=10,
                      linestyle='none', color='black', clip_on=False)
    else:
        figure, axis = plt.subplots(figsize=(8, 5), layout='constrained')
        axes = [axis]
        axis.set_xlim(0, 90)
    for name, selected, label, colour, linestyle in groups:
        counts, _ = np.histogram(abs_latitude_deg[selected], bins=latitude_bin_edges_deg)
        density_per_1000_deg2 = counts / latitude_bin_area_deg2 * 1000
        latitude_area_summary[f'{name}_count'] = counts
        latitude_area_summary[f'{name}_per_1000_deg2'] = density_per_1000_deg2
        for axis in axes:
            axis.stairs(density_per_1000_deg2, latitude_bin_edges_deg,
                        color=colour, linestyle=linestyle, linewidth=2.2, label=label)
        print(label + ':', int(selected.sum()))
    for tier in ('main', 'additional'):
        counts, _ = np.histogram(abs_latitude_deg[candidates.tier.eq(tier)], bins=latitude_bin_edges_deg)
        latitude_area_summary[f'{tier}_unique_count'] = counts
    for axis in axes:
        axis.axvline(25, color='0.45', linestyle=':', linewidth=1.5,
                     label=r'$|b|=25^\circ$ selection limit')
        axis.set_ylim(0, None)
        axis.minorticks_on()
        axis.tick_params(axis='both', which='both', top=True, right=True)
        axis.grid(alpha=0.15)
    if axis_break:
        axes[0].tick_params(axis='y', which='both', right=False)
        axes[1].tick_params(axis='y', which='both', left=False, labelleft=False)
    axes[0].set_ylabel(r'$N$ systems per $1000\,\mathrm{deg}^2$')
    figure.supxlabel(r'Absolute Galactic latitude, $|b|$ [deg]', fontsize=16)
    axes[-1].legend(loc='upper right', frameon=False, fontsize=9.5,
                    labelspacing=.25, borderaxespad=.35)
    figure_root.mkdir(parents=True, exist_ok=True)
    for extension in ('png', 'pdf'):
        figure_path = figure_root / f'galactic_latitude.{extension}'
        figure.savefig(figure_path, dpi=220, bbox_inches='tight')
        print('Saved:', figure_path.resolve())
    plt.close(figure)
    summary_dir.mkdir(parents=True, exist_ok=True)
    latitude_area_summary.to_csv(summary_dir / 'galactic_latitude.csv', index=False)
    shared = int(boolean_flag(candidates, 'selected_by_both_searches').sum())
    print(f'{len(candidates)} unique systems; {shared} shared systems appear in both applicable route curves.')
    print('Five-degree bins use full geometric sky area, without footprint or parent-density corrections.')
    if axis_break:
        print('Display-only x-axis break: 1--24 degrees; no candidates or histogram bins were removed.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidates', type=Path, default=FINAL_CANDIDATES_PATH)
    parser.add_argument('--figure-dir', type=Path, default=FIGURE_ROOT / 'joint')
    parser.add_argument('--summary-dir', type=Path, default=WORK_ROOT)
    parser.add_argument('--no-axis-break', action='store_true', help='Display the full 0--90 degree range')
    args = parser.parse_args()
    candidates = pd.read_csv(args.candidates, dtype=IDENTIFIER_DTYPES, float_precision='round_trip')
    print('Input:', args.candidates.resolve())
    plot_latitude(candidates, args.figure_dir, args.summary_dir, axis_break=not args.no_axis_break)


if __name__ == '__main__':
    main()
