"""Plot the unique final candidates in projected separation versus redshift."""

import argparse
from pathlib import Path

from astropy.cosmology import Planck18
import numpy as np
import pandas as pd

from dual_agn.plot_style import plt
from dual_agn.settings import FIGURE_ROOT, FINAL_CANDIDATES_PATH


CANDIDATE_STYLES = {
    'Imaging: one Gaia': ('Imaging only: one Gaia', '#00BFD8', 'o', 48),
    'Imaging: two Gaia': ('Imaging only: two Gaia', '#E832A2', 'D', 52),
    'One-arcsecond resolved': ('One-arcsecond only: resolved', 'tab:red', '^', 42),
    'One-arcsecond unresolved': ('One-arcsecond only: unresolved', 'tab:orange', 'v', 42),
    'Both searches': ('Recovered by both searches', 'black', 'P', 52),
}
PAPER_STYLE = {
    'font.serif': ['Computer Modern Roman', 'CMU Serif', 'DejaVu Serif'],
    'axes.unicode_minus': False,
    'axes.labelsize': 18,
    'xtick.top': True,
    'ytick.right': True,
    'xtick.major.width': 1.6,
    'ytick.major.width': 1.6,
    'xtick.minor.width': 1.2,
    'ytick.minor.width': 1.2,
    'xtick.major.size': 7,
    'ytick.major.size': 7,
    'xtick.minor.size': 4,
    'ytick.minor.size': 4,
}


def read_candidates(candidates_path):
    candidates = pd.read_csv(
        candidates_path,
        usecols=['candidate_id', 'plot_group', 'redshift', 'projected_separation_kpc'],
        dtype={'candidate_id': 'string', 'plot_group': 'string'},
        float_precision='round_trip',
    )
    if candidates.candidate_id.isna().any() or not candidates.candidate_id.is_unique:
        raise ValueError('The final catalogue must have a unique, non-missing candidate_id per system.')
    unknown_groups = ~candidates.plot_group.isin(CANDIDATE_STYLES)
    if unknown_groups.any():
        raise ValueError(f'Unrecognised plot_group values: {candidates.loc[unknown_groups, "plot_group"].unique().tolist()}')
    finite_coordinates = (
        np.isfinite(candidates.redshift)
        & np.isfinite(candidates.projected_separation_kpc)
    )
    positive_coordinates = (
        candidates.redshift.gt(0)
        & candidates.projected_separation_kpc.gt(0)
    )
    plot_mask = finite_coordinates & positive_coordinates
    print(f'Input: {candidates_path.resolve()}')
    print(f'Unique systems: {len(candidates)}; plotted: {plot_mask.sum()}; '
          f'excluded nonfinite or nonpositive coordinates: {(~plot_mask).sum()}')
    if not plot_mask.any():
        raise ValueError('No candidates have finite, positive redshift and projected separation.')
    return candidates.loc[plot_mask].copy()


def plot_candidates(candidates, figure_directory, dpi):
    # Projected separations in the combined catalogue already use Planck18.
    redshift_max = max(4.0, 1.05 * candidates.redshift.max())
    projected_max_kpc = max(14.0, 1.05 * candidates.projected_separation_kpc.max())
    with plt.rc_context(PAPER_STYLE):
        figure, axis = plt.subplots(figsize=(12.4, 7.4), layout='constrained')
        axis.set_box_aspect(1 / 1.65)
        for group, (label, colour, marker, size) in CANDIDATE_STYLES.items():
            selected = candidates.loc[candidates.plot_group.eq(group)]
            print(f'{label}: {len(selected)}')
            if selected.empty:
                continue
            axis.scatter(
                selected.redshift, selected.projected_separation_kpc,
                s=size, color=colour, marker=marker, edgecolors='black',
                linewidths=0.55, alpha=0.82, zorder=3, label=label,
            )

        redshift_grid = np.linspace(0.001, redshift_max, 600)
        kpc_per_arcsec = Planck18.kpc_proper_per_arcmin(redshift_grid).value / 60
        for separation_arcsec, colour in ((0.3, '0.70'), (0.5, '0.58'),
                                          (1.0, '0.43'), (1.5, '0.28')):
            projected_kpc = separation_arcsec * kpc_per_arcsec
            axis.plot(redshift_grid, projected_kpc, color=colour,
                      linestyle=':', linewidth=1.6, zorder=0)
            label_index = int(0.87 * len(redshift_grid))
            axis.text(
                redshift_grid[label_index], projected_kpc[label_index],
                f'{separation_arcsec:g}' + r'$^{\prime\prime}$', color=colour,
                fontsize=12, ha='left', va='bottom', rotation=12,
            )

        axis.set(xlim=(0, redshift_max), ylim=(0, projected_max_kpc),
                 xlabel=r'Redshift, $z$', ylabel='Projected separation [kpc]')
        axis.minorticks_on()
        axis.grid(False)
        axis.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
                    frameon=False, borderaxespad=0, handletextpad=0.45,
                    labelspacing=0.38, fontsize=11)
        figure_directory.mkdir(parents=True, exist_ok=True)
        for extension in ('png', 'pdf'):
            figure_path = figure_directory / f'separation_redshift.{extension}'
            figure.savefig(figure_path, dpi=dpi, bbox_inches='tight')
            print(f'Saved: {figure_path.resolve()}')
        plt.close(figure)
    print('Coordinates: adopted catalogue redshift and projected proper separation (Planck18).')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidates', type=Path, default=FINAL_CANDIDATES_PATH,
                        help='Unique final candidate CSV from combine_candidate_catalogues')
    parser.add_argument('--figure-dir', type=Path, default=FIGURE_ROOT / 'joint')
    parser.add_argument('--dpi', type=int, default=250)
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error('--dpi must be positive')
    candidates = read_candidates(args.candidates.expanduser())
    plot_candidates(candidates, args.figure_dir.expanduser(), args.dpi)


if __name__ == '__main__':
    main()
