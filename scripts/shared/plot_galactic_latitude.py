"""Plot Galactic latitude for the joint vetted sample."""

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from dual_agn.plot_style import plt
from dual_agn.settings import FIGURE_ROOT, IDENTIFIER_DTYPES, FINAL_CANDIDATES_PATH, WORK_ROOT

figure_root = FIGURE_ROOT / 'joint'
figure_root.mkdir(parents=True, exist_ok=True)
WORK_ROOT.mkdir(parents=True, exist_ok=True)
joint_candidates = pd.read_csv(
    FINAL_CANDIDATES_PATH,
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip',
)
candidate_sky = SkyCoord(
    joint_candidates.primary_ra_deg.to_numpy() * u.deg,
    joint_candidates.primary_dec_deg.to_numpy() * u.deg
)
abs_latitude_deg = np.abs(candidate_sky.galactic.b.deg)
finite_latitude = np.isfinite(abs_latitude_deg)
latitude_bin_edges_deg = np.arange(0, 95, 5)
latitude_bin_area_deg2 = 4 * np.pi * (180 / np.pi) ** 2 * np.diff(np.sin(np.deg2rad(latitude_bin_edges_deg)))
latitude_area_summary = pd.DataFrame({
    'abs_b_min_deg': latitude_bin_edges_deg[:-1],
    'abs_b_max_deg': latitude_bin_edges_deg[1:],
    'geometric_area_deg2': latitude_bin_area_deg2
})
figure, axis = plt.subplots(figsize=(8, 5), layout='constrained')
for name, selection, label, colour, linestyle in [
    (
        'resolved',
        joint_candidates.selected_by_resolved,
        'One-arcsecond resolved',
        'tab:blue',
        '-'
    ),
    (
        'unresolved',
        joint_candidates.selected_by_unresolved,
        'One-arcsecond unresolved',
        'tab:green',
        '--'
    ),
    ('imaging', joint_candidates.selected_by_imaging, 'Imaging', 'tab:red', '-.')
]:
    selected = finite_latitude & selection.to_numpy()
    counts, _ = np.histogram(abs_latitude_deg[selected], bins=latitude_bin_edges_deg)
    density_per_1000_deg2 = counts / latitude_bin_area_deg2 * 1000
    latitude_area_summary[f'{name}_count'] = counts
    latitude_area_summary[f'{name}_per_1000_deg2'] = density_per_1000_deg2
    axis.stairs(
        density_per_1000_deg2,
        latitude_bin_edges_deg,
        color=colour,
        linestyle=linestyle,
        linewidth=2.2,
        label=f'{label} ($N={selected.sum()}$)'
    )
axis.axvline(25, color='0.45', linestyle=':', linewidth=1.5, label='$|b|=25^\\circ$ selection limit')
axis.set(
    xlabel='Absolute Galactic latitude, $|b|$ [deg]',
    ylabel='$N$ systems per $1000\\,\\mathrm{deg}^2$',
    xlim=(0, 90),
    ylim=(0, None)
)
axis.minorticks_on()
axis.tick_params(axis='both', which='both', top=True, right=True)
axis.grid(alpha=0.15)
axis.legend(loc='upper left', frameon=False, fontsize=11)
for extension in ('png', 'pdf'):
    figure_path = figure_root / f'galactic_latitude.{extension}'
    figure.savefig(figure_path, dpi=180, bbox_inches='tight')
    print('Saved:', figure_path)
plt.close("all")
latitude_area_summary.to_csv(WORK_ROOT / 'galactic_latitude.csv', index=False)
print(f'{len(joint_candidates)} unique final systems; {(~finite_latitude).sum()} missing latitudes. {joint_candidates.selected_by_both_searches.sum()} shared systems appear in both method curves.')
print('Counts use 5-degree bins, divided by full geometric latitude-band area; the effective survey footprint and parent-quasar density are not corrected.')
