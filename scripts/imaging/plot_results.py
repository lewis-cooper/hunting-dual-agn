import json
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18
from dual_agn.settings import *
from dual_agn.photometry import configure_dustmaps
from qso_pcolor import Photometry
from qso_pcolor.extinction import deredden
from qso_pcolor.unified import UnifiedPSFModel
import dual_agn.colour_projection as colour_tools
from dual_agn.plot_style import plt, save_figure
configure_dustmaps()
log_r_per_unit_z_threshold = json.loads(THRESHOLD_CALIBRATION_PATH.read_text())['adopted_log_r_threshold']

MODEL_POINTERS = MODEL_BUNDLES
IMAGING_PLOT_MODEL = 'either'
paper_figure_root = FIGURE_ROOT / 'imaging'
models = {name: UnifiedPSFModel.load(path) for name, path in MODEL_BUNDLES.items()}
system_scores = pd.read_csv(
    IMAGING_SCORES_PATH,
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)
component_photometry = pd.read_csv(
    IMAGING_ROOT / 'component_photometry.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)

def save_paper_figure(figure, stem):
    save_figure(figure, stem, 'imaging')
fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)
for axis, model in zip(axes[:2], MODEL_POINTERS):
    eligible = system_scores[f'both_components_eligible_{model}']
    axis.scatter(
        system_scores.loc[eligible, f'primary_log_r_per_unit_z_{model}'],
        system_scores.loc[eligible, f'secondary_log_r_per_unit_z_{model}'],
        s=18,
        alpha=0.6
    )
    axis.set(
        xlabel='Primary log R per unit z',
        ylabel='Secondary log R per unit z',
        title=f'{model}: {eligible.sum()} eligible systems'
    )
both_models_eligible = system_scores.both_components_eligible_magdep & system_scores.both_components_eligible_magindep
axes[2].scatter(
    system_scores.loc[both_models_eligible, 'minimum_component_log_r_magdep'],
    system_scores.loc[both_models_eligible, 'minimum_component_log_r_magindep'],
    s=18,
    alpha=0.6
)
axes[2].set(
    xlabel='Minimum component log R: magdep',
    ylabel='Minimum component log R: magindep',
    title=f'Both variants eligible: {both_models_eligible.sum()}'
)
save_paper_figure(fig, 'model_comparison')
plt.close("all")


imaging_score_selected = system_scores.loc[system_scores.joint_score_pass.eq(True)].copy()
imaging_plot_label = f'Both component log R ≥ {log_r_per_unit_z_threshold:g} in the same model'

print(imaging_plot_label, '| score-selected systems:', len(imaging_score_selected))
print('Selection stage: after the upstream one-Gaia-only DCR cut; before secondary PM and morphology cuts.')
print('Of these, passing the secondary PM rule:', int(imaging_score_selected.secondary_pm_pass.sum()))
print(
    'Saved morphology:',
    int(imaging_score_selected.morphology_available.sum()),
    '| missing:',
    int((~imaging_score_selected.morphology_available).sum())
)
print('Plot outputs:', paper_figure_root)
imaging_colour_model = models['magdep'].base.model  # Reference locus overlaid on the either-model sample.
colour_pairs = [('g', 'r'), ('r', 'z')]
imaging_colour_components = component_photometry.loc[component_photometry.qso_parent_index.isin(imaging_score_selected.qso_parent_index)].copy()
for first_band, second_band in colour_pairs:
    imaging_colour_components[f'{first_band}_minus_{second_band}'] = np.nan
for native_system, group in imaging_colour_components.groupby('photometric_system'):
    bands = tuple((f'{native_system}:{band}' for band in 'grz'))
    photometry = Photometry(
        group[[f'flux_{band}' for band in 'grz']].to_numpy(float),
        group[[f'variance_{band}' for band in 'grz']].to_numpy(float),
        bands
    )
    coordinates = SkyCoord(
        group.component_ra.to_numpy(float) * u.deg,
        group.component_dec.to_numpy(float) * u.deg
    )
    corrected_photometry, _, _ = deredden(
        photometry,
        coordinates.galactic.l.deg,
        coordinates.galactic.b.deg,
        imaging_colour_model.meta['extinction']
    )
    features = imaging_colour_model.transform(corrected_photometry)
    for first_band, second_band in colour_pairs:
        first_index = features.labels.index(f'{native_system}:{first_band}')
        second_index = features.labels.index(f'{native_system}:{second_band}')
        usable_colour = features.observed[:, first_index] & features.observed[:, second_index] & np.isfinite(features.x[:, first_index]) & np.isfinite(features.x[:, second_index])
        imaging_colour_components.loc[group.index[usable_colour], f'{first_band}_minus_{second_band}'] = features.x[usable_colour, first_index] - features.x[usable_colour, second_index]
imaging_colour_pairs = imaging_score_selected.copy()
for role in ('primary', 'secondary'):
    component_colours = imaging_colour_components.loc[
        imaging_colour_components.component_role.eq(role),
        ['qso_parent_index', 'g_minus_r', 'r_minus_z']
    ].rename(columns={name: f'{role}_luptitude_{name}' for name in ('g_minus_r', 'r_minus_z')})
    imaging_colour_pairs = imaging_colour_pairs.merge(component_colours, on='qso_parent_index', how='left', validate='one_to_one')
print('Colour coordinates: extinction-corrected native Legacy luptitude differences; weak/negative fluxes retained.')
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
imaging_colour_log_r_minimum = log_r_per_unit_z_threshold
minimum_score_column = f'minimum_component_log_r_{IMAGING_PLOT_MODEL}'
imaging_colour_plot_rows = imaging_colour_pairs.loc[imaging_colour_pairs[minimum_score_column].ge(imaging_colour_log_r_minimum)]
score_values = imaging_colour_plot_rows[minimum_score_column]
score_scale = Normalize(score_values.min(), score_values.max())
native_weights = imaging_colour_components.loc[imaging_colour_components.qso_parent_index.isin(imaging_colour_plot_rows.qso_parent_index)].photometric_system.value_counts(normalize=True).sort_index()
model_redshift = np.unique(np.r_[
    np.linspace(*imaging_colour_model.qso.support, 200),
    imaging_colour_model.qso.z_centres
])
redshift_min = imaging_colour_plot_rows.adopted_redshift.min()
redshift_max = imaging_colour_plot_rows.adopted_redshift.max()
# Leave space for markers centred on the lowest and highest redshifts.
redshift_padding = 0.02 * (redshift_max - redshift_min)
redshift_limits = (redshift_min - redshift_padding, redshift_max + redshift_padding)
figure, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
imaging_colour_coverage = []
for axis, (first_band, second_band) in zip(axes, colour_pairs):
    colour_name = f'{first_band}_minus_{second_band}'
    primary_colour, secondary_colour = [f'{role}_luptitude_{colour_name}' for role in ('primary', 'secondary')]
    finite_redshift = np.isfinite(imaging_colour_plot_rows.adopted_redshift)
    finite_colours = np.isfinite(imaging_colour_plot_rows[[primary_colour, secondary_colour]]).all(axis=1)
    panel = imaging_colour_plot_rows.loc[finite_redshift & finite_colours]
    projected_models = {native_system: colour_tools.colour_model(
        imaging_colour_model.qso,
        (f'{native_system}:{first_band}', f'{native_system}:{second_band}')
    ) for native_system in ('decals_dr9_south', 'decals_dr9_north')}
    median_colour = np.array([colour_tools.colour_quantile(colour_tools.mixture_at_z(projected_models['decals_dr9_north'], redshift), 0.5) for redshift in model_redshift])
    all_colours = np.r_[panel[[primary_colour, secondary_colour]].to_numpy().ravel(), median_colour]
    padding = max(0.06 * np.ptp(all_colours), 0.1)
    colour_limits = (all_colours.min() - padding, all_colours.max() + padding)
    colour_grid = np.linspace(*colour_limits, 250)
    model_density = sum((
        weight * np.exp(projected_models[native_system].log_p_colour_given_z(colour_grid[:, None], None, model_redshift)) for native_system,
        weight in native_weights.items()
    ))
    axis.pcolormesh(
        colour_grid,
        model_redshift,
        model_density.T,
        norm=LogNorm(model_density.max() * 0.001, model_density.max()),
        cmap='Greys',
        alpha=0.22,
        shading='auto',
        rasterized=True
    )
    axis.plot(median_colour, model_redshift, color='firebrick', linewidth=1.8, linestyle='--', zorder=4)
    for row in panel.itertuples():
        axis.plot(
            [getattr(row, primary_colour), getattr(row, secondary_colour)],
            [row.adopted_redshift] * 2,
            color='black',
            linewidth=0.6,
            alpha=0.4
        )
    axis.scatter(
        panel[primary_colour],
        panel.adopted_redshift,
        c=panel[minimum_score_column],
        cmap='viridis',
        norm=score_scale,
        marker='o',
        s=36,
        alpha=0.75,
        edgecolors='black',
        linewidths=0.35,
        zorder=3
    )
    axis.scatter(
        panel[secondary_colour],
        panel.adopted_redshift,
        c=panel[minimum_score_column],
        cmap='viridis',
        norm=score_scale,
        marker='^',
        s=36,
        alpha=0.75,
        edgecolors='black',
        linewidths=0.35,
        zorder=3
    )
    axis.set(xlim=colour_limits, ylim=redshift_limits)
    axis.set_box_aspect(1)
    axis.set_xlabel(f'${first_band}-{second_band}$ (mag)', fontsize=18)
    axis.tick_params(axis='both', labelsize=15)
    axis.grid(alpha=0.15)
    imaging_colour_coverage.append({
        'colour': colour_name,
        'systems_plotted': len(panel),
        'missing_colour_or_redshift': len(imaging_colour_plot_rows) - len(panel)
    })
axes[0].set_ylabel('Primary redshift, $z_0$', fontsize=18)
colourbar = figure.colorbar(
    ScalarMappable(norm=score_scale, cmap='viridis'),
    cax=axes[-1].inset_axes([1.04, 0, 0.05, 1])
)
colourbar.set_label('Pair score, $\\ln\\mathcal{R}$', fontsize=18)
colourbar.ax.tick_params(labelsize=15)
figure.legend(
    handles=[
        Line2D([], [], color='firebrick', linestyle='--', linewidth=1.8, label='Model'),
        Line2D(
            [],
            [],
            marker='o',
            linestyle='none',
            markerfacecolor='0.5',
            markeredgecolor='black',
            markeredgewidth=0.5,
            markersize=6,
            label='Primary'
        ),
        Line2D(
            [],
            [],
            marker='^',
            linestyle='none',
            markerfacecolor='0.5',
            markeredgecolor='black',
            markeredgewidth=0.5,
            markersize=6,
            label='Secondary'
        )
    ],
    loc='lower center',
    bbox_to_anchor=(0.49, 0.02),
    ncol=3,
    fontsize=15,
    handletextpad=0.6,
    columnspacing=1.8
)
figure.subplots_adjust(left=0.1, right=0.88, bottom=0.22, top=0.96, wspace=0.1)
save_paper_figure(figure, 'colour_redshift')
plt.close("all")
print(f'Both-component colour-plot cut: log_r_per_unit_z >= {imaging_colour_log_r_minimum:.2f}; {len(imaging_colour_plot_rows)}/{len(imaging_colour_pairs)} systems retained.')
for coverage in imaging_colour_coverage:
    print(f"{coverage['colour']}: {coverage['systems_plotted']}/{len(imaging_colour_plot_rows)} systems plotted; "
          f"{coverage['missing_colour_or_redshift']} missing colours or redshift.")


pre_morphology_pass = system_scores.joint_pre_morphology_pass.fillna(False)
finite_projection = np.isfinite(system_scores.adopted_redshift) & system_scores.adopted_redshift.gt(0) & np.isfinite(system_scores.component_sep_arcsec) & system_scores.component_sep_arcsec.gt(0)
imaging_separation_rows = system_scores.loc[
    pre_morphology_pass & finite_projection,
    [
        'qso_parent_index',
        'candidate_tier',
        'adopted_redshift',
        'component_sep_arcsec',
        'r_perp_sigma_ratio'
    ]
].copy()
imaging_separation_rows['projected_separation_kpc'] = imaging_separation_rows.component_sep_arcsec.to_numpy() * Planck18.kpc_proper_per_arcmin(imaging_separation_rows.adopted_redshift.to_numpy()).value / 60
separation_bin_width_kpc = 0.5
separation_maximum_kpc = max(
    14.0,
    separation_bin_width_kpc * np.ceil(imaging_separation_rows.projected_separation_kpc.max() / separation_bin_width_kpc)
)
separation_bins_kpc = np.arange(0, separation_maximum_kpc + separation_bin_width_kpc, separation_bin_width_kpc)
figure, axis = plt.subplots(figsize=(8, 5), layout='constrained')
for tier, label, colour, linestyle in [('one_gaia', 'One Gaia', 'tab:blue', '--'), ('two_gaia', 'Two Gaia', 'tab:red', '-')]:
    group = imaging_separation_rows.loc[imaging_separation_rows.candidate_tier.eq(tier)]
    axis.hist(
        group.projected_separation_kpc,
        bins=separation_bins_kpc,
        histtype='step',
        color=colour,
        linestyle=linestyle,
        linewidth=2.2,
        label=f'{label} ($N={len(group)}$)'
    )
axis.set(
    xlabel='Projected separation, $r_{\\mathrm{p}}$ [kpc]',
    ylabel=f'$N$ systems per {separation_bin_width_kpc:g} kpc bin',
    xlim=(0, separation_maximum_kpc),
    ylim=(0, None)
)
axis.minorticks_on()
axis.tick_params(axis='both', which='both', top=True, right=True)
axis.grid(axis='x', alpha=0.18, linewidth=0.6)
axis.grid(axis='y', alpha=0.1, linewidth=0.5)
axis.legend(loc='upper left', frameon=False)
save_paper_figure(figure, 'pre_morphology_separation')
plt.close("all")
print(f'Imaging after score + secondary PM, before morphology: {pre_morphology_pass.sum()} systems; {len(imaging_separation_rows)} plotted; {(pre_morphology_pass & ~finite_projection).sum()} lack a valid projection. DCR was applied upstream only to one-Gaia pairs. Bin width: {separation_bin_width_kpc:g} kpc; cosmology: Planck18.')
print(f'Morphology is not applied in this plot: {imaging_separation_rows.r_perp_sigma_ratio.lt(R_PERP_MIN).sum()} plotted systems fall below R_perp = {R_PERP_MIN:g}.')
print('Paper figure folder:', paper_figure_root)
