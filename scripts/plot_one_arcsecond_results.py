import json
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from settings import *
from photometry import configure_dustmaps
from qso_pcolor import Photometry
from qso_pcolor.extinction import deredden
from qso_pcolor.unified import UnifiedPSFModel
import colour_projection as colour_tools
from plot_style import plt, save_figure
from score_selection import eligible_model_scores
configure_dustmaps()
calibration = json.loads(THRESHOLD_CALIBRATION_PATH.read_text())
log_r_per_unit_z_threshold = calibration['adopted_log_r_threshold']

ranked = pd.read_csv(ONE_ROOT / 'ranked_resolved.csv', dtype=IDENTIFIER_DTYPES, float_precision='round_trip')
resolved_models = {name: UnifiedPSFModel.load(path) for name, path in MODEL_BUNDLES.items()}
model_bundles = {name: {'bundle': path} for name, path in MODEL_BUNDLES.items()}

def save_paper_figure(figure, stem):
    save_figure(figure, stem, 'one_arcsecond')
ranked['log_r_per_unit_z'] = eligible_model_scores(ranked, 'eligible', 'log_r_per_unit_z').max(axis=1)
finite_score = np.isfinite(ranked.log_r_per_unit_z)
finite_redshift = np.isfinite(ranked['primary_qso_z'])
plot_rows = ranked.loc[finite_score & finite_redshift].copy()
model = resolved_models['magdep'].base.model
resolved_colour_panels = []
if plot_rows.empty:
    print('No eligible current-model resolved ranks; colour/rank plots are omitted.')
else:
    colour_pairs = [('g', 'r'), ('r', 'z')]
    for system_name, group in plot_rows.groupby('photometric_system'):
        bands = tuple((f'{system_name}:{band}' for band in 'grz'))
        for component in ('primary', 'secondary'):
            flux = group[[f'{component}_decals_flux_{band}' for band in 'grz']].to_numpy(float)
            inverse_variance = group[[f'{component}_decals_flux_ivar_{band}' for band in 'grz']].to_numpy(float)
            positive_ivar = np.isfinite(inverse_variance) & (inverse_variance > 0)
            variance = np.divide(
                1.0,
                inverse_variance,
                out=np.full_like(inverse_variance, np.inf),
                where=positive_ivar
            )
            component_sky = SkyCoord(
                group[f'{component}_component_ra'].to_numpy(float) * u.deg,
                group[f'{component}_component_dec'].to_numpy(float) * u.deg
            )
            corrected_photometry, _, _ = deredden(
                Photometry(flux, variance, bands),
                component_sky.galactic.l.deg,
                component_sky.galactic.b.deg,
                model.meta['extinction']
            )
            features = model.transform(corrected_photometry)
            for first_band, second_band in colour_pairs:
                first_index = features.labels.index(f'{system_name}:{first_band}')
                second_index = features.labels.index(f'{system_name}:{second_band}')
                plot_rows.loc[group.index, f'{component}_{first_band}_minus_{second_band}'] = features.x[:, first_index] - features.x[:, second_index]
    redshift_limits = (0, 3.5)
    redshift_grid_points = 300
    colour_grid_points = 300
    model_redshift_min = max(redshift_limits[0], model.qso.support[0])
    model_redshift_max = min(redshift_limits[1], model.qso.support[1])
    model_redshift = np.unique(np.r_[
        np.linspace(model_redshift_min, model_redshift_max, redshift_grid_points),
        model.qso.z_centres[(model.qso.z_centres >= model_redshift_min) & (model.qso.z_centres <= model_redshift_max)]
    ])
    region_weights = plot_rows['photometric_system'].value_counts(normalize=True).sort_index()
    for first_band, second_band in colour_pairs:
        colour_name = f'{first_band}_minus_{second_band}'
        component_columns = [f'{component}_{colour_name}' for component in ('primary', 'secondary')]
        finite_colours = np.isfinite(plot_rows[component_columns]).all(axis=1)
        within_redshift_limits = plot_rows['primary_qso_z'].between(*redshift_limits)
        panel = plot_rows.loc[finite_colours & within_redshift_limits]
        projected_models = {}
        model_medians = {}
        for system_name in region_weights.index:
            projected = colour_tools.colour_model(model.qso, (f'{system_name}:{first_band}', f'{system_name}:{second_band}'))
            projected_models[system_name] = projected
            model_medians[system_name] = np.array([colour_tools.colour_quantile(colour_tools.mixture_at_z(projected, redshift), 0.5) for redshift in model_redshift])
        x_values = np.r_[panel[component_columns].to_numpy().ravel(), *model_medians.values()]
        padding = max(0.06 * np.ptp(x_values), 0.1)
        colour_limits = (x_values.min() - padding, x_values.max() + padding)
        model_colour = np.linspace(*colour_limits, colour_grid_points)
        model_density = sum((
            weight * np.exp(projected_models[system_name].log_p_colour_given_z(model_colour[:, None], None, model_redshift)) for system_name,
            weight in region_weights.items()
        ))
        resolved_colour_panels.append({
            'colour_name': colour_name,
            'axis_label': f'${first_band}-{second_band}$ (mag)',
            'component_columns': component_columns,
            'candidates': panel,
            'colour_limits': colour_limits,
            'colour_grid': model_colour,
            'redshift_grid': model_redshift,
            'model_density': model_density,
            'model_medians': model_medians
        })
        print(f'{first_band}−{second_band}: {len(panel)}/{len(plot_rows)} scored pairs plotted; {(~finite_colours).sum()} without both component colours; {(finite_colours & ~within_redshift_limits).sum()} outside redshift limits.')
    print(f"Model: {model_bundles['magdep']['bundle']}; curve support: {model_redshift_min:.2f}–{model_redshift_max:.2f}; {finite_score.sum()}/{len(ranked)} finite scores; {(finite_score & ~finite_redshift).sum()} missing redshifts.")

pm_score_candidates = pd.read_csv(ONE_ROOT / 'threshold_fit_candidates.csv', float_precision='round_trip')
fit = np.load(ONE_ROOT / 'threshold_gaussians.npz')
gaussian_means, gaussian_covariances = fit['means'], fit['covariances']
gaussian_component_order = fit['order']
contour_probability = calibration['contour_probability']
beta = calibration['beta']
threshold_grid_points = calibration['grid_points']
tradeoff = pd.read_csv(ONE_ROOT / 'threshold_tradeoff.csv', float_precision='round_trip')
score_thresholds = tradeoff.log_r_threshold.to_numpy()
f_beta_score = tradeoff.f_beta_score.to_numpy()
best_threshold_index = int(np.argmax(f_beta_score))
recommended_log_r_threshold = float(score_thresholds[best_threshold_index])
from scipy.stats import chi2
from matplotlib.lines import Line2D

contour_radius = np.sqrt(chi2.ppf(contour_probability, df=2))
green_component, red_component = gaussian_component_order
score_sigma = np.sqrt(gaussian_covariances[:, 1, 1])
green_score_minimum = gaussian_means[green_component, 1] - contour_radius * score_sigma[green_component]
red_score_maximum = gaussian_means[red_component, 1] + contour_radius * score_sigma[red_component]

# Evaluate responsibilities from the saved full Gaussians; no refitting in plotting.
from scipy.stats import multivariate_normal
pm_score_coordinates = pm_score_candidates[['secondary_pm_significance', 'either_model_log_r_per_unit_z']].to_numpy(float)
gaussian_densities = np.column_stack([
    weight * multivariate_normal.pdf(pm_score_coordinates, mean=mean, cov=covariance)
    for weight, mean, covariance in zip(fit['weights'], gaussian_means, gaussian_covariances)])
gaussian_membership_probabilities = gaussian_densities / gaussian_densities.sum(axis=1, keepdims=True)
# Source cell resolved-score-pm-tradeoff-multiplot
from matplotlib.patches import Ellipse
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import norm
figure_size = (9, 7)
curve_panel_width = 0.28
panel_gap = 0.01
interest_colour = 'gold'
interest_alpha = 0.12
component_colours = ('tab:blue', 'tab:red')
threshold_colour = 'black'
marker_size = 32
pm_score_membership_probabilities = gaussian_membership_probabilities
membership_cmap = LinearSegmentedColormap.from_list('pm_group_membership', [component_colours[1], '0.85', component_colours[0]])
fig, (pm_axis, curve_axis) = plt.subplots(
    1,
    2,
    figsize=figure_size,
    sharey=True,
    layout='constrained',
    gridspec_kw={'width_ratios': (1, curve_panel_width), 'wspace': panel_gap}
)
fig.get_layout_engine().set(w_pad=0.02, wspace=0)
pm_axis.set_box_aspect(1)
pm_axis.set_anchor('E')
curve_axis.set_box_aspect(1 / curve_panel_width)
curve_axis.set_anchor('W')
points = pm_axis.scatter(
    pm_score_candidates.secondary_pm_significance,
    pm_score_candidates.either_model_log_r_per_unit_z,
    c=pm_score_membership_probabilities[:, gaussian_component_order[0]],
    cmap=membership_cmap,
    vmin=0,
    vmax=1,
    s=marker_size,
    edgecolors='black',
    linewidths=0.25,
    label=f'Candidate quasar sample',
    zorder=3
)
for component, colour in zip(gaussian_component_order, component_colours):
    mean = gaussian_means[component]
    eigenvalues, eigenvectors = np.linalg.eigh(gaussian_covariances[component])
    radii = contour_radius * np.sqrt(eigenvalues)
    angle_deg = np.degrees(np.arctan2(eigenvectors[1, 1], eigenvectors[0, 1]))
    pm_axis.add_patch(Ellipse(
        mean,
        width=2 * radii[1],
        height=2 * radii[0],
        angle=angle_deg,
        fill=False,
        edgecolor=colour,
        linewidth=1.8
    ))
    pm_axis.scatter(*mean, marker='x', color=colour, s=80, linewidths=2, zorder=4)
displayed_score_limits = pm_axis.get_ylim()
panel_score_thresholds = np.linspace(*displayed_score_limits, threshold_grid_points)
panel_completeness = norm.sf(
    panel_score_thresholds,
    loc=gaussian_means[green_component, 1],
    scale=score_sigma[green_component]
)
panel_red_retained = norm.sf(
    panel_score_thresholds,
    loc=gaussian_means[red_component, 1],
    scale=score_sigma[red_component]
)
panel_purity = fit['weights'][green_component] * panel_completeness / (fit['weights'][green_component] * panel_completeness + fit['weights'][red_component] * panel_red_retained)
panel_f_beta = (1 + beta ** 2) * panel_purity * panel_completeness / (beta ** 2 * panel_purity + panel_completeness)
curve_axis.plot(
    panel_completeness,
    panel_score_thresholds,
    color='black',
    linestyle='-',
    linewidth=2,
    label=f'$C\\:(t)$'
)
curve_axis.plot(
    panel_purity,
    panel_score_thresholds,
    color='black',
    linestyle='--',
    linewidth=2,
    label=f'$P\\:(t)$'
)
curve_axis.plot(
    panel_f_beta,
    panel_score_thresholds,
    color='black',
    linestyle=':',
    linewidth=2,
    label=f'$S\\:(t)$'
)
curve_axis.set_ylim(displayed_score_limits)
curve_axis.scatter(
    f_beta_score[best_threshold_index],
    recommended_log_r_threshold,
    color=threshold_colour,
    s=35,
    zorder=5
)
for axis in (pm_axis, curve_axis):
    axis.axhspan(
        green_score_minimum,
        red_score_maximum,
        color=interest_colour,
        alpha=interest_alpha,
        zorder=0
    )
    axis.axhline(
        green_score_minimum,
        color=component_colours[0],
        linestyle='--',
        linewidth=1.2,
        label=f'Low-PM contour min. {green_score_minimum:.2f}'
    )
    axis.axhline(
        red_score_maximum,
        color=component_colours[1],
        linestyle='--',
        linewidth=1.2,
        label=f'High-PM contour max. {red_score_maximum:.2f}'
    )
    axis.axhline(
        recommended_log_r_threshold,
        color=threshold_colour,
        linestyle=':',
        linewidth=1.8,
        label=f'Maximum $F_1$ at {recommended_log_r_threshold:.2f}'
    )
    axis.grid(alpha=0.15)
pm_axis.set(xlabel='Secondary $S_{PM}$', ylabel='Maximum eligible $\\ln\\mathcal{R}$')
pm_axis.legend(loc='upper right', fontsize=11)
curve_axis.set(xlabel='Model\nestimate', xlim=(0, 1.05), xticks=[0, 0.5, 1])
curve_axis.tick_params(axis='y', left=False, labelleft=False)
curve_axis.legend(
    handles=curve_axis.lines[:3],
    loc='upper left',
    fontsize=11,
    handlelength=2,
    labelspacing=0.3
)
membership_colourbar = fig.colorbar(
    points,
    cax=curve_axis.inset_axes([1.12, 0, 0.12, 1]),
    orientation='vertical',
    ticks=[0, 0.25, 0.5, 0.75, 1],
    label='Low-PM component probability'
)
save_paper_figure(fig, 'resolved_pm_calibration')
plt.close("all")


# Source cell def252f9
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
resolved_colour_log_r_minimum = log_r_per_unit_z_threshold
resolved_colour_rows = plot_rows.loc[plot_rows.log_r_per_unit_z.ge(resolved_colour_log_r_minimum)]
score_scale = Normalize(
    resolved_colour_rows.log_r_per_unit_z.min(),
    resolved_colour_rows.log_r_per_unit_z.max()
)
redshift_min = resolved_colour_rows.primary_qso_z.min()
redshift_max = resolved_colour_rows.primary_qso_z.max()
# Leave space for markers centred on the lowest and highest redshifts.
redshift_padding = 0.02 * (redshift_max - redshift_min)
redshift_limits = (redshift_min - redshift_padding, redshift_max + redshift_padding)
fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
for ax, colour_panel in zip(axes, resolved_colour_panels):
    panel = colour_panel['candidates']
    panel = panel.loc[panel.index.isin(resolved_colour_rows.index)]
    primary_colour, secondary_colour = colour_panel['component_columns']
    density = colour_panel['model_density']
    ax.pcolormesh(
        colour_panel['colour_grid'],
        colour_panel['redshift_grid'],
        density.T,
        norm=LogNorm(density.max() * 0.001, density.max()),
        cmap='Greys',
        alpha=0.3,
        rasterized=True
    )
    ax.plot(
        colour_panel['model_medians']['decals_dr9_north'],
        colour_panel['redshift_grid'],
        color='firebrick',
        linewidth=1.8,
        linestyle='--',
        zorder=4
    )
    for row in panel.itertuples():
        ax.plot(
            [getattr(row, primary_colour), getattr(row, secondary_colour)],
            [row.primary_qso_z, row.primary_qso_z],
            color='black',
            linewidth=0.6,
            alpha=0.4
        )
    ax.scatter(
        panel[primary_colour],
        panel.primary_qso_z,
        c=panel.log_r_per_unit_z,
        cmap='viridis',
        norm=score_scale,
        marker='o',
        s=36,
        edgecolors='black',
        linewidths=0.5,
        zorder=3
    )
    ax.scatter(
        panel[secondary_colour],
        panel.primary_qso_z,
        c=panel.log_r_per_unit_z,
        cmap='viridis',
        norm=score_scale,
        marker='^',
        s=36,
        edgecolors='black',
        linewidths=0.5,
        zorder=3
    )
    ax.set(xlim=colour_panel['colour_limits'], ylim=redshift_limits)
    ax.set_box_aspect(1)
    first_band, second_band = colour_panel['colour_name'].split('_minus_')
    ax.set_xlabel(f'${first_band}-{second_band}$ (mag)', fontsize=18)
    ax.tick_params(axis='both', labelsize=15)
    ax.grid(alpha=0.15)
axes[0].set_ylabel('Primary redshift, $z_0$', fontsize=18)
colourbar = fig.colorbar(
    ScalarMappable(norm=score_scale, cmap='viridis'),
    cax=axes[-1].inset_axes([1.04, 0, 0.05, 1])
)
colourbar.set_label('Maximum eligible $\\ln\\mathcal{R}$', fontsize=18)
colourbar.ax.tick_params(labelsize=15)
fig.legend(
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
fig.subplots_adjust(left=0.1, right=0.88, bottom=0.22, top=0.96, wspace=0.1)
save_paper_figure(fig, 'resolved_colour_redshift')
plt.close("all")
print(f'Original-subset colour diagnostic, before secondary PM and visual exclusions: log_r_per_unit_z >= {resolved_colour_log_r_minimum:.2f}; {len(resolved_colour_rows)}/{len(plot_rows)} scored systems retained.')


# Source cell resolved-log-r-threshold-diagnostics
if plot_rows.empty:
    print('No eligible resolved candidates for this diagnostic.')
else:
    panel_columns = [
        ('primary_qso_z', 'primary_abs_galactic_b_deg'),
        ('primary_abs_galactic_b_deg', 'ref_mag'),
        ('primary_qso_z', 'ref_mag'),
        ('secondary_r_snr', 'ref_mag'),
        ('primary_r_snr', 'ref_mag'),
        ('qso_ood_sigma', 'primary_qso_z'),
        ('secondary_pm_significance', 'primary_abs_galactic_b_deg'),
        ('component_separation_arcsec', 'ref_mag'),
        ('abs_delta_g_minus_r', 'primary_qso_z'),
        ('abs_delta_r_minus_z', 'primary_qso_z')
    ]
    panels_per_row = 3
    log_r_reference = -20
    x_scales = {'secondary_r_snr': 'log', 'primary_r_snr': 'log', 'secondary_pm_significance': 'log'}
    from matplotlib.colors import Normalize
    threshold_candidates = plot_rows.copy()
    finite_primary_position = np.isfinite(threshold_candidates[['primary_component_ra', 'primary_component_dec']]).all(axis=1)
    primary_positions = SkyCoord(
        ra=threshold_candidates.loc[finite_primary_position, 'primary_component_ra'].to_numpy(float) * u.deg,
        dec=threshold_candidates.loc[finite_primary_position, 'primary_component_dec'].to_numpy(float) * u.deg
    )
    threshold_candidates['primary_abs_galactic_b_deg'] = np.nan
    threshold_candidates.loc[finite_primary_position, 'primary_abs_galactic_b_deg'] = np.abs(primary_positions.galactic.b.deg)
    for component in ('primary', 'secondary'):
        r_flux = threshold_candidates[f'{component}_decals_flux_r']
        r_ivar = threshold_candidates[f'{component}_decals_flux_ivar_r']
        measured_r = np.isfinite(r_flux) & np.isfinite(r_ivar) & r_ivar.gt(0)
        threshold_candidates[f'{component}_r_snr'] = r_flux * np.sqrt(r_ivar.where(measured_r))
    finite_component_positions = np.isfinite(threshold_candidates[[
        'primary_component_ra',
        'primary_component_dec',
        'secondary_component_ra',
        'secondary_component_dec'
    ]]).all(axis=1)
    primary_pair_positions = SkyCoord(
        ra=threshold_candidates.loc[finite_component_positions, 'primary_component_ra'].to_numpy(float) * u.deg,
        dec=threshold_candidates.loc[finite_component_positions, 'primary_component_dec'].to_numpy(float) * u.deg
    )
    secondary_pair_positions = SkyCoord(
        ra=threshold_candidates.loc[finite_component_positions, 'secondary_component_ra'].to_numpy(float) * u.deg,
        dec=threshold_candidates.loc[finite_component_positions, 'secondary_component_dec'].to_numpy(float) * u.deg
    )
    threshold_candidates['component_separation_arcsec'] = np.nan
    threshold_candidates.loc[finite_component_positions, 'component_separation_arcsec'] = primary_pair_positions.separation(secondary_pair_positions).arcsec
    for colour_name in ('g_minus_r', 'r_minus_z'):
        threshold_candidates[f'abs_delta_{colour_name}'] = np.abs(threshold_candidates[f'primary_{colour_name}'] - threshold_candidates[f'secondary_{colour_name}'])
    threshold_labels = {
        'primary_qso_z': 'Primary redshift',
        'primary_abs_galactic_b_deg': 'Primary Galactic |b| (deg)',
        'ref_mag': 'Secondary r luptitude (mag)',
        'secondary_r_snr': 'Secondary r S/N',
        'primary_r_snr': 'Primary r S/N',
        'qso_ood_sigma': 'Quasar model colour distance',
        'secondary_pm_significance': 'Secondary PM significance',
        'component_separation_arcsec': 'Component separation (arcsec)',
        'abs_delta_g_minus_r': '$|\\Delta(g-r)_{\\mathrm{p-s}}|$ (mag)',
        'abs_delta_r_minus_z': '$|\\Delta(r-z)_{\\mathrm{p-s}}|$ (mag)'
    }
    threshold_panel_specs = [
        (
            x_column,
            colour_column,
            x_scales.get(x_column, 'linear'),
            threshold_labels.get(x_column, x_column)
        ) for x_column,
        colour_column in panel_columns
    ]
    threshold_colour_scales = {}
    for colour_column in dict.fromkeys((spec[1] for spec in threshold_panel_specs)):
        colour_values = threshold_candidates[colour_column].to_numpy(float)
        finite_colour_values = colour_values[np.isfinite(colour_values)]
        threshold_colour_scales[colour_column] = Normalize(vmin=finite_colour_values.min(), vmax=finite_colour_values.max())
    finite_threshold_score = np.isfinite(threshold_candidates['log_r_per_unit_z'])
    threshold_score_values = threshold_candidates.loc[finite_threshold_score, 'log_r_per_unit_z']
    score_padding = max(0.05 * np.ptp(threshold_score_values), 0.5)
    threshold_score_limits = (
        threshold_score_values.min() - score_padding,
        threshold_score_values.max() + score_padding
    )
    number_columns = min(panels_per_row, len(threshold_panel_specs))
    number_rows = int(np.ceil(len(threshold_panel_specs) / number_columns))
    fig, axes = plt.subplots(
        number_rows,
        number_columns,
        figsize=(6.0 * number_columns, 5 * number_rows),
        sharey=True,
        squeeze=False,
        layout='constrained'
    )
    for ax, (x_column, colour_column, x_scale, title) in zip(axes.flat, threshold_panel_specs):
        finite_parameters = np.isfinite(threshold_candidates[[x_column, 'log_r_per_unit_z', colour_column]].to_numpy(float)).all(axis=1)
        if x_scale == 'log':
            finite_parameters &= threshold_candidates[x_column].gt(0).to_numpy()
        panel = threshold_candidates.loc[finite_parameters]
        points = ax.scatter(
            panel[x_column],
            panel['log_r_per_unit_z'],
            c=panel[colour_column],
            cmap='viridis',
            norm=threshold_colour_scales[colour_column],
            s=25,
            alpha=0.75,
            edgecolors='black',
            linewidths=0.25
        )
        ax.axhline(log_r_reference, color='black', linestyle='--', linewidth=1)
        ax.set_xscale(x_scale)
        ax.set(
            xlabel=threshold_labels.get(x_column, x_column),
            ylabel='log_r_per_unit_z',
            title=f'{title} (n={len(panel)}/{len(threshold_candidates)})',
            ylim=threshold_score_limits
        )
        ax.set_box_aspect(1)
        ax.grid(alpha=0.2)
        colourbar = fig.colorbar(points, ax=ax, fraction=0.045, pad=0.025)
        colourbar.set_label(threshold_labels.get(colour_column, colour_column))
    for unused_ax in axes.flat[len(threshold_panel_specs):]:
        unused_ax.set_visible(False)
    fig.suptitle(f'Resolved-candidate cutoff diagnostics — dashed line: {log_r_reference:g}')
    save_paper_figure(fig, 'resolved_diagnostics')
    plt.close("all")
