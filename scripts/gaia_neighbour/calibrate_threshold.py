import json
import numpy as np
import pandas as pd
from dual_agn.settings import ONE_ROOT, OUTPUT_ROOT, THRESHOLD_CALIBRATION_PATH, IDENTIFIER_DTYPES
from dual_agn.score_selection import eligible_model_scores
from sklearn.mixture import GaussianMixture
from scipy.stats import chi2, norm
import argparse
parser = argparse.ArgumentParser(description='Calibrate the shared resolved-score threshold.')
parser.add_argument('--previous-threshold', type=float, default=-6.70)
arguments = parser.parse_args()
ranked = pd.read_csv(ONE_ROOT / 'ranked_resolved.csv', dtype=IDENTIFIER_DTYPES, float_precision='round_trip')
model_scores = eligible_model_scores(ranked, 'eligible', 'log_r_per_unit_z')
ranked['either_model_log_r_per_unit_z'] = model_scores.max(axis=1)
finite_score = np.isfinite(ranked.either_model_log_r_per_unit_z)
finite_redshift = np.isfinite(ranked.primary_qso_z)
plot_rows = ranked.loc[finite_score & finite_redshift].copy()
pm_score_columns = ['secondary_pm_significance', 'either_model_log_r_per_unit_z', 'primary_qso_z']
secondary_pm_significance_max = 70.0
finite_pm_score_redshift = np.isfinite(plot_rows[pm_score_columns]).all(axis=1)
pm_significance_pass = plot_rows.secondary_pm_significance.lt(secondary_pm_significance_max)
pm_score_candidates = plot_rows.loc[finite_pm_score_redshift & pm_significance_pass].copy()
print(f'{len(pm_score_candidates)} systems plotted and fitted; {(~finite_pm_score_redshift).sum()} missing measurements; {(finite_pm_score_redshift & ~pm_significance_pass).sum()} at/above PM significance {secondary_pm_significance_max:g} excluded.')
pm_score_coordinates = pm_score_candidates[pm_score_columns[:2]].to_numpy(float)
coordinate_mean = pm_score_coordinates.mean(axis=0)
coordinate_scale = pm_score_coordinates.std(axis=0)
if len(pm_score_candidates) < 6 or not (np.isfinite(coordinate_scale) & (coordinate_scale > 0)).all():
    raise ValueError('Threshold calibration needs at least six finite eligible resolved systems with varying score and PM significance below 70; the supplied catalogue cannot support this fit')
pm_score_gaussian_fit = GaussianMixture(
    n_components=2,
    covariance_type='full',
    n_init=32,
    random_state=42,
    reg_covar=1e-06,
    tol=1e-07,
    max_iter=1000
).fit((pm_score_coordinates - coordinate_mean) / coordinate_scale)
gaussian_means = coordinate_mean + pm_score_gaussian_fit.means_ * coordinate_scale
gaussian_covariances = pm_score_gaussian_fit.covariances_ * coordinate_scale[None, :, None] * coordinate_scale[None, None, :]
gaussian_component_order = np.argsort(gaussian_means[:, 0])
pm_score_component_rows = []
for component in gaussian_component_order:
    mean = gaussian_means[component]
    covariance = gaussian_covariances[component]
    sigma = np.sqrt(np.diag(covariance))
    pm_score_component_rows.append({
        'weight': pm_score_gaussian_fit.weights_[component],
        'mean_pm_significance': mean[0],
        'mean_log_r_per_unit_z': mean[1],
        'sigma_pm_significance': sigma[0],
        'sigma_log_r_per_unit_z': sigma[1],
        'correlation': covariance[0, 1] / np.prod(sigma)
    })
pm_score_gaussian_components = pd.DataFrame(pm_score_component_rows)
contour_probability = 0.68
contour_radius = np.sqrt(chi2.ppf(contour_probability, df=2))
green_component, red_component = gaussian_component_order
score_sigma = np.sqrt(gaussian_covariances[:, 1, 1])
green_score_minimum = gaussian_means[green_component, 1] - contour_radius * score_sigma[green_component]
red_score_maximum = gaussian_means[red_component, 1] + contour_radius * score_sigma[red_component]
beta = 1.0
threshold_grid_points = 2001
score_thresholds = np.linspace(green_score_minimum, red_score_maximum, threshold_grid_points)
green_component, red_component = gaussian_component_order
score_sigma = np.sqrt(gaussian_covariances[:, 1, 1])
green_retained = norm.sf(
    score_thresholds,
    loc=gaussian_means[green_component, 1],
    scale=score_sigma[green_component]
)
red_retained = norm.sf(
    score_thresholds,
    loc=gaussian_means[red_component, 1],
    scale=score_sigma[red_component]
)
green_selected_mass = pm_score_gaussian_fit.weights_[green_component] * green_retained
red_selected_mass = pm_score_gaussian_fit.weights_[red_component] * red_retained
estimated_completeness = green_retained
estimated_purity = green_selected_mass / (green_selected_mass + red_selected_mass)
f_beta_score = (1 + beta ** 2) * estimated_purity * estimated_completeness / (beta ** 2 * estimated_purity + estimated_completeness)
best_threshold_index = int(np.argmax(f_beta_score))
recommended_log_r_threshold = float(score_thresholds[best_threshold_index])
threshold_tradeoff = pd.DataFrame({
    'log_r_threshold': score_thresholds,
    'estimated_completeness': estimated_completeness,
    'estimated_purity': estimated_purity,
    'f_beta_score': f_beta_score,
    'observed_systems_retained': (pm_score_candidates.either_model_log_r_per_unit_z.to_numpy()[:, None] >= score_thresholds[None, :]).sum(axis=0)
})
score_limit_summary = threshold_tradeoff.iloc[[0, best_threshold_index, len(score_thresholds) - 1]].copy()
score_limit_summary.insert(0, 'choice', ['Completeness boundary', f'Maximum F{beta:g}', 'Purity boundary'])
print(f'For the fitted {len(pm_score_candidates)}-system sample, maximum F{beta:g} recommends either_model_log_r_per_unit_z >= {recommended_log_r_threshold:.2f}.')
print(f'Search interval: {green_score_minimum:.2f} to {red_score_maximum:.2f}; grid spacing {score_thresholds[1] - score_thresholds[0]:.4f}.')
if best_threshold_index in (0, len(score_thresholds) - 1):
    print('The optimum is at a search boundary; it is constrained by the requested interval.')
print('These estimates assume green = desired objects and red = contaminants; there are no truth labels here.')
print('They apply to the fitted PM-measured sample and have not been validated as catalogue purity/completeness.')
log_r_per_unit_z_threshold = round(recommended_log_r_threshold, 2)
adopted_retained_fractions = norm.sf(log_r_per_unit_z_threshold, loc=gaussian_means[:, 1], scale=score_sigma)
adopted_selected_mass = pm_score_gaussian_fit.weights_ * adopted_retained_fractions
threshold_calibration = {
    'recommended_log_r_threshold': recommended_log_r_threshold,
    'adopted_log_r_threshold': log_r_per_unit_z_threshold,
    'previous_adopted_log_r_threshold': arguments.previous_threshold,
    'fit_systems': len(pm_score_candidates),
    'pm_significance_upper_limit': secondary_pm_significance_max,
    'contour_probability': contour_probability,
    'beta': beta,
    'maximum_f_beta': float(f_beta_score[best_threshold_index]),
    'adopted_low_pm_proxy_retention': float(adopted_retained_fractions[green_component]),
    'adopted_low_pm_proxy_share': float(adopted_selected_mass[green_component] / adopted_selected_mass.sum()),
    'adopted_fit_systems_retained': int(pm_score_candidates.either_model_log_r_per_unit_z.ge(log_r_per_unit_z_threshold).sum()),
    'search_interval': [float(green_score_minimum), float(red_score_maximum)],
    'grid_points': threshold_grid_points,
    'model': 'either',
    'calibration_statistic': 'maximum finite eligible log R per unit z across both models',
    'sample_unit': 'anchor-entry',
    'gaussian_population_labels_are_proxies': True
}
ONE_ROOT.mkdir(parents=True, exist_ok=True)
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
THRESHOLD_CALIBRATION_PATH.write_text(json.dumps(threshold_calibration, indent=2) + '\n')
threshold_tradeoff.to_csv(ONE_ROOT / 'threshold_tradeoff.csv', index=False)
pm_score_candidates.to_csv(ONE_ROOT / 'threshold_fit_candidates.csv', index=False)
print(f'Adopted threshold (rounded to two decimals): {log_r_per_unit_z_threshold:.2f}')
np.savez(
    ONE_ROOT / 'threshold_gaussians.npz',
    means=gaussian_means,
    covariances=gaussian_covariances,
    weights=pm_score_gaussian_fit.weights_,
    order=gaussian_component_order
)
pm_score_gaussian_components.to_csv(ONE_ROOT / 'threshold_gaussian_components.csv', index=False)
print('Output:', THRESHOLD_CALIBRATION_PATH)
