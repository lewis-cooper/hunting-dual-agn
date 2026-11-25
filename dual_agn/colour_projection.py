from __future__ import annotations
import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr
from qso_pcolor.gaussmix import GaussianMixture
from qso_pcolor.qso_model import SlicedColourRedshiftModel

def project_colour(mixture: GaussianMixture, bands: tuple[str, str]) -> GaussianMixture:
    """Intrinsic density of u_a-u_b, normalised over colour in magnitudes.

    The projected variance in mag squared is V_aa+V_bb-2 V_ab. All other
    luptitudes and the common brightness coordinate are marginalised exactly.
    """
    a = np.zeros(mixture.n_dim)
    a[mixture.labels.index(bands[0])] = 1.0
    a[mixture.labels.index(bands[1])] = -1.0
    mean = mixture.means @ a
    variance = np.einsum('i,kij,j->k', a, mixture.covs, a)
    return GaussianMixture(mixture.weights, mean[:, None], variance[:, None, None], (f'{bands[0]}-{bands[1]}',))

def colour_model(model: SlicedColourRedshiftModel, bands: tuple[str, str]) -> SlicedColourRedshiftModel:
    """Project every slice, preserving the fitted model's density interpolation."""
    mixtures = [project_colour(m, bands) for m in model.mixtures]
    return SlicedColourRedshiftModel(model.z_centres, mixtures, model.n_train, model.system, mixtures[0].labels)

def mixture_at_z(model: SlicedColourRedshiftModel, z: float) -> GaussianMixture:
    """The normalised interpolated colour mixture at one supported redshift."""
    centres = model.z_centres
    i = int(np.clip(np.searchsorted(centres, z) - 1, 0, len(centres) - 2))
    t = float(np.clip((z - centres[i]) / (centres[i + 1] - centres[i]), 0, 1))
    if t == 0:
        return model.mixtures[i]
    if t == 1:
        return model.mixtures[i + 1]
    return model.mixtures[i].mixture_with(model.mixtures[i + 1], 1 - t)

def colour_cdf(mixture: GaussianMixture, colour: float) -> float:
    """CDF of a univariate mixture; colour is in magnitudes."""
    sigma = np.sqrt(mixture.covs[:, 0, 0])
    return float(mixture.weights @ ndtr((colour - mixture.means[:, 0]) / sigma))

def colour_quantile(mixture: GaussianMixture, probability: float) -> float:
    """Mixture quantile in magnitudes, including multimodal colour densities."""
    mean, sigma = (mixture.means[:, 0], np.sqrt(mixture.covs[:, 0, 0]))
    return brentq(
        lambda c: colour_cdf(mixture, c) - probability,
        float(np.min(mean - 12 * sigma)),
        float(np.max(mean + 12 * sigma))
    )
