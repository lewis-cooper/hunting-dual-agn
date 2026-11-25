"""Score components and blended hypotheses with per-object redshift uncertainty."""

from dataclasses import replace

import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm


def score_components_with_redshift_errors(model, photometry, redshift_errors, **arguments):
    """Use a scalar RedshiftMatch error for each uncertain object."""
    errors = np.asarray(redshift_errors, dtype=float)
    if not np.isfinite(errors).all() or np.any(errors < 0):
        raise ValueError('Redshift errors must be finite and non-negative')
    scores = [None] * len(errors)
    eligible = np.zeros(len(errors), dtype=bool)
    reasons = np.full(len(errors), "", dtype=object)
    percentiles = np.full(len(errors), np.nan)
    groups = [(0.0, np.flatnonzero(errors <= 0))]
    groups.extend((float(error), np.array([index]))
                  for index, error in enumerate(errors) if error > 0)
    row_arguments = ("ra_deg", "dec_deg", "z_primary", "morphology",
                     "separation_arcsec", "fracflux", "candidate_id", "primary_id")
    for error, indices in groups:
        if not len(indices):
            continue
        subset_arguments = arguments.copy()
        for name in row_arguments:
            subset_arguments[name] = np.asarray(arguments[name])[indices]
        subset_arguments["match"] = replace(arguments["match"], z_primary_err=error)
        subset_scores, decision = model.score(photometry.subset(indices), **subset_arguments)
        for index, score in zip(indices, subset_scores, strict=True):
            scores[index] = score
        eligible[indices] = decision["eligible"]
        reasons[indices] = decision["reason"]
        percentiles[indices] = decision["qso_support"]["percentile"]
    return scores, {"eligible": eligible, "reason": reasons,
                    "qso_support": {"percentile": percentiles}}


def blend_score_with_redshift_error(model, photometry, ra, dec, z0, zerr, n=9):
    """Average QQ and QS likelihoods over a Gaussian-weighted clipped redshift grid."""
    central = model.score(photometry, ra_deg=[ra], dec_deg=[dec], z_qso=[z0])[0]
    if not zerr > 0:
        return central, {
            "redshift_error": 0.0, "grid_nodes": 1, "finite_nodes": int(
                np.isfinite(central.log_like_qq) and np.isfinite(central.log_like_qs)),
            "unique_nodes": 1, "clipped_nodes": 0, "retained_weight": 1.0,
        }, []

    unbounded_z = np.linspace(z0 - 3 * zerr, z0 + 3 * zerr, n)
    redshifts = np.clip(unbounded_z, 0.15, 4.35)
    logweights = norm.logpdf(redshifts, z0, zerr)
    logweights -= logsumexp(logweights)
    grid_scores = []
    for redshift in redshifts:
        if redshift == z0:
            score = central
        else:
            score = model.score(photometry, ra_deg=[ra], dec_deg=[dec], z_qso=[redshift])[0]
        grid_scores.append(score)
    logqq = np.array([score.log_like_qq for score in grid_scores])
    logqs = np.array([score.log_like_qs for score in grid_scores])
    finite = np.isfinite(logqq) & np.isfinite(logqs)
    qq_mean = float(logsumexp(logqq[finite] + logweights[finite]))
    qs_mean = float(logsumexp(logqs[finite] + logweights[finite]))
    with np.errstate(invalid="ignore"):
        logbf = qq_mean - qs_mean
    with np.errstate(invalid='ignore'):
        posterior = float(np.exp(-np.logaddexp(0.0, -logbf)))
    marginal = replace(
        central, log_like_qq=qq_mean, log_like_qs=qs_mean, log_bf_qq_qs=logbf, p_qq=posterior,
        notes=tuple(central.notes) + ("quaia_redshift_error_nine_point_average",),
    )
    audit = {
        "redshift_error": float(zerr), "grid_nodes": n, "finite_nodes": int(finite.sum()),
        "unique_nodes": int(len(np.unique(redshifts))),
        "clipped_nodes": int(np.count_nonzero(redshifts != unbounded_z)),
        "retained_weight": float(np.exp(logweights[finite]).sum()),
    }
    grid = [
        {"redshift": float(redshift), "weight": float(np.exp(logweight)),
         "log_like_qq": score.log_like_qq, "log_like_qs": score.log_like_qs,
         "status": score.status, "finite_in_both_hypotheses": bool(valid)}
        for redshift, logweight, score, valid in zip(redshifts, logweights, grid_scores, finite, strict=True)
    ]
    return marginal, audit, grid
