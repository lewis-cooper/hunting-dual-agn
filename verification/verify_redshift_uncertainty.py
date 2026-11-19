"""Check Gaussian likelihood averaging and per-object redshift-error dispatch."""

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
import sys
import warnings

import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from redshift_uncertainty import (blend_score_with_redshift_error,
                                  score_components_with_redshift_errors)


@dataclass
class Score:
    log_like_qq: float
    log_like_qs: float
    log_bf_qq_qs: float
    p_qq: float = 0.5
    notes: tuple = ()
    status: str = 'ok'


class BlendExample:
    def __init__(self, finite=True):
        self.finite = finite
        self.calls = []

    def score(self, photometry, *, ra_deg, dec_deg, z_qso):
        redshift = z_qso[0]
        self.calls.append(redshift)
        qq = -3 * (redshift - 0.8) ** 2
        qs = -2 * (redshift - 1.7) ** 2
        if not self.finite or redshift < 0.4:
            qq = -np.inf
        return [Score(qq, qs, qq - qs)]


@dataclass
class Match:
    half_width_kms: float = 2000.0
    z_primary_err: float = 0.0


class PhotometryExample:
    def __init__(self, indices):
        self.indices = np.asarray(indices)

    def subset(self, indices):
        return PhotometryExample(self.indices[indices])


class ComponentExample:
    def __init__(self):
        self.calls = []

    def score(self, photometry, **arguments):
        indices = photometry.indices
        error = arguments['match'].z_primary_err
        self.calls.append((indices.tolist(), error))
        np.testing.assert_array_equal(arguments['candidate_id'], indices)
        return [SimpleNamespace(value=int(index) + error) for index in indices], {
            'eligible': indices % 2 == 0,
            'reason': np.array([str(index) for index in indices]),
            'qso_support': {'percentile': indices / 10},
        }


def main():
    model = BlendExample()
    exact, audit, grid = blend_score_with_redshift_error(model, None, 1, 2, 1.0, 0)
    assert len(model.calls) == 1 and not grid
    assert audit['redshift_error'] == 0
    assert exact.log_bf_qq_qs == exact.log_like_qq - exact.log_like_qs

    # This case includes clipped duplicate nodes and nodes invalid in only QQ.
    averaged, audit, grid = blend_score_with_redshift_error(BlendExample(), None, 1, 2, 0.3, 0.3)
    redshifts = np.clip(np.linspace(0.3 - 0.9, 0.3 + 0.9, 9), 0.15, 4.35)
    logweights = norm.logpdf(redshifts, 0.3, 0.3)
    logweights -= logsumexp(logweights)
    qq = -3 * (redshifts - 0.8) ** 2
    qs = -2 * (redshifts - 1.7) ** 2
    valid = redshifts >= 0.4
    expected = logsumexp(qq[valid] + logweights[valid]) - logsumexp(qs[valid] + logweights[valid])
    np.testing.assert_allclose(averaged.log_bf_qq_qs, expected, rtol=0, atol=1e-14)
    assert len(grid) == 9 and audit['unique_nodes'] < 9
    assert audit['finite_nodes'] == int(valid.sum())
    np.testing.assert_allclose(audit['retained_weight'], np.exp(logweights[valid]).sum())

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        invalid, audit, _ = blend_score_with_redshift_error(BlendExample(False), None, 1, 2, 1.0, 0.1)
    assert audit['finite_nodes'] == 0 and np.isnan(invalid.log_bf_qq_qs)

    indices = np.arange(4)
    arguments = {name: indices for name in ('ra_deg', 'dec_deg', 'z_primary', 'morphology',
                 'separation_arcsec', 'fracflux', 'candidate_id', 'primary_id')}
    model = ComponentExample()
    scores, decisions = score_components_with_redshift_errors(
        model, PhotometryExample(indices), [0, 0.2, 0, 0.1], match=Match(), **arguments)
    assert model.calls == [([0, 2], 0.0), ([1], 0.2), ([3], 0.1)]
    np.testing.assert_allclose([score.value for score in scores], [0, 1.2, 2, 3.1])
    np.testing.assert_array_equal(decisions['eligible'], [True, False, True, False])
    np.testing.assert_array_equal(decisions['reason'], ['0', '1', '2', '3'])
    np.testing.assert_allclose(decisions['qso_support']['percentile'], indices / 10)
    print('Passed: exact-redshift dispatch, separate QQ/QS averaging, clipped duplicates,')
    print('finite-node intersection, no-finite-node result, per-object errors and row ordering.')


if __name__ == '__main__':
    main()
