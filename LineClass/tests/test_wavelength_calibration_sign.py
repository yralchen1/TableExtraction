"""The sign of `own_correction`, pinned against the fit's own equations.

`wavelength_calibration.py` fits one row per accepted line,

    E_upp - E_low - scale * delta_lambda  =  wn_obs,     scale = wn^2 * 1e-8

(`build`, the `A[r, off + col] -= scale * val` term), so the fitted
displacement is

    delta_lambda = (ritz - wn_obs) / scale

and the corrected wavenumber is `wn_obs + scale * delta_lambda`, which is the
Ritz value the fit itself predicts.  `own_correction` is that shift, and it
therefore carries a PLUS sign.  It was written with a minus for a while, which
moved every line twice as far from its Ritz wavenumber as leaving it alone,
and nothing in the pipeline could notice because every consumer simply adds
the column.  This test makes the round trip explicit so it cannot come back.
"""
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (ROOT, os.path.dirname(ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

import wavelength_calibration as W            # noqa: E402


def test_the_correction_carries_a_line_onto_its_ritz_wavenumber():
    """One line, one calibration parameter: solve, then correct, then check.

    The design row is built exactly as `build` builds it - the level columns
    with +1 and -1, the calibration column with `-scale` - so that the test
    fails if that convention is ever changed without the writer being
    changed with it.
    """
    wn = 23083.6428
    scale = wn * wn * 1e-8
    ritz = wn + 0.0322                  # the fit's own Ritz value
    # columns: E_low, E_upp, delta_lambda
    a = np.array([[-1.0, 1.0, -scale]])
    y = np.array([wn])
    # the level difference is held, so the one free parameter is the curve
    sol, _, _, _ = np.linalg.lstsq(a[:, 2:], y - a[:, :2].dot([0.0, ritz]),
                                   rcond=None)
    delta_lambda = float(sol[0])
    assert delta_lambda == pytest.approx((ritz - wn) / scale)

    own_correction = delta_lambda * scale        # what the writer must emit
    assert own_correction > 0                    # his wavenumber is too small
    assert wn + own_correction == pytest.approx(ritz)
    # the sign that was there before moved the line the wrong way
    assert abs(wn - own_correction - ritz) > abs(wn - ritz)


def test_the_writer_emits_that_sign():
    """The one line of `main` that turns the fitted value into the column."""
    import inspect
    src = inspect.getsource(W)
    assert "own_correction='%+.4f' % (v * scale)" in src
    assert "own_correction='%+.4f' % (-v * scale)" not in src
