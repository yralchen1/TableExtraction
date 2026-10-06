"""Tests for the kappa safeguard of the plate calibration (2026-10-06): the
jackknife uncertainty of each class's kappa, and the deadband inside which
a class keeps the kappa [hfs.kappa] holds now.

Run from the LineClass directory:  python -m pytest tests -q

The tests fix:

  * `kappa_jackknife`: each line's leave-one-out change of kappa is the one
    a refit without the line gives, and the uncertainty is their spread;
  * `kappa_holds`: a class moves only beyond the deadband, and a class the
    configuration lacks always moves;
  * `build` with a held class: no column, and kappa * D off the right-hand
    side, so the fit with the class held at its own fitted value is the
    free fit.
"""
import math
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import hfs_kappa                           # noqa: E402
import wavelength_calibration as W         # noqa: E402

LEVELS = ['059003.%06d' % k for k in (10, 20, 30, 40, 50)]
E = dict(zip(LEVELS, (0.0, 1000.0, 2500.0, 4000.0, 5200.0)))
KAPPA = {'plain_1974': 0.8, 'c': 0.2}


def network(seed=1):
    """Every pair of five levels, twice: once plain, once c, each with its
    own D, measured at its class's kappa with a little noise."""
    rng = np.random.default_rng(seed)
    lines = []
    for i, low in enumerate(LEVELS):
        for upp in LEVELS[i + 1:]:
            for cls in ('plain', 'c'):
                D = float(rng.uniform(-0.8, 0.8))
                k = KAPPA['c' if cls == 'c' else 'plain_1974']
                wn = (E[upp] - E[low] + 30000.0 + k * D
                      + float(rng.normal(0.0, 0.01)))
                lines.append(hfs_kappa.Line(
                    wn=wn, char='c' if cls == 'c' else '', era=1974,
                    cls=cls, ucls=('', 1974), low=low, upp=upp, D=D,
                    dJ=None))
    # the first level is 30000 lower for every line: shift it once
    return lines


def fit(lines, keep=None, fixed=None):
    n = len(lines)
    keep = [True] * n if keep is None else keep
    parts = W.build(lines, [3000.0] * n, [0] * n, [0] * n, [0.01] * n, keep,
                    0, lambda *a: [], fixed)
    return parts, W.solve(*parts[5:])


def test_the_jackknife_is_the_refit_without_the_line():
    lines = network()
    parts, f = fit(lines)
    jack = W.kappa_jackknife(parts, f, lines, top=len(lines))
    # plain_1969 is a required class with no line here: it gets a column
    # that nothing votes on, and no spread
    assert jack['plain_1969'][:2] == (0.0, 0)
    for c in ('plain_1974', 'c'):
        j = parts[4][c]
        u, n, top = jack[c]
        assert n == len(lines) // 2               # the lines that vote
        d = {wn: dk for wn, dk in top}
        for i, ln in enumerate(lines):
            keep = [k != i for k in range(len(lines))]
            parts2, f2 = fit(lines, keep)
            want = f2['sol'][parts2[4][c]] - f['sol'][j]
            assert d[ln.wn] == pytest.approx(want, abs=1e-9)
        dk = np.array([d[ln.wn] for ln in lines])
        n = len(lines)
        assert u == pytest.approx(
            math.sqrt((n - 1) / n * np.sum((dk - dk.mean()) ** 2)))


def test_the_lines_that_move_a_class_most_come_first():
    lines = network()
    parts, f = fit(lines)
    for c in ('plain_1974', 'c'):
        top = W.kappa_jackknife(parts, f, lines, top=3)[c][2]
        assert len(top) == 3
        assert abs(top[0][1]) >= abs(top[1][1]) >= abs(top[2][1])


def test_a_class_moves_only_beyond_the_deadband():
    previous = {'plain_1974': (0.803, 0.013), 'c': (0.146, 0.051)}
    # 23529.3109 of 2026-10-06: 0.024 against a realistic 0.032 - held
    fitted = {'plain_1974': (0.779, 0.032), 'c': (0.110, 0.045)}
    assert W.kappa_holds(fitted, previous) == {'plain_1974': 0.803,
                                               'c': 0.146}
    # beyond it the class moves (c: 0.056 off against 0.045); a narrower
    # deadband moves a class sooner
    fitted = {'plain_1974': (0.779, 0.032), 'c': (0.090, 0.045)}
    assert W.kappa_holds(fitted, previous) == {'plain_1974': 0.803}
    assert W.kappa_holds(fitted, previous, deadband=0.5) == {}
    # a deadband of 0 takes every fitted value
    assert W.kappa_holds({'c': (0.146, 0.05)}, previous, deadband=0.0) == {}
    # a class the configuration does not have yet always moves
    assert W.kappa_holds({'flag_unlisted': (0.9, 0.025)}, previous) == {}


def test_a_held_class_has_no_column_and_its_shift_on_the_right():
    lines = network()
    parts, f = fit(lines)
    k_c = float(f['sol'][parts[4]['c']])
    held, fh = fit(lines, fixed={'c': k_c})
    assert 'c' not in held[4] and 'plain_1974' in held[4]
    assert held[5].shape[1] == parts[5].shape[1] - 1
    # held at its own fitted value, the rest of the fit does not change
    assert fh['sol'][held[4]['plain_1974']] == pytest.approx(
        f['sol'][parts[4]['plain_1974']], abs=1e-9)
    c_rows = [r for r, i in enumerate(held[0]) if lines[i].cls == 'c']
    for r in c_rows:
        ln = lines[held[0][r]]
        assert held[6][r] == pytest.approx(ln.wn - k_c * ln.D)
