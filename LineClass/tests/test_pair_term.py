"""The close-pair term of the line model (hfs_kappa.pair_uncertainties,
2026-11-21): the pull of a partly resolved neighbor, as an uncertainty."""
import math

import pytest

import hfs_kappa


def width(_wn):
    return 0.2


def test_pull_falls_with_separation_and_vanishes_beyond_xmax():
    near = hfs_kappa.pull_term(0.1, 0.2, 0.5)
    far = hfs_kappa.pull_term(0.3, 0.2, 0.5)
    assert near == pytest.approx(0.1 * 0.5 * hfs_kappa.PULL_G0
                                 * math.exp(-0.5 / hfs_kappa.PULL_X0))
    assert 0.0 < far < near
    assert hfs_kappa.pull_term(hfs_kappa.PULL_XMAX * 0.2, 0.2, 0.5) == 0.0


def test_the_weaker_line_of_a_pair_takes_the_larger_term():
    u = hfs_kappa.pair_uncertainties([(1000.0, 1.0), (1000.1, 3.0)], width)
    assert u[1000.0] == pytest.approx(3 * u[1000.1])
    assert u[1000.0] == pytest.approx(hfs_kappa.pull_term(0.1, 0.2, 0.75))


def test_both_neighbors_count_and_a_lone_line_has_none():
    u = hfs_kappa.pair_uncertainties(
        [(1000.0, 2.0), (1000.1, 2.0), (1000.2, 2.0), (1005.0, 2.0)], width)
    one = hfs_kappa.pull_term(0.1, 0.2, 0.5)
    assert u[1000.1] == pytest.approx(math.sqrt(2) * one)
    assert u[1000.0] == pytest.approx(one)
    assert u[1005.0] == 0.0


def test_missing_intensities_count_as_equal_shares():
    u = hfs_kappa.pair_uncertainties([(1000.0, 0.0), (1000.1, 5.0)], width)
    assert u[1000.0] == pytest.approx(u[1000.1])


def test_line_sigma_adds_the_term_in_quadrature():
    ln = hfs_kappa.Line(wn=20000.0, char='', era=1974, cls='plain',
                        ucls=('', 1974), low='a', upp='b', D=0.0, dJ=1,
                        key=20000.0, pair=0.03)
    bare = ln._replace(pair=0.0)
    ab = (0.002, 0.015)
    u0 = hfs_kappa.line_sigma(bare, ab, ab, 0.0)
    assert hfs_kappa.line_sigma(ln, ab, ab, 0.0) == pytest.approx(
        math.hypot(u0, 0.03))


def test_the_class_fit_takes_the_pair_variance_as_fixed():
    """A class whose residuals are all pair variance fits small terms."""
    import numpy as np
    rng = np.random.default_rng(1)
    vals = [(float(r), 0.0, 20000.0, 0.0, None, 0.05 ** 2)
            for r in rng.normal(0.0, 0.05, 400)]
    a, b = hfs_kappa.fit_two_term(vals, fix_a=True)
    assert b < 0.02
