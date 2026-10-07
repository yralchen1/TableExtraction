"""Tests for eigenvector_check.py: the recursive check of the Cowan
eigenvectors by hyperfine constants and line intensities (2026-10-07).

Run from the LineClass directory:  python -m pytest tests -q
"""
import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import eigenvector_check as E          # noqa: E402


def line(low, upp, x, u=0.1):
    return {'low': low, 'upp': upp, 'x': x, 'u_calc': u}


def branch(upp, lows, const, offsets=None):
    offsets = offsets or [0.0] * len(lows)
    return [line(lo, upp, const + d) for lo, d in zip(lows, offsets)]


def test_population_cancels_within_a_branch():
    """An upper level whose lines all sit one constant off Icalc - its
    population - is good; one whose ratios are wrong is not."""
    lines = (branch('U1', ['G1', 'G2', 'G3'], +1.7)
             + branch('U2', ['G1', 'G2', 'G3'], -0.9, [0.0, +2.5, -2.5]))
    chk = E.Check(lines, 0.3)
    good = chk.run({'G1', 'G2', 'G3'}, ['U1', 'U2'])
    assert good['U1'] == 1 and 'U2' not in good
    n_down, n_up, dof, chi2, _ = chk.test('U2', good)
    assert (n_down, n_up, dof) == (3, 0, 2)
    assert chi2 / dof > E.VETO


def test_a_lower_level_is_tested_through_good_upper_levels():
    """X has no lines of its own going down: it is judged by its line in
    the branch of each good upper level, against that branch's other
    lines."""
    lines = []
    for k, c in enumerate((0.4, -1.2, 2.0)):
        lines += branch('U%d' % k, ['G1', 'G2'], c)
        lines.append(line('X', 'U%d' % k, c + 0.05))
        lines.append(line('Y', 'U%d' % k, c + (3.0 if k % 2 else -3.0)))
    chk = E.Check(lines, 0.3)
    good = chk.run({'G1', 'G2', 'U0', 'U1', 'U2'}, ['X', 'Y'])
    assert 'X' in good and 'Y' not in good
    _, n_up, dof, chi2, mean_up = chk.test('X', good)
    assert n_up == 3 and dof == 3
    assert mean_up == pytest.approx(0.05 / math.hypot(0.1, 0.3) * 1.0,
                                    rel=0.2)


def test_the_recursion_climbs():
    """A level made good in one round lets the next round test another."""
    lines = (branch('A', ['G1', 'G2'], 0.0)
             + branch('B', ['A', 'G1', 'G2'], 1.0))
    good = E.Check(lines, 0.3).run({'G1', 'G2'}, ['A', 'B'])
    # A has one degree of freedom from its own branch: not enough alone,
    # but B's branch tests it from above only once B is good - neither is
    # reached, so nothing is made good from too little
    assert good == {'G1': 0, 'G2': 0}
    lines += [line('G3', 'A', 0.0)]
    good = E.Check(lines, 0.3).run({'G1', 'G2', 'G3'}, ['A', 'B'])
    assert good['A'] == 1 and good['B'] == 2


def test_status():
    assert E.status_of(True, 0, 5, 9.0) == 'good_A'
    assert E.status_of(True, 2, 5, 1.0) == 'good_I'
    assert E.status_of(False, None, 1, 9.0) == 'untested'
    assert E.status_of(False, None, 30, 2.0) == 'doubtful'
    assert E.status_of(False, None, 4, E.VETO + 0.1) == 'fail'


def test_observational_scatter_removes_the_calculated_part():
    lines = []
    for k in range(40):
        lines += branch('U%d' % k, ['L1', 'L2', 'L3', 'L4'], k * 0.1,
                        [-0.6, -0.2, +0.2, +0.6])
    s_obs, sd, u = E.observational_scatter(lines)
    assert u == pytest.approx(0.1)
    assert s_obs == pytest.approx(math.sqrt(sd * sd - 0.01))
