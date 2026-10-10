"""Tests for `lopt_perturb`: the systematic uncertainty by perturbation runs.

Run from the LineClass directory:  python -m pytest tests -q

LOPT itself is not run here; a stand-in solver, linear like LOPT, takes
its place where the runs are tested.
"""
import math
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (ROOT, os.path.dirname(ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

import lopt_perturb as L                 # noqa: E402
import make_LOPT_input as M              # noqa: E402

LOW, UPP, OTHER = '059003.000001', '059003.000300', '059003.000200'


def record(wn, low=LOW, upp=UPP, flag='', weight=1.0):
    return M.format_line(wn, 0.05, 100.0, low, upp, flag, weight)


def test_a_record_takes_every_decimal_the_field_holds():
    text = record(32947.173)
    moved = L.with_wn(text, 32947.173 + 0.1234567)
    assert moved[:12] == '32947.296457'
    assert moved[12:] == text[12:]
    assert L.format_wn(121547.155) == '121547.15500'
    assert L.format_wn(9328.777) == '9328.7770000'


def test_the_records_are_read_back(tmp_path):
    path = str(tmp_path / 'in.txt')
    with open(path, 'w', newline='') as fh:
        fh.write(record(20000.0) + '\r\n' + record(19000.5, OTHER, UPP, 'P',
                                                   0.0) + '\r\n')
    recs = L.read_records(path)
    assert [(r[1], r[2], r[3]) for r in recs] == [
        (20000.0, LOW, UPP), (19000.5, OTHER, UPP)]


def test_the_patterns_rebuild_the_covariance():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(5, 5))
    C = X.dot(X.T) * 1e-6
    a, ev = L.patterns(C)
    assert np.allclose(a.T.dot(a), C)
    assert list(ev) == sorted(ev, reverse=True)


def test_a_record_moves_against_its_wavelength_shift():
    recs = [(None, 20000.0, LOW, UPP), (None, 30000.0, OTHER, UPP),
            (None, 25000.0, OTHER, LOW)]
    d = L.record_shifts(recs, ['g1', 'g2', None], {'g1': 0, 'g2': 1},
                        np.array([0.001, -0.002]))
    assert d[0] == pytest.approx(-0.001 * 20000.0 ** 2 * 1e-8)
    assert d[1] == pytest.approx(+0.002 * 30000.0 ** 2 * 1e-8)
    assert d[2] == 0.0


def test_a_record_finds_its_line_through_the_rounding():
    keys = {(LOW, UPP): [(32947.1729, 32947.1485), (32950.0, 32950.02)]}
    groups = {'32947.1485': ('3025-3050', 0.00125)}
    recs = [(None, 32947.173, LOW, UPP), (None, 1.0, OTHER, UPP),
            (None, 32950.0, LOW, UPP)]
    got, missing = L.record_groups(recs, keys, {}, groups)
    assert got == ['3025-3050', None, None]
    assert missing == {'not in the table': 1, 'no group': 1}


def test_files_of_two_calibration_runs_are_refused():
    names = ['a', 'b']
    C = np.diag([1e-6, 4e-6])
    assert L.check_covariance({'1': ('a', 0.001), '2': ('b', 0.002)},
                              names, C) == pytest.approx(0.0)
    with pytest.raises(SystemExit):
        L.check_covariance({'1': ('a', 0.0012)}, names, C)


def test_the_runs_give_each_pattern_s_effect_divided_back():
    """A linear stand-in for LOPT: two levels, each measured by one line,
    E(upper) = the line.  A pattern scaled by f moves the line by f*d, and
    the response divided by f is d whatever f was."""
    def solve(shifts):
        return {LOW: (0.0, 0.0), UPP: (20000.0 + shifts, 0.01)}
    base = solve(0.0)
    d = [-0.003, 0.0004]
    factors = [10.0, 2500.0]
    runs = [solve(f * x) for f, x in zip(factors, d)]
    modes = L.responses(base, runs, factors)
    assert modes[UPP] == pytest.approx(d)
    assert L.ritz_unc(modes, LOW, UPP) == pytest.approx(math.hypot(*d))


def test_the_runs_go_in_parallel_and_all_come_back():
    got = L.run_all(3, {k: k for k in range(7)},
                    runner=lambda d: {'x': (float(d), 0.0)})
    assert sorted(got) == list(range(7))
    assert got[4]['x'][0] == 4.0
