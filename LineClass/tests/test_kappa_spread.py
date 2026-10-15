"""Tests for the spread of kappa (2026-10-15): the third term of the line
uncertainty model, s * |D| per era, fitted in the plate calibration
(hfs_kappa.fit_classes) and carried into u_hfs_shift (hfs_correction).

Run from the LineClass directory:  python -m pytest tests -q

The tests fix:

  * the joint fit recovers a known (a, b) and s from synthetic residuals,
    and returns s = 0 where the residuals do not grow with |D|;
  * which lines carry the term: the fitted kappa classes, not the anchor;
  * the term in a line's u_hfs_shift and in Set 3's center-of-gravity
    uncertainty, and that it is not counted as owed to undetermined levels;
  * the configuration table [hfs.kappa_spread] and its writer.
"""
import math
import os
import sys
import tomllib

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import config                              # noqa: E402
import hfs_correction as H                 # noqa: E402
import hfs_kappa                           # noqa: E402

from test_set3 import (KAPPA, LOW, UPP, OTHER, UNDET, S, a_table, tsv,  # noqa: E402,E501
                       EXC_HEADER, SAT_HEADER)


# ---------------------------------------------------------------------------
# The fit
# ---------------------------------------------------------------------------
def synthetic(n, a, b, s, seed, era_wn=(20000.0, 47000.0)):
    rng = np.random.default_rng(seed)
    wn = rng.uniform(*era_wn, n)
    d = np.abs(rng.normal(0.0, 0.2, n))
    h = rng.uniform(0.0, 0.4, n)
    u = np.sqrt((a * 1e-8 * wn * wn) ** 2 + b * b + (s * d) ** 2)
    r = rng.normal(0.0, 1.0, n) * u * np.sqrt(1.0 - h)
    return [(float(r[i]), float(h[i]), float(wn[i]), float(d[i]))
            for i in range(n)]


def test_the_joint_fit_recovers_the_three_terms():
    groups = {('', 1974): synthetic(4000, 0.0022, 0.016, 0.07, 1)}
    fitted, s = hfs_kappa.fit_classes(groups)
    a, b = fitted[('', 1974)]
    assert a == pytest.approx(0.0022, rel=0.1)
    assert b == pytest.approx(0.016, rel=0.15)
    assert s[1974] == pytest.approx(0.07, rel=0.15)


def test_no_spread_is_fitted_where_the_residuals_do_not_grow_with_D():
    groups = {('', 1974): synthetic(3000, 0.0022, 0.02, 0.0, 2)}
    _, s = hfs_kappa.fit_classes(groups)
    assert s[1974] < 0.02


def test_the_two_term_fit_takes_up_the_spread_in_b():
    """Without the term the constant b absorbs it, which is the symptom."""
    groups = {('', 1974): synthetic(4000, 0.0022, 0.016, 0.07, 3)}
    two, s = hfs_kappa.fit_classes(groups, spread=False)
    three, _ = hfs_kappa.fit_classes(groups)
    assert s == {}
    assert two[('', 1974)][1] > three[('', 1974)][1] * 1.1


def test_the_anchor_class_carries_no_spread():
    ln = hfs_kappa.Line(wn=30000.0, char='*r', era=1974, cls='flag',
                        ucls=('*r', 1974), low='a', upp='b', D=0.5, dJ=1)
    assert hfs_kappa.spread_D(ln) == 0.0
    ln = ln._replace(char='', cls='plain', D=-0.3)
    assert hfs_kappa.spread_D(ln) == pytest.approx(0.3)


def test_a_line_sigma_adds_the_spread_in_quadrature():
    ln = hfs_kappa.Line(wn=30000.0, char='', era=1974, cls='plain',
                        ucls=('', 1974), low='a', upp='b', D=0.4, dJ=1)
    base = hfs_kappa.two_term(0.002, 0.015, 30000.0)
    got = hfs_kappa.line_sigma(ln, (0.002, 0.015), (0.002, 0.015), 0.05)
    assert got == pytest.approx(math.hypot(base, 0.05 * 0.4))


# ---------------------------------------------------------------------------
# The term in the correction
# ---------------------------------------------------------------------------
SPREAD = ((1969, 0.2), (1974, 0.05))


def make_model(tmp_path, spread=SPREAD, exceptions=(), satellites=()):
    return H.Model(config.HfsSettings(
        apply=True, A_levels=a_table(tmp_path), kappa=KAPPA,
        kappa_spread=spread,
        satellites=tsv(tmp_path, 'sat.txt', SAT_HEADER, satellites),
        kappa_exceptions=tsv(tmp_path, 'exc.txt', EXC_HEADER, exceptions)))


def test_a_plain_line_gains_s_D_in_u_hfs_shift(tmp_path):
    with_s = make_model(tmp_path)
    without = make_model(tmp_path, spread=())
    D = S(UPP) - S(LOW)
    comp = [(LOW, UPP, 1.0)]
    _, u0, _ = without.line_shift('', 30000.0, comp)
    shift, u1, _ = with_s.line_shift('', 30000.0, comp)
    assert u1 == pytest.approx(math.hypot(u0, 0.05 * D))
    assert with_s.line_spread('', 30000.0, comp) == pytest.approx(
        abs(0.05 * D))
    # the era decides s
    assert with_s.line_spread('', 60000.0, comp) == pytest.approx(
        abs(0.2 * D))
    # the shift itself is unchanged
    assert shift == pytest.approx(without.line_shift('', 30000.0, comp)[0])


def test_a_listed_flag_and_a_registry_head_carry_none(tmp_path):
    m = make_model(tmp_path)
    assert m.line_spread('*v', 30000.0, [(LOW, UPP, 1.0)]) == 0.0
    m = make_model(tmp_path, satellites=[
        ('29999.7', '30000', LOW, UPP, '1', '', '', '')])
    assert m.line_spread('', 30000.0, [(LOW, UPP, 1.0)]) == 0.0
    m = make_model(tmp_path, exceptions=[
        ('30000', LOW, UPP, 'cg', '', '', '')])
    assert m.line_spread('', 30000.0, [(LOW, UPP, 1.0)]) == 0.0


def test_in_a_flagged_blend_only_the_plain_component_carries_it(tmp_path):
    m = make_model(tmp_path)
    d2 = S(UPP) - S(OTHER)
    got = m.line_spread('*r', 30000.0, [(LOW, UPP, 0.7), (OTHER, UPP, 0.3)])
    assert got == pytest.approx(abs(0.3 * 0.05 * d2))


def test_the_cg_uncertainty_carries_it_and_it_is_not_undetermined(tmp_path):
    with_s = make_model(tmp_path)
    without = make_model(tmp_path, spread=())
    D = S(UPP) - S(LOW)
    comp = [(LOW, UPP, 1.0)]
    g1 = with_s.cg_shift('', 30000.0, comp)
    g0 = without.cg_shift('', 30000.0, comp)
    assert g1.u == pytest.approx(math.hypot(g0.u, 0.05 * D))
    assert g1.u_undet == pytest.approx(g0.u_undet)
    # a line touching an undetermined level: its u_undet is unchanged too
    comp = [(UNDET, UPP, 1.0)]
    assert with_s.cg_shift('', 30000.0, comp).u_undet == pytest.approx(
        without.cg_shift('', 30000.0, comp).u_undet)


def test_combine_takes_five_or_six_tuples():
    five = [(1.0, 0.75, 0.03, 0.2, 0.01)]
    six = [five[0] + (0.1,)]
    s5, u5, k5 = H.combine(five)
    s6, u6, k6 = H.combine(six)
    assert s5 == s6 and k5 == k6
    assert u6 == pytest.approx(math.hypot(u5, 0.1 * 0.2))
    assert H.combine_terms(six) == H.combine_terms(five)


# ---------------------------------------------------------------------------
# The configuration
# ---------------------------------------------------------------------------
PARENT = """\
[hfs]
apply = false

[hfs.kappa]
flag          = [1.0, 0.0]
plain_1974    = [0.751, 0.034]
plain_1969    = [0.440, 0.049]
c             = [0.039, 0.041]

[range]
wn_min = 0.0
"""


def test_the_writer_adds_the_table_and_then_replaces_its_values(tmp_path):
    path = tmp_path / 'k.toml'
    path.write_text(PARENT, encoding='utf-8', newline='\n')
    _, before = config.set_kappa_spread(str(path), {1974: 0.067, 1969: 0.21},
                                        'first')
    assert before == {}
    with open(path, 'rb') as fh:
        raw = tomllib.load(fh)
    assert raw['hfs']['kappa_spread'] == {'1974': 0.067, '1969': 0.21}
    assert raw['range'] == {'wn_min': 0.0}          # nothing else touched
    _, before = config.set_kappa_spread(str(path), {1974: 0.06, 1969: 0.2},
                                        'second')
    assert before == {1974: 0.067, 1969: 0.21}
    text = path.read_text(encoding='utf-8')
    assert text.count('[hfs.kappa_spread]') == 1
    assert 'first' not in text and 'second' in text
    assert text.count('1974 = ') == 1
    # the kappa writer still finds its own table
    config.set_hfs_kappa(str(path), {'plain_1974': (0.75, 0.03),
                                     'plain_1969': (0.44, 0.05),
                                     'c': (0.04, 0.04)})
    with open(path, 'rb') as fh:
        raw = tomllib.load(fh)
    assert raw['hfs']['kappa_spread'] == {'1974': 0.06, '1969': 0.2}


def test_the_project_configuration_loads_the_table():
    cfg = config.load()
    spread = dict(cfg.hfs.kappa_spread)
    assert set(spread) <= {1974, 1969}
    assert all(v >= 0.0 for v in spread.values())
