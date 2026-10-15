"""Tests for the position prior of `level_positions`.

Run from the LineClass directory:  python -m pytest tests -q

ln R weighs the recorded lines only.  The position prior adds what the
calculation says about where the level lies: a Student t in
z = (E - E_calc)/W, fitted on the found levels with W taken leave-one-out, and
written as ln_P = ln t(z) - ln t(0), zero at E_calc and never positive.  The
tests fix the leave-one-out z, the fit, the shape of ln_P, and the three
places it enters: the scan of a known level, the table of an unknown one and
the disposition of an alternate.
"""
import math
import os
import sys

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import level_positions as lp          # noqa: E402
import test_level_positions as tlp    # noqa: E402


def enlev(groups, unknown=0):
    """A level list: {cfg: [E_obs - E_calc of its found levels]}."""
    rows = []
    for cfg, omc in groups.items():
        for o in omc:
            rows.append(dict(cfg=cfg, omc=float(o), known=True))
    for _ in range(unknown):
        rows.append(dict(cfg='x', omc=0.0, known=False))
    return pd.DataFrame(rows)


def test_z_is_taken_against_the_other_levels_of_the_configuration():
    en = enlev({'a': [3.0, 4.0, 0.0, -5.0], 'b': [1.0, 2.0]}, unknown=5)
    z = lp.prior_z(en)
    # 'b' has too few found levels; the unfound ones never enter
    assert len(z) == 4
    others = np.array([4.0, 0.0, -5.0])
    assert z[0] == pytest.approx(3.0 / math.sqrt(np.mean(others ** 2)))


def test_the_fit_finds_gaussian_tails_gaussian_and_heavy_tails_heavy():
    rng = np.random.default_rng(1)
    gauss = enlev({f'c{i}': rng.normal(0.0, 50.0, 100) for i in range(10)})
    heavy = enlev({f'c{i}': 50.0 * rng.standard_t(3.0, 100)
                   for i in range(10)})
    pg = lp.fit_position_prior(gauss, log=lambda *a: None)
    ph = lp.fit_position_prior(heavy, log=lambda *a: None)
    assert pg.nu > 20.0
    assert pg.tau == pytest.approx(1.0, abs=0.1)
    assert 2.0 < ph.nu < 5.0


def test_too_few_levels_give_no_prior():
    en = enlev({'a': [1.0, 2.0, 3.0]})
    assert lp.fit_position_prior(en, log=lambda *a: None) is None


def test_ln_p_is_zero_at_e_calc_and_falls_symmetrically():
    pr = lp.PositionPrior(nu=4.0, tau=0.9, n=100)
    e = np.array([1000.0, 1020.0, 980.0, 1060.0])
    v = lp.ln_prior(e, 1000.0, 20.0, pr)
    assert v[0] == 0.0
    assert v[1] == pytest.approx(v[2])
    assert v[3] < v[1] < 0.0
    t = 1.0 / 0.9
    assert v[1] == pytest.approx(-2.5 * math.log1p(t * t / 4.0))


def test_ln_p_is_zero_without_a_prior_or_a_window():
    pr = lp.PositionPrior(nu=4.0, tau=1.0, n=100)
    assert lp.ln_prior([1050.0], 1000.0, 20.0, None)[0] == 0.0
    assert lp.ln_prior([1050.0], np.nan, 20.0, pr)[0] == 0.0
    assert lp.ln_prior([1050.0], 1000.0, np.nan, pr)[0] == 0.0


def test_the_prior_raises_an_alternate_near_e_calc():
    """The position the lines prefer by 34 nats sits ten W from E_calc; the
    one at E_calc becomes its rival only once the prior is counted."""
    ctx = tlp.j_ctx()
    r = lp.scan_level(ctx, lp.UNKNOWN_ID, 100000.0, 100010.0, 1.0, step=0.05)
    assert r['alternates'] == []
    ctx.prior = lp.PositionPrior(nu=1000.0, tau=1.0, n=100)
    r = lp.scan_level(ctx, lp.UNKNOWN_ID, 100000.0, 100010.0, 1.0, step=0.05)
    assert r['ln_P'] < -40.0
    e, v, p = r['alternates'][0]
    assert e == pytest.approx(100010.0, abs=0.05)
    assert p == pytest.approx(0.0, abs=0.01)
    assert v + p > r['ln_R'] + r['ln_P']


def test_the_unknown_table_counts_ln_p_in_ln_r_j():
    ctx = tlp.j_ctx()
    r = dict(E_calc=100000.0, W=20.0,
             positions=[(100000.0, 8.0, -3.0), (100010.0, 6.0, -0.5)])
    tab = lp.unknown_table(ctx, r)
    clean = tab[tab['E'] == 100000.0].iloc[0]
    dirty = tab[tab['E'] == 100010.0].iloc[0]
    assert clean['ln_P'] == pytest.approx(-3.0)
    assert clean['ln_R_J'] == pytest.approx(5.0)
    assert dirty['ln_R_J'] == pytest.approx(5.5 - dirty['ln_J'])
    assert clean['look'] == pytest.approx(5.0 - math.log(2))


def test_a_relocation_must_convince_after_the_prior():
    row = dict(near_dE=np.nan, dE_alt=40.0, top_share=0.3, n_free=6,
               ln_R_alt=4.0, ln_R_alt_J=4.0, ln_R_alt_post=4.0,
               free_gain=12.0, look=5.0, n_own=7, n_kept=0)
    assert lp.disposition(row) == 'relocate'
    row['ln_R_alt_post'] = -0.5
    assert lp.disposition(row) == 'weak'
