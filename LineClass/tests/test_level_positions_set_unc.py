"""Tests for the three level_positions amendments of 2026-10-15.

Run from the LineClass directory:  python -m pytest tests -q

1. One uncertainty model.  The score weighs a line as LOPT does: an accepted
   line by its record in the set's LOPT input (placed through the centroid
   when it is a blend component), any other line by unc_wn_obs and the
   terms make_LOPT_input.hfs_line_unc would add for the candidate.
2. Soft blend credit.  Claims that already account for the whole feature
   are believed only up to about I_obs, smoothly over one intensity width.
3. No action label where the adopted position holds (gain <= 0).
"""
import math
import os
import sys
import tempfile

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import hfs_correction                      # noqa: E402
import level_positions as lp               # noqa: E402
import make_LOPT_input                     # noqa: E402
from test_level_positions import fake_ctx  # noqa: E402


class _Model(object):
    """Stands for hfs_correction.Model where only its presence is asked."""
    kappa_spread = {}


def _hfs_ctx(kap, ukap, head=False, spread=0.0, unc=0.03):
    ctx = fake_ctx([30000.0])
    ctx.unc_o = np.array([unc])
    ctx.hfs = _Model()
    ctx.kap_o, ctx.ukap_o = np.array([kap]), np.array([ukap])
    ctx.head_o = np.array([head])
    ctx.spread_o = np.array([spread])
    return ctx


def _tt(D, u2k, u2u, w2a, w2u, acc_i=-1, acc_wn=0.0, acc_unc=0.0):
    return dict(D=np.array([D]), u2k=np.array([u2k]), u2u=np.array([u2u]),
                w2a=np.array([w2a]), w2u=np.array([w2u]),
                acc_i=np.array([acc_i]), acc_wn=np.array([acc_wn]),
                acc_unc=np.array([acc_unc]))


# ---------------------------------------------------------------------------
# 1. The set's uncertainties
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('kap,ukap,D,uk,uu,wa,wu', [
    (0.75, 0.034, 0.40, 0.03, 0.00, 0.00, 0.00),
    (0.44, 0.049, -0.25, 0.02, 0.05, 0.03, 0.04),
    (0.75, 0.034, 0.10, 0.00, 0.08, 0.00, 0.02),
    (0.90, 0.025, 0.60, 0.01, 0.00, 0.05, 0.00),
])
def test_a_free_line_is_weighed_as_make_LOPT_input_would_weigh_it(
        kap, ukap, D, uk, uu, wa, wu):
    """candidate_terms against make_LOPT_input.hfs_line_unc fed what
    hfs_correction.combine_terms gives for one component of BF = 1."""
    unc = 0.03
    ctx = _hfs_ctx(kap, ukap, unc=unc)
    var, d = lp.candidate_terms(ctx, _tt(D, uk ** 2, uu ** 2, wa ** 2,
                                         wu ** 2),
                                np.array([[0]]), np.array([[30000.0]]))
    shift, u_A, u_k, _ = hfs_correction.combine_terms(
        [(1.0, kap, ukap, D, math.hypot(uk, uu))])
    _, u_A0, u_k0, _ = hfs_correction.combine_terms(
        [(1.0, kap, ukap, D, uk)])
    u_shift = math.hypot(u_A, u_k)
    u_undet = math.sqrt(max(u_shift ** 2 - u_A0 ** 2 - u_k0 ** 2, 0.0))
    ref = make_LOPT_input.hfs_line_unc(unc, [wa] if wa else [],
                                       [wu] if wu else [], u_shift, u_undet,
                                       False)
    assert math.sqrt(float(var[0, 0])) == pytest.approx(ref, rel=1e-9)
    assert float(d[0, 0]) == pytest.approx(shift, abs=1e-8)


def test_a_line_at_the_head_of_its_pattern_takes_neither_shift_nor_widths():
    ctx = _hfs_ctx(1.0, 0.0, head=True)
    var, d = lp.candidate_terms(ctx, _tt(0.5, 0.01, 0.0, 0.04, 0.0),
                                np.array([[0]]), np.array([[30000.0]]))
    assert float(var[0, 0]) == pytest.approx(0.03 ** 2)
    assert float(d[0, 0]) == pytest.approx(0.0)


def test_the_spread_of_kappa_widens_a_wide_pattern():
    ctx = _hfs_ctx(0.75, 0.0, spread=0.2)
    var, _ = lp.candidate_terms(ctx, _tt(0.5, 0.0, 0.0, 0.0, 0.0),
                                np.array([[0]]), np.array([[30000.0]]))
    assert float(var[0, 0]) == pytest.approx(0.03 ** 2 + (0.2 * 0.5) ** 2)


def test_without_the_correction_a_line_takes_the_widths_of_both_levels():
    ctx = fake_ctx([30000.0])
    ctx.unc_o = np.array([0.03])
    var, d = lp.candidate_terms(ctx, _tt(0.0, 0.0, 0.0, 0.05 ** 2, 0.0),
                                np.array([[0]]), np.array([[30000.1]]))
    assert float(var[0, 0]) == pytest.approx(0.03 ** 2 + 0.05 ** 2)
    assert float(d[0, 0]) == pytest.approx(-0.1)


def test_the_accepted_line_takes_its_lopt_record():
    ctx = _hfs_ctx(0.75, 0.034)
    var, d = lp.candidate_terms(
        ctx, _tt(0.4, 0.01, 0.0, 0.0, 0.0, acc_i=0, acc_wn=30000.25,
                 acc_unc=0.233),
        np.array([[0]]), np.array([[30000.0]]))
    assert float(var[0, 0]) == pytest.approx(0.233 ** 2)
    assert float(d[0, 0]) == pytest.approx(0.25)


def _record(wn, unc, low, upp, weight, flag=''):
    """One record of a LOPT input, field by field as make_LOPT_input writes
    them."""
    rec = [' '] * 93
    for name, text in (('wavenumber', '%.3f' % wn), ('uncertainty',
                                                      '%.3f' % unc),
                       ('intensity', '100.0'), ('lower_level', low),
                       ('upper_level', upp), ('flags', flag),
                       ('weight', '%.4f' % weight), ('units', 'cm-1')):
        a, b = make_LOPT_input.FIELDS[name]
        rec[a - 1:a - 1 + len(text)] = list(text)
    return ''.join(rec).rstrip() + '\n'


def test_a_blend_component_is_placed_through_the_centroid():
    """Two components of one line, weights 0.75 and 0.25.  With the other's
    Ritz value held, the record puts each at (wn - w_o Ritz_o) / w_t, with
    the uncertainty sigma / w_t."""
    ctx = fake_ctx([1000.0])
    ctx.hfs = None
    ctx.e_final = {'A': 0.0, 'B': 0.0, 'U': 1000.2, 'V': 999.6}
    acc = pd.DataFrame({'low_id': ['A', 'B'], 'upp_id': ['U', 'V'],
                        'wn_obs': [1000.0, 1000.0], 'line_idx': [0, 0]})
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, 'LOPT_input_lines.txt')
        with open(path, 'w', encoding='ascii', newline='\n') as fh:
            fh.write(_record(1000.0, 0.02, 'A', 'U', 0.75))
            fh.write(_record(1000.0, 0.02, 'B', 'V', 0.25))
        old = lp.LOPT_INPUT_FILE
        lp.LOPT_INPUT_FILE = path
        try:
            lp.accepted_fit(acc, ctx, log=lambda *a: None)
        finally:
            lp.LOPT_INPUT_FILE = old
    i, nu_a, u_a = ctx.acc_rec[('A', 'U')]
    assert nu_a == pytest.approx((1000.0 - 0.25 * 999.6) / 0.75)
    assert u_a == pytest.approx(0.02 / 0.75)
    _, nu_b, u_b = ctx.acc_rec[('B', 'V')]
    assert nu_b == pytest.approx((1000.0 - 0.75 * 1000.2) / 0.25)
    assert u_b == pytest.approx(0.02 / 0.25)
    # the centroid of the two placed components is the line itself
    assert list(acc['wn_fit']) == [1000.0, 1000.0]


def test_a_flagged_record_is_not_an_accepted_one():
    ctx = fake_ctx([1000.0])
    ctx.hfs = None
    ctx.e_final = {}
    acc = pd.DataFrame({'low_id': ['A'], 'upp_id': ['U'], 'wn_obs': [1000.0],
                        'line_idx': [0]})
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, 'LOPT_input_lines.txt')
        with open(path, 'w', encoding='ascii', newline='\n') as fh:
            fh.write(_record(1000.0, 0.5, 'A', 'U', 0.0, flag='P'))
        old = lp.LOPT_INPUT_FILE
        lp.LOPT_INPUT_FILE = path
        try:
            lp.accepted_fit(acc, ctx, log=lambda *a: None)
        finally:
            lp.LOPT_INPUT_FILE = old
    # the candidate formula: the line's own 0.05 of fake_ctx
    assert ctx.acc_rec[('A', 'U')][2] == pytest.approx(0.05)


# ---------------------------------------------------------------------------
# 2. Soft blend credit
# ---------------------------------------------------------------------------
def test_claims_that_leave_room_are_believed_as_they_stand():
    x = math.log(1.0e4)
    assert float(lp.softened_claim(x, 10.0, 1.0)) == pytest.approx(
        math.log(10.0), abs=1e-3)


def test_claims_far_above_the_feature_are_believed_up_to_its_light():
    x = math.log(1.0e4)
    assert float(lp.softened_claim(x, 1.0e7, 1.0)) == pytest.approx(
        x, abs=1e-3)
    # and where they equal it, one width's worth of doubt
    assert float(lp.softened_claim(x, 1.0e4, 1.0)) == pytest.approx(
        x - math.log(2.0))


def test_a_free_feature_has_no_claims_to_soften():
    assert float(lp.softened_claim(5.0, 0.0, 1.0)) == -np.inf


def _claimed_ln_g(i_obs, claimed, i_pred):
    ctx = fake_ctx([50000.0], i_obs=[i_obs], claimed=[claimed])
    ctx.by_level['X'] = dict(
        partner=np.array(['A']), sign=np.array([1.0]),
        i_pred=np.array([float(i_pred)]), e_m=np.array([0.0]),
        u_m=np.array([0.0]), degenerate=np.array([False]))
    _, tab = lp.ln_ratio(ctx, 'X', np.array([50000.0]), detail=True)
    return float(tab['ln_G'][0])


def _hard_ln_g(i_obs, claimed, i_pred, s=1.2):
    x = math.log(i_obs)
    ln_bg = lp.ln_norm(x, math.log(1.0e4), 1.2)
    return float(lp.ln_intensity(x, claimed + i_pred, s, ln_bg)
                 - lp.ln_intensity(x, claimed, s, ln_bg))


def test_a_feature_already_explained_still_lends_a_claimant_some_light():
    """Claims of three times the feature's light: the old comparison charged
    a claimant of the feature's own size heavily; the soft one charges it
    less, though still something."""
    soft = _claimed_ln_g(1.0e4, 3.0e4, 1.0e4)
    hard = _hard_ln_g(1.0e4, 3.0e4, 1.0e4)
    assert hard < soft < 0.0


def test_a_feature_with_room_is_judged_as_before():
    soft = _claimed_ln_g(1.0e4, 1.0e2, 5.0e3)
    hard = _hard_ln_g(1.0e4, 1.0e2, 5.0e3)
    assert soft == pytest.approx(hard, abs=0.02)


# ---------------------------------------------------------------------------
# 3. The action label
# ---------------------------------------------------------------------------
def _row(gain):
    return {'near_dE': np.nan, 'dE_alt': -0.75, 'n_kept': 8, 'n_own': 8,
            'kept_light': 1.0, 'top_share': 0.2, 'n_free': 0,
            'ln_R_alt': 7.9, 'free_gain': 0.0, 'look': gain, 'gain': gain}


def test_an_alternate_the_evidence_does_not_prefer_carries_no_action():
    assert lp.action_of(_row(-0.13)) == ''
    assert lp.action_of(_row(0.0)) == ''


def test_a_preferred_alternate_is_still_named():
    assert lp.action_of(_row(0.5)) == 'refit'
