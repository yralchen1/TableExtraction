"""Tests for the three things `level_positions` learned to do.

Run from the LineClass directory:  python -m pytest tests -q

The first is a correction.  The audit used to call every alternate position
within two wavenumbers a `refit` - the level stays where it is, and one bad
assignment is dragging the least-squares fit off the position its own lines
want.  That is true only when the level KEEPS the recorded lines it is
assigned now, and distance does not decide it: 059003.000617 keeps none of its
three at a position 1.3 cm^-1 away, so every one of its assignments has to be
dropped, another set made, and the level entered in a ledger at a new
position, which is the opposite piece of work.  The tests below fix the two
numbers that separate the cases, n_own and n_kept, and the rule that reads
them.

The second is a search.  A level of the calculation that has never been found
has no adopted energy and no accepted lines, only a row in enlev.dat and, in
Cowan's transition list, transitions to levels that HAVE been found.  Those
are enough to evaluate ln R at any trial energy, so the window the calculation
allows can be scanned and the positions where ln R > 0 reported.  The tests
build a miniature run - eight recorded lines, three found partners, one
calculated level that nobody has placed - and check that the predictions are
assembled from the right transitions with the pipeline's own intensity
formula, that the scan finds a level put there on purpose, and that a position
resting on one line is called what it is.

The third is a second correction, to the likelihood itself and to the same
disposition.  A match on a free feature has two readings - the feature IS the
transition, or the transition is hidden in a feature that would have been
recorded there in any case - and a row is worth the better of them.  Reading
the brightness under the second and the position under the first, which is
what the intensity floor alone did, let a prediction of 1 unit landing on a
feature of 10000 collect the whole of the positional credit for a
coincidence.  And a refit is allowed to give up half a level's lines but not
the greater part of its light, since the line a short move drops is usually
the brightest one and its Ritz mismatch is the thing to explain rather than
the thing to discard.

None of them reads the pipeline's output; one, marked `real_file`, opens the
real IDEN2 directory.
"""
import math
import os
import sys
import tempfile

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import classify_lines as cl                # noqa: E402
import cowan_gA                            # noqa: E402
import level_interchange as li             # noqa: E402
import level_positions as lp               # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REAL_ENLEV = os.path.join(HERE, 'IDEN2', 'enlev.dat')

real_file = pytest.mark.skipif(not os.path.exists(REAL_ENLEV),
                               reason='the real IDEN2 directory is absent')

C = float(cl.CFG.intensity_model['C'])
KT = float(cl.CFG.intensity_model['kT'])


# ---------------------------------------------------------------------------
# A miniature run
# ---------------------------------------------------------------------------
def fake_ctx(wn_obs, i_obs=None, claimed=None):
    """A Context carrying only what ln_ratio and support ask of a run.

    `calib` is None, which makes `p_obs_vec` return the coverage alone, so
    nothing here depends on the noise threshold fitted to any real run.  The
    recorded lines are quoted to 0.05 cm^-1 in a neighbourhood 0.1 lines per
    wavenumber thick; `claimed` is the predicted intensity that accepted
    assignments already put on each of them, zero everywhere by default.
    """
    ctx = lp.Context()
    wn = np.asarray(wn_obs, dtype=float)
    ctx.wn_o = wn
    ctx.int_o = (np.full(len(wn), 1.0e4) if i_obs is None
                 else np.asarray(i_obs, dtype=float))
    ctx.ln_int_o = np.log(ctx.int_o)
    ctx.unc_o = np.full(len(wn), 0.05)
    ctx.meas_o = np.full(len(wn), 0.05)
    ctx.sig_loc_o = np.full(len(wn), 0.05)
    ctx.rho_o = np.full(len(wn), 0.1)
    ctx.char_o = np.array([''] * len(wn))
    ctx.bg_mu_o = np.full(len(wn), math.log(1.0e4))
    ctx.bg_sd_o = np.full(len(wn), 1.2)
    ctx.claimed_tot = (np.zeros(len(wn)) if claimed is None
                       else np.asarray(claimed, dtype=float))
    ctx.n_acc_line = (ctx.claimed_tot > 0).astype(int)
    ctx.k_n = {1: 1.0, 2: 1.0}
    ctx.w_hfs = {}
    ctx.eta = 0.0
    ctx.s, ctx.s_L = 1.2, 0.2
    ctx.calib = None
    ctx.bias = {}
    ctx.own_claim = {}
    ctx.by_level = {}
    ctx.e_final = {}
    ctx.u_M = {}
    ctx.n_acc_level = {}
    ctx.shared = {}
    # no free/identified distinction in a synthetic run: the background of
    # unrelated lines is the whole list, as it was before the calculated
    # transitions of the unfound levels were brought into it
    lp.plain_background(ctx)
    return ctx


def theory(rows, trans_rows):
    """`(en, trans, mapping, e_meas, id_of)` from a handful of levels.

    `rows` are (idx, E_calc, E_obs, J, cfg, known); `trans_rows` are
    (idx_a, idx_b, gA).  The Cowan level numbering is taken to be the IDEN2
    row numbering, so `mapping` is the identity - the real correspondence is
    `cowan_gA.match_to_enlev`'s business and is tested there.
    """
    en = pd.DataFrame(
        [dict(idx=i, E_calc=ec, E_obs=eo, omc=eo - ec, J=j, cfg=cfg,
              label=f'{cfg} _test', term='test', u_obs=0.01, known=k)
         for i, ec, eo, j, cfg, k in rows]).set_index('idx')
    trans = pd.DataFrame([dict(lid1=a, lid2=b, gA=g)
                          for a, b, g in trans_rows])
    mapping = {int(i): int(i) for i in en.index}
    e_meas = {int(i): float(en['E_obs'][i])
              for i in en.index if en['known'][i]}
    id_of = {int(i): f'L{int(i):03d}' for i in en.index}
    return en, trans, mapping, e_meas, id_of


# ---------------------------------------------------------------------------
# The correction: a refit is a level that keeps its lines
# ---------------------------------------------------------------------------
def test_observed_index_finds_the_line_a_wavenumber_names():
    ctx = fake_ctx([100.0, 200.0, 300.0])
    got = lp.observed_index(ctx, [300.0, 100.0, 199.9])
    assert list(got) == [2, 0, 1]


def test_a_near_alternate_that_keeps_its_lines_is_a_refit():
    row = dict(near_dE=np.nan, dE_alt=0.5, top_share=0.3, n_free=6,
               ln_R_alt=4.0, free_gain=12.0, look=5.0, n_own=7, n_kept=6)
    assert lp.disposition(row) == 'refit'


def test_a_near_alternate_that_keeps_none_of_them_is_not_a_refit():
    """059003.000617: 1.3 cm^-1 away, and not one of its three lines left."""
    row = dict(near_dE=np.nan, dE_alt=1.32, top_share=0.45, n_free=3,
               ln_R_alt=4.93, free_gain=8.04, look=-0.04, n_own=3, n_kept=0)
    assert lp.disposition(row) != 'refit'


def test_keeping_exactly_half_is_not_enough():
    """059003.000602: two of its four lines survive, which is a move."""
    row = dict(near_dE=np.nan, dE_alt=-0.79, top_share=0.32, n_free=3,
               ln_R_alt=4.97, free_gain=6.01, look=-0.99, n_own=4, n_kept=2)
    assert lp.disposition(row) != 'refit'


def test_a_level_with_no_lines_at_all_cannot_be_refitted():
    row = dict(near_dE=np.nan, dE_alt=0.1, top_share=0.3, n_free=5,
               ln_R_alt=4.0, free_gain=12.0, look=5.0, n_own=0, n_kept=0)
    assert lp.disposition(row) != 'refit'


def test_an_interchange_is_still_decided_first():
    """The alternate lands on another level: whose lines they are is not
    this scan's question, however many of them the level would keep."""
    row = dict(near_dE=0.2, dE_alt=0.3, top_share=0.3, n_free=6,
               ln_R_alt=4.0, free_gain=12.0, look=5.0, n_own=7, n_kept=7)
    assert lp.disposition(row) == 'interchange'


def test_support_counts_the_lines_the_level_would_give_up():
    """Two recorded lines, both assigned to the level, and a position two
    wavenumbers away that matches neither of them."""
    ctx = fake_ctx([50000.0, 60000.0])
    ctx.own_claim = {'X': {0: 1.0e4, 1: 1.0e4}}
    ctx.by_level['X'] = dict(
        partner=np.array(['A', 'B']), sign=np.array([1.0, 1.0]),
        i_pred=np.array([1.0e4, 1.0e4]),
        e_m=np.array([60000.0, 50000.0]), u_m=np.array([0.0, 0.0]),
        degenerate=np.array([False, False]))
    here = lp.support(ctx, 'X', 110000.0)
    away = lp.support(ctx, 'X', 110002.0)
    assert here['n_own'] == 2 and here['n_kept'] == 2
    assert away['n_own'] == 2 and away['n_kept'] == 0


# ---------------------------------------------------------------------------
# The search: a level nobody has found
# ---------------------------------------------------------------------------
LEVELS = [(1, 0.0, 0.0, 0.5, 'f25d', True),
          (2, 10000.0, 10000.0, 1.5, 'f25d', True),
          (3, 20000.0, 20000.0, 2.5, 'f25d', True),
          (4, 90000.0, 90000.0, 2.5, 'f25f', False),
          (5, 95000.0, 95000.0, 1.5, 'f25f', False)]
# the unfound level 4 reaches all three found levels, and level 5; the last
# is a transition between two levels neither of which has been placed
TRANS = [(4, 1, 1.0e10), (4, 2, 1.0e10), (4, 3, 1.0e10), (4, 5, 1.0e10)]


def test_only_transitions_to_a_found_level_can_be_searched_on():
    """The other end of the rest is itself unplaced, so they predict no
    wavenumber at all."""
    ctx = fake_ctx([70000.0, 80000.0, 90000.0])
    en, trans, mapping, e_meas, id_of = theory(LEVELS, TRANS)
    for i in (1, 2, 3):
        ctx.e_final[f'L{i:03d}'] = float(en['E_obs'][i])
        ctx.n_acc_level[f'L{i:03d}'] = 5
    n = lp.register_unknown(ctx, 4, en, trans, mapping, e_meas, id_of,
                            log=lambda *a: None)
    assert n == 3
    assert sorted(ctx.by_level[lp.UNKNOWN_ID]['partner']) == ['L001', 'L002',
                                                              'L003']


def test_the_predicted_intensity_is_the_pipelines_own_formula():
    ctx = fake_ctx([70000.0, 80000.0, 90000.0])
    en, trans, mapping, e_meas, id_of = theory(LEVELS, TRANS)
    for i in (1, 2, 3):
        ctx.e_final[f'L{i:03d}'] = float(en['E_obs'][i])
        ctx.n_acc_level[f'L{i:03d}'] = 5
    lp.register_unknown(ctx, 4, en, trans, mapping, e_meas, id_of,
                        log=lambda *a: None)
    g = ctx.by_level[lp.UNKNOWN_ID]
    for p, e_m, i_pred in zip(g['partner'], g['e_m'], g['i_pred']):
        nu = 90000.0 - e_m
        assert i_pred == pytest.approx(
            C * 1.0e10 * (nu / 1.0e8) * math.exp(-90000.0 / KT))
    # the unfound level is the upper partner of all three
    assert list(g['sign']) == [1.0, 1.0, 1.0]


def test_the_unfound_level_is_registered_at_its_calculated_energy():
    ctx = fake_ctx([70000.0, 80000.0, 90000.0])
    en, trans, mapping, e_meas, id_of = theory(LEVELS, TRANS)
    for i in (1, 2, 3):
        ctx.e_final[f'L{i:03d}'] = float(en['E_obs'][i])
        ctx.n_acc_level[f'L{i:03d}'] = 5
    lp.register_unknown(ctx, 4, en, trans, mapping, e_meas, id_of,
                        log=lambda *a: None)
    assert ctx.e_final[lp.UNKNOWN_ID] == 90000.0


def test_the_window_is_three_configuration_widths_either_side():
    ctx = fake_ctx([70000.0, 80000.0, 90000.0])
    en, trans, mapping, e_meas, id_of = theory(LEVELS, TRANS)
    for i in (1, 2, 3):
        ctx.e_final[f'L{i:03d}'] = float(en['E_obs'][i])
        ctx.n_acc_level[f'L{i:03d}'] = 5
    lp.register_unknown(ctx, 4, en, trans, mapping, e_meas, id_of,
                        log=lambda *a: None)
    r = lp.scan_unknown(ctx, 4, en, {'f25f': 20.0}, 132.0, step=0.5)
    assert r['lo'] == pytest.approx(90000.0 - 60.0)
    assert r['hi'] == pytest.approx(90000.0 + 60.0)
    assert r['W'] == 20.0


def test_a_configuration_with_no_found_level_takes_the_list_wide_width():
    ctx = fake_ctx([70000.0, 80000.0, 90000.0])
    en, trans, mapping, e_meas, id_of = theory(LEVELS, TRANS)
    for i in (1, 2, 3):
        ctx.e_final[f'L{i:03d}'] = float(en['E_obs'][i])
        ctx.n_acc_level[f'L{i:03d}'] = 5
    lp.register_unknown(ctx, 4, en, trans, mapping, e_meas, id_of,
                        log=lambda *a: None)
    r = lp.scan_unknown(ctx, 4, en, {}, 132.0, step=2.0)
    assert r['W'] == 132.0


def test_the_scan_finds_a_level_put_there_on_purpose():
    """Three recorded lines placed at the wavenumbers a level 30 cm^-1 below
    the calculated position would predict, with the intensities it predicts
    for them, and nowhere else."""
    e_true = 89970.0
    wn = sorted(e_true - np.array([0.0, 10000.0, 20000.0]))
    i_obs = [C * 1.0e10 * (nu / 1.0e8) * math.exp(-90000.0 / KT)
             for nu in (70000.0, 80000.0, 90000.0)]
    ctx = fake_ctx(wn, i_obs=i_obs)
    en, trans, mapping, e_meas, id_of = theory(LEVELS, TRANS)
    for i in (1, 2, 3):
        ctx.e_final[f'L{i:03d}'] = float(en['E_obs'][i])
        ctx.n_acc_level[f'L{i:03d}'] = 5
    lp.register_unknown(ctx, 4, en, trans, mapping, e_meas, id_of,
                        log=lambda *a: None)
    r = lp.scan_unknown(ctx, 4, en, {'f25f': 20.0}, 132.0, step=0.02)
    assert r['positions'], 'the scan found no position at all'
    assert r['positions'][0][0] == pytest.approx(e_true, abs=0.05)
    assert r['positions'][0][1] > 0.0


def test_a_position_resting_on_one_line_is_called_no_support():
    assert lp.position_verdict(1, 4.0, 1.0, 3.0) == 'no support'
    assert lp.position_verdict(5, 12.0, 0.9, 5.0) == 'no support'


def test_a_position_the_audit_would_call_firm_is_called_firm():
    assert lp.position_verdict(lp.AUDIT_MIN_FREE, lp.AUDIT_MIN_GAIN,
                               lp.AUDIT_MAX_SHARE, lp.AUDIT_MIN_LOOK) == 'firm'


def test_a_position_that_does_not_survive_the_look_elsewhere_is_weak():
    """The same support, found by a scan with more chances to find it."""
    assert lp.position_verdict(6, 15.0, 0.2, lp.AUDIT_MIN_LOOK - 0.1) == 'weak'


def test_the_look_elsewhere_correction_pays_for_the_whole_scan():
    """Two positions, so each of them is worth ln 2 less than its ln R."""
    ctx = fake_ctx([70000.0, 80000.0])
    r = dict(E_calc=90000.0, W=20.0,
             positions=[(89990.0, 6.0), (90010.0, 5.0)])
    ctx.by_level[lp.UNKNOWN_ID] = dict(
        partner=np.array(['A']), sign=np.array([1.0]),
        i_pred=np.array([1.0e4]), e_m=np.array([70000.0]),
        u_m=np.array([0.0]), degenerate=np.array([False]))
    tab = lp.unknown_table(ctx, r)
    assert list(tab['look']) == pytest.approx([6.0 - math.log(2),
                                               5.0 - math.log(2)])


@real_file
def test_the_row_numbers_are_the_ones_iden2_shows():
    """--unknown takes the running number in column 1 of enlev.dat, and every
    unfound row of the real file can be asked for by it."""
    en = li.read_enlev(REAL_ENLEV).set_index('idx')
    unfound = en[~en['known']]
    assert len(unfound) > 100
    for i in list(unfound.index[:5]) + list(unfound.index[-5:]):
        assert int(i) in en.index
        assert np.isfinite(float(en['E_calc'][i]))


# ---------------------------------------------------------------------------
# The second correction: what a strong line on a weak prediction is worth,
# and which line a refit is allowed to give up
# ---------------------------------------------------------------------------
def test_a_feature_the_prediction_cannot_account_for_is_marked():
    """10000 units of light where 1 was predicted: an unidentified line
    explains that better than the prediction does."""
    x = np.log([1.0e4])
    ln_bg = lp.ln_norm(x, math.log(1.0e4), 1.2)
    assert bool(lp.intensity_floored(x, [1.0], 1.2, ln_bg)[0])


def test_a_feature_as_bright_as_predicted_is_not_marked():
    x = np.log([1.0e4])
    ln_bg = lp.ln_norm(x, math.log(1.0e4), 1.2)
    assert not bool(lp.intensity_floored(x, [1.0e4], 1.2, ln_bg)[0])


def test_a_feature_dimmer_than_predicted_is_not_marked():
    """The floor is one-sided: a deficit is evidence, and is charged."""
    x = np.log([1.0e2])
    ln_bg = lp.ln_norm(x, math.log(1.0e4), 1.2)
    assert not bool(lp.intensity_floored(x, [1.0e4], 1.2, ln_bg)[0])


def _one_prediction(ctx, i_pred):
    ctx.by_level['X'] = dict(
        partner=np.array(['A']), sign=np.array([1.0]),
        i_pred=np.array([float(i_pred)]), e_m=np.array([0.0]),
        u_m=np.array([0.0]), degenerate=np.array([False]))
    return float(lp.ln_ratio(ctx, 'X', np.array([50000.0]))[0][0])


def test_a_strong_line_credits_a_weak_prediction_with_almost_nothing():
    """A recorded feature of 10000 exactly where a prediction of 1 falls.

    The line would have been recorded there whether or not the level is at
    this energy, so the coincidence is worth almost nothing - at most what a
    line known to be somewhere in the window is worth, not what one drawn
    from the local line density is.  This is the row that carried four levels
    to alternate positions assigning a very strong line to a very weak
    transition.
    """
    ctx = fake_ctx([50000.0], i_obs=[1.0e4])
    coincidence = _one_prediction(ctx, 1.0)
    genuine = _one_prediction(ctx, 1.0e4)
    assert coincidence < 1.7
    assert genuine > coincidence + 2.0


def test_a_line_merely_brighter_than_predicted_keeps_its_position():
    """The correction is not a rule against under-predicted lines.

    In a neighbourhood thin enough for a coincidence to be unlikely, a line
    thirty times brighter than the calculation says is still much better
    explained as the transition than as a line that was going to be there
    anyway, and the first reading goes on winning for it.  One thirty
    thousand times brighter is not.
    """
    ctx = fake_ctx([50000.0], i_obs=[3.0e4])
    ctx.rho_o = np.full(1, 0.01)
    lp.plain_background(ctx)
    ordinary = _one_prediction(ctx, 1.0e3)
    assert lp.support(ctx, 'X', 50000.0)['n_free'] == 1
    wild = _one_prediction(ctx, 1.0)
    assert lp.support(ctx, 'X', 50000.0)['n_free'] == 0
    assert ordinary > wild + 1.0


def test_such_a_row_is_not_free_support():
    ctx = fake_ctx([50000.0], i_obs=[1.0e4])
    _one_prediction(ctx, 1.0)
    assert lp.support(ctx, 'X', 50000.0)['n_free'] == 0
    _one_prediction(ctx, 1.0e4)
    assert lp.support(ctx, 'X', 50000.0)['n_free'] == 1


def test_the_light_of_a_blended_line_is_split_by_branching_fraction():
    """A feature of 1000 units carrying two accepted components, this
    level's worth a quarter of the prediction: 250 units are its."""
    ctx = fake_ctx([50000.0], i_obs=[1000.0], claimed=[400.0])
    assert lp.own_light(ctx, {0: 100.0}) == pytest.approx({0: 250.0})


def test_a_line_nothing_else_claims_brings_the_level_all_of_it():
    ctx = fake_ctx([50000.0], i_obs=[1000.0])
    assert lp.own_light(ctx, {0: 100.0}) == pytest.approx({0: 1000.0})


def test_support_weighs_the_lines_a_move_gives_up_as_well_as_counting_them():
    """Two lines, one of them a hundred times the brighter, and a position
    that matches only the faint one: six lines out of seven would say the
    level keeps what matters, and the light says it does not."""
    ctx = fake_ctx([50000.0, 60000.0], i_obs=[1.0e2, 1.0e4])
    ctx.own_claim = {'X': {0: 1.0e2, 1: 1.0e4}}
    ctx.by_level['X'] = dict(
        partner=np.array(['A', 'B']), sign=np.array([1.0, 1.0]),
        i_pred=np.array([1.0e2, 1.0e4]),
        e_m=np.array([60000.0, 50000.0]), u_m=np.array([0.0, 0.0]),
        degenerate=np.array([False, False]))
    both = lp.support(ctx, 'X', 110000.0)
    assert both['n_kept'] == 2 and both['kept_light'] == pytest.approx(1.0)
    # the same level with only the faint line within reach: half the lines
    # by count, a hundredth of the light
    ctx.by_level['X'] = dict(
        partner=np.array(['A']), sign=np.array([1.0]),
        i_pred=np.array([1.0e2]), e_m=np.array([60000.0]),
        u_m=np.array([0.0]), degenerate=np.array([False]))
    faint = lp.support(ctx, 'X', 110000.0)
    assert faint['n_own'] == 2 and faint['n_kept'] == 1
    assert faint['kept_light'] == pytest.approx(1.0e2 / 1.01e4)


def test_a_refit_that_drops_the_levels_brightest_line_is_not_a_refit():
    """059003.000457: six of its seven lines survive the move, but the one
    it drops is its strongest by an order of magnitude."""
    row = dict(near_dE=np.nan, dE_alt=-0.51, top_share=0.39, n_free=2,
               ln_R_alt=0.85, free_gain=4.59, look=-1.59, n_own=7, n_kept=6,
               kept_light=0.155)
    assert lp.disposition(row) == 'top line'


def test_a_refit_that_keeps_the_light_is_still_a_refit():
    row = dict(near_dE=np.nan, dE_alt=-0.54, top_share=0.30, n_free=3,
               ln_R_alt=5.81, free_gain=8.53, look=0.27, n_own=4, n_kept=4,
               kept_light=1.0)
    assert lp.disposition(row) == 'refit'


# ---------------------------------------------------------------------------
# The background of lines that belong to nothing yet known (section 4a)
# ---------------------------------------------------------------------------
def test_a_free_feature_is_scored_against_the_unrelated_lines_alone():
    """rho_t is the rate of lines the level list does not account for.

    Of the recorded lines, two thirds already carry an accepted transition:
    they are explained without the level under test and cannot be what a
    coincidence is drawn from.  Scoring a free feature against all of them
    understates the row by exactly the ratio of the two rates.
    """
    ctx = fake_ctx([50000.0], i_obs=[1.0e4])
    whole_list = _one_prediction(ctx, 1.0e4)
    ctx.rho_bg_o = ctx.rho_o / 4.0
    unrelated_only = _one_prediction(ctx, 1.0e4)
    assert unrelated_only == pytest.approx(whole_list + math.log(4.0), abs=0.02)


def test_the_unrelated_rate_does_not_touch_a_claimed_feature():
    """A blended row is already scored against 1/(2W), not against rho.

    The question there is whether adding the prediction improves an account of
    a brightness that other transitions already give, and no density of
    unrelated lines enters it.
    """
    ctx = fake_ctx([50000.0], i_obs=[1.0e4], claimed=[1.0e4])
    before = _one_prediction(ctx, 1.0e4)
    ctx.rho_bg_o = ctx.rho_o / 100.0
    assert _one_prediction(ctx, 1.0e4) == pytest.approx(before, abs=1e-9)


def test_a_brightness_an_unfound_level_could_produce_earns_less():
    """g for a free feature includes what the calculation predicts.

    Where a level nobody has found is expected to put a line of just this
    brightness, a feature of that brightness is no longer surprising under H0,
    and the row must be worth less than where nothing is expected at all.
    """
    ctx = fake_ctx([50000.0], i_obs=[1.0e4])
    ctx.bg_mu_f = np.full(1, math.log(10.0))       # free lines are faint here
    ctx.bg_sd_f = np.full(1, 0.5)
    ctx.unk_mu_o = np.full(1, math.log(10.0))
    ctx.unk_sd_o = np.full(1, 0.5)
    nothing_expected = _one_prediction(ctx, 1.0e4)
    ctx.unk_mu_o = np.full(1, math.log(1.0e4))     # an unfound level would
    ctx.unk_sd_o = np.full(1, 0.5)                 # put a line of 10000 here
    expected = _one_prediction(ctx, 1.0e4)
    assert expected < nothing_expected - 1.0


def test_the_free_background_is_capped_by_the_whole_list_and_floored_below_it():
    """Neither estimate may run away where there is no free line near.

    The free lines are a subset of the recorded ones, so their rate cannot
    exceed the all-lines rate; and where both estimates extrapolate to nearly
    nothing the rate is held at RHO_BG_FLOOR of it, which caps what the
    extrapolation can buy at ln(1/RHO_BG_FLOOR) per matched row.
    """
    class A:
        plain_background = False

    wn = np.linspace(20000.0, 60000.0, 400)
    ctx = fake_ctx(wn)
    ctx.n_acc_line = np.ones(len(wn), dtype=int)   # every line identified:
    ctx.n_acc_line[:5] = 0                         # too few free ones to fit
    lp.free_line_background(ctx, A(), log=lambda *a, **k: None)
    assert np.all(ctx.rho_bg_o <= ctx.rho_o + 1e-12)
    assert np.all(ctx.rho_bg_o >= lp.RHO_BG_FLOOR * ctx.rho_o - 1e-12)


# ---------------------------------------------------------------------------
# The width model of section 4
# ---------------------------------------------------------------------------
def test_the_character_column_is_reduced_to_the_published_vocabulary():
    """'**' is multiplicity, not width, and '*v'/'*r' are hyperfine."""
    assert lp.normalize_char('**') == ''
    assert lp.normalize_char('w**') == 'w'
    assert lp.normalize_char('cl**') == 'cl'
    assert lp.normalize_char('*v') == 'hfs'
    assert lp.normalize_char('*r') == 'hfs'
    assert lp.normalize_char(' cl ') == 'cl'
    assert lp.normalize_char(None) == ''
    # a code that carries a width keeps it; 'hfs' and '' carry none
    assert lp.CHAR_DLAM.get(('1969', 'cl'), 0.0) > 0
    assert lp.CHAR_DLAM.get(('1974', 'hfs'), 0.0) == 0.0


def test_a_character_width_grows_as_the_square_of_the_wavenumber():
    """The blend and resolution codes are constant in WAVELENGTH.

    A width of dlam angstrom is dlam*1e-8*nu^2 cm^-1, so doubling the
    wavenumber quadruples it - while the precision floor does not move.  Both
    must be visible in sigma_meas.
    """
    lo = float(lp.meas_sigma([50000.0], ['cl'])[0])     # 'cl' is a 1969 code,
    hi = float(lp.meas_sigma([100000.0], ['cl'])[0])    # so both are of 1969
    plain_lo = float(lp.meas_sigma([50000.0], [''])[0])
    plain_hi = float(lp.meas_sigma([100000.0], [''])[0])
    excess_lo = math.sqrt(lo ** 2 - plain_lo ** 2)
    excess_hi = math.sqrt(hi ** 2 - plain_hi ** 2)
    assert abs(excess_hi / excess_lo - 4.0) < 1e-6
    # and the floor is there: at low wavenumber a plain line is all floor
    assert float(lp.meas_sigma([9000.0], [''])[0]) >= lp.UNC_FLOOR


def test_the_wide_code_is_constant_in_wavenumber():
    """'w' is a width of the LINE, not of the reading of the plate.

    Its excess is the same number of cm^-1 wherever the line falls, so it does
    NOT quadruple when the wavenumber doubles the way a wavelength-constant
    code does.  It is also the only code of its era in CHAR_DWN, so no other
    one may pick the term up.
    """
    out = []
    for wn in (15000.0, 30000.0):
        tot = float(lp.meas_sigma([wn], ['w'])[0])
        plain = float(lp.meas_sigma([wn], [''])[0])
        out.append(math.sqrt(tot ** 2 - plain ** 2))
    assert abs(out[1] - out[0]) < 1e-9
    assert abs(out[0] - lp.CHAR_DWN[('1974', 'w')]) < 1e-9
    assert ('1974', 'w') not in lp.CHAR_DLAM
    assert set(lp.CHAR_DWN) == {('1974', 'w')}


def test_a_line_widened_by_hand_keeps_its_quoted_uncertainty():
    """Within a class the quoted value is the rule times one fixed factor.

    So a line quoted far above the usual value of its own class was widened
    deliberately, and the quoted value is honoured; one quoted the usual
    multiple of its class is left to the model, which knows the class.
    """
    wn = np.full(40, 30000.0)
    char = np.array(['w'] * 40)
    rule = float(lp.meas_sigma(wn[:1], [''])[0])
    unc = np.full(40, 2.3 * rule)      # the usual factor of the class
    unc[0] = 0.9                        # one line widened by hand
    s = lp.meas_sigma(wn, char, unc)
    assert s[0] == 0.9
    assert s[1] < unc[1]                # the class factor buys nothing extra
    assert abs(s[1] - float(lp.meas_sigma(wn[:1], ['w'])[0])) < 1e-12


def test_a_hyperfine_width_belongs_to_both_levels_of_the_line():
    w = {'A': 0.2, 'B': 0.1}
    v = lp.hfs_pair(w, ['A', 'A', 'C'], ['B', 'C', 'C'])
    assert abs(v[0] - (0.04 + 0.01)) < 1e-12
    assert abs(v[1] - 0.04) < 1e-12
    assert v[2] == 0.0


# ---------------------------------------------------------------------------
# The hyperfine width of a level the fit cannot measure: section 4.
#
# A level with fewer than HFS_MIN_LINES lines of 1974 used to be given a width
# of zero, which says it has no hyperfine structure.  What it really means is
# that nobody has measured its hyperfine structure, and the best estimate of it
# is what the other levels of the same CONFIGURATION show - the configuration
# being what decides how strongly the outer electron feels the nuclear magnetic
# moment.  Where the configuration itself has too few measured levels the
# outermost electron answers for it instead, its own orbital first and then its
# orbital letter.
# ---------------------------------------------------------------------------
def test_the_outermost_orbital_is_read_off_the_iden2_label():
    # a digit before a letter is a principal quantum number, one after it is
    # an occupation, and the last shell written is the outermost
    assert lp.outer_orbital('f26s') == '6s'        # 4f^2.6s
    assert lp.outer_orbital('f5d2') == '5d'        # 4f.5d^2, the 2 is a count
    assert lp.outer_orbital('4f3') == '4f'         # 4f^3
    assert lp.outer_orbital('fd6p') == '6p'        # 4f.5d.6p
    assert lp.outer_orbital('f5d6s') == '6s'
    assert lp.outer_orbital('') == ''


def test_a_configuration_width_is_the_rms_of_the_widths_it_can_measure():
    # three measured levels of one configuration, 0.1, 0.0 and 0.0: the
    # quantity that adds in quadrature is the mean VARIANCE, not the mean width
    w = [0.10, 0.00, 0.00, 0.30]
    n = [10, 10, 10, 1]                 # the fourth is not measured
    cfg = ['f26s'] * 4
    widths, wide = lp.configuration_widths(w, n, cfg)
    assert abs(widths['f26s'] - math.sqrt(0.01 / 3)) < 1e-12
    # and the unmeasured level takes it, its own 0.30 being unsupported
    assert abs(lp.width_for('f26s', widths, wide) - widths['f26s']) < 1e-12


def test_a_width_below_the_gate_counts_as_the_zero_it_is():
    # HFS_APPLY is noise of the fit, so it must not be averaged in as signal
    w = [lp.HFS_APPLY * 0.9, 0.0, 0.0]
    widths, _ = lp.configuration_widths(w, [10, 10, 10], ['f25g'] * 3)
    assert widths['f25g'] == 0.0


def test_a_configuration_with_too_few_levels_asks_the_outer_electron():
    # 4f.5d.6s has one measured level and cannot speak; its outer electron is
    # the same 6s as 4f^2.6s's, and that is the width it gets
    n_ok = [10] * lp.HFS_CFG_MIN
    w = [0.10] * lp.HFS_CFG_MIN + [0.0]
    cfg = ['f26s'] * lp.HFS_CFG_MIN + ['f5d6s']
    widths, wide = lp.configuration_widths(w, n_ok + [1], cfg)
    assert 'f5d6s' not in widths
    assert abs(lp.width_for('f5d6s', widths, wide) - widths['6s']) < 1e-12
    assert abs(widths['6s'] - 0.10) < 1e-12


def test_failing_that_the_orbital_letter_answers_and_then_the_whole_list():
    # 4f^2.6f has no 6f anywhere else, but an f electron is an f electron
    n_ok = [10] * lp.HFS_CFG_MIN
    w = [0.02001] * lp.HFS_CFG_MIN
    widths, wide = lp.configuration_widths(
        w + [0.0], n_ok + [1], ['f25f'] * lp.HFS_CFG_MIN + ['f26f'])
    assert abs(lp.width_for('f26f', widths, wide) - widths['f']) < 1e-12
    # a level with no configuration at all has nothing to ask, and takes the
    # width of the whole list
    assert lp.width_for('', widths, wide) == wide


def test_a_level_the_fit_can_measure_is_not_touched_by_any_of_this():
    # the whole rule applies to under-determined levels only
    widths, wide = lp.configuration_widths(
        [0.10] * 4, [10] * 4, ['f26s'] * 4)
    assert widths['f26s'] == 0.10          # and the level keeps its own 0.10


def test_the_applied_width_is_what_both_modules_read():
    # level_positions writes w_applied; make_LOPT_input must read the same
    # column, or the uncertainties LOPT is given and the ones the scan uses
    # would disagree
    import make_LOPT_input as mli
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, 'w.csv')
        with open(p, 'w', newline='\n') as fh:
            fh.write('level_id,w_hfs,n_1974,w_applied,w_source\n')
            fh.write('A,0.0,1,0.0300,configuration\n')
            fh.write('B,0.0500,9,0.0500,level\n')
            fh.write('C,0.0,9,0.0,\n')
        got = mli.read_hfs_widths(p)
    assert got == {'A': 0.03, 'B': 0.05}


def test_an_old_width_file_without_the_column_still_works():
    import make_LOPT_input as mli
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, 'w.csv')
        with open(p, 'w', newline='\n') as fh:
            fh.write('level_id,w_hfs,n_1974\n')
            fh.write('A,0.0500,9\n')
            fh.write('B,0.0500,2\n')          # too few lines, the old rule
            fh.write('C,0.0100,9\n')          # below the old gate
        got = mli.read_hfs_widths(p)
    assert got == {'A': 0.05}


# The configuration of a level is looked up through IDEN2/IDEN_level_ids.txt,
# the table giving every level identifier the number of its row in enlev.dat.
# That table is what survives a level being moved in IDEN2: the energy in
# enlev.dat changes, the row number does not.

ENLEV_ROWS = (
    "    1  10000.000   5000.000   10000.500 *        0.500  3.5 / f26s _3H4G/\n"
    "    2  20000.000   5000.000   20000.000 *        0.000  2.5 / f25f _1G2F/\n"
)


def _iden2(tmp):
    """A miniature IDEN2 pair: two levels, the second one moved by 3 cm^-1."""
    en = os.path.join(tmp, 'enlev.dat')
    ids = os.path.join(tmp, 'IDEN_level_ids.txt')
    with open(en, 'w', encoding='latin-1', newline='\n') as fh:
        fh.write(ENLEV_ROWS)
    with open(ids, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('level_id\tIDEN_id\nL1\t1\nL2\t2\n')
    return en, ids


def test_a_level_moved_in_iden2_is_still_found_by_its_identifier():
    with tempfile.TemporaryDirectory() as tmp:
        en, ids = _iden2(tmp)
        # L2 sits 3 cm^-1 from the energy enlev.dat holds - far outside the
        # window that matching on energy allows
        lv = pd.DataFrame({'level_id': ['L1', 'L2'],
                           'E_input': [10000.5, 20003.0]})
        got = lp.read_level_configs(lv, enlev=en, ids=ids)
    assert got == {'L1': 'f26s', 'L2': 'f25f'}


def test_a_level_with_no_energy_yet_takes_no_configuration_by_accident():
    with tempfile.TemporaryDirectory() as tmp:
        en, ids = _iden2(tmp)
        lv = pd.DataFrame({'level_id': ['L2'], 'E_input': [float('nan')]})
        got = lp.read_level_configs(lv, enlev=en, ids=ids)
    # the table answers; nearest-energy matching would have handed it the
    # first row of the file
    assert got == {'L2': 'f25f'}


def test_a_level_the_table_does_not_list_is_reported_not_hidden():
    said = []
    with tempfile.TemporaryDirectory() as tmp:
        en, ids = _iden2(tmp)
        lv = pd.DataFrame({'level_id': ['L9'], 'E_input': [99999.0]})
        got = lp.read_level_configs(lv, log=said.append, enlev=en, ids=ids)
    assert got == {}
    assert any('wants checking' in s for s in said)


# ---------------------------------------------------------------------------
# An interchange is an exchange of identities, so only a level of the same J
# and parity can be the other half of one
# ---------------------------------------------------------------------------
def _kinds():
    per = pd.DataFrame({
        'level_id': ['L1', 'L2', 'L3', 'L4'],
        'E_final': [10000.0, 10000.2, 10040.0, 10000.1],
        'J': ['5/2', '17/2', '5/2', '5/2'],
        'parity': ['e', 'e', 'e', 'o']})
    return lp.level_kinds(per)


def test_a_level_of_another_j_is_not_the_other_half_of_an_interchange():
    """L2 lies 0.2 cm^-1 from the alternate and L3 forty; L3 is the answer,
    because L2 is J = 17/2 and has none of L1's lines to exchange."""
    kinds, kind_of = _kinds()
    nid, dE = lp.nearest_same_kind(kinds, kind_of, 'L1', 10040.0)
    assert nid == 'L3'
    assert dE == pytest.approx(0.0)


def test_a_level_of_the_other_parity_is_not_either():
    kinds, kind_of = _kinds()
    nid, _ = lp.nearest_same_kind(kinds, kind_of, 'L1', 10000.1)
    assert nid == 'L3'


def test_the_level_itself_is_never_its_own_near_level():
    kinds, kind_of = _kinds()
    nid, _ = lp.nearest_same_kind(kinds, kind_of, 'L1', 10000.05)
    assert nid == 'L3'


def test_no_other_level_of_that_j_names_none_at_all():
    kinds, kind_of = _kinds()
    nid, dE = lp.nearest_same_kind(kinds, kind_of, 'L2', 10000.0)
    assert nid == ''
    assert math.isnan(dE)


def test_an_interchange_needs_the_near_level_of_the_same_kind():
    """disposition() reads near_dE, and near_dE is now only ever filled in
    from a level that could actually be swapped with this one."""
    kinds, kind_of = _kinds()
    _, dE = lp.nearest_same_kind(kinds, kind_of, 'L1', 10000.2)
    assert abs(dE) > lp.INTERCHANGE_DE


# ---------------------------------------------------------------------------
# --unknown on a row the run has already found
# ---------------------------------------------------------------------------
def _found_ctx():
    """The miniature run with level 4 already found as L004, holding one
    accepted line - the 90000 cm^-1 feature, which is its transition to L001
    and which nothing else claims.  The recorded intensities are those the
    level's own predictions ask for, so that a match is read as the line
    rather than as a transition hidden in a much brighter feature."""
    ctx = fake_ctx([70000.0, 80000.0, 90000.0],
                   i_obs=[1.2e6, 1.38e6, 1.55e6],
                   claimed=[0.0, 0.0, 1.55e6])
    en, trans, mapping, e_meas, id_of = theory(LEVELS, TRANS)
    for i in (1, 2, 3):
        ctx.e_final[f'L{i:03d}'] = float(en['E_obs'][i])
        ctx.n_acc_level[f'L{i:03d}'] = 5
    ctx.e_final['L004'] = 90000.0
    ctx.own_claim['L004'] = {2: 1.55e6}
    ctx.n_acc_level['L004'] = 1
    return ctx, en, trans, mapping, e_meas, id_of


def test_a_found_rows_own_lines_are_released_to_the_search():
    ctx, en, trans, mapping, e_meas, id_of = _found_ctx()
    lp.register_unknown(ctx, 4, en, trans, mapping, e_meas, id_of,
                        log=lambda *a: None)
    # the search reads the level's own claims as its own, exactly as the scan
    # of a known level does, so its own lines come out free rather than as a
    # blend of the level with itself
    assert ctx.own_claim[lp.UNKNOWN_ID] == {2: 1.55e6}
    assert ctx.unknown_level['level_id'] == 'L004'
    assert ctx.unknown_level['in_run'] is True
    assert ctx.unknown_level['n_own'] == 1


def test_the_released_lines_make_the_position_the_level_holds_score():
    ctx, en, trans, mapping, e_meas, id_of = _found_ctx()
    lp.register_unknown(ctx, 4, en, trans, mapping, e_meas, id_of,
                        log=lambda *a: None)
    with_release = lp.support(ctx, lp.UNKNOWN_ID, 90000.0)
    ctx.own_claim.pop(lp.UNKNOWN_ID)          # as it was before the fix
    without = lp.support(ctx, lp.UNKNOWN_ID, 90000.0)
    assert with_release['n_free'] > without['n_free']
    assert with_release['ln_R'] > without['ln_R']
    # and the cost of a move is reported: the level has one line of its own
    assert with_release['n_own'] == 1
    assert with_release['n_kept'] == 1


def test_a_row_found_in_enlev_but_absent_from_the_run_releases_nothing():
    ctx = fake_ctx([70000.0, 80000.0, 90000.0])
    en, trans, mapping, e_meas, id_of = theory(LEVELS, TRANS)
    for i in (1, 2, 3):
        ctx.e_final[f'L{i:03d}'] = float(en['E_obs'][i])
        ctx.n_acc_level[f'L{i:03d}'] = 5
    lp.register_unknown(ctx, 4, en, trans, mapping, e_meas, id_of,
                        log=lambda *a: None)
    assert lp.UNKNOWN_ID not in ctx.own_claim
    assert ctx.unknown_level['in_run'] is False


def test_one_rows_release_does_not_carry_over_to_the_next():
    ctx, en, trans, mapping, e_meas, id_of = _found_ctx()
    lp.register_unknown(ctx, 4, en, trans, mapping, e_meas, id_of,
                        log=lambda *a: None)
    del ctx.e_final['L004']                   # row 5 is not a level of the run
    lp.register_unknown(ctx, 5, en, trans, mapping, e_meas, id_of,
                        log=lambda *a: None)
    assert lp.UNKNOWN_ID not in ctx.own_claim


# ---------------------------------------------------------------------------
# --drop-all-questionable
# ---------------------------------------------------------------------------
class _Args:
    plain_background = True


def _claimed_ctx():
    """Three recorded lines, two of them claimed: one by a questionable level
    Q alone, one by Q sharing the feature with K."""
    ctx = fake_ctx([70000.0, 80000.0, 90000.0], claimed=[0.0, 3.0, 7.0])
    ctx.acc = pd.DataFrame([
        dict(low_id='K', upp_id='Q', line_idx=1, calc_intens=3.0),
        dict(low_id='K', upp_id='Q', line_idx=2, calc_intens=4.0),
        dict(low_id='K', upp_id='M', line_idx=2, calc_intens=3.0)])
    ctx.own_claim = {'Q': {1: 3.0, 2: 4.0}, 'K': {1: 3.0, 2: 7.0},
                     'M': {2: 3.0}}
    ctx.n_acc_line = np.array([0, 1, 2])
    ctx.n_acc_level = {'Q': 2, 'K': 2, 'M': 1}
    return ctx


def test_releasing_a_level_frees_the_feature_it_held_alone():
    ctx = _claimed_ctx()
    n = lp.release_levels(ctx, ['Q'], _Args(), log=lambda *a: None)
    assert n == 2
    assert ctx.claimed_tot[1] == pytest.approx(0.0)
    assert ctx.n_acc_line[1] == 0
    assert 'Q' not in ctx.own_claim


def test_releasing_a_level_leaves_the_other_claims_on_a_shared_feature():
    ctx = _claimed_ctx()
    lp.release_levels(ctx, ['Q'], _Args(), log=lambda *a: None)
    # the 4.0 of the Q-K row goes; the 3.0 of the K-M row stays, and so does
    # the feature's second claimant
    assert ctx.claimed_tot[2] == pytest.approx(3.0)
    assert ctx.n_acc_line[2] == 1
    assert ctx.own_claim['M'] == {2: 3.0}
    assert ctx.own_claim['K'] == {2: pytest.approx(3.0)}


def test_a_released_level_still_serves_as_a_partner():
    """Its energy is what it was, and n_acc_level - which is what marks a
    partner degenerate - is untouched, or the release would strike out the
    very predictions that reach the lines it has just freed."""
    ctx = _claimed_ctx()
    lp.release_levels(ctx, ['Q'], _Args(), log=lambda *a: None)
    assert ctx.n_acc_level['Q'] == 2


def test_releasing_levels_that_hold_nothing_changes_nothing():
    ctx = _claimed_ctx()
    before = ctx.claimed_tot.copy()
    n = lp.release_levels(ctx, ['Z'], _Args(), log=lambda *a: None)
    assert n == 0
    assert np.array_equal(ctx.claimed_tot, before)


def test_the_questionable_levels_are_the_ones_the_report_marks():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, 'level_positions.csv')
        pd.DataFrame({'level_id': ['L1', 'L2', 'L3'],
                      'question': ['?', '', '?']}).to_csv(path, index=False)
        assert lp.questionable_levels(path) == ['L1', 'L3']


def test_a_report_without_the_column_is_refused_not_guessed_at():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, 'level_positions.csv')
        pd.DataFrame({'level_id': ['L1']}).to_csv(path, index=False)
        with pytest.raises(SystemExit):
            lp.questionable_levels(path)


def test_the_option_is_refused_without_unknown():
    with pytest.raises(SystemExit):
        lp.main(['--drop-all-questionable'])


# ---------------------------------------------------------------------------
# What a match is: the hidden reading needs the light, n_match needs the
# brightness, and the window can be widened for a mixed level
# ---------------------------------------------------------------------------
def faint_background(ctx):
    """The unrelated lines of a real list are faint: 100, not the 10000 of
    fake_ctx.  Against a background as bright as the prediction a faint line
    is as unlikely under either hypothesis, and brightness decides nothing."""
    ctx.bg_mu_o = np.full(len(ctx.wn_o), math.log(1.0e2))
    lp.plain_background(ctx)
    return ctx


def test_a_transition_cannot_hide_in_a_feature_fainter_than_itself():
    """A free feature of 100 where 10000 is predicted.

    It could be the transition, under-recorded, and that reading pays for
    the deficit.  It cannot be a transition hidden in a feature that is
    there anyway, because the feature does not have the light.  Before the
    condition the row took the hidden reading's positional credit and was
    tagged `bright`.
    """
    ctx = faint_background(fake_ctx([50000.0], i_obs=[1.0e2]))
    ln_r = _one_prediction(ctx, 1.0e4)
    _, tab = lp.ln_ratio(ctx, 'X', np.array([50000.0]), detail=True)
    assert not bool(tab['over'][0])
    assert ln_r < 0.0


def test_a_feature_brighter_than_predicted_may_still_hide_it():
    ctx = fake_ctx([50000.0], i_obs=[1.0e4])
    _one_prediction(ctx, 1.0)
    _, tab = lp.ln_ratio(ctx, 'X', np.array([50000.0]), detail=True)
    assert bool(tab['over'][0])


def test_a_match_is_a_line_of_the_right_brightness_not_any_line_nearby():
    """One line as bright as predicted, one a hundred times fainter, one a
    wavenumber off: all three are in the matching window, one is a match."""
    ctx = faint_background(fake_ctx([50000.0, 60000.0, 70000.25],
                                    i_obs=[1.0e4, 1.0e2, 1.0e4]))
    ctx.by_level['X'] = dict(
        partner=np.array(['A', 'B', 'C']), sign=np.array([1.0, 1.0, 1.0]),
        i_pred=np.array([1.0e4, 1.0e4, 1.0e4]),
        e_m=np.array([50000.0, 40000.0, 30000.0]),
        u_m=np.array([0.0, 0.0, 0.0]),
        degenerate=np.array([False, False, False]))
    sup = lp.support(ctx, 'X', 100000.0)
    assert sup['n_seen_alt'] == 3
    assert sup['n_real'] == 1
    assert sup['n_free'] == 1


def test_the_unknown_table_splits_the_observable_predictions_three_ways():
    ctx = faint_background(fake_ctx([50000.0, 60000.0],
                                    i_obs=[1.0e4, 1.0e2]))
    ctx.by_level[lp.UNKNOWN_ID] = dict(
        partner=np.array(['A', 'B', 'C']), sign=np.array([1.0, 1.0, 1.0]),
        i_pred=np.array([1.0e4, 1.0e4, 1.0e4]),
        e_m=np.array([50000.0, 40000.0, 25000.0]),
        u_m=np.array([0.0, 0.0, 0.0]),
        degenerate=np.array([False, False, False]))
    r = dict(E_calc=100000.0, W=20.0, positions=[(100000.0, 1.0)])
    row = lp.unknown_table(ctx, r).iloc[0]
    assert (row['n_obs'], row['n_match'], row['n_poor'], row['n_miss']) == \
        (3, 1, 1, 1)


def test_the_window_can_be_widened_for_a_mixed_level():
    ctx = fake_ctx([70000.0, 80000.0, 90000.0])
    en, trans, mapping, e_meas, id_of = theory(LEVELS, TRANS)
    for i in (1, 2, 3):
        ctx.e_final[f'L{i:03d}'] = float(en['E_obs'][i])
        ctx.n_acc_level[f'L{i:03d}'] = 5
    lp.register_unknown(ctx, 4, en, trans, mapping, e_meas, id_of,
                        log=lambda *a: None)
    r = lp.scan_unknown(ctx, 4, en, {'f25f': 20.0}, 132.0, step=0.5, k=4.5)
    assert r['lo'] == pytest.approx(90000.0 - 90.0)
    assert r['hi'] == pytest.approx(90000.0 + 90.0)
    assert r['W'] == 20.0


# ---------------------------------------------------------------------------
# The delta J fingerprint
# ---------------------------------------------------------------------------
def j_tab(rows):
    """A detail table of (partner, P_obs, found) triples.

    Only the columns match_kinds() and ln_j_pattern() read are filled: a found
    row is a well-centred match on a free feature with a positive ln R_t, a
    missing row has no line at all.
    """
    return pd.DataFrame([
        dict(partner=p, P_obs=float(pv), matched=bool(f),
             d=0.0 if f else np.nan, C=0.0 if f else np.nan,
             over=False, W=0.4, ln_R=1.0 if f else -0.5)
        for p, pv, f in rows])


def test_a_class_found_as_often_as_predicted_costs_nothing():
    """The detection rate is fitted, and where it comes out at 1 it is free.

    Three predictions at P_obs = 0.9, 0.8 and 0.7, two of them recorded: 2.4
    were expected and 2 turned up, which is what a level of this J looks like.
    """
    lam, gain = lp.fit_detection_rate([0.9, 0.8, 0.7], [True, True, False])
    assert 0.8 < lam <= 1.0
    assert gain < 0.35


def test_a_class_that_is_never_recorded_is_charged_its_whole_absence():
    """None of them found: the rate goes to zero and the gain is the product
    of the absences, which is the plain probability of missing them all."""
    p = [0.7, 0.6, 0.5]
    lam, gain = lp.fit_detection_rate(p, [False, False, False])
    assert lam == 0.0
    assert gain == pytest.approx(-sum(math.log(1.0 - x) for x in p), abs=1e-9)


def test_a_surplus_of_matches_is_not_evidence_about_j():
    """Found oftener than predicted: the rate is capped at 1 and gains
    nothing.  A class cannot support a position by being over-observed."""
    lam, gain = lp.fit_detection_rate([0.3, 0.3], [True, True])
    assert lam == 1.0
    assert gain == 0.0


def test_a_missing_bright_class_accuses_the_position():
    """Every dJ = -1 partner recorded, every dJ = 0 partner missing.

    That is the fingerprint of a level of a different J sitting at this
    energy, and ln_j_pattern charges for it.
    """
    j_of = {'a': 4.5, 'b': 4.5, 'c': 4.5, 'd': 5.5, 'e': 5.5, 'f': 5.5}
    tab = j_tab([('a', 0.9, True), ('b', 0.8, True), ('c', 0.7, True),
                 ('d', 0.7, False), ('e', 0.7, False), ('f', 0.6, False)])
    r = lp.ln_j_pattern(tab, 5.5, j_of)
    assert r['ln_J'] > 1.5
    got = {c['dJ']: c for c in r['classes']}
    assert got[-1]['n_match'] == 3 and got[-1]['lam'] == 1.0
    assert got[0]['n_match'] == 0 and got[0]['lam'] == 0.0


def test_a_class_too_faint_to_record_cannot_accuse_anything():
    """THE POINT OF WEIGHING THE TEST BY P_obs.

    The same missing dJ = 0 class, but its transitions are predicted far too
    weak to have been recorded.  Nothing is expected of them, so their absence
    says nothing about J and ln_J stays at zero - where a test made of
    selection rules alone would have accused the position exactly as hard as
    it does above.
    """
    j_of = {'a': 4.5, 'b': 4.5, 'c': 4.5, 'd': 5.5, 'e': 5.5, 'f': 5.5}
    tab = j_tab([('a', 0.9, True), ('b', 0.8, True), ('c', 0.7, True),
                 ('d', 0.02, False), ('e', 0.05, False), ('f', 0.01, False)])
    r = lp.ln_j_pattern(tab, 5.5, j_of)
    assert r['ln_J'] == 0.0
    assert [c['dJ'] for c in r['classes']] == [-1]


def test_one_absence_is_not_a_fingerprint():
    """A single missing prediction, however bright, is one line and not a
    pattern: J_CLASS_COST takes the whole of what it gains."""
    j_of = {'a': 4.5, 'd': 5.5}
    tab = j_tab([('a', 0.9, True), ('d', 0.9, False)])
    assert lp.ln_j_pattern(tab, 5.5, j_of)['ln_J'] == 0.0


def test_a_partner_with_no_j_is_left_out_of_the_test():
    j_of = {'a': 4.5, 'b': 4.5}
    tab = j_tab([('a', 0.9, True), ('b', 0.8, True), ('z', 0.7, False)])
    r = lp.ln_j_pattern(tab, 5.5, j_of)
    assert r['n_no_j'] == 1
    assert [c['n_obs'] for c in r['classes']] == [2]


def test_the_level_whose_j_is_unknown_has_no_fingerprint():
    tab = j_tab([('a', 0.9, True), ('d', 0.9, False)])
    r = lp.ln_j_pattern(tab, float('nan'), {'a': 4.5, 'd': 5.5})
    assert r['ln_J'] == 0.0


def test_support_reports_the_fingerprint():
    """support() carries ln_J, so the audit and the unknown table get it."""
    ctx = fake_ctx([50000.0])
    ctx.j_of = {'X': 5.5, 'A': 4.5}
    _one_prediction(ctx, 1.0e4)
    s = lp.support(ctx, 'X', 50000.0)
    assert 'ln_J' in s and s['ln_J'] == 0.0


# ---------------------------------------------------------------------------
# Stage 2: the fingerprint enters ln R, and only between rival positions
# ---------------------------------------------------------------------------
def test_the_cleanest_rival_keeps_its_ln_r_untouched():
    """Nothing is taken from the position whose delta J pattern is best.

    The fingerprint of a single position is dominated by how well the
    calculated gA divides the level's strength among its branches, which is a
    property of the wavefunction and the same at every energy, so it must not
    be charged to anybody's whereabouts.
    """
    assert lp.fold_j([10.0, 8.0], [3.0, 5.0])[0] == pytest.approx(10.0)


def test_a_rival_pays_what_its_fingerprint_is_worse_by():
    assert lp.fold_j([10.0, 8.0], [3.0, 5.0])[1] == pytest.approx(8.0 - 2.0)


def test_equally_dirty_rivals_are_left_exactly_as_they_were():
    """Both positions score 5: the branch strengths are wrong, but they are
    wrong in the same way at both, so there is nothing to choose between."""
    assert lp.fold_j([10.0, 8.0], [5.0, 5.0]) == pytest.approx([10.0, 8.0])


def test_a_position_with_no_rival_pays_nothing():
    assert lp.fold_j([146.0], [6.0]) == pytest.approx([146.0])
    assert lp.fold_j([146.0], []) == pytest.approx([146.0])


def j_ctx():
    """A level of J = 4 with two partners of J = 3 and two of its own J.

    Its four predicted transitions are recorded at 100000 cm^-1; ten
    wavenumbers higher only the two dJ = -1 partners are, so the whole dJ = 0
    class is missing there and nowhere else.
    """
    ctx = faint_background(fake_ctx([50000.0, 50010.0, 60000.0, 60010.0,
                                     70000.0, 80000.0]))
    ctx.by_level[lp.UNKNOWN_ID] = dict(
        partner=np.array(['A', 'B', 'C', 'D']),
        sign=np.ones(4), i_pred=np.full(4, 1.0e4),
        e_m=np.array([50000.0, 40000.0, 30000.0, 20000.0]),
        u_m=np.zeros(4), degenerate=np.zeros(4, dtype=bool))
    ctx.j_of = {str(lp.UNKNOWN_ID): 4.0, 'A': 3.0, 'B': 3.0,
                'C': 4.0, 'D': 4.0}
    return ctx


def test_the_position_that_loses_a_whole_dj_class_scores_it():
    ctx = j_ctx()
    assert lp.support(ctx, lp.UNKNOWN_ID, 100000.0)['ln_J'] == \
        pytest.approx(0.0)
    assert lp.support(ctx, lp.UNKNOWN_ID, 100010.0)['ln_J'] > 1.0


def test_the_unknown_table_judges_the_candidates_on_ln_r_j():
    """Two candidates worth the same in lines; the dJ pattern separates
    them, and it is ln_R_J - not ln_R - that the look-elsewhere correction
    and the order are taken on."""
    ctx = j_ctx()
    r = dict(E_calc=100000.0, W=20.0,
             positions=[(100000.0, 6.0), (100010.0, 6.0)])
    tab = lp.unknown_table(ctx, r)
    clean = tab[tab['E'] == 100000.0].iloc[0]
    dirty = tab[tab['E'] == 100010.0].iloc[0]
    assert clean['ln_R_J'] == pytest.approx(6.0)
    assert dirty['ln_R_J'] == pytest.approx(6.0 - dirty['ln_J'])
    assert dirty['look'] == pytest.approx(dirty['ln_R_J'] - math.log(2))
    assert tab.iloc[0]['E'] == 100000.0


def test_the_fingerprint_can_overturn_the_order_of_the_candidates():
    """The position the lines prefer by two nats is not the one to take if
    it is missing a class of transitions that should have been recorded."""
    ctx = j_ctx()
    r = dict(E_calc=100000.0, W=20.0,
             positions=[(100010.0, 8.0), (100000.0, 6.0)])
    tab = lp.unknown_table(ctx, r)
    assert sorted(tab['ln_R']) == pytest.approx([6.0, 8.0])
    assert tab.iloc[0]['E'] == 100000.0


def test_a_relocation_must_convince_after_the_fingerprint():
    """ln_R_alt is what the lines say; ln_R_alt_J is what is left once the
    alternate has paid for a dJ pattern worse than the adopted position's."""
    row = dict(near_dE=np.nan, dE_alt=40.0, top_share=0.3, n_free=6,
               ln_R_alt=4.0, ln_R_alt_J=4.0, free_gain=12.0, look=5.0,
               n_own=7, n_kept=0)
    assert lp.disposition(row) == 'relocate'
    row['ln_R_alt_J'] = -0.5
    assert lp.disposition(row) == 'weak'


# ---------------------------------------------------------------------------
# The per-line intensity width: sigma^2 = s0^2 + (k u_calc)^2 with
# s0 = a0 + a1 log10 P and k = c log10 S
# ---------------------------------------------------------------------------
def _width_sample(a0, a1, c, n=5000, seed=7, cut=True):
    """Residuals drawn with exactly the width the model assumes, and - when
    `cut` - only those above each line's detection cut-off kept, as the
    recorded lines are."""
    rng = np.random.default_rng(seed)
    m = 3 * n
    u = rng.uniform(0.2, 1.5, m)
    p10 = rng.uniform(-0.5, 4.0, m)
    ls10 = rng.uniform(-6.5, -0.5, m)
    s0 = a0 + a1 * p10
    k = np.maximum(c * ls10, lp.K_FLOOR)
    r = rng.normal(0.0, np.sqrt(s0 ** 2 + (k * u) ** 2))
    rthr = -p10 * math.log(10.0)      # the residual of a line at threshold
    keep = r > rthr if cut else np.ones(m, dtype=bool)
    idx = np.flatnonzero(keep)[:n]
    return r[idx], u[idx], p10[idx], ls10[idx], rthr[idx]


def _fit(r, u, p10, ls10, rthr, s=1.0, fit=None):
    return lp.fit_intensity_width(r, u, p10, ls10, rthr, fit, s,
                                  lambda *_: None)


def test_the_width_model_is_recovered_from_the_run():
    """a0, a1 and c are fitted, not assumed: planted, they come back."""
    w = _fit(*_width_sample(0.90, -0.06, -0.22))
    assert w.a0 == pytest.approx(0.90, abs=0.06)
    assert w.a1 == pytest.approx(-0.06, abs=0.03)
    assert w.c == pytest.approx(-0.22, abs=0.03)


def test_the_detection_cut_off_is_what_keeps_the_fit_honest():
    """The same sample fitted as though nothing below the threshold had been
    lost reads the selection as a narrower, shifted distribution; with the
    cut-off in the likelihood the planted values come back."""
    r, u, p10, ls10, rthr = _width_sample(0.90, -0.06, -0.22, seed=11)
    honest = _fit(r, u, p10, ls10, rthr)
    naive = _fit(r, u, p10, ls10, np.full(len(r), -50.0))
    assert abs(honest.c + 0.22) < abs(naive.c + 0.22)


def test_k_and_s0_never_fall_below_their_floors():
    """k stops at K_FLOOR for the strongest transitions: even their
    calculation is not taken as exact."""
    w = lp.IntensityWidth(0.9, -0.5, -0.22, 1.0, -3.5)
    s0, k = lp.width_terms(w, np.array([0.0, 10.0]), np.array([-5.0, 0.5]))
    assert k[0] == pytest.approx(1.1)
    assert k[1] == pytest.approx(lp.K_FLOOR)
    assert s0[0] == pytest.approx(0.9)
    assert s0[1] == pytest.approx(lp.S0_FLOOR)


def test_a_width_that_carries_nothing_leaves_the_single_width():
    """Residuals of one width whatever u, P and S say: a1 falls to 0, k to
    its floor, and every line gets close to the one width.  (c itself is not
    defined then: any c with c log10 S below K_FLOOR everywhere is the same
    model.)"""
    rng = np.random.default_rng(3)
    n = 4000
    u = rng.uniform(0.2, 1.5, n)
    p10 = rng.uniform(0, 4, n)
    ls10 = rng.uniform(-6.5, -0.5, n)
    r = rng.normal(0.0, 0.9, n)
    w = _fit(r, u, p10, ls10, np.full(n, -50.0), s=0.9)
    assert abs(w.a1) < 0.03
    s0, k = lp.width_terms(w, p10, ls10)
    assert np.all(k < lp.K_FLOOR + 0.05)
    sig = np.sqrt(s0 ** 2 + (k * u) ** 2)
    assert np.all(np.abs(sig - 0.9) < 0.1)


def test_the_fit_is_declined_when_there_is_too_little_to_fit_it_on():
    r, u, p10, ls10, rthr = _width_sample(0.8, 0.0, -0.2, n=50)
    assert _fit(r, u, p10, ls10, rthr, s=1.15) == lp.flat_width(1.15)
    assert (lp.fit_intensity_width(r, None, p10, ls10, rthr, None, 1.15,
                                   lambda *_: None) == lp.flat_width(1.15))


def test_only_the_rows_marked_for_the_fit_are_used():
    """Blended and multiply classified lines are kept out: their residual is
    not theirs alone."""
    r, u, p10, ls10, rthr = _width_sample(0.90, -0.06, -0.22)
    bad = np.zeros(len(r), dtype=bool)
    bad[:1500] = True
    r = np.where(bad, r + 3.0 * u, r)
    w = _fit(r, u, p10, ls10, rthr, fit=~bad)
    assert w.c == pytest.approx(-0.22, abs=0.04)


def test_a_run_with_no_fitted_width_uses_the_pooled_one():
    ctx = lp.Context()
    ctx.s = 0.9
    assert lp.sigma_intensity(ctx, [0.5, 1.5]) == pytest.approx([0.9, 0.9])
    ctx.width = lp.flat_width(0.9)
    assert lp.sigma_intensity(ctx, [0.5, 1.5]) == pytest.approx([0.9, 0.9])


def test_a_transition_with_no_u_calc_is_judged_by_the_pooled_width():
    ctx = lp.Context()
    ctx.s = 1.163
    ctx.width = lp.IntensityWidth(0.902, -0.062, -0.218, 1.0, -3.5)
    got = lp.sigma_intensity(ctx, [np.nan, 0.0, 0.58], [2.0, 2.0, 2.0],
                             [-4.0, -4.0, -4.0])
    assert got[0] == pytest.approx(1.163)
    assert got[1] == pytest.approx(1.163)
    s0 = 0.902 - 0.062 * 2.0
    k = 0.218 * 4.0
    assert got[2] == pytest.approx(math.sqrt(s0 ** 2 + (k * 0.58) ** 2))


def test_a_missing_p_or_s_is_replaced_by_the_median_of_the_fit():
    ctx = lp.Context()
    ctx.s = 1.163
    ctx.width = lp.IntensityWidth(0.902, -0.062, -0.218, 1.0, -3.5)
    a = lp.sigma_intensity(ctx, [0.58], [np.nan], [np.nan])
    b = lp.sigma_intensity(ctx, [0.58], [1.0], [-3.5])
    assert a == pytest.approx(b)


def test_a_weak_transition_believes_its_u_calc_more_than_a_strong_one():
    ctx = lp.Context()
    ctx.s = 1.163
    ctx.width = lp.IntensityWidth(0.902, -0.062, -0.218, 1.0, -3.5)
    weak, strong = lp.sigma_intensity(ctx, [0.6, 0.6], [2.0, 2.0],
                                      [-5.5, -1.0])
    assert weak > strong


def test_a_well_calculated_line_is_held_to_a_tighter_intensity():
    """The same intensity mismatch costs more when the calculation of that
    one transition is a good one, and less when it is a poor one."""
    ctx = fake_ctx([50000.0], i_obs=[1.0e4])
    ctx.s = 1.163
    ctx.width = lp.IntensityWidth(0.902, -0.062, -0.218, 1.0, -3.5)

    def score(u):
        ctx.by_level['X'] = dict(
            partner=np.array(['A']), sign=np.array([1.0]),
            i_pred=np.array([1.0e5]), e_m=np.array([0.0]),
            u_m=np.array([0.0]), u_calc=np.array([float(u)]),
            ls10=np.array([-4.0]), degenerate=np.array([False]))
        return float(lp.ln_ratio(ctx, 'X', np.array([50000.0]))[0][0])

    assert score(0.25) < score(2.0)


def test_line_strength_in_atomic_units():
    """S = 3.0376e-6 lambda gf, gf = 1.499e-16 lambda^2 gA, lambda in A."""
    lam = 2000.0
    gA = 1.0e8
    got = lp.ls.line_strength(np.array([gA, np.nan, 0.0]),
                               np.array([1e8 / lam, 5e4, 5e4]))
    assert got[0] == pytest.approx(3.0376e-6 * lam * 1.499e-16 * lam ** 2 * gA)
    assert np.isnan(got[1]) and np.isnan(got[2])


def test_with_k_at_zero_the_level_offset_is_the_old_correction():
    """The marginalised offset is written in A = sum w/sigma^2 and
    B = sum w r/sigma^2; at one width for every row it must reduce to the
    -0.5 ln(1 + n s_L^2/s^2) + ... of section 3."""
    w = np.array([1.0, 0.6, 0.9])
    r = np.array([0.4, -0.3, 0.8])
    s, s_l = 1.163, 0.217
    A, B = (w / s ** 2).sum(), (w * r / s ** 2).sum()
    new = (-0.5 * math.log(1.0 + s_l ** 2 * A)
           + s_l ** 2 * B ** 2 / (2.0 * (1.0 + s_l ** 2 * A)))
    n, sum_r = w.sum(), (w * r).sum()
    old = (-0.5 * math.log(1.0 + n * s_l ** 2 / s ** 2)
           + s_l ** 2 * sum_r ** 2
           / (2.0 * s ** 2 * (s ** 2 + n * s_l ** 2)))
    assert new == pytest.approx(old)
