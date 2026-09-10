"""Tests for `unfound_levels`.

Run from the LineClass directory:  python -m pytest tests -q

Most of these build a miniature level list - a handful of calculated levels,
some found and some not - and check the pieces that decide the answer: that
only the transitions joining an unfound level to a found one are kept, that the
predicted intensity is the pipeline's own formula, that a transition outside
the recorded wavenumber range can never be promising, and that the probability
really is averaged over the width the level's configuration allows rather than
evaluated at the calculated energy alone.  They use a stand-in for the run
carrying only the four things the probability needs, so they do not read the
pipeline's output.

Two tests marked `real_file` open the real IDEN2 directory and the real
transition-probability workbook; they are skipped when those are absent.
"""
import math
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import classify_lines as cl                # noqa: E402
import cowan_gA                            # noqa: E402
import level_positions as lp               # noqa: E402
import unfound_levels as uf                # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REAL_IDEN2 = os.path.join(HERE, 'IDEN2')

real_file = pytest.mark.skipif(
    not (os.path.exists(os.path.join(REAL_IDEN2, 'enlev.dat'))
         and os.path.exists(cowan_gA.TP_FILE)),
    reason='the real IDEN2 directory or the transition workbook is absent')

C = float(cl.CFG.intensity_model['C'])
KT = float(cl.CFG.intensity_model['kT'])


class FakeRun:
    """The four things `observation_probabilities` asks of the run.

    `calib` is None, which makes `level_positions.p_obs_vec` return the
    coverage alone - a definite function of wavelength and nothing else - so
    these tests do not depend on the fitted noise threshold of any particular
    run.
    """
    calib = None
    bias = {}
    wn_o = np.array([1.0e4, 5.0e4, 1.0e5])
    rho_o = np.array([0.1, 0.1, 0.1])


def enlev_frame(rows):
    """`(en, e_meas)` from (idx, E_calc, E_obs, J, cfg, term, known) rows.

    `en` is what `level_interchange.read_enlev` returns, indexed by the IDEN2
    row number, and `e_meas` the measured energies of the found levels, which
    is what the transitions are built against.
    """
    recs = []
    for idx, e_calc, e_obs, j, cfg, term, known in rows:
        recs.append(dict(idx=idx, E_calc=e_calc, u_obs=0.01, E_obs=e_obs,
                         omc=e_obs - e_calc, J=j, label=cfg + term,
                         cfg=cfg, term=term, known=known))
    en = pd.DataFrame(recs).set_index('idx')
    e_meas = {i: float(en['E_obs'][i]) for i in en.index if en['known'][i]}
    return en, e_meas


def trans_frame(rows):
    """Cowan's transition list from (lid1, lid2, gA) rows."""
    return pd.DataFrame([dict(lid1=a, lid2=b, gA=g) for a, b, g in rows])


# ---------------------------------------------------------------------------
# Which transitions are the level's own
# ---------------------------------------------------------------------------
def test_only_transitions_with_one_end_found_are_kept():
    """A transition between two found levels is not this module's business,
    and one between two unfound levels cannot be given a wavenumber at all -
    neither end is known - so it cannot be searched for."""
    en, e_meas = enlev_frame([
        (1, 1000.0, 1000.5, 2.5, 'f25g', '_3H4H', True),
        (2, 20000.0, 20000.0, 3.5, 'f25g', '_3H4G', False),
        (3, 40000.0, 40000.0, 2.5, 'f26p', '_3F2F', False),
        (4, 60000.0, 60000.7, 3.5, 'f26p', '_3F2G', True),
    ])
    mapping = {10: 1, 20: 2, 30: 3, 40: 4}
    trans = trans_frame([(10, 20, 1.0e8),      # found - unfound: kept
                         (10, 40, 1.0e8),      # found - found: dropped
                         (20, 30, 1.0e8),      # unfound - unfound: dropped
                         (30, 40, 1.0e8)])     # unfound - found: kept
    t = uf.unfound_transitions(en, trans, mapping, e_meas, log=lambda *a: None)
    assert sorted(zip(t['idx'], t['partner'])) == [(2, 1), (3, 4)]


def test_the_unfound_level_may_be_the_upper_or_the_lower_partner():
    """`upper` says which of the two carries the Boltzmann factor."""
    en, e_meas = enlev_frame([
        (1, 1000.0, 1000.0, 2.5, 'f25g', '_3H4H', True),
        (2, 20000.0, 20000.0, 3.5, 'f25g', '_3H4G', False),
        (3, 90000.0, 90000.0, 2.5, 'f25g', '_3H4F', True),
    ])
    t = uf.unfound_transitions(en, trans_frame([(10, 20, 1.0e8),
                                                (20, 30, 1.0e8)]),
                               {10: 1, 20: 2, 30: 3}, e_meas,
                               log=lambda *a: None)
    role = dict(zip(t['partner'], t['upper']))
    assert role[1] == 'L'      # the unfound level at 20000 is above level 1
    assert role[3] == 'M'      # and below level 3


def test_a_transition_with_no_gA_is_dropped():
    """Cowan prints gA only above a cutoff; a transition he does not print is
    not one of unknown strength but one known to be weak."""
    en, e_meas = enlev_frame([
        (1, 1000.0, 1000.0, 2.5, 'f25g', '_3H4H', True),
        (2, 20000.0, 20000.0, 3.5, 'f25g', '_3H4G', False),
    ])
    t = uf.unfound_transitions(en, trans_frame([(10, 20, 0.0)]),
                               {10: 1, 20: 2}, e_meas, log=lambda *a: None)
    assert len(t) == 0


# ---------------------------------------------------------------------------
# How far the level may be from where the calculation puts it
# ---------------------------------------------------------------------------
def test_the_window_is_the_rms_of_its_own_configuration():
    """W is measured per configuration, over the found levels of that
    configuration alone: here 4f^2.5g is out by 3, 4 and 5 cm^-1, so
    W = sqrt((9 + 16 + 25)/3)."""
    en, _ = enlev_frame([
        (1, 1000.0, 1003.0, 2.5, 'f25g', '_3H4H', True),
        (2, 2000.0, 2004.0, 3.5, 'f25g', '_3H4G', True),
        (3, 3000.0, 3005.0, 4.5, 'f25g', '_3H4F', True),
        (4, 4000.0, 4100.0, 2.5, 'f26p', '_3F2F', True),
        (5, 5000.0, 5000.0, 2.5, 'f26p', '_3F2G', False),
    ])
    per_cfg, whole = uf.windows(en, log=lambda *a: None)
    assert per_cfg['f25g'] == pytest.approx(math.sqrt((9 + 16 + 25) / 3.0))
    assert per_cfg['f26p'] == pytest.approx(100.0)
    assert whole == pytest.approx(math.sqrt((9 + 16 + 25 + 10000) / 4.0))


def test_a_configuration_with_no_found_level_takes_the_list_wide_width():
    """It has nothing of its own to average, and the row says so."""
    en, _ = enlev_frame([
        (1, 1000.0, 1003.0, 2.5, 'f25g', '_3H4H', True),
        (2, 2000.0, 2004.0, 3.5, 'f25g', '_3H4G', True),
        (3, 9000.0, 9000.0, 2.5, 'f5d6d', '_3P2P', False),
    ])
    per_cfg, whole = uf.windows(en, log=lambda *a: None)
    t = pd.DataFrame(dict(idx=[3], partner=[1], gA=[1.0], E_calc=[9000.0],
                          E_partner=[1003.0], upper=['L']))
    t = uf.attach_windows(t, en, per_cfg, whole)
    assert t['W_src'].iloc[0] == 'list'
    assert t['W'].iloc[0] == pytest.approx(whole)


# ---------------------------------------------------------------------------
# The predicted intensity and the probability of having been recorded
# ---------------------------------------------------------------------------
def test_the_predicted_intensity_is_the_pipeline_s_own_formula():
    """I = C gA (nu/1e8) exp(-E_up/kT) with the C and kT of the configuration
    file - the same numbers Icalc.xlsx carries for the transitions between two
    found levels."""
    e_lo, e_up, gA = 12000.0, 52000.0, 3.0e7
    t = pd.DataFrame(dict(idx=[2], partner=[1], gA=[gA], E_calc=[e_up],
                          E_partner=[e_lo], upper=['L'], W=[10.0],
                          W_src=['cfg']))
    out = uf.observation_probabilities(t, FakeRun())
    nu = e_up - e_lo
    assert out['wn'].iloc[0] == pytest.approx(nu)
    assert out['lam'].iloc[0] == pytest.approx(1.0e8 / nu)
    assert out['I_pred'].iloc[0] == pytest.approx(
        C * gA * (nu / 1.0e8) * math.exp(-e_up / KT))


def test_a_transition_outside_the_recorded_range_is_never_promising():
    """Sugar's list runs from WN_MIN to WN_MAX; a transition outside it could
    not have been recorded whatever its strength."""
    e_lo = 1000.0
    e_up = e_lo + cl.WN_MAX + 5000.0
    t = pd.DataFrame(dict(idx=[2], partner=[1], gA=[1.0e10], E_calc=[e_up],
                          E_partner=[e_lo], upper=['L'], W=[1.0],
                          W_src=['cfg']))
    out = uf.observation_probabilities(t, FakeRun())
    assert out['P_obs'].iloc[0] == 0.0
    assert out['P_at_calc'].iloc[0] == 0.0


def test_the_probability_is_averaged_over_the_window():
    """The level's energy is not known, so P_obs is the weighted mean of the
    probability over the trial positions E_calc + z W, not its value at E_calc.

    The mean is checked against the same function evaluated by hand at each of
    those positions - a level given no window at all (W = 0) reports its
    probability at the energy it is put at - and then against the value at
    E_calc, which it must differ from wherever the neighbourhood is not
    uniform.  Here the level straddles the lower edge of the recorded range, so
    some of its trial positions have no transition to record at all.
    """
    e_lo = 1000.0
    e_up = e_lo + cl.WN_MIN + 100.0        # 100 cm^-1 inside the lower edge
    W = 400.0

    def one(E, w):
        t = pd.DataFrame(dict(idx=[2], partner=[1], gA=[1.0e10], E_calc=[E],
                              E_partner=[e_lo], upper=['L'], W=[w],
                              W_src=['cfg']))
        return uf.observation_probabilities(t, FakeRun())

    out = one(e_up, W)
    by_hand = sum(w * one(e_up + z * W, 0.0)['P_at_calc'].iloc[0]
                  for z, w in zip(uf.Z_GRID, uf.Z_WEIGHT))
    assert out['P_obs'].iloc[0] == pytest.approx(by_hand)
    assert abs(out['P_obs'].iloc[0] - out['P_at_calc'].iloc[0]) > 1e-6


# ---------------------------------------------------------------------------
# The table
# ---------------------------------------------------------------------------
def _mini_table(p_values, cfg='f25g'):
    """A one-level table whose transitions carry the given probabilities."""
    n = len(p_values)
    t = pd.DataFrame(dict(idx=[7] * n, partner=list(range(1, n + 1)),
                          gA=[1.0] * n, E_calc=[5000.0] * n,
                          E_partner=[100.0] * n, upper=['L'] * n,
                          W=[20.0] * n, W_src=['cfg'] * n,
                          wn=[4900.0] * n, lam=[2.0e4] * n,
                          I_pred=list(range(1, n + 1)),
                          P_at_calc=list(p_values), P_obs=list(p_values),
                          rho=[0.1] * n))
    en = pd.DataFrame([dict(idx=7, E_calc=5000.0, u_obs=5000.0, E_obs=5000.0,
                            omc=0.0, J=2.5, label=cfg + '_3H4H', cfg=cfg,
                            term='_3H4H', known=False)]).set_index('idx')
    return t, en


def test_n_prom_counts_only_the_promising_transitions():
    """n_prom is the count at P_PROMISING, n_obs the weaker bar
    level_positions.py uses, and sum_P the expected number recorded."""
    p = [0.9, 0.6, 0.5, 0.3, 0.21, 0.1, 0.0]
    t, en = _mini_table(p)
    tab = uf.level_table(t, en)
    row = tab.iloc[0]
    assert row['n_prom'] == 3                       # 0.9, 0.6, 0.5
    assert row['n_obs'] == 5                        # + 0.3, 0.21
    assert row['n_pred'] == len(p)
    assert row['sum_P'] == pytest.approx(sum(p), abs=5e-3)
    assert row['parity'] == 'e'                     # 4f^2.5g


def test_the_table_is_sorted_by_the_count_that_matters():
    """Best first, and ties broken by the expected number of lines."""
    t1, en1 = _mini_table([0.9, 0.9], cfg='f25g')
    t2, en2 = _mini_table([0.9, 0.9, 0.9], cfg='f26p')
    t3, en3 = _mini_table([0.6, 0.6], cfg='f27s')
    t2['idx'] = 8
    t3['idx'] = 9
    en2 = en2.rename(index={7: 8})
    en3 = en3.rename(index={7: 9})
    tab = uf.level_table(pd.concat([t1, t2, t3], ignore_index=True),
                         pd.concat([en1, en2, en3]))
    assert list(tab['idx']) == [8, 7, 9]


# ---------------------------------------------------------------------------
# The real files
# ---------------------------------------------------------------------------
@real_file
def test_the_calculated_levels_line_up_with_the_real_enlev():
    """Every level of enlev.dat is matched to a level of Cowan's transition
    list, and every one whose identifier both files know is checked."""
    en, trans, mapping = uf.read_theory(log=lambda *a: None)
    assert len(en) == len(mapping)
    assert int(en['known'].sum()) + int((~en['known']).sum()) == len(en)


@real_file
def test_the_report_covers_the_unfound_levels_and_only_those():
    """A full run of the real files: every level in the table is one that has
    not been found, and the transitions it is judged on all lead to levels
    that have."""
    args = uf.parse_args(['--quiet-run'])
    _, en, t, tab, _ = uf.build(args, log=lambda *a: None)
    assert len(tab) > 0
    assert not en.loc[tab['idx'], 'known'].any()
    assert en.loc[t['partner'].unique(), 'known'].all()
    assert (tab['n_prom'] <= tab['n_obs']).all()
    assert (tab['n_obs'] <= tab['n_pred']).all()
    assert tab['n_prom'].iloc[0] == tab['n_prom'].max()
