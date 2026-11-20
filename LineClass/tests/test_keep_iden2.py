"""Tests for classify_lines.py --keep-IDEN2-ids (keep_iden2.py, 2026-11-20).

The mode holds a classification run to the identifications of the set's
IDEN2/trans.dat: every one is accepted, every other candidate rejected, and
each enforced verdict that changes anything is written into the decision
ledger and the log.  The tests fix:

  * the join of IDEN2's identifications to the run's lines and levels, and
    that a broken join stops;
  * how IDEN2's verdicts are laid over the ledger's, including a ledger row
    that contradicts IDEN2;
  * the weeding's last word under the mode (apply_line_decisions);
  * which enforced verdicts give a ledger row, and that a failure to obey
    IDEN2 stops before anything is written;
  * the ledger file: rows appended and reversed, its line ending kept;
  * a forced offset IDEN2 holds is warned of, not stopped on;
  * the option on the command line.

Run from the LineClass directory:  python -m pytest tests -q
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import classify_lines as cl                                     # noqa: E402
import keep_iden2 as K                                          # noqa: E402
from models import EnergyLevel, SpectralLine, Transition, UNASSIGNED  # noqa: E402

A, R = cl.DECISIONS_ACCEPT, cl.DECISIONS_REJECT


def level(lid, energy, parity, j=2.5):
    return EnergyLevel(level_id=lid, energy=energy, parity=parity,
                       J_str=str(j), J_val=float(j))


def obs(wn, key=None):
    return SpectralLine(wavenumber=wn, wn_uncertainty=0.01, intensity=10.0,
                        line_character='', wn_key=key if key else wn)


def cand(lo, up, line, accepted=None, forced=False):
    t = Transition(lower_level=lo, upper_level=up, calc_intensity=5.0,
                   u_calc=0.5, assigned_to=line, accepted=accepted)
    t.forced = forced
    line.assigned_transitions.append(t)
    return t


@pytest.fixture
def world():
    """Levels 000100 (1000, e), 000200 (1100, o), 000300 (1250, o), and two
    lines: 100.003 (named 100.0) and 250.0."""
    lv = {x.level_id: x for x in (level('000100', 1000.0, 'e'),
                                  level('000200', 1100.0, 'o'),
                                  level('000300', 1250.0, 'o'))}
    l1, l2 = obs(100.003, 100.0), obs(250.0)
    return lv, l1, l2


# ---- the join ----------------------------------------------------------------

def test_a_line_is_found_by_its_key_and_the_levels_ordered_by_energy(world):
    lv, l1, l2 = world
    got = K.join([(7, 100.0, 99.9, '000200', '000100')], [l1, l2], lv, 0.01)
    assert got == {id(l1): (l1, {('000100', '000200')})}


def test_a_row_without_a_key_is_found_by_its_wavenumber(world):
    lv, l1, l2 = world
    got = K.join([(8, None, 250.004, '000100', '000300')], [l1, l2], lv, 0.01)
    assert got[id(l2)][1] == {('000100', '000300')}


@pytest.mark.parametrize('row, what', [
    ((9, None, 180.0, '000100', '000200'), 'no observed line'),
    ((9, 100.0, 100.0, '000100', '000999'), 'not in the level list')])
def test_a_broken_join_stops(world, row, what):
    lv, l1, l2 = world
    with pytest.raises(ValueError, match=what):
        K.join([row], [l1, l2], lv, 0.01)


# ---- laying IDEN2 over the ledger -------------------------------------------

def lay(monkeypatch, world, iden2, ledger):
    lv, l1, l2 = world
    monkeypatch.setattr(K, 'read_assignments', lambda d: iden2)
    sites = []
    for (line, low, upp, verdict) in ledger:
        line.decisions[(low, upp)] = verdict
        sites.append(((line.wn_key, low, upp), line, verdict))
    st = K.lay_on([l1, l2], lv, sites, 'IDEN2', A, R, 0.01)
    return st, sites


def test_iden2_verdicts_are_laid_over_the_ledger(monkeypatch, world):
    lv, l1, l2 = world
    st, sites = lay(monkeypatch, world,
                    [(1, 100.0, 100.0, '000100', '000200')],
                    [(l1, '000100', '000300', (A, 'by hand')),     # not in IDEN2
                     (l2, '000100', '000300', (R, 'too weak'))])    # agrees
    assert l1.iden2_pairs == {('000100', '000200')}
    assert l2.iden2_pairs == set()
    # the IDEN2 assignment is accepted, provisionally
    assert l1.decisions[('000100', '000200')] == (A, K.ACCEPT_REASON)
    assert st.provisional == {(id(l1), ('000100', '000200'))}
    # the contradicting ledger row is reversed, in the sites as well
    assert l1.decisions[('000100', '000300')] == (R, K.REJECT_REASON)
    assert st.reversed == [((100.0, '000100', '000300'), l1, R, A, 'by hand')]
    assert ((100.0, '000100', '000300'), l1, (R, K.REJECT_REASON)) in sites
    # an agreeing row is left alone
    assert l2.decisions[('000100', '000300')] == (R, 'too weak')
    assert len(sites) == 3


def test_a_ledger_rejection_of_an_iden2_assignment_is_reversed(monkeypatch,
                                                               world):
    lv, l1, l2 = world
    st, _ = lay(monkeypatch, world, [(1, 100.0, 100.0, '000100', '000200')],
                [(l1, '000100', '000200', (R, 'old'))])
    assert l1.decisions[('000100', '000200')] == (A, K.ACCEPT_REASON)
    assert st.provisional == set()
    assert st.reversed[0][2:] == (A, R, 'old')


# ---- the weeding's last word ------------------------------------------------

def test_under_the_mode_every_other_candidate_is_rejected(world):
    lv, l1, _ = world
    keep = cand(lv['000100'], lv['000200'], l1, accepted=0)
    other = cand(lv['000100'], lv['000300'], l1, accepted=1)
    other.notes2 = 'Step1 Accepted'
    l1.iden2_pairs = {('000100', '000200')}
    l1.decisions[('000100', '000200')] = (A, K.ACCEPT_REASON)
    cl.apply_line_decisions(l1)
    assert keep.accepted == 1 and keep.algo_accepted == 0
    assert other.accepted == 0 and other.algo_accepted == 1
    assert other.algo_note == 'Step1 Accepted'
    assert other.notes2.startswith('Kept IDEN2 rejection (overrides)')


def test_a_line_iden2_leaves_empty_loses_all_its_candidates(world):
    lv, l1, _ = world
    t = cand(lv['000100'], lv['000200'], l1, accepted=1)
    l1.iden2_pairs = set()
    cl.apply_line_decisions(l1)
    assert t.accepted == 0


def test_the_manual_column_says_where_a_verdict_comes_from(world):
    lv, l1, _ = world
    lo, up, up3 = lv['000100'], lv['000200'], lv['000300']
    held = cand(lo, up, l1, accepted=0)             # laid on from IDEN2
    hand = cand(lo, up3, l1, accepted=1)            # a hand-written row
    l1.iden2_pairs = {('000100', '000200'), ('000100', '000300')}
    l1.decisions[('000100', '000200')] = (A, K.ACCEPT_REASON)
    l1.decisions[('000100', '000300')] = (A, 'by eye in IDEN2')
    cl.apply_line_decisions(l1)
    assert held.manual == cl.IDEN2_MARK and hand.manual == A
    # a ledger row the mode wrote, read in a run without the mode
    _, l2, _ = world
    l2.decisions[('000100', '000200')] = (
        R, K.REJECT_REASON + ' (was accept: old)')
    t = cand(lo, up, l2, accepted=1)
    cl.apply_line_decisions(l2)
    assert t.manual == cl.IDEN2_MARK and t.accepted == 0
    assert cl.origin_of(t) == 'IDEN2'


def test_a_rejection_the_mode_enforces_is_marked_only_if_it_overrides(world):
    lv, l1, _ = world
    was_on = cand(lv['000100'], lv['000200'], l1, accepted=1)
    was_off = cand(lv['000100'], lv['000300'], l1, accepted=0)
    l1.iden2_pairs = set()
    cl.apply_line_decisions(l1)
    assert was_on.manual == cl.IDEN2_MARK
    assert was_off.manual == ''
    assert cl.origin_of(was_off) == 'new'


def test_without_the_mode_nothing_changes(world):
    lv, l1, _ = world
    t = cand(lv['000100'], lv['000200'], l1, accepted=1)
    cl.apply_line_decisions(l1)
    assert t.accepted == 1 and t.algo_accepted is None


# ---- what is written --------------------------------------------------------

def state(l1, l2, provisional=()):
    st = K.KeepState()
    st.provisional = set(provisional)
    for line in (l1, l2):
        if line.iden2_pairs is None:
            line.iden2_pairs = set()
    return st


def test_the_ledger_rows_follow_what_the_mode_changed(world):
    lv, l1, l2 = world
    lo, up, up3 = lv['000100'], lv['000200'], lv['000300']
    l1.iden2_pairs = {('000100', '000200'), ('000100', '000300')}
    forced = cand(lo, up, l1, accepted=1, forced=True)
    forced.algo_accepted = 0
    rejected = cand(lo, up3, l1, accepted=1)
    rejected.algo_accepted, rejected.algo_note = 0, 'Step3 Rejected: new'
    l2.iden2_pairs = set()
    extra = cand(lo, up3, l2, accepted=0)
    extra.algo_accepted = 1
    agreed = cand(lo, up, l2, accepted=0)
    agreed.algo_accepted = 0
    st = state(l1, l2, {(id(l1), ('000100', '000200')),
                        (id(l1), ('000100', '000300'))})
    rows, log = K.finish([l1, l2], st, A, R, 'offset-ok', UNASSIGNED)
    assert sorted((r[0].wn_key, r[1], r[2], r[3], r[4]) for r in rows) == [
        (100.0, '000100', '000200', A, K.ACCEPT_REASON),
        (100.0, '000100', '000300', A, K.ACCEPT_REASON),
        (250.0, '000100', '000300', R, K.REJECT_REASON)]
    cats = sorted(r[0] for r in log)
    assert cats == sorted([K.MISSING, K.REJECTED, K.EXTRA])


def test_an_assignment_the_run_accepts_anyway_needs_no_row(world):
    lv, l1, l2 = world
    l1.iden2_pairs = {('000100', '000200')}
    t = cand(lv['000100'], lv['000200'], l1, accepted=1)
    t.algo_accepted = 1
    st = state(l1, l2, {(id(l1), ('000100', '000200'))})
    assert K.finish([l1, l2], st, A, R, 'offset-ok', UNASSIGNED) == ([], [])


def test_an_offset_iden2_holds_is_written_offset_ok(world):
    lv, l1, l2 = world
    l1.iden2_pairs = {('000100', '000200')}
    cand(lv['000100'], lv['000200'], l1, accepted=1, forced=True)
    st = state(l1, l2, {(id(l1), ('000100', '000200'))})
    st.offset = [(l1, ('000100', '000200'), 0.9)]
    rows, log = K.finish([l1, l2], st, A, R, 'offset-ok', UNASSIGNED)
    assert rows[0][4] == K.ACCEPT_REASON + ' (offset-ok)'
    assert K.OFFSET in [r[0] for r in log]


@pytest.mark.parametrize('accepted, what', [(0, 'not accepted'),
                                            (None, 'not a candidate')])
def test_an_iden2_assignment_not_accepted_stops(world, accepted, what):
    lv, l1, l2 = world
    l1.iden2_pairs = {('000100', '000200')}
    if accepted is not None:
        cand(lv['000100'], lv['000200'], l1, accepted=accepted)
    with pytest.raises(ValueError, match=what):
        K.finish([l1, l2], state(l1, l2), A, R, 'offset-ok', UNASSIGNED)


def test_an_accepted_assignment_iden2_lacks_stops(world):
    lv, l1, l2 = world
    l1.iden2_pairs = set()
    cand(lv['000100'], lv['000200'], l1, accepted=1)
    with pytest.raises(ValueError, match='not identified in IDEN2'):
        K.finish([l1, l2], state(l1, l2), A, R, 'offset-ok', UNASSIGNED)


# ---- the ledger file ----------------------------------------------------------

def test_the_ledger_is_appended_and_reversed_in_its_own_line_ending(tmp_path,
                                                                    world):
    lv, l1, l2 = world
    path = tmp_path / 'line_decisions.csv'
    path.write_bytes(b'wn_key,low_id,upp_id,decision,date,reason\r\n'
                     b'100.0,000100,000300,accept,9/3/2026,"by hand, once"\r\n'
                     b'250.0,000100,000200,reject,9/3/2026,weak\r\n')
    st = K.KeepState()
    st.reversed = [((100.0, '000100', '000300'), l1, R, A, 'by hand, once')]
    rows = [(l2, '000100', '000300', A, K.ACCEPT_REASON)]
    assert K.update_ledger(str(path), rows, st, A, 'offset-ok') == (1, 1)
    text = path.read_bytes().decode('utf-8')
    assert '\n' not in text.replace('\r\n', '')
    lines = text.split('\r\n')
    assert lines[1].startswith('100.0,000100,000300,reject,' + K.today())
    assert 'kept IDEN2 rejection (was accept: by hand, once)' in lines[1]
    assert lines[2] == '250.0,000100,000200,reject,9/3/2026,weak'
    assert lines[3] == ('250.0000,000100,000300,accept,%s,kept IDEN2 '
                        'assignment' % K.today())
    # and the ledger reads back
    dmap = cl.read_line_decisions(str(path))
    assert dmap[(250.0, '000100', '000300')] == (A, K.ACCEPT_REASON)


def test_nothing_to_write_leaves_the_ledger_alone(tmp_path):
    path = tmp_path / 'line_decisions.csv'
    path.write_bytes(b'wn_key,low_id,upp_id,decision\n')
    assert K.update_ledger(str(path), [], K.KeepState(), A, 'x') == (0, 0)
    assert path.read_bytes() == b'wn_key,low_id,upp_id,decision\n'


def test_the_log_lists_every_enforced_decision(tmp_path, world):
    lv, l1, l2 = world
    st = K.KeepState()
    st.n_iden2 = 2
    log = [(K.EXTRA, l2, '000100', '000300', 'Step1 Accepted'),
           (K.MISSING, l1, '000100', '000200', '')]
    path = tmp_path / K.LOG_NAME
    count = K.write_log(str(path), log, st, 'line_decisions.csv', 2, 0)
    assert count[K.EXTRA] == 1 and count[K.MISSING] == 1
    text = path.read_text(encoding='utf-8')
    assert '2 identifications read from IDEN2' in text
    assert '000100 - 000300 Step1 Accepted' in text


# ---- the forced offset ----------------------------------------------------------

def test_an_offset_iden2_holds_is_warned_of_not_stopped_on(monkeypatch, world):
    lv, l1, _ = world
    monkeypatch.setattr(cl, 'WN_MIN', 0.0)
    monkeypatch.setattr(cl, 'WN_MAX', 1.0e6)
    monkeypatch.setattr(cl, 'MAX_FORCED_OFFSET', 0.5)
    key = (99.0, '000100', '000200')         # Ritz 100.0, 1.0 off
    sites = [(key, l1, (A, K.ACCEPT_REASON))]
    with pytest.raises(ValueError, match='misses this line'):
        cl.check_forced_decisions(sites, lv)
    assert cl.check_forced_decisions(sites, lv, {key}) == [(key, 1.0)]


# ---- the command line -----------------------------------------------------------

@pytest.mark.parametrize('flag', ['--keep-IDEN2-ids', '--keep_IDEN2_ids'])
def test_the_option_reaches_main(monkeypatch, flag):
    got = {}
    monkeypatch.setattr(cl, 'apply_config', lambda *a, **k: None)
    monkeypatch.setattr(cl.config, 'load', lambda *a, **k: None)
    monkeypatch.setattr(cl, 'main', lambda **k: got.update(k))
    cl.cli([flag])
    assert got['keep_iden2'] is True
    cl.cli([])
    assert got['keep_iden2'] is False


def test_the_mode_refuses_a_calibration_run():
    with pytest.raises(ValueError, match='real classification'):
        cl.main(wn_shift=1.0, keep_iden2=True, write_files=False)
