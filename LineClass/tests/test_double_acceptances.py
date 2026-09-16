"""Tests for the two safeguards against a level fitted away from its lines.

A pair of levels emits at one wavenumber, so the same pair accepted on two
different observed lines is never a physical result.  It arises when a level is
dragged off its position by a heavy line the ledger holds on it, while its
published identifications stay accepted where they were written.  Two checks
stand against that:

  * check_forced_decisions() refuses a ledger "accept" row whose Ritz
    wavenumber, at the energies the run starts from, misses the line it names
    by more than decisions.max_forced_offset - the stale row, written before
    one of its levels was re-positioned;
  * check_double_acceptances() refuses the outcome itself, whatever produced
    it, and reports the weights that moved the level.

Run from the LineClass directory:  python -m pytest tests -q
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import classify_lines as cl
from models import EnergyLevel, SpectralLine, Transition


def level(lid, energy, parity, j):
    return EnergyLevel(level_id=lid, energy=energy, parity=parity,
                       J_str=str(j), J_val=float(j))


@pytest.fixture(autouse=True)
def plain_world(monkeypatch):
    """A run with no level moved and the shipped offset limit in force."""
    monkeypatch.setattr(cl, 'WN_MIN', 0.0)
    monkeypatch.setattr(cl, 'WN_MAX', 1.0e6)
    monkeypatch.setattr(cl, 'MAX_FORCED_OFFSET', 5.0)
    monkeypatch.setattr(cl, 'LEVEL_REVISIONS', {})
    monkeypatch.setattr(cl, 'LEVEL_OVERRIDES', 'revised_level_energies.csv')


def two_levels(e_low=1000.0, e_up=1100.0):
    lo = level('000100', e_low, 'e', 2.5)
    up = level('000200', e_up, 'o', 2.5)
    return {lo.level_id: lo, up.level_id: up}, lo, up


def site(wn, low, upp, line, decision, reason=''):
    """One entry of the list attach_line_decisions returns."""
    line.decisions[(low, upp)] = (decision, reason)
    return ((wn, low, upp), line, (decision, reason))


# --------------------------------------------------------------------------
# the stale ledger row
# --------------------------------------------------------------------------
def test_an_accept_row_far_from_its_line_stops_the_run():
    levels, lo, up = two_levels()
    line = SpectralLine(wavenumber=80.0, wn_uncertainty=0.01, intensity=50.0,
                        line_character='')
    sites = [site(80.0, '000100', '000200', line, cl.DECISIONS_ACCEPT)]

    with pytest.raises(ValueError) as exc:
        cl.check_forced_decisions(sites, levels)

    msg = str(exc.value)
    assert 'misses this line by +20.0000 cm^-1' in msg
    assert 'max_forced_offset' in msg
    assert 'offset-ok' in msg


def test_the_message_names_the_move_that_made_the_row_stale(monkeypatch):
    levels, lo, up = two_levels()
    monkeypatch.setattr(cl, 'LEVEL_REVISIONS',
                        {'000200': (1080.0, 1100.0,
                                    're-positioned; taken from IDEN2')})
    line = SpectralLine(wavenumber=80.0, wn_uncertainty=0.01, intensity=50.0,
                        line_character='')
    sites = [site(80.0, '000100', '000200', line, cl.DECISIONS_ACCEPT)]

    with pytest.raises(ValueError) as exc:
        cl.check_forced_decisions(sites, levels)

    msg = str(exc.value)
    assert '000200 was moved by revised_level_energies.csv' in msg
    assert '1080.0000 -> 1100.0000 (+20.0000 cm^-1)' in msg
    assert 're-positioned; taken from IDEN2' in msg


def test_a_row_naming_no_moved_level_says_so():
    levels, lo, up = two_levels()
    line = SpectralLine(wavenumber=80.0, wn_uncertainty=0.01, intensity=50.0,
                        line_character='')
    sites = [site(80.0, '000100', '000200', line, cl.DECISIONS_ACCEPT)]

    with pytest.raises(ValueError) as exc:
        cl.check_forced_decisions(sites, levels)

    assert 'Neither level has been moved' in str(exc.value)


def test_offset_ok_in_the_reason_exempts_that_one_row():
    levels, lo, up = two_levels()
    line = SpectralLine(wavenumber=80.0, wn_uncertainty=0.01, intensity=50.0,
                        line_character='')
    sites = [site(80.0, '000100', '000200', line, cl.DECISIONS_ACCEPT,
                  reason='pulls the level to its new place, offset-ok')]

    cl.check_forced_decisions(sites, levels)     # no exception


def test_an_offset_within_the_limit_passes():
    levels, lo, up = two_levels()
    line = SpectralLine(wavenumber=97.0, wn_uncertainty=0.01, intensity=50.0,
                        line_character='')
    sites = [site(97.0, '000100', '000200', line, cl.DECISIONS_ACCEPT)]

    cl.check_forced_decisions(sites, levels)     # 3 cm^-1 out, allowed


def test_a_reject_row_is_not_measured_against_the_limit():
    levels, lo, up = two_levels()
    line = SpectralLine(wavenumber=80.0, wn_uncertainty=0.01, intensity=50.0,
                        line_character='')
    sites = [site(80.0, '000100', '000200', line, cl.DECISIONS_REJECT)]

    cl.check_forced_decisions(sites, levels)     # a refusal rules on nothing


# --------------------------------------------------------------------------
# the outcome: one transition on two lines
# --------------------------------------------------------------------------
def accepted(lo, up, line, new=1, notes2='', manual='', u=0.3):
    """An accepted assignment of the pair (lo, up) to `line`, wired both ways."""
    line.wn_uncertainty = u
    t = Transition(lower_level=lo, upper_level=up, calc_intensity=10.0,
                   u_calc=0.5, assigned_to=line, accepted=1, new=new,
                   grade='5A', notes2=notes2, manual=manual)
    line.assigned_transitions.append(t)
    up.from_transitions.append(t)
    lo.to_transitions.append(t)
    return t


def doubled_world():
    """000389 moved +8.8 cm^-1, its pair accepted at the old and new places."""
    lo = level('000054', 0.0, 'e', 4.0)
    up = level('000389', 116462.04, 'o', 3.0)
    partner = level('000578', 134406.98, 'e', 3.0)
    old = SpectralLine(wavenumber=99937.020, wn_uncertainty=0.3,
                       intensity=1000.0, line_character='')
    new = SpectralLine(wavenumber=99947.514, wn_uncertainty=0.3,
                       intensity=35.0, line_character='b')
    ir = SpectralLine(wavenumber=17944.636, wn_uncertainty=0.0111,
                      intensity=500.0, line_character='')
    accepted(lo, up, old, new=0,
             notes2='Step3 Accepted: old, no solid evidence for rejection')
    accepted(lo, up, new, new=1, notes2='Step1 Accepted New: Ritz Ok')
    heavy = accepted(up, partner, ir, new=1, manual='accept', u=0.0111)
    lines = [old, new, ir]
    weights = {id(t): 1.0 / t.assigned_to.wn_uncertainty ** 2
               for line in lines for t in line.assigned_transitions}
    return lines, weights, heavy


def test_a_pair_on_two_lines_stops_the_run():
    lines, weights, _heavy = doubled_world()
    inputs = {'000054': 0.0, '000389': 116453.20, '000578': 134406.98}

    with pytest.raises(ValueError) as exc:
        cl.check_double_acceptances(lines, inputs, weights)

    msg = str(exc.value)
    assert '1 transition(s) accepted on more than one observed line' in msg
    assert '000054-000389  accepted on 2 lines' in msg
    assert '99937.0200' in msg and '99947.5140' in msg
    assert 'Step3 Accepted: old, no solid evidence for rejection' in msg


def test_the_report_shows_the_shift_and_the_weights_behind_it():
    lines, weights, _heavy = doubled_world()
    inputs = {'000054': 0.0, '000389': 116453.20, '000578': 134406.98}

    with pytest.raises(ValueError) as exc:
        cl.check_double_acceptances(lines, inputs, weights)

    msg = str(exc.value)
    assert '000389  116453.2000 -> 116462.0400 (+8.8400 cm^-1)' in msg
    # the infrared line, 27 times lighter in uncertainty and 700 times heavier
    # in weight, is at the top of the list and is marked as the ledger's
    assert 'heaviest first' in msg
    first = [ln for ln in msg.splitlines() if '17944.6360' in ln]
    assert first and 'ledger: accept' in first[0]
    # a level the fit did not move gets no weight table
    assert '000054  0.0000 -> 0.0000 (+0.0000 cm^-1)' in msg


def test_a_calibration_run_reports_but_does_not_stop(capsys):
    lines, weights, _heavy = doubled_world()
    inputs = {'000054': 0.0, '000389': 116453.20, '000578': 134406.98}

    n = cl.check_double_acceptances(lines, inputs, weights, strict=False)

    assert n == 1
    assert 'Warning:' in capsys.readouterr().out


def test_a_clean_run_passes():
    lo = level('000054', 0.0, 'e', 4.0)
    up = level('000389', 116453.20, 'o', 3.0)
    line = SpectralLine(wavenumber=116453.20, wn_uncertainty=0.3,
                        intensity=1000.0, line_character='')
    accepted(lo, up, line, new=0)

    assert cl.check_double_acceptances([line], {}, {}) == 0


def test_a_rejected_second_assignment_is_not_a_double():
    lo = level('000054', 0.0, 'e', 4.0)
    up = level('000389', 116453.20, 'o', 3.0)
    keep = SpectralLine(wavenumber=116453.20, wn_uncertainty=0.3,
                        intensity=1000.0, line_character='')
    drop = SpectralLine(wavenumber=116462.04, wn_uncertainty=0.3,
                        intensity=35.0, line_character='')
    accepted(lo, up, keep, new=0)
    t = accepted(lo, up, drop, new=1)
    t.accepted = 0

    assert cl.check_double_acceptances([keep, drop], {}, {}) == 0
