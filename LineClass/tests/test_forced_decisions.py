"""Tests for the decision ledger as a source of orders, not votes.

An "accept" row in line_decisions.csv names an identification the analyst has
made by eye.  Some of those sit further from the Ritz wavenumber than the
automatic matching tolerance reaches, so the run would never propose them; the
ledger therefore has to *create* the assignment, and a row that cannot be
carried out has to stop the run instead of being silently dropped.  These
tests cover both halves.

Run from the LineClass directory:  python -m pytest tests -q
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import classify_lines as cl
from models import EnergyLevel, SpectralLine, Transition, UNASSIGNED


def level(lid, energy, parity, j):
    return EnergyLevel(level_id=lid, energy=energy, parity=parity,
                       J_str=str(j), J_val=float(j))


@pytest.fixture
def world(monkeypatch):
    """Two opposite-parity levels 100 cm^-1 apart and the line that joins them."""
    monkeypatch.setattr(cl, 'WN_MIN', 0.0)
    monkeypatch.setattr(cl, 'WN_MAX', 1.0e6)
    lo = level('000100', 1000.0, 'e', 2.5)
    up = level('000200', 1100.0, 'o', 2.5)
    levels = {lo.level_id: lo, up.level_id: up}
    line = SpectralLine(wavenumber=100.3, wn_uncertainty=0.01, intensity=50.0,
                        line_character='')
    cand = Transition(lower_level=lo, upper_level=up, calc_intensity=42.0,
                      u_calc=0.7, orig_calc_intensity=42.0, orig_u_calc=0.7)
    return levels, line, cand


def site(wn, low, upp, line, decision, reason=''):
    """One entry of the list attach_line_decisions returns."""
    line.decisions[(low, upp)] = (decision, reason)
    return ((wn, low, upp), line, (decision, reason))


# --------------------------------------------------------------------------
# creating the assignment
# --------------------------------------------------------------------------
def test_an_accepted_pair_the_matching_never_proposes_is_created(world):
    levels, line, cand = world
    site(100.3, '000100', '000200', line, cl.DECISIONS_ACCEPT)

    proto = cl.forced_pair_prototypes([line], [cand])
    assert set(proto) == {('000100', '000200')}

    matches = cl.force_ledger_assignments(line, [], proto)

    assert len(matches) == 1
    forced = matches[0]
    assert forced is not cand                      # a candidate of its own
    assert forced.assigned_to is line
    assert forced in line.assigned_transitions
    # the predicted intensity travels with it, so the grading and the
    # intensity-pattern check see the same numbers as for any other candidate
    assert forced.calc_intensity == cand.calc_intensity
    assert forced.u_calc == cand.u_calc

    cl.apply_line_decisions(line)
    assert forced.accepted == 1
    assert forced.manual == cl.DECISIONS_ACCEPT


def test_a_pair_already_on_the_line_is_not_duplicated(world):
    levels, line, cand = world
    site(100.3, '000100', '000200', line, cl.DECISIONS_ACCEPT)
    proto = cl.forced_pair_prototypes([line], [cand])

    already = Transition(lower_level=cand.lower_level,
                         upper_level=cand.upper_level, assigned_to=line)
    line.assigned_transitions.append(already)

    matches = cl.force_ledger_assignments(line, [already], proto)
    assert matches == [already]
    assert len(line.assigned_transitions) == 1


def test_a_rejected_pair_is_never_created(world):
    levels, line, cand = world
    site(100.3, '000100', '000200', line, cl.DECISIONS_REJECT)

    assert cl.forced_pair_prototypes([line], [cand]) == {}
    assert cl.force_ledger_assignments(line, [], {}) == []
    assert line.assigned_transitions == []


# --------------------------------------------------------------------------
# refusing an order that cannot be carried out
# --------------------------------------------------------------------------
def test_a_row_naming_an_unknown_level_stops_the_run(world):
    levels, line, cand = world
    sites = [site(100.3, '000100', '000999', line, cl.DECISIONS_ACCEPT)]
    with pytest.raises(ValueError, match='000999 is not in the level list'):
        cl.check_forced_decisions(sites, levels)


def test_an_accepted_pair_of_the_same_parity_stops_the_run(world):
    levels, line, cand = world
    levels['000200'].parity = 'e'
    sites = [site(100.3, '000100', '000200', line, cl.DECISIONS_ACCEPT)]
    with pytest.raises(ValueError, match='opposite parities'):
        cl.check_forced_decisions(sites, levels)


def test_an_accepted_pair_breaking_the_J_rule_stops_the_run(world):
    levels, line, cand = world
    levels['000200'].J_val = 5.5
    sites = [site(100.3, '000100', '000200', line, cl.DECISIONS_ACCEPT)]
    with pytest.raises(ValueError, match=r'J\(upper\) - J\(lower\)'):
        cl.check_forced_decisions(sites, levels)


def test_the_same_pair_accepted_on_two_lines_stops_the_run(world):
    levels, line, cand = world
    other = SpectralLine(wavenumber=100.9, wn_uncertainty=0.01, intensity=5.0,
                         line_character='')
    sites = [site(100.3, '000100', '000200', line, cl.DECISIONS_ACCEPT),
             site(100.9, '000100', '000200', other, cl.DECISIONS_ACCEPT)]
    with pytest.raises(ValueError, match='cannot belong to two lines'):
        cl.check_forced_decisions(sites, levels)


def test_accepting_on_one_line_and_rejecting_on_another_is_allowed(world):
    """That is how an assignment is moved from one line to another."""
    levels, line, cand = world
    other = SpectralLine(wavenumber=100.9, wn_uncertainty=0.01, intensity=5.0,
                         line_character='')
    sites = [site(100.3, '000100', '000200', line, cl.DECISIONS_ACCEPT),
             site(100.9, '000100', '000200', other, cl.DECISIONS_REJECT)]
    cl.check_forced_decisions(sites, levels)      # no exception


def test_every_fault_is_reported_at_once(world):
    levels, line, cand = world
    other = SpectralLine(wavenumber=100.9, wn_uncertainty=0.01, intensity=5.0,
                         line_character='')
    sites = [site(100.3, '000100', '000999', line, cl.DECISIONS_ACCEPT),
             site(100.9, '000777', '000200', other, cl.DECISIONS_REJECT)]
    with pytest.raises(ValueError) as exc:
        cl.check_forced_decisions(sites, levels)
    assert '000999' in str(exc.value) and '000777' in str(exc.value)


def test_a_ritz_wavenumber_outside_the_range_stops_the_run(world, monkeypatch):
    levels, line, cand = world
    monkeypatch.setattr(cl, 'WN_MIN', 200.0)       # the run covers 200-300
    monkeypatch.setattr(cl, 'WN_MAX', 300.0)
    sites = [site(100.3, '000100', '000200', line, cl.DECISIONS_ACCEPT)]
    with pytest.raises(ValueError, match='outside the range'):
        cl.check_forced_decisions(sites, levels)
