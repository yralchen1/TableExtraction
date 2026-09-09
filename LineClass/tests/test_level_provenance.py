"""Tests for what the revised-energies comment column does to old status.

An identification taken from the published line list names its two levels in
the line workbook, which is never edited, so it goes on naming a level that
the identification work has since moved.  The comment column of
revised_level_energies.csv says what was done, and these tests cover the three
readings of it: a level re-positioned, whose published identifications lose
their old status; two levels whose measured positions were exchanged, whose
published identifications keep it under the other identifier; and a row that
says neither, which changes nothing.

Run from the LineClass directory:  python -m pytest tests -q
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import classify_lines as cl
from models import EnergyLevel, SpectralLine, Transition


def level(lid, energy, parity='e', j=2.5):
    return EnergyLevel(level_id=lid, energy=energy, parity=parity,
                       J_str=str(j), J_val=float(j))


def overrides(tmp_path, rows):
    """Write a revised-energies file of (level_id, E_input, comment) rows."""
    path = tmp_path / 'revised_level_energies.csv'
    text = ['level_id,E_input,comment']
    for lid, e, comment in rows:
        text.append('%s,%s,"%s"' % (lid, e, comment))
    path.write_text('\n'.join(text) + '\n', encoding='utf-8')
    return str(path)


@pytest.fixture(autouse=True)
def clean(monkeypatch):
    """Neither kind of move is in force unless the test puts one there."""
    monkeypatch.setattr(cl, 'LEGACY_MOVED', set())
    monkeypatch.setattr(cl, 'LEGACY_SWAPPED', {})


# ---------------------------------------------------------------------------
# Reading the comment
# ---------------------------------------------------------------------------
def test_repositioned_and_moved_are_read(tmp_path):
    path = overrides(tmp_path, [
        ('000565', 118967.3776, 're-positioned from 119225.52; weighted mean'),
        ('000424', 133078.09, 'repositioned from 133352.11 (dE = -274.02)'),
        ('000371', 113823.12, 'moved on the evidence of nine new lines'),
    ])
    moved, swapped = cl.read_level_provenance(path)
    assert moved == {'000565', '000424', '000371'}
    assert swapped == {}


def test_an_exchange_is_read_from_either_wording(tmp_path):
    path = overrides(tmp_path, [
        ('000483', 117686.46, 'exchanged with 000398 on 9/6/2026 (levels '
                              '000398 and 000483 were swapped)'),
        ('000447', 84409.903, 'exchanged with 000298'),
    ])
    moved, swapped = cl.read_level_provenance(path)
    assert moved == set()
    assert swapped == {'000483': '000398', '000447': '000298'}


def test_a_plain_revision_is_neither(tmp_path):
    path = overrides(tmp_path, [
        ('000600', 136638.4947, 'from iter22/LOPT_output_levels.txt'),
    ])
    assert cl.read_level_provenance(path) == (set(), {})


def test_an_exchange_is_an_exchange_even_if_it_also_says_moved(tmp_path):
    path = overrides(tmp_path, [
        ('000483', 117686.46, 'moved: levels 000398 and 000483 were swapped'),
    ])
    moved, swapped = cl.read_level_provenance(path)
    assert moved == set()
    assert swapped == {'000483': '000398'}


# ---------------------------------------------------------------------------
# Adopting it
# ---------------------------------------------------------------------------
def test_an_exchange_recorded_once_is_recorded_from_both_sides(tmp_path):
    path = overrides(tmp_path, [
        ('000483', 117686.46, 'levels 000398 and 000483 were swapped'),
    ])
    levels = {lid: level(lid, 0.0) for lid in ('000398', '000483')}
    cl.set_level_provenance(levels, path)
    assert cl.LEGACY_SWAPPED == {'000483': '000398', '000398': '000483'}


def test_an_unknown_partner_stops_the_run(tmp_path):
    path = overrides(tmp_path, [
        ('000483', 117686.46, 'exchanged with 000999'),
    ])
    with pytest.raises(ValueError, match='not in the level list'):
        cl.set_level_provenance({'000483': level('000483', 0.0)}, path)


def test_a_half_finished_exchange_stops_the_run(tmp_path):
    path = overrides(tmp_path, [
        ('000483', 117686.46, 'exchanged with 000398'),
        ('000398', 117592.99, 'exchanged with 000447'),
    ])
    levels = {lid: level(lid, 0.0)
              for lid in ('000398', '000447', '000483')}
    with pytest.raises(ValueError, match='does not agree'):
        cl.set_level_provenance(levels, path)


# ---------------------------------------------------------------------------
# What it does to a published identification
# ---------------------------------------------------------------------------
def test_an_untouched_pair_keeps_its_identification():
    assert cl.legacy_identification('000187', '000298') == ('000187',
                                                            '000298')


def test_a_repositioned_level_has_no_published_identification(monkeypatch):
    monkeypatch.setattr(cl, 'LEGACY_MOVED', {'000424'})
    assert cl.legacy_identification('000164', '000424') is None


def test_an_exchange_carries_the_identification_over(monkeypatch):
    monkeypatch.setattr(cl, 'LEGACY_SWAPPED',
                        {'000298': '000447', '000447': '000298'})
    assert cl.legacy_identification('000187', '000298') == ('000187',
                                                            '000447')


def line_with_legacy(pairs, energies=None):
    """An observed line carrying one published identification per pair.

    The candidate is seeded on the line exactly as read_observed_lines seeds
    it: on assigned_transitions, where the weeding sees it, as well as on
    original_assignments.  `energies` gives the upper levels a Ritz
    wavenumber of the caller's choosing; the default puts every one of them
    on the line.
    """
    line = SpectralLine(wavenumber=100.0, wn_uncertainty=0.01,
                        intensity=50.0, line_character='')
    for lo, up in pairs:
        tr = Transition(lower_level=level(lo, 0.0),
                        upper_level=level(up, (energies or {}).get(up, 100.0)),
                        assigned_to=line)
        line.original_assignments.append(tr)
        line.assigned_transitions.append(tr)
    return line


def pairs_on(line):
    """The level pairs the line holds a candidate for."""
    return {(t.lower_level.level_id, t.upper_level.level_id)
            for t in line.assigned_transitions}


def test_retagging_keeps_carries_and_withdraws(monkeypatch, capsys):
    monkeypatch.setattr(cl, 'LEGACY_MOVED', {'000424'})
    monkeypatch.setattr(cl, 'LEGACY_SWAPPED',
                        {'000298': '000447', '000447': '000298'})
    line = line_with_legacy([('000100', '000200'),    # untouched
                             ('000187', '000298'),    # exchanged
                             ('000164', '000424')])   # re-positioned
    cl.retag_legacy_identifications([line])
    assert line.legacy_keys == {('000100', '000200'), ('000187', '000447')}
    said = capsys.readouterr().out
    assert '1 stand as written' in said
    assert '1 carried over' in said
    assert '1 withdrawn' in said


def test_old_status_follows_the_retagged_pair(monkeypatch):
    """The pair the line now names is old; the pair it used to name is not."""
    monkeypatch.setattr(cl, 'LEGACY_SWAPPED',
                        {'000298': '000447', '000447': '000298'})
    line = line_with_legacy([('000187', '000298')])
    cl.retag_legacy_identifications([line])
    stale = Transition(lower_level=level('000187', 0.0),
                       upper_level=level('000298', 100.0),
                       calc_intensity=42.0, u_calc=0.7, assigned_to=line)
    carried = Transition(lower_level=level('000187', 0.0),
                         upper_level=level('000447', 100.0),
                         calc_intensity=42.0, u_calc=0.7, assigned_to=line)
    cl.assign_grades(line, [stale, carried])
    assert stale.new == 1
    assert carried.new == 0


# ---------------------------------------------------------------------------
# What it does to the candidate the workbook seeded
# ---------------------------------------------------------------------------
def test_the_stale_pair_is_no_longer_a_candidate(monkeypatch):
    """The exchanged-away pair leaves the line instead of standing there.

    Its identifier now sits far from this line - that is what the exchange
    said - so a seeded candidate for it is a row the analyst has to read and
    dismiss, with a Ritz mismatch of the whole distance the level moved.
    """
    monkeypatch.setattr(cl, 'LEGACY_SWAPPED',
                        {'000298': '000447', '000447': '000298'})
    line = line_with_legacy([('000187', '000298')],
                            energies={'000298': 192.9})
    cl.retag_legacy_identifications([line])
    assert pairs_on(line) == set()
    assert line.original_assignments == []
    assert line.legacy_keys == {('000187', '000447')}


def test_a_repositioned_level_takes_its_candidate_with_it(monkeypatch):
    monkeypatch.setattr(cl, 'LEGACY_MOVED', {'000424'})
    line = line_with_legacy([('000164', '000424'), ('000100', '000200')])
    cl.retag_legacy_identifications([line])
    assert pairs_on(line) == {('000100', '000200')}


def test_the_exchange_re_seeds_under_the_new_identifier(monkeypatch):
    """Given the level list, the identification is put back, re-keyed."""
    monkeypatch.setattr(cl, 'LEGACY_SWAPPED',
                        {'000298': '000447', '000447': '000298'})
    line = line_with_legacy([('000187', '000298')])
    levels = {'000187': level('000187', 0.0),
              '000447': level('000447', 100.0)}
    calc = {('000187', '000447'): {'calc_intensity': 42.0, 'u_calc': 0.7,
                                   'assigned_to': None}}
    cl.retag_legacy_identifications([line], levels, calc)
    assert pairs_on(line) == {('000187', '000447')}
    seeded = line.original_assignments[0]
    assert seeded.calc_intensity == 42.0
    assert calc[('000187', '000447')]['assigned_to'] is line


def test_un_seeding_releases_the_calculated_transition(monkeypatch):
    """The record of which line holds a pair goes with the candidate."""
    monkeypatch.setattr(cl, 'LEGACY_MOVED', {'000424'})
    line = line_with_legacy([('000164', '000424')])
    calc = {('000164', '000424'): {'calc_intensity': 1.0, 'u_calc': 0.5,
                                   'assigned_to': line}}
    cl.retag_legacy_identifications([line], {}, calc)
    assert calc[('000164', '000424')]['assigned_to'] is None


def test_without_the_step_the_seeded_transitions_are_the_old_ones():
    """A run that never retagged keeps the behaviour it had before."""
    line = line_with_legacy([('000187', '000298')])
    assert line.legacy_keys is None
    stale = line.original_assignments[0]
    other = Transition(lower_level=level('000187', 0.0),
                       upper_level=level('000447', 100.0),
                       calc_intensity=42.0, u_calc=0.7, assigned_to=line)
    cl.assign_grades(line, [stale, other])
    assert stale.new == 0
    assert other.new == 1
