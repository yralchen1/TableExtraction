"""Tests for `move_level`.

Run from the LineClass directory:  python -m pytest tests -q

`move_level.py` runs LOPT and three other programs, so the tests here do not
run the chain.  They fix the decisions the chain is asked to carry out - which
of an already assigned level's lines survive a move and which are released,
and with what reason - and the file handling that has to be right before any
of it is safe:

  * the four release tests, each in isolation and in the order they are applied:
    a line too far from the Ritz wavenumber the new position gives, a line too
    strong for the transition to account for, a component too small a share of
    its blend to matter, and a pair with no calculated transition at all;
  * a line that fails none of them survives, and its record, its mark and its
    ledger row are left alone;
  * the level's energy is written to `new_levels.txt` when that is where it
    lives and to `revised_level_energies.csv` when it is not, and to exactly
    one of the two;
  * the ledger rows a move takes out are removed from `line_decisions.csv` and
    kept in `line_decisions_removed.csv` with the reason;
  * `order_pair` puts the two levels of a transition the right way round,
    whichever way the move leaves their energies.
"""
import csv
import io
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (ROOT, os.path.dirname(ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

import move_level as ML                  # noqa: E402
import insert_new_level as INL           # noqa: E402
from models import EnergyLevel, SpectralLine   # noqa: E402


# ---------------------------------------------------------------------------
# A stand-in for the command line, carrying only what the release tests read
# ---------------------------------------------------------------------------
class Args(object):
    def __init__(self, **kw):
        self.window = INL.DEF_WINDOW
        self.strong_sigma = INL.DEF_STRONG_SIGMA
        self.release_sigma = ML.DEF_RELEASE_SIGMA
        self.masked_share = ML.DEF_MASKED_SHARE
        for k, v in kw.items():
            setattr(self, k, v)


LOW = '059003.000100'
UPP = '059003.000200'
OTHER = '059003.000300'


def a_model(e_low=1000.0, e_upp=21000.0, wn=20000.0, unc=0.05,
            i_obs=100.0, i_calc=100.0, u_calc=0.6):
    """One transition, one observed line, and the indexes the tests need."""
    levels = {
        LOW: EnergyLevel(level_id=LOW, energy=e_low, parity='e', J_str='2',
                         J_val=2.0),
        UPP: EnergyLevel(level_id=UPP, energy=e_upp, parity='o', J_str='3',
                         J_val=3.0),
        OTHER: EnergyLevel(level_id=OTHER, energy=e_upp + 5.0, parity='o',
                           J_str='3', J_val=3.0),
    }
    calc_index = {(LOW, UPP): {'calc_intensity': i_calc, 'u_calc': u_calc},
                  (LOW, OTHER): {'calc_intensity': 500.0, 'u_calc': u_calc}}
    line = SpectralLine(wavenumber=wn, wn_uncertainty=unc, intensity=i_obs,
                        line_character='')
    wns, ordered = ML.line_index([line])
    return levels, calc_index, wns, ordered


def one_lopt_row(wn, low_id=LOW, upp_id=UPP, flag=''):
    return {'wn': wn, 'unc': 0.05, 'intens': 100.0, 'low_id': low_id,
            'upp_id': upp_id, 'flag': flag, 'weight': 1.0, 'raw': None}


# ---------------------------------------------------------------------------
# D. The release tests
# ---------------------------------------------------------------------------
def test_a_line_that_fails_nothing_survives():
    levels, calc_index, wns, ordered = a_model()
    why = ML.release_reason(20000.0, LOW, UPP, levels, calc_index, wns,
                            ordered, [one_lopt_row(20000.0)],
                            {(LOW, UPP): 100.0}, Args())
    assert why == ''


def test_too_far_from_ritz():
    # The level has moved 1 cm^-1, which is 20 times the line's 0.05 cm^-1
    # uncertainty: the line is not this transition's any more.
    levels, calc_index, wns, ordered = a_model(e_upp=21001.0)
    why = ML.release_reason(20000.0, LOW, UPP, levels, calc_index, wns,
                            ordered, [one_lopt_row(20000.0)],
                            {(LOW, UPP): 100.0}, Args())
    assert why.startswith(ML.REASON_FAR)
    assert '20.0 sigma' in why


def test_a_small_move_inside_the_uncertainty_is_not_far():
    levels, calc_index, wns, ordered = a_model(e_upp=21000.1)
    why = ML.release_reason(20000.0, LOW, UPP, levels, calc_index, wns,
                            ordered, [one_lopt_row(20000.0)],
                            {(LOW, UPP): 100.0}, Args())
    assert why == ''


def test_too_weak_to_explain_iobs():
    # The line is 1000 times stronger than the transition can make it.
    levels, calc_index, wns, ordered = a_model(i_obs=100000.0, i_calc=100.0)
    why = ML.release_reason(20000.0, LOW, UPP, levels, calc_index, wns,
                            ordered, [one_lopt_row(20000.0)],
                            {(LOW, UPP): 100.0}, Args())
    assert why.startswith(ML.REASON_WEAK)


def test_far_from_ritz_is_tested_before_intensity():
    """A line that fails both tests is released for the wavenumber.

    A line nowhere near the new Ritz wavenumber is not this transition's
    whatever its intensity, so that is the reason to record.
    """
    levels, calc_index, wns, ordered = a_model(e_upp=21001.0, i_obs=100000.0)
    why = ML.release_reason(20000.0, LOW, UPP, levels, calc_index, wns,
                            ordered, [one_lopt_row(20000.0)],
                            {(LOW, UPP): 100.0}, Args())
    assert why.startswith(ML.REASON_FAR)


def test_too_weak_to_contribute_to_blend():
    # Two components on one line: this one is 1% of the predicted intensity.
    levels, calc_index, wns, ordered = a_model()
    rows = [one_lopt_row(20000.0, LOW, UPP),
            one_lopt_row(20000.0, LOW, OTHER)]
    calc_of = {(LOW, UPP): 1.0, (LOW, OTHER): 99.0}
    why = ML.release_reason(20000.0, LOW, UPP, levels, calc_index, wns,
                            ordered, rows, calc_of, Args())
    assert why.startswith('Too weak to contribute to blend')
    assert '1.0%' in why


def test_a_transition_alone_on_its_line_is_never_too_weak_for_it():
    """The blend share is not a test a single component can fail."""
    levels, calc_index, wns, ordered = a_model()
    why = ML.release_reason(20000.0, LOW, UPP, levels, calc_index, wns,
                            ordered, [one_lopt_row(20000.0)],
                            {(LOW, UPP): 1.0}, Args())
    assert why == ''


def test_a_flagged_record_is_not_a_component_of_the_blend():
    """A `P` record is a candidate LOPT is shown and does not fit."""
    levels, calc_index, wns, ordered = a_model()
    rows = [one_lopt_row(20000.0, LOW, UPP),
            one_lopt_row(20000.0, LOW, OTHER, flag='P')]
    calc_of = {(LOW, UPP): 1.0, (LOW, OTHER): 99.0}
    why = ML.release_reason(20000.0, LOW, UPP, levels, calc_index, wns,
                            ordered, rows, calc_of, Args())
    assert why == ''


def test_no_calculated_transition_at_the_new_position():
    levels, calc_index, wns, ordered = a_model()
    why = ML.release_reason(20000.0, LOW, OTHER + '9', levels, calc_index,
                            wns, ordered, [], {}, Args())
    assert why == ML.REASON_NONE


def test_no_observed_line_at_that_wavenumber():
    levels, calc_index, wns, ordered = a_model()
    why = ML.release_reason(12345.678, LOW, UPP, levels, calc_index, wns,
                            ordered, [], {}, Args())
    assert why == ML.REASON_NONE


def test_release_sigma_is_the_knob():
    levels, calc_index, wns, ordered = a_model(e_upp=21000.2)   # 4 sigma
    assert ML.release_reason(20000.0, LOW, UPP, levels, calc_index, wns,
                             ordered, [one_lopt_row(20000.0)],
                             {(LOW, UPP): 100.0},
                             Args(release_sigma=3.0)) != ''
    assert ML.release_reason(20000.0, LOW, UPP, levels, calc_index, wns,
                             ordered, [one_lopt_row(20000.0)],
                             {(LOW, UPP): 100.0},
                             Args(release_sigma=5.0)) == ''


# ---------------------------------------------------------------------------
# Small pieces
# ---------------------------------------------------------------------------
def test_order_pair_follows_the_energies():
    levels, _c, _w, _o = a_model()
    assert ML.order_pair(UPP, LOW, levels) == (LOW, UPP)
    assert ML.order_pair(LOW, UPP, levels) == (LOW, UPP)


def test_order_pair_leaves_an_unknown_level_alone():
    levels, _c, _w, _o = a_model()
    assert ML.order_pair(LOW, 'nobody', levels) == (LOW, 'nobody')


def test_line_at_takes_the_nearest_within_the_tolerance():
    a = SpectralLine(wavenumber=20000.000, wn_uncertainty=0.05,
                     intensity=1.0, line_character='')
    b = SpectralLine(wavenumber=20000.004, wn_uncertainty=0.05,
                     intensity=1.0, line_character='')
    wns, ordered = ML.line_index([a, b])
    assert ML.line_at(wns, ordered, 20000.0045) is b
    assert ML.line_at(wns, ordered, 20000.0) is a
    assert ML.line_at(wns, ordered, 20001.0) is None


def test_holding_says_where_it_is_recorded():
    h = ML.Holding(20000.0, LOW, UPP)
    assert h.where() == '-'
    h.in_lopt = True
    h.partner_row = 742
    h.ledger.append({'decision': 'accept'})
    assert h.where() == 'LOPT+IDEN2+ledger x1'


# ---------------------------------------------------------------------------
# B/F. Where the energy is written
# ---------------------------------------------------------------------------
@pytest.fixture
def sandbox(tmp_path, monkeypatch):
    """`new_levels.txt` and `revised_level_energies.csv` in a temporary copy."""
    new_levels = tmp_path / 'new_levels.txt'
    new_levels.write_text(
        'level_id\tE\tJ\tparity\tiden2_row\tcowan_lid\tcomment\n'
        '059003.000624\t116327.560\t5/2\to\t742\t294\tfound with --unknown\n',
        encoding='utf-8')
    overrides = tmp_path / 'revised_level_energies.csv'
    overrides.write_text(
        'level_id,E_input,comment\n'
        '059003.000565,118967.3776,"re-positioned from 119225.52"\n',
        encoding='utf-8')
    monkeypatch.setattr(ML, 'NEW_LEVELS', str(new_levels))
    monkeypatch.setattr(INL, 'NEW_LEVELS', str(new_levels))
    monkeypatch.setattr(ML, 'OVERRIDES', str(overrides))
    return tmp_path


def test_a_new_level_keeps_its_energy_in_new_levels(sandbox):
    assert ML.where_the_energy_lives('059003.000624') == 'new_levels'


def test_a_published_level_is_moved_by_an_override(sandbox):
    assert ML.where_the_energy_lives('059003.000565') == 'overrides'
    assert ML.where_the_energy_lives('059003.000001') == 'overrides'


def test_setting_a_new_level_energy_writes_only_new_levels(sandbox):
    before = (sandbox / 'revised_level_energies.csv').read_text(
        encoding='utf-8')
    ML.set_energy('059003.000624', 116330.1234, 'new_levels', 'moved',
                  lambda t: None)
    _fields, rows = INL.read_new_levels(str(sandbox / 'new_levels.txt'))
    assert rows[0]['E'] == '116330.1234'
    assert rows[0]['comment'] == 'moved'
    assert (sandbox / 'revised_level_energies.csv').read_text(
        encoding='utf-8') == before


def test_setting_a_published_level_energy_updates_its_override_row(sandbox):
    before = (sandbox / 'new_levels.txt').read_text(encoding='utf-8')
    ML.set_energy('059003.000565', 118970.0, 'overrides', 'moved again',
                  lambda t: None)
    fields, rows = ML.read_overrides(str(sandbox / 'revised_level_energies.csv'))
    assert len(rows) == 1
    assert rows[0]['E_input'] == '118970.0000'
    assert rows[0]['comment'] == 'moved again'
    assert (sandbox / 'new_levels.txt').read_text(encoding='utf-8') == before


def test_a_published_level_with_no_override_row_yet_gets_one(sandbox):
    ML.set_energy('059003.000424', 133078.09, 'overrides', 'first move',
                  lambda t: None)
    _fields, rows = ML.read_overrides(
        str(sandbox / 'revised_level_energies.csv'))
    assert [r['level_id'] for r in rows] == ['059003.000565', '059003.000424']
    assert rows[1]['E_input'] == '133078.0900'


def test_the_overrides_file_keeps_lf_endings(sandbox):
    ML.set_energy('059003.000424', 133078.09, 'overrides', 'first move',
                  lambda t: None)
    raw = (sandbox / 'revised_level_energies.csv').read_bytes()
    assert b'\r' not in raw


# ---------------------------------------------------------------------------
# E. The ledger rows a move takes out
# ---------------------------------------------------------------------------
def test_removed_ledger_rows_are_kept_with_their_reason(tmp_path):
    removed = tmp_path / 'line_decisions_removed.csv'
    rec = {'wn_obs': '20000.0000', 'low_id': LOW, 'upp_id': UPP,
           'decision': 'accept', 'date': '9/3/2026', 'reason': 'hand-made'}
    ML.archive_ledger_rows(str(removed), [(rec, ML.REASON_FAR)],
                           UPP, '9/15/2026')
    with io.open(str(removed), encoding='utf-8-sig', newline='') as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 1
    assert rows[0]['removed_on'] == '9/15/2026'
    assert rows[0]['removed_for'] == UPP
    assert rows[0]['removed_because'] == ML.REASON_FAR
    assert rows[0]['reason'] == 'hand-made'          # the original, kept
    assert rows[0]['wn_obs'] == '20000.0000'


def test_archiving_appends_and_writes_one_header(tmp_path):
    removed = tmp_path / 'line_decisions_removed.csv'
    rec = {'wn_obs': '20000.0000', 'low_id': LOW, 'upp_id': UPP,
           'decision': 'accept', 'date': '9/3/2026', 'reason': 'hand-made'}
    ML.archive_ledger_rows(str(removed), [(rec, 'one')], UPP, '9/15/2026')
    ML.archive_ledger_rows(str(removed), [(rec, 'two')], UPP, '9/15/2026')
    text = removed.read_text(encoding='utf-8')
    assert text.count('removed_on') == 1
    with io.open(str(removed), encoding='utf-8-sig', newline='') as fh:
        assert len(list(csv.DictReader(fh))) == 2


def test_archiving_nothing_creates_no_file(tmp_path):
    removed = tmp_path / 'line_decisions_removed.csv'
    ML.archive_ledger_rows(str(removed), [], UPP, '9/15/2026')
    assert not removed.exists()


def test_the_ledger_is_rewritten_without_the_removed_rows(tmp_path):
    ledger = tmp_path / 'line_decisions.csv'
    ledger.write_text(
        'wn_obs,low_id,upp_id,decision,date,reason\n'
        '20000.0000,%s,%s,accept,9/3/2026,kept\n'
        '30000.0000,%s,%s,accept,9/3/2026,goes\n' % (LOW, OTHER, LOW, UPP),
        encoding='utf-8')
    fields, rows = ML.read_ledger_rows(str(ledger))
    keep = [r for r in rows if r['upp_id'] != UPP]
    ML.write_ledger(str(ledger), fields, keep)
    _f, again = ML.read_ledger_rows(str(ledger))
    assert [r['reason'] for r in again] == ['kept']
    assert b'\r' not in ledger.read_bytes()


# ---------------------------------------------------------------------------
# The two tools stay in step
# ---------------------------------------------------------------------------
def test_move_level_backs_up_everything_insert_new_level_does():
    """A move writes every file an insertion writes, and two more.

    If insert_new_level.py ever learns to write a file that move_level.py does
    not back up, a failed move would leave that one file changed and every
    other put back - which is the silent disagreement both tools exist to
    prevent.
    """
    missing = [p for p in INL.WRITABLE
               if p != INL.LOGFILE and p not in ML.WRITABLE]
    assert missing == []
    assert ML.OVERRIDES in ML.WRITABLE
    assert ML.REMOVED in ML.WRITABLE
    assert INL.LOGFILE not in ML.WRITABLE


def test_the_two_tools_use_their_own_backup_directories():
    assert ML.BACKUP_DIR != INL.BACKUP_DIR


def test_preflight_takes_a_backup_directory(tmp_path):
    """The generalization move_level.py needs of insert_new_level.preflight."""
    src = tmp_path / 'a.txt'
    src.write_text('x', encoding='utf-8')
    where = tmp_path / 'backup'
    saved = INL.preflight([str(src)], lambda t=None: None, str(where))
    assert saved == {str(src): str(where / 'a.txt')}
    src.write_text('y', encoding='utf-8')
    INL.restore(saved, lambda t=None: None)
    assert src.read_text(encoding='utf-8') == 'x'


def test_the_command_line_insists_on_a_level_and_a_destination():
    with pytest.raises(SystemExit):
        ML.parse_args(['--energy', '100.0'])          # no level
    with pytest.raises(SystemExit):
        ML.parse_args(['059003.000565'])              # nowhere to go
    with pytest.raises(SystemExit):
        ML.parse_args(['059003.000565', '--energy', '1.0', '--from-lopt'])
    args = ML.parse_args(['059003.000565', '--energy', '118967.3776'])
    assert args.level_id == '059003.000565'
    assert args.energy == pytest.approx(118967.3776)
    assert args.propose is True
    assert ML.parse_args(['--undo']).undo is True
