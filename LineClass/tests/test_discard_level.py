"""Tests for `discard_level`.

Run from the LineClass directory:  python -m pytest tests -q

`discard_level.py` runs LOPT and three other programs, so the tests here do not
run the chain.  They fix the decisions the chain is asked to carry out, and the
file handling that has to be right before any of it is safe:

  * the readers: a discarded level leaves both `levels_dict` and `levels_list`;
    a level named in the ledger that is in no level list raises, and so does
    one named in both the ledger and the revised-energies file;
  * `legacy_identification` withdraws a published identification naming a
    discarded level, and carries none of it over;
  * the identifier-map row is taken out and its number kept in the ledger,
    which is what `--undo` restores from;
  * the unfound shape on disk: the adopted energy becomes the calculated one,
    the O-C becomes zero, the star goes out, the uncertainty is NOT touched -
    it is sync_IDEN2.py's to write - and every copy of the energy `trans.dat`
    keeps in other levels' blocks comes with it;
  * the per-configuration search window narrows when a badly placed level
    stops counting as found, and the MIN_CFG_LEVELS threshold is reported:
    discarding one of three found levels leaves two, too few to measure
    anything, and those unfound levels keep the value they had;
  * the evidence guard: a level with positive ln_R, a level whose alternate
    the audit calls firm, and a level with no audit row at all are all refused
    without --force.
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

import classify_lines as CL              # noqa: E402
import discard_level as DL               # noqa: E402
import insert_new_level as INL           # noqa: E402
import swap_line_assignments_IDEN as IDEN  # noqa: E402
from models import EnergyLevel           # noqa: E402


LOW = '059003.000100'
UPP = '059003.000200'
GONE = '059003.000300'


class Args(object):
    """A stand-in for the command line, carrying only what a test reads."""

    def __init__(self, **kw):
        self.audit = ''
        self.force = False
        for k, v in kw.items():
            setattr(self, k, v)


def write_csv(path, fields, rows):
    with io.open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fields, lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow(r)


def a_level_list():
    levels = {}
    order = []
    for lid, e in ((LOW, 1000.0), (UPP, 21000.0), (GONE, 30000.0)):
        lev = EnergyLevel(level_id=lid, energy=e, parity='e', J_str='2',
                          J_val=2.0)
        levels[lid] = lev
        order.append(lev)
    return levels, order


# ---------------------------------------------------------------------------
# The reader: classify_lines.drop_discarded_levels
# ---------------------------------------------------------------------------
def test_a_discarded_level_leaves_both_the_dict_and_the_list(tmp_path,
                                                             monkeypatch):
    path = str(tmp_path / 'discarded_levels.csv')
    write_csv(path, DL.DISCARD_COLUMNS,
              [{'level_id': GONE, 'iden2_row': '812', 'date': '2026-09-20',
                'ln_R': '-5.492', 'reason': 'no support'}])
    monkeypatch.setattr(CL, 'LEVEL_OVERRIDES', '')
    levels, order = a_level_list()
    n = CL.drop_discarded_levels(levels, order, path)
    assert n == 1
    assert GONE not in levels
    assert [lv.level_id for lv in order] == [LOW, UPP]
    assert CL.DISCARDED == {GONE}


def test_a_level_named_in_no_list_raises(tmp_path, monkeypatch):
    """A name that matches nothing is a typo, not a silent no-op."""
    path = str(tmp_path / 'discarded_levels.csv')
    write_csv(path, DL.DISCARD_COLUMNS,
              [{'level_id': '059003.999999', 'iden2_row': '1',
                'date': '2026-09-20', 'ln_R': '', 'reason': 'x'}])
    monkeypatch.setattr(CL, 'LEVEL_OVERRIDES', '')
    levels, order = a_level_list()
    with pytest.raises(ValueError, match='not in the level list'):
        CL.drop_discarded_levels(levels, order, path)


def test_discarded_and_revised_at_once_raises(tmp_path, monkeypatch):
    """One file says where the level is, the other says nobody knows."""
    path = str(tmp_path / 'discarded_levels.csv')
    write_csv(path, DL.DISCARD_COLUMNS,
              [{'level_id': GONE, 'iden2_row': '812', 'date': '2026-09-20',
                'ln_R': '-5.5', 'reason': 'x'}])
    overrides = str(tmp_path / 'revised_level_energies.csv')
    write_csv(overrides, ['level_id', 'E_input', 'comment'],
              [{'level_id': GONE, 'E_input': '30001.0', 'comment': ''}])
    monkeypatch.setattr(CL, 'LEVEL_OVERRIDES', overrides)
    levels, order = a_level_list()
    with pytest.raises(ValueError, match='in both'):
        CL.drop_discarded_levels(levels, order, path)


def test_a_ledger_without_a_level_id_column_raises(tmp_path):
    path = str(tmp_path / 'discarded_levels.csv')
    with io.open(path, 'w', encoding='utf-8', newline='') as fh:
        fh.write('level\treason\n812\tno support\n')
    with pytest.raises(ValueError, match='level_id'):
        CL.read_discarded_records(path)


def test_a_published_identification_of_a_discarded_level_is_withdrawn(
        monkeypatch):
    """Withdrawn, and carried nowhere: there is no position to carry it to."""
    monkeypatch.setattr(CL, 'LEGACY_MOVED', set())
    monkeypatch.setattr(CL, 'LEGACY_SWAPPED', {})
    monkeypatch.setattr(CL, 'DISCARDED', {GONE})
    assert CL.legacy_identification(LOW, GONE) is None
    assert CL.legacy_identification(GONE, UPP) is None
    assert CL.legacy_identification(LOW, UPP) == (LOW, UPP)


# ---------------------------------------------------------------------------
# F. The identifier map
# ---------------------------------------------------------------------------
def test_the_map_row_is_removed_and_its_number_returned(tmp_path):
    """The number is what `iden2_row` keeps, and what --undo restores from."""
    path = str(tmp_path / 'IDEN_level_ids.txt')
    with io.open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('level_id\tIDEN_id\n%s\t100\n%s\t812\n%s\t900\n'
                 % (LOW, GONE, UPP))
    said = []
    row = DL.remove_id_map_row(path, GONE, said.append)
    assert row == 812
    left = INL.read_id_map_rows(path)
    assert GONE not in left
    assert left == {LOW: 100, UPP: 900}


def test_removing_a_map_row_that_is_not_there_says_so(tmp_path):
    path = str(tmp_path / 'IDEN_level_ids.txt')
    with io.open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('level_id\tIDEN_id\n%s\t100\n' % LOW)
    said = []
    assert DL.remove_id_map_row(path, GONE, said.append) == 0
    assert 'no row for' in said[0]


def test_the_ledger_round_trips(tmp_path):
    path = str(tmp_path / 'discarded_levels.csv')
    fields, rows = DL.read_discarded(path)
    rows.append({'level_id': GONE, 'iden2_row': '812', 'date': '2026-09-20',
                 'ln_R': '-5.492', 'reason': 'no support, 2 lines'})
    DL.write_discarded(path, fields, rows)
    back = DL.discarded_row(path, GONE)
    assert back['iden2_row'] == '812'
    assert back['reason'] == 'no support, 2 lines'
    assert DL.discarded_row(path, LOW) is None
    with io.open(path, 'rb') as fh:
        assert b'\r\n' not in fh.read()


# ---------------------------------------------------------------------------
# I. The unfound shape on disk
# ---------------------------------------------------------------------------
def enlev_record(index, e_calc, unc, e_obs, known, j, label):
    """One row of enlev.dat, at the fixed widths the file is read by."""
    rec = ('%4d%12.3f%10.3f%12.3f%s%12.3f%5s /%s/'
           % (index, e_calc, unc, e_obs,
              IDEN.STAR_ON if known else IDEN.STAR_OFF,
              e_obs - e_calc, j, label))
    assert len(rec) > IDEN.EN_LABEL
    return rec


def an_enlev(tmp_path, rows):
    path = str(tmp_path / 'enlev.dat')
    with io.open(path, 'w', encoding='latin-1', newline='\n') as fh:
        for r in rows:
            fh.write(r + '\n')
    return path


def a_trans(tmp_path, records):
    path = str(tmp_path / 'trans.dat')
    with io.open(path, 'w', encoding='latin-1', newline='\n') as fh:
        for r in records:
            fh.write(r + '\n')
    return path


def header_record(index, e_obs, known, j='5.5', label='f26d _1G2H'):
    rec = ('$%4d    J=%4s  %-11.11s%13.3f%s '
           % (index, j, label, e_obs,
              IDEN.STAR_ON if known else IDEN.STAR_OFF))
    return rec


def transition_record(partner, e_partner, partner_known, rwn):
    rec = '+%4d%5d%12.3f%s%13.3f' % (partner, -10, e_partner,
                                     IDEN.STAR_ON if partner_known
                                     else IDEN.STAR_OFF, rwn)
    return rec + IDEN.BLANK_OBS


def test_the_unfound_shape(tmp_path, monkeypatch):
    """E_obs becomes E_calc, O-C becomes zero, the star goes out - and the
    uncertainty is left exactly as it was, for sync_IDEN2.py to rewrite."""
    enlev = an_enlev(tmp_path, [
        enlev_record(812, 110903.800, 0.077, 111039.696, True, '5.5',
                     ' f26d _1G2H'),
        enlev_record(700, 100000.000, 0.050, 100010.000, True, '4.5',
                     ' f26d _1G2G'),
    ])
    trans = a_trans(tmp_path, [
        header_record(812, 111039.696, True),
        transition_record(700, 100010.000, True, 11029.696),
        header_record(700, 100010.000, True, j='4.5'),
        transition_record(812, 111039.696, True, 11029.696),
    ])
    monkeypatch.setattr(DL, 'ENLEV', enlev)
    monkeypatch.setattr(DL, 'TRANS', trans)
    said = []
    e_before, e_calc = DL.unfind_in_iden2(812, said.append)
    assert e_before == pytest.approx(111039.696)
    assert e_calc == pytest.approx(110903.800)

    back = IDEN.Enlev(enlev)
    assert back.known(812) is False
    assert back.e_obs(812) == pytest.approx(110903.800)
    assert IDEN.number(back.record(812), IDEN.EN_OMC) == pytest.approx(0.0)
    # The uncertainty is a prediction on an unfound row and a measurement on a
    # found one, and what belongs there is the configuration's rms; this tool
    # does not guess at it.
    assert IDEN.number(back.record(812), IDEN.EN_UNC) == pytest.approx(0.077)
    # The level next door is untouched.
    assert back.known(700) is True
    assert back.e_obs(700) == pytest.approx(100010.000)


def test_every_copy_of_the_energy_comes_with_it(tmp_path, monkeypatch):
    """check_sync reports an ERROR when trans.dat no longer follows from
    enlev.dat, so the copies and the predicted wavenumbers are brought along."""
    enlev = an_enlev(tmp_path, [
        enlev_record(812, 110903.800, 0.077, 111039.696, True, '5.5',
                     ' f26d _1G2H'),
        enlev_record(700, 100000.000, 0.050, 100010.000, True, '4.5',
                     ' f26d _1G2G'),
    ])
    trans = a_trans(tmp_path, [
        header_record(812, 111039.696, True),
        transition_record(700, 100010.000, True, 11029.696),
        header_record(700, 100010.000, True, j='4.5'),
        transition_record(812, 111039.696, True, 11029.696),
    ])
    monkeypatch.setattr(DL, 'ENLEV', enlev)
    monkeypatch.setattr(DL, 'TRANS', trans)
    DL.unfind_in_iden2(812, lambda *_a: None)

    back = IDEN.Trans(trans)
    head = back.records[back.header_of[812]]
    assert '*' not in IDEN.field(head, IDEN.TH_STAR)
    assert IDEN.number(head, IDEN.TH_EOBS) == pytest.approx(110903.800)
    # 700's block holds a copy of 812's energy and flag; both follow it.
    rec = back.records[back.row_of[(700, 812)]]
    assert IDEN.number(rec, IDEN.TR_EPART) == pytest.approx(110903.800)
    assert '*' not in IDEN.field(rec, IDEN.TR_STAR)
    assert IDEN.number(rec, IDEN.TR_WN) == pytest.approx(10893.800, abs=0.002)
    # And so does the predicted wavenumber in 812's own block.
    own = back.records[back.row_of[(812, 700)]]
    assert IDEN.number(own, IDEN.TR_WN) == pytest.approx(10893.800, abs=0.002)


# ---------------------------------------------------------------------------
# The per-configuration search window
# ---------------------------------------------------------------------------
def a_configuration(tmp_path, found_omc, n_unfound, label=' f26d _1G2H'):
    """A configuration with `found_omc` found levels and `n_unfound` unfound.

    The level to be discarded is row 1, and its O-C is the first of the list.
    """
    rows = []
    for k, omc in enumerate(found_omc, start=1):
        rows.append(enlev_record(k, 100000.0 + 10 * k, 0.05,
                                 100000.0 + 10 * k + omc, True, '5.5', label))
    for m in range(n_unfound):
        n = len(found_omc) + 1 + m
        e = 100000.0 + 10 * n
        rows.append(enlev_record(n, e, 5000.0, e, False, '5.5', label))
    return an_enlev(tmp_path, rows)


def test_the_window_narrows_when_a_bad_level_stops_counting(tmp_path):
    """The point of unfinding rather than deleting: the level was inflating
    the window every other level of its configuration is searched in."""
    enlev = a_configuration(tmp_path, [136.0, 10.0, -8.0, 6.0, -12.0], 3)
    said = []
    cfg, n_before, n_after, rms_before, rms_after, n_window = \
        DL.configuration_window(enlev, 1, said.append)
    assert cfg == 'f26d'
    assert (n_before, n_after) == (5, 4)
    assert n_window == 4          # the three unfound, plus the one being given up
    assert rms_after < rms_before
    # sqrt((136^2 + 10^2 + 8^2 + 6^2 + 12^2)/5) against sqrt((10^2 + 8^2 +
    # 6^2 + 12^2)/4): one level 136 cm^-1 out was setting the window for all
    # of them, and without it the window is six times narrower.
    assert rms_before == pytest.approx(61.384, abs=0.002)
    assert rms_after == pytest.approx(9.274, abs=0.002)
    assert any('search window' in s for s in said)


def test_below_the_threshold_the_placeholder_survives_and_is_reported(tmp_path):
    """An rms over two levels measures nothing, so sync_IDEN2.py leaves those
    rows alone and their unfound levels keep the 5000.000 placeholder.  That
    is the truth, and the run says so rather than passing it over."""
    import sync_IDEN2
    assert sync_IDEN2.MIN_CFG_LEVELS == 3
    enlev = a_configuration(tmp_path, [136.0, 10.0, -8.0], 2)
    said = []
    _cfg, n_before, n_after, _rb, _ra, _n = \
        DL.configuration_window(enlev, 1, said.append)
    assert (n_before, n_after) == (3, 2)
    assert any('WARNING' in s and '5000.000' in s for s in said)


def test_a_configuration_still_above_the_threshold_is_not_warned_about(
        tmp_path):
    enlev = a_configuration(tmp_path, [136.0, 10.0, -8.0, 6.0], 2)
    said = []
    DL.configuration_window(enlev, 1, said.append)
    assert not any('WARNING' in s for s in said)


def test_after_the_discard_only_one_state_is_reported(tmp_path):
    """Called again once the level is unfound there is nothing to compare, and
    the report says the window rather than pretending to a change."""
    enlev = a_configuration(tmp_path, [10.0, -8.0, 6.0, -12.0], 3)
    said = []
    # Row 5 is one of the unfound ones.
    _cfg, n_before, n_after, _rb, _ra, n_window = \
        DL.configuration_window(enlev, 5, said.append)
    assert n_before == n_after == 4
    assert n_window == 3
    assert not any('->' in s for s in said)


# ---------------------------------------------------------------------------
# D. The evidence guard
# ---------------------------------------------------------------------------
def an_audit(tmp_path, **cells):
    path = str(tmp_path / 'level_positions.csv')
    row = {'level_id': GONE, 'E': '111039.66', 'ln_R': '-5.492',
           'n_alt': '405', 'action': 'no support'}
    row.update({k: str(v) for k, v in cells.items()})
    write_csv(path, list(row), [row])
    return path


def test_a_negative_ln_r_passes(tmp_path):
    args = Args(audit=an_audit(tmp_path))
    assert DL.check_the_evidence(GONE, args, lambda *_a: None) == \
        pytest.approx(-5.492)


def test_a_positive_ln_r_is_refused(tmp_path):
    """The lines are MORE likely with the level there; this is not a discard."""
    args = Args(audit=an_audit(tmp_path, ln_R='3.2'))
    with pytest.raises(INL.Abort, match='move_level.py is the tool'):
        DL.check_the_evidence(GONE, args, lambda *_a: None)


def test_a_firm_alternate_is_refused(tmp_path):
    """A level with somewhere to go wants move_level.py."""
    args = Args(audit=an_audit(tmp_path, action='firm'))
    with pytest.raises(INL.Abort, match='somewhere to go'):
        DL.check_the_evidence(GONE, args, lambda *_a: None)


def test_no_audit_row_is_refused(tmp_path):
    """The audit is the evidence; without it the discard rests on nothing."""
    args = Args(audit=an_audit(tmp_path, level_id='059003.000999'))
    with pytest.raises(INL.Abort, match='run  python level_positions.py'):
        DL.check_the_evidence(GONE, args, lambda *_a: None)


def test_force_overrides_every_one_of_them(tmp_path):
    for cells in ({'ln_R': '3.2'}, {'action': 'firm'},
                  {'level_id': '059003.000999'}):
        args = Args(audit=an_audit(tmp_path, **cells), force=True)
        DL.check_the_evidence(GONE, args, lambda *_a: None)
