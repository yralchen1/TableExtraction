"""Tests for `insert_new_level` and for the two changes it needed elsewhere.

Run from the LineClass directory:  python -m pytest tests -q

`insert_new_level.py` runs three other programs and LOPT, so the tests here do
not run the chain.  They fix the pieces that decide what the chain is asked to
do, and the two behaviours that make it safe to run at all:

  * the delimiter of `new_levels` follows the extension - tab for the `.txt`
    the analyst can open in Excel without `7/2` turning into a date, comma for
    the `.csv` that came before it - and the two columns saying which
    calculated level a new level is are read;
  * a new level's calculated transitions, derived on the fly from the Cowan
    transition list, are the ones `icalc_new.xlsx` used to carry by hand, on
    the same intensity scale and cut off at the same gA;
  * a line much stronger than predicted is never proposed, and one too far
    from the prediction is not proposed either;
  * an observed line that gains a second component has ALL its unflagged
    weights divided again, the flagged ones left at zero, and the records the
    run did not touch come out of the file byte for byte as they went in;
  * the backups put every file back exactly as it was.
"""
import io
import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (ROOT, os.path.dirname(ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

import insert_new_level as INL           # noqa: E402
import classify_lines as CL              # noqa: E402
import make_LOPT_input as MLI            # noqa: E402
from models import EnergyLevel, SpectralLine   # noqa: E402


# ---------------------------------------------------------------------------
# new_levels: the delimiter, and the two new columns
# ---------------------------------------------------------------------------
TAB_FILE = ('level_id\tE\tJ\tparity\tiden2_row\tcowan_lid\tcomment\n'
            '059003.000900\t100000.5\t7/2\to\t742\t294\tfound in IDEN2\n')
CSV_FILE = ('level_id,E,J,parity,comment\n'
            '059003.000900,100000.5,7/2,o,found in IDEN2\n')


def _levels():
    """A level list of one published level, to add to."""
    lev = EnergyLevel(level_id='059003.000001', energy=0.0, parity='e',
                      J_str='4', J_val=4.0)
    return {lev.level_id: lev}, [lev]


@pytest.mark.parametrize('name, text, row, lid', [
    ('new_levels.txt', TAB_FILE, 742, 294),
    ('new_levels.csv', CSV_FILE, 0, 0),
])
def test_delimiter_follows_the_extension(tmp_path, name, text, row, lid):
    path = tmp_path / name
    path.write_text(text, encoding='utf-8')
    d, lst = _levels()
    assert CL.add_new_levels(d, lst, str(path)) == 1
    added = d['059003.000900']
    assert added.energy == pytest.approx(100000.5)
    assert added.J_str == '7/2' and added.parity == 'o'
    assert added.is_new == 1 and added.is_added == 1
    # the two optional columns, absent from the .csv, are zero when absent
    assert added.iden2_row == row
    assert added.cowan_lid == lid


def test_a_tab_file_read_as_csv_would_lose_the_columns(tmp_path):
    """The point of the extension rule: the wrong delimiter reads one column.

    A tab-separated row read with the comma reader is a single field, so the
    level_id would carry the whole line.  add_new_levels must not do that.
    """
    path = tmp_path / 'new_levels.txt'
    path.write_text(TAB_FILE, encoding='utf-8')
    d, lst = _levels()
    CL.add_new_levels(d, lst, str(path))
    assert '059003.000900' in d


# ---------------------------------------------------------------------------
# The derived calculated transitions
# ---------------------------------------------------------------------------
class _FakeTrans(object):
    """The three columns read_cowan_transitions selects on, as a DataFrame
    would answer them."""

    def __init__(self, rows):
        self.rows = rows

    def __getitem__(self, key):
        if isinstance(key, str):
            return _Col([r[key] for r in self.rows])
        return _FakeTrans([r for r, keep in zip(self.rows, key.values)
                           if keep])

    def iterrows(self):
        return enumerate(self.rows)


class _Col(object):
    def __init__(self, values):
        self.values = values

    def __eq__(self, other):
        return _Col([v == other for v in self.values])

    def __or__(self, other):
        return _Col([a or b for a, b in zip(self.values, other.values)])


def test_derived_transitions_land_on_the_pipeline_intensity_scale(monkeypatch):
    """Icalc is C*gA*(rwn/1e8)*exp(-Eup/kT) with the ADOPTED energies.

    Not the calculated ones, and not a column copied from a workbook that a
    later refit of C and kT would leave stale.
    """
    new = EnergyLevel(level_id='059003.000900', energy=100000.0, parity='o',
                      J_str='5/2', J_val=2.5, is_new=1, is_added=1,
                      cowan_lid=294)
    low = EnergyLevel(level_id='059003.000047', energy=20000.0, parity='e',
                      J_str='7/2', J_val=3.5)
    levels = {new.level_id: new, low.level_id: low}

    rows = [
        # taken: partner known, gA above the cutoff
        {'lid1': 294, 'lid2': 7, 'id1': '', 'id2': '059003.000047',
         'gA': 5.0e4, 'u_gA_pct': 50.0},
        # dropped: gA below the printing cutoff Icalc.xlsx itself obeys
        {'lid1': 294, 'lid2': 8, 'id1': '', 'id2': '059003.000047',
         'gA': 5.0e2, 'u_gA_pct': 50.0},
        # dropped: the partner is not in the level list
        {'lid1': 294, 'lid2': 9, 'id1': '', 'id2': '',
         'gA': 5.0e4, 'u_gA_pct': 50.0},
        # not this level at all
        {'lid1': 3, 'lid2': 9, 'id1': '', 'id2': '059003.000047',
         'gA': 5.0e4, 'u_gA_pct': 50.0},
    ]
    monkeypatch.setattr(CL.cowan_gA if hasattr(CL, 'cowan_gA') else CL,
                        'read_transitions', lambda **kw: _FakeTrans(rows),
                        raising=False)
    import cowan_gA
    monkeypatch.setattr(cowan_gA, 'read_transitions',
                        lambda **kw: _FakeTrans(rows))

    index = {}
    n = CL.read_cowan_transitions(levels, index)
    assert n == 1
    key = ('059003.000047', '059003.000900')
    assert key in index

    C = float(CL.CFG.intensity_model['C'])
    kT = float(CL.CFG.intensity_model['kT'])
    rwn = 80000.0
    want = C * 5.0e4 * (rwn / 1e8) * math.exp(-100000.0 / kT)
    assert index[key]['calc_intensity'] == pytest.approx(want, rel=1e-12)
    assert index[key]['u_calc'] == pytest.approx(math.log(1.5), rel=1e-12)


def test_a_row_already_in_the_table_is_not_replaced(monkeypatch):
    """files.icalc_extra keeps the last word: a hand-made row survives."""
    new = EnergyLevel(level_id='059003.000900', energy=100000.0, parity='o',
                      J_str='5/2', J_val=2.5, cowan_lid=294)
    low = EnergyLevel(level_id='059003.000047', energy=20000.0, parity='e',
                      J_str='7/2', J_val=3.5)
    levels = {new.level_id: new, low.level_id: low}
    rows = [{'lid1': 294, 'lid2': 7, 'id1': '', 'id2': '059003.000047',
             'gA': 5.0e4, 'u_gA_pct': 50.0}]
    import cowan_gA
    monkeypatch.setattr(cowan_gA, 'read_transitions',
                        lambda **kw: _FakeTrans(rows))
    index = {('059003.000047', '059003.000900'):
             {'calc_intensity': 1.25, 'u_calc': 0.5, 'assigned_to': None}}
    assert CL.read_cowan_transitions(levels, index) == 0
    assert index[('059003.000047', '059003.000900')]['calc_intensity'] == 1.25


# ---------------------------------------------------------------------------
# Choosing the lines
# ---------------------------------------------------------------------------
def _candidate(rwn, wn, i_obs, i_calc, u_calc=0.6):
    line = SpectralLine(wavenumber=wn, wn_uncertainty=0.3, intensity=i_obs,
                        line_character='')
    return INL.Candidate('059003.000047', '059003.000900', rwn, line,
                         i_calc, u_calc)


def test_a_line_much_stronger_than_predicted_is_left_free():
    """The 91856.116 case: 58 times the predicted intensity, 4.5 sigma."""
    strong = _candidate(91856.713, 91856.116, 3.74e5, 6413.0, 0.59)
    INL.choose([strong], marked=[], strong_sigma=3.0, propose_window=1.0)
    assert strong.source == ''
    assert strong.verdict.startswith('left free')
    assert strong.z_intensity > 4.0


def test_a_line_the_prediction_accounts_for_is_proposed():
    ok = _candidate(90917.873, 90917.831, 2912.0, 1882.0, 0.68)
    INL.choose([ok], marked=[], strong_sigma=3.0, propose_window=1.0)
    assert ok.source == 'proposed'


def test_a_line_too_far_from_the_prediction_is_not_proposed():
    far = _candidate(22556.197, 22554.249, 177.5, 192.6, 0.60)
    INL.choose([far], marked=[], strong_sigma=3.0, propose_window=1.0)
    assert far.source == ''
    assert 'too far' in far.verdict


def test_a_hand_mark_is_taken_however_strong_the_line_is():
    """What the analyst marked in IDEN2 is a decision already taken."""
    strong = _candidate(94049.953, 94050.809, 7.882e4, 3734.0, 0.49)
    INL.choose([strong], marked=[94050.809], strong_sigma=3.0,
               propose_window=1.0)
    assert strong.source == 'IDEN2'
    assert strong.z_intensity > 3.0


def test_one_transition_takes_at_most_one_line():
    a = _candidate(1000.0, 1000.1, 100.0, 100.0)
    b = _candidate(1000.0, 1000.4, 100.0, 100.0)
    INL.choose([a, b], marked=[], strong_sigma=3.0, propose_window=1.0)
    assert a.source == 'proposed'
    assert b.source == '' and 'already has a line' in b.verdict


# ---------------------------------------------------------------------------
# The LOPT input: reading it, dividing the weights, writing it back
# ---------------------------------------------------------------------------
def _lopt_file(tmp_path, rows):
    path = tmp_path / 'LOPT_input_lines.txt'
    with io.open(str(path), 'w', encoding='ascii', newline='') as fh:
        for wn, unc, intens, low, upp, flag, weight in rows:
            fh.write(MLI.format_line(wn, unc, intens, low, upp, flag, weight)
                     + MLI.EOL_LINES)
    return path


def test_untouched_records_come_out_byte_for_byte(tmp_path):
    """Including the blends whose weight column holds a raw intensity.

    A few records of the real file carry `35318.` rather than a fraction, from
    a run made before the weights were normalised.  Rewriting one of those in
    the fraction format overflows the six-character column, so a record this
    run does not change is passed through exactly as it stands.
    """
    path = _lopt_file(tmp_path, [
        (1000.0, 0.5, 10.0, '059003.000001', '059003.000100', '', 1.0),
        (900.0, 0.5, 20.0, '059003.000002', '059003.000200', 'P', 0.0),
    ])
    raw = io.open(str(path), encoding='ascii', newline='').read()
    # a weight that no fraction format could reproduce
    raw = raw.replace('1.0000  cm-1', '35318.  cm-1')
    io.open(str(path), 'w', encoding='ascii', newline='').write(raw)

    rows = INL.read_lopt_input(str(path))
    assert [r['weight'] for r in rows] == [35318.0, 0.0]
    INL.write_lopt_input(str(path), rows)
    assert io.open(str(path), encoding='ascii', newline='').read() == raw


def test_a_new_component_divides_the_weight_and_leaves_P_alone(tmp_path):
    """The rule asked for: every unflagged record of a touched line is
    re-weighted, in proportion to the calculated intensities; the P records,
    which LOPT carries but does not fit, stay at zero."""
    path = _lopt_file(tmp_path, [
        (1000.0, 0.5, 10.0, '059003.000001', '059003.000100', '', 1.0),
        (1000.0, 0.5, 10.0, '059003.000003', '059003.000300', 'P', 0.0),
        (800.0, 0.5, 10.0, '059003.000004', '059003.000400', '', 1.0),
    ])
    rows = INL.read_lopt_input(str(path))
    rows.append({'wn': 1000.0, 'unc': 0.5, 'intens': 10.0,
                 'low_id': '059003.000002', 'upp_id': '059003.000200',
                 'flag': '', 'weight': 1.0, 'raw': None})
    calc = {('059003.000001', '059003.000100'): 3.0,
            ('059003.000002', '059003.000200'): 1.0,
            ('059003.000003', '059003.000300'): 99.0,
            ('059003.000004', '059003.000400'): 5.0}
    INL.reweigh(rows, [1000.0], calc, lambda *_a: None)
    got = {(r['low_id'], r['flag']): r['weight'] for r in rows}
    assert got[('059003.000001', '')] == pytest.approx(0.75)
    assert got[('059003.000002', '')] == pytest.approx(0.25)
    assert got[('059003.000003', 'P')] == 0.0
    # the line the run did not touch keeps its weight
    assert got[('059003.000004', '')] == 1.0


def test_the_weights_of_a_touched_line_add_up_to_one(tmp_path):
    path = _lopt_file(tmp_path, [
        (1000.0, 0.5, 10.0, '059003.000001', '059003.000100', '', 1.0),
    ])
    rows = INL.read_lopt_input(str(path))
    for n, lid in enumerate(('059003.000002', '059003.000005'), start=2):
        rows.append({'wn': 1000.0, 'unc': 0.5, 'intens': 10.0,
                     'low_id': lid, 'upp_id': '059003.00030%d' % n,
                     'flag': '', 'weight': 1.0, 'raw': None})
    calc = {(r['low_id'], r['upp_id']): 1.0 + i
            for i, r in enumerate(rows)}
    INL.reweigh(rows, [1000.0], calc, lambda *_a: None)
    assert sum(r['weight'] for r in rows) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# The backups
# ---------------------------------------------------------------------------
def test_restore_puts_every_file_back(tmp_path, monkeypatch):
    a = tmp_path / 'a.txt'
    b = tmp_path / 'b.txt'
    a.write_bytes(b'one\r\ntwo\n')          # mixed endings, kept as they are
    b.write_bytes(b'\x00\x01\x02')
    monkeypatch.setattr(INL, 'BACKUP_DIR', str(tmp_path / 'backup'))
    saved = INL.preflight([str(a), str(b)], lambda *_a: None)
    assert len(saved) == 2
    a.write_bytes(b'ruined')
    b.write_bytes(b'')
    INL.restore(saved, lambda *_a: None)
    assert a.read_bytes() == b'one\r\ntwo\n'
    assert b.read_bytes() == b'\x00\x01\x02'


def test_the_preflight_names_a_file_it_cannot_write(tmp_path, monkeypatch):
    """Excel's lock has to stop the run at the start, not at the end."""
    a = tmp_path / 'a.txt'
    a.write_text('x', encoding='ascii')
    monkeypatch.setattr(INL, 'BACKUP_DIR', str(tmp_path / 'backup'))
    monkeypatch.setattr(INL.output_files, 'why_unwritable',
                        lambda p: 'it is open in Excel'
                        if p.endswith('a.txt') else None)
    with pytest.raises(SystemExit) as exc:
        INL.preflight([str(a)], lambda *_a: None)
    assert 'a.txt' in str(exc.value)


# ---------------------------------------------------------------------------
# The identifier
# ---------------------------------------------------------------------------
def test_the_next_identifier_is_the_largest_plus_one():
    assert INL.next_free_id(['059003.000001', '059003.000623',
                             '059003.000042']) == '059003.000624'


def test_j_is_written_the_way_the_level_list_writes_it():
    assert INL.j_string(2.5) == '5/2'
    assert INL.j_string(3.0) == '3'
    assert INL.j_string(0.5) == '1/2'


def test_the_assignment_tail_keeps_the_width_trans_dat_expects():
    import swap_line_assignments_IDEN as IDEN
    tail = INL.make_assignment(53, 71647.505, -0.164, 1511)
    assert len(tail) == IDEN.TR_OBS_ROW[1]
    assert IDEN.obs_wavenumber(tail) == pytest.approx(71647.505)
    assert IDEN.obs_omc(tail) == pytest.approx(-0.164)
    assert IDEN.obs_row(tail) == 1511
    assert IDEN.has_line(tail)


# ---------------------------------------------------------------------------
# A level whose lines are already marked in IDEN2 is not added to
# ---------------------------------------------------------------------------
def test_hand_marks_are_the_whole_answer_when_there_are_any():
    """A level gone through on the screen gets nothing proposed beside it.

    The analyst passed over the other lines for reasons - the branch, the
    plate, a stronger owner elsewhere - that none of the arithmetic here can
    see, so passing them over is part of the answer and not a gap in it.
    """
    marked = _candidate(90917.873, 90917.831, 2912.0, 1882.0, 0.68)
    also_good = _candidate(85072.903, 85073.215, 567.6, 2485.0, 0.60)
    INL.choose([marked, also_good], marked=[90917.831], strong_sigma=3.0,
               propose_window=1.0)
    assert marked.source == 'IDEN2'
    assert also_good.source == ''
    assert 'marked in IDEN2' in also_good.verdict


def test_proposals_are_made_when_no_line_of_the_level_is_marked():
    a = _candidate(90917.873, 90917.831, 2912.0, 1882.0, 0.68)
    b = _candidate(85072.903, 85073.215, 567.6, 2485.0, 0.60)
    b.low_id = '059003.000135'          # a different transition of the level
    INL.choose([a, b], marked=[], strong_sigma=3.0, propose_window=1.0)
    assert a.source == 'proposed' and b.source == 'proposed'


def test_propose_true_proposes_beside_the_hand_marks():
    marked = _candidate(90917.873, 90917.831, 2912.0, 1882.0, 0.68)
    other = _candidate(85072.903, 85073.215, 567.6, 2485.0, 0.60)
    other.low_id = '059003.000135'
    INL.choose([marked, other], marked=[90917.831], strong_sigma=3.0,
               propose_window=1.0, propose=True)
    assert marked.source == 'IDEN2' and other.source == 'proposed'


def test_a_mark_is_recognised_though_trans_dat_rounds_the_wavenumber():
    """trans.dat carries three decimals; the line workbook carries more.

    92677.294 in Pr3_lines.xlsx is 92677.290 in IDEN2, four thousandths away,
    and a mark that is not recognised as one would be quietly re-decided.
    """
    c = _candidate(92675.733, 92677.294, 7642.0, 3422.0, 0.60)
    INL.choose([c], marked=[92677.290], strong_sigma=3.0,
               propose_window=1.0)
    assert c.source == 'IDEN2'


# ---------------------------------------------------------------------------
# What LOPT says about the fit as a whole
# ---------------------------------------------------------------------------
def test_the_residual_sum_of_squares_is_read_out_of_what_lopt_printed():
    out = ('Adding preliminary values...Done.\n'
           '\nRSS/degrees_of_freedom = 1.16 (5261 degrees_of_freedom)\n\n'
           'Rounding levels...Done.\n')
    assert INL.lopt_rss(out) == (1.16, 5261)


def test_no_residual_sum_of_squares_is_not_an_error():
    assert INL.lopt_rss('Ok.\n') == (None, None)


# ---------------------------------------------------------------------------
# A blended component has no residual of its own, and is not silently dropped
# ---------------------------------------------------------------------------
def _lopt_output(tmp_path, rows):
    head = ['int', 'wn_o', 'uWnO1', 'uWnO2', 'uWnOTot', 'dWO-C', 'F',
            'L1', 'L2']
    path = tmp_path / 'LOPT_output_lines.txt'
    with io.open(str(path), 'w', encoding='ascii', newline='') as fh:
        fh.write('\t'.join(head) + '\n')
        for r in rows:
            fh.write('\t'.join(r) + '\n')
    return str(path)


def test_a_shared_line_is_reported_as_having_no_residual_of_its_own(tmp_path):
    path = _lopt_output(tmp_path, [
        ['1', '90917.831', '0', '0', '0.30', '-0.042', '', 'A', 'B'],
        ['1', '93276.982', '0', '0', '0.30', '_', 'c', 'C', 'B'],
    ])
    bad, worst, every, without = INL.ritz_culprits(
        path, {('A', 'B'), ('C', 'B')}, 4.0)
    assert not bad
    assert [t[0] for t in every] == ['A']
    assert [t[0] for t in without] == ['C']
    assert worst == pytest.approx(0.14, abs=0.01)


# ---------------------------------------------------------------------------
# A level enlev.dat has not been given a measured energy for
# ---------------------------------------------------------------------------
def test_a_row_with_no_measured_energy_reads_as_none(tmp_path):
    import swap_line_assignments_IDEN as IDEN
    rec = ('%4d' % 742) + ('%12.3f' % 116300.0) + (' ' * 10) + (' ' * 12)
    rec = rec + '  ' + (' ' * 12) + ('%5.1f' % 2.5) + ' /f25f/'
    path = tmp_path / 'enlev.dat'
    with io.open(str(path), 'w', encoding='ascii', newline='') as fh:
        fh.write(rec + '\n')
    e_obs, e_calc, j_val, _label = INL.enlev_row(str(path), 742)
    assert e_obs is None
    assert e_calc == pytest.approx(116300.0)
    assert j_val == pytest.approx(2.5)


def _classification_file(path, rows):
    """A minimal `line_classifications.csv` holding `rows`."""
    fields = ['wn_obs', 'unc_wn_obs', 'obs_intens', 'char', 'low_id', 'upp_id',
              'calc_intens', 'dif_wn_O-C', 'grade', 'notes2', 'accepted']
    with io.open(str(path), 'w', encoding='utf-8', newline='') as fh:
        fh.write(','.join(fields) + '\n')
        for r in rows:
            fh.write(','.join(str(r.get(k, '')) for k in fields) + '\n')
    return str(path)


def test_the_rows_of_one_observed_line_are_kept_together(tmp_path):
    path = _classification_file(tmp_path / 'c.csv', [
        {'wn_obs': 95033.381, 'low_id': 'A', 'upp_id': 'B'},
        {'wn_obs': 95033.3812, 'low_id': 'C', 'upp_id': 'D'},
        {'wn_obs': 90917.831, 'low_id': 'E', 'upp_id': 'F'},
    ])
    by_wn = INL.classification_rows(path)
    assert sorted(by_wn) == [90917.83, 95033.38]
    assert [r['low_id'] for r in by_wn[95033.38]] == ['A', 'C']
    # the workbook has a header, so the first transition is its row 2
    assert [r['xlsx_row'] for r in by_wn[95033.38]] == [2, 3]
    assert by_wn[90917.83][0]['xlsx_row'] == 4


def test_the_dossier_names_every_transition_on_the_line(tmp_path):
    path = _classification_file(tmp_path / 'now.csv', [
        {'wn_obs': 95033.381, 'unc_wn_obs': 0.361, 'obs_intens': 2409.0,
         'low_id': '059003.000099', 'upp_id': '059003.000407',
         'calc_intens': 6817.0, 'dif_wn_O-C': -1.025, 'grade': '3B',
         'notes2': 'Step2 Accepted (pair)', 'accepted': 1},
        {'wn_obs': 95033.381, 'unc_wn_obs': 0.361, 'obs_intens': 2409.0,
         'low_id': '059003.000080', 'upp_id': '059003.000624',
         'calc_intens': 13510.0, 'dif_wn_O-C': 0.707, 'grade': '2G',
         'notes2': 'Step1 Accepted New', 'accepted': 1},
    ])
    before = _classification_file(tmp_path / 'was.csv', [
        {'wn_obs': 95033.381, 'low_id': '059003.000099',
         'upp_id': '059003.000407', 'dif_wn_O-C': -1.311, 'accepted': 0},
    ])
    said = []
    now = INL.classification_rows(path)[95033.38]
    was = {(r['low_id'], r['upp_id']): r
           for r in INL.classification_rows(before)[95033.38]}
    INL.line_dossier(95033.381, now, was, '059003.000624', said.append)
    text = '\n'.join(said)
    # both components are there, the new level's marked with a star
    assert '059003.000407' in text and '059003.000624' in text
    assert '* 059003.000080' in text
    # the verdict this run changed is shown as the change it is
    assert 'accepted (was rejected)' in text
    assert 'accepted (new)' in text
    # and so is the wavenumber residual that moved
    assert 'dif_O-C was -1.311' in text
    assert 'Step2 Accepted (pair)' in text


def test_the_dossier_is_silent_about_what_did_not_change(tmp_path):
    path = _classification_file(tmp_path / 'now.csv', [
        {'wn_obs': 90917.831, 'low_id': 'A', 'upp_id': 'B',
         'calc_intens': 421.1, 'dif_wn_O-C': 0.913, 'grade': '3C',
         'notes2': 'kept', 'accepted': 1},
    ])
    before = _classification_file(tmp_path / 'was.csv', [
        {'wn_obs': 90917.831, 'low_id': 'A', 'upp_id': 'B',
         'dif_wn_O-C': 0.9132, 'accepted': 1},
    ])
    said = []
    was = {(r['low_id'], r['upp_id']): r
           for r in INL.classification_rows(before)[90917.83]}
    INL.line_dossier(90917.831, INL.classification_rows(path)[90917.83],
                     was, '059003.000624', said.append)
    text = '\n'.join(said)
    assert 'was' not in text.replace('lower', '')  # no "(was ...)" note at all
    assert 'accepted' in text


# --- --accept: adopting the assignments a previous run only named ----------
def test_accept_without_a_reason_records_that_it_was_looked_at():
    assert INL.parse_accept(['95033.381']) == \
        [(95033.381, '', INL.REASON_ADOPT)]


def test_accept_keeps_the_reason_it_is_given():
    assert INL.parse_accept(['93276.982=better CoG']) == \
        [(93276.982, '', 'better CoG')]


def test_accept_wants_a_wavenumber():
    with pytest.raises(SystemExit):
        INL.parse_accept(['the second one=good'])


def test_accept_may_be_repeated_and_keeps_its_lines_apart():
    args = INL.parse_args(['--iden2-row', '742', '--accept', '95033.381',
                           '--accept', '93276.982=better CoG'])
    assert INL.parse_accept(args.accept) == [
        (95033.381, '', INL.REASON_ADOPT), (93276.982, '', 'better CoG')]


# ---------------------------------------------------------------------------
# Which transition a mark is on, and which one an option names
# ---------------------------------------------------------------------------
def test_a_mark_belongs_to_a_transition_not_to_a_wavenumber():
    """Two transitions of the level on one blended feature, one of them marked.

    39785.512 carries 059003.000433-059003.000625 in IDEN2.  The level's other
    candidate on the same feature, 059003.000243-059003.000625, was never
    marked, and calling it "marked in IDEN2" would both overstate what the
    analyst decided and, since a mark is never second-guessed, assign it.
    """
    marked = _candidate(39785.500, 39785.512, 1000.0, 900.0)
    marked.low_id = '059003.000433'
    other = _candidate(39785.480, 39785.512, 1000.0, 5.0)
    other.low_id = '059003.000243'
    INL.choose([marked, other],
               marked={('059003.000433', '059003.000900'): 39785.512},
               strong_sigma=3.0, propose_window=1.0)
    assert marked.source == 'IDEN2'
    assert other.source == ''
    assert other.verdict != 'marked in IDEN2'


def test_the_same_transition_on_another_line_says_which_line_took_it():
    marked = _candidate(1000.0, 1000.05, 100.0, 100.0)
    elsewhere = _candidate(1000.0, 1000.40, 100.0, 100.0)
    INL.choose([marked, elsewhere],
               marked={('059003.000047', '059003.000900'): 1000.05},
               strong_sigma=3.0, propose_window=1.0, propose=True)
    assert marked.source == 'IDEN2'
    assert elsewhere.source == '' and '1000.050' in elsewhere.verdict


def test_reject_names_one_component_of_a_blend():
    spec, = INL.parse_reject(['39785.512/000243=Too weak'])
    assert spec.wn == 39785.512 and spec.reason == 'Too weak'
    assert spec.matches(39785.512, '059003.000243', '059003.000625')
    assert not spec.matches(39785.512, '059003.000433', '059003.000625')
    # the partner may be written out in full
    full, = INL.parse_reject(['39785.512/059003.000243=Too weak'])
    assert full.matches(39785.512, '059003.000243', '059003.000625')


def test_reject_without_a_partner_names_every_transition_on_the_line():
    spec, = INL.parse_reject(['39785.512=Too weak'])
    assert spec.matches(39785.512, '059003.000243', '059003.000625')
    assert spec.matches(39785.512, '059003.000433', '059003.000625')
    assert not spec.matches(39785.400, '059003.000433', '059003.000625')


def test_reject_still_wants_a_reason():
    with pytest.raises(SystemExit):
        INL.parse_reject(['39785.512'])


def test_accept_may_name_one_component_too():
    spec, = INL.parse_accept(['39785.512/000433'])
    assert spec.partner == '000433' and spec.reason == INL.REASON_ADOPT


# ---------------------------------------------------------------------------
# A component too faint to matter to the blend it would join
# ---------------------------------------------------------------------------
def _lopt_row(wn, low, upp, flag=''):
    return {'wn': wn, 'unc': 0.005, 'intens': 1.0, 'low_id': low,
            'upp_id': upp, 'flag': flag, 'weight': 1.0, 'raw': ''}


def test_a_hundredth_of_the_blend_is_not_proposed():
    # 5 against the 900 of the strong component: the intensity the line shows
    # is its own, so nothing here makes the candidate look too strong - the
    # share is the only test it fails.
    weak = _candidate(39785.500, 39785.512, 5.0, 5.0)
    weak.low_id = '059003.000243'
    share = INL.blend_share_of(
        [weak],
        {('059003.000433', '059003.000900'): {'calc_intensity': 900.0},
         ('059003.000243', '059003.000900'): {'calc_intensity': 5.0}},
        [_lopt_row(39785.512, '059003.000433', '059003.000900')])
    assert abs(share(weak) - 5.0 / 905.0) < 1e-9
    INL.choose([weak], marked={}, strong_sigma=3.0, propose_window=1.0,
               share=share, min_share=INL.DEF_MIN_SHARE)
    assert weak.source == ''
    assert 'too weak to contribute' in weak.verdict


def test_a_transition_alone_on_its_line_has_no_share_to_fail():
    lone = _candidate(90917.873, 90917.831, 2912.0, 1882.0, 0.68)
    share = INL.blend_share_of([lone], {}, [])
    assert share(lone) is None
    INL.choose([lone], marked={}, strong_sigma=3.0, propose_window=1.0,
               share=share, min_share=0.5)
    assert lone.source == 'proposed'


def test_a_rejected_candidate_of_the_line_is_not_part_of_the_blend():
    """LOPT's ``P`` records are shown to the fit but not fitted, so they are
    not components and do not dilute anyone's share."""
    c = _candidate(39785.500, 39785.512, 1000.0, 5.0)
    share = INL.blend_share_of(
        [c], {('x', 'y'): {'calc_intensity': 900.0}},
        [_lopt_row(39785.512, 'x', 'y', flag='P')])
    assert share(c) is None


# ---------------------------------------------------------------------------
# The candidate table
# ---------------------------------------------------------------------------
def test_the_table_tells_an_old_mark_from_a_new_one():
    old = _candidate(1000.0, 1000.05, 100.0, 100.0)
    new = _candidate(2000.0, 2000.05, 100.0, 100.0)
    new.low_id = '059003.000135'
    INL.choose([old, new],
               marked={('059003.000047', '059003.000900'): 1000.05},
               strong_sigma=3.0, propose_window=1.0, propose=True)
    out = []
    INL.candidate_table([old, new], out.append)
    body = '\n'.join(out[1:])
    assert ' old ' in body and ' new ' in body


def test_the_table_counts_what_it_passes_over_and_lists_it_on_request():
    taken = _candidate(1000.0, 1000.05, 100.0, 100.0)
    skipped = _candidate(1000.0, 1000.40, 100.0, 100.0)
    INL.choose([taken, skipped], marked={}, strong_sigma=3.0,
               propose_window=1.0)
    out = []
    INL.candidate_table([taken, skipped], out.append)
    assert any('1 further candidate' in t for t in out)
    assert not any('1000.400' in t for t in out)
    out = []
    INL.candidate_table([taken, skipped], out.append, show_skipped=True)
    assert any('1000.400' in t for t in out)


# ---------------------------------------------------------------------------
# The precision LOPT_input_lines.txt keeps
# ---------------------------------------------------------------------------
def test_the_file_is_matched_at_its_own_three_decimals():
    """81027.41941532 out of Pr3_lines.xlsx is 81027.419 in the file.

    make_LOPT_input.format_line writes the wavenumber as '%.3f', so a
    candidate and the record of it can never be equal at four decimals.  The
    run that took 059003.000625 into the fit rounded to four, decided that the
    five records the analyst had put in by hand were absent, and inserted them
    a second time; LOPT fitted each line twice and gave each copy half the
    weight.
    """
    assert INL.lopt_key(81027.41941532436) == 81027.419
    assert INL.lopt_key('81027.419') == 81027.419
    assert round(81027.41941532436, 4) != 81027.419


def test_a_record_already_in_the_file_is_recognised_at_full_precision(tmp_path):
    path = _lopt_file(tmp_path, [
        (81027.41941532436, 0.263, 1296.888,
         '059003.000118', '059003.000625', '', 1.0),
    ])
    rows = INL.read_lopt_input(str(path))
    present = {(INL.lopt_key(r['wn']), r['low_id'], r['upp_id'])
               for r in rows}
    assert (INL.lopt_key(81027.41941532436),
            '059003.000118', '059003.000625') in present


def test_a_touched_line_is_reweighed_at_full_precision(tmp_path):
    """The same rounding kept reweigh() from finding the line's records.

    With the wavenumber given at the precision of the line list, none of the
    file's records matched, the group came out empty and no weight was ever
    recomputed - the run reported '0 weight(s) changed' however many
    components the line had.
    """
    path = _lopt_file(tmp_path, [
        (39785.512, 0.048, 15946.9, '059003.000433', '059003.000625', '', 1.0),
        (39785.512, 0.048, 15946.9, '059003.000500', '059003.000700', '', 1.0),
    ])
    rows = INL.read_lopt_input(str(path))
    calc = {('059003.000433', '059003.000625'): 3.0,
            ('059003.000500', '059003.000700'): 1.0}
    n = INL.reweigh(rows, [39785.51150016838], calc, lambda *_a: None)
    assert n == 2
    got = {r['low_id']: r['weight'] for r in rows}
    assert got['059003.000433'] == pytest.approx(0.75)
    assert got['059003.000500'] == pytest.approx(0.25)


# ---------------------------------------------------------------------------
# --accept across the rounds: the component whose line this run takes a share of
#
# 059003.000629 was inserted on 59485.815, a line 059003.000121-059003.000453
# already held alone.  The new component's predicted intensity is 12 times the
# old one's, so step D cut the old one's weight from 1.0000 to 0.0746 and its
# share of the observed intensity with it.  Round 1 still accepted it and it
# was marked in IDEN2, so the run left it alone as a published identification
# nobody questioned; round 2, with the new level's accept row in the ledger,
# rejected it, and nothing recorded that.  check_sync.py then found it
# identified in trans.dat and weighted in the fit while the classification
# denied it, and the whole run was rolled back.
# ---------------------------------------------------------------------------
_OLD = (59485.815, '059003.000121', '059003.000453')
_KEY = (59485.815, '059003.000121', '059003.000453')
_ID_ROWS = {'059003.000121': 929, '059003.000453': 1161}
_ON_SCREEN = {(59485.82, frozenset((929, 1161)))}


def _adopted():
    spec = INL.Target(59485.815, '000453', 'IDEN2/LOPT; may exclude from LOPT')
    got = {}
    INL.register_adoptions([spec], {_KEY: True}, '059003.000629', got)
    return got


def test_accept_registers_the_component_the_classification_accepts():
    assert _adopted() == {
        _KEY: (59485.815, '059003.000121', '059003.000453',
               'IDEN2/LOPT; may exclude from LOPT')}


def test_accept_leaves_the_new_levels_own_assignments_to_the_other_loop():
    """A pair the level being inserted is part of is never adopted here: the
    candidate loop above writes those, with the reason REASON_ACCEPT."""
    spec = INL.Target(59485.815, '', 'looked at')
    got = {}
    INL.register_adoptions(
        [spec], {(59485.815, '059003.000134', '059003.000629'): True},
        '059003.000629', got)
    assert got == {}


def test_a_spec_that_matches_nothing_is_reported_back():
    spec = INL.Target(12345.678, '', 'looked at')
    assert INL.register_adoptions([spec], {_KEY: True}, '059003.000629', {}) \
        == [spec]


def test_a_mark_on_a_line_this_run_leaves_alone_needs_no_ledger_row():
    rows, left = INL.adoption_rows(_adopted(), {_KEY: True}, set(), set(),
                                   _ON_SCREEN, set(), _ID_ROWS, '9/16/2026')
    assert rows == []
    assert left == [_OLD]


def test_a_mark_on_a_line_this_run_takes_a_share_of_is_written():
    """The fix.  touched_wn holds the wavenumbers step D put a record on, and
    a component of one of them has had its share of the blend changed by this
    run: it is no longer an assignment nobody questioned."""
    rows, left = INL.adoption_rows(_adopted(), {_KEY: True}, set(), set(),
                                   _ON_SCREEN, {59485.82}, _ID_ROWS,
                                   '9/16/2026')
    assert left == []
    (key, row, still), = rows
    assert key == _KEY and still is True
    assert row == {'wn_obs': '59485.8150', 'low_id': '059003.000121',
                   'upp_id': '059003.000453', 'decision': 'accept',
                   'date': '9/16/2026',
                   'reason': 'IDEN2/LOPT; may exclude from LOPT'}


def test_a_component_a_later_round_rejects_is_written_however_it_is_marked():
    """The backstop.  Registered in round 1 while it was accepted, rejected in
    round 2 by this level's own rows; the screen must not exempt it, or the
    fit and IDEN2 keep an assignment the classification denies."""
    rows, left = INL.adoption_rows(_adopted(), {_KEY: False}, set(), set(),
                                   _ON_SCREEN, set(), _ID_ROWS, '9/16/2026')
    assert left == []
    (_key, _row, still), = rows
    assert still is False


def test_a_row_the_ledger_already_holds_is_not_written_again():
    rows, left = INL.adoption_rows(_adopted(), {_KEY: False}, {_KEY}, set(),
                                   _ON_SCREEN, {59485.82}, _ID_ROWS,
                                   '9/16/2026')
    assert (rows, left) == ([], [])


def test_a_row_written_earlier_in_this_run_is_not_written_again():
    rows, left = INL.adoption_rows(_adopted(), {_KEY: False}, set(), {_KEY},
                                   _ON_SCREEN, {59485.82}, _ID_ROWS,
                                   '9/16/2026')
    assert (rows, left) == ([], [])


def test_an_adopted_component_not_on_the_screen_is_written():
    rows, _left = INL.adoption_rows(_adopted(), {_KEY: True}, set(), set(),
                                    set(), set(), _ID_ROWS, '9/16/2026')
    assert len(rows) == 1


# ---------------------------------------------------------------------------
# --accept on a component the classification has ALREADY rejected
#
# The 11476.092 case.  That observed line carried 059003.000058-059003.000127
# at full weight.  The new level 059003.000630 was put on it as well, with a
# predicted intensity 55 times the sitting component's, so step D cut that
# component's share of the blend to 0.0180 and the very first round of
# classify_lines.py rejected it - "calc contribution too weak vs I_cum".  It
# was therefore never accepted at any point of the run, so nothing adopted it,
# and --accept 11476.092/000127 reported that it matched no assignment.  The
# component stayed weighted in LOPT_input_lines.txt and marked in trans.dat
# while the classification denied it, which is the two errors check_sync.py
# stops the run for, and every file was rolled back.
#
# A rejected component is adopted when this run is what rejected it: its line
# is one step D put a record on, and the fit or the screen still holds it.
# ---------------------------------------------------------------------------
_WEAK = (11476.092, '059003.000058', '059003.000127')
_WEAK_ROWS = {'059003.000058': 411, '059003.000127': 505}
_WEAK_SCREEN = {(11476.09, frozenset((411, 505)))}
_WEAK_REASON = 'IDEN2/LOPT; may as well exclude from LOPT'


def _weak_spec():
    return INL.Target(11476.092, '000127', _WEAK_REASON)


def test_a_component_this_run_rejected_and_the_fit_holds_is_adopted():
    got = {}
    missed = INL.register_adoptions(
        [_weak_spec()], {_WEAK: False}, '059003.000630', got,
        touched_wn={11476.09}, in_fit={_WEAK})
    assert missed == []
    assert got == {_WEAK: (11476.092, '059003.000058', '059003.000127',
                           _WEAK_REASON)}


def test_a_component_this_run_rejected_that_only_iden2_holds_is_adopted():
    got = {}
    INL.register_adoptions([_weak_spec()], {_WEAK: False}, '059003.000630',
                           got, touched_wn={11476.09}, in_fit=set(),
                           on_screen=_WEAK_SCREEN, id_rows=_WEAK_ROWS)
    assert list(got) == [_WEAK]


def test_a_component_rejected_on_a_line_this_run_leaves_alone_is_left():
    got = {}
    missed = INL.register_adoptions(
        [_weak_spec()], {_WEAK: False}, '059003.000630', got,
        touched_wn=set(), in_fit={_WEAK}, on_screen=_WEAK_SCREEN,
        id_rows=_WEAK_ROWS)
    assert got == {}
    assert [s.wn for s in missed] == [11476.092]


def test_a_component_neither_the_fit_nor_iden2_holds_is_left_rejected():
    # 059003.000557-059003.000407 sits on the same line as a P record at
    # weight zero: rejected long before this run, so nothing contradicts the
    # classification and --accept does not revive it.
    got = {}
    INL.register_adoptions(
        [INL.Target(11476.092, '', _WEAK_REASON)], {_WEAK: False},
        '059003.000630', got, touched_wn={11476.09}, in_fit=set(),
        on_screen=set(), id_rows=_WEAK_ROWS)
    assert got == {}


def test_an_adopted_rejected_component_gets_its_row_with_its_own_reason():
    got = {}
    INL.register_adoptions([_weak_spec()], {_WEAK: False}, '059003.000630',
                           got, touched_wn={11476.09}, in_fit={_WEAK})
    rows, left = INL.adoption_rows(got, {_WEAK: False}, set(), set(),
                                   _WEAK_SCREEN, {11476.09}, _WEAK_ROWS,
                                   '9/16/2026')
    assert left == []
    (key, row, still), = rows
    assert key == _WEAK and still is False
    assert row == {'wn_obs': '11476.0920', 'low_id': '059003.000058',
                   'upp_id': '059003.000127', 'decision': 'accept',
                   'date': '9/16/2026', 'reason': _WEAK_REASON}


def test_standing_needs_both_a_touched_line_and_a_holder():
    assert INL.standing(11476.092, _WEAK[1], _WEAK[2], {11476.09}, {_WEAK},
                        set()) is True
    assert INL.standing(11476.092, _WEAK[1], _WEAK[2], set(), {_WEAK},
                        set()) is False
    assert INL.standing(11476.092, _WEAK[1], _WEAK[2], {11476.09}, set(),
                        set()) is False


# ---------------------------------------------------------------------------
# Every program reads the added levels, and predicts them
# ---------------------------------------------------------------------------
# level_positions.py --audit once considered 593 levels and none of the 16 of
# new_levels.txt.  chance_mc.read_input_levels() - which level_shifts.py,
# level_positions.py and level_interchange.py all take their level list from -
# read the file with a comma, so its tab-separated header was a single column,
# no row had a level_id, and every one was skipped without a word.  And had
# they been read, level_shifts.load_predictions() gave them no predicted
# transition (Icalc.xlsx lacks them and icalc_extra is retired), and the audit
# leaves out a level with none.
def test_a_file_without_the_columns_is_an_error_not_an_empty_list(tmp_path):
    path = tmp_path / 'new_levels.txt'
    path.write_text(TAB_FILE.replace('\t', ';'), encoding='utf-8')
    with pytest.raises(ValueError, match='missing column'):
        CL.read_new_level_records(str(path))


def test_read_input_levels_takes_the_levels_of_a_tab_file(tmp_path,
                                                          monkeypatch):
    import chance_mc as MC
    path = tmp_path / 'new_levels.txt'
    path.write_text(TAB_FILE, encoding='utf-8')
    monkeypatch.setattr(CL, 'NEW_LEVELS', str(path))
    monkeypatch.setattr(CL, 'LEVEL_OVERRIDES', '')
    df = MC.read_input_levels()
    row = df[df['level_id'] == '059003.000900']
    assert len(row) == 1
    assert row['E_input'].iloc[0] == pytest.approx(100000.5)
    assert row['is_new_level'].iloc[0] == 1
    assert row['J'].iloc[0] == '7/2'


def test_new_level_cowan_lids(tmp_path):
    path = tmp_path / 'new_levels.txt'
    path.write_text(TAB_FILE, encoding='utf-8')
    assert CL.new_level_cowan_lids(str(path)) == {'059003.000900': 294}


def test_an_added_level_gets_predictions_from_the_cowan_list(monkeypatch):
    """Same selection, same gA cutoff, same intensity relation as the
    classification; the Ritz wavenumber from the run's own energies."""
    import level_shifts as LS
    rows = [
        {'lid1': 294, 'lid2': 7, 'id1': '', 'id2': '059003.000047',
         'gA': 5.0e4, 'u_gA_pct': 50.0},
        {'lid1': 294, 'lid2': 8, 'id1': '', 'id2': '059003.000047',
         'gA': 5.0e2, 'u_gA_pct': 50.0},         # below the cutoff
    ]
    import cowan_gA
    monkeypatch.setattr(cowan_gA, 'read_transitions',
                        lambda **kw: _FakeTrans(rows))
    monkeypatch.setattr(CL, 'new_level_cowan_lids',
                        lambda path='': {'059003.000900': 294})
    e = {'059003.000900': 100000.0, '059003.000047': 20000.0}
    out = LS._cowan_predictions(set(e), e, [])
    assert len(out) == 1
    lo, up, i_pred = out[0]
    assert (lo, up) == ('059003.000047', '059003.000900')
    C = float(CL.CFG.intensity_model['C'])
    kT = float(CL.CFG.intensity_model['kT'])
    want = C * 5.0e4 * (80000.0 / 1e8) * math.exp(-100000.0 / kT)
    assert i_pred == pytest.approx(want, rel=1e-12)
    # a pair Icalc.xlsx already predicts is not predicted twice
    assert LS._cowan_predictions(set(e), e, [(lo, up, 1.0)]) == []
