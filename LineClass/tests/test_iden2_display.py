"""Tests for [hfs] iden2_display: IDEN2 showing what LOPT is given.

Run from the LineClass directory:  python -m pytest tests -q

With the head-frame hyperfine correction on, LOPT is given each assigned line
moved by its hfs shift, with the shift's uncertainty and any hfs width added
to the line's.  IDEN2 draws dlv.dat, so by default it shows the measured line
and the departure it shows is not LOPT's O-C.  With iden2_display = 'lopt',
sync_IDEN2.py writes every line the fit uses into dlv.dat as LOPT is given
it, lists those rows in IDEN2/dlv_shown.txt, and check_sync.py checks them
against LOPT_input_lines.txt.

The miniature set used throughout:

    line  wn_key      measured    u      LOPT record
    1     30000.100   30000.000   0.10   29999.700  0.130   accepted, shifted
    2     29999.900   29999.800   0.10   none                unassigned
    3     25000.100   25000.000   0.05   25000.000  0.050   accepted, no shift
    4     20000.100   20000.000   0.08   19999.900  0.090   flagged P, shifted

Line 1, moved by -0.300, passes line 2, so the two rows change place.
"""
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import check_sync                          # noqa: E402
import config                              # noqa: E402
import hfs_correction                      # noqa: E402
import sync_IDEN2 as sync                  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

LINES = [(30000.100, 30000.000, 0.10),
         (29999.900, 29999.800, 0.10),
         (25000.100, 25000.000, 0.05),
         (20000.100, 20000.000, 0.08)]
LOW, UPP = '059003.000004', '059003.000003'


def quiet(*_a):
    pass


def dlv_record(wn, row, u_lam=0.0100, code=62, character='c'):
    """One row of dlv.dat at its fixed width of 66 characters."""
    rec = ('%5d%14.3f%14.4f  /%10s/%13.4f%6d'
           % (code, wn, sync.wavelength_of(wn), character, u_lam, row))
    assert len(rec) == sync.DLV_WIDTH, len(rec)
    return rec


def make_lopt_record(wn, unc, low, upp, flag):
    """A record of LOPT_input_lines.txt in make_LOPT_input's own columns."""
    buf = [' '] * 93

    def put(span, text):
        buf[span[0]:span[0] + len(text)] = text
    put(sync.FIELD_WN, '%.3f' % wn)
    put(sync.FIELD_UNC, '%.3f' % unc)
    put(sync.FIELD_LOW, low)
    put(sync.FIELD_UPP, upp)
    put(sync.FIELD_FLAGS, flag)
    put((89, 93), 'cm-1')
    return ''.join(buf).rstrip() + '\r\n'


def lopt_input(tmp_path, records=None, shifts=None):
    """LOPT_input_lines.txt and LOPT_hfs_shifts.txt for the miniature set."""
    if records is None:
        records = [(29999.700, 0.130, '059003.000010', '059003.000020', ' '),
                   (25000.000, 0.050, '059003.000011', '059003.000021', ' '),
                   (19999.900, 0.090, '059003.000012', '059003.000022', 'P')]
    if shifts is None:
        shifts = [('059003.000010', '059003.000020', '30000.0', 29999.700,
                   -0.300, 0.040),
                  ('059003.000012', '059003.000022', '20000.0', 19999.900,
                   -0.100, 0.030)]
    path = tmp_path / 'LOPT_input_lines.txt'
    with open(path, 'w', encoding='latin-1', newline='') as fh:
        for rec in records:
            fh.write(make_lopt_record(*rec))
    hfs_correction.write_shifts(hfs_correction.shifts_path(str(path)), shifts)
    return str(path)


def measured_rows():
    """dlv.dat as a 'measured' sync leaves it."""
    raw = [dlv_record(wn, k + 1) for k, (_key, wn, _u) in enumerate(LINES)]
    rows, _rep = sync.rewrite_dlv(raw, LINES, quiet)
    return rows


def wn_of(rec):
    return float(rec[sync.DLV_WN[0]:sync.DLV_WN[1]])


def row_of(rec):
    return int(rec[sync.DLV_ROW[0]:sync.DLV_ROW[1]])


def u_of(rec):
    """The uncertainty a row states, on the wavenumber scale."""
    lam = float(rec[sync.DLV_LAMBDA[0]:sync.DLV_LAMBDA[1]])
    return float(rec[sync.DLV_UNC[0]:sync.DLV_UNC[1]]) * wn_of(rec) / lam


# ---------------------------------------------------------------------------
# The configuration
# ---------------------------------------------------------------------------
def write_config(tmp_path, hfs_text):
    path = tmp_path / 'lineclass_config.toml'
    base = os.path.join(HERE, 'lineclass_config.toml').replace('\\', '/')
    path.write_text('inherit = "%s"\n\n[files]\noutput = "o.xlsx"\n'
                    'output_csv = "o.csv"\n\n[hfs]\n%s' % (base, hfs_text),
                    encoding='utf-8')
    return str(path)


def test_the_default_is_measured():
    assert config.HfsSettings().iden2_display == 'measured'
    assert config.load(os.path.join(HERE, 'lineclass_config.toml')) \
        .hfs.iden2_display == 'measured'


def test_lopt_is_accepted_with_the_correction_on(tmp_path):
    path = write_config(tmp_path, 'apply = true\niden2_display = "lopt"\n')
    assert config.load(path).hfs.iden2_display == 'lopt'


def test_lopt_without_the_correction_is_refused(tmp_path):
    """With the correction off LOPT is given the measured lines, so 'lopt'
    would show nothing different and is taken for a mistake."""
    path = write_config(tmp_path, 'iden2_display = "lopt"\n')
    with pytest.raises(config.ConfigError, match='needs apply'):
        config.load(path)


def test_an_unknown_display_is_refused(tmp_path):
    path = write_config(tmp_path, 'apply = true\niden2_display = "head"\n')
    with pytest.raises(config.ConfigError, match='iden2_display'):
        config.load(path)


# ---------------------------------------------------------------------------
# What LOPT is given, line by line
# ---------------------------------------------------------------------------
def test_a_line_the_fit_uses_takes_lopts_values(tmp_path):
    shown = sync.lopt_view(LINES, lopt_input(tmp_path), quiet)
    assert shown[0] == (30000.100, 29999.700, 0.130)


def test_an_unassigned_line_and_a_rejected_one_stay_as_measured(tmp_path):
    shown = sync.lopt_view(LINES, lopt_input(tmp_path), quiet)
    assert shown[1] == LINES[1]
    assert shown[3] == LINES[3]


def test_lopts_print_of_the_measured_value_keeps_the_measured_value(tmp_path):
    """LOPT is given three decimals; the measured line, printed, is not a
    difference worth showing."""
    lines = list(LINES)
    lines[2] = (25000.100, 25000.0004, 0.0503)
    shown = sync.lopt_view(lines, lopt_input(tmp_path), quiet)
    assert shown[2] == lines[2]


def test_one_line_given_to_lopt_two_ways_stops_the_run(tmp_path):
    records = [(25000.000, 0.050, '059003.000011', '059003.000021', ' '),
               (25000.000, 0.070, '059003.000013', '059003.000023', ' ')]
    with pytest.raises(sync.SyncError, match='two uncertainties'):
        sync.lopt_view(LINES, lopt_input(tmp_path, records, []), quiet)


# ---------------------------------------------------------------------------
# dlv.dat
# ---------------------------------------------------------------------------
def lopt_sync(tmp_path, rows=None, previous=None):
    shown = sync.lopt_view(LINES, lopt_input(tmp_path), quiet)
    return sync.rewrite_dlv(rows or measured_rows(), LINES, quiet, shown,
                            previous)


def test_the_row_shows_lopts_wavenumber_and_uncertainty(tmp_path):
    rows, rep = lopt_sync(tmp_path)
    got = [r for r in rows if row_of(r) == 1][0]
    assert wn_of(got) == 29999.700
    assert abs(u_of(got) - 0.130) < 5e-4
    # the wavelength is the standard one for the wavenumber shown, so IDEN2
    # rebuilds the same wavenumber from it
    assert float(got[sync.DLV_LAMBDA[0]:sync.DLV_LAMBDA[1]]) == \
        round(sync.wavelength_of(29999.700), 4)
    assert rep['shown'] == [(1, 30000.100, 30000.000, 0.10, 29999.700, 0.130)]


def test_a_row_that_passes_its_neighbor_changes_place(tmp_path):
    """The rows stay in decreasing order of what they show, each with its own
    line number, as IDEN2 places a line inserted on its screen."""
    rows, rep = lopt_sync(tmp_path)
    assert [row_of(r) for r in rows] == [2, 1, 3, 4]
    assert rep['n_reordered'] == 2


def test_a_second_lopt_sync_changes_nothing(tmp_path):
    once, rep = lopt_sync(tmp_path)
    previous = {m[0]: m[1:] for m in rep['shown']}
    twice, rep2 = lopt_sync(tmp_path, once, previous)
    assert twice == once and rep2['changed'] == []


def test_a_shown_row_is_still_its_line_without_the_record(tmp_path):
    once, _rep = lopt_sync(tmp_path)
    twice, rep = lopt_sync(tmp_path, once, {})
    assert twice == once and not rep['unmatched'] and not rep['absent']


def test_going_back_to_measured_restores_the_file(tmp_path):
    once, rep = lopt_sync(tmp_path)
    previous = {m[0]: m[1:] for m in rep['shown']}
    back, rep2 = sync.rewrite_dlv(once, LINES, quiet, None, previous)
    assert back == measured_rows()
    assert rep2['shown'] == []


def test_a_new_shift_finds_the_row_through_the_record(tmp_path):
    """LOPT given the line somewhere else since the last sync: the row shows
    neither the measured nor the new value, and the record says whose it
    is."""
    once, rep = lopt_sync(tmp_path)
    previous = {m[0]: m[1:] for m in rep['shown']}
    records = [(29999.400, 0.150, '059003.000010', '059003.000020', ' ')]
    shifts = [('059003.000010', '059003.000020', '30000.0', 29999.400,
               -0.600, 0.050)]
    shown = sync.lopt_view(LINES, lopt_input(tmp_path, records, shifts),
                           quiet)
    again, rep2 = sync.rewrite_dlv(once, LINES, quiet, shown, previous)
    assert not rep2['unmatched']
    assert [wn_of(r) for r in again if row_of(r) == 1] == [29999.400]


def test_a_shown_value_on_another_lines_wavenumber_stops_the_run(tmp_path):
    lines = LINES + [(29999.750, 29999.700, 0.10)]
    raw = [dlv_record(wn, k + 1) for k, (_key, wn, _u) in enumerate(lines)]
    shown = sync.lopt_view(lines, lopt_input(tmp_path), quiet)
    with pytest.raises(sync.SyncError, match='another row'):
        sync.rewrite_dlv(raw, lines, quiet, shown)


def test_the_record_reads_back_what_was_written(tmp_path):
    _rows, rep = lopt_sync(tmp_path)
    sync.write_shown(sync.shown_path(str(tmp_path)), rep['shown'])
    assert sync.read_shown(str(tmp_path)) == {
        1: (30000.100, 30000.000, 0.10, 29999.700, 0.130)}
    assert sync.read_shown(str(tmp_path / 'nowhere')) == {}


def test_an_accepted_transition_finds_its_row_at_the_shown_value(tmp_path):
    path = lopt_input(tmp_path)
    shown = sync.lopt_view(LINES, path, quiet)
    rows, _rep = sync.rewrite_dlv(measured_rows(), LINES, quiet, shown)
    lopt_lines = sync.as_shown(sync.read_lopt_transitions(path), LINES, shown)
    by_wn = sync.dlv_rows_by_wavenumber(rows)
    accepted = {p: wn for p, (ok, wn) in lopt_lines.items() if ok}
    assert accepted == {('059003.000010', '059003.000020'): 29999.700,
                        ('059003.000011', '059003.000021'): 25000.000}
    assert [by_wn[round(wn, 3)][0] for wn in sorted(accepted.values())] \
        == [3, 1]


# ---------------------------------------------------------------------------
# check_sync.py
# ---------------------------------------------------------------------------
class FakeLayout(object):
    def __init__(self, columns):
        self.columns = columns
        self.sheet = None


class FakeConfig(object):
    """Only what read_line_list asks of a configuration."""

    def __init__(self, lines_file):
        self.lines_file = lines_file
        self.lines = FakeLayout({'wn': 'own_corr', 'u_wn': 'unc_own_corr',
                                 'wn_key': 'own'})
        self.inflated_unc = ''


def line_list(tmp_path):
    path = tmp_path / 'lines.xlsx'
    pd.DataFrame({'own': [t[0] for t in LINES],
                  'own_corr': [t[1] for t in LINES],
                  'unc_own_corr': [t[2] for t in LINES]}).to_excel(
        path, index=False)
    return FakeConfig(str(path))


def check(tmp_path, records, lopt=True, previous=None):
    rows = check_sync._dlv_internal(records, check_sync.Report(), 12)
    rep = check_sync.Report()
    check_sync._dlv_vs_lines(rows, line_list(tmp_path), rep,
                             check_sync.DEF_DLV_TOL, 12,
                             lopt_input(tmp_path) if lopt else None,
                             previous)
    return rep


def severities(rep):
    return [f['severity'] for f in rep.findings
            if f['check'] == 'IDEN2/dlv.dat']


def items(rep, severity):
    out = []
    for f in rep.findings:
        if f['severity'] == severity:
            out.extend(f['items'])
    return out


def test_a_file_synced_in_lopt_mode_passes(tmp_path):
    rows, _rep = lopt_sync(tmp_path)
    rep = check(tmp_path, rows)
    assert severities(rep) == [check_sync.OK]


def test_a_file_still_measured_in_lopt_mode_is_out_of_step(tmp_path):
    rep = check(tmp_path, measured_rows())
    assert check_sync.ERROR not in severities(rep)
    assert items(rep, check_sync.WARN) == [
        'line      1     30000.000 ->    29999.700  (-0.300)']


def test_a_file_left_in_lopt_mode_by_a_measured_set_is_out_of_step(tmp_path):
    rows, rep = lopt_sync(tmp_path)
    previous = {m[0]: m[1:] for m in rep['shown']}
    rep = check(tmp_path, rows, lopt=False, previous=previous)
    assert check_sync.ERROR not in severities(rep)
    assert items(rep, check_sync.WARN) == [
        'line      1     29999.700 ->    30000.000  (+0.300)']


def test_an_identification_stands_for_the_measured_line(tmp_path):
    """Checks 9 and 11 compare an identified line with the classification,
    which holds the measured wavenumber."""
    rows, rep = lopt_sync(tmp_path)
    previous = {m[0]: m[1:] for m in rep['shown']}
    iden = check_sync.Iden2(None, None, {}, rows, previous)
    assert iden.wn_of_line[1] == 30000.000
    assert iden.wn_of_line[3] == 25000.000
