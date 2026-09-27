"""Tests for the dlv.dat checks of `check_sync`.

Run from the LineClass directory:  python -m pytest tests -q

`dlv.dat` is IDEN2's copy of the observed line list: one fixed-width row per
measured line, carrying its wavenumber, its standard wavelength and - as a
wavelength uncertainty in angstroms - how well it is known.  Nothing in IDEN2
writes it, so it is the one file that can sit on another set's wavenumbers
while every screen goes on looking right, and these tests are about saying so.

Most of them build a handful of rows in memory and read the findings out of a
Report.  The ones marked `real_file` open the project's own IDEN2 directory and
check what the analysis actually uses; they are skipped where it is absent.
"""
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import check_sync                          # noqa: E402
import sync_IDEN2 as sync                  # noqa: E402
import swap_line_assignments_IDEN as IDEN  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REAL_IDEN2 = os.path.join(HERE, 'IDEN2')

real_file = pytest.mark.skipif(
    not os.path.exists(os.path.join(REAL_IDEN2, 'dlv.dat')),
    reason='the real IDEN2 directory is absent')


def dlv_record(wn, lam, u_lam, row, code=62, character='c'):
    """One row of dlv.dat at its fixed width of 66 characters."""
    rec = ('%5d%14.3f%14.4f  /%10s/%13.4f%6d'
           % (code, wn, lam, character, u_lam, row))
    assert len(rec) == sync.DLV_WIDTH, len(rec)
    return rec


def severities(rep, check='IDEN2/dlv.dat'):
    return [f['severity'] for f in rep.findings if f['check'] == check]


def messages(rep, severity):
    return [f['message'] for f in rep.findings if f['severity'] == severity]


def items(rep, severity):
    out = []
    for f in rep.findings:
        if f['severity'] == severity:
            out.extend(f['items'])
    return out


class FakeLayout(object):
    def __init__(self, columns):
        self.columns = columns
        self.sheet = None


class FakeConfig(object):
    """Only what read_line_list asks of a configuration."""

    def __init__(self, lines_file, columns):
        self.lines_file = lines_file
        self.lines = FakeLayout(columns)


def line_list(tmp_path, own, own_corr, unc, corrected=True):
    """A miniature line workbook and the configuration that names it."""
    path = tmp_path / 'lines.xlsx'
    pd.DataFrame({'own': own, 'own_corr': own_corr,
                  'unc_own_corr': unc}).to_excel(path, index=False)
    columns = {'wn': 'own_corr', 'u_wn': 'unc_own_corr'}
    if corrected:
        columns['wn_key'] = 'own'
    return FakeConfig(str(path), columns)


# ---------------------------------------------------------------------------
# The file's own arithmetic
# ---------------------------------------------------------------------------
def test_a_well_formed_file_passes():
    rep = check_sync.Report()
    rows = check_sync._dlv_internal(
        [dlv_record(20000.100, 4999.9750, 0.0500, 1),
         dlv_record(10000.200, 9999.7500, 0.1000, 2)], rep, 12)
    assert [r.row for r in rows] == [1, 2]
    assert severities(rep) == [check_sync.OK]


def test_a_row_number_that_is_not_its_position_is_an_error():
    """Every assignment in trans.dat names its line by this number."""
    rep = check_sync.Report()
    check_sync._dlv_internal([dlv_record(20000.100, 4999.9750, 0.05, 1),
                              dlv_record(10000.200, 9999.7500, 0.10, 7)],
                             rep, 12)
    assert check_sync.ERROR in severities(rep)
    assert 'row number 7' in ' '.join(items(rep, check_sync.ERROR))


def test_two_rows_carrying_one_wavenumber_are_an_error():
    rep = check_sync.Report()
    check_sync._dlv_internal([dlv_record(20000.100, 4999.9750, 0.05, 1),
                              dlv_record(20000.100, 4999.9750, 0.05, 2)],
                             rep, 12)
    assert check_sync.ERROR in severities(rep)
    assert '20000.100' in ' '.join(items(rep, check_sync.ERROR))


def test_a_row_out_of_order_is_a_warning():
    """IDEN2 lists the lines in the order the file is in, decreasing."""
    rep = check_sync.Report()
    check_sync._dlv_internal([dlv_record(10000.200, 9999.7500, 0.10, 1),
                              dlv_record(20000.100, 4999.9750, 0.05, 2)],
                             rep, 12)
    assert severities(rep) == [check_sync.WARN]


def test_a_record_of_the_wrong_width_is_an_error():
    rep = check_sync.Report()
    check_sync._dlv_internal([dlv_record(20000.100, 4999.9750, 0.05, 1),
                              '   62     19999.000'], rep, 12)
    assert check_sync.ERROR in severities(rep)
    assert '19 characters' in ' '.join(items(rep, check_sync.ERROR))


# ---------------------------------------------------------------------------
# Against the set's own line list
# ---------------------------------------------------------------------------
def compare(tmp_path, records, own, own_corr, unc, corrected=True,
            tol=check_sync.DEF_DLV_TOL):
    rows = check_sync._dlv_internal(records, check_sync.Report(), 12)
    cfg = line_list(tmp_path, own, own_corr, unc, corrected)
    rep = check_sync.Report()
    check_sync._dlv_vs_lines(rows, cfg, rep, tol, 12)
    return rep


def test_a_file_in_step_with_the_line_list_passes(tmp_path):
    """The uncertainty is a wavelength uncertainty: a line at 20000 cm^-1 and
    5000 A known to 0.2 cm^-1 is known to 0.2 * 5000 / 20000 = 0.05 A."""
    rep = compare(tmp_path, [dlv_record(20000.000, 5000.0000, 0.0500, 1)],
                  [20000.100], [20000.000], [0.2])
    assert severities(rep) == [check_sync.OK]


def test_the_uncorrected_scale_is_named_as_what_it_is(tmp_path):
    """A set seeded by copying another set's IDEN2 keeps its wavenumbers."""
    rep = compare(tmp_path, [dlv_record(20000.100, 5000.0000, 0.0500, 1)],
                  [20000.100], [20000.000], [0.2])
    assert check_sync.ERROR in severities(rep)
    said = ' '.join(messages(rep, check_sync.ERROR))
    assert 'uncorrected own scale' in said
    listed = ' '.join(items(rep, check_sync.ERROR))
    assert '20000.100' in listed and '20000.000' in listed
    # and it is not reported as a line nobody has heard of
    assert not [f for f in rep.findings
                if f['severity'] == check_sync.WARN
                and 'match no line' in f['message']]


def test_a_wavenumber_belonging_to_no_line_is_a_warning(tmp_path):
    rep = compare(tmp_path, [dlv_record(19000.000, 5263.1000, 0.0500, 1)],
                  [20000.100], [20000.000], [0.2])
    assert check_sync.ERROR not in severities(rep)
    assert 'match no line' in ' '.join(messages(rep, check_sync.WARN))


def test_a_line_with_no_row_is_a_warning(tmp_path):
    """A row cannot be added without renumbering the file that every
    assignment refers to, so this is reported and not treated as a fault."""
    rep = compare(tmp_path, [dlv_record(20000.000, 5000.0000, 0.0500, 1)],
                  [20000.100, 31604.539], [20000.000, 31604.529], [0.2, 0.18])
    assert check_sync.ERROR not in severities(rep)
    assert 'no row in dlv.dat' in ' '.join(messages(rep, check_sync.WARN))
    assert '31604.529' in ' '.join(items(rep, check_sync.WARN))


def test_a_last_digit_difference_is_counted_and_not_listed(tmp_path):
    """Both files carry three decimals, rounded from different intermediate
    values, so one in the last digit is printing and nothing else."""
    rep = compare(tmp_path, [dlv_record(20000.001, 5000.0000, 0.0500, 1)],
                  [20000.100], [20000.000], [0.2])
    assert severities(rep) == [check_sync.OK]
    assert '1 of them differ in the last printed digit' \
        in ' '.join(messages(rep, check_sync.OK))


def test_an_uncertainty_out_of_step_is_an_error(tmp_path):
    """The user's own case: an uncertainty inflated by hand in the line list
    and never written back, so the window on the screen is half the width the
    fit is weighted by."""
    rep = compare(tmp_path, [dlv_record(20000.000, 5000.0000, 0.0500, 1)],
                  [20000.000], [20000.000], [0.4])
    assert check_sync.ERROR in severities(rep)
    listed = ' '.join(items(rep, check_sync.ERROR))
    assert '0.2000' in listed and '0.4000' in listed


def test_the_last_digit_of_the_angstrom_field_is_not_an_error(tmp_path):
    """0.0001 A is all the file can state, which at 100000 cm^-1 is 0.005
    cm^-1, and the check must not call that a disagreement."""
    rep = compare(tmp_path, [dlv_record(100000.000, 1000.0000, 0.0030, 1)],
                  [100000.000], [100000.000], [0.3025])
    assert severities(rep) == [check_sync.OK]


def test_a_line_list_that_contradicts_itself_stops_the_comparison(tmp_path):
    """Two rows of one observed line with different corrected wavenumbers:
    dlv.dat has one row for the line and no way to hold both."""
    rep = check_sync.Report()
    cfg = line_list(tmp_path, [100.0, 100.0], [100.5, 100.9], [0.01, 0.01])
    check_sync._dlv_vs_lines([], cfg, rep, check_sync.DEF_DLV_TOL, 12)
    assert check_sync.ERROR in severities(rep)
    assert 'cannot be read' in ' '.join(messages(rep, check_sync.ERROR))


# ---------------------------------------------------------------------------
# Against the identifications
# ---------------------------------------------------------------------------
class FakeTrans(object):
    """Just enough of IDEN.Trans: the records and where each transition is."""

    def __init__(self, assignments):
        self.records = []
        self.row_of = {}
        for i, (owner, partner, tail) in enumerate(assignments):
            self.records.append(sync.transition_row(partner, 40, 1000.0, True,
                                                    5000.0, tail))
            self.row_of[(owner, partner)] = i


def test_an_identification_naming_a_row_that_is_not_there_is_an_error():
    rep = check_sync.Report()
    rows = check_sync._dlv_internal(
        [dlv_record(20000.000, 5000.0000, 0.0500, 1)], check_sync.Report(), 12)
    trans = FakeTrans([(4, 3, sync.make_assignment(62, 20000.000, 0.0, 77))])
    check_sync._dlv_vs_trans(rows, trans, rep, 12)
    assert check_sync.ERROR in severities(rep)
    assert 'name row 77' in ' '.join(items(rep, check_sync.ERROR))


def test_an_identification_on_another_wavenumber_than_its_row_is_an_error():
    """IDEN2 reads the row number, so the assignment is read as one line and
    reported as another."""
    rep = check_sync.Report()
    rows = check_sync._dlv_internal(
        [dlv_record(20000.000, 5000.0000, 0.0500, 1)], check_sync.Report(), 12)
    trans = FakeTrans([(4, 3, sync.make_assignment(62, 20000.100, 0.0, 1))])
    check_sync._dlv_vs_trans(rows, trans, rep, 12)
    assert check_sync.ERROR in severities(rep)
    listed = ' '.join(items(rep, check_sync.ERROR))
    assert '20000.100' in listed and '20000.000' in listed


def test_identifications_that_agree_pass():
    rep = check_sync.Report()
    rows = check_sync._dlv_internal(
        [dlv_record(20000.000, 5000.0000, 0.0500, 1)], check_sync.Report(), 12)
    trans = FakeTrans([(4, 3, sync.make_assignment(62, 20000.000, 0.0, 1)),
                       (4, 2, IDEN.BLANK_OBS)])
    check_sync._dlv_vs_trans(rows, trans, rep, 12)
    assert severities(rep) == [check_sync.OK]
    assert 'all 1 identifications' in ' '.join(messages(rep, check_sync.OK))


# ---------------------------------------------------------------------------
# The real files
# ---------------------------------------------------------------------------
@real_file
def test_the_real_dlv_is_numbered_by_its_rows():
    records, _ends = IDEN.read_records(os.path.join(REAL_IDEN2, 'dlv.dat'))
    rep = check_sync.Report()
    rows = check_sync._dlv_internal(records, rep, 12)
    assert len(rows) > 6000
    assert check_sync.ERROR not in severities(rep)


@real_file
def test_every_real_identification_names_a_row_that_exists():
    trans_path = os.path.join(REAL_IDEN2, 'trans.dat')
    if not os.path.exists(trans_path):
        pytest.skip('trans.dat is absent')
    records, _ends = IDEN.read_records(os.path.join(REAL_IDEN2, 'dlv.dat'))
    rows = check_sync._dlv_internal(records, check_sync.Report(), 12)
    rep = check_sync.Report()
    check_sync._dlv_vs_trans(rows, IDEN.Trans(trans_path), rep, 12)
    assert severities(rep) == [check_sync.OK]
