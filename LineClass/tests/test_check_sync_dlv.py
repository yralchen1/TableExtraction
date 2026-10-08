"""Tests for the dlv.dat checks of `check_sync`.

Run from the LineClass directory:  python -m pytest tests -q

`dlv.dat` is IDEN2's copy of the observed line list: one fixed-width row per
measured line, carrying its wavenumber, its standard wavelength and - as a
wavelength uncertainty in angstroms - how well it is known.  IDEN2 rewrites it
only when a line is edited on its screen, and then rebuilds each air row's
wavenumber from its wavelength, so it is the one file that can sit on another
set's wavenumbers, or drift off its own, while every screen goes on looking
right, and these tests are about saying so.

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

    def __init__(self, lines_file, columns, inflated_unc=''):
        self.lines_file = lines_file
        self.lines = FakeLayout(columns)
        self.inflated_unc = inflated_unc


def line_list(tmp_path, own, own_corr, unc, corrected=True, registry=None):
    """A miniature line workbook and the configuration that names it.

    ``registry`` is the text of an inflated_unc_lines.txt to go with it."""
    path = tmp_path / 'lines.xlsx'
    pd.DataFrame({'own': own, 'own_corr': own_corr,
                  'unc_own_corr': unc}).to_excel(path, index=False)
    columns = {'wn': 'own_corr', 'u_wn': 'unc_own_corr'}
    if corrected:
        columns['wn_key'] = 'own'
    inflated = ''
    if registry is not None:
        inflated = str(tmp_path / 'inflated_unc_lines.txt')
        with open(inflated, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(registry)
    return FakeConfig(str(path), columns, inflated)


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


def test_a_line_inserted_in_iden2_under_a_later_number_passes():
    """IDEN2 gives an inserted line the next free number and places it by
    wavenumber, so the line numbers stop following the rows."""
    rep = check_sync.Report()
    rows = check_sync._dlv_internal(
        [dlv_record(31635.240, 3160.1170, 0.0074, 1),
         dlv_record(31604.528, 3163.1880, 0.0183, 4),
         dlv_record(31597.667, 3163.8749, 0.0082, 2),
         dlv_record(31564.068, 3167.2429, 0.0082, 3)], rep, 12)
    assert [r.row for r in rows] == [1, 4, 2, 3]
    assert severities(rep) == [check_sync.OK]
    assert '1 of them inserted' in ' '.join(messages(rep, check_sync.OK))
    assert items(rep, check_sync.OK) == ['line 4 at record 2 (31604.528)']


def test_two_rows_carrying_one_line_number_are_an_error():
    """Every assignment in trans.dat names its line by this number."""
    rep = check_sync.Report()
    check_sync._dlv_internal([dlv_record(20000.100, 4999.9750, 0.05, 1),
                              dlv_record(10000.200, 9999.7500, 0.10, 1)],
                             rep, 12)
    assert check_sync.ERROR in severities(rep)
    assert 'line 1 is on records 1 (20000.100), 2 (10000.200)' \
        in items(rep, check_sync.ERROR)


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
            tol=check_sync.DEF_DLV_TOL, registry=None):
    rows = check_sync._dlv_internal(records, check_sync.Report(), 12)
    cfg = line_list(tmp_path, own, own_corr, unc, corrected, registry)
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


def test_a_drifted_row_is_its_line_and_is_not_unknown(tmp_path):
    """The case of 2026-09-30: IDEN2 had rebuilt the wavenumbers from
    wavelengths the old sync had nudged, and 53 rows stood 0.0015 to 0.0017
    from their lines.  Each is still its line, and is reported as drifted."""
    rep = compare(tmp_path, [dlv_record(29710.207, 3364.8841, 0.0011, 5026)],
                  [29710.2085], [29710.2085], [0.01])
    said = ' '.join(messages(rep, check_sync.WARN))
    assert 'drifted' in said and 'match no line' not in said
    assert 'no row in dlv.dat' not in said
    assert '29710.2085' in ' '.join(items(rep, check_sync.WARN))
    assert any('sync_IDEN2.py' in text for _p, text in rep.actions)
    assert check_sync.OK not in severities(rep)


def test_the_last_digit_of_an_air_wavelength_near_2000_a_is_printing(
        tmp_path):
    """0.00005 A at 2001 A is 0.00125 cm^-1, and IDEN2 rebuilds the
    wavenumber from that wavelength, so 0.0016 there is still printing."""
    wn = 49961.4424
    rep = compare(tmp_path, [dlv_record(49961.444, 2000.8953, 0.0041, 1)],
                  [wn], [wn], [0.1])
    assert severities(rep) == [check_sync.OK]


def test_a_wavelength_not_standard_for_its_wavenumber_is_a_warning():
    """IDEN2 would move the wavenumber to the wavelength at its next save."""
    rep = check_sync.Report()
    w = 29710.209
    check_sync._dlv_dispersion(
        [dlv_record(w, sync.wavelength_of(w), 0.0011, 1),
         dlv_record(20000.000, 5000.0000, 0.05, 2)], rep)
    assert severities(rep) == [check_sync.WARN]
    assert 'line      2' in ' '.join(items(rep, check_sync.WARN))
    rep = check_sync.Report()
    check_sync._dlv_dispersion(
        [dlv_record(w, sync.wavelength_of(w), 0.0011, 1)], rep)
    assert severities(rep) == [check_sync.OK]


def test_an_uncertainty_out_of_step_is_an_error(tmp_path):
    """The user's own case: an uncertainty inflated by hand in the line list
    and never written back, so the window on the screen is half the width the
    fit is weighted by."""
    rep = compare(tmp_path, [dlv_record(20000.000, 5000.0000, 0.0500, 1)],
                  [20000.000], [20000.000], [0.4])
    assert check_sync.ERROR in severities(rep)
    assert 'narrower' in ' '.join(messages(rep, check_sync.ERROR))
    listed = ' '.join(items(rep, check_sync.ERROR))
    assert '0.2000' in listed and '0.4000' in listed


REGISTRY = 'wn_key\tunc_wn\tdate\treason\n'


def test_an_uncertainty_widened_in_iden2_and_registered_passes(tmp_path):
    """0.4 cm^-1 at 20000 cm^-1 and 5000 A is 0.1 A; the registry names the
    line by its wn_key, here with two decimals as the early entries were."""
    rep = compare(tmp_path, [dlv_record(20000.000, 5000.0000, 0.1000, 1)],
                  [20000.1049], [20000.000], [0.2],
                  registry=REGISTRY + '20000.10\t0.4\t\twidened\n')
    assert severities(rep) == [check_sync.OK]
    assert '1 of the uncertainties are the wider values' \
        in ' '.join(messages(rep, check_sync.OK))


def test_the_registry_never_narrows_a_line(tmp_path):
    """The rule classify_lines.py applies: a registry value smaller than the
    line list's own leaves the line list's."""
    rep = compare(tmp_path, [dlv_record(20000.000, 5000.0000, 0.0500, 1)],
                  [20000.100], [20000.000], [0.2],
                  registry=REGISTRY + '20000.1000\t0.1\t\told\n')
    assert severities(rep) == [check_sync.OK]


def test_an_uncertainty_widened_in_iden2_only_is_its_own_error(tmp_path):
    """The cure is the registry, and sync_IDEN2.py would undo the widening."""
    rep = compare(tmp_path, [dlv_record(20000.000, 5000.0000, 0.1000, 1)],
                  [20000.100], [20000.000], [0.2], registry=REGISTRY)
    assert check_sync.ERROR in severities(rep)
    said = ' '.join(messages(rep, check_sync.ERROR))
    assert 'wider' in said and 'inflated_unc_lines.txt' in said
    assert 'narrower' not in said
    assert any('inflated_unc_lines.txt' in text for _p, text in rep.actions)


def test_a_registered_line_narrower_in_iden2_is_listed_with_the_value(
        tmp_path):
    rep = compare(tmp_path, [dlv_record(20000.000, 5000.0000, 0.0500, 1)],
                  [20000.100], [20000.000], [0.2],
                  registry=REGISTRY + '20000.1000\t0.4\t\twidened\n')
    assert 'narrower' in ' '.join(messages(rep, check_sync.ERROR))
    listed = ' '.join(items(rep, check_sync.ERROR))
    assert 'pipeline    0.4000' in listed
    assert 'inflated_unc_lines.txt 0.4000' in listed


def test_a_registry_entry_naming_no_line_is_a_warning(tmp_path):
    rep = compare(tmp_path, [dlv_record(20000.000, 5000.0000, 0.0500, 1)],
                  [20000.100], [20000.000], [0.2],
                  registry=REGISTRY + '30000.5\t0.4\t\tgone\n')
    assert check_sync.ERROR not in severities(rep)
    assert items(rep, check_sync.WARN) == ['30000.5']


def test_a_registry_that_cannot_be_applied_stops_the_comparison(tmp_path):
    """An entry too short to tell two lines apart: the pipeline refuses it,
    so no uncertainty can be said to be the pipeline's."""
    rep = compare(tmp_path, [dlv_record(20000.000, 5000.0000, 0.0500, 1),
                             dlv_record(19999.990, 5000.0025, 0.0500, 2)],
                  [20000.100, 20000.090], [20000.000, 19999.990], [0.4, 0.4],
                  registry=REGISTRY + '20000.1\t0.5\t\tshort\n')
    said = ' '.join(messages(rep, check_sync.ERROR))
    assert 'cannot be applied' in said and 'more than one line' in said
    assert 'narrower' not in said and 'wider' not in said


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
    assert 'name line 77' in ' '.join(items(rep, check_sync.ERROR))


def test_the_wavenumber_an_identification_carries_is_ignored():
    """IDEN2 reads only the line number; the wavenumber beside it in
    trans.dat means nothing."""
    rep = check_sync.Report()
    rows = check_sync._dlv_internal(
        [dlv_record(20000.000, 5000.0000, 0.0500, 1)], check_sync.Report(), 12)
    trans = FakeTrans([(4, 3, sync.make_assignment(62, 20000.100, 0.0, 1))])
    check_sync._dlv_vs_trans(rows, trans, rep, 12)
    assert severities(rep) == [check_sync.OK]


def test_an_identified_line_has_the_wavenumber_of_its_dlv_row():
    """Check 9 compares the classification with the line the line number
    names, not with the wavenumber trans.dat happens to carry."""
    trans = FakeTrans([(4, 3, sync.make_assignment(62, 20000.100, 0.0, 1)),
                       (4, 2, sync.make_assignment(62, 30000.000, 0.0, 9))])
    id_of_row = {2: 'B', 3: 'C', 4: 'D'}
    found, _ = check_sync.read_iden_assignments(None, trans, id_of_row,
                                                {1: 20000.000})
    assert found[('C', 'D')] == 20000.000
    # a line number dlv.dat has not got falls back on trans.dat's value
    assert found[('B', 'D')] == 30000.000


def test_identifications_that_agree_pass():
    rep = check_sync.Report()
    rows = check_sync._dlv_internal(
        [dlv_record(20000.000, 5000.0000, 0.0500, 1)], check_sync.Report(), 12)
    trans = FakeTrans([(4, 3, sync.make_assignment(62, 20000.000, 0.0, 1)),
                       (4, 2, IDEN.BLANK_OBS)])
    check_sync._dlv_vs_trans(rows, trans, rep, 12)
    assert severities(rep) == [check_sync.OK]
    assert 'all 1 identifications' in ' '.join(messages(rep, check_sync.OK))


def keyed_iden2(wn_shown, wn_written, key=70349.27026812005,
                wn_measured=70349.1541):
    """An Iden2 of one dlv.dat row, line 1600, written by sync_IDEN2.py for
    the line `key` (dlv_keys.txt) when it measured `wn_measured`
    (dlv_shown.txt); `wn_shown` is what the row shows now."""
    rec = dlv_record(wn_shown, 1e8 / wn_shown, 0.0046, 1600)
    return check_sync.Iden2(
        None, None, {}, [rec],
        {1600: (key, wn_measured, 0.1374, wn_written, 0.2260)},
        {1600: (key, wn_written)})


def test_a_recalibrated_line_is_still_its_line():
    """2026-10-08: a calibration run since the last sync moved 70349.154 to
    70349.067.  Compared by the old number the identification on it read as
    "elsewhere" in all three files; by its wn_key it is the table's line."""
    iden = keyed_iden2(70349.070, 70349.070)
    assert iden.wn_of_line[1600] == pytest.approx(70349.1541)
    cls = pd.DataFrame({'wn_key': [70349.27026812005, 50000.0],
                        'wn_obs': [70349.067, 50000.0]})
    assert iden.follow(cls) == 1
    assert iden.wn_of_line[1600] == pytest.approx(70349.067)
    assert iden.follow(cls) == 0                    # nothing left to rename


def test_a_row_edited_in_iden2_keeps_its_own_wavenumber():
    """A row that no longer shows what sync_IDEN2.py wrote there is a line
    edited on the screen, and its key may no longer name it."""
    iden = keyed_iden2(70349.200, 70349.070)
    cls = pd.DataFrame({'wn_key': [70349.27026812005], 'wn_obs': [70349.067]})
    assert iden.follow(cls) == 0
    assert iden.wn_of_line[1600] == pytest.approx(70349.200)


def test_a_key_the_table_has_not_got_keeps_the_recorded_wavenumber():
    iden = keyed_iden2(70349.070, 70349.070)
    assert iden.follow(pd.DataFrame({'wn_key': [1.0], 'wn_obs': [1.0]})) == 0
    assert iden.follow(None) == 0
    assert iden.wn_of_line[1600] == pytest.approx(70349.1541)


# ---------------------------------------------------------------------------
# The real files
# ---------------------------------------------------------------------------
@real_file
def test_the_real_dlv_has_one_row_per_line_number():
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
