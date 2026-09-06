"""Tests for swap_line_assignments_LOPT.py and swap_line_assignments_IDEN.py.

Every fixture here is built from scratch, so nothing in this file needs the
real line lists or the real IDEN2 files.  The three row builders at the top
reproduce the layouts of the real files byte for byte - each was checked
against a row taken from LOPT_input_lines.txt, IDEN2/enlev.dat and
IDEN2/trans.dat - and they are written out in full rather than assembled
from the modules' own column constants, so that a mistake in those constants
shows up as a failing test instead of being copied into the fixture.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import swap_line_assignments_IDEN as iden
import swap_line_assignments_LOPT as lopt


# ---------------------------------------------------------------------------
# Row builders
# ---------------------------------------------------------------------------
def lopt_row(wn, low, upp, flag='', weight=1.0):
    """One record of the LOPT transitions file, in the sample layout."""
    buf = [' '] * 93
    def put(first, text):
        buf[first - 1:first - 1 + len(text)] = text
    put(1, '%.3f' % wn)
    put(14, '0.005')
    put(25, '100.000')
    put(43, low)
    put(59, upp)
    put(73, flag)
    put(82, '%.4f' % weight)
    put(90, 'cm-1')
    return ''.join(buf).rstrip()


PAR = '\n'.join([
    'lines.txt   ; TRANSITIONS input file name',
    '43     ;1st column of LOWER LEVEL LABEL in the transitions-input file',
    '55     ;LAST column of LOWER LEVEL LABEL in the transitions-input file',
    ' 59    ;1st column of UPPER level LABEL in the transitions-input file',
    '  82    ;1st column of line WEIGHT in the transitions-input file',
]) + '\n'


def en_row(index, e_calc, unc, e_obs, known, j, label):
    """One row of enlev.dat."""
    return ('%4d%12.3f%10.3f%12.3f%2s%12.3f%5.1f /%-11s/'
            % (index, e_calc, unc, e_obs, ' *' if known else '  ',
               e_obs - e_calc, j, label))


def tr_header(index, j, label, e_obs, known):
    """A block header of trans.dat."""
    return ('$%4d    J=%4.1f   %-10s%13.3f%2s '
            % (index, j, label, e_obs, ' *' if known else '  '))


def tr_row(partner, v, e_partner, known, wn, v_obs=0, wn_obs=0.0,
           omc=0.0, dlv=0):
    """One predicted transition of trans.dat."""
    return ('+%4d%5d%12.3f%2s%13.3f%6d%12.3f%11.3f%6d'
            % (partner, v, e_partner, ' *' if known else '  ', wn,
               v_obs, wn_obs, omc, dlv))


def write(path, rows, eol='\r\n'):
    with open(path, 'w', encoding='latin-1', newline='') as fh:
        for row in rows:
            fh.write(row + eol)
    return str(path)


# ---------------------------------------------------------------------------
# A pair of levels to exchange, and their neighbourhood
# ---------------------------------------------------------------------------
ID1, ID2 = '059003.000483', '059003.000398'
N1, N2 = 721, 723
E1, E2 = 117592.990, 117686.460
U1, U2 = 0.170, 0.130
CALC1, CALC2 = 117545.000, 117478.900

# Three lower levels the pair connects to, and their energies.
LOW = {800: 20000.000, 801: 21000.000, 802: 22000.000}


def make_enlev(tmp_path, **kw):
    rows = [
        en_row(N1, CALC1, U1, E1, True, 2.5, ' f25f ~3F2F'),
        en_row(722, 117500.000, 5000.0, 117500.000, False, 2.5, ' f25f ~1D2D'),
        en_row(N2, CALC2, U2, E2, True, 2.5, ' f25f ~3F4D'),
    ]
    for n in sorted(LOW):
        rows.append(en_row(n, LOW[n] - 50.0, 0.004, LOW[n], True, 3.5,
                           ' f25d _3H4H'))
    return write(tmp_path / 'enlev.dat', rows)


def make_trans(tmp_path, drop=(), extra_between=False):
    """A trans.dat in which both levels reach all three lower levels.

    Every transition of N1 carries an identified observed line and none of
    N2's does, which makes it easy to see where the assignments end up.
    ``drop`` names (owner, partner) rows to leave out, so that a level can
    be given a transition its partner does not have.
    """
    rows = []
    for n, e in ((N1, E1), (722, 117500.000), (N2, E2)):
        rows.append(tr_header(n, 2.5, 'f25f ~3F2F', e, n != 722))
        if n == 722:
            continue
        for p in sorted(LOW):
            if (n, p) in drop:
                continue
            wn = e - LOW[p]
            if n == N1:
                rows.append(tr_row(p, 90, LOW[p], True, wn,
                                   v_obs=50, wn_obs=wn + 0.100,
                                   omc=0.100, dlv=1000 + p))
            else:
                rows.append(tr_row(p, 80, LOW[p], True, wn))
        if extra_between and n == N1:
            rows.append(tr_row(N2, 70, E2, True, abs(E1 - E2)))
    for p in sorted(LOW):
        rows.append(tr_header(p, 3.5, 'f25d _3H4H', LOW[p], True))
    return write(tmp_path / 'trans.dat', rows)


def make_levels(tmp_path, e1=E1, e2=E2):
    rows = ['Designation\tEnergy\tID1',
            '%s\t%.4f\t0.005' % (ID1, e1),
            '%s\t%.4f\t0.005' % (ID2, e2)]
    return write(tmp_path / 'levels.txt', rows, eol='\n')


class Args(object):
    """The handful of settings resolve() and crossed() look at."""

    def __init__(self, **kw):
        self.index = kw.get('index')
        self.map = kw.get('map')
        self.levels = kw.get('levels', '')
        self.tol = kw.get('tol', 0.5)
        self.force = kw.get('force', False)


def quiet(*_a, **_k):
    pass


# ---------------------------------------------------------------------------
# Reading and writing without disturbing the line endings
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('module', [lopt, iden])
@pytest.mark.parametrize('text', [
    'a\r\nb\r\n',        # DOS, final newline
    'a\nb\n',            # Unix, final newline
    'a\r\nb',            # DOS, no final newline
    'a\nb\r\nc',         # mixed
    '',                  # empty
])
def test_round_trip_is_byte_exact(tmp_path, module, text):
    path = tmp_path / 'f.txt'
    with open(path, 'w', encoding='latin-1', newline='') as fh:
        fh.write(text)
    records, endings = module.read_records(str(path))
    out = tmp_path / 'g.txt'
    module.write_records(str(out), records, endings)
    assert open(out, 'rb').read() == open(path, 'rb').read()


# ---------------------------------------------------------------------------
# swap_line_assignments_LOPT.py
# ---------------------------------------------------------------------------
def test_columns_come_from_the_parameter_file(tmp_path):
    par = write(tmp_path / 'LOPT.par', PAR.splitlines(), eol='\n')
    assert lopt.read_columns(par) == (43, 55, 59)


def test_columns_missing_from_the_parameter_file(tmp_path):
    par = write(tmp_path / 'p.par', ['x ; nothing useful'], eol='\n')
    with pytest.raises(ValueError, match='does not say where'):
        lopt.read_columns(par)


def test_layout_spans_the_two_identifier_fields():
    layout = lopt.Layout((43, 55, 59))
    assert layout.width == 13
    assert layout.low == (42, 55)
    assert layout.upp == (58, 71)
    row = lopt_row(1000.0, ID1, ID2)
    assert layout.get(row, layout.low) == ID1
    assert layout.get(row, layout.upp) == ID2


def test_an_identifier_too_wide_for_its_field_is_refused():
    layout = lopt.Layout((43, 55, 59))
    with pytest.raises(ValueError, match='does not fit'):
        layout.put(lopt_row(1000.0, ID1, ID2), layout.low, 'x' * 14)


def test_swap_exchanges_both_ends_and_leaves_the_rest_alone():
    layout = lopt.Layout((43, 55, 59))
    rows = [lopt_row(1000.0, '059003.000063', ID1),
            lopt_row(900.0, ID2, '059003.000600', flag='P', weight=0.0),
            lopt_row(800.0, '059003.000063', '059003.000064')]
    out, info = lopt.swap_records(rows, layout, ID1, ID2)
    assert layout.get(out[0], layout.upp) == ID2
    assert layout.get(out[1], layout.low) == ID1
    assert out[2] == rows[2]
    # Nothing outside the two fields moved.
    for old, new in zip(rows, out):
        assert len(new.rstrip()) == len(old.rstrip())
        assert new[:42] == old[:42]
        assert new[71:] == old[71:]
    assert info['counts'][(ID1, 'upper')] == 1
    assert info['counts'][(ID2, 'lower')] == 1
    assert len(info['moved']) == 2


def test_an_absent_identifier_stops_the_swap():
    layout = lopt.Layout((43, 55, 59))
    rows = [lopt_row(1000.0, '059003.000063', ID1)]
    _, info = lopt.swap_records(rows, layout, ID1, ID2)
    with pytest.raises(ValueError, match='does not occur'):
        lopt.check(info, ID1, ID2, 'lines.txt')


def test_a_line_joining_the_two_levels_stops_the_swap():
    layout = lopt.Layout((43, 55, 59))
    rows = [lopt_row(93.5, ID1, ID2)]
    _, info = lopt.swap_records(rows, layout, ID1, ID2)
    with pytest.raises(ValueError, match='level with itself'):
        lopt.check(info, ID1, ID2, 'lines.txt')


def test_a_duplicate_made_by_the_swap_stops_it():
    layout = lopt.Layout((43, 55, 59))
    rows = [lopt_row(1000.0, '059003.000063', ID1),
            lopt_row(1000.0, '059003.000063', ID2),
            lopt_row(1000.0, '059003.000063', ID1)]
    _, info = lopt.swap_records(rows, layout, ID1, ID2)
    with pytest.raises(ValueError, match='counted more than once'):
        lopt.check(info, ID1, ID2, 'lines.txt')


def test_fixed_levels_are_exchanged_when_both_are_fixed(tmp_path):
    path = write(tmp_path / 'fix.txt',
                 ['%s       1.000        0.010' % ID1,
                  '%s       2.000        0.020' % ID2], eol='\n')
    records, _, note = lopt.swap_fixed_levels(path, ID1, ID2)
    assert 'exchanged' in note
    assert records[0].split() == [ID1, '2.000', '0.020']
    assert records[1].split() == [ID2, '1.000', '0.010']


def test_one_fixed_level_of_the_two_stops_the_swap(tmp_path):
    path = write(tmp_path / 'fix.txt', ['%s   0   0' % ID1], eol='\n')
    with pytest.raises(ValueError, match='pinning the wrong level'):
        lopt.swap_fixed_levels(path, ID1, ID2)


def test_no_fixed_level_of_the_two_is_left_alone(tmp_path):
    path = write(tmp_path / 'fix.txt', ['059003.000001   0   0'], eol='\n')
    records, _, note = lopt.swap_fixed_levels(path, ID1, ID2)
    assert 'nothing to do' in note
    assert records == ['059003.000001   0   0']


def test_the_ledger_is_reported_with_the_replacement_identifier(tmp_path):
    path = write(tmp_path / 'ledger.csv',
                 ['wn_obs,low_id,upp_id,decision,date,reason',
                  '94058.9743,059003.000095,%s,accept,9/3/2026,by hand' % ID1,
                  '100.0,059003.000001,059003.000002,reject,9/3/2026,no'],
                 eol='\n')
    rows = lopt.ledger_rows(path, ID1, ID2)
    assert len(rows) == 1
    n, wn, low, low_new, upp, upp_new, decision, reason = rows[0]
    assert (n, low, low_new, upp, upp_new) == (2, '059003.000095',
                                               '059003.000095', ID1, ID2)
    assert decision == 'accept'


def test_lopt_end_to_end_keeps_every_other_byte(tmp_path):
    lines = write(tmp_path / 'LOPT_input_lines.txt',
                  [lopt_row(1000.0, '059003.000063', ID1),
                   lopt_row(900.0, '059003.000064', ID2),
                   lopt_row(800.0, '059003.000063', '059003.000064')])
    par = write(tmp_path / 'LOPT.par', PAR.splitlines(), eol='\n')
    out = str(tmp_path / 'out.txt')
    code = lopt.main([ID1, ID2, '--lines', lines, '--par', par,
                      '--no-fixlev', '--ledger', str(tmp_path / 'none.csv'),
                      '--out', out])
    assert code == 0
    before = open(lines, 'rb').read().split(b'\r\n')
    after = open(out, 'rb').read().split(b'\r\n')
    assert len(before) == len(after)
    assert [len(x) for x in before] == [len(x) for x in after]
    assert sum(1 for a, b in zip(before, after) if a != b) == 2


def test_lopt_dry_run_writes_nothing(tmp_path):
    lines = write(tmp_path / 'lines.txt',
                  [lopt_row(1000.0, '059003.000063', ID1),
                   lopt_row(900.0, '059003.000064', ID2)])
    before = open(lines, 'rb').read()
    lopt.main([ID1, ID2, '--lines', lines, '--par', str(tmp_path / 'x.par'),
               '--no-fixlev', '--ledger', str(tmp_path / 'none.csv'),
               '--dry-run'])
    assert open(lines, 'rb').read() == before
    assert not os.path.exists(lines + '.bak')


# ---------------------------------------------------------------------------
# swap_line_assignments_IDEN.py - reading the files
# ---------------------------------------------------------------------------
def test_enlev_fields_are_read_at_the_right_columns(tmp_path):
    en = iden.Enlev(make_enlev(tmp_path))
    assert en.e_calc(N1) == pytest.approx(CALC1)
    assert en.e_obs(N1) == pytest.approx(E1)
    assert en.j(N1) == '2.5'
    assert en.label(N1) == 'f25f ~3F2F'
    assert en.known(N1) is True
    assert en.known(722) is False


def test_setting_a_measurement_keeps_the_row_width_and_redoes_o_minus_c(tmp_path):
    en = iden.Enlev(make_enlev(tmp_path))
    before = en.record(N1)
    en.set_measurement(N1, U2, E2, True)
    after = en.record(N1)
    assert len(after) == len(before)
    assert after[:16] == before[:16]          # the calculated energy is untouched
    assert after[52:] == before[52:]          # so are J and the label
    assert en.e_obs(N1) == pytest.approx(E2)
    assert iden.number(after, iden.EN_OMC) == pytest.approx(E2 - CALC1, abs=1e-3)


def test_level_energies_read_from_both_table_shapes(tmp_path):
    tab = make_levels(tmp_path)
    assert iden.read_level_energies(tab)[ID1] == pytest.approx(E1)
    csv_path = write(tmp_path / 'e.csv',
                     ['level_id,E_input,comment', '%s,123.5,x' % ID1],
                     eol='\n')
    assert iden.read_level_energies(csv_path)[ID1] == pytest.approx(123.5)


def test_a_table_without_the_two_columns_is_refused(tmp_path):
    path = write(tmp_path / 'e.csv', ['a,b', '1,2'], eol='\n')
    with pytest.raises(ValueError, match='no identifier and energy columns'):
        iden.read_level_energies(path)


def test_lookup_table_with_and_without_a_header(tmp_path):
    with_header = write(tmp_path / 'm1.csv',
                        ['level_id,iden2_number', '%s,721' % ID1], eol='\n')
    assert iden.read_map(with_header) == {ID1: N1}
    bare = write(tmp_path / 'm2.csv', ['%s,721' % ID1, '%s,723' % ID2],
                 eol='\n')
    assert iden.read_map(bare) == {ID1: N1, ID2: N2}


def test_resolving_by_energy(tmp_path):
    en = iden.Enlev(make_enlev(tmp_path))
    args = Args(levels=make_levels(tmp_path))
    assert iden.resolve((ID1, ID2), args, en, quiet) == (N1, N2)


def test_resolving_by_energy_finds_nothing(tmp_path):
    en = iden.Enlev(make_enlev(tmp_path))
    args = Args(levels=make_levels(tmp_path, e1=1.0))
    with pytest.raises(ValueError, match='within 0.5'):
        iden.resolve((ID1, ID2), args, en, quiet)


def test_resolving_by_energy_is_ambiguous(tmp_path):
    en = iden.Enlev(make_enlev(tmp_path))
    args = Args(levels=make_levels(tmp_path), tol=200.0)
    with pytest.raises(ValueError, match='name the right one'):
        iden.resolve((ID1, ID2), args, en, quiet)


def test_two_row_numbers_need_no_lookup(tmp_path):
    """A bare run of digits is the level's row number in enlev.dat, which is
    what IDEN2 itself calls the level; nothing has to be translated."""
    en = iden.Enlev(make_enlev(tmp_path))
    assert iden.resolve((str(N1), str(N2)), Args(), en, quiet) == (N1, N2)


def test_a_row_number_and_an_identifier_may_be_mixed(tmp_path):
    """Each argument is taken on its own: the number is already the answer,
    the identifier is looked up by energy."""
    en = iden.Enlev(make_enlev(tmp_path))
    args = Args(levels=make_levels(tmp_path))
    assert iden.resolve((str(N1), ID2), args, en, quiet) == (N1, N2)
    assert iden.resolve((ID1, str(N2)), args, en, quiet) == (N1, N2)


def test_a_row_number_needs_no_level_table(tmp_path):
    """The energy route is what wants LOPT's output; naming the rows by
    number avoids it, and so avoids the question of which state that table
    is in."""
    en = iden.Enlev(make_enlev(tmp_path))
    assert iden.resolve((str(N2), str(N1)), Args(levels='no_such_file'),
                        en, quiet) == (N2, N1)


def test_a_mixed_pair_still_needs_the_table_for_its_identifier(tmp_path):
    en = iden.Enlev(make_enlev(tmp_path))
    with pytest.raises(ValueError, match='no level-energy table'):
        iden.resolve((str(N1), ID2), Args(levels='no_such_file'), en, quiet)


def test_index_and_map_take_precedence(tmp_path):
    en = iden.Enlev(make_enlev(tmp_path))
    assert iden.resolve((ID1, ID2), Args(index=[N2, N1]), en, quiet) == (N2, N1)
    m = write(tmp_path / 'm.csv', ['%s,723' % ID1, '%s,721' % ID2], eol='\n')
    assert iden.resolve((ID1, ID2), Args(map=m), en, quiet) == (N2, N1)


def test_trans_is_indexed_by_owner_and_partner(tmp_path):
    tr = iden.Trans(make_trans(tmp_path))
    assert set(tr.header_of) == {N1, 722, N2} | set(LOW)
    assert tr.partners_of(N1) == set(LOW)
    assert tr.row(N1, 800) == tr.row(800, N1)
    assert tr.row(N1, 999) is None


# ---------------------------------------------------------------------------
# swap_line_assignments_IDEN.py - the exchange
# ---------------------------------------------------------------------------
def test_the_exchange_moves_energies_lines_and_nothing_else(tmp_path):
    en = iden.Enlev(make_enlev(tmp_path))
    tr = iden.Trans(make_trans(tmp_path))
    before = {k: v for k, v in enumerate(tr.records)}
    info = iden.swap(en, tr, N1, N2, quiet)

    # The measurements changed places, the calculation did not.
    assert en.e_obs(N1) == pytest.approx(E2)
    assert en.e_obs(N2) == pytest.approx(E1)
    assert iden.number(en.record(N1), iden.EN_UNC) == pytest.approx(U2)
    assert en.e_calc(N1) == pytest.approx(CALC1)
    assert en.label(N1) == 'f25f ~3F2F'

    # All three lines moved from N1 to N2 and none the other way.
    assert len(info['moved']) == 3
    assert {m['to'] for m in info['moved']} == {N2}
    assert not info['orphans'] and not info['mismatch']

    for p in LOW:
        row = tr.records[tr.row(N2, p)]
        assert iden.obs_row(iden.assignment(row)) == 1000 + p
        assert iden.obs_omc(iden.assignment(row)) == pytest.approx(0.100)
        # The predicted wavenumber follows the energy that came with them.
        assert iden.number(row, iden.TR_WN) == pytest.approx(E1 - LOW[p],
                                                             abs=1e-3)
        assert not iden.has_line(iden.assignment(tr.records[tr.row(N1, p)]))

    # Every row keeps its width; only the rows naming the pair changed.
    assert all(len(tr.records[k]) == len(before[k]) for k in before)
    changed = [k for k in before if tr.records[k] != before[k]]
    assert len(changed) == info['rows'] + 2   # the two block headers as well


def test_the_block_headers_carry_the_new_energies(tmp_path):
    en = iden.Enlev(make_enlev(tmp_path))
    tr = iden.Trans(make_trans(tmp_path))
    iden.swap(en, tr, N1, N2, quiet)
    head = tr.records[tr.header_of[N1]]
    assert iden.number(head, iden.TH_EOBS) == pytest.approx(E2, abs=1e-3)
    assert len(head) == 44


def test_a_line_with_nowhere_to_go_is_reported_and_removed(tmp_path):
    """N1 reaches level 802 and N2 does not, so that line cannot move."""
    en = iden.Enlev(make_enlev(tmp_path))
    tr = iden.Trans(make_trans(tmp_path, drop=((N2, 802),)))
    info = iden.swap(en, tr, N1, N2, quiet)
    assert len(info['moved']) == 2
    assert len(info['orphans']) == 1
    orphan = info['orphans'][0]
    assert (orphan['from'], orphan['to'], orphan['partner']) == (N1, N2, 802)
    assert orphan['dlv'] == 1802
    assert not iden.has_line(iden.assignment(tr.records[tr.row(N1, 802)]))


def test_a_transition_between_the_two_levels_stops_the_exchange(tmp_path):
    en = iden.Enlev(make_enlev(tmp_path))
    tr = iden.Trans(make_trans(tmp_path, extra_between=True))
    with pytest.raises(ValueError, match='level with itself'):
        iden.swap(en, tr, N1, N2, quiet)


def test_o_minus_c_of_every_transferred_line_is_preserved(tmp_path):
    en = iden.Enlev(make_enlev(tmp_path))
    tr = iden.Trans(make_trans(tmp_path))
    info = iden.swap(en, tr, N1, N2, quiet)
    assert [m['omc'] for m in info['moved']] == [pytest.approx(0.100)] * 3
    for p in LOW:
        row = tr.records[tr.row(N2, p)]
        obs = iden.assignment(row)
        assert (iden.obs_wavenumber(obs)
                - iden.number(row, iden.TR_WN)) == pytest.approx(0.100,
                                                                 abs=1e-3)


def test_a_disturbed_file_is_caught_by_the_o_minus_c_check(tmp_path):
    """An assignment whose O-C does not match its own row cannot be moved
    safely, and the mismatch is what says so."""
    en = iden.Enlev(make_enlev(tmp_path))
    tr = iden.Trans(make_trans(tmp_path))
    k = tr.row(N1, 800)
    row = tr.records[k]
    tr.records[k] = row[:iden.TR_OBS] + iden.set_omc(iden.assignment(row), 9.0)
    info = iden.swap(en, tr, N1, N2, quiet)
    assert len(info['mismatch']) == 1
    assert info['mismatch'][0]['was'] == pytest.approx(9.0)


def test_the_crossed_guard_sees_an_exchange_already_made(tmp_path):
    en = iden.Enlev(make_enlev(tmp_path))
    args = Args(index=[N1, N2], levels=make_levels(tmp_path, e1=E2, e2=E1))
    with pytest.raises(ValueError, match='already'):
        iden.crossed(en, (N1, N2), (ID1, ID2), args, quiet)


def test_the_crossed_guard_passes_a_file_not_yet_exchanged(tmp_path):
    en = iden.Enlev(make_enlev(tmp_path))
    args = Args(index=[N1, N2], levels=make_levels(tmp_path))
    iden.crossed(en, (N1, N2), (ID1, ID2), args, quiet)


# ---------------------------------------------------------------------------
# swap_line_assignments_IDEN.py - the command line
# ---------------------------------------------------------------------------
def test_iden_end_to_end_writes_both_files(tmp_path):
    make_enlev(tmp_path)
    make_trans(tmp_path)
    levels = make_levels(tmp_path)
    out = tmp_path / 'out'
    code = iden.main([ID1, ID2, '--iden2-dir', str(tmp_path),
                      '--levels', levels, '--out-dir', str(out)])
    assert code == 0
    for name in ('enlev.dat', 'trans.dat'):
        before = open(tmp_path / name, 'rb').read().split(b'\r\n')
        after = open(out / name, 'rb').read().split(b'\r\n')
        assert [len(x) for x in before] == [len(x) for x in after]
    en = iden.Enlev(str(out / 'enlev.dat'))
    assert en.e_obs(N1) == pytest.approx(E2)


def test_iden_dry_run_writes_nothing(tmp_path):
    path = make_enlev(tmp_path)
    make_trans(tmp_path)
    before = open(path, 'rb').read()
    iden.main([ID1, ID2, '--iden2-dir', str(tmp_path),
               '--levels', make_levels(tmp_path), '--dry-run'])
    assert open(path, 'rb').read() == before
    assert not os.path.exists(path + '.bak')


def test_iden_refuses_two_levels_of_different_j(tmp_path):
    make_enlev(tmp_path)
    make_trans(tmp_path)
    levels = write(tmp_path / 'lv.txt',
                   ['Designation\tEnergy',
                    '%s\t%.4f' % (ID1, E1),
                    '%s\t%.4f' % (ID2, LOW[800])], eol='\n')
    with pytest.raises(ValueError, match='different J'):
        iden.main([ID1, ID2, '--iden2-dir', str(tmp_path),
                   '--levels', levels, '--dry-run'])


def test_iden_refuses_a_level_never_observed(tmp_path):
    make_enlev(tmp_path)
    make_trans(tmp_path)
    with pytest.raises(ValueError, match='must both have an observed energy'):
        iden.main([ID1, ID2, '--iden2-dir', str(tmp_path),
                   '--levels', make_levels(tmp_path),
                   '--index', str(N1), '722', '--dry-run'])


def test_iden_refuses_the_same_identifier_twice(tmp_path):
    make_enlev(tmp_path)
    assert iden.main([ID1, ID1, '--iden2-dir', str(tmp_path)]) == 2
