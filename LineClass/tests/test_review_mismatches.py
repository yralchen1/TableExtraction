"""Tests for `review_mismatches`.

Run from the LineClass directory:  python -m pytest tests -q

The tool compares two sets of LOPT files, one transition at a time, and builds
a worksheet of the places where they disagree about whether an assignment is
accepted.  Then, on a second run, it carries out whatever verdict has been
written into that worksheet, in the four files a verdict has to reach: the
shared ledger, the inflated-uncertainty registry, and IDEN2's own `dlv.dat`
and `trans.dat`.

Most of these tests build a miniature set of files in a temporary directory -
two transitions, three observed lines - and read the worksheet or the rewritten
files back.  The point of doing it that small is that every number in the
expected result can be worked out by hand.
"""
import io
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import make_LOPT_input                     # noqa: E402
import review_mismatches as rm             # noqa: E402
import sync_IDEN2 as sync                  # noqa: E402
import swap_line_assignments_IDEN as IDEN  # noqa: E402


# ---------------------------------------------------------------------------
# the fixed-column records the tool reads
# ---------------------------------------------------------------------------
def lopt_record(wn, unc, intens, low, upp, flag, weight):
    """One record of LOPT's transitions file.

    Built by the writer the pipeline itself uses, so that a column moving in
    one place and not the other cannot make these tests pass.
    """
    return make_LOPT_input.format_line(wn, unc, intens, low, upp, flag, weight)


def dlv_record(wn, lam, u_lam, row, code=62, character='c'):
    """One row of dlv.dat at its fixed width of 66 characters."""
    rec = ('%5d%14.3f%14.4f  /%10s/%13.4f%6d'
           % (code, wn, lam, character, u_lam, row))
    assert len(rec) == sync.DLV_WIDTH, len(rec)
    return rec


def trans_header(index, energy):
    rec = '$%4d    J= 1.5   f5d6d~3D4D %13.3f   ' % (index, energy)
    return rec


def trans_row(partner, code, e_partner, rwn, assignment=None):
    rec = sync.transition_row(partner, code, e_partner, True, rwn,
                              assignment or IDEN.BLANK_OBS)
    assert len(rec) == IDEN.ROW_WIDTH if hasattr(IDEN, 'ROW_WIDTH') else True
    return rec


CLS_HEADER = ('wn_obs,wn_key,unc_wn_obs,obs_intens,char,low_id,upp_id,'
              'calc_intens,dif_wn_O-C,grade,notes1,notes2,manual,new,'
              'accepted,n_accepted,rwn\n')


def cls_row(wn_obs, wn_key, unc, intens, char, low, upp, icalc, omc,
            accepted, n_accepted, grade='5A', new='1'):
    return ('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,,,,%s,%s,%s,0.0\n'
            % (wn_obs, wn_key, unc, intens, char, low, upp, icalc, omc,
               grade, new, accepted, n_accepted))


LOW_A, UPP_A = '059003.000104', '059003.000224'
LOW_B, UPP_B = '059003.000199', '059003.000327'
IDEN_A, IDEN_B = 11, 22
IDEN_LOW_A, IDEN_LOW_B = 33, 44


@pytest.fixture
def miniature(tmp_path):
    """A baseline directory and an `iter` set inside it, both complete.

    One observed line at 41473.8 carries two candidate transitions.  The
    baseline accepts both - it is a blend there - and the set accepts only the
    second, so the first is the disagreement the tool is about.
    """
    base = tmp_path / 'LineClass'
    iter_dir = base / 'iter'
    iden2 = iter_dir / 'IDEN2'
    iden2.mkdir(parents=True)

    def write(path, text, crlf=False):
        newline = '\r\n' if crlf else '\n'
        with io.open(str(path), 'w', encoding='utf-8', newline='') as fh:
            fh.write(text.replace('\n', newline))

    # the baseline's LOPT input: both components accepted
    write(base / 'LOPT_input_lines.txt',
          lopt_record(41473.799, 0.099, 16761.0, LOW_A, UPP_A, '', 0.5) + '\n'
          + lopt_record(41473.799, 0.099, 16761.0, LOW_B, UPP_B, '', 0.5)
          + '\n')
    # the set's: the first component is flagged P, the second is not
    write(iter_dir / 'LOPT_input_lines.txt',
          lopt_record(41473.820, 0.050, 16761.0, LOW_A, UPP_A, 'P', 0.0) + '\n'
          + lopt_record(41473.820, 0.050, 16761.0, LOW_B, UPP_B, '', 1.0)
          + '\n')

    write(iter_dir / 'line_classifications.csv',
          CLS_HEADER
          + cls_row('41473.8197', '41473.79906049907', '0.0499', '16761.3',
                    'w', LOW_A, UPP_A, '17122.3', '-0.1600', '0.0', '1')
          + cls_row('41473.8197', '41473.79906049907', '0.0499', '16761.3',
                    'w', LOW_B, UPP_B, '16267.6', '0.3836', '1.0', '1'))

    write(base / 'line_decisions.csv',
          'wn_key,low_id,upp_id,decision,date,reason\n')
    write(base / 'inflated_unc_lines.txt', 'wn_key\tunc_wn\tdate\treason\n')

    write(iden2 / 'IDEN_level_ids.txt',
          'level_id\tIDEN_id\n'
          + '%s\t%d\n' % (LOW_A, IDEN_LOW_A)
          + '%s\t%d\n' % (UPP_A, IDEN_A)
          + '%s\t%d\n' % (LOW_B, IDEN_LOW_B)
          + '%s\t%d\n' % (UPP_B, IDEN_B))

    # dlv.dat: three rows, numbered by position, decreasing in wavenumber
    write(iden2 / 'dlv.dat',
          dlv_record(50000.000, 1999.9200, 0.0100, 1) + '\n'
          + dlv_record(41473.820, 2411.1700, 0.0029, 2) + '\n'
          + dlv_record(30000.000, 3332.8000, 0.0100, 3) + '\n', crlf=True)

    # trans.dat: a block for each upper level, neither assignment made
    write(iden2 / 'trans.dat',
          trans_header(IDEN_A, 102830.460) + '\n'
          + trans_row(IDEN_LOW_A, 97, 61357.018, 41473.980) + '\n'
          + trans_header(IDEN_B, 102830.460) + '\n'
          + trans_row(IDEN_LOW_B, 96, 61357.018, 41473.436) + '\n', crlf=True)

    return base


@pytest.fixture
def paths(miniature, monkeypatch):
    monkeypatch.setattr(rm, 'HERE', str(miniature))
    import swap_paths
    monkeypatch.setattr(swap_paths, 'PROJECT', str(miniature))
    return rm.resolve('iter')


# ---------------------------------------------------------------------------
# the arithmetic
# ---------------------------------------------------------------------------
def test_center_of_gravity_weights_by_predicted_intensity():
    rows = [{'calc_intens': '3000', 'dif_wn_O-C': '0.4'},
            {'calc_intens': '1000', 'dif_wn_O-C': '-0.4'}]
    assert rm.center_of_gravity(rows) == pytest.approx(0.2)


def test_center_of_gravity_is_none_without_predicted_intensity():
    assert rm.center_of_gravity([{'calc_intens': '0',
                                  'dif_wn_O-C': '0.4'}]) is None
    assert rm.center_of_gravity([]) is None


def test_code_is_ten_times_the_natural_log():
    assert rm.code_of(1.0) == 0
    assert rm.code_of(100.0) == 46          # 10*ln(100) = 46.05
    assert rm.code_of(0.0) is None


def test_twenty_code_units_is_a_factor_of_e_squared():
    import math
    assert rm.code_of(100.0 * math.e ** 2) - rm.code_of(100.0) == 20


def test_preferred_uncertainty_prefers_the_half_wavenumber_value():
    # within a factor of 1.5 of 0.5, the standing far-ultraviolet value
    assert rm.preferred_uncertainty(0.42) == 0.5
    assert rm.preferred_uncertainty(-0.70) == 0.5
    # outside it, the center of gravity itself
    assert rm.preferred_uncertainty(0.105) == pytest.approx(0.105)
    assert rm.preferred_uncertainty(0.9) == pytest.approx(0.9)


def test_much_fainter_than_the_line_holder_stays_rejected():
    verdict, why = rm.suggest({'code_behind': 41, 'cg_with': 0.1,
                               'n_accepted': 1, 'n_candidates': 2})
    assert verdict == 'stays rejected'
    assert 'fainter' in why


def test_a_center_of_gravity_beyond_the_limit_stays_rejected():
    verdict, why = rm.suggest({'code_behind': 2, 'cg_with': -0.9,
                               'n_accepted': 1, 'n_candidates': 2})
    assert verdict == 'stays rejected'
    assert 'center of gravity' in why


def test_nothing_to_blend_with_stays_rejected():
    verdict, why = rm.suggest({'code_behind': None, 'cg_with': -0.6,
                               'n_accepted': 0, 'n_candidates': 1})
    assert verdict == 'stays rejected'
    assert 'no blend' in why


def test_everything_else_is_left_for_the_analyst():
    verdict, why = rm.suggest({'code_behind': 5, 'cg_with': 0.1,
                               'n_accepted': 1, 'n_candidates': 2})
    assert verdict == 'review'
    assert 'IDEN2' in why


# ---------------------------------------------------------------------------
# reading
# ---------------------------------------------------------------------------
def test_a_transition_named_twice_in_one_lopt_file_is_an_error(tmp_path):
    p = tmp_path / 'LOPT_input_lines.txt'
    with io.open(str(p), 'w', encoding='utf-8', newline='') as fh:
        fh.write(lopt_record(100.0, 0.1, 1.0, LOW_A, UPP_A, '', 1.0) + '\n')
        fh.write(lopt_record(200.0, 0.1, 1.0, LOW_A, UPP_A, '', 1.0) + '\n')
    with pytest.raises(rm.ReviewError) as exc:
        rm.read_lopt(str(p))
    assert 'twice' in str(exc.value)


def test_the_id_table_must_carry_its_own_header(tmp_path):
    p = tmp_path / 'IDEN_level_ids.txt'
    with io.open(str(p), 'w', encoding='utf-8', newline='') as fh:
        fh.write('%s\t7\n' % LOW_A)
    with pytest.raises(rm.ReviewError):
        rm.read_ids(str(p))


def test_the_nearest_dlv_row_is_found_and_only_within_the_tolerance():
    wavenumbers, numbers = [30000.0, 41473.82, 50000.0], [3, 2, 1]
    assert rm.dlv_row_of(wavenumbers, numbers, 41473.8195) == 2
    assert rm.dlv_row_of(wavenumbers, numbers, 41473.9) is None


# ---------------------------------------------------------------------------
# the worksheet
# ---------------------------------------------------------------------------
def test_the_disagreement_becomes_one_row_and_the_agreement_none(paths):
    cases, decided = rm.build_cases(paths, 0, False)
    assert decided == 0
    assert [(c['level_id1'], c['level_id2']) for c in cases] \
        == [(LOW_A, UPP_A)]


def test_the_worksheet_carries_the_blend_arithmetic(paths):
    case, = rm.build_cases(paths, 0, False)[0]
    assert case['lc_flag'] == '(blank)'
    assert case['iter_flag'] == 'P'
    assert case['iden2_id1'] == IDEN_LOW_A
    assert case['iden2_id2'] == IDEN_A
    # the accepted component alone stands +0.3836 from Ritz; with this one
    # added, weighted by 16267.6 and 17122.3, it stands +0.1049
    assert case['cg_without'] == '+0.3836'
    assert case['cg_with'] == '+0.1048'
    assert case['cg_gain'] == '+0.2788'
    assert case['cg_over_unc'] == '2.10'
    # and because that is wider than the line's own 0.0499, an inflated
    # uncertainty is proposed
    assert case['unc_wn_inflated'] == '0.1048'
    assert case['decision'] == ''


def test_the_line_holder_in_trans_dat_is_reported(paths):
    # nothing is assigned in the miniature trans.dat
    case, = rm.build_cases(paths, 0, False)[0]
    assert case['held_by'] == '(unassigned)'
    assert case['held_code'] == ''
    assert case['code_behind'] == ''


def test_a_disagreement_the_ledger_rules_on_is_left_out(paths):
    with io.open(paths['decisions'], 'a', encoding='utf-8', newline='') as fh:
        fh.write('41473.79906049907,%s,%s,accept,9/26/2026,tested\n'
                 % (LOW_A, UPP_A))
    cases, decided = rm.build_cases(paths, 0, False)
    assert cases == [] and decided == 1
    kept, decided = rm.build_cases(paths, 0, True)
    assert decided == 1
    assert 'the ledger already rules on this' in kept[0]['Notes']


def test_from_iden2_id_skips_the_levels_already_reviewed(paths):
    assert rm.build_cases(paths, IDEN_A, False)[0]        # 11, kept
    assert rm.build_cases(paths, IDEN_A + 1, False)[0] == []


def test_a_level_missing_from_the_id_table_is_an_error(paths):
    text = io.open(paths['ids'], encoding='utf-8', newline='').read()
    kept = [ln for ln in text.splitlines() if LOW_A not in ln]
    with io.open(paths['ids'], 'w', encoding='utf-8', newline='') as fh:
        fh.write('\n'.join(kept) + '\n')
    with pytest.raises(rm.ReviewError) as exc:
        rm.build_cases(paths, 0, False)
    assert 'the lookup is broken' in str(exc.value)


def test_the_worksheet_round_trips(paths, tmp_path):
    cases, _ = rm.build_cases(paths, 0, False)
    out = str(tmp_path / 'decisions.txt')
    rm.write_worksheet(out, cases)
    back = rm.read_worksheet(out)
    assert len(back) == 1
    assert back[0]['level_id1'] == LOW_A
    assert back[0]['cg_with'] == '+0.1048'


# ---------------------------------------------------------------------------
# the free lines
# ---------------------------------------------------------------------------
def test_a_free_line_off_its_ritz_value_is_not_offered(paths):
    # the rejected component stands 0.1600 from Ritz on an uncertainty of
    # 0.0499, which is 3.2 sigma
    assert rm.build_free_lines(paths, 0, 2.0, 0.04) == []


def test_a_free_line_within_the_window_is_offered(paths):
    case, = rm.build_free_lines(paths, 0, 5.0, 0.04)
    assert case['level_id1'] == LOW_A
    assert case['iden2_id2'] == IDEN_A
    assert case['held_by'] == '(unassigned)'
    assert case['cg_over_unc'] == '3.21'
    assert 'the line is free' in case['why']
    assert case['decision'] == ''


def test_a_line_another_assignment_holds_is_not_offered(paths):
    rm.update_trans(paths['trans'], paths['dlv'],
                    [(LOW_B, UPP_B, 41473.8197)],
                    rm.read_ids(paths['ids']), lambda *a: None, False)
    assert rm.build_free_lines(paths, 0, 5.0, 0.04) == []


def test_a_transition_predicted_too_faint_is_not_offered(paths):
    # the component carries 17122 against the line's 16761, so a floor of
    # twice the observed intensity rules it out
    assert rm.build_free_lines(paths, 0, 5.0, 2.0) == []
    assert rm.build_free_lines(paths, 0, 5.0, 1.0) != []


def test_a_free_line_the_ledger_rules_on_is_not_offered(paths):
    with io.open(paths['decisions'], 'a', encoding='utf-8', newline='') as fh:
        fh.write('41473.79906049907,%s,%s,reject,9/26/2026,tested\n'
                 % (LOW_A, UPP_A))
    assert rm.build_free_lines(paths, 0, 5.0, 0.04) == []


def test_a_free_line_below_the_starting_level_is_not_offered(paths):
    assert rm.build_free_lines(paths, IDEN_A + 1, 5.0, 0.04) == []


def test_a_free_line_can_be_applied_like_any_other_row(paths, tmp_path):
    cases = rm.build_free_lines(paths, 0, 5.0, 0.04)
    cases[0]['decision'] = 'accepted'
    cases[0]['Notes'] = 'recovered at the level in IDEN2'
    out = str(tmp_path / 'free.txt')
    rm.write_worksheet(out, cases)
    rm.apply_worksheet(paths, out, lambda *a: None, False)
    ledger = io.open(paths['decisions'], encoding='utf-8', newline='').read()
    assert ',accept,' in ledger
    trans = IDEN.Trans(paths['trans'])
    tail = IDEN.assignment(trans.records[trans.row(IDEN_A, IDEN_LOW_A)])
    assert IDEN.obs_row(tail) == 2


# ---------------------------------------------------------------------------
# carrying a decision out
# ---------------------------------------------------------------------------
def worksheet_with(paths, tmp_path, decision, notes='tested', unc=None):
    cases, _ = rm.build_cases(paths, 0, False)
    cases[0]['decision'] = decision
    cases[0]['Notes'] = notes
    if unc is not None:
        cases[0]['unc_wn_inflated'] = unc
    out = str(tmp_path / 'decisions.txt')
    rm.write_worksheet(out, cases)
    return out


def test_an_unknown_word_in_the_decision_column_stops_the_run(paths, tmp_path):
    out = worksheet_with(paths, tmp_path, 'maybe')
    with pytest.raises(rm.ReviewError) as exc:
        rm.apply_worksheet(paths, out, lambda *a: None, True)
    assert "'maybe'" in str(exc.value)


def test_an_empty_decision_changes_nothing(paths, tmp_path):
    out = worksheet_with(paths, tmp_path, '')
    before = io.open(paths['decisions'], encoding='utf-8', newline='').read()
    rm.apply_worksheet(paths, out, lambda *a: None, False)
    assert io.open(paths['decisions'], encoding='utf-8',
                   newline='').read() == before


def test_accepting_writes_the_ledger_row_and_the_registry(paths, tmp_path):
    out = worksheet_with(paths, tmp_path, 'accepted, inflated',
                         notes='blend restored')
    rm.apply_worksheet(paths, out, lambda *a: None, False)
    ledger = io.open(paths['decisions'], encoding='utf-8', newline='').read()
    assert '41473.79906049907,%s,%s,accept' % (LOW_A, UPP_A) in ledger
    assert 'blend restored' in ledger
    registry = io.open(paths['inflated'], encoding='utf-8', newline='').read()
    assert '41473.7991\t0.1048' in registry


def test_rejecting_writes_a_reject_row_and_no_assignment(paths, tmp_path):
    out = worksheet_with(paths, tmp_path, 'rejected', unc='')
    rm.apply_worksheet(paths, out, lambda *a: None, False)
    ledger = io.open(paths['decisions'], encoding='utf-8', newline='').read()
    assert ',reject,' in ledger
    trans = io.open(paths['trans'], encoding='utf-8', newline='').read()
    assert IDEN.BLANK_OBS in trans


def test_the_dlv_uncertainty_is_converted_with_the_rows_own_wavelength(paths):
    # row 2 is at 41473.820 cm-1 and 2411.1700 A, so 0.1049 cm-1 is
    # 0.1049 * 2411.17 / 41473.820 = 0.0061 A
    changed = rm.update_dlv(paths['dlv'], [(41473.8197, 0.1049)],
                            lambda *a: None, False)
    assert len(changed) == 1
    row, wn, u_old, u_lam, u_wn = changed[0]
    assert row == 2
    assert u_lam == pytest.approx(0.1049 * 2411.17 / 41473.820, abs=5e-5)
    text = io.open(paths['dlv'], encoding='utf-8', newline='').read()
    assert '\r\n' in text                  # the file's own line endings kept
    assert '%13.4f' % u_lam in text


def test_a_line_with_no_dlv_row_is_reported_and_not_written(paths):
    said = []
    changed = rm.update_dlv(paths['dlv'], [(12345.678, 0.5)], said.append,
                            False)
    assert changed == []
    assert any('has no row' in s for s in said)


def test_the_assignment_names_the_dlv_row(paths):
    rm.update_trans(paths['trans'], paths['dlv'],
                    [(LOW_A, UPP_A, 41473.8197)],
                    rm.read_ids(paths['ids']), lambda *a: None, False)
    trans = IDEN.Trans(paths['trans'])
    tail = IDEN.assignment(trans.records[trans.row(IDEN_A, IDEN_LOW_A)])
    assert IDEN.has_line(tail)
    assert IDEN.obs_row(tail) == 2
    assert IDEN.obs_wavenumber(tail) == pytest.approx(41473.820, abs=0.001)


def test_a_row_that_already_holds_a_line_is_left_alone(paths):
    ids = rm.read_ids(paths['ids'])
    rm.update_trans(paths['trans'], paths['dlv'],
                    [(LOW_A, UPP_A, 41473.8197)], ids, lambda *a: None, False)
    said = []
    made = rm.update_trans(paths['trans'], paths['dlv'],
                           [(LOW_A, UPP_A, 41473.8197)], ids, said.append,
                           False)
    assert made == []
    assert any('left as it stands' in s for s in said)


def test_the_registry_key_is_the_one_the_reader_matches_on(paths, tmp_path):
    import hfs_kappa
    out = worksheet_with(paths, tmp_path, 'accepted, inflated')
    rm.apply_worksheet(paths, out, lambda *a: None, False)
    # four decimals, not the full-precision wn_key, or nothing would match it
    assert hfs_kappa.read_inflated(paths['inflated']) \
        == {'41473.7991': 0.1048}


def test_accepting_a_component_shares_the_weight_by_predicted_intensity(paths,
                                                                        tmp_path):
    out = worksheet_with(paths, tmp_path, 'accepted, inflated')
    rm.apply_worksheet(paths, out, lambda *a: None, False)
    lopt = rm.read_lopt(paths['set_lopt'])
    # 17122.3 and 16267.6 of 33389.9
    assert lopt[(LOW_A, UPP_A)][0] == ''
    assert lopt[(LOW_B, UPP_B)][0] == ''
    records = {k: v for k, v in
               (((rm.lopt_field(r, 'lower_level'),
                  rm.lopt_field(r, 'upper_level')), r)
                for r in rm.records_of(paths['set_lopt']) if r.strip())}
    assert rm.lopt_field(records[(LOW_A, UPP_A)], 'flags') == ''
    assert records[(LOW_A, UPP_A)].endswith('0.5128  cm-1')
    assert records[(LOW_B, UPP_B)].endswith('0.4872  cm-1')


def test_the_inflated_uncertainty_reaches_every_record_of_the_line(paths,
                                                                  tmp_path):
    out = worksheet_with(paths, tmp_path, 'accepted, inflated')
    rm.apply_worksheet(paths, out, lambda *a: None, False)
    lopt = rm.read_lopt(paths['set_lopt'])
    # both components of one observed line share its uncertainty
    assert lopt[(LOW_A, UPP_A)][2] == pytest.approx(0.105)
    assert lopt[(LOW_B, UPP_B)][2] == pytest.approx(0.105)


def test_rejecting_flags_the_record_and_gives_the_line_back(paths, tmp_path):
    # reject the component the set accepts, leaving the other still rejected
    cases, _ = rm.build_cases(paths, 0, False)
    case = dict(cases[0])
    case['level_id1'], case['level_id2'] = LOW_B, UPP_B
    case['decision'], case['Notes'] = 'rejected', 'tested'
    case['unc_wn_inflated'] = ''
    out = str(tmp_path / 'decisions.txt')
    rm.write_worksheet(out, [case])
    rm.apply_worksheet(paths, out, lambda *a: None, False)
    lopt = rm.read_lopt(paths['set_lopt'])
    assert lopt[(LOW_B, UPP_B)][0] == 'P'
    records = {(rm.lopt_field(r, 'lower_level'),
                rm.lopt_field(r, 'upper_level')): r
               for r in rm.records_of(paths['set_lopt']) if r.strip()}
    assert records[(LOW_B, UPP_B)].endswith('0.0000  cm-1')


def test_two_inflated_uncertainties_for_one_line_stop_the_run(paths, tmp_path):
    cases, _ = rm.build_cases(paths, 0, False)
    first = dict(cases[0])
    first['decision'], first['Notes'] = 'accepted, inflated', 'a'
    first['unc_wn_inflated'] = '0.1048'
    second = dict(first)
    second['level_id1'], second['level_id2'] = LOW_B, UPP_B
    second['unc_wn_inflated'] = '0.2000'
    out = str(tmp_path / 'decisions.txt')
    rm.write_worksheet(out, [first, second])
    with pytest.raises(rm.ReviewError) as exc:
        rm.apply_worksheet(paths, out, lambda *a: None, True)
    assert 'two different' in str(exc.value)


def test_the_lopt_input_keeps_its_line_endings(paths, tmp_path):
    out = worksheet_with(paths, tmp_path, 'accepted, inflated')
    rm.apply_worksheet(paths, out, lambda *a: None, False)
    text = io.open(paths['set_lopt'], encoding='utf-8', newline='').read()
    assert '\r\n' not in text          # the miniature set was written with LF


def test_the_dry_run_writes_nothing(paths, tmp_path):
    out = worksheet_with(paths, tmp_path, 'accepted, inflated')
    before = [io.open(p, encoding='utf-8', newline='').read()
              for p in (paths['decisions'], paths['inflated'], paths['dlv'],
                        paths['trans'])]
    rm.apply_worksheet(paths, out, lambda *a: None, True)
    after = [io.open(p, encoding='utf-8', newline='').read()
             for p in (paths['decisions'], paths['inflated'], paths['dlv'],
                       paths['trans'])]
    assert before == after


def test_a_second_apply_does_not_double_the_ledger_row(paths, tmp_path):
    out = worksheet_with(paths, tmp_path, 'accepted, inflated')
    rm.apply_worksheet(paths, out, lambda *a: None, False)
    rm.apply_worksheet(paths, out, lambda *a: None, False)
    ledger = io.open(paths['decisions'], encoding='utf-8', newline='').read()
    assert ledger.count('41473.79906049907,%s' % LOW_A) == 1
    registry = io.open(paths['inflated'], encoding='utf-8', newline='').read()
    assert registry.count('41473.7991') == 1


# ---------------------------------------------------------------------------
# the other direction, and the transition the set never proposed
# ---------------------------------------------------------------------------
LOW_C, UPP_C = '059003.000141', '059003.000454'
IDEN_C, IDEN_LOW_C = 55, 66


def append(path, text):
    with io.open(str(path), 'a', encoding='utf-8', newline='') as fh:
        fh.write(text)


def test_direction_of_names_which_set_disagrees():
    # P means printed but not fitted, and a transition absent from a set's
    # LOPT input is not fitted there either: the same answer to the only
    # question asked, so neither is a disagreement with the other
    assert rm.direction_of('P', '(absent)') is None
    assert rm.direction_of('(absent)', 'P') is None
    assert rm.direction_of('', '') is None
    assert rm.direction_of('P', 'P') is None
    assert rm.direction_of('', 'P') == 'rejected in %s'
    assert rm.direction_of('', '(absent)') == 'rejected in %s'
    assert rm.direction_of('P', '') == 'accepted in %s'
    assert rm.direction_of('(absent)', '') == 'accepted in %s'


def test_an_assignment_only_the_set_fits_is_also_a_disagreement(paths,
                                                               miniature):
    # the baseline flags the second component P as well, so now the set fits
    # an assignment the baseline does not
    text = io.open(paths['base_lopt'], encoding='utf-8', newline='').read()
    lines = text.splitlines()
    lines[1] = lopt_record(41473.799, 0.099, 16761.0, LOW_B, UPP_B, 'P', 0.0)
    with io.open(paths['base_lopt'], 'w', encoding='utf-8',
                 newline='') as fh:
        fh.write('\n'.join(lines) + '\n')
    cases, _ = rm.build_cases(paths, 0, False)
    by_pair = dict(((c['level_id1'], c['level_id2']), c) for c in cases)
    case = by_pair[(LOW_B, UPP_B)]
    assert case['direction'] == 'accepted in iter'
    assert case['lc_flag'] == 'P' and case['iter_flag'] == '(blank)'
    # its own departure is +0.3836 on an uncertainty of 0.0499, and it is the
    # only accepted component, so the line does not sit where it should
    assert case['cg_with'] == '+0.3836'
    assert case['suggested'] == 'review'


def test_a_transition_the_baseline_does_not_carry_at_all_is_reported(paths):
    append(paths['set_lopt'],
           lopt_record(30000.000, 0.030, 900.0, LOW_C, UPP_C, '', 1.0) + '\n')
    append(paths['ids'], '%s\t%d\n%s\t%d\n'
           % (LOW_C, IDEN_LOW_C, UPP_C, IDEN_C))
    cases, _ = rm.build_cases(paths, 0, False)
    case, = [c for c in cases if c['level_id2'] == UPP_C]
    assert case['lc_flag'] == '(absent)'
    assert case['direction'] == 'accepted in iter'


def test_a_pair_neither_set_fits_is_no_disagreement(paths):
    append(paths['base_lopt'],
           lopt_record(30000.000, 0.030, 900.0, LOW_C, UPP_C, 'P', 0.0) + '\n')
    append(paths['ids'], '%s\t%d\n%s\t%d\n'
           % (LOW_C, IDEN_LOW_C, UPP_C, IDEN_C))
    cases, _ = rm.build_cases(paths, 0, False)
    assert [c for c in cases if c['level_id2'] == UPP_C] == []


def test_a_classification_without_a_wn_key_column_names_lines_by_wn_obs(
        tmp_path):
    path = str(tmp_path / 'baseline.csv')
    with io.open(path, 'w', encoding='utf-8', newline='') as fh:
        fh.write('wn_obs,unc_wn_obs,low_id,upp_id,calc_intens,dif_wn_O-C,'
                 'accepted\n')
        fh.write('41473.7991,0.099,%s,%s,17122.3,-0.1600,1\n' % (LOW_A, UPP_A))
    rows, by_line, by_pair = rm.read_classifications(path)
    assert rows[0]['wn_key'] == '41473.7991'
    assert list(by_line) == ['41473.7991']
    assert (LOW_A, UPP_A) in by_pair


def test_level_energies_come_from_the_table_and_then_from_lopt(tmp_path):
    rows = [{'low_id': LOW_A, 'upp_id': UPP_A,
             'low_E': '61357.0180', 'upp_E': '102830.4600'}]
    levels = str(tmp_path / 'LOPT_output_levels.txt')
    with io.open(levels, 'w', encoding='utf-8', newline='') as fh:
        fh.write('Designation\tEnergy\tD1\n')
        fh.write('%s\t99999.0000\t0.01\n' % UPP_A)      # the table wins
        fh.write('%s\t12345.6789\t0.01\n' % LOW_C)      # only here
    energies = rm.read_level_energies(rows, levels)
    assert energies[UPP_A] == pytest.approx(102830.4600)
    assert energies[LOW_C] == pytest.approx(12345.6789)


def baseline_classification(paths, third=True):
    """The baseline's own table: the two components, and optionally a third.

    The baseline names its lines by the observed wavenumber - its `wn` column
    is the column a line is named by - so it writes no `wn_key` column, which
    is what the third component's row has to be found through.
    """
    text = ('wn_obs,unc_wn_obs,obs_intens,char,low_id,upp_id,calc_intens,'
            'dif_wn_O-C,grade,new,accepted,n_accepted\n')
    text += ('41473.79906049907,0.099,16761.3,w,%s,%s,17122.3,-0.1100,5A,1,'
             '1,2\n' % (LOW_A, UPP_A))
    text += ('41473.79906049907,0.099,16761.3,w,%s,%s,16267.6,0.4336,5A,1,'
             '1,2\n' % (LOW_B, UPP_B))
    if third:
        text += ('41473.79906049907,0.099,16761.3,w,%s,%s,4000.0,0.2000,3G,0,'
                 '1,3\n' % (LOW_C, UPP_C))
    with io.open(paths['base_classifications'], 'w', encoding='utf-8',
                 newline='') as fh:
        fh.write(text)


def a_third_component_absent_from_the_set(paths):
    """A transition the baseline fits and the set's table does not carry.

    The corrected wavenumber put the line outside this transition's Ritz
    window, so the set's classification says nothing about it at all - the
    case that has to be assembled from elsewhere.
    """
    append(paths['base_lopt'],
           lopt_record(41473.799, 0.099, 16761.0, LOW_C, UPP_C, '', 0.3)
           + '\n')
    append(paths['ids'], '%s\t%d\n%s\t%d\n'
           % (LOW_C, IDEN_LOW_C, UPP_C, IDEN_C))
    baseline_classification(paths)


def test_a_transition_absent_from_the_sets_table_is_filled_in(paths):
    a_third_component_absent_from_the_set(paths)
    # the set's level energies, which no row of its table carries
    with io.open(paths['set_levels'], 'w', encoding='utf-8', newline='') as fh:
        fh.write('Designation\tEnergy\tD1\n')
        fh.write('%s\t20000.0000\t0.01\n' % LOW_C)
        fh.write('%s\t61473.5000\t0.01\n' % UPP_C)
    cases, _ = rm.build_cases(paths, 0, False)
    case, = [c for c in cases if c['level_id2'] == UPP_C]
    # the line itself from the other component of the blend in the set's table
    assert case['obs_wn'] == '41473.8197'
    assert case['unc_obs_wn'] == '0.0499'
    assert case['wn_key'] == '41473.79906049907'
    assert case['char'] == 'w'
    # the predicted intensity from the baseline's table
    assert case['Icalc'] == '4000'
    # and O-C worked out here: 41473.8197 - (61473.5000 - 20000.0000)
    assert case['dif_wn_O-C'] == '+0.3197'
    assert 'absent from this set' in case['why']
    assert 'other component of the blend' in case['why']
    assert 'against this set\'s level energies' in case['why']


def test_the_filled_in_row_joins_the_blend_arithmetic(paths):
    a_third_component_absent_from_the_set(paths)
    with io.open(paths['set_levels'], 'w', encoding='utf-8', newline='') as fh:
        fh.write('Designation\tEnergy\tD1\n')
        fh.write('%s\t20000.0000\t0.01\n' % LOW_C)
        fh.write('%s\t61473.5000\t0.01\n' % UPP_C)
    cases, _ = rm.build_cases(paths, 0, False)
    case, = [c for c in cases if c['level_id2'] == UPP_C]
    # the accepted component alone stands +0.3836 from Ritz; this one at
    # +0.3197 with a predicted intensity of 4000 against its 16267.6 moves
    # the center of gravity to (0.3836*16267.6 + 0.3197*4000)/20267.6
    assert case['cg_without'] == '+0.3836'
    assert case['cg_with'] == '+0.3710'
    assert case['n_candidates'] == 3
    assert case['n_accepted'] == 1


def test_the_line_list_fills_a_row_no_other_component_can(paths, tmp_path):
    a_third_component_absent_from_the_set(paths)
    # the set's table no longer carries the line at all, so the only place
    # left to find its corrected wavenumber is the set's line list
    with io.open(paths['classifications'], 'w', encoding='utf-8',
                 newline='') as fh:
        fh.write(CLS_HEADER)
    import openpyxl
    book = openpyxl.Workbook()
    sheet = book.active
    sheet.append([rm.COL_OWN, rm.COL_UNC_CORR, rm.COL_INTENS, rm.COL_CHAR,
                  rm.COL_OWN_CORR])
    sheet.append([41473.79906049907, 0.0499, 16761.3, 'w', 41473.8197])
    book.save(paths['corrected_lines'])
    cases, _ = rm.build_cases(paths, 0, False)
    case, = [c for c in cases if c['level_id2'] == UPP_C]
    assert case['obs_wn'] == '41473.8197'
    assert case['unc_obs_wn'] == '0.0499'
    assert case['wn_key'] == '41473.79906049907'
    assert rm.NAME_CORRECTED in case['why']


def test_a_row_that_cannot_be_filled_in_says_so_and_is_kept(paths):
    a_third_component_absent_from_the_set(paths)
    with io.open(paths['classifications'], 'w', encoding='utf-8',
                 newline='') as fh:
        fh.write(CLS_HEADER)
    if os.path.exists(paths['corrected_lines']):
        os.remove(paths['corrected_lines'])
    cases, _ = rm.build_cases(paths, 0, False)
    case, = [c for c in cases if c['level_id2'] == UPP_C]
    assert case['obs_wn'] == ''
    assert case['suggested'] == 'review'
    assert 'could not be found' in case['why']


def test_an_acceptance_inside_the_lines_uncertainty_is_left_alone():
    verdict, why = rm.suggest_accepted({'cg_with': 0.03, 'unc': 0.05,
                                        'code_behind': None, 'n_accepted': 2,
                                        'n_candidates': 2})
    assert verdict == 'stays accepted'
    assert 'within its own uncertainty' in why


def test_an_acceptance_outside_it_is_left_for_the_analyst():
    verdict, why = rm.suggest_accepted({'cg_with': 0.30, 'unc': 0.05,
                                        'code_behind': None, 'n_accepted': 2,
                                        'n_candidates': 2})
    assert verdict == 'review'
    assert 'IDEN2' in why
