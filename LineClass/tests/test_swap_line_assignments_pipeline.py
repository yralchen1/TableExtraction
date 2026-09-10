"""Tests for swap_line_assignments_pipeline.py, swap_line_assignments.py and
the shared path resolution of swap_paths.py.

Every fixture is built from scratch, so nothing here needs the real level
table, the real classification table or the real ledger.  The two levels are
called A and B throughout, A sitting 100 cm^-1 below B before the exchange
and taking B's position after it.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import swap_line_assignments as wrapper
import swap_line_assignments_IDEN as iden
import swap_line_assignments_pipeline as pipe
import swap_paths


A = '059003.000483'
B = '059003.000398'
LOW = '059003.000095'
DATE = '9/6/2026'

E_A = 117592.99          # where A sits before the exchange
E_B = 117686.46          # where B sits before it


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
def write(path, text):
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        fh.write(text)
    return str(path)


def levels_file(tmp_path, e_a=E_A, e_b=E_B, name='LOPT_output_levels.txt'):
    """LOPT's output level table: tab-separated, one header row."""
    rows = ['Designation\tEnergy\tD1\tD2stat\tD2sys\tD2tot\tN_lines\tComments',
            '%s\t%.4f\t0.17\t0.11\t0.00\t0.11\t12\t' % (A, e_a),
            '%s\t%.4f\t0.13\t0.11\t0.00\t0.11\t13\t' % (B, e_b),
            '%s\t1398.3437\t0.0035\t0.0036\t0.0000\t0.0036\t21\t' % LOW]
    return write(tmp_path / name, '\n'.join(rows) + '\n')


CLASS_HEAD = 'wn_obs,low_id,upp_id,obs_intens,new,accepted,low_E,upp_E'


def classifications(tmp_path, rows=None, e_a=E_A, e_b=E_B,
                    name='line_classifications.csv'):
    """The classification table, cut down to the columns this script reads."""
    if rows is None:
        rows = [
            # a legacy identification of B, at B's own position
            ('94153.5202543,%s,%s,120,0,1,1000.0,%.4f' % (LOW, B, e_b)),
            # a new identification of A, which follows the energies on its own
            ('94058.9743210,%s,%s,80,1,1,1000.0,%.4f' % (LOW, A, e_a)),
            # a legacy identification of B that was NOT accepted
            ('90000.0000000,%s,%s,10,0,0,1000.0,%.4f' % (LOW, B, e_b)),
            # somebody else's line entirely
            ('80000.0000000,%s,%s,10,0,1,0.0,80000.0' % (LOW, '059003.000999')),
        ]
    return write(tmp_path / name, CLASS_HEAD + '\n' + '\n'.join(rows) + '\n')


LEDGER_HEAD = 'wn_obs,low_id,upp_id,decision,date,reason'


def ledger(tmp_path, rows=None, name='line_decisions.csv'):
    if rows is None:
        rows = ['94058.9743,%s,%s,accept,9/3/2026,"hand-made identification, '
                'in the LOPT input since iter20"' % (LOW, A)]
    return write(tmp_path / name, LEDGER_HEAD + '\n' + '\n'.join(rows) + '\n')


def overrides(tmp_path, rows=(), name='revised_level_energies.csv'):
    return write(tmp_path / name,
                 'level_id,E_input,comment\n' + ''.join(r + '\n' for r in rows))


# ---------------------------------------------------------------------------
# swap_paths: the directory the command was run from wins
# ---------------------------------------------------------------------------
def test_a_file_in_the_current_directory_is_the_one_meant(tmp_path):
    """Standing in an iteration folder and asking for the transitions file
    can only mean that folder's transitions file."""
    write(tmp_path / 'LOPT.par', 'x\n')
    assert (swap_paths.working_path('LOPT.par', cwd=str(tmp_path))
            == os.path.join(str(tmp_path), 'LOPT.par'))


def test_a_file_absent_there_falls_back_to_the_project(tmp_path):
    """Which is what finds the files that exist in only one place - the
    decision ledger, the level overrides."""
    assert (swap_paths.working_path('line_decisions.csv', cwd=str(tmp_path))
            == os.path.join(swap_paths.PROJECT, 'line_decisions.csv'))


def test_a_directory_resolves_like_a_file(tmp_path):
    os.mkdir(os.path.join(str(tmp_path), 'IDEN2'))
    assert (swap_paths.working_path('IDEN2', cwd=str(tmp_path))
            == os.path.join(str(tmp_path), 'IDEN2'))


def test_a_name_in_neither_place_comes_back_as_the_project_s(tmp_path):
    """So that the caller's own error message names the file it expected."""
    got = swap_paths.working_path('no_such_file.txt', cwd=str(tmp_path))
    assert got == os.path.join(swap_paths.PROJECT, 'no_such_file.txt')


# ---------------------------------------------------------------------------
# The IDEN2 script: a positional may be a row number
# ---------------------------------------------------------------------------
def test_a_bare_number_is_a_row_number_and_an_identifier_is_not():
    assert iden.row_number('721') == 721
    assert iden.row_number(' 721 ') == 721
    assert iden.row_number(A) is None
    assert iden.row_number('') is None


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------
def test_the_level_table_is_read_by_identifier(tmp_path):
    got = pipe.read_lopt_levels(levels_file(tmp_path))
    assert got[A] == pytest.approx(E_A)
    assert got[B] == pytest.approx(E_B)
    assert 'Designation' not in got          # the header row is not a level


def test_a_wavenumber_is_written_to_four_decimals():
    """The classification table carries the full binary expansion; the ledger
    matches within 0.01 cm^-1, so four decimals identify a line many times
    over and are what the file is written in."""
    assert pipe.wavenumber('94153.52025432966') == '94153.5203'
    assert pipe.wavenumber('') == ''


def test_the_energies_the_pipeline_has_come_from_both_ends(tmp_path):
    got = pipe.pipeline_energies(classifications(tmp_path), (A, B))
    assert got[A] == pytest.approx(E_A)
    assert got[B] == pytest.approx(E_B)


def test_only_accepted_identifications_of_the_published_list_are_collected(tmp_path):
    """new = 0 and accepted = 1: an identification that comes back on every
    run whatever the energies are, which is exactly the dangerous kind."""
    got = pipe.legacy_assignments(classifications(tmp_path), (A, B))
    assert [(g['wn'], g['level']) for g in got] == [('94153.5203', B)]


# ---------------------------------------------------------------------------
# The energies
# ---------------------------------------------------------------------------
def test_a_table_already_exchanged_is_taken_as_it_stands(tmp_path):
    """A LOPT run made after the exchange: its two energies are the answer."""
    lopt = pipe.read_lopt_levels(levels_file(tmp_path, e_a=E_B, e_b=E_A))
    current = {A: E_A, B: E_B}
    targets, whence = pipe.target_energies((A, B), lopt, current,
                                           False, False, 0.5)
    assert targets[A] == pytest.approx(E_B)
    assert targets[B] == pytest.approx(E_A)
    assert 'already been exchanged' in whence


def test_a_table_not_yet_exchanged_is_crossed_here(tmp_path):
    """No LOPT run yet: the same exchange, done by arithmetic instead."""
    lopt = pipe.read_lopt_levels(levels_file(tmp_path))
    targets, whence = pipe.target_energies((A, B), lopt, {A: E_A, B: E_B},
                                           False, False, 0.5)
    assert targets[A] == pytest.approx(E_B)
    assert targets[B] == pytest.approx(E_A)
    assert 'exchanged here' in whence


def test_exchange_does_the_crossing_itself(tmp_path):
    lopt = pipe.read_lopt_levels(levels_file(tmp_path))
    targets, _ = pipe.target_energies((A, B), lopt, {A: E_A, B: E_B},
                                      True, False, 0.5)
    assert targets[A] == pytest.approx(E_B)
    assert targets[B] == pytest.approx(E_A)


def test_exchange_on_a_table_already_exchanged_is_refused(tmp_path):
    """It would put the two levels back."""
    lopt = pipe.read_lopt_levels(levels_file(tmp_path, e_a=E_B, e_b=E_A))
    with pytest.raises(ValueError, match='would put them back'):
        pipe.target_energies((A, B), lopt, {A: E_A, B: E_B}, True, False, 0.5)


def test_as_fitted_on_a_table_not_yet_exchanged_is_refused(tmp_path):
    """It would write both levels back where they already are."""
    lopt = pipe.read_lopt_levels(levels_file(tmp_path))
    with pytest.raises(ValueError, match='record no exchange'):
        pipe.target_energies((A, B), lopt, {A: E_A, B: E_B}, False, True, 0.5)


def test_a_table_agreeing_with_neither_is_refused(tmp_path):
    """Something other than the exchange has moved these levels."""
    lopt = pipe.read_lopt_levels(levels_file(tmp_path, e_a=1.0, e_b=2.0))
    with pytest.raises(ValueError, match='neither what'):
        pipe.target_energies((A, B), lopt, {A: E_A, B: E_B}, False, False, 0.5)


def test_without_the_pipeline_energies_one_of_the_two_must_be_named(tmp_path):
    """Nothing to compare the table with, so it cannot be read either way."""
    lopt = pipe.read_lopt_levels(levels_file(tmp_path))
    with pytest.raises(ValueError, match='no energy for'):
        pipe.target_energies((A, B), lopt, {}, False, False, 0.5)
    targets, _ = pipe.target_energies((A, B), lopt, {}, True, False, 0.5)
    assert targets[A] == pytest.approx(E_B)
    targets, _ = pipe.target_energies((A, B), lopt, {}, False, True, 0.5)
    assert targets[A] == pytest.approx(E_A)


def test_two_levels_at_the_same_energy_have_nothing_to_exchange(tmp_path):
    lopt = pipe.read_lopt_levels(levels_file(tmp_path, e_a=E_A, e_b=E_A))
    with pytest.raises(ValueError, match='nothing to exchange'):
        pipe.target_energies((A, B), lopt, {}, False, False, 0.5)


# ---------------------------------------------------------------------------
# The overlay
# ---------------------------------------------------------------------------
def test_a_level_not_yet_overridden_gets_a_row():
    rows, report = pipe.override_rows([], {A: E_B, B: E_A}, (A, B), DATE)
    assert [r['level_id'] for r in rows] == [A, B]
    assert rows[0]['E_input'] == '%.4f' % E_B
    assert pipe.mark(A, B) in rows[0]['comment']
    assert len(report) == 2


def test_a_level_already_overridden_keeps_its_row_and_its_comment():
    """A second row for the same level would be read as a contradiction."""
    old = [{'level_id': A, 'E_input': '117000.0000',
            'comment': 're-positioned from somewhere'}]
    rows, _ = pipe.override_rows(old, {A: E_B, B: E_A}, (A, B), DATE)
    assert len(rows) == 2
    assert rows[0]['E_input'] == '%.4f' % E_B
    assert rows[0]['comment'].startswith('re-positioned from somewhere')
    assert pipe.mark(A, B) in rows[0]['comment']


def test_recording_the_same_exchange_twice_refreshes_only_the_energy():
    """The comment says the exchange has been made; saying it again would
    read as two exchanges."""
    rows, _ = pipe.override_rows([], {A: E_B, B: E_A}, (A, B), DATE)
    once = rows[0]['comment']
    rows, report = pipe.override_rows(rows, {A: E_B + 0.01, B: E_A},
                                      (A, B), '9/9/2026')
    assert rows[0]['comment'] == once
    assert rows[0]['E_input'] == '%.4f' % (E_B + 0.01)
    assert 'already recorded' in report[0]


# ---------------------------------------------------------------------------
# The ledger
# ---------------------------------------------------------------------------
def test_an_order_naming_one_level_is_pointed_at_the_other():
    rows = [{'wn_obs': '94058.9743', 'low_id': LOW, 'upp_id': A,
             'decision': 'accept', 'date': '9/3/2026', 'reason': 'hand-made'}]
    touched, skipped = pipe.rekey_ledger(rows, A, B, DATE)
    assert rows[0]['upp_id'] == B and rows[0]['low_id'] == LOW
    assert rows[0]['reason'] == (
        'hand-made; Lines belonging to levels %s and %s were swapped on %s.'
        % (B, A, DATE))
    assert len(touched) == 1 and skipped == []


def test_an_order_naming_neither_level_is_untouched():
    rows = [{'wn_obs': '1.0', 'low_id': LOW, 'upp_id': '059003.000999',
             'decision': 'reject', 'date': '', 'reason': 'no'}]
    before = dict(rows[0])
    pipe.rekey_ledger(rows, A, B, DATE)
    assert rows[0] == before


def test_an_order_already_marked_is_not_re_keyed_a_second_time():
    """Running the script twice must not point the order back at the level it
    started from."""
    rows = [{'wn_obs': '94058.9743', 'low_id': LOW, 'upp_id': A,
             'decision': 'accept', 'date': '', 'reason': 'hand-made'}]
    pipe.rekey_ledger(rows, A, B, DATE)
    once = dict(rows[0])
    touched, skipped = pipe.rekey_ledger(rows, A, B, '9/9/2026')
    assert rows[0] == once
    assert touched == [] and len(skipped) == 1


def test_the_mark_does_not_depend_on_the_order_of_the_two_levels():
    assert pipe.mark(A, B) == pipe.mark(B, A)


def test_a_legacy_identification_gets_a_reject_and_an_accept():
    legacy = [{'wn': '94153.5203', 'low': LOW, 'upp': B, 'level': B,
               'intens': '120'}]
    out = pipe.legacy_orders(legacy, A, B, DATE, True)
    assert [(o['decision'], o['low_id'], o['upp_id']) for o in out] == [
        ('reject', LOW, B), ('accept', LOW, A)]
    assert all(pipe.mark(A, B) in o['reason'] for o in out)


def test_rejects_only_leaves_the_accept_to_the_classification():
    legacy = [{'wn': '94153.5203', 'low': LOW, 'upp': B, 'level': B,
               'intens': '120'}]
    out = pipe.legacy_orders(legacy, A, B, DATE, False)
    assert [o['decision'] for o in out] == ['reject']


def test_an_order_that_says_the_same_thing_is_not_written_twice():
    rows = [{'wn_obs': '94153.5203', 'low_id': LOW, 'upp_id': B,
             'decision': 'reject', 'date': '', 'reason': 'already'}]
    new = [{'wn_obs': '94153.5203', 'low_id': LOW, 'upp_id': B,
            'decision': 'reject', 'date': DATE, 'reason': 'again'}]
    added, already = pipe.merge_orders(rows, new)
    assert added == [] and len(already) == 1
    assert len(rows) == 1


def test_an_order_that_contradicts_one_already_there_stops_the_script():
    """classify_lines.py aborts on two rows ruling differently on the same
    assignment, so the clash has to be caught here, not on the next run."""
    rows = [{'wn_obs': '94153.5203', 'low_id': LOW, 'upp_id': B,
             'decision': 'accept', 'date': '', 'reason': 'kept'}]
    new = [{'wn_obs': '94153.5203', 'low_id': LOW, 'upp_id': B,
            'decision': 'reject', 'date': DATE, 'reason': 'moved'}]
    with pytest.raises(ValueError, match='already orders accept'):
        pipe.merge_orders(rows, new)


# ---------------------------------------------------------------------------
# End to end
# ---------------------------------------------------------------------------
def run(tmp_path, extra=(), **kw):
    """Run the script on the four fixture files, making any that is missing.

    A file already in tmp_path is left as it is, so that a second run sees
    what the first one wrote - which is what the idempotency test needs.
    """
    def have(name, make):
        path = str(tmp_path / name)
        return path if os.path.exists(path) else make(tmp_path)

    argv = [A, B,
            '--levels', kw.get('levels')
            or have('LOPT_output_levels.txt', levels_file),
            '--overrides', kw.get('overrides')
            or have('revised_level_energies.csv', overrides),
            '--ledger', kw.get('ledger') or have('line_decisions.csv', ledger),
            '--classifications', kw.get('classifications')
            or have('line_classifications.csv', classifications),
            '--date', DATE, '--no-backup'] + list(extra)
    return pipe.main(argv)


def test_end_to_end_writes_both_files(tmp_path, capsys):
    assert run(tmp_path, ['--exchange']) == 0
    fields, rows = pipe.read_table(str(tmp_path / 'revised_level_energies.csv'))
    assert {r['level_id']: r['E_input'] for r in rows} == {
        A: '%.4f' % E_B, B: '%.4f' % E_A}

    _fields, orders = pipe.read_table(str(tmp_path / 'line_decisions.csv'))
    keyed = {(r['wn_obs'], r['low_id'], r['upp_id']): r['decision']
             for r in orders}
    # the order that named A now names B
    assert keyed[('94058.9743', LOW, B)] == 'accept'
    # the legacy identification is left to classify_lines.py, which reads the
    # exchange from the comment written into the overrides above and moves it
    # to the other identifier itself
    assert not [k for k in keyed if k[0] == '94153.5203']
    # the unaccepted legacy row and the unrelated line are not ruled on
    assert not [k for k in keyed if k[0] == '90000.0000']
    assert not [k for k in keyed if '059003.000999' in k]


def test_legacy_orders_writes_the_pair_of_orders(tmp_path):
    """--legacy-orders is the old behaviour, kept for a run that records the
    exchange nowhere else."""
    assert run(tmp_path, ['--exchange', '--legacy-orders']) == 0
    _fields, orders = pipe.read_table(str(tmp_path / 'line_decisions.csv'))
    keyed = {(r['wn_obs'], r['low_id'], r['upp_id']): r['decision']
             for r in orders}
    # the legacy identification of B is rejected there and accepted at A
    assert keyed[('94153.5203', LOW, B)] == 'reject'
    assert keyed[('94153.5203', LOW, A)] == 'accept'


def test_no_overrides_still_writes_the_legacy_orders(tmp_path):
    """Nothing then records the exchange, so the ledger has to."""
    assert run(tmp_path, ['--exchange', '--no-overrides']) == 0
    _fields, orders = pipe.read_table(str(tmp_path / 'line_decisions.csv'))
    keyed = {(r['wn_obs'], r['low_id'], r['upp_id']): r['decision']
             for r in orders}
    assert keyed[('94153.5203', LOW, B)] == 'reject'
    assert keyed[('94153.5203', LOW, A)] == 'accept'


def test_end_to_end_is_idempotent(tmp_path):
    """The second run must not undo the first: no order may be re-keyed back
    and no order duplicated."""
    run(tmp_path, ['--exchange'])
    first = open(str(tmp_path / 'line_decisions.csv'), encoding='utf-8').read()
    # the level table has to be the exchanged one now, as it would be after
    # a LOPT run; the classification table still has the old energies.
    levels_file(tmp_path, e_a=E_B, e_b=E_A)
    assert run(tmp_path) == 0
    assert open(str(tmp_path / 'line_decisions.csv'),
                encoding='utf-8').read() == first


def test_the_files_are_written_with_lf_endings(tmp_path):
    """Every text file tracked here is LF; a csv module default of \\r\\n
    would put this file out of step with the repository."""
    run(tmp_path, ['--exchange'])
    for name in ('revised_level_energies.csv', 'line_decisions.csv'):
        raw = open(str(tmp_path / name), 'rb').read()
        assert b'\r\n' not in raw


def test_a_dry_run_writes_nothing(tmp_path):
    led = ledger(tmp_path)
    before = open(led, encoding='utf-8').read()
    assert run(tmp_path, ['--exchange', '--dry-run'], ledger=led) == 0
    assert open(led, encoding='utf-8').read() == before
    assert not os.path.exists(str(tmp_path / 'revised_level_energies.csv.bak'))


def test_a_row_number_is_refused_here(tmp_path):
    """IDEN2's numbering means nothing to the pipeline."""
    with pytest.raises(ValueError, match='row number of enlev.dat'):
        pipe.main(['721', '723', '--levels', levels_file(tmp_path)])


def test_the_same_identifier_twice_is_refused(tmp_path):
    assert pipe.main([A, A, '--levels', levels_file(tmp_path)]) == 2


# ---------------------------------------------------------------------------
# The wrapper
# ---------------------------------------------------------------------------
class WArgs(object):
    def __init__(self, **kw):
        self.id1, self.id2 = A, B
        self.only = 'lopt,iden,pipeline'
        self.after_lopt = self.dry_run = self.no_backup = False
        self.rejects_only = self.force = self.no_lopt_run = False
        self.legacy_orders = self.no_reweight = False
        self.reweight_tol = None
        self.index = self.date = None
        self.lopt_command = 'lopt'
        self.__dict__.update(kw)


def test_the_steps_run_in_the_order_they_have_to():
    """LOPT first, so that a pair it refuses never reaches the other two."""
    assert wrapper.wanted(WArgs()) == ['lopt', 'iden', 'pipeline']
    assert wrapper.wanted(WArgs(only='pipeline,lopt')) == ['lopt', 'pipeline']


def test_after_lopt_runs_only_the_pipeline_step():
    assert wrapper.wanted(WArgs(after_lopt=True)) == ['pipeline']


def test_an_unknown_step_is_refused():
    with pytest.raises(ValueError, match='no step called'):
        wrapper.wanted(WArgs(only='icalc'))


def test_the_pipeline_step_is_left_to_read_the_level_table_itself():
    """Whether LOPT ran or not, the step works out which table it has."""
    assert wrapper.step_argv('pipeline', WArgs()) == [A, B]
    assert wrapper.step_argv('pipeline', WArgs(after_lopt=True)) == [A, B]


def test_lopt_runs_in_the_directory_its_parameter_file_was_found_in(tmp_path,
                                                                   monkeypatch):
    """LOPT.par names the transitions file, so the directory decides the fit."""
    (tmp_path / 'LOPT.par').write_text('x\n')
    monkeypatch.chdir(str(tmp_path))
    assert wrapper.lopt_directory() == str(tmp_path)


def test_a_missing_lopt_command_is_a_warning_not_a_failure(tmp_path):
    """The run has to go on: the pipeline step can do the exchange alone."""
    trouble = wrapper.run_lopt('there_is_no_such_command_here',
                               str(tmp_path))
    assert trouble and 'status' in trouble


def test_lopt_that_works_reports_nothing(tmp_path):
    assert wrapper.run_lopt('python -c "pass"', str(tmp_path)) is None


def test_each_step_is_given_only_the_options_it_has():
    lopt_argv = wrapper.step_argv('lopt', WArgs(index=[721, 723], force=True,
                                                rejects_only=True))
    assert lopt_argv == [A, B]
    iden_argv = wrapper.step_argv('iden', WArgs(index=[721, 723], force=True))
    assert iden_argv == [A, B, '--index', '721', '723', '--force']


def test_dry_run_reaches_every_step():
    for step in wrapper.STEPS:
        assert '--dry-run' in wrapper.step_argv(step, WArgs(dry_run=True))


def test_the_reweighting_options_reach_the_lopt_step_only():
    assert '--no-reweight' in wrapper.step_argv('lopt', WArgs(no_reweight=True))
    assert wrapper.step_argv('lopt', WArgs(reweight_tol=0.01))[-2:] \
        == ['--reweight-tol', '0.01']
    assert '--no-reweight' not in wrapper.step_argv(
        'pipeline', WArgs(no_reweight=True))
    assert '--no-reweight' not in wrapper.step_argv('lopt', WArgs())


def test_legacy_orders_reaches_the_pipeline_step_only():
    assert '--legacy-orders' in wrapper.step_argv(
        'pipeline', WArgs(legacy_orders=True))
    assert '--legacy-orders' not in wrapper.step_argv(
        'lopt', WArgs(legacy_orders=True))
    assert '--legacy-orders' not in wrapper.step_argv('pipeline', WArgs())
