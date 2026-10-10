"""Tests for `make_LOPT_input`.

Run from the LineClass directory:  python -m pytest tests -q

One observed line is one measurement, so every record LOPT reads for it has to
carry the same wavenumber and the same uncertainty; two uncertainties on one
wavenumber are two weights on one measurement.  The hyperfine term of
`total_unc` is a property of the transition, not of the line, so it is what
pulls the components of a blend apart, and `blend_uncertainty` is what puts
them back together.  The tests here fix:

  * the weighted mean itself, including the degenerate cases - nothing
    fitted, and one component carrying the whole line;
  * that the written file gives one uncertainty per observed wavenumber,
    counting the records flagged `P` that share it;
  * that the shared value lies between the components' own, and equals the
    single component's value where a line has only one;
  * that the weights are untouched by all of this.
"""
import io
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (ROOT, os.path.dirname(ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

import make_LOPT_input as M              # noqa: E402


def row(wn, unc, low, upp, accepted=1, calc=1.0, intens=100.0):
    return {'wn_obs': wn, 'unc_wn_obs': unc, 'obs_intens': intens,
            'low_id': low, 'upp_id': upp, 'accepted': accepted,
            'calc_intens': calc}


def read_back(path):
    """{wavenumber: [(uncertainty, flags, weight)]} of a written file."""
    out = {}
    for rec in io.open(path, encoding='ascii', newline=''):
        rec = rec.rstrip('\r\n')
        if not rec.strip():
            continue
        got = {k: rec[a - 1:b].strip() for k, (a, b) in M.FIELDS.items()}
        out.setdefault(got['wavenumber'], []).append(
            (got['uncertainty'], got['flags'], got['weight']))
    return out


# ---------------------------------------------------------------------------
# The mean
# ---------------------------------------------------------------------------
def test_weighted_mean():
    assert M.blend_uncertainty([0.10, 0.20], [0.75, 0.25]) == pytest.approx(0.125)


def test_whole_line_on_one_component():
    """A component carrying the line governs its width."""
    assert M.blend_uncertainty([0.10, 0.90], [1.0, 0.0]) == pytest.approx(0.10)


def test_nothing_fitted_takes_the_plain_mean():
    """Records that are all flagged have no weights to average with."""
    assert M.blend_uncertainty([0.10, 0.30], [0.0, 0.0]) == pytest.approx(0.20)


# ---------------------------------------------------------------------------
# The file
# ---------------------------------------------------------------------------
def test_one_uncertainty_per_observed_line(tmp_path):
    w_hfs = {'059003.000100': 0.30, '059003.000400': 0.0}
    rows = [
        # a blend whose two components sit on levels of different widths
        row('20000.000', '0.050', '059003.000100', '059003.000300', calc=3.0),
        row('20000.000', '0.050', '059003.000200', '059003.000400', calc=1.0),
        # a line with one accepted component and one rejected candidate.
        # The pairs of levels differ from the ones above: a transition gets
        # one record, so reusing a pair here would drop a row.
        row('19000.000', '0.050', '059003.000100', '059003.000301'),
        row('19000.000', '0.050', '059003.000200', '059003.000401',
            accepted=0),
        # a line assigned once
        row('18000.000', '0.050', '059003.000200', '059003.000402'),
    ]
    path = str(tmp_path / 'lines.txt')
    M.write_lines_file(rows, path, w_hfs)
    got = read_back(path)

    for wn, recs in got.items():
        assert len({u for u, _, _ in recs}) == 1, wn

    wide = (0.050 ** 2 + 0.30 ** 2) ** 0.5     # the component on level ...100
    narrow = 0.050

    # 3:1 in calculated intensity, so three quarters of the wide component
    assert float(got['20000.000'][0][0]) == pytest.approx(
        0.75 * wide + 0.25 * narrow, abs=5e-4)
    # the flagged record does not dilute the one that is fitted
    assert float(got['19000.000'][0][0]) == pytest.approx(wide, abs=5e-4)
    # nothing to share: the line keeps its own value
    assert float(got['18000.000'][0][0]) == pytest.approx(narrow, abs=5e-4)


def test_the_shared_value_lies_between_the_components(tmp_path):
    w_hfs = {'059003.000100': 0.30}
    rows = [
        row('20000.000', '0.050', '059003.000100', '059003.000300', calc=1.0),
        row('20000.000', '0.050', '059003.000200', '059003.000400', calc=1.0),
    ]
    path = str(tmp_path / 'lines.txt')
    M.write_lines_file(rows, path, w_hfs)
    shared = float(read_back(path)['20000.000'][0][0])
    assert 0.050 < shared < (0.050 ** 2 + 0.30 ** 2) ** 0.5


def test_the_weights_are_untouched(tmp_path):
    w_hfs = {'059003.000100': 0.30}
    rows = [
        row('20000.000', '0.050', '059003.000100', '059003.000300', calc=3.0),
        row('20000.000', '0.050', '059003.000200', '059003.000400', calc=1.0),
        row('20000.000', '0.050', '059003.000200', '059003.000500',
            accepted=0, calc=5.0),
    ]
    path = str(tmp_path / 'lines.txt')
    M.write_lines_file(rows, path, w_hfs)
    recs = read_back(path)['20000.000']
    assert sorted(w for _, _, w in recs) == ['0.0000', '0.2500', '0.7500']
    assert sorted(f for _, f, _ in recs) == ['', '', 'P']


# ---------------------------------------------------------------------------
# The guard
# ---------------------------------------------------------------------------
def test_check_sync_sees_a_line_with_two_uncertainties(tmp_path):
    """`check_sync.one_uncertainty_per_line` is what keeps this from coming
    back silently, whichever program wrote the file."""
    import check_sync as CS

    rows = [
        M.format_line(20000.0, 0.050, 10.0, '059003.000100',
                      '059003.000300', '', 0.75),
        M.format_line(20000.0, 0.300, 10.0, '059003.000200',
                      '059003.000400', '', 0.25),
        M.format_line(19000.0, 0.050, 10.0, '059003.000100',
                      '059003.000300', '', 1.0),
    ]
    path = str(tmp_path / 'lines.txt')
    with io.open(path, 'w', encoding='ascii', newline='') as fh:
        for r in rows:
            fh.write(r + M.EOL_LINES)

    split = CS.one_uncertainty_per_line(path)
    assert list(split) == ['20000.000']
    assert {u for u, _, _, _ in split['20000.000']} == {'0.050', '0.300'}


# ---------------------------------------------------------------------------
# One record per transition
# ---------------------------------------------------------------------------
def test_the_accepted_record_of_a_transition_is_the_one_kept():
    """Two observed lines assigned to one pair of levels: one energy
    difference, so one record, and the accepted one is it."""
    rejected = row('31391.380', '0.050', '059003.000141', '059003.000218',
                   accepted=0)
    taken = row('31391.653', '0.050', '059003.000141', '059003.000218')
    for rows in ([rejected, taken], [taken, rejected]):
        kept, dropped = M.one_record_per_transition(rows)
        assert [r['wn_obs'] for r in kept] == ['31391.653']
        assert [r['wn_obs'] for r in dropped] == ['31391.380']


def test_two_flagged_records_keep_the_first():
    """Neither is in the fit, so which wavenumber is printed does not
    matter; what matters is that only one is."""
    rows = [row('31391.653', '0.050', '059003.000141', '059003.000218',
                accepted=0),
            row('31391.380', '0.050', '059003.000141', '059003.000218',
                accepted=0)]
    kept, dropped = M.one_record_per_transition(rows)
    assert [r['wn_obs'] for r in kept] == ['31391.653']
    assert len(dropped) == 1


def test_a_transition_accepted_twice_stops_the_run():
    rows = [row('31391.653', '0.050', '059003.000141', '059003.000218'),
            row('31391.380', '0.050', '059003.000141', '059003.000218')]
    with pytest.raises(SystemExit) as err:
        M.one_record_per_transition(rows)
    assert '059003.000141' in str(err.value)


def test_the_written_file_holds_each_transition_once(tmp_path):
    rows = [
        row('31391.653', '0.050', '059003.000141', '059003.000218'),
        row('31391.380', '0.050', '059003.000141', '059003.000218',
            accepted=0),
        row('20000.000', '0.050', '059003.000200', '059003.000400'),
    ]
    path = str(tmp_path / 'lines.txt')
    written, accepted, flagged, widened, dropped = M.write_lines_file(
        rows, path, {})
    assert (written, accepted, flagged, dropped) == (2, 2, 0, 1)
    got = read_back(path)
    assert sorted(got) == ['20000.000', '31391.653']


def test_a_dropped_record_does_not_weigh_on_the_line_it_left(tmp_path):
    """The repeat goes before the weights are worked out, so the blend it
    was a component of is shared out between the records that remain."""
    rows = [
        # the transition this record names is already accepted at 31391.653
        row('20000.000', '0.050', '059003.000141', '059003.000218',
            accepted=0),
        row('31391.653', '0.050', '059003.000141', '059003.000218'),
        row('20000.000', '0.050', '059003.000200', '059003.000400',
            calc=3.0),
        row('20000.000', '0.050', '059003.000300', '059003.000500',
            calc=1.0),
    ]
    path = str(tmp_path / 'lines.txt')
    M.write_lines_file(rows, path, {})
    recs = read_back(path)['20000.000']
    assert sorted(w for _, _, w in recs) == ['0.2500', '0.7500']


# ---------------------------------------------------------------------------
# Where the inputs are looked for
# ---------------------------------------------------------------------------
def write_classifications(path, rows):
    cols = ['wn_obs', 'unc_wn_obs', 'obs_intens', 'low_id', 'upp_id',
            'accepted', 'calc_intens']
    with io.open(path, 'w', encoding='utf-8', newline='') as fh:
        fh.write(','.join(cols) + '\n')
        for rec in rows:
            fh.write(','.join(str(rec[c]) for c in cols) + '\n')


def test_a_set_without_the_samples_uses_the_project_copies(tmp_path, capsys):
    """Run for an iteration folder, the bare-named inputs come from the
    project directory instead of stopping the run."""
    set_dir = tmp_path / 'iter'
    set_dir.mkdir()
    write_classifications(
        str(set_dir / 'line_classifications.csv'),
        [row('20000.000', '0.050', '059003.000100', '059003.000300')])

    assert M.main(['--classifications',
                   str(set_dir / 'line_classifications.csv')]) == 0
    out = capsys.readouterr().out

    # the three files were written beside the classification table
    for name in ('LOPT_input_lines.txt', 'LOPT_fixlev.txt', 'LOPT.par'):
        assert (set_dir / name).exists(), name
    # and the inputs the set has not got came from the project directory
    assert M.SCRIPT_DIR in out
    assert 'no hyperfine widths applied' not in out


def test_the_set_own_copy_wins(tmp_path):
    """A file the set does have is the one used."""
    set_dir = tmp_path / 'iter'
    set_dir.mkdir()
    write_classifications(
        str(set_dir / 'line_classifications.csv'),
        [row('20000.000', '0.050', '059003.000100', '059003.000300')])
    (set_dir / 'Pr3_line_class5_fixlev.txt').write_text(
        'this set only\n', encoding='utf-8')

    assert M.main(['--classifications',
                   str(set_dir / 'line_classifications.csv')]) == 0
    assert (set_dir / 'LOPT_fixlev.txt').read_text(
        encoding='utf-8') == 'this set only\n'


# ---------------------------------------------------------------------------
# The parameter file says the line uncertainties are statistical
# ---------------------------------------------------------------------------
def par_sample(tmp_path, extra=()):
    text = ['a.txt  ; TRANSITIONS input file name',
            'b.txt  ; FIXED LEVELS input file name',
            'c.txt  ; Levels OUTPUT file name',
            'd.txt  ; Transitions OUTPUT file name',
            'Y      ; TUNE single-line levels (Y/N)?'] + list(extra)
    path = tmp_path / 'sample.par'
    path.write_text('\n'.join(text) + '\n', encoding='utf-8')
    return str(path)


def kind_lines(path):
    with open(path, encoding='utf-8') as fh:
        return [t.rstrip('\n') for t in fh
                if M.PAR_UNC_KIND_KEY in t.upper()]


def test_the_par_file_gets_stat_where_the_sample_has_no_kind(tmp_path):
    out = str(tmp_path / 'LOPT.par')
    M.write_par_file(par_sample(tmp_path), out, 'l', 'f', 'lo', 'li')
    kinds = kind_lines(out)
    assert len(kinds) == 1 and kinds[0].split(';')[0].strip() == 'stat'


def test_the_par_file_turns_a_sample_tot_into_stat(tmp_path):
    sample = par_sample(tmp_path, [
        'tot     ; Kind of measurement uncertainty given in the lines input '
        'file [tot/stat] (default tot)'])
    out = str(tmp_path / 'LOPT.par')
    M.write_par_file(sample, out, 'l', 'f', 'lo', 'li')
    kinds = kind_lines(out)
    assert len(kinds) == 1 and kinds[0].split(';')[0].strip() == 'stat'
    assert kinds[0].endswith('(default tot)')
