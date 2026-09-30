"""Tests for `--set` in level_positions.py and find_unknown_levels.py.

Run from the LineClass directory:  python -m pytest tests -q

A working set (`iter/`, `final/`) holds its own configuration, classification
table, LOPT output and IDEN2; a file it does not hold is the project
directory's.  `--set` must point every one of those readers at the set, and
put what the run writes into the set, without touching any file.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

import classify_lines as cl                # noqa: E402
import config                              # noqa: E402
import find_unknown_levels as fu           # noqa: E402
import level_positions as lp               # noqa: E402

ITER = os.path.join(HERE, 'iter')
real_set = pytest.mark.skipif(
    not os.path.exists(os.path.join(ITER, 'lineclass_config.toml')),
    reason='needs the iter working set')


@pytest.fixture
def keep_globals(monkeypatch):
    """Restore what use_set rebinds: the module paths and cl's configuration."""
    for name in ('HFS_FILE', 'LEVEL_IDS', 'ENLEV', 'LOPT_LEVELS_FILE'):
        monkeypatch.setattr(lp, name, getattr(lp, name))
    for name in ('ENLEV', 'ID_MAP'):
        monkeypatch.setattr(fu, name, getattr(fu, name))
    yield
    cl.apply_config(config.load())


@real_set
def test_level_positions_reads_the_set(keep_globals):
    args = lp.parse_args(['--set', ITER, '--scan'])
    lp.use_set(args, log=lambda *a: None)
    assert lp.ENLEV == os.path.join(ITER, 'IDEN2', 'enlev.dat')
    assert lp.LEVEL_IDS == os.path.join(ITER, 'IDEN2', 'IDEN_level_ids.txt')
    assert lp.LOPT_LEVELS_FILE == os.path.join(ITER, 'LOPT_output_levels.txt')
    assert os.path.dirname(cl.OUTPUT_CSV) == ITER
    assert os.path.dirname(cl.LINES_FILE) == ITER
    # what the run writes describes the set
    assert args.out == os.path.join(ITER, 'level_positions.csv')
    assert args.firm == os.path.join(ITER, lp.FIRM_FILE)


@real_set
def test_level_positions_widths_fall_back_but_are_written_in_the_set(
        keep_globals):
    lp.use_set(lp.parse_args(['--set', ITER]), log=lambda *a: None)
    read = lp.HFS_FILE
    lp.use_set(lp.parse_args(['--set', ITER, '--fit-hfs']),
               log=lambda *a: None)
    name = 'level_hfs_widths.csv'
    if not os.path.exists(os.path.join(ITER, name)):
        assert read == os.path.join(HERE, name)
    assert lp.HFS_FILE == os.path.join(ITER, name)


def test_level_positions_without_set_is_the_baseline():
    args = lp.parse_args(['--scan'])
    assert args.set_dir is None
    assert args.out == 'level_positions.csv'
    assert lp.ENLEV == os.path.join(HERE, 'IDEN2', 'enlev.dat')


def test_level_positions_rejects_a_missing_set(tmp_path):
    with pytest.raises(SystemExit):
        lp.use_set(lp.parse_args(['--set', str(tmp_path / 'nowhere')]),
                   log=lambda *a: None)


def test_find_unknown_levels_baseline_defaults(keep_globals):
    args = fu.parse_args([])
    assert fu.use_set(args) == []
    assert (args.src, args.out, args.log) == (fu.IN_CSV, fu.OUT_CSV,
                                              fu.LOGFILE)
    assert fu.ENLEV == os.path.join(HERE, 'IDEN2', 'enlev.dat')


def test_find_unknown_levels_in_a_set(tmp_path, keep_globals, capsys):
    (tmp_path / 'IDEN2').mkdir()
    args = fu.parse_args(['--set', str(tmp_path)])
    extra = fu.use_set(args)
    assert extra == ['--set', str(tmp_path)]
    assert fu.ENLEV == str(tmp_path / 'IDEN2' / 'enlev.dat')
    assert fu.ID_MAP == str(tmp_path / 'IDEN2' / 'IDEN_level_ids.txt')
    # the ranked list is the project's unless the set has its own copy
    assert args.src == os.path.join(HERE, 'unfound_levels.csv')
    assert args.out == str(tmp_path / 'found_levels.csv')
    assert args.log == str(tmp_path / 'find_unknown_levels.log')


def test_find_unknown_levels_prefers_the_sets_own_list(tmp_path,
                                                      keep_globals, capsys):
    (tmp_path / 'unfound_levels.csv').write_text('idx\n', newline='\n')
    args = fu.parse_args(['--set', str(tmp_path)])
    fu.use_set(args)
    assert args.src == str(tmp_path / 'unfound_levels.csv')


def test_find_unknown_levels_writes_a_table_of_failed_searches(tmp_path):
    """Every search failed: every measured column is blank, and the table
    must still be written, since it is what names the rows to run again."""
    import pandas as pd
    lev = dict(idx=986, label='x', J=1.5, parity='o', cfg='f5d2',
               E_calc=75159.9, W=226.8)
    row = fu.collect(lev, 1, None, False)
    cols = (fu.KEEP + fu.CORE + ['cowan_lid'] + fu.FLAGS + ['Note'])
    path = tmp_path / 'found.csv'
    fu.write(pd.DataFrame([row], columns=cols), str(path))
    back = pd.read_csv(path)
    assert back.loc[0, 'verdict'] == 'error'
