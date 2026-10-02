"""The lock on a working set (config.require_unlocked).

A set is locked by `locked = true` in its own lineclass_config.toml.  These
tests check that:

  * the lock belongs to the file that says it - a configuration that
    inherits a locked one is not locked;
  * a file belongs to the nearest directory above it that holds a
    configuration, so the baseline's IDEN2 is the baseline's and iter's is
    iter's;
  * a run that would write a locked set stops before it writes, naming the
    set and the files, and --unlock (or the environment variable that passes
    it on to the chain) lets it through;
  * classify_lines.py writes no table of another set than its
    configuration's, so a set that forgot to name its own output cannot
    overwrite its parent's;
  * every program that writes a set's files asks: classify_lines.py,
    make_LOPT_input.py and sync_IDEN2.py, and insert_new_level.py,
    move_level.py and discard_level.py before they start - and these three
    leave the shared hand-kept files out of the question.

Nothing here writes the real baseline: each program is stopped by a stub of
output_files.require_writable, which every one of them calls only after the
lock, and the run is expected to stop at the lock before reaching it.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import config  # noqa: E402
import output_files  # noqa: E402


class Reached(Exception):
    """The run got past the lock to the writability check."""


@pytest.fixture(autouse=True)
def no_override(monkeypatch):
    monkeypatch.delenv(config.UNLOCK_ENV, raising=False)


@pytest.fixture
def stop_at_writable(monkeypatch):
    def stub(*a, **k):
        raise Reached()
    monkeypatch.setattr(output_files, 'require_writable', stub)


@pytest.fixture
def sets(tmp_path):
    """A locked base and an unlocked set inheriting it, as LineClass/ and
    iter/ are."""
    base = tmp_path / 'base'
    (base / 'IDEN2').mkdir(parents=True)
    (base / 'lineclass_config.toml').write_text('locked = true\n')
    child = base / 'iter'
    (child / 'IDEN2').mkdir(parents=True)
    (child / 'lineclass_config.toml').write_text(
        'inherit = "../lineclass_config.toml"\n')
    return base, child


def test_the_baseline_is_locked_and_iter_is_not():
    assert config.is_locked(ROOT)
    assert not config.is_locked(os.path.join(ROOT, 'iter'))
    # and iter still reads the rest of the baseline's configuration
    assert config.load(os.path.join(ROOT, 'iter', 'lineclass_config.toml'))


def test_a_lock_is_not_inherited(sets):
    base, child = sets
    assert config.is_locked(str(base))
    assert not config.is_locked(str(child))
    assert 'locked' not in config._read(str(child / 'lineclass_config.toml'))


def test_a_file_belongs_to_the_nearest_set_above_it(sets, tmp_path):
    base, child = sets
    assert config.set_of(str(base / 'IDEN2' / 'enlev.dat')) == str(base)
    assert config.set_of(str(child / 'IDEN2' / 'enlev.dat')) == str(child)
    assert config.set_of(str(child / 'LOPT.par')) == str(child)
    assert config.set_of(str(tmp_path / 'x.csv')) is None


def test_a_run_on_a_locked_set_stops_and_names_it(sets):
    base, child = sets
    target = str(base / 'IDEN2' / 'enlev.dat')
    with pytest.raises(SystemExit) as exc:
        config.require_unlocked([str(child / 'LOPT.par'), target], 'prog.py')
    text = str(exc.value)
    assert 'prog.py: this run would write a locked set' in text
    assert str(base) in text and target in text
    assert str(child / 'LOPT.par') not in text
    assert config.require_unlocked([str(child / 'LOPT.par')], 'prog.py') == []


def test_unlock_lets_it_through(sets, monkeypatch, capsys):
    base, _child = sets
    target = [str(base / 'line_classifications.csv')]
    assert config.require_unlocked(target, 'prog.py', unlock=True) \
        == [str(base)]
    assert 'writing the locked set' in capsys.readouterr().out
    monkeypatch.setenv(config.UNLOCK_ENV, '1')
    assert config.require_unlocked(target, 'prog.py') == [str(base)]


def test_a_lock_that_is_not_a_boolean_is_an_error(tmp_path):
    (tmp_path / 'lineclass_config.toml').write_text('locked = "yes"\n')
    with pytest.raises(config.ConfigError, match='true or false'):
        config.is_locked(str(tmp_path))


# --- the programs ---------------------------------------------------------------
def test_classify_lines_refuses_the_baseline(stop_at_writable, monkeypatch):
    import classify_lines as CL
    monkeypatch.setattr(CL, 'CFG', config.load(
        os.path.join(ROOT, 'lineclass_config.toml')))
    monkeypatch.setattr(CL, 'OUTPUT_FILE',
                        os.path.join(ROOT, 'line_classifications.xlsx'))
    monkeypatch.setattr(CL, 'OUTPUT_CSV',
                        os.path.join(ROOT, 'line_classifications.csv'))
    monkeypatch.setattr(CL, 'UNLOCK', False)
    with pytest.raises(SystemExit, match='locked set'):
        CL.main()
    monkeypatch.setattr(CL, 'UNLOCK', True)
    with pytest.raises(Reached):
        CL.main()


def test_classify_lines_writes_iter(stop_at_writable, monkeypatch):
    import classify_lines as CL
    it = os.path.join(ROOT, 'iter')
    monkeypatch.setattr(CL, 'CFG', config.load(
        os.path.join(it, 'lineclass_config.toml')))
    monkeypatch.setattr(CL, 'OUTPUT_FILE',
                        os.path.join(it, 'line_classifications.xlsx'))
    monkeypatch.setattr(CL, 'OUTPUT_CSV',
                        os.path.join(it, 'line_classifications.csv'))
    monkeypatch.setattr(CL, 'UNLOCK', False)
    with pytest.raises(Reached):
        CL.main()


def test_own_output_follows_the_nearest_set(sets, tmp_path):
    base, child = sets
    cfg = str(child / 'lineclass_config.toml')
    config.require_own_output(cfg, [str(child / 'out.csv')], 'prog.py')
    # the parent's table, inherited by a child that did not name its own
    with pytest.raises(SystemExit) as exc:
        config.require_own_output(cfg, [str(base / 'out.csv')], 'prog.py')
    assert "would write another set's files" in str(exc.value)
    assert str(base / 'out.csv') in str(exc.value)
    # and the other way round: the parent writing into a child's directory
    with pytest.raises(SystemExit):
        config.require_own_output(str(base / 'lineclass_config.toml'),
                                  [str(child / 'out.csv')], 'prog.py')
    # a configuration outside any set, as a test's temporary one is, writes
    # beside itself
    other = tmp_path / 'loose'
    other.mkdir()
    config.require_own_output(str(other / 'cfg.toml'),
                              [str(other / 'out.csv')], 'prog.py')


def test_a_set_that_forgot_its_output_cannot_write_its_parents(
        stop_at_writable, monkeypatch, tmp_path):
    """An iter_hfs/ beside iter/ whose configuration does not name its own
    output inherits iter's table; classify_lines.py stops before writing."""
    import classify_lines as CL
    d = tmp_path / 'iter_hfs'
    d.mkdir()
    it = os.path.join(ROOT, 'iter', 'lineclass_config.toml').replace('\\', '/')
    (d / 'lineclass_config.toml').write_text(
        f'inherit = "{it}"\n[hfs]\napply = true\n')
    cfg = config.load(str(d / 'lineclass_config.toml'))
    assert config.set_of(cfg.output_csv) == os.path.join(ROOT, 'iter')
    monkeypatch.setattr(CL, 'CFG', cfg)
    monkeypatch.setattr(CL, 'OUTPUT_FILE', cfg.output_file)
    monkeypatch.setattr(CL, 'OUTPUT_CSV', cfg.output_csv)
    with pytest.raises(SystemExit, match="another set's files"):
        CL.main()
    # named, it goes through to the writability check
    (d / 'lineclass_config.toml').write_text(
        f'inherit = "{it}"\n[files]\noutput = "t.xlsx"\n'
        'output_csv = "t.csv"\n[hfs]\napply = true\n')
    cfg = config.load(str(d / 'lineclass_config.toml'))
    monkeypatch.setattr(CL, 'CFG', cfg)
    monkeypatch.setattr(CL, 'OUTPUT_FILE', cfg.output_file)
    monkeypatch.setattr(CL, 'OUTPUT_CSV', cfg.output_csv)
    with pytest.raises(Reached):
        CL.main()


def test_make_lopt_input_refuses_the_baseline(stop_at_writable):
    import make_LOPT_input as MLI
    table = os.path.join(ROOT, 'line_classifications.csv')
    with pytest.raises(SystemExit, match='locked set'):
        MLI.main(['--classifications', table])
    with pytest.raises(Reached):
        MLI.main(['--classifications', table, '--unlock'])
    with pytest.raises(Reached):
        MLI.main(['--classifications',
                  os.path.join(ROOT, 'iter', 'line_classifications.csv')])


def test_sync_iden2_refuses_the_baseline(stop_at_writable):
    import sync_IDEN2 as S
    with pytest.raises(SystemExit, match='locked set'):
        S.main(['--set', ROOT])
    with pytest.raises(Reached):
        S.main(['--set', ROOT, '--unlock'])
    with pytest.raises(Reached):
        S.main(['--set', os.path.join(ROOT, 'iter')])


def test_a_set_without_its_own_iden2_cannot_write_the_baselines(
        stop_at_writable, tmp_path):
    """sync_IDEN2 looks for a set's IDEN2 in the set and then in the project:
    the lock is what keeps such a set from writing the baseline's."""
    import shutil
    import sync_IDEN2 as S
    d = tmp_path / 'noiden2'
    d.mkdir()
    for name in ('LOPT_output_levels.txt', 'LOPT_input_lines.txt'):
        shutil.copy(os.path.join(ROOT, 'iter', name), d / name)
    with pytest.raises(SystemExit, match='locked set'):
        S.main(['--set', str(d), '--allow-mixed'])


@pytest.mark.parametrize('module, argv', [
    ('insert_new_level', ['--iden2-row', '658']),
    ('insert_new_level', ['--undo']),
    ('move_level', ['--undo']),
    ('discard_level', ['--undo', '--iden2-row', '658']),
])
def test_the_level_tools_refuse_the_baseline(module, argv, monkeypatch):
    import importlib
    import insert_new_level as INL
    monkeypatch.setattr(INL, 'UNLOCK', False)
    tool = importlib.import_module(module)
    with pytest.raises(SystemExit, match='locked set'):
        tool.main(argv)
    assert INL.UNLOCK is False


def test_the_level_tools_leave_the_shared_files_out(sets, monkeypatch):
    """With --set iter the ledger and new_levels.txt, which live in the
    baseline's directory, are written; they are not the baseline's own."""
    import insert_new_level as INL
    monkeypatch.setattr(INL, 'UNLOCK', False)
    INL.require_unlocked([INL.NEW_LEVELS, INL.LINE_DECISIONS],
                         'insert_new_level.py', False)
    with pytest.raises(SystemExit, match='locked set'):
        INL.require_unlocked([INL.NEW_LEVELS, INL.CLASSIFICATIONS],
                             'insert_new_level.py', False)


def test_unlock_reaches_the_chain(monkeypatch):
    import subprocess
    import insert_new_level as INL
    seen = {}

    class Done:
        returncode, stdout = 0, ''

    def fake(cmd, cwd, env, **k):
        seen.update(env)
        return Done()
    monkeypatch.setattr(subprocess, 'run', fake)
    monkeypatch.setattr(INL, 'UNLOCK', False)
    INL.require_unlocked([INL.CLASSIFICATIONS], 'insert_new_level.py', True)
    INL.run(['x'], lambda s: None)
    assert seen.get(config.UNLOCK_ENV) == '1'
