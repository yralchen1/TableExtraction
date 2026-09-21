"""Tests for the check that an output file can be written before it is.

A workbook open in Excel is locked against writing on Windows, so a run whose
output is open fails at the end, after all of the work.  output_files makes
the same test at the start of the run instead.  Excel cannot be started in a
test, so a refused write is produced instead by taking away the file's write
permission: the check goes by the refusal, which is the same either way.

Run from the LineClass directory:  python -m pytest tests -q
"""
import os
import stat
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import output_files


def test_a_free_file_is_writable(tmp_path):
    path = tmp_path / 'line_classifications.csv'
    path.write_text('wn_obs\n100.0\n', encoding='utf-8')
    assert output_files.why_unwritable(str(path)) is None
    assert output_files.unwritable([str(path)]) == []
    output_files.require_writable([str(path)])       # does not raise


def test_a_file_that_does_not_exist_yet_is_writable(tmp_path):
    """The first run of all writes files that are not there yet."""
    assert output_files.why_unwritable(str(tmp_path / 'new.xlsx')) is None


def test_a_missing_folder_is_reported(tmp_path):
    path = tmp_path / 'no' / 'such' / 'place' / 'out.csv'
    reason = output_files.why_unwritable(str(path))
    assert reason is not None and 'does not exist' in reason


def test_a_folder_in_place_of_a_file_is_reported(tmp_path):
    (tmp_path / 'out.csv').mkdir()
    reason = output_files.why_unwritable(str(tmp_path / 'out.csv'))
    assert reason is not None and 'folder' in reason


def test_a_file_that_cannot_be_written_is_reported(tmp_path):
    """A refused write is what Excel's lock looks like from here.

    Excel cannot be started in a test, so the same refusal is produced by
    taking the write permission away: in both cases the open for update
    raises PermissionError, which is all the check goes by.
    """
    path = tmp_path / 'line_classifications.xlsx'
    path.write_bytes(b'PK')
    os.chmod(str(path), stat.S_IREAD)
    try:
        reason = output_files.why_unwritable(str(path))
    finally:
        os.chmod(str(path), stat.S_IREAD | stat.S_IWRITE)
    assert reason is not None and 'open in another program' in reason


def test_excel_s_owner_file_is_named_when_it_is_there(tmp_path):
    """`~$name.xlsx` beside a held workbook says which program holds it."""
    path = tmp_path / 'line_classifications.xlsx'
    path.write_bytes(b'PK')
    (tmp_path / '~$line_classifications.xlsx').write_bytes(b'')
    os.chmod(str(path), stat.S_IREAD)
    try:
        reason = output_files.why_unwritable(str(path))
    finally:
        os.chmod(str(path), stat.S_IREAD | stat.S_IWRITE)
    assert 'Excel' in reason and '~$line_classifications.xlsx' in reason


def test_a_stale_owner_file_alone_does_not_stop_a_run(tmp_path):
    """Excel leaves `~$` files behind when it is killed; they are not a lock."""
    path = tmp_path / 'line_classifications.xlsx'
    path.write_bytes(b'PK')
    (tmp_path / '~$line_classifications.xlsx').write_bytes(b'')
    assert output_files.why_unwritable(str(path)) is None


def test_an_unreadable_path_stops_the_run(tmp_path):
    """require_writable raises SystemExit, and names the file in the message."""
    (tmp_path / 'out.csv').mkdir()
    with pytest.raises(SystemExit) as exc:
        output_files.require_writable([str(tmp_path / 'out.csv')])
    assert 'out.csv' in str(exc.value)


def test_empty_paths_are_ignored(tmp_path):
    """A run that writes no csv passes '' for it and is not stopped."""
    output_files.require_writable(['', None, str(tmp_path / 'fresh.csv')])


def test_a_path_named_twice_is_tested_once(tmp_path):
    (tmp_path / 'out.csv').mkdir()
    bad = output_files.unwritable([str(tmp_path / 'out.csv'),
                                   str(tmp_path / 'out.csv')])
    assert len(bad) == 1


def test_with_twin_names_the_xlsx_beside_the_csv():
    assert output_files.with_twin('report.csv') == ['report.csv', 'report.xlsx']
    assert output_files.with_twin('') == []


def test_the_message_names_every_file_that_is_held(tmp_path):
    for name in ('a.csv', 'b.xlsx'):
        (tmp_path / name).mkdir()
    with pytest.raises(SystemExit) as exc:
        output_files.require_writable([str(tmp_path / 'a.csv'),
                                       str(tmp_path / 'b.xlsx')])
    message = str(exc.value)
    assert 'a.csv' in message and 'b.xlsx' in message
    assert 'run again' in message


# ---------------------------------------------------------------------------
# read_retry: a read of a file another program holds open for a moment
# ---------------------------------------------------------------------------

def test_read_retry_returns_the_value_on_the_first_attempt():
    calls = []
    got = output_files.read_retry(lambda p: calls.append(p) or 'ok', 'f.txt')
    assert got == 'ok' and calls == ['f.txt']


def test_read_retry_waits_out_a_lock_that_lets_go(monkeypatch):
    """Two refusals then success: the value comes back and nothing is raised."""
    monkeypatch.setattr(output_files.time, 'sleep', lambda _s: None)
    state = {'n': 0}

    def read(path):
        state['n'] += 1
        if state['n'] < 3:
            raise PermissionError(13, 'Permission denied')
        return 'contents'

    notes = []
    assert output_files.read_retry(read, 'IDEN_level_ids.txt',
                                   log=notes.append) == 'contents'
    assert state['n'] == 3
    # One note per repeat, each naming the file.
    assert len(notes) == 2
    assert all('IDEN_level_ids.txt' in note for note in notes)


def test_read_retry_gives_up_and_names_the_file(monkeypatch):
    monkeypatch.setattr(output_files.time, 'sleep', lambda _s: None)

    def read(_path):
        raise PermissionError(13, 'Permission denied')

    with pytest.raises(PermissionError) as exc:
        output_files.read_retry(read, os.path.join('IDEN2', 'enlev.dat'),
                                tries=3)
    message = str(exc.value)
    assert 'enlev.dat' in message and '3 attempts' in message


def test_read_retry_does_not_wait_for_an_error_waiting_cannot_cure(monkeypatch):
    """A missing file or a bad format is raised at once, not retried."""
    slept = []
    monkeypatch.setattr(output_files.time, 'sleep', slept.append)

    def read(_path):
        raise FileNotFoundError('no such file')

    with pytest.raises(FileNotFoundError):
        output_files.read_retry(read, 'gone.csv')
    assert slept == []
