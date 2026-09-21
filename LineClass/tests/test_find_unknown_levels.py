"""Tests for how find_unknown_levels.py reports a search that failed.

A walk over several hundred levels runs each search as its own subprocess.
When one of them dies - a file held open by another program is the usual
reason - the table gets the verdict ``error``, which says only that the
subprocess exited non-zero.  These tests cover the two things that make such a
failure legible: the exception line kept for the console, and the verdict
itself.

Run from the LineClass directory:  python -m pytest tests -q
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import find_unknown_levels as ful

TRACEBACK = """\
reading tp_E1_no_trials.cache.csv
Traceback (most recent call last):
  File "level_positions.py", line 4672, in main
    en, trans, mapping = uf.read_theory(log=log)
PermissionError: [Errno 13] Permission denied
"""

LEVEL = {'idx': 951, 'label': 'f5d2 ~3P2F', 'J': 2.5, 'parity': 'o',
         'cfg': 'f5d2', 'E_calc': 83140.8, 'W': 228.2}


def test_failure_line_is_the_exception():
    assert ful.failure_line(TRACEBACK) == \
        'PermissionError: [Errno 13] Permission denied'


def test_failure_line_survives_an_empty_output():
    assert 'printed nothing' in ful.failure_line('')
    assert 'printed nothing' in ful.failure_line('\n  \n')


def test_a_failed_search_is_an_error_not_a_level_without_positions():
    """Exit non-zero with no counts line: the level was never searched."""
    row = ful.collect(LEVEL, 1, ful.block_of(TRACEBACK, 951), False)
    assert row['verdict'] == 'error'
    assert row['n_found'] is None          # not 0: nothing was counted
    assert row['not_found'] is False
    assert row['unique_found'] is False


def test_a_search_that_exited_cleanly_without_counting_is_not_an_error():
    row = ful.collect(LEVEL, 0, '', False)
    assert row['verdict'] == 'not searched'
    assert row['n_found'] is None


def test_an_empty_window_is_a_count_of_zero():
    block = ('\nIDEN2 row 951  f5d2 ~3P2F  J = 2.5  not found\n'
             '  ' + ful.EMPTY + '\n')
    row = ful.collect(LEVEL, 0, ful.block_of(block, 951), False)
    assert row['n_found'] == 0
    assert row['not_found'] is True
    assert row['verdict'] != 'error'


# ---------------------------------------------------------------------------
# --merge: a re-run of a few levels corrects the table instead of replacing it
# ---------------------------------------------------------------------------

COLUMNS = ['idx', 'verdict', 'n_found', 'Note']


def _tab(rows):
    return pd.DataFrame(rows, columns=COLUMNS)


def test_merge_replaces_the_row_in_place_and_keeps_the_order():
    old = _tab([{'idx': 951, 'verdict': 'error', 'n_found': None, 'Note': ''},
                {'idx': 513, 'verdict': 'weak', 'n_found': 2, 'Note': ''},
                {'idx': 691, 'verdict': 'error', 'n_found': None, 'Note': ''}])
    new = _tab([{'idx': 691, 'verdict': 'firm', 'n_found': 1, 'Note': ''}])
    got = ful.merge_into(old, new, COLUMNS)
    assert list(got['idx']) == [951, 513, 691]
    assert list(got['verdict']) == ['error', 'weak', 'firm']


def test_merge_adds_a_level_the_table_did_not_have():
    old = _tab([{'idx': 513, 'verdict': 'weak', 'n_found': 2, 'Note': ''}])
    new = _tab([{'idx': 209, 'verdict': 'firm', 'n_found': 1, 'Note': ''}])
    got = ful.merge_into(old, new, COLUMNS)
    assert list(got['idx']) == [513, 209]
