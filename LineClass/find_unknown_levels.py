#!/usr/bin/env python
"""Search for every promising unfound level, and tabulate what turns up.

Run from inside LineClass/:

    python find_unknown_levels.py                 # the whole promising list
    python find_unknown_levels.py --min-n-prom 8  # a shorter, safer list
    python find_unknown_levels.py --idx 742 913   # just these two rows


1.  What this does
------------------
``unfound_levels.py`` ranks the levels of Cowan's calculation that have never
been placed in the spectrum by ``n_prom``: the number of their calculated
transitions that would have been recorded on the plates if the level is where
the calculation puts it.  A level with one such transition cannot be found at
all - one line can be made to fit any energy - and a level with ten is one
whose position, if it is there, is heavily overdetermined.  Its table is
``unfound_levels.csv``, sorted by ``n_prom``, most promising first.

``level_positions.py --unknown ROW`` then does the search for one of them: it
scans the energy window the calculation allows, and reports every position in
it where the recorded lines are more likely with a level there than with none
(ln R > 0), with the audit's own verdict on each - ``firm``, ``weak`` or ``no
support``.

This module is the loop between the two.  It walks ``unfound_levels.csv`` from
the top, runs the search on each row, stops at the first row whose ``n_prom``
falls below ``--min-n-prom``, and collects the best position of each search
into one table, ``found_levels.csv``.  That table is the worklist: the rows
whose verdict is ``firm`` are the positions worth putting into IDEN2 and
looking at by eye.

Nothing here is written back into the pipeline.  The searches are read-only,
and the only files this module writes are its own output and its log.


2.  What is taken from each search
----------------------------------
Of the counts line -

    2 positions where ln R > 0: 1 firm, 1 weak, 0 no support

- the three numbers, whose sum is ``n_found``.  Two flags follow from them:

``unique_found``
    exactly one ``firm`` position and no ``weak`` one.  The search has found
    one place for the level and offers no competitor to it, which is the case
    worth looking at first: where two positions both carry support they are
    fed by the same free lines and only one of them can be right.

``not_found``
    no position at all where ln R > 0.  The window holds nothing: the recorded
    lines in it are no more likely with a level there than with none.

Of the table that follows the counts line, the top row only - the searches
sort their own table by verdict first and by ln R within a verdict, so the top
row is the best position the window offers:

``E_found``  where it is, cm^-1
``verdict``  ``firm`` / ``weak`` / ``no support``, by the audit's constants
``ln_R``     how much more likely the recorded lines are with a level there
``n_match``  transitions of the level that a recorded line falls on
``n_free``   how many of those lines no accepted transition already claims -
             the lines the level can take without taking them from a level
             that has already been found.  This is the number that decides a
             verdict: a position fed entirely by lines another level holds is
             not a position, it is a collision
``dE_o_c``   E_found minus the calculated energy, cm^-1

``--all-columns`` adds the rest of what the search prints: ``look`` (ln R after
the look-elsewhere correction for the width of the window scanned), ``z``
(dE_o_c in units of the configuration's own rms W), ``n_obs``, ``n_miss``,
``free_gain`` (what the free lines alone are worth in ln R) and ``top_share``
(the largest single line's share of the positive evidence - a position where
one line is most of the case is not one to trust).


3.  Note - positions that are somebody else's position
------------------------------------------------------
Each level is searched for on its own, and a search has no way of knowing that
the position it likes is one another level has already been given, or is being
given in the same run.  Two levels of the same parity but different J draw on
overlapping sets of lines, so the same free lines can carry both, and their
searches then land within a few hundredths of a wavenumber of each other.
Only one of them can be the level that is really there: the lines cannot be
spent twice.

The ``Note`` column says so, for every position within ``--dup-tol`` cm^-1
(0.1 by default) of another one OF THE SAME PARITY - two levels of opposite
parity at one energy are two different levels, not one collision.  It names
the other row of this table, or, worse, the level IDEN2 has already accepted
there, whose lines the new position would have to take.  A row whose own IDEN2
entry already sits at the position its search found is not a collision, and
the note says that instead.


4.  cowan_lid
-------------
The last column is the level number of the same level in Cowan's calculation,
``tp_E1_no_trials.xlsx``.  It is what ``new_levels.txt`` needs in order to give
a newly accepted level its calculated transitions, and it is resolved here the
way ``insert_new_level.py`` resolves it: the calculated level list is matched
to ``IDEN2/enlev.dat`` through ``cowan_gA.match_to_enlev``, which joins the two
lists by ``IDEN2/IDEN_level_ids.txt`` where a level has an identifier.  A row
that matches nothing leaves the column blank; it is never guessed at.
"""

import argparse
import os
import subprocess
import sys
import time

import pandas as pd

import cowan_gA
import output_files

HERE = os.path.dirname(os.path.abspath(__file__))
ENLEV = os.path.join(HERE, 'IDEN2', 'enlev.dat')
ID_MAP = os.path.join(HERE, 'IDEN2', 'IDEN_level_ids.txt')
IN_CSV = os.path.join(HERE, 'unfound_levels.csv')
OUT_CSV = os.path.join(HERE, 'found_levels.csv')
LOGFILE = os.path.join(HERE, 'find_unknown_levels.log')

MIN_N_PROM = 5       # where the list stops being worth the search
DUP_TOL = 0.1        # cm^-1 within which two positions are one position.  The
                     # scans step in 0.020 cm^-1 and two searches fed by
                     # overlapping line sets land a few hundredths apart.

# The seven columns of unfound_levels.csv that identify the level.
KEEP = ['idx', 'label', 'J', 'parity', 'cfg', 'E_calc', 'W']

# The columns taken from the top row of the search, in the order asked for.
CORE = ['E_found', 'verdict', 'n_found', 'ln_R', 'n_match', 'n_free',
        'dE_o_c']
EXTRA = ['look', 'z', 'n_obs', 'n_miss', 'free_gain', 'top_share']
FLAGS = ['unique_found', 'not_found']

COUNTS = 'positions where ln R > 0:'
EMPTY = 'no position in the window where ln R > 0'


# --------------------------------------------------------------------------
# reading one search
# --------------------------------------------------------------------------

def run_search(idx, extra_args=(), python=None):
    """Run level_positions.py --unknown IDX and hand back what it printed.

    The search is a plain subprocess rather than an import so that every level
    is searched by a model built from scratch: --unknown registers the level
    being looked for as a pseudo-level of the run, and a hundred registrations
    in one process is a hundred chances for one search to colour the next.  A
    few seconds per level is a cheap price for that independence.
    """
    cmd = [python or sys.executable, os.path.join(HERE, 'level_positions.py'),
           '--unknown', str(idx), '--drop-all-questionable']
    cmd.extend(extra_args)
    proc = subprocess.run(cmd, cwd=HERE, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True)
    return proc.returncode, proc.stdout


def block_of(text, idx):
    """The part of the output that belongs to IDEN2 row idx.

    Everything before ``IDEN2 row N`` is the model the search builds first -
    the width model, the background, the level list - which is the same for
    every level and is of no interest here.
    """
    at = text.find('\nIDEN2 row %d' % idx)
    if at < 0:
        return ''
    return text[at + 1:]


def parse_counts(block):
    """(n_firm, n_weak, n_no_support) from the counts line, or None.

    None means the search never got as far as counting: the row is not in
    enlev.dat, or it has no calculated transition to any level that has been
    found.  That is not the same as a search that counted zero, which is what
    ``no position in the window`` is, and the two must not be run together.
    """
    for line in block.splitlines():
        s = line.strip()
        if s.startswith(EMPTY):
            return (0, 0, 0)
        if COUNTS not in s:
            continue
        tail = s.split(COUNTS, 1)[1].split(';')[0]
        got = {}
        for part in tail.split(','):
            n, _, name = part.strip().partition(' ')
            got[name.strip()] = int(n)
        return (got.get('firm', 0), got.get('weak', 0),
                got.get('no support', 0))
    return None


def parse_top_row(block):
    """The first data row of the printed table, as {column: value}.

    The table is printed by DataFrame.to_string(index=False), so the header
    names the columns and every value is whitespace-separated.  The verdict is
    the last column and the only one that can hold a space ('no support'), so
    the other fields are split off the front and whatever is left is it.
    """
    lines = block.splitlines()
    for i, line in enumerate(lines):
        if COUNTS not in line:
            continue
        for j in range(i + 1, len(lines)):
            head = lines[j].split()
            if head[:2] == ['E', 'ln_R']:
                break
        else:
            return None
        if j + 1 >= len(lines):
            return None
        data = lines[j + 1].split()
        if len(data) < len(head):
            return None
        row = dict(zip(head[:-1], data[:len(head) - 1]))
        row[head[-1]] = ' '.join(data[len(head) - 1:])
        return row
    return None


def as_float(row, name):
    """One value of the top row as a number, or None where it is '-'."""
    if row is None:
        return None
    value = row.get(name)
    if value is None or value == '-':
        return None
    try:
        return float(value)
    except ValueError:
        return None


# --------------------------------------------------------------------------
# cowan_lid
# --------------------------------------------------------------------------

def cowan_lids(log=print):
    """{enlev.dat row: Cowan level number}, by insert_new_level.py's route."""
    trans = cowan_gA.read_transitions(log=lambda m: log('   ' + m))
    calc_levels = cowan_gA.levels(trans)
    enlev_levels = cowan_gA.read_enlev_levels(ENLEV)
    id_of_row = cowan_gA.read_id_map(ID_MAP)
    mapping, _report = cowan_gA.match_to_enlev(calc_levels, enlev_levels,
                                               id_of_row)
    return {row: lid for lid, row in mapping.items()}


# --------------------------------------------------------------------------
# positions that are the same position
# --------------------------------------------------------------------------

def accepted_levels():
    """Every level IDEN2 already holds a position for, with its parity.

    ``level_interchange.read_enlev`` reads IDEN2/enlev.dat into one row per
    calculated level: ``known`` says whether a position has been accepted for
    it, ``E_obs`` is that position, and the parity follows from the
    configuration through ``level_interchange.configuration_parity`` - the
    calculation's own label for the row, not a guess.
    """
    import level_interchange as li
    en = li.read_enlev(ENLEV)
    en = en[en['known']].copy()
    en['parity'] = [li.configuration_parity(c) for c in en['cfg']]
    return en


def add_notes(tab, tol, id_of_row, log=print):
    """Fill the Note column: positions that are somebody else's position.

    The search looks for one level at a time, and it has no way of knowing
    that the position it likes is one another level has already been given or
    is being given in the same run.  Two levels of the same parity but
    different J draw on overlapping sets of lines, so the same free lines can
    carry both, and the two searches then land within a few hundredths of a
    wavenumber of each other.  Only one of them can be the level that is
    really there: the lines cannot be spent twice.

    Two kinds of collision are reported, both within ``tol`` cm^-1 and both
    only between levels of the SAME parity - two levels of opposite parity at
    the same energy are two different levels, not one:

    * with another row of this table - the searches of this run;
    * with a level IDEN2 has already accepted, which is the worse of the two:
      the position is taken, and the lines that support it are lines that
      level already holds.

    A row whose own IDEN2 entry is already at the position it found is not a
    collision at all - the search has simply confirmed where the level sits -
    and says so.
    """
    accepted = accepted_levels()
    notes = []
    for i, row in tab.iterrows():
        e = row['E_found']
        found = []
        if pd.notna(e):
            for j, other in tab.iterrows():
                if j == i or pd.isna(other['E_found']):
                    continue
                if other['parity'] != row['parity']:
                    continue
                if abs(other['E_found'] - e) <= tol:
                    found.append('duplicate with row %d (J %s), %.3f cm^-1 '
                                 'away' % (other['idx'], other['J'],
                                           abs(other['E_found'] - e)))
            near = accepted[(accepted['parity'] == row['parity'])
                            & ((accepted['E_obs'] - e).abs() <= tol)]
            for _, lev in near.iterrows():
                # A row enlev.dat calls found need not be a level of the
                # pipeline: a position marked in IDEN2 by hand has no
                # level_id until it is brought in through new_levels.txt.
                # That is a state, not a failed lookup, and it is said so.
                got = id_of_row.get(int(lev['idx']))
                what = ('as %s' % got if got else
                        'which has no level_id yet: it is marked in IDEN2 '
                        'but is not a level of the pipeline')
                if int(lev['idx']) == int(row['idx']):
                    found.append('this is the position IDEN2 already holds '
                                 'for this row, %s' % what)
                else:
                    found.append('duplicate with the position IDEN2 holds '
                                 'for row %d (J %s), %s, %.3f cm^-1 away'
                                 % (int(lev['idx']), lev['J'], what,
                                    abs(lev['E_obs'] - e)))
        notes.append('; '.join(found))
    tab['Note'] = notes
    n = int(sum(1 for x in notes if x))
    log('%d of the %d positions found are somebody else\'s position too '
        '(within %.3f cm^-1)' % (n, int(tab['E_found'].notna().sum()), tol))
    return tab


# --------------------------------------------------------------------------
# the loop
# --------------------------------------------------------------------------

def shortest(path):
    """The path as this directory sees it, or in full if it is elsewhere.

    os.path.relpath raises on Windows when the two paths are on different
    drives, which an --out in a temporary directory readily is.
    """
    try:
        return os.path.relpath(path, HERE)
    except ValueError:
        return path


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--in', dest='src', default=IN_CSV,
                   help='the ranked list to walk (default %(default)s)')
    p.add_argument('--out', default=OUT_CSV,
                   help='where the table goes (default %(default)s)')
    p.add_argument('--min-n-prom', '--min_n_prom', dest='min_n_prom',
                   type=int, default=MIN_N_PROM, metavar='N',
                   help='stop at the first level with fewer than N promising '
                        'transitions (default %(default)s).  The list is '
                        'sorted by that count, so this is where it stops '
                        'being worth the search')
    p.add_argument('--idx', nargs='+', type=int, default=None,
                   metavar='IDEN2_ROW',
                   help='search these rows only, in this order, whatever '
                        'their n_prom')
    p.add_argument('--limit', type=int, default=None, metavar='N',
                   help='stop after N levels, whatever their n_prom')
    p.add_argument('--dup-tol', '--dup_tol', dest='dup_tol', type=float,
                   default=DUP_TOL, metavar='X',
                   help='two positions of the same parity within X cm^-1 of '
                        'each other are the same position, and each is noted '
                        'against the other (default %(default)s)')
    p.add_argument('--all-columns', action='store_true',
                   help='also keep look, z, n_obs, n_miss, free_gain and '
                        'top_share from the search')
    p.add_argument('--notes-only', action='store_true',
                   help='do no searching: read the table already written to '
                        '--out, work the Note column out again and write it '
                        'back.  For trying another --dup-tol on a run that '
                        'took ten minutes')
    p.add_argument('--log', default=LOGFILE,
                   help='where the full output of every search is kept '
                        '(default %(default)s); "-" keeps none')
    p.add_argument('--pass-through', nargs=argparse.REMAINDER, default=[],
                   metavar='ARG',
                   help='everything after this is handed to '
                        'level_positions.py unchanged')
    return p.parse_args(argv)


def select(src, args):
    """The rows to search, in the order they are to be searched."""
    if args.idx:
        known = set(src['idx'])
        missing = [i for i in args.idx if i not in known]
        if missing:
            raise SystemExit('%s has no row for %s'
                             % (args.src, ', '.join(map(str, missing))))
        todo = src.set_index('idx').loc[args.idx].reset_index()
    else:
        # The list is sorted by n_prom, so the first level below the threshold
        # ends it: this is a walk down a ranked list, not a filter.
        below = src.index[src['n_prom'] < args.min_n_prom]
        todo = src if not len(below) else src.loc[:below[0] - 1]
    if args.limit is not None:
        todo = todo.head(args.limit)
    return todo


def collect(lev, code, block, all_columns):
    """One output row, from the level and from what its search printed."""
    row = {name: lev[name] for name in KEEP}
    counts = parse_counts(block) if block else None
    top = parse_top_row(block) if block else None

    if counts is None:
        # The search never got as far as counting.  Say so in the verdict
        # rather than write a zero, which would read as a window searched and
        # found empty - the opposite of a search that did not happen.
        row['verdict'] = 'not searched' if code == 0 else 'error'
        row['n_found'] = None
        row['unique_found'] = False
        row['not_found'] = False
    else:
        n_firm, n_weak, n_none = counts
        row['n_found'] = n_firm + n_weak + n_none
        row['unique_found'] = bool(n_firm == 1 and n_weak == 0)
        row['not_found'] = row['n_found'] == 0
        row['verdict'] = (top or {}).get('verdict', '')
    row['E_found'] = as_float(top, 'E')
    row['ln_R'] = as_float(top, 'ln_R')
    row['n_match'] = as_float(top, 'n_match')
    row['n_free'] = as_float(top, 'n_free')
    row['dE_o_c'] = as_float(top, 'dE')
    if all_columns:
        for name in EXTRA:
            row[name] = as_float(top, name)
    return row


def write(tab, path):
    """The table to disk, with counts written as counts and LF endings."""
    # A count written as 2.0 reads as a measurement, and Int64 is the dtype
    # that also holds the blank of a level whose search never got as far as
    # counting.
    for name in ('n_found', 'n_match', 'n_free', 'n_obs', 'n_miss',
                 'cowan_lid'):
        if name in tab.columns:
            tab[name] = tab[name].astype('Int64')
    for name in ('E_calc', 'W', 'E_found', 'ln_R', 'dE_o_c', 'look', 'z',
                 'free_gain', 'top_share'):
        if name in tab.columns:
            tab[name] = tab[name].round(3)
    # newline='' keeps the file LF-only on Windows, as every text file of this
    # repository is.
    with open(path, 'w', newline='') as fh:
        tab.to_csv(fh, index=False, lineterminator='\n')


def main(argv=None):
    args = parse_args(argv)
    output_files.require_writable([args.out], 'output file')

    if args.notes_only:
        tab = pd.read_csv(args.out)
        add_notes(tab, args.dup_tol, cowan_gA.read_id_map(ID_MAP))
        write(tab, args.out)
        print('rewrote the Note column of %s' % shortest(args.out))
        return 0

    src = pd.read_csv(args.src)
    for name in KEEP + ['n_prom']:
        if name not in src.columns:
            raise SystemExit('%s has no %s column' % (args.src, name))
    todo = select(src, args)

    print('%d of the %d levels of %s to search'
          % (len(todo), len(src), os.path.basename(args.src)))
    if not len(todo):
        return 1
    print('resolving the Cowan level numbers')
    lid_of_row = cowan_lids()

    columns = (KEEP + CORE + (EXTRA if args.all_columns else [])
               + ['cowan_lid'] + FLAGS + ['Note'])
    out_rows = []
    transcript = []
    t0 = time.time()

    for n, (_, lev) in enumerate(todo.iterrows(), 1):
        idx = int(lev['idx'])
        code, text = run_search(idx, args.pass_through)
        block = block_of(text, idx)
        transcript.append('=' * 78)
        transcript.append('IDEN2 row %d   exit %d' % (idx, code))
        transcript.append('=' * 78)
        transcript.append(block or text)

        row = collect(lev, code, block, args.all_columns)
        row['cowan_lid'] = lid_of_row.get(idx)
        out_rows.append(row)
        print('%4d/%d  row %-5d n_prom %2d  %-12s n_found %-3s%s'
              % (n, len(todo), idx, int(lev['n_prom']), row['verdict'] or '-',
                 '-' if row['n_found'] is None else int(row['n_found']),
                 '  E = %.3f' % row['E_found']
                 if row['E_found'] is not None else ''))

    tab = pd.DataFrame(out_rows, columns=columns)
    print('')
    add_notes(tab, args.dup_tol, cowan_gA.read_id_map(ID_MAP))
    write(tab, args.out)
    if args.log != '-':
        with open(args.log, 'w', newline='') as fh:
            fh.write('\n'.join(transcript) + '\n')

    firm = tab['verdict'] == 'firm'
    print('')
    print('%d searched in %.0f s: %d firm, %d of them the only position with '
          'any support, %d with nothing in the window'
          % (len(tab), time.time() - t0, int(firm.sum()),
             int((firm & tab['unique_found']).sum()),
             int(tab['not_found'].sum())))
    print('wrote %s' % shortest(args.out))
    if args.log != '-':
        print('every search in full is in %s' % shortest(args.log))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
