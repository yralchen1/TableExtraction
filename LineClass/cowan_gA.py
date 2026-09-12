"""The Cowan transition probabilities, and how their levels line up with IDEN2.

WHAT THE FILE IS
================
``tp_E1_no_trials.xlsx`` is the electric-dipole (E1) transition list of the
Cowan calculation for Pr III: one row per calculated transition, 120274 of
them, among 1254 calculated levels.  It is the *source* of ``Icalc.xlsx`` -
every one of the 30206 rows of ``Icalc.xlsx`` is one of these rows, with the
same gA to five significant figures and the same percentage uncertainty - but
it is much larger, because ``Icalc.xlsx`` keeps only the transitions whose two
levels are both experimentally known, and most calculated levels are not.

The columns this module reads (numbered as in the spreadsheet, counting from
one) are

    2, 3, 4    conf1, term1, J1    the lower level; a ``*`` in term1 means the
                                   level has odd parity, otherwise even
    5, 6, 7    conf2, term2, J2    the upper level, marked the same way
    10, 11     E_low, E_up         calculated energies in kK (1 kK = 1000
                                   cm^-1); this module converts them to cm^-1
    12, 13     lid1, lid2          the calculation's own level numbers
    14, 15     id1, id2            the experimental level identifiers of
                                   Wyart's list, where the level is known
    17, 18     E_low_exp, E_up_exp Wyart's measured energies, in kK
    20         A(s-1)              the transition probability A
    29         std_A%              the percentage uncertainty of A, taken as
                                   the standard deviation over Monte-Carlo
                                   trials in which the Slater parameters and
                                   the dipole matrix elements were varied

Every other column belongs to unrelated research and is ignored.

The quantity the rest of the pipeline works in is not A but

    gA = (2*J2 + 1) * A ,

the transition probability weighted by the statistical weight of the upper
level, and that is what :func:`read_transitions` returns.

WHY A MODULE OF ITS OWN
=======================
Two programs need this table: ``sync_IDEN2.py``, which writes IDEN2's files
from it, and ``level_positions.py``, which needs the calculated strength of
transitions from levels that have never been found and therefore appear
nowhere in ``Icalc.xlsx``.  Reading a 29 MB spreadsheet takes about eleven
seconds, so :func:`read_transitions` keeps a comma-separated cache beside the
workbook and re-reads the workbook only when it is newer than the cache.

THE LEVEL NUMBERS, AND HOW THEY REACH IDEN2
===========================================
Three different numberings of the same levels are in play:

    lid          the Cowan calculation's level number, used in this
                 spreadsheet and nowhere else
    IDEN2 index  the row number of ``IDEN2/enlev.dat``, which is what
                 ``IDEN2/trans.dat`` refers to
    level_id     Wyart's experimental identifier, ``059003.000042`` and the
                 like, used by every other file of the pipeline

``IDEN2/IDEN_level_ids.txt`` maps the IDEN2 index to the level_id, but only
for the 594 levels that have been found; the other 659 have no level_id at
all, so for them there is nothing to join on.  :func:`match_to_enlev` supplies
the missing half.

It works by alignment rather than by matching.  Both lists are the same set of
levels in the same calculated order, but the two files were written from runs
that differ slightly - the calculated energies disagree by up to 5 cm^-1, and
in ten places two neighbouring levels are half a wavenumber apart and come out
in opposite order.  So the two energy-ordered lists are aligned the way two
versions of a text are aligned, by dynamic programming: the cost of pairing
two levels is the difference of their calculated energies, a pairing of two
different J values is forbidden outright, and leaving a level unpaired costs
``GAP``.  An order reversal shows up as two adjacent unpaired levels, and a
second pass pairs those off against each other by energy.

The result is then checked, not trusted: every one of the 594 levels whose
IDEN2 index and level_id are both known must come out paired with the
spreadsheet row carrying that same level_id.  A single disagreement raises
:class:`MatchError` rather than returning a mapping that is wrong somewhere.
"""

import io
import os
import re

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TP_FILE = os.path.join(HERE, 'tp_E1_no_trials.xlsx')
KK = 1000.0          # the spreadsheet's energies are in kK = 1000 cm^-1
GAP = 25.0           # cm^-1: the cost of leaving a level unpaired
CROSS_TOL = 20.0     # cm^-1: how far apart two levels left over from an order
                     # reversal may be and still be paired with each other

# The spreadsheet columns, by position (counting from zero), and the names
# this module gives them.  Positions rather than headers: the headers of the
# columns belonging to the unrelated research are not ours to depend on.
_COLS = [(1, 'conf1'), (2, 'term1'), (3, 'J1'),
         (4, 'conf2'), (5, 'term2'), (6, 'J2'),
         (9, 'E_low'), (10, 'E_up'),
         (11, 'lid1'), (12, 'lid2'),
         (13, 'id1'), (14, 'id2'),
         (16, 'E_low_exp'), (17, 'E_up_exp'),
         (19, 'A'), (28, 'u_gA_pct')]


class MatchError(Exception):
    """The calculated levels could not be lined up with ``enlev.dat``."""


def level_id(value):
    """Wyart's identifier as the rest of the pipeline spells it.

    The spreadsheet holds it as the number 59003.000042; every other file
    holds the string ``059003.000042``.  An empty cell gives ``''``.
    """
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ''
    return '0%012.6f' % float(value)


# ---------------------------------------------------------------------------
# The transitions
# ---------------------------------------------------------------------------
def read_transitions(path=TP_FILE, cache=True, log=None):
    """The E1 transition list, one row per transition.

    Columns: ``lid1, conf1, term1, J1, parity1`` and the same five for the
    upper level; ``E_low, E_up`` (cm^-1, calculated); ``E_low_exp, E_up_exp``
    (cm^-1, measured, NaN where the level is not known); ``id1, id2`` (Wyart
    identifiers as strings, ``''`` where the level is not known); ``A``, ``gA``
    and ``u_gA_pct``.

    ``parity`` is ``'odd'`` or ``'even'``, read off the asterisk the
    spreadsheet puts on the term of an odd level; the term itself is returned
    without it.

    With ``cache`` the table is kept beside the workbook as
    ``<name>.cache.csv`` and the workbook is re-read only when it is the newer
    of the two.  Reading the workbook takes about eleven seconds and the cache
    about half a second.
    """
    say = log or (lambda *_: None)
    cache_path = os.path.splitext(path)[0] + '.cache.csv'
    if cache and os.path.exists(cache_path) \
            and os.path.getmtime(cache_path) >= os.path.getmtime(path):
        say(f"reading {os.path.basename(cache_path)}")
        t = pd.read_csv(cache_path, dtype={'id1': str, 'id2': str})
        return t.fillna({'id1': '', 'id2': ''})

    say(f"reading {os.path.basename(path)} (this takes a few seconds)")
    raw = pd.read_excel(path, usecols=[c for c, _ in _COLS])
    # usecols returns the columns in spreadsheet order whatever order they
    # were asked for, so name them in that order.
    raw.columns = [n for _, n in sorted(_COLS)]

    t = pd.DataFrame({
        'lid1': raw['lid1'].astype(int),
        'conf1': raw['conf1'].astype(str),
        'term1': [s.rstrip('*') for s in raw['term1'].astype(str)],
        'parity1': ['odd' if s.endswith('*') else 'even'
                    for s in raw['term1'].astype(str)],
        'J1': raw['J1'].astype(float),
        'lid2': raw['lid2'].astype(int),
        'conf2': raw['conf2'].astype(str),
        'term2': [s.rstrip('*') for s in raw['term2'].astype(str)],
        'parity2': ['odd' if s.endswith('*') else 'even'
                    for s in raw['term2'].astype(str)],
        'J2': raw['J2'].astype(float),
        'E_low': raw['E_low'].astype(float) * KK,
        'E_up': raw['E_up'].astype(float) * KK,
        'E_low_exp': raw['E_low_exp'].astype(float) * KK,
        'E_up_exp': raw['E_up_exp'].astype(float) * KK,
        'id1': [level_id(v) for v in raw['id1']],
        'id2': [level_id(v) for v in raw['id2']],
        'A': raw['A'].astype(float),
        'u_gA_pct': raw['u_gA_pct'].astype(float),
    })
    # gA is A weighted by the statistical weight of the upper level.  The
    # uncertainty is a percentage of A, and multiplying by a constant leaves a
    # percentage unchanged, so u_gA_pct is std_A% as it stands.
    t['gA'] = (2.0 * t['J2'] + 1.0) * t['A']

    if cache:
        try:
            t.to_csv(cache_path, index=False, lineterminator='\n')
            say(f"  cached as {os.path.basename(cache_path)}")
        except OSError as exc:                       # a read-only directory
            say(f"  could not write the cache: {exc}")
    return t


def levels(trans):
    """One row per calculated level, in order of calculated energy.

    Columns ``lid, conf, term, parity, J, E_calc, E_exp, level_id``.  Every
    level of the calculation appears in some transition, as the lower level of
    one or the upper level of another, so the two ends of the transition list
    between them name them all.
    """
    def side(n):
        return pd.DataFrame({
            'lid': trans['lid%d' % n], 'conf': trans['conf%d' % n],
            'term': trans['term%d' % n], 'parity': trans['parity%d' % n],
            'J': trans['J%d' % n],
            'E_calc': trans['E_low' if n == 1 else 'E_up'],
            'E_exp': trans['E_low_exp' if n == 1 else 'E_up_exp'],
            'level_id': trans['id%d' % n]})

    out = pd.concat([side(1), side(2)]).drop_duplicates('lid')
    return out.sort_values(['E_calc', 'lid']).reset_index(drop=True)


# ---------------------------------------------------------------------------
# enlev.dat, read only for what the alignment needs
# ---------------------------------------------------------------------------
def read_enlev_levels(path):
    """``idx, E_calc, E_obs, J, label`` for every row of ``enlev.dat``.

    Only what the alignment and the report need; ``sync_IDEN2.py`` reads the
    file properly, through ``swap_line_assignments_IDEN.Enlev``, when it comes
    to writing it.
    """
    rows = []
    with io.open(path, encoding='latin-1', newline='') as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            head, _, tail = line.partition('/')
            label = tail.partition('/')[0]
            v = head.replace('*', ' ').split()
            rows.append((int(v[0]), float(v[1]), float(v[3]), float(v[5]),
                         label))
    return pd.DataFrame(rows, columns=['idx', 'E_calc', 'E_obs', 'J', 'label'])


def read_id_map(path):
    """``{IDEN2 index: level_id}`` from ``IDEN2/IDEN_level_ids.txt``."""
    t = pd.read_csv(path, sep='\t', dtype={'level_id': str})
    return dict(zip(t['IDEN_id'].astype(int), t['level_id']))


# ---------------------------------------------------------------------------
# Lining the two level lists up
# ---------------------------------------------------------------------------
def _align(e_left, j_left, e_right, j_right, gap=GAP):
    """Pair two energy-ordered level lists without crossing.

    Returns ``(pairs, left_over_left, left_over_right)``, the pairs as index
    pairs into the two lists.  Needleman-Wunsch with the cost of a pair the
    difference of the two energies, a pair of unequal J forbidden, and an
    unpaired level costing ``gap``.
    """
    n, m = len(e_left), len(e_right)
    dist = np.full((n + 1, m + 1), np.inf)
    move = np.zeros((n + 1, m + 1), dtype=np.int8)
    dist[0, 0] = 0.0
    dist[1:, 0] = gap * np.arange(1, n + 1)
    move[1:, 0] = 1
    dist[0, 1:] = gap * np.arange(1, m + 1)
    move[0, 1:] = 2
    for i in range(1, n + 1):
        cost = np.abs(e_right - e_left[i - 1])
        cost[j_right != j_left[i - 1]] = np.inf
        row, prev = dist[i], dist[i - 1]
        for j in range(1, m + 1):
            pair = prev[j - 1] + cost[j - 1]
            skip_l = prev[j] + gap
            skip_r = row[j - 1] + gap
            if pair <= skip_l and pair <= skip_r:
                row[j], move[i, j] = pair, 0
            elif skip_l <= skip_r:
                row[j], move[i, j] = skip_l, 1
            else:
                row[j], move[i, j] = skip_r, 2

    i, j = n, m
    pairs, spare_l, spare_r = [], [], []
    while i > 0 or j > 0:
        step = move[i, j]
        if step == 0:
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif step == 1:
            spare_l.append(i - 1)
            i -= 1
        else:
            spare_r.append(j - 1)
            j -= 1
    pairs.reverse()
    return pairs, spare_l, spare_r


def match_to_enlev(calc_levels, enlev_levels, id_map=None, tol=CROSS_TOL):
    """``{lid: IDEN2 index}``, checked against every level_id both files know.

    ``calc_levels`` is what :func:`levels` returns and ``enlev_levels`` what
    :func:`read_enlev_levels` returns.  ``id_map`` is
    ``{IDEN2 index: level_id}`` from :func:`read_id_map`; when it is given,
    every level it names must come out matched to the spreadsheet row carrying
    the same level_id, and :class:`MatchError` is raised if one does not.  A
    row whose spreadsheet identifier is blank is passed over: that is a level
    identified since the identifiers of the calculation were written down, so
    there is nothing there to contradict anything.

    Returns ``(mapping, report)``.  ``report`` is a dict with ``n_matched``,
    ``max_dE`` (the largest energy disagreement of a matched pair, in cm^-1),
    ``n_reordered`` (pairs recovered from an order reversal), ``unmatched_lid``
    and ``unmatched_idx``.
    """
    left = enlev_levels.sort_values('E_calc').reset_index(drop=True)
    right = calc_levels.reset_index(drop=True)      # already energy-ordered
    e_l = left['E_calc'].to_numpy(float)
    e_r = right['E_calc'].to_numpy(float)
    j_l = left['J'].to_numpy(float)
    j_r = right['J'].to_numpy(float)

    pairs, spare_l, spare_r = _align(e_l, j_l, e_r, j_r)

    # Two levels half a wavenumber apart can come out in opposite order in the
    # two files.  A no-crossing alignment cannot pair both, and leaves them as
    # two neighbouring unpaired levels; pair those off here, closest first.
    free_l, free_r = set(spare_l), set(spare_r)
    candidates = sorted((abs(e_l[a] - e_r[b]), a, b)
                        for a in spare_l for b in spare_r if j_l[a] == j_r[b])
    n_reordered = 0
    for gapsize, a, b in candidates:
        if a in free_l and b in free_r and gapsize <= tol:
            pairs.append((a, b))
            free_l.discard(a)
            free_r.discard(b)
            n_reordered += 1

    idx = left['idx'].to_numpy()
    lid = right['lid'].to_numpy()
    mapping = {int(lid[b]): int(idx[a]) for a, b in pairs}

    if id_map:
        wanted = right['level_id'].to_numpy()
        wrong = []
        for a, b in pairs:
            expect = id_map.get(int(idx[a]))
            # An identifier the spreadsheet does not carry belongs to a level
            # found since the calculation's identifiers were written down:
            # there is nothing to check it against, and its absence is not a
            # disagreement.
            if expect and wanted[b] and expect != wanted[b]:
                wrong.append('enlev.dat row %d is %s, but it was matched to '
                             'the calculated level %d, which the spreadsheet '
                             'calls %r'
                             % (idx[a], expect, lid[b], wanted[b] or None))
        missed = [idx[a] for a in free_l if id_map.get(int(idx[a]))]
        if missed:
            wrong.append('no calculated level was found for enlev.dat rows '
                         + ', '.join(str(v) for v in sorted(missed)))
        if wrong:
            raise MatchError('the calculated levels do not line up with '
                             'enlev.dat:\n  ' + '\n  '.join(wrong[:20]))

    dE = np.array([abs(e_l[a] - e_r[b]) for a, b in pairs]) if pairs \
        else np.zeros(0)
    report = {
        'n_matched': len(pairs),
        'max_dE': float(dE.max()) if len(dE) else 0.0,
        'n_reordered': n_reordered,
        'n_checked': sum(1 for a, _ in pairs if id_map and
                         id_map.get(int(idx[a]))) if id_map else 0,
        'unmatched_lid': sorted(int(lid[b]) for b in free_r),
        'unmatched_idx': sorted(int(idx[a]) for a in free_l),
    }
    return mapping, report
