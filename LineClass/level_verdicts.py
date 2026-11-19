"""The user's verdicts on the levels the audit questions.

`LineClass/level_verdicts.txt` (files.level_verdicts; 2026-11-19) is the one
place where the decisions on questionable levels are kept.  It is hand-kept
and tab-delimited, with the columns

    level_id   the level
    verdict    `firm`: the position stands; `?`: the level stays where it
               is, questioned
    E_alt      the energy of an alternate position looked at and denied
               (one row per alternate; a level has one verdict)
    date       when
    note       why

and may be saved by Excel.  level_positions.py applies the verdicts over the
audit's question marks (`mark`), and so does --drop-all-questionable;
prepare_audit.py shows them in the review copy beside the position the audit
now suggests.
"""
import csv
import io
import os
from typing import NamedTuple

import numpy as np

import output_files

#: cm^-1: an alternate position within this of an E_alt is the one looked at
VERDICT_TOL = 0.5
#: the verdicts the file may give
VERDICTS = ('firm', '?')
#: the file's columns, in the order a new file is written with
COLUMNS = ('level_id', 'verdict', 'E_alt', 'date', 'note')


def _read_bytes(path):
    with open(path, 'rb') as fh:
        return fh.read()


class Verdict(NamedTuple):
    verdict: str          # 'firm' or '?'
    e_alt: tuple          # the alternate positions looked at and denied
    note: str


def read(path):
    """`{level_id: Verdict}` from the user's verdicts file; {} if `path` is
    '' or the file is absent.

    The file is hand-kept and tab-delimited, with a header naming at least
    level_id and verdict, and optionally E_alt, date and note.  Excel may
    have saved it: cells are unquoted first, blank cells at a row's end are
    dropped, a row whose first cell starts with '#' is a comment.  A level
    may have several rows, one per alternate position denied, but one
    verdict.  `firm`: the position stands; the audit's question mark is
    lifted while the alternate it finds is one of E_alt (VERDICT_TOL), or
    there is none.  `?`: the level stays where it is, questioned, whatever
    the audit finds."""
    if not path or not os.path.exists(path):
        return {}
    raw = output_files.read_retry(_read_bytes, path)
    try:
        text = raw.decode('utf-8-sig')
    except UnicodeDecodeError:
        text = raw.decode('cp1252')
    head, out = None, {}
    for n, cells in enumerate(csv.reader(io.StringIO(text, newline=''),
                                         delimiter='\t'), start=1):
        cells = [c.strip() for c in cells]
        while cells and cells[-1] == '':
            cells.pop()
        if not cells or cells[0].startswith('#'):
            continue
        if head is None:
            head = cells
            if 'level_id' not in head or 'verdict' not in head:
                raise SystemExit(f"{path}: the header must name the columns "
                                 f"level_id and verdict; it reads {head}")
            continue
        row = dict(zip(head, cells + [''] * (len(head) - len(cells))))
        lid, verdict = row['level_id'], row['verdict'].lower()
        if verdict not in VERDICTS:
            raise SystemExit(f"{path}, line {n}: verdict {row['verdict']!r} "
                             f"for {lid}; it must be one of {VERDICTS}")
        e_alt = ()
        if row.get('E_alt', ''):
            try:
                e_alt = (float(row['E_alt']),)
            except ValueError:
                raise SystemExit(f"{path}, line {n}: E_alt "
                                 f"{row['E_alt']!r} of {lid} is not a "
                                 f"number")
        old = out.get(lid)
        if old is not None and old.verdict != verdict:
            raise SystemExit(f"{path}, line {n}: {lid} is given the verdicts "
                             f"{old.verdict!r} and {verdict!r}")
        note = row.get('note', '')
        out[lid] = (Verdict(verdict, old.e_alt + e_alt,
                            '; '.join(x for x in (old.note, note) if x))
                    if old else Verdict(verdict, e_alt, note))
    return out


def mark(row, rec, tol=VERDICT_TOL):
    """`(question, words)`: the question mark ('?' or '') the user's verdict
    `rec` gives a report row, and what q_why says of it.

    A `?` verdict always marks.  A `firm` one lifts the mark unless the row
    has an alternate position (n_alt > 0) at an energy E + dE_alt that is no
    E_alt of the verdict: an alternate not yet looked at."""
    if rec.verdict == '?':
        return '?', 'user ?'

    def num(k):
        try:
            return float(row.get(k))
        except (TypeError, ValueError):
            return float('nan')
    n_alt, E, dE = num('n_alt'), num('E'), num('dE_alt')
    if not np.isfinite(n_alt):
        n_alt = 1.0 if np.isfinite(dE) else 0.0
    if n_alt <= 0:
        return '', 'firm (user)'
    e_new = E + dE
    if np.isfinite(e_new) and any(abs(e_new - e) < tol for e in rec.e_alt):
        return '', 'firm (user, alternate %.4f denied)' % e_new
    return '?', ('firm (user), but the alternate at %.4f is not one denied'
                 % e_new if np.isfinite(e_new) else
                 'firm (user), but an alternate position is not placed')


def nearest(rec, e):
    """The E_alt of `rec` nearest the energy `e`, or the first one if `e`
    is not a number; None if the verdict lists none."""
    if not rec.e_alt:
        return None
    try:
        e = float(e)
    except (TypeError, ValueError):
        e = float('nan')
    if not np.isfinite(e):
        return rec.e_alt[0]
    return min(rec.e_alt, key=lambda x: abs(x - e))
