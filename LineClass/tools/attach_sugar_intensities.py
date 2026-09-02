#!/usr/bin/env python3
"""Copy the intensities Sugar printed into a column of the working line list.

Why this tool exists
--------------------
The working line list ``Pr3_lines.xlsx`` carries a column ``Icor``: an
intensity that has already been put on a common scale by a wavelength-dependent
*correction*, because the original measurements were made on several
photographic plates whose sensitivity differed and varied along the plate.  The
numbers Sugar himself printed in his two papers - small whole numbers such as
1, 3, 20, 500, 9000, one per line - were the starting point of that column, but
they were not kept anywhere in the working file.  Everything downstream
therefore had to *undo* the correction to get back to them, and could only do
so as accurately as the correction file was.

This tool puts the printed numbers themselves into the working file, in a new
column (``Iorig`` by default), so that the intensity scale can afterwards be
rebuilt from the measurements rather than from a reconstruction of them.

How lines are matched
---------------------
Not by wavelength.  Sugar's papers print, for each line, both a wavelength and
a wavenumber (a wavenumber is 1e8 divided by the vacuum wavelength in
angstroms, and is measured in reciprocal centimetres, cm^-1), and the two do
not always agree to the last digit he printed.  The extraction workbooks
therefore carry a reconciled value - a weighted mean of the two - and it is
that value, and only that, which the working line list was built from.  Each
source file is told which column holds it:

    --source "Sugar1969=file.xlsm:Table 1:wn_mean:Intensity"
             ^name       ^workbook ^sheet   ^wavenumber ^intensity

And not by the ``ref`` column of the working line list either.  That column
says whose *classification* a row carries, which is a different question from
who *measured* the line: the rows credited to Wyart were measured by Sugar like
all the others, and their intensities stand in his workbooks.  Which workbook
holds a given line is decided by the line's wavenumber, because the two
recordings cover disjoint stretches of the spectrum - in this project
everything above about 47500 cm^-1 was measured in 1969 and everything below it
in 1974.  The tool needs no boundary to be configured: it looks a line up in
every source given and keeps the nearest wavenumber found in any of them, which
comes out the same as the boundary rule whenever the sources are disjoint, and
picks the better of the two where they are not.  The name on the left of a
``--source`` is then only a label, used in the report and in the audit file.

A column may be named by its header (letter case, surrounding blanks and
embedded line breaks are ignored) or by its 1-based number, which is the safer
choice when the header contains a line break.

The match is rejected if the nearest source wavenumber is farther away than
``--tol`` (default 1e-4 cm^-1, far below any real difference between two lines:
the agreement in practice is at the level of 1e-10, the last binary digit of a
stored number).  One source line may legitimately serve several rows of the
working list, because a single observed line that has two candidate
identifications appears there once per candidate.

What is written
---------------
Only the one column, and only into a copy of the workbook made after the
original has been backed up (``*.bak_iorig.xlsx``).  The column is created at
the right-hand end if it is not there yet and overwritten in place if it is, so
running the tool again never adds a second copy of it.  Rows that find no
source line within ``--tol`` are left as they were; nothing else in the
workbook is touched.  A CSV audit listing every row, its wavenumber, the source
it was matched to, the source row, the distance between the two wavenumbers and
the value copied is written beside it.

Caution
-------
Rewriting an .xlsx with openpyxl discards the cached results of any array
formula in it (the formulas survive and Excel recomputes them when the file is
next opened) and can change the last binary digit of a stored number.
"""

import argparse
import bisect
import csv
import os
import shutil
import sys

import openpyxl


# ----------------------------------------------------------------- utilities

def norm(s):
    """A header reduced to what is worth comparing: lower case, no blanks."""
    return ''.join(str(s).split()).lower() if s is not None else ''


def column_index(header, token):
    """Position (1-based) of the column named `token` in the header row."""
    if token.isdigit():
        return int(token)
    want = norm(token)
    for i, h in enumerate(header):
        if norm(h) == want:
            return i + 1
    raise SystemExit('no column %r in %s' % (token, [h for h in header if h]))


def read_sheet(path, sheet):
    """Header row and data rows of one sheet, values only."""
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    if sheet not in wb.sheetnames:
        raise SystemExit('no sheet %r in %s (have %s)' % (sheet, path, wb.sheetnames))
    ws = wb[sheet]
    rows = list(ws.iter_rows(values_only=True))
    wb.close()
    if not rows:
        raise SystemExit('%s!%s is empty' % (path, sheet))
    return rows[0], rows[1:]


def parse_source(spec):
    """`ref=path:sheet:wncol:icol`, splitting from the right so paths may
    contain colons (as `F:\\...` does on Windows)."""
    if '=' not in spec:
        raise SystemExit('--source needs REF=FILE:SHEET:WNCOL:ICOL, got %r' % spec)
    ref, rest = spec.split('=', 1)
    parts = rest.rsplit(':', 3)
    if len(parts) != 4:
        raise SystemExit('--source needs REF=FILE:SHEET:WNCOL:ICOL, got %r' % spec)
    path, sheet, wncol, icol = parts
    return ref.strip(), path, sheet, wncol, icol


# ---------------------------------------------------------------------- main

def main(argv=None):
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p = argparse.ArgumentParser(
        description='Copy Sugar\'s printed intensities into the working line list, '
                    'matched by reconciled observed wavenumber.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--lines', default=os.path.join(here, 'Pr3_lines.xlsx'),
                   help='working line list to add the column to')
    p.add_argument('--sheet', default='Sheet1')
    p.add_argument('--col-wn', default='own',
                   help='observed wavenumber in the working list, cm^-1')
    p.add_argument('--col-ref', default='ref',
                   help='column naming whose classification each row carries; '
                        'not used to choose a source, only reported in the '
                        'audit file (give "" to leave it out)')
    p.add_argument('--col-out', default='Iorig',
                   help='column to write; created at the right end if absent, '
                        'overwritten in place if already there')
    p.add_argument('--source', action='append', default=[], metavar='SPEC',
                   help='REF=FILE:SHEET:WNCOL:ICOL, repeatable; WNCOL and ICOL '
                        'are header names or 1-based column numbers')
    p.add_argument('--tol', type=float, default=1e-4,
                   help='largest wavenumber difference, cm^-1, still called a match')
    p.add_argument('--out-csv', default=os.path.join(here, 'iorig_match.csv'),
                   help='audit of every row and what it was matched to')
    p.add_argument('--no-backup', action='store_true')
    p.add_argument('--dry-run', action='store_true',
                   help='report the matching, write nothing')
    a = p.parse_args(argv)

    if not a.source:
        raise SystemExit('at least one --source is required')

    # ---- the source tables, each sorted by wavenumber so a match is a bisect
    src = {}
    for spec in a.source:
        ref, path, sheet, wncol, icol = parse_source(spec)
        header, rows = read_sheet(path, sheet)
        jw = column_index(header, wncol) - 1
        ji = column_index(header, icol) - 1
        table = []
        for k, r in enumerate(rows):
            wn = r[jw] if jw < len(r) else None
            iv = r[ji] if ji < len(r) else None
            if isinstance(wn, (int, float)) and isinstance(iv, (int, float)):
                table.append((float(wn), float(iv), k + 2))
        table.sort()
        src[ref] = table
        print('%-12s %d usable lines from %s!%s (columns %d and %d)'
              % (ref, len(table), os.path.basename(path), sheet, jw + 1, ji + 1))
        if not table:
            raise SystemExit('no usable lines in %s' % path)

    keys = {ref: [t[0] for t in tab] for ref, tab in src.items()}

    # ---- the working list
    wb = openpyxl.load_workbook(a.lines)
    if a.sheet not in wb.sheetnames:
        raise SystemExit('no sheet %r in %s' % (a.sheet, a.lines))
    ws = wb[a.sheet]
    header = [ws.cell(1, c).value for c in range(1, ws.max_column + 1)]
    cw = column_index(header, a.col_wn)
    cr = column_index(header, a.col_ref) if a.col_ref else None
    try:
        co = column_index(header, a.col_out)
        print('column %r already there at %d; it is updated in place' % (a.col_out, co))
    except SystemExit:
        co = ws.max_column + 1
        print('column %r created at %d' % (a.col_out, co))

    matched = missing = 0
    worst = 0.0
    used = {}
    audit = []
    for row in range(2, ws.max_row + 1):
        wn = ws.cell(row, cw).value
        if not isinstance(wn, (int, float)):
            continue
        wn = float(wn)
        ref = ''
        if cr is not None:
            v = ws.cell(row, cr).value
            ref = str(v).strip() if v is not None else ''
        # the nearest source line, looked for in every source given
        best = None
        for name, tab in src.items():
            ks = keys[name]
            j = bisect.bisect_left(ks, wn)
            for jj in (j - 1, j, j + 1):
                if 0 <= jj < len(ks):
                    d = abs(ks[jj] - wn)
                    if best is None or d < best[0]:
                        best = (d, name, jj)
        if best is None or best[0] > a.tol:
            missing += 1
            audit.append((row, wn, ref, '', '', '', ''))
            continue
        d, name, jj = best
        wn_s, iv, srow = src[name][jj]
        matched += 1
        worst = max(worst, d)
        used[(name, jj)] = used.get((name, jj), 0) + 1
        audit.append((row, wn, ref, name, srow, '%.12g' % d, '%.12g' % iv))
        if not a.dry_run:
            ws.cell(row, co).value = iv

    shared = sum(1 for n in used.values() if n > 1)
    print('matched %d rows (worst wavenumber difference %.3g cm^-1)' % (matched, worst))
    if missing:
        print('%d rows found no source line within %g cm^-1; left as they were'
              % (missing, a.tol))
    if shared:
        print('%d source lines serve more than one row (alternative identifications)' % shared)
    for ref, tab in src.items():
        u = sum(1 for k in used if k[0] == ref)
        print('%-12s %d of its %d lines were used' % (ref, u, len(tab)))

    if a.dry_run:
        print('dry run: nothing written')
        return 0

    if not a.no_backup:
        bak = os.path.splitext(a.lines)[0] + '.bak_iorig.xlsx'
        shutil.copyfile(a.lines, bak)
        print('backup written to %s' % bak)

    ws.cell(1, co).value = a.col_out
    wb.save(a.lines)
    print('%s updated' % a.lines)

    with open(a.out_csv, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['row', 'wn', 'ref', 'source', 'source_row', 'dwn', a.col_out])
        w.writerows(audit)
    print('audit written to %s' % a.out_csv)
    return 0


if __name__ == '__main__':
    sys.exit(main())
