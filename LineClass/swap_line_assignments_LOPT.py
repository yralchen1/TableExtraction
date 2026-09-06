"""
swap_line_assignments_LOPT.py
=============================
Exchange two levels' observed lines in the LOPT transitions file, leaving
the level identifiers themselves untouched.

What this is for
----------------
Every energy level in this analysis is known by two different things.

  * Its **identifier** - a string such as ``059003.000483``.  It is only a
    name.  It appears in the LOPT input, in the classification tables, in
    the Cowan-code bookkeeping, and in lookup tables that live in other
    folders altogether, so it is expensive to change and must not be.

  * Its **energy** - a number in cm^-1.  Nothing writes it down directly:
    LOPT (the least-squares level optimiser) computes it from the observed
    lines that are assigned to that identifier.  Change which lines carry
    the identifier and you change the energy it comes out with.

``level_interchange.py`` looks for a specific mistake: two levels of the
same parity and the same J, close enough in energy that the calculation
cannot tell them apart, whose **theoretical identities were interchanged**.
The observed energies are right; what is wrong is which theoretical level -
which configuration, which term, and with them the whole set of calculated
transition probabilities - each measured energy was given.

The repair could be made by exchanging the two identifiers.  That is the
one thing that cannot be done cheaply here, because those identifiers are
scattered over files in several folders and in copies inside this
repository, and IDEN2's own level numbering is fixed by the *calculated*
energy and so cannot move at all.  The cheap repair is the mirror image of
it: **keep both identifiers where they are and exchange everything the
measurement gave them** - their observed lines, and hence their energies.

That is all this script does.  In the LOPT transitions file every line
assigned to ``id1`` becomes a line of ``id2`` and every line of ``id2``
becomes a line of ``id1``.  The next LOPT run then returns the two energies
exchanged, with the same uncertainties and the same residuals as before,
because nothing about the lines themselves has changed - only which of two
names each group of them is filed under.

``swap_line_assignments_IDEN.py`` is the companion that does the
corresponding thing to the IDEN2 working files.

What the file looks like
------------------------
The LOPT transitions file is fixed-column: LOPT is told in its parameter
file (``LOPT.par``) at which character each field begins and ends.  One
record is one observed line, and it carries

    the wavenumber, its uncertainty, the intensity, the identifier of the
    LOWER level, the identifier of the UPPER level, flags, a weight, and
    the unit of the first field.

Only the two identifier fields are touched here.  The columns they occupy
are read from the parameter file, so a re-arranged layout is followed
automatically; ``--par`` names it and ``--columns`` overrides it.

Every other character of every record - the wavenumber, the weight, the
flags, the spacing, the DOS line endings - is preserved byte for byte, so
a difference of the before and after files shows exactly the two columns
and nothing else.

What it checks before writing
-----------------------------
  * Both identifiers must actually occur in the file.  A typo in an
    identifier would otherwise produce a file that is silently unchanged.
  * No single record may join the two levels to each other.  Such a record
    would become a line from a level to itself, which is meaningless; it
    also cannot arise from a real interchange, since two levels of the same
    parity have no electric-dipole transition between them.
  * After the exchange no two records may name the same wavenumber and the
    same pair of levels, which would be one observed line counted twice.

Any of these stops the script before anything is written.

The fixed-levels file
---------------------
LOPT also reads a file of levels whose energies are held fixed.  A fixed
energy belongs to a *position in the spectrum*, not to a name, so when the
two names exchange positions the two fixed values must exchange with them.
If both identifiers appear there, their energies and uncertainties are
exchanged.  If only one appears the script stops: one of the two levels is
then pinned and the other is not, and what should happen is a decision
about the analysis, not something to guess.

The decision ledger, and what this script does not do
----------------------------------------------------
``line_decisions.csv`` records hand-made orders about individual lines -
accept this classification, reject that one - each keyed by a wavenumber
and the two level identifiers, with a written reason.  An exchange makes
every such order that names one of the two levels point at the wrong level.
This script prints every affected row and the identifier it would have to
become, and leaves the file alone;
``swap_line_assignments_pipeline.py`` is what rewrites it.

That is part of a larger division of labour.  The LOPT transitions file is
*generated* from ``line_classifications.csv`` by ``make_LOPT_input.py``, so
re-running the classification pipeline rebuilds it and undoes this
exchange.  What the pipeline would restore is the two **energies**: it gives
every level the energy of the published list, overridden by
``revised_level_energies.csv``, and generates that level's candidate
transitions there.  Recording the exchange so that it survives is therefore
a matter of writing the two exchanged energies into that overlay and
re-keying the ledger - both of which are
``swap_line_assignments_pipeline.py``'s work, not this script's.

Nothing in ``Icalc.xlsx`` changes.  It is keyed by the pair of level
identifiers, and an identifier keeps its configuration, its term and its
calculated transition probabilities through the exchange: it is the
measurement that moves, not the calculation.

Where the files are looked for
------------------------------
Each default file name is looked for first in the directory the command was
run from and then in the directory holding this script, so that working in
an iteration folder acts on that folder's copies while the files that exist
in only one place - the decision ledger above all - are still found.  The
report names the file it settled on.

Usage
-----
    python swap_line_assignments_LOPT.py 059003.000483 059003.000398
    python swap_line_assignments_LOPT.py id1 id2 --dry-run
    python swap_line_assignments_LOPT.py id1 id2 --lines iter22/LOPT_input_lines.txt
    python swap_line_assignments_LOPT.py id1 id2 --out swapped_lines.txt

Without ``--out`` the file is rewritten in place and the previous contents
are kept in ``<file>.bak`` (``--no-backup`` turns that off).  ``--dry-run``
reports and writes nothing.

``swap_line_assignments.py`` runs this script, its IDEN2 companion and
``swap_line_assignments_pipeline.py`` in one command.
"""

import argparse
import csv
import os
import sys

import swap_paths

# The names of the files, not their places: where each of them is looked for
# is decided when the command line is parsed, by swap_paths.working_path -
# the directory the command was run from if the file is there, otherwise the
# directory holding this script.  Standing in an iteration folder and asking
# for the transitions file therefore means that folder's transitions file.
NAME_LINES = 'LOPT_input_lines.txt'
NAME_PAR = 'LOPT.par'
NAME_FIXLEV = 'LOPT_fixlev.txt'
NAME_LEDGER = 'line_decisions.csv'

# The layout of the sample transitions file, used when no parameter file is
# available.  Columns are 1-based and inclusive, exactly as LOPT is given
# them in LOPT.par and as make_LOPT_input.py writes them.
DEF_COLUMNS = (43, 55, 59)


# ---------------------------------------------------------------------------
# Reading files without disturbing their line endings
# ---------------------------------------------------------------------------
def read_records(path):
    """Return (records, endings, trailer) for a fixed-column text file.

    ``records[i]`` is the i-th line with no end-of-line characters on it and
    ``endings[i]`` is the exact end-of-line that followed it ('\\r\\n', '\\n'
    or '' for a last line that had none).  Writing ``records[i] + endings[i]``
    back out reproduces the file byte for byte, whatever mixture of DOS and
    Unix endings it had.
    """
    with open(path, 'r', encoding='latin-1', newline='') as fh:
        text = fh.read()
    pieces = text.split('\n')
    trailing_newline = pieces[-1] == ''
    if trailing_newline:
        pieces = pieces[:-1]
    records, endings = [], []
    for k, piece in enumerate(pieces):
        last = (k == len(pieces) - 1)
        eol = '' if (last and not trailing_newline) else '\n'
        if piece.endswith('\r'):
            piece = piece[:-1]
            eol = '\r' + eol
        records.append(piece)
        endings.append(eol)
    return records, endings


def write_records(path, records, endings):
    with open(path, 'w', encoding='latin-1', newline='') as fh:
        for rec, eol in zip(records, endings):
            fh.write(rec + eol)


# ---------------------------------------------------------------------------
# Where the level identifiers sit
# ---------------------------------------------------------------------------
def read_columns(par_path):
    """The (first, last, first_upper) columns of the two identifier fields.

    LOPT's parameter file gives one setting per line as a value, a semicolon
    and a comment.  The three settings wanted here are found by their
    comments rather than by counting lines, so that an edited parameter file
    with extra or re-ordered settings still works.  Columns are 1-based and
    inclusive; the upper-level field has the same width as the lower one,
    which is why LOPT is given only its first column.
    """
    wanted = {
        'lower_first': ('1ST COLUMN', 'LOWER LEVEL'),
        'lower_last': ('LAST COLUMN', 'LOWER LEVEL'),
        'upper_first': ('1ST COLUMN', 'UPPER LEVEL'),
    }
    found = {}
    with open(par_path, 'r', encoding='latin-1') as fh:
        for raw in fh:
            if ';' not in raw:
                continue
            value, comment = raw.split(';', 1)
            comment = comment.upper()
            for key, needles in wanted.items():
                if key in found:
                    continue
                if all(n in comment for n in needles):
                    try:
                        found[key] = int(value.strip())
                    except ValueError:
                        pass
    missing = [k for k in wanted if k not in found]
    if missing:
        raise ValueError(
            '%s does not say where the level identifiers are: no setting '
            'for %s' % (par_path, ', '.join(sorted(missing))))
    return found['lower_first'], found['lower_last'], found['upper_first']


class Layout(object):
    """The character positions of the two identifier fields of a record."""

    def __init__(self, columns):
        first, last, upper_first = columns
        if last < first:
            raise ValueError('the lower-level field ends (column %d) before '
                             'it begins (column %d)' % (last, first))
        self.width = last - first + 1
        self.low = (first - 1, last)
        self.upp = (upper_first - 1, upper_first - 1 + self.width)

    def get(self, record, span):
        return record[span[0]:span[1]].strip()

    def put(self, record, span, value):
        """``record`` with ``value`` written left-justified into ``span``."""
        if len(value) > self.width:
            raise ValueError(
                'the identifier %r is %d characters and does not fit the '
                '%d-character field at columns %d-%d'
                % (value, len(value), self.width, span[0] + 1, span[1]))
        padded = record.ljust(span[1])
        return (padded[:span[0]] + value.ljust(self.width)
                + padded[span[1]:])


# ---------------------------------------------------------------------------
# The exchange itself
# ---------------------------------------------------------------------------
def swap_records(records, layout, id1, id2):
    """Exchange id1 and id2 in the two identifier fields of every record.

    Returns the new records and a report: how many times each identifier
    was found at each end of a transition, the wavenumbers involved, and
    any record that joins the two levels to each other.
    """
    out = []
    counts = {(id1, 'lower'): 0, (id1, 'upper'): 0,
              (id2, 'lower'): 0, (id2, 'upper'): 0}
    moved, self_lines, seen = [], [], {}
    for k, rec in enumerate(records):
        if not rec.strip():
            out.append(rec)
            continue
        low = layout.get(rec, layout.low)
        upp = layout.get(rec, layout.upp)
        if {low, upp} == {id1, id2}:
            self_lines.append((k, rec))
        new = rec
        if low in (id1, id2):
            counts[(low, 'lower')] += 1
            new = layout.put(new, layout.low, id2 if low == id1 else id1)
        if upp in (id1, id2):
            counts[(upp, 'upper')] += 1
            new = layout.put(new, layout.upp, id2 if upp == id1 else id1)
        if new != rec:
            moved.append((k, rec.split()[0] if rec.split() else '', low, upp))
        out.append(new)
        key = (rec.split()[0] if rec.split() else '',
               layout.get(new, layout.low), layout.get(new, layout.upp))
        seen.setdefault(key, []).append(k)
    duplicates = {k: v for k, v in seen.items() if len(v) > 1}
    return out, {'counts': counts, 'moved': moved,
                 'self_lines': self_lines, 'duplicates': duplicates}


# ---------------------------------------------------------------------------
# The fixed-levels file
# ---------------------------------------------------------------------------
def swap_fixed_levels(path, id1, id2):
    """Exchange the two levels' rows in the fixed-levels file.

    Returns (records, endings, note).  ``note`` says what was done.  A file
    that names neither level is left exactly as it is; one that names only
    one of them raises, because the two levels then differ in a way that
    the exchange cannot carry over on its own.
    """
    records, endings = read_records(path)
    hit = {}
    for k, rec in enumerate(records):
        fields = rec.split()
        if fields and fields[0] in (id1, id2):
            hit.setdefault(fields[0], []).append(k)
    if not hit:
        return records, endings, 'neither level is fixed; nothing to do'
    if len(hit) == 1:
        only = list(hit)[0]
        raise ValueError(
            '%s fixes the energy of %s but not of %s.  A fixed energy '
            'belongs to a position in the spectrum, so after the exchange it '
            'would be pinning the wrong level.  Decide by hand whether the '
            'other level should be fixed instead, or neither, and edit the '
            'file before running this.' % (path, only, id2 if only == id1 else id1))
    if any(len(v) > 1 for v in hit.values()):
        raise ValueError('%s names one of the two levels more than once'
                         % path)
    k1, k2 = hit[id1][0], hit[id2][0]
    # Keep each identifier where it is and give it the other's numbers.
    tail1 = records[k1][len(id1):]
    tail2 = records[k2][len(id2):]
    records[k1] = id1 + tail2
    records[k2] = id2 + tail1
    return records, endings, 'the two fixed energies were exchanged'


# ---------------------------------------------------------------------------
# The decision ledger
# ---------------------------------------------------------------------------
def ledger_rows(path, id1, id2):
    """Rows of the decision ledger that name either level, as printed text."""
    if not os.path.exists(path):
        return []
    rows = []
    with open(path, 'r', encoding='utf-8-sig', newline='') as fh:
        for n, row in enumerate(csv.DictReader(fh), start=2):
            low = (row.get('low_id') or '').strip()
            upp = (row.get('upp_id') or '').strip()
            if low not in (id1, id2) and upp not in (id1, id2):
                continue
            other = {id1: id2, id2: id1}
            rows.append((n, row.get('wn_obs', ''),
                         low, other.get(low, low),
                         upp, other.get(upp, upp),
                         (row.get('decision') or '').strip(),
                         (row.get('reason') or '').strip()))
    return rows


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def report(log, id1, id2, info, records, layout):
    counts = info['counts']
    log('')
    log('Lines found:')
    for lid in (id1, id2):
        low = counts[(lid, 'lower')]
        upp = counts[(lid, 'upper')]
        log('  %s   %3d as the lower level, %3d as the upper level, %3d in all'
            % (lid, low, upp, low + upp))
    log('  %d records rewritten.' % len(info['moved']))
    if info['moved']:
        log('')
        log('  %-14s  %-13s -> %-13s  %-13s -> %-13s'
            % ('wavenumber', 'lower was', 'becomes', 'upper was', 'becomes'))
        for k, wn, low, upp in info['moved']:
            other = {id1: id2, id2: id1}
            log('  %-14s  %-13s -> %-13s  %-13s -> %-13s'
                % (wn, low, other.get(low, low), upp, other.get(upp, upp)))


def check(info, id1, id2, path):
    """Raise on anything that makes the exchange meaningless."""
    counts = info['counts']
    for lid in (id1, id2):
        if counts[(lid, 'lower')] + counts[(lid, 'upper')] == 0:
            raise ValueError(
                '%s does not occur in %s.  Nothing would be exchanged; check '
                'the identifier.' % (lid, path))
    if info['self_lines']:
        k, rec = info['self_lines'][0]
        raise ValueError(
            'line %d of %s joins %s to %s, so exchanging them would make a '
            'transition of a level with itself:\n    %s'
            % (k + 1, path, id1, id2, rec.rstrip()))
    if info['duplicates']:
        key, where = sorted(info['duplicates'].items())[0]
        raise ValueError(
            'after the exchange the wavenumber %s would be assigned to '
            '%s - %s by %d different records (lines %s of %s), i.e. the same '
            'observed line counted more than once'
            % (key[0], key[1], key[2], len(where),
               ', '.join(str(w + 1) for w in where), path))


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def parse_args(argv):
    p = argparse.ArgumentParser(
        description='Exchange two levels\' observed lines in the LOPT '
                    'transitions file, keeping both level identifiers where '
                    'they are.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('id1', help='identifier of the first level')
    p.add_argument('id2', help='identifier of the second level')
    p.add_argument('--lines', default=swap_paths.working_path(NAME_LINES),
                   help='the LOPT transitions file to rewrite')
    p.add_argument('--par', default=swap_paths.working_path(NAME_PAR),
                   help='the LOPT parameter file, read for the column '
                        'positions of the two identifier fields')
    p.add_argument('--columns', nargs=3, type=int, metavar=('LOW1', 'LOW2', 'UPP1'),
                   help='override the parameter file: first and last column '
                        'of the lower-level field and first column of the '
                        'upper-level field, 1-based and inclusive')
    p.add_argument('--fixlev', default=swap_paths.working_path(NAME_FIXLEV),
                   help='the LOPT fixed-levels file')
    p.add_argument('--no-fixlev', action='store_true',
                   help='leave the fixed-levels file alone')
    p.add_argument('--ledger', default=swap_paths.working_path(NAME_LEDGER),
                   help='the decision ledger, reported on but never rewritten')
    p.add_argument('--out',
                   help='write the rewritten transitions file here instead '
                        'of over the original')
    p.add_argument('--no-backup', action='store_true',
                   help='do not keep the previous contents in <file>.bak')
    p.add_argument('--dry-run', action='store_true',
                   help='report what would change and write nothing')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    log = print
    id1, id2 = args.id1.strip(), args.id2.strip()
    if id1 == id2:
        log('The two identifiers are the same; there is nothing to exchange.')
        return 2

    if args.columns:
        columns = tuple(args.columns)
        source = 'the --columns option'
    elif os.path.exists(args.par):
        columns = read_columns(args.par)
        source = args.par
    else:
        columns = DEF_COLUMNS
        source = 'the built-in default layout (%s was not found)' % args.par
    layout = Layout(columns)

    log('swap_line_assignments_LOPT.py')
    log('  transitions file : %s'
        % swap_paths.describe(NAME_LINES, args.lines))
    log('  level fields     : lower at columns %d-%d, upper at columns %d-%d'
        % (layout.low[0] + 1, layout.low[1],
           layout.upp[0] + 1, layout.upp[1]))
    log('                     (from %s)' % source)
    log('  exchanging       : %s  <->  %s' % (id1, id2))

    records, endings = read_records(args.lines)
    new_records, info = swap_records(records, layout, id1, id2)
    report(log, id1, id2, info, records, layout)
    check(info, id1, id2, args.lines)

    fix_records = fix_endings = None
    fix_note = 'not touched (--no-fixlev)'
    if not args.no_fixlev and os.path.exists(args.fixlev):
        fix_records, fix_endings, fix_note = swap_fixed_levels(
            args.fixlev, id1, id2)
        if 'nothing to do' in fix_note:
            fix_records = None
    elif not args.no_fixlev:
        fix_note = 'not found: %s' % args.fixlev
    log('')
    log('Fixed levels: %s' % fix_note)

    rows = ledger_rows(args.ledger, id1, id2)
    log('')
    if rows:
        log('The decision ledger %s has %d order(s) naming one of these two '
            'levels.  This script does not touch them; '
            'swap_line_assignments_pipeline.py re-keys them, together with '
            'the level energies the pipeline reads:'
            % (args.ledger, len(rows)))
        for n, wn, low, low_new, upp, upp_new, decision, reason in rows:
            log('  row %d: %s  %s -> %s   %s -> %s   %s   %s'
                % (n, wn, low, low_new, upp, upp_new, decision, reason))
    else:
        log('The decision ledger names neither level; no order is affected.')

    if args.dry_run:
        log('')
        log('--dry-run: nothing was written.')
        return 0

    target = args.out or args.lines
    if not args.out and not args.no_backup:
        backup = args.lines + '.bak'
        write_records(backup, records, endings)
        log('')
        log('Previous contents kept in %s' % backup)
    write_records(target, new_records, endings)
    log('Written: %s' % target)

    if fix_records is not None:
        if not args.no_backup:
            saved, saved_eol = read_records(args.fixlev)
            write_records(args.fixlev + '.bak', saved, saved_eol)
            log('Previous contents kept in %s' % (args.fixlev + '.bak'))
        write_records(args.fixlev, fix_records, fix_endings)
        log('Written: %s' % args.fixlev)

    log('')
    log('Run LOPT again to get the two exchanged energies.  This file is '
        'still only a working copy: make_LOPT_input.py rebuilds it from '
        'line_classifications.csv, so until swap_line_assignments_pipeline.py '
        'has written the exchanged energies into the level overrides and '
        're-keyed the decision ledger, the next run of the pipeline puts both '
        'levels back where they were.')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except ValueError as exc:
        sys.stderr.write('\nswap_line_assignments_LOPT.py: %s\n' % exc)
        sys.exit(1)
