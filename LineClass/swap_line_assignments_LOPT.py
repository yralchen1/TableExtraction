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

Every other character of every record - the wavenumber, the flags, the
spacing, the DOS line endings - is preserved byte for byte, with the one
exception described next, so a difference of the before and after files
shows the two identifier columns, the weights of the disturbed blends, and
nothing else.

The shares of a blend
---------------------
Two observed transitions can fall on one measured line and be impossible to
separate: a **blend**.  ``make_LOPT_input.py`` gives such a line one record
per component and divides its weight, 1, between them in proportion to their
**calculated** intensities - the ``Icalc`` column of ``Icalc.xlsx``, the
intensity Cowan's codes predict for that pair of levels - so that LOPT
attributes to each level the share of the measured position that the
calculation says belongs to it.

An exchange invalidates those shares.  A component that named one of the two
levels now names the other, and the calculation predicts a different
intensity for the new pair; the share worked out for the old pair means
nothing.  So this script recomputes them.  For every measured wavenumber
where at least one component moved it reads the calculated intensity of each
component's new pair and divides the line's weight in proportion to them
again, exactly as ``make_LOPT_input.py`` would, and writes the new shares
into the weight column.  ``--no-reweight`` leaves the old shares in place.

It refuses to guess.  A group is recomputed only when

  * the calculated file gives an intensity for every component of it, before
    and after the exchange, and
  * recomputing the group's **old** shares reproduces the weights the file
    actually carries, to within ``--reweight-tol`` (0.002 by default).

A share too small to be written as a four-decimal number is written in
exponential form instead - ``25e-06`` rather than ``0.0000`` - which LOPT
reads just as well.  Nothing is ever written as an exact zero: LOPT stops
with an error message when a record that is not flagged "P" (predicted) or
"M" (masked) carries the weight 0, so a group in which the calculation gives
a component no intensity at all is left alone instead.

Faint shares are reported as a warning.  A component that comes out with
less than a tenth of its measured line is one the calculation does not
really place on that transition; fitting it with a negligible weight is
worth less than rejecting the assignment, and the report says which
components those are so that the decision can be made.

The second test is what makes the first safe.  It checks, on the file
itself, that this script is reading the same intensities and applying the
same rule as the run that wrote the weights.  When it fails - a transitions
file left over from an older classification, an intensity model refitted
since, a pair whose intensity had to be imputed rather than read - the group
is left exactly as it is and the report says which group and why.  Such a
group is then stale until ``classify_lines.py`` and ``make_LOPT_input.py``
rebuild the file.

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

Nothing in ``Icalc.xlsx`` changes; it is only read.  It is keyed by the
pair of level identifiers, and an identifier keeps its configuration, its
term and its calculated transition probabilities through the exchange: it is
the measurement that moves, not the calculation.  That is precisely why the
shares of a blend have to be recomputed - the same measured line is now
being divided between two different predictions.

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
    python swap_line_assignments_LOPT.py id1 id2 --no-reweight

Without ``--out`` the file is rewritten in place and the previous contents
are kept in ``<file>.bak`` (``--no-backup`` turns that off).  ``--dry-run``
reports and writes nothing.

``swap_line_assignments.py`` runs this script, its IDEN2 companion and
``swap_line_assignments_pipeline.py`` in one command.
"""

import argparse
import csv
import io
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

# The weight field of the same sample layout, 1-based and inclusive.  It is
# read from the parameter file like the others when that file names it; only
# the recomputed blend shares are written into it.
DEF_WEIGHT_COLUMNS = (82, 87)

# How far a recomputed old share may sit from the weight the file carries
# before the group is left alone; see "The shares of a blend" above.
DEF_REWEIGHT_TOL = 0.002
FAINT_SHARE = 0.10


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
        'weight_first': ('1ST COLUMN', 'WEIGHT'),
        'weight_last': ('LAST COLUMN', 'WEIGHT'),
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
    # The weight is wanted only for the recomputed blend shares; a parameter
    # file that does not name its columns falls back to the sample layout.
    missing = [k for k in wanted if k not in found and 'weight' not in k]
    if missing:
        raise ValueError(
            '%s does not say where the level identifiers are: no setting '
            'for %s' % (par_path, ', '.join(sorted(missing))))
    return (found['lower_first'], found['lower_last'], found['upper_first'],
            found.get('weight_first', DEF_WEIGHT_COLUMNS[0]),
            found.get('weight_last', DEF_WEIGHT_COLUMNS[1]))


class Layout(object):
    """The character positions of the two identifier fields of a record."""

    def __init__(self, columns):
        columns = tuple(columns)
        if len(columns) == 3:
            columns += DEF_WEIGHT_COLUMNS
        first, last, upper_first, weight_first, weight_last = columns
        if last < first:
            raise ValueError('the lower-level field ends (column %d) before '
                             'it begins (column %d)' % (last, first))
        self.width = last - first + 1
        self.low = (first - 1, last)
        self.upp = (upper_first - 1, upper_first - 1 + self.width)
        self.wt = (weight_first - 1, weight_last)

    def get(self, record, span):
        return record[span[0]:span[1]].strip()

    def put(self, record, span, value):
        """``record`` with ``value`` written left-justified into ``span``."""
        width = span[1] - span[0]
        if len(value) > width:
            raise ValueError(
                'the value %r is %d characters and does not fit the '
                '%d-character field at columns %d-%d'
                % (value, len(value), width, span[0] + 1, span[1]))
        padded = record.ljust(span[1])
        return (padded[:span[0]] + value.ljust(width) + padded[span[1]:])


# ---------------------------------------------------------------------------
# The exchange itself
# ---------------------------------------------------------------------------
def record_weight(record):
    """The weight LOPT gives this record, or None if the record has none.

    The weight is the last number of the record - the field after the flags
    and before the unit.  A record flagged "P" carries the weight 0: LOPT
    prints it and fits nothing to it, which is how ``make_LOPT_input.py``
    writes a classification the pipeline did not accept (``accepted = 0``).
    A weight strictly between 0 and 1 means the observed line is a blend and
    this record is one component of it, holding the share of the measured
    intensity that the calculated intensities gave it.
    """
    for token in reversed(record.split()):
        try:
            return float(token)
        except ValueError:
            continue
    return None


def swap_records(records, layout, id1, id2):
    """Exchange id1 and id2 in the two identifier fields of every record.

    Returns the new records and a report: how many times each identifier
    was found at each end of a transition, the wavenumbers involved, and
    any record that joins the two levels to each other.
    """
    out = []
    counts = {(id1, 'lower'): 0, (id1, 'upper'): 0,
              (id2, 'lower'): 0, (id2, 'upper'): 0}
    fitted = dict.fromkeys(counts, 0)
    moved, self_lines, seen = [], [], {}
    weight_of, by_wavenumber, about = {}, {}, {}
    for k, rec in enumerate(records):
        if not rec.strip():
            out.append(rec)
            continue
        low = layout.get(rec, layout.low)
        upp = layout.get(rec, layout.upp)
        if {low, upp} == {id1, id2}:
            self_lines.append((k, rec))
        weight = record_weight(rec)
        weight_of[k] = weight
        wn_text = rec.split()[0] if rec.split() else ''
        if weight:
            by_wavenumber.setdefault(wn_text, []).append(k)
        new = rec
        if low in (id1, id2):
            counts[(low, 'lower')] += 1
            if weight:
                fitted[(low, 'lower')] += 1
            new = layout.put(new, layout.low, id2 if low == id1 else id1)
        if upp in (id1, id2):
            counts[(upp, 'upper')] += 1
            if weight:
                fitted[(upp, 'upper')] += 1
            new = layout.put(new, layout.upp, id2 if upp == id1 else id1)
        if new != rec:
            moved.append((k, wn_text, low, upp))
        about[k] = {'wn': wn_text, 'weight': weight, 'moved': new != rec,
                    'old': (low, upp),
                    'new': (layout.get(new, layout.low),
                            layout.get(new, layout.upp))}
        out.append(new)
        key = (rec.split()[0] if rec.split() else '',
               layout.get(new, layout.low), layout.get(new, layout.upp))
        seen.setdefault(key, []).append(k)
    duplicates = {k: v for k, v in seen.items() if len(v) > 1}
    # A moved record that shares its wavenumber with another weighted record
    # is one component of a blend, and its weight is a share worked out from
    # the calculated intensities of the pair it named BEFORE the exchange.
    blends = []
    for k, wn, low, upp in moved:
        if not weight_of.get(k):
            continue
        company = [j for j in by_wavenumber.get(wn, []) if j != k]
        if company:
            blends.append((wn, low, upp, weight_of[k], len(company) + 1))
    return out, {'counts': counts, 'fitted': fitted, 'moved': moved,
                 'self_lines': self_lines, 'duplicates': duplicates,
                 'blends': blends, 'about': about,
                 'by_wavenumber': by_wavenumber}


# ---------------------------------------------------------------------------
# The calculated intensities, and the shares of a blend
# ---------------------------------------------------------------------------
def read_calc_intens(log=None):
    """{(lower id, upper id): calculated intensity} for every predicted pair.

    The calculated intensity is what divides a blend between its components,
    so it has to be the same number ``classify_lines.py`` used.  It is
    therefore read by that script's own reader, from the files its
    configuration names: ``Icalc.xlsx`` and, for levels added since that
    workbook was made, the supplementary ``icalc_new.xlsx``, whose intensity
    is recomputed from gA there.  Pairs the reader drops - below the printing
    cutoff of Cowan's codes - are absent here too, and a group containing one
    of them is left alone rather than guessed at.

    Reading it prints a couple of lines of its own; they are passed on to
    ``log`` indented, so that the report says where the numbers came from.
    """
    import contextlib
    import classify_lines as cl

    class AnyLevel(object):
        """Stands in for the level list: every identifier is wanted here."""

        def __contains__(self, item):
            return True

    said = io.StringIO()
    with contextlib.redirect_stdout(said):
        index = cl.read_transitions(AnyLevel())
    if log is not None:
        for line in said.getvalue().splitlines():
            if line.strip() and not line.startswith('Step '):
                log('  %s' % line.strip())
    return {pair: data['calc_intensity'] for pair, data in index.items()
            if data['calc_intensity'] is not None}


def shares(intensities):
    """The intensities as fractions of their sum, or None if they say nothing.

    This is ``make_LOPT_input.blend_weights`` for a group whose calculated
    intensities are all known: each component takes the share of the measured
    line that its predicted intensity is of the total.
    """
    total = sum(intensities)
    if total <= 0:
        return None
    return [v / total for v in intensities]


def weight_text(value, width):
    """``value`` as at most ``width`` characters of a number LOPT can read.

    A weight field is narrow - six characters in the sample file - and a
    share can be very small.  The plain four-decimal form is used whenever it
    says something; when it would collapse to ``0.0000`` the value is written
    in exponential form instead, with as many figures as the field allows and
    the leading zeros of the exponent dropped, so that 2.5e-05 is written
    ``2.5e-5`` rather than being lost.  LOPT reads either form.
    """
    plain = '%.*f' % (4, value)
    if len(plain) <= width and (float(plain) != 0 or value == 0):
        return plain
    for digits in range(width, -1, -1):
        mantissa, _, exponent = ('%.*e' % (digits, value)).partition('e')
        sign, figures = exponent[0], exponent[1:].lstrip('0') or '0'
        text = '%se%s%s' % (mantissa, '' if sign == '+' else sign, figures)
        if len(text) <= width:
            return text
    raise ValueError('%r does not fit a %d-character weight field'
                     % (value, width))


def reweight_blends(records, layout, info, calc, tol=DEF_REWEIGHT_TOL):
    """Recompute the share of every blend the exchange disturbed.

    ``records`` are the records as the exchange left them and ``calc`` is
    ``read_calc_intens()``.  Returns ``(records, report)``; the records are
    returned changed only where a group passed both tests described in "The
    shares of a blend" at the top of this file.  ``report`` holds

        ``changed``  one entry per rewritten record: wavenumber, the pair it
                     names now, whether that pair is new, the weight it had,
                     the weight it has been given and how that weight was
                     written into the field;
        ``kept``     one entry per group left alone: wavenumber and why;
        ``faint``    the entries of ``changed`` whose new share is below
                     ``FAINT_SHARE``, which are reported as a warning.
    """
    out = list(records)
    changed, kept, faint = [], [], []
    about = info['about']
    for wn, ks in sorted(info['by_wavenumber'].items(),
                         key=lambda t: -_number(t[0])):
        if len(ks) < 2 or not any(about[k]['moved'] for k in ks):
            continue                      # not a blend, or nothing moved in it
        old = [calc.get(about[k]['old']) for k in ks]
        new = [calc.get(about[k]['new']) for k in ks]
        missing = [about[k]['new'] for k, v in zip(ks, new) if v is None]
        missing += [about[k]['old'] for k, v in zip(ks, old) if v is None]
        if missing:
            kept.append((wn, 'no calculated intensity for %s'
                         % ', '.join('%s - %s' % p for p in missing)))
            continue
        was, now = shares(old), shares(new)
        if was is None or now is None:
            kept.append((wn, 'the calculated intensities of the group add up '
                             'to nothing'))
            continue
        off = max(abs(a - about[k]['weight']) for a, k in zip(was, ks))
        if off > tol:
            kept.append((wn, 'the weights in the file are not the ones these '
                             'calculated intensities give (largest difference '
                             '%.4f, tolerance %.4f): the file was written by '
                             'another run' % (off, tol)))
            continue
        if min(now) <= 0:
            kept.append((wn, 'the calculation gives a component no intensity '
                             'at all, and LOPT stops on a zero weight unless '
                             'the record is flagged "P" or "M"'))
            continue
        width = layout.wt[1] - layout.wt[0]
        for k, w in zip(ks, now):
            text = weight_text(w, width)
            out[k] = layout.put(out[k], layout.wt, text)
            entry = (wn, about[k]['new'][0], about[k]['new'][1],
                     about[k]['moved'], about[k]['weight'], w, text)
            changed.append(entry)
            if w < FAINT_SHARE:
                faint.append(entry)
    return out, {'changed': changed, 'kept': kept, 'faint': faint, 'tol': tol}


def _number(text):
    try:
        return float(text)
    except ValueError:
        return 0.0


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
    counts, fitted = info['counts'], info['fitted']
    log('')
    log('Lines found  (every record of the file, whether LOPT fits it or not):')
    for lid in (id1, id2):
        low = counts[(lid, 'lower')]
        upp = counts[(lid, 'upper')]
        keep = fitted[(lid, 'lower')] + fitted[(lid, 'upper')]
        log('  %s   %3d as the lower level, %3d as the upper level, %3d in all'
            % (lid, low, upp, low + upp))
        log('  %s   of those, %d carry a weight and are fitted; %d are flagged '
            '"P"' % (' ' * len(lid), keep, low + upp - keep))
        log('  %s   (accepted = 0 in the classification: LOPT prints them and '
            'fits nothing to them)' % (' ' * len(lid)))
    log('  %d records rewritten.  Every one of them moves: an exchange changes '
        'only which' % len(info['moved']))
    log('  of the two names each record carries, so no assignment can be lost '
        'here.')
    if info['blends']:
        log('')
        log('  %d of the moved records are components of a blend: one measured '
            'line divided' % len(info['blends']))
        log('  between the transitions that fall on it, in proportion to their '
            'calculated')
        log('  intensities.  The exchange gives the moved component a '
            'different pair of')
        log('  levels, for which the calculation predicts a different '
            'intensity, so the')
        log('  share it was given no longer applies.')
        rew = info.get('reweight')
        if rew is None:
            log('  --no-reweight: the old shares are carried over unchanged '
                'and stand until')
            log('  make_LOPT_input.py rebuilds the file from a fresh '
                'classify_lines.py run.')
        elif rew.get('trouble'):
            log('  The calculated intensities could not be read, so the old '
                'shares are carried')
            log('  over unchanged: %s' % rew['trouble'])
        if rew and rew.get('changed'):
            log('')
            log('  Shares recomputed from the calculated intensities of the '
                'new pairs:')
            log('    %-14s  %-13s  %-13s  %8s  %8s  %s'
                % ('wavenumber', 'lower', 'upper', 'was', 'becomes', 'pair'))
            for wn, low, upp, moved, was, _now, text in rew['changed']:
                log('    %-14s  %-13s  %-13s  %8.4f  %8s  %s'
                    % (wn, low, upp, was, text,
                       'exchanged' if moved else 'unchanged'))
        if rew and rew.get('faint'):
            log('')
            log('  WARNING: %d of them now take less than %g per cent of the '
                'measured line.' % (len(rew['faint']), FAINT_SHARE * 100))
            log('  A share that small says the calculation does not really '
                'place the line on')
            log('  that transition.  Rejecting the assignment is worth more '
                'than fitting it')
            log('  with a negligible weight; these are the ones to look at:')
            log('    %-14s  %-13s  %-13s  %8s  %s'
                % ('wavenumber', 'lower', 'upper', 'share', 'pair'))
            for wn, low, upp, moved, _was, _now, text in rew['faint']:
                log('    %-14s  %-13s  %-13s  %8s  %s'
                    % (wn, low, upp, text,
                       'exchanged' if moved else 'unchanged'))
        if rew and rew.get('kept'):
            log('')
            log('  Left as they are, and stale until make_LOPT_input.py '
                'rebuilds the file:')
            for wn, why in rew['kept']:
                log('    %-14s  %s' % (wn, why))
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
    p.add_argument('--no-reweight', action='store_true',
                   help='do not recompute the shares of the blends the '
                        'exchange disturbs; carry the old shares over, as '
                        'every other field is carried over')
    p.add_argument('--reweight-tol', type=float, default=DEF_REWEIGHT_TOL,
                   metavar='W',
                   help='how far a recomputed old share may sit from the '
                        'weight the file carries before the group is left '
                        'alone as written by another run')
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
    if info['blends'] and not args.no_reweight:
        log('')
        log('Calculated intensities, for the shares of the %d disturbed blend '
            'component(s):' % len(info['blends']))
        try:
            calc = read_calc_intens(log)
        except Exception as exc:                  # a missing or unreadable file
            info['reweight'] = {'changed': [], 'kept': [],
                                'trouble': '%s: %s'
                                           % (type(exc).__name__, exc)}
            log('  could not be read: %s' % exc)
        else:
            log('  %d predicted pairs read.' % len(calc))
            new_records, rew = reweight_blends(new_records, layout, info,
                                               calc, args.reweight_tol)
            info['reweight'] = rew
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
