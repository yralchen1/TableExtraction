"""
swap_line_assignments_IDEN.py
=============================
Exchange two levels' observed energies and observed lines in the IDEN2
working files, leaving IDEN2's level numbering untouched.

What this is for
----------------
This is the companion of ``swap_line_assignments_LOPT.py``.  Both repair
the same mistake in the same way; they differ only in which files they
rewrite.

The mistake is the one ``level_interchange.py`` looks for: two levels of
the same parity and the same J, close enough in energy that the calculation
cannot tell them apart, whose **theoretical identities were interchanged**.
Both measured energies are right; what is wrong is which calculated level -
which configuration, which term, and with them the whole set of calculated
transition probabilities - each measured energy was attached to.

The obvious repair, exchanging the two levels' names, is impossible here.
IDEN2 knows a level by its **row number in ``enlev.dat``**, and that file is
sorted by *calculated* energy.  A row number is therefore a statement about
the calculation, and it cannot move without re-sorting the file and
invalidating every other file that refers to a row by number.  So the row
numbers stay exactly where they are and everything the *measurement*
contributed is exchanged instead:

  * the observed energy of each of the two levels, and its uncertainty;
  * every observed line assigned to them.

Afterwards row ``n1`` carries what row ``n2`` used to carry and the other
way round, and IDEN2 shows each calculated level against the observed lines
that now belong to it.

Which files are rewritten
-------------------------
``enlev.dat`` - one row per calculated level, fixed-column::

    columns  1-4    the row number, which is what IDEN2 calls the level
             5-16   the calculated energy      (never touched)
            17-26   the uncertainty of the observed energy
            27-38   the observed energy
            39-40   " *" when the level has been observed at all
            41-52   observed minus calculated
            53-57   J
            58-     the label: configuration, a separator, the term

The uncertainty, the observed energy and the asterisk belong to the
measurement and are exchanged between the two rows.  The calculated energy
and the label belong to the calculation and stay.  Observed minus calculated
is then recomputed for each row from its own calculated energy, so it comes
out as the O-C the exchanged assignment implies - which is the whole point
of the exercise, and worth looking at in the report.

``trans.dat`` - the predicted transitions.  It is written in blocks: a
header row beginning with ``$`` names a level, and the rows under it are
that level's predicted transitions to levels **lower down the file**, one
row each::

    columns  1      always "+"
             2-5    the row number of the level at the other end
             6-10   the predicted intensity of the transition
            11-22   the observed energy of the level at the other end
            23-24   " *" when that level has been observed
            25-37   the predicted wavenumber, |E(this) - E(other)|
            38-43   the intensity of the observed line assigned to it, 0 if none
            44-55   the wavenumber of that observed line
            56-66   observed minus predicted
            67-72   the row of that line in dlv.dat, 0 if none

Each predicted transition appears exactly once, in the block of whichever
of its two levels stands higher in the file.  The last four fields together
are the **assignment**: the observed line that has been identified with the
predicted transition.  Those four fields are what moves.

For every level ``p`` the assignment of the transition ``(n1, p)`` and the
assignment of the transition ``(n2, p)`` are exchanged, and likewise for
every block ``b`` that contains a row to ``n1`` and a row to ``n2``.  Then
every row that mentions either level has its observed-energy field, its
asterisk and its predicted wavenumber recomputed from the exchanged
energies, and its observed-minus-predicted recomputed from that.

``dlv.dat`` (the observed line list) and ``numset.dat`` carry no level
information at all and are not touched.

A self-check that comes for free
--------------------------------
An observed line assigned to the transition ``(n1, p)`` sat at some
observed-minus-predicted value.  After the exchange that line is assigned
to ``(n2, p)``, and ``n2`` now holds the energy ``n1`` used to hold, so the
predicted wavenumber of its new row is the same number as the predicted
wavenumber of its old one.  **Every transferred line must therefore keep its
observed-minus-predicted value exactly.**  The script verifies this for
every line it moves and stops if any of them changes: that would mean the
files were not in the state this rewrite assumes.

Lines that cannot be transferred
--------------------------------
A predicted transition is in ``trans.dat`` only if the calculation gave it
an intensity above the printing threshold.  Two levels that are being
exchanged do not have the same set of predicted transitions, so a level
``p`` may be connected to ``n1`` in the file and not to ``n2``.  If an
observed line is assigned to such a transition it has nowhere to go: after
the exchange the level it belongs to has no row for it.

Those lines are reported one by one - the predicted transition they were on,
their wavenumber and their row in ``dlv.dat`` - and their assignment is
removed, because leaving it in place would file the line under the wrong
level, which is the very mistake being repaired.  They are not lost: they
are still in ``dlv.dat`` and can be re-identified in IDEN2 by hand.  Read
this part of the report; a level that loses several strong lines this way is
a sign that the exchange needs a second look.

Naming the two levels
---------------------
Each of the two mandatory arguments is **either** the experimental level
identifier used everywhere else in this analysis, such as
``059003.000483``, **or** the level's row number in ``enlev.dat``, which is
what IDEN2 itself calls the level.  A bare run of digits is read as a row
number and anything else as an identifier, so the two forms need no flag to
tell them apart and may be mixed.

``--index`` is a shortcut for the *lookup* and for nothing else.  What the
script does to the two rows does not depend on how they were named: their
observed energies, their uncertainties and their observed-line assignments
are exchanged in ``enlev.dat`` and ``trans.dat`` in every case.  Without
``--index`` the two identifiers are first translated into row numbers by one
of the routes below, and the exchange then proceeds exactly as it would have
done had the numbers been typed.

A row number needs no translation and is the shortest way to name a level
when IDEN2 is already open in front of you.  Identifiers are accepted
because this script is the companion of
``swap_line_assignments_LOPT.py``, which can be given nothing else - the
LOPT files know levels only by identifier - and because
``swap_line_assignments.py`` drives both from one pair of arguments.  An
identifier has to be translated into a row number, and there are three ways,
in order of precedence:

  ``--index N1 N2``   give the two row numbers directly.
  ``--map FILE``      a two-column table of identifier and row number, i.e.
                      one of the lookup tables kept alongside the analysis.
                      Its columns may be named or not; a header naming them
                      ``level_id`` and something containing ``iden``, ``num``
                      or ``index`` is used if present.
  by energy           the default.  The identifiers are looked up in a table
                      of experimental energies (``--levels``, by default
                      LOPT's output level list) and each is matched to the
                      row of ``enlev.dat`` whose observed energy agrees to
                      within ``--tol`` (0.5 cm^-1).  The match must be
                      unique.

The energy route reads the state **before** the exchange, so run this script
before re-running LOPT, or give it the level table that was current when the
IDEN2 files were made.  Naming the rows by number avoids the question
entirely.

Two guards catch the commonest mistakes.  The script stops if the two rows
do not have the same J, since levels of different J cannot have had their
identities interchanged.  And when the rows were named with ``--index`` or
``--map`` it also stops if they already hold each other's energies - that is
what a second run looks like, and a second run would quietly undo the first.
The energy route cannot see this, because it finds whichever row holds the
identifier's energy; run it once, on the files as IDEN2 last left them.

Usage
-----
    python swap_line_assignments_IDEN.py 059003.000483 059003.000398
    python swap_line_assignments_IDEN.py 721 723
    python swap_line_assignments_IDEN.py id1 id2 --dry-run
    python swap_line_assignments_IDEN.py id1 id2 --index 721 723
    python swap_line_assignments_IDEN.py id1 id2 --iden2-dir IDEN2_snr
    python swap_line_assignments_IDEN.py id1 id2 --out-dir IDEN2_swapped

Without ``--out-dir`` the files are rewritten in place and their previous
contents are kept in ``<file>.bak`` (``--no-backup`` turns that off).
``--dry-run`` reports and writes nothing.

``--iden2-dir`` and ``--levels`` default to the directory the command was
run from if the file is there and to the directory holding this script
otherwise, so that working in an iteration folder acts on that folder's
copies.

This script rewrites only IDEN2's own working files.  What makes an
exchange survive a re-run of the classification pipeline is
``swap_line_assignments_pipeline.py``, which writes the two exchanged
energies into ``revised_level_energies.csv`` and re-keys
``line_decisions.csv``; ``swap_line_assignments.py`` runs all three scripts
in one command.
"""

import argparse
import csv
import os
import sys

import swap_paths

# The names of the two things the script needs to find, not their places:
# where each is looked for is decided when the command line is parsed, by
# swap_paths.working_path - the directory the command was run from if it is
# there, otherwise the directory holding this script.
NAME_IDEN2 = 'IDEN2'
NAME_LEVELS = 'LOPT_output_levels.txt'
DEF_TOL = 0.5

# enlev.dat, as character positions (0-based, end-exclusive).
EN_INDEX = (0, 4)
EN_ECALC = (4, 16)
EN_UNC = (16, 26)
EN_EOBS = (26, 38)
EN_STAR = (38, 40)
EN_OMC = (40, 52)
EN_J = (52, 57)
EN_LABEL = 57

# trans.dat header rows.
TH_INDEX = (1, 5)
TH_EOBS = (28, 41)
TH_STAR = (41, 43)

# trans.dat transition rows.
TR_PARTNER = (1, 5)
TR_EPART = (10, 22)
TR_STAR = (22, 24)
TR_WN = (24, 37)
TR_OBS = 37              # the assignment: everything from here to the end
TR_OBS_WN = (6, 18)      # inside the assignment
TR_OBS_OMC = (18, 29)    # inside the assignment
TR_OBS_ROW = (29, 35)    # inside the assignment
BLANK_OBS = '     0       0.000      0.000     0'

STAR_ON = ' *'
STAR_OFF = '  '


# ---------------------------------------------------------------------------
# Reading and writing without disturbing the line endings
# ---------------------------------------------------------------------------
def read_records(path):
    """Return (records, endings): the lines without their end-of-line
    characters, and the exact end-of-line that followed each, so that the
    file can be written back byte for byte."""
    with open(path, 'r', encoding='latin-1', newline='') as fh:
        text = fh.read()
    pieces = text.split('\n')
    trailing = pieces[-1] == ''
    if trailing:
        pieces = pieces[:-1]
    records, endings = [], []
    for k, piece in enumerate(pieces):
        eol = '' if (k == len(pieces) - 1 and not trailing) else '\n'
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


def field(record, span):
    return record[span[0]:span[1]]


def put(record, span, text):
    """``record`` with ``text`` written into ``span``, which it must fill."""
    width = span[1] - span[0]
    if len(text) != width:
        raise ValueError('%r is %d characters, not the %d of the field at '
                         '%d-%d' % (text, len(text), width,
                                    span[0] + 1, span[1]))
    padded = record.ljust(span[1])
    return padded[:span[0]] + text + padded[span[1]:]


def number(record, span):
    text = field(record, span).strip()
    return float(text) if text else 0.0


# ---------------------------------------------------------------------------
# enlev.dat
# ---------------------------------------------------------------------------
class Enlev(object):
    """The rows of enlev.dat, addressable by IDEN2's level number."""

    def __init__(self, path):
        self.path = path
        self.records, self.endings = read_records(path)
        self.row_of = {}
        for k, rec in enumerate(self.records):
            if not rec.strip():
                continue
            text = field(rec, EN_INDEX).strip()
            if not text.isdigit():
                raise ValueError('%s line %d does not begin with a level '
                                 'number:\n    %s' % (path, k + 1, rec))
            self.row_of[int(text)] = k
        if not self.row_of:
            raise ValueError('%s holds no level rows' % path)

    def record(self, n):
        try:
            return self.records[self.row_of[n]]
        except KeyError:
            raise ValueError('%s has no level number %d' % (self.path, n))

    def e_obs(self, n):
        return number(self.record(n), EN_EOBS)

    def e_calc(self, n):
        return number(self.record(n), EN_ECALC)

    def known(self, n):
        return '*' in field(self.record(n), EN_STAR)

    def j(self, n):
        return field(self.record(n), EN_J).strip()

    def label(self, n):
        return self.record(n)[EN_LABEL:].strip().strip('/').strip()

    def set_measurement(self, n, unc, e_obs, known):
        """Give level ``n`` an observed energy, its uncertainty and the
        observed-minus-calculated that follows from its own calculated
        energy."""
        k = self.row_of[n]
        rec = self.records[k]
        rec = put(rec, EN_UNC, '%10.3f' % unc)
        rec = put(rec, EN_EOBS, '%12.3f' % e_obs)
        rec = put(rec, EN_STAR, STAR_ON if known else STAR_OFF)
        rec = put(rec, EN_OMC, '%12.3f' % (e_obs - number(rec, EN_ECALC)))
        self.records[k] = rec


# ---------------------------------------------------------------------------
# Translating an experimental identifier into an IDEN2 level number
# ---------------------------------------------------------------------------
def read_level_energies(path):
    """{identifier: energy} from a table of experimental level energies.

    Both shapes in use here are accepted: LOPT's tab-delimited output level
    list, whose columns are named ``Designation`` and ``Energy``, and a
    comma-separated table with a ``level_id`` column and an energy column
    called ``E``, ``E_input`` or ``Energy``.
    """
    with open(path, 'r', encoding='utf-8-sig', newline='') as fh:
        sample = fh.readline()
        fh.seek(0)
        delim = '\t' if sample.count('\t') > sample.count(',') else ','
        reader = csv.DictReader(fh, delimiter=delim)
        names = reader.fieldnames or []
        id_col = _pick(names, ('designation', 'level_id', 'id', 'level'))
        e_col = _pick(names, ('energy', 'e_input', 'e'))
        if id_col is None or e_col is None:
            raise ValueError(
                '%s has no identifier and energy columns; its columns are %s'
                % (path, ', '.join(names)))
        out = {}
        for row in reader:
            key = (row.get(id_col) or '').strip()
            text = (row.get(e_col) or '').strip()
            if not key or not text:
                continue
            try:
                out[key] = float(text)
            except ValueError:
                continue
    return out


def _pick(names, wanted):
    lowered = {n.strip().lower(): n for n in names if n}
    for w in wanted:
        if w in lowered:
            return lowered[w]
    return None


def read_map(path):
    """{identifier: IDEN2 level number} from a lookup table.

    A header naming the columns is used when present; otherwise the first
    two columns are taken as the identifier and the number, in that order.
    """
    out = {}
    with open(path, 'r', encoding='utf-8-sig', newline='') as fh:
        sample = fh.readline()
        fh.seek(0)
        delim = '\t' if sample.count('\t') > sample.count(',') else ','
        rows = list(csv.reader(fh, delimiter=delim))
    if not rows:
        raise ValueError('%s is empty' % path)
    start, id_col, n_col = 0, 0, 1
    head = [c.strip().lower() for c in rows[0]]
    if not any(c.replace('.', '', 1).isdigit() for c in rows[0] if c.strip()):
        start = 1
        for k, name in enumerate(head):
            if name in ('level_id', 'designation', 'id', 'level'):
                id_col = k
            elif any(w in name for w in ('iden', 'num', 'index', 'row')):
                n_col = k
    for row in rows[start:]:
        if len(row) <= max(id_col, n_col):
            continue
        key = row[id_col].strip()
        text = row[n_col].strip()
        if not key or not text:
            continue
        try:
            out[key] = int(float(text))
        except ValueError:
            continue
    if not out:
        raise ValueError('%s holds no identifier-to-number pairs' % path)
    return out


def row_number(text):
    """The enlev.dat row number a positional argument names, or None.

    The two mandatory arguments may be given either way round the same
    translation this script otherwise has to make: an experimental
    identifier such as ``059003.000483``, which IDEN2 does not know and
    which therefore has to be looked up, or the level's row number in
    ``enlev.dat``, which is what IDEN2 itself calls the level and needs no
    lookup at all.  A bare run of digits is a row number; anything else -
    an identifier always carries a dot - is an identifier.  The two forms
    may be mixed, one of each.
    """
    text = text.strip()
    return int(text) if text.isdigit() else None


def resolve(ids, args, enlev, log):
    """The two IDEN2 level numbers of the two levels named on the command line.

    ``--index`` names them outright and wins over everything.  Otherwise each
    of the two arguments is taken on its own: one given as a row number is
    already the answer, and one given as an identifier is looked up, in the
    table of ``--map`` if there is one and by energy otherwise.
    """
    if args.index:
        log('  level numbers    : %d and %d (from --index)' % tuple(args.index))
        return tuple(args.index)

    direct = [row_number(lid) for lid in ids]
    if all(n is not None for n in direct):
        log('  level numbers    : %d and %d (given as row numbers of %s)'
            % (direct[0], direct[1], os.path.basename(enlev.path)))
        return tuple(direct)

    if args.map:
        table = read_map(args.map)
        out = []
        for lid, n in zip(ids, direct):
            if n is not None:
                out.append(n)
                continue
            if lid not in table:
                raise ValueError('%s does not list %s' % (args.map, lid))
            out.append(table[lid])
        log('  level numbers    : %d and %d (from %s)'
            % (out[0], out[1], args.map))
        return tuple(out)

    if not os.path.exists(args.levels):
        raise ValueError(
            'no level-energy table at %s, so the identifiers cannot be '
            'translated into IDEN2 level numbers.  Give --levels, --map or '
            '--index, or name the levels by their row number in enlev.dat.'
            % args.levels)
    energies = read_level_energies(args.levels)
    out = []
    for lid, given in zip(ids, direct):
        if given is not None:
            out.append(given)
            continue
        if lid not in energies:
            raise ValueError('%s does not list %s' % (args.levels, lid))
        want = energies[lid]
        hits = [n for n in enlev.row_of
                if enlev.known(n) and abs(enlev.e_obs(n) - want) <= args.tol]
        if not hits:
            raise ValueError(
                'no row of %s has an observed energy within %g cm^-1 of '
                '%.4f, the energy of %s.  If the files have already been '
                'changed, name the rows with --index.'
                % (enlev.path, args.tol, want, lid))
        if len(hits) > 1:
            raise ValueError(
                'rows %s of %s all lie within %g cm^-1 of %.4f, the energy '
                'of %s; name the right one with --index.'
                % (', '.join(str(h) for h in sorted(hits)), enlev.path,
                   args.tol, want, lid))
        out.append(hits[0])
    log('  level numbers    : %d and %d (%s)'
        % (out[0], out[1],
           'matched by energy against %s' % args.levels
           if all(n is None for n in direct) else
           'one given as a row number, the other matched by energy against %s'
           % args.levels))
    return tuple(out)


def crossed(enlev, rows, ids, args, log):
    """Stop if the two rows already hold each other's energies.

    Only ``--index`` and ``--map`` can produce this, since matching by
    energy finds the row that holds the identifier's energy by definition.
    It is what running the script twice looks like, and it is worth
    catching, because the second run would quietly put everything back.
    """
    named_directly = (args.index or args.map
                      or any(row_number(i) is not None for i in ids))
    if not named_directly or not os.path.exists(args.levels):
        return
    try:
        energies = read_level_energies(args.levels)
    except ValueError:
        return
    if not all(i in energies for i in ids):
        return
    n1, n2 = rows
    id1, id2 = ids
    straight = (abs(enlev.e_obs(n1) - energies[id1]) <= args.tol
                and abs(enlev.e_obs(n2) - energies[id2]) <= args.tol)
    swapped = (abs(enlev.e_obs(n1) - energies[id2]) <= args.tol
               and abs(enlev.e_obs(n2) - energies[id1]) <= args.tol)
    if swapped and not straight:
        if not args.force:
            raise ValueError(
                'level %d already holds the energy %s has in %s, and level %d '
                'holds %s\'s.  The exchange looks as if it has been made '
                'already; running it again would put everything back.  Use '
                '--force if that is what you want.'
                % (n1, id2, args.levels, n2, id1))
        log('  --force: the two levels already hold each other\'s energies')


# ---------------------------------------------------------------------------
# trans.dat
# ---------------------------------------------------------------------------
class Trans(object):
    """The rows of trans.dat, indexed by (owner, partner)."""

    def __init__(self, path):
        self.path = path
        self.records, self.endings = read_records(path)
        self.header_of = {}
        self.row_of = {}
        owner = None
        for k, rec in enumerate(self.records):
            if not rec.strip():
                continue
            if rec.startswith('$'):
                owner = int(field(rec, TH_INDEX))
                self.header_of[owner] = k
            else:
                if owner is None:
                    raise ValueError('%s line %d comes before any level '
                                     'header' % (path, k + 1))
                partner = int(field(rec, TR_PARTNER))
                self.row_of[(owner, partner)] = k

    def partners_of(self, n):
        """Every level that ``n`` has a row for, whichever side it is on."""
        out = set()
        for owner, partner in self.row_of:
            if owner == n:
                out.add(partner)
            elif partner == n:
                out.add(owner)
        return out

    def row(self, a, b):
        """The row index of the transition between ``a`` and ``b``, or None."""
        if (a, b) in self.row_of:
            return self.row_of[(a, b)]
        return self.row_of.get((b, a))


def assignment(record):
    return record[TR_OBS:]


def has_line(obs):
    text = obs[TR_OBS_ROW[0]:TR_OBS_ROW[1]].strip()
    return bool(text) and text != '0'


def obs_wavenumber(obs):
    text = obs[TR_OBS_WN[0]:TR_OBS_WN[1]].strip()
    return float(text) if text else 0.0


def obs_omc(obs):
    text = obs[TR_OBS_OMC[0]:TR_OBS_OMC[1]].strip()
    return float(text) if text else 0.0


def obs_row(obs):
    text = obs[TR_OBS_ROW[0]:TR_OBS_ROW[1]].strip()
    return int(text) if text else 0


def set_omc(obs, value):
    return (obs[:TR_OBS_OMC[0]] + '%11.3f' % value
            + obs[TR_OBS_OMC[1]:])


# ---------------------------------------------------------------------------
# The exchange
# ---------------------------------------------------------------------------
def swap(enlev, trans, n1, n2, log):
    """Exchange the measurements of levels n1 and n2 in both files.

    Returns a report: the transferred lines, the ones that had nowhere to
    go, and anything that failed the observed-minus-predicted check.
    """
    old = {}
    for n in (n1, n2):
        old[n] = {'unc': number(enlev.record(n), EN_UNC),
                  'e': enlev.e_obs(n),
                  'known': enlev.known(n)}
    enlev.set_measurement(n1, old[n2]['unc'], old[n2]['e'], old[n2]['known'])
    enlev.set_measurement(n2, old[n1]['unc'], old[n1]['e'], old[n1]['known'])

    if trans is None:
        return {'moved': [], 'orphans': [], 'mismatch': [], 'rows': 0}

    if trans.row(n1, n2) is not None:
        raise ValueError(
            'trans.dat has a predicted transition between level %d and level '
            '%d.  Exchanging them would turn it into a transition of a level '
            'with itself, so this is not an interchange of the kind this '
            'script repairs.' % (n1, n2))

    # Pair up the rows that have to hand their assignment to each other.
    # A row of the pair is addressed by the level at its other end.
    pairs = []
    for p in sorted(trans.partners_of(n1) | trans.partners_of(n2)):
        if p in (n1, n2):
            continue
        pairs.append((p, trans.row(n1, p), trans.row(n2, p)))

    moved, orphans = [], []
    new_assignment = {}
    for p, k1, k2 in pairs:
        obs1 = assignment(trans.records[k1]) if k1 is not None else None
        obs2 = assignment(trans.records[k2]) if k2 is not None else None
        for source, target, from_level, to_level in (
                (obs1, k2, n1, n2), (obs2, k1, n2, n1)):
            if source is None:
                continue
            if target is None:
                if has_line(source):
                    orphans.append({'from': from_level, 'to': to_level,
                                    'partner': p,
                                    'wn': obs_wavenumber(source),
                                    'dlv': obs_row(source)})
                continue
            new_assignment[target] = source
            if has_line(source):
                moved.append({'from': from_level, 'to': to_level,
                              'partner': p,
                              'wn': obs_wavenumber(source),
                              'omc': obs_omc(source),
                              'dlv': obs_row(source)})
        for k in (k1, k2):
            if k is not None and k not in new_assignment:
                new_assignment[k] = BLANK_OBS

    # Every row that mentions either level: new energies, new predicted
    # wavenumbers, and the assignment it has just been handed.
    touched = set(new_assignment)
    for owner, partner in trans.row_of:
        if owner in (n1, n2) or partner in (n1, n2):
            touched.add(trans.row_of[(owner, partner)])

    owner_of = {}
    for (owner, partner), k in trans.row_of.items():
        owner_of[k] = (owner, partner)

    # The block headers of the two levels carry their observed energy.
    for n in (n1, n2):
        k = trans.header_of.get(n)
        if k is None:
            continue
        rec = trans.records[k]
        rec = put(rec, TH_EOBS, '%13.3f' % enlev.e_obs(n))
        rec = put(rec, TH_STAR, STAR_ON if enlev.known(n) else STAR_OFF)
        trans.records[k] = rec

    mismatch = []
    for k in sorted(touched):
        owner, partner = owner_of[k]
        rec = trans.records[k]
        wn_pred = abs(enlev.e_obs(owner) - enlev.e_obs(partner))
        rec = put(rec, TR_EPART, '%12.3f' % enlev.e_obs(partner))
        rec = put(rec, TR_STAR,
                  STAR_ON if enlev.known(partner) else STAR_OFF)
        rec = put(rec, TR_WN, '%13.3f' % wn_pred)
        obs = new_assignment.get(k, assignment(rec))
        if has_line(obs):
            was = obs_omc(obs)
            now = obs_wavenumber(obs) - wn_pred
            if abs(now - was) > 5e-4:
                mismatch.append({'row': k + 1, 'owner': owner,
                                 'partner': partner, 'was': was, 'now': now})
            obs = set_omc(obs, now)
        else:
            obs = BLANK_OBS
        trans.records[k] = rec[:TR_OBS] + obs

    return {'moved': moved, 'orphans': orphans, 'mismatch': mismatch,
            'rows': len(touched)}


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def parse_args(argv):
    p = argparse.ArgumentParser(
        description='Exchange two levels\' observed energies and observed '
                    'lines in the IDEN2 working files, keeping IDEN2\'s level '
                    'numbering as it is.  The two observed energies in '
                    'enlev.dat are always exchanged; --index and --map only '
                    'save the script the work of finding the two rows.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('id1', help='the first level: its experimental '
                                'identifier, or its row number in enlev.dat')
    p.add_argument('id2', help='the second level: its experimental '
                               'identifier, or its row number in enlev.dat')
    p.add_argument('--iden2-dir', default=swap_paths.working_path(NAME_IDEN2),
                   help='the directory holding enlev.dat and trans.dat')
    p.add_argument('--index', nargs=2, type=int, metavar=('N1', 'N2'),
                   help='the two levels\' row numbers in enlev.dat, given '
                        'directly instead of being looked up.  This only '
                        'saves the lookup: the two observed energies, their '
                        'uncertainties and their line assignments are '
                        'exchanged in enlev.dat and trans.dat either way')
    p.add_argument('--map',
                   help='a lookup table of experimental identifier and '
                        'IDEN2 level number')
    p.add_argument('--levels', default=swap_paths.working_path(NAME_LEVELS),
                   help='a table of experimental level energies, used to '
                        'find each identifier in enlev.dat by its energy')
    p.add_argument('--tol', type=float, default=DEF_TOL,
                   help='how far, in cm^-1, an energy in --levels may lie '
                        'from the one in enlev.dat and still be the same level')
    p.add_argument('--force', action='store_true',
                   help='exchange two levels of different J, which is '
                        'normally refused')
    p.add_argument('--out-dir',
                   help='write the rewritten files here instead of over the '
                        'originals')
    p.add_argument('--no-backup', action='store_true',
                   help='do not keep the previous contents in <file>.bak')
    p.add_argument('--dry-run', action='store_true',
                   help='report what would change and write nothing')
    return p.parse_args(argv)


def find_file(directory, name):
    """The file called ``name`` in ``directory``, whatever its case."""
    if not os.path.isdir(directory):
        raise ValueError('there is no directory %s' % directory)
    for entry in os.listdir(directory):
        if entry.lower() == name:
            return os.path.join(directory, entry)
    return None


def main(argv=None):
    args = parse_args(argv)
    log = print
    id1, id2 = args.id1.strip(), args.id2.strip()
    if id1 == id2:
        log('The two identifiers are the same; there is nothing to exchange.')
        return 2

    enlev_path = find_file(args.iden2_dir, 'enlev.dat')
    if enlev_path is None:
        raise ValueError('there is no enlev.dat in %s' % args.iden2_dir)
    trans_path = find_file(args.iden2_dir, 'trans.dat')

    log('swap_line_assignments_IDEN.py')
    log('  IDEN2 directory  : %s'
        % swap_paths.describe(NAME_IDEN2, args.iden2_dir))
    log('  exchanging       : %s  <->  %s' % (id1, id2))

    enlev = Enlev(enlev_path)
    n1, n2 = resolve((id1, id2), args, enlev, log)
    if n1 == n2:
        raise ValueError('both identifiers resolve to level %d' % n1)

    for n, lid in ((n1, id1), (n2, id2)):
        log('  %-13s = level %4d  %-11s J = %-4s  E_obs = %12.3f  '
            'E_calc = %11.3f  O-C = %9.3f'
            % (lid, n, enlev.label(n), enlev.j(n), enlev.e_obs(n),
               enlev.e_calc(n), enlev.e_obs(n) - enlev.e_calc(n)))

    if enlev.j(n1) != enlev.j(n2):
        if not args.force:
            raise ValueError(
                'level %d has J = %s and level %d has J = %s.  Two levels of '
                'different J cannot have had their identities interchanged, '
                'so one of the two identifiers is probably not the level you '
                'meant.  Use --force if you really mean it.'
                % (n1, enlev.j(n1), n2, enlev.j(n2)))
        log('  --force: exchanging levels of different J')
    crossed(enlev, (n1, n2), (id1, id2), args, log)
    if not (enlev.known(n1) and enlev.known(n2)):
        raise ValueError(
            'level %d and level %d must both have an observed energy; '
            'level %d has none.'
            % (n1, n2, n1 if not enlev.known(n1) else n2))

    trans = Trans(trans_path) if trans_path else None
    if trans is None:
        log('  no trans.dat in %s; only enlev.dat will be rewritten'
            % args.iden2_dir)

    info = swap(enlev, trans, n1, n2, log)

    log('')
    log('After the exchange:')
    for n, lid in ((n1, id1), (n2, id2)):
        log('  level %4d  %-11s E_obs = %12.3f   O-C = %9.3f'
            % (n, enlev.label(n), enlev.e_obs(n),
               enlev.e_obs(n) - enlev.e_calc(n)))

    log('')
    log('trans.dat: %d rows rewritten, %d observed lines transferred.'
        % (info['rows'], len(info['moved'])))
    if info['moved']:
        by_level = {}
        for m in info['moved']:
            by_level.setdefault(m['to'], []).append(m)
        for n in sorted(by_level):
            lines = by_level[n]
            log('  level %d receives %d line(s), O-C from %+.3f to %+.3f'
                % (n, len(lines),
                   min(m['omc'] for m in lines),
                   max(m['omc'] for m in lines)))
    if info['orphans']:
        log('')
        log('%d observed line(s) could not be transferred: the level they '
            'move to has no predicted transition to the level at the other '
            'end, so trans.dat has no row to hold them.  Their assignment '
            'has been removed; they are still in dlv.dat and can be '
            're-identified by hand in IDEN2.' % len(info['orphans']))
        log('  %-8s %-8s %-8s %14s %8s'
            % ('was on', 'would be', 'other end', 'wavenumber', 'dlv row'))
        for o in info['orphans']:
            log('  %-8d %-8d %-8d %14.3f %8d'
                % (o['from'], o['to'], o['partner'], o['wn'], o['dlv']))
    if info['mismatch']:
        raise ValueError(
            'the observed-minus-predicted value of %d transferred line(s) '
            'changed, which cannot happen when two levels simply exchange '
            'their energies and their lines.  The first is row %d, the '
            'transition %d-%d, where it went from %+.3f to %+.3f.  The files '
            'are not in the state this rewrite assumes; nothing has been '
            'written.'
            % (len(info['mismatch']), info['mismatch'][0]['row'],
               info['mismatch'][0]['owner'], info['mismatch'][0]['partner'],
               info['mismatch'][0]['was'], info['mismatch'][0]['now']))

    if args.dry_run:
        log('')
        log('--dry-run: nothing was written.')
        return 0

    out_dir = args.out_dir or args.iden2_dir
    if args.out_dir and not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)
    log('')
    for path, obj in ((enlev_path, enlev), (trans_path, trans)):
        if obj is None:
            continue
        target = os.path.join(out_dir, os.path.basename(path))
        if not args.out_dir and not args.no_backup:
            saved, saved_eol = read_records(path)
            # The originals are still on disk: read them back before writing.
            write_records(path + '.bak', saved, saved_eol)
            log('Previous contents kept in %s' % (path + '.bak'))
        write_records(target, obj.records, obj.endings)
        log('Written: %s' % target)

    log('')
    log('dlv.dat and numset.dat carry no level information and were not '
        'touched.  Open the two levels in IDEN2 and check the branching '
        'patterns before accepting the exchange.')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except ValueError as exc:
        sys.stderr.write('\nswap_line_assignments_IDEN.py: %s\n' % exc)
        sys.exit(1)
