"""
swap_line_assignments_pipeline.py
=================================
Record an exchange of two levels' measurements in the files the
classification pipeline reads, so that re-running the pipeline reproduces it
instead of undoing it.

What this is for
----------------
``swap_line_assignments_LOPT.py`` and ``swap_line_assignments_IDEN.py``
exchange two levels' observed lines - and with them their observed energies -
in the LOPT working files and in the IDEN2 working files.  Neither of those
is an input of the classification pipeline.  ``make_LOPT_input.py`` writes
the LOPT transitions file out of ``line_classifications.csv``, which
``classify_lines.py`` produces from the published level list, the observed
line list and the calculated transitions.  Run the pipeline again and the
transitions file is rebuilt from those, and the exchange is gone.

This script is what makes it stay.  Two things have to be written.

**The energies.**  ``classify_lines.py`` gives every level the energy of the
published list, overridden by ``revised_level_energies.csv``, and generates
that level's candidate transitions at that energy.  The exchange moved both
energies, so both belong in that overlay: one row each, the identifier
keeping its own configuration and its own calculated transition
probabilities and taking the other's measured position.

**The decision ledger.**  ``line_decisions.csv`` holds the verdicts reached
by hand - accept this assignment, reject that one - each keyed by an
observed wavenumber and the two level identifiers.  Those identifiers are
written out, so every order naming one of the exchanged levels now points at
the wrong one; each is re-keyed here, and the reason it carries is annotated
with the date of the exchange rather than being rewritten, because the
reason is prose and only its subject has changed.

**And the legacy identifications.**  This is the part that used to be easy
to miss and expensive to get wrong.  An identification taken from the
published line list - Sugar's own, the ones ``line_classifications.csv``
marks ``new`` = 0 - names its two levels explicitly in the line workbook, an
input file that is never edited.  Step 3 of ``classify_lines.py`` keeps such
an identification unless something rejects it ("old, no solid evidence for
rejection"), and it used to go on naming its level after the exchange had
moved that level a hundred wavenumbers away.  The level then ended up fitted
between its true lines and its stale legacy ones.  That is what pulled level
``059003.000424`` 81 cm^-1 off its LOPT position on 2026-09-05.  So this
script wrote two orders for every one of them: a ``reject`` where the line
stood and an ``accept`` under the other identifier.

``classify_lines.py`` now does that itself.  It reads the ``comment`` column
of ``revised_level_energies.csv`` - the same comment this script writes,
which says in plain words that the two levels were swapped - and moves every
published identification of either level to the other identifier before the
matching begins (``retag_legacy_identifications``): the pair the workbook
names is taken off the observed line, and the same identification is seeded
under the partner, which is the level that now sits where the line is.  The
ledger orders would only say a second time what the overlay already says,
and their ``reject`` halves would then match no candidate at all and be
reported, one line each, as unapplied decisions on every future run.  So
they are **no longer written** when this script writes the overlay.  They
are still written when it does not (``--no-overrides``), because then
nothing else records the exchange, and ``--legacy-orders`` asks for them
outright.

What is NOT written
-------------------
``Icalc.xlsx`` and ``icalc_new.xlsx`` are left exactly as they are, and so
is the published line workbook.  The calculated transitions are keyed by the
pair of level identifiers, and an identifier keeps its configuration, its
term and its calculated transition probabilities through the exchange - what
moves is the measurement, not the calculation.  ``classify_lines.py`` reads
only the ``Icalc`` and ``u%gA`` columns of the main file, so a level's
calculated intensities do not even depend on where the level sits; the
predicted wavenumbers are computed from the level energies, which is exactly
what the overlay changes.

Where the new energies come from
--------------------------------
From ``LOPT_output_levels.txt``, LOPT's output level table, and it does not
matter whether the LOPT run behind it was made before or after the exchange.
No LOPT run is needed to record the exchange at all.

An exchange does not invent energies.  Both were already measured; what
changes is which identifier each belongs to.  So the two numbers this script
has to write are the two the level table already holds, the other way round.
A LOPT run made after the exchange writes them that way itself - each
identifier refitted with the other's lines, which differs from the plain
exchange of the two old numbers by a few thousandths of a wavenumber - and
then they are taken as they stand.

The script works out which of the two cases it is in, by comparing the table
with the energies the pipeline currently has (the ``low_E`` and ``upp_E``
columns of ``line_classifications.csv``):

* the table still agrees with the pipeline - the LOPT run has not been made
  yet, and the two energies are exchanged here;
* the table has them **crossed** with respect to the pipeline - the LOPT run
  has been made, and each identifier keeps the energy the table gives it.

Either way the result is the same exchange, to within those few thousandths,
so running the script before the LOPT run and again after it is the ordinary
way to work: the second run replaces the two energies with the fitted ones
and recognises everything else it wrote the first time.  It stops only in the
third case - a table that is neither, meaning something other than the
exchange has moved these levels.  ``--exchange`` and ``--as-fitted`` force
one of the two readings when the pipeline has no energy to compare against.

Usage
-----
    python swap_line_assignments_pipeline.py 059003.000483 059003.000398
    python swap_line_assignments_pipeline.py id1 id2 --dry-run
    python swap_line_assignments_pipeline.py id1 id2 --legacy-orders

Every default file name is looked for first in the directory the command was
run from and then in the directory holding this script, so that working in an
iteration folder acts on that folder's copies while the files that exist in
only one place are still found.  Files are rewritten in place with their
previous contents kept in ``<file>.bak`` (``--no-backup`` turns that off);
``--dry-run`` reports and writes nothing.

``swap_line_assignments.py`` runs this script together with its two
companions in one command.
"""

import argparse
import csv
import datetime
import os
import sys

import output_files
import swap_paths

NAME_LEVELS = 'LOPT_output_levels.txt'
NAME_OVERRIDES = 'revised_level_energies.csv'
NAME_LEDGER = 'line_decisions.csv'
NAME_CLASSIFICATIONS = 'line_classifications.csv'

# How far, in cm^-1, two energies may differ and still be called the same.
# The levels an exchange can be made between are separated by far more than
# this - they have to be resolvable as two levels at all - and a LOPT refit
# moves an energy by a few thousandths.
DEF_TOL = 0.5

ACCEPT = 'accept'
REJECT = 'reject'


def mark(id1, id2):
    """The phrase every row this script writes carries, and reads back.

    Running the script twice on the same pair must not undo the first run:
    re-keying a ledger order a second time would point it back at the level
    it started from.  So each row it writes or annotates says, in plain
    words, that these two levels have been exchanged, and each row that
    already says so is left alone.  The two identifiers are put in a fixed
    order here, so that the phrase does not depend on which of them was
    given first, and the date is not part of it, so that a later run on a
    later date still recognises the row.
    """
    a, b = sorted((id1, id2))
    return 'levels %s and %s were swapped' % (a, b)


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------
def read_lopt_levels(path):
    """{level_id: energy} from LOPT's output level table.

    The file is tab-separated with a header row; the first column is the
    level identifier and the second its fitted energy in cm^-1.
    """
    out = {}
    with open(path, 'r', encoding='latin-1', newline='') as fh:
        for n, raw in enumerate(fh, start=1):
            fields = raw.rstrip('\r\n').split('\t')
            if len(fields) < 2:
                continue
            key, text = fields[0].strip(), fields[1].strip()
            if not key or n == 1:
                continue
            try:
                out[key] = float(text)
            except ValueError:
                continue
    if not out:
        raise ValueError('%s holds no level energies' % path)
    return out


def read_table(path):
    """(fieldnames, rows) of a csv, the rows as dicts in the order written."""
    with open(path, 'r', encoding='utf-8-sig', newline='') as fh:
        rdr = csv.DictReader(fh)
        return list(rdr.fieldnames or []), [dict(r) for r in rdr]


def write_table(path, fieldnames, rows):
    """Write a csv with LF line endings, which is what this repository keeps."""
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, '') for k in fieldnames})


def pipeline_energies(path, ids):
    """The energy the pipeline currently has for each of the two levels.

    Taken from ``line_classifications.csv``, whose ``low_E`` and ``upp_E``
    columns are the energies every candidate transition was generated at, so
    they are the published energy already overridden by
    ``revised_level_energies.csv`` - which is what has to be compared with a
    LOPT table to see whether the exchange has been made.  A level that
    appears in no row of the file is simply absent from the answer.
    """
    want = set(ids)
    out = {}
    with open(path, 'r', encoding='utf-8-sig', newline='') as fh:
        for row in csv.DictReader(fh):
            for end, col in (('low_id', 'low_E'), ('upp_id', 'upp_E')):
                lid = (row.get(end) or '').strip()
                if lid in want and lid not in out:
                    try:
                        out[lid] = float(row.get(col) or '')
                    except ValueError:
                        pass
        return out


def wavenumber(text):
    """An observed wavenumber as the ledger writes it: four decimals.

    The classification table carries the full binary expansion of the
    number, which is unreadable and pointless here: classify_lines.py finds
    a ledger row's line by taking the nearest observed line within
    0.01 cm^-1, and the closest two lines of this spectrum are 0.10 cm^-1
    apart, so four decimals identify a line a thousand times over.
    """
    text = (text or '').strip()
    try:
        return '%.4f' % float(text)
    except ValueError:
        return text


def legacy_assignments(path, ids):
    """The accepted legacy identifications of the two levels.

    "Legacy" is ``new`` = 0: the pair of levels stands against this line in
    the published line workbook, so the identification is Sugar's and comes
    back on every run whatever the level energies are.  Returns one dict per
    identification with its wavenumber as written, its two identifiers, and
    which of the two exchanged levels it is an identification of.
    """
    out = []
    with open(path, 'r', encoding='utf-8-sig', newline='') as fh:
        for row in csv.DictReader(fh):
            low = (row.get('low_id') or '').strip()
            upp = (row.get('upp_id') or '').strip()
            mine = [i for i in ids if i in (low, upp)]
            if not mine:
                continue
            if (row.get('accepted') or '').strip() not in ('1', '1.0'):
                continue
            if (row.get('new') or '').strip() not in ('0', '0.0'):
                continue
            out.append({'wn': wavenumber(row.get('wn_key')
                                         or row.get('wn_obs')),
                        'low': low, 'upp': upp, 'level': mine[0],
                        'intens': (row.get('obs_intens') or '').strip()})
    return out


# ---------------------------------------------------------------------------
# The energies
# ---------------------------------------------------------------------------
def target_energies(ids, lopt, current, exchange, as_fitted, tol):
    """The energy each identifier must be given, and how it was arrived at.

    The two energies to be written are always the two the level table holds,
    the only question being whether the table already has them the right way
    round.  It is answered by comparing the table with the energies the
    pipeline currently has: if the two have crossed, the table is from a LOPT
    run made after the exchange and each identifier keeps what the table
    gives it; if they still agree, no such run has been made and the two are
    exchanged here.  A table that is neither stops the script, because
    something other than the exchange has then moved these levels.

    ``exchange`` and ``as_fitted`` force one of the two readings, for the case
    where the pipeline has no energy for either level and the comparison
    cannot be made.  Each is refused if the table plainly says the other.
    """
    id1, id2 = ids
    for lid in ids:
        if lid not in lopt:
            raise ValueError('the level table does not list %s' % lid)
    e1, e2 = lopt[id1], lopt[id2]
    if abs(e1 - e2) <= tol:
        raise ValueError(
            '%s and %s are given the same energy (%.4f and %.4f cm^-1) in the '
            'level table; there is nothing to exchange.' % (id1, id2, e1, e2))

    known = all(lid in current for lid in ids)
    crossed = known and (abs(e1 - current[id2]) <= tol
                         and abs(e2 - current[id1]) <= tol)
    straight = known and (abs(e1 - current[id1]) <= tol
                          and abs(e2 - current[id2]) <= tol)

    HERE = 'the two energies of the level table, exchanged here'
    FITTED = 'the level table, in which they have already been exchanged'

    if exchange and as_fitted:
        raise ValueError('--exchange and --as-fitted contradict each other.')
    if exchange:
        if crossed:
            raise ValueError(
                'the level table already holds the two energies exchanged '
                '(%s at %.4f, %s at %.4f), so --exchange would put them back. '
                'Drop --exchange to record the table as it stands.'
                % (id1, e1, id2, e2))
        return {id1: e2, id2: e1}, HERE
    if as_fitted:
        if straight:
            raise ValueError(
                'the level table still holds the energies the pipeline has '
                '(%s at %.4f, %s at %.4f), so --as-fitted would write them '
                'back unchanged and record no exchange.'
                % (id1, e1, id2, e2))
        return {id1: e1, id2: e2}, FITTED

    if crossed:
        return {id1: e1, id2: e2}, FITTED
    if straight:
        return {id1: e2, id2: e1}, HERE
    if known:
        raise ValueError(
            'the level table gives %s %.4f and %s %.4f, which is neither what '
            'the pipeline has for them (%.4f and %.4f) nor those two values '
            'exchanged.  Something other than the exchange has moved these '
            'levels; decide by hand what their energies should be.'
            % (id1, e1, id2, e2, current[id1], current[id2]))
    raise ValueError(
        'the pipeline has no energy for %s or for %s, so the level table '
        'cannot be checked and there is no telling whether its two energies '
        'still have to be exchanged.  Pass --exchange to exchange them here, '
        'or --as-fitted to record the table as it stands.' % (id1, id2))


def override_rows(rows, targets, ids, date):
    """Insert or update the two rows of the level-override file.

    Returns the new list of rows and one report line per level.  A level
    already listed keeps its place in the file and has its energy replaced,
    with the old comment kept and the exchange added to it; the alternative -
    a second row for the same level - would be read as a contradiction.
    """
    id1, id2 = ids
    note = ('exchanged with %%s on %%s (%s): this identifier keeps its own '
            'configuration and its own calculated transition probabilities '
            'and takes the other\'s measured position'
            % mark(id1, id2))
    by_id = {}
    for row in rows:
        by_id.setdefault((row.get('level_id') or '').strip(), row)

    out, report = list(rows), []
    for lid in ids:
        other = id2 if lid == id1 else id1
        energy = '%.4f' % targets[lid]
        comment = note % (other, date)
        row = by_id.get(lid)
        if row is None:
            out.append({'level_id': lid, 'E_input': energy,
                        'comment': comment})
            report.append('  %s  new row at %s cm^-1' % (lid, energy))
            continue
        was = (row.get('E_input') or '').strip()
        old_comment = (row.get('comment') or '').strip()
        recorded = mark(id1, id2) in old_comment
        row['E_input'] = energy
        if not recorded:
            row['comment'] = ((old_comment + '. ' + comment) if old_comment
                              else comment)
        report.append('  %s  %s -> %s cm^-1%s'
                      % (lid, was or '(blank)', energy,
                         '  (the exchange was already recorded here; only '
                         'the energy is refreshed)' if recorded else ''))
    return out, report


# ---------------------------------------------------------------------------
# The decision ledger
# ---------------------------------------------------------------------------
def swap_id(lid, id1, id2):
    return id2 if lid == id1 else (id1 if lid == id2 else lid)


def ledger_key(row):
    """The key classify_lines.py reads a ledger row by."""
    try:
        wn = float((row.get('wn_key') or '').strip())
    except ValueError:
        wn = None
    return (wn, (row.get('low_id') or '').strip(),
            (row.get('upp_id') or '').strip())


def rekey_ledger(rows, id1, id2, date):
    """Point every order naming one of the two levels at the other one.

    The reason is left as it was and the exchange appended to it, in the
    words the analysis uses for it, because the reason is prose about a
    measured line and the line has not changed - only which identifier it
    is filed under.
    """
    marker = mark(id1, id2)
    sentence = '; Lines belonging to %s on %s.' % (marker, date)
    touched, skipped = [], []
    for n, row in enumerate(rows, start=2):
        low = (row.get('low_id') or '').strip()
        upp = (row.get('upp_id') or '').strip()
        if low not in (id1, id2) and upp not in (id1, id2):
            continue
        reason = (row.get('reason') or '').rstrip()
        if marker in reason:
            skipped.append((n, (row.get('wn_key') or '').strip(), low, upp))
            continue
        was = (low, upp)
        row['low_id'] = swap_id(low, id1, id2)
        row['upp_id'] = swap_id(upp, id1, id2)
        row['reason'] = (reason + sentence) if reason else sentence[2:].strip()
        touched.append((n, (row.get('wn_key') or '').strip(), was,
                        (row['low_id'], row['upp_id']),
                        (row.get('decision') or '').strip()))
    return touched, skipped


def legacy_orders(legacy, id1, id2, date, with_accepts):
    """The orders that move each legacy identification to the other level.

    One ``reject`` where the line stands - the level it names is about to be
    a hundred wavenumbers away from it - and one ``accept`` under the other
    identifier, which is the level that now sits where the line is.
    """
    out = []
    for item in legacy:
        lid = item['level']
        other = swap_id(lid, id1, id2)
        low_new = swap_id(item['low'], id1, id2)
        upp_new = swap_id(item['upp'], id1, id2)
        out.append({
            'wn_key': item['wn'], 'low_id': item['low'],
            'upp_id': item['upp'], 'decision': REJECT, 'date': date,
            'reason': ('identification of the published line list, at the '
                       'position %s held before the exchange; lines belonging '
                       'to %s on %s, so this line now belongs to %s'
                       % (lid, mark(id1, id2), date, other))})
        if with_accepts:
            out.append({
                'wn_key': item['wn'], 'low_id': low_new,
                'upp_id': upp_new, 'decision': ACCEPT, 'date': date,
                'reason': ('the same identification of the published line '
                           'list, re-keyed: lines belonging to %s on %s, and '
                           '%s is the level that now sits where this line is'
                           % (mark(id1, id2), date, other))})
    return out


def merge_orders(rows, new_rows):
    """Add the new orders, refusing to contradict one already in the file.

    ``classify_lines.py`` stops on two rows that rule differently on the same
    assignment, so a contradiction made here would abort the next run rather
    than being noticed at the end of this one.  A row that says the same
    thing as one already there is dropped instead of being duplicated.
    """
    have = {}
    for row in rows:
        have.setdefault(ledger_key(row), (row.get('decision') or '').strip())
    added, already, clashes = [], [], []
    for row in new_rows:
        key = ledger_key(row)
        if key in have:
            if have[key] == row['decision']:
                already.append(row)
            else:
                clashes.append((row, have[key]))
            continue
        have[key] = row['decision']
        rows.append(row)
        added.append(row)
    if clashes:
        row, standing = clashes[0]
        raise ValueError(
            'the ledger already orders %s for the assignment %s - %s at %s, '
            'and this exchange would order %s.  %d order(s) clash like this. '
            'The exchange cannot be recorded until they are settled by hand: '
            'classify_lines.py stops on two rows that rule differently on the '
            'same assignment.'
            % (standing, row['low_id'], row['upp_id'], row['wn_key'],
               row['decision'], len(clashes)))
    return added, already


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def parse_args(argv):
    p = argparse.ArgumentParser(
        description='Record an exchange of two levels\' measurements in the '
                    'files the classification pipeline reads, so that a '
                    're-run reproduces it instead of undoing it.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('id1', help='identifier of the first level')
    p.add_argument('id2', help='identifier of the second level')
    p.add_argument('--levels', default=swap_paths.working_path(NAME_LEVELS),
                   help='LOPT\'s output level table, read for the two '
                        'energies; a run made before or after the exchange '
                        'will do, and which it is is worked out here')
    p.add_argument('--overrides',
                   default=swap_paths.working_path(NAME_OVERRIDES),
                   help='the revised-energies overlay the pipeline reads')
    p.add_argument('--ledger', default=swap_paths.working_path(NAME_LEDGER),
                   help='the decision ledger')
    p.add_argument('--classifications',
                   default=swap_paths.working_path(NAME_CLASSIFICATIONS),
                   help='the classification table, read for the energies the '
                        'pipeline currently has and for the legacy '
                        'identifications of the two levels')
    p.add_argument('--exchange', action='store_true',
                   help='take the level table to be from a run made before '
                        'the exchange and exchange its two energies here, '
                        'instead of working out which it is')
    p.add_argument('--as-fitted', action='store_true',
                   help='take the level table to be from a run made after '
                        'the exchange and record its two energies as they '
                        'stand, instead of working out which it is')
    p.add_argument('--legacy-orders', action='store_true',
                   help='write a reject and an accept into the ledger for '
                        'every identification of the published line list '
                        'held by either level.  Not needed since '
                        'classify_lines.py began carrying those '
                        'identifications over itself, from the exchange '
                        'recorded in the comment of the level overrides; '
                        'they are written anyway when --no-overrides means '
                        'nothing records the exchange')
    p.add_argument('--rejects-only', action='store_true',
                   help='with --legacy-orders, write only the reject at the '
                        'level the identification names, leaving the accept '
                        'at the other level to the classification')
    p.add_argument('--no-ledger', action='store_true',
                   help='leave the decision ledger alone')
    p.add_argument('--no-overrides', action='store_true',
                   help='leave the revised-energies overlay alone')
    p.add_argument('--date',
                   help='the date written into the new rows')
    p.add_argument('--tol', type=float, default=DEF_TOL,
                   help='how far, in cm^-1, two energies may differ and still '
                        'be called the same')
    p.add_argument('--no-backup', action='store_true',
                   help='do not keep the previous contents in <file>.bak')
    p.add_argument('--dry-run', action='store_true',
                   help='report what would change and write nothing')
    return p.parse_args(argv)


def today():
    now = datetime.date.today()
    return '%d/%d/%d' % (now.month, now.day, now.year)


def save(path, rows, fieldnames, backup, log):
    if backup and os.path.exists(path):
        with open(path, 'r', encoding='utf-8-sig', newline='') as fh:
            text = fh.read()
        with open(path + '.bak', 'w', encoding='utf-8', newline='') as fh:
            fh.write(text)
        log('Previous contents kept in %s' % (path + '.bak'))
    write_table(path, fieldnames, rows)
    log('Written: %s' % path)


def main(argv=None):
    args = parse_args(argv)
    log = print
    id1, id2 = args.id1.strip(), args.id2.strip()
    if id1 == id2:
        log('The two identifiers are the same; there is nothing to exchange.')
        return 2
    for lid in (id1, id2):
        if lid.isdigit():
            raise ValueError(
                '%s is a row number of enlev.dat.  The pipeline does not know '
                'IDEN2\'s numbering; name both levels by their experimental '
                'identifier here.' % lid)
    date = args.date or today()

    # Both files are csv, and both are read in Excel; a run that found one of
    # them locked would stop half way, having written the other.
    if not args.dry_run:
        output_files.require_writable([args.overrides, args.ledger], 'ledger')

    log('swap_line_assignments_pipeline.py')
    log('  exchanging       : %s  <->  %s' % (id1, id2))
    log('  level table      : %s'
        % swap_paths.describe(NAME_LEVELS, args.levels))
    log('  level overrides  : %s'
        % swap_paths.describe(NAME_OVERRIDES, args.overrides))
    log('  decision ledger  : %s'
        % swap_paths.describe(NAME_LEDGER, args.ledger))
    log('  classifications  : %s'
        % swap_paths.describe(NAME_CLASSIFICATIONS, args.classifications))
    log('  date written     : %s' % date)

    if not os.path.exists(args.levels):
        raise ValueError('there is no level table at %s' % args.levels)
    lopt = read_lopt_levels(args.levels)
    current = ({} if not os.path.exists(args.classifications)
               else pipeline_energies(args.classifications, (id1, id2)))
    if not current:
        log('')
        log('  %s gives no energy for either level, so the level table cannot '
            'be checked against what the pipeline has.'
            % os.path.basename(args.classifications))

    targets, whence = target_energies((id1, id2), lopt, current,
                                      args.exchange, args.as_fitted, args.tol)
    log('')
    log('Energies (%s):' % whence)
    for lid in (id1, id2):
        was = current.get(lid)
        log('  %-13s  %s -> %12.4f cm^-1  (%s)'
            % (lid, ('%12.4f' % was) if was is not None else ' ' * 8 + 'n/a',
               targets[lid],
               ('%+.4f' % (targets[lid] - was)) if was is not None
               else 'not in the classification table'))

    # --- the overlay ------------------------------------------------------
    ov_rows = ov_fields = None
    if args.no_overrides:
        log('')
        log('The level overrides are not touched (--no-overrides).')
    else:
        if not os.path.exists(args.overrides):
            ov_fields, ov_rows = ['level_id', 'E_input', 'comment'], []
            log('')
            log('%s does not exist; it will be created.' % args.overrides)
        else:
            ov_fields, ov_rows = read_table(args.overrides)
            for col in ('level_id', 'E_input'):
                if col not in ov_fields:
                    raise ValueError('%s has no %s column'
                                     % (args.overrides, col))
            if 'comment' not in ov_fields:
                ov_fields.append('comment')
        ov_rows, ov_report = override_rows(ov_rows, targets, (id1, id2),
                                           date)
        log('')
        log('%s:' % os.path.basename(args.overrides))
        for line in ov_report:
            log(line)

    # --- the ledger -------------------------------------------------------
    led_rows = led_fields = None
    if args.no_ledger:
        log('')
        log('The decision ledger is not touched (--no-ledger).')
    else:
        if not os.path.exists(args.ledger):
            raise ValueError('there is no decision ledger at %s'
                             % args.ledger)
        led_fields, led_rows = read_table(args.ledger)
        for col in ('wn_key', 'low_id', 'upp_id', 'decision'):
            if col not in led_fields:
                raise ValueError('%s has no %s column' % (args.ledger, col))
        for col in ('date', 'reason'):
            if col not in led_fields:
                led_fields.append(col)

        touched, skipped = rekey_ledger(led_rows, id1, id2, date)
        log('')
        log('%s: %d existing order(s) re-keyed.'
            % (os.path.basename(args.ledger), len(touched)))
        for n, wn, was, now, decision in touched:
            log('  row %-4d %-14s %s - %s  ->  %s - %s   %s'
                % (n, wn, was[0], was[1], now[0], now[1], decision))
        if skipped:
            log('  %d order(s) already say this exchange has been made and '
                'were left as they are - re-keying them a second time would '
                'point them back at the level they started from:'
                % len(skipped))
            for n, wn, low, upp in skipped:
                log('    row %-4d %-14s %s - %s' % (n, wn, low, upp))

        legacy = []
        if not os.path.exists(args.classifications):
            log('  %s is not there, so the legacy identifications of the two '
                'levels could not be looked up.  Check by hand that neither '
                'level keeps an identification of the published line list at '
                'its old position.' % args.classifications)
        else:
            legacy = legacy_assignments(args.classifications, (id1, id2))
        # The overlay records the exchange, and classify_lines.py reads it
        # and moves the published identifications itself; an order here would
        # only repeat that, and its reject half would match no candidate.
        write_orders = args.legacy_orders or args.no_overrides
        log('')
        log('Identifications of the published line list held by the two '
            'levels: %d.' % len(legacy))
        for item in legacy:
            log('  %-14s %s - %s   (of %s)'
                % (item['wn'], item['low'], item['upp'], item['level']))
        if not write_orders:
            log('  No order is written for them.  classify_lines.py reads the '
                'exchange from the')
            log('  comment written into %s above, and moves'
                % os.path.basename(args.overrides))
            log('  every one of them to the other identifier by itself, so a '
                'ledger order would')
            log('  say the same thing twice, and its')
            log('  reject half would match no candidate of the run.  '
                '--legacy-orders writes them.')
        else:
            orders = legacy_orders(legacy, id1, id2, date,
                                   not args.rejects_only)
            added, already = merge_orders(led_rows, orders)
            if args.no_overrides:
                log('  --no-overrides: nothing else records the exchange, so '
                    'each identification')
                log('  needs an order of its own.')
            log('  %d order(s) added, %d already in the ledger.'
                % (len(added), len(already)))
            for row in added:
                log('    %-7s %-14s %s - %s'
                    % (row['decision'], row['wn_key'], row['low_id'],
                       row['upp_id']))

    if args.dry_run:
        log('')
        log('--dry-run: nothing was written.')
        return 0

    log('')
    if ov_rows is not None:
        save(args.overrides, ov_rows, ov_fields, not args.no_backup, log)
    if led_rows is not None:
        save(args.ledger, led_rows, led_fields, not args.no_backup, log)

    log('')
    log('Run classify_lines.py next.  Check afterwards that each of the two '
        'levels came out at the energy written here - lopt_vs_classify_levels'
        '.csv compares them - and that neither kept a line at its old '
        'position; a level fitted between two groups of lines is the sign '
        'that an identification was left behind.')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except ValueError as exc:
        sys.stderr.write('\nswap_line_assignments_pipeline.py: %s\n' % exc)
        sys.exit(1)
