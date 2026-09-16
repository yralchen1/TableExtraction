"""Move an already assigned energy level to a new position, from end to end.

WHAT THIS IS FOR
================
``insert_new_level.py`` puts into the pipeline a level that was never there.
This tool is for the other case: a level that IS there, whose lines are in the
fit and marked in IDEN2 and ruled on in the ledger, and whose position turns
out to be wrong.  ``level_positions.py --scan`` offers the alternate positions;
accepting one is a decision made in IDEN2, by eye; everything that has to
follow it is mechanical and touches the same six files in both directions at
once - the old position has to be taken apart before the new one can be built.

A level can be moved whichever list it came from:

* a level of the published list (Wyart's) is moved by a row of
  ``revised_level_energies.csv`` - the published workbook is an external input
  and is never edited;
* a level found since, one of ``new_levels.txt``, is moved by its own ``E``.

The tool works out which of the two the level is and writes the right file.

WHAT IT DOES, IN ORDER
======================
A. **Preflight.**  Every file the run can write is tested for writability -
   Excel locks the workbook it has open - and copied byte for byte into
   ``move_level_backup``, so that any later step can put everything back.
   ``--undo`` puts them back afterwards.

B. **The level and the two positions.**  The level is named by its identifier
   or by its row of ``enlev.dat``; the two are joined through
   ``IDEN2/IDEN_level_ids.txt`` and never by energy, because a level being
   moved is exactly the level whose energy has just changed.  The position it
   holds now is the pipeline's own adopted energy, overrides and all.  Where it
   is to go is ``--energy E``, or the position IDEN2 already holds for it
   (``--from-enlev``), or the one the current fit gives it (``--from-lopt``).

C. **What the level holds at the old position.**  Three files are read for
   everything that pins the level where it is: the records of
   ``LOPT_input_lines.txt`` that name it, the assignments marked for its block
   in ``IDEN2/trans.dat``, and the rows of ``line_decisions.csv`` that rule on
   its transitions.  A level that holds nothing anywhere has never been
   assigned and is ``insert_new_level.py``'s business, not this tool's.

D. **What survives the move.**  The whole model is rebuilt with the level at
   its new position and every one of those assignments is put to it again.  One
   that is still a candidate there survives and is left exactly as it is - the
   same record, the same mark, the same ledger row.  One that is not is
   RELEASED, and the reason is the test it failed::

       Too far from Ritz                 the observed wavenumber is more than
                                         --release-sigma times the line's own
                                         uncertainty from the Ritz wavenumber
                                         the new position gives
       Too weak to explain Iobs          the line is more than --strong-sigma
                                         stronger than this transition could
                                         make it
       Too weak to contribute to blend;  the line has other components and this
       may retain in pub list as masked. one's predicted share of it is below
                                         --masked-share
       No calculated transition at the   the pair is not in the calculated
       new position                      transition list at all

   Then the lines the new position opens up are proposed, on the rules
   ``insert_new_level.py`` proposes on: inside the Ritz window, never a line
   much stronger than predicted, never a component below ``--masked-share`` of
   the blend it would join, one line to a transition and one transition to a
   line.  Assignments already marked in IDEN2 at the new position are taken as
   given and never second-guessed - and, as there, a mark counts on the
   transition it was made on and not on every transition of this level that
   happens to reach the same observed line.  The ``mark`` column of the table
   is ``old`` for what IDEN2 shows already and ``new`` for what this run would
   add; ``--show-skipped`` lists the rows it passed over.  ``--reject`` and
   ``--accept`` take the observed wavenumber, with ``/PARTNER`` after it to
   name one component of a blend.

   Without ``--yes`` the run stops here, having written nothing, and prints both
   tables: what is released and why, and what is assigned.

E. **The old position is taken apart.**  The released records come out of
   ``LOPT_input_lines.txt``, the released marks out of ``IDEN2/trans.dat``, and
   the released ledger rows out of ``line_decisions.csv`` - the last kept, not
   thrown away, in ``line_decisions_removed.csv`` with the date and the reason,
   so that a move can be read back afterwards.

F. **The new position is built.**  The new energy goes into whichever of
   ``new_levels.txt`` and ``revised_level_energies.csv`` the level belongs to,
   one record per new assignment goes into the LOPT input, and the weights of
   EVERY observed wavenumber the move touched - the ones that lost a component
   as much as the ones that gained one - are divided again among their
   unflagged records in proportion to their calculated intensities.

G. **LOPT.**  ``lopt.bat LOPT.par`` (the Perl v5 build; never the Java jar,
   which has no centroid blend model), run once before anything is changed so
   that there is a before to compare with, and once after.
   ``RSS/degrees_of_freedom`` - the sum over every fitted line of its squared
   observed-minus-Ritz difference divided by its uncertainty, per degree of
   freedom, about 1 when the uncertainties are honest and the identifications
   right - is reported both ways, because a level moved at the cost of the rest
   of the fit has to show itself.  Then every line of the moved level: a
   residual above ``--ritz-sigma`` times its own uncertainty puts every file
   back and stops the run.  Otherwise the level's energy is replaced by the
   optimized one.

H. **classify_lines.py and the ledger.**  The classification is run and its
   verdicts on THIS LEVEL's transitions compared with what the move intends.
   One it will not make becomes an ``accept`` row; one it makes that the move
   does not intend becomes a ``reject`` row, with the reason chosen by the tests
   of D.  This repeats until nothing new has to be written, because a ledger row
   changes the classification.  ``make_LOPT_input.py`` IS NOT RUN: it would
   rebuild the fit's input out of the classification and so put into the fit
   every assignment the classification accepts, including ones nobody has
   looked at.  ``--rebuild`` asks for it explicitly.

I. **IDEN2.**  What the classification finally accepts for this level is what
   ``trans.dat`` is made to show: marks are added and released ones taken out.
   ANY OTHER LEVEL'S component on a line the move touched is named and left
   alone - releasing a line hands it back to whatever else can have it, and
   which of them gets it is a judgement to make on the screen, not one to make
   silently here.  The run is then HELD, with what it wrote left in place and
   ``--undo`` to put it back, until those are decided with ``--accept`` and
   ``--reject``.

J. **check_sync.py**, which must report no errors, then **sync_IDEN2.py**, which
   brings ``IDEN2/enlev.dat`` and the rest of IDEN2 into step with the fit.
   Then the report, on screen and in ``move_level.log``.

USAGE
=====
::

    python move_level.py 059003.000565 --energy 118967.3776
    python move_level.py --iden2-row 828 --from-lopt
    python move_level.py 059003.000565 --energy 118967.3776 --yes

The first form writes nothing: it prints the release and assignment tables and
stops.  ``--dry-run`` is the opposite kind of rehearsal - it runs the whole
sequence for real, programs and all, and then puts every file back, which is
the only way to see what the fit and the classification will say.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import datetime
import io
import math
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import insert_new_level as INL
import make_LOPT_input
import swap_line_assignments_IDEN as IDEN

from insert_new_level import Abort, Held, Log

# ---------------------------------------------------------------------------
# Where everything is
# ---------------------------------------------------------------------------
# The files the two tools share are insert_new_level.py's own, so that there
# can never be two spellings of one path.
ENLEV = INL.ENLEV
TRANS = INL.TRANS
ID_MAP = INL.ID_MAP
DLV = INL.DLV
NEW_LEVELS = INL.NEW_LEVELS
LINE_DECISIONS = INL.LINE_DECISIONS
CLASSIFICATIONS = INL.CLASSIFICATIONS
LOPT_INPUT = INL.LOPT_INPUT
LOPT_LINES = INL.LOPT_LINES
LOPT_LEVELS = INL.LOPT_LEVELS
SYNC_REPORT = INL.SYNC_REPORT

# The two this tool has of its own: the overrides file, which is how a level of
# the published list is moved, and the removed-rows file, where the ledger rows
# a move takes out are kept.
OVERRIDES = os.path.join(HERE, 'revised_level_energies.csv')
REMOVED = os.path.join(HERE, 'line_decisions_removed.csv')

LOGFILE = os.path.join(HERE, 'move_level.log')
BACKUP_DIR = os.path.join(HERE, 'move_level_backup')

WRITABLE = [p for p in INL.WRITABLE if p != INL.LOGFILE]
WRITABLE += [OVERRIDES, REMOVED, LOGFILE]

# ---------------------------------------------------------------------------
# The rules the release tests obey
# ---------------------------------------------------------------------------
# How many times its own uncertainty an observed wavenumber may be from the
# Ritz wavenumber the NEW position gives and still be kept.  Beyond it the line
# is not this transition's any more.
DEF_RELEASE_SIGMA = 3.0

# The smallest share of a blended feature's predicted intensity a component may
# have and still be kept on it.  Below this the component contributes nothing
# the fit can see; it is released, and the reason says it may still be worth
# printing in the published line list as a masked component.
DEF_MASKED_SHARE = 0.05

REASON_FAR = 'Too far from Ritz'
REASON_WEAK = 'Too weak to explain Iobs'
REASON_MASKED = ('Too weak to contribute to blend; may retain in pub list as '
                 'masked.')
REASON_NONE = 'No calculated transition at the new position'
REASON_ACCEPT = INL.REASON_ACCEPT


# ---------------------------------------------------------------------------
# B. The level, and where its energy lives
# ---------------------------------------------------------------------------
def resolve_level(args, log):
    """``(level_id, iden2_row)`` from whichever of the two was given.

    The join is always through ``IDEN2/IDEN_level_ids.txt`` and never by
    energy.  Every level the pipeline and LOPT use is in ``enlev.dat``, so a
    miss here means the lookup is broken, not that the level is absent, and the
    run says so rather than falling back to anything.
    """
    id_rows = INL.read_id_map_rows(ID_MAP)          # {level_id: row}
    by_row = {v: k for k, v in id_rows.items()}
    if args.level_id:
        level_id = args.level_id.strip()
        row = id_rows.get(level_id)
        if row is None:
            raise Abort('%s has no row for %s.  Every level the pipeline uses '
                        'is in enlev.dat, so this is a broken lookup, not a '
                        'missing level: check the table.'
                        % (os.path.basename(ID_MAP), level_id))
    else:
        row = args.iden2_row
        level_id = by_row.get(row)
        if level_id is None:
            raise Abort('%s gives no identifier to row %d'
                        % (os.path.basename(ID_MAP), row))
    log('   %s is row %d of enlev.dat' % (level_id, row))
    return level_id, row


def read_overrides(path):
    """``(fieldnames, rows)`` of ``revised_level_energies.csv``."""
    if not os.path.exists(path):
        return (['level_id', 'E_input', 'comment'], [])
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        rdr = csv.DictReader(fh)
        return (list(rdr.fieldnames or ['level_id', 'E_input', 'comment']),
                list(rdr))


def write_overrides(path, fields, rows):
    with io.open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fields, lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in fields})


def where_the_energy_lives(level_id):
    """``'new_levels'`` or ``'overrides'``: which file holds this level's
    adopted energy, and so which one a move has to write.

    A level that came from ``new_levels.txt`` carries its energy there and has
    no row in the published workbook at all, so an override row for it would be
    a second place to look and a second thing to forget.  Every other level's
    energy comes from the workbook, which is never edited, and is moved by an
    override row.
    """
    _fields, rows = INL.read_new_levels(NEW_LEVELS)
    for r in rows:
        if (r.get('level_id') or '').strip() == level_id:
            return 'new_levels'
    return 'overrides'


def set_energy(level_id, energy, where, comment, log, apply=True):
    """Put `energy` into whichever file holds this level's position."""
    if where == 'new_levels':
        fields, rows = INL.read_new_levels(NEW_LEVELS)
        for r in rows:
            if (r.get('level_id') or '').strip() != level_id:
                continue
            log('   %s: E %s -> %.4f cm^-1 in %s'
                % (level_id, r['E'], energy, os.path.basename(NEW_LEVELS)))
            r['E'] = '%.4f' % energy
            if comment:
                r['comment'] = comment
        if apply:
            INL.write_new_levels(NEW_LEVELS, fields, rows)
        return
    fields, rows = read_overrides(OVERRIDES)
    for name in ('level_id', 'E_input', 'comment'):
        if name not in fields:
            fields.append(name)
    for r in rows:
        if (r.get('level_id') or '').strip() == level_id:
            log('   %s: E_input %s -> %.4f cm^-1 in %s'
                % (level_id, r.get('E_input'), energy,
                   os.path.basename(OVERRIDES)))
            r['E_input'] = '%.4f' % energy
            if comment:
                r['comment'] = comment
            break
    else:
        log('   %s: a new row at %.4f cm^-1 in %s'
            % (level_id, energy, os.path.basename(OVERRIDES)))
        rows.append({'level_id': level_id, 'E_input': '%.4f' % energy,
                     'comment': comment or ''})
    if apply:
        write_overrides(OVERRIDES, fields, rows)


def refresh_energy(level_id, where, log):
    """Replace the level's energy by LOPT's optimized one, wherever it lives.

    The next classification then starts from where the fit put the level and
    not from the position it was moved to by hand, so its candidate windows and
    its Ritz tests are the fit's own.
    """
    optimized = INL.lopt_level_energy(LOPT_LEVELS, level_id)
    if optimized is None:
        log('   %s is not in %s; its energy is left as it was'
            % (level_id, os.path.basename(LOPT_LEVELS)))
        return None
    set_energy(level_id, optimized, where, None, log)
    return optimized


# ---------------------------------------------------------------------------
# C. What the level holds at the old position
# ---------------------------------------------------------------------------
class Holding(object):
    """One assignment the level holds at the position it is being moved from.

    The same assignment can be recorded in three places and this is all three
    of them in one object: the record of the LOPT input that puts it into the
    fit, the mark of ``trans.dat`` that shows it on the screen, and the rows of
    the ledger that rule on it.  A move that took it out of one and left it in
    another would be exactly the silent disagreement ``check_sync.py`` exists to
    find.
    """

    def __init__(self, wn, low_id, upp_id):
        self.wn = wn
        self.low_id = low_id
        self.upp_id = upp_id
        self.in_lopt = False
        self.partner_row = None     # the other level's row of enlev.dat
        self.ledger = []            # the rows of line_decisions.csv for it
        self.reason = ''            # why it is released; '' while it survives

    def key(self):
        return (round(self.wn, 3), self.low_id, self.upp_id)

    def where(self):
        bits = []
        if self.in_lopt:
            bits.append('LOPT')
        if self.partner_row is not None:
            bits.append('IDEN2')
        if self.ledger:
            bits.append('ledger x%d' % len(self.ledger))
        return '+'.join(bits) or '-'


def read_ledger_rows(path):
    """``(fieldnames, rows)`` of ``line_decisions.csv``."""
    if not os.path.exists(path):
        return (['wn_obs', 'low_id', 'upp_id', 'decision', 'date', 'reason'],
                [])
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        rdr = csv.DictReader(fh)
        return (list(rdr.fieldnames or []), list(rdr))


def write_ledger(path, fields, rows):
    with io.open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fields, lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in fields})


def archive_ledger_rows(path, rows, level_id, today):
    """Keep the ledger rows a move took out, rather than losing them.

    A ledger row is a verdict somebody reached on the screen.  The move makes
    it inapplicable - the transition it rules on no longer exists at the level's
    new position - but it is still the record of what was decided and why, and
    if the move is undone it is the thing to put back.
    """
    if not rows:
        return
    fields = ['removed_on', 'removed_for', 'removed_because',
              'wn_obs', 'low_id', 'upp_id', 'decision', 'date', 'reason']
    exists = os.path.exists(path)
    with io.open(path, 'a' if exists else 'w', encoding='utf-8',
                 newline='') as fh:
        w = csv.DictWriter(fh, fields, lineterminator='\n')
        if not exists:
            w.writeheader()
        for rec, why in rows:
            out = {k: rec.get(k, '') for k in fields}
            out['removed_on'] = today
            out['removed_for'] = level_id
            out['removed_because'] = why
            w.writerow(out)


def gather_holdings(level_id, index, levels_dict, log):
    """Everything that pins the level where it is, from all three files."""
    held = {}

    def get(wn, low_id, upp_id):
        key = (round(wn, 3), low_id, upp_id)
        if key not in held:
            held[key] = Holding(wn, low_id, upp_id)
        return held[key]

    lopt_rows = INL.read_lopt_input(LOPT_INPUT)
    for r in lopt_rows:
        if level_id in (r['low_id'], r['upp_id']):
            get(r['wn'], r['low_id'], r['upp_id']).in_lopt = True

    id_rows = INL.read_id_map_rows(ID_MAP)
    by_row = {v: k for k, v in id_rows.items()}
    trans = IDEN.Trans(TRANS)
    for partner_row, (wn, _dlv_row, _code) in \
            INL.marked_assignments(trans, index).items():
        other = by_row.get(partner_row)
        if other is None:
            log('   IDEN2 row %d, marked at %.3f cm^-1, has no identifier in '
                '%s; it cannot be moved with the level'
                % (partner_row, wn, os.path.basename(ID_MAP)))
            continue
        low_id, upp_id = order_pair(level_id, other, levels_dict)
        h = get(wn, low_id, upp_id)
        h.partner_row = partner_row

    for rec in read_ledger_rows(LINE_DECISIONS)[1]:
        low = (rec.get('low_id') or '').strip()
        upp = (rec.get('upp_id') or '').strip()
        if level_id not in (low, upp):
            continue
        try:
            wn = float(rec['wn_obs'])
        except (KeyError, TypeError, ValueError):
            continue
        get(wn, low, upp).ledger.append(rec)

    return [held[k] for k in sorted(held, reverse=True)], lopt_rows


def order_pair(a_id, b_id, levels_dict):
    """``(lower level, upper level)`` of two identifiers, by their energies."""
    a, b = levels_dict.get(a_id), levels_dict.get(b_id)
    if a is None or b is None:
        return (a_id, b_id)
    return (a_id, b_id) if a.energy <= b.energy else (b_id, a_id)


# ---------------------------------------------------------------------------
# D. What survives the move
# ---------------------------------------------------------------------------
def line_index(lines):
    """``(sorted wavenumbers, lines)`` for looking a line up by wavenumber."""
    ordered = sorted(lines, key=lambda l: l.wavenumber)
    return [l.wavenumber for l in ordered], ordered


def line_at(wns, ordered, wn, tol=5e-3):
    """The observed line at `wn`, or None."""
    j = bisect.bisect_left(wns, wn - tol)
    best, gap = None, tol
    while j < len(wns) and wns[j] <= wn + tol:
        d = abs(wns[j] - wn)
        if d <= gap:
            best, gap = ordered[j], d
        j += 1
    return best


def blend_share(wn, low_id, upp_id, lopt_rows, calc_of):
    """This component's share of the observed line's predicted intensity.

    A line several transitions share is fitted by LOPT as one blended feature,
    the components dividing it in proportion to their calculated intensities.  A
    component with a small enough share moves the feature's centroid by less
    than the fit can see, and keeping it on the line says more than the data
    support.  ``None`` when the line has no other component: a transition that
    is the whole of its line is never too weak for it.
    """
    group = [r for r in lopt_rows
             if round(r['wn'], 3) == round(wn, 3)
             and 'P' not in r['flag'].upper()]
    if len(group) < 2:
        return None
    total, mine = 0.0, 0.0
    for r in group:
        v = calc_of.get((r['low_id'], r['upp_id']), 0.0) or 0.0
        total += v
        if (r['low_id'], r['upp_id']) == (low_id, upp_id):
            mine = v
    if total <= 0:
        return None
    return mine / total


def release_reason(wn, low_id, upp_id, levels_dict, calc_index, wns, ordered,
                   lopt_rows, calc_of, args):
    """Why this assignment cannot be kept at the new position, or ``''``.

    The tests are in the order of how badly they fail: a line that is nowhere
    near the new Ritz wavenumber is not this transition's whatever its
    intensity, and a line whose intensity the transition cannot account for is
    not its however well the wavenumbers agree.
    """
    d = calc_index.get((low_id, upp_id))
    if d is None:
        return REASON_NONE
    low, upp = levels_dict.get(low_id), levels_dict.get(upp_id)
    if low is None or upp is None:
        return REASON_NONE
    rwn = upp.energy - low.energy
    line = line_at(wns, ordered, wn)
    if line is None:
        return REASON_NONE
    unc = getattr(line, 'wn_uncertainty', None) or 0.0
    resid = line.wavenumber - rwn
    if unc > 0 and abs(resid) > args.release_sigma * unc:
        return '%s (%+.3f cm^-1, %.1f sigma)' % (REASON_FAR, resid,
                                                 abs(resid) / unc)
    if abs(resid) > args.window:
        return '%s (%+.3f cm^-1)' % (REASON_FAR, resid)
    calc = d.get('calc_intensity')
    if calc and line.intensity > 0:
        u = math.sqrt((d.get('u_calc') or 1.0) ** 2 + INL.U_OBS_LN ** 2)
        z = math.log(line.intensity / calc) / u
        if z > args.strong_sigma:
            return '%s (%.0f times stronger than predicted, %.1f sigma)' \
                % (REASON_WEAK, line.intensity / calc, z)
    share = blend_share(wn, low_id, upp_id, lopt_rows, calc_of)
    if share is not None and share < args.masked_share:
        return '%s (%.1f%% of the blend)' % (REASON_MASKED, 100.0 * share)
    return ''


def holdings_table(holdings, log):
    log('   %-13s %-13s %11s %-16s %s'
        % ('lower', 'upper', 'obs wn', 'recorded in', 'what happens to it'))
    for h in holdings:
        log('   %-13s %-13s %11.3f %-16s %s'
            % (h.low_id, h.upp_id, h.wn, h.where(),
               h.reason or 'kept: still a candidate at the new position'))


# ---------------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description='Move an assigned level to a new position, and bring the '
                    'fit, the classification, the ledger and IDEN2 with it.')
    p.add_argument('level_id', nargs='?', default=None,
                   help='identifier of the level to move, 059003.000565')
    p.add_argument('--iden2-row', type=int, default=None, metavar='N',
                   help='name the level by its row of enlev.dat instead')
    g = p.add_mutually_exclusive_group()
    g.add_argument('--energy', type=float, default=None, metavar='E',
                   help='the new position, cm^-1')
    g.add_argument('--from-enlev', action='store_true',
                   help='take the new position from the measured energy '
                        'IDEN2 holds for the level in enlev.dat')
    g.add_argument('--from-lopt', action='store_true',
                   help='take the new position from the current fit, '
                        'LOPT_output_levels.txt')
    p.add_argument('--comment', default=None,
                   help='what to write in the comment column of the file that '
                        'records the move')
    p.add_argument('--reject', action='append',
                   metavar='WN[/PARTNER]=REASON', default=None,
                   help='never make this assignment; write a reject verdict '
                        'with this reason.  WN is the OBSERVED wavenumber, as '
                        'the obs wn column prints it, not the Ritz one; WN '
                        'alone names every transition of this level on that '
                        'line, WN/PARTNER (the other level, in full or by its '
                        'tail) names the one')
    p.add_argument('--accept', action='append',
                   metavar='WN[/PARTNER][=REASON]', default=None,
                   help='adopt the other components of the observed line at '
                        'WN that the classification accepts - the ones a '
                        'previous run listed as nobody\'s decision yet.  '
                        'WN/PARTNER names a single component')
    p.add_argument('--show-skipped', action='store_true',
                   help='list in the candidate table the rows the run passed '
                        'over, not only the ones it acts on')
    p.add_argument('--window', type=float, default=INL.DEF_WINDOW,
                   metavar='CM',
                   help='half-width of the Ritz window looked in, cm^-1 '
                        '(default %(default)s)')
    p.add_argument('--propose-window', type=float,
                   default=INL.DEF_PROPOSE_WINDOW, metavar='CM',
                   help='half-width of the window a line has to fall in '
                        'before it is proposed, cm^-1 (default %(default)s)')
    p.add_argument('--strong-sigma', type=float, default=INL.DEF_STRONG_SIGMA,
                   metavar='S',
                   help='how much stronger than predicted a line may be and '
                        'still be taken (default %(default)s)')
    p.add_argument('--release-sigma', type=float, default=DEF_RELEASE_SIGMA,
                   metavar='S',
                   help='how far from the new Ritz wavenumber, in the line\'s '
                        'own uncertainties, an existing assignment may be and '
                        'still be kept (default %(default)s)')
    p.add_argument('--masked-share', type=float, default=DEF_MASKED_SHARE,
                   metavar='F',
                   help='the smallest share of a blend a component may have '
                        'and still be kept on it (default %(default)s)')
    p.add_argument('--ritz-sigma', type=float, default=INL.DEF_RITZ_SIGMA,
                   metavar='S',
                   help='a residual above this many times the line\'s own '
                        'uncertainty in the fit stops the run '
                        '(default %(default)s)')
    p.add_argument('--no-propose', dest='propose', action='store_false',
                   default=True,
                   help='keep what survives the move and propose nothing new')
    p.add_argument('--max-rss', type=float, default=None, metavar='R',
                   help='stop if RSS/degrees_of_freedom comes out above this')
    p.add_argument('--rebuild', action='store_true',
                   help='run make_LOPT_input.py after the classification, so '
                        'that everything it accepts goes into the fit')
    p.add_argument('--undo', action='store_true',
                   help='put back everything the last run wrote and stop')
    p.add_argument('--yes', action='store_true',
                   help='apply; without it the run writes nothing')
    p.add_argument('--dry-run', action='store_true',
                   help='run the whole sequence for real and then put every '
                        'file back')
    p.add_argument('--no-sync', action='store_true',
                   help='skip check_sync.py and sync_IDEN2.py')
    args = p.parse_args(argv)
    if not args.undo and not args.level_id and args.iden2_row is None:
        p.error('name the level, either as an identifier or with --iden2-row')
    if not args.undo and args.energy is None and not args.from_enlev \
            and not args.from_lopt:
        p.error('say where the level goes: --energy E, --from-enlev or '
                '--from-lopt')
    return args


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------
def main(argv=None):
    args = parse_args(argv)
    rejects = INL.parse_reject(args.reject)
    accepts = INL.parse_accept(args.accept)
    log = Log()
    if args.undo:
        saved = {p: os.path.join(BACKUP_DIR, os.path.basename(p))
                 for p in WRITABLE
                 if os.path.exists(os.path.join(BACKUP_DIR,
                                                os.path.basename(p)))}
        if not saved:
            log('nothing to undo: %s holds no backup'
                % os.path.basename(BACKUP_DIR))
            return 1
        INL.restore(saved, log)
        return 0
    today = datetime.date.today().strftime('%-m/%-d/%Y'
                                           if os.name != 'nt'
                                           else '%#m/%#d/%Y')
    log('move_level.py - %s' % (args.level_id or 'IDEN2 row %d'
                                % args.iden2_row))
    log('')
    saved = INL.preflight(WRITABLE, log, BACKUP_DIR)
    try:
        return _run(args, rejects, accepts, saved, log, today)
    except Held as exc:
        log('')
        log('HELD: %s' % exc)
        if args.dry_run:
            log('')
            log('--dry-run: putting every file back')
            INL.restore(saved, log)
        else:
            log('   what this run wrote is left in place; %s --undo puts it '
                'back' % os.path.basename(__file__))
        return 2
    except Abort as exc:
        log('')
        log('STOPPED: %s' % exc)
        INL.restore(saved, log)
        return 1
    except BaseException as exc:
        log('')
        log('STOPPED by an unexpected %s: %s' % (type(exc).__name__, exc))
        INL.restore(saved, log)
        raise
    finally:
        log.save(LOGFILE)


def _run(args, rejects, accepts, saved, log, today):
    import classify_lines as CL

    # --- B. the level and the two positions ---------------------------------
    log('B. the level and the two positions')
    level_id, index = resolve_level(args, log)
    where = where_the_energy_lives(level_id)
    log('   its position is held in %s'
        % os.path.basename(NEW_LEVELS if where == 'new_levels' else OVERRIDES))

    levels_dict, _levels_list = CL.read_energy_levels()
    if level_id not in levels_dict:
        raise Abort('%s is not in the pipeline\'s level list.  A level that is '
                    'not there has nothing to move; insert_new_level.py is '
                    'what puts one in.' % level_id)
    e_old = levels_dict[level_id].energy

    if args.from_enlev:
        e_new, e_calc, _j, _label = INL.enlev_row(ENLEV, index)
        if e_new is None:
            raise Abort('enlev.dat row %d carries no measured energy (the '
                        'calculated position is %.3f cm^-1, which is not a '
                        'position to move to)' % (index, e_calc))
    elif args.from_lopt:
        e_new = INL.lopt_level_energy(LOPT_LEVELS, level_id)
        if e_new is None:
            raise Abort('%s has no row for %s'
                        % (os.path.basename(LOPT_LEVELS), level_id))
    else:
        e_new = args.energy
    log('   %.4f -> %.4f cm^-1   (dE = %+.4f)' % (e_old, e_new, e_new - e_old))
    if abs(e_new - e_old) < 5e-4:
        raise Abort('the level is already at %.4f cm^-1; there is nothing to '
                    'move' % e_new)

    # --- C. what it holds at the old position -------------------------------
    log('')
    log('C. what the level holds at the position it is leaving')
    holdings, lopt_rows = gather_holdings(level_id, index, levels_dict, log)
    if not holdings:
        raise Abort('%s holds no line in the fit, no mark in IDEN2 and no row '
                    'of the ledger.  It has never been assigned, so there is '
                    'nothing to move: insert_new_level.py is the tool for a '
                    'level being assigned for the first time.' % level_id)
    log('   %d assignment(s): %d in the LOPT input, %d marked in IDEN2, %d '
        'ruled on in the ledger'
        % (len(holdings),
           sum(1 for h in holdings if h.in_lopt),
           sum(1 for h in holdings if h.partner_row is not None),
           sum(1 for h in holdings if h.ledger)))

    # --- D. the model at the new position -----------------------------------
    log('')
    log('D. the level at its new position')
    levels_dict[level_id].energy = e_new
    calc_index = CL.read_transitions(levels_dict)
    lines = CL.read_observed_lines(levels_dict, calc_index)
    wns, ordered = line_index(lines)
    calc_of = {(k[0], k[1]): (v.get('calc_intensity') or 0.0)
               for k, v in calc_index.items()}
    calc_of.update(INL.classification_intensities(CLASSIFICATIONS))

    for h in holdings:
        h.reason = release_reason(h.wn, h.low_id, h.upp_id, levels_dict,
                                  calc_index, wns, ordered, lopt_rows,
                                  calc_of, args)
    for h in holdings:
        if not h.reason and any(s.matches(h.wn, h.low_id, h.upp_id)
                                for s in rejects):
            h.reason = 'rejected on the command line'
    kept = [h for h in holdings if not h.reason]
    log('   %d of the %d assignment(s) survive the move'
        % (len(kept), len(holdings)))
    holdings_table(holdings, log)

    # What IDEN2 still shows for this level after the released marks come out:
    # those are the decisions already taken, and they are what the proposals
    # are allowed to stand beside.
    kept_marked = {(h.low_id, h.upp_id): h.wn
                   for h in kept if h.partner_row is not None}

    level = levels_dict[level_id]
    candidates = INL.gather_candidates(level, levels_dict, calc_index, lines,
                                       args.window)
    released_keys = {h.key() for h in holdings if h.reason}
    free, held_back = [], []
    for c in candidates:
        if (round(c.wn, 3), c.low_id, c.upp_id) in released_keys:
            held_back.append(c)
        else:
            free.append(c)
    INL.choose(free, kept_marked, args.strong_sigma, args.propose_window,
               args.propose,
               share=INL.blend_share_of(free, calc_index, lopt_rows,
                                        kept_marked),
               min_share=args.masked_share)
    for c in held_back:
        c.source = ''
        c.verdict = 'released by the move'
    for c in free:
        if any(s.matches(c.wn, c.low_id, c.upp_id) for s in rejects):
            c.source = ''
            c.verdict = 'rejected on the command line'

    log('')
    log('   the lines the new position offers:')
    INL.candidate_table(sorted(candidates, key=lambda c: -c.rwn), log,
                        args.show_skipped)

    accepted = [c for c in free if c.source]
    log('')
    log('   %d assignment(s) at the new position: %d kept from IDEN2, %d '
        'proposed here' % (len(accepted),
                           sum(1 for c in accepted if c.source == 'IDEN2'),
                           sum(1 for c in accepted if c.source == 'proposed')))

    if not args.yes:
        log('')
        log('Nothing was written.  Run again with --yes to apply.')
        INL.restore(saved, log)
        return 0

    # --- G0. the fit as it stands, before anything is changed ---------------
    # The comparison a move has to answer for is against the fit as it was with
    # the level where it was.  It is taken before a single file is touched, so
    # that the RSS after the move has something to be measured against.
    log('')
    log('G0. the fit as it stands')
    rss_before, _dof = INL.run_lopt(log)

    # --- E. the old position is taken apart ---------------------------------
    log('')
    log('E. taking the old position apart')
    released = [h for h in holdings if h.reason]
    touched_wn = [h.wn for h in released] + [c.wn for c in accepted]

    ledger_fields, ledger_all = read_ledger_rows(LINE_DECISIONS)
    archived = []
    for h in released:
        for rec in h.ledger:
            archived.append((rec, h.reason))
    if archived:
        keep_rows = []
        dropped = 0
        want = {(round(float(r.get('wn_obs') or 0), 3),
                 (r.get('low_id') or '').strip(),
                 (r.get('upp_id') or '').strip()) for r, _w in archived}
        for rec in ledger_all:
            try:
                key = (round(float(rec['wn_obs']), 3),
                       (rec.get('low_id') or '').strip(),
                       (rec.get('upp_id') or '').strip())
            except (KeyError, TypeError, ValueError):
                keep_rows.append(rec)
                continue
            if key in want:
                dropped += 1
                continue
            keep_rows.append(rec)
        write_ledger(LINE_DECISIONS, ledger_fields, keep_rows)
        archive_ledger_rows(REMOVED, archived, level_id, today)
        log('   %d ledger row(s) taken out of %s and kept in %s'
            % (dropped, os.path.basename(LINE_DECISIONS),
               os.path.basename(REMOVED)))
    else:
        log('   no ledger row rules on a released assignment')

    trans = IDEN.Trans(TRANS)
    n_cleared = 0
    for h in released:
        if h.partner_row is None:
            continue
        k = trans.row(index, h.partner_row)
        if k is None:
            log('   %11.3f: no row of trans.dat to clear' % h.wn)
            continue
        trans.records[k] = trans.records[k][:IDEN.TR_OBS] + \
            INL.make_assignment(0, 0.0, 0.0, 0)
        log('   %11.3f  %s - %s  mark taken out of trans.dat'
            % (h.wn, h.low_id, h.upp_id))
        n_cleared += 1
    if n_cleared:
        IDEN.write_records(TRANS, trans.records, trans.endings)
    log('   %d mark(s) cleared in %s' % (n_cleared, os.path.basename(TRANS)))

    # --- F. the new position is built ---------------------------------------
    log('')
    log('F. building the new position')
    set_energy(level_id, e_new, where,
               args.comment or 're-positioned from %.4f (dE = %+.4f) by '
                               'move_level.py on %s' % (e_old, e_new - e_old,
                                                        today),
               log)

    lopt_rows = INL.read_lopt_input(LOPT_INPUT)
    gone = {(INL.lopt_key(h.wn), h.low_id, h.upp_id) for h in released}
    before_n = len(lopt_rows)
    lopt_rows = [r for r in lopt_rows
                 if (INL.lopt_key(r['wn']), r['low_id'], r['upp_id'])
                 not in gone]
    n_out = before_n - len(lopt_rows)
    # INL.lopt_key on both sides: the file holds three decimals and a
    # candidate the full precision of Pr3_lines.xlsx, so a finer rounding
    # makes a record that IS there look absent and inserts it twice.
    present = {(INL.lopt_key(r['wn']), r['low_id'], r['upp_id'])
               for r in lopt_rows}
    n_new = 0
    for c in accepted:
        key = (INL.lopt_key(c.wn), c.low_id, c.upp_id)
        if key in present:
            continue
        lopt_rows.append({'wn': c.wn, 'unc': c.line.wn_uncertainty,
                          'intens': c.line.intensity, 'low_id': c.low_id,
                          'upp_id': c.upp_id, 'flag': '', 'weight': 1.0,
                          'raw': None})
        present.add(key)
        n_new += 1
    n_w = INL.reweigh(lopt_rows, touched_wn, calc_of, log)
    INL.write_lopt_input(LOPT_INPUT, lopt_rows)
    log('   %d record(s) taken out, %d put in, %d weight(s) changed, %d '
        'record(s) in all' % (n_out, n_new, n_w, len(lopt_rows)))

    # --- G. LOPT ------------------------------------------------------------
    log('')
    log('G. LOPT')
    rss_after, _dof = INL.run_lopt(log)
    if rss_before is not None and rss_after is not None:
        log('   RSS/degrees_of_freedom %.2f -> %.2f (%+.2f)'
            % (rss_before, rss_after, rss_after - rss_before))
        if rss_after > rss_before:
            log('   the fit as a whole is worse than it was; the move is '
                'being paid for somewhere and every line of this level is '
                'worth looking at again')
    if args.max_rss is not None and rss_after is not None \
            and rss_after > args.max_rss:
        raise Abort('RSS/degrees_of_freedom came out at %.2f, above the %.2f '
                    'of --max-rss' % (rss_after, args.max_rss))
    pairs = {c.key() for c in accepted}
    intended = {(round(c.wn, 3), c.low_id, c.upp_id) for c in accepted}
    bad, worst, every, without = INL.ritz_culprits(LOPT_LINES, pairs,
                                                   args.ritz_sigma)
    log('   %d assignment(s) went into the fit; %d of them have a residual of '
        'their own, worst first:' % (len(accepted), len(every)))
    for low, upp, wn, d, z in every[:8]:
        log('     %11.3f  %s - %s  O-C = %+.3f cm^-1 (%.1f sigma)'
            % (wn, low, upp, d, z))
    if without:
        log('   %d share their line with another transition, so the fit gives '
            'the feature one residual and the component none:' % len(without))
        for low, upp, wn, flag in sorted(without, key=lambda t: -t[2]):
            log('     %11.3f  %s - %s  (%s)' % (wn, low, upp, flag or 'blended'))
    log('   the largest Ritz residual of the moved level: %.2f sigma' % worst)
    if bad:
        log('   ABNORMAL Ritz mismatches:')
        for low, upp, wn, d, z in sorted(bad, key=lambda t: -t[4]):
            log('     %11.3f  %s - %s  O-C = %+.3f cm^-1 (%.1f sigma)'
                % (wn, low, upp, d, z))
        raise Abort('%d line(s) of the moved level are further than %g sigma '
                    'from the fit' % (len(bad), args.ritz_sigma))

    refresh_energy(level_id, where, log)

    # --- H. classify_lines and the ledger -----------------------------------
    # Two kinds of disagreement, and a row for each.  An assignment the move
    # intends that the classification will not make needs an accept row, as it
    # does for a new level.  An assignment the classification makes that the
    # move does NOT intend needs a reject row - and that one is new here: a
    # level that has moved leaves behind transitions the classification is
    # still perfectly willing to propose at the old wavenumbers, and without a
    # verdict they would simply come back.  The reason is the test the
    # assignment failed in D, so the ledger says why and not merely that.
    log('')
    log('H. classify_lines.py and the ledger')
    reason_of = {h.key(): h.reason for h in holdings if h.reason}
    on_screen = set()
    if accepts:
        _tr = IDEN.Trans(TRANS)
        for (_owner, _partner), _k in _tr.row_of.items():
            _obs = IDEN.assignment(_tr.records[_k])
            if IDEN.has_line(_obs):
                on_screen.add((round(IDEN.obs_wavenumber(_obs), 2),
                               frozenset((_owner, _partner))))
        del _tr
    id_rows = INL.read_id_map_rows(ID_MAP)
    ledger_rows = []
    written = set()
    settled = False
    for attempt in range(1, 4):
        log('   round %d' % attempt)
        if INL.run([sys.executable, 'classify_lines.py'], log)[0] != 0:
            raise Abort('classify_lines.py failed')
        verdicts = INL.classification_verdicts(CLASSIFICATIONS)
        have = INL.read_ledger_keys(LINE_DECISIONS)
        fresh = []
        for c in accepted:
            key = (round(c.wn, 3), c.low_id, c.upp_id)
            if key in have or key in written or verdicts.get(key) is True:
                continue
            why = ('classify_lines did not propose it' if key not in verdicts
                   else 'classify_lines rejected it')
            fresh.append({'wn_obs': '%.4f' % c.wn, 'low_id': c.low_id,
                          'upp_id': c.upp_id, 'decision': 'accept',
                          'date': today, 'reason': REASON_ACCEPT})
            written.add(key)
            log('     accept %11.3f  %s - %s  (%s)'
                % (c.wn, c.low_id, c.upp_id, why))
        for (wn, low, upp), ok in sorted(verdicts.items(), reverse=True):
            if not ok or level_id not in (low, upp):
                continue
            key = (round(wn, 3), low, upp)
            if key in intended or key in have or key in written:
                continue
            reason = reason_of.get(key)
            if reason is None:
                reason = next((s.reason for s in rejects
                               if s.matches(wn, low, upp)), None)
            if reason is None:
                reason = ('not part of the move: this level was moved to '
                          '%.4f cm^-1 and this assignment was not carried '
                          'over' % e_new)
            fresh.append({'wn_obs': '%.4f' % wn, 'low_id': low,
                          'upp_id': upp, 'decision': 'reject',
                          'date': today, 'reason': reason})
            written.add(key)
            log('     reject %11.3f  %s - %s  (%s)' % (wn, low, upp, reason))
        if attempt == 1:
            for spec in accepts:
                wn, reason = spec.wn, spec.reason
                hit = 0
                for (kwn, low, upp), ok in sorted(verdicts.items()):
                    if not ok or level_id in (low, upp)                             or not spec.matches(kwn, low, upp):
                        continue
                    hit += 1
                    key = (round(kwn, 3), low, upp)
                    if key in have or key in written:
                        continue
                    rows = frozenset((id_rows.get(low), id_rows.get(upp)))
                    if (round(kwn, 2), rows) in on_screen:
                        log('     %11.3f  %s - %s is on the screen already; '
                            'left as it is' % (kwn, low, upp))
                        continue
                    fresh.append({'wn_obs': '%.4f' % kwn, 'low_id': low,
                                  'upp_id': upp, 'decision': 'accept',
                                  'date': today, 'reason': reason})
                    written.add(key)
                    log('     adopt  %11.3f  %s - %s  (%s)'
                        % (kwn, low, upp, reason))
                if not hit:
                    log('     --accept %.3f matches no assignment the '
                        'classification accepts on that line' % wn)
        if fresh:
            INL.append_ledger(LINE_DECISIONS, fresh)
            ledger_rows.extend(fresh)
            log('     %d row(s) added to %s'
                % (len(fresh), os.path.basename(LINE_DECISIONS)))
        else:
            log('     the ledger already says everything this run needs it to '
                'say')

        if args.rebuild:
            if INL.run([sys.executable, 'make_LOPT_input.py'], log)[0] != 0:
                raise Abort('make_LOPT_input.py failed')
            rss_r, _dof = INL.run_lopt(log)
            if rss_before is not None and rss_r is not None:
                log('     RSS/degrees_of_freedom %.2f -> %.2f (%+.2f)'
                    % (rss_before, rss_r, rss_r - rss_before))
            if args.max_rss is not None and rss_r is not None \
                    and rss_r > args.max_rss:
                raise Abort('RSS/degrees_of_freedom came out at %.2f after '
                            'the rebuild, above the %.2f of --max-rss'
                            % (rss_r, args.max_rss))
            bad, worst, _every, _w = INL.ritz_culprits(LOPT_LINES, pairs,
                                                       args.ritz_sigma)
            log('     the largest Ritz residual after the rebuild: %.2f sigma'
                % worst)
            if bad:
                for low, upp, wn, d, z in sorted(bad, key=lambda t: -t[4]):
                    log('       %11.3f  %s - %s  O-C = %+.3f cm^-1 '
                        '(%.1f sigma)' % (wn, low, upp, d, z))
                raise Abort('%d line(s) are further than %g sigma from the '
                            'rebuilt fit' % (len(bad), args.ritz_sigma))
            refresh_energy(level_id, where, log)
        if not fresh:
            settled = True
            break
    if not settled:
        raise Abort('the classification and the ledger did not settle in three '
                    'rounds; look at the assignments by hand')

    # --- I. IDEN2 -----------------------------------------------------------
    log('')
    log('I. IDEN2')
    verdicts_final = INL.classification_verdicts(CLASSIFICATIONS)
    adopt_wn = {round(s.wn, 2) for s in accepts}
    final = {k: v for k, v in verdicts_final.items()
             if level_id in (k[1], k[2])
             or (v and round(k[0], 2) in adopt_wn)}
    trans = IDEN.Trans(TRANS)
    marked = INL.marked_assignments(trans, index)
    marked_wn = {round(wn, 2) for wn, _row, _code in marked.values()}
    orphan = [wn for wn in marked_wn
              if not any(round(k[0], 2) == wn and ok and level_id in (k[1], k[2])
                         for k, ok in final.items())]
    if orphan:
        raise Abort('%d line(s) still marked in IDEN2 for row %d are not '
                    'accepted by the classification even with the ledger rows '
                    'this run wrote: %s.  Decide them by hand before running '
                    'this again.'
                    % (len(orphan), index,
                       ', '.join('%.3f' % w for w in sorted(orphan))))

    dlv = INL.read_dlv(DLV)
    by_wn = sorted(dlv.items(), key=lambda kv: kv[1][1])
    dlv_wns = [v[1] for _k, v in by_wn]
    n_marked, n_already, n_missed = 0, 0, 0
    for (wn, low_id, upp_id), ok in sorted(final.items(), reverse=True):
        if not ok:
            continue
        row_low, row_upp = id_rows.get(low_id), id_rows.get(upp_id)
        if row_low is None or row_upp is None:
            log('   %11.3f  %s - %s: a level has no IDEN2 row; not marked'
                % (wn, low_id, upp_id))
            n_missed += 1
            continue
        k = trans.row(row_upp, row_low)
        if k is None:
            log('   %11.3f  %s - %s: no row in trans.dat; not marked'
                % (wn, low_id, upp_id))
            n_missed += 1
            continue
        obs = IDEN.assignment(trans.records[k])
        if IDEN.has_line(obs) and abs(IDEN.obs_wavenumber(obs) - wn) < 0.02:
            n_already += 1
            continue
        j = bisect.bisect_left(dlv_wns, wn)
        hit = None
        for m in (j - 1, j):
            if 0 <= m < len(by_wn) and abs(dlv_wns[m] - wn) < 0.02:
                hit = by_wn[m]
        if hit is None:
            log('   %11.3f: no row of dlv.dat; not marked' % wn)
            n_missed += 1
            continue
        dlv_row, (code, dlv_wn) = hit
        rwn = levels_dict[upp_id].energy - levels_dict[low_id].energy
        trans.records[k] = trans.records[k][:IDEN.TR_OBS] + \
            INL.make_assignment(code, dlv_wn, dlv_wn - rwn, dlv_row)
        log('   %11.3f  %s - %s  marked (IDEN2 rows %d - %d)'
            % (dlv_wn, low_id, upp_id, row_low, row_upp))
        n_marked += 1
    if n_marked:
        IDEN.write_records(TRANS, trans.records, trans.endings)
    log('   %d assignment(s) written into %s, %d were already there, %d could '
        'not be placed'
        % (n_marked, os.path.basename(TRANS), n_already, n_missed))

    # Every other level's component on a line this move touched that the
    # classification now accepts and IDEN2 does not show.  A released line is
    # handed back to whatever else can have it, and which of them gets it is
    # the judgement this tool does not make.
    touched = {round(w, 2) for w in touched_wn}
    marked_all = set()
    for (owner, partner), k in trans.row_of.items():
        obs = IDEN.assignment(trans.records[k])
        if IDEN.has_line(obs):
            marked_all.add((round(IDEN.obs_wavenumber(obs), 2), owner, partner))
    unchecked = []
    for (wn, low_id, upp_id), ok in sorted(verdicts_final.items(),
                                           reverse=True):
        if not ok or level_id in (low_id, upp_id):
            continue
        if round(wn, 2) not in touched:
            continue
        row_low, row_upp = id_rows.get(low_id), id_rows.get(upp_id)
        if row_low is None or row_upp is None:
            continue
        if (round(wn, 2), row_upp, row_low) in marked_all \
                or (round(wn, 2), row_low, row_upp) in marked_all:
            continue
        unchecked.append((wn, low_id, upp_id, row_low, row_upp))
    if unchecked:
        log('')
        log('   %d assignment(s) of other levels on the lines this move '
            'touched are now accepted by the classification and are not '
            'marked in IDEN2:' % len(unchecked))
        now_rows = INL.classification_rows(CLASSIFICATIONS)
        before_all = INL.classification_rows(saved.get(CLASSIFICATIONS)
                                             or CLASSIFICATIONS)
        for wn, low_id, upp_id, row_low, row_upp in unchecked:
            log('')
            log('     %11.3f  %s - %s  (IDEN2 rows %d - %d)'
                % (wn, low_id, upp_id, row_low, row_upp))
            key = round(wn, 2)
            was_rows = {((r.get('low_id') or '').strip(),
                         (r.get('upp_id') or '').strip()): r
                        for r in before_all.get(key, [])}
            INL.line_dossier(wn, now_rows.get(key, []), was_rows, level_id, log)
        log('')
        log('   They are not this level\'s and this move did not propose them, '
            'so it has not written them.  A line this move released is free '
            'again and something else may now be able to have it - which is '
            'exactly the judgement to make on the screen.  Look at each one; '
            'then run this tool again, naming the lines you found good with '
            '--accept and the rest with --reject:')
        log('')
        log('     python %s %s --energy %.4f --rebuild --yes%s'
            % (os.path.basename(__file__), level_id, e_new,
               ''.join(' --accept %.3f' % wn
                       for wn in sorted({round(u[0], 3) for u in unchecked}))))

    # --- J. check_sync and sync_IDEN2 ---------------------------------------
    if args.no_sync:
        log('')
        log('J. skipped (--no-sync)')
    else:
        log('')
        log('J. check_sync.py')
        status, _out = INL.run([sys.executable, 'check_sync.py'], log)
        if status > 1:
            log('   the ERROR findings of %s:' % os.path.basename(SYNC_REPORT))
            shutil.copy2(SYNC_REPORT,
                         os.path.join(BACKUP_DIR, 'sync_report.failed.txt'))
            for rec in io.open(SYNC_REPORT, encoding='utf-8',
                               errors='replace'):
                if rec.strip().startswith('ERROR'):
                    log('     ' + rec.rstrip())
            if unchecked:
                raise Held('check_sync.py reports the %d assignment(s) above '
                           'as findings.  Decide them, then run this tool '
                           'again with --accept/--reject as printed above.'
                           % len(unchecked))
            raise Abort('check_sync.py reported errors')
        log('   no errors%s' % (' (warnings only)' if status == 1 else ''))

        if unchecked:
            raise Held('%d assignment(s) on the lines this move touched are '
                       'nobody\'s decision yet; sync_IDEN2.py is not run until '
                       'they are made.' % len(unchecked))

        log('')
        log('K. sync_IDEN2.py')
        if INL.run([sys.executable, 'sync_IDEN2.py'], log)[0] != 0:
            raise Abort('sync_IDEN2.py failed')

    log('')
    log('L. done')
    log('   %s moved %.4f -> %.4f cm^-1 (dE = %+.4f): %d assignment(s) '
        'released, %d kept, %d proposed, %d ledger row(s) written, %d taken '
        'out' % (level_id, e_old, e_new, e_new - e_old, len(released),
                 sum(1 for c in accepted if c.source == 'IDEN2'),
                 sum(1 for c in accepted if c.source == 'proposed'),
                 len(ledger_rows), len(archived)))
    if args.dry_run:
        log('')
        log('--dry-run: putting every file back')
        INL.restore(saved, log)
    return 0


if __name__ == '__main__':
    sys.exit(main())
