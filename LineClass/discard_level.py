"""Give up a level's measured position, from end to end.

WHAT THIS IS FOR
================
``level_positions.py --audit`` sometimes reports a level whose present
position is not merely second best but worse than nothing: the observed lines
are better explained with no level there at all.  Level 812 is one - two lines,
both blends, neither of them its own - and level 436 another, its four lines
wanting four different energies spread over 1.9 cm^-1.  Neither has an
alternate position the audit supports, so there is nowhere to move them to.
This tool carries out the decision to stop claiming the position.

WHAT IS BEING DENIED IS THE POSITION, NOT THE LEVEL
===================================================
Cowan's calculation predicts both levels and will go on predicting them.  Of
the 1253 calculated levels of Pr III, 636 have been found and 617 have never
been placed, and what decides which pool a level is in is a single character -
the ``*`` in columns 39-40 of ``IDEN2/enlev.dat``.  So this tool does not
delete a level.  It UNFINDS it: puts it back among the levels nobody has
found, where ``unfound_levels.py`` ranks it by ``n_prom`` and
``level_positions.py --unknown`` can search for it again.

That choice is what keeps the change reversible and the statistics honest.  A
deleted level would quietly remove a data point from the per-configuration
window - the rms of ``E_obs - E_calc`` over the found levels of a
configuration, which sets how wide every future search of that configuration
scans - and the found and unfound counts would stop adding to 1253.

The recomputation of that window is the part of a discard most easily
overlooked, and it is a gain rather than a cost.  A level held at a position
the evidence does not support is usually one of the worst-placed levels of its
configuration, and while it counts as found it is inflating the search window
of every OTHER level of that configuration.  Taking it out of the found set
narrows the window for all of them: dropping 812 takes ``f26d`` from 60.4 to
57.6 cm^-1 for its 61 unfound levels, and dropping 436 takes ``f25f`` from
106.1 to 103.4 for its 15.  ``sync_IDEN2.configuration_uncertainties`` does the
arithmetic, and step J runs it, so this tool has only to clear the star.

CLEARING THE STAR IS NOT ENOUGH ON ITS OWN
==========================================
The pipeline does not read IDEN2's star when it builds its level list; it
reads Wyart's published workbook, which is an external input and is never
edited.  A level unfound in IDEN2 alone would be generated with its full set of
calculated transitions on the next ``classify_lines.py`` run and would simply
re-acquire lines - and a published identification naming it is protected by
Step 3's "old, no solid evidence for rejection" rule and would come straight
back.  Hence ``discarded_levels.csv``, which the pipeline's own readers obey.

WHAT IT DOES, IN ORDER
======================
A. **Preflight.**  Every file the run can write is tested for writability -
   Excel locks the workbook it has open - and copied byte for byte into
   ``discard_level_backup``.  ``--undo`` puts them back.

B. **The level**, named by its identifier or by its row of ``enlev.dat``, the
   two joined through ``IDEN2/IDEN_level_ids.txt`` and never by energy.

C. **What it holds**: the records of ``LOPT_input_lines.txt`` that name it, the
   assignments marked for its block in ``IDEN2/trans.dat``, and the rows of
   ``line_decisions.csv`` that rule on its transitions.

D. **Two guards, before anything is written.**

   *The no-op guard*, at the top of the run so that starting the last run again
   by mistake - the commonest way of losing several minutes to PyCharm's re-run
   icon - costs nothing: a level already absent from the level list, with no
   map row, no holdings anywhere and a row already in ``discarded_levels.csv``,
   has been discarded already.

   *The evidence guard*: the level's row of the audit table is read, and the
   run refuses unless ``ln_R`` is negative and no alternate position of the
   level is firm.  A level with real support that merely sits in the wrong
   place is ``move_level.py``'s business, and the two tools must not become
   interchangeable by accident.  ``--force`` overrides and says so in the log
   and in the ledger.

   Without ``--yes`` the run stops here having written nothing, printing the
   release table: every held line, its observed intensity, this level's share
   of it and what else claims it.

E. **The old position is taken apart.**  Every holding is released - there is
   no surviving set to work out, which is the one place a discard is simpler
   than a move.  The records come out of the LOPT input; the weights of every
   observed wavenumber that lost a component are divided again among its
   unflagged records in proportion to their calculated intensities; the marks
   come out of ``trans.dat``; and the ledger rows are moved, not thrown away,
   into ``line_decisions_removed.csv`` with the date and the reason.

F. **The level leaves the list.**  The ``discarded_levels.csv`` row is written,
   the ``IDEN_level_ids.txt`` row removed - which is what makes the rest of the
   change cost nothing, see ``remove_id_map_row`` - and any ``new_levels.txt``
   or ``revised_level_energies.csv`` row of the level removed with it.

G. **LOPT**, and the Ritz check INVERTED.  A discard leaves the level with no
   line of its own to test, so the test moves to the blend partners: every
   observed wavenumber that lost a component is re-examined, and a SURVIVING
   component whose residual now exceeds ``--ritz-sigma`` times its own
   uncertainty puts every file back and stops the run.  A feature whose two
   Ritz values straddle the observed wavenumber is a real blend, and taking one
   component away leaves the other with a residual it cannot explain.  It is
   the one way a discard can damage the fit, and the run refuses rather than
   reporting it afterward.

H. **classify_lines.py and the ledger.**  A discard writes NO ledger row for
   the level, and must not: ``check_forced_decisions`` treats a row naming a
   level that is not in the level list as a fault and stops the run, for a
   ``reject`` row as much as an ``accept`` one.  A move needs reject rows
   because the level is still in the list and its old lines would be proposed
   again; a discarded level generates no candidate at all, so there is nothing
   for a verdict to rule on, and what was decided is kept in
   ``line_decisions_removed.csv``.  What this step does instead is check: the
   classification must propose nothing for the level, and no ledger row may
   still name it.  Anything it now accepts for ANOTHER level on a freed line is
   reported and left alone.  ``make_LOPT_input.py`` IS NOT RUN - it would
   rebuild the fit's input out of the classification and put back every
   assignment the classification accepts, including ones nobody has looked at.
   ``--rebuild`` asks for it explicitly.

I. **IDEN2.**  The level's ``enlev.dat`` row and its ``trans.dat`` block are
   returned to the unfound shape, and every copy of its energy that ``trans.dat``
   keeps in other levels' blocks is brought with them.  Any other level's
   component on a freed line is NAMED AND LEFT ALONE: releasing a line hands it
   back to whatever else can have it, and which of them gets it is a judgment
   made on the screen.  The run is then HELD until those are settled.

J. **check_sync.py**, which must exit no worse than warnings, then
   **sync_IDEN2.py**, whose ``configuration_uncertainties`` step gives the
   discarded level and every other unfound level of its configuration the
   recomputed search window.  Then the report, on screen and in
   ``discard_level.log``.

USAGE
=====
::

    python discard_level.py --iden2-row 812 --reason "no support; 2 lines"
    python discard_level.py 059003.000561 --reason "..." --yes
    python discard_level.py --iden2-row 812 --undo

The first form writes nothing.  ``--dry-run`` is the opposite kind of
rehearsal: it runs the whole sequence for real, programs and all, and then puts
every file back, which is the only way to see what the fit and the
classification will say.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import io
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import insert_new_level as INL
import level_interchange
import move_level as ML
import swap_line_assignments_IDEN as IDEN

from insert_new_level import Abort, Held, Log

# ---------------------------------------------------------------------------
# Where everything is
# ---------------------------------------------------------------------------
# Shared paths come from insert_new_level.py through move_level.py, so that no
# path can ever have two spellings.
ENLEV = INL.ENLEV
TRANS = INL.TRANS
ID_MAP = INL.ID_MAP
NEW_LEVELS = INL.NEW_LEVELS
LINE_DECISIONS = INL.LINE_DECISIONS
CLASSIFICATIONS = INL.CLASSIFICATIONS
LOPT_INPUT = INL.LOPT_INPUT
LOPT_LINES = INL.LOPT_LINES
LOPT_LEVELS = INL.LOPT_LEVELS
SYNC_REPORT = INL.SYNC_REPORT
OVERRIDES = ML.OVERRIDES
REMOVED = ML.REMOVED

DISCARDED = os.path.join(HERE, 'discarded_levels.csv')
AUDIT = os.path.join(HERE, 'level_positions.csv')

LOGFILE = os.path.join(HERE, 'discard_level.log')
BACKUP_DIR = os.path.join(HERE, 'discard_level_backup')

WRITABLE = [p for p in ML.WRITABLE if p != ML.LOGFILE]
WRITABLE += [DISCARDED, LOGFILE]

DISCARD_COLUMNS = ['level_id', 'iden2_row', 'date', 'ln_R', 'reason']

# The reason written into every ledger row a discard makes.  It names the tool
# rather than a test, because unlike a move nothing here failed a test: the
# assignment is released because the level it belonged to has no position any
# more.
REASON_DISCARDED = 'the level has been discarded; it has no position'


# ---------------------------------------------------------------------------
# The ledger of discards
# ---------------------------------------------------------------------------
def read_discarded(path: str) -> tuple:
    """``(fieldnames, rows)`` of ``discarded_levels.csv``."""
    if not os.path.exists(path):
        return (list(DISCARD_COLUMNS), [])
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        rdr = csv.DictReader(fh)
        return (list(rdr.fieldnames or DISCARD_COLUMNS), list(rdr))


def write_discarded(path: str, fields, rows) -> None:
    for name in DISCARD_COLUMNS:
        if name not in fields:
            fields.append(name)
    with io.open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fields, lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in fields})


def discarded_row(path: str, level_id: str):
    """The row of the ledger for this level, or ``None``."""
    for r in read_discarded(path)[1]:
        if (r.get('level_id') or '').strip() == level_id:
            return r
    return None


# ---------------------------------------------------------------------------
# F. The identifier map
# ---------------------------------------------------------------------------
def remove_id_map_row(path: str, level_id: str, log) -> int:
    """Take this level's line out of ``IDEN_level_ids.txt``; return its row.

    THIS IS WHAT MAKES THE REST OF THE CHANGE COST NOTHING, and it is worth
    being explicit about why.  Two programs would otherwise have to learn what
    a discarded level is:

    * ``sync_IDEN2.adopted_energies()`` raises for a level the map names and
      the LOPT output does not - which a discarded level is, its records having
      left the fit.  With the map row gone the level is simply not starred, and
      that function's own existing branch leaves it exactly as it stands.
    * ``check_sync._iden_vs_fit()`` raises an ERROR for every mapped level
      whose ``enlev.dat`` row carries no star, and lists it again as being in
      the map and not in ``LOPT_output_levels``.  Both vanish with the row.

    The association is not lost: it is in the ``iden2_row`` column of
    ``discarded_levels.csv``, which is where ``--undo`` reads it back from.
    """
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        text = fh.read()
    lines = text.replace('\r\n', '\n').rstrip('\n').split('\n')
    head, body = lines[0], lines[1:]
    kept, row = [], None
    for rec in body:
        if rec.split('\t')[0].strip() == level_id:
            try:
                row = int(rec.split('\t')[1])
            except (IndexError, ValueError):
                row = None
            continue
        kept.append(rec)
    if row is None:
        log('   %s has no row for %s; nothing to take out'
            % (os.path.basename(path), level_id))
        return 0
    with io.open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join([head] + kept) + '\n')
    log('   %s -> row %d taken out of %s'
        % (level_id, row, os.path.basename(path)))
    return row


# ---------------------------------------------------------------------------
# D. The evidence guard
# ---------------------------------------------------------------------------
def audit_row(path: str, level_id: str):
    """The level's row of ``level_positions.csv``, or ``None``."""
    if not os.path.exists(path):
        return None
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        for rec in csv.DictReader(fh):
            if (rec.get('level_id') or '').strip() == level_id:
                return rec
    return None


def _float(rec, key):
    try:
        return float((rec.get(key) or '').strip())
    except (AttributeError, TypeError, ValueError):
        return None


def check_the_evidence(level_id, args, log):
    """Refuse a discard the audit does not support; return ``ln_R``.

    Two tests, and both are about keeping this tool and ``move_level.py``
    distinct.  ``ln_R`` is the log of how much more likely the observed lines
    are with the level at its adopted energy than with no level there at all:
    a discard is only honest when it is negative, which is to say when the
    lines are better explained without the level.  And a level with a firm
    alternate position has somewhere to go, so what it wants is a move.

    A level with no row in the audit table is refused too.  The audit is the
    evidence, and running this tool on a level nobody has audited would be
    acting on no evidence whatever.
    """
    rec = audit_row(args.audit, level_id)
    if rec is None:
        if args.force:
            log('   %s has no row in %s; --force goes on anyway'
                % (level_id, os.path.basename(args.audit)))
            return None
        raise Abort('%s has no row in %s.  That table is the evidence a '
                    'discard rests on, so run  python level_positions.py '
                    '--audit  first.' % (level_id,
                                         os.path.basename(args.audit)))
    ln_r = _float(rec, 'ln_R')
    action = (rec.get('action') or '').strip()
    n_alt = _float(rec, 'n_alt')
    log('   the audit says: ln_R = %s, %s alternate position(s), action %r'
        % ('%.3f' % ln_r if ln_r is not None else '?',
           '%d' % n_alt if n_alt is not None else '?', action or '-'))
    if ln_r is not None and ln_r >= 0 and not args.force:
        raise Abort('%s has ln_R = %+.3f: the observed lines are MORE likely '
                    'with the level at this position than without it, so the '
                    'position has support and this is not a level to discard.  '
                    'If it is in the wrong place, move_level.py is the tool.  '
                    '--force overrides.' % (level_id, ln_r))
    if action in ('firm', 'relocate') and not args.force:
        raise Abort('the audit calls %s\'s best alternate position %r, so the '
                    'level has somewhere to go and wants move_level.py, not '
                    'this tool.  --force overrides.' % (level_id, action))
    return ln_r


# ---------------------------------------------------------------------------
# D. The release table
# ---------------------------------------------------------------------------
def mark_fitted(holdings, lopt_rows):
    """Say which holdings are IN the fit and which are only printed beside it.

    A record of ``LOPT_input_lines.txt`` flagged ``P`` is one the
    classification REJECTED: make_LOPT_input.py writes it with weight 0 so
    that LOPT prints a predicted wavenumber for the pair, and the fit takes no
    notice of it.  Level 812 has twenty-two records and two weights, which is
    why the audit counts two lines for it and not twenty-two.

    The distinction runs through the whole run.  Only a fitted record can move
    a level, so only a fitted record's wavenumber needs its weights divided
    again, and only a fitted record's blend partners can be left with a
    residual they cannot explain.  Every record goes out of the file either
    way: the level is leaving the level list, and make_LOPT_input.py would
    write none of them at the next rebuild.

    Sets ``fitted`` and ``weight`` on each holding; returns how many are fitted.
    """
    weight_of = {}
    for r in lopt_rows:
        key = (INL.lopt_key(r['wn']), r['low_id'], r['upp_id'])
        if 'P' in (r['flag'] or '').upper():
            weight_of.setdefault(key, None)
        else:
            weight_of[key] = r['weight']
    n = 0
    for h in holdings:
        w = weight_of.get((INL.lopt_key(h.wn), h.low_id, h.upp_id))
        h.weight = w
        h.fitted = w is not None
        n += 1 if h.fitted else 0
    return n


def release_table(holdings, lopt_rows, calc_of, lines_at, log):
    """Every line the level holds, and what the feature looks like without it.

    This is the table the analyst reads before saying ``--yes``.  For a level
    like 812, whose two fitted lines are both blends it owns barely half of,
    what matters is not that the lines are lost but who else is on them: a line
    the level owns outright goes back to being unidentified, and a line it
    shares goes back to its other components, whose residuals step G then has
    to re-examine.

    The rows the fit ignores are shown too, and marked, rather than left out:
    they are going out of the file with the rest, and a table that showed only
    two of twenty-two records would leave the other twenty to be discovered in
    the diff.
    """
    log('   %-13s %-13s %11s %9s %9s %7s %-16s %s'
        % ('lower', 'upper', 'obs wn', 'I obs', 'weight', 'share',
           'recorded in', 'what else claims the line'))
    for h in holdings:
        share = (ML.blend_share(h.wn, h.low_id, h.upp_id, lopt_rows, calc_of)
                 if h.fitted else None)
        others = [r for r in lopt_rows
                  if round(r['wn'], 3) == round(h.wn, 3)
                  and (r['low_id'], r['upp_id']) != (h.low_id, h.upp_id)
                  and 'P' not in (r['flag'] or '').upper()]
        line = lines_at(h.wn)
        log('   %-13s %-13s %11.3f %9s %9s %7s %-16s %s'
            % (h.low_id, h.upp_id, h.wn,
               '%.2f' % line.intensity if line is not None else '-',
               '%.4f' % h.weight if h.fitted else 'rejected',
               ('%.0f%%' % (100.0 * share) if share is not None
                else ('all' if h.fitted else '-')),
               h.where(),
               ', '.join('%s - %s' % (r['low_id'], r['upp_id'])
                         for r in others)
               or ('nothing: the line is its own' if h.fitted
                   else 'nothing the fit uses')))


# ---------------------------------------------------------------------------
# G. The inverted Ritz check
# ---------------------------------------------------------------------------
def partner_culprits(path, wavenumbers, level_id, sigma):
    """Surviving components of the freed lines that the fit can no longer place.

    The test a move makes on the moved level's own lines has nothing to work on
    here - a discarded level has no lines left - so it is turned round and made
    on the lines the discard touched.  A component that was sharing its feature
    with the discarded level's transition is now alone on it, and if the
    observed wavenumber was really the centroid of a blend, the survivor is
    left with a residual it cannot explain.  That is the 42117.349 case: a
    feature whose two Ritz values straddle the observed wavenumber is a real
    blend, and taking one component away is not free.

    Returns ``(culprits, worst, every)`` in the shape ``INL.ritz_culprits``
    returns them, over every row of ``LOPT_output_lines.txt`` at one of
    `wavenumbers` that does not name `level_id`.
    """
    want = {INL.lopt_key(w) for w in wavenumbers}
    bad, worst, every = [], 0.0, []
    for r in INL.lopt_lines(path):
        low, upp = r.get('L1', ''), r.get('L2', '')
        if level_id in (low, upp):
            continue
        try:
            wn = float(r['wn_o'])
        except (KeyError, ValueError):
            continue
        if INL.lopt_key(wn) not in want:
            continue
        try:
            d = float(r['dWO-C'])
            u = float(r['uWnOTot'])
        except (KeyError, ValueError):
            continue
        if u <= 0:
            continue
        z = abs(d) / u
        worst = max(worst, z)
        every.append((low, upp, wn, d, z))
        if z > sigma:
            bad.append((low, upp, wn, d, z))
    return bad, worst, sorted(every, key=lambda t: -t[4])


# ---------------------------------------------------------------------------
# I. The unfound shape on disk
# ---------------------------------------------------------------------------
def unfind_in_iden2(index, log, apply=True):
    """Return level `index` to the shape of a level nobody has found.

    An unfound row of ``enlev.dat`` has an exact shape, and this restores it
    character for character: the adopted energy becomes the calculated one, the
    observed-minus-calculated becomes zero, and the star goes out.

    THE UNCERTAINTY COLUMN IS LEFT ALONE HERE, deliberately.  It means two
    different things on the two kinds of row - a measurement's uncertainty on a
    found level, and on an unfound one a prediction, the half-width of the
    window IDEN2 will search - and what belongs there is the rms of
    ``E_obs - E_calc`` over the levels of this level's configuration that are
    STILL found.  ``sync_IDEN2.configuration_uncertainties`` computes exactly
    that, for this level and every other unfound level of the configuration,
    from the energies it is about to write, and step J runs it.  Writing
    anything here would only be overwritten by something better.

    ``trans.dat`` keeps a copy of every level's energy and found flag in every
    block that has a row for it, and ``check_sync._iden_internal`` reports an
    ERROR when a copy no longer follows from ``enlev.dat``.  So every copy of
    this level's energy is brought with it, and every predicted wavenumber that
    was built from it recomputed.  Returns ``(e_obs before, e_calc)``.
    """
    enlev = IDEN.Enlev(ENLEV)
    e_before = enlev.e_obs(index)
    e_calc = round(enlev.e_calc(index), 3)
    unc = IDEN.number(enlev.record(index), IDEN.EN_UNC)
    enlev.set_measurement(index, unc, e_calc, False)
    if apply:
        IDEN.write_records(ENLEV, enlev.records, enlev.endings)
    log('   enlev.dat row %d: E_obs %.3f -> %.3f (the calculated position), '
        'O-C -> 0.000, star cleared' % (index, e_before, e_calc))
    log('   its uncertainty is left for sync_IDEN2.py, which writes the rms '
        'of this configuration\'s still-found levels there')

    trans = IDEN.Trans(TRANS)
    k = trans.header_of.get(index)
    if k is not None:
        trans.records[k] = IDEN.put(
            IDEN.put(trans.records[k], IDEN.TH_EOBS, '%13.3f' % e_calc),
            IDEN.TH_STAR, IDEN.STAR_OFF)
    n_part, n_wn = 0, 0
    for (owner, partner), j in trans.row_of.items():
        if index not in (owner, partner):
            continue
        rec = trans.records[j]
        if partner == index:
            rec = IDEN.put(rec, IDEN.TR_EPART, '%12.3f' % e_calc)
            rec = IDEN.put(rec, IDEN.TR_STAR, IDEN.STAR_OFF)
            n_part += 1
        try:
            pred = abs(enlev.e_obs(owner) - enlev.e_obs(partner))
        except ValueError:
            trans.records[j] = rec
            continue
        rec = IDEN.put(rec, IDEN.TR_WN, '%13.3f' % pred)
        trans.records[j] = rec
        n_wn += 1
    if apply:
        IDEN.write_records(TRANS, trans.records, trans.endings)
    log('   trans.dat: the block header unstarred, %d partner energ(ies) and '
        '%d predicted wavenumber(s) brought with it' % (n_part, n_wn))
    return e_before, e_calc


# ---------------------------------------------------------------------------
# The per-configuration search window, before and after
# ---------------------------------------------------------------------------
def configuration_window(path, index, log):
    """What the discard does to this configuration's search window.

    The same quantity ``sync_IDEN2.configuration_uncertainties`` writes: the
    rms of ``E_obs - E_calc`` over the FOUND levels of a configuration.  It is
    recomputed here only so that the report can quote it both ways, before and
    after, without waiting for sync_IDEN2.py to run - and so that the one case
    where the discard makes things worse can be caught and said out loud.

    That case is a configuration left with fewer than ``MIN_CFG_LEVELS`` found
    levels.  An rms over one or two levels measures nothing, so
    ``configuration_uncertainties`` leaves those rows exactly as they stand,
    carrying whatever ``enlev.dat`` already holds - which for a level nobody has
    ever touched is the placeholder 5000.000, the conversion from Cowan's
    output saying that nobody has filled the column in.  That is the truth and
    not a fault, but it is worth being told about, because a discard can push a
    configuration over that edge and silently return its unfound levels to the
    placeholder.

    Called before the discard it quotes both numbers; called after it, when the
    level is already unfound, there is only one state left to report and it
    says so rather than printing the same figure twice as a change.

    Returns ``(cfg, n_before, n_after, rms_before, rms_after, n_unfound)``.
    """
    import sync_IDEN2

    enlev = IDEN.Enlev(path)
    cfg = level_interchange.split_label(enlev.label(index))[0]
    still_found = enlev.known(index)
    omc, n_unfound = [], 0
    for other in sorted(enlev.row_of):
        if level_interchange.split_label(enlev.label(other))[0] != cfg:
            continue
        if not enlev.known(other):
            n_unfound += 1
            continue
        omc.append((other, enlev.e_obs(other) - round(enlev.e_calc(other), 3)))

    def rms(values):
        if not values:
            return float('nan')
        return round(math.sqrt(sum(v * v for v in values) / len(values)), 3)

    before = [v for _n, v in omc]
    after = [v for n, v in omc if n != index]
    mine = next((v for n, v in omc if n == index), None)
    # The unfound levels this window is written for are the ones unfound now
    # plus, while the discard is still to come, this level itself.
    n_window = n_unfound + (1 if still_found else 0)
    if still_found:
        log('   configuration %s: %d found level(s) before, %d after; the '
            'search window its %d unfound level(s) are given goes '
            '%.3f -> %.3f cm^-1'
            % (cfg, len(before), len(after), n_window, rms(before),
               rms(after)))
    else:
        log('   configuration %s: %d found level(s), and the %d unfound ones '
            'are searched in a window of %.3f cm^-1'
            % (cfg, len(after), n_window, rms(after)))
    if mine is not None:
        log('   this level\'s own O-C is %+.3f cm^-1, %.1f times that window: '
            'while it counted as found it was widening the search for every '
            'other level of %s' % (mine, abs(mine) / (rms(before) or 1.0), cfg))
    if len(after) < sync_IDEN2.MIN_CFG_LEVELS:
        log('   WARNING: %s is left with %d found level(s), below the %d '
            'sync_IDEN2.py needs to measure a window at all.  Its unfound '
            'levels will keep the 5000.000 placeholder, which says only that '
            'nobody knows - true, but a real loss of information that this '
            'discard caused.'
            % (cfg, len(after), sync_IDEN2.MIN_CFG_LEVELS))
    return cfg, len(before), len(after), rms(before), rms(after), n_window


# ---------------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description='Give up a level\'s measured position: release its lines, '
                    'take it out of the pipeline\'s level list and return it '
                    'to the levels nobody has found.')
    p.add_argument('level_id', nargs='?', default=None,
                   help='identifier of the level to discard, 059003.000561')
    p.add_argument('--iden2-row', type=int, default=None, metavar='N',
                   help='name the level by its row of enlev.dat instead')
    p.add_argument('--reason', default=None,
                   help='why the position is being given up; it goes into the '
                        'reason column of discarded_levels.csv and is what the '
                        'decision can be read back from')
    p.add_argument('--audit', default=AUDIT, metavar='PATH',
                   help='the table level_positions.py --audit writes, which is '
                        'the evidence the discard rests on '
                        '(default %(default)s)')
    p.add_argument('--ritz-sigma', type=float, default=INL.DEF_RITZ_SIGMA,
                   metavar='S',
                   help='a SURVIVING component of a freed line whose residual '
                        'exceeds this many times its own uncertainty stops the '
                        'run (default %(default)s)')
    p.add_argument('--max-rss', type=float, default=None, metavar='R',
                   help='stop if RSS/degrees_of_freedom comes out above this')
    p.add_argument('--rebuild', action='store_true',
                   help='run make_LOPT_input.py after the classification, so '
                        'that everything it accepts goes into the fit')
    p.add_argument('--force', action='store_true',
                   help='go on although the audit does not support the '
                        'discard, and run LOPT although its files are as it '
                        'left them; what was overridden is said in the log and '
                        'in the ledger')
    p.add_argument('--undo', action='store_true',
                   help='put back everything the last run wrote and stop.  '
                        'With no backup left - a discard made in an earlier '
                        'session - it takes the ledger row out and restores '
                        'the IDEN_level_ids.txt row from it')
    p.add_argument('--yes', action='store_true',
                   help='apply; without it the run writes nothing')
    p.add_argument('--dry-run', action='store_true',
                   help='run the whole sequence for real and then put every '
                        'file back')
    p.add_argument('--no-sync', action='store_true',
                   help='skip check_sync.py and sync_IDEN2.py.  The discarded '
                        'level then keeps the uncertainty it had in enlev.dat '
                        'instead of its configuration\'s search window')
    args = p.parse_args(argv)
    if not args.level_id and args.iden2_row is None:
        p.error('name the level, either as an identifier or with --iden2-row')
    if not args.undo and not args.reason:
        p.error('--reason is required: a discard that says nothing about why '
                'cannot be read back or argued with')
    return args


# ---------------------------------------------------------------------------
# Undoing a discard made in an earlier session
# ---------------------------------------------------------------------------
def undo_from_ledger(args, log):
    """Put the level back in the list, from the ledger alone.

    The run's own backups are gone - this is a discard made in an earlier
    session - so what there is to work from is the ``discarded_levels.csv`` row
    itself, and the ``iden2_row`` column of it is exactly why that column is
    there.  Deleting the row returns the level to the pipeline's level list,
    and restoring its ``IDEN_level_ids.txt`` line joins it to IDEN2 again.

    THE LEVEL IS NOT PUT BACK AT THE POSITION IT WAS DISCARDED FROM.  There was
    never any evidence for that one; that is what the discard said.  It comes
    back with no lines, which is ``insert_new_level.py``'s starting state, and
    is entered afresh at whatever position the search has since found.
    """
    fields, rows = read_discarded(DISCARDED)
    hit = None
    for r in rows:
        if args.level_id and (r.get('level_id') or '').strip() == args.level_id:
            hit = r
        elif args.iden2_row is not None \
                and (r.get('iden2_row') or '').strip() == str(args.iden2_row):
            hit = r
    if hit is None:
        log('nothing to undo: no backup in %s and no row in %s for %s'
            % (os.path.basename(BACKUP_DIR), os.path.basename(DISCARDED),
               args.level_id or 'IDEN2 row %d' % args.iden2_row))
        return 1
    level_id = (hit.get('level_id') or '').strip()
    try:
        index = int((hit.get('iden2_row') or '').strip())
    except ValueError:
        raise Abort('the %s row for %s carries no iden2_row, so there is '
                    'nothing to restore the map from.  Put the row back by '
                    'hand.' % (os.path.basename(DISCARDED), level_id))
    write_discarded(DISCARDED, fields,
                    [r for r in rows if r is not hit])
    log('   %s taken out of %s (discarded %s: %s)'
        % (level_id, os.path.basename(DISCARDED), hit.get('date') or '?',
           hit.get('reason') or ''))
    if level_id not in INL.read_id_map_rows(ID_MAP):
        INL.add_id_map_row(ID_MAP, level_id, index)
        log('   %s -> row %d put back into %s'
            % (level_id, index, os.path.basename(ID_MAP)))
    log('')
    log('%s is in the level list again, holding no line.  It is not back at '
        'the position it was discarded from - there was never any evidence for '
        'that one - so enter it where the search now puts it:' % level_id)
    log('   python level_positions.py --unknown %d' % index)
    log('   python insert_new_level.py --iden2-row %d' % index)
    return 0


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------
def main(argv=None):
    args = parse_args(argv)
    log = Log()
    try:
        if args.undo:
            saved = {p: os.path.join(BACKUP_DIR, os.path.basename(p))
                     for p in WRITABLE
                     if os.path.exists(os.path.join(BACKUP_DIR,
                                                    os.path.basename(p)))}
            if saved:
                INL.restore(saved, log)
                return 0
            return undo_from_ledger(args, log)
        today = datetime.date.today().strftime('%Y-%m-%d')
        log('discard_level.py - %s' % (args.level_id or 'IDEN2 row %d'
                                       % args.iden2_row))
        log('')
        saved = INL.preflight(WRITABLE, log, BACKUP_DIR)
        try:
            return _run(args, saved, log, today)
        except Held as exc:
            log('')
            log('HELD: %s' % exc)
            if args.dry_run:
                log('')
                log('--dry-run: putting every file back')
                INL.restore(saved, log)
            else:
                log('   what this run wrote is left in place; %s --undo puts '
                    'it back' % os.path.basename(__file__))
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
    except Abort as exc:
        log('')
        log('STOPPED: %s' % exc)
        return 1
    finally:
        log.save(LOGFILE)


def _run(args, saved, log, today):
    import classify_lines as CL

    # --- B. the level -------------------------------------------------------
    log('B. the level')
    id_rows = INL.read_id_map_rows(ID_MAP)
    already = None
    if args.level_id and args.level_id.strip() not in id_rows:
        # The map row is the first thing a discard takes out, so a level that
        # is not in the map may well be one this tool has already dealt with.
        # That is the no-op guard's business, not a broken lookup.
        already = discarded_row(DISCARDED, args.level_id.strip())
    elif args.iden2_row is not None \
            and args.iden2_row not in set(id_rows.values()):
        for rec in read_discarded(DISCARDED)[1]:
            if (rec.get('iden2_row') or '').strip() == str(args.iden2_row):
                already = rec
    if already is not None:
        log('')
        log('%s was discarded on %s already: %s'
            % ((already.get('level_id') or '').strip(),
               already.get('date') or '?', already.get('reason') or ''))
        log('There is nothing left for this run to do, and it has spent '
            'nothing: no LOPT, no classification.')
        log('%s --undo puts the level back into the list.'
            % os.path.basename(__file__))
        INL.restore(saved, log)
        return 0

    level_id, index = ML.resolve_level(args, log)
    levels_dict, _levels_list = CL.read_energy_levels()
    if level_id not in levels_dict:
        raise Abort('%s is not in the pipeline\'s level list, yet it is still '
                    'in %s.  Those two disagree, and the disagreement has to '
                    'be settled by hand before anything is written.'
                    % (level_id, os.path.basename(ID_MAP)))
    e_old = levels_dict[level_id].energy
    where = ML.where_the_energy_lives(level_id)
    log('   at %.4f cm^-1; its position is held in %s'
        % (e_old, os.path.basename(NEW_LEVELS if where == 'new_levels'
                                   else OVERRIDES)))

    # --- C. what it holds ---------------------------------------------------
    log('')
    log('C. what the level holds')
    holdings, lopt_rows = ML.gather_holdings(level_id, index, levels_dict, log)
    n_fitted = mark_fitted(holdings, lopt_rows)
    log('   %d assignment(s): %d in the LOPT input and %d of those carrying a '
        'weight, %d marked in IDEN2, %d ruled on in the ledger'
        % (len(holdings),
           sum(1 for h in holdings if h.in_lopt), n_fitted,
           sum(1 for h in holdings if h.partner_row is not None),
           sum(1 for h in holdings if h.ledger)))
    log('   the other %d record(s) are flagged "P": classifications already '
        'rejected, which LOPT prints and does not fit'
        % (sum(1 for h in holdings if h.in_lopt) - n_fitted))
    if not holdings and discarded_row(DISCARDED, level_id) is not None:
        log('')
        log('%s holds nothing anywhere and is in %s already; nothing to do.'
            % (level_id, os.path.basename(DISCARDED)))
        INL.restore(saved, log)
        return 0

    # --- D. the guards ------------------------------------------------------
    log('')
    log('D. the evidence')
    ln_r = check_the_evidence(level_id, args, log)

    calc_index = CL.read_transitions(levels_dict)
    lines = CL.read_observed_lines(levels_dict, calc_index)
    wns, ordered = ML.line_index(lines)
    calc_of = {(k[0], k[1]): (v.get('calc_intensity') or 0.0)
               for k, v in calc_index.items()}
    calc_of.update(INL.classification_intensities(CLASSIFICATIONS))

    log('')
    log('   what the discard releases:')
    release_table(holdings, lopt_rows, calc_of,
                  lambda wn: ML.line_at(wns, ordered, wn), log)
    log('')
    configuration_window(ENLEV, index, log)

    if not args.yes:
        log('')
        log('Nothing was written.  Run again with --yes to apply.')
        INL.restore(saved, log)
        return 0

    # --- G0. the fit as it stands -------------------------------------------
    log('')
    log('G0. the fit as it stands')
    rss_before, _dof = INL.call_LOPT(log, force=args.force)

    # --- E. the old position is taken apart ---------------------------------
    # Everything, with no surviving set to work out: that is the one place a
    # discard is simpler than a move.
    log('')
    log('E. releasing everything the level holds')
    for h in holdings:
        h.reason = REASON_DISCARDED
    # Only a record the fit uses can change a weight or leave a partner with a
    # residual; a "P" record going out of the file changes nothing at all.
    touched_wn = [h.wn for h in holdings if h.fitted]

    ledger_fields, ledger_all = ML.read_ledger_rows(LINE_DECISIONS)
    archived = [(rec, REASON_DISCARDED) for h in holdings for rec in h.ledger]
    if archived:
        want = {(round(float(r.get('wn_obs') or 0), 3),
                 (r.get('low_id') or '').strip(),
                 (r.get('upp_id') or '').strip()) for r, _w in archived}
        keep_rows, dropped = [], 0
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
        ML.write_ledger(LINE_DECISIONS, ledger_fields, keep_rows)
        ML.archive_ledger_rows(REMOVED, archived, level_id, today)
        log('   %d ledger row(s) taken out of %s and kept in %s'
            % (dropped, os.path.basename(LINE_DECISIONS),
               os.path.basename(REMOVED)))
    else:
        log('   no ledger row rules on any of its assignments')

    trans = IDEN.Trans(TRANS)
    n_cleared = 0
    for h in holdings:
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

    lopt_rows = INL.read_lopt_input(LOPT_INPUT)
    gone = {(INL.lopt_key(h.wn), h.low_id, h.upp_id) for h in holdings}
    before_n = len(lopt_rows)
    lopt_rows = [r for r in lopt_rows
                 if (INL.lopt_key(r['wn']), r['low_id'], r['upp_id'])
                 not in gone]
    n_out = before_n - len(lopt_rows)
    n_w = INL.reweigh(lopt_rows, touched_wn, calc_of, log)
    INL.write_lopt_input(LOPT_INPUT, lopt_rows)
    log('   %d record(s) taken out, %d weight(s) changed, %d record(s) in all'
        % (n_out, n_w, len(lopt_rows)))

    # --- F. the level leaves the list ---------------------------------------
    log('')
    log('F. the level leaves the level list')
    fields, rows = read_discarded(DISCARDED)
    reason = args.reason
    if args.force:
        reason += '  [--force: the audit did not support this discard]'
    rows = [r for r in rows if (r.get('level_id') or '').strip() != level_id]
    rows.append({'level_id': level_id, 'iden2_row': str(index),
                 'date': today,
                 'ln_R': '' if ln_r is None else '%.3f' % ln_r,
                 'reason': reason})
    rows.sort(key=lambda r: (r.get('level_id') or ''))
    write_discarded(DISCARDED, fields, rows)
    log('   %s written into %s' % (level_id, os.path.basename(DISCARDED)))
    remove_id_map_row(ID_MAP, level_id, log)

    if where == 'new_levels':
        nl_fields, nl_rows = INL.read_new_levels(NEW_LEVELS)
        kept = [r for r in nl_rows
                if (r.get('level_id') or '').strip() != level_id]
        if len(kept) != len(nl_rows):
            INL.write_new_levels(NEW_LEVELS, nl_fields, kept)
            log('   its row taken out of %s: the level was one of ours, and a '
                'level with no position has nothing to enter the pipeline with'
                % os.path.basename(NEW_LEVELS))
    ov_fields, ov_rows = ML.read_overrides(OVERRIDES)
    ov_kept = [r for r in ov_rows
               if (r.get('level_id') or '').strip() != level_id]
    if len(ov_kept) != len(ov_rows):
        ML.write_overrides(OVERRIDES, ov_fields, ov_kept)
        log('   its row taken out of %s: one file may not say where the level '
            'is while another says nobody knows'
            % os.path.basename(OVERRIDES))

    # --- G. LOPT, and the Ritz check inverted -------------------------------
    log('')
    log('G. LOPT')
    rss_after, _dof = INL.call_LOPT(log, force=args.force)
    if rss_before is not None and rss_after is not None:
        log('   RSS/degrees_of_freedom %.2f -> %.2f (%+.2f)'
            % (rss_before, rss_after, rss_after - rss_before))
    if args.max_rss is not None and rss_after is not None \
            and rss_after > args.max_rss:
        raise Abort('RSS/degrees_of_freedom came out at %.2f, above the %.2f '
                    'of --max-rss' % (rss_after, args.max_rss))
    bad, worst, every = partner_culprits(LOPT_LINES, touched_wn, level_id,
                                         args.ritz_sigma)
    log('   the components left on the %d line(s) this discard freed, worst '
        'first:' % len(set(INL.lopt_key(w) for w in touched_wn)))
    for low, upp, wn, d, z in every[:8]:
        log('     %11.3f  %s - %s  O-C = %+.3f cm^-1 (%.1f sigma)'
            % (wn, low, upp, d, z))
    if not every:
        log('     none: every freed line was this level\'s alone and is '
            'unidentified again')
    log('   the largest residual among them: %.2f sigma' % worst)
    if bad:
        log('   ABNORMAL: taking this level off these lines has left a '
            'surviving component where the fit cannot place it:')
        for low, upp, wn, d, z in sorted(bad, key=lambda t: -t[4]):
            log('     %11.3f  %s - %s  O-C = %+.3f cm^-1 (%.1f sigma)'
                % (wn, low, upp, d, z))
        raise Abort('%d surviving component(s) of the freed lines are further '
                    'than %g sigma from the fit.  A feature whose components '
                    'straddle the observed wavenumber is a real blend, and '
                    'this discard has taken one side of it away; decide those '
                    'lines before discarding the level.'
                    % (len(bad), args.ritz_sigma))

    # --- H. classify_lines and the ledger -----------------------------------
    log('')
    log('H. classify_lines.py and the ledger')
    # A DISCARD WRITES NO LEDGER ROW FOR THE LEVEL, and it must not.
    # classify_lines.check_forced_decisions() treats a ledger row naming a
    # level that is not in the level list as a fault and stops the run - for a
    # "reject" row as much as an "accept" one, since such a row can never rule
    # on anything.  A discarded level is exactly that, so a reject row written
    # here would abort the very classification run below, and every run after
    # it, until somebody found the row and deleted it.
    #
    # Nothing is lost by not writing one.  A move needs reject rows because the
    # level is still in the list and the classification would cheerfully
    # propose its old lines again; a discarded level generates no candidate at
    # all, so there is nothing for a verdict to rule on.  What was decided is
    # kept where step E put it: the rows the discard took out are in
    # line_decisions_removed.csv with the date and the reason, and the discard
    # itself is in discarded_levels.csv.
    if INL.run([sys.executable, 'classify_lines.py'], log)[0] != 0:
        raise Abort('classify_lines.py failed')
    verdicts = INL.classification_verdicts(CLASSIFICATIONS)
    # A level the classification has been told to drop cannot appear in its
    # output at all.  If it does, the reader and this tool disagree about what
    # a discard is, and that is a fault in the code rather than something a
    # ledger row should paper over.
    still = sorted(k for k, ok in verdicts.items()
                   if ok and level_id in (k[1], k[2]))
    if still:
        raise Abort('classify_lines.py still accepts %d assignment(s) for %s '
                    'although %s names it: %s.  The reader and this tool '
                    'disagree about what a discard is; this is a bug, not '
                    'something to write a ledger row for.'
                    % (len(still), level_id, os.path.basename(DISCARDED),
                       ', '.join('%.3f' % k[0] for k in still[:5])))
    left = sorted(k for k in INL.read_ledger_keys(LINE_DECISIONS)
                  if level_id in (k[1], k[2]))
    if left:
        raise Abort('%d row(s) of %s still name %s, which is no longer in the '
                    'level list; the next classification run would stop on '
                    'them: %s'
                    % (len(left), os.path.basename(LINE_DECISIONS), level_id,
                       ', '.join('%.3f' % k[0] for k in left[:5])))
    log('   the classification proposes nothing for %s, and no ledger row '
        'names it; what it used to hold is in %s'
        % (level_id, os.path.basename(REMOVED)))
    # Nothing is appended to line_decisions.csv by a discard; see above.

    if args.rebuild:
        if INL.run([sys.executable, 'make_LOPT_input.py'], log)[0] != 0:
            raise Abort('make_LOPT_input.py failed')
        rss_r, _dof = INL.call_LOPT(log, force=args.force)
        if rss_before is not None and rss_r is not None:
            log('   RSS/degrees_of_freedom %.2f -> %.2f (%+.2f)'
                % (rss_before, rss_r, rss_r - rss_before))
        bad, worst, _every = partner_culprits(LOPT_LINES, touched_wn,
                                              level_id, args.ritz_sigma)
        log('   the largest residual on the freed lines after the rebuild: '
            '%.2f sigma' % worst)
        if bad:
            raise Abort('%d surviving component(s) of the freed lines are '
                        'further than %g sigma from the rebuilt fit'
                        % (len(bad), args.ritz_sigma))

    # --- I. IDEN2 -----------------------------------------------------------
    log('')
    log('I. IDEN2')
    e_before, e_calc = unfind_in_iden2(index, log)

    # Every other level's component on a line this discard freed that the
    # classification now accepts and IDEN2 does not show.  A released line is
    # handed back to whatever else can have it, and which of them gets it is
    # the judgement this tool does not make.
    verdicts_final = INL.classification_verdicts(CLASSIFICATIONS)
    id_rows = INL.read_id_map_rows(ID_MAP)
    trans = IDEN.Trans(TRANS)
    marked_all = set()
    for (owner, partner), k in trans.row_of.items():
        obs = IDEN.assignment(trans.records[k])
        if IDEN.has_line(obs):
            marked_all.add((round(IDEN.obs_wavenumber(obs), 2), owner, partner))
    touched = {round(w, 2) for w in touched_wn}
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
        log('   %d assignment(s) of other levels on the lines this discard '
            'freed are now accepted by the classification and are not marked '
            'in IDEN2:' % len(unchecked))
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
        log('   This discard did not propose them and has not written them.  A '
            'line it released is free again and something else may now be able '
            'to have it - which is exactly the judgement to make on the '
            'screen.  Look at each one, then mark it in IDEN2 or rule on it in '
            '%s.' % os.path.basename(LINE_DECISIONS))

    # --- J. check_sync, then sync_IDEN2 -------------------------------------
    if args.no_sync:
        log('')
        log('J. skipped (--no-sync).  enlev.dat row %d still carries the '
            'uncertainty it had as a found level; only sync_IDEN2.py replaces '
            'it by the configuration\'s search window.' % index)
    else:
        log('')
        log('J. check_sync.py')
        status, _out = INL.run([sys.executable, 'check_sync.py'], log)
        if status > 1:
            log('   the ERROR findings of %s:' % os.path.basename(SYNC_REPORT))
            for rec in io.open(SYNC_REPORT, encoding='utf-8',
                               errors='replace'):
                if rec.strip().startswith('ERROR'):
                    log('     ' + rec.rstrip())
            INL.keep_failed(BACKUP_DIR, log)
            raise Abort('check_sync.py reported errors')
        log('   no errors%s' % (' (warnings only)' if status == 1 else ''))

        if unchecked:
            raise Held('%d assignment(s) on the lines this discard freed are '
                       'nobody\'s decision yet; sync_IDEN2.py is not run until '
                       'they are made.' % len(unchecked))

        log('')
        log('K. sync_IDEN2.py')
        if INL.run([sys.executable, 'sync_IDEN2.py'], log)[0] != 0:
            raise Abort('sync_IDEN2.py failed')
        log('')
        log('   the search window sync_IDEN2.py has just written for this '
            'configuration:')
        configuration_window(ENLEV, index, log)

    log('')
    log('L. done')
    log('   %s discarded: %d record(s) released, %d of them carrying a weight '
        'in the fit, and %d ledger row(s) taken out and kept in %s.  Its '
        'position at %.4f cm^-1 is given up and enlev.dat row %d is back at '
        'the calculated %.3f, unstarred.'
        % (level_id, len(holdings), n_fitted, len(archived),
           os.path.basename(REMOVED), e_old, index, e_calc))
    log('   The level is not gone: it is among the ones nobody has found, '
        'where  python level_positions.py --unknown %d  can search for it '
        'again.' % index)
    if args.dry_run:
        log('')
        log('--dry-run: putting every file back')
        INL.restore(saved, log)
    return 0


if __name__ == '__main__':
    sys.exit(main())
