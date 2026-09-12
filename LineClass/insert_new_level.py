"""Put a newly found energy level into the pipeline, from end to end.

WHAT THIS IS FOR
================
``level_positions.py --unknown`` finds a position where a calculated level
that nobody has ever found would explain a group of observed lines.  Accepting
that position is a decision the analyst makes in IDEN2, by eye; everything that
has to follow it is mechanical, and until now it was a dozen edits made by hand
across six files:

    new_levels.txt              the level itself, with the next free identifier
    IDEN2/IDEN_level_ids.txt    that identifier against its row of enlev.dat
    LOPT_input_lines.txt        one record per accepted line, and the weights
                                of every other component of the lines it shares
    line_decisions.csv          the verdicts that make the assignments stick
    IDEN2/trans.dat             the assignments as IDEN2 shows them
    IDEN2/enlev.dat             the level's measured energy, once LOPT has it

with three programs to be run in between, in the right order, and a check
after each of them.  A step forgotten leaves two files describing different
identifications, and the disagreement is silent.

This tool does the whole sequence, checks it, and can undo it.

WHAT IT DOES, IN ORDER
======================
A. **Preflight.**  Every file the run can write is tested for writability -
   Excel locks the workbook it has open - and then copied byte for byte into a
   run directory, so that any later step can put everything back.

B. **The level.**  A row is added to ``new_levels.txt`` if there is none for
   this row of ``enlev.dat``: the next free identifier, the adopted energy and
   J of ``enlev.dat``, the parity of the calculated level, and the two columns
   that say which calculated level it is (``iden2_row``, ``cowan_lid``).  An
   existing row is checked and left alone.  ``IDEN2/IDEN_level_ids.txt`` gets
   the identifier against the row if it lacks it.

C. **The lines.**  The assignments already marked in ``IDEN2/trans.dat`` for
   the level's block are taken as given, and WHEN THERE ARE ANY THEY ARE ALL
   THERE IS: a level whose lines have been gone through on the screen has been
   decided, and the tool adds nothing of its own to it.  Only for a level with
   no mark anywhere does it propose lines itself - each an observed line inside
   the Ritz window whose intensity the prediction can account for - and then A
   LINE MUCH STRONGER THAN PREDICTED IS NEVER PROPOSED: it is listed as left
   free, because a line whose strength the new level cannot explain belongs to
   some other transition, and taking it would both misplace this level and hide
   the real owner.  ``--propose`` asks for proposals beside the hand marks,
   ``--no-propose`` for none ever.  ``--reject WN=reason`` writes a verdict and
   assigns nothing.  Without ``--yes`` the run stops here, having written
   nothing, and prints the table.

D. **The LOPT input.**  One record per accepted assignment is inserted into
   ``LOPT_input_lines.txt`` in its own fixed-column layout.  Then, for every
   observed wavenumber the run touches, the weights of ALL its records without
   the ``P`` flag are recomputed together, in proportion to their calculated
   intensities: a line that was one transition's alone and is now shared has to
   be divided again, and the components that were already there are as much
   affected as the new one.  The ``P`` records - the rejected candidates LOPT
   carries but does not fit - keep weight zero.

E. **LOPT.**  ``lopt.bat LOPT.par`` (the Perl v5 build; never the Java jar,
   which has no centroid blend model), run once on the untouched input first
   so that there is a before to compare with.  Two things are then read out of
   the fit.

   ``RSS/degrees_of_freedom`` is LOPT's own measure of how well the whole fit
   holds together - the sum over every fitted line of the squared difference
   between its observed wavenumber and the one the levels imply, each divided
   by that line's uncertainty, per degree of freedom.  It is about 1 when the
   uncertainties are honest and the identifications right.  The before and
   after are both reported, so that a level bought by making the rest of the
   fit worse shows itself; ``--max-rss`` stops the run above a value.

   Then every line of the new level: an observed-minus-Ritz residual above
   four times the line's own uncertainty means the position is wrong or one of
   the lines is not its.  The culprits are named, every file is put back as it
   was, and the run stops.  A line shared with another transition has no
   residual of its own - the fit tests the blended feature, not the component -
   and those are listed separately rather than passed over in silence.
   Otherwise the level's energy in ``new_levels.txt`` is replaced by the
   optimized one.

F. **classify_lines.py and the ledger.**  The classification is run and its
   verdict on each of the run's lines compared with the intended one.  Every
   disagreement - an assignment rejected, or never even proposed - becomes a
   row of ``line_decisions.csv`` dated today with the reason ``IDEN2/LOPT``,
   and each ``--reject`` becomes a reject row with the reason given.  Rows the
   ledger already holds are left alone.  This repeats until nothing new has to
   be written, because a ledger row changes the classification and a line the
   first pass accepted of its own accord can be rejected by the next.  Three
   rounds without settling is an oscillation and stops the run.
   ``make_LOPT_input.py`` IS NOT RUN: it would rebuild the fit's input out of
   the classification and so put into the fit every assignment the
   classification accepts, including ones nobody has looked at.  ``--rebuild``
   asks for it when that is what is wanted.  ``--accept WN`` writes an accept
   row for every component of the observed line at ``WN`` that the
   classification accepts and this level is not part of - the ones a previous
   run listed as nobody's decision yet, once they have been looked at.

G. **IDEN2.**  What the classification finally accepts FOR THIS LEVEL is what
   ``IDEN2/trans.dat`` is made to show: marks made by hand are kept, missing
   ones are added.  A hand mark the classification will not accept even with
   the ledger rows this run wrote stops the run instead of being overwritten -
   that is a judgement to make by hand, not one to make silently.  Any other
   component of a touched line that the classification has newly accepted is
   named and left alone, for the same reason - unless ``--accept WN`` names
   its line, which says it has been looked at and is good, and then it is
   written like this level's own and holds the run no longer.

H. **check_sync.py.**  It must report no errors; warnings are printed and
   ignored.  An error puts every file back and stops the run - unless the
   findings are those unauthorized assignments, in which case the run is held
   where it is, with everything it wrote left in place to be looked at, and
   ``--undo`` puts it back.

I. **sync_IDEN2.py.**  Only when H is clean: it brings ``IDEN2/enlev.dat`` and
   the rest of IDEN2 into step with the fit.

J. **The report**, on screen and in ``insert_new_level.log``.

USAGE
=====
::

    python insert_new_level.py --iden2-row 742
    python insert_new_level.py --iden2-row 742 --yes \\
        --reject 91856.116="May add to pub line list as masked"
    python insert_new_level.py --iden2-row 742 --yes --rebuild \
        --reject 91856.116="May add to pub line list as masked" \
        --accept 95033.381 --accept 93276.982 --accept 90917.831="better CoG"

The first form writes nothing: it prints the proposal table and stops.
``--dry-run`` is the opposite kind of rehearsal - it runs the whole sequence
for real, programs and all, and then puts every file back as it was, which is
the only way to see what the classification and the fit will say.
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
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import cowan_gA
import make_LOPT_input
import output_files
import swap_line_assignments_IDEN as IDEN

# ---------------------------------------------------------------------------
# Where everything is
# ---------------------------------------------------------------------------
IDEN2_DIR = os.path.join(HERE, 'IDEN2')
ENLEV = os.path.join(IDEN2_DIR, 'enlev.dat')
TRANS = os.path.join(IDEN2_DIR, 'trans.dat')
GLUELEV = os.path.join(IDEN2_DIR, 'gluelev.dat')
ID_MAP = os.path.join(IDEN2_DIR, 'IDEN_level_ids.txt')
DLV = os.path.join(IDEN2_DIR, 'dlv.dat')

NEW_LEVELS = os.path.join(HERE, 'new_levels.txt')
LINE_DECISIONS = os.path.join(HERE, 'line_decisions.csv')
CLASSIFICATIONS = os.path.join(HERE, 'line_classifications.csv')
CLASSIFICATIONS_XLSX = os.path.join(HERE, 'line_classifications.xlsx')
LOPT_PAR = os.path.join(HERE, 'LOPT.par')
LOPT_INPUT = os.path.join(HERE, 'LOPT_input_lines.txt')
LOPT_LINES = os.path.join(HERE, 'LOPT_output_lines.txt')
LOPT_LEVELS = os.path.join(HERE, 'LOPT_output_levels.txt')
LOPT_FIXLEV = os.path.join(HERE, 'LOPT_fixlev.txt')
SYNC_REPORT = os.path.join(HERE, 'sync_report.txt')
LOGFILE = os.path.join(HERE, 'insert_new_level.log')
BACKUP_DIR = os.path.join(HERE, 'insert_new_level_backup')

# Every file the run can write.  The preflight tests all of them before any
# work is done, and each is copied into the run directory so that an abort can
# put it back.
WRITABLE = [NEW_LEVELS, LINE_DECISIONS, CLASSIFICATIONS, CLASSIFICATIONS_XLSX,
            LOPT_INPUT, LOPT_LINES, LOPT_LEVELS, LOPT_FIXLEV, LOPT_PAR,
            SYNC_REPORT, LOGFILE, TRANS, ENLEV, GLUELEV, ID_MAP]

# ---------------------------------------------------------------------------
# The rules the proposals obey
# ---------------------------------------------------------------------------
# Half-width of the Ritz window, in cm^-1: an observed line further than this
# from the predicted wavenumber is not a candidate at all.  Two wavenumbers
# out of the seventeen assigned to level 742 by hand differ by 1.6 cm^-1, so
# the window has to be wider than that to see them.
DEF_WINDOW = 2.0

# Half-width of the window a line has to fall in before the tool will propose
# it ITSELF, in cm^-1.  It is tighter than the window it looks in, and
# deliberately so: everything inside the wider window is shown, because the
# analyst may recognise a line the arithmetic cannot, but a line more than
# this far from the prediction is not a case simple enough to decide without
# looking.  The seventeen lines marked by hand for level 742 lie within
# 1.6 cm^-1 of their predictions and all but two within 0.9.
DEF_PROPOSE_WINDOW = 1.0

# How much stronger than predicted a line may be and still be proposed,
# measured in standard deviations of ln(I_obs/I_calc).  The uncertainty is
# the calculated one (u_calc = ln(1 + u%gA/100), often around 0.6, sometimes
# above 2) added in quadrature to the observed one, taken as ln 2 - the same
# figure classify_lines.py weeds with.  Three sigma is deliberately generous:
# the tool is not trying to reproduce the analyst's judgement, only to refuse
# the cases where a line is so much stronger than the new level could make it
# that it must belong to something else.  Nothing here overrules a hand mark
# in IDEN2, which is read as a decision already taken.
DEF_STRONG_SIGMA = 3.0

# The observed intensity's own uncertainty on the logarithmic scale.
U_OBS_LN = math.log(2.0)

# A Ritz residual above this many times the line's own uncertainty, in LOPT's
# fit, means the position or one of the lines is wrong.
DEF_RITZ_SIGMA = 4.0

REASON_ACCEPT = 'IDEN2/LOPT'
REASON_ADOPT = 'checked on the screen; adopted with the new level'


class Abort(Exception):
    """Something the run cannot go on from; every file is put back first."""


class Held(Exception):
    """The run did all it was authorized to do and stops there, with what it
    wrote left in place.

    This is not a failure and nothing is put back.  It is what happens when
    the sequence reaches something that is the analyst's to decide - an
    assignment nobody has looked at, a finding of check_sync.py - and the
    files have to stay as they are for that decision to be made on them.
    ``--undo`` puts them back if the decision goes the other way.
    """


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------
class Log(object):
    """Print, and keep a copy for the log file."""

    def __init__(self):
        self.lines = []

    def __call__(self, text=''):
        print(text)
        self.lines.append(text)

    def save(self, path):
        try:
            with io.open(path, 'w', encoding='utf-8', newline='\n') as fh:
                fh.write('\n'.join(self.lines) + '\n')
        except OSError as exc:
            print('could not write %s: %s' % (os.path.basename(path), exc))


def j_string(j_val: float) -> str:
    """``2.5`` as ``"5/2"``, ``3.0`` as ``"3"``."""
    twice = int(round(2 * j_val))
    return str(twice // 2) if twice % 2 == 0 else '%d/2' % twice


def next_free_id(level_ids) -> str:
    """The next identifier free: the largest in use, plus one.

    An identifier is ``059003.000042``: a spectrum part that never changes and
    a six-digit level number.  Taking one past the largest number in use is
    what makes it impossible for a level found now to collide with one of the
    published list, whose numbers are all smaller.
    """
    prefix, biggest = None, 0
    for lid in level_ids:
        head, _, tail = str(lid).partition('.')
        if not tail.isdigit():
            continue
        prefix = head
        biggest = max(biggest, int(tail))
    if prefix is None:
        raise Abort('no level identifier of the form 059003.000042 was found')
    return '%s.%06d' % (prefix, biggest + 1)


# ---------------------------------------------------------------------------
# A. Preflight and the backups
# ---------------------------------------------------------------------------
def preflight(paths, log) -> dict:
    """Test every file for writability, then copy it into the run directory.

    Returns ``{path: backup path}``.  ``require_writable`` raises SystemExit
    naming the files Excel is holding, which is the whole point of doing this
    first: the run stops before any of the work rather than after all of it.
    """
    output_files.require_writable(paths, 'file of this run')
    if os.path.isdir(BACKUP_DIR):
        shutil.rmtree(BACKUP_DIR)
    os.makedirs(BACKUP_DIR)
    saved = {}
    for path in paths:
        if not os.path.exists(path):
            continue
        dest = os.path.join(BACKUP_DIR, os.path.basename(path))
        shutil.copy2(path, dest)
        saved[path] = dest
    log('A. %d file(s) writable and backed up in %s'
        % (len(saved), os.path.basename(BACKUP_DIR)))
    return saved


def restore(saved: dict, log) -> None:
    """Put every backed-up file back, byte for byte."""
    for path, dest in saved.items():
        shutil.copy2(dest, path)
    log('   every file was put back as it was (%d file(s))' % len(saved))


# ---------------------------------------------------------------------------
# IDEN2: enlev.dat, the identifier map, the assignments already marked
# ---------------------------------------------------------------------------
def enlev_row(path: str, index: int):
    """``(E_obs, E_calc, J, label)`` of one row of ``enlev.dat``.

    ``E_obs`` is ``None`` when the row carries no measured energy: the level
    is one the calculation predicts and nobody has found yet, and where to put
    it is exactly the decision that has not been made.
    """
    with io.open(path, encoding='latin-1', newline='') as fh:
        for rec in fh:
            rec = rec.rstrip('\r\n')
            if not rec.strip():
                continue
            if int(IDEN.field(rec, IDEN.EN_INDEX)) != index:
                continue
            label = rec.partition('/')[2].partition('/')[0]
            try:
                e_obs = float(IDEN.field(rec, IDEN.EN_EOBS))
            except ValueError:
                e_obs = None
            if e_obs == 0.0 and index != 1:
                e_obs = None
            return (e_obs, float(IDEN.field(rec, IDEN.EN_ECALC)),
                    float(IDEN.field(rec, IDEN.EN_J)), label)
    raise Abort('%s has no row %d' % (os.path.basename(path), index))


def read_id_map_rows(path: str) -> dict:
    """``{level_id: IDEN2 row}`` from ``IDEN_level_ids.txt``."""
    out = {}
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        rdr = csv.DictReader(fh, delimiter='\t')
        for rec in rdr:
            lid = (rec.get('level_id') or '').strip()
            row = (rec.get('IDEN_id') or '').strip()
            if lid and row:
                out[lid] = int(row)
    return out


def add_id_map_row(path: str, level_id: str, index: int) -> None:
    """Append one ``level_id -> IDEN2 row`` line, keeping the file sorted."""
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        text = fh.read()
    lines = text.replace('\r\n', '\n').rstrip('\n').split('\n')
    head, body = lines[0], lines[1:]
    body.append('%s\t%d' % (level_id, index))
    body.sort(key=lambda s: s.split('\t')[0])
    with io.open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join([head] + body) + '\n')


def marked_assignments(trans: 'IDEN.Trans', index: int) -> dict:
    """``{partner row: (obs wavenumber, dlv row, intensity code)}`` for a block.

    Every ``+`` row of ``trans.dat`` that names ``index`` on either side and
    carries a line - the last field of its assignment is the row of
    ``dlv.dat``, and zero there means no line.
    """
    out = {}
    for (owner, partner), k in trans.row_of.items():
        if index not in (owner, partner):
            continue
        obs = IDEN.assignment(trans.records[k])
        if not IDEN.has_line(obs):
            continue
        other = partner if owner == index else owner
        code = obs[:IDEN.TR_OBS_WN[0]].strip()
        out[other] = (IDEN.obs_wavenumber(obs), IDEN.obs_row(obs),
                      int(code) if code else 0)
    return out


def read_dlv(path: str) -> dict:
    """``{row: (intensity code, wavenumber)}`` from IDEN2's observed line list."""
    out = {}
    with io.open(path, encoding='latin-1', newline='') as fh:
        for rec in fh:
            v = rec.split()
            if len(v) < 2:
                continue
            try:
                out[int(v[-1])] = (int(float(v[0])), float(v[1]))
            except ValueError:
                continue
    return out


def make_assignment(code: int, wn: float, omc: float, row: int) -> str:
    """The assignment tail of a ``trans.dat`` row, at its fixed widths."""
    text = '%6d%12.3f%11.3f%6d' % (code, wn, omc, row)
    width = IDEN.TR_OBS_ROW[1]
    if len(text) != width:
        raise Abort('built an assignment of %d characters, not %d: %r'
                    % (len(text), width, text))
    return text


# ---------------------------------------------------------------------------
# B. The level in new_levels.txt
# ---------------------------------------------------------------------------
NEW_LEVEL_COLUMNS = ['level_id', 'E', 'J', 'parity', 'iden2_row', 'cowan_lid',
                     'comment']


def read_new_levels(path: str) -> tuple:
    """``(fieldnames, rows)`` of ``new_levels.txt``; the delimiter is the
    extension's - tab for ``.txt``, comma for ``.csv``, as classify_lines.py
    reads it."""
    if not os.path.exists(path):
        return list(NEW_LEVEL_COLUMNS), []
    delim = '\t' if os.path.splitext(path)[1].lower() == '.txt' else ','
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        rdr = csv.DictReader(fh, delimiter=delim)
        return list(rdr.fieldnames or NEW_LEVEL_COLUMNS), list(rdr)


def write_new_levels(path: str, fields, rows) -> None:
    delim = '\t' if os.path.splitext(path)[1].lower() == '.txt' else ','
    with io.open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fields, delimiter=delim, lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in fields})


# ---------------------------------------------------------------------------
# C. The candidate lines
# ---------------------------------------------------------------------------
class Candidate(object):
    """One observed line offered to one transition of the new level."""

    def __init__(self, low_id, upp_id, rwn, line, calc_intens, u_calc):
        self.low_id = low_id
        self.upp_id = upp_id
        self.rwn = rwn
        self.line = line
        self.calc_intens = calc_intens
        self.u_calc = u_calc
        self.source = ''        # 'IDEN2' | 'proposed' | ''
        self.verdict = ''       # why it was taken or left

    @property
    def wn(self):
        return self.line.wavenumber

    @property
    def residual(self):
        return self.line.wavenumber - self.rwn

    @property
    def ln_ratio(self):
        if not self.calc_intens or self.line.intensity <= 0:
            return None
        return math.log(self.line.intensity / self.calc_intens)

    @property
    def z_intensity(self):
        """How many standard deviations the line is stronger than predicted."""
        r = self.ln_ratio
        if r is None:
            return None
        u = math.sqrt((self.u_calc or 1.0) ** 2 + U_OBS_LN ** 2)
        return r / u

    def key(self):
        return (self.low_id, self.upp_id)


def gather_candidates(level, levels_dict, calc_index, lines, window):
    """Every observed line within `window` of a predicted transition's Ritz
    wavenumber, as a list of Candidate."""
    ordered = sorted(lines, key=lambda l: l.wavenumber)
    wns = [l.wavenumber for l in ordered]
    out = []
    for (low_id, upp_id), d in calc_index.items():
        if level.level_id not in (low_id, upp_id):
            continue
        low, upp = levels_dict[low_id], levels_dict[upp_id]
        rwn = upp.energy - low.energy
        if rwn <= 0:
            continue
        a = bisect.bisect_left(wns, rwn - window)
        b = bisect.bisect_right(wns, rwn + window)
        for line in ordered[a:b]:
            out.append(Candidate(low_id, upp_id, rwn, line,
                                 d.get('calc_intensity'), d.get('u_calc')))
    out.sort(key=lambda c: -c.rwn)
    return out


def choose(candidates, marked_wn, strong_sigma, propose_window, propose=None):
    """Decide what each candidate is: taken from IDEN2, proposed, or left.

    The hand marks in IDEN2 are decisions already taken and are never
    second-guessed.

    WHETHER THE TOOL PROPOSES ANYTHING AT ALL depends on whether it finds any
    such mark.  A level whose lines are already marked in ``trans.dat`` has
    been gone through by eye, line by line, on the screen where the plates and
    the branch structure can be seen; the marks are the whole of the answer,
    and a line the analyst passed over was passed over for a reason the
    arithmetic here cannot see.  So when there are marks, they are taken and
    nothing else is - the tool adds no assignment of its own.  Only for a
    level with no mark anywhere, one straight out of ``level_positions.py
    --unknown``, does it propose lines itself.  `propose` overrides that:
    ``True`` proposes even beside hand marks, ``False`` never proposes, and
    ``None`` (the default) is the rule just described.

    Among the lines it does propose, one transition may take at most one line
    and one line at most one transition of this level - a level cannot have
    two of its transitions on the same feature - and a line much stronger than
    predicted is refused outright, whatever else is on offer.
    """
    for c in candidates:
        if any(abs(c.wn - w) < 5e-3 for w in marked_wn):
            c.source = 'IDEN2'
            c.verdict = 'marked in IDEN2'

    if propose is None:
        propose = not any(c.source == 'IDEN2' for c in candidates)
    if not propose:
        for c in candidates:
            if not c.source:
                c.verdict = ("left: this level's lines are the ones marked in "
                             'IDEN2')
        return candidates

    free = [c for c in candidates if not c.source]
    # Refuse the too-strong ones first, so that they cannot win a transition
    # away from a line the prediction can account for.
    rest = []
    for c in free:
        z = c.z_intensity
        if z is not None and z > strong_sigma:
            c.verdict = ('left free: %.0f times stronger than predicted (%.1f '
                         'sigma)' % (math.exp(c.ln_ratio), z))
        elif abs(c.residual) > propose_window:
            c.verdict = ('left: %.3f cm^-1 from the prediction, too far to '
                         'take without looking' % abs(c.residual))
        else:
            rest.append(c)

    taken_trans, taken_wn = set(), set()
    for c in candidates:
        if c.source == 'IDEN2':
            taken_trans.add(c.key())
            taken_wn.add(round(c.wn, 4))
    # Closest in wavenumber first: the best explanation of a line gets it.
    for c in sorted(rest, key=lambda c: abs(c.residual)):
        if c.key() in taken_trans:
            c.verdict = 'left: the transition already has a line'
            continue
        if round(c.wn, 4) in taken_wn:
            c.verdict = 'left: the line already has a transition of this level'
            continue
        c.source = 'proposed'
        c.verdict = 'proposed'
        taken_trans.add(c.key())
        taken_wn.add(round(c.wn, 4))
    return candidates


def candidate_table(candidates, log):
    log('   %-13s %-13s %11s %11s %7s %10s %10s %6s  %s'
        % ('lower', 'upper', 'Ritz wn', 'obs wn', 'O-C', 'I_obs', 'I_calc',
           'sigma', 'verdict'))
    for c in candidates:
        z = c.z_intensity
        log('   %-13s %-13s %11.3f %11.3f %+7.3f %10.4g %10.4g %6s  %s'
            % (c.low_id, c.upp_id, c.rwn, c.wn, c.residual,
               c.line.intensity, c.calc_intens or 0.0,
               '%+.1f' % z if z is not None else '-', c.verdict))


# ---------------------------------------------------------------------------
# D. The LOPT input
# ---------------------------------------------------------------------------
def read_lopt_input(path: str) -> list:
    """The records of ``LOPT_input_lines.txt``, parsed by their fixed columns.

    Each is a dict with the numbers the weighting needs and, under ``raw``, the
    record exactly as it stands.  Reading it back rather than rebuilding it
    from the classification table is what lets one level be added to a fit
    without disturbing a single other record - and a record is written out
    again from its fields only when this run actually changed it.  That matters
    because the file is not uniform: most weights are the fractions
    make_LOPT_input.py writes, but a few blends carry the raw calculated
    intensities of an older run instead (``35318.``), which the fraction format
    cannot even express.  Those records are passed through untouched.
    """
    fields = make_LOPT_input.FIELDS
    rows = []
    with io.open(path, encoding='ascii', newline='') as fh:
        for rec in fh:
            rec = rec.rstrip('\r\n')
            if not rec.strip():
                continue
            got = {}
            for name, (first, last) in fields.items():
                got[name] = rec[first - 1:last].strip()
            rows.append({
                'wn': float(got['wavenumber']),
                'unc': float(got['uncertainty']),
                'intens': float(got['intensity']),
                'low_id': got['lower_level'],
                'upp_id': got['upper_level'],
                'flag': got['flags'],
                'weight': float(got['weight'] or 0.0),
                'raw': rec,
            })
    return rows


def record_text(r) -> str:
    """One record: as it stands if this run did not change it, rebuilt if it
    did."""
    if r.get('raw') is not None:
        return r['raw']
    return make_LOPT_input.format_line(
        r['wn'], r['unc'], r['intens'], r['low_id'], r['upp_id'],
        r['flag'], r['weight'])


def write_lopt_input(path: str, rows) -> None:
    rows = sorted(rows, key=lambda r: -r['wn'])
    with io.open(path, 'w', encoding='ascii', newline='') as fh:
        for r in rows:
            fh.write(record_text(r) + make_LOPT_input.EOL_LINES)


def reweigh(rows, wavenumbers, calc_of, log) -> int:
    """Divide each touched observed line's weight among its unflagged records.

    This is the rule the analyst asked for and the one make_LOPT_input.py
    applies when it builds the file from scratch: the components of a blend
    share the line in proportion to their calculated intensities, and a record
    carrying the ``P`` flag - a candidate LOPT is shown but does not fit - is
    not a component and keeps weight zero.  Adding one transition to a line
    that already had one therefore changes the OTHER component's weight too,
    which is why every record of every touched wavenumber is rewritten and not
    only the new ones.
    """
    touched = set(round(w, 4) for w in wavenumbers)
    n = 0
    for wn in sorted(touched, reverse=True):
        group = [r for r in rows
                 if round(r['wn'], 4) == wn and 'P' not in r['flag'].upper()]
        if not group:
            continue
        proxy = [{'calc_intens': calc_of.get((r['low_id'], r['upp_id']), 0.0)}
                 for r in group]
        weights = make_LOPT_input.blend_weights(proxy)
        for r, p in zip(group, proxy):
            new = weights[id(p)]
            if abs(new - r['weight']) > 5e-5:
                log('   weight %11.3f  %s - %s  %.4f -> %.4f'
                    % (r['wn'], r['low_id'], r['upp_id'], r['weight'], new))
                n += 1
                r['weight'] = new
                r['raw'] = None          # rebuild this one, keep the others
    return n


# ---------------------------------------------------------------------------
# E. LOPT, and the Ritz check
# ---------------------------------------------------------------------------
def run(cmd, log, cwd=HERE):
    """Run one program of the chain, echoing its last lines.

    Returns ``(status, output)``: the whole of what the program wrote is
    handed back, because some of what matters - LOPT's residual sum of
    squares among it - is printed long before the program's last line and
    would otherwise fall off the end of the echo.
    """
    log('   $ ' + ' '.join(cmd))
    env = dict(os.environ)
    parent = os.path.dirname(HERE)
    env['PYTHONPATH'] = parent + os.pathsep + env.get('PYTHONPATH', '')
    proc = subprocess.run(cmd, cwd=cwd, env=env, shell=False,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, errors='replace')
    out = proc.stdout or ''
    tail = [s for s in out.splitlines() if s.strip()][-12:]
    for s in tail:
        log('     | ' + s)
    return proc.returncode, out


# LOPT prints one line that says how well the whole fit holds together:
#
#     RSS/degrees_of_freedom = 1.16 (5261 degrees_of_freedom)
#
# RSS is the sum, over every observed line the fit uses, of the squared
# difference between the observed wavenumber and the one the fitted levels
# imply, each divided by that line's own statistical uncertainty; degrees of
# freedom is the number of lines less the number of levels the fit determines.
# The ratio is therefore about 1 when the wavenumber uncertainties are honest
# and the identifications are right, and it grows when a line is put where it
# does not belong.  It is the one number that says whether a new level has
# been paid for by making everything else fit worse.
RSS_RE = 'RSS/degrees_of_freedom'


def lopt_rss(output: str):
    """``(RSS/dof, degrees of freedom)`` out of what LOPT printed."""
    for rec in (output or '').splitlines():
        if RSS_RE not in rec:
            continue
        _head, _sep, tail = rec.partition('=')
        parts = tail.replace('(', ' ').split()
        try:
            return float(parts[0]), int(parts[1])
        except (IndexError, ValueError):
            return None, None
    return None, None


def run_lopt(log):
    """Run LOPT and return ``(RSS/dof, degrees of freedom)``."""
    # A .bat file is not an executable image, so it is handed to the command
    # interpreter rather than started directly.  lopt.bat is the Perl v5 build
    # of LOPT; the Java jar under F:/J/LOPT has no centroid blend model and
    # renames its outputs, and must never be used here.
    cmd = (['cmd', '/c', 'lopt.bat', os.path.basename(LOPT_PAR)]
           if os.name == 'nt' else ['lopt.bat', os.path.basename(LOPT_PAR)])
    status, out = run(cmd, log)
    if status != 0:
        raise Abort('lopt.bat returned %d' % status)
    rss, dof = lopt_rss(out)
    if rss is None:
        log('   LOPT printed no %s line' % RSS_RE)
    else:
        log('   %s = %.2f (%d degrees of freedom)' % (RSS_RE, rss, dof))
    return rss, dof


def lopt_lines(path: str) -> list:
    """The rows of ``LOPT_output_lines.txt`` as dicts, tab-separated with a
    header; ``_`` is LOPT's way of writing an empty cell."""
    with io.open(path, encoding='latin-1', newline='') as fh:
        text = fh.read().replace('\r\n', '\n')
    rows = [r for r in text.split('\n') if r.strip()]
    head = rows[0].split('\t')
    out = []
    for rec in rows[1:]:
        d = dict(zip(head, rec.split('\t')))
        out.append({k: ('' if v.strip() in ('_', '') else v.strip())
                    for k, v in d.items()})
    return out


def ritz_culprits(path, pairs, sigma):
    """How each of the run's lines sits in the fit.

    Returns ``(culprits, worst, all, without)``: the assignments whose
    observed-minus-Ritz residual exceeds `sigma` times the line's own
    uncertainty, the largest such ratio, every one of them as
    ``(low, upp, wn, residual, sigmas)`` so that the report can show the
    picture and not only the failures, and the ones the fit gives no residual
    of at all, as ``(low, upp, wn, flag)``.

    THE LAST OF THOSE IS NOT A LINE THAT WAS DROPPED.  A line several
    transitions share is fitted by LOPT as one blended feature, its centroid
    against the intensity-weighted mean of the components' Ritz wavenumbers;
    there is then one residual for the feature and none for any single
    component, and LOPT writes ``_`` in the component's O-C column.  The
    assignment is in the fit and is carrying its weight; it simply has no
    number of its own to test, so this check cannot say anything about it and
    the report says so rather than leaving it out in silence.
    """
    bad, worst, every, without = [], 0.0, [], []
    for r in lopt_lines(path):
        key = (r.get('L1', ''), r.get('L2', ''))
        if key not in pairs:
            continue
        try:
            wn = float(r['wn_o'])
        except (KeyError, ValueError):
            wn = 0.0
        if r.get('F', ''):
            without.append((key[0], key[1], wn, r.get('F', '')))
            continue
        try:
            d = float(r['dWO-C'])
            u = float(r['uWnOTot'])
        except (KeyError, ValueError):
            without.append((key[0], key[1], wn, 'blended'))
            continue
        if u <= 0:
            without.append((key[0], key[1], wn, 'no uncertainty'))
            continue
        z = abs(d) / u
        worst = max(worst, z)
        every.append((key[0], key[1], wn, d, z))
        if z > sigma:
            bad.append((key[0], key[1], wn, d, z))
    return bad, worst, sorted(every, key=lambda t: -t[4]), without


def lopt_level_energy(path: str, level_id: str):
    """The optimized energy of one level in ``LOPT_output_levels.txt``."""
    with io.open(path, encoding='latin-1', newline='') as fh:
        for rec in fh:
            if level_id in rec:
                for tok in rec.replace('\t', ' ').split():
                    try:
                        val = float(tok)
                    except ValueError:
                        continue
                    if val > 1000.0 and tok != level_id:
                        return val
    return None


def refresh_level_energy(level_id: str, log) -> None:
    """Put LOPT's optimized energy for the level into ``new_levels.txt``.

    The next classification then starts from where the fit put the level and
    not from the position it was guessed at, so its candidate windows and its
    Ritz tests are the fit's own.
    """
    optimized = lopt_level_energy(LOPT_LEVELS, level_id)
    if optimized is None:
        return
    fields, rows = read_new_levels(NEW_LEVELS)
    for r in rows:
        if (r.get('level_id') or '').strip() != level_id:
            continue
        if abs(float(r['E']) - optimized) < 5e-4:
            return
        log('   %s: E %s -> %.3f cm^-1 (LOPT)' % (level_id, r['E'], optimized))
        r['E'] = '%.3f' % optimized
    write_new_levels(NEW_LEVELS, fields, rows)


# ---------------------------------------------------------------------------
# F. The ledger
# ---------------------------------------------------------------------------
def read_ledger_keys(path: str) -> set:
    """The ``(wn, low, upp)`` the ledger already rules on, wavenumber rounded
    to three decimals so that a row written to fewer digits still matches."""
    keys = set()
    if not os.path.exists(path):
        return keys
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        for rec in csv.DictReader(fh):
            try:
                wn = round(float(rec['wn_obs']), 3)
            except (KeyError, TypeError, ValueError):
                continue
            keys.add((wn, (rec.get('low_id') or '').strip(),
                      (rec.get('upp_id') or '').strip()))
    return keys


def append_ledger(path: str, rows) -> None:
    """Append verdicts to ``line_decisions.csv``, header and all if it is new."""
    exists = os.path.exists(path)
    fields = ['wn_obs', 'low_id', 'upp_id', 'decision', 'date', 'reason']
    if exists:
        with io.open(path, encoding='utf-8-sig', newline='') as fh:
            fields = csv.DictReader(fh).fieldnames or fields
    with io.open(path, 'a' if exists else 'w', encoding='utf-8',
                 newline='') as fh:
        w = csv.DictWriter(fh, fields, lineterminator='\n')
        if not exists:
            w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in fields})


def classification_intensities(path: str) -> dict:
    """``{(low, upp): calc_intens}`` from ``line_classifications.csv``."""
    out = {}
    if not os.path.exists(path):
        return out
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        for rec in csv.DictReader(fh):
            low = (rec.get('low_id') or '').strip()
            upp = (rec.get('upp_id') or '').strip()
            try:
                out[(low, upp)] = float(rec.get('calc_intens') or 0.0)
            except ValueError:
                continue
    return out


def classification_verdicts(path: str) -> dict:
    """``{(wn rounded, low, upp): accepted}`` from ``line_classifications.csv``."""
    out = {}
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        for rec in csv.DictReader(fh):
            low = (rec.get('low_id') or '').strip()
            upp = (rec.get('upp_id') or '').strip()
            if not (low and upp):
                continue
            try:
                wn = round(float(rec['wn_obs']), 3)
            except (KeyError, TypeError, ValueError):
                continue
            out[(wn, low, upp)] = float(rec.get('accepted') or 0) == 1
    return out


def classification_rows(path: str) -> dict:
    """``{wavenumber rounded to 0.01: [row, ...]}`` of ``line_classifications``.

    Every transition the classification put on one observed line, in the order
    the file holds them, each row carrying the extra key ``xlsx_row``: the row
    of ``line_classifications.xlsx`` it is, so that the analyst can open the
    workbook and go straight there.  The header is row 1 of the workbook and
    the first transition row 2, which is what the ``+ 2`` below counts.
    """
    out = {}
    if not os.path.exists(path):
        return out
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        for n, rec in enumerate(csv.DictReader(fh)):
            try:
                wn = round(float(rec['wn_obs']), 2)
            except (KeyError, TypeError, ValueError):
                continue
            rec['xlsx_row'] = n + 2
            out.setdefault(wn, []).append(rec)
    return out


def _num(rec, key, default=None):
    """One numeric field of a classification row, ``default`` if unreadable."""
    try:
        return float(rec.get(key) or 0.0)
    except (TypeError, ValueError):
        return default


def line_dossier(wn, rows, before, level_id, log) -> None:
    """Print every transition the classification puts on one observed line.

    ``rows`` are that line's rows of the classification as it stands now and
    ``before`` the same line's rows as the classification stood before this
    run, keyed by ``(low_id, upp_id)``, so that a verdict this run changed can
    be shown as the change it is.  ``level_id`` is the new level, whose own
    transition on the line is the reason any of the others may have moved.

    The columns are the ones a decision is made on: the observed intensity of
    the line as a whole, then per transition the predicted intensity, the
    difference between the observed wavenumber and the Ritz wavenumber the two
    levels give (``dif_wn_O-C``, cm^-1), the grade classify_lines gave it, the
    verdict, and the reason behind that verdict.
    """
    if not rows:
        return
    head = rows[0]
    log('     %11.3f +/- %.3f cm^-1   I_obs %10.4g   %s   '
        'line_classifications.xlsx rows %d-%d'
        % (wn, _num(head, 'unc_wn_obs', 0.0) or 0.0,
           _num(head, 'obs_intens', 0.0) or 0.0,
           (head.get('char') or '').strip() or '-',
           rows[0]['xlsx_row'], rows[-1]['xlsx_row']))
    log('       %-13s %-13s %10s %9s %5s %-19s %s'
        % ('lower', 'upper', 'I_calc', 'dif_O-C', 'grade', 'verdict',
           'reason'))
    for rec in rows:
        low = (rec.get('low_id') or '').strip()
        upp = (rec.get('upp_id') or '').strip()
        ok = _num(rec, 'accepted', 0.0) == 1
        mark = '*' if level_id in (low, upp) else ' '
        was = before.get((low, upp))
        verdict = 'accepted' if ok else 'rejected'
        if was is None:
            verdict += ' (new)'
        elif (_num(was, 'accepted', 0.0) == 1) != ok:
            verdict += ' (was %s)' % ('accepted' if not ok else 'rejected')
        omc = _num(rec, 'dif_wn_O-C')
        omc_txt = '%+9.3f' % omc if omc is not None else '        -'
        reason = (rec.get('notes2') or '').strip()
        if was is not None and omc is not None:
            omc_was = _num(was, 'dif_wn_O-C')
            if omc_was is not None and '%+.3f' % omc_was != '%+.3f' % omc:
                reason += '   [dif_O-C was %+.3f]' % omc_was
        log('     %s %-13s %-13s %10.4g %s %5s %-19s %s'
            % (mark, low, upp, _num(rec, 'calc_intens', 0.0) or 0.0,
               omc_txt, (rec.get('grade') or '').strip() or '-', verdict,
               reason))
    log('       * = a transition of the new level, the reason the rest of the '
        'line may have moved')


# ---------------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------------
def parse_reject(values):
    """``--reject WN=reason`` into ``[(wavenumber, reason), ...]``."""
    out = []
    for text in values or []:
        wn, sep, reason = text.partition('=')
        if not sep:
            raise SystemExit('--reject wants WN=reason, not %r' % text)
        try:
            out.append((float(wn), reason.strip()))
        except ValueError:
            raise SystemExit('--reject: %r is not a wavenumber' % wn)
    return out


def parse_accept(values):
    """``--accept WN`` or ``--accept WN=reason`` into ``[(wn, reason), ...]``.

    The reason is what goes into ``line_decisions.csv``; left out, it is
    ``REASON_ADOPT``, which says only that the assignment was looked at.
    """
    out = []
    for text in values or []:
        wn, sep, reason = text.partition('=')
        try:
            out.append((float(wn), reason.strip() if sep and reason.strip()
                        else REASON_ADOPT))
        except ValueError:
            raise SystemExit('--accept: %r is not a wavenumber' % wn)
    return out


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description='Put a newly found level into the pipeline, from the row '
                    'of IDEN2/enlev.dat it was found at to a fit and a set of '
                    'files that agree with each other.')
    p.add_argument('--iden2-row', type=int, default=None, metavar='N',
                   help='the row of IDEN2/enlev.dat the level is')
    p.add_argument('--level-id', metavar='ID', default=None,
                   help='the identifier to give it; the next one free by '
                        'default')
    p.add_argument('--energy', type=float, default=None, metavar='E',
                   help='where to put the level, cm^-1.  Needed when the row '
                        'of enlev.dat has no measured energy yet - the '
                        'position is the analyst\'s decision and there is '
                        'nowhere else to read it from.  Given as well as one '
                        'in enlev.dat, it overrides it.')
    p.add_argument('--reject', action='append', metavar='WN=REASON',
                   help='record a reject verdict for the line at WN, with the '
                        'reason given; may be repeated')
    p.add_argument('--accept', action='append', metavar='WN[=REASON]',
                   help='adopt the other assignment(s) the classification has '
                        'newly accepted on the observed line at WN - the ones '
                        "a previous run listed as nobody's decision yet.  "
                        'They are written into IDEN2 and into '
                        "line_decisions.csv like this level's own, and no "
                        'longer hold the run.  May be repeated')
    p.add_argument('--window', type=float, default=DEF_WINDOW, metavar='CM',
                   help='half-width of the Ritz window, cm^-1 (default %(default)s)')
    p.add_argument('--propose-window', type=float,
                   default=DEF_PROPOSE_WINDOW, metavar='CM',
                   help='propose a line only if it is within this of the '
                        'prediction, cm^-1 (default %(default)s)')
    p.add_argument('--strong-sigma', type=float, default=DEF_STRONG_SIGMA,
                   metavar='S',
                   help='refuse a line more than S sigma stronger than '
                        'predicted (default %(default)s)')
    p.add_argument('--ritz-sigma', type=float, default=DEF_RITZ_SIGMA,
                   metavar='S',
                   help='a LOPT residual above S times the line uncertainty '
                        'stops the run (default %(default)s)')
    p.add_argument('--propose', dest='propose', action='store_true',
                   default=None,
                   help='propose lines of its own even for a level whose '
                        'lines are already marked in IDEN2 (by default it '
                        'proposes only when there is no mark at all)')
    p.add_argument('--no-propose', dest='propose', action='store_false',
                   help='take only what IDEN2 marks, whatever there is')
    p.add_argument('--max-rss', type=float, default=None, metavar='R',
                   help='stop if LOPT\'s RSS/degrees_of_freedom comes out '
                        'above R; by default it is only reported, before and '
                        'after')
    p.add_argument('--rebuild', action='store_true',
                   help='rebuild LOPT_input_lines.txt with make_LOPT_input.py '
                        'after the ledger rows are written.  Off by default: '
                        'a rebuild puts into the fit every assignment the '
                        'classification accepts, including ones this run '
                        'never proposed and nobody has looked at')
    p.add_argument('--undo', action='store_true',
                   help='put back every file of the last run from '
                        + os.path.basename(BACKUP_DIR) + ' and do nothing '
                        'else')
    p.add_argument('--yes', action='store_true',
                   help='apply; without it the run stops after the proposal '
                        'table and writes nothing')
    p.add_argument('--dry-run', action='store_true',
                   help='go through the whole sequence, programs and all, and '
                        'then put every file back, so that nothing is left '
                        'changed; implies --yes')
    p.add_argument('--no-sync', action='store_true',
                   help='skip check_sync.py and sync_IDEN2.py')
    args = p.parse_args(argv)
    if args.dry_run:
        args.yes = True          # a rehearsal has to get past the proposals
    if args.iden2_row is None and not args.undo:
        p.error('--iden2-row is required')
    return args


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------
def main(argv=None):
    args = parse_args(argv)
    rejects = parse_reject(args.reject)
    accepts = parse_accept(args.accept)
    log = Log()
    if args.undo:
        saved = {p: os.path.join(BACKUP_DIR, os.path.basename(p))
                 for p in WRITABLE
                 if os.path.exists(os.path.join(BACKUP_DIR,
                                                os.path.basename(p)))}
        if not saved:
            log('nothing to undo: %s holds no backup' % BACKUP_DIR)
            return 1
        restore(saved, log)
        return 0
    today = datetime.date.today().strftime('%-m/%-d/%Y'
                                           if os.name != 'nt'
                                           else '%#m/%#d/%Y')
    log('insert_new_level.py - IDEN2 row %d' % args.iden2_row)
    log('')

    saved = preflight(WRITABLE, log)
    try:
        return _run(args, rejects, accepts, saved, log, today)
    except Held as exc:
        log('')
        log('HELD: %s' % exc)
        if args.dry_run:
            log('')
            log('--dry-run: putting every file back')
            restore(saved, log)
        else:
            log('   what this run wrote is left in place; %s --undo puts it '
                'back' % os.path.basename(__file__))
        return 2
    except Abort as exc:
        log('')
        log('STOPPED: %s' % exc)
        restore(saved, log)
        return 1
    except BaseException as exc:
        # Anything at all - a bug here, a keyboard interrupt - must leave the
        # files as they were.  Half-written pipeline files are worse than no
        # run at all, and they are silent.
        log('')
        log('STOPPED by an unexpected %s: %s' % (type(exc).__name__, exc))
        restore(saved, log)
        raise
    finally:
        log.save(LOGFILE)


def _run(args, rejects, accepts, saved, log, today):
    import classify_lines as CL

    index = args.iden2_row
    log('')
    log('B. the level')

    # --- which calculated level this row is ---------------------------------
    trans_table = cowan_gA.read_transitions(log=lambda m: log('   ' + m))
    calc_levels = cowan_gA.levels(trans_table)
    enlev_levels = cowan_gA.read_enlev_levels(ENLEV)
    id_of_row = cowan_gA.read_id_map(ID_MAP)
    mapping, _report = cowan_gA.match_to_enlev(calc_levels, enlev_levels,
                                               id_of_row)
    lid_of_row = {v: k for k, v in mapping.items()}
    cowan_lid = lid_of_row.get(index)
    if cowan_lid is None:
        raise Abort('enlev.dat row %d matches no calculated level' % index)
    row_calc = calc_levels[calc_levels['lid'] == cowan_lid].iloc[0]
    e_obs, e_calc, j_val, label = enlev_row(ENLEV, index)
    parity = 'o' if row_calc['parity'] == 'odd' else 'e'
    # Where the level goes.  enlev.dat carries a measured energy once the
    # position has been accepted in IDEN2 and written there; until then the
    # row holds only the calculated position, which is nowhere near good
    # enough to search on, and the energy has to be given.
    if args.energy is not None:
        if e_obs is not None and abs(e_obs - args.energy) > 5e-4:
            log('   --energy %.3f overrides the %.3f cm^-1 of enlev.dat row %d'
                % (args.energy, e_obs, index))
        e_obs = args.energy
    elif e_obs is None:
        raise Abort('enlev.dat row %d has no measured energy - the level has '
                    'not been fixed at a position there.  Say where to put it '
                    'with --energy E (the calculated position is %.3f cm^-1, '
                    'which is not a position to search on).' % (index, e_calc))
    log('   enlev.dat row %d: E = %.3f cm^-1, J = %s, parity %s, %s '
        '(Cowan level %d)'
        % (index, e_obs, j_string(j_val), parity, label.strip(), cowan_lid))

    # --- the row of new_levels.txt ------------------------------------------
    fields, rows = read_new_levels(NEW_LEVELS)
    for name in NEW_LEVEL_COLUMNS:
        if name not in fields:
            fields.append(name)
    existing = [r for r in rows
                if str(r.get('iden2_row') or '').strip() == str(index)]
    id_rows = read_id_map_rows(ID_MAP)
    by_row = {v: k for k, v in id_rows.items()}
    if existing:
        level_id = existing[0]['level_id'].strip()
        log('   %s is already in %s at %s cm^-1; left as it stands'
            % (level_id, os.path.basename(NEW_LEVELS), existing[0]['E']))
        if args.energy is not None and abs(float(existing[0]['E'])
                                           - args.energy) > 5e-4:
            log('   --energy moves it to %.3f cm^-1' % args.energy)
            existing[0]['E'] = '%.3f' % args.energy
            if args.yes:
                write_new_levels(NEW_LEVELS, fields, rows)
        new_row = None
    else:
        level_id = args.level_id or next_free_id(
            list(by_row.values())
            + [r['level_id'].strip() for r in rows if r.get('level_id')])
        new_row = {'level_id': level_id, 'E': '%.3f' % e_obs,
                   'J': j_string(j_val), 'parity': parity,
                   'iden2_row': str(index), 'cowan_lid': str(cowan_lid),
                   'comment': 'Found with level_positions --unknown %d.' % index}
        log('   %s to be added at %.3f cm^-1, J = %s, parity %s'
            % (level_id, e_obs, j_string(j_val), parity))

    map_row_needed = by_row.get(index) != level_id
    if map_row_needed and index in by_row:
        raise Abort('%s already gives row %d to %s, not to %s'
                    % (os.path.basename(ID_MAP), index, by_row[index],
                       level_id))

    # --- the state of the pipeline -----------------------------------------
    if new_row is not None and args.yes:
        write_new_levels(NEW_LEVELS, fields, rows + [new_row])
    if map_row_needed and args.yes:
        add_id_map_row(ID_MAP, level_id, index)

    log('')
    log('C. the lines')
    levels_dict, _levels_list = CL.read_energy_levels()
    if level_id not in levels_dict:      # not applied yet: add it in memory
        from models import EnergyLevel
        levels_dict[level_id] = EnergyLevel(
            level_id=level_id, energy=e_obs, parity=parity,
            J_str=j_string(j_val), J_val=j_val, is_new=1, is_added=1,
            iden2_row=index, cowan_lid=cowan_lid)
    calc_index = CL.read_transitions(levels_dict)
    lines = CL.read_observed_lines(levels_dict, calc_index)

    trans = IDEN.Trans(TRANS)
    marked = marked_assignments(trans, index)
    marked_wn = [wn for wn, _row, _code in marked.values()]
    log('   %d assignment(s) already marked in %s for row %d'
        % (len(marked), os.path.basename(TRANS), index))

    level = levels_dict[level_id]
    candidates = gather_candidates(level, levels_dict, calc_index, lines,
                                   args.window)
    choose(candidates, marked_wn, args.strong_sigma, args.propose_window,
           args.propose)
    if marked and args.propose is None:
        log('   the lines of this level are taken from IDEN2 alone; nothing '
            'is proposed here (--propose to propose as well)')

    # a --reject wavenumber is never assigned, whatever else was decided
    reject_wn = [wn for wn, _r in rejects]
    for c in candidates:
        if any(abs(c.wn - wn) < 5e-3 for wn in reject_wn):
            c.source = ''
            c.verdict = 'rejected on the command line'
    candidate_table(candidates, log)

    accepted = [c for c in candidates if c.source]
    left = [c for c in candidates
            if not c.source and c.verdict.startswith('left free')]
    log('')
    log('   %d assignment(s): %d marked in IDEN2, %d proposed here'
        % (len(accepted),
           sum(1 for c in accepted if c.source == 'IDEN2'),
           sum(1 for c in accepted if c.source == 'proposed')))
    if left:
        log('   %d line(s) left free for another transition:' % len(left))
        for c in left:
            log('     %11.3f  %s' % (c.wn, c.verdict))

    if not args.yes:
        log('')
        log('Nothing was written.  Run again with --yes to apply.')
        restore(saved, log)
        return 0

    unmarked = [wn for wn, _r, _c in marked.values()
                if not any(abs(c.wn - wn) < 5e-3 for c in candidates)]
    if unmarked:
        raise Abort('%d line(s) marked in IDEN2 for row %d are outside the '
                    'Ritz window or have no calculated transition: %s'
                    % (len(unmarked), index,
                       ', '.join('%.3f' % w for w in sorted(unmarked))))

    # --- the fit as it stands, before anything is added ---------------------
    # LOPT is run once on the untouched input so that there is something to
    # compare the new fit with.  Without it the RSS/degrees_of_freedom the run
    # ends with is a number with nothing to measure it against, and the whole
    # question - has this level been paid for by making the rest of the fit
    # worse? - cannot be answered.
    log('')
    log('D0. the fit as it stands')
    rss_before, _dof = run_lopt(log)

    # --- D. the LOPT input --------------------------------------------------
    log('')
    log('D. the LOPT input')
    # What the weights are divided in proportion to.  make_LOPT_input.py takes
    # this from the calc_intens column of line_classifications.csv - the
    # predicted intensity as the classification finally left it, adjusted by
    # the levels' intensity factors - so the same column is used here, and the
    # calculated table is the fallback for the pairs that are not in it yet.
    calc_of = {(k[0], k[1]): (v.get('calc_intensity') or 0.0)
               for k, v in calc_index.items()}
    calc_of.update(classification_intensities(CLASSIFICATIONS))
    lopt_rows = read_lopt_input(LOPT_INPUT)
    present = {(round(r['wn'], 4), r['low_id'], r['upp_id'])
               for r in lopt_rows}
    n_new = 0
    for c in accepted:
        key = (round(c.wn, 4), c.low_id, c.upp_id)
        if key in present:
            continue
        lopt_rows.append({'wn': c.wn, 'unc': c.line.wn_uncertainty,
                          'intens': c.line.intensity, 'low_id': c.low_id,
                          'upp_id': c.upp_id, 'flag': '', 'weight': 1.0,
                          'raw': None})
        n_new += 1
    n_w = reweigh(lopt_rows, [c.wn for c in accepted], calc_of, log)
    write_lopt_input(LOPT_INPUT, lopt_rows)
    log('   %d record(s) inserted, %d weight(s) changed, %d record(s) in all'
        % (n_new, n_w, len(lopt_rows)))

    # --- E. LOPT ------------------------------------------------------------
    log('')
    log('E. LOPT')
    rss_after, _dof = run_lopt(log)
    if rss_before is not None and rss_after is not None:
        log('   RSS/degrees_of_freedom %.2f -> %.2f (%+.2f)'
            % (rss_before, rss_after, rss_after - rss_before))
        if rss_after > rss_before:
            log('   the fit as a whole is worse than it was; every line of '
                'the new level is worth looking at again')
    if args.max_rss is not None and rss_after is not None \
            and rss_after > args.max_rss:
        raise Abort('RSS/degrees_of_freedom came out at %.2f, above the %.2f '
                    'of --max-rss' % (rss_after, args.max_rss))
    pairs = {c.key() for c in accepted}
    bad, worst, every, without = ritz_culprits(LOPT_LINES, pairs,
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
    log('   the largest Ritz residual of the new level: %.2f sigma' % worst)
    if bad:
        log('   ABNORMAL Ritz mismatches:')
        for low, upp, wn, d, z in sorted(bad, key=lambda t: -t[4]):
            log('     %11.3f  %s - %s  O-C = %+.3f cm^-1 (%.1f sigma)'
                % (wn, low, upp, d, z))
        raise Abort('%d line(s) of the new level are further than %g sigma '
                    'from the fit' % (len(bad), args.ritz_sigma))

    refresh_level_energy(level_id, log)

    # --- F. classify_lines and the ledger -----------------------------------
    # The classification is run and every verdict of its that differs from
    # what this run intends becomes a ledger row.  It is a loop rather than a
    # single pass because the ledger changes the classification: a line the
    # first run accepted of its own accord can be rejected by the next, once
    # the other assignments around it have moved, and that one needs a ledger
    # row of its own.  Two rounds settle it in practice; three without
    # settling means something is oscillating and is worth looking at by hand.
    #
    # WHAT THIS DELIBERATELY DOES NOT DO is run make_LOPT_input.py.  That
    # program builds the whole of LOPT_input_lines.txt afresh out of the
    # classification, and so puts into the fit every assignment the
    # classification accepts - including ones this run never proposed, that
    # nobody has looked at on the screen, and that the analyst may well mean
    # to deny.  The fit this run makes holds exactly the records step D put in
    # it, and nothing else.  --rebuild asks for the rebuild explicitly, for
    # when the classification's own verdicts have already been gone through.
    log('')
    log('F. classify_lines.py and the ledger')
    # What IDEN2 already shows, as {(wn to 0.01, {row, row})}: an assignment
    # --accept names that is on the screen already needs no ledger row, and
    # writing one would turn a published identification nobody questioned into
    # a decision this run claims to have taken.
    on_screen = set()
    if accepts:
        _tr = IDEN.Trans(TRANS)
        for (_owner, _partner), _k in _tr.row_of.items():
            _obs = IDEN.assignment(_tr.records[_k])
            if IDEN.has_line(_obs):
                on_screen.add((round(IDEN.obs_wavenumber(_obs), 2),
                               frozenset((_owner, _partner))))
        del _tr
    ledger_rows = []
    written = set()
    settled = False
    for attempt in range(1, 4):
        log('   round %d' % attempt)
        if run([sys.executable, 'classify_lines.py'], log)[0] != 0:
            raise Abort('classify_lines.py failed')
        verdicts = classification_verdicts(CLASSIFICATIONS)
        have = read_ledger_keys(LINE_DECISIONS)
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
        if attempt == 1:
            # --accept names an observed line, not a pair: what is adopted is
            # every component of that line the classification accepts and this
            # level is not part of - which is exactly the list a previous run
            # printed as nobody's decision yet.  A component the
            # classification rejects is left rejected: --accept says "the ones
            # you showed me are good", never "take everything on this line".
            for wn, reason in accepts:
                hit = 0
                for (kwn, low, upp), ok in sorted(verdicts.items()):
                    if not ok or abs(kwn - wn) > 5e-3                             or level_id in (low, upp):
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
            for wn, reason in rejects:
                for c in candidates:
                    if abs(c.wn - wn) > 5e-3:
                        continue
                    key = (round(c.wn, 3), c.low_id, c.upp_id)
                    if key in have or key in written:
                        continue
                    fresh.append({'wn_obs': '%.4f' % c.wn, 'low_id': c.low_id,
                                  'upp_id': c.upp_id, 'decision': 'reject',
                                  'date': today, 'reason': reason})
                    written.add(key)
                    log('     reject %11.3f  %s - %s  (%s)'
                        % (c.wn, c.low_id, c.upp_id, reason))
        if fresh:
            append_ledger(LINE_DECISIONS, fresh)
            ledger_rows.extend(fresh)
            log('     %d row(s) added to %s'
                % (len(fresh), os.path.basename(LINE_DECISIONS)))
        else:
            log('     the ledger already says everything this run needs it to '
                'say')

        if args.rebuild:
            if run([sys.executable, 'make_LOPT_input.py'], log)[0] != 0:
                raise Abort('make_LOPT_input.py failed')
            rss_r, _dof = run_lopt(log)
            if rss_before is not None and rss_r is not None:
                log('     RSS/degrees_of_freedom %.2f -> %.2f (%+.2f)'
                    % (rss_before, rss_r, rss_r - rss_before))
            if args.max_rss is not None and rss_r is not None \
                    and rss_r > args.max_rss:
                raise Abort('RSS/degrees_of_freedom came out at %.2f after '
                            'the rebuild, above the %.2f of --max-rss'
                            % (rss_r, args.max_rss))
            bad, worst, _every, _w = ritz_culprits(LOPT_LINES, pairs,
                                                   args.ritz_sigma)
            log('     the largest Ritz residual after the rebuild: %.2f sigma'
                % worst)
            if bad:
                for low, upp, wn, d, z in sorted(bad, key=lambda t: -t[4]):
                    log('       %11.3f  %s - %s  O-C = %+.3f cm^-1 (%.1f sigma)'
                        % (wn, low, upp, d, z))
                raise Abort('%d line(s) are further than %g sigma from the '
                            'rebuilt fit' % (len(bad), args.ritz_sigma))
            refresh_level_energy(level_id, log)
        if not fresh:
            settled = True
            break
    if not settled:
        raise Abort('the classification and the ledger did not settle in three '
                    'rounds; look at the assignments by hand')

    # --- G. IDEN2 -----------------------------------------------------------
    log('')
    log('G. IDEN2')
    # What this run assigned is what IDEN2 has to show: that is the very thing
    # check_sync.py compares, and a line in the fit that the screen does not
    # show is what it calls an error.
    #
    # ONLY THE NEW LEVEL'S OWN ASSIGNMENTS ARE WRITTEN.  Adding a level to a
    # line that already had components changes the intensity accounting on
    # that line, and the classification can turn one of the other components
    # from rejected into accepted on the strength of it.  Those are real
    # assignments to real levels and they may well be right, but they are not
    # what this run was authorized to make, and writing them into trans.dat
    # would put them on the screen as settled without anyone having looked at
    # the branch, the plate or the intensity.  They are listed below instead,
    # and the run stops before sync_IDEN2.py so that they can be looked at.
    verdicts_final = classification_verdicts(CLASSIFICATIONS)
    adopt_wn = {round(wn, 2) for wn, _r in accepts}
    # An assignment named by --accept has been looked at, so it is written
    # like this level's own instead of being listed and held on.
    final = {k: v for k, v in verdicts_final.items()
             if level_id in (k[1], k[2])
             or (v and round(k[0], 2) in adopt_wn)}
    if adopt_wn:
        log('   the other assignments on %s are adopted on the command line '
            'and are written as well'
            % ', '.join('%.2f' % w for w in sorted(adopt_wn, reverse=True)))
    trans = IDEN.Trans(TRANS)
    marked = marked_assignments(trans, index)
    marked_wn = {round(wn, 2) for wn, _row, _code in marked.values()}
    orphan = [wn for wn in marked_wn
              if not any(round(k[0], 2) == wn and ok and level_id in (k[1], k[2])
                         for k, ok in final.items())]
    if orphan:
        raise Abort('%d line(s) marked in IDEN2 for row %d are not accepted by '
                    'the classification even with the ledger rows this run '
                    'wrote: %s.  Decide them by hand - accept them in '
                    'line_decisions.csv, or take the mark out in IDEN2 - '
                    'before running this again.'
                    % (len(orphan), index,
                       ', '.join('%.3f' % w for w in sorted(orphan))))

    dlv = read_dlv(DLV)
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
        trans.records[k] = trans.records[k][:IDEN.TR_OBS] + make_assignment(
            code, dlv_wn, dlv_wn - rwn, dlv_row)
        log('   %11.3f  %s - %s  marked (IDEN2 rows %d - %d)'
            % (dlv_wn, low_id, upp_id, row_low, row_upp))
        n_marked += 1
    if n_marked:
        IDEN.write_records(TRANS, trans.records, trans.endings)
    log('   %d assignment(s) written into %s, %d were already there, %d could '
        'not be placed'
        % (n_marked, os.path.basename(TRANS), n_already, n_missed))

    # Every other component of a line this run touched that the classification
    # now accepts and IDEN2 does not show.  These are the ones nobody has
    # looked at; they are named, not written.
    touched = {round(c.wn, 2) for c in accepted}
    touched.update(round(wn, 2) for wn, _r in rejects)
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
        log('   %d other assignment(s) on the lines this run touched are now '
            'accepted by the classification and are not marked in IDEN2:'
            % len(unchecked))
        now_rows = classification_rows(CLASSIFICATIONS)
        before_all = classification_rows(saved.get(CLASSIFICATIONS)
                                         or CLASSIFICATIONS)
        for wn, low_id, upp_id, row_low, row_upp in unchecked:
            log('')
            log('     %11.3f  %s - %s  (IDEN2 rows %d - %d)'
                % (wn, low_id, upp_id, row_low, row_upp))
            key = round(wn, 2)
            was_rows = {((r.get('low_id') or '').strip(),
                         (r.get('upp_id') or '').strip()): r
                        for r in before_all.get(key, [])}
            line_dossier(wn, now_rows.get(key, []), was_rows, level_id, log)
        log('')
        log('   They are not this level\'s and this run did not propose them, '
            'so it has not written them.  Look at each one; then run this '
            'tool again, naming the lines you found good with --accept and '
            'the rest with --reject, and it will write them and go on to '
            'check_sync.py and sync_IDEN2.py itself:')
        log('')
        log('     python %s --iden2-row %d --rebuild --yes%s'
            % (os.path.basename(__file__), index,
               ''.join(' --accept %.3f' % wn
                       for wn in sorted({round(u[0], 3) for u in unchecked}))))
        log('')
        log('   --rebuild is what puts them into the fit: without it the LOPT '
            'input keeps only the records this run wrote.  A line you do not '
            'want goes in as --reject WN=reason instead, or straight into %s.'
            % os.path.basename(LINE_DECISIONS))

    # --- H. check_sync ------------------------------------------------------
    if args.no_sync:
        log('')
        log('H/I. skipped (--no-sync)')
    else:
        log('')
        log('H. check_sync.py')
        status, _out = run([sys.executable, 'check_sync.py'], log)
        if status > 1:
            log('   the ERROR findings of %s:' % os.path.basename(SYNC_REPORT))
            shutil.copy2(SYNC_REPORT,
                         os.path.join(BACKUP_DIR, 'sync_report.failed.txt'))
            for rec in io.open(SYNC_REPORT, encoding='utf-8',
                               errors='replace'):
                if rec.strip().startswith('ERROR'):
                    log('     ' + rec.rstrip())
            log('   the report of the failed run is kept as %s'
                % os.path.join(os.path.basename(BACKUP_DIR),
                               'sync_report.failed.txt'))
            if unchecked:
                # The findings are the assignments named just above: they are
                # errors only in the sense that the fit and the screen do not
                # yet agree, and making them agree is the decision this run
                # deliberately did not take.  Nothing is put back - the files
                # have to stay as they are for the decision to be made on
                # them.
                raise Held('check_sync.py reports the %d assignment(s) above '
                           'as findings.  Decide them, then run this tool '
                           'again with --accept/--reject as printed above; '
                           '--yes and --rebuild alone will not decide them.'
                           % len(unchecked))
            raise Abort('check_sync.py reported errors')
        log('   no errors%s' % (' (warnings only)' if status == 1 else ''))

        if unchecked:
            raise Held('%d assignment(s) on the lines this run touched are '
                       'nobody\'s decision yet; sync_IDEN2.py is not run '
                       'until they are made.  Name their lines with --accept '
                       'or --reject as printed above.' % len(unchecked))

        log('')
        log('I. sync_IDEN2.py')
        if run([sys.executable, 'sync_IDEN2.py'], log)[0] != 0:
            raise Abort('sync_IDEN2.py failed')

    # --- J. the report ------------------------------------------------------
    log('')
    log('J. done')
    log('   level %s at IDEN2 row %d, %d line(s) assigned '
        '(%d of them proposed here), %d left free, %d ledger row(s) written'
        % (level_id, index, len(accepted),
           sum(1 for c in accepted if c.source == 'proposed'),
           len(left), len(ledger_rows)))
    if args.dry_run:
        log('')
        log('--dry-run: putting every file back')
        restore(saved, log)
    return 0


if __name__ == '__main__':
    sys.exit(main())
