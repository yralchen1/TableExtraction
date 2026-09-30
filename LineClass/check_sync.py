#!/usr/bin/env python
"""Is everything in this directory talking about the same identification?

The work passes through five files that each hold a copy of the same facts -
which observed line belongs to which pair of levels, and where each level sits
- and nothing keeps those copies in step automatically:

    Pr3_lines.xlsx            the measured lines (the only real input)
        |   classify_lines.py, obeying line_decisions.csv and
        |   revised_level_energies.csv
        v
    line_classifications.csv / .xlsx   every candidate transition, accepted
        |                              or not
        |   make_LOPT_input.py
        v
    LOPT_input_lines.txt      the accepted ones, weighted
        |   lopt.bat
        v
    LOPT_output_lines.txt     the fit
    LOPT_output_levels.txt    the optimized level energies
        |   by hand, rarely
        v
    IDEN2/enlev.dat           the level list IDEN2 shows on screen
    IDEN2/trans.dat           its predicted transitions and their assignments
    IDEN2/dlv.dat             the observed lines it shows them against, which
                              are Pr3_lines.xlsx again, in IDEN2's own form

A step skipped anywhere leaves two files disagreeing, and the disagreement is
silent: every tool downstream keeps working, on stale numbers.  The usual
symptom is a report that recommends something already done - level_positions.py
proposing a move to an energy the level was moved to yesterday - which is
indistinguishable, when reading the report, from a recommendation worth acting
on.

This script reads all of them and says where they differ.  It changes nothing.

WHAT IS A MISMATCH AND WHAT IS ONLY DRIFT
-----------------------------------------
The IDEN2 files are brought up to date by hand and only when there is reason
to, so their level energies are expected to lag the current LOPT fit by a
little: a hundredth or two of a wavenumber is the ordinary state of affairs and
is not worth reporting.  Two thresholds separate that from a real difference:

    --drift  (default 0.5 cm^-1)  below this, an energy difference is the
                                  expected lag and is only counted, not listed
    --gross  (default 5.0 cm^-1)  at or above this, the two files are
                                  describing different levels, not the same
                                  level known to different precision

Everything else - a line accepted in one file and absent from another, a
decision in the ledger the classification does not obey, a transition assigned
in IDEN2 that the classification has withdrawn - has no tolerance at all and is
reported item by item.

THE THREE SEVERITIES
--------------------
    ERROR  the files contradict each other; a report built on them can be
           wrong in a way that no reading of it would reveal
    WARN   the files disagree in a way that is expected, or is confined to one
           artefact that can simply be rebuilt
    ok     checked, nothing to say

The exit status is 2 if anything is an ERROR, 1 if anything is a WARN, 0 if
everything is ok, so the script can gate a pipeline.

WHAT IT CHECKS
--------------
 1. Freshness.  Every artefact against the files it is built from, by
    modification time.  This is the cheapest check and catches most of it.
 2. line_classifications.csv against line_classifications.xlsx.  They are
    written from one table by classify_lines.py, so a difference means one of
    them has been edited by hand or one is left over from an earlier run.
 3. line_classifications against LOPT_input_lines.txt.  What has to agree is
    the fit: every accepted classification carrying weight in the transitions
    file, nothing else carrying weight, and each of them carrying the right
    share of its observed line.  The share, not the number in the weight
    column: LOPT normalises the weights of one observed line to sum to one, so
    a blend assigned by hand with the calculated intensities pasted straight
    out of Icalc.xlsx fits exactly as the pipeline's own fractions would, and
    only a share that is out by more than one per cent is an error.
 4. LOPT_input_lines.txt against LOPT_output_lines.txt: LOPT must have been
    run on the transitions file that is on disk now.
 5. LOPT_output_lines.txt against LOPT_output_levels.txt: the level energies
    the line output quotes against the level output's own.
 6. revised_level_energies.csv against the LOPT levels: a revision recorded in
    the ledger but not visible in the fit was never carried out.
 7. IDEN2/trans.dat against IDEN2/enlev.dat: the partner energies, the found
    flags and the predicted wavenumbers trans.dat carries are copies of what
    enlev.dat holds, and IDEN2 rewrites them together.
 8. IDEN2/enlev.dat against LOPT_output_levels.txt, joined by
    IDEN2/IDEN_level_ids.txt: the energy of every level, and the membership of
    the two lists both ways.
 9. IDEN2/trans.dat against the accepted classifications: which transitions
    IDEN2 shows as identified and which the classification actually accepts.
10. IDEN2/dlv.dat against the set's own line list, joined on the wavenumber
    every row should be carrying.  This is the file the eye reads a
    measurement off - the wavenumber a transition is looked for near, and the
    uncertainty that says how far from the prediction a line may stand and
    still be the line - and IDEN2 rewrites only what is edited on its screen,
    so a set seeded by copying another set's IDEN2 directory keeps that set's
    wavenumbers and uncertainties until sync_IDEN2.py is run.  A row that
    instead matches the set's wn_key column, the wavenumber that never moves,
    is named as what it is: a row the wavelength calibration correction has
    never been applied to.  The uncertainty there is a wavelength uncertainty
    in angstroms, stated to four decimals, and it is converted with the row's
    own wavelength before being compared with the uncertainty the pipeline
    uses: the line list's own, or the value inflated_unc_lines.txt sets for
    the line where that is larger.  The file's own arithmetic is checked
    first: one line number per row and one wavenumber per line, because
    trans.dat names every observed line by its line number.  The line numbers
    need not follow the order of the rows: a line inserted in IDEN2 takes the
    next free number and stands in its place by wavenumber.
11. IDEN2/trans.dat against IDEN2/dlv.dat: every identification names its
    line by the line number, and dlv.dat must have that line.  An
    identification also carries a wavenumber, but IDEN2 ignores it, and so
    does this script: the wavenumber of an identified line - here, and in
    check 9 - is the one on the dlv.dat row its line number names.
12. Stability.  unstable_candidates.csv - the assignments classify_lines.py
    withdrew because they oscillate - and line_decisions.csv, whose verdicts
    the classification must obey exactly.

THE PROBLEM LISTS
-----------------
Every list of transitions is written so that it can be worked through at the
screen without a second lookup.  The rows are ordered by the energy of the
upper level and then of the lower one - the order the same transitions appear
in in IDEN2 and in the workbook - and each row carries both levels with
IDEN2's own level numbers beside them, both energies, and three words saying
where the transition stands:

    .xlsx   accepted      line_classifications.xlsx accepts it here
            rejected      it is a candidate on this line and was turned down
            elsewhere     it is a candidate, but on other lines only
            missing       it is not a candidate anywhere
    LOPT    included      on this line in LOPT_input_lines.txt with a weight
                          and no P flag
            not incl. P   in the file, but flagged predicted or of zero
                          weight: it does not pull on the levels
            elsewhere     in the file on other lines only
            not incl. --  not in the file at all
    IDEN2   assigned      trans.dat identifies it with this line
            elsewhere     trans.dat identifies it with another line
            not assigned  IDEN2 predicts it and has no line for it
            no trans row  IDEN2 does not predict it at all, which for a real
                          pair of levels means trans.dat is stale

Every one of the three is answered for the line the row names and not for the
pair of levels in general.  The same pair can be a candidate on several
observed lines - moving an assignment leaves the old row in place to be
rejected, which is how the pipeline is made to turn down an assignment it
would otherwise keep - so a status read off another line would print a row
that is in perfect order as a problem.

The terminal shows the first --list rows of each list, enough to see the shape
of the trouble; the whole of every list is written to sync_report.txt, because
a list cut off after a dozen rows cannot be worked through and the rows that
were cut are then invisible.

What the lists leave out is the difference that is on paper only.  A
transition that carries no weight in the fit and no identification in IDEN2,
and that the classification does not accept, is saying the same thing in all
three files whichever of them happens to hold a row for it: it changes
neither the level optimisation, nor the next pass of classify_lines.py, nor
the view in IDEN2.  It arises in both directions - make_LOPT_input.py writes
every classified candidate, the rejected ones at zero weight, so settling a
classification after the transitions file was built leaves rejected
candidates in the workbook alone, and withdrawing one leaves a P-flagged row
in the transitions file alone - and in both directions it is counted in a
warning and left out of the lists.  The same transition weighted in the fit,
or identified in IDEN2, is a real disagreement and stays in the list.

Usage:

    python check_sync.py                     the report, and sync_report.txt
    python check_sync.py --quiet             only what is wrong
    python check_sync.py --list 40           show up to 40 items per finding
    python check_sync.py --out mine.txt      write the full report elsewhere
    python check_sync.py --drift 0.1         a stricter idea of "up to date"
    python check_sync.py --wn-tol 0.02       a stricter idea of "one line"
    python check_sync.py --dlv-tol 0.001    a stricter idea of one wavenumber
                                            printed twice, in dlv.dat
    python check_sync.py --csv sync.csv      the findings as a table too

It reads the files in the directory it is run from, falling back to the
directory holding the scripts for the ones that exist in only one place - the
same rule the swap tools use, so that it can be run inside an iteration folder
against that folder's own LOPT files.
"""

import argparse
import bisect
import csv
import os
import sys

import numpy as np
import pandas as pd

import config
import hfs_kappa
import output_files
import swap_line_assignments_IDEN as IDEN
import sync_IDEN2
from swap_paths import working_path


# ---------------------------------------------------------------------------
# What is built from what.  Only the modification times are used here; the
# contents are compared by the checks further down.
# ---------------------------------------------------------------------------
CLASSIFICATIONS_CSV = 'line_classifications.csv'
CLASSIFICATIONS_XLSX = 'line_classifications.xlsx'
LOPT_IN = 'LOPT_input_lines.txt'
LOPT_LINES = 'LOPT_output_lines.txt'
LOPT_LEVELS = 'LOPT_output_levels.txt'
DECISIONS = 'line_decisions.csv'
REVISIONS = 'revised_level_energies.csv'
UNSTABLE = 'unstable_candidates.csv'
POSITIONS = 'level_positions.csv'
LINES_XLSX = 'Pr3_lines.xlsx'
INFLATED = 'inflated_unc_lines.txt'
NAME_IDEN2 = 'IDEN2'
NAME_CONFIG = 'lineclass_config.toml'

# artefact -> the files it is built from.  A missing name is skipped.
DEPENDS = [
    (CLASSIFICATIONS_CSV, [LINES_XLSX, DECISIONS, REVISIONS]),
    (CLASSIFICATIONS_XLSX, [LINES_XLSX, DECISIONS, REVISIONS]),
    (LOPT_IN, [CLASSIFICATIONS_CSV]),
    (LOPT_LINES, [LOPT_IN]),
    (LOPT_LEVELS, [LOPT_IN]),
    (POSITIONS, [LOPT_LINES, LOPT_LEVELS, CLASSIFICATIONS_XLSX]),
]

# The columns of LOPT_input_lines.txt, as make_LOPT_input.py writes them:
# 1-based first and last column of each field.
IN_FIELDS = {
    'wavenumber': (1, 12),
    'uncertainty': (14, 19),
    'intensity': (25, 42),
    'lower_level': (43, 55),
    'upper_level': (59, 71),
    'flags': (73, 77),
    'weight': (82, 87),
}

# Where the complete report is written.  The terminal shows the first few
# items of each list; a list cut off after a dozen rows cannot be worked
# through, and the rows that were cut are then invisible.
REPORT = 'sync_report.txt'

DEF_DRIFT = 0.5
DEF_GROSS = 5.0
DEF_LIST = 12

# How far two files may put one transition's observed wavenumber apart before
# it counts as a different line rather than a difference of printing
# precision.  LOPT echoes each wavenumber to the precision its uncertainty
# warrants, which costs up to 0.03 cm^-1 on the widest lines.
DEF_WN_TOL = 0.05

# How far a wavenumber in IDEN2/dlv.dat may stand from the line list's own
# before the two are describing different numbers rather than one number
# printed twice.  Both carry three decimals and the roundings are made from
# different intermediate values, so a difference of one in the last digit is
# printing and nothing else; a wavelength calibration correction is tens of
# times larger than that.  An air row is allowed more where its wavelength's
# last printed digit is worth more (see _dlv_allowance), because IDEN2
# rebuilds the wavenumbers of the air rows from their wavelengths.  A row
# further off than this but within sync_IDEN2.DLV_MATCH is still its line,
# and is reported as having drifted from it.
DEF_DLV_TOL = 0.0015

# dlv.dat states the uncertainty of a line as a wavelength uncertainty in
# angstroms, to four decimals, so the coarsest it can state a wavenumber
# uncertainty is half of that converted, and it is allowed this fraction of
# the value on top: a hundredth of an uncertainty changes no judgment made by
# eye, and the line list's own column is itself a rounded number.
DLV_U_REL = 0.02

# classify_lines.py finds the line a ledger row rules on within this distance
# (its DECISIONS_WN_MATCH), so a ledger wavenumber that misses by less is
# obeyed and is not a mismatch.
DECISIONS_WN_MATCH = 0.01

ERROR, WARN, OK = 'ERROR', 'WARN', 'ok'
RANK = {ERROR: 0, WARN: 1, OK: 2}


class Report(object):
    """The findings, in the order they were made, and the actions they imply.

    A finding is a severity, the name of the check, one line of English, and
    the items it applies to - level identifiers, line wavenumbers, pairs.  The
    items are what makes a finding actionable, so they are kept whole and only
    abbreviated when printed.
    """

    def __init__(self, list_n=DEF_LIST):
        self.findings = []
        self.actions = []
        self.list_n = list_n

    def add(self, severity, check, message, items=None, head=None):
        self.findings.append(dict(severity=severity, check=check,
                                  message=message, items=list(items or []),
                                  head=head))

    def error(self, check, message, items=None, head=None):
        self.add(ERROR, check, message, items, head)

    def warn(self, check, message, items=None, head=None):
        self.add(WARN, check, message, items, head)

    def ok(self, check, message, items=None, head=None):
        self.add(OK, check, message, items, head)

    def act(self, priority, text):
        """Something to do.  Lower priority numbers are done first, because
        the pipeline runs one way and a fix applied out of order is undone by
        the next step."""
        self.actions.append((priority, text))

    @property
    def worst(self):
        if any(f['severity'] == ERROR for f in self.findings):
            return ERROR
        if any(f['severity'] == WARN for f in self.findings):
            return WARN
        return OK

    def print(self, quiet=False, out=None, full=False):
        """Write the findings.

        ``full`` writes every item of every list, which is what the
        report file is for: the terminal shows the first few rows so
        that the shape of the trouble is visible at once, and the file
        holds the whole of it so that it can be worked through.
        """
        out = out or sys.stdout
        limit = None if full else self.list_n
        section = None
        for f in self.findings:
            if quiet and f['severity'] == OK:
                continue
            if f['check'].split(':')[0] != section:
                section = f['check'].split(':')[0]
                out.write('\n' + section + '\n')
                out.write('-' * len(section) + '\n')
            tag = {ERROR: 'ERROR', WARN: 'warn ', OK: '  ok '}[f['severity']]
            out.write('  %s %s\n' % (tag, f['message']))
            shown = f['items'] if limit is None else f['items'][:limit]
            if shown and f['head']:
                out.write('         %s\n' % f['head'])
            for item in shown:
                out.write('         %s\n' % item)
            if limit is not None and len(f['items']) > limit:
                out.write('         ... and %d more; all of them are in %s\n'
                          % (len(f['items']) - limit, REPORT))

    def print_actions(self, out=None):
        out = out or sys.stdout
        out.write('\nWhat to do\n')
        out.write('----------\n')
        if not self.actions:
            out.write('  Nothing: the files agree.\n')
            return
        seen = set()
        n = 0
        for _, text in sorted(self.actions, key=lambda t: t[0]):
            if text in seen:
                continue
            seen.add(text)
            n += 1
            out.write('  %d. %s\n' % (n, text))
        out.write('\n  They are in the order the pipeline runs.  Doing a'
                  ' later\n  one first only means doing it again after the'
                  ' earlier one.\n')

    def to_frame(self):
        return pd.DataFrame([dict(severity=f['severity'], check=f['check'],
                                  message=f['message'],
                                  n_items=len(f['items']),
                                  items='; '.join(str(i) for i in f['items']))
                             for f in self.findings])


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------
def read_classifications(path):
    """line_classifications.csv, with the identifiers kept as text.

    Read as numbers, "059003.000456" becomes 59003.000456 and loses its
    leading zero and its last digit; every join in this script would then be
    empty and would report a total mismatch.
    """
    return pd.read_csv(path, low_memory=False,
                       dtype={'low_id': str, 'upp_id': str})


def classified(df):
    """The rows that name both levels: the ones make_LOPT_input.py writes."""
    lo = df['low_id'].fillna('').str.strip()
    up = df['upp_id'].fillna('').str.strip()
    return df[(lo != '') & (up != '')]


def pair_key(low, upp):
    """The identity of a candidate transition: the pair of levels it joins.

    Not the wavenumber as well, although that is the obvious key, because the
    files write it to different precisions: make_LOPT_input.py gives three
    decimals, LOPT echoes each line to the precision its own uncertainty
    warrants (121279.416 comes back as 121279.42), and the decision ledger is
    written by hand from a printed list.  Keying on the wavenumber therefore
    reports thousands of differences that are only rounding.  The pair is a
    key in its own right - one transition can belong to only one observed
    line, which check_duplicate_pairs() verifies rather than assumes - so the
    wavenumbers are compared afterwards, with a tolerance.
    """
    return (str(low).strip(), str(upp).strip())


def fmt_pair(pair, wn=None):
    if wn is None:
        return '%s - %s' % pair
    return '%12.3f  %s - %s' % (float(wn), pair[0], pair[1])


def align(a, b, tol):
    """Match two files' records transition by transition.

    Both arguments are {pair: [record, ...]} with the observed wavenumber
    first in every record.  Returns (matched, a_only, b_only), where matched
    is a list of (pair, record_a, record_b).

    Within one pair of levels the records are matched by wavenumber, nearest
    first, because the pair alone is not a key: the same two levels can be
    candidates on several observed lines, and moving an assignment from one
    line to another - which is how a correction is made - leaves the pair on
    both.  Matching on the wavenumber alone is no good either, since the files
    write it to different precisions; hence the pair, then the wavenumber
    within `tol`.
    """
    matched, a_only, b_only = [], [], []
    for pair in sorted(set(a) | set(b)):
        left = sorted(a.get(pair, []), key=lambda r: r[0])
        right = list(sorted(b.get(pair, []), key=lambda r: r[0]))
        for rec in left:
            best, best_d = None, tol
            for other in right:
                d = abs(float(rec[0]) - float(other[0]))
                if d <= best_d:
                    best, best_d = other, d
            if best is None:
                a_only.append((pair, rec))
            else:
                right.remove(best)
                matched.append((pair, rec, best))
        for rec in right:
            b_only.append((pair, rec))
    return matched, a_only, b_only


class Context(object):
    """Where each transition stands in every file, and how to print it.

    Every list of problem transitions answers the same three questions, and
    the answers are what turn a list into something that can be worked
    through at the screen.  Each is asked about the transition ON THE LINE
    the row names, not about the pair of levels in general: one pair can be a
    candidate on several observed lines, and reading the pair's status off
    another line is how a row that is in perfect order comes to be printed as
    a problem.

    * ``.xlsx``  - what line_classifications.xlsx says: ``accepted``, or
      ``rejected`` (the pair is a candidate on this line and was turned
      down), or ``elsewhere`` (it is a candidate, but on other lines only),
      or ``missing`` (the pair is not a candidate anywhere).
    * ``LOPT``   - whether the transition carries weight in the fit.
      ``included`` means it is in LOPT_input_lines.txt on this line with a
      weight above zero and no ``P`` flag; ``not incl. P`` means it is in the
      file but flagged predicted or given zero weight; ``elsewhere`` means
      the file has the pair on other lines only, and ``not incl. --`` that it
      is not in the file at all.  All but the first are the same thing to
      LOPT: this line does not pull on the levels.
    * ``IDEN2``  - what trans.dat says: ``assigned`` (the transition is
      identified with this line), ``elsewhere`` (identified with another
      line), ``not assigned`` (IDEN2 predicts it and has no line for it), or
      ``no trans row`` (IDEN2 does not predict it at all, which for a real
      pair of levels means trans.dat is stale).

    The two levels are printed as one field, lower first, in the form a
    transition is written and searched for; their energies follow in the same
    order, and then IDEN2's own level numbers - the numbers typed into IDEN2 -
    so that a row can be looked up without a search through
    IDEN_level_ids.txt.  IDEN_id1 is the lower level's number and IDEN_id2 the
    upper level's.  The rows are still sorted by upper energy and then by
    lower energy, which is the order the same transitions are met in IDEN2.
    """

    FMT = '%-10s %-27s %10s %10s  %-8s %-8s %-9s %-12s %-12s'
    HEAD = FMT % ('wn_obs', 'lower_id - upper_id', 'E_lower', 'E_upper',
                  'IDEN_id1', 'IDEN_id2', '.xlsx', 'LOPT', 'IDEN2')

    def __init__(self, cls=None, levels=None, lopt_in=None, iden=None,
                 wn_tol=DEF_WN_TOL):
        self.wn_tol = wn_tol
        self.levels = dict(levels or {})
        self.iden = iden
        # The transitions file as it was read, keyed the way the file writes
        # the pair, for align(); and again keyed on the sorted pair, for the
        # questions the item lists ask.
        self.lopt_in = lopt_in
        self.lopt = {}
        for pair, recs in (lopt_in or {}).items():
            self.lopt.setdefault(tuple(sorted(pair)), []).extend(recs)
        self.rows = {}
        if cls is not None:
            for pair, recs in classified_rows(classified(cls)).items():
                self.rows.setdefault(tuple(sorted(pair)), []).extend(recs)
        # The transitions some earlier check has already printed in full,
        # as {sorted pair: [wavenumber, ...]}.  A list is worth reading only
        # for what it adds: the same transition met again for the same reason
        # further down the report is noise, and the count in the message says
        # it is there without spending a line on it.  See listed() and
        # already_listed().
        self.shown = {}
        self.iden_row, self.assigned, self.predicted = {}, {}, set()
        if iden is not None:
            self.iden_row = dict(iden.row_of)
            self.assigned, _ = read_iden_assignments(
                iden.enlev, iden.trans, iden.id_of_row, iden.wn_of_line)
            for owner, partner in iden.trans.row_of:
                a = iden.id_of_row.get(owner)
                b = iden.id_of_row.get(partner)
                if a and b:
                    self.predicted.add(tuple(sorted((a, b))))

    # -- the level ------------------------------------------------------
    def energy(self, lid):
        """The fitted energy, or IDEN2's if the fit does not know the level.

        A level absent from both gives a nan, which sorts last and prints as
        a question mark rather than as a number that is not there.
        """
        if lid in self.levels:
            return float(self.levels[lid])
        if self.iden is not None and lid in self.iden_row:
            try:
                return float(self.iden.enlev.e_obs(self.iden_row[lid]))
            except (ValueError, KeyError):
                pass
        return float('nan')

    def iden_of(self, lid):
        n = self.iden_row.get(lid)
        return '-' if n is None else str(n)

    def order(self, pair):
        """(upper, lower), decided by energy rather than by identifier."""
        a, b = tuple(pair)
        ea, eb = self.energy(a), self.energy(b)
        if ea == ea and eb == eb and ea != eb:
            return (a, b) if ea > eb else (b, a)
        return (max(a, b), min(a, b))

    def sort_key(self, pair, wn=None):
        """Upper energy first, then lower energy: the order the same
        transitions are looked at in IDEN2 and in the workbook."""
        upp, low = self.order(pair)
        eu, el = self.energy(upp), self.energy(low)
        far = 1e12
        return (eu if eu == eu else far, el if el == el else far,
                float(wn) if wn is not None else 0.0)

    # -- the transition -------------------------------------------------
    def _here(self, wn, other):
        """Is that record's wavenumber the line being asked about?

        With no wavenumber to ask about - a question put about the pair of
        levels alone - every record counts.
        """
        return (wn is None or other is None
                or abs(float(wn) - float(other)) <= self.wn_tol)

    def in_xlsx(self, pair, wn=None):
        """What the classification says about this transition on this line.

        One pair of levels can be a candidate on several observed lines, and
        that is not an accident: moving an assignment from one line to
        another leaves the old row in place to be rejected, which is how the
        pipeline is made to turn down an assignment it would otherwise keep.
        The verdict therefore belongs to the line, not to the pair, and a
        pair classified only on other lines is `elsewhere` - nothing is wrong
        with it here, it simply belongs to a different observed line.
        """
        recs = self.rows.get(tuple(sorted(pair)), [])
        if not recs:
            return 'missing'
        here = [acc for w, acc in recs if self._here(wn, w)]
        if not here:
            return 'elsewhere'
        return 'accepted' if any(here) else 'rejected'

    def in_lopt(self, pair, wn=None):
        recs = [r for r in self.lopt.get(tuple(sorted(pair)), [])
                if self._here(wn, r[0])]
        if not recs:
            return ('elsewhere' if self.lopt.get(tuple(sorted(pair)))
                    else 'not incl. --')
        if wn is not None:
            recs.sort(key=lambda r: abs(float(r[0]) - float(wn)))
        return 'included' if _weighted(recs[0]) else 'not incl. P'

    def in_iden(self, pair, wn=None):
        if self.iden is None:
            return '?'
        key = tuple(sorted(pair))
        got = self.assigned.get(key)
        if got is not None:
            return 'assigned' if self._here(wn, got) else 'elsewhere'
        return 'not assigned' if key in self.predicted else 'no trans row'

    def item(self, pair, wn=None, note=''):
        """One line of a problem list, with everything needed to look it
        up: the observed wavenumber, both levels with IDEN2's numbers and
        their energies, and where the transition stands in the three
        files."""
        upp, low = self.order(pair)
        text = self.FMT % (
            '' if wn is None else '%.3f' % float(wn),
            '%s-%s' % (low, upp),
            self._e_text(low), self._e_text(upp),
            self.iden_of(low), self.iden_of(upp),
            self.in_xlsx(pair, wn), self.in_lopt(pair, wn),
            self.in_iden(pair, wn))
        return (text + '  ' + note).rstrip()

    def listed(self, pairs_wn):
        """Remember that these transitions have been printed."""
        for pair, wn in pairs_wn:
            self.shown.setdefault(tuple(sorted(pair)), []).append(
                None if wn is None else float(wn))

    def already_listed(self, pair, wn=None):
        """Has this transition, on this line, been printed already?"""
        return any(self._here(wn, other)
                   for other in self.shown.get(tuple(sorted(pair)), []))

    def items(self, pairs_wn, note_of=None, by_wn=False):
        """A whole list, sorted by upper energy and then by lower energy.

        With by_wn the rows are sorted by decreasing wavenumber instead.  That
        is the order for a list about the observed lines themselves rather
        than about the levels - a wrongly divided blend is a fault of one
        observed line, and the strongest end of the spectrum is where the
        blends are, so the lines worth looking at first come first.
        """
        if by_wn:
            def order(t):
                wn = t[1]
                return (0, -float(wn)) if wn is not None else (1, 0.0)
        else:
            def order(t):
                return self.sort_key(t[0], t[1])
        rows = sorted(pairs_wn, key=order)
        return [self.item(p, wn, (note_of or {}).get((p, wn), ''))
                for p, wn in rows]

    def _e_text(self, lid):
        e = self.energy(lid)
        return '?' if e != e else '%.3f' % e


def read_lopt_input(path):
    """LOPT_input_lines.txt as {pair: (wavenumber, flag, weight)}, by columns.

    The fields cannot be found by splitting on blanks: an intensity of 128400
    fills its column and runs into the next one, so a whitespace split silently
    merges two fields on exactly the brightest lines.
    """
    def cut(rec, name):
        first, last = IN_FIELDS[name]
        return rec[first - 1:last]

    out = {}
    with open(path, 'r', encoding='latin-1', newline='') as fh:
        for raw in fh:
            rec = raw.rstrip('\r\n')
            if not rec.strip():
                continue
            wn = cut(rec, 'wavenumber').strip()
            low = cut(rec, 'lower_level').strip()
            upp = cut(rec, 'upper_level').strip()
            if not wn or not low or not upp:
                continue
            weight = cut(rec, 'weight').strip()
            try:
                value = (float(wn), cut(rec, 'flags').strip(),
                         float(weight) if weight else 0.0)
            except ValueError:
                continue
            out.setdefault(pair_key(low, upp), []).append(value)
    return out


def read_lopt_output_lines(path):
    """LOPT_output_lines.txt as {pair: (wavenumber, weight, E_low, E_upp)}."""
    out = {}
    with open(path, 'r', encoding='latin-1', newline='') as fh:
        rd = csv.reader(fh, delimiter='\t')
        header = next(rd)
        ix = {h.strip(): i for i, h in enumerate(header)}
        for need in ('wn_o', 'L1', 'L2', 'E1', 'E2', 'Weight'):
            if need not in ix:
                raise SystemExit(
                    '%s has no %r column - switch its printing on in the LOPT '
                    'parameter file' % (os.path.basename(path), need))
        for r in rd:
            if len(r) <= ix['Weight']:
                continue
            low, upp = r[ix['L1']].strip(), r[ix['L2']].strip()
            if not low or not upp:
                continue
            out.setdefault(pair_key(low, upp), []).append(
                (_num(r[ix['wn_o']]), _num(r[ix['Weight']]),
                 _num(r[ix['E1']]), _num(r[ix['E2']])))
    return out


def _num(text):
    text = (text or '').strip()
    if text in ('', '_'):
        return np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


def read_lopt_levels(path):
    """{level_id: energy} from LOPT_output_levels.txt."""
    df = pd.read_csv(path, sep='\t', dtype={'Designation': str})
    return dict(zip(df['Designation'].str.strip(),
                    df['Energy'].astype(float)))


def classified_rows(df, column='wn_obs'):
    """{pair: [(wavenumber, accepted), ...]} for the rows naming both levels.

    A pair may legitimately appear more than once - the same two levels can be
    candidates on several observed lines, and only one of those rows is
    accepted.  Two ACCEPTED rows for one pair would be a contradiction, since
    one transition can be on one line only; check_double_acceptance() looks
    for that.

    `column` is the wavenumber the rows are listed under: `wn_obs` for the
    wavenumber the set works with, `wn_key` for Sugar's value, which is
    what the ledger names a line by.
    """
    out = {}
    for wn, low, upp, acc in zip(df[column], df['low_id'], df['upp_id'],
                                 df['accepted']):
        out.setdefault(pair_key(low, upp), []).append(
            (float(wn), float(acc or 0) == 1.0))
    return out


def check_double_acceptance(rows, rep, list_n):
    """One transition on two observed lines at once.

    Nothing downstream can represent it: the transitions file would carry the
    pair twice, LOPT would fit it to both lines, and trans.dat has one row per
    pair and would keep whichever came last.
    """
    bad = ['%s on %s' % (fmt_pair(pair),
                         ' and '.join('%.3f' % w for w, a in recs if a))
           for pair, recs in sorted(rows.items())
           if sum(1 for _, a in recs if a) > 1]
    if bad:
        rep.error('Classification',
                  '%d transitions are accepted on more than one observed line'
                  % len(bad), bad)
        rep.act(10, 'Decide which line each of those %d transitions belongs '
                    'to and record it in %s.' % (len(bad), DECISIONS))
    return bad


def read_iden_assignments(enlev, trans, id_of_row, wn_of_line=None):
    """{(low_id, upp_id) sorted: wavenumber} for every identified transition
    of trans.dat, and the rows whose levels are not in the map.

    A trans.dat transition row carries, after its predicted wavenumber, the
    observed line it has been identified with: a wavenumber, the
    observed-minus-predicted, and the line number of the line in dlv.dat.  A
    line number of zero means the transition has not been identified.  IDEN2
    reads only the line number, so the wavenumber is taken from dlv.dat,
    ``wn_of_line`` ({line number: wavenumber}); the one in trans.dat is used
    only for a line number dlv.dat has not got, or when there is no dlv.dat.
    """
    wn_of_line = wn_of_line or {}
    found, unmapped = {}, []
    for (owner, partner), k in trans.row_of.items():
        obs = IDEN.assignment(trans.records[k])
        if not IDEN.has_line(obs):
            continue
        wn = wn_of_line.get(IDEN.obs_row(obs))
        if wn is None:
            wn = IDEN.obs_wavenumber(obs)
        if owner not in id_of_row or partner not in id_of_row:
            unmapped.append((owner, partner, wn))
            continue
        pair = tuple(sorted((id_of_row[owner], id_of_row[partner])))
        found[pair] = wn
    return found, unmapped


# ---------------------------------------------------------------------------
# 1. Freshness
# ---------------------------------------------------------------------------
def check_freshness(paths, rep):
    """An artefact older than something it was built from is out of date.

    This is a statement about the files, not about their contents: a rebuild
    that changed nothing still moves the timestamp, so a stale timestamp is
    only a reason to look, and the checks that follow say whether the
    staleness matters.  It is put first because it is the one check that can
    explain all the others at once.
    """
    stale = False
    for name, sources in DEPENDS:
        path = paths.get(name)
        if not path or not os.path.exists(path):
            continue
        t = os.path.getmtime(path)
        older = []
        for src in sources:
            spath = paths.get(src)
            if spath and os.path.exists(spath) and os.path.getmtime(spath) > t:
                older.append('%s (%s newer by %s)'
                             % (src, os.path.basename(spath),
                                _age(os.path.getmtime(spath) - t)))
        if older:
            stale = True
            rep.warn('Freshness', '%s is older than what it is built from:'
                     % name, older)
    if not stale:
        rep.ok('Freshness', 'every artefact is newer than its inputs')


def _age(seconds):
    if seconds < 90:
        return '%.0f s' % seconds
    if seconds < 5400:
        return '%.0f min' % (seconds / 60.0)
    if seconds < 172800:
        return '%.1f h' % (seconds / 3600.0)
    return '%.1f days' % (seconds / 86400.0)


# ---------------------------------------------------------------------------
# 2. The two copies of the classification
# ---------------------------------------------------------------------------
def check_csv_xlsx(csv_path, xlsx_path, rep, list_n):
    """classify_lines.py writes the csv and the workbook from one table, so
    any difference between them was made afterwards - by an edit in Excel, or
    by one of the two being left over from an earlier run."""
    if not (csv_path and os.path.exists(csv_path)):
        rep.error('Classification', '%s is missing' % CLASSIFICATIONS_CSV)
        return None, None
    a = read_classifications(csv_path)
    if not (xlsx_path and os.path.exists(xlsx_path)):
        rep.warn('Classification', '%s is missing; the csv is used alone'
                 % CLASSIFICATIONS_XLSX)
        return a, None
    b = pd.read_excel(xlsx_path, dtype={'low_id': str, 'upp_id': str})
    if len(a) != len(b):
        rep.error('Classification',
                  'the csv has %d rows and the workbook %d - they are from '
                  'different runs' % (len(a), len(b)))
        rep.act(20, 'Re-run classify_lines.py so that %s and %s are written '
                    'together.' % (CLASSIFICATIONS_CSV, CLASSIFICATIONS_XLSX))
        return a, b

    pa = classified_rows(classified(a))
    pb = classified_rows(classified(b))
    check_double_acceptance(pa, rep, list_n)

    matched, a_only, b_only = align(pa, pb, 0.001)
    if a_only or b_only:
        rep.error('Classification',
                  '%d classified transitions are in the csv only and %d in '
                  'the workbook only' % (len(a_only), len(b_only)),
                  [('csv only: ' + fmt_pair(p, r[0])) for p, r in a_only[:20]]
                  + [('workbook only: ' + fmt_pair(p, r[0]))
                     for p, r in b_only[:20]])
        rep.act(20, 'Re-run classify_lines.py; the csv and the workbook hold '
                    'different assignments.')
        return a, b

    bad = [fmt_pair(p, ra[0]) for p, ra, rb in matched if ra[1] != rb[1]]
    if bad:
        rep.error('Classification',
                  '%d transitions are accepted in one copy and not the other'
                  % len(bad), bad)
        rep.act(20, 'Re-run classify_lines.py; the two copies disagree about '
                    'which assignments are accepted.')
    else:
        rep.ok('Classification',
               'the csv and the workbook hold the same %d classified '
               'transitions, on the same lines and with the same verdicts'
               % len(matched))
    return a, b


# ---------------------------------------------------------------------------
# 3-5. The LOPT chain
# ---------------------------------------------------------------------------
def _weighted(rec):
    """Does this row of LOPT_input_lines.txt pull on the levels?

    Only if it has a weight above zero and is not flagged P.  The two are
    separate fields and either one on its own takes the line out of the fit,
    so both have to be read: a P-flagged row often carries its old weight,
    and reading the weight alone would call it fitted.
    """
    return rec[2] > 0 and 'P' not in (rec[1] or '').upper()


def _split_idle(missing, ctx):
    """Split the classified-but-absent transitions into the ones that matter
    and the ones that do not.

    make_LOPT_input.py writes every classified candidate, the rejected ones
    with zero weight, so a classification settled after the transitions file
    was built leaves rejected candidates in the workbook and nowhere else.
    That difference is on paper only: a rejected candidate that carries no
    weight in the fit and no identification in IDEN2 is saying the same thing
    in all three files, and rebuilding LOPT_input_lines.txt at the end of the
    analysis - when the fully synchronised files are wanted for the
    publication tables - is what puts it in writing.  Such a transition is
    counted, not listed.

    A rejected candidate that IS weighted in the fit, or IS identified in
    IDEN2, is a different matter and stays in the list: there the files
    disagree about the identification itself.

    Returns (real, idle).
    """
    real, idle = [], []
    for pair, rec in missing:
        wn, accepted = rec[0], rec[1]
        if (not accepted and ctx.in_lopt(pair, wn) != 'included'
                and ctx.in_iden(pair, wn) != 'assigned'):
            idle.append((pair, rec))
        else:
            real.append((pair, rec))
    return real, idle


def _split_dead(extra, ctx):
    """The same test the other way round, for the transitions the file holds
    and the classification does not.

    A row of LOPT_input_lines.txt flagged P or given zero weight is not in
    the fit; if the classification does not accept it and IDEN2 does not
    identify it, it is a leftover of an earlier classification that is doing
    nothing anywhere, and the next rebuild of the file will drop it.  A row
    that IS weighted, or IS identified in IDEN2, is a real disagreement and
    stays in the list.

    Returns (real, dead).
    """
    real, dead = [], []
    for pair, rec in extra:
        wn = rec[0]
        if (not _weighted(rec) and ctx.in_xlsx(pair, wn) != 'accepted'
                and ctx.in_iden(pair, wn) != 'assigned'):
            dead.append((pair, rec))
        else:
            real.append((pair, rec))
    return real, dead


def expected_shares(df):
    """{(pair, wavenumber): the share of its observed line each accepted
    classification should carry in LOPT_input_lines.txt}.

    This is make_LOPT_input.py's own arithmetic, repeated so that a weight
    written by hand can be compared with the one the pipeline would write: a
    line with a single accepted classification gives it the whole weight, and
    the components of a blend divide it in proportion to their calculated
    intensities, sharing it equally if those intensities say nothing.
    """
    groups = {}
    for wn, low, upp, acc, calc in zip(df['wn_obs'], df['low_id'],
                                       df['upp_id'], df['accepted'],
                                       df['calc_intens']):
        if float(acc or 0) != 1.0:
            continue
        c = float(calc) if calc == calc else 0.0
        groups.setdefault(round(float(wn), 3), []).append(
            (tuple(sorted(pair_key(low, upp))), float(wn), max(c, 0.0)))
    out = {}
    for rows in groups.values():
        if len(rows) == 1:
            pair, wn, _c = rows[0]
            out[(pair, round(wn, 3))] = 1.0
            continue
        total = sum(c for _p, _w, c in rows)
        if total <= 0:
            rows = [(p, w, 1.0) for p, w, _c in rows]
            total = float(len(rows))
        for pair, wn, c in rows:
            out[(pair, round(wn, 3))] = c / total
    return out


def line_weights(lopt_in):
    """{wavenumber: the sum of the weights on that observed line}.

    LOPT divides each weight by this sum, so what a row is worth in the fit
    is its share of its line and never the number in the column.
    """
    total = {}
    for _pair, recs in lopt_in.items():
        for rec in recs:
            if _weighted(rec):
                key = round(float(rec[0]), 3)
                total[key] = total.get(key, 0.0) + float(rec[2])
    return total


def same_share(want, got):
    """Do the two files give a transition the same share of its line?

    Not the same number: a blend assigned by hand carries the calculated
    intensities pasted straight out of Icalc.xlsx, which do not sum to one,
    and LOPT normalises them itself, so the columns are expected to look
    unlike each other and to fit identically.  What has to agree is the
    share, to one per cent.  The 1e-4 is the rounding of the four-decimal
    weight column, which is the whole of the difference for the small
    component of a lopsided blend.
    """
    d = abs(want - got)
    return d <= 1e-4 or d <= 0.01 * max(want, got)


def one_record_per_transition(path):
    """The transitions of the file that are written more than once.

    Returns `{(lower, upper): [(wavenumber, flags)]}` for the pairs of levels
    that carry several records.  A transition is one energy difference, so a
    second record of it adds nothing to the fit; and if two of them are
    weighted, LOPT fits the pair once per record and splits the evidence
    between them.
    """
    def cut(rec, name):
        first, last = IN_FIELDS[name]
        return rec[first - 1:last].strip()

    pairs = {}
    with open(path, 'r', encoding='latin-1', newline='') as fh:
        for raw in fh:
            rec = raw.rstrip('\r\n')
            if not rec.strip():
                continue
            low, upp = cut(rec, 'lower_level'), cut(rec, 'upper_level')
            if not low or not upp:
                continue
            pairs.setdefault((low, upp), []).append(
                (cut(rec, 'wavenumber'), cut(rec, 'flags')))
    return {k: v for k, v in pairs.items() if len(v) > 1}


def one_uncertainty_per_line(path):
    """The observed lines of the transitions file that carry two uncertainties.

    Returns `{wavenumber: [(uncertainty, flags, lower, upper)]}` for the
    wavenumbers whose records disagree.  One observed line is one measurement,
    so LOPT has to be given one uncertainty for it however many transitions it
    is assigned to; two values on one wavenumber are two weights on one
    measured position, which is a statement the measurement cannot make.  The
    file is read here by columns rather than through `read_lopt_input`, which
    keys its rows by the pair of levels and does not carry the uncertainty.
    """
    def cut(rec, name):
        first, last = IN_FIELDS[name]
        return rec[first - 1:last].strip()

    groups = {}
    with open(path, 'r', encoding='latin-1', newline='') as fh:
        for raw in fh:
            rec = raw.rstrip('\r\n')
            if not rec.strip():
                continue
            wn, unc = cut(rec, 'wavenumber'), cut(rec, 'uncertainty')
            if not wn or not unc:
                continue
            groups.setdefault(wn, []).append(
                (unc, cut(rec, 'flags'),
                 cut(rec, 'lower_level'), cut(rec, 'upper_level')))
    return {wn: recs for wn, recs in groups.items()
            if len({r[0] for r in recs}) > 1}


def check_lopt_chain(cls, paths, ctx, rep, wn_tol):
    """The classification, the transitions file LOPT was given, and the two
    files LOPT wrote, are four views of one fit; each is made from the one
    before it by a command that has to be run by hand."""
    inp = paths.get(LOPT_IN)
    out = paths.get(LOPT_LINES)
    lev = paths.get(LOPT_LEVELS)
    lopt_in = ctx.lopt_in
    levels = ctx.levels

    if inp:
        split = one_uncertainty_per_line(inp)
        if split:
            fitted = sum(1 for recs in split.values()
                         if len({r[0] for r in recs
                                 if 'P' not in (r[1] or '').upper()}) > 1)
            items = []
            for wn in sorted(split, key=lambda t: -float(t)):
                for unc, flag, low, upp in split[wn]:
                    items.append('%12s  %7s  %-2s %s - %s'
                                 % (wn, unc, flag, low, upp))
            rep.warn('LOPT',
                     '%d observed lines in %s carry more than one '
                     'uncertainty, %d of them among the records the fit '
                     'uses. LOPT reads the uncertainty of each record as the '
                     'uncertainty of the measurement it constrains, so those '
                     'lines are weighted differently in the components of one '
                     'blend. Rebuilding the file cures it: make_LOPT_input.py '
                     'gives an observed line one uncertainty, the mean of its '
                     'transitions\' values weighted by the weights LOPT is '
                     'given.' % (len(split), LOPT_IN, fitted), items)
        else:
            rep.ok('LOPT', 'every observed line in %s carries one '
                           'uncertainty' % LOPT_IN)

        twice = one_record_per_transition(inp)
        if twice:
            weighted = sum(1 for recs in twice.values()
                           if sum(1 for _, f in recs
                                  if 'P' not in (f or '').upper()) > 1)
            items = []
            for (low, upp) in sorted(twice):
                for wn, flag in twice[(low, upp)]:
                    items.append('%s - %s  %12s  %s' % (low, upp, wn, flag))
            rep.warn('LOPT',
                     '%d transitions in %s are written more than once, %d of '
                     'them with more than one weighted record. A transition '
                     'is one energy difference, so a second record of it adds '
                     'nothing to the fit, and two weighted records make LOPT '
                     'fit the pair twice and split its evidence. Rebuilding '
                     'the file cures it: make_LOPT_input.py keeps the '
                     'accepted record of a transition and drops the rest.'
                     % (len(twice), LOPT_IN, weighted), items)
        else:
            rep.ok('LOPT', 'every transition in %s is written once' % LOPT_IN)

    if cls is not None and lopt_in is not None:
        want = classified_rows(classified(cls))
        matched, missing, extra = align(want, lopt_in, wn_tol)
        shares = expected_shares(classified(cls))
        total = line_weights(lopt_in)
        wrong, note = [], {}
        # A row that both files hold can still say opposite things.  Being in
        # the transitions file is not being in the fit - a P flag or a zero
        # weight takes it out - so an accepted classification that is in the
        # file without weight belongs with the ones that are not in the file
        # at all, and a weighted row the classification does not accept
        # belongs with the ones the classification has never heard of.  Only
        # when both files put the line in the fit is there a weight to
        # compare.
        for pair, ra, rb in matched:
            fitted = _weighted(rb)
            if ra[1] and not fitted:
                missing.append((pair, ra))
            elif fitted and not ra[1]:
                extra.append((pair, rb))
            elif fitted:
                key = round(float(rb[0]), 3)
                got = float(rb[2]) / total[key] if total.get(key) else 0.0
                due = shares.get((tuple(sorted(pair)),
                                  round(float(ra[0]), 3)), 1.0)
                if not same_share(due, got):
                    wrong.append((pair, ra))
                    note[(pair, ra[0])] = ('%.4f of the line here, %.4f in %s'
                                           % (due, got, LOPT_IN))
        missing, idle = _split_idle(missing, ctx)
        extra, dead = _split_dead(extra, ctx)
        if idle:
            rep.warn('LOPT',
                     '%d rejected candidates are in the classification and '
                     'nowhere else: absent from %s and unassigned in IDEN2. '
                     'They are counted, not listed - carrying no weight and '
                     'no line, they say the same thing in all three files, '
                     'and the rebuild that puts them on paper belongs at the '
                     'end of the analysis.' % (len(idle), LOPT_IN))
        if dead:
            rep.warn('LOPT',
                     '%d transitions in %s are flagged P or of zero weight, '
                     'are not accepted in the classification, and are not '
                     'identified in IDEN2: leftovers of an earlier '
                     'classification that pull on nothing.  They are counted, '
                     'not listed; the next rebuild of the file drops them.'
                     % (len(dead), LOPT_IN))
        if missing or extra or wrong:
            # The last flag is the sort: the first two groups are about
            # levels and keep the level order, the weights are about the
            # observed lines and run from the highest wavenumber down.
            groups = [
                ('in %s, not in %s:' % (CLASSIFICATIONS_XLSX, LOPT_IN),
                 missing, False),
                ('in %s, not in %s:' % (LOPT_IN, CLASSIFICATIONS_XLSX),
                 extra, False),
                ('the share of the observed line differs by more than '
                 '1 per cent:', wrong, True),
            ]
            items = []
            for label, rows, by_wn in groups:
                if not rows:
                    continue
                items.append(label)
                items.extend(ctx.items([(p, r[0]) for p, r in rows], note,
                                       by_wn=by_wn))
            # What the IDEN2 check finds next is largely these same
            # transitions seen from the other end - an assignment absent from
            # the fit is usually the one IDEN2 still shows - so tell it what
            # has been printed here.  The weights are deliberately included:
            # a wrongly weighted blend is a different fault from an unmade
            # identification, and it is worth meeting twice.
            ctx.listed([(p, r[0]) for _label, rows, _by_wn in groups
                        for p, r in rows])
            said = [n for n in
                    ('%d accepted transitions are not in the fit it makes'
                     % len(missing) if missing else '',
                     '%d are weighted in it and not accepted in the '
                     'classification' % len(extra) if extra else '',
                     '%d are weighted differently' % len(wrong) if wrong
                     else '') if n]
            rep.error('LOPT', '%s does not match the classification: %s'
                      % (LOPT_IN, ', '.join(said)), items, head=Context.HEAD)
            rep.act(30, 'Look at each of those %d transitions in IDEN2 and '
                        'settle it, recording the verdict in %s; where it is '
                        'the share of a blend that differs, correct the '
                        'weight in %s by hand.  Do not simply re-run '
                        'make_LOPT_input.py: that rebuilds %s from the '
                        'current classification, so the list above '
                        'disappears while every question in it stays open '
                        'and is now untracked.'
                        % (len(missing) + len(extra) + len(wrong),
                           DECISIONS, LOPT_IN, LOPT_IN))
        else:
            rep.ok('LOPT', '%s holds every classified transition that bears '
                           'on the fit (%d of them), weighted as accepted'
                   % (LOPT_IN, len(matched)))

    lopt_out = None
    if out and os.path.exists(out):
        lopt_out = read_lopt_output_lines(out)
        if lopt_in is not None:
            matched, missing, extra = align(lopt_in, lopt_out, wn_tol)
            if missing or extra:
                rep.error('LOPT',
                          'LOPT was run on a different transitions file: %d '
                          'lines of %s are not in the output, %d lines of the '
                          'output are not in %s'
                          % (len(missing), LOPT_IN, len(extra), LOPT_IN),
                          [fmt_pair(p, r[0])
                           for p, r in (missing + extra)[:20]])
                rep.act(40, 'Re-run lopt.bat on the current %s.' % LOPT_IN)
            else:
                rep.ok('LOPT', 'the fit was made on the transitions file that '
                               'is on disk now (%d lines)' % len(matched))

    if lev and os.path.exists(lev):
        if lopt_out is not None:
            bad = []
            for (low, upp), recs in lopt_out.items():
                _wn, _w, e1, e2 = recs[0]
                for lid, e in ((low, e1), (upp, e2)):
                    if lid in levels and not np.isnan(e) \
                            and abs(levels[lid] - e) > 0.5:
                        bad.append('%s  %.3f in the line output, %.3f in the '
                                   'level output' % (lid, e, levels[lid]))
            bad = sorted(set(bad))
            if bad:
                rep.error('LOPT',
                          'the line output and the level output give %d '
                          'levels different energies - they are from '
                          'different runs' % len(bad), bad)
                rep.act(40, 'Re-run lopt.bat; %s and %s are from different '
                            'runs.' % (LOPT_LINES, LOPT_LEVELS))
            else:
                rep.ok('LOPT', 'the line output and the level output agree '
                               'about every level energy')
    return lopt_out


def check_revisions(paths, levels, rep, gross):
    """A level revision is recorded in revised_level_energies.csv and carried
    out by re-running the pipeline.  A row whose energy the fit does not show
    is a revision that was written down and never made - or, just as often,
    one that was made and then undone by a later run started from an older
    classification."""
    path = paths.get(REVISIONS)
    if not (path and os.path.exists(path)) or not levels:
        return
    rev = pd.read_csv(path, dtype={'level_id': str})
    if 'E_input' not in rev.columns:
        return
    late, absent = [], []
    for lid, e in zip(rev['level_id'].str.strip(), rev['E_input']):
        if lid not in levels:
            absent.append(lid)
            continue
        d = float(levels[lid]) - float(e)
        if abs(d) >= gross:
            late.append('%s  ledger %.3f, fit %.3f  (%+.3f)'
                        % (lid, float(e), levels[lid], d))
    if absent:
        rep.error('Revisions', '%d revised levels are not in %s'
                  % (len(absent), LOPT_LEVELS), absent)
    if late:
        rep.error('Revisions',
                  '%d recorded revisions are not in the fit - the energy the '
                  'ledger gives and the energy LOPT settled on differ by more '
                  'than %g cm^-1' % (len(late), gross), late)
        rep.act(20, 'Re-run classify_lines.py and then the LOPT step: %s '
                    'records level energies the fit does not have.'
                    % REVISIONS)
    elif not absent:
        rep.ok('Revisions', 'all %d recorded revisions are in the fit'
               % len(rev))


# ---------------------------------------------------------------------------
# 7-9. IDEN2
# ---------------------------------------------------------------------------
class Iden2(object):
    """The IDEN2 files, read once.

    They are read before any check is made rather than inside the IDEN2
    check, because every problem list - including the LOPT ones - names
    IDEN2's level numbers and says whether the transition is assigned there.

    ``dlv`` is the records of ``dlv.dat``, the observed line list, or None
    when the set has not got the file.  It is the one file of the three that
    IDEN2 never writes: it is the measurement, and every assignment in
    ``trans.dat`` names a line by the number of its row in it.
    """

    def __init__(self, enlev, trans, row_of, dlv=None):
        self.enlev = enlev
        self.trans = trans
        self.row_of = row_of
        self.id_of_row = {n: lid for lid, n in row_of.items()}
        self.dlv = dlv
        # {line number: wavenumber} - the wavenumber an identification in
        # trans.dat stands for.  A malformed row is left out here and
        # reported by check_dlv.
        self.wn_of_line = {}
        for rec in dlv or []:
            try:
                self.wn_of_line[int(rec[_span(sync_IDEN2.DLV_ROW)])] = \
                    float(rec[_span(sync_IDEN2.DLV_WN)])
            except ValueError:
                continue


def load_iden2(iden2_dir, rep):
    """Iden2, or None with a warning if a file is missing."""
    enlev_path = os.path.join(iden2_dir, 'enlev.dat')
    trans_path = os.path.join(iden2_dir, 'trans.dat')
    map_path = os.path.join(iden2_dir, 'IDEN_level_ids.txt')
    dlv_path = os.path.join(iden2_dir, 'dlv.dat')
    for p in (enlev_path, trans_path, map_path):
        if not os.path.exists(p):
            rep.warn('IDEN2', '%s is missing; the IDEN2 checks are skipped'
                     % p)
            return None
    dlv = None
    if os.path.exists(dlv_path):
        dlv, _ends = IDEN.read_records(dlv_path)
    else:
        rep.warn('IDEN2/dlv.dat', '%s is missing; the observed lines are not '
                                  'checked' % dlv_path)
    return Iden2(IDEN.Enlev(enlev_path), IDEN.Trans(trans_path),
                 IDEN.read_map(map_path), dlv)


def check_iden2(iden, cls, ctx, levels, rep, drift, gross, list_n):
    """The two IDEN2 files against each other, against the fit, and against
    the accepted assignments.

    IDEN2 is where the identifications are looked at, so a difference here is
    the one the eye will meet: a transition it shows as identified that the
    classification has since withdrawn looks, on the screen, exactly like one
    that is still accepted.
    """
    if iden is None:
        return
    _iden_internal(iden.enlev, iden.trans, rep, list_n)
    _iden_vs_fit(iden.enlev, iden.row_of, levels, rep, drift, gross, list_n)
    _iden_vs_classification(iden.enlev, iden.trans, iden.id_of_row, cls, ctx,
                            rep, list_n, iden.wn_of_line)


def _iden_internal(enlev, trans, rep, list_n):
    """trans.dat carries, for every transition, a copy of the partner's
    observed energy and found flag and the wavenumber that follows from them.
    IDEN2 rewrites all of it whenever a level energy changes, so a copy that
    no longer agrees with enlev.dat means trans.dat was written before the
    last change to enlev.dat and everything it predicts is stale."""
    bad_head, bad_part, bad_wn = [], [], []
    for owner, k in trans.header_of.items():
        rec = trans.records[k]
        if owner not in enlev.row_of:
            bad_head.append('level %d is in trans.dat and not in enlev.dat'
                            % owner)
            continue
        starred = '*' in IDEN.field(rec, IDEN.TH_STAR)
        if abs(IDEN.number(rec, IDEN.TH_EOBS) - enlev.e_obs(owner)) > 0.0005 \
                or starred != enlev.known(owner):
            bad_head.append('level %d: %.3f in trans.dat, %.3f in enlev.dat'
                            % (owner, IDEN.number(rec, IDEN.TH_EOBS),
                               enlev.e_obs(owner)))
    for (owner, partner), k in trans.row_of.items():
        rec = trans.records[k]
        if partner not in enlev.row_of or owner not in enlev.row_of:
            continue
        e_part = IDEN.number(rec, IDEN.TR_EPART)
        if abs(e_part - enlev.e_obs(partner)) > 0.0005:
            bad_part.append('%d - %d: partner at %.3f in trans.dat, %.3f in '
                            'enlev.dat' % (owner, partner, e_part,
                                           enlev.e_obs(partner)))
        pred = abs(enlev.e_obs(owner) - enlev.e_obs(partner))
        if abs(IDEN.number(rec, IDEN.TR_WN) - pred) > 0.0015:
            bad_wn.append('%d - %d: %.3f in trans.dat, %.3f from the two '
                          'level energies' % (owner, partner,
                                        IDEN.number(rec, IDEN.TR_WN), pred))
    items = (bad_head + bad_part + bad_wn)
    if items:
        rep.error('IDEN2',
                  'trans.dat is older than enlev.dat: %d level headers, %d '
                  'partner energies and %d predicted wavenumbers no longer '
                  'follow from the level energies'
                  % (len(bad_head), len(bad_part), len(bad_wn)), items)
        rep.act(50, 'Let IDEN2 rewrite trans.dat from the current enlev.dat; '
                    'every predicted wavenumber it shows is stale.')
    else:
        rep.ok('IDEN2', 'trans.dat follows from enlev.dat: %d level headers '
                        'and %d transitions all consistent'
               % (len(trans.header_of), len(trans.row_of)))


def _iden_vs_fit(enlev, row_of, levels, rep, drift, gross, list_n):
    """enlev.dat holds a copy of every level energy; the fit holds the current
    one.  The copy is updated by hand and is expected to lag."""
    if not levels:
        return
    joined, lag, big, absent = 0, [], [], []
    for lid, n in sorted(row_of.items()):
        if n not in enlev.row_of:
            absent.append('%s is mapped to enlev.dat row %d, which does not '
                          'exist' % (lid, n))
            continue
        if lid not in levels:
            absent.append('%s is in the IDEN2 map and not in %s'
                          % (lid, LOPT_LEVELS))
            continue
        joined += 1
        d = float(levels[lid]) - enlev.e_obs(n)
        if abs(d) >= gross:
            big.append('%s (row %d)  enlev %.3f, fit %.3f  (%+.3f)'
                       % (lid, n, enlev.e_obs(n), levels[lid], d))
        elif abs(d) > drift:
            lag.append('%s (row %d)  %+.3f' % (lid, n, d))
    if absent:
        rep.error('IDEN2', '%d levels cannot be joined between the IDEN2 map '
                           'and the fit' % len(absent), absent)
    if big:
        rep.error('IDEN2',
                  '%d levels sit at least %g cm^-1 from where the fit puts '
                  'them - enlev.dat and the fit are describing different '
                  'positions' % (len(big), gross), big)
        rep.act(50, 'Update the energies of those levels in IDEN2 (enlev.dat, '
                    'and let IDEN2 rewrite trans.dat) before looking at them '
                    'on the screen.')
    if lag:
        rep.warn('IDEN2', '%d levels are more than %g cm^-1 from the fit but '
                          'less than %g - the ordinary lag, worth an update '
                          'the next time IDEN2 is opened'
                 % (len(lag), drift, gross), lag)
    if not (absent or big or lag):
        rep.ok('IDEN2', 'all %d levels agree with the fit to within %g cm^-1'
               % (joined, drift))

    # levels the fit knows and IDEN2 does not, and the reverse
    unmapped_star = sorted(n for n in enlev.row_of
                           if enlev.known(n) and n not in row_of.values())
    unlisted = sorted(set(levels) - set(row_of))
    unstarred = sorted(lid for lid, n in row_of.items()
                       if n in enlev.row_of and not enlev.known(n))
    if unmapped_star:
        rep.warn('IDEN2', '%d levels are marked found in enlev.dat and '
                          'have no entry in IDEN_level_ids.txt'
                 % len(unmapped_star),
                 ['enlev.dat row %d at %.3f' % (n, enlev.e_obs(n))
                  for n in unmapped_star])
    if unlisted:
        rep.warn('IDEN2', '%d levels of the fit have no entry in '
                          'IDEN_level_ids.txt' % len(unlisted),
                 unlisted)
    if unstarred:
        rep.error('IDEN2', '%d levels of the fit are not marked found in '
                           'enlev.dat' % len(unstarred), unstarred)
        rep.act(50, 'Mark those levels found in enlev.dat; IDEN2 will not '
                    'offer their transitions until it is done.')


def _already(n):
    """The clause that stands in place of a list already printed above."""
    if not n:
        return ''
    return ('; %d of them are the transitions already listed under LOPT '
            'above and are not repeated here' % n)


def _iden_vs_classification(enlev, trans, id_of_row, cls, ctx, rep,
                            list_n, wn_of_line=None):
    """Which transitions IDEN2 shows as identified, against which the
    classification accepts.  This is the check the rest of the script exists
    to make."""
    if cls is None:
        return
    found, unmapped = read_iden_assignments(enlev, trans, id_of_row,
                                            wn_of_line)
    accepted = {}
    for pair, recs in classified_rows(classified(cls)).items():
        for wn, acc in recs:
            if acc:
                accepted[tuple(sorted(pair))] = wn

    only_iden = sorted(set(found) - set(accepted))
    only_cls = sorted(set(accepted) - set(found))
    moved, moved_note = [], {}
    for pair in sorted(set(found) & set(accepted)):
        if abs(found[pair] - accepted[pair]) > 0.01:
            moved.append((pair, accepted[pair]))
            moved_note[(pair, accepted[pair])] = (
                'IDEN2 has it on %.3f' % found[pair])
    if unmapped:
        rep.warn('IDEN2', '%d identified transitions of trans.dat join levels '
                          'that are not in IDEN_level_ids.txt' % len(unmapped),
                 ['rows %d - %d at %.3f' % u for u in unmapped])
    if only_iden:
        fresh = [p for p in only_iden if not ctx.already_listed(p, found[p])]
        rep.error('IDEN2',
                  '%d transitions are identified in trans.dat and are no '
                  'longer accepted - IDEN2 will show them as good%s'
                  % (len(only_iden), _already(len(only_iden) - len(fresh))),
                  ctx.items([(p, found[p]) for p in fresh]),
                  head=Context.HEAD)
        rep.act(60, 'Withdraw those %d assignments in IDEN2, or accept them '
                    'in the classification: at present the screen and the fit '
                    'disagree.' % len(only_iden))
    if only_cls:
        fresh = [p for p in only_cls
                 if not ctx.already_listed(p, accepted[p])]
        rep.error('IDEN2',
                  '%d accepted assignments are not identified in trans.dat - '
                  'IDEN2 will show those lines as unassigned%s'
                  % (len(only_cls), _already(len(only_cls) - len(fresh))),
                  ctx.items([(p, accepted[p]) for p in fresh]),
                  head=Context.HEAD)
        rep.act(60, 'Carry those %d new assignments into IDEN2 so that the '
                    'screen shows the fit that is actually being used.'
                    % len(only_cls))
    if moved:
        rep.error('IDEN2', '%d transitions are on a different observed '
                           'line in the two files' % len(moved),
                  ctx.items(moved, moved_note), head=Context.HEAD)
        rep.act(60, 'Re-identify those transitions in IDEN2: they have been '
                    'moved to another line since IDEN2 last saw them.')
    if not (only_iden or only_cls or moved):
        rep.ok('IDEN2', 'IDEN2 and the classification agree about all %d '
                        'identified transitions' % len(found))


# ---------------------------------------------------------------------------
# 10-11. IDEN2/dlv.dat
# ---------------------------------------------------------------------------
class DlvRow(object):
    """One row of dlv.dat: an observed line as IDEN2 holds it.

    ``index`` is the row's place in the file, counted from 0; ``row`` is the
    line number in columns 61-66, which is how trans.dat names the line.  The
    two agree in a file as it was first built and part company as soon as a
    line is inserted in IDEN2.
    """

    def __init__(self, index, record):
        self.index = index
        self.row = int(record[_span(sync_IDEN2.DLV_ROW)])
        self.wn = float(record[_span(sync_IDEN2.DLV_WN)])
        self.lam = float(record[_span(sync_IDEN2.DLV_LAMBDA)])
        self.u_lam = float(record[_span(sync_IDEN2.DLV_UNC)])

    @property
    def u_wn(self):
        """The uncertainty on the wavenumber scale.

        The file states it as a wavelength uncertainty in angstroms, so it is
        converted by u_wn = u_lambda * wn / lambda.  The row's own wavelength
        is used and not 1e8 / wn, because the wavelength there is the standard
        one - vacuum in the ultraviolet, air above 2000 A - and the two differ
        by up to three angstroms.
        """
        return self.u_lam * self.wn / self.lam


def _span(field):
    return slice(field[0], field[1])


def check_dlv(iden, cfg_path, rep, dlv_tol, list_n):
    """The observed line list IDEN2 shows, against the set's own line list.

    This is the file the eye reads a measurement off: the wavenumber it is
    looking for a transition near, and the uncertainty that says how far from
    the prediction a line may be and still be the line.  IDEN2 rewrites it
    when a line is edited on its screen, rebuilding every air row's
    wavenumber from its wavelength, and no earlier check touches it, so a set
    seeded by copying another set's IDEN2 directory keeps the other set's
    wavenumbers - which is silent, because every screen goes on working and
    the numbers on it are of the right size.
    """
    if iden is None or iden.dlv is None:
        return
    rows = _dlv_internal(iden.dlv, rep, list_n)
    if not rows:
        return
    _dlv_vs_trans(rows, iden.trans, rep, list_n)
    _dlv_dispersion(iden.dlv, rep)
    if not os.path.exists(cfg_path):
        rep.warn('IDEN2/dlv.dat', '%s does not exist; dlv.dat is not checked '
                                  'against the line list' % cfg_path)
        return
    cfg = config.load(cfg_path)
    if not os.path.exists(cfg.lines_file):
        rep.warn('IDEN2/dlv.dat', '%s does not exist; dlv.dat is not checked '
                                  'against the line list' % cfg.lines_file)
        return
    _dlv_vs_lines(rows, cfg, rep, dlv_tol, list_n)


def _dlv_internal(records, rep, list_n):
    """The file's own arithmetic: fixed width, one line number per row, and
    one wavenumber per line.

    The line number, in columns 61-66, is how every assignment in trans.dat
    names its line, so two rows carrying one line number, or two rows
    carrying one wavenumber, make an assignment mean something other than
    what it says.  The line numbers need not follow the order of the rows.
    The file is first built with each row numbered by its place, but a line
    inserted in IDEN2 later - one the line list has and dlv.dat was built
    without - takes the next free number and is placed among the others by
    its wavenumber, and every row keeps the number it had.
    """
    rows, malformed = [], []
    for i, rec in enumerate(records):
        if not rec.strip():
            continue
        if len(rec) < sync_IDEN2.DLV_WIDTH:
            malformed.append('record %d is %d characters, not %d'
                             % (i + 1, len(rec), sync_IDEN2.DLV_WIDTH))
            continue
        try:
            rows.append(DlvRow(i, rec))
        except ValueError:
            malformed.append('record %d does not parse: %r'
                             % (i + 1, rec[:40]))
    if malformed:
        rep.error('IDEN2/dlv.dat',
                  '%d records of dlv.dat are not of the fixed 66-character '
                  'form the file is read by' % len(malformed), malformed)
        return rows
    records_of = {}
    for r in rows:
        records_of.setdefault(r.row, []).append(r)
    shared = ['line %d is on records %s' % (n, ', '.join(
                  '%d (%.3f)' % (r.index + 1, r.wn) for r in rs))
              for n, rs in sorted(records_of.items()) if len(rs) > 1]
    bad = ['record %d carries line number %d' % (r.index + 1, r.row)
           for r in rows if r.row < 1]
    if shared or bad:
        rep.error('IDEN2/dlv.dat',
                  '%d line numbers of dlv.dat are not one line each, and '
                  'every assignment in trans.dat names its line by that '
                  'number' % (len(shared) + len(bad)), shared + bad)
        rep.act(10, 'Restore dlv.dat: two rows carry one line number, so no '
                    'assignment in trans.dat that names it can be trusted to '
                    'mean the line it means.')
        return rows
    # A line inserted in IDEN2 carries a larger number than the row after it.
    inserted = [a for a, b in zip(rows, rows[1:]) if b.row < a.row]
    seen, repeated = {}, []
    for r in rows:
        key = round(r.wn, 3)
        if key in seen:
            repeated.append('lines %d and %d both carry %.3f'
                            % (seen[key], r.row, r.wn))
        seen[key] = r.row
    if repeated:
        rep.error('IDEN2/dlv.dat',
                  '%d wavenumbers appear on more than one row of dlv.dat; the '
                  'pipeline names an observed line by its wavenumber and '
                  'cannot tell which row is meant' % len(repeated), repeated)
    out_of_order = sum(1 for a, b in zip(rows, rows[1:]) if b.wn > a.wn)
    if out_of_order:
        rep.warn('IDEN2/dlv.dat',
                 'dlv.dat is not in decreasing order of wavenumber: %d rows '
                 'stand above a larger wavenumber. IDEN2 lists the lines in '
                 'the order the file is in' % out_of_order)
    if not (repeated or out_of_order):
        rep.ok('IDEN2/dlv.dat', 'dlv.dat holds %d rows, one line number and '
                                'one wavenumber each, in decreasing order%s'
               % (len(rows), '; %d of them inserted in IDEN2 under a later '
                             'line number' % len(inserted) if inserted
                  else ''),
               ['line %d at record %d (%.3f)' % (r.row, r.index + 1, r.wn)
                for r in inserted])
    return rows


def _dlv_dispersion(records, rep):
    """Every row's wavelength against its wavenumber.

    IDEN2 treats the wavelength as the measurement: when it saves dlv.dat it
    rebuilds the wavenumber of every air row from the row's wavelength.  A
    row whose two numbers disagree - which is what the rescaling of the old
    sync_IDEN2.py left behind - therefore has its wavenumber moved the next
    time any line is edited on IDEN2's screen, with nothing on the screen to
    say so.
    """
    bad, n = sync_IDEN2.dispersion_disagreements(records)
    if not n:
        return
    if bad:
        rep.warn('IDEN2/dlv.dat',
                 '%d rows of dlv.dat carry a wavelength that is not the '
                 'standard one for their wavenumber; IDEN2 will move their '
                 'wavenumbers to match the wavelengths the next time it saves '
                 'the file' % len(bad),
                 ['line %6d  %12.3f  %11.4f A, should be %11.4f A'
                  % b for b in bad],
                 head='line          wavenumber   wavelength')
        rep.act(66, 'Run sync_IDEN2.py in this set: it writes every '
                    "wavelength from the line list's wavenumber.")
    else:
        rep.ok('IDEN2/dlv.dat', 'the wavelength of each of the %d rows of '
                                'dlv.dat is the standard one for its '
                                'wavenumber' % n)


def _dlv_allowance(r, dlv_tol):
    """How far row ``r`` of dlv.dat may stand from its line's wavenumber and
    still be that wavenumber printed.

    ``dlv_tol``, or for an air row the rounding of both printed fields, if
    that is more: the wavenumber to 0.0005 cm^-1, and the wavelength IDEN2
    rebuilds it from to 0.00005 A, which is ``0.00005 * wn / lambda`` cm^-1
    and reaches 0.0013 cm^-1 near 2000 A.
    """
    if r.wn > sync_IDEN2.VACUUM_ABOVE:
        return dlv_tol
    return max(dlv_tol, 0.0005 + 0.00005 * r.wn / r.lam + 1e-6)


def _dlv_vs_lines(rows, cfg, rep, dlv_tol, list_n):
    """Every wavenumber and every uncertainty in dlv.dat against the line list
    the set is classified and fitted on.

    The set's own line list, read on the columns its configuration names, so
    that the baseline set is checked against Sugar's wavenumbers and a
    corrected set against the corrected ones.  The join is on the wavenumber
    the row should be carrying; a row that instead matches the set's
    ``wn_key`` column - the wavenumber that never moves, which is the scale
    dlv.dat was built on - is a row the correction has never been applied to,
    and that is reported as what it is rather than as an unknown line.

    Both joins accept a row within ``sync_IDEN2.DLV_MATCH`` of its line, the
    tolerance sync_IDEN2.py matches on.  A row further from its line than
    printing accounts for (``_dlv_allowance``) has drifted: it is still that
    line, and sync_IDEN2.py puts it back, so it is reported as drifted and
    not as a line of unknown origin.

    The uncertainty a row should carry is the one the pipeline classifies and
    fits the line on: the line list's own, raised to the value
    inflated_unc_lines.txt sets for the line where that is larger - the rule
    classify_lines.py applies.  An uncertainty widened in IDEN2 is therefore
    in step once it is entered in the registry, and one widened in IDEN2 and
    never entered is reported on its own, since the cure for it is the
    registry and not sync_IDEN2.py, which writes the pipeline's value.
    """
    try:
        lines = sync_IDEN2.read_line_list(cfg, lambda *a: None)
    except sync_IDEN2.SyncError as exc:
        rep.error('IDEN2/dlv.dat', 'the line list cannot be read: %s' % exc)
        return
    key_name = cfg.lines.columns.get('wn_key') or cfg.lines.columns['wn']
    wn_name = cfg.lines.columns['wn']
    registry = _dlv_registry(cfg, [key for key, _wn, _u in lines], rep)
    want = sorted((wn, u, key) for key, wn, u in lines)
    want_wn = [t[0] for t in want]
    keyed = sorted((key, wn) for key, wn, _u in lines)
    keyed_wn = [t[0] for t in keyed]

    stale, unknown, narrow, wide, moved = [], [], [], [], []
    drifted = 0
    worst = worst_moved = 0.0
    n_registry = 0
    matched = set()
    join = max(dlv_tol, sync_IDEN2.DLV_MATCH)
    for r in rows:
        allow = _dlv_allowance(r, dlv_tol)
        k = _nearest(want_wn, r.wn, join)
        if k is not None and abs(want_wn[k] - r.wn) > allow:
            j = _nearest(keyed_wn, r.wn, allow)
            if j is not None and key_name != wn_name:
                k = None            # it is the key of a line, exactly
        if k is None:
            j = _nearest(keyed_wn, r.wn, join)
            if j is not None and key_name != wn_name:
                stale.append('line %6d  %12.3f -> %12.3f  (%+.3f)'
                             % (r.row, r.wn, keyed[j][1],
                                keyed[j][1] - r.wn))
            else:
                unknown.append('line %6d  %12.3f' % (r.row, r.wn))
            continue
        matched.add(k)
        if abs(want_wn[k] - r.wn) > allow:
            worst_moved = max(worst_moved, abs(want_wn[k] - r.wn))
            moved.append('line %6d  %12.3f  %12.4f  (%+.4f)'
                         % (r.row, r.wn, want_wn[k], r.wn - want_wn[k]))
        elif abs(want_wn[k] - r.wn) > 0.0005:
            drifted += 1
            worst = max(worst, abs(want_wn[k] - r.wn))
        if registry is None:
            continue
        u_list, key = want[k][1], want[k][2]
        u_set = registry.lookup(key)
        u_want = u_list if u_set is None else max(u_list, u_set)
        if u_set is not None and u_set > u_list:
            n_registry += 1
        quantum = 0.00005 * r.wn / r.lam
        if abs(r.u_wn - u_want) <= quantum + DLV_U_REL * u_want:
            continue
        item = ('line %6d  %12.3f  dlv %9.4f  pipeline %9.4f  cm-1  (list '
                '%.4f%s)' % (r.row, r.wn, r.u_wn, u_want, u_list,
                             '' if u_set is None
                             else ', %s %.4f' % (INFLATED, u_set)))
        (wide if r.u_wn > u_want else narrow).append(item)
    absent = ['%12.3f  u %9.4f  cm-1' % want[k][:2] for k in range(len(want))
              if k not in matched]

    head = 'line          wavenumber'
    if stale:
        rep.error('IDEN2/dlv.dat',
                  '%d rows of dlv.dat are still on the uncorrected %s scale: '
                  'IDEN2 is showing wavenumbers the set no longer uses'
                  % (len(stale), key_name),
                  stale, head='line     dlv.dat now      %s' % wn_name)
        rep.act(65, 'Run sync_IDEN2.py in this set: dlv.dat holds the '
                    'wavenumbers of another set, so every window drawn around '
                    'a prediction is drawn in the wrong place.')
    if narrow:
        rep.error('IDEN2/dlv.dat',
                  '%d uncertainties in dlv.dat are narrower than the ones the '
                  'pipeline classifies and fits on, by more than their '
                  'printing precision; IDEN2 judges every identification by '
                  'this number' % len(narrow), narrow)
        rep.act(65, 'Run sync_IDEN2.py in this set: it writes the '
                    'pipeline\'s uncertainty - the line list\'s, or the '
                    'value in %s where that is larger - into dlv.dat.'
                    % INFLATED)
    if wide:
        rep.error('IDEN2/dlv.dat',
                  '%d uncertainties in dlv.dat are wider than the ones the '
                  'pipeline classifies and fits on, by more than their '
                  'printing precision: widened in IDEN2 and not entered in '
                  '%s, or left behind when the line list narrowed them'
                  % (len(wide), INFLATED), wide)
        rep.act(65, 'Enter in %s each wider uncertainty in dlv.dat that was '
                    'set on purpose, so that classify_lines.py and LOPT use '
                    'it too.  Do not run sync_IDEN2.py on the others until '
                    'that is done: it would narrow every one of them to the '
                    'pipeline\'s value.' % INFLATED)
    if moved:
        rep.warn('IDEN2/dlv.dat',
                 "%d rows of dlv.dat have drifted from their line's %s by "
                 'more than printing accounts for (up to %.4f cm^-1): each is '
                 'still its line, but IDEN2 shows it slightly off'
                 % (len(moved), wn_name, worst_moved),
                 moved, head='line     dlv.dat now      %s' % wn_name)
        rep.act(66, 'Run sync_IDEN2.py in this set: it writes the line '
                    "list's wavenumber, and a wavelength consistent with it, "
                    'back into every drifted row of dlv.dat.')
    if unknown:
        rep.warn('IDEN2/dlv.dat',
                 '%d rows of dlv.dat match no line of %s within %.4f cm^-1'
                 % (len(unknown), os.path.basename(cfg.lines_file), join),
                 unknown, head=head)
    if absent:
        rep.warn('IDEN2/dlv.dat',
                 '%d lines of %s have no row in dlv.dat, so they cannot be '
                 'identified on the screen until they are inserted in IDEN2, '
                 'where each takes the next free line number'
                 % (len(absent), os.path.basename(cfg.lines_file)),
                 absent, head='  wavenumber  uncertainty')
    if not (stale or narrow or wide or moved):
        if registry is None:
            said = 'the %s' % wn_name
        else:
            said = 'the %s and the uncertainty the pipeline uses' % wn_name
        rep.ok('IDEN2/dlv.dat',
               'all %d matched rows of dlv.dat carry %s (%d of them differ '
               'in the last printed digit, by up to %.4f cm^-1)%s'
               % (len(matched), said, drifted, worst,
                  '; %d of the uncertainties are the wider values of %s'
                  % (n_registry, INFLATED) if n_registry else ''))


def _dlv_registry(cfg, wn_keys, rep):
    """The registry of hand-set uncertainties the set's configuration names,
    checked against the line list's ``wn_keys``.

    An empty registry if the configuration names none or the file is absent.
    None if it cannot be used - an entry too short to tell two lines apart,
    or one line entered twice with two values - in which case the pipeline
    refuses it as well and the uncertainties are not compared.
    """
    path = getattr(cfg, 'inflated_unc', '') or ''
    if not path or not os.path.exists(path):
        return hfs_kappa.Registry()
    try:
        registry = hfs_kappa.read_inflated(path)
        unused = registry.check(wn_keys, os.path.basename(path))
    except ValueError as exc:
        rep.error('IDEN2/dlv.dat', '%s cannot be applied, so the '
                                   'uncertainties of dlv.dat are not compared:'
                                   ' %s' % (INFLATED, exc))
        return None
    if unused:
        rep.warn('IDEN2/dlv.dat',
                 '%d entries of %s name no line of %s'
                 % (len(unused), INFLATED, os.path.basename(cfg.lines_file)),
                 unused)
    return registry


def _nearest(values, wn, tol):
    """The index of the value nearest ``wn``, or None if none is within
    ``tol``.  ``values`` is sorted."""
    k = bisect.bisect_left(values, wn)
    best = None
    for m in (k - 1, k, k + 1):
        if 0 <= m < len(values) and abs(values[m] - wn) <= tol:
            if best is None or abs(values[m] - wn) < abs(values[best] - wn):
                best = m
    return best


def _dlv_vs_trans(rows, trans, rep, list_n):
    """Every identification in trans.dat against the row of dlv.dat it names.

    An identification names its observed line by the line number in columns
    61-66 of dlv.dat, and IDEN2 shows the row that carries it.  The
    wavenumber an identification also carries is ignored, by IDEN2 and here:
    the line number is the identification.
    """
    have = set(r.row for r in rows)
    missing = []
    n = 0
    for (owner, partner), k in trans.row_of.items():
        obs = IDEN.assignment(trans.records[k])
        if not IDEN.has_line(obs):
            continue
        n += 1
        row = IDEN.obs_row(obs)
        if row not in have:
            missing.append('levels %d - %d at %.3f name line %d'
                           % (partner, owner, IDEN.obs_wavenumber(obs), row))
    if missing:
        rep.error('IDEN2/dlv.dat',
                  '%d identifications in trans.dat name a line number that '
                  'dlv.dat has not got' % len(missing), missing)
    else:
        rep.ok('IDEN2/dlv.dat', 'all %d identifications in trans.dat name a '
                                'line that dlv.dat has' % n)


# ---------------------------------------------------------------------------
# 12. Stability
# ---------------------------------------------------------------------------
def check_stability(paths, cls, ctx, rep, list_n):
    """Two ways the classification can be unsettled rather than merely stale.

    An oscillating assignment is one classify_lines.py accepted on one pass
    and rejected on the next until it blacklisted it; the run converged, but
    only because that assignment was taken out of the argument, and which way
    it should have gone is undecided.  An unobeyed ledger row is worse: the
    analyst has ruled on a line and the classification does not show the
    ruling.
    """
    path = paths.get(UNSTABLE)
    if path and os.path.exists(path):
        un = pd.read_csv(path, dtype={'low_id': str, 'upp_id': str})
        if len(un):
            rows, notes = [], {}
            for _, r in un.iterrows():
                key = (pair_key(r['low_id'], r['upp_id']),
                       float(r['wn_obs']))
                rows.append(key)
                notes[key] = '%s -> %s' % (r.get('oscillated_between', ''),
                                           r.get('states', ''))
            rep.warn('Stability',
                     '%d assignments were withdrawn for oscillating: the run '
                     'converged only because they were blacklisted'
                     % len(un), ctx.items(rows, notes), head=Context.HEAD)
            rep.act(10, 'Decide the %d oscillating assignment(s) in %s by '
                        'hand and record the verdict in %s, so that the next '
                        'run is stable rather than merely converged.'
                        % (len(un), UNSTABLE, DECISIONS))
        else:
            rep.ok('Stability', 'no assignment oscillated in the last run')
    else:
        rep.warn('Stability', '%s is missing; classify_lines.py has not been '
                              'run, or was run with --no-write' % UNSTABLE)

    path = paths.get(DECISIONS)
    if not (path and os.path.exists(path)) or cls is None:
        return
    dec = pd.read_csv(path, dtype={'low_id': str, 'upp_id': str})
    if 'wn_key' not in dec.columns:
        rep.error('Stability', '%s has no wn_key column' % DECISIONS)
        return
    # The ledger names a line by Sugar's wavenumber (wn_key); a calibrated
    # set's wn_obs is the corrected one, so the join is on its wn_key.
    state = classified_rows(classified(cls), 'wn_key' if 'wn_key'
                            in cls.columns else 'wn_obs')
    # The lines on which the ledger orders each pair accepted.  A rejection
    # with an accept of the same pair on another line is how an assignment is
    # moved; the transition settling there is the ruling obeyed, not a
    # transition the ledger says nothing about.
    ruled = {}
    for low, upp, wn, d in zip(dec['low_id'], dec['upp_id'], dec['wn_key'],
                               dec['decision']):
        if str(d).strip().lower() == 'accept':
            ruled.setdefault(pair_key(low, upp), []).append(float(wn))
    unapplied, unknown, elsewhere = [], [], []
    notes = {}
    for _, r in dec.iterrows():
        key = pair_key(r['low_id'], r['upp_id'])
        want = str(r['decision']).strip().lower() == 'accept'
        wn = float(r['wn_key'])
        # A ledger row rules on ONE observed line, and classify_lines.py finds
        # it within DECISIONS_WN_MATCH.  The same pair on another line is a
        # different question, and is how an assignment is moved from one line
        # to another: both rows are then obeyed and neither contradicts the
        # other.
        here = [(w, a) for w, a in state.get(key, [])
                if abs(w - wn) <= DECISIONS_WN_MATCH]
        if not here:
            if want:
                unknown.append((key, wn))
                notes[(key, wn)] = 'the ledger orders it accepted'
            else:
                unruled = [x for x, a in state.get(key, []) if a
                           and not any(abs(x - y) <= DECISIONS_WN_MATCH
                                       for y in ruled.get(key, []))]
                if unruled:
                    w = unruled[0]
                    elsewhere.append((key, w))
                    notes[(key, w)] = ('rejected on %.3f as ordered, now '
                                       'accepted on this line' % wn)
            continue
        got_wn, got = here[0]
        if got != want:
            unapplied.append((key, wn))
            notes[(key, wn)] = ('the ledger says %s, the classification says '
                                '%s' % ('accept' if want else 'reject',
                                        'accept' if got else 'reject'))
    if unknown:
        rep.error('Stability',
                  '%d ledger rows order an assignment the classification does '
                  'not hold at all' % len(unknown),
                  ctx.items(unknown, notes), head=Context.HEAD)
    if unapplied:
        rep.error('Stability',
                  '%d ledger verdicts are not obeyed by the classification - '
                  'it was written before those rows, or by a run that did not '
                  'read the ledger' % len(unapplied),
                  ctx.items(unapplied, notes), head=Context.HEAD)
        rep.act(10, 'Re-run classify_lines.py: %s has verdicts the '
                    'classification does not obey.' % DECISIONS)
    if elsewhere:
        rep.warn('Stability',
                 '%d ledger rejections were obeyed on the line they name, but '
                 'the transition has settled on another line - the ledger '
                 'says nothing about that one' % len(elsewhere),
                 ctx.items(elsewhere, notes), head=Context.HEAD)
    if not (unknown or unapplied or elsewhere):
        rep.ok('Stability', 'all %d ledger verdicts are obeyed' % len(dec))


# ---------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description='Check that the classification, the LOPT files and the '
                    'IDEN2 files describe the same identification, and that '
                    'the last classification run was stable.')
    ap.add_argument('--drift', type=float, default=DEF_DRIFT, metavar='CM',
                    help='an IDEN2 level energy this far from the fit is the '
                         'expected lag, not a mismatch (default: %(default)s)')
    ap.add_argument('--gross', type=float, default=DEF_GROSS, metavar='CM',
                    help='at this distance the two files are describing '
                         'different positions (default: %(default)s)')
    ap.add_argument('--wn-tol', type=float, default=DEF_WN_TOL, metavar='CM',
                    dest='wn_tol',
                    help='how far two files may put one transition apart in '
                         'wavenumber and still mean the same observed line '
                         '(default: %(default)s)')
    ap.add_argument('--dlv-tol', type=float, default=DEF_DLV_TOL,
                    metavar='CM', dest='dlv_tol',
                    help='how far a wavenumber in IDEN2/dlv.dat may stand '
                         'from the line list and still be the same number '
                         'printed twice (default: %(default)s)')
    ap.add_argument('--list', type=int, default=DEF_LIST, metavar='N',
                    dest='list_n',
                    help='how many items to print under each finding '
                         '(default: %(default)s)')
    ap.add_argument('--quiet', action='store_true',
                    help='print only the checks that found something')
    ap.add_argument('--out', metavar='FILE', default=REPORT,
                    help='where the complete report is written, every item '
                         'of every list (default: %(default)s)')
    ap.add_argument('--csv', metavar='FILE',
                    help='also write the findings as a table')
    ap.add_argument('--iden2', metavar='DIR', default=None,
                    help='the IDEN2 directory (default: IDEN2, resolved like '
                         'every other file)')
    ap.add_argument('--set', metavar='DIR', default=None, dest='set_dir',
                    help='the working set to check: the directory holding its '
                         'own classification table, LOPT files and IDEN2 '
                         '(iter, final). A file that directory does not hold '
                         'is taken from the project directory, which is what '
                         'makes the hand-kept ledgers shared between the sets '
                         '(default: the current directory)')
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    # `--set iter` is `working_path` asked to look in `iter` first: the files
    # that set has of its own - its classification table, its LOPT files, its
    # IDEN2 - are its own, and the ones it has not got, the decision ledger
    # and the level overrides, come from the project directory, which is
    # exactly the arrangement the sets are meant to have.  The report is
    # written into the set, since it describes that set.
    set_dir = os.path.abspath(args.set_dir or os.getcwd())
    if not os.path.isdir(set_dir):
        raise SystemExit('check_sync.py: --set %s is not a directory'
                         % args.set_dir)
    out = args.out if os.path.isabs(args.out) else os.path.join(set_dir,
                                                                args.out)
    csv_out = args.csv
    if csv_out and not os.path.isabs(csv_out):
        csv_out = os.path.join(set_dir, csv_out)
    args.out, args.csv = out, csv_out
    output_files.require_writable([args.out, args.csv], 'report file')
    names = [CLASSIFICATIONS_CSV, CLASSIFICATIONS_XLSX, LOPT_IN, LOPT_LINES,
             LOPT_LEVELS, DECISIONS, REVISIONS, UNSTABLE, POSITIONS,
             LINES_XLSX, NAME_CONFIG]
    paths = {n: working_path(n, cwd=set_dir) for n in names}
    iden2 = args.iden2 or working_path(NAME_IDEN2, cwd=set_dir)

    head = ['check_sync.py - are the files describing the same '
            'identification?',
            'working set: %s' % set_dir,
            '%-32s %s' % ('file', 'last written')]
    for n in names + ['IDEN2/enlev.dat', 'IDEN2/trans.dat',
                      'IDEN2/dlv.dat']:
        p = paths.get(n) or os.path.join(iden2, os.path.basename(n))
        head.append('%-32s %s'
                    % (n, _stamp(p) if os.path.exists(p) else 'not found'))
    print('\n'.join(head))

    rep = Report(args.list_n)
    check_freshness(paths, rep)
    cls, _ = check_csv_xlsx(paths[CLASSIFICATIONS_CSV],
                            paths[CLASSIFICATIONS_XLSX], rep, args.list_n)

    levels = (read_lopt_levels(paths[LOPT_LEVELS])
              if os.path.exists(paths[LOPT_LEVELS]) else None)
    lopt_in = (read_lopt_input(paths[LOPT_IN])
               if os.path.exists(paths[LOPT_IN]) else None)
    iden = load_iden2(iden2, rep)
    ctx = Context(cls, levels, lopt_in, iden, args.wn_tol)

    check_lopt_chain(cls, paths, ctx, rep, args.wn_tol)
    check_revisions(paths, levels, rep, args.gross)
    check_iden2(iden, cls, ctx, levels, rep, args.drift, args.gross,
                args.list_n)
    check_dlv(iden, paths[NAME_CONFIG], rep, args.dlv_tol, args.list_n)
    check_stability(paths, cls, ctx, rep, args.list_n)

    if rep.worst != OK:
        rep.act(70, 'Re-run level_positions.py (and level_shifts.py) last: '
                    'their reports are built on everything above, and a '
                    'report made before the fixes will recommend things that '
                    'have already been done.')

    n_err = sum(1 for f in rep.findings if f['severity'] == ERROR)
    n_warn = sum(1 for f in rep.findings if f['severity'] == WARN)
    summary = ('%d error(s), %d warning(s), %d check(s) clean.'
               % (n_err, n_warn,
                  sum(1 for f in rep.findings if f['severity'] == OK)))

    rep.print(quiet=args.quiet)
    rep.print_actions()
    print('\n' + summary)

    if args.out:
        with open(args.out, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write('\n'.join(head) + '\n')
            rep.print(quiet=False, out=fh, full=True)
            rep.print_actions(out=fh)
            fh.write('\n' + summary + '\n')
        print('the complete report, with every item of every list, is in %s'
              % args.out)

    if args.csv:
        rep.to_frame().to_csv(args.csv, index=False, lineterminator='\n')
        print('findings written to %s' % args.csv)

    return 2 if n_err else (1 if n_warn else 0)


def _stamp(path):
    return (pd.Timestamp(os.path.getmtime(path), unit='s', tz='UTC')
            .tz_convert(None).strftime('%Y-%m-%d %H:%M'))


if __name__ == '__main__':
    sys.exit(main())
