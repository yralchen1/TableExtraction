"""Work through the accept/reject disagreements between two sets of LOPT files.

A "set" here is one of the three directories the README describes: the
baseline `LineClass/`, the calibrated working set `LineClass/iter/`, and the
publication set `LineClass/final/`.  Each set has its own line list, its own
classification table and its own LOPT files, and the three share one verdict
ledger, `line_decisions.csv`, keyed on the wavenumber a line is *named* by
(`wn_key` - Sugar's own value, which no calibration moves).

When a set is re-measured, the classification changes its mind about some
assignments.  A line that the baseline accepted as one component of a blend
can fail in the calibrated set, because the corrected wavenumbers carry much
tighter uncertainties: an assignment that stood 1.1 sigma from its Ritz value
on Sugar's uncertainty stands 5 sigma from it on the fitted one, and the
classification rejects it.  Each of those disagreements has to be looked at
by hand, because the question - is this observed feature one line or two? -
is answered on the IDEN2 screen and not by any number in these files.

This tool does the mechanical half of that review twice over.

**Without ``--apply``** it writes a worksheet, one row per disagreement, with
every quantity the judgment needs already worked out: the observed line and
its uncertainty, every transition the classification proposes for it with its
predicted intensity and its departure from its Ritz value, the center of
gravity of those departures with and without the disputed component, how much
that center of gravity improves or worsens, which transition holds the line
in IDEN2 today and how much brighter it is predicted to be, and the two
levels with their IDEN2 row numbers.  The rows are ordered by the IDEN2 row
of the upper level, which is the order the levels are reviewed in.

Both directions of disagreement are reported, and the ``direction`` column
says which one a row is.  The common one is an assignment the baseline fits
and the set does not.  The other happens too: a corrected wavenumber can
bring a line *into* a Ritz window it missed before, and then the set fits an
assignment the baseline never did, which is a decision nobody has taken
either.  A transition flagged ``P`` in a set's LOPT input and one absent from
it altogether are the same answer to the only question asked here - does this
set fit the assignment? - so a pair that neither set fits is not a
disagreement, and a pair both fit is not one.

A disagreement in which the set's classification never proposed the
transition at all - the corrected wavenumber put the line outside its Ritz
window, so the table says nothing whatever about it - is the hardest of all
to judge by hand, and is assembled here rather than left blank.  The line
itself, its corrected wavenumber and uncertainty, its intensity and its
character, comes from another component of the same blend in the set's own
table, and failing that from the set's line list, matched on the wavenumber
the line is named by; the predicted intensity comes from the baseline's
table, which did propose the transition; and the departure from the Ritz
value is worked out here against the set's own level energies.  The ``why``
column says which of those sources each row used, and the row then joins the
blend arithmetic like any other.

**With ``--apply``** it reads the worksheet back and carries out whatever has
been written in its ``decision`` column, in all five places a decision has to
land:

* ``line_decisions.csv`` - an ``accept`` or ``reject`` row, so that the next
  classification run keeps the verdict;
* ``inflated_unc_lines.txt`` - the widened uncertainty, when one is given,
  under the line's wn_key to four decimals;
* ``LOPT_input_lines.txt`` of the set - the flag, the weight and the
  uncertainty of every record of the observed line.  All three are properties
  of the line as a whole, so all of its records are rewritten together: the
  flag is empty for an accepted classification and ``P`` for one that is not;
  the accepted classifications divide the line's weight in proportion to their
  calculated intensities, so accepting a second component changes the weight
  of the first; and the records share one uncertainty.  The rules are
  ``make_LOPT_input``'s own - ``blend_weights``, ``total_unc`` and
  ``blend_uncertainty`` - so the file says what a full regeneration would say;
* ``IDEN2/dlv.dat`` - the same uncertainty converted to the wavelength scale
  that file states it on, ``u_lambda = u_wn * lambda / wn``, with the row's
  own wavelength, so that the window IDEN2 judges the identification against
  on the screen is the window the ledger means;
* ``IDEN2/trans.dat`` - the assignment itself, for a transition that becomes
  accepted.

Nothing is written without a decision written by hand.  The ``suggested``
column is advisory, it is never copied into ``decision``, and ``--apply``
ignores a row whose ``decision`` is empty.

Why the suggestion is only ever advisory
----------------------------------------
The recipe the suggestions follow was checked against the half of the list
already reviewed by hand - 39 disagreements at upper levels up to IDEN2 row
830, of which 3 were accepted, 3 rejected on the ledger and 33 left rejected.
Three rules reproduce a verdict reliably:

* the disputed transition is predicted at least 20 IDEN2 intensity units
  fainter than the transition holding the line (6 of 6 correct) - the units
  are ``round(10*ln(Icalc))``, so 20 units is a factor of e^2;
* there is nothing at the line to blend it with (5 of 6);
* the center of gravity of the blend would stand more than 0.75 cm^-1 from
  the Ritz value (6 of 7).

Together they account for 22 of the 35 rejections - and they would also have
rejected 2 of the 3 assignments that were in fact accepted.  The remaining 13
rejections are not distinguishable, by any quantity in these files, from the
3 acceptances: the center of gravity improves, the intensities are
comparable, and the verdict went the other way.  The judgment is being made
on the IDEN2 screen.  So the tool lays out the evidence and leaves the
verdict to the analyst; it does not guess.

Usage
-----
::

    python review_mismatches.py --set iter
    python review_mismatches.py --set iter --from-iden2-id 831
    python review_mismatches.py --set iter --apply
    python review_mismatches.py --set iter --apply --dry-run

``--from-iden2-id`` starts at an upper level, for picking the review up where
it was left.  ``--all`` keeps the disagreements that already carry a ledger
row, which are otherwise left out as decided.
"""
import argparse
import bisect
import collections
import csv
import datetime
import io
import math
import os
import sys

import config
import hfs_correction
import hfs_kappa
import make_LOPT_input
import output_files
import swap_line_assignments_IDEN as IDEN
import sync_IDEN2 as sync
from swap_paths import working_path

HERE = os.path.dirname(os.path.abspath(__file__))

NAME_CLASSIFICATIONS = 'line_classifications.csv'
NAME_LOPT_IN = 'LOPT_input_lines.txt'
NAME_DECISIONS = 'line_decisions.csv'
NAME_INFLATED = 'inflated_unc_lines.txt'
NAME_IDEN2 = 'IDEN2'
NAME_WORKSHEET = 'decisions.txt'
NAME_LOPT_LEVELS = 'LOPT_output_levels.txt'

# The set's line list, read only to fill in a line the set's classification
# never proposed the transition for: the corrected wavenumber and its
# uncertainty live in it whether any Ritz window caught the line or not.
NAME_CORRECTED = 'Pr3_lines_corrected.xlsx'
COL_OWN = 'own'
COL_OWN_CORR = 'own_corr'
COL_UNC_CORR = 'unc_own_corr'
COL_INTENS = 'Icor'
COL_CHAR = 'Ch.'

# The columns of LOPT's transitions file, from make_LOPT_input.FIELDS.
LOPT_FIELDS = {'wavenumber': (1, 12), 'uncertainty': (14, 19),
               'lower_level': (43, 55), 'upper_level': (59, 71),
               'flags': (73, 77)}

# A wavenumber in the ledger names a line to within this much; the value
# classify_lines.attach_line_decisions matches on.
WN_MATCH = 0.01

# How far the center of gravity of a blend may stand from its Ritz value and
# the blend still be worth proposing.  Beyond this the assignment stays
# rejected whatever else is true of it.
CG_LIMIT = 0.75

# The inflated uncertainty preferred when the center of gravity is within a
# factor PREFER_FACTOR of it: the far-ultraviolet blends reconcile with Ritz
# on one constant value rather than on a value each.
PREFER_UNC = 0.5
PREFER_FACTOR = 1.5

# An IDEN2 intensity code is round(10*ln(Icalc)), so a difference of 20 is a
# factor of e^2 in predicted intensity.
CODE_MUCH_BRIGHTER = 20

WORKSHEET_COLUMNS = [
    # what was asked for, in the order it was asked for
    'obs_wn', 'unc_obs_wn', 'level_id1', 'level_id2', 'iden2_id1', 'iden2_id2',
    'lc_flag', 'iter_flag', 'unc_wn_inflated', 'decision', 'Notes',
    # the evidence the decision is made on
    'direction', 'wn_key', 'char', 'obs_intens', 'grade', 'new', 'Icalc', 'code',
    'dif_wn_O-C', 'n_candidates', 'n_accepted', 'cg_without', 'cg_with',
    'cg_gain', 'cg_over_unc', 'held_by', 'held_code', 'code_behind',
    'suggested', 'why',
]

DECISIONS_APPLIED = ('accepted', 'accepted, inflated', 'rejected')
DECISIONS_QUIET = ('stays rejected', 'stays accepted', 'review', '')


class ReviewError(Exception):
    """Something in the files contradicts itself; the run stops."""


# ---------------------------------------------------------------------------
# reading
# ---------------------------------------------------------------------------
def records_of(path):
    """The lines of a fixed-column file, without their line endings."""
    with io.open(path, encoding='utf-8', newline='') as fh:
        return [ln.rstrip('\r\n') for ln in fh]


def lopt_field(record, name):
    first, last = LOPT_FIELDS[name]
    return record[first - 1:last].strip()


def read_lopt(path):
    """``{(low_id, upp_id): (flag, wavenumber, uncertainty)}``."""
    out = {}
    for rec in records_of(path):
        if not rec.strip():
            continue
        key = (lopt_field(rec, 'lower_level'), lopt_field(rec, 'upper_level'))
        if key in out:
            raise ReviewError('%s names the transition %s - %s twice'
                              % (path, key[0], key[1]))
        out[key] = (lopt_field(rec, 'flags'),
                    float(lopt_field(rec, 'wavenumber')),
                    float(lopt_field(rec, 'uncertainty')))
    return out


def read_ids(path):
    """``{level_id: IDEN2 row}`` from IDEN_level_ids.txt.

    This is the only way the two level lists are joined.  An energy match
    loses exactly the levels that have just been worked on, because editing a
    level's position in IDEN2 changes its energy and not its row number.
    """
    out = {}
    with io.open(path, encoding='utf-8', newline='') as fh:
        header = fh.readline()
        if 'IDEN_id' not in header:
            raise ReviewError('%s does not start with the level_id/IDEN_id '
                              'header' % path)
        for ln in fh:
            parts = ln.rstrip('\r\n').split('\t')
            if len(parts) >= 2 and parts[0].strip():
                out[parts[0].strip()] = int(parts[1])
    return out


def read_classifications(path):
    """The classified rows, and an index of them by line and by transition.

    A set whose ``wn`` column *is* the column a line is named by writes no
    ``wn_key`` column - the baseline is such a set - so where the column is
    missing the observed wavenumber is the name, which is what
    ``lineclass_config.toml`` says the default is.

    The row of a resolved hfs companion (grade ``hfs``) names the transition
    its main line is, so it stays out of the index by transition: it is in
    ``by_line`` only.
    """
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        rows = [r for r in csv.DictReader(fh)
                if (r.get('low_id') or '').strip()
                and (r.get('upp_id') or '').strip()]
    by_line = collections.defaultdict(list)
    by_pair = {}
    for r in rows:
        if not (r.get('wn_key') or '').strip():
            r['wn_key'] = r['wn_obs']
        by_line[r['wn_key']].append(r)
        if not is_companion(r):
            by_pair[(r['low_id'], r['upp_id'])] = r
    return rows, by_line, by_pair


def is_companion(row):
    """True for the row of a resolved hfs companion: it is the hfs
    component of its main line's transition, which LOPT has under the main
    line, and is never a record of its own."""
    return ((row.get('grade') or '').strip()
            == hfs_correction.COMPANION_GRADE)


def read_level_energies(rows, levels_path):
    """``{level_id: energy}`` as the classification itself used them.

    The energies are taken from the classification's own ``low_E``/``upp_E``
    columns, so that a Ritz wavenumber worked out here is on exactly the scale
    the table's other rows were judged on.  A level named by no row at all -
    every one of its transitions fell outside every window - is looked up in
    the set's ``LOPT_output_levels.txt`` instead.
    """
    out = {}
    for r in rows:
        for which, column in (('low_id', 'low_E'), ('upp_id', 'upp_E')):
            level, energy = r[which].strip(), (r.get(column) or '').strip()
            if level and energy and level not in out:
                try:
                    out[level] = float(energy)
                except ValueError:
                    pass
    if os.path.exists(levels_path):
        with io.open(levels_path, encoding='utf-8', newline='') as fh:
            fh.readline()
            for ln in fh:
                parts = ln.rstrip('\r\n').split('\t')
                if len(parts) >= 2 and parts[0].strip() not in out:
                    try:
                        out[parts[0].strip()] = float(parts[1])
                    except ValueError:
                        pass
    return out


class CorrectedLines(object):
    """The set's line list, searchable by the wavenumber a line is named by.

    ``Pr3_lines_corrected.xlsx`` holds Sugar's own wavenumber (``own``) beside
    the calibrated one (``own_corr``) and its uncertainty, so a line whose
    transition no Ritz window caught can still be given its corrected
    position.  The search is by the ``own`` value, to within ``WN_MATCH``,
    because that is the only quantity the two sets agree on digit for digit.
    """

    def __init__(self, path):
        import openpyxl
        book = openpyxl.load_workbook(path, read_only=True, data_only=True)
        sheet = book[book.sheetnames[0]]
        source = sheet.iter_rows(values_only=True)
        header = [('' if c is None else str(c).strip()) for c in next(source)]
        for column in (COL_OWN, COL_OWN_CORR, COL_UNC_CORR):
            if column not in header:
                raise ReviewError('%s has no %s column' % (path, column))
        index = dict((name, n) for n, name in enumerate(header))
        self.lines = []
        for record in source:
            own = record[index[COL_OWN]]
            if own is None:
                continue
            self.lines.append({
                'wn_key': str(float(own)),
                'own': float(own),
                'wn_obs': record[index[COL_OWN_CORR]],
                'unc_wn_obs': record[index[COL_UNC_CORR]],
                'obs_intens': record[index[COL_INTENS]]
                if COL_INTENS in index else '',
                'char': record[index[COL_CHAR]]
                if COL_CHAR in index else '',
            })
        self.lines.sort(key=lambda r: r['own'])
        self.wavenumbers = [r['own'] for r in self.lines]
        book.close()

    def match(self, wn, tol=WN_MATCH):
        """The line whose own wavenumber is nearest ``wn``, or None."""
        k = bisect.bisect_left(self.wavenumbers, wn)
        best, distance = None, tol
        for m in (k - 1, k, k + 1):
            if 0 <= m < len(self.lines):
                d = abs(self.wavenumbers[m] - wn)
                if d <= distance:
                    best, distance = self.lines[m], d
        return best


def read_ledger(path):
    """``{(rounded wn_key, low_id, upp_id): row}`` of line_decisions.csv."""
    out = {}
    if not os.path.exists(path):
        return out
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        for r in csv.DictReader(fh):
            try:
                wn = float(r['wn_key'])
            except (KeyError, TypeError, ValueError):
                continue
            out[(round(wn, 2), r['low_id'], r['upp_id'])] = r
    return out


def read_dlv(path):
    """``(records, endings, rows_by_number, wavenumbers, row numbers)``.

    Read through IDEN's own reader, which keeps every line ending exactly as
    it is: these files are CRLF in the working tree and a rewrite must not
    quietly convert them.
    """
    records, endings = IDEN.read_records(path)
    rows = {}
    for index, rec in enumerate(records):
        if len(rec) < sync.DLV_WIDTH or not rec.strip():
            continue
        number = int(rec[sync.DLV_ROW[0]:sync.DLV_ROW[1]])
        rows[number] = (index,
                        float(rec[sync.DLV_WN[0]:sync.DLV_WN[1]]),
                        float(rec[sync.DLV_LAMBDA[0]:sync.DLV_LAMBDA[1]]),
                        float(rec[sync.DLV_UNC[0]:sync.DLV_UNC[1]]))
    order = sorted(rows, key=lambda n: rows[n][1])
    return records, endings, rows, [rows[n][1] for n in order], order


def dlv_row_of(wavenumbers, numbers, wn, tol=WN_MATCH):
    """The row of dlv.dat nearest ``wn``, or None if none is within ``tol``."""
    k = bisect.bisect_left(wavenumbers, wn)
    best, distance = None, tol
    for m in (k - 1, k, k + 1):
        if 0 <= m < len(wavenumbers):
            d = abs(wavenumbers[m] - wn)
            if d <= distance:
                best, distance = numbers[m], d
    return best


def read_trans_holders(path):
    """``{dlv row: [(code, owner, partner)]}`` - what holds each line today."""
    held = collections.defaultdict(list)
    owner = None
    records, _endings = IDEN.read_records(path)
    for rec in records:
        if not rec.strip():
            continue
        if rec.startswith('$'):
            owner = int(IDEN.field(rec, IDEN.TH_INDEX))
            continue
        if owner is None:
            raise ReviewError('%s has a transition row before any level '
                              'header' % path)
        tail = IDEN.assignment(rec)
        if not IDEN.has_line(tail):
            continue
        held[IDEN.obs_row(tail)].append(
            (int(IDEN.field(rec, sync.TR_ICALC)), owner,
             int(IDEN.field(rec, IDEN.TR_PARTNER))))
    return held


# ---------------------------------------------------------------------------
# arithmetic
# ---------------------------------------------------------------------------
def number(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return 0.0


def is_accepted(row):
    return number(row.get('accepted')) >= 0.5


def code_of(icalc):
    """IDEN2's logarithmic intensity, ``round(10*ln(Icalc))``."""
    return int(round(10.0 * math.log(icalc))) if icalc > 0 else None


def center_of_gravity(rows):
    """The predicted-intensity-weighted mean departure from the Ritz value.

    This is what a blend's measured position means: one feature whose center
    lies between the components, weighted by how bright each is predicted to
    be.  Returns None when nothing carries any predicted intensity.
    """
    weight = sum(number(r['calc_intens']) for r in rows)
    if weight <= 0:
        return None
    return sum(number(r['dif_wn_O-C']) * number(r['calc_intens'])
               for r in rows) / weight


def preferred_uncertainty(cg):
    """The inflated uncertainty a center of gravity ``cg`` asks for.

    A value within a factor of 1.5 of 0.5 cm^-1 is reported as 0.5: the
    far-ultraviolet blends reconcile with their Ritz values on one constant
    uncertainty, and a separate value per line would be a fit to noise.
    """
    magnitude = abs(cg)
    if PREFER_UNC / PREFER_FACTOR <= magnitude <= PREFER_UNC * PREFER_FACTOR:
        return PREFER_UNC
    return round(magnitude, 4)


# ---------------------------------------------------------------------------
# the worksheet
# ---------------------------------------------------------------------------
UNUSED_FLAGS = ('P', '(absent)')


def direction_of(base_flag, set_flag):
    """Which way the two sets disagree about a transition, or None.

    ``P`` means printed but not fitted, and a transition absent from a set's
    LOPT input is not fitted there either, so the two are the same answer to
    the only question asked here - does this set use the assignment?  A pair
    that neither set uses, and a pair both use, are not disagreements.
    """
    if base_flag == set_flag:
        return None
    if base_flag in UNUSED_FLAGS and set_flag in UNUSED_FLAGS:
        return None
    return 'rejected in %s' if set_flag in UNUSED_FLAGS \
        else 'accepted in %s'


def suggest(case):
    """The advisory verdict, and the rule that produced it.

    Only the three rules the reviewed half confirms are applied, and only in
    the direction they were confirmed in - towards leaving an assignment
    rejected.  Everything else comes back as ``review``: the numbers do not
    decide it.
    """
    if case['code_behind'] is not None \
            and case['code_behind'] >= CODE_MUCH_BRIGHTER:
        return ('stays rejected',
                'predicted %d IDEN2 units fainter than the transition holding '
                'the line' % case['code_behind'])
    if case['cg_with'] is not None and abs(case['cg_with']) > CG_LIMIT:
        return ('stays rejected',
                'the blend center of gravity would stand %.3f cm-1 from Ritz, '
                'beyond %.2f' % (abs(case['cg_with']), CG_LIMIT))
    if case['n_accepted'] == 0 and case['n_candidates'] <= 1:
        return ('stays rejected',
                'nothing else is classified at this line, so there is no '
                'blend to reconcile it with')
    return ('review', 'the numbers do not decide it; look at the level in '
                      'IDEN2')


def suggest_accepted(case):
    """The advisory verdict for the other direction of disagreement.

    Here the set fits an assignment the baseline does not, so the question is
    whether the new acceptance holds.  Only one thing is said with any
    confidence: if the line's center of gravity, this component included,
    lies within the line's own uncertainty, the acceptance is consistent with
    the fit and nothing needs doing.  Everything else is for the eye.
    """
    if case['cg_with'] is not None and case['unc'] \
            and abs(case['cg_with']) <= case['unc']:
        return ('stays accepted',
                'the center of gravity of the line lies within its own '
                'uncertainty with this component included')
    return ('review', 'the set fits an assignment the baseline does not; '
                      'look at the level in IDEN2')


def borrowed_row(key, base_row, by_line, energies, corrected):
    """A classification row for a transition the set never proposed.

    The corrected wavenumber can put a line outside every Ritz window, and
    then the set's table says nothing at all about the transition - which is
    the hardest case to judge by hand and the one most worth assembling.
    Everything needed is elsewhere:

    * the line itself - its corrected wavenumber and uncertainty, its
      intensity and its character - from another component of the same blend
      in this set's own table, and failing that from the set's line list,
      matched on the wavenumber the line is named by;
    * the predicted intensity from the baseline's table, which proposed the
      transition;
    * the departure from the Ritz value worked out here, against this set's
      level energies, since no row of this set carries it.

    Returns ``(row, provenance)``, or ``(None, why)`` when the line cannot be
    found at all.
    """
    low, upp = key
    where = []
    line = None
    if base_row is not None:
        siblings = by_line.get(base_row['wn_key'])
        if siblings:
            line = dict(siblings[0])
            where.append('the line from this set\'s other component of the '
                         'blend, %s - %s' % (line['low_id'], line['upp_id']))
    if line is None:
        own = None
        if base_row is not None:
            own = number(base_row['wn_obs'])
        if own and corrected is not None:
            found = corrected.match(own)
            if found is not None:
                line = dict(found)
                where.append('the line from %s' % NAME_CORRECTED)
    if line is None:
        return None, ('the transition is absent from this set\'s '
                      'classification and the line could not be found in it '
                      'or in %s either' % NAME_CORRECTED)

    row = {
        'wn_key': line['wn_key'],
        'wn_obs': line['wn_obs'],
        'unc_wn_obs': line['unc_wn_obs'],
        'obs_intens': line.get('obs_intens') or '',
        'char': line.get('char') or '',
        'grade': (base_row or {}).get('grade', ''),
        'new': (base_row or {}).get('new', ''),
        'low_id': low, 'upp_id': upp,
        'accepted': '0',
        'calc_intens': (base_row or {}).get('calc_intens', '') or '0',
    }
    if base_row is not None:
        where.append('the predicted intensity from the baseline\'s table')
    if low in energies and upp in energies:
        ritz = energies[upp] - energies[low]
        row['dif_wn_O-C'] = '%.4f' % (number(row['wn_obs']) - ritz)
        where.append('O-C worked out against this set\'s level energies '
                     '(Ritz %.4f)' % ritz)
    else:
        row['dif_wn_O-C'] = ''
        where.append('no O-C: %s is in no row of this set and in no level of '
                     'its LOPT output'
                     % (low if low not in energies else upp))
    return row, ('the transition is absent from this set\'s classification; '
                 + '; '.join(where))


def build_cases(paths, from_iden2_id, keep_decided):
    """One record per accept/reject disagreement, in review order.

    Both directions are reported.  A transition the baseline fits and this set
    does not is the common case - the corrected uncertainties are tighter, and
    a blend component that passed on Sugar's uncertainty fails on the fitted
    one - but the reverse happens too: the corrected wavenumber can bring a
    line into a window it missed before, and then this set fits an assignment
    the baseline never did.  Both are decisions nobody has taken yet.
    """
    base = read_lopt(paths['base_lopt'])
    theirs = read_lopt(paths['set_lopt'])
    ids = read_ids(paths['ids'])
    rows, by_line, by_pair = read_classifications(paths['classifications'])
    energies = read_level_energies(rows, paths['set_levels'])
    base_by_pair = {}
    if os.path.exists(paths['base_classifications']):
        _r, _l, base_by_pair = read_classifications(
            paths['base_classifications'])
    corrected = None
    ledger = read_ledger(paths['decisions'])
    _recs, _ends, _dlv_rows, dlv_wn, dlv_numbers = read_dlv(paths['dlv'])
    held = read_trans_holders(paths['trans'])
    set_name = os.path.basename(paths['set_dir']) or 'the set'

    cases = []
    for key in set(base) | set(theirs):
        mine, base_record = theirs.get(key), base.get(key)
        flag = '(absent)' if base_record is None else (base_record[0] or '')
        set_flag = '(absent)' if mine is None else (mine[0] or '')
        direction = direction_of(flag, set_flag)
        if direction is None:
            continue
        low, upp = key
        row = by_pair.get(key)
        iden1, iden2 = ids.get(low), ids.get(upp)
        if iden1 is None or iden2 is None:
            raise ReviewError(
                'no IDEN2 row for %s in %s.  Every level the pipeline uses is '
                'in that table; a miss means the lookup is broken, not that '
                'the level is new.'
                % (low if iden1 is None else upp, paths['ids']))

        case = dict.fromkeys(WORKSHEET_COLUMNS, '')
        case.update(level_id1=low, level_id2=upp, iden2_id1=iden1,
                    iden2_id2=iden2, lc_flag=flag or '(blank)',
                    iter_flag=set_flag or '(blank)',
                    direction=direction % set_name)
        borrowed = ''

        if row is None:
            # The corrected wavenumber moved the line out of every Ritz
            # window, so this set's classification never proposed the
            # transition.  Everything the judgment needs is still recoverable
            # from the other component of the blend, the baseline's table and
            # this set's level energies.
            if corrected is None and os.path.exists(paths['corrected_lines']):
                corrected = CorrectedLines(paths['corrected_lines'])
            row, borrowed = borrowed_row(key, base_by_pair.get(key), by_line,
                                        energies, corrected)
            if row is None:
                case.update(n_candidates=0, n_accepted=0,
                            suggested='review', why=borrowed)
                cases.append(case)
                continue

        candidates = by_line.get(row['wn_key'], [])
        if row not in candidates:
            candidates = candidates + [row]
        others = [c for c in candidates if c is not row]
        accepted = [c for c in others if is_accepted(c)]
        cg_without = center_of_gravity(accepted)
        cg_with = center_of_gravity(accepted + [row])
        unc = number(row['unc_wn_obs'])
        icalc = number(row['calc_intens'])
        code = code_of(icalc)

        number_of_row = dlv_row_of(dlv_wn, dlv_numbers, number(row['wn_obs']))
        holders = held.get(number_of_row, []) if number_of_row else []
        held_code = max((c for c, _o, _p in holders), default=None)
        held_by = ' '.join('%d-%d(%d)' % (o, p, c) for c, o, p in
                           sorted(holders, reverse=True))

        case.update(
            obs_wn='%.4f' % number(row['wn_obs']),
            unc_obs_wn='%.4f' % unc,
            wn_key=row['wn_key'],
            char=row['char'], obs_intens='%.1f' % number(row['obs_intens']),
            grade=row['grade'], new=row['new'],
            Icalc='%.0f' % icalc, code='' if code is None else code,
            **{'dif_wn_O-C': '%+.4f' % number(row['dif_wn_O-C'])})
        case.update(
            n_candidates=len(candidates), n_accepted=len(accepted),
            cg_without='' if cg_without is None else '%+.4f' % cg_without,
            cg_with='' if cg_with is None else '%+.4f' % cg_with,
            cg_gain='' if (cg_without is None or cg_with is None)
                    else '%+.4f' % (abs(cg_without) - abs(cg_with)),
            cg_over_unc='' if (cg_with is None or not unc)
                        else '%.2f' % (abs(cg_with) / unc),
            held_by=held_by or '(unassigned)',
            held_code='' if held_code is None else held_code,
            code_behind='' if (held_code is None or code is None)
                        else held_code - code)

        # the numeric forms the suggestion needs, kept off the worksheet
        case['_cg_with'] = cg_with
        case['_unc'] = unc
        evidence = {
            'code_behind': None if (held_code is None or code is None)
                           else held_code - code,
            'cg_with': cg_with, 'unc': unc, 'n_accepted': len(accepted),
            'n_candidates': len(candidates)}
        verdict, why = suggest(evidence) if direction.startswith('rejected') \
            else suggest_accepted(evidence)
        case['suggested'] = verdict
        case['why'] = ('%s; %s' % (borrowed, why)) if borrowed else why
        if verdict == 'review' and cg_with is not None and unc \
                and abs(cg_with) > unc:
            case['unc_wn_inflated'] = '%.4f' % preferred_uncertainty(cg_with)
        cases.append(case)

    decided = 0
    kept = []
    for case in cases:
        if case['iden2_id2'] < from_iden2_id:
            continue
        found = None
        if case['wn_key']:
            found = ledger.get((round(float(case['wn_key']), 2),
                                case['level_id1'], case['level_id2']))
        if found is not None:
            decided += 1
            if not keep_decided:
                continue
            case['Notes'] = 'the ledger already rules on this: %s %s, %s' % (
                found['decision'], found['date'], found['reason'])
        kept.append(case)
    kept.sort(key=lambda c: (c['iden2_id2'], c['level_id1']))
    return kept, decided


def build_free_lines(paths, from_iden2_id, z_max, ic_min):
    """Rejected transitions whose observed line nothing holds in IDEN2.

    A line that no assignment claims is a line still on offer.  Where the
    classification proposed a transition for it, rejected it, and the
    rejection was not about the Ritz agreement, the transition is worth a
    second look at the level in IDEN2 - which is how the assignments already
    recovered by hand were found.

    Two filters keep the list to what is worth looking at.  ``z_max`` is how
    many of the line's own uncertainties its departure from the Ritz value may
    be.  ``ic_min`` is the smallest share of the observed intensity the
    transition may be predicted to carry: the classification accepts a legacy
    Ritz match at 0.04, and a coincidence with a transition a hundred times
    too faint is a coincidence, not evidence.

    The rows come back in the worksheet's own format, so ``--apply`` carries
    them out exactly as it carries out the disagreements.
    """
    ids = read_ids(paths['ids'])
    rows, by_line, _by_pair = read_classifications(paths['classifications'])
    _recs, _ends, _dlv_rows, dlv_wn, dlv_numbers = read_dlv(paths['dlv'])
    held = read_trans_holders(paths['trans'])
    ledger = read_ledger(paths['decisions'])

    out = []
    for r in rows:
        if is_accepted(r):
            continue
        upper = ids.get(r['upp_id'])
        lower = ids.get(r['low_id'])
        if upper is None or lower is None or upper < from_iden2_id:
            continue
        if (round(float(r['wn_key']), 2), r['low_id'], r['upp_id']) in ledger:
            continue
        unc = number(r['unc_wn_obs'])
        icalc = number(r['calc_intens'])
        iobs = number(r['obs_intens'])
        if not unc or icalc <= 0 or iobs <= 0:
            continue
        omc = number(r['dif_wn_O-C'])
        if abs(omc) > z_max * unc or icalc < ic_min * iobs:
            continue
        number_of_row = dlv_row_of(dlv_wn, dlv_numbers, number(r['wn_obs']))
        if number_of_row is None or held.get(number_of_row):
            continue                 # another assignment already holds it

        candidates = by_line[r['wn_key']]
        case = dict.fromkeys(WORKSHEET_COLUMNS, '')
        case.update(
            obs_wn='%.4f' % number(r['wn_obs']),
            unc_obs_wn='%.4f' % unc, wn_key=r['wn_key'],
            level_id1=r['low_id'], level_id2=r['upp_id'],
            iden2_id1=lower, iden2_id2=upper,
            lc_flag='', iter_flag='P',
            char=r['char'], obs_intens='%.1f' % iobs, grade=r['grade'],
            new=r['new'], Icalc='%.0f' % icalc, code=code_of(icalc),
            n_candidates=len(candidates), n_accepted=0,
            cg_with='%+.4f' % omc,
            cg_over_unc='%.2f' % (abs(omc) / unc),
            held_by='(unassigned)', suggested='review',
            why='the line is free and this transition is a Ritz match at '
                '%.2f sigma carrying %.0f%% of its intensity; the pipeline '
                'rejected it - %s'
                % (abs(omc) / unc, 100.0 * icalc / iobs,
                   (r['notes2'] or r['notes1'] or 'no reason recorded')))
        case.update(**{'dif_wn_O-C': '%+.4f' % omc})
        out.append(case)
    out.sort(key=lambda c: (c['iden2_id2'], c['obs_wn']))
    return out


def write_worksheet(path, cases):
    with io.open(path, 'w', encoding='utf-8', newline='') as fh:
        fh.write('\t'.join(WORKSHEET_COLUMNS) + '\n')
        for case in cases:
            fh.write('\t'.join(str(case[c]) for c in WORKSHEET_COLUMNS) + '\n')


# ---------------------------------------------------------------------------
# applying what the worksheet says
# ---------------------------------------------------------------------------
def read_worksheet(path):
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    missing = [c for c in ('wn_key', 'level_id1', 'level_id2', 'decision',
                           'unc_wn_inflated', 'Notes') if rows
               and c not in rows[0]]
    if missing:
        raise ReviewError('%s has no %s column' % (path, ', '.join(missing)))
    return rows


def append_ledger(path, wanted, log):
    """Add the accept and reject rows, leaving any row already there alone."""
    existing = read_ledger(path)
    today = datetime.date.today().strftime('%m/%d/%Y').lstrip('0').replace(
        '/0', '/')
    added = []
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        header = fh.readline().rstrip('\r\n')
    lines = []
    for wn_key, low, upp, decision, reason in wanted:
        if (round(float(wn_key), 2), low, upp) in existing:
            log('  %s %s - %s: the ledger already rules on it' %
                (wn_key, low, upp))
            continue
        lines.append([wn_key, low, upp, decision, today, reason])
        added.append((wn_key, low, upp, decision))
    if lines:
        with io.open(path, 'a', encoding='utf-8', newline='') as fh:
            writer = csv.writer(fh, lineterminator='\n')
            for row in lines:
                writer.writerow(row)
    log('  %s: %d rows added (header %r)' % (os.path.basename(path),
                                             len(lines), header))
    return added


def append_inflated(path, wanted, log):
    """Add the widened uncertainties, in the registry's own tab-delimited form.

    The wavenumber is written to four decimals, the precision
    ``hfs_kappa.inflated_key`` gives.  A line the registry already lists -
    under any precision, an older entry of two or three decimals included -
    is left as it stands.
    """
    registry = hfs_kappa.read_inflated(path)
    today = datetime.date.today().strftime('%m/%d/%Y').lstrip('0').replace(
        '/0', '/')
    added = []
    with io.open(path, 'a', encoding='utf-8', newline='') as fh:
        for wn_key, unc, reason in wanted:
            key = hfs_kappa.inflated_key(float(wn_key))
            there = registry.key_of(float(wn_key))
            if there is not None:
                log('  %s is already in the registry%s; left as it stands'
                    % (key, '' if there == key else ' as ' + there))
                continue
            fh.write('%s\t%s\t%s\t%s\n' % (key, unc, today, reason))
            registry[key] = float(unc)
            added.append((key, unc))
    log('  %s: %d rows added' % (os.path.basename(path), len(added)))
    return added


def update_lopt_input(path, by_line, by_pair, verdicts, inflations, w_hfs,
                      log, dry_run, u_ln_blend=0.0):
    """Carry the decisions into the set's own LOPT transitions file.

    Every record of an observed line is rewritten together, because the three
    things a decision changes are all properties of the line as a whole:

    * the **flag**, which is empty for an accepted classification and ``P``
      for one that is not.  ``P`` tells LOPT to print the line and not fit it;
    * the **weight**, which the accepted classifications of one observed line
      divide between them in proportion to their calculated intensities.
      Accepting a second component therefore changes the weight of the first,
      and rejecting one gives the whole line back to what is left.  The rule
      is ``make_LOPT_input.blend_weights``, so the file says what a full
      regeneration would say;
    * the **uncertainty**, which the records of one line share.  It is the
      widened value where the decision widens it, the level widths added in
      quadrature by ``make_LOPT_input.total_unc``, and then the
      weight-weighted mean across the line by
      ``make_LOPT_input.blend_uncertainty``, widened in quadrature by the
      uncertainty of the blend's centroid, ``make_LOPT_input.centroid_unc``
      with ``u_ln_blend`` ([blends] u_ln_intensity) - the same steps, in the
      same order, that wrote the file in the first place.

    A transition that has no record here cannot be given one, because the
    classification has none either: there is no observed intensity and no
    predicted wavenumber to write.  It is reported instead, and a
    classification run will put it in.
    """
    records, endings = IDEN.read_records(path)
    index = {}
    for k, rec in enumerate(records):
        if rec.strip():
            index[(lopt_field(rec, 'lower_level'),
                   lopt_field(rec, 'upper_level'))] = k

    # which observed lines a decision touches, and what it says about each
    # transition of them
    lines, absent = {}, []
    for (low, upp), verdict in verdicts.items():
        row = by_pair.get((low, upp))
        if row is None:
            absent.append((low, upp, 'no row in the classification table'))
            continue
        lines.setdefault(row['wn_obs'], {})[(low, upp)] = verdict

    changed = []
    for wn_obs, decided in sorted(lines.items(), key=lambda kv: -float(kv[0])):
        group = [r for r in by_line[by_pair[list(decided)[0]]['wn_key']]
                 if r['wn_obs'] == wn_obs]
        # a blended companion's hfs component keeps its share of the line
        other = sum(number(r['calc_intens'] or 0) for r in group
                    if is_companion(r))
        group = [r for r in group if not is_companion(r)]
        accepted = []
        for r in group:
            verdict = decided.get((r['low_id'], r['upp_id']))
            if verdict == 'accept' or (verdict is None and is_accepted(r)):
                accepted.append(r)
        weights = make_LOPT_input.blend_weights(accepted) if accepted else {}

        base = inflations.get(wn_obs)
        per_row = {}
        for r in group:
            quoted = base if base is not None else number(r['unc_wn_obs'])
            per_row[id(r)] = make_LOPT_input.total_unc(
                quoted, r['low_id'].strip(), r['upp_id'].strip(), w_hfs)
        shared = None
        if len(group) > 1:
            shared = make_LOPT_input.blend_uncertainty(
                [per_row[id(r)] for r in group],
                [weights.get(id(r), 0.0) for r in group])
            shared = math.hypot(shared, make_LOPT_input.centroid_unc(
                group, weights, u_ln_blend))
        share = make_LOPT_input.accepted_share(accepted, other)

        for r in group:
            key = (r['low_id'], r['upp_id'])
            k = index.get(key)
            if k is None:
                absent.append((key[0], key[1], 'no record in %s'
                               % os.path.basename(path)))
                continue
            unc = (shared if shared is not None else per_row[id(r)]) / share
            if id(r) in weights:
                flag, weight = '', weights[id(r)]
            else:
                flag, weight = 'P', 0.0
            new = make_LOPT_input.format_line(
                number(r['wn_obs']), unc, number(r['obs_intens']),
                r['low_id'].strip(), r['upp_id'].strip(), flag, weight)
            old = records[k]
            if new != old:
                records[k] = new
                changed.append((number(r['wn_obs']), key, old, new))

    for wn, key, old, new in changed:
        log('  %12.3f  %s - %s' % (wn, key[0], key[1]))
        log('      was  %s' % old)
        log('      now  %s' % new)
    for low, upp, why in absent:
        log('  %s - %s not changed: %s' % (low, upp, why))
    if changed and not dry_run:
        IDEN.write_records(path, records, endings)
    log('  %s: %d records %s' % (os.path.basename(path), len(changed),
                                 'would change' if dry_run else 'changed'))
    return changed


def update_dlv(path, wanted, log, dry_run):
    """Put the widened uncertainties into dlv.dat, on its wavelength scale.

    ``wanted`` gives wavenumbers in the set's own scale and uncertainties in
    cm^-1; the file states an uncertainty in angstroms, to four decimals, so
    each is converted with the row's own wavelength.
    """
    records, endings, rows, wavenumbers, numbers = read_dlv(path)
    changed = []
    for wn, u_wn in wanted:
        number_of_row = dlv_row_of(wavenumbers, numbers, wn)
        if number_of_row is None:
            log('  %.4f has no row in %s; its uncertainty is not updated'
                % (wn, os.path.basename(path)))
            continue
        index, row_wn, row_lam, u_old = rows[number_of_row]
        u_lam = u_wn * row_lam / row_wn
        records[index] = IDEN.put(records[index], sync.DLV_UNC,
                                  '%13.4f' % u_lam)
        changed.append((number_of_row, row_wn, u_old, u_lam, u_wn))
    for number_of_row, row_wn, u_old, u_lam, u_wn in changed:
        log('  row %6d  %12.3f  %9.4f -> %9.4f A  (%.4f cm-1)'
            % (number_of_row, row_wn, u_old, u_lam, u_wn))
    if changed and not dry_run:
        IDEN.write_records(path, records, endings)
    log('  %s: %d uncertainties %s' % (os.path.basename(path), len(changed),
                                       'would change' if dry_run
                                       else 'changed'))
    return changed


def update_trans(path, dlv_path, wanted, ids, log, dry_run):
    """Assign the newly accepted transitions in trans.dat.

    ``wanted`` gives ``(low_id, upp_id, wavenumber)``.  The assignment names
    the line by its row in dlv.dat, which is what IDEN2 reads, and carries
    the observed wavenumber and its departure from the row's own predicted
    value.
    """
    _recs, _ends, rows, wavenumbers, numbers = read_dlv(dlv_path)
    trans = IDEN.Trans(path)
    assigned, already, absent = [], [], []
    for low, upp, wn in wanted:
        a, b = ids[low], ids[upp]
        k = trans.row(a, b)
        if k is None:
            absent.append((low, upp, 'no row in trans.dat for %d - %d'
                           % (a, b)))
            continue
        number_of_row = dlv_row_of(wavenumbers, numbers, wn)
        if number_of_row is None:
            absent.append((low, upp, '%.4f has no row in dlv.dat' % wn))
            continue
        record = trans.records[k]
        if IDEN.has_line(IDEN.assignment(record)):
            already.append((low, upp, IDEN.obs_row(IDEN.assignment(record))))
            continue
        _index, row_wn, _lam, _u = rows[number_of_row]
        predicted = float(IDEN.field(record, IDEN.TR_WN))
        code = int(IDEN.field(record, sync.TR_ICALC))
        tail = sync.make_assignment(code, row_wn, row_wn - predicted,
                                    number_of_row)
        trans.records[k] = record[:IDEN.TR_OBS] + tail
        assigned.append((low, upp, number_of_row, row_wn, row_wn - predicted))
    for low, upp, row, wn, omc in assigned:
        log('  %s - %s  ->  dlv row %d at %.3f  (O-C %+.3f)'
            % (low, upp, row, wn, omc))
    for low, upp, row in already:
        log('  %s - %s already holds dlv row %d; left as it stands'
            % (low, upp, row))
    for low, upp, why in absent:
        log('  %s - %s not assigned: %s' % (low, upp, why))
    if assigned and not dry_run:
        IDEN.write_records(path, trans.records, trans.endings)
    log('  %s: %d assignments %s' % (os.path.basename(path), len(assigned),
                                     'would be made' if dry_run else 'made'))
    return assigned


def apply_worksheet(paths, worksheet, log, dry_run):
    """Carry out every decision the worksheet states."""
    rows = read_worksheet(worksheet)
    ids = read_ids(paths['ids'])
    _cls, by_line, by_pair = read_classifications(paths['classifications'])
    w_hfs = ({} if not os.path.exists(paths['hfs_widths'])
             else make_LOPT_input.read_hfs_widths(paths['hfs_widths']))

    unknown = sorted({r['decision'].strip() for r in rows} -
                     set(DECISIONS_APPLIED) - set(DECISIONS_QUIET))
    if unknown:
        raise ReviewError(
            '%s: the decision column says %s, which is not one of %s'
            % (worksheet, ', '.join(repr(u) for u in unknown),
               ', '.join(repr(d) for d in DECISIONS_APPLIED + DECISIONS_QUIET)))

    ledger_rows, inflations, dlv_changes, assignments = [], [], [], []
    verdicts, line_unc = {}, {}
    for r in rows:
        decision = r['decision'].strip()
        if decision not in DECISIONS_APPLIED:
            continue
        if not r['wn_key'].strip():
            raise ReviewError(
                '%s - %s is marked %r but carries no wn_key, so no ledger row '
                'can be keyed to it' % (r['level_id1'], r['level_id2'],
                                        decision))
        reason = r['Notes'].strip() or 'IDEN2/LOPT'
        verdict = 'reject' if decision == 'rejected' else 'accept'
        ledger_rows.append((r['wn_key'], r['level_id1'], r['level_id2'],
                            verdict, reason))
        verdicts[(r['level_id1'], r['level_id2'])] = verdict
        if verdict == 'accept':
            assignments.append((r['level_id1'], r['level_id2'],
                                float(r['obs_wn'])))
        unc = r['unc_wn_inflated'].strip()
        if unc:
            inflations.append((r['wn_key'], unc, reason))
            dlv_changes.append((float(r['obs_wn']), float(unc)))
            classified = by_pair.get((r['level_id1'], r['level_id2']))
            if classified is not None:
                held = line_unc.setdefault(classified['wn_obs'], float(unc))
                if held != float(unc):
                    raise ReviewError(
                        'the worksheet gives the line at %s two different '
                        'inflated uncertainties, %s and %s.  One observed '
                        'line has one uncertainty, which its components '
                        'share.' % (r['obs_wn'], held, unc))

    log('%d decisions to carry out: %d ledger rows, %d inflations, '
        '%d assignments' % (len(ledger_rows), len(ledger_rows),
                            len(inflations), len(assignments)))
    if not ledger_rows:
        return

    targets = [paths['decisions'], paths['inflated'], paths['set_lopt'],
               paths['dlv'], paths['trans']]
    if not dry_run:
        output_files.require_writable(targets, 'file this decision changes')

    log('line_decisions.csv')
    if dry_run:
        already = read_ledger(paths['decisions'])
        for wn, low, upp, verdict, reason in ledger_rows:
            there = (round(float(wn), 2), low, upp) in already
            log('  %s  %s - %s  %s  %s'
                % (wn, low, upp, verdict,
                   'the ledger already rules on it' if there else reason))
    else:
        append_ledger(paths['decisions'], ledger_rows, log)

    log('inflated_unc_lines.txt')
    if dry_run:
        for wn, unc, reason in inflations:
            log('  %s  %s  %s' % (hfs_kappa.inflated_key(float(wn)), unc,
                                  reason))
    else:
        append_inflated(paths['inflated'], inflations, log)

    log('LOPT_input_lines.txt')
    # a directory without a configuration (a test's) has no [blends] term
    u_ln_blend = (config.load(paths['config']).blend_u_ln
                  if os.path.isfile(paths['config']) else 0.0)
    update_lopt_input(paths['set_lopt'], by_line, by_pair, verdicts, line_unc,
                      w_hfs, log, dry_run, u_ln_blend=u_ln_blend)

    log('IDEN2/dlv.dat')
    update_dlv(paths['dlv'], dlv_changes, log, dry_run)

    log('IDEN2/trans.dat')
    update_trans(paths['trans'], paths['dlv'], assignments, ids, log, dry_run)


# ---------------------------------------------------------------------------
# the command line
# ---------------------------------------------------------------------------
def resolve(set_name):
    """Where every file this tool reads and writes lives."""
    set_dir = os.path.join(HERE, set_name) if set_name else HERE
    if not os.path.isdir(set_dir):
        raise ReviewError('%s is not a directory' % set_dir)
    iden2 = working_path(NAME_IDEN2, cwd=set_dir)
    return {
        'set_dir': set_dir,
        'base_lopt': os.path.join(HERE, NAME_LOPT_IN),
        'set_lopt': os.path.join(set_dir, NAME_LOPT_IN),
        'classifications': working_path(NAME_CLASSIFICATIONS, cwd=set_dir),
        'decisions': working_path(NAME_DECISIONS, cwd=set_dir),
        'inflated': working_path(NAME_INFLATED, cwd=set_dir),
        'hfs_widths': working_path(make_LOPT_input.DEF_HFS_WIDTHS,
                                   cwd=set_dir),
        'base_classifications': os.path.join(HERE, NAME_CLASSIFICATIONS),
        'set_levels': os.path.join(set_dir, NAME_LOPT_LEVELS),
        'corrected_lines': working_path(NAME_CORRECTED, cwd=set_dir),
        'ids': os.path.join(iden2, 'IDEN_level_ids.txt'),
        'dlv': os.path.join(iden2, 'dlv.dat'),
        'trans': os.path.join(iden2, 'trans.dat'),
        'worksheet': os.path.join(set_dir, NAME_WORKSHEET),
        'config': working_path(config.CONFIG_NAME, cwd=set_dir),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Work through the accept/reject disagreements between the '
                    'baseline LOPT input and another set\'s.')
    ap.add_argument('--set', default='iter', metavar='NAME',
                    help='the set to review against the baseline '
                         '(default: %(default)s)')
    ap.add_argument('--out', default=None, metavar='FILE',
                    help='where the worksheet goes (default: '
                         '<set>/decisions.txt)')
    ap.add_argument('--from-iden2-id', type=int, default=0, metavar='N',
                    dest='from_iden2_id',
                    help='start at this row of IDEN2/enlev.dat, for picking '
                         'the review up where it was left')
    ap.add_argument('--all', action='store_true',
                    help='keep the disagreements the ledger already rules on, '
                         'which are otherwise left out as decided')
    ap.add_argument('--free-lines', action='store_true', dest='free_lines',
                    help='instead of the disagreements, list the rejected '
                         'transitions whose observed line no assignment holds '
                         'in IDEN2 - the ones worth a second look at the level')
    ap.add_argument('--z-max', type=float, default=2.0, metavar='Z',
                    dest='z_max',
                    help='with --free-lines, how many of the line\'s own '
                         'uncertainties its departure from the Ritz value may '
                         'be (default: %(default)s)')
    ap.add_argument('--ic-min', type=float, default=0.04, metavar='R',
                    dest='ic_min',
                    help='with --free-lines, the smallest share of the '
                         'observed intensity the transition may be predicted '
                         'to carry (default: %(default)s)')
    ap.add_argument('--apply', action='store_true',
                    help='read the worksheet back and carry out its decision '
                         'column')
    ap.add_argument('--dry-run', action='store_true', dest='dry_run',
                    help='with --apply, say what would be written and write '
                         'nothing')
    args = ap.parse_args(argv)

    def log(text=''):
        print(text)

    try:
        paths = resolve(args.set)
        # dlv.dat is read and written here by the measured wavenumber of
        # each line, so a file showing lines where LOPT is given them
        # ([hfs] iden2_display = 'lopt') would lose them.
        if sync.read_shown(os.path.dirname(paths['dlv'])):
            raise ReviewError(
                '%s shows lines as LOPT is given them (%s lists them).  Set '
                "[hfs] iden2_display = 'measured' and run sync_IDEN2.py in "
                'the set first.' % (paths['dlv'], sync.SHOWN_FILE))
        worksheet = args.out or paths['worksheet']
        if args.apply:
            if not os.path.exists(worksheet):
                raise ReviewError('%s does not exist; run without --apply '
                                  'first' % worksheet)
            log('reading the decisions of %s' % worksheet)
            apply_worksheet(paths, worksheet, log, args.dry_run)
            return 0

        for what in ('base_lopt', 'set_lopt', 'classifications', 'ids',
                     'dlv', 'trans'):
            if not os.path.exists(paths[what]):
                raise ReviewError('%s does not exist' % paths[what])
        if args.free_lines:
            cases, decided = build_free_lines(
                paths, args.from_iden2_id, args.z_max, args.ic_min), 0
            what = 'free lines to look at'
        else:
            cases, decided = build_cases(paths, args.from_iden2_id, args.all)
            what = 'disagreements to review'
        write_worksheet(worksheet, cases)

        counts = collections.Counter(c['suggested'] for c in cases)
        log('%s: %d %s' % (worksheet, len(cases), what))
        if decided:
            log('  %d more already carry a ledger row and are left out%s'
                % (decided, ' (--all keeps them)' if not args.all else ''))
        for verdict, n in counts.most_common():
            log('  %-16s %d' % (verdict, n))
        log()
        log('Fill the decision column with one of: %s.'
            % ', '.join(DECISIONS_APPLIED))
        log('Leave it empty, or write "stays rejected", to change nothing.')
        log('Then: python %s --set %s --apply'
            % (os.path.basename(__file__), args.set))
    except ReviewError as exc:
        sys.stderr.write('review_mismatches.py: %s\n' % exc)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
