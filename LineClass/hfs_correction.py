"""The head-frame hyperfine correction of the observed wavenumbers.

This is Step 5 of `Work_on_hfs_plan.md`, as agreed with the user on
2026-09-29.  It is switched on by `apply = true` in the `[hfs]` section of a
working set's `lineclass_config.toml`; with it off, nothing in the pipeline
changes.

THE PHYSICS, IN ONE PARAGRAPH
=============================
141Pr has nuclear spin I = 5/2.  A level of angular momentum J is therefore a
group of hyperfine sublevels F = |J-I| ... J+I, spread about the level's center
of gravity by its magnetic-dipole constant A (cm^-1).  The outermost sublevel,
F = I+J, lies

    S = I * A * J

from the center of gravity (above it for A > 0, below it for A < 0).  A line
joining two levels is a pattern of components; its strongest component joins
the two F = I+J sublevels and lies

    D = S(upper) - S(lower)

from the pattern's center of gravity.

THE HEAD FRAME
==============
The level energies of the pipeline, of LOPT and of IDEN2 are taken to be the
F = I+J sublevels, E_head = E_cg + S.  That is what Sugar's resolved lines
measure - he tabulated the strongest component (plan, D7) - and what the
pipeline has always nearly fitted.  The Ritz wavenumber E_head(upper) -
E_head(lower) is then the position of the strongest component.

A line whose pattern Sugar did not resolve was measured at a point between the
center of gravity and the strongest component.  Step 3 (`hfs_kappa.py`)
measured where, per class of line, as the convention factor kappa:

    wn_measured = wn_cg + kappa * D = Ritz_head - (1 - kappa) * D

kappa = 1 for the lines he flagged `*r` or `*v` - those sit exactly on the
head-frame Ritz value - 0.944 for the plain lines of 1974, 0.633 for the plain
lines of 1969 and 0.20 for the lines he called complex.  The classes and
their values live in `[hfs.kappa]`.

So a candidate transition i of an observed line is predicted at

    Ritz_i - (1 - kappa) * D_i

by `classify_lines.py`, and an observed line whose accepted transitions have
the branching fractions BF_i (their shares of the calculated intensity, the
`BF` column, the weights LOPT is given) is corrected into the head frame by

    hfs_shift = (1 - kappa) * sum_i BF_i * D_i

which `make_LOPT_input.py` adds to the wavenumber LOPT is given.  One shift per
observed line, because LOPT's centroid blend model compares the one measured
wavenumber of a blend with the BF-weighted mean of its components' Ritz
values, and it recognizes the components of one blend by their sharing that
wavenumber.  A line with no accepted transition is not corrected: an
unidentified line cannot be, and need not be, since nothing is fitted to it.

A flagged blend is the one line whose components do not all share its kappa
(user's decision of 2026-09-30).  Sugar's flag describes the resolved pattern
of one transition, taken to be the blend's strongest (largest BF), which
takes kappa = 1; the others take the kappa the line would have without the
flag, the plain kappa of its era.  Each component then carries its own
(1 - kappa_i) * D_i:

    hfs_shift = sum_i BF_i * (1 - kappa_i) * D_i

since the measured wavenumber is the BF-weighted mean of where each component
was measured.  The line's kappa is reported as the BF-weighted mean,
sum_i BF_i * kappa_i; the two forms agree whenever the D_i do.

IDEN2 needs nothing: `enlev.dat` holds the head-frame energies LOPT fits, and
`dlv.dat` the lines as measured.  The centers of gravity, E_cg = E_head - S,
are computed once at the end, for the published level list and the fit of the
Cowan parameters; that is the only place the absolute A values enter.

THE A CONSTANTS, AND WHAT IS MISSING FROM THEM
==============================================
`A_hfs_levels.csv` (files.hfs_A_levels) gives A, its uncertainty u_A and the
source of the value for each level it lists.  A level falls in one of three
groups:

* **determined** - A comes from Reader and Sugar's calculated values, from the
  level's composition, or from its flag intervals.  S = I*A*J, with the
  uncertainty I*J*u_A;
* **undetermined** - listed, with the source "not determined" (the 41 4f5d2
  and 4f3 levels of the plan's section 3.1), or with a source reporting a
  conflict.  S = 0, and the unknown A is carried as the uncertainty I*J*u_A
  (plan, D10);
* **absent** - not listed.  Nothing is known, S = 0 and no uncertainty is
  added; the empirical hyperfine width of `level_hfs_widths.csv` still
  covers such a level in `make_LOPT_input.py`.

A level in `[hfs] resolved_levels` counts as S = 0 whatever the table says:
its sublevels are resolved, so every line of it ends on one sublevel, the one
its energy stands at, and there is no unresolved pattern of its own to
correct for.  Two J = 1/2 levels of 4f2.6s so far, each standing at its
F = 3 sublevel: 059003.000642 (IDEN2 row 1093, placed there on 2026-09-29)
and 059003.000190 (IDEN2 row 1089).

A level whose lines show a resolved second component (the next rung of the
pattern, F_low - 1 -> F_upp - 1) does not belong here: both levels change F
along that ladder, so the tabulated line is the strongest component of the
whole pattern - a property of the line, as Sugar's flags are, not of the
level.

RESOLVED COMPANIONS (files.hfs_satellites)
==========================================
Where the line list itself resolves a pattern - a weaker line a few tenths of
a cm^-1 from a classified one, at the offset of the next rung of that
transition's pattern (hfs_patterns.rung) - both lines are the one transition.
`hfs_satellites.txt`, kept by hand, names them: the companion line (wn_key),
the line it belongs to (main_wn_key), the transition (low_id, upp_id) and the
rung k, the companion joining F_low - k and F_upp - k.  It rules in two ways:

* the companion is classified as that transition's hfs component, grade
  `hfs` (COMPANION_GRADE), never accepted: it is not a second measurement of
  the energy difference, so it is kept out of LOPT and out of IDEN2, and no
  other candidate is sought for it.  This holds whether `apply` is on or off;
* with `apply` on, the named transition of the main line takes kappa = 1,
  as a flagged line's does: the pattern is resolved, so the main line is its
  strongest component, on the head-frame Ritz value.  In a blend the other
  components keep the kappa they would have without it (the plain one of the
  era, if the line is also flagged).  The main line gives up its level
  widths in make_LOPT_input.py as a flagged line does.

The head of a pattern may be missing from the line list while weaker rungs
are there (main_wn_key left empty, 2026-10-02).  The companions are then
classified as above, grade hfs and never accepted, and no line takes
kappa = 1 for the transition.

A companion cannot also be an identification of its own (attach_hfs_satellites
stops the run), except a published one the ledger rejects, which is withdrawn.

A companion may be blended with a transition of its own: the weaker part of
a line Sugar identified, or one the classification accepts (column `blend` =
1, 2026-10-02).  Such a line is classified as any other, published and
ledger identifications included, and its hfs component is written beside
them as an extra row, grade hfs, never accepted.  The component takes its
share of the line's calculated intensity - the main transition's, times the
fraction of the pattern's strength in that rung (`rung_share`) - so the
line's other transitions carry only the rest of the BF: in the
classification's weights, and in LOPT's as an uncertainty divided by that
rest (make_LOPT_input.accepted_share).  Their hfs shift is that of
the line without the component: it does not move them.  The main line takes
kappa = 1 as above.

WHAT A CORRECTED LINE NO LONGER NEEDS (plan, D13)
=================================================
The corrections replace two allowances that were made for hfs before it could
be computed, and counting both would count it twice:

* a level whose A is determined no longer carries its hyperfine width from
  `level_hfs_widths.csv` into LOPT's uncertainty (`make_LOPT_input.py`);
* an observed line whose accepted transitions touch only such levels no
  longer carries the extra uncertainty of an `inflated_unc_lines.txt` row
  tagged hfs (`classify_lines.py`).

A level without a determined A keeps both, but an *undetermined* level's
unknown hfs is then carried once (2026-10-01).  Its u_S = I*J*u_A already
enters u_hfs_shift, so `classify_lines.py` writes the part of u_hfs_shift owed
to such levels as `u_hfs_undet`, and `make_LOPT_input.py` adds the larger of
that part and the levels' widths, not both.  A line that keeps a registry
allowance tagged hfs (`hfs_allowance` = 1) adds neither: the hand-set value is
the allowance for the hfs of those levels.  An absent level keeps its width in
every case, since it has no u_S to stand in for it.
"""
import collections
import csv
import math
import os
import re

import hfs_kappa
import hfs_patterns

#: nuclear spin of 141Pr.
I_SPIN = hfs_kappa.I_SPIN

#: the source of a listed level whose A is not known (A_hfs_levels.csv).
UNDETERMINED = 'not determined'

#: a registry row (inflated_unc_lines.txt) whose reason contains this word is
#: an allowance for hyperfine structure.
HFS_TAG = 'hfs'

#: the record make_LOPT_input.py leaves beside LOPT_input_lines.txt of what it
#: added to each observed wavenumber.  The readers that match LOPT's records
#: to the observed lines by wavenumber (sync_IDEN2.py, check_sync.py) take the
#: measured value back from it.
SHIFTS_FILE = 'LOPT_hfs_shifts.txt'
SHIFTS_COLUMNS = ('low_id', 'upp_id', 'wn_obs', 'wn_lopt', 'hfs_shift',
                  'u_hfs_shift')

#: the grade classify_lines.py gives the row of a resolved hfs companion,
#: which make_LOPT_input.py leaves out of the LOPT input.
COMPANION_GRADE = 'hfs'

#: the columns of the registry of resolved companions (files.hfs_satellites).
SATELLITE_COLUMNS = ('wn_key', 'main_wn_key', 'low_id', 'upp_id', 'rung',
                     'blend', 'date', 'reason')


def is_determined(source):
    """True if a row of A_hfs_levels.csv with this source gives a usable A."""
    return source != UNDETERMINED and 'CONFLICT' not in source


class Model(object):
    """The correction of one working set: A constants, kappas, resolved levels.

    `settings` is the `hfs` attribute of a `config.Config`.  The model is
    only built when `settings.apply` is on; building it reads the table of A
    constants.
    """

    def __init__(self, settings):
        self.kappa = dict(settings.kappa)
        self.resolved = set(settings.resolved_levels)
        self.satellites = read_satellites(settings.satellites)
        self.levels = {}     # level_id -> (J, A, u_A, determined)
        with open(settings.A_levels, encoding='utf-8', newline='') as fh:
            for row in csv.DictReader(fh):
                src = row['source']
                self.levels[row['level_id'].strip()] = (
                    float(row['J']), float(row['A_cm-1']), float(row['u_A']),
                    is_determined(src))

    # --- levels -------------------------------------------------------------
    def S(self, level_id):
        """`(S, u_S)` of a level, cm^-1: the head sublevel's offset from the
        center of gravity, and its uncertainty.  See the module docstring for
        the three groups of levels."""
        if level_id in self.resolved:
            return 0.0, 0.0
        rec = self.levels.get(level_id)
        if rec is None:
            return 0.0, 0.0
        J, A, u_A, determined = rec
        return (I_SPIN * A * J if determined else 0.0), I_SPIN * J * u_A

    def is_corrected(self, level_id):
        """True if the level's hyperfine structure is accounted for by the
        correction: its A is determined, or its sublevels are resolved."""
        if level_id in self.resolved:
            return True
        rec = self.levels.get(level_id)
        return rec is not None and rec[3]

    def is_undetermined(self, level_id):
        """True if the level is listed without a usable A: its S is 0 and
        its u_S the unknown A (see the module docstring's three groups)."""
        if level_id in self.resolved:
            return False
        rec = self.levels.get(level_id)
        return rec is not None and not rec[3]

    # --- lines --------------------------------------------------------------
    def kappa_of(self, char, wn_key):
        """`(kappa, u_kappa)` of an observed line with Sugar's character
        `char`, named by `wn_key` (his own wavenumber, which says the era)."""
        cls = hfs_kappa.kappa_class(char or '')
        if cls == 'plain':
            cls = 'plain_%d' % hfs_kappa.era_of(wn_key)
        return self.kappa[cls]

    def component_kappas(self, char, wn_key, bfs, heads=None):
        """`[(kappa_i, u_kappa_i)]`, one per accepted component of an observed
        line, `bfs` being their shares of the line (BF) and `heads` saying,
        per component, whether the registry of resolved companions makes the
        line that transition's strongest component.

        Every component takes the kappa of the line's class, with two
        exceptions.  A component `heads` marks takes the flag's kappa, and
        the others of its line the kappa the line would have without that
        mark (the plain kappa of its era if the line is also flagged).
        Otherwise, in a flagged blend the strongest component (largest BF,
        the first of equals) takes the flag's kappa and the others the plain
        kappa of the line's era, the one the line would have without the
        flag.
        """
        k = self.kappa_of(char, wn_key)
        flagged = hfs_kappa.kappa_class(char or '') == 'flag'
        plain = self.kappa['plain_%d' % hfs_kappa.era_of(wn_key)]
        if heads and any(heads):
            other = plain if flagged else k
            return [self.kappa['flag'] if h else other for h in heads]
        if len(bfs) < 2 or not flagged:
            return [k] * len(bfs)
        top = max(range(len(bfs)), key=lambda i: bfs[i])
        return [k if i == top else plain for i in range(len(bfs))]

    def is_head_line(self, wn_key):
        """True if the registry of resolved companions names the line
        `wn_key` as the one its companions belong to."""
        return bool(self.satellites.head_pairs(wn_key))

    def D(self, low_id, upp_id):
        """`(D, u_D)` of the transition, cm^-1."""
        s_l, u_l = self.S(low_id)
        s_u, u_u = self.S(upp_id)
        return s_u - s_l, math.hypot(u_u, u_l)

    def line_shift(self, char, wn_key, components):
        """`(hfs_shift, u_hfs_shift, kappa)` of one observed line, cm^-1.

        `components` are `(low_id, upp_id, BF)` of the line's accepted
        transitions; see `combine` for the three values.  A line with none
        is not shifted and has the kappa of its class.
        """
        if not components:
            return 0.0, 0.0, self.kappa_of(char, wn_key)[0]
        heads = self.satellites.head_pairs(wn_key)
        kappas = self.component_kappas(
            char, wn_key, [bf for _, _, bf in components],
            [(low, upp) in heads for low, upp, _ in components])
        return combine([(bf,) + k + self.D(low, upp)
                        for (low, upp, bf), k in zip(components, kappas)])


def combine(parts):
    """`(hfs_shift, u_hfs_shift, kappa)` of an observed line from its accepted
    components, given as `(BF, kappa, u_kappa, D, u_D)`.

    The shift is sum BF_i * (1 - kappa_i) * D_i, summed class by class of
    kappa so that a line of one class gets exactly (1 - kappa) * sum BF_i*D_i.
    Its uncertainty has two parts, added in quadrature: the A constants,
    sum BF_i * (1 - kappa_i) * u_D_i, summed linearly since the components of
    one line share levels and configurations; and kappa, |sum BF_i * u_kappa_i
    * D_i|, linear too, the kappa of one class being one number.  The kappa
    reported is the BF-weighted mean of the components', or the common one.
    """
    by_kappa = {}
    for bf, kappa, u_kappa, d, u_d in parts:
        acc = by_kappa.setdefault((kappa, u_kappa), [0.0, 0.0, 0.0])
        acc[0] += bf * d
        acc[1] += bf * u_d
        acc[2] += bf
    shift = u_A = u_k = 0.0
    for (kappa, u_kappa), (sum_D, sum_u, _) in by_kappa.items():
        shift += (1.0 - kappa) * sum_D
        u_A += (1.0 - kappa) * sum_u
        u_k += u_kappa * sum_D
    if len(by_kappa) == 1:
        mean_kappa = next(iter(by_kappa))[0]
    else:
        total = sum(acc[2] for acc in by_kappa.values())
        mean_kappa = sum(k * acc[2] for (k, _), acc in by_kappa.items()) / total
    return shift, math.hypot(u_A, u_k), mean_kappa


def is_hfs_reason(reason):
    """True if a registry row with this reason is an allowance for hfs: the
    word HFS_TAG appears in it ("hfs of f26s")."""
    return re.search(r'\b%s\b' % HFS_TAG, reason or '', re.I) is not None


def hfs_registry_keys(path):
    """The registry keys (as `hfs_kappa.Registry` files them) of the rows of
    `inflated_unc_lines.txt` whose reason is an allowance for hfs."""
    out = set()
    if not path or not os.path.exists(path):
        return out
    with open(path, encoding='utf-8', newline='') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            wn = (row.get('wn_key') or '').strip()
            if not wn or wn.startswith('#'):
                continue
            if is_hfs_reason(row.get('reason')):
                out.add(hfs_kappa.registry_key(wn))
    return out


# --- the registry of resolved companions -------------------------------------
Satellite = collections.namedtuple(
    'Satellite', 'key main_key low_id upp_id rung blend reason')


def rung_share(J_low, J_upp, k):
    """The fraction of a line's strength in the k-th rung of its pattern
    (F_low = I + J_low - k -> F_upp = I + J_upp - k; k = 0 is the strongest
    component), 0 beyond the ladder.  It depends on the two J alone."""
    comps = hfs_patterns.components(J_low, 0.0, J_upp, 0.0)
    total = sum(s for _, s, _, _ in comps)
    f1, f2 = I_SPIN + J_low - k, I_SPIN + J_upp - k
    part = sum(s for _, s, a, b in comps
               if abs(a - f1) < 1e-9 and abs(b - f2) < 1e-9)
    return part / total if total else 0.0


class Satellites(object):
    """`read_satellites`'s result: the rows of the registry of resolved hfs
    companions, looked up by the wn_key of a line.

    An entry names a line as `hfs_kappa.Registry` matches them: when the
    line's wn_key, rounded to the number of decimals the entry was written
    with, is the entry.
    """

    def __init__(self, rows=()):
        self.rows = list(rows)
        self._companion = {r.key: r for r in self.rows}
        self._heads = {}
        for r in self.rows:
            if r.main_key is None:          # the head is not observed
                continue
            self._heads.setdefault(r.main_key, set()).add((r.low_id,
                                                           r.upp_id))

    def __len__(self):
        return len(self.rows)

    @staticmethod
    def _find(table, wn):
        for d in range(hfs_kappa.KEY_DECIMALS, -1, -1):
            hit = table.get('%.*f' % (d, wn))
            if hit is not None:
                return hit
        return None

    def companion(self, wn_key):
        """The row naming the line `wn_key` as a companion, or None."""
        return self._find(self._companion, wn_key)

    def head_pairs(self, wn_key):
        """The `(low_id, upp_id)` pairs whose strongest component the line
        `wn_key` is, by the registry; empty if it names none."""
        return self._find(self._heads, wn_key) or frozenset()

    def check(self, wn_keys, source):
        """Raise unless every entry, companion and main line alike, names
        exactly one of the lines `wn_keys`."""
        keys = hfs_kappa.Registry.fromkeys(
            set(self._companion) | set(self._heads), 0.0)
        lines_of, unknown = keys.match(wn_keys)
        if unknown:
            raise ValueError('%s: no observed line has the wn_key %s cm^-1'
                             % (source, ', '.join(unknown)))
        wide = sorted(k for k, v in lines_of.items() if len(v) > 1)
        if wide:
            raise ValueError('%s: an entry names more than one line; write it '
                             'with the full-precision wn_key: %s'
                             % (source, ', '.join(wide)))


def read_satellites(path):
    """The registry of resolved hfs companions (files.hfs_satellites), as
    `Satellites`; see the module docstring.

    Tab-delimited, with the columns of SATELLITE_COLUMNS; a row whose
    wn_key starts with '#' is a comment.  `blend` is 1 for a companion
    blended with a transition of its own, empty or 0 otherwise; a file
    without the column has none.  An empty main_wn_key means the head of the
    pattern is not observed: the row's main_key is None.  A missing file, or
    none named, is an empty registry; a row that cannot mean what it says - a
    line that is its own companion, a companion entered twice, a line entered
    both as a companion and as a main line, a rung that is not a positive
    integer, a blend that is not 0 or 1, a blended companion without a head -
    raises.
    """
    if not path or not os.path.exists(path):
        return Satellites()
    name = os.path.basename(path)
    rows = []
    with open(path, encoding='utf-8', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        missing = [c for c in SATELLITE_COLUMNS[:5]
                   if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError('%s lacks the column(s) %s (header: %s)'
                             % (name, ', '.join(missing), reader.fieldnames))
        for rec in reader:
            wn = (rec.get('wn_key') or '').strip()
            if not wn or wn.startswith('#'):
                continue
            main = (rec.get('main_wn_key') or '').strip()
            low = (rec.get('low_id') or '').strip()
            upp = (rec.get('upp_id') or '').strip()
            rung = (rec.get('rung') or '').strip()
            if not (low and upp) or not rung.isdigit() or int(rung) < 1:
                raise ValueError('%s: the row of %s needs low_id, upp_id '
                                 'and a rung of 1 or more' % (name, wn))
            blend = (rec.get('blend') or '').strip()
            if blend not in ('', '0', '1'):
                raise ValueError('%s: the row of %s has blend = %r; it is 1 '
                                 'for a blended companion, empty or 0 '
                                 'otherwise' % (name, wn, blend))
            if not main and blend == '1':
                raise ValueError('%s: the row of %s has no main_wn_key and '
                                 'blend = 1; a companion whose head is not '
                                 'observed cannot be blended (its share of '
                                 'the line is taken from the head\'s '
                                 'calculated intensity)' % (name, wn))
            rows.append(Satellite(hfs_kappa.registry_key(wn),
                                  hfs_kappa.registry_key(main) if main
                                  else None, low, upp,
                                  int(rung), blend == '1',
                                  (rec.get('reason') or '').strip()))
    seen, mains = set(), {r.main_key for r in rows if r.main_key}
    for r in rows:
        if r.key == r.main_key:
            raise ValueError('%s: %s is entered as its own companion'
                             % (name, r.key))
        if r.key in seen:
            raise ValueError('%s: the companion %s is entered twice'
                             % (name, r.key))
        if r.key in mains:
            raise ValueError('%s: %s is entered both as a companion and as '
                             'the line companions belong to' % (name, r.key))
        seen.add(r.key)
    return Satellites(rows)


# --- the record of what LOPT was given ---------------------------------------
def shifts_path(lopt_input):
    """Where the record of the shifts belongs: beside the LOPT input file."""
    return os.path.join(os.path.dirname(os.path.abspath(lopt_input)),
                        SHIFTS_FILE)


def write_shifts(path, rows):
    """Write the record: `rows` are `(low_id, upp_id, wn_obs, wn_lopt,
    shift, u_shift)`, one per record of the LOPT input whose wavenumber was
    changed.  `wn_obs` is written exactly as given - the classification
    table's own text, when it comes from there - so that reading it back
    gives the very number the table holds.  A transition has one record in
    the LOPT input, so the pair of levels names the record; the wavenumber
    alone would not, since a shifted value can fall on another line's
    measured one."""
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        fh.write('\t'.join(SHIFTS_COLUMNS) + '\n')
        for low, upp, wn, wn_lopt, d, u in sorted(
                rows, key=lambda r: (-float(r[2]), r[0], r[1])):
            wn = wn if isinstance(wn, str) else repr(float(wn))
            fh.write('%s\t%s\t%s\t%.3f\t%+.4f\t%.4f\n'
                     % (low, upp, wn, wn_lopt, d, u))


def read_shifts(lopt_input):
    """`{(low_id, upp_id): (wn_lopt, wn_obs)}` from the record beside
    `lopt_input`, or `{}` if there is none - a LOPT input written with the
    correction off, whose wavenumbers are the measured ones."""
    path = shifts_path(lopt_input)
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, encoding='utf-8', newline='') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            out[(row['low_id'], row['upp_id'])] = (
                round(float(row['wn_lopt']), 3), float(row['wn_obs']))
    return out


def measured(low_id, upp_id, wn_lopt, shifts, tol=0.0005):
    """The measured wavenumber of the LOPT record of `low_id - upp_id` at
    `wn_lopt`.

    A record the file does not list, or lists at a wavenumber more than
    `tol` (cm^-1) away, was not shifted and is returned as it is.  `tol`
    allows for the rounding of a printed value: LOPT prints its output
    wavenumbers to the precision of their uncertainties, so a reader of its
    output passes a wider one.  A record found at the very wavenumber that
    was written gets the measured value exactly, as the table holds it; one
    found at a rounded value gets that value less the shift, so that the
    rounding stays where LOPT put it.
    """
    rec = shifts.get((low_id, upp_id))
    if rec is None or abs(rec[0] - wn_lopt) > tol:
        return wn_lopt
    if abs(rec[0] - wn_lopt) < 1e-6:
        return rec[1]
    return wn_lopt - (rec[0] - rec[1])
