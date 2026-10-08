#!/usr/bin/env python
"""Which lines were measured otherwise than their kappa class says.

WHAT THIS IS FOR
================
`kappa_exceptions.txt` (hfs_correction.py) fixes the measurement convention
kappa of single lines: `head` (kappa = 1, the strongest hyperfine component),
`cg` (kappa = 0, the center of gravity of the pattern) or a class of
[hfs.kappa].  This program looks through a working set for the lines that
belong there, and checks the rows already entered.  It writes nothing but
its own report; every row of the registry is still the analyst's decision.

THE TEST
========
A line of a class measured at kappa sits (1 - kappa) * D below the head-frame
Ritz wavenumber, D being the distance of the pattern's strongest component
from its center of gravity (hfs_correction.py).  The working set's LOPT
fit gives the line in the head frame, at the kappa it was given, and its
residual r = dWO-C.  Had it been given kappa' instead, the residual would be

    r(kappa') = r(kappa) + (kappa - kappa') * D

so three hypotheses are compared on one number: the line's class, the head
and the center of gravity.

The residual of a line the fit has used is partly absorbed into the level
energies it shares.  LOPT reports how much: the uncertainty of the Ritz
wavenumber uWnCstat against the line's own uWnOStat give its leverage

    h = (uWnCstat / uWnOStat) ** 2,

the fraction of the line's own residual the fit has already taken up.  The
residual the line would have if the fit had not used it is r / (1 - h), and
the uncertainty of the Ritz value the other lines give alone is

    u_rest = uWnCstat * uWnOStat / sqrt(uWnOStat**2 - uWnCstat**2).

This is the level support of the test, and it counts precision, not lines.
A level with many assigned lines, only one of them precise, gives that
line h close to 1 and an unbounded u_rest: the line is untestable. Counting
the lines would have called it well supported.

The line's own uncertainty in the test is the calibration's model value
for its character and era (the class table of wavelength_calibration.txt),
never the hand-set value of inflated_unc_lines.txt.  A line was usually
inflated because of the very residual tested here, and testing it on the
inflated value would hide every line the test is for.  The uncertainty of
the hfs shift (the A constants) is added in quadrature.  So

    sigma = sqrt(u_model**2 + u_hfs_shift**2 + u_rest**2),
    z(kappa') = (r(kappa') / (1 - h)) / sigma,

the residual without the line, under kappa', in units of sigma.

A line is testable when sigma / |D|, the uncertainty of the kappa it
measures, is at most --max-u-kappa (0.3).  It is a candidate when its class
is rejected against the better of head and cg by chi^2(class) -
chi^2(best) >= --dchi2 (9, three sigma) and the better one fits within
--fit-sigma (2).  A row already in the registry is audited: it fails when
its own kappa is rejected as a class is, beyond sqrt(--dchi2) (three
sigma), not at --fit-sigma - among some 20 audited rows one or two lie
beyond 2 sigma by chance alone (2026-10-05).  The kappa audited is the
registry's (kappa_held), which differs from the one the fit used
(kappa_now) for a row added after the last chain run.  A held line whose
inflation in inflated_unc_lines.txt already covers its residual, to the
two decimals the advice prints, passes the audit too (2026-10-06): the
pair of rows says all the test can.  The advice itself is left blank when
the value a line needs is within 10 per cent of the registry's.

PARTIAL PATTERNS: INFORMATION, NEVER A VERDICT
==============================================
A line may have been measured on a part of its pattern, one end of it lost:
the head side (`partial:k-`, rungs k to the end) or the tail side
(`partial:0-k`, the head to rung k); see hfs_patterns.partial_pattern.  For
a line that does not fit where it is now (beyond --fit-sigma at its held or
fitted kappa) the part that fits it best is reported with the evidence for
it: where the part left out would lie, a line of the list there
(`line_at_missing`, with what it is assigned to) and any candidate
transition of the classification predicted there (`blend_at_missing`).
About ten parts per pattern, spread 0.2 - 0.3 apart in kappa, mean that
nearly any line fits one by chance, so the residual alone justifies none.
The user's rule (2026-10-05): a part left out must be accounted for by a
blending transition, or the line stays at the cg with an inflated
uncertainty - the plates cannot be inspected and IDEN2 cannot show it.

WHAT IS LEFT OUT
================
Blends (lines with more than one accepted transition): one measured
position cannot give each component its own kappa.  Lines touching a level
whose A is not determined: their D is not known.  Lines LOPT did not fit
(flag P).  They are counted in the summary.

OUTPUT
======
`<set>/kappa_candidates.csv`, one row per tested line, the candidates and
the audited rows first, with both IDEN2 row numbers (IDEN_level_ids.txt)
and the evidence the decision should rest on besides the residual: Sugar's
character, the intensity, the nearest other line of the line list, and the
registry row if the line has one.  A summary by character and era goes to
the console.

    python kappa_candidates.py                  # iter_hfs/
    python kappa_candidates.py --set iter_hfs --dchi2 16
"""
import argparse
import bisect
import collections
import csv
import math
import os
import re
import sys

import config
import cowan_gA
import hfs_correction
import hfs_kappa
import hfs_patterns
import output_files

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(HERE, 'wavelength_calibration.txt')
OUT = 'kappa_candidates.csv'

DCHI2 = 9.0
FIT_SIGMA = 2.0
AGREE = 0.10
MAX_U_KAPPA = 0.3

COLUMNS = ('verdict', 'wn_key', 'char', 'era', 'low_id', 'upp_id',
           'iden_low', 'iden_upp', 'obs_intens', 'D', 'kappa_class',
           'kappa_now', 'held', 'kappa_held', 'kappa_line', 'u_kappa_line',
           'leverage', 'u_rest', 'u_model', 'u_hfs_shift', 'r_loo', 'z_now',
           'z_held', 'z_class', 'z_head', 'z_cg', 'dchi2', 'registry_unc',
           'registry_advice', 'registry_reason',
           'nearest_dwn', 'nearest_intens', 'rung1_dwn', 'line_at_rung1',
           'partial_best', 'kappa_partial', 'z_partial', 'missing_side',
           'missing_share', 'missing_dwn', 'line_at_missing',
           'blend_at_missing')

#: the columns of the best partial pattern and the evidence for it: for
#: information only, never a verdict (see the module docstring).
PARTIAL_COLUMNS = COLUMNS[COLUMNS.index('partial_best'):]

#: how close a line of the list must lie to the predicted rung 1 to be
#: named in `line_at_rung1`, cm^-1.
RUNG_WINDOW = 0.1

#: how close a line of the list, or the predicted position of a candidate
#: transition, must lie to the part of a pattern a partial measurement left
#: out to be named in `line_at_missing` or `blend_at_missing`, cm^-1.
MISSING_WINDOW = 0.1


# --- the model uncertainty -------------------------------------------------
UncClass = collections.namedtuple('UncClass', 'a b')


def read_unc_table(path=REPORT):
    """`{(char, era): UncClass}` from the class table of the calibration
    report (the lines under `char  era  n  a_A  b_cm1 ...`); '(plain)' is
    the empty character.  A class whose constants are nan is left out."""
    with open(path, encoding='utf-8') as fh:
        lines = fh.read().splitlines()
    try:
        start = next(i for i, ln in enumerate(lines)
                     if re.match(r'char\s+era\s+n\s+a_A\s+b_cm1', ln))
    except StopIteration:
        raise SystemExit('%s has no table of uncertainty classes; run '
                         'wavelength_calibration.py first' % path) from None
    table = {}
    for ln in lines[start + 1:]:
        parts = ln.split()
        if len(parts) < 5 or not parts[1].isdigit():
            break
        char = '' if parts[0] == '(plain)' else parts[0]
        a, b = float(parts[3]), float(parts[4])
        if not (math.isnan(a) or math.isnan(b)):
            table[(char, int(parts[1]))] = UncClass(a, b)
    return table


def model_unc(table, char, wn):
    """The calibration's statistical uncertainty of a line of character
    `char` at Sugar's wavenumber `wn`, ignoring the inflation registry: its
    class's two-term model, the pooled class `other` of its era for a
    character without one, Sugar's stated value if neither is there; never
    below `hfs_kappa.character_floor`."""
    era = hfs_kappa.era_of(wn)
    plain = table.get(('', era))
    floor = hfs_kappa.character_floor(char, era, wn,
                                      (plain.a, plain.b) if plain else None)
    u = table.get((char, era)) or table.get(('other', era))
    if u is None:
        return max(hfs_kappa.stated_uncertainty(char, era, wn),
                   hfs_kappa.FLOOR, floor)
    return max(hfs_kappa.two_term(u.a, u.b, wn), hfs_kappa.FLOOR, floor)


# --- the test --------------------------------------------------------------
def leave_one_out(r, u_obs, u_calc):
    """`(h, r_loo, u_rest)` of a line LOPT fitted with the uncertainty
    `u_obs`, its Ritz value then known to `u_calc`, leaving the residual
    `r`.  u_rest is None when the line alone fixes its Ritz value (h of 1
    within rounding): nothing else measures it."""
    h = (u_calc / u_obs) ** 2 if u_obs > 0 else 1.0
    if h >= 0.999:
        return h, None, None
    return (h, r / (1.0 - h),
            u_calc * u_obs / math.sqrt(u_obs ** 2 - u_calc ** 2))


def median(values):
    """The median of a non-empty list.  Not the standard library's: a
    `statistics.py` of the repository root (Mandel-Paule, imported by
    classify_lines.py) shadows it whenever the root is on the path."""
    v = sorted(values)
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def hypotheses(r_loo, kappa_now, kappa_class, D, sigma, kappa_held=None):
    """`{'now', 'class', 'head', 'cg'}` -> z, the residual without the line
    under each kappa in units of sigma; 'now' is the kappa the line was
    fitted at.  With `kappa_held` (a line the registries hold) there is a
    'held' key too: it differs from 'now' when the registry row is newer
    than the fit."""
    def z(k):
        return (r_loo + (kappa_now - k) * D) / sigma
    out = {'now': z(kappa_now), 'class': z(kappa_class), 'head': z(1.0),
           'cg': z(0.0)}
    if kappa_held is not None:
        out['held'] = z(kappa_held)
    return out


def needed_unc(zz, sigma, u_shift, u_rest):
    """The own uncertainty that makes the residual zz * sigma a 1-sigma one
    together with the hfs shift and the rest of the network:
    sqrt(r^2 - u_hfs_shift^2 - u_rest^2), or 0 when those two cover it."""
    r = zz * sigma
    return math.sqrt(max(r ** 2 - u_shift ** 2 - u_rest ** 2, 0.0))


def inflation_covers(reg_u, z, sigma, u_shift, u_rest):
    """Whether the registry value `reg_u` already covers the residual of a
    held line at its held kappa (z['held'], else z['now']): the need, rounded
    to the two decimals registry_advice prints, is at most reg_u."""
    if reg_u is None:
        return False
    need = needed_unc(z.get('held', z['now']), sigma, u_shift, u_rest)
    return round(need, 2) <= reg_u + 1e-9


def registry_advice(reg_u, z, sigma, u_shift, u_rest, fit_sigma=FIT_SIGMA,
                    at=None, agree=AGREE):
    """What the inflation registry's value `reg_u` of a line should become,
    judged at the kappa the line is held at (z['held']) or else fitted at
    (z['now']), with `sigma` the model's own total uncertainty.  `at`
    ('head' or 'cg') judges a candidate where it would be held instead, and
    says so: 'at cg: delete'.

    'delete' when the residual fits the model within `fit_sigma`, the same
    bar a head or cg alternative must clear: the line then needs no own
    value.  Otherwise 'needs X (now Y)', X being the own uncertainty that
    makes the residual a 1-sigma one together with the hfs shift and the
    rest of the network: sqrt(r^2 - u_hfs_shift^2 - u_rest^2) - or nothing
    when X is within `agree` (10 per cent) of Y, which has nothing to change
    (2026-10-06)."""
    if reg_u is None:
        return ''
    zz = z[at] if at else z.get('held', z['now'])
    prefix = 'at %s: ' % at if at else ''
    if abs(zz) <= fit_sigma:
        return prefix + 'delete'
    need = needed_unc(zz, sigma, u_shift, u_rest)
    if abs(need - reg_u) <= agree * reg_u:
        return ''
    return prefix + 'needs %.2f (now %g)' % (need, reg_u)


def advice_of(v, reg_u, z, sigma, u_shift, u_rest, fit_sigma=FIT_SIGMA,
              kappa_class=None):
    """The `registry_advice` of a line whose verdict is `v`: a candidate is
    judged where it would be held; a line that fails its class is judged at
    the class and, after '; ', at the better of head and cg, where a row of
    kappa_exceptions.txt would hold it - its failure is then a measurement
    one, and that value is what the inflation should be."""
    if v.startswith('candidate-'):
        return registry_advice(reg_u, z, sigma, u_shift, u_rest, fit_sigma,
                               at=v[len('candidate-'):])
    out = [registry_advice(reg_u, z, sigma, u_shift, u_rest, fit_sigma)]
    if v == 'class-fails':
        out.append(registry_advice(reg_u, z, sigma, u_shift, u_rest,
                                   fit_sigma,
                                   at=best_alternative(z, kappa_class)))
    return '; '.join(a for a in out if a)


def best_alternative(z, kappa_class=None):
    """'head' or 'cg', whichever fits the line better; a class whose kappa is
    itself 1 (a flagged line) has no head alternative."""
    alts = ['cg'] if kappa_class is not None and kappa_class >= 1.0 \
        else ['head', 'cg']
    return min(alts, key=lambda k: abs(z[k]))


def verdict(z, held, dchi2_min=DCHI2, fit_sigma=FIT_SIGMA, kappa_class=None,
            covered=False):
    """`(verdict, dchi2)` of a testable line.

    `held` is the registry's class name for the line ('' if none).  A held
    line is audited against the kappa the registry holds it at (z['held'],
    or z['now'] when absent): 'audit-fails' when that is rejected beyond
    sqrt(dchi2_min), as a class is, and its inflation does not cover the
    residual (`covered`, inflation_covers), else 'audit-ok'.
    Otherwise 'candidate-head' or 'candidate-cg' when the class loses to the
    better alternative by at least dchi2_min and that one fits within
    fit_sigma; 'class-fails' when the class is rejected and neither
    alternative fits; 'fits-class' else.  A class whose kappa is itself 1
    (a flagged line) has no head alternative."""
    best = best_alternative(z, kappa_class)
    dchi2 = z['class'] ** 2 - z[best] ** 2
    if held:
        return ('audit-ok' if covered or abs(z.get('held', z['now']))
                <= math.sqrt(dchi2_min)
                else 'audit-fails', dchi2)
    if dchi2 >= dchi2_min and abs(z[best]) <= fit_sigma:
        return 'candidate-' + best, dchi2
    if abs(z['class']) > math.sqrt(dchi2_min):
        return 'class-fails', dchi2
    return 'fits-class', dchi2


# --- reading the set -------------------------------------------------------
def _num(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def read_lopt_lines(path):
    """`{(low_id, upp_id): row}` of the lines LOPT fitted (no flag P), with
    the floats this program needs."""
    out = {}
    with open(path, encoding='utf-8', newline='') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if (row.get('F') or '').strip() == 'P':
                continue
            vals = {k: _num(row.get(k)) for k in
                    ('wn_o', 'uWnOStat', 'Wn_c', 'uWnCstat', 'dWO-C',
                     'Weight')}
            if None in vals.values():
                continue
            out[(row['L1'].strip(), row['L2'].strip())] = vals
    return out


def nearest_other(keys, intens, i):
    """`(offset, intensity)` of the line of the list nearest the i-th."""
    best = None
    for j in (i - 1, i + 1):
        if 0 <= j < len(keys):
            d = keys[j] - keys[i]
            if best is None or abs(d) < abs(best[0]):
                best = (d, intens[j])
    return best or (None, None)


def rung1(model, low, upp):
    """Where the first rung of the line's pattern lies from its strongest
    component (hfs_patterns.rung), or None when a level of the transition
    is resolved, whose lines end on one sublevel and have no ladder."""
    if low in model.resolved or upp in model.resolved:
        return None
    J1, A1 = model.levels[low][:2]
    J2, A2 = model.levels[upp][:2]
    return hfs_patterns.rung(1, J1, A1, J2, A2)


def partial_family(model, low, upp, unresolved=()):
    """`[(name, Partial, u_kappa)]`: the parts of the transition's pattern
    a line can be measured on when one end of it is lost - the head side
    (`partial:k-`, k = 1 ... the end of the ladder) or the tail side
    (`partial:0-k`, k = 1 ... one short of it; `partial:0-0`, the head and
    its nearest off-ladder components, is left to the class head).  Empty
    when the transition has no pattern of its own."""
    try:
        J1, A1, _, J2, A2, _ = model.pattern_of(low, upp, unresolved)
    except ValueError:
        return []
    n = hfs_patterns.ladder_length(J1, J2)
    out = []
    for first, last in ([(k, None) for k in range(1, n + 1)]
                        + [(0, k) for k in range(1, n)]):
        try:
            part = hfs_patterns.partial_pattern(J1, A1, J2, A2, first, last)
        except ValueError:
            continue
        u = model.partial_kappa(low, upp, first, last, unresolved)[1]
        out.append((hfs_correction.partial_name(first, last), part, u))
    return out


def best_partial(family, r_loo, kappa_now, D, sigma):
    """`(name, Partial, z)` of the part of `family` that fits the line best,
    the uncertainty of its kappa added to `sigma`; None if there is none."""
    best = None
    for name, part, u in family:
        z = ((r_loo + (kappa_now - part.kappa) * D)
             / math.hypot(sigma, u * D))
        if best is None or abs(z) < abs(best[2]):
            best = (name, part, z)
    return best


def missing_offset(part, kappa_line, D):
    """Where the part of the pattern a partial measurement left out lies,
    from the line as measured, cm^-1: its center of gravity less the
    measured point, both from the pattern's center of gravity."""
    return part.missing[1] - kappa_line * D


def line_near(keys, wn, window=RUNG_WINDOW, exclude=None):
    """The line of the sorted list `keys` nearest `wn`, if within
    `window`; never the line `exclude`."""
    i = bisect.bisect_left(keys, wn)
    near = [keys[j] for j in (i - 2, i - 1, i, i + 1)
            if 0 <= j < len(keys) and keys[j] != exclude]
    best = min(near, key=lambda k: abs(k - wn), default=None)
    return best if best is not None and abs(best - wn) <= window else None


def analyse(set_dir, dchi2_min=DCHI2, fit_sigma=FIT_SIGMA,
            max_u_kappa=MAX_U_KAPPA, report=REPORT):
    """The rows of the output and the counts of the summary."""
    cfg = config.load(os.path.join(set_dir, 'lineclass_config.toml'))
    model = hfs_correction.Model(cfg.hfs)
    unc = read_unc_table(report)
    registry = hfs_kappa.read_inflated()
    reasons = {}
    with open(os.path.join(HERE, hfs_kappa.INFLATED), encoding='utf-8',
              newline='') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            wn = (row.get('wn_key') or '').strip()
            if wn and not wn.startswith('#'):
                reasons[hfs_kappa.registry_key(wn)] = (row.get('reason')
                                                       or '').strip()
    iden = {lid: n for n, lid in cowan_gA.read_id_map(
        os.path.join(set_dir, 'IDEN2', 'IDEN_level_ids.txt')).items()}
    lopt = read_lopt_lines(os.path.join(set_dir, 'LOPT_output_lines.txt'))

    table = list(csv.DictReader(open(
        os.path.join(set_dir, 'line_classifications.csv'), encoding='utf-8',
        newline='')))
    intens_of = {}
    for r in table:
        intens_of.setdefault(float(r['wn_key']), float(r['obs_intens']))
    keys = sorted(intens_of)
    pos = {k: i for i, k in enumerate(keys)}
    intens = [intens_of[k] for k in keys]
    rows_of = collections.defaultdict(list)
    for r in table:
        rows_of[float(r['wn_key'])].append(r)

    def name(low, upp):
        return '%s-%s' % (iden.get(low, low[-6:]), iden.get(upp, upp[-6:]))

    count = collections.Counter()
    rows = []
    for r in table:
        if r['accepted'] != '1.0':
            continue
        if r['n_accepted'] != '1':
            count['blend'] += 1
            continue
        low, upp = r['low_id'], r['upp_id']
        key = float(r['wn_key'])
        if not (model.is_corrected(low) and model.is_corrected(upp)):
            count['A unknown'] += 1
            continue
        lp = lopt.get((low, upp))
        if lp is None:
            count['not fitted by LOPT'] += 1
            continue
        D = float(r['hfs_D'])
        char = r['char']
        h, r_loo, u_rest = leave_one_out(lp['dWO-C'], lp['uWnOStat'],
                                         lp['uWnCstat'])
        if r_loo is None or D == 0.0:
            count['untestable'] += 1
            continue
        u_mod = model_unc(unc, char, key)
        u_shift = float(r['u_hfs_shift'] or 0.0)
        sigma = math.sqrt(u_mod ** 2 + u_shift ** 2 + u_rest ** 2)
        if sigma / abs(D) > max_u_kappa:
            count['untestable'] += 1
            continue
        kappa_now = float(r['kappa'])
        kappa_class = model.kappa_of(char, key)[0]
        exc = model.exception(key, low, upp)
        held = exc.cls if exc else ('head' if (low, upp) in
                                    model.satellites.head_pairs(key) else '')
        kappa_held = (model.row_kappa(exc)[0] if exc else
                      1.0 if held else None)
        if held and abs(kappa_held - kappa_now) > 1e-6:
            # the registry row is newer than the fit, which still has the
            # line at kappa_now: audit it at the kappa the row asks for
            count['held, not yet in the fit'] += 1
        z = hypotheses(r_loo, kappa_now, kappa_class, D, sigma, kappa_held)
        reg_u = registry.lookup(key)
        v, dchi2 = verdict(z, held, dchi2_min, fit_sigma, kappa_class,
                           covered=bool(held) and inflation_covers(
                               reg_u, z, sigma, u_shift, u_rest))
        count[v] += 1
        near = nearest_other(keys, intens, pos[key])
        r1 = rung1(model, low, upp)
        at_r1 = None if r1 is None else line_near(keys, key + r1)
        kappa_line = kappa_now + r_loo / D
        part = best_partial(
            partial_family(model, low, upp, exc.unresolved if exc else ()),
            r_loo, kappa_now, D, sigma)
        miss = {}
        # only a line that does not fit where it is now is worth a part
        if part is not None and abs(z.get('held', z['now'])) > fit_sigma:
            count['partial tested'] += 1
            off = missing_offset(part[1], kappa_line, D)
            at = key + off
            there = line_near(keys, at, MISSING_WINDOW, exclude=key)
            blends = []
            for k in [key] + ([there] if there is not None else []):
                for c in rows_of[k]:
                    if not (c['low_id'] and c['upp_id']) or (
                            k == key and (c['low_id'], c['upp_id'])
                            == (low, upp)):
                        continue
                    pred = k - float(c['dif_wn_O-C'] or 0.0)
                    if abs(pred - at) <= MISSING_WINDOW:
                        blends.append('%s at %+.3f, Icalc %.0f%s' % (
                            name(c['low_id'], c['upp_id']), pred - at,
                            float(c['calc_intens'] or 0.0),
                            ', accepted' if c['accepted'] == '1.0' else ''))
            assigned = [name(c['low_id'], c['upp_id'])
                        for c in rows_of.get(there, ())
                        if c['accepted'] == '1.0']
            miss = {
                'partial_best': part[0],
                'kappa_partial': '%.2f' % part[1].kappa,
                'z_partial': '%+.1f' % part[2],
                'missing_side': part[1].missing[0],
                'missing_share': '%.0f%%' % (100 * part[1].missing[2]),
                'missing_dwn': '%+.3f' % off,
                'line_at_missing': '' if there is None else
                '%.4f (I %.0f%s)' % (there, intens_of[there],
                                     '; ' + ', '.join(assigned)
                                     if assigned else '; unassigned'),
                'blend_at_missing': '; '.join(blends),
            }
            if abs(part[2]) <= fit_sigma and (
                    v in ('class-fails', 'audit-fails')):
                count['partial fits'] += 1
                if there is not None or blends:
                    count['partial fits, evidence'] += 1
        rows.append({
            'verdict': v, 'wn_key': r['wn_key'], 'char': char,
            'era': hfs_kappa.era_of(key), 'low_id': low, 'upp_id': upp,
            'iden_low': iden.get(low, ''), 'iden_upp': iden.get(upp, ''),
            'obs_intens': '%.0f' % float(r['obs_intens']),
            'D': '%+.4f' % D, 'kappa_class': '%.3f' % kappa_class,
            'kappa_now': '%.3f' % kappa_now, 'held': held,
            'kappa_held': '' if kappa_held is None else '%.3f' % kappa_held,
            'kappa_line': '%.2f' % kappa_line,
            'u_kappa_line': '%.2f' % (sigma / abs(D)),
            'leverage': '%.3f' % h, 'u_rest': '%.4f' % u_rest,
            'u_model': '%.4f' % u_mod, 'u_hfs_shift': '%.4f' % u_shift,
            'r_loo': '%+.4f' % r_loo,
            'z_now': '%+.1f' % z['now'],
            'z_held': '%+.1f' % z['held'] if 'held' in z else '',
            'z_class': '%+.1f' % z['class'], 'z_head': '%+.1f' % z['head'],
            'z_cg': '%+.1f' % z['cg'], 'dchi2': '%.1f' % dchi2,
            'registry_unc': '' if reg_u is None else '%g' % reg_u,
            'registry_advice': advice_of(v, reg_u, z, sigma, u_shift,
                                         u_rest, fit_sigma, kappa_class),
            'registry_reason': '' if reg_u is None else reasons.get(
                _registry_name(reasons, key), ''),
            'nearest_dwn': '' if near[0] is None else '%+.3f' % near[0],
            'nearest_intens': '' if near[1] is None else '%.0f' % near[1],
            'rung1_dwn': '' if r1 is None else '%+.3f' % r1,
            'line_at_rung1': '' if at_r1 is None else '%.4f (I %.0f)' % (
                at_r1, intens_of[at_r1]),
            **{c: miss.get(c, '') for c in PARTIAL_COLUMNS},
        })
    order = {'audit-fails': 0, 'candidate-cg': 1, 'candidate-head': 1,
             'class-fails': 2, 'audit-ok': 3, 'fits-class': 4}
    rows.sort(key=lambda x: (order[x['verdict']], -float(x['dchi2'])))
    return rows, count


def _registry_name(reasons, key):
    """The registry key, as written, that names the line `key`."""
    for d in range(hfs_kappa.KEY_DECIMALS, -1, -1):
        name = '%.*f' % (d, key)
        if name in reasons:
            return name
    return None


def summary(rows, count):
    lines = ['tested %d line(s); left out: %s' % (
        len(rows), ', '.join('%s %d' % (k, count[k]) for k in
                             ('blend', 'A unknown',
                              'not fitted by LOPT', 'untestable')))]
    lines.append('verdicts: ' + ', '.join(
        '%s %d' % (k, count[k]) for k in
        ('candidate-cg', 'candidate-head', 'class-fails', 'audit-ok',
         'audit-fails', 'fits-class')))
    if count['partial tested']:
        lines.append('partial patterns (information only, never a verdict): '
                     '%d failing line(s) fit one within %s sigma; %d of them '
                     'have a line or a candidate transition where the part '
                     'left out would lie' % (
                         count['partial fits'], '%g' % FIT_SIGMA,
                         count['partial fits, evidence']))
    if count['held, not yet in the fit']:
        lines.append('%d held line(s) not yet in the fit (registry row newer '
                     'than the LOPT run; audited at kappa_held): %s' % (
                         count['held, not yet in the fit'], ', '.join(
                             r['wn_key'][:10] for r in rows if r['held']
                             and r['kappa_held'] != r['kappa_now'])))
    lines.append('')
    lines.append('char  era     n  median kappa_line  candidates (cg/head)'
                 '  class fails')
    groups = collections.defaultdict(list)
    for r in rows:
        if not r['held']:
            groups[(r['char'] or '-', r['era'])].append(r)
    for (char, era), rs in sorted(groups.items(),
                                  key=lambda kv: (kv[0][1], -len(kv[1]))):
        lines.append('%-5s %d %5d  %17.2f  %10d/%d  %11d' % (
            char, era, len(rs),
            median([float(r['kappa_line']) for r in rs]),
            sum(r['verdict'] == 'candidate-cg' for r in rs),
            sum(r['verdict'] == 'candidate-head' for r in rs),
            sum(r['verdict'] == 'class-fails' for r in rs)))
    # A level whose lines fail in different directions is more likely wrong
    # itself - its A or its energy - than each of its lines.
    bad = collections.defaultdict(list)
    for r in rows:
        if r['verdict'] in ('candidate-cg', 'candidate-head', 'class-fails',
                            'audit-fails'):
            for lid, n in ((r['low_id'], r['iden_low']),
                           (r['upp_id'], r['iden_upp'])):
                bad[(lid, n)].append('%s(%s)' % (r['wn_key'][:9],
                                                 r['kappa_line']))
    shared = [(k, v) for k, v in bad.items() if len(v) > 1]
    if shared:
        lines.append('')
        lines.append('levels with more than one line off its class '
                     '(kappa_line in brackets):')
        for (lid, n), v in sorted(shared, key=lambda kv: -len(kv[1])):
            lines.append('  %s (IDEN2 %s): %s' % (lid[-6:], n, ', '.join(v)))
    return '\n'.join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--set', default='iter_hfs', dest='set_dir',
                    help='the working set (default: %(default)s); it needs '
                         'the hfs correction on and a LOPT fit')
    ap.add_argument('--dchi2', type=float, default=DCHI2,
                    help='how much better (in chi^2) the head or cg must fit '
                         'than the class for a candidate; its square root '
                         'is also where a class fails and where an audit of '
                         'a held line fails (default: %(default)s, i.e. 3 '
                         'sigma)')
    ap.add_argument('--fit-sigma', type=float, default=FIT_SIGMA,
                    help='for a candidate only: how close the head or cg '
                         'must itself fit, in sigma; not used by the audit '
                         '(default: %(default)s)')
    ap.add_argument('--max-u-kappa', type=float, default=MAX_U_KAPPA,
                    help='the largest uncertainty of a line\'s own kappa '
                         'that is still tested (default: %(default)s)')
    args = ap.parse_args(argv)
    set_dir = os.path.join(HERE, args.set_dir)
    out = os.path.join(set_dir, OUT)
    output_files.require_writable([out], 'output file')
    rows, count = analyse(set_dir, args.dchi2, args.fit_sigma,
                          args.max_u_kappa)
    with open(out, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, COLUMNS, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print(summary(rows, count))
    print('\nwritten: %s' % os.path.relpath(out, HERE))
    return 0


if __name__ == '__main__':
    sys.exit(main())
