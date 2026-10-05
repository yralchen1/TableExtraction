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
--fit-sigma (2).  A row already in the registry is audited: it should still
fit its own kappa.

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
import statistics
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
MAX_U_KAPPA = 0.3

COLUMNS = ('verdict', 'wn_key', 'char', 'era', 'low_id', 'upp_id',
           'iden_low', 'iden_upp', 'obs_intens', 'D', 'kappa_class',
           'kappa_now', 'held', 'kappa_line', 'u_kappa_line', 'leverage',
           'u_rest', 'u_model', 'u_hfs_shift', 'r_loo', 'z_class', 'z_head',
           'z_cg', 'dchi2', 'registry_unc', 'registry_reason',
           'nearest_dwn', 'nearest_intens', 'rung1_dwn', 'line_at_rung1')

#: how close a line of the list must lie to the predicted rung 1 to be
#: named in `line_at_rung1`, cm^-1.
RUNG_WINDOW = 0.1


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
    character without one, Sugar's stated value if neither is there."""
    era = hfs_kappa.era_of(wn)
    u = table.get((char, era)) or table.get(('other', era))
    if u is None:
        return max(hfs_kappa.stated_uncertainty(char, era, wn),
                   hfs_kappa.FLOOR)
    return max(hfs_kappa.two_term(u.a, u.b, wn), hfs_kappa.FLOOR)


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


def hypotheses(r_loo, kappa_now, kappa_class, D, sigma):
    """`{'class', 'head', 'cg'}` -> z, the residual without the line under
    each kappa in units of sigma."""
    def z(k):
        return (r_loo + (kappa_now - k) * D) / sigma
    return {'class': z(kappa_class), 'head': z(1.0), 'cg': z(0.0)}


def verdict(z, held, dchi2_min=DCHI2, fit_sigma=FIT_SIGMA, kappa_class=None):
    """`(verdict, dchi2)` of a testable line.

    `held` is the registry's class name for the line ('' if none).  A held
    line is audited against its own kappa: 'audit-ok' or 'audit-fails'.
    Otherwise 'candidate-head' or 'candidate-cg' when the class loses to the
    better alternative by at least dchi2_min and that one fits within
    fit_sigma; 'class-fails' when the class is rejected and neither
    alternative fits; 'fits-class' else.  A class whose kappa is itself 1
    (a flagged line) has no head alternative."""
    alts = ['cg'] if kappa_class is not None and kappa_class >= 1.0 \
        else ['head', 'cg']
    best = min(alts, key=lambda k: abs(z[k]))
    dchi2 = z['class'] ** 2 - z[best] ** 2
    if held:
        own = 'head' if held in ('head', 'flag') else (
            'cg' if held == 'cg' else 'class')
        return ('audit-ok' if abs(z[own]) <= fit_sigma else 'audit-fails',
                dchi2)
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


def line_near(keys, wn, window=RUNG_WINDOW):
    """The line of the sorted list `keys` nearest `wn`, if within
    `window`."""
    i = bisect.bisect_left(keys, wn)
    near = [keys[j] for j in (i - 1, i) if 0 <= j < len(keys)]
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
        z = hypotheses(r_loo, kappa_now, kappa_class, D, sigma)
        v, dchi2 = verdict(z, held, dchi2_min, fit_sigma, kappa_class)
        count[v] += 1
        near = nearest_other(keys, intens, pos[key])
        reg_u = registry.lookup(key)
        r1 = rung1(model, low, upp)
        at_r1 = None if r1 is None else line_near(keys, key + r1)
        rows.append({
            'verdict': v, 'wn_key': r['wn_key'], 'char': char,
            'era': hfs_kappa.era_of(key), 'low_id': low, 'upp_id': upp,
            'iden_low': iden.get(low, ''), 'iden_upp': iden.get(upp, ''),
            'obs_intens': '%.0f' % float(r['obs_intens']),
            'D': '%+.4f' % D, 'kappa_class': '%.3f' % kappa_class,
            'kappa_now': '%.3f' % kappa_now, 'held': held,
            'kappa_line': '%.2f' % (kappa_now + r_loo / D),
            'u_kappa_line': '%.2f' % (sigma / abs(D)),
            'leverage': '%.3f' % h, 'u_rest': '%.4f' % u_rest,
            'u_model': '%.4f' % u_mod, 'u_hfs_shift': '%.4f' % u_shift,
            'r_loo': '%+.4f' % r_loo,
            'z_class': '%+.1f' % z['class'], 'z_head': '%+.1f' % z['head'],
            'z_cg': '%+.1f' % z['cg'], 'dchi2': '%.1f' % dchi2,
            'registry_unc': '' if reg_u is None else '%g' % reg_u,
            'registry_reason': '' if reg_u is None else reasons.get(
                _registry_name(reasons, key), ''),
            'nearest_dwn': '' if near[0] is None else '%+.3f' % near[0],
            'nearest_intens': '' if near[1] is None else '%.0f' % near[1],
            'rung1_dwn': '' if r1 is None else '%+.3f' % r1,
            'line_at_rung1': '' if at_r1 is None else '%.4f (I %.0f)' % (
                at_r1, intens_of[at_r1]),
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
            statistics.median(float(r['kappa_line']) for r in rs),
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
                    help='how much better the alternative must fit '
                         '(default: %(default)s)')
    ap.add_argument('--fit-sigma', type=float, default=FIT_SIGMA,
                    help='how close the alternative must fit, in sigma '
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
