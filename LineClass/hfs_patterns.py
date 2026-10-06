"""What Sugar's resolved hfs components say about the A constants.

WHAT THIS IS FOR
================
`hfs_components.py` recovers 500 measured positions of individual hyperfine
components from Sugar's 1974 table and pairs each with the line it belongs to.
This module is what those measurements are for.  It answers the three
questions of `Work_on_hfs_plan.md` section 3.5, in order:

1. what feature of the pattern the tabulated wavelength of a flagged line is;
2. whether the semiempirical A constants from the level compositions reproduce
   the measured component spacings;
3. whether the assumption kappa(flag) = 1 - the single external assumption the
   corrected level list rests on - survives.

THE PATTERN, IN THE NOTATION USED HERE
======================================
The nucleus of 141Pr has spin I = 5/2.  A level of electronic angular momentum
J splits into components F = |I - J| ... I + J, and the component F lies at

    W(F) = (A/2) * K(F),    K(F) = F(F+1) - I(I+1) - J(J+1)

above the level's center of gravity, where A is the magnetic dipole hyperfine
constant of that level, in cm^-1.  (The electric quadrupole term is left out:
nothing in these data calls for it.)  A line between a lower level (J1, A1)
and an upper level (J2, A2) therefore breaks into components

    sigma(F1 -> F2) = sigma_cg + (A2/2) K(F2) - (A1/2) K(F1),   |F2 - F1| <= 1

whose relative strengths are the standard

    I(F1 -> F2) = (2F1+1)(2F2+1) * { J1  F1  I ; F2  J2  1 }^2

with the braces a Wigner 6-j symbol.  Those strengths concentrate in the six
"principal" components F1 -> F1 + (J2 - J1), which carry 80-90 % of the line
and form a ladder starting at F1 = I + J1.  The extreme member of that ladder,
F1 = I + J1 -> F2 = I + J2, is always the strongest single component, and it
sits at

    D = I * (A2 * J2 - A1 * J1)

from the center of gravity, because K(I+J) = 2 I J.  D is exactly the quantity
the plan calls D = S_upper - S_lower with S = I * A * J.

WHAT THE MEASUREMENTS SAY
=========================
Taking Sugar's tabulated wavelength to be that extreme component, and his k-th
listed companion to be the k-th rung below it, every listed position is
predicted with no free parameter at all from the A constants of
`hfs_level_shifts.csv`:

    d(k) = (k/2) * [ (2F1max + 1 - k) A1 - (2F2max + 1 - k) A2 ]

over 313 measured positions on 158 patterns whose two levels both have a
semiempirical A.  The result is an rms of 0.090 cm^-1, against a spread of
displacements that runs to 2.2 cm^-1 and a measurement precision of about
0.03 cm^-1 per component; the sign is right in 500 cases out of 500.  With one
free parameter - a single scale factor applied to every A constant in the
project - the rms falls to 0.084 cm^-1 at

    s = 0.958 +- 0.006

that is, the measured intervals are 4.2 % smaller than the semiempirical ones.
An alternative one-parameter model, a constant displacement of the tabulated
line towards the degraded side, does worse (0.085 cm^-1 at -0.030 cm^-1), and
when both are free the displacement collapses to +0.003 +- 0.005 cm^-1 while
the scale factor stays at 0.955.  So the tabulated wavenumber of a flagged
line is the strongest component of the pattern, with no measurable offset, and
not its center of gravity: the center of gravity is at sigma_tab - D.

Remeasured 2026-10-03 on the current assignments and the levels table
`A_hfs_levels.csv`, with only the *semiempirical* constants tested (sources
"Reader & Sugar 1965 (semiempirical)" and "composition..."; the 7 flag-interval
levels the first measurement included are measurements, not semiempirical values):
135 patterns, 271 positions, rms 0.072 -> 0.065 cm^-1 at

    s = 0.9638 +- 0.0079   (jackknife; the formal uncertainty is 0.0048)

which is `hfs_kappa.A_SCALE`.  It is applied to semiempirical constants only
(`hfs_kappa.scaled_A`), in the pipeline and in the fits of kappa alike.  The
first spacing has the right sign in 135 of 135 patterns.  Fitted on the *r
lines alone the scale is 0.981 +- 0.009, on the *v lines alone 0.958 +-
0.006; with scale and offset both free, 0.981 +- 0.011 and -0.016 +- 0.009
cm^-1 - 1.7 sigma, not adopted.

THE THREE ANSWERS
=================
1. The tabulated wavelength of a flagged 1974 line is the extreme principal
   component F = I+J1 -> I+J2, the strongest one.  This is decided, not
   assumed: were it the center of gravity the components would straddle it,
   and not one of the 500 does.
2. The semiempirical A constants reproduce the measured spacings to 4 % in scale
   and 0.084 cm^-1 in scatter.  The 4 % is a measurement about the radial
   parameters behind those constants and belongs in the error budget of every
   computed A; it is not absorbed silently here.
3. kappa(flag) = 1 is confirmed, and is no longer an assumption.  The
   displacement the pipeline already applies to a flagged line,
   `delta_cm-1` of `hfs_line_corrections.csv`, is identically
   I*(A2 J2 - A1 J1) with kappa = 1, and that is exactly the displacement of
   the strongest component from the center of gravity that these measurements
   confirm.

WHAT IS STILL OPEN
==================
The global fit that frees every level's A constant and solves for all of them
together (`fit_levels`) has a nearly singular direction: adding the same
constant to every A changes the predicted spacings only through the
J_low - J_upp difference, so the absolute scale of the A constants is weakly
determined even though their differences are not.  The fitted constants are
therefore reported with that caveat and are NOT written back into
`A_hfs_levels.csv`.  `hfs_A_fit.py` (Step 4, 2026-10-04) removes the
singular direction by holding the confirmed semiempirical constants at their
scaled values as priors, and writes the values it determines.

Run as ``python hfs_patterns.py``.  Nothing in the assignment pipeline is
touched, and no file is written unless ``--out`` is given.
"""
import argparse
import bisect
import collections
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))

COMPONENTS = 'hfs_components.csv'
LEVELS = 'hfs_level_shifts.csv'

# The nuclear spin of 141Pr, the only isotope in the spectrum.
I_SPIN = 2.5

# Provisional weight on one measured component position, cm^-1: Sugar's
# wavelengths carry three decimals in angstroms, which is 0.03 cm^-1 near
# 3000 A.  It sets the scale of the reported chi^2 and nothing else.
SIGMA = 0.03


# --- angular momentum -------------------------------------------------------
def _triangle(a, b, c):
    return math.sqrt(math.factorial(int(a + b - c))
                     * math.factorial(int(a - b + c))
                     * math.factorial(int(-a + b + c))
                     / math.factorial(int(a + b + c + 1)))


def sixj(a, b, c, d, e, f):
    """The Wigner 6-j symbol { a b c ; d e f }, by Racah's sum."""
    def ok(x, y, z):
        return abs(x - y) <= z <= x + y and (x + y + z) % 1 == 0
    if not (ok(a, b, c) and ok(a, e, f) and ok(d, b, f) and ok(d, e, c)):
        return 0.0
    pre = (_triangle(a, b, c) * _triangle(a, e, f)
           * _triangle(d, b, f) * _triangle(d, e, c))
    total = 0.0
    t = max(a + b + c, a + e + f, d + b + f, d + e + c)
    hi = min(a + b + d + e, b + c + e + f, a + c + d + f)
    while t <= hi:
        total += ((-1) ** int(t) * math.factorial(int(t + 1))
                  / (math.factorial(int(t - a - b - c))
                     * math.factorial(int(t - a - e - f))
                     * math.factorial(int(t - d - b - f))
                     * math.factorial(int(t - d - e - c))
                     * math.factorial(int(a + b + d + e - t))
                     * math.factorial(int(b + c + e + f - t))
                     * math.factorial(int(a + c + d + f - t))))
        t += 1
    return pre * total


def components(J1, A1, J2, A2, spin=I_SPIN):
    """Every hfs component of the line, as (position, strength, F1, F2).

    The position is measured from the center of gravity of the line and the
    strengths are relative, summing over the pattern to the line's strength.
    """
    def fs(J):
        n = int(round(2 * min(spin, J)))
        return [abs(spin - J) + k for k in range(n + 1)]

    out = []
    for f1 in fs(J1):
        k1 = f1 * (f1 + 1) - spin * (spin + 1) - J1 * (J1 + 1)
        for f2 in fs(J2):
            if abs(f1 - f2) > 1 or (f1 == 0 and f2 == 0):
                continue
            k2 = f2 * (f2 + 1) - spin * (spin + 1) - J2 * (J2 + 1)
            w = sixj(J1, f1, spin, f2, J2, 1)
            out.append((0.5 * A2 * k2 - 0.5 * A1 * k1,
                        (2 * f1 + 1) * (2 * f2 + 1) * w * w, f1, f2))
    return out


def head_displacement(J1, A1, J2, A2, spin=I_SPIN):
    """D, the strongest component's distance from the center of gravity."""
    return spin * (A2 * J2 - A1 * J1)


def rung(k, J1, A1, J2, A2, spin=I_SPIN):
    """Where the k-th principal component lies relative to the strongest one.

    k = 0 is the strongest component itself.  The ladder runs down in F on
    both levels together, F1 = I + J1 - k and F2 = I + J2 - k, and is defined
    while both stay above |I - J|.
    """
    f1, f2 = spin + J1, spin + J2
    return 0.5 * k * ((2 * f1 + 1 - k) * A1 - (2 * f2 + 1 - k) * A2)


def ladder_length(J1, J2, spin=I_SPIN):
    """How many rungs below the strongest component the ladder has."""
    return int(round(min(2 * min(spin, J1), 2 * min(spin, J2))))


#: what partial_pattern returns: the kappa of the measured part, its share
#: of the pattern's strength, and the part left out (None if nothing is) as
#: (side, position, share) - the side 'head' or 'tail', the position of its
#: center of gravity measured from the whole pattern's, cm^-1, and its share.
Partial = collections.namedtuple('Partial', 'kappa share missing')


def rung_of(J1, A1, J2, A2, spin=I_SPIN):
    """`[(position, strength, k)]`: every component of the pattern with the
    rung k of the principal ladder it lies nearest (the first of equals),
    positions measured from the center of gravity.  A component off the
    ladder thus joins the part of the pattern it sits in."""
    D = head_displacement(J1, A1, J2, A2, spin)
    rungs = [D + rung(k, J1, A1, J2, A2, spin)
             for k in range(ladder_length(J1, J2, spin) + 1)]
    return [(p, s, min(range(len(rungs)), key=lambda k: abs(p - rungs[k])))
            for p, s, _, _ in components(J1, A1, J2, A2, spin)]


def partial_pattern(J1, A1, J2, A2, first, last=None, spin=I_SPIN):
    """Where a line is measured that saw only the rungs `first` ... `last` of
    its pattern (`last` None: to the end of the ladder), each rung with the
    components nearest it; see `rung_of`.

    Rung 0 is the strongest component, the head; a part left out lies on one
    side of what is kept, or on both if `first` > 0 and `last` < the end -
    `missing` then describes the head side.  kappa is the measured part's
    center of gravity over D, the head's: 1 at the head, 0 at the whole
    pattern's center of gravity.  Raises ValueError for a range outside the
    ladder, for the whole ladder (that is the class cg) and for D = 0.
    """
    n = ladder_length(J1, J2, spin)
    last = n if last is None else last
    if not 0 <= first <= last <= n:
        raise ValueError('rungs %d-%d are outside the ladder 0-%d of this '
                         'pattern' % (first, last, n))
    if first == 0 and last == n:
        raise ValueError('rungs 0-%d are the whole pattern: that is the '
                         'class cg' % n)
    D = head_displacement(J1, A1, J2, A2, spin)
    if D == 0.0:
        raise ValueError('the pattern has D = 0: a part of it has no kappa')
    comps = [c for c in rung_of(J1, A1, J2, A2, spin) if c[1] > 1e-12]
    total = sum(s for _, s, _ in comps)
    if not any(first <= k <= last for _, _, k in comps):
        raise ValueError('rungs %d-%d carry no strength in this pattern'
                         % (first, last))
    if all(first <= k <= last for _, _, k in comps):
        raise ValueError('rungs %d-%d carry the whole strength of this '
                         'pattern: that is the class cg' % (first, last))
    kept = [(p, s) for p, s, k in comps if first <= k <= last]
    head = [(p, s) for p, s, k in comps if k < first]
    tail = [(p, s) for p, s, k in comps if k > last]

    def centre(part):
        w = sum(s for _, s in part)
        return sum(p * s for p, s in part) / w, w

    p, w = centre(kept)
    out = head or tail
    missing = None
    if out:
        pm, wm = centre(out)
        missing = ('head' if head else 'tail', pm, wm / total)
    return Partial(p / D, w / total, missing)


# --- the measurements -------------------------------------------------------
Pattern = collections.namedtuple(
    'Pattern', 'wn char low_id upp_id J1 J2 A1 A2 u1 u2 src1 src2 dwn')


def read_patterns(components_csv=None, levels_csv=None):
    """Group `hfs_components.csv` into patterns and attach the two levels.

    Only the patterns of a line with exactly one accepted classification can
    be used: a blend does not belong to one pair of levels.  Returns
    (patterns, skipped) with `skipped` counting why each was left out.
    """
    path = components_csv or os.path.join(HERE, COMPONENTS)
    levels_path = levels_csv or os.path.join(HERE, LEVELS)
    with open(levels_path, encoding='utf-8', newline='') as fh:
        levels = {r['level_id']: r for r in csv.DictReader(fh)}
    with open(path, encoding='utf-8', newline='') as fh:
        rows = list(csv.DictReader(fh))
    return patterns_from(rows, levels, os.path.basename(levels_path))


def patterns_from(rows, levels, levels_name='the level table'):
    """`read_patterns` on rows already read: `rows` as `hfs_components.csv`
    gives them, `levels` as `{level_id: row}` with the columns J, A_cm-1,
    u_A and A_source."""
    grouped = collections.OrderedDict()
    for r in rows:
        grouped.setdefault(r['wn_line'], []).append(r)

    out, skipped = [], collections.Counter()
    for wn, group in grouped.items():
        first = group[0]
        if first['n_class'] != '1':
            skipped['line unclassified or a blend'] += 1
            continue
        low, upp = levels.get(first['low_id']), levels.get(first['upp_id'])
        if low is None or upp is None:
            skipped['level absent from %s' % levels_name] += 1
            continue
        J1, J2 = float(low['J']), float(upp['J'])
        dwn = sorted((float(r['dwn']) for r in group), key=abs)
        if len(dwn) > ladder_length(J1, J2):
            skipped['more components than the ladder has rungs'] += 1
            continue
        out.append(Pattern(float(wn), first['char'], first['low_id'],
                           first['upp_id'], J1, J2,
                           float(low['A_cm-1']), float(upp['A_cm-1']),
                           float(low['u_A']), float(upp['u_A']),
                           low['A_source'], upp['A_source'], dwn))
    return out, skipped


def with_both_A(patterns):
    """The patterns whose two levels both have an A constant to test."""
    return [p for p in patterns
            if 'not determined' not in (p.src1, p.src2)]


def residuals(patterns):
    """(predicted, observed, sign) for every measured component position."""
    pred, obs, sign = [], [], []
    for p in patterns:
        s = 1.0 if p.char == '*v' else -1.0
        for k, d in enumerate(p.dwn, 1):
            pred.append(rung(k, p.J1, p.A1, p.J2, p.A2))
            obs.append(d)
            sign.append(s)
    return pred, obs, sign


def _rms(v):
    return math.sqrt(sum(x * x for x in v) / len(v)) if v else float('nan')


def scale_tests(patterns):
    """Fit the three one- and two-parameter models of the docstring."""
    pred, obs, sign = residuals(patterns)
    n = len(obs)
    pp = sum(x * x for x in pred)
    s = sum(a * b for a, b in zip(pred, obs)) / pp
    r_s = [b - s * a for a, b in zip(pred, obs)]
    u_s = _rms(r_s) / math.sqrt(pp)
    c = sum(g * (b - a) for a, b, g in zip(pred, obs, sign)) / n
    r_c = [b - a - c * g for a, b, g in zip(pred, obs, sign)]
    u_c = _rms(r_c) / math.sqrt(n)
    # both free: two normal equations in (s, c)
    sc = sum(a * g for a, g in zip(pred, sign))
    det = pp * n - sc * sc
    rhs1 = sum(a * b for a, b in zip(pred, obs))
    rhs2 = sum(g * b for g, b in zip(sign, obs))
    s2 = (rhs1 * n - sc * rhs2) / det
    c2 = (pp * rhs2 - sc * rhs1) / det
    r_2 = [b - s2 * a - c2 * g for a, b, g in zip(pred, obs, sign)]
    return {'n': n,
            'rms_none': _rms([b - a for a, b in zip(pred, obs)]),
            'scale': s, 'u_scale': u_s, 'rms_scale': _rms(r_s),
            'offset': c, 'u_offset': u_c, 'rms_offset': _rms(r_c),
            'both': (s2, c2), 'rms_both': _rms(r_2)}


def fit_levels(patterns, sigma=SIGMA):
    """Solve for every level's A constant from the measured spacings.

    Returns (ids, A, u_A, n_patterns, diagnostics).  See the docstring for why
    the absolute scale of the solution is weakly determined; `diagnostics`
    carries the rank and the chi^2 so that the caller can say so.
    """
    import numpy as np

    ids = sorted({p.low_id for p in patterns} | {p.upp_id for p in patterns})
    index = {x: i for i, x in enumerate(ids)}
    counts = collections.Counter()
    rows, y = [], []
    for p in patterns:
        counts[p.low_id] += 1
        counts[p.upp_id] += 1
        f1, f2 = I_SPIN + p.J1, I_SPIN + p.J2
        for k, d in enumerate(p.dwn, 1):
            row = [0.0] * len(ids)
            row[index[p.low_id]] += 0.5 * k * (2 * f1 + 1 - k)
            row[index[p.upp_id]] -= 0.5 * k * (2 * f2 + 1 - k)
            rows.append(row)
            y.append(d)
    m = np.array(rows) / sigma
    rhs = np.array(y) / sigma
    sol, _, rank, _ = np.linalg.lstsq(m, rhs, rcond=None)
    resid = (m @ sol - rhs) * sigma
    cov = np.linalg.pinv(m.T @ m)
    dof = len(y) - rank
    diag = {'equations': len(y), 'unknowns': len(ids), 'rank': int(rank),
            'rms': float(np.sqrt(np.mean(resid ** 2))),
            'chi2_per_dof': float(np.sum((m @ sol - rhs) ** 2) / dof)}
    return ids, sol, np.sqrt(np.diag(cov)), counts, diag


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--components', default=None)
    ap.add_argument('--levels', default=None)
    ap.add_argument('--out', default=None,
                    help='write the fitted A constants to this csv')
    args = ap.parse_args(argv)

    patterns, skipped = read_patterns(args.components, args.levels)
    print('patterns usable                  %4d' % len(patterns))
    for k, v in sorted(skipped.items()):
        print('   set aside, %-38s %3d' % (k + ':', v))
    tested = with_both_A(patterns)
    print('patterns whose two levels both have an A  %4d' % len(tested))

    t = scale_tests(tested)
    print('\nthe tabulated line as the strongest component, %d measured '
          'positions' % t['n'])
    print('   no free parameter                      rms %.4f cm^-1'
          % t['rms_none'])
    print('   one scale factor on every A: %.4f +- %.4f   rms %.4f'
          % (t['scale'], t['u_scale'], t['rms_scale']))
    print('   one offset towards the degraded side: %+.4f +- %.4f  rms %.4f'
          % (t['offset'], t['u_offset'], t['rms_offset']))
    print('   both free: scale %.4f, offset %+.4f          rms %.4f'
          % (t['both'][0], t['both'][1], t['rms_both']))

    ids, sol, err, counts, diag = fit_levels(patterns)
    print('\nglobal fit: %(equations)d equations, %(unknowns)d level constants,'
          ' rank %(rank)d' % diag)
    print('   rms %(rms).4f cm^-1, chi2/dof %(chi2_per_dof).2f' % diag)
    with open(args.levels or os.path.join(HERE, LEVELS),
              encoding='utf-8', newline='') as fh:
        levels = {r['level_id']: r for r in csv.DictReader(fh)}
    known = [(x, sol[i], err[i], levels[x]) for i, x in enumerate(ids)
             if levels[x]['A_source'] != 'not determined' and err[i] < 0.02]
    d = [a - float(r['A_cm-1']) for _, a, _, r in known]
    print('   against the semiempirical constants, %d levels with u < 0.02:'
          % len(known))
    print('      mean difference %+.4f, scatter about it %.4f cm^-1'
          % (sum(d) / len(d), _rms([x - sum(d) / len(d) for x in d])))
    print('   levels whose A these measurements give for the first time: %d'
          % sum(1 for i, x in enumerate(ids)
                if levels[x]['A_source'] == 'not determined' and err[i] < 0.02))

    if args.out:
        with open(args.out, 'w', encoding='utf-8', newline='\n') as fh:
            w = csv.writer(fh, lineterminator='\n')
            w.writerow(['level_id', 'J', 'n_patterns', 'A_fit', 'u_fit',
                        'A_calc', 'u_calc', 'A_source'])
            for i, x in enumerate(ids):
                r = levels[x]
                w.writerow([x, r['J'], counts[x], '%+.4f' % sol[i],
                            '%.4f' % err[i], r['A_cm-1'], r['u_A'],
                            r['A_source']])
        print('\nwrote %s (%d levels)' % (args.out, len(ids)))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
