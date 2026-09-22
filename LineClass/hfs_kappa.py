"""The measurement convention of Sugar's line lists, fitted on honest weights.

WHAT THIS IS FOR
================
This is Step 3 of `Work_on_hfs_plan.md`, together with the part of Step 1 that
Step 3 cannot move without.  It answers one question:

    when Sugar tabulated the wavelength of a line whose hyperfine pattern he
    could NOT resolve, did he measure the pattern's center of gravity, or did
    he measure its strongest component, as he demonstrably did for the lines
    he flagged?

THE MODEL, AND WHY ONLY DIFFERENCES ARE MEASURABLE
==================================================
141Pr has nuclear spin I = 5/2, so every level of electronic angular momentum
J is a group of sublevels and every line is a pattern of components.  Write

    A       the magnetic dipole hyperfine constant of a level, in cm^-1
    S       = I * A * J, the displacement of the outermost sublevel (F = I+J)
            of that level from the level's own center of gravity
    D       = S_upper - S_lower, the displacement of the pattern's extreme
            component from the center of gravity of the LINE

and model the tabulated wavenumber of a line as

    wn_tabulated = wn_center_of_gravity + kappa * D

with one convention factor `kappa` per class of line.  kappa = 1 means the
measurement was made at the extreme component, kappa = 0 at the center of
gravity.  The three classes are the ones Sugar's own characters define:

    flag    the lines he marked `*r` or `*v`, whose pattern he saw
    c       the lines he marked `c` for complex
    plain   everything else

**The line list alone cannot measure the absolute scale of kappa.**  This is
the degeneracy theorem of the plan's section 2.1, and it is exact: adding the
same constant c to the kappa of every class, and c * S_k to the fitted energy
of every level k, changes no predicted wavenumber whatever, because a line's
prediction contains the two only in the combination
`(E_u + c S_u) - (E_l + c S_l) - (kappa + c)(S_u - S_l)`.  A fit of all three
kappas at once is therefore rank-deficient by exactly one, and what the data
measure are the two differences.  One class must be anchored from outside, and
Step 2 is what anchors it: Sugar's own resolved components put his tabulated
flagged wavelength at the extreme component of the pattern, with no measurable
offset, so `kappa(flag) = 1` is now a measurement rather than an assumption.

Step 2 measured one thing more that belongs here.  Fitted against the observed
component spacings, the calculated A constants are 4.2 +- 0.6 per cent too
large.  Since `D` is computed from those constants, every kappa fitted against
`D` absorbs that error, and the honest anchor is not `kappa(flag) = 1` but
`kappa(flag) = 0.958 +- 0.006` with `D` as computed.  Both anchors are
reported below; the difference between them is smaller than the uncertainty of
the result, which is worth knowing.

THE WEIGHTS
===========
The weights are the reason this fit is being done again.  The `unc_own` column
of `Pr3_lines.xlsx` inflates the uncertainty of a line that departs from its
Ritz value, which is right for Set 1 and fatal here: a trend fitted through
those weights is fitted through a function of the residuals it is trying to
explain.  The plan's section 2.6 prescribes what to do instead, and this
module does it:

1. start every line at the uncertainty Sugar himself states for its character
   and era - 0.0030 angstrom for a plain 1974 line, 0.0040 angstrom for a
   plain 1969 line, 0.0070 angstrom, his stated average deviation from Ritz,
   for everything else, converted to cm^-1 through `u_wn = u_lambda * wn^2`;
2. fit;
3. set aside the lines that depart from the fit by more than four of their own
   sigma as probably misassigned - they are candidates for review, not
   deletions, and nothing is written to `line_decisions.csv` here;
4. adopt the rms of the surviving residuals, per character per era, as the
   uncertainty of that class, and repeat from 2 until the adopted values stop
   moving.

Step 1 then has its deliverable - `hfs_uncertainties.csv`, a number per
character per era with the count behind it, reproducible by a reader from the
published line table - and Step 3 has its weights.

WHAT THIS MODULE DOES NOT DO
============================
It does not touch the pipeline, `A_hfs_levels.csv` or `hfs_line_corrections.csv`,
and it writes no ledger row.  The line set it fits is the one
`line_classifications.csv` currently accepts, singly assigned; a line assigned
to more than one transition is a blend whose measured wavenumber belongs to no
single pair of levels, and is left out.
"""
import argparse
import collections
import contextlib
import csv
import io
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

LINES = 'line_classifications.csv'
A_LEVELS = 'A_hfs_levels.csv'
OUT_UNCERTAINTIES = 'hfs_uncertainties.csv'
OUT_MISASSIGNED = 'hfs_misassigned_candidates.csv'

#: nuclear spin of 141Pr, the only stable isotope.
I_SPIN = 2.5

#: Sugar measured the region below this wavenumber in 1974 and the region
#: above it in 1969; the plan calls it the era boundary.
ERA_SPLIT = 47500.0

#: what Sugar states, in angstrom.  A plain line of 1974 is "~0.003", a line
#: of 1969 "+-0.004"; for everything else he states only an average deviation
#: from Ritz values of "+-0.007", and that is what the other characters start
#: from.  These are starting values only: the adopted uncertainties come out
#: of the residuals.
STATED_PLAIN = {1974: 0.0030, 1969: 0.0040}
STATED_OTHER = 0.0070

#: the scale error of the calculated A constants measured in Step 2, and its
#: uncertainty.  `hfs_patterns.scale_tests` is where these come from.
A_SCALE = 0.9582
U_A_SCALE = 0.0058

#: a line this far from the fit in units of its own uncertainty is set aside
#: as probably misassigned.
OUTLIER_SIGMA = 4.0

#: no class is given an uncertainty below this, the precision floor of the
#: project's own uncertainty model (README, `UNC_FLOOR`): a measurement read
#: off a photographic plate is not more precise than this however small the
#: residuals of a handful of lines happen to be.
FLOOR = 0.0055

#: an uncertainty class needs this many lines before the residuals of that
#: character and era are allowed to speak for themselves; below it the
#: characters are pooled into one "other" class per era.
MIN_CLASS = 20

Line = collections.namedtuple(
    'Line', 'wn char era cls ucls low upp D dJ')


def era_of(wn):
    """1974 below the boundary, 1969 above it."""
    return 1974 if wn < ERA_SPLIT else 1969


def kappa_class(char):
    """The three classes of the convention model."""
    if char in ('*r', '*v'):
        return 'flag'
    if char == 'c':
        return 'c'
    return 'plain'


def stated_uncertainty(char, era, wn):
    """Sugar's own stated uncertainty for this line, in cm^-1.

    `u_lambda` in angstrom becomes `u_wn = u_lambda * 1e-8 * wn^2` in cm^-1,
    because wn = 1e8 / lambda.
    """
    ang = STATED_PLAIN[era] if char == '' else STATED_OTHER
    return ang * 1e-8 * wn * wn


def read_A_constants(path=None):
    """`level_id -> (J, A, u_A, source)` from `A_hfs_levels.csv`.

    A level whose constant is undetermined, and a level whose calculated
    constant contradicts its own flags (section 3.2), are both given A = 0
    here - which is what `hfs_line_corrections.csv` does, so that the two
    agree on every line they share.  A level absent from the file has no
    constant at all and is likewise given 0.
    """
    path = path or os.path.join(HERE, A_LEVELS)
    out = {}
    with open(path, encoding='utf-8', newline='') as fh:
        for row in csv.DictReader(fh):
            src = row['source']
            A = float(row['A_cm-1'])
            if src == 'not determined' or 'CONFLICT' in src:
                A = 0.0
            out[row['level_id']] = (float(row['J']), A, float(row['u_A']), src)
    return out


def read_levels():
    """`level_id -> J`, through the pipeline's own reader.

    `classify_lines.read_energy_levels` is read-only: it opens the published
    level workbook, adds the levels found since it was published, drops the
    ones given up and applies the revised energies, which is exactly the level
    set `line_classifications.csv` was produced with.  Its progress report is
    swallowed here.
    """
    import classify_lines
    with contextlib.redirect_stdout(io.StringIO()):
        levels, _ = classify_lines.read_energy_levels()
    return {lid: lev.J_val for lid, lev in levels.items()}


def read_lines(path=None, constants=None, J_of=None):
    """The accepted, singly assigned lines, with their displacement D.

    A line accepted for more than one transition is a blend: its measured
    wavenumber is the position of a feature made of several patterns and
    belongs to no single pair of levels, so it is not used here.
    """
    path = path or os.path.join(HERE, LINES)
    constants = read_A_constants() if constants is None else constants
    J_of = read_levels() if J_of is None else J_of

    def S(level_id):
        rec = constants.get(level_id)
        return 0.0 if rec is None else I_SPIN * rec[1] * rec[0]

    lines = []
    with open(path, encoding='utf-8', newline='') as fh:
        for row in csv.DictReader(fh):
            if row['accepted'] != '1.0' or row['n_accepted'] != '1':
                continue
            low, upp = row['low_id'], row['upp_id']
            if not low or not upp:
                continue
            wn = float(row['wn_obs'])
            char = row['char']
            era = era_of(wn)
            jl, ju = J_of.get(low), J_of.get(upp)
            dJ = None if jl is None or ju is None else ju - jl
            lines.append(Line(wn=wn, char=char, era=era,
                              cls=kappa_class(char), ucls=(char, era),
                              low=low, upp=upp, D=S(upp) - S(low), dJ=dJ))
    return lines


def uncertainty_classes(lines, min_class=MIN_CLASS):
    """Group the lines by character and era, pooling the rare characters.

    Returns `{line index: class key}`.  A character with fewer than
    `min_class` lines in an era cannot measure its own rms, and joins the
    pooled class `('other', era)` of that era.
    """
    counts = collections.Counter(ln.ucls for ln in lines)
    keys = {}
    for i, ln in enumerate(lines):
        keys[i] = ln.ucls if counts[ln.ucls] >= min_class else ('other', ln.era)
    return keys


def solve(lines, class_of, anchor, anchor_value, scale=1.0, sigma=None,
          use=None):
    """One weighted least-squares solve for the level values and the kappas.

    `class_of(line) -> key` says which convention class a line belongs to;
    the class named by `anchor` is held at `anchor_value` and the rest are
    free, which is what makes the system full rank (see the degeneracy
    theorem in the module docstring).  `scale` multiplies every displacement
    D, and is how Step 2's measured scale error of the A constants is carried
    in.  `use` is an optional boolean mask of the lines to fit.

    Returns a dict with the fitted kappas, their covariance, the residual of
    every line (including the ones masked out), and the fit diagnostics.
    """
    sigma = [stated_uncertainty(ln.char, ln.era, ln.wn) for ln in lines] \
        if sigma is None else sigma
    use = [True] * len(lines) if use is None else use

    keys = [class_of(ln) for ln in lines]
    free = sorted({k for k, u in zip(keys, use) if u and k != anchor},
                  key=str)
    levels = sorted({ln.low for ln, u in zip(lines, use) if u}
                    | {ln.upp for ln, u in zip(lines, use) if u})
    index = {lid: i for i, lid in enumerate(levels)}
    kcol = {k: len(levels) + i for i, k in enumerate(free)}

    rows = [i for i, u in enumerate(use) if u]
    M = np.zeros((len(rows), len(levels) + len(free)))
    b = np.zeros(len(rows))
    for r, i in enumerate(rows):
        ln = lines[i]
        w = 1.0 / sigma[i]
        M[r, index[ln.upp]] += w
        M[r, index[ln.low]] -= w
        rhs = ln.wn
        if keys[i] == anchor:
            rhs -= anchor_value * scale * ln.D
        else:
            M[r, kcol[keys[i]]] += w * scale * ln.D
        b[r] = w * rhs

    x, _, rank, _ = np.linalg.lstsq(M, b, rcond=None)
    cov = np.linalg.pinv(M.T @ M)

    # The leverage of a line: how much of its own residual the fit has
    # already absorbed into the level values it shares.  A line that is the
    # only one determining a level has h = 1 and a residual of zero however
    # badly it was measured, so the residuals of a class cannot be read as
    # its uncertainty without this.  Each row has at most three non-zero
    # entries, so the quadratic form is cheap.
    leverage = [0.0] * len(lines)
    for r, i in enumerate(rows):
        cols = np.nonzero(M[r])[0]
        v = M[r, cols]
        leverage[i] = float(v @ cov[np.ix_(cols, cols)] @ v)

    energy = {lid: x[index[lid]] for lid in levels}
    kappa = {anchor: anchor_value}
    u_kappa = {anchor: 0.0}
    for k in free:
        kappa[k] = x[kcol[k]]
        u_kappa[k] = math.sqrt(max(cov[kcol[k], kcol[k]], 0.0))

    residual = []
    for i, ln in enumerate(lines):
        e_l = energy.get(ln.low)
        e_u = energy.get(ln.upp)
        if e_l is None or e_u is None:
            residual.append(None)
            continue
        k = kappa.get(keys[i])
        if k is None:                      # a class present only among the
            k = anchor_value               # lines this solve left out
        residual.append(ln.wn - (e_u - e_l + k * scale * ln.D))

    chi2 = sum((residual[i] / sigma[i]) ** 2 for i in rows)
    dof = len(rows) - rank
    return {
        'kappa': kappa, 'u_kappa': u_kappa, 'free': free, 'anchor': anchor,
        'cov': cov, 'kcol': kcol, 'energy': energy, 'residual': residual,
        'sigma': sigma, 'used': use, 'rows': rows,
        'chi2': chi2, 'dof': dof, 'rank': rank, 'leverage': leverage,
        'unknowns': len(levels) + len(free), 'n': len(rows),
        'rms': math.sqrt(sum(residual[i] ** 2 for i in rows) / len(rows)),
    }


def difference(fit, a, b):
    """kappa(a) - kappa(b) and its uncertainty, correlations included."""
    ka = fit['kappa'][a] - fit['kappa'][b]
    cov, kcol = fit['cov'], fit['kcol']
    va = cov[kcol[a], kcol[a]] if a in kcol else 0.0
    vb = cov[kcol[b], kcol[b]] if b in kcol else 0.0
    cab = cov[kcol[a], kcol[b]] if a in kcol and b in kcol else 0.0
    return ka, math.sqrt(max(va + vb - 2 * cab, 0.0))


def adopt_uncertainties(lines, scale=1.0, outlier=OUTLIER_SIGMA,
                        max_rounds=10, tol=0.01):
    """The Step 1 loop: stated uncertainties in, adopted uncertainties out.

    Returns `(sigma, keep, table, fit)`: the adopted uncertainty of every
    line, the mask of the lines that survived the outlier filter, the adopted
    value per uncertainty class with the count behind it, and the last fit.

    The rms of the residuals of a fit understates the uncertainty, because the
    fitted level values have already followed the residuals; the factor
    `sqrt(n / (n - rank))` puts that back, and is applied to every class
    alike since the level values are shared.
    """
    ucls = uncertainty_classes(lines)
    sigma = [stated_uncertainty(ln.char, ln.era, ln.wn) for ln in lines]
    keep = [True] * len(lines)
    stated = dict(sigma=list(sigma))
    table, fit = {}, None
    for _ in range(max_rounds):
        fit = solve(lines, lambda ln: ln.cls, 'flag', 1.0, scale=scale,
                    sigma=sigma, use=keep)
        res = fit['residual']
        lev = fit['leverage']
        keep = [res[i] is not None and abs(res[i]) <= outlier * sigma[i]
                for i in range(len(lines))]
        groups = collections.defaultdict(list)
        for i, ln in enumerate(lines):
            if keep[i]:
                groups[ucls[i]].append((res[i], lev[i]))
        new = {}
        for key, values in groups.items():
            # E[r^2] = sigma^2 (1 - h), so the sum of (1 - h) is the number
            # of residuals the class really has left to measure itself with.
            free_n = sum(1.0 - h for _, h in values)
            rms = math.sqrt(sum(r * r for r, _ in values) / max(free_n, 1.0))
            new[key] = (max(rms, FLOOR), len(values))
        moved = max(abs(new[k][0] - table.get(k, (0.0,))[0])
                    / max(new[k][0], 1e-12) for k in new)
        table = new
        sigma = [table[ucls[i]][0] if ucls[i] in table else sigma[i]
                 for i in range(len(lines))]
        if moved < tol:
            break
    # the stated value each class started from, for the report
    for key in table:
        lo = [stated['sigma'][i] for i in range(len(lines))
              if ucls[i] == key]
        table[key] = (table[key][0], table[key][1],
                      sum(lo) / len(lo) if lo else float('nan'))
    return sigma, keep, table, fit


def write_uncertainties(table, path=None):
    """The Step 1 deliverable, one row per character per era."""
    path = path or os.path.join(HERE, OUT_UNCERTAINTIES)
    rows = sorted(table.items(), key=lambda kv: (kv[0][1], kv[0][0]))
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        w = csv.writer(fh, lineterminator='\n')
        w.writerow(['char', 'era', 'n_lines', 'u_stated_mean_cm-1',
                    'u_adopted_cm-1'])
        for (char, era), (adopted, n, started) in rows:
            w.writerow([char, era, n, '%.4f' % started, '%.4f' % adopted])
    return path


def write_misassigned(lines, fit, sigma, keep, path=None):
    """The lines the trend fit rejects, as candidates for review by hand."""
    path = path or os.path.join(HERE, OUT_MISASSIGNED)
    res = fit['residual']
    bad = [i for i in range(len(lines))
           if not keep[i] and res[i] is not None]
    bad.sort(key=lambda i: -abs(res[i] / sigma[i]))
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        w = csv.writer(fh, lineterminator='\n')
        w.writerow(['wn_obs', 'char', 'class', 'low_id', 'upp_id',
                    'D_cm-1', 'residual_cm-1', 'u_adopted_cm-1', 'n_sigma'])
        for i in bad:
            ln = lines[i]
            w.writerow(['%.4f' % ln.wn, ln.char, ln.cls, ln.low, ln.upp,
                        '%+.4f' % ln.D, '%+.4f' % res[i], '%.4f' % sigma[i],
                        '%+.1f' % (res[i] / sigma[i])])
    return path, len(bad)


def era_test(lines, sigma, keep, scale=1.0):
    """Section 3.3, first test: does kappa(plain) drift between the eras?

    A pattern resolved at 4000 angstrom need not be resolved at 2000, so the
    effective convention factor of an unflagged line could drift with
    wavelength.  The two eras are also two different measurement campaigns,
    which makes this the natural place to split.
    """
    def class_of(ln):
        return (ln.cls, ln.era) if ln.cls == 'plain' else ln.cls
    return solve(lines, class_of, 'flag', 1.0, scale=scale, sigma=sigma,
                 use=keep)


def branch_test(lines, sigma, keep, scale=1.0):
    """Section 3.3, second test: does kappa(flag) depend on the branch?

    Step 2 has already answered the physics - the extreme component is the
    strongest one for every J pair the line list contains - but the residuals
    can be asked independently.  The flag class is split by
    dJ = J_upper - J_lower and the dJ = +1 branch is held at 1.
    """
    def class_of(ln):
        if ln.cls != 'flag' or ln.dJ is None:
            return ln.cls
        return ('flag', int(round(ln.dJ)))
    return solve(lines, class_of, ('flag', 1), 1.0, scale=scale, sigma=sigma,
                 use=keep)


def drift_scan(lines, sigma, keep, scale=1.0, nbins=6):
    """How kappa(plain) moves across the spectrum, in bins of wavenumber.

    The era test only asks whether the two campaigns differ.  If the cause is
    the one section 3.3 names - a pattern that is resolved at 4000 angstrom
    and not at 2000 - then kappa should drift smoothly with wavenumber rather
    than step at the boundary, and the bins say which it does.  The bins hold
    equal numbers of the plain lines that carry a displacement; a plain line
    with D = 0 tells nothing about kappa and only helps fix the level values.
    """
    informative = sorted(ln.wn for i, ln in enumerate(lines)
                         if keep[i] and ln.cls == 'plain' and ln.D)
    edges = [informative[int(round(k * len(informative) / nbins))]
             for k in range(1, nbins)]

    def bin_of(wn):
        return sum(1 for e in edges if wn >= e)

    def class_of(ln):
        return ('plain', bin_of(ln.wn)) if ln.cls == 'plain' else ln.cls

    fit = solve(lines, class_of, 'flag', 1.0, scale=scale, sigma=sigma,
                use=keep)
    rows = []
    for k in range(nbins):
        wns = [ln.wn for i, ln in enumerate(lines)
               if keep[i] and ln.cls == 'plain' and ln.D
               and bin_of(ln.wn) == k]
        key = ('plain', k)
        rows.append((min(wns), max(wns), len(wns),
                     fit['kappa'][key], fit['u_kappa'][key]))
    return fit, rows


def null_test(lines, sigma, keep, scale=1.0, seed=20260921):
    """Permute the displacements and refit: kappa must come out at zero.

    The displacements D are shuffled among the lines of each class, which
    destroys the pairing of a line with its own hyperfine geometry and leaves
    everything else alone - the level values, the weights, the class sizes.
    A kappa that survives this is being produced by the fit's own machinery
    rather than by hyperfine structure.
    """
    rng = np.random.default_rng(seed)
    byclass = collections.defaultdict(list)
    for i, ln in enumerate(lines):
        byclass[ln.cls].append(i)
    shuffled = list(lines)
    for cls, idx in byclass.items():
        order = list(idx)
        rng.shuffle(order)
        for i, j in zip(idx, order):
            shuffled[i] = lines[i]._replace(D=lines[j].D)
    return solve(shuffled, lambda ln: ln.cls, 'flag', 1.0, scale=scale,
                 sigma=sigma, use=keep)


def shared_level_test(lines, sigma, keep, constants=None, scale=1.0):
    """The era test again, on the levels the two eras have in common.

    If kappa(plain) came out different in the two eras only because the two
    eras see different levels - and so different A constants, of different
    provenance and different reliability - then confining the test to the
    levels that carry a displacement in BOTH eras must remove the difference.
    Whatever survives this is a property of the measurement, not of the
    constants.
    """
    constants = read_A_constants() if constants is None else constants

    def carriers(ln):
        out = set()
        for lid in (ln.low, ln.upp):
            rec = constants.get(lid)
            if rec is not None and rec[1]:
                out.add(lid)
        return out

    seen = {1974: set(), 1969: set()}
    for i, ln in enumerate(lines):
        if keep[i] and ln.cls == 'plain' and ln.D:
            seen[ln.era] |= carriers(ln)
    shared = seen[1974] & seen[1969]
    use = [keep[i] and (lines[i].cls != 'plain' or not lines[i].D
                        or carriers(lines[i]) <= shared)
           for i in range(len(lines))]

    def class_of(ln):
        return (ln.cls, ln.era) if ln.cls == 'plain' else ln.cls

    fit = solve(lines, class_of, 'flag', 1.0, scale=scale, sigma=sigma,
                use=use)
    counts = {}
    for era in (1974, 1969):
        counts[era] = sum(1 for i, ln in enumerate(lines)
                          if use[i] and ln.cls == 'plain' and ln.D
                          and ln.era == era)
    return fit, shared, counts


def set1_uncertainties(lines, path=None):
    """The `unc_wn_obs` of the line list, for comparison only.

    These are the Set 1 weights the plan's section 2.6 rules out for this
    purpose, because they are inflated by a line's own departure from its
    Ritz value.  They are read here so that the report can show what using
    them does to the answer, which is the evidence for that ruling.
    """
    path = path or os.path.join(HERE, LINES)
    quoted = {}
    with open(path, encoding='utf-8', newline='') as fh:
        for row in csv.DictReader(fh):
            if row['accepted'] == '1.0' and row['n_accepted'] == '1':
                try:
                    quoted['%.6f' % float(row['wn_obs'])] = float(
                        row['unc_wn_obs'])
                except ValueError:
                    pass
    out = []
    for ln in lines:
        u = quoted.get('%.6f' % ln.wn)
        out.append(u if u else stated_uncertainty(ln.char, ln.era, ln.wn))
    return out


def _fmt(value, u):
    return '%+.3f +- %.3f' % (value, u)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--lines', default=None)
    ap.add_argument('--scale', type=float, default=A_SCALE,
                    help='scale applied to every A constant; Step 2 measured '
                         '%.4f, and 1.0 uses the constants as calculated'
                         % A_SCALE)
    ap.add_argument('--outlier', type=float, default=OUTLIER_SIGMA)
    ap.add_argument('--nulls', type=int, default=20,
                    help='permutations of the displacements in control 1')
    ap.add_argument('--dry-run', action='store_true',
                    help='report without writing the two csv files')
    args = ap.parse_args(argv)

    constants = read_A_constants()
    J_of = read_levels()
    lines = read_lines(args.lines, constants, J_of)
    counts = collections.Counter(ln.cls for ln in lines)
    with_D = collections.Counter(ln.cls for ln in lines if ln.D)
    print('lines fitted %d, accepted and singly assigned' % len(lines))
    for cls in ('plain', 'flag', 'c'):
        print('   %-6s %5d   of them with a computed displacement %5d'
              % (cls, counts[cls], with_D[cls]))

    sigma, keep, table, fit = adopt_uncertainties(
        lines, scale=args.scale, outlier=args.outlier)
    print('\nadopted uncertainties, per character per era (Step 1)')
    print('   char  era   n_lines   stated    adopted   (cm^-1)')
    for (char, era), (adopted, n, started) in sorted(
            table.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        print('   %-5s %4d %8d   %8.4f %10.4f'
              % (char or "''", era, n, started, adopted))

    n_out = sum(1 for k in keep if not k)
    print('\nset aside as probably misassigned, beyond %.0f sigma: %d of %d'
          % (args.outlier, n_out, len(lines)))

    print('\nthe convention factors, A constants scaled by %.4f' % args.scale)
    print('   fit: %d lines, %d unknowns, rank %d, chi2/dof %.2f, rms %.4f'
          % (fit['n'], fit['unknowns'], fit['rank'],
             fit['chi2'] / max(fit['dof'], 1), fit['rms']))
    for a, b in (('plain', 'flag'), ('c', 'flag'), ('c', 'plain')):
        d, u = difference(fit, a, b)
        print('   kappa(%s) - kappa(%s)   %s' % (a, b, _fmt(d, u)))
    print('   anchored at kappa(flag) = 1:')
    for cls in ('plain', 'c'):
        print('      kappa(%s) = %s' % (cls, _fmt(fit['kappa'][cls],
                                                  fit['u_kappa'][cls])))

    print('\ntest 1 (section 3.3): kappa(plain) by era')
    fe = era_test(lines, sigma, keep, scale=args.scale)
    for key in (('plain', 1974), ('plain', 1969)):
        print('      %-14s %s' % ('%s %d' % key,
                                  _fmt(fe['kappa'][key], fe['u_kappa'][key])))
    d, u = difference(fe, ('plain', 1974), ('plain', 1969))
    print('      difference     %s   (%.1f sigma)' % (_fmt(d, u), abs(d / u)))

    print('      the pooled value above must lie between these two, and the')
    print('      test is whether they agree with it, not with each other')

    print('\nthe same in bins holding equal numbers of informative plain lines')
    fd, rows = drift_scan(lines, sigma, keep, scale=args.scale)
    print('      wavenumber range       lines   kappa(plain)')
    for lo, hi, n, k, u in rows:
        print('      %9.1f %9.1f %7d   %s' % (lo, hi, n, _fmt(k, u)))

    print('\ntest 2 (section 3.3): kappa(flag) by branch, dJ = +1 held at 1')
    fb = branch_test(lines, sigma, keep, scale=args.scale)
    for key in sorted(k for k in fb['kappa'] if isinstance(k, tuple)):
        n = sum(1 for i, ln in enumerate(lines)
                if keep[i] and ln.cls == 'flag' and ln.dJ is not None
                and int(round(ln.dJ)) == key[1])
        print('      dJ = %+d  %4d lines   %s'
              % (key[1], n, _fmt(fb['kappa'][key], fb['u_kappa'][key])))

    print('\ntest 1 again, on the 74 levels both eras have in common')
    fs, shared, counts = shared_level_test(lines, sigma, keep, constants,
                                           scale=args.scale)
    for era in (1974, 1969):
        key = ('plain', era)
        print('      plain %d  %5d lines   %s'
              % (era, counts[era], _fmt(fs['kappa'][key], fs['u_kappa'][key])))
    print('      %d levels carry a displacement in both eras' % len(shared))

    print('\ncontrol 1: the displacements permuted within each class, %d times' % args.nulls)
    draws = collections.defaultdict(list)
    for seed in range(args.nulls):
        fn = null_test(lines, sigma, keep, scale=args.scale, seed=20260921 + seed)
        for cls in ('plain', 'c'):
            draws[cls].append(fn['kappa'][cls])
    for cls in ('plain', 'c'):
        d = draws[cls]
        mean = sum(d) / len(d)
        sd = math.sqrt(sum((v - mean) ** 2 for v in d) / max(len(d) - 1, 1))
        print('      kappa(%s) = %s   against %s measured'
              % (cls, _fmt(mean, sd),
                 _fmt(fit['kappa'][cls], fit['u_kappa'][cls])))

    print('\ncontrol 2: the same fit on the Set 1 weights, which section 2.6 rules out')
    f1 = solve(lines, lambda ln: ln.cls, 'flag', 1.0, scale=args.scale,
               sigma=set1_uncertainties(lines), use=keep)
    for cls in ('plain', 'c'):
        print('      kappa(%s) = %s' % (cls, _fmt(f1['kappa'][cls],
                                                  f1['u_kappa'][cls])))

    if not args.dry_run:
        p1 = write_uncertainties(table)
        p2, n_bad = write_misassigned(lines, fit, sigma, keep)
        print('\nwritten: %s, %s (%d rows)'
              % (os.path.basename(p1), os.path.basename(p2), n_bad))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
