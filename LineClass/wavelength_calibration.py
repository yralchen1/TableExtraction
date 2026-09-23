#!/usr/bin/env python
"""The wavelength calibration of Sugar's observed line list.

What is being measured
----------------------
Sugar's wavenumbers were obtained from wavelengths measured on photographic
plates and reduced against comparison lines recorded on the same plate.  If
the comparison scale of a plate is slightly wrong, or the plate sat slightly
off its nominal position, every line of that plate is displaced by the same
amount *in wavelength* - that is what a plate calibration error is.  Call it
`delta_lambda`, in angstrom, positive when the reported wavelength is too
long.

A wavelength too long is a wavenumber too small, because `wn = 1e8 / lambda`
with `lambda` in angstrom and `wn` in cm^-1.  Differentiating,

    delta_wn = -delta_lambda * wn^2 * 1e-8

so the model this program fits to every accepted, singly assigned line is

    wn_measured = (E_upp - E_low + kappa * D) - delta_lambda * wn^2 * 1e-8

where `E_upp` and `E_low` are the two level energies, free parameters of the
same fit, and `kappa * D` is the hyperfine displacement of the line: `D` is
the calculated separation between the center of gravity of the hyperfine
pattern and the point of it that Sugar's measurement actually followed, and
`kappa` is the fraction of that separation his measuring convention picked up
(`hfs_kappa.py`, plan Step 3).  One `kappa` is carried for the plain lines of
1974, one for the plain lines of 1969, and one for the lines Sugar marked
complex, `c`; the lines he flagged `*r` and `*v` have their measured
displacement removed directly and carry no `kappa`.

Where a plate ends
------------------
Sugar's exposure boundaries were never published, but the line list maps
them: at a join no line of any species was recorded, so the list goes blind
over a stretch of wavelength.  `tools/coverage_map.py` has measured those
stretches into `coverage_gaps.txt`, and `tools/calibrate_intensities.py`
already imposes them on the intensity scale.  They are imposed here on the
wavelength scale for the same reason, and a stronger one: across a join no
line was recorded twice, so nothing whatever ties the two calibrations
together.  Carrying a correction across a join is not a smoothing
approximation, it is an invention.

Each blind stretch is a block in its own right, not a hole.  It is thin, not
empty, and the few lines that were recorded inside it came from whatever
exposure did reach there - which is a plate of its own, and certainly not
either of its neighbors.

How `delta_lambda` runs inside a block
--------------------------------------
Two models are carried, and both are fitted on every run so that they can be
compared block by block.

`--model groups`, the default, is piecewise constant.  Inside a block,
`delta_lambda` is taken to be constant over `BIN_A` angstrom; walking a block
from its blue end, those bins are accumulated until the group holds
`MIN_IN_BIN` lines, and a short group left over at the red end is merged back
into its predecessor.  A group therefore never spans a join, and a block that
cannot raise `MIN_IN_BIN` lines in total is given no parameter at all: its
calibration is unmeasurable, which is an honest answer, where the neighboring
plate's value would be a wrong one.

`--model poly` puts one polynomial in wavelength on each block, with the
degree chosen by the data: a degree is added while it lowers the total
chi-square by more than `DCHI2`, up to `DEG_MAX`, and a block may carry one
coefficient per `MIN_IN_BIN` of its lines.  The basis is Legendre in the
block's own wavelength range, which keeps the coefficients nearly
uncorrelated.  This is the form a calibration error is usually written in,
and it would be the better model if a block were a single smooth reduction.

Which one the spectrum actually asks for is a measurement, not an assumption,
and the report prints it: the two chi-squares of every block side by side
against the number of parameters between them.  A block whose plate is smooth
gives back about one unit of chi-square per parameter the staircase spends; a
block that gives back much more is telling us that its displacement moves
faster than any low polynomial can follow - an undetected join inside it, or
a reduction done in pieces - and there the staircase is the honest model.

What comes out
--------------
`wavelength_calibration.txt`  the report: the blocks, the fitted curve, the
                              diagnostics, and the calibration uncertainty.
`wavelength_calibration.csv`  the same curve, one row per group, machine
                              readable.
`wavelength_calibration_corrections.csv`
                              the deliverable: the correction to add to
                              every observed wavenumber of the list, with
                              its uncertainty from the covariance matrix of
                              this fit.  Every line, not only the ones the
                              fit rests on.
`wavelength_calibration_poly.csv`
                              (`--model poly` only) the fitted coefficients
                              and their full covariance.
`wavelength_calibration_points.csv`
                              the input of the fit, one row per line: the
                              wavelength, the displacement in angstrom that
                              line asks for, and its uncertainty.
`wavelength_calibration_fit/` the same points as tab-delimited input files
                              for `fit_power.py`, one per block, weighted by
                              `u_eff_A`, so that the wavelength dependence
                              inside a plate can be fitted as a smooth
                              function rather than as a staircase.  A block
                              too thin to carry a fitted group here gets a
                              file all the same: its points
                              are measurements of its own plate, and a
                              one-parameter fit to them is worth more than
                              nothing.

The per-line displacement is

    dl_i = (E_upp - E_low + kappa * D - wn_i) / (wn_i^2 * 1e-8)

with the level energies and the kappas of this fit, and its uncertainty is
the line's adopted uncertainty on the same scale,
`u_i = sigma_i / (wn_i^2 * 1e-8)`.  These are the numbers the fitted
`delta_lambda` of a group is the weighted mean of.

One caution about them.  The level energies are fitted to these same lines,
so a point is not independent of the level values it is compared with; the
leverage `h_i` of the fit says how much of it the levels have already
absorbed, the expected square of a residual being `sigma_i^2 * (1 - h_i)`
rather than `sigma_i^2`.  The leverage is written out with each point, and
so is the effective uncertainty

    u_eff_i = u_i / (1 - h_i)

which is what an outside fit to these points should weight them by.  The
two effects it has to cover pull opposite ways: the scatter of the points
about their own group is smaller than `u_i` because the levels have followed
them, while the uncertainty of the group's mean is larger than an
independent average of `u_i` would give because the points share those
levels.  No per-point uncertainty can reproduce a correlated covariance
exactly, but this one comes close: taken over the 151 groups of this fit,
the weighted mean of `u_eff_i` reproduces the group uncertainty of the full
covariance matrix to a median of 0.97, where the plain `u_i` gives 0.74.
It also disposes of the `h_i = 1` lines, which carry nothing.

The anchor
----------
The fit is reported twice, once with every group free and once with the
best-populated group in the 1974 region held at zero.  The comparison is
printed: if freeing the anchor raises the rank of the design matrix by one,
the curve is determined only up to an additive constant and the anchor is
needed to name it; if the rank does not move, the anchor is an unnecessary
constraint and the free solution is the one to quote.  The energy zero is a
degeneracy of the level system in any case - the whole odd parity can be
displaced against the whole even parity without changing a single predicted
wavenumber - so a rank deficiency of exactly one is expected and harmless.
"""
import argparse
import collections
import csv
import io
import math
import os

import numpy as np

import hfs_kappa

HERE = os.path.dirname(os.path.abspath(__file__))

#: the width of the finest bin, in angstrom, and the fewest lines a fitted
#: group may rest on.  25 angstrom is well inside a plate and 12 lines is
#: about where the weighted mean of a group stops being dominated by its
#: single best line.
BIN_A = 25.0
MIN_IN_BIN = 12

#: a block that cannot raise MIN_IN_BIN lines still gets one constant over
#: the whole block, provided it has at least this many lines.  Two lines
#: measure a mean; one line measures nothing, since the fit can satisfy it
#: by moving a level instead.  Such a block is one group and never more.
MIN_THIN = 2

#: a line whose leverage reaches this carries no information about the
#: calibration at all: one of its levels rests on that line alone, so the fit
#: satisfies it exactly whatever `delta_lambda` is.  Such a line is written
#: into `wavelength_calibration_points.csv` with an empty `u_eff_A` and is
#: left out of the `fit_power.py` input files.
H_MAX = 0.999

GAPS = os.path.join(HERE, 'coverage_gaps.txt')
OUT_REPORT = os.path.join(HERE, 'wavelength_calibration.txt')
OUT_CURVE = os.path.join(HERE, 'wavelength_calibration.csv')
OUT_POINTS = os.path.join(HERE, 'wavelength_calibration_points.csv')
OUT_CORR = os.path.join(HERE, 'wavelength_calibration_corrections.csv')
OUT_POLY = os.path.join(HERE, 'wavelength_calibration_poly.csv')
OUT_FIT = os.path.join(HERE, 'wavelength_calibration_fit')

#: the highest degree a block's polynomial may reach, and the drop in the
#: total chi-square that buys one more degree.  9 is 3 sigma on the one
#: parameter added; a block may carry one parameter per MIN_IN_BIN lines.
DEG_MAX = 5
DCHI2 = 9.0

#: the anchor is looked for above this wavelength, in the middle of the 1974
#: region where the lines are densest and Sugar's stated uncertainty smallest.
ANCHOR_ABOVE = 2190.0


def read_blocks(path=GAPS):
    """The stretches of wavelength the blind ones cut the spectrum into.

    Returns a list of `(lam_lo, lam_hi)` covering the whole axis without
    overlap: the observed stretches and the blind stretches alternate, and
    each blind stretch is a block in its own right.
    """
    cuts = []
    with open(path, encoding='utf-8') as fh:
        for row in fh:
            if row.startswith('#'):
                continue
            field = row.split()
            if field:
                cuts.append((float(field[0]), float(field[1])))
    cuts.sort()
    edges, prev = [], 0.0
    for lo, hi in cuts:
        edges.append((prev, lo))
        edges.append((lo, hi))
        prev = hi
    edges.append((prev, 1e9))
    return edges


def block_of(blocks, value):
    """The block a wavelength falls in."""
    for k, (lo, hi) in enumerate(blocks):
        if lo <= value < hi:
            return k
    return len(blocks) - 1


def observed_wavenumbers(path=None):
    """Every distinct observed wavenumber of the line list, with its character.

    The fit rests on the accepted, singly assigned lines only, but the
    correction is a property of the plate and applies to every line recorded
    on it - the unclassified ones above all, which is where the next
    identification has to be made.
    """
    path = path or os.path.join(hfs_kappa.HERE, hfs_kappa.LINES)
    seen = {}
    with open(path, encoding='utf-8', newline='') as fh:
        for row in csv.DictReader(fh):
            wn = float(row['wn_obs'])
            seen.setdefault(round(wn, 6), (wn, row['char']))
    return [seen[k] for k in sorted(seen)]


def group_bins(blocks, lam, keep):
    """Assign every line to a fitted group, or to none.

    Returns `(block_of, group_of, counts, span)`: the block index of each
    line, the group key of each line (`None` where the block carries no
    parameter), the number of kept lines behind each group, and the
    wavelength range each group covers, clipped to its own block.
    """
    def block_of(value):
        for k, (lo, hi) in enumerate(blocks):
            if lo <= value < hi:
                return k
        return len(blocks) - 1

    blk = [block_of(value) for value in lam]
    raw = [(blk[i], int(lam[i] // BIN_A)) for i in range(len(lam))]
    raw_counts = collections.Counter(raw[i] for i in range(len(lam))
                                     if keep[i])
    per_block = collections.Counter()
    for (b, _), n in raw_counts.items():
        per_block[b] += n

    group, counts, span = {}, collections.Counter(), {}
    for b in sorted(per_block):
        steps = sorted(s for (bb, s) in raw_counts if bb == b)
        if per_block[b] < MIN_IN_BIN:
            if per_block[b] < MIN_THIN:
                continue
            key = (b, steps[0])
            for st in steps:
                group[(b, st)] = key
            counts[key] = per_block[b]
            span[key] = (max(steps[0] * BIN_A, blocks[b][0]),
                         min((steps[-1] + 1) * BIN_A, blocks[b][1]))
            continue
        made, cur, n = [], [], 0
        for s in steps:
            cur.append(s)
            n += raw_counts[(b, s)]
            if n >= MIN_IN_BIN:
                made.append((cur, n))
                cur, n = [], 0
        if cur:
            if made:
                made[-1] = (made[-1][0] + cur, made[-1][1] + n)
            else:
                made.append((cur, n))
        for cur, n in made:
            key = (b, cur[0])
            for s in cur:
                group[(b, s)] = key
            counts[key] = n
            span[key] = (max(cur[0] * BIN_A, blocks[b][0]),
                         min((cur[-1] + 1) * BIN_A, blocks[b][1]))
    return blk, [group.get(raw[i]) for i in range(len(lam))], counts, span



def legendre_row(deg, t):
    """`[P_0(t), ..., P_deg(t)]`, the Legendre polynomials on `[-1, 1]`.

    They are used instead of the plain powers `1, t, t^2, ...` because they
    are nearly uncorrelated over a well filled block, which keeps the normal
    matrix well conditioned and lets a fitted coefficient, and its
    uncertainty, be read one at a time.
    """
    out = [1.0]
    if deg >= 1:
        out.append(t)
    for k in range(1, deg):
        out.append(((2 * k + 1) * t * out[k] - k * out[k - 1]) / (k + 1))
    return out


def model_groups(good):
    """The staircase: one constant `delta_lambda` per group.

    Returns `(n_columns, terms)`, where `terms(block, group, lambda)` gives
    the `(column, value)` pairs a line contributes to the calibration part of
    the design matrix.  A line in a block that carries no parameter
    contributes nothing and so receives no correction.
    """
    idx = {g: k for k, g in enumerate(good)}

    def terms(b, g, lam_value):
        return [(idx[g], 1.0)] if g in idx else []

    return len(good), terms


def poly_keys(degrees, omit=None):
    """The calibration columns of the polynomial model, in order."""
    return [(b, d) for b in sorted(degrees)
            for d in range(degrees[b] + 1) if (b, d) != omit]


def model_poly(degrees, ranges, omit=None):
    """A polynomial in wavelength per block, piecewise across the joins.

    `degrees[b]` is the degree carried by block `b` and `ranges[b]` the
    wavelength interval its lines cover, which the Legendre argument is
    scaled to.  `omit` drops one coefficient; it is used only to test whether
    the model is determined absolutely or only up to that coefficient.
    """
    idx = {kk: k for k, kk in enumerate(poly_keys(degrees, omit))}

    def terms(b, g, lam_value):
        if b not in degrees:
            return []
        lo, hi = ranges[b]
        t = 0.0 if hi <= lo else 2.0 * (lam_value - lo) / (hi - lo) - 1.0
        p = legendre_row(degrees[b], t)
        return [(idx[(b, d)], p[d]) for d in range(degrees[b] + 1)
                if (b, d) != omit]

    return len(idx), terms


def build(lines, lam, blk, binof, sigma, keep, ncal, terms):
    """The design matrix, the right-hand side and the weights.

    The columns are the level energies, then the `ncal` calibration
    parameters laid out by `terms`, then the three `kappa` of the hyperfine
    convention.
    """
    rows = [i for i in range(len(lines)) if keep[i]]
    levels = sorted(set(lines[i].low for i in rows)
                    | set(lines[i].upp for i in rows))
    lidx = {v: k for k, v in enumerate(levels)}
    off = len(levels)
    kcls = ['plain_1974', 'plain_1969', 'c']
    kidx = {c: off + ncal + k for k, c in enumerate(kcls)}

    A = np.zeros((len(rows), off + ncal + len(kcls)))
    y = np.zeros(len(rows))
    w = np.zeros(len(rows))
    for r, i in enumerate(rows):
        ln = lines[i]
        A[r, lidx[ln.upp]] += 1.0
        A[r, lidx[ln.low]] -= 1.0
        scale = ln.wn * ln.wn * 1e-8
        for col, val in terms(blk[i], binof[i], lam[i]):
            A[r, off + col] -= scale * val
        D = hfs_kappa.A_SCALE * ln.D
        if ln.cls == 'flag':
            y[r] = ln.wn - D
        else:
            A[r, kidx['c' if ln.cls == 'c' else 'plain_%d' % ln.era]] += D
            y[r] = ln.wn
        w[r] = 1.0 / sigma[i]
    return rows, levels, lidx, off, kidx, A, y, w


def solve(A, y, w):
    """One weighted least-squares solve, with leverage and covariance.

    The covariance is the pseudo-inverse of the weighted normal matrix, which
    is the right thing here because the level system is rank deficient by the
    energy zero: the pseudo-inverse puts no variance along that direction and
    gives the correct variance of every estimable function, which every
    calibration quantity in this program is.  It is symmetric by
    construction; it is symmetrized explicitly all the same, and how far it
    had drifted is reported.
    """
    Aw = A * w[:, None]
    sol, _, rank, _ = np.linalg.lstsq(Aw, y * w, rcond=None)
    pinv = np.linalg.pinv(Aw.T.dot(Aw))
    asym = float(np.max(np.abs(pinv - pinv.T)))
    big = float(np.max(np.abs(pinv)))
    pinv = 0.5 * (pinv + pinv.T)
    lev = np.diag(Aw.dot(pinv).dot(Aw.T))
    res = y - A.dot(sol)
    chi2 = float(np.sum((res * w) ** 2))
    return dict(sol=sol, rank=rank, cov=pinv, leverage=lev, residual=res,
                chi2=chi2, dof=len(y) - rank,
                asym=(asym / big if big else 0.0))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--model', choices=('groups', 'poly'), default='groups',
                    help='groups: one constant per group of at least %d '
                         'lines, piecewise in wavelength (default).  poly: '
                         'one polynomial in wavelength per block, its degree '
                         'chosen by the data.' % MIN_IN_BIN)
    ap.add_argument('--degree', type=int, default=None,
                    help='force this degree on every parameterized block '
                         'instead of letting the data choose it')
    ap.add_argument('--no-write', action='store_true',
                    help='print the report without writing any file')
    args = ap.parse_args(argv)

    blocks = read_blocks()
    lines = hfs_kappa.read_lines()
    sigma, keep, table, _ = hfs_kappa.adopt_uncertainties(
        lines, scale=hfs_kappa.A_SCALE)
    lam = [1e8 / ln.wn for ln in lines]
    blk, binof, counts, span = group_bins(blocks, lam, keep)
    good = sorted(counts)

    def run(ncal, terms):
        parts = build(lines, lam, blk, binof, sigma, keep, ncal, terms)
        return parts, solve(*parts[5:])

    def extent(good):
        """The blocks that carry a parameter, their line counts and the
        wavelength interval their own lines cover - which is what the
        Legendre argument is scaled to, not the nominal block edges, which
        run into the middle of the blind stretches on either side."""
        pb = sorted(set(g[0] for g in good))
        n = collections.Counter()
        lo_of, hi_of = {}, {}
        for i in range(len(lines)):
            b = blk[i]
            if not keep[i] or b not in pb:
                continue
            n[b] += 1
            lo_of[b] = min(lo_of.get(b, lam[i]), lam[i])
            hi_of[b] = max(hi_of.get(b, lam[i]), lam[i])
        return pb, n, {b: (lo_of[b], hi_of[b]) for b in pb}

    # ---- the thin blocks, on trial ---------------------------------------
    # A block too thin for a group of its own is given one constant over the
    # whole block and then asked whether that constant is worth keeping.  The
    # test is the parameter's own significance: if the uncertainty of the
    # constant is larger than the constant itself, the block has not measured
    # anything and the parameter is fixed at zero.  Keeping such a parameter
    # is not free - it trades almost exactly against the energies of the few
    # levels its lines touch, and that near-degeneracy leaks through the
    # common mode into the uncertainty of every other group in the spectrum.
    # The lines of a fixed block still get an uncertainty - the one measured
    # here - against a correction of zero.
    thin = [g for g in good if counts[g] < MIN_IN_BIN]
    thin_u = {}
    if thin:
        tcal, tterms = model_groups(good)
        tparts, tfit = run(tcal, tterms)
        toff = tparts[3]
        tidx = {g: k for k, g in enumerate(good)}
        drop = []
        for g in thin:
            k = toff + tidx[g]
            v = float(tfit['sol'][k])
            u = math.sqrt(max(float(tfit['cov'][k, k]), 0.0))
            thin_u[g[0]] = (v, u, counts[g])
            if not v or abs(u / v) >= 1.0:
                drop.append(g)
        if drop:
            good = [g for g in good if g not in set(drop)]

    parblocks, nblk, ranges = extent(good)

    # ---- the staircase, always fitted: it is the reference the polynomial
    # ---- has to match, block by block, before it may replace it ----------
    gcal, gterms = model_groups(good)
    gparts, gfit = run(gcal, gterms)

    # ---- the polynomial, one degree at a time ---------------------------
    search = []
    degrees = {b: 0 for b in parblocks}
    if args.degree is not None:
        for b in parblocks:
            degrees[b] = min(args.degree, DEG_MAX,
                             max(0, nblk[b] // MIN_IN_BIN - 1))
    ncal, terms = model_poly(degrees, ranges)
    parts, cur = run(ncal, terms)
    if args.degree is None:
        # Every degree up to the block's ceiling is tried, not only degrees
        # up to the first one that fails.  A plate whose error is a symmetric
        # bow gains nothing from a slope and a great deal from a curve, so a
        # search that stops at the first disappointment never finds it.  The
        # price of a degree is DCHI2 per parameter it adds, and the winner is
        # the degree that beats that price by the widest margin.  The blocks
        # are coupled through the level energies they share, so the scan is
        # repeated until no block changes its mind.
        for _pass in range(3):
            moved = False
            for b in parblocks:
                top = min(DEG_MAX, max(0, nblk[b] // MIN_IN_BIN - 1))
                best, gain = degrees[b], 0.0
                for d in range(degrees[b] + 1, top + 1):
                    trial = dict(degrees)
                    trial[b] = d
                    nc2, t2 = model_poly(trial, ranges)
                    p2, f2 = run(nc2, t2)
                    drop = cur['chi2'] - f2['chi2']
                    need = DCHI2 * (d - degrees[b])
                    search.append((b, d, drop, drop > need))
                    if drop - need > gain:
                        best, gain = d, drop - need
                        keep_ = (nc2, t2, p2, f2)
                if best != degrees[b]:
                    degrees[b] = best
                    ncal, terms, parts, cur = keep_
                    moved = True
            if not moved:
                break

    if args.model == 'groups':
        ncal, terms, parts, fit = gcal, gterms, gparts, gfit
    else:
        fit = cur
    rows, levels, lidx, off, kidx, A, y, w = parts
    sol, cov = fit['sol'], fit['cov']
    ncol = A.shape[1]

    def cvec(b, g, lam_value):
        """The row of the design matrix that reads off `delta_lambda`.

        `delta_lambda = c . sol` and `u^2 = c . cov . c`, which is the only
        correct way to take the uncertainty of a polynomial's value: the
        coefficients are correlated with one another, and the calibration is
        correlated with the level energies.
        """
        c = np.zeros(ncol)
        for col, val in terms(b, g, lam_value):
            c[off + col] = val
        return c

    def value(c):
        return (float(c.dot(sol)),
                math.sqrt(max(float(c.dot(cov).dot(c)), 0.0)))

    goodset = set(good)
    has_par = [(blk[i] in degrees) if args.model == 'poly'
               else (binof[i] in goodset) for i in range(len(lines))]

    # ---- is the curve determined absolutely, or only up to a constant? ---
    anchor = max((g for g in good if span[g][0] >= ANCHOR_ABOVE),
                 key=lambda g: counts[g])
    if args.model == 'groups':
        held = run(*model_groups([g for g in good if g != anchor]))[1]
    else:
        held = run(*model_poly(degrees, ranges, omit=(anchor[0], 0)))[1]
    anchored = fit['rank'] == held['rank']

    # ---- what the one rank deficiency actually is ------------------------
    _, sv, vt = np.linalg.svd(A * w[:, None], full_matrices=False)
    nullv = vt[int(np.argmin(sv))]
    null = (max(abs(nullv[off + k]) for k in range(ncal)) if ncal else 0.0,
            max(abs(nullv[k]) for k in range(len(levels))))

    out = io.StringIO()

    def say(text=''):
        out.write(text + '\n')

    say('wavelength_calibration.py - the plate calibration of the observed '
        'line list')
    say('the displacement delta_lambda is in angstrom, positive where '
        "Sugar's wavelength is too long;")
    say('the correction to add to his wavenumber is '
        '-delta_lambda * wn^2 * 1e-8, so a positive')
    say('delta_lambda means his wavenumber is too small.')
    say()
    say('model: %s' % ('one polynomial in wavelength per block, its degree '
                       'chosen by the data'
                       if args.model == 'poly'
                       else 'one constant per group of at least %d lines, '
                            'and one over any thinner block that has at '
                            'least %d' % (MIN_IN_BIN, MIN_THIN)))
    say()
    if thin_u:
        say('The thin blocks, on trial')
        say('-------------------------')
        say('A block that cannot raise %d lines for a group of its own was '
            'given one' % MIN_IN_BIN)
        say('constant over the whole block and then asked whether that '
            'constant is worth')
        say('keeping.  The test is the parameter against its own '
            'uncertainty: a constant')
        say('smaller than its uncertainty has measured nothing and is fixed '
            'at zero.  Such')
        say('a block trades its constant almost exactly against the energies '
            'of the few')
        say('levels its lines touch, and keeping it lets that near-degeneracy '
            'leak through')
        say('the common mode into every other group in the spectrum: with '
            'all six kept the')
        say('largest eigenvalue of the calibration covariance rises by a '
            'factor of 66.')
        say('The lines of a fixed block still carry the uncertainty measured '
            'here, against')
        say('a correction of zero.')
        say('%5s %6s %13s %11s %8s   %s'
            % ('blk', 'lines', 'd_lambda_A', 'u_A', '|u/d|', 'verdict'))
        for b in sorted(thin_u):
            v, u, n = thin_u[b]
            r = abs(u / v) if v else float('inf')
            say('%5d %6d %13.5f %11.5f %8.2f   %s'
                % (b, n, v, u, r,
                   'kept' if r < 1.0 else 'fixed at zero'))
        say()
    say('The blocks, cut at the measured plate joins of coverage_gaps.txt')
    say('-----------------------------------------------------------------')
    say('%10s %10s %8s %6s %6s' % ('from_A', 'to_A', 'lines', 'groups',
                                   'degree'))
    per_block = collections.Counter()
    for i in range(len(lines)):
        if keep[i]:
            per_block[blk[i]] += 1
    for k, (lo, hi) in enumerate(blocks):
        n = per_block.get(k, 0)
        if not n:
            continue
        nb = len([g for g in good if g[0] == k])
        say('%10.2f %10.2f %8d %6d %6s%s'
            % (lo, min(hi, 99999.0), n, nb,
               degrees[k] if k in degrees else '-',
               '' if nb else '   no parameter'))

    say()
    say('The fit')
    say('-------')
    say('lines %d, unknowns %d, rank %d, chi2/dof %.3f, rms %.4f cm-1'
        % (len(rows), ncol, fit['rank'], fit['chi2'] / fit['dof'],
           math.sqrt(float(np.mean(fit['residual'] ** 2)))))
    say('calibration parameters %d, levels %d, kappa 3'
        % (ncal, len(levels)))
    say('with every calibration parameter free: rank %d' % fit['rank'])
    say('with one of them held at zero:          rank %d' % held['rank'])
    if anchored:
        say('the two ranks agree, so the curve is determined only up to an '
            'additive constant')
        say('and one parameter would have to be held to name it.')
    else:
        say('the rank falls by one when a parameter is held, so every '
            'parameter is determined')
        say('by the data: the curve is absolute and nothing is held.')
    say('kappa: ' + '  '.join(
        '%s %+.3f+-%.3f' % (c, sol[j], math.sqrt(cov[j, j]))
        for c, j in kidx.items()))
    say('the one remaining rank deficiency is the energy zero: every level '
        'can be raised by')
    say('the same amount without changing a single level difference, so '
        'nothing observable')
    say('moves along it.  The direction was found from the design matrix '
        'itself, and the')
    say('largest weight it puts on any calibration parameter is %.2e '
        'against %.2e on the' % null)
    say('levels, so the curve is not touched by it and its uncertainties '
        'are not inflated by it.')
    say()
    say('The covariance matrix is the pseudo-inverse of the weighted normal '
        'matrix, formed')
    say('once for the whole fit, so every uncertainty below carries the '
        'correlation between')
    say('the calibration and the level energies, and between one '
        'coefficient and the next.')
    say('Its largest departure from symmetry before it was symmetrized was '
        '%.2e of its' % fit['asym'])
    say('largest element - rounding, as it must be.  The leverage inflation '
        'u/(1-h) is not')
    say('used here and cannot distort it: the fit weights each line by that '
        "line's own")
    say('uncertainty, and the leverage is computed from the solution, never '
        'fed back into it.')
    calcov = cov[off:off + ncal, off:off + ncal]
    ev = np.linalg.eigvalsh(calcov)
    say('The calibration block of it, %d by %d, has eigenvalues from %.3e to '
        '%.3e: all' % (ncal, ncal, float(ev.min()), float(ev.max())))
    say('positive, so it is a proper covariance matrix and any linear '
        'combination of the')
    say('calibration parameters has a positive variance.')

    say()
    say('Is a plate smooth?  The polynomial against the staircase')
    say('-------------------------------------------------------')
    say('Every block is fitted both ways.  chi2 is that block\'s share of '
        'the total, and the')
    say('difference is what the staircase buys with its extra parameters.  '
        'If the plate is')
    say('smooth, the staircase is only fitting noise and d_chi2 is about '
        'd_par; if a block')
    say('needs more freedom than its polynomial has, d_chi2 is much the '
        'larger of the two.')
    say('%4s %8s %6s %7s %11s %11s %8s %6s'
        % ('blk', 'lines', 'deg', 'groups', 'chi2_poly', 'chi2_step',
           'd_chi2', 'd_par'))
    grows = {i: r for r, i in enumerate(gparts[0])}
    prows = {i: r for r, i in enumerate(rows)}
    tot = [0.0, 0.0, 0]
    for b in parblocks:
        sel = [i for i in range(len(lines)) if keep[i] and blk[i] == b]
        c_p = sum((cur['residual'][prows[i]] / sigma[i]) ** 2
                  for i in sel if i in prows)
        c_g = sum((gfit['residual'][grows[i]] / sigma[i]) ** 2
                  for i in sel if i in grows)
        ng = len([g for g in good if g[0] == b])
        dpar = ng - (degrees[b] + 1)
        tot[0] += c_p
        tot[1] += c_g
        tot[2] += dpar
        say('%4d %8d %6d %7d %11.1f %11.1f %8.1f %6d'
            % (b, nblk[b], degrees[b], ng, c_p, c_g, c_p - c_g, dpar))
    say('%4s %8d %6s %7d %11.1f %11.1f %8.1f %6d'
        % ('all', sum(nblk.values()), '', len(good), tot[0], tot[1],
           tot[0] - tot[1], tot[2]))
    if search:
        say()
        say('The degrees were chosen one block at a time, and up '
            'to the ceiling of the block:')
        say('every degree is tried, not only the degrees up to the first one that '
            'fails: a plate')
        say('whose error is a symmetric bow gains nothing from a slope and a '
            'great deal from')
        say('a curve, and a search that stopped at the first disappointment '
            'would never find')
        say('it.  A degree costs %.1f in chi2 per parameter it adds, which '
            'is 3 sigma each,' % DCHI2)
        say('and the winner is the degree that beats that price by the '
            'widest margin; a block')
        say('may carry a parameter only per %d of its lines.  "pays" below '
            'is whether that one' % MIN_IN_BIN)
        say('degree covers its own cost, not whether it won: the fitted '
            'degree is in the')
        say('table above.  The scan is repeated until no block changes its '
            'mind, since the')
        say('blocks are coupled through the level energies they share.')
        say('%4s %6s %10s %6s' % ('blk', 'degree', 'd_chi2', 'pays'))
        for b, d, drop, kept in search:
            say('%4d %6d %10.1f %6s' % (b, d, -drop, 'yes' if kept else 'no'))

    say()
    say('The curve, at the midpoint of each group')
    say('----------------------------------------')
    say('%4s %12s %6s %12s %10s %12s %10s %7s'
        % ('blk', 'lambda_A', 'n', 'd_lambda_A', 'u_A', 'shift_cm-1',
           'u_cm-1', 'sigma'))
    curve, gvec = [], []
    for g in good:
        lo_a, hi_a = span[g]
        mid = 0.5 * (lo_a + hi_a)
        c = cvec(g[0], g, mid)
        v, u = value(c)
        shift = -v * 1e8 / (mid * mid)
        ushift = u * 1e8 / (mid * mid)
        curve.append((g[0], lo_a, hi_a, mid, counts[g], v, u, shift, ushift))
        gvec.append(c)
        say('%4d %5.0f-%6.0f %6d %+12.5f %10.5f %+12.4f %10.4f %7.1f'
            % (g[0], lo_a, hi_a, counts[g], v, u, shift, ushift,
               (v / u if u else 0.0)))

    say()
    say('How large the calibration error is, and how well it is known')
    say('------------------------------------------------------------')
    say('%-24s %6s %8s %11s %11s %11s %11s'
        % ('region', 'groups', 'at 2sig', 'mean_d_lam', 'rms_d_lam',
           'median_u', 'max|d|'))

    def common_mode(sel):
        """The mean displacement over a set of group midpoints.

        Two of these values are not independent - they share the level values
        the whole fit rests on, and inside a block they share the
        polynomial's coefficients - so the uncertainty of their mean has to
        come from the covariance matrix, not from an individual uncertainty
        divided by a square root.
        """
        if not sel:
            return 0.0, 0.0
        return value(sum(sel) / len(sel))

    eras = [('1969  below 2105 A', lambda a: a < 2105.0),
            ('1974  2105 - 4500 A', lambda a: 2105.0 <= a < 4500.0),
            ('1974  above 4500 A', lambda a: a >= 4500.0)]
    for name, test in eras:
        sel = [c for c in curve if test(c[3])]
        if not sel:
            continue
        rms = math.sqrt(sum(c[5] ** 2 for c in sel) / len(sel))
        med = sorted(c[6] for c in sel)[len(sel) // 2]
        big = max(abs(c[5]) for c in sel)
        n2 = len([c for c in sel if c[6] and abs(c[5]) > 2.0 * c[6]])
        say('%-24s %6d %8d %+11.5f %11.5f %11.5f %11.5f'
            % (name, len(sel), n2, sum(c[5] for c in sel) / len(sel),
               rms, med, big))
    allsel = [c for c in curve if c[6]]
    rms_all = math.sqrt(sum(c[5] ** 2 for c in allsel) / len(allsel))
    say()
    say('Over the whole spectrum the calibration error has an rms of '
        '%.5f angstrom,' % rms_all)
    say('measured at %d group midpoints, of which %d differ from zero by '
        'more than twice'
        % (len(allsel), len([c for c in allsel
                             if abs(c[5]) > 2.0 * c[6]])))
    say('their own uncertainty.')
    say('A calibration error is not noise: at a given wavelength it '
        'displaces every line the')
    say('same way, so it does not average down over the lines of a level '
        'and belongs in a')
    say('LOPT run as a systematic uncertainty, not as a widening of the '
        'individual line')
    say('uncertainties.')

    say()
    say('The common mode: the zero point of the whole wavelength scale')
    say('-------------------------------------------------------------')
    say('%-24s %6s %13s %13s %7s'
        % ('region', 'points', 'mean_d_lam_A', 'u_A', 'sigma'))
    for name, test in eras + [('all', lambda a: True)]:
        sel = [gvec[k] for k, c in enumerate(curve) if test(c[3])]
        m, um = common_mode(sel)
        say('%-24s %6d %+13.5f %13.5f %7.1f'
            % (name, len(sel), m, um, (m / um if um else 0.0)))
    say('The mean over all groups is what the fit says about the wavelength '
        'scale as a')
    say('whole, and it is not held to zero by anything: the shape of the '
        'displacement in')
    say('wavenumber, wn^2, is not a shape the level values can imitate, so '
        'the network of')
    say('Ritz combinations measures it.  It is measured weakly, though, '
        'which is what the')
    say('uncertainty in this table says; the differences between one part '
        'of the curve and')
    say('another are much better determined than their common level.')
    say()
    say('The statistical uncertainty of a corrected wavenumber')
    say('-----------------------------------------------------')
    say('Once the calibration error is taken out of a line, what is left is '
        "Sugar's measuring")
    say('scatter, and it is that - not his stated uncertainty, and not the '
        'uncertainty adopted')
    say('before the curve was known - that belongs in the LOPT input as the '
        'uncertainty of a')
    say('corrected wavenumber.  The table is over the lines of '
        'parameterized blocks only, the')
    say('rest having no correction to remove.  "after" is the rms of the '
        'residuals of this')
    say('fit over the effective number of degrees of freedom of the class, '
        'sum(1 - h), which')
    say('is what is left of its lines once the fitted level values have '
        'taken their share;')
    say('the floor of %.4f cm-1 is hfs_kappa.FLOOR.' % hfs_kappa.FLOOR)
    say()
    say('%-10s %6s %7s %11s %11s %11s %8s'
        % ('char', 'era', 'n', 'Sugar_A', 'adopted', 'after', 'factor'))
    ucls = hfs_kappa.uncertainty_classes(lines)
    classes = collections.defaultdict(list)
    for r, i in enumerate(rows):
        if has_par[i]:
            classes[ucls[i]].append((fit['residual'][r], fit['leverage'][r],
                                     lines[i]))
    for key in sorted(classes, key=lambda k: (-k[1], str(k[0]))):
        v = classes[key]
        dof = max(sum(1.0 - h for _, h, _ in v), 1.0)
        rms = max(math.sqrt(sum(e * e for e, _, _ in v) / dof),
                  hfs_kappa.FLOOR)
        ch, era = key
        stated = sum(hfs_kappa.stated_uncertainty(
            ch if ch != 'other' else 'w', era, ln.wn) for _, _, ln in v)
        adopted = table[key][0] if key in table else float('nan')
        say('%-10s %6d %7d %11.4f %11.4f %11.4f %8.2f'
            % (ch or '(plain)', era, len(v), stated / len(v), adopted, rms,
               rms / adopted))

    say()
    say('The lines that carry it')
    say('-----------------------')
    nodata = int(np.sum(fit['leverage'] >= H_MAX))
    say('mean leverage %.3f over %d lines: that fraction of each residual '
        'has already been'
        % (float(np.mean(fit['leverage'])), len(rows)))
    say('absorbed by the fitted level values, and a fit made to '
        'wavelength_calibration_points.csv')
    say('without allowing for it will report a chi-square low by about as '
        'much.  The column')
    say('u_eff_A = u_d_lambda_A / (1 - leverage) is the weight to fit by '
        'instead, and it is')
    say('what the files in wavelength_calibration_fit/ carry.  It is a '
        'weight for an outside')
    say('fit only: this program never uses it, so it cannot reach the '
        'covariance matrix.')
    say('%d line%s reach%s a leverage of %.3f or more - one of their levels '
        'rests on that'
        % (nodata, '' if nodata == 1 else 's',
           'es' if nodata == 1 else '', H_MAX))
    say('line alone, so the fit satisfies them exactly whatever the '
        'calibration is; they')
    say('measure nothing and are left out of those files.')

    # ---- the per-line points --------------------------------------------
    E = {v: sol[k] for k, v in enumerate(levels)}
    K = {c: sol[j] for c, j in kidx.items()}
    points = []
    for r, i in enumerate(rows):
        ln = lines[i]
        k = 1.0 if ln.cls == 'flag' else K['c' if ln.cls == 'c'
                                           else 'plain_%d' % ln.era]
        ritz = E[ln.upp] - E[ln.low] + k * hfs_kappa.A_SCALE * ln.D
        scale = ln.wn * ln.wn * 1e-8
        g = binof[i]
        points.append(dict(
            wn_obs='%.4f' % ln.wn, lambda_A='%.4f' % lam[i],
            block=blk[i],
            group=('%.0f-%.0f' % span[g]) if g in span else '',
            char=ln.char, era=ln.era, low_id=ln.low, upp_id=ln.upp,
            hfs_shift_cm1='%+.4f' % (k * hfs_kappa.A_SCALE * ln.D),
            ritz_cm1='%.4f' % ritz,
            d_lambda_A='%+.5f' % ((ritz - ln.wn) / scale),
            u_d_lambda_A='%.5f' % (sigma[i] / scale),
            u_line_cm1='%.4f' % sigma[i],
            leverage='%.4f' % fit['leverage'][r],
            u_eff_A=('' if fit['leverage'][r] >= H_MAX
                     else '%.5f' % (sigma[i] / scale
                                    / (1.0 - fit['leverage'][r])))))

    # ---- the correction of every observed wavenumber --------------------
    # Every line of the list, not only the ones the fit rests on: an
    # unclassified line needs the correction too, for dlv.dat and for the
    # next round of identification.  With a polynomial this costs nothing -
    # it is defined everywhere inside its block - where a staircase has it
    # only where a group happens to fall.
    infit = set(rows)
    mid_of = {g: 0.5 * (span[g][0] + span[g][1]) for g in good}
    in_block = collections.defaultdict(list)
    for g in good:
        in_block[g[0]].append(g)

    # A line the fit never saw - an unclassified one, or one in a 25 A bin
    # that holds no accepted line - has no group of its own, but it was
    # recorded on the same plate as the groups around it.  It is given the
    # nearest group of its own block, and none at all if that block carries
    # no parameter.
    fitted = {round(lines[i].wn, 6) for i in infit}

    def nearest_group(b, lam_value):
        here = in_block.get(b)
        if not here:
            return None
        return min(here, key=lambda g: abs(mid_of[g] - lam_value))

    corr = []
    for wn, char in observed_wavenumbers():
        lam_value = 1e8 / wn
        b = block_of(blocks, lam_value)
        scale = wn * wn * 1e-8
        g = nearest_group(b, lam_value)
        row = dict(wn_obs='%.4f' % wn, lambda_A='%.4f' % lam_value,
                   block=b, char=char, era=hfs_kappa.era_of(wn),
                   group=('%.0f-%.0f' % span[g]) if g in span else '',
                   in_fit=1 if round(wn, 6) in fitted else 0)
        if (b in degrees) if args.model == 'poly' else (g is not None):
            v, u = value(cvec(b, g, lam_value))
            row.update(d_lambda_A='%+.5f' % v, u_d_lambda_A='%.5f' % u,
                       own_correction='%+.4f' % (-v * scale),
                       u_own_correction='%.4f' % (u * scale))
        elif b in thin_u:
            # The block was tried and its constant fixed at zero because it
            # was smaller than its own uncertainty.  Zero is then the
            # correction, and the uncertainty the trial measured is what the
            # line carries - it is not zero, and dlv.dat needs it.
            u = thin_u[b][1]
            row.update(d_lambda_A='+0.00000', u_d_lambda_A='%.5f' % u,
                       own_correction='+0.0000',
                       u_own_correction='%.4f' % (u * scale))
        else:
            row.update(d_lambda_A='', u_d_lambda_A='',
                       own_correction='', u_own_correction='')
        corr.append(row)

    if args.no_write:
        print(out.getvalue(), end='')
        return 0

    with open(OUT_REPORT, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(out.getvalue())
    with open(OUT_CURVE, 'w', encoding='utf-8', newline='') as fh:
        w2 = csv.writer(fh, lineterminator='\n')
        w2.writerow(['block', 'lambda_lo_A', 'lambda_hi_A', 'lambda_mid_A',
                     'n_lines', 'd_lambda_A', 'u_d_lambda_A',
                     'shift_cm-1', 'u_shift_cm-1'])
        for c in curve:
            w2.writerow([c[0], '%.3f' % c[1], '%.3f' % c[2], '%.3f' % c[3],
                         c[4], '%+.5f' % c[5], '%.5f' % c[6],
                         '%+.4f' % c[7], '%.4f' % c[8]])
    with open(OUT_POINTS, 'w', encoding='utf-8', newline='') as fh:
        w2 = csv.DictWriter(fh, list(points[0]), lineterminator='\n')
        w2.writeheader()
        w2.writerows(points)
    with open(OUT_CORR, 'w', encoding='utf-8', newline='') as fh:
        w2 = csv.DictWriter(fh, list(corr[0]), lineterminator='\n')
        w2.writeheader()
        w2.writerows(corr)

    # ---- the parameters themselves, with their covariance ---------------
    if args.model == 'poly':
        keys = poly_keys(degrees)
        vec = {}
        for k, kk in enumerate(keys):
            c = np.zeros(ncol)
            c[off + k] = 1.0
            vec[kk] = c
        with open(OUT_POLY, 'w', encoding='utf-8', newline='') as fh:
            w2 = csv.writer(fh, lineterminator='\n')
            w2.writerow(['block', 'degree', 'lambda_lo_A', 'lambda_hi_A',
                         'coefficient_A', 'u_coefficient_A']
                        + ['cov_%d_%d' % kk for kk in keys])
            for kk in keys:
                c = vec[kk]
                v, u = value(c)
                w2.writerow([kk[0], kk[1], '%.3f' % ranges[kk[0]][0],
                             '%.3f' % ranges[kk[0]][1],
                             '%+.6e' % v, '%.6e' % u]
                            + ['%+.6e' % float(c.dot(cov).dot(vec[j]))
                               for j in keys])

    if not os.path.isdir(OUT_FIT):
        os.mkdir(OUT_FIT)
    made = []
    for b in sorted(set(p['block'] for p in points)):
        sel = [p for p in points if p['block'] == b and p['u_eff_A']]
        sel.sort(key=lambda p: float(p['lambda_A']))
        name = os.path.join(OUT_FIT, 'block_%02d.txt' % b)
        with open(name, 'w', encoding='utf-8', newline='\n') as fh:
            for p in sel:
                fh.write('%s\t%s\t%s\n'
                         % (p['lambda_A'], p['d_lambda_A'], p['u_eff_A']))
            fh.write('\n')
            fh.write('c0\t0.0\tvary\n')
            fh.write('c1\t0.0\tvary\n')
        made.append((b, len(sel)))

    print(out.getvalue(), end='')
    print('written: %s' % os.path.basename(OUT_REPORT))
    print('written: %s (%d groups)'
          % (os.path.basename(OUT_CURVE), len(curve)))
    print('written: %s (%d lines)'
          % (os.path.basename(OUT_POINTS), len(points)))
    print('written: %s (%d lines, %d corrected)'
          % (os.path.basename(OUT_CORR), len(corr),
             len([c for c in corr if c['own_correction']])))
    if args.model == 'poly':
        print('written: %s (%d coefficients)'
              % (os.path.basename(OUT_POLY), ncal))
    print('written: %s/ (%d blocks, %d points)'
          % (os.path.basename(OUT_FIT), len(made),
             sum(n for _, n in made)))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
