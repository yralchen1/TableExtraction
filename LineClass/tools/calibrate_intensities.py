#!/usr/bin/env python
"""Find the intensity-scale correction of an observed line list, automatically.

This is a standalone tool.  It imports nothing from the rest of the project:
only the standard library, numpy, openpyxl (to read .xlsx workbooks) and, if
--plot is used, matplotlib.  Every file name, sheet name and column name it
needs is a command-line option.


1.  The problem
---------------
An old line list gives, for every observed spectral line, an "intensity" that
is only a relative number: it was read off a photographic plate, and the
sensitivity of the plate, of the grating and of the optics all vary with
wavelength.  Two lines of equal true brightness recorded in different parts of
the spectrum, or on different plates, therefore carry different numbers.  In
Sugar's Pr III list the numbers run from 1 to about 9000, and the spectrum was
covered in several pieces, so the scale changes SMOOTHLY inside each piece and
JUMPS at the joins between pieces.

The correction we want is a function of wavelength, written in logarithms,

    Icor(lambda)  =  Iobs * exp( P(lambda) ) ,

where

    lambda = the vacuum wavelength of the line in angstroms (1 A = 1e-8 cm),
             obtained from the wavenumber wn in cm^-1 as lambda = 1e8/wn,
    Iobs   = the intensity as originally reported,
    Icor   = the intensity on one uniform scale, valid across the whole
             spectrum,
    P      = a PIECEWISE polynomial: the wavelength range is cut into a few
             contiguous regions, and inside each region P is an ordinary
             polynomial of low degree.  The cuts are where the scale jumps.

Both the positions of the cuts and the polynomials are unknown, and this tool
determines both from the data.


2.  What tells us the correction
--------------------------------
The reference against which the observed intensities are judged is the
CALCULATED intensity of the same transition.  In a plasma in local
thermodynamic equilibrium the number of atoms sitting in the upper level of a
transition falls off with the energy of that level as exp(-Eup/kT), where

    Eup = energy of the upper level of the transition, in cm^-1,
    kT  = an effective excitation temperature, written as an energy in cm^-1
          (kT in cm^-1 times 1.2398e-4 = kT in eV).

The number of photons emitted per second in the line is that population times
gA, the statistical weight of the upper level times the transition probability
per second (a quantity computed by an atomic-structure code).  What a
photographic plate or a photomultiplier records, however, is not the number of
photons but the ENERGY they carry: the reported intensity is the energy flux
under the line contour.  One photon of wavenumber rwn carries an energy
proportional to rwn, so the calculated intensity of a line is

    Icalc  =  C * gA * (rwn/1e8) * exp(-Eup/kT) ,                        (*)

with rwn the Ritz wavenumber of the transition in cm^-1 (the "Ritz" wavenumber
is the one computed from the two level energies, which is more precise than
the measured one) and C a single constant that absorbs everything that does
not vary from line to line - the number of emitting atoms, the solid angle,
the exposure.  The factor rwn/1e8 is 1/lambda; writing it that way keeps the
formula proportional to the wavenumber, as an energy flux must be.

    THE ENERGY-FLUX CONVENTION IS USED EVERYWHERE IN THIS TOOL.  If some
    particular line list did report photon numbers instead, its intensities
    would differ from energy fluxes by a factor lambda, and a rescaling of
    intensities by a smooth function of lambda is exactly what P(lambda)
    absorbs - so nothing here has to change to accommodate that rare case.

Taking the logarithm of (*) turns it into a straight line, the Boltzmann plot:

    y  =  ln( Iobs * 1e8 / (gA * rwn) )  =  ln C  -  Eup/kT .

Every identified line - one whose two levels are known, so that Eup, rwn and
gA are known - gives one point (Eup, y).  Fitting a straight line through them
gives kT from the slope (kT = -1/slope) and C from the intercept
(C = exp(intercept)).

The quantity that drives everything below is the disagreement of one line with
that fit,

    dlnI  =  ln( Icalc / Iobs ) ,

the logarithm of the ratio of the calculated to the observed intensity.  If
the intensity scale were already uniform, dlnI would scatter about zero with
no systematic dependence on wavelength.  It does not: it drifts and jumps.
That drift, plotted against lambda, IS the correction we are after, so the
correction polynomial is fitted to dlnI, and after applying it the Boltzmann
fit is repeated.  Because the fit itself moves when the intensities change,
the two steps have to be repeated until they stop moving.


3.  The sequence this tool performs
-----------------------------------
  step 1  Restore the original intensities, if asked (--restore-from): if the
          line list on disk already carries corrected intensities made with an
          earlier version of the piecewise function, dividing them by
          exp(P_old(lambda)) gives the original numbers back.
  step 2  Boltzmann fit of the current intensities -> C, kT, dlnI per line.
  step 3  Cut the wavelength range into regions and fit a low-degree
          polynomial in each; this is P(lambda).  See section 4.
  step 4  Icor = Iobs*exp(P); go back to step 2.  Stop when no corrected
          intensity moves by more than --tol between two rounds.
  step 5  Remove outliers, in three deliberate stages, repeating steps 2-4
          after each.  See section 5.  Nothing is ever removed by a blunt cut
          on |dlnI| alone.
  step 6  Write the region boundaries and polynomial coefficients out in the
          same plain-text format the file is read in, one region per line:

              lambda_lo  lambda_hi  c0;c1;c2;...

          meaning P(lambda) = c0 + c1*lambda + c2*lambda^2 + ... inside
          [lambda_lo, lambda_hi].  The coefficients are in ASCENDING order and
          the number of them differs from region to region - a reader must
          take the degree from how many are written, not assume one.  The
          regions cover the spectrum contiguously except at the coverage gaps
          (--gap), which no region covers; a wavelength outside every region
          is given the value of the nearest region at that region's own end,
          never an extrapolation.  The terms of a high-degree polynomial in
          lambda ~ 1e4 are large and cancel, so all the digits written matter.


4.  How the regions are found
-----------------------------
Recursive binary splitting.  A region is described by the polynomial degree
that fits it best; "best" is decided by the Bayesian information criterion

    BIC  =  n*ln(SSE/n)  +  k*ln(n) ,

with n the number of points in the region, SSE the sum of squared residuals
and k the number of fitted coefficients.  BIC rewards a smaller residual and
charges for every coefficient, so it will not add a coefficient that buys
nothing.  Then every possible place to cut the region in two is tried, and the
cut is kept if it lowers the total BIC of the two halves by more than one
further penalty of --split-penalty*ln(n) (a breakpoint is itself a fitted
quantity and must earn its place), and if the jump of P across the cut is at
least --min-jump.  The two halves are then examined the same way, and so on.

A region is also split, penalty or not, when the polynomial it needs is too
stiff: if BIC still wants a degree above --max-degree (default 5), the region
is broken in two and each half is fitted with a low-degree polynomial instead.
That is the rule "if the required power is too high, break the region into
smaller pieces".

Two kinds of breakpoint must not be confused:

  * the ones just described, FOUND in the data.  They are placed where the
    lines say the scale changes, and where the true scale is continuous the
    fitted polynomials of two neighbouring regions meet at nearly the same
    value, because the lines on both sides of the cut demand it.
  * the ones GIVEN in advance by --gap: the wavelength intervals no exposure
    covered.  Where the spectrum was recorded in pieces that do not overlap,
    no line was ever registered on both plates, so nothing ties the two
    intensity scales together and the correction is genuinely discontinuous
    there, by an amount no data can reveal.  Such a place is a boundary before
    any fitting is attempted: the lines on the two sides are segmented and
    fitted independently, and neither influences the other.  Fitting one
    polynomial across a gap would smear a real, arbitrarily large step over
    the whole neighbourhood and corrupt the corrected intensities on both
    sides of it.  The gap interval itself belongs to no region, since nothing
    could have been observed inside it.

Every region is kept at no fewer than --min-points points, so a handful of
stray lines cannot manufacture a region of their own.

The correction is fitted only on the identified lines, but it has to be applied
to every line that carries an intensity, and some of those lie beyond the ends
of the fitted range.  The outermost regions are therefore stretched to cover
the whole wavelength range of the line list, and their degree is lowered until
the polynomial no longer runs away over the stretched part, where no line holds
it down.

Two different measures of "runs away" are applied there, and both must pass:

  --max-swing   how far P goes outside the band the region's own lines occupy
                (their 1st to 99th percentile of dlnI).
  --max-extrap  how far P MOVES, over the stretched part, away from the value
                it takes at the outermost line that was actually fitted.

The second is the one that bites at the ends of the spectrum.  The scatter of
dlnI is a unit or two in natural logarithms - gA being a calculated quantity -
so a curved fit can dive by two units past its last point and still sit inside
the band that --max-swing compares it with.  That dive is not scatter: it is
the same systematic factor applied to every corrected intensity out there.
--max-extrap measures it directly and lowers the degree until it is small,
which in practice means a straight line over an edge region, as it should be:
where there are no lines there is nothing to justify a curve.


5.  How outliers are removed
----------------------------
A line whose dlnI is large is not automatically wrong, and throwing away
everything beyond a fixed |dlnI| would bias the correction towards whatever
survives.  The removals are therefore made one physically motivated class at a
time, with the whole fit repeated after each:

  stage 1  No removals.  A first, rough correction.
  stage 2  BAD CALCULATED INTENSITIES.  gA comes from a calculation and has
           its own uncertainty, quoted in percent; in logarithms that
           uncertainty is u_ln = ln(1 + u/100), directly comparable with dlnI.
           Among the lines with |dlnI| > the cut of this pass, those are
           dropped for which the calculation is a sufficient explanation of
           the disagreement, namely
           |dlnI| <= --u-sigma * u_ln (default 2 uncertainties).  Such lines
           say nothing about the plate.  A flat threshold on the uncertainty
           itself - "drop everything beyond the cut whose gA is uncertain by
           more than 50 per cent" - is available as --u-rule percent, and both
           counts are printed at every pass whichever rule is in force; in
           Pr III the flat rule is far too permissive, because nine gA values
           in ten are uncertain by more than 50 per cent, so it would throw
           away almost every outlier without asking whether the uncertainty is
           large enough to matter.
           The stage looks again after the fit has moved, up to --max-passes
           times, since a line can become an outlier only once the correction
           has changed shape.
           THE CUT IS STAGED.  --dln-cut takes a list, by default 4 3 2: the
           first pass condemns only lines beyond a factor e^4 = 55, the second
           beyond e^3 = 20, and every later pass - and stages 3 and 4 - beyond
           e^2 = 7.4.  The reason is that on the first pass the correction is
           still the rough stage-1 one, and a line whose only fault is that it
           sits where the correction has not yet found its shape would be
           thrown out on the strength of an error the fit was about to
           remove.  Loosening the early passes lets the shape settle first and
           leaves the tight cut to judge lines against a correction that is
           nearly final.  Give a single value for a fixed cut.
  stage 3  SELF-ABSORBED LINES.  A photon emitted deep in the source can be
           reabsorbed by an atom in the same lower level before it leaves,
           which makes the line look weaker than it is.  This bites hardest
           for lines that end on a heavily populated low level - above all the
           ground level (a "resonance line").  Self-absorption can only make a
           line look too WEAK, i.e. can only make dlnI = ln(Icalc/Iobs) too
           LARGE.  So among the lines still left with dlnI > +--dln-cut, those
           whose lower level lies at or below --self-abs-elow cm^-1 are
           dropped.  Lines that are too STRONG are never dropped by this test,
           because self-absorption cannot explain them.
  stage 4  WHAT IS LEFT.  If, after stages 2 and 3, only a few outliers per
           region remain and every region still has enough points for a
           confident fit, they are all dropped and the whole sequence is run
           once more with the fine convergence threshold --tol-final; this is
           the final answer.  If instead some region would lose too large a
           fraction of its points (--max-drop-frac) or fall below
           --min-points, the tool STOPS and reports the problem rather than
           producing a correction that rests on a handful of survivors.


6.  A caution about what the correction can and cannot separate
---------------------------------------------------------------
The Boltzmann fit describes the intensities with the upper-level energy Eup,
and the correction describes them with the wavelength lambda.  These two are
not independent - a line of high upper energy tends to be a short-wavelength
line - so part of a genuine temperature effect can be absorbed into P(lambda)
and vice versa.  The excitation temperature that comes out of this procedure
is therefore an effective number, not a measurement of the plasma; what the
procedure IS reliable for is the relative intensity scale, which is what it is
used for.


Usage
-----
    python tools/calibrate_intensities.py                       # full run
    python tools/calibrate_intensities.py --restore-only        # step 1 only
    python tools/calibrate_intensities.py --max-degree 4 --plot cal.png

The defaults name the files of this project, but nothing in the code knows
about them.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys

import numpy as np

# ---------------------------------------------------------------------------
# coverage gaps
# ---------------------------------------------------------------------------
# Wavelength intervals, in vacuum angstroms, that the spectrograph never
# recorded: they fall between two photographic exposures whose ranges do not
# overlap.  No line can have been observed inside such an interval, and - what
# matters for the correction - the intensity scales of the two exposures on
# either side were never tied to each other, because no line was recorded on
# both plates.  The correction may therefore JUMP by an arbitrary amount across
# a gap, and a single polynomial must never be fitted through one: the fit
# would smear a real, unknown step over the whole neighbourhood.
#
# The defaults are the two gaps of Sugar's Pr III exposures; they are the same
# numbers as COVERAGE_GAPS_A in level_shifts.py, which uses them to drop
# predicted transitions that could not have been observed.  Use --gap to give
# other ones, --no-gaps for a spectrum recorded in one continuous piece.
DEFAULT_GAPS_A = ((1522.49, 1529.85), (2103.46, 2107.92))


# ---------------------------------------------------------------------------
# reading tables
# ---------------------------------------------------------------------------


def read_table(path, sheet=None):
    """Read a .xlsx worksheet or a .csv file as (header, rows-of-strings).

    Cells are returned exactly as they are stored, except that numbers are
    kept as numbers.  Level identifiers such as '059003.000042' must not lose
    their leading zero, which is why openpyxl is used directly rather than a
    dataframe reader that would guess the type of the column.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.csv', '.txt'):
        with open(path, newline='', encoding='utf-8-sig') as fh:
            rd = list(csv.reader(fh))
        if not rd:
            raise SystemExit('%s is empty' % path)
        return rd[0], rd[1:]
    import openpyxl
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    ws = wb[sheet] if sheet else wb[wb.sheetnames[0]]
    it = ws.iter_rows(values_only=True)
    header = ['' if c is None else str(c).strip() for c in next(it)]
    rows = [list(r) for r in it]
    wb.close()
    return header, rows


def column(header, rows, name, path):
    """Return one column of a table, by the name written in its header row."""
    try:
        j = header.index(name)
    except ValueError:
        raise SystemExit("column '%s' is not in %s; it has: %s"
                         % (name, path, ', '.join(h for h in header if h)))
    return [r[j] if j < len(r) else None for r in rows]


def as_float(v):
    """A cell as a float, or None if it is empty or not a number."""
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return float(v)
    s = str(v).strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def as_id(v):
    """A level identifier as text, with no numeric mangling."""
    if v is None:
        return None
    s = str(v).strip()
    return s or None


# ---------------------------------------------------------------------------
# the piecewise correction: file format, evaluation
# ---------------------------------------------------------------------------


class Region(object):
    """One contiguous wavelength region and its polynomial.

    lo, hi  = the ends of the region in angstroms
    coef    = polynomial coefficients in ASCENDING order, in plain lambda:
              P(lambda) = coef[0] + coef[1]*lambda + coef[2]*lambda^2 + ...
    poly    = the same polynomial as a numpy object fitted on a rescaled
              variable; it is what the code evaluates, because evaluating a
              high-degree polynomial directly in lambda ~ 1e4 loses digits.
    n       = how many points were fitted in it
    rms     = root-mean-square residual of dlnI in it
    """

    def __init__(self, lo, hi, poly, n=0, rms=float('nan')):
        self.lo = float(lo)
        self.hi = float(hi)
        self.poly = poly
        self.coef = np.asarray(poly.convert().coef, float)
        self.n = int(n)
        self.rms = float(rms)

    @property
    def degree(self):
        return len(self.poly.coef) - 1

    def __call__(self, lam):
        return self.poly(lam)


def read_correction(path):
    """Read a piecewise-correction file: 'lo hi c0;c1;c2;...' per line."""
    regions = []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            lo, hi = float(parts[0]), float(parts[1])
            coef = [float(c) for c in parts[2].split(';') if c.strip()]
            poly = np.polynomial.Polynomial(coef)
            regions.append(Region(lo, hi, poly))
    if not regions:
        raise SystemExit('no regions found in %s' % path)
    regions.sort(key=lambda r: r.lo)
    return regions


def write_correction(path, regions):
    """Write the regions in the same plain-text format they are read in."""
    with open(path, 'w', encoding='utf-8') as fh:
        for r in regions:
            coef = ';'.join('%.15g' % c for c in r.coef)
            fh.write('%-10.4f %-11.4f %s\n' % (r.lo, r.hi, coef))


def apply_correction(regions, lam):
    """Evaluate the piecewise polynomial at every wavelength in lam.

    A wavelength outside every region is given the value of the nearest
    region, evaluated at that region's own end - never extrapolated, because a
    polynomial leaves its fitted range explosively.
    """
    lam = np.asarray(lam, float)
    out = np.full(lam.shape, np.nan)
    for r in regions:
        m = (lam >= r.lo) & (lam <= r.hi) & np.isnan(out)
        if m.any():
            out[m] = r(lam[m])
    miss = np.isnan(out)
    if miss.any():
        los = np.array([r.lo for r in regions])
        his = np.array([r.hi for r in regions])
        for i in np.nonzero(miss)[0]:
            x = lam[i]
            j = int(np.argmin(np.minimum(np.abs(los - x), np.abs(his - x))))
            out[i] = regions[j](min(max(x, regions[j].lo), regions[j].hi))
    return out


# ---------------------------------------------------------------------------
# the Boltzmann fit
# ---------------------------------------------------------------------------


def boltzmann_y(intens, gA, rwn):
    """The ordinate of the Boltzmann plot, ln(Iobs*1e8/(gA*rwn))."""
    return np.log(intens * 1e8 / (gA * rwn))


def joint_fit(d, use, regions):
    """Fit the Boltzmann line and all the polynomials in one least squares.

    The model of section 2, with the correction written in, is

        ln(I0 * 1e8/(gA*rwn))  =  ln C  -  Eup/kT  -  P(lambda) ,

    I0 being the original (uncorrected) intensity.  Every unknown in it - ln C,
    1/kT and every polynomial coefficient of every region - enters linearly, so
    ONE least-squares solution gives them all at once.  Fitting them separately
    and going back and forth instead converges very slowly, because the upper
    energy Eup and the wavelength lambda are correlated (section 6): a change
    of temperature can be partly mimicked by a tilt of P, and the two keep
    handing that tilt back to each other.  Solving for them together settles
    the split in one step, in the way least squares always settles correlated
    parameters.

    The constant term of P and ln C are indistinguishable, so the solution is
    the minimum-norm one and is then re-levelled to make P average to zero over
    the fitted lines.

    Returns (C, kT, P, regions), the regions being the same ones with their
    polynomials replaced by the jointly fitted ones.
    """
    y = boltzmann_y(d['I0'], d['gA'], d['rwn'])
    n = len(y)
    cols = [np.ones(n), d['Eup']]
    spans = []
    if regions:
        where = region_of(regions, d['lam'])
        for k, r in enumerate(regions):
            m = (where == k).astype(float)
            t = (2 * d['lam'] - (r.lo + r.hi)) / (r.hi - r.lo)
            j0 = len(cols)
            for j in range(r.degree + 1):
                cols.append(-m * t ** j)
            spans.append((j0, len(cols)))
    A = np.vstack(cols).T
    coef, _, _, _ = np.linalg.lstsq(A[use], y[use], rcond=None)
    C = math.exp(coef[0])
    kT = -1.0 / coef[1] if coef[1] != 0 else float('inf')
    out = []
    for r, (j0, j1) in zip(regions or [], spans):
        poly = np.polynomial.Polynomial(coef[j0:j1], domain=[r.lo, r.hi],
                                        window=[-1.0, 1.0])
        out.append(Region(r.lo, r.hi, poly, r.n, r.rms))
    P = apply_correction(out, d['lam']) if out else np.zeros(n)
    shift = float(np.mean(P[use])) if out else 0.0
    if shift:
        for r in out:
            r.poly = r.poly - shift
            r.coef = np.asarray(r.poly.convert().coef, float)
        P = P - shift
        C = C * math.exp(-shift)
    return C, kT, P, out


def calc_intensity(C, kT, gA, rwn, Eup):
    """The calculated intensity, energy-flux convention: eq. (*) above."""
    return C * gA * (rwn / 1e8) * np.exp(-Eup / kT)


# ---------------------------------------------------------------------------
# fitting one region, choosing its degree
# ---------------------------------------------------------------------------


def _fit(x, y, deg):
    """Least-squares polynomial of the given degree; returns (poly, SSE)."""
    poly = np.polynomial.Polynomial.fit(x, y, deg)
    resid = y - poly(x)
    return poly, float(resid @ resid)


def _bic(n, sse, k):
    """Bayesian information criterion of a fit with k free coefficients."""
    return n * math.log(max(sse, 1e-300) / n) + k * math.log(n)


def swing(poly, x, y):
    """How far a fitted polynomial goes beyond what its own points ask for.

    A least-squares polynomial is held down only where there are points; near
    the ends of a region, where they thin out, a high-degree fit can shoot off.
    Since the correction multiplies the intensities by exp(P), such an
    excursion would be ruinous, so it is measured and refused.  The comparison
    is with the 1st and 99th percentile of the points, not their extremes, so
    that one stray line does not licence an excursion.
    """
    g = np.linspace(x.min(), x.max(), 200)
    p = poly(g)
    lo, hi = np.percentile(y, 1), np.percentile(y, 99)
    return max(float(p.max()) - hi, lo - float(p.min()), 0.0)


def choose_degree(x, y, max_degree, pts_per_coef, max_swing=float('inf')):
    """Best polynomial for one region.

    Returns (poly, sse, bic, wants_more), where wants_more is True when the
    criterion would still like a degree above max_degree - the signal that the
    region should be broken into smaller pieces instead.  A degree whose
    polynomial runs away (see swing) is not eligible, however good its BIC.
    """
    n = len(x)
    cap = min(max_degree, max(0, n // max(pts_per_coef, 1) - 1))
    best = None
    for d in range(cap + 1):
        poly, sse = _fit(x, y, d)
        b = _bic(n, sse, d + 1)
        if d and swing(poly, x, y) > max_swing:
            continue
        if best is None or b < best[2]:
            best = (poly, sse, b, d)
    if best is None:
        poly, sse = _fit(x, y, 0)
        best = (poly, sse, _bic(n, sse, 1), 0)
    wants_more = False
    if cap == max_degree:
        for d in (max_degree + 1, max_degree + 2):
            if n // max(pts_per_coef, 1) - 1 < d:
                break
            _, sse = _fit(x, y, d)
            if _bic(n, sse, d + 1) < best[2]:
                wants_more = True
                break
    return best[0], best[1], best[2], wants_more


def _split_positions(x, min_points, max_tries):
    """Indices at which a region may be cut, at most max_tries of them."""
    n = len(x)
    first, last = min_points, n - min_points
    if last <= first:
        return []
    idx = list(range(first, last + 1))
    # a cut may only be made between two different wavelengths
    idx = [i for i in idx if x[i] > x[i - 1]]
    if len(idx) > max_tries:
        step = len(idx) / float(max_tries)
        idx = [idx[int(k * step)] for k in range(max_tries)]
    return idx


def best_split(x, y, opt):
    """The cut that lowers the total BIC of the two halves the most.

    Returns (index, bic_total, jump) or None.  jump is the size of the step of
    the fitted correction across the cut, in natural logarithms.
    """
    cand = _split_positions(x, opt.min_points, opt.max_split_tries)
    best = None
    for i in cand:
        xl, yl, xr, yr = x[:i], y[:i], x[i:], y[i:]
        try:
            pl, sl = _fit(xl, yl, min(opt.split_degree, len(xl) // opt.pts_per_coef - 1))
            pr, sr = _fit(xr, yr, min(opt.split_degree, len(xr) // opt.pts_per_coef - 1))
        except Exception:
            continue
        b = (_bic(len(xl), sl, opt.split_degree + 1)
             + _bic(len(xr), sr, opt.split_degree + 1))
        if best is None or b < best[1]:
            cut = 0.5 * (x[i - 1] + x[i])
            best = (i, b, abs(float(pr(cut) - pl(cut))))
    return best


def segment(x, y, opt, depth=0):
    """Cut [x] into regions and fit each; returns a list of (i0, i1, poly).

    x must be sorted.  Indices are into the arrays passed in.
    """
    n = len(x)
    poly, sse, bic, wants_more = choose_degree(
        x, y, opt.max_degree, opt.pts_per_coef, opt.max_swing)
    if n < 2 * opt.min_points or depth >= opt.max_depth:
        return [(0, n, poly)]
    sp = best_split(x, y, opt)
    if sp is None:
        return [(0, n, poly)]
    i, bic_split, jump = sp
    penalty = opt.split_penalty * math.log(n)
    worth_it = bic_split + penalty < bic and jump >= opt.min_jump
    if not (worth_it or wants_more):
        return [(0, n, poly)]
    left = segment(x[:i], y[:i], opt, depth + 1)
    right = segment(x[i:], y[i:], opt, depth + 1)
    return left + [(a + i, b + i, p) for (a, b, p) in right]


def gap_blocks(x, gaps, log=None):
    """Cut the sorted wavelengths x into the blocks separated by the gaps.

    Returns a list of (i0, i1, lo_edge, hi_edge): the lines x[i0:i1] of one
    block, and the wavelengths its outermost regions must reach - the near
    side of the gap that bounds the block, or None where the block is bounded
    by the end of the spectrum instead and the caller decides how far to go.
    Lines that fall inside a gap (there should be none) are reported and left
    out of every block; they still get corrected, by the clamping rule of
    apply_correction.
    """
    blocks = []
    start = 0
    lo_edge = None
    for glo, ghi in gaps:
        i = int(np.searchsorted(x, glo, side='left'))
        j = int(np.searchsorted(x, ghi, side='right'))
        if j > i and log is not None:
            log('    WARNING: %d fitted lines fall inside the coverage gap '
                '%.2f - %.2f A, where nothing could have been observed; they '
                'are left out of the region fits' % (j - i, glo, ghi))
        if i > start:
            blocks.append((start, i, lo_edge, glo))
        elif log is not None and lo_edge is not None:
            log('    note: no lines between the gaps ending at %.2f A and '
                'starting at %.2f A' % (lo_edge, glo))
        lo_edge = ghi
        start = max(start, j)
    if start < len(x):
        blocks.append((start, len(x), lo_edge, None))
    return blocks


def _extrapolation(region, x):
    """How far P moves over the parts of a region beyond its outermost lines.

    A region stretched to a gap edge or to the end of the spectrum is
    evaluated where no line holds the polynomial down.  What matters there is
    not how far P goes beyond the SCATTER of the region (that is swing, and a
    region whose lines scatter by two units in ln can hide a two-unit
    excursion inside it) but how far P MOVES away from the last value the
    lines actually determined.  That distance is returned, in natural
    logarithms; it is 0 for a region that is not stretched at all.
    """
    lo_fit, hi_fit = float(np.min(x)), float(np.max(x))
    worst = 0.0
    for edge, g in ((lo_fit, np.linspace(region.lo, lo_fit, 200)),
                    (hi_fit, np.linspace(hi_fit, region.hi, 200))):
        if g[-1] - g[0] <= 0:
            continue
        ref = float(region.poly(edge))
        worst = max(worst, float(np.abs(region.poly(g) - ref).max()))
    return worst


def _tame(region, x, y, max_swing, max_extrap=float('inf')):
    """Lower the degree of a stretched region until it stops running away.

    A region whose end was pushed out to a gap edge or to the end of the
    spectrum is fitted over a wider interval than its lines cover, and a
    polynomial is held down only where there are points.  The degree is
    lowered until BOTH of the following hold over the whole stretched
    interval:

      * the polynomial stays within max_swing of what its own lines ask for
        (the 1st to 99th percentile band of their dlnI), and
      * it moves by no more than max_extrap away from the value it takes at
        the outermost line of the region (see _extrapolation).

    The second test is the one that matters at the ends of the spectrum.  The
    first compares the polynomial with the SCATTER of the region, and where
    that scatter is a unit or two in ln - which it is, gA being what it is -
    an excursion of the same size passes it unnoticed even though it is a
    systematic error in every corrected intensity out there, not scatter.
    """
    g = np.linspace(region.lo, region.hi, 400)
    while region.degree > 0:
        p = region.poly(g)
        within = max(p.max() - np.percentile(y, 99),
                     np.percentile(y, 1) - p.min()) <= max_swing
        if within and _extrapolation(region, x) <= max_extrap:
            break
        region.poly = _fit(x, y, region.degree - 1)[0]
        region.coef = np.asarray(region.poly.convert().coef, float)


def build_regions(lam, dln, opt, full_range=None, log=None):
    """Segment the wavelength axis and return a list of Region objects.

    lam, dln are the lines the correction is fitted to.  full_range, when
    given, is the wavelength interval the finished correction must cover -
    wider than the fitted lines, because lines thrown out of the fit still
    have to be corrected.  The outermost regions are stretched to reach it,
    and their degree is lowered until they no longer run away over the
    stretched interval, where nothing holds them down.

    The coverage gaps in opt.gaps are respected absolutely: they are cuts
    before any fitting is attempted, so no region ever spans one, and the two
    regions that meet at a gap are fitted independently - the step between
    them is whatever the lines on the two sides ask for, however large.
    Inside a block between two gaps the cuts are found from the data, as
    described in section 4 of the module docstring.  The gaps themselves are
    covered by no region: nothing was observed there, and a wavelength asked
    for inside one is given the value at the near edge of the neighbouring
    region (apply_correction).
    """
    order = np.argsort(lam)
    x, y = lam[order], dln[order]
    lo_end, hi_end = full_range if full_range else (x[0], x[-1])
    gaps = getattr(opt, 'gaps', ()) or ()
    regions = []
    for i0, i1, lo_edge, hi_edge in gap_blocks(x, gaps, log):
        xb, yb = x[i0:i1], y[i0:i1]
        pieces = segment(xb, yb, opt)
        first = len(regions)
        for k, (a, b, poly) in enumerate(pieces):
            lo = xb[a] if k == 0 else 0.5 * (xb[a - 1] + xb[a])
            hi = xb[b - 1] if k == len(pieces) - 1 else 0.5 * (xb[b - 1] + xb[b])
            resid = yb[a:b] - poly(xb[a:b])
            regions.append(Region(lo, hi, poly, b - a,
                                  math.sqrt(float(resid @ resid) / (b - a))))
        # the ends of the block: a gap edge, or the end of the spectrum
        r_first, r_last = regions[first], regions[-1]
        r_first.lo = (float(lo_edge) if lo_edge is not None
                      else min(r_first.lo, float(lo_end)) - 1e-6)
        r_last.hi = (float(hi_edge) if hi_edge is not None
                     else max(r_last.hi, float(hi_end)) + 1e-6)
        for r, (a, b) in ((r_first, pieces[0][:2]), (r_last, pieces[-1][:2])):
            _tame(r, xb[a:b], yb[a:b], opt.max_swing,
                  getattr(opt, 'max_extrap', float('inf')))
    if not regions:
        raise SystemExit('no lines left to fit the correction to')
    return regions


# ---------------------------------------------------------------------------
# the iteration: Boltzmann fit <-> correction
# ---------------------------------------------------------------------------


class Fit(object):
    """Everything one converged calibration produced."""

    def __init__(self, regions, C, kT, P, dln, rounds, converged, dev):
        self.regions = regions
        self.C = C
        self.kT = kT
        self.P = P            # the correction at every line
        self.dln = dln        # ln(Icalc/Icor) at every line, after the fit
        self.rounds = rounds
        self.converged = converged
        self.dev = dev


def calibrate(d, use, opt, tol, log):
    """Iterate 'Boltzmann fit -> correction polynomials' to a fixed point.

    d    = the data (a dict of equal-length arrays; see load_data)
    use  = boolean mask, True for the lines allowed into the fits
    tol  = stop when no corrected intensity moves by more than this fraction

    Returns a Fit.  The correction is always defined relative to the ORIGINAL
    intensities, so it never has to be composed with itself between rounds.
    """
    P = np.zeros(len(d['lam']))
    regions = []
    dev = float('inf')
    converged = False
    for rnd in range(1, opt.max_rounds + 1):
        C, kT, Pnew, regions = joint_fit(d, use, regions)
        if not np.isfinite(kT) or kT <= 0:
            raise SystemExit('the Boltzmann fit gave kT = %s; the input '
                             'intensities or gA values cannot be right' % kT)
        dev = float(np.max(np.abs(np.expm1(Pnew - P))))
        P = Pnew
        # dlnI measured against the ORIGINAL intensities: this is what the
        # correction has to reproduce; what is left of it after P is taken out
        # is the disagreement of the corrected intensity with the fit.
        Icalc = calc_intensity(C, kT, d['gA'], d['rwn'], d['Eup'])
        target = np.log(Icalc / d['I0'])
        dln = target - P
        log('    round %2d: %2d regions, C = %.4g, kT = %.1f cm^-1, '
            'rms dlnI = %.4f, max |dIcor| = %s'
            % (rnd, len(regions), C, kT, float(np.std(dln[use])),
               '-' if rnd == 1 else '%.4g' % dev))
        if dev < tol and rnd > 1:
            converged = True
            break
        # where the correction should be cut, and what degree each piece
        # needs, is decided anew from the current disagreement
        regions = build_regions(d['lam'][use], target[use], opt,
                                d.get('lam_cover',
                                      (d['lam'].min(), d['lam'].max())),
                                log if rnd == 1 else None)
    where = region_of(regions, d['lam'])
    for k, r in enumerate(regions):
        m = use & (where == k)
        r.n = int(m.sum())
        r.rms = float(np.sqrt(np.mean(dln[m] ** 2))) if r.n else float('nan')
    return Fit(regions, C, kT, P, dln, rnd, converged, dev)


# ---------------------------------------------------------------------------
# loading the observed lines and the calculated transitions
# ---------------------------------------------------------------------------


def load_data(opt, log):
    """Join the observed lines to the calculated transitions.

    Produces arrays, one entry per identified line that has a calculated
    transition:
        wn    observed wavenumber, cm^-1
        lam   vacuum wavelength, angstroms (1e8/wn)
        Iraw  the intensity as it stands in the line list
        I0    the original intensity (Iraw, or Iraw with an earlier
              correction divided out - see --restore-from)
        gA    statistical weight x transition probability, s^-1
        u     uncertainty of gA, percent
        rwn   Ritz wavenumber of the transition, cm^-1
        Eup   energy of the upper level, cm^-1
        Elow  energy of the lower level, cm^-1 (= Eup - rwn)
    """
    lh, lr = read_table(opt.lines, opt.lines_sheet)
    wn_c = [as_float(v) for v in column(lh, lr, opt.col_wn, opt.lines)]
    in_c = [as_float(v) for v in column(lh, lr, opt.col_intensity, opt.lines)]
    i1_c = [as_id(v) for v in column(lh, lr, opt.col_id1, opt.lines)]
    i2_c = [as_id(v) for v in column(lh, lr, opt.col_id2, opt.lines)]

    th, tr = read_table(opt.transitions, opt.transitions_sheet)
    t1 = [as_id(v) for v in column(th, tr, opt.col_t_id1, opt.transitions)]
    t2 = [as_id(v) for v in column(th, tr, opt.col_t_id2, opt.transitions)]
    trwn = [as_float(v) for v in column(th, tr, opt.col_rwn, opt.transitions)]
    tgA = [as_float(v) for v in column(th, tr, opt.col_gA, opt.transitions)]
    tu = [as_float(v) for v in column(th, tr, opt.col_u, opt.transitions)]
    tEup = [as_float(v) for v in column(th, tr, opt.col_Eup, opt.transitions)]

    trans = {}
    for k in range(len(t1)):
        if t1[k] and t2[k] and tgA[k] and trwn[k] and tEup[k] is not None:
            trans[(t1[k], t2[k])] = (trwn[k], tgA[k], tu[k], tEup[k])

    cols = {n: [] for n in ('wn', 'Iraw', 'gA', 'u', 'rwn', 'Eup')}
    n_lines = n_ided = 0
    cover = []           # every line that will have to be corrected
    for k in range(len(wn_c)):
        if wn_c[k] is None or in_c[k] is None or in_c[k] <= 0:
            continue
        n_lines += 1
        cover.append(1e8 / wn_c[k])
        if not (i1_c[k] and i2_c[k]):
            continue
        n_ided += 1
        t = trans.get((i1_c[k], i2_c[k])) or trans.get((i2_c[k], i1_c[k]))
        if t is None:
            continue
        cols['wn'].append(wn_c[k])
        cols['Iraw'].append(in_c[k])
        cols['rwn'].append(t[0])
        cols['gA'].append(t[1])
        cols['u'].append(t[2] if t[2] is not None else float('nan'))
        cols['Eup'].append(t[3])
    d = {k: np.asarray(v, float) for k, v in cols.items()}
    d['lam'] = 1e8 / d['wn']
    d['Elow'] = d['Eup'] - d['rwn']
    # The finished correction has to cover every line whose intensity will be
    # corrected, not only the identified ones it is fitted to.
    d['lam_cover'] = (min(cover), max(cover))
    if opt.cover:
        d['lam_cover'] = (min(d['lam_cover'][0], opt.cover[0]),
                          max(d['lam_cover'][1], opt.cover[1]))
    log('lines read: %d with an intensity, %d of them identified, '
        '%d with a calculated transition' % (n_lines, n_ided, len(d['wn'])))
    log('    the correction must cover %.2f - %.2f A (all lines with an '
        'intensity); it is fitted on %.2f - %.2f A'
        % (d['lam_cover'][0], d['lam_cover'][1],
           float(d['lam'].min()), float(d['lam'].max())))
    if len(d['wn']) < 50:
        raise SystemExit('too few usable lines (%d) to calibrate anything'
                         % len(d['wn']))
    return d


def restore_originals(d, opt, log):
    """Undo an earlier piecewise correction, recovering the original scale.

    The stored intensity is Iraw = K * I0 * exp(P_old(lambda)), where K is
    whatever overall factor the earlier work carried (a constant in P_old is
    indistinguishable from it).  Dividing by exp(P_old) leaves K*I0; K is then
    taken from --restore-scale, which with the default 'auto' is the power of
    ten nearest to the weakest restored intensity - i.e. the assumption that
    the faintest line of the list was recorded as intensity 1.
    """
    old = read_correction(opt.restore_from)
    log('restoring the original intensities with %d regions from %s'
        % (len(old), opt.restore_from))
    v = d['Iraw'] / np.exp(apply_correction(old, d['lam']))
    if opt.restore_scale == 'auto':
        K = 10.0 ** round(math.log10(float(np.min(v))))
    else:
        K = float(opt.restore_scale)
    v = v / K
    log('    scale factor divided out: %g   ->  restored intensities run '
        'from %.3g to %.3g' % (K, float(v.min()), float(v.max())))
    if opt.restore_round:
        r = np.round(v)
        r[r < 1] = 1.0
        off = np.abs(v / r - 1.0)
        log('    snapped to whole numbers (the original list reports whole '
            'numbers): largest departure %.1f%%, 99th percentile %.1f%%'
            % (100 * off.max(), 100 * np.percentile(off, 99)))
        if off.max() > 0.35:
            log('    WARNING: some restored values are far from a whole '
                'number; --no-restore-round may be safer here')
        v = r
    return v


# ---------------------------------------------------------------------------
# the staged removal of outliers
# ---------------------------------------------------------------------------


def region_of(regions, lam):
    """Index of the region each wavelength falls in (nearest one if outside)."""
    out = np.full(len(lam), -1, int)
    for k, r in enumerate(regions):
        m = (lam >= r.lo) & (lam <= r.hi) & (out < 0)
        out[m] = k
    for i in np.nonzero(out < 0)[0]:
        out[i] = int(np.argmin([min(abs(r.lo - lam[i]), abs(r.hi - lam[i]))
                                for r in regions]))
    return out


def run(opt, log):
    d = load_data(opt, log)
    d['I0'] = restore_originals(d, opt, log) if opt.restore_from else d['Iraw'].copy()
    if opt.restore_only:
        return d, None, None

    n = len(d['I0'])
    use = np.ones(n, bool)
    dropped = np.zeros(n, int)   # 0 = kept; 2,3,4 = the stage that dropped it

    log('')
    log('stage 1 - first correction, nothing removed (%d lines)' % use.sum())
    fit = calibrate(d, use, opt, opt.tol, log)

    cuts = ([float(c) for c in opt.dln_cut]
            if isinstance(opt.dln_cut, (list, tuple)) else [float(opt.dln_cut)])
    final_cut = cuts[-1]

    def cut_of(npass):
        return cuts[min(npass, len(cuts)) - 1]

    log('')
    sched = ('|dlnI| > %s in successive passes' % ', '.join('%.2f' % c for c in cuts)
             if len(cuts) > 1 else '|dlnI| > %.2f' % final_cut)
    if opt.u_rule == 'percent':
        log('stage 2 - dropping lines with %s whose calculated gA '
            'is uncertain by more than %.0f%%' % (sched, opt.u_cut))
    else:
        log('stage 2 - dropping lines with %s whose disagreement '
            'is within %.1f uncertainties of their own calculated gA'
            % (sched, opt.u_sigma))
    for npass in range(1, opt.max_passes + 1):
        cut = cut_of(npass)
        beyond = use & (np.abs(fit.dln) > cut)
        u_ln = np.log1p(d['u'] / 100.0)          # gA uncertainty in logarithms
        by_pct = beyond & (d['u'] > opt.u_cut)
        by_sig = beyond & (np.abs(fit.dln) <= opt.u_sigma * u_ln)
        bad = by_pct if opt.u_rule == 'percent' else by_sig
        log('    pass %d (cut %.2f): %d lines beyond it; the percentage rule '
            'would drop %d of them, the uncertainty rule %d; %d dropped'
            % (npass, cut, int(beyond.sum()), int(by_pct.sum()),
               int(by_sig.sum()), int(bad.sum())))
        if not bad.any() and cut == final_cut:
            break
        if not bad.any():
            continue
        use &= ~bad
        dropped[bad] = 2
        fit = calibrate(d, use, opt, opt.tol, log)

    log('')
    log('stage 3 - dropping apparently too-weak lines that end on a level at '
        'or below %.0f cm^-1 (self-absorption)' % opt.self_abs_elow)
    log('    lines still too weak by more than the cut, by the energy of the '
        'level they end on:')
    weak = use & (fit.dln > final_cut)
    for e in (0.0, 1000.0, 2000.0, 5000.0, 10000.0):
        log('        lower level <= %6.0f cm^-1 : %4d of the %d such lines '
            'are that weak' % (e, int((weak & (d['Elow'] <= e)).sum()),
                               int((use & (d['Elow'] <= e)).sum())))
    for npass in range(1, opt.max_passes + 1):
        sa = use & (fit.dln > final_cut) & (d['Elow'] <= opt.self_abs_elow)
        log('    pass %d: %d lines dropped as self-absorbed' % (npass, int(sa.sum())))
        if not sa.any():
            break
        use &= ~sa
        dropped[sa] = 3
        fit = calibrate(d, use, opt, opt.tol, log)

    log('')
    log('stage 4 - what is left beyond |dlnI| > %.2f' % final_cut)
    out = use & (np.abs(fit.dln) > final_cut)
    reg = region_of(fit.regions, d['lam'])
    safe = True
    for k, r in enumerate(fit.regions):
        inreg = use & (reg == k)
        nout = int((out & (reg == k)).sum())
        nin = int(inreg.sum())
        frac = nout / float(nin) if nin else 1.0
        left = nin - nout
        flag = ''
        if frac > opt.max_drop_frac or left < opt.min_points:
            flag = '   <-- TOO MANY'
            safe = False
        log('    region %d  %8.2f - %8.2f A : %4d lines, %3d outliers '
            '(%4.1f%%), %4d would remain%s'
            % (k + 1, r.lo, r.hi, nin, nout, 100 * frac, left, flag))
    if not safe:
        log('')
        log('STOPPING: at least one region would lose more than %.0f%% of its '
            'lines or fall below %d lines.  The remaining disagreement is not '
            'a matter of a few bad lines, and dropping it would make the '
            'correction rest on whatever happened to survive.  Look at the '
            'diagnostics file, or relax --dln-cut / --max-drop-frac / '
            '--min-points, before trusting a correction from this data.'
            % (100 * opt.max_drop_frac, opt.min_points))
        return d, fit, dict(use=use, dropped=dropped, stopped=True)

    use &= ~out
    dropped[out] = 4
    log('    %d outliers dropped; final iteration to %.3g' % (int(out.sum()), opt.tol_final))
    fit = calibrate(d, use, opt, opt.tol_final, log)
    if not fit.converged:
        log('    WARNING: the final iteration did not converge: the largest '
            'change of a corrected intensity in the last round was %.3g, '
            'against the threshold %.3g' % (fit.dev, opt.tol_final))
    return d, fit, dict(use=use, dropped=dropped, stopped=False)


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------


def write_lines_csv(path, d, fit, state):
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['wn', 'lambda', 'I_reported', 'I_original', 'I_corrected',
                    'P', 'dlnI', 'gA', 'u_pct_gA', 'Eup', 'Elow',
                    'used_in_fit', 'dropped_at_stage'])
        Icor = d['I0'] * np.exp(fit.P)
        for i in range(len(d['wn'])):
            w.writerow(['%.6f' % d['wn'][i], '%.4f' % d['lam'][i],
                        '%.6g' % d['Iraw'][i], '%.6g' % d['I0'][i],
                        '%.6g' % Icor[i], '%.6f' % fit.P[i],
                        '%.6f' % fit.dln[i], '%.6g' % d['gA'][i],
                        '' if np.isnan(d['u'][i]) else '%.2f' % d['u'][i],
                        '%.2f' % d['Eup'][i], '%.2f' % d['Elow'][i],
                        int(state['use'][i]), state['dropped'][i]])


def make_plot(path, d, fit, state, opt):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    use = state['use']
    fig, ax = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    ax[0].plot(d['lam'][~use], np.log(calc_intensity(fit.C, fit.kT, d['gA'][~use],
               d['rwn'][~use], d['Eup'][~use]) / d['I0'][~use]), '.',
               ms=3, color='0.75', label='dropped')
    ax[0].plot(d['lam'][use], (fit.dln + fit.P)[use], '.', ms=3, color='C0',
               label='kept')
    # one curve per region, so that the jumps between them are not bridged
    for k, r in enumerate(fit.regions):
        g = np.linspace(r.lo, r.hi, 400)
        ax[0].plot(g, r(g), '-', color='C3', lw=1.5,
                   label='P(lambda)' if k == 0 else None)
    for r in fit.regions[1:]:
        for a in ax:
            a.axvline(r.lo, color='0.6', lw=0.8, ls='--')
    for lo, hi in (getattr(opt, 'gaps', ()) or ()):
        for a in ax:
            a.axvspan(lo, hi, color='0.85', zorder=0)
    ax[0].set_ylabel('ln(Icalc / Ioriginal)')
    ax[0].legend(markerscale=3, fontsize=8)
    ax[0].set_title('intensity-scale correction: C = %.4g, kT = %.1f cm^-1, '
                    '%d regions' % (fit.C, fit.kT, len(fit.regions)))
    ax[1].plot(d['lam'][use], fit.dln[use], '.', ms=3, color='C2')
    ax[1].axhline(0, color='k', lw=0.8)
    ax[1].set_ylabel('ln(Icalc / Icorrected)')
    ax[1].set_xlabel('vacuum wavelength (A)')
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def report(fit, d, state, opt, log):
    """Print the regions, and check that no polynomial runs away inside one.

    A least-squares polynomial can swing violently near the end of its region
    where there are few points to hold it down; since the correction is a
    multiplication by exp(P), such a swing would ruin the intensities there.
    The check compares the range P covers inside a region with the range the
    lines in that region actually ask for.
    """
    use = state['use']
    where = region_of(fit.regions, d['lam'])
    target = fit.dln + fit.P
    log('')
    log('final correction: %d regions, C = %.6g, kT = %.2f cm^-1 (%.4f eV)'
        % (len(fit.regions), fit.C, fit.kT, fit.kT * 1.23984e-4))
    gaps = getattr(opt, 'gaps', ()) or ()
    if gaps:
        log('  coverage gaps honoured (no region crosses one, and the '
            'correction may jump across it): '
            + ', '.join('%.2f - %.2f A' % (lo, hi) for lo, hi in gaps))
    log('  region       from          to    deg     n   rms dlnI    P(from)    '
        'P(to)   swing')
    for k, r in enumerate(fit.regions):
        for lo, hi in gaps:
            if abs(r.lo - hi) < 1e-6:
                log('         %10.2f  %10.2f      coverage gap: no exposures, '
                    'the scales on the two sides are unrelated' % (lo, hi))
        g = np.linspace(r.lo, r.hi, 400)
        p = r(g)
        m = use & (where == k)
        span = (float(np.percentile(target[m], 1)),
                float(np.percentile(target[m], 99))) if m.any() else (0.0, 0.0)
        swing = max(p.max() - span[1], span[0] - p.min(), 0.0)
        log('  %5d  %10.2f  %10.2f  %5d %5d     %6.3f  %9.3f %8.3f  %6.2f%s'
            % (k + 1, r.lo, r.hi, r.degree, r.n, r.rms, r(r.lo), r(r.hi),
               swing, '   <-- runs away' if swing > opt.max_swing else ''))
    log('  (swing = how far P goes beyond what the lines of the region ask '
        'for, in natural logarithms; a large value at an end of a region '
        'means the polynomial is unconstrained there - lower --max-degree or '
        'raise --pts-per-coef)')


# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------


def parse_args(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    p = argparse.ArgumentParser(
        description='Find the piecewise intensity-scale correction of an '
                    'observed line list from a Boltzmann plot.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    g = p.add_argument_group('input files')
    g.add_argument('--lines', default=os.path.join(root, 'Pr3_lines.xlsx'),
                   help='observed line list (.xlsx or .csv)')
    g.add_argument('--lines-sheet', default='Sheet1')
    g.add_argument('--transitions', default=os.path.join(root, 'Icalc.xlsx'),
                   help='calculated transitions (.xlsx or .csv)')
    g.add_argument('--transitions-sheet', default='Icalc')

    g = p.add_argument_group('column names')
    g.add_argument('--col-wn', default='own', help='observed wavenumber, cm^-1')
    g.add_argument('--col-intensity', default='Icor', help='reported intensity')
    g.add_argument('--col-id1', default='id1', help='lower-level identifier')
    g.add_argument('--col-id2', default='id2', help='upper-level identifier')
    g.add_argument('--col-t-id1', default='id1')
    g.add_argument('--col-t-id2', default='id2')
    g.add_argument('--col-rwn', default='rwn', help='Ritz wavenumber, cm^-1')
    g.add_argument('--col-gA', default='gA', help='g times A, s^-1')
    g.add_argument('--col-u', default='u%gA', help='uncertainty of gA, percent')
    g.add_argument('--col-Eup', default='Eup', help='upper-level energy, cm^-1')

    g = p.add_argument_group('restoring the original intensities')
    g.add_argument('--restore-from', default=None, metavar='FILE',
                   help='piecewise correction already applied to the '
                        'intensity column; divide it out first')
    g.add_argument('--restore-scale', default='auto',
                   help="overall factor to divide out as well; 'auto' takes "
                        'the power of ten nearest the weakest restored line')
    g.add_argument('--no-restore-round', dest='restore_round',
                   action='store_false',
                   help='do not snap the restored intensities to whole numbers')
    g.add_argument('--restore-only', action='store_true',
                   help='restore and report, then stop')

    g = p.add_argument_group('regions and polynomials')
    g.add_argument('--cover', type=float, nargs=2, default=None,
                   metavar=('LO', 'HI'),
                   help='wavelength interval, in angstroms, that the finished '
                        'correction must cover, if wider than the line list '
                        'itself; use the wavelength range of the spectrum the '
                        'rest of the work scans, so that nothing asks the '
                        'correction for a wavelength it does not have')
    g.add_argument('--max-degree', type=int, default=5,
                   help='highest polynomial degree allowed in one region; a '
                        'region that wants more is split instead')
    g.add_argument('--min-points', type=int, default=40,
                   help='fewest lines a region may contain')
    g.add_argument('--pts-per-coef', type=int, default=8,
                   help='fewest lines per fitted coefficient')
    g.add_argument('--split-penalty', type=float, default=2.0,
                   help='how many coefficients a new breakpoint must be worth')
    g.add_argument('--min-jump', type=float, default=0.0,
                   help='smallest step of P across a breakpoint, in natural '
                        'logarithms, for the breakpoint to be kept')
    g.add_argument('--split-degree', type=int, default=3,
                   help='degree used while searching for the best breakpoint')
    g.add_argument('--max-split-tries', type=int, default=150,
                   help='candidate breakpoints examined per region')
    g.add_argument('--max-depth', type=int, default=6,
                   help='deepest recursion of the splitting')
    g.add_argument('--gap', dest='gaps', action='append', nargs=2, type=float,
                   metavar=('LO', 'HI'), default=None,
                   help='a wavelength interval, in angstroms, that no exposure '
                        'covered: a forced region boundary, across which the '
                        'correction is free to jump.  Repeat for several; the '
                        'default is the two gaps of this spectrum, %s'
                        % ', '.join('%.2f-%.2f' % g for g in DEFAULT_GAPS_A))
    g.add_argument('--no-gaps', action='store_true',
                   help='the spectrum was recorded in one continuous piece: '
                        'no forced boundaries')

    g = p.add_argument_group('iteration and outliers')
    g.add_argument('--tol', type=float, default=0.02,
                   help='rough convergence: largest allowed relative change '
                        'of a corrected intensity between two rounds')
    g.add_argument('--tol-final', type=float, default=0.002,
                   help='the same, for the final iteration')
    g.add_argument('--max-rounds', type=int, default=25)
    g.add_argument('--dln-cut', type=float, nargs='+', default=[4.0, 3.0, 2.0],
                   metavar='CUT',
                   help='|ln(Icalc/Iobs)| beyond which a line is an '
                        'outlier.  Several values are a SCHEDULE, one '
                        'per pass of stage 2: the first pass uses the '
                        'first, the second the second, and every later '
                        'pass and stages 3 and 4 use the last.  The '
                        'default 4 3 2 removes only the wildest lines '
                        'while the correction is still rough and '
                        'tightens as it settles, so that a line is not '
                        'condemned by a shape the fit has not found '
                        'yet.  Give one value for a fixed cut.')
    g.add_argument('--u-rule', choices=('percent', 'sigma'), default='percent',
                   help="which outliers are blamed on the calculated gA: "
                        "'sigma' drops those whose |dlnI| is no larger than "
                        "--u-sigma times the uncertainty of their own gA "
                        "(in logarithms, ln(1+u/100)); 'percent' drops those "
                        'whose gA is uncertain by more than --u-cut percent')
    g.add_argument('--u-sigma', type=float, default=2.0,
                   help='how many of its own uncertainties a calculated '
                        'intensity may be away before the calculation stops '
                        'being a sufficient explanation')
    g.add_argument('--u-cut', type=float, default=50.0,
                   help="uncertainty of gA, in percent, for --u-rule percent")
    g.add_argument('--max-passes', type=int, default=6,
                   help='how many times a removal stage may look again after '
                        'the fit has moved')
    g.add_argument('--self-abs-elow', type=float, default=0.0,
                   help='lower-level energy, cm^-1, at or below which a '
                        'too-weak line is taken to be self-absorbed; 0 means '
                        'the ground level only (resonance lines)')
    g.add_argument('--max-extrap', type=float, default=0.5,
                   help='how far P may move, in natural logarithms, '
                        'over the part of a region that reaches '
                        'beyond its outermost line - a gap edge, or '
                        'the end of the spectrum.  Out there nothing '
                        'holds the polynomial down, so a curved fit '
                        'is free to dive or climb; the degree is '
                        'lowered until it does not.  0.5 in ln is a '
                        'factor 1.6 in intensity.  Use inf to switch '
                        'the test off and keep the older behaviour, '
                        'which policed only --max-swing.')
    g.add_argument('--max-swing', type=float, default=1.0,
                   help='how far, in natural logarithms, a fitted polynomial '
                        'may go beyond what the lines of its own region ask '
                        'for before it is reported as running away')
    g.add_argument('--max-drop-frac', type=float, default=0.05,
                   help='largest fraction of a region that stage 4 may drop '
                        'before the tool refuses to go on')

    g = p.add_argument_group('output')
    g.add_argument('--out-functions',
                   default=os.path.join(root, 'intensity_correction_auto.txt'),
                   help='where to write the regions and coefficients')
    g.add_argument('--out-lines',
                   default=os.path.join(root, 'intensity_calibration_lines.csv'),
                   help='per-line diagnostics')
    g.add_argument('--plot', default=None, metavar='PNG')
    g.add_argument('--quiet', action='store_true')

    opt = p.parse_args(argv)
    if opt.min_points < opt.pts_per_coef * (opt.split_degree + 1):
        opt.min_points = opt.pts_per_coef * (opt.split_degree + 1)
    if opt.no_gaps:
        opt.gaps = []
    elif opt.gaps is None:
        opt.gaps = [tuple(g) for g in DEFAULT_GAPS_A]
    else:
        opt.gaps = [tuple(g) for g in opt.gaps]
    opt.gaps.sort()
    last = None
    for lo, hi in opt.gaps:
        if not hi > lo:
            raise SystemExit('--gap %g %g: the second wavelength must be the '
                             'larger one' % (lo, hi))
        if last is not None and lo <= last:
            raise SystemExit('--gap: the gaps overlap at %g A' % lo)
        last = hi
    return opt


def main(argv=None):
    opt = parse_args(argv)
    out = []

    def log(msg):
        out.append(msg)
        if not opt.quiet:
            print(msg)
            sys.stdout.flush()

    d, fit, state = run(opt, log)
    if fit is None:
        return 0
    report(fit, d, state, opt, log)
    if not state['stopped']:
        write_correction(opt.out_functions, fit.regions)
        log('')
        log('correction written to %s' % opt.out_functions)
    write_lines_csv(opt.out_lines, d, fit, state)
    log('per-line diagnostics written to %s' % opt.out_lines)
    if opt.plot:
        make_plot(opt.plot, d, fit, state, opt)
        log('plot written to %s' % opt.plot)
    return 1 if state['stopped'] else 0


if __name__ == '__main__':
    sys.exit(main())
