#!/usr/bin/env python
"""Map the wavelength coverage of a photographic line list.

This is a standalone tool.  It imports nothing from the pipeline except the
table reader of its sibling `calibrate_intensities.py`, and takes every file
and column name as a command-line option.


1.  What the tool is for
------------------------
A predicted transition that theory says should be strong, but which the
observed line list does not contain, is evidence against the identification
that predicted it - but only if the line COULD have been recorded where it
falls.  Three things can stop that, and none of them has anything to do with
the atom:

  * the wavelength fell between two photographic exposures, so no plate was
    looking there at all;
  * a defect of the emulsion - a scratch, a fog patch, a grain flaw - sat on
    top of it;
  * a line of another species (an impurity, or another ion of the same
    element) sat on top of it.

For this spectrum we have no record of any of the three.  The plates were
hundreds of separate exposures with varying settings, and their boundaries
were never published; the defects were never catalogued; the impurity spectra
were never listed.  Any list of "known gaps" is therefore guesswork, and this
project's current list - two intervals, hard-coded - has to be treated as a
guess that may be wrong or incomplete.

The point of this tool is that we do not need the record.  All three causes
have the same observable consequence: WHEREVER THEY ACT, THE LINE LIST IS
EMPTY.  Not just empty of the one transition we care about - empty of every
line, of every species, that would otherwise have been recorded there.  So the
line list maps its own blind spots.  A stretch of wavelength that holds far
fewer recorded lines than its neighbours held is a stretch where something was
in the way, and we do not have to know which of the three it was in order to
stop counting a missing line there as evidence.

What the tool produces is a COVERAGE FUNCTION

    c(lambda)  =  the probability that a line which the spectrum really
                  contains, at vacuum wavelength lambda, and which is bright
                  enough to be recorded, did in fact get recorded,
                  RELATIVE to how well that is done in a typical part of
                  the spectrum.

c = 1 means "as well recorded as anywhere else"; c = 0 means "nothing there
could have been seen".  It is deliberately a RELATIVE measure.  The absolute
part - a line's brightness against the detection threshold of the plate -
is a separate question, already answered by the noise level N(lambda) that
`estimate_snr.py` fits.  The two multiply, and each covers what the other
cannot: N(lambda) is a smooth function of wavelength and describes how faint
a line may be before it is lost; c(lambda) is sharp and local and describes
where a line of ANY brightness is lost.


2.  Wavelength, and the coordinate the work is done in
------------------------------------------------------
lambda is the VACUUM wavelength in angstroms (1 A = 1e-8 cm), obtained from a
line's wavenumber wn in cm^-1 as lambda = 1e8 / wn.

All of the fitting is done not in lambda but in

    u  =  ln lambda ,

its natural logarithm.  The reason is dynamic range: this list runs from 822 A
to 10720 A, a factor of thirteen, and the density of recorded lines falls by
more than a factor of ten across it.  A bin of fixed width in lambda would
hold a hundred lines at one end of the list and none at the other.  A bin of
fixed width in u holds a roughly comparable number everywhere, because a fixed
step in u is a fixed step in RELATIVE wavelength: du = dlambda / lambda.  The
bin width --bin is therefore given as a fraction: --bin 0.0002 means each bin
spans 0.02 % in wavelength, which is 0.16 A at 822 A and 2.1 A at 10700 A.


3.  The model
-------------
Every bin b of the grid is in one of two hidden states:

    covered   the plates were looking, and nothing was in the way;
    blocked   they were not, or something was.

Let

    n_b       the number of lines the list records in bin b (an integer);
    rho(u)    the BASELINE: the density of recordable lines per unit u that
              a covered stretch of spectrum at that wavelength shows.  This
              is a property of the atom and of the plate sensitivity, and it
              varies smoothly and slowly with wavelength;
    R_b       = rho(u_b) * (bin width in u), the number of lines bin b would
              hold if it were covered;
    r0        the leak-through: the fraction of R_b that still gets recorded
              in a BLOCKED bin.  Not zero, because obscuration is rarely
              total - a scratch may hide part of a bin, an impurity line
              hides only what it overlaps - and because leaving it at exactly
              zero would let a single stray line veto a hole that is otherwise
              unmistakable.

The number of lines in a bin is Poisson: a bin with an expected count m holds
n lines with probability exp(-m) m^n / n!.  So

    covered bin:   n_b ~ Poisson(R_b)
    blocked bin:   n_b ~ Poisson(r0 * R_b)

and the states are not independent from bin to bin - obscuration comes in
runs, a plate gap or a scratch being many bins wide.  That is modelled with
two transition probabilities,

    a  =  P(the next bin is blocked  | this one is covered)
    b  =  P(the next bin is covered  | this one is blocked)

so that a run of blocked bins lasts 1/b bins on average, and the fraction of
the spectrum that is blocked is a/(a+b).  This is a hidden Markov model, and
the standard forward-backward recursion turns the whole line list into

    c_b  =  P(bin b is covered | every count in the list) ,

which is the coverage function asked for.  Because the recursion looks at the
whole list at once, a bin gets the benefit of the doubt from its neighbours: a
single empty bin inside a well-populated stretch is a fluctuation and keeps
c near 1, while an empty bin in the middle of a long desert is condemned by
the company it keeps.


4.  Where the baseline comes from
---------------------------------
rho cannot be estimated from the observed list alone.  The density of recorded
lines in a spectrum like this one swings by a factor of twenty between
neighbouring stretches, because transition arrays cluster in wavelength, and
that is real structure of the atom, not of the plates.  A smooth curve fitted
through the observed list either follows the structure - and then bends into a
gap and hides it - or does not, and then hands the structure to the blocked
state and calls half the spectrum blind.  Both failures happen in practice,
and which one you get depends on how flexible the curve is allowed to be.

The way out is that the CALCULATED transitions know the structure and know
nothing about the plates.  Their density per unit u, read from TRANS.DAT and
blurred by --pred-blur to allow for the error of a calculated wavelength,
supplies the SHAPE of the baseline:

    R_b  =  (calculated transitions in bin b)  *  exp(spline(u_b))

so that all the fitted curve has left to describe is the slow EFFICIENCY -
what fraction of the calculated lines a covered plate actually recorded, which
also absorbs whatever part of the list belongs to other species.  That is a
smooth, slowly varying thing, and it can be kept stiff without blaming the
atom's structure on the plates.  (--no-pred falls back to a spline-only
baseline, for comparison.)

The efficiency is fitted by Poisson regression of the counts on a penalized
cubic spline in u - --knots evenly spaced basis functions and a penalty
--smooth on the second differences of their coefficients - with

    expected exposure of bin b   e_b = c_b + (1 - c_b) * r0

as an offset, so that bins which currently look blocked do not drag it down.

How flexible that curve may be is not left to taste, because the answer
depends on it.  --knots is set generously; it is the PENALTY that controls the
flexibility, and it is chosen by the Bayesian information criterion over
--smooth-grid: twice the negative log-likelihood plus ln(number of bins) for
every EFFECTIVE parameter the spline spends.  The effective count is the trace
of the matrix that maps the data to the fit, not the number of coefficients,
because the penalty ties neighbouring coefficients together and a penalized
fit with two hundred coefficients may be spending thirty.

Baseline, leak-through, transition probabilities and coverage are all unknown
together, so they are found by iterating: guess the coverage, fit the
baseline, re-derive the coverage, and so on (the EM algorithm).  Twenty
rounds is more than enough; the log-likelihood printed each round shows it
settling.


5.  What comes out
------------------
  --out-map      one row per bin: its wavelength range, the lines recorded in
                 it, the lines the baseline expects, and c.
  --out-gaps     the runs of bins whose c falls below --report-threshold,
                 merged into segments, with the number of lines each segment
                 should have held and did not.  This is the table to read.
  --out-log      the running commentary, including the fitted r0, a and b.
  --plot         a picture of the recorded density, the baseline and c.

The gap table is a MEASUREMENT, not a list of plate boundaries: it says where
the record is thin, not why.  A segment in it may be a plate gap, a defect, an
impurity band - or, at the sparse red end of the list, a genuinely empty piece
of spectrum.  It is meant to be used as the continuous function c, which
weakens the evidence of a missing line in proportion to how thin the record
is; reading it as a set of hard boundaries would throw that proportion away.
"""

import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from calibrate_intensities import read_table, column, as_float   # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------------------
# input
# ---------------------------------------------------------------------------


def load_lines(opt, log):
    """The vacuum wavelengths of the observed lines, in angstroms, sorted."""
    header, rows = read_table(opt.lines, opt.lines_sheet)
    wn = [as_float(v) for v in column(header, rows, opt.col_wn, opt.lines)]
    wn = np.array(sorted(w for w in wn if w is not None and w > 0), float)
    if wn.size < 100:
        raise SystemExit('only %d usable wavenumbers in %s' % (wn.size,
                                                               opt.lines))
    lam = 1e8 / wn[::-1]
    log('lines read: %d, %.4f - %.4f A' % (lam.size, lam[0], lam[-1]))
    return lam


def load_predicted(opt, log):
    """Vacuum wavelengths of the CALCULATED transitions, with their strength.

    TRANS.DAT holds one row per predicted transition, each starting with '+'.
    Columns 6-10 are the predicted intensity on the file's logarithmic scale
    (v = 10*ln(intensity), so v is negative for a line weaker than 1) and
    columns 25-38 the predicted wavenumber in cm^-1.
    """
    lam, v = [], []
    with open(opt.pred, encoding='latin-1') as fh:
        for row in fh:
            if row[:1] != '+':
                continue
            try:
                wn, vv = float(row[24:38]), float(row[5:10])
            except ValueError:
                continue
            if wn > 0 and vv > opt.pred_vcut:
                lam.append(1e8 / wn)
                v.append(vv)
    lam = np.array(lam)
    log('predicted transitions read: %d with v > %g, %.1f - %.1f A'
        % (lam.size, opt.pred_vcut, lam.min(), lam.max()))
    return np.sort(lam)


def predicted_per_bin(edges, lam_p, opt, log):
    """The calculated transitions per bin, blurred by their own error.

    A calculated wavelength is only as good as the calculated energies behind
    it - for a level that has not been found experimentally the error can
    reach a few hundred cm^-1, which is a few angstroms in the ultraviolet.
    The count is therefore smoothed along ln(lambda) with a Gaussian of width
    --pred-blur before it is used, so that the baseline carries the SHAPE of
    the predicted density and not its accidental graininess.
    """
    n = np.bincount(np.clip(np.searchsorted(edges, np.log(lam_p), 'right') - 1,
                            0, edges.size - 2), minlength=edges.size - 1
                    ).astype(float)
    sig = opt.pred_blur / opt.bin
    if sig > 0.3:
        half = int(np.ceil(4 * sig))
        k = np.exp(-0.5 * (np.arange(-half, half + 1) / sig) ** 2)
        n = np.convolve(n, k / k.sum(), mode='same')
    log('predicted density: %.2f per bin on average, blurred by %.4f in '
        'ln(lambda)' % (n.mean(), opt.pred_blur))
    return np.maximum(n, opt.pred_floor)


def make_grid(lam, opt, log):
    """Bins of equal width in u = ln lambda, and the count in each."""
    u = np.log(lam)
    lo, hi = u[0], u[-1]
    nbin = int(np.ceil((hi - lo) / opt.bin))
    edges = lo + opt.bin * np.arange(nbin + 1)
    edges[-1] = max(edges[-1], hi + 1e-12)
    n = np.bincount(np.clip(np.searchsorted(edges, u, 'right') - 1,
                            0, nbin - 1), minlength=nbin).astype(float)
    mid = 0.5 * (edges[:-1] + edges[1:])
    log('grid: %d bins of %.5f in ln(lambda) (%.3f A at %.0f A, %.2f A at '
        '%.0f A); %.2f lines per bin on average'
        % (nbin, opt.bin, opt.bin * lam[0], lam[0], opt.bin * lam[-1],
           lam[-1], n.mean()))
    return edges, mid, n


# ---------------------------------------------------------------------------
# the smooth baseline: a penalized-spline Poisson fit of ln(rho)
# ---------------------------------------------------------------------------


def spline_basis(x, nknots, degree=3):
    """Cubic B-spline design matrix on an even knot grid spanning x."""
    from scipy.interpolate import BSpline
    lo, hi = float(x[0]), float(x[-1])
    step = (hi - lo) / (nknots - 1)
    knots = np.concatenate([lo - step * np.arange(degree, 0, -1),
                            np.linspace(lo, hi, nknots),
                            hi + step * np.arange(1, degree + 1)])
    return np.asarray(BSpline.design_matrix(np.clip(x, lo, hi), knots,
                                            degree, extrapolate=False)
                      .todense())


def diff_penalty(ncoef, order=2):
    """Second-difference penalty matrix: it charges for curvature."""
    d = np.eye(ncoef)
    for _ in range(order):
        d = np.diff(d, axis=0)
    return d.T @ d


def effective_df(basis, penalty, m, smooth):
    """How many free parameters the penalized spline really spends.

    A penalized fit with many coefficients is not as flexible as their number
    suggests, because the penalty ties neighbouring ones together.  The honest
    count is the trace of the matrix that maps the data to the fit,
    trace[(B'WB + s*P)^-1 B'WB], and it is what the information criterion that
    chooses the penalty has to charge for.
    """
    bwb = basis.T @ (basis * m[:, None])
    return float(np.trace(np.linalg.solve(bwb + smooth * penalty, bwb)))


def fit_baseline(basis, penalty, n, log_offset, beta, smooth, iters=40):
    """Poisson regression of the counts on the spline basis, penalized.

    Maximizes  sum_b [ n_b*(o_b + B*beta)_b - exp((o_b + B*beta)_b) ]
               - 0.5*smooth*beta'*penalty*beta
    by Newton's method (iteratively reweighted least squares).
    """
    for _ in range(iters):
        m = np.exp(log_offset + basis @ beta)
        grad = basis.T @ (n - m) - smooth * (penalty @ beta)
        hess = basis.T @ (basis * m[:, None]) + smooth * penalty
        step = np.linalg.solve(hess, grad)
        beta = beta + step
        if np.max(np.abs(step)) < 1e-9:
            break
    return beta


# ---------------------------------------------------------------------------
# the hidden Markov model over the bins
# ---------------------------------------------------------------------------


def poisson_logpmf(n, m):
    from scipy.special import gammaln
    return n * np.log(m) - m - gammaln(n + 1.0)


def forward_backward(logem, a, b):
    """Posterior state probabilities and expected transition counts.

    logem   (nbin, 2) log P(count in bin | state), state 0 = blocked,
            state 1 = covered
    a       P(covered -> blocked),   b  P(blocked -> covered)

    Returns (gamma, n10, n01, loglik) where gamma[:, 1] is the coverage,
    n10 the expected number of covered -> blocked steps and n01 the expected
    number of blocked -> covered ones.
    """
    nbin = logem.shape[0]
    trans = np.array([[1.0 - b, b], [a, 1.0 - a]])
    pi = np.array([a / (a + b), b / (a + b)])

    scale = logem.max(axis=1)
    em = np.exp(logem - scale[:, None])

    alpha = np.empty((nbin, 2))
    norm = np.empty(nbin)
    v = pi * em[0]
    norm[0] = v.sum()
    alpha[0] = v / norm[0]
    for i in range(1, nbin):
        v = (alpha[i - 1] @ trans) * em[i]
        norm[i] = v.sum()
        alpha[i] = v / norm[i]

    beta = np.empty((nbin, 2))
    beta[-1] = 1.0
    for i in range(nbin - 2, -1, -1):
        beta[i] = (trans @ (em[i + 1] * beta[i + 1])) / norm[i + 1]

    gamma = alpha * beta
    gamma /= gamma.sum(axis=1, keepdims=True)

    xi = np.zeros((2, 2))
    for i in range(nbin - 1):
        xi += (alpha[i][:, None] * trans
               * (em[i + 1] * beta[i + 1])[None, :]) / norm[i + 1]

    loglik = float(np.sum(np.log(norm)) + scale.sum())
    return gamma, float(xi[1, 0]), float(xi[0, 1]), loglik


# ---------------------------------------------------------------------------
# the whole estimate
# ---------------------------------------------------------------------------


def estimate(mid, n, pred, opt, log):
    """Iterate baseline, leak-through, transitions and coverage to a fixed
    point.  Returns (coverage, expected covered count per bin, r0, a, b).

    `pred` is the calculated transitions per bin, or None.  When it is given,
    the baseline is  pred_b * exp(spline)  and the spline has only the slow
    EFFICIENCY left to describe - what fraction of the calculated lines a
    covered plate actually recorded - so it can be kept stiff.  When it is
    None the spline has to carry the line density itself.
    """
    log_pred = np.log(pred) if pred is not None else np.zeros(mid.size)
    basis = spline_basis(mid, opt.knots)
    penalty = diff_penalty(basis.shape[1])
    log('baseline: %d spline coefficients, knot spacing %.4f in ln(lambda) '
        '(%.1f %% in wavelength), curvature penalty %.3g'
        % (basis.shape[1], (mid[-1] - mid[0]) / (opt.knots - 1),
           100.0 * (mid[-1] - mid[0]) / (opt.knots - 1), opt.smooth))

    beta = (np.zeros(basis.shape[1])
            + np.log(max(n.sum(), 1.0) / max(np.exp(log_pred).sum(), 1e-9)))
    cov = np.ones(mid.size)
    rate = np.full(mid.size, max(n.mean(), 1e-6))
    r0, a, b = opt.leak, opt.p_block, opt.p_unblock
    prev = None

    for it in range(1, opt.iters + 1):
        expo = cov + (1.0 - cov) * r0
        beta = fit_baseline(basis, penalty, n, log_pred + np.log(expo), beta,
                            opt.smooth)
        rate = np.exp(log_pred + basis @ beta)         # covered count per bin

        logem = np.column_stack([
            poisson_logpmf(n, np.maximum(r0 * rate, 1e-12)),
            poisson_logpmf(n, np.maximum(rate, 1e-12))])
        gamma, n10, n01, loglik = forward_backward(logem, a, b)
        cov = gamma[:, 1]

        a = float(np.clip(n10 / max(cov[:-1].sum(), 1e-9),
                          1e-7, opt.max_block))
        b = float(np.clip(n01 / max((1.0 - cov[:-1]).sum(), 1e-9),
                          1.0 / opt.max_run, 0.9))
        if not opt.fix_leak:
            r0 = float(np.clip(np.sum((1.0 - cov) * n)
                               / max(np.sum((1.0 - cov) * rate), 1e-9),
                               opt.leak_min, opt.leak_max))

        log('  round %2d  lnL %14.2f  leak-through %.4f  blocked fraction '
            '%.4f  mean run %.1f bins' % (it, loglik, r0, a / (a + b), 1.0 / b))
        if prev is not None and abs(loglik - prev) < opt.tol:
            log('  converged')
            break
        prev = loglik

    edf = effective_df(basis, penalty, np.maximum(rate * expo, 1e-12),
                       opt.smooth)
    bic = -2.0 * loglik + np.log(float(mid.size)) * (edf + 3.0)
    log('  effective parameters in the baseline %.1f; BIC %.1f' % (edf, bic))
    return cov, rate, r0, a, b, loglik, edf, bic


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------


def segments(edges, n, rate, cov, opt):
    """Merge the runs of bins with c below the threshold into segments."""
    below = cov < opt.report_threshold
    out, i = [], 0
    while i < below.size:
        if not below[i]:
            i += 1
            continue
        j = i
        while j + 1 < below.size and below[j + 1]:
            j += 1
        sl = slice(i, j + 1)
        lam_lo = float(np.exp(edges[i]))
        lam_hi = float(np.exp(edges[j + 1]))
        out.append(dict(lam_lo=lam_lo, lam_hi=lam_hi, width=lam_hi - lam_lo,
                        nbin=j - i + 1,
                        expected=float(rate[sl].sum()),
                        observed=float(n[sl].sum()),
                        lost=float(np.sum((rate[sl] - n[sl])
                                          * (1.0 - cov[sl]))),
                        c_min=float(cov[sl].min()),
                        c_mean=float(cov[sl].mean())))
        i = j + 1
    return out


def write_map(path, edges, n, rate, cov):
    with open(path, 'w', newline='\n') as fh:
        w = csv.writer(fh, lineterminator='\n')
        w.writerow(['lam_lo', 'lam_hi', 'n_obs', 'n_expected', 'coverage'])
        for i in range(n.size):
            w.writerow(['%.4f' % np.exp(edges[i]),
                        '%.4f' % np.exp(edges[i + 1]),
                        '%d' % int(n[i]), '%.4f' % rate[i], '%.5f' % cov[i]])


def write_gaps(path, segs, r0, a, b, opt):
    with open(path, 'w', newline='\n') as fh:
        fh.write('# wavelength stretches where the observed line list is '
                 'thinner than its surroundings.\n')
        fh.write('# produced by tools/coverage_map.py; coverage below %.2f; '
                 'leak-through %.4f, blocked fraction %.4f, mean run %.1f '
                 'bins\n' % (opt.report_threshold, r0, a / (a + b), 1.0 / b))
        fh.write('# lam_lo   lam_hi    width  expected  observed  lost  '
                 'c_min  c_mean\n')
        for s in segs:
            fh.write('%9.3f %9.3f %8.3f %9.1f %9d %7.1f %7.3f %7.3f\n'
                     % (s['lam_lo'], s['lam_hi'], s['width'], s['expected'],
                        int(round(s['observed'])), s['lost'], s['c_min'],
                        s['c_mean']))


def plot(path, mid, n, rate, cov):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    lam = np.exp(mid)
    nrow = 5
    cuts = np.exp(np.linspace(mid[0], mid[-1], nrow + 1))
    fig, axes = plt.subplots(nrow, 1, figsize=(13, 2.6 * nrow))
    for k, ax in enumerate(axes):
        lo, hi = cuts[k], cuts[k + 1]
        m = (lam >= lo) & (lam <= hi)
        ax.plot(lam[m], n[m], lw=0.4, color='0.65', label='lines recorded')
        ax.plot(lam[m], rate[m], lw=1.3, color='C0', label='baseline expected')
        ax2 = ax.twinx()
        ax2.fill_between(lam[m], 0, cov[m], color='C3', alpha=0.15, lw=0)
        ax2.plot(lam[m], cov[m], lw=0.8, color='C3')
        ax2.set_ylim(0, 1.05)
        ax2.set_ylabel('coverage c', color='C3')
        ax.set_xlim(lo, hi)
        ax.set_ylabel('lines per bin')
        if k == 0:
            ax.legend(loc='upper right', fontsize=8)
        if k == nrow - 1:
            ax.set_xlabel('vacuum wavelength, A')
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def parse_args(argv):
    p = argparse.ArgumentParser(
        description='map where a photographic line list is blind, from the '
                    'emptiness of the list itself',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    g = p.add_argument_group('input')
    g.add_argument('--lines', default=os.path.join(ROOT, 'Pr3_lines.xlsx'))
    g.add_argument('--lines-sheet', default='Sheet1')
    g.add_argument('--col-wn', default='own', help='observed wavenumber, cm^-1')
    g.add_argument('--pred', default=os.path.join(ROOT, 'IDEN2', 'trans.dat'),
                   help='calculated transitions, whose density gives the '
                        'baseline its shape')
    g.add_argument('--no-pred', action='store_true',
                   help='do not use the calculated transitions; let the '
                        'spline carry the line density on its own')
    g.add_argument('--pred-vcut', type=float, default=-20.0,
                   help='ignore calculated transitions weaker than this on '
                        "TRANS.DAT's logarithmic intensity scale")
    g.add_argument('--pred-blur', type=float, default=0.002,
                   help='smoothing of the calculated density along '
                        'ln(lambda), for the error of the calculated '
                        'wavelengths')
    g.add_argument('--pred-floor', type=float, default=0.02,
                   help='smallest calculated count a bin may be given')

    g = p.add_argument_group('the grid')
    g.add_argument('--bin', type=float, default=2e-4,
                   help='bin width in ln(lambda), i.e. as a fraction of the '
                        'wavelength')

    g = p.add_argument_group('the baseline')
    g.add_argument('--knots', type=int, default=200,
                   help='spline knots across the whole range.  It is the '
                        'curvature penalty, not this number, that sets how '
                        'flexible the baseline really is, so this only has '
                        'to be generous')
    g.add_argument('--smooth', type=float, default=1.0,
                   help='penalty on the curvature of ln(baseline)')
    g.add_argument('--select-smooth', action='store_true', default=True,
                   help='choose --smooth from --smooth-grid by BIC')
    g.add_argument('--no-select-smooth', dest='select_smooth',
                   action='store_false')
    g.add_argument('--smooth-grid', type=float, nargs='+',
                   default=[0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0,
                            3000.0, 10000.0, 30000.0])

    g = p.add_argument_group('the hidden Markov model')
    g.add_argument('--leak', type=float, default=0.10,
                   help='starting value of the leak-through r0')
    g.add_argument('--leak-min', type=float, default=0.02)
    g.add_argument('--leak-max', type=float, default=0.35)
    g.add_argument('--fix-leak', action='store_true',
                   help='keep r0 at --leak instead of estimating it')
    g.add_argument('--p-block', type=float, default=0.002,
                   help='starting P(covered -> blocked) per bin')
    g.add_argument('--p-unblock', type=float, default=0.05,
                   help='starting P(blocked -> covered) per bin')
    g.add_argument('--max-block', type=float, default=0.2,
                   help='ceiling on P(covered -> blocked)')
    g.add_argument('--max-run', type=float, default=1000.0,
                   help='longest mean run of blocked bins the fit may claim')
    g.add_argument('--iters', type=int, default=25)
    g.add_argument('--tol', type=float, default=1e-3)

    g = p.add_argument_group('output')
    g.add_argument('--report-threshold', type=float, default=0.5,
                   help='a bin is listed in the gap table below this coverage')
    g.add_argument('--out-map', default=os.path.join(ROOT, 'coverage_map.csv'))
    g.add_argument('--out-gaps',
                   default=os.path.join(ROOT, 'coverage_gaps.txt'))
    g.add_argument('--out-log', default=os.path.join(ROOT, 'coverage_map.log'))
    g.add_argument('--plot', default=os.path.join(ROOT, 'coverage_map.png'))
    g.add_argument('--no-plot', action='store_true')
    return p.parse_args(argv)


def main(argv=None):
    opt = parse_args(argv)
    fh = open(opt.out_log, 'w', newline='\n')

    def log(msg=''):
        print(msg)
        fh.write(msg + '\n')

    lam = load_lines(opt, log)
    edges, mid, n = make_grid(lam, opt, log)
    pred = None
    if not opt.no_pred:
        pred = predicted_per_bin(edges, load_predicted(opt, log), opt, log)
    if opt.select_smooth:
        log()
        log('choosing the curvature penalty by BIC (Bayesian information '
            'criterion: twice the negative log-likelihood, plus ln(number of '
            'bins) for every effective parameter the baseline spends; the '
            'smallest wins):')
        log('     penalty       lnL      params       BIC')
        best = None
        for sm in opt.smooth_grid:
            opt.smooth = sm
            res = estimate(mid, n, pred, opt, lambda m='': None)
            log('  %10.4g %10.1f %10.1f %10.1f' % (sm, res[5], res[6], res[7]))
            if best is None or res[7] < best[0]:
                best = (res[7], sm)
        opt.smooth = best[1]
        log('  chosen: penalty %g' % opt.smooth)
        log()
    cov, rate, r0, a, b = estimate(mid, n, pred, opt, log)[:5]

    write_map(opt.out_map, edges, n, rate, cov)
    segs = segments(edges, n, rate, cov, opt)
    write_gaps(opt.out_gaps, segs, r0, a, b, opt)

    log()
    log('blocked stretches (coverage below %.2f), worst first:'
        % opt.report_threshold)
    log('   lam_lo    lam_hi    width  expected  observed    lost   c_min')
    for s in sorted(segs, key=lambda s: -s['lost'])[:40]:
        log('%9.3f %9.3f %8.3f %9.1f %9d %7.1f %7.3f'
            % (s['lam_lo'], s['lam_hi'], s['width'], s['expected'],
               int(round(s['observed'])), s['lost'], s['c_min']))
    log()
    log('%d segments; %.1f lines lost in total, %.2f %% of the %.0f the '
        'baseline expects over the whole range'
        % (len(segs), sum(s['lost'] for s in segs),
           100.0 * sum(s['lost'] for s in segs) / rate.sum(), rate.sum()))
    log('written: %s, %s' % (opt.out_map, opt.out_gaps))

    if not opt.no_plot:
        plot(opt.plot, mid, n, rate, cov)
        log('written: %s' % opt.plot)
    fh.close()


if __name__ == '__main__':
    main()
