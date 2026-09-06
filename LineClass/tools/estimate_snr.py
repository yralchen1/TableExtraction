#!/usr/bin/env python
"""Estimate the noise level of an observed line list, and turn it into SNR.

This is a standalone tool.  It imports nothing from the pipeline except the
region-fitting machinery of its sibling `calibrate_intensities.py`, and takes
every file, sheet and column name as a command-line option.


1.  What the tool is for
------------------------
An old photographic line list gives, for every line it records, an intensity:
a number read off the plate.  What it does NOT give is the NOISE - the
brightness a line must reach before it appears on the plate at all.  Without
that number, one cannot say whether a transition that theory predicts but the
list does not contain is genuinely absent or merely too weak to have been
recorded, and one cannot say how solid a line that IS recorded is.

The ratio that answers both questions is the signal-to-noise ratio

    SNR(line)  =  I(line) / N(lambda) ,

with

    I(line)   = the intensity of the line, on the scale the plate recorded it,
    lambda    = its vacuum wavelength in angstroms (1 A = 1e-8 cm), obtained
                from the wavenumber wn in cm^-1 as lambda = 1e8/wn,
    N(lambda) = the noise level at that wavelength, on the SAME scale.

This tool estimates N, computes SNR for every observed line, and - if asked -
writes the SNR into the two working files of the identification program IDEN2
in place of the intensity, so that the display can be thinned out: a predicted
transition whose SNR is below 1 could not have been seen, and there is no
point in it cluttering the screen or occupying a column.

The scale used throughout is the ORIGINAL reported intensity, not the
intensity after the plate calibration has been applied.  The reason is that
noise is a property of the plate: it is the reported number that is censored
at the detection threshold, so it is on the reported scale that the threshold
is a smooth function of wavelength.  SNR is a ratio, so it does not matter for
the answer which scale it is computed on, as long as numerator and denominator
are on the same one; but the ESTIMATION of N is only clean on the original.


2.  How the noise level is estimated
------------------------------------
A line list is a censored sample: every line in it was above the detection
threshold, and everything below the threshold is missing.  The weakest lines
present are therefore the ones sitting just above the threshold, and the
threshold can be read off as the lower envelope of the recorded intensities.

    step 1  Sort the lines by wavelength and cut them at the COVERAGE GAPS -
            the wavelength intervals between two photographic exposures whose
            ranges do not overlap (--gap; the defaults are Sugar's two).  No
            line was recorded on both of the plates that meet at a gap, so
            nothing ties their scales - or their noise levels - together, and
            the estimate must not be carried across.
    step 2  Slide a window of --window consecutive lines along each block and
            take the WEAKEST line in the window as the noise level there,
            attributing it to the median wavelength of the window.  Twenty
            lines is enough that at least one of them is likely to sit near
            the threshold, and few enough that the window is short compared
            with the scale on which the plate sensitivity changes.
            The windows do not overlap by default (--step defaults to the
            window size).  Overlapping windows would repeat the same minimum
            many times over, and the region-finding of step 3 charges for
            every breakpoint against the NUMBER of points, so repeating them
            would make it believe it has far more evidence than it does.
    step 3  Fit ln N against lambda with a piecewise polynomial of low degree,
            using the machinery of calibrate_intensities.py: the wavelength
            range is cut into contiguous regions, the cuts being chosen by the
            Bayesian information criterion, and inside each region ln N is an
            ordinary polynomial.  The fit is in logarithms because the noise
            spans decades and because a polynomial in ln is what keeps N
            positive everywhere.
    step 4  Discard outliers and refit, up to --max-passes times.  The
            rejection is ASYMMETRIC, and deliberately so.  A window whose
            weakest line sits far ABOVE the fitted curve is a window where the
            lines happened to be sparse, so its minimum overestimates the
            threshold; it is dropped when the excess exceeds --high-cut.  A
            window whose minimum sits BELOW the curve is a window that
            happened to catch a line very close to the threshold - which is
            the quantity being estimated - so it is kept unless it is absurd,
            beyond --low-cut, which is set much looser.  Rejecting the two
            sides symmetrically would pull the estimate up towards the middle
            of the intensity distribution, which is not what a threshold is.

The result is written in the same plain-text format as the intensity
correction: one region per line,

    lambda_lo  lambda_hi  c0;c1;c2;...

meaning ln N(lambda) = c0 + c1*lambda + c2*lambda^2 + ... inside that region.


3.  Writing the SNR into the IDEN2 files
----------------------------------------
IDEN2 (Azarov, Kramida & Vokhmentsev, Comput. Phys. Commun. 225 (2018) 149)
reads two working files that carry an intensity in a fixed-width integer
field:

    dlv.dat    one row per OBSERVED line; the field is columns 1-5, and the
               wavenumber of the line is columns 6-18.
    TRANS.DAT  one row per PREDICTED transition, grouped under a header row
               for each upper level.  The row carries TWO intensities and TWO
               wavenumbers, and they must not be confused:

                   columns  6-10   the PREDICTED intensity
                   columns 11-22   the energy of the level at the other end
                   columns 25-38   the wavenumber of the PREDICTED transition
                   columns 39-43   the intensity of the OBSERVED line the
                                   prediction is identified with, if any
                   columns 45-57   the wavenumber of that observed line
                   columns 69-74   its row number in dlv.dat, 0 if none

               Columns 23-24 hold an asterisk when the level at the other end
               has an experimentally known energy, so the wavenumber field
               cannot be found by splitting the row on spaces.  The wavelength
               at which the plate calibration and the noise are evaluated is
               1e8 / (the wavenumber in columns 25-38) - NOT 1e8 / (the energy
               in columns 11-22), which is a different line of the spectrum
               altogether and, for the transitions to low-lying levels that
               make up the short-wavelength end of this list, is out of the
               calibrated range entirely.

Both hold the intensity on a logarithmic scale,

    v  =  k * ln(Icor) ,     k = --iden2-k, default 10,

where Icor is the intensity after the plate calibration in force when the
files were made has been applied, Icor = S * Iorig * exp(P(lambda)) with
S = --trans-scale and P the piecewise correction named by --trans-correction.
(That the relation is exactly this one, for this project's files, was checked
by regressing the dlv.dat integers on ln(Icor): slope 10.016, intercept -0.19,
scatter 0.10 in ln, i.e. rounding.)

  * In dlv.dat the observed line is looked up in the line list by wavenumber
    and its SNR written in place of v.  Nothing else on the row is touched.
  * In TRANS.DAT the predicted intensity is turned back into an intensity on
    the plate's own scale by undoing that relation,

        Icor_pred  = exp(v/k)
        Iorig_pred = Icor_pred / (S * exp(P(lambda)))
        SNR_pred   = Iorig_pred / N(lambda) ,

    and the result written in its place.  This is why --trans-correction must
    name the calibration THE FILE WAS MADE WITH, not the one now in force: the
    number in the file has that one baked into it, and it has to come out
    again.
  * The observed intensity in the same TRANS.DAT row is replaced by the SNR of
    that same line, taken from the line list exactly as dlv.dat is.  The two
    columns sit next to each other on the screen and are read against each
    other, so leaving one of them on the plate scale would show a difference
    of some sixty units at the short-wavelength end that is entirely an
    artefact of the two scales, and none of it a disagreement between the
    prediction and the measurement.

The scale of the number written is chosen by --scale:

    linear   SNR itself, rounded to an integer.
    log      round(k * ln(SNR)), which is negative below SNR = 1.
    auto     linear if every observed line falls between 1 and 999, so that
             the column reads directly as "how many times the noise"; log
             otherwise.  (default)

The rewritten files go to --out-dir, never over the originals.


4.  What the estimate is and is not
-----------------------------------
N is a DETECTION THRESHOLD read off the line list itself, not a measurement of
the noise on the plate.  It inherits every property of the list: where the
compiler of the list stopped measuring weak lines, N rises, whether or not the
plate was any noisier there; where a region of the spectrum is crowded, the
weakest recorded line is closer to the true threshold than where it is empty.
It is the right quantity for the question "could this predicted transition
have appeared in THIS list", which is the question being asked of it, and it
should not be quoted as anything else.

Usage
-----
    python tools/estimate_snr.py
    python tools/estimate_snr.py --window 30 --plot snr.png
    python tools/estimate_snr.py --iden2-dir IDEN2 --out-dir IDEN2_snr
"""
from __future__ import annotations

import argparse
import bisect
import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from calibrate_intensities import (          # noqa: E402
    DEFAULT_GAPS_A, Region, apply_correction, build_regions, read_correction,
    read_table, column, as_float, write_correction,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------------------
# the observed lines
# ---------------------------------------------------------------------------


def load_lines(opt, log):
    """The observed lines as (wavenumber, original intensity), sorted by wn."""
    header, rows = read_table(opt.lines, opt.lines_sheet)
    wn = [as_float(v) for v in column(header, rows, opt.col_wn, opt.lines)]
    io = [as_float(v) for v in column(header, rows, opt.col_orig, opt.lines)]
    keep = [(w, i) for w, i in zip(wn, io)
            if w is not None and w > 0 and i is not None and i > 0]
    keep.sort()
    if not keep:
        raise SystemExit('no line has both a wavenumber and a positive %s'
                         % opt.col_orig)
    w = np.array([k[0] for k in keep], float)
    i = np.array([k[1] for k in keep], float)
    log('lines read: %d with a wavenumber and a positive %s, %.2f - %.2f A'
        % (len(w), opt.col_orig, 1e8 / w[-1], 1e8 / w[0]))
    return w, i


def gap_split(lam, gaps):
    """Index ranges of the blocks the coverage gaps cut the sorted lam into."""
    blocks, start = [], 0
    for glo, ghi in gaps:
        i = int(np.searchsorted(lam, glo, side='left'))
        j = int(np.searchsorted(lam, ghi, side='right'))
        if i > start:
            blocks.append((start, i))
        start = max(start, j)
    if start < len(lam):
        blocks.append((start, len(lam)))
    return blocks


def noise_points(lam, inten, opt, log):
    """One (wavelength, weakest intensity) point per window, block by block."""
    step = opt.step if opt.step else opt.window
    xs, ys, ns = [], [], []
    for a, b in gap_split(lam, opt.gaps):
        n = b - a
        if n < opt.window:
            log('    block %.2f - %.2f A holds only %d lines, fewer than the '
                'window of %d; its single weakest line is used'
                % (lam[a], lam[b - 1], n, opt.window))
            xs.append(float(np.median(lam[a:b])))
            ys.append(float(inten[a:b].min()))
            ns.append(n)
            continue
        k = a
        while k + opt.window <= b:
            sl = slice(k, k + opt.window)
            xs.append(float(np.median(lam[sl])))
            ys.append(float(inten[sl].min()))
            ns.append(opt.window)
            k += step
        if k < b:                     # the tail, folded into the last window
            sl = slice(b - opt.window, b)
            xs.append(float(np.median(lam[sl])))
            ys.append(float(inten[sl].min()))
            ns.append(opt.window)
    x = np.array(xs, float)
    order = np.argsort(x)
    log('    %d windows of %d lines (step %d) over %d blocks'
        % (len(xs), opt.window, step, len(gap_split(lam, opt.gaps))))
    return x[order], np.array(ys, float)[order]


# ---------------------------------------------------------------------------
# the fit
# ---------------------------------------------------------------------------


class _Seg(object):
    """The few attributes build_regions() reads off its options object."""

    def __init__(self, opt, lo, hi):
        self.max_degree = opt.max_degree
        self.min_points = opt.min_points
        self.pts_per_coef = opt.pts_per_coef
        self.split_penalty = opt.split_penalty
        self.min_jump = 0.0
        self.split_degree = min(opt.max_degree, 3)
        self.max_split_tries = 150
        self.max_depth = 6
        self.max_swing = opt.max_swing
        self.max_extrap = opt.max_extrap
        self.gaps = opt.gaps
        self.cover = (lo, hi)


def fit_noise(x, y, opt, cover, log):
    """Fit ln(noise) against wavelength; returns (regions, keep-mask)."""
    ly = np.log(y)
    keep = np.ones(len(x), bool)
    seg = _Seg(opt, *cover)
    regions = build_regions(x[keep], ly[keep], seg, cover, log)
    for npass in range(1, opt.max_passes + 1):
        resid = ly - apply_correction(regions, x)
        high = keep & (resid > opt.high_cut)
        low = keep & (resid < -opt.low_cut)
        log('    pass %d: %d windows above the curve by more than %.2f, '
            '%d below it by more than %.2f'
            % (npass, int(high.sum()), opt.high_cut,
               int(low.sum()), opt.low_cut))
        drop = high | low
        if not drop.any():
            break
        keep &= ~drop
        regions = build_regions(x[keep], ly[keep], seg, cover, log)
    resid = (ly - apply_correction(regions, x))[keep]
    log('    %d of %d windows kept; rms of ln N about the curve %.3f'
        % (int(keep.sum()), len(x), float(np.sqrt((resid ** 2).mean()))))
    return regions, keep


# ---------------------------------------------------------------------------
# the IDEN2 files
# ---------------------------------------------------------------------------


class Matcher(object):
    """Nearest observed line by wavenumber, within a tolerance."""

    def __init__(self, wn, val):
        order = np.argsort(wn)
        self.wn = list(np.asarray(wn, float)[order])
        self.val = list(np.asarray(val, float)[order])

    def __call__(self, w, tol):
        i = bisect.bisect_left(self.wn, w)
        best, bd = None, float('inf')
        for j in (i - 1, i, i + 1):
            if 0 <= j < len(self.wn):
                d = abs(self.wn[j] - w)
                if d < bd:
                    best, bd = j, d
        return (self.val[best], bd) if best is not None and bd <= tol else (None, bd)


def encode(snr, scale, k):
    """The integer IDEN2 is to carry for this signal-to-noise ratio."""
    if snr is None or not np.isfinite(snr) or snr <= 0:
        return -9999 if scale == 'log' else 0
    if scale == 'linear':
        return int(round(min(snr, 99999.0)))
    return int(round(max(min(k * math.log(snr), 99999.0), -9999.0)))


def rewrite_dlv(src, dst, matcher, scale, k, tol, log):
    n, miss = 0, 0
    out = []
    for raw in open(src, encoding='latin-1', newline=''):
        line = raw.rstrip('\r\n')
        if not line.strip():
            out.append(raw)
            continue
        try:
            wn = float(line[5:18])
        except ValueError:
            out.append(raw)
            continue
        snr, _ = matcher(wn, tol)
        if snr is None:
            miss += 1
            out.append(raw)
            continue
        eol = raw[len(line):]
        out.append('%5d%s%s' % (encode(snr, scale, k), line[5:], eol))
        n += 1
    with open(dst, 'w', encoding='latin-1', newline='') as fh:
        fh.writelines(out)
    log('  %s: %d rows rewritten, %d left as they were (no observed line '
        'within %.3f cm^-1)' % (os.path.basename(dst), n, miss, tol))


TRANS_V     = (5, 10)     # the predicted intensity, k*ln(Icor)
TRANS_EPART = (10, 22)    # the energy of the level at the other end
TRANS_WN    = (24, 38)    # the wavenumber of the PREDICTED transition
TRANS_VOBS  = (38, 43)    # the intensity of the observed line assigned to it
TRANS_WNOBS = (44, 57)    # the wavenumber of that observed line
TRANS_LINE  = (68, 74)    # its row number in dlv.dat; 0 when none is assigned


def trans_field(line, span):
    """The number in a fixed-width TRANS.DAT field, or None if it is blank."""
    try:
        return float(line[span[0]:span[1]])
    except ValueError:
        return None


def trans_zero(src, regions, corr, opt, matcher, log):
    """The constant that puts the predicted intensities on the observed scale.

    The predicted intensities in TRANS.DAT carry an arbitrary overall
    normalisation - they come from calculated transition probabilities times a
    population factor, with no absolute calibration - so dividing them by the
    noise gives a number proportional to a signal-to-noise ratio but not equal
    to one, and the useful threshold "SNR = 1" would sit anywhere.  The offset
    is measured on the rows TRANS.DAT itself marks as identified: each of them
    names the observed line assigned to the prediction, so the two sides can
    be compared without any matching of our own, and the median difference in
    logarithms is the offset.  Rows with no identification are compared by
    wavenumber instead, which is what happens when the file carries none.

    That median is biased upward - a prediction is more likely to have found
    an observed line when it is strong - so it is a conservative zero point:
    it makes the tool call a transition unobservable slightly more readily
    than the truth warrants.  Give --trans-zero a number to override it, or 0
    to leave the predicted scale exactly as it is.
    """
    diffs = []
    assigned = 0
    for raw in open(src, encoding='latin-1', newline=''):
        line = raw.rstrip('\r\n')
        if not line or line[0] != '+':
            continue
        v = trans_field(line, TRANS_V)
        wn = trans_field(line, TRANS_WN)
        if v is None or wn is None or wn <= 0:
            continue
        wn_obs = trans_field(line, TRANS_WNOBS)
        n_obs = trans_field(line, TRANS_LINE)
        if n_obs and wn_obs:
            obs, _ = matcher(wn_obs, opt.match_tol)
            assigned += 1
        else:
            obs, _ = matcher(wn, opt.match_tol)
        if obs is None or obs <= 0:
            continue
        lam = np.array([1e8 / wn])
        icor = math.exp(v / opt.iden2_k)
        iorig = icor / (opt.trans_scale * math.exp(apply_correction(corr, lam)[0]))
        noise = math.exp(apply_correction(regions, lam)[0])
        diffs.append(math.log(iorig / noise) - math.log(obs))
    if not diffs:
        log('    no predicted transition matched an observed line; the '
            'predicted scale is left as it is')
        return 0.0
    z = -float(np.median(diffs))
    log('    %d predicted transitions have an observed line (%d of them '
        'identified as such in the file itself); their predicted SNR is a '
        'factor %.3g away from the observed one (median), and the predicted '
        'scale is shifted by that much so that SNR = 1 means "at the '
        'threshold"' % (len(diffs), assigned, math.exp(-z)))
    return z


def rewrite_trans(src, dst, regions, corr, opt, matcher, log):
    """Write the SNR into both intensity columns of TRANS.DAT.

    The predicted intensity is converted; the intensity of the observed line
    the row is identified with is looked up in the line list, exactly as
    dlv.dat is, so that the two columns the user reads side by side are on one
    scale.  Leaving the observed column alone would show an SNR next to a
    plate-calibrated intensity and invite the comparison to be read as a
    disagreement between prediction and measurement.
    """
    n, n_obs, miss = 0, 0, 0
    out = []
    for raw in open(src, encoding='latin-1', newline=''):
        line = raw.rstrip('\r\n')
        if not line or line[0] != '+':
            out.append(raw)
            continue
        v = trans_field(line, TRANS_V)
        wn = trans_field(line, TRANS_WN)
        if v is None or wn is None or wn <= 0:
            out.append(raw)
            continue
        eol = raw[len(line):]
        lam = np.array([1e8 / wn])
        icor = math.exp(v / opt.iden2_k)
        iorig = icor / (opt.trans_scale * math.exp(apply_correction(corr, lam)[0]))
        noise = math.exp(apply_correction(regions, lam)[0])
        snr = iorig / noise * math.exp(opt.trans_zero_value)
        line = ('%s%5d%s'
                % (line[:TRANS_V[0]], encode(snr, opt.scale, opt.iden2_k),
                   line[TRANS_V[1]:]))
        wn_obs = trans_field(line, TRANS_WNOBS)
        if trans_field(line, TRANS_LINE) and wn_obs:
            obs, _ = matcher(wn_obs, opt.match_tol)
            if obs is None:
                miss += 1
            else:
                line = ('%s%5d%s'
                        % (line[:TRANS_VOBS[0]],
                           encode(obs, opt.scale, opt.iden2_k),
                           line[TRANS_VOBS[1]:]))
                n_obs += 1
        out.append(line + eol)
        n += 1
    with open(dst, 'w', encoding='latin-1', newline='') as fh:
        fh.writelines(out)
    log('  %s: %d predicted transitions rewritten, and the observed intensity '
        'of %d identified rows (%d had no line within %.3f cm^-1)'
        % (os.path.basename(dst), n, n_obs, miss, opt.match_tol))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def parse_args(argv):
    p = argparse.ArgumentParser(
        description='estimate the noise level of a line list and express the '
                    'intensities as signal-to-noise ratios',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    g = p.add_argument_group('input')
    g.add_argument('--lines', default=os.path.join(ROOT, 'Pr3_lines.xlsx'))
    g.add_argument('--lines-sheet', default='Sheet1')
    g.add_argument('--col-wn', default='own', help='observed wavenumber, cm^-1')
    g.add_argument('--col-orig', default='Iorig',
                   help='the intensity AS REPORTED, before the plate '
                        'calibration; the noise is estimated on this scale')
    g.add_argument('--cover', type=float, nargs=2, default=None,
                   metavar=('LO', 'HI'),
                   help='wavelength interval, in angstroms, the finished noise '
                        'function must cover; the range of the lines by default')
    g.add_argument('--gap', dest='gaps', action='append', nargs=2, type=float,
                   metavar=('LO', 'HI'), default=None,
                   help='a wavelength interval no exposure covered: a forced '
                        'boundary, across which the noise may jump.  Repeat '
                        'for several; the default is this spectrum\'s two, %s'
                        % ', '.join('%.2f-%.2f' % g for g in DEFAULT_GAPS_A))
    g.add_argument('--no-gaps', action='store_true')

    g = p.add_argument_group('the estimate')
    g.add_argument('--window', type=int, default=20,
                   help='lines per window; the weakest of them is the noise')
    g.add_argument('--step', type=int, default=0,
                   help='lines between two windows; 0 = the window size, i.e. '
                        'windows that do not overlap')
    g.add_argument('--max-degree', type=int, default=3)
    g.add_argument('--min-points', type=int, default=25,
                   help='fewest windows a region may contain')
    g.add_argument('--pts-per-coef', type=int, default=8)
    g.add_argument('--split-penalty', type=float, default=2.0)
    g.add_argument('--max-swing', type=float, default=1.0)
    g.add_argument('--max-extrap', type=float, default=0.5)
    g.add_argument('--high-cut', type=float, default=1.0,
                   help='drop a window whose weakest line sits this far ABOVE '
                        'the fitted curve, in natural logarithms: its lines '
                        'were too sparse to reach the threshold')
    g.add_argument('--low-cut', type=float, default=2.5,
                   help='drop a window this far BELOW the curve.  Set much '
                        'looser than --high-cut on purpose: a low point is '
                        'evidence about the threshold, a high one is not')
    g.add_argument('--max-passes', type=int, default=5)

    g = p.add_argument_group('IDEN2')
    g.add_argument('--iden2-dir', default=None, metavar='DIR',
                   help='directory holding dlv.dat and TRANS.DAT; without it '
                        'the tool only estimates the noise')
    g.add_argument('--out-dir', default=None, metavar='DIR',
                   help='where the rewritten copies go (never over the '
                        'originals); DIR_snr by default')
    g.add_argument('--scale', choices=('auto', 'linear', 'log'), default='auto',
                   help='how the SNR is written: linear, log = k*ln(SNR), or '
                        'auto = linear if every observed line lies between 1 '
                        'and 999')
    g.add_argument('--iden2-k', type=float, default=10.0,
                   help='the k of the logarithmic scale v = k*ln(I) these '
                        'files carry')
    g.add_argument('--trans-correction',
                   default=os.path.join(ROOT, 'intensity_correction_manual.txt'),
                   help='the piecewise plate calibration the IDEN2 files were '
                        'MADE with - it is baked into their intensities and '
                        'has to be undone, so this is not necessarily the '
                        'calibration now in force')
    g.add_argument('--trans-scale', type=float, default=1000.0,
                   help='the overall factor of that calibration, '
                        'Icor = S*Iorig*exp(P)')
    g.add_argument('--trans-zero', default='auto',
                   help='zero point of the predicted intensity scale, in '
                        'natural logarithms: "auto" measures it against the '
                        'observed lines (see trans_zero), a number sets it, '
                        '0 leaves the predicted scale untouched')
    g.add_argument('--match-tol', type=float, default=0.02,
                   help='how far a dlv.dat wavenumber may miss its line, cm^-1')

    g = p.add_argument_group('output')
    g.add_argument('--out-noise', default=os.path.join(ROOT, 'noise_level.txt'),
                   help='the fitted ln N(lambda), region by region')
    g.add_argument('--out-lines', default=os.path.join(ROOT, 'line_snr.csv'),
                   help='per-line wavenumber, wavelength, intensity, noise, SNR')
    g.add_argument('--out-windows', default=os.path.join(ROOT, 'noise_windows.csv'),
                   help='the window minima the fit rests on')
    g.add_argument('--plot', default=None, metavar='PNG')
    g.add_argument('--quiet', action='store_true')

    opt = p.parse_args(argv)
    opt.gaps = () if opt.no_gaps else tuple(
        tuple(sorted(g)) for g in (opt.gaps or [list(g) for g in DEFAULT_GAPS_A]))
    opt.gaps = tuple(sorted(opt.gaps))
    return opt


def main(argv=None):
    opt = parse_args(argv)

    def log(msg):
        if not opt.quiet:
            print(msg)

    wn, inten = load_lines(opt, log)
    lam = 1e8 / wn
    order = np.argsort(lam)
    lam, wn_s, inten_s = lam[order], wn[order], inten[order]
    cover = tuple(opt.cover) if opt.cover else (float(lam[0]), float(lam[-1]))

    log('')
    log('noise windows')
    x, y = noise_points(lam, inten_s, opt, log)

    floor = float(inten_s.min())
    at_floor = float((y <= floor).mean())
    log('    %.0f%% of the windows have their weakest line at the smallest '
        'intensity the list uses at all (%g)' % (100 * at_floor, floor))
    if at_floor > 0.5:
        log('    NOTE: the reported intensities are themselves a per-plate '
            'scale - the compiler of the list graded each exposure from its '
            'own weakest visible line upward - so the threshold on THIS scale '
            'is nearly flat, and the SNR of an observed line is close to its '
            'reported intensity.  That is not a failure of the estimate: it '
            'means the whole wavelength dependence sits in the plate '
            'calibration, and the useful half of the exercise is the other '
            'one - dividing a PREDICTED intensity by exp(P(lambda)) to bring '
            'it down onto this scale, which is what the TRANS.DAT rewrite '
            'below does.  Where the estimate does rise above the floor, the '
            'lines are sparse there and the weakest one recorded is not '
            'necessarily near the threshold.')

    log('')
    log('fitting ln N(lambda)')
    regions, keep = fit_noise(x, y, opt, cover, log)

    noise = np.exp(apply_correction(regions, lam))
    snr = inten_s / noise
    lo, hi = float(snr.min()), float(snr.max())
    log('')
    log('signal-to-noise ratio of the observed lines: %.3g to %.3g, '
        'median %.3g' % (lo, hi, float(np.median(snr))))
    for q in (1, 5, 25, 50, 75, 95, 99):
        log('    %2d%% of the lines are below SNR %8.2f'
            % (q, float(np.percentile(snr, q))))
    below = int((snr < 1).sum())
    if below:
        log('    %d lines fall below SNR 1; the estimate is a threshold read '
            'off the list itself, so a few of the weakest must' % below)

    scale = opt.scale
    if scale == 'auto':
        scale = 'linear' if (lo >= 1.0 and hi <= 999.0) else 'log'
        log('    --scale auto -> %s (%s)' % (
            scale, 'the whole range fits between 1 and 999'
            if scale == 'linear' else
            'the range does not fit between 1 and 999'))
    opt.scale = scale

    log('')
    log('the fitted noise level')
    log('  region       from          to    deg      n      N(from)      N(to)')
    for i, r in enumerate(regions):
        log('  %6d %10.2f %11.2f %6d %6d %12.4g %10.4g'
            % (i + 1, r.lo, r.hi, r.degree, r.n,
               math.exp(float(r(r.lo))), math.exp(float(r(r.hi)))))
    write_correction(opt.out_noise, regions)
    log('  written to %s' % opt.out_noise)

    with open(opt.out_lines, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['wn', 'lambda_A', opt.col_orig, 'noise', 'snr', 'iden2'])
        for a, b, c, d, e in zip(wn_s, lam, inten_s, noise, snr):
            w.writerow(['%.4f' % a, '%.4f' % b, '%g' % c, '%.4g' % d,
                        '%.4g' % e, encode(e, scale, opt.iden2_k)])
    log('  per-line SNR written to %s' % opt.out_lines)

    with open(opt.out_windows, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['lambda_A', 'weakest_line', 'fitted_noise', 'used'])
        fitted = np.exp(apply_correction(regions, x))
        for a, b, c, u in zip(x, y, fitted, keep):
            w.writerow(['%.4f' % a, '%g' % b, '%.4g' % c, int(u)])
    log('  window minima written to %s' % opt.out_windows)

    if opt.iden2_dir:
        out_dir = opt.out_dir or (opt.iden2_dir.rstrip('/\\') + '_snr')
        os.makedirs(out_dir, exist_ok=True)
        log('')
        log('IDEN2 files rewritten with SNR on the %s scale, into %s'
            % (scale, out_dir))
        corr = read_correction(opt.trans_correction)
        matcher = Matcher(wn_s, snr)
        tsrc = None
        for name in os.listdir(opt.iden2_dir):
            if name.lower() == 'trans.dat':
                tsrc = os.path.join(opt.iden2_dir, name)
        if opt.trans_zero == 'auto':
            opt.trans_zero_value = (
                trans_zero(tsrc, regions, corr, opt, matcher, log)
                if tsrc else 0.0)
        else:
            opt.trans_zero_value = float(opt.trans_zero)
            log('    predicted scale shifted by %+.3f in ln (--trans-zero)'
                % opt.trans_zero_value)
        for name in sorted(os.listdir(opt.iden2_dir)):
            src = os.path.join(opt.iden2_dir, name)
            dst = os.path.join(out_dir, name)
            low = name.lower()
            if low == 'dlv.dat':
                rewrite_dlv(src, dst, matcher, scale, opt.iden2_k,
                            opt.match_tol, log)
            elif low == 'trans.dat':
                rewrite_trans(src, dst, regions, corr, opt, matcher, log)
            elif os.path.isfile(src):
                with open(src, 'rb') as a, open(dst, 'wb') as b:
                    b.write(a.read())
                log('  %s: copied unchanged' % name)

    if opt.plot:
        make_plot(opt.plot, lam, inten_s, x, y, keep, regions, log)
    return 0


def make_plot(path, lam, inten, x, y, keep, regions, log):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(lam, inten, '.', ms=1.5, color='0.75', label='observed lines')
    ax.plot(x[keep], y[keep], 'o', ms=3, color='tab:blue',
            label='window minima used')
    if (~keep).any():
        ax.plot(x[~keep], y[~keep], 'x', ms=5, color='tab:red',
                label='window minima dropped')
    for r in regions:
        g = np.linspace(r.lo, r.hi, 400)
        ax.plot(g, np.exp(r(g)), '-', color='tab:orange', lw=2)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('vacuum wavelength, A')
    ax.set_ylabel('reported intensity')
    ax.set_title('noise level of the line list (orange) under the recorded lines')
    ax.legend(loc='best', fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    log('  plot written to %s' % path)


if __name__ == '__main__':
    sys.exit(main())
