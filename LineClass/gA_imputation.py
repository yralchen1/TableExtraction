#!/usr/bin/env python
"""Value to give a transition that is absent from the calculated-intensity file.

Terms used throughout
---------------------
gA
    The statistical weight of the upper level times the probability that the
    atom makes the transition per second, in s^-1.  It is the quantity
    Cowan's atomic-structure codes print for every transition they compute,
    and it measures how strong the line is expected to be.
u%gA
    The uncertainty of gA quoted in the same file, in percent of gA.
u_ln
    That same uncertainty on a logarithmic scale, u_ln = ln(1 + u%gA/100).
    A logarithmic scale is the natural one here because the uncertainty of a
    calculated gA is a factor, not an amount: u_ln = 0.7 means "uncertain by
    a factor of about two", whichever decade gA lies in.
cutoff
    Cowan's codes printed a transition only when gA reached this value
    (1e3 s^-1 for this calculation).  A transition that does not appear in the
    file is therefore not a transition of unknown strength: it is one whose
    strength is *known to be below the cutoff*.  Statisticians call such data
    left-censored.  Treating it as unknown throws away that information and
    lets a candidate identification escape the intensity test entirely.

What this module computes
-------------------------
A single pair (gA_missing, u_ln_missing) to stand for every censored
transition, chosen so that

    gA_missing * exp(u_ln_missing) = cutoff,

that is: the upper edge of the one-standard-deviation interval of the imputed
value sits exactly on the censoring threshold.  This is the natural choice -
the imputed transition is as strong as it could be while still being
consistent with having been left out of the file.

u_ln_missing itself is read off the surviving data.  u_ln grows as gA falls
(weak transitions are the badly calculated ones), so the value appropriate
just below the cutoff is the value seen just above it.  Two ways of doing
that are offered:

  * the plain estimate (`self_consistent=False`): u_ln_missing is the root
    mean square of u_ln over the rows in the lowest window of the file,
    cutoff <= gA <= window*cutoff.  With the default window of 10 this is one
    decade of gA, the decade adjacent to the censored region.
  * the extrapolated estimate (`self_consistent=True`): u_ln is fitted as a
    straight line in log10(gA) over the decades named by `fit_range_decades`,
    and the pair of equations

        u_ln_missing = fit(log10(gA_missing))
        gA_missing   = cutoff * exp(-u_ln_missing)

    is solved together, so that the uncertainty used is the one belonging to
    the imputed gA itself rather than to the decade above it.  The solution is
    found by iteration; the fit has a negative slope of small magnitude, so
    the iteration is a contraction and converges geometrically.

The plain estimate is the conservative one (a larger gA_missing, hence a
smaller penalty on the identification); the extrapolated one is the more
faithful to the trend of the data.

The second half of the module turns an imputed gA into an imputed intensity.
The intensities in the file follow

    Icalc = C * gA * exp(-Eup/kT) / rwn

exactly (to eight significant figures over all 30206 rows), where Eup is the
energy of the upper level of the transition and rwn its Ritz wavenumber, both
in cm^-1, and kT is an effective excitation temperature expressed as an
energy in the same unit.  `fit_intensity_model` recovers C and kT from the
file, and `impute_intensity` applies the relation.

Run `python gA_imputation.py --report` for the calibration report.
"""
import argparse
import math
import os
import sys

import numpy as np
import pandas as pd

import config

# Logical column names this module needs from the calculated-transition file.
NEEDED = ('gA', 'u_ln', 'Icalc', 'Eup', 'rwn')


class ImputationError(Exception):
    """Raised when the data cannot support the estimate that was asked for."""


# --- reading ----------------------------------------------------------------

def load_icalc(cfg=None) -> pd.DataFrame:
    """The calculated-transition file as a frame with logical column names.

    Columns are found by the names given in the configuration, then renamed to
    the logical names used here, so that nothing downstream depends on the
    spelling used in the workbook.
    """
    cfg = cfg or config.load()
    df = pd.read_excel(cfg.icalc_file, sheet_name=cfg.icalc.sheet)
    have = {str(c).strip(): c for c in df.columns}
    rename = {}
    for logical in NEEDED:
        physical = cfg.icalc.columns[logical]
        if physical not in have:
            raise ImputationError(
                f"{os.path.basename(cfg.icalc_file)}[{cfg.icalc.sheet}]: no "
                f"column named {physical!r} (wanted for {logical!r}); the "
                f"header row has: " + ', '.join(repr(h) for h in have))
        rename[have[physical]] = logical
    return df.rename(columns=rename)[list(NEEDED)].astype(float)


# --- the uncertainty just above the cutoff ----------------------------------

def _estimate(values: np.ndarray, estimator: str) -> float:
    """One summary number for a set of u_ln values."""
    if estimator == 'rms':
        return float(np.sqrt(np.mean(values ** 2)))
    if estimator == 'mean':
        return float(np.mean(values))
    if estimator == 'median':
        return float(np.median(values))
    raise ImputationError(f"unknown estimator {estimator!r}; "
                          "expected 'rms', 'mean' or 'median'")


def window_u_ln(df: pd.DataFrame, cutoff: float, window: float,
                estimator: str, min_rows: int = 30) -> float:
    """u_ln summarised over the rows with cutoff <= gA <= window*cutoff."""
    gA, u_ln = df['gA'].to_numpy(float), df['u_ln'].to_numpy(float)
    sel = (gA >= cutoff) & (gA <= window * cutoff) & np.isfinite(u_ln)
    if int(sel.sum()) < min_rows:
        raise ImputationError(
            f"only {int(sel.sum())} rows with {cutoff:g} <= gA <= "
            f"{window * cutoff:g}; too few to estimate u_ln from")
    return _estimate(u_ln[sel], estimator)


def decade_profile(df: pd.DataFrame, estimator: str,
                   lo: float = None, hi: float = None, min_rows: int = 30):
    """u_ln summarised decade by decade of gA.

    Returns three parallel lists: the centre of each decade in log10(gA), the
    summary of u_ln in it, and the number of rows it holds.  `lo` and `hi`
    bound log10(gA); by default every populated decade of the file is used.
    """
    gA, u_ln = df['gA'].to_numpy(float), df['u_ln'].to_numpy(float)
    ok = np.isfinite(gA) & (gA > 0) & np.isfinite(u_ln)
    x = np.log10(gA[ok])
    u_ln = u_ln[ok]
    lo = math.floor(x.min()) if lo is None else lo
    hi = math.ceil(x.max()) if hi is None else hi
    centres, values, counts = [], [], []
    edge = lo
    while edge < hi - 1e-9:
        sel = (x >= edge) & (x < edge + 1.0)
        n = int(sel.sum())
        if n >= min_rows:
            centres.append(edge + 0.5)
            values.append(_estimate(u_ln[sel], estimator))
            counts.append(n)
        edge += 1.0
    return centres, values, counts


def fit_u_ln(df: pd.DataFrame, estimator: str, fit_range_decades):
    """Straight-line fit of u_ln against log10(gA), decade by decade.

    Returns (slope, intercept) of  u_ln = slope*log10(gA) + intercept, fitted
    to one point per decade of gA within `fit_range_decades`.  Fitting the
    decade summaries rather than the individual rows keeps the estimator
    ('rms', 'mean' or 'median') meaningful: a least-squares fit to the rows
    themselves would always return the mean, whatever was asked for.
    """
    lo, hi = float(fit_range_decades[0]), float(fit_range_decades[1])
    centres, values, _ = decade_profile(df, estimator, lo, hi)
    if len(centres) < 2:
        raise ImputationError(
            f"only {len(centres)} populated decades between 10^{lo:g} and "
            f"10^{hi:g}; at least two are needed for the fit")
    slope, intercept = np.polyfit(centres, values, 1)
    return float(slope), float(intercept)


def estimate_missing_gA(df: pd.DataFrame, cutoff: float, window: float = 10.0,
                        estimator: str = 'rms', self_consistent: bool = False,
                        fit_range_decades=(3.0, 6.0), tol: float = 1e-14,
                        max_iter: int = 200):
    """The gA and the u_ln to stand for a transition absent from the file.

    `df` needs the columns 'gA' and 'u_ln' (see `load_icalc`).  Rows with
    gA below `cutoff` play no part: they are the handful that slipped below
    the printing threshold when the wavenumbers were rescaled, and they are
    not a sample of the censored population.

    Returns (gA_missing, u_ln_missing), always satisfying
    gA_missing*exp(u_ln_missing) = cutoff.

    Why the upper one-standard-deviation bound, and not something else
    -----------------------------------------------------------------
    All that is known about a censored transition is "gA is somewhere below
    the cutoff".  The imputed value is placed so that the top of its
    one-standard-deviation interval, gA_missing*exp(u_ln_missing), lands
    exactly on the cutoff: the transition is made as strong as it could be
    while still being consistent with having been left out of the file.  Any
    larger value would contradict the censoring; any smaller one would assert
    knowledge that is not there, and would punish the identification harder
    than the data warrant.

    WARNING - an earlier justification of this same formula was wrong, and the
    wrong version must not be reinstated.  It read: "subtract the uncertainty
    from gA and require the result to be zero".  On a logarithmic uncertainty
    the standard deviation in linear units is sigma = gA*(e^u_ln - 1), so

        gA - sigma = gA*(2 - e^u_ln),

    which vanishes only for the single value u_ln = ln 2 = 0.693, not for the
    1.738 measured here.  The reading is unsound in principle as well: gA is
    log-normally distributed, so it cannot reach zero at any finite number of
    standard deviations.  The defining relation implemented above is
    gA_missing*exp(u_ln_missing) = cutoff and nothing else.  For orientation,
    the lower one-standard-deviation bound of the imputed value is
    cutoff*exp(-2*u_ln_missing), about 3% of the cutoff.
    """
    if not self_consistent:
        u = window_u_ln(df, cutoff, window, estimator)
        return cutoff * math.exp(-u), u

    slope, intercept = fit_u_ln(df, estimator, fit_range_decades)
    if abs(slope) >= math.log(10.0):
        raise ImputationError(
            f"u_ln falls by {abs(slope):.3f} per decade of gA, too steeply "
            "for the self-consistent estimate to have a unique solution")
    # Iterate  log10(gA) = log10(cutoff) - u_ln(log10 gA)/ln(10),  starting at
    # the cutoff itself.  |slope|/ln(10) < 1 makes this a contraction.
    lg_cut = math.log10(cutoff)
    lg = lg_cut
    for _ in range(max_iter):
        nxt = lg_cut - (slope * lg + intercept) / math.log(10.0)
        if abs(nxt - lg) < tol:
            lg = nxt
            break
        lg = nxt
    else:
        raise ImputationError('the self-consistent estimate did not converge')
    u = slope * lg + intercept
    # Return the gA implied by u exactly, so that the defining relation holds
    # to machine precision whichever branch produced u.
    return cutoff * math.exp(-u), u


# --- the intensity that goes with an imputed gA -----------------------------

def fit_intensity_model(df: pd.DataFrame):
    """Recover C and kT of  Icalc = C*gA*exp(-Eup/kT)/rwn  from the file.

    Taking logarithms turns the relation into a straight line,
    ln(Icalc*rwn/gA) = ln(C) - Eup/kT, which is fitted by least squares.
    Returns (C, kT), kT in cm^-1.
    """
    d = df[(df['Icalc'] > 0) & (df['gA'] > 0) & (df['rwn'] > 0)]
    if len(d) < 2:
        raise ImputationError('too few usable rows to fit the intensity model')
    y = np.log(d['Icalc'].to_numpy(float) * d['rwn'].to_numpy(float)
               / d['gA'].to_numpy(float))
    slope, intercept = np.polyfit(d['Eup'].to_numpy(float), y, 1)
    if slope >= 0:
        raise ImputationError('the fitted excitation temperature is not '
                              'positive; the file does not follow the model')
    return float(np.exp(intercept)), float(-1.0 / slope)


def impute_intensity(gA, Eup, rwn, C: float, kT: float):
    """The intensity the model gives for a transition of strength `gA`."""
    return C * np.asarray(gA, float) * np.exp(-np.asarray(Eup, float) / kT) \
        / np.asarray(rwn, float)


def intensity_model_error(df: pd.DataFrame, C: float, kT: float) -> float:
    """Largest relative departure of the file from the model, as a fraction."""
    d = df[(df['Icalc'] > 0) & (df['gA'] > 0) & (df['rwn'] > 0)]
    pred = impute_intensity(d['gA'], d['Eup'], d['rwn'], C, kT)
    return float(np.max(np.abs(pred / d['Icalc'].to_numpy(float) - 1.0)))


def check_intensity_model(df: pd.DataFrame, cfg, stream=sys.stderr):
    """Compare the configured C and kT with a fresh fit; warn on disagreement.

    The configured values are returned whatever the outcome: a re-fit is a
    diagnostic, not a silent substitution.  Returns (C, kT, fitted_C,
    fitted_kT, agrees).
    """
    C_cfg = float(cfg.intensity_model['C'])
    kT_cfg = float(cfg.intensity_model['kT'])
    tol = float(cfg.intensity_model['verify_tolerance'])
    C_fit, kT_fit = fit_intensity_model(df)
    dC = abs(C_fit - C_cfg) / C_cfg
    dkT = abs(kT_fit - kT_cfg) / kT_cfg
    agrees = dC <= tol and dkT <= tol
    if not agrees and stream is not None:
        print(f"WARNING: the intensity model fitted to "
              f"{os.path.basename(cfg.icalc_file)} disagrees with the "
              f"configuration by more than {tol:.3%}:\n"
              f"    C  = {C_cfg:.6g} in the configuration, {C_fit:.6g} fitted "
              f"({dC:.3%} apart)\n"
              f"    kT = {kT_cfg:.6g} in the configuration, {kT_fit:.6g} "
              f"fitted ({dkT:.3%} apart)\n"
              "    The configured values are being used; update "
              "[intensity_model] if the fitted ones are right.",
              file=stream)
    return C_cfg, kT_cfg, C_fit, kT_fit, agrees


# --- report -----------------------------------------------------------------

def report(cfg=None, stream=sys.stdout) -> None:
    """Print the calibration of the imputation on the real file."""
    cfg = cfg or config.load()
    df = load_icalc(cfg)
    cutoff = cfg.gA_cutoff
    mg = cfg.missing_gA

    def p(*a):
        print(*a, file=stream)

    p(f"Calculated transitions: {os.path.basename(cfg.icalc_file)}"
      f"[{cfg.icalc.sheet}], {len(df)} rows")
    p(f"Printing cutoff of the calculation: gA >= {cutoff:g} s^-1"
      f"  ({int((df['gA'] < cutoff).sum())} rows fall below it)")
    p('')

    p("u_ln = ln(1 + u%gA/100), the uncertainty of gA as a factor, by decade:")
    p(f"    {'decade of gA':>16}  {'rows':>6}  {'rms':>7}  {'mean':>7}"
      f"  {'median':>7}")
    centres, rms, counts = decade_profile(df, 'rms')
    _, means, _ = decade_profile(df, 'mean')
    _, meds, _ = decade_profile(df, 'median')
    for c, n, r, m, d in zip(centres, counts, rms, means, meds):
        lo = int(round(c - 0.5))
        label = f"10^{lo}..10^{lo + 1}"
        p(f"    {label:>16}  {n:6d}  {r:7.4f}  {m:7.4f}  {d:7.4f}")
    p('')

    est = mg['u_ln_estimator']
    win = float(mg['u_ln_window'])
    gA0, u0 = estimate_missing_gA(df, cutoff, win, est, False)
    p(f"Plain estimate ({est} of u_ln over {cutoff:g} <= gA <= "
      f"{win * cutoff:g}):")
    p(f"    u_ln_missing = {u0:.4f}")
    p(f"    gA_missing   = {cutoff:g} * exp(-{u0:.4f}) = {gA0:.1f} s^-1")
    p('')

    lo, hi = mg['fit_range_decades']
    p(f"Extrapolated estimate (self-consistent), from a straight-line fit of "
      f"the\nper-decade summary against log10(gA) over "
      f"10^{lo:g}..10^{hi:g}:")
    for e in ('rms', 'mean', 'median'):
        try:
            slope, intercept = fit_u_ln(df, e, (lo, hi))
            gA1, u1 = estimate_missing_gA(df, cutoff, win, e, True, (lo, hi))
        except ImputationError as exc:
            p(f"    {e:6s}  not available: {exc}")
            continue
        mark = ' <- configured' if e == est else ''
        p(f"    {e:6s}  u_ln = {slope:+.4f}*log10(gA) {intercept:+.4f}"
          f"   ->  gA_missing = {gA1:7.1f} s^-1,  u_ln_missing = {u1:.4f}"
          f"{mark}")
    p('')

    C, kT, C_fit, kT_fit, agrees = check_intensity_model(df, cfg, stream=None)
    p("Intensity model  Icalc = C * gA * exp(-Eup/kT) / rwn  "
      "(Eup, rwn, kT in cm^-1):")
    p(f"    configured  C = {C:.6g}   kT = {kT:.6g}")
    p(f"    fitted      C = {C_fit:.6g}   kT = {kT_fit:.6g}"
      f"   ({'agree' if agrees else 'DISAGREE'} within "
      f"{cfg.intensity_model['verify_tolerance']:.3%})")
    p(f"    largest relative departure of the file from the fit: "
      f"{intensity_model_error(df, C_fit, kT_fit):.2e}")
    p('')

    gA_use, u_use = estimate_missing_gA(df, cutoff, win, est,
                                        bool(mg['self_consistent']), (lo, hi))
    p(f"With the configured settings (estimator={est!r}, window={win:g}, "
      f"self_consistent={bool(mg['self_consistent'])}):")
    p(f"    gA_missing = {gA_use:.1f} s^-1,  u_ln_missing = {u_use:.4f}")
    p("    imputed intensity for a few upper levels and wavenumbers:")
    p(f"    {'Eup (cm^-1)':>12}  {'rwn (cm^-1)':>12}  {'Icalc imputed':>14}")
    for Eup in (20000.0, 40000.0, 60000.0):
        for rwn in (10000.0, 30000.0):
            p(f"    {Eup:12.0f}  {rwn:12.0f}  "
              f"{float(impute_intensity(gA_use, Eup, rwn, C, kT)):14.4g}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description='Calibrate the value given to a transition that is absent '
                    'from the calculated-intensity file.')
    ap.add_argument('--config', default=config.DEFAULT_PATH, metavar='FILE',
                    help='configuration file (default: %(default)s)')
    ap.add_argument('--report', action='store_true',
                    help='print the calibration report')
    args = ap.parse_args(argv)
    cfg = config.load(args.config)
    if args.report:
        report(cfg)
    else:
        df = load_icalc(cfg)
        mg = cfg.missing_gA
        gA, u = estimate_missing_gA(df, cfg.gA_cutoff,
                                    float(mg['u_ln_window']),
                                    mg['u_ln_estimator'],
                                    bool(mg['self_consistent']),
                                    mg['fit_range_decades'])
        print(f"gA_missing = {gA:.4f} s^-1   u_ln_missing = {u:.6f}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
