#!/usr/bin/env python
"""Measure epsilon(lambda), the residual obscuration rate of Sugar's plates.

Run from inside LineClass/:

    python tools/obscuration_rate.py


1.  What is being measured, and why it is not the coverage function
-------------------------------------------------------------------
A predicted transition that theory says should be strong, and that the
observed line list does not contain, is evidence against whatever
identification predicted it - but only if a line COULD have been recorded
where it falls.  Three things stop that, and none of them has anything to do
with the atom: the wavelength fell between two photographic exposures; a
defect of the emulsion sat on it; a line of another species sat on it.

`tools/coverage_map.py` measures the first of those, and the large instances
of the other two, by asking where the line list is EMPTY over a stretch of
wavelength.  It has to work on stretches, because that is the only way a
statistical model can tell a hole from a fluctuation: one empty bin proves
nothing, a hundred consecutive ones prove a great deal.  Its output is the
coverage function

    c(lambda) = P(a line bright enough to be recorded at vacuum wavelength
                  lambda was in fact recorded), relative to a typical part
                  of the spectrum,

and it is blind, by construction, to obscuration that acts on ONE line: a
grain flaw a few tenths of an angstrom across, a single impurity line, a
scratch narrower than the model's own resolution.  Such events do not thin
out a stretch of the record enough to be detected as a stretch.  They are
also, individually, invisible: nothing in the line list marks the place where
a line should have been and was not.

What they do have is a RATE, and the rate is measurable.  Call it

    epsilon(lambda) = P(a line that was certainly bright enough to be
                        recorded, at a wavelength the coverage map calls
                        covered, was nevertheless not recorded),

so that the full observability of a predicted transition is

    P(recorded) = c(lambda) * (1 - epsilon(lambda)) * D(I_pred / N(lambda))

with N the noise level fitted by `tools/estimate_snr.py` and D the
probability that a line of that predicted strength clears it.  c is measured
where the record is thin; epsilon is measured where it is not.  Neither can
substitute for the other, and epsilon is the one the likelihood needs most,
because a level's strongest predicted branch almost always falls in a fully
covered region - of the 14 levels whose questionable mark rests on an absent
strongest branch, every one has that branch at c >= 0.997.


2.  The measurement
-------------------
The rate cannot be measured on the levels under test, because for those the
question of whether a predicted line is really missing is exactly what is in
dispute.  It is measured instead on the ESTABLISHED levels: those with more
old than new supporting identifications and present in the ASD compilation
(384 of the 594 in the current run).  Their energies are certain to a few
thousandths of a wavenumber, so their predicted transitions fall at known
places, and an absence there is the plate's doing rather than the
identification's.

For every predicted transition between two established levels the tool asks:

  * WHERE would it fall?  The Ritz wavelength lambda = 1e8 / (E_upper -
    E_lower) in vacuum angstroms, from the run's optimized energies.

  * HOW BRIGHT would it be?  I_pred, the calculated intensity from
    Icalc.xlsx on the pipeline's linear scale, against the noise level at
    that wavelength - the linear intensity corresponding to Sugar's plate
    intensity 1, which `level_shifts.read_intensity_calibration` reconstructs
    from the intensity-correction polynomials.  Their ratio r = I_pred /
    N(lambda) is the only strength measure used; it already absorbs the
    wavelength dependence of the plate's sensitivity.

  * WAS ANYTHING RECORDED THERE?  Three outcomes:

      recorded - a measured line lies within the matching window
                 W = 5.5 * sigma_lambda of the Ritz wavelength, where
                 sigma_lambda is the local median wavelength uncertainty of
                 the measured lines (0.003 - 0.009 A over the whole range; in
                 wavenumber terms 0.5 cm^-1 at 900 A and 0.007 cm^-1 at
                 1 um).  5.5 is the pipeline's own matching factor.  The line
                 need not be classified as this transition, or as Pr III at
                 all: the question is whether the plate registered anything
                 at that place.

      masked   - nothing within W, but a recorded line close enough and
                 strong enough to swallow the prediction, by the same test
                 `level_shifts` uses for its own masking check: within one
                 effective line width (the instrumental 0.035 A convolved
                 with the Doppler width) a line of at least half the
                 predicted intensity hides it, and beyond one width a much
                 stronger line hides it under its Gaussian wing.

      absent   - neither.

    "recorded" and "masked" together are counted as EXPLAINED; epsilon is
    measured on the rest.

  * WOULD ANYTHING HAVE BEEN RECORDED THERE ANYWAY?  This is the correction
    the measurement lives or dies by.  The matching window is 0.02 A wide and
    the recorded lines are, in the crowded ultraviolet, 0.13 A apart, so about
    one predicted position in eight lands on a recorded line by pure
    coincidence.  Left uncorrected that would hide a fifth of the obscuration.
    The chance rate is measured, not modelled: every prediction is re-tested
    at 18 CONTROL POSITIONS, displaced by +-6, +-9, ... +-30 window widths,
    and the fraction of those that come out explained is that prediction's
    chance rate f.  Displacing in units of the window keeps the controls
    inside the same local line density at every wavelength.  (As a check, the
    Poisson estimate 1 - exp(-2*W*rho) from the local line density rho agrees
    with the measured f to about 0.005 in every wavelength band.)

    With N predictions, y_i = 1 if explained, and c_i the coverage:

        epsilon = 1 - (sum y_i - sum f_i) / sum (1 - f_i) * c_i

    and its confidence interval comes from the profile likelihood of the
    binomial model p_i = f_i + (1 - f_i) * c_i * (1 - epsilon).


3.  Why the measurement is confined to the strongest predictions
----------------------------------------------------------------
The absence of a WEAK predicted line means nothing: it may simply have been
too faint.  Absence becomes informative only where D, the probability of
clearing the noise, is 1, and the honest way to find that place is not to
model D but to look at where the measured absence rate stops falling:

    I_pred / noise      predictions   epsilon
    10^1.5 - 10^2.0         526        0.110
    10^2.0 - 10^2.5         232        0.028
    10^2.5 - 10^3.0          73        0.058
    10^3.0 - 10^3.5          27        0.037

It stops at a hundred times the noise level, and everything above that is one
population.  That cut is the tool's `--strong` option.

A hundred times the noise sounds extravagant, and the reason it is needed is
worth stating, because it is a result in its own right.  If one models D as a
log-normal - the predicted intensity times a scatter factor whose spread is
measured on the accepted identifications, sd of ln(I_obs*BF / I_pred) = 1.30
over 3689 established-level lines, with BF the branching fraction that splits
a blended feature among its components - then D reaches 0.99 already at ten
times the noise, and the absence rate there should be 1 - 0.99 = 0.01.  It is
in fact 0.14.  The log-normal fitted to the accepted lines is far too
optimistic at the faint end, for the obvious reason: it is fitted to lines
that WERE recorded, so its lower tail has been cut off by the very effect it
is being asked to predict.  Fitting the tail and epsilon together does not
help - the two are degenerate, and a free fit puts the sd at 2.4 and epsilon
at 0, attributing every absence to the tail.  What breaks the degeneracy is
that at a hundred times the noise a log-normal of ANY plausible width predicts
essentially no loss, so whatever is left there is not the intensity model.

The consequence for the likelihood is that the missing-line term must use a
detection curve MEASURED this way (the D(z) table this tool writes) and not a
log-normal around the noise level.


4.  The result
--------------
Of the 335 predicted transitions between established levels that are at least
a hundred times the noise level and fall where the coverage map reports
c >= 0.95:  317 recorded, 7 masked, 11 absent, against a chance-match rate of
0.116.  That gives

    epsilon = 0.035,  95 per cent interval 0.018 - 0.061

and the rate is the same everywhere: per wavelength band the values run
0.000 - 0.072 and a single constant fits them (likelihood ratio 7.8 on 5
degrees of freedom, p = 0.17).  Roughly one predicted line in thirty is lost
to something too small for the coverage map to see, wherever in the spectrum
it falls.

Do not read the per-band numbers off a fit that uses the weak predictions as
well.  Such a fit has to assume the detection curve has the same shape in
every band, and it does not - the noise level is a fitted polynomial, and a
small error in it at one end shifts the whole curve there.  Fitted that way
the 7000-11000 A band comes out at 0.26 against 0.072 measured on its own
strong predictions.  The strong subsample is the measurement; everything else
is the diagnostic that justifies it.


5.  What epsilon is not
-----------------------
It is a RATE, not a map: it says how often a line is lost, never where.  Used
in the likelihood it weakens the evidence of one absent branch by the factor
1 - epsilon = 0.965, which is small for a single line and decisive only when
several of a level's branches are absent at once - and there it has to be
combined with the correlation length of the obscuration, which is measured by
tools/obscuration_length.py and comes out at L = 0.067 A (0.027 - 0.132), the
width of a resolution element.  Two predictions closer than that are hidden by
one grain flaw and their joint absence costs epsilon rather than epsilon
squared; farther apart than that they are independent, and since the closest
two branches of any established level ever come, among predictions strong
enough for epsilon to apply to them, is 0.446 A, the correction is in practice
never needed.

Unlike tools/coverage_map.py this is not a standalone tool: it reads the run
through level_shifts.py and chance_mc.py, because it needs the optimized level
energies, the old/new status of every level and the calculated intensities.

Outputs (in the working directory unless --out-dir is given):
  obscuration_rate.txt   the headline number, the band table, the detection
                         curve, and the list of absent strong predictions
  obscuration_rate.csv   one row per prediction, with every quantity above
  obscuration_rate.png   the detection curve, epsilon per band, and where the
                         absent strong predictions fall
  obscuration_rate.log   the console output
"""
import argparse
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import chance_mc as mc                      # noqa: E402
import level_shifts as ls                   # noqa: E402
import lopt_lines                           # noqa: E402


# --- the population ---------------------------------------------------------
def established_levels(per: pd.DataFrame) -> set:
    """Level ids whose position is not in question.

    Two conditions, both taken from the run itself: the level is supported by
    more old than new identifications (new_star = 0 in level_shifts' sense),
    and it is present in the ASD compilation (is_new_level = 0).  A level that
    fails either is a level whose energy the current work may still move, and
    a prediction made from a moving energy cannot be asked whether it landed
    on a recorded line.
    """
    return set(per.loc[(per['new_star'] == 0) & (per['is_new_level'] == 0),
                       'level_id'])


def build_predictions(opt, log):
    """Every predicted transition between two established levels.

    The returned frame carries the Ritz wavelength lam (vacuum angstroms), the
    calculated intensity I_pred on the pipeline's linear scale, the noise level
    thr at that wavelength on the same scale, and z = ln(I_pred / thr).
    """
    if opt.lopt:
        ls.LOPT_LINES = opt.lopt
        ls.LOPT_LEVELS = opt.lopt_levels
        ls.ENERGIES = opt.energies
    levels, real, _ = ls.read_run()
    per = mc.per_level_table(real, levels)
    per['new_star'] = (per['n_new'] > per['n_old']).astype(int)
    e_final = dict(zip(per['level_id'], per['E_final']))
    est = established_levels(per)
    log(f"levels in the run: {len(per)}; established (old and in ASD): "
        f"{len(est)}")

    preds = ls.load_predictions(set(per['level_id']), e_final)
    keep = preds['lo_id'].isin(est) & preds['up_id'].isin(est)
    P = preds[keep].reset_index(drop=True).copy()
    log(f"predicted transitions: {len(preds)} in all, {len(P)} between two "
        f"established levels")

    calib = ls.read_intensity_calibration()
    if calib is None:
        raise SystemExit("intensity_correction_functions.txt not found - the "
                         "noise level cannot be reconstructed")
    P['wn'] = [e_final[u] - e_final[l]
               for l, u, _ in P.itertuples(index=False)]
    P['lam'] = 1.0e8 / P['wn']
    P['thr'] = [ls.noise_threshold_linear(w, calib) for w in P['wn']]
    P = P[P['thr'].notna() & (P['I_pred'] > 0)].reset_index(drop=True)
    P['z'] = np.log(P['I_pred'] / P['thr'])
    return P, preds, real, est


# --- where a line could have been seen --------------------------------------
def matching_window(P: pd.DataFrame, obs: pd.DataFrame, k: float, log):
    """Half-width W of the position-matching window, in angstroms.

    The measured lines carry a wavenumber uncertainty; converted to wavelength
    (sigma_lambda = sigma_wn * lambda^2 / 1e8) it is nearly constant over the
    whole range, 0.003 - 0.009 A.  A predicted line has no uncertainty of its
    own to quote, so it borrows the local one: the running median over 101
    neighbouring measured lines.  W = k * sigma_lambda with k = 5.5, the factor
    classify_lines.py uses for its own candidate matching.
    """
    o = obs.sort_values('lam').reset_index(drop=True)
    u = (o['unc_wn_obs'] * o['lam'] ** 2 * 1.0e-8).to_numpy()
    s = (pd.Series(u).rolling(101, center=True, min_periods=11)
         .median())
    s = s.ffill().bfill().to_numpy()
    lam_o = o['lam'].to_numpy()
    idx = np.clip(np.searchsorted(lam_o, P['lam'].to_numpy()), 0, len(s) - 1)
    W = k * s[idx]
    log(f"matching window: {k} * sigma_lambda, median half-width "
        f"{np.median(W):.4f} A ({np.min(W):.4f} - {np.max(W):.4f})")
    return W


def local_density(P: pd.DataFrame, lam_o: np.ndarray, k: int = 41):
    """Recorded lines per angstrom near each prediction, from the k nearest."""
    i = np.searchsorted(lam_o, P['lam'].to_numpy())
    lo = np.clip(i - k // 2, 0, len(lam_o) - k)
    return k / (lam_o[lo + k - 1] - lam_o[lo])


def is_recorded(pos: np.ndarray, W: np.ndarray, lam_o: np.ndarray):
    """True where a measured line lies within W of the position."""
    i = np.searchsorted(lam_o, pos)
    d = np.full(len(pos), np.inf)
    for j in (-1, 0):
        k = np.clip(i + j, 0, len(lam_o) - 1)
        d = np.minimum(d, np.abs(lam_o[k] - pos))
    return d <= W


def is_masked(pos: np.ndarray, i_pred: np.ndarray, wn_o, int_o, dop):
    """True where a recorded line is close and strong enough to swallow the
    prediction, by level_shifts' own masking test."""
    return np.array([ls._hidden_by_a_stronger_line(1.0e8 / p, v, wn_o, int_o,
                                                   dop)
                     for p, v in zip(pos, i_pred)], dtype=bool)


def chance_rate(P, W, lam_o, wn_o, int_o, dop, offsets, log):
    """Probability that a position is explained by coincidence alone.

    Every prediction is re-tested at control positions displaced by the given
    multiples of its own window.  The controls sample the same local line
    density and the same masking environment, so their explained fraction is
    the chance rate the real position faced.

    Returns (f_explained, f_recorded): the second counts only the control
    positions that landed on a recorded line, which is the part the Poisson
    estimate from the local line density can be compared with.
    """
    lam = P['lam'].to_numpy()
    ip = P['I_pred'].to_numpy()
    hits = np.zeros(len(P))
    hits_rec = np.zeros(len(P))
    for m in offsets:
        pos = lam + m * W
        rec = is_recorded(pos, W, lam_o)
        hits_rec += rec
        hits += (rec | is_masked(pos, ip, wn_o, int_o, dop))
    f = hits / len(offsets)
    log(f"chance rate from {len(offsets)} control positions per prediction: "
        f"mean {f.mean():.4f} explained, of which "
        f"{hits_rec.sum() / hits.sum():.0%} by landing on a recorded line")
    return f, hits_rec / len(offsets)


def attach_coverage(P: pd.DataFrame, path: str, log):
    """c(lambda) at each prediction, interpolated from coverage_map.csv."""
    if not os.path.exists(path):
        log(f"WARNING: {os.path.basename(path)} not found - coverage taken as "
            f"1 everywhere, so epsilon will absorb the mapped blind stretches "
            f"as well")
        return np.ones(len(P))
    cov = pd.read_csv(path)
    mid = 0.5 * (cov['lam_lo'] + cov['lam_hi']).to_numpy()
    return np.interp(P['lam'].to_numpy(), mid, cov['coverage'].to_numpy())


# --- the estimator ----------------------------------------------------------
def nll(eps: float, y, f, c) -> float:
    """Negative log-likelihood of the binomial model
    p_i = f_i + (1 - f_i) * c_i * (1 - eps)."""
    p = np.clip(f + (1.0 - f) * c * (1.0 - eps), 1e-12, 1.0 - 1e-12)
    return -float((y * np.log(p) + (1.0 - y) * np.log(1.0 - p)).sum())


def fit_eps(T: pd.DataFrame):
    """(eps, low, high, nll) - maximum likelihood and its 95 per cent profile
    interval, the values where the log-likelihood has fallen by 1.92."""
    from scipy.optimize import brentq, minimize_scalar
    y = T['explained'].to_numpy().astype(float)
    f = T['f_chance'].to_numpy()
    c = T['c'].to_numpy()
    r = minimize_scalar(lambda e: nll(e, y, f, c), bounds=(0.0, 0.95),
                        method='bounded', options=dict(xatol=1e-8))
    e, f0 = float(r.x), float(r.fun)

    def drop(x):
        return nll(x, y, f, c) - f0 - 1.92

    lo = brentq(drop, 0.0, e) if drop(0.0) > 0 else 0.0
    hi = brentq(drop, e, 0.95) if drop(0.95) > 0 else 0.95
    return e, lo, hi, f0


def detection_curve(P: pd.DataFrame, edges):
    """Chance- and coverage-corrected fraction explained, by strength bin.

    This is D(z) of the module docstring, measured rather than modelled: the
    probability that a predicted transition of a given strength relative to the
    noise level was recorded, with coincidences removed and the mapped blind
    stretches divided out.  Its plateau is 1 - epsilon.
    """
    rows = []
    for a, b in zip(edges[:-1], edges[1:]):
        T = P[(P['z'] >= a) & (P['z'] < b)]
        if len(T) < 10:
            continue
        y = T['explained'].to_numpy().astype(float)
        f = T['f_chance'].to_numpy()
        c = T['c'].to_numpy()
        den = ((1.0 - f) * c).sum()
        d = (y.sum() - f.sum()) / den
        se = math.sqrt(float(((y - f - (1.0 - f) * c * d) ** 2).sum())) / den
        rows.append((a, b, len(T), d, se))
    return rows


# --- diagnostics ------------------------------------------------------------
def intensity_residuals(real, preds, est, log):
    """sd of ln(I_obs * BF / I_pred) over the accepted lines of established
    pairs - the width of the intensity model that section 3 of the module
    docstring shows to be too optimistic at the faint end.

    BF is the branching fraction: the share of a blended feature's intensity
    that belongs to this one transition.  Comparing a whole blended feature
    with one component's prediction would put the residual above the scale by
    the reciprocal of that share.
    """
    ip = {tuple(sorted((l, u))): v
          for l, u, v in preds.itertuples(index=False)}
    acc = real[real['accepted'] == 1]
    bf = (pd.to_numeric(acc['BF'], errors='coerce').fillna(1.0)
          if 'BF' in acc.columns else pd.Series(1.0, index=acc.index))
    lam, res = [], []
    for lo, up, wn, io, b in zip(acc['low_id'], acc['upp_id'], acc['wn_obs'],
                                 acc['obs_intens'], bf):
        lo, up = str(lo), str(up)
        if lo not in est or up not in est:
            continue
        v = ip.get(tuple(sorted((lo, up))))
        if v is None or not (v > 0) or not (io > 0):
            continue
        b = float(b) if float(b) > 0 else 1.0
        lam.append(1.0e8 / float(wn))
        res.append(math.log(float(io) * b / float(v)))
    if len(res) < 20:
        return None
    r = np.asarray(res)
    log(f"intensity model on {len(r)} accepted lines of established pairs: "
        f"median ln(I_obs*BF/I_pred) = {np.median(r):+.3f}, sd = "
        f"{r.std(ddof=1):.3f}")
    return pd.DataFrame({'lam': lam, 'res': res})


# --- output -----------------------------------------------------------------
def write_report(path, strong, res, curve, bands, opt, log):
    e, lo, hi, _ = fit_eps(strong)
    out = []
    w = out.append
    w("# residual obscuration rate of the plates, epsilon")
    w("# produced by tools/obscuration_rate.py")
    w("#")
    w("# epsilon = P(a predicted line certainly bright enough to be recorded,")
    w("#           at a wavelength the coverage map calls covered, was")
    w("#           nevertheless not recorded).")
    w("#")
    w("# population: predicted transitions between two established levels,")
    w(f"#             I_pred >= {opt.strong:g} x the noise level, coverage "
      f">= {opt.cmin:g}")
    w("")
    w(f"n_predictions   {len(strong)}")
    w(f"n_recorded      {int(strong['recorded'].sum())}")
    w(f"n_masked_only   {int((strong['masked'] & ~strong['recorded']).sum())}")
    w(f"n_absent        {int((~strong['explained']).sum())}")
    w(f"chance_rate     {strong['f_chance'].mean():.4f}")
    w(f"epsilon         {e:.4f}")
    w(f"epsilon_lo95    {lo:.4f}")
    w(f"epsilon_hi95    {hi:.4f}")
    if res is not None:
        w(f"intensity_sd    {res['res'].std(ddof=1):.4f}   "
          f"# sd of ln(I_obs*BF/I_pred); see section 3")
    w("")
    w("# epsilon by wavelength band, on the same strong subsample")
    w("# lam_lo   lam_hi      n  absent   epsilon    lo95    hi95")
    for a, b, n, na, be, blo, bhi in bands:
        w(f"{a:9.1f}{b:10.1f}{n:7d}{na:8d}{be:10.4f}{blo:8.4f}{bhi:8.4f}")
    w("")
    w("# the detection curve D(z), z = ln(I_pred / noise level): the chance-")
    w("# and coverage-corrected fraction of predictions recorded or masked.")
    w("# Its plateau is 1 - epsilon.  Use this, not a log-normal around the")
    w("# noise level, for the missing-line term of the likelihood.")
    w("#    z_lo     z_hi     I_pred/noise range        n       D      sd")
    for a, b, n, d, se in curve:
        z_lo = "-inf" if not np.isfinite(a) else f"{a:.2f}"
        z_hi = "+inf" if not np.isfinite(b) else f"{b:.2f}"
        rng = (f"{math.exp(a) if np.isfinite(a) else 0:.3g} - "
               f"{math.exp(b) if np.isfinite(b) else float('inf'):.3g}")
        w(f"{z_lo:>9s}{z_hi:>9s}{rng:>23s}{n:9d}{d:8.4f}{se:8.4f}")
    w("")
    w("# the absent strong predictions, strongest first")
    w("#   lambda      I_pred     noise   I/noise   coverage")
    M = strong[~strong['explained']].sort_values('I_pred', ascending=False)
    for r in M.itertuples(index=False):
        w(f"{r.lam:10.3f}{r.I_pred:12.4g}{r.thr:10.4g}"
          f"{r.I_pred / r.thr:10.4g}{r.c:11.4f}")
    with open(path, 'w', newline='\n') as fh:
        fh.write("\n".join(out) + "\n")
    log(f"wrote {os.path.basename(path)}")
    return e, lo, hi


def plot(path, strong, curve, bands, e, lo, hi, opt, log):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 1, figsize=(11, 12))

    a0 = ax[0]
    # the two end bins are open, so give them the width of their neighbours
    zc, d, se = [], [], []
    for a, b, n, dd, ss in curve:
        aa = b - 1.15 if not np.isfinite(a) else a
        bb = a + 1.15 if not np.isfinite(b) else b
        zc.append(0.5 * (aa + bb) / math.log(10))
        d.append(dd)
        se.append(ss)
    a0.errorbar(zc, d, yerr=se, fmt='o-', color='tab:blue',
                label='measured D')
    a0.axhline(1 - e, color='tab:red', ls='--',
               label=f'plateau 1 - eps = {1 - e:.3f}')
    a0.axvspan(math.log10(opt.strong), max(zc) + 0.6, color='0.9', zorder=0)
    a0.set_xlabel('log10( predicted intensity / noise level )')
    a0.set_ylabel('P(recorded or masked)')
    a0.set_title('the detection curve, measured on the established levels; '
                 'shaded = the subsample epsilon is measured on')
    a0.set_ylim(0, 1.05)
    a0.legend(loc='upper left')
    a0.grid(alpha=0.3)

    a1 = ax[1]
    xs = [math.sqrt(x[0] * x[1]) for x in bands]
    xe = [[math.sqrt(x[0] * x[1]) - x[0] for x in bands],
          [x[1] - math.sqrt(x[0] * x[1]) for x in bands]]
    ys = [x[4] for x in bands]
    ye = [[x[4] - x[5] for x in bands], [x[6] - x[4] for x in bands]]
    a1.errorbar(xs, ys, xerr=xe, yerr=ye, fmt='o', color='tab:blue')
    a1.axhline(e, color='tab:red', ls='--', label=f'epsilon = {e:.3f}')
    a1.axhspan(lo, hi, color='tab:red', alpha=0.12)
    a1.set_xscale('log')
    a1.set_xlabel('vacuum wavelength, A')
    a1.set_ylabel('epsilon')
    a1.set_title('epsilon by wavelength band, at matched predicted strength')
    a1.legend(loc='upper left')
    a1.grid(alpha=0.3)

    a2 = ax[2]
    a2.semilogx(strong['lam'], np.log10(strong['I_pred'] / strong['thr']),
                '.', ms=4, color='0.7', label='recorded or masked')
    M = strong[~strong['explained']]
    a2.semilogx(M['lam'], np.log10(M['I_pred'] / M['thr']), 'o', ms=8,
                mfc='none', color='tab:red', label='absent')
    a2.set_xlabel('vacuum wavelength, A')
    a2.set_ylabel('log10( I_pred / noise level )')
    a2.set_title(f'the {len(strong)} strong predictions, and the {len(M)} '
                 f'of them that are absent')
    a2.legend(loc='upper left')
    a2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    log(f"wrote {os.path.basename(path)}")


# --- driver -----------------------------------------------------------------
def parse_args(argv):
    p = argparse.ArgumentParser(
        prog='obscuration_rate.py',
        description=__doc__.split('\n')[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    g = p.add_argument_group('the run to measure')
    g.add_argument('--lopt', default=None,
                   help='LOPT line-output file; without it the '
                        'classify_lines.py run is used')
    g.add_argument('--lopt-levels', default=None,
                   help='LOPT level-output file, for --energies lopt')
    g.add_argument('--energies', default='refit', choices=('refit', 'lopt'),
                   help='which energies a --lopt run uses')
    g = p.add_argument_group('the population')
    g.add_argument('--strong', type=float, default=100.0,
                   help='I_pred / noise level above which an absence is '
                        'informative')
    g.add_argument('--cmin', type=float, default=0.95,
                   help='minimum coverage; below it the coverage map, not '
                        'epsilon, explains the absence')
    g.add_argument('--k-window', type=float, default=5.5,
                   help='matching window in units of the local wavelength '
                        'uncertainty')
    g.add_argument('--offsets', type=float, nargs='+',
                   default=[6, 9, 12, 15, 18, 21, 24, 27, 30],
                   help='control displacements, in window widths, used both '
                        'ways to measure the chance rate')
    g.add_argument('--bands', type=float, nargs='+',
                   default=[800, 1200, 1700, 2300, 3000, 4500, 7000, 11000],
                   help='wavelength band edges for epsilon(lambda), A')
    g.add_argument('--min-band', type=int, default=10,
                   help='predictions a band needs before it is reported')
    g = p.add_argument_group('input and output')
    g.add_argument('--coverage',
                   default=os.path.join(ROOT, 'coverage_map.csv'))
    g.add_argument('--out-dir', default=ROOT)
    g.add_argument('--prefix', default='obscuration_rate')
    g.add_argument('--no-plot', action='store_true')
    return p.parse_args(argv)


def main(argv=None):
    opt = parse_args(argv)
    lines = []

    def log(msg):
        print(msg)
        lines.append(str(msg))

    P, preds, real, est = build_predictions(opt, log)
    res = intensity_residuals(real, preds, est, log)

    # The complete observed-line workbook, not the run's candidate table: the
    # question here is whether the plate registered anything at a wavelength,
    # and a line that no candidate transition ever matched answers it just as
    # well as a classified one.
    obs = lopt_lines.read_observed_lines()
    obs['lam'] = 1.0e8 / obs['wn_obs']
    obs = obs.sort_values('wn_obs').reset_index(drop=True)
    wn_o = obs['wn_obs'].to_numpy()
    int_o = np.nan_to_num(obs['obs_intens'].to_numpy())
    lam_o = np.sort(1.0e8 / wn_o)
    log(f"measured lines: {len(obs)}, {lam_o[0]:.2f} - {lam_o[-1]:.2f} A")

    W = matching_window(P, obs, opt.k_window, log)
    P['W'] = W
    P['rho'] = local_density(P, lam_o)
    P['f_chance_poisson'] = 1.0 - np.exp(-2.0 * W * P['rho'])
    dop = ls.doppler_fwhm_factor()
    P['recorded'] = is_recorded(P['lam'].to_numpy(), W, lam_o)
    P['masked'] = is_masked(P['lam'].to_numpy(), P['I_pred'].to_numpy(),
                            wn_o, int_o, dop)
    P['explained'] = P['recorded'] | P['masked']
    offsets = np.concatenate([np.asarray(opt.offsets, dtype=float),
                              -np.asarray(opt.offsets, dtype=float)])
    P['f_chance'], P['f_chance_recorded'] = chance_rate(
        P, W, lam_o, wn_o, int_o, dop, offsets, log)
    log(f"  landing on a recorded line: measured "
        f"{P['f_chance_recorded'].mean():.4f}, Poisson estimate from the "
        f"local line density {P['f_chance_poisson'].mean():.4f} - a check on "
        f"the controls, not an input")
    P['c'] = attach_coverage(P, opt.coverage, log)

    curve = detection_curve(P[P['c'] >= opt.cmin],
                            [-np.inf, -2.3, -1.15, 0.0, 1.15, 2.3, 3.45, 4.6,
                             np.inf])
    strong = P[(P['c'] >= opt.cmin)
               & (P['I_pred'] >= opt.strong * P['thr'])].reset_index(drop=True)
    log(f"strong subsample: {len(strong)} predictions "
        f"(I_pred >= {opt.strong:g} x noise, coverage >= {opt.cmin:g})")

    bands, nll_bands = [], 0.0
    edges = list(opt.bands)
    for a, b in zip(edges[:-1], edges[1:]):
        T = strong[(strong['lam'] >= a) & (strong['lam'] < b)]
        if len(T) < opt.min_band:
            log(f"  band {a:.0f}-{b:.0f} A: only {len(T)} strong "
                f"prediction(s) - not reported")
            continue
        be, blo, bhi, f0 = fit_eps(T)
        nll_bands += f0
        bands.append((a, b, len(T), int((~T['explained']).sum()), be, blo,
                      bhi))

    e, lo, hi, _ = fit_eps(strong)
    log("")
    log(f"epsilon = {e:.4f}   95% interval {lo:.4f} - {hi:.4f}   "
        f"({int((~strong['explained']).sum())} absent of {len(strong)})")
    if len(bands) > 1:
        from scipy.stats import chi2
        used = pd.concat([strong[(strong['lam'] >= a) & (strong['lam'] < b)]
                          for a, b, *_ in bands])
        _, _, _, f_used = fit_eps(used)
        lr = 2.0 * (f_used - nll_bands)
        dof = len(bands) - 1
        log(f"constant across the bands: likelihood ratio {lr:.1f} on {dof} "
            f"degrees of freedom, p = {1 - chi2.cdf(lr, dof):.2f}")
    for a, b, n, na, be, blo, bhi in bands:
        log(f"  {a:6.0f} - {b:6.0f} A   n = {n:4d}   absent {na:3d}   "
            f"epsilon = {be:.3f}  ({blo:.3f} - {bhi:.3f})")

    os.makedirs(opt.out_dir, exist_ok=True)
    base = os.path.join(opt.out_dir, opt.prefix)
    P.to_csv(base + '.csv', index=False, lineterminator='\n')
    log(f"wrote {opt.prefix}.csv ({len(P)} predictions)")
    write_report(base + '.txt', strong, res, curve, bands, opt, log)
    if not opt.no_plot:
        plot(base + '.png', strong, curve, bands, e, lo, hi, opt, log)
    with open(base + '.log', 'w', newline='\n') as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == '__main__':
    main()
