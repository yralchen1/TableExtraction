#!/usr/bin/env python
"""How strongly the observed lines say that a level sits where it is put.

Run from inside LineClass/:

    python level_positions.py                       # every level, adopted E
    python level_positions.py --detail 059003.000271
    python level_positions.py --scan                # + the alternate-position scan
    python level_positions.py --audit               # + what to do about them


1.  The quantity
----------------
For a level L and a candidate energy E (cm^-1) this module computes

    ln R(L, E) = ln [ P(what was recorded | a real level of Pr III sits at E)
                    / P(what was recorded | nothing sits at E) ]

H1 says a real level sits at E and the lines assigned to its transitions are
those transitions.  H0 says nothing sits there: every observed line that
happens to lie near one of the predicted positions lies there by coincidence,
and every predicted transition with no line near it is absent because there was
never anything to see.

ln R is a sum over the level's PREDICTED transitions, not over its accepted
lines.  That one choice merges two of the four things a level has to answer
for - the improbability of the alignments (a prediction with a line on it) and
the penalty for absences (a prediction without one) - into a single formula,
because a prediction with no line near it is simply a prediction whose datum is
"no line".

The other two follow from the same expression without any term being added by
hand: the intensity test is a second factor inside the branch that says the
line is genuine, and the blend penalty is what happens when the observed
feature is already explained by transitions that have nothing to do with L.


2.  The per-prediction formula
------------------------------
Take one predicted transition t of L.  Its partner level M sits at the fixed
energy E_M, so the transition is expected at the Ritz wavenumber

    nu_t(E) = E - E_M   (L the upper level)   or   E_M - E   (L the lower)

- the two move in opposite directions as E is scanned, which is what makes the
maximum of ln R sharp rather than flat.  Write

    sigma_t   the position uncertainty of a line found there (section 4)
    W_t       = 4 sigma_t, the half-width of the window searched for a line
    P_t       = c(lambda) D(z), the probability that the line would have been
                recorded at all: c the coverage function measured by
                tools/coverage_map.py, D the detection curve measured by
                tools/obscuration_rate.py, both reached through
                level_shifts.observation_probability
    rho_t     the density of observed lines near nu_t, per cm^-1
    eta       the probability that a genuine recorded line sits at an
              anomalous position - a corrupted measurement, a misread blend -
              so that its position carries no information (section 4)
    G_t       the intensity likelihood ratio of section 3

Model the observed line list near nu_t as a Poisson process of unrelated lines
of rate rho_t, plus, under H1, at most one genuine line, present with
probability P_t and placed Gaussian about nu_t.  The Poisson background factor
is identical under the two hypotheses and cancels, leaving, with
p = P_t (1 - eta):

    a line found at residual d = wn_obs - nu_t:

        R_t = (1 - p) + p * N(d; 0, sigma_t) / rho_t * G_t

    no line inside the window:

        R_t = 1 - p (1 - q_out),     q_out = 2 Phi(-4) = 6.3e-5

The first branch of the matched case is the world in which the line is there by
chance even though the level is real; that is why the intensity ratio G_t
multiplies only the second branch, and why nothing here can be worse than
ln(1 - p) however wrong the position or the intensity turns out to be.  A
single corrupted measurement therefore cannot kill a real level, which matters:
the residuals of this run are measurably heavy-tailed, |t| > 4 occurring twenty
times as often as a Gaussian allows.

Two consequences worth stating in the other direction.  A well-centred line on
an open plate has N(0, sigma)/rho of about 25 for the typical sigma = 0.1 cm^-1
and rho = 0.04 lines per cm^-1 of this spectrum, so it is worth ln R_t = +3.2,
and ten of them e^32.  A missing prediction costs ln(1 - p), which is bounded,
small where the plate is blind or the line is at the noise, and near ln(eta)
only for a prediction that was certain to be seen.


3.  The intensity ratio G_t
---------------------------
A feature can be free - nothing else accepted on it - or claimed, meaning that
other accepted transitions, between levels that do not involve L, already
account for it.  Write I_obs for its measured intensity, f(lambda) for the
far-ultraviolet scale correction of the predicted intensities, I_t for the
predicted intensity of t, and C for the summed predicted intensity of the
claimed components (C = 0 for a free feature).  With s the scatter of
ln(I_obs / I_pred) measured on the accepted lines,

    free feature      G_t = N(ln I_obs; ln(I_t f), s) / g(ln I_obs | lambda)
    claimed feature   G_t = N(ln I_obs; ln((C + I_t) f), s)
                          / N(ln I_obs; ln(C f), s)

with g the local distribution of ln I_obs over the observed line list, which is
what an unrelated line would have been drawn from.  Both numerator and
denominator are densities of the same variable, so no Jacobian survives.

Each N above is floored at g on ONE side - where the feature is brighter than
predicted, but not where it is dimmer.  The asymmetry is the physics: a feature
can always carry light from transitions nobody has identified, so an excess
over the prediction is not evidence against anything, while a deficit cannot be
explained away because light already emitted cannot be taken back.  Without the
floor the denominator of a blended feature is evaluated at the summed
prediction of its other components, and where those happen to be weak that
density is astronomically small - so a level collects enormous credit for
"explaining" a brightness the components it is compared with never claimed to
explain.  One such row was worth ln G = +22, and it alone moved a solidly
established level 219 cm^-1 onto the position of a different real level.

This is the blend penalty (factor 3), and it is not an added term but the same
formula with two densities changed.  Under H0 a claimed feature is present with
certainty, so its presence is no evidence for L: in the positional part rho_t
is replaced by 1/(2 W_t), the density of a line known to be somewhere in the
window, which drops the best attainable R_t from about 25 to 3.2.  And in the
intensity part the question stops being "is a line of this brightness here" and
becomes "does adding the predicted intensity of t improve the account of a
brightness already explained" - so a strong prediction dumped on a feature that
is fully explained without it is charged for the excess, exactly as it should
be.

Note that ln(I_obs (C+I_t)^-1 ... ) is the same number as the branching-fraction
form ln(I_obs BF / (I_t f)) used elsewhere in the pipeline, since BF is by
definition I_t / (C + I_t): one transition's prediction is never compared with a
whole blended feature.  The C form is used here because it is defined at any
scanned position, where BF is not.

The level offset.  All of a level's lines share whatever the calculation got
wrong about that level, so the residuals r_i = ln(I_obs BF / (I_t f)) are
written r_i = mu_L + e_i with mu_L ~ N(0, s_L).  mu_L is marginalised, not
fitted: the intensity model is absolutely calibrated (mean residual -0.02) and
the level offset explains only a ninth of the variance, so fitting it freely
would throw away the calibration to buy very little.  Marginalising n free rows
adds the single correction

    -0.5 ln(1 + n s_L^2/s^2) + s_L^2 (sum r)^2 / (2 s^2 (s^2 + n s_L^2))

to the total, with each row weighted by its posterior of being a genuine line
rather than a coincidence.  The weighting is what keeps the correction bounded:
the positive part of that expression is cancelled by the -r^2/2s^2 the same
rows carry in the base term, but a row whose line is a coincidence has its
ln R_t floored at ln(1 - p) and so carries no such term - an unweighted sum
would collect the gain without ever having paid for it, and a displaced
position with a few wildly mismatched intensities would score hundreds of units
it has not earned.


4.  What is measured, and what it is measured on
------------------------------------------------
Everything the formula needs is fitted to the run itself, by build(), and
printed at the head of the report so that no number in it is a guess.

sigma_t.  The quoted wavenumber uncertainty of the matched line, inflated by
two measured factors and combined with the partner energy uncertainty:

    sigma_t^2 = (k(n) k(char) unc_wn_obs)^2 + u_M^2

k(n) is the rms of t = (wn_obs - rwn)/unc_wn_obs among features carrying n
accepted transitions - about 1.01 unblended and 1.47 for a two-component blend,
so a blend is less well placed than its already inflated uncertainty admits.
k(char) is the same rms by line-character code, on residuals already divided by
k(n).  Both are shrunk toward 1 by a pseudo-count so that a code with five
lines does not acquire a factor of its own.  u_M is the partner level D1 from
LOPT_output_levels.txt, 0.002 - 0.03 cm^-1; adding it in quadrature slightly
overstates sigma, since the two are positively correlated, and the direction is
conservative.

eta is fitted by maximum likelihood on the accepted lines, as the mixing weight
that makes (1 - eta) N(d; 0, sigma) + eta rho fit their residuals best.

rho(nu) is the local density of the observed line list, from the K nearest
recorded lines; g(ln I | lambda) is the local mean and spread of ln I_obs over
a wider neighbourhood; s and s_L are the within- and between-level spreads of
the intensity residual over levels with at least five accepted lines.

Transitions that carry no positional information are dropped, which is the
caveat the alignment factor needs: a prediction is dropped when its partner
level M would have no accepted line left after removing the lines M shares
with L, because M's energy was then fitted to those lines and the residual is
zero by construction.  In this run the seven such lines have residual 0.000
exactly.  (A hat-matrix treatment generalising this continuously was tried and
rejected: LOPT fits a blended feature through a centroid model, so its normal
equations are not those of independent observations, and the leverages that
come out of assuming they are make the residuals worse rather than flatter.)


5.  The scan, and the question marks
------------------------------------
The immediate use.  Scan E across the interval where the calculation allows the
level to be - E_calc from IDEN2/enlev.dat, plus or minus three times the rms of
E_obs - E_calc over the known levels of the same dominant configuration, which
is level_interchange.configuration_windows() - and ask whether the adopted
position is the only maximum of ln R.  Every local maximum within `--alt-drop`
of the adopted one is an alternate position, and a level with any is a level
that has earned a question mark.  The count of them is reported as a column of
its own: it is the criterion, not any single p-value.

The predicted intensities are held fixed while E moves.  Over the few hundred
cm^-1 of a scan window the Boltzmann factor exp(-E/kT) with kT = 12900 cm^-1
changes by about two per cent, far below s.

6.  The audit: what to do about a level that has one
----------------------------------------------------
ln R says which of two energies the lines prefer.  It does not say the
preference is worth acting on, and four quite different situations produce the
same number, so --audit separates them before anything is called a revision.

Three corrections stand between "the lines prefer that energy" and "move the
level".

THE LOOK-ELSEWHERE CORRECTION.  A level whose configuration is 4f.5d.6p is
scanned over +/- 7000 cm^-1 and a well-determined one over +/- 130, so the
first gets hundreds of chances at a good maximum and the second a handful.
Their gains cannot be compared raw.  The scan already counts the chances -
n_alt, the local maxima within --alt-drop of the adopted peak - and under H0
the heights of those maxima are close to exponentially distributed, so the
best of n_alt of them stands about ln(n_alt) above a typical one.  The
corrected gain is

    look = [ln R(alternate) - ln R(adopted)] - ln(n_alt) .

THE FREE-LINE TEST.  Count only the alternate's matched rows that are
well-centred (within FREE_SIGMA) and sit on features that no accepted
transition claims (C = 0).  Those are the lines the level could take without
robbing another level.  Support that is entirely blended is not support: a
prediction can be dumped on an already-explained feature almost anywhere, and
the blend branch of the formula gives it credit for doing so.  Alongside it,
top_share - the largest single row's share of all the positive evidence -
throws out the positions whose case is one lucky line.

THE VACANCY TEST, which is the only test in this module that comes from the
theory rather than from the line list, and the only one that can rule a
position out outright.  A preferred alternate can be read two ways: the level
moves there, or an unknown level sits there and the lines are its.  The second
reading needs the calculation to still have a level of the right J and parity
spare at that energy, and n_vacant counts them - unfound rows of enlev.dat, of
the level's own J and parity, within one configuration window of the
alternate.  Where n_vacant is zero the new-level reading is dead however good
the lines look.  It is a real filter: at the alternate position of
059003.000305, eight unassigned, well-centred lines with coherent intensities
looked like an undiscovered level, and the nearest unfound calculated level of
J = 9/2 is 6900 cm^-1 away in the odd system and 18900 in the even one.  There
is nothing there to find.

What survives is sorted into five dispositions, named rather than ranked
because they are five different pieces of work: `interchange` (the alternate
lands on another level of the run - level_interchange.py's question, on
evidence this scan does not use), `refit` (the alternate is a fraction of a
wavenumber away, so the level stays and some accepted line is dragging the
LOPT fit), `relocate` (a free position, broadly supported, surviving the
look-elsewhere correction), `weak` (a free position with thin support) and
`no support`.

One caution the report repeats, because it is the easiest mistake to make with
this tool: LEVELS SHARE LINES, so every ln R is conditional on the rest of the
level list, and a revision changes the verdict on its neighbours.  Accept them
one at a time and re-scan.  059003.000604 is the demonstration - with
059003.000371 at its old position the adopted energy of 000604 wins, and only
after 000371 moves does the alternate at 136728.75 become the best position in
the whole window.

7.  Reading the report
----------------------
WHAT ln R MEANS.  ln R is a log odds, so its sign is the whole of its meaning.

    ln R  >  0   the recorded lines are more likely if a real level sits at
                 this energy than if nothing does - the position is supported,
                 and the larger the number the stronger the support.  Adding
                 10 to ln R means the data are e^10, about 20000 times, more
                 likely under the level than without it.  A well-centred line
                 on an open plate is worth about +3.2, so ln R = +30 is a
                 level carrying ten of them.
    ln R  =  0   the lines say nothing either way.
    ln R  <  0   the recorded lines are LESS likely if a level sits there:
                 the position predicts transitions that should have been seen
                 and were not, or lines whose intensities contradict the
                 prediction.  A strongly negative ln R is the signature of a
                 spurious position.

The same number is used for a candidate energy as for the adopted one, which
is what makes the two comparable: gain = ln R(alternate) - ln R(adopted) is a
log odds ratio between two positions of the same level.

WHAT COUNTS AS A PREDICTION.  n_pred is the number of transitions the
calculation gives the level anywhere in the recorded range, and on its own it
says very little: most of them are far too faint to have been recorded at all,
and their number is set by how many partner levels are known, not by the
experiment.  n_obs is the count that matters - the predictions whose
probability of having been recorded, P_obs = c(lambda) D(z) of section 2, is
at least P_SEEN.  Those are the transitions the level can actually be judged
on.  n_seen is how many of the n_obs were found, and n_miss = n_obs - n_seen
how many were not.  A level with n_obs = 9, n_seen = 8 is in good order
whatever n_pred says; n_obs = 9 with n_seen = 2 is not.  (n_match, kept in the
csv, counts every prediction that found a line including the faint ones whose
line is almost certainly a coincidence, and is the number NOT to quote.)

THE COLUMNS, in the order they are written.

    E            the adopted energy, cm^-1
    n_drop       predictions discarded because the partner level would have no
                 accepted line left without the ones it shares with this level
                 (section 4); they carry no positional information
    n_pred       predicted transitions in range
    n_obs        of those, the ones that could have been recorded
    n_seen       observable predictions that found a line
    n_miss       observable predictions that did not - the absences
    n_match      every prediction that found a line, faint coincidences and
                 all; n_seen is the meaningful count
    n_claimed    matched predictions whose feature other accepted transitions
                 already explain
    ln_R         the total, at the adopted energy
    ln_R_match   the part of it contributed by the matched predictions
    ln_R_miss    the part contributed by the absences (never positive)
    sum_lnG      the part contributed by intensity agreement alone

with --scan, for the best alternate position found:

    scan_width   the width of the interval scanned, cm^-1
    n_alt        local maxima of ln R within --alt-drop of the adopted peak.
                 This is the count of alternate positions and it, not any
                 single probability, is what earns the level its question mark
    ln_R_alt     ln R at the best of them
    d_ln_R       ln_R - ln_R_alt.  NEGATIVE means the lines prefer the
                 alternate.  (The audit's `gain` is the same quantity with the
                 clearer sign: gain = -d_ln_R, positive = the alternate wins.)
    dE_alt       E_alternate - E_adopted, cm^-1: how far the move would be
    near_level   the level of the run whose energy is CLOSEST TO THE
                 ALTERNATE.  It is a neighbourhood label, not a rival: it is
                 filled in for every alternate, and for most of them it names
                 a level that merely happens to lie nearby and wants nothing.
                 Only when near_dE is small - under INTERCHANGE_DE, half a
                 wavenumber - does it mean anything, and then it means the
                 alternate IS that level's position, so the two identities may
                 need exchanging rather than either of them moving.  A
                 near_level several cm^-1 away is noise; read dE_alt instead.
    near_dE      E_alternate - E(near_level), cm^-1
    question     '?' when n_alt > 0

with --audit, for that same alternate:

    n_free       matched lines at the alternate that no accepted transition
                 claims and that sit within FREE_SIGMA of prediction
    free_gain    what those free lines are worth in ln R
    top_share    the largest single line's share of the positive evidence
    gain         ln_R_alt - ln_R
    look         gain - ln(n_alt), the look-elsewhere correction of section 6
    n_obs_alt    observable predictions at the alternate
    n_seen_alt   how many of them would find a line there
    n_miss_alt   how many would be absences there
    n_vacant     unfound calculated levels of the same J and parity within one
                 configuration window of the alternate.  Zero kills the
                 reading "an unknown level sits there"; it does not kill "this
                 level moves there"
    z_alt        (E_alternate - E_calc)/W, the distance from where the
                 calculation puts the level in units of that configuration's
                 own scatter
    action       the disposition of section 6

THE ORDER OF THE ROWS.  The report is written best-first, so it can be read
from the top: with --audit, relocations before interchanges before refits
before the weak ones, and within each group the largest look-elsewhere
corrected gain; with --scan alone, the largest gain; with neither, the weakest
positions first, since a plain report is a list of complaints.  The csv keeps
that order.  The tables printed under "the twenty weakest positions" are
sorted by ln R regardless, because that is what they are for.

"""

import argparse
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import chance_mc as mc                        # noqa: E402
import classify_lines as cl                   # noqa: E402
import level_interchange as li                # noqa: E402
import level_shifts as ls                     # noqa: E402
import lopt_lines                             # noqa: E402


# --- constants --------------------------------------------------------------
N_SIGMA = 4.0        # half-width of the matching window, in sigma
P_SEEN = 0.2         # P_obs at which a prediction counts observable
Q_OUT = 2.0 * 0.5 * math.erfc(N_SIGMA / math.sqrt(2.0))   # 6.3e-5
K_RHO = 41           # recorded lines averaged for the local line density
K_BG = 201           # recorded lines averaged for the local intensity spread
K_UNC = 101          # recorded lines averaged for the local uncertainty
SHRINK_N = 25.0      # pseudo-count pulling k(n) toward 1
SHRINK_CHAR = 50.0   # pseudo-count pulling k(char) toward 1
K_WORST = 1.6        # the largest k(n) k(char) the run produces, used
                     # to widen the matching window (see ln_ratio)
ETA_MAX = 0.25       # upper bound of the anomalous-position rate
MIN_LEVEL_LINES = 5  # accepted lines a level needs to enter the s_L estimate
GRID_STEP = 0.02     # cm^-1, the scan step
GRID_MAX = 300000    # cap on the number of scan points per level
BLOCK = 4000         # scan points evaluated in one array operation
ALT_SEP = 0.5        # cm^-1: two maxima closer than this are one maximum
ALT_DROP = 5.0       # ln R below the adopted peak that still counts as an
                     # alternate position
LOPT_LEVELS_FILE = os.path.join(HERE, 'LOPT_output_levels.txt')

_LN_2PI = math.log(2.0 * math.pi)


# ---------------------------------------------------------------------------
# Vectorised forms of the level_shifts scalar helpers
# ---------------------------------------------------------------------------
def noise_threshold_vec(wn, calib):
    """level_shifts.noise_threshold_linear over an array; nan outside range."""
    lam = 1.0e8 / np.asarray(wn, dtype=float)
    his = np.array([hi for _, hi, _ in calib])
    out = np.full(lam.shape, np.nan)
    inside = (lam >= calib[0][0]) & (lam <= calib[-1][1])
    j = np.clip(np.searchsorted(his, lam, side='left'), 0, len(calib) - 1)
    for r in np.unique(j[inside]) if inside.any() else ():
        m = inside & (j == r)
        p = np.zeros(int(m.sum()))
        x = lam[m]
        for k, ck in enumerate(calib[r][2]):
            p += ck * x ** k
        out[m] = np.exp(p) * 1000.0
    return out


def coverage_vec(wn):
    """level_shifts.coverage over an array."""
    lam = 1.0e8 / np.asarray(wn, dtype=float)
    m = ls.read_coverage_map()
    if not m:
        c = np.ones_like(lam)
        for lo, hi in ls.COVERAGE_GAPS_A:
            c[(lam >= lo) & (lam <= hi)] = 0.0
        return c
    return np.interp(lam, m[0], m[1])


def scale_factor_vec(wn, bias):
    """level_shifts.scale_bias_factor over an array."""
    lam = 1.0e8 / np.asarray(wn, dtype=float)
    base = bias.get(None, (1.0, 1.0))[0] if bias else 1.0
    f = np.full(lam.shape, float(base))
    for band, v in (bias or {}).items():
        if band is None:
            continue
        f[(lam >= band[0]) & (lam < band[1])] = v[0]
    return f


def p_obs_vec(wn, i_pred, calib, bias):
    """level_shifts.observation_probability over an array: c(lambda) D(z)."""
    wn = np.asarray(wn, dtype=float)
    i_pred = np.asarray(i_pred, dtype=float)
    c = coverage_vec(wn)
    if calib is None:
        return c
    thr = noise_threshold_vec(wn, calib)
    f = scale_factor_vec(wn, bias)
    ok = np.isfinite(thr) & (thr > 0) & (i_pred > 0)
    z = np.zeros(wn.shape)
    z[ok] = np.log(i_pred[ok] * f[ok] / thr[ok])
    cur = ls.read_detection_curve()
    if cur:
        d = np.interp(z, cur[0], cur[1])
    else:
        d = 0.5 * (1.0 + np.vectorize(math.erf)(z / math.sqrt(2.0)))
    d = np.where(ok, d, 1.0)
    return c * d


def ln_intensity(x, predicted, s, ln_bg):
    """ln of the density of x = ln I_obs given a predicted total `predicted`.

    N(x; ln predicted, s) where the feature is dimmer than predicted, and
    max of that and the background density where it is brighter.  The
    asymmetry is the physics: a feature can always carry light from
    transitions nobody has identified, so an excess over the prediction is
    not evidence against anything and must not be charged as though it were;
    a deficit cannot be explained away, because light already emitted cannot
    be taken back.

    Without the floor the denominator of a blended feature is evaluated at
    the summed prediction of its other components, and where those are weak
    the density is astronomically small - so a level acquires enormous credit
    for "explaining" a brightness the components it is being compared with
    never claimed to explain.  That defect alone moved a solidly established
    level 219 cm^-1, onto the position of a different real level.
    """
    m = np.log(np.maximum(np.asarray(predicted, dtype=float), 1e-300))
    ln_n = ln_norm(x, m, s)
    return np.where((x > m) & (ln_bg > ln_n), ln_bg, ln_n)


def ln_norm(x, mu, sd):
    """ln of the normal density, elementwise."""
    z = (np.asarray(x, dtype=float) - mu) / sd
    return -0.5 * (z * z) - np.log(sd) - 0.5 * _LN_2PI


# ---------------------------------------------------------------------------
# The run, and everything measured on it
# ---------------------------------------------------------------------------
class Context:
    """The run under test with every ingredient of the likelihood attached."""


def running(values, k, fn):
    """fn over a centred window of k neighbours, as an array of the same
    length; used for the local density, uncertainty and intensity spread."""
    s = pd.Series(values).rolling(k, center=True, min_periods=max(3, k // 8))
    return getattr(s, fn)().ffill().bfill().to_numpy()


def read_partner_uncertainties(path=LOPT_LEVELS_FILE):
    """{level_id: D1}, the partner-energy uncertainty of LOPT's level output.

    D1 is what a Ritz wavenumber inherits from the partner level, 0.002 to
    0.03 cm^-1 here.  Returns {} if the file is absent, in which case the
    partner contributes nothing to sigma_t and the report says so.
    """
    if not os.path.exists(path):
        return {}
    d = pd.read_csv(path, sep='\t', dtype={'Designation': str})
    return {str(k): float(v) for k, v in zip(d['Designation'], d['D1'])
            if pd.notna(v)}


def fit_k_blend(t, n, log):
    """k(n): how much the quoted uncertainty understates a blend of n parts."""
    out, rows = {}, []
    for nb in sorted(set(int(x) for x in n if x >= 1)):
        m = (n == nb)
        cnt = int(m.sum())
        if cnt == 0:
            continue
        ms = float(np.mean(t[m] ** 2))
        k = math.sqrt((cnt * ms + SHRINK_N) / (cnt + SHRINK_N))
        out[nb] = k
        rows.append((nb, cnt, math.sqrt(ms), k))
    log("  blend inflation of the position uncertainty, k(n):")
    for nb, cnt, rms, k in rows:
        log(f"    n={nb:2d}  N={cnt:5d}  rms(t)={rms:.3f}  k={k:.3f}")
    return out


def fit_k_char(t, char, log):
    """k(char): the same by line-character code, on t already divided by k(n)."""
    out, rows = {}, []
    for c in sorted(set(char)):
        m = (char == c)
        cnt = int(m.sum())
        ms = float(np.mean(t[m] ** 2))
        k = math.sqrt((cnt * ms + SHRINK_CHAR) / (cnt + SHRINK_CHAR))
        out[c] = k
        rows.append((c, cnt, math.sqrt(ms), k))
    rows.sort(key=lambda r: -r[1])
    log("  character-code inflation, k(char) (codes with 20 lines or more):")
    for c, cnt, rms, k in rows:
        if cnt >= 20:
            log(f"    {c if c else '(none)':<8s} N={cnt:5d}  rms={rms:.3f}  "
                f"k={k:.3f}")
    return out


def fit_eta(d, sigma, rho, log):
    """The anomalous-position rate, by maximum likelihood on the accepted
    lines: the weight eta that best mixes the background density rho into the
    Gaussian.  A grid on ln eta, refined; the surface is unimodal."""
    good = np.isfinite(d) & np.isfinite(sigma) & (sigma > 0) & (rho > 0)
    d, sigma, rho = d[good], sigma[good], rho[good]
    g = np.exp(-0.5 * (d / sigma) ** 2) / (sigma * math.sqrt(2.0 * math.pi))

    def nll(e):
        v = (1.0 - e) * g + e * rho
        return -float(np.sum(np.log(np.maximum(v, 1e-300))))

    lo, hi = 1e-5, ETA_MAX
    for _ in range(60):
        a = lo + (hi - lo) / 3.0
        b = hi - (hi - lo) / 3.0
        if nll(a) < nll(b):
            hi = b
        else:
            lo = a
    eta = 0.5 * (lo + hi)
    log(f"  anomalous-position rate eta = {eta:.4f} "
        f"(2 dNLL against eta=0: {2 * (nll(0.0) - nll(eta)):.1f})")
    return eta


def fit_intensity_scatter(r, level_of, log):
    """(s, s_L): the within-level and between-level spread of the intensity
    residual r = ln(I_obs BF / (I_pred f))."""
    ok = np.isfinite(r)
    r, level_of = r[ok], np.asarray(level_of)[ok]
    df = pd.DataFrame({'r': r, 'lid': level_of})
    g = df.groupby('lid')['r']
    big = g.count()[g.count() >= MIN_LEVEL_LINES].index
    sub = df[df['lid'].isin(big)]
    if len(sub) > 10:
        within = sub.groupby('lid')['r'].transform('mean')
        s = float(np.std(sub['r'] - within, ddof=1))
        means = sub.groupby('lid')['r'].mean()
        sd_means = float(np.std(means, ddof=1))
        n_bar = float(sub.groupby('lid')['r'].count().mean())
        s_l = math.sqrt(max(sd_means ** 2 - s ** 2 / n_bar, 1e-4))
    else:
        s, s_l = float(np.std(r, ddof=1)), 0.1
    log(f"  intensity residual: mean {float(np.mean(r)):+.3f}, sd "
        f"{float(np.std(r, ddof=1)):.3f} over {len(r)} accepted lines")
    log(f"  within-level s = {s:.3f}, between-level s_L = {s_l:.3f} "
        f"({len(big)} levels with {MIN_LEVEL_LINES} lines or more)")
    return s, s_l


def complete_accepted(acc, ctx, log=print):
    """Fill in the columns a LOPT run does not carry, in place.

    lopt_lines.read_run gives the accepted set but leaves rwn (the Ritz
    wavenumber), calc_intens (the predicted intensity of the component) and
    therefore BF empty, because LOPT knows nothing of the theoretical
    intensities.  Both are recoverable and both are needed here - rwn for the
    position residuals that measure k(n), k(char) and eta, calc_intens for the
    claimed intensity C of a blended feature.  The Ritz wavenumber is the
    difference of the two optimized energies; the predicted intensity is the
    Icalc entry for that pair, which is the same number classify_lines.py
    would have written; BF follows from the shares within each feature.
    """
    n = len(acc)
    if acc['rwn'].isna().any():
        miss = acc['rwn'].isna()
        acc.loc[miss, 'rwn'] = (acc.loc[miss, 'upp_E']
                                - acc.loc[miss, 'low_E'])
        log(f"  rwn reconstructed from the optimized energies for "
            f"{int(miss.sum())} of {n} accepted rows")
    if acc['calc_intens'].isna().any():
        ip = {(str(a), str(b)): float(c)
              for a, b, c in ctx.preds_all.itertuples(index=False)}
        miss = acc['calc_intens'].isna()
        acc.loc[miss, 'calc_intens'] = [
            ip.get((a, b), np.nan)
            for a, b in zip(acc.loc[miss, 'low_id'], acc.loc[miss, 'upp_id'])]
        got = int(acc['calc_intens'].notna().sum())
        log(f"  calc_intens taken from Icalc for {int(miss.sum())} of {n} "
            f"accepted rows ({n - got} still without a prediction)")
        tot = acc.groupby('wn_obs')['calc_intens'].transform('sum')
        bf = acc['calc_intens'] / tot.where(tot > 0)
        acc['BF'] = bf.fillna(1.0 / acc['n_accepted'].clip(lower=1))
    return acc


def build(args, log=print):
    """Read the run and measure every ingredient of the likelihood."""
    ctx = Context()
    if args.lopt:
        ls.LOPT_LINES = args.lopt
        ls.LOPT_LEVELS = args.lopt_levels
        ls.ENERGIES = args.energies
    if args.energies_csv:
        ls.E_INPUT_CSV = args.energies_csv

    levels, real, _ = ls.read_run()
    per = mc.per_level_table(real, levels)
    ctx.per = per
    ctx.e_final = dict(zip(per['level_id'], per['E_final']))
    ctx.real = real
    log(f"levels in the run: {len(per)}")

    # --- the complete observed line list --------------------------------
    obs = lopt_lines.read_observed_lines()
    obs = obs.dropna(subset=['wn_obs']).sort_values('wn_obs')
    obs = obs.drop_duplicates('wn_obs').reset_index(drop=True)
    obs['char'] = obs['char'].fillna('').astype(str)
    ctx.obs = obs
    ctx.wn_o = obs['wn_obs'].to_numpy()
    ctx.int_o = obs['obs_intens'].to_numpy(dtype=float)
    ctx.unc_o = obs['unc_wn_obs'].to_numpy(dtype=float)
    ctx.char_o = obs['char'].to_numpy()
    with np.errstate(divide='ignore', invalid='ignore'):
        ctx.ln_int_o = np.log(np.where(ctx.int_o > 0, ctx.int_o, np.nan))
    log(f"observed lines: {len(obs)}, {ctx.wn_o[0]:.1f} - {ctx.wn_o[-1]:.1f} "
        f"cm^-1")

    # local density, local uncertainty, local intensity distribution
    ctx.rho_o = K_RHO / (
        ctx.wn_o[np.clip(np.arange(len(ctx.wn_o)) + K_RHO // 2, 0,
                         len(ctx.wn_o) - 1)]
        - ctx.wn_o[np.clip(np.arange(len(ctx.wn_o)) - K_RHO // 2, 0,
                           len(ctx.wn_o) - 1)])
    ctx.sig_loc_o = running(ctx.unc_o, K_UNC, 'median')
    ctx.bg_mu_o = running(ctx.ln_int_o, K_BG, 'mean')
    ctx.bg_sd_o = running(ctx.ln_int_o, K_BG, 'std')
    ctx.bg_sd_o = np.maximum(ctx.bg_sd_o, 0.2)
    log(f"  local line density {np.min(ctx.rho_o):.3f} - "
        f"{np.max(ctx.rho_o):.3f} per cm^-1 (median "
        f"{np.median(ctx.rho_o):.3f})")
    log(f"  local ln I_obs: mean {np.nanmedian(ctx.bg_mu_o):.2f}, sd "
        f"{np.nanmedian(ctx.bg_sd_o):.2f}")

    # --- the calibration, the coverage and the detection curve ----------
    ctx.calib = ls.read_intensity_calibration()
    if ctx.calib is None:
        raise SystemExit("intensity_correction_functions.txt not found")
    ctx.preds_all = ls.load_predictions(set(per['level_id']), ctx.e_final)
    ctx.bias = ls.intensity_scale_bias(real, ctx.preds_all, ctx.e_final)
    log(f"predicted transitions: {len(ctx.preds_all)}")

    # --- the accepted set: claimed features, partners, ingredients ------
    acc = real[real['accepted'] == 1].copy()
    acc['low_id'] = acc['low_id'].astype(str)
    acc['upp_id'] = acc['upp_id'].astype(str)
    complete_accepted(acc, ctx, log)
    ctx.acc = acc
    # every accepted row placed on its observed line
    idx = np.searchsorted(ctx.wn_o, acc['wn_obs'].to_numpy())
    idx = np.clip(idx, 0, len(ctx.wn_o) - 1)
    left = np.clip(idx - 1, 0, len(ctx.wn_o) - 1)
    take_left = (np.abs(ctx.wn_o[left] - acc['wn_obs'].to_numpy())
                 < np.abs(ctx.wn_o[idx] - acc['wn_obs'].to_numpy()))
    line_idx = np.where(take_left, left, idx)
    acc['line_idx'] = line_idx
    miss = np.abs(ctx.wn_o[line_idx] - acc['wn_obs'].to_numpy()) > 1e-6
    if miss.any():
        log(f"  WARNING: {int(miss.sum())} accepted rows have no exact line in "
            f"the observed workbook")

    # C, the predicted intensity claimed on each observed line, and who claims
    ctx.claimed_tot = np.zeros(len(ctx.wn_o))
    np.add.at(ctx.claimed_tot, line_idx,
              acc['calc_intens'].fillna(0.0).to_numpy(dtype=float))
    ctx.n_acc_line = np.zeros(len(ctx.wn_o), dtype=int)
    np.add.at(ctx.n_acc_line, line_idx, 1)
    own = {}
    for lid, li_, ic in zip(acc['low_id'], acc['line_idx'],
                            acc['calc_intens'].fillna(0.0)):
        own.setdefault(lid, {})
        own[lid][int(li_)] = own[lid].get(int(li_), 0.0) + float(ic)
    for lid, li_, ic in zip(acc['upp_id'], acc['line_idx'],
                            acc['calc_intens'].fillna(0.0)):
        own.setdefault(lid, {})
        own[lid][int(li_)] = own[lid].get(int(li_), 0.0) + float(ic)
    ctx.own_claim = own          # {level_id: {line index: its own I_pred}}

    # how many accepted lines each level has, and how many it shares with each
    # partner - the degeneracy rule of section 4
    ctx.n_acc_level = {}
    for lid, m in own.items():
        ctx.n_acc_level[lid] = len(m)
    shared = {}
    for lo, up in zip(acc['low_id'], acc['upp_id']):
        shared[(lo, up)] = shared.get((lo, up), 0) + 1
    ctx.shared = shared

    # --- the measured ingredients ---------------------------------------
    log("ingredients measured on this run:")
    a = acc.dropna(subset=['rwn', 'unc_wn_obs']).copy()
    a = a[a['unc_wn_obs'] > 0]
    t = ((a['wn_obs'] - a['rwn']) / a['unc_wn_obs']).to_numpy()
    nb = a['n_accepted'].fillna(1).to_numpy(dtype=float)
    ctx.k_n = fit_k_blend(t, nb, log)
    kn = np.array([ctx.k_n.get(int(x), 1.0) for x in nb])
    ctx.k_char = fit_k_char(t / kn, a['char'].fillna('').astype(str).to_numpy(),
                            log)
    kc = np.array([ctx.k_char.get(c, 1.0)
                   for c in a['char'].fillna('').astype(str)])
    sigma = kn * kc * a['unc_wn_obs'].to_numpy()
    rho_a = np.interp(a['wn_obs'].to_numpy(), ctx.wn_o, ctx.rho_o)
    ctx.eta = fit_eta((a['wn_obs'] - a['rwn']).to_numpy(), sigma, rho_a, log)

    f = np.array([ls.scale_bias_factor(w, ctx.bias) for w in a['wn_obs']])
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.log(a['obs_intens'].to_numpy(dtype=float)
                   * a['BF'].fillna(1.0).to_numpy(dtype=float)
                   / (a['calc_intens'].to_numpy(dtype=float) * f))
    lids = np.where(a['low_id'].astype(str) != '', a['low_id'].astype(str),
                    a['upp_id'].astype(str))
    ctx.s, ctx.s_L = fit_intensity_scatter(r, lids, log)

    ctx.u_M = read_partner_uncertainties()
    if ctx.u_M:
        v = np.array(list(ctx.u_M.values()))
        log(f"  partner energy uncertainty D1: median {np.median(v):.4f}, "
            f"max {np.max(v):.4f} cm^-1 over {len(v)} levels")
    else:
        log("  LOPT_output_levels.txt not found - the partner energy "
            "contributes nothing to sigma_t")

    # --- predictions grouped by level ------------------------------------
    ctx.by_level = group_predictions(ctx)
    return ctx


def group_predictions(ctx):
    """{level_id: frame of its predicted transitions}.

    sign is +1 where the level is the upper partner, so that nu = sign*(E-E_M)
    with E_M the partner energy; drop marks the predictions whose partner would
    have no accepted line left without the ones it shares with this level.
    """
    p = ctx.preds_all
    out = {}
    lo = p['lo_id'].astype(str).to_numpy()
    up = p['up_id'].astype(str).to_numpy()
    ip = p['I_pred'].to_numpy(dtype=float)
    rows = {}
    for i in range(len(p)):
        rows.setdefault(up[i], []).append((lo[i], +1.0, ip[i]))
        rows.setdefault(lo[i], []).append((up[i], -1.0, ip[i]))
    for lid, rr in rows.items():
        partner = np.array([x[0] for x in rr])
        sign = np.array([x[1] for x in rr])
        i_pred = np.array([x[2] for x in rr])
        e_m = np.array([ctx.e_final.get(x, np.nan) for x in partner])
        u_m = np.array([ctx.u_M.get(x, 0.0) for x in partner])
        n_p = np.array([ctx.n_acc_level.get(x, 0) for x in partner])
        sh = np.array([ctx.shared.get((x, lid), 0) + ctx.shared.get((lid, x), 0)
                       for x in partner])
        out[lid] = dict(partner=partner, sign=sign, i_pred=i_pred, e_m=e_m,
                        u_m=u_m, degenerate=(n_p - sh) <= 0)
    return out


# ---------------------------------------------------------------------------
# The likelihood ratio
# ---------------------------------------------------------------------------
def ln_ratio(ctx, level_id, e_grid, detail=False):
    """ln R(level, E) for every E in e_grid.

    Returns the array of ln R, and - when detail is asked for a single E - a
    frame with one row per predicted transition showing what it contributed.
    """
    g = ctx.by_level.get(level_id)
    e_grid = np.atleast_1d(np.asarray(e_grid, dtype=float))
    if g is None:
        return np.zeros(len(e_grid)), None
    keep = ~g['degenerate'] & np.isfinite(g['e_m']) & (g['i_pred'] > 0)
    if not keep.any():
        return np.zeros(len(e_grid)), None
    sign, e_m = g['sign'][keep], g['e_m'][keep]
    i_pred, u_m = g['i_pred'][keep], g['u_m'][keep]
    partner = g['partner'][keep]
    own = ctx.own_claim.get(level_id, {})

    total = np.zeros(len(e_grid))
    rows = None
    for a in range(0, len(e_grid), BLOCK):
        e = e_grid[a:a + BLOCK]
        # nu[t, k]: the Ritz wavenumber of transition t at candidate energy k
        nu = sign[:, None] * (e[None, :] - e_m[:, None])
        ok = (nu > cl.WN_MIN) & (nu < cl.WN_MAX)
        nu_s = np.where(ok, nu, cl.WN_MIN + 1.0)

        # what would have been seen there, and where a line could hide
        p_t = p_obs_vec(nu_s.ravel(),
                        np.repeat(i_pred, nu_s.shape[1]),
                        ctx.calib, ctx.bias).reshape(nu_s.shape)
        sig_loc = np.interp(nu_s, ctx.wn_o, ctx.sig_loc_o)

        # the nearest recorded line
        j = np.searchsorted(ctx.wn_o, nu_s)
        j1 = np.clip(j, 0, len(ctx.wn_o) - 1)
        j0 = np.clip(j - 1, 0, len(ctx.wn_o) - 1)
        d1 = ctx.wn_o[j1] - nu_s
        d0 = ctx.wn_o[j0] - nu_s
        pick = np.where(np.abs(d0) < np.abs(d1), j0, j1)
        d = ctx.wn_o[pick] - nu_s

        # The window is set by whichever is worse, the uncertainty typical of
        # the neighbourhood or the one the candidate line itself carries: a
        # line quoted to 0.6 cm^-1 where its neighbours are quoted to 0.1 must
        # not be ruled out for lying 0.4 away.  K_WORST covers the largest
        # k(n) k(char) the run produces.  Widening the window costs nothing:
        # a line far out has N(d; 0, sigma)/rho far below 1 and contributes
        # ln(1 - p), the same as an absence.
        unc_near = ctx.unc_o[pick]
        unc_near = np.where(np.isfinite(unc_near) & (unc_near > 0), unc_near,
                            sig_loc)
        sig_w = np.sqrt(np.maximum(sig_loc, K_WORST * unc_near) ** 2
                        + (u_m ** 2)[:, None])
        w = N_SIGMA * sig_w
        matched = (np.abs(d) <= w) & ok

        p = p_t * (1.0 - ctx.eta)
        r_t = np.where(ok, 1.0 - p * (1.0 - Q_OUT), 1.0)

        if matched.any():
            mi = pick[matched]
            unc = ctx.unc_o[mi]
            unc = np.where(np.isfinite(unc) & (unc > 0), unc,
                           sig_loc[matched])
            # C: what other levels already claim on that feature
            own_here = np.array([own.get(int(x), 0.0) for x in mi])
            claimed = np.maximum(ctx.claimed_tot[mi] - own_here, 0.0)
            n_blend = ctx.n_acc_line[mi] - (own_here > 0).astype(int) + 1
            kn = np.array([ctx.k_n.get(int(x), ctx.k_n.get(1, 1.0))
                           for x in n_blend])
            kc = np.array([ctx.k_char.get(c, 1.0) for c in ctx.char_o[mi]])
            u_here = np.repeat(u_m, nu_s.shape[1]).reshape(nu_s.shape)[matched]
            sig = np.sqrt((kn * kc * unc) ** 2 + u_here ** 2)

            # positional density under H0
            free = claimed <= 0
            rho = np.where(free, ctx.rho_o[mi], 1.0 / (2.0 * w[matched]))
            dens = (np.exp(-0.5 * (d[matched] / sig) ** 2)
                    / (sig * math.sqrt(2.0 * math.pi)))

            # the intensity ratio
            i_here = np.repeat(i_pred, nu_s.shape[1]).reshape(nu_s.shape)[matched]
            f = scale_factor_vec(nu_s[matched], ctx.bias)
            x = ctx.ln_int_o[mi]
            ln_bg = ln_norm(x, ctx.bg_mu_o[mi], ctx.bg_sd_o[mi])
            with np.errstate(divide='ignore', invalid='ignore'):
                ln_p1 = ln_intensity(x, (claimed + i_here) * f, ctx.s, ln_bg)
                ln_p0 = np.where(free, ln_bg,
                                 ln_intensity(x, claimed * f, ctx.s, ln_bg))
            ln_g = ln_p1 - ln_p0
            ln_g = np.where(np.isfinite(ln_g), ln_g, 0.0)
            gain = np.exp(np.clip(ln_g, -50.0, 50.0))

            pm = p[matched]
            genuine = pm * dens / rho * gain
            r_t[matched] = (1.0 - pm) + genuine
            # posterior that the matched line really is the transition, not a
            # coincidence: the weight this row carries in the level offset
            w_gen = np.clip(genuine / np.maximum(r_t[matched], 1e-300),
                            0.0, 1.0)

        ln_r = np.log(np.maximum(r_t, 1e-300))
        block_total = ln_r.sum(axis=0)

        # The marginalised level offset, over the free matched rows, each
        # weighted by its posterior of being genuine.  The weight is what
        # keeps the correction bounded: written for a plain Gaussian product
        # its positive part is cancelled by the -r^2/2s^2 the same rows
        # already carry, but a row whose line is a coincidence carries no such
        # term - its ln R_t is floored at ln(1 - p) - so an unweighted sum
        # would collect the gain without ever having paid for it, and a
        # displaced position with a few wildly mismatched intensities would
        # score hundreds of units it has not earned.
        if matched.any():
            weight = np.zeros(matched.shape)
            weight[matched] = w_gen * free
            with np.errstate(divide='ignore', invalid='ignore'):
                res = np.zeros(matched.shape)
                res[matched] = x - np.log(np.maximum((claimed + i_here) * f,
                                                     1e-300))
            res = np.where(np.isfinite(res), res, 0.0)
            n_w = weight.sum(axis=0)
            sum_r = (weight * res).sum(axis=0)
            s2, sl2 = ctx.s ** 2, ctx.s_L ** 2
            block_total += (-0.5 * np.log(1.0 + n_w * sl2 / s2)
                            + sl2 * sum_r ** 2 / (2.0 * s2 * (s2 + n_w * sl2)))

        total[a:a + BLOCK] = block_total

        if detail and len(e_grid) == 1:
            tab = pd.DataFrame({
                'partner': partner,
                'nu': nu_s[:, 0],
                'I_pred': i_pred,
                'P_obs': p_t[:, 0],
                'W': w[:, 0],
                'matched': matched[:, 0],
                'ln_R': ln_r[:, 0],
            })
            tab['wn_obs'] = np.where(matched[:, 0], ctx.wn_o[pick[:, 0]],
                                     np.nan)
            tab['d'] = np.where(matched[:, 0], d[:, 0], np.nan)
            tab['I_obs'] = np.where(matched[:, 0], ctx.int_o[pick[:, 0]],
                                    np.nan)
            claim_col = np.full(len(partner), np.nan)
            gain_col = np.full(len(partner), np.nan)
            if matched.any():
                claim_col[matched[:, 0]] = claimed
                gain_col[matched[:, 0]] = ln_g
            tab['C'] = claim_col
            tab['ln_G'] = gain_col
            tab = tab[ok[:, 0]]
            rows = tab.sort_values('ln_R', ascending=False).reset_index(
                drop=True)
    return total, rows


# ---------------------------------------------------------------------------
# The scan
# ---------------------------------------------------------------------------
def local_maxima(y):
    """Indices of the strict local maxima of y, ends included."""
    m = np.zeros(len(y), dtype=bool)
    m[1:-1] = (y[1:-1] >= y[:-2]) & (y[1:-1] > y[2:])
    if len(y) > 1:
        m[0] = y[0] > y[1]
        m[-1] = y[-1] > y[-2]
    return np.flatnonzero(m)


def scan_level(ctx, level_id, e_adopted, e_calc, window, step=GRID_STEP,
               alt_drop=ALT_DROP):
    """Scan ln R across the interval the calculation allows.

    Returns a dict with the adopted value, the alternate maxima ranked by
    ln R, and the grid itself.
    """
    if not np.isfinite(e_calc) or not np.isfinite(window) or window <= 0:
        lo, hi = e_adopted - 50.0, e_adopted + 50.0
    else:
        lo, hi = e_calc - 3.0 * window, e_calc + 3.0 * window
    lo = min(lo, e_adopted - 5.0)
    hi = max(hi, e_adopted + 5.0)
    n = int((hi - lo) / step) + 1
    if n > GRID_MAX:
        step = (hi - lo) / (GRID_MAX - 1)
        n = GRID_MAX
    grid = lo + step * np.arange(n)
    y, _ = ln_ratio(ctx, level_id, grid)
    y_adopted, _ = ln_ratio(ctx, level_id, np.array([e_adopted]))
    y_adopted = float(y_adopted[0])

    peaks = local_maxima(y)
    alts = []
    for i in peaks:
        if abs(grid[i] - e_adopted) < ALT_SEP:
            continue
        if y[i] >= y_adopted - alt_drop:
            alts.append((float(grid[i]), float(y[i])))
    alts.sort(key=lambda z: -z[1])
    # keep only maxima ALT_SEP apart from each other
    kept = []
    for e, v in alts:
        if all(abs(e - k[0]) >= ALT_SEP for k in kept):
            kept.append((e, v))
    i_best = int(np.argmax(y))
    return dict(level_id=level_id, e_adopted=e_adopted, ln_R=y_adopted,
                e_best=float(grid[i_best]), ln_R_best=float(y[i_best]),
                lo=lo, hi=hi, step=step, grid=grid, y=y, alternates=kept)


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------
def report_table(ctx, level_ids):
    """ln R of every level at its adopted position, with its parts."""
    rows = []
    for lid in level_ids:
        e = ctx.e_final[lid]
        v, tab = ln_ratio(ctx, lid, np.array([e]), detail=True)
        if tab is None:
            rows.append(dict(level_id=lid, E=e, n_drop=0, n_pred=0,
                             n_obs=0, n_seen=0, n_miss=0, n_match=0,
                             n_claimed=0, ln_R=0.0, ln_R_match=0.0,
                             ln_R_miss=0.0, sum_lnG=0.0))
            continue
        m = tab['matched'].to_numpy(dtype=bool)
        seen = tab['P_obs'].to_numpy(dtype=float) >= P_SEEN
        rows.append(dict(level_id=lid, E=e,
                         n_drop=int(ctx.by_level[lid]['degenerate'].sum()),
                         n_pred=len(tab),
                         n_obs=int(seen.sum()),
                         n_seen=int((seen & m).sum()),
                         n_miss=int((seen & ~m).sum()),
                         n_match=int(m.sum()),
                         n_claimed=int((tab['C'].fillna(0) > 0).sum()),
                         ln_R=float(v[0]),
                         ln_R_match=float(tab['ln_R'][m].sum()),
                         ln_R_miss=float(tab['ln_R'][~m].sum()),
                         sum_lnG=float(tab['ln_G'][m].sum())))
    return pd.DataFrame(rows)


def print_detail(ctx, level_id):
    e = ctx.e_final[level_id]
    v, tab = ln_ratio(ctx, level_id, np.array([e]), detail=True)
    print(f"\nlevel {level_id} at E = {e:.4f} cm^-1")
    print(f"ln R = {float(v[0]):+.2f} over {0 if tab is None else len(tab)} "
          f"predicted transitions")
    if tab is None:
        return
    g = ctx.by_level[level_id]
    print(f"  {int(g['degenerate'].sum())} predictions dropped: the partner "
          f"has no accepted line of its own")
    m = tab['matched'].to_numpy(dtype=bool)
    print(f"  {int(m.sum())} matched, contributing {tab['ln_R'][m].sum():+.2f}; "
          f"{int((~m).sum())} absent, costing {tab['ln_R'][~m].sum():+.2f}")
    print()
    hdr = (f"{'partner':<14}{'nu':>12}{'I_pred':>10}{'P_obs':>7}"
           f"{'wn_obs':>12}{'d':>8}{'I_obs':>10}{'C':>10}{'lnG':>7}{'lnR':>8}")
    print(hdr)
    print('-' * len(hdr))
    for _, r in tab.iterrows():
        if r['matched']:
            print(f"{r['partner']:<14}{r['nu']:>12.3f}{r['I_pred']:>10.1f}"
                  f"{r['P_obs']:>7.3f}{r['wn_obs']:>12.3f}{r['d']:>8.3f}"
                  f"{r['I_obs']:>10.1f}{r['C']:>10.1f}{r['ln_G']:>7.2f}"
                  f"{r['ln_R']:>8.2f}")
        else:
            print(f"{r['partner']:<14}{r['nu']:>12.3f}{r['I_pred']:>10.1f}"
                  f"{r['P_obs']:>7.3f}{'-':>12}{'-':>8}{'-':>10}{'-':>10}"
                  f"{'-':>7}{r['ln_R']:>8.2f}")


# ---------------------------------------------------------------------------
# The audit: what a preferred alternate position actually rests on
# ---------------------------------------------------------------------------
# ln R says which of two energies the lines prefer.  It does not say whether
# the preference is worth acting on, and four different things can produce the
# same number.  The audit separates them.
FREE_SIGMA = 2.0        # a supporting line must sit within this many sigma
AUDIT_MIN_FREE = 4      # free supporting lines a firm relocation needs
AUDIT_MIN_GAIN = 8.0    # ln R those free lines have to be worth
AUDIT_MAX_SHARE = 0.45  # largest share of the evidence one line may carry
AUDIT_ANY_SHARE = 0.60  # above this one line IS the case, and there is none
AUDIT_MIN_FEW = 2       # fewer free lines than this is not support at all
AUDIT_MIN_LOOK = 3.0    # gain a relocation must keep after look-elsewhere
REFIT_DE = 2.0          # an alternate this close is the same level, refitted
INTERCHANGE_DE = 0.5    # an alternate this close to another level is a swap


def j_value(j):
    """4.5 from the '9/2' of the level table; nan from anything unreadable."""
    try:
        s = str(j).strip()
        if '/' in s:
            a, b = s.split('/')
            return float(a) / float(b)
        return float(s)
    except (ValueError, TypeError):
        return float('nan')


def support(ctx, level_id, e):
    """What the position E rests on: its free lines, their weight, its spread.

    A matched prediction is FREE support when the observed feature is claimed
    by no other accepted transition (C = 0) and the line sits within
    FREE_SIGMA of where it is predicted.  Those are the lines the level can
    take without robbing another level of its evidence; support that is
    entirely blended is not support, because a prediction can be dumped on an
    already-explained feature almost anywhere.

    top_share is the largest single row's share of all the positive evidence.
    A position whose case is one line is not a case: one line can be a
    coincidence, and the scan looked at tens of thousands of positions.
    """
    v, tab = ln_ratio(ctx, level_id, np.array([float(e)]), detail=True)
    if tab is None:
        return dict(ln_R=0.0, n_free=0, free_gain=0.0,
                    top_share=np.nan, n_obs_alt=0, n_seen_alt=0,
                    n_miss_alt=0)
    m = tab['matched'].to_numpy(dtype=bool)
    d = np.nan_to_num(tab['d'].to_numpy(dtype=float), nan=np.inf)
    c = np.nan_to_num(tab['C'].to_numpy(dtype=float), nan=-1.0)
    w = tab['W'].to_numpy(dtype=float)
    r = tab['ln_R'].to_numpy(dtype=float)
    free = m & (c == 0.0) & (np.abs(d) <= FREE_SIGMA * w / N_SIGMA)
    pos = r[m & (r > 0.0)]
    seen = tab['P_obs'].to_numpy(dtype=float) >= P_SEEN
    return dict(ln_R=float(v[0]), n_free=int(free.sum()),
                free_gain=float(r[free].sum()),
                n_obs_alt=int(seen.sum()),
                n_seen_alt=int((seen & m).sum()),
                n_miss_alt=int((seen & ~m).sum()),
                top_share=(float(pos.max() / pos.sum()) if pos.size
                           else np.nan))


def vacancies(en):
    """The calculated levels of enlev.dat that have not been found.

    A preferred alternate can be read two ways: the level moves there, or an
    unknown level sits there and the coincidence is with its lines.  The
    second reading needs the calculation to still have a level of the right J
    and parity spare at that energy.  Where it has none the reading is dead
    however good the lines look - which is the one test in this whole module
    that comes from the theory rather than from the line list.
    """
    v = en[~en['known']].copy()
    v['par'] = v['cfg'].map(li.configuration_parity)
    return v


def count_vacancies(vac, j, parity, e, window):
    """Unfound calculated levels of this J and parity within +/- window of E.

    -1 when the level has no J, parity or window to test it against.
    """
    j = j_value(j)
    if not (np.isfinite(j) and np.isfinite(e) and np.isfinite(window)
            and window > 0 and parity in ('e', 'o')):
        return -1
    sel = ((vac['J'].to_numpy() == j) & (vac['par'].to_numpy() == parity)
           & (np.abs(vac['E_calc'].to_numpy() - e) <= window))
    return int(sel.sum())


def disposition(row):
    """What kind of problem a level with a preferred alternate actually is.

    Five different pieces of work, which is why they are named rather than
    ranked:

      interchange  the alternate lands on another level of the run.  Nothing
                   moves anywhere new; the two identities may be swapped, and
                   level_interchange.py decides that on evidence this scan
                   does not look at.
      refit        the alternate is a fraction of a wavenumber away.  The
                   level stays where it is; some accepted line is dragging
                   the LOPT fit off the position its own lines want.
      relocate     a free position, broadly supported by lines nobody is
                   using, surviving the look-elsewhere correction, and
                   convincing in its own right - ln R must be positive there,
                   not merely better than where the level is now.
      weak         a free position whose support is thin.  Leave it until the
                   region around it is settled - the neighbours will change
                   the answer.
      no support   one line, or none that is free.  No case at all.
    """
    if np.isfinite(row['near_dE']) and abs(
            row['near_dE']) < INTERCHANGE_DE:
        return 'interchange'
    if np.isfinite(row['dE_alt']) and abs(row['dE_alt']) < REFIT_DE:
        return 'refit'
    share = row['top_share']
    if row['n_free'] < AUDIT_MIN_FEW or (np.isfinite(share)
                                         and share > AUDIT_ANY_SHARE):
        return 'no support'
    if (row['n_free'] >= AUDIT_MIN_FREE
            and row['ln_R_alt'] > 0.0
            and row['free_gain'] >= AUDIT_MIN_GAIN
            and np.isfinite(share) and share <= AUDIT_MAX_SHARE
            and row['look'] >= AUDIT_MIN_LOOK):
        return 'relocate'
    return 'weak'


def audit_table(ctx, tab, e_calc, w_of, vac):
    """Add the audit columns to a scanned report table.

    gain     ln R the alternate wins by (positive: the lines prefer it)
    look     that gain after the look-elsewhere correction.  A level allowed
             a 7000 cm^-1 window gets hundreds of chances at a good maximum
             and a level allowed 130 gets a handful, so the two cannot be
             compared raw.  Under H0 the heights of the local maxima are
             roughly exponential, so the best of n_alt of them stands about
             ln(n_alt) above a typical one; requiring the gain to beat that
             puts every level on the same footing.
    n_free   lines supporting the alternate that no other level is using
    free_gain  what those lines are worth
    top_share  the largest row's share of the positive evidence
    n_vacant   unfound calculated levels of the same J and parity within one
             configuration window of the alternate - can an UNKNOWN level be
             there instead?
    z_alt    (E_alt - E_calc)/W: how far the alternate sits from where the
             calculation puts the level, in units of that configuration's own
             scatter
    action   see disposition()
    """
    cols = dict(n_free=[], free_gain=[], top_share=[], gain=[], look=[],
                n_obs_alt=[], n_seen_alt=[], n_miss_alt=[],
                n_vacant=[], z_alt=[], action=[])
    lev = ctx.per.set_index('level_id')
    for _, r in tab.iterrows():
        lid = r['level_id']
        if not r.get('n_alt', 0) or not np.isfinite(r.get('dE_alt', np.nan)):
            for k in cols:
                cols[k].append('' if k == 'action' else np.nan)
            continue
        e_a = float(r['E']) + float(r['dE_alt'])
        s = support(ctx, lid, e_a)
        gain = float(r['ln_R_alt']) - float(r['ln_R'])
        look = gain - math.log(max(int(r['n_alt']), 1))
        w = float(w_of.get(lid, np.nan))
        ec = float(e_calc.get(lid, np.nan))
        row = dict(r)
        row.update(s)
        row['gain'] = gain
        row['look'] = look
        cols['n_free'].append(s['n_free'])
        cols['free_gain'].append(s['free_gain'])
        cols['top_share'].append(s['top_share'])
        cols['n_obs_alt'].append(s['n_obs_alt'])
        cols['n_seen_alt'].append(s['n_seen_alt'])
        cols['n_miss_alt'].append(s['n_miss_alt'])
        cols['gain'].append(gain)
        cols['look'].append(look)
        cols['n_vacant'].append(count_vacancies(
            vac, lev.at[lid, 'J'] if lid in lev.index else None,
            lev.at[lid, 'parity'] if lid in lev.index else None, e_a, w))
        cols['z_alt'].append((e_a - ec) / w if (np.isfinite(ec)
                                                and np.isfinite(w) and w > 0)
                             else np.nan)
        cols['action'].append(disposition(row))
    for k, v in cols.items():
        tab[k] = v
    return tab


ACTION_ORDER = ['relocate', 'interchange', 'refit', 'weak',
                'no support', '']


def order_report(tab):
    """Put the rows worth acting on at the top.

    The report is read from the top down, so the first row should be the move
    with the best case for it, not the level with the worst ln R.  With
    --audit that is the disposition order of ACTION_ORDER, and within each
    disposition the look-elsewhere-corrected gain, largest first.  With --scan
    alone there is no disposition, so it is simply the size of the gain the
    alternate offers.  Without either, the weakest positions come first, which
    is the only ordering a plain report can have.
    """
    tab = tab.copy()
    if 'action' in tab.columns:
        rank = {a: i for i, a in enumerate(ACTION_ORDER)}
        tab['_r'] = [rank.get(str(a), len(ACTION_ORDER))
                     for a in tab['action']]
        tab = tab.sort_values(['_r', 'look', 'gain'],
                              ascending=[True, False, False],
                              na_position='last')
        tab = tab.drop(columns='_r')
    elif 'd_ln_R' in tab.columns:
        tab = tab.sort_values(['d_ln_R', 'ln_R'], ascending=[True, True],
                              na_position='last')
    else:
        tab = tab.sort_values('ln_R')
    return tab.reset_index(drop=True)


def print_audit(tab, alt_drop):
    """The audit report: what to do about the levels that have an alternate."""
    q = tab[tab['action'].astype(str) != ''].copy()
    if q.empty:
        print('\nno level has an alternate position to audit')
        return
    pref = q[q['gain'] > 0]
    print(f"\naudit of the {len(q)} levels with an alternate position within "
          f"{alt_drop:g} of the adopted one")
    print(f"  {len(pref)} of them have an alternate the lines actually prefer")
    print('\n  what kind of problem each one is:')
    for name, g in pref.groupby('action'):
        print(f"    {name:<12} {len(g):>4}")
    firm = pref[pref['action'] == 'relocate'].sort_values(
        'look', ascending=False)
    cols = ['level_id', 'E', 'n_obs', 'n_seen', 'ln_R', 'dE_alt',
            'ln_R_alt', 'n_seen_alt', 'n_miss_alt', 'gain', 'n_alt', 'look',
            'n_free', 'free_gain', 'top_share', 'n_vacant', 'z_alt']
    with pd.option_context('display.width', 220, 'display.max_columns', 24):
        print(f"\n  firm grounds for relocation ({len(firm)}):")
        print(firm[cols].to_string(index=False) if len(firm) else '    none')
        swap = pref[pref['action'] == 'interchange']
        if len(swap):
            print(f"\n  interchanges, for level_interchange.py ({len(swap)}):")
            print(swap[cols + ['near_level', 'near_dE']].to_string(
                index=False))
        rf = pref[pref['action'] == 'refit']
        if len(rf):
            print(f"\n  the level stays; an accepted line is dragging the fit "
                  f"({len(rf)}):")
            print(rf[cols].to_string(index=False))
    hole = pref[(pref['action'] == 'relocate') & (pref['n_vacant'] > 0)]
    if len(hole):
        print(f"\n  {len(hole)} of the relocations sit where the calculation "
              f"still has an unfound level of the same J and parity, so the "
              f"lines may belong to THAT level rather than to this one: "
              f"{', '.join(hole['level_id'])}")
    print('\n  A relocation changes the lines available to its neighbours, so '
          'accept them\n  one at a time and re-scan: the verdict on every '
          'level is conditional on the\n  rest of the list.')


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--detail', metavar='LEVEL_ID',
                   help='per-transition table for one level')
    p.add_argument('--scan', action='store_true',
                   help='scan the alternate-position window of every level')
    p.add_argument('--audit', action='store_true',
                   help='what each alternate position rests on and what to '
                        'do about it (implies --scan)')
    p.add_argument('--levels', nargs='*', default=None,
                   help='restrict the report to these level ids')
    p.add_argument('--alt-drop', type=float, default=ALT_DROP,
                   help='ln R below the adopted peak that still counts as an '
                        'alternate position (default %(default)s)')
    p.add_argument('--step', type=float, default=GRID_STEP,
                   help='scan step in cm^-1 (default %(default)s)')
    p.add_argument('--lopt', default=None,
                   help='build the run from a LOPT line-output file')
    p.add_argument('--lopt-levels', default=None)
    p.add_argument('--energies', default='refit', choices=['refit', 'lopt'])
    p.add_argument('--energies-csv', default=None,
                   help='csv of revised adopted energies (level_id,E_input)')
    p.add_argument('--out', default='level_positions.csv')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    ctx = build(args)

    if args.detail:
        print_detail(ctx, args.detail)
        return 0

    ids = args.levels if args.levels else list(ctx.per['level_id'])
    ids = [i for i in ids if i in ctx.by_level]
    tab = report_table(ctx, ids)

    if args.scan or args.audit:
        en = li.read_enlev()
        win = li.configuration_windows(en)
        lv, unmatched = li.attach_identities(ctx.per, en, win)
        e_calc = dict(zip(lv['level_id'], lv['E_calc']))
        w_of = dict(zip(lv['level_id'], lv['W']))
        print(f"\nscanning {len(ids)} levels "
              f"({len(unmatched)} matched no row of enlev.dat)")
        n_alt, best_alt, d_alt, sep, width = [], [], [], [], []
        # a maximum that falls on another level of the run is a different
        # question from one that falls on empty energy: the first is an
        # interchange, which level_interchange.py judges on evidence this
        # scan does not use, the second a position nobody has claimed
        lev_e = ctx.per['E_final'].to_numpy(dtype=float)
        lev_id = ctx.per['level_id'].to_numpy()
        o = np.argsort(lev_e)
        lev_e, lev_id = lev_e[o], lev_id[o]
        alt_lev, alt_lev_dE = [], []
        for lid in tab['level_id']:
            r = scan_level(ctx, lid, ctx.e_final[lid],
                           e_calc.get(lid, np.nan), w_of.get(lid, np.nan),
                           step=args.step, alt_drop=args.alt_drop)
            n_alt.append(len(r['alternates']))
            width.append(r['hi'] - r['lo'])
            if r['alternates']:
                e_a, v_a = r['alternates'][0]
                best_alt.append(v_a)
                d_alt.append(r['ln_R'] - v_a)
                sep.append(e_a - r['e_adopted'])
                k = int(np.clip(np.searchsorted(lev_e, e_a), 1,
                                len(lev_e) - 1))
                if abs(lev_e[k] - e_a) >= abs(lev_e[k - 1] - e_a):
                    k -= 1
                alt_lev.append(lev_id[k])
                alt_lev_dE.append(float(e_a - lev_e[k]))
            else:
                best_alt.append(np.nan)
                d_alt.append(np.nan)
                sep.append(np.nan)
                alt_lev.append('')
                alt_lev_dE.append(np.nan)
        tab['scan_width'] = width
        tab['n_alt'] = n_alt
        tab['ln_R_alt'] = best_alt
        tab['d_ln_R'] = d_alt
        tab['dE_alt'] = sep
        tab['near_level'] = alt_lev
        tab['near_dE'] = alt_lev_dE
        tab['question'] = np.where(np.asarray(n_alt) > 0, '?', '')
        if args.audit:
            tab = audit_table(ctx, tab, e_calc, w_of, vacancies(en))

    tab = order_report(tab)
    mc.save_table(tab, args.out, decimals={
        'E': 4, 'ln_R': 3, 'ln_R_match': 3, 'ln_R_miss': 3,
        'sum_lnG': 3, 'ln_R_alt': 3, 'd_ln_R': 3, 'dE_alt': 4,
        'near_dE': 4, 'scan_width': 1, 'gain': 3, 'look': 3,
        'free_gain': 3, 'top_share': 3, 'z_alt': 2})
    print(f"ln R at the adopted position: median {tab['ln_R'].median():.1f}, "
          f"{int((tab['ln_R'] <= 0).sum())} at or below zero, "
          f"{int((tab['ln_R'] < 10).sum())} below 10")
    print(f"predictions dropped for partner degeneracy: "
          f"{int(tab['n_drop'].sum())}")
    if 'n_alt' in tab:
        q = tab[tab['n_alt'] > 0]
        print(f"\nlevels with an alternate position within "
              f"{args.alt_drop:g} of the adopted one: {len(q)} of {len(tab)}")
        print(f"  of these, {int((q['d_ln_R'] < 0).sum())} have an alternate "
              f"the lines prefer to the adopted position")
        onlev = (q['near_dE'].abs() < 0.5).sum()
        print(f"  {int(onlev)} of the best alternates fall on another level "
              f"of the run - an interchange, not a free position")
        cols = ['level_id', 'E', 'n_obs', 'n_seen', 'n_miss', 'ln_R',
                'n_alt', 'ln_R_alt', 'd_ln_R', 'dE_alt', 'near_level',
                'near_dE']
        with pd.option_context('display.width', 200,
                               'display.max_columns', 24):
            print(q.sort_values('d_ln_R')[cols].head(40).to_string(index=False))
    if 'action' in tab.columns:
        print_audit(tab, args.alt_drop)
    with pd.option_context('display.width', 200,
                           'display.max_columns', 24):
        print("\nthe twenty weakest positions:")
        print(tab.sort_values('ln_R').head(20).to_string(index=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
