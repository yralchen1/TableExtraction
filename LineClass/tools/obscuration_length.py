"""How wide is the thing that hides a line: the correlation length of the
residual obscuration.

tools/obscuration_rate.py measures epsilon, the probability that a predicted
line certainly bright enough to be recorded, at a wavelength the coverage map
calls covered, was nevertheless not recorded: about one line in thirty, the
same everywhere.  That number is a rate and not a map - it says how often a
line is lost, never where.  For a single absent branch a rate is all the
likelihood needs.  For two it is not, and the reason is the whole subject of
this tool.

1. THE QUESTION

A level whose strongest predicted branch is missing is suspicious; a level
with two strong branches missing is more suspicious still - but how much more?
If the two absences are independent events the penalty multiplies, epsilon
squared, one chance in eight hundred.  If one grain flaw took both, the
penalty is paid once, one chance in thirty.  Between those two readings lies a
factor of thirty in the evidence, which is enough to decide a level.

What separates them is a length.  An obscuring agent - an impurity line, a
scratch, a flaw in the emulsion, a plate seam - occupies a stretch of the
plate, and two predicted lines that fall inside the same stretch are lost
together.  Call that stretch's width L.  Two branches closer than L are one
absence; two branches farther apart than L are two.  So the quantity to
measure is

    phi(d) = P(both obscured | separation d) - epsilon^2, in units of
             epsilon(1 - epsilon)

the correlation coefficient of the obscuration along the spectrum: 1 when the
two positions always share the fate of a single event, 0 when their fates are
independent.  L is the width of phi.

2. THE MEASUREMENT

The population is the one obscuration_rate.py uses, and for the same reason:
the 13709 predicted transitions between two ESTABLISHED levels - supported by
more old than new identifications and present in the ASD compilation - whose
energies are certain to a few thousandths of a wavenumber, so that an absence
at their Ritz wavelength is the plate's doing and not the identification's.
Every pair of predictions closer than --dmax angstroms is examined, 276567 of
them out to 5 A, and asked whether BOTH are absent, where absent means what it
means in obscuration_rate.py: no measured line within the matching window and
no recorded line close and strong enough to have swallowed it.

Joint absence has to be compared with something, and three separate effects
make nearby predictions look jointly absent when nothing is hiding them.

*The line is faint.*  Most predicted transitions are far below the noise, and
two faint predictions are jointly absent for no interesting reason.  So each
prediction carries a modelled absence probability q, fitted to the run itself:
a logistic model of the outcome in the predicted strength ln(I_pred / noise)
and the wavelength, of the form P(explained) = f + (1 - f) c D, with f the
prediction's own chance rate and c its coverage, and then raked so that the
modelled and observed absence rates agree in every one of forty wavelength
blocks and ten strength bins.  The comparison is between the observed joint
absences and the sum of q_i q_j over the same pairs.

*The two positions see the same piece of plate by coincidence.*  With a
matching window of 0.022 A and, in the crowded ultraviolet, a recorded line
every 0.13 A, a pair of predictions a hundredth of an angstrom apart is
explained or not explained together about one time in six by pure chance.
That is not obscuration, and it is removed the way obscuration_rate.py removes
the same effect from epsilon: the whole pair analysis is repeated at the 18
control positions, each prediction displaced by +-6, 9, ... 30 window widths,
which preserves every pair's separation, both members' predicted strengths and
the local line density, and changes only that the position is no longer one
where a line was expected.  The controls show the artefact plainly - joint
absence 1.18 times the independent value below 0.02 A, falling to 1.00 beyond
0.05 A - and dividing by it removes it.

*The baseline is imperfect.*  A logistic model of a photographic plate is not
exact, and any error in it that varies slowly with wavelength makes nearby
predictions look correlated at every separation.  The remedy is to let the fit
follow the wavelength closely enough (forty blocks, about 250 A each) that the
long-range correlation disappears, and to check that it has: the raw
correlation of the absences falls from 0.09 at every separation out to 20 A,
with a baseline that has no wavelength dependence at all, to 0.00 with the
fitted one.  What is left of it is reported as a floor and fitted as a
nuisance parameter, never assumed to be zero.

3. WHY THE ESTIMATOR LOOKS LIKE THIS

Write A_i for absent, O_i for obscured (probability epsilon) and F_i for
'too faint or otherwise lost, independently' (probability d_i), so that
A = O or F and q_i = epsilon + (1 - epsilon) d_i.  Then

    P(A_i A_j) - q_i q_j = epsilon (1 - epsilon) phi(d) (1 - d_i)(1 - d_j)
                           + O(epsilon^2)

and since (1 - d_i) = (1 - q_i)/(1 - epsilon),

    phi(d) = (1 - epsilon)/epsilon * [ sum A_iA_j - Rc(d) sum q_iq_j ]
                                   / [ sum (1 - q_i)(1 - q_j) ]

with Rc(d) the chance factor from the controls.  Two things follow from the
shape of this formula and both matter.

The first is that no strength cut is needed.  A pair of faint predictions has
(1 - q_i)(1 - q_j) near zero, so it carries almost no weight in the answer -
the formula down-weights it automatically, and correctly, because a faint
prediction's absence says nothing about the plate.  Pairs are weighted by the
inverse of their variance, (1 - q_i)(1 - q_j) / q_i q_j, which is another way
of saying that the measurement is made on the strong pairs without having to
choose a threshold for them.

The second is that the factor (1 - epsilon)/epsilon = 28 amplifies everything,
the systematic errors of the baseline included.  A residual error of half a
percent in the modelled joint absence rate appears as phi = 0.14.  That, and
not the counting statistics, is what limits this measurement, and it is why
the floor is fitted rather than assumed.

4. THE RESULT

Fitting phi(d) = (1 - d/L)+ , the overlap of two positions with one obscuring
stretch of width L, over all pairs closer than 1 A, with the floor free:

    L = 0.07 A, 95 % interval 0.03 - 0.15 A

against no correlation at all by 2 dNLL = 14, and unchanged if the pair range
is 0.5, 3 or 5 A or if epsilon is set anywhere in its own interval.  In the
crowded ultraviolet where nearly all the close pairs lie, 0.07 A is about two
resolution elements: the effective line width there, the instrumental 0.035 A
convolved with the Doppler width, is 0.036 A.  The obscuration is, as far as
this measurement can see, the width of a line and no wider.

The amplitude of the clustering is itself a second measurement of epsilon,
made from a different statistic than the rate: fitting epsilon and L together
to the pairs alone gives epsilon = 0.069 (0.037 - 0.145) against the rate
measurement's 0.035 (0.018 - 0.061).  The intervals overlap; the central
values differ by a factor of two, in the direction of MORE obscuration among
the close pairs, which is where it would be expected, since almost all of them
lie between 1000 and 2800 A and the rate measurement puts epsilon there at
0.054 rather than 0.022.

5. WHAT IT MEANS FOR THE LIKELIHOOD

Almost always, nothing - and that is the useful part of the answer.  The
branches of one level go to different lower levels, and different lower levels
are hundreds or thousands of wavenumbers apart, so the branches are far apart
on the plate.  Among predictions at least a hundred times the noise level, the
closest two branches of any established level ever come is 0.446 A, six times
L; the median separation is 206 A.  Absences of one level's branches may
therefore be multiplied as independent events, which is what the missing-line
term already does.

The exception is worth a warning in the code that uses it, because when it
does occur it is large.  Two predicted branches of one level closer than about
0.1 A are one test, not two: their joint absence costs epsilon, not epsilon
squared, and treating them as independent overstates the case against the
level by a factor of 1/epsilon = 29.  Among predictions ten times the noise or
stronger there are five such within-level pairs, and this tool lists them.

6. WHAT THIS IS NOT

It is not a measurement of how obscuration is distributed along the plate at
scales above an angstrom.  The systematic floor described in section 3 sets
the sensitivity at |phi| of about 0.15, so a weak large-scale modulation - a
plate that was slightly fogged over tens of angstroms, say - would not be seen
here, and the coverage function in coverage_map.csv is the only handle on
structure at that scale.

Nor does it separate the agents.  An impurity line, a grain flaw and a scratch
all produce the same statistic, and the measured width is close enough to the
resolution element that the data cannot distinguish an emulsion defect from a
line of some other species too weak to have been measured.

Like tools/obscuration_rate.py, and unlike tools/coverage_map.py, this is not
a standalone tool: it reads the run through level_shifts.py and chance_mc.py,
and it reuses obscuration_rate.py for the population, the matching window and
the outcomes, so that the two measurements are made on exactly the same
predictions.

Outputs (in the working directory unless --out-dir is given):
  obscuration_length.txt   the fitted length, phi(d) in bins, the control
                           check, the within-level separations, and the
                           jointly absent close pairs
  obscuration_length.csv   one row per close pair of predictions
  obscuration_length.png   phi(d) with the fit, the control artefact, and the
                           within-level separations against L
  obscuration_length.log   the console output
"""
import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import level_shifts as ls                   # noqa: E402
import lopt_lines                           # noqa: E402
import obscuration_rate as orate            # noqa: E402


# --- the outcomes, real and control -----------------------------------------
def outcomes(opt, log):
    """The predictions, their real outcome, and their outcome at every control.

    Everything here is obscuration_rate.py's, called function by function so
    that the two tools cannot drift apart: the same established levels, the
    same matching window, the same masking test, the same control offsets.
    The one thing kept that obscuration_rate.py throws away is the outcome at
    each individual control position, which is what the pair analysis needs.
    """
    P, _, _, _ = orate.build_predictions(opt, log)

    obs = lopt_lines.read_observed_lines()
    obs['lam'] = 1.0e8 / obs['wn_obs']
    obs = obs.sort_values('wn_obs').reset_index(drop=True)
    wn_o = obs['wn_obs'].to_numpy()
    int_o = np.nan_to_num(obs['obs_intens'].to_numpy())
    lam_o = np.sort(1.0e8 / wn_o)
    log(f"measured lines: {len(obs)}, {lam_o[0]:.2f} - {lam_o[-1]:.2f} A")

    W = orate.matching_window(P, obs, opt.k_window, log)
    P['W'] = W
    dop = ls.doppler_fwhm_factor()
    lam = P['lam'].to_numpy()
    ip = P['I_pred'].to_numpy()
    P['recorded'] = orate.is_recorded(lam, W, lam_o)
    P['masked'] = orate.is_masked(lam, ip, wn_o, int_o, dop)
    P['explained'] = P['recorded'] | P['masked']
    P['c'] = orate.attach_coverage(P, opt.coverage, log)

    offsets = np.concatenate([np.asarray(opt.offsets, dtype=float),
                              -np.asarray(opt.offsets, dtype=float)])
    Y = np.zeros((len(offsets), len(P)), dtype=bool)
    for k, m in enumerate(offsets):
        pos = lam + m * W
        Y[k] = (orate.is_recorded(pos, W, lam_o)
                | orate.is_masked(pos, ip, wn_o, int_o, dop))
    P['f_chance'] = Y.mean(axis=0)
    log(f"{len(offsets)} control positions per prediction: explained "
        f"{Y.mean():.4f} of the time, against {P['explained'].mean():.4f} at "
        f"the real positions")

    # everything downstream walks the predictions in wavelength order, and the
    # control outcomes have to be carried along in the same order
    order = np.argsort(P['lam'].to_numpy())
    return P.iloc[order].reset_index(drop=True), Y[:, order], offsets


# --- the modelled absence probability ---------------------------------------
def hat_basis(x: np.ndarray, k: int) -> np.ndarray:
    """Piecewise-linear bumps on k quantile knots of x.

    A spline basis with nothing clever about it: each column rises from zero
    at its neighbouring knots to one at its own.  Quantile knots put the
    resolution where the predictions are.
    """
    kn = np.unique(np.quantile(x, np.linspace(0.0, 1.0, k)))
    B = np.zeros((len(x), len(kn)))
    for i, t in enumerate(kn):
        if i == 0:
            w = kn[1] - kn[0]
        elif i == len(kn) - 1:
            w = kn[-1] - kn[-2]
        else:
            w = max(kn[i + 1] - kn[i], kn[i] - kn[i - 1])
        B[:, i] = np.clip(1.0 - np.abs(x - t) / w, 0.0, None)
    return B


def fit_baseline(P: pd.DataFrame, opt, log) -> np.ndarray:
    """q_i, the probability that prediction i would be absent on its own.

    The link is the structural one of obscuration_rate.py,
    P(explained) = f + (1 - f) c D, with f the prediction's chance rate and c
    its coverage, and D a logistic function of the predicted strength and the
    wavelength.  Nothing here is meant to be a physical model of detection -
    the measured detection curve in obscuration_rate.txt is that - only a
    smooth interpolation good enough that the residual absences carry no
    structure of their own.
    """
    from scipy.optimize import minimize
    y = P['explained'].to_numpy().astype(float)
    z = P['z'].to_numpy()
    f = P['f_chance'].to_numpy()
    c = P['c'].to_numpy()
    lam = P['lam'].to_numpy()

    Bz = hat_basis(z, opt.knots_z)
    Bl = hat_basis(np.log(lam), opt.knots_lam)
    cols = [np.ones((len(z), 1)), Bz, Bl]
    edges = np.quantile(lam, np.linspace(0.0, 1.0, opt.blocks_z + 1))[1:-1]
    lb = np.digitize(lam, edges)
    for b in range(opt.blocks_z):
        cols.append(Bz * (lb == b)[:, None])
    X = np.hstack(cols)

    def negll(beta):
        g = 1.0 / (1.0 + np.exp(-(X @ beta)))
        p = np.clip(f + (1.0 - f) * c * g, 1e-9, 1.0 - 1e-9)
        d = (y / p - (1.0 - y) / (1.0 - p)) * (1.0 - f) * c * g * (1.0 - g)
        return (-(y * np.log(p) + (1.0 - y) * np.log(1.0 - p)).sum()
                + 1e-3 * beta @ beta,
                -(X * d[:, None]).sum(axis=0) + 2e-3 * beta)

    r = minimize(negll, np.zeros(X.shape[1]), jac=True, method='L-BFGS-B',
                 options=dict(maxiter=8000, maxfun=40000))
    g = 1.0 / (1.0 + np.exp(-(X @ r.x)))
    q = 1.0 - np.clip(f + (1.0 - f) * c * g, 1e-9, 1.0 - 1e-9)
    log(f"baseline: {X.shape[1]} parameters, mean modelled absence "
        f"{q.mean():.4f} against {1 - y.mean():.4f} observed")
    return q


def rake(q, a, lam, s, nlam, nstr, rounds=8):
    """Scale q inside wavelength blocks and strength bins until its mean
    matches the observed absence rate in each.

    The fitted model is close but not exact - the structural link does not
    force the mean to come out right - and the correlation estimator amplifies
    whatever is left by 1/epsilon.  Raking removes it in the two directions
    that matter.  Each cell holds a few hundred predictions, so forcing its
    mean removes at most a percent of any real correlation with it.
    """
    q = q.copy()
    lb = np.digitize(lam, np.quantile(lam, np.linspace(0, 1, nlam + 1))[1:-1])
    sb = np.digitize(s, np.quantile(s, np.linspace(0, 1, nstr + 1))[1:-1])
    for _ in range(rounds):
        for cell in (lb, sb):
            for v in np.unique(cell):
                m = cell == v
                q[m] = np.clip(q[m] * a[m].mean() / max(q[m].mean(), 1e-9),
                               1e-6, 1.0 - 1e-6)
    return q


# --- pairs ------------------------------------------------------------------
def build_pairs(lam: np.ndarray, dmax: float, log):
    """Every pair of predictions closer than dmax angstroms, lam sorted."""
    I, J = [], []
    n = len(lam)
    for i in range(n):
        j = i + 1
        while j < n and lam[j] - lam[i] <= dmax:
            I.append(i)
            J.append(j)
            j += 1
    I = np.asarray(I)
    J = np.asarray(J)
    log(f"pairs of predictions within {dmax:g} A: {len(I)}")
    return I, J, lam[J] - lam[I]


def chance_factor(D, AC, QC, edges):
    """Rc(d): how much more often two control positions are jointly absent
    than their own marginals imply.  Everything not obscuration that couples
    two nearby tests is in here."""
    Rc = np.ones(len(D))
    k = np.searchsorted(edges, D, side='right') - 1
    table = []
    for b in range(len(edges) - 1):
        m = k == b
        if not m.any():
            continue
        v = AC[m].sum() / max(QC[m].sum(), 1e-12)
        Rc[m] = v
        table.append((edges[b], edges[b + 1], int(m.sum()), v))
    return Rc, table


# --- the estimator ----------------------------------------------------------
def phi_binned(D, AA, QQ, UU, Rc, edges, eps, sel=None):
    """phi(d) in bins, inverse-variance weighted within each bin."""
    k = np.searchsorted(edges, D, side='right') - 1
    out = np.full(len(edges) - 1, np.nan)
    w = UU / np.maximum(QQ, 1e-9)
    base = np.ones(len(D), dtype=bool) if sel is None else sel
    for b in range(len(edges) - 1):
        m = base & (k == b)
        if not m.any():
            continue
        num = (w[m] * (AA[m] - Rc[m] * QQ[m])).sum()
        den = (w[m] * UU[m]).sum()
        out[b] = (1.0 - eps) / eps * num / max(den, 1e-12)
    return out


def fit_length(D, AA, QQ, UU, Rc, eps, grid_L, grid_f, dfit, sel=None):
    """Profile the top-hat width L over the floor, on pairs closer than dfit.

    phi(d) = (1 - d/L)+ is the probability that one obscuring stretch of width
    L, dropped at random, covers both positions - the simplest model with a
    length in it and the only one the data can distinguish from no correlation
    at all.  The floor is a constant added to phi: whatever correlation the
    baseline has failed to remove, fitted rather than assumed away.
    """
    m = D <= dfit if sel is None else (sel & (D <= dfit))
    a, b, u, dd = AA[m], (Rc * QQ)[m], UU[m], D[m]
    K = eps / (1.0 - eps)
    S = np.empty((len(grid_L), len(grid_f)))
    for i, L in enumerate(grid_L):
        phi = np.clip(1.0 - dd / L, 0.0, None)
        for j, f0 in enumerate(grid_f):
            p = np.clip(b + K * (phi + f0) * u, 1e-9, 1.0 - 1e-9)
            S[i, j] = -(a * np.log(p) + (1.0 - a) * np.log(1.0 - p)).sum()
    i, j = np.unravel_index(S.argmin(), S.shape)
    prof = S.min(axis=1)
    ok = grid_L[prof <= prof.min() + 1.92]
    null = min(-(a * np.log(np.clip(b + K * f0 * u, 1e-9, 1 - 1e-9))
                 + (1.0 - a) * np.log(1.0 - np.clip(b + K * f0 * u, 1e-9,
                                                    1 - 1e-9))).sum()
               for f0 in grid_f)
    return grid_L[i], ok.min(), ok.max(), grid_f[j], 2.0 * (null - S.min())


def fit_eps_from_clustering(D, AA, QQ, UU, Rc, grid_e, grid_L, grid_f, dfit):
    """epsilon and L together, from the pairs alone.

    The amplitude of the excess joint absence is proportional to epsilon, so
    the clustering measures it a second time, by a statistic that has nothing
    in common with the rate measurement but the outcomes themselves.
    """
    m = D <= dfit
    a, b, u, dd = AA[m], (Rc * QQ)[m], UU[m], D[m]
    S = np.full((len(grid_e), len(grid_L)), np.inf)
    for ie, eps in enumerate(grid_e):
        K = eps / (1.0 - eps)
        for il, L in enumerate(grid_L):
            phi = np.clip(1.0 - dd / L, 0.0, None)
            core = K * phi * u
            for f0 in grid_f:
                p = np.clip(b + core + f0 * u, 1e-9, 1.0 - 1e-9)
                v = -(a * np.log(p) + (1.0 - a) * np.log(1.0 - p)).sum()
                if v < S[ie, il]:
                    S[ie, il] = v
    ie, il = np.unravel_index(S.argmin(), S.shape)
    pe = S.min(axis=1)
    ok = grid_e[pe <= pe.min() + 1.92]
    return grid_e[ie], ok.min(), ok.max(), grid_L[il]


# --- what it means in practice ----------------------------------------------
def within_level_separations(P: pd.DataFrame, cuts, cmin):
    """How far apart the predicted branches of one level fall.

    The answer decides whether the correlation length is ever reached in
    practice.  Both the lower and the upper level of a transition define a set
    of branches, and both are counted.
    """
    rows = []
    for cut in cuts:
        Q = P[(P['I_pred'] >= cut * P['thr']) & (P['c'] >= cmin)]
        seps = []
        for _, g in itertools.chain(Q.groupby('lo_id'), Q.groupby('up_id')):
            v = np.sort(g['lam'].to_numpy())
            if len(v) > 1:
                seps.append(np.diff(v))
        seps = np.concatenate(seps) if seps else np.array([np.inf])
        rows.append((cut, len(Q), len(seps) if np.isfinite(seps).all() else 0,
                     seps.min(), float(np.median(seps)), seps))
    return rows


def close_branch_pairs(P, cut, cmin, L):
    """Branches of one level closer than L: the pairs the likelihood must
    treat as one test rather than two."""
    Q = P[(P['I_pred'] >= cut * P['thr']) & (P['c'] >= cmin)]
    rows = []
    for kind, key in (('lower', 'lo_id'), ('upper', 'up_id')):
        for lev, g in Q.groupby(key):
            g = g.sort_values('lam')
            lam = g['lam'].to_numpy()
            for i in range(len(lam) - 1):
                if lam[i + 1] - lam[i] < L:
                    r1, r2 = g.iloc[i], g.iloc[i + 1]
                    rows.append((str(lev), kind, r1['lam'], r2['lam'],
                                 lam[i + 1] - lam[i],
                                 r1['I_pred'] / r1['thr'],
                                 r2['I_pred'] / r2['thr'],
                                 bool(r1['explained']), bool(r2['explained'])))
    return rows


def close_absent_pairs(P, I, J, D, dmax, smin):
    """Pairs of absent predictions, both above smin x noise, within dmax."""
    s = (P['I_pred'] / P['thr']).to_numpy()
    a = (~P['explained'].to_numpy())
    lam = P['lam'].to_numpy()
    m = (D <= dmax) & a[I] & a[J] & (s[I] >= smin) & (s[J] >= smin)
    T = pd.DataFrame({
        'lam_1': lam[I[m]], 'lam_2': lam[J[m]], 'dlam': D[m],
        'strength_1': s[I[m]], 'strength_2': s[J[m]],
        'coverage': P['c'].to_numpy()[I[m]],
        'same_level': ((P['lo_id'].to_numpy()[I[m]]
                        == P['lo_id'].to_numpy()[J[m]])
                       | (P['up_id'].to_numpy()[I[m]]
                          == P['up_id'].to_numpy()[J[m]])),
    })
    return T.sort_values('dlam').reset_index(drop=True)


# --- output -----------------------------------------------------------------
def write_report(path, opt, eps, res, phi, phi_err, ctrl, wl, pairs, extra,
                 log):
    L, Llo, Lhi, floor, sig = res
    with open(path, 'w', newline='\n') as fh:
        w = fh.write
        w("# correlation length of the residual obscuration\n")
        w("# produced by tools/obscuration_length.py\n#\n")
        w("# phi(d) = correlation coefficient of the obscuration at\n")
        w("#          separation d: 1 if two positions d apart always share\n")
        w("#          the fate of one obscuring event, 0 if their fates are\n")
        w("#          independent.  L is the width of phi.\n#\n")
        w(f"# population: {extra['n_pred']} predicted transitions between\n")
        w("#             two established levels, all pairs closer than "
          f"{opt.dmax:g} A\n\n")
        w(f"epsilon_used     {eps:.4f}   # from the same run, "
          "tools/obscuration_rate.py\n")
        w(f"n_pairs          {extra['n_pairs']}\n")
        w(f"n_pairs_fitted   {extra['n_fit']}   # separation <= "
          f"{opt.dfit:g} A\n")
        w(f"L                {L:.4f}   # angstrom, top-hat width\n")
        w(f"L_lo95           {Llo:.4f}\n")
        w(f"L_hi95           {Lhi:.4f}\n")
        w(f"L_leave_one_out  {extra['L_lo']:.4f} - {extra['L_hi']:.4f}"
          "   # refitting without any one wavelength block\n")
        w(f"floor_phi0       {floor:+.3f}   # unremoved baseline "
          "correlation; the method's systematic\n")
        w(f"significance     {sig:.1f}   # 2 dNLL against no correlation, "
          "1 degree of freedom\n")
        w(f"fwhm_at_1500A    {extra['fwhm']:.4f}   # effective resolution "
          "there, for comparison with L\n\n")
        w("# epsilon measured a second time, from the clustering amplitude\n")
        w(f"epsilon_pairs    {extra['eps_pair']:.4f}\n")
        w(f"epsilon_pairs_lo {extra['eps_lo']:.4f}\n")
        w(f"epsilon_pairs_hi {extra['eps_hi']:.4f}\n\n")

        w("# phi(d) in bins, inverse-variance weighted; the error is a\n")
        w("# jackknife over the wavelength blocks.  The floor above has NOT\n")
        w("# been subtracted here.\n")
        w("#    d_lo     d_hi    pairs       phi     error\n")
        for (lo, hi), p, e, n in zip(zip(opt.bins[:-1], opt.bins[1:]),
                                     phi, phi_err, extra['bin_n']):
            w(f"  {lo:8.3f} {hi:8.3f} {n:8d} {p:9.3f} {e:9.3f}\n")

        w("\n# the control check: the same statistic at the displaced\n")
        w("# positions, where nothing was expected.  Its excess below one\n")
        w("# window width is the coincidence artefact, and it is divided\n")
        w("# out of phi above.\n")
        w("#    d_lo     d_hi    pairs   R_control\n")
        for lo, hi, n, v in ctrl:
            w(f"  {lo:8.3f} {hi:8.3f} {n:8d} {v:11.3f}\n")

        w("\n# how far apart the predicted branches of one level fall, on\n")
        w("# which the whole practical consequence turns\n")
        w("#  I/noise   predictions   branch pairs   closest   median   "
          f"closer than L={L:.3f} A\n")
        for cut, npred, nsep, mn, med, seps in wl:
            w(f"  {cut:8.0f} {npred:13d} {nsep:14d} {mn:9.3f} {med:8.1f}"
              f" {int((seps < L).sum()):12d}\n")

        w("\n# branches of one level closer than L, at "
          f"{extra['branch_cut']:g} x the noise level or stronger:\n")
        w("# these are one test and not two, and their joint absence costs\n")
        w("# epsilon, not epsilon squared\n")
        if extra['branch_pairs']:
            w("#      level    via     lam_1      lam_2     dlam   "
              "I/noise_1   I/noise_2   outcome\n")
            for lev, kind, l1, l2, dd, s1, s2, e1, e2 in extra['branch_pairs']:
                w(f"  {lev:>12s} {kind:>6s} {l1:10.4f} {l2:10.4f} {dd:8.4f} "
                  f"{s1:11.1f} {s2:11.1f}   "
                  f"{'seen' if e1 else 'ABSENT'}/"
                  f"{'seen' if e2 else 'ABSENT'}\n")
        else:
            w("#   none\n")

        w("\n# jointly absent pairs, both at least "
          f"{opt.pair_strength:g} x the noise level, within "
          f"{opt.pair_window:g} A\n")
        w("#     lam_1      lam_2     dlam   I/noise_1   I/noise_2  cov  "
          "same level\n")
        for r in pairs.itertuples(index=False):
            w(f"  {r.lam_1:10.4f} {r.lam_2:10.4f} {r.dlam:8.4f} "
              f"{r.strength_1:11.1f} {r.strength_2:11.1f} "
              f"{r.coverage:5.3f}  {'yes' if r.same_level else 'no'}\n")
    log(f"wrote {os.path.basename(path)}")


def plot(path, opt, res, phi, phi_err, ctrl, wl, log):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    L, Llo, Lhi, floor, _ = res
    mid = 0.5 * (np.asarray(opt.bins[:-1]) + np.asarray(opt.bins[1:]))
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))

    a = ax[0]
    a.axhline(0.0, color='0.7', lw=0.8)
    a.axhline(floor, color='tab:orange', lw=0.9, ls=':',
              label=f'fitted floor {floor:+.2f}')
    a.errorbar(mid, phi, yerr=phi_err, fmt='o', ms=4, color='tab:blue',
               capsize=2, label='measured')
    d = np.concatenate([[0.0], np.exp(np.linspace(np.log(1e-3),
                                                  np.log(opt.bins[-1]), 400))])
    a.plot(d, np.clip(1.0 - d / L, 0.0, None) + floor, color='tab:red',
           lw=1.4, label=f'top hat, L = {L:.3f} A')
    a.axvspan(Llo, Lhi, color='tab:red', alpha=0.10, lw=0)
    a.set_xscale('log')
    a.set_xlabel('separation of two predictions, A')
    a.set_ylabel('phi, obscuration correlation')
    a.set_title('the obscuration is one resolution element wide')
    a.set_ylim(-0.8, 2.2)
    a.legend(fontsize=8)

    a = ax[1]
    lo = np.array([c[0] for c in ctrl])
    hi = np.array([c[1] for c in ctrl])
    v = np.array([c[3] for c in ctrl])
    a.axhline(1.0, color='0.7', lw=0.8)
    a.step(np.append(lo, hi[-1]), np.append(v, v[-1]), where='post',
           color='tab:green', lw=1.4)
    a.set_xscale('log')
    a.set_xlabel('separation, A')
    a.set_ylabel('joint absence / independent, at the controls')
    a.set_title('what coincidence alone couples: removed from phi')

    a = ax[2]
    for (cut, _, _, _, _, seps), col in zip(wl, ['0.75', 'tab:blue',
                                                 'tab:red', 'tab:green']):
        s = seps[np.isfinite(seps) & (seps > 0)]
        a.hist(s, bins=np.logspace(np.log10(max(s.min(), 1e-3)),
                                   np.log10(s.max()), 60),
               histtype='step', color=col, label=f'I >= {cut:g} x noise')
    a.axvline(L, color='k', lw=1.2, ls='--', label=f'L = {L:.3f} A')
    a.set_xscale('log')
    a.set_xlabel('separation of two branches of one level, A')
    a.set_ylabel('branch pairs')
    a.set_title('the correlation is essentially never reached')
    a.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    log(f"wrote {os.path.basename(path)}")


# --- driver -----------------------------------------------------------------
def parse_args(argv):
    p = argparse.ArgumentParser(
        prog='obscuration_length.py',
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
    g = p.add_argument_group('the measurement')
    g.add_argument('--dmax', type=float, default=5.0,
                   help='largest separation, A, at which pairs are formed')
    g.add_argument('--dfit', type=float, default=1.0,
                   help='separation, A, out to which the length is fitted')
    g.add_argument('--epsilon', type=float, default=None,
                   help='obscuration rate; by default measured on this run '
                        'the way tools/obscuration_rate.py measures it')
    g.add_argument('--bins', type=float, nargs='+',
                   default=[0.0, 0.015, 0.03, 0.05, 0.08, 0.12, 0.18, 0.28,
                            0.45, 0.8, 1.5, 3.0, 5.0],
                   help='separation bins for phi(d), A')
    g.add_argument('--offsets', type=float, nargs='+',
                   default=[6, 9, 12, 15, 18, 21, 24, 27, 30],
                   help='control displacements, in window widths, used both '
                        'ways')
    g.add_argument('--k-window', type=float, default=5.5,
                   help='matching window in units of the local wavelength '
                        'uncertainty')
    g = p.add_argument_group('the baseline')
    g.add_argument('--knots-z', type=int, default=16,
                   help='spline knots in ln(I_pred / noise)')
    g.add_argument('--knots-lam', type=int, default=40,
                   help='spline knots in ln(lambda); enough of them that no '
                        'correlation survives at large separation')
    g.add_argument('--blocks-z', type=int, default=4,
                   help='wavelength blocks with a strength spline of their '
                        'own')
    g.add_argument('--rake-lam', type=int, default=40,
                   help='wavelength blocks the raking matches exactly')
    g.add_argument('--rake-str', type=int, default=10,
                   help='strength bins the raking matches exactly')
    g.add_argument('--jackknife', type=int, default=40,
                   help='wavelength blocks left out one at a time for the '
                        'errors')
    g = p.add_argument_group('reporting')
    g.add_argument('--strong', type=float, default=100.0,
                   help='I_pred / noise for the epsilon measurement and the '
                        'branch-separation table')
    g.add_argument('--cmin', type=float, default=0.95,
                   help='minimum coverage for the epsilon measurement')
    g.add_argument('--branch-cuts', type=float, nargs='+',
                   default=[1.0, 10.0, 100.0],
                   help='strength cuts for the branch-separation table')
    g.add_argument('--pair-strength', type=float, default=5.0,
                   help='both members of a listed absent pair must be this '
                        'many times the noise level')
    g.add_argument('--pair-window', type=float, default=0.25,
                   help='largest separation, A, of a listed absent pair')
    g = p.add_argument_group('input and output')
    g.add_argument('--coverage',
                   default=os.path.join(ROOT, 'coverage_map.csv'))
    g.add_argument('--out-dir', default=ROOT)
    g.add_argument('--prefix', default='obscuration_length')
    g.add_argument('--no-plot', action='store_true')
    return p.parse_args(argv)


def main(argv=None):
    opt = parse_args(argv)
    lines = []

    def log(msg):
        print(msg)
        lines.append(str(msg))

    P, Y, offsets = outcomes(opt, log)

    strong = P[(P['c'] >= opt.cmin)
               & (P['I_pred'] >= opt.strong * P['thr'])]
    eps, elo, ehi, _ = orate.fit_eps(strong)
    if opt.epsilon is not None:
        log(f"epsilon on this run: {eps:.4f} ({elo:.4f} - {ehi:.4f}) from "
            f"{len(strong)} strong predictions; overridden by --epsilon "
            f"{opt.epsilon:g}")
        eps = opt.epsilon
    else:
        log(f"epsilon on this run: {eps:.4f} ({elo:.4f} - {ehi:.4f}) from "
            f"{len(strong)} strong predictions")

    lam = P['lam'].to_numpy()
    a = (~P['explained'].to_numpy()).astype(float)
    s = (P['I_pred'] / P['thr']).to_numpy()
    q = fit_baseline(P, opt, log)
    q = rake(q, a, lam, s, opt.rake_lam, opt.rake_str)
    qc = (1.0 - P['f_chance'].to_numpy())
    Ac = (~Y).astype(float)
    qc = qc * (Ac.mean() / qc.mean())

    I, J, D = build_pairs(lam, opt.dmax, log)
    AA = a[I] * a[J]
    QQ = q[I] * q[J]
    UU = (1.0 - q[I]) * (1.0 - q[J])
    AC = (Ac[:, I] * Ac[:, J]).mean(axis=0)
    QC = qc[I] * qc[J]

    cedges = np.array([0.0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.15, 0.4,
                       max(0.5, opt.dmax)])
    Rc, ctrl = chance_factor(D, AC, QC, cedges)
    log("chance factor from the controls: "
        + " ".join(f"{lo:g}-{hi:g}:{v:.3f}" for lo, hi, _, v in ctrl))

    bins = np.asarray(opt.bins, dtype=float)
    phi = phi_binned(D, AA, QQ, UU, Rc, bins, eps)
    blk = np.digitize(lam, np.quantile(lam,
                                       np.linspace(0, 1, opt.jackknife + 1)
                                       )[1:-1])
    pb = np.maximum(blk[I], blk[J])
    jk = np.array([phi_binned(D, AA, QQ, UU, Rc, bins, eps, sel=(pb != v))
                   for v in np.unique(pb)])
    nb = len(jk)
    phi_err = np.sqrt((nb - 1) / nb
                      * np.nansum((jk - np.nanmean(jk, axis=0)) ** 2, axis=0))

    grid_L = np.exp(np.linspace(np.log(0.004), np.log(1.0), 50))
    grid_f = np.linspace(-0.2, 1.5, 60)
    res = fit_length(D, AA, QQ, UU, Rc, eps, grid_L, grid_f, opt.dfit)
    L, Llo, Lhi, floor, sig = res
    # a stability check, not an error bar: the uncertainty on L is the
    # profile interval above, and a leave-one-block-out fit lands on the same
    # coarse grid of trial widths nearly every time
    jkL = np.array([fit_length(D, AA, QQ, UU, Rc, eps, grid_L, grid_f,
                               opt.dfit, sel=(pb != v))[0]
                    for v in np.unique(pb)])

    grid_e = np.exp(np.linspace(np.log(0.004), np.log(0.5), 40))
    ep, ep_lo, ep_hi, _ = fit_eps_from_clustering(
        D, AA, QQ, UU, Rc, grid_e, np.exp(np.linspace(np.log(0.004),
                                                      np.log(1.0), 40)),
        np.linspace(-0.02, 0.15, 25), opt.dfit)

    log("")
    log(f"L = {L:.4f} A   95% profile {Llo:.4f} - {Lhi:.4f}")
    log(f"  leaving out any one of the {nb} wavelength blocks: "
        f"{jkL.min():.4f} - {jkL.max():.4f} A")
    log(f"  against no correlation at all: 2 dNLL = {sig:.1f} on 1 degree of "
        f"freedom")
    log(f"  unremoved baseline correlation, fitted as a floor: "
        f"phi0 = {floor:+.3f}")
    fwhm = float(np.hypot(ls.INSTR_FWHM_A, ls.doppler_fwhm_factor() * 1500.0))
    log(f"  effective resolution at 1500 A, where nearly all the close pairs "
        f"lie: {fwhm:.4f} A")
    log(f"epsilon from the clustering amplitude alone: {ep:.4f} "
        f"({ep_lo:.4f} - {ep_hi:.4f}), against {eps:.4f} from the rate")
    log("")
    log(f"{'d_lo':>8}{'d_hi':>8}{'pairs':>8}{'phi':>9}{'error':>9}")
    k = np.searchsorted(bins, D, side='right') - 1
    bin_n = [int((k == b).sum()) for b in range(len(bins) - 1)]
    for b in range(len(bins) - 1):
        log(f"{bins[b]:8.3f}{bins[b+1]:8.3f}{bin_n[b]:8d}{phi[b]:9.3f}"
            f"{phi_err[b]:9.3f}")

    wl = within_level_separations(P, opt.branch_cuts, opt.cmin)
    log("")
    log("branches of one level, and how close they ever come:")
    for cut, npred, nsep, mn, med, seps in wl:
        log(f"  I >= {cut:6.0f} x noise: {npred:5d} predictions, {nsep:6d} "
            f"branch pairs, closest {mn:.3f} A, median {med:.0f} A, "
            f"{int((seps < L).sum())} closer than L")

    # the middle strength cut: low enough that such pairs exist at all, high
    # enough that an absence at that strength is worth explaining
    bp_cut = sorted(opt.branch_cuts)[len(opt.branch_cuts) // 2]
    bp = close_branch_pairs(P, bp_cut, opt.cmin, L)
    log(f"branches of one level closer than L, at {bp_cut:g} x the noise or "
        f"stronger: {len(bp)}"
        + ("" if not bp else "  ("
           + ", ".join(f"{l1:.3f}/{l2:.3f}"
                       + ("" if (e1 and e2) else " - one is ABSENT")
                       for _, _, l1, l2, _, _, _, e1, e2 in bp) + ")"))

    pairs = close_absent_pairs(P, I, J, D, opt.pair_window, opt.pair_strength)
    log(f"jointly absent pairs within {opt.pair_window:g} A, both at least "
        f"{opt.pair_strength:g} x the noise level: {len(pairs)}")

    os.makedirs(opt.out_dir, exist_ok=True)
    base = os.path.join(opt.out_dir, opt.prefix)
    pairs.to_csv(base + '.csv', index=False, lineterminator='\n')
    log(f"wrote {opt.prefix}.csv ({len(pairs)} pairs)")
    extra = dict(n_pred=len(P), n_pairs=len(I),
                 n_fit=int((D <= opt.dfit).sum()),
                 L_lo=jkL.min(), L_hi=jkL.max(), fwhm=fwhm,
                 eps_pair=ep, eps_lo=ep_lo, eps_hi=ep_hi, bin_n=bin_n,
                 branch_cut=bp_cut, branch_pairs=bp)
    write_report(base + '.txt', opt, eps, res, phi, phi_err, ctrl, wl, pairs,
                 extra, log)
    if not opt.no_plot:
        plot(base + '.png', opt, res, phi, phi_err, ctrl, wl, log)
    with open(base + '.log', 'w', newline='\n') as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == '__main__':
    main()
