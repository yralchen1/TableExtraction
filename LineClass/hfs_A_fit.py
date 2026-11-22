"""A constants from Sugar's resolved hfs components, anchored on the confirmed
semiempirical constants of the level compositions (Work_on_hfs_plan.md, Step
4, item 1).

WHAT THIS IS FOR
================
A level's magnetic dipole hyperfine constant A (cm^-1) sets how far its
hyperfine components F = |I - J| ... I + J lie from its center of gravity;
I = 5/2 is the nuclear spin of 141Pr and J the level's angular momentum.  The
pipeline moves a line by the A constants of its two levels only where A rests
on measurement, or on a calculation that measurements confirm (the user's
rule, plan section 3.1).  The levels whose A is "not determined" in
`A_hfs_levels.csv` - 39 of them on 2026-10-03, nearly all 4f5d2 - have their
lines left unmoved, with the unknown A as an allowance in the uncertainty.

Sugar printed, for many of his flagged (`*r`, `*v`) lines, the positions of the
resolved hyperfine components.  `hfs_patterns.py` showed that the tabulated
line is the strongest component and that the k-th listed companion is the k-th
rung of the ladder below it,

    d(k) = (k/2) [ (2 F1max + 1 - k) A1 - (2 F2max + 1 - k) A2 ],

with F1max = I + J1 and F2max = I + J2 for the lower (1) and upper (2) level.
Each listed position is one linear equation in the two A constants.

THE FIT
=======
Solving for every A at once (`hfs_patterns.fit_levels`) leaves a nearly
singular direction: adding one constant to every A changes the spacings only
through J1 - J2.  Here that direction is fixed by the levels whose
semiempirical A from the composition the spacings confirm (source
"composition ..."): each enters as a prior, A = A_SCALE * A_semi within its
uncertainty (`hfs_kappa.scaled_A`), so the fit can pull it only as far as the
spacings demand.  Every other level met by a pattern is free: the
undetermined ones, those whose A is a flag interval, the semiempirical
constants the spacings contradict (`CONTRADICTED`), and Reader and Sugar's
(`FREED_SOURCES`).

Reader and Sugar's 20 values of 4f2.6s were priors too until 2026-10-05.
Fitted with their own priors removed, the components put them at 0.76 to 1.10
times the published values, a scatter far beyond the uncertainties that no
single scale describes (the mean, 0.968 +- 0.005, is A_SCALE): a fault of
their treatment level by level, not of their radial parameters.  The user
had the component values adopted instead.

The weight of a measured position is Sugar's precision, 0.003 A in
wavelength, taken twice because a spacing is the difference of two measured
positions: sigma = sqrt(2) * 0.003e-8 * sigma_line^2 cm^-1.

A free level tied to no prior, even through other free levels, sits on a
singular direction of its own and is reported as not determined.

THE PRIORS COMBINED WITH THE SPACINGS (--combine-priors, 2026-10-09)
====================================================================
The fit's value for a level with a prior is the prior and the spacings
combined, and is better than either: on the 52 composition levels the
patterns reach, it brings the median u_A from 0.0062 to 0.0024 cm^-1.  The
table held the prior alone until the user had the combined values adopted
(2026-10-09); a leave-one-out test had shown the priors' u_A realistic (rms
z 0.82 against the spacings without them), so the combination is a fair one.

`--combine-priors` writes the fitted value and `fitted_u` into the level's
row of `A_hfs_levels.csv`, under a source that says so (COMBINED_PREFIX,
hfs_kappa.COMBINED_TAG), and moves the composition row it replaces into
`A_hfs_priors.csv` (PRIORS) the first time.  From then on the prior is read
from there (`anchors_of`): the fit must go on using the prior alone, or the
spacings would be counted twice, once in the table's value and once as
positions.  A second run recomputes every combined value from the priors and
the current patterns.  A combined row is on the measured scale - it is not
multiplied by A_SCALE again - and is no pure measurement: hfs_A_theory.py
leaves it out of the fit of its radial parameters.  A level whose prior the
spacings contradict can be left as it is (NOT_COMBINED) for the user to
decide; none is now.

CALCULATED PRIORS WHERE NO COMPOSITION PRIOR REACHES (2026-11-20)
=================================================================
A pattern between two levels that no composition prior reaches, even
through other patterns, leaves both on a singular direction: a pattern of a
transition between equal J measures only the difference of the two A, and
others hardly more.  Where such levels have a calculated A (hfs_A_theory),
those values become priors for that part of the fit only (`full_fit`), so
that the spacings can improve on them; the rest of the fit, and every level
a composition prior reaches, is unchanged.  The calculated constants are no
priors anywhere else (user, 2026-11-20: not a global use).

`--combine-priors` writes the combination for a level whose calculated
uncertainty it brings down to CALC_GAIN or less (it at least halves the
variance), whose prior it pulls by no more than AGREE_SIGMA, and for a
"not determined" or "flag interval" level of such a part, which has only
the spacings and its partners' priors; its source then begins
CALC_COMBINED_PREFIX (hfs_kappa.CALC_COMBINED_TAG), and the calculated row
it replaces is kept in `A_hfs_priors.csv`, as a composition prior is.  A
level whose combination gains less keeps its calculated row and still acts
as the prior of its partner.  First case, 24593.418 *r (000310 - 000655,
J 11/2 - 11/2): A(000310) -0.0051(109) calculated, -0.0107(42) combined.

THE CHECKS
==========
* per pattern: each free level's A from each of its patterns alone, the
  partner at its fitted value, against the joint value;
* flags: every flagged line with one accepted classification bounds A from
  one side (an `*r` line has D > 0, a `*v` line D < 0, with
  D = I (A2 J2 - A1 J1)); the bounds are compared with the fit;
* lines: the level's own lines in the classification table, of every class
  of character, give A through the kappa model, independently and weakly
  (plan section 3.1, route 3).  The table must have been made with the A
  table given by --levels: its O-C carry that table's constants.  The route
  is circular for a level whose semiempirical D the kappas were fitted with
  (000127).

Usage
    python hfs_A_fit.py                       # report only
    python hfs_A_fit.py --out hfs_A_anchored.csv
    python hfs_A_fit.py --write-levels        # adopt the values (see adoptable)
    python hfs_A_fit.py --combine-priors      # priors combined with the spacings

The classification table defaults to the hfs working set's,
`iter_hfs/line_classifications.csv`; it is read, never written.
"""
import argparse
import collections
import csv
import datetime
import math
import os
import re

import numpy as np

import hfs_A_theory
import hfs_components
import hfs_kappa
import hfs_patterns

HERE = os.path.dirname(os.path.abspath(__file__))
LINES = os.path.join('iter_hfs', 'line_classifications.csv')
A_LEVELS = 'A_hfs_levels.csv'
OUT = 'hfs_A_anchored.csv'
PRIORS = 'A_hfs_priors.csv'

#: the source of a prior combined with the spacings (--combine-priors)
COMBINED_PREFIX = 'resolved components %s, ' % hfs_kappa.COMBINED_TAG

#: the source of a calculated prior combined with the spacings, used only
#: where no composition prior reaches the pattern (full_fit)
CALC_COMBINED_PREFIX = 'resolved components %s, ' % \
    hfs_kappa.CALC_COMBINED_TAG

#: a calculated prior's combination is written only if it brings u_A down
#: to this fraction of the calculated one or less: at least half the
#: variance gone (user, 2026-11-20: calculated priors only where they
#: improve the value)
CALC_GAIN = 0.7

#: priors the spacings contradict, not combined (--combine-priors) but left
#: in the table as they are, with the reason.  Empty since 000229 was
#: combined (2026-10-09): its prior had been held because 30879.98 *r gave
#: +0.0568(43) against +0.0286(61), until the user dismissed that line's two
#: printed components as misidentified (hfs_components.DISMISSED); its other
#: pattern, 33021.15 *v, agrees with the prior.
NOT_COMBINED = {}

#: Sugar's precision on one tabulated wavelength, in cm (0.003 A).
DLAM = 0.003e-8

#: semiempirical constants the measured spacings contradict (plan section 3.1):
#: freed rather than anchored.
CONTRADICTED = {'059003.000127'}

#: sources of semiempirical constants that are not priors but free levels of
#: the fit, and replaced by `--write-levels` (user, 2026-10-05; see the
#: module docstring).
FREED_SOURCES = ('Reader & Sugar',)

#: patterns left out of the fit, by the line's wavenumber to 0.001 cm^-1.
SET_ASIDE = {
    29753.592: 'the one listed companion (+1.188) is rung 2 of the pattern '
               '(+1.06 to +1.20 predicted); rung 1 (+0.62) is not listed, so '
               'the ladder read as listed puts A(000274) 11 sigma off its '
               'other pattern',
    29151.128: 'the listed companion (-0.579) would make A(000151) = -0.068, '
               'a width no other line of this 4f2 5d level shows: its strong '
               '18418.856 is unflagged and narrow, and the calculated A is '
               '+0.023; the companion is taken to be another line (user, '
               '2026-10-07)',
}

#: printed components that are two or more rungs read as one, by the line's
#: wavenumber to 0.001 cm^-1: {wn: {k: (rungs, reason)}}, the k-th printed
#: component (in order of distance from the line) being the strength-weighted
#: mean of the listed rungs (their strengths depend on J alone, so the
#: equation stays linear in A).  Widths of 0.20-0.47 cm^-1 merge them.
MERGED_RUNGS = {
    33185.849: {1: ((1, 2), 'rung 1 (-0.270) is too close to the head to be '
                            'resolved and merges with rung 2 (-0.495); the '
                            'broadened maximum lies at -0.29 to -0.37 against '
                            'the printed -0.378 (Contradictions_for_paper.md '
                            'item 2f; user, 2026-11-21)')},
}


def rung_strength(k, J1, J2, spin=hfs_patterns.I_SPIN):
    """The relative strength of rung k (F1 = I + J1 - k -> F2 = I + J2 - k),
    in the units of hfs_patterns.components."""
    f1, f2 = spin + J1 - k, spin + J2 - k
    for _p, s, a, b in hfs_patterns.components(J1, 0.0, J2, 0.0, spin):
        if abs(a - f1) < 1e-9 and abs(b - f2) < 1e-9:
            return s
    return 0.0


def rung_coefficients(p, k, merged=None):
    """(c_low, c_upp): the k-th printed component of pattern `p` as
    c_low * A_low + c_upp * A_upp - rung k (hfs_patterns.rung), or the
    strength-weighted mean of the rungs MERGED_RUNGS names for it."""
    f1, f2 = hfs_patterns.I_SPIN + p.J1, hfs_patterns.I_SPIN + p.J2
    merged = MERGED_RUNGS if merged is None else merged
    rungs = merged.get(round(p.wn, 3), {}).get(k, ((k,), ''))[0]
    w = [rung_strength(j, p.J1, p.J2) for j in rungs] if len(rungs) > 1 \
        else [1.0]
    c1 = sum(x * 0.5 * j * (2 * f1 + 1 - j) for x, j in zip(w, rungs)) / sum(w)
    c2 = -sum(x * 0.5 * j * (2 * f2 + 1 - j) for x, j in zip(w, rungs)) / sum(w)
    return c1, c2

#: levels the fit determines but whose value is not adopted, with the reason.
HOLD = {
    '059003.000448': 'one component position, checked by nothing else; its '
                     'own lines give +0.013 +- 0.018',
}

#: the sources whose rows `--write-levels` may replace; its own earlier
#: values (SOURCE_PREFIX) are refreshed too, and FREED_SOURCES replaced.
REPLACEABLE = ('not determined', 'flag interval')
SOURCE_PREFIX = 'resolved components, '

#: a calculated A (`hfs_A_theory.py`) and a measured one replace each other
#: only when the newcomer is much more accurate - its uncertainty at most
#: MUCH_MORE_ACCURATE times the other's - and agrees with it within
#: AGREE_SIGMA combined standard uncertainties (user, 2026-10-07: neither
#: always wins).  Between the two ratios the row in the table stays.
MUCH_MORE_ACCURATE = 1.0 / 3.0
AGREE_SIGMA = 3.0

Fit = collections.namedtuple(
    'Fit', 'ids A u chi2 dof n_positions n_priors null pulls resid_rms')


def sigma_of(wn):
    """The uncertainty of one measured spacing on a line at wn, cm^-1."""
    return math.sqrt(2.0) * DLAM * wn * wn


def read_A_table(path=None):
    """`{level_id: row}` of `A_hfs_levels.csv`, in file order."""
    path = path or os.path.join(HERE, A_LEVELS)
    with open(path, encoding='utf-8', newline='') as fh:
        return collections.OrderedDict(
            (r['level_id'], r) for r in csv.DictReader(fh))


def read_priors(path=None):
    """`{level_id: row}` of `A_hfs_priors.csv` (PRIORS): the composition rows
    `--combine-priors` replaced in the table, as they were.  Empty if there
    is no such file."""
    path = path or os.path.join(HERE, PRIORS)
    if not os.path.exists(path):
        return collections.OrderedDict()
    with open(path, encoding='utf-8', newline='') as fh:
        return collections.OrderedDict(
            (r['level_id'], r) for r in csv.DictReader(fh))


def anchors_of(table, contradicted=CONTRADICTED, priors=None):
    """`{level_id: (A, u_A)}` of the confirmed semiempirical constants, on the
    measured scale: the priors of the fit.  A level whose table row is a
    prior combined with the spacings takes its prior from `priors`
    (`read_priors`), never the combined value; such a row without its prior
    stops the run."""
    priors = priors or {}
    out = {}
    for lid, r in table.items():
        src = r['source']
        if lid in priors and (hfs_kappa.is_combined(src)
                              or hfs_kappa.is_calc_combined(src)
                              or hfs_kappa.is_semiempirical(src)):
            r = priors[lid]
            src = r['source']
        elif hfs_kappa.is_combined(src) or hfs_kappa.is_calc_combined(src):
            raise SystemExit('hfs_A_fit.py: %s is a combined value (%s) but '
                             '%s has no prior for it' % (lid, src, PRIORS))
        if (lid in contradicted
                or not hfs_kappa.is_semiempirical(src)
                or src.startswith(FREED_SOURCES)):
            continue
        out[lid] = hfs_kappa.scaled_A(float(r['A_cm-1']), float(r['u_A']),
                                      src)
    return out


def calc_priors_of(table, priors=None):
    """`{level_id: (A, u_A)}` of the calculated constants: every row taken
    from hfs_A_theory (its source begins hfs_A_theory.CALC_SOURCE), and for a
    row that is already a calculated prior combined with the spacings, the
    prior kept in `priors` - never the combined value, which would count
    the spacings twice.  Such a row without its prior stops the run.  The
    calculated constants of a configuration hfs_kappa.CALC_CFG_SCALE names
    are taken on its scale (hfs_kappa.scaled_A)."""
    priors = priors or {}
    out = {}
    for lid, r in table.items():
        src = r['source']
        if hfs_kappa.is_calc_combined(src):
            if lid not in priors:
                raise SystemExit('hfs_A_fit.py: %s is a combined value (%s) '
                                 'but %s has no prior for it'
                                 % (lid, src, PRIORS))
            r = priors[lid]
        elif not src.startswith(hfs_A_theory.CALC_SOURCE):
            continue
        # on the pipeline's scale (hfs_kappa.CALC_CFG_SCALE, 2026-11-21);
        # the cfg is the table's, which a prior row may lack
        out[lid] = hfs_kappa.scaled_A(float(r['A_cm-1']), float(r['u_A']),
                                      r['source'],
                                      table[lid].get('cfg') or r.get('cfg'))
    return out


def full_fit(patterns, anchors, calc):
    """`(result, extra, dependent)`: the fit on the composition priors, and
    where it leaves levels undetermined, again with the calculated constants
    of those levels (`calc`, calc_priors_of) as priors too.

    `extra` holds the calculated priors used; `dependent` the levels
    determined only through them.  The levels the composition priors reach
    are untouched: an undetermined part of the fit is joined to the rest by
    no pattern, so adding priors to it moves nothing else."""
    res0 = fit(patterns, anchors)
    extra = {x: calc[x] for x in res0.null if x in calc}
    if not extra:
        return res0, {}, set()
    res = fit(patterns, dict(anchors, **extra))
    return res, extra, res0.null - res.null


def calc_combinable(lid, table, result, extra, dependent):
    """True if `--combine-priors` writes the level as the spacings combined
    with calculated priors: a level determined only through them
    (`dependent`) whose own calculated prior the combination brings down to
    CALC_GAIN or less of its u_A without pulling it beyond AGREE_SIGMA, or
    that was combined so before, or whose row is "not determined" or a flag
    interval (it has no value of its own to improve on)."""
    r = table.get(lid)
    if (r is None or lid not in dependent or lid in HOLD
            or lid in NOT_COMBINED):
        return False
    if lid in extra:
        return (hfs_kappa.is_calc_combined(r['source'])
                or (result.u[lid] <= CALC_GAIN * extra[lid][1]
                    and abs(result.pulls[lid]) <= AGREE_SIGMA))
    return r['source'] in REPLACEABLE


def calc_combined_source(n_patterns, n_positions, date):
    return (CALC_COMBINED_PREFIX + '%d pattern%s, %d position%s '
            '(Step 4 anchored fit, %s)'
            % (n_patterns, '' if n_patterns == 1 else 's', n_positions,
               '' if n_positions == 1 else 's', date))


def read_patterns(line_list, table, J_of, set_aside=SET_ASIDE):
    """The measured patterns on the current assignments.

    Returns (patterns, skipped).  `J_of` gives every level's J, so that a
    level absent from the A table can enter the fit too."""
    records, _ = hfs_components.build(line_list=line_list)
    rows = [{k: str(v) for k, v in r.items()} for r in records]
    levels = {}
    for lid, J in J_of.items():
        r = table.get(lid)
        levels[lid] = {'J': J,
                       'A_cm-1': r['A_cm-1'] if r else '0',
                       'u_A': r['u_A'] if r else '0',
                       'A_source': r['source'] if r else 'absent'}
    patterns, skipped = hfs_patterns.patterns_from(rows, levels,
                                                   'the level list')
    kept = []
    for p in patterns:
        if round(p.wn, 3) in set_aside:
            skipped['set aside by hand (SET_ASIDE)'] += 1
        else:
            kept.append(p)
    return kept, skipped


def _equations(patterns, anchors, ids):
    ix = {x: i for i, x in enumerate(ids)}
    rows, y, w = [], [], []
    for p in patterns:
        s = sigma_of(p.wn)
        for k, d in enumerate(p.dwn, 1):
            row = np.zeros(len(ids))
            c1, c2 = rung_coefficients(p, k)
            row[ix[p.low_id]] += c1
            row[ix[p.upp_id]] += c2
            rows.append(row)
            y.append(d)
            w.append(1.0 / s)
    n_pos = len(y)
    for x in ids:
        if x in anchors:
            row = np.zeros(len(ids))
            row[ix[x]] = 1.0
            rows.append(row)
            y.append(anchors[x][0])
            w.append(1.0 / anchors[x][1])
    w = np.array(w)
    return np.array(rows) * w[:, None], np.array(y) * w, n_pos


def fit(patterns, anchors):
    """The weighted least-squares solution with the anchors as priors.

    `u` is the formal uncertainty, not scaled by chi^2/dof (which comes out
    below 1 on Sugar's precision).  `null` is the set of levels on a singular
    direction, whose values mean nothing."""
    ids = sorted({p.low_id for p in patterns} | {p.upp_id for p in patterns})
    M, b, n_pos = _equations(patterns, anchors, ids)
    sol, *_ = np.linalg.lstsq(M, b, rcond=None)
    _, S, Vt = np.linalg.svd(M, full_matrices=False)
    null = set()
    for j in np.where(S < 1e-6 * S[0])[0]:
        null.update(ids[i] for i in np.where(np.abs(Vt[j]) > 1e-3)[0])
    cov = np.linalg.pinv(M.T @ M, rcond=1e-12)
    r = M @ sol - b
    pulls = {x: float(r[n_pos + i]) for i, x in
             enumerate(x for x in ids if x in anchors)}
    resid = [r[i] * sigma_of(p.wn) for i, p in
             enumerate(p for p in patterns for _ in p.dwn)]
    dof = len(b) - int(np.sum(S >= 1e-6 * S[0]))
    return Fit(ids, dict(zip(ids, sol)),
               dict(zip(ids, np.sqrt(np.clip(np.diag(cov), 0.0, None)))),
               float(r @ r), dof, n_pos, len(pulls), null, pulls,
               float(np.sqrt(np.mean(np.square(resid)))))


def single_pattern_values(patterns, result, lid):
    """[(wn, A, u)] - A of `lid` from each of its patterns alone, the partner
    held at its fitted value."""
    out = []
    for p in patterns:
        if lid not in (p.low_id, p.upp_id):
            continue

        def r(a):
            A1 = a if p.low_id == lid else result.A[p.low_id]
            A2 = a if p.upp_id == lid else result.A[p.upp_id]
            return np.array([sum(c * a for c, a in
                                 zip(rung_coefficients(p, k), (A1, A2))) - d
                             for k, d in enumerate(p.dwn, 1)])
        r0, g = r(0.0), r(1.0) - r(0.0)
        out.append((p.wn, -float(g @ r0) / float(g @ g),
                    sigma_of(p.wn) / math.sqrt(float(g @ g))))
    return out


def read_lines(path):
    """The accepted rows of a classification table."""
    with open(path, encoding='utf-8', newline='') as fh:
        return [r for r in csv.DictReader(fh)
                if r.get('accepted') in ('1', '1.0', 'True')]


def flag_bounds(lines, lid, A_of, J_of, known):
    """(lower, upper, n) bounds on A(lid) from its flagged lines with one
    accepted classification whose partner has a known A."""
    n_acc = collections.Counter(r['wn_obs'] for r in lines)
    lo, hi, n = -math.inf, math.inf, 0
    for r in lines:
        if r['char'] not in ('*r', '*v') or n_acc[r['wn_obs']] != 1:
            continue
        if lid not in (r['low_id'], r['upp_id']):
            continue
        other = r['upp_id'] if r['low_id'] == lid else r['low_id']
        if not known(other):
            continue
        a0 = A_of(other) * J_of[other] / J_of[lid]       # D = 0 here
        if (r['char'] == '*r') == (r['upp_id'] == lid):
            lo = max(lo, a0)
        else:
            hi = min(hi, a0)
        n += 1
    return lo, hi, n


def lines_value(lines, lid, A_now, J):
    """(A, u, n): A of `lid` from its own lines through the kappa model.

    A line of class kappa sits (1 - kappa) * BF * D below the head-frame
    Ritz value (BF its share of a blend), so with every other A held, its
    O-C moves by s (1 - kappa) BF I J dA, s = +1 if `lid` is the upper
    level, while a shift delta of the level's energy moves it by -s delta.
    Both are fitted, weighted by the line's uncertainty; A_now is the A the
    table was made with.  None for fewer than 3 lines."""
    X, y, w = [], [], []
    for r in lines:
        if lid not in (r['low_id'], r['upp_id']):
            continue
        s = 1.0 if r['upp_id'] == lid else -1.0
        f = ((1.0 - float(r['kappa'])) * float(r['BF'] or 1.0)
             * hfs_patterns.I_SPIN * J)
        oc0 = float(r['dif_wn_O-C']) - s * f * A_now
        X.append([s * f, -s])
        y.append(-oc0)
        w.append(1.0 / float(r['unc_wn_obs']) ** 2)
    if len(y) < 3:
        return None
    X, y, W = np.array(X), np.array(y), np.diag(w)
    try:
        cov = np.linalg.inv(X.T @ W @ X)
    except np.linalg.LinAlgError:
        return None
    p = cov @ X.T @ W @ y
    return float(p[0]), math.sqrt(cov[0, 0]), len(y)


def source_of(n_patterns, n_positions, date):
    return (SOURCE_PREFIX + '%d pattern%s, %d position%s '
            '(Step 4 anchored fit, %s)'
            % (n_patterns, '' if n_patterns == 1 else 's', n_positions,
               '' if n_positions == 1 else 's', date))


def much_more_accurate(A_new, u_new, A_old, u_old):
    """True if (A_new, u_new) may replace (A_old, u_old) in the table."""
    return (u_new <= MUCH_MORE_ACCURATE * u_old
            and abs(A_new - A_old) <= AGREE_SIGMA * math.hypot(u_new, u_old))


def fitted_u(result, lid):
    """The uncertainty `--write-levels` writes: the fitted one with the
    scale's added, since the values stand on the measured scale only
    through the anchors."""
    return math.hypot(result.u[lid], result.A[lid] * hfs_kappa.U_A_SCALE)


def adoptable(lid, table, result, dependent=()):
    """True if `--write-levels` replaces the row of `lid`.  A level
    determined only through calculated priors (`dependent`, full_fit) is not
    a measurement and is never written here; `--combine-priors` writes it
    (calc_combinable)."""
    r = table.get(lid)
    if (r is None or lid not in result.A or lid in result.null or lid in HOLD
            or lid in dependent):
        return False
    if (r['source'].startswith(hfs_A_theory.CALC_SOURCE)
            or hfs_kappa.is_calc_combined(r['source'])):
        return much_more_accurate(result.A[lid], fitted_u(result, lid),
                                  float(r['A_cm-1']), float(r['u_A']))
    return (r['source'] in REPLACEABLE or lid in CONTRADICTED
            or r['source'].startswith(SOURCE_PREFIX)
            or r['source'].startswith(FREED_SOURCES))


def combinable(lid, table, result, anchors):
    """True if `--combine-priors` writes the level's combined value: a
    composition prior (or an earlier combination of it) that the patterns
    reach and determine, not in NOT_COMBINED or HOLD."""
    r = table.get(lid)
    if (r is None or lid not in anchors or lid not in result.A
            or lid in result.null or lid in NOT_COMBINED or lid in HOLD):
        return False
    return (hfs_kappa.is_combined(r['source'])
            or r['source'].startswith('composition'))


def combined_source(n_patterns, n_positions, date):
    return (COMBINED_PREFIX + '%d pattern%s, %d position%s '
            '(Step 4 anchored fit, %s)'
            % (n_patterns, '' if n_patterns == 1 else 's', n_positions,
               '' if n_positions == 1 else 's', date))


def _keep_date(r, A, u, source_fn, n_patterns, n_positions, date):
    """The source to write for a row given the values `A`, `u` (as written):
    the row's own source if neither value changed and the source, with its own
    date, would read the same; else `source_fn` with today's `date`.  A rerun
    that changes nothing then leaves the table's dates alone (2026-11-21; it
    used to redate every combined row)."""
    if r.get('A_cm-1') == A and r.get('u_A') == u:
        m = re.search(r'(\d{4}-\d{2}-\d{2})\)\s*$', r.get('source', ''))
        if m and source_fn(n_patterns, n_positions, m.group(1)) == r['source']:
            return r['source']
    return source_fn(n_patterns, n_positions, date)


def write_combined(table, priors, result, anchors, counts, path, priors_path,
                   date, extra=None, dependent=()):
    """`--combine-priors`: write the combined values into `path` and keep the
    rows they replace in `priors_path` - the composition priors, and the
    calculated ones of `extra` (calc_combinable).  Returns the level ids
    written.  The priors file is written first, so that no prior is ever
    lost."""
    extra = extra or {}
    with open(path, encoding='utf-8', newline='') as fh:
        fields = csv.DictReader(fh).fieldnames
    changed, calc = [], set()
    for lid, r in table.items():
        if calc_combinable(lid, table, result, extra, dependent):
            calc.add(lid)
        elif not combinable(lid, table, result, anchors):
            continue
        # a "not determined" row combined through its partners' priors has
        # no prior of its own to keep
        if lid not in priors and (lid not in calc or lid in extra):
            priors[lid] = dict(r)
            priors[lid]['moved'] = date
        changed.append(lid)
    pfields = list(fields) + ['moved']
    with open(priors_path, 'w', encoding='utf-8', newline='\n') as fh:
        w = csv.DictWriter(fh, fieldnames=pfields, lineterminator='\n')
        w.writeheader()
        w.writerows(priors.values())
    for lid in changed:
        r = table[lid]
        A = '%+.4f' % result.A[lid]
        if lid in calc:
            # on the measured scale through the calculated priors, not
            # through the composition scale
            u, fn = '%.4f' % result.u[lid], calc_combined_source
        else:
            u, fn = '%.4f' % fitted_u(result, lid), combined_source
        r['source'] = _keep_date(r, A, u, fn, counts[lid][0], counts[lid][1],
                                 date)
        r['A_cm-1'], r['u_A'] = A, u
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator='\n')
        w.writeheader()
        w.writerows(table.values())
    return changed


def write_levels(table, result, counts, path, date, dependent=()):
    """Rewrite `A_hfs_levels.csv` with the adoptable fitted values.

    The value is the fitted A; its uncertainty `fitted_u`.  Every other
    row is written back unchanged.  Returns the list of level ids
    replaced."""
    with open(path, encoding='utf-8', newline='') as fh:
        fields = csv.DictReader(fh).fieldnames
    changed = []
    for lid, r in table.items():
        if not adoptable(lid, table, result, dependent):
            continue
        A, u = '%+.4f' % result.A[lid], '%.4f' % fitted_u(result, lid)
        r['source'] = _keep_date(r, A, u, source_of, counts[lid][0],
                                 counts[lid][1], date)
        r['A_cm-1'], r['u_A'] = A, u
        changed.append(lid)
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator='\n')
        w.writeheader()
        w.writerows(table.values())
    return changed


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--lines', default=os.path.join(HERE, LINES),
                    help='classification table (default iter_hfs)')
    ap.add_argument('--levels', default=os.path.join(HERE, A_LEVELS))
    ap.add_argument('--out', default=None,
                    help='write every fitted level to this csv')
    ap.add_argument('--write-levels', action='store_true',
                    help='adopt the values into the A table')
    ap.add_argument('--priors', default=None,
                    help='the composition priors replaced in the table '
                         '(default: %s beside --levels)' % PRIORS)
    ap.add_argument('--combine-priors', action='store_true',
                    help='write each composition prior the patterns reach, '
                         'combined with the spacings, into the A table, and '
                         'keep the prior in --priors')
    args = ap.parse_args(argv)

    table = read_A_table(args.levels)
    priors_path = args.priors or os.path.join(
        os.path.dirname(os.path.abspath(args.levels)), PRIORS)
    priors = read_priors(priors_path)
    J_of = hfs_kappa.read_levels()
    anchors = anchors_of(table, priors=priors)
    patterns, skipped = read_patterns(args.lines, table, J_of)
    print('patterns used %d' % len(patterns))
    for k, v in sorted(skipped.items()):
        print('   set aside, %-44s %3d' % (k + ':', v))
    for wn, why in SET_ASIDE.items():
        print('      %.3f: %s' % (wn, why))

    res, extra, dependent = full_fit(patterns, anchors,
                                     calc_priors_of(table, priors))
    print('\n%d positions, %d priors (%d of them calculated), %d levels; '
          'chi2/dof %.3f; rms of the positions %.4f cm^-1'
          % (res.n_positions, res.n_priors, len(extra), len(res.ids),
             res.chi2 / res.dof, res.resid_rms))
    pulled = sorted((abs(v), x) for x, v in res.pulls.items()
                    if abs(v) > 2 and x in anchors)
    for _, x in reversed(pulled):
        print('   prior pulled %+.1f sigma: %s, %+.4f -> %+.4f'
              % (res.pulls[x], x, anchors[x][0], res.A[x]))
    if res.null:
        print('   not determined (no prior reaches them): %s'
              % ', '.join(sorted(res.null)))

    counts = {}
    for p in patterns:
        for x in (p.low_id, p.upp_id):
            n, m = counts.get(x, (0, 0))
            counts[x] = (n + 1, m + len(p.dwn))

    lines = read_lines(args.lines)
    now = {k: v[1] for k, v in hfs_kappa.read_A_constants(args.levels).items()}

    def known(x):
        return ((x in table and table[x]['source'] != 'not determined')
                or (x in res.A and x not in res.null))

    def A_of(x):
        return res.A[x] if x in res.A and x not in res.null else now.get(x, 0.0)

    free = [x for x in res.ids if x not in anchors]
    print('\n%-14s %-6s %4s %3s %3s %9s %7s  %9s %7s  %-17s %-17s %-17s %s'
          % ('level', 'cfg', 'J', 'pat', 'pos', 'A_fit', 'u', 'A_table', 'u',
             'flag bounds', 'own lines', 'worst pattern', 'table source'))
    rows_out = []
    for x in free:
        r = table.get(x)
        lo, hi, nf = flag_bounds(lines, x, A_of, J_of, known)
        own = lines_value(lines, x, now.get(x, 0.0), J_of[x])
        single = single_pattern_values(patterns, res, x)
        worst = max(((a - res.A[x]) / ua for _, a, ua in single),
                    key=abs) if len(single) > 1 else None
        z_own = ((own[0] - res.A[x]) / math.hypot(own[1], res.u[x])
                 if own else None)
        status = ('not determined' if x in res.null else
                  'held' if x in HOLD else
                  'adopt' if adoptable(x, table, res, dependent) else
                  'calc prior' if x in dependent else '')
        print('%-14s %-6s %4.1f %3d %3d %+9.4f %7.4f  %9s %7s  %-17s %-17s '
              '%-17s %s%s'
              % (x, r['cfg'] if r else '', J_of[x], counts[x][0],
                 counts[x][1], res.A[x], res.u[x],
                 r['A_cm-1'] if r else '', r['u_A'] if r else '',
                 '%s..%s' % ('%+.3f' % lo if lo > -math.inf else '-inf',
                             '%+.3f' % hi if hi < math.inf else '+inf')
                 if nf else '',
                 '%+.4f(%.4f) %+.1f' % (own[0], own[1], z_own) if own else '',
                 '%+.1f sigma' % worst if worst is not None else '',
                 (r['source'][:30] if r else 'absent from the table'),
                 '   [%s]' % status if status else ''))
        rows_out.append([x, r['cfg'] if r else '', J_of[x], counts[x][0],
                         counts[x][1], '%+.4f' % res.A[x],
                         '%.4f' % res.u[x],
                         r['A_cm-1'] if r else '', r['u_A'] if r else '',
                         r['source'] if r else 'absent',
                         '%+.4f' % lo if nf and lo > -math.inf else '',
                         '%+.4f' % hi if nf and hi < math.inf else '',
                         '%+.4f' % own[0] if own else '',
                         '%.4f' % own[1] if own else '',
                         '%+.2f' % worst if worst is not None else '',
                         status])
    for x, why in HOLD.items():
        print('held: %s - %s' % (x, why))

    if args.out:
        with open(args.out, 'w', encoding='utf-8', newline='\n') as fh:
            w = csv.writer(fh, lineterminator='\n')
            w.writerow(['level_id', 'cfg', 'J', 'n_patterns', 'n_positions',
                        'A_fit', 'u_fit', 'A_table', 'u_table',
                        'source_table', 'flag_lower', 'flag_upper',
                        'A_lines', 'u_lines', 'worst_pattern_sigma',
                        'status'])
            w.writerows(rows_out)
        print('\nwrote %s (%d levels)' % (args.out, len(rows_out)))
    if args.write_levels:
        date = datetime.date.today().isoformat()
        changed = write_levels(table, res, counts, args.levels, date,
                               dependent)
        print('\nwrote %s: %d levels replaced' % (args.levels, len(changed)))
    combine = [x for x in res.ids if combinable(x, table, res, anchors)]
    print('\npriors the spacings tighten: %d' % len(combine))
    for x in combine:
        p = anchors[x]
        print('   %s J %.1f  prior %+.4f(%.4f) -> combined %+.4f(%.4f)  '
              'u_S %.4f -> %.4f'
              % (x, J_of[x], p[0], p[1], res.A[x], fitted_u(res, x),
                 2.5 * J_of[x] * p[1], 2.5 * J_of[x] * fitted_u(res, x)))
    for x, why in NOT_COMBINED.items():
        print('not combined: %s - %s' % (x, why))
    print('\nlevels no composition prior reaches, determined through '
          'calculated priors: %d (written if u_A falls to %.2f of the '
          'calculated one or less)' % (len(dependent), CALC_GAIN))
    for x in sorted(dependent):
        r = table.get(x)
        p = extra.get(x)
        print('   %s J %.1f  %s -> combined %+.4f(%.4f)%s  [%s]'
              % (x, J_of[x],
                 'calculated %+.4f(%.4f), pulled %+.1f sigma'
                 % (p[0], p[1], res.pulls[x]) if p else
                 'no calculated A (%s)' % (r['source'] if r else 'absent'),
                 res.A[x], res.u[x],
                 '  u ratio %.2f' % (res.u[x] / p[1]) if p else '',
                 'write' if calc_combinable(x, table, res, extra, dependent)
                 else 'keep'))
    if args.combine_priors:
        import output_files
        output_files.require_writable([args.levels, priors_path])
        date = datetime.date.today().isoformat()
        changed = write_combined(table, priors, res, anchors, counts,
                                 args.levels, priors_path, date, extra,
                                 dependent)
        print('\nwrote %s: %d priors combined with the spacings; the priors '
              'are kept in %s' % (args.levels, len(changed), priors_path))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
