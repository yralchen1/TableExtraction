"""Detection of interchanged level identifications.

The failure this looks for
--------------------------
Every energy level of Pr III carries two things that are established
separately:

  * its ENERGY, a number in cm^-1 measured from the observed lines.  A level
    is found by noticing that several observed lines have wavenumber
    differences that fit one common position, and the least-squares
    optimization then fixes that position to a hundredth of a cm^-1.

  * its THEORETICAL IDENTITY: which level of the parametric (Cowan-code)
    calculation it is.  That identity is a label - a configuration such as
    4f^2.5d and a term such as ^3H - and it carries with it the whole set of
    CALCULATED transition probabilities gA of the level, i.e. the predicted
    intensity of every line the level can emit.  In this project those
    predictions live in Icalc.xlsx, keyed by the level id, and the label
    itself lives in the IDEN2 working file enlev.dat.

The identity is assigned by matching the observed energy to a calculated
energy, and the calculated energies of this ion are not accurate.  How
inaccurate depends strongly on the configuration: over the levels whose
energies are experimentally known, the root-mean-square of E_obs - E_calc is
22 cm^-1 for 4f^2.5g but 404 cm^-1 for 4f.5d.6p (the whole table is printed
by this tool).  Wherever two calculated levels of the SAME parity and the
SAME total angular momentum J lie closer together than that error, the energy
match cannot tell which is which, and the two identities can have been
INTERCHANGED: level A wears B's label and B's gA values, and level B wears
A's.

Why it matters, and why nothing else here catches it
----------------------------------------------------
The energies stay right under an interchange - they were measured from the
lines and no line moves - so the energy-shift test of level_shifts.py (dE,
the distance the optimizer drags a level from its input energy) says nothing
at all.  What goes wrong is the INTENSITIES: each level is then judged
against the predicted branching pattern of the other level, so its observed
lines and its predicted ones stop agreeing.  level_shifts.py sees this only
through pattern_V, one ingredient of a folded probability, and a level can
keep a small p_spur with a bad pattern.  In the run of 2026-09-05 the worst
pattern of all 208 tested levels belonged to 059003.000572, whose
configuration 4f.5d.6p is exactly the one with the largest calculated-energy
error - the signature of this failure - and it carried no questionable mark.

The test
--------
Take a candidate pair (A, B): same parity, same J, and close enough in energy
that the calculation cannot separate them (see "The search window" below).
Neither level moves.  What is exchanged is only which set of gA values
belongs to which level.  For each of the two hypotheses

    ASSIGNED  - A keeps A's gA values, B keeps B's;
    SWAPPED   - A is scored against B's gA values, B against A's,

every accepted line of A and of B is given the intensity that hypothesis
predicts for it,

    I_calc  =  C * gA(identity, partner) * (wn/1e8) * exp(-E_up/kT) ,

with the observed wavenumber wn of the line and the observed energy E_up of
its upper level (both unchanged by the swap - only gA changes), and the
constants C and kT of the pipeline's Boltzmann intensity model.  A level pair
absent from Icalc.xlsx is not a transition of unknown strength but one known
to be weaker than the printing cutoff of Cowan's codes, and gets the imputed
gA of gA_imputation.py, exactly as in the classification itself; so a
hypothesis that predicts nothing where a strong line is observed is penalized
rather than excused.

Two numbers then compare the hypotheses.

1.  THE SLOPE.  For each level, ln(I_obs) is fitted against ln(I_calc) by
    ordinary least squares.  A correct identity gives a slope near 1: the
    observed intensities follow the predicted ones over their whole range.  A
    wrong identity gives a slope near 0 - the predicted pattern carries no
    information about the observed one - and often a negative one.  The tool
    reports the four slopes (each level under each hypothesis) and prefers
    the interpretation whose slopes are consistently nearer 1.

2.  THE LIKELIHOOD RATIO, which is what actually decides.  The slope alone
    ignores the SCALE: a set of predictions ten times too large can still
    have slope 1.  Under a correct identity the residual

        r = ln(I_obs) - ln(I_calc * f(lambda))

    is centred on zero with a known spread - f and the spread are the band
    scale factor and the standard deviation of ln(I_obs/I_pred) measured on
    the run itself by level_shifts.intensity_scale_bias, which is how the
    far-ultraviolet scale error of the calculated intensities is taken out.
    Summing -r^2/(2*sd^2) over the accepted lines of both levels under each
    hypothesis and subtracting gives

        lnR = ln L(assigned) - ln L(swapped) ,

    positive when the assignment as it stands fits the observed intensities
    better, negative when the swap does.  lnR is a difference of two fits to
    the SAME observed intensities, so everything common to the two - the
    plate calibration, the Boltzmann factor, the observed intensity scale -
    cancels.

The search window
-----------------
Two levels are candidates for an interchange when the calculation cannot
place them in order, i.e. when their separation is no larger than the error
of the calculated energies THERE.  That error is not a global number: it is
taken per configuration, as the root-mean-square of E_obs - E_calc over the
experimentally known levels of the level's own dominant configuration, read
from enlev.dat.  The window of a pair is `--window` times the larger of the
two levels' configuration rms values (default 1.0).

The report also gives, for each pair, how far the swap would move each level
from its calculated position, in units of that same rms (the columns
omc_swap_A_sigma and omc_swap_B_sigma).  A swap that puts both levels within
about one rms of their new calculated positions is energetically as good as
the assignment; one that pushes a level to three rms is not, whatever the
intensities say.

How large an lnR is large
-------------------------
lnR is calibrated against swaps that are known to be wrong.  For every level
pair of the same J and parity that is NOT a candidate - separated by more
than the window but by less than `--null-window` times the configuration rms
(default 10), so that the two levels are still in the same part of the
spectrum and their gA values are of comparable magnitude - the same lnR is
computed.  Every one of those swaps is false by construction, so their lnR
distribution is the null: it shows how strongly the evidence normally favours
the true assignment.  Each candidate pair gets p_null, the fraction of null
pairs whose lnR is as low as or lower than its own.  A small p_null means the
swap is far better supported than a swap of this kind ever is.

This is the same device the rest of the validation uses: decoy_mc.py plants
levels that cannot be real and measures what the algorithm does with them;
here the null swaps cannot be real and measure what the intensities do with
them.

Reading the output
------------------
One row per candidate pair in level_interchange.csv (+ .xlsx twin), sorted by
lnR, most swap-favouring first.  A pair is FLAGGED when all three hold: lnR
is negative (the swap fits the intensities better), both slopes are nearer 1
under the swap than as assigned, and p_null is at or below `--p-flag`
(default 0.05).  Applying that same rule to the calibration swaps, all of
which are false, gives the expected number of false flags among the
candidates, which the report prints - the analogue of the expected false
passes in the criterion grids of level_shifts.py.

Read lnR_A and lnR_B, not only lnR: a genuine interchange improves both sides.
A large gain on one level alone says that level's label is wrong while the
other's true identity may be a third level, very possibly one not yet
observed, which this test cannot reach.

A flag is a request for a look in IDEN2, not a verdict.  Relabelling changes
no energy and no identification - only which theoretical level the published
table names - so it is the analyst's decision, taken with the term structure
and the g-factors in hand.

Usage
-----
    python level_interchange.py
    python level_interchange.py --lopt LOPT_output_lines_revised.txt \
                                --e-input revised_level_energies.csv
    python level_interchange.py --detail 059003.000572 059003.000560
        the line-by-line table of one pair: every accepted line of either
        level, its observed intensity, and the intensity each hypothesis
        predicts for it.
    python level_interchange.py --window 1.5 --min-lines 2
"""
import argparse
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import chance_mc as mc
import output_files
import classify_lines as cl
import gA_imputation
import level_shifts as ls

ENLEV = os.path.join(HERE, 'IDEN2', 'enlev.dat')
REPORT_CSV = os.path.join(HERE, 'level_interchange.csv')

# A level of the run is tied to its row in enlev.dat by energy: the two lists
# come from the same identification work and agree far more closely than this,
# so no match within it means the files are out of step.
ENLEV_MATCH_TOL = 0.5   # cm^-1

DECIMALS = {'E_A': 4, 'E_B': 4, 'dE_AB': 3, 'window': 1,
            'Ecalc_A': 3, 'Ecalc_B': 3,
            'omc_swap_A_sigma': 2, 'omc_swap_B_sigma': 2,
            'slope_A_assigned': 3, 'slope_A_swapped': 3,
            'slope_B_assigned': 3, 'slope_B_swapped': 3,
            'rms_assigned': 3, 'rms_swapped': 3,
            'lnR': 2, 'lnR_A': 2, 'lnR_B': 2, 'p_null': 4}


# ---------------------------------------------------------------------------
# The theoretical identities: enlev.dat
# ---------------------------------------------------------------------------
def split_label(label: str):
    """("configuration", "separator + term") of an enlev.dat label.

    The separator is "_" or "~" and never occurs inside either part, so the
    first one found divides them; the term keeps its separator because that is
    how IDEN2 itself writes and sorts these labels.  A label with no separator
    is all configuration.
    """
    for i, ch in enumerate(label):
        if ch in '_~':
            return label[:i].strip(), label[i:].strip()
    return label.strip(), ''


_ORBITAL_L = dict(s=0, p=1, d=2, f=3, g=4, h=5, i=6)


def configuration_parity(cfg: str) -> str:
    """'e' or 'o' for an enlev.dat configuration abbreviation.

    IDEN2 writes a configuration as its open shells run together, each an
    orbital letter optionally preceded by the principal quantum number and
    followed by its occupation: "f26p" is 4f^2.6p, "fd6p" is 4f.5d.6p, "5d3"
    is 5d^3.  The parity is the parity of the sum of l over the electrons, so
    it needs only the letters and the occupations, and the principal quantum
    numbers - which are digits too - must not be counted.  A digit is an
    occupation exactly when it FOLLOWS its orbital letter; a digit before a
    letter is that letter's principal quantum number.
    """
    total = 0
    i = 0
    while i < len(cfg):
        ch = cfg[i].lower()
        if ch in _ORBITAL_L:
            n = 1
            if i + 1 < len(cfg) and cfg[i + 1].isdigit():
                n = int(cfg[i + 1])
                i += 1
            total += _ORBITAL_L[ch] * n
        i += 1
    return 'o' if total % 2 else 'e'


def read_enlev(path: str = ENLEV) -> pd.DataFrame:
    """The theoretical level list of IDEN2, with its calculated energies.

    enlev.dat has one fixed-width row per level of the calculation:

        columns   1-5   running index
                  6-16  E_calc, the calculated energy (cm^-1)
                 17-27  the uncertainty of the observed energy
                 28-39  E_obs, the observed energy (a copy of E_calc when the
                        level has not been found)
                 39-40  ' *' when the level has been found experimentally
                 42-53  E_obs - E_calc
                 54-    J, then the label between slashes, e.g. "f25g _3H4G"
                        or "p5f3d_4I6Ga": the configuration in IDEN2's own
                        abbreviation (f25g = 4f^2.5g, fd6p = 4f.5d.6p), then a
                        separator - "_" for LS coupling, "~" for jj - then the
                        term.  The field is padded to a fixed width but not
                        always by the same amount, so the configuration is
                        taken as everything before the first separator, not by
                        a column cut.

    Returns the columns idx, E_calc, u_obs, E_obs, omc, J, label, cfg, term,
    known.

    `known` - has this level been found experimentally? - is read from the STAR
    ALONE.  That is the flag IDEN2 itself sets and reads, and the flag its own
    tools write when a level is fixed.  The uncertainty column is not a second
    opinion on it: 5000 cm^-1 there is the placeholder the conversion from
    Cowan's output writes for every level, and IDEN2 changes it only when the
    user orders that level by level, so a level that has been found may still
    be carrying the placeholder.  Reading the uncertainty instead of the star
    silently loses such a level (059003.000623 in the run of 2026-09-05).
    """
    recs = []
    with open(path, encoding='latin-1') as fh:
        for raw in fh:
            line = raw.rstrip('\r\n')
            if len(line) < 54 or '/' not in line:
                continue
            try:
                idx = int(line[0:5])
                e_calc = float(line[5:16])
                u_obs = float(line[16:27])
                e_obs = float(line[27:39])
                omc = float(line[41:53])
            except ValueError:
                continue
            rest = line[53:]
            j = float(rest.split()[0])
            label = rest.split('/')[1].strip()
            cfg, term = split_label(label)
            recs.append((idx, e_calc, u_obs, e_obs, omc, j, label,
                         cfg, term, '*' in line[38:40]))
    df = pd.DataFrame(recs, columns=['idx', 'E_calc', 'u_obs', 'E_obs', 'omc',
                                     'J', 'label', 'cfg', 'term', 'known'])
    if df.empty:
        raise ValueError(f"{path}: no level rows could be read")
    return df


def configuration_windows(en: pd.DataFrame) -> pd.Series:
    """{configuration: rms of E_obs - E_calc over its known levels}.

    The measure of how far the calculation can be wrong in that part of the
    spectrum, and therefore of how far apart two levels may lie and still be
    interchangeable.  Only experimentally known levels enter: for the others
    E_obs is a copy of E_calc and the difference is zero by construction.
    """
    known = en[en['known']]
    return known.groupby('cfg')['omc'].apply(
        lambda s: float(np.sqrt(float((s ** 2).mean()))))


def attach_identities(levels: pd.DataFrame, en: pd.DataFrame,
                      windows: pd.Series, tol: float = ENLEV_MATCH_TOL):
    """Give every level of the run its configuration, E_calc and window.

    The two lists share no key, so they are tied by energy: enlev.dat holds
    the observed energies of the same identification work.  Only the rows
    marked found - the starred ones - carry a real observed energy; in the rest
    E_obs is a copy of E_calc, and matching against those would tie a level of
    the run to a position nobody has measured.  Returns the levels with the
    columns cfg, term, E_calc, omc and W added, and the list of level ids that
    no row of enlev.dat matches.
    """
    known = en[en['known']].sort_values('E_obs').reset_index(drop=True)
    e_arr = known['E_obs'].to_numpy()
    cfg, term, e_calc, omc, unmatched = [], [], [], [], []
    for lid, e in zip(levels['level_id'], levels['E_final']):
        i = int(np.argmin(np.abs(e_arr - float(e))))
        if abs(e_arr[i] - float(e)) > tol:
            unmatched.append(lid)
            cfg.append('')
            term.append('')
            e_calc.append(float('nan'))
            omc.append(float('nan'))
            continue
        row = known.iloc[i]
        cfg.append(row['cfg'])
        term.append(row['term'])
        e_calc.append(float(row['E_calc']))
        omc.append(float(row['omc']))
    out = levels.copy()
    out['cfg'] = cfg
    out['term'] = term
    out['E_calc'] = e_calc
    out['omc'] = omc
    out['W'] = out['cfg'].map(windows)
    return out, unmatched


def nearest_enlev_row(en: pd.DataFrame, e: float):
    """(row, |difference|) of the enlev.dat row nearest in E_obs to `e`.

    Starred or not - which is the point.  A level of the run that matched
    nothing is reported with this row beside it, because the two files coming
    out of step is a far likelier cause than a level the calculation does not
    have at all: a row sitting at the right energy but not starred means IDEN2
    has not been told the level was found.
    """
    d = (en['E_obs'] - float(e)).abs().to_numpy()
    i = int(np.argmin(d))
    return en.iloc[i], float(d[i])


# ---------------------------------------------------------------------------
# The calculated transition strengths
# ---------------------------------------------------------------------------
def read_gA() -> dict:
    """{(lower id, upper id) sorted: gA} from Icalc.xlsx and icalc_new.xlsx.

    gA (statistical weight times transition probability, s^-1) is the quantity
    the swap exchanges: it belongs to the theoretical level, not to the
    observed one.  The supplementary file is read after the main one for the
    same reason as in level_shifts.load_predictions - a level found since the
    calculated table was made has its transitions only there.
    """
    import openpyxl
    out = {}
    files = [cl.ICALC_FILE]
    if getattr(cl, 'ICALC_EXTRA', ''):
        files.append(cl.ICALC_EXTRA)
    for path in files:
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        ws = wb[cl.CFG.icalc.sheet]
        col = cl.column_index(ws, cl.CFG.icalc, path)
        for row in ws.iter_rows(min_row=2):
            a = cl.to_str_id(row[col['id1']].value)
            b = cl.to_str_id(row[col['id2']].value)
            g = row[col['gA']].value
            if not a or not b or g is None:
                continue
            try:
                g = float(g)
            except (TypeError, ValueError):
                continue
            if g > 0:
                out[tuple(sorted((a, b)))] = g
        wb.close()
    return out


def imputed_gA() -> float:
    """The gA standing for a transition absent from the calculated file.

    Absence means gA below the printing cutoff of Cowan's codes, not an
    unknown strength; gA_imputation.py derives the value whose upper
    one-standard-deviation bound sits on the cutoff.  Returns 0.0 under
    missing_gA.policy = "none", in which case such lines are dropped from both
    hypotheses instead (see common_lines).
    """
    if cl.MISSING_POLICY != 'impute':
        return 0.0
    mg = cl.CFG.missing_gA
    df = gA_imputation.load_icalc(cl.CFG)
    g, _u = gA_imputation.estimate_missing_gA(
        df, cl.CFG.gA_cutoff, mg['u_ln_window'], mg['u_ln_estimator'],
        mg['self_consistent'], mg['fit_range_decades'])
    return float(g)


# ---------------------------------------------------------------------------
# The observed side
# ---------------------------------------------------------------------------
def accepted_lines_by_level(real: pd.DataFrame, e_final: dict) -> dict:
    """{level_id: [(partner id, wn_obs, I_obs, E_up, n_blend, sib), ...]}.

    One entry per accepted assignment of one level, holding what the intensity
    comparison needs about it.  E_up is the observed energy of the upper of the
    two levels, which the intensity model needs and which an interchange does
    not touch.

    The last two entries describe the BLEND the assignment belongs to.  One
    measured line is often produced by several predicted transitions falling at
    the same wavenumber; the classification then assigns all of them to it, and
    every one of those rows carries the SAME measured intensity, that of the
    whole feature.  Crediting each of them with the whole of it would say that
    every component is as strong as the feature - which, for a component the
    calculation makes a twentieth of it, overstates the measurement twentyfold.

    `n_blend` is how many accepted transitions share the measured line (1 when
    it is not blended), and `sib` is the sum of the CALCULATED intensities of
    the OTHER components, or None when one of them has no calculated intensity
    at all.  From those two, level_residuals works out the share of the measured
    intensity the component being scored should be given - separately under each
    hypothesis, because a hypothesis that changes a component's calculated
    intensity changes its share of the blend with it.  See blend_share.
    """
    acc = real[real['accepted'] == 1]
    # The components of one blend are the accepted rows carrying the same
    # printed wavenumber; that grouping reproduces the n_accepted and BF
    # columns classify_lines.py writes, exactly.
    key = acc['wn_obs'].astype(str)
    ci = pd.to_numeric(acc['calc_intens'], errors='coerce')
    grp = ci.groupby(key)
    n_bl = key.map(key.value_counts())
    tot = grp.transform('sum')
    bad = grp.transform(lambda c: (~(c > 0)).any())
    out = {}
    for lo, up, wn, i_obs, own, nb, tt, ms in zip(
            acc['low_id'], acc['upp_id'], acc['wn_obs'], acc['obs_intens'],
            ci, n_bl, tot, bad):
        lo, up = str(lo), str(up)
        if lo not in e_final or up not in e_final:
            continue
        try:
            wn = float(wn)
            i_obs = float(i_obs)
        except (TypeError, ValueError):
            continue
        if not (wn > 0 and i_obs > 0):
            continue
        nb = int(nb)
        sib = None if (nb > 1 and bool(ms)) else float(tt) - float(own)
        e_up = max(e_final[lo], e_final[up])
        out.setdefault(lo, []).append((up, wn, i_obs, e_up, nb, sib))
        out.setdefault(up, []).append((lo, wn, i_obs, e_up, nb, sib))
    return out


# ---------------------------------------------------------------------------
# The comparison
# ---------------------------------------------------------------------------
def predicted_intensity(g: float, wn: float, e_up: float,
                        C: float, kT: float) -> float:
    """The Boltzmann intensity model of the pipeline, on the observed scale.

    I_calc = C * gA * (wn/1e8) * exp(-E_up/kT); the factor wn/1e8 is the
    reciprocal of the vacuum wavelength in angstroms, and makes the predicted
    intensity an energy flux rather than a photon number.
    """
    return C * g * (wn / 1.0e8) * math.exp(-e_up / kT)


def blend_share(i_calc: float, n_blend: int, sib) -> float:
    """The share of a blended line's measured intensity one component carries.

    A measured line to which several transitions are assigned has one measured
    intensity for the lot of them.  The share given to one component is what
    the calculation says it contributes:

        share = I_calc(this component) / sum of I_calc over all components

    which is the quantity the pipeline stores in the BF column of
    line_classifications.csv and squares to weight blend components in the
    least-squares fit.  `sib` is the sum over the OTHER components, so the
    denominator here is i_calc + sib.  When one of the others has no calculated
    intensity there is nothing to compare it with and the line is split evenly,
    1/n_blend each - the same fallback the pipeline uses.

    `i_calc` is the value under the hypothesis being scored, not the one in the
    file: a hypothesis that makes a component fainter also makes it take a
    smaller piece of the blend, and that is part of what is being tested.  The
    scale correction of the wavelength band is common to all components and
    cancels in the ratio, so it is left out.
    """
    if n_blend <= 1:
        return 1.0
    if sib is None:
        return 1.0 / float(n_blend)
    den = i_calc + sib
    return (i_calc / den) if den > 0 else 1.0 / float(n_blend)


def level_residuals(lines, identity: str, gA: dict, g_imp: float,
                    bias: dict, C: float, kT: float):
    """The intensity evidence of one level's lines under one identity.

    `lines` are that level's accepted lines as (partner, wn, I_obs, E_up,
    n_blend, sib); `identity` is the level id whose gA values are used, which is
    the level itself under the ASSIGNED hypothesis and its partner in the pair
    under the SWAPPED one.

    A line that several transitions share contributes only the share of its
    measured intensity this component is calculated to carry, worked out under
    THIS hypothesis - see blend_share.

    Returns (x, y, r, w): the logarithm of the predicted intensity, the
    logarithm of the observed one, the residual y - x, and 1/sd^2 of the
    wavelength band, one entry per line.
    """
    x, y, r, w = [], [], [], []
    for partner, wn, i_obs, e_up, n_bl, sib in lines:
        g = gA.get(tuple(sorted((identity, partner))), g_imp)
        if not (g > 0):
            continue
        i_calc = predicted_intensity(g, wn, e_up, C, kT)
        if not (i_calc > 0):
            continue
        i_own = i_obs * blend_share(i_calc, n_bl, sib)
        if not (i_own > 0):
            continue
        f, sd = ls.scale_bias(wn, bias)
        lx = math.log(i_calc * f)
        ly = math.log(i_own)
        x.append(lx)
        y.append(ly)
        r.append(ly - lx)
        w.append(1.0 / (sd * sd) if sd > 0 else 1.0)
    return (np.asarray(x), np.asarray(y), np.asarray(r), np.asarray(w))


def common_lines(lines, id_a: str, id_b: str, gA: dict, g_imp: float):
    """The lines both hypotheses can score, so that lnR is a paired figure.

    A line to a partner for which one of the two identities has no predicted
    intensity at all is left out of both sides rather than counted against
    one; under the adopted policy missing_gA = "impute" every pair has one and
    nothing is dropped.  A line joining the two levels of the pair is dropped
    too - it cannot exist, since they share a parity - as is a line to a level
    that IS one of the pair, whose own gA would change under the swap.
    """
    def usable(identity):
        return {t[0] for t in lines
                if t[0] not in (id_a, id_b)
                and gA.get(tuple(sorted((identity, t[0]))), g_imp) > 0}

    ok = usable(id_a) & usable(id_b)
    return [t for t in lines if t[0] in ok]


def slope(x, y) -> float:
    """Ordinary least-squares slope of y on x; nan when it is undefined.

    Fewer than three points do not determine a slope worth reading (two always
    give an exact fit), and neither does a set of predictions that are all
    equal.
    """
    if len(x) < 3:
        return float('nan')
    var = float(np.var(x, ddof=1))
    if not (var > 0):
        return float('nan')
    return float(np.cov(x, y, ddof=1)[0, 1] / var)


def loglike(r, w) -> float:
    """-sum(r^2 / (2 sd^2)): the fit of one hypothesis to one level's lines."""
    return float(-0.5 * np.sum(w * r * r)) if len(r) else 0.0


def pair_evidence(id_a: str, id_b: str, by_level: dict, gA: dict,
                  g_imp: float, bias: dict, C: float, kT: float):
    """Compare the assignment of a level pair with its interchange.

    Returns None when either level has no line both hypotheses can score;
    otherwise a dict with the four slopes, the residual rms of each hypothesis
    and lnR = ln L(assigned) - ln L(swapped).

    lnR_A and lnR_B split lnR into the two levels' own contributions, and are
    worth reading separately: an interchange is symmetric and should improve
    both sides, whereas a large improvement on one level alone says that THAT
    level's label is wrong while the other's true identity may be a third
    level - very possibly one not yet observed, which this test cannot reach.
    """
    la = common_lines(by_level.get(id_a, []), id_a, id_b, gA, g_imp)
    lb = common_lines(by_level.get(id_b, []), id_a, id_b, gA, g_imp)
    if not la or not lb:
        return None
    aa = level_residuals(la, id_a, gA, g_imp, bias, C, kT)   # A as assigned
    ab = level_residuals(la, id_b, gA, g_imp, bias, C, kT)   # A under the swap
    bb = level_residuals(lb, id_b, gA, g_imp, bias, C, kT)   # B as assigned
    ba = level_residuals(lb, id_a, gA, g_imp, bias, C, kT)   # B under the swap
    if not (len(aa[0]) and len(ab[0]) and len(bb[0]) and len(ba[0])):
        return None
    r_as = np.concatenate((aa[2], bb[2]))
    r_sw = np.concatenate((ab[2], ba[2]))
    lnR_A = loglike(aa[2], aa[3]) - loglike(ab[2], ab[3])
    lnR_B = loglike(bb[2], bb[3]) - loglike(ba[2], ba[3])
    return {
        'n_A': int(len(aa[0])), 'n_B': int(len(bb[0])),
        'lnR_A': lnR_A, 'lnR_B': lnR_B,
        'slope_A_assigned': slope(aa[0], aa[1]),
        'slope_A_swapped': slope(ab[0], ab[1]),
        'slope_B_assigned': slope(bb[0], bb[1]),
        'slope_B_swapped': slope(ba[0], ba[1]),
        'rms_assigned': float(np.sqrt(np.mean(r_as ** 2))),
        'rms_swapped': float(np.sqrt(np.mean(r_sw ** 2))),
        'lnR': lnR_A + lnR_B,
    }


def slopes_prefer_swap(ev: dict) -> bool:
    """Are BOTH levels' slopes nearer 1 under the swap than as assigned?

    Undefined for a level with fewer than three lines, whose slope is nan;
    such a pair is left to the likelihood ratio and is not flagged
    automatically, so this returns False.
    """
    vals = (ev['slope_A_assigned'], ev['slope_A_swapped'],
            ev['slope_B_assigned'], ev['slope_B_swapped'])
    if any(math.isnan(v) for v in vals):
        return False
    return (abs(ev['slope_A_swapped'] - 1.0)
            < abs(ev['slope_A_assigned'] - 1.0)
            and abs(ev['slope_B_swapped'] - 1.0)
            < abs(ev['slope_B_assigned'] - 1.0))


# ---------------------------------------------------------------------------
# The search
# ---------------------------------------------------------------------------
def same_jp_blocks(levels: pd.DataFrame):
    """The level table split into blocks of one parity and one J.

    Only levels inside such a block can be interchanged: a swap that changed
    either quantity would change which lines the level may emit at all, and
    would be caught by the ordinary classification long before this test.
    """
    for _key, g in levels.groupby(['parity', 'J'], sort=True):
        yield g.reset_index(drop=True)


def pair_window(a, b, factor: float) -> float:
    """How far apart two levels may lie and still be interchangeable.

    The larger of the two configurations' rms of E_obs - E_calc, times
    `factor`: the calculation has to be able to misplace at least one of the
    two by the separation for the identities to be confusable.
    """
    w = [x for x in (a['W'], b['W']) if x == x]
    return factor * max(w) if w else float('nan')


def search_pairs(levels: pd.DataFrame, min_lines: int, window: float,
                 null_window: float):
    """(candidates, null pairs) as lists of (row A, row B, window).

    A candidate is a same-parity, same-J pair separated by no more than the
    pair window and supported by at least `min_lines` accepted lines on each
    side.  A null pair is the same but separated by MORE than the window and
    at most `null_window` times the configuration rms: a swap that the
    energies rule out, made between levels of the same part of the spectrum so
    that the two sets of gA values are of comparable magnitude.
    """
    cand, null = [], []
    sup = levels[levels['n_tot'] >= min_lines]
    for block in same_jp_blocks(sup):
        for i in range(len(block)):
            for j in range(i + 1, len(block)):
                a, b = block.iloc[i], block.iloc[j]
                w = pair_window(a, b, window)
                if not (w == w):
                    continue
                sep = abs(float(a['E_final']) - float(b['E_final']))
                if sep <= w:
                    cand.append((a, b, w))
                elif sep <= pair_window(a, b, null_window):
                    null.append((a, b, w))
    return cand, null


def mark_rows(rows: list, null_lnR, p_flag: float) -> list:
    """Add p_null, slopes_favour_swap and flagged to each evidence row.

    p_null is the share of the calibration swaps whose lnR is at or below this
    row's own; a pair is flagged when the swap fits better (lnR < 0), both
    slopes move nearer 1 under it, and p_null is at or below `p_flag`.  With no
    calibration sample p_null is nan and the first two conditions decide.
    """
    for r in rows:
        r['p_null'] = (float(np.mean(null_lnR <= r['lnR'])) if len(null_lnR)
                       else float('nan'))
        r['slopes_favour_swap'] = int(slopes_prefer_swap(r))
        r['flagged'] = int(r['lnR'] < 0 and r['slopes_favour_swap']
                           and (r['p_null'] <= p_flag
                                if r['p_null'] == r['p_null'] else True))
    return rows


def evidence_rows(pairs, by_level, gA, g_imp, bias, C, kT) -> list:
    """pair_evidence for each pair, with the identifying columns in front."""
    rows = []
    for a, b, w in pairs:
        ev = pair_evidence(a['level_id'], b['level_id'], by_level, gA, g_imp,
                           bias, C, kT)
        if ev is None:
            continue
        row = {
            'level_A': a['level_id'], 'level_B': b['level_id'],
            'J': a['J'], 'parity': a['parity'],
            'E_A': float(a['E_final']), 'E_B': float(b['E_final']),
            'dE_AB': float(b['E_final']) - float(a['E_final']),
            'window': w,
            'cfg_A': a['cfg'], 'term_A': a['term'],
            'cfg_B': b['cfg'], 'term_B': b['term'],
            'Ecalc_A': float(a['E_calc']), 'Ecalc_B': float(b['E_calc']),
            # How far the swap would move each level from its NEW calculated
            # position, in units of the rms of that configuration.
            'omc_swap_A_sigma': ((float(a['E_final']) - float(b['E_calc']))
                                 / float(b['W']) if b['W'] else float('nan')),
            'omc_swap_B_sigma': ((float(b['E_final']) - float(a['E_calc']))
                                 / float(a['W']) if a['W'] else float('nan')),
            'new_star_A': int(a.get('new_star', 0)),
            'new_star_B': int(b.get('new_star', 0)),
        }
        row.update(ev)
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def print_configuration_table(en: pd.DataFrame, windows: pd.Series,
                              factor: float):
    print('=' * 74)
    print('CALCULATED-ENERGY ERROR PER CONFIGURATION (from enlev.dat)')
    print('=' * 74)
    print('The rms of E_obs - E_calc over the experimentally known levels of')
    print('each dominant configuration.  Two levels of one parity and one J')
    print(f'closer than {factor:g} x this rms cannot be told apart by their')
    print('energies, and are candidates for an interchange.')
    print(f"  {'config':8s} {'levels':>7s} {'rms':>9s} {'max|O-C|':>10s} "
          f"{'window':>9s}")
    known = en[en['known']]
    for cfg, w in windows.sort_values(ascending=False).items():
        g = known[known['cfg'] == cfg]['omc']
        print(f'  {cfg:8s} {len(g):7d} {w:9.1f} {g.abs().max():10.1f} '
              f'{factor * w:9.1f}')


def print_pair(row: dict, flagged: bool):
    mark = ' <== FLAGGED' if flagged else ''
    print(f"\n  {row['level_A']} ({row['cfg_A']} {row['term_A']}, "
          f"E = {row['E_A']:.3f}, {row['n_A']} lines)")
    print(f"  {row['level_B']} ({row['cfg_B']} {row['term_B']}, "
          f"E = {row['E_B']:.3f}, {row['n_B']} lines)")
    print(f"    J = {row['J']}{row['parity']}, separation "
          f"{abs(row['dE_AB']):.2f} cm^-1, window {row['window']:.1f} cm^-1"
          f"{mark}")
    print(f"    slope of ln I_obs on ln I_calc  "
          f"A: {row['slope_A_assigned']:+.2f} -> "
          f"{row['slope_A_swapped']:+.2f}   "
          f"B: {row['slope_B_assigned']:+.2f} -> "
          f"{row['slope_B_swapped']:+.2f}   (assigned -> swapped)")
    print(f"    rms of ln(I_obs/I_calc)  {row['rms_assigned']:.2f} -> "
          f"{row['rms_swapped']:.2f};  lnR = {row['lnR']:+.1f}"
          f"  ({row['lnR_A']:+.1f} from A, {row['lnR_B']:+.1f} from B; "
          f"negative favours the swap);  p_null = {row['p_null']:.3f}")
    print(f"    the swap would put A at {row['omc_swap_A_sigma']:+.1f} and "
          f"B at {row['omc_swap_B_sigma']:+.1f} rms of its calculated "
          f"position")


def print_detail(id_a: str, id_b: str, by_level, gA, g_imp, bias, C, kT):
    """The line-by-line table behind one pair."""
    for lid, other in ((id_a, id_b), (id_b, id_a)):
        lines = common_lines(by_level.get(lid, []), id_a, id_b, gA, g_imp)
        print(f'\n{lid}: {len(lines)} accepted lines usable by both '
              f'hypotheses')
        if not lines:
            continue
        print(f"    {'partner':>14s} {'wn_obs':>12s} {'lambda':>9s} "
              f"{'I_obs':>10s} {'n':>2s} {'share as':>9s} {'share sw':>9s} "
              f"{'I_calc(as)':>11s} {'I_calc(sw)':>11s} "
              f"{'ln ratio as':>12s} {'ln ratio sw':>12s}")
        rows = sorted(lines, key=lambda t: -t[2])
        for partner, wn, i_obs, e_up, n_bl, sib in rows:
            g_as = gA.get(tuple(sorted((lid, partner))), g_imp)
            g_sw = gA.get(tuple(sorted((other, partner))), g_imp)
            f, _sd = ls.scale_bias(wn, bias)
            p_as = predicted_intensity(g_as, wn, e_up, C, kT)
            p_sw = predicted_intensity(g_sw, wn, e_up, C, kT)
            b_as = blend_share(p_as, n_bl, sib)
            b_sw = blend_share(p_sw, n_bl, sib)
            i_as, i_sw = p_as * f, p_sw * f
            print(f'    {partner:>14s} {wn:12.3f} {1.0e8 / wn:9.1f} '
                  f'{i_obs:10.2f} {n_bl:2d} {b_as:9.3f} {b_sw:9.3f} '
                  f'{i_as:11.2f} {i_sw:11.2f} '
                  f'{math.log(i_obs * b_as / i_as):12.2f} '
                  f'{math.log(i_obs * b_sw / i_sw):12.2f}')


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def build_context(args, verbose: bool = True):
    """Everything the comparison needs, read once.

    Returns (levels with identities, by_level, gA, g_imp, bias, C, kT, en,
    windows).
    """
    ls.LOPT_LINES = args.lopt
    ls.LOPT_LEVELS = args.lopt_levels
    ls.ENERGIES = args.energies
    ls.E_INPUT_CSV = args.e_input
    levels, real, _obs = ls.read_run()
    per = mc.per_level_table(real, levels)
    per['n_tot'] = per['n_old'] + per['n_new']
    per['new_star'] = (per['n_new'] > per['n_old']).astype(int)
    e_final = dict(zip(per['level_id'], per['E_final']))

    en = read_enlev(args.enlev)
    windows = configuration_windows(en)
    per, unmatched = attach_identities(per, en, windows)

    gA = read_gA()
    g_imp = imputed_gA()
    by_level = accepted_lines_by_level(real, e_final)

    # The predicted intensities the scale bias is measured on are the same
    # ones level_shifts.py uses, so the correction and the scatter that enter
    # lnR are exactly the numbers its report quotes.
    preds = ls.load_predictions(set(per['level_id']), e_final)
    bias = ls.intensity_scale_bias(real, preds, e_final, verbose)
    C = float(cl.CFG.intensity_model['C'])
    kT = float(cl.CFG.intensity_model['kT'])
    return per, unmatched, by_level, gA, g_imp, bias, C, kT, en, windows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description='Find level pairs whose theoretical identities may have '
                    'been interchanged.')
    ap.add_argument('--window', type=float, default=1.0,
                    help='search window, in units of the rms of E_obs - E_calc '
                         'of the configuration (default 1.0)')
    ap.add_argument('--null-window', type=float, default=10.0,
                    help='outer edge of the calibration sample, in the same '
                         'units (default 10)')
    ap.add_argument('--min-lines', type=int, default=3,
                    help='accepted lines each level of a pair must have '
                         '(default 3)')
    ap.add_argument('--p-flag', type=float, default=0.05,
                    help='flag a swap-favouring pair at or below this p_null '
                         '(default 0.05)')
    ap.add_argument('--enlev', default=ENLEV,
                    help='IDEN2 level file holding the labels and E_calc')
    ap.add_argument('--report', default=REPORT_CSV,
                    help='output csv (an .xlsx twin is written beside it)')
    ap.add_argument('--detail', nargs=2, metavar=('LEVEL_A', 'LEVEL_B'),
                    help='print the line-by-line table of one pair and stop')
    ap.add_argument('--lopt', default=None,
                    help='build the run from a LOPT line-output file instead '
                         'of the pipeline output')
    ap.add_argument('--lopt-levels', default=None,
                    help='LOPT level-output file, for --energies lopt')
    ap.add_argument('--energies', default='refit', choices=('refit', 'lopt'),
                    help="whose energies to use in --lopt mode "
                         "(default refit; see level_shifts.read_run)")
    ap.add_argument('--e-input', default=None,
                    help='csv of revised adopted energies (level_id,E_input)')
    args = ap.parse_args(argv)

    if not args.detail:
        output_files.require_writable(output_files.with_twin(args.report),
                                      'report file')
        print('Interchanged-identity check: are two levels of the same parity '
              'and J wearing')
        print("each other's calculated intensities?  Scale of the predicted "
              'intensities,')
        print('as level_shifts.py measures it on this run:')
    (per, unmatched, by_level, gA, g_imp, bias, C, kT,
     en, windows) = build_context(args, verbose=not args.detail)

    if args.detail:
        print_detail(args.detail[0], args.detail[1], by_level, gA, g_imp,
                     bias, C, kT)
        return 0

    print_configuration_table(en, windows, args.window)
    print()
    print(f"Levels in the run: {len(per)}; matched to a known level of "
          f"{os.path.basename(args.enlev)}: {len(per) - len(unmatched)}")
    if unmatched:
        base = os.path.basename(args.enlev)
        print(f"  no level marked found (*) in {base} lies within "
              f"{ENLEV_MATCH_TOL} cm^-1, so no label and no calculated "
              f"position:")
        e_of = dict(zip(per['level_id'], per['E_final']))
        for lid in unmatched:
            e = float(e_of[lid])
            row, d = nearest_enlev_row(en, e)
            mark = 'found (*)' if row['known'] else 'NOT marked found (*)'
            print(f"    {lid}  E = {e:.3f};  nearest row of {base} is "
                  f"{int(row['idx'])} {row['label']} at E_obs = "
                  f"{row['E_obs']:.3f} ({row['E_obs'] - e:+.3f}), {mark}")
        print("    a row at the right energy that is NOT marked found means "
              "the two files are out of step - fix enlev.dat rather than "
              "reading past this")
        print("    such a level still takes part - the swap exchanges gA "
              "values, which it has - but it borrows its partner's window and "
              "its calculated-position column is blank")
    if g_imp > 0:
        print(f"  a level pair absent from the calculated file is given the "
              f"imputed gA = {g_imp:.1f} s^-1")
    else:
        print("  missing_gA.policy = 'none': lines whose predicted intensity "
              "one hypothesis lacks are left out of both")

    cand, null = search_pairs(per, args.min_lines, args.window,
                              args.null_window)
    print(f"\nCandidate pairs (same parity and J, separated by at most "
          f"{args.window:g} x the configuration rms,")
    print(f"  both levels with at least {args.min_lines} accepted lines): "
          f"{len(cand)}")
    print(f"Calibration pairs (the same, but separated by more than the "
          f"window and at most {args.null_window:g} x the rms - swaps the "
          f"energies rule out): {len(null)}")

    null_rows = evidence_rows(null, by_level, gA, g_imp, bias, C, kT)
    null_lnR = np.array([r['lnR'] for r in null_rows], dtype=float)
    if len(null_lnR):
        print(f"\n  lnR of the {len(null_lnR)} calibration swaps: "
              f"median {np.median(null_lnR):+.1f}, "
              f"p05 {np.percentile(null_lnR, 5):+.1f}, "
              f"p01 {np.percentile(null_lnR, 1):+.1f}, "
              f"min {null_lnR.min():+.1f}")
        print(f"  {int((null_lnR < 0).sum())} of them ("
              f"{100.0 * (null_lnR < 0).mean():.1f}%) favour the swap at all, "
              f"which is the rate of accidental preference this test carries.")
    else:
        print("\n  no calibration swaps: p_null cannot be measured, and the "
              "flagging falls back on lnR and the slopes alone.")

    rows = mark_rows(evidence_rows(cand, by_level, gA, g_imp, bias, C, kT),
                     null_lnR, args.p_flag)
    rows.sort(key=lambda r: r['lnR'])

    # What the flagging rule does to swaps that are false by construction is
    # the expected number of false flags among the candidates, exactly as the
    # criterion grids of level_shifts.py quote expected false passes.  The
    # null pairs are ranked against their own distribution, which makes the
    # rate slightly conservative (each pair contributes to the sample it is
    # judged against).
    if len(null_lnR):
        false_rate = float(np.mean([r['flagged'] for r in
                                    mark_rows(null_rows, null_lnR,
                                              args.p_flag)]))
        print(f"  the same rule applied to those false swaps flags "
              f"{false_rate * 100:.2f}% of them, so among {len(rows)} "
              f"candidates {false_rate * len(rows):.2f} false flags are "
              f"expected.")

    print('\n' + '=' * 74)
    print(f'CANDIDATE PAIRS ({len(rows)}), most swap-favouring first')
    print('=' * 74)
    if not rows:
        print('  none: no two levels of the same parity and J lie close '
              'enough, with enough accepted lines on both sides.')
    for r in rows:
        print_pair(r, bool(r['flagged']))

    n_flag = sum(r['flagged'] for r in rows)
    print(f"\nFlagged: {n_flag} of {len(rows)} candidate pairs "
          f"(lnR < 0, both slopes nearer 1 under the swap, "
          f"p_null <= {args.p_flag:g}).")
    if n_flag:
        print("A flag asks for a look in IDEN2; it does not relabel a level.")
        print("To carry out an exchange once it has been decided on, use "
              "swap_line_assignments_LOPT.py")
        print("and swap_line_assignments_IDEN.py, which keep both level "
              "identifiers where they are.")

    if rows:
        out = pd.DataFrame(rows)
        mc.save_table(out, args.report, DECIMALS)
    return 0


if __name__ == '__main__':
    sys.exit(main())
