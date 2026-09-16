"""Validation of level confirmations: energy shifts, decoy calibration, and
per-level spurious probabilities.

Background. Wyart derived his new Pr III levels from real, previously
unclassified lines of the very line list used here, but his identifications
were never published and are lost. When the pipeline re-finds those lines,
the least-squares level optimization must reproduce the input energy of the
level almost exactly: the tested levels connect predominantly to old levels,
which are anchored by the legacy identifications and move very little. A
level supported by chance coincidences behaves differently: the optimizer
drags it to the weighted mean of arbitrary line positions, and its final
energy wanders away from the input value.

Which levels are tested. A level is "new*" if, in the real run, it is
supported by more new than old identifications (n_new > n_old). This is the
set whose reality needs validation; it comprises the levels absent from the
ASD compilation plus a few levels that are formally in ASD but whose
energies rest mainly on the new identifications. The remaining levels are
"old" and serve as the reference of genuine behavior.

Per-level quantities (from the real run):
    n  = n_old + n_new = number of accepted lines supporting the level;
    dE = E_final - E_input = optimized minus input energy;
    d  = dE * n, the support-normalized shift.
The two populations differ sharply in the size of the shift. Genuine levels
(the old ones) moved by |dE| of about 0.015 cm^-1 (median, nearly
independent of n; in the aggregate the old levels follow
rms(dE*n) ~ 0.7 cm^-1, a number dominated by the many-line levels). Levels
supported by false matches moved by about 0.2 cm^-1 (median, measured on
the decoys) with a wide spread. For given n, the |d| distribution of each
population is modeled as a log-normal whose median depends on n (a
straight-line fit of ln|d| against ln n), fitted to the old levels and to
the decoys, respectively.

Reference of spurious behavior: the decoy runs (decoy_mc.py). Each decoy is
a copy of a tested level whose energy is displaced so that none of its
transitions can coincide with a true wavenumber; every line accepted for a
decoy is a false match by construction. The decoys compete in the real run,
one per tested level.

Reference of genuine behavior: the old levels (anchored by legacy
identifications); and a second, independent reference of spurious behavior
is provided by the shifted-wavenumber runs (chance_mc.py), in which all
observed wavenumbers are displaced so that every acceptance is false. The
shifted runs leave the whole line list free for false matching (nothing is
consumed by true identifications), so their false-pass counts are upper
bounds; agreement between the two independent references shows that neither
is an artifact of how it was measured.

Per-level probability of being spurious. For each tested level two
independent pieces of evidence are weighed between the two hypotheses
("genuine" against "spurious", the latter meaning: behaves like a decoy):
  1. the energy shift: |d| follows the log-normal fitted to the old levels
     (genuine) or the log-normal fitted to the decoys (spurious);
  2. the intensity pattern: pattern_V falls into one of five ranges, whose
     probabilities for genuine and for spurious levels are measured on the
     old levels and on the decoys, respectively.
The only unknown is S, the total number of spurious levels among the tested
ones. For a trial S, the prior probability that a tested level with support
n is spurious is pi_n = S * q_n / m_n, where q_n is the decoy probability of
support n (given at least one line) and m_n is the observed number of tested
levels with support n; S is estimated by maximizing the likelihood of the
observed evidence. Three posterior probabilities are written to the report:
    p_spur_decoys  - from the energy shift alone;
    p_spur_pattern - from the intensity pattern alone;
    p_spur         - from both together (the folded value; the two pieces
                     of evidence are treated as independent for given
                     hypothesis and support n).
The sum of the folded p_spur over the tested levels reproduces the joint
estimate of S, and the expected number of false validations under any cut
equals the sum of p_spur over the levels passing the cut.

Adjudication of the questionable marks. Two things put a tested level up for
examination: the probabilities question it (folded p_spur >= P_QUESTIONABLE),
or its STRONGEST predicted transition is absent from the accepted lines with
nothing to explain the absence and no reason it would have escaped the plates
(report column top1 = 'missing'). Both are then
put through the masking check (question_status): a strong predicted transition
may be absent only because a nearby stronger observed line hides it on the
photographic plates. If the dominant part of the level's missing predicted
intensity is hidden in this way, the energy-shift evidence alone is not
suspect, and the strongest branch is not among the unexplained absences, the
questionable mark is cleared. The report columns question_status
('removed' = mark cleared, 'retained' = mark kept) and reason record the
outcome.

Why the strongest branch is judged separately. A level can score well on both
the energy shift and pattern_V and still be wrong, if the lines it explains
are weaker ones while its brightest predicted branch is simply not there.
pattern_V does not catch this: V asks whether the ACCEPTED lines are the
strongest of the predictions, so a level whose two accepted lines are ranked 2
and 3 can reach V near 1 with rank 1 missing. A missing, unmasked strongest
branch is the most direct evidence there is against a level, and it keeps the
questionable mark whatever the probabilities say. The top1 column reads
'found' (the strongest predicted transition is accepted), 'masked' (absent,
but an observed line hides it), 'faint' (absent, nothing hides it, but the
line would not reliably have been recorded anyway - see below), 'missing'
(absent, unexplained), or blank (the level has no observable prediction).

Absence is evidence only where a line would have been seen. Icalc is on the
same linear scale as the observed intensities from 1000 A upward, but not in
the far ultraviolet: below 900 A the median of I_obs/I_pred over the accepted
identifications is 0.09, and the ratio of the band totals - all observed
intensities over all predicted ones, which does not depend on which lines
were classified - agrees at 0.16. A prediction there is therefore about ten
times fainter on the plate than its I_pred says. intensity_scale_bias
measures that factor band by band, together with the scatter of
ln(I_obs/I_pred), which is of order 1 (a factor e either way) even where the
scale is right. The corrected intensity is what the observability filter
compares with the noise level. Sugar's plates do reach his intensity 1
everywhere - lines of intensity 1 are recorded and classified down to the
short-wavelength edge at 821.9 A - so this is a fault of the calculated
intensities, not of the plates.

Whether a predicted line would have been seen at all is observation_probability
= c(lambda) x D(z), both factors measured on this run rather than assumed:
c from tools/coverage_map.py, the probability that the plate was exposed and
measured at that wavelength, and D from tools/obscuration_rate.py, the observed
fraction of predictions of that strength that were recorded where the plate is
open. D has a ceiling at 1 - epsilon = 0.965, because about one predicted line
in thirty is lost to a grain flaw or a single impurity line wherever it falls.
A missing strongest branch counts as evidence (top1 = 'missing') only above
DETECT_CONFIDENCE, and reads 'faint' below it.

Intensity-pattern check. Independently of the energy shifts, a real level
must reproduce the PATTERN of theoretically predicted intensities: its
strongest predicted transitions should be present among the accepted lines
(the classical "square array" argument used by Wyart and his
contemporaries). For every level, all theoretically predicted transitions
to/from it (from Icalc.xlsx) with both partner levels known and the Ritz
wavenumber inside the observed range are collected and sorted by predicted
intensity, and two scores are computed:
    top10_found = how many of the 10 strongest predicted transitions are
                  among the accepted lines (n_top10 = how many were
                  available, if fewer than 10 are predicted);
    pattern_C   = sum of predicted intensities of the accepted transitions
                  divided by the sum over all predicted observable ones
                  (intensity-weighted completeness: 1 means every predicted
                  unit of intensity was found, 0 means none; missing a
                  strong predicted line lowers it strongly, missing a weak
                  one barely);
    pattern_V   = C divided by the largest C the level could have reached
                  with its number of matched lines (reached if they were
                  exactly the strongest predictions). V = 1 means the
                  accepted lines ARE the strongest predicted ones; V = 0
                  means that predictions exist but none of the accepted
                  lines is among them. Unlike C, V does not punish a level
                  for accepted lines that theory does not cover, and it is
                  nearly independent of the number of lines, which makes it
                  usable as a second, independent piece of evidence.
When the intensity-calibration file intensity_correction_functions.txt is
present (piecewise polynomials converting between Sugar's plate intensity
scale and the uniform linear scale), predicted transitions whose Sugar-scale
intensity falls below the noise level (Sugar intensity 1) are excluded from
the predicted-observable set before any score is computed - a prediction
that could not have been seen must not count as "missing". Accepted lines
whose predictions are below the noise level stay in the support count n
(the decoys are treated identically, so the calibration stays fair); the
report column n_acc_below_noise makes them visible.
The same scores computed for the decoys (false by construction) and for the
old levels (genuine) calibrate what fake and real patterns look like. Note
that the weeding already checks each individual line's intensity against
theory; the new information in these scores is the COVERAGE of the strong
predictions - the weeding never penalized a level for strong predicted
lines that are absent.

Where the run comes from. By default the accepted identifications and the
optimized level energies are read from the classify_lines.py output table
(line_classifications.csv), whose energies come from that program's own
internal weighted least squares. Identifications are also revised by hand -
in IDEN2, or by editing the LOPT input file - after which LOPT itself is
re-run; --lopt builds the same report from LOPT's line-output file instead
(see lopt_lines.py). For the same set of identifications the reconstruction
reproduces the classify table row for row.

Which energies a --lopt report uses, and why it is not LOPT's own by default.
The report judges a tested level by its energy shift dE against the same
quantity measured on levels that are false by construction: the decoys of
decoy_mc.py and the chance levels of chance_mc.py. Those reference populations
are made by re-running the pipeline, so their E_final comes from the
pipeline's optimizer. Taking E_final for the tested levels from LOPT instead
would measure the two sides of that comparison with different instruments:
LOPT tunes single-line levels, treats blends by its centroid model and rounds
its output, and its energies differ from the pipeline's by up to 0.2 cm^-1 -
the same size as the dE under test, and enough to widen the fitted scatter of
the genuine population by half again. So --lopt takes the ACCEPTED SET from
LOPT and re-solves the energies from it with the pipeline's own least squares
(lopt_lines.refit_energies). Fed the pipeline's own identifications the refit
reproduces the pipeline's energies to 1e-6 cm^-1 (median) and 0.0012 cm^-1
(maximum), so a --lopt report and the ordinary one differ only in the thing
being studied: which lines are assigned to which levels, which is what makes
a before/after comparison of a revision fair. Pass --energies lopt (with
--lopt-levels LOPT_output_levels_revised.txt, whose four-decimal energies are
better than the rounded E1/E2 of the line output) to report LOPT's own
optimized energies instead; those are the energies to publish, but they are
not comparable with the calibrations.

One thing does NOT follow a --lopt run automatically. E_input still comes from
the adopted-level workbook. dE asks whether the optimizer stays where the
identifications put the level, so for a level that a revision deliberately
MOVED, the workbook energy is the wrong reference and dE would report the size
of the move rather than test it. Give the new starting energies in a
two-column csv (level_id,E_input) via --e-input; the levels listed there are
reported, the rest keep the workbook value. Note that dE is then not an
independent test for those levels either, if the new E_input was itself
computed from the very lines the revision assigns to them: their evidence is
the intensity pattern and the residual spread, not dE.

Outputs: level_shift_report.csv/.xlsx (one row per level: energies, dE, d,
support counts, new* flag, p_spur, pattern scores, top1, and the two
observability columns top1_pobs and n_obs_expected) and the printed
calibration tables. With --lopt the report is named after the LOPT file
(LOPT_output_lines_revised.txt -> level_shift_report_revised.csv) so that it
never overwrites the pipeline's own; --report overrides the name.

Usage:
    python level_shifts.py                          # the full report
    python level_shifts.py --detail <level_id> ...  # inspect single levels:
        prints the level's predicted transitions (strongest first) and what
        happened to each in the real run - accepted, rejected (with the
        weeding note), or not matched by any observed line.  In --lopt mode
        the calculated intensity, grade and weeding note of a candidate are
        not available and print as blanks: LOPT does not record them.
    python level_shifts.py --lopt LOPT_output_lines_revised.txt \
                           --e-input revised_level_energies.csv
        the report of a hand-revised, LOPT-optimized run, written to
        level_shift_report_revised.csv.  Add
            --energies lopt --lopt-levels LOPT_output_levels_revised.txt
        to report LOPT's own energies instead of the calibration-consistent
        refit, and --report NAME to choose the output name.
"""
import math
import os
import sys

import numpy as np
import pandas as pd

import chance_mc as mc
import output_files
import classify_lines as cl
import decoy_mc as dmc
import lopt_lines

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT_CSV = os.path.join(HERE, 'level_shift_report.csv')

# Set from the command line by --lopt / --report / --e-input; see read_run()
# and the usage block of the module docstring.
LOPT_LINES = None    # LOPT line-output file to build the run from
LOPT_LEVELS = None   # LOPT level-output file, for --energies lopt
ENERGIES = 'refit'   # 'refit' (calibration-consistent) or 'lopt'; see read_run
E_INPUT_CSV = None   # csv of revised adopted energies (level_id,E_input)

SUPPORT_MAX_BIN = 8   # last support bin is "8 or more"
N_WORST = 12          # how many largest-|d| tested levels to list on the console
MAX_CONTAMINATION = 0.05  # acceptable expected-false share among validated levels
K_TOP = 10            # strongest predicted transitions defining top10_found
P_QUESTIONABLE = 0.1  # tested levels with p_spur >= this are listed in detail
INSTR_FWHM_A = 0.035  # angstrom: instrumental full width at half maximum of the
                      # grating spectrographs (nearly constant in wavelength;
                      # estimated at 0.03-0.04 A from the closest measured line
                      # pairs). The effective line width convolves this with
                      # the Doppler width, which grows with wavelength.
MASK_STRENGTH = 0.5   # within one effective line width the two lines are not
                      # resolved: a line of at least this fraction of the
                      # predicted intensity hides (absorbs) the prediction
MASK_CLEAR_FRACTION = 0.9  # share of the MISSING predicted intensity that must be
                           # masked for the questionable mark to be cleared
CALIB_FILE = os.path.join(HERE, 'intensity_correction_functions.txt')
# --- where a predicted line could have been seen ---------------------------
# Two separate things keep a predicted transition out of the line list, and
# both are now measured rather than assumed.
#
# COVERAGE, c(lambda): whether the plate was exposed and measured at that
# wavelength at all.  tools/coverage_map.py fits a hidden Markov model to the
# density of the observed lines against a smooth baseline and returns, bin by
# bin, the posterior probability that the record is open there; the result is
# coverage_map.csv, and coverage() below reads it.  It replaces the two
# hand-entered intervals that used to stand here, which came from the seams of
# the intensity calibration and are far too narrow for this purpose: the map
# finds ten blind stretches totalling 1078 A, and the two longest of them run
# out to 1663.82 A and 2190.91 A rather than stopping at 1529.85 A and
# 2107.92 A.  The old pair survives below only as the fallback used when
# coverage_map.csv is absent.  The same two numbers are also the forced region
# boundaries of the intensity calibration (DEFAULT_GAPS_A in
# tools/calibrate_intensities.py), and for THAT purpose - where the intensity
# scales of the bordering exposures could not be reconciled - they are right as
# they stand.  The two lists answer different questions and must not be kept
# equal.
#
# DETECTION, D(z): whether a line of that strength would have been recorded
# where the plate is open, with z = ln(I_pred x scale factor / noise level).
# tools/obscuration_rate.py measures it on the predicted transitions between
# two established levels, removing chance coincidences and dividing out the
# mapped blind stretches, and writes it into obscuration_rate.txt.  Its plateau
# is 1 - epsilon = 0.965: about one predicted line in thirty is lost even when
# it is bright and the plate is open, to a grain flaw or a single impurity
# line, and no map of any kind can say which one.
#
# The probability that a prediction would have been observed at all is the
# product
#     P_obs = c(lambda) x D(z)
# which is what observation_probability() returns.  The two factors do not
# double-count: D is measured only where c >= 0.95, so it describes the open
# plate and c carries everything else.
#
# When several branches of one level are absent, P_obs applies once per branch
# and not once for the level: tools/obscuration_length.py measures the
# obscuration to be 0.067 A wide, one resolution element, while the closest two
# branches of any established level ever come is 0.446 A.  The exception the
# same tool lists - two branches of one level closer together than 0.067 A - is
# one test and not two.
COVERAGE_FILE = os.path.join(HERE, 'coverage_map.csv')
DETECTION_FILE = os.path.join(HERE, 'obscuration_rate.txt')
COVERAGE_MIN = 0.5    # below this posterior the record is blind and a
                      # prediction there is dropped from the observable set.
                      # It is the threshold tools/coverage_map.py itself uses
                      # to list the blind stretches in coverage_gaps.txt, and
                      # the map is decisive nearly everywhere: three quarters
                      # of its bins are above 0.98 and half above 0.999, so
                      # little rides on the exact value.
COVERAGE_GAPS_A = [(1522.49, 1529.85), (2103.46, 2107.92)]  # fallback only

# Wavelength bands in which the predicted intensities are checked against the
# observed ones for a scale error (see intensity_scale_bias).  Sugar's plates
# reach his intensity 1 everywhere - lines of intensity 1 are recorded and
# classified down to the short-wavelength edge at 821.9 A - but the CALCULATED
# intensities in the far ultraviolet are an order of magnitude too large, so a
# prediction there is far weaker on the plate than its I_pred suggests.
SCALE_BIAS_BANDS_A = [(0.0, 900.0), (900.0, 1000.0)]
SCALE_BIAS_MIN_N = 10          # accepted lines needed to measure a band
SCALE_BIAS_MAX_FACTOR = 0.5    # only a bias this strong is worth correcting
# A missing strongest branch counts as evidence against a level only if the
# line would have been observed with at least this probability (top1_status).
#
# The probability is observation_probability below - the coverage of the plate
# times the MEASURED detection curve - and no longer the log-normal around the
# noise level that used to stand here.  That raises the bar a long way, which
# is the point of the change.  The log-normal was fitted to the ACCEPTED lines,
# whose lower tail the noise has already removed, so it was optimistic exactly
# where it mattered: it put the recording probability at 0.99 by ten times the
# noise level, where the measured rate is 0.81.  The measured D does not reach
# 0.9 until about sixty times the noise level, so this threshold now selects
# what it always claimed to select - branches whose absence really is
# unexplained - instead of admitting lines that go unrecorded one time in five.
# It cannot be raised much further: D has a ceiling at 1 - epsilon = 0.965, and
# a threshold above that would make the test unpassable at any brightness.
DETECT_CONFIDENCE = 0.9


def q(x, p):
    """Percentile that tolerates empty input: the value below which p percent
    of the sample falls."""
    x = np.asarray(x, dtype=float)
    return float(np.percentile(x, p)) if x.size else float('nan')


def support_bin(n):
    return np.minimum(np.asarray(n, dtype=int), SUPPORT_MAX_BIN)


def dE_by_n_table(sup: pd.DataFrame, ref: pd.DataFrame, count_col: str,
                  ref_name: str):
    """|dE| of real tested levels vs the false-by-construction reference.

    `ref` holds the reference levels (decoys, or levels of the shifted-line
    runs), whose accepted lines are all false matches; `count_col` names the
    column with their accepted-line counts. For each support n the table
    lists the median |dE| (half the cases are below it) and the p90 value
    (90 percent of the cases are below it).
    """
    print(f"\n|dE| by support n: real tested levels vs {ref_name} (runs pooled):")
    print(f"  {'n':>3}   {'real: N':>8} {'median':>8} {'p90':>8}   |"
          f"   {ref_name + ': N':>12} {'median':>8} {'p90':>8}")
    for n in range(1, SUPPORT_MAX_BIN + 1):
        if n < SUPPORT_MAX_BIN:
            r = sup.loc[sup['n_tot'] == n, 'abs_dE']
            c = ref.loc[ref[count_col] == n, 'abs_dE']
            label = f"{n}"
        else:
            r = sup.loc[sup['n_tot'] >= n, 'abs_dE']
            c = ref.loc[ref[count_col] >= n, 'abs_dE']
            label = f"{n}+"
        print(f"  {label:>3}   {len(r):8d} {q(r, 50):8.4f} {q(r, 90):8.4f}   |"
              f"   {len(c):12d} {q(c, 50):8.4f} {q(c, 90):8.4f}")


def criterion_grid(new: pd.DataFrame, ref: pd.DataFrame, count_col: str,
                   n_runs: int, Ts: list, n_tested: int, ref_name: str):
    """Tested levels passing (n >= N, |dE| <= T) vs expected false passes.

    For every combination of N (minimum number of supporting lines) and T
    (largest allowed |dE|), prints two numbers: how many tested levels pass
    the cut, and how many of the false-by-construction reference levels pass
    it per run - i.e. how many false validations to expect. Returns the cut
    with the largest real yield whose expected share of false validations
    stays below MAX_CONTAMINATION, or None if no cut reaches that.
    """
    print(f"\nCriterion (n >= N and |dE| <= T): tested levels passing / "
          f"expected false passes per run ({ref_name}):")
    print("  N \\ T " + "".join(f"{name:>16}" for name, _ in Ts))
    best = None
    for N in range(1, SUPPORT_MAX_BIN + 1):
        cells = []
        for name, t in Ts:
            rp = int(((new['n_tot'] >= N) & (new['abs_dE'] <= t)).sum())
            ce = float(((ref[count_col] >= N) &
                        (ref['abs_dE'] <= t)).sum()) / n_runs
            cells.append(f"{rp:6d} /{ce:6.2f}")
            if rp > 0 and ce / rp <= MAX_CONTAMINATION and (best is None or rp > best[2]):
                best = (N, name, rp, ce, t)
        print(f"  >= {N} " + "".join(f"{c:>16}" for c in cells))
    if best is not None:
        N, name, rp, ce, t = best
        print(f"Best cut at <= {MAX_CONTAMINATION:.0%} expected contamination: "
              f"n >= {N} and |dE| <= {name} = {t:.4f} cm^-1 "
              f"-> {rp} of {n_tested} levels validated, expected false: {ce:.2f}")
    else:
        print(f"No cut reaches <= {MAX_CONTAMINATION:.0%} expected contamination.")
    return best


def read_intensity_calibration(path: str = CALIB_FILE):
    """Piecewise-polynomial intensity calibration of Sugar's plates.

    Each line of the file holds three whitespace-separated fields:
    lambda_start, lambda_end (vacuum wavelength, angstroms) and the
    polynomial coefficients c0;c1;...;cn in ASCENDING powers of the vacuum
    wavelength. Sugar's reported plate intensity converts to the uniform
    linear scale as I_linear = 1000 * I_Sugar * exp(P(lambda)); the factor
    1000 was applied when the polynomials were derived, so that the
    linearized intensities in Pr3_lines.xlsx are integers (smallest 21).
    Consequently a theoretical intensity (linear scale) converts back to
    Sugar's scale as I_calc / (1000 * exp(P(lambda))), and the linear
    intensity corresponding to Sugar's noise level (intensity 1) is
    1000 * exp(P(lambda)).

    Returns a sorted list of (lambda_start, lambda_end, coefficients), or
    None if the file is absent.
    """
    if not os.path.exists(path):
        return None
    regions = []
    with open(path) as f:
        for raw in f:
            raw = raw.strip()
            if not raw:
                continue
            parts = raw.split()
            regions.append((float(parts[0]), float(parts[1]),
                            [float(x) for x in parts[2].split(';')]))
    regions.sort()
    return regions


def noise_threshold_linear(wn: float, calib) -> float | None:
    """Linear-scale intensity equal to Sugar's noise level (intensity 1) at
    the vacuum wavelength corresponding to wavenumber wn; None outside the
    calibrated wavelength range."""
    lam = 1.0e8 / wn
    if lam < calib[0][0] or lam > calib[-1][1]:
        return None
    coeffs = None
    for lo, hi, cf in calib:
        if lam <= hi:      # first region whose end is not passed; the tiny
            coeffs = cf    # gaps between regions fall to the next region
            break
    if coeffs is None:
        coeffs = calib[-1][2]
    p = 0.0
    for k, ck in enumerate(coeffs):
        p += ck * lam ** k
    return float(np.exp(p))*1000


_coverage_map = None    # (bin centres, coverage) of coverage_map.csv, cached


def read_coverage_map(path: str = COVERAGE_FILE):
    """The coverage function c(lambda) as (bin centres, coverage).

    coverage_map.csv is written by tools/coverage_map.py: one row per
    wavelength bin, with the posterior probability that the record is open
    there.  Returns () if the file is absent, which sends coverage() back to
    the two hand-entered intervals of COVERAGE_GAPS_A.
    """
    global _coverage_map
    if _coverage_map is None:
        if os.path.exists(path):
            m = pd.read_csv(path)
            _coverage_map = (0.5 * (m['lam_lo'] + m['lam_hi']).to_numpy(),
                             m['coverage'].to_numpy(dtype=float))
        else:
            _coverage_map = ()
            print(f"WARNING: {os.path.basename(path)} not found - falling back "
                  f"on COVERAGE_GAPS_A, which understates the blind stretches "
                  f"by a factor of ninety in total width; run "
                  f"tools/coverage_map.py")
    return _coverage_map


def coverage(wn: float) -> float:
    """c(lambda): the probability that the record is open at the vacuum
    wavelength of this wavenumber.

    Linearly interpolated between the bin centres of the map and held at the
    end values outside its range, which runs 822 - 10719 A, the extremes of the
    observed line list.  Without the map this degrades to 0 inside
    COVERAGE_GAPS_A and 1 outside it.
    """
    lam = 1.0e8 / wn
    m = read_coverage_map()
    if not m:
        return 0.0 if any(lo <= lam <= hi for lo, hi in COVERAGE_GAPS_A) else 1.0
    return float(np.interp(lam, m[0], m[1]))


def in_coverage_gap(wn: float) -> bool:
    """True where the record is blind - coverage below COVERAGE_MIN - so that
    nothing at this wavelength could have been observed and the absence of a
    predicted line there says nothing at all."""
    return coverage(wn) < COVERAGE_MIN


def intensity_scale_bias(real: pd.DataFrame, preds: pd.DataFrame,
                         e_final: dict, verbose: bool = True) -> dict:
    """Scale error and scatter of the predicted intensities, band by band.

    Icalc is meant to be on the same linear scale as the observed
    intensities, and from 1000 A upward it is: the median of I_obs / I_pred
    over the accepted identifications is within about 30 percent of 1.  Below
    900 A it is not - the median falls to about 0.09, and the ratio of the
    band totals (all observed intensities over all predicted ones, which does
    not depend on which lines were classified) agrees at 0.16.  A predicted
    transition there is therefore about ten times fainter on the plate than
    its I_pred says, and judging it against the noise level on the
    uncorrected scale calls lines observable that are in fact at the
    detection limit.

    The scatter matters as much as the factor: ln(I_obs/I_pred) has a
    standard deviation of order 1 (a factor e either way) even where the
    scale is right, because gA values of this size carry uncertainties of
    tens of percent and Sugar's eye estimates are coarse.  It is measured
    here so that the detectability of a single predicted line can be stated
    as a probability rather than a yes or no (see top1_status).

    Returns {(lambda_low, lambda_high): (factor, sd)} for the bands of
    SCALE_BIAS_BANDS_A measured on at least SCALE_BIAS_MIN_N accepted lines,
    plus the entry None: (1.0, sd) for every other wavelength.  A band whose
    factor is not below SCALE_BIAS_MAX_FACTOR gets factor 1.0: the intensity
    scale there is right and nothing is corrected.
    """
    ip = {tuple(sorted((lo, up))): v
          for lo, up, v in preds.itertuples(index=False)}
    acc = real[real['accepted'] == 1]
    # A measured line that several transitions share carries one intensity for
    # all of them.  What is compared with ONE transition's predicted intensity
    # is that transition's share of it - the branching fraction the pipeline
    # stores in the BF column (1 for an unblended line, and where the column is
    # missing, as in a run rebuilt from LOPT output alone).  Comparing the whole
    # feature with one component would put a blended line above the scale by the
    # reciprocal of its share, which for the faint components of a blend is a
    # factor of ten or more.
    bf_col = (pd.to_numeric(acc['BF'], errors='coerce').fillna(1.0)
              if 'BF' in acc.columns else pd.Series(1.0, index=acc.index))
    lam, res = [], []
    for lo, up, wn, io_, bf in zip(acc['low_id'], acc['upp_id'],
                                   acc['wn_obs'], acc['obs_intens'], bf_col):
        v = ip.get(tuple(sorted((str(lo), str(up)))))
        if v is None or not (v > 0) or not (io_ > 0):
            continue
        bf = float(bf)
        if not (bf > 0):
            bf = 1.0
        lam.append(1.0e8 / float(wn))
        res.append(math.log(float(io_) * bf / float(v)))
    lam = np.asarray(lam, dtype=float)
    res = np.asarray(res, dtype=float)
    out, covered = {}, np.zeros(len(lam), dtype=bool)
    for lo_a, hi_a in SCALE_BIAS_BANDS_A:
        m = (lam >= lo_a) & (lam < hi_a)
        if int(m.sum()) < SCALE_BIAS_MIN_N:
            continue
        covered |= m
        f = float(np.exp(np.median(res[m])))
        sd = float(np.std(res[m], ddof=1))
        if f >= SCALE_BIAS_MAX_FACTOR:
            f = 1.0
        out[(lo_a, hi_a)] = (f, sd)
        if verbose:
            print(f"  intensity scale {lo_a:.0f}-{hi_a:.0f} A: median "
                  f"I_obs/I_pred = {np.exp(np.median(res[m])):.3f} "
                  f"(sd of ln = {sd:.2f}) over {int(m.sum())} accepted lines"
                  + ("" if f != 1.0 else "  - no correction"))
    rest = res[~covered]
    sd0 = float(np.std(rest, ddof=1)) if len(rest) > 2 else 1.0
    out[None] = (1.0, sd0)
    if verbose:
        print(f"  intensity scale elsewhere: no correction "
              f"(sd of ln = {sd0:.2f} over {len(rest)} accepted lines)")
    return out


def scale_bias(wn: float, bias: dict) -> tuple:
    """(factor, sd of ln) of intensity_scale_bias applying at this wavenumber."""
    if not bias:
        return 1.0, 1.0
    lam = 1.0e8 / wn
    for band, v in bias.items():
        if band is not None and band[0] <= lam < band[1]:
            return v
    return bias.get(None, (1.0, 1.0))


def scale_bias_factor(wn: float, bias: dict) -> float:
    """The intensity-scale factor applying at this wavenumber."""
    return scale_bias(wn, bias)[0]


_detection_curve = None    # (z anchors, D anchors) of obscuration_rate.txt


def read_detection_curve(path: str = DETECTION_FILE):
    """The measured detection curve D(z) as (z anchors, D anchors).

    obscuration_rate.txt gives D in eight bins of z = ln(I_pred / noise level),
    the chance- and coverage-corrected fraction of predicted transitions that
    were recorded or masked.  Each bounded bin is anchored at its midpoint and
    the two open end bins at their finite edge; detection_probability
    interpolates linearly between the anchors and holds the end values beyond
    them, so D is flat at the measured plateau above a hundred times the noise
    level and flat at the measured floor a tenth of it.  Returns () if the file
    is absent or carries no curve.
    """
    global _detection_curve
    if _detection_curve is not None:
        return _detection_curve
    _detection_curve = ()
    if not os.path.exists(path):
        print(f"WARNING: {os.path.basename(path)} not found - falling back on "
              f"the log-normal detection model, which is optimistic at the "
              f"faint end; run tools/obscuration_rate.py")
        return _detection_curve
    z, d, inside = [], [], False
    with open(path) as fh:
        for line in fh:
            t = line.strip()
            if t.startswith('# the detection curve'):
                inside = True
            elif inside and not t:
                break
            elif inside and not t.startswith('#'):
                w = t.split()
                if len(w) < 5:
                    continue
                lo, hi, dz = float(w[0]), float(w[1]), float(w[-2])
                # the strength range in the middle carries spaces of its own,
                # so the numbers are taken from the two ends of the row
                z.append(hi if lo == -math.inf else
                         lo if hi == math.inf else 0.5 * (lo + hi))
                d.append(min(max(dz, 0.0), 1.0))
    if len(z) >= 2:
        o = np.argsort(z)
        _detection_curve = (np.asarray(z)[o], np.asarray(d)[o])
    return _detection_curve


def detection_probability(wn: float, i_pred: float, calib, bias: dict) -> float:
    """Probability that a line this bright would have been recorded WHERE THE
    PLATE IS OPEN - the second factor of observation_probability.

    The line reaches the plate with intensity I_pred times the scale factor of
    its band, and z = ln(that / the noise level at the same wavelength).  The
    answer is the measured curve of read_detection_curve, whose ceiling is
    1 - epsilon = 0.965 rather than 1: even a line a thousand times the noise
    level is lost about one time in thirty.

    The scale factor is applied here but was not applied when the curve was
    measured, and the two agree because the factor is 1 above 1000 A and the
    curve is a curve for that region: the 243 predictions below 1000 A in the
    measurement move no bin of it by as much as one standard error.  Correcting
    z here therefore maps a far-ultraviolet prediction onto the same curve
    instead of reading it at a z that is too high by ln(10).

    Falls back on the old log-normal around the noise level, of width the
    measured sd of ln(I_obs/I_pred), when obscuration_rate.txt is not there.
    Returns 1.0 where there is no calibration to compare with.
    """
    if calib is None:
        return 1.0
    thr = noise_threshold_linear(wn, calib)
    if thr is None or thr <= 0 or i_pred <= 0:
        return 1.0
    f, sd = scale_bias(wn, bias)
    z = math.log(i_pred * f / thr)
    cur = read_detection_curve()
    if cur:
        return float(np.interp(z, cur[0], cur[1]))
    if not (sd > 0):
        return 1.0 if z >= 0.0 else 0.0
    return 0.5 * (1.0 + math.erf(z / (sd * math.sqrt(2.0))))


def observation_probability(wn: float, i_pred: float, calib,
                            bias: dict) -> float:
    """P_obs = c(lambda) x D(z): the probability that this predicted transition
    would have appeared in the observed line list at all.

    Both factors are measured - the coverage by tools/coverage_map.py, the
    detection curve by tools/obscuration_rate.py - and they are independent by
    construction, the curve having been measured only where the coverage is at
    least 0.95.  This is the quantity a missing-line term needs: the absence of
    a prediction is evidence against a level in proportion to how surely it
    would have been seen, and it is no evidence at all where the plate is blind
    or the line is at the noise.
    """
    return coverage(wn) * detection_probability(wn, i_pred, calib, bias)


def filter_observable(preds: pd.DataFrame, e_final: dict, calib,
                      bias: dict | None = None):
    """Drop predicted transitions that could not have been observed.

    Two causes: (a) the Ritz wavelength falls where the record is blind, that
    is where the measured coverage is below COVERAGE_MIN; (b) the Sugar-scale
    intensity
    I_pred/exp(P(lambda)) is below the noise level 1, i.e. I_pred is below
    the linear noise threshold at the Ritz wavelength (requires the
    calibration; predictions outside the calibrated range are kept).

    `bias`, from intensity_scale_bias, corrects the known scale error of the
    predicted intensities in the far ultraviolet before that comparison: what
    has to clear the noise level is the intensity the line would actually
    have on the plate, I_pred times the factor of its band.  The pattern
    scores keep using the uncorrected I_pred, which is a relative weight
    among a level's branches, not a plate intensity.
    Returns (kept, dropped_below_noise, n_in_gaps).
    """
    gap = np.array([in_coverage_gap(e_final[up] - e_final[lo])
                    for lo, up, _ in preds.itertuples(index=False)], dtype=bool)
    kept = preds[~gap].reset_index(drop=True)
    if calib is None:
        return kept, preds.iloc[0:0], int(gap.sum())
    keep = []
    for lo, up, ip in kept.itertuples(index=False):
        wn = e_final[up] - e_final[lo]
        thr = noise_threshold_linear(wn, calib)
        keep.append(thr is None or ip * scale_bias_factor(wn, bias) >= thr)
    keep = np.array(keep, dtype=bool)
    return (kept[keep].reset_index(drop=True),
            kept[~keep].reset_index(drop=True), int(gap.sum()))


def load_predictions(level_ids: set, e_final: dict) -> pd.DataFrame:
    """Theoretically predicted observable transitions, from Icalc.xlsx.

    Keeps entries where both partner levels are in the level list, the
    predicted intensity is present, and the Ritz wavenumber (from the final
    optimized energies) lies inside the observed range [WN_MIN, WN_MAX].
    Returns a DataFrame with columns lo_id, up_id, I_pred.

    A level added by files.new_levels has predicted transitions here too: from
    the Cowan transition list through its cowan_lid (_cowan_predictions), and
    from the supplementary file files.icalc_extra where one is still named.
    Without it such a level would arrive at the intensity-pattern test with no
    predictions at all, and the test - which asks how much of what theory
    expects to see was in fact found - would return nothing for precisely the
    level whose reality is least established.  Its Icalc is recomputed from gA
    for the same reason as in classify_lines.read_extra_transitions(): the
    column of that file is not rewritten when C and kT move.
    """
    import openpyxl
    rows = []
    files = [(cl.ICALC_FILE, False)]
    if getattr(cl, 'ICALC_EXTRA', ''):
        files.append((cl.ICALC_EXTRA, True))
    for path, recompute in files:
        rows.extend(_read_predictions(path, level_ids, e_final, recompute))
    rows.extend(_cowan_predictions(level_ids, e_final, rows))
    return pd.DataFrame(rows, columns=['lo_id', 'up_id', 'I_pred'])


def _cowan_predictions(level_ids, e_final, have):
    """The predictions of the levels files.new_levels adds, from the Cowan list.

    Such a level has no row in Icalc.xlsx, and files.icalc_extra - where its
    transitions used to be written by hand - is retired: classify_lines.py
    takes them from tp_E1_no_trials.xlsx through the cowan_lid column instead
    (classify_lines.cowan_transitions_of).  This reads them the same way, with
    the same gA cutoff and the same intensity relation, so that the level is
    predicted here exactly as it was when its lines were classified.  Without
    it every such level had no predictions at all, and level_positions.py,
    which scans only levels that have some, left all of them out.

    The Ritz wavenumber is taken from e_final, as for every other row, and a
    pair already read from Icalc.xlsx or files.icalc_extra is not added again.
    """
    lids = {k: v for k, v in cl.new_level_cowan_lids().items() if k in level_ids}
    if not lids:
        return []
    C = float(cl.CFG.intensity_model['C'])
    kT = float(cl.CFG.intensity_model['kT'])
    seen = {tuple(sorted((a, b))) for a, b, _ in have}
    found, _n_cut, _n_unknown = cl.cowan_transitions_of(lids, set(level_ids))
    out = []
    for new_id, partner_id, gA, _u in found:
        lo, up = sorted((new_id, partner_id), key=lambda x: e_final[x])
        if (tuple(sorted((lo, up))) in seen or lo not in e_final
                or up not in e_final):
            continue
        wn = e_final[up] - e_final[lo]
        if wn < cl.WN_MIN or wn > cl.WN_MAX:
            continue
        i_pred = C * gA * (wn / 1e8) * math.exp(-e_final[up] / kT)
        if i_pred <= 0:
            continue
        seen.add(tuple(sorted((lo, up))))
        out.append((lo, up, i_pred))
    return out


def _read_predictions(path, level_ids, e_final, recompute):
    import openpyxl
    C = float(cl.CFG.intensity_model['C'])
    kT = float(cl.CFG.intensity_model['kT'])
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[cl.CFG.icalc.sheet]
    col = cl.column_index(ws, cl.CFG.icalc, path)
    rows = []
    for row in ws.iter_rows(min_row=2):
        id1 = cl.to_str_id(row[col['id1']].value)  # lower level
        id2 = cl.to_str_id(row[col['id2']].value)  # upper level
        icalc = row[col['Icalc']].value            # predicted intensity, observed scale
        if not id1 or not id2 or icalc is None:
            continue
        if id1 not in level_ids or id2 not in level_ids:
            continue
        try:
            if recompute:
                i_pred = (C * float(row[col['gA']].value)
                          * (float(row[col['rwn']].value) / 1e8)
                          * math.exp(-float(row[col['Eup']].value) / kT))
            else:
                i_pred = float(icalc)
        except (TypeError, ValueError):
            continue
        if i_pred <= 0:
            continue
        wn = e_final[id2] - e_final[id1]
        if wn < cl.WN_MIN or wn > cl.WN_MAX:
            continue
        rows.append((id1, id2, i_pred))
    wb.close()
    return rows


def pattern_scores(level_ids, accepted_pairs: set, preds: pd.DataFrame) -> dict:
    """Intensity-pattern scores of each level against the predictions.

    accepted_pairs holds the accepted transitions of one run as unordered
    level-id pairs (tuple(sorted((lo, up)))). For each level the function
    returns (n_pred, n_topK, topK_found, pattern_C, pattern_V):
      n_pred      = number of predicted observable transitions of the level;
      n_topK      = min(K_TOP, n_pred), the number of strongest predictions
                    considered;
      topK_found  = how many of those n_topK are among the accepted lines;
      pattern_C   = intensity-weighted completeness: the sum of predicted
                    intensities over the accepted transitions divided by the
                    sum over all n_pred predicted ones;
      pattern_V   = C divided by the largest C achievable with the level's
                    number of matched predictions (V = 1: the accepted lines
                    are exactly the strongest predictions; V = 0: none of
                    the accepted lines is among the predictions).
    """
    by_level = {}
    for lo, up, ip in preds.itertuples(index=False):
        pair = tuple(sorted((lo, up)))
        by_level.setdefault(lo, []).append((ip, pair))
        by_level.setdefault(up, []).append((ip, pair))
    out = {}
    for lid in level_ids:
        plist = sorted(by_level.get(lid, []), key=lambda t: -t[0])
        if not plist:
            out[lid] = (0, 0, 0, np.nan, np.nan)
            continue
        k = min(K_TOP, len(plist))
        top_found = sum(1 for ip, pair in plist[:k] if pair in accepted_pairs)
        i_all = sum(ip for ip, _ in plist)
        matched = [ip for ip, pair in plist if pair in accepted_pairs]
        c = sum(matched) / i_all
        if matched:
            c_max = sum(ip for ip, _ in plist[:len(matched)]) / i_all
            v = c / c_max
        else:
            v = 0.0
        out[lid] = (len(plist), k, top_found, c, v)
    return out


def accepted_pair_set(df: pd.DataFrame, strip_decoy_prefix: bool = False) -> set:
    """Unordered level-id pairs of the accepted rows of a run output."""
    acc = df[df['accepted'] == 1] if 'accepted' in df.columns else df
    pairs = set()
    for lo, up in zip(acc['low_id'], acc['upp_id']):
        lo, up = str(lo), str(up)
        if strip_decoy_prefix:
            lo = lo[len(cl.DECOY_PREFIX):] if lo.startswith(cl.DECOY_PREFIX) else lo
            up = up[len(cl.DECOY_PREFIX):] if up.startswith(cl.DECOY_PREFIX) else up
        pairs.add(tuple(sorted((lo, up))))
    return pairs


def fit_lognormal_d(n_arr, d_abs_arr, b_fixed=None):
    """Fit ln|d| = a + b*ln(n) with Gaussian scatter s around the line.

    Returns (a, b, s, n_points) and thereby defines a log-normal probability
    density of |d| at given n. Points with |d| = 0 are dropped. If b_fixed is
    given, the slope is held at that value and only a and s are fitted
    (b_fixed = 1 corresponds to |dE| independent of n).
    """
    n_arr = np.asarray(n_arr, dtype=float)
    d_abs_arr = np.asarray(d_abs_arr, dtype=float)
    keep = d_abs_arr > 0
    x = np.log(n_arr[keep])
    y = np.log(d_abs_arr[keep])
    if b_fixed is None:
        b, a = np.polyfit(x, y, 1)
    else:
        b = float(b_fixed)
        a = float((y - b * x).mean())
    s = float((y - (a + b * x)).std(ddof=2))
    return float(a), float(b), s, int(keep.sum())


def lognormal_density(d_abs, n, a, b, s):
    d_abs = max(float(d_abs), 1e-6)
    mu = a + b * np.log(float(n))
    return (np.exp(-0.5 * ((np.log(d_abs) - mu) / s) ** 2)
            / (d_abs * s * np.sqrt(2 * np.pi)))


V_EDGES = [0.0, 0.1, 0.3, 0.6, 0.9, 1.0000001]
V_LABELS = ['0-0.1', '0.1-0.3', '0.3-0.6', '0.6-0.9', '0.9-1']


def v_bin_index(v):
    return int(np.clip(np.digitize([float(v)], V_EDGES)[0] - 1,
                       0, len(V_EDGES) - 2))


def v_bin_probs(v_values):
    """Probability of each pattern_V range, with 0.5 counts added to every
    range so that none has zero probability."""
    v = np.asarray(v_values, dtype=float)
    v = v[np.isfinite(v)]
    cnt = np.histogram(v, bins=V_EDGES)[0].astype(float) + 0.5
    return cnt / cnt.sum()


def fit_spurious_probabilities(tested: pd.DataFrame, decoys: pd.DataFrame,
                               old_cal: pd.DataFrame, decoy_pattern=None,
                               gen_b_fixed=None):
    """Per-level spurious probabilities from the energy shifts and,
    optionally, the intensity patterns.

    tested: the new* levels of the real run (columns n_tot, abs_d and, for
            the pattern evidence, pattern_V, n_pred).
    decoys: all decoy rows of decoy_mc_levels.csv (n_decoy_lines, dE).
    old_cal: old levels with >= 1 accepted line - the reference of genuine
            behavior (columns n_tot, abs_d, pattern_V, n_pred).
    decoy_pattern: per-trial pattern scores of the supported decoys (column
            V); if None, only the energy-shift evidence is used.
    gen_b_fixed: if given, hold the slope of the genuine ln|d| vs ln n line
            at this value instead of fitting it (robustness check of the
            extrapolation to small n, where old calibration levels are
            scarce).

    Returns a dict with one entry per evidence stream - 'dE' and, when
    decoy_pattern is given, 'pattern' and 'joint' - each holding
    (p_spur array aligned with tested, S_best, (S_lo95, S_hi95)); the key
    'info' holds the fitted model parameters. Each stream estimates its own
    S by maximizing its likelihood, so the three posteriors are separately
    self-consistent.
    """
    # --- spurious model, part 1: support distribution of a fake level ---
    # probability that a decoy ends with support in bin b, given >= 1 line
    dsup = decoys[decoys['n_decoy_lines'] >= 1].copy()
    dsup_bins = support_bin(dsup['n_decoy_lines'])
    q_bin = {b: (dsup_bins == b).sum() / len(dsup)
             for b in range(1, SUPPORT_MAX_BIN + 1)}

    # --- spurious model, part 2: |d| at given n (log-normal, fitted) ---
    dsup['abs_d'] = dsup['dE'].abs() * dsup['n_decoy_lines']
    a_s, b_s, s_s, n_s = fit_lognormal_d(dsup['n_decoy_lines'], dsup['abs_d'])

    # --- genuine model: |d| at given n, fitted to the old levels ---
    a_g, b_g, s_g, n_g = fit_lognormal_d(old_cal['n_tot'], old_cal['abs_d'],
                                         b_fixed=gen_b_fixed)

    # --- evidence factors of the tested levels ---
    m_bin = {b: int((support_bin(tested['n_tot']) == b).sum())
             for b in range(1, SUPPORT_MAX_BIN + 1)}
    n_arr = tested['n_tot'].to_numpy(dtype=float)
    bin_arr = np.asarray(support_bin(tested['n_tot']))
    fs_d = np.array([lognormal_density(d, n, a_s, b_s, s_s)
                     for d, n in zip(tested['abs_d'], n_arr)])
    fg_d = np.array([lognormal_density(d, n, a_g, b_g, s_g)
                     for d, n in zip(tested['abs_d'], n_arr)])

    def scan(fs, fg):
        """Maximum-likelihood S and posteriors for one evidence stream."""
        def loglik(S):
            ll = 0.0
            for i in range(len(tested)):
                b = int(bin_arr[i])
                pi = min(1.0, S * q_bin[b] / m_bin[b]) if m_bin[b] else 0.0
                ll += np.log(max(pi * fs[i] + (1.0 - pi) * fg[i], 1e-300))
            return ll

        S_grid = np.arange(0.0, len(tested) + 0.001, 0.25)
        ll_grid = np.array([loglik(S) for S in S_grid])
        i_best = int(np.argmax(ll_grid))
        S_best = float(S_grid[i_best])
        ok95 = S_grid[ll_grid >= ll_grid[i_best] - 1.92]
        p = np.empty(len(tested))
        for i in range(len(tested)):
            b = int(bin_arr[i])
            pi = min(1.0, S_best * q_bin[b] / m_bin[b]) if m_bin[b] else 0.0
            p[i] = pi * fs[i] / max(pi * fs[i] + (1.0 - pi) * fg[i], 1e-300)
        return p, S_best, (float(ok95.min()), float(ok95.max()))

    info = {'a_spur': a_s, 'b_spur': b_s, 's_spur': s_s, 'n_spur_fit': n_s,
            'a_gen': a_g, 'b_gen': b_g, 's_gen': s_g, 'n_gen_fit': n_g,
            'q_bin': q_bin, 'm_bin': m_bin}
    res = {'dE': scan(fs_d, fg_d)}

    if decoy_pattern is not None:
        # Pattern evidence: probability of the level's pattern_V range under
        # each hypothesis; levels without predictions carry no information
        pg = v_bin_probs(old_cal.loc[old_cal['n_pred'] > 0, 'pattern_V'])
        ps = v_bin_probs(decoy_pattern['V'])
        fs_p = np.ones(len(tested))
        fg_p = np.ones(len(tested))
        v_arr = tested['pattern_V'].to_numpy(dtype=float)
        npred_arr = tested['n_pred'].to_numpy(dtype=int)
        for i in range(len(tested)):
            if npred_arr[i] > 0 and np.isfinite(v_arr[i]):
                k = v_bin_index(v_arr[i])
                fs_p[i] = ps[k]
                fg_p[i] = pg[k]
        info['v_probs_gen'] = pg
        info['v_probs_spur'] = ps
        res['pattern'] = scan(fs_p, fg_p)
        res['joint'] = scan(fs_d * fs_p, fg_d * fg_p)

    res['info'] = info
    return res


def read_energy_overrides(path: str, levels: pd.DataFrame) -> pd.DataFrame:
    """Replace E_input for the levels listed in a two-column csv.

    The file has a header and the columns level_id and E_input.  It is needed
    when identifications were revised in a way that MOVES a level: dE tests
    whether the optimizer stays at the energy the identifications were built
    for, so for a level deliberately re-positioned the starting energy is the
    new hypothesised one, not the value in the adopted-level workbook.  Levels
    absent from the file keep their workbook energy.

    Once the same file is declared as files.level_overrides in the
    configuration, chance_mc.read_input_levels() has already applied it and
    --e-input only repeats the substitution (and prints it); it is then the
    way to test a position that is NOT the adopted one.
    """
    emap = cl.read_energy_overrides(path)
    known = set(levels['level_id'])
    unknown = sorted(set(emap) - known)
    if unknown:
        raise ValueError(f"{os.path.basename(path)}: level id(s) not in the "
                         f"level list: {', '.join(unknown)}")
    out = levels.copy()
    old = dict(zip(out['level_id'], out['E_input']))
    out['E_input'] = out['level_id'].map(lambda l: emap.get(l, old[l]))
    print(f"E_input overridden for {len(emap)} level(s) from "
          f"{os.path.basename(path)}:")
    for lid, e in sorted(emap.items()):
        print(f"  {lid}  {old[lid]:.4f} -> {e:.4f}  "
              f"({e - old[lid]:+.4f} cm^-1)")
    return out


def read_run():
    """The run the report describes: (levels, accepted-line table, obs_lines).

    Without --lopt this is the classify_lines.py run: its output table, whose
    low_E / upp_E columns hold the energies of its own internal least-squares
    optimization.  With --lopt the accepted set is read from a LOPT
    line-output file instead (see lopt_lines.py), which is how a hand-revised
    set of identifications enters the report.

    Which energies the LOPT path uses is set by ENERGIES (--energies):
      'refit' (the default) re-solves the level energies from the accepted
              lines with the pipeline's own weighted least squares.  This is
              what the report needs.  The energy shift dE of a tested level is
              judged against the dE of levels that are false by construction -
              the decoys of decoy_mc.py and the chance levels of chance_mc.py -
              and those are produced by re-running the pipeline, so their
              E_final comes from the pipeline's optimizer.  Measuring the
              tested levels with LOPT instead would compare the two sides with
              different instruments: LOPT's energies differ from the
              pipeline's by up to 0.2 cm^-1, which is the same size as the dE
              under test.  Fed the pipeline's own identifications the refit
              reproduces its energies to 1e-6 cm^-1 (median), so a --lopt
              report and the ordinary one differ only in the identifications.
      'lopt'  uses LOPT's own optimized energies, from --lopt-levels when a
              level-output file is given and from the E1/E2 columns of the
              line output otherwise.  Right for reporting the adopted
              energies; not comparable with the calibrations.

    obs_lines (wavenumber and intensity of every measured line) is used only by
    the masking check.  In LOPT mode it comes from the observed-line workbook,
    which is the complete list; the classify path keeps using the candidate
    table it always used, so the existing report is unchanged.
    """
    levels = mc.read_input_levels()
    if E_INPUT_CSV:
        levels = read_energy_overrides(E_INPUT_CSV, levels)
    if LOPT_LINES:
        print(f"Run read from LOPT output instead of "
              f"{os.path.basename(cl.OUTPUT_CSV)}:")
        real = lopt_lines.read_run(
            LOPT_LINES, energies=ENERGIES, levels_path=LOPT_LEVELS,
            e_input=dict(zip(levels['level_id'], levels['E_input'])))
        obs_lines = lopt_lines.read_observed_lines()[['wn_obs', 'obs_intens']]
    else:
        real = pd.read_csv(cl.OUTPUT_CSV, dtype={'low_id': str, 'upp_id': str})
        obs_lines = real.drop_duplicates('wn_obs')[['wn_obs', 'obs_intens']]
    return levels, real, obs_lines


def main():
    output_files.require_writable(output_files.with_twin(REPORT_CSV),
                                  'report file')
    # ------------------------------------------------------------------
    # Real run, per level
    # ------------------------------------------------------------------
    levels, real, obs_lines = read_run()
    per = mc.per_level_table(real, levels)
    per['n_tot'] = per['n_old'] + per['n_new']
    per['abs_dE'] = per['dE'].abs()
    per['d'] = per['dE'] * per['n_tot']
    per['abs_d'] = per['d'].abs()
    # new* = supported by more new than old identifications: the tested set
    per['new_star'] = (per['n_new'] > per['n_old']).astype(int)

    # Intensity-pattern scores of every level (see module docstring)
    e_final_map = dict(zip(per['level_id'], per['E_final']))
    preds = load_predictions(set(per['level_id']), e_final_map)
    calib = read_intensity_calibration()
    n_pred_all = len(preds)
    bias = intensity_scale_bias(real, preds, e_final_map)
    preds, preds_below, n_gap = filter_observable(preds, e_final_map, calib,
                                                  bias)
    real_pairs = accepted_pair_set(real)
    sc = pattern_scores(per['level_id'], real_pairs, preds)
    per['n_pred'] = [sc[l][0] for l in per['level_id']]
    per['n_top10'] = [sc[l][1] for l in per['level_id']]
    per['top10_found'] = [sc[l][2] for l in per['level_id']]
    per['pattern_C'] = [sc[l][3] for l in per['level_id']]
    per['pattern_V'] = [sc[l][4] for l in per['level_id']]
    # fate of each level's strongest predicted observable transition
    by_level = predictions_by_level(preds, e_final_map)
    per['top1'] = per['level_id'].map(
        top1_status(per['level_id'], by_level, real_pairs, obs_lines,
                    calib, bias))
    # how surely the level's branches would have been seen: the observation
    # probability of the strongest one (what makes top1 read 'faint') and the
    # expected number over all of them (what a missing-line term weighs the
    # accepted count against).  n_pred counts the same branches without asking
    # how visible they were, so n_obs_expected is always the smaller.
    obs = observability(per['level_id'], by_level, calib, bias)
    per['top1_pobs'] = [obs[l][0] for l in per['level_id']]
    per['n_obs_expected'] = [obs[l][1] for l in per['level_id']]
    # accepted lines whose predicted intensity is below Sugar's noise level
    # (kept in the support count n; this column only makes them visible)
    below_cnt = {}
    for lo, up, ip in preds_below.itertuples(index=False):
        if tuple(sorted((lo, up))) in real_pairs:
            below_cnt[lo] = below_cnt.get(lo, 0) + 1
            below_cnt[up] = below_cnt.get(up, 0) + 1
    per['n_acc_below_noise'] = per['level_id'].map(lambda l: below_cnt.get(l, 0))

    n_tested = int(per['new_star'].sum())
    print(f"Levels read: {len(per)} total; tested (new*, n_new > n_old): {n_tested}; "
          f"old (reference): {len(per) - n_tested}")
    if calib is not None:
        print(f"Predicted transitions loaded: {n_pred_all}; observable: "
              f"{len(preds)} (excluded: {n_gap} in blind stretches, "
              f"{len(preds_below)} below Sugar's noise level)")
    else:
        print(f"Predicted transitions loaded: {n_pred_all}; observable: "
              f"{len(preds)} ({n_gap} in blind stretches excluded; "
              f"intensity_correction_functions.txt not found - no noise-level "
              f"cut applied)")
    tab = per.groupby(['is_new_level', 'new_star']).size()
    print("cross-tab vs ASD status (is_new_level = 1 means absent from ASD):")
    for (asd, star), cnt in tab.items():
        print(f"  absent_from_ASD={asd}, new_star={star}: {cnt}")

    old = per[per['new_star'] == 0]
    new = per[per['new_star'] == 1]

    # ------------------------------------------------------------------
    # Old-level calibration: legitimate |dE| range and the d statistic
    # ------------------------------------------------------------------
    cal = old[old['n_tot'] >= 1]
    a = cal['abs_dE'].to_numpy()
    dsig = cal['d'].to_numpy()
    print("\n" + "=" * 60)
    print(f"OLD-LEVEL CALIBRATION ({len(cal)} old levels with >= 1 accepted line)")
    print("=" * 60)
    pcts = {'p50': q(a, 50), 'p90': q(a, 90), 'p95': q(a, 95),
            'p99': q(a, 99), 'max': float(a.max()) if a.size else float('nan')}
    print("|dE| percentiles (pNN = value below which NN% of the levels fall):")
    print("  " + ",  ".join(f"{k} = {v:.4f}" for k, v in pcts.items()) + "  cm^-1")
    Ts = [(k, pcts[k]) for k in ('p90', 'p95', 'p99', 'max')]
    print(f"d = dE*n: aggregate rms = {np.sqrt((dsig ** 2).mean()):.3f} cm^-1 "
          f"(dominated by the many-line levels);")
    print("per given n the genuine shifts are much tighter:")
    print(f"  {'support bin':>12} {'count':>6} {'median |dE|':>12} {'median |d|':>11}")
    for lo, hi in ((1, 2), (3, 4), (5, 7), (8, 11), (12, 17), (18, 25), (26, 10 ** 6)):
        sbin = cal[(cal['n_tot'] >= lo) & (cal['n_tot'] <= hi)]
        if len(sbin):
            label = f"{lo}-{hi}" if hi < 10 ** 6 else f"{lo}+"
            print(f"  {label:>12} {len(sbin):6d} {sbin['abs_dE'].median():12.4f} "
                  f"{sbin['abs_d'].median():11.3f}")

    # ------------------------------------------------------------------
    # Tested (new*) levels in the real run
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print(f"TESTED (NEW*) LEVELS ({len(new)}) IN THE REAL RUN")
    print("=" * 60)
    counts = new['n_tot'].clip(upper=SUPPORT_MAX_BIN).value_counts().sort_index()
    dist = '  '.join(f"{('%d+' % n) if n == SUPPORT_MAX_BIN else n}: {c}"
                     for n, c in counts.items())
    print(f"support distribution (n = accepted lines):  {dist}")
    print("|dE|: " + ",  ".join(f"p{p} = {q(new['abs_dE'], p):.4f}"
                                for p in (50, 90, 95)) +
          f",  max = {new['abs_dE'].max():.4f}  cm^-1")
    for name, t in Ts:
        k = int((new['abs_dE'] <= t).sum())
        print(f"  within old-level {name} ({t:.4f} cm^-1): {k} of {len(new)}")

    worst = new.sort_values('abs_d', ascending=False).head(N_WORST)
    print(f"\nLargest |d| = |dE*n| among tested levels (top {len(worst)}):")
    print(f"  {'level_id':>16} {'J':>4} par {'E_input':>11} {'dE':>9} {'n':>4} {'d':>7}")
    for _, r in worst.iterrows():
        print(f"  {r['level_id']:>16} {r['J']:>4} {r['parity']:>3} "
              f"{r['E_input']:11.2f} {r['dE']:+9.4f} {r['n_tot']:4d} {r['d']:+7.2f}")

    # ------------------------------------------------------------------
    # Primary calibration: in-situ decoys
    # ------------------------------------------------------------------
    p_spur = None
    p_de_arr = None
    p_pat_arr = None
    qstat = {}
    if os.path.exists(dmc.LEVELS_CSV):
        dc = pd.read_csv(dmc.LEVELS_CSV, dtype={'level_id': str})
        dc['abs_dE'] = dc['dE'].abs()
        n_runs_d = int(dc['delta'].nunique())
        n_planted = int(dc.groupby('delta').size().iloc[0])
        print("\n" + "=" * 60)
        print(f"IN-SITU DECOY CALIBRATION ({n_runs_d} runs x {n_planted} decoys) "
              f"- primary")
        print("=" * 60)
        if n_planted != n_tested:
            print(f"WARNING: {n_planted} decoys per run vs {n_tested} tested levels "
                  f"- rerun decoy_mc.py after changing the tested set.")
        dE_by_n_table(new, dc, 'n_decoy_lines', 'decoys')
        criterion_grid(new, dc, 'n_decoy_lines', n_runs_d, Ts, n_tested,
                       'in-situ decoys, worst case')

        # ------------------------------------------------------------------
        # Pattern scores of the decoys (evidence for the spurious hypothesis)
        # ------------------------------------------------------------------
        dpat = None
        if os.path.exists(dmc.LINES_CSV):
            dl = pd.read_csv(dmc.LINES_CSV, dtype={'low_id': str, 'upp_id': str})
            drows = []
            for delta, grp in dl.groupby('delta'):
                pairs_d = accepted_pair_set(grp, strip_decoy_prefix=True)
                scd = pattern_scores(new['level_id'], pairs_d, preds)
                nsup = dc[dc['delta'] == delta].set_index('level_id')['n_decoy_lines']
                for lid in new['level_id']:
                    if int(nsup.get(cl.DECOY_PREFIX + lid, 0)) >= 1:
                        n_pred_d, n_topk, found, c_val, v_val = scd[lid]
                        if n_pred_d > 0:
                            drows.append((delta, lid, n_topk, found, c_val, v_val))
            dpat = pd.DataFrame(drows, columns=['delta', 'level_id', 'n_topk',
                                                'top_found', 'C', 'V'])
        else:
            print("\ndecoy_mc_lines.csv not found - pattern evidence "
                  "unavailable; probabilities use the energy shifts only.")

        # ------------------------------------------------------------------
        # Spurious probabilities: energy-shift, pattern, and folded
        # ------------------------------------------------------------------
        res = fit_spurious_probabilities(new, dc, cal, decoy_pattern=dpat)
        info = res['info']
        p_de_arr, S_de, (S_de_lo, S_de_hi) = res['dE']
        # Robustness of the genuine-model extrapolation to small n, where
        # old calibration levels are scarce: repeat with the slope of the
        # genuine ln|d| vs ln n line held at 1 (|dE| independent of n, as
        # the per-bin medians of the old levels suggest).
        res_alt = fit_spurious_probabilities(new, dc, cal, decoy_pattern=dpat,
                                             gen_b_fixed=1.0)
        _, S_alt, (S_alt_lo, S_alt_hi) = res_alt['dE']
        if dpat is not None:
            p_pat_arr, S_pat, (S_pat_lo, S_pat_hi) = res['pattern']
            p_spur, S_joint, (S_j_lo, S_j_hi) = res['joint']
            _, S_jalt, (S_jalt_lo, S_jalt_hi) = res_alt['joint']
        else:
            p_pat_arr = None
            p_spur, S_joint, (S_j_lo, S_j_hi) = p_de_arr, S_de, (S_de_lo, S_de_hi)
        newp = new.assign(
            p_spur=p_spur, p_spur_decoys=p_de_arr,
            p_spur_pattern=(p_pat_arr if p_pat_arr is not None else np.nan))

        # ------------------------------------------------------------------
        # Intensity-pattern check (the square-array argument)
        # ------------------------------------------------------------------
        if dpat is not None:
            print("\n" + "=" * 60)
            print("INTENSITY-PATTERN CHECK")
            print("=" * 60)
            print("top10 = share of the 10 strongest predicted transitions of the "
                  "level found among\nthe accepted lines; C = intensity-weighted "
                  "completeness; V = C relative to the\nbest C achievable with the "
                  "level's matched lines. Tested levels are split by the\n"
                  "energy-shift evidence alone to show that the pattern is an "
                  "independent check.")

            def pat_line(label, top_found, n_topk, c_vals, v_vals):
                frac = np.asarray(top_found, dtype=float) / np.asarray(n_topk, dtype=float)
                c_vals = np.asarray(c_vals, dtype=float)
                v_vals = np.asarray(v_vals, dtype=float)
                print(f"  {label:<38} {len(c_vals):5d}  {np.mean(frac):6.2f}  "
                      f"{q(c_vals, 50):6.2f}  {q(v_vals, 25):6.2f} "
                      f"{q(v_vals, 50):6.2f} {q(v_vals, 75):6.2f}")

            print(f"\n  {'group':<38} {'N':>5}  {'top10':>6}  {'med C':>6}  "
                  f"{'V: p25':>6} {'p50':>6} {'p75':>6}")
            gsup = cal[cal['n_pred'] > 0]
            pat_line("old levels (genuine)", gsup['top10_found'], gsup['n_top10'],
                     gsup['pattern_C'], gsup['pattern_V'])
            t_lo = newp[(newp['p_spur_decoys'] < P_QUESTIONABLE) & (newp['n_pred'] > 0)]
            pat_line(f"tested, p_spur_decoys < {P_QUESTIONABLE}", t_lo['top10_found'],
                     t_lo['n_top10'], t_lo['pattern_C'], t_lo['pattern_V'])
            t_hi = newp[(newp['p_spur_decoys'] >= P_QUESTIONABLE) & (newp['n_pred'] > 0)]
            pat_line(f"tested, p_spur_decoys >= {P_QUESTIONABLE}", t_hi['top10_found'],
                     t_hi['n_top10'], t_hi['pattern_C'], t_hi['pattern_V'])
            pat_line("decoys with >= 1 line (fake)", dpat['top_found'],
                     dpat['n_topk'], dpat['C'], dpat['V'])
            good = dpat[dpat['top_found'] >= 3]
            print(f"\nDecoy trials with >= 3 of the top-{K_TOP} predictions accepted: "
                  f"{len(good)} of {len(dpat)} ({100 * len(good) / len(dpat):.1f}%)"
                  f" - an upper bound on the chance\nthat a displaced level "
                  f"accidentally sits on a real, previously unknown level.")
            print("\npattern_V range probabilities used as evidence "
                  "(genuine / spurious):")
            print("  " + "  ".join(f"{lab}: {g:.2f}/{s:.2f}"
                                   for lab, g, s in zip(V_LABELS,
                                                        info['v_probs_gen'],
                                                        info['v_probs_spur'])))

        # ------------------------------------------------------------------
        # Probability results
        # ------------------------------------------------------------------
        print("\n" + "=" * 60)
        print("PER-LEVEL SPURIOUS PROBABILITIES")
        print("=" * 60)
        print("log-normal |d| models (ln|d| = a + b*ln n, scatter s):")
        print(f"  spurious (decoys):    a = {info['a_spur']:+.3f}, "
              f"b = {info['b_spur']:+.3f}, s = {info['s_spur']:.3f}   "
              f"({info['n_spur_fit']} supported decoys)")
        print(f"  genuine (old levels): a = {info['a_gen']:+.3f}, "
              f"b = {info['b_gen']:+.3f}, s = {info['s_gen']:.3f}   "
              f"({info['n_gen_fit']} old levels)")
        print(f"\nMost probable number of spurious levels among the {n_tested} tested"
              " (with 95% confidence ranges):")
        print(f"  energy shifts alone:        S = {S_de:5.1f}  [{S_de_lo:.1f} .. {S_de_hi:.1f}]")
        print(f"    robustness check (genuine |dE| independent of n): "
              f"S = {S_alt:.1f}  [{S_alt_lo:.1f} .. {S_alt_hi:.1f}]")
        if dpat is not None:
            print(f"  intensity pattern alone:    S = {S_pat:5.1f}  [{S_pat_lo:.1f} .. {S_pat_hi:.1f}]")
            print(f"  both folded (final p_spur): S = {S_joint:5.1f}  [{S_j_lo:.1f} .. {S_j_hi:.1f}]")
            print(f"    robustness check, folded: S = {S_jalt:.1f}  "
                  f"[{S_jalt_lo:.1f} .. {S_jalt_hi:.1f}]")
        print(f"consistency: sum of folded p_spur = {np.asarray(p_spur).sum():.1f}")
        edges = [0.05, 0.5, 0.9]
        c0 = int((p_spur < edges[0]).sum())
        c1 = int(((p_spur >= edges[0]) & (p_spur < edges[1])).sum())
        c2 = int(((p_spur >= edges[1]) & (p_spur < edges[2])).sum())
        c3 = int((p_spur >= edges[2]).sum())
        print(f"levels by folded p_spur: <5%: {c0},  5-50%: {c1},  "
              f"50-90%: {c2},  >90%: {c3}")

        quest = newp[newp['p_spur'] >= P_QUESTIONABLE].sort_values(
            'p_spur', ascending=False)
        print(f"\nTested levels with folded p_spur >= {P_QUESTIONABLE} "
              f"({len(quest)}):")
        print(f"  {'level_id':>16} {'J':>4} par {'E_input':>11} {'dE':>9} "
              f"{'n':>3} {'top10':>6} {'V':>5} {'p_dE':>6} {'p_pat':>6} {'p_spur':>7}")
        for _, r in quest.iterrows():
            top = f"{int(r['top10_found'])}/{int(r['n_top10'])}"
            vval = f"{r['pattern_V']:.2f}" if np.isfinite(r['pattern_V']) else '  --'
            ppat = (f"{r['p_spur_pattern']:6.3f}"
                    if np.isfinite(r['p_spur_pattern']) else '    --')
            print(f"  {r['level_id']:>16} {r['J']:>4} {r['parity']:>3} "
                  f"{r['E_input']:11.2f} {r['dE']:+9.4f} {r['n_tot']:3d} "
                  f"{top:>6} {vval:>5} {r['p_spur_decoys']:6.3f} {ppat} "
                  f"{r['p_spur']:7.3f}")

        # ------------------------------------------------------------------
        # Masking check: adjudication of the questionable marks
        # ------------------------------------------------------------------
        qstat = question_status(newp, by_level, real_pairs, obs_lines,
                                dict(zip(per['level_id'], per['top1'])))
        t1 = per.loc[per['new_star'] == 1, 'top1'].value_counts()
        print(f"\nStrongest predicted branch of the {n_tested} tested levels: "
              f"{int(t1.get('found', 0))} accepted, "
              f"{int(t1.get('masked', 0))} absent but masked\nby a stronger "
              f"line, {int(t1.get('faint', 0))} absent but too faint on the "
              f"plate to have been\nrecorded, {int(t1.get('missing', 0))} "
              f"absent with nothing hiding them, {int(t1.get('', 0))} without "
              f"an observable prediction.")
        print(f"\nMasking check of the levels the probabilities question "
              f"(p_spur >= {P_QUESTIONABLE}) or whose strongest\npredicted "
              f"branch is missing (a missing predicted transition counts as "
              f"hidden if, within\none effective line width - instrumental "
              f"FWHM {INSTR_FWHM_A:g} angstrom convolved with the Doppler "
              f"width -\nan observed line reaches "
              f"{MASK_STRENGTH:g} x its predicted intensity, or if the "
              f"Gaussian wing of a stronger\nline still exceeds the predicted "
              f"intensity at its position; the mark is cleared when >= "
              f"\n{MASK_CLEAR_FRACTION:.0%} of the missing predicted intensity "
              f"is hidden, the energy shift alone is not suspect,\nand the "
              f"strongest branch is not among the unexplained absences):")
        removed = sorted((lid, f) for lid, (st, _, f) in qstat.items()
                         if st == 'removed')
        print(f"  question_status = removed (questionable mark cleared): "
              f"{len(removed)}")
        for lid, f in removed:
            print(f"    {lid}  ({f:.0%} of the missing predicted intensity is masked)")
        reasons = {}
        for lid, (st, why, _) in qstat.items():
            if st == 'retained':
                reasons[why] = reasons.get(why, 0) + 1
        print(f"  question_status = retained: {sum(reasons.values())}")
        for why, cnt in sorted(reasons.items(), key=lambda kv: -kv[1]):
            print(f"    {cnt:3d}  {why}")
    else:
        print("\ndecoy_mc_levels.csv not found - run decoy_mc.py for the "
              "primary calibration and the spurious probabilities.")

    # ------------------------------------------------------------------
    # Independent cross-check: shifted-wavenumber chance runs (upper bound)
    # ------------------------------------------------------------------
    if os.path.exists(mc.LEVELS_CSV):
        ch = pd.read_csv(mc.LEVELS_CSV, dtype={'level_id': str})
        ch['abs_dE'] = ch['dE'].abs()
        n_runs_c = int(ch['shift'].nunique())
        star_set = set(new['level_id'])
        chn = ch[ch['level_id'].isin(star_set)]
        print("\n" + "=" * 60)
        print(f"SHIFTED-WAVENUMBER CHANCE CALIBRATION ({n_runs_c} runs) - upper bound")
        print("=" * 60)
        print("Counts restricted to the tested levels; the whole line list is "
              "available for false matching\nhere (nothing is consumed by true "
              "identifications), so these false-pass counts are upper bounds.")
        dE_by_n_table(new, chn, 'n_chance_lines', 'chance')
        criterion_grid(new, chn, 'n_chance_lines', n_runs_c, Ts, n_tested,
                       'shifted-wavenumber chance, worst case')
    else:
        print("\nchance_mc_levels.csv not found - run chance_mc.py for the "
              "upper-bound cross-check.")

    # ------------------------------------------------------------------
    # Report file
    # ------------------------------------------------------------------
    out = per.copy()
    out['p_spur_decoys'] = np.nan
    out['p_spur_pattern'] = np.nan
    out['p_spur'] = np.nan
    if p_spur is not None:
        out.loc[new.index, 'p_spur'] = p_spur
        out.loc[new.index, 'p_spur_decoys'] = p_de_arr
        if p_pat_arr is not None:
            out.loc[new.index, 'p_spur_pattern'] = p_pat_arr
    out['question_status'] = out['level_id'].map(
        lambda lid: qstat.get(lid, ('', ''))[0])
    out['reason'] = out['level_id'].map(
        lambda lid: qstat.get(lid, ('', ''))[1])
    cols = ['level_id', 'is_new_level', 'new_star', 'note', 'J', 'parity',
            'E_input', 'E_final', 'dE', 'n_old', 'n_new', 'n_tot', 'd',
            'n_pred', 'n_obs_expected', 'n_top10', 'top10_found', 'top1',
            'top1_pobs', 'pattern_C', 'pattern_V',
            'n_acc_below_noise', 'p_spur_decoys', 'p_spur_pattern', 'p_spur',
            'question_status', 'reason']
    print()
    mc.save_table(out[cols], REPORT_CSV,
                  decimals={'E_input': 4, 'E_final': 4, 'dE': 4, 'd': 3,
                            'n_obs_expected': 2, 'top1_pobs': 3,
                            'pattern_C': 3, 'pattern_V': 3, 'p_spur_decoys': 3,
                            'p_spur_pattern': 3, 'p_spur': 3})


def doppler_fwhm_factor() -> float:
    """Doppler full width at half maximum divided by the wavelength
    (dimensionless): sqrt(8*ln2*kT/m)/c, with the pipeline's plasma
    temperature and atomic mass (8.22e-6 for Pr at T = 1.6 eV)."""
    eV_to_J = 1.602176634e-19
    u_to_kg = 1.66053906660e-27
    c = 2.99792458e8  # m/s
    return float(np.sqrt(8.0 * np.log(2.0) * cl.T * eV_to_J
                         / (cl.ATOMIC_MASS * u_to_kg)) / c)


def _hidden_by_a_stronger_line(ritz: float, i_pred: float, wn_arr, int_arr,
                               dop: float) -> bool:
    """Is a predicted transition hidden under an observed line?

    ritz is the predicted position in cm^-1 and i_pred its predicted
    intensity; wn_arr and int_arr are the positions and intensities of all
    observed lines, wn_arr sorted; dop is the Doppler width factor returned by
    doppler_fwhm_factor().

    The test uses the effective line width (full width at half maximum): the
    instrumental width INSTR_FWHM_A convolved with the Doppler width, both in
    wavelength units, where the spectrograph resolution is nearly constant.
    An observed line hides the prediction in either of two regimes:
      - within one effective width the two lines are not resolved, so a line
        of at least MASK_STRENGTH times the predicted intensity absorbs it;
      - beyond one width the prediction is hidden under the WING of a much
        stronger line: a Gaussian profile of the effective width is applied to
        the observed line, and the prediction counts as hidden if the profile
        intensity at the predicted position still exceeds i_pred.  The masking
        distance therefore grows with the intensity ratio.
    """
    lam = 1.0e8 / ritz                      # angstrom
    fwhm_wn = float(np.hypot(INSTR_FWHM_A, dop * lam)) * ritz * ritz * 1.0e-8
    # candidates out to 3 widths (a Gaussian wing at 3 FWHM requires an
    # intensity ratio above e^25 - far beyond any real line pair)
    win = 3.0 * fwhm_wn
    i0 = np.searchsorted(wn_arr, ritz - win)
    i1 = np.searchsorted(wn_arr, ritz + win)
    four_ln2 = 4.0 * np.log(2.0)
    for j in range(i0, i1):
        sep = abs(wn_arr[j] - ritz)
        if sep <= fwhm_wn:
            if int_arr[j] >= MASK_STRENGTH * i_pred:     # unresolved
                return True
        elif (int_arr[j] * np.exp(-four_ln2 * (sep / fwhm_wn) ** 2)) >= i_pred:
            return True
    return False


def predictions_by_level(preds: pd.DataFrame, e_final_map: dict) -> dict:
    """{level_id: [(I_pred, pair, ritz), ...]} sorted strongest first.

    pair is the unordered level-id pair of the predicted transition and ritz
    its Ritz position, E(upper) - E(lower), from the optimized energies.
    """
    by_level = {}
    for lo, up, ip in preds.itertuples(index=False):
        pair = tuple(sorted((lo, up)))
        ritz = e_final_map[up] - e_final_map[lo]
        by_level.setdefault(lo, []).append((ip, pair, ritz))
        by_level.setdefault(up, []).append((ip, pair, ritz))
    for lid in by_level:
        by_level[lid].sort(key=lambda t: -t[0])
    return by_level


def observability(level_ids, by_level: dict, calib=None,
                  bias: dict | None = None) -> dict:
    """{level_id: (P_obs of the strongest prediction, sum of P_obs)}.

    The second number is how many of the level's predicted branches should
    have appeared in the line list at all - the expected count, against which
    the number actually accepted is what a missing-line term has to weigh.  It
    is well below n_pred for most levels, because most predicted branches sit
    near the noise: of the 7391 predictions that survive the observability
    filter of the run, the great majority are between one and ten times the
    noise level, where the measured detection curve stands at 0.4 to 0.8.

    The first number says why a level's strongest branch was called 'faint':
    it is the quantity compared with DETECT_CONFIDENCE in top1_status.
    """
    out = {}
    for lid in level_ids:
        plist = by_level.get(lid, [])
        if not plist:
            out[lid] = (np.nan, 0.0)
            continue
        po = [observation_probability(ritz, ip, calib, bias)
              for ip, _, ritz in plist]
        out[lid] = (po[0], float(sum(po)))
    return out


def top1_status(level_ids, by_level: dict, accepted_pairs: set,
                obs_lines: pd.DataFrame, calib=None,
                bias: dict | None = None) -> dict:
    """Fate of each level's STRONGEST predicted observable transition.

    'found'   - it is among the accepted identifications;
    'masked'  - it is not, but an observed line strong enough and close enough
                to hide it stands at its predicted position, so its absence
                from the line list is explained;
    'faint'   - it is not, nothing hides it, but the line would not reliably
                have been observed in the first place: the measured chance of
                seeing it, the coverage of the plate there times the measured
                detection curve at its intensity, is below DETECT_CONFIDENCE
                (observation_probability).  That covers both the far
                ultraviolet, where the calculated intensities are an order of
                magnitude too large and the line is at the noise, and the
                half-blind stretches, where the plate is thinly exposed;
    'missing' - it is not, nothing hides it, and it should have been
                recorded: the branch the theory makes the level's brightest
                is simply not there;
    ''        - the level has no observable prediction, so the test does not
                apply.

    Why this is worth a column of its own.  A level can score well on both the
    energy shift and the intensity pattern and still be wrong, if what it
    explains is a set of weaker lines while its brightest predicted branch is
    absent for no reason.  pattern_V does not catch that on its own: V measures
    whether the ACCEPTED lines are the strongest of the predictions, so a level
    whose two accepted lines happen to be ranked 2 and 3 can reach V near 1
    while rank 1 is missing.  A missing, unmasked strongest branch is the
    single most direct piece of evidence against a level, and it keeps the
    questionable mark whatever the probabilities say.
    """
    obs_sorted = obs_lines.sort_values('wn_obs')
    wn_arr = obs_sorted['wn_obs'].to_numpy(dtype=float)
    int_arr = obs_sorted['obs_intens'].to_numpy(dtype=float)
    dop = doppler_fwhm_factor()
    out = {}
    for lid in level_ids:
        plist = by_level.get(lid, [])
        if not plist:
            out[lid] = ''
            continue
        ip, pair, ritz = plist[0]
        if pair in accepted_pairs:
            out[lid] = 'found'
        elif _hidden_by_a_stronger_line(ritz, ip, wn_arr, int_arr, dop):
            out[lid] = 'masked'
        elif observation_probability(ritz, ip, calib, bias) < DETECT_CONFIDENCE:
            out[lid] = 'faint'
        else:
            out[lid] = 'missing'
    return out


def question_status(newp: pd.DataFrame, by_level: dict, accepted_pairs: set,
                    obs_lines: pd.DataFrame, top1: dict) -> dict:
    """Adjudication of the questionable levels: can the bad score be excused?

    A level is examined here when either the probabilities call it into
    question - folded p_spur >= P_QUESTIONABLE - or its strongest predicted
    branch is missing and unmasked (top1_status == 'missing'), which is
    evidence against the level no matter what the probabilities say.  A
    strongest branch that is absent but too faint to have been recorded
    (top1_status == 'faint') is not evidence and does not put the level up
    for examination on its own.

    The examination is the one an expert would make: a strong predicted
    transition may be absent from the accepted lines simply because a nearby
    stronger line hides it on the photographic plates
    (_hidden_by_a_stronger_line).  If at least MASK_CLEAR_FRACTION of the
    missing predicted intensity among the K_TOP strongest predictions is
    hidden - i.e. the DOMINANT missing predictions are masked, a small
    unmasked residue of weak predictions being tolerated because the theory is
    least reliable for weak transitions - then the bad pattern score carries
    no evidence against the level and the questionable mark is cleared
    (question_status = 'removed').  Clearing requires two further things: the
    energy-shift evidence alone must not remain suspect
    (p_spur_decoys < P_QUESTIONABLE), and the strongest predicted branch must
    not be missing and unmasked.

    Returns {level_id: (question_status, reason, masked_fraction)} for the
    levels examined; masked_fraction is NaN when the test does not apply.
    """
    obs_sorted = obs_lines.sort_values('wn_obs')
    wn_arr = obs_sorted['wn_obs'].to_numpy(dtype=float)
    int_arr = obs_sorted['obs_intens'].to_numpy(dtype=float)
    dop = doppler_fwhm_factor()

    out = {}
    for _, r in newp.iterrows():
        lid = r['level_id']
        suspect = bool(r['p_spur'] >= P_QUESTIONABLE)
        no_top1 = top1.get(lid) == 'missing'
        if not (suspect or no_top1):
            continue
        plist = by_level.get(lid, [])
        if not plist:
            out[lid] = ('retained',
                        'no theoretical predictions (pattern check not applicable)',
                        np.nan)
            continue
        top = plist[:min(K_TOP, len(plist))]
        missing = [(ip, pair, ritz) for ip, pair, ritz in top
                   if pair not in accepted_pairs]
        if not missing:
            out[lid] = ('retained', 'no reason to clear', np.nan)
            continue
        i_missing = sum(ip for ip, _, _ in missing)
        i_masked = sum(ip for ip, _, ritz in missing
                       if _hidden_by_a_stronger_line(ritz, ip, wn_arr, int_arr, dop))
        frac = i_masked / i_missing
        if no_top1:
            out[lid] = ('retained', 'the strongest predicted transition of the '
                                    'level is absent and nothing hides it', frac)
        elif frac < MASK_CLEAR_FRACTION:
            out[lid] = ('retained', 'no reason to clear', frac)
        elif r['p_spur_decoys'] >= P_QUESTIONABLE:
            out[lid] = ('retained', 'masking found, but the energy shift '
                                    'alone remains suspect', frac)
        else:
            out[lid] = ('removed', 'strong missing transitions are masked '
                                   'by nearby stronger lines', frac)
    return out


def _opt(value, prefix: str = '', width: int = 0) -> str:
    """`prefix + value`, or nothing when the value is missing.

    In --lopt mode the grade and the weeding note of a candidate do not exist
    (LOPT records neither), and printing them as "nan" only clutters the line.
    """
    text = '' if value is None else str(value).strip()
    if text in ('', 'nan', 'None'):
        return ''
    return prefix + (text[:width] if width else text)


def print_level_detail(level_id: str):
    """Inspection aid: one level's predicted spectrum against the real run.

    Lists every predicted observable transition of the level (strongest
    first) with the partner level, the Ritz wavenumber, the predicted
    intensity, and its fate in the real run: ACCEPTED (with the observed
    line), REJECTED (with the weeding note), or NO MATCH (no observed line
    was matched to it). Accepted lines without a theoretical prediction are
    listed at the end.
    """
    levels, real, _ = read_run()
    per = mc.per_level_table(real, levels)
    row = per[per['level_id'] == level_id]
    if row.empty:
        print(f"level {level_id} not found in the level list")
        return
    row = row.iloc[0]
    print(f"Level {level_id}   J={row['J']}  parity={row['parity']}  "
          f"E_input={row['E_input']:.2f}  E_final={row['E_final']:.4f}  "
          f"dE={row['dE']:+.4f}  n_old={int(row['n_old'])}  n_new={int(row['n_new'])}")
    if os.path.exists(REPORT_CSV):
        rep = pd.read_csv(REPORT_CSV, dtype={'level_id': str})
        rr = rep[rep['level_id'] == level_id]
        if len(rr) and np.isfinite(rr.iloc[0]['p_spur']):
            rr = rr.iloc[0]
            print(f"p_spur_decoys={rr['p_spur_decoys']:.3f}  "
                  f"p_spur_pattern={rr['p_spur_pattern']:.3f}  "
                  f"p_spur={rr['p_spur']:.3f}")

    e_final_map = dict(zip(per['level_id'], per['E_final']))
    preds = load_predictions(set(per['level_id']), e_final_map)
    calib = read_intensity_calibration()
    bias = intensity_scale_bias(real, preds, e_final_map, verbose=False)
    mine = preds[(preds['lo_id'] == level_id) |
                 (preds['up_id'] == level_id)].sort_values('I_pred',
                                                           ascending=False)
    d = real.copy()
    for c in ('low_id', 'upp_id'):
        d[c] = d[c].fillna('').astype(str)
    cand = d[(d['low_id'] == level_id) | (d['upp_id'] == level_id)]
    by_pair = {}
    for _, cr in cand.iterrows():
        by_pair.setdefault(tuple(sorted((cr['low_id'], cr['upp_id']))),
                           []).append(cr)

    print(f"\n{len(mine)} predicted observable transitions (strongest first):")
    print(f"  {'partner':>16} {'Ritz wn':>11} {'I_pred':>10}  fate in the run")
    for _, p in mine.iterrows():
        partner = p['up_id'] if p['lo_id'] == level_id else p['lo_id']
        wn_ritz = e_final_map[p['up_id']] - e_final_map[p['lo_id']]
        rows = by_pair.get(tuple(sorted((p['lo_id'], p['up_id']))), [])
        if not rows:
            fate = "NO MATCH (no observed line matched)"
        else:
            cr = next((r for r in rows if r['accepted'] == 1), rows[0])
            base = (f"wn_obs={cr['wn_obs']:.3f}  I_obs={cr['obs_intens']:.0f}"
                    + _opt(cr['grade'], '  grade='))
            if cr['accepted'] == 1:
                fate = "ACCEPTED  " + base
            else:
                fate = "REJECTED  " + base + _opt(cr['notes2'], '  ', 55)
        cov = coverage(wn_ritz)
        if cov < COVERAGE_MIN:
            fate += f"  [blind stretch, coverage {cov:.2f} - not observable]"
        elif calib is not None:
            thr = noise_threshold_linear(wn_ritz, calib)
            f = scale_bias_factor(wn_ritz, bias)
            pobs = observation_probability(wn_ritz, p['I_pred'], calib, bias)
            if thr is not None and p['I_pred'] * f < thr:
                fate += ("  [below noise level]" if f == 1.0 else
                         f"  [below noise level on the plate: I_pred x {f:.2f}"
                         f" = {p['I_pred'] * f:.0f} vs {thr:.0f}]")
            elif pobs < DETECT_CONFIDENCE:
                # printed for accepted lines too - that a line which WAS seen
                # was only marginally likely to be is worth knowing - so the
                # wording has to read as a probability, not as a verdict
                fate += (f"  [P(seen) = {pobs:.2f}: "
                         + (f"I_pred x {f:.2f} = {p['I_pred'] * f:.0f}"
                            if f != 1.0 else f"I_pred = {p['I_pred']:.0f}")
                         + f" vs noise {thr:.0f}"
                         + (f", coverage {cov:.2f}" if cov < 0.99 else "")
                         + "]")
        print(f"  {partner:>16} {wn_ritz:11.3f} {p['I_pred']:10.1f}  {fate}")

    pred_pairs = {tuple(sorted((p['lo_id'], p['up_id'])))
                  for _, p in mine.iterrows()}
    extra = [(pr, next(r for r in rows if r['accepted'] == 1))
             for pr, rows in by_pair.items()
             if pr not in pred_pairs and any(r['accepted'] == 1 for r in rows)]
    if extra:
        print(f"\naccepted lines without a theoretical prediction ({len(extra)}):")
        for pr, cr in extra:
            partner = pr[0] if pr[1] == level_id else pr[1]
            print(f"  {partner:>16}  wn_obs={cr['wn_obs']:.3f}  "
                  f"I_obs={cr['obs_intens']:.0f}"
                  + _opt(cr['grade'], '  grade='))


def _take_option(argv: list, name: str):
    """Remove `--name VALUE` from argv and return VALUE (None if absent)."""
    if name not in argv:
        return None
    i = argv.index(name)
    if i + 1 >= len(argv):
        raise SystemExit(f"{name} needs a file name")
    value = argv[i + 1]
    del argv[i:i + 2]
    return value


def _default_report_name(lopt_path: str) -> str:
    """level_shift_report<tag>.csv, where <tag> follows LOPT_output_lines.

    LOPT_output_lines_revised.txt -> level_shift_report_revised.csv, so that a
    LOPT-driven report never silently overwrites the pipeline's own.
    """
    base = os.path.splitext(os.path.basename(lopt_path))[0]
    stem = 'LOPT_output_lines'
    tag = base[len(stem):] if base.startswith(stem) else '_' + base
    return os.path.join(HERE, f"level_shift_report{tag}.csv")


if __name__ == '__main__':
    _argv = sys.argv[1:]
    LOPT_LINES = _take_option(_argv, '--lopt')
    LOPT_LEVELS = _take_option(_argv, '--lopt-levels')
    ENERGIES = _take_option(_argv, '--energies') or 'refit'
    if ENERGIES not in ('refit', 'lopt'):
        raise SystemExit("--energies takes 'refit' or 'lopt'")
    if LOPT_LEVELS and ENERGIES != 'lopt':
        raise SystemExit("--lopt-levels is only used with --energies lopt")
    E_INPUT_CSV = _take_option(_argv, '--e-input')
    _report = _take_option(_argv, '--report')
    if LOPT_LINES:
        REPORT_CSV = _report or _default_report_name(LOPT_LINES)
    elif _report:
        REPORT_CSV = _report
    if _argv and _argv[0] == '--detail':
        if len(_argv) < 2:
            raise SystemExit("--detail needs at least one level id")
        for _lid in _argv[1:]:
            print_level_detail(_lid)
            print()
    elif _argv:
        raise SystemExit(f"unrecognised argument(s): {' '.join(_argv)}")
    else:
        main()
