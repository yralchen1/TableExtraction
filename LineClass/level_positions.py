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
    rho_t     the density near nu_t, per cm^-1, of the lines that are NOT
                accounted for by the level list - the unrelated ones, which
                are what a coincidence is drawn from (section 4a)
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
an open plate has N(0, sigma)/rho of about 40 for the typical sigma = 0.1 cm^-1
and rho = 0.039 unrelated lines per cm^-1 of this spectrum, so it is worth
ln R_t = +3.7, and ten of them e^37.  A missing prediction costs ln(1 - p), which is bounded,
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

with g the local distribution of ln I_obs over the UNRELATED lines - not over
the observed list as a whole, which is dominated by the identified and
therefore bright ones - since that is what a coincidence would have been drawn
from (section 4a).  Both numerator and
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

Two readings of a match, and the row is worth the better of them.  A matched
free feature can be read as being the transition - position drawn from
N(d; 0, sigma_t) against the local line density rho, brightness from
N(ln I_obs; ln(I_t f), s) against g - or as a feature that is there in any
case, carrying light nobody has identified, with the transition hidden in it.
Under the second reading H0 has the line present with certainty, so rho_t is
1/(2 W_t) as for a blend and the brightness says nothing either way; its
attainable value is 3.2 against the first reading's 25.  R_t takes whichever
is larger.

Taking the floor of the intensity term without also replacing rho was the
error: it read the brightness under the second hypothesis and the position
under the first, so a prediction of I = 0.2 landing on a feature of 13372 kept
the whole of the positional credit - ln R_t = +1.4 for a transition with a
1.5 per cent chance of having been recorded at all.  Rows of that kind are the
strong observed line assigned to the very weak transition, the assignment
classify_lines refuses and the analyst leaves free for a better one, and half
a dozen of them carried four levels of the audit to alternate positions.  A
genuine line under-predicted by a factor of thirty is unharmed: the first
reading is still much the better one for it, and it goes on being read that
way.  A row the second reading wins counts as neither free support nor
evidence in the level offset.

The second reading is only open to a feature at least as bright as the
transition it would be hiding: ln I_obs >= ln(I_t f) - s, one scatter width
of allowance.  A feature cannot carry light it does not have, and without the
condition a line nine times FAINTER than predicted - 7801 recorded where
70072 was predicted, for 059003.000391 at 116617.18 - was credited +1.5 as a
transition hidden in it and tagged `bright`, which is the opposite of what it
is.  Such a row is read as the transition or not at all, and pays for the
deficit.

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

sigma_t.  The wavenumber uncertainty of the matched line, built from the three
things that are known to widen a feature, and combined with the partner energy
uncertainty:

    sigma_t^2 = (k(n) sigma_meas)^2 + w_hfs(low)^2 + w_hfs(upp)^2 + u_M^2

    sigma_meas^2 = (dlam(era) 1e-8 nu^2)^2 + u_F^2 + (dlam(char) 1e-8 nu^2)^2
                                                     + k(char)^2

The first two terms are the measurement itself: a wavelength-constant reading
error dlam - 0.0030 A for the 1974 measurements, 0.0040 A for the 1969 ones -
converted to wavenumber by 1e-8 nu^2, plus a constant precision floor
u_F = 0.0055 cm^-1.  A wavelength-constant width grows as nu^2, so the two
cross at about 43000 cm^-1.  Both are measured, not assumed: among the plain
lines whose quoted uncertainty follows that rule the residual scatter needs
nothing whatever beyond it (1579 lines of 1974, rms of residual/sigma = 1.00).

The last two terms are the line's recorded character.  Every blend and
resolution code - w wide, d double, bl blended, ch complex and hazy, c complex
(1974), cl perturbed by a close neighbouring line (1969, the NIST code 'p'),
h hazy - widens the feature by an amount of its own, fitted by maximum
likelihood on the accepted lines of the run.

Which form that amount takes is a question about where the width comes from,
and the two possible answers are physically different.  A width that comes
from READING the plate - deciding where the middle of a smeared or doubled
image lies - is a length on the plate, so it is constant in wavelength and
grows as nu^2 in wavenumber; these are CHAR_DLAM, in A.  A width that belongs
to the LINE - an unresolved structure inside the feature itself - is a
constant in wavenumber, the same number of cm^-1 wherever the line falls;
these are CHAR_DWN, in cm^-1.  Each code is fitted both ways and given the
form its own residuals prefer.

Most codes choose the wavelength form, which is what they should choose: c is
a complex feature and cl a line pulled off centre by a close neighbour, and
both are difficulties of deciding where on the plate the line lies.  (They are
not the same code differently spelled - c is complex, cl is perturbed - and
their widths, 0.034 A and 0.0044 A, have no reason to agree.)  w (wide) chooses
the other form, by 12 units of ln L over its 312
lines, and says so twice over: split into four bands of wavenumber the fitted
excess is 0.023, 0.030, 0.019, 0.045 cm^-1 - flat within its errors - while
the same four bands expressed as a wavelength fall 0.0060, 0.0027, 0.0013,
0.0025 A, a factor of five across the range.  That is what 'wide' should mean
if the mark records a feature that is genuinely broader than the instrument
rather than one that was awkward to read, and it halves the width the old
wavelength fit gave at the blue end while nearly doubling it at the red.

These fits exclude the lines LOPT marks S = '*' (whose observed value IS the
Ritz value, by construction), F = 'P' (predicted, not observed), and those
whose Ritz wavenumber is not firm enough to test the observation
(uWnCtot > 0.5 uWnOTot); a blended feature enters once, through LOPT's own
centroid deviation dEcent, never once per component.

w_hfs is hyperfine structure, and belongs to the LEVEL, not to the line: an
unresolved hyperfine pattern displaces every line that touches that level, by
an amount set by the level's hyperfine splitting and therefore constant in
WAVENUMBER, not in wavelength.  The residuals say so plainly - over the band
holding 255 of the 285 flagged 1974 lines the excess is flat at 0.133 cm^-1
while its wavelength equivalent moves by 26 %, and the likelihood prefers the
wavenumber form by 138 units of ln L.  One width per level is fitted to the
1974 lines alone (fit_hfs_widths, --fit-hfs) and stored in
level_hfs_widths.csv; the levels that need more than HFS_MARK are marked.

A level with fewer than HFS_MIN_LINES lines of 1974 has no width the fit can
measure, and used to be given zero.  That is the one place the model was
plainly wrong: such a level is not a level without hyperfine structure, it is
a level whose hyperfine structure nobody has measured.  It is now given the
typical width of its CONFIGURATION instead (configuration_widths), which is
what the physics says governs a hyperfine splitting - an s electron reaches
the nucleus and feels its magnetic moment, a 5f or 5g electron does not - and
which the fitted widths bear out without being told: over the levels the fit
can constrain, 4f2.6s has a median width of 0.063 cm^-1 with 82 % of its
levels above the gate, while every configuration with no penetrating outer
electron has a median of zero; a permutation test on the configuration labels
gives p < 5e-5.  The rule is worth 445 units of held-out ln L over the 2057
accepted unblended 1974 lines, five-fold, against 22 (5-95 %: -42 to +129)
when the configuration labels are shuffled among the levels, and 462 of those
445 are earned by 31 lines alone: the ones Sugar marked *r or *v, which are
given no character width precisely because w_hfs is meant to carry them, and
whose level the fit could not reach.  It changes nothing for a level the fit
can constrain, and u_F is unmoved by it - re-fitted on the clean plain 1974
lines with the configuration widths in place it comes out 0.0054 (0.0049 -
0.0060) cm^-1, against 0.0050 (0.0045 - 0.0056) without them.
The fit rediscovers Sugar's own flags without being shown them - 89 % of the
*r lines and 61 % of the *v lines land on a marked level, against 4 % of the
unflagged ones - and then predicts the 1969 region out of sample: unflagged
1969 lines on marked levels go from rms 1.60 to 1.05 when the width is carried
across unchanged, and lines on unmarked levels do not move.  So the flag itself
is not used: a line gets the hyperfine width of its two levels whatever its
character, which is the only way the unflagged ones can be reached.

k(n) is the rms of t = (wn_obs - rwn)/sigma_meas among features carrying n
accepted transitions - about 1.01 unblended and 1.47 for a two-component blend,
so a blend is less well placed than the model for its width admits.  It is
shrunk toward 1 by a pseudo-count so that a class of five lines does not
acquire a factor of its own.  u_M is the partner level D1 from
LOPT_output_levels.txt, 0.002 - 0.03 cm^-1; adding it in quadrature slightly
overstates sigma, since the two are positively correlated, and the direction is
conservative.

The quoted unc_wn_obs of the line list is used only where the model cannot
speak for it.  A flagged line's quoted value is the reading rule times a fixed
factor of its code and carries nothing the code does not already say, so the
model replaces it; a line with no character flag whose quoted value exceeds the
rule can only have been widened by a deliberate judgement made when the list
was compiled, so there the larger of the two is kept.  Two codes are dropped
outright: ** (multiply classified) duplicates k(n), which prices multiplicity
already, and *v / *r (hyperfine, shaded violet or red) are replaced by w_hfs,
which needs no character excess on top of it once the per-level widths are in
(rms 1.015 and 0.941 with nothing added).

eta is fitted by maximum likelihood on the accepted lines, as the mixing weight
that makes (1 - eta) N(d; 0, sigma) + eta rho fit their residuals best.

rho(nu) and g(ln I | lambda) are the rate and the brightness distribution of
the unrelated lines, and section 4a is about nothing else; s and s_L are the
within- and between-level spreads of the intensity residual over levels with at
least five accepted lines.

s_i, the width of the intensity term, is not one number for every line.  Each
predicted transition carries its own uncertainty on the calculated intensity -
u_calc, the rms over Cowan Monte Carlo trials, on the natural-log scale - and

    s_i = sqrt(s0^2 + (k u_calc_i)^2),   s0 = a0 + a1 log10 P,   k = c log10 S

is used in its place (fit_intensity_width).  P = I_pred f / I_thr is the
intensity the calculation expects on the plate, in units of the detection
threshold at that wavelength, and S the calculated line strength in atomic
units.  s0 is the floor no per-line uncertainty explains, and it falls slowly
as the expected plate intensity grows; k is how much of u_calc to believe, and
it falls as the transition grows stronger: Cowan's u_calc is about right for
the weakest lines (log10 S near -5, k near 1) and overstated for the strong
ones.  k is never allowed below K_FLOOR = 0.2 - even the strongest
transition's calculation is not taken as exact - and s0 never below
S0_FLOOR.  The two linear
forms were chosen on the run of September 2026 from free fits of s0 and k in
bins: a quadratic term in k was not defined (87 per cent error once the bin
errors were scaled to a reduced chi^2 of 1), and a constant term beside c
log10 S was consistent with zero.  That run gave a0 = 0.902(27), a1 =
-0.062(16) per dex and c = -0.218(11), on 4245 lines.

a0, a1, c and a free mean are fitted together by maximum likelihood on the
single, non-bl accepted lines, with each residual's normal density cut off at
the detection threshold, so the faint predictions - accepted only where a line
happens to be bright enough to have been recorded - enter without their
selection being read as scatter.  Everything is refitted on each run, so a new
calculation behind u_calc needs no maintenance here.  Where a transition has
no u_calc the pooled s is used.


4a.  rho and g: the lines that belong to nothing yet known
----------------------------------------------------------
rho_t is the rate at which a line the level under test has nothing to do with
turns up near nu_t, and g is the brightness such a line would have.  Both were
first measured on the whole recorded list, from the K nearest lines.  That is
the wrong population, and it is wrong in both halves.  Of the 6668 recorded
lines, 4526 already carry an accepted transition: they are identified, they are
by construction explained without the level under test, and they are also the
bright ones, because a line was identified in the first place partly by being
strong enough to measure well.  Only the remaining 2142 are candidates for a
coincidence.  Using all 6668 puts rho about two and a half times too high, and
puts the mean of g about 0.6 in ln I too high.

The free lines are not a mystery.  A recorded line that carries no accepted
transition is, in this spectrum, almost always a transition of a level that the
calculation predicts and nobody has yet found: 658 of the 1253 calculated
levels of IDEN2/enlev.dat are in that state, 61290 of their E1 transitions have
a partner that HAS been found, and summing the probability P_t that each of
those would have been recorded gives about 2970 expected lines - of the same
order as the 2142 free ones actually there.  Cowan's full transition list,
tp_E1_no_trials.xlsx, therefore predicts the very population rho and g are
supposed to describe, and it predicts its structure as well as its size: the
expected rate runs from 0.001 per cm^-1 below 10000 to 0.052 near 55000, more
than an order of magnitude, which no single number and no 41-nearest-neighbour
smoothing of 2142 sparse lines can express.

An unfound level's energy is known only to W, the rms of E_obs - E_calc over
the found levels of its configuration, 40 to 800 cm^-1, so a single one of its
transitions predicts no position.  Smeared over W it predicts a RATE, which is
exactly what a Poisson background needs:

    rho_unk(nu) = sum over unfound-level transitions of
                  P_t * N(nu; nu_calc, W_cfg)                 per cm^-1

and, with the same weights, the mean and spread of ln(I_t f) give the
brightness such a line would have.  unknown_transition_background() builds
both, as one convolution per configuration width.

Two estimates of one population, then, and they are combined by taking the
LARGER of them - the larger rate, and the larger of the two brightness
densities.  That is the conservative direction, since a bigger background is a
smaller ln R, and it is the same device as the two readings of a match in
section 3: an unrelated line is allowed whichever account of itself is the
better, the one the free lines actually recorded suggest or the one the
calculation predicts for a level still missing.  The result is capped above by
the all-lines density, which the free lines are a subset of, and floored at
RHO_BG_FLOOR of it, because both estimates extrapolate where no free line is
near and an extrapolated density must not buy unbounded credit; the brightness
half is floored the same way, at LN_BG_DROP below the all-lines value.

Where this does NOT apply.  A CLAIMED feature keeps the all-lines rate and the
all-lines g: there the question is the blend one, rho_t has already been
replaced by 1/(2 W_t), and the one-sided floor of ln_intensity is a floor - a
guard against an astronomically small denominator - for which the more generous
distribution is the safer one.  Nor does rho_unk enter the "hidden in a feature
that is there anyway" reading.  That reading compares a world in which the
observed feature is an unrelated line with one in which the same unrelated line
has the transition blended into it; the feature is observed in both, so its
rate cancels and only the probability that the transition falls inside it,
2 W N(d; 0, sigma), survives.

What it costs and what it changes.  Reading the transition list adds about two
seconds to build().  Because the correction removes a background that was too
large, every level's ln R rises - the median of this run from 22.4 to 31.6, a
little under one unit per matched free line - so FIRM_LN_R and any threshold
read off an absolute ln R mean something slightly different than they did
before.  Differences between two positions of the SAME level, which is what the
scan and the audit are about, are much less affected: the shift is common to
both positions wherever they match the same number of free lines.

The effect worth knowing about is on the alternate positions.  Summed over the
593 levels, the local maxima within --alt-drop of the adopted peak fall from
2823 to 1912, and the levels carrying a question mark from 85 to 68: a third of
the alternate positions were manufactured by the background itself, a
coincidence scored against a rate two and a half times too low looking like a
match, and enough of those in one window raising a maximum.  That also loosens
the look-elsewhere correction, which is ln(n_alt) - 059003.000228 becomes firm
grounds for relocation not because its alternate gained but because the rival
maxima in its window fell from 107 to 10 - so the corrected gains want reading
with that in mind.  --plain-background restores the old behaviour exactly.

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

THE REGISTRY OF SETTLED POSITIONS.  The scan is what makes a run slow, and
most of it is spent re-proving what was proved last time: a level with ln R
above FIRM_LN_R and no alternate anywhere in its window is not in doubt, and
the great majority of the levels are in that state.  Those levels are written
to a registry file, one line each, with a FINGERPRINT of everything their scan
depends on: the interval scanned and its step, `--alt-drop`, the level's own
energy, every prediction it makes (partner, partner energy, partner
uncertainty, predicted intensity), the measured line list, the constants
fitted to the run, and the accepted assignments lying anywhere the level's
predictions can reach as it moves across its window.  With `--use-firm` a
level whose fingerprint is unchanged is not scanned again.

The fingerprint is deliberately LOCAL.  A run in which one relocation has been
accepted differs from the previous one in that level, in its partners, and in
the lines near it; a fingerprint over the run as a whole would invalidate all
594 entries and the registry would never save anything.  As written, accepting
a relocation costs a re-scan of the levels that share lines or partners with
it and of nothing else.

Two things guard against a stale entry being honoured.  The fingerprint is the
first.  The second is that ln R at the adopted position is computed for every
level on every run whether it is scanned or not - it is one evaluation, not a
scan - so a registry entry is used only if that value still comes out within
FIRM_TOL of the one recorded beside the fingerprint.  A level taken from the
registry is marked `registry` in the `scanned` column; one that was scanned is
marked `this run`.  `--no-firm` ignores the registry entirely, which is what to
use when a result must not depend on any earlier run.

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
OBSERVABLE (P_obs >= P_SEEN, section 2), well-centred (within FREE_SIGMA) and
sit on features that no accepted transition claims (C = 0).  Those are the
lines the level could take without robbing another level.  The observability
condition is the same one that separates n_seen from n_match, and it matters
as much here: a prediction with a 1.5 per cent chance of having been recorded
that nevertheless finds a line has found a coincidence, not support, and
counting such rows inflates n_free without adding anything to free_gain -
they are worth a few hundredths of a unit each.  Support that is entirely
blended is not support: a
prediction can be dumped on an already-explained feature almost anywhere, and
the blend branch of the formula gives it credit for doing so.  Alongside it,
top_share - the largest single row's share of all the positive evidence -
throws out the positions whose case is one lucky line.

THE DELTA J FINGERPRINT (ln_J at the adopted position, ln_J_alt at the
alternate).  For an E1 transition J changes by 0 or 1, so which partners a
level is actually recorded with is a fingerprint of its own J - and ln R cannot
read it, because the two energies being compared are the same level and predict
the same transitions to the same partners with the same gA.  The J information
survives only as a correlation among the absences, which ln R treats as
independent.  ln_j_pattern scores it against the rival hypothesis that the
lines belong to a level of another J, by fitting a detection multiplier to each
delta J class; see the block above ln_j_pattern for what it does and does not
claim.  It is weighed throughout by P_obs, so a class too faint to have been
recorded cannot accuse a position.

It enters ln R as a DIFFERENCE BETWEEN TWO POSITIONS OF THE SAME LEVEL and
never on its own.  A single position's ln_J is dominated by how well the
calculated gA divides that level's strength among its branches - a property of
the wavefunction, the same at every energy - so charging it to one position
would charge the intensity model's error to the level's whereabouts.  Between
two candidate energies of one level the partners and the gA are identical and
that fault cancels, leaving only what the test is for.  So every set of rival
positions is offset by the cleanest of them: with b = min ln_J over the set,
the position at energy E is judged on

    ln_R_J = ln_R - (ln_J(E) - b)

which leaves the best-fitting position's ln R untouched, charges each rival
what its fingerprint is worse by, and gives a level with no rival nothing to
pay.  With --audit the set is the adopted position and its alternate, so
gain = (ln_R_alt - ln_R) + d_ln_J; with --unknown it is the candidate positions
of the scan, which are then ranked, corrected for look-elsewhere and given
their verdicts on ln_R_J.

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

What survives is sorted into six dispositions, named rather than ranked
because they are six different pieces of work: `interchange` (the alternate
lands on another level of the run - level_interchange.py's question, on
evidence this scan does not use), `refit` (the level stays where it is and
some accepted line is dragging the LOPT fit), `top line` (the same short move,
but the line it gives up is the level's strongest), `relocate` (the level
moves: a free position, broadly supported, surviving the look-elsewhere
correction), `weak` (a move whose support is thin) and `no support`.

WHAT SEPARATES A REFIT FROM A RELOCATION is not how far the alternate lies -
it is whether the level keeps the recorded lines it is assigned now.  n_own
counts those lines and n_kept how many of them the alternate still matches.
A level that keeps most of them stays where it is: one bad assignment is
pulling the LOPT fit a fraction of a wavenumber off the position its own
lines want, and correcting it is a refit and nothing more (059003.000457
keeps six of its seven at a position 0.5 cm^-1 away).  A level that keeps
none of them has to be re-identified however short the move: every one of its
present assignments is dropped, another set is made, and the level is entered
in a ledger at a new position (059003.000617 keeps NONE of its three at a
position 1.3 cm^-1 away).  Sorting the two by distance put them side by side
under the same heading and told the reader to do the opposite of the work.

AND COUNTING THE LINES IS NOT ENOUGH, because the one being given up is
usually the one that matters.  kept_light weighs what n_kept counts: the share
of the level's own observed light the alternate still matches, each feature's
measured intensity split among its accepted components by branching fraction
so that a level is never credited with the whole of a blend.  059003.000457 is
the case - six of its seven lines survive a move of 0.51 cm^-1, which reads as
one bad line dragging the fit, but the line it drops is its strongest by an
order of magnitude (kept_light = 0.155) and the Ritz mismatch the move removes
is the mismatch of the level's brightest branch.  Removing it improves the fit
by construction and leaves the level deprived of the transition that most
defines it, which is a thing to explain and not a thing to discard.  That row
is called `top line`, and what wants opening in IDEN2 is the line, not the
position.

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
                 ALTERNATE, AMONG THE LEVELS OF THE SAME J AND PARITY.  Only
                 those can be the other half of an interchange: an interchange
                 exchanges two identities, and a level of J = 5/2 cannot take
                 over the lines of a level of J = 17/2 however close the two
                 energies are - the selection rules give the two of them
                 different transitions altogether.  A level of another J
                 sitting on the alternate is not a rival and is not named
                 here; the column is empty when the run has no other level of
                 that J and parity.  Even within one J it is a neighbourhood
                 label rather than a rival: it is filled in for every
                 alternate, and for most of them it names a level that merely
                 happens to lie nearby and wants nothing.  Only when near_dE
                 is small - under INTERCHANGE_DE, half a wavenumber - does it
                 mean anything, and then it means the alternate IS that
                 level's position, so the two identities may need exchanging
                 rather than either of them moving.  A near_level several
                 cm^-1 away is noise; read dE_alt instead.
    near_dE      E_alternate - E(near_level), cm^-1
    question     '?' when n_alt > 0
    scanned      `this run` when the window was scanned, `registry` when the
                 level was taken from the registry of settled positions
                 unchanged (section 5); the other scan columns then hold what
                 the registry recorded

with --audit, for that same alternate:

    n_free       OBSERVABLE matched lines at the alternate that no accepted
                 transition claims and that sit within FREE_SIGMA of
                 prediction.  Faint matches are coincidences and are not
                 counted, for the reason n_match is not quoted
    free_gain    what those free lines are worth in ln R
    top_share    the largest single line's share of the positive evidence
    gain_R       ln_R_alt - ln_R, the line evidence alone
    gain         gain_R + d_ln_J, what the alternate wins once the delta J
                 fingerprint is counted with it.  THE one that decides
    look         gain - ln(n_alt), the look-elsewhere correction of section 6
    n_obs_alt    observable predictions at the alternate
    n_seen_alt   how many of them would find a line there
    n_miss_alt   how many would be absences there
    ln_J         the delta J fingerprint at the adopted position, ln_J_alt at
                 the alternate: how much better the pattern of matches and
                 absences fits a level of a DIFFERENT J.  Zero is clean
    d_ln_J       ln_J - ln_J_alt, which is what gain counts.  A single
                 position's ln_J is dominated by how well the calculated gA
                 divides the level's strength among its branches, which is the
                 same at every energy and cancels in the difference.  Positive
                 means the alternate fits this level's J better
    ln_R_J       ln_R and ln_R_alt each offset by the cleaner fingerprint of
    ln_R_alt_J   the two, so that their difference is gain
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
before the levels giving up their strongest line before the weak ones, and within each group the largest look-elsewhere
corrected gain; with --scan alone, the largest gain; with neither, the weakest
positions first, since a plain report is a list of complaints.  The csv keeps
that order.  The tables printed under "the twenty weakest positions" are
sorted by ln R regardless, because that is what they are for.

"""

import argparse
import hashlib
import math
import os
import sys
from typing import NamedTuple

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import chance_mc as mc                        # noqa: E402
import output_files                          # noqa: E402
import classify_lines as cl                   # noqa: E402
import level_interchange as li                # noqa: E402
import cowan_gA                               # noqa: E402
import level_shifts as ls                     # noqa: E402
import lopt_lines                             # noqa: E402


# --- constants --------------------------------------------------------------
N_SIGMA = 4.0        # half-width of the matching window, in sigma
P_SEEN = 0.2         # P_obs at which a prediction counts observable
Q_OUT = 2.0 * 0.5 * math.erfc(N_SIGMA / math.sqrt(2.0))   # 6.3e-5
K_RHO = 41           # recorded lines averaged for the local line density
K_BG = 201           # recorded lines averaged for the local intensity spread
K_UNC = 101          # recorded lines averaged for the local uncertainty
K_RHO_FREE = 21      # free lines averaged for the local free-line density
K_BG_FREE = 101      # free lines averaged for their local intensity spread
UNK_STEP = 20.0      # cm^-1, the grid the modelled unknown-level background is
                     # accumulated on; far finer than the configuration windows
                     # that smear it, which are 40 - 800 cm^-1
UNK_SMEAR_MIN = 10.0 # cm^-1, the narrowest smearing used for an unfound level,
                     # so that no configuration ever contributes a spike
RHO_BG_FLOOR = 0.05  # the background rate of unrelated lines is never taken
                     # below this fraction of the all-lines density: both
                     # estimates of it extrapolate where there is no free line
                     # nearby, and this caps the credit that extrapolation can
                     # buy at ln 20 = 3.0 per matched row
LN_BG_DROP = 3.0     # and the same cap on the intensity half: the H0 density
                     # of a free feature's brightness is never more than this
                     # far below the all-lines value
SHRINK_N = 25.0      # pseudo-count pulling k(n) toward 1
K_WORST = 1.6        # the largest k(n) the run produces, used to widen the
                     # matching window (see ln_ratio)
ETA_MAX = 0.25       # upper bound of the anomalous-position rate
MIN_LEVEL_LINES = 5  # accepted lines a level needs to enter the s_L estimate
S0_FLOOR = 0.05      # the smallest s0 the intensity-width model may take
K_FLOOR = 0.2        # the smallest k: even the strongest transition's
                     # u_calc is believed to this extent

# --- the wavenumber-uncertainty model of section 4 -------------------------
UNC_FLOOR = 0.0055   # cm^-1, the precision floor every measurement carries
ERA_SPLIT = 47500.0  # cm^-1: above it the 1969 measurements, below the 1974
DLAM_ERA = {'1974': 0.0030, '1969': 0.0040}   # A, the plain reading error
CHAR_DLAM = {        # A, the wavelength-constant excess of each character,
                     # fitted jointly with the per-level hyperfine widths
    ('1974', 'd'): 0.0073,     # double
    ('1974', 'ch'): 0.0103,    # complex and hazy
    ('1974', 'c'): 0.0337,     # complex
    ('1974', 'bl'): 0.0085,    # blended (2 lines; the 1969 value)
    ('1974', 'h'): 0.0082,     # hazy (3 lines; the 1969 value)
    ('1969', 'cl'): 0.0044,    # perturbed by a close neighbouring line
                               # (Sugar's own words); the NIST code 'p'
    ('1969', 'h'): 0.0082,     # hazy
    ('1969', 'bl'): 0.0085,    # blended
    ('1969', 'w'): 0.0104,     # wide: 21 lines in the list, one of which
                     # survives the Ritz cut, so the 1974 refit below could
                     # not be repeated here; the wavelength value is kept
                     # because it is the larger of the two
}
CHAR_DWN = {         # cm^-1, the WAVENUMBER-constant excess of a character:
                     # a width belonging to the line itself rather than to the
                     # reading of the plate (see section 4)
    ('1974', 'w'): 0.0251,     # wide (312 lines)
}
HFS_FILE = os.path.join(HERE, 'level_hfs_widths.csv')
LEVEL_IDS = os.path.join(HERE, 'IDEN2', 'IDEN_level_ids.txt')
HFS_APPLY = 0.02     # cm^-1: a fitted width below this is noise of the fit
                     # and is not applied to any line
HFS_MARK = 0.08      # cm^-1: above this the level is reported as having a
                     # large hyperfine structure
HFS_MIN_LINES = 4    # 1974 lines a level needs before its width is believed
HFS_CFG_MIN = 3      # levels of its own the configuration needs before its
                     # typical width is used for a level the fit cannot
                     # constrain; below that the list-wide width is used
HFS_PRIOR = 0.02     # cm^2: the scale of the exponential prior on the fitted
                     # variances, which pulls an unsupported level to zero
HFS_MAX = 4.0        # cm^2, the upper bound of one level's fitted variance
HAND_FACTOR = 1.5    # how far above the usual quoted value of its own class a
                     # line must be quoted before the quoted value is read as a
                     # deliberate widening and honoured
GRID_STEP = 0.02     # cm^-1, the scan step
GRID_MAX = 300000    # cap on the number of scan points per level
BLOCK = 4000         # scan points evaluated in one array operation
ALT_SEP = 0.5        # cm^-1: two maxima closer than this are one maximum
WINDOW_K = 3.0       # configuration widths either side of E_calc a scan covers
ALT_DROP = 5.0       # ln R below the adopted peak that still counts as an
                     # alternate position
LOPT_LEVELS_FILE = os.path.join(HERE, 'LOPT_output_levels.txt')
FIRM_LN_R = 3.0      # ln R at or above which a level with no alternate
                     # position is settled enough to go into the registry
FIRM_TOL = 0.5       # a registry entry is honoured only if ln R at the
                     # adopted position still comes out within this of the
                     # value that was recorded with it
FIRM_FILE = 'level_positions_firm.csv'
REACH_PAD = 5.0      # cm^-1 added to the wavenumber range a level's
                     # predictions sweep, to cover the matching window

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
    return np.where(intensity_floored(x, predicted, s, ln_bg), ln_bg, ln_n)


def intensity_floored(x, predicted, s, ln_bg):
    """Where the floor of ln_intensity binds: the feature is already light.

    True marks a feature so much brighter than everything predicted on it
    that an unrecognised line accounts for its brightness better than the
    prediction does.  Charging that excess against the position would be
    wrong, which is what the floor is for; crediting the position with having
    found the feature is the same mistake with its sign reversed, and that is
    what the second reading of a match in ln_ratio is for.
    """
    m = np.log(np.maximum(np.asarray(predicted, dtype=float), 1e-300))
    return (x > m) & (ln_bg > ln_norm(x, m, s))


def ln_norm(x, mu, sd):
    """ln of the normal density, elementwise."""
    z = (np.asarray(x, dtype=float) - mu) / sd
    return -0.5 * (z * z) - np.log(sd) - 0.5 * _LN_2PI


# ---------------------------------------------------------------------------
# The run, and everything measured on it
# ---------------------------------------------------------------------------
class Context:
    """The run under test with every ingredient of the likelihood attached."""

    def __init__(self):
        # J of every level, for the delta J fingerprint; a context built
        # without it simply has no fingerprint to report
        self.j_of = {}


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
        # never below 1: the width model of section 4 is what the measurement
        # itself is worth, and a class whose residuals happen to fall inside it
        # buys no licence to quote the position better than it was measured
        k = max(k, 1.0)
        out[nb] = k
        rows.append((nb, cnt, math.sqrt(ms), k))
    log("  blend inflation of the position uncertainty, k(n):")
    for nb, cnt, rms, k in rows:
        log(f"    n={nb:2d}  N={cnt:5d}  rms(t)={rms:.3f}  k={k:.3f}")
    return out


# ---------------------------------------------------------------------------
# The width of an observed feature (section 4)
# ---------------------------------------------------------------------------
def normalize_char(c):
    """The character code of a line, reduced to the published vocabulary.

    The raw column mixes three kinds of mark.  A trailing '**' says the line is
    multiply classified, which is not a property of the feature at all and is
    priced by k(n); it is stripped.  '*v' and '*r' say the feature is shaded
    violet or red by an unresolved hyperfine pattern; they become 'hfs', which
    carries no width of its own because the per-level widths carry it.  What
    is left is one of the resolution codes - w, d, bl, ch, c, cl, h - or the
    empty string for a line with no mark.
    """
    s = str(c or '').strip()
    if s in ('*v', '*r'):
        return 'hfs'
    while s.endswith('*'):
        s = s[:-1]
    return s.strip()


def era_of(wn):
    """'1969' above ERA_SPLIT, '1974' below: which measurement a line is from."""
    return np.where(np.asarray(wn, dtype=float) > ERA_SPLIT, '1969', '1974')


def meas_sigma(wn, char, unc=None):
    """sigma_meas: the reading error, the precision floor and the character.

    Vectorized over arrays of wavenumber and of raw character code.  `unc` is
    the quoted unc_wn_obs of the line list; it is honoured only where it says
    something the model cannot - see the last paragraph of section 4.  The test
    is made against the line's own class: within a class the quoted value is
    the reading rule times one fixed factor, so a line quoted more than
    HAND_FACTOR times the usual value of its class was widened by hand when the
    list was compiled, on knowledge no flag records, and its quoted value is
    then used as a floor.  The medians are taken from the lines passed in, so
    nothing about the factors has to be assumed.
    """
    wn = np.asarray(wn, dtype=float)
    code = np.array([normalize_char(c) for c in np.atleast_1d(char)])
    era = era_of(wn)
    dl = np.array([DLAM_ERA[e] for e in np.atleast_1d(era)])
    dc = np.array([CHAR_DLAM.get((e, c), 0.0)
                   for e, c in zip(np.atleast_1d(era), code)])
    kc = np.array([CHAR_DWN.get((e, c), 0.0)
                   for e, c in zip(np.atleast_1d(era), code)])
    rule = np.sqrt(((dl * wn ** 2) * 1e-8) ** 2 + UNC_FLOOR ** 2)
    s = np.sqrt(rule ** 2 + ((dc * wn ** 2) * 1e-8) ** 2 + kc ** 2)
    if unc is None:
        return s
    unc = np.asarray(unc, dtype=float)
    good = np.isfinite(unc) & (unc > 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(good, unc / rule, np.nan)
    typical = np.ones(len(wn))
    for key in set(zip(era, code)):
        m = (era == key[0]) & (code == key[1])
        if m.sum() >= 5 and np.isfinite(ratio[m]).any():
            typical[m] = np.nanmedian(ratio[m])
    hand = good & (ratio > HAND_FACTOR * np.maximum(typical, 1.0))
    return np.where(hand, np.maximum(unc, s), s)


def configuration_widths(w_hfs, n_1974, cfg):
    """The typical hyperfine width of each configuration, in cm^-1.

    A level's hyperfine splitting is set by how strongly its outer electron
    feels the nuclear magnetic moment, and that is a property of which
    orbitals the electron occupies - the level's CONFIGURATION - far more than
    of anything else about it.  An s electron has a non-zero probability
    density at the nucleus and feels the moment directly (the Fermi contact
    term); a 5f or 5g electron does not come near it.  141Pr is the only
    isotope, so there is no isotope shift to confuse this with, and the
    ordering the fitted widths come out in is the ordering that argument
    predicts: 4f2.6s at 0.111 cm^-1, then the configurations with a 6p or a
    5d2, then the pure 4f3 and the non-penetrating 4f2.5f and 4f2.5g at zero.

    So a level with too few lines of its own for a width to be fitted is not
    a level with no hyperfine structure; it is a level whose hyperfine
    structure has not been measured, and the best estimate of it is what the
    levels of the same configuration show.  That is what this returns: for
    each configuration the root mean square of the widths APPLIED to the
    levels of that configuration the fit can constrain - the mean of the
    variances, which is the quantity that adds in quadrature and therefore
    the one to predict with.  Levels below the gate count as the zeros they
    are, so a configuration whose levels are all narrow gets a narrow width.

    A configuration with fewer than HFS_CFG_MIN determined levels cannot speak
    for itself, and is answered for by the electron that decides the question
    anyway - the outermost one.  Its orbital is taken from the label and the
    levels of every configuration with that same outer orbital are pooled: 4f
    5d.6s has one determined level of its own and would say nothing, but its
    outer electron is the same 6s as 4f2.6s's, and that is a 6s electron's
    answer to give.  Failing that the levels are pooled by the orbital LETTER
    alone - every s, every d - which still separates a penetrating outer
    electron from one that never reaches the nucleus, and only failing that
    does the list-wide value apply.  Each step is wider than the truth more
    often than narrower, which is the safe direction: a width can only make a
    sigma larger and a position less certain.

    Returns (dict key -> width, the list-wide width).  The keys are
    configurations, then outer orbitals ("6s"), then orbital letters ("s");
    width_for() walks them in that order.
    """
    w = np.asarray(w_hfs, dtype=float)
    n = np.asarray(n_1974, dtype=float)
    cfg = np.asarray([str(c or '') for c in cfg])
    sup = n >= HFS_MIN_LINES
    var = np.where(w > HFS_APPLY, w, 0.0) ** 2
    wide = math.sqrt(var[sup].mean()) if sup.any() else 0.0
    out = {}
    if not sup.any():
        return out, wide
    orb = np.array([outer_orbital(c) for c in cfg])
    for keys in (cfg, orb, np.array([o[-1:] for o in orb])):
        for k in sorted(set(keys[sup])):
            m = sup & (keys == k)
            if k and k not in out and int(m.sum()) >= HFS_CFG_MIN:
                out[k] = math.sqrt(var[m].mean())
    return out, wide


def outer_orbital(cfg):
    """The outermost orbital of an enlev.dat configuration, as "6s" or "5d".

    IDEN2 runs the open shells together, each an orbital letter optionally
    preceded by its principal quantum number and followed by its occupation -
    "f26p" is 4f^2.6p, "fd6p" is 4f.5d.6p, "4f3" is 4f^3 - so a digit BEFORE a
    letter is that letter's principal quantum number and a digit after it is
    an occupation.  The last shell written is the outermost one, which is the
    electron whose penetration of the nucleus decides the hyperfine splitting.
    A principal quantum number IDEN2 leaves out is left out here too, so "f"
    and "d" are their own keys; that costs nothing, because the letter-only
    pooling catches them.
    """
    cfg = str(cfg or '')
    last = ''
    i = 0
    while i < len(cfg):
        ch = cfg[i].lower()
        if ch in li._ORBITAL_L:
            pre = cfg[i - 1] if i and cfg[i - 1].isdigit() else ''
            last = pre + ch
            if i + 1 < len(cfg) and cfg[i + 1].isdigit():
                i += 1
        i += 1
    return last


def width_for(cfg, widths, wide):
    """The width to give a level of this configuration that the fit could not
    measure: its configuration's, else its outer orbital's, else that
    orbital's letter, else the list-wide value."""
    orb = outer_orbital(cfg)
    for k in (str(cfg or ''), orb, orb[-1:]):
        if k and k in widths:
            return widths[k]
    return wide


def read_hfs_widths(path=HFS_FILE, log=None):
    """The hyperfine width applied to each level, from level_hfs_widths.csv.

    The file's w_applied column already holds the rule of section 4: a level
    the fit can constrain - HFS_MIN_LINES lines of 1974 or more - carries its
    own fitted width if that width clears HFS_APPLY, and a level the fit
    cannot constrain carries the typical width of its configuration instead of
    zero.  Files written before that column existed are read the old way, so
    an old level_hfs_widths.csv still works.
    """
    if not os.path.exists(path):
        if log:
            log(f"  no {os.path.basename(path)}: no hyperfine widths applied")
        return {}
    d = pd.read_csv(path, dtype={'level_id': str})
    if 'w_applied' in d.columns:
        d = d[d['w_applied'] > 0]
        col, kind = 'w_applied', d.get('w_source')
        n_cfg = int((kind == 'configuration').sum()) if kind is not None else 0
    else:
        d = d[(d['w_hfs'] > HFS_APPLY) & (d['n_1974'] >= HFS_MIN_LINES)]
        col, n_cfg = 'w_hfs', 0
    out = {str(k): float(v) for k, v in zip(d['level_id'], d[col])}
    big = sum(1 for v in out.values() if v > HFS_MARK)
    if log:
        log(f"  hyperfine widths: {len(out)} levels carry one, "
            f"{min(out.values(), default=0):.3f} - "
            f"{max(out.values(), default=0):.3f} cm^-1; {big} of them above "
            f"{HFS_MARK} cm^-1 are marked as having a large hyperfine "
            f"structure; {n_cfg} of them are the typical width of the level's "
            f"configuration, the fit having too few lines to measure its own")
    return out


def hfs_pair(w_hfs, low, upp):
    """w_hfs(low)^2 + w_hfs(upp)^2 for arrays of level identifiers."""
    a = np.array([w_hfs.get(str(x), 0.0) for x in np.atleast_1d(low)])
    b = np.array([w_hfs.get(str(x), 0.0) for x in np.atleast_1d(upp)])
    return a ** 2 + b ** 2


def read_level_configs(levels, log=None, enlev=None, ids=None):
    """level_id -> dominant configuration, from IDEN2/enlev.dat.

    The two lists are tied by IDEN2/IDEN_level_ids.txt, the table that gives
    every level identifier the number of its row in enlev.dat.  That table is
    the authority, and it is the one thing that survives a level being moved:
    when a position is edited in IDEN2 the energy in enlev.dat changes and the
    row number does not.  Matching on energy instead loses exactly those
    levels - four of them in the run of 2026-09-11, each shifted 0.5 - 1.1
    cm^-1 beyond level_interchange.ENLEV_MATCH_TOL.

    Energy matching is kept only for a level the table does not list, and a
    level that neither route finds is reported: every level of the set the
    pipeline and LOPT use is in enlev.dat, so a miss means the lookup is
    broken and wants a human, not a silent fallback width.
    """
    try:
        en = li.read_enlev(enlev or li.ENLEV)
        cfg_of_row = dict(zip(en['idx'].astype(int), en['cfg']))
    except Exception as exc:                      # no IDEN2 files, say
        if log:
            log(f"  no configurations ({exc}): every level takes the "
                f"list-wide hyperfine width")
        return {}

    want = [str(k) for k in levels['level_id']]
    out = {}
    try:
        row_of_id = {v: k for k, v in
                     cowan_gA.read_id_map(ids or LEVEL_IDS).items()}
    except Exception as exc:
        row_of_id = {}
        if log:
            log(f"  IDEN_level_ids.txt could not be read ({exc}); "
                f"falling back to matching on energy")
    for lid in want:
        row = row_of_id.get(lid)
        if row in cfg_of_row:
            out[lid] = str(cfg_of_row[row] or '')

    left = [k for k in want if k not in out]
    if left:                                      # not in the table: by energy
        try:
            win = li.configuration_windows(en)
            sub = levels[[str(k) in set(left) for k in levels['level_id']]]
            lv, unmatched = li.attach_identities(
                sub.assign(E_final=sub['E_input']), en, win)
            for k, c in zip(lv['level_id'], lv['cfg']):
                if str(c or ''):
                    out[str(k)] = str(c)
        except Exception as exc:
            if log:
                log(f"  energy matching failed ({exc})")

    missing = [k for k in want if not out.get(k)]
    if log:
        log(f"  configurations: {len(want) - len(missing)} of {len(want)} "
            f"levels read from enlev.dat through IDEN_level_ids.txt")
        if missing:
            log(f"  {len(missing)} levels are in neither the table nor the "
                f"energy window and take the list-wide hyperfine width - "
                f"every level of the fit should be in enlev.dat, so this "
                f"wants checking: {', '.join(missing[:8])}"
                + (' ...' if len(missing) > 8 else ''))
    return out


def fit_hfs_widths(acc, levels, log=print):
    """One hyperfine width per level, by maximum likelihood on the 1974 lines.

    The model is sigma^2 = base^2 + v(low) + v(upp), with base the measurement
    and character part of section 4 and v >= 0 one variance per level.  Only
    the 1974 lines are used: below ERA_SPLIT the measurement base is 0.006 -
    0.01 cm^-1 and a hyperfine width of 0.1 - 0.3 cm^-1 stands out of it, while
    above it the base alone is 0.29 cm^-1 at 85000 cm^-1 and the same width is
    invisible.  Carrying the widths across that boundary unchanged is then a
    genuine out-of-sample prediction, and it is what makes the treatment
    believable; see section 4.

    An exponential prior of scale HFS_PRIOR on the sum of the variances keeps a
    level with two lines and one large residual from acquiring a width.

    k(n) and the widths are fitted together, by alternating: a k(n) measured
    with no hyperfine term in the base comes out too large, because a blend
    that happens to sit on a wide level charges its width to the blending, and
    the inflated k(n) then leaves nothing for the level to explain.  Four
    passes settle it.

    Returns a frame with one row per level, ready to be written to
    level_hfs_widths.csv.
    """
    from scipy.optimize import minimize

    a = acc.dropna(subset=['rwn', 'wn_obs', 'low_id', 'upp_id']).copy()
    meas = meas_sigma(a['wn_obs'].to_numpy(),
                      a['char'].fillna('').astype(str).to_numpy(),
                      a['unc_wn_obs'].to_numpy(dtype=float))
    nb = a['n_accepted'].fillna(1).to_numpy(dtype=float)
    r_all = (a['wn_obs'] - a['rwn']).to_numpy()

    ids = sorted(set(a['low_id'].astype(str)) | set(a['upp_id'].astype(str)))
    idx = {k: i for i, k in enumerate(ids)}
    # every level the accepted lines touch, not only those of the level list:
    # a level added by this run carries lines before it carries a position
    e_of = dict(zip(levels['level_id'].astype(str), levels['E_input']))
    level_cfg = read_level_configs(
        pd.DataFrame({'level_id': ids,
                      'E_input': [e_of.get(k, np.nan) for k in ids]}), log)
    lo_all = a['low_id'].astype(str).map(idx).to_numpy()
    up_all = a['upp_id'].astype(str).map(idx).to_numpy()
    n = len(ids)
    prior = 1.0 / HFS_PRIOR

    # the fit of the widths themselves sees the 1974 lines only
    f74 = (a['wn_obs'] <= ERA_SPLIT).to_numpy()
    lo, up, r = lo_all[f74], up_all[f74], r_all[f74]

    v = np.zeros(n)
    for _ in range(4):
        # the widths the last pass could support, as the production model
        # would apply them
        cnt74 = pd.concat([a.loc[f74, 'low_id'].astype(str),
                           a.loc[f74, 'upp_id'].astype(str)]).value_counts()
        vt = np.array([v[idx[k]] if (math.sqrt(v[idx[k]]) > HFS_APPLY
                                     and cnt74.get(k, 0) >= HFS_MIN_LINES)
                       else 0.0 for k in ids])
        kn = fit_k_blend(r_all / np.sqrt(meas ** 2 + vt[lo_all] + vt[up_all]),
                         nb, lambda *_: None)
        b2 = (np.array([kn.get(int(x), 1.0) for x in nb]) * meas)[f74] ** 2

        def nll(x):
            s2 = b2 + x[lo] + x[up]
            val = 0.5 * np.sum(np.log(s2) + r ** 2 / s2) + prior * np.sum(x)
            g0 = 0.5 * (1.0 / s2 - r ** 2 / s2 ** 2)
            g = np.full(n, prior)
            np.add.at(g, lo, g0)
            np.add.at(g, up, g0)
            return val, g

        v = minimize(nll, np.full(n, 0.005), jac=True, method='L-BFGS-B',
                     bounds=[(0.0, HFS_MAX)] * n,
                     options={'maxiter': 2000}).x
    w = np.sqrt(np.maximum(v, 0.0))
    a = a[f74]

    cnt = pd.concat([a['low_id'].astype(str),
                     a['upp_id'].astype(str)]).value_counts()
    lv = levels.set_index(levels['level_id'].astype(str))
    out = pd.DataFrame({
        'level_id': ids,
        'w_hfs': np.round(w, 4),
        'n_1974': [int(cnt.get(k, 0)) for k in ids],
        'E_input': [lv['E_input'].get(k, np.nan) for k in ids],
        'J': [lv['J'].get(k, '') for k in ids],
        'parity': [lv['parity'].get(k, '') for k in ids],
        'cfg': [level_cfg.get(k, '') for k in ids],
    })

    # the width each level will actually be given: its own where the fit can
    # measure one, its configuration's where it cannot - see
    # configuration_widths
    wc, wide = configuration_widths(out['w_hfs'], out['n_1974'], out['cfg'])
    sup = out['n_1974'] >= HFS_MIN_LINES
    own = np.where(out['w_hfs'] > HFS_APPLY, out['w_hfs'], 0.0)
    fall = np.array([width_for(c, wc, wide) for c in out['cfg']])
    out['w_applied'] = np.round(np.where(sup, own, fall), 4)
    out['w_source'] = np.where(sup, 'level', 'configuration')
    out.loc[out['w_applied'] <= 0, 'w_source'] = ''

    used = out[out['w_applied'] > 0]
    marked = used[used['w_applied'] > HFS_MARK]
    n_own = int(((out['w_source'] == 'level') & (out['w_applied'] > 0)).sum())
    n_cfg = int((out['w_source'] == 'configuration').sum())
    log(f"  hyperfine widths fitted on {len(a)} lines of 1974 over "
        f"{n} levels; {len(used)} will be applied - {n_own} the level's own "
        f"fitted width above {HFS_APPLY} cm^-1 on {HFS_MIN_LINES} lines or "
        f"more, {n_cfg} the typical width of the level's configuration where "
        f"the fit has too few lines to measure one; {len(marked)} above "
        f"{HFS_MARK} cm^-1 are marked as large")
    log('  configuration widths, cm^-1: '
        + ', '.join(f'{c} {v:.4f}'
                    for c, v in sorted(wc.items(), key=lambda t: -t[1]))
        + f'; anything else {wide:.4f}')
    pooled = out.loc[out['w_source'] == 'configuration', 'cfg']
    pooled = sorted(set(c for c in pooled if c not in wc))
    if pooled:
        log('  answered for by the outer electron: '
            + ', '.join(f'{c} -> {outer_orbital(c) or "?"} '
                        f'{width_for(c, wc, wide):.4f}' for c in pooled))
    return out.sort_values('level_id').reset_index(drop=True)


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


def fit_intensity_scatter(r, level_of, log, u=None, lp=None, ls10=None,
                          rthr=None, fit=None):
    """(s, s_L, width): the spread of the intensity residual
    r = ln(I_obs BF / (I_pred f)).

    s is the pooled within-level spread and s_L the between-level one, both
    measured on every accepted line.  `width` is the per-line model of
    section 3, fitted by fit_intensity_width on the rows `fit` marks, with
    lp, ls10 and rthr the quantities it needs for each line (see there).
    """
    ok = np.isfinite(r)
    df = pd.DataFrame({'r': r[ok], 'lid': np.asarray(level_of)[ok]})
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
        s, s_l = float(np.std(df['r'], ddof=1)), 0.1
    log(f"  intensity residual: mean {float(df['r'].mean()):+.3f}, sd "
        f"{float(df['r'].std(ddof=1)):.3f} over {len(df)} accepted lines")
    log(f"  within-level s = {s:.3f}, between-level s_L = {s_l:.3f} "
        f"({len(big)} levels with {MIN_LEVEL_LINES} lines or more)")
    width = fit_intensity_width(r, u, lp, ls10, rthr, fit, s, log)
    return s, s_l, width


class IntensityWidth(NamedTuple):
    """sigma^2 = s0(P)^2 + (k(S) u_calc)^2, with

        s0 = a0 + a1 log10 P        (never below S0_FLOOR)
        k  = c log10 S              (never below k_min)

    p_ref and s_ref are the medians of log10 P and log10 S over the lines the
    fit was made on, used for a transition that lacks one of them; mu is the
    fitted mean residual, kept for the report only; k_min is K_FLOOR for a
    fitted model and 0 for the flat one."""
    a0: float
    a1: float
    c: float
    p_ref: float
    s_ref: float
    mu: float = 0.0
    k_min: float = K_FLOOR


def flat_width(s):
    """The model with no per-line information: sigma = s for every line."""
    return IntensityWidth(float(s), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def width_terms(w, lp, ls10):
    """(s0, k) of the model w at log10 P = lp and log10 S = ls10, arrays."""
    lp = np.asarray(lp, dtype=float)
    ls10 = np.asarray(ls10, dtype=float)
    lp = np.where(np.isfinite(lp), lp, w.p_ref)
    ls10 = np.where(np.isfinite(ls10), ls10, w.s_ref)
    s0 = np.maximum(w.a0 + w.a1 * lp, S0_FLOOR)
    k = np.maximum(w.c * ls10, w.k_min)
    return s0, k


def fit_intensity_width(r, u, lp, ls10, rthr, fit, s, log):
    """The per-line intensity width of section 3, by maximum likelihood.

    For each line, P = I_pred f / I_thr is the intensity the calculation
    expects on the plate in units of the detection threshold I_thr at that
    wavelength, and S the calculated line strength in atomic units; lp and
    ls10 are their log10.  rthr = ln(I_thr BF / (I_pred f)) is the residual a
    line exactly at the threshold would have.  A line fainter than the
    threshold is never recorded, so the residuals of the recorded lines are
    drawn from the normal density N(r; mu, sigma) cut off below rthr, and
    each line's likelihood is

        N(r; mu, sigma) / Phi((mu - rthr) / sigma)

    with Phi the standard normal cumulative distribution.  The cut-off is
    what makes the faint predictions usable: without it their selection - a
    faint prediction is accepted only where a line happens to be bright
    enough - would be read as scatter.

    a0, a1, c and mu are fitted together on the rows `fit` marks (single
    classified, not bl).  The 1-sigma errors of a0, a1 and c from the
    curvature of the likelihood are reported, with 2 dlnL against the flat
    model (a1 = 0, k = 0).  Falls back to flat_width(s) - one width for
    every line - when the ingredients are missing or too few rows remain.
    """
    if u is None or lp is None or ls10 is None or rthr is None:
        return flat_width(s)
    u, lp, ls10, rthr = (np.asarray(x, dtype=float)
                         for x in (u, lp, ls10, rthr))
    m = (np.isfinite(r) & np.isfinite(u) & (u > 0) & np.isfinite(lp)
         & np.isfinite(ls10) & np.isfinite(rthr))
    if fit is not None:
        m &= np.asarray(fit, dtype=bool)
    n = int(m.sum())
    if n < 200:
        log(f"  intensity width not fitted: only {n} usable lines; the single "
            f"width s = {s:.3f} is used")
        return flat_width(s)
    from scipy.optimize import minimize
    from scipy.special import log_ndtr

    rr, uu, pp, ss, tt = r[m], u[m], lp[m], ls10[m], rthr[m]

    def nll(q, k_min=K_FLOOR):
        s0 = q[0] + q[1] * pp
        if np.any(s0 < S0_FLOOR):
            return 1e12
        k = np.maximum(q[2] * ss, k_min)
        v = s0 ** 2 + (k * uu) ** 2
        return float(np.sum(0.5 * np.log(2.0 * np.pi * v)
                            + (rr - q[3]) ** 2 / (2.0 * v)
                            + log_ndtr((q[3] - tt) / np.sqrt(v))))

    opts = {'maxiter': 20000, 'maxfev': 20000, 'xatol': 1e-7, 'fatol': 1e-7}
    best = None
    for c0 in (-0.1, -0.2, -0.3):
        f = minimize(nll, [s, 0.0, c0, 0.0], method='Nelder-Mead',
                     options=opts)
        if best is None or f.fun < best.fun:
            best = f
    flat = minimize(lambda q: nll([q[0], 0.0, 0.0, q[1]], 0.0), [s, 0.0],
                    method='Nelder-Mead', options=opts)
    a0, a1, c, mu = (float(x) for x in best.x)
    gain = 2.0 * (float(flat.fun) - float(best.fun))
    if not all(np.isfinite([a0, a1, c, mu])) or gain <= 0.0:
        log(f"  the per-line width buys nothing (2 dlnL = {gain:.1f}); the "
            f"single width s = {s:.3f} is used")
        return flat_width(s)
    err = curvature_errors(nll, best.x)
    w = IntensityWidth(a0, a1, c, float(np.median(pp)), float(np.median(ss)),
                       mu)
    s0, k = width_terms(w, pp, ss)
    log(f"  per-line intensity width sigma^2 = s0^2 + (k u_calc)^2, fitted on "
        f"{n} single, non-bl lines cut off at the detection threshold:")
    log(f"    s0 = {a0:.4f}({err[0]:.4f}) {a1:+.4f}({err[1]:.4f}) log10 P, "
        f"{s0.min():.2f} to {s0.max():.2f} over those lines")
    log(f"    k  = max({c:+.4f}({err[2]:.4f}) log10 S, {K_FLOOR}), "
        f"{k.min():.2f} to "
        f"{k.max():.2f}; mean residual {mu:+.3f}")
    log(f"    2 dlnL against one width for every line: {gain:.1f}")
    return w


def curvature_errors(nll, x, h=1e-3):
    """1-sigma errors from the inverse of the numerical Hessian of nll at x;
    nan where it is not positive definite."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    e = np.eye(n) * h
    f0 = nll(x)
    H = np.zeros((n, n))
    for i in range(n):
        H[i, i] = (nll(x + e[i]) - 2.0 * f0 + nll(x - e[i])) / h ** 2
        for j in range(i + 1, n):
            H[i, j] = H[j, i] = (nll(x + e[i] + e[j]) - nll(x + e[i] - e[j])
                                 - nll(x - e[i] + e[j])
                                 + nll(x - e[i] - e[j])) / (4.0 * h ** 2)
    try:
        d = np.diag(np.linalg.inv(H))
    except np.linalg.LinAlgError:
        return np.full(n, np.nan)
    return np.where(d > 0, np.sqrt(np.abs(d)), np.nan)


def sigma_intensity(ctx, u, lp=None, ls10=None):
    """The per-line intensity width sqrt(s0(P)^2 + (k(S) u)^2).

    lp and ls10 are log10 P and log10 S of each transition (see
    fit_intensity_width); either may be omitted or NaN, and the median of the
    fitted lines then stands in.  Falls back to the pooled s wherever u is
    missing, so a transition whose calculation carries no uncertainty is
    judged exactly as before, and everywhere when the run has no fitted
    width.
    """
    u = np.atleast_1d(np.asarray(u, dtype=float))
    w = getattr(ctx, 'width', None)
    if w is None:
        return np.full(u.shape, ctx.s)
    lp = np.full(u.shape, np.nan) if lp is None else lp
    ls10 = np.full(u.shape, np.nan) if ls10 is None else ls10
    s0, k = width_terms(w, lp, ls10)
    v = np.sqrt(s0 ** 2 + (k * np.where(np.isfinite(u), u, 0.0)) ** 2)
    return np.where(np.isfinite(u) & (u > 0), v, ctx.s)


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
        ip = {(str(r[0]), str(r[1])): float(r[2])
              for r in ctx.preds_all.itertuples(index=False)}
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


def unknown_transition_background(ctx, log=print):
    """The modelled density of recorded lines that belong to unfound levels.

    Every level of the calculation that has not yet been found - 658 of the
    1253 rows of IDEN2/enlev.dat in this run - still has a place in Cowan's
    E1 transition list, and every one of its transitions to a level that HAS
    been found predicts a wavenumber, an intensity, and, through the coverage
    and detection curves, a probability P_obs that such a line would have been
    recorded.  Those are exactly the lines that make up the free part of the
    observed list: features nobody has identified because the level at one end
    of them is not yet known.

    The wavenumber is predicted only to the width W of the level's
    configuration - the rms of E_obs - E_calc over the found levels of that
    configuration, 40 to 800 cm^-1 - so a single transition does not predict a
    position.  Smeared over W it predicts a RATE, which is the object a Poisson
    background needs.  Each transition contributes P_obs of a line spread as a
    normal density of width W about its calculated wavenumber, and the
    accumulated rate is the density of unrelated lines that the calculation
    says should be there, per cm^-1.

    Returns ``(grid, rate, mu, sd)``: the wavenumber grid, the rate on it, and
    the mean and spread of ln I of those same transitions - their predicted
    intensities put on the observed scale by the far-ultraviolet correction
    f(lambda) - which is what an unrelated line's BRIGHTNESS would have been
    drawn from if it belongs to a level nobody has found.  Returns None if the
    transition list or the calculated level list is not available, in which
    case only the empirical estimate is used and the report says so.
    """
    try:
        import unfound_levels as uf
    except Exception as exc:                                  # pragma: no cover
        log(f"  unfound-level background unavailable: {exc}")
        return None
    quiet = (lambda *a, **k: None)
    try:
        en, trans, mapping = uf.read_theory(log=quiet)
        per_cfg, whole = uf.windows(en, log=quiet)
        e_meas = uf.measured_energies(en, e_final=ctx.e_final, log=quiet)
        t = uf.unfound_transitions(en, trans, mapping, e_meas, log=quiet)
        t = uf.attach_windows(t, en, per_cfg, whole)
        t = uf.observation_probabilities(t, ctx)
    except Exception as exc:                                  # pragma: no cover
        log(f"  unfound-level background unavailable: {exc}")
        return None
    if not len(t):
        return None

    wn = t['wn'].to_numpy(dtype=float)
    p = t['P_obs'].to_numpy(dtype=float)
    ipr = t['I_pred'].to_numpy(dtype=float)
    wid = np.maximum(t['W'].to_numpy(dtype=float), UNK_SMEAR_MIN)
    f = scale_factor_vec(wn, ctx.bias)
    keep = (np.isfinite(wn) & (wn > cl.WN_MIN) & (wn < cl.WN_MAX)
            & np.isfinite(p) & (p > 0) & (ipr > 0) & np.isfinite(f) & (f > 0))
    wn, p, ipr, wid = wn[keep], p[keep], ipr[keep], wid[keep]
    x = np.log(ipr * f[keep])

    edges = np.arange(cl.WN_MIN, cl.WN_MAX + UNK_STEP, UNK_STEP)
    grid = 0.5 * (edges[:-1] + edges[1:])
    rate = np.zeros(len(grid))
    s1 = np.zeros(len(grid))
    s2 = np.zeros(len(grid))
    # One convolution per distinct configuration width, which is a few dozen
    # kernels rather than 61290 of them.
    for w0 in np.unique(wid):
        m = wid == w0
        k = np.arange(-int(4.0 * w0 / UNK_STEP) - 1,
                      int(4.0 * w0 / UNK_STEP) + 2) * UNK_STEP
        ker = np.exp(-0.5 * (k / w0) ** 2)
        ker /= ker.sum()
        for acc, wt in ((rate, p[m]), (s1, p[m] * x[m]), (s2, p[m] * x[m] ** 2)):
            h, _ = np.histogram(wn[m], bins=edges, weights=wt)
            acc += np.convolve(h, ker, mode='same')
    with np.errstate(divide='ignore', invalid='ignore'):
        mu = np.where(rate > 0, s1 / np.maximum(rate, 1e-300), np.nan)
        var = np.where(rate > 0, s2 / np.maximum(rate, 1e-300) - mu ** 2,
                       np.nan)
    sd = np.sqrt(np.maximum(var, 0.04))          # floored at 0.2, as bg_sd_o is
    rate = rate / UNK_STEP                       # counts per bin -> per cm^-1
    log(f"  unfound levels: {int((~en['known']).sum())} of {len(en)} "
        f"calculated; {len(t)} of their transitions have a found partner and "
        f"{p.sum():.0f} of those would have been recorded")
    return grid, rate, mu, sd


def plain_background(ctx):
    """Make the unrelated-line background the whole observed list.

    The state of every run made before the transition-probability file was
    available: the free lines are not separated from the identified ones and
    the levels that have not been found contribute nothing.  Used by
    --plain-background, as the fallback when the calculated level list cannot
    be read, and by the tests, whose synthetic runs have no such distinction
    to make.
    """
    ctx.rho_bg_o = np.asarray(ctx.rho_o, dtype=float).copy()
    ctx.rho_unk_o = np.zeros(len(ctx.rho_o))
    ctx.bg_mu_f = np.asarray(ctx.bg_mu_o, dtype=float).copy()
    ctx.bg_sd_f = np.asarray(ctx.bg_sd_o, dtype=float).copy()
    ctx.unk_mu_o = np.asarray(ctx.bg_mu_o, dtype=float).copy()
    ctx.unk_sd_o = np.asarray(ctx.bg_sd_o, dtype=float).copy()


def free_line_background(ctx, args, log=print):
    """rho_bg(nu) and g_free(ln I | nu): the population of UNRELATED lines.

    The likelihood asks, of a feature that no accepted transition claims, how
    surprising it is to find a line of that brightness at that place if the
    level under test is not there.  Both halves of that question are about the
    lines nobody has identified, not about the observed list as a whole - the
    identified lines are, by construction, the ones already explained, and they
    are also systematically the bright ones.  Two independent estimates of the
    same population are formed and combined:

      empirical - the 1468 recorded lines that carry no accepted transition,
        their local density from the K_RHO_FREE nearest of them and their local
        ln I from the K_BG_FREE nearest;
      modelled - unknown_transition_background() above, the rate and brightness
        the calculation predicts for the transitions of the levels that have
        not been found.

    The rate used is the LARGER of the two, capped by the all-lines density
    (free lines are a subset of the recorded ones, so no estimate of their rate
    may exceed it) and floored at RHO_BG_FLOOR of it.  The brightness density
    used is likewise the larger of the two normal densities.  Taking the larger
    in both places is the conservative choice - a bigger background is a
    smaller ln R - and it is the same device as the two readings of a match: an
    unrelated line is offered whichever account of itself is the better, the
    one the free lines actually recorded suggest or the one the calculation
    predicts for a level still missing.

    Sets ctx.rho_bg_o, ctx.bg_mu_f, ctx.bg_sd_f and ctx.rho_unk_o on the
    observed-line grid.  With --plain-background all three fall back to the
    all-lines values, which is the behaviour of the runs made before this was
    written.
    """
    n = len(ctx.wn_o)
    if getattr(args, 'plain_background', False):
        plain_background(ctx)
        log("  --plain-background: the unrelated-line background is the whole "
            "observed list, as before")
        return
    ctx.rho_unk_o = np.zeros(n)

    free = ctx.n_acc_line == 0
    n_free = int(free.sum())
    if n_free >= K_RHO_FREE:
        wf = ctx.wn_o[free]
        i = np.arange(len(wf))
        hi = np.clip(i + K_RHO_FREE // 2, 0, len(wf) - 1)
        lo = np.clip(i - K_RHO_FREE // 2, 0, len(wf) - 1)
        span = wf[hi] - wf[lo]
        rho_f = np.where(span > 0, K_RHO_FREE / np.where(span > 0, span, 1.0),
                         np.nan)
        rho_emp = np.interp(ctx.wn_o, wf, rho_f)
        mu_f = running(ctx.ln_int_o[free], K_BG_FREE, 'mean')
        sd_f = np.maximum(running(ctx.ln_int_o[free], K_BG_FREE, 'std'), 0.2)
        ctx.bg_mu_f = np.interp(ctx.wn_o, wf, mu_f)
        ctx.bg_sd_f = np.interp(ctx.wn_o, wf, sd_f)
    else:                                                     # pragma: no cover
        rho_emp = ctx.rho_o.copy()
        ctx.bg_mu_f = ctx.bg_mu_o.copy()
        ctx.bg_sd_f = ctx.bg_sd_o.copy()

    # The modelled half depends on the calculated transition list and the
    # partner energies, not on which observed features are claimed, so it
    # survives a release of assignments unchanged and is computed once.  The
    # empirical half above does depend on it, and has just been re-measured.
    if hasattr(ctx, 'unk_bg'):
        unk = ctx.unk_bg
    else:
        unk = unknown_transition_background(ctx, log=log)
        ctx.unk_bg = unk
    if unk is not None:
        grid, rate, mu, sd = unk
        ctx.rho_unk_o = np.interp(ctx.wn_o, grid, rate)
        ok = np.isfinite(mu) & np.isfinite(sd)
        if ok.any():
            ctx.unk_mu_o = np.interp(ctx.wn_o, grid[ok], mu[ok])
            ctx.unk_sd_o = np.maximum(np.interp(ctx.wn_o, grid[ok], sd[ok]),
                                      0.2)
        else:                                                 # pragma: no cover
            ctx.unk_mu_o = ctx.bg_mu_f.copy()
            ctx.unk_sd_o = ctx.bg_sd_f.copy()
    else:                                                     # pragma: no cover
        ctx.unk_mu_o = ctx.bg_mu_f.copy()
        ctx.unk_sd_o = ctx.bg_sd_f.copy()
    bad = ~np.isfinite(ctx.unk_mu_o) | ~np.isfinite(ctx.unk_sd_o)
    ctx.unk_mu_o = np.where(bad, ctx.bg_mu_f, ctx.unk_mu_o)
    ctx.unk_sd_o = np.where(bad, ctx.bg_sd_f, ctx.unk_sd_o)

    rho = np.maximum(rho_emp, ctx.rho_unk_o)
    rho = np.minimum(rho, ctx.rho_o)
    ctx.rho_bg_o = np.maximum(rho, RHO_BG_FLOOR * ctx.rho_o)

    log(f"  unrelated-line background: {n_free} of {n} recorded lines carry "
        f"no accepted transition")
    log(f"    empirical density {np.nanmedian(rho_emp):.4f}, modelled "
        f"{np.median(ctx.rho_unk_o):.4f}, used {np.median(ctx.rho_bg_o):.4f} "
        f"per cm^-1 (median; all lines {np.median(ctx.rho_o):.4f})")
    log(f"    ln I of a free line: mean {np.nanmedian(ctx.bg_mu_f):.2f} "
        f"(all lines {np.nanmedian(ctx.bg_mu_o):.2f}), of a predicted "
        f"unfound-level line {np.nanmedian(ctx.unk_mu_o):.2f}")


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
    # J of every level of the run, for the delta J fingerprint
    ctx.j_of.update({str(k): j_value(v)
                     for k, v in zip(per['level_id'], per['J'])})
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
    # the modelled width of each recorded feature, section 4: this, and not
    # the quoted unc_wn_obs, is what a candidate position is judged against
    ctx.meas_o = meas_sigma(ctx.wn_o, ctx.char_o, ctx.unc_o)
    ctx.w_hfs = read_hfs_widths(log=log)
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
    ctx.sig_loc_o = running(ctx.meas_o, K_UNC, 'median')
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
    ctx.preds_all = ls.load_predictions(set(per['level_id']), ctx.e_final,
                                        with_u=True)
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
    meas = meas_sigma(a['wn_obs'].to_numpy(),
                      a['char'].fillna('').astype(str).to_numpy(),
                      a['unc_wn_obs'].to_numpy(dtype=float))
    v_hfs = hfs_pair(ctx.w_hfs, a['low_id'].astype(str),
                     a['upp_id'].astype(str))
    # k(n) is what is left over once the width of the feature is accounted
    # for, hyperfine part included: measured without it, it would charge the
    # blending for a width that belongs to the level
    t = ((a['wn_obs'] - a['rwn']).to_numpy() / np.sqrt(meas ** 2 + v_hfs))
    nb = a['n_accepted'].fillna(1).to_numpy(dtype=float)
    ctx.k_n = fit_k_blend(t, nb, log)
    kn = np.array([ctx.k_n.get(int(x), 1.0) for x in nb])
    sigma = np.sqrt((kn * meas) ** 2 + v_hfs)
    res = (a['wn_obs'] - a['rwn']).to_numpy()
    z = res / sigma
    log(f"  the width model: rms of residual/sigma over the {len(a)} accepted "
        f"lines = {math.sqrt(float(np.mean(z ** 2))):.3f}; a class that comes "
        f"out far from 1 is a class whose width is mis-stated")
    code = np.array([normalize_char(c) for c in a['char'].fillna('')])
    era = era_of(a['wn_obs'].to_numpy())
    seen = sorted(set(zip(era, code)))
    log("    " + "  ".join(
        f"{e[2:]}/{c if c else '-'}:{int(m.sum())}@"
        f"{math.sqrt(float(np.mean(z[m] ** 2))):.2f}"
        for e, c in seen
        for m in [(era == e) & (code == c)] if m.sum() >= 10))
    # the all-lines density here, not the unrelated one of section 4a: eta is
    # the rate at which a GENUINE line sits at an anomalous position, and the
    # background it is mixed against is every line it could have been taken for
    rho_a = np.interp(a['wn_obs'].to_numpy(), ctx.wn_o, ctx.rho_o)
    ctx.eta = fit_eta((a['wn_obs'] - a['rwn']).to_numpy(), sigma, rho_a, log)

    f = np.array([ls.scale_bias_factor(w, ctx.bias) for w in a['wn_obs']])
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.log(a['obs_intens'].to_numpy(dtype=float)
                   * a['BF'].fillna(1.0).to_numpy(dtype=float)
                   / (a['calc_intens'].to_numpy(dtype=float) * f))
    lids = np.where(a['low_id'].astype(str) != '', a['low_id'].astype(str),
                    a['upp_id'].astype(str))
    pairs = [tuple(sorted((str(lo), str(up))))
             for lo, up in zip(a['low_id'], a['upp_id'])]
    pa = ctx.preds_all
    by_pair = {tuple(sorted((str(x.lo_id), str(x.up_id)))): x
               for x in pa.itertuples(index=False)}
    u_line = (a['u_calc'].to_numpy(dtype=float)
              if 'u_calc' in a.columns else None)
    if u_line is None:
        u_line = np.array([getattr(by_pair.get(k), 'u_calc', np.nan)
                           for k in pairs], dtype=float)
    # the two coordinates of the width model and the detection cut-off
    # (fit_intensity_width): log10 P, the predicted plate intensity over the
    # threshold; log10 S, the calculated line strength; and rthr, the residual
    # a line exactly at the threshold would have
    s_line = np.array([getattr(by_pair.get(k), 'S', np.nan) for k in pairs],
                      dtype=float)
    thr = (noise_threshold_vec(a['wn_obs'].to_numpy(dtype=float), ctx.calib)
           if ctx.calib is not None else np.full(len(a), np.nan))
    i_c = a['calc_intens'].to_numpy(dtype=float)
    bf = a['BF'].fillna(1.0).to_numpy(dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        lp_line = np.log10(i_c * f / thr)
        ls_line = np.log10(s_line)
        rthr = np.log(thr * bf / (i_c * f))
    # the fit is made on the lines whose intensity is theirs alone: one
    # accepted transition, not flagged bl
    single = ((a['n_accepted'].fillna(1).to_numpy(dtype=float) <= 1)
              & ~a['char'].fillna('').astype(str).str.contains('bl')
              .to_numpy())
    ctx.s, ctx.s_L, ctx.width = fit_intensity_scatter(
        r, lids, log, u=u_line, lp=lp_line, ls10=ls_line, rthr=rthr,
        fit=single)

    ctx.u_M = read_partner_uncertainties()
    if ctx.u_M:
        v = np.array(list(ctx.u_M.values()))
        log(f"  partner energy uncertainty D1: median {np.median(v):.4f}, "
            f"max {np.max(v):.4f} cm^-1 over {len(v)} levels")
    else:
        log("  LOPT_output_levels.txt not found - the partner energy "
            "contributes nothing to sigma_t")

    # --- the background of lines that belong to nothing yet known --------
    free_line_background(ctx, args, log)

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
    uc = (p['u_calc'].to_numpy(dtype=float) if 'u_calc' in p.columns
          else np.full(len(p), np.nan))
    with np.errstate(divide='ignore', invalid='ignore'):
        sl = (np.log10(p['S'].to_numpy(dtype=float)) if 'S' in p.columns
              else np.full(len(p), np.nan))
    rows = {}
    for i in range(len(p)):
        rows.setdefault(up[i], []).append((lo[i], +1.0, ip[i], uc[i], sl[i]))
        rows.setdefault(lo[i], []).append((up[i], -1.0, ip[i], uc[i], sl[i]))
    for lid, rr in rows.items():
        partner = np.array([x[0] for x in rr])
        sign = np.array([x[1] for x in rr])
        i_pred = np.array([x[2] for x in rr])
        u_calc = np.array([x[3] for x in rr])
        ls10 = np.array([x[4] for x in rr])
        e_m = np.array([ctx.e_final.get(x, np.nan) for x in partner])
        u_m = np.array([ctx.u_M.get(x, 0.0) for x in partner])
        n_p = np.array([ctx.n_acc_level.get(x, 0) for x in partner])
        sh = np.array([ctx.shared.get((x, lid), 0) + ctx.shared.get((lid, x), 0)
                       for x in partner])
        out[lid] = dict(partner=partner, sign=sign, i_pred=i_pred, e_m=e_m,
                        u_m=u_m, u_calc=u_calc, ls10=ls10,
                        degenerate=(n_p - sh) <= 0)
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
    # what the intensity width of each predicted transition is made of: its
    # own calculated-intensity uncertainty and line strength; the plate
    # intensity P depends on where the candidate energy puts the line, and
    # enters below
    uc_pred = g.get('u_calc', np.full(len(keep), np.nan))[keep]
    ls_pred = g.get('ls10', np.full(len(keep), np.nan))[keep]
    partner = g['partner'][keep]
    own = ctx.own_claim.get(level_id, {})
    # the hyperfine width of the pair of levels the transition joins: the
    # level being scanned contributes its own to every one of its lines
    v_hfs = (ctx.w_hfs.get(str(level_id), 0.0) ** 2
             + np.array([ctx.w_hfs.get(str(x), 0.0) ** 2 for x in partner]))

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

        # The window is set by whichever is worse, the width typical of the
        # neighbourhood or the one the candidate line itself is modelled at: a
        # feature flagged complex where its neighbours are plain must not be
        # ruled out for lying half its own width away.  K_WORST covers the
        # largest k(n) the run produces, and the hyperfine widths of the two
        # levels widen it further.  Widening the window costs nothing: a line
        # far out has N(d; 0, sigma)/rho far below 1 and contributes
        # ln(1 - p), the same as an absence.
        unc_near = ctx.meas_o[pick]
        unc_near = np.where(np.isfinite(unc_near) & (unc_near > 0), unc_near,
                            sig_loc)
        sig_w = np.sqrt(np.maximum(sig_loc, K_WORST * unc_near) ** 2
                        + (u_m ** 2 + v_hfs)[:, None])
        w = N_SIGMA * sig_w
        matched = (np.abs(d) <= w) & ok

        p = p_t * (1.0 - ctx.eta)
        r_t = np.where(ok, 1.0 - p * (1.0 - Q_OUT), 1.0)

        if matched.any():
            mi = pick[matched]
            unc = ctx.meas_o[mi]
            unc = np.where(np.isfinite(unc) & (unc > 0), unc,
                           sig_loc[matched])
            # C: what other levels already claim on that feature
            own_here = np.array([own.get(int(x), 0.0) for x in mi])
            claimed = np.maximum(ctx.claimed_tot[mi] - own_here, 0.0)
            n_blend = ctx.n_acc_line[mi] - (own_here > 0).astype(int) + 1
            kn = np.array([ctx.k_n.get(int(x), ctx.k_n.get(1, 1.0))
                           for x in n_blend])
            u_here = np.repeat(u_m, nu_s.shape[1]).reshape(nu_s.shape)[matched]
            v_here = np.repeat(v_hfs,
                               nu_s.shape[1]).reshape(nu_s.shape)[matched]
            sig = np.sqrt((kn * unc) ** 2 + u_here ** 2 + v_here)

            # the intensity ratio
            i_here = np.repeat(i_pred, nu_s.shape[1]).reshape(nu_s.shape)[matched]
            f = scale_factor_vec(nu_s[matched], ctx.bias)
            thr = (noise_threshold_vec(nu_s[matched], ctx.calib)
                   if ctx.calib is not None
                   else np.full(int(matched.sum()), np.nan))
            with np.errstate(divide='ignore', invalid='ignore'):
                lp_here = np.log10(i_here * f / thr)
            s_here = sigma_intensity(
                ctx,
                np.repeat(uc_pred, nu_s.shape[1]).reshape(nu_s.shape)[matched],
                lp_here,
                np.repeat(ls_pred, nu_s.shape[1]).reshape(nu_s.shape)[matched])
            x = ctx.ln_int_o[mi]
            ln_bg = ln_norm(x, ctx.bg_mu_o[mi], ctx.bg_sd_o[mi])
            free = claimed <= 0
            # g for a FREE feature is the brightness distribution of the
            # unrelated population - the recorded lines that carry no accepted
            # transition, and the lines the calculation predicts for the levels
            # nobody has found - not that of the observed list as a whole,
            # which is dominated by the identified and therefore bright lines.
            # The better of the two accounts is offered, and the result is
            # never allowed more than LN_BG_DROP below the all-lines value:
            # both estimates extrapolate where free lines are sparse, and an
            # extrapolated density must not buy unbounded credit.
            ln_bg_f = np.maximum(
                np.maximum(ln_norm(x, ctx.bg_mu_f[mi], ctx.bg_sd_f[mi]),
                           ln_norm(x, ctx.unk_mu_o[mi], ctx.unk_sd_o[mi])),
                ln_bg - LN_BG_DROP)
            ln_bg_f = np.where(np.isfinite(ln_bg_f), ln_bg_f, ln_bg)
            with np.errstate(divide='ignore', invalid='ignore'):
                ln_p1 = ln_intensity(x, (claimed + i_here) * f, s_here, ln_bg)
                ln_p0 = np.where(free, ln_bg_f,
                                 ln_intensity(x, claimed * f, s_here, ln_bg))
            ln_g = np.where(np.isfinite(ln_p1 - ln_p0), ln_p1 - ln_p0, 0.0)

            dens = (np.exp(-0.5 * (d[matched] / sig) ** 2)
                    / (sig * math.sqrt(2.0 * math.pi)))
            # TWO READINGS OF A MATCH ON A FREE FEATURE, and the row is worth
            # whichever is the better of them.
            #
            #   the feature IS the transition:  its position is drawn from
            #     N(d; 0, sigma) against the local line density rho, and its
            #     brightness from N(x; ln I_t f, s) against the local
            #     distribution g of ln I_obs.
            #   the transition is hidden IN the feature:  a line is there in
            #     any case and carries light nobody has identified, so under
            #     H0 the feature is present with certainty - the density of a
            #     line known to be somewhere in the window, 1/(2 W) - and its
            #     brightness says nothing either way.
            #
            # The second reading is what the one-sided floor of ln_intensity
            # was reaching for, and taking the floor without also replacing
            # rho was the error: a prediction of I = 0.2 landing on a feature
            # of 13372 kept the whole of the positional credit, ln R_t = +1.4
            # for a transition with a 1.5 per cent chance of having been
            # recorded at all.  Rows of that kind - the strong observed line
            # assigned to the very weak transition, the assignment
            # classify_lines refuses and the analyst leaves free for a better
            # one - carried four levels of the audit to alternate positions.
            # Taking the better of the two readings rather than switching on
            # the floor is what keeps a genuine line: one under-predicted by
            # a factor of thirty is still far better explained as the
            # transition than as a coincidence, and it goes on being read
            # that way.
            # rho for the "it IS the transition" reading of a FREE feature is
            # the rate of UNRELATED lines, not of all recorded lines: a feature
            # that no accepted transition claims cannot have been drawn from
            # the identified ones.  A claimed feature keeps the all-lines rate,
            # because there the question is the blend one and rho_t is replaced
            # by 1/(2W) in any case.
            rho_here = np.where(free, ctx.rho_bg_o[mi], ctx.rho_o[mi])
            with np.errstate(divide='ignore', invalid='ignore'):
                ln_here = (np.log(dens / rho_here)
                           + ln_norm(x, np.log(np.maximum(i_here * f, 1e-300)),
                                     s_here) - np.where(free, ln_bg_f, ln_bg))
                ln_hidden = np.log(dens * 2.0 * w[matched])
                ln_blend = np.log(dens * 2.0 * w[matched]) + ln_g
                # a transition can only hide in a feature that has its light:
                # a feature fainter than predicted by more than the scatter
                # is read as the transition or as nothing
                ln_hidden = np.where(
                    x >= np.log(np.maximum(i_here * f, 1e-300)) - s_here,
                    ln_hidden, -np.inf)
            over = free & (ln_hidden >= ln_here)
            ln_t = np.where(free, np.maximum(ln_here, ln_hidden), ln_blend)
            ln_t = np.where(np.isfinite(ln_t), ln_t, -50.0)
            pm = p[matched]
            genuine = pm * np.exp(np.clip(ln_t, -50.0, 50.0))
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
            weight[matched] = w_gen * (free & ~over)
            with np.errstate(divide='ignore', invalid='ignore'):
                res = np.zeros(matched.shape)
                res[matched] = x - np.log(np.maximum((claimed + i_here) * f,
                                                     1e-300))
            res = np.where(np.isfinite(res), res, 0.0)
            # with a width of its own per row the correction is written in
            # the two sufficient statistics A = sum w/sigma^2 and
            # B = sum w r/sigma^2, which reduce to n/s^2 and (sum r)/s^2 when
            # every sigma is the same s - the form of section 3
            sv = np.zeros(matched.shape)
            sv[matched] = s_here ** 2
            sv = np.where(sv > 0, sv, ctx.s ** 2)
            sl2 = ctx.s_L ** 2
            A = (weight / sv).sum(axis=0)
            B = (weight * res / sv).sum(axis=0)
            block_total += (-0.5 * np.log(1.0 + sl2 * A)
                            + sl2 * B ** 2 / (2.0 * (1.0 + sl2 * A)))

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
            over_col = np.zeros(len(partner), dtype=bool)
            if matched.any():
                claim_col[matched[:, 0]] = claimed
                gain_col[matched[:, 0]] = ln_g
                over_col[matched[:, 0]] = over
            tab['C'] = claim_col
            tab['ln_G'] = gain_col
            tab['over'] = over_col
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


def scan_interval(e_adopted, e_calc, window, k=WINDOW_K):
    """The energies the scan covers: k (three) configuration windows either
    side of where the calculation puts the level, always including the adopted
    position with five wavenumbers to spare.  Fifty wavenumbers either side of
    the adopted position when the calculation says nothing."""
    if not np.isfinite(e_calc) or not np.isfinite(window) or window <= 0:
        lo, hi = e_adopted - 50.0, e_adopted + 50.0
    else:
        lo, hi = e_calc - k * window, e_calc + k * window
    return min(lo, e_adopted - 5.0), max(hi, e_adopted + 5.0)


def scan_level(ctx, level_id, e_adopted, e_calc, window, step=GRID_STEP,
               alt_drop=ALT_DROP):
    """Scan ln R across the interval the calculation allows.

    Returns a dict with the adopted value, the alternate maxima ranked by
    ln R, and the grid itself.
    """
    lo, hi = scan_interval(e_adopted, e_calc, window)
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
# The registry of settled positions
# ---------------------------------------------------------------------------
def run_digest(ctx):
    """A hash of the ingredients of ln R that no single level owns.

    The list of measured lines and the constants fitted to the run enter
    every level's scan, so when they change no registry entry can be trusted.
    The constants are rounded to three decimals first: they are re-fitted on
    every run and move in the last figures even when nothing of substance has
    changed, and a registry that threw itself away on that would never be
    used.
    """
    h = hashlib.sha1()
    for a in (ctx.wn_o, ctx.int_o, ctx.unc_o):
        h.update(np.ascontiguousarray(a, dtype=float).tobytes())
    h.update('\x00'.join(ctx.char_o).encode('utf-8'))
    w = getattr(ctx, 'width', None) or flat_width(ctx.s)
    h.update(('%.3f|%.3f|%.3f|%.3f|%.3f|%.3f|%.2f|%.2f|%.3f'
              % (ctx.eta, ctx.s, ctx.s_L, w.a0, w.a1, w.c, w.p_ref,
                 w.s_ref, w.k_min)).encode())
    # the width model: a changed constant changes every level's sigma
    h.update(repr((UNC_FLOOR, ERA_SPLIT, sorted(DLAM_ERA.items()),
                   sorted(CHAR_DLAM.items()),
                   sorted(CHAR_DWN.items()))).encode())
    # the unrelated-line background: it enters every matched free row
    for a in (ctx.rho_bg_o, ctx.bg_mu_f, ctx.bg_sd_f, ctx.unk_mu_o,
              ctx.unk_sd_o):
        h.update(np.round(np.asarray(a, dtype=float), 3).tobytes())
    for d in (ctx.k_n, ctx.w_hfs, ctx.bias):
        h.update(repr(sorted((str(k), np.round(np.asarray(v, dtype=float), 3)
                              .tolist())
                             for k, v in d.items())).encode())
    return h.hexdigest()[:16]


def accepted_index(ctx):
    """The accepted assignments, sorted by wavenumber, each as a short string.

    A level's scan sees the rest of the run only through the assignments that
    sit on the lines its own predictions can reach: those are what set C, the
    intensity another transition already claims on a feature, and the number
    of components blended into it.  Keeping them in wavenumber order lets the
    fingerprint of one level pick out exactly the ones that could affect it,
    so that accepting a relocation invalidates the neighbours of that level
    and nothing else.
    """
    a = ctx.acc
    o = np.argsort(a['wn_obs'].to_numpy(dtype=float))
    wn = a['wn_obs'].to_numpy(dtype=float)[o]
    key = np.array(['%.4f|%s|%s|%.4g' % (w, lo, up, ic) for w, lo, up, ic in
                    zip(wn, a['low_id'].to_numpy()[o],
                        a['upp_id'].to_numpy()[o],
                        a['calc_intens'].fillna(0.0)
                        .to_numpy(dtype=float)[o])], dtype=object)
    return wn, key


def scan_fingerprint(ctx, level_id, e_adopted, e_calc, window, step, alt_drop,
                     glob, acc_wn, acc_key):
    """A hash of everything the scan of this one level depends on.

    Two runs that give this the same value give the same scan: the interval
    scanned and its step, the level's own energy, every prediction it makes
    (partner, partner energy, partner uncertainty, predicted intensity), and
    the accepted assignments lying anywhere the predictions can reach.  A
    level is taken from the registry only when its fingerprint is unchanged
    AND ln R at the adopted position, which is computed for every level
    anyway, still comes out what the registry recorded.
    """
    lo, hi = scan_interval(e_adopted, e_calc, window)
    h = hashlib.sha1(glob.encode())
    h.update(('%s|%.4f|%.3f|%.3f|%.4f|%.2f'
              % (level_id, e_adopted, lo, hi, step, alt_drop)).encode())
    g = ctx.by_level.get(level_id)
    if g is None:
        return h.hexdigest()[:16]
    keep = ~g['degenerate'] & np.isfinite(g['e_m']) & (g['i_pred'] > 0)
    e_m, sign = g['e_m'][keep], g['sign'][keep]
    uc = g.get('u_calc', np.full(len(keep), np.nan))[keep]
    h.update(repr(sorted(
        '%s|%.3f|%.3f|%.4g|%.3f' % (p, e, u, i, c if np.isfinite(c) else -1.0)
        for p, e, u, i, c in
        zip(g['partner'][keep], e_m, g['u_m'][keep],
            g['i_pred'][keep], uc))).encode())
    # the wavenumbers this level's predictions sweep as it moves over lo..hi
    a = sign * (lo - e_m)
    b = sign * (hi - e_m)
    reach = np.zeros(len(acc_wn), dtype=bool)
    for x, y in zip(np.minimum(a, b) - REACH_PAD, np.maximum(a, b) + REACH_PAD):
        reach[np.searchsorted(acc_wn, x):np.searchsorted(acc_wn, y)] = True
    h.update('\x00'.join(acc_key[reach]).encode('utf-8'))
    return h.hexdigest()[:16]


def read_firm(path):
    """{level_id: row} of the registry, or {} when there is no registry."""
    if not path or not os.path.exists(path):
        return {}
    t = pd.read_csv(path, dtype={'level_id': str, 'fingerprint': str})
    return {str(r['level_id']): r for _, r in t.iterrows()}


def write_firm(path, rows, log=print):
    """Write the registry, newest verdict per level, in energy order."""
    cols = ['level_id', 'E', 'ln_R', 'scan_width', 'fingerprint']
    t = pd.DataFrame(rows, columns=cols)
    # an empty registry is still written: the alternative is to leave the
    # entries of a previous run standing for levels this run has just found
    # alternates for
    t = t.sort_values('E').reset_index(drop=True)
    t['E'] = t['E'].round(4)
    t['ln_R'] = t['ln_R'].round(3)
    t['scan_width'] = t['scan_width'].round(1)
    output_files.require_writable([path], 'registry file')
    t.to_csv(path, index=False, lineterminator='\n')
    log(f"registry: {len(t)} settled levels written to {path}")


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
                             ln_R_miss=0.0, sum_lnG=0.0, ln_J=0.0))
            continue
        m = tab['matched'].to_numpy(dtype=bool)
        seen, real, _ = match_kinds(tab)
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
                         sum_lnG=float(tab['ln_G'][m].sum()),
                         ln_J=ln_j_pattern(
                             tab, ctx.j_of.get(str(lid), float('nan')),
                             ctx.j_of, seen, real)['ln_J']))
    return pd.DataFrame(rows)


def iden2_identities(per, en, tol=li.ENLEV_MATCH_TOL):
    """{level_id: (index in enlev.dat, label)} for the levels of the run.

    The two lists share no key, so they are tied by energy exactly as
    level_interchange.attach_identities ties them: a level of the run is the
    starred - experimentally found - row of enlev.dat whose observed energy is
    nearest, if that row is within `tol`.  The index is the running number in
    column 1 of enlev.dat, which is the level id IDEN2 shows on the screen,
    and it is what a detail table has to print if the table is to be read
    beside IDEN2 without looking every partner up by hand.
    """
    known = en[en['known']].sort_values('E_obs').reset_index(drop=True)
    e_arr = known['E_obs'].to_numpy(dtype=float)
    out = {}
    for lid, e in zip(per['level_id'], per['E_final']):
        if not len(e_arr):
            break
        i = int(np.argmin(np.abs(e_arr - float(e))))
        if abs(e_arr[i] - float(e)) <= tol:
            out[str(lid)] = (int(known['idx'][i]), str(known['label'][i]))
    return out


def detail_energy(ctx, level_id, at, alt_drop=ALT_DROP):
    """The energy --detail is to be evaluated at, and how it was arrived at.

    ``at`` is None for the adopted position, the word ``alt`` for the best
    alternate the scan finds - which is the position an audit row proposes,
    and the one there is any point inspecting when the report suggests a move
    - or an energy in cm^-1.
    """
    e_adopted = ctx.e_final[level_id]
    if at is None:
        return e_adopted, 'the adopted position'
    if str(at).strip().lower() != 'alt':
        return float(at), 'the energy asked for'
    en = li.read_enlev()
    win = li.configuration_windows(en)
    lv, _ = li.attach_identities(ctx.per, en, win, row_of_id=li.id_rows())
    e_calc = dict(zip(lv['level_id'], lv['E_calc']))
    w_of = dict(zip(lv['level_id'], lv['W']))
    r = scan_level(ctx, level_id, e_adopted, e_calc.get(level_id, np.nan),
                   w_of.get(level_id, np.nan), alt_drop=alt_drop)
    if not r['alternates']:
        return e_adopted, 'the adopted position - the scan finds no alternate'
    e_a = r['alternates'][0][0]
    return e_a, ('the best of the %d alternate positions, %+.4f cm^-1 from '
                 'the adopted one' % (len(r['alternates']), e_a - e_adopted))


def print_detail(ctx, level_id, at=None, show_all=False,
                 alt_drop=ALT_DROP):
    """The per-transition table for one level, at one energy.

    Only the rows worth looking at are printed: the observable predictions
    (P_obs >= P_SEEN) and any row that found a line, whether observable or
    not.  The rest - and there are usually a hundred or more of them, faint
    predictions to distant partners that were never going to be recorded -
    are counted and left out, since they say nothing about the position and
    burying the handful that matter is what makes the full table useless.
    ``--detail-all`` prints them anyway.

    The `what` column is the one to read against the audit row: `free` marks
    the real matches (support()) on features no accepted transition claims,
    which are exactly the n_free of the audit; `blend` a real match on a
    feature that is already explained; `bright` a match on a feature so much
    stronger than the prediction that an unidentified line explains it better
    - it would have been recorded there in any case, so it counts for
    nothing; `poor` a line in the matching window that is not a match, too
    far out or of the wrong brightness (its ln R_t is not positive); `faint`
    a match to a prediction that could not have been recorded, which is a
    coincidence; `-` an absence.  free + blend is the n_match of --unknown.
    """
    e, how = detail_energy(ctx, level_id, at, alt_drop)
    v, tab = ln_ratio(ctx, level_id, np.array([e]), detail=True)
    ident = iden2_identities(ctx.per, li.read_enlev())
    own = ident.get(str(level_id))
    print(f"\nlevel {level_id} at E = {e:.4f} cm^-1 - {how}")
    if own:
        print(f"  IDEN2 level {own[0]}  {own[1]}")
    print(f"ln R = {float(v[0]):+.2f} over {0 if tab is None else len(tab)} "
          f"predicted transitions")
    if tab is None:
        return
    g = ctx.by_level[level_id]
    print(f"  {int(g['degenerate'].sum())} predictions dropped: the partner "
          f"has no accepted line of its own")
    m = tab['matched'].to_numpy(dtype=bool)
    o = tab['over'].to_numpy(dtype=bool)
    r = tab['ln_R'].to_numpy(dtype=float)
    seen, real, free = match_kinds(tab)
    what = np.where(~m, '-',
                    np.where(o, 'bright',
                             np.where(~seen, 'faint',
                                      np.where(~real, 'poor',
                                               np.where(free, 'free',
                                                        'blend')))))
    print(f"  {int(m.sum())} matched, contributing {tab['ln_R'][m].sum():+.2f}; "
          f"{int((~m).sum())} absent, costing {tab['ln_R'][~m].sum():+.2f}")
    print(f"  {int(seen.sum())} of the predictions could have been recorded "
          f"(P_obs >= {P_SEEN:g}); {int(real.sum())} of them match a "
          f"recorded line, {int(free.sum())} of those free and worth "
          f"{r[free].sum():+.2f}")
    j_level = ctx.j_of.get(str(level_id), float('nan'))
    jf = ln_j_pattern(tab, j_level, ctx.j_of, seen, real)
    if jf['classes']:
        print(f"  the delta J fingerprint, over those observable "
              f"predictions (J of this level = {j_level:g}):")
        for cl_ in jf['classes']:
            tag_j = '%+d' % cl_['dJ'] if cl_['dJ'] else ' 0'
            print(f"    dJ = {tag_j}:  {cl_['n_match']} of "
                  f"{cl_['n_obs']} found, {cl_['sum_P']:.2f} expected"
                  f"   detection rate {cl_['lam']:.2f}"
                  f"   worth {cl_['gain']:.2f}"
                  + ('' if cl_['counted'] else
                     '  (one prediction: not a pattern, not counted)'))
        if jf['n_no_j']:
            print(f"    {jf['n_no_j']} observable prediction(s) left out: "
                  f"the partner has no J")
        print(f"  ln_J = {jf['ln_J']:.2f}"
              + ("  - the matches are spread over the dJ classes as the "
                 "predicted intensities say they should be"
                 if jf['ln_J'] <= 0.0 else
                 "  - nothing much: no dJ class is missing by more than the "
                 "line list alone would explain"
                 if jf['ln_J'] < J_MARK else
                 "  - a dJ class that should have been recorded is not "
                 "there.  Read against the same number at the other "
                 "candidate energies before concluding anything about the "
                 "position: a level whose calculated gA divides its strength "
                 "wrongly among the branches scores this at EVERY energy, "
                 "and only the DIFFERENCE between two candidate energies "
                 "enters ln R"))
    show = np.ones(len(tab), dtype=bool) if show_all else (seen | m)
    hidden = int((~show).sum())
    if hidden:
        print(f"  {hidden} faint predictions that found no line are not "
              f"listed (--detail-all lists them)")
    print()
    j_p = np.array([ctx.j_of.get(str(x), float('nan'))
                    for x in tab['partner']])
    dj = (np.round(j_p - j_level) if np.isfinite(j_level)
          else np.full(len(tab), np.nan))
    hdr = (f"{'partner':<14}{'IDEN2':>6} {'label':<12}{'E_partner':>12}"
           f"{'dJ':>4}{'nu':>12}{'I_pred':>10}{'P_obs':>7}"
           f"{'wn_obs':>12}{'d':>8}"
           f"{'I_obs':>10}{'C':>10}{'lnG':>7}{'lnR':>8}  what")
    print(hdr)
    print('-' * len(hdr))
    for (_, r), tag, keep, dd in zip(tab.iterrows(), what, show, dj):
        if not keep:
            continue
        pid = str(r['partner'])
        idx, label = ident.get(pid, ('', ''))
        e_m = ctx.e_final.get(pid, float('nan'))
        head = (f"{pid:<14}{idx:>6} {label:<12}{e_m:>12.3f}"
                f"{('%+d' % dd) if np.isfinite(dd) else '-':>4}"
                f"{r['nu']:>12.3f}{r['I_pred']:>10.1f}{r['P_obs']:>7.3f}")
        if r['matched']:
            print(f"{head}{r['wn_obs']:>12.3f}{r['d']:>8.3f}"
                  f"{r['I_obs']:>10.1f}{r['C']:>10.1f}{r['ln_G']:>7.2f}"
                  f"{r['ln_R']:>8.2f}  {tag}")
        else:
            print(f"{head}{'-':>12}{'-':>8}{'-':>10}{'-':>10}"
                  f"{'-':>7}{r['ln_R']:>8.2f}  {tag}")


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
REFIT_DE = 2.0          # an alternate this close may be the same level
REFIT_KEEP = 0.5        # ... but only if it keeps more than this share of
                        # the recorded lines the level is assigned now
REFIT_LIGHT = 0.5       # ... and more than this share of their light, so
                        # that dropping the level's brightest branch is never
                        # a refit however many faint lines survive
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



def match_kinds(tab):
    """seen, real and free of a detail table - the one definition of each.

    ``seen``  the prediction could have been recorded at all (P_obs >= P_SEEN)
    ``real``  and it matches a recorded line that counts as evidence for the
              level: within FREE_SIGMA of where it is predicted, not a row the
              `bright` reading won, and with its own ln R_t positive
    ``free``  and the feature is claimed by no other accepted transition

    support(), print_detail() and the J fingerprint all ask the same question
    of the same table and must get the same answer, so they all ask it here.
    """
    m = tab['matched'].to_numpy(dtype=bool)
    seen = tab['P_obs'].to_numpy(dtype=float) >= P_SEEN
    d = np.nan_to_num(tab['d'].to_numpy(dtype=float), nan=np.inf)
    c = np.nan_to_num(tab['C'].to_numpy(dtype=float), nan=-1.0)
    o = tab['over'].to_numpy(dtype=bool)
    w = tab['W'].to_numpy(dtype=float)
    r = tab['ln_R'].to_numpy(dtype=float)
    real = (seen & m & ~o & (r > 0.0)
            & (np.abs(d) <= FREE_SIGMA * w / N_SIGMA))
    return seen, real, real & (c == 0.0)


# ---------------------------------------------------------------------------
# The J fingerprint
# ---------------------------------------------------------------------------
# For an electric dipole transition J changes by 0 or 1, so WHICH partners a
# level is seen with is a fingerprint of its own J.  ln R cannot read it.  The
# two energies being compared are the same level, so they predict the same
# transitions to the same partners with the same gA and the same intensities;
# only the wavenumbers shift.  The J information therefore never appears as a
# difference in predicted intensity - it appears as a CORRELATION AMONG THE
# ABSENCES, and ln R multiplies the absences as if they were independent, which
# is exactly the assumption that makes a correlation cost nothing.
#
# To charge for it there has to be a rival hypothesis under which the
# correlation is expected, and the physical one is: the lines at this position
# are not this level's, but those of a level of some OTHER J - in which case a
# whole delta J class of the predictions was never going to be there.  Give
# each class its own detection multiplier lambda_g, so that a prediction of
# that class is recorded with probability lambda_g * P_obs instead of P_obs,
# and fit the lambda to the found/missing pattern.  ln_J is how much better the
# pattern is explained that way.
#
# NOTHING HERE IS A SELECTION RULE.  A prediction enters only if it could have
# been recorded at all (P_obs >= P_SEEN) and it enters weighted by its own
# P_obs: a missing partner at P_obs = 0.29 contributes ln(1 - 0.29) = -0.34 and
# one at 0.72 contributes -1.27.  A class of transitions that is simply too
# weak to see reaches lambda = 0 at almost no gain in likelihood and cannot
# accuse a position; only a class that SHOULD have been recorded and was not,
# again and again, can.  That is the whole point of doing this with P_obs
# rather than with delta J alone.
J_CLASS_COST = 1.0   # nats charged for each detection rate the test fits,
                     # so that a class cannot buy credit by being small
J_MARK = 1.0         # ln_J above which the pattern is worth remarking on
J_MIN_CLASS = 2      # observable predictions a class needs before its
                     # absences can be read as a pattern rather than as the
                     # one absence ln R has already charged for


def fit_detection_rate(p, found):
    """The detection multiplier of one delta J class, and what it gains.

    ``p`` is P_obs of each observable prediction of the class and ``found``
    says which of them match a recorded line.  Under the plain model each is
    recorded with probability p; under the rival it is lambda * p, with
    lambda in [0, 1] - a class can be missed more often than predicted, never
    more often found, since a surplus of matches is not evidence about J.

        ln L(lambda) = n ln lambda + sum_missing ln(1 - lambda p) + const

    Returns ``(lambda_hat, ln L(lambda_hat) - ln L(1))``, the second never
    negative.
    """
    p = np.clip(np.asarray(p, dtype=float), 0.0, 0.999)
    found = np.asarray(found, dtype=bool)
    if p.size == 0:
        return 1.0, 0.0
    q = p[~found]                      # the ones that were not recorded
    n = int(found.sum())
    base = float(np.sum(np.log1p(-q))) if q.size else 0.0   # ln L(1)

    def ln_l(lam):
        if lam <= 0.0:
            return -np.inf if n else float(np.sum(np.log1p(-0.0 * q)))
        return n * math.log(lam) + float(np.sum(np.log1p(-lam * q)))

    if q.size == 0:                    # every one of them was found
        return 1.0, 0.0
    if n == 0:                         # none was: the class is extinguished
        return 0.0, -base
    # d/dlambda of ln L, decreasing in lambda; lambda_hat = 1 when it is still
    # positive there, which is a class found as often as predicted or oftener
    def slope(lam):
        return n / lam - float(np.sum(q / (1.0 - lam * q)))
    if slope(1.0) >= 0.0:
        return 1.0, 0.0
    lo, hi = 1e-9, 1.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if slope(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    lam = 0.5 * (lo + hi)
    return lam, max(ln_l(lam) - base, 0.0)


def ln_j_pattern(tab, j_level, j_of, seen=None, real=None):
    """How much better the pattern of absences fits a level of another J.

    Groups the observable predictions by delta J = J(partner) - J(level) and
    fits fit_detection_rate to each.  ``ln_J`` is what the rival hypothesis
    wins, after J_CLASS_COST is charged for every rate fitted, floored at
    zero: a position with ln_J = 0 says the recorded lines are spread over the
    delta J classes as the intensities predict, and a large ln_J says one
    class is systematically absent, which is what a level of a different J
    sitting there would look like.

    Partners whose J cannot be read are left out of the test (they are still
    in ln R), and so is every prediction too faint to have been recorded.
    Returns a dict with ln_J and one row per class.
    """
    out = dict(ln_J=0.0, classes=[], n_no_j=0)
    if tab is None or not np.isfinite(j_level):
        return out
    if seen is None or real is None:
        seen, real, _ = match_kinds(tab)
    jp = np.array([j_of.get(str(x), float('nan')) for x in tab['partner']])
    p = tab['P_obs'].to_numpy(dtype=float)
    ok = seen & np.isfinite(jp)
    out['n_no_j'] = int((seen & ~np.isfinite(jp)).sum())
    dj = np.round(jp - j_level)
    total = 0.0
    for d in (-1.0, 0.0, 1.0):
        m = ok & (dj == d)
        if not m.any():
            continue
        lam, gain = fit_detection_rate(p[m], real[m])
        # A CLASS OF ONE IS NOT A PATTERN.  Its whole content is that one
        # prediction was missed, and ln R has charged ln(1 - P_obs) for that
        # already; reading it as a fingerprint of J would be charging the
        # same absence twice under another name.  The fingerprint is a
        # CORRELATION among absences, and a correlation needs two.
        counted = int(m.sum()) >= J_MIN_CLASS
        if counted:
            total += max(gain - J_CLASS_COST, 0.0)
        out['classes'].append(dict(dJ=int(d), n_obs=int(m.sum()),
                                   sum_P=float(p[m].sum()),
                                   n_match=int(real[m].sum()),
                                   lam=float(lam), gain=float(gain),
                                   counted=counted))
    out['ln_J'] = float(total)
    return out


def fold_j(ln_r, ln_j):
    """ln R of rival positions of ONE level, with the fingerprint folded in.

    THE ONLY PLACE ln_J IS ALLOWED TO CHANGE ln R, and it changes it only
    between rivals.  The positions handed in are candidate energies of a
    single level, so they predict the same transitions to the same partners
    with the same gA: whatever the calculated branch strengths get wrong about
    how that level divides its strength among its delta J classes is common to
    all of them, and only the difference between their fingerprints says
    anything about where the level is.  Subtracting each ln_J outright would
    charge the intensity model's error to the level's whereabouts - over this
    run the largest ln_J belongs to 059003.000069, whose position is worth
    ln R = +146 on 82 observable predictions and is not in any doubt.

    So the cleanest of the rivals sets the zero: with b = min ln_J,

        ln_R_J = ln_R - (ln_J - b)

    The best-fitting fingerprint keeps its ln R exactly, each other position
    pays what its own fingerprint is worse by, and a level whose rivals are
    equally clean - or which has no rival at all - pays nothing.

    Returns a list of the offset values, in the order given.
    """
    if not len(ln_j):
        return [float(x) for x in ln_r]
    b = min(float(x) for x in ln_j)
    return [float(r) - (float(j) - b) for r, j in zip(ln_r, ln_j)]


def level_kinds(per):
    """The levels of the run grouped by (J, parity), energies sorted.

    Returns ``(kinds, kind_of)``: ``{(J, parity): (energies, level ids)}`` and
    ``{level_id: (J, parity)}``.  A level whose J cannot be read is in
    neither.
    """
    j = np.array([j_value(x) for x in per['J']])
    par = np.array([str(x) for x in per['parity']])
    e = per['E_final'].to_numpy(dtype=float)
    lid = per['level_id'].to_numpy()
    kind_of = {str(i): (jj, pp) for i, jj, pp in zip(lid, j, par)
               if np.isfinite(jj)}
    kinds = {}
    for key in set(kind_of.values()):
        m = (j == key[0]) & (par == key[1])
        o = np.argsort(e[m])
        kinds[key] = (e[m][o], lid[m][o])
    return kinds, kind_of


def nearest_same_kind(kinds, kind_of, level_id, e):
    """The level of the same J and parity whose energy is nearest to E.

    Returns ``(level_id, E - E(that level))``, or ``('', nan)`` when the run
    holds no OTHER level of that J and parity - and then the alternate cannot
    be an interchange whatever sits near it, so there is nothing to name.

    Restricting the search to one J and one parity is the point of the column.
    An interchange is an exchange of two identities: each level takes the
    other's energy and, with it, the other's lines.  Two levels of different J
    do not have the same lines to take - their transitions go to different
    partners entirely - so a near-coincidence of their energies says nothing,
    and a swap of them is not a thing that can be carried out.  The level
    itself is excluded as well: an alternate is always at least ALT_SEP from
    where the level stands, so the level's own position is never the answer,
    and naming it would turn a plain move into a self-interchange.
    """
    key = kind_of.get(str(level_id))
    if key is None:
        return '', float('nan')
    ev, iv = kinds.get(key, (np.empty(0), np.empty(0, dtype=object)))
    m = np.array([str(x) != str(level_id) for x in iv], dtype=bool)
    if not m.any():
        return '', float('nan')
    ev, iv = ev[m], iv[m]
    k = int(np.argmin(np.abs(ev - float(e))))
    return str(iv[k]), float(e - ev[k])


def observed_index(ctx, wn):
    """Where these recorded wavenumbers sit in the observed line list.

    The detail table names a matched line by its wavenumber; the accepted
    assignments name it by its position in ctx.wn_o.  This turns the first
    into the second, so that the two can be compared.
    """
    wn = np.asarray(wn, dtype=float)
    i = np.clip(np.searchsorted(ctx.wn_o, wn), 0, len(ctx.wn_o) - 1)
    j = np.clip(i - 1, 0, len(ctx.wn_o) - 1)
    return np.where(np.abs(ctx.wn_o[j] - wn) < np.abs(ctx.wn_o[i] - wn), j, i)


def own_light(ctx, own):
    """How much observed light each of the level's own lines brings it.

    ``own`` is ctx.own_claim's entry for the level: {index of the observed
    line -> the predicted intensity of this level's component on it}.  What
    the level can claim of a feature it shares is the feature's measured
    intensity times the branching fraction of its own component, never the
    whole feature - the same rule the rest of the pipeline applies to a
    blend.  Returns {index -> the light that is the level's}.
    """
    out = {}
    for i, i_pred in own.items():
        i = int(i)
        tot = float(ctx.claimed_tot[i]) if i < len(ctx.claimed_tot) else 0.0
        share = (float(i_pred) / tot) if tot > 0 else 1.0
        out[i] = float(ctx.int_o[i]) * min(max(share, 0.0), 1.0)
    return out


def support(ctx, level_id, e):
    """What the position E rests on: its free lines, their weight, its spread.

    A matched prediction is a REAL match (n_real) when it could have been
    recorded at all (P_obs >= P_SEEN), the line sits within FREE_SIGMA of
    where it is predicted, the row is not one the `bright` reading won (a
    feature far brighter than the prediction, which classify_lines leaves free
    for a better transition), and the row's own ln R_t is positive - the
    recorded line, position and brightness together, is more likely with the
    level there than as a coincidence.  On a free feature that last condition
    is the intensity test: a line nine times fainter than predicted fails it.
    On a feature other transitions already claim it is the blend test: the
    component's light must improve, or at least not spoil, the account of the
    feature's brightness.  It is FREE support when, in addition, the observed
    feature is claimed by no other accepted transition (C = 0).  Without
    the first condition a prediction with a 1.5 per cent chance of being
    recorded counts as support the moment any line happens to lie within the
    matching window - a coincidence, and the same one that makes n_match the
    number not to quote.  Those are the lines the level can take without
    robbing another level of its evidence; support that is entirely blended
    is not support, because a prediction can be dumped on an
    already-explained feature almost anywhere.

    top_share is the largest single row's share of all the positive evidence.
    A position whose case is one line is not a case: one line can be a
    coincidence, and the scan looked at tens of thousands of positions.

    n_own and n_kept say what the move would cost.  n_own is the number of
    recorded lines the level is assigned NOW; n_kept is how many of those the
    position E still matches.  The two together are what separates a level
    that stays put from a level that moves, and the distance between the two
    energies does not do it: 059003.000617 keeps NONE of its three lines at a
    position 1.3 cm^-1 away, so every one of them would have to be given up
    and another set assigned.

    kept_light is the same cost weighed rather than counted: the share of the
    level's own observed light - each feature's intensity times the branching
    fraction of the level's component in it - that the position still matches.
    Counting alone cannot see which line is being given up, and the one being
    given up is usually the one that matters.  059003.000457 keeps six of its
    seven lines at a position 0.5 cm^-1 away, which reads as one bad line
    dragging the fit; but the line it drops is its strongest by an order of
    magnitude, and a position bought by discarding a level's brightest branch
    is not a refit however many faint lines it keeps.
    """
    v, tab = ln_ratio(ctx, level_id, np.array([float(e)]), detail=True)
    light = own_light(ctx, ctx.own_claim.get(level_id, {}))
    own = set(light)
    if tab is None:
        return dict(ln_R=0.0, n_free=0, free_gain=0.0,
                    top_share=np.nan, n_obs_alt=0, n_seen_alt=0, n_real=0,
                    n_miss_alt=0, n_own=len(own), n_kept=0,
                    kept_light=0.0 if own else np.nan, ln_J=0.0, j_classes=[])
    m = tab['matched'].to_numpy(dtype=bool)
    r = tab['ln_R'].to_numpy(dtype=float)
    seen, real, free = match_kinds(tab)
    j = ln_j_pattern(tab, ctx.j_of.get(str(level_id), float('nan')),
                     ctx.j_of, seen, real)
    pos = r[m & (r > 0.0)]
    here = set(int(x) for x in
               observed_index(ctx, tab.loc[m, 'wn_obs'].to_numpy()))
    whole = sum(light.values())
    return dict(ln_R=float(v[0]), n_free=int(free.sum()),
                n_own=len(own), n_kept=len(own & here),
                kept_light=(sum(light[i] for i in own & here) / whole
                            if whole > 0 else (np.nan if not own else 0.0)),
                free_gain=float(r[free].sum()),
                n_obs_alt=int(seen.sum()),
                n_seen_alt=int((seen & m).sum()),
                n_real=int(real.sum()),
                n_miss_alt=int((seen & ~m).sum()),
                ln_J=j['ln_J'], j_classes=j['classes'],
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

      interchange  the alternate lands on another level of the run OF THE
                   SAME J AND PARITY - the only kind of level whose identity
                   this one could be exchanged with.  Nothing moves anywhere
                   new; the two identities may be swapped, and
                   level_interchange.py decides that on evidence this scan
                   does not look at.
      refit        the alternate is a fraction of a wavenumber away AND
                   the level keeps more than half of the recorded lines it
                   is assigned now.  Nothing is re-identified: the level
                   stays, and some accepted line is dragging the LOPT fit off
                   the position its own lines want.  Both halves are needed.
                   A short move that keeps none of the lines is not a refit
                   at all - every one of its assignments has to be dropped,
                   another set made, and the level entered in a ledger at a
                   new position - and calling it one told the reader to do
                   the opposite of the work the level actually needs.
      top line     the same short move, keeping most of the level's lines by
                   count, but giving up the greater part of its light: the
                   line it drops is the level's strongest.  The alternate is
                   then bought with the one assignment the level can least
                   afford to lose, and the Ritz mismatch it removes is the
                   mismatch of the brightest branch, which is a thing to
                   explain and not a thing to discard.  Look at that line -
                   is it a blend, is its wavenumber right, is the
                   identification right - rather than at the position.
      relocate     a free position, broadly supported by lines nobody is
                   using, surviving the look-elsewhere correction, and
                   convincing in its own right - ln R must be positive there,
                   after the delta J fingerprint has been charged against it
                   (ln_R_alt_J), and not merely better than where the level is
                   now.  How far
                   away it is does not enter: a position 1 cm^-1 away that
                   takes a different set of lines is as much a relocation as
                   one 300 cm^-1 away, and needs the same ledger record.
      weak         a free position whose support is thin.  Leave it until the
                   region around it is settled - the neighbours will change
                   the answer.
      no support   one line, or none that is free.  No case at all.
    """
    if np.isfinite(row['near_dE']) and abs(
            row['near_dE']) < INTERCHANGE_DE:
        return 'interchange'
    if (np.isfinite(row['dE_alt']) and abs(row['dE_alt']) < REFIT_DE
            and row['n_kept'] > REFIT_KEEP * row['n_own']):
        kept = row.get('kept_light', np.nan)
        if np.isfinite(kept) and kept <= REFIT_LIGHT:
            return 'top line'
        return 'refit'
    share = row['top_share']
    if row['n_free'] < AUDIT_MIN_FEW or (np.isfinite(share)
                                         and share > AUDIT_ANY_SHARE):
        return 'no support'
    if (row['n_free'] >= AUDIT_MIN_FREE
            and row.get('ln_R_alt_J', row['ln_R_alt']) > 0.0
            and row['free_gain'] >= AUDIT_MIN_GAIN
            and np.isfinite(share) and share <= AUDIT_MAX_SHARE
            and row['look'] >= AUDIT_MIN_LOOK):
        return 'relocate'
    return 'weak'


def audit_table(ctx, tab, e_calc, w_of, vac):
    """Add the audit columns to a scanned report table.

    gain_R   ln R the alternate wins by on the line evidence alone
             (positive: the lines prefer it)
    gain     gain_R + d_ln_J: the same, once the delta J fingerprint is
             counted with it.  This is the column the dispositions are taken
             on, and the one to read
    look     that gain after the look-elsewhere correction.  A level allowed
             a 7000 cm^-1 window gets hundreds of chances at a good maximum
             and a level allowed 130 gets a handful, so the two cannot be
             compared raw.  Under H0 the heights of the local maxima are
             roughly exponential, so the best of n_alt of them stands about
             ln(n_alt) above a typical one; requiring the gain to beat that
             puts every level on the same footing.
    n_free   lines supporting the alternate that no other level is using
    n_own    recorded lines the level is assigned now
    n_kept   how many of those the alternate still matches - the test that
             separates a refit from a relocation
    kept_light  the share of the level's own observed light those kept lines
             carry: the test that separates a refit from a move that survives
             by dropping the level's brightest line
    free_gain  what those lines are worth
    top_share  the largest row's share of the positive evidence
    ln_J_alt   the delta J fingerprint of the alternate (ln_j_pattern): how
             much better the pattern of absences THERE fits a level of another
             J than this one.  The column beside it, ln_J, is the same test at
             the adopted position.  A move that raises it is a move onto lines
             that do not belong to this level's J
    d_ln_J   ln_J - ln_J_alt, and THE COLUMN TO READ.  ln_J at a single
             position mixes two things: a position that belongs to a level of
             another J, and a calculated gA that misjudges how this level's
             strength is divided among its dJ branches.  The second is a
             property of the level's wavefunction, not of its energy, and it
             is the larger of the two over this run - 059003.000069 scores
             ln_J = 6.0 at a position worth ln R = +146 on 82 observable
             predictions, which is a branch-strength fault and not a doubt
             about where the level is.  The difference between two candidate
             energies FOR THE SAME LEVEL cancels it: the same partners, the
             same gA, only the wavenumbers differ.  Positive d_ln_J means the
             alternate fits this level's J better than the adopted position
             does, and gain counts it
    ln_R_J   ln_R and ln_R_alt with the fingerprint folded in: each is offset
    ln_R_alt_J  by the cleaner of the two, b = min(ln_J, ln_J_alt), so that
             the better-fitting position keeps its ln R untouched, the other
             pays what its fingerprint is worse by, and the difference between
             them is gain.  A level whose two positions have the same
             fingerprint pays nothing, which is right: there is then nothing
             in the pattern of absences to choose between them
    n_vacant   unfound calculated levels of the same J and parity within one
             configuration window of the alternate - can an UNKNOWN level be
             there instead?
    z_alt    (E_alt - E_calc)/W: how far the alternate sits from where the
             calculation puts the level, in units of that configuration's own
             scatter
    action   see disposition()
    """
    cols = dict(n_free=[], free_gain=[], top_share=[], gain=[], gain_R=[],
                look=[], n_obs_alt=[], n_seen_alt=[], n_miss_alt=[], n_own=[],
                n_kept=[], kept_light=[], n_vacant=[], z_alt=[], ln_J_alt=[],
                d_ln_J=[], ln_R_J=[], ln_R_alt_J=[], action=[])
    lev = ctx.per.set_index('level_id')
    for _, r in tab.iterrows():
        lid = r['level_id']
        if not r.get('n_alt', 0) or not np.isfinite(r.get('dE_alt', np.nan)):
            for k in cols:
                cols[k].append('' if k == 'action' else np.nan)
            continue
        e_a = float(r['E']) + float(r['dE_alt'])
        s = support(ctx, lid, e_a)
        # THE FINGERPRINT ENTERS HERE, AND ONLY AS A DIFFERENCE.  The same
        # level at two energies predicts the same transitions to the same
        # partners with the same gA, so whatever its calculated branch
        # strengths get wrong is the same at both and cancels; what is left
        # is the pattern of absences that belongs to one J rather than the
        # other.  Offsetting both positions by the cleaner of the two leaves
        # the better fingerprint's ln R alone and charges the worse one the
        # difference, so nothing is taken from a level whose two positions
        # are equally clean.
        d_j = float(r['ln_J']) - s['ln_J']
        r_j, a_j = fold_j([r['ln_R'], r['ln_R_alt']],
                          [r['ln_J'], s['ln_J']])
        gain_R = float(r['ln_R_alt']) - float(r['ln_R'])
        gain = a_j - r_j
        look = gain - math.log(max(int(r['n_alt']), 1))
        w = float(w_of.get(lid, np.nan))
        ec = float(e_calc.get(lid, np.nan))
        row = dict(r)
        row.update(s)
        row['gain'] = gain
        row['gain_R'] = gain_R
        row['look'] = look
        row['ln_R_alt_J'] = a_j
        cols['n_free'].append(s['n_free'])
        cols['free_gain'].append(s['free_gain'])
        cols['top_share'].append(s['top_share'])
        cols['n_obs_alt'].append(s['n_obs_alt'])
        cols['n_seen_alt'].append(s['n_seen_alt'])
        cols['n_miss_alt'].append(s['n_miss_alt'])
        cols['n_own'].append(s['n_own'])
        cols['n_kept'].append(s['n_kept'])
        cols['kept_light'].append(s['kept_light'])
        cols['ln_J_alt'].append(s['ln_J'])
        cols['d_ln_J'].append(d_j)
        cols['ln_R_J'].append(r_j)
        cols['ln_R_alt_J'].append(a_j)
        cols['gain'].append(gain)
        cols['gain_R'].append(gain_R)
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


ACTION_ORDER = ['relocate', 'interchange', 'refit', 'top line', 'weak',
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
            'ln_R_alt', 'n_seen_alt', 'n_miss_alt', 'gain_R', 'd_ln_J',
            'gain', 'n_alt', 'look',
            'n_free', 'free_gain', 'top_share', 'n_own', 'n_kept',
            'kept_light', 'n_vacant', 'z_alt', 'ln_J', 'ln_J_alt']
    with pd.option_context('display.width', 220, 'display.max_columns', 24):
        print(f"\n  firm grounds for relocation ({len(firm)}):")
        print(firm[cols].to_string(index=False, na_rep='-')
              if len(firm) else '    none')
        swap = pref[pref['action'] == 'interchange']
        if len(swap):
            print(f"\n  interchanges, for level_interchange.py ({len(swap)}):")
            print(swap[cols + ['near_level', 'near_dE']].to_string(
                index=False, na_rep='-'))
        rf = pref[pref['action'] == 'refit']
        if len(rf):
            print(f"\n  the level stays where it is and keeps its lines; an "
                  f"accepted line is dragging the fit ({len(rf)}):")
            print(rf[cols].to_string(index=False, na_rep='-'))
        tl = pref[pref['action'] == 'top line']
        if len(tl):
            print(f"\n  the alternate is bought by giving up the level's "
                  f"strongest line ({len(tl)}).  It keeps most of the "
                  f"assignments by count\n  and less than half their light: "
                  f"read kept_light against n_kept/n_own.  What wants "
                  f"looking at is that line -\n  its wavenumber, whether it "
                  f"is a blend, whether the identification is right - and "
                  f"not the position:")
            print(tl[cols].to_string(index=False, na_rep='-'))
        wk = pref[pref['action'] == 'weak'].sort_values(
            'look', ascending=False)
        if len(wk):
            print(f"\n  the level would move onto a different set of lines - "
                  f"a relocation, and a ledger record - but the case for the "
                  f"new position is thin ({len(wk)}).\n  Read n_kept against "
                  f"n_own: that is how many of its present assignments it "
                  f"would have to give up:")
            print(wk[cols].to_string(index=False, na_rep='-'))
    hole = pref[(pref['action'] == 'relocate') & (pref['n_vacant'] > 0)]
    if len(hole):
        print(f"\n  {len(hole)} of the relocations sit where the calculation "
              f"still has an unfound level of the same J and parity, so the "
              f"lines may belong to THAT level rather than to this one: "
              f"{', '.join(hole['level_id'])}")
    print('\n  A relocation changes the lines available to its neighbours, so '
          'accept them\n  one at a time and re-scan: the verdict on every '
          'level is conditional on the\n  rest of the list.')


# ---------------------------------------------------------------------------
# Levels nobody has found: searching by IDEN2 row
# ---------------------------------------------------------------------------
# Everything above scans a level the run already has: an adopted energy, a set
# of accepted lines, a level_id.  A level of the CALCULATION that has never
# been found has none of those.  What it has is a row in IDEN2/enlev.dat - a
# calculated energy E_calc, a J, a configuration label - and a place in
# Cowan's transition list, which gives it transitions to every other
# calculated level.  The ones whose other end HAS been found are predictions
# with a known partner energy: at any trial energy E they predict a
# wavenumber, and a wavenumber is all ln R needs.  So the same likelihood
# ratio can be evaluated for a level that does not exist yet, and scanning it
# across the window the calculation allows IS the search.  Where ln R > 0 the
# recorded lines are more likely with a level at that energy than with nothing
# there; where it is largest they are most likely.
#
# This is the scan the user does by hand in IDEN2, with the level list scrolled
# to a trial position and the predicted transitions checked against the plate
# list one by one.  The arithmetic is the same arithmetic; what it adds is
# that every position in the window is tried instead of the ones a person has
# the patience for, and that the absences count against a position as well as
# the coincidences count for it.
#
# THREE THINGS DIFFER from the scan of a level the run already has, and all
# three make the answer harsher rather than kinder:
#
#   - The level owns no accepted line.  Every feature that an existing
#     identification already claims therefore counts against it in full
#     through C (section 5): to be believed, an unfound level must explain
#     lines that the levels already found have left alone.  That is the right
#     test, and it is stricter than the one a known level faces.
#   - There is no adopted position to compare against, so there is no gain.
#     ln R itself is the verdict.
#   - The look-elsewhere problem is the whole window rather than the distance
#     between two candidates: a scan that returns n positive maxima has had n
#     chances at every one of them, so look = ln R - ln n is what a position
#     is worth after the search that found it is paid for.  It is the same
#     correction the audit applies to a relocation.
#
# THE PREDICTED INTENSITIES ARE FIXED AT E_calc and not recomputed as the scan
# moves.  I = C gA (nu/1e8) exp(-E_up/kT) does depend on the trial energy
# through both nu and the Boltzmann factor, but over a window of a few hundred
# wavenumbers out of a hundred thousand that is a few per cent, against an
# observed-to-predicted scatter of a factor of three.  It is also exactly what
# the scan of a known level already does, so the two are comparable.
UNKNOWN_ID = '(unknown)'   # the level_id an unfound level is scanned under


def unfound_theory(ctx, log=print):
    """The calculated level list, Cowan's transitions and the partner energies.

    Returned as ``(en, trans, mapping, e_meas, per_cfg, whole, id_of)``: the
    rows of enlev.dat indexed by their IDEN2 row number, Cowan's E1 transition
    list, the ``{Cowan level number: IDEN2 row}`` correspondence, the measured
    energy of every found row, the rms of E_obs - E_calc per configuration and
    over the whole list, and the ``{IDEN2 row: level_id}`` map.

    unfound_levels.py imports THIS module, so it can only be imported from
    inside a function here, once this module is built.  It is worth the
    awkwardness: the reading of Cowan's table, the matching of its level
    numbering to IDEN2's and the per-configuration windows are all measured
    and tested there, and a second copy of them would be a second answer.
    """
    import cowan_gA
    import unfound_levels as uf
    en, trans, mapping = uf.read_theory(log=log)
    e_meas = uf.measured_energies(en, e_final=ctx.e_final, log=log)
    per_cfg, whole = uf.windows(en, log=log)
    return en, trans, mapping, e_meas, per_cfg, whole, cowan_gA.read_id_map(
        uf.IDS)


def register_unknown(ctx, idx, en, trans, mapping, e_meas, id_of, log=print):
    """Give the calculated level of IDEN2 row ``idx`` a set of predictions.

    They are entered in ctx under UNKNOWN_ID, so that ln_ratio, support and
    print_detail work on it unchanged.  Only transitions to levels that have
    been found are kept: the other end of the rest is itself unplaced, so no
    wavenumber can be predicted for them at all.  Returns the number of
    predictions the level is left with.

    WHEN THE ROW IS A LEVEL THE RUN HAS ALREADY FOUND - and it may be: the
    search is run on found rows to ask whether the position they hold is the
    one the lines want - that level's own accepted assignments are handed to
    the search.  ctx.own_claim[UNKNOWN_ID] is set to the level's own claims,
    which is exactly what the scan of a known level does with them: through C
    in ln_ratio, a feature the level itself claims is then read as FREE, and
    only the light OTHER levels have put on it counts against the position.

    Without this the search is asked an impossible question.  Every line the
    level is assigned now is a line some level has claimed, so at the level's
    own position each of its own assignments is charged to it as a blend with
    itself: 059003.000538, which holds 139081.51 cm^-1 on four accepted
    lines and scores ln R = +12.9 as a known level, came out at +0.4 with no
    free line at all, ranked below three positions in the window that rest on
    nothing.  The level was being made to compete with itself, and it lost.

    n_own and n_kept in the support columns then read as they do in the
    audit: how many recorded lines the level is assigned now, and how many of
    them a candidate position still matches.
    """
    a = np.array([mapping.get(int(v), -1) for v in trans['lid1']])
    b = np.array([mapping.get(int(v), -1) for v in trans['lid2']])
    gA = trans['gA'].to_numpy(dtype=float)
    up = (trans['u_gA_pct'].to_numpy(dtype=float)
          if 'u_gA_pct' in trans.columns else np.full(len(trans), np.nan))
    hit = ((a == idx) | (b == idx)) & (a > 0) & (b > 0) & (gA > 0)
    other = np.where(a[hit] == idx, b[hit], a[hit])
    gA, up = gA[hit], up[hit]
    known = en['known'].to_dict()
    keep = np.array([bool(known.get(int(i), False)) for i in other],
                    dtype=bool)
    other, gA, up = other[keep], gA[keep], up[keep]

    pid = np.array([str(id_of.get(int(i), '')) for i in other])
    in_run = np.array([p in ctx.e_final for p in pid], dtype=bool)
    if not in_run.all():
        log(f"  {int((~in_run).sum())} partners are found in enlev.dat but "
            f"are not levels of this run; their transitions are dropped")
    other, gA, pid, up = other[in_run], gA[in_run], pid[in_run], up[in_run]

    E_c = float(en['E_calc'][idx])
    E_p = np.array([float(e_meas[int(i)]) for i in other])
    C = float(cl.CFG.intensity_model['C'])
    kT = float(cl.CFG.intensity_model['kT'])
    nu = np.abs(E_c - E_p)
    i_pred = C * gA * (nu / 1.0e8) * np.exp(-np.maximum(E_c, E_p) / kT)

    ctx.by_level[UNKNOWN_ID] = dict(
        partner=pid,
        sign=np.where(E_c > E_p, 1.0, -1.0),
        i_pred=i_pred,
        e_m=E_p,
        u_m=np.array([ctx.u_M.get(p, 0.0) for p in pid]),
        u_calc=np.array([ls._u_ln(v) for v in up]),
        degenerate=np.array([ctx.n_acc_level.get(p, 0) <= 0 for p in pid]))
    ctx.e_final[UNKNOWN_ID] = E_c
    ctx.j_of[UNKNOWN_ID] = j_value(en['J'][idx])

    lid = str(id_of.get(int(idx), ''))
    ctx.unknown_level = None
    ctx.own_claim.pop(UNKNOWN_ID, None)
    if lid:
        own = ctx.own_claim.get(lid, {}) if lid in ctx.e_final else {}
        ctx.unknown_level = dict(
            level_id=lid, in_run=lid in ctx.e_final,
            E=float(ctx.e_final[lid]) if lid in ctx.e_final else float('nan'),
            n_own=len(own), released=lid in getattr(ctx, 'released', set()))
        if own:
            ctx.own_claim[UNKNOWN_ID] = dict(own)
    return len(pid)


def questionable_levels(path, log=print):
    """The levels the last report of ``path`` marks with a question mark.

    That is the `question` column of level_positions.csv: '?' wherever the
    scan found an alternate position within --alt-drop of the adopted one.  It
    is the run's own record of which positions are not settled, written by the
    last --scan or --audit.
    """
    if not os.path.exists(path):
        raise SystemExit(f"{path} not found: --drop-all-questionable needs "
                         f"the report of a previous --scan or --audit run")
    t = pd.read_csv(path, dtype={'level_id': str})
    if 'question' not in t.columns:
        raise SystemExit(f"{path} has no `question` column: it was not "
                         f"written by a --scan or --audit run")
    q = t[t['question'].astype(str).str.strip() == '?']
    return [str(x) for x in q['level_id']]


def release_levels(ctx, ids, args, log=print):
    """Drop every accepted assignment of these levels; free their lines.

    What a released level loses is its CLAIM on the recorded features: its
    accepted rows are subtracted from C, from the count of accepted
    transitions on each feature, and from the own_claim of both levels the row
    joined.  A feature that no other accepted transition claims is then free,
    and the unrelated-line background is re-measured over the enlarged free
    set - a release makes free lines commoner, and a commoner background is a
    smaller ln R, so the position being searched for is not handed the lines
    for nothing.

    What a released level KEEPS is its adopted energy.  It goes on serving as
    a partner: a transition to it is still predicted at a definite wavenumber,
    and that is the whole point of releasing it - the lines it was holding are
    offered to the level under test while the predictions that reach them
    survive.  Its energy is of course conditional on the very assignments that
    have just been dropped, which is the standing caveat on a questionable
    level and the reason it is marked one.

    Returns the number of accepted rows dropped.
    """
    ids = set(str(i) for i in ids)
    ctx.released = ids
    acc = ctx.acc
    lo = acc['low_id'].astype(str).to_numpy()
    up = acc['upp_id'].astype(str).to_numpy()
    sel = np.array([(a in ids) or (b in ids) for a, b in zip(lo, up)],
                   dtype=bool)
    if not sel.any():
        log(f"  none of the {len(ids)} levels named carries an accepted "
            f"line; nothing is released")
        return 0
    li_ = acc['line_idx'].to_numpy(dtype=int)[sel]
    ic = acc['calc_intens'].fillna(0.0).to_numpy(dtype=float)[sel]
    np.add.at(ctx.claimed_tot, li_, -ic)
    ctx.claimed_tot = np.maximum(ctx.claimed_tot, 0.0)
    np.add.at(ctx.n_acc_line, li_, -1)
    ctx.n_acc_line = np.maximum(ctx.n_acc_line, 0)
    for a, b, i, c in zip(lo[sel], up[sel], li_, ic):
        for who in (a, b):
            m = ctx.own_claim.get(who)
            if m is None or int(i) not in m:
                continue
            m[int(i)] -= float(c)
            if m[int(i)] <= 1e-12:
                del m[int(i)]
            if not m:
                ctx.own_claim.pop(who, None)
    # ctx.n_acc_level is NOT recomputed: it is what marks a partner degenerate
    # - a level with no accepted line of its own has no independently measured
    # energy - and a released level's energy is still the one the run adopted.
    # Recomputing it here would strike out exactly the predictions that reach
    # the lines the release has just freed, which is the opposite of the
    # intent.
    n_free_now = int((ctx.n_acc_line == 0).sum())
    free_line_background(ctx, args, log=lambda *a, **k: None)
    log(f"  released {int(sel.sum())} accepted assignments of "
        f"{len(ids)} questionable levels; {n_free_now} of "
        f"{len(ctx.wn_o)} recorded lines now carry no accepted transition")
    log(f"    the unrelated-line background is re-measured over them: "
        f"density {np.median(ctx.rho_bg_o):.4f} per cm^-1 (median), ln I of "
        f"a free line {np.nanmedian(ctx.bg_mu_f):.2f}")
    return int(sel.sum())


def scan_unknown(ctx, idx, en, per_cfg, whole, step=GRID_STEP, min_ln_r=0.0,
                 k=WINDOW_K):
    """Every position in the calculation's window where the lines want a level.

    The window is E_calc +/- k (three) times the rms of E_obs - E_calc over
    the FOUND levels of the same configuration - the same interval a known
    level of that configuration is scanned over, and for the same reason: how
    far the calculation can be wrong is a property of the configuration, not
    of the level.  That reason fails for a level strongly mixed with a level
    of another configuration, whose label says little about which
    configuration's scatter it shares; --window-sigmas widens the window for
    it.  Maxima closer together than ALT_SEP are one maximum.
    """
    E_c = float(en['E_calc'][idx])
    W = float(per_cfg.get(en['cfg'][idx], whole))
    lo, hi = scan_interval(E_c, E_c, W, k=k)
    n = int((hi - lo) / step) + 1
    if n > GRID_MAX:
        step = (hi - lo) / (GRID_MAX - 1)
        n = GRID_MAX
    grid = lo + step * np.arange(n)
    y, _ = ln_ratio(ctx, UNKNOWN_ID, grid)
    hits = [(float(grid[i]), float(y[i])) for i in local_maxima(y)
            if y[i] > min_ln_r]
    hits.sort(key=lambda z: -z[1])
    kept = []
    for e, v in hits:
        if all(abs(e - k[0]) >= ALT_SEP for k in kept):
            kept.append((e, v))
    return dict(idx=idx, E_calc=E_c, W=W, lo=lo, hi=hi, step=step,
                min_ln_r=min_ln_r, grid=grid, y=y, positions=kept,
                ln_R_best=float(y.max()) if len(y) else float('nan'),
                E_best=float(grid[int(np.argmax(y))]) if len(y) else
                float('nan'))


def unknown_table(ctx, r):
    """One row per candidate position, best first.

    ln_R    what the recorded lines alone say about that energy
    ln_R_J  the same with the delta J fingerprint folded in, and THE COLUMN
            THE TABLE IS RANKED AND JUDGED ON.  The candidates are positions
            of ONE level, so they predict the same transitions to the same
            partners with the same gA and whatever the calculated branch
            strengths get wrong is common to them all; only the difference
            between their fingerprints says anything about where the level
            is.  Every candidate is therefore offset by the cleanest of them,
            ln_R_J = ln_R - (ln_J - min ln_J): the best-fitting candidate
            keeps its ln R exactly, each other pays what its fingerprint is
            worse by, and a scan whose candidates are all equally clean is
            left as it was
    look    ln_R_J after the look-elsewhere correction for the whole scan
    dE, z   how far it is from the calculated energy, in cm^-1 and in units
            of that configuration's own scatter.  A position at z = 2.8 is
            asking the calculation to be wrong by nearly three times as much
            as it usually is for that configuration.
    n_obs   predictions that could have been recorded there (P_obs >= P_SEEN)
    n_match how many of those REALLY match a recorded line (n_real of
            support(): within FREE_SIGMA, and a line the row's own ln R_t
            counts as evidence for the level - brightness agreeing with the
            prediction on a free feature, admissible as a component on a
            claimed one) - the count to read against n_obs, since a position
            that predicts twenty observable lines and matches three is
            contradicted, not supported
    n_poor  observable predictions with a recorded line in the matching
            window that is NOT a match: off-centre, or of the wrong
            brightness.  Each counts against the position about as much as
            an absence does
    n_miss  observable predictions with no recorded line in the window at
            all; n_obs = n_match + n_poor + n_miss
    n_free  of the matches, the well-centred ones on features no accepted
            transition already claims: the lines the level could take without
            taking them from a level already found
    free_gain  what those free lines are worth in ln R
    top_share  the largest single line's share of the positive evidence
    ln_J    the delta J fingerprint: how much better the pattern of matches
            and absences fits a level of a DIFFERENT J than this one.  Zero
            where the lines are spread over the delta J classes as their
            predicted intensities say they should be; a few nats where a whole
            class that should have been recorded is missing.  ln_R - ln_R_J
            is what it cost this candidate against the cleanest one
    verdict how the audit above would read that support, by exactly the same
            constants: `firm` where a relocation would be called firm, `no
            support` where one line is the whole case or there are fewer than
            two free ones, `weak` in between.  A level nobody has found is
            held to the standard a level of the run is held to, because the
            evidence for it is the same kind of evidence and the scan that
            found it had the same freedom to look.
    """
    rows = []
    n = max(len(r['positions']), 1)
    sups = [support(ctx, UNKNOWN_ID, e) for e, _ in r['positions']]
    # the cleanest fingerprint among the candidates carries no charge; fold_j
    # says why only the difference between them means anything
    folded = fold_j([v for _, v in r['positions']],
                    [sp['ln_J'] for sp in sups])
    for (e, v), sup, v_j in zip(r['positions'], sups, folded):
        rows.append(dict(E=e, ln_R=v, ln_R_J=v_j, look=v_j - math.log(n),
                         dE=e - r['E_calc'],
                         z=(e - r['E_calc']) / r['W'] if r['W'] > 0
                         else np.nan,
                         n_obs=sup['n_obs_alt'], n_match=sup['n_real'],
                         n_poor=sup['n_seen_alt'] - sup['n_real'],
                         n_miss=sup['n_miss_alt'], n_free=sup['n_free'],
                         free_gain=sup['free_gain'],
                         top_share=sup['top_share'], ln_J=sup['ln_J'],
                         verdict=position_verdict(
                             sup['n_free'], sup['free_gain'],
                             sup['top_share'], v_j - math.log(n))))
    tab = pd.DataFrame(rows, columns=['E', 'ln_R', 'ln_R_J', 'look', 'dE',
                                      'z',
                                      'n_obs', 'n_match', 'n_poor',
                                      'n_miss',
                                      'n_free', 'free_gain', 'top_share',
                                      'ln_J', 'verdict'])
    if tab.empty:
        return tab
    rank = {'firm': 0, 'weak': 1, 'no support': 2}
    tab['_r'] = [rank[x] for x in tab['verdict']]
    return tab.sort_values(['_r', 'ln_R_J'], ascending=[True, False]).drop(
        columns='_r').reset_index(drop=True)


def position_verdict(n_free, free_gain, top_share, look):
    """What the audit's own constants make of one candidate position.

    The tests are those of disposition(): AUDIT_MIN_FEW free lines before
    there is any support at all, AUDIT_ANY_SHARE above which one line IS the
    case, and the four a firm relocation has to pass.  There is no ln R of an
    adopted position to beat here - the level has no adopted position - so ln
    R itself, after the look-elsewhere correction, plays the part the gain
    plays there.  What is handed in as `look` is the look-elsewhere correction
    of ln_R_J, not of ln_R, so a candidate whose delta J pattern is worse than
    another candidate's has to make that up in line evidence before it can be
    called firm.
    """
    share = top_share
    if n_free < AUDIT_MIN_FEW or (np.isfinite(share)
                                  and share > AUDIT_ANY_SHARE):
        return 'no support'
    if (n_free >= AUDIT_MIN_FREE and free_gain >= AUDIT_MIN_GAIN
            and np.isfinite(share) and share <= AUDIT_MAX_SHARE
            and look >= AUDIT_MIN_LOOK):
        return 'firm'
    return 'weak'


def print_unknown(ctx, idx, en, r, n_pred, top=0):
    """The report for one unfound level."""
    row = en.loc[idx]
    print(f"\nIDEN2 row {idx}  {row['label']}  J = {row['J']}  "
          f"{'FOUND' if row['known'] else 'not found'}")
    held = getattr(ctx, 'unknown_level', None)
    if held is not None and not held['in_run']:
        print(f"  enlev.dat calls this row found, as {held['level_id']}, but "
              f"that level is not one of this run's: it has nothing of its "
              f"own here, and the window is searched as for a level nobody "
              f"has found")
    elif held is not None:
        n_own = held['n_own']
        print(f"  this row IS a level of the run: {held['level_id']}, which "
              f"holds {held['E']:.3f} cm^-1")
        if n_own:
            print(f"  its {n_own} accepted line"
                  f"{' is' if n_own == 1 else 's are'} released for the "
                  f"search, so the position it holds now is scanned on the "
                  f"same footing as every other: it does not have to compete "
                  f"with itself for its own lines")
        elif held['released']:
            print(f"  it is itself one of the questionable levels, so its "
                  f"assignments have already been freed for every level in "
                  f"this run of the search, this one included")
        else:
            print(f"  it carries no accepted line, so there is nothing of "
                  f"its own to release")
    print(f"  calculated at {r['E_calc']:.1f} cm^-1; levels of {row['cfg']} "
          f"turn out to be {r['W']:.1f} cm^-1 from where the calculation "
          f"puts them (rms)")
    at_calc = support(ctx, UNKNOWN_ID, r['E_calc'])
    print(f"  {n_pred} calculated transitions to levels that HAVE been found, "
          f"{at_calc['n_obs_alt']} of which could have been recorded at the "
          f"calculated position")
    print(f"  scanned {r['lo']:.1f} - {r['hi']:.1f} cm^-1 in steps of "
          f"{r['step']:.3f}")
    if at_calc['n_obs_alt'] < 2:
        print(f"  With fewer than two transitions that could have been "
              f"recorded, nothing in this window can be believed: one line "
              f"can be made to fit any energy, and a scan over hundreds of "
              f"wavenumbers will always find one.  unfound_levels.py ranks "
              f"the levels by how many such transitions they have, and this "
              f"one is not worth a search.")
    cut = r.get('min_ln_r', 0.0)
    if not r['positions']:
        print(f"  no position in the window where ln R > {cut:g}.  The best "
              f"the window offers is {r['ln_R_best']:+.2f} at "
              f"{r['E_best']:.3f} cm^-1"
              + (", which is not a candidate: the recorded lines are no more "
                 "likely with a level there than with none."
                 if r['ln_R_best'] <= 0.0 else "."))
        print(f"  A position just outside the window is not seen: "
              f"--unknown {idx} --at E scores any energy, and "
              f"--window-sigmas widens the search.")
        return
    tab = unknown_table(ctx, r)
    counts = tab['verdict'].value_counts()
    shown = tab if not top else tab.head(top)
    print(f"  {len(tab)} positions where ln R > {cut:g}: "
          + ', '.join(f"{int(counts.get(k, 0))} {k}"
                      for k in ('firm', 'weak', 'no support'))
          + (f"; the best {len(shown)} shown" if len(shown) < len(tab)
             else ''))
    with pd.option_context('display.width', 200, 'display.max_columns', 16):
        print()
        print(shown.to_string(index=False, na_rep='-',
                              float_format=lambda x: f'{x:.3f}'))
    print(f"\n  --unknown {idx} --at E lists the transitions at any one of "
          f"these energies.")
    if not int(counts.get('firm', 0)) and not int(counts.get('weak', 0)):
        print(f"  Every one of them rests on one line or on none that is "
              f"free, which is what a scan of {r['hi'] - r['lo']:.0f} cm^-1 "
              f"finds when there is nothing there.")
    print(f"  A position here is conditional on the rest of the level list, "
          f"exactly as an alternate position is: the lines it takes are "
          f"lines its neighbours could take instead.")


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--detail', metavar='LEVEL_ID',
                   help='per-transition table for one level')
    p.add_argument('--at', metavar='E|alt', default=None,
                   help='energy at which --detail is evaluated: a value in '
                        'cm^-1, or "alt" for the best alternate position the '
                        'scan finds - the one an audit row proposes moving '
                        'to (default: the adopted position)')
    p.add_argument('--detail-all', action='store_true',
                   help='--detail lists every prediction, including the faint '
                        'ones that could not have been recorded and found no '
                        'line (default: only the observable and the matched)')
    p.add_argument('--unknown', nargs='+', metavar='IDEN2_ROW', type=int,
                   default=None,
                   help='search for a level of the CALCULATION by its row '
                        'number in IDEN2/enlev.dat, starting from the E_calc '
                        'of that row: scan the window the calculation allows '
                        'and report every position where ln R > 0, with what '
                        'each one rests on.  The row need not be one of the '
                        'levels that have been found')
    p.add_argument('--drop-all-questionable', '--drop_all_questionable',
                   dest='drop_questionable', action='store_true',
                   help='before searching with --unknown, drop every accepted '
                        'assignment of every level the last report (--out, '
                        'default level_positions.csv) marks questionable - '
                        'the levels whose `question` column is `?`, those the '
                        'scan found an alternate position for - so that their '
                        'lines are free for the level being searched for.  '
                        'The released levels keep their adopted energies and '
                        'go on serving as partners; only their claim on the '
                        'recorded features is dropped, and the unrelated-line '
                        'background is re-measured over the enlarged free set')
    p.add_argument('--min-ln-r', type=float, default=0.0, metavar='X',
                   help='--unknown reports the positions with ln R above this '
                        '(default %(default)s: the energies at which the '
                        'recorded lines are more likely with a level there '
                        'than without one)')
    p.add_argument('--window-sigmas', type=float, default=WINDOW_K,
                   metavar='K',
                   help='--unknown scans E_calc +/- K times the rms of '
                        'E_obs - E_calc of the configuration the level is '
                        'labelled with (default %(default)s).  Widen it for a '
                        'level strongly mixed with a level of another '
                        'configuration: its label, and so the scatter its '
                        'window is taken from, can be the wrong one.  z in '
                        'the table stays in units of the labelled '
                        'configuration scatter')
    p.add_argument('--top', type=int, default=0, metavar='N',
                   help='--unknown lists only the best N positions '
                        '(default: all of them)')
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
    p.add_argument('--firm', default=FIRM_FILE, metavar='PATH',
                   help='registry of the levels whose position is settled - '
                        'ln R at least %g and no alternate anywhere in the '
                        'window.  It is rewritten by every --scan or --audit '
                        'run (default %%(default)s)' % FIRM_LN_R)
    p.add_argument('--no-firm', dest='firm', action='store_const', const=None,
                   help='neither read nor write the registry')
    p.add_argument('--use-firm', action='store_true',
                   help='do not re-scan the levels the registry calls '
                        'settled, when nothing their scan depends on has '
                        'changed since it was written.  This is what makes a '
                        'second run fast; a level is re-scanned the moment '
                        'its energy, one of its partners, or an accepted '
                        'assignment within reach of its predictions moves')
    p.add_argument('--plain-background', action='store_true',
                   help='take the background of unrelated lines from the whole '
                        'observed list, as the runs made before the '
                        'transition-probability file was available did, '
                        'instead of from the free lines and the predicted '
                        'transitions of the levels that have not been found')
    p.add_argument('--lopt', default=None,
                   help='build the run from a LOPT line-output file')
    p.add_argument('--lopt-levels', default=None)
    p.add_argument('--energies', default='refit', choices=['refit', 'lopt'])
    p.add_argument('--energies-csv', default=None,
                   help='csv of revised adopted energies (level_id,E_input)')
    p.add_argument('--out', default='level_positions.csv')
    p.add_argument('--fit-hfs', action='store_true',
                   help='re-fit one hyperfine width per level on the 1974 '
                        'lines of the run, write level_hfs_widths.csv and '
                        'stop; nothing else is scanned')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.drop_questionable and not args.unknown:
        raise SystemExit('--drop-all-questionable only has a meaning with '
                         '--unknown: it frees the lines of the questionable '
                         'levels for a level being searched for, and a report '
                         'written with those lines free would not be a report '
                         'of this run')
    if args.fit_hfs:
        output_files.require_writable([HFS_FILE], 'hyperfine width file')
        if args.lopt:
            ls.LOPT_LINES = args.lopt
            ls.LOPT_LEVELS = args.lopt_levels
            ls.ENERGIES = args.energies
        if args.energies_csv:
            ls.E_INPUT_CSV = args.energies_csv
        levels, real, _ = ls.read_run()
        acc = real[real['accepted'] == 1].copy()
        acc['low_id'] = acc['low_id'].astype(str)
        acc['upp_id'] = acc['upp_id'].astype(str)
        w = fit_hfs_widths(acc, levels)
        w.to_csv(HFS_FILE, index=False, lineterminator='\n')
        used = w[w['w_applied'] > 0]
        marked = used[used['w_applied'] > HFS_MARK]
        n_cfg = int((w['w_source'] == 'configuration').sum())
        print(f"\n{HFS_FILE}: {len(w)} levels, {len(used)} carry a width "
              f"({n_cfg} of them their configuration's, the fit having too "
              f"few lines of their own), {len(marked)} marked as having "
              f"a large hyperfine structure")
        print(marked.sort_values('E_input').to_string(index=False))
        return 0
    if not args.detail and not args.unknown:
        output_files.require_writable(output_files.with_twin(args.out),
                                      'report file')
        if args.firm and (args.scan or args.audit):
            output_files.require_writable([args.firm], 'registry file')
    ctx = build(args)

    if args.detail:
        print_detail(ctx, args.detail, args.at, args.detail_all,
                     args.alt_drop)
        return 0

    if args.unknown:
        if args.drop_questionable:
            release_levels(ctx, questionable_levels(args.out), args)
        en, trans, mapping, e_meas, per_cfg, whole, id_of = unfound_theory(ctx)
        for idx in args.unknown:
            if idx not in en.index:
                print(f"\nIDEN2 row {idx} is not in enlev.dat")
                continue
            n_pred = register_unknown(ctx, idx, en, trans, mapping, e_meas,
                                      id_of)
            if not n_pred:
                print(f"\nIDEN2 row {idx} has no calculated transition to "
                      f"any level that has been found: there is nothing to "
                      f"search on")
                continue
            if args.at is not None:
                print_detail(ctx, UNKNOWN_ID, args.at, args.detail_all,
                             args.alt_drop)
                continue
            r = scan_unknown(ctx, idx, en, per_cfg, whole, step=args.step,
                             min_ln_r=args.min_ln_r, k=args.window_sigmas)
            print_unknown(ctx, idx, en, r, n_pred, top=args.top)
        return 0

    ids = args.levels if args.levels else list(ctx.per['level_id'])
    # a level with no predicted transition cannot be scanned; it is named,
    # not dropped in silence - that is how every level of files.new_levels
    # once went missing from the audit without a word
    dropped = [i for i in ids if i not in ctx.by_level]
    ids = [i for i in ids if i in ctx.by_level]
    if dropped:
        print(f"\n{len(dropped)} level(s) have no predicted transition and "
              f"cannot be considered: {', '.join(map(str, dropped))}")
    tab = report_table(ctx, ids)

    if args.scan or args.audit:
        en = li.read_enlev()
        win = li.configuration_windows(en)
        lv, unmatched = li.attach_identities(ctx.per, en, win,
                                             row_of_id=li.id_rows())
        e_calc = dict(zip(lv['level_id'], lv['E_calc']))
        w_of = dict(zip(lv['level_id'], lv['W']))
        registry = read_firm(args.firm)
        known = registry if args.use_firm else {}
        if registry and not args.use_firm:
            print(f"\nthe registry {args.firm} calls {len(registry)} levels "
                  f"settled; --use-firm skips their scans")
        glob = run_digest(ctx)
        acc_wn, acc_key = accepted_index(ctx)
        print(f"\n{len(ids)} levels to consider "
              f"({len(unmatched)} matched no row of enlev.dat)"
              + (f"; the registry may spare some of them their scan"
                 if known else ''))
        n_alt, best_alt, d_alt, sep, width = [], [], [], [], []
        scanned, fresh = [], {k: dict(v) for k, v in registry.items()}
        # a maximum that falls on another level of the run OF THE SAME J
        # AND PARITY is a different question from one that falls on empty
        # energy: the first is an interchange, which level_interchange.py
        # judges on evidence this scan does not use, the second a position
        # nobody has claimed.  A level of a different J or parity sitting at
        # the same energy is neither - it has different transitions, so it is
        # not competing for these lines and the two identities cannot be
        # exchanged - and it is not named at all.
        kinds, kind_of = level_kinds(ctx.per)
        alt_lev, alt_lev_dE = [], []
        for lid, ln_r in zip(tab['level_id'], tab['ln_R']):
            e_a, ln_r = ctx.e_final[lid], float(ln_r)
            ec, wd = e_calc.get(lid, np.nan), w_of.get(lid, np.nan)
            fp = scan_fingerprint(ctx, lid, e_a, ec, wd, args.step,
                                  args.alt_drop, glob, acc_wn, acc_key)
            was = known.get(str(lid))
            if (was is not None and str(was['fingerprint']) == fp
                    and abs(ln_r - float(was['ln_R'])) <= FIRM_TOL):
                # nothing this level's scan depends on has moved, and ln R at
                # the adopted position still comes out what it did: the scan
                # would find the same empty window it found last time
                scanned.append('registry')
                n_alt.append(0)
                width.append(float(was['scan_width']))
                best_alt.append(np.nan)
                d_alt.append(np.nan)
                sep.append(np.nan)
                alt_lev.append('')
                alt_lev_dE.append(np.nan)
                fresh[str(lid)] = dict(level_id=lid, E=e_a, ln_R=ln_r,
                                       scan_width=float(was['scan_width']),
                                       fingerprint=fp)
                continue
            scanned.append('this run')
            r = scan_level(ctx, lid, e_a, ec, wd,
                           step=args.step, alt_drop=args.alt_drop)
            if not r['alternates'] and ln_r >= FIRM_LN_R:
                fresh[str(lid)] = dict(level_id=lid, E=e_a, ln_R=ln_r,
                                       scan_width=r['hi'] - r['lo'],
                                       fingerprint=fp)
            else:
                fresh.pop(str(lid), None)
            n_alt.append(len(r['alternates']))
            width.append(r['hi'] - r['lo'])
            if r['alternates']:
                e_a, v_a = r['alternates'][0]
                best_alt.append(v_a)
                d_alt.append(r['ln_R'] - v_a)
                sep.append(e_a - r['e_adopted'])
                nid, ndE = nearest_same_kind(kinds, kind_of, lid, e_a)
                alt_lev.append(nid)
                alt_lev_dE.append(ndE)
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
        tab['scanned'] = scanned
        n_cached = scanned.count('registry')
        print(f"  windows actually scanned: {len(scanned) - n_cached}"
              + (f"; the other {n_cached} were taken from the registry "
                 f"unchanged" if n_cached else ''))
        if args.firm:
            write_firm(args.firm, list(fresh.values()))
        if args.audit:
            tab = audit_table(ctx, tab, e_calc, w_of, vacancies(en))

    tab = order_report(tab)
    mc.save_table(tab, args.out, decimals={
        'E': 4, 'ln_R': 3, 'ln_R_match': 3, 'ln_R_miss': 3,
        'sum_lnG': 3, 'ln_R_alt': 3, 'd_ln_R': 3, 'dE_alt': 4,
        'near_dE': 4, 'scan_width': 1, 'gain': 3, 'look': 3,
        'free_gain': 3, 'top_share': 3, 'kept_light': 3, 'z_alt': 2,
        'ln_J': 3, 'ln_J_alt': 3, 'd_ln_J': 3, 'gain_R': 3,
        'ln_R_J': 3, 'ln_R_alt_J': 3})
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
        onlev = (q['near_dE'].abs() < INTERCHANGE_DE).sum()
        print(f"  {int(onlev)} of the best alternates fall on another level "
              f"of the run of the same J and parity - an interchange, not a "
              f"free position")
        cols = ['level_id', 'E', 'n_obs', 'n_seen', 'n_miss', 'ln_R',
                'n_alt', 'ln_R_alt', 'd_ln_R', 'dE_alt', 'near_level',
                'near_dE']
        with pd.option_context('display.width', 200,
                               'display.max_columns', 24):
            print(q.sort_values('d_ln_R')[cols].head(40)
                  .to_string(index=False, na_rep='-'))
        # An empty action is not a clean bill of health.  It means the scan
        # found no alternate, which is reassuring only when the adopted
        # position is itself well supported; when it is not, the level has
        # nothing holding it where it is AND nowhere better to go, and it is
        # the hardest case in the report rather than the safest.  The audit
        # cannot show these - it lists levels that have an alternate - so
        # they are listed here.
        stuck = tab[(tab['n_alt'] == 0) & (tab['ln_R'] < FIRM_LN_R)]
        print(f"\n{len(stuck)} levels have no alternate anywhere in the "
              f"window AND ln R below {FIRM_LN_R:g} where they stand.  They "
              f"carry no action, because there is nothing to propose moving "
              f"them to;\n  that is not a clean bill of health.  The lines do "
              f"not support them where they are and the scan has nothing "
              f"better to offer, so what is wrong with them - if\n  anything "
              f"is - lies outside the window the calculation allows, or in "
              f"the lines assigned to them:")
        if len(stuck):
            with pd.option_context('display.width', 200,
                                   'display.max_columns', 24):
                print(stuck.sort_values('ln_R')[
                    ['level_id', 'E', 'n_pred', 'n_obs', 'n_seen', 'n_miss',
                     'ln_R', 'scan_width']].to_string(index=False,
                                                      na_rep='-'))
    if 'action' in tab.columns:
        print_audit(tab, args.alt_drop)
    with pd.option_context('display.width', 200,
                           'display.max_columns', 24):
        print("\nthe twenty weakest positions:")
        print(tab.sort_values('ln_R').head(20)
              .to_string(index=False, na_rep='-'))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
