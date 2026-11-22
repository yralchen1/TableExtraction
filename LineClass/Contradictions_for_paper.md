# Contradictions and remaining deficiencies - draft for the paper's discussion

Draft of 2026-11-21, on the state of finalization round 10 (the 4f²5d scale
of the calculated A constants applied; Set 3 in `final/`). Nothing listed here
is acted on in the line or level lists unless an item says so. The user
reviews this draft; it is not yet text for the paper.

**The user's review (2026-11-21).** Items 1 and 2 to be acted on: D
uncertainties added to item 1, whose two lines leave the hyperfine fits;
item 2's seven cases analyzed (below), with proposals awaiting decision.
Items 3-6, 8 and 9 need no action beyond the discussion in the paper. Item 7:
the treatment of the doublet near 15505 cm⁻¹ is accepted as correct.

Terms used throughout:

- **A**: the magnetic-dipole hyperfine constant of a level, cm⁻¹. ¹⁴¹Pr has
  nuclear spin I = 5/2, so a level of angular momentum J splits into
  sublevels F = |J − I| ... J + I.
- **Center of gravity (cg)**: the intensity-weighted mean position of a
  line's hyperfine pattern, or of a level's sublevels. The published level
  energies are cg energies.
- **Head, D**: the strongest component of a pattern, at the extreme F; D is
  its displacement from the cg, D = I(A_upp J_upp − A_low J_low).
- **κ (kappa)**: where between the cg (κ = 0) and the head (κ = 1) Sugar's
  reading of an unresolved line lies, fitted per class of line.
- **Effective line width**: the instrumental width (0.035 Å) and the Doppler
  width combined; FWHM, full width at half maximum.
- **Broadened pattern**: the calculated components, each given the effective
  line width, summed. Its maxima are where a reader would set the line.
- **z**: a difference divided by its combined uncertainty.

---

## 1. Flags shaded against the sign of D

Sugar's flag `*r` says the pattern is shaded to the red (longer wavelengths),
`*v` to the violet. A pattern shaded to the red has its head at the violet
end and trails off to the red, so with the adopted A constants a `*r` line
should have D > 0 (head above the cg in wavenumber) and a `*v` line D < 0.
Two flags disagree with the sign of D. The uncertainty of D, u_D, is
I·sqrt((J_upp u_A,upp)² + (J_low u_A,low)²) from the two levels' A
uncertainties (as the model scales them):

| Line (cm⁻¹) | Flag | Transition | D (cm⁻¹) | u_D (cm⁻¹) | D/u_D | Line width (FWHM) |
|---|---|---|---|---|---|---|
| 36667.3188 | `*r` | 000130-000224 | −0.0019 | 0.027 | −0.1 | 0.56 |
| 29151.1276 | `*r` | 000151-000220 | −0.0768 | 0.048 | −1.6 | 0.38 |

Neither sign is established: 36667.32's D is zero within 0.03 cm⁻¹, and
29151.13's is 1.6 u_D from zero. Both patterns are also far narrower than
the line width, so neither transition can produce a visible shading; the
flag must come from something else on the line. For 29151.13 that is found
(item 2, case a): the head of another transition, 000187-000444, lies 0.57
cm⁻¹ to the red, where Sugar printed his one component. For 36667.32 nothing
is known.

**To be acted on (user, 2026-11-21): both lines are taken out of the
hyperfine fits.** Their flag is no evidence of where on the pattern they
were read. Proposed means: a `kappa_exceptions.txt` row of class `cg` for
each (a line narrower than its width is read at its cg), which also removes
it from every kappa fit. 29151.13's printed component is already left out
of the A fit (`hfs_A_fit.SET_ASIDE`, 2026-10-07); 36667.32 has none.

## 2. Printed hyperfine components that disagree with the head position

Sugar printed component positions for 228 flagged lines. For each, the
pattern's cg found from the printed components is compared with the cg found
from the head and D (`make_set3.component_check`). It assumes that the
tabulated line is the head of the transition with the largest share and that
the k-th printed component is the k-th rung. The median difference is 0.009
cm⁻¹; 7 lines disagree beyond 3σ (there were 20 in October 2026, before Step
4's A constants):

| Line | Flag | n comp. | D | s | difference (cm⁻¹) | z |
|---|---|---|---|---|---|---|
| 29151.1276 | `*r` | 1 | −0.0768 | −15.19 | −1.243 | 20.6 |
| 29753.5924 | `*v` | 1 | −0.7906 | 1.88 | +0.695 | −18.5 |
| 32105.4086 | `*v` | 3 | −0.8529 | 0.90 | −0.090 | 5.6 |
| 32294.9185 | `*r` | 1 | +0.6184 | 0.61 | +0.239 | −5.3 |
| 28477.5676 | `*v` | 5 | −0.7556 | 0.95 | −0.039 | 4.4 |
| 33185.8488 | `*r` | 1 | +0.3929 | 1.40 | −0.158 | 3.6 |
| 31945.1479 | `*v` | 3 | −0.6336 | 1.07 | +0.045 | −3.0 |

s is the scale the printed spacings ask for relative to the calculated pattern
(s = 1: agreement); z uses only the reading precision of a component (0.03
cm⁻¹), not the uncertainty of the A constants.

**Analyzed 2026-11-21** (round 15; `Devel/review_2026-11-21/contra_12.py`,
`blend_profile.py`, `alt_search.py`). Each line was modeled as Sugar saw it:
every transition accepted on it, at its own Ritz cg (Set 3's energies) and
weighted by its calculated intensity, each component given a Gaussian width,
summed; the maxima of that profile, relative to the measured line, were
compared with the printed positions. Widths of 0.20 and 0.28 cm⁻¹ were used
(Sugar resolved these components, so his effective width was below the
model's 0.37-0.47 cm⁻¹, which merges them); the conclusions hold for both.
For the single components, every transition of `Icalc.xlsx` between known
levels within 2.5 cm⁻¹ was searched.

None of the seven asks for a different A constant. Two of them (cases a, b)
are not components of the assigned transition; three are blends, where the
check's premises fail (c, d, e); in two the reading of the line itself, or of
a merged pair of rungs, is the cause (f, g).

a. **29151.1276 (000151-000220): the component is the head of another
   transition.** The assigned pattern is 0.19 cm⁻¹ wide in all, so a
   component 0.58 below the line cannot belong to it. 000187-000444
   (J 1.5-0.5, calculated intensity 3552, assigned to no line) has its cg at
   29151.207 and D = −0.640(47) cm⁻¹, so its head lies at 29150.567(50);
   the printed component, calibrated, is at 29150.559. In the summed profile
   it makes a separate maximum at −0.571 (height 0.11 of the line's). It
   would also explain the `*r` flag of item 1: the red side of the feature
   carries a head. 000444 is well established (10 accepted lines, observed
   to calculated intensity 0.9-2.4, all within 0.06 cm⁻¹ but one blend), and
   000187-000444 is among its strongest unobserved transitions.
   Already left out of the A fit since 2026-10-07.
b. **29753.5924 (000182-000274): no explanation found.** The assigned
   pattern's profile has maxima at −0.06 (head), +0.58 (0.70-0.75 of the
   head), +1.09 to +1.10 (0.53-0.59) and +1.49; the printed +1.188 is 0.09
   above the third, and the second, nearly as strong as the head, is not
   printed. The line is also 14 times stronger than the transition's
   calculated intensity (6451 against 466). No other transition between
   known levels lies within 2.5 cm⁻¹ of the line or 1.5 of the component.
   The reading of 10/7 (the component is rung 2, rung 1 unlisted) is
   possible but not convincing. Already left out of the A fit.
c. **32105.4086, a blend (000139-000219 with 0.70 of the intensity,
   000142-000432 with 0.30).** Their heads lie 0.44 and 0.25 cm⁻¹ below the
   line. The summed profile has its strongest maximum at −0.40 to −0.43 and
   further maxima at +0.04, +0.45, +0.81 to +0.82 and +1.12 to +1.13: the
   tabulated line is the +0.04 maximum, not a head, and the printed +0.442,
   +0.792, +1.102 match the next three maxima within 0.003-0.03. The A fit
   skips blends. The line is inflated to 0.35 cm⁻¹ already.
d. **32294.9185, a blend (000149-000238, 000165-000258, about equal).** The
   check took 000165-000258 (share 0.508 against 0.491). The printed −0.250
   is rung 1 of 000149-000238 read from the line (−0.2505, s = 0.999); in the
   summed profile it is the strongest maximum (−0.23 to −0.25), where the
   head of 000165-000258 (−0.241) and that rung meet. Skipped by the A fit.
e. **28477.5676, a blend (000189-000284 with 0.81, 000148-000211 with
   0.19).** The profile's maxima at +0.37, +0.72 to +0.73, +1.05, +1.34 to
   +1.36, +1.57 to +1.60, measured from its head maximum (−0.03 to −0.04),
   reproduce the printed +0.332, +0.695, +1.053, +1.321, +1.595 within
   0.01-0.03. Skipped by the A fit.
f. **33185.8488 (000150-000250): merged rungs.** Rung 1 (−0.270) is closer to
   the head than any width allows to resolve, and rungs 2-5 cluster at
   −0.49 to −0.67. The profile's first maximum below the head lies at −0.29
   to −0.37 depending on the width, against the printed −0.378: the
   component is rungs 1 and 2 read together, not rung 1. Read as rung 1, it
   would make A(000250) = +0.042(8) against the calculated +0.024(3); the
   table keeps the calculated value, so no adopted constant rests on it.
g. **31945.1479 (000133-000206): the line, not the components.** From the
   Ritz head (0.072 above the measured line) the rungs fall at +0.454,
   +0.784, +1.063, against the printed +0.455, +0.753, +1.051 (within 0.03).
   The components agree with the adopted A constants; it is the line, the
   head, that lies 0.075 cm⁻¹ below its Ritz value (LOPT O−C −0.076). Since
   the A fit measures every component from the line, this pattern pulls
   A(000206) down: left out, A(000206) (combined with its composition prior)
   goes from +0.0287(13) to +0.0306(14), and A(000133) from +0.0800 to
   +0.0794(12).

**Decided (user, 2026-11-21).**
- 36667.32 and 29151.13 (item 1): `kappa_exceptions.txt` class `cg`.
- 33185.85: its component enters the A fit as rungs 1 and 2 read as one,
  their strength-weighted mean (`hfs_A_fit.MERGED_RUNGS`; −0.358 predicted
  against the printed −0.378). No adopted constant changes.
- 31945.15: first left out of the A fit, then put back. Left out, A(000206)
  rose to +0.0306(14), and the calibration's spread of kappa for 1974 rose
  from 0.079 to 0.090 (with that constant alone reverted, 0.079): the
  positions of 000206's plain lines favor about +0.0287, the value the fit
  gives with 31945.15 in. Its single-pattern values: 29212.40 +0.0299 to
  +0.0304, 29918.67 +0.0311, 31945.15 +0.0257. The reading of case g (the
  line read low) assumed A(000206) = +0.0287 for the Ritz head, so it does not
  decide between them. A(000206) stays +0.0287(13).
- 29150.559 (case a): a line of `Pr3_lines.xlsx` (Icor 2470, Iorig 1.57, the
  rung-1 rule: 0.7864 × the main line's), assigned to 000187-000444 by the
  ledger for the user's review in IDEN2, held at the head (`kappa_exceptions`
  `head`; at the wide class's kappa 0.875 its O−C was −0.083). In iter_hfs's
  LOPT: grade 3A, Iobs/Icalc 0.70; 000444 and 000187 move by less than 0.005.
  29151.13 then fits at O−C +0.029 (was +0.099).
- 29754.785 (case b): an unidentified line of `Pr3_lines.xlsx` (Icor 3519,
  Iorig 1.64; 0.5455 × the main line's).
- Both printed components are dropped from `hfs_components.csv`
  (`hfs_components.DISMISSED`); `final/merged_lines_log.txt` section C.
- The blends (c, d, e): no action. `make_set3.component_check` could model
  them as above, so that the report lists only real disagreements.

## 3. 000229 and the printed components of 30879.98

Level 000229 (J = 5/2) had two patterns that disagreed. The two components
Sugar printed for 30879.98 `*r` (000150–000229), at −0.430 and −0.820 cm⁻¹
from the head, asked for A(000229) = +0.063(4). The one component of the
seven times stronger 33021.15 `*v` (000229 with another level) gave
+0.024(9), in line with the composition value +0.0286(61) and the eigenvector
calculation +0.0280(40): 3.8σ apart. No identified level has a transition at
the two component positions. The user dismissed them as misidentified
(2026-10-09); A(000229) is now +0.0270(51). They remain lines of the list:
30879.549 became a companion of 30879.98 (Sugar printed it as rung 1 of the
pattern; it lies on the broadened pattern's maximum), and 30879.159 stays a
line of its own.

## 4. Lines of the resolved levels 000190 and 000642 without an F = 2 partner

000190 and 000642 have J = 1/2, so each splits into two sublevels, F = 3 and
F = 2, far enough apart (3A = 1.5 cm⁻¹ for 000642) that each of their
transitions should appear twice. For 14 lines of these levels the F = 2
partner is not in Sugar's list (3 more have their partner among the
companions). For 000190 the F = 2 groups are among Sugar's printed hfs
components, which are not lines of the list. Part of the general discussion of
predicted strong lines that are missing.

## 5. The calculated A constants

The A constants of 527 levels are calculated from the eigenvectors of the
Cowan least-squares fit (LSF) and radial parameters fitted to measured
constants. Two deficiencies:

- **The accuracy of the LSF limits them.** Its eigenvectors set every
  calculated A and every predicted intensity; the user will not refine the LSF
  in this paper.
- **The 4f²5d constants were 13% too large.** No measured A tests the
  radial parameters of 4f²5d, the leading configuration of 104 calculated
  levels. Measured on the line positions (each calculated 4f²5d A multiplied
  by a trial factor t and the calibration refitted): t = 0.872(39), 3.3σ
  from 1, while the other 423 calculated levels give 1.044(34). Adopted
  2026-11-21 as a scale of 0.87(4) on those constants, its uncertainty added
  to each u_A. With it, the difference between the κ of plain 1974 and 1969
  lines on the 410 levels both eras see is again significant (5.6σ against
  2.3σ without the scale), as it was on the measured constants of September
  2026. Not fully explained: the 1969 lines alone ask for t = 0.66(10), the
  1974 lines for 0.905(42).

## 6. 34396.02: two halves of a pattern, possibly with a blend

Sugar listed two lines, 34395.48 (I 2) and 34396.56 (I 7), symmetrically
placed at −0.54 and +0.54 cm⁻¹ about the cg of 000148–000255, whose
broadened pattern is flat-topped (four maxima of nearly equal height, 2.2
cm⁻¹ wide). They were merged into one line at their mean position
(`final/merged_lines_log.txt`, entry 4). The pattern is nearly symmetric, but
the intensities are not: the stronger half may carry a blend with an
unidentified line. No candidate transition was found. Documented, not acted
on (user, 2026-11-20).

## 7. The doublet at 15505 cm⁻¹ (and 13595): agreement after all

*The user asked for this case to be explained; it was first entered on this
list as evidence that a calculated A is too large, and that reading is
withdrawn.*

**The lines.** Sugar listed 15504.946 (I 5) and 15505.047 (I 4), 0.101 cm⁻¹
apart (0.042 Å at 6449 Å). His classification 000065–000143 was on 15505.05;
the user moved it to 15504.95, which agrees better with Sugar's own Ritz
wavenumber (the note "err in Sugar"). On 2026-11-20 15505.05 was made a
companion of 15504.95 (`final/merged_lines_log.txt`, entry 5).

**Why a doublet is expected.** Both levels have J = 1/2, so each has two
sublevels, F = 2 and F = 3. The lower level 000065 (4f³) has a large
calculated A, +0.0531(34) cm⁻¹; the upper 000143 a small one, +0.0069
(calculated, 4f²5d, scaled). The pattern is therefore two groups of
components, one from each sublevel of 000065, about 3·A(000065) ≈ 0.16
cm⁻¹ apart:

| Component (F_low → F_up) | Position from cg (cm⁻¹) | Share |
|---|---|---|
| 3 → 2 | −0.078 | 0.32 |
| 3 → 3 | −0.058 | 0.26 |
| 2 → 2 | +0.081 | 0.09 |
| 2 → 3 | +0.102 | 0.32 |

**What was claimed, and why it was wrong.** Taking Sugar's two lines as the
centers of the two groups, their separation 0.101 asks for A(000065) ≈ +0.034,
a third below the calculated value. But the groups overlap: the effective line
width at 6449 Å is 0.153 cm⁻¹, more than their separation. Broadened to that
width the pattern shows two maxima only if A(000065) ≥ 0.050, and their
separation then runs from 0.118 cm⁻¹ (A = 0.050, just resolved) to 0.145
(A = 0.053, the calculated value; maxima at −0.063 and +0.082, heights 1.00
and 0.75). No A gives a resolved doublet only 0.101 apart at this width. A =
0.034 would need a width of 0.09 cm⁻¹, below the Doppler width alone (about
0.13 cm⁻¹ at this wavelength).

**Conclusion.** Sugar's readings are 0.017 cm⁻¹ closer than the nearest
resolvable doublet, about what the readings of a barely resolved pair allow
(see item 9, where close pairs are found measured too close together). The
doublet agrees with the calculated A(000065) = +0.053(3); it also confirms
that the identification belongs on the stronger, red group (calculated height
1.00 against 0.75; observed intensities 5 and 4), where the user put it. The
summed intensity of the two lines matches the calculated one (ratio 0.93).

**13595.21/13595.34**, a transition from the same lower level (000065–000138,
upper J = 3/2), has a second line 0.131 cm⁻¹ above the classified one, also
with the note "err in Sugar" on both, as at 15505. At the modeled width and
the adopted A values the calculated pattern has a single maximum (the upper
level's structure, A(000138) = +0.0112, scaled, tier 3, fills the gap between
the groups). A second maximum, the F = 2 group of 000065, appears only with
A(000138) near the low end of its uncertainty and a width below the modeled
one, and then 0.09-0.13 cm⁻¹ from the first against 0.131 observed; its
height, 0.74 of the first, matches Sugar's 15/20 = 0.75, as the 15505
doublet's 4/5 matches its 0.75. A possible companion, not a proven one.

## 8. Shifts of the level energies from their published values

A level's energy differs from the published one for three reasons, separated
by three LOPT fits of Set 3's input: (a) the wavelength calibration of the
lines, (b) the hyperfine correction to the cg, (c) the rest: changed
assignments, new lines, blends, and the published fit itself
(`Devel/calib_level_shifts_2026-11-21`; the line-by-line check in
`Devel/review_2026-11-21/level_review.py`). Over the 562 published levels
that are neither new, revised nor discarded, rms 0.20 cm⁻¹ in all:
calibration 0.08, hyperfine 0.14, rest 0.13.

**Whose value.** Only three of the ten levels below are Sugar's (`L3974` in
`TableExtraction/Pr3_lev_Wyart_1999.xlsm`, the levels in NIST ASD): 000300,
000305, 000310. The other seven are Wyart's unpublished levels (`u#`), whose
energies were never finalized; that is why none of them has a classified line
in Sugar's lists. None of the ten has a verdict in `level_verdicts.txt` or a
mark of the audit (2026-11-19).

"Calibration", "hyperfine" and "rest" are measured in that order, so a strong
flagged line read at its head counts in "rest" until the hyperfine step; the
last column says what the lines show in Set 3's own terms (calibrated
wavenumbers at the cg, the other levels at their Set 3 energies).

| Level | Source | Shift | Calibration | Hyperfine (cg) | Rest | Account |
|---|---|---|---|---|---|---|
| 000493 | Wyart | +1.453 | +0.165 | +0.056 | +1.232 | No classified lines in Sugar's lists. Wyart's position was probably based on a few lines, all now found heavily blended by transitions from other levels; only 4 of its 15 assigned lines belong to it alone. Nothing alarming (user). |
| 000300 | Sugar | −0.914 | +0.118 | +0.043 | −1.075 | His 5 lines, all accepted now, and the 2 new ones (53175.40, 47789.63) cannot move it by more than about 0.05 (user). With Sugar's own uncalibrated wavenumbers all 7 agree on 84429.91(10), 1.07 below his 84430.98: probably a misprint by 1 cm⁻¹ (84429.98). |
| 000457 | Wyart | +0.321 | +0.153 | −0.521 | +0.688 | No classified lines in Sugar's lists; the position cannot be explained by the observed lines (user). |
| 000455 | Wyart | +0.042 | +0.189 | −0.567 | +0.421 | No classified lines in Sugar's lists; the position cannot be explained by the observed lines (user). |
| 000488 | Wyart | −0.408 | +0.072 | +0.043 | −0.524 | No classified lines in Sugar's lists. Probably based on four lines, all heavily blended by transitions from other levels (user). |
| 000461 | Wyart | −0.393 | −0.079 | +0.057 | −0.371 | No classified lines in Sugar's lists; the position cannot be explained by the observed lines (user). |
| 000305 | Sugar | −0.826 | +0.037 | −0.118 | −0.746 | Inherited from its lower levels. Sugar classified 3 lines (1969: 50500.18 h, 53657.57, 58215.01), all accepted now and unblended; they go to 000133, 000148 and 000161, 4f²6s levels with large hyperfine structure (S 1.05-1.10) whose energies moved −0.69 to −0.76 to their cg. In Set 3's terms his 3 lines put 000305 at −0.69(12) from his value and its 12 new lines (to other levels; 24906.60 the most precise) at −0.85(3): they agree. With uncorrected wavenumbers they would disagree by 0.7. |
| 000506 | Wyart | +0.531 | −0.164 | +0.013 | +0.683 | No classified lines in Sugar's lists. Three lines now (96636.63, 93436.48 blended, 86524.12), implying +0.42(10) in Set 3's terms (spread −0.03 to +0.56); the position cannot be explained by the observed lines, and the level rests on few lines. |
| 000450 | Wyart | +0.770 | −0.019 | +0.248 | +0.541 | No classified lines in Sugar's lists. Its 9 lines agree on +0.74(4) in Set 3's terms. Wyart's value matches the uncorrected wavenumber of one of them, 52321.75 (to 000142), within 0.03: a line held at the head of a wide hyperfine pattern (000142, S = 0.81; its cg shift +1.04). Read without the hyperfine correction, that line gives Wyart's value. |
| 000310 | Sugar | −0.091 | +0.103 | +0.275 | −0.468 | No real shift: in Set 3 it is 0.09 below Sugar's value. His 9 lines, all accepted, put it at −0.05(4), the 6 new lines at −0.09(2). It appears here because the most precise new line, 24593.42 `*r`, is read at its head, 0.57 from its cg, in the uncorrected fit. |

The calibration accounts for 0.04-0.19 cm⁻¹ of these shifts, not for their
size. Two of the larger ones come from the hyperfine structure (000305 through
its lower levels, 000450 through a head reading); one is a probable misprint
(000300); 000310 is not shifted; the Wyart levels without lines of Sugar's
differ from values that were never finalized.

## 9. Close pairs of lines are measured too close together

66 pairs of Sugar's lines closer than one effective line width have both lines
identified, with different transitions (51 of 1969, 15 of 1974). Each line
fits its own Ritz wavenumber (rms z 1.01 over the 132 lines, none beyond 3σ),
so they are genuinely two lines each, resolved by Sugar down to about half the
modeled width. But in 53 of the 66 pairs the measured separation is smaller
than the Ritz separation: by 0.20(3) cm⁻¹ on average, 7σ. This is the pull of
two partly resolved images toward each other.

Measured on all 617 pairs closer than 2.5 widths with both lines identified
(`Devel/close_pair_pull_2026-11-21`), each line reads toward its partner by
on average sep · share · g(x), with x = sep/FWHM, g(x) = 0.98 exp(−x/0.37)
(0.25 at x = 0.5, 0.09 at 0.9, 0.02 at 1.5), and share the partner's part of
the two lines' summed intensity: the weaker line of a pair is pulled about
three times as far as the stronger. The model lowers χ² from 1068 to 995 over
the 1234 lines.

**Acted on (user, 2026-11-21: as an uncertainty, not a correction).** Every
line of the list now carries this amount, from its neighbors below and above
in quadrature, as a fixed part of its variance: in the calibration's weights
and class fit, and in the uncertainty it writes for the line, which the
classification, LOPT and IDEN2 take (`hfs_kappa.pair_uncertainties`). 1660 of
the 6669 lines have such a term, 795 above 0.01 cm⁻¹, the largest 0.47. On
the lines below 1.2 widths it brings the rms z from 0.97-1.42 to 0.79-1.17,
the other lines' being about 0.92.

---

*Sources:* `final/set3_report.txt` (items 1, 2, 4); README, "000229 and the
components of 30879.98" (3); Work_on_hfs_plan.md item 3c and
`Devel/calc_A_scale_2026-11-21` (5); `final/merged_lines_log.txt` (6, 7);
`Devel/calib_level_shifts_2026-11-21/levels_decomposed.csv` (8);
`Devel/sugar_duplicates_2026-11-20/close_pairs.csv` with
`final/LOPT_output_lines.txt` (9).
