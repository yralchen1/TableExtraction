# Contradictions and remaining deficiencies - draft for the paper's discussion

Draft of 2026-11-21, on the state of finalization round 10 (the 4f²5d scale
of the calculated A constants applied; Set 3 in `final/`). Nothing listed here
is acted on in the line or level lists unless an item says so. The user
reviews this draft; it is not yet text for the paper.

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
Two flags disagree with the sign of D:

| Line (cm⁻¹) | Flag | D (cm⁻¹) |
|---|---|---|
| 36667.3188 | `*r` | −0.0019 |
| 29151.1276 | `*r` | −0.0768 |

36667.32's pattern is negligibly small (|D| = 0.002 cm⁻¹, far below the line
width), so its flag cannot come from this transition's structure; a blend or a
misreading of the shading is the likelier cause. 29151.13 is also the worst
case of item 2 (component scale s = −15).

## 2. Printed hyperfine components that disagree with the head position

Sugar printed component positions for 228 flagged lines. For each, the
pattern's cg found from the printed components is compared with the cg found
from the head and D. The median difference is 0.009 cm⁻¹; 7 lines disagree
beyond 3σ (there were 20 in October 2026, before Step 4's A constants):

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
(s = 1: agreement). To be written for each: the suspected level(s), whether
the A constant or the assignment is in doubt, and how the contradiction could
be resolved. The two single-component cases with |s| far from 1 (29151.13,
29753.59) are most simply printed components of some other transition.

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
with the note "err in Sugar". Here the calculated pattern has a single maximum
(the upper level's structure fills the gap between the groups), and the
F = 2 group lies at +0.06, not +0.125. It is left as it is; whether it is a
companion is for the user.

## 8. Shifts of the level energies from Sugar's published values

A level's energy differs from Sugar's for three reasons, separated by three
LOPT fits of Set 3's input: (a) the wavelength calibration of the lines, (b)
the hyperfine correction to the cg, (c) the rest: changed assignments, new
lines, blends, and Sugar's own fit (`Devel/calib_level_shifts_2026-11-21`).
Over the 562 levels Sugar published that are neither new, revised nor
discarded, rms 0.20 cm⁻¹ in all: calibration 0.08, hyperfine 0.14, rest 0.13.

| Level | Shift | Calibration | Hyperfine (cg) | Rest | The user's account |
|---|---|---|---|---|---|
| 000493 | +1.453 | +0.165 | +0.056 | +1.232 | No classified lines in Sugar's lists. His position was probably based on a few lines, all now found heavily blended by transitions from other levels; only 4 of its 15 assigned lines belong to it alone. Nothing alarming. |
| 000300 | −0.914 | +0.118 | +0.043 | −1.075 | Unexplained by the lines: his 5 lines, all accepted now, and the 2 new ones (53175.40, 47789.63) cannot move it by more than about 0.05. With Sugar's own uncalibrated wavenumbers all 7 agree on 84429.91(10), 1.07 below his 84430.98. Probably a misprint by 1 cm⁻¹ (84429.98). |
| 000457 | +0.321 | +0.153 | −0.521 | +0.688 | No classified lines in Sugar's lists; his position cannot be explained by the observed lines. |
| 000455 | +0.042 | +0.189 | −0.567 | +0.421 | No classified lines in Sugar's lists; his position cannot be explained by the observed lines. |
| 000488 | −0.408 | +0.072 | +0.043 | −0.524 | No classified lines in Sugar's lists. Probably based on four lines, all heavily blended by transitions from other levels. |
| 000461 | −0.393 | −0.079 | +0.057 | −0.371 | No classified lines in Sugar's lists; his position cannot be explained by the observed lines. |

The present positions are well supported by the assigned lines and are
considered firm (user, 2026-11-21). The calibration accounts for 0.07–0.19
cm⁻¹ of these shifts, not for their size. Four more unmarked levels have a
"rest" above 0.45 cm⁻¹ and are not yet reviewed: 000305 (−0.75), 000506
(+0.68), 000450 (+0.54), 000310 (−0.47).

## 9. Close pairs of lines are measured too close together

66 pairs of Sugar's lines closer than one effective line width have both lines
identified, with different transitions (51 of 1969, 15 of 1974). Each line
fits its own Ritz wavenumber (rms z 1.01 over the 132 lines, none beyond 3σ),
so they are genuinely two lines each, resolved by Sugar down to about half the
modeled width. But in 53 of the 66 pairs the measured separation is smaller
than the Ritz separation: by 0.20(3) cm⁻¹ on average, 7σ. The shortfall grows
as the lines approach each other: 27% of the separation below 0.7 FWHM, 17%
at 0.7–0.85, 7% at 0.85–1.0, in both eras (1974 0.096(23), 1969 0.229(37)
cm⁻¹). This is the pull of two partly resolved images toward each other. Each
line of such a pair is displaced toward its partner by about half the
shortfall, 0.05–0.18 cm⁻¹, which its uncertainty does not include. Not yet
acted on.

---

*Sources:* `final/set3_report.txt` (items 1, 2, 4); README, "000229 and the
components of 30879.98" (3); Work_on_hfs_plan.md item 3c and
`Devel/calc_A_scale_2026-11-21` (5); `final/merged_lines_log.txt` (6, 7);
`Devel/calib_level_shifts_2026-11-21/levels_decomposed.csv` (8);
`Devel/sugar_duplicates_2026-11-20/close_pairs.csv` with
`final/LOPT_output_lines.txt` (9).
