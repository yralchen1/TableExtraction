# Plan for the hyperfine-structure (hfs) work in Pr III

Status: written 2026-09-12, revised 2026-09-13 after the user's review.
Nothing in the analysis pipeline has been changed yet. Three result files
exist (`A_hfs_levels.csv`, `hfs_line_corrections.csv`,
`hfs_level_shifts.csv`); they are additive, and no stage of the pipeline has
been re-run.

**The revision of 2026-09-13 changes the foundations of this work, not its
details.** Three things came in:

1. Sugar's own statement, for a neighboring spectrum taken with comparable
   instrumentation, that he measured *the centers of unresolved hfs blends*
   (section 1.1). The working hypothesis is therefore center-of-gravity
   measurement everywhere except the flagged 1974 lines - not measurement at
   the strongest component everywhere.
2. The hfs trends must be derived from **Sugar's own stated uncertainties**,
   not from the inflated uncertainties now in `Pr3_lines.xlsx`, which were
   made for Set 1 and bias any trend fitted through them (section 2.6).
3. `Pr3_lines.xlsx` is missing 501 measured wavelengths of the resolved hfs
   components of the 1974 flagged lines. They survive in the checked
   extraction workbook and have never been used (section 3.5). They are
   direct measurements of the quantity the kappa model only infers.

Until these three are worked through, **the kappa values quoted in section
2.2 and the level shifts of section 2.5 are provisional. They are not to be
published and not to be put into the pipeline.**

---

## 0. Terms used in this document

Every symbol is defined here once, and the definitions are not repeated below.

* **hfs — hyperfine structure.** The nucleus of the only stable praseodymium
  isotope, 141-Pr, has spin **I = 5/2**. Each atomic energy level of total
  electronic angular momentum **J** is therefore not a single level but a
  group of sublevels labeled by **F = |J − I| … J + I**. A spectral line
  between two such levels is not a single line but a pattern of components —
  the *hyperfine pattern*.
* **A (the magnetic-dipole hyperfine constant)** of a level, in cm^-1. The
  sublevel of quantum number F lies at
  `W(F) = (A/2)·[F(F+1) − J(J+1) − I(I+1)]`
  above the level's **center of gravity (cg)** — the intensity-weighted mean
  position of all its sublevels, which is what a level value in a published
  table is supposed to be.
* **S ≡ I·A·J.** The displacement of the outermost hyperfine sublevel
  (`F = J + I`) of a level from that level's own cg. Positive A puts it above
  the cg, negative A below.
* **D ≡ S_upper − S_lower.** The displacement, from the cg wavenumber of a
  line, of the extreme hyperfine component
  `F_low = J_low + I → F_up = J_up + I` — which for a well-separated pattern
  is the strongest component. `D > 0` means the pattern's strong end lies on
  the high-wavenumber (violet) side; `D < 0` on the low-wavenumber (red) side.
* **Sugar's flags `*r` and `*v`.** Characters in Sugar's 1974 line list
  marking a line whose hfs pattern was visibly resolved or broadened, quoted
  from the original paper: "*r or *v = strongest component of a hfs flag
  pattern where r means the pattern is shaded to the red and v means the
  pattern is shaded to the violet."
  **`*r` means `D > 0` and `*v` means `D < 0`**, and the tabulated wavenumber
  of such a line is the extreme component, not the cg. The argument: the
  tabulated wavelength is the strongest component, and "shaded to the red"
  means the rest of the pattern trails away from it toward the red, i.e. to
  lower wavenumber. The cg therefore lies below the tabulated component, so
  `measured - cg = D > 0`. `*v` is the mirror image. Verified 2026-09-21 by
  three independent routes, which is why this is no longer called an
  interpretation:
  1. *Sugar's own resolved components* (section 3.5): the companions of an
     `*r` line all lie redward of it and those of a `*v` line all violetward,
     500 rows out of 500, with mean displacements -0.56 and +0.82 cm^-1.
  2. *The computed A constants* (section 2.3), through
     `hfs_line_corrections.csv`: of the flagged lines whose `D` is non-zero,
     136 of 137 `*r` lines have `D > 0` (mean +0.47 cm^-1) and 129 of 131
     `*v` lines have `D < 0` (mean -0.79 cm^-1). The magnitudes match item 1.
  3. *One directly resolved pair* (`4f26s_2P_0.5_hfs.md` section 1): the line
     29294.708, tagged `*r`, has its partner recorded at 29294.209, 0.499
     cm^-1 redward. Its lower level 1089 has a measured `A = -0.157`, hence
     `S_low < 0` and `D > 0` with an upper level of small `A` - which is what
     `*r` asserts. This is a single line, seen, with both components on the
     plate, and it agrees.

* **κ (kappa), the measurement-convention factor.** The model
  `wavenumber measured = wavenumber of cg + κ · D`.
  `κ = 0` means the measurement was made at the center of gravity; `κ = 1`
  means it was made at the extreme (strongest) component. κ is fitted per
  class of line.
* **Line class.** `flag` = carries `*r` or `*v`; `c` = carries Sugar's
  character `c` ("complex", 1974 list); `plain` = neither.
* **obs − Ritz.** The difference between a measured line wavenumber and the
  difference of the two fitted level values it connects. The *Ritz* value is
  the prediction from the level system.
* **Set 1 and Set 2.** The two result sets proposed earlier: Set 1 reports
  Sugar's measured wavenumbers unchanged with honest uncertainties, and
  therefore levels that are *effective* (hfs-displaced, not physical); Set 2
  reports hfs-corrected wavenumbers and therefore levels that are physical
  centers of gravity.
* **θ₁, θ₂ (theta).** The two angular coefficients in Reader & Sugar's
  expression `A = θ₁·a_4f + θ₂·a_6s` for a level of the 4f²6s configuration.
* **a_nl (radial hyperfine parameters)**, in cm^-1: `a_4f = 0.029 ± 0.003`,
  `a_6s = 0.639 ± 0.002` (Reader & Sugar 1965); adopted here in addition
  `a_6p(1/2) = 0.12 ± 0.05` and `a_6p(3/2) = 0.020 ± 0.010`, both estimated,
  not measured.
* **Component wavelength.** A wavelength printed by Sugar for one resolved
  component of an hfs pattern, as distinct from the wavelength of the line
  itself. Section 3.5: there are 501 of them, in the 1974 extraction workbook,
  as the rows with `Intens = 0`.
* **`unc_own`.** The wavenumber uncertainty column of `Pr3_lines.xlsx` as it
  now stands. It was built for Set 1 and deliberately inflates the uncertainty
  of a line that deviates far from its Ritz value. Section 2.6: it must not be
  used as a weight when a trend is fitted to those same deviations.
* **Sugar's stated uncertainties.** The values Sugar himself gives in the two
  papers - 0.004 A (1969), 0.003 A for plain lines and an average deviation of
  0.007 A from Ritz for blended and complex lines (1974). These, not
  `unc_own`, are the starting weights for the hfs trend work.
* **LOPT.** The least-squares level-optimisation program that produces the
  level values and the Ritz wavenumbers from the classified line list.

---

## 1. Why this work exists

Pr III level values in this project are obtained by a weighted least-squares
fit of measured line wavenumbers. If the measured wavenumber of a line is not
the wavenumber of its center of gravity — because the measurer pointed at the
brightest hyperfine component instead — then every level value carries a
systematic displacement, and the published level list is not the physical
level list. The displacement is not small: it reaches **1.2 cm^-1** for the
worst level, against a line measurement uncertainty near 0.006 cm^-1.

The work has three goals, in order of increasing ambition:

1. **Honest uncertainties.** Whatever else is done, a line whose position is
   uncertain because of unresolved hfs must carry an uncertainty that says so,
   and that uncertainty must be reproducible by a reader from the published
   tabulation alone, with no hidden inflation factor.
2. **A correct level list.** If the hfs displacement of each line can be
   computed, it can be subtracted, and the fit then returns physical centers
   of gravity. This is Set 2.
3. **New spectroscopic information.** The flags `*r` and `*v`, read as
   inequalities on `D`, are genuine measurements of hyperfine constants of
   levels for which none were published.

---

### 1.1 What Sugar said about his own measurements, and what he did not

In none of his papers on Pr III does Sugar state whether a tabulated
wavelength is the center of gravity of an hfs pattern or the position of one
of its components. The single exception is the `*r`/`*v` lines of the 1974
list, for which he identifies the strongest components of the flag patterns.

The best external evidence is his work on the neighboring spectrum Pr IV,
made with comparable instrumentation [R3]:

> "Nearly all the lines above 2100 angstrom are unresolved hfs patterns. The
> wavelengths of the centers of these blends were measured. The uncertainty
> in these data is about +-0.01 angstrom."

**The adopted working hypothesis is therefore: Sugar measured centers of
gravity for every line except those he flagged.** It is a hypothesis, not a
law - it is evidence from a different spectrum in a different paper - and it
may have to be revised, possibly per line, possibly per era (1969 or 1974),
if the Pr III data say otherwise. Any such revision must be presented in the
publication as a revision, with the evidence that forced it.

Note what this does to the kappa model. Under this hypothesis kappa = 0 is
the *expected* value for the unflagged lines, where the provisional fit of
section 2.2 returned kappa(plain) = 0.96. Either that fit is being driven by
something other than hfs - the inflated uncertainties of section 2.6 are the
leading suspect - or the hypothesis is wrong. Resolving the disagreement is
the first scientific task of this work.

**Resolved, 2026-09-21 (Step 3, `hfs_kappa.py`).** The hypothesis is wrong for
the 1974 lines, and the weights were not what produced the 0.96. Refitted on
Sugar's own stated uncertainties, with the displacements permuted as a control
and the flagged lines anchored by the measurement of section 3.5,

    kappa(plain, 1974) = 0.944 +- 0.014      (wavenumbers below 47500 cm^-1)
    kappa(plain, 1969) = 0.633 +- 0.035      (above it)

The Pr IV sentence quoted above covers wavelengths **above 2100 angstrom**,
which is 47600 cm^-1 - the 1974 region almost exactly. There the Pr III data
contradict it: an unflagged 1974 line sits within six per cent of its
strongest component, not at the center of gravity. Below 2100 angstrom, which
is the 1969 material, the measurement does move a third of the way towards the
center of gravity, so the effect Sugar describes is real and is a
short-wavelength effect. Either way the correction applies to the unflagged
lines as well, which is what the rest of this document assumed and what the
user conjectured from the start.

## 2. What has already been established

These are results, not plans. They are the foundation the rest rests on.

### 2.1 The degeneracy theorem (proved, and verified numerically)

Write the measured wavenumber of a line as `cg + κ·D = cg + κ(S_u − S_l)`.
Then

```
obs − Ritz = (κ·S_u − ε_u) − (κ·S_l − ε_l)
```

where `ε_k` is whatever displacement the fit has already put into level *k*.
The fit absorbs a *common* part of κ exactly: adding the same constant to κ
for every class of line leaves the fit's χ² completely unchanged. This was
confirmed to the last printed decimal — χ² = 11360.77 for common shifts of
0.00, 0.25, 0.50 and 1.00.

**Consequence, and it must be stated in any publication:** the line list alone
cannot determine the absolute scale of the hfs correction. Only *differences*
of κ between classes are measurable. The absolute scale rests on one external
physical statement — that a line whose pattern Sugar could see and flag was
measured at its resolved strongest component, i.e. **κ(flag) = 1**. That is
the single assumption Set 2 adds beyond Reader & Sugar's published constants.

Section 3.5 has now removed that assumption (2026-09-21). Sugar printed the
wavelengths of the resolved components of the flagged patterns; measured
against them, the tabulated wavenumber of a flagged line is the strongest
component of the pattern with no measurable offset, so **kappa(flag) = 1 is
measured, not assumed**. What remains external is the scale of the A constants
themselves, and section 3.5 measures that too: the observed intervals are
4.2 +- 0.6 per cent smaller than the semiempirical ones.

### 2.2 The measurement convention, as fitted provisionally

**Superseded on 2026-09-21 by the refit of Step 3, which is at the end of this
section. The provisional numbers immediately below were obtained with the Set 1
uncertainties and are kept only to show what changed.**

Solving for one κ per class jointly and linearly with the level values gives

```
kappa(flag) - kappa(c)      = +0.756 +- 0.021
kappa(flag) - kappa(plain)  = +0.039 +- 0.014
kappa(c)    - kappa(plain)  = -0.717 +- 0.023
```

Anchoring `κ(flag) = 1`: **κ(plain) = 0.96 ± 0.01, κ(c) = 0.24 ± 0.02.**
An independent, slower grid scan reproduced 0.960 and 0.240.

This confirms the user's own conjecture: **Sugar measured unflagged lines at
the same feature as flagged ones.** The flags are a warning label about the
appearance of the line, not a change of measurement method. Only the `c`
("complex") lines were measured near the center of gravity, and their large
residuals in the uncorrected fit are the shadow of the other lines dragging
the levels away from cg.

**The refit of 2026-09-21 (Step 3, `hfs_kappa.py`).** Same method, new
weights: every line enters at the uncertainty adopted in section 2.6 rather
than at its `unc_own`, the ten lines the fit rejects beyond 4 sigma are set
aside, the displacements are computed with the A constants scaled by the
0.9582 that section 3.5 measured, and `kappa(flag) = 1` is now an anchor
measured in section 3.5 rather than an assumption. 4365 lines, 650 unknowns,
rank 649 - the one missing rank is the zero of the energy scale, the
degeneracy of section 2.1 having been removed by the anchor - chi2/dof 1.00 by
construction, rms 0.174 cm^-1.

```
kappa(plain) - kappa(flag)  = -0.071 +- 0.014
kappa(c)     - kappa(flag)  = -0.802 +- 0.051
anchored at kappa(flag) = 1:  kappa(plain) = 0.929 +- 0.014
                              kappa(c)     = 0.198 +- 0.051
```

**The conclusion of section 2.2 stands, and it now stands on honest weights.**
Three things were checked before saying so.

* *The displacements carry the signal.* Permuting D among the lines of a class
  - which leaves the level values, the weights and the class sizes alone and
  breaks only the pairing of a line with its own pattern - gives
  kappa(plain) = -0.003 +- 0.006 over 20 permutations, against the 0.929
  measured. Nothing in the machinery manufactures a kappa.
* *The old weights were not what produced the old answer.* The same fit on the
  Set 1 `unc_own` weights gives kappa(plain) = 0.951 +- 0.012, against
  0.929 +- 0.014 on the honest ones. Section 2.6 was right to rule them out
  and right that they bias the answer, but the bias is 0.02, not the 0.9 that
  would have been needed to explain the result away.
* *The anchor's own scale barely matters.* Taking the A constants as
  calculated instead of scaled by 0.9582 moves kappa(plain) to 0.932 +- 0.013.

**kappa(c) is the one value that must not be quoted with its formal error.**
The same permutation control gives +0.026 +- 0.115 for the c class - a scatter
twice the formal 0.051, because 50 informative lines of one class cannot
average anything down. The honest statement is kappa(c) = 0.20 +- 0.12: the
complex lines are measured much nearer the center of gravity than any other
class, and that is all these data say about them.

**And a single kappa per class does not survive section 3.3's first test**;
the plain class needs one value per era. See section 3.3.

### 2.3 Theoretical A values from the level compositions

`levels_pub.xlsx` gives the eigenvector composition of each calculated level.
For a level written in J1-j coupling as `4f2.(3H<4>).5p6.6p (4,1/2)*` — parent
term 3H, parent J1 = 4, outer electron j = 1/2 — the A constant of each
component is

```
A = a_core * theta(parent) * P(J1, j, J)  +  a_outer * P(j, J1, J)
P(X, Y, J) = [ J(J+1) + X(X+1) - Y(Y+1) ] / [ 2 J (J+1) ]
```

with `P` the Landé projection factor and `theta(parent)` a per-parent-term
coefficient calibrated empirically against Reader & Sugar's published pure-θ₁
values (the two J rows of every parent agree to better than 1%).

**Validation against Reader and Sugar's calculation.** The method reproduces
all 22 of Reader and Sugar's 4f²6s A constants with **rms 0.0048 cm^-1** over
a range −0.106 to +0.205 cm^-1. Freeing the radial constants returns
`a_4f = 0.0282 ± 0.0011` and `a_6s = 0.6357 ± 0.0072` against the published
0.029 ± 0.003 and 0.639 ± 0.002 — both inside 1 sigma.

**What that check proves, and what it does not (corrected 2026-09-21).** Those
22 constants are **not measurements**, and this document called them measured
until this correction. Reader and Sugar [R1] publish the angular coefficients
θ₁ and θ₂ — their Table III, headed *"Calculated values of θ₁ and θ₂
appearing in the expression A = θ₁a_4f + θ₂a_6s"* — together with the radial
constants a_4f and a_6s least-squares fitted to 23 line-width equations. They
publish no per-level A value anywhere. Every row marked
`Reader & Sugar 1965 (semiempirical)` in `A_hfs_levels.csv` is θ₁a_4f + θ₂a_6s
evaluated from that table by this project; for `059003.000190` (J = 1/2),
−0.1381×0.639 + 0.33×0.029 = −0.0787, the tabulated value exactly. Their
measurements also excluded low J outright: *"In order to avoid complex
patterns, levels with J < 7/2 were not used."*

The rms of 0.0048 cm^-1 therefore measures how well this project's angular
algebra agrees with Reader and Sugar's angular algebra. That is a real and
necessary check — it proves the θ computation here is right, and the returned
radial constants confirm it independently — but it is **not** a test against
experiment, and 0.0048 cm^-1 must not be quoted as the uncertainty of a
computed A constant against the true one.

**The measured calibration set.** The only measured A constants of these
levels in existence are this project's own, from resolved patterns:

| level | measured A (cm^-1) | semiempirical A (cm^-1) | source of the measurement |
|---|---|---|---|
| 1093, 4f²6s `_3P4P`, J = 1/2 | +0.505 ± 0.007 | +0.500 ± 0.005 | `4f26s_4P_0.5_hfs.md` sec. 2 |
| 1089, 4f²6s `_3P2P`, J = 1/2 | −0.157 ± 0.010 | −0.0787 ± 0.0022 | `4f26s_2P_0.5_hfs.md` sec. 1 |
| 964 (upper partner) | +0.084 ± 0.009 | +0.087 | `4f26s_4P_0.5_hfs.md` sec. 2 |
| 972 (upper partner, sign undetermined) | 0.132 ± 0.012 | +0.030 ± 0.012 | `4f26s_4P_0.5_hfs.md` sec. 2 |

Two of the four agree within about 1 %, and two disagree by a factor of two
and by 8 sigma. The pattern in the failures is physical, not random: 1093's A
is dominated by a single large term, while 1089's is the small residue of a
cancellation between −a_6s/3 = −0.213 and a 4f contribution near +0.13, so a
few per cent of composition error moves it by a factor of two. **The honest
statement is that a computed A constant is reliable to about 1 % where one
term dominates and to no better than a factor of two where the value is a
residue of cancellation**, and every level whose computed A is small relative
to its individual terms must be flagged as such wherever these constants are
used or published. Enlarging this measured set is one of the returns expected
from Step 2.

`levels_pub.xlsx` also independently confirms all 22 parent assignments that
had been made by the energy-ladder argument, including the 3P2 pair that was
the source of an earlier transcription error.

### 2.4 The flags as inequalities

Of 162 flag constraints whose partner level has a computable composition, 154
are satisfied everywhere in the plausible range of the 6p radial constants,
5 give one-sided lower bounds only, and 3 are unsatisfiable. **The flags
cannot measure a_6p** — Sugar wrote a flag only when the pattern *looked like*
a flag pattern, which presupposes small partner hfs: a selection effect. The
6p constants were therefore adopted from the 6s-to-6p ratio trend in
comparable systems (Cs I 0.127/0.022, Ba II 0.185/0.028). Across the whole
plausible box this choice never changes a predicted flag direction and moves
any line's correction by less than 0.09 cm^-1.

Flag-direction agreement between prediction and Sugar's own flags is
**261 of 264**.

### 2.5 Current results

**The fit quality and level shifts below inherit the provisional kappa values
of section 2.2 and are provisional with them.** The A constants and their
provenance do not depend on kappa and stand on their own.

* `A_hfs_levels.csv` — 129 levels: 22 evaluated from Reader and Sugar's
  published θ coefficients [R1] — semiempirical (calculated θ, fitted radial parameters), not measured (section 2.3) —
  57 computed from
  composition, 7 bounded by flags alone, 41 undetermined, 2 in conflict.
* One genuinely new constant measured by the flags alone:
  **A(059003.000233, 4f5d², J = 5/2) = −0.0462 ± 0.0040 cm^-1** from 3 lines.
* Fit quality, standalone weighted fit, 3834 lines:

```
                  chi2/dof   rms residual   |res| > 3 sigma
  Set 1             4.765       0.2358 cm-1       227
  Set 2            3.415       0.2303 cm-1       174
     c lines                    0.383 -> 0.191
     flagged                    0.143 -> 0.137
     plain                      0.237 -> 0.237
```

* Level shifts between the two sets: rms 0.177 cm^-1 over 535 levels, 43
  levels beyond 0.3 cm^-1, largest −1.230 cm^-1 (`059003.000189`). For the 22
  known levels the shift equals `−I·A·J` exactly, as it must.

---

### 2.6 Why section 2.2 has to be done again

Two defects, both fatal to the numbers as they stand.

**The uncertainties used were the wrong ones.** The weights in that fit came
from `Pr3_lines.xlsx`, whose `unc_own` column was built *for Set 1*: it
inflates the uncertainty of a line that deviates strongly from its Ritz
value, precisely so that a fit using measured (hfs-displaced) wavenumbers
stays statistically honest. A trend fitted through those weights is fitted
through a function of the very residuals it is trying to explain. It will
find structure whether or not any is there.

**The trend derivation must start from Sugar's own stated uncertainties**,
with only the few modifications the data force:

| source | what Sugar states | applies to |
|---|---|---|
| Sugar 1969 [R4] | "From the variation of wavelengths obtained with different exposures the accuracy of the final values is estimated to be +-0.004 angstrom" | believed to be the plain (blank-character) lines only |
| Sugar 1974 [R5] | plain-line uncertainty "~0.003 angstrom" | plain lines |
| Sugar 1974 [R5] | "average deviation of +-0.007 angstrom" from Ritz values | blended and complex lines, for which no uncertainty is stated |

The procedure:

1. Weight every line by the value Sugar states for its character and era.
2. Fit the trend.
3. Mark the lines that deviate clearly from the established trend as
   **probably misassigned** and exclude them from the trend derivation. They
   are not thereby deleted from the line list: they become candidates for
   review by hand, and each decision goes in `line_decisions.csv` with its
   reason, per the project's normal ledger discipline.
4. Adopt the **rms of the filtered residuals, per character per era**, as
   Sugar's uncertainty for that class. That number goes into the publication,
   and a reader can reproduce it from the published table - which is the
   standing requirement on every uncertainty in this work.

The inflated `unc_own` values are Set 1 machinery and play no part in Set 2.

**Done, 2026-09-21, in `hfs_kappa.py`; the table is `hfs_uncertainties.csv`.**
The procedure above is carried out as written, with one addition the data
forced. The rms of the residuals of a *fitted* quantity understates its
uncertainty, because the fitted level values have already followed the
residuals; a class of 29 lines came out at 0.0023 cm^-1, which is four times
better than Sugar's best plain line and is an artifact of that. What the class
really has left is the sum of `1 - h` over its lines, where `h` is the line's
leverage in the fit - the share of its own residual the level values have
already absorbed - and the adopted value is the rms taken over that sum. No
class is allowed below the project's precision floor of 0.0055 cm^-1, and a
character with fewer than 20 lines in an era is pooled into an `other` class
of that era rather than measuring itself.

| char | era | lines | Sugar's stated value | adopted |
|---|---|---|---|---|
| (plain) | 1974 | 1614 | 0.0360 | **0.0459** |
| `w` | 1974 | 420 | 0.0767 | 0.0485 |
| `*r` | 1974 | 133 | 0.0668 | 0.0858 |
| `*v` | 1974 | 124 | 0.0704 | 0.0725 |
| `c` | 1974 | 63 | 0.0734 | 0.1866 |
| `d` | 1974 | 41 | 0.0460 | 0.0581 |
| `ch` | 1974 | 29 | 0.0112 | 0.0186 |
| other | 1974 | 18 | 0.0981 | 0.0234 |
| (plain) | 1969 | 1755 | 0.2457 | **0.2454** |
| `cl` | 1969 | 76 | 0.5004 | 0.2751 |
| `h` | 1969 | 49 | 0.3639 | 0.4798 |
| other | 1969 | 22 | 0.3474 | 0.3789 |
| `bl` | 1969 | 21 | 0.2960 | 0.6371 |

All values are cm^-1, and the stated column is the mean over the lines of that
class of Sugar's own figure converted line by line through
`u_wn = u_lambda * wn^2`. The 1969 values are larger than the 1974 ones by a
factor of five not because the plates were worse but because the same error in
angstrom is worth `wn^2` more in cm^-1 at short wavelength; a reader can
reproduce every number in the table from the published line list and these two
stated figures.

**What the table says.** Sugar's stated uncertainties are close to right where
he stated them: a plain 1974 line comes out at 0.046 cm^-1 against his 0.036,
and a plain 1969 line at 0.2454 against 0.2457. The characters he said nothing
about spread either way - `w` and `cl` are *better* than the 0.007 angstrom
average deviation he quotes, `c` and `bl` two to three times worse.

**Ten lines are set aside as probably misassigned**, listed in
`hfs_misassigned_candidates.csv` with their residual in units of their own
adopted uncertainty. They are candidates for review in IDEN2, not deletions:
no row was written to `line_decisions.csv`, and each still needs the plate
evidence looked at by hand. Four of the ten are flagged lines, which is where
a misread flag or a wrong A constant would show first - 37136.190 `*v` at
-13.9 sigma is the largest residual in the whole line list - and two involve
`059003.000233`, the level whose A the flags alone measure.

**The kappa model may not survive this.** If the unflagged classes come out
at kappa = 0 once honest weights are used, the center-of-gravity hypothesis
of section 1.1 is confirmed, the correction applies to the flagged lines
only, and the level shifts of section 2.5 shrink by roughly an order of
magnitude. That would be a *better* outcome than the present one: fewer
assumptions, a smaller correction, and one that rests on Sugar's own words
rather than on a fitted convention factor.

## 3. Open scientific problems

### 3.1 The 41 undetermined levels — the main limitation

**Largely resolved 2026-10-04 (Step 4):** Sugar's resolved components, fitted
with the confirmed semiempirical constants as priors, give A to 35 of the 39
levels; 4 remain "not determined" (000258, 000295, 000453, and 000448 held
back). What follows is the state before that.

*Rewritten 2026-10-03 with the user.* The levels in `A_hfs_levels.csv` with
the source "not determined" have no A constant. There were 41 when this
section was first written (37 4f5d², 4 4f³). Hand measurements have since
given values to two of them, so 39 remain: 35 4f5d², 2 4f³, one 4f5d6s and
one 4f²6p. A line touching such a level gets no displacement. The unknown
A, taken as 0.05 cm^-1, times I·J, is added to the line's uncertainty instead
(decision D10). Here I = 5/2 is the nuclear spin of ¹⁴¹Pr and J the level's
angular momentum.

**The user's rule, 2026-10-03: a line is moved only by an A constant that
rests on measurement, or on a calculation the measurements have confirmed.**
Until such a value exists, the lines of these levels stay at their observed
positions, in Set 3 as in the head frame, with the allowance above.

**Why a calculation alone is not enough for these levels.**

* *a_5d is unknown.* a_5d is the radial hyperfine constant of a 5d
  electron, the 5d counterpart of the a_4f and a_6s of section 2.3. Nobody
  has measured or published it for Pr III. Any value would be scaled from
  other spectra or from the radial integrals of the Cowan calculation.
* *The 4f5d² levels are strongly mixed.* Their A would be a sum of
  contributions of both signs. That is the cancellation case in which
  section 2.3 found the calculation wrong by a factor of 2 (IDEN2 row 1089).
* *There is almost nothing to test it against.* Only one of them,
  `059003.000233`, has an A measured independently: −0.0462 ± 0.0040 cm^-1,
  from its flags.

The coupling scheme of the labels is not the obstacle. `levels_pub.xlsx`
writes the levels in mixed schemes: LS, jj, JK, and some with different
components in different schemes. The full eigenvector composition in any one
scheme gives A. LS is the natural scheme for 4f5d², where the 4f and 5d
one-electron contributions are projected onto J through the intermediate LS
quantum numbers. The user holds the complete composition tables in every
scheme. The LS tables go into the project when the check calculation of
Step 4 is made, and not before.

**The semiempirical 4f²6s and 4f²6p constants are a different case: they are
confirmed by measurement.** Sugar measured the positions of the resolved
hyperfine components of his flagged lines (section 3.5). A component
*spacing* is the distance between two neighboring components of one line's
pattern. The semiempirical A constants of the line's two levels predict every
spacing, with no free parameter except Step 2's scale. The test of
2026-10-03 used 135 patterns (271 spacings) whose two levels both have a
semiempirical A, after scaling by 0.9638:

```
rms of the measured spacings themselves        0.795 cm-1   (what A = 0 predicts)
rms of measured minus predicted                0.065 cm-1
rms expected from Sugar's 0.003 A alone        0.041 cm-1
share of the measured variation explained      99.3 %
first spacing with the right sign              135 of 135
measured/predicted, first spacing              median 0.991; 124 of 135 within 20 %, all within 50 %
```

Those 135 patterns test 74 of the 81 levels with a semiempirical A, and 41 of them in
three or more patterns. The calculation is worse than Sugar's measurement
by a factor of 1.6, and better than no correction by a factor of 12. One
level is contradicted: `059003.000127`. It is in 4 of the 10 worst spacings,
all 25 to 35 per cent narrower than predicted, and it was among the levels
beyond 3 sigma in Step 6. Seven levels with a semiempirical A, all 4f²6p, are tested by
no pattern: `059003.000232`, `000249`, `000273`, `000290`, `000299`,
`000301` and `000444`. Two of them have sizeable I·A·J: 0.35 cm^-1 for
000299 and 0.51 cm^-1 for 000301. Step 4 deals with both cases.

**The routes to an A from the observed data,** for Step 4:

1. *Flags.* Sugar wrote `*r` or `*v` only where he saw the pattern shaded,
   so a flag bounds the A of the line's levels (section 2.4). This is done:
   it gave 7 levels an interval and 000233 a value.
2. *Resolved components.* The global fit of section 3.5
   (`hfs_A_from_patterns.csv`) already gives a value for 35 of the 39
   levels: 32 4f5d², both 4f³ and the 4f²6p one. 18 of them rest on a
   single pattern. That fit freed every A at once, so it left a common
   offset undetermined (+0.03 cm^-1 against the calculation). Holding the
   confirmed semiempirical constants fixed at their scaled values removes the
   offset, and with it the reason the values were not adopted.
3. *Lines of different classes.* Each class of line sits a different
   fraction κ of the way from the pattern's center of gravity to its
   strongest component: 1 for flags, 0.93 for plain 1974 lines, 0.71 for
   plain 1969 lines, 0.27 for `c` lines. With the level's own A set to 0,
   its lines of class κ are left off the head-frame Ritz value by
   (1−κ)·I·A·J. Comparing its classes therefore measures I·A·J. This is a
   weak route: it needs a level with lines in two classes and an I·A·J of
   the order of the 1969 uncertainty (0.25 cm^-1) or more.
4. *Line widths.* The per-level widths of `level_hfs_widths.csv` measure
   the size of A but not its sign.


### 3.2 The two conflicting levels

`059003.000223` (theory +0.0227 ± 0.0054, flags demand A < −0.017) and
`059003.000248` (theory +0.0251 ± 0.0055, flags demand simultaneously
A > +0.125 and A < −0.031 — no value can satisfy both). Closing the
`000223` gap would need `a_6p(3/2)` above 2 cm^-1, about 100 times the
plausible value, so the theory is not the suspect.

**Work required:** treat these as defective *lines*, not defective levels.
Re-examine each flagged line touching them, with the plate evidence in IDEN2.
The leading candidate is the line at 39917.569 cm^-1 (residual +1.010 cm^-1).
Decide line by line whether the flag is a misreading, the assignment is wrong,
or the line is a blend. Record each decision in `line_decisions.csv` with its
reason, per the project's normal ledger discipline. Do **not** silently drop
constraints to make the fit agree.

**Withdrawn 2026-09-29.** Neither conflict survives a line-by-line check, and
both levels keep their semiempirical A. The label "CONFLICTS with own flags" was
removed from `A_hfs_levels.csv`.

* `000223`: the bound A < −0.017 came from one line alone, 34926.912 `*v`
  (to `000134`, A = −0.0206), which the user rejected on 9/12. Its four
  remaining flags (35188.811 `*v`, 33162.180 `*v`, 32455.968 `*r`,
  30631.472 `*v`) all go to 4f2.6s levels whose A is three to four times
  larger, and they bound A only to −0.022 < A < +0.078. The calculated
  +0.0227 lies inside that range.
* `000248`: the bound A > +0.125 came from 39917.569 `*r` (to `000127`),
  also rejected on 9/12. The bound A < −0.031 came from 33417.350 `*v` (to
  `000149`). That line is a blend: its stronger component is
  `000133`–`000213`, with a calculated intensity six times larger, and the
  semiempirical A predicts `*v` for that component, so the flag belongs to it.
  The three other flags (33510.631, 33000.490, 30353.320) agree with the
  semiempirical +0.0251. On 9/30 the user rejected the `000149`–`000248`
  component of 33417.350 outright ("worsens cg; opposite shift to hfs
  calculation"), so the line no longer constrains `000248` at all.

`hfs_line_corrections.csv` still gives these two levels S = 0, because it
was written while the labels stood. It is left as it is until Step 5 rewrites
it (the user's decision of 2026-09-29), and `tests/test_hfs_components.py`
counts its seven affected flagged lines separately. Until then
`hfs_kappa.py`, which reads the labels, uses the semiempirical A for both levels
while the correction file gives them S = 0.

### 3.3 Whether the κ model is the right physical model at all

The model assumes a single κ per class, constant across the spectrum. Two
things could break it:

* **Resolution varies with wavelength.** A pattern that is resolved at 4000 Å
  may be unresolved at 2000 Å for the same A values, so the effective κ of a
  plain line should drift. **Test:** split the plain class in two by
  wavenumber at the 1974/1969 era boundary (47500 cm^-1) and fit κ separately.
  If the two agree within their uncertainties, the single-κ model stands.
* **Pattern intensity structure.** The extreme component is the strongest only
  for a well-separated pattern of a transition with ΔJ = ΔF; for other
  branches the strongest component is an interior one. **Test:** fit κ
  separately for ΔJ = +1, 0, −1 within the flag class. A significant spread
  means κ should be replaced by a per-line geometric factor computed from the
  relative intensities of the pattern components, which is a standard
  Racah-algebra calculation and is entirely doable.

Both tests are cheap — each is a re-run of the existing joint linear solve
with one extra column — and both should be done **before** anything is adopted
into the pipeline.

**Both were run on 2026-09-21 (Step 3, `hfs_kappa.py`). The first fails and
the second passes.**

**Test 1 fails: the plain class needs one kappa per era.**

```
kappa(plain, 1974)  +0.944 +- 0.014     below 47500 cm^-1
kappa(plain, 1969)  +0.633 +- 0.035     above it
difference          +0.311 +- 0.034     9.2 sigma
```

The plan's own acceptance test was agreement within 2 sigma, so the single-κ
model is refused and **Step 5 must carry one kappa per era for the unflagged
lines.** Three things were done before accepting that.

* *Is it a drift or a step?* Section 3.3 predicts a smooth drift with
  wavelength. It is not one: in six bins holding equal numbers of informative
  plain lines, κ is flat at 0.92-0.99 from 14800 to 52800 cm^-1 and only the
  last bin falls. Splitting at the era boundary instead gives 0.945 +- 0.014
  below it and 0.652, 0.585 and 0.657 in three bins above it - a step at the
  boundary, and flat on both sides of it. The boundary is 2105 angstrom, which
  is where Sugar's own Pr IV statement (section 1.1) puts the change in his
  measuring practice, and the two campaigns are also different plates and a
  different instrument.
* *Is it the A constants rather than the measurement?* The two eras see
  different levels, and a systematic error in the A constants of the levels
  peculiar to one of them would look the same. It is not that: confined to the
  74 levels that carry a displacement in **both** eras, the split is 0.938 +-
  0.015 against 0.653 +- 0.037, unchanged.
* *Is it real at all?* Permuting the displacements gives -0.003 +- 0.006.

**Test 2 passes: no per-line geometric factor is needed.** Splitting the flag
class by branch, with dJ = +1 held at 1, gives dJ = -1: 1.006 +- 0.019 (67
lines) and dJ = 0: 0.957 +- 0.019 (94 lines) - a spread of 0.05 over the three
branches. Section 3.5 had already settled the physics behind it: the extreme
component is the strongest one for every J pair the line list contains, the
single exception in the whole table being 1/2 -> 1/2, which no flagged line
is. The Racah-algebra per-line factor this test was to trigger is therefore
not needed, and κ stays a number per class.

### 3.4 Lines whose correction is estimated rather than computed

7 levels are constrained by flags to an interval only; 41 are undetermined.
Every line touching one of them carries an inflated correction uncertainty. As
3.1 progresses this shrinks. The final paper must report the exact count of
lines whose correction is computed, estimated, and absent.

---

### 3.5 The 501 measured hfs component wavelengths - the highest priority

`Pr3_lines.xlsx` was built for Set 1 and carries one row per line. Sugar's
1974 Table 1 carries more than that: for the flagged lines it also prints the
wavelengths of the resolved hfs components, which the user identifies as the
strongest components of the flag patterns. They survive in the checked
extraction `Pr_3_Sugar74_Table1_extracted_v3_gemini-3-flash-preview.xlsm`
(sheet `Table 1`) as the rows carrying `Intens = 0` in column 3; the
wavenumber to use is column 31, `wn adopted`.

What is in there, counted directly from the workbook on 2026-09-13:

* **501 such rows**, out of 4334 rows in the sheet.
* Each pairs with a flagged line: 500 of the 501 lie within 3 cm^-1 of an
  `*r`/`*v` line. The one that does not must be looked at by hand.
* **The side is right in 500 cases out of 500.** Every companion of an `*r`
  line lies at lower wavenumber (redward of the tabulated line); every
  companion of a `*v` line lies at higher wavenumber (violetward). Not one
  exception. Both the extraction and the reading of the flags are sound.
* The displacements from the tabulated line run 0.10 to 1.38 cm^-1 for `*r`
  (193 rows, mean -0.56) and 0.14 to 2.20 cm^-1 for `*v` (307 rows, mean
  +0.82) - the scale the hfs model predicts for D.
* 248 of the 296 flagged lines carry at least one companion: 89 carry one,
  96 two, 41 three, 14 four, 8 five. A line with several companions has
  several resolved components of one pattern listed, and - this is the part
  that needs explaining - all of them on the flagged side of the tabulated
  wavelength.

**Why this outranks everything else in this document.** The kappa model of
section 2.2 *infers* the displacement D of a pattern from residuals. These
rows *measure* it, for about 250 lines, at Sugar's own measurement precision.
They can settle, with no fitted convention factor anywhere:

1. **Which feature the tabulated wavelength is.** If the tabulated line were
   the center of gravity of the pattern, its components would straddle it.
   They do not: every listed component lies on one side. Taken at face value
   that says the tabulated wavelength of a flagged line is one end of the
   pattern, which is in tension with the center-of-gravity hypothesis of
   section 1.1 and must be reconciled before anything else is done. The
   possibilities to weigh are that Sugar tabulated the extreme component for
   flagged lines only (consistent with both, since section 1.1 exempts the
   flagged lines), or that the components listed are a partial list - only
   those he could resolve - which for a pattern whose strong end is resolved
   and whose weak end is not would look exactly like this.
2. **A constants, measured directly.** In the resolved limit the spacings of
   consecutive components of one pattern are the hfs intervals of the two
   levels; a pattern with three or more listed components over-determines the
   two A constants of the transition. These are measurements, independent of
   the `levels_pub.xlsx` compositions of section 2.3 and of the
   flags-as-inequalities argument of section 2.4.
3. **The kappa(flag) = 1 assumption**, currently the single external
   assumption Set 2 rests on (section 2.1), becomes testable, and if it holds
   it becomes measured rather than assumed.

**What was done, and what it settled (2026-09-21).** All four items below are
carried out in `hfs_components.py` (the extraction) and `hfs_patterns.py` (the
physics), with `tests/test_hfs_components.py` holding every number quoted here
as a regression test.

1. **`hfs_components.csv`**, 500 rows: `wn_line, char, wn_component, dwn,
   n_components, low_id, upp_id` and seven audit columns, joined to
   `line_classifications.csv` by nearest wavenumber within 0.006 cm^-1. The
   one unpaired row is resolved: `Intens = 0` is overloaded, and three of the
   501 rows were zeroed by the curation of the workbook rather than printed
   without an intensity by Sugar. Two of those three (sheet rows 80 and 2389)
   are lines whose classification duplicated the neighboring flagged line's
   and sit 0.1-0.3 cm^-1 on its flagged side - resolved components in all but
   name - and are kept, marked `origin = curated_duplicate`. The third (sheet
   row 2053) is an error in Sugar's table, 6.9 cm^-1 from any flagged line,
   and is dropped. All 296 flagged lines find their row in the pipeline's
   line list; 248 of them carry components, 89 with one, 96 two, 41 three, 14
   four and 8 five.
2. **The interval rule holds with no free parameter.** Taking the tabulated
   line to be the extreme component F = I+J_low -> I+J_upp and the k-th
   listed companion to be the k-th rung below it, the position of every
   companion is
   `d(k) = (k/2)[(2F1max+1-k) A_low - (2F2max+1-k) A_upp]`,
   which is linear in the two A constants and needs no matching of
   intensities. Over the 313 measured positions of the 155 patterns whose two
   levels both have a semiempirical A, that predicts the measurement with an
   **rms of 0.090 cm^-1** against displacements running to 2.2 cm^-1 and a
   measurement precision near 0.03 cm^-1. One free parameter - a single scale
   factor on every A constant in the project - brings it to 0.084 cm^-1 at
   **s = 0.958 +- 0.006**: the measured intervals are 4.2 per cent smaller
   than the semiempirical ones. That 4.2 per cent is a measurement about the
   radial parameters of section 2.3 and belongs in the error budget of every
   computed A.
3. **What the tabulated wavelength is: the strongest component of the
   pattern.** It is not the center of gravity - not one of the 500 components
   straddles it. It is not displaced from the strongest component either: the
   competing one-parameter model, a constant shift towards the degraded side,
   fits worse (0.085 cm^-1 at -0.030 cm^-1), and with both parameters free
   the shift collapses to +0.003 +- 0.005 cm^-1 while the scale factor stands.
   Decision D7 is therefore answered from the data and is no longer a
   decision.
4. **kappa(flag) = 1 is confirmed as exact.** The strongest component lies at
   `D = I (A_upp J_upp - A_low J_low)` from the center of gravity, because
   `K(I+J) = 2 I J`; D is precisely the quantity of section 0, and
   `delta_cm-1` of `hfs_line_corrections.csv` is already identically that,
   with kappa = 1, on 148 of the 155 patterns. The seven exceptions are the
   lines of the two conflicting levels of section 3.2, for which the
   correction file deliberately writes S = 0. So the correction the pipeline
   is built to apply is the one the measurements confirm; what changes is its
   scale, by the 4.2 per cent of item 2, and its uncertainty.

**What the same data say about the level constants.** Solving for every
level's A constant at once from the 435 measured positions (216 patterns, 121
levels) reproduces the semiempirical constants of 71 levels with a scatter of
**0.020 cm^-1** about a common offset of +0.022 cm^-1, and gives a first A
constant for 38 levels that had none. The common offset is the fit's own
near-singular direction - adding the same constant to every A changes the
predicted spacings only through `J_low - J_upp`, and not at all when they are
equal - so the *differences* of these constants are measured and their
absolute scale is not. They are therefore reported in
`hfs_A_from_patterns.csv` and are **not** written back into
`A_hfs_levels.csv`; that is Step 4's business, with the LS-coupled calculation
in hand. Two results worth carrying there: level `059003.000190`, the 4f2.6s
2P_1/2 whose semiempirical A is the residue of a cancellation (section 2.3),
comes out at -0.117 +- 0.005 against the semiempirical -0.0787 and the directly
measured -0.157 +- 0.010, again too small in the calculation; and the two
conflicting levels of section 3.2, `059003.000223` and `059003.000248`, both
come out larger than semiempirical (+0.043 +- 0.011 against +0.0227,
+0.057 +- 0.016 against +0.0251).

**Refreshed 2026-09-29.** The numbers above are those of 2026-09-21. Since
then the workbook's curation gave row 2053 back its intensity (so 500 rows
carry `Intens = 0` and none is dropped), and the classification changed on
five flagged lines: three of the ten Step 1 candidates (25042.950, 31391.380,
31558.079) are rejected, 11956.920 is unclassified and 11777.630 is no longer
a blend. On the regenerated `hfs_components.csv`: 309 positions of 152
patterns, rms 0.087 cm^-1 with no free parameter and 0.080 at
s = 0.958 +- 0.006; the offset model again loses (0.081 at -0.031), and with
both free the offset is -0.003 +- 0.005. kappa = 1 holds on 145 of the 152,
the same seven conflicted lines excepted. The global fit (430 positions, 213
patterns, 120 levels) drops from chi2/dof 2.92 to 1.47, rms 0.031 cm^-1: the
rejected lines were the worst-fitting ones. Scatter against the calculated
constants 0.023 cm^-1 about +0.030; 37 first A constants. Level 190 is
unchanged at -0.117 +- 0.005; 223 and 248 come out at +0.050 +- 0.011 and
+0.065 +- 0.016. `tests/test_hfs_components.py` holds the new values.

**Remeasured 2026-10-03, on semiempirical constants only.** The scale of item 2
had been fitted on every level with a determined A. That included 7 levels
whose A is a flag interval, which is a measurement and not a calculation.
Fitted on the semiempirical constants alone ("Reader & Sugar 1965 (semiempirical)"
and the sources beginning "composition"), on the current iter_hfs
assignments: 135 patterns, 271 positions, rms 0.072 -> 0.065 cm^-1 at
**s = 0.9638 +- 0.0079**. The uncertainty is the jackknife one, from
leaving out one pattern at a time; the formal one is 0.0048. The old method
on the current inputs still gives 0.958, so the assignments changed nothing.
The new value lies within 1 sigma of the old one.

The scale is now applied to semiempirical constants only, and in the pipeline
and the fits of kappa alike: `hfs_kappa.scaled_A`, read by
`hfs_correction.Model` and by `hfs_kappa.read_A_constants`.
`wavelength_calibration.py` no longer multiplies D by a scale of its own.
The constants from flags, resolved components and hand estimates are used
as written.

Two results are recorded and not acted on. Fitted on the `*r` lines alone
the scale is 0.981 +- 0.009, on the `*v` lines alone 0.958 +- 0.006, which
is 2.2 sigma apart. With scale and offset both free, the fit gives
0.981 +- 0.011 and an offset of −0.016 +- 0.009 cm^-1 toward the shaded
side, 1.7 sigma; the first measurement had found −0.003 +- 0.005. The
one-parameter model is kept. Section 3.1 gives the same test pattern by
pattern.

**Left for later.** Three patterns list more components than the principal
ladder has rungs - 29514.531 (J 3/2 -> 1/2), 29584.911 and 29602.342 (both
J 5/2 -> 3/2). At such low J the satellite components carry 6-9 per cent each
and are resolvable, so these three need the full component list rather than
the ladder and are set aside by `read_patterns`.

### 3.6 Every external datum must be referenced

Every external number used in the analysis - a radial parameter, a quoted
uncertainty, an hfs constant borrowed from an analogous system, a nuclear
moment - must be cited to a source the user can check, and the citation must
be in section 9 of this document so that it can be carried into the paper
without being reconstructed.

Two items now in use are **not** yet traceable to a verified source and must
be settled before publication: the Cs I and Ba II 6s/6p hfs constants used in
section 2.4 to adopt `a_6p`, and the exact bibliographic data of the two
Sugar line lists (section 9).

### 3.7 Which inflated-uncertainty lines hfs explains

`inflated_unc_lines.txt` holds 157 lines whose uncertainty the user raised by
hand (140 until 2026-09-29; 17 added on 9/30, see below). Many are firmly established identifications made by Sugar, and no
single cause explains them all. The user's guess is unrecognized blends with
impurity lines or with Pr II and Pr IV lines, plus some measurement or
typesetting errors that went undetected. The user's notes there ("hfs of
f26s", "hfs of fds") were written from IDEN2, where only one end of a
transition is in view at a time, so they are incomplete and not always right.
**Only the lines that can be traced to the hfs of specific levels are to be
marked and used**, for instance in the Step 5 correction or as tests of a
revised A. The rest stay unexplained.

`inflated_unc_hfs.csv` (2026-09-29, iter classification set) is the condensed
result of checking every one of the 140 lines. It gives one row per line:
the classification with the strongest hfs case (the line's number of
classifications is in `n_class`), `D`, the calculated pattern width, the iter
LOPT residual, and an attribution:

* `hfs` (36 lines): |D| >= 0.3 cm^-1, or a calculated pattern at least
  1.0 cm^-1 wide, or the line involves one of the two suspect levels below.
* `hfs possible` (33): one level lacks an A, and the other has
  |I·J·A| >= 0.3 cm^-1.
* `no hfs evidence` (71): neither test is met in any classification of the
  line.

For 4f5d6s levels, `A_estimated = yes` marks an A estimated from the 6s
contact term alone in pure LS coupling (a_6s = 0.639 cm^-1), good to perhaps
±50 %. Level `000642` (IDEN2 1093) is given its measured A = +0.505.

Two levels stand out:
* `000150` (IDEN2 1129, 4f2.6s (3F)2F, J = 7/2) is under 13 inflated lines,
  ten of them tagged "hfs" by the user. Its semiempirical A = −0.0329 is the
  residue of a cancellation, as for `000190`, where the calculation is off by
  a factor of two. It is the first candidate for a measured A.
* `000448` (IDEN2 938, 4f5d6s (1G)2G) is under four inflated lines. The two
  unblended ones have O−C of +0.58 and +0.71 cm^-1.

A large D alone leaves no LOPT residual: a line read at its strongest
component carries the shift into the level energies. Residuals come from
lines read elsewhere, such as `c` lines, the 1969 lines, and very wide
patterns. The script that wrote the file is not in the repository. It must
be rebuilt from `A_hfs_levels.csv` when that file changes.

**Additions of 2026-09-30.** The user compared the transitions accepted in the
baseline with those accepted in `iter/` and made sure that `iter/` keeps
every one IDEN2 and LOPT support. The check found missed classifications, and
lines that sit at the center of gravity of the unresolved hfs of a 4f2.6s
level. Seventeen lines were given inflated uncertainties and accepted, each
with a ledger row reading `IDEN2/LOPT; unc inflated`:
* 14 tagged "hfs of f26s" (0.3-0.5 cm^-1): 59829.1143 and 30126.3092
  (`000141`); 33722.4702, 28508.8216, 44134.0401, 27443.3325, 33934.5817
  (`000190`, IDEN2 1089, whose A is about twice the semiempirical value);
  51071.138, 31478.9216 (`000140`); 33044.2299, 36647.5418, 38425.7687
  (`000183`); 61680.8807 and 36914.5611, accepted by the classification
  without a ledger row. For 33044.2299 and 36647.5418 the competing
  candidates `000189`–`000297` and `000139`–`000252` were rejected.
* 3 tagged "outlier in IDEN2/LOPT, OK with inflated unc" (0.12-0.2 cm^-1):
  34406.941, 31664.7593, and 38052.9185, whose reject of 9/5 was replaced by
  an accept.
30126.3092 (`000141`–`000210`), rejected on 9/11 as "looks wrong", is
accepted again. 22561.9597 (`000388`–`000617`) was accepted with its own
uncertainty. `inflated_unc_hfs.csv` covers only the first 140 lines.

## 4. Open decisions for the user

These are decisions, not tasks. Nothing below is done without a word.

| # | Decision | Recommendation |
|---|---|---|
| D1 | Apply the correction inside `make_LOPT_input.py`, so that `Pr3_lines.xlsx` remains the record of what was *measured* and the correction is a pipeline transformation, and refit with LOPT | **Answered, 2026-09-29: yes**, behind a switch in `lineclass_config.toml` (off in the baseline; on in the set that is corrected), together with the plate calibration of Step 3c. Whether the correction reaches `classify_lines.py` and IDEN2 as well is D11 |
| D2 | Publish both sets, or only Set 2 | Both, with Set 1 the primary line list (measured quantities) and Set 2 the primary level list (physical quantities) |
| D3 | Whether flagged lines receive the full D correction | Yes in principle, but the size of D for a flagged line now comes from the measured components of section 3.5 wherever one exists, not from theory and not from a fitted kappa |
| D4 | Update `CHAR_DLAM[('1974','c')]` from 0.0337 A to the leftover 0.0203 A once the correction is adopted | Yes, but only after the correction is in the pipeline, and as part of the second calibration pass (section 6). The value itself will change once section 2.6 is done |
| D5 | Whether unflagged lines need an extra wavenumber-constant term for unresolved hfs, which re-opens `UNC_FLOOR = 0.0055` | Defer until after sections 2.6 and 3.3. Under the center-of-gravity hypothesis the right form is a symmetric widening of the unflagged lines, not a shift, and its size follows from the pattern widths rather than from a fit to residuals |
| D6 | Raise the `*r`/`*v` extra uncertainty: the flat fit wants about 0.147 cm^-1 per 6s level against the 0.05-0.08 applied now | Moot in Set 2, where the displacement is subtracted rather than covered by an uncertainty; keep as is for Set 1 |
| D7 | What the tabulated wavenumber of a flagged 1974 line is - the pattern center of gravity, or one end of the pattern | **Answered, 2026-09-21.** It is the strongest component, F = I+J_low -> I+J_upp, which is the extreme one; it lies D = I(A_upp J_upp - A_low J_low) from the center of gravity and shows no further offset (+0.003 +- 0.005 cm^-1). See section 3.5 |
| D8 | Whether `hfs_components.csv` (section 3.5) becomes a tracked input of the pipeline or an analysis-only file | A tracked input. It is measured data taken from the source paper and belongs beside `Pr3_lines.xlsx`, not in a scratch directory |
| D9 | Whether the unflagged 1969 lines get Sugar's 0.004 A as stated, or a per-character value refitted as in section 2.6 | **Answered by the data, 2026-09-21.** It makes no difference for the plain 1969 lines: the refit returns 0.2454 cm^-1 against the 0.2457 his 0.004 A gives, so his stated value is adopted and the refit is reported beside it. It matters for the characters he states nothing about, where the refit moves the value by a factor of two either way (section 2.6) |
| D10 | Whether Step 5 goes ahead before Step 4 has given the 41 levels without an A constant (37 4f5d2, 4 4f3) their values | **Answered, 2026-09-29: yes.** A line touching such a level gets no displacement, and the uncertainty of the missing A (about 0.05 cm^-1 per level, times I J) is added to its uncertainty. Step 4 refines it later |
| D11 | Where the correction applies: LOPT's input only, or `classify_lines.py` and IDEN2 as well | **Answered, 2026-09-29: all three**, in a set of its own. LOPT: observed wavenumber - kappa*D. `classify_lines.py`: Ritz + kappa*D per candidate. IDEN2: shifted lines in the set's `dlv.dat`, levels at E_cg (revised the same day: level shifts cannot carry kappa). A multiply classified line gets one shift, kappa * sum BF_i*D_i. **Reopened the same day and replaced by the head frame, agreed 2026-09-29 (Step 5):** levels at E_cg + I*A*J, only assigned lines corrected, by (1-kappa)*sum BF_i*D_i, and the centers of gravity computed at the end; implemented 2026-09-30, behind `[hfs] apply`. The plan had proposed LOPT only; the user's objection was that the three tools have been kept in step so that an outlier shows up in all three. Details, under Step 5 below |
| D12 | Which convention factors kappa | **Answered, 2026-09-29:** the Step 3 values - 1 for `*r`/`*v`, 0.944 for plain 1974, 0.633 for plain 1969, 0.20 for `c` - unchanged, because the plate calibration was fitted jointly with them. The displacements are computed fresh; `hfs_line_corrections.csv` (pooled 0.961, an older classification) is not read. **Amended 2026-10-03:** the kappas follow the calibration they were fitted with; after the second pass (Step 7, on the iter_hfs assignments) they are 0.925(12) plain 1974, 0.708(21) plain 1969, 0.271(44) `c`. **Amended again 2026-10-03:** pass 4, with the semiempirical A constants on the measured scale (D14), gives 0.925(12), 0.712(20), 0.276(44). **Amended 2026-10-04:** pass 5, after 83315.6988 was inflated to 0.5 cm^-1, gives 0.926(12), 0.717(21), 0.276(44); these are in `[hfs.kappa]` |
| D13 | How to avoid counting hfs twice | **Answered, 2026-09-29:** for a line whose displacement is corrected, the level widths of `level_hfs_widths.csv` and the extra uncertainty of the registry rows tagged hfs are both turned off; for a line touching a level without A both are kept. `inflated_unc_lines.txt` itself stays as it is: the exemption is made in code |
| D14 | Which scale of the A constants D is computed with, and where | **Answered, 2026-10-03:** Step 2's scale, remeasured on semiempirical constants only, 0.9638 +- 0.0079 (section 3.5). It is applied to semiempirical constants only, the same way in `hfs_correction.py` (the pipeline) and in `hfs_kappa.py` and `wavelength_calibration.py` (the kappas). It replaces 0.9582 on every D in the kappa fit and 1 in the pipeline. The kappas are refitted in calibration pass 4 |
| D15 | Whether a semiempirical A constant may move a line | **Answered, 2026-10-03:** only where measurements confirm it. The 4f²6s and 4f²6p constants are confirmed by Sugar's component spacings (section 3.1) and are used. 000127, which its patterns contradict, gets its value from its own patterns in Step 4. The undetermined levels get no calculated value: their lines stay at their measured positions until Step 4 gives them an A from the data. **Amended 2026-10-04:** Step 4 did, for 35 of them, from Sugar's resolved components with the confirmed constants as priors; 000127 has +0.0775 from its 13 patterns. Those values move lines; 4 levels remain undetermined |

---

## 5. Implementation plan, in order

Each step names its deliverable and its test. Steps 0-4 change nothing in the
pipeline.

**Step 0 - finish the work already in progress, before any of this starts.**
The hfs analysis rests on the identifications, and the identifications are not
final. Two things close first, and neither belongs to the hfs project:

* the **revision of the questionable levels** - every question mark raised by
  `level_positions.py --audit` accepted or rejected by hand, with its ledger
  row;
* the **assignment of the new levels** - the `insert_new_level.py` work,
  including the level at IDEN2 row 742 already accepted by hand, and the
  candidates that `level_positions.py --unknown` and `unfound_levels.py`
  point to.

Every number in this document is computed from the line-to-level assignment
those two change. Running the hfs analysis before they are done means running
it twice.

**Closed 2026-09-29, by the user.** Every question mark of `--audit` has been
reviewed and every firmly found new level inserted. Reviewed does not mean
cleared: about 50 levels remain questionable (`level_positions.csv` of
2026-09-23 marks 58 with `?`). A few of those are 4f2.5g levels, very pure in
JK coupling, whose position is well supported by the very small calculated
separations within their pairs even though they rest on few, often blended
lines; the rest are questionable in fact. **The lines of the questionable
levels are excluded from every analysis that relies on obs-Ritz residuals** -
the kappa fits, the pattern fits and the Step 6 comparison among them. They
stay in the line list and in the LOPT fit.
Of the 58, the user counts 16 as firm (2026-09-29): the twelve 4f2.5g levels
000585, 000586, 000587, 000591, 000592, 000599, 000605, 000611, 000613,
000616, 000617, 000620, and 000642 (the alternate is its other hfs
component), 000538, 000376 and 000347 (the alternate loses lines). The other
42 are questionable. The list is `questionable_levels.txt`.
*Test:* `check_sync.py` reports no ERROR findings, and `level_positions.py
--audit` leaves no unreviewed question mark.

**Step 1 - rebuild the uncertainty basis for the trend work (section 2.6).
DONE 2026-09-21, except for the ledger rows, which need the user.**
Sugar's stated uncertainties per character per era, the trend fit on those
weights, the outlier filter, and the per-class rms adopted as the final
uncertainty.
*Deliverable:* `hfs_uncertainties.csv`, 13 classes with the count behind each
(the table is reproduced in section 2.6), and
`hfs_misassigned_candidates.csv`, the 10 lines the fit rejects beyond 4 sigma.
Both are written by `hfs_kappa.py`.
*Test:* the adopted uncertainties are reproducible from the published line
table by a reader with no access to this project - they are: the input is
Sugar's two stated figures and the published wavenumber of each line - and no
line is excluded without a row in `line_decisions.csv` saying why. **That
second half is not done and is not this step's to do:** each of the ten lines
needs its plate evidence looked at in IDEN2 first, and until then they are
excluded from the *trend derivation* only, not from the line list.

**Step 2 - the measured hfs components (section 3.5). DONE 2026-09-21.**
*Deliverable:* `hfs_components.csv` (500 rows), `hfs_components.py`,
`hfs_patterns.py`, `hfs_A_from_patterns.csv` and
`tests/test_hfs_components.py`; the interval-rule test at rms 0.090 cm^-1
with no free parameter and 0.084 at a scale factor 0.958 +- 0.006; the
tabulated wavelength identified as the strongest component; kappa(flag) = 1
confirmed. Nothing in the pipeline was run or changed.
*Test:* the A constants fitted from the multi-component patterns agree with
the 22 semiempirical values computed from Reader and Sugar's θ coefficients [R1] wherever
a pattern touches one of those levels. Read this test with section 2.3 in
hand: those values are semiempirical, so a disagreement on a level whose A is a
residue of cancellation is a result about the calculation, not a failure of
the pattern — level 1089 is exactly that case. It is a disagreement on a level
whose A is dominated by a single term, as 1093's is, that condemns the pattern
interpretation, and there nothing built on it is adopted.

**Step 3 - redo the convention fit (section 2.2) on the Step 1 weights, and
run the two physics tests of section 3.3. DONE 2026-09-21.**
*Deliverable:* `hfs_kappa.py` and `tests/test_hfs_kappa.py`. kappa is **not**
zero for the unflagged classes: anchored at the measured kappa(flag) = 1,
kappa(plain) = 0.929 +- 0.014 and kappa(c) = 0.20 +- 0.12, with a permutation
control at -0.003 +- 0.006. The center-of-gravity hypothesis of section 1.1 is
refused for the 1974 lines and the correction applies to every class.
*Test:* the split fits did **not** agree with the pooled value, so the model
is revised as the plan requires: the plain class carries one kappa per era,
0.944 +- 0.014 for 1974 and 0.633 +- 0.035 for 1969, and Step 5 must be
written that way. The branch test passed. The result is consistent with
Step 2: Step 2 anchors this fit rather than being tested by it, and the two
agree on the one thing they share - the size of the flagged lines'
displacement - to within the 0.02 that the anchor's own scale moves it.

**Step 3b - the plate wavelength calibration. DONE 2026-09-22.**
Sugar reduced each plate against comparison lines recorded on that same plate,
so a calibration error is a displacement in *wavelength* shared by the lines of
one plate and by nothing else. Fitted here as one `delta_lambda` per group of
lines inside a plate, jointly with the level energies and the Step 3 kappas -
one weighted least-squares solve, 4373 lines against 803 unknowns, never one
group at a time. The plate joins are the blind stretches of `coverage_gaps.txt`
(`tools/coverage_map.py`), which cut the spectrum into 19 blocks, together with
the hand breaks of `calibration_breaks.txt` (2026-09-23), which cut it into 30;
a group is 25 angstrom bins accumulated until it holds 12 lines, never spanning
a break, and a block that cannot raise 12 lines is given no parameter rather
than its neighbor's value. 152 groups over 24 blocks.
*Deliverable:* `wavelength_calibration.py`, `wavelength_calibration.txt`,
`wavelength_calibration.csv` (the curve) and `wavelength_calibration_points.csv`
(the fit's input, one row per line; the per-block `fit_power.py` files were
dropped on 2026-09-23 at the user's request - the points paste into Excel from
`wavelength_calibration_points.xlsm`). Nothing in the pipeline was run or
changed. See README, "The wavelength calibration of the observed line
list".
*Result:* chi^2/dof 0.701, rms 0.1398 cm^-1; the calibration error has an rms of
0.0048 angstrom over the 151 groups, 48 of them beyond twice their own
uncertainty, and `shift/u_shift` has an rms of 2.1 and runs from -5.9 to +7.0.
It is a real, measured effect, not a fitting artifact.
*Model, settled 2026-09-22 with the user:* the correction is **piecewise**, one
constant per group, and not a smooth polynomial over a plate. Both models are
now fitted on every run and compared block by block (`--model poly` puts one
Legendre polynomial per block, degree chosen by the data, 21 coefficients
against 151 steps). The staircase gives back 241 in chi-square for 124 extra
parameters where noise alone would give 124 +- 16, a 7-sigma preference, so
the staircase is the default and the one to apply. Block by block: 8, 10 and
12 are strongly not smooth (block 10 gives back 106 for 20), 14 is marginal,
and 2, 3, 4, 6, 7, 9, 16 and 18 are smooth.
*Corrected 2026-09-22, on the user's report that fit_power.py fits block 6
with a 4th-degree polynomial at chi2/dof 0.60:* the degree search was greedy
and stopped at the first degree that did not pay, so block 6 - whose degree 1
buys 0.6 and whose degree 4 buys 109 - was reported as degree 0 and as the
worst block in the spectrum. Every degree up to the ceiling is now tried.
*Thin blocks:* a block too thin for a group of its own is given one constant
over the whole block and judged by |u/value|, the user's own criterion from
fit_power.py. All six such blocks fail it (ratios 1.17 to 3.85) and are fixed
at zero; their lines keep the uncertainty the trial measured. Keeping them is
not free - the constant of a 3-line block trades against the energies of the
levels its lines touch, and that near-degeneracy raises the largest
eigenvalue of the calibration covariance by a factor of 66 and the median
u_own_correction from 0.020 to 0.184 cm^-1.
*Deliverable added:* `wavelength_calibration_corrections.csv` - one row per
distinct observed wavenumber of the whole list (6669, of which 6667 are
corrected), with `own_correction` in cm^-1 to add to Sugar's wavenumber and
`u_own_correction` from the full covariance matrix of the fit. rms correction
0.106 cm^-1, median uncertainty 0.017 cm^-1, 2567 lines corrected by more than
twice their own uncertainty. A line the fit never saw takes the nearest group
of its own block; a line in an unparameterized block gets nothing.
*Covariance, checked:* the pseudo-inverse of the weighted normal matrix -
pseudo because the level system is rank deficient by the energy zero. Largest
asymmetry before symmetrizing 5.7e-12 of the largest element; the 151x151
calibration block has eigenvalues 6.1e-08 to 7.4e-04, all positive. The
`u_eff = u/(1-h)` inflation is never used inside the fit and cannot reach the
covariance: it is only a weight for an outside fit to the points file.
*Breaks inside a plate, entered 2026-09-23 on the user's account of the
instrument:* the plates are rigid flat glass bent onto the Rowland circle by two
to four mechanical holders, which pinch the plate and deform the dispersion
curve at each holder without leaving any gap in the line list. Eleven such
breaks, read off `wavelength_calibration_points.xlsm` by the user, are now
carried in `calibration_breaks.txt`: 1750, 2225, 2250, 2283, 2343, 2498, 2543.1,
2631.6, 5370, 7385 and 8190 angstrom. They settle the piecewise-against-smooth
question. With them in, no block shows the staircase ahead of the best
polynomial it can carry by as much as 1 sigma, and over the spectrum the
polynomials are 109 in chi-square better with 37 fewer parameters; the
2191-2781 angstrom stretch, which alone had put the staircase 100 ahead for 18
parameters (13.6 sigma), is entirely accounted for by the seven breaks inside
it. **`--model poly` is therefore the model to apply**, 47 coefficients against
152 steps, chi2/dof 0.690 against 0.698, and the deliverables of 2026-09-23 are
written from it. The plates are smooth between the holders.
*Corrected 2026-09-23:* the smoothness test itself was wrong as well. It
compared the staircase with the polynomial the price rule had *selected*, which
convicts a block of structure whenever its curve is real but too shallow to buy
a degree at 9 chi-square per parameter. It is now made against the best
polynomial the block can carry, at its ceiling degree. That alone took the
2809-3211 angstrom block from "not smooth" to 1.2 sigma, which is where the
user's own parabola fit (chi2/dof 0.426) had put it.
*Test:* holding one group at `delta_lambda = 0` was found to be an unnecessary
constraint - freeing every group raises the rank of the design matrix by one, so
the curve is determined absolutely. The one remaining rank deficiency is the
energy zero, and the null direction taken from the design matrix itself puts
2.65e-16 on any `delta_lambda` against 3.93e-02 on the levels, so the curve is
untouched by it. The weighted mean of the per-line points of each group
reproduces that group's fitted `delta_lambda` to 1e-5 angstrom, over all 151
groups. The common mode - the mean over all groups, with its uncertainty from
the covariance matrix, because the groups share the levels - is
+0.00117 +- 0.00140 angstrom, so the wavelength scale as a whole is not in
question; what is measured is the difference between one plate and the next.

**Step 3c - apply the calibration to the observed wavenumbers (decided by the
user, 2026-09-22; not yet written).** The correction is significant and is to be
applied, not merely carried as an uncertainty. Three things follow, in order:

* **The corrected wavenumber.** `wn_corrected = wn_obs + shift(lambda)`, with
  `shift = -delta_lambda * wn^2 * 1e-8`. Whether `shift` is the staircase of
  `wavelength_calibration.csv` or a smooth function of wavelength fitted per
  block with `fit_power.py` was the user's to decide, and it is decided: the
  per-block polynomial, for the reason measured under Step 3b once the plate
  holders were allowed for. `own_correction` and
  `u_own_correction` of `wavelength_calibration_corrections.csv` are the columns
  to apply. A line in a block with no fitted group gets no correction.
* **The statistical uncertainty of a corrected wavenumber**, re-derived per era
  per character, since removing the calibration error removes part of what the
  Step 1 uncertainties were absorbing. The table is in
  `wavelength_calibration.txt`, "The statistical uncertainty of a corrected
  wavenumber": the rms of the fit's residuals over the class's effective degrees
  of freedom `sum(1 - h)`, floored at `hfs_kappa.FLOOR`. It runs 0.78 to 1.02 of
  the presently adopted value, the 1969 plain class falling the most
  (0.2450 -> 0.1900 cm^-1). This is what goes in the LOPT input as the input
  wavenumber uncertainty.
* **The systematic uncertainty per group**, `u_shift_cm-1` of
  `wavelength_calibration.csv`, which LOPT can carry as a function of wavelength
  for each assigned group of lines. It must stay separate from the statistical
  part: a calibration error displaces every line of a group the same way, so it
  does not average down over the lines of a level.
* **For IDEN2**, `dlv.dat` wants one total uncertainty per wavelength - the
  statistical and the systematic parts in quadrature - as the basis for visual
  identification. The rewrite of `dlv.dat` is already on the section 6 list
  under `unc_own`.

*Deliverable (proposed):* the correction and both uncertainties applied in
`make_LOPT_input.py` behind a config switch, beside the hfs switch of Step 5, so
that the two corrections are turned on together and `Pr3_lines.xlsx` stays the
record of measurement.
*Test (proposed):* with the switch off, `LOPT_input_lines.txt` is byte-identical
to the current file; with it on, only the wavenumber and uncertainty columns
differ, and the fit's chi^2/dof moves toward 1 rather than away from it.

**Step 4 - A constants from the observed data, with a calculation as a
check (section 3.1); resolve the two conflicts (section 3.2) by hand with
IDEN2.** *Rewritten 2026-10-03.* The earlier title, "extend the A
calculation to LS-coupled levels", made the calculation the source of the
values. It is now the check. A line is moved only by an A constant that
rests on measurement, or on a calculation the measurements have confirmed
(the user's rule, section 3.1).

*Scope:*

1. *The 39 undetermined levels.* Their A comes from the data, by the routes
   of section 3.1. First is the global pattern fit, rerun on the current
   assignments with the confirmed semiempirical constants held at their scaled
   values; then the flags, the classes and the widths as checks. A value is
   adopted only where the data bound it from both sides. The rest stay "not
   determined", and their lines stay unmoved with the D10 allowance.
2. *The semiempirical constants the data contradict or do not test.*
   `059003.000127` gets its value from its own patterns. The 7 untested
   4f²6p levels of section 3.1 keep their semiempirical values with a wider
   u_A; how much wider is the user's to choose.
3. *The hand-measured constants that used a partner's semiempirical A.* Four
   sources in `A_hfs_levels.csv` name the partner's value, A(000187),
   A(000140) or A(000161), unscaled. They are recomputed with the partner on
   the measured scale of Step 2 (0.9638).
4. *The check calculation.* It is made only once 1 is done. It needs the LS
   composition tables from the user, and a_5d estimated two ways: (a) scaled
   from a_4f and a_6s by the ratios of the radial integrals of the Cowan
   calculation; (b) fitted to the values of 1. It is a check on 1. It does
   not replace a value from 1, and it moves no line unless the user decides
   otherwise.

*Deliverable:* `A_hfs_levels.csv` with the new values, each with a source
naming the data it comes from; the generating script moved out of the
scratchpad into `LineClass/`, with a test file; ledger rows in
`line_decisions.csv`, one per conflicting line examined.

*Test:*

* every adopted value is reproduced from its own data, within its
  uncertainty;
* `059003.000233` is reproduced within the combined uncertainty;
* flag-direction agreement is no worse than 261 of 264;
* afterwards no level has flag constraints that contradict its A by more
  than its uncertainty;
* for the check calculation, the 22 values from Reader and Sugar's θ
  coefficients are still reproduced at rms 0.005 cm^-1 or better. That is a
  check of algebra against algebra (section 2.3), not a test against
  experiment.

*Items 1 to 3 done 2026-10-04* (`hfs_A_fit.py`, `tests/test_hfs_A_fit.py`,
the per-level record `hfs_A_anchored.csv`).

*The fit.* Sugar's components were joined to the iter_hfs assignments on
wn_key (`hfs_components.read_line_list` now does that wherever a table has
wn_key). Each listed position is one equation of section 3.5; the 73
semiempirical constants the spacings confirm enter as priors at their scaled
values and uncertainties; every other level met by a pattern is free. The
weight of a position is Sugar's 0.003 A, taken twice for a spacing.

```
patterns used                     218   (25 blends or unclassified, 4 too long for the ladder,
                                          1 set aside by hand)
positions                         440
chi2/dof                          0.61  (1.04 before the one pattern was set aside)
rms of the positions              0.027 cm^-1
priors pulled beyond 2 sigma      3 of 73, the largest 2.9 (000229)
levels on a singular direction    000034, 000076, 000078 (tied only to one another)
```

The pattern set aside is 29753.592 `*v` (000182-000274). Its one listed
companion, +1.188, is rung 2 of the pattern (+1.06 to +1.20 predicted for
the A its other pattern gives), and rung 1 (+0.62) is not listed. Read as
rung 1 it put A(000274) 11 sigma off its other pattern and pulled the priors
of 000182 and 000183 by 3 sigma. No other pattern fits better with its
rungs relabeled.

*Adopted into `A_hfs_levels.csv`:* 46 levels, with the source "resolved
components, N patterns, M positions (Step 4 anchored fit, 2026-10-04)" and
u_A the fit's uncertainty with A * 0.0079 (the scale's) added:

* 35 of the 39 undetermined levels. They come out between -0.001 and
  +0.037 cm^-1, with uncertainties of 0.001 to 0.011. The unanchored fit of
  section 3.5 had them 0.03 higher, which was its free offset;
* the 7 flag-interval levels, whose intervals the values all fall inside,
  except 000233 (below);
* `059003.000127`: +0.0775 +- 0.0014 from 13 patterns, against the
  semiempirical +0.0935 after scaling (17 per cent smaller). Its own patterns
  agree within 1.4 sigma.

*Not adopted:* `059003.000448` (f5d6s, J = 9/2) is held (`HOLD`). The fit
gives +0.088 +- 0.007 from a single position, 34042.649 `*r`, which nothing
else checks, and its own 15 lines give +0.013 +- 0.018. 000258, 000295 and
000453 are reached by no usable pattern. These 4 stay "not determined".

*000233 is not reproduced, and the old value was wrong.* Its flag value,
-0.0462 +- 0.0040, came from three flags. The bound A < -0.043 came from
31558.079 `*v` (to 000149), which was rejected on 2026-09-29, after the
value was set. Its two remaining flags allow -0.052 < A < +0.130. Its three
component positions, through two partners (000140 and 000142), give +0.021
to +0.024 each; the fit gives +0.0193 +- 0.0044. Until now the pipeline
moved 000233's lines with the stale value.

*The checks:*

```
flag direction, flagged lines with one accepted classification
   with the constants before Step 4      181 of 181 testable agree
   with the fitted constants             253 of 253 testable agree
flag bounds against the fitted values    none outside (2 sigma)
each pattern alone against the joint     worst 1.9 sigma (000293), 24 levels with 2+ patterns
own lines (route 3), 47 levels           44 within 3 sigma; 000127 +3.4, 000243 -4.7,
                                          000448 -3.9
satellites (hfs_satellites.txt)          mostly predicted closer than before
```

The 261 of 264 of section 2.4 was counted on an older classification; on
the current one no flag disagrees, before or after. Route 3 is weak and
assumes the kappa model. For 000127 it is circular, since the kappas were
fitted with 000127's semiempirical D. For 000243 it is outweighed by two long
ladders (5 and 4 rungs) that fit to 0.03 cm^-1. 000448 is the case held
back.

*Item 3.* The four hand-measured constants that used a partner's calculated
A were recomputed with the partner scaled by 0.9638. The old values are
reproduced with the old partners:

```
level    partner A, before -> scaled     A before    A now      u_A
000260   000187 +0.2053 -> +0.1979       +0.0089     +0.0013    0.0080
000432   000140 -0.0826 -> -0.0796       +0.0383     +0.0416    0.0180
000439   000161 +0.0942 -> +0.0908       +0.0200     +0.0172    0.0120
000629   000161 +0.0942 -> +0.0908       +0.0199     +0.0155    0.0240
```

*Found on the way, not acted on:*

* Six levels met by patterns were absent from `A_hfs_levels.csv`, so the
  pipeline gave them S = 0 and kept their widths. The fit determines three
  of them, and on the user's decision they were added to the table on
  2026-10-04 with the Step 4 source (n_flagged: their *r and *v lines in
  iter_hfs): 000151 (4f2.5d, -0.0683 +- 0.0057, one position), 000654
  (+0.0173 +- 0.0036) and 000656 (+0.0140 +- 0.0043). 000151 rests on one
  position like 000448, but its own lines disagree less (-0.023 +- 0.016,
  2.7 sigma against 4.2 for 000448); a first-order estimate on its 36 lines
  takes chi2 from 34.1 to 40.0, the worst line 83051.0561 from 2.4 to 3.3
  sigma. The other three are the 000034-000076-000078 cluster, which no
  prior reaches.
* Eleven rows of `inflated_unc_lines.txt` with the reason "hfs of f26s"
  now have every level determined: 37577.0988, 39960.6491, 30260.0502,
  35327.7691, 39719.0502, 33934.5817, 31478.9216, 36647.5418, 44494.0809,
  35266.469 and 44570.63 (first counted as nine). The constants found are
  small, and the hfs explains none of these residuals. Corrected on
  2026-10-04: in iter and iter_hfs the D13 release does not reach them.
  `wavelength_calibration.py` writes every registry value, tagged hfs or
  not, into `unc_own_corr`, so classify_lines finds nothing to widen and
  records nothing to give back (all 34 tagged rows carry the registry
  value and hfs_allowance = 0 in iter_hfs, the kappa = 1 line 29056.1776
  included). The release acts only on a row added after the last
  calibration run, as 36302.172 was on 2026-10-03. The tag is harmless to
  these rows as the code stands, but no longer true.
  On 2026-10-05 the user reworded the eleven reasons without "hfs". The
  two c lines, 30260.0502 and 35327.7691, left the registry, since the
  c class fits them once A(000250) is known. The other nine were then
  checked against the full component patterns (the scratchpad's
  step4/peak13.py and pair190.py):
  - 37577.0988 sits where a blend of head and rung 1 would (-0.16 from
    the head);
  - 31478.9216 sits between the head and the center of gravity;
  - the six w/d lines sit at or a little past the center of gravity
    (kappa_line -0.4 to 0.2). So do all 8 other w/d lines of 1974 with
    |D| > 0.3 and both A known;
  - 33934.5817 and 27443.3325 are the centroids of the F = 3 and F = 2
    groups of 000190 (-0.213 and -0.234 predicted, -0.203 and -0.230
    observed).

  The user decided that w/d/h lines do not form one class (some are
  measured near the head). Instead, `kappa_exceptions.txt`
  (hfs_correction.py) fixes kappa per transition of a line, at head,
  cg or an existing class's value, and can count a resolved level with
  its whole pattern for that line (`unresolved`). The eight rows above
  were entered as cg. In the same change `wavelength_calibration.py`
  began reading lines through the set's hfs model: resolved levels
  count S = 0, the satellite heads and the exceptions are held, not
  fitted. A(000190) became the measured -0.157 +- 0.010.

*Item 4*, the check calculation, waits for the LS composition tables.


**Step 5 - put the correction in the pipeline (decisions D1, D11). Agreed
2026-09-29: the head frame. Implemented 2026-09-30 in `classify_lines.py` and
`make_LOPT_input.py`, switched on by `[hfs] apply`; off in every set until
the user turns it on.**

*The agreed design.* The user's objection to every earlier version (history
below): the hfs shift belongs to a transition, so an unidentified line cannot
be shifted, and IDEN2, the tool new identifications are made with, can only
apply per-level shifts. So the pipeline works in the head frame:

* The working level energies are the F = I+J sublevels, E_head = E_cg +
  I*A*J. This is nearly what the pipeline fits already, because the strongest
  component of every transition joins the F = I+J sublevels of both levels and
  Sugar measured the strongest component (D7, kappa = 1 for flags).
* A line of class kappa then lies (1-kappa)*D below the head-frame Ritz value:
  exactly on it for flags, within 0.056*D for 1974 plain lines, 0.367*D for
  1969 lines, 0.80*D for `c` lines. In `iter/line_classifications.csv` (A
  unknown counted as 0), 440 of 4965 accepted lines need more than 0.1 cm^-1
  in this frame, against 1981 in the center-of-gravity frame; 336 need more
  than their own uncertainty, against 1964.
* The calibration stays what it is: an assignment-independent correction of
  every observed wavenumber (`own_corr`), already carried into `dlv.dat` by
  `sync_IDEN2.py`.
* `classify_lines.py`: each candidate is predicted at Ritz_head -
  (1-kappa)*D_i, and the line gets hfs_shift = (1-kappa)*sum BF_i*D_i.
  `make_LOPT_input.py` adds hfs_shift to the wavenumber it writes.
  `sync_IDEN2.py` is unchanged in what it writes: `enlev.dat` holds the
  head-frame energies and `dlv.dat` the calibrated, unshifted lines.
* The centers of gravity are computed once, at the end, E_cg = E_head - I*A*J
  with I*J*u_A added in quadrature, for the published level list and the
  Cowan fit. That is the only place the absolute A values enter.
* One kappa per observed line, since kappa is a property of the line (flag,
  era, character), with one exception decided by the user on 2026-09-30: a
  flagged blend takes kappa as the BF-weighted mean of its components'. The
  flag describes the pattern of one transition, taken to be the strongest
  (largest BF; Sugar's measured components in `hfs_components.csv` agree
  wherever both components have an A to test, except 28477.539, where the
  weaker 000148-000211 ladder fits them better but the stronger one's A,
  000284's, is not determined). That one takes kappa = 1, the others the
  plain kappa of the era. Each component carries its own (1 - kappa_i)*D_i,
  shift = sum BF_i*(1 - kappa_i)*D_i, which is the BF-weighted mean of where
  the components were measured; the `kappa` column shows sum BF_i*kappa_i.
  In the iter/ scratch run the 8 flagged blends move by -0.019 to +0.012
  cm^-1. The classification itself still compares every candidate of a
  flagged line with its Ritz value (kappa = 1): which component is the
  strongest is only known once it has converged.
* Consequence for 000642 (IDEN2 1093, J = 1/2, A = +0.505): the head frame
  puts it at F = 3, 1.53 cm^-1 above its former F = 2 position (the
  alternate `level_positions.py` prefers). Done 2026-09-29: `new_levels.txt`
  and the ledger (7 F = 3 lines accepted, the F = 2 lines rejected); the
  iter LOPT fit of 9/30 gives 49653.872 +- 0.015 from 9 lines. Its own A
  counts as 0 in D, since each line ends on one sublevel; the A enters only
  at the c.g. step, E_cg = E(F=3) - 0.631. 000190 (IDEN2 1089, J = 1/2,
  A about -0.157) is at F = 3 already; added to `resolved_levels` on
  2026-09-30 at the user's request, which takes 0.0055 cm^-1 off the
  shifts of seven 1974 plain lines and 0.079 off two `c` lines. Sugar's 500 printed components are not rows of the line list, so
  they consume nothing; resolved companions listed as separate lines would be
  classified as hfs companions of their transition, kept out of LOPT.

*What was written, 2026-09-30.*

* `lineclass_config.toml` gains `files.hfs_A_levels` (`A_hfs_levels.csv`) and
  an `[hfs]` section: `apply` (false in the baseline), `resolved_levels`
  (`["059003.000642", "059003.000190"]`), and `[hfs.kappa]` with `[kappa, u_kappa]` per class -
  `flag` [1.0, 0], `plain_1974` [0.944, 0.014], `plain_1969` [0.633, 0.035],
  `c` [0.20, 0.12] (D12). A set turns the correction on with `[hfs]
  apply = true` in its own configuration; everything else is inherited.
* `hfs_correction.py` (new) holds the model. A level is *determined* (A from
  Reader and Sugar, the composition or the flag intervals: S = I*A*J, u_S =
  I*J*u_A), *undetermined* (the 41 levels "not determined": S = 0, and the
  unknown A carried as u_S = I*J*u_A, D10) or *absent* from the table (S = 0,
  no uncertainty; its hyperfine width still covers it). A level in
  `resolved_levels` counts S = 0. The uncertainty of a line's shift is
  (1-kappa)*sum BF_i*u_D_i, summed linearly because the components of one line
  share levels, in quadrature with u_kappa*|sum BF_i*D_i|.
* `classify_lines.py`: every comparison of a line with a Ritz value (matching
  window, grades, the blend tests of the weeding, O-C) uses Ritz - (1-kappa)*D,
  and every energy determination (level optimization, u_energy) uses the
  measured wavenumber + (1-kappa)*D. With the switch on, the table gains four
  columns: `kappa`, `hfs_D` (per row), `hfs_shift` and `u_hfs_shift` (per
  line, 0 for a line with nothing accepted); `dif_wn_O-C` is the residual
  against the corrected prediction. With it off the table is byte-identical
  to what the code wrote before (checked on `iter/`).
* `make_LOPT_input.py` reads the switch from the configuration beside the
  table, writes every record of a line at wn_obs + hfs_shift, adds
  u_hfs_shift to the line's uncertainty in quadrature, and writes
  `LOPT_hfs_shifts.txt` beside `LOPT_input_lines.txt`: one row per record
  whose wavenumber changed. A table written under the other setting stops the
  run. With the switch off the LOPT input is byte-identical to before, and a
  shift record left from an earlier run is removed.
* `sync_IDEN2.py` and `check_sync.py` find the observed line of a LOPT record
  by its wavenumber (the `dlv.dat` row; the classification row), so both take
  the measured value back from `LOPT_hfs_shifts.txt`. Without that, 1930
  accepted transitions of the scratch test below would have found no row in
  `dlv.dat`.
* D13, as implemented: a level with a determined A (or resolved sublevels)
  no longer carries its hyperfine width from `level_hfs_widths.csv` into
  LOPT's uncertainty, and a flagged line carries none; a level without a
  determined A keeps its width. A line widened by an `inflated_unc_lines.txt`
  row whose reason says hfs gets the line list's own uncertainty back, at the
  end of the classification run, when every accepted transition of it touches
  only such levels or the line is flagged; the other lines keep the
  allowance. The registry itself is not touched.
* D13 amended 2026-10-01 (user's approval): an *undetermined* level's hfs is
  carried once. `classify_lines.py` writes `u_hfs_undet` (the part of
  u_hfs_shift owed to such levels) and `hfs_allowance` (the line keeps a
  registry allowance for hfs); `make_LOPT_input.py` adds the larger of
  u_hfs_undet and the level's width, or neither under a kept allowance. An
  absent level keeps its width. Scratch estimate on the iter_hfs table of
  2026-10-01: 124 fitted records lose up to 20% of their uncertainty,
  RSS/dof stays 0.77, no level moves more than 0.005 cm^-1.
* `tests/test_hfs_correction.py`, 22 tests.

*Test, run 2026-09-30 in the scratchpad (nothing in `iter/` or the project was
written).* With the switch off, `classify_lines.py` on the iter inputs
reproduces `iter/line_classifications.csv` byte for byte, and
`make_LOPT_input.py` gives the same file as the code before the change. With
it on: 88 levels carry an S; 5554 accepted assignments against 5549 (14 lost,
19 gained, to be reviewed before any of it is adopted); rms O-C of the singly
assigned lines 0.0792 -> 0.0732 cm^-1 for flags, 0.0467 -> 0.0464 for 1974
plain, 0.2136 -> 0.2024 for 1969 plain; 1727 observed lines shifted in the
LOPT input, 402 by more than 0.1 cm^-1, the largest by 0.910; 8 hfs registry
allowances withdrawn, 25 kept; levels moved by 0.032 cm^-1 rms in
`classify_lines.py`'s own fit, most by 0.41 (000304). `check_sync.py` finds
every classified transition in the shifted LOPT input, and `sync_IDEN2.py
--dry-run` finds a `dlv.dat` row for every accepted one. LOPT itself was not
run on it: that is Step 6.

*Satellites of five 4f2.6s levels, examined 2026-09-30.* The user asked
whether the weaker (or equally strong) unassigned line beside one or more
transitions of IDEN2 rows 1138, 1131, 1129, 1118 and 1091 is a resolved hfs
component, and if so whether the level belongs in `resolved_levels`. The
pattern of each transition was computed with `hfs_patterns.components`
(A from `A_hfs_levels.csv`; for a partner whose A is not determined, the A
that would put the second component on the satellite is given instead).
"rung 1" is the second principal component, F_low - 1 -> F_upp - 1, relative
to the strongest, F_low = I+J_low -> F_upp = I+J_upp; its predicted strength
relative to the strongest is 0.75-0.83 for these J. Offsets in cm^-1, positive
to the violet, observed from `iter/line_classifications.csv`; intensities are
Sugar's own.

```
row   transition (low -> upp)   line        satellite   offset   rung 1 predicted            I sat/line
1138  000040 -> 000141          20971.778   20972.130   +0.352   needs A(000040) = +0.022    2/3
1131  000148 -> 000253          34139.985   34140.416   +0.431   needs A(000253) = +0.039    3/4
1129  000150 -> 000247          32742.923   32742.240   -0.683   needs A(000247) = +0.069    8/8
1118  000161 -> 000215          25787.266   25787.757   +0.492   +0.461 (A both known)       3/3
1118  000161 -> 000205          24110.254   24110.736   +0.482   +0.481; rung 2 +0.890       4/8
                                            24111.152   +0.91      against +0.91 observed    3/8
1118  000161 -> 000204          24087.303   24087.714   +0.412   needs A(000204) = +0.035    3/6
1118  000161 -> 000203          23792.532   23792.956   +0.424   +0.452 (A both known)       1/1
1091  000124 -> 000188          22545.679   22546.015   +0.336   needs A(000124) = +0.024    1/1
```

All are resolved hfs components. The four rows of 1118 whose partner A is
known match rung 1 to 0.03 cm^-1, and at 24110 rung 2 as well. The needed A
values of the undetermined partners, +0.022 to +0.039, lie where the other
4f5d2 and 4f3 levels do. 1129 is the one doubtful case: rung 1 would need
A(000247) = +0.069, whereas with the composition estimate +0.025 the
satellite falls on rung 2 (-0.687). The main line would then be a blend of
the strongest component and rung 1, not the strongest component alone.

None of the five levels belongs in `resolved_levels`. Along this ladder both
levels change F together. So the listed main line is the strongest component
of the whole pattern, exactly what a Sugar `*r`/`*v` line is, and needs no
shift (kappa = 1). Making the level "resolved" would instead set its own S to
0 and leave (1 - kappa)*S(partner). For 000161 -> 000205 that is +0.017
cm^-1, where the correct value is 0 and the plain-line value is -0.042.
Worse, the level's other, unresolved lines would lose their correct shift.
What the observation calls for is a per-line marker: kappa = 1 for these
main lines, and the satellites classed as hfs companions of their
transition, kept out of LOPT.

*The marker, agreed and written 2026-09-30.* The user agreed to both, noting
that the satellites, once accounted for, can also improve the derived A
constants. `hfs_satellites.txt` (`files.hfs_satellites`, kept by hand; the
columns wn_key, main_wn_key, low_id, upp_id, rung, date, reason) holds the
eight companions of the table above; 1129's is a comment row until its rung
is settled. `classify_lines.attach_hfs_satellites` reads it in the real run
only, as the ledger is read. The companion gets a row of its transition with
grade `hfs`, accepted = 0, its offset from the Ritz value in `dif_wn_O-C`,
and no candidate of its own; `make_LOPT_input.read_classifications` drops
that row. The main line's named transition is compared with the head-frame
Ritz value (`Transition.hfs_offset` = 0 for a pair in
`SpectralLine.hfs_head_pairs`), takes kappa = 1 in the line's shift
(`Model.component_kappas`, the other components of a blend keeping the kappa
they would have without the mark), gives up its level widths in
`make_LOPT_input.py` and its hfs registry allowance, as a flagged line does.

Scratch runs on the iter inputs (nothing in `iter/` written): with the switch
off, the only change is the eight companion rows; accepted assignments,
level energies and the LOPT input are unchanged. With it on, accepted
assignments are unchanged (5556), and the seven main lines lose shifts of
-0.014 to -0.060 cm^-1 and their level widths (34139.985: 0.043 -> 0.032
cm^-1 in the LOPT input). The levels move by at most 0.031 cm^-1 (000161).
How the lines of the levels concerned fit, chi^2/n of the accepted lines
(O-C over the line's uncertainty) before -> after:

```
level    n    before  after   main lines' O-C after
000148   32   1.09    1.02    -0.001  (was -0.049)
000253   23   0.73    0.65    -0.001
000040   15   0.74    0.70    -0.009  (was -0.018)
000141   31   1.38    1.35    -0.009
000124   10   0.42    0.48    +0.019  (was +0.012)
000188   13   0.47    0.53    +0.019
000161   27   2.11    2.58    -0.034 -0.059 +0.041 +0.050
```

000161 is the exception. Its four main lines (0.022-0.040 cm^-1) are more
precise than its 13 flagged `*v` lines (0.074), and at kappa = 1 they pull
the level down by 0.031, which leaves the `*v` lines about 0.03 low; much of
the rise in chi^2 is one line already 5.5 sigma out, 40015.2007 (000161 -
000281, `**`); without it, 1.05 -> 1.20. Putting every line of 000161 at
kappa = 1 is worse still (3.35), so the rest of its lines do want the plain
kappa. Either the four main lines were measured a little below their heads,
or the upper levels 000203, 000204, 000205 and 000215 are misplaced by a few
hundredths; the rungs themselves agree with the A constants (000205 needs
A = +0.0222 against +0.0223 listed). Left for the user to judge.

The A constants. Each rung-k interval gives k/2 * ((2F_low + 1 - k) A_low -
(2F_upp + 1 - k) A_upp), one linear equation in the two A values
(`hfs_patterns.rung`). The registry holds what is needed to add these to the
pattern fit of `hfs_A_from_patterns.csv` and to Step 4; not done yet.

*Not yet hfs-aware.* `insert_new_level.py`, `move_level.py` and
`discard_level.py` edit the baseline's `LOPT_input_lines.txt` record by
record and do not read the shift record; they must not be used on a set with
the switch on. `level_positions.py` and `review_mismatches.py` recompute
residuals as wn_obs - Ritz from the table and do not yet subtract the
correction.

*History: the center-of-gravity frame, superseded.* The first version of this
step (D1, D11 of 2026-09-29) put the levels at their centers of gravity and
corrected the observed lines by kappa*D - in LOPT's input, in
`classify_lines.py` (Ritz + kappa*D per candidate) and in the set's
`dlv.dat` - one shift per line, kappa * sum BF_i*D_i. LOPT alone was not
enough: the displacement kappa*D exceeds 0.1 cm^-1 for 1507 of 3834 lines and
reaches 1.23 cm^-1, against 0.01-0.04 cm^-1 for a sharp 1974 line, and a fit
on corrected input would have put centers of gravity into `enlev.dat` and
`new_levels.txt`, against which IDEN2 and `classify_lines.py` would have seen
every 4f2.6s line off by its own kappa*D. Putting the correction into the
level energies instead (E_cg + I*A*J in `enlev.dat`, lines as measured) was
withdrawn because a level shift cannot carry kappa, which belongs to the line:
plain lines would have stayed (1-kappa)*D off. The user then objected that an
unidentified line cannot be shifted either, and the head frame above replaced
both.

**Step 6 - run LOPT both ways and compare.** Per CLAUDE.md, `lopt.bat
LOPT.par` from `LineClass/`, never the Java jar.
*Deliverable:* the two level lists, and a difference table.
*Test:* the standalone fit numbers are reproduced by LOPT within rounding. If
LOPT disagrees materially, the standalone fit is wrong and the discrepancy is
resolved before going on.

*Done 2026-10-03, in the scratchpad (`step6/`).* Both fits were built from
the iter_hfs table and the pass-3 line list, so the assignments (5601
weighted lines) and the calibration are the same and only the correction
differs. The "on" build reproduces iter_hfs's own LOPT input and output
exactly. The "off" build is the same table with its six hfs columns removed,
so make_LOPT_input.py treats it as uncorrected: measured wavenumbers, with
the level widths of `level_hfs_widths.csv` on every line. It differs from a
true measured-frame classification in one way. The 8 registry allowances
the correction withdraws are missing from it as well.

```
                           on (head frame)   off (measured)
RSS/dof, 4373 dof               0.80              0.91

singly assigned lines (4498), O-C rms in cm^-1 and in units of its uncertainty
class            n      on       off      on/u    off/u   >3 sigma on  off
flagged *r/*v   260   0.0716   0.0924   0.919   0.878        0      0
plain 1974     1672   0.0493   0.0499   0.781   0.798        2      2
plain 1969     1798   0.1816   0.1881   0.873   0.933       10     10
c                68   0.1881   0.3711   0.821   1.810        0      7
other chars     700   0.2018   0.2154   0.840   0.871        0      0
```

The levels hardly move between the two frames: E_on - E_off has rms
0.025 cm^-1 over 653 levels, 8 levels beyond 0.1 and none beyond 0.3. The
largest is 000451 at +0.20 +- 0.13; its A is only estimated. The direction
follows the physics. Of the 78 levels with a nonzero A that move by more
than 0.02, 77 move toward their head, by the sign of A. The ten that move by
more than 3 times the combined D1 are 000127 (4f2.6s, +0.14) and nine
4f2.6p levels (+0.033 to +0.050). The size of the move is roughly
(1-kappa)*I*A*J, the part of the pattern the plain lines were not measured at.
The combined D1 ignores the correlation of the two fits, so 3 sigma here
overstates the significance.

*The test against the standalone fit.* The calibration's
center-of-gravity Ritz value (ritz_cm1 - hfs_shift_cm1 of
`wavelength_calibration_points.csv`) was compared with LOPT's head-frame
Ritz value minus 0.9582*D on the 4494 lines of the calibration fit. The
mean difference is -0.0003 cm^-1 and the rms 0.018. The median is 0.22 of
LOPT's Ritz uncertainty, and 11 lines lie beyond 3. Those 11 sit on a few
levels (000569, 000504, 000508) whose LOPT positions also rest on blended
lines the calibration leaves out. The standalone fit is reproduced. It is
not reproduced to rounding, because the two fits use different line sets
and weights.

*An inconsistency found by the test.* The kappas were fitted with D scaled
by A_SCALE = 0.9582 (`hfs_kappa.py`), but `hfs_correction.py` applies
(1-kappa)*D with D unscaled. The difference is (1-kappa)*0.042*D. That is
at most about 0.004 cm^-1 for a plain line, but up to about 0.03 cm^-1 for
a `c` line. One of the two should change; that is for the user to decide.
**Resolved 2026-10-03 (decision D14).** Both programs now read A through
`hfs_kappa.scaled_A`. It scales only the semiempirical constants, by 0.9638
(section 3.5), and adds the scale's uncertainty to u_A. The kappas must be
refitted with D on that scale. That is calibration pass 4 of Step 7, and
the iter_hfs chain follows it.

**Step 7 - second calibration pass** (see section 6). The character constants
of `level_positions.py` were fitted against the uncorrected residuals and must
be iterated to convergence against the corrected fit.

*Started 2026-10-03.* `wavelength_calibration.py --classifications PATH` reads
a working set's table, on its wn_key (Sugar's published wavenumber), never on
its calibrated wn_obs; without the option it still reads the frozen baseline
table of 2026-09-24. A preview on `iter_hfs/line_classifications.csv`, written
to the scratchpad only: 4494 lines against 4370, 152 groups, chi^2/dof 0.969,
rms 0.1471 cm^-1; kappa(plain 1974) 0.925 +- 0.012 (was 0.947), kappa(plain
1969) 0.708 +- 0.021 (0.718), kappa(c) 0.271 +- 0.044 (0.185). No line's
correction moves by more than twice the combined systematic uncertainty of the
two passes; the largest change is +0.085 cm^-1 (1969, near 55500 cm^-1). The
systematic uncertainties fall by a fifth to a third. The 1974 `d` class widens from
0.042 to 0.071 cm^-1 on the lines it measures itself (registry lines left
out), driven by lines of the newly found levels. Not yet written: the user
runs it.

*Calibration converged 2026-10-03, by the user's judgment.* Passes 2 and 3
were written on the iter_hfs table, and the iter_hfs chain was rerun after
each one. Pass 3 changed no line's correction by 0.01 cm^-1 or more; 96 lines
moved by more than 0.005, all of them in the 3075-3100 A group. That group
moved by +0.008 cm^-1, 0.6 of its own uncertainty, because one newly accepted
line joined it. chi^2/dof stayed at 0.969, and the kappas held at 0.924(12)
for plain 1974, 0.708(21) for plain 1969 and 0.271(44) for `c`. LOPT's RSS/dof
on iter_hfs is 0.80. The other half of this step, refitting the character
constants of `level_positions.py`, is still open. Making
`level_positions.py` subtract the correction would allow it on iter_hfs (see
*Not yet hfs-aware* under Step 5). The user's Set 3 (Step 7b) removes that
need: its lines are already at their centers of gravity. So the constants
are refitted on Set 3. (Withdrawn 2026-10-11: no level analysis in Set 3; see the end
of Step 7b.)

*Pass 4, 2026-10-03, after decision D14.* Run by the user on
`iter_hfs/line_classifications.csv`, with only the semiempirical A constants
scaled, by 0.9638, instead of every D by 0.9582. 4496 lines, chi^2/dof 0.969,
rms 0.1472 cm^-1. The kappas are 0.925(12) for plain 1974, 0.712(20) for
plain 1969 and 0.276(44) for `c`; they moved by 0.2 sigma or less, and were
copied into `[hfs.kappa]`. No line's calibration correction changed by more
than 0.0025 cm^-1 (43 lines by more than 0.002, near 33050 cm^-1), a fifth of
its uncertainty. The iter_hfs chain has to be rerun on the new kappas and the
new scale.

*Pass 5, 2026-10-04.* Run by the user after the iter_hfs chain
was rerun on pass 4, with 83315.6988 (000139-000630) inflated to 0.5 cm^-1:
it is the center of gravity of its pattern, not its head. 4496 lines,
chi^2/dof 0.969, rms 0.1473 cm^-1. The kappas are 0.926(12), 0.717(21) and
0.276(44), copied into `[hfs.kappa]`. Only the 1200-1225 A group moved: its
correction went from +0.0586(425) to +0.0420(430) cm^-1, 0.4 of its
uncertainty; no other group moved by more than 0.0011 cm^-1.

**Step 7b - Set 3, the published set (added by the user 2026-10-03).** This
is the final set: calibration corrected, levels at their centers of gravity,
and every line reduced to its center of gravity. It replaces Set 2 as the
set to publish. Its line list differs from Sugar's in four ways:

* *The hfs components are merged.* The newly identified companions of
  `hfs_satellites.txt` and Sugar's printed components (`hfs_components.csv`)
  do not appear as lines of their own; they are merged into one line at the
  center of gravity of the pattern. The same goes for the widely split
  levels carried as two positions, such as 000642 (IDEN2 row 1093, F=3 and
  F=2). The level becomes one level at its center of gravity, and the lines
  of both components become single lines at their centers of gravity.
* *Flagged lines are moved to their centers of gravity.* A `*r` or `*v`
  line, measured at its head, is moved to head - D. Its character becomes
  `hfs,l` or `hfs,s`: the pattern shaded toward longer or shorter
  wavelengths. Sugar's r and v (red, violet) give these directly, and the
  sign of D gives them independently. The few lines where the two disagree
  (section 2.4: 261 of 264 agree) are listed for the user.
* *Other characters become `hfs` where hyperfine structure is shown to be
  their cause.* Some `c`, `cl`, `w`, `d`, `h` and `ch` lines are
  hyperfine patterns. The criterion is the user's to choose. A candidate is
  a computed pattern width that is large compared with the plate resolution,
  together with a residual that improves when the correction is applied.
* *Unflagged lines with a known mean offset are unbiased.* A line of a class
  whose mean position in its pattern is measured (kappa) is moved by
  -kappa*D. Its uncertainty gains u(kappa)*D and the uncertainty of D, which
  comes from u_A.

Points to settle before building it:

* the lines of the levels without a determined A: **settled 2026-10-03.**
  They stay at their measured positions and keep the allowance they carry
  now, until Step 4 gives a level an A from the data (section 3.1);
* which D: **settled 2026-10-03 (D14).** Semiempirical constants are scaled by
  0.9638, measured ones are not, the same as in the pipeline and the kappa
  fit;
* whether the center of gravity of a flagged line with measured components
  comes from those components (intensity-weighted with the calculated
  relative intensities) or from head - D, and what is done where the two
  disagree;
* the intensity of a merged line: the sum of its components;
* the directory: `final/` exists, holding only an IDEN2 copy.

*Built 2026-10-09: `make_set3.py`, README "Set 3, the published set".* `final/` is a
derived set (`derived_from = "make_set3.py"`): iter_hfs's table and fit re-expressed at the
centers of gravity, never classified. The open points were settled with the user:

* **the cg of a flagged line is head − D for every line.** Sugar's printed components are
  only a check. Over 228 lines their per-line cg differs from head − D by 0.021 cm^-1 at the
  median, and 20 lines beyond 3σ are listed in `final/set3_report.txt`;
* **a `c`, `cl`, `w`, `d`, `h` or `ch` line becomes `hfs,l`/`hfs,s` when |D| ≥ one effective
  line width** (0.035 Å with Doppler). 78 lines do now;
* **the F = 2 partners of a resolved level's lines are paired in `make_set3.py`**, not
  registered. 4 partners of 000642 are merged.

The A uncertainty is a level's, not a line's. The lines keep the head frame's u_hfs_shift, and
`levels_set3.csv` adds u_S to each level. The fit of Set 3 then reproduces iter_hfs's moved by
−S: within 0.001 cm^-1 on all 653 levels in a scratch LOPT run. The user runs the chain in
`final/`. The refit of the width model on Set 3 (the rest of Step 7) is next.

Set 3 has no displacement left between line and level, so `level_positions.py`,
`level_shifts.py` and `review_mismatches.py` can run on it unchanged. The
lines still have their widths, though: a center of gravity estimated from an
unresolved pattern is not a sharp line. The width model has to be refitted
on Set 3 (the rest of Step 7) before their verdicts mean anything. Most of
the per-level widths of `level_hfs_widths.csv` should then shrink toward zero.
(Superseded 2026-10-11, next paragraph: none of these runs on Set 3.)

*2026-10-11: no width model for Set 3; level analyses only in iter_hfs.* A width model was
refitted on final/'s LOPT residuals (per-set `[widths]` constants for `level_positions.py`) and
withdrawn the same day at the user's call. final/ is the container of the reported quantities,
for judging overall agreement and internal consistency only. Its accepted lines were moved by
their own assignments, through D and kappa, which carry large systematic errors, and the other
lines were not moved, so its residuals mix the measurement errors with the hfs correction's. No
level is searched for, and no assignment judged, in final/. `level_positions.py --set` refuses a
derived set. What the fit found is kept for the calibration's refit (README, "No level analysis
in `final/`"):
- a floor of about 0.0125 cm^-1, flat over 20000–47500 cm^-1, of unknown cause;
- bands off by up to 20 % in sigma, the same in iter_hfs against the table's own uncertainties;
- a per-line spread of kappa about its class of about 0.05·|D| (1974) and 0.26·|D| (1969).

*After the user's conference (user, 2026-10-11): refit the kappas in the calibration
accounting for the newly found trend.* The per-line spread of kappa grows with |D|, strongly
in 1969, and the floor is not understood; both are to be explained as empirical trends. The
calibration and its uncertainties are not to be changed before then. Still open, the user's
decision: whether the 1969 kappa spread, an independent per-line error, also goes into
final/'s LOPT line uncertainties. That would change the weights, and so every rerun after it.
If so, the spread is to be measured on iter_hfs's head-frame residuals, not on final/'s.

**Step 8 - regenerate the tracked outputs and the documentation.**
`level_positions.csv` and `.xlsx`, `line_classifications.*`, `README.md`.
Excel must be closed; budget about 10 minutes.

**Step 9 - the publication tables.** The Set 1 and Set 2 line lists, the A
constant table with provenance per level, the adopted uncertainties of Step 1
with the counts behind them, and a methods section stating every assumption
that survives - including whatever replaces kappa(flag) = 1 - and carrying the
reference list of section 9.

*Open item (user, 2026-10-09): the uncertainties of the Set 3 energies and Ritz
wavenumbers.* A level's cg energy carries several systematic contributions besides
LOPT's statistical one. Each is correlated in its own way, so energies and Ritz
wavenumbers do not take them the same way:

* **The ground level's S, 0.343 ± 0.021 cm^-1 (A = +0.0305(19)).** It is additive and
  common to every excited level: their energies are measured from the ground level's
  cg. It cancels in every Ritz wavenumber between excited levels, and enters a
  transition to the ground level once.
* **The scale 0.9638(79) of the semiempirical A constants.** It is common to every level
  whose A is semiempirical, in proportion to its S.
* **The radial parameters, per configuration.** Correlated within a configuration.
* **The rest of u_A per level.** u_S = I·J·u_A, independent from level to level. It is
  in `final/levels_set3.csv` as u_S, which today lumps all three A terms together.
* **LOPT's own systematic uncertainty** (D2sys, the calibration's group functions). Not
  entered yet; the user will enter it.

A Ritz wavenumber's uncertainty then needs the covariance of its two levels, not
hypot(u_E1, u_E2). This is to be worked out before the tables are made.

*Plan agreed 2026-10-10.*
1. Write the calibration groups' covariance. **Done:** `wavelength_calibration_cov.csv`;
   `make_LOPT_input.py` always writes `stat` in LOPT.par.
2. Put each calibrated line's group number in its LOPT flag, and write the group
   functions into LOPT.par (Å, power 0, `u_d_lambda_A`), for the final set only.
3. Split u_S into its correlated parts. Check the shortcut (dE_cg = kappa_bar·dS, with
   1 as the bound) against one perturbation run of the A scale.
4. Run four kappa perturbations. The user then decides on the kappa term of
   `u_hfs_shift`.
5. Assemble the totals in make_set3.

*Finding of step 1 (scratch run on iter_hfs, 2026-10-10).*
- The 151 groups are positively correlated throughout: median r +0.32, mean +0.34.
- One mode carries 23% of the variance, and it is d_lambda ∝ lambda (R² 0.93): a scale
  error, u_s = 3.1e-7, which gives u_s·E per energy (0.031 cm^-1 at 100000).
- The fitted curve holds s = +5.3e-7 of it.
- Without that mode, the correlations are moderate: median |r| 0.16, 90th percentile 0.41.
- LOPT's independent groups cannot carry this covariance. The exact way is one LOPT
  perturbation run per eigenvector of the covariance: 151 runs at 28 s, about 70 min.
- The user decided to keep the calibration's scale (2026-10-10).
- Step 3 is coded as `lopt_perturb.py --source hfs` (2026-10-10). The untested-configuration
  term gets one pattern per configuration (user's choice). Each pattern is a change of the
  model and a recomputation of every line's cg_shift, followed by LOPT.
- **Open: kappa held against the fit.** Since the composition priors were combined (commit
  c5b9ea6), the calibration fits plain_1974 = 0.542(60) and plain_1969 = 0.347(49). The
  config holds 0.751(38) and 0.440(47), kept by `--keep-kappa`. The difference is 3.5 and
  1.9 of the new uncertainties. To be settled before step 4.
- Plan step 2 is replaced by `lopt_perturb.py`: one LOPT run per eigenvector of the
  covariance, 16 runs at a time (README, "The calibration's systematic uncertainty"). LOPT's
  own groups are not used. The same driver is to carry the kappa and A-scale perturbations of
  steps 3 and 4.

*The u_A check (2026-10-09, scratch, on Sugar's resolved components; nothing changed).*
z is a difference divided by its combined uncertainty, so rms z = 1 means the
uncertainties are right.

* **Resolved components (437 positions, 121 levels).** The anchored fit has chi2/dof 0.54;
  the positions alone give chi2 0.40 per position. The spacing precision √2·0.003 Å is
  conservative, so the fitted u_A are too, by about a factor 1.6.
* **Composition priors (52 of 59 met by a pattern).** Leave-one-out: each prior freed in
  turn and compared with what the spacings give without it. rms z 0.82, so u_A is
  realistic and slightly conservative. One outlier: 000229, prior +0.0286(61) against
  spacings +0.0568(43), z = +3.8. The scale's share of u_A is about 3 %.
* **Calculated A (67 levels the spacings determine without a prior).** rms z 0.81, with no
  |z| > 3. This is a consistency check rather than an independent test, since these
  levels were among those hfs_A_theory was tested on.
* **Composition against calculated.** rms z 0.39. This is no test: both come from the same
  Cowan eigenvectors, so their errors are correlated.
* **What the spacings could add.** Combining each composition prior with the spacings
  (the anchored fit's own value for that level) would cut its median u_A from 0.0062 to
  0.0024, and u_S from 0.064 to 0.021 cm^-1. 22 of the 52 meet the table's 1/3 rule.
  **Done 2026-10-09** (`hfs_A_fit.py --combine-priors`): 51 levels were combined, with
  the priors kept in `A_hfs_priors.csv`. 000229 was held at first and then combined: its
  conflict came from the two printed components of 30879.98 *r, which the user dismissed
  as misidentified (`hfs_components.DISMISSED`). They are now lines of `Pr3_lines.xlsx`,
  with their intensities estimated from the rung shares; see the README. This belongs on
  the paper's list of contradictions.

*For the paper's discussion, not to be acted on (user, 2026-10-09):*

* **Flags whose shading disagrees with the sign of D.** 36667.32 and 29151.13, both
  `*r` with D < 0. 29151.13 also has the per-line component scale s = −6.7.
* **The 20 lines whose printed components disagree with head − D beyond 3σ**
  (`final/set3_report.txt`). For each: the suspected level(s), whether the A or the
  assignment is in doubt, and how the contradiction could be resolved.
* **The inaccuracy of the LSF**, which limits the calculated A constants and the
  predictions. The user will not refine it in this paper.
* **The 15 lines of 000190/000642 without a listed F = 2 partner.** Part of the
  general discussion of predicted strong lines missing from the list.

*Status 2026-10-11: the totals are done.* `make_set3.py --levels` writes
`final/levels_set3.csv` (u_stat = max(D1, D2stat), u_cal, u_hfs, u_kappa, D2tot*),
`final/ritz_set3.csv` and `final/levels_set3_report.txt`; see the README, "The totals".
The ground level's common part, 0.0160 cm^-1, is the ground level's own D2tot* and is
in no excited level's; the paper must say that an excited level's total uncertainty is
its D2tot* and the ground level's in quadrature (user, 2026-10-11). Before that run,
52989.1149 was rejected and 49077.8771 and 66477.4764 inflated (user, 2026-10-11): no
line of the fit is now beyond 3.5 sigma.

**Step 10 - the final line list (user, 2026-10-11; the approach is to be explained by
the user).** The set of transitions in `final/line_classifications.*` is not the one
to publish:
- it holds many transitions that lie close to a Ritz wavenumber but are too weak to be
  observable;
- it omits many potentially observable transitions that are far from any line Sugar
  observed, among them transitions masked by strong lines nearby. Some of these, not
  all, are marked as masked in `line_decisions.csv`.

The list for the database is to contain:
- every potentially observable transition, with its Ritz wavenumber and uncertainty;
- every transition Sugar identified that is now deemed unobservable, with its gA and
  calculated intensity. This is the one exception to the observability rule;
- every line Sugar observed, identified or not.

Columns still to be specified include Sugar's line characters (the list given earlier
was incomplete) and a column of coded notes: one-letter codes, comma-delimited where
several apply. For example:
- D: a transition Sugar identified with an observed line, now discarded;
- R: a reclassified line;
- N: a new identification.

The first stage is only the list of transitions. Those missing from the table have to
be added to LOPT's input, so that LOPT computes their Ritz wavenumbers and
uncertainties.

*Note on the mechanics (2026-10-11).*
- LOPT reads a record flagged `P` (predicted, lopt5.0.pl `$FPredict`) as outside the
  fit: `$FMP` excludes it from the level optimization, and its Ritz wavenumber and
  uncertainty are still computed. The table's rejected candidates already reach LOPT
  this way.
- Transitions added as `P` records therefore leave every level energy unchanged. The
  per-level shifts of `lopt_perturb.py` are unchanged too, so the two passes need not
  be rerun: `make_set3.py --levels` takes each Ritz wavenumber's systematic parts from
  the modes files, for any pair of levels, and its `SYS_TOL` check of the energies
  confirms that the fit is the same.
- LOPT needs a value in the wavelength field of a `P` record. Its Ritz output does not
  depend on that value.

---

## 6. Housekeeping the hfs work has accumulated

These are small and unrelated to the physics, but they must not be lost: they
belong to the uncertainty model the hfs work sits on.

* Second calibration pass: the character constants in `level_positions.py`
  were fitted against the *old* fit's residuals. They must be iterated to
  convergence against the corrected fit.
* Repair stale text in `level_positions.py`: the "(312 lines)" comment (371
  now); the §4 `w` evidence sentence's band numbers and its "12 units of ln L";
  the `('1969','w')` comment's false claim that one line survives the Ritz cut
  (13 do). Same for `README.md:3086`.
* Decide `CHAR_DLAM[('1969','cl')] = 0.0044` Å against the 87 lines that fit
  0.0012–0.0015 Å (upper limit about 0.0023).
* Cap or refit the `c` width above 38356 cm^-1 instead of extrapolating the
  lambda-squared law further.
* `Unsolved_problems.md` entry 1, lines 143–149, claims a `reject` row in
  `line_decisions.csv` for `(10682.6425, 059003.000328, 059003.000369)` that
  does not exist. Fix the text.
* Decide whether `HFS_APPLY` should be lowered — the held-out evidence is
  −58.2, i.e. against.
* After any further change to `unc_own`: rewrite it in `Pr3_lines.xlsx` and
  rewrite `dlv.dat` (the script is in the scratchpad, about one minute).
* Keep-or-roll-back decision on the uncertainty-model trial as a whole, after
  the user's `check_sync.py` inspection.

---

## 7. Files

### Existing, additive; nothing else modified

| file | contents |
|---|---|
| `A_hfs_levels.csv` | `level_id, cfg, J, A_cm-1, u_A, source, n_flagged` — 129 levels |
| `hfs_line_corrections.csv` | `wn_obs, char, class, low_id, upp_id, S_low, S_upp, kappa, delta_cm-1, u_delta, wn_cg` — 3834 lines |
| `hfs_level_shifts.csv` | `level_id, J, A_cm-1, u_A, A_source, E_setI, E_setII, shift` — 535 levels |
| `level_hfs_widths.csv` | the per-level hfs width used by the uncertainty model |
| `levels_pub.xlsx` | the calculated level compositions (an input, not generated here) |

The column names `E_setI` / `E_setII` in `hfs_level_shifts.csv` predate the
Set 1 / Set 2 naming and are to be renamed `E_set1` / `E_set2` when that file
is next regenerated.

### To be created

| file | purpose |
|---|---|
| `hfs_components.csv` | **written** - the 500 measured hfs component positions of section 3.5, joined to their parent lines; a tracked input, extracted from the 1974 workbook |
| `hfs_components.py` | **written** - the extraction, the pairing and the side check |
| `hfs_patterns.py` | **written** - the interval rule, the 6-j intensities, the scale tests and the global fit of section 3.5 |
| `hfs_A_from_patterns.csv` | **written** - the A constants the measured patterns give, with the caveat on their absolute scale; not an input to anything yet |
| `hfs_kappa.py` | **written** - the adopted uncertainties of section 2.6, the convention fit of section 2.2 and the two tests of section 3.3 |
| `hfs_uncertainties.csv` | **written** - the adopted uncertainty per character per era with the count behind it; the Step 1 deliverable |
| `hfs_misassigned_candidates.csv` | **written** - the lines the trend fit rejects beyond 4 sigma, for review in IDEN2; no ledger row is written from it |
| `tests/test_hfs_kappa.py` | **written** - the degeneracy theorem, the recovery of a kappa put in by hand, and the conclusions of Steps 1 and 3 as ranges |
| `wavelength_calibration.py` | **written** - the plate wavelength calibration of Step 3b, fitted jointly with the level energies and the kappas |
| `wavelength_calibration.csv` | **written** - the fitted curve, one row per group, with the correction in cm^-1 and its systematic uncertainty |
| `wavelength_calibration_points.csv` | **written** - the input of that fit, one row per line: the wavelength, the displacement it asks for, its uncertainty, its effective uncertainty and its leverage |
| `calibration_breaks.txt` | **written** - the breaks inside a plate, one wavelength per line with the evidence for it; read by `wavelength_calibration.py` |
| `inflated_unc_hfs.csv` | **written** - the hfs attribution of each line of `inflated_unc_lines.txt` (section 3.7); a 2026-09-29 snapshot, not an input to anything |
| `hfs_theory.py` | the composition-to-A calculation, moved out of the scratchpad, with the LS extension |
| `hfs_A_fit.py` | **written** 2026-10-04 - the anchored fit of Step 4: A constants from Sugar's resolved components with the confirmed semiempirical constants as priors, its checks, and `--write-levels`; tests in `tests/test_hfs_A_fit.py` |
| `hfs_A_anchored.csv` | **written** 2026-10-04 - every level the Step 4 fit meets, with its value, the table's value before, its flag bounds, its value from its own lines, and whether it was adopted |
| `hfs_correction.py` | **written** - the head-frame model of Step 5: S = I*A*J per level, kappa per class, the shift of a line and its uncertainty, and the record `LOPT_hfs_shifts.txt` |
| `tests/test_hfs_theory.py` | the 22-level reproduction test and the parser tests |
| `tests/test_hfs_correction.py` | **written** - the configuration, the groups of levels, the arithmetic of a shift, the allowances D13 withdraws, the LOPT input and the readers downstream |
| `tests/test_hfs_components.py` | **written** - the extraction of the 501 rows, the pairing to flagged lines, the 500-of-500 side agreement, and every number section 3.5 quotes |

### Modified for Step 5 (2026-09-30)

`classify_lines.py`, `make_LOPT_input.py`, `models.py`, `config.py`,
`lineclass_config.toml`, `sync_IDEN2.py`, `check_sync.py`, `README.md`.
Still to be made hfs-aware: `level_positions.py`, `review_mismatches.py`,
`insert_new_level.py`, `move_level.py`, `discard_level.py`.

### Added for the per-line kappa exceptions (2026-10-05)

`kappa_exceptions.txt` (new, hand-kept, `files.kappa_exceptions`) and
`tests/test_kappa_exceptions.py` (new). Modified: `hfs_correction.py`,
`classify_lines.py`, `models.py`, `hfs_kappa.py` (`read_lines(model=)`),
`wavelength_calibration.py` (`--config`, the set's hfs model), `config.py`,
`lineclass_config.toml`, `A_hfs_levels.csv` (000190), `README.md`.

### The kappas written by the calibration (2026-10-05)

`wavelength_calibration.py` now writes the kappas it fits into
`[hfs.kappa]` itself. It writes the file of the configuration chain that
holds the table, today the baseline's (`config.kappa_file`,
`config.set_hfs_kappa`), and `--keep-kappa` opts out. So D12, "change them
together with the calibration", no longer needs a hand step.
`tests/test_set_hfs_kappa.py` is new. The pin of plain_1969 in
`tests/test_hfs_correction.py` became a range check.

### The search for kappa exceptions (2026-10-05)

The kappas converged at plain_1974 0.817(13), plain_1969 0.632(20) and
c 0.174(50).

`kappa_candidates.py` is new, with `tests/test_kappa_candidates.py`. It
is a read-only search of a working set for registry candidates. For each
line it compares three hypotheses for where the line was measured: its
class's kappa, the head (kappa = 1) and the center of gravity (kappa = 0).
The residual tested is the one the line would have without itself in the
fit, taken from LOPT's leverage (uWnCstat/uWnOStat)^2 (the user's choice
over counting lines per level). The uncertainty tested is the
calibration's model value, never the registry's. The method is described
in the README.

First run on iter_hfs, after convergence:
- 6 candidates at the cg, 3 at the head, and 15 lines whose class fails
  with neither alternative fitting;
- all 18 held rows pass their audit. 44494.0809 (2.6 sigma) and the
  satellite head 23792.4971 (2.3 sigma) failed the first rule, which
  cut at 2 sigma. On the user's decision an audit now fails only
  beyond 3 sigma, as a class does, and is made at the held kappa
  itself (z_now);
- levels 000150 (1129) and 000165 (1116) carry most of the failures, in
  both directions.

Next step: the user reviews the candidates in IDEN2 and records the
accepted rows. The calibration is rerun only after that.

### Partial patterns (2026-10-05)

A registry class `partial:<first>-<last>` holds a line at the center of
gravity of rungs first to last of its own pattern. The κ comes from the two
levels' J and A (`hfs_patterns.partial_pattern`, `Model.row_kappa`) and
follows them. `kappa_candidates.py` reports the best such part for every
line that does not fit where it is, as information only, together with the
evidence: a line, or a candidate transition, where the left-out part would
lie. User's rule: a left-out part must be accounted for by a blending
transition; otherwise the line stays at the cg with an inflated
uncertainty. On the first run 8 failing lines fit a part and none had
evidence. Tests: `tests/test_partial_kappa.py`.

### Reader and Sugar's A constants replaced (2026-10-05)

Terminology, on the user's correction: Reader and Sugar's (1965) values are
semiempirical, not calculated - only the angular factors theta are
calculated; a_4f and a_6s are fitted to measured hyperfine structure. The
word is changed throughout (`hfs_kappa.is_semiempirical`, the source label
"Reader & Sugar 1965 (semiempirical)").

Fitted from Sugar's resolved components with their own priors removed, the
20 levels come out at 0.76 - 1.10 times the published values; the mean,
0.968 +- 0.005, reproduces A_SCALE, but the scatter is level by level
(000149, 000150, 000145, 000142 low by 2.3 - 4.6 sigma; 000165, 000191,
000140 high by 2.5 - 3.5 sigma). 000149 and 000145 are the levels whose
4f2.5d admixture the intensity analysis of 2026-09-27 found larger than
Cowan's. The user had the component values adopted: the 20 are free levels
of `hfs_A_fit.py` (`FREED_SOURCES`), only the 59 composition values remain
priors, and `--write-levels` replaced 66 rows (the 20, and the 46 earlier
fitted values refreshed; 13 moved by more than 0.0005, at most +0.0079 for
000656). For 000165 the new A, -0.0431 +- 0.0015, does not explain its
off-class lines (30074.3893 would need about -0.064). A refit of Reader and
Sugar's parameters with the present compositions waits for the redone
level-structure fit. Next: recalibration and the iter_hfs chain (the user).

### Head rows for wide patterns tried and withdrawn (2026-10-06)

The unflagged 1974 lines of 000142, 000148 and 000182 (4f2.6s, the widest
patterns) sit near their heads, and with kappa held at 0.811 a scratch LOPT
fit favored head rows for them (chi2 -95; their flagged lines came in). The
9 rows were written, and the next calibration refitted plain_1974 from 0.811
to 0.697 (plain_1969 0.628 -> 0.549, c 0.157 -> 0.100). Eight lines IDEN2
accepts were rejected, among them 24792.516 (000035-000145, measured at kappa
0.85). Scratch calibrations on one table gave 0.810 with the old registry and
0.697 with the new one: seven lines move the class by 0.11. Lines of narrow
patterns carry next to no information on kappa, because the level energies
absorb it. The class value rests on a few dozen wide-pattern lines on levels
that have flagged lines too, and that population is mixed. Registering its
head-measured members as exceptions pulls kappa down, and the next members
fail: a runaway. The rows were removed and the three inflations restored
the same day, and the recalibration gave plain_1974 0.809(13).

### The flags split (2026-10-06)

kappa(flag) = 1 is the anchor (the degeneracy theorem, section 2.1), and its
evidence is the printed components, which exist for 248 of Sugar's 296
flagged lines. On the user's question the other 48 were tested. The 38
of them in the fit sit 0.043 +- 0.010 cm^-1 from the head towards the
center of gravity (kappa 0.92 +- 0.02, against 1.01 +- 0.01 for those with
printed components and 0.81 for plain_1974). A main-line offset fitted
with the A constants on the printed components is +0.010 +- 0.006 cm^-1,
independent of the first rung's separation. So an unlisted flagged pattern
was not resolved, and Sugar set on a peak pulled towards the center of
gravity, though not as far as on a plain line.

Implemented: `files.hfs_components` (the baseline names hfs_components.csv)
splits the flags; a line within 0.006 cm^-1 of a `wn_line` there is `flag`,
the others `flag_unlisted`, a fitted class of `[hfs.kappa]`
(`hfs_kappa.kappa_class`, `hfs_correction.Model.line_class`,
`wavelength_calibration.kappa_name`, `config.FITTED_KAPPA_CLASSES`).
Measured in scratch on iter_hfs/'s table of 10:21: plain_1974 0.783 with
and without the split, flag_unlisted 0.895 +- 0.024, written into the
configuration as the first value. The same scratch run showed plain_1974
following the accepted set it is fitted on: 0.809 from the 09:14 table,
0.783 from the 10:21 one. Next: recalibration and the iter_hfs chain (the
user).

---

## 8. Standing rules this work obeys

* No stage of the assignment pipeline (`classify_lines.py`,
  `make_LOPT_input.py`, `sync_IDEN2.py`, `insert_new_level.py`) is run without
  the user asking for that specific run. Changes are proposed, not applied.
* Every uncertainty that appears in a publication must be reproducible by a
  reader from the published tabulation. No hidden inflation factors.
* LOPT's own reported uncertainties are model-dependent and are never
  presented as an absolute sigma scale.
* All text files are written with LF line endings only.
* Every changed assignment is reviewed by hand before it is accepted.
* Every external datum is cited to a source the user can check, in section 9,
  so that the paper's reference list is assembled as the work goes rather than
  reconstructed at the end. A number whose source cannot be produced is not
  used.
* Sugar's own stated uncertainties are the starting point for anything fitted
  to residuals. The Set 1 inflation of `unc_own` is Set 1 machinery and is
  never used as a weight in the derivation of hfs trends.

---

## 9. References

Every external datum in this document traces to one of these, and the paper
will carry the same list. Where a reference is not yet verified, it says so.

* **[R1] J. Reader and J. Sugar**, "Nuclear Magnetic Moment of Pr-141 from the
  Hyperfine Structure of Doubly Ionized Praseodymium", *Phys. Rev.* **137**,
  B784-B789 (1965) (Volume 137, Number 4B, 22 February 1965). Local copy:
  `LineClass/Reader_el_3278_1965.pdf`; the volume, issue, date and page range
  were read from that copy. Source of a_4f = 0.029 +- 0.003 and
  a_6s = 0.639 +- 0.002 cm^-1, of the angular coefficients θ₁ and θ₂ of its
  Table III from which the 22 semiempirical 4f²6s A constants in
  `A_hfs_levels.csv` are evaluated here, and of the total hfs widths in the
  `hfs_width` sheet of the 1974 extraction workbook. **Reader and Sugar
  publish no per-level A constant**, and their measurements excluded levels
  with J < 7/2; see section 2.3.
* **[R2] J. Sugar**, *Phys. Rev. Lett.* **14**, 731-732 (1965) - the Pr IV
  analysis; source of the quotation in section 1.1 about measuring the centers
  of unresolved hfs blends. Citation as supplied by the user; to be checked
  against the original before publication.
* **[R3]** = [R2]; section 1.1 cites it under that label.
* **[R4] J. Sugar**, "The Third Spectrum of Praseodymium (Pr III) in the 
  Vacuum Ultraviolet", *J. Res. Natl. Bur. Stand. (U.S.), Sect. A* 73, 333–381 
	(1969) - the 1969 Pr III line list. Source of the 0.004 angstrom accuracy 
	statement.
* **[R5] J. Sugar**, "Spectrum of Doubly Ionized Praseodymium from 2107 Å to 
  10716 Å", *J. Res. Natl. Bur. Stand. (U.S.), Sect. A" **78**, 555–593 
	(1974); DOI: 10.6028/jres.078A.036. Source of the 0.003 angstrom
  plain-line uncertainty, the 0.007 angstrom average deviation from Ritz, the
  `*r`/`*v` flags, and the 501 component wavelengths of section 3.5.
* **[R6] Cs I and Ba II 6s and 6p hfs constants** (0.127/0.022 and
  0.185/0.028 cm^-1) used in section 2.4 to adopt a_6p(1/2) and a_6p(3/2).
  **Not traceable to a verified source at present.** Either a citation is
  found and checked, or the adopted values are refitted from Pr III data alone
  and the analogy is dropped from the paper.
* **[R7] A. Kramida**, *Eur. Phys. J. D* **78**, 36 (2024) - the method by
  which the DREAM transition-rate uncertainties were determined (cited in
  `README.md`).
* **[R8] A. Kramida, A. N. Ryabtsev and P. R. Young**, *Astrophys. J. Suppl.
  Ser.* **258**, 37 (2022) - the reduction of Sugar's plate intensities to a
  common linear scale (cited in `README.md`).
