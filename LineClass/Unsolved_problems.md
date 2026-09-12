# Unsolved problems

A running record of cases where the analysis of the Pr III line list reaches a
dead end: the evidence is complete, the alternatives have been enumerated, and
none of them accounts for everything. These are kept deliberately, as examples
of the residual uncertainty that remains in a level optimization of this kind,
and are intended for inclusion in the article describing this project.

Each entry states what is observed, what the level network says, what was tried,
and what would be needed to settle it.

---

## 1. The line at 10682.643 cm-1 - a perfect identification that the level network refuses

### The observation

`Pr3_lines.xlsx` carries an observed line at **10682.6425 cm^-1**, measured
uncertainty **0.0034 cm^-1**, from Sugar (1974). Its recorded plate intensity is
`Iorig = 25`, the faintest step of the visual estimation scale used on those
photographic plates. In the pipeline it is currently **unassigned**
(`line_classifications.csv`, no `low_id`/`upp_id`).

### Why it looks like a solved case

*It is the only possible identification.* A search of all 595 established
(observed) levels of `IDEN2/enlev.dat` for a pair whose separation falls within
+/-0.25 cm^-1 of the observed wavenumber, with |dJ| <= 1 and opposite parity (the
electric-dipole selection rules), returns exactly one candidate - and no second
candidate anywhere near:

| IDEN2 row | level id | E (cm^-1) | J | dominant configuration | parity |
|---|---|---|---|---|---|
| 878 | `059003.000328` | 102 981.887 | 6.5 | 4f2.6d `_3H4I` | even |
| 778 | `059003.000369` | 113 664.597 | 5.5 | 4f2.5f `~3H4I` | odd |

The transition 878 - 778 is electric-dipole allowed. Its Ritz wavenumber - the
difference of the two adopted level energies - is **10682.7100 cm^-1**.

*The predicted intensity agrees.* `Icalc.xlsx` gives this transition
`gA = 1.66e8 s^-1` and a predicted intensity `Icalc = 441`, with a logarithmic
uncertainty `u_ln = 0.40` - a one-standard-deviation band of about a factor 1.5
either way. The observed calibrated intensity is `Icor = 725`, a factor 1.6
above prediction, i.e. agreement at 1.3 sigma.

### Why it cannot be accepted

**The wavenumber disagrees by far more than the data allow.**
O - C = 10682.6425 - 10682.7100 = **-0.0675 cm^-1**. Three independent ways of
scaling that:

* against the line's own uncertainty (0.0034 cm^-1): **-20 sigma**;
* against the uncertainty of the level separation itself. Forming the normal
  matrix of the full weighted least-squares problem over all 5003 fitted lines
  and solving for the variance of the single combination E(778) - E(878) gives
  **+/-0.0083 cm^-1**; combined with the line's own error the discrepancy is
  **-7.5 sigma**;
* against its peers. Among the 130 accepted, full-weight lines of the same
  precision class (uncertainty <= 0.006 cm^-1, all from the same source) the
  root-mean-square of residual/uncertainty is **1.10** - the quoted
  uncertainties are honest, not optimistic - and the single largest residual
  anywhere in that set is **0.0166 cm^-1**. The disputed line would be four
  times worse than the worst of its 130 peers.

**It is not one weak link.** Independent Ritz paths from 878 to 778 that avoid
the disputed line, each progressively banning the weakest link of the previous
one:

| route (IDEN2 rows) | E(778) - E(878) | +/- |
|---|---|---|
| 469 - 331 - 470 - 879 | 10682.741 | 0.013 |
| 471 - 324 - 469 - 331 - 470 - 879 | 10682.678 | 0.018 |
| 379 - 324 - 469 - 331 - 470 - 879 | 10682.664 | 0.019 |
| 471 - 420 - 467 - 421 - 470 - 879 | 10682.721 | 0.025 |
| 379 - 609 - 375 - 595 - 367 - 350 (no 879, no infrared link) | 10682.770 | 0.039 |

Five routes through different levels, one of which never touches level 879 or
the near-infrared region at all, all land between 10682.66 and 10682.77. None
reaches 10682.643.

**A calibration offset cannot explain it.** The disputed line and every
constraint it fights - 10800.473, 10834.136, 11064.270, 11749.721, 11831.808,
12112.214, 12236.800, 12345.696, 12426.524, 12499.263 - come from the same
source and the same spectral region, and its immediate neighbours 10800.473 and
10834.136 fit the network to 0.0006 and 0.0050 cm^-1 respectively.

**Forcing it in spreads damage rather than absorbing it.** A trial refit of the
whole network with the line added at full weight raises chi^2 by **+54**, and
the line still ends at +2.8 sigma. Two rigid, disjoint blocks of levels move
against each other - 878, 523, 471, 332, 379 all rise by about +0.03 cm^-1 while
778 and 320 fall by about -0.02 cm^-1 - and seven previously good lines degrade,
notably 12112.214 (-1.7 -> -4.7 sigma), 10834.136 (-1.2 -> -3.1 sigma) and
12345.696 (+1.1 -> +2.9 sigma). By the standing rollback rule - one weak
casualty may be paid for a strong branch, but widespread residual disturbance
means roll back - the trial is rejected.

### Why the intensity agreement carries less weight than it appears to

`Icor` for this line is not an independent measurement. It is the plate estimate
`Iorig = 25` - the faintest available step - put through the smooth,
wavenumber-dependent intensity calibration, so it is a deterministic function of
that step and the wavenumber. Every `Iorig = 25` feature in this window comes
out with the same value:

```
10606.75  Icor 680     10665.32  Icor 715     10682.64  Icor 725  <- the line
10653.21  Icor 708     10676.11  Icor 721     10686.57  Icor 727
10657.46  Icor 710                            10729.41  Icor 751
```

So "the intensity matches perfectly" would be equally true of any of the six
other unidentified threshold-level features in the same window. Worse for the
assignment, `Icalc = 441` lies **below** the faintest feature actually recorded
in the region (the plate floor is at `Icor` about 680): the transition is
predicted at or under the detection limit, so there is no expectation that it
was recorded at all.

The chance arithmetic is consistent with a coincidence. Between 10400 and
11000 cm^-1 there are 35 recorded lines, 19 of them unidentified - a density of
0.032 unidentified lines per cm^-1. The probability that one falls within
+/-0.0675 cm^-1 of an arbitrary Ritz value is about 0.4 %. Small for a single
transition, but thousands of predicted transitions are scanned, so coincidences
at this level are expected to occur, and this is what one looks like.

### The readings that survive

1. **Most likely:** the 878 - 778 transition is below the detection limit of the
   available plates and was not recorded. The feature at 10682.643 is a
   different, still unidentified line - a Pr III transition involving a level
   not yet found, a line of another species, or a plate artefact - that happens
   to lie 0.067 cm^-1 from the Ritz value.
2. **A blend.** If 878 - 778 contributes part of the light and a brighter
   component dominates, the measured centroid is set by the other component. The
   intensity agreement would then also be accidental, and the wavenumber still
   must not enter the fit. Note that no width flag (`Ch. = w`) is recorded for
   10682.643, although one is recorded for 10789.69 and 12112.214 in the same
   region, so there is no measurement note supporting this.
3. **A wavelength blunder** in the original measurement of this one faint line.
   0.0675 cm^-1 at 10682 cm^-1 is 0.059 A at 9361 A - not impossible on a
   photographic plate, but no level may be moved on that hypothesis without
   independent evidence, and the outcome for the fit is the same.

### Current disposition

The line is left unassigned and out of the LOPT input. A `reject` row in
`line_decisions.csv` for `(10682.6425, 059003.000328, 059003.000369)` records
the reason so the candidate does not resurface on the next classification run.
In a published line list the honest note is "possibly 878 - 778,
O - C = -0.068 cm^-1, not used".

### What would settle it

**A new recording of the near-infrared region at higher exposure.** Everything
hangs on features at or below the plate threshold: the predicted transition is
fainter than anything recorded there, and the feature that was recorded carries
the coarsest possible intensity estimate. A deeper exposure of the region around
9360 A (10600-10800 cm^-1), with a modern detector and wavenumber calibration,
would decide the case directly:

* if a line appears at 10682.71 +/- 0.01 cm^-1, the assignment is correct and
  10682.643 is a separate, unrelated feature;
* if 10682.643 is confirmed at its measured position with no companion at the
  Ritz value, the transition is genuinely absent and the feature belongs to
  something else - and, being then well measured, it becomes a useful constraint
  on whatever level is still missing;
* if the feature at 10682.643 is not reproduced at all, the original measurement
  carried a blunder and the case closes.

The same argument applies to the other threshold-level features listed above,
several of which are unidentified for the same reason. A higher-exposure
near-infrared survey is the single measurement most likely to resolve this class
of problem.

### Secondary lead

Level 878's weakest precise anchor is the line at **12112.214 cm^-1**
(878 - 469): it is flagged wide (`Ch. = w`), has the largest uncertainty of
878's four infrared constraints (0.0103 cm^-1) and already the largest residual
(+0.017 cm^-1, 1.7 sigma). It is the only constraint whose re-examination could
move 878 downward. However, moving 878 by 0.067 cm^-1 breaks 33052.338,
12426.524 and 12236.800 - the last of which currently fits to 0.0003 cm^-1 - so
this is not expected to survive scrutiny.
