# Pr III Line Classification Pipeline (`LineClass`)

A fully-automated, deterministic pipeline for classifying observed spectral lines of **Pr III** against every theoretically allowed transition between known energy levels.

This code was designed to solve a problem specific to the current state of knowledge on the Pr III spectrum:
- Two separately published line lists represent a complete set of fairly precisely measured lines (Sugar 1969, 1974), but they are **incomplete** in the sense that many lines are unclassified.
- J.-F. Wyart (1999) constructed an unpublished set of **adopted energy levels** for Pr III, but he never published a complete set of observed lines defining each of these levels. 
- An extensive list of calculated transition rates exists in a public database DREAM. This list was derived with the use of most but not all of Wyart's energy levels and covers only a limited range of wavelengths.

The main purpose of this code is to find all observed lines defining each known energy level. Prior to running it, the input data were prepared in the following way:
- Each observed wavenumber has been assigned an estimated uncertainty based on original Sugar's general statements and a number of other criteria.
- A complete list of all previously known energy levels of Pr III has been constructed, including the data from NIST ASD and Wyart's new levels communicated by him to the DREAM team.
- Uncertainties of the DREAM transition rates for Pr III have been determined by comparison with other available data using the method described by A. Kramida,
Eur. Phys. J. D 78, 36 (2024).
- The observed line intensities reported by Sugar (1969, 1974) have been reduced to a common linear scale by using the method described by A. Kramida, A. N. Ryabtsev, and P. R. Young,
Astrophys. J., Suppl. Ser. 258, 37 (2022). This procedure simultaneously produced a set of calculated intensities corresponding to the same intensity scale. These calculated intensities are based on a simplified Boltzmann model of level populations with effective temperature of 1.6 eV and the transition probabilities from DREAM. 

Given these data for observed lines and energy levels, the pipeline performs the following actions:
- enumerates all electric-dipole (E1) transitions in range;
- matches them to observed lines;
- grades each match on wavenumber and intensity agreement;
- resolves conflicts;
- iteratively "weeds" out spurious classifications using per-level intensity calibration;
- approximately refines the level energies by a simplified weighted least squares method. 

The output is a classification table ready to be fed to the **LOPT** level-optimization program.

The classification is accompanied by a **statistical validation suite** (`decoy_mc.py`, `level_shifts.py`, `level_interchange.py`, with shared utilities in `chance_mc.py`) that measures how often the pipeline would confirm a level that is *not* real, assigns to every tested level a probability of being spurious, and looks for pairs of levels whose theoretical identities may have been interchanged. See **Validation of the classification** below.

This module is the deterministic successor to the LLM-based post-processing that previously followed the [`TableExtraction`](../TableExtraction/README.md) OCR pipeline. It shares one input file (the Wyart 1999 energy levels) with `TableExtraction`.

---

## Prerequisites

- **Python 3.10+** — the code uses PEP 604 union annotations (`dict | None`) and builtin-generic hints (`dict[str, EnergyLevel]`, `tuple[...]`).
- Python packages (all already listed in `../TableExtraction/requirements.txt`):

  | Package    | Purpose                                                     |
  |------------|-------------------------------------------------------------|
  | `numpy`    | Numerical work in weeding & optimization                    |
  | `pandas`   | Building and writing the output table                       |
  | `openpyxl` | Reading `.xlsm`/`.xlsx` inputs and writing the Excel output |

- **Bundled `statistics.py`** — the repository root contains a small `statistics.py` that provides `mandel_paule()` (Paule–Mandel consensus weighted mean, NBS 1982). It is imported (inside `_wm_mandel_paule`) **only** when the per-level intensity-factor machinery is enabled with `FACTOR_WEIGHTING = 'mandel_paule'`; with the current defaults (`USE_INTENSITY_ADJUSTMENT = 0`) the module is never imported and no `PYTHONPATH` setup is needed. If you re-enable that option, the repository root must be importable (see *Running* below), because the local module **shadows the standard-library `statistics`**, which does not provide `mandel_paule`.

> **Note:** `LineClass` uses no LLM, no network access, and no API key. It is purely numerical.

---

## Setup

`LineClass` has no separate `requirements.txt`; it reuses the environment created for `TableExtraction`.

```bash
# From the repository root
python3 -m venv venv
source venv/bin/activate            # macOS/Linux
# venv\Scripts\activate             # Windows
pip install -r TableExtraction/requirements.txt
```

---

## Running

The script imports two modules that live in **two different directories**:

- `models` — in `LineClass/` (the script's own directory)
- `statistics` — the bundled `statistics.py` in the **repository root**

Both directories must therefore be on `sys.path`. The script's own directory is always added automatically, but the repository root is not, so add it via `PYTHONPATH`:

```bash
# From the repository root:
PYTHONPATH=. python LineClass/classify_lines.py

# …or from inside LineClass:
cd LineClass
PYTHONPATH=.. python classify_lines.py
```

> **Why `PYTHONPATH`?** Without the repository root on the path, `from statistics import mandel_paule` resolves to the standard-library `statistics` (which has no `mandel_paule`) and the run fails during weeding. Running the file directly with `python LineClass/classify_lines.py` — with no `PYTHONPATH` — will fail for this reason.

### Running from PyCharm

PyCharm puts the repository root on `sys.path` for you, so no `PYTHONPATH` is needed — provided the root directory is registered as a content root or a sources root. Either of the following is sufficient:

1. **Open the project at the repository root.** When the repo root is the project's content root, PyCharm's run configurations add it to `PYTHONPATH` automatically (the **Add content roots to PYTHONPATH** option, enabled by default). The bundled `statistics.py` then correctly shadows the standard library.

2. **Or mark the repository root as a Sources Root.** In the *Project* tool window, right-click the repository-root folder → **Mark Directory as → Sources Root** (it turns blue). This adds it to `PYTHONPATH` via the **Add source roots to PYTHONPATH** option.

To confirm the settings for your run configuration: **Run → Edit Configurations… →** select the `classify_lines` configuration, and under the environment/options section verify that **Add content roots to PYTHONPATH** and **Add source roots to PYTHONPATH** are both checked (both are on by default). `models.py` resolves regardless, since PyCharm runs the script with its own directory (`LineClass/`) on the path.

> If you ever see `ImportError: cannot import name 'mandel_paule' from 'statistics'`, the repository root is not on the path — apply one of the two options above (or use the `PYTHONPATH=.` command line shown earlier).

`classify_lines.py` has no command-line flags. Behavior is controlled by module-level constants and function defaults (see *Configuration*). `main(max_cycles=20)` runs the full pipeline and writes the output files; the optional keyword arguments `wn_shift`, `decoy_shift`, `decoy_ids`, and `write_files` exist for the validation suite (see *Validation of the classification*).

The validation scripts are run from the same directory:

```bash
python decoy_mc.py             # decoy runs (add --smoke for a quick plumbing test)
python level_shifts.py         # the validation report and level_shift_report.csv/.xlsx
python level_shifts.py --detail 059003.000483   # inspect one level in depth
python level_interchange.py    # interchanged theoretical identities
python swap_line_assignments.py id1 id2        # repair one: all three files below at once
python swap_line_assignments_LOPT.py id1 id2   #   the LOPT transitions file
python swap_line_assignments_IDEN.py id1 id2   #   the IDEN2 working files
python swap_line_assignments_pipeline.py id1 id2  #   the overlay the pipeline reads
python chance_mc.py            # optional: shifted-wavenumber cross-check runs
```

---

## Input Files

All paths are resolved relative to the script's own directory (`SCRIPT_DIR`).

| File                                | Location (constant)                   | Sheet       | Description                                                               |
|-------------------------------------|---------------------------------------|-------------|---------------------------------------------------------------------------|
| `Pr3_lev_Wyart_1999.xlsm`           | `../TableExtraction/` (`LEVELS_FILE`) | `Wyart2000` | Energy levels (Wyart 1999 adopted values). Shared with `TableExtraction`. |
| `Icalc.xlsx`                        | `LineClass/` (`files.icalc`)          | `Icalc`     | Calculated (theoretical) transition intensities and their uncertainties.  |
| `Pr3_lines.xlsx`                    | `LineClass/` (`LINES_FILE`)           | `Sheet1`    | Observed spectral lines with existing (Sugar 1969/1974) classifications.  |
| `revised_level_energies.csv`        | `LineClass/` (`files.level_overrides`, optional) | csv | Adopted energies revised by the identification work, as `level_id,E_input` plus a comment column. A level listed here is read at that energy instead of the workbook value, and the comment is read too — it is what says whether the level was re-positioned or exchanged with another, and hence what becomes of the identifications the published line list makes with it. See [Revised level energies](#revised-level-energies). |
| `line_decisions.csv`                | `LineClass/` (`files.line_decisions`, optional) | csv | The verdicts reached by hand on individual assignments, as `wn_obs,low_id,upp_id,decision` plus free `date`/`reason` columns. Applied after every automatic step, so it has the last word; an `accept` row is an order and *creates* its assignment when the matching never proposes it, and a row that cannot be carried out stops the run. See [The decision ledger](#the-decision-ledger-identifications-ruled-on-by-hand). |
| `intensity_correction_functions.txt`| `LineClass/` (`CALIB_FILE`, optional) | text        | Validation only: piecewise polynomials `P(λ_vac)` per wavelength region (`λ_start λ_end c0;c1;…;cn`, ascending powers, λ in Å). Sugar's plate intensity converts to the linear scale as `I_linear = 1000·I_Sugar·exp(P(λ))` (the factor 1000 makes the linearized intensities in `Pr3_lines.xlsx` integers, the smallest being 21); used by `level_shifts.py` to drop predictions below Sugar's noise level, which is `1000·exp(P(λ))` on the linear scale. |

### Column mapping (as read by the code)

The file names, worksheet names and column names below are the ones configured in
`lineclass_config.toml`; the readers locate each column by its **name** in the header
row, so the `#` position columns are informative only. See [Configuration](#configuration).

**Energy levels** — `read_energy_levels()`, sheet `Wyart2000` (header row skipped):

| Column    | #  | Field             | Notes                                                         |
|-----------|----|-------------------|---------------------------------------------------------------|
| `J_adopt` | 11 | `J_str` → `J_val` | Parsed by `parse_J` (`'5/2'` → `2.5`)                         |
| `E_ASD`   | 12 | `is_new`          | Blank ⇒ `is_new = 1`: level absent from the ASD compilation   |
| `E_adopt` | 13 | `energy`          | cm⁻¹                                                          |
| `par_ASD` | 14 | `parity`          | `e` / `o`                                                     |
| `ASD_id`  | 23 | `level_id`        | Unique level identifier (string)                              |

Rows missing energy, J, or level id are skipped.

**Calculated transitions** — `read_transitions()`, sheet `Sheet1`:

| Column    | # | Field            | Notes                                                |
|-----------|---|------------------|------------------------------------------------------|
| `id1_W99` | 1 | lower level id   | Must exist in the levels dictionary                  |
| `id2_W99` | 2 | upper level id   | Must exist in the levels dictionary                  |
| `u%gA`    | 5 | `u_calc`         | Stored as `ln(u%gA/100 + 1)` (log-scale uncertainty) |
| `Icalc`   | 6 | `calc_intensity` | Theoretical intensity, same scale as observed        |

Indexed by `(lower_id, upper_id)`.

**Observed lines** — `read_observed_lines()`, sheet `Sheet1`:

| Column    | # | Field                   | Notes                                                |
|-----------|---|-------------------------|------------------------------------------------------|
| `own`     | 1 | `wavenumber`            | cm⁻¹                                                 |
| `unc_own` | 2 | `wn_uncertainty`        | cm⁻¹ (σ)                                             |
| `Icor`    | 3 | `intensity`             | Observed intensity, reduced to a uniform scale (> 0) |
| `Ch.`     | 4 | `line_character`        | `h`, `w`, `bl`, …                                    |
| `id1`     | 5 | original lower level id | Original Sugar classification                        |
| `id2`     | 6 | original upper level id | Original Sugar classification                        |

Column 17, `Iorig`, is not read by the pipeline: it holds the intensity Sugar printed for
the line, from which `Icor` is built, and is used only by the calibration tools (see
[Sugar's own intensities](#sugars-own-intensities-toolsattach_sugar_intensitiespy)). Every
line has one, including the 23 whose `ref` is `Wyart1974`: that column names whose
classification the row carries, not who measured the line, and all the lines in the list
were measured by Sugar.

Consecutive rows with the **same wavenumber** are folded into a single `SpectralLine` whose `original_assignments` collect all listed classifications (blended components).

---


### Revised level energies

Wyart's workbook is an external, published list and is never edited. When the
identification work concludes that one of his levels sits at a different
energy, the new value is declared instead in the csv named by
`files.level_overrides` in `lineclass_config.toml`:

```
level_id,E_input,comment
059003.000565,118967.3776,"re-positioned from 119225.52; weighted mean of ..."
```

Three readers consult it, so the whole pipeline sees one energy for the level:

* `classify_lines.read_energy_levels()` generates the level's candidate
  transitions at the new energy;
* `chance_mc.read_input_levels()` supplies `E_input`, which `decoy_mc.py`
  displaces to plant the level's decoys and from which `level_shifts.py`
  measures every `dE`;
* `level_shifts.py --e-input <csv>` still exists and now merely repeats the
  same substitution; use it only to test a position that is *not* the adopted
  one.

A level in the file that is absent from the workbook raises an error rather
than passing silently as "nothing was moved". Comment the `level_overrides`
line out of the configuration to run on Wyart's energies unchanged.

#### The comment column, and what becomes of the published identifications

An identification taken from the published line list — the rows
`line_classifications.csv` marks `new` = 0, and which this README calls *old* —
names its two levels in the line workbook, and that workbook is an external
input and is never edited. So an identification made while a level sat at its
published energy goes on naming that level after the identification work has
moved it, and Step 3 keeps such an identification unless something rejects it
("old, no solid evidence for rejection"). The level is then fitted between its
true lines and its stale legacy ones. That is what pulled `059003.000424`
81 cm⁻¹ off its LOPT position on 2026-09-05.

The comment column is what tells the pipeline which of two things was done, and
`read_level_provenance()` reads it:

| the comment says | what was done | what happens to the published identifications naming the level |
|---|---|---|
| `re-positioned`, or `moved` | the level was found at a different energy | they lose their old status and are weighed as new candidates are, on their own evidence |
| `levels A and B were swapped`, or `exchanged with B` | two levels of the same parity and *J* exchanged their measured positions, each identifier keeping its own calculated identity | they keep their old status, **under the other identifier** — the observed lines did not move, only the name written over them |
| anything else | the energy was refined and nothing else | they stand exactly as written |

The phrase an exchange is recognised by is the one
`swap_line_assignments_pipeline.py` stamps on every row it writes, so an
exchange recorded by that script is recognised without anything further being
written by hand. A row that names an exchange partner absent from the level
list, or one whose partner's own row names a third level, stops the run: that
is a typo or a half-finished exchange, and either way there is no telling which
of two positions an identification belongs to.

`retag_legacy_identifications()` then settles, for every observed line, which
pairs of levels it holds a published identification for, and reports the count:

```
  Published identifications: 3462 stand as written, 20 carried over to the
  level that now holds the measured position, 10 withdrawn because a level of
  theirs has been re-positioned; the pairs they used to name are no longer
  candidates on their lines.
```

What this settles is which of a line's candidates are entitled to be treated
as old, and that entitlement belongs to the observed line and its measured
position rather than to the identifier written beside it in a workbook that
predates the move.

The candidate seeded from the workbook goes wherever its entitlement goes. A
line's published identifications are put on it as candidates whatever their
Ritz wavenumber — that is what makes them published identifications rather
than proposals, and it is why they can be weighed at all — so a pair the line
no longer holds has to be taken off it. Otherwise the row survives as a
candidate exempt from the matching tolerance, at a Ritz mismatch of the whole
distance the level moved: at 95158.229 cm⁻¹ the pair `000087-000398`, whose
upper level swapped its measured position away on 2026-09-06, stood in the
output with `dif_wn_O-C = 92.898` cm⁻¹, an identification nobody holds and
that nothing in the run can accept. For an exchange the identification is
re-seeded under the partner's identifier, which is where it now lives; a
sound exchange puts that pair within the matching tolerance anyway, so this
usually adds no row that was not there already.

The verdicts these rows carried in `line_decisions.csv` become moot when the
pair leaves the run, and are listed by `report_unapplied_decisions()` at the
end of the run — which is what that report is for. They can be deleted from
the ledger or left as a record of the decision; nothing reads them again.

The point of routing the revision through this file rather than through the
workbook is that the level then re-enters the probabilistic validation: its
decoys are planted around the position actually being claimed, so it is tested
on the same footing as every other new level. Note the limit this leaves: when
the revised energy was derived as the mean of the Ritz energies of the level's
own lines, its `dE` mainly measures internal consistency. That is the same
circularity Wyart's own new levels carry, which is why the two are comparable —
but the intensity-pattern test and the decoy tail remain the real evidence.

Do **not** feed LOPT's optimized energies back in as the next run's input
energies and iterate: that loop has no fixed point to defend. Set the input
energies once, from the Ritz-derived values, and report LOPT's optimized output
as the result.

### Levels found since the level list was published

`revised_level_energies.csv` can move a level Wyart already has. It cannot
create one. A level found in IDEN2 after the adopted list was published exists
in no input file of this pipeline at all: it has no row in the level workbook,
which is an external published list and is never edited, and it has no
calculated transitions in `Icalc.xlsx`. Nothing could therefore ever propose an
observed line for it, and a `line_decisions.csv` row accepting one of its lines
would be ruling on a candidate that is never generated — a silent no-op.

`files.new_levels`, named in `[files]`, is how such a level enters — and
`insert_new_level.py` (below) is what writes it, along with everything else the
new level needs.

**`files.new_levels`** (`new_levels.txt`) — one row per level, **tab**-separated:

```
level_id	E	J	parity	iden2_row	cowan_lid	comment
059003.000623	132958.176	7/2	o	474	396	Found in IDEN2; 7 lines …
```

The file is tab-separated and named `.txt` for one reason: Excel turns the J
value `7/2` into the date 7 February when it opens a `.csv`, and offers no way
to stop it, whereas opening a `.txt` gives the import dialogue, where the J
column can be declared Text. A `.csv` is still read — the delimiter follows the
extension — so the older form goes on working.

`iden2_row` and `cowan_lid` are optional and are written by
`insert_new_level.py`: the row of `IDEN2/enlev.dat` the level is, and its level
number in the Cowan calculation. `cowan_lid` is what its calculated transitions
are found by, so a level that carries one needs no `icalc_extra` rows at all
(see below).

`level_id` is the next one free: the last six digits of the largest identifier
in use, plus one. The level enters the run marked `is_new = 1`, exactly like a
level of the published list that is absent from the ASD — it is a level whose
reality the work is establishing, so the decoy calibration must test it on the
same footing as the others — and `is_added = 1` records that it came from here
rather than from the workbook. Reusing an identifier that is already in the
level list is an error, not a silent overwrite.

**The calculated transitions of a new level** are derived on the fly, by
`classify_lines.read_cowan_transitions()`, from the Cowan transition list
`tp_E1_no_trials.xlsx` — the very file `Icalc.xlsx` was made from, which holds
every calculated transition of the ion whether its levels are known or not.
The level's `cowan_lid` selects its rows; every one of them whose partner
carries a Wyart identifier the level list knows becomes a transition, with
`Icalc` computed from `gA` as below and the same `gA` cutoff applied that
`Icalc.xlsx` itself obeys, so that a new level's transition list stops where
every other level's stops and the completeness rule keeps meaning the same
thing. **`icalc_new.xlsx` is therefore no longer needed**, and the
`icalc_extra` line of `lineclass_config.toml` is commented out.

**`files.icalc_extra`** (formerly `icalc_new.xlsx`) — still supported, and
still first in precedence: a pair listed there overrides the derived row, so a
calculated transition can be corrected by hand without touching the
calculation. In the layout of `Icalc.xlsx` (same worksheet name, same column
names); only the pairs involving a new level need be listed. The complete
calculated table for the ion, including every level not yet identified, is far
too large to carry here; and a pair absent from *both* files is treated by the
completeness rule exactly as before, as one whose `gA` falls below the printing
cutoff of Cowan's codes (see *Transitions missing from `Icalc.xlsx`*).

**The `Icalc` column of the supplementary file is not read.** The predicted
intensity is recomputed from `gA` by the same relation the main file obeys,

```
Icalc = C * gA * (rwn/1e8) * exp(-Eup/kT)
```

with `C` and `kT` taken from `[intensity_model]`. The reason is that
`tools/fit_boltzmann.py` rewrites the `Icalc` column of the main file every
time `C` and `kT` move, and does not touch this one: reading its column would
leave a handful of transitions on a stale intensity scale, and the intensity
tests would then blame the resulting disagreement on the identification.
Recomputing costs nothing and cannot drift.

With both files in place, the new level's transitions are ordinary candidates:
they are matched, graded, weeded and ruled on by the ledger like any others.

Comment either line out of the configuration to run without them.

### The decision ledger: identifications ruled on by hand

The classification is re-derived from scratch on every run, so a verdict
reached by eye — in IDEN2, or by watching what a trial LOPT run does to the
residuals — is lost unless it is written down. Without such a record a
candidate that was refused last round is proposed again this round, hand-made
identifications disappear from a freshly generated LOPT input, and the accepted
set keeps moving for reasons that have nothing to do with any new evidence.
The **decision ledger** is that record: the csv named by `files.line_decisions`
in `lineclass_config.toml`.

```
wn_obs,low_id,upp_id,decision,date,reason
40047.2789,059003.000218,059003.000337,accept,2026-09-03,"well supported by the pattern of observed lines seen in IDEN2"
44022.4591,059003.000105,059003.000245,reject,2026-09-03,"firmly rejected, stays so"
```

* `wn_obs` — the observed wavenumber of the line, cm⁻¹. It is matched to the
  nearest observed line within `DECISIONS_WN_MATCH` = 0.01 cm⁻¹, so it may be
  written to fewer decimals than the line list carries; the closest two
  observed lines of this spectrum are 0.10 cm⁻¹ apart, so the match cannot be
  ambiguous. A wavenumber that matches no line stops the run.
* `low_id`, `upp_id` — the two levels of the assignment being ruled on. One
  line may carry several rows, one per assignment.
* `decision` — `accept` or `reject`. Anything else stops the run, as does one
  row contradicting another.
* `date`, `reason` — free columns. `reason` is echoed in the reports.

**An `accept` row is an order, not a vote.** The automatic matching proposes a
pair only when the observed wavenumber lies within 5.5 combined standard
deviations of the Ritz wavenumber *E*(upper) − *E*(lower). An identification
made by eye can be further out than that: the analyst has evidence the matching
does not use — the appearance of the line on the plate, the branch structure of
the level, what a trial LOPT fit does with it — and a hand-made identification
is frequently the very thing that pulls a level back to where it belongs. Six of
the present ones sit at *O*−*C* ≈ 0.3 cm⁻¹, several times the matching window.
Such a pair used to be unreachable: with no candidate on the line there was
nothing for the verdict to be applied to, and the row was merely reported as
unapplied at the end of the run.

So `force_ledger_assignments()` now **creates** the assignment when the matching
does not propose it. The calculated intensity and its uncertainty are copied
from the generated candidate for that pair, and the new transition is graded,
takes part in the conflict resolution and is accepted at the end of the weeding
exactly like any other — the only difference is where it came from. It is
rebuilt at the start of every cycle, so it survives the level optimization
moving its levels around.

Because an order that cannot be obeyed must not be silently dropped,
`check_forced_decisions()` runs before any of the work and stops the run,
listing **every** fault it found rather than the first, if:

* a row names a level id that is not in the level list (a typo, or a level that
  has been renumbered) — an error for a `reject` row too, since such a row can
  never rule on anything;
* an accepted pair breaks the electric-dipole selection rules used to generate
  candidates (opposite parity, |Δ*J*| ≤ 1, not both *J* = 0), or its Ritz
  wavenumber falls outside `[wn_min, wn_max]` — no such transition can exist;
* the same pair of levels is accepted on two different observed lines — one
  transition can belong to only one line, so the two orders contradict each
  other. (Accepting a pair on one line while rejecting it on another is *not* a
  contradiction: that is how an assignment is moved from one line to another,
  and both rows are obeyed. Two rows that accept and reject the same pair on the
  *same* line are caught earlier, when the file is read.)

With this in force, only `reject` rows should ever be reported as unapplied at
the end of a run — a rejected pair can still fall out of the matching window
when its levels move, which makes the verdict moot rather than lost.

`apply_line_decisions()` runs at the end of the weeding of each line, after
Steps 1–3 and after the oscillation blacklist, so the ledger has the **last
word** on every assignment it names. Two further places respect it:
`resolve_conflicts()` gives an accepted pair the transition outright (a verdict
could not be honoured later if the transition were taken away from the line
first), and the automatic verdict it replaces is kept in `notes2`, prefixed
`Manual Accepted (overrides)` / `(agrees)`, so the disagreements can be read
off the output. The `manual` output column holds `accept`/`reject` for the
rows the ledger rules on.

This is what makes a rerun reproducible in the only sense that matters here:
**only assignments the ledger does not name can move between two rounds.**
When the remaining movement is adjudicated and written into the ledger, the
output stops changing and the data files can be frozen.

Chance-coincidence and decoy runs ignore the ledger (`classify_lines.main()`
reads it only when `wn_shift` and `decoy_shift` are both zero): those runs
measure what the algorithm alone does with an input it should find nothing in,
and laying human verdicts over them would corrupt the false-positive rate.

A row whose two levels never come up as a candidate for that line is legal —
the levels may have moved far enough for the assignment to fall outside the
matching tolerance — and is listed at the end of the run by
`report_unapplied_decisions()`. A stale `reject` is harmless; a stale `accept`
means the identification is **not** in this run's output, which is why it is
printed.

## Output

Two files are written to the `LineClass/` directory (`OUTPUT_FILE`, `OUTPUT_CSV`):

| File                        | Description                                            |
|-----------------------------|--------------------------------------------------------|
| `line_classifications.xlsx` | Classification table with per-column number formatting |
| `line_classifications.csv`  | Identical data in CSV                                  |

#### A file open in Excel stops the run at the start, not at the end

Excel holds the workbook or comma-separated file it has open locked against
writing, so a run whose output file is open used to die at the last moment —
after every level had been read, every candidate weighed and every cycle
turned — with `PermissionError: [Errno 13] Permission denied`, and the whole
run was lost.

Every script that writes a file Excel can open now asks, before it starts
work, whether that file can be written, and stops with the name of the file
if it cannot:

```
Cannot write the output file of this run:
    line_classifications.xlsx - it is open in another program (Excel locks the
    file it has open) or is read-only

Close the file and run again.  The run stops here, before any of the work is
done, rather than after all of it.
```

The test is a real one, not a guess from the file name: the file is opened for
update — read and write, no truncation — which asks the operating system for
exactly the access the eventual write will need and touches not a single byte.
It lives in `output_files.py` (`require_writable`), and is made by
`classify_lines.py`, `make_LOPT_input.py`, `check_sync.py`, `chance_mc.py`,
`decoy_mc.py`, `level_shifts.py`, `level_positions.py` and
`level_interchange.py`. A run that writes no files (`main(write_files=False)`,
used by the Monte-Carlo drivers) skips it, as does any `--detail` mode, which
only prints.

The gap between the check and the write is the length of the run: nothing stops
a file being opened in Excel while the run is going. What the check settles is
the common case — the file was already open when the run began.

Rows are sorted by **decreasing observed wavenumber** (`wn_obs`), then decreasing Ritz wavenumber (`rwn`), then increasing `grade` for ties (`build_output`). Unclassified observed lines still appear, as a single blank-classification row.

### Output Columns (24)

| #  | Column             | Description                                                                |
|----|--------------------|----------------------------------------------------------------------------|
| 1  | `wn_obs`           | Observed wavenumber (cm⁻¹)                                                 |
| 2  | `unc_wn_obs`       | Wavenumber uncertainty σ (cm⁻¹)                                            |
| 3  | `obs_intens`       | Observed intensity (`Icor`)                                                |
| 4  | `char`             | Line character (`h`, `w`, `bl`, …)                                         |
| 5  | `low_id`           | Lower energy level id                                                      |
| 6  | `upp_id`           | Upper energy level id                                                      |
| 7  | `calc_intens`      | Adjusted calculated intensity (after per-level factors)                    |
| 8  | `orig_calc_intens` | Original theoretical intensity from `Icalc.xlsx`                           |
| 9  | `u_calc`           | Adjusted log-uncertainty of `calc_intens`                                  |
| 10 | `imputed`          | `1` = the pair is absent from `Icalc.xlsx` and its intensity was imputed   |
| 11 | `intens_from_f`    | Upper-level intensity correction factor (ln scale)                         |
| 12 | `intens_to_f`      | Lower-level intensity correction factor (ln scale)                         |
| 13 | `dif_wn_O-C`       | Observed − Ritz wavenumber (cm⁻¹)                                          |
| 14 | `grade`            | 2D grade: tier (`2`–`5`) + subgrade (`A`–`E`, `G`); see *Grading*          |
| 15 | `notes1`           | Conflict flags: `F` (conflicting), `R` (revised)                           |
| 16 | `notes2`           | Per-transition decision trace from weeding (`Step1`/`Step2`/`Step3` label) |
| 17 | `manual`           | `accept`/`reject` if the decision ledger rules on this assignment          |
| 18 | `new`              | `1` = new classification, `0` = original Sugar, blank = unclassified       |
| 19 | `accepted`         | `1` = accepted by weeding, `0` = rejected, blank = unclassified line       |
| 20 | `n_accepted`       | Number of accepted classifications of this observed line (0, 1, 2, …)      |
| 21 | `low_E`            | Lower level energy (cm⁻¹)                                                  |
| 22 | `upp_E`            | Upper level energy (cm⁻¹)                                                  |
| 23 | `rwn`              | Ritz wavenumber = `upp_E − low_E` (cm⁻¹)                                   |
| 24 | `BF`               | Branching fraction of this component (0 if not accepted)                   |

> **`BF` and why it is not called "weight".** The output column is the **branching fraction** `BF` of the component, which is exactly what LOPT takes as input: LOPT (centroid model) squares `BF` internally to form the weight, so you feed it `BF`, not `BF²`. The column is named `BF` rather than "weight" to avoid confusion with the script's own internal `weights` dictionary (`calc_weights()`), which stores the *squared* branching fractions `BF² / σ²` used for the internal level optimization. For a single accepted transition on a line, `BF = 1`. Rows with `accepted = 1` and their `BF` values form the input for LOPT. See `calc_weights()` and the `bf` computation in `build_output()`.

---

## Data Flow

```
Pr3_lev_Wyart_1999.xlsm ┐
Icalc.xlsx              ├─▶ read_* ─▶ EnergyLevel / Transition / SpectralLine  (models.py)
Pr3_lines.xlsx          ┘
                                 │
                                 ▼
                generate_all_possible_transitions   (Step 4)
                                 │
             ┌───────────────────┴────────────────────────────┐
             │            main() outer cycle (≤ max_cycles)    │
             │                                                 │
             │   compute_level_uncertainties ─▶ match_and_grade│
             │   ─▶ resolve_conflicts ─▶                       │
             │   weed_assignments (iterative) ─▶ calc_weights  │
             │             ─▶ optimize_levels                  │
             └───────────────────┬────────────────────────────┘
                                 │  converged?
                                 ▼
                 build_output ─▶ write_output (xlsx + csv)
```

The **outer cycle** in `main()` repeats *match → resolve → weed → weight → optimize levels* until the number of accepted transitions is unchanged **and** the largest level-energy change in the first optimization iteration is `< 0.001 cm⁻¹` (`max_cycles = 20`; with the current inputs the cycle converges in 8–9 cycles).

---

## The 6-Step Pipeline

The docstring and `main()` organize the work into six steps.

### Steps 1–3 — Read inputs
`read_energy_levels`, `read_transitions`, `read_observed_lines` load the three files into the dataclasses in `models.py`:

- **`EnergyLevel`** — id, energy, parity, J (string and float), plus per-level intensity factors, the estimated energy uncertainty `u_energy` (see *Level-Energy Uncertainty Estimation* below), and the lists `from_transitions` / `to_transitions`.
- **`Transition`** — lower/upper level, `calc_intensity`, `u_calc`, grade, `notes1`/`notes2`, `new`, `accepted`, and the *immutable originals* `orig_calc_intensity` / `orig_u_calc`. `calculated_wavenumber` is a property = `upper.energy − lower.energy`.
- **`SpectralLine`** — wavenumber, σ, intensity, character, `assigned_transitions`, and `original_assignments`.
- A module-level sentinel `UNASSIGNED = Transition(notes1='R')` marks a line whose only classification was revised away.

### Step 4 — Generate all possible transitions
`generate_all_possible_transitions` sorts levels by energy and forms every level pair that satisfies the **E1 selection rules**:

- opposite parity (`lev_up.parity != lev_lo.parity`),
- `|ΔJ| ≤ 1` (implemented as `abs(J_up − J_lo) > 1.5` → skip, to tolerate float noise),
- not both `J = 0`,
- Ritz wavenumber within **`[WN_MIN, WN_MAX] = [9327.0, 121665.0]` cm⁻¹**.

Each transition inherits `calc_intensity`/`u_calc` from `Icalc.xlsx` when the pair is present there, otherwise `None`. The list is sorted by `calculated_wavenumber` (ascending) to enable binary search.

### Step 5 — Assignment cycle
Before `assignment_cycle` runs, `main()` calls `compute_level_uncertainties(levels_dict)` (see *Level-Energy Uncertainty Estimation* below) to (re-)estimate each level's `u_energy` from the *previous* cycle's accepted transitions — `assignment_cycle` itself then runs three sub-steps:

**5.1 Match & grade** (`match_and_grade`, `assign_grades`).
For each observed line, a **bisection** over the sorted transition wavenumbers selects candidates within a tolerance of **`5.5 σ_comb`**, where `σ_comb = sqrt(wn_uncertainty² + u_energy(lower)² + u_energy(upper)²)` — the line's own wavenumber uncertainty combined in quadrature with the estimated uncertainty of each candidate's two level energies, so a transition to/from a not-yet-precisely-known level is not excluded just because the *line* itself was measured very precisely. Each candidate is reused if already attached to the line, else copied and attached. `assign_grades` sets the 2D grade and the `new` flag (`0` if the pair matches an `original_assignments` entry, else `1`); the grade's **tier** still uses the line's own `wn_uncertainty` alone (see *Grading* below) — only the matching window and the Step-1 accept/reject test (next) use `σ_comb`. Grading covers not only the fresh matches but also any seeded original assignment whose Ritz wavenumber fell outside the matching window in the current cycle — such an assignment still takes part in the weeding (e.g. in the Step-2 blend logic) and must not keep an undefined grade or `new` flag.

**5.2 Resolve conflicts** (`resolve_conflicts`).
When one transition is matched to several observed lines, the candidates are sorted and the best is kept. Four scoring approaches are present in the source; **Approach D — "Conservative Sort" is ACTIVE**, the other three (Lexicographic, chi², Rank-based) are commented out. The conservative key is:

```
(manual, is_new, |Obs − Ritz|, tier, z, int_err)
```

`manual` is `−1` when the decision ledger accepts this (line, transition) pair, `+1` when it rejects it, `0` otherwise, so a hand-ruled pair wins or gives up the transition outright — a verdict could not be honoured later in the weeding if the transition were taken away from the line here. Below it, legacy assignments are favored, then the smallest wavenumber mismatch, then tier, then z (σ residual), then intensity mismatch. The winner is flagged **`F`** (conflicting); if any *original* (legacy) assignment lost, the winner is also flagged **`R`** (revised). Losing transitions are removed from their lines; a line that loses its only classification (and had originals) receives the `UNASSIGNED` sentinel.

**5.3 Iterative intensity weeding** (`weed_assignments`) — see *Weeding* below.

### Step 6 — Energy level optimization
`optimize_levels` refines level energies — see *Energy Level Optimization* below.

---

## Grading (`assign_grades`)

Each match gets a two-character grade, **tier + subgrade**.

**Tier — wavenumber agreement** (`wn_diff_abs = |wn_obs − wn_calc|`; `σ` here is the line's own `wn_uncertainty` only — unlike the matching window and the Step-1 Ritz test, the tier does **not** include level-energy uncertainty, so a transition accepted thanks to a sizeable `u_energy` can still display a poor (`4`/`5`) tier):

| Tier | Condition                                              |
|------|--------------------------------------------------------|
| `2`  | `wn_diff_abs ≤ 2.0 σ`                                  |
| `3`  | `wn_diff_abs ≤ 3.0 σ`                                  |
| `4`  | `wn_diff_abs ≤ 4.0 σ`                                  |
| `5`  | otherwise (`> 4.0 σ`, up to the 5.5 σ match tolerance) |

**Subgrade — intensity consistency.** With `u = u_calc` and `d = |ln(I_calc / I_obs)| − u` (`ln2 = ln 2`, `ln5 = ln 5`):

| Subgrade | Condition                                                                               |
|----------|-----------------------------------------------------------------------------------------|
| `G`      | no theoretical intensity (`calc_intensity` or `u_calc` is `None`), or none of the below |
| `A`      | `d ≤ 0` and `u ≤ ln2`                                                                   |
| `B`      | `d ≤ ln2` and `u ≤ ln5`                                                                 |
| `C`      | `ln2 < d ≤ ln5` and `ln2 < u ≤ ln5`                                                     |
| `D`      | `d ≤ ln5 < u`                                                                           |
| `E`      | `d > ln5` and `u > ln5`                                                                 |

There is no subgrade `F`; `F` and `R` are `notes1` flags produced by conflict resolution.

---

## Iterative Weeding (`weed_assignments`)

Weeding decides, for each observed line, which candidate transitions to **accept** (`accepted = 1`) or **reject** (`accepted = 0`). It runs an **outer factor-refinement loop** around a **per-line 3-step decision process**.

### Outer factor loop
1. **Iteration 0 (baseline):** weed once with raw theoretical intensities (`calc_intensity == orig_calc_intensity`).
2. **Iterations 1…N** (`compute_intensity_factors`, `apply_intensity_adjustments`) — **disabled by default**: with `USE_INTENSITY_ADJUSTMENT = 0`, no factors are computed and every iteration uses the raw theoretical intensities. The machinery was evaluated and switched off because it added only ~8 accepted classifications (5080 vs 5072) at the cost of ~240 extra fitted parameters. When enabled (`USE_INTENSITY_ADJUSTMENT = 1`), it works as follows:
   - For every level compute a **per-level intensity factor** from its accepted, *unblended* transitions (`_compute_factor_for_transition_list`): qualifying transitions need `accepted == 1`, positive `orig_calc_intensity`/`orig_u_calc`, a parent line with exactly **one** accepted transition (`n_accepted == 1`), and `|ln(I_obs/I_calc_orig)| ≤ 3`. With **≥ `min_n` = 10** qualifiers, the factor is the plain (unweighted) mean of `ln(I_obs / I_calc_orig)` (`FACTOR_WEIGHTING = 'unweighted'`; the Mandel–Paule weighted mean is available as an option but is biased here, because the `1/u_calc²` weights correlate with the residuals). Factors are **clamped to ±5** ln-units.
   - The upper-level (`intens_from_factor`) and lower-level (`intens_to_factor`) factors are fitted by **backfitting** (alternating fits on partial residuals), so that shared structure is not absorbed twice by the two marginal fits.
   - **Damped update:** `applied = alpha × raw + (1 − alpha) × prev` with **`alpha = 0.5`** (no damping on the first factor computation).
   - **Apply** to every transition: `I_adj = I_orig × exp(f_from + f_to)`; uncertainty adds in quadrature. Originals are never mutated, preventing feedback amplification.
   - **Reset** all decisions (`reset_weeding_state`) and re-weed.
3. **Convergence:** stop when no transition changes acceptance state between passes, else after `max_iterations` (`assignment_cycle` passes `max_iterations = 100`; the function's own default is 59).
4. **Oscillation blacklist:** the blacklist is created once in `main()` and passed into `weed_assignments`, so it **persists across the outer cycles**. A transition enters it when its acceptance goes **A→B→A** over the last three snapshots — either three consecutive weeding passes within a cycle, or three end-of-cycle snapshots (`blacklist_cycle_oscillations` in `main()`); both call the one detector, `_find_oscillations`. Blacklisted transitions are **rejected** in Step 3 ("decisions oscillate between weeding iterations or cycles"). Blacklist keys are `trans_key()` tuples (lower id, upper id, identity of the assigned line), which remain stable when transition objects are recreated between cycles.

   **A→B→A is the whole criterion, and it was not always so.** Until 2026-09-03 a *single* flip between two consecutive weeding passes was enough, and absence from a snapshot counted as a state of its own. Both are wrong, and both were destabilizing:

   - the weeding loop is a damped fixed-point iteration on the per-level intensity factors, so a single flip is what convergence looks like from close up. Blacklisting on it withdrew whatever candidate happened to sit near a decision threshold at pass 3 — a set that depends on every other transition in the run, so a small change in the input (one level moved by a fraction of a cm⁻¹) reshuffled it and moved identifications that had nothing to do with the change;
   - absence from a snapshot means the transition was not a candidate that cycle at all, because the levels had moved and it fell outside the matching tolerance. That is a missing decision, not a wavering one. Counting it blacklisted **168** transitions in the third cycle of a normal run, while the energies were still moving by of the order of 1 cm⁻¹.

   With the genuine criterion the same run blacklists **18** transitions in total, of which 17 are still candidates at the end.

   Because a blacklisting is the one rejection the pipeline makes on **no evidence about the identification itself** — the transition was dropped because its acceptance would not settle, which says nothing about whether it is right — it must not be silent. `write_unstable_report()` writes every one of them to **`unstable_candidates.csv`** (observed wavenumber, both level ids, intensities, grade, which loop it oscillated in, and the three states) so they can be looked at in IDEN2 and settled for good in the decision ledger, which overrides the blacklist.

### Per-line 3-step decision (`weed_assignments_line`)
Candidates are sorted by decreasing `calc_intensity`, then by Ritz σ.

- **Step 1 — clear-cut decisions.** New transitions with Ritz mismatch `> 3 σ_comb` are rejected outright, and (for legacy candidates with no theoretical intensity) mismatch `> 5 σ_comb` is likewise an outright rejection — `σ_comb` is the same combined line+level-uncertainty sigma used for matching (`_ritz_sigma`, see *Match & grade* above and *Level-Energy Uncertainty Estimation* below). A **relative-intensity filter** (`_apply_relative_intensity_filter`) rejects new candidates contributing `< 10 %` of the accepted cumulative intensity, and applies conservative tests to legacy candidates around the `4 %` / `2 %` thresholds (`check_internal_significance` guarantees rejection of statistically-insignificant legacy lines). Surviving candidates are tested by **asymmetric z-score** against `I_obs` — thresholds are looser when the prediction is *too strong* and for **resonance** lower levels (ids ending `001`, via `_fudge_factor_for_asym_intensity`), stricter when *too weak*. Decisions are tagged `Step1…` in `notes2`.
- **Step 2 — grouping & late arbitration.** Undecided candidates are tested for **pair** and **triple** acceptance (`_try_pair_stage`, `_try_triple_stage`): a group passes if its intensity-weighted **center of gravity** matches the observed wavenumber and its **effective spread** — the Doppler width (full width at half maximum, `sqrt(8·ln2·kT/m)/c` = 8.22×10⁻⁶ of the wavenumber for Pr at `T = 1.6 eV`, `ATOMIC_MASS = 140.90765 u`) convolved with the measurement σ (`_effective_spread_combined`) — is consistent with the line profile. For `h`/`w` (and other broadened) characters the thresholds are relaxed by a broadening multiplier (`_broadening_multiplier`). A CoG-sensitivity test rejects candidates that *worsen* the group CoG by more than 1 σ. Finally, still-undecided **new** candidates carrying `F`/`R` in `notes1` are rejected. Decisions are tagged `Step2…`.
- **Step 3 — defaults.** Any remaining undecided **new** candidate is rejected ("no solid evidence for acceptance"); any remaining **legacy** candidate is accepted ("no solid evidence for rejection"). Blacklisted transitions land here and are annotated "decisions oscillate between weeding iterations or cycles". Decisions are tagged `Step3…`.
- **Last word — the decision ledger.** `apply_line_decisions(line)` then sets every assignment named in `files.line_decisions` to the verdict written there, whatever Steps 1–3 and the blacklist concluded, and records the overruled verdict in `notes2`. See *The decision ledger* above.

---

## Level-Energy Uncertainty Estimation (`compute_level_uncertainties`)

Both the matching tolerance and the Step-1 Ritz-mismatch test (above) compare an observed-minus-calculated wavenumber residual against a σ. Using the observed line's own `wn_uncertainty` alone silently assumes the two level energies that define the Ritz wavenumber are exact — reasonable for old, well-established levels, but not for a level newly proposed by Wyart, whose adopted energy carries its own (unstated) uncertainty. `compute_level_uncertainties` estimates that missing per-level uncertainty, `EnergyLevel.u_energy`, so it can be folded into the σ used everywhere above.

**When it runs.** Once per outer Step-5 cycle in `main()`, immediately after `clear_assignments` and *before* `assignment_cycle` (hence before `match_and_grade`). At that point `level.from_transitions`/`to_transitions` still hold the *previous* cycle's final, accepted transitions (`clear_assignments` resets line-level and transition-level state but does not touch these level-side lists), so the estimate always uses the most recent settled classification. On the very first cycle these lists are empty, so every level starts at `u_energy = 0.0` — consistent with "no assignments yet, no basis for an uncertainty."

**The estimate.** For a level `L` with an accepted transition to a partner level `P` (line wavenumber uncertainty `u_obs`), that transition implies `E_L = E_P(adopted) ∓ wn_obs`, with combined uncertainty `sqrt(u_obs² + u_energy(P)²)`. Collecting one such implied energy per accepted transition of `L` (both `from_transitions`, where `L` is the upper level, and `to_transitions`, where it is the lower level), `u_energy(L)` is set to the chi²-inflated weighted-mean uncertainty of that set — the same statistic (`_wm_red_chi`) already used for the per-level intensity factors (see *Iterative Weeding* above). A level with zero accepted transitions keeps `u_energy = 0.0` (unknown) until later cycles bring in assignments.

**Propagation through chains.** `u_energy(P)` in the formula above is itself an estimate being solved for, so the whole level system is updated by fixed-point (Gauss-Seidel-style) sweeps: every level is recomputed from the *previous* sweep's `u_energy` values, repeated until the largest per-level change drops below `tol` (or `max_sweeps` is reached — see *Configuration*). This lets uncertainty propagate outward from the well-determined old levels through however many links a chain of new levels has, without the algorithm ever needing to identify a single "defining transition" for a level in the middle of a ladder: each level simply combines whatever accepted transitions it has, weighted by their own (possibly still-converging) partner uncertainties, and a few sweeps are enough for the chain to settle.

**Leaving the candidate out (why the test is not circular).** Taken at face value, `u_energy` would make every acceptance test judge a transition against an error bar the transition itself helped to inflate: a line that disagrees with its level widens the scatter of that level's implied energies, the chi²-inflation in `_wm_red_chi` turns that scatter into a larger `u_energy`, and the larger `u_energy` shrinks the very `σ_comb` ratio the line has to pass. Three mutually incompatible lines could therefore hold each other up. The test is made non-circular by **leave-one-out**: `_ritz_sigma(candidate)` asks `u_energy_excluding(level, candidate)` for each of the candidate's two levels, which recomputes that level's chi²-inflated weighted-mean uncertainty from its *other* accepted determinations only (they are kept for this purpose in `EnergyLevel.u_contrib`, keyed by `trans_key`, which is the transition identity that survives the re-creation of `Transition` objects each cycle). A level whose only support is the candidate itself falls back to `u_energy = 0`, i.e. the candidate is judged on the observed wavenumber uncertainty alone — the honest answer, since without that transition the level has no measured energy. Nothing is capped: the *remaining* determinations still inflate normally when they genuinely disagree. No separate consistency-weeding step is needed, because each Step-5 cycle re-decides every assignment from scratch: the worst line of an inconsistent level now fails its own test, is dropped, and the level re-forms from what is left.

**Consequence for convergence.** Because `u_energy` is recomputed from scratch every cycle from that cycle's starting (previous-cycle) acceptance state, and the matching/Step-1 tolerances depend on it, a transition that is rejected in one cycle purely because a partner level's uncertainty was still underestimated can be picked up in a later cycle, once more of that level's other transitions have been accepted and its `u_energy` has grown to reflect the real scatter.

---

## Energy Level Optimization (`optimize_levels`)

After weeding, each non-ground level (`energy != 0`) is re-estimated as the **weighted mean** of the energies implied by its accepted transitions:

- from a transition *from* the level (it is the upper level): `E = E_lower + wn_obs`,
- from a transition *to* the level (it is the lower level): `E = E_upper − wn_obs`.

Weights come from `calc_weights`: for each accepted transition on a line, `BF = I_calc / Σ I_calc` over the line's accepted set (or `1/n` when any intensity is missing), and `w = BF² / σ_wn²` — the LOPT "centroid" weighting that shares a blend's weight among its components. Lines with a single accepted transition get `BF = 1` → `w = 1/σ²`. Iteration runs up to `MAX_IT = 50`, stopping when the maximum energy change `< TOL = 0.001 cm⁻¹`. The function returns the *first-iteration* maximum change, which the outer cycle uses as its convergence signal.

> This optimization is intentionally minimalistic — it exists only to sharpen the Ritz wavenumbers that drive assignment decisions. **Final** level optimization is deferred to **LOPT (v ≥ 5)**, fed by the accepted rows and their `BF` values.

---

## The intensity model (Boltzmann plot)

### The relation

`gA` is the statistical weight of the upper level times the probability per second that the
atom makes the transition; it is what Cowan's codes compute. In a plasma in local
thermodynamic equilibrium the number of atoms sitting in an upper level of energy `Eup`
falls off as `exp(−Eup/kT)`, where `kT` is an effective excitation temperature written as an
energy in cm⁻¹, the same unit as `Eup`. The observed intensity is the **energy flux under
the line contour**, so it carries one factor of the photon energy — that is, one factor of
the wavenumber `rwn`. Hence

```
Icalc = C · gA · (rwn / 1e8) · exp(−Eup / kT)
```

with `rwn` the Ritz wavenumber in cm⁻¹ and `1e8/rwn` the vacuum wavelength in ångström.
`C` absorbs everything constant across the spectrum: the number of emitters, the solid
angle, and the plate response, which the intensity calibration has already removed from the
observed intensities. Taking logarithms turns the relation into a straight line — the
**Boltzmann plot**:

```
y = ln[ Iobs · (1e8/rwn) / gA ]  =  ln C  −  Eup / kT
```

so a straight-line fit of `y` against `Eup` gives `kT` from the slope (`kT = −1/slope`) and
`C` from the intercept (`C = exp(intercept)`).

### Fitting it: `tools/fit_boltzmann.py`

```
python tools/fit_boltzmann.py                                  # fit and report
python tools/fit_boltzmann.py --write                          # ... and rewrite Icalc.xlsx
python tools/fit_boltzmann.py --plot fit.png --points pts.csv  # ... with diagnostics
```

It takes one point per accepted identification (`accepted = 1` in
`line_classifications.csv`), reading `obs_intens` and `rwn` from that file and `gA` and
`Eup` for the same level pair from `Icalc.xlsx`; a pair absent from `Icalc.xlsx` has no `gA`
and cannot enter the plot. It then fits twice: once on everything, then again after
dropping the points whose predicted intensity disagrees with the observed one by more than
a factor `e²` (`|ln(Icalc/Iobs)| > 2`, the `--clip` value). With `--write` it recomputes the
`Icalc` column of every row of `Icalc.xlsx` from the fitted constants and the file's own
`gA`, `Eup` and `rwn`, keeping a copy of the previous workbook as
`Icalc_before_boltzmann.xlsx`. Because openpyxl stores no cached result for a formula, the
`u_ln` column — written in the workbook as `=LN(1+u%gA/100)` — is materialised as the
number that formula gives, so that the file reads back correctly for every consumer.

The first fit, made on the 3300 identifications the old model had accepted (3296 of them
have a `gA`), gave

| pass | points | C | kT | rms of ln I |
|---|---|---|---|---|
| first  | 3296 | 193.05 | 12140.6 cm⁻¹ | 1.381 |
| second (after dropping 319) | 2977 | 175.0916 | 11953.797 cm⁻¹ (1.482 eV) | 0.822 |

but those constants are not the ones in force: they were superseded by the iteration
described below, which ended at `C` = 131.092, `kT` = 13009.1 cm⁻¹ (1.613 eV) — and those
in turn by the refit of the plate calibration, which is described further down and leaves
`C` = 252.641, `kT` = 12365.9 cm⁻¹ (1.533 eV), and finally by the refit that followed the
removal of the circularity from the Ritz acceptance test (see *Leaving the candidate out*
above), which leaves **`C` = 251.606, `kT` = 12375.4 cm⁻¹** (1.534 eV) in force.

`C` and `kT` live in `[intensity_model]` of `lineclass_config.toml`; the single
implementation of the relation is `gA_imputation.impute_intensity`, which the imputation of
censored transitions and `tools/fit_boltzmann.py` both call.
`gA_imputation.check_intensity_model` re-fits the constants on the file at every run and
warns if they have drifted apart by more than `verify_tolerance` (1%).

> **Superseded on 2026-08-28.** The `Icalc` column formerly obeyed
> `Icalc = C·gA·exp(−Eup/kT)/rwn` with `C = 135.8`, `kT = 12905 cm⁻¹` — proportional to the
> wavelength instead of to the wavenumber, a factor `rwn²` different in shape across the
> spectrum. Those constants came from a fit made with the DREAM `gA` values; with the
> present `gA` values they disagree badly with the observed intensities. Recomputing the
> column raised `Icalc` by a factor of 15.6 in the median (range 0.20 to 90.8).

The first run made with this model (un-iterated constants) is archived under
`baseline/boltzmann_model/` and the converged one under `baseline/boltzmann_converged/`
(classification, level-shift report, both Monte Carlos, run logs), beside the two policy
archives.

### Making the fit self-consistent: `tools/iterate_boltzmann.py`

The constants are fitted on the accepted identifications, and which identifications are
accepted depends on the calculated intensities, hence on the constants. A single fit
therefore leaves the model and the line list disagreeing with each other. The loop is closed
by

```
python tools/iterate_boltzmann.py            # to convergence, tol = 1%
python tools/iterate_boltzmann.py --dry-run  # one fit, print it, change nothing
```

which repeats: fit `C` and `kT` on the accepted lines → recompute the whole `Icalc` column
with them → write them into `[intensity_model]` → reclassify every observed line against the
new `Icalc` → fit again. The convergence test is on the calculated intensities themselves,
not on the constants:

```
max over all 30206 transitions of | Icalc_now / Icalc_previous − 1 |  <  0.01
```

i.e. no calculated intensity moved by as much as one per cent in the last round. Testing the
intensities rather than `C` and `kT` matters because the two constants trade off against each
other: a larger `C` with a smaller `kT` reproduces almost the same intensities over the
energy range where most lines lie, so the constants can still be moving while nothing the
pipeline uses is. The criterion above cannot be satisfied that way, since it looks at every
transition, including those at the ends of the energy range where the trade-off fails.

**A trap to know about.** The first fit of a run reads the observed intensities out of
`line_classifications.csv`, i.e. out of the *previous* classification. If the `Icor` column
of `Pr3_lines.xlsx` has just been changed, that file still holds the intensities from before
the change, so round 1 returns exactly the constants that were already in force, `Icalc`
does not move, and the loop declares itself converged after one round - on the old answer.
The reclassification it performs at the end of that round does refresh the file, so simply
running the tool a second time starts it correctly. Whenever the intensity column has been
rewritten, run `iterate_boltzmann.py` twice and take the second run as the real one. The
giveaway in the log is that the number of accepted identifications the round *reclassified*
differs from the number it *fitted on*: nothing moved, and yet the answer changed.

It converged in three rounds (2026-08-28):

| round | C | kT (cm⁻¹) | max \|ΔIcalc\| | accepted identifications after |
|---|---|---|---|---|
| 1 | 134.651 | 12915.77 | 0.855 | 4890 |
| 2 | 131.111 | 13005.43 | 0.0500 | 4894 |
| 3 | **131.092** | **13009.11** | **0.00293** | **4894** |

Almost all of the change is in the first round; the second moves every intensity by at most
5%, the third by at most 0.3%. Against the un-iterated fit the classification barely moves
(4839 → 4894 accepted identifications), which is the reassuring outcome: the fixed point is
close to the first fit, so the answer does not depend on where the iteration started.
The converged run is archived under `baseline/boltzmann_converged/` together with the round
logs.

That table is history: the observed intensities themselves were afterwards put on a
refitted plate calibration (see [The adopted calibration](#the-adopted-calibration)), and
the same iteration was run again on them - twice, in fact, the second time after Sugar's
own printed intensities replaced the reconstructed ones. The constants in force are the ones
that run ended at, `C` = 252.641 and `kT` = 12365.9 cm⁻¹. They were re-fitted once more on
2026-08-30, after the leave-one-out change to the Ritz test withdrew 6 identifications, and
the loop closed in a single round: `C` = 251.606, `kT` = 12375.4 cm⁻¹, with no calculated
intensity moving by as much as 0.5% (max 0.46%) and the re-fit on the reclassified list
returning the same two constants to six figures. On 2026-08-31 the plate calibration was
refitted once again, this time with the two coverage gaps of Sugar's exposures imposed as
region boundaries (see [The coverage gaps, and the refit that imposed them](#the-coverage-gaps-and-the-refit-that-imposed-them-2026-08-31)),
and the loop closed in two rounds at `C` = 268.527, `kT` = 12241.8 cm⁻¹. Those are the
values now in force.

Where the validation stood after the Ritz acceptance test was made non-circular, and before
the coverage gaps were imposed on the plate calibration (all figures from
`baseline/noncircular_ritz/`; for the state now in force see
[The coverage gaps, and the refit that imposed them](#the-coverage-gaps-and-the-refit-that-imposed-them-2026-08-31)): 4874 accepted identifications; 208 levels are
tested (a tested level is one supported by more new than previously published
identifications). The in-situ decoys — shadow copies of the tested levels displaced in
energy, so every line they accept is false by construction — put the best validation cut at
`n >= 4` accepted lines with `|dE| <= 0.1374 cm^-1`, giving **113 of the 208 levels
validated with 4.75 expected false among them** (4.2%). The shifted-wavenumber cross-check,
a looser upper bound, passes 113 of 208 at the same support with 2.62 expected false.
Folding the energy shift and the intensity pattern into a per-level spurious probability
leaves **199 levels below 5%, 8 between 5% and 50%, one above 50%**, for a most probable
3.0 spurious levels among the 208 (95% range 0.0–10.0). Seven levels come out at
`p_spur >= 0.1`, and the masking check clears one of them (059003.000457, 96% of its missing
predicted intensity hidden under observed lines), leaving six questionable levels retained.

Almost all of the change from the previous state is one level. `059003.000622` had been held
up by three mutually incompatible lines that the circular test let through; with the
circularity gone it keeps a single line, and its spurious probability moves from 0.005 to
0.873 — the level is not refuted, but it is no longer supported, and the report now says so.
`059003.000268` fell out of the tested set for the neutral reason that it lost two new lines
and its previously published support (6) now outnumbers its new support (5), so it counts as
a reference level instead. Every other tested level moved by less than 0.06 in `p_spur`.
The rise of the most probable number of spurious levels from 1.8 to 3.0 is therefore not a
loss of quality: it is one falsely validated level being correctly re-labelled as doubtful.

For comparison, three earlier states of the same validation. **Immediately before the
circularity was removed** (`baseline/sugar_intensities_complete/`): 4880 accepted
identifications, 209 tested levels, the cut at `n >= 4`, `|dE| <= 0.1289 cm^-1` validating
112 of 209 with 5.50 expected false, and 202 / 7 / 0 levels below 5%, between 5% and 50%,
and above 50% `p_spur`. With Sugar's intensities used for
the 6760 lines credited to him and the remaining 23 still reconstructed
(`baseline/sugar_intensities_converged/`): 4884 accepted identifications, the same cut
validating 112 of 209 with 5.25 expected false, and 204 / 5 / 0 levels below 5%, between 5%
and 50%, and above 50% `p_spur`. **Before the plate calibration was refitted at all**
(`baseline/boltzmann_converged/`): 4894 accepted identifications, the cut at `n >= 4`,
`|dE| <= 0.1299 cm^-1` validating 115 of 209 with 5.75 expected false (5.0%), and 196 / 7 / 6.
The recalibration therefore left the number of validated levels essentially unchanged while
removing every level that had stood above 50% spurious probability; the last step, closing
the 23 remaining gaps, moved nothing beyond the ordinary jitter of the two Monte Carlos.


### The plate calibration behind the observed intensities

The linearized observed intensities in `Pr3_lines.xlsx`, and the plate calibration
`intensity_correction_functions.txt` that produced them, were originally derived by hand
from a similar Boltzmann model but with the DREAM `gA` values. They have since been
refitted with the new `gA` values by `tools/calibrate_intensities.py` and adopted into
`Pr3_lines.xlsx` by `tools/apply_calibration.py`, both starting from the intensities Sugar
actually printed, which `tools/attach_sugar_intensities.py` brought into the line list as
the column `Iorig`; all three are described below. Note
that the noise-level formula of `level_shifts.py` (`1000·exp(P(λ))`, the linear-scale
intensity of Sugar's faintest recordable line) needs no adjustment for the new model: it
concerns only the observed intensities being brought to a common linear scale, and the
refitted `Icalc` is directly comparable with those linearized observed intensities.

### Sugar's own intensities: `tools/attach_sugar_intensities.py`

Everything below rests on the numbers Sugar printed beside each line - small whole numbers
(1, 3, 20, 500, … 9000) expressing how black the line came out on the photographic plate.
Those numbers were the starting point of the `Icor` column, but for a long time they were
not kept anywhere in `Pr3_lines.xlsx`, so the calibration tools had to *reconstruct* them
by undoing whatever correction was in force - and could only do so as accurately as that
correction file happened to be.

They are now in the line list, in the column **`Iorig`**, copied straight from the two
manually checked extractions of Sugar's tables:

| source label | source workbook, sheet `Table 1` | wavenumber range measured | wavenumber column | intensity column |
|---|---|---|---|---|
| `Sugar1969` | `Pr3_Sugar69_extracted_gemini-3-flash-preview_v1.xlsm` | 47541 - 121663 cm^-1 | `wn_mean` (24) | `Intensity` (3) |
| `Sugar1974` | `Pr_3_Sugar74_Table1_extracted_v3_gemini-3-flash-preview.xlsm` | 9329 - 47440 cm^-1 | `wn adopted` (31) | `Intens` (3) |

Lines are matched **by wavenumber, not by wavelength**. Sugar printed both for every line,
and the two do not always agree to his last digit; each extraction workbook therefore
carries a reconciled value - a weighted mean of the wavenumber he printed and the one
implied by the wavelength he printed - and it is that value, and only that, from which
`Pr3_lines.xlsx` was built. The tool is told which column holds it, takes the nearest
wavenumber, and refuses a match farther away than `--tol` (default 1e-4 cm^-1).

They are **not** matched through the `ref` column of the line list. That column says whose
*classification* a row carries, which is a different question from who *measured* the line:
the rows credited to `Wyart1974` were measured by Sugar like all the others, and their
intensities stand in his workbooks. Which of the two workbooks holds a given line is
settled by the line's wavenumber, because the two recordings cover disjoint stretches of
the spectrum - everything above about 47500 cm^-1 was recorded in 1969 and everything below
it in 1974 (the ranges in the table leave a 100 cm^-1 gap between them). No boundary has to
be configured: the tool looks each line up in every source given and keeps the nearest
wavenumber found in any of them, which is the same thing whenever the sources are disjoint
and the better choice where they are not. The label on the left of a `--source` is
therefore only a name, used in the report and in the audit file.

```bash
python tools/attach_sugar_intensities.py \
  --source "Sugar1969=Pr3_Sugar69_extracted_gemini-3-flash-preview_v1.xlsm:Table 1:wn_mean:Intensity" \
  --source "Sugar1974=Pr_3_Sugar74_Table1_extracted_v3_gemini-3-flash-preview.xlsm:Table 1:wn adopted:Intens"
```

A column may be named by its header (letter case, blanks and embedded line breaks ignored)
or by its 1-based number, which is safer when the header contains a line break. The
workbook is copied to `Pr3_lines.bak_iorig.xlsx` first, only the one column is written -
created at the right-hand end if it is not there yet, updated in place if it is, so a
second run never adds a second copy of it - and a row-by-row audit goes to
`iorig_match.csv`.

**How it came out.** All 6783 lines of the working list matched, the worst wavenumber
difference being 4.4e-11 cm^-1 - the last binary digit of a stored number, i.e. exact.
115 source lines serve two rows each, which is correct: an observed line with two candidate
identifications appears twice in the line list. 2883 of the 2946 lines of the 1969 paper
and 3785 of the 4333 of the 1974 paper were used; the rest belong to lines the working list
does not carry. `Iorig` runs from 1 to 9000 and has no empty cell.

The reconstruction it replaces was in fact good: over the 6760 lines matched at the first
attempt, the previously restored intensities and the true ones differ by a standard
deviation of 8% (in natural logarithms), 99% of them within ±15%, with one line off by a
factor 245. The refit below therefore moves the answer only slightly - but it now rests on
the measurements instead of on an inversion of a hand-made file.

The 23 rows credited to `Wyart1974` were left empty at that first attempt, when the source
file was chosen by `ref`; they were filled in afterwards, once it became clear that `ref`
names the author of the classification and not the observer. Their reconstructed
intensities turn out to have been good too - the true numbers differ from them by between
-7% and +21%.

Its tests are in `tests/test_attach_sugar_intensities.py`.

### Refitting the plate calibration: `tools/calibrate_intensities.py`

`intensity_correction_functions.txt` holds the correction that turns Sugar's reported
plate intensities into the uniform linear scale of the `Icor` column of `Pr3_lines.xlsx`:

```
Icor(lambda) = 1000 * I_Sugar * exp(P(lambda))
```

with `lambda` the vacuum wavelength in angstroms (`1e8/wn`) and `P` a piecewise
polynomial - the file gives, one region per text line, `lambda_lo lambda_hi c0;c1;c2;...`
meaning `P = c0 + c1*lambda + c2*lambda^2 + ...` inside that region. The regions are the
pieces the spectrum was recorded in: inside one, the sensitivity of plate, grating and
optics varies smoothly with wavelength; between two, it jumps. That file was built by
hand. `tools/calibrate_intensities.py` builds it automatically, boundaries included.

It is a **standalone** tool: it imports nothing from the pipeline (only numpy, openpyxl
and, for `--plot`, matplotlib), and takes every file, sheet and column name as an option.

```bash
# the run that produced the correction now in force: it fits Sugar's own
# intensities, so nothing has to be undone first
python tools/calibrate_intensities.py --col-intensity Iorig --cover 821.92 10721.57 \
       --out-functions intensity_correction_sugar.txt --plot intensity_calibration.png

# the older way, when the originals are not written down: undo the correction
# that is in force and fit what is left
python tools/calibrate_intensities.py --restore-from intensity_correction_functions.txt --no-restore-round

# just undo the correction and report the range of the restored intensities
python tools/calibrate_intensities.py --restore-from intensity_correction_functions.txt --restore-only
```

**What it does.** It joins the identified lines of `Pr3_lines.xlsx` to `Icalc.xlsx` on
their two level ids, so that each has a calculated intensity
`Icalc = C*gA*(rwn/1e8)*exp(-Eup/kT)` - the energy-flux form, proportional to the
wavenumber, used everywhere in this project. The disagreement `dlnI = ln(Icalc/Iobs)`
plotted against wavelength is the correction to be found. Then:

1. **Regions.** Two kinds of boundary. The **coverage gaps** are given in advance
   (`--gap LO HI`, repeatable; the defaults are Sugar's two, `1522.49-1529.85` and
   `2103.46-2107.92` A, the same numbers as `COVERAGE_GAPS_A` in `level_shifts.py`;
   `--no-gaps` switches them off). No exposure covered those intervals, so no line was
   ever recorded on both of the plates that meet there, nothing tied their intensity
   scales together, and the correction is genuinely discontinuous by an amount no data
   can reveal. They are therefore cuts *before* any fitting: the lines on the two sides
   are segmented and fitted independently, and the gap interval itself is covered by no
   region (a wavelength inside it, or beyond the ends, takes the value of the nearest
   region at that region's own end - never an extrapolation). Fitting one polynomial
   across a gap would smear a real, arbitrarily large step over the whole neighbourhood.
   Everything else is found in the data: recursive binary splitting. Each region's
   polynomial degree is chosen by
   the Bayesian information criterion `BIC = n*ln(SSE/n) + k*ln(n)` (`n` points, `SSE`
   sum of squared residuals, `k` fitted coefficients); a cut is kept when it lowers the
   total BIC of the two halves by more than a further `--split-penalty*ln(n)`. A region
   that still wants a degree above `--max-degree` (default 5) is split instead of fitted
   stiffly - the rule "if the required power is too high, break the region into smaller
   pieces". No region may fall below `--min-points`, and a polynomial that would run away
   at a thinly populated end of its region is refused in favour of a lower degree
   (`--max-swing`).
2. **The fit.** `ln C`, `1/kT` and every polynomial coefficient enter the model linearly,
   so all of them are found in one least-squares solution rather than by alternating
   between the Boltzmann line and the polynomials. That matters: `Eup` and `lambda` are
   correlated, a change of temperature can be partly mimicked by a tilt of `P`, and
   alternation passes that tilt back and forth for dozens of rounds. Solved jointly, the
   whole calibration converges in 3 or 4 rounds - the rounds being needed only because
   the region boundaries themselves are decided anew from each fit.
3. **Convergence** is measured as in `iterate_boltzmann.py`: the largest relative change
   of any corrected intensity between two rounds, `--tol` (2%) for the preliminary stages
   and `--tol-final` (0.2%) for the last.
4. **Outliers**, in physically motivated stages, the whole fit repeated after each and
   each stage looking again up to `--max-passes` times, never by a blunt cut on `|dlnI|`:
   - *stage 1* nothing removed;
   - *stage 2* of the lines with `|dlnI| > --dln-cut` (2), those whose calculated `gA` is
     uncertain by more than `--u-cut` percent (50) - for them the disagreement is
     plausibly the calculation's fault. The alternative `--u-rule sigma` instead removes
     those whose `|dlnI|` is within `--u-sigma` uncertainties of their own `gA`
     (`u_ln = ln(1+u%gA/100)`), and both counts are printed whichever rule is in force;
   - *stage 3* of the lines still too **weak** (`dlnI > +cut`, the only sign
     self-absorption can produce), those ending on a level at or below
     `--self-abs-elow` (default 0, the ground level: resonance lines). The tool also
     prints how many anomalously weak lines end on levels below 1000, 2000, 5000 and
     10000 cm^-1, so the threshold can be reconsidered from the data;
   - *stage 4* what is left: dropped, and the fit repeated at `--tol-final`, **only if**
     every region loses at most `--max-drop-frac` (5%) of its lines and keeps at least
     `--min-points`. Otherwise the tool stops with a message and writes no correction
     file, because the residual disagreement is then not a matter of a few bad lines.

**Where it lands on the Pr III data.** The run in force fits the `Iorig` column, i.e.
Sugar's own numbers, so nothing is restored and nothing can be lost in restoring. From the
3486 identified lines that have a `gA`, it converges to 9 regions and `kT` = 12253 cm^-1
(1.52 eV), with an rms `dlnI` of 0.807 (a factor 2.2), after 509 lines removed in stage 2,
none in stage 3 (no resonance line is anomalously weak) and 27 in stage 4. The region
boundaries are

```
1431.09  2098.04  2176.43  2450.71  2766.63  3979.61  4875.46  6877.19  A
```

which land near the hand-made ones (1530, 2108, 2752, 4455, 6050 A) - about as close as two
different segmentations of the same jumps can be expected to agree - and within a few
angstroms of the earlier automatic fit that worked from restored rather than true
intensities. The last region carries a `swing` of 0.97 against the 1.0 the tool allows: its
polynomial is only weakly held down at the long-wavelength end, where the lines thin out.

**Why one pass is now enough.** Fitting `Iorig` makes the calibration independent of
everything the pipeline does downstream. It reads Sugar's intensities, Sugar's own
identifications (the `id1`/`id2` columns of `Pr3_lines.xlsx`) and the `gA`, `Eup` and `rwn`
of `Icalc.xlsx` - never `Icor`, never our classification, never the fitted `C` and `kT`. So
the correction cannot drift as the classification settles, and re-running the tool after
the whole chain has converged reproduces the file byte for byte. Before `Iorig` existed the
tool had to undo the correction that was in force, which fed its own output back into its
input and made the fixed point something to be checked rather than guaranteed.

**A note on snapping.** When the originals must be restored rather than read
(`--restore-from` without `Iorig`), they are *not* snapped to whole numbers here
(`--no-restore-round`), although Sugar reported whole numbers. The reason is in the report
the snapping prints: over all 6783 lines the nearest whole number was up to 16% away, and
one line in a hundred more than 14% away, concentrated below 1000 A. The stored
coefficients of the hand-made correction were not accurate enough to land back exactly on
the original integers there, so snapping would have invented errors of that size instead of
removing them.

Two cautions, both printed by the tool. The excitation temperature it returns is an
*effective* number, not a measurement of the plasma: `Eup` and `lambda` are correlated,
so part of a real temperature effect can be absorbed into `P(lambda)` and the other way
round. And the stage-2 rule is decisive here - with the flat 50% threshold the procedure
completes as above, while `--u-rule sigma` removes only 133 lines, leaves 336 outliers,
and the tool refuses at stage 4 because 13 of its 15 regions would lose more than 5% of
their lines. Nine `gA` values in ten in this file are uncertain by more than 50%, so the
flat threshold removes nearly every outlier while a 2-sigma test against those same
uncertainties explains almost nothing. The output files are
`intensity_correction_auto.txt` (same format as the input) and
`intensity_calibration_lines.csv` (per-line `P`, `dlnI`, and which stage dropped it).

The correction is fitted on the identified lines only (879 A upwards here), but it has to
be applied to **every** line that carries an intensity, so it is made to cover the full
wavelength range of the line list: the outermost regions are stretched to reach the ends,
and their degree is lowered until the polynomial no longer runs away over the stretched
part, where no line holds it down. `--cover LO HI` widens that further; here it is set to
the wavelength range of the configured spectral interval (`[range]` of
`lineclass_config.toml`, 9327 - 121665 cm^-1, i.e. 821.92 - 10721.57 A), a shade wider
than the observed lines themselves, so that the noise-level formula of `level_shifts.py`
finds a correction at every wavenumber it is asked about instead of returning "no value"
at the two ends.

One trap worth naming, for the `--restore-from` route only: it must be given the correction
**currently applied to the column**, not the hand-made one it descends from. Once a new
correction has been adopted, that file is the new one; restoring with the old one instead
silently leaves the new correction inside the "original" intensities and the fit wanders
off (in this data it ended by refusing at stage 4). Fitting `Iorig` avoids the question
entirely.

Its own tests are in `tests/test_calibrate_intensities.py`: a synthetic spectrum whose
scale is spoiled by a known two-region function, which the tool must find back, breakpoint
included.

### Adopting a new correction: `tools/apply_calibration.py`

The intensity column of `Pr3_lines.xlsx` does not hold Sugar's numbers; it holds them
already corrected, `Icor = 1000 * I_Sugar * exp(P_old(lambda))`. A new correction can
therefore not simply be multiplied in: the old one has to be divided out first.
`tools/apply_calibration.py` does exactly that, for every line that has a wavenumber and
a positive intensity:

```
I_Sugar  = Icor / ( 1000 * exp(P_old(lambda)) )
Icor_new = 1000 * I_Sugar * exp( P_new(lambda) )
```

Undoing is only necessary while `I_Sugar` is not written down anywhere. Now that it is,
`--orig-col Iorig` takes it from the column instead, with no undoing and no snapping to
whole numbers, so the first line above is skipped and the second is exact. Rows whose
`Iorig` is empty would still go the long way round through `--old`; since every line of
`Pr3_lines.xlsx` now has one, none do. The overall factor `K` can then no longer be guessed
from the restored values, so it is given explicitly; `--scale 1000` is the scale this
project has always used.

```bash
# the run that adopted the correction now in force
python tools/apply_calibration.py --orig-col Iorig --scale 1000 \
       --old intensity_correction_functions.txt --new intensity_correction_sugar.txt

python tools/apply_calibration.py --old ... --new ... --dry-run   # report only
```

The workbook is copied to `Pr3_lines.bak.xlsx` before it is touched, only the one column
is written, and a per-line record of the change goes to `intensity_rescale.csv`
(`wn`, `lambda`, the old intensity, `P_old`, the restored original, `P_new`, the new
intensity). Because the operation is a change of scale and not an accumulation, applying
`old -> new` and then `new -> old` returns the column exactly as it was; that round trip
is one of the tests in `tests/test_apply_calibration.py`.

### The adopted calibration

The automatic correction, fitted to Sugar's own intensities, is the one in force. The files
now stand as follows:

| file | what it is |
|------|------------|
| `intensity_correction_functions.txt` | the correction in force; `level_shifts.py` reads this name for its noise level |
| `intensity_correction_sugar.txt` | the same file under the name it was produced with (the fit to `Iorig`) |
| `intensity_correction_auto.txt` | what `calibrate_intensities.py` last produced - identical to the file in force, which is the fixed-point check |
| `baseline/before_coverage_gaps/` | line list, `Icalc.xlsx`, classification, level-shift report, both Monte Carlos, LOPT output, config and correction as they were before the coverage gaps were imposed |
| `intensity_correction_manual.txt` | the original hand-made correction, kept for reference |
| `Pr3_lines.bak.xlsx` | the line list with the previous intensities (before this adoption) |
| `Pr3_lines.bak_iorig.xlsx` | the line list before the `Iorig` column was added |
| `baseline/sugar_intensities_complete/` | the run in force: line list, `Icalc.xlsx`, classification, level-shift report, both Monte Carlos, run logs, config and correction |
| `baseline/sugar_intensities_converged/` | the previous run, with the 23 lines credited to Wyart still on reconstructed intensities |
| `baseline/before_sugar_intensities/` | line list, `Icalc.xlsx`, classification, level-shift report, config and correction as they were before Sugar's own intensities were used |
| `baseline/before_intensity_recalibration/` | the same, one step further back: before the plate calibration was refitted at all |

Adopting it changed every intensity, but only slightly, because the correction it replaced
was itself a fit to nearly the same numbers: the median line moved by a factor 0.99, the
middle 90% of them by factors between 0.92 and 1.08, and the largest single change was a
factor 4.4. Every line in the column now satisfies
`Icor = 1000 * Iorig * exp(P(lambda))` to within 6e-16, i.e. exactly.

Since the classification weighs an identification partly by how well the observed intensity
agrees with the calculated one, the fit and the classification then had to be brought back
into agreement with `tools/iterate_boltzmann.py`. That took two rounds (`--tol 0.005`,
after the first, stale-input run described above) and ended at

```
C = 252.641,  kT = 12365.92 cm^-1,  rms |ln(Icalc/Iobs)| = 0.867
```

with no calculated intensity moving by more than 0.31% in the last round. `C` is not
comparable with its value under the earlier hand-made correction: `P(lambda)` is fixed only
up to an additive constant and `C` absorbs it. `kT` is comparable, and the sequence of
values it has taken is 13009 (hand-made correction) -> 12358 (refit from restored
originals) -> 12425 (refit from Sugar's own numbers for the lines he was credited with) ->
**12366** cm^-1 (the same with his numbers for every line in the list) -> **12375** cm^-1
(after the leave-one-out change to the Ritz test, 2026-08-30) -> **12242** cm^-1 (after
the coverage gaps were imposed on the calibration, 2026-08-31; see the next subsection).

On the classification the change is small, as it should be: **4880** accepted
identifications, the same number as under the previous correction, with the grades barely
moving (main grade 2/3/4/5 = 4500/335/37/8 against 4505/334/37/8; subgrade A/B/C/D/E/G =
1633/1963/109/211/69/895 against 1633/1966/117/210/71/887).

**The loop closes, and this time by construction.** Running `calibrate_intensities.py`
again after the whole chain had converged reproduced `intensity_correction_functions.txt`
byte for byte, which is no longer a coincidence to be checked but a consequence of the fit
reading only inputs that the pipeline does not write (see *Why one pass is now enough*,
above).

### The coverage gaps, and the refit that imposed them (2026-08-31)

Everything above found its region boundaries in the data alone. That was wrong at two
wavelengths. Sugar's spectrum was not recorded in one continuous sweep: no exposure covers
1522.49-1529.85 A, and none covers 2103.46-2107.92 A. Across such a gap no line was ever
recorded on both of the plates that meet there, so nothing tied their intensity scales
together and Sugar had no way of stitching them; the correction is therefore genuinely
discontinuous at each gap, by an amount that no amount of data can reveal. The automatic
segmentation had run single polynomials straight across both - the worst case was the old
region `2098.04 - 2176.43 A`, which straddled the second gap and tried to absorb the step
with a cubic (`-816142; 1141.6; -0.532; 8.3e-5`), smearing a true discontinuity over its
whole neighbourhood.

The gaps are now given to `calibrate_intensities.py` in advance (`--gap LO HI`, defaults
`DEFAULT_GAPS_A`; they were also the starting numbers of `COVERAGE_GAPS_A` in
`level_shifts.py`, which now uses the measured coverage function instead and keeps that list
only as a fallback — the two answer different questions and are no longer the same,
see [Where the plates were blind](#where-the-plates-were-blind-toolscoverage_mappy)) and cut the
lines into blocks *before* any fitting, so the two sides of a gap are segmented and fitted
independently and the step between them is unconstrained. The refit gives 8 regions,
`C` = 0.2622, `kT` = 12226.9 cm^-1, rms `dlnI` = 0.8175:

| region | from (A) | to (A) | degree | lines | rms dlnI |
|---|---|---|---|---|---|
| 1 | 821.92 | 1522.49 | 3 | 787 | 0.829 |
| | *coverage gap* | | | | *P jumps by +3.45* |
| 2 | 1529.85 | 2103.46 | 1 | 399 | 0.861 |
| | *coverage gap* | | | | *P jumps by +6.76* |
| 3 | 2107.92 | 2476.72 | 3 | 560 | 0.785 |
| 4 | 2476.72 | 2883.39 | 3 | 400 | 0.798 |
| 5 | 2883.39 | 3889.18 | 3 | 382 | 0.889 |
| 6 | 3889.18 | 4875.46 | 2 | 140 | 0.835 |
| 7 | 4875.46 | 6877.19 | 2 | 123 | 0.628 |
| 8 | 6877.19 | 10721.57 | 3 | 168 | 0.739 |

The second jump is large - a factor ~860 in Sugar's scale - and it is what the lines ask
for: the four lines between 2107.92 and 2125 A have `ln(Icalc/Iobs)` = 5.4 +- 0.3 against
about -0.5 just below the gap. Only its first few angstroms are thinly populated (the
region as a whole holds 560 lines), so the value 6.29 exactly at the edge is a short
extrapolation; the swing diagnostic there is 0.89, just under the `--max-swing` limit of 1.

Adopting it (`apply_calibration.py --orig-col Iorig --scale 1000`) moved the median line by
a factor 1.01, the middle 90% by factors 0.70 to 1.42, and the largest single line by 37;
the column again satisfies `Icor = 1000 * Iorig * exp(P(lambda))` to 6e-16. The Boltzmann
iteration then closed in two rounds (`--tol 0.005`, after the stale-input round) at

```
C = 268.527,  kT = 12241.81 cm^-1,  rms |ln(Icalc/Iobs)| = 0.876
```

with no calculated intensity moving by more than 0.47% in the last round, and
`classify_lines.py` re-run on that model reproducing the same answer. `kT` from the
independent joint fit inside the calibration (12227 cm^-1) and `kT` from the pipeline
(12242 cm^-1) now agree to 0.1%.

On the classification: **4891** accepted identifications against 4874 before (main grade
2/3/4/5 = 4535/316/37/3 against 4522/312/35/5; subgrade A/B/C/D/E/G =
1609/1929/130/212/69/942 against 1628/1962/109/210/69/896). On the validation: 208 tested
levels, best decoy-calibrated cut at `n >= 5` and `|dE| <= 0.2335 cm^-1` validating **114
of 208** with 4.75 expected false, and 198 / 9 / 1 levels below 5%, between 5% and 50%, and
between 50% and 90% spurious probability, for a most probable **3.5** spurious levels among
the 208 (95% range 0.0-10.8). The state before this change is archived in
`baseline/before_coverage_gaps/`.

### The short-wavelength failure, and the refit that cured it (2026-09-05)

The correction above was wrong below about 1000 A, badly, and the way it showed up is
worth recording because the diagnostic is general.

**The symptom.** Level `059003.000533` has three strong lines that IDEN2 shows aligning
cleanly with it — 121229.786, 120101.207 and 118982.601 cm⁻¹ (824.9, 832.6 and 840.5 A).
On the log scale IDEN2 displays, `I_log = 10*ln(I_lin)`, their predicted intensities are
97, 93 and 97 and their observed ones 91, 86 and 97: agreement to a few tenths in `ln`.
But `line_classifications.csv` gave observed intensities of 211, 169 and 604 against
calculated 11143, 7591 and 11816 — the observed values too small by factors of 20 to 50.

**The cause.** Region 1 ran 821.92 to 1522.49 A and was fitted with a **cubic**,

```
P = -93.735889568457 + 0.227980354340891*lam - 1.80299782559081e-4*lam^2
    + 4.57997976538457e-8*lam^3
```

whose value at 824.9 A is **-2.654**, against **+0.938** for the hand-made linear fit of
the same region — a factor `exp(3.59) = 36` in the corrected intensity, which is exactly
the shortfall observed. Two things went wrong together:

* The region was fitted on 879.44-1522.49 A but had to *cover* 821.92-1522.49 A, so the
  last 57 A were pure extrapolation, and a cubic extrapolates violently.
* Nothing caught it. The existing guard, `--max-swing`, asks how far `P` goes outside the
  band the region's own lines occupy — their 1st to 99th percentile of `dlnI`. That band
  is about two units wide in `ln`, because `gA` is a calculated quantity and scatters, so
  a three-unit dive at the edge sat comfortably inside it. Swing was 0.00.

**The fix, in two parts.**

`--max-extrap` (new, default 0.5 in natural logarithms) measures the *right* thing: how far
`P` moves, over the stretched part of a region, away from the value it takes at the
outermost line that was actually fitted. That is not scatter — it is one systematic factor
applied to every corrected intensity out there — and `_tame()` now lowers the degree of an
edge region until both this and `--max-swing` are satisfied. Region 1 became a parabola.

`--dln-cut` (changed) now takes a **schedule**, by default `4 3 2`, one value per pass of
outlier stage 2. On the first pass the correction is still the rough stage-1 shape, and a
line whose only fault is that it sits where the shape has not been found yet was being
condemned on the strength of an error the fit was about to remove. Loosening the early
passes and tightening as the fit settles is what was always done by hand. It keeps more
lines: 432 dropped at stage 2 out of 3486, where the fixed cut of 2 had dropped more, and
region 1 is now fitted on 808 lines instead of 641.

**The result.** 8 regions, `C` = 0.2671, `kT` = 12081.4 cm⁻¹, rms `dlnI` = 0.8356:

| region | from (A) | to (A) | degree | lines | rms dlnI |
|---|---|---|---|---|---|
| 1 | 821.92 | 1522.49 | 2 | 808 | 0.850 |
| | *coverage gap* | | | | |
| 2 | 1529.85 | 2103.46 | 1 | 416 | 0.890 |
| | *coverage gap* | | | | |
| 3 | 2107.92 | 2450.71 | 3 | 523 | 0.797 |
| 4 | 2450.71 | 2801.04 | 2 | 413 | 0.814 |
| 5 | 2801.04 | 3833.82 | 3 | 429 | 0.906 |
| 6 | 3833.82 | 4865.93 | 2 | 145 | 0.830 |
| 7 | 4865.93 | 6877.19 | 2 | 124 | 0.636 |
| 8 | 6877.19 | 10721.57 | 4 | 169 | 0.741 |

At 824.9 A the new `P` is **-0.378** against the old **-2.654**: the three lines of
`000533` gained a factor 9.7, 8.3 and 7.1 and now read

| line (cm⁻¹) | I_obs before | I_obs now | I_calc | grade |
|---|---|---|---|---|
| 121229.786 | 211 | 2056 | 11113 | 2G |
| 120101.207 | 169 | 1402 | 7571 | 2C |
| 118982.601 | 604 | 4297 | 11785 | 2B |

Everywhere else the change is small: adopting it moved the median line by a factor 0.998,
the middle 90% by 0.86 to 1.35, and the largest single line by 10.4. The Boltzmann
iteration closed in **one** round (`--tol 0.005`) at

```
C = 269.980,  kT = 12233.07 cm^-1,  rms |ln(Icalc/Iobs)| = 0.876
```

with no calculated intensity moving by more than 0.47%.

**What is still not settled below 1000 A, and cannot be.** The fitted region-1 lines are
almost all longward of 1000 A: 5 lines lie below 950 A and 52 below 1000 A, out of 808.
The remaining disagreement with the hand-made linear correction — about 1.3 in `ln`, a
factor 4, at 825 A, where elsewhere the two agree to 0.2-0.5 — is therefore not a defect
of either fit but the honest width of what five lines can determine. The residual
`ln(Icalc/Iobs)` of the three `000533` lines is now 1.69, 1.68 and 1.01 against a region
rms of 0.85, i.e. about two standard deviations: no longer an anomaly, but a hint that the
straight line may still be the better shape out there. Anyone who needs intensities below
1000 A should treat them as uncertain by a factor of a few, whichever correction is in
force.

---

### The noise level and signal-to-noise ratios: `tools/estimate_snr.py`

The plate calibration says how the intensity *scale* varies with wavelength. It does not
say how faint a line could be and still have been recorded. That second number — the
detection threshold — is what decides whether a predicted transition that the line list
does not contain is genuinely absent or merely too weak to have appeared, and IDEN2 has no
way to show it: its display carries intensities, so a screen full of predicted transitions
that could never have been seen looks exactly like a screen full of real candidates.

`tools/estimate_snr.py` estimates the threshold and turns both files into signal-to-noise
ratios. It is standalone, like `calibrate_intensities.py`, whose region-fitting machinery
it reuses.

```bash
python tools/estimate_snr.py                                  # estimate only
python tools/estimate_snr.py --iden2-dir IDEN2 --out-dir IDEN2_snr --plot noise_level.png
```

**How the threshold is estimated.** A line list is a censored sample: everything in it was
above the threshold and everything below it is missing, so the weakest lines present sit
just above the threshold and it can be read off as their lower envelope. The lines are
sorted by wavelength and cut at the coverage gaps (the same `1522.49-1529.85` and
`2103.46-2107.92` A: nothing ties the noise levels of two plates that never shared a line);
a window of `--window` (default 20) consecutive lines is slid along each block, and the
weakest line in each window is the noise there, attributed to the window's median
wavelength. Windows do not overlap by default, because the region-finding charges for each
breakpoint against the number of points and repeated minima would make it believe it has
far more evidence than it does. `ln N` is then fitted against wavelength by the same
piecewise low-degree polynomials, written out in the same `lo hi c0;c1;...` format.

Outlier rejection is **asymmetric**, and that is the one place the algorithm departs from
the intensity calibration. A window whose weakest line sits far *above* the curve is a
window where the lines happened to be sparse, so its minimum overestimates the threshold:
it is dropped beyond `--high-cut` (1.0). A window whose minimum sits *below* the curve
caught a line very close to the threshold, which is the quantity being estimated: it is
kept unless it is absurd, beyond `--low-cut` (2.5). Rejecting both sides equally would pull
the estimate up toward the middle of the intensity distribution, which is not what a
threshold is.

**What it found for this list, and why the answer is not what one expects.** On Sugar's
*reported* scale the threshold is nearly flat: 84% of the windows have their weakest line
at intensity 1, the smallest value the list uses at all, and the fitted `N` runs between
1.0 and 2.3 across the whole spectrum. That is not a failure of the estimate — it is the
correct answer, and it says something about the line list. Sugar graded each exposure from
its own weakest visible line upward, so his reported intensity *already is* a per-plate
scale whose floor is the local threshold. The whole wavelength dependence of detectability
therefore sits in the plate calibration, and

> the SNR of an *observed* line is essentially its reported intensity; the useful half of
> the exercise is the other one — bringing a *predicted* intensity down onto that scale.

**Rewriting the IDEN2 files.** Both `dlv.dat` (one row per observed line, intensity in
columns 1-5, wavenumber in columns 6-18) and `TRANS.DAT` (one row per predicted transition)
carry `v = 10*ln(Icor)`, with `Icor = 1000 * Iorig * exp(P(lambda))` and `P` the
calibration in force *when the files were made*. That this is exactly the relation was
checked by regressing the `dlv.dat` integers on `ln(Icor)`: slope 10.016, intercept -0.19,
scatter 0.10 in `ln`, i.e. rounding. So `--trans-correction` must name the correction the
files were made with — `intensity_correction_manual.txt` here — not the one now in force:
the old one is baked into the numbers and has to come out again.

A `TRANS.DAT` row carries **two** intensities and **two** wavenumbers, and the conversion
depends on telling them apart:

| columns | contents |
|---|---|
| 6-10 | the **predicted** intensity of the transition |
| 11-22 | the energy of the level at the other end of it |
| 23-24 | an asterisk when that level's energy is experimentally known |
| 25-38 | the wavenumber of the **predicted transition** |
| 39-43 | the intensity of the **observed line** identified with it, if any |
| 45-57 | the wavenumber of that observed line |
| 69-74 | its row number in `dlv.dat`, `0` when the prediction is unidentified |

The asterisk in columns 23-24 means the row cannot be taken apart by splitting on spaces;
the fields have to be cut at their fixed positions. The wavelength at which the plate
calibration and the noise are evaluated is `1e8 / wn` from **columns 25-38**. Both
intensity columns are rewritten: the predicted one by

```
Icor_pred  = exp(v/10)
Iorig_pred = Icor_pred / (1000 * exp(P(lambda)))
SNR_pred   = Iorig_pred / N(lambda) * exp(z)
```

and the observed one by looking up its own line in the list, exactly as `dlv.dat` is. The
two columns sit side by side on the screen and are read against each other, so leaving one
of them on the plate scale would show a difference of some sixty units at the
short-wavelength end that is an artefact of the two scales and no part of any disagreement
between the prediction and the measurement.

The zero point `z` is there because calculated intensities carry an arbitrary overall
normalisation, so without it "SNR = 1" would sit anywhere. `--trans-zero auto` (the
default) measures it on the rows `TRANS.DAT` itself marks as identified — each names the
observed line assigned to it, so no matching of our own is needed — and takes the median
difference in logarithms; for the current files it rests on 5197 comparisons, 4948 of them
identified in the file, and the predicted scale ran a factor 1.39 high. That median is
biased upward — a
prediction is likelier to have found a line when it is strong — so the zero point is
conservative: it calls a transition unobservable slightly more readily than the truth
warrants. Pass a number to override it, or `0` to leave the predicted scale alone.

The number written is linear if every observed line falls between 1 and 999, so the column
reads directly as "how many times the noise", and `round(10*ln(SNR))` otherwise
(`--scale`). For this list the observed SNR runs from 0.46 to 4550, so the log scale is
chosen. Nothing else on a row is touched — the fixed-width layout is preserved byte for
byte — and the rewritten files go to `--out-dir`, never over the originals.

**What this buys.** Of the 104271 predicted transitions in `TRANS.DAT`, **70% fall below
SNR 0.1** and a further 18% below SNR 1; only 12% could have been recorded at all. Those
are the rows that can be dropped from the display and from the column ordering. On the
4948 rows that *are* identified, the predicted SNR now sits a median 1 unit (0.1 in `ln`)
above the observed SNR of the line, in every wavelength decade from 822 to 16500 A.

**What the estimate is not.** `N` is a detection threshold read off the line list, not a
measurement of the noise on the plate. Where the compiler of the list stopped measuring
weak lines, `N` rises whether or not the plate was noisier there; where the spectrum is
crowded, the weakest recorded line is closer to the true threshold than where it is empty.
It is the right quantity for "could this predicted transition have appeared in *this*
list", and it should not be quoted as anything else.

Outputs: `noise_level.txt` (the fitted `ln N` region by region), `line_snr.csv` (per line:
wavenumber, wavelength, reported intensity, noise, SNR, and the integer IDEN2 would carry),
`noise_windows.csv` (the window minima the fit rests on, and which were dropped), and
optionally a plot.

### Where the plates were blind: `tools/coverage_map.py`

**The question.** A predicted transition that theory says should be strong, and that the
line list does not contain, is evidence against whatever identification predicted it — but
only if a line *could* have been recorded where it falls. Three things stop that, and none
of them has anything to do with the atom: the wavelength fell between two exposures; a
defect of the emulsion sat on it; a line of another species sat on it. For this spectrum we
have a record of none of the three. Sugar made hundreds of exposures with varying plate
positions, sources and spectrographs, and their boundaries were never published; the
defects were never catalogued; the impurity spectra were never listed. The two intervals
this project used to carry hard-coded (`COVERAGE_GAPS_A` in `level_shifts.py`,
`DEFAULT_GAPS_A` in `tools/calibrate_intensities.py`) were an informed guess.

**The idea.** The record is not needed. All three causes have the same observable
consequence: *wherever they act, the line list is empty* — not merely of the transition we
care about, but of every line of every species that would otherwise have been recorded
there. So the line list maps its own blind spots, and it does so without our having to know
which of the three causes was at work.

What the tool produces is a **coverage function**

    c(lambda) = P(a line the spectrum really contains at vacuum wavelength lambda, bright
                  enough to be recorded, was in fact recorded), relative to how well that
                  is done in a typical part of the spectrum.

It is deliberately *relative*. The absolute part — how bright a line must be — is the noise
level `N(lambda)` of the previous section. The two multiply and neither can replace the
other: `N` is smooth in wavelength and says how faint a line may be before it is lost; `c`
is sharp and local and says where a line of *any* brightness is lost.

**The model.** Work is done in `u = ln lambda`, so that a bin of fixed width is a fixed
*fraction* of the wavelength — necessary because the list spans 822–10720 Å and the density
of recorded lines falls by more than a factor of ten across it. The default bin is
`--bin 2e-4`, 0.02 % in wavelength: 0.16 Å at 822 Å, 2.1 Å at 10700 Å. Each bin `b` is in
one of two hidden states, *covered* or *blocked*, and the number of lines it holds is
Poisson:

    covered:   n_b ~ Poisson(R_b)
    blocked:   n_b ~ Poisson(r0 * R_b)

with `R_b` the lines the bin would hold if covered and `r0` the **leak-through**, the
fraction that still gets recorded in a blocked bin — not zero, because obscuration is
rarely total (a scratch hides part of a bin, an impurity line hides only what it overlaps)
and because an exact zero would let one stray line veto an otherwise unmistakable hole.
Blocking comes in runs, so the states follow a two-state Markov chain with
`a = P(covered→blocked)` and `b = P(blocked→covered)` per bin; the standard
forward-backward recursion then returns `c_b = P(bin b is covered | the whole line list)`.
Because the recursion sees the whole list at once, an empty bin inside a populated stretch
is a fluctuation and keeps `c ≈ 1`, while an empty bin inside a desert is condemned by the
company it keeps.

**Why the baseline comes from theory.** `R_b` cannot be estimated from the observed list
alone. The density of recorded lines in this spectrum swings by a factor of twenty between
neighbouring stretches — real structure, since transition arrays cluster in wavelength — and
a smooth curve fitted through it either follows that structure (and then bends into a gap,
hiding it) or does not (and then hands the real structure to the blocked state). Both
failures were observed in practice. The way out is that the **calculated** transitions know
the structure and know nothing about the plates. `TRANS.DAT` holds 104271 predicted
transitions; their density per unit `u`, blurred by `--pred-blur` (default 0.002 in
`ln lambda`, a few Å in the ultraviolet, for the error of a calculated wavelength) gives the
baseline its *shape*. What is left for a fitted curve is only the slow **efficiency** — what
fraction of the calculated lines a covered plate actually recorded, which also absorbs the
non-Pr III part of the list — and that is fitted as a penalized cubic spline in `u` by
Poisson regression with the current coverage as an exposure offset, so that blocked bins do
not drag it down. `--no-pred` falls back to a spline-only baseline.

The check that justifies this: across every candidate blind stretch the predicted density is
flat while the observed density collapses. In 1530–1570 Å theory predicts 28.4 transitions
per Å against 34.0 per Å in the dense 1450–1516 Å next door, a difference of 16 %; the
observed density falls from 3.92 per Å to 0.05, a factor of **65**. Nothing atomic happens
there. The same holds at 2110–2150 Å (predicted 21.1 vs 21.9 per Å; observed 0.80 vs 2.64)
and at 890–960 Å (predicted 57.4 vs 64.0; observed 0.27 vs 2.32).

**How flexible the efficiency curve is allowed to be** is not left to taste, because the
answer depends on it: a stiff curve blames the atom's structure on the plates, a flexible
one bends into the holes and hides them. `--knots` is set generously (200) and the curvature
penalty `--smooth` — which is what really controls the flexibility — is chosen by BIC over
`--smooth-grid`, charging `ln(number of bins)` for every *effective* parameter the spline
spends (the trace of the smoother matrix, not the coefficient count). With the current
inputs the criterion picks penalty 1000, an effective **28.7 parameters** over the whole
range, and the leak-through settles at `r0 = 0.23`.

**What it finds.** Ten stretches with `c < 0.5`, holding 13.2 % of the wavelength range but
only 2.6 % of the observed lines and 2.7 % of the accepted classifications:

| λ range, Å (vac) | width, Å | lines expected | observed | lost | mean c |
|---|---:|---:|---:|---:|---:|
| 821.94 – 827.88 | 5.9 | 48.6 | 9 | 38.6 | 0.075 |
| 887.73 – 966.49 | 78.8 | 130.4 | 24 | 104.2 | 0.030 |
| 1164.29 – 1174.82 | 10.5 | 24.4 | 7 | 13.5 | 0.233 |
| **1522.75 – 1663.82** | **141.1** | 206.2 | 25 | 179.9 | 0.011 |
| **2103.32 – 2190.91** | **87.6** | 271.2 | 71 | 198.8 | 0.009 |
| 2781.29 – 2808.68 | 27.4 | 41.5 | 15 | 22.6 | 0.178 |
| 3211.42 – 3256.05 | 44.6 | 32.3 | 9 | 19.9 | 0.142 |
| 4935.80 – 5067.84 | 132.0 | 19.0 | 3 | 13.8 | 0.132 |
| 8755.78 – 9033.27 | 277.5 | 32.4 | 10 | 20.0 | 0.119 |
| 10447.00 – 10720.04 | 273.0 | 8.0 | 1 | 5.9 | 0.260 |

**Both hard-coded gaps are real, and both are truncated.** 1522.49–1529.85 Å and
2103.46–2107.92 Å come out blind (mean `c` = 0.06 in each), and their *lower* edges are
right to a fraction of an angstrom — 1522.49 against a measured 1522.75, 2103.46 against
2103.32. Their upper edges are not: each declared gap covers about 5 % of the blind stretch
it begins. Where the record actually resumes is 1663.8 Å and 2190.9 Å. The pattern is the
same in both cases and suggests how the numbers were arrived at: the start of each gap was
identified correctly and the recovery point was read off far too early.

Two blind stretches of comparable weight were not declared at all: 887.7–966.5 Å (130
expected lines, 24 recorded) and 1164.3–1174.8 Å. The stretch at the very short-wavelength
end, 821.9–827.9 Å, is the edge of the survey itself rather than a hole inside it.

**The intensity calibration keeps the narrow gaps, and that is a result, not an oversight.**
`tools/calibrate_intensities.py` forces a region boundary at each of these same two places
(`DEFAULT_GAPS_A`), for a different reason: across a seam between two exposures no line was
recorded on both plates, so nothing ties the two intensity scales together, and a polynomial
fitted through the seam would smear a real, unknown step over the whole neighbourhood. The
obvious move, once the blind stretches were measured, was to widen those two boundaries from
1522.49–1529.85 Å and 2103.46–2107.92 Å to the measured 1522.75–1663.82 and
2103.32–2190.91 Å. It was tried, and it is wrong.

A blind stretch is where *few lines were recorded*; a calibration gap is where *the
intensity scale breaks*; the two have different lengths. Past a seam the new plate is far
less sensitive near its own short-wavelength edge, so it registers few lines for some way —
which is exactly why the coverage map calls the stretch blind — but the lines it does
register are on the new plate's scale, and a polynomial that starts at the seam corrects
them correctly. The direct test is the residual left by the correction now in force:

| blind stretch, Å | median ln(Icalc/Icorrected), below / inside / above |
|---|---|
| 887.7 – 966.5 | — / +0.16 / +0.03 |
| 1522.7 – 1663.8 | +0.03 / +0.51 / +0.04 |
| 2103.3 – 2190.9 | 0.00 / −0.05 / +0.01 |
| 2781.3 – 2808.7 | +0.11 / −0.49 / −0.08 |
| 3211.4 – 3256.0 | −0.06 / −0.24 / −0.13 |
| 8755.8 – 9033.3 | +0.27 / −0.10 / −0.27 |

Inside every one of them the correction is already right, to better than 0.5 in natural
logarithms against a line-to-line scatter of about 1.2. Widening the two gaps instead throws
the 27 lines the fit was using inside them out of the fit and replaces their correction by
the clamped value of the neighbouring region: the lines at 2103–2191 Å come out a median
factor of 3 too faint, and the first of them, at 2110.6 Å, a factor of 700. Forcing all
eight measured stretches as boundaries fails the tool's own safety check outright — the cut
at 8755.8–9033.3 Å leaves 35 lines beyond it, below `--min-points`.

**The intensity scale is what tells the two causes of blindness apart.** Across these two
seams the median of ln(Icalc/Iobserved) steps by +4.4 and +3.7 — from −0.45 to +3.99 at
2103 Å, from −2.52 to +1.18 at 1523 Å — while the scatter inside the stretch stays no larger
than outside it. That is a new plate on its own scale, not a scatter of randomly obscured
lines. Across the other six measured stretches the median does not step at all. So the two
stretches the coverage map confirms are the two exposure seams, and the six it adds are
emulsion defects and impurity lines *inside* one exposure: they hide lines, which is what
the observability term needs to know, and they break no scale, which is what the calibration
needs to know. The declared edges are not guesses either — each is exactly the empty
interval between the last line recorded before the seam and the first one after it
(1522.489 → 1529.851, 2103.455 → 2107.918).

**What it is not.** The table is a measurement of where the record is thin, not a list of
plate boundaries: it says *that* something was in the way, never *which* of the three causes
it was. At the red end especially, a low-density stretch may be a genuinely empty piece of
spectrum. The table is meant to be consumed as the continuous function `c`, which weakens
the evidence of a missing line in proportion to how thin the record is; read as a set of
hard on/off boundaries it would throw that proportion away.

**Weight of the correction.** About **15 %** of the predicted transitions above the noise
fall where `c < 0.5` — roughly one predicted branch in seven that is currently counted as
missing evidence is missing for instrumental reasons. Against the 14 levels whose
questionable mark rests on an absent strongest branch (`top1 = missing`,
`question_status = retained`), the map changes nothing: every one of the 14 has its
strongest branch in a fully covered region (`c` ≥ 0.997). It does move weaker branches of
those levels — for instance 059003.000447 has a predicted intensity 2744 branch at
1642.66 Å with `c = 0.000`, and 059003.000575 one of intensity 2070 at 952.62 Å, also
`c = 0.000` — which used to be counted against them and no longer are.

`level_shifts.py` now reads this file. `coverage()` interpolates `c` between the bin centres,
`in_coverage_gap()` is `c < COVERAGE_MIN` = 0.5, and `observation_probability()` multiplies `c`
by the measured detection curve of the next section. The count of predictions dropped as
unobservable rises from **111** under the two hand-entered intervals to **3813** of 29260, and
1271 more sit in the half-lit band 0.5 ≤ `c` < 0.95, kept but discounted by their `c` wherever
a probability rather than a filter is wanted.

Usage and outputs:

```bash
python tools/coverage_map.py            # defaults; from inside LineClass/
python tools/coverage_map.py --no-pred  # spline-only baseline, for comparison
```

`coverage_map.csv` (one row per bin: wavelength range, lines recorded, lines expected, `c`),
`coverage_gaps.txt` (the segment table above), `coverage_map.log` (the BIC scan and the
fitted `r0`, `a`, `b`) and `coverage_map.png` (recorded density, baseline and `c`, in five
panels).

### How often a line is lost where the map says the plate was looking: `tools/obscuration_rate.py`

**What the coverage function cannot see.** `c(λ)` is measured by asking where the line list
is *empty over a stretch of wavelength*, and it has to be: one empty bin proves nothing, a
hundred consecutive ones prove a great deal. That makes it blind, by construction, to
obscuration that acts on a single line — a grain flaw a few tenths of an angstrom across, one
impurity line, a scratch narrower than the model's own resolution. Such an event thins the
record by one line, which no statistical model can distinguish from a fluctuation, and it
leaves no mark: nothing in the line list says where a line should have been and was not.

That matters more than it sounds, because the missing branches the validation actually argues
about are almost never in a mapped blind stretch. Of the 14 levels whose questionable mark
rests on an absent strongest branch, **every one has that branch at `c` ≥ 0.997** — in a part
of the spectrum the map calls fully covered. For those levels `c` contributes nothing, and
what is needed instead is a *rate*:

    ε(λ) = P(a line that was certainly bright enough to be recorded, at a wavelength
             the coverage map calls covered, was nevertheless not recorded)

so that the observability of a predicted transition factorizes as

    P(recorded) = c(λ) · (1 − ε(λ)) · D(I_pred / N(λ))

with `N(λ)` the noise level of `tools/estimate_snr.py` and `D` the probability that a line of
that predicted strength clears it. `c` is measured where the record is thin, ε where it is
not; neither substitutes for the other.

**How ε is measured.** Not on the levels under test — for those, whether a predicted line is
really missing is the very question in dispute. It is measured on the **established levels**:
supported by more old than new identifications *and* present in the ASD compilation, 384 of
the 594 in the current run. Their energies are certain to a few thousandths of a wavenumber,
so their predicted transitions fall at known places, and an absence there is the plate's doing
and not the identification's. Each of the 13709 predicted transitions between two such levels
is asked three questions.

*Where would it fall?* At the Ritz wavelength λ = 10⁸/(E_upper − E_lower), vacuum angstroms.

*How bright would it be?* `I_pred` from `Icalc.xlsx`, against the noise level at that
wavelength — the linear intensity corresponding to Sugar's plate intensity 1. Only the ratio
`I_pred / N(λ)` is used; it already absorbs the wavelength dependence of the plate's
sensitivity.

*Was anything recorded there?* Three outcomes. **Recorded** — a measured line lies within
`W = 5.5·σ_λ` of the Ritz wavelength, where σ_λ is the local median wavelength uncertainty of
the measured lines (0.003–0.009 Å over the whole range, so W ≈ 0.022 Å; 5.5 is the pipeline's
own matching factor). The line need not be classified as this transition, or as Pr III at all
— the question is only whether the plate registered anything there. **Masked** — nothing
within W, but a recorded line close enough and strong enough to swallow the prediction, by the
same test `level_shifts.py` uses for its own masking check. **Absent** — neither. Recorded and
masked together count as *explained*; ε is measured on the rest.

**The correction the measurement lives or dies by.** The window is 0.022 Å wide and in the
crowded ultraviolet the recorded lines are 0.13 Å apart, so about one predicted position in
eight lands on a recorded line by pure coincidence; left uncorrected that would hide a fifth of
the obscuration. The chance rate is measured, not modelled: every prediction is re-tested at
18 control positions displaced by ±6, ±9, … ±30 window widths, and the fraction of *those* that
come out explained is that prediction's chance rate `f`. Displacing in units of the window
keeps the controls inside the same local line density at every wavelength. The check: the
share of controls that land on a recorded line, 0.0831, matches the Poisson estimate
`1 − exp(−2Wρ)` from the local line density ρ, 0.0819. With `y_i = 1` for explained,

    ε = 1 − (Σ y_i − Σ f_i) / Σ (1 − f_i)·c_i

and the interval comes from the profile likelihood of `p_i = f_i + (1 − f_i)·c_i·(1 − ε)`.

**Why only the strongest predictions count, and why that is itself a result.** The absence of a
*weak* predicted line means nothing; absence becomes informative only where `D` — the
probability of clearing the noise — is 1. The honest way to find that place is not to model `D`
but to watch where the measured absence rate stops falling:

| I_pred / noise | predictions | ε |
|---|---:|---:|
| 10¹·⁵ – 10² | 526 | 0.110 |
| 10² – 10²·⁵ | 232 | 0.028 |
| 10²·⁵ – 10³ | 73 | 0.058 |
| 10³ – 10³·⁵ | 27 | 0.037 |

It stops at a **hundred times** the noise level. That sounds extravagant, and the reason it is
needed is worth stating. Model `D` as a log-normal — the predicted intensity times a scatter
factor whose spread is measured on the accepted identifications, sd of ln(I_obs·BF / I_pred)
= 1.30 over 3689 established-level lines, with BF the branching fraction that splits a blended
feature among its components — and `D` reaches 0.99 already at *ten* times the noise, so the
absence rate there should be 0.01. It is 0.14. The log-normal is far too optimistic at the
faint end for the obvious reason: it is fitted to lines that *were* recorded, so its lower tail
has been cut off by the very effect it is being asked to predict. Fitting the tail and ε
together does not rescue it — the two are degenerate, and a free fit puts the sd at 2.4 and ε
at 0, attributing every absence to the tail. What breaks the degeneracy is that at a hundred
times the noise a log-normal of *any* plausible width predicts essentially no loss, so whatever
is left there is not the intensity model.

The consequence for the likelihood is that its missing-line term must use the **measured**
detection curve and not a log-normal around the noise level, and `level_shifts.py` now does:
`detection_probability()` reads this table out of `obscuration_rate.txt` and interpolates it,
anchoring each bounded bin at its midpoint in `z = ln(I_pred/noise)` and the two open end bins
at their finite edge. The scale factor `f(λ)` of the far ultraviolet is applied to `I_pred`
before the lookup although it was not applied when the curve was measured; the two are
consistent because `f` is 1 above 1000 Å and the curve is a curve for that region — the 243
predictions below 1000 Å in the measurement move no bin of it by as much as one standard
error — so correcting `z` maps a far-ultraviolet prediction onto the same curve instead of
reading it at a `z` too high by ln 10.

| I_pred / noise | predictions | D = P(recorded or masked) |
|---|---:|---:|
| below 0.1 | 3763 | 0.015 ± 0.008 |
| 0.1 – 0.32 | 1360 | 0.043 ± 0.013 |
| 0.32 – 1 | 1693 | 0.155 ± 0.013 |
| 1 – 3.2 | 1690 | 0.376 ± 0.015 |
| 3.2 – 10 | 1348 | 0.630 ± 0.015 |
| 10 – 32 | 979 | 0.812 ± 0.014 |
| 32 – 100 | 526 | 0.890 ± 0.015 |
| above 100 | 336 | 0.965 ± 0.011 |

(chance- and coverage-corrected; the plateau is 1 − ε).

**The result.** Of the 335 predicted transitions between established levels that are at least a
hundred times the noise level and fall where `c ≥ 0.95`: 317 recorded, 7 masked, 11 absent,
against a chance-match rate of 0.116. That gives

**ε = 0.035, 95 % interval 0.018 – 0.061** — about one predicted line in thirty is lost to
something too small for the coverage map to see.

And the rate is the same everywhere. Per wavelength band, at matched predicted strength:

| λ range, Å | strong predictions | absent | ε | 95 % interval |
|---|---:|---:|---:|---|
| 800 – 1200 | 66 | 1 | 0.022 | 0.001 – 0.095 |
| 1700 – 2300 | 85 | 4 | 0.054 | 0.017 – 0.122 |
| 2300 – 3000 | 18 | 0 | 0.000 | 0.000 – 0.115 |
| 3000 – 4500 | 49 | 1 | 0.022 | 0.001 – 0.093 |
| 4500 – 7000 | 51 | 0 | 0.000 | 0.000 – 0.037 |
| 7000 – 11000 | 65 | 5 | 0.072 | 0.024 – 0.154 |

A single constant fits them (likelihood ratio 7.8 on 5 degrees of freedom, p = 0.17). The
1200–1700 Å band holds only one strong prediction and is not reported. Nor does the answer
depend much on the cuts: ε is 0.054 at a 50× strength threshold, 0.035 at 100× and at 200×,
0.061 at 400× (on 83 predictions), and 0.039, 0.035, 0.035 for `c ≥ 0.90`, 0.95, 0.99.

**Do not read ε(λ) off a fit that uses the weak predictions too.** Such a fit has to assume the
detection curve has the same *shape* in every band, and it does not — the noise level is a
fitted polynomial, and a small error in it at one end shifts the whole curve there. Fitted that
way the 7000–11000 Å band comes out at 0.26 against the 0.072 measured on its own strong
predictions. The strong subsample is the measurement; everything else is the diagnostic that
justifies it.

**What ε is not.** It is a rate, not a map: it says how often a line is lost, never where. In
the likelihood it weakens the evidence of one absent branch by the factor 1 − ε = 0.965, which
is small for a single line and decisive only when several of a level's branches are absent at
once — and there it has to be combined with the **correlation length** of the obscuration, which
`tools/obscuration_length.py` measures in the next section: L = 0.067 Å, one resolution element,
so that two absences of one level are independent events unless the two predictions fall closer
together than that, which among predictions strong enough for ε to apply to them never happens.

Usage and outputs:

```bash
python tools/obscuration_rate.py                 # defaults; from inside LineClass/
python tools/obscuration_rate.py --strong 200    # a stricter definition of "certainly visible"
python tools/obscuration_rate.py --lopt LOPT_output_lines.txt   # measure a LOPT run instead
```

`obscuration_rate.txt` (the headline number, the band table, the detection curve, and the list
of absent strong predictions), `obscuration_rate.csv` (one row per prediction, with position,
strength, noise level, coverage, matching window, the three outcomes and the chance rate),
`obscuration_rate.png` (the detection curve, ε per band, and where the absent strong
predictions fall) and `obscuration_rate.log`.

### How wide is the thing that hides a line: `tools/obscuration_length.py`

**Why a rate is not enough for two absences.** ε says that about one predicted line in thirty
is lost where the coverage map calls the plate covered. For a single absent branch that is the
whole story. For two it is not. If the two absences are independent events, the missing-line
term multiplies them: ε² ≈ 1 chance in 800, and the level is in serious trouble. If one grain
flaw took both, the penalty is paid once: ε ≈ 1 chance in 30, and the level is merely
unlucky. Between those readings lies a factor of thirty in the evidence, which is enough to
decide a level on its own.

What separates them is a length. An obscuring agent — an impurity line, a scratch, a flaw in
the emulsion — occupies a stretch of plate, and two predicted lines inside the same stretch are
lost together. So the quantity to measure is

    φ(d) = the correlation coefficient of the obscuration at separation d
         = 1 if two positions d apart always share the fate of one event,
           0 if their fates are independent

and the correlation length **L** is the width of φ. Two branches closer than L are one absence;
two farther apart are two.

**How it is measured.** On exactly the population ε was measured on — the 13709 predicted
transitions between two established levels — but now pairwise: every pair of predictions closer
than 5 Å, 276 567 of them, asked whether *both* are absent, with absent meaning what it means
in `obscuration_rate.py` (no measured line within the matching window, and no recorded line
close and strong enough to have swallowed it). Joint absence has to be compared with something,
and three separate effects make nearby predictions look jointly absent when nothing hid them.

*The lines are faint.* Most predicted transitions are far below the noise, and two faint
predictions are jointly absent for no interesting reason. Each prediction therefore carries a
modelled absence probability `q`, fitted to the run: a logistic model of the outcome in the
predicted strength ln(I_pred/noise) and the wavelength, of the structural form
`P(explained) = f + (1 − f)·c·D` with `f` the prediction's own chance rate and `c` its coverage,
then raked so the modelled and observed absence rates agree in each of 40 wavelength blocks and
10 strength bins. The comparison is between the observed joint absences and Σ q_i q_j.

*The two positions see the same piece of plate by coincidence.* With a 0.022 Å window and, in
the crowded ultraviolet, a recorded line every 0.13 Å, two predictions a hundredth of an
angstrom apart are explained or not explained together about one time in six by chance alone.
That is removed the way ε removes it: the entire pair analysis is repeated at the 18 control
positions, each prediction displaced by ±6, 9, … 30 window widths — which preserves every
pair's separation, both members' predicted strengths and the local line density, and changes
only that the position is no longer one where a line was expected. The controls show the
artefact plainly, and it is confined to the window:

| separation, Å | joint absence ÷ independent, at the controls |
|---|---:|
| 0 – 0.01 | 1.176 |
| 0.01 – 0.02 | 1.165 |
| 0.02 – 0.03 | 1.128 |
| 0.03 – 0.05 | 1.091 |
| 0.05 – 0.08 | 1.029 |
| beyond 0.08 | 0.99 – 1.00 |

*The baseline is imperfect.* Any error in the model that varies slowly with wavelength makes
nearby predictions look correlated at every separation. With a baseline having no wavelength
dependence at all the raw correlation of the absences sits at 0.09 flat out to 20 Å; letting
the fit follow the wavelength on a 300 Å scale drives it to 0.00, and finer than that begins to
subtract real signal. Whatever survives is fitted as a floor, never assumed to be zero.

**The estimator, and what limits it.** Writing `A` for absent, `O` for obscured (probability ε)
and `F` for lost independently (probability d), so that `A = O or F` and `q = ε + (1 − ε)d`,

    P(A_i A_j) − q_i q_j = ε(1 − ε)·φ(d)·(1 − d_i)(1 − d_j) + O(ε²)

    φ(d) = (1 − ε)/ε · [ Σ A_iA_j − R_c(d)·Σ q_iq_j ] / Σ (1 − q_i)(1 − q_j)

with `R_c(d)` the chance factor from the controls. Two consequences. First, **no strength cut is
needed**: a pair of faint predictions has (1 − q_i)(1 − q_j) near zero and so carries almost no
weight, which is correct, because a faint prediction's absence says nothing about the plate. The
measurement is made on the strong pairs without anyone having to choose a threshold for them.
Second, the factor (1 − ε)/ε = 28 amplifies the systematic errors along with the signal: half a
percent of error in the modelled joint absence rate shows up as φ = 0.14. That, and not counting
statistics, is the limit of the method.

**The result.** Fitting the overlap of two positions with one obscuring stretch of width L,
φ(d) = (1 − d/L)₊, over all pairs closer than 1 Å with the floor free:

**L = 0.067 Å, 95 % interval 0.027 – 0.132 Å** — favoured over no correlation at all by
2 ΔNLL = 14.0 on one degree of freedom, unchanged if the pair range is 0.5, 3 or 5 Å or if ε is
set anywhere in its own interval, and 0.060 – 0.067 Å when any one of the 40 wavelength blocks
is left out. Nearly all the close pairs lie between 1000 and 2800 Å, where the effective line
width — the instrumental 0.035 Å convolved with the Doppler width — is 0.037 Å. **The
obscuration is the width of a line and no wider.** In bins, before the floor (+0.29) is
subtracted:

| separation, Å | pairs | φ | jackknife error |
|---|---:|---:|---:|
| 0 – 0.015 | 775 | 2.16 | 1.20 |
| 0.015 – 0.03 | 786 | 0.09 | 0.10 |
| 0.03 – 0.05 | 1098 | 0.35 | 0.29 |
| 0.05 – 0.08 | 1596 | 0.02 | 0.10 |
| 0.08 – 0.12 | 2255 | 0.11 | 0.12 |
| 0.12 – 0.18 | 3318 | 0.00 | 0.10 |
| 0.18 – 0.28 | 5447 | 0.22 | 0.22 |
| 0.28 – 5 | 261 292 | 0.08 | 0.04 |

The single bin below one window width carries the measurement, as it must: that is the only
separation at which two predictions are inside the same defect often enough to show.

The clustering has a second use. Its *amplitude* is proportional to ε, so fitting ε and L
together to the pairs alone measures ε again, by a statistic sharing nothing with the rate
measurement but the outcomes themselves: **ε = 0.069 (0.037 – 0.145)** against the rate's
0.035 (0.018 – 0.061). The intervals overlap and the central values differ by two, in the
direction of *more* obscuration among the close pairs — which is where more would be expected,
since almost all of them lie in the 1000–2800 Å region where the rate measurement itself puts
ε at 0.054 rather than 0.022.

The individual cases are worth looking at, and the tool lists them. Of the eleven strong
predictions that ε is measured on and that are absent, four have another absent prediction
within 0.2 Å — including the strongest absence in the whole set, 2055.265 Å at 1418 times the
noise level, which has a companion at 2055.274 Å, 0.009 Å away and 70 times the noise, absent
with it. Thirty such jointly absent pairs exist within 0.25 Å with both members at least five
times the noise, and only one of the thirty is a pair of branches of the same level.

**What it means in practice: almost always nothing, and that is the useful part.** The branches
of one level go to different lower levels, and different lower levels are hundreds or thousands
of wavenumbers apart, so the branches land far apart on the plate:

| strength cut | predictions | branch pairs | closest | median | closer than L |
|---|---:|---:|---:|---:|---:|
| I ≥ 1 × noise | 4879 | 9157 | 0.044 Å | 29 Å | 11 |
| I ≥ 10 × noise | 1840 | 3163 | 0.055 Å | 50 Å | 4 |
| I ≥ 100 × noise | 335 | 379 | 0.446 Å | 206 Å | 0 |

Among the predictions strong enough for ε to apply to them at all, the closest two branches of
any established level ever come is 0.446 Å, six times L. **Absences of one level's branches may
therefore be multiplied as independent events**, which is what the missing-line term already
does — the correction this measurement was made to supply turns out not to be needed.

The exception deserves the warning that goes with it, because where it does apply it is large.
Two predicted branches of one level closer than about 0.1 Å are one test and not two: their
joint absence costs ε and not ε², and treating them as independent overstates the case against
the level by a factor of 1/ε = 29. Four such pairs exist among predictions ten times the noise
or stronger, and `obscuration_length.txt` names them — as it happens both members of all four
were recorded, so none of them costs anything today, but the check has to be made rather than
assumed.

**What this is not.** It is not a measurement of how obscuration is distributed above an
angstrom. The systematic floor sets the sensitivity at |φ| ≈ 0.15, so a weak large-scale
modulation — a plate slightly fogged over tens of angstroms — would not be seen here, and
`coverage_map.csv` remains the only handle on structure at that scale. Nor does it separate the
agents: an impurity line, a grain flaw and a scratch all produce the same statistic, and the
measured width is close enough to the resolution element that the data cannot tell an emulsion
defect from a line of some other species too weak to have been measured.

Usage and outputs:

```bash
python tools/obscuration_length.py                # defaults; from inside LineClass/
python tools/obscuration_length.py --dfit 0.5     # fit the length on closer pairs only
python tools/obscuration_length.py --epsilon 0.05 # impose a rate instead of measuring it
```

`obscuration_length.txt` (the fitted length, φ(d) in bins, the control check, the branch
separations and the jointly absent close pairs), `obscuration_length.csv` (one row per jointly
absent close pair), `obscuration_length.png` (φ(d) with the fit, the coincidence artefact, and
the branch separations against L) and `obscuration_length.log`.

## Transitions missing from `Icalc.xlsx`: the censoring correction

### What the absence of a transition means

`Icalc.xlsx` holds the calculated transitions: for each pair of energy levels it gives
`gA` — the statistical weight of the upper level multiplied by the probability, per
second, that the atom makes that transition — together with the predicted intensity
`Icalc` and the uncertainty `u%gA` of `gA` in percent. `gA` measures how strong the line
is expected to be.

The file is **not** a list of every transition that exists. Cowan's atomic-structure
codes, which produced it, printed a transition only when `gA` reached **1000 s⁻¹**. So a
level pair that does not appear in the file is not a transition of unknown strength: it
is one **known to be weaker than the cutoff**. Statisticians call data of this kind
*left-censored* — the value was measured, found to be below a threshold, and only the
fact of being below it was recorded.

Checked against the data (593 levels of the `Wyart2000` sheet):

|                                                                                  | pairs         |
|----------------------------------------------------------------------------------|---------------|
| electric-dipole-allowed level pairs within the wavenumber range of the file       | 31646         |
| present in `Icalc.xlsx`                                                           | 30206 (95.4%) |
| absent, i.e. `gA` < 1000 s⁻¹                                                      | 1442 (4.6%)   |
| present but *not* dipole-allowed                                                  | 0             |
| rows dropped because a level id is unknown                                        | 0             |

Every pair in the file is allowed and every allowed pair is either in the file or below
the cutoff, so the censoring reading holds exactly, and it applies to 4.6% of the
candidate pairs.

Until 2026-08 the pipeline treated such a pair as having *no* predicted intensity, which
meant that both intensity tests — the relative-intensity filter on a blended line and the
z-test of predicted against observed intensity on a line with a single candidate — were
simply skipped for it. An identification resting on a censored pair therefore escaped the
strongest evidence against it: the theory says this line should be far too weak to see,
and it was seen. The censoring correction gives such a pair the intensity implied by a
`gA` just below the cutoff, so that its predicted weakness counts against the
identification instead of exempting it.

### The imputed value

Write `u_ln = ln(1 + u%gA/100)` for the uncertainty of `gA` on a logarithmic scale (the
`u_ln` column of the file). A logarithmic scale is the natural one because the uncertainty
of a calculated `gA` is a factor, not an amount: `u_ln = 0.69` means "uncertain by about a
factor of two", in whatever decade `gA` happens to lie.

1. **`u_ln_missing`** = root mean square of `u_ln` over the decade just above the cutoff,
   1000 ≤ `gA` ≤ 10000 s⁻¹. Measured: **1.7383**. (`u_ln` grows as `gA` falls — the weak
   transitions are the badly calculated ones — so the uncertainty appropriate just below
   the cutoff is the one seen just above it.)
2. **`gA_missing` = cutoff · exp(−`u_ln_missing`) = 176 s⁻¹**, i.e. the value whose
   **upper one-standard-deviation bound sits exactly on the cutoff**,
   `gA_missing · exp(u_ln_missing) = 1000 s⁻¹`. The censored transition is thereby made as
   strong as it could be while still being consistent with having been left out of the
   file. A larger value would contradict the censoring; a smaller one would assert
   knowledge that is not there and punish the identification harder than the data warrant.

   > **Not** "gA minus its uncertainty is zero." That reading appeared in early notes and
   > is wrong. With a logarithmic uncertainty the linear standard deviation is
   > `σ = gA·(e^u_ln − 1)`, so `gA − σ = gA·(2 − e^u_ln)` vanishes only at
   > `u_ln = ln 2 = 0.693`, not at the 1.738 measured here — and on a log-normally
   > distributed quantity zero is unreachable at any finite number of standard deviations.
   > The lower one-standard-deviation bound of the imputed value is
   > `cutoff·exp(−2·u_ln_missing)` = 31 s⁻¹, 3% of the cutoff.

3. The intensity follows from the relation every row of the file obeys (see
   [The intensity model](#the-intensity-model-boltzmann-plot) below):

   ```
   Icalc = C · gA · (rwn / 1e8) · exp(−Eup / kT) ,   C = 251.606 ,  kT = 12375.4 cm⁻¹
   ```

   with `Eup` the energy of the upper level and `rwn` the Ritz wavenumber, both in cm⁻¹.
   `kT` is an effective excitation temperature written as an energy in the same unit.
   `gA_imputation.check_intensity_model` re-fits `C` and `kT` on each run and warns if the
   file disagrees with the configured values by more than 1%.

`python gA_imputation.py --report` prints the whole calibration: the decade-by-decade
profile of `u_ln`, the estimate in force, and the model check.

### Switching the policy

`missing_gA.policy` in `lineclass_config.toml` selects the treatment, and
`--missing-gA {none,impute}` overrides it on the command line of `classify_lines.py`,
`decoy_mc.py` and `chance_mc.py`:

| value                | meaning                                                                             |
|----------------------|-------------------------------------------------------------------------------------|
| `"impute"` (default) | a censored pair gets the imputed intensity above and faces both intensity tests      |
| `"none"`             | a censored pair has no predicted intensity and the intensity tests are skipped for it (the behaviour before 2026-08) |

The switch is on the Monte-Carlo drivers as well because a false-positive calibration made
under one policy must never be used to judge a run made under the other: the decoy runs are
the yardstick, and the yardstick has to be marked in the same units as the thing measured.

The imputed intensity is computed at the point of use (`classify_lines.imputed_values`),
not when the file is read, so a decoy level's transitions are imputed exactly when the real
transitions they shadow would be. That is what keeps the false-positive calibration fair.

### What the correction does to the result

Full comparison, evidence and verdict: `PLAN_missing_gA.md`, "Stage 5 results" (a working
document kept outside version control). Both runs
are archived under `baseline/policy_none/` and `baseline/policy_impute/` and are reproduced
exactly by the current code; `tools/analyze_policy_diff.py OLD NEW [--list]` regenerates the
comparison.

Accepted identifications fall from 3336 to 3300 — **36 lost, 2 gained, 1.1% of the total**;
33 of the 36 rest on a censored pair. The identifications removed behave like false ones on
three independent measures:

|                                                                       | `none`         | `impute`      |
|-----------------------------------------------------------------------|----------------|---------------|
| accepted **decoy** identifications per run (false by construction)     | 91.2           | **50.0**      |
| expected false level confirmations at the best criterion               | 2.12           | **1.00**      |
| tested levels validated at that criterion                              | 50/123 (40.6%) | 46/113 (40.7%)|
| rms level shift &#124;dE&#124; over the 192 new levels                 | 0.1318         | **0.0941** cm⁻¹ |
| rms &#124;dE&#124; over the 64 levels that lost an identification      | 0.1837         | **0.0845** cm⁻¹ |
| new levels with spurious probability `p_spur` < 0.1                    | 49             | **57**        |

The policy removes 45% of the demonstrably false identifications at a cost of 1.1% of all
of them — about 40 false for each real one — while the same proportion of tested levels
survives validation and the expected contamination among them halves. The levels that lose
an identification are precisely the ones that were wandering: their rms energy shift was
twice the average before the correction and is below it after, which is the signature of
removing a false line rather than a true one.

One number carries the argument without depending on the imputation recipe at all: under
`"none"`, identifications resting on a censored pair are **25.8%** of what the decoys accept
(188 of 730 pooled decoy identifications) but only **1.2%** of what the real run accepts (40
of 3336). They are 21 times over-represented among identifications that are false by
construction, so exempting them from the intensity tests was letting through the very
population the decoys mark as spurious.

Nine new levels lost every identification they had. All nine were supported *only* by
identifications on censored pairs, and all nine already carried `p_spur` ≥ 0.40 under the
old policy — not one was a level the old run considered established. One level already in the
NIST ASD, 059003.000304, also lost both of its identifications; both were new proposals on
censored pairs, so its published status is untouched. The sole-candidate
z-thresholds were examined and deliberately **left unchanged**: for the 33 rejected imputed
identifications the z of the intensity test runs from 3.34 to 5.88 (median 4.79) against a
rejection threshold of 1.25, so every one of these rejections would stand at any threshold
up to 3.3 and the thresholds have no leverage on the outcome.

One caution when comparing the two archived runs: `p_spur` is **not** comparable one-to-one
between them, because each run carries its own decoy calibration and the `"impute"`
calibration is the tighter one. Two levels cross the 0.1 mark while keeping exactly the same
identifications and the same energy shift.

---

## Validation of the classification (`decoy_mc.py`, `level_shifts.py`, `level_interchange.py`)

### The problem

With the current inputs the pipeline accepts 5155 identifications: 3456 previously published and 1699 new ones. The new identifications are what re-establishes the levels of Wyart's list that were never confirmed in print. Two questions must be answered before those levels can be adopted: *how often would this pipeline confirm a level that is not real*, and *which individual levels are doubtful*. The validation suite answers both. All numbers quoted below refer to the current input data and are reproduced by running the suite.

### Which levels are tested

A level is **tested** ("new\*") if, in the real run, it is supported by more new than old identifications (`n_new > n_old`): 214 of the 593 levels. This set comprises the 192 levels absent from the ASD compilation (blank `E_ASD` in the levels file) plus 22 levels that are formally in ASD but whose energies rest mainly on the new identifications. The remaining **379 old levels**, anchored by their previously published identifications, serve as the reference of genuine behavior.

### The quantities used to judge a level

After the run, each level is characterized by:

- **n** = `n_old + n_new` — the number of accepted lines supporting the level;
- **ΔE** (`dE`) = `E_final − E_input` — how far the level's energy moved during the run, where `E_input` is the value in the input file (Wyart's) and `E_final` is the optimized value;
- **d** = `ΔE · n` — the support-normalized shift.

Why ΔE discriminates: Wyart derived each new level from real lines of this same line list. If the pipeline re-finds those lines, the least-squares optimization must return the level's energy to almost exactly Wyart's value, because those lines tie it to old levels that barely move. If instead the accepted lines are accidental wavenumber coincidences, they drag the level toward unrelated positions. With the current inputs, old levels move by a median |ΔE| of 0.015 cm⁻¹ (90% below 0.070, 95% below 0.091, 99% below 0.188, maximum 0.379 — these percentiles are the candidate thresholds `T` used below), while levels supported by false matches move by ~0.2 cm⁻¹ (median) with a long tail to above 1 cm⁻¹.

### The decoy experiment (`decoy_mc.py`)

To measure how often a level that is certainly *not* real would pass any acceptance rule, the classification is repeated with one **decoy** added per tested level. A decoy (created by `add_decoy_levels()` in `classify_lines.py`, activated by `main(decoy_shift=δ, decoy_ids=…)`) is an exact copy of its level — same J, same parity, same set of possible transitions, same calculated intensities — with the energy displaced by δ, the sign alternating with parity so that the displacement can never cancel out of any transition wavenumber (transitions connect opposite parities; decoy–real pairs are displaced by δ, decoy–decoy pairs by 2δ). δ exceeds the widest matching window the pipeline ever uses (5.5 × max σ = 6.3 cm⁻¹), so no decoy transition can land on a true line of its original; eight values between 6.9 and 11.3 cm⁻¹ are used. Everything else about the run is fully realistic — real lines, legacy identifications, real levels all in place, decoys competing in every step — so every line accepted for a decoy is a false match obtained under real-run conditions. With one decoy per tested level, the mean number of decoys passing a rule equals the number of false confirmations to expect **if all tested levels were fake** (a deliberate worst case).

Findings worth stating: mere support proves nothing — 89% of decoys collect at least one accepted line and 49% collect four or more; it is the energy shift and the intensity pattern that separate fake from real. The decoy competition perturbs the real solution only mildly (mean accepted count 5152 vs 5155; 73% of level×run combinations keep their support exactly, 24% change by ±1 line).

`decoy_mc.py` writes `decoy_mc_levels.csv` (per decoy and run: support and ΔE), `decoy_mc_real_levels.csv` (the real levels in each decoy run, for the perturbation check), `decoy_mc_lines.csv` (the accepted decoy lines), and `decoy_mc_summary.csv` — each with an `.xlsx` twin. `--smoke` runs one δ with 2 cycles as a plumbing test.

### The shared utilities and the shifted-wavenumber cross-check (`chance_mc.py`)

`chance_mc.py` now serves mainly as a library for the other two scripts: `read_input_levels()` (input energies and ASD status of every level), `per_level_table()` (per-level support counts, final energies and ΔE from a run output), and `save_table()` (Excel-friendly output writing, see below). Run directly, it performs an independent cross-check: all observed wavenumbers are shifted by a constant (again 6.9–11.3 cm⁻¹, beyond the widest matching window), so that no line can coincide with any true transition and, once more, everything accepted is a false match. Because the shift destroys the true matches, the whole line list is left free for false matching — these rates are therefore upper bounds. The two experiments measure the same rates by very different means and agree within ~25% at every setting, which shows that neither result is an artifact of how it was measured. `level_shifts.py` includes the cross-check section when `chance_mc_levels.csv` is present and skips it otherwise.

### The validation report (`level_shifts.py`)

Running `python level_shifts.py` produces the full analysis:

1. **Old-level calibration** — the |ΔE| percentiles quoted above, and the per-support behavior of d.
2. **Criterion grids** — for every combination of `N` (minimum number of supporting lines) and `T` (largest allowed |ΔE|), the number of tested levels passing the cut `n ≥ N and |ΔE| ≤ T` against the expected number of false passes (decoy calibration, and the shifted-wavenumber upper bound).
3. **Intensity-pattern check** — the modern form of the classical "square-array" argument: a real level must reproduce the *pattern* of theoretically predicted intensities. For every level, all predicted transitions to/from it (from `Icalc.xlsx`, both partners known, Ritz wavenumber inside the observed range) are sorted by predicted intensity, and three scores are computed. Two classes of predictions are first excluded, because a transition that could not have been seen must not count as missing: (a) predictions whose Ritz wavelength falls where **the record is blind** — where the measured coverage `c(λ)` of `coverage_map.csv` is below `COVERAGE_MIN` = 0.5, the same threshold `tools/coverage_map.py` uses to list its blind stretches. This replaced the two hand-entered intervals (`COVERAGE_GAPS_A`: 1522.49–1529.85 and 2103.46–2107.92 Å vacuum, the seams of the intensity-calibration regions), which were far too narrow for this purpose: the map finds ten blind stretches totalling 1078 Å, and the two longest run out to 1663.82 and 2190.91 Å rather than stopping at 1529.85 and 2107.92 Å. The old pair survives only as the fallback used when `coverage_map.csv` is absent; (b) when the intensity-calibration file `intensity_correction_functions.txt` is present (piecewise polynomials `P(λ)` in vacuum wavelength converting Sugar's plate intensity to the uniform linear scale via `I_linear = 1000·I_Sugar·exp(P(λ))`), predictions whose Sugar-scale intensity `I_pred/(1000·exp(P(λ)))` falls below the noise level (Sugar intensity 1). A third correction enters that comparison. `Icalc` is meant to be on the same linear scale as the observed intensities, and from 1000 Å upward it is — the median of `I_obs/I_pred` over the accepted identifications is within about 30 % of 1. Below 900 Å it is not: the median falls to 0.09 (83 accepted lines), and the ratio of the band totals — all observed intensities over all predicted ones, which does not depend on which lines were classified — agrees at 0.16. A prediction there is therefore about ten times fainter on the plate than its `I_pred` says. `intensity_scale_bias` measures that factor band by band (`SCALE_BIAS_BANDS_A`, at least `SCALE_BIAS_MIN_N` = 10 accepted lines, correction applied only where the factor is below `SCALE_BIAS_MAX_FACTOR` = 0.5) together with the scatter of `ln(I_obs/I_pred)`. A blended line enters this comparison with the share of its measured intensity the component being compared carries — its `BF` — for the reason given under `level_interchange.py` below: the whole feature measured against one component's prediction would put a blended line above the scale by the reciprocal of its share. Splitting the blends moves the measured factors of the run of 2026-09-05 from 0.482 to **0.414** below 900 Å and from 0.735 to **0.567** in 900–1000 Å, and tightens the scatter elsewhere from sd 1.32 to **1.25**. The corrected intensity `I_pred·f(λ)` is what has to clear the noise level. This is a fault of the calculated intensities, not of the plates: Sugar's plates reach his intensity 1 everywhere, and Pr III lines of intensity 1 are recorded and classified down to the short-wavelength edge at 821.9 Å. With the current inputs 3813 predictions are excluded as falling in the blind stretches (111 under the old hand-entered pair) and 18056 below the noise level, of 29260 loaded — most theoretically predicted transitions of these high levels are simply too faint for Sugar's plates. Accepted lines whose predictions are below the noise level stay in the support count n (the decoys are treated identically, so the calibration stays fair); the report column `n_acc_below_noise` makes them visible.
   - `top10_found` — how many of the 10 strongest predictions are among the accepted lines (`n_top10` = how many were available);
   - `pattern_C` — intensity-weighted completeness: the sum of predicted intensities over the accepted transitions divided by the sum over all predicted observable ones;
   - `pattern_V` — C divided by the largest C achievable with the level's number of matched lines. V = 1 means the accepted lines are exactly the strongest predictions; V = 0 means predictions exist but none was found. Unlike C, V does not punish a level for accepted lines that theory does not cover, and it is nearly independent of n.
   - `top1` — the fate of the level's **strongest** predicted observable transition on its own: `found` (it is among the accepted lines), `masked` (absent, but an observed line hides it — see the masking test in item 5), `faint` (absent, nothing hides it, but the line would not reliably have been recorded anyway — see below), `missing` (absent and unexplained), or blank (no observable prediction). This is scored separately because `pattern_V` cannot catch it: V asks whether the *accepted* lines are the strongest of the predictions, so a level whose two accepted lines are ranked 2 and 3 reaches V ≈ 1 with rank 1 missing. A missing, unmasked strongest branch is the most direct evidence there is against a level — but only where the line would have been seen, and how surely it would have been seen is now measured rather than modelled. `observation_probability` returns **P_obs = c(λ)·D(z)**: the coverage of the plate at that wavelength times the measured detection curve at that strength, `z = ln(I_pred·f(λ) / noise level)`. `D` comes from `obscuration_rate.txt` — the observed fraction of predictions of each strength that were recorded where the plate is open — and has a ceiling at `1 − ε` = 0.965, because about one predicted line in thirty is lost to a grain flaw or a single impurity line wherever it falls. Below `DETECT_CONFIDENCE` = 0.9 the absence carries no information and `top1` reads `faint` instead of `missing`. Two report columns expose the quantity: `top1_pobs`, the P_obs of the strongest branch, which is what the threshold is applied to, and `n_obs_expected`, the sum of P_obs over all of a level's observable predictions — how many of its branches should have appeared at all, which is the number a missing-line term has to weigh the accepted count against. It is well below `n_pred`: the median ratio over the 594 levels is 0.60, because most predicted branches sit between one and ten times the noise, where `D` stands at 0.4 to 0.8.

     This replaced a log-normal centred on `I_pred·f(λ)` with the measured scatter of `ln(I_obs/I_pred)`, integrated above the noise level. That model was fitted to the **accepted** lines, whose lower tail the noise has already removed, so it was optimistic exactly where it mattered: it put the recording probability at 0.99 by ten times the noise level, where the measured rate is 0.81, and it reached 0.9 at about five times the noise where the measured curve does not reach it until sixty. The bar therefore rose by a factor of twelve in strength, and the effect is large. Running the pipeline twice on identical inputs, once with `coverage_map.csv` and `obscuration_rate.txt` in place and once with them hidden so that the old model is used, `top1 = missing` falls over all 594 levels from **27 to 8**: 17 of the 27 become `faint` — their strongest branch would not surely have been seen — and 3 become `found`, their previously strongest branch having fallen in a blind stretch the coverage map now excludes, so that the next one down, which is accepted, takes its place. One goes the other way, `found` to `missing`, for the same re-ranking reason. Among the 210 tested levels the count goes 13 → 2, and `question_status = retained` goes 13 → 3: ten levels leave the questionable examination altogether. The eight that remain missing are strong cases, `top1_pobs` between 0.90 and the 0.965 ceiling.

     Two side effects are worth naming. The shrunken observable set (7391 predictions against 7776) moves `pattern_V` for 110 levels, by up to 0.53 — the denominator of the pattern score no longer contains branches that were never on a plate. `p_spur` barely moves at all (4 levels by more than 0.005, at most 0.15), because the V distributions it is weighed against are refitted on the same run and absorb the shift. And `pattern_C`/`pattern_V` still treat every surviving prediction as equally observable: the coverage function enters them only through the hard `c < 0.5` cut, not as a weight. Weighting the pattern denominator by `P_obs` is the natural next step and is left for the likelihood proper, since it would change what `C` means — the share of predicted intensity found becomes the share of *expected observable* intensity found, and `C` would no longer be bounded by 1.

     The threshold is the one knob: at `DETECT_CONFIDENCE` = 0.85 thirteen levels read `missing`, at 0.8 twenty, at 0.7 twenty-five. It cannot be raised much above 0.9, since `D` can never exceed 0.965 and a threshold above that would make the test unpassable at any brightness.
   
   Calibration: genuine (old) levels have a median C of 0.91 and median V of 0.94; decoys have 0.03 and 0.06. The information here is genuinely new — the weeding checks each line's intensity individually but never penalizes a level for strong predicted lines that are *absent*. Two known weaknesses: a single grossly wrong predicted rate can wreck C and V of a real level (theory faults of this kind are common), and levels whose transitions are not covered by the theory get no pattern verdict at all.
4. **Per-level spurious probabilities.** For each tested level, two independent pieces of evidence are weighed between "genuine" and "spurious" (= behaves like a decoy): the energy shift (|d| modeled as a log-normal distribution whose median depends on n, fitted separately to the old levels and to the decoys) and the intensity pattern (the probability of the level's V range under each hypothesis, measured on the old levels and on the decoys). The prior probability that a tested level with support n is spurious is `S·q(n)/m(n)`, where `q(n)` is the decoy probability of ending with support n and `m(n)` is the observed number of tested levels with support n; the single unknown S — the total number of spurious levels — is estimated by maximizing the likelihood. Three posterior probabilities go into the report:
   - `p_spur_decoys` — from the energy shift alone;
   - `p_spur_pattern` — from the intensity pattern alone;
   - `p_spur` — both folded together (the final value).
   
   With the current inputs (converged model, 208 tested levels, run of 2026-09-05): S = 0.0 [0.0–10.8 at 95% confidence] from the energy shifts alone (26.5 [17.0–42.8] under the alternative small-n assumption — poorly determined), 2.8 [0.0–9.0] from the pattern alone, and **2.0 [0.0–7.2] folded** (4.2 in the robustness variant) — the pattern evidence removes the model sensitivity. The sum of `p_spur` over the tested levels reproduces S, and the expected number of false confirmations under any cut equals the sum of `p_spur` over the levels passing it. Marking levels with `p_spur ≥ 0.1` as questionable flags 2 levels, both of them supported by a single line; the 206 unmarked levels then carry a summed probability of only 0.45 expected spurious levels, while the marked pair is expected to contain ~1.55 spurious and ~0.45 genuine members (a "questionable" mark is a caution, not a verdict).
5. **Masking check (adjudication of the questionable marks).** Two things put a tested level up for examination: the probabilities question it (`p_spur ≥ 0.1`), or its strongest predicted transition is absent from the accepted lines and would have been recorded had it been there (`top1` reads `missing`; a branch too faint to have been recorded reads `faint` and does not put the level up for examination). A strong predicted transition may be absent from the accepted lines simply because a nearby **stronger** observed line hides it — on the photographic plates, a weaker line close to a much stronger one cannot be measured. Every such level is therefore re-examined. The test uses the **effective line width** (full width at half maximum): the instrumental width `INSTR_FWHM_A` (0.035 Å, estimated from the closest measured line pairs; the spectrograph resolution is nearly constant in wavelength) convolved with the Doppler width, which grows with wavelength. For each missing prediction among the level's strongest ten, an observed line hides it in either of two regimes: **within one effective width** the two lines are not resolved, so a line of at least `MASK_STRENGTH` (0.5) times the predicted intensity absorbs it; **beyond one width**, a Gaussian profile of the effective width is applied to the observed line, and the prediction counts as hidden if the profile intensity at its position still exceeds the predicted intensity — so the masking reach of a much stronger line extends into its wings, growing with the intensity ratio. (The manually verified masking cases show separations up to 0.029 Å, inside one width.) If at least `MASK_CLEAR_FRACTION` (90%) of the *missing predicted intensity* is hidden — i.e. the dominant missing predictions are masked, while a small residue of weak ones is tolerated (the theory is least reliable for weak transitions) — the bad pattern score carries no evidence, and the questionable mark is cleared — provided the energy-shift evidence alone is not suspect (`p_spur_decoys < 0.1`) **and** the strongest branch is not itself among the unexplained absences. A missing, unmasked strongest branch keeps the mark whatever the probabilities say. The outcome goes into two report columns: `question_status` (`removed` = mark cleared, `retained` = mark kept) and `reason` (e.g. "strong missing transitions are masked by nearby stronger lines", "no reason to clear", "no theoretical predictions (pattern check not applicable)", "masking found, but the energy shift alone remains suspect"). With the current inputs (2026-09-05), 15 levels are examined (2 flagged by `p_spur`, 14 by a missing strongest branch, 059003.000589 by both); 1 is cleared (059003.000615, 96% of its missing predicted intensity masked), leaving **14** questionable, all of them because their strongest predicted transition is absent and nothing hides it.
6. **Report file** — `level_shift_report.csv/.xlsx`, one row per level: `level_id`, `is_new_level` (absent from ASD), `new_star` (tested), `note`, `J`, `parity`, `E_input`, `E_final`, `dE`, `n_old`, `n_new`, `n_tot`, `d`, `n_pred`, `n_top10`, `top10_found`, `top1`, `pattern_C`, `pattern_V`, `n_acc_below_noise`, `p_spur_decoys`, `p_spur_pattern`, `p_spur`, `question_status`, `reason` (the probabilities and the adjudication columns are filled only where applicable).

**Single-level inspection:** `python level_shifts.py --detail <level_id> …` prints one level's predicted transitions (strongest first) with the fate of each in the real run — accepted (with the observed line), rejected (with the weeding note), or not matched by any observed line — plus the accepted lines that have no theoretical prediction. This is the working tool for judging individual questionable levels.

### Revised identifications: building the report from LOPT (`lopt_lines.py`)

Identifications do not stay as the pipeline left them. They are revised by hand — in
IDEN2, or by editing the LOPT input file — and **LOPT** (the least-squares level
optimizer) is then re-run on the revised list. `--lopt` builds the same validation
report from LOPT's own line-output file, so a hand-revised run can be validated with
the identical machinery:

```
python level_shifts.py --lopt LOPT_output_lines_revised.txt \
                       --e-input revised_level_energies.csv
```

The report is named after the LOPT file (`LOPT_output_lines_revised.txt` →
`level_shift_report_revised.csv`) so it never overwrites the pipeline's own; `--report`
overrides the name.

**What `lopt_lines.py` recovers, and from where.** LOPT's line output is tab-separated
with one header row; which columns it holds depends on the print options of the
parameter file, so they are looked up by name, and a number LOPT does not give prints
as `_`.

| the report needs | where it comes from |
|---|---|
| `accepted` | `Weight > 0`. Rows flagged `P` carry weight 0 — predicted (Ritz) positions LOPT was asked to print but not to fit. Rows flagged `c` are the components of a resolved blend and carry their branching fraction; they are accepted. |
| `low_E`, `upp_E` | re-solved from the accepted lines (see below), or LOPT's own energies with `--energies lopt` |
| `new` | not in LOPT's output. An identification is *old* when the same pair of levels stands in the `id1`/`id2` columns of the line workbook — the published identifications — against a line at the same wavenumber. Level pairs are unique across the line list, so the pair alone is a safe key and the wavenumber is only a guard. |
| `wn_obs`, `obs_intens` | LOPT prints both, but rounded to the precision of the line's uncertainty; the exact values are taken from the line workbook by nearest-wavenumber match within 0.05 cm⁻¹. |

Columns LOPT cannot supply — a candidate's calculated intensity, its grade, the weeding
notes — are written as blanks, and `--detail` says so.

**Why the energies are re-solved rather than taken from LOPT.** The report judges a
tested level by its energy shift ΔE against the same quantity measured on levels that
are false by construction: the decoys of `decoy_mc.py` and the chance levels of
`chance_mc.py`. Those reference populations are made by re-running the pipeline, so
their `E_final` comes from the pipeline's optimizer. Taking `E_final` for the tested
levels from LOPT instead would measure the two sides of the comparison with different
instruments: LOPT tunes single-line levels, treats blends by its centroid model and
rounds its output, and its energies differ from the pipeline's by up to 0.2 cm⁻¹ — the
same size as the ΔE under test, and enough to widen the fitted scatter of the genuine
population from s = 1.34 to s = 1.98, which quietly makes every level look more genuine.

So `--lopt` takes the **accepted set** from LOPT and re-solves the energies from it with
the pipeline's own weighted least squares (`lopt_lines.refit_energies`): each accepted
transition is one observation equation `E(upper) − E(lower) = wn_obs` weighted by
`BF²/u²`, the ground level is held fixed, and one level of any group not connected to it
is anchored at its input energy. Fed the pipeline's own identifications, the refit
reproduces `line_classifications.csv` to a median of 1.3·10⁻⁶ cm⁻¹ and a maximum of
0.0012 cm⁻¹, and the resulting report reproduces every integer column of
`level_shift_report.csv` exactly and `p_spur` to 0.001. A `--lopt` report and the
ordinary one therefore differ only in the thing being studied — which lines are assigned
to which levels — which is what makes a before/after comparison of a revision fair.

Add `--energies lopt --lopt-levels LOPT_output_levels_revised.txt` to report LOPT's own
optimized energies instead. Those are the energies to publish; they are not comparable
with the calibrations. The level-output file is worth passing because its four decimals do not
depend on how much precision LOPT chose to print for each level, unlike the `E1`/`E2`
columns of the line output; for `LOPT_output_lines_revised.txt` the two agree exactly on
all 593 levels.

**`--e-input`: levels a revision deliberately moved.** ΔE asks whether the optimizer
stays where the identifications put the level. For a level that was *re-positioned*, the
adopted-level workbook is the wrong reference and ΔE would report the size of the move
rather than test it. Give the new starting energies in a two-column CSV
(`level_id,E_input`); the levels listed there are reported on the console, the rest keep
the workbook value. Note that ΔE is then not an independent test for those levels
either, if the new `E_input` was itself computed from the very lines the revision
assigns to them — their evidence is the intensity pattern and the residual spread.

**What the calibrations do and do not need.** `decoy_mc.py` and `chance_mc.py` take no
LOPT input and need none: they insert fake levels into the pipeline, or shift every
wavenumber, and measure what the pipeline does with them. Their output is a *null*
distribution — how a level that is false by construction behaves — and it is the real
side of the comparison that has to be brought into the same units, which is what the
refit does. Both are deterministic given the pipeline, so re-running them after a
revision of the identifications reproduces the same files. The one residual
approximation is that the null was measured in an environment where the pipeline's own
assignments consumed the observed lines; a revision that changes a handful of the ~4900
accepted lines changes that environment by well under a percent.

### Interchanged identities: `level_interchange.py`

Every level here carries two things established separately. Its **energy** is a number in
cm⁻¹ measured from the observed lines and fixed by the least-squares optimization to about
a hundredth of a cm⁻¹. Its **theoretical identity** — which level of the Cowan-code
calculation it is — is a label (a configuration such as `4f².5d` and a term such as `³H`),
and that label carries with it the whole set of calculated transition probabilities *gA*
of the level, i.e. the predicted intensity of every line it can emit. The predictions live
in `Icalc.xlsx` keyed by level id; the label itself lives in IDEN2's `enlev.dat`.

The identity is assigned by matching the observed energy to a calculated one, and the
calculated energies of Pr III are not accurate. **How inaccurate depends strongly on the
configuration.** Over the levels whose energies are experimentally known, the rms of
`E_obs − E_calc` runs from 22 cm⁻¹ for `4f².5g` to 404 cm⁻¹ for `4f.5d.6p`:

| config | levels | rms (cm⁻¹) | max ǀO−Cǀ | | config | levels | rms (cm⁻¹) | max ǀO−Cǀ |
|---|---:|---:|---:|---|---|---:|---:|---:|
| `fd6p` | 9 | 404.2 | 827.4 | | `f27d` | 5 | 87.9 | 122.3 |
| `f5d2` | 66 | 207.6 | 623.8 | | `f28s` | 8 | 74.5 | 125.6 |
| `f5d6s` | 19 | 184.1 | 607.5 | | `f26d` | 45 | 61.0 | 152.5 |
| `f25d` | 105 | 147.8 | 532.8 | | `f27s` | 12 | 58.6 | 88.0 |
| `f26p` | 62 | 114.3 | 287.2 | | `f26f` | 23 | 45.4 | 93.2 |
| `f25f` | 118 | 104.1 | 362.4 | | `f27p` | 15 | 32.2 | 68.6 |
| `f26s` | 22 | 98.6 | 207.5 | | `f25g` | 45 | 22.2 | 60.0 |
| `4f3` | 40 | 97.3 | 234.0 | | | | | |

(`enlev.dat`'s own abbreviations: `f25g` = 4f².5g, `fd6p` = 4f.5d.6p, `f5d6s` = 4f.5d.6s.)

Wherever two calculated levels of the **same parity** and the **same J** lie closer
together than that error, the energy match cannot tell which is which, and the two
identities may have been **interchanged**: level A wears B's label and B's *gA* values, and
B wears A's.

**Why no other test sees it.** Under an interchange the energies stay right — they were
measured from the lines, and no line moves — so ΔE says nothing at all. What breaks is the
intensities: each level is judged against the branching pattern of the other.
`level_shifts.py` sees that only through `pattern_V`, one ingredient of a folded
probability, and a level can keep a small `p_spur` with a bad pattern.

#### The test

For a candidate pair (A, B) nothing moves; only which set of *gA* values belongs to which
level is exchanged. Under the two hypotheses — **assigned** (A keeps A's *gA*, B keeps B's)
and **swapped** (A is scored against B's *gA*, B against A's) — every accepted line of both
levels is given the intensity that hypothesis predicts for it,

```
I_calc = C · gA(identity, partner) · (wn/1e8) · exp(−E_up/kT)
```

with the observed wavenumber and the observed upper-level energy of the line — both
untouched by the swap — and the pipeline's own `C` and `kT`. A level pair absent from
`Icalc.xlsx` gets the imputed *gA* of `gA_imputation.py` (175.8 s⁻¹ here), so a hypothesis
that predicts nothing where a strong line is observed is penalized rather than excused.

**Blended lines are split, and split under each hypothesis separately.** One measured line
is often produced by several predicted transitions falling together, and the classification
then assigns all of them to it; every one of those rows carries the same measured
intensity, that of the whole feature. What is compared with *one* transition's prediction
is that transition's share of the feature,

```
share = I_calc(this component) / Σ I_calc(all accepted components of the line)
```

which is the branching fraction the pipeline already stores in the `BF` column of
`line_classifications.csv` and squares to weight blend components in the least squares.
Giving a component the whole feature instead would overstate its observed intensity by
`1/share` — a factor of 25 for a component calculated to carry 4 % of it. In the run of
2026-09-05, **906 of the 4939 accepted lines are blend components**, 478 of them with a
share below ½ and 211 below ⅕, and the levels of the two flagged pairs are among the worst
affected: 8 of the 14 lines of `059003.000483` are blend components (shares down to 0.09)
and 5 of the 13 of `059003.000398` (down to 0.03).

The share is recomputed **inside each hypothesis**, not read from the file: the stored `BF`
was worked out under the assigned identities, and a hypothesis that makes a component
fainter also makes it take a smaller piece of the blend. Using the stored value on both
sides would divide the two hypotheses by the same number and cancel out. So the component
whose identity the swap changes gets its `I_calc` from that hypothesis's *gA* while its
neighbours in the blend, untouched by the swap, keep theirs. (The band scale factor `f(λ)`
below is common to a blend and cancels in the ratio.) A component of a blend one of whose
members has no calculated intensity at all is given an even split, `1/n`, which is the same
fallback `classify_lines.calc_weights` uses. See `blend_share()` and
`accepted_lines_by_level()`.

Two numbers compare the hypotheses.

- **The slope.** `ln(I_obs)` fitted against `ln(I_calc)` by least squares, per level. A
  correct identity gives a slope near 1; a wrong one gives a slope near 0, because the
  predicted pattern then carries no information about the observed one. All four slopes
  are reported (each level under each hypothesis).
- **The likelihood ratio**, which decides. The slope ignores the *scale* — predictions ten
  times too large can still have slope 1. Under a correct identity the residual
  `r = ln(I_obs) − ln(I_calc · f(λ))` is centred on zero with a known spread, where `f` and
  the spread are the band scale factor and the sd of `ln(I_obs/I_pred)` that
  `level_shifts.intensity_scale_bias` measures on the same run (this is what removes the
  far-ultraviolet scale error of the calculated intensities). Summing `−r²/(2·sd²)` over
  both levels' lines under each hypothesis gives

  ```
  lnR = ln L(assigned) − ln L(swapped)
  ```

  positive when the assignment fits better, negative when the swap does. It is a difference
  of two fits to the *same* observed intensities, so the plate calibration, the Boltzmann
  factor and the observed scale all cancel. `lnR_A` and `lnR_B` split it into the two
  levels' own contributions: a genuine interchange should improve **both** sides, whereas a
  large gain on one level alone says that level's label is wrong while the other's true
  identity may be a third level — possibly one not yet observed, which this test cannot
  reach.

**Which row of `enlev.dat` a level is.** The run's level list and `enlev.dat` share no
identifier, so a level is tied to its row by **energy**: the nearest row within
`ENLEV_MATCH_TOL` = 0.5 cm⁻¹, which is far coarser than the agreement the two files
actually have. Only the rows IDEN2 marks **found** — the asterisk in columns 39–40 — are
searched, because in an unfound row the observed-energy field holds a copy of the
calculated energy, a position nobody has measured, and matching against it would be
meaningless. The asterisk is the *only* thing consulted. In particular the uncertainty
column is **not** a second opinion on it: the value there starts life as a placeholder
written for every level by the conversion from the Cowan-code output, and IDEN2 changes it
only when the user orders that change for that level, so a level that has been found may
still be carrying the placeholder. Reading the uncertainty instead lost `059003.000623` in
the run of 2026-09-05 — the level was in the file, row 474, starred, at the right energy,
and was reported as having no calculated counterpart. A level of the run that matches
nothing is now printed together with the nearest row of `enlev.dat` whether that row is
starred or not (`nearest_enlev_row()`), because a row sitting at the right energy but
unstarred means the two files are out of step, which is a likelier cause than a level the
calculation does not contain.

**The search window.** A pair is a candidate when its separation is at most `--window`
(default 1.0) times the larger of the two levels' configuration rms above — the error of
the calculated energies *there*, not a global number. The report also gives how far the
swap would move each level from its **new** calculated position, in units of that same rms
(`omc_swap_A_sigma`, `omc_swap_B_sigma`): a swap that leaves both within about one rms is
energetically as good as the assignment; one that pushes a level to three rms is not,
whatever the intensities say.

**How large an `lnR` is large.** It is calibrated the way the rest of this validation
calibrates everything — against swaps that are false by construction. Every same-parity,
same-J pair separated by *more* than the window but less than `--null-window` (default 10)
times the rms gets the same `lnR`; those levels are still in the same part of the spectrum,
so their *gA* values are of comparable magnitude, but the energies rule the swap out. Each
candidate's `p_null` is the fraction of that sample lying at or below its own `lnR`.

#### The run of 2026-09-05

594 levels, **all 594** matched to a labelled level of `enlev.dat` by energy. **18 candidate
pairs**, against **395 calibration swaps** whose `lnR` has median +71, 5th percentile +7.7,
1st percentile −5.8 — only **8 of 395 (2.0%)** favour the swap at all, which is this test's
rate of accidental preference. Applying the whole flagging rule to those false swaps flags
**1.27%** of them, so among 18 candidates **0.23 false flags** are expected.

Two pairs are flagged (`lnR` < 0, both slopes nearer 1 under the swap, `p_null` ≤ 0.05):

| A | B | J, par | sep (cm⁻¹) | window | slopes A | slopes B | rms | `lnR` | `p_null` |
|---|---|---|---:|---:|---|---|---|---:|---:|
| `059003.000483` `f25f ~3F2F` | `059003.000398` `f25f ~3F4D` | 5/2 o | 93.5 | 104.1 | +0.46 → +0.92 | +0.54 → +0.94 | 2.10 → 0.91 | −30.9 | 0.003 |
| `059003.000447` `f26p ~3P2D` | `059003.000298` `f5d6s ~3F4F` | 3/2 o | 139.4 | 184.1 | +0.48 → +0.61 | +0.30 → +0.34 | 1.63 → 1.38 | −4.1 | 0.015 |

The first is the stronger case by an order of magnitude, and it corroborates independently:
`059003.000483` is the tested level with the **largest** energy shift of the whole run
(ΔE = +0.77 cm⁻¹, `d` = +10.7, top of the `level_shifts.py` worst-shift list), and the swap
would leave both levels within 1.1 and 1.4 rms of their calculated positions. Its
line-by-line table (`--detail`) shows the reason plainly. Under the swap the 27 lines of
the two levels fall onto their predicted intensities: the rms of `ln(I_obs/I_calc)` drops
from 2.10 to **0.91**, which is the scatter of the intensity scale itself (sd 1.25 over the
whole run), and both slopes go from about ½ to about 1. The individual outliers go with it
— for `059003.000398` the two worst lines read `ln(I_obs/I_calc)` = +4.11 and +4.03 as
assigned and +0.46 and +0.17 under the swap; for `059003.000483`, +6.67 and +4.02 as
assigned and −0.56 and −2.13 under the swap.

`059003.000572` — the `4f.5d.6p` level with the worst `pattern_V` of the run, the case that
motivated this tool — appears among the 18 candidates twice (paired with `059003.000419`
and with `059003.000573`) and is **not** flagged: `lnR` = +27.4 and +41.5, and either swap
would push it to 4.9 or 4.4 rms of its calculated position. Its bad pattern is therefore
not an interchange with an observed neighbour. That leaves the two possibilities this test
cannot separate: its true partner is a calculated level not yet found experimentally, or
the `4f.5d.6p` eigenvectors are simply too poor there to predict branching at all.

A flag asks for a look in IDEN2; it does not relabel a level. Relabelling changes no energy
and no identification — only which theoretical level the published table names — so it is a
decision for the analyst, taken with the term structure and *g*-factors in hand. When the
decision is made, `swap_line_assignments.py` (next section) carries it out.

#### Usage

```
python level_interchange.py                       # the report, over the pipeline run
python level_interchange.py --lopt LOPT_output_lines_revised.txt --e-input revised_level_energies.csv
python level_interchange.py --detail 059003.000483 059003.000398
python level_interchange.py --window 1.5 --min-lines 2
```

`--detail` prints the two levels' accepted lines with the intensity each hypothesis
predicts for each — the table to read before accepting or dismissing a flag. Its `n` column
is how many accepted transitions share the measured line, and `share as` / `share sw` are
the fraction of the measured intensity the line contributes under each hypothesis, so a
blended line can be read for what it is. Output:
`level_interchange.csv` (+ `.xlsx` twin), one row per candidate pair.

### Making the exchange: `swap_line_assignments.py` and its three steps

A flag from `level_interchange.py` says that two levels' **theoretical identities** were
interchanged: both measured energies are right, but each was attached to the wrong
calculated level — the wrong configuration, the wrong term, and with them the wrong set of
calculated transition probabilities.

The obvious repair is to exchange the two **level identifiers** (`059003.000483` and the
like). That is the one thing that cannot be done cheaply. An identifier is only a name, but
it appears in the classification tables, in the LOPT input, in the Cowan-code bookkeeping
kept in a different folder from the experimental analysis, in lookup tables that translate
between the two, and in copies of some of those files in this repository. Nothing keeps
them in step. And IDEN2 cannot follow at all: it knows a level by its **row number in
`enlev.dat`**, that file is sorted by *calculated* energy, so a row number is a statement
about the calculation and cannot move without re-sorting the file and invalidating every
reference to a row by number.

The cheap repair is the mirror image. **Keep both identifiers exactly where they are and
exchange everything the measurement gave them:** their observed energies, and every
observed line assigned to them. The result is the same table — level A is now the
theoretical level B used to be, with B's calculated intensities and B's old energy — and
nothing outside these files has to be told about it.

Three scripts do it, one per place the exchange has to be made, and none of the three can
see the other two. `swap_line_assignments.py` runs them in order and is the ordinary way in.
All four take the two identifiers as their only mandatory arguments, report before they
write, keep the previous contents of every file in `<file>.bak`, and have `--dry-run`.

```
python swap_line_assignments.py 059003.000483 059003.000398 --dry-run
python swap_line_assignments.py 059003.000483 059003.000398
```

| script | what it rewrites |
|---|---|
| `swap_line_assignments_LOPT.py` | `LOPT_input_lines.txt`, `LOPT_fixlev.txt` |
| `swap_line_assignments_IDEN.py` | `IDEN2/enlev.dat`, `IDEN2/trans.dat` |
| `swap_line_assignments_pipeline.py` | `revised_level_energies.csv`, `line_decisions.csv` |

Nothing in `Icalc.xlsx`, in `icalc_new.xlsx` or in the published line workbook changes. The
calculated transitions are keyed by the **pair of level identifiers**, and an identifier
keeps its configuration, its term and its calculated transition probabilities through the
exchange — what moves is the measurement. The calculated intensities are *read*, by the
first script, to divide the blends again (below). `classify_lines.py` reads only the `Icalc` and
`u%gA` columns of that file, so a level's calculated intensities do not depend on where the
level sits at all; the predicted wavenumbers are computed from the level energies, which is
exactly what the third script changes.

**Where the files are looked for.** Every default file name is looked for first in the
directory the command was run from and then in the directory holding the scripts
(`swap_paths.working_path`). Working in an iteration folder therefore acts on that folder's
copies, while the files that exist in only one place — the decision ledger, the level
overrides — are still found. Each script names the file it settled on, and says which of the
two directories it came from, before it writes anything.

#### `swap_line_assignments_LOPT.py`

LOPT does not read level energies; it computes them from the lines assigned to each
identifier. Exchanging the lines therefore exchanges the energies, and nothing else needs
saying. In `LOPT_input_lines.txt` every record naming `id1` at either end is rewritten to
name `id2`, and the other way round. The columns of the two identifier fields are read
from `LOPT.par` (`--par`, or `--columns` to override), so a re-arranged layout is followed
automatically; every other character of every record — wavenumber, flags, spacing, DOS line
endings — is preserved byte for byte, with the single exception of the blend weights
described next, and a difference of the before and after files shows those two columns, the
weights of the disturbed blends, and nothing else.

**The shares of a blend are recomputed.** Where two transitions fall on one measured line,
`make_LOPT_input.py` writes one record per component and divides the line's weight, 1,
between them in proportion to their **calculated** intensities, so that LOPT gives each
level the share of the measured position the calculation says belongs to it. An exchange
gives the moved component a different pair of levels, for which the calculation predicts a
different intensity, so the share it was given no longer applies. The script therefore reads
`Icalc.xlsx` and `icalc_new.xlsx` — through `classify_lines.read_transitions`, so that the
numbers are the ones the pipeline itself would use — and divides each disturbed line again
by the same rule, reporting every weight it changes.

It recomputes a blend only when two things hold: the calculated file gives an intensity for
every component of it, before and after the exchange, and recomputing the shares the file
*already* carries reproduces them to within `--reweight-tol` (0.002). The second test is
what makes the first safe — it checks, on the file itself, that the script is reading the
same intensities and applying the same rule as the run that wrote the weights. A transitions
file left over from an older classification, an intensity model refitted since, or a
component whose intensity had to be imputed rather than read therefore leaves that blend
exactly as it is, with the reason printed. `--no-reweight` turns the whole of it off.

A share too small for the four-decimal weight field is written in exponential form —
`2.5e-5` rather than `0.0000` — which LOPT reads just as well. Nothing is written as an
exact zero: LOPT stops with an error when a record that is not flagged `P` (predicted) or
`M` (masked) carries the weight 0, so a blend in which the calculation gives a component
no intensity at all is left alone instead. Any component that comes out below a tenth of
its measured line is listed as a **warning**: a share that small says the calculation does
not really place the line on that transition, and rejecting the assignment is worth more
than fitting it with a negligible weight. For `059003.000237` ↔ `059003.000234` three of
the sixteen recomputed components fall below that mark — 0.0400 at 45293.419, 0.0034 at
43640.810 and 0.0717 at 38164.609 cm⁻¹.

This matters for the LOPT refit that `swap_line_assignments.py` runs between its first and
third steps: those are the energies the third step writes into `revised_level_energies.csv`.
The next full `classify_lines.py` → `make_LOPT_input.py` rebuild works the same shares out
again from scratch.

Three things stop it before anything is written: an identifier that does not occur in the
file at all (a typo would otherwise produce a silently unchanged file); a record joining
the two levels to each other, which the exchange would turn into a transition of a level
with itself; and a duplicate created by the exchange, i.e. one observed line ending up
assigned twice to the same pair of levels.

`LOPT_fixlev.txt` holds the levels whose energies are held fixed. A fixed energy belongs to
a position in the spectrum, not to a name, so if both levels are fixed their values are
exchanged with them. If only one of the two is fixed the script stops: that is a decision
about the analysis, not something to guess. (Neither level of either flagged pair is
fixed — only the ground level is.)

`line_decisions.csv`, the ledger of hand-made orders, is **reported on but never
rewritten**. Every order naming one of the two levels is printed together with the
identifier it would have to become. They are left alone deliberately: each order carries a
reason written in prose, and after an exchange that reason usually has to be re-worded, not
merely re-keyed.

One more thing the script says at the end and that is worth repeating here:
`LOPT_input_lines.txt` is *generated* from `line_classifications.csv` by
`make_LOPT_input.py`, so re-running the classification pipeline rebuilds it and undoes the
exchange. The exchange is a statement about which theoretical level each measured energy
belongs to, and to survive a rebuild it has to be recorded where the pipeline reads it —
in the calculated-intensity bookkeeping.

#### `swap_line_assignments_IDEN.py`

Here the energies are written down, so they have to be moved by hand along with the lines.

`enlev.dat` — one row per calculated level, fixed-column:

| columns | contents |
|---|---|
| 1–4 | the row number, which is what IDEN2 calls the level |
| 5–16 | the calculated energy |
| 17–26 | the uncertainty of the observed energy |
| 27–38 | the observed energy |
| 39–40 | `" *"` when the level has been observed at all — the **only** marker of that; the uncertainty column is a placeholder until someone changes it |
| 41–52 | observed minus calculated |
| 53–57 | J |
| 58– | the label: configuration, separator, term |

The uncertainty, the observed energy and the asterisk belong to the measurement and change
places. The calculated energy and the label belong to the calculation and stay. Observed
minus calculated is then recomputed for each row from **its own** calculated energy, and
what comes out is worth looking at — it is the O−C the exchange implies.

`trans.dat` — the predicted transitions, written in blocks. A header row beginning with
`$` names a level; the rows under it are that level's predicted transitions to levels lower
down the file, one row each:

| columns | contents |
|---|---|
| 1 | always `+` |
| 2–5 | the row number of the level at the other end |
| 6–10 | the predicted intensity of the transition |
| 11–22 | the observed energy of that other level |
| 23–24 | `" *"` when it has been observed |
| 25–37 | the predicted wavenumber, \|E(this) − E(other)\| |
| 38–43 | the intensity of the observed line assigned to it, `0` if none |
| 44–55 | the wavenumber of that observed line |
| 56–66 | observed minus predicted |
| 67–72 | the row of that line in `dlv.dat`, `0` if none |

Each predicted transition appears exactly once, in the block of whichever of its two levels
stands higher in the file. The last four fields are the **assignment**, and they are what
moves: for every level `p`, the assignment of the transition `(n1, p)` and that of `(n2, p)`
are exchanged. Then every row that mentions either level has its observed-energy field, its
asterisk and its predicted wavenumber recomputed, and its observed-minus-predicted with
them. `dlv.dat` and `numset.dat` carry no level information and are not touched.

**A self-check that comes for free.** A line assigned to `(n1, p)` sat at some
observed-minus-predicted value. After the exchange it is assigned to `(n2, p)`, and `n2`
now holds the energy `n1` used to hold, so its new row predicts the same wavenumber as its
old one. **Every transferred line must therefore keep its O−C exactly**, and the script
verifies this line by line and stops if any of them changes: that would mean the files were
not in the state the rewrite assumes.

**Lines that cannot be transferred.** A predicted transition is in `trans.dat` only if the
calculation gave it an intensity above the printing threshold, and the two levels being
exchanged do not have the same set of predicted transitions. A level `p` may therefore be
connected to `n1` in the file and not to `n2`; an observed line assigned to such a
transition has nowhere to go. Those lines are listed one by one — the transition they were
on, their wavenumber, their row in `dlv.dat` — and their assignment is removed, because
leaving it would file the line under the wrong level, which is the mistake being repaired.
They are not lost: they are still in `dlv.dat` and can be re-identified by hand. Read this
part of the report. A level that loses several strong lines this way is a sign that the
exchange needs a second look. (Neither of the two flagged pairs loses any.)

**Naming the two levels.** Each of the two arguments is **either** an experimental
identifier **or** the level's row number in `enlev.dat`, which is what IDEN2 itself calls
the level. A bare run of digits is read as a row number and anything else as an identifier,
so the two forms need no flag to tell them apart and may be mixed:

```
python swap_line_assignments_IDEN.py 721 723                        # IDEN2's own numbers
python swap_line_assignments_IDEN.py 059003.000483 059003.000398    # identifiers
```

A row number needs no translation, and is the shortest way to name a level with IDEN2 open
in front of you. Identifiers are accepted because this script is the companion of the other
two, which can be given nothing else — the LOPT files and the pipeline know levels only by
identifier — and because `swap_line_assignments.py` drives all three from one pair of
arguments. An identifier has to be translated, and there are three ways, in order of
precedence: `--index N1 N2` gives both row numbers directly; `--map FILE` reads one of the
lookup tables kept alongside the analysis (a column of identifiers and a column of numbers,
header optional); and by default the identifier is looked up in a table of experimental
energies (`--levels`, LOPT's output level list) and matched to the row of `enlev.dat` whose
observed energy agrees to within `--tol` (0.5 cm⁻¹), the match having to be unique.

The energy route reads the state **before** the exchange, so run the script once, on the
files as IDEN2 last left them; naming the rows by number avoids the question entirely. Two
guards catch the usual mistakes: it stops if the two rows do not have the same J, and — when
the rows were named directly rather than found by energy — if they already hold each other's
energies, which is what a second run looks like.

`--iden2-dir` points at the directory (so `IDEN2_snr`, the copy carrying signal-to-noise
ratios in place of intensities, can be given the same exchange), and `--out-dir` writes the
rewritten files somewhere else instead of over the originals.

#### `swap_line_assignments_pipeline.py`

The other two scripts rewrite **working files**. Neither is an input of the classification.
`make_LOPT_input.py` writes the LOPT transitions file out of `line_classifications.csv`,
which `classify_lines.py` produces from the published level list, the observed line list and
the calculated transitions — so re-running the pipeline rebuilds the transitions file from
those and the exchange is gone. This script is what makes it stay, and it writes three
things into the two small files that are under the analyst's hand.

**The two energies**, into `revised_level_energies.csv`. `classify_lines.py` gives every
level the energy of the published list, overridden by this file, and generates that level's
candidate transitions there. The exchange moved both energies, so both belong here: one row
each, the identifier keeping its own calculated identity and taking the other's measured
position. A level already listed keeps its row and has its energy replaced — a second row
for the same level would be read as a contradiction.

**The orders already in the ledger.** `line_decisions.csv` keys each verdict by an observed
wavenumber and the two level identifiers, written out, so every order naming one of the
exchanged levels now points at the wrong one. Each is re-keyed, and the reason it carries is
kept and annotated — `; Lines belonging to levels A and B were swapped on <date>.` — rather
than rewritten, because the reason is prose about a measured line and the line has not
changed, only which identifier it is filed under.

**The identifications of the published line list.** This is the part that is easy to miss.
An identification that is Sugar's own — the rows `line_classifications.csv` marks `new` = 0
— names its two levels explicitly in the line workbook, an input file that is never edited.
Step 3 keeps such an identification unless something rejects it ("old, no solid evidence for
rejection"), and it goes on naming its level after the exchange has moved that level a
hundred wavenumbers away; the level is then fitted between its true lines and its stale
legacy ones. That is what pulled `059003.000424` 81 cm⁻¹ off its LOPT position on
2026-09-05. So for each of them the script writes two orders: a `reject` where the line
stands and an `accept` under the other identifier, which is the level that now sits where
the line is. `--rejects-only` writes only the first and leaves the rest to the
classification. `classify_lines.py` now reads the exchange out of the comment column of
`revised_level_energies.csv` and carries the old status over by itself (see
[The comment column](#the-comment-column-and-what-becomes-of-the-published-identifications)),
so these orders are no longer what stands between the exchange and a level fitted to its
stale legacy lines; they remain a written record of the two verdicts and are still worth
having, and where they exist already the two agree. An order that would contradict one already in the ledger stops the script,
because `classify_lines.py` aborts on two rows ruling differently on the same assignment —
better here than on the next run.

**Where the new energies come from.** From `LOPT_output_levels.txt`, and a LOPT run made
either before or after the exchange will do. An exchange does not invent energies: both were
measured already, and what changes is which identifier each belongs to, so the two numbers to
be written are the two the level table holds, the other way round. A LOPT run made after the
exchange has already written them that way — each identifier refitted with the other's lines,
which differs from the plain exchange of the two old numbers by a few thousandths of a
wavenumber — and then they are taken as they stand.

Which of the two cases it is in, the script works out, by comparing the table with the
energies the pipeline currently has (the `low_E` / `upp_E` columns of
`line_classifications.csv`):

| the level table | means | what is written |
|---|---|---|
| still agrees with the pipeline | LOPT has not been re-run | those two energies, exchanged here |
| has them **crossed** | LOPT has been re-run | those two energies as they stand |
| agrees with neither | something else moved these levels | nothing; the script stops |

`--exchange` and `--as-fitted` force one of the two readings, which is needed only when the
pipeline has no energy for either level and the comparison cannot be made.

Running it twice does not undo the first run. Every row it writes or annotates says in plain
words that these two levels have been exchanged, and a row that already says so is left
alone. So running it once before the LOPT run and again after it is the ordinary way to
work: the second run replaces the two override energies with the fitted ones and changes
nothing else.

#### `swap_line_assignments.py`

Runs the three in the order above and stops at the first that refuses, so that a pair the
LOPT step rejects — a typo in an identifier, a line joining the two levels, a fixed energy
on one of them only — never reaches the other two. `--dry-run` runs all three as reports.

Between the first step and the third it also tries to run LOPT itself: the command `lopt`
(`--lopt-command`, `--no-lopt-run` to skip it), in the directory the LOPT files it has just
rewritten were found in — LOPT reads `LOPT.par` from the directory it is started in, and that
file names the transitions file. The refit is not needed to record the exchange, but it turns
the two exchanged numbers into least-squares energies fitted with every other level, and it
is the check that the rewritten transitions file is one LOPT will still read.

If the command is not there or returns an error, that is a **warning** and the run continues:
the third step then sees a level table that does not yet know about the exchange, says so,
and exchanges the two energies itself. Running LOPT by hand afterwards and then
`swap_line_assignments.py id1 id2 --after-lopt` replaces those two energies with the fitted
ones and leaves everything else alone — worth running after any later LOPT run that has moved
the two levels.

`--only lopt,iden` and the like select steps, and `--index N1 N2` passes IDEN2's row numbers
straight to the second step.

```
python swap_line_assignments.py 059003.000483 059003.000398 --dry-run
python swap_line_assignments.py 059003.000483 059003.000398   # runs lopt itself
python swap_line_assignments.py 059003.000483 059003.000398 --after-lopt   # if it did not
python classify_lines.py
```

Afterwards check in `lopt_vs_classify_levels.csv` that each of the two levels came out at
the energy that was written for it. A level fitted between two groups of lines has kept an
identification at its old position.

#### The two flagged pairs

Both exchanges were checked with `--dry-run` against the files of 2026-09-05. In each the
calculated order and the observed order of the two levels are **inverted** as assigned, and
the exchange puts them back in the same order — an argument independent of the intensities,
and the reason the inversion had to be retained in the least-squares fit as long as the
published compositions were trusted.

| | `059003.000483` / `059003.000398` | `059003.000447` / `059003.000298` |
|---|---|---|
| IDEN2 rows | 721, 723 | 943, 945 |
| labels | `f25f ~3F2F`, `f25f ~3F4D` | `f26p ~3P2D`, `f5d6s~3F4F` |
| E(calc) | 117545.0, 117478.9 | 84331.8, 84184.2 |
| E(obs) as assigned | 117592.99, 117686.46 | 84270.53, 84409.90 |
| O−C as assigned | +47.99, +207.56 | −61.27, +225.70 |
| O−C after the exchange | +141.46, +114.09 | +78.10, +86.33 |
| records rewritten in `LOPT_input_lines.txt` | 45 | 29 |
| observed lines moved in `trans.dat` | 27 (13 ↔ 14) | 17 (10 ↔ 7) |
| lines that could not be transferred | none | none |
| orders in `line_decisions.csv` to re-word | 2 | 1 |

The exchange does not make either level's O−C small — these are `4f².5f`, `4f².6p` and
`4f.5d.6s` levels, whose calculated energies are good to 104, 114 and 184 cm⁻¹ rms — but it
takes the pair from one level near its calculated position and one far from it to two levels
about equally close, which is what a correct pair of identities looks like. In the first
pair the two O−C values become 1.4 and 1.1 times the rms of `4f².5f`; in the second, 0.7
of the `4f².6p` rms and 0.5 of the `4f.5d.6s` rms.


### Is the level where it is put: `level_positions.py`

**What it answers.** Every other test in this chapter asks whether a level is *real*. This one
asks the sharper question the question marks need: granted that something is there, is the
adopted energy the only energy the observed lines allow? The answer is a likelihood ratio

    ln R(L, E) = ln [ P(what was recorded | a real level sits at E)
                    / P(what was recorded | nothing sits at E) ]

evaluated at the adopted energy and then scanned across the interval Cowan's calculation allows.
Every local maximum within a few units of the adopted one is an alternate position, and a level
with any has earned its question mark. **The count of them is the verdict**, reported as a column
of its own; no single p-value stands in for it.

**One sum, not four factors.** The four things a level answers for — the improbability of its
alignments, the intensity pattern, the blends, the absences — are usually four penalties
multiplied together, and the arbitrary part is then how to weigh them against one another. Here
they are not separate terms. ln R is summed over the level's **predicted** transitions rather
than over its accepted lines, and that one choice does the work: a prediction with no line near
it is simply a prediction whose datum is "no line", so the alignment evidence and the
missing-line penalty become the two branches of a single formula. The intensity test is a second
factor inside the branch that says the line is genuine. The blend penalty is what happens to the
same formula when the feature is already explained by transitions that have nothing to do with L.

**The per-prediction formula.** For a predicted transition t of L whose partner M sits at E_M,
the Ritz wavenumber is nu = E − E_M when L is the upper level and E_M − E when it is the lower —
the two move in *opposite* directions as E is scanned, which is what makes the maximum sharp
rather than flat. Model the line list near nu as a Poisson process of unrelated lines of local
density rho — *unrelated*, so the density is that of the lines the level list does not already
account for, not of the whole recorded list; see **The background of unrelated lines** below —
plus, under H1, at most one genuine line, present with probability P = c(lambda)·D(z)
and placed Gaussian about nu. The Poisson background is identical under both hypotheses and
cancels, leaving, with p = P(1 − eta):

| the datum at nu | R_t |
|---|---|
| a line at residual d | (1 − p) + p · N(d; 0, sigma) / rho · G |
| no line in the window | 1 − p (1 − q_out),  q_out = 2Φ(−4) = 6.3e-5 |

c is the coverage function of `tools/coverage_map.py` and D the detection curve of
`tools/obscuration_rate.py`, reached through `level_shifts.observation_probability`. The absence
of a prediction is therefore charged in exact proportion to how surely it would have been seen,
and charged nothing where the plate was blind or the line sat at the noise.

The first branch of the matched row is the world in which the line is there by chance even though
the level is real. That is why the intensity ratio G multiplies only the second branch, and why
**no single row can cost more than ln(1 − p)** however badly placed or however wrong in intensity.
The bound is not a nicety: the residuals of this run are measurably heavy-tailed — |t| > 4 occurs
twenty times as often as a Gaussian allows — and under a pure Gaussian one corrupted measurement
would charge e^-8 and kill a level that is perfectly real. Here a bad line merely fails to help.

**The intensity ratio, and where the blend penalty comes from.** A feature is *free* when nothing
else is accepted on it and *claimed* when other accepted transitions, between levels that do not
involve L, already account for it. With I_t the predicted intensity of t, C the summed predicted
intensity of the claimants, f(lambda) the far-ultraviolet scale correction, s the measured scatter
of ln(I_obs/I_pred) and g the local distribution of ln I_obs over the whole observed list — what
an unrelated line would have been drawn from —

| feature | G |
|---|---|
| free | p(ln I_obs given I_t f) / g(ln I_obs given lambda) |
| claimed | p(ln I_obs given (C + I_t) f) / p(ln I_obs given C f) |

Both numerator and denominator are densities of the same variable, so no Jacobian survives.

The density p is a Gaussian of width s about the prediction, **floored at g on one side only** —
where the feature is brighter than predicted, not where it is dimmer. The asymmetry is the
physics: a feature can always carry light from transitions nobody has identified, so an excess
over the prediction is no evidence against anything, while a deficit cannot be explained away,
because light already emitted cannot be taken back. The floor is not cosmetic. Without it the
denominator of a blended feature is evaluated at the summed prediction of its other components,
and where those happen to be weak that density is astronomically small — so a level collects
enormous credit for "explaining" a brightness the components it is compared with never claimed to
explain. One such row, an observed feature of intensity 7874 whose other claimant is predicted at
2.3, was worth ln G = +22, and it alone moved a solidly established level 219 cm^-1 onto the
position of a different real level.

The blend penalty is then not an added term but the same formula with two densities changed.
Under H0 a claimed feature is present with *certainty*, so its presence is no evidence for L: in
the positional part rho is replaced by 1/(2W), the density of a line known only to be somewhere
in the window, and the best attainable R_t falls from about 25 to 3.2. In the intensity part the
question stops being "is a line of this brightness here" and becomes "does adding the predicted
intensity of t improve an account of that brightness which already works".

**Two readings of a match, and the row is worth the better of them.** A matched free feature can
be read two ways. Either *the feature is the transition* — its position drawn from N(d; 0, sigma)
against the local line density rho, its brightness from N(ln I_obs; ln(I_t·f), s) against g — or
*the transition is hidden in a feature that is there anyway*, carrying light nobody has identified.
Under the second reading H0 has the line present with certainty, so rho is 1/(2W) as for a blend
and the brightness says nothing either way; the best it can be worth is 3.2 against the first
reading's 25. R_t takes whichever is larger.

Taking the intensity floor without also replacing rho was the error: it read the brightness under
the second hypothesis and the position under the first. A prediction of I = 0.2 landing on a
feature of 13372 then kept the whole of the positional credit — ln R_t = +1.4 for a transition with
a 1.5 % chance of having been recorded at all. Rows of that kind are the strong observed line
assigned to the very weak transition, the assignment `classify_lines.py` refuses and the analyst
leaves free for a more adequate one, and half a dozen of them were carrying four levels of the
2026-09-10 audit (`059003.000521`, `.000522`, `.000534`, `.000536`) to alternate positions. Taking
the better of the two readings, rather than switching on the floor, is what leaves a genuine line
alone: one under-predicted by a factor of thirty is still far better explained as the transition
than as a coincidence, and it goes on being read that way. A row the second reading wins counts as
neither free support (`n_free`) nor evidence in the level offset, and `--detail` marks it `bright`
in the `what` column. In the 2026-09-10 run the correction leaves the median `ln R` at 22.4, takes
the preferred alternates from 57 to 30 and the `refit` rows from 5 to 1, and — because it lowers
`ln R` at an adopted position as readily as at an alternate — promotes one level, `059003.000477`,
from `weak` to firm grounds for relocation: the four free lines 100.5 cm⁻¹ away were always worth
`ln R` = +2.5, and what changed is that the position it sits on now scores −5.4 instead of −3.3.

**The background of unrelated lines.** rho is the rate at which a line the level under test has
nothing to do with turns up near nu, and g is the brightness such a line would have. Both were
first measured on the whole recorded list, from the 41 and 201 nearest lines. That is the wrong
population, in both halves. Of the 6668 recorded lines, 4526 already carry an accepted transition:
they are explained without the level under test, and they are systematically the *bright* ones,
because a line was identified in the first place partly by being strong enough to measure well.
Only the remaining 2142 can be what a coincidence is drawn from. Using all 6668 puts rho about two
and a half times too high and the mean of g about 0.6 in ln I too high — a bias that understates
every matched free row, by ln(0.094/0.039) ≈ 0.9 in the positional half alone.

The free lines are not a mystery, and this is where `tp_E1_no_trials.xlsx` — Cowan's complete E1
transition list, which was not in the project when the likelihood was written — earns its place in
the null model. A recorded line that carries no accepted transition is, in this spectrum, almost
always a transition of a level the calculation predicts and nobody has found: **658** of the 1253
calculated levels of `IDEN2/enlev.dat` are in that state, **61290** of their E1 transitions have a
partner that *has* been found, and summing the probability P that each would have been recorded
gives **≈ 2970** expected lines — the same order as the 2142 free ones actually there. The file
therefore predicts the very population rho and g describe, and predicts its *structure* as well as
its size: the expected rate runs from 0.001 per cm⁻¹ below 10000 to 0.052 near 55000, more than an
order of magnitude, which no single number and no 21-nearest-neighbour smoothing of sparse data can
express.

An unfound level's energy is known only to W, the rms of E_obs − E_calc over the found levels of
its configuration (40 to 800 cm⁻¹), so one of its transitions predicts no position. Smeared over W
it predicts a *rate*, which is exactly what a Poisson background needs:

    rho_unk(nu) = Σ over unfound-level transitions of  P · N(nu; nu_calc, W_cfg)   per cm^-1

and the same weights give the mean and spread of ln(I·f) — the brightness such a line would have.
`unknown_transition_background()` builds both as one convolution per configuration width, about two
seconds on top of `build()`.

Two estimates of one population, then — the free lines actually recorded, and the lines the
calculation says are missing — and they are combined by taking the **larger** of them, the larger
rate and the larger of the two brightness densities. That is the conservative direction, since a
bigger background is a smaller `ln R`, and it is the same device as the two readings of a match: an
unrelated line is allowed whichever account of itself is the better. The result is capped above by
the all-lines density (the free lines are a subset of the recorded ones) and floored at
`RHO_BG_FLOOR` = 0.05 of it, because both estimates extrapolate where no free line is near and an
extrapolated density must not buy unbounded credit; the brightness half is floored the same way at
`LN_BG_DROP` = 3.0 below the all-lines value.

Three places this deliberately does **not** reach. A *claimed* feature keeps the all-lines rate and
the all-lines g: there rho has already been replaced by 1/(2W), and the one-sided floor of
`ln_intensity` is a guard against an astronomically small denominator, for which the more generous
distribution is the safer one. `eta` keeps the all-lines rate too: it is the rate at which a
*genuine* line sits at an anomalous position, and the background it is mixed against is every line
it could have been taken for. And rho_unk does not enter the *hidden in a feature that is there
anyway* reading — that reading compares a world in which the observed feature is an unrelated line
with one in which the same unrelated line has the transition blended into it, so the feature is
present in both, its rate cancels, and only 2W·N(d; 0, sigma), the probability that the transition
falls inside it, survives.

**What it does to the run** (`--audit --no-firm`, both ways, 593 levels). Because the correction
removes a background that was too large, every level's `ln R` **rises** — the median goes from 22.4
to 31.6, a little under one unit per matched free line — so `FIRM_LN_R` and any threshold read off
an absolute `ln R` mean slightly less than they did before. Differences between two positions of
the *same* level, which is what the scan and the audit are about, are far less affected: the shift
is common to both wherever they match the same number of free lines.

| | plain background | unrelated-line background |
|---|---|---|
| median `ln R` at the adopted position | 22.4 | 31.6 |
| levels at `ln R` ≤ 0 | 21 | 7 |
| local maxima within `--alt-drop`, summed over all levels | 2823 | **1912** |
| levels carrying a question mark | 85 | **68** |
| of those, alternate preferred to the adopted position | 30 | 27 |
| `relocate` / `interchange` / `refit` / `top line` / `weak` | 1 / 3 / 0 / 1 / 17 | 2 / 1 / 2 / 0 / 16 |

The line that matters is the third. A third of the alternate positions were **manufactured by the
background itself**: scoring a coincidence against a rate two and a half times too low made it look
like a match, and enough such rows within a scan window raise a local maximum. With the right rate
the surface is sharper and 911 of those maxima are gone, and with them 17 question marks. That also
loosens the look-elsewhere correction, which is ln(n_alt): `059003.000228` moves from `weak` to
firm grounds for relocation not because its alternate gained (it did, +5.9) but because the count of
rival maxima in its window fell from 107 to 10, and ln(n_alt) is a proxy sensitive to exactly that.
Read the new `look` values with that in mind. `--plain-background` restores the old behaviour
exactly, for comparison.

The C form is the branching-fraction form: ln(I_obs/((C + I_t)f)) is identically
ln(I_obs·BF/(I_t·f)), since BF is by definition I_t/(C + I_t). One transition's prediction is
never compared with a whole blended feature. C is used because it is defined at any scanned
position, where BF is not.

**The level offset.** All of a level's lines share whatever the calculation got wrong about that
level, so the intensity residuals are written r_i = mu_L + e_i with mu_L ~ N(0, s_L). mu_L is
marginalised, not fitted, and the reason is a measurement: the within-level spread is s = 1.215
against a between-level spread of level means of s_L = 0.280, so the offset explains about a
twentieth of the variance, while the intensity model is absolutely calibrated (mean residual
−0.006). Fitting it freely would throw the calibration away to buy very little. Marginalising adds
one scalar correction, −½ln(1 + n s_L²/s²) + s_L²(Σr)²/(2s²(s² + n s_L²)), **with each row
weighted by its posterior of being genuine rather than a coincidence**. The weighting is what
keeps the correction bounded: its positive part is meant to be cancelled by the −r²/2s² the same
rows carry in the base term, but a row whose line is a coincidence has its ln R_t floored at
ln(1 − p) and carries no such term. Unweighted, a displaced position with a few wildly mismatched
intensities collects the gain without ever having paid for it — in testing, +192 where the honest
answer was −8.

**Every ingredient is measured on the run**, and printed at the head of every report so that no
number in it is a guess (values of the run of 2026-09-07, 594 levels, 29260 predicted transitions,
6668 recorded lines):

| ingredient | what it is | value |
|---|---|---|
| sigma_meas | the width of the feature: a wavelength-constant reading error (0.0030 Å for 1974, 0.0040 Å for 1969), a precision floor of 0.0055 cm⁻¹, and an excess for the line's character code, wavelength-constant where the width comes from reading the plate and wavenumber-constant where it belongs to the line itself | 0.006 – 0.01 cm⁻¹ plain below 47500 cm⁻¹; `c` (1974) 0.0337 Å, `cl` (1969) 0.0044 Å, `d` 0.0073 Å, `ch` 0.0103 Å, `h` 0.0082 Å, `bl` 0.0085 Å, and `w` (1974) 0.0251 cm⁻¹ — the one code the residuals put in wavenumber |
| w_hfs | hyperfine width of a LEVEL, constant in wavenumber, fitted on the 1974 lines alone and carried unchanged into 1969; a level with fewer than 4 lines of its own takes the typical width of its CONFIGURATION instead of zero | 274 levels carry one — 61 their own fitted width above 0.020 cm⁻¹, 213 their configuration's; 5 above 0.08 and marked as large; largest 0.102 cm⁻¹ |
| w_cfg | the typical width of a configuration: the rms of the widths applied to the levels of it the fit can measure | 4f².6s 0.0642, 4f².6p 0.0237, 4f.5d² 0.0211, 4f².7s 0.0142, 4f².6d 0.0123, 4f².5g 0.0089, 4f².5f 0.0069, 4f³ 0.0047, 4f².5d 0.0043 cm⁻¹ |
| k(n) | what is left over for a blend of n parts once that width is accounted for | 1.000, 1.483, 1.242, 1.116 for n = 1…4 |
| eta | rate at which a genuine recorded line sits at an anomalous position | 0.0002 (2ΔNLL = 51 against eta = 0) |
| rho | UNRELATED lines per cm^-1: the free lines and what the calculation predicts for the unfound levels | median 0.039 (all recorded lines: 0.004 – 0.336, median 0.094) |
| s, s_L | within-level and between-level spread of ln(I_obs·BF/I_pred·f) | 1.215, 0.280 |
| u_M | partner energy uncertainty, D1 of `LOPT_output_levels.txt` | median 0.014, max 0.410 cm^-1 |

sigma_t² = (k(n)·sigma_meas)² + w_hfs(low)² + w_hfs(upp)² + u_M², and the matching window is four
sigma — set by whichever is worse, the width typical of the neighbourhood or the one the candidate
line itself is modelled at. The quoted `unc_wn_obs` enters only as a floor, and only for a line
quoted more than 1.5 times the usual value of its own character class: within a class the quoted
value is the reading rule times one fixed factor and says nothing the class does not, but a line
quoted far above its class was widened by hand on knowledge no flag records. The codes `**`
(multiply classified) and `*v`/`*r` (hyperfine) carry no width of their own — the first is priced
by k(n), the second by w_hfs. `level_hfs_widths.csv` holds the fitted widths and is rewritten by
`python level_positions.py --fit-hfs`; its `w_applied` column is the width actually used, and
`make_LOPT_input.py` reads the same column, so the uncertainties LOPT is given and the ones the
scan uses cannot drift apart.

**A level the fit cannot measure takes its configuration's width.** A level with fewer than four
lines of 1974 has no width of its own, and used to be given zero — which asserts it has no
hyperfine structure, when all that is known is that nobody has measured it. Its configuration is
the right thing to ask, because that is what decides how strongly the outer electron feels the
nuclear magnetic moment: an s electron has a non-zero probability density at the nucleus and feels
it directly, a 5f or 5g electron never comes near it, and ¹⁴¹Pr is the only isotope, so there is
no isotope shift mixed in. The fitted widths bear that out without being told — pooled by the
outermost electron's orbital letter they come out s 0.0495, p 0.0214, d 0.0127, g 0.0089,
f 0.0055 cm⁻¹, which is the order of penetration — and a permutation test on the configuration
labels of the measured levels gives p < 5×10⁻⁵. A configuration with fewer than `HFS_CFG_MIN`
measured levels of its own is answered for by its outermost orbital (4f.5d.6s takes the 6s width),
then by that orbital's letter (4f².6f takes the f width), then by the list-wide value.

The configuration itself is read from `IDEN2/enlev.dat`, and a level is tied to its row there
through `IDEN2/IDEN_level_ids.txt`, the table giving every level identifier the number of its
row — never by energy. The row number is what survives a level being moved: when a position is
edited in IDEN2 the energy in `enlev.dat` changes and the row number does not, so an energy
match loses exactly the levels that have just been worked on, and hands a level that has no
energy yet whatever row lies nearest. All 595 levels of the run resolve through the table.

Out of sample the rule is worth **+434 units of ln L** over the 2057 accepted unblended 1974
lines, five-fold, against +23 (5–95 %: −42 to +124) when the configuration labels are shuffled
among the levels — no shuffle of forty reached it. 462 of that is earned by 31 lines alone: the
ones Sugar marked `*r` or `*v`, which are given no character width precisely because w_hfs is
meant to carry them, and whose level the fit could not reach, so the model gave them nothing for
hyperfine structure at all. Nothing changes for a level the fit can measure, and the precision
floor is unmoved: re-fitted on the clean plain 1974 lines with the configuration widths in place
it is 0.0054 (0.0049 – 0.0060) cm⁻¹, against 0.0059 (0.0054 – 0.0065) without them, so
`UNC_FLOOR` stays at 0.0055.

Taking the worse of the two widths for the window matters: the neighbourhood alone left 204 of the 4937 accepted
lines outside their own window, because a line quoted to 0.6 cm^-1 among neighbours quoted to 0.1
was being ruled out for lying 0.4 away. Taking the worse of the two leaves 1. Widening costs
nothing, since a line far out has N(d; 0, sigma)/rho well below 1 and contributes the same
ln(1 − p) as an absence.

eta deserves a word, because the measurement contradicted the guess that prompted it. Reading the
residual tails suggested 0.01; maximum likelihood puts it at 0.0002, and the reason is that mixing
in a *flat* background at rho ≈ 0.09 per cm^-1 only helps beyond |t| ≈ 5 — at t = 3.5 the Gaussian
is still ten times rho. The protection against a corrupted measurement does not come from eta at
all. It comes from the (1 − p) branch, which floors every row whatever its residual.

**What is dropped.** A prediction is dropped when its partner M would have no accepted line left
after removing the lines M shares with L: M's energy was then fitted to those very lines, and the
residual is zero by construction rather than by agreement. Seven predictions in this run. (A
hat-matrix treatment generalising this continuously was tried and rejected — LOPT fits a blended
feature through a centroid model, so its normal equations are not those of independent
observations, and leverages computed as though they were make the residuals worse, not flatter:
the rms of t/sqrt(1−h) climbs to 6.0 in the top leverage band. The explicit rule is exact and
needs no matrix. It is also not needed for the scan proper, where E is scanned rather than fitted
and only the partner energies can be degenerate.)

**The scan.** E runs over E_calc ± 3W, with E_calc from `IDEN2/enlev.dat` and W the rms of
E_obs − E_calc over the known levels of the same dominant configuration
(`level_interchange.configuration_windows`) — 133 cm^-1 wide for the best-determined
configurations, 2425 for the worst. The predicted intensities are held fixed while E moves: over
a few hundred cm^-1 the Boltzmann factor with kT = 12900 cm^-1 changes by about two per cent, far
below s. A maximum that lands on another level of the run is reported as such (`near_level`,
`near_dE`), because that is an interchange — `level_interchange.py`'s question, judged on
evidence this scan does not use — and not a free position nobody has claimed.

**What it finds.** Over the 594 levels, ln R at the adopted position runs from −10 to +168 with a
median of 22.6; 21 levels are at or below zero. 122 have at least one alternate position within 5
units, 57 have one the lines actually prefer, and 5 of those fall on another level of the run.
(These figures, and the table below, are from the run before the two-readings correction described
two sections down. That correction moves the median by −0.2, to 22.4, and costs the line-rich low
levels between 5 and 16 units each — it withdraws the positional credit their strongly
under-predicted matches were collecting twice. It bites hardest where it should: the levels with
an alternate within 5 units fall from 122 to 85 and those with a *preferred* alternate from 57 to
30, because a large part of what made an alternate attractive was rows of that kind.)

The test knows nothing about which levels are Wyart's established ones and which are under
review, so the split between them is a check on the method rather than an input to it:

| | levels | median ln R | ln R ≤ 0 | with an alternate | with a preferred alternate |
|---|---:|---:|---:|---:|---:|
| legacy (n_old ≥ n_new) | 384 | 39.0 | 3 | 20 (5 %) | 8 |
| new* (n_new > n_old) | 210 | 8.2 | 18 | 102 (49 %) | 49 |

A ten-to-one difference in the rate of alternate positions between the set that should be safe and
the set that is being tested, from a statistic that was never told which was which.

**Limits.** Three, all in the same direction — they flatter the adopted position:

- The partner energies E_M are held fixed while E is scanned, and they were fitted using L's own
  accepted lines. The degeneracy rule removes the extreme case; the residual effect is that the
  adopted position is favoured over its alternates by a little more than it should be.
- The alternate at a given energy is scored with L's predicted intensities, which belong to L at
  its adopted position. Where an alternate is far away this is only approximately right.
- A maximum found at another real level says the two are interchangeable **on this evidence**;
  it does not say they should be swapped. `level_interchange.py` weighs slopes, common lines and
  branch structure that this scan does not look at.

**How to read `ln R`.** It is a log odds, so the sign is the whole of the meaning. `ln R > 0`: the
recorded lines are more likely if a real level sits at that energy than if nothing does — the
position is supported, and the bigger the number the stronger the support, each +1 being a factor
of e in the odds. `ln R = 0`: the lines say nothing either way. `ln R < 0`: the lines are *less*
likely with a level there than without one — the position predicts transitions that should have
been recorded and were not, or lines whose brightnesses contradict the prediction. **A strongly
negative `ln R` is the signature of a spurious position.** For scale, a well-centred line on an
open plate is worth about +3.2, so `ln R` = +30 is a level carrying ten of them, and the median
level of this run sits near +22. The same quantity is computed at a candidate energy as at the
adopted one, which is what makes the two comparable: `gain` = `ln R`(alternate) − `ln R`(adopted)
is a log odds ratio between two positions of the same level.

**Predicted is not observable.** `n_pred` is how many transitions the calculation gives the level
inside the recorded range, and on its own it says very little: most are far too faint to have been
recorded at all, and their number is set by how many partner levels happen to be known rather than
by the experiment. The count that means something is **`n_obs`** — the predictions whose
probability of having been recorded, P = c(lambda)·D(z), is at least `P_SEEN` (0.2). Of those,
`n_seen` found a line and `n_miss` did not. A level with `n_obs` = 9 and `n_seen` = 8 is in good
order however large `n_pred` is; `n_obs` = 9 with `n_seen` = 2 is not. `n_match` counts every
prediction that found a line, faint coincidences included, and is kept in the csv only — it is the
number *not* to quote.

**`near_level` is a neighbourhood label, not a rival.** It names the level of the run whose energy
is closest to the alternate position, and `near_dE` is the gap between them. It is filled in for
every alternate, so for most rows it names a level that merely happens to lie nearby and wants
nothing: read `dE_alt`, the size of the move, instead. Only when `|near_dE|` is under half a
wavenumber does it carry information, and then it says the alternate *is* that level's position —
an interchange, where the two identities may need exchanging rather than either level moving.

**The audit (`--audit`).** `ln R` says which of two energies the lines prefer. It does not say the
preference is worth acting on, and several quite different situations produce the same number, so
the audit separates them before anything is called a revision. Three corrections stand between
"the lines prefer that energy" and "move the level":

- **Look-elsewhere.** A level whose configuration is 4f.5d.6p is scanned over ±7000 cm⁻¹ and a
  well-determined one over ±130, so the first gets hundreds of chances at a good maximum and the
  second a handful; the raw gains are not comparable. Under H0 the heights of the local maxima are
  close to exponentially distributed, so the best of `n_alt` of them stands about ln(`n_alt`) above
  a typical one, and `look` = `gain` − ln(`n_alt`) puts every level on the same footing. It is not
  a small correction: the median `n_alt` among the levels with a preferred alternate is about 50.
- **Free lines.** `n_free` counts only the alternate's matched rows that are **observable**
  (`P_obs` ≥ 0.2 — the prediction had at least a one-in-five chance of having been recorded),
  well-centred, and on features **no accepted transition claims**; `free_gain` is what they are
  worth. The observability condition is the same one that separates `n_seen` from `n_match`, and
  it matters as much here: a prediction with a 1.5 per cent chance of being recorded that finds a
  line has found a coincidence. Such rows are worth a few hundredths of a unit each, so leaving
  them in inflated `n_free` without moving `free_gain` — for `059003.000613` the alternate at
  +28.87 cm⁻¹ was reported with `n_free` = 9 when only 4 of the 9 could have been recorded at all,
  and that inflation alone was carrying the row over the `n_free` ≥ 4 bar. Support
  that is entirely blended is not support: a prediction can be dumped on an already-explained
  feature almost anywhere, and the blend branch of the formula gives it credit for doing so.
  `top_share`, the largest single row's share of the positive evidence, throws out the positions
  whose whole case is one lucky line.
- **Vacancy.** A preferred alternate can be read two ways — the level moves there, or an *unknown*
  level sits there and the lines are its. The second reading needs the calculation still to have a
  level of the right J and parity spare at that energy, and `n_vacant` counts them: unfound rows of
  `IDEN2/enlev.dat` (no `*` in columns 39–40), of the level's own J and parity, within one
  configuration window of the alternate. **Where `n_vacant` = 0 the new-level reading is dead
  however good the lines look.** This is the only test in the module that comes from the theory
  rather than from the line list, and the only one that can rule a position out outright. It is a
  real filter: at the alternate position of `059003.000305`, eight unassigned, well-centred lines
  with coherent intensities looked like an undiscovered level, and the nearest unfound calculated
  level of J = 9/2 is 6900 cm⁻¹ away in the odd system and 18900 in the even one. Across the whole
  run, 50 of the 57 preferred alternates have `n_vacant` = 0.

What survives is sorted into six **dispositions**, named rather than ranked because they are six
different pieces of work:

| `action` | what it means | what to do |
|---|---|---|
| `interchange` | the alternate lands on another level of the run (`\|near_dE\|` < 0.5) | `level_interchange.py`, which weighs evidence this scan does not use |
| `refit` | the alternate is a fraction of a wavenumber away (`\|dE_alt\|` < 2), **the level keeps more than half of the recorded lines it is assigned now** (`n_kept` > `n_own`/2) **and more than half of their light** (`kept_light` > 0.5) | nothing is re-identified; some accepted line is dragging the LOPT fit off the position the level's own lines want |
| `top line` | the same short move, keeping most of the lines by count but **half or less of the light** (`kept_light` ≤ 0.5): the line it drops is the level's strongest | look at that line — its wavenumber, whether it is a blend, whether the identification is right — not at the position |
| `relocate` | a free position, ≥ 4 free lines worth ≥ 8, no line carrying more than 45 % of the case, `look` ≥ 3, **and `ln R` positive there** — convincing in its own right, not merely better than where the level is now. How far away it is does not enter | a candidate revision |
| `weak` | a move whose support is thin | leave it until the region around it is settled |
| `no support` | one free line, or none | no case at all |

**What separates a refit from a relocation is not distance — it is whether the level keeps its
lines.** `n_own` counts the recorded lines the level is assigned now and `n_kept` how many of them
the alternate still matches, and the two decide which of two opposite pieces of work the level
needs. A level that keeps most of them stays where it is: one bad assignment is pulling the
least-squares fit a fraction of a wavenumber off the position its own lines want, and correcting
that is a refit and nothing more — `059003.000457` keeps six of its seven at a position 0.51 cm⁻¹
away. A level that keeps none of them has to be re-identified however short the move: every one of
its present assignments is dropped, another set is made, and the level is entered in a ledger at a
new position — `059003.000617` keeps **none** of its three at a position 1.32 cm⁻¹ away, and
`059003.000602` two of its four at 0.79. Sorting the two by distance alone put them under the same
heading and told the reader to do the opposite of the work required. With the test applied, the 7
`refit` rows of the 2026-09-10 run become 5 refits and 2 relocations, both of which land in `weak`
— they are moves, and the report now prints that group so they are not lost.

**Counting the lines is not enough: the one being given up is usually the one that matters.**
`kept_light` weighs what `n_kept` counts — the share of the level's own observed light the
alternate still matches, each feature's measured intensity split among its accepted components by
branching fraction, so that a level is never credited with the whole of a blend. `059003.000457`
is the case that forced it: six of its seven lines survive a move of 0.51 cm⁻¹, which reads as one
bad line dragging the fit — but the line it drops is 49280.699, its strongest by an order of
magnitude (`kept_light` = 0.155), and the Ritz mismatch the move removes is the mismatch of the
level's brightest branch. Removing it improves the fit by construction and leaves the level
deprived of the transition that most defines it, which is a thing to explain and not a thing to
discard. Such a row is now called `top line` instead of `refit`, and what wants opening in IDEN2
is that line, not the position.

**The rows are ordered best-first**, so the report can be read from the top: with `--audit`,
relocations before interchanges before refits before the weak ones, and within each group the
largest `look`; with `--scan` alone, the largest gain; with neither, the weakest positions first.
The csv keeps that order.

**One caution, because it is the easiest mistake to make with this tool.** Levels share lines, so
every `ln R` is conditional on the rest of the level list, and a revision changes the verdict on
its neighbours. **Accept them one at a time and re-scan.** `059003.000604` is the demonstration:
with `059003.000371` at its old position the adopted energy of 000604 wins, and only after 000371
moves does the alternate at 136728.75 become the best position in the whole window.

Usage and outputs:

```bash
python level_positions.py                        # every level at its adopted energy
python level_positions.py --detail 059003.000271 # the per-transition table for one level
python level_positions.py --detail 059003.000613 --at alt   # ... at the alternate the audit proposes
python level_positions.py --scan                 # + the alternate-position scan
python level_positions.py --audit                # + what each alternate rests on and what to do
python level_positions.py --scan --alt-drop 3    # a stricter definition of "alternate"
python level_positions.py --audit --use-firm     # + skip the levels already settled (fast)
python level_positions.py --lopt LOPT_output_lines.txt   # judge a hand-revised run
```

**An empty `action` is not a clean bill of health.** `action` answers one question only — is
there somewhere better to put this level? — and it is blank whenever the scan found no alternate.
Whether the level is supported *where it is* is a separate question, and `ln R` answers it. The
two come apart in both directions. A level with `ln R` = 12 and `action` = `refit` is in good
order where it stands; the alternate is a worse position (`gain` negative) and the row is only
telling you that an accepted line is pulling on the fit. A level with `ln R` = −0.54 and a blank
`action` — `059003.000322` is the current example — is the hardest case in the report: the lines
do not support it where it is, and the scan finds nothing better anywhere in the 366 cm⁻¹ the
calculation allows, so whatever is wrong lies outside that window or in the lines assigned to it.
The report prints those levels in a section of their own, after the scan table, for exactly this
reason. The levels worth opening in IDEN2 are the union of the two lists: everything the audit
prints (it already restricts itself to alternates the lines actually prefer) and everything with
`ln R` below about 3, with or without an action.

**The registry of settled positions.** Most of a run is spent re-proving what was proved last
time: a level with `ln R` above 3 and no alternate anywhere in its window is not in doubt, and
about 485 of the 594 levels are in that state. Every `--scan` or `--audit` run writes them to
`level_positions_firm.csv` — `level_id`, `E`, `ln_R`, `scan_width` and a **fingerprint** of
everything that level's scan depends on: the interval scanned and its step, `--alt-drop`, the
level's own energy, every prediction it makes (partner, partner energy, partner uncertainty,
predicted intensity), the measured line list, the constants fitted to the run, and the accepted
assignments lying anywhere its predictions can reach as it moves across the window. `--use-firm`
skips the scan for a level whose fingerprint is unchanged.

The fingerprint is deliberately **local**. A run in which one relocation has been accepted differs
from the previous one in that level, in its partners and in the lines near it; a fingerprint taken
over the run as a whole would invalidate all 594 entries and the registry would never save
anything. As written, accepting a relocation costs a re-scan of the levels that share lines or
partners with it, and of nothing else. A second guard sits behind the fingerprint: `ln R` at the
adopted position is computed for every level on every run whether it is scanned or not — that is
one evaluation, not a scan — so an entry is honoured only if that value still comes out within 0.5
of the one recorded with it. The `scanned` column says which levels were taken from the registry
(`registry`) and which were scanned (`this run`). `--no-firm` ignores the registry altogether, for
a result that must not depend on any earlier run.

`--detail` prints the level's IDEN2 index and label, and the same for every partner, so the table
can be read beside IDEN2 without looking each partner up by hand. It lists only the rows worth
looking at — the observable predictions and any row that found a line — and says how many faint
predictions it left out; `--detail-all` lists those too. The `what` column names what each match
is: `free` (observable, well-centred, on an unclaimed feature — these are exactly the `n_free` of
the audit), `blend` (the feature is already explained), `off` (matched but too far out to count as
free), `faint` (the prediction could not have been recorded, so the match is a coincidence). By
default the table is evaluated at the adopted energy; `--at alt` evaluates it at the best
alternate the scan finds, which is the position an audit row proposes moving to and the only one
worth inspecting when the report suggests a move. `--at 138851.2` takes any energy.

`level_positions.csv`, one row per level:

| column | meaning |
|---|---|
| `E` | the adopted energy, cm⁻¹ |
| `n_drop` | predictions discarded because the partner would have no accepted line left without the ones it shares with this level |
| `n_pred` | predicted transitions in range |
| `n_obs` | of those, the ones that could have been recorded (P ≥ `P_SEEN`) |
| `n_seen` / `n_miss` | observable predictions that did / did not find a line |
| `n_match` | every prediction that found a line, faint coincidences included |
| `n_claimed` | matched predictions on features something else already explains |
| `ln_R` | the total (see "How to read `ln R`") |
| `ln_R_match` / `ln_R_miss` | its split between the matched rows and the absences |
| `sum_lnG` | the intensity evidence alone — a useful diagnostic on its own: a large negative value means the level's lines are there but their brightnesses are not what theory expects |
| `scan_width` | width of the interval scanned, cm⁻¹ (`--scan`) |
| `n_alt` | local maxima within `--alt-drop` of the adopted peak. **This count is the verdict** behind the `question` flag |
| `ln_R_alt` | `ln R` at the best of them |
| `d_ln_R` | `ln_R` − `ln_R_alt`; **negative means the lines prefer the alternate** (`gain` is the same quantity with the clearer sign) |
| `dE_alt` | E(alternate) − E(adopted): the size of the move |
| `near_level` / `near_dE` | the run's level nearest the alternate, and the gap — see above |
| `question` | `?` when `n_alt` > 0 |
| `scanned` | `this run` when the window was scanned, `registry` when the level was taken unchanged from the registry of settled positions |
| `n_free` / `free_gain` / `top_share` | the free-line test (`--audit`) |
| `n_own` / `n_kept` | recorded lines the level is assigned now, and how many of them the alternate still matches — the test that separates a refit from a relocation |
| `gain` / `look` | the gain and its look-elsewhere correction |
| `n_obs_alt` / `n_seen_alt` / `n_miss_alt` | the observable-prediction counts the level would have at the alternate |
| `n_vacant` | unfound calculated levels of the same J and parity within one configuration window of the alternate |
| `z_alt` | (E(alternate) − E_calc)/W: how far the alternate is from where the calculation puts the level, in units of that configuration's own scatter |
| `action` | the disposition |

### Searching for a level nobody has found: `level_positions.py --unknown`

Everything above scans a level the run already has: an adopted energy, a set of accepted lines, a
`level_id`. A level of the **calculation** that has never been found has none of those. What it has
is a row in `IDEN2/enlev.dat` — a calculated energy `E_calc`, a `J`, a configuration label — and a
place in Cowan's transition list, which gives it transitions to every other calculated level. The
ones whose other end **has** been found are predictions with a known partner energy: at any trial
energy E they predict a wavenumber, and a wavenumber is all `ln R` needs. So the same likelihood
ratio can be evaluated for a level that does not exist yet, and **scanning it across the window the
calculation allows is the search**. Where `ln R` > 0 the recorded lines are more likely with a level
at that energy than with nothing there; where it is largest they are most likely.

This is the scan done by hand in IDEN2, with the level list scrolled to a trial position and the
predicted transitions checked against the plate list one at a time. The arithmetic is the same
arithmetic. What it adds is that every position in the window is tried instead of the ones a person
has the patience for, and that the **absences** count against a position as well as the coincidences
count for it.

**Three things differ from the scan of a level the run already has, and all three make the answer
harsher rather than kinder.**

- **The level owns no accepted line.** Every feature that an existing identification already claims
  therefore counts against it in full through `C`. To be believed, an unfound level must explain
  lines the levels already found have left alone — a stricter test than a known level faces.
- **There is no adopted position to compare against**, so there is no `gain`. `ln R` itself is the
  verdict.
- **The look-elsewhere problem is the whole window**, not the distance between two candidates: a
  scan returning *n* positive maxima has had *n* chances at every one of them, so
  `look` = `ln R` − ln *n* is what a position is worth once the search that found it is paid for.

The window is `E_calc` ± 3W, W being the rms of E_obs − E_calc over the **found** levels of the same
configuration — the same interval a known level of that configuration is scanned over, and for the
same reason: how far the calculation can be wrong is a property of the configuration, not of the
level. The six configurations with no found level take the list-wide 132 cm⁻¹.

The predicted intensities are computed once at `E_calc` and held fixed as the scan moves.
`I = C·gA·(ν/1e8)·exp(−E_up/kT)` does depend on the trial energy through both ν and the Boltzmann
factor, but over a few hundred wavenumbers out of a hundred thousand that is a few per cent, against
an observed-to-predicted scatter of a factor of three — and it is exactly what the scan of a known
level already does, so the two are comparable. The `gA` values, the matching of Cowan's level
numbering to IDEN2's, and the per-configuration widths are all read through `unfound_levels.py`, so
there is one answer to each of those questions and not two.

**Each candidate position is graded by the audit's own constants**, because the evidence is the same
kind of evidence and the scan that found it had the same freedom to look: `no support` when fewer
than two free lines support it or one line carries more than 60 % of the case, `firm` when it would
pass every test a firm relocation passes, `weak` in between.

**What it finds.** IDEN2 row 742 — `f25f ~3F4G`, J = 5/2, calculated at 116406.7 and the top of
`unfound_levels.py`'s list — has 95 calculated transitions to found levels, 40 of which could have
been recorded. The scan of its 621 cm⁻¹ window returns two positions with `ln R` > 0:

```
         E  ln_R   look      dE      z  n_obs  n_match  n_miss  n_free  free_gain  top_share verdict
116327.610 9.180  8.487 -79.090 -0.765     40       20      20      10     23.200      0.132    firm
116327.070 0.156 -0.537 -79.630 -0.770     40       19      21       9     14.832      0.196    weak
```

The first rests on ten free lines worth 23.2 in `ln R`, no one of them carrying more than 13 % of
the case, and sits 0.77 configuration widths below the calculated energy. `--at 116327.610` lists
them: all ten are transitions to levels of 4f².5d, spread from 987 to 1396 Å, each matched to a
recorded line within 0.7 cm⁻¹ on a feature nothing else claims.

The counter-example is as useful. IDEN2 row 1 — `f5d6d`, one of the six configurations with no found
level — returns **269** positions with `ln R` > 0, and every one of them is `no support`: `n_obs` is
zero at all of them, so not one of the level's 45 transitions could have been recorded, and the
positive `ln R` comes entirely from faint coincidences: the six best carry `top_share` between 0.50
and 0.92, which is to say one line is the whole of the case for each of them. That
is what a 793 cm⁻¹ scan finds when there is nothing there, and the report says so in as many words.
`look` is negative at every one of them.

Usage:

```bash
python level_positions.py --unknown 742             # the candidate positions for one row
python level_positions.py --unknown 742 913 986     # several at once
python level_positions.py --unknown 742 --at 116327.610   # the transitions at one of them
python level_positions.py --unknown 742 --top 10    # only the best ten positions
python level_positions.py --unknown 742 --min-ln-r 3      # a higher bar than "better than nothing"
```

A run takes about three seconds after the 1.6 s the likelihood takes to build.

**What this does not do.** It says which energies the recorded lines want a level at; it does not
say the level is there. A position here is conditional on the rest of the level list exactly as an
alternate position is — the lines it takes are lines its neighbours could take instead — so a
candidate has to be checked against IDEN2 and accepted one at a time. And it searches only where
`unfound_levels.py` says there is something to search on: a level with fewer than two transitions
that could have been recorded is not searchable at all, one line being enough to fit any energy, and
the report refuses to pretend otherwise.

### Limits of the validation (to be stated alongside the results)

- **Recovery, not physical proof.** A small ΔE certifies that the accepted lines reproduce the energy encoded in Wyart's input value — i.e. that his identifications were recovered. If Wyart himself was misled by chance coincidences, our run re-finds the same coincidences with a small ΔE; only the intensity pattern (and physics arguments: theory, g-factors, term structure) can catch that case.
- **Theory-fault sensitivity.** A bad pattern score can mean a spurious level *or* a level whose theoretical description is wrong; levels whose accepted lines are mostly uncovered by theory (intensity grade `G`) get weak-quality pattern verdicts. In the converged run all 5 marked levels do have observable predictions (in earlier runs a sizeable part of them did not, and there `p_spur` rests on the energy shift and the support count alone). A related caution: for a few marked levels every observable prediction lies within a factor ~3 of the local noise — there a bad pattern score means little, because the predicted intensities themselves scatter by a factor ~2.5 against the observed ones.
- **Decoys can accidentally be real.** A displaced decoy may land on a real, previously unknown level (most likely a neighboring-J member of the same term at high energies) or on the true position of a level misassigned in the underlying Cowan-code fit. Both effects make some decoy "false positives" actually real, so the decoy-based false rates are slight overestimates — the bias is in the safe direction. The fraction of supported decoy trials showing three or more of their top-10 predictions accepted (~6%) is an upper bound on this contamination.
- **A bad intensity pattern can be hidden by the folding.** `p_spur` folds the energy-shift evidence with the pattern evidence, and a level whose strongest predicted branch is *masked* by a nearby stronger line has that half of the pattern evidence withdrawn — correctly, since a masked branch says nothing. But a level can then keep a very poor pattern score and still come out with a small `p_spur` and no questionable mark. In the run of 2026-09-05 the worst intensity pattern of all 208 tested levels belongs to `059003.000572` (`pattern_V` = 0.122, 4 of its top 10 predictions accepted), which carries `p_spur` = 0.033 and is not marked; its configuration, `4f.5d.6p`, is the one whose calculated positions deviate most from the observed ones in the Cowan-code fit, so a misidentification of the level with the wrong theoretical partner is exactly what one would expect there. **`pattern_V` is worth reading directly, not only through `p_spur`.** The specific failure — two levels of the same parity and `J`, close in energy, whose theoretical identities have been interchanged — is what `level_interchange.py` looks for (next section); it is a separate test because an interchange leaves every energy right and so is invisible to ΔE.
- **Lines omitted as blends with other ionization stages.** Sugar assigned lines to Pr II, III, or IV by comparing exposures at different degrees of excitation; a Pr III line nearly coinciding with a stronger Pr II or Pr IV line could not be assigned confidently and was omitted from his Pr III list. Such omissions are invisible to the masking check (the search covers only Sugar's Pr III lines), so some "missing" strong predictions are excusable in a way the automation cannot see. Statistically the effect is absorbed — the old (genuine) reference levels suffer the same omissions, so the pattern-score comparison between the classes stays fair — but the per-level adjudication is conservative: a questionable mark that an expert would clear on this ground stays retained. Automating this excuse would require Pr II and Pr IV line lists.

### Are the files in step: `check_sync.py`

The same facts — which observed line belongs to which pair of levels, and where each level
sits — are held in five places, and nothing keeps them in step automatically:

```
Pr3_lines.xlsx                    the measured lines (the only real input)
    |  classify_lines.py, obeying line_decisions.csv and revised_level_energies.csv
    v
line_classifications.csv / .xlsx  every candidate transition, accepted or not
    |  make_LOPT_input.py
    v
LOPT_input_lines.txt              the accepted ones, weighted
    |  lopt.bat
    v
LOPT_output_lines.txt             the fit
LOPT_output_levels.txt            the optimized level energies
    |  by hand, rarely
    v
IDEN2/enlev.dat                   the level list IDEN2 shows on screen
IDEN2/trans.dat                   its predicted transitions and their assignments
```

A step skipped anywhere leaves two of them disagreeing, and the disagreement is silent:
every tool downstream keeps working, on stale numbers. The symptom is a report that
recommends something already done — `level_positions.py` proposing a move to the energy the
level was moved to yesterday — which, read on its own, looks exactly like a recommendation
worth acting on. `check_sync.py` compares all five and says where they differ. It changes
nothing.

**What is a mismatch and what is only drift.** The IDEN2 files are brought up to date by
hand and only when there is reason to, so their level energies are expected to lag the
current fit a little; a hundredth or two of a wavenumber is the ordinary state of affairs.
Two thresholds separate that from a real difference: `--drift` (0.5 cm⁻¹ by default) below
which an energy difference is the expected lag and is only counted, and `--gross` (5.0)
at which the two files are describing different positions rather than one position known to
different precision. A third, `--wn-tol` (0.05), is not a tolerance on the physics but on
the printing: `make_LOPT_input.py` writes three decimals, LOPT echoes each wavenumber to
the precision its own uncertainty warrants (121279.416 comes back as 121279.42), and the
decision ledger is written by hand from a printed list. Everything else — a line accepted
in one file and absent from another, a ledger verdict the classification does not obey, a
transition assigned in IDEN2 that the classification has withdrawn — has no tolerance at
all and is reported item by item.

**Why the level pair is the key and the wavenumber is not.** Two files are compared
transition by transition, matched first on the pair of levels and then, within that pair, on
the nearest wavenumber. Keying on the wavenumber alone reports thousands of differences
that are only rounding; keying on the pair alone is wrong the other way, because the same
two levels can be candidates on several observed lines — moving an assignment from one line
to another, which is how a correction is made, leaves the pair on both, and both ledger
rows are then obeyed. Two *accepted* rows for one pair would be a real contradiction, since
one transition can be on one line only, and that is reported.

**The three severities.** `ERROR` — the files contradict each other, and a report built on
them can be wrong in a way no reading of it would reveal. `WARN` — they disagree in a way
that is expected, or in one artefact that can simply be rebuilt. `ok` — checked, nothing to
say. The exit status is 2 if anything is an error, 1 if anything is a warning, 0 if
everything is clean, so the script can gate a pipeline.

The ten checks, in the order they are made:

| # | what is compared | what a difference means |
|---|---|---|
| 1 | every artefact against the files it is built from, by modification time | the cheapest check, and it usually explains all the others at once |
| 2 | `line_classifications.csv` against `.xlsx` | they are written together by `classify_lines.py`, so a difference was made afterwards — an edit in Excel, or one file left over from an earlier run |
| 3 | the classification against `LOPT_input_lines.txt` | what has to agree is the fit: every accepted classification weighted, nothing else weighted, and each carrying the right share of its observed line |
| 4 | `LOPT_input_lines.txt` against `LOPT_output_lines.txt` | LOPT must have been run on the transitions file that is on disk now |
| 5 | LOPT's line output against its level output | different level energies in the two mean they are from different runs |
| 6 | `revised_level_energies.csv` against the fitted levels | a revision written down and never carried out — or one made and then undone by a later run started from an older classification |
| 7 | `IDEN2/trans.dat` against `IDEN2/enlev.dat` | trans.dat carries a copy of every partner energy, found flag and predicted wavenumber, and IDEN2 rewrites them together; a copy that no longer follows means every wavenumber trans.dat shows is stale |
| 8 | `IDEN2/enlev.dat` against `LOPT_output_levels.txt`, joined through `IDEN2/IDEN_level_ids.txt` | the energy of every level, and the membership of the two lists both ways: levels starred in enlev.dat with no entry in the map, levels of the fit not marked found |
| 9 | `IDEN2/trans.dat` against the accepted classifications | **the check the rest exists to make**: a transition IDEN2 shows as identified but the classification has withdrawn looks, on the screen, exactly like one that is still accepted |
| 10 | `unstable_candidates.csv` and `line_decisions.csv` | stability rather than staleness — see below |

**Stability.** Two things can be unsettled rather than merely out of date. An *oscillating*
assignment is one `classify_lines.py` accepted on one pass and rejected on the next until it
blacklisted the pair; the run converged, but only because that assignment was taken out of
the argument, and which way it should have gone is still undecided. It is reported as a
warning with the acceptance states that made it an oscillator, and the action is to settle
it by hand and record the verdict in `line_decisions.csv`. An *unobeyed ledger row* is
worse and is an error: the analyst has ruled on a line and the classification does not show
the ruling, which means it was written before the row, or by a run that did not read the
ledger. A rejection obeyed on the line it names, while the same pair is accepted on another
line, is neither — that is how an assignment is moved, and it is reported as a warning only
so that the move is visible.

**The problem lists.** Every list of transitions is written so that it can be worked
through at the screen without a second lookup. The rows are ordered by the energy of the
upper level and then of the lower one — the order the same transitions come in in IDEN2 and
in the workbook — and each row carries both levels with IDEN2's own level numbers beside
them (the numbers typed into IDEN2, from `IDEN_level_ids.txt`), both energies, and three
words saying where the transition stands:

| column | value | meaning |
|---|---|---|
| `.xlsx` | `accepted` | `line_classifications.xlsx` accepts it on this line |
| | `rejected` | it is a candidate on this line and was turned down |
| | `elsewhere` | it is a candidate, but on other observed lines only |
| | `missing` | it is not a candidate anywhere |
| `LOPT` | `included` | on this line in `LOPT_input_lines.txt` with a weight and no `P` flag |
| | `not incl. P` | in the file, but flagged predicted or given zero weight — either way it does not pull on the levels |
| | `elsewhere` | in the file on other observed lines only |
| | `not incl. --` | not in the file at all |
| `IDEN2` | `assigned` | `trans.dat` identifies it with this line |
| | `elsewhere` | `trans.dat` identifies it with another line |
| | `not assigned` | IDEN2 predicts it and has no line for it |
| | `no trans row` | IDEN2 does not predict it at all, which for a real pair of levels means `trans.dat` is stale |

Every one of the three is answered for the observed line the row names, never for the pair
of levels in general. The same pair can be a candidate on several lines — moving an
assignment leaves the old row in place to be rejected, which is how the pipeline is made to
turn down an assignment it would otherwise keep — so a status read off another line would
print a row that is in perfect order as a problem.

**The weights.** A row that both files hold can still say opposite things, because being in
the transitions file is not being in the fit: a `P` flag or a zero weight takes it out. An
accepted classification that is in the file without weight is therefore listed with the ones
that are not in the file at all, and a weighted row the classification does not accept with
the ones the classification has never heard of. Only when both files put the line in the fit
is there a weight to compare, and what is compared is the *share of the observed line*, not
the number in the weight column: LOPT normalises the weights of one line to sum to one, so a
blend assigned by hand with the calculated intensities pasted straight out of `Icalc.xlsx`
fits exactly as `make_LOPT_input.py`'s own fractions would and looks nothing like them on
paper. The tool repeats that arithmetic — the whole line to a single accepted
classification, a blend divided between its components in proportion to `calc_intens` — and
reports a share that is out by more than one per cent, printing both.

**What the lists leave out.** A transition that carries no weight in the fit and no
identification in IDEN2, and that the classification does not accept, says the same thing
in all three files whichever of them happens to hold a row for it: it changes neither the
level optimisation, nor the next pass of `classify_lines.py`, nor the view in IDEN2. It
arises in both directions — `make_LOPT_input.py` writes every classified candidate, the
rejected ones at zero weight, so settling a classification after the transitions file was
built leaves rejected candidates in the workbook alone, and withdrawing one leaves a
`P`-flagged row in the transitions file alone — and in both directions it is counted in a
warning and kept out of the detailed lists. The rebuild that puts the two files literally
in step belongs at the end of the analysis, when fully synchronised files are wanted for
the publication tables. The same transition weighted in the fit, or identified in IDEN2, is
a real disagreement and stays in the list.

**Where the report goes.** The terminal shows the first `--list` rows of each list, enough
to see the shape of the trouble at once; the complete report, every item of every list, is
written to `sync_report.txt`. A list cut off after a dozen rows cannot be worked through,
and the rows that were cut are then invisible.

**The suggested actions** are printed last, in the order the pipeline runs rather than by
severity, because doing a later one first only means doing it again after the earlier one:
settle the unstable assignments and the ledger, re-run `classify_lines.py`, rebuild the
LOPT input, re-run LOPT, carry the changes into IDEN2, and re-run the validation reports
last of all. One thing the tool deliberately does *not* advise is rebuilding
`LOPT_input_lines.txt` to make a mismatch with the classification go away. Rebuilding it
does exactly that — makes the list go away — while every question in the list stays open
and is now untracked: the transitions have to be looked at in IDEN2 one by one and settled
in `line_decisions.csv` first.

Usage:

```bash
python check_sync.py                     # the report, and sync_report.txt
python check_sync.py --quiet             # only the checks that found something
python check_sync.py --list 40           # up to 40 items under each finding on screen
python check_sync.py --out mine.txt      # write the complete report elsewhere
python check_sync.py --drift 0.1         # a stricter idea of "up to date"
python check_sync.py --csv sync.csv      # the findings as a table as well
```

Like the swap tools, it looks for each file in the directory it was run from and falls back
to the directory holding the scripts, so it can be run inside an iteration folder
(`iter22/`) against that folder's own LOPT files while still reading the one copy of the
ledgers.

### Bringing IDEN2 up to date: `sync_IDEN2.py`

`check_sync.py` says when IDEN2's files have fallen behind. `sync_IDEN2.py` catches them
up. It rewrites both of them from the current fit and the current calculated intensities:

* **`IDEN2/enlev.dat`** — the measured energy of every found level and its uncertainty are
  taken from `LOPT_output_levels.txt`. The uncertainty is the larger of LOPT's `D1` and
  `D2tot`, rounded to a thousandth of a wavenumber, which is the rule the file already
  follows for all 594 levels. The observed-minus-calculated column is recomputed from each
  level's own calculated energy. A level that has **not** been found keeps its calculated
  energy as its adopted energy — this run knows nothing better — but its uncertainty is
  rewritten too, and there the column means something else entirely.
* **`IDEN2/trans.dat`** — every row is rewritten with the partner's new energy, the new
  predicted wavenumber (the difference of the two energies exactly as `enlev.dat` writes
  them) and a new intensity code.

**The uncertainty column holds two different quantities.** On a found level it is the
uncertainty of a measurement: how well the fit knows where the level is, a few thousandths
of a wavenumber. On a level nobody has found there is no measurement to be uncertain about,
and the only honest reading of the column is a prediction — how far from its **calculated**
position the level is likely to turn out to be. IDEN2 uses it that way: it is the half-width
of the window it searches. That distance is a property of the configuration rather than of
the level, because the calculation describes some configurations far better than others, so
it is measured per configuration as the rms of E_obs − E_calc over the levels of that
configuration that have been found. It ranges from 24 cm⁻¹ for `f25g` to 404 cm⁻¹ for
`fd6p`, against 132 cm⁻¹ for the level list as a whole — a factor of seventeen, which is why
one number for all of them will not do. Six configurations (`5d3`, `f26g`, `f5d6d`,
`f5d7s`, `f6s2`, `p5f3d`) have not one level found in them; their 149 levels are left
exactly as they are, still carrying the placeholder 5000, because the list-wide rms would be
a claim about them that nothing supports. Both kinds of change — the found levels' `max(D1,
D2tot)` and the unfound levels' configuration rms — are listed in the run's report.

The intensity code `trans.dat` carries is

```
Icalc_IDEN2 = round(10 * ln(Icalc)) ,   Icalc = C * gA * (rwn/1e8) * exp(-Eup/kT)
```

a small signed integer — one unit is a factor 1.105 in intensity. `gA` comes from
**`tp_E1_no_trials.xlsx`**, the Cowan E1 transition list, which is the file `Icalc.xlsx`
itself was drawn from; `C` and `kT` come from `[intensity_model]` of the configuration and
are checked against `Icalc.xlsx` before anything is written. The codes now in the file were
computed long ago at kT ≈ 13100 cm⁻¹ with a normalisation to match: they rank the
transitions almost perfectly (rank correlation 0.9994 with the current `Icalc`) but they are
a factor 0.70 and a Boltzmann tilt away from the numbers the pipeline works in.

**An identification is never dropped.** A transition that carries an observed line stays in
`trans.dat` however weak its new intensity code makes it, and whether or not the Cowan table
still has a `gA` for it — which one of them does not: the hand-added row for Sugar's
identification of a line below the `gA` = 10³ s⁻¹ floor of the Cowan run, the one imputed
`Icalc` in `line_classifications.csv`. Its row is carried through with the code it has.

Which of the other transitions are listed is set by `--cutoff`, on the intensity code. The
default, −41, is the value that leaves the present list most nearly alone: on the first run
it added 420 transitions the old scale had cut off and dropped 47, out of 104272. Each step
of one removes about 700 more.

A level that ends up with no transition above the cutoff heads no block, and IDEN2 then
shows nothing below it. That is not a loss of information: IDEN2 displays only the
transitions with a code of 0 or more anyway, so a level whose whole block falls below the
cutoff had nothing on screen to lose. The report names those levels.

```bash
python sync_IDEN2.py                 # rewrite both files, keeping .presync copies
python sync_IDEN2.py --dry-run       # say what would change, write nothing
python sync_IDEN2.py --cutoff -45    # a longer list, down to weaker lines
python sync_IDEN2.py --report sync_IDEN2_report.txt
```

**Close IDEN2 first.** It holds its files open and rewrites them from memory when it exits,
which would undo everything the run has done.

The reader of the Cowan table is `cowan_gA.py`, a module of its own, because
`level_positions.py` needs it too: it is the only source of a calculated strength for a
transition of a level that has never been found, and so the only way to say which of the
659 unfound levels are worth looking for. It also carries the join between the three
numberings of the same levels — the calculation's `lid`, IDEN2's row number, and Wyart's
`level_id` — which it works out by aligning the two energy-ordered level lists and then
checking the result against all 594 levels whose row number and `level_id` are both known.

### Which unfound levels are worth searching for: `unfound_levels.py`

**What it answers.** Cowan's calculation gives Pr III 1253 levels. 594 have been found; the
other 659 have never been placed, and each of them is an energy the calculation predicts and
the line list has never been searched for. Searching for one is an afternoon's work in IDEN2,
so the question is which of the 659 to spend it on. A level can only be found through its
lines, so the answer is a count: **how many of the transitions the calculation gives it would
have been recorded on the plates**. One is never enough — a single line can be made to fit any
energy, so one coincidence is no evidence — and ten is a position that, if the level is there
at all, is heavily overdetermined.

**What makes a transition promising.** For an unfound level L at a trial energy E and a partner
M that *has* been found, at its measured energy E_M, the line would sit at ν = |E − E_M| with
the intensity `I = C·gA·(ν/1e8)·exp(−E_up/kT)`, the same formula and the same C and kT as
everything else in the pipeline. gA comes from Cowan's own transition list through
`cowan_gA.py`, which is the only place it can come from: `Icalc.xlsx` holds only transitions
whose two levels are both known, and by definition none of these are. That the two agree is
worth checking and was checked — recomputing the 29260 predictions of `Icalc.xlsx` from gA this
way reproduces the file to an rms of 0.2 % in ln I. The probability that such a line reached
Sugar's list is then the same `P_obs = c(λ)·D(z)` the rest of the analysis uses, and a
transition is **promising** when P_obs ≥ 0.5.

**The level's energy is not known, and the count says so.** How far a level turns out to be
from its calculated position is a property of the *configuration*, not of the level: W, the rms
of E_obs − E_calc over the found levels of the same configuration, runs from 24 cm⁻¹ for
4f².5g to 404 cm⁻¹ for 4f.5d.6p. So ν is uncertain by that much, and with it λ, the local noise
level and therefore P_obs. Every probability here is averaged over a Gaussian of width W
centred on E_calc — nine trial positions from −2W to +2W, weighted by the normal density — so
what comes out is a probability in the ordinary sense: the chance the line is on the plates
given both the physics and the ignorance about E. Summing it over a level's transitions gives
`sum_P`, the expected **number** of its lines that were recorded. The six configurations with
no found level of their own (`5d3`, `f26g`, `f5d6d`, `f5d7s`, `f6s2`, `p5f3d`) have no W; they
are listed anyway, with the list-wide 132 cm⁻¹ and a `list` in the `W_src` column, because a
ranking is a suggestion about where to spend an afternoon and not a number written into a file.

**What it finds.** Of the 659, every one has at least one calculated transition to a found
level, but

| promising transitions | levels |
|---|---:|
| none — nothing to search on | 232 |
| exactly one — not enough | 97 |
| two or more — worth a search | **330** |
| five or more | 106 |

and the search list is concentrated: `f27d` 62, `f25g` 50, `fd6p` 49, `f26d` 32, `f27p` 29,
`f26f` 21. The best of them is IDEN2 row 742, `f25f ~3F4G`, J = 5/2, calculated at 116406.7:
23 promising transitions and an expected 22.9 recorded lines, 21 of them to levels of 4f².5d
and 18 of those between 985 and 1125 Å.

**The count is blind, so it can be scored on the levels already found** (`--validate`). Nothing
in it looks at whether a line was ever assigned to a level — only calculated energies,
calculated gA, and the measured coverage and noise of the plates — so every found level can be
put back at its *calculated* energy and counted as though nobody knew where it was, and the
answer compared with the number of accepted lines it really has:

| `n_prom` | levels | median accepted lines | mean |
|---|---:|---:|---:|
| 0 | 1 | 5 | 5.0 |
| 1 | 11 | 3 | 3.0 |
| 2–4 | 88 | 4 | 4.2 |
| 5–9 | 138 | 7 | 7.6 |
| 10+ | 356 | 19 | 23.8 |

Rank correlation 0.925 for `n_prom` and 0.943 for `sum_P`. And 582 of the 594 clear the bar of
two promising transitions, which is that bar being checked from the other side: a level the
count would have called hopeless is hardly ever a level that was in fact found. The found
levels are not a fair sample — they were found *because* they had lines — so this says the
count ranks levels correctly, not that a level with `n_prom` = 10 will turn up.

**Limits.** Three, and they are all of the same kind: the count says the evidence would exist,
not that it can be recognised.

- **It does not know whether the lines are still free.** A promising transition may fall on a
  feature another level already claims. While the level could be anywhere in a window hundreds
  of wavenumbers wide there is no way to say which lines it would reach; `rho_med`, the density
  of recorded lines around its promising transitions, is the nearest thing to a warning — where
  it is high, a trial energy will find lines whether the level is there or not.
- **It does not know whether the partners are sound.** Every partner is a found level, but
  found levels differ in how firmly; `level_positions.py` is what says which is which.
- **It does not know whether the calculated strengths are right.** Observed intensities scatter
  against predicted ones by a factor of about three. P_obs carries that spread through the
  detection curve, but a transition predicted at ten times the noise can still be absent.

Usage:

```bash
python unfound_levels.py                   # the ranked table + unfound_levels.csv/.xlsx
python unfound_levels.py --detail 742      # the transitions of one level, by IDEN2 row
python unfound_levels.py --detail 742 --detail-all   # ... including the hopeless ones
python unfound_levels.py --validate        # score the count on the levels already found
python unfound_levels.py --min-prom 5 --top 0        # print the whole short list
```

`unfound_levels.csv`, one row per unfound level, best first:

| column | meaning |
|---|---|
| `idx` | the level's row number in `IDEN2/enlev.dat` — the key IDEN2 itself uses |
| `label` / `J` / `parity` / `cfg` | as `enlev.dat` writes them |
| `E_calc` | the calculated energy, cm⁻¹ |
| `W` / `W_src` | how far levels of that configuration turn out to be from their calculated position, and whether that came from the configuration itself or from the list as a whole |
| `n_pred` | transitions to found levels the calculation gives it |
| `n_obs` | of those, the ones with P_obs ≥ 0.2, the bar `level_positions.py` uses |
| `n_prom` | **the verdict**: transitions with P_obs ≥ 0.5 |
| `sum_P` | the expected number of its lines actually on the plates |
| `I_max` / `wn_I_max` | its strongest predicted transition and where it would fall |
| `rho_med` | recorded lines per cm⁻¹ around the promising transitions |

The `--detail` table names each partner by IDEN2 row, label and `level_id`, says whether the
unfound level is the upper or the lower of the two, and gives ν and λ at E_calc, the predicted
intensity, P_obs and the local line density. It is meant to be read beside IDEN2.

**This is the first half of the search.** It says where to point a search, not what a search
finds. Taking a level off this list and asking whether the lines actually put it somewhere — the
scan of `ln R` over the window, at trial energies rather than at an adopted one — is the other half,
and it is `level_positions.py --unknown`
([Searching for a level nobody has found](#searching-for-a-level-nobody-has-found-level_positionspy---unknown)).

### Putting a newly found level in: `insert_new_level.py`

`level_positions.py --unknown` finds a position where a calculated level nobody
has ever found would explain a group of observed lines. Accepting that position
is a decision made by eye, in IDEN2. Everything that has to follow it is
mechanical, and it used to be a dozen edits by hand across six files —
`new_levels.txt`, `IDEN2/IDEN_level_ids.txt`, `LOPT_input_lines.txt`,
`line_decisions.csv`, `IDEN2/trans.dat`, `IDEN2/enlev.dat` — with three programs
and LOPT to be run in between, in the right order. A step forgotten leaves two
of them describing different identifications, and the disagreement is silent.

```
python insert_new_level.py --iden2-row 742
python insert_new_level.py --iden2-row 742 --yes \
    --reject 91856.116="May add to pub line list as masked"
python insert_new_level.py --iden2-row 742 --yes --rebuild \
    --reject 91856.116="May add to pub line list as masked" \
    --accept 95033.381 --accept 93276.982 --accept 90917.831="better CoG"
```

The first form writes nothing: it prints the proposal table and stops. What it
does with `--yes`, in order:

| | |
|---|---|
| **A** | every file the run can write is tested — `output_files.require_writable()`, so that a workbook open in Excel stops the run at the start rather than at the end — and then copied byte for byte into `insert_new_level_backup/` |
| **B** | the level is added to `new_levels.txt` if it is not there: the next free identifier, the adopted energy and J of `enlev.dat`, the parity of the calculated level, `iden2_row` and `cowan_lid`. `IDEN2/IDEN_level_ids.txt` gets the identifier against the row |
| **C** | the lines: what is already marked in `IDEN2/trans.dat` for that block, taken as given — and when there is any such mark, that is **all** the level gets; the tool's own proposals, only for a level with no mark anywhere; and `--reject WN=reason`, which writes a verdict and assigns nothing |
| **D0** | `lopt.bat LOPT.par` on the untouched input, to have a fit to compare with |
| **D** | one record per accepted assignment goes into `LOPT_input_lines.txt`, and every unflagged record of every observed line the run touches is re-weighted together |
| **E** | `lopt.bat LOPT.par` again, then `RSS/degrees_of_freedom` before against after, then the four-sigma Ritz check |
| **F** | `classify_lines.py` and the ledger rows, repeated until nothing new has to be written. **`make_LOPT_input.py` is not run** unless `--rebuild` asks for it |
| **G** | `IDEN2/trans.dat` is made to show what the classification accepts **for this level**; anything else it has newly accepted is named and left alone, unless `--accept WN` names its line, and then it is written too |
| **H** | `check_sync.py` must report **no errors**; warnings are printed and ignored |
| **I** | `sync_IDEN2.py` finishes the job |
| **J** | the report, on screen and in `insert_new_level.log` |

Six things about it are worth knowing.

**A level already marked up in IDEN2 gets nothing added to it.** The lines of a
level that has been gone through on the screen, line by line, against the
plates and the branch structure, are the whole of the answer: a line the
analyst passed over was passed over for a reason none of the arithmetic here
can see. So when `trans.dat` holds any mark for the level, the tool takes those
marks and proposes nothing of its own; only a level straight out of
`level_positions.py --unknown`, with no mark anywhere, is assigned
automatically. `--propose` asks for proposals beside the hand marks and
`--no-propose` for none ever.

**`make_LOPT_input.py` is not run.** It rebuilds `LOPT_input_lines.txt` out of
the whole classification, and so puts into the fit every assignment the
classification accepts — including ones this run never proposed and nobody has
looked at in IDEN2. The fit this run makes holds exactly the records step D put
into it. The same applies in step G: adding a level to a line that already had
components changes the intensity accounting on that line, and the
classification can turn one of the other components from rejected into
accepted on the strength of it. Those are named, not written, and the run is
**held** there — everything it wrote stays in place so the assignments can be
looked at in IDEN2, `check_sync.py`'s findings are printed, `sync_IDEN2.py` is
not run, and `insert_new_level.py --undo` puts the whole run back if the
verdict goes the other way. `--rebuild` asks for the rebuild explicitly, for
when those verdicts have already been given.

Naming such an assignment is not enough to decide it, so each one is printed
together with **every** transition the classification puts on the same observed
line — the new level's own among them, marked `*`, since it is the reason the
rest of the line may have moved. For each: the predicted intensity, `dif_O-C`
(the observed wavenumber less the Ritz wavenumber the two levels imply, in
cm⁻¹) with its previous value if this run changed it, the grade, the verdict
with its previous value if this run changed that, and the reason
`classify_lines.py` gave. The head of the block carries the line's own observed
wavenumber, uncertainty, intensity and character, and the rows of
`line_classifications.xlsx` the block covers, so the workbook can be opened at
the right place for anything the printout does not settle. This is what shows,
for instance, that 95033.381 turned from rejected to accepted with its
`dif_O-C` unchanged at −1.025 cm⁻¹: what changed was the centre of gravity of
the blend it now forms with the new level's line, which is exactly the kind of
reasoning the verdict has to be judged on.

**`--accept WN` is how the verdict is given.** Neither `--yes` nor `--rebuild`
decides these assignments — `--yes` authorizes the run to write what it
proposed itself, `--rebuild` only says the fit's input may be built from the
whole classification — so a run that repeats them holds again, in the same
place, with the same list. `--accept 95033.381`, repeatable and taking an
optional `=reason`, says that the assignments printed for that observed line
have been looked at and are good: each one the classification accepts on that
line and this level is not part of gets an accept row in `line_decisions.csv`
and a mark in `IDEN2/trans.dat`, exactly as the level's own do, and the run
goes on through `check_sync.py` to `sync_IDEN2.py`. A component of that line
the classification *rejects* is left rejected — `--accept` adopts what was
shown, never everything on the line — and one IDEN2 already shows is left
alone rather than given a ledger row, since a published identification nobody
questioned is not this run's decision to claim. A line that should not be taken
goes in as `--reject WN=reason` instead. Because these assignments enter the
fit only through the classification, the run that adopts them wants `--rebuild`
as well; the held run prints the exact command line to use.

**`RSS/degrees_of_freedom` is reported before and after.** This is LOPT's own
measure of how well the whole fit holds together: the sum, over every observed
line the fit uses, of the squared difference between the observed wavenumber
and the one the fitted levels imply, each divided by that line's own
statistical uncertainty, per degree of freedom (lines less levels determined).
It is about 1 when the uncertainties are honest and the identifications are
right, and it grows when a line is put where it does not belong — so it is the
one number that says whether a new level has been paid for by making
everything else fit worse. It stands at 1.16. `--max-rss R` stops the run if
it comes out above `R`; by default it is only reported.

**A line shared with another transition has no residual of its own.** LOPT
fits such a line as one blended feature — its centroid against the
intensity-weighted mean of the components' Ritz wavenumbers — so there is one
residual for the feature and `_` in each component's O−C column. Those
assignments are in the fit and carrying their weight; the four-sigma check
simply has nothing to test them with, and they are listed separately rather
than passed over in silence. Six of level 742's sixteen lines are of this kind.

**Where the level goes.** `enlev.dat` carries a measured energy for a row only
once the position has been accepted in IDEN2 and written there. Until then the
row holds the calculated position, which is nowhere near good enough to search
on, and `--energy E` says where to put the level. Given as well as a measured
energy, it overrides it.

**A line much stronger than predicted is never proposed.** It is listed as left
free, with the factor and how many standard deviations of `ln(I_obs/I_calc)` it
amounts to. A line whose strength the new level cannot account for belongs to
some other transition; taking it would both misplace this level and hide the
real owner. The 91856.116 cm⁻¹ line of level 742 — 58 times its predicted
intensity, 4.5 σ — is the case the rule was written for. The threshold
(`--strong-sigma`, 3 σ) is deliberately generous, and so is the window a line
must fall in before the tool will propose it at all (`--propose-window`,
1 cm⁻¹, against the 2 cm⁻¹ it looks in): the tool is not trying to reproduce
the analyst's judgement, only to decline the cases that are not simple. Nothing
here overrules a mark made by hand in IDEN2, which is read as a decision
already taken.

**Anything that goes wrong puts every file back.** A Ritz residual above four
times a line's own uncertainty, an unexplained error from `check_sync.py`, a
hand mark the classification will not accept even with the ledger rows the run
wrote, a bug, a Ctrl-C: the culprits are named and every backed-up file is
restored byte for byte. (A *hold* is the one stop that does not restore, because
the files have to stay as they are for the decision to be made on them;
`--undo` restores them.) `--dry-run` does the same thing on purpose — it runs the whole sequence,
programs and all, and then puts everything back, which is the only way to see
what the classification and the fit will say before committing to them.

**The weights of a shared line are all recomputed, not just the new one.** When
a new transition joins an observed line that another transition already had,
the line has to be divided again, in proportion to the calculated intensities;
the component that was there alone loses exactly what the new one gains. The
records carrying the `P` flag — candidates LOPT is shown but does not fit — are
not components and stay at weight zero. This is
`make_LOPT_input.blend_weights()`, the same rule that a full rebuild applies.

### Excel-friendly output files

All validation tables are written by `save_table()` (in `chance_mc.py`): every CSV gets an `.xlsx` twin, and floating-point columns are rounded to physically meaningful decimals. The rounding matters for CSVs: Python prints a 64-bit float with up to 17 significant digits (the number needed to reproduce the binary value exactly — not extra precision), while Excel reads at most 15 and converts longer numbers to text; in the `.xlsx` twins, J values such as `3/2` stay text instead of being converted to dates. If a target file is locked (open in Excel), the writer falls back to a `_new`-suffixed name instead of aborting the run.

### Workflow

```
python classify_lines.py     # 1. the real classification → line_classifications.csv/.xlsx
python decoy_mc.py           # 2. eight decoy runs → decoy_mc_*.csv/.xlsx
python chance_mc.py          # 3. optional: shifted-wavenumber cross-check → chance_mc_*.csv/.xlsx
python level_shifts.py       # 4. calibrations, probabilities → level_shift_report.csv/.xlsx
python level_interchange.py  # 5. interchanged identities → level_interchange.csv/.xlsx
python level_positions.py --scan  # 6. alternate positions, question marks → level_positions.csv
python check_sync.py         # 7. do all the files still describe the same identification?
python unfound_levels.py     # 8. which levels nobody has found are worth searching for
python level_positions.py --unknown 742   # 9. and where the lines want one of them
python insert_new_level.py --iden2-row 742   # 10. and how one of them gets in
```

Repairing an interchange that step 5 flags is a separate act, done once and by hand:

```
python swap_line_assignments.py id1 id2        # LOPT files, IDEN2 files, lopt, pipeline overlay
python swap_line_assignments.py id1 id2 --after-lopt   # only if lopt did not run above
```

When the **observed intensities** have to be recalibrated — a new intensity column, or a
correction found to be wrong — three more steps come first, and they invalidate everything
downstream, so the four above must all be repeated afterwards:

```
python tools/calibrate_intensities.py --col-intensity Iorig --cover 821.92 10721.57        --out-functions intensity_correction_sugar.txt --plot intensity_calibration.png
python tools/apply_calibration.py --orig-col Iorig --scale 1000        --old intensity_correction_functions.txt --new intensity_correction_sugar.txt
cp intensity_correction_sugar.txt intensity_correction_functions.txt   # adopt it
python tools/iterate_boltzmann.py --tol 0.005     # re-converge C, kT; rewrites Icalc.xlsx
```

`tools/estimate_snr.py` stands outside this chain: it reads the line list and the
correction and writes only its own outputs, so it can be run whenever, but its numbers are
only as good as the calibration in force when it ran.

After identifications have been revised by hand and LOPT re-run on them, step 4 is
repeated against the revised run, reusing the calibrations of steps 2–3:

```
python make_LOPT_input.py                     # optional: pipeline table → LOPT input files
#   … revise the identifications, run LOPT …
python level_shifts.py --lopt LOPT_output_lines_revised.txt \
                       --e-input revised_level_energies.csv
```

The first three accept `--missing-gA {none,impute}`, which overrides `missing_gA.policy` of
the configuration file. Steps 1 and 2 must be run under the **same** policy, since step 2
calibrates step 1.

Step 2 reads the baseline output of step 1 (to define the tested levels and the perturbation reference), so the order matters. Each full validation run takes a few minutes. Step 5 needs only step 1 and `IDEN2/enlev.dat` — it uses no decoys — and takes a few seconds; it accepts the same `--lopt` / `--e-input` switches as step 4.

---

## Configuration

**`lineclass_config.toml`**, read by `config.py` (`config.load()`) at the start of every
script, holds everything that describes the *inputs*: file names, worksheet names, and the
column names the readers look for. Columns are found by their **name** in the header row of
the worksheet, not by position, so a rearranged input workbook needs no change in the code.
Relative paths are taken relative to the directory holding the configuration file.

| section                     | what it fixes                                                                     |
|-----------------------------|-----------------------------------------------------------------------------------|
| `[files]`                   | the input/output workbook names, plus four optional overlays: `level_overrides` (revised adopted energies), `line_decisions` (the decision ledger), `new_levels` (levels found since the level list was published) and `icalc_extra` (calculated transitions supplied by hand; commented out, since a new level's are now derived from `tp_E1_no_trials.xlsx`) |
| `[range]`                   | `wn_min`, `wn_max`: the wavenumber interval (cm⁻¹) in which candidate transitions are generated |
| `[levels.layout]`, `[lines.layout]`, `[icalc.layout]` | worksheet name and column names of each input file            |
| `[icalc.completeness]`      | `gA_cutoff`: the printing threshold of Cowan's codes, 1000 s⁻¹ — the basis of the censoring correction above |
| `[missing_gA]`              | `policy` (`"impute"` or `"none"`) and the settings of the imputation recipe          |
| `[intensity_model]`         | `C` and `kT` of `Icalc = C·gA·exp(−Eup/kT)/rwn`, plus the tolerance of the re-fit check |

The remaining tunables are decisions about the *method* rather than the data, and stay as
module constants and function defaults in `classify_lines.py`:

| Name                       | Value                | Meaning                                                        |
|----------------------------|----------------------|----------------------------------------------------------------|
| `WN_MIN`, `WN_MAX`         | `9327.0`, `121665.0` | Allowed Ritz-wavenumber range (cm⁻¹); set from `[range]` of the configuration file |
| `MISSING_POLICY`           | `'impute'`           | Treatment of a level pair absent from `Icalc.xlsx`; set from `[missing_gA]`, overridden by `--missing-gA` |
| `ATOMIC_MASS`              | `140.90765` u        | Pr mass for Doppler width                                      |
| `T`                        | `1.6` eV             | Plasma temperature for Doppler width                           |
| match tolerance            | `5.5 σ_comb`          | Candidate window in `match_and_grade`; `σ_comb` combines `wn_uncertainty` with `u_energy` of both levels |
| level-uncertainty sweeps   | `max_sweeps=50`, `tol=1e-4` | Convergence cap/tolerance for `compute_level_uncertainties`'s propagation sweep |
| `USE_INTENSITY_ADJUSTMENT` | `0`                  | Per-level intensity factors off (`1` enables them)             |
| `FACTOR_WEIGHTING`         | `'unweighted'`       | Factor averaging when enabled (`'mandel_paule'` optional)      |
| `alpha`                    | `0.5`                | Factor-update damping in `weed_assignments`                    |
| `DECISIONS_WN_MATCH`       | `0.01` cm⁻¹          | How far a decision-ledger wavenumber may miss its observed line |
| `min_n`                    | `10`                 | Min. qualifying transitions per level factor                   |
| `CLAMP`                    | `5.0`                | Max. magnitude of a per-level ln-factor                        |
| `max_iterations`           | `100` (weeding)      | Inner weeding iteration cap                                    |
| `max_cycles`               | `20` (`main`)        | Outer match→weed→optimize cycles                               |
| `DECOY_PREFIX`             | `'D_'`               | Level-id prefix of decoy levels (validation runs only)         |

---

## Project Structure

```
LineClass/
├── README.md                     ← this file
├── classify_lines.py             # Full pipeline: read → generate → match → grade → resolve → weed → optimize → output
│                                 #   (+ decoy-level support for the validation suite)
├── models.py                     # Dataclasses: EnergyLevel, SpectralLine, Transition, UNASSIGNED
├── decoy_mc.py                   # Validation: decoy (shadow-level) runs → false-confirmation rates in situ
├── level_shifts.py               # Validation: calibrations, criterion grids, pattern scores, p_spur; --detail mode
├── level_interchange.py          # Validation: are two levels of one parity and J wearing each
│                                 #   other's calculated intensities?  --detail mode
├── level_positions.py             # Validation: the likelihood of a level position, and the
│                                 #   scan that puts a question mark on it; --detail, --scan
├── swap_line_assignments.py      # Repair: run the three scripts below in order
├── swap_line_assignments_LOPT.py # Repair: exchange two levels' lines in the LOPT transitions
│                                 #   file, keeping both level identifiers where they are
├── swap_line_assignments_IDEN.py # Repair: exchange two levels' observed energies and lines in
│                                 #   enlev.dat and trans.dat, keeping IDEN2's numbering
├── swap_line_assignments_pipeline.py  # Repair: record the exchange in the level overrides and
│                                 #   the decision ledger, so a pipeline re-run keeps it
├── swap_paths.py                 # Where the four scripts above look for their files:
│                                 #   the current directory first, then the project
├── chance_mc.py                  # Shared utilities (read_input_levels, per_level_table, save_table,
│                                 #   apply_policy_option); run directly for the optional
│                                 #   shifted-wavenumber cross-check
├── config.py                     # Reads lineclass_config.toml; file, sheet and column names
├── lineclass_config.toml         # Configuration: inputs, layouts, gA cutoff, missing-gA policy
├── gA_imputation.py              # The censoring correction: gA_missing, u_ln_missing, Icalc model
│                                 #   (`--report` prints the calibration)
├── tests/                        # pytest suite (python -m pytest -q)
├── tools/                        # compare_runs.py, analyze_policy_diff.py, fit_boltzmann.py
│                                 #   (Boltzmann fit of C and kT), iterate_boltzmann.py
│                                 #   (fit → rewrite Icalc → reclassify, to convergence),
│                                 #   calibrate_intensities.py (standalone: finds the
│                                 #   piecewise plate-intensity correction and its regions),
│                                 #   apply_calibration.py (takes the intensity column
│                                 #   from one such correction to another),
│                                 #   attach_sugar_intensities.py (copies Sugar's printed
│                                 #   intensities into the line list as the column Iorig),
│                                 #   estimate_snr.py (standalone: the detection threshold
│                                 #   of the line list, and the IDEN2 files rewritten with
│                                 #   signal-to-noise ratios in place of intensities)
├── baseline/                     # Archived reference runs (policy_none/, policy_impute/,
│                                 #   boltzmann_model/, boltzmann_converged/,
│                                 #   before_intensity_recalibration/,
│                                 #   before_sugar_intensities/,
│                                 #   sugar_intensities_converged/,
│                                 #   sugar_intensities_complete/ — the run in force)
├── PLAN_missing_gA.md            # The censoring correction: plan, evidence, Stage-5 verdict
├── Icalc.xlsx                    # Input: calculated transition intensities & uncertainties
├── Pr3_lines.xlsx                # Input: observed spectral lines (Sugar 1969/1974);
│                                 #   column Iorig holds Sugar's own printed intensities
├── Pr3_Sugar69_extracted_*.xlsm  # Source: the checked extraction of Sugar 1969, Table 1
├── Pr_3_Sugar74_Table1_*.xlsm    # Source: the checked extraction of Sugar 1974, Table 1
├── intensity_correction_functions.txt  # Input (optional): Sugar plate-intensity calibration
├── new_levels.txt                # Input (optional): levels found since the level list
│                                 #   was published, tab-separated (level_id, E, J,
│                                 #   parity, iden2_row, cowan_lid, comment)
├── new_levels.csv                # The comma-separated form it replaced, kept as a backup
├── icalc_new.xlsx                # Retired: the calculated transitions of those levels,
│                                 #   now derived from tp_E1_no_trials.xlsx instead
├── insert_new_level.py           # Puts a newly found level into the pipeline, end to end
├── revised_level_energies.csv    # Input (optional): revised adopted energies
├── line_decisions.csv            # Input (optional): the decision ledger
├── IDEN2/                        # The IDEN2 working files as last saved (dlv.dat,
│                                 #   TRANS.DAT, enlev.dat, numset.dat)
├── IDEN2_snr/                    # The same, with signal-to-noise ratios in place of
│                                 #   intensities (tools/estimate_snr.py)
├── noise_level.txt               # Output: the fitted detection threshold ln N(lambda)
├── line_snr.csv                  # Output: per-line noise level and SNR
├── noise_windows.csv             # Output: the window minima the noise fit rests on
├── line_classifications.xlsx     # Output: classification table (LOPT-ready)
├── line_classifications.csv      # Output: same data as CSV
├── level_shift_report.csv/.xlsx  # Output: per-level validation table (dE, pattern scores, p_spur)
├── decoy_mc_*.csv/.xlsx          # Output: decoy-run tables (levels, real levels, lines, summary)
└── chance_mc_*.csv/.xlsx         # Output: shifted-wavenumber cross-check tables (optional)

# Shared / external to this directory:
../TableExtraction/Pr3_lev_Wyart_1999.xlsm   # Input: energy levels (sheet 'Wyart2000')
../statistics.py                             # Bundled mandel_paule() (used only if FACTOR_WEIGHTING='mandel_paule')
```

---

## Pipeline Architecture

```
┌─────────────────────┐    ┌────────────────────┐     ┌───────────────────────┐    ┌─────────────────────────┐
│  Read Levels,       │───→│  Generate & Match  │────→│  Resolve Conflicts    │───→│ Iterative Weeding       │
│  Lines, Icalc       │    │  (binary search,   │     │  (conservative sort,  │    │ per-level factors,      │
│  (Steps 1–4)        │    │   2D grading)      │     │   F/R flags)          │    │ 3-step per-line logic,  │
└─────────────────────┘    └─────────↑──────────┘     └───────────────────────┘    │ blacklist oscillations) │
                                     │                                             └──────────┬──────────────┘
                                     │                                                        │
                                     │                                         ┌──────────────▼────────────────┐
                                     │                                         │  Energy Level Optimization    │
                                     │                                         │  (weighted mean, BF weights)  │
                                     │                                         └──────────┬────────────────────┘
                                     │                                                    ▼
                                     │                                                   / \
                                     │                                                  /   \
                                     │                                                 /     \
                                     │                                                /       \
                                     │                                        No     /         \
                                     └ ──────────────────────────────────────────── / Converged?\
                                                                                    \           /
                                                                                     \         /
                                                                                      \       /
                                                                                       \     /
                                                                                        \   /
                                                                                         \ /
                                                                                          │ Yes
                                                                               ┌──────────▼──────────────┐           
                                                                               │  LOPT-ready output      │
                                                                               │  (xlsx + csv)           │
                                                                               └─────────────────────────┘
```

The outer loop repeats until both the number of accepted transitions and the level energies are stable across cycles (first-iteration max energy change `< 0.001 cm⁻¹`).

---

## License

For academic and research purposes, as part of the parent `TableExtraction` project.
