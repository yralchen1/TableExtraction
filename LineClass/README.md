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

The classification is accompanied by a **statistical validation suite** (`decoy_mc.py`, `level_shifts.py`, with shared utilities in `chance_mc.py`) that measures how often the pipeline would confirm a level that is *not* real, and assigns to every tested level a probability of being spurious. See **Validation of the classification** below.

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

## Output

Two files are written to the `LineClass/` directory (`OUTPUT_FILE`, `OUTPUT_CSV`):

| File                        | Description                                            |
|-----------------------------|--------------------------------------------------------|
| `line_classifications.xlsx` | Classification table with per-column number formatting |
| `line_classifications.csv`  | Identical data in CSV                                  |

Rows are sorted by **decreasing observed wavenumber** (`wn_obs`), then decreasing Ritz wavenumber (`rwn`), then increasing `grade` for ties (`build_output`). Unclassified observed lines still appear, as a single blank-classification row.

### Output Columns (22)

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
| 10 | `intens_from_f`    | Upper-level intensity correction factor (ln scale)                         |
| 11 | `intens_to_f`      | Lower-level intensity correction factor (ln scale)                         |
| 12 | `dif_wn_O-C`       | Observed − Ritz wavenumber (cm⁻¹)                                          |
| 13 | `grade`            | 2D grade: tier (`2`–`5`) + subgrade (`A`–`E`, `G`); see *Grading*          |
| 14 | `notes1`           | Conflict flags: `F` (conflicting), `R` (revised)                           |
| 15 | `notes2`           | Per-transition decision trace from weeding (`Step1`/`Step2`/`Step3` label) |
| 16 | `new`              | `1` = new classification, `0` = original Sugar, blank = unclassified       |
| 17 | `accepted`         | `1` = accepted by weeding, `0` = rejected, blank = unclassified line       |
| 18 | `n_accepted`       | Number of accepted classifications of this observed line (0, 1, 2, …)      |
| 19 | `low_E`            | Lower level energy (cm⁻¹)                                                  |
| 20 | `upp_E`            | Upper level energy (cm⁻¹)                                                  |
| 21 | `rwn`              | Ritz wavenumber = `upp_E − low_E` (cm⁻¹)                                   |
| 22 | `BF`               | Branching fraction of this component (0 if not accepted)                   |

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
(is_new, |Obs − Ritz|, tier, z, int_err)
```

so legacy assignments are favored, then the smallest wavenumber mismatch, then tier, then z (σ residual), then intensity mismatch. The winner is flagged **`F`** (conflicting); if any *original* (legacy) assignment lost, the winner is also flagged **`R`** (revised). Losing transitions are removed from their lines; a line that loses its only classification (and had originals) receives the `UNASSIGNED` sentinel.

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
4. **Oscillation blacklist:** the blacklist is created once in `main()` and passed into `weed_assignments`, so it **persists across the outer cycles**. A transition enters it in two ways: (a) it flips acceptance state between consecutive weeding passes within a cycle; (b) its end-of-cycle acceptance state alternates between outer cycles (A→B→A over the last three cycle snapshots, detected by `blacklist_cycle_oscillations` in `main()`). Blacklisted transitions are permanently **rejected** in Step 3 ("decisions oscillate between weeding iterations or cycles") — an oscillating decision is treated as evidence of an unreliable assignment. Blacklist keys are `trans_key()` tuples (lower id, upper id, identity of the assigned line), which remain stable when transition objects are recreated between cycles.

### Per-line 3-step decision (`weed_assignments_line`)
Candidates are sorted by decreasing `calc_intensity`, then by Ritz σ.

- **Step 1 — clear-cut decisions.** New transitions with Ritz mismatch `> 3 σ_comb` are rejected outright, and (for legacy candidates with no theoretical intensity) mismatch `> 5 σ_comb` is likewise an outright rejection — `σ_comb` is the same combined line+level-uncertainty sigma used for matching (`_ritz_sigma`, see *Match & grade* above and *Level-Energy Uncertainty Estimation* below). A **relative-intensity filter** (`_apply_relative_intensity_filter`) rejects new candidates contributing `< 10 %` of the accepted cumulative intensity, and applies conservative tests to legacy candidates around the `4 %` / `2 %` thresholds (`check_internal_significance` guarantees rejection of statistically-insignificant legacy lines). Surviving candidates are tested by **asymmetric z-score** against `I_obs` — thresholds are looser when the prediction is *too strong* and for **resonance** lower levels (ids ending `001`, via `_fudge_factor_for_asym_intensity`), stricter when *too weak*. Decisions are tagged `Step1…` in `notes2`.
- **Step 2 — grouping & late arbitration.** Undecided candidates are tested for **pair** and **triple** acceptance (`_try_pair_stage`, `_try_triple_stage`): a group passes if its intensity-weighted **center of gravity** matches the observed wavenumber and its **effective spread** — the Doppler width (full width at half maximum, `sqrt(8·ln2·kT/m)/c` = 8.22×10⁻⁶ of the wavenumber for Pr at `T = 1.6 eV`, `ATOMIC_MASS = 140.90765 u`) convolved with the measurement σ (`_effective_spread_combined`) — is consistent with the line profile. For `h`/`w` (and other broadened) characters the thresholds are relaxed by a broadening multiplier (`_broadening_multiplier`). A CoG-sensitivity test rejects candidates that *worsen* the group CoG by more than 1 σ. Finally, still-undecided **new** candidates carrying `F`/`R` in `notes1` are rejected. Decisions are tagged `Step2…`.
- **Step 3 — defaults.** Any remaining undecided **new** candidate is rejected ("no solid evidence for acceptance"); any remaining **legacy** candidate is accepted ("no solid evidence for rejection"). Blacklisted transitions land here and are annotated "decisions oscillate in iterations". Decisions are tagged `Step3…`.

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
`DEFAULT_GAPS_A`, the same numbers as `COVERAGE_GAPS_A` in `level_shifts.py`) and cut the
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

---

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

## Validation of the classification (`decoy_mc.py`, `level_shifts.py`)

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
3. **Intensity-pattern check** — the modern form of the classical "square-array" argument: a real level must reproduce the *pattern* of theoretically predicted intensities. For every level, all predicted transitions to/from it (from `Icalc.xlsx`, both partners known, Ritz wavenumber inside the observed range) are sorted by predicted intensity, and three scores are computed. Two classes of predictions are first excluded, because a transition that could not have been seen must not count as missing: (a) predictions whose Ritz wavelength falls inside one of **Sugar's coverage gaps** (`COVERAGE_GAPS_A`: 1522.49–1529.85 and 2103.46–2107.92 Å vacuum — intervals with no exposures, whose edges coincide with the seams of the intensity-calibration regions); (b) when the intensity-calibration file `intensity_correction_functions.txt` is present (piecewise polynomials `P(λ)` in vacuum wavelength converting Sugar's plate intensity to the uniform linear scale via `I_linear = 1000·I_Sugar·exp(P(λ))`), predictions whose Sugar-scale intensity `I_pred/(1000·exp(P(λ)))` falls below the noise level (Sugar intensity 1). A third correction enters that comparison. `Icalc` is meant to be on the same linear scale as the observed intensities, and from 1000 Å upward it is — the median of `I_obs/I_pred` over the accepted identifications is within about 30 % of 1. Below 900 Å it is not: the median falls to 0.09 (83 accepted lines), and the ratio of the band totals — all observed intensities over all predicted ones, which does not depend on which lines were classified — agrees at 0.16. A prediction there is therefore about ten times fainter on the plate than its `I_pred` says. `intensity_scale_bias` measures that factor band by band (`SCALE_BIAS_BANDS_A`, at least `SCALE_BIAS_MIN_N` = 10 accepted lines, correction applied only where the factor is below `SCALE_BIAS_MAX_FACTOR` = 0.5) together with the scatter of `ln(I_obs/I_pred)`, and the corrected intensity `I_pred·f(λ)` is what has to clear the noise level. This is a fault of the calculated intensities, not of the plates: Sugar's plates reach his intensity 1 everywhere, and Pr III lines of intensity 1 are recorded and classified down to the short-wavelength edge at 821.9 Å. With the current inputs 113 predictions are excluded in the gaps and 21178 below the noise level, of 29162 loaded — most theoretically predicted transitions of these high levels are simply too faint for Sugar's plates. Accepted lines whose predictions are below the noise level stay in the support count n (the decoys are treated identically, so the calibration stays fair); the report column `n_acc_below_noise` makes them visible.
   - `top10_found` — how many of the 10 strongest predictions are among the accepted lines (`n_top10` = how many were available);
   - `pattern_C` — intensity-weighted completeness: the sum of predicted intensities over the accepted transitions divided by the sum over all predicted observable ones;
   - `pattern_V` — C divided by the largest C achievable with the level's number of matched lines. V = 1 means the accepted lines are exactly the strongest predictions; V = 0 means predictions exist but none was found. Unlike C, V does not punish a level for accepted lines that theory does not cover, and it is nearly independent of n.
   - `top1` — the fate of the level's **strongest** predicted observable transition on its own: `found` (it is among the accepted lines), `masked` (absent, but an observed line hides it — see the masking test in item 5), `faint` (absent, nothing hides it, but the line would not reliably have been recorded anyway — see below), `missing` (absent and unexplained), or blank (no observable prediction). This is scored separately because `pattern_V` cannot catch it: V asks whether the *accepted* lines are the strongest of the predictions, so a level whose two accepted lines are ranked 2 and 3 reaches V ≈ 1 with rank 1 missing. A missing, unmasked strongest branch is the most direct evidence there is against a level — but only where the line would have been seen. `detection_probability` turns the corrected intensity and the measured scatter of `ln(I_obs/I_pred)` into the probability that the branch would have been recorded at all (a log-normal centred on `I_pred·f(λ)` with that scatter, integrated above the noise level); below `DETECT_CONFIDENCE` = 0.9 the absence carries no information and `top1` reads `faint` instead of `missing`. The two classes separate cleanly — with the current inputs every branch still called `missing` would have been recorded with probability ≥ 0.94, and every `faint` one with 0.85 or less — so the exact cut does not matter. Of the 208 tested levels, 180 have their strongest branch accepted, 10 have it masked, 2 have it faint and 16 have it missing.
   
   Calibration: genuine (old) levels have a median C of 0.91 and median V of 0.94; decoys have 0.03 and 0.06. The information here is genuinely new — the weeding checks each line's intensity individually but never penalizes a level for strong predicted lines that are *absent*. Two known weaknesses: a single grossly wrong predicted rate can wreck C and V of a real level (theory faults of this kind are common), and levels whose transitions are not covered by the theory get no pattern verdict at all.
4. **Per-level spurious probabilities.** For each tested level, two independent pieces of evidence are weighed between "genuine" and "spurious" (= behaves like a decoy): the energy shift (|d| modeled as a log-normal distribution whose median depends on n, fitted separately to the old levels and to the decoys) and the intensity pattern (the probability of the level's V range under each hypothesis, measured on the old levels and on the decoys). The prior probability that a tested level with support n is spurious is `S·q(n)/m(n)`, where `q(n)` is the decoy probability of ending with support n and `m(n)` is the observed number of tested levels with support n; the single unknown S — the total number of spurious levels — is estimated by maximizing the likelihood. Three posterior probabilities go into the report:
   - `p_spur_decoys` — from the energy shift alone;
   - `p_spur_pattern` — from the intensity pattern alone;
   - `p_spur` — both folded together (the final value).
   
   With the current inputs (converged model, 208 tested levels): S = 0.0 [0.0–8.0 at 95% confidence] from the energy shifts alone (36.8 [19.5–53.8] under the alternative small-n assumption — poorly determined), 4.8 [0.0–13.0] from the pattern alone, and **3.5 [0.0–10.8] folded** (9.0 in the robustness variant) — the pattern evidence removes the model sensitivity. The sum of `p_spur` over the tested levels reproduces S, and the expected number of false confirmations under any cut equals the sum of `p_spur` over the levels passing it. Marking levels with `p_spur ≥ 0.1` as questionable flags 8 levels; the 200 unmarked levels then carry a summed probability of only 0.7 expected spurious levels, while the marked group is expected to contain ~2.9 spurious and ~5.2 genuine members (a "questionable" mark is a caution, not a verdict).
5. **Masking check (adjudication of the questionable marks).** Two things put a tested level up for examination: the probabilities question it (`p_spur ≥ 0.1`), or its strongest predicted transition is absent from the accepted lines and would have been recorded had it been there (`top1` reads `missing`; a branch too faint to have been recorded reads `faint` and does not put the level up for examination). A strong predicted transition may be absent from the accepted lines simply because a nearby **stronger** observed line hides it — on the photographic plates, a weaker line close to a much stronger one cannot be measured. Every such level is therefore re-examined. The test uses the **effective line width** (full width at half maximum): the instrumental width `INSTR_FWHM_A` (0.035 Å, estimated from the closest measured line pairs; the spectrograph resolution is nearly constant in wavelength) convolved with the Doppler width, which grows with wavelength. For each missing prediction among the level's strongest ten, an observed line hides it in either of two regimes: **within one effective width** the two lines are not resolved, so a line of at least `MASK_STRENGTH` (0.5) times the predicted intensity absorbs it; **beyond one width**, a Gaussian profile of the effective width is applied to the observed line, and the prediction counts as hidden if the profile intensity at its position still exceeds the predicted intensity — so the masking reach of a much stronger line extends into its wings, growing with the intensity ratio. (The manually verified masking cases show separations up to 0.029 Å, inside one width.) If at least `MASK_CLEAR_FRACTION` (90%) of the *missing predicted intensity* is hidden — i.e. the dominant missing predictions are masked, while a small residue of weak ones is tolerated (the theory is least reliable for weak transitions) — the bad pattern score carries no evidence, and the questionable mark is cleared — provided the energy-shift evidence alone is not suspect (`p_spur_decoys < 0.1`) **and** the strongest branch is not itself among the unexplained absences. A missing, unmasked strongest branch keeps the mark whatever the probabilities say. The outcome goes into two report columns: `question_status` (`removed` = mark cleared, `retained` = mark kept) and `reason` (e.g. "strong missing transitions are masked by nearby stronger lines", "no reason to clear", "no theoretical predictions (pattern check not applicable)", "masking found, but the energy shift alone remains suspect"). With the current inputs, 20 levels are examined (8 flagged by `p_spur`, 16 by a missing strongest branch, 4 by both); 1 is cleared (059003.000457, 96% of its missing predicted intensity masked), leaving 19 questionable — 16 because their strongest predicted transition is absent and nothing hides it, 3 with no reason to clear.
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

### Limits of the validation (to be stated alongside the results)

- **Recovery, not physical proof.** A small ΔE certifies that the accepted lines reproduce the energy encoded in Wyart's input value — i.e. that his identifications were recovered. If Wyart himself was misled by chance coincidences, our run re-finds the same coincidences with a small ΔE; only the intensity pattern (and physics arguments: theory, g-factors, term structure) can catch that case.
- **Theory-fault sensitivity.** A bad pattern score can mean a spurious level *or* a level whose theoretical description is wrong; levels whose accepted lines are mostly uncovered by theory (intensity grade `G`) get weak-quality pattern verdicts. In the converged run all 5 marked levels do have observable predictions (in earlier runs a sizeable part of them did not, and there `p_spur` rests on the energy shift and the support count alone). A related caution: for a few marked levels every observable prediction lies within a factor ~3 of the local noise — there a bad pattern score means little, because the predicted intensities themselves scatter by a factor ~2.5 against the observed ones.
- **Decoys can accidentally be real.** A displaced decoy may land on a real, previously unknown level (most likely a neighboring-J member of the same term at high energies) or on the true position of a level misassigned in the underlying Cowan-code fit. Both effects make some decoy "false positives" actually real, so the decoy-based false rates are slight overestimates — the bias is in the safe direction. The fraction of supported decoy trials showing three or more of their top-10 predictions accepted (~6%) is an upper bound on this contamination.
- **Lines omitted as blends with other ionization stages.** Sugar assigned lines to Pr II, III, or IV by comparing exposures at different degrees of excitation; a Pr III line nearly coinciding with a stronger Pr II or Pr IV line could not be assigned confidently and was omitted from his Pr III list. Such omissions are invisible to the masking check (the search covers only Sugar's Pr III lines), so some "missing" strong predictions are excusable in a way the automation cannot see. Statistically the effect is absorbed — the old (genuine) reference levels suffer the same omissions, so the pattern-score comparison between the classes stays fair — but the per-level adjudication is conservative: a questionable mark that an expert would clear on this ground stays retained. Automating this excuse would require Pr II and Pr IV line lists.

### Excel-friendly output files

All validation tables are written by `save_table()` (in `chance_mc.py`): every CSV gets an `.xlsx` twin, and floating-point columns are rounded to physically meaningful decimals. The rounding matters for CSVs: Python prints a 64-bit float with up to 17 significant digits (the number needed to reproduce the binary value exactly — not extra precision), while Excel reads at most 15 and converts longer numbers to text; in the `.xlsx` twins, J values such as `3/2` stay text instead of being converted to dates. If a target file is locked (open in Excel), the writer falls back to a `_new`-suffixed name instead of aborting the run.

### Workflow

```
python classify_lines.py     # 1. the real classification → line_classifications.csv/.xlsx
python decoy_mc.py           # 2. eight decoy runs → decoy_mc_*.csv/.xlsx
python chance_mc.py          # 3. optional: shifted-wavenumber cross-check → chance_mc_*.csv/.xlsx
python level_shifts.py       # 4. calibrations, probabilities → level_shift_report.csv/.xlsx
```

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

Step 2 reads the baseline output of step 1 (to define the tested levels and the perturbation reference), so the order matters. Each full validation run takes a few minutes.

---

## Configuration

**`lineclass_config.toml`**, read by `config.py` (`config.load()`) at the start of every
script, holds everything that describes the *inputs*: file names, worksheet names, and the
column names the readers look for. Columns are found by their **name** in the header row of
the worksheet, not by position, so a rearranged input workbook needs no change in the code.
Relative paths are taken relative to the directory holding the configuration file.

| section                     | what it fixes                                                                     |
|-----------------------------|-----------------------------------------------------------------------------------|
| `[files]`                   | the four input/output workbook names                                                |
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
│                                 #   intensities into the line list as the column Iorig)
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
