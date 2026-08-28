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

**Consequence for convergence.** Because `u_energy` is recomputed from scratch every cycle from that cycle's starting (previous-cycle) acceptance state, and the matching/Step-1 tolerances depend on it, a transition that is rejected in one cycle purely because a partner level's uncertainty was still underestimated can be picked up in a later cycle, once more of that level's other transitions have been accepted and its `u_energy` has grown to reflect the real scatter.

---

## Energy Level Optimization (`optimize_levels`)

After weeding, each non-ground level (`energy != 0`) is re-estimated as the **weighted mean** of the energies implied by its accepted transitions:

- from a transition *from* the level (it is the upper level): `E = E_lower + wn_obs`,
- from a transition *to* the level (it is the lower level): `E = E_upper − wn_obs`.

Weights come from `calc_weights`: for each accepted transition on a line, `BF = I_calc / Σ I_calc` over the line's accepted set (or `1/n` when any intensity is missing), and `w = BF² / σ_wn²` — the LOPT "centroid" weighting that shares a blend's weight among its components. Lines with a single accepted transition get `BF = 1` → `w = 1/σ²`. Iteration runs up to `MAX_IT = 50`, stopping when the maximum energy change `< TOL = 0.001 cm⁻¹`. The function returns the *first-iteration* maximum change, which the outer cycle uses as its convergence signal.

> This optimization is intentionally minimalistic — it exists only to sharpen the Ritz wavenumbers that drive assignment decisions. **Final** level optimization is deferred to **LOPT (v ≥ 5)**, fed by the accepted rows and their `BF` values.

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

3. The intensity follows from the relation the file itself obeys, to eight significant
   figures over all 30206 rows:

   ```
   Icalc = C · gA · exp(−Eup / kT) / rwn ,      C = 135.8 ,  kT = 12905 cm⁻¹
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
3. **Intensity-pattern check** — the modern form of the classical "square-array" argument: a real level must reproduce the *pattern* of theoretically predicted intensities. For every level, all predicted transitions to/from it (from `Icalc.xlsx`, both partners known, Ritz wavenumber inside the observed range) are sorted by predicted intensity, and three scores are computed. Two classes of predictions are first excluded, because a transition that could not have been seen must not count as missing: (a) predictions whose Ritz wavelength falls inside one of **Sugar's coverage gaps** (`COVERAGE_GAPS_A`: 1522.49–1529.85 and 2103.46–2107.92 Å vacuum — intervals with no exposures, whose edges coincide with the seams of the intensity-calibration regions); (b) when the intensity-calibration file `intensity_correction_functions.txt` is present (piecewise polynomials `P(λ)` in vacuum wavelength converting Sugar's plate intensity to the uniform linear scale via `I_linear = 1000·I_Sugar·exp(P(λ))`), predictions whose Sugar-scale intensity `I_pred/(1000·exp(P(λ)))` falls below the noise level (Sugar intensity 1). With the current inputs 56 predictions are excluded in the gaps and 11905 below the noise level, of 18331 loaded — most theoretically predicted transitions of these high levels are simply too faint for Sugar's plates. Accepted lines whose predictions are below the noise level stay in the support count n (the decoys are treated identically, so the calibration stays fair); the report column `n_acc_below_noise` makes them visible.
   - `top10_found` — how many of the 10 strongest predictions are among the accepted lines (`n_top10` = how many were available);
   - `pattern_C` — intensity-weighted completeness: the sum of predicted intensities over the accepted transitions divided by the sum over all predicted observable ones;
   - `pattern_V` — C divided by the largest C achievable with the level's number of matched lines. V = 1 means the accepted lines are exactly the strongest predictions; V = 0 means predictions exist but none was found. Unlike C, V does not punish a level for accepted lines that theory does not cover, and it is nearly independent of n.
   
   Calibration: genuine (old) levels have a median C of 0.92 and median V of 0.95; decoys have 0.03 and 0.06. The information here is genuinely new — the weeding checks each line's intensity individually but never penalizes a level for strong predicted lines that are *absent*. Two known weaknesses: a single grossly wrong predicted rate can wreck C and V of a real level (theory faults of this kind are common), and levels whose transitions are not covered by the theory get no pattern verdict at all.
4. **Per-level spurious probabilities.** For each tested level, two independent pieces of evidence are weighed between "genuine" and "spurious" (= behaves like a decoy): the energy shift (|d| modeled as a log-normal distribution whose median depends on n, fitted separately to the old levels and to the decoys) and the intensity pattern (the probability of the level's V range under each hypothesis, measured on the old levels and on the decoys). The prior probability that a tested level with support n is spurious is `S·q(n)/m(n)`, where `q(n)` is the decoy probability of ending with support n and `m(n)` is the observed number of tested levels with support n; the single unknown S — the total number of spurious levels — is estimated by maximizing the likelihood. Three posterior probabilities go into the report:
   - `p_spur_decoys` — from the energy shift alone;
   - `p_spur_pattern` — from the intensity pattern alone;
   - `p_spur` — both folded together (the final value).
   
   With the current inputs: S = 18 [0–43 at 95% confidence] from the energy shifts alone (61 under the alternative small-n assumption — poorly determined), 37 [24–51] from the pattern alone, and **33 [22–46] folded** (42 in the robustness variant) — the pattern evidence removes the model sensitivity. The sum of `p_spur` over the tested levels reproduces S, and the expected number of false confirmations under any cut equals the sum of `p_spur` over the levels passing it. Marking levels with `p_spur ≥ 0.1` as questionable flags 56 levels; the 158 unmarked levels then carry a summed probability of only 3.7 expected spurious levels, while the marked group is expected to contain ~30 spurious and ~26 genuine members (a "questionable" mark is a caution, not a verdict).
5. **Masking check (adjudication of the questionable marks).** A strong predicted transition may be absent from the accepted lines simply because a nearby **stronger** observed line hides it — on the photographic plates, a weaker line close to a much stronger one cannot be measured. Every questionable level is therefore re-examined. The test uses the **effective line width** (full width at half maximum): the instrumental width `INSTR_FWHM_A` (0.035 Å, estimated from the closest measured line pairs; the spectrograph resolution is nearly constant in wavelength) convolved with the Doppler width, which grows with wavelength. For each missing prediction among the level's strongest ten, an observed line hides it in either of two regimes: **within one effective width** the two lines are not resolved, so a line of at least `MASK_STRENGTH` (0.5) times the predicted intensity absorbs it; **beyond one width**, a Gaussian profile of the effective width is applied to the observed line, and the prediction counts as hidden if the profile intensity at its position still exceeds the predicted intensity — so the masking reach of a much stronger line extends into its wings, growing with the intensity ratio. (The manually verified masking cases show separations up to 0.029 Å, inside one width.) If at least `MASK_CLEAR_FRACTION` (90%) of the *missing predicted intensity* is hidden — i.e. the dominant missing predictions are masked, while a small residue of weak ones is tolerated (the theory is least reliable for weak transitions) — the bad pattern score carries no evidence, and the questionable mark is cleared, provided the energy-shift evidence alone is not suspect (`p_spur_decoys < 0.1`). The outcome goes into two report columns: `question_status` (`removed` = mark cleared, `retained` = mark kept) and `reason` (e.g. "strong missing transitions are masked by nearby stronger lines", "no reason to clear", "no theoretical predictions (pattern check not applicable)", "masking found, but the energy shift alone remains suspect"). With the current inputs, 6 of the 56 flagged levels are cleared, leaving 50 questionable (28 with no reason to clear, 16 with no observable predictions, 6 where masking was found but the energy shift alone remains suspect).
6. **Report file** — `level_shift_report.csv/.xlsx`, one row per level: `level_id`, `is_new_level` (absent from ASD), `new_star` (tested), `note`, `J`, `parity`, `E_input`, `E_final`, `dE`, `n_old`, `n_new`, `n_tot`, `d`, `n_pred`, `n_top10`, `top10_found`, `pattern_C`, `pattern_V`, `n_acc_below_noise`, `p_spur_decoys`, `p_spur_pattern`, `p_spur`, `question_status`, `reason` (the probabilities and the adjudication columns are filled only where applicable).

**Single-level inspection:** `python level_shifts.py --detail <level_id> …` prints one level's predicted transitions (strongest first) with the fate of each in the real run — accepted (with the observed line), rejected (with the weeding note), or not matched by any observed line — plus the accepted lines that have no theoretical prediction. This is the working tool for judging individual questionable levels.

### Limits of the validation (to be stated alongside the results)

- **Recovery, not physical proof.** A small ΔE certifies that the accepted lines reproduce the energy encoded in Wyart's input value — i.e. that his identifications were recovered. If Wyart himself was misled by chance coincidences, our run re-finds the same coincidences with a small ΔE; only the intensity pattern (and physics arguments: theory, g-factors, term structure) can catch that case.
- **Theory-fault sensitivity.** A bad pattern score can mean a spurious level *or* a level whose theoretical description is wrong; levels whose accepted lines are mostly uncovered by theory (intensity grade `G`) get weak-quality pattern verdicts. Sixteen of the 56 marked levels have no observable predictions at all; their `p_spur` rests on the energy shift and the support count alone. A related caution: for a few marked levels every observable prediction lies within a factor ~3 of the local noise — there a bad pattern score means little, because the predicted intensities themselves scatter by a factor ~2.5 against the observed ones.
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
├── tools/                        # compare_runs.py, analyze_policy_diff.py
├── baseline/                     # Archived reference runs (policy_none/, policy_impute/)
├── PLAN_missing_gA.md            # The censoring correction: plan, evidence, Stage-5 verdict
├── Icalc.xlsx                    # Input: calculated transition intensities & uncertainties
├── Pr3_lines.xlsx                # Input: observed spectral lines (Sugar 1969/1974)
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
