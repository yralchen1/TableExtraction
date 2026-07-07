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

- **Bundled `statistics.py`** — the repository root contains a small `statistics.py` that provides `mandel_paule()` (Paule–Mandel consensus weighted mean, NBS 1982). `weed_assignments` imports it via `from statistics import mandel_paule`. This local module **shadows the standard-library `statistics` module**, so the repository root must be importable at run time (see *Running* below). The standard library's `statistics` does **not** provide `mandel_paule`.

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

There are no command-line flags. Behavior is controlled by module-level constants and function defaults (see *Configuration*). `main(max_cycles=10)` runs the full pipeline and writes the output files.

---

## Input Files

All paths are resolved relative to the script's own directory (`SCRIPT_DIR`).

| File                      | Location (constant)                   | Sheet       | Description                                                               |
|---------------------------|---------------------------------------|-------------|---------------------------------------------------------------------------|
| `Pr3_lev_Wyart_1999.xlsm` | `../TableExtraction/` (`LEVELS_FILE`) | `Wyart2000` | Energy levels (Wyart 1999 adopted values). Shared with `TableExtraction`. |
| `Icalc.xlsx`              | `LineClass/` (`ICALC_FILE`)           | `Sheet1`    | Calculated (theoretical) transition intensities and their uncertainties.  |
| `Pr3_lines.xlsx`          | `LineClass/` (`LINES_FILE`)           | `Sheet1`    | Observed spectral lines with existing (Sugar 1969/1974) classifications.  |

### Column mapping (as read by the code)

**Energy levels** — `read_energy_levels()`, sheet `Wyart2000` (header row skipped):

| Column    | #  | Field             | Notes                                 |
|-----------|----|-------------------|---------------------------------------|
| `J_adopt` | 11 | `J_str` → `J_val` | Parsed by `parse_J` (`'5/2'` → `2.5`) |
| `E_adopt` | 13 | `energy`          | cm⁻¹                                  |
| `par_ASD` | 14 | `parity`          | `e` / `o`                             |
| `ASD_id`  | 23 | `level_id`        | Unique level identifier (string)      |

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

### Output Columns (21)

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
| 18 | `low_E`            | Lower level energy (cm⁻¹)                                                  |
| 19 | `upp_E`            | Upper level energy (cm⁻¹)                                                  |
| 20 | `rwn`              | Ritz wavenumber = `upp_E − low_E` (cm⁻¹)                                   |
| 21 | `BF`               | Branching fraction of this component (0 if not accepted)                   |

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
             │   match_and_grade ─▶ resolve_conflicts ─▶       │
             │   weed_assignments (iterative) ─▶ calc_weights  │
             │             ─▶ optimize_levels                  │
             └───────────────────┬────────────────────────────┘
                                 │  converged?
                                 ▼
                 build_output ─▶ write_output (xlsx + csv)
```

The **outer cycle** in `main()` repeats *match → resolve → weed → weight → optimize levels* until the number of accepted transitions is unchanged **and** the largest level-energy change in the first optimization iteration is `< 0.001 cm⁻¹` (`max_cycles = 10`).

---

## The 6-Step Pipeline

The docstring and `main()` organize the work into six steps.

### Steps 1–3 — Read inputs
`read_energy_levels`, `read_transitions`, `read_observed_lines` load the three files into the dataclasses in `models.py`:

- **`EnergyLevel`** — id, energy, parity, J (string and float), plus per-level intensity factors and the lists `from_transitions` / `to_transitions`.
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
`assignment_cycle` runs three sub-steps:

**5.1 Match & grade** (`match_and_grade`, `assign_grades`).
For each observed line, a **bisection** over the sorted transition wavenumbers selects candidates within a tolerance of **`5.5 σ`** (`5.5 × wn_uncertainty`). Each candidate is reused if already attached to the line, else copied and attached. `assign_grades` sets the 2D grade and the `new` flag (`0` if the pair matches an `original_assignments` entry, else `1`).

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

**Tier — wavenumber agreement** (`wn_diff_abs = |wn_obs − wn_calc|`):

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
2. **Iterations 1…N** (`compute_intensity_factors`, `apply_intensity_adjustments`):
   - For every level compute a **per-level intensity factor** from its accepted, *unblended* transitions (`_compute_factor_for_transition_list`): qualifying transitions need `accepted == 1`, positive `orig_calc_intensity`/`orig_u_calc`, and a parent line with exactly **one** accepted transition. With **≥ `min_n` = 5** qualifiers, the factor is the **Mandel–Paule** weighted mean of `ln(I_obs / I_calc_orig)` (`_wm_mandel_paule`; a chi²-inflated weighted mean `_wm_red_chi` is the alternative, currently disabled). Factors are **clamped to ±5** ln-units.
   - Two factors per level are maintained: `intens_from_factor` (transitions *from* this level, i.e. it is the upper level) and `intens_to_factor` (transitions *to* this level, i.e. it is the lower level). See the docstring of `compute_intensity_factors` for the physical justification (population deviations, self-absorption, systematic errors in transition probabilities).
   - **Damped update:** `applied = alpha × raw + (1 − alpha) × prev` with **`alpha = 0.5`** (no damping on the first factor computation).
   - **Apply** to every transition: `I_adj = I_orig × exp(f_from + f_to)`; uncertainty adds in quadrature `u_adj = sqrt(u_orig² + u_from² + u_to²)`. Originals are never mutated, preventing feedback amplification.
   - **Reset** all decisions (`reset_weeding_state`) and re-weed.
3. **Convergence:** stop when no transition changes acceptance state between passes, else after `max_iterations` (`assignment_cycle` passes `max_iterations = 100`; the function's own default is 59).
4. **Oscillation blacklist:** after iteration 2, any transition that flips acceptance state between consecutive passes is added to a `blacklist` and thereafter skipped in Steps 1–2, falling through to its Step-3 default (`detect_oscillations` provides diagnostics; `snapshot_accepted`/`count_changes` track state).

### Per-line 3-step decision (`weed_assignments_line`)
Candidates are sorted by decreasing `calc_intensity`, then by Ritz σ.

- **Step 1 — clear-cut decisions.** New transitions with Ritz mismatch `> 3 σ` are rejected outright. A **relative-intensity filter** (`_apply_relative_intensity_filter`) rejects new candidates contributing `< 10 %` of the accepted cumulative intensity, and applies conservative tests to legacy candidates around the `4 %` / `2 %` thresholds (`check_internal_significance` guarantees rejection of statistically-insignificant legacy lines). Surviving candidates are tested by **asymmetric z-score** against `I_obs` — thresholds are looser when the prediction is *too strong* and for **resonance** lower levels (ids ending `001`, via `_fudge_factor_for_asym_intensity`), stricter when *too weak*. Decisions are tagged `Step1…` in `notes2`.
- **Step 2 — grouping & late arbitration.** Undecided candidates are tested for **pair** and **triple** acceptance (`_try_pair_stage`, `_try_triple_stage`): a group passes if its intensity-weighted **center of gravity** matches the observed wavenumber and its **effective spread** — the Doppler width (from `T = 1.6 eV`, `ATOMIC_MASS = 140.90765 u`) convolved with the measurement σ (`_effective_spread_combined`) — is consistent with the line profile. For `h`/`w` (and other broadened) characters the thresholds are relaxed by a broadening multiplier (`_broadening_multiplier`). A CoG-sensitivity test rejects candidates that *worsen* the group CoG by more than 1 σ. Finally, still-undecided **new** candidates carrying `F`/`R` in `notes1` are rejected. Decisions are tagged `Step2…`.
- **Step 3 — defaults.** Any remaining undecided **new** candidate is rejected ("no solid evidence for acceptance"); any remaining **legacy** candidate is accepted ("no solid evidence for rejection"). Blacklisted transitions land here and are annotated "decisions oscillate in iterations". Decisions are tagged `Step3…`.

---

## Energy Level Optimization (`optimize_levels`)

After weeding, each non-ground level (`energy != 0`) is re-estimated as the **weighted mean** of the energies implied by its accepted transitions:

- from a transition *from* the level (it is the upper level): `E = E_lower + wn_obs`,
- from a transition *to* the level (it is the lower level): `E = E_upper − wn_obs`.

Weights come from `calc_weights`: for each accepted transition on a line, `BF = I_calc / Σ I_calc` over the line's accepted set (or `1/n` when any intensity is missing), and `w = BF² / σ_wn²` — the LOPT "centroid" weighting that shares a blend's weight among its components. Lines with a single accepted transition get `BF = 1` → `w = 1/σ²`. Iteration runs up to `MAX_IT = 50`, stopping when the maximum energy change `< TOL = 0.001 cm⁻¹`. The function returns the *first-iteration* maximum change, which the outer cycle uses as its convergence signal.

> This optimization is intentionally minimalistic — it exists only to sharpen the Ritz wavenumbers that drive assignment decisions. **Final** level optimization is deferred to **LOPT (v ≥ 5)**, fed by the accepted rows and their `BF` values.

---

## Configuration

There is no config file; tunables are module constants and function defaults in `classify_lines.py`:

| Name               | Value                | Meaning                                      |
|--------------------|----------------------|----------------------------------------------|
| `WN_MIN`, `WN_MAX` | `9327.0`, `121665.0` | Allowed Ritz-wavenumber range (cm⁻¹)         |
| `ATOMIC_MASS`      | `140.90765` u        | Pr mass for Doppler width                    |
| `T`                | `1.6` eV             | Plasma temperature for Doppler width         |
| match tolerance    | `5.5 σ`              | Candidate window in `match_and_grade`        |
| `alpha`            | `0.5`                | Factor-update damping in `weed_assignments`  |
| `min_n`            | `5`                  | Min. qualifying transitions per level factor |
| `CLAMP`            | `5.0`                | Max. magnitude of a per-level ln-factor      |
| `max_iterations`   | `100` (weeding)      | Inner weeding iteration cap                  |
| `max_cycles`       | `10` (`main`)        | Outer match→weed→optimize cycles             |

---

## Project Structure

```
LineClass/
├── README.md                     ← this file
├── classify_lines.py             # Full pipeline: read → generate → match → grade → resolve → weed → optimize → output
├── models.py                     # Dataclasses: EnergyLevel, SpectralLine, Transition, UNASSIGNED
├── Icalc.xlsx                    # Input: calculated transition intensities & uncertainties
├── Pr3_lines.xlsx                # Input: observed spectral lines (Sugar 1969/1974)
├── line_classifications.xlsx     # Output: classification table (LOPT-ready)
└── line_classifications.csv      # Output: same data as CSV

# Shared / external to this directory:
../TableExtraction/Pr3_lev_Wyart_1999.xlsm   # Input: energy levels (sheet 'Wyart2000')
../statistics.py                             # Bundled mandel_paule() (shadows stdlib 'statistics')
```

---

## Pipeline Architecture

```
┌─────────────────────┐    ┌────────────────────┐    ┌───────────────────────┐
│  Read Levels,       │───▶│  Generate & Match  │───▶│  Resolve Conflicts    │
│  Lines, Icalc       │    │  (binary search,   │    │  (conservative sort,  │
│  (Steps 1–4)        │    │   2D grading)      │    │   F/R flags)          │
└─────────────────────┘    └────────────────────┘    └──────────┬────────────┘
                                                                 │
                                                  ┌──────────────▼────────────┐
                                                  │  Iterative Weeding        │
                                                  │  (per-level factors,      │
                                                  │   3-step per-line logic,  │
                                                  │   blacklist oscillations) │
                                                  └──────────┬────────────────┘
                                                             │
                                              ┌──────────────▼────────────────┐
                                              │  Energy Level Optimization    │
                                              │  (weighted mean, BF weights)  │
                                              └──────────┬────────────────────┘
                                                         │ converged?
                                              ┌──────────▼──────────────┐
                                              │  LOPT-ready output      │
                                              │  (xlsx + csv)           │
                                              └─────────────────────────┘
```

The outer loop repeats until both the number of accepted transitions and the level energies are stable across cycles (first-iteration max energy change `< 0.001 cm⁻¹`).

---

## License

For academic and research purposes, as part of the parent `TableExtraction` project.
