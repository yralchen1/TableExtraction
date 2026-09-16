# TableExtraction

An AI-powered pipeline for extracting spectral line tables from scanned scientific PDFs using **Google Gemini** multimodal vision and **LangChain/LangGraph**.

Originally built for digitizing tables from **Sugar's Pr III** spectroscopy papers ([1969](https://doi.org/10.6028/jres.073A.029) & [1974](https://doi.org/10.6028/jres.078A.036)), the pipeline converts PDF pages to images, sends them to Gemini for OCR/extraction, validates against physics formulas and [NIST ASD](https://physics.nist.gov/asd) reference energy levels, corrects page-level misalignments, and outputs structured Excel/CSV files. A companion `LineClass` module provides a fully-automated, deterministic pipeline for classifying observed spectral lines — including iterative weeding of spurious classifications and energy level optimization — and produces output ready for use with the LOPT level-optimization program. A statistical **validation suite** in the same module measures how often the classifier would confirm an energy level that is not real (by planting displaced "decoy" copies of the tested levels into the real run) and assigns each re-established level a probability of being spurious, combining the level's energy stability with the agreement between its observed and theoretically predicted intensity pattern.

---

## Features

- **Multimodal OCR** — Extracts structured data from scanned two-column table layouts using Gemini vision models (Pro & Flash)
- **Physics-aware validation** — Validates air↔vacuum wavelength/wavenumber via Edlén's 1966 formula
- **NIST ASD cross-reference** — Looks up energy levels against NIST Atomic Spectra Database with auto-correction for typos, parity flips, and J-value mismatches
- **Selection rule checks** — Enforces E1 transition rules (parity change, ΔJ ≤ 1)
- **Misalignment correction** — Detects and fixes OCR-induced row shifts using Ritz wavenumber matching
- **Line Classification** — Automated matching of observed lines to Pr III energy levels with 2D grading (tier by wavenumber agreement in σ units; subgrade by log intensity-ratio consistency)
- **Iterative Intensity Weeding** — A 3-step per-line decision process (clear-cut tests, blend grouping by center-of-gravity, conservative defaults) iterated to convergence; assignments whose accept/reject decision oscillates between iterations or cycles are blacklisted and rejected. Optional per-level intensity calibration factors are available but disabled by default
- **Level-Energy Uncertainty Estimation** — Each level's energy uncertainty is estimated once per outer cycle from the scatter of energies implied by its currently accepted transitions (propagated through chains of new levels by fixed-point iteration), and folded into the matching window and Ritz-mismatch test, so lines are not rejected just because a not-yet-precisely-known level's adopted energy has a small, unaccounted-for error
- **Energy Level Optimization** — Weighted least-squares refinement of level energies from accepted transitions, iterated until both accepted count and energy values converge; output is LOPT-ready
- **Classification Validation** — Decoy (shadow-level) runs measure the false-confirmation rate of the pipeline under fully realistic conditions; each re-established level receives a spurious probability `p_spur` that folds two independent pieces of evidence: the stability of its optimized energy and its observed-vs-predicted intensity pattern
- **Robust JSON parsing** — 6-step fallback parser for malformed LLM output
- **Configurable** — Year-specific prompts, column contexts, layout definitions, and environment-driven settings

---

## Prerequisites

- **Python 3.10+** (LineClass uses PEP 604 `X | Y` type annotations; the bundled `statistics.py` in the repository root, which shadows the standard-library `statistics` module, is needed only if the optional Mandel–Paule intensity-factor weighting is enabled in LineClass)
- **poppler** (for PDF → image conversion):
  ```bash
  # macOS
  brew install poppler

  # Ubuntu/Debian
  sudo apt-get install poppler-utils

  # Windows — download from https://github.com/oschwartz10612/poppler-windows/releases
  ```
- **Google Gemini API key** — get one free at [https://aistudio.google.com/apikey](https://aistudio.google.com/apikey)

---

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/yourusername/TableExtraction.git
cd TableExtraction/TableExtraction
```

### 2. Create & activate a virtual environment

```bash
python3 -m venv venv
source venv/bin/activate        # macOS / Linux
# venv\Scripts\activate         # Windows
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Configure your API key

```bash
cp .env.example .env
```

Edit `.env` and add your Gemini API key:
```
GOOGLE_API_KEY=your_api_key_here
```

### 5. (Optional) Choose a Gemini model

Edit `.env`:
```
GEMINI_MODEL=gemini-2.5-pro            # Best accuracy, higher cost / stricter quota
# GEMINI_MODEL=gemini-2.5-flash        # Fast and cheap
# GEMINI_MODEL=gemini-2.0-flash        # Alternative
# GEMINI_MODEL=gemini-3-flash-preview  # High accuracy, moderate cost
```

Or pass as a CLI flag: `--model gemini-2.5-pro`

> **Note:** `gemini-2.5-pro` produces significantly better table alignment than the flash models, especially for two-column layouts.

---

## Usage

```bash
# Full pipeline (all pages):
python extract_table.py

# Specific page range:
python extract_table.py --pages 3 39

# Single page (for testing):
python extract_table.py --single-page 5

# Use a specific model:
python extract_table.py --pages 3 39 --model gemini-2.5-pro

# Reuse existing page images (skip PDF conversion):
python extract_table.py --skip-images --model gemini-2.5-pro

# Lower DPI for faster processing:
python extract_table.py --dpi 200 --pages 3 39

# Dry run (show config, no processing):
python extract_table.py --dry-run --pages 3 39

# Convert PDF to images only:
python pdf_to_images.py --year 1974 --pages 3 39

# Validate a previously extracted raw table:
python validate_extracted_table.py --input output/raw_extracted_table.xlsx

# --- Spectral Line Classification ---

# Run the full classification pipeline:
python LineClass/classify_lines.py

# --- Classification validation (run from LineClass/, after classify_lines.py) ---
python LineClass/decoy_mc.py       # decoy runs → in-situ false-confirmation rates
python LineClass/chance_mc.py      # optional shifted-wavenumber cross-check
python LineClass/level_shifts.py   # validation report + per-level spurious probabilities
python LineClass/level_shifts.py --detail 059003.000483   # inspect one level

# --- Revising levels and finding new ones (run from LineClass/) ---
python level_interchange.py        # levels whose theoretical identities may be interchanged
python level_positions.py --scan   # alternate positions for the questionable levels
python check_sync.py               # do all the files still describe the same identification?
python sync_IDEN2.py               # bring the IDEN2 working files back into step
python unfound_levels.py           # which never-found levels are worth searching for
python level_positions.py --unknown 742      # search for one of them
python find_unknown_levels.py                # ... or for the whole promising list at once
python insert_new_level.py --iden2-row 742   # put an accepted position into the pipeline
python move_level.py 059003.000565 --energy 118967.3776   # ... or move one already in
```

---

## Output

Results are saved to `output/`:

| File | Description |
|------|-------------|
| `*_extracted.csv` | CSV with all output columns |
| `*_extracted.xlsx` | Excel with auto-sized columns |

> **CSV J-values:** Half-integer J values (e.g., `5/2`) are written as `="5/2"` in CSV to prevent Excel from misinterpreting them as dates. The Excel file stores them natively as text.

### LineClass Output

| File | Description |
|------|-------------|
| `LineClass/line_classifications.xlsx` | Classification results, sorted by decreasing wavenumber |
| `LineClass/line_classifications.csv` | Same data in CSV format |
| `LineClass/level_shift_report.csv/.xlsx` | Validation: one row per level — energy shift, support counts, intensity-pattern scores, the spurious probabilities `p_spur_decoys`, `p_spur_pattern`, `p_spur`, and the adjudication columns `question_status`/`reason` (questionable marks cleared when the missing predicted lines are masked by nearby stronger lines) |
| `LineClass/decoy_mc_*.csv/.xlsx` | Validation: decoy-run tables (per-decoy support and energy wander, perturbation check, accepted decoy lines, per-run summary) |
| `LineClass/chance_mc_*.csv/.xlsx` | Validation (optional): shifted-wavenumber cross-check tables |
| `LineClass/level_interchange.csv/.xlsx` | Pairs of levels whose theoretical identities may be interchanged |
| `LineClass/level_positions.csv` | Alternate positions offered to each questionable level by `level_positions.py --scan` |
| `LineClass/unfound_levels.csv/.xlsx` | The calculated levels nobody has found, ranked by how many recordable transitions they would have |
| `LineClass/found_levels.csv` | One row per searched level: where `find_unknown_levels.py` found a position, and how well supported it is |
| `LineClass/line_decisions.csv` | The decision ledger — hand-ruled `accept`/`reject` verdicts on individual assignments, with dates and reasons |
| `LineClass/line_decisions_removed.csv` | The ledger rows a level move took out, kept with the date, the level and the reason they no longer apply |
| `LineClass/revised_level_energies.csv` | The levels the identification work has moved, as `level_id, E_input` — how a level of the published list is re-positioned without editing the published workbook |
| `LineClass/LOPT_input_lines.txt` | The accepted transitions in LOPT's fixed-column input format, with their branching-fraction weights |
| `LineClass/sync_report.txt` | `check_sync.py`'s findings: whether the pipeline, LOPT and IDEN2 files still describe the same identification |

All validation tables are Excel-friendly: floats are rounded to meaningful decimals and each CSV has an `.xlsx` twin in which text cells such as J = `3/2` are not converted to dates.

#### Output Columns

| Column | Description |
|--------|-------------|
| `wn_obs` | Observed wavenumber (cm⁻¹), 3 decimal places |
| `unc_wn_obs` | Wavenumber uncertainty (cm⁻¹) |
| `obs_intens` | Observed intensity |
| `char` | Line character (`h`, `w`, `bl`, etc.) |
| `low_id` | Lower energy level ID |
| `upp_id` | Upper energy level ID |
| `calc_intens` | Adjusted calculated intensity (after per-level factors) |
| `orig_calc_intens` | Original theoretical intensity from `Icalc.xlsx` |
| `u_calc` | Adjusted log-uncertainty of `calc_intens` |
| `imputed` | `1` = the pair is absent from `Icalc.xlsx` and its intensity was imputed |
| `intens_from_f` | Upper-level intensity correction factor (ln scale) |
| `intens_to_f` | Lower-level intensity correction factor (ln scale) |
| `dif_wn_O-C` | Observed − Ritz wavenumber difference (cm⁻¹) |
| `grade` | 2D grade: tier (2–5) + subgrade (`A`–`E`, `G`); see grading scheme below |
| `notes1` | Conflict flags: `F` (conflicting), `R` (revised from original) |
| `notes2` | Per-transition decision trace from weeding (Step1/Step2/Step3 label) |
| `manual` | `accept`/`reject` if the decision ledger (`line_decisions.csv`) rules on this assignment |
| `new` | `1` = new classification, `0` = original Sugar 1969/1974 classification, blank = unclassified |
| `accepted` | `1` = accepted by weeding, `0` = rejected, blank = unclassified line |
| `n_accepted` | Number of accepted classifications of this observed line |
| `low_E` | Lower level energy (cm⁻¹) |
| `upp_E` | Upper level energy (cm⁻¹) |
| `rwn` | Ritz wavenumber = upp_E − low_E (cm⁻¹) |
| `BF` | Branching fraction of this component (0 if not accepted); the LOPT centroid weight is `BF²` |

> Rows where `accepted = 1` with their `BF` values form direct input for LOPT (which squares `BF` internally to obtain the centroid weight).

### Extraction Output Columns

| # | Column | Description |
|---|--------|-------------|
| 1 | Wavelength (Å) | Air wavelength, 3 decimal places |
| 2 | Intensity | Integer intensity value |
| 3 | Line Character | Symbols (h, w, bl, etc.) |
| 4 | Wavenumber (cm⁻¹) | Vacuum wavenumber, 2 decimal places |
| 5 | Lower Level Int | Integer part of lower energy level |
| 6 | Lower Parity | `e` (even) or `o` (odd) |
| 7 | Lower J | J-value (e.g., "9/2") |
| 8 | Upper Level Int | Integer part of upper energy level |
| 9 | Upper Parity | `e` (even) or `o` (odd) |
| 10 | Upper J | J-value of upper level |
| 11 | Lower Level NIST (cm⁻¹) | Exact level from NIST ASD |
| 12 | Upper Level NIST (cm⁻¹) | Exact level from NIST ASD |
| 13 | Page | Source PDF page number |
| 14 | Notes | All corrections, errors, and flags |

---

## Validation Pipeline

### Validation Order

Each row is validated in a specific order to avoid incorrect corrections:

1. **NIST ASD lookup** — Looks up both energy levels with auto-correction
2. **Selection rule checks** — Parity change + ΔJ ≤ 1 (only when both levels are found)
3. **Wavelength/wavenumber check** — Only corrects if NIST lookup passed
4. **Page-level misalignment correction** — Relocates misplaced classifications

> **Key insight:** If the NIST lookup fails, the classification is likely misaligned by OCR, so the wavelength/wavenumber is probably correct as-is. The pipeline flags "suspect misalignment" without modifying values, then attempts a Ritz-based fix.

### NIST ASD Auto-Correction

When a level isn't found exactly, the pipeline tries (in order):
- ±1 in integer part with exact parity + J
- Flip parity (`e` ↔ `o`) — catches ° misread as digit 0
- Different J from NIST — for levels where NIST revised the J value
- Single-digit typo in the integer part
- Adjacent-digit swap in the integer part

### Selection Rules (E1 Transitions)

- **Parity must change:** lower and upper levels must have different parity
- **ΔJ ≤ 1:** `|J_upper - J_lower| ≤ 1`
- Only checked when **both** levels are found in NIST

### Wavelength ↔ Wavenumber (Edlén's 1966 Formula)

Converts air wavelength to vacuum wavenumber using Edlén's formula. On mismatch:
- **NIST OK** → tries typo correction (single digit, adjacent swap)
- **NIST failed** → marks as "suspect misalignment" — no correction applied

### Page-Level Misalignment Correction

When a page has **> 2 error rows**, the pipeline hypothesizes OCR shifted some classifications to wrong rows. For each error row with NIST-confirmed levels:
1. Compute **Ritz wavenumber** = upper_exact − lower_exact
2. Search same page for the row whose observed wavenumber is closest
3. If the target row is unclassified → relocate the classification ✅
4. If already classified → flag as unresolvable ❌

---

## Project Structure

```
TableExtraction/
├── README.md                              ← this file
├── .env.example                           # API key template
├── .gitignore
├── astools.py                             # Atomic spectroscopy tools (ASTools library)
├── test_astools.py                        # Unit tests for ASTools
├── prompts/
│   ├── column_contexts/                   # Year-specific column definitions (JSON)
│   └── layout_definitions/                # Table layout prompt templates
├── output/                                # CSV/Excel output
├── TableExtraction/                       # Core extraction pipeline
│   ├── config.py                          # Configuration (env-driven)
│   ├── extract_table.py                   # CLI entry point (extraction)
│   ├── validate_extracted_table.py        # CLI entry point (validation only)
│   ├── extract_page_direct.py             # Quick single-page raw OCR
│   ├── extraction_prompt.py               # Gemini prompt builder
│   ├── graph.py                           # LangGraph pipeline (nodes & state)
│   ├── tools.py                           # Validation tools (Edlén, NIST, selection rules)
│   ├── pdf_to_images.py                   # PDF → PNG conversion
│   ├── Pr3_lev_Wyart_1999.xlsm           # Energy levels for Pr III (Wyart 1999) — shared input for LineClass
│   └── requirements.txt                   # Python dependencies
└── LineClass/                             # Spectral Line Classification + validation suite
    ├── README.md                          # LineClass component documentation (incl. validation)
    ├── classify_lines.py                  # Full automated pipeline: match → weed → optimize → LOPT output
    ├── models.py                          # Data models (EnergyLevel, SpectralLine, Transition)
    ├── config.py                          # Reads lineclass_config.toml; resolves every path relative to it
    ├── output_files.py                    # Writable-file preflight (catches workbooks held open by Excel)
    ├── cowan_gA.py                        # Cowan calculation: gA values, level matching to IDEN2
    ├── gA_imputation.py                   # Imputed intensity for pairs absent from Icalc.xlsx
    ├── make_LOPT_input.py                 # LOPT input rows and their branching-fraction weights
    ├── decoy_mc.py                        # Validation: decoy (shadow-level) runs
    ├── level_shifts.py                    # Validation: calibrations, pattern scores, p_spur; --detail mode
    ├── chance_mc.py                       # Shared utilities + optional shifted-wavenumber cross-check
    ├── level_interchange.py               # Levels whose theoretical identities may be interchanged
    ├── swap_line_assignments*.py          # Repair an accepted interchange across all three file sets
    ├── level_positions.py                 # --scan: alternate positions; --unknown: search for a lost level
    ├── unfound_levels.py                  # Never-found levels, ranked by recordable transitions
    ├── find_unknown_levels.py             # The --unknown search run down that whole list → found_levels.csv
    ├── insert_new_level.py                # Put an accepted position into the pipeline, LOPT and IDEN2
    ├── move_level.py                      # Move an assigned level, releasing the lines it leaves behind
    ├── check_sync.py                      # Do the pipeline, LOPT and IDEN2 files still agree?
    ├── sync_IDEN2.py                      # Bring the IDEN2 working files back into step
    ├── lineclass_config.toml              # Paths, thresholds and model constants
    ├── Pr3_lines.xlsx                     # Observed spectral lines (Sugar 1969/1974)
    ├── Icalc.xlsx                         # Calculated transition intensities and uncertainties
    ├── line_decisions.csv                 # Decision ledger: hand-ruled accept/reject verdicts
    ├── new_levels.txt                     # Levels added by hand, tab-delimited
    ├── IDEN2/                             # The IDEN2 working files (enlev.dat, trans.dat, IDEN_level_ids.txt)
    ├── tools/                             # Calibration and diagnostic scripts (intensities, coverage, SNR)
    ├── line_classifications.xlsx          # Classification output (LOPT-ready)
    ├── line_classifications.csv           # Classification output (CSV)
    ├── LOPT_input_lines.txt               # LOPT input; LOPT_output_levels.txt / _lines.txt come back
    ├── level_shift_report.csv/.xlsx       # Validation output: per-level table with p_spur
    ├── decoy_mc_*.csv/.xlsx               # Validation output: decoy-run tables
    ├── chance_mc_*.csv/.xlsx              # Validation output: cross-check tables (optional)
    ├── level_interchange.csv/.xlsx        # Interchange candidates
    ├── level_positions.csv                # Alternate positions for the questionable levels
    ├── unfound_levels.csv/.xlsx           # The ranked search list
    ├── found_levels.csv                   # What the searches found
    └── sync_report.txt                    # check_sync.py findings
```

---

## Pipeline Architecture

```
┌─────────────┐    ┌─────────────────┐    ┌──────────────────────┐
│ convert_pdf  │───▶│ extract_tables  │───▶│   validate_rows      │
│ (pdf2image)  │    │ (Gemini Vision) │    │ (NIST + rules + λ/σ) │
└─────────────┘    └─────────────────┘    └──────────┬───────────┘
                                                      │
                                          ┌───────────▼──────────┐
                                          │  fix_misalignment    │
                                          │  (Ritz-based move)   │
                                          └──────────┬───────────┘
                                                      │
                                          ┌───────────▼──────────┐
                                          │   compile_output     │
                                          │   (CSV + Excel)      │
                                          └──────────────────────┘
```

Built with [LangGraph](https://github.com/langchain-ai/langgraph) `StateGraph`.

### LineClass Pipeline Architecture

```
┌─────────────────────┐    ┌────────────────────┐    ┌───────────────────────┐
│  Read Levels,       │───▶│  Generate & Match  │───▶│  Resolve Conflicts    │
│  Lines, Icalc       │    │  (binary search,   │    │  (conservative sort,  │
│                     │    │  2D grading)       │    │  F/R flags)           │
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
                                              │  (weighted least squares)     │
                                              └──────────┬────────────────────┘
                                                         │ converged?
                                              ┌──────────▼──────────────┐
                                              │  LOPT-ready output      │
                                              │  (xlsx + csv)           │
                                              └─────────────────────────┘
```

The outer loop repeats until both the number of accepted transitions and the maximum level energy change (< 0.001 cm⁻¹) are stable across cycles. At the start of each cycle (not shown), `compute_level_uncertainties` re-estimates every level's `u_energy` from the previous cycle's accepted transitions, before `Generate & Match` runs.

---

## Line Classification Pipeline (`LineClass`)

The `LineClass` module provides a fully-automated, deterministic pipeline for classifying observed spectral lines of Pr III against all theoretically allowed transitions. The pipeline supersedes the previous LLM-based post-processing step with a convergent iterative algorithm. A main outer loop repeats the full assignment cycle until both the accepted transition count and the refined energy levels have stabilized.

### Input Files

| File | Description |
|------|-------------|
| `Pr3_lev_Wyart_1999.xlsm` | Energy levels (Wyart 1999) |
| `Icalc.xlsx` | Theoretical transition intensities and uncertainties |
| `Pr3_lines.xlsx` | Observed spectral lines with existing Sugar classifications |

### 6-Step Workflow

**Steps 1–4** prepare the input data:

1. **Read Energy Levels** — Loads J, energy, and parity from `Pr3_lev_Wyart_1999.xlsm`
2. **Read Calculated Transitions** — Loads intensities and log-uncertainties from `Icalc.xlsx`
3. **Read Observed Lines** — Loads wavenumbers, intensities, and original Sugar assignments from `Pr3_lines.xlsx`
4. **Generate Transitions** — Enumerates all possible transitions satisfying E1 selection rules (|ΔJ| ≤ 1, parity change, wavenumber in [9327, 121665] cm⁻¹)

**Step 5** runs the assignment cycle:

**5.1 Match & Grade** — Binary search finds all possible transitions within `5.5 σ_comb` of each observed line, where `σ_comb` combines the line's own wavenumber uncertainty with the estimated energy uncertainty of both candidate levels (`u_energy`, re-estimated once per outer cycle — see *Level-Energy Uncertainty Estimation* below). Each candidate is assigned a 2D grade: **tier** (2–5) by wavenumber agreement in the line's own σ (not `σ_comb`); **subgrade** (A–G) by intensity consistency (log-ratio of I_obs to I_calc vs. `u_calc`).

**5.2 Resolve Conflicts** — When a transition is matched to multiple observed lines, a conservative sort selects the best match (prefer legacy over new, minimize |Obs−Ritz|, then tier). Losers are removed. The `notes1` field records `F` (conflicting) and `R` (revised).

**5.3 Iterative Intensity Weeding** — See subsection below.

**Step 6** runs energy level optimization — see subsection below.

### Iterative Weeding Algorithm

`weed_assignments()` decides, for each observed line, which of its candidate transitions to accept or reject. It uses per-level intensity calibration factors derived from already-accepted transitions, updated iteratively.

**Outer iteration loop:**
- **Iteration 0:** Baseline pass using raw theoretical intensities from `Icalc.xlsx`
- **Iterations 1–N:** re-weed until no transition changes its acceptance state (or `max_iterations` is reached). An optional per-level intensity calibration (factors computed from accepted unblended transitions by backfitting, with damped updates) can modify the theoretical intensities between iterations; it is **disabled by default** (`USE_INTENSITY_ADJUSTMENT = 0`), having been found to add only ~8 accepted classifications at the cost of ~240 fitted parameters
- **Oscillation blacklist:** transitions whose accept/reject decision flips between weeding iterations, or whose end-of-cycle state alternates between outer cycles, are blacklisted and permanently rejected — an oscillating decision is treated as evidence of an unreliable assignment. The blacklist persists across the outer cycles

**Per-line 3-step decision process (`weed_assignments_line`):**

Each line's candidate transitions are sorted by decreasing `calc_intensity`, then by Ritz σ.

- **Step 1 — Clear-cut decisions:** New transitions with Ritz mismatch > 3σ_comb are rejected (legacy transitions with no theoretical intensity are rejected only beyond 5σ_comb); `σ_comb` is the same combined line+level-uncertainty sigma used for matching. A relative-intensity filter removes new transitions contributing < 10% of total `calc_intensity` (< 4% for legacy). Surviving candidates are tested by z-score against I_obs: asymmetric thresholds apply (looser for resonance lines and for cases where `I_calc > I_obs`). Decisions are labeled "Step1" in the `notes2` field.

- **Step 2 — Grouping and late arbitration:** Undecided candidates are tested for pair/triple acceptance: groups pass if their joint center-of-gravity matches the observed wavenumber and their combined effective spread (Doppler width convolved with measurement uncertainty) is consistent with the line profile. For lines flagged `h` or `w`, broadening thresholds are relaxed. Remaining undecided new candidates with `F` or `R` in `notes1` are rejected.

- **Step 3 — Default decisions:** Undecided **new** candidates: rejected ("no solid evidence for acceptance"). Undecided **legacy** candidates: accepted ("no solid evidence for rejection").

### Level-Energy Uncertainty Estimation

Before each outer cycle's match & grade step, `compute_level_uncertainties()` (re-)estimates every level's energy uncertainty `u_energy` from the scatter of energies implied by its currently accepted transitions (using the previous cycle's accepted set — empty on the first cycle, so `u_energy` starts at 0 for every level). Because a partner level's own `u_energy` is itself part of the unknown, the whole level system is solved by repeated sweeps (Gauss-Seidel-style fixed-point iteration) until it stabilizes, letting uncertainty propagate outward from well-established old levels through however many links a chain of new levels has, without needing to identify a single transition that "defines" each level. `u_energy` then widens the matching window and the Step-1 Ritz-mismatch test (above) via `σ_comb`, so a transition is not rejected merely because a not-yet-precisely-known level's adopted energy carries an unstated error. See [`LineClass/README.md`](LineClass/README.md) for the full algorithm.

### Energy Level Optimization

After weeding, `optimize_levels()` refines each non-ground energy level by a weighted mean of all implied energies from accepted transitions (both upper-from-lower and lower-from-upper). Weights are `1/σ²` for unblended lines; for blended lines the weight is shared proportionally to `calc_intensity`. The iteration runs until max energy change < 0.001 cm⁻¹. The `BF` column in the output is each accepted component's branching fraction; LOPT squares it internally to obtain the centroid weight.

### Classification Validation

Many of the accepted classifications re-establish energy levels that were never confirmed in print (Wyart's unpublished levels), so the reliability of the result must be measured, not assumed. Three scripts in `LineClass/` do this (details in [`LineClass/README.md`](LineClass/README.md)):

- **`decoy_mc.py`** repeats the real classification with one **decoy** added per tested level: an exact copy (same J, parity, possible transitions, and predicted intensities) whose energy is displaced far enough that none of its transitions can coincide with a true line. The decoys compete with the real levels in every step, so every line accepted for a decoy is a false match obtained under fully realistic conditions; with one decoy per tested level, the number of decoys passing any acceptance rule equals the number of false confirmations to expect if all tested levels were fake.
- **`chance_mc.py`** provides shared utilities for the suite and, run directly, an independent cross-check in which all observed wavenumbers are shifted so that every accepted match is false; its rates are upper bounds and agree with the decoy rates within ~25%.
- **`level_shifts.py`** turns the calibrations into per-level verdicts. Each tested level is judged by two independent pieces of evidence: the stability of its optimized energy (a genuine level returns almost exactly to its input value, because its lines tie it to well-anchored known levels; a level built from chance coincidences drifts away) and its intensity pattern (a real level's strongest theoretically predicted transitions must be present among its accepted lines — the classical "square-array" argument). Both are calibrated on the known-genuine old levels and on the known-fake decoys, and folded into a per-level **probability of being spurious** (`p_spur`), written to `level_shift_report.csv/.xlsx` together with the separate energy-only and pattern-only probabilities. A `--detail` mode prints any single level's predicted transitions with the fate of each in the run, for case-by-case inspection.

A level the suite puts in doubt then has to be dealt with, and a level nobody has ever found has to be looked for. Four further scripts do that work, all of them documented in [`LineClass/README.md`](LineClass/README.md):

- **`level_interchange.py`** looks for pairs of levels of the same J and parity whose theoretical identities may simply have been swapped — the energies are right and the labels are on the wrong rows. `swap_line_assignments.py` repairs an accepted swap in all three files at once.
- **`level_positions.py --scan`** offers each questionable level the alternate positions its lines would also allow, so that a level can be moved rather than merely doubted; `--unknown ROW` turns the same machinery outwards and searches the window a calculated level could occupy for a position the recorded lines support.
- **`unfound_levels.py`** ranks the calculated levels that have never been placed by the number of their transitions that would have been recorded on the plates, and **`find_unknown_levels.py`** runs the `--unknown` search down that ranked list in one pass, reducing each search to one row of `found_levels.csv`.
- **`insert_new_level.py`** takes an accepted position through everything that must follow it — the calculated transitions, the LOPT input rows and their weights, the LOPT run, the classification run, the ledger rows and the IDEN2 files — and restores every file if the fit rejects the level. **`move_level.py`** does the same for a level that is already in and has to move: it takes the old position apart first, releasing every line the new position cannot account for, with the reason it failed on written into the ledger, and keeping the ledger rows it removes in `line_decisions_removed.csv`.

Every one of these proposes; none of them decides. Each changed assignment is reviewed by hand in IDEN2 against the plates and the branch structure before it is accepted.

### Grading Scheme

Matches are graded on two dimensions:

-   **Tier (Wavenumber Agreement):**
    -   `2`: Wavenumber residual ≤ 2.0σ (uncertainty)
    -   `3`: Wavenumber residual ≤ 3.0σ
    -   `4`: Wavenumber residual ≤ 4.0σ
    -   `5`: Wavenumber residual > 4.0σ (up to the 5.5σ match tolerance)
-   **Subgrade (Intensity Consistency):**
    Assigned based on agreement between observed and calculated intensities (I_obs, I_calc) and uncertainty `u_calc = ln(uA_pcnt/100 + 1)`. Let `ln_I_ratio = |ln(I_calc/I_obs)|`.
    -   `G`: No theoretical intensity available
    -   `A`: (ln_I_ratio - u_calc ≤ 0) AND (u_calc ≤ ln(2))
    -   `B`: (ln_I_ratio - u_calc ≤ ln(2)) AND (ln(2) < u_calc ≤ ln(5))
    -   `C`: (ln(2) < ln_I_ratio - u_calc ≤ ln(5)) AND (ln(2) < u_calc ≤ ln(5))
    -   `D`: (ln_I_ratio - u_calc ≤ ln(2)) AND (u_calc > ln(5))
    -   `E`: (ln_I_ratio - u_calc > ln(5)) AND (u_calc > ln(5))

**Notes:**
-   `N`: Newly assigned (not in original source)
-   `F`: Conflicting assignment (multiple lines match one transition)
-   `R`: Revised (original classification was changed or line is now unassigned)
-   `notes2`: Per-transition decision trace from weeding (e.g., "Step1: intensity ratio too low", "Step3: legacy default accepted")

---

## ASTools Library

The repository includes `astools.py` — a Python translation of the VBA module `ASTools.bas` for atomic spectroscopy computations:

- **Physical constants** across CODATA vintages (1932–2022)
- **Refractive index of air** (Edlén 1953/1966, Peck & Reeder 1972, Meggers & Peters 1919, Barrell 1951, Ciddor 1996)
- **Air↔vacuum wavelength conversions** using multiple formulations
- **Water vapor pressure** calculations (Stone-Zimmerman, Buck, Saul & Wagner)

---

## Dependencies

| Package | Purpose |
|---------|---------|
| `langchain-google-genai` | Gemini API integration |
| `langchain-core` | LangChain core |
| `langgraph` | Pipeline orchestration |
| `pdf2image` | PDF → PNG conversion |
| `Pillow` | Image handling |
| `openpyxl` | Excel I/O |
| `python-dotenv` | Environment config |
| `pandas` | Data manipulation |
| `numpy` | Numerical computations in LineClass weeding and optimization |

---

## License

This project is for academic and research purposes.
