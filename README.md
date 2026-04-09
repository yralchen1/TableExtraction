# TableExtraction

An AI-powered pipeline for extracting spectral line tables from scanned scientific PDFs using **Google Gemini** multimodal vision and **LangChain/LangGraph**.

Originally built for digitizing tables from **Sugar's Pr III** spectroscopy papers ([1969](https://doi.org/10.6028/jres.073A.029) & [1974](https://doi.org/10.6028/jres.078A.036)), the pipeline converts PDF pages to images, sends them to Gemini for OCR/extraction, validates against physics formulas and [NIST ASD](https://physics.nist.gov/asd) reference energy levels, corrects page-level misalignments, and outputs structured Excel/CSV files. A companion `LineClass` module provides a fully-automated, deterministic pipeline for classifying observed spectral lines — including iterative intensity weeding and energy level optimization — and produces output ready for use with the LOPT level-optimization program.

---

## Features

- **Multimodal OCR** — Extracts structured data from scanned two-column table layouts using Gemini vision models (Pro & Flash)
- **Physics-aware validation** — Validates air↔vacuum wavelength/wavenumber via Edlén's 1966 formula
- **NIST ASD cross-reference** — Looks up energy levels against NIST Atomic Spectra Database with auto-correction for typos, parity flips, and J-value mismatches
- **Selection rule checks** — Enforces E1 transition rules (parity change, ΔJ ≤ 1)
- **Misalignment correction** — Detects and fixes OCR-induced row shifts using Ritz wavenumber matching
- **Line Classification** — Automated matching of observed lines to Pr III energy levels with 2D grading (tier by wavenumber agreement in σ units; subgrade by log intensity-ratio consistency)
- **Iterative Intensity Weeding** — Per-level intensity adjustment factors computed from accepted transitions (Mandel-Paule weighted mean with damped updates) drive iterative re-weeding until convergence; oscillating assignments are blacklisted
- **Energy Level Optimization** — Weighted least-squares refinement of level energies from accepted transitions, iterated until both accepted count and energy values converge; output is LOPT-ready
- **Robust JSON parsing** — 6-step fallback parser for malformed LLM output
- **Configurable** — Year-specific prompts, column contexts, layout definitions, and environment-driven settings

---

## Prerequisites

- **Python 3.11+** (required for `statistics.mandel_paule` used by the LineClass weeding algorithm)
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
| `intens_from_f` | Upper-level intensity correction factor (ln scale) |
| `intens_to_f` | Lower-level intensity correction factor (ln scale) |
| `dif_wn_O-C` | Observed − Ritz wavenumber difference (cm⁻¹) |
| `grade` | 2D grade: tier (2–5) + subgrade (A–G); see grading scheme below |
| `notes1` | Conflict flags: `F` (conflicting), `R` (revised from original) |
| `notes2` | Per-transition decision trace from weeding (Step1/Step2/Step3 label) |
| `new` | `1` = new classification, `0` = original Sugar 1969/1974 classification |
| `accepted` | `1` = accepted by weeding, `0` = rejected, blank = unclassified line |
| `low_E` | Lower level energy (cm⁻¹) |
| `upp_E` | Upper level energy (cm⁻¹) |
| `rwn` | Ritz wavenumber = upp_E − low_E (cm⁻¹) |
| `weight` | LOPT-compatible weight for this transition |

> Rows where `accepted = 1` with their `weight` values form direct input for LOPT.

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
└── LineClass/                             # Spectral Line Classification
    ├── classify_lines.py                  # Full automated pipeline: match → weed → optimize → LOPT output
    ├── models.py                          # Data models (EnergyLevel, SpectralLine, Transition)
    ├── Pr3_lines.xlsx                     # Observed spectral lines (Sugar 1969/1974)
    ├── Icalc.xlsx                         # Calculated transition intensities and uncertainties
    └── line_classifications.xlsx          # Classification output (LOPT-ready)
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

The outer loop repeats until both the number of accepted transitions and the maximum level energy change (< 0.001 cm⁻¹) are stable across cycles.

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

**5.1 Match & Grade** — Binary search finds all possible transitions within 5.5σ of each observed line. Each candidate is assigned a 2D grade: **tier** (2–5) by wavenumber agreement in σ units; **subgrade** (A–G) by intensity consistency (log-ratio of I_obs to I_calc vs. `u_calc`).

**5.2 Resolve Conflicts** — When a transition is matched to multiple observed lines, a conservative sort selects the best match (prefer legacy over new, minimize |Obs−Ritz|, then tier). Losers are removed. The `notes1` field records `F` (conflicting) and `R` (revised).

**5.3 Iterative Intensity Weeding** — See subsection below.

**Step 6** runs energy level optimization — see subsection below.

### Iterative Weeding Algorithm

`weed_assignments()` decides, for each observed line, which of its candidate transitions to accept or reject. It uses per-level intensity calibration factors derived from already-accepted transitions, updated iteratively.

**Outer iteration loop:**
- **Iteration 0:** Baseline pass using raw theoretical intensities from `Icalc.xlsx`
- **Iterations 1–N:**
  1. Compute per-level intensity factors from currently accepted, unblended transitions (minimum 5 qualifying entries per level; Mandel-Paule weighted mean with chi²-inflated fallback; result clamped to ±5 ln-units)
  2. Apply damped update: `factor_new = alpha × factor_computed + (1−alpha) × factor_prev` (default `alpha = 0.5`)
  3. Adjust all theoretical intensities: `I_adj = I_orig × exp(f_upper + f_lower)`; uncertainty propagates in quadrature
  4. Reset all acceptance decisions and re-run weeding with adjusted intensities
- **Convergence:** stops when no transition changes acceptance state, or after `max_iterations`
- **Oscillation blacklist:** after iteration 2, transitions that flip acceptance state are blacklisted and fall through to their Step 3 default on all subsequent passes

**Per-line 3-step decision process (`weed_assignments_line`):**

Each line's candidate transitions are sorted by decreasing `calc_intensity`, then by Ritz σ.

- **Step 1 — Clear-cut decisions:** New transitions with Ritz mismatch > 4σ are rejected. A relative-intensity filter removes new transitions contributing < 10% of total `calc_intensity` (< 4% for legacy). Surviving candidates are tested by z-score against I_obs: asymmetric thresholds apply (looser for resonance lines and for cases where `I_calc > I_obs`). Decisions are labeled "Step1" in the `notes2` field.

- **Step 2 — Grouping and late arbitration:** Undecided candidates are tested for pair/triple acceptance: groups pass if their joint center-of-gravity matches the observed wavenumber and their combined effective spread (Doppler width convolved with measurement uncertainty) is consistent with the line profile. For lines flagged `h` or `w`, broadening thresholds are relaxed. Remaining undecided new candidates with `F` or `R` in `notes1` are rejected.

- **Step 3 — Default decisions:** Undecided **new** candidates: rejected ("no solid evidence for acceptance"). Undecided **legacy** candidates: accepted ("no solid evidence for rejection").

### Energy Level Optimization

After weeding, `optimize_levels()` refines each non-ground energy level by a weighted mean of all implied energies from accepted transitions (both upper-from-lower and lower-from-upper). Weights are `1/σ²` for unblended lines; for blended lines the weight is shared proportionally to `calc_intensity`. The iteration runs until max energy change < 0.001 cm⁻¹. The `weight` column in the output is the LOPT-compatible weight for each accepted transition.

### Grading Scheme

Matches are graded on two dimensions:

-   **Tier (Wavenumber Agreement):**
    -   `2`: Wavenumber residual ≤ 2.0σ (uncertainty)
    -   `3`: Wavenumber residual ≤ 3.0σ
    -   `4`: Wavenumber residual ≤ 4.0σ
    -   `5`: Wavenumber residual > 5.0σ
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
