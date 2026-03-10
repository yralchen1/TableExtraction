# TableExtraction

An AI-powered pipeline for extracting spectral line tables from scanned scientific PDFs using **Google Gemini** multimodal vision and **LangChain/LangGraph**.

Originally built for digitizing tables from **Sugar's Pr III** spectroscopy papers ([1969](https://doi.org/10.6028/jres.073A.029) & [1974](https://doi.org/10.6028/jres.078A.036)), the pipeline converts PDF pages to images, sends them to Gemini for OCR/extraction, validates against physics formulas and [NIST ASD](https://physics.nist.gov/asd) reference energy levels, corrects page-level misalignments, and outputs structured Excel/CSV files.

---

## Features

- **Multimodal OCR** — Extracts structured data from scanned two-column table layouts using Gemini vision models (Pro & Flash)
- **Physics-aware validation** — Validates air↔vacuum wavelength/wavenumber via Edlén's 1966 formula
- **NIST ASD cross-reference** — Looks up energy levels against NIST Atomic Spectra Database with auto-correction for typos, parity flips, and J-value mismatches
- **Selection rule checks** — Enforces E1 transition rules (parity change, ΔJ ≤ 1)
- **Misalignment correction** — Detects and fixes OCR-induced row shifts using Ritz wavenumber matching
- **Line Classification** — Automated matching of observed lines to energy levels with multi-tier grading (A1-D6)
- **Robust JSON parsing** — 6-step fallback parser for malformed LLM output
- **Configurable** — Year-specific prompts, column contexts, layout definitions, and environment-driven settings

---

## Prerequisites

- **Python 3.10+**
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

### Output Columns

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
│   └── requirements.txt                   # Python dependencies
└── LineClass/                             # Spectral Line Classification
    ├── classify_lines.py                  # CLI entry point (classification)
    ├── models.py                          # Data models (EnergyLevel, SpectralLine, Transition)
    ├── Pr3_lines.xlsx                     # Observed lines input
    ├── Pr3_tp_Dream.xlsm                  # DREAM calculated transitions
    └── line_classifications.xlsx          # Classification output
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

---

## Line Classification Pipeline (`LineClass`)

This module classifies observed spectral lines by matching them against all possible transitions between known energy levels.

### 6-Step Workflow

1.  **Read Energy Levels** — Loads data from `Pr3_lev_Wyart_1999.xlsm`
2.  **Read Observed Lines** — Loads data from `Pr3_lines.xlsx`
3.  **Read Theoretical Transitions** — Loads DREAM calculated intensities/CFs from `Pr3_tp_Dream.xlsm`
4.  **Generate Transitions** — Calculates all possible transitions satisfying selection rules (ΔJ ≤ 1, parity change)
5.  **Match & Grade** — Matches observed lines to transitions and assigns a grade (A1 to D6) based on wavenumber residual and intensity consistency
6.  **Resolve Conflicts** — Handles transitions assigned to multiple lines, selecting the best match and tagging revisions

### Grading Scheme

Matches are graded on two dimensions:

-   **Tier (Wavenumber Agreement):**
    -   `A`: Wavenumber residual ≤ 1.0σ (uncertainty)
    -   `B`: Wavenumber residual ≤ 2.0σ
    -   `C`: Wavenumber residual ≤ 3.0σ
    -   `D`: Wavenumber residual > 3.0σ
-   **Subgrade (Intensity Consistency):**
    -   `1` to `3`: Good to poor agreement with DREAM calculated intensity (CF ≥ 0.1)
    -   `4` to `5`: Agreement with low-CF transitions
    -   `6`: No theoretical intensity available

**Suffixes:**
-   `N`: Newly assigned (not in original source)
-   `F`: Conflicting assignment (multiple lines match one transition)
-   `R`: Revised (original classification was moved to a better row)

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

---

## License

This project is for academic and research purposes.
