# Pr III Table Extraction Pipeline (v3)

Extract **Table 1** from the scanned Sugar Pr III paper ([jresv78An5p555_A1b.pdf](https://doi.org/10.6028/jres.078A.032)) using **Google Gemini** multimodal vision via a **LangChain/LangGraph** pipeline.

The pipeline converts PDF pages to images, sends them to Gemini for OCR/extraction, validates against physics formulas and NIST ASD energy levels, corrects page-level misalignments, and outputs a 14-column Excel/CSV sorted by wavelength (descending).

---

## Prerequisites

- **Python 3.10+**
- **poppler** (for PDF→image conversion):
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

### 1. Create & activate a virtual environment

```bash
cd "/Users/themanaspandey/Documents/Table Extraction"
python3 -m venv venv
source venv/bin/activate        # macOS/Linux
# venv\Scripts\activate         # Windows
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

### 3. Configure your API key

```bash
cp .env.example .env
nano .env    # add your Gemini API key
```

Your `.env` file should contain:
```
GOOGLE_API_KEY=AIza...your_key_here
```

### 4. (Optional) Choose a Gemini model

Edit `.env`:
```
GEMINI_MODEL=gemini-2.5-pro            # Best accuracy, higher cost / stricter quota
# GEMINI_MODEL=gemini-2.5-flash        # Default: fast and cheap
# GEMINI_MODEL=gemini-2.0-flash        # Alternative
# GEMINI_MODEL=gemini-3-flash-preview  # High accuracy, moderate cost
```

Or pass as CLI flag: `--model gemini-2.5-pro`

> **Note:** `gemini-2.5-pro` produces significantly better table alignment than `gemini-2.5-flash`, especially for the two-column layout. The flash model tends to misalign classification columns across rows.

---

## Usage

```bash
source venv/bin/activate

# Full pipeline (all pages):
python extract_table_Sugar1974.py

# Specific page range:
python extract_table_Sugar1974.py --pages 3 39

# Single page (for testing):
python extract_table_Sugar1974.py --single-page 5

# Use gemini-2.5-pro:
python extract_table_Sugar1974.py --pages 3 39 --model gemini-2.5-pro

# Reuse existing page images (skip PDF conversion):
python extract_table_Sugar1974.py --skip-images --model gemini-2.5-pro

# Lower DPI for faster processing:
python extract_table_Sugar1974.py --dpi 200 --pages 3 39

# Dry run (show config, no processing):
python extract_table_Sugar1974.py --dry-run --pages 3 39

# Convert PDF to images only:
python pdf_to_images.py --pages 3 39
```

---

## Output

Results in `output/`:

| File | Description |
|------|-------------|
| `Pr_III_Table1_extracted.csv` | CSV with all 14 columns |
| `Pr_III_Table1_extracted.xlsx` | Excel with auto-sized columns |

> **CSV J-values:** Half-integer J values (e.g., `5/2`) are written as `="5/2"` in CSV to prevent Excel from misinterpreting them as dates or fractions. The Excel file stores them natively as text.

### 14 Output Columns

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

## Validation Pipeline (v3)

### Validation Order

The pipeline validates each row in a specific order to avoid incorrect corrections:

1. **NIST ASD lookup first** — looks up both levels with auto-correction
2. **Selection rule checks** — parity change + ΔJ ≤ 1 — **only when BOTH levels are found** (v3 fix)
3. **Wavelength/wavenumber check** — only corrects if NIST confirmed OK
4. **Page-level misalignment correction** (v3) — relocates misplaced classifications

> **Key insight:** If the NIST lookup fails (levels not found or Ritz mismatch), the row's classification is likely misaligned by OCR, meaning the wavelength/wavenumber is probably correct as-is. The pipeline marks it as "suspect misalignment" without modifying the wavelength, then attempts to fix the misalignment in step 4.

### 1. NIST ASD Auto-Correction

When a level isn't found exactly, the pipeline tries (in order):
- **±1 in integer part** with exact parity + J
- **Flip parity** (`e` ↔ `o`) — catches ° misread as digit 0
- **Different J from NIST** — for levels where NIST revised the J value
- **Single-digit typo** in the integer part
- **Adjacent-digit swap** in the integer part

### 2. Selection Rules (E1 Transitions)

- **Parity must change**: lower and upper levels must have different parity
- **ΔJ ≤ 1**: `|J_upper - J_lower| ≤ 1`
- **Only checked when BOTH levels found in NIST** (v3 fix — previously could fire with one corrected + one missing level, producing false positives)

If selection rules fail after correction, the line is flagged `RITZ_OK_SELECTION_FAIL`.

### 3. Wavelength ↔ Wavenumber (Edlen's 1966 Formula)

Converts air wavelength to vacuum wavenumber using Edlen's formula. If mismatch detected:
- **NIST OK** → tries typo correction (single digit, adjacent swap) — correction is safe
- **NIST failed** → marks as "suspect misalignment" — NO correction applied

### 4. Page-Level Misalignment Correction (v3)

When a page has **> 2 error rows** (NIST failures or Ritz mismatches), the pipeline hypothesizes that OCR shifted some classifications to wrong rows (typically shifted upward by 1–10 rows).

For each error row with NIST-confirmed levels:
1. Compute **Ritz wavenumber** = upper_exact − lower_exact
2. Search all rows on the **same page** for the one whose observed wavenumber is closest
3. If that row is **unclassified** → move the classification there ✅
4. If that row is **already classified** → flag as unresolvable ❌

After relocation, the target row's wavelength/wavenumber is re-validated.

### 5. Robust JSON Parsing

Handles malformed LLM output with 6-step fallback:
1. Direct parse
2. Strip markdown code fences
3. Fix trailing commas
4. Fix unquoted property names
5. Extract JSON array from surrounding text
6. Parse individual `{...}` objects

On total parse failure, saves raw response to `page_images/{page}_raw_response.txt` for debugging.

### 6. Rate-Limit Handling

The pipeline catches 429 (rate limit) errors from the Gemini API. If you hit quota limits:
- Wait for the retry period shown in the error message
- Consider using `gemini-3-flash-preview` (much better OCR than with gemini-2.5-flash)
- Process fewer pages at a time with `--single-page` or `--pages`

---

## Project Structure

```
Table Extraction/
├── extract_table_Sugar1974.py     # CLI entry point for Table 1 of Sugar 1974
├── config_Sugar1974.py            # Configuration & env loading
├── pdf_to_images_Sugar1974.py     # PDF → PNG conversion
├── extraction_prompt_Sugar1974.py # Gemini extraction prompt (v3)
├── tools_Sugar1974.py             # Validation (Edlen, selection rules, NIST)
├── graph_Sugar1974.py             # LangGraph pipeline (v3)
├── extract_table_Sugar1969.py     # CLI entry point for Table 1 of Sugar 1969
├── config_Sugar1969.py            # Configuration & env loading
├── pdf_to_images_Sugar1969.py     # PDF → PNG conversion
├── extraction_prompt_Sugar1969.py # Gemini extraction prompt (v1)
├── tools_Sugar1969.py             # Validation (Edlen, selection rules, NIST) (v1)
├── graph_Sugar1969.py             # LangGraph pipeline (v1)
├── requirements.txt               # Python dependencies
├── .env.example                   # API key template
├── .env                           # Your API key (create from template)
├── jresv78An5p555_A1b.pdf         # Source PDF for Sugar 1974
├── jresv73An3p333_A1b.pdf         # Complete PDF for Sugar 1969
├── jresv73An3p333_A1b_TableX.pdf  # Source PDF for Table X of Sugar 1969
├── Pr3_lev_ASD512.xlsx            # NIST ASD reference data
├── page_images/                   # Generated page images
└── output/                        # CSV/Excel output
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
