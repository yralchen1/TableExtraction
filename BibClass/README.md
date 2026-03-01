# Benchmarking PDFs — Automated Keyword Extraction for Atomic Spectroscopy

This project automates the extraction of **`keywords_el`** (Atomic Energy Levels & Spectra bibliographic keywords) from scientific papers in PDF format. It uses **Google Gemini** (via LangChain) to read each paper and produce BibTeX-format keyword annotations following the NIST ASD bibliographic database conventions.

---

## Table of Contents

1. [Overview](#overview)
2. [Folder Structure](#folder-structure)
3. [Prerequisites](#prerequisites)
4. [Getting a Gemini API Key](#getting-a-gemini-api-key)
5. [Setting Up the Environment](#setting-up-the-environment)
6. [Running the PDF Processing Script](#running-the-pdf-processing-script)
7. [Understanding the Output](#understanding-the-output)
8. [Changing the Prompt](#changing-the-prompt)
9. [Changing the Gemini Model](#changing-the-gemini-model)
10. [Other Scripts](#other-scripts)
11. [Troubleshooting](#troubleshooting)

---

## Overview

The core workflow is:

1. **Organize** raw PDFs (named like `Author_el_ID_YEAR.pdf`) into structured folders using `organize_pdfs.py`.
2. **Distribute** existing BibTeX catalogue entries into each folder using `catalogue_bibtex.py`.
3. **Process** each paper through Google Gemini to automatically extract atomic-physics keywords using `process_pdfs_langchain.py`.

The AI-generated keywords are saved alongside each paper as `bibtex_AI_Generated.txt`, which can then be compared with the human-catalogued `bibtex_catalogued.txt` for benchmarking.

---

## Folder Structure

```
Benchmarking PDFs/
│
├── .env                          # Your Gemini API key (DO NOT SHARE)
├── .env.example                  # Template for the .env file
├── requirements.txt              # Python dependencies
├── process_pdfs_langchain.py     # Main AI processing script
├── organize_pdfs.py              # Organizes raw PDFs into folders
├── catalogue_bibtex.py           # Distributes BibTeX entries to folders
├── bibtex_Sr1_1976-2025.txt      # Master BibTeX catalogue file
├── README.md                     # This file
│
├── 2025_Cheung_el_23861/         # Example paper folder
│   ├── main_article.pdf          # The primary PDF paper
│   ├── bibtex_catalogued.txt     # Human-written BibTeX keywords (ground truth)
│   ├── bibtex_AI_Generated.txt   # AI-generated BibTeX keywords (Gemini output)
│   └── suppl/                    # Supplementary materials (if any)
│
├── 2024_Bothwell_el_22893/       # Another paper folder
│   ├── main_article.pdf
│   ├── bibtex_catalogued.txt
│   ├── bibtex_AI_Generated.txt
│   └── suppl/
│
└── ... (90+ paper folders)
```

---

## Prerequisites

- **Python 3.10+** (tested with Python 3.14)
- **pip** (Python package manager)
- **A Google Gemini API key** (see [Getting a Gemini API Key](#getting-a-gemini-api-key))
- **Operating System:** macOS, Linux, or Windows 10/11

---

## Getting a Gemini API Key

You need a Google Gemini API key to run the processing script. Here's how to get one:

1. **Go to Google AI Studio**: [https://aistudio.google.com/apikey](https://aistudio.google.com/apikey)
2. **Sign in** with your Google account.
3. **Click "Create API Key"**.
4. **Select a Google Cloud project** (or create a new one).
5. **Copy the generated API key** — you'll need it in the next step.

> **⚠️ Important:** Keep your API key secret. Never commit it to version control or share it publicly.

### API Key Pricing

- The **Gemini Flash Lite** model (used by default) is part of Google's **free tier** for low-volume usage.
- For high-volume processing (100+ papers), check the [Gemini API pricing page](https://ai.google.dev/pricing) for current limits and costs.

---

## Setting Up the Environment

### Step 1 — Clone/Navigate to the project

**macOS / Linux:**
```bash
cd "/path/to/BibClass project folder"
```

**Windows (Command Prompt):**
```cmd
cd "C:\path\to\BibClass project folder"
```

**Windows (PowerShell):**
```powershell
cd "C:\path\to\BibClass project folder"
```

### Step 2 — Create a virtual environment (recommended)

**macOS / Linux:**
```bash
python3 -m venv .venv
source .venv/bin/activate
```

**Windows (Command Prompt):**
```cmd
python -m venv .venv
.venv\Scripts\activate.bat
```

**Windows (PowerShell):**
```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
```

> **Note (Windows PowerShell):** If you get an execution policy error, run:
> ```powershell
> Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
> ```

### Step 3 — Install dependencies

```bash
pip install -r requirements.txt
```

> This command is the same on all platforms once the virtual environment is activated.

This installs:
| Package | Purpose |
|---|---|
| `langchain` | LLM orchestration framework |
| `langchain-google-genai` | Google Gemini integration for LangChain |
| `pypdf` | PDF text extraction |
| `python-dotenv` | Load environment variables from `.env` |

### Step 4 — Set up your API key

Create a `.env` file in the `Benchmarking PDFs/` directory:

**macOS / Linux:**
```bash
cp .env.example .env
```

**Windows (Command Prompt):**
```cmd
copy .env.example .env
```

**Windows (PowerShell):**
```powershell
Copy-Item .env.example .env
```

Then open `.env` in any text editor and replace the placeholder with your actual API key:

```env
# Gemini API Key
GOOGLE_API_KEY=your_actual_gemini_api_key_here
```

> **Tip:** You can also set the key as a shell environment variable instead:
>
> **macOS / Linux:**
> ```bash
> export GOOGLE_API_KEY="your_actual_gemini_api_key_here"
> ```
>
> **Windows (Command Prompt):**
> ```cmd
> set GOOGLE_API_KEY=your_actual_gemini_api_key_here
> ```
>
> **Windows (PowerShell):**
> ```powershell
> $env:GOOGLE_API_KEY = "your_actual_gemini_api_key_here"
> ```

---

## Running the PDF Processing Script

The main script is **`process_pdfs_langchain.py`**. It reads each paper's PDF, sends the extracted text to Google Gemini, and saves the AI-generated keywords.

### Process All Papers

```bash
python process_pdfs_langchain_with_suppl.py
```

This will:
1. Scan all subfolders for `main_article.pdf` files (the initial distribution includes only one such subfolder named `tmp`). If such folder contains a subfolder `suppl`, its contents are treated as supplementary files for `main_article.pdf` and will be processed together with it.
2. Extract text from each PDF using `pypdf`.
3. Send the text (up to ~50,000 characters) to Gemini with the keyword-extraction prompt.
4. Save the response as `bibtex_AI_Generated.txt` inside each paper's folder.
5. Print a progress summary showing successes and errors.

**Example terminal output:**
```
🚀 Initializing Google Gemini...
📁 Found 1 paper(s) to process

[1/1] Processing: tmp
  📄 Extracting PDF text...
  📝 Extracted 42,317 characters
  🤖 Sending to Gemini...
  ✅ Saved to: bibtex_AI_Generated.txt

...

==================================================
Summary:
  ✅ Successful: 1
  ❌ Errors: 0
  📁 Total: 1
```

### Process a Single Paper

To process only one specific paper folder:

```bash
python process_pdfs_langchain_with_suppl.py --single tmp
```

The `--single` flag accepts the **folder name** (not the full path). The folder must exist inside the `Benchmarking PDFs/` directory and contain a `main_article.pdf`.

### Dry Run (Preview Without Processing)

To see which papers would be processed without actually calling the API:

```bash
python process_pdfs_langchain.py --dry-run
```

**Output:**
```
📁 Found 1 paper(s) to process

DRY RUN - Would process:
  • tmp
  ...
```

### Combine Flags

You can combine `--single` and `--dry-run`:

```bash
python process_pdfs_langchain.py --single tmp --dry-run
```

---

## Understanding the Output

Each processed paper gets a `bibtex_AI_Generated.txt` file saved in its folder. The file contains the extracted keywords in BibTeX format:

**Example output** (`tmp/bibtex_AI_Generated.txt`):
```
keywords_el={Ca I; 43Ca I: Hfs: E
Ca I: IS: E
Ca I; 43Ca I: CL: E
Ca I; 43Ca I: W: E},
keywords_tp={},
keywords_lb={}
```

### How to Read the Keywords

Each line inside `keywords_el={...}` follows one of these formats:

- **General Interest:** `GENINT: [Code]: [Method]`
  - Example: `GENINT: 1.8: T` = Atomic codes (Theory)
  - Example: `GENINT: 1.13: E` = Atomic Clocks (Experiment)

- **Element-Specific:** `[Spectrum]: [Subject Code]: [Method]`
  - Example: `Sr I: SE: T` = Strontium I, Stark Effect, Theory
  - Example: `Sr I; 87Sr I: EL: E` = Strontium I & Strontium-87 I, Energy Levels, Experiment

### Method Types
| Code | Meaning |
|---|---|
| `E` | Experimental data |
| `T` | Theoretical calculations |
| `O` | Other / semi-empirical |

### Common Subject Codes
| Code | Full Name |
|---|---|
| `EL` | Energy Levels |
| `W` | Wavelengths / Frequencies |
| `CL` | Classified Lines |
| `TE` | Theoretical Energies |
| `AT` | Ab Initio Theory |
| `SE` | Stark Effect / Polarizability / BBR shifts |
| `ZE` | Zeeman Effect / g-factors |
| `Hfs` | Hyperfine Structure |
| `IS` | Isotope Shifts |
| `IP` | Ionization Potential |
| `QF` | QED / Lamb Shifts |

---

## Changing the Prompt

The extraction prompt is defined as the `SYSTEM_PROMPT` variable at the top of `process_pdfs_langchain_with_suppl.py` (starting at **line 36**). This is the core instruction that tells Gemini what to extract and how to format the output.

### How to Modify the Prompt

1. **Open the script** in your editor:

   ```bash
   open process_pdfs_langchain_with_suppl.py
   ```

2. **Find the `SYSTEM_PROMPT` variable** — it starts at line 36 with a triple-quoted string:

   ```python
   SYSTEM_PROMPT = """**Role:** Your sole task is to extract and format `keywords_el`, `keywords_tp`, and `keywords_lb` from the attached PDF and its supplementary files (if any).
   ...
   """
   ```

3. **Edit the prompt text** between the triple quotes (`"""`). The prompt is organized into these sections:

   | Section | What It Controls | Lines |
   |---|---|---|
   | **Role & Output Format** | Tells Gemini what role to play and the exact output format | 36–44 |
   | **Critical Rules** | Core extraction rules (be conservative, one method per line, etc.) | 46–76 |
   | **Real Examples** | Worked examples showing correct keyword formatting | 78–89 |
   | **SPECS** | Full specification document for the keyword system | 91–286 |

4. **Save the file** — the next time you run `process_pdfs_langchain.py`, it will use your updated prompt.

### Common Prompt Modifications

#### Add a New Subject Code

If you need to add a new keyword category, add it to the **Subject Codes** list in the prompt (around line 53):

```python
# Add after the existing codes:
   - `NEW` = New Code Description
```

And add it to the SPECS table section as well.

#### Change the Output Format

If you want JSON output instead of BibTeX format, modify the **Output** section (lines 39–44):

```python
# Change from:
**Output:** Output ONLY the keywords_el field in exact BibTeX format.

# To something like:
**Output:** Output a JSON object with a single key "keywords_el" containing an array of keyword strings.
```

> **⚠️ Warning:** If you change the output format, you may also need to update the `process_paper()` function (line 306) to handle the new format correctly.

#### Adjust the Extraction Strictness

The prompt currently says "Be CONSERVATIVE". To make it more or less strict:

```python
# More aggressive (may produce false positives):
1. **Be THOROUGH** - Assign keywords for data that is explicitly or implicitly present in the paper.

# More conservative (may miss valid keywords):
1. **Be EXTREMELY CONSERVATIVE** - Only assign keywords when you are 100% certain the data is explicitly reported in numerical tables.
```

#### Add More Examples

Adding more examples improves Gemini's accuracy. Add them in the **REAL EXAMPLES** section (after line 78):

```python
Paper about isotope shift measurements:
keywords_el={Ca I; Ca II: IS: E
41Ca I; 43Ca I; 45Ca I: Hfs: E}
```

### Where the Prompt Is Used

The prompt flows through the code as follows:

```
SYSTEM_PROMPT (line 36)
    ↓
process_paper() function (line 306)
    ↓
SystemMessage(content=SYSTEM_PROMPT)  →  sent to Gemini as the system instruction
HumanMessage(content=pdf_text)        →  the paper text is sent as the user message
    ↓
Gemini returns the keywords_el string
    ↓
Saved to bibtex_AI_Generated.txt
```

---

## Changing the Gemini Model

The Gemini model is configured on **line 358** of `process_pdfs_langchain.py`:

```python
llm = ChatGoogleGenerativeAI(
    model="gemini-flash-lite-latest",    # The model to use
    google_api_key=GOOGLE_API_KEY,
    temperature=0.1                      # Low = more deterministic
)
```

### Available Models

| Model Name | Speed | Quality | Cost |
|---|---|---|---|
| `gemini-flash-lite-latest` | ⚡ Fastest | Good | Lowest |
| `gemini-2.0-flash` | Fast | Better | Low |
| `gemini-2.5-pro-preview-06-05` | Slower | Better | Higher |
| `gemini-3-flash-preview` | Slower | Best | Higher |

To switch models, change the `model=` parameter:

```python
llm = ChatGoogleGenerativeAI(
    model="gemini-2.0-flash",      # Switch to a more capable model
    google_api_key=GOOGLE_API_KEY,
    temperature=0.1
)
```

### Adjusting Temperature

- `temperature=0.0` — Most deterministic, always picks the highest-probability token.
- `temperature=0.1` — **(Current)** Slightly creative but mostly consistent.
- `temperature=0.5` — More varied outputs; not recommended for structured extraction.

---

## Other Scripts

### `organize_pdfs.py` — Organize Raw PDFs into Folders

Takes flat PDF files named `Author_el_ID_YEAR.pdf` and organizes them into structured folders.

```bash
# Preview what would happen (dry run):
python organize_pdfs.py --dry-run

# Actually organize the PDFs:
python organize_pdfs.py

# Organize PDFs in a different directory:
python organize_pdfs.py /path/to/pdf/directory
```

**Input:** `Cheung_el_23861_2025.pdf`
**Output:** `2025_Cheung_el_23861/main_article.pdf` + `suppl/` subdirectory

### `catalogue_bibtex.py` — Distribute BibTeX Entries

Reads the master BibTeX file (`bibtex_Sr1_1976-2025.txt`) and places the relevant entry into each paper's folder as `bibtex_catalogued.txt`.

```bash
python catalogue_bibtex.py
```

---

## Troubleshooting

### ❌ "GOOGLE_API_KEY not found!"

Your `.env` file is missing or doesn't contain the key.

**macOS / Linux:**
```bash
cat .env
# Should contain: GOOGLE_API_KEY=AIza...your_key_here
```

**Windows:**
```cmd
type .env
```

### ❌ "ModuleNotFoundError: No module named 'langchain'"

Dependencies are not installed. Make sure you're in the virtual environment:

**macOS / Linux:**
```bash
source .venv/bin/activate
pip install -r requirements.txt
```

**Windows (Command Prompt):**
```cmd
.venv\Scripts\activate.bat
pip install -r requirements.txt
```

**Windows (PowerShell):**
```powershell
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### ❌ "Failed to extract text from PDF"

The PDF may be image-based (scanned) rather than text-based. `pypdf` can only extract text from text-based PDFs. For scanned PDFs, you would need OCR (e.g., `pytesseract`).

### ❌ API Rate Limit Errors

If processing many papers at once, you may hit Gemini's rate limits. Try:
- Using `--single` to process one paper at a time.
- Waiting a few minutes between batches.
- Upgrading to a paid API tier.

### ❌ Truncated Output

The script limits input to ~50,000 characters per paper (line 310). Very long papers may have their later sections cut off. To increase this limit:

```python
# In process_paper(), line 310, change:
HumanMessage(content=f"...{pdf_text[:50000]}")

# To a larger value:
HumanMessage(content=f"...{pdf_text[:100000]}")
```

> **Note:** Larger inputs consume more tokens and may increase API costs.

---

## Quick Reference

**macOS / Linux:**
```bash
# Setup (one-time):
cd "Benchmarking PDFs"
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cp .env.example .env
# Edit .env with your API key

# Run:
python process_pdfs_langchain.py              # Process all papers
python process_pdfs_langchain.py --single X   # Process folder X only
python process_pdfs_langchain.py --dry-run    # Preview only
```

**Windows (Command Prompt):**
```cmd
:: Setup (one-time):
cd "Benchmarking PDFs"
python -m venv .venv
.venv\Scripts\activate.bat
pip install -r requirements.txt
copy .env.example .env
:: Edit .env with your API key

:: Run:
python process_pdfs_langchain.py
python process_pdfs_langchain.py --single X
python process_pdfs_langchain.py --dry-run
```

**Windows (PowerShell):**
```powershell
# Setup (one-time):
cd "Benchmarking PDFs"
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
Copy-Item .env.example .env
# Edit .env with your API key

# Run:
python process_pdfs_langchain.py
python process_pdfs_langchain.py --single X
python process_pdfs_langchain.py --dry-run
```
