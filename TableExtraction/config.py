"""
Configuration for the Table Extraction Pipeline.
Loads settings from .env file and provides defaults.
"""

import os
import sys
import datetime
from pathlib import Path
from dotenv import load_dotenv

# Load environment variables from .env
load_dotenv()

# --- Paths ---
BASE_DIR = Path(__file__).parent
PAGE_IMAGES_DIR = BASE_DIR / "page_images"
OUTPUT_DIR = BASE_DIR / "output"

# --- API Key ---
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY", "")

def validate_api_key():
    """Check that the API key is set."""
    if not GOOGLE_API_KEY or GOOGLE_API_KEY == "your_api_key_here":
        print("❌ Error: GOOGLE_API_KEY not configured!")
        print("   1. Copy .env.example to .env:")
        print("      cp .env.example .env")
        print("   2. Edit .env and add your Gemini API key")
        print("   Get a key at: https://aistudio.google.com/apikey")
        sys.exit(1)

# --- Model ---
GEMINI_MODEL = os.getenv("GEMINI_MODEL", "gemini-3-flash-preview")

# --- PDF Conversion ---
PDF_DPI = int(os.getenv("PDF_DPI", "600"))

# --- Table Config ---
PDF_PATH = Path(os.getenv("PDF_PATH", str(BASE_DIR / "jresv73An3p333_A1b_TableX.pdf")))
REF_LEVELS_PATH = Path(os.getenv("REF_LEVELS_PATH", str(BASE_DIR / "Pr3_lev_ASD512.xlsx")))

WAVENUMBER_TOLERANCE_CM1 = float(os.getenv("WAVENUMBER_TOLERANCE_CM1", "0.2"))
REF_RITZ_TOLERANCE_CM1 = float(os.getenv("REF_RITZ_TOLERANCE_CM1", "1.5"))
WAVELENGTH_TYPE = os.getenv("WAVELENGTH_TYPE", "vac")

LAYOUT_PROMPT_PATH = Path(os.getenv("LAYOUT_PROMPT_PATH", str(BASE_DIR / "prompts/layout_definitions/split_two_columns.txt")))
COLUMN_CONTEXT_PATH = Path(os.getenv("COLUMN_CONTEXT_PATH", str(BASE_DIR / "prompts/column_contexts/sugar_1969.json")))

CSV_NAME = os.getenv("CSV_NAME", "Table_extracted.csv")
EXCEL_NAME = os.getenv("EXCEL_NAME", "Table_extracted.xlsx")
EXCEL_SHEET_NAME = os.getenv("EXCEL_SHEET_NAME", "Table 1")

# --- New Options ---
AIR_VAC_CONVERSION_METHOD = os.getenv("AIR_VAC_CONVERSION_METHOD", "Edlen1966")
current_year = datetime.datetime.now().year
PUBLICATION_YEAR = int(os.getenv("PUBLICATION_YEAR", str(current_year)))

def get_config() -> dict:
    """Returns the configuration as a dictionary."""
    return {
        "pdf_path": PDF_PATH,
        "ref_levels_path": REF_LEVELS_PATH,
        "wavenumber_tolerance_cm1": WAVENUMBER_TOLERANCE_CM1,
        "ref_ritz_tolerance_cm1": REF_RITZ_TOLERANCE_CM1,
        "wavelength_type": WAVELENGTH_TYPE,
        "csv_name": CSV_NAME,
        "excel_name": EXCEL_NAME,
        "excel_sheet_name": EXCEL_SHEET_NAME,
        "air_vac_conversion_method": AIR_VAC_CONVERSION_METHOD,
        "publication_year": PUBLICATION_YEAR,
        "layout_prompt_path": LAYOUT_PROMPT_PATH,
        "column_context_path": COLUMN_CONTEXT_PATH,
    }
