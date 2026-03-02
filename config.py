"""
Configuration for the Pr III Table Extraction Pipeline.
Loads settings from .env file and provides defaults.
"""

import os
import sys
from pathlib import Path
from typing import Dict, Any
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

# --- Year-Specific Configs ---
CONFIG_1969 = {
    "pdf_path": BASE_DIR / "jresv73An3p333_A1b_TableX.pdf",
    "nist_asd_path": BASE_DIR / "Pr3_lev_ASD512.xlsx",
    "wavenumber_tolerance_cm1": 0.2,
    "nist_ritz_tolerance_cm1": 1.5,
    "wavelength_type": "vac",
    "csv_name": "Pr_III_TableX_extracted.csv",
    "excel_name": "Pr_III_TableX_extracted.xlsx",
    "excel_sheet_name": "Table X",
}

CONFIG_1974 = {
    "pdf_path": BASE_DIR / "jresv78An5p555_A1b.pdf",
    "nist_asd_path": BASE_DIR / "Pr3_lev_ASD512.xlsx",
    "wavenumber_tolerance_cm1": 0.01,
    "nist_ritz_tolerance_cm1": 1.5,
    "wavelength_type": "air",
    "csv_name": "Pr_III_Table1_extracted.csv",
    "excel_name": "Pr_III_Table1_extracted.xlsx",
    "excel_sheet_name": "Table 1",
}

CONFIGS = {
    "1969": CONFIG_1969,
    "1974": CONFIG_1974,
}

_current_year = "1969"

def set_year(year: str):
    global _current_year
    if str(year) not in CONFIGS:
        raise ValueError(f"Invalid year: {year}. Must be '1969' or '1974'.")
    _current_year = str(year)

def get_config() -> Dict[str, Any]:
    return CONFIGS[_current_year]
