"""
Configuration for the Pr III Table Extraction Pipeline.
Loads settings from .env file and provides defaults.
"""

import os
import sys
from pathlib import Path
from dotenv import load_dotenv

# Load environment variables from .env
load_dotenv()

# --- Paths ---
BASE_DIR = Path(__file__).parent
PDF_PATH = BASE_DIR / "jresv73An3p333_A1b_TableX.pdf"
NIST_ASD_PATH = BASE_DIR / "Pr3_lev_ASD512.xlsx"
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

# --- Validation Thresholds ---
WAVENUMBER_TOLERANCE_CM1 = 0.2  # Tolerance for wavelength↔wavenumber rounding check
NIST_RITZ_TOLERANCE_CM1 = 1.5   # Tolerance for Ritz wavenumber vs observed
