#!/usr/bin/env python3
"""
Simple script to extract a single PDF page directly to a CSV file.
It skips all validation and misalignment correction – just raw OCR.
"""

import sys
import base64
from pathlib import Path
import pandas as pd
from pdf2image import convert_from_path
import google.generativeai as genai
from dotenv import load_dotenv
import os
import json
import re

# Load environment variables
load_dotenv()
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")
if not GOOGLE_API_KEY:
    print("Error: GOOGLE_API_KEY environment variable not set.")
    sys.exit(1)

genai.configure(api_key=GOOGLE_API_KEY)

# Use the same prompt structure from the main pipeline
SYSTEM_PROMPT = """You are an expert spectroscopist digitizing a scanned table of spectral lines of doubly ionized praseodymium (Pr III) from a 1970s scientific paper by Sugar.

Your task: Extract EVERY row from the table image and return structured JSON.

TABLE LAYOUT — READ CAREFULLY:
Each page is split into TWO completely distinct vertical halves (a left section and a right section).
You MUST treat the left section and right section as two separate tables that happen to be printed side-by-side.
Each section independently contains FOUR main columns in this exact order:
  Wavelength (Å)  |  Intensity  |  Wavenumber (cm⁻¹)  |  Classification

The Classification column shows transitions like: "12847₉/₂ — 22176₁₁/₂"
- The left part (before the dash) is the LOWER energy level
- The right part (after the dash) is the UPPER energy level

Column details:
1. Wavelength (Å) — Usually a decimal number with EXACTLY 3 decimal places.
   ⚠️ IMPORTANT: If a spectral line has multiple classifications, the Wavelength on the continuation row may be a ditto mark (like `"` or `〃` or `''`) or completely blank. Extract the ditto mark EXACTLY as printed into the JSON string.
2. Intensity — An integer, optionally followed by character symbols. If the intensity is exactly "0", extract it as the integer 0.
3. Wavenumber (cm⁻¹) — A decimal number with EXACTLY 2 decimal places.
4. Classification — Wait for the integer, parity (° for odd, nothing for even), and J value.

EXTRACTION RULES:
1. You MUST process the LEFT section of the page entirely (from top to bottom), and ONLY THEN process the RIGHT section (from top to bottom). Do not read horizontally across the middle gap between the two sections.
2. Read each row strictly as a complete horizontal line within its section.
3. If a row has no classification, set all 6 classification fields to null.

Output format:
Return a JSON array where each element is an object:
[
  {
    "wavelength": 10716.061,
    "intensity": 50,
    "line_character": "h",
    "wavenumber": 9329.40,
    "lower_level_int": 12847,
    "lower_parity": "e",
    "lower_j": "9/2",
    "upper_level_int": 22176,
    "upper_parity": "o",
    "upper_j": "11/2"
  }
]
Return ONLY the JSON array. No markdown, no explanations."""

USER_PROMPT = "Extract all rows from this spectral table image into structured JSON as instructed. Read each row strictly horizontally — do not misalign columns between rows."


def _try_parse_json(text: str) -> list[dict] | None:
    text = text.strip()
    match = re.search(r'\[.*\]', text, re.DOTALL)
    if match:
        try:
            return json.loads(match.group(0))
        except json.JSONDecodeError:
            pass
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return None


def main():
    if len(sys.argv) < 3:
        print("Usage: python extract_page_direct.py <pdf_file> <page_number>")
        sys.exit(1)

    pdf_path = Path(sys.argv[1])
    if not pdf_path.exists():
        print(f"Error: PDF not found at {pdf_path}")
        sys.exit(1)

    try:
        page_num = int(sys.argv[2])
    except ValueError:
        print("Error: page_number must be an integer.")
        sys.exit(1)

    print(f"📄 Converting page {page_num} of {pdf_path.name} to image...")
    # Convert exactly one page (pdf2image uses 1-based indexing for first_page and last_page)
    images = convert_from_path(pdf_path, dpi=300, first_page=page_num, last_page=page_num)
    if not images:
        print("Error: Could not extract that page.")
        sys.exit(1)
    
    img = images[0]
    # Save a temporary image file needed for the SDK upload
    tmp_img_path = f"tmp_page_{page_num}.png"
    img.save(tmp_img_path, "PNG")

    print("🤖 Uploading to Gemini and extracting table... (this may take a minute)")
    
    # Upload to Gemini File API (recommended for multimodal)
    myfile = genai.upload_file(tmp_img_path)
    
    model = genai.GenerativeModel(
        model_name="gemini-2.5-pro",
        system_instruction=SYSTEM_PROMPT,
        generation_config={"temperature": 0.1}
    )

    response = model.generate_content([USER_PROMPT, myfile])
    content = response.text

    # Clean up the uploaded file and local temp file
    genai.delete_file(myfile.name)
    os.remove(tmp_img_path)

    rows = _try_parse_json(content)
    
    if not rows:
        print("\n❌ Failed to parse JSON from response. Raw response was:")
        print("---")
        print(content)
        print("---")
        sys.exit(1)

    print(f"\n✅ Successfully extracted {len(rows)} rows.")
    
    # Format and save to CSV
    # Ensure J columns are stored as text by padding them with Excel's literal string marker
    df = pd.DataFrame(rows)
    for j_col in ["lower_j", "upper_j"]:
        if j_col in df.columns:
            df[j_col] = df[j_col].apply(
                lambda v: f'="{v}"' if pd.notna(v) and str(v).strip() != "" and "/" in str(v) else v
            )
            
    out_file = f"page_{page_num}_raw_extraction.csv"
    df.to_csv(out_file, index=False, encoding="utf-8-sig")
    print(f"💾 Saved raw output to: {out_file}")


if __name__ == "__main__":
    main()
