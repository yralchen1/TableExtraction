#!/usr/bin/env python3
"""
Pr III Table Validation Pipeline — CLI Entry Point.

Reads an intermediate raw OCR Excel file (.xlsx) containing `all_rows`,
runs the validation pipeline (NIST checks, selection rules, wavelength/wavenumber
conversions, page misalignment fixes), and outputs the final CSV and structured Excel files.

Usage:
    python validate_extracted_table.py --input output/raw_extracted_table_1969.xlsx
"""

import argparse
import sys
import os
from pathlib import Path

os.environ["PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION"] = "python"

from config import get_config
from graph import run_validation_pipeline


def main():
    parser = argparse.ArgumentParser(
        description="Validate raw extracted Pr III tables (from Excel) using LangGraph",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python validate_extracted_table.py --input output/raw_extracted_table_1969.xlsx
  python validate_extracted_table.py --input output/raw_extracted_table_1969.xlsx --dry-run
        """,
    )

    parser.add_argument(
        "--input", "-i", type=str, required=True,
        help="Path to the raw extracted Excel file (.xlsx) output by extract_table.py.",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Show configuration and file paths without processing.",
    )

    args = parser.parse_args()

    cfg = get_config()
    input_path = Path(args.input)

    if not input_path.exists():
        print(f"❌ Error: Input file not found: {input_path}")
        sys.exit(1)

    # Dry run
    if args.dry_run:
        print(f"🔍 Dry Run — Validation Configuration:")
        print(f"   Input Raw Excel: {input_path}")
        print(f"   CSV Out:         {cfg['csv_name']}")
        print(f"   Excel Out:       {cfg['excel_name']}")
        print(f"   Wavelength:      {cfg['wavelength_type']}")
        print(f"   Pub Year:        {cfg['publication_year']}")
        print(f"   Air/Vac Method:  {cfg['air_vac_conversion_method']}")
        print(f"   Tolerance (WL):  {cfg['wavenumber_tolerance_cm1']} cm⁻¹")
        print(f"   Tolerance (Ritz):{cfg['ref_ritz_tolerance_cm1']} cm⁻¹")
        sys.exit(0)

    # Run the validation pipeline
    result = run_validation_pipeline(raw_xlsx_path=str(input_path))

    # Exit code based on errors
    if result.get("errors"):
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
