#!/usr/bin/env python3
"""
Pr III Table Extraction Pipeline — CLI Entry Point.

Extracts tables from the Sugar Pr III PDF using Gemini multimodal vision,
validates against Edlen's formula and NIST ASD, and outputs CSV/Excel.

Usage:
    python extract_table.py                           # Process all pages
    python extract_table.py --pages 4 36              # Process pages 4–36
    python extract_table.py --single-page 5           # Process page 5 only
    python extract_table.py --skip-images             # Skip PDF→image conversion
    python extract_table.py --model gemini-2.5-pro    # Use specific model
"""

import argparse
import sys
import os
os.environ["PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION"] = "python"

from config import validate_api_key, get_config, GEMINI_MODEL, PDF_DPI
from graph import run_extraction_pipeline


def main():
    parser = argparse.ArgumentParser(
        description="Extract Pr III Tables from Sugar PDFs using Gemini + LangGraph",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python extract_table.py --pages 4 36              Process pages 4 through 36
  python extract_table.py --single-page 5           Process only page 5
  python extract_table.py --skip-images             Reuse existing page images
  python extract_table.py --model gemini-2.5-pro    Use a different Gemini model
  python extract_table.py --dpi 200                 Lower DPI for faster processing
  python extract_table.py --dry-run                 Show config without processing
        """,
    )

    parser.add_argument(
        "--pages", type=int, nargs=2, metavar=("START", "END"),
        help="Page range to process (1-indexed, inclusive). Default: all pages.",
    )
    parser.add_argument(
        "--single-page", type=int, metavar="N",
        help="Process a single page only (1-indexed). Overrides --pages.",
    )
    parser.add_argument(
        "--model", type=str, default=GEMINI_MODEL,
        help=f"Gemini model to use (default: {GEMINI_MODEL}). Options: gemini-2.5-flash, gemini-2.5-pro, gemini-2.0-flash",
    )
    parser.add_argument(
        "--dpi", type=int, default=PDF_DPI,
        help=f"Image DPI for PDF conversion (default: {PDF_DPI}). Higher = better quality but more tokens.",
    )
    parser.add_argument(
        "--pdf", type=str, default=None,
        help="Path to the PDF file (default: from config based on year).",
    )
    parser.add_argument(
        "--skip-images", action="store_true",
        help="Skip PDF-to-image conversion; use existing images in page_images/.",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Show configuration and exit without processing.",
    )

    args = parser.parse_args()

    cfg = get_config()

    # Resolve page range
    first_page = None
    last_page = None

    if args.single_page:
        first_page = args.single_page
        last_page = args.single_page
    elif args.pages:
        first_page = args.pages[0]
        last_page = args.pages[1]

    actual_pdf = args.pdf if args.pdf else str(cfg["pdf_path"])

    # Dry run (no API key needed)
    if args.dry_run:
        print(f"🔍 Dry Run — Configuration:")
        print(f"   PDF:       {actual_pdf}")
        print(f"   Model:     {args.model}")
        print(f"   Pages:     {first_page or 'first'} – {last_page or 'last'}")
        print(f"   DPI:       {args.dpi}")
        print(f"   Skip imgs: {args.skip_images}")
        print(f"   CSV Out:   {cfg['csv_name']}")
        print(f"   Excel Out: {cfg['excel_name']}")
        print(f"   Wavelength:{cfg['wavelength_type']}")
        print(f"   Pub Year:  {cfg['publication_year']}")
        print(f"   Air/Vac:   {cfg['air_vac_conversion_method']}")
        sys.exit(0)

    # Validate API key (needed for actual processing)
    validate_api_key()

    # Run the extraction pipeline
    result = run_extraction_pipeline(
        pdf_path=actual_pdf,
        first_page=first_page,
        last_page=last_page,
        dpi=args.dpi,
        model_name=args.model,
        skip_images=args.skip_images,
    )

    # Exit code based on errors
    if result.get("errors"):
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
