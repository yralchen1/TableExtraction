#!/usr/bin/env python3
"""
Script to organize PDFs into a structured folder format:

root_directory/
│
├── [YEAR_Author_Keyword]/      # Folder dedicated to each article
│   ├── main_article.pdf        # The primary PDF
│   └── suppl/                  # Sub-directory for supplements
│       ├── supplementary_note_1.pdf
│       ├── raw_data_table.csv
│       └── image_S1.png

Expected input filename format: Author_el_XXXXX_YEAR.pdf
Output folder format: YEAR_Author_el_XXXXX
"""

import os
import re
import shutil
import argparse


def parse_filename(filename):
    """
    Parse PDF filename to extract author, identifier, and year.
    Expected format: Author_el_XXXXX_YEAR.pdf
    
    Returns:
        tuple: (author, identifier, year) or None if pattern doesn't match
    """
    pattern = r'^(.+?)_el_(\d+)_(\d{4})\.pdf$'
    match = re.match(pattern, filename)
    
    if match:
        author = match.group(1)
        identifier = match.group(2)
        year = match.group(3)
        return author, identifier, year
    return None


def organize_pdfs(source_dir, dry_run=False, verbose=True):
    """
    Organize PDFs in the source directory into structured folders.
    
    Args:
        source_dir: Directory containing PDF files to organize
        dry_run: If True, only print what would be done without making changes
        verbose: If True, print detailed progress
    
    Returns:
        dict: Statistics about the organization process
    """
    stats = {
        'organized': 0,
        'skipped': [],
        'already_organized': 0,
        'errors': []
    }
    
    # Get all PDF files in the source directory (not in subdirectories)
    pdf_files = [f for f in os.listdir(source_dir) 
                 if f.endswith('.pdf') and os.path.isfile(os.path.join(source_dir, f))]
    
    if not pdf_files:
        print("No PDF files found in the source directory.")
        return stats
    
    print(f"Found {len(pdf_files)} PDF files to organize.\n")
    
    for pdf in pdf_files:
        parsed = parse_filename(pdf)
        
        if parsed:
            author, identifier, year = parsed
            
            # Create folder name: YEAR_Author_el_ID
            folder_name = f"{year}_{author}_el_{identifier}"
            folder_path = os.path.join(source_dir, folder_name)
            suppl_path = os.path.join(folder_path, "suppl")
            
            src = os.path.join(source_dir, pdf)
            dst = os.path.join(folder_path, "main_article.pdf")
            
            if dry_run:
                print(f"[DRY RUN] Would create: {folder_name}/")
                print(f"[DRY RUN] Would create: {folder_name}/suppl/")
                print(f"[DRY RUN] Would move: {pdf} -> {folder_name}/main_article.pdf")
                stats['organized'] += 1
            else:
                try:
                    # Create directory if it doesn't exist
                    os.makedirs(folder_path, exist_ok=True)
                    
                    # Create suppl subdirectory
                    os.makedirs(suppl_path, exist_ok=True)
                    
                    # Move and rename PDF to main_article.pdf
                    shutil.move(src, dst)
                    
                    if verbose:
                        print(f"✓ {pdf} -> {folder_name}/main_article.pdf")
                    
                    stats['organized'] += 1
                    
                except Exception as e:
                    error_msg = f"Error processing {pdf}: {str(e)}"
                    stats['errors'].append(error_msg)
                    print(f"✗ {error_msg}")
        else:
            stats['skipped'].append(pdf)
            if verbose:
                print(f"⊘ Skipped (pattern mismatch): {pdf}")
    
    return stats


def print_summary(stats):
    """Print a summary of the organization process."""
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"  Successfully organized: {stats['organized']}")
    print(f"  Already organized:      {stats['already_organized']}")
    print(f"  Skipped (no match):     {len(stats['skipped'])}")
    print(f"  Errors:                 {len(stats['errors'])}")
    
    if stats['skipped']:
        print(f"\nSkipped files (didn't match expected pattern):")
        for f in stats['skipped']:
            print(f"  - {f}")
    
    if stats['errors']:
        print(f"\nErrors encountered:")
        for e in stats['errors']:
            print(f"  - {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Organize PDFs into structured folders",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Expected input filename format: Author_el_XXXXX_YEAR.pdf
Output folder structure:
  YEAR_Author_el_XXXXX/
    ├── main_article.pdf
    └── suppl/
        """)
    
    parser.add_argument(
        'source_dir',
        nargs='?',
        default=os.path.dirname(os.path.abspath(__file__)),
        help='Directory containing PDFs to organize (default: script directory)'
    )
    
    parser.add_argument(
        '--dry-run', '-n',
        action='store_true',
        help='Show what would be done without making changes'
    )
    
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Only show summary, not individual file operations'
    )
    
    args = parser.parse_args()
    
    source_dir = os.path.abspath(args.source_dir)
    
    if not os.path.isdir(source_dir):
        print(f"Error: {source_dir} is not a valid directory")
        return 1
    
    print(f"Organizing PDFs in: {source_dir}")
    if args.dry_run:
        print("(DRY RUN - no changes will be made)\n")
    else:
        print()
    
    stats = organize_pdfs(source_dir, dry_run=args.dry_run, verbose=not args.quiet)
    print_summary(stats)
    
    return 0


if __name__ == "__main__":
    exit(main())
