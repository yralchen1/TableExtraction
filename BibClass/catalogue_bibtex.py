#!/usr/bin/env python3
"""
Script to extract bibtex entries from bibtex_Sr1_1976-2025.txt and place them 
in their corresponding paper folders as bibtex_catalogued.txt
"""

import os
import re

# Paths
BASE_DIR = "/Users/themanaspandey/Documents/Alexander Marianna Yuri Automations/Benchmarking PDFs"
BIBTEX_FILE = os.path.join(BASE_DIR, "bibtex_Sr1_1976-2025.txt")

def parse_bibtex_entries(filepath):
    """Parse bibtex file and return a dictionary of {id: full_entry}"""
    entries = {}
    
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Split by @ to get individual entries
    # Pattern: @article{ID, ... } or @incollection{ID, ...} or any other entry type
    pattern = r'(@\w+\{(\d+)EL,.*?(?=\n@|\Z))'
    matches = re.findall(pattern, content, re.DOTALL)
    
    for full_entry, entry_id in matches:
        entries[entry_id] = full_entry.strip()
    
    return entries

def get_folder_ids(base_dir):
    """Get all paper folders and extract their IDs"""
    folders = {}
    
    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        if os.path.isdir(item_path):
            # Extract ID from folder name pattern: YEAR_Author_el_ID
            match = re.search(r'_el_(\d+)$', item)
            if match:
                folder_id = match.group(1)
                folders[folder_id] = item_path
    
    return folders

def main():
    print("Parsing bibtex entries...")
    entries = parse_bibtex_entries(BIBTEX_FILE)
    print(f"Found {len(entries)} bibtex entries")
    
    print("\nScanning folders...")
    folders = get_folder_ids(BASE_DIR)
    print(f"Found {len(folders)} paper folders")
    
    # Match and write
    matched = 0
    unmatched_folders = []
    
    for folder_id, folder_path in folders.items():
        if folder_id in entries:
            output_file = os.path.join(folder_path, "bibtex_catalogued.txt")
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(entries[folder_id])
            print(f"✓ Wrote bibtex for ID {folder_id} to {os.path.basename(folder_path)}")
            matched += 1
        else:
            unmatched_folders.append((folder_id, os.path.basename(folder_path)))
    
    print(f"\n{'='*50}")
    print(f"Summary:")
    print(f"  Total bibtex entries: {len(entries)}")
    print(f"  Total paper folders: {len(folders)}")
    print(f"  Successfully matched: {matched}")
    
    if unmatched_folders:
        print(f"\n  Unmatched folders ({len(unmatched_folders)}):")
        for fid, fname in unmatched_folders:
            print(f"    - {fname} (ID: {fid})")

if __name__ == "__main__":
    main()
