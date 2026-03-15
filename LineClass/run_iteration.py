#!/usr/bin/env python3
"""
run_iteration.py
================
Orchestrates one iteration of the classify → weed → export pipeline.

Usage:
  python run_iteration.py --iteration 1
  python run_iteration.py --iteration 2 --levels-file /path/to/optimized_levels.xlsm
  python run_iteration.py --iteration 1 --skip-weeder     # classify only
  python run_iteration.py --iteration 2 --skip-classify    # weeder only (reuse previous classify output)

Each iteration creates an iterations/N/ directory containing:
  - metadata.json       — records which levels file and parameters were used
  - weeder_input.csv    — input for the LLM weeder
  - weeder_output.csv   — output from the LLM weeder (if weeder was run)

After the weeder output is reviewed, the user re-optimizes energy levels
and runs the next iteration with --levels-file pointing to the new levels.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ITERATIONS_DIR = os.path.join(SCRIPT_DIR, 'iterations')
DEFAULT_LEVELS_FILE = os.path.join(SCRIPT_DIR, '..', 'TableExtraction', 'Pr3_lev_Wyart_1999.xlsm')


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def run_step(cmd, description):
    """Run a subprocess and print its output in real time."""
    print(f"\n{'='*60}")
    print(f"  {description}")
    print(f"{'='*60}\n")
    result = subprocess.run(cmd, cwd=SCRIPT_DIR)
    if result.returncode != 0:
        print(f"\n❌ FAILED: {description} (exit code {result.returncode})")
        sys.exit(result.returncode)
    print(f"\n✅ {description} completed.")


def main():
    ap = argparse.ArgumentParser(
        description="Run one iteration of the classify → weed → export pipeline."
    )
    ap.add_argument("--iteration", "-n", type=int, required=True,
                    help="Iteration number (1, 2, 3, ...)")
    ap.add_argument("--levels-file", dest="levels_file", default=None,
                    help="Path to the energy levels file. Defaults to the original Wyart levels.")
    ap.add_argument("--skip-classify", dest="skip_classify", action="store_true",
                    help="Skip classify_lines.py (reuse existing output)")
    ap.add_argument("--skip-weeder", dest="skip_weeder", action="store_true",
                    help="Skip llm_weeder.py (classify only)")
    ap.add_argument("--weeder-args", dest="weeder_args", nargs=argparse.REMAINDER,
                    default=[], help="Additional arguments to pass to llm_weeder.py")
    args = ap.parse_args()

    iteration_n = args.iteration
    levels_file = os.path.abspath(args.levels_file) if args.levels_file else DEFAULT_LEVELS_FILE

    # Create iteration directory
    iter_dir = ensure_dir(os.path.join(ITERATIONS_DIR, str(iteration_n)))

    # --- Metadata ---
    metadata = {
        "iteration": iteration_n,
        "levels_file": levels_file,
        "timestamp": datetime.now().isoformat(),
        "skipped_classify": args.skip_classify,
        "skipped_weeder": args.skip_weeder,
    }

    meta_path = os.path.join(iter_dir, "metadata.json")
    with open(meta_path, "w") as f:
        json.dump(metadata, f, indent=2)
    print(f"Iteration {iteration_n} — metadata written to {meta_path}")
    print(f"  Levels file: {levels_file}")

    # --- Step 1: Classify ---
    if not args.skip_classify:
        # If a custom levels file is provided, we need to tell classify_lines.py
        # Currently classify_lines.py uses a hardcoded LEVELS_FILE path.
        # For iterations with custom levels, we temporarily symlink/copy.
        if args.levels_file:
            dest = os.path.join(SCRIPT_DIR, '..', 'TableExtraction', 'Pr3_lev_Wyart_1999.xlsm')
            backup = dest + '.bak'
            if os.path.exists(dest) and not os.path.exists(backup):
                shutil.copy2(dest, backup)
                print(f"  Backed up original levels to {backup}")
            shutil.copy2(levels_file, dest)
            print(f"  Installed levels file: {levels_file} -> {dest}")

        run_step([sys.executable, "classify_lines.py"], "classify_lines.py")
    else:
        print("\n⏭ Skipping classify_lines.py (--skip-classify)")

    # --- Copy classify output to iteration dir ---
    classify_csv = os.path.join(SCRIPT_DIR, "line_classifications.csv")
    weeder_input_src = os.path.join(SCRIPT_DIR, "weeder_input.csv")

    if os.path.exists(classify_csv):
        shutil.copy2(classify_csv, os.path.join(iter_dir, "line_classifications.csv"))
    if os.path.exists(weeder_input_src):
        iter_weeder_input = os.path.join(iter_dir, "weeder_input.csv")
        shutil.copy2(weeder_input_src, iter_weeder_input)
        print(f"  Copied weeder_input.csv -> {iter_weeder_input}")

    # --- Step 2: LLM Weeder ---
    if not args.skip_weeder:
        iter_weeder_input = os.path.join(iter_dir, "weeder_input.csv")
        iter_weeder_output = os.path.join(iter_dir, "weeder_output.csv")

        if not os.path.exists(iter_weeder_input):
            print(f"\n❌ Cannot run weeder: {iter_weeder_input} not found.")
            sys.exit(1)

        weeder_cmd = [
            sys.executable, "llm_weeder.py",
            "--in", iter_weeder_input,
            "--out", iter_weeder_output,
        ] + args.weeder_args

        run_step(weeder_cmd, "llm_weeder.py")
    else:
        print("\n⏭ Skipping llm_weeder.py (--skip-weeder)")

    # --- Summary ---
    print(f"\n{'='*60}")
    print(f"  Iteration {iteration_n} complete")
    print(f"{'='*60}")
    print(f"  Output directory: {iter_dir}")
    print(f"  Files:")
    for f in sorted(os.listdir(iter_dir)):
        fpath = os.path.join(iter_dir, f)
        size = os.path.getsize(fpath)
        print(f"    {f} ({size:,} bytes)")
    print()
    print("Next steps:")
    print(f"  1. Review {os.path.join(iter_dir, 'weeder_output.csv')} (if weeder was run)")
    print(f"  2. Re-optimize energy levels based on accepted assignments")
    print(f"  3. Run: python run_iteration.py --iteration {iteration_n + 1} --levels-file /path/to/new_levels.xlsm")


if __name__ == "__main__":
    main()
