#!/usr/bin/env python3
"""
llm_weeder.py
=================
LLM-based weeding of improbable assignments (per observed line).

Input: CSV where each row is a candidate assignment for an observed line.
Output: same rows + columns: accepted (0/1), reason_code, reasoning

Requires:
  pip install -U google-genai pandas

Usage:
  python llm_weeder.py --in candidates.csv --out candidates_scored.csv --model models/gemini-3-flash-preview --api-key YOUR_KEY

Notes:
- The script batches by observed-line key and caps total candidates per request.
- It validates JSON strictly and retries once on transient failures.
"""

import argparse
import json
import math
import os
import time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import pandas as pd
from google import genai
from google.genai import errors as genai_errors
from google.genai import types as genai_types

from dotenv import load_dotenv
# from pathlib import Path

# Load environment variables from .env
load_dotenv()

# --- API Key ---
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY", "")

# ---------------------------
# CONFIG DEFAULTS
# ---------------------------

GEMINI_MODEL = os.getenv("GEMINI_MODEL", "gemini-3-flash-preview")
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Input/output files
INPUT_PATH = str(SCRIPT_DIR + "/weeder_input.csv")
OUTPUT_PATH = str(SCRIPT_DIR + "/weeder_output.csv")


# Batch controls
DEFAULT_MAX_LINES_PER_BATCH = 500         # unique observed lines
DEFAULT_MAX_CANDIDATES_PER_BATCH = 700    # total candidate rows

# Throttling (avoid token/min bursts; adjust as needed)
DEFAULT_MIN_SECONDS_BETWEEN_CALLS = 60.0

# Output control
ACCEPT_COL = "accepted"
REASON_CODE_COL = "reason_code"
REASON_COL = "reasoning"


CANONICAL_REASON_CODES = [
    "STRONG_Z",
    "MODERATE_Z",
    "INT_MATCH",
    "Z_MISMATCH",
    "ABS_MISM",
    "I_UNCERT",
    "I_MISMATCH",
    "ORIG_MATCH",
    "MULTI_SUM",
    "UNCERTAIN",
    "MISSING_DATA",
]

# --- Structured output schema (array of objects) for Gemini ---
WEEDER_OUTPUT_SCHEMA = {
    "type": "array",
    "items": {
        "type": "object",
        "properties": {
            "row_index": {"type": "integer"},
            "accepted": {"type": "integer", "enum": [0, 1]},
            "reason_code": {"type": "string", "enum": CANONICAL_REASON_CODES},
            "reasoning": {"type": "string"}
        },
        "required": ["row_index", "accepted", "reason_code", "reasoning"]
    }
}

# --- Local pre-filter thresholds (adjust as needed) ---
PREFILTER_Z_REJECT = 8.0       # Auto-reject candidates with z > this
PREFILTER_ACCEPT_GRADES = {'2A', '2B'}  # Auto-accept sole originals with these grades

# Keep your long “physics caveats” here (once per request).
# You can edit this text freely.
PROMPT_PREAMBLE = """Context (read once, apply to all rows in this batch):

You are assisting with classification of atomic spectral lines.

Your task is to evaluate candidate transition assignments for observed spectral lines and decide whether each candidate should be accepted for further analysis.

You must follow the rules below strictly. Do not invent new rules or external knowledge.

Definitions:
- obs_wn: observed wavenumber of the line
- rwn: Ritz wavenumber predicted from energy levels (E_upper - E_lower).
- Δwn = obs_wn - rwn; wn_diff_abs = |Δwn|.
- sigma: σ = uncertainty of obs_wn.
- z = wn_diff_abs / σ.

Intensity information:
- obs_I: observed intensity (approximate, derived from photographic plate blackening)
- calc_I: calculated transition intensity (approximate, may be missing)
- u_calc: uncertainty estimate for calc_I when available; it is given on a natural logarithmic scale, converted from percentage uncertainty of calc_I, calc_I_unc_percent: u_calc = ln(calc_I_unc_percent / 100 + 1). 

Intensity guidance:
- intensity agreement is weak evidence
- |ln(calc_I/obs_I)| ≤ ln(2) → good
- |ln(calc_I/obs_I)| ≥ ln(5) → strong mismatch
- missing calc_I means intensity cannot be used → rely on wavenumber consistency, flags, line character, and whether the classification is new or not.

Line character (line_char) guidance:
- Codes 'cl', 'bl', 'h', 'w', 'ch', 'd', '*r', or '*v' denote lines that have complex profile and are likely significantly broader than those without such characters. These lines are more likely to be blends of several transitions.
- Characters '**' denote lines that had multiple original classifications.     

Grades (heuristic summary from earlier code):
- first symbol encodes z:
  2: z ≤ 2
  3: z ≤ 3
  4: z ≤ 4
  5: z > 4
- second symbol summarizes intensity agreement and reliability:
  A = very good match with dependable calc_I
  B = satisfactory match with fairly dependable calc_I
  C = moderately poor match with moderately uncertain calc_I
  D = poor match with poorly known calc_I
  E = very poor match with poorly known calc_I
  G = extremely poor match or calc_I is unknown

`new` flag:
- 1 means this is a new classification and should be treated with extra caution

Important:
- Prefer previously established classifications unless evidence strongly contradicts them.
- New classifications require stronger evidence.
- If uncertain between two candidates, prefer rejecting the weaker one rather than accepting both.
- Decisions must rely only on the provided data.
- Use grade only as a secondary heuristic. Prefer the raw fields (obs_wn, sigma, obs_I, calc_I, u_calc, rwn, `new`) if they suggest a different conclusion.

Data caveats:
- Observed intensities are approximated from photographic plate blackening, scaled roughly to a linear scale; they can have systematic errors (often within a factor ~2–3).
- The scale of observed intensities (obs_I) is approximately the same as the scale of calculated intensities (calc_I). 
- Calculated intensities are incomplete and based on a crude Boltzmann population model; they may be missing and have uncertainties, which are themselves unreliable, more so for greater uncertainties.
- Observed wavenumber uncertainties are approximate; energy levels may be refined later, which can reduce uncertainty estimates for some true assignments.
- Use intensity as weak evidence unless mismatch is large.

Decision procedure (apply in this order):

1. Wavenumber test
   - If z > 8 → reject (Z_MISMATCH)
   - If 5 < z ≤ 8 → usually reject unless other evidence strongly supports it

2. Absolute residual comparison
   - If several candidates exist for one line, prefer smaller wn_diff_abs

3. Intensity compatibility
   - If calc_I and obs_I exist and mismatch ≥ factor 5 → reject unless wavenumber agreement is excellent
   - If mismatch ≤ factor 2 → weak support

4. Multiple-candidate pruning
   - Sort candidates by calc_I descending
   - Maintain running sum S of accepted calc_I
   - If S ≥ 1.5 × obs_I:
       reject candidates with calc_I < 0.1 × S
       unless their z is significantly smaller

5. Sign consistency test
   - If multiple accepted candidates have Δwn values with opposite signs, this supports a blended line explanation
   - If all Δwn have the same sign and offsets are large, treat with suspicion

6. Conservative rule
   - If `new` is 1 (new classification), require stronger evidence than for existing ones

7. If evidence is inconclusive → mark UNCERTAIN

Input rows format (CSV example):
id, obs_wn, sigma, I_obs, line_char, rwn, I_calc, u_calc, grade, new
3,121547.2990863,0.5909498,1.0264487,,1,102644.18,,,2G,1
5,121279.4157051,0.5883480,0.5195639,,2,121280.90,,,3G,0
...

Meaning:
id     = sequential row identifier
obs_wn = observed wavenumber (unique for an observed line)
sigma  = uncertainty of obs_wn
I_obs  = observed intensity
char   = line character string
rwn    = lower level energy
I_calc = calculated intensity (may be null)
u_calc = uncertainty of I_calc (may be null)
grade  = grading code
new    = 1 if new classification else 0

"""

# ---------------------------
# Helpers
# ---------------------------

@dataclass
class BatchItem:
    row_index: int
    row: Dict[str, Any]


def pick_obs_key(df: pd.DataFrame, user_key: Optional[str]) -> str:
    if user_key:
        if user_key not in df.columns:
            raise ValueError(f"--obs-key '{user_key}' not found in CSV columns: {list(df.columns)}")
        return user_key
    # Sensible defaults:
    for k in ["obs_id", "obs_wn", "observed_wavenumber", "wn_obs"]:
        if k in df.columns:
            return k
    raise ValueError("Could not infer observed-line key. Provide --obs-key (e.g., obs_id or obs_wn).")


def safe_float(x: Any) -> Optional[float]:
    try:
        if x is None or (isinstance(x, float) and math.isnan(x)):
            return None
        return float(x)
    except (ValueError, TypeError, OverflowError):
        return None


def build_rows_payload(df_batch: pd.DataFrame, required_cols: List[str]) -> List[Dict[str, Any]]:
    payload = []
    for _, r in df_batch.iterrows():
        obj = {c: (None if (pd.isna(r[c]) if c in df_batch.columns else True) else r[c]) for c in required_cols}
        payload.append(obj)
    return payload


def pre_filter(df: pd.DataFrame, obs_key: str) -> tuple:
    """Auto-decide obvious accept/reject cases without calling the LLM.

    Returns:
        decided_indices: dict mapping row index -> (accepted, reason_code, reasoning)
        undecided_df: DataFrame of rows that still need LLM evaluation
    """
    decided = {}

    for obs_wn, group in df.groupby(obs_key):
        for idx, row in group.iterrows():
            sigma = safe_float(row.get('sigma'))
            rwn = safe_float(row.get('rwn'))
            obs_wn_val = safe_float(row.get('obs_wn'))
            grade = str(row.get('grade', '')) if pd.notna(row.get('grade')) else ''
            new_flag = int(row.get('new', 1)) if pd.notna(row.get('new')) else 1

            if sigma and sigma > 0 and obs_wn_val is not None and rwn is not None:
                z = abs(obs_wn_val - rwn) / sigma
            else:
                z = None

            # Auto-reject: z > threshold
            if z is not None and z > PREFILTER_Z_REJECT:
                decided[idx] = (0, 'Z_MISMATCH', f'z={z:.1f}; auto-reject (pre-filter)')
                continue

            # Auto-accept: sole original candidate with excellent grade
            if (len(group) == 1
                    and grade in PREFILTER_ACCEPT_GRADES
                    and new_flag == 0
                    and z is not None and z <= 2.0):
                decided[idx] = (1, 'STRONG_Z', f'z={z:.1f}; sole orig {grade}; auto-accept (pre-filter)')
                continue

    undecided_df = df.drop(index=list(decided.keys()))
    return decided, undecided_df


def build_prompt(rows_payload: List[Dict[str, Any]]) -> str:
    # Output format: array of objects (enforced by structured output schema)

    schema_text = f"""
Output format:
Return exactly one JSON array of objects.
Each object must have these fields:
- row_index: integer, 0-based index of the input row in this batch
- accepted: integer, either 0 or 1
- reason_code: one of {CANONICAL_REASON_CODES}
- reasoning: very short string, at most 80 characters
  - prefer semicolon-separated fragments instead of full sentences
  - do not repeat field names unnecessarily

Rules:
- Return exactly {len(rows_payload)} output objects.
- Preserve input row order.
- Do not omit any row.
- Do not add any extra rows.
- If critical data are missing, use MISSING_DATA, I_UNCERT, or UNCERTAIN as appropriate.

Example output:
[
  {{"row_index": 0, "accepted": 1, "reason_code": "STRONG_Z", "reasoning": "z=1.2; old cls; no conflict"}},
  {{"row_index": 1, "accepted": 0, "reason_code": "Z_MISMATCH", "reasoning": "z=7.8; reject"}}
]
""".strip()

    instructions = f"""
Input rows are given below as one JSON array.
Each input row already contains its own row_index field.

Use only the provided data.
Return only the JSON output array in the required format.

Input rows:
""".strip()

    prompt = (
        PROMPT_PREAMBLE.strip()
        + "\n\n"
        + schema_text
        + "\n\n"
        + instructions
        + "\n"
        + json.dumps(rows_payload, ensure_ascii=False, separators=(",", ":"))
    )
    return prompt


def call_model(client: genai.Client, model: str, prompt: str,
               use_structured_output: bool = True, max_retries: int = 1) -> str:
    """Call Gemini with optional structured JSON output enforcement."""
    last_err = None
    for attempt in range(max_retries + 1):
        try:
            config = None
            if use_structured_output:
                config = genai_types.GenerateContentConfig(
                    response_mime_type='application/json',
                    response_json_schema=WEEDER_OUTPUT_SCHEMA,
                )
            resp = client.models.generate_content(
                model=model,
                contents=prompt,
                config=config,
            )
            # google-genai returns response.text
            return resp.text or ""
        except (genai_errors.ClientError, genai_errors.ServerError) as e:
            last_err = e
            # Retry only on typical transient conditions (429/503)
            status = getattr(e, "status_code", None)
            msg = getattr(e, "message", str(e))
            if attempt < max_retries and (status in (429, 503) or "UNAVAILABLE" in msg or "RESOURCE_EXHAUSTED" in msg):
                # basic backoff
                time.sleep(5 + 2 * attempt)
                continue
            raise
    raise last_err  # pragma: no cover


def parse_model_json(text: str, expected_n: int) -> List[Dict[str, Any]]:
    """Parse and validate the model's JSON output (array of objects)."""
    text = text.strip()

    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        # Try to recover if the model wrapped the JSON with accidental extra text
        start = text.find("[")
        end = text.rfind("]")
        if start != -1 and end != -1 and end > start:
            data = json.loads(text[start:end + 1])
        else:
            raise ValueError("Model did not return a valid JSON array.")

    if not isinstance(data, list):
        raise ValueError(f"Expected top-level JSON array, got {type(data)}")

    if len(data) != expected_n:
        raise ValueError(
            f"Expected JSON array of length {expected_n}, got {len(data)}"
        )

    seen_row_indices = set()

    for i, item in enumerate(data):
        if not isinstance(item, dict):
            raise ValueError(f"Output element {i} is not an object (got {type(item).__name__}).")

        for key in ('row_index', 'accepted', 'reason_code', 'reasoning'):
            if key not in item:
                raise ValueError(f"Output element {i} missing required field '{key}'.")

        row_index = item['row_index']
        accepted = item['accepted']
        reason_code = item['reason_code']
        reasoning = item['reasoning']

        if not isinstance(row_index, int):
            raise ValueError(f"Output element {i} has non-integer row_index: {row_index!r}")

        if row_index in seen_row_indices:
            raise ValueError(f"Duplicate row_index in output: {row_index}")
        seen_row_indices.add(row_index)

        if not isinstance(accepted, int) or accepted not in (0, 1):
            raise ValueError(f"Output element {i} has invalid accepted value: {accepted!r}")

        if not isinstance(reason_code, str) or reason_code not in CANONICAL_REASON_CODES:
            raise ValueError(
                f"Output element {i} has invalid reason_code: {reason_code!r}"
            )

        if not isinstance(reasoning, str):
            raise ValueError(f"Output element {i} has non-string reasoning: {reasoning!r}")

        if len(reasoning) > 120:
            # Warn but don't fail — structured output may produce slightly longer text
            print(f"  Warning: element {i} reasoning is {len(reasoning)} chars (truncating to 120)")
            item['reasoning'] = reasoning[:120]

    missing = set(range(expected_n)) - seen_row_indices
    extra = seen_row_indices - set(range(expected_n))
    if missing:
        raise ValueError(f"Missing row_index values in output: {sorted(missing)[:10]}")
    if extra:
        raise ValueError(f"Out-of-range row_index values in output: {sorted(extra)[:10]}")

    return data


def make_batches(df: pd.DataFrame, obs_key: str, max_lines: int, max_candidates: int) -> List[pd.DataFrame]:
    batches: List[pd.DataFrame] = []
    # Group by observed line, preserve order by first appearance
    # (stable ordering helps reproducibility)
    seen = []
    seen_set = set()
    for v in df[obs_key].tolist():
        if v not in seen_set:
            seen_set.add(v)
            seen.append(v)

    i = 0
    while i < len(seen):
        lines = []
        cand_count = 0
        while i < len(seen) and len(lines) < max_lines:
            key_val = seen[i]
            group = df[df[obs_key] == key_val]
            gsz = len(group)
            # If adding this group would exceed candidate cap, stop batch (unless empty batch)
            if lines and (cand_count + gsz) > max_candidates:
                break
            lines.append(key_val)
            cand_count += gsz
            i += 1

        batch_df = df[df[obs_key].isin(lines)].copy()
        # Keep original row order
        batch_df = batch_df.sort_index()
        batches.append(batch_df)

    return batches


# noinspection PyArgumentList
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_csv", default=INPUT_PATH, help="Input candidates CSV")
    ap.add_argument("--out", dest="out_csv", default=OUTPUT_PATH, help="Output CSV path")
    ap.add_argument("--api-key", dest="api_key", default=GOOGLE_API_KEY, help="Gemini API key")
    ap.add_argument("--model", dest="model", default=GEMINI_MODEL, help=f"Model name (default: {GEMINI_MODEL})")
    ap.add_argument("--obs-key", dest="obs_key", default='obs_wn', help="Column used to group candidates per observed line (e.g., obs_id or obs_wn)")
    ap.add_argument("--max-lines", dest="max_lines", type=int, default=DEFAULT_MAX_LINES_PER_BATCH)
    ap.add_argument("--max-candidates", dest="max_candidates", type=int, default=DEFAULT_MAX_CANDIDATES_PER_BATCH)
    ap.add_argument("--min-delay", dest="min_delay", type=float, default=DEFAULT_MIN_SECONDS_BETWEEN_CALLS)
    ap.add_argument("--dry-run", dest="dry_run", action="store_true", help="Build batches and prompts but do not call model")
    ap.add_argument("--no-prefilter", dest="no_prefilter", action="store_true", help="Disable local pre-filter (send all rows to LLM)")
    ap.add_argument("--no-structured-output", dest="no_structured", action="store_true", help="Disable structured JSON output mode")
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv)
    obs_key = pick_obs_key(df, args.obs_key)

    # Initialize output cols
    for c in [ACCEPT_COL, REASON_CODE_COL, REASON_COL]:
        if c not in df.columns:
            df[c] = None

    # --- Local pre-filter ---
    n_prefiltered = 0
    if not args.no_prefilter:
        decided, undecided_df = pre_filter(df, obs_key)
        n_prefiltered = len(decided)
        for idx, (accepted, reason_code, reasoning) in decided.items():
            df.at[idx, ACCEPT_COL] = accepted
            df.at[idx, REASON_CODE_COL] = reason_code
            df.at[idx, REASON_COL] = reasoning
        print(f"Pre-filter: {n_prefiltered} rows auto-decided, {len(undecided_df)} rows remain for LLM.")
    else:
        undecided_df = df
        print("Pre-filter disabled.")

    # Choose which columns to send to the model:
    # Adjust this list to match your CSV headers.
    # Keep it minimal to save tokens.
    required_cols = [
        obs_key, "id", "obs_wn", "sigma", "obs_I", "line_char",
        "rwn", "calc_I", "u_calc", "grade", "new"
    ]
    # Keep only columns that actually exist; the model will see None for missing fields
    required_cols = [c for c in required_cols if c in df.columns]

    batches = make_batches(undecided_df, obs_key, args.max_lines, args.max_candidates)
    print(f"Using obs_key='{obs_key}'. Undecided rows={len(undecided_df)}. Batches={len(batches)}.")
    print(f"Batch caps: max_lines={args.max_lines}, max_candidates={args.max_candidates}")

    use_structured = not args.no_structured
    if use_structured:
        print("Structured JSON output mode: ENABLED")
    else:
        print("Structured JSON output mode: DISABLED")

    client = None
    if not args.dry_run:
        client = genai.Client(api_key=args.api_key)

    last_call_t = 0.0
    for bi, df_batch in enumerate(batches, start=1):
        idxs = df_batch.index.tolist()

        # --- Checkpoint/resume: skip batches already completed ---
        if not pd.isna(df.at[idxs[0], ACCEPT_COL]):
            print(f"Batch {bi}/{len(batches)} already completed, skipping.")
            continue

        rows_payload = []
        for local_i, (row_idx, row) in enumerate(df_batch.iterrows()):
            obj = {c: (None if (pd.isna(row[c]) if c in df_batch.columns else True) else row[c]) for c in required_cols}
            obj["row_index"] = local_i  # batch-local index
            rows_payload.append(obj)

        prompt = build_prompt(rows_payload)

        if args.dry_run:
            print(f"\n--- DRY RUN batch {bi}/{len(batches)} rows={len(df_batch)} ---")
            print(prompt[:1500] + ("\n...[truncated]..." if len(prompt) > 1500 else ""))
            continue

        # Throttle
        now = time.time()
        wait = (last_call_t + args.min_delay) - now
        if wait > 0:
            time.sleep(wait)

        print(f"\nCalling model for batch {bi}/{len(batches)}: rows={len(df_batch)} unique_lines={df_batch[obs_key].nunique()}")
        text = call_model(client, args.model, prompt,
                          use_structured_output=use_structured, max_retries=1)
        last_call_t = time.time()

        out = parse_model_json(text, expected_n=len(df_batch))

        # Apply outputs back to df by original row order
        # Each output item is: {"row_index": ..., "accepted": ..., "reason_code": ..., "reasoning": ...}
        seen_local_indices = set()

        for item in out:
            local_i = int(item['row_index'])
            accepted = item['accepted']
            reason_code = item['reason_code']
            reasoning = item['reasoning']

            if local_i in seen_local_indices:
                raise ValueError(f"Duplicate row_index in model output: {local_i}")
            seen_local_indices.add(local_i)

            if not (0 <= local_i < len(idxs)):
                raise ValueError(f"row_index out of range in model output: {local_i}")

            row_idx = idxs[local_i]
            df.at[row_idx, ACCEPT_COL] = int(accepted)
            df.at[row_idx, REASON_CODE_COL] = reason_code
            df.at[row_idx, REASON_COL] = reasoning

        if len(seen_local_indices) != len(df_batch):
            missing = sorted(set(range(len(df_batch))) - seen_local_indices)
            raise ValueError(f"Model output is missing row_index values: {missing[:10]}")

        # Save incremental checkpoint
        df.to_csv(args.out_csv, index=False)
        print(f"Saved checkpoint: {args.out_csv}")

    # Final save (includes pre-filtered + LLM-decided rows)
    df.to_csv(args.out_csv, index=False)
    print(f"\nDONE. Final output: {args.out_csv}")
    print(f"  Pre-filtered: {n_prefiltered} rows")
    print(f"  LLM-decided:  {len(df) - n_prefiltered} rows")
    accepted_count = (df[ACCEPT_COL] == 1).sum() if ACCEPT_COL in df.columns else 0
    rejected_count = (df[ACCEPT_COL] == 0).sum() if ACCEPT_COL in df.columns else 0
    print(f"  Accepted: {accepted_count}, Rejected: {rejected_count}")


if __name__ == "__main__":
    main()