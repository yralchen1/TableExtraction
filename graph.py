"""
LangGraph Pipeline for Pr III Table Extraction.

v4: Row structure fixes (doubly-classified lines, split-row merging),
    enhanced misalignment correction (all-row Ritz check).

Pipeline:
  convert_pdf → extract_tables → fix_row_structure → validate_rows
  → fix_misalignment → compile_output
"""

import json
import re
import base64
import traceback
from pathlib import Path
from typing import TypedDict
from collections import defaultdict

import pandas as pd
from langgraph.graph import StateGraph, START, END

from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_core.messages import HumanMessage, SystemMessage

from config import (
    get_config, GOOGLE_API_KEY, GEMINI_MODEL, PDF_DPI,
    PAGE_IMAGES_DIR, OUTPUT_DIR,
)
from pdf_to_images import convert_pdf_to_images, get_existing_images
from extraction_prompt import get_system_prompt, USER_PROMPT
from tools import (
    validate_full_row,
    wavenumber_from_vac_wavelength,
    wavelength_vac_from_wavenumber,
    wavenumber_from_air_wavelength,
    wavelength_air_from_wavenumber,
)


# ──────────────────────────────────────────────
# State definition
# ──────────────────────────────────────────────

class PipelineState(TypedDict):
    """State passed between graph nodes."""
    pdf_path: str
    first_page: int | None
    last_page: int | None
    dpi: int
    model_name: str
    skip_images: bool

    page_images: list[str]
    raw_extractions: list[dict]
    all_rows: list[dict]
    final_rows: list[dict]
    errors: list[str]
    output_csv: str
    output_excel: str


# ──────────────────────────────────────────────
# Robust JSON parser
# ──────────────────────────────────────────────

def _try_parse_json(content: str) -> list[dict] | None:
    """
    Attempt to parse JSON with fallbacks for common LLM output issues.
    """
    content = content.strip()

    # Remove markdown code fences
    if content.startswith("```"):
        first_nl = content.find("\n")
        if first_nl >= 0:
            content = content[first_nl + 1:]
    if content.endswith("```"):
        content = content[:-3].strip()
    content = content.strip()

    # Direct parse
    try:
        result = json.loads(content)
        if isinstance(result, list):
            return result
        if isinstance(result, dict):
            return [result]
    except json.JSONDecodeError:
        pass

    # Fix trailing commas
    fixed = re.sub(r',\s*([}\]])', r'\1', content)
    try:
        result = json.loads(fixed)
        if isinstance(result, list):
            return result
        if isinstance(result, dict):
            return [result]
    except json.JSONDecodeError:
        pass

    # Fix unquoted property names
    fixed2 = re.sub(r'(?<=[{,])\s*(\w+)\s*:', r' "\1":', fixed)
    try:
        result = json.loads(fixed2)
        if isinstance(result, list):
            return result
        if isinstance(result, dict):
            return [result]
    except json.JSONDecodeError:
        pass

    # Extract JSON array from surrounding text
    bracket_start = content.find('[')
    bracket_end = content.rfind(']')
    if bracket_start >= 0 and bracket_end > bracket_start:
        json_substr = content[bracket_start:bracket_end + 1]
        json_substr_fixed = re.sub(r',\s*([}\]])', r'\1', json_substr)
        try:
            result = json.loads(json_substr_fixed)
            if isinstance(result, list):
                return result
        except json.JSONDecodeError:
            pass

    # Parse individual JSON objects
    objects = []
    depth = 0
    start_idx = None
    for i, ch in enumerate(content):
        if ch == '{':
            if depth == 0:
                start_idx = i
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0 and start_idx is not None:
                obj_str = content[start_idx:i + 1]
                obj_str_fixed = re.sub(r',\s*([}\]])', r'\1', obj_str)
                obj_str_fixed = re.sub(r'(?<=[{,])\s*(\w+)\s*:', r' "\1":', obj_str_fixed)
                try:
                    obj = json.loads(obj_str_fixed)
                    if isinstance(obj, dict):
                        objects.append(obj)
                except json.JSONDecodeError:
                    pass
                start_idx = None

    if objects:
        return objects

    return None


# ──────────────────────────────────────────────
# Helper: check if a value is "empty" / missing
# ──────────────────────────────────────────────

def _is_empty(val) -> bool:
    """Check if a value is empty/missing/null."""
    if val is None:
        return True
    s = str(val).strip()
    return s == "" or s.lower() == "none" or s.lower() == "null"


def _is_ditto_mark(val) -> bool:
    """Check if a value is a ditto mark indicating doubly-classified line."""
    if val is None:
        return False
    s = str(val).strip()
    return s in ('"', '″', '〃', "''", '``', 'ditto', 'do.', 'do')


# Extraction fields — the columns from the PDF table
_OBSERVED_FIELDS = ["wavelength", "intensity", "line_character", "wavenumber"]
_CLASSIFICATION_FIELDS = [
    "lower_level_int", "lower_parity", "lower_j",
    "upper_level_int", "upper_parity", "upper_j",
]


# ──────────────────────────────────────────────
# Node 1: Convert PDF to Images
# ──────────────────────────────────────────────

def convert_pdf_node(state: PipelineState) -> dict:
    """Convert PDF pages to PNG images."""
    if state.get("skip_images"):
        print("⏭️  Skipping image conversion — using existing images")
        images = get_existing_images(
            PAGE_IMAGES_DIR,
            state.get("first_page"),
            state.get("last_page"),
        )
        if not images:
            return {"errors": ["No existing page images found in page_images/"], "page_images": []}
        print(f"   Found {len(images)} existing page images")
        return {"page_images": [str(p) for p in images]}

    images = convert_pdf_to_images(
        pdf_path=Path(state["pdf_path"]),
        first_page=state.get("first_page"),
        last_page=state.get("last_page"),
        dpi=state.get("dpi", PDF_DPI),
    )
    return {"page_images": [str(p) for p in images]}


# ──────────────────────────────────────────────
# Node 2: Extract Tables via Gemini
# ──────────────────────────────────────────────

def extract_tables_node(state: PipelineState) -> dict:
    """Send each page image to Gemini and extract table rows."""
    page_images = state.get("page_images", [])
    if not page_images:
        return {"errors": ["No page images to process"], "raw_extractions": [], "all_rows": []}

    model_name = state.get("model_name", GEMINI_MODEL)
    print(f"\n🤖 Extracting tables using {model_name}...")

    llm = ChatGoogleGenerativeAI(
        model=model_name,
        google_api_key=GOOGLE_API_KEY,
        temperature=1.0,
        # max_output_tokens=8192
    )

    raw_extractions = []
    all_rows = []
    errors = []

    for i, image_path in enumerate(page_images):
        page_name = Path(image_path).stem
        print(f"\n   [{i + 1}/{len(page_images)}] Processing {page_name}...")

        try:
            with open(image_path, "rb") as f:
                image_bytes = f.read()
            image_b64 = base64.b64encode(image_bytes).decode("utf-8")
            image_data_uri = f"data:image/png;base64,{image_b64}"

            messages = [
                SystemMessage(content=get_system_prompt()),
                HumanMessage(content=[
                    {"type": "text", "text": USER_PROMPT},
                    {"type": "image_url", "image_url": image_data_uri},
                ]),
            ]

            response = llm.invoke(messages)
            content = response.content

            if isinstance(content, list):
                text_parts = []
                for part in content:
                    if isinstance(part, dict) and "text" in part:
                        text_parts.append(part["text"])
                    elif isinstance(part, str):
                        text_parts.append(part)
                content = "\n".join(text_parts)

            rows = _try_parse_json(content)

            if rows is not None:
                for row in rows:
                    row["_source_page"] = page_name
                raw_extractions.append({"page": page_name, "row_count": len(rows), "status": "OK"})
                all_rows.extend(rows)
                print(f"   ✅ Extracted {len(rows)} rows from {page_name}")
            else:
                debug_path = PAGE_IMAGES_DIR / f"{page_name}_raw_response.txt"
                with open(str(debug_path), "w", encoding="utf-8") as f:
                    f.write(content)
                raw_extractions.append({"page": page_name, "row_count": 0, "status": "JSON parse error"})
                errors.append(f"{page_name}: JSON parse error — saved to {debug_path.name}")
                print(f"   ❌ JSON parse error on {page_name}")

        except Exception as e:
            raw_extractions.append({"page": page_name, "row_count": 0, "status": f"Error: {e}"})
            errors.append(f"{page_name}: {e}")
            print(f"   ❌ Error on {page_name}: {e}")
            traceback.print_exc()

    print(f"\n📊 Total rows extracted: {len(all_rows)} from {len(page_images)} pages")
    return {"raw_extractions": raw_extractions, "all_rows": all_rows, "errors": errors}


# ──────────────────────────────────────────────
# Node 3: Fix Row Structure
# ──────────────────────────────────────────────

def fix_row_structure_node(state: PipelineState) -> dict:
    """
    Fix structural issues before validation:
      Pass 1: Handle doubly-classified lines (ditto marks, duplicate rows)
      Pass 2: Merge rows with missing wavelength (has wavenumber)
      Pass 3: Merge rows with missing wavenumber (has wavelength)
    All operations are page-local.
    """
    all_rows = state.get("all_rows", [])
    if not all_rows:
        return {"all_rows": []}

    print(f"\n🔧 Fixing row structure ({len(all_rows)} rows)...")

    # Group rows by page, preserving order
    pages = defaultdict(list)
    for row in all_rows:
        page = row.get("_source_page", "unknown")
        pages[page].append(row)

    fixed_rows = []
    stats = {"doubly_classified": 0, "duplicate_rows": 0,
             "wl_merged": 0, "wn_merged": 0,
             "merge_fail_occupied": 0, "merge_fail_distant": 0,
             "rows_deleted": 0}

    for page, page_rows in sorted(pages.items()):
        # ── Pass 1: Doubly-classified lines ──
        _pass1_doubly_classified(page_rows, stats)

        # ── Pass 2: Missing wavelength (has wavenumber) ──
        _pass2_missing_wavelength(page_rows, stats, page)

        # ── Pass 3: Missing wavenumber (has wavelength) ──
        _pass3_missing_wavenumber(page_rows, stats, page)

        # Remove rows marked for deletion (all fields empty after merge)
        for row in page_rows:
            if not row.get("_delete"):
                fixed_rows.append(row)
            else:
                stats["rows_deleted"] += 1

    print(f"\n   ── Row Structure Summary ──")
    print(f"   🔄 Doubly-classified lines detected: {stats['doubly_classified']}")
    print(f"   🔄 Duplicate consecutive rows: {stats['duplicate_rows']}")
    print(f"   🔗 Rows merged (missing wavelength): {stats['wl_merged']}")
    print(f"   🔗 Rows merged (missing wavenumber): {stats['wn_merged']}")
    print(f"   ❌ Merge failed (target occupied): {stats['merge_fail_occupied']}")
    print(f"   ❌ Merge failed (too distant): {stats['merge_fail_distant']}")
    print(f"   🗑️  Rows deleted: {stats['rows_deleted']}")
    print(f"   📊 Rows after fixing: {len(fixed_rows)}")

    return {"all_rows": fixed_rows}


def _pass1_doubly_classified(page_rows: list[dict], stats: dict):
    """
    Detect doubly-classified lines:
    1. Rows with ditto mark in wavelength (", ″, etc.) — Sugar's convention
       for continuation rows of multiply-classified lines.
    2. Consecutive rows with identical wavelength+wavenumber — both have
       observed data, just mark with **.

    For case 1: copy wavelength/intensity/character/wavenumber from previous
    observed row, append ** to character on ALL rows of the group.

    For case 2: append ** to character on all matching rows.
    """
    # Find the "parent" row for each ditto row by scanning backward
    i = 0
    while i < len(page_rows):
        row = page_rows[i]
        wl = row.get("wavelength")

        # Case 1: missing/ditto observed data (it has a classification but no explicit wavelength/wavenumber)
        # We classify it as a continuation row if:
        #  - It HAS a classification AND
        #  - (wavelength is empty/ditto/"0") OR (intensity, character, and wavenumber are ALL empty)
        
        has_class = _row_has_classification(row)
        wl_is_missing = _is_empty(wl) or _is_ditto_mark(wl) or str(wl).strip() == "0"
        
        int_empty = _is_empty(row.get("intensity"))
        char_empty = _is_empty(row.get("line_character"))
        wn_empty = _is_empty(row.get("wavenumber"))
        obs_all_empty = int_empty and char_empty and wn_empty

        # if has_class and (wl_is_missing or obs_all_empty):
        if (wl_is_missing or obs_all_empty):
            # Find the parent row (closest previous row with observed wavelength)
            parent_idx = None
            for j in range(i - 1, -1, -1):
                if not _is_empty(page_rows[j].get("wavelength")) and not _is_ditto_mark(page_rows[j].get("wavelength")):
                    parent_idx = j
                    break

            if parent_idx is not None:
                parent = page_rows[parent_idx]
                # Append ** to parent's character (only once)
                parent_char = str(parent.get("line_character", "") or "")
                if "**" not in parent_char:
                    parent["line_character"] = (parent_char + "**").strip()

                # Copy observed data from parent to ditto row, with **
                row["wavelength"] = parent.get("wavelength")
                row["intensity"] = parent.get("intensity")
                row["wavenumber"] = parent.get("wavenumber")
                char_val = str(parent.get("line_character", "") or "")
                if "**" not in char_val:
                    char_val = (char_val + "**").strip()
                row["line_character"] = char_val

                stats["doubly_classified"] += 1
            i += 1
            continue

        # Case 2: identical consecutive wavelength
        if i > 0 and not _is_empty(wl):
            prev = page_rows[i - 1]
            prev_wl = prev.get("wavelength")
            prev_wn = prev.get("wavenumber")
            cur_wn = row.get("wavenumber")

            # if (not _is_empty(prev_wl) and not _is_empty(cur_wn)
            #         and not _is_empty(prev_wn)):
            if (not _is_empty(prev_wl)):
                try:
                    # if (abs(float(wl) - float(prev_wl)) < 0.0005
                    #         and abs(float(cur_wn) - float(prev_wn)) < 0.005):
                    if (abs(float(wl) - float(prev_wl)) < 0.0005):
                        # Mark both with **
                        for r in [prev, row]:
                            char_val = str(r.get("line_character", "") or "")
                            if "**" not in char_val:
                                r["line_character"] = (char_val + "**").strip()
                        if (not _is_empty(prev_wn)):
                            row["wavenumber"] = prev_wn
                        else:
                            prev["wavenumber"] = cur_wn
                        prev_intens = prev.get("intensity")
                        if (not _is_empty(prev_intens)):
                            row["intensity"] = prev_intens
                        else:
                            prev["intensity"] = row.get("intensity")
                        stats["duplicate_rows"] += 1
                except (ValueError, TypeError):
                    pass

        i += 1


def _pass2_missing_wavelength(page_rows: list[dict], stats: dict, page: str):
    """
    Handle rows with wavenumber but no wavelength.

    1. Compute vacuum wavelength from wavenumber.
    2. Find closest wavelength match on same page (|dwl| ≤ 0.005).
    3. If target fields are free → merge; otherwise → flag.
    """
    for row in page_rows:
        if row.get("_delete"):
            continue

        wl = row.get("wavelength")
        wn = row.get("wavenumber")

        # Only process: has wavenumber, missing wavelength
        if _is_empty(wl) and not _is_empty(wn):
            try:
                wn_val = float(wn)
                if get_config()["wavelength_type"] == "air":
                    computed_wl = wavelength_air_from_wavenumber(wn_val)
                else:
                    computed_wl = wavelength_vac_from_wavenumber(wn_val)
            except (ValueError, TypeError, ZeroDivisionError):
                continue

            # Find closest wavelength match on this page
            best_target = None
            best_diff = float("inf")
            for cand in page_rows:
                if cand is row or cand.get("_delete"):
                    continue
                cand_wl = cand.get("wavelength")
                if _is_empty(cand_wl):
                    continue
                try:
                    diff = abs(float(cand_wl) - computed_wl)
                    if diff < best_diff:
                        best_diff = diff
                        best_target = cand
                except (ValueError, TypeError):
                    continue

            if best_target is None or best_diff > 0.005:
                row["_structure_note"] = (
                    f"Missing wavelength: computed λ={computed_wl:.3f} from σ={wn}, "
                    f"no matching row on {page} within 0.005 Å"
                )
                stats["merge_fail_distant"] += 1
                continue

            # Check if any non-empty field in current row conflicts with target
            conflict = False
            non_empty_fields = []
            for field in _OBSERVED_FIELDS + _CLASSIFICATION_FIELDS:
                val = row.get(field)
                if not _is_empty(val) and not _is_ditto_mark(val):
                    non_empty_fields.append(field)
                    target_val = best_target.get(field)
                    if not _is_empty(target_val):
                        conflict = True
                        break

            if conflict:
                row["_structure_note"] = (
                    f"Missing wavelength: computed λ={computed_wl:.3f}, "
                    f"target row (λ={best_target.get('wavelength')}) has conflicting field '{field}'"
                )
                stats["merge_fail_occupied"] += 1
                continue

            # Merge: move all non-empty values to target
            for field in non_empty_fields:
                best_target[field] = row[field]

            # Track the merge in target's notes
            merge_note = f"Merged from split row (missing λ, had σ={wn})"
            existing_note = best_target.get("_structure_note", "")
            best_target["_structure_note"] = f"{existing_note}; {merge_note}" if existing_note else merge_note

            # Mark current row for deletion
            row["_delete"] = True
            stats["wl_merged"] += 1


def _pass3_missing_wavenumber(page_rows: list[dict], stats: dict, page: str):
    """
    Handle rows with wavelength but no wavenumber.

    1. Compute vacuum wavenumber from vacuum wavelength.
    2. Find closest wavenumber match on same page (|dwn| ≤ 0.2).
    3. If target fields are free → merge; otherwise → flag.
    """
    for row in page_rows:
        if row.get("_delete"):
            continue

        wl = row.get("wavelength")
        wn = row.get("wavenumber")

        # Only process: has wavelength, missing wavenumber
        if not _is_empty(wl) and _is_empty(wn):
            try:
                wl_val = float(wl)
                if get_config()["wavelength_type"] == "air":
                    computed_wn = wavenumber_from_air_wavelength(wl_val)
                else:
                    computed_wn = wavenumber_from_vac_wavelength(wl_val)
            except (ValueError, TypeError, ZeroDivisionError):
                continue

            # Find closest wavenumber match on this page
            best_target = None
            best_diff = float("inf")
            for cand in page_rows:
                if cand is row or cand.get("_delete"):
                    continue
                cand_wn = cand.get("wavenumber")
                if _is_empty(cand_wn):
                    continue
                try:
                    diff = abs(float(cand_wn) - computed_wn)
                    if diff < best_diff:
                        best_diff = diff
                        best_target = cand
                except (ValueError, TypeError):
                    continue

            if best_target is None or best_diff > 0.2:
                row["_structure_note"] = (
                    f"Missing wavenumber: computed σ={computed_wn:.2f} from λ={wl}, "
                    f"no matching row on {page} within 0.2 cm⁻¹"
                )
                stats["merge_fail_distant"] += 1
                continue

            # Check for conflicts
            conflict = False
            non_empty_fields = []
            for field in _OBSERVED_FIELDS + _CLASSIFICATION_FIELDS:
                val = row.get(field)
                if not _is_empty(val) and not _is_ditto_mark(val):
                    non_empty_fields.append(field)
                    target_val = best_target.get(field)
                    if not _is_empty(target_val):
                        conflict = True
                        break

            if conflict:
                row["_structure_note"] = (
                    f"Missing wavenumber: computed σ={computed_wn:.2f}, "
                    f"target row (σ={best_target.get('wavenumber')}) has conflicting field '{field}'"
                )
                stats["merge_fail_occupied"] += 1
                continue

            # Merge
            for field in non_empty_fields:
                best_target[field] = row[field]

            merge_note = f"Merged from split row (missing σ, had λ={wl})"
            existing_note = best_target.get("_structure_note", "")
            best_target["_structure_note"] = f"{existing_note}; {merge_note}" if existing_note else merge_note

            row["_delete"] = True
            stats["wn_merged"] += 1


# ──────────────────────────────────────────────
# Node 4: Validate All Rows
# ──────────────────────────────────────────────

def validate_rows_node(state: PipelineState) -> dict:
    """
    Run full validation on each row:
    1. NIST ASD lookup (with parity/J/int auto-correction)
    2. Selection rule checks (only when both levels found)
    3. Wavelength↔wavenumber validation (only correct if NIST is OK)
    """
    all_rows = state.get("all_rows", [])
    if not all_rows:
        return {"final_rows": [], "errors": ["No rows to validate"]}

    print(f"\n🔬 Validating {len(all_rows)} rows (NIST + selection rules + wavelength)...")

    final_rows = []
    stats = {
        "nist_ok": 0, "nist_ok_sel_fail": 0, "nist_not_found": 0,
        "nist_ritz_mismatch": 0, "unclassified": 0,
        "wl_ok": 0, "wl_corrected": 0, "wl_error": 0, "wl_suspect": 0,
        "parity_corrected": 0, "j_corrected": 0, "int_corrected": 0,
    }

    for row in all_rows:
        # User requested: If intensity is 0, NEVER assign a classification
        intensity_val = str(row.get("intensity", "")).strip()
        if intensity_val == "0" and _row_has_classification(row):
            for field in ["lower_level_int", "lower_parity", "lower_j",
                          "upper_level_int", "upper_parity", "upper_j"]:
                row[field] = None
            row["_structure_note"] = "Classification deleted (lines with intensity 0 must not be classified)"

        validated = validate_full_row(row)

        nist_status = validated.get("nist_status", "")
        wl_status = validated.get("wl_wn_status", "")

        if nist_status == "OK":
            stats["nist_ok"] += 1
        elif nist_status == "RITZ_OK_SELECTION_FAIL":
            stats["nist_ok_sel_fail"] += 1
        elif nist_status == "LEVEL_NOT_FOUND":
            stats["nist_not_found"] += 1
        elif nist_status == "RITZ_MISMATCH":
            stats["nist_ritz_mismatch"] += 1
        elif nist_status == "UNCLASSIFIED":
            stats["unclassified"] += 1

        if wl_status == "OK":
            stats["wl_ok"] += 1
        elif wl_status in ("CORRECTED_WAVELENGTH", "CORRECTED_WAVENUMBER"):
            stats["wl_corrected"] += 1
        elif wl_status == "MISMATCH_SUSPECT_MISALIGN":
            stats["wl_suspect"] += 1
        elif wl_status == "ERROR":
            stats["wl_error"] += 1

        if validated.get("lower_parity_original") or validated.get("upper_parity_original"):
            stats["parity_corrected"] += 1
        if validated.get("lower_j_original") or validated.get("upper_j_original"):
            stats["j_corrected"] += 1
        if validated.get("lower_level_int_original") or validated.get("upper_level_int_original"):
            stats["int_corrected"] += 1

        final_rows.append(validated)

    print(f"\n   ── NIST ASD Validation ──")
    print(f"   ✅ OK: {stats['nist_ok']}")
    print(f"   ⚠️  Ritz OK but selection rule fail: {stats['nist_ok_sel_fail']}")
    print(f"   ❌ Level not found: {stats['nist_not_found']}")
    print(f"   ❌ Ritz mismatch: {stats['nist_ritz_mismatch']}")
    print(f"   ➖ Unclassified: {stats['unclassified']}")

    print(f"\n   ── Auto-Corrections Applied ──")
    print(f"   🔧 Parity corrected: {stats['parity_corrected']}")
    print(f"   🔧 J-value corrected: {stats['j_corrected']}")
    print(f"   🔧 Energy level int corrected: {stats['int_corrected']}")

    print(f"\n   ── Wavelength/Wavenumber ──")
    print(f"   ✅ OK: {stats['wl_ok']}")
    print(f"   🔧 Corrected (NIST-confirmed): {stats['wl_corrected']}")
    print(f"   🟡 Suspect misalignment (not corrected): {stats['wl_suspect']}")
    print(f"   ❌ Error: {stats['wl_error']}")

    return {"final_rows": final_rows}


# ──────────────────────────────────────────────
# Node 5: Page-Level Misalignment Correction
# ──────────────────────────────────────────────

def _row_has_classification(row: dict) -> bool:
    return row.get("lower_level_int") is not None and row.get("upper_level_int") is not None


def _row_is_unclassified(row: dict) -> bool:
    return row.get("lower_level_int") is None and row.get("upper_level_int") is None


def _run_one_misalignment_pass(final_rows, pages, iteration):
    """
    Single pass of misalignment correction across all pages.

    For EVERY classified row with known NIST levels, find the row on the
    same page (within ±10 positions) whose observed wavenumber is closest
    to the Ritz wavenumber. If the closest is a different row and is
    unclassified, move the classification there.

    Returns (moves_made, flags_made) counts.
    """
    from tools import validate_wavelength_wavenumber

    total_moved = 0
    total_flagged = 0

    SEARCH_WINDOW = 10  # Look ±10 rows from current position

    for page, indices in sorted(pages.items()):
        # Build list of ALL classified rows with known NIST levels
        ritz_candidates = []
        for pos, idx in enumerate(indices):
            row = final_rows[idx]
            if not _row_has_classification(row):
                continue
            lower_exact = row.get("lower_exact_nist")
            upper_exact = row.get("upper_exact_nist")
            if lower_exact is None or upper_exact is None:
                continue
            ritz_wn = upper_exact - lower_exact
            ritz_candidates.append((idx, pos, ritz_wn))

        if not ritz_candidates:
            continue

        moved_this_page = set()

        for err_idx, err_pos, ritz_wn in ritz_candidates:
            if err_idx in moved_this_page:
                continue
            err_row = final_rows[err_idx]
            if not _row_has_classification(err_row):
                continue  # May have been cleared by an earlier move this pass

            # Search ±SEARCH_WINDOW rows from this row's position on the page
            search_start = max(0, err_pos - SEARCH_WINDOW)
            search_end = min(len(indices), err_pos + SEARCH_WINDOW + 1)

            best_match_idx = None
            best_diff = float("inf")

            for search_pos in range(search_start, search_end):
                cand_idx = indices[search_pos]
                cand_row = final_rows[cand_idx]
                cand_wn = cand_row.get("wavenumber")
                if cand_wn is None:
                    continue
                try:
                    diff = abs(float(cand_wn) - ritz_wn)
                    if diff < best_diff:
                        best_diff = diff
                        best_match_idx = cand_idx
                except (ValueError, TypeError):
                    continue

            if best_match_idx is None or best_diff > get_config()["nist_ritz_tolerance_cm1"]:
                continue

            # If best match is current row, classification is in the right place
            if best_match_idx == err_idx:
                continue

            target_row = final_rows[best_match_idx]

            # Check if target row already has the EXACT SAME classification
            is_exact_duplicate = False
            if _row_has_classification(target_row):
                if (target_row.get("lower_exact_nist") == err_row.get("lower_exact_nist") and
                    target_row.get("upper_exact_nist") == err_row.get("upper_exact_nist")):
                    is_exact_duplicate = True

            if _row_is_unclassified(target_row) or is_exact_duplicate:
                if is_exact_duplicate:
                    # Target already has this classification. Just delete from current row.
                    move_note = f"[iter {iteration}] Exact duplicate classification deleted (merged into σ_obs={float(target_row.get('wavenumber', 0)):.2f})"
                    
                    # Clear classification from error row
                    for field in ["lower_level_int", "lower_parity", "lower_j",
                                  "upper_level_int", "upper_parity", "upper_j",
                                  "lower_exact_nist", "upper_exact_nist"]:
                        err_row[field] = None
                    err_row["nist_status"] = "UNCLASSIFIED"
                    err_note = err_row.get("nist_note", "")
                    err_row["nist_note"] = f"{err_note}; {move_note}" if err_note else move_note
                    err_row["wl_wn_status"] = "OK"
                    err_row["wl_wn_note"] = ""

                    moved_this_page.add(err_idx)
                    total_moved += 1
                    print(f"      [{page} iter {iteration}] Deleted exact duplicate classification (Ritz σ={ritz_wn:.2f})")
                    
                else:
                    # Move classification from err_row to target_row
                    classification_fields = [
                        "lower_level_int", "lower_parity", "lower_j",
                        "upper_level_int", "upper_parity", "upper_j",
                        "lower_exact_nist", "upper_exact_nist",
                        "nist_status",
                        "lower_parity_original", "upper_parity_original",
                        "lower_j_original", "upper_j_original",
                        "lower_level_int_original", "upper_level_int_original",
                    ]
                    for field in classification_fields:
                        if field in err_row and err_row[field] is not None:
                            target_row[field] = err_row[field]

                    # Re-validate wavelength/wavenumber on target
                    target_wl = target_row.get("wavelength")
                    target_wn = target_row.get("wavenumber")
                    if target_wl is not None and target_wn is not None:
                        wl_check = validate_wavelength_wavenumber(float(target_wl), float(target_wn))
                        target_row["wl_wn_status"] = wl_check["status"]
                        target_row["wl_wn_note"] = wl_check.get("note", "")

                    if best_diff <= get_config()["nist_ritz_tolerance_cm1"]:
                        target_row["nist_status"] = "OK"
                    move_note = f"[iter {iteration}] Classification moved from misaligned row (Ritz σ={ritz_wn:.2f}, Δ={best_diff:.4f})"
                    existing = target_row.get("nist_note", "")
                    target_row["nist_note"] = f"{existing}; {move_note}" if existing else move_note

                    # Clear classification from error row
                    for field in ["lower_level_int", "lower_parity", "lower_j",
                                  "upper_level_int", "upper_parity", "upper_j",
                                  "lower_exact_nist", "upper_exact_nist"]:
                        err_row[field] = None
                    err_row["nist_status"] = "UNCLASSIFIED"
                    err_note = err_row.get("nist_note", "")
                    clear_note = f"[iter {iteration}] Classification moved to row with σ={float(target_row.get('wavenumber', 0)):.2f}"
                    err_row["nist_note"] = f"{err_note}; {clear_note}" if err_note else clear_note
                    err_row["wl_wn_status"] = "OK"
                    err_row["wl_wn_note"] = ""

                    moved_this_page.add(err_idx)
                    total_moved += 1

                    print(f"      [{page} iter {iteration}] Moved classification "
                          f"(Ritz σ={ritz_wn:.2f}) to row with σ_obs={float(target_row.get('wavenumber', 0)):.2f}")
            else:
                # Target classified — will be re-evaluated in next iteration
                # (the target's classification might itself be moved away)
                total_flagged += 1

    return total_moved, total_flagged


def fix_misalignment_node(state: PipelineState) -> dict:
    """
    Iterating page-level misalignment correction.

    Checks ALL classified rows on every page. For each classified row whose
    Ritz wavenumber doesn't match its own observed wavenumber:
      - Find the row on the same page with closest observed wavenumber
      - If unclassified → move classification there
      - If classified → skip (it may be freed in a later iteration)

    Iterates until no more moves are made (convergence), up to max 10 iterations.
    This handles cascading moves: when row A's classification is moved away,
    row A becomes free for another classification in the next iteration.
    """
    final_rows = state.get("final_rows", [])
    if not final_rows:
        return {"final_rows": []}

    # Group by page
    pages = defaultdict(list)
    for idx, row in enumerate(final_rows):
        page = row.get("_source_page", "unknown")
        pages[page].append(idx)

    print(f"\n🔧 Page-level misalignment correction (iterating until convergence)...")

    grand_total_moved = 0
    max_iterations = 10

    for iteration in range(1, max_iterations + 1):
        moved, flagged = _run_one_misalignment_pass(final_rows, pages, iteration)
        grand_total_moved += moved

        print(f"   Iteration {iteration}: moved {moved}, blocked {flagged}")

        if moved == 0:
            break

    # Final pass: flag remaining unresolvable cases
    remaining_flagged = 0
    for page, indices in sorted(pages.items()):
        for idx in indices:
            row = final_rows[idx]
            if not _row_has_classification(row):
                continue
            lower_exact = row.get("lower_exact_nist")
            upper_exact = row.get("upper_exact_nist")
            if lower_exact is None or upper_exact is None:
                continue
            ritz_wn = upper_exact - lower_exact
            obs_wn = row.get("wavenumber")
            if obs_wn is None:
                continue
            try:
                obs_diff = abs(float(obs_wn) - ritz_wn)
            except (ValueError, TypeError):
                continue
            if obs_diff > get_config()["nist_ritz_tolerance_cm1"]:
                # Still misaligned after all iterations
                # Find where it should be
                best_idx = None
                best_diff = float("inf")
                for cand_idx in indices:
                    cand_wn = final_rows[cand_idx].get("wavenumber")
                    if cand_wn is None:
                        continue
                    try:
                        d = abs(float(cand_wn) - ritz_wn)
                        if d < best_diff:
                            best_diff = d
                            best_idx = cand_idx
                    except (ValueError, TypeError):
                        continue
                if best_idx is not None and best_diff <= get_config()["nist_ritz_tolerance_cm1"] and best_idx != idx:
                    target = final_rows[best_idx]
                    flag = (f"Misalignment unresolved: Ritz σ={ritz_wn:.2f} matches row with "
                            f"σ_obs={float(target.get('wavenumber', 0)):.2f} (Δ={best_diff:.4f}), "
                            f"but target is classified — needs manual review")
                    existing = row.get("nist_note", "")
                    # Don't duplicate flags
                    if "Misalignment unresolved" not in existing:
                        row["nist_note"] = f"{existing}; {flag}" if existing else flag
                    remaining_flagged += 1

    print(f"\n   ── Misalignment Summary ──")
    print(f"   🔧 Classifications moved (total across iterations): {grand_total_moved}")
    print(f"   ❌ Unresolvable after all iterations: {remaining_flagged}")

    return {"final_rows": final_rows}


# ──────────────────────────────────────────────
# Node 6: Compile Output
# ──────────────────────────────────────────────

def compile_output_node(state: PipelineState) -> dict:
    """Sort rows by wavelength (descending) and write CSV + Excel."""
    final_rows = state.get("final_rows", [])
    if not final_rows:
        return {"errors": ["No rows to output"]}

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"\n📁 Compiling output for {len(final_rows)} rows...")

    output_records = []
    for row in final_rows:
        notes = []

        # Structure notes (from fix_row_structure)
        struct_note = row.get("_structure_note", "")
        if struct_note:
            notes.append(struct_note)

        # Wavelength/wavenumber notes
        wl_note = row.get("wl_wn_note", "")
        if wl_note:
            notes.append(wl_note)

        # NIST notes
        nist_note = row.get("nist_note", "")
        if nist_note:
            notes.append(nist_note)

        # Correction tracking
        corrections = []
        if row.get("wavelength_original"):
            corrections.append(f"λ_orig={row['wavelength_original']}")
        if row.get("wavenumber_original"):
            corrections.append(f"σ_orig={row['wavenumber_original']}")
        if row.get("lower_parity_original"):
            corrections.append(f"lower_par_orig={row['lower_parity_original']}")
        if row.get("upper_parity_original"):
            corrections.append(f"upper_par_orig={row['upper_parity_original']}")
        if row.get("lower_j_original"):
            corrections.append(f"lower_J_orig={row['lower_j_original']}")
        if row.get("upper_j_original"):
            corrections.append(f"upper_J_orig={row['upper_j_original']}")
        if row.get("lower_level_int_original"):
            corrections.append(f"lower_int_orig={row['lower_level_int_original']}")
        if row.get("upper_level_int_original"):
            corrections.append(f"upper_int_orig={row['upper_level_int_original']}")
        if corrections:
            notes.append("Corrections: " + ", ".join(corrections))

        # Page number
        source_page = row.get("_source_page", "")
        page_num = ""
        if source_page:
            match = re.search(r'(\d+)', source_page)
            if match:
                page_num = int(match.group(1))

        record = {
            "Wavelength (Å)": row.get("wavelength"),
            "Intensity": row.get("intensity"),
            "Line Character": row.get("line_character", ""),
            "Wavenumber (cm⁻¹)": row.get("wavenumber"),
            "Lower Level Int": row.get("lower_level_int"),
            "Lower Parity": row.get("lower_parity"),
            "Lower J": row.get("lower_j"),
            "Upper Level Int": row.get("upper_level_int"),
            "Upper Parity": row.get("upper_parity"),
            "Upper J": row.get("upper_j"),
            "Lower Level NIST (cm⁻¹)": row.get("lower_exact_nist"),
            "Upper Level NIST (cm⁻¹)": row.get("upper_exact_nist"),
            "Page": page_num,
            "Notes": "; ".join(notes) if notes else "",
        }
        output_records.append(record)

    df = pd.DataFrame(output_records)

    # Replaced: NO SORTING by wavelength! We must preserve Author's original order
    # (Left column top-to-bottom, then right column)
    # df["_wl_sort"] = pd.to_numeric(df["Wavelength (Å)"], errors="coerce")
    # df = df.sort_values("_wl_sort", ascending=False).drop(columns=["_wl_sort"])
    # df = df.reset_index(drop=True)

    # ── Save CSV with ="5/2" formatting for J columns ──
    csv_path = OUTPUT_DIR / get_config()["csv_name"]
    df_csv = df.copy()
    for j_col in ["Lower J", "Upper J"]:
        df_csv[j_col] = df_csv[j_col].apply(
            lambda v: f'="{v}"' if pd.notna(v) and v != "" and "/" in str(v) else v
        )
    df_csv.to_csv(str(csv_path), index=False, encoding="utf-8-sig")
    print(f"   ✅ CSV saved: {csv_path}")

    # ── Save Excel ──
    excel_path = OUTPUT_DIR / get_config()["excel_name"]
    sheet_name = get_config()["excel_sheet_name"]
    with pd.ExcelWriter(str(excel_path), engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)
        ws = writer.sheets[sheet_name]
        for col_idx, col_name in enumerate(df.columns, 1):
            max_len = max(
                len(str(col_name)),
                df[col_name].astype(str).str.len().max() if len(df) > 0 else 0,
            )
            ws.column_dimensions[ws.cell(row=1, column=col_idx).column_letter].width = min(max_len + 2, 60)

    print(f"   ✅ Excel saved: {excel_path}")

    total = len(df)
    classified = df["Lower Level Int"].notna().sum()
    notes_count = (df["Notes"] != "").sum()
    print(f"\n📊 Summary:")
    print(f"   Total rows: {total}")
    print(f"   Classified: {classified}")
    print(f"   Unclassified: {total - classified}")
    print(f"   Rows with notes/warnings: {notes_count}")

    return {"output_csv": str(csv_path), "output_excel": str(excel_path)}


# ──────────────────────────────────────────────
# Build the Graph
# ──────────────────────────────────────────────

def build_pipeline() -> StateGraph:
    """Build and compile the LangGraph pipeline."""
    workflow = StateGraph(PipelineState)

    workflow.add_node("convert_pdf", convert_pdf_node)
    workflow.add_node("extract_tables", extract_tables_node)
    workflow.add_node("fix_row_structure", fix_row_structure_node)
    workflow.add_node("validate_rows", validate_rows_node)
    workflow.add_node("fix_misalignment", fix_misalignment_node)
    workflow.add_node("compile_output", compile_output_node)

    workflow.add_edge(START, "convert_pdf")
    workflow.add_edge("convert_pdf", "extract_tables")
    workflow.add_edge("extract_tables", "fix_row_structure")
    workflow.add_edge("fix_row_structure", "validate_rows")
    workflow.add_edge("validate_rows", "fix_misalignment")
    workflow.add_edge("fix_misalignment", "compile_output")
    workflow.add_edge("compile_output", END)

    return workflow.compile()


def run_pipeline(
    pdf_path: str | None = None,
    first_page: int | None = None,
    last_page: int | None = None,
    dpi: int | None = None,
    model_name: str | None = None,
    skip_images: bool = False,
) -> dict:
    """Run the complete extraction pipeline."""
    pipeline = build_pipeline()

    if pdf_path is None:
        pdf_path = str(get_config()["pdf_path"])
    if dpi is None:
        dpi = PDF_DPI
    if model_name is None:
        model_name = GEMINI_MODEL

    initial_state: PipelineState = {
        "pdf_path": pdf_path,
        "first_page": first_page,
        "last_page": last_page,
        "dpi": dpi,
        "model_name": model_name,
        "skip_images": skip_images,
        "page_images": [],
        "raw_extractions": [],
        "all_rows": [],
        "final_rows": [],
        "errors": [],
        "output_csv": "",
        "output_excel": "",
    }

    print("=" * 60)
    print("  Pr III Table Extraction Pipeline (v4)")
    print("=" * 60)
    print(f"  Model: {model_name}")
    print(f"  Pages: {first_page or 'first'} – {last_page or 'last'}")
    print(f"  DPI: {dpi}")
    print(f"  Skip images: {skip_images}")
    print("=" * 60)

    result = pipeline.invoke(initial_state)

    print("\n" + "=" * 60)
    print("  Pipeline Complete!")
    print("=" * 60)
    if result.get("output_csv"):
        print(f"  CSV:   {result['output_csv']}")
    if result.get("output_excel"):
        print(f"  Excel: {result['output_excel']}")
    if result.get("errors"):
        print(f"  ⚠️  Errors: {len(result['errors'])}")
        for err in result["errors"][:10]:
            print(f"     • {err}")
    print("=" * 60)

    return result
