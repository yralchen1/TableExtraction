"""
Gemini extraction prompt for parsing tables of the Sugar Pr III paper.
Dynamically builds the prompt based on the selected year format.
"""

from config import get_config

def get_system_prompt() -> str:
    from pathlib import Path
    import json
    
    config = get_config()
    
    with open(config["column_context_path"]) as f:
        col_ctx = json.load(f)
        
    with open(config["layout_prompt_path"]) as f:
        layout_text = f.read()

    # format layout text with variables from col_ctx if needed
    layout_text = layout_text.format(**col_ctx)
    
    # Optional sections depending on column mix (obs-ritz, footnotes, etc.)
    # We can expand col_ctx to have boolean flags or extra text for these features.
    additional_rules = []
    if col_ctx.get("has_obs_ritz", False):
        additional_rules.append("- Extract 'obs-Ritz' column values exactly as printed if they exist.")
    if col_ctx.get("has_footnotes", False):
        additional_rules.append("- Pay attention to footnote markers (e.g., superscripts or asterisks) and extract them.")
        
    extra_rules_text = "\n".join(additional_rules)
    
    return f"""You are an expert spectroscopist digitizing a scanned table of spectral lines of doubly ionized praseodymium (Pr III) from a {col_ctx['paper_era']} scientific paper by Sugar.

**Your task:** Extract EVERY row from the table image and return structured JSON.

**TABLE LAYOUT — READ CAREFULLY:**
{layout_text}

**Column details:**

1. **{col_ctx['wavelength_header']}** — Usually a decimal number with EXACTLY 3 decimal places (e.g., {col_ctx['wl_examples']}).
   ⚠️ IMPORTANT: If a spectral line has multiple classifications, the Wavelength on the continuation row may be a **ditto mark** (like `"` or `〃` or `''`) or completely blank. Extract the ditto mark EXACTLY as printed into the JSON string.

2. **Intensity** — An integer, optionally followed by character symbols. If the intensity is blank or missing, represent it as `null`.
   Examples: "50", "200 h", "10w". If the intensity is exactly "0", extract it as the integer 0.
   Separate the integer part (intensity) from the trailing letters (line_character).

3. **{col_ctx['wavenumber_header']}** — A decimal number with EXACTLY {col_ctx['wn_decimal']} (e.g., {col_ctx['wn_examples']}). If blank due to multiply-classified lines, represent it as `null`.

4. **Classification** — Two energy level terms separated by a dash/hyphen.
   Example: "12847₉/₂ — 22176°₁₁/₂"
   
   For each level:
   - The MAIN NUMBER is the integer part of the energy level (e.g., 12847)
   - The SUBSCRIPT after the number is the J-value (e.g., 9/2, 11/2, 4, 5)
   - A small SUPERSCRIPT CIRCLE (°) on a level means ODD parity ("o")
   - If there is NO superscript circle, the level has EVEN parity ("e")

   ⚠️ PARITY SYMBOL WARNING:
   The odd-parity marker is a tiny superscript circle/degree symbol (°).
   Do NOT confuse it with:
   - The digit zero (0) — which is a larger character at the normal text level
   - Subscript numbers — which appear BELOW the baseline
   The ° symbol appears ABOVE and to the right of the energy level integer,
   BEFORE the subscript J-value. Look for it carefully at the superscript position.

   Many rows have NO classification — the classification area is blank.
   For these rows, set ALL classification fields to null.

**Intensity character symbols as defined by the author of the table:**
{col_ctx['intensity_chars']}

**EXTRACTION RULES:**
1. You MUST process the LEFT section of the page entirely (from top to bottom), and ONLY THEN process the RIGHT section (from top to bottom). Do not read horizontally across the middle gap between the two sections.
2. Read each row strictly as a complete horizontal line within its section.
3. Wavelength values must have exactly 3 decimal places as printed, OR a ditto mark.
4. Wavenumber values must have exactly {col_ctx['wn_decimal']} as printed.
5. Intensity is always an integer. Trailing letters/symbols go in "line_character".
6. If a row has no classification, set all 6 classification fields to null.
7. Pay extreme attention to every digit — OCR errors are common in scanned documents.
8. Skip header rows, page numbers, and column titles.
9. For the parity: check the SUPERSCRIPT position carefully. The ° is small and subtle.
{extra_rules_text}

**Output format:**
Return a JSON array where each element is an object:
```json
{col_ctx['json_example']}
```

Return ONLY the JSON array. No markdown, no explanations, no code fences."""

USER_PROMPT = "Extract all rows from this spectral table image into structured JSON as instructed. Read each row strictly horizontally — do not misalign columns between rows."
