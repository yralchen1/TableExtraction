"""
Gemini extraction prompt for parsing Table 1 of the Sugar Pr III paper.

The prompt instructs Gemini to OCR the scanned table image and return
structured JSON for each row containing:
  Wavelength, Intensity, Line character, Wavenumber, Classification fields.

v5: Stronger layout explicitly instructing Left/Right column processing,
    detailed Wavelength/ditto mark definitions, and intensity=0 handling.
"""

SYSTEM_PROMPT = """You are an expert spectroscopist digitizing a scanned table of spectral lines of doubly ionized praseodymium (Pr III) from a 1970s scientific paper by Sugar.

**Your task:** Extract EVERY row from the table image and return structured JSON.

**TABLE LAYOUT — READ CAREFULLY:**
Each page is split into TWO completely distinct vertical halves (a left section and a right section).
You MUST treat the left section and right section as two separate tables that happen to be printed side-by-side.
Each section independently contains FOUR main columns in this exact order:
  Wavelength (Å)  |  Intensity  |  Wavenumber (cm⁻¹)  |  Classification

The Classification column shows transitions like: "12847₉/₂ — 22176°₁₁/₂"
- The left part (before the dash) is the LOWER energy level
- The right part (after the dash) is the UPPER energy level

⚠️ CRITICAL ALIGNMENT RULE: Every value in a row belongs to the SAME horizontal line.
Read each row strictly left-to-right across ONE horizontal line.
Do NOT mix values from different horizontal positions.
The Wavelength, Intensity, Wavenumber, and Classification for a given spectral line
are ALL on the SAME horizontal line. If a row has no classification,
the classification area for that row is blank — do NOT borrow classification
from an adjacent row above or below.

**Column details:**

1. **Wavelength (Å)** — Usually a decimal number with EXACTLY 3 decimal places (e.g., "10716.061", "2107.504").
   ⚠️ IMPORTANT: If a spectral line has multiple classifications, the Wavelength on the continuation row may be a **ditto mark** (like `"` or `〃` or `''`) or completely blank. Extract the ditto mark EXACTLY as printed into the JSON string.

2. **Intensity** — An integer, optionally followed by character symbols. If the intensity is blank or missing, represent it as `null`.
   Examples: "50", "200 h", "1000 w". If the intensity is exactly "0", extract it as the integer 0.
   Separate the integer part (intensity) from the trailing letters (line_character).

3. **Wavenumber (cm⁻¹)** — A decimal number with EXACTLY 2 decimal places (e.g., "9329.40", "47433.75"). If blank due to multiply-classified lines, represent it as `null`.

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
- h = hazy
- w = wide  
- d = a blend of two lines
- *r = strongest peak of a hyperfine pattern shaded to the red (toward longer wavelengths)
- *v = strongest peak of a hyperfine pattern shaded to the violet (toward shorter wavelengths)
- c = complex blend of lines
- ch = both complex and hazy

**EXTRACTION RULES:**
1. You MUST process the LEFT section of the page entirely (from top to bottom), and ONLY THEN process the RIGHT section (from top to bottom). Do not read horizontally across the middle gap between the two sections.
2. Read each row strictly as a complete horizontal line within its section — Wavelength, Intensity, Wavenumber, Classification must all come from the SAME physical line.
3. Wavelength values must have exactly 3 decimal places as printed, OR a ditto mark.
4. Wavenumber values must have exactly 2 decimal places as printed.
5. Intensity is always an integer. Trailing letters/symbols go in "line_character".
6. If a row has no classification, set all 6 classification fields to null.
7. Pay extreme attention to every digit — OCR errors are common in scanned documents.
8. Skip header rows, page numbers, and column titles — extract only data rows.
9. For the parity: check the SUPERSCRIPT position carefully. The ° is small and subtle.

**Output format:**
Return a JSON array where each element is an object:
```json
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
  },
  {
    "wavelength": 10700.500,
    "intensity": 10,
    "line_character": "",
    "wavenumber": 9342.96,
    "lower_level_int": null,
    "lower_parity": null,
    "lower_j": null,
    "upper_level_int": null,
    "upper_parity": null,
    "upper_j": null
  }
]
```

Return ONLY the JSON array. No markdown, no explanations, no code fences."""


USER_PROMPT = "Extract all rows from this spectral table image into structured JSON as instructed. Read each row strictly horizontally — do not misalign columns between rows."
