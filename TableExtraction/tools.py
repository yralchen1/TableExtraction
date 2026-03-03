"""
Validation tools for the Pr III Table Extraction Pipeline.

v2: Added selection rule checks, NIST-guided parity/J correction,
    energy level integer typo correction, and smarter wavelength validation.

Contains:
1. Edlen's 1966 formula for air↔vacuum wavelength conversion
2. Wavelength ↔ Wavenumber cross-validation with typo detection
3. Selection rule validation (parity change, ΔJ ≤ 1)
4. NIST ASD energy level lookup with auto-correction
5. Full row validation with proper ordering
"""

import math
from pathlib import Path
from functools import lru_cache

import openpyxl

from config import get_config
import astools


# ──────────────────────────────────────────────
# Edlen's 1966 formula  (air → vacuum)
# ──────────────────────────────────────────────

def edlen_n_air(sigma_cm1: float) -> float:
    """
    Compute the refractive index of standard air using Edlen's 1966 formula.

    The formula (Edlen, Metrologia 2, 71, 1966) for the refractive index of
    standard air (15°C, 760 mmHg, 0.03% CO₂) as a function of wavenumber:

        (n - 1) × 10⁸ = 8342.13 + 2406030 / (130 - σ²) + 15997 / (38.9 - σ²)

    where σ is in μm⁻¹ (= cm⁻¹ × 10⁻⁴).
    """
    sigma_um = sigma_cm1 * 1e-4
    sigma_sq = sigma_um ** 2
    n_minus_1_times_1e8 = (
        8342.13
        + 2406030.0 / (130.0 - sigma_sq)
        + 15997.0 / (38.9 - sigma_sq)
    )
    return 1.0 + n_minus_1_times_1e8 * 1e-8


def wavelength_air_to_vacuum(wavelength_air_angstrom: float) -> float:
    """Convert air wavelength (Å) to vacuum wavelength (Å) using Edlen's formula."""
    lambda_air_cm = wavelength_air_angstrom * 1e-8
    sigma_approx = 1.0 / lambda_air_cm

    for _ in range(3):
        n = edlen_n_air(sigma_approx)
        lambda_vac_cm = lambda_air_cm * n
        sigma_approx = 1.0 / lambda_vac_cm

    return lambda_vac_cm * 1e8


def wavenumber_from_air_wavelength(wavelength_air_angstrom: float) -> float:
    """Compute vacuum wavenumber (cm⁻¹) from air wavelength (Å)."""
    lambda_vac_angstrom = wavelength_air_to_vacuum(wavelength_air_angstrom)
    lambda_vac_cm = lambda_vac_angstrom * 1e-8
    return 1.0 / lambda_vac_cm


def wavenumber_from_vac_wavelength(wavelength_vac_angstrom: float) -> float:
    """Compute vacuum wavenumber (cm⁻¹) from vacuum wavelength (Å)."""
    if wavelength_vac_angstrom <= 0:
        raise ValueError(f"Wavelength must be positive, got {wavelength_vac_angstrom}")
    return 1.0e8 / wavelength_vac_angstrom


def wavelength_vac_from_wavenumber(sigma_cm1: float) -> float:
    """
    Compute vacuum wavelength (Å) from vacuum wavenumber (cm⁻¹).
    """
    if sigma_cm1 <= 0:
        raise ValueError(f"Wavenumber must be positive, got {sigma_cm1}")
    return 1.0e8 / sigma_cm1


def wavelength_air_from_wavenumber(sigma_cm1: float) -> float:
    """
    Compute air wavelength (Å) from vacuum wavenumber (cm⁻¹).

    Direct computation — no iteration needed because Edlen's formula
    gives n as a function of sigma (vacuum wavenumber) directly:
        lambda_vac = 1 / sigma
        lambda_air = lambda_vac / n_air(sigma)
    """
    if sigma_cm1 <= 0:
        raise ValueError(f"Wavenumber must be positive, got {sigma_cm1}")
    lambda_vac_cm = 1.0 / sigma_cm1
    n = edlen_n_air(sigma_cm1)
    lambda_air_cm = lambda_vac_cm / n
    return lambda_air_cm * 1e8


# ──────────────────────────────────────────────
# Typo Detection / Correction helpers
# ──────────────────────────────────────────────

def _try_single_digit_typo(value_str: str, compute_target, target: float, tolerance: float) -> str | None:
    """Try replacing each digit with all other digits to find a match."""
    for i, ch in enumerate(value_str):
        if ch.isdigit():
            for d in "0123456789":
                if d == ch:
                    continue
                candidate = value_str[:i] + d + value_str[i + 1:]
                try:
                    candidate_val = float(candidate)
                    computed = compute_target(candidate_val)
                    if abs(computed - target) <= tolerance:
                        return candidate
                except (ValueError, ZeroDivisionError):
                    continue
    return None


def _try_swap_adjacent_digits(value_str: str, compute_target, target: float, tolerance: float) -> str | None:
    """Try swapping each pair of adjacent digits."""
    chars = list(value_str)
    for i in range(len(chars) - 1):
        if chars[i].isdigit() and chars[i + 1].isdigit() and chars[i] != chars[i + 1]:
            swapped = chars.copy()
            swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
            candidate = "".join(swapped)
            try:
                candidate_val = float(candidate)
                computed = compute_target(candidate_val)
                if abs(computed - target) <= tolerance:
                    return candidate
            except (ValueError, ZeroDivisionError):
                continue
    return None


# ──────────────────────────────────────────────
# Wavelength ↔ Wavenumber validation
# ──────────────────────────────────────────────

def validate_wavelength_wavenumber(
    wavelength: float,
    wavenumber: float,
    tolerance: float = None,
) -> dict:
    """
    Validate that wavelength and wavenumber are consistent via Edlen's formula.
    If they disagree, attempt typo correction on each value.
    """
    if tolerance is None:
        tolerance = get_config()["wavenumber_tolerance_cm1"]

    try:
        wl_type = get_config()["wavelength_type"]
        if wl_type == "air":
            method_name = get_config().get("air_vac_conversion_method", "Edlen1966")
            if method_name == "Edlen1966":
                computed_sigma = wavenumber_from_air_wavelength(wavelength)
            elif method_name == "Edlen1953":
                computed_sigma = 1.0 / (astools.LvacE53(wavelength) * 1e-8)
            elif method_name == "Peck&Reeder1973":
                computed_sigma = 1.0 / (astools.Lvac(wavelength) * 1e-8)
            elif method_name == "Meggers&Peters1919":
                computed_sigma = 1.0 / (astools.LvacMP(wavelength) * 1e-8)
            else:
                computed_sigma = wavenumber_from_air_wavelength(wavelength)
        else:
            computed_sigma = wavenumber_from_vac_wavelength(wavelength)
    except (ZeroDivisionError, ValueError):
        return {
            "status": "ERROR",
            "computed_wavenumber": None,
            "difference": None,
            "corrected_wavelength": None,
            "corrected_wavenumber": None,
            "note": f"Cannot compute wavenumber from wavelength {wavelength}",
        }

    diff = abs(computed_sigma - wavenumber)

    if diff <= tolerance:
        return {
            "status": "OK",
            "computed_wavenumber": round(computed_sigma, 2),
            "difference": round(diff, 4),
            "corrected_wavelength": None,
            "corrected_wavenumber": None,
            "note": "",
        }

    # --- Attempt typo correction ---
    wl_str = f"{wavelength:.3f}"
    wn_str = f"{wavenumber:.2f}"

    wl_type = get_config()["wavelength_type"]
    if wl_type == "air":
        method_name = get_config().get("air_vac_conversion_method", "Edlen1966")
        if method_name == "Edlen1966":
            compute_fn = wavenumber_from_air_wavelength
        elif method_name == "Edlen1953":
            compute_fn = lambda wl: 1.0 / (astools.LvacE53(wl) * 1e-8)
        elif method_name == "Peck&Reeder1973":
            compute_fn = lambda wl: 1.0 / (astools.Lvac(wl) * 1e-8)
        elif method_name == "Meggers&Peters1919":
            compute_fn = lambda wl: 1.0 / (astools.LvacMP(wl) * 1e-8)
        else:
            compute_fn = wavenumber_from_air_wavelength
    else:
        compute_fn = wavenumber_from_vac_wavelength

    # Try fixing wavelength (single digit)
    fixed_wl = _try_single_digit_typo(
        wl_str, compute_fn, wavenumber, tolerance,
    )
    if fixed_wl:
        return {
            "status": "CORRECTED_WAVELENGTH",
            "computed_wavenumber": round(computed_sigma, 2),
            "difference": round(diff, 4),
            "corrected_wavelength": float(fixed_wl),
            "corrected_wavenumber": None,
            "note": f"Wavelength typo: {wl_str} → {fixed_wl}",
        }

    # Try fixing wavelength (swap adjacent digits)
    fixed_wl = _try_swap_adjacent_digits(
        wl_str, compute_fn, wavenumber, tolerance,
    )
    if fixed_wl:
        return {
            "status": "CORRECTED_WAVELENGTH",
            "computed_wavenumber": round(computed_sigma, 2),
            "difference": round(diff, 4),
            "corrected_wavelength": float(fixed_wl),
            "corrected_wavenumber": None,
            "note": f"Wavelength digit swap: {wl_str} → {fixed_wl}",
        }

    # Try fixing wavenumber (single digit)
    fixed_wn = _try_single_digit_typo(wn_str, lambda wn: wn, computed_sigma, tolerance)
    if fixed_wn:
        return {
            "status": "CORRECTED_WAVENUMBER",
            "computed_wavenumber": round(computed_sigma, 2),
            "difference": round(diff, 4),
            "corrected_wavelength": None,
            "corrected_wavenumber": float(fixed_wn),
            "note": f"Wavenumber typo: {wn_str} → {fixed_wn}",
        }

    # Try fixing wavenumber (swap adjacent digits)
    fixed_wn = _try_swap_adjacent_digits(wn_str, lambda wn: wn, computed_sigma, tolerance)
    if fixed_wn:
        return {
            "status": "CORRECTED_WAVENUMBER",
            "computed_wavenumber": round(computed_sigma, 2),
            "difference": round(diff, 4),
            "corrected_wavelength": None,
            "corrected_wavenumber": float(fixed_wn),
            "note": f"Wavenumber digit swap: {wn_str} → {fixed_wn}",
        }

    # Cannot fix
    return {
        "status": "ERROR",
        "computed_wavenumber": round(computed_sigma, 2),
        "difference": round(diff, 4),
        "corrected_wavelength": None,
        "corrected_wavenumber": None,
        "note": f"Wavelength/wavenumber mismatch: Δσ={diff:.4f} cm⁻¹ (expected ≤{tolerance})",
    }


# ──────────────────────────────────────────────
# Selection Rules
# ──────────────────────────────────────────────

def _parse_j_numeric(j_str: str) -> float | None:
    """Parse J value string to numeric. '9/2' → 4.5, '4' → 4.0."""
    if j_str is None:
        return None
    j_str = str(j_str).strip()
    if "/" in j_str:
        parts = j_str.split("/")
        try:
            return float(parts[0]) / float(parts[1])
        except (ValueError, ZeroDivisionError):
            return None
    try:
        return float(j_str)
    except ValueError:
        return None


def check_selection_rules(
    lower_parity: str | None,
    lower_j: str | None,
    upper_parity: str | None,
    upper_j: str | None,
) -> dict:
    """
    Check electric dipole (E1) selection rules for a transition.

    Rules:
    1. Parity must CHANGE: lower and upper must have different parity
    2. ΔJ = |J_upper - J_lower| ≤ 1
    3. J=0 → J=0 is forbidden (not relevant for half-integer J)

    Returns:
        dict with:
            - parity_ok: bool or None (if parity unknown)
            - delta_j_ok: bool or None (if J unknown)
            - delta_j: float or None
            - notes: list of violation descriptions
    """
    notes = []

    # Parity check
    parity_ok = None
    if lower_parity and upper_parity:
        lp = lower_parity.strip().lower()
        up = upper_parity.strip().lower()
        if lp in ("e", "o") and up in ("e", "o"):
            parity_ok = (lp != up)
            if not parity_ok:
                notes.append(f"Parity violation: both levels have parity '{lp}' (must change for E1)")

    # ΔJ check
    delta_j_ok = None
    delta_j = None
    j_lower = _parse_j_numeric(lower_j)
    j_upper = _parse_j_numeric(upper_j)
    if j_lower is not None and j_upper is not None:
        delta_j = abs(j_upper - j_lower)
        delta_j_ok = delta_j <= 1.0
        if not delta_j_ok:
            notes.append(f"ΔJ violation: |{upper_j} - {lower_j}| = {delta_j} (must be ≤ 1)")
        # J=0→J=0 forbidden (unlikely for half-integer, but check anyway)
        if j_lower == 0 and j_upper == 0:
            delta_j_ok = False
            notes.append("J=0 → J=0 transition is forbidden")

    return {
        "parity_ok": parity_ok,
        "delta_j_ok": delta_j_ok,
        "delta_j": delta_j,
        "notes": notes,
    }


# ──────────────────────────────────────────────
# Reference Level Lookup
# ──────────────────────────────────────────────

@lru_cache(maxsize=1)
def _load_reference_levels(ref_path: str | None = None) -> dict:
    """
    Load reference energy levels from the Excel file into lookup structures.

    Returns dict with:
        - exact: {(int_part, parity, J) → level_float}
        - by_int_parity: {(int_part, parity) → [(J, level_float), ...]}
        - by_int: {int_part → [(parity, J, level_float), ...]}
        - all_levels: [(int_part, parity, J, level_float), ...]
    """
    path = Path(ref_path) if ref_path else get_config()["ref_levels_path"]
    wb = openpyxl.load_workbook(str(path), read_only=True)
    ws = wb.active

    exact = {}
    by_int_parity = {}
    by_int = {}
    all_levels = []

    for row in ws.iter_rows(min_row=2, values_only=True):
        j_value = row[2]
        parity = row[4]
        level_str = row[5]

        if j_value is None or level_str is None or parity is None:
            continue

        j_str = str(j_value).strip()
        if j_str == "---":
            continue

        try:
            level_float = float(str(level_str).strip())
        except ValueError:
            continue
            
        # Auto-convert eV to cm-1 if the value looks like it's in eV (e.g. max level ~15-20 eV vs 100,000 cm-1. Ground state is 0.0)
        # We can loosely guess it's eV if the level is >0 but very small, however it's safer to just provide 
        # the conversion capability. To be robust, we'll check the cell header if possible, or just convert
        # if the max level in the document seems to be < 1000.
        # For now, we will expose the conversion based on publication year as required by the plan.
        if level_float > 0 and level_float < 1000.0:  
            # Assuming values < 1000 cm-1 are actually eV for this ion (Pr III levels go up to 100,000s).
            # This handles older papers that might have published in eV.
            pub_year = get_config().get("publication_year", 2022)
            level_float = level_float * astools.evcm(pub_year)

        int_part = int(math.floor(level_float))
        parity_str = str(parity).strip().lower()
        if parity_str not in ("e", "o"):
            continue

        key = (int_part, parity_str, j_str)
        exact[key] = level_float

        key2 = (int_part, parity_str)
        if key2 not in by_int_parity:
            by_int_parity[key2] = []
        by_int_parity[key2].append((j_str, level_float))

        if int_part not in by_int:
            by_int[int_part] = []
        by_int[int_part].append((parity_str, j_str, level_float))

        all_levels.append((int_part, parity_str, j_str, level_float))

    wb.close()
    return {
        "exact": exact,
        "by_int_parity": by_int_parity,
        "by_int": by_int,
        "all_levels": all_levels,
    }


def lookup_reference_level(level_int: int, parity: str, j_value: str) -> dict:
    """
    Look up an energy level in Reference Levels, with auto-correction attempts.

    Tries in order:
    1. Exact match (int_part, parity, J)
    2. ±1 in integer part with exact parity+J
    3. Flipped parity (° misread as 0 or vice versa) with same J
    4. Same int_part+parity, different J from Reference
    5. Int_part typo correction (single digit, adjacent swap)

    Returns:
        dict with found, exact_level, corrected_parity, corrected_j, corrected_int, note
    """
    data = _load_reference_levels()
    j_str = str(j_value).strip()
    par = parity.strip().lower() if parity else "e"

    result_base = {
        "found": False,
        "exact_level": None,
        "corrected_parity": None,
        "corrected_j": None,
        "corrected_int": None,
        "note": "",
    }

    # 1. Exact match
    key = (level_int, par, j_str)
    if key in data["exact"]:
        return {**result_base, "found": True, "exact_level": data["exact"][key]}

    # 2. ±1 in integer part with exact parity+J
    for delta in [-1, 1]:
        alt_key = (level_int + delta, par, j_str)
        if alt_key in data["exact"]:
            return {
                **result_base,
                "found": True,
                "exact_level": data["exact"][alt_key],
                "corrected_int": level_int + delta,
                "note": f"Int part ±1: {level_int}→{level_int + delta}",
            }

    # 3. Flipped parity with same J (° misread as 0 or vice versa)
    flipped_par = "o" if par == "e" else "e"
    key_flip = (level_int, flipped_par, j_str)
    if key_flip in data["exact"]:
        return {
            **result_base,
            "found": True,
            "exact_level": data["exact"][key_flip],
            "corrected_parity": flipped_par,
            "note": f"Parity corrected: {par}→{flipped_par} (likely ° misread)",
        }
    # Also try flipped parity with ±1 int
    for delta in [-1, 1]:
        key_flip_delta = (level_int + delta, flipped_par, j_str)
        if key_flip_delta in data["exact"]:
            return {
                **result_base,
                "found": True,
                "exact_level": data["exact"][key_flip_delta],
                "corrected_parity": flipped_par,
                "corrected_int": level_int + delta,
                "note": f"Parity corrected: {par}→{flipped_par}, int ±1: {level_int}→{level_int + delta}",
            }

    # 4. Same int_part, any parity, different J from Reference
    #    Check with original parity first, then flipped
    for check_par in [par, flipped_par]:
        key2 = (level_int, check_par)
        if key2 in data["by_int_parity"]:
            candidates = data["by_int_parity"][key2]
            if len(candidates) == 1:
                # Only one J available → likely Reference revised the J value
                ref_j, ref_level = candidates[0]
                corrections = []
                if check_par != par:
                    corrections.append(f"parity {par}→{check_par}")
                corrections.append(f"J {j_str}→{ref_j} (Reference)")
                return {
                    **result_base,
                    "found": True,
                    "exact_level": ref_level,
                    "corrected_parity": check_par if check_par != par else None,
                    "corrected_j": ref_j,
                    "note": "Corrected: " + ", ".join(corrections),
                }
        # Also try ±1 int with different J
        for delta in [-1, 1]:
            key2_d = (level_int + delta, check_par)
            if key2_d in data["by_int_parity"]:
                candidates = data["by_int_parity"][key2_d]
                if len(candidates) == 1:
                    ref_j, ref_level = candidates[0]
                    corrections = []
                    if check_par != par:
                        corrections.append(f"parity {par}→{check_par}")
                    corrections.append(f"int {level_int}→{level_int + delta}")
                    corrections.append(f"J {j_str}→{ref_j} (Reference)")
                    return {
                        **result_base,
                        "found": True,
                        "exact_level": ref_level,
                        "corrected_parity": check_par if check_par != par else None,
                        "corrected_j": ref_j,
                        "corrected_int": level_int + delta,
                        "note": "Corrected: " + ", ".join(corrections),
                    }

    # 5. Integer part typo correction (single digit or adjacent swap)
    int_str = str(level_int)
    for check_par_inner in [par, flipped_par]:
        # Single digit typo
        for i, ch in enumerate(int_str):
            if ch.isdigit():
                for d in "0123456789":
                    if d == ch:
                        continue
                    candidate_str = int_str[:i] + d + int_str[i + 1:]
                    try:
                        candidate_int = int(candidate_str)
                    except ValueError:
                        continue
                    # Check if this integer part exists with the right parity and J
                    cand_key = (candidate_int, check_par_inner, j_str)
                    if cand_key in data["exact"]:
                        corrections = [f"int typo: {level_int}→{candidate_int}"]
                        if check_par_inner != par:
                            corrections.append(f"parity {par}→{check_par_inner}")
                        return {
                            **result_base,
                            "found": True,
                            "exact_level": data["exact"][cand_key],
                            "corrected_parity": check_par_inner if check_par_inner != par else None,
                            "corrected_int": candidate_int,
                            "note": "Corrected: " + ", ".join(corrections),
                        }

        # Adjacent digit swap
        chars = list(int_str)
        for i in range(len(chars) - 1):
            if chars[i] != chars[i + 1]:
                swapped = chars.copy()
                swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
                candidate_str = "".join(swapped)
                try:
                    candidate_int = int(candidate_str)
                except ValueError:
                    continue
                cand_key = (candidate_int, check_par_inner, j_str)
                if cand_key in data["exact"]:
                    corrections = [f"int digit swap: {level_int}→{candidate_int}"]
                    if check_par_inner != par:
                        corrections.append(f"parity {par}→{check_par_inner}")
                    return {
                        **result_base,
                        "found": True,
                        "exact_level": data["exact"][cand_key],
                        "corrected_parity": check_par_inner if check_par_inner != par else None,
                        "corrected_int": candidate_int,
                        "note": "Corrected: " + ", ".join(corrections),
                    }

    # Nothing worked
    # Provide helpful info about what IS available
    info_parts = []
    if level_int in data["by_int"]:
        available = data["by_int"][level_int]
        info_parts.append(f"Level {level_int} exists as: " +
                          ", ".join(f"{p} J={j}" for p, j, _ in available))
    return {
        **result_base,
        "note": f"Level {level_int} parity={par} J={j_str} not found in Reference Levels" +
                (f". {'; '.join(info_parts)}" if info_parts else ""),
    }


def validate_reference_levels(
    lower_level_int: int | None,
    lower_parity: str | None,
    lower_j: str | None,
    upper_level_int: int | None,
    upper_parity: str | None,
    upper_j: str | None,
    observed_wavenumber: float,
    tolerance: float = None,
) -> dict:
    """
    Validate a spectral line's classification against Reference levels.

    Steps:
    1. Look up both levels (with auto-correction)
    2. Check selection rules
    3. If selection rules fail with raw values but pass with Reference-corrected values, use corrected
    4. Check Ritz wavenumber agreement

    Returns dict with:
        lower_exact, upper_exact, corrected_lower_parity, corrected_lower_j,
        corrected_upper_parity, corrected_upper_j, ritz_wavenumber, ritz_difference,
        status, note
    """
    if tolerance is None:
        tolerance = get_config()["ref_ritz_tolerance_cm1"]

    # Unclassified line
    if lower_level_int is None or upper_level_int is None:
        return {
            "lower_exact": None, "upper_exact": None,
            "corrected_lower_parity": None, "corrected_lower_j": None,
            "corrected_upper_parity": None, "corrected_upper_j": None,
            "corrected_lower_int": None, "corrected_upper_int": None,
            "ritz_wavenumber": None, "ritz_difference": None,
            "status": "UNCLASSIFIED", "note": "Unclassified line",
        }

    # Look up both levels (with auto-correction)
    lower = lookup_reference_level(lower_level_int, lower_parity or "e", lower_j or "0")
    upper = lookup_reference_level(upper_level_int, upper_parity or "e", upper_j or "0")

    notes = []
    if lower["note"]:
        notes.append(f"Lower: {lower['note']}")
    if upper["note"]:
        notes.append(f"Upper: {upper['note']}")

    result = {
        "lower_exact": lower["exact_level"],
        "upper_exact": upper["exact_level"],
        "corrected_lower_parity": lower.get("corrected_parity"),
        "corrected_lower_j": lower.get("corrected_j"),
        "corrected_lower_int": lower.get("corrected_int"),
        "corrected_upper_parity": upper.get("corrected_parity"),
        "corrected_upper_j": upper.get("corrected_j"),
        "corrected_upper_int": upper.get("corrected_int"),
        "ritz_wavenumber": None,
        "ritz_difference": None,
        "status": "LEVEL_NOT_FOUND",
        "note": "; ".join(notes) if notes else "",
    }

    # If either level not found, return early — do NOT check selection rules
    # (the uncorrected values would be misleading)
    if not lower["found"] or not upper["found"]:
        return result

    # Both levels found — NOW check selection rules with effective (corrected) values
    eff_lower_parity = lower.get("corrected_parity") or lower_parity
    eff_lower_j = lower.get("corrected_j") or lower_j
    eff_upper_parity = upper.get("corrected_parity") or upper_parity
    eff_upper_j = upper.get("corrected_j") or upper_j

    sel = check_selection_rules(eff_lower_parity, eff_lower_j, eff_upper_parity, eff_upper_j)
    if sel["notes"]:
        notes.extend([f"Selection rule: {n}" for n in sel["notes"]])

    # Compute Ritz wavenumber
    ritz = upper["exact_level"] - lower["exact_level"]
    ritz_diff = abs(ritz - observed_wavenumber)

    result["ritz_wavenumber"] = round(ritz, 4)
    result["ritz_difference"] = round(ritz_diff, 4)

    if ritz_diff <= tolerance:
        # Check if selection rules are also satisfied
        if sel["parity_ok"] is False or sel["delta_j_ok"] is False:
            result["status"] = "RITZ_OK_SELECTION_FAIL"
            result["note"] = "; ".join(notes) + f"; Ritz OK (Δ={ritz_diff:.4f}) but selection rules violated"
        else:
            result["status"] = "OK"
    else:
        result["status"] = "RITZ_MISMATCH"
        mismatch_note = f"Ritz mismatch: |{ritz:.4f} - {observed_wavenumber:.2f}| = {ritz_diff:.4f} cm⁻¹ (> {tolerance})"
        notes.append(mismatch_note)
        result["note"] = "; ".join(notes)

    return result


def validate_full_row(row: dict) -> dict:
    """
    Full validation of a single extracted row.
    
    Strategy:
    1. First do Reference lookup (with auto-correction of parity, J, energy level int)
    2. If Reference lookup OK → do wavelength/wavenumber validation → correct typos if needed
    3. If Reference lookup FAILS → do wavelength/wavenumber validation but DO NOT correct
       (because the mismatch is likely due to misaligned classification, not a wavelength typo)
    4. Check selection rules on the corrected values

    Returns the row dict with added validation fields.
    """
    wl = row.get("wavelength")
    wn = row.get("wavenumber")
    lower_int = row.get("lower_level_int")
    lower_par = row.get("lower_parity")
    lower_j = row.get("lower_j")
    upper_int = row.get("upper_level_int")
    upper_par = row.get("upper_parity")
    upper_j = row.get("upper_j")

    # Parse numeric values
    try:
        wl = float(wl) if wl is not None else None
    except (ValueError, TypeError):
        wl = None
    try:
        wn = float(wn) if wn is not None else None
    except (ValueError, TypeError):
        wn = None
        
    obs_ritz = row.get("obs_ritz", None)
    if obs_ritz is not None:
        try:
            obs_ritz = float(obs_ritz)
        except (ValueError, TypeError):
            obs_ritz = None

    is_classified = lower_int is not None and upper_int is not None

    # ── Step 1: Reference validation ──
    ref_result = validate_reference_levels(
        lower_int, lower_par, lower_j,
        upper_int, upper_par, upper_j,
        wn if wn else 0.0,
    )

    row["lower_exact_ref"] = ref_result["lower_exact"]
    row["upper_exact_ref"] = ref_result["upper_exact"]
    row["ref_status"] = ref_result["status"]
    row["ref_note"] = ref_result["note"]

    # Apply corrections from Reference lookup
    if ref_result.get("corrected_lower_parity"):
        row["lower_parity_original"] = lower_par
        row["lower_parity"] = ref_result["corrected_lower_parity"]
    if ref_result.get("corrected_lower_j"):
        row["lower_j_original"] = lower_j
        row["lower_j"] = ref_result["corrected_lower_j"]
    if ref_result.get("corrected_lower_int"):
        row["lower_level_int_original"] = lower_int
        row["lower_level_int"] = ref_result["corrected_lower_int"]
    if ref_result.get("corrected_upper_parity"):
        row["upper_parity_original"] = upper_par
        row["upper_parity"] = ref_result["corrected_upper_parity"]
    if ref_result.get("corrected_upper_j"):
        row["upper_j_original"] = upper_j
        row["upper_j"] = ref_result["corrected_upper_j"]
    if ref_result.get("corrected_upper_int"):
        row["upper_level_int_original"] = upper_int
        row["upper_level_int"] = ref_result["corrected_upper_int"]

    # ── Step 2: Wavelength ↔ Wavenumber validation ──
    if wl is not None and wn is not None:
        ref_ok = ref_result["status"] in ("OK", "RITZ_OK_SELECTION_FAIL", "UNCLASSIFIED")

        wl_wn_result = validate_wavelength_wavenumber(wl, wn)

        if wl_wn_result["status"] == "OK":
            row["wl_wn_status"] = "OK"
            row["wl_wn_note"] = ""
        elif ref_ok or not is_classified:
            # Reference is OK (or unclassified) → wavelength/wavenumber typo correction IS safe
            row["wl_wn_status"] = wl_wn_result["status"]
            row["wl_wn_note"] = wl_wn_result["note"]
            if wl_wn_result["status"] == "CORRECTED_WAVELENGTH" and wl_wn_result["corrected_wavelength"]:
                row["wavelength_original"] = wl
                row["wavelength"] = wl_wn_result["corrected_wavelength"]
            elif wl_wn_result["status"] == "CORRECTED_WAVENUMBER" and wl_wn_result["corrected_wavenumber"]:
                row["wavenumber_original"] = wn
                row["wavenumber"] = wl_wn_result["corrected_wavenumber"]
        else:
            # Reference failed → likely misaligned classification → do NOT correct wl/wn
            # Just report the mismatch without applying corrections
            row["wl_wn_status"] = "MISMATCH_SUSPECT_MISALIGN"
            row["wl_wn_note"] = (
                f"Wl/Wn mismatch (Δσ={wl_wn_result['difference']:.4f} cm⁻¹) "
                f"NOT corrected — Reference validation also failed, likely misaligned classification"
            )
    else:
        row["wl_wn_status"] = "MISSING_DATA"
        row["wl_wn_note"] = "Missing wavelength or wavenumber"

    # ── Step 3: obs-Ritz validation ──
    if obs_ritz is not None and ref_result.get("ritz_wavenumber") is not None and wn is not None:
        expected_obs_ritz = wn - ref_result["ritz_wavenumber"]
        row["expected_obs_ritz"] = round(expected_obs_ritz, 4)
        diff = abs(obs_ritz - expected_obs_ritz)
        if diff <= 0.02: # Tolerance for printed obs-ritz
            row["obs_ritz_status"] = "OK"
        else:
            row["obs_ritz_status"] = "MISMATCH"
            row["obs_ritz_note"] = f"Printed obs-Ritz {obs_ritz} != Expected {expected_obs_ritz:.4f}"

    row["computed_wavenumber"] = wl_wn_result["computed_wavenumber"] if wl is not None and wn is not None else None

    return row
