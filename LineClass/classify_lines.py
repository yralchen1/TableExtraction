"""
classify_lines.py
=================
Classify observed spectral lines of Pr III by matching them against all
possible transitions between known energy levels. Grade each match by
wavenumber agreement and intensity consistency, resolve conflicts, and
output a sorted classification table to Excel.

Steps:
  1. Read energy levels from Pr3_lev_Wyart_1999.xlsm
  2. Read calculated transitions from Icalc.xlsx
  3. Read observed spectral lines from Pr3_lines.xlsx
  4. Generate all possible transitions satisfying selection rules
  5. Match observed lines to possible transitions & grade
  6. Resolve conflicts & write output
"""

import os
import math
from typing import List
import numpy as np
import bisect
import pandas as pd
import openpyxl
from models import EnergyLevel, SpectralLine, Transition, UNASSIGNED

# ---------------------------------------------------------------------------
# Paths (relative to this script's directory)
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LEVELS_FILE = os.path.join(SCRIPT_DIR, '..', 'TableExtraction', 'Pr3_lev_Wyart_1999.xlsm')
LINES_FILE = os.path.join(SCRIPT_DIR, 'Pr3_lines.xlsx')
ICALC_FILE = os.path.join(SCRIPT_DIR, 'Icalc.xlsx')
OUTPUT_FILE = os.path.join(SCRIPT_DIR, 'line_classifications.xlsx')
OUTPUT_CSV = os.path.join(SCRIPT_DIR, 'line_classifications.csv')

# Wavenumber range for possible transitions (cm^{-1})
WN_MIN = 9327.0
WN_MAX = 121665.0

# Atomic mass and plasma temperature for Doppler width calculation
ATOMIC_MASS = 140.90765  # u, standard mass unit
T = 1.6   # eV, plasma temperature


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def parse_J(j_str: str) -> float:
    """Parse a J quantum number string to a float.
    '5/2' -> 2.5, '3' -> 3.0
    """
    j_str = str(j_str).strip()
    if '/' in j_str:
        num, den = j_str.split('/')
        return float(num) / float(den)
    return float(j_str)


# noinspection PyUnresolvedReferences
def same_transitions(t1: Transition, t2: Transition) -> bool:
    """Return True if two transitions connect the same pair of levels."""
    return (t1.lower_level.level_id == t2.lower_level.level_id and
            t1.upper_level.level_id == t2.upper_level.level_id)




def to_str_id(val) -> str:
    """Convert a level_id value (possibly read as float from Excel) to string.
    E.g., 123.0 -> '123', 'abc' -> 'abc', None -> ''
    """
    if val is None:
        return ''
    if isinstance(val, float):
        if val != val:  # NaN check
            return ''
        return str(int(val))
    return str(val).strip()


# ===========================================================================
# STEP 1: Read energy levels
# ===========================================================================
def read_energy_levels() -> tuple[dict[str, EnergyLevel], list[EnergyLevel]]:
    """Read energy levels from the Wyart2000 worksheet.

    Returns:
        levels_dict: dict keyed by level_id (str) -> EnergyLevel
        levels_list: ordered list of EnergyLevel objects
    """
    print("Step 1: Reading energy levels...")
    wb = openpyxl.load_workbook(LEVELS_FILE, read_only=True, data_only=True)
    ws = wb['Wyart2000']

    levels_dict = {}
    levels_list = []

    for row in ws.iter_rows(min_row=2):  # skip header
        # Column indices are 0-based in openpyxl rows, but spec uses 1-based
        # J_adopt = col 11 (index 10), E_adopt = col 13 (index 12),
        # par_ASD = col 14 (index 13), ASD_id = col 23 (index 22)
        j_str_val = row[10].value   # J_adopt (col 11)
        energy_val = row[12].value  # E_adopt (col 13)
        parity_val = row[13].value  # par_ASD (col 14)
        level_id_val = row[22].value  # ASD_id (col 23)

        # Skip rows with missing essential data
        if energy_val is None or j_str_val is None or level_id_val is None:
            continue

        level_id = to_str_id(level_id_val)
        if level_id == '':
            continue

        j_str = str(j_str_val).strip()
        parity = str(parity_val).strip() if parity_val is not None else ''

        try:
            # noinspection PyTypeChecker
            energy = float(energy_val)
            j_val = parse_J(j_str)
        except (ValueError, TypeError):
            continue

        lev = EnergyLevel(
            level_id=level_id,
            energy=energy,
            parity=parity,
            J_str=j_str,
            J_val=j_val
        )
        levels_dict[level_id] = lev
        levels_list.append(lev)

    wb.close()
    print(f"  Read {len(levels_list)} energy levels.")
    return levels_dict, levels_list


# ===========================================================================
# STEP 2: Read calculated transitions
# ===========================================================================
# noinspection PyTypeChecker
def read_transitions(levels_dict: dict) -> dict:
    """Read calculated transitions from Icalc.xlsx.

    Returns:
        calc_trans_index: dict keyed by (lower_id, upper_id) -> dict with
                          keys 'calc_intensity', 'u_calc', 'assigned_to'
    """
    print("Step 2: Reading calculated transitions...")
    wb = openpyxl.load_workbook(ICALC_FILE, read_only=True, data_only=True)
    ws = wb['Sheet1']

    calc_trans_index = {}

    for row in ws.iter_rows(min_row=2):  # skip header
        # Col 1: id1_W99 (idx 0), Col 2: id2_W99 (idx 1),
        # Col 5: u%gA (idx 4), Col 6: Icalc (idx 5)
        id1_val = to_str_id(row[0].value)
        id2_val = to_str_id(row[1].value)

        if not id1_val or not id2_val:
            continue

        if id1_val not in levels_dict or id2_val not in levels_dict:
            continue

        ua_pcnt_val = row[4].value
        icalc_val = row[5].value

        u_calc = None
        if ua_pcnt_val is not None:
            try:
                u_calc = math.log(float(ua_pcnt_val) / 100.0 + 1.0)
            except (ValueError, TypeError):
                pass

        calc_intensity = None
        if icalc_val is not None:
            try:
                calc_intensity = float(icalc_val)
            except (ValueError, TypeError):
                pass

        calc_trans_index[(id1_val, id2_val)] = {
            'calc_intensity': calc_intensity,
            'u_calc': u_calc,
            'assigned_to': None
        }

    wb.close()
    print(f"  Read {len(calc_trans_index)} calculated transitions.")
    return calc_trans_index


# ===========================================================================
# STEP 3: Read observed spectral lines
# ===========================================================================
# noinspection PyTypeChecker,PyUnresolvedReferences
def read_observed_lines(levels_dict: dict, calc_trans_index: dict):
    """Read observed spectral lines from Pr3_lines.xlsx.

    Returns:
        observed_lines: list of SpectralLine objects
    """
    print("Step 3: Reading observed spectral lines...")
    wb = openpyxl.load_workbook(LINES_FILE, read_only=True, data_only=True)
    ws = wb['Sheet1']

    observed_lines = []
    prev_line = None

    for row in ws.iter_rows(min_row=2):  # skip header
        wn_val = row[0].value       # wavenumber
        unc_val = row[1].value      # wn_uncertainty
        intens_val = row[2].value   # intensity
        char_val = row[3].value     # line_character

        if wn_val is None or unc_val is None or intens_val is None:
            continue

        try:
            wavenumber = float(wn_val)
            uncertainty = float(unc_val)
            intensity = float(intens_val)
        except (ValueError, TypeError):
            continue

        line_char = str(char_val).strip() if char_val is not None else ''

        id1_val = to_str_id(row[4].value)  # lower level_id
        id2_val = to_str_id(row[5].value)  # upper level_id

        if prev_line is not None and wavenumber == prev_line.wavenumber:
            current_line = prev_line
        else:
            current_line = SpectralLine(
                wavenumber=wavenumber,
                wn_uncertainty=uncertainty,
                intensity=intensity,
                line_character=line_char
            )
            observed_lines.append(current_line)
            prev_line = current_line

        if id1_val and id2_val:
            lower = levels_dict.get(id1_val)
            upper = levels_dict.get(id2_val)
            if lower is not None and upper is not None:
                calc_data = calc_trans_index.get((id1_val, id2_val))
                tr = Transition(
                    lower_level=lower,
                    upper_level=upper,
                    calc_intensity=calc_data['calc_intensity'] if calc_data else None,
                    u_calc=calc_data['u_calc'] if calc_data else None,
                    assigned_to=current_line,
                    orig_calc_intensity=calc_data['calc_intensity'] if calc_data else None,
                    orig_u_calc=calc_data['u_calc'] if calc_data else None,
                )
                current_line.assigned_transitions.append(tr)
                current_line.original_assignments.append(tr)
                if calc_data:
                    calc_data['assigned_to'] = current_line

    wb.close()
    print(f"  Read {len(observed_lines)} observed lines.")
    return observed_lines


# ===========================================================================
# STEP 4: Generate all possible transitions
# ===========================================================================
# noinspection PyUnresolvedReferences
def generate_all_possible_transitions(levels_list: list,
                                       calc_trans_index: dict) -> list:
    """Generate all possible transitions between known levels satisfying
    selection rules, within the wavenumber range [WN_MIN, WN_MAX].

    Selection rules (J-based only, per spec):
      - |J_up - J_lo| <= 1
      - Not both J = 0

    Returns:
        all_possible_transitions: list of Transition objects sorted by
                                   calculated_wavenumber (ascending)
    """
    print("Step 4: Generating all possible transitions...")
    all_possible = []
    n_levels = len(levels_list)
    levels_list = sorted(levels_list, key=lambda l: l.energy)  # sort by energy for efficiency

    for i in range(n_levels):
        lev_lo = levels_list[i]
        for j in range(i+1, n_levels):
            lev_up = levels_list[j]
            if lev_up.parity == lev_lo.parity:
                continue

            wn = lev_up.energy - lev_lo.energy

            # Wavenumber range check
            if wn < WN_MIN or wn > WN_MAX:
                continue

            # Selection rules
            delta_j = abs(lev_up.J_val - lev_lo.J_val)
            if delta_j > 1.5:  # Allowing for some numerical imprecision in J values
                continue
            if lev_up.J_val == 0.0 and lev_lo.J_val == 0.0:
                continue

            # Look up in calculated transitions
            key = (lev_lo.level_id, lev_up.level_id)
            calc_data = calc_trans_index.get(key)

            if calc_data is not None:
                tr = Transition(
                    lower_level=lev_lo,
                    upper_level=lev_up,
                    calc_intensity=calc_data['calc_intensity'],
                    u_calc=calc_data['u_calc'],
                    assigned_to=calc_data['assigned_to'],
                    orig_calc_intensity=calc_data['calc_intensity'],
                    orig_u_calc=calc_data['u_calc'],
                )
            else:
                tr = Transition(
                    lower_level=lev_lo,
                    upper_level=lev_up
                )

            all_possible.append(tr)

    # Sort by calculated wavenumber
    all_possible.sort(key=lambda t: t.calculated_wavenumber)
    print(f"  Generated {len(all_possible)} possible transitions.")
    return all_possible


# ===========================================================================
# STEP 5: Match observed lines & grade assignments
# ===========================================================================
def assign_grades(line: SpectralLine, transitions: list):
    """Implement new 2D grading scheme based on WN residual and intensity consistency."""
    ln2 = math.log(2)
    ln5 = math.log(5)

    for tr in transitions:
        wn_diff_abs = abs(line.wavenumber - tr.calculated_wavenumber)
        # Tier (Wavenumber Agreement)
        if wn_diff_abs <= 2.0 * line.wn_uncertainty:
            tier = '2'
        elif wn_diff_abs <= 3.0 * line.wn_uncertainty:
            tier = '3'
        elif wn_diff_abs <= 4.0 * line.wn_uncertainty:
            tier = '4'
        else:
            tier = '5'

        # Subgrade (Intensity Consistency)
        if tr.calc_intensity is None or tr.u_calc is None:
            subgrade = 'G'
        else:
            ln_u_I_calc = tr.u_calc
            ln_I_ratio = abs(math.log(tr.calc_intensity / line.intensity)) if line.intensity > 0 else float('inf')
            diff = ln_I_ratio - ln_u_I_calc

            if diff <= 0 and ln_u_I_calc <= ln2:
                subgrade = 'A'
            elif diff <= ln2 and ln_u_I_calc <= ln5:
                subgrade = 'B'
            elif ln2 < diff <= ln5 and ln2 < ln_u_I_calc <= ln5:
                subgrade = 'C'
            elif diff <= ln5 < ln_u_I_calc:
                subgrade = 'D'
            elif diff > ln5 and ln_u_I_calc > ln5:
                subgrade = 'E'
            else:
                subgrade = 'G' # Fallback

        tr.grade = tier + subgrade

        # Handle UNASSIGNED removal
        if UNASSIGNED in line.assigned_transitions:
            line.assigned_transitions.remove(UNASSIGNED)

        # Set new flag: 1 = new classification, 0 = original
        is_original = any(same_transitions(tr, orig) for orig in line.original_assignments)
        tr.new = 0 if is_original else 1
        tr.notes1 = ""  # notes1 holds F/R flags from conflict resolution


def match_and_grade(observed_lines: list, all_possible_transitions: list, verbose: bool = False):
    """For each observed line, find matching transitions and grade them.

    Modifies observed_lines in place (updates assigned_transitions and grades).
    Also tracks all new assignments for conflict resolution.

    Returns:
        transition_assignments: dict keyed by (lower_id, upper_id) ->
                                list of (Transition, SpectralLine, wn_diff_abs)
                                for conflict resolution
    """
    if verbose: print("Step 5.1: Matching observed lines to transitions & grading...")

    # Build a list of wavenumbers for binary search
    all_wn = [t.calculated_wavenumber for t in all_possible_transitions]

    # Track all assignments: (lower_id, upper_id) -> [(transition, line, wn_diff)]
    transition_assignments = {}

    for line_idx, obs_line in enumerate(observed_lines):
        # if (line_idx + 1) % 1000 == 0:
        #     print(f"  Processing line {line_idx + 1}/{len(observed_lines)}...")

        tolerance = 5.5 * obs_line.wn_uncertainty
        wn_lo = obs_line.wavenumber - tolerance
        wn_hi = obs_line.wavenumber + tolerance

        idx_lo = bisect.bisect_left(all_wn, wn_lo)
        idx_hi = bisect.bisect_right(all_wn, wn_hi)

        matches = []
        for t in all_possible_transitions[idx_lo:idx_hi]:
            wn_diff_abs = abs(obs_line.wavenumber - t.calculated_wavenumber)
            if wn_diff_abs <= tolerance:
                # Reuse or create transition
                existing = None
                for otr in obs_line.assigned_transitions:
                    if same_transitions(t, otr):
                        existing = otr
                        break
                
                if existing:
                    target_tr = existing
                else:
                    target_tr = Transition(
                        lower_level=t.lower_level,
                        upper_level=t.upper_level,
                        calc_intensity=t.calc_intensity,
                        u_calc=t.u_calc,
                        assigned_to=obs_line,
                        orig_calc_intensity=t.orig_calc_intensity,
                        orig_u_calc=t.orig_u_calc,
                    )
                    obs_line.assigned_transitions.append(target_tr)
                
                matches.append(target_tr)

        if matches:
            assign_grades(obs_line, matches)
            for m in matches:
                key = (m.lower_level.level_id, m.upper_level.level_id)
                if key not in transition_assignments:
                    transition_assignments[key] = []
                transition_assignments[key].append((m, obs_line, abs(obs_line.wavenumber - m.calculated_wavenumber)))

    if verbose: print(f"  Matching complete.")
    return transition_assignments


# ===========================================================================
# STEP 5.2: Conflict resolution, "R" tagging, and output
# ===========================================================================
def resolve_conflicts(transition_assignments: dict, verbose: bool = False):
    """Resolve conflicts and tag 'R' (Revised) notes1."""
    if verbose: print("Step 5.2: Resolving conflicts...")
    
    # 1. Selection for multiple lines pointing to one transition
    max_num_assignments = 0
    num_conflicts = 0
    for key, assignments in transition_assignments.items():
        if len(assignments) > max_num_assignments:
            max_num_assignments = len(assignments)
        if len(assignments) <= 1: continue

        num_conflicts += 1

        # --- Shared helpers for all scoring approaches ---
        # noinspection PyTypeChecker
        def _get_z_and_int_err(t, lin, wn_diff_abs):
            """Compute z (sigma-normalized residual) and int_err for a candidate."""
            sigma = getattr(lin, "wn_uncertainty", None)
            if not sigma or sigma <= 0:
                sigma = 1e-6
            z = wn_diff_abs / sigma

            Iobs = getattr(lin, "intensity", None)
            Icalc = getattr(t, "calc_intensity", None)
            u_Icalc = getattr(t, "u_calc", None)

            if Iobs and Iobs > 0 and Icalc and Icalc > 0 and u_Icalc is not None:
                # noinspection PyUnresolvedReferences
                int_err = max(abs(math.log(Icalc / Iobs)) - u_Icalc, 0.0)
            else:
                int_err = 0.0

            return z, int_err

        # =====================================================================
        # APPROACH A: Lexicographic Sort (INACTIVE)
        # Sorts by (tier, z, int_err, is_new). No weight tuning needed.
        # Tier dominates, then z breaks ties, then intensity, then new-vs-old.
        # =====================================================================
        # def sort_key_lexicographic(item):
        #     t, lin, wn_diff_abs = item
        #     z, int_err = _get_z_and_int_err(t, lin, wn_diff_abs)
        #     tier = 2 if z <= 2 else 3 if z <= 3 else 4 if z <= 4 else 5
        #     is_new = t.new if t.new is not None else 0
        #     return tier, z, int_err, is_new
        #
        # sort_key = sort_key_lexicographic

        # =====================================================================
        # APPROACH B: Normalized chi-squared Score (INACTIVE)
        # Combines z and intensity into a single chi² with one tunable scale.
        # Uncomment this block and comment out Approach A to activate.
        # =====================================================================
        # INT_ERR_SCALE = 1.5  # Normalizes int_err to z-equivalent units
        #
        # def sort_key_chi2(item):
        #     t, lin, wn_diff_abs = item
        #     z, int_err = _get_z_and_int_err(t, lin, wn_diff_abs)
        #     is_new = t.new if t.new is not None else 0
        #     chi2 = z * z + (int_err / INT_ERR_SCALE) ** 2 + 0.5 * is_new
        #     return chi2
        #
        # sort_key = sort_key_chi2

        # =====================================================================
        # APPROACH C: Rank-Based Scoring (INACTIVE)
        # Uses ranks within the conflict group, robust to outliers.
        # Uncomment this block and comment out Approach A to activate.
        # =====================================================================
        # def sort_key_rank(item):
        #     """Compute rank-based score. Must be called after computing all
        #     z and int_err values for the conflict group."""
        #     t, lin, wn_diff_abs = item
        #     z, int_err = _get_z_and_int_err(t, lin, wn_diff_abs)
        #     is_new = t.new if t.new is not None else 0
        #     # Note: ranks are computed below after building the values list
        #     return (z, int_err, is_new)  # placeholder; see rank sort below
        #
        # # To use rank-based scoring, replace the assignments.sort() call below:
        # # vals = [(z, ie, is_new) for (t, l, d) in assignments
        # #         for z, ie in [_get_z_and_int_err(t, l, d)]
        # #         for is_new in [t.new if t.new is not None else 0]]
        # # z_ranks = {v: r for r, v in enumerate(sorted(set(v[0] for v in vals)))}
        # # ie_ranks = {v: r for r, v in enumerate(sorted(set(v[1] for v in vals)))}
        # # ranked = [(z_ranks[v[0]] + 0.3 * ie_ranks[v[1]] + 0.2 * v[2], i)
        # #           for i, v in enumerate(vals)]
        # # ranked.sort()
        # # assignments = [assignments[i] for _, i in ranked]

        # =====================================================================
        # APPROACH D: Conservative Sort (ACTIVE)
        # Sorts by (is_new, |dwn_O-R|, tier, z, int_err ). No weight tuning needed.
        # Old assignments are favored, then the absolute "Obs-Ritz" mismatch,
        # then tier, then z breaks ties, then intensity mismatch.
        # =====================================================================
        def sort_key_conservative(item):
            t, lin, wn_diff_abs = item
            z, int_err = _get_z_and_int_err(t, lin, wn_diff_abs)
            tier = 2 if z <= 2 else 3 if z <= 3 else 4 if z <= 4 else 5
            is_new = t.new if t.new is not None else 0
            return is_new, wn_diff_abs, tier, z, int_err

        sort_key = sort_key_conservative

        assignments.sort(key=sort_key)
        winner_tr, winner_line, _ = assignments[0]
        winner_tr.notes1 += 'F'

        # Check if any original classification is being moved
        # (i.e., an original classification loses to a different line)
        original_moved = False
        for tr, line, _ in assignments[1:]:
            if tr.new == 0:
                original_moved = True

        if original_moved:
            winner_tr.notes1 += 'R'

        for tr, line, _ in assignments[1:]:
            if tr in line.assigned_transitions:
                line.assigned_transitions.remove(tr)
                if line.original_assignments and not line.assigned_transitions:
                    line.assigned_transitions.append(UNASSIGNED)

    return max_num_assignments, num_conflicts


def build_output(observed_lines: list, weights: dict) -> pd.DataFrame:
    """Build the output DataFrame from observed lines and their assignments.

    Numeric columns (calc_intens, CF, dif_wn_O-C) are stored as actual floats
    with NaN for missing values, so they appear as numbers in Excel/CSV.
    """
    print("  Building output...")
    output_rows = []
    for obs_line in observed_lines:
        unassigned_row = {
            'wn_obs': obs_line.wavenumber,
            'unc_wn_obs': obs_line.wn_uncertainty,
            'obs_intens': obs_line.intensity,
            'char': obs_line.line_character,
            'low_id': '',
            'upp_id': '',
            'calc_intens': np.nan,
            'orig_calc_intens': np.nan,
            'u_calc': np.nan,
            'intens_from_f': np.nan,
            'intens_to_f': np.nan,
            'dif_wn_O-C': np.nan,
            'grade': '',
            'notes1': '',
            'notes2': '',
            'new': '',
            'accepted': np.nan,
            'low_E': np.nan,
            'upp_E': np.nan,
            'rwn': np.nan,
            'weight': np.nan
        }
        if not obs_line.assigned_transitions:
            output_rows.append(unassigned_row)
        else:
            for tr in obs_line.assigned_transitions:
                if tr is UNASSIGNED:
                    output_rows.append(unassigned_row)
                    continue

                calc_wn = tr.calculated_wavenumber
                wn_diff = obs_line.wavenumber - calc_wn
                u_own = obs_line.wn_uncertainty
                # Weight for LOPT input is the factor by which LOPT will multiply 1/u_own**2
                w = weights[id(tr)] * (u_own**2) if tr.accepted is not None and tr.accepted == 1 else 0.0

                output_rows.append({
                    'wn_obs': obs_line.wavenumber,
                    'unc_wn_obs': u_own,
                    'obs_intens': obs_line.intensity,
                    'char': obs_line.line_character,
                    'low_id': tr.lower_level.level_id if tr.lower_level else '',
                    'upp_id': tr.upper_level.level_id if tr.upper_level else '',
                    'calc_intens': tr.calc_intensity,
                    'orig_calc_intens': tr.orig_calc_intensity,
                    'u_calc': tr.u_calc,
                    'intens_from_f': tr.upper_level.intens_from_factor if tr.upper_level else np.nan,
                    'intens_to_f': tr.lower_level.intens_to_factor if tr.lower_level else np.nan,
                    'dif_wn_O-C': wn_diff,
                    'grade': tr.grade if tr.grade else '',
                    'notes1': tr.notes1,
                    'notes2': tr.notes2 if tr.notes2 else '',
                    'new': tr.new if tr.new is not None else '',
                    'accepted': tr.accepted if tr.accepted is not None else np.nan,
                    'low_E': tr.lower_level.energy if tr.lower_level else np.nan,
                    'upp_E': tr.upper_level.energy if tr.upper_level else np.nan,
                    'rwn': tr.upper_level.energy-tr.lower_level.energy if tr.upper_level and tr.lower_level else np.nan,
                    'weight': w
                })

    df = pd.DataFrame(output_rows)

    # Sort: decreasing wavenumber, then increasing grade for ties
    df['_sort_grade'] = df['grade'].apply(lambda g: g if g else 'ZZZZZ')
    df['_sort_rwn'] = df['rwn'].apply(lambda r: r if r else 1.0e10)
    df = df.sort_values(
        by=['wn_obs', '_sort_rwn', '_sort_grade'],
        ascending=[False, False, True]
    ).drop(columns=['_sort_grade','_sort_rwn']).reset_index(drop=True)

    return df


def write_output(df: pd.DataFrame):
    """Write the output DataFrame to Excel and CSV files.

    Excel output gets number formatting for precision display.
    """
    # --- Excel ---
    print(f"  Writing output to {OUTPUT_FILE}...")
    df.to_excel(OUTPUT_FILE, index=False, engine='openpyxl')

    # Apply number formats for display precision in Excel
    wb = openpyxl.load_workbook(OUTPUT_FILE)
    ws = wb.active
    # Map column names to Excel number formats
    col_formats = {
        'wn_obs':  '0.000',
        'unc_wn_obs': '0.000',
        'obs_intens': '0.000',
        'calc_intens': '0.000',
        'orig_calc_intens': '0.000',
        'u_calc': '0.000',
        'intens_from_f': '0.000',
        'intens_to_f': '0.000',
        'dif_wn_O-C': '0.000',
        'low_E': '0.000',
        'upp_E': '0.000',
        'rwn': '0.000',
        'weight': '0.000'
    }
    # Find column indices from header row
    # noinspection PyUnresolvedReferences
    header = {str(cell.value): cell.column for cell in ws[1]}
    for col_name, fmt in col_formats.items():
        if col_name in header:
            col_idx = header[col_name]
            # noinspection PyUnresolvedReferences
            for row in ws.iter_rows(min_row=2, min_col=col_idx,
                                     max_col=col_idx):
                for cell in row:
                    if cell.value is not None:
                        cell.number_format = fmt
    wb.save(OUTPUT_FILE)
    wb.close()

    # --- CSV ---
    print(f"  Writing output to {OUTPUT_CSV}...")
    df.to_csv(OUTPUT_CSV, index=False)

    print(f"  Output written: {len(df)} rows.")


def build_level_transition_lists(levels_dict: dict, observed_lines: list) -> None:
    """Populate from_transitions and to_transitions on each EnergyLevel.

    Uses the same Transition instances attached to lines — no copying.
    Called once before the iterative weeding loop.
    """
    for level in levels_dict.values():
        level.from_transitions = []
        level.to_transitions = []

    for line in observed_lines:
        for t in line.assigned_transitions:
            if t is UNASSIGNED or t.calculated_wavenumber == 0.0:
                continue
            t.upper_level.from_transitions.append(t)
            t.lower_level.to_transitions.append(t)


def reset_weeding_state(observed_lines: list) -> None:
    """Clear decisions from a previous weeding pass so the next iteration starts clean."""
    for line in observed_lines:
        for t in line.assigned_transitions:
            if t is UNASSIGNED:
                continue
            t.accepted = None
            t.notes2 = ""


def snapshot_accepted(observed_lines: list) -> dict:
    """Return {id(transition): accepted} for all non-UNASSIGNED transitions."""
    snap = {}
    for line in observed_lines:
        for t in line.assigned_transitions:
            if t is UNASSIGNED:
                continue
            snap[id(t)] = t.accepted
    return snap


def snapshot_levels(levels_list: list) -> dict:
    """Return {id(level): energy} for all levels."""
    snap = {}
    for lev in levels_list:
        snap[id(lev)] = lev.energy
    return snap


def count_changes(prev: dict, curr: dict) -> int:
    """Count transitions whose accepted value differs between two snapshots."""
    n = 0
    for key in prev:
        if key in curr and prev[key] != curr[key]:
            n += 1
    return n


def _wm_red_chi(values: list, uncertainties: list) -> tuple:
    """Weighted mean of a set of values with uncertainties, with chi-squared inflation.
       Returns: (weighted mean, wm uncertainty). If sum of weights is too small, returns (None, None).
    """
    # Weighted mean : w_i = 1/u_calc_i^2, r_i = ln(I_obs / I_calc_orig)
    sum_w = 0.0
    sum_wr = 0.0
    n_val = len(values)
    for i in range(n_val):
        w = 1.0 / (uncertainties[i] ** 2)
        r = values[i]
        sum_w += w
        sum_wr += w * r

    if sum_w < 1e-30:
        return None, None

    wm = sum_wr / sum_w
    u_stat = 1.0 / math.sqrt(sum_w)

    # Chi-squared inflation
    chi2 = 0.0
    for i in range(n_val):
        w = 1.0 / (uncertainties[i] ** 2)
        r = values[i]
        chi2 += w * (r - wm) ** 2

    nu = n_val - 1
    chi2_nu = chi2 / nu if nu > 0 else 1.0
    u_wm = u_stat * math.sqrt(max(1.0, chi2_nu))
    return wm, u_wm

def _wm_mandel_paule(values: list, uncertainties: list) -> tuple:
    from statistics import mandel_paule
    wm, u_wm, _ = mandel_paule(values, uncertainties)
    return wm, u_wm


def _compute_factor_for_transition_list(transitions: list, min_n: int=5) -> tuple:
    """Compute a weighted-mean intensity adjustment factor from a list of transitions.

    Qualifying transitions: accepted == 1, orig_calc_intensity and orig_u_calc > 0,
    and the parent line has exactly one accepted transition (unblended).

    Returns (factor, u_factor, n_qualifying).
    If n_qualifying < min_n, returns (0.0, 0.0, n_qualifying).
    """
    CLAMP = 5.0

    qualifying = []
    for t in transitions:
        if t.accepted != 1:
            continue
        if t.orig_calc_intensity is None or t.orig_u_calc is None:
            continue
        if t.orig_calc_intensity <= 0 or t.orig_u_calc <= 0:
            continue
        line = t.assigned_to
        if line is None:
            continue
        n_accepted_on_line = sum(1 for tt in line.assigned_transitions if tt.accepted == 1)
        if n_accepted_on_line != 1:
            continue
        qualifying.append(t)

    n = len(qualifying)
    if n < min_n:
        return 0.0, 0.0, n

    # Weighted mean: w_i = 1/u_calc_i^2, r_i = ln(I_obs / I_calc_orig)
    r_values = []
    r_uncertainties = []
    for t in qualifying:
        r_uncertainties += [t.orig_u_calc]
        r = math.log(t.assigned_to.intensity / t.orig_calc_intensity)
        r_values += [r]

    use_mandel_paule = True
    if use_mandel_paule:
        factor, u_factor = _wm_mandel_paule(r_values, r_uncertainties)
    else:
        factor, u_factor = _wm_red_chi(r_values, r_uncertainties)
    if factor is None:
        return 0.0, 0.0, n

    # Magnitude clamping
    factor = max(-CLAMP, min(CLAMP, factor))

    return factor, u_factor, n


def compute_intensity_factors(levels_dict: dict, alpha: float = 0.5,
                              prev_factors: dict = None, min_n: int=5) -> dict:
    """Compute per-level intensity adjustment factors from accepted transitions.

    The same procedure is used for adjusting calculated intensities for transitions from each level
    and for transitions to each level.

    For transitions from a given level to lower levels, the adjustment is physically justified
    by possible deviations of level populations from the Boltzmann distribution used in the intensity calculation.
    These adjustments can be either positive or negative.

    For transitions to a given level from upper levels, one part of the adjustment is physically explained by
    self-absorption, which effectively decreases the observed intensity down to highly populated metastable levels.
    Another effect that could explain under- or over-prediction of the intensity is the possible presence
    of systematic errors in the calculated transition probabilities.

    Uses orig_calc_intensity/orig_u_calc (raw theoretical values) to prevent
    feedback amplification across iterations.

    Args:
        levels_dict: dict of level_id -> EnergyLevel
        alpha: damping coefficient (0..1). 1.0 = no damping, 0.5 = blend
               equally with previous. Ignored on first iteration (prev_factors=None).
        prev_factors: {level_id: (from_factor, to_factor)} from previous iteration.
                      None on first iteration.
        min_n: Minimum number of qualifying transitions required to compute a factor.

    Returns a diagnostics dictionary including per-level details.
    """
    n_from_computed = 0
    n_from_skipped = 0
    n_to_computed = 0
    n_to_skipped = 0
    from_factors = []
    to_factors = []
    level_details = {}  # level_id -> (n_from, from_factor, n_to, to_factor)

    for lid, level in levels_dict.items():
        # Upper-level (emission) correction
        raw_from, u_from, n_from = _compute_factor_for_transition_list(level.from_transitions, min_n)
        if prev_factors is not None and lid in prev_factors:
            prev_from = prev_factors[lid][0]
            applied_from = alpha * raw_from + (1.0 - alpha) * prev_from
        else:
            applied_from = raw_from
        level.intens_from_factor = applied_from
        level.u_intens_from_factor = u_from
        if n_from >= 5:
            n_from_computed += 1
            if applied_from != 0.0:
                from_factors.append(applied_from)
        else:
            n_from_skipped += 1

        # Lower-level correction
        raw_to, u_to, n_to = _compute_factor_for_transition_list(level.to_transitions, min_n)
        if prev_factors is not None and lid in prev_factors:
            prev_to = prev_factors[lid][1]
            applied_to = alpha * raw_to + (1.0 - alpha) * prev_to
        else:
            applied_to = raw_to
        level.intens_to_factor = applied_to
        level.u_intens_to_factor = u_to
        if n_to >= 5:
            n_to_computed += 1
            if applied_to != 0.0:
                to_factors.append(applied_to)
        else:
            n_to_skipped += 1

        level_details[lid] = (n_from, applied_from, n_to, applied_to)

    return {
        'n_from_computed': n_from_computed,
        'n_from_skipped': n_from_skipped,
        'n_to_computed': n_to_computed,
        'n_to_skipped': n_to_skipped,
        'from_factors': from_factors,
        'to_factors': to_factors,
        'level_details': level_details,
    }


def print_factor_diagnostics(iteration: int, diag: dict):
    """Print summary statistics of intensity adjustment factors."""
    print(f"  --- Factor diagnostics (iteration {iteration}) ---")
    for label, key_comp, key_skip, factors in [
        ("FROM (upper)", 'n_from_computed', 'n_from_skipped', diag['from_factors']),
        ("TO   (lower)", 'n_to_computed', 'n_to_skipped', diag['to_factors']),
    ]:
        n_comp = diag[key_comp]
        n_skip = diag[key_skip]
        print(f"    {label}: {n_comp} levels with factors, {n_skip} skipped (n<5)")
        if factors:
            arr = np.array(factors)
            print(f"      mean={arr.mean():.3f}  median={np.median(arr):.3f}  "
                  f"std={arr.std():.3f}  min={arr.min():.3f}  max={arr.max():.3f}")
            if label.startswith("TO"):
                n_pos = int(np.sum(arr > 0))
                n_neg = int(np.sum(arr < 0))
                print(f"      sign distribution: {n_pos} positive, {n_neg} negative")


def print_factor_stability(iteration: int, level_details: dict, prev_factors: dict|None):
    """Log per-level factor changes between consecutive iterations.

    Prints only levels where either from_factor or to_factor shifted by more
    than 0.05 (log units) relative to the previous iteration.
    """
    if prev_factors is None:
        return
    THRESHOLD = 0.05
    changed = []
    for lid, (n_from, from_f, n_to, to_f) in level_details.items():
        if lid not in prev_factors:
            continue
        prev_from, prev_to = prev_factors[lid]
        delta_from = abs(from_f - prev_from)
        delta_to = abs(to_f - prev_to)
        if delta_from > THRESHOLD or delta_to > THRESHOLD:
            changed.append((lid, n_from, prev_from, from_f, delta_from,
                            n_to, prev_to, to_f, delta_to))
    if changed:
        print(f"  --- Factor stability (iter {iteration}): "
              f"{len(changed)} levels shifted >0.05 ---")
        changed.sort(key=lambda x: max(x[4], x[8]), reverse=True)
        for row in changed[:20]:
            lid, nf, pf, cf, df, nt, pt, ct, dt = row
            print(f"    {lid:12s}  FROM n={nf:3d} {pf:+.3f}->{cf:+.3f} Δ{df:.3f}"
                  f"   TO n={nt:3d} {pt:+.3f}->{ct:+.3f} Δ{dt:.3f}")
        if len(changed) > 20:
            print(f"    ... and {len(changed) - 20} more")
    else:
        print(f"  --- Factor stability (iter {iteration}): all levels stable (<0.05) ---")


def detect_oscillations(accepted_history: list, trans_map: dict):
    """Identify transitions that flip acceptance state A→B→A across last 3 snapshots.

    Args:
        accepted_history: list of snapshot dicts [{id(t): accepted}, ...], oldest first.
        trans_map: {id(t): Transition} for reverse lookup.
    """
    if len(accepted_history) < 3:
        return
    s1, s2, s3 = accepted_history[-3], accepted_history[-2], accepted_history[-1]
    oscillating = []
    for tid in s1:
        if tid not in s2 or tid not in s3:
            continue
        v1, v2, v3 = s1[tid], s2[tid], s3[tid]
        if v1 == v3 and v1 != v2:
            oscillating.append((tid, v1, v2, v3))
    if oscillating:
        print(f"  --- Oscillation detection: {len(oscillating)} transitions flip A→B→A in last 3 iterations ---")
        for tid, v1, v2, v3 in oscillating[:20]:
            t = trans_map.get(tid)
            if t is None:
                continue
            upper = t.upper_level.level_id if t.upper_level else "?"
            lower = t.lower_level.level_id if t.lower_level else "?"
            wn = getattr(t, 'calculated_wavenumber', '?')
            print(f"    {upper} -> {lower}  wn={wn}  "
                  f"accepted: {v1}->{v2}->{v3}")
        if len(oscillating) > 20:
            print(f"    ... and {len(oscillating) - 20} more")
    else:
        print(f"  --- Oscillation detection: no A→B→A flips in last 3 iterations ---")


def apply_intensity_adjustments(observed_lines: list) -> None:
    """Overwrite calc_intensity and u_calc on every Transition with adjusted values.

    I_calc_adj = orig_calc_intensity * exp(intens_from_factor + intens_to_factor)
    u_calc_adj = sqrt(orig_u_calc² + u_from² + u_to²)

    weed_assignments_line() sees adjusted values transparently.
    """
    for line in observed_lines:
        for t in line.assigned_transitions:
            if t is UNASSIGNED:
                continue
            if t.orig_calc_intensity is None:
                t.calc_intensity = None
                t.u_calc = None
                continue
            from_f = t.upper_level.intens_from_factor
            to_f = t.lower_level.intens_to_factor
            t.calc_intensity = t.orig_calc_intensity * math.exp(from_f + to_f)

            u_from = t.upper_level.u_intens_from_factor
            u_to = t.lower_level.u_intens_to_factor
            t.u_calc = math.sqrt(t.orig_u_calc ** 2 + u_from ** 2 + u_to ** 2)


# noinspection PyTypeChecker,PyUnresolvedReferences
def calculate_sum_uncertainty(transitions: List[Transition]):
    """Calculates the propagated log-uncertainty of a sum of log-normal intensities."""
    valid = [t for t in transitions if t.calc_intensity is not None and t.u_calc is not None]
    if not valid: return 0.0

    total_i = sum(t.calc_intensity for t in valid)
    if total_i == 0: return 0.0

    # Square root of sum of (I * u_ln)^2 divided by Total_I
    sum_sq_errors = sum((t.calc_intensity * t.u_calc) ** 2 for t in valid)
    return np.sqrt(sum_sq_errors) / total_i


def check_internal_significance(trans: Transition, accepted_already: List[Transition], S: float = 2.0):
    """
    Checks if a legacy transition is statistically guaranteed to be < 4% of I_cum.
    Returns True if it should be rejected as insignificant.
    """
    if not accepted_already or trans.calc_intensity is None:
        return False

    i_cum = sum(t.calc_intensity for t in accepted_already if t.calc_intensity)
    if i_cum == 0: return False

    # 1. Propagated uncertainty of the already accepted sum
    u_ln_cum = calculate_sum_uncertainty(accepted_already)

    # 2. Combined uncertainty of the ratio (Candidate / I_cum)
    u_ratio = np.sqrt((trans.u_calc or 0) ** 2 + u_ln_cum ** 2)

    # 3. Log-difference from the 4% (0.04) threshold
    # L = ln(I_trans) - ln(0.04 * I_cum)
    l_val = np.log(trans.calc_intensity) - np.log(0.04 * i_cum)

    # 4. Z-score (using legacy Z=3.0)
    z_internal = l_val / (S * u_ratio)

    if z_internal < -3.0:
        return True  # Guarantee reached: Reject
    return False  # Not guaranteed: Keep/Retain


# noinspection PyTypeChecker
def weed_assignments_line(line: SpectralLine, u_obs_ln: float=np.log(2), S: float = 2.0, blacklist: set = None):
    """
    Evaluate assignments for a line in steps.
    Step 1: Only make a decision if the case is statistically 'clear-cut'.
    Step 2: Decide unclear cases with additional center-of-gravity (cog) and line character considerations.
    Step 3: Reject all undecided candidates.
    """
    # Reference for the statistical theory:
    # Fisher, R.A. (1992). Statistical Methods for Research Workers. In: Kotz, S., Johnson, N.L. (eds)
    # Breakthroughs in Statistics. Springer Series in Statistics. Springer, New York, NY.
    # https://doi.org/10.1007/978-1-4612-4380-9_6

    # Sorting by I_calc (desc), then Ritz sigma (asc)
    line.assigned_transitions.sort(key=lambda t: (
        t.calc_intensity is None,
        -(t.calc_intensity or 0),
        abs(line.wavenumber - t.calculated_wavenumber) / line.wn_uncertainty
    ))

    # --- tiny local helpers (behavior-preserving) ---
    def _accepted_already() -> list[Transition]:
        return [t for t in line.assigned_transitions if t.accepted == 1]

    def _i_cum(transitions_list: list[Transition]) -> float:
        # NOTE: preserve "is not None" semantics (0.0 contributes, None doesn't)
        return sum(t.calc_intensity for t in transitions_list if t.calc_intensity is not None)

    def _current_group_stats(accepted_list: list[Transition], candidate: Transition):
        # NOTE: preserve truthiness semantics (0.0 is excluded, None excluded)
        cur_group = accepted_list
        if not candidate in accepted_list:
            cur_group = cur_group + [candidate]
        i_tot_theory = _i_cum(cur_group)
        u_ln_theo_group = calculate_sum_uncertainty(cur_group)
        return cur_group, i_tot_theory, u_ln_theo_group

    def _ritz_sigma(candidate: Transition) -> float:
        return abs(line.wavenumber - candidate.calculated_wavenumber) / line.wn_uncertainty

    def _u_sys(candidate: Transition):
        # NOTE: preserve truthiness semantics (0.0 behaves as missing intensity)
        return np.sqrt((candidate.u_calc or 0) ** 2 + u_obs_ln ** 2) if candidate.calc_intensity else None

    def _lower_id_resonance(candidate: Transition) -> bool:
        return (candidate.lower_level is not None and
                str(candidate.lower_level.level_id).endswith('001'))

    def _any_i_none(group: list[Transition]):
        return any(t.calc_intensity is None for t in group)
    
    def _all_decided(group: list[Transition]):
        not_decided = any(t.accepted is None for t in group)
        return not not_decided

    def _any_new(group: list[Transition]):
        return any(t.new == 1 for t in group)

    def _cog(group: list[Transition]) -> float:
        if not group: return 999.0
        any_none = _any_i_none(group)
        i_cum_group = _i_cum(group)
        if not any_none and i_cum_group > 0:
            return sum(t.calculated_wavenumber * t.calc_intensity for t in group) / i_cum_group
        else:
            # Fallback to unweighted average if any I_calc is missing
            return sum(t.calculated_wavenumber for t in group) / len(group)

    def _effective_spread_combined(group: list[Transition],
                                   obs_wn_unc: float,
                                   temperature_eV: float,
                                   atomic_mass_u: float) -> float:
        """
        Effective spread normalized by combined resolution width.

        The observed line profile is a convolution of the Doppler (thermal)
        profile and the instrumental profile.  The effective resolution is
        sqrt(doppler_sigma² + wn_uncertainty²), so at low wavenumber where
        Doppler widths are small the measurement uncertainty dominates
        (recovering sigma-spread behavior), and at high wavenumber the
        Doppler width dominates.

        Returns dimensionless value:
            effective_spread / sqrt(doppler_sigma² + wn_uncertainty²)
        """

        if not group: return 999.0
        if temperature_eV <= 0: temperature_eV = 1.0
        if atomic_mass_u <= 0: atomic_mass_u = 100.0

        # --- constants ---
        eV_to_J = 1.602176634e-19
        u_to_kg = 1.66053906660e-27
        c = 2.99792458e10  # cm/s

        # --- Doppler sigma in cm^-1 ---
        T_J = temperature_eV * eV_to_J
        m_kg = atomic_mass_u * u_to_kg

        # sqrt(kT / m) / c
        doppler_sigma_factor = (T_J / m_kg) ** 0.5 / c

        # Use representative wavenumber (CoG is best)
        cog = _cog(group)
        doppler_sigma = cog * doppler_sigma_factor

        # --- combined resolution width ---
        effective_resolution = math.sqrt(doppler_sigma ** 2 + obs_wn_unc ** 2)

        # --- effective spread in cm^-1 ---
        any_none = _any_i_none(group)
        sum_i = _i_cum(group)

        if not any_none and sum_i > 0:
            # noinspection PyUnresolvedReferences
            eff_spread = sum(
                t.calc_intensity * abs(t.calculated_wavenumber - cog)
                for t in group
            ) / sum_i
        else:
            # If no theoretical intensities, apply twice greater effective spread to avoid
            # insufficiently grounded decisions
            max_wn = max(t.calculated_wavenumber for t in group)
            min_wn = min(t.calculated_wavenumber for t in group)
            eff_spread = max_wn - min_wn

        # --- normalize ---
        if effective_resolution > 0:
            return eff_spread / effective_resolution
        else:
            return 999.0

    def _fudge_factor_for_asym_intensity(candidates: list[Transition],
                                         i_theo: float,
                                         line_intensity: float,
                                         i_threshold) -> float:
        if i_theo is None or i_theo == 0 or line_intensity == 0:
            return 1.0
        ratio = i_theo / line_intensity
        if ratio > i_threshold:
            # Assume that candidates[0] is the strongest of the candidates (has max calc_intensity)
            return 2.5 if _lower_id_resonance(candidates[0]) else 2.0
        return 1.0

    def _apply_relative_intensity_filter(step_label: str,
                                         candidate: Transition,
                                         accepted_list: list,
                                         i_cum_sum: float,
                                         i_tot_theo: float,
                                         u_ln_theo_gr: float,
                                         allow_insufficient_evidence_new_when_i_cum_zero: bool):
        """
        Returns (decided: bool, accepted_value: Optional[int], notes2: Optional[str]).
        Caller is responsible for counter updates and `continue`.
        """
        # --- Pre-check: Relative Intensity Filter (The 10% Rule) ---
        if candidate.calc_intensity is not None:
            if i_cum_sum > 0:
                rel_contribution = candidate.calc_intensity / i_cum_sum

                # If new and contribution is tiny, reject immediately
                if candidate.new == 1 and rel_contribution < 0.1:
                    return True, 0, f"{step_label} Rejected New: Insignificant intensity (<10% of I_cum)"

                # If legacy and tiny, apply a conservative uncertainty test
                if candidate.new == 0:
                    if rel_contribution >= 0.04:
                        return True, 1, f"{step_label} Accepted Legacy: Significant intensity (>=4% of I_cum)"

                    # Tightened: very small share of I_cum — reject before internal-significance / z-test
                    if rel_contribution < 0.02:
                        return True, 0, (
                            f"{step_label} Rejected Legacy: calc contribution too weak vs I_cum (<2% of I_cum)"
                        )

                    if check_internal_significance(candidate, accepted_list):
                        return True, 0, f"{step_label} Rejected Legacy: Statistically guaranteed <4% calc contribution"

                    # Relative contribution checked; now compare with observed intensity.
                    # Combine theoretical uncertainty with observed intensity uncertainty
                    u_sys_tot = np.sqrt(u_ln_theo_gr ** 2 + u_obs_ln ** 2)
                    # Z-test on the TOTAL sum vs the observation (asymmetric: too-weak vs too-strong)
                    ln_ratio_tot = np.log(i_tot_theo / line.intensity)
                    d_ln_abs = abs(ln_ratio_tot)
                    z_score = d_ln_abs / (S * u_sys_tot)
                    if ln_ratio_tot > 0:
                        # Predicted sum too strong vs Iobs: relax rejection threshold
                        z_reject = 2.0 * _fudge_factor_for_asym_intensity([candidate],
                                                                          i_tot_theo, line.intensity, 1.0)
                        if z_score <= 2.0:
                            note = f"{step_label} Accepted Legacy: Ic stat significant & compatible with Iobs (asym: sum too strong)"
                            if _lower_id_resonance(candidate):
                                note += ", resonance lower"
                            return True, 1, note
                        if z_score > z_reject:
                            note = f"{step_label} Rejected Legacy: Ic statistically incompatible with Iobs (asym: sum too strong, z>{z_reject:.1f})"
                            return True, 0, note
                    elif ln_ratio_tot < 0:
                        # Predicted sum too weak: stricter rejection
                        if z_score <= 2.0:
                            return True, 1, f"{step_label} Accepted Legacy: Ic stat significant & compatible with Iobs (asym: sum too weak)"
                        if z_score > 2.5:
                            return True, 0, f"{step_label} Rejected Legacy: Ic statistically incompatible with Iobs (asym: sum too weak, z>2.5)"
                    else:
                        if z_score <= 2.0:
                            return True, 1, f"{step_label} Accepted Legacy: Ic stat significant & compatible with Iobs"
                        if z_score > 3.0:
                            return True, 0, f"{step_label} Rejected Legacy: Ic statistically incompatible with Iobs"
                    return False, None, None
            elif candidate.new == 1 and i_cum_sum == 0.0: # This is the top-calc-intensity candidate. Step must be 1
                # Compare Icalc with Iobs
                # noinspection PyUnresolvedReferences
                u_sys_tot = np.sqrt(candidate.u_calc ** 2 + u_obs_ln ** 2)
                # Z-test on the TOTAL sum vs the observation (asymmetric: too-weak vs too-strong)
                ln_ratio_tot = np.log(candidate.calc_intensity / line.intensity)
                d_ln_abs = abs(ln_ratio_tot)
                z_score = d_ln_abs / (S * u_sys_tot)
                if ln_ratio_tot > 0:
                    # Predicted sum too strong vs Iobs: relax rejection threshold
                    z_reject = 1.0 * _fudge_factor_for_asym_intensity([candidate],
                                                                      i_tot_theo, line.intensity, 1.0)
                    if z_score <= 1.0:
                        note = f"{step_label} Accepted New: Ic stat significant & compatible with Iobs (asym: sum too strong)"
                        if _lower_id_resonance(candidate):
                            note += ", resonance lower"
                        return True, 1, note
                    if z_score > z_reject:
                        note = f"{step_label} Rejected new: Ic statistically incompatible with Iobs (asym: sum too strong, z>{z_reject:.1f})"
                        return True, 0, note
                elif ln_ratio_tot < 0:
                    # Predicted sum too weak: stricter rejection
                    if z_score <= 1.0:
                        return True, 1, f"{step_label} Accepted New: Ic stat significant & compatible with Iobs (asym: sum too weak)"
                    if z_score > 1.25:
                        return True, 0, f"{step_label} Rejected New: Ic statistically incompatible with Iobs (asym: sum too weak, z>1.25)"
                else:
                    if z_score <= 1.0:
                        return True, 1, f"{step_label} Accepted New: Ic stat significant & compatible with Iobs"
                    if z_score > 1.5:
                        return True, 0, f"{step_label} Rejected New: Ic statistically incompatible with Iobs"
                return False, None, None
            elif allow_insufficient_evidence_new_when_i_cum_zero and candidate.new == 1:  # i_cum == 0, new assignment
                return True, 0, f"{step_label} Rejected New: Insufficient evidence for acceptance"

        return False, None, None

    def _broadening_multiplier() -> float:
        ch = (line.line_character or "").strip()
        if ch == "" or ch == "**":
            return 1.0
        if ch == "h":
            return 1.5
        return 2.0

    def _try_pair_stage(include_broadening: bool,
                        n_decided_step2: int,
                        n_undecided: int) -> tuple[bool, int, int]:
        pool_not_rejected = [t for t in line.assigned_transitions
                             if t.accepted != 0 and not (blacklist and id(t) in blacklist)]
        n_pool = len(pool_not_rejected)
        pair_pass_made = False
        progress = False

        for ii in range(n_pool):
            if pair_pass_made:
                break

            for jj in range(ii + 1, n_pool):
                ta, tb = pool_not_rejected[ii], pool_not_rejected[jj]

                # Ignore pairs that cannot change anything
                if _all_decided([ta, tb]):
                    continue

                any_new = _any_new([ta, tb])
                any_int_none = _any_i_none([ta, tb])

                # Use pair thresholds aligned with single-candidate policy:
                # - any new in pair -> stricter
                # - all legacy -> looser
                joint_cog_max_sig = 2.5 if any_new else 4.0
                spread_max_sig = 2.0 if any_new else 4.0
                sum_int_z_max = 1.0 if any_new else 2.0

                m_broadening = _broadening_multiplier() if include_broadening else 1.0

                # Joint CoG and effective spread
                joint_wavenumber = _cog([ta, tb])
                sum_i = _i_cum([ta, tb])
                effective_spread_sig = _effective_spread_combined([ta, tb],
                                                                 line.wn_uncertainty, T, ATOMIC_MASS)

                joint_cog_sig = abs(line.wavenumber - joint_wavenumber) / line.wn_uncertainty

                if joint_cog_sig > joint_cog_max_sig:
                    relaxed_joint_cog_sig = joint_cog_sig / m_broadening
                    if relaxed_joint_cog_sig > joint_cog_max_sig:
                        continue
                    broadening_helped = True
                else:
                    broadening_helped = False

                # Spread penalty: relaxed by broadening multiplier for broadened lines
                if effective_spread_sig > spread_max_sig * m_broadening:
                    continue

                # Joint intensity check only when both intensities are available
                if not any_int_none:
                    u_ln_pair = calculate_sum_uncertainty([ta, tb])
                    u_sys_tot = np.sqrt(u_ln_pair ** 2 + u_obs_ln ** 2)
                    if u_sys_tot <= 0:
                        continue

                    # Fudge factor uses full group (accepted + pair) for total predicted intensity
                    acc_already = _accepted_already()
                    full_group = acc_already + [t for t in (ta, tb) if t.accepted != 1]
                    i_full = _i_cum(full_group)
                    full_group_sorted = sorted(full_group, key=lambda t: -(t.calc_intensity or 0))

                    z_sum_intens = abs(np.log(sum_i / line.intensity)) / (S * u_sys_tot)
                    FF = _fudge_factor_for_asym_intensity(full_group_sorted,
                                                         i_full, line.intensity, 1.0)
                    if z_sum_intens > sum_int_z_max * FF:
                        continue
                else:
                    z_sum_intens = None

                for t in (ta, tb):
                    if t.accepted is None:
                        t.accepted = 1
                        n_decided_step2 += 1
                        n_undecided -= 1
                        progress = True

                        pair_type = "mixed/new" if any_new else "legacy"
                        msg = (
                            f"Step2 Accepted (pair, {pair_type}): "
                            f"joint CoG sigma={joint_cog_sig:.2f}; "
                            f"eff. spread sigma={effective_spread_sig:.2f}"
                        )
                        if z_sum_intens is not None:
                            msg += f"; sum I z={z_sum_intens:.2f}"
                            if z_sum_intens > sum_int_z_max:
                                msg += "; Ok for too-strong Ic"
                        if broadening_helped:
                            msg += "; accepted: Ritz Ok for broadened line"
                        t.notes2 = msg

                pair_pass_made = True
                break
        return progress, n_decided_step2, n_undecided

    def _try_triple_stage(include_broadening: bool,
                        n_decided_step2: int,
                        n_undecided: int) -> tuple[bool, int, int]:
        pool_nr = [t for t in line.assigned_transitions
                   if t.accepted != 0 and not (blacklist and id(t) in blacklist)]
        npool = len(pool_nr)
        triple_pass_made = False
        progress = False

        for ii in range(npool):
            if triple_pass_made:
                break
            for jj in range(ii + 1, npool):
                if triple_pass_made:
                    break
                for kk in range(jj + 1, npool):
                    ta, tb, tc = pool_nr[ii], pool_nr[jj], pool_nr[kk]

                    if _all_decided([ta, tb, tc]):  # Ignore triplets that cannot change anything
                        continue

                    any_new = _any_new([ta, tb, tc])
                    any_int_none = _any_i_none([ta, tb, tc])

                    joint_cog_max_sigma = 2.5 if any_new else 4.0
                    spread_max_sigma = 2.0 if any_new else 4.0
                    sum_i_z_max = 1.0 if any_new else 2.0

                    m_broad = _broadening_multiplier() if include_broadening else 1.0

                    # Joint CoG and effective spread
                    joint_wn = _cog([ta, tb, tc])
                    sum_i = _i_cum([ta, tb, tc])
                    effective_spread_sig = _effective_spread_combined([ta, tb, tc],
                                                                     line.wn_uncertainty, T, ATOMIC_MASS)

                    joint_cog_sigma = abs(line.wavenumber - joint_wn) / line.wn_uncertainty

                    if joint_cog_sigma > joint_cog_max_sigma:
                        relaxed_joint_cog_sigma = joint_cog_sigma / m_broad
                        if relaxed_joint_cog_sigma > joint_cog_max_sigma:
                            continue
                        broadened_helped = True
                    else:
                        broadened_helped = False

                    # Spread penalty: relaxed by broadening multiplier for broadened lines
                    if effective_spread_sig > spread_max_sigma * m_broad:
                        continue

                    if not any_int_none:
                        u_ln_triple = calculate_sum_uncertainty([ta, tb, tc])
                        u_sys_total = np.sqrt(u_ln_triple ** 2 + u_obs_ln ** 2)
                        if u_sys_total <= 0:
                            continue

                        # Fudge factor uses full group (accepted + triple) for total predicted intensity
                        acc_already = _accepted_already()
                        full_group = acc_already + [t for t in (ta, tb, tc) if t.accepted != 1]
                        i_full = _i_cum(full_group)
                        full_group_sorted = sorted(full_group, key=lambda t: -(t.calc_intensity or 0))

                        z_sum_i = abs(np.log(sum_i / line.intensity)) / (S * u_sys_total)
                        FF = _fudge_factor_for_asym_intensity(full_group_sorted,
                                                             i_full, line.intensity, 1.0)
                        if z_sum_i > sum_i_z_max * FF:
                            continue
                    else:
                        z_sum_i = None

                    for t in (ta, tb, tc):
                        if t.accepted is None:
                            t.accepted = 1
                            n_decided_step2 += 1
                            n_undecided -= 1
                            progress = True

                            triple_type = "mixed/new" if any_new else "legacy"
                            msg = (
                                f"Step2 Accepted (triple, {triple_type}): "
                                f"joint CoG sigma={joint_cog_sigma:.2f}; "
                                f"eff. spread sigma={effective_spread_sig:.2f}"
                            )
                            if z_sum_i is not None:
                                msg += f"; sum I z={z_sum_i:.2f}"
                                if z_sum_i > sum_i_z_max:
                                    msg += "; Ok for too-strong Ic"
                            if broadened_helped:
                                msg += "; accepted: Ritz Ok for broadened line"
                            t.notes2 = msg
                    triple_pass_made = True
                    break
        return progress, n_decided_step2, n_undecided

    num_undecided = len(line.assigned_transitions)
    step1_done = False
    step2_done = False
    while num_undecided > 0:
        undecided = [t for t in line.assigned_transitions if t.accepted is None]
        num_undecided = len(undecided)

        if not step1_done:
            num_decided_step1 = 0
            for rank, trans in enumerate(undecided, start=1):
                if blacklist and id(trans) in blacklist:
                    continue
                # 1. Initialize cumulative values from already accepted (accepted == 1) transitions
                accepted_already = _accepted_already()

                # Track cumulative intensity and weighted wavenumber sum
                i_cum = _i_cum(accepted_already)

                # Ritz mismatch in sigma units
                ritz_sigma = _ritz_sigma(trans)
                u_sys = _u_sys(trans)

                if trans.new == 1:  # Test ritz_sigma first
                    if ritz_sigma > 4.0:
                        trans.accepted = 0
                        num_undecided = num_undecided - 1
                        trans.notes2 = "Step1 Rejected New: Strong Ritz mismatch"
                        continue


                # Calculate total intensity and its uncertainty of what's already accepted + the current candidate
                _current_group, i_total_theory, u_ln_theory_group = _current_group_stats(accepted_already, trans)

                decided, accepted_value, notes2 = _apply_relative_intensity_filter(
                    "Step1",
                    trans,
                    accepted_already,
                    i_cum,
                    i_total_theory,
                    u_ln_theory_group,
                    allow_insufficient_evidence_new_when_i_cum_zero=False
                )

                if decided:
                    trans.accepted = accepted_value
                    num_decided_step1 = num_decided_step1 + 1
                    num_undecided = num_undecided - 1
                    trans.notes2 = notes2
                    continue

                # --- BRANCH 1: Legacy Assignments (is_new == 0) ---
                if trans.new == 0:
                    if trans.calc_intensity is not None:
                        d_signed = np.log(trans.calc_intensity/line.intensity)
                        # F is fudge factor to relax rejection for strong predictions
                        F = 2.0 if d_signed > np.log(1.0/5) else 1.0
                        d_abs = abs(d_signed)
                        resonance = _lower_id_resonance(trans)
                        if d_signed > 0:
                            thresh_reject = (4.0 if resonance else 3.5) * S * F * u_sys
                        elif d_signed < 0:
                            thresh_reject = 2.5 * S * F * u_sys
                        else:
                            thresh_reject = 3.0 * S * F * u_sys
                        if d_abs > thresh_reject:
                            trans.accepted = 0
                            num_undecided = num_undecided - 1
                            if d_signed < 0:
                                tag = "asym: Ic too weak vs Iobs"
                            elif d_signed > 0:
                                tag = "asym: Ic too strong vs Iobs"
                                if resonance:
                                    tag += ", resonance lower"
                            else:
                                tag = "Ic vs Iobs mismatch"
                            trans.notes2 = (
                                f"Step1 Rejected Legacy: Int mismatch ({d_abs / u_sys:.1f} sigma; {tag})"
                            )
                        elif ritz_sigma <= 4.0:
                            if trans.calc_intensity / line.intensity >= 0.04:
                                trans.accepted = 1
                                num_undecided = num_undecided - 1
                                trans.notes2 = "Step1 Accepted Legacy: Ritz match, Ic/Iobs >= 0.04"
                            else:
                                compat_lim = (1.5 * S * F * u_sys) if d_signed < 0 else (2.0 * S * F * u_sys)
                                if d_abs <= compat_lim:
                                    trans.accepted = 1
                                    num_undecided = num_undecided - 1
                                    suf = " (asym: tighter compat when Ic too weak)" if d_signed < 0 else ""
                                    trans.notes2 = "Step1 Accepted Legacy: Ritz Ok, Ic stat compatible with Iobs" + suf
                                else:
                                    trans.notes2 = "Retained Legacy: Pending further checks"
                        else:
                            # We don't accept yet; we leave it as None (Retained/Undecided)
                            # to see if a CoG or a better 'New' transition competes with it later.
                            trans.notes2 = "Retained Legacy: Pending further checks"
                    else:
                        # No intensity for legacy? Generally keep unless Ritz is impossible.
                        if ritz_sigma > 5.0:
                            trans.accepted = 0
                            num_undecided = num_undecided - 1
                            trans.notes2 = "Step1 Rejected Legacy: Strong Ritz mismatch"
                        elif ritz_sigma <= 2.0:
                            trans.accepted = 1
                            num_undecided = num_undecided - 1
                            trans.notes2 = "Step1 Accepted Legacy: Good Ritz match"
                        elif ritz_sigma <= 4.0:
                            trans.accepted = 1
                            num_undecided = num_undecided - 1
                            trans.notes2 = "Step1 Accepted Legacy: Weak Ritz match"
                        else:
                            trans.notes2 = "Retained Legacy: No Ic, weak Ritz match"

                # --- BRANCH 2: New Assignments (new == 1) ---
                elif trans.new == 1:
                    if ritz_sigma > 4.0:
                        trans.accepted = 0
                        num_undecided = num_undecided - 1
                        trans.notes2 = "Step1 Rejected New: Strong Ritz mismatch"
                    else:
                        if trans.calc_intensity is not None:
                            d_signed = np.log(trans.calc_intensity/line.intensity)
                            z_thresh = 1.0 if ritz_sigma <= 2.5 else 0.5
                            F = _fudge_factor_for_asym_intensity([trans],trans.calc_intensity,line.intensity, 
                                                                 0.2)
                            # F = 2.0 if d_signed > np.log(1.0/5) else 1.0
                            d_abs = abs(d_signed)
                            resonance = _lower_id_resonance(trans)
                            if d_signed > 0:
                                thresh_reject = (4.0 if resonance else 3.5) * S * F * u_sys
                            elif d_signed < 0:
                                thresh_reject = 2.5 * S * F * u_sys
                            else:
                                thresh_reject = 3.0 * S * F * u_sys

                            if d_abs <= (z_thresh * S * F * u_sys):
                                trans.accepted = 1
                                num_undecided = num_undecided - 1
                                trans.notes2 = f"Step1 Accepted New: Ritz Ok, Fair Ic match ({d_abs / u_sys:.1f} sigma)"
                            elif d_abs > thresh_reject:
                                trans.accepted = 0
                                num_undecided = num_undecided - 1
                                if d_signed < 0:
                                    tag = "asym: Ic too weak vs Iobs"
                                elif d_signed > 0:
                                    tag = "asym: Ic too strong vs Iobs"
                                    if resonance:
                                        tag += ", resonance line"
                                else:
                                    tag = "Ic vs Iobs mismatch"
                                trans.notes2 = f"Step1 Rejected New: Strong Ic mismatch ({tag})"
                            # Else: Leave accepted as None for Step 2 (Doubtful/Borderline)
                        else:
                            # New with no I_calc: Only accept if Ritz is near-perfect and rank 1;
                            # rank here is residual rank for the undecided group, not original rank.
                            # The original sorting order is used here.
                            # Only the best of the undecided candidates is considered here.
                            if rank == 1:
                                if ritz_sigma <= 2.0:
                                    trans.accepted = 1
                                    num_undecided = num_undecided - 1
                                    trans.notes2 = "Step1 Accepted New: Excellent Ritz (no I_calc)"
                                else:
                                    trans.accepted = 0
                                    num_undecided = num_undecided - 1
                                    trans.notes2 = "Step1 Rejected New: Poor Ritz and no I_calc"
            step1_done = True
        elif not step2_done:
            num_decided_step2 = 0
            for trans in undecided:
                if blacklist and id(trans) in blacklist:
                    continue
                # 1. Initialize cumulative values from already accepted (accepted == 1) transitions
                accepted_already = _accepted_already()
                # Track cumulative intensity and weighted wavenumber sum
                i_cum = _i_cum(accepted_already)
                any_i_none = _any_i_none(accepted_already + [trans])
                # Calculate initial CoG
                if not any_i_none and i_cum > 0:
                    current_cog_wn = sum(t.calculated_wavenumber * t.calc_intensity for t in accepted_already) / i_cum
                else:
                    # Fallback to unweighted average if any I_calc is missing
                    current_cog_wn = np.mean(
                        [t.calculated_wavenumber for t in accepted_already]) if accepted_already else None
                current_cog_sigma = abs(
                    line.wavenumber - current_cog_wn) / line.wn_uncertainty if current_cog_wn else 999.0

                # Calculate total intensity and its uncertainty of what's already accepted + the current candidate
                current_group, i_total_theory, u_ln_theory_group = _current_group_stats(accepted_already, trans)

                decided, accepted_value, notes2 = _apply_relative_intensity_filter(
                    "Step2",
                    trans,
                    accepted_already,
                    i_cum,
                    i_total_theory,
                    u_ln_theory_group,
                    allow_insufficient_evidence_new_when_i_cum_zero=True
                )
                if decided:
                    trans.accepted = accepted_value
                    num_decided_step2 = num_decided_step2 + 1
                    num_undecided = num_undecided - 1
                    trans.notes2 = notes2
                    continue

                # --- CoG Sensitivity Test ---
                # See if accepting this candidate improves the match
                # test_group = accepted_already + [trans]
                # test_i_cum = sum(t.calc_intensity for t in test_group if t.calc_intensity is not None)
                test_any_none = any(t.calc_intensity is None for t in current_group)

                if not test_any_none and i_total_theory > 0:
                    # Calculated cog wn as weighted mean (calc intensities as weights)
                    test_cog_wn = sum(t.calculated_wavenumber * t.calc_intensity for t in current_group)/i_total_theory
                else:
                    # Calculated cog wn as unweighted mean
                    test_cog_wn = np.mean([t.calculated_wavenumber for t in current_group])

                test_cog_sigma = abs(line.wavenumber - test_cog_wn) / line.wn_uncertainty

                # --- Final Decision Logic for Step 2 ---
                improvement = current_cog_sigma - test_cog_sigma

                # Defer when nothing accepted yet (no last-rank forced legacy/new; no immediate CoG accept)
                if not accepted_already:
                    continue
                if improvement < -1.0:
                    trans.accepted = 0
                    num_decided_step2 = num_decided_step2 + 1
                    num_undecided = num_undecided - 1
                    trans.notes2 = f"Step2 Rejected: Worsened CoG sigma by {abs(improvement):.2f}"
                    continue

                # leave undecided for later Step2 substages; if still unresolved, Step3 will decide

            # --- Step2 pair stage: unordered pairs from non-rejected pool (accepted is None or 1) ---
            pair_progress, num_decided_step2, num_undecided = _try_pair_stage(include_broadening=False,
                                                                            n_decided_step2=num_decided_step2,
                                                                            n_undecided=num_undecided)
            if pair_progress:
                continue
            pair_progress, num_decided_step2, num_undecided = _try_pair_stage(include_broadening=True,
                                                                            n_decided_step2=num_decided_step2,
                                                                            n_undecided=num_undecided)
            if pair_progress:
                continue


            # --- Step2 triple stage: unordered triples from non-rejected pool (accepted is None or 1) ---
            triple_progress, num_decided_step2, num_undecided = _try_triple_stage(include_broadening=False,
                                                                            n_decided_step2=num_decided_step2,
                                                                            n_undecided=num_undecided)
            if triple_progress:
                continue
            triple_progress, num_decided_step2, num_undecided = _try_triple_stage(include_broadening=True,
                                                                            n_decided_step2=num_decided_step2,
                                                                            n_undecided=num_undecided)
            if triple_progress:
                continue

            # --- Step2 late arbitration: conflict flags in notes1 ---
            for trans in [t for t in line.assigned_transitions
                          if t.accepted is None and not (blacklist and id(t) in blacklist)]:
                notes1 = trans.notes1 or ""

                if trans.new == 1:
                    if 'R' in notes1:
                        trans.accepted = 0
                        num_decided_step2 += 1
                        num_undecided -= 1
                        trans.notes2 = (
                            "Step2 Rejected New: Conflict with legacy (notes1: R)"
                        )
                    elif 'F' in notes1:
                        trans.accepted = 0
                        num_decided_step2 += 1
                        num_undecided -= 1
                        trans.notes2 = (
                            "Step2 Rejected New: Conflicting/suspicious (notes1: F)"
                        )

            if num_decided_step2 == 0:
                step2_done = True
        else:
            # Step 3: Reject all undecided new candidates and accept all undecided old candidates
            for rank, trans in enumerate(undecided, start=1):
                if trans.new == 1:
                    trans.accepted = 0
                    trans.notes2 = "Step3 Rejected: new, no solid evidence for acceptance"
                else:
                    trans.accepted = 1
                    trans.notes2 = "Step3 Accepted: old, no solid evidence for rejection"
                if blacklist and id(trans) in blacklist:
                    trans.notes2 += "; decisions oscillate in iterations"
                num_undecided = num_undecided - 1

    num_decided_step1 = 0
    num_decided_step2 = 0
    num_decided_step3 = 0
    num_accepted = 0
    for t in line.assigned_transitions:
        if 'Step1' in t.notes2:
            num_decided_step1 = num_decided_step1 + 1
        elif 'Step2' in t.notes2:
            num_decided_step2 = num_decided_step2 + 1
        elif 'Step3' in t.notes2:
            num_decided_step3 = num_decided_step3 + 1
        if t.accepted == 1:
            num_accepted = num_accepted + 1
    return num_decided_step1, num_decided_step2, num_decided_step3, num_accepted

def weed_single_pass(observed_lines: list, blacklist: set = None):
    """Run one pass of the weeding algorithm over all observed lines."""
    num_decided_step1 = 0
    num_decided_step2 = 0
    num_decided_step3 = 0
    num_accepted = 0
    for line_no, line in enumerate(observed_lines, start=1):
        if not line.assigned_transitions or (len(line.assigned_transitions) == 1 and
                                             line.assigned_transitions[0].calculated_wavenumber == 0.0):
            continue
        nd1, nd2, nd3, na = weed_assignments_line(line, blacklist=blacklist)
        num_decided_step1 = num_decided_step1 + nd1
        num_decided_step2 = num_decided_step2 + nd2
        num_decided_step3 = num_decided_step3 + nd3
        num_accepted = num_accepted + na
    return num_decided_step1, num_decided_step2, num_decided_step3, num_accepted


def print_iteration_summary(iteration: int, stats: tuple):
    """Print weeding pass summary."""
    nd1, nd2, nd3, na = stats
    label = "Baseline" if iteration == 0 else f"Iteration {iteration}"
    print(f"  {label}: Step1={nd1}, Step2={nd2}, Step3={nd3}, Accepted={na}")


def weed_assignments(observed_lines: list, levels_dict: dict,
                     max_iterations: int = 59, alpha: float = 0.5, min_n: int = 5, verbose: bool = False):
    """Iterative weeding with per-level intensity adjustment factors.

    1. Run baseline weeding (no adjustments)
    2. Compute per-level factors from accepted transitions (with damping)
    3. Reset state, apply adjusted intensities, re-run weeding
    4. Repeat until convergence or max_iterations reached

    Args:
        observed_lines: list of SpectralLine objects.
        levels_dict: dictionary of energy levels, keyed by level ID.
        max_iterations: maximum number of weeding iterations to perform.
        alpha: damping coefficient for factor updates (0..1).
               1.0 = no damping, 0.5 = equal blend with previous iteration.
        min_n: minimum number of transitions for calculation of intensity adjustment factors
        verbose: True for printing diagnostics
    """
    if verbose: print("Step 5.3: Weeding assignments (iterative)...")
    build_level_transition_lists(levels_dict, observed_lines)

    # Build transition map once for oscillation detection
    trans_map = {
        id(t): t
        for line in observed_lines
        for t in line.assigned_transitions
        if t is not UNASSIGNED
    }

    blacklist = set()

    # Iteration 0: baseline (calc_intensity == orig_calc_intensity)
    stats = weed_single_pass(observed_lines, blacklist=blacklist)
    if verbose: print_iteration_summary(0, stats)
    prev_snapshot = snapshot_accepted(observed_lines)
    accepted_history = [prev_snapshot]
    prev_factors = None  # no damping on first factor computation

    for iteration in range(1, max_iterations + 1):
        diag = compute_intensity_factors(levels_dict, alpha=alpha,
                                         prev_factors=prev_factors, min_n=min_n)
        print_diagnostics = False
        if verbose and print_diagnostics:
            print_factor_diagnostics(iteration, diag)
            print_factor_stability(iteration, diag['level_details'], prev_factors)

        # Save applied factors for next iteration's damping
        prev_factors = {
            lid: (vals[1], vals[3])
            for lid, vals in diag['level_details'].items()
        }

        reset_weeding_state(observed_lines)
        apply_intensity_adjustments(observed_lines)

        stats = weed_single_pass(observed_lines, blacklist=blacklist)
        if verbose: print_iteration_summary(iteration, stats)

        curr_snapshot = snapshot_accepted(observed_lines)
        accepted_history.append(curr_snapshot)
        if verbose and print_diagnostics: detect_oscillations(accepted_history, trans_map)

        n_changed = count_changes(prev_snapshot, curr_snapshot)

        # Breakdown: accepted→rejected vs rejected→accepted
        n_newly_rejected = sum(
            1 for k in prev_snapshot
            if k in curr_snapshot and prev_snapshot[k] == 1 and curr_snapshot[k] != 1
        )
        n_newly_accepted = sum(
            1 for k in prev_snapshot
            if k in curr_snapshot and prev_snapshot[k] != 1 and curr_snapshot[k] == 1
        )
        if verbose:
            print(f"  Iteration {iteration}: {n_changed} changes "
                  f"(+{n_newly_accepted} accepted, -{n_newly_rejected} rejected)")

        # Blacklist transitions that flipped acceptance state
        if iteration > 2:
            for tid in prev_snapshot:
                if tid in curr_snapshot and prev_snapshot[tid] != curr_snapshot[tid]:
                    if tid not in blacklist:
                        blacklist.add(tid)

        if n_changed == 0:
            if verbose: print(f"  Converged after {iteration} iteration(s).")
            break

        prev_snapshot = curr_snapshot

    if blacklist and verbose:
        print(f"  Blacklisted {len(blacklist)} oscillating transition(s).")
    return stats

def assignment_cycle(observed_lines: list, all_possible: list, levels_dict: dict, verbose: bool = False) -> int:

    # Step 5.1: Match & grade
    transition_assignments = match_and_grade(observed_lines, all_possible, verbose=verbose)

    # Step 5.2: Resolve conflicts & output
    max_num_assignments, num_conflicts = resolve_conflicts(transition_assignments, verbose=verbose)
    if verbose:
        print(f'  Max number of conflicting assignments: {max_num_assignments}')
        print(f'  Total number of conflicting assignments: {num_conflicts}')

    # Step 5.3: Weed assignments (iterative with per-level intensity adjustments)
    stats = weed_assignments(observed_lines, levels_dict, max_iterations=100, alpha=0.5, min_n=5, verbose=verbose)
    nd1, nd2, nd3, na = stats
    if verbose: print(f"Weeding complete. Decisions made: {nd1} in Step1, {nd2} in Step2, {nd3} in Step3")
    print(f"  Total number of accepted assignments: {na}")
    return na


def calc_weights(lines: dict) -> dict:
    """Calculate weights for energy levels based on accepted transitions."""
    weights = {}
    for line in lines:
        all_accepted = [t for t in line.assigned_transitions if t.accepted == 1]
        if not all_accepted: continue
        any_none = any(t.calc_intensity is None for t in line.assigned_transitions)
        sum_i = float(len(all_accepted)) if any_none else sum(t.calc_intensity for t in all_accepted)
        for t in all_accepted:
            intens = 1.0 if any_none else t.calc_intensity
            w = 1.0 / (t.assigned_to.wn_uncertainty ** 2) * (intens / sum_i)
            weights[id(t)] = w
    return weights


def optimize_levels(levels_list: list, levels_history: list, weights: dict, verbose: bool = False) -> float:
    MAX_IT = 50  # Max. number of iterations
    TOL = 0.001  # Tolerance on change in level energies
    levels_history.append([(lev.level_id, lev.energy) for lev in levels_list])
    max_dif_it1 = 0.0
    for it in range (MAX_IT):
        max_change = 0.0
        for lev in levels_list:
            if lev.energy == 0.0: continue    # Skip the ground level
            # Compute weighted mean of energy differences for accepted transitions from/to this level
            e_sum = 0.0
            w_sum = 0.0

            for t in lev.from_transitions:
                if t.accepted != 1: continue
                # w = 1.0 / (t.assigned_to.wn_uncertainty ** 2)
                w = weights[id(t)] if id(t) in weights else 1.0 / (t.assigned_to.wn_uncertainty ** 2)
                e_sum += (t.lower_level.energy + t.assigned_to.wavenumber) * w
                w_sum += w
            for t in lev.to_transitions:
                if t.accepted != 1: continue
                # w = 1.0 / (t.assigned_to.wn_uncertainty ** 2)
                w = weights[id(t)] if id(t) in weights else 1.0 / (t.assigned_to.wn_uncertainty ** 2)
                e_sum += (t.upper_level.energy - t.assigned_to.wavenumber) * w
                w_sum += w

            if w_sum > 0:
                new_energy = e_sum / w_sum
                change = abs(new_energy - lev.energy)
                max_change = max(max_change, change)
                lev.energy = new_energy
        if verbose: print(f"Iteration {it+1}: max energy change = {max_change:.6f} cm^-1")
        if it == 0: max_dif_it1 = max_change
        if max_change < TOL:
            if verbose: print("Convergence achieved.")
            break
    return max_dif_it1


def clear_assignments(all_possible: list, lines: list):
    for t in all_possible:
        t.accepted = None
        t.notes1 = ""
        t.notes2 = ""
    for line in lines:
        line.assigned_transitions = []
        for tr in line.original_assignments:
            line.assigned_transitions.append(tr)


# ===========================================================================
# MAIN
# ===========================================================================
def main(max_cycles: int = 10):
    print("=" * 60)
    print("classify_lines.py — Spectral Line Classification for Pr III")
    print("=" * 60)

    # Step 1: Read energy levels
    levels_dict, levels_list = read_energy_levels()

    # Step 2: Read calculated transitions
    calc_trans_index = read_transitions(levels_dict)

    # Step 3: Read observed spectral lines
    observed_lines = read_observed_lines(levels_dict, calc_trans_index)

    # Step 4: Generate all possible transitions
    all_possible = generate_all_possible_transitions(levels_list, calc_trans_index)

    levels_history = []  # List of level snapshots per cycle
    i, na_prev, num_accepted, max_dif_it1 = 0, 0, 0, 0.0
    weights = {}
    # Step 5 - Main cycle: match & grade, resolve conflicts, weed assignments, optimize levels
    for i in range(max_cycles):
        print(f"Step 5 cycle {i+1}:")
        if i > 0: clear_assignments(all_possible, observed_lines)
        num_accepted = assignment_cycle(observed_lines, all_possible, levels_dict, verbose=False)
        weights = calc_weights(observed_lines)
        max_dif_it1 = optimize_levels(levels_list, levels_history, weights, verbose=False)
        print(f"  Accepted assignments: {num_accepted}, max energy change in optimization: {max_dif_it1:.6f} cm^-1")
        if num_accepted == na_prev and max_dif_it1 < 0.001: break
        na_prev = num_accepted
    if i >= max_cycles and (num_accepted != na_prev or max_dif_it1 >= 0.001):
        print(f"Warning: Step 5 iterations did not converge.")
    levels_history.append([(lev.level_id, lev.energy) for lev in levels_list])
    df = build_output(observed_lines, weights)
    write_output(df)

    print("=" * 60)
    print("Done.")
    print("=" * 60)


if __name__ == '__main__':
    main()
