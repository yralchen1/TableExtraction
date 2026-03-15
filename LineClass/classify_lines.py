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
                    assigned_to=current_line
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
                    assigned_to=calc_data['assigned_to']
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
        tr.notes = ""  # notes now only holds F/R flags from conflict resolution


def match_and_grade(observed_lines: list, all_possible_transitions: list):
    """For each observed line, find matching transitions and grade them.

    Modifies observed_lines in place (updates assigned_transitions and grades).
    Also tracks all new assignments for conflict resolution.

    Returns:
        transition_assignments: dict keyed by (lower_id, upper_id) ->
                                list of (Transition, SpectralLine, wn_diff_abs)
                                for conflict resolution
    """
    print("Step 5: Matching observed lines to transitions & grading...")

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
                        assigned_to=obs_line
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

    print(f"  Matching complete.")
    return transition_assignments


# ===========================================================================
# STEP 6: Conflict resolution, "R" tagging, and output
# ===========================================================================
# def resolve_conflicts(observed_lines: list, transition_assignments: dict):
def resolve_conflicts(transition_assignments: dict):
    """Resolve conflicts and tag 'R' (Revised) notes."""
    print("Step 6: Resolving conflicts...")
    
    # 1. Selection for multiple lines pointing to one transition
    max_num_assignments = 0
    num_conflicts = 0
    for key, assignments in transition_assignments.items():
        if len(assignments) > max_num_assignments:
            max_num_assignments = len(assignments)
        if len(assignments) <= 1: continue

        num_conflicts += 1

        # if len(assignments) == 3:
        #     print("3-way conflict for transition", key, "candidates:",
        #           [(a[0].grade, a[1].wavenumber, a[2]) for a in assignments])

        # --- Shared helpers for all scoring approaches ---
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
                int_err = max(abs(math.log(Icalc / Iobs)) - u_Icalc, 0.0)
            else:
                int_err = 0.0

            return z, int_err

        # =====================================================================
        # APPROACH A: Lexicographic Sort (ACTIVE)
        # Sorts by (tier, z, int_err, is_new). No weight tuning needed.
        # Tier dominates, then z breaks ties, then intensity, then new-vs-old.
        # =====================================================================
        def sort_key_lexicographic(item):
            t, lin, wn_diff_abs = item
            z, int_err = _get_z_and_int_err(t, lin, wn_diff_abs)
            tier = 2 if z <= 2 else 3 if z <= 3 else 4 if z <= 4 else 5
            is_new = t.new if t.new is not None else 0
            return (tier, z, int_err, is_new)

        sort_key = sort_key_lexicographic

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

        assignments.sort(key=sort_key)
        winner_tr, winner_line, _ = assignments[0]
        winner_tr.notes += 'F'

        # Check if any original classification is being moved
        # (i.e., an original classification loses to a different line)
        original_moved = False
        for tr, line, _ in assignments[1:]:
            if tr.new == 0:
                original_moved = True

        if original_moved:
            winner_tr.notes += 'R'

        for tr, line, _ in assignments[1:]:
            if tr in line.assigned_transitions:
                line.assigned_transitions.remove(tr)
                if line.original_assignments and not line.assigned_transitions:
                    line.assigned_transitions.append(UNASSIGNED)

    return max_num_assignments, num_conflicts


def build_output(observed_lines: list) -> pd.DataFrame:
    """Build the output DataFrame from observed lines and their assignments.

    Numeric columns (calc_intens, CF, dif_wn_O-C) are stored as actual floats
    with NaN for missing values, so they appear as numbers in Excel/CSV.
    """
    print("  Building output...")
    import numpy as np
    output_rows = []

    for obs_line in observed_lines:
        if not obs_line.assigned_transitions:
            output_rows.append({
                'wn_obs': obs_line.wavenumber,
                'unc_wn_obs': obs_line.wn_uncertainty,
                'obs_intens': obs_line.intensity,
                'char': obs_line.line_character,
                'low_id': '',
                'upp_id': '',
                'calc_intens': np.nan,
                'u_calc': np.nan,
                'dif_wn_O-C': np.nan,
                'grade': '',
                'notes': '',
                'new': '',
                'low_E': np.nan,
                'upp_E': np.nan,
                'rwn': np.nan
            })
        else:
            for tr in obs_line.assigned_transitions:
                if tr is UNASSIGNED:
                    output_rows.append({
                        'wn_obs': obs_line.wavenumber,
                        'unc_wn_obs': obs_line.wn_uncertainty,
                        'obs_intens': obs_line.intensity,
                        'char': obs_line.line_character,
                        'low_id': '',
                        'upp_id': '',
                        'calc_intens': np.nan,
                        'u_calc': np.nan,
                        'dif_wn_O-C': np.nan,
                        'grade': '',
                        'notes': tr.notes,
                        'new': '',
                        'low_E': np.nan,
                        'upp_E': np.nan,
                        'rwn': np.nan
                    })
                    continue

                calc_wn = tr.calculated_wavenumber
                wn_diff = obs_line.wavenumber - calc_wn

                output_rows.append({
                    'wn_obs': obs_line.wavenumber,
                    'unc_wn_obs': obs_line.wn_uncertainty,
                    'obs_intens': obs_line.intensity,
                    'char': obs_line.line_character,
                    'low_id': tr.lower_level.level_id if tr.lower_level else '',
                    'upp_id': tr.upper_level.level_id if tr.upper_level else '',
                    'calc_intens': tr.calc_intensity,
                    'u_calc': tr.u_calc,
                    'dif_wn_O-C': wn_diff,
                    'grade': tr.grade if tr.grade else '',
                    'notes': tr.notes,
                    'new': tr.new if tr.new is not None else '',
                    'low_E': tr.lower_level.energy if tr.lower_level else np.nan,
                    'upp_E': tr.upper_level.energy if tr.upper_level else np.nan,
                    'rwn': tr.upper_level.energy-tr.lower_level.energy if tr.upper_level and tr.lower_level else np.nan
                })

    df = pd.DataFrame(output_rows)

    # Sort: decreasing wavenumber, then increasing grade for ties
    df['_sort_grade'] = df['grade'].apply(lambda g: g if g else 'ZZZZZ')
    df = df.sort_values(
        by=['wn_obs', '_sort_grade'],
        ascending=[False, True]
    ).drop(columns=['_sort_grade']).reset_index(drop=True)

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
        'u_calc': '0.000',
        'dif_wn_O-C': '0.000',
        'low_E': '0.00',
        'upp_E': '0.00',
        'rwn': '0.00'
    }
    # Find column indices from header row
    header = {cell.value: cell.column for cell in ws[1]}
    for col_name, fmt in col_formats.items():
        if col_name in header:
            col_idx = header[col_name]
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


WEEDER_INPUT_CSV = os.path.join(SCRIPT_DIR, 'weeder_input.csv')


def prepare_weeder_input(df: pd.DataFrame):
    """Generate weeder_input.csv from classification output.

    Filters to classified lines only and renames columns to match the
    expected input format for llm_weeder.py.
    """
    print("Step 7: Preparing weeder input...")
    import numpy as np

    # Only rows with a grade (i.e., classified lines)
    weeder = df[df['grade'].astype(str).str.strip() != ''].copy()

    # Rename columns to match llm_weeder expected format
    weeder = weeder.rename(columns={
        'wn_obs': 'obs_wn',
        'unc_wn_obs': 'sigma',
        'obs_intens': 'obs_I',
        'char': 'line_char',
        'calc_intens': 'calc_I',
    })

    # Add sequential id
    weeder.insert(0, 'id', range(1, len(weeder) + 1))

    # Select and order columns for the weeder
    cols = ['id', 'obs_wn', 'sigma', 'obs_I', 'line_char',
            'rwn', 'calc_I', 'u_calc', 'grade', 'new']
    # Only keep columns that exist
    cols = [c for c in cols if c in weeder.columns]
    weeder = weeder[cols]

    weeder.to_csv(WEEDER_INPUT_CSV, index=False)
    print(f"  Weeder input written: {len(weeder)} rows -> {WEEDER_INPUT_CSV}")


# ===========================================================================
# MAIN
# ===========================================================================
def main():
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

    # Step 5: Match & grade
    transition_assignments = match_and_grade(observed_lines, all_possible)

    # Step 6: Resolve conflicts & output
    max_num_assignments, num_conflicts = resolve_conflicts(transition_assignments)
    df = build_output(observed_lines)
    write_output(df)

    # Step 7: Prepare weeder input
    prepare_weeder_input(df)

    print("=" * 60)
    print("Done.")
    print(f'Max number of conflicting assignments: {max_num_assignments}')
    print(f'Total number of conflicting assignments: {num_conflicts}')
    print("=" * 60)


if __name__ == '__main__':
    main()
