"""
classify_lines.py
=================
Classify observed spectral lines of Pr III by matching them against all
possible transitions between known energy levels. Grade each match by
wavenumber agreement and intensity consistency, resolve conflicts, and
output a sorted classification table to Excel.

Steps:
  1. Read energy levels from Pr3_lev_Wyart_1999.xlsm
  2. Read observed spectral lines from Pr3_lines.xlsx
  3. Read DREAM calculated transitions from Pr3_tp_Dream.xlsm
  4. Generate all possible transitions satisfying selection rules
  5. Match observed lines to possible transitions & grade
  6. Resolve conflicts & write output
"""

import os
import bisect
import pandas as pd
import openpyxl
from models import EnergyLevel, SpectralLine, Transition

# ---------------------------------------------------------------------------
# Paths (relative to this script's directory)
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LEVELS_FILE = os.path.join(SCRIPT_DIR, '..', 'TableExtraction', 'Pr3_lev_Wyart_1999.xlsm')
LINES_FILE = os.path.join(SCRIPT_DIR, 'Pr3_lines.xlsx')
DREAM_FILE = os.path.join(SCRIPT_DIR, 'Pr3_tp_Dream.xlsm')
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


def grade_base(grade_str: str) -> str:
    """Strip suffix letters N, F, R from a grade for comparison purposes."""
    if grade_str is None:
        return ''
    base = grade_str
    for suffix in ('N', 'F', 'R'):
        base = base.replace(suffix, '')
    return base


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
def read_energy_levels() -> (dict, list):
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
# STEP 2: Read observed spectral lines
# ===========================================================================
def read_observed_lines(levels_dict: dict):
    """Read observed spectral lines from Pr3_lines.xlsx.

    Returns:
        observed_lines: list of SpectralLine objects
        all_assigned_transitions: list of Transition objects (Sugar's classifications)
        assigned_index: dict keyed by (lower_id, upper_id) -> Transition
    """
    print("Step 2: Reading observed spectral lines...")
    wb = openpyxl.load_workbook(LINES_FILE, read_only=True, data_only=True)
    ws = wb['Sheet1']

    observed_lines = []
    all_assigned_transitions = []
    assigned_index = {}  # (lower_id, upper_id) -> Transition

    prev_line = None

    for row in ws.iter_rows(min_row=2):  # skip header
        # own (col 1, idx 0), unc_own (col 2, idx 1), Icor (col 3, idx 2),
        # Ch. (col 4, idx 3), id1 (col 5, idx 4), id2 (col 6, idx 5)
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

        # Check if this is a continuation of the previous line
        # (same wavenumber = multiply-assigned line)
        if prev_line is not None and wavenumber == prev_line.wavenumber:
            # Same line, additional classification
            current_line = prev_line
        else:
            # New spectral line
            current_line = SpectralLine(
                wavenumber=wavenumber,
                wn_uncertainty=uncertainty,
                intensity=intensity,
                line_character=line_char
            )
            observed_lines.append(current_line)
            prev_line = current_line

        # If both level IDs are present, create a Transition (Sugar's classification)
        if id1_val and id2_val:
            lower = levels_dict.get(id1_val)
            upper = levels_dict.get(id2_val)
            if lower is not None and upper is not None:
                tr = Transition(
                    lower_level=lower,
                    upper_level=upper,
                    assigned_to=current_line
                )
                current_line.assigned_transitions.append(tr)
                all_assigned_transitions.append(tr)
                assigned_index[(id1_val, id2_val)] = tr

    wb.close()
    print(f"  Read {len(observed_lines)} observed lines with "
          f"{len(all_assigned_transitions)} existing classifications.")
    return observed_lines, all_assigned_transitions, assigned_index


# ===========================================================================
# STEP 3: Read DREAM calculated transitions
# ===========================================================================
def read_dream_transitions(levels_dict: dict, assigned_index: dict) -> dict:
    """Read DREAM calculated transitions from Pr3_tp_Dream.xlsm.

    Returns:
        calc_trans_index: dict keyed by (lower_id, upper_id) -> dict with
                          keys 'calc_intensity', 'CF', 'assigned_to'
    """
    print("Step 3: Reading DREAM calculated transitions...")
    wb = openpyxl.load_workbook(DREAM_FILE, read_only=True, data_only=True)
    ws = wb['trans']

    calc_trans_index = {}

    for row in ws.iter_rows(min_row=2):  # skip header
        # id1_W99 (col 15, idx 14), id2_W99 (col 16, idx 15),
        # CF (col 10, idx 9), Icalc (col 45, idx 44)
        id1_val = to_str_id(row[14].value)  # id1_W99
        id2_val = to_str_id(row[15].value)  # id2_W99

        if not id1_val or not id2_val:
            continue

        # Check that these levels exist in our levels list
        if id1_val not in levels_dict or id2_val not in levels_dict:
            continue

        cf_val = row[9].value     # CF (col 10)
        icalc_val = row[44].value  # Icalc (col 45)

        cf = None
        if cf_val is not None:
            try:
                cf = float(cf_val)
            except (ValueError, TypeError):
                pass

        calc_intensity = None
        if icalc_val is not None:
            try:
                calc_intensity = float(icalc_val)
            except (ValueError, TypeError):
                pass

        # Look up if this transition is already assigned to an observed line
        assigned_tr = assigned_index.get((id1_val, id2_val))
        assigned_to = None
        if assigned_tr is not None:
            assigned_to = assigned_tr.assigned_to
            for tr in assigned_to.assigned_transitions:
                if same_transitions(tr, assigned_tr):
                    # Update CF and calc_intensity in the assigned_to SpectralLine object,
                    # for the transition that matches this calculated transition
                    tr.calc_intensity = calc_intensity
                    tr.CF = cf
                    break

        calc_trans_index[(id1_val, id2_val)] = {
            'calc_intensity': calc_intensity,
            'CF': cf,
            'assigned_to': assigned_to
        }

    wb.close()
    print(f"  Read {len(calc_trans_index)} DREAM calculated transitions.")
    return calc_trans_index


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
                    CF=calc_data['CF'],
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
def compute_grade(wn_diff_abs: float, uncertainty: float,
                  calc_intensity: float, obs_intensity: float,
                  cf: float) -> str:
    """Compute the grade for a transition assignment.

    wn_diff_abs: |obs_wavenumber - calc_wavenumber|
    uncertainty: obs_line.wn_uncertainty
    calc_intensity: from DREAM (can be None)
    obs_intensity: observed intensity
    cf: cancellation factor (can be None)
    """
    # Determine tier
    if wn_diff_abs <= 1.0 * uncertainty:
        tier = 'A'
    elif wn_diff_abs <= 2.0 * uncertainty:
        tier = 'B'
    elif wn_diff_abs <= 3.0 * uncertainty:
        tier = 'C'
    else:
        tier = 'D'

    # Determine subgrade
    if cf is None or calc_intensity is None:
        subgrade = '6'
    else:
        intensity_ratio = calc_intensity / obs_intensity if obs_intensity != 0 else float('inf')
        if cf < 0.1:
            # Calculated intensity is almost meaningless
            if 0.1 <= intensity_ratio <= 10:
                subgrade = '4'
            else:
                subgrade = '5'
        else:
            # Compare observed and calculated intensities
            if 0.5 <= intensity_ratio <= 2:
                subgrade = '1'
            elif 0.1 <= intensity_ratio <= 10:
                subgrade = '2'
            else:
                subgrade = '3'

    return tier + subgrade


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
        if (line_idx + 1) % 1000 == 0:
            print(f"  Processing line {line_idx + 1}/{len(observed_lines)}...")

        tolerance = 4.5 * obs_line.wn_uncertainty
        wn_lo = obs_line.wavenumber - tolerance
        wn_hi = obs_line.wavenumber + tolerance

        # Binary search for transitions in range
        idx_lo = bisect.bisect_left(all_wn, wn_lo)
        idx_hi = bisect.bisect_right(all_wn, wn_hi)

        this_line_possible = all_possible_transitions[idx_lo:idx_hi]

        if not this_line_possible:
            # No possible transitions found
            if obs_line.assigned_transitions:
                for tr in obs_line.assigned_transitions:
                    tr.grade = 'X'
            continue

        for t in this_line_possible:
            wn_diff = obs_line.wavenumber - t.calculated_wavenumber
            wn_diff_abs = abs(wn_diff)

            # Only consider transitions within tolerance (should be guaranteed by search)
            if wn_diff_abs > tolerance:
                continue

            grade = compute_grade(
                wn_diff_abs, obs_line.wn_uncertainty,
                t.calc_intensity, obs_line.intensity, t.CF
            )

            # Check if this transition is already assigned to this line
            t_assigned = False
            for assigned_tr in obs_line.assigned_transitions:
                if same_transitions(t, assigned_tr):
                    t_assigned = True
                    assigned_tr.grade = grade
                    # Record for conflict resolution
                    key = (t.lower_level.level_id, t.upper_level.level_id)
                    if key not in transition_assignments:
                        transition_assignments[key] = []
                    transition_assignments[key].append(
                        (assigned_tr, obs_line, wn_diff_abs)
                    )
                    break

            if not t_assigned:
                # New assignment — create a copy of the transition for this line
                new_tr = Transition(
                    lower_level=t.lower_level,
                    upper_level=t.upper_level,
                    calc_intensity=t.calc_intensity,
                    CF=t.CF,
                    assigned_to=obs_line,
                    grade=grade + 'N'  # Tag as newly assigned
                )
                obs_line.assigned_transitions.append(new_tr)

                # Record for conflict resolution
                key = (t.lower_level.level_id, t.upper_level.level_id)
                if key not in transition_assignments:
                    transition_assignments[key] = []
                transition_assignments[key].append(
                    (new_tr, obs_line, wn_diff_abs)
                )

    print(f"  Matching complete.")
    return transition_assignments


# ===========================================================================
# STEP 6: Conflict resolution, "R" tagging, and output
# ===========================================================================
def resolve_conflicts(transition_assignments: dict):
    """Resolve conflicts where a transition is assigned to Multiple lines.

    Rules:
    - Compare base grades (ignore N, C, R suffixes); keep the better one.
    - If equal: prefer non-N (original) over N (new).
    - If still equal: prefer smaller wn residual.
    - Winner gets "C" appended (conflicting).
    - If an original (non-N) classification is moved from its original line
      to a different one, append "R" (revised) to the winner's grade.
    """
    print("Step 6: Resolving conflicts...")
    conflict_count = 0

    for key, assignments in transition_assignments.items():
        if len(assignments) <= 1:
            continue

        # Sort by: base grade (ascending), then prefer non-N, then by wn residual
        def sort_key(item):
            tr, line, wn_diff = item
            base = grade_base(tr.grade)
            is_new = 1 if (tr.grade and 'N' in tr.grade) else 0
            return (base, is_new, wn_diff)

        assignments.sort(key=sort_key)

        # The winner is the first one (best grade)
        winner_tr, winner_line, winner_diff = assignments[0]
        losers = assignments[1:]

        # Check if any original (non-N) classification is being moved
        # (i.e., an original classification loses to a different line)
        original_moved = False
        for tr, line, diff in losers:
            if tr.grade and 'N' not in tr.grade:
                # This is an original Sugar classification being displaced
                original_moved = True

        # Remove the losing transitions from their lines
        for tr, line, diff in losers:
            if tr in line.assigned_transitions:
                line.assigned_transitions.remove(tr)

        # Tag the winner
        base = grade_base(winner_tr.grade)
        suffixes = ''
        if winner_tr.grade and 'N' in winner_tr.grade:
            suffixes += 'N'
        suffixes += 'F'  # Conflicting assignment
        if original_moved:
            suffixes += 'R'  # Revised from original
        winner_tr.grade = base + suffixes

        conflict_count += 1

    print(f"  Resolved {conflict_count} conflicts.")


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
            # Unassigned line — transition-related fields are NaN/empty
            output_rows.append({
                'wn_obs': obs_line.wavenumber,
                'unc_wn_obs': obs_line.wn_uncertainty,
                'obs_intens': obs_line.intensity,
                'char': obs_line.line_character,
                'low_id': '',
                'upp_id': '',
                'calc_intens': np.nan,
                'CF': np.nan,
                'dif_wn_O-C': np.nan,
                'grade': ''
            })
        else:
            for tr in obs_line.assigned_transitions:
                calc_wn = tr.calculated_wavenumber
                wn_diff = obs_line.wavenumber - calc_wn

                # Round to the specified precision, keep as float
                calc_intens_val = round(tr.calc_intensity, 2) if tr.calc_intensity is not None else np.nan
                cf_val = round(tr.CF, 3) if tr.CF is not None else np.nan
                dif_val = round(wn_diff, 3)

                output_rows.append({
                    'wn_obs': obs_line.wavenumber,
                    'unc_wn_obs': obs_line.wn_uncertainty,
                    'obs_intens': obs_line.intensity,
                    'char': obs_line.line_character,
                    'low_id': tr.lower_level.level_id,
                    'upp_id': tr.upper_level.level_id,
                    'calc_intens': calc_intens_val,
                    'CF': cf_val,
                    'dif_wn_O-C': dif_val,
                    'grade': tr.grade if tr.grade else ''
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
        'calc_intens': '0.00',
        'CF': '0.000',
        'dif_wn_O-C': '0.000',
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


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    print("=" * 60)
    print("classify_lines.py — Spectral Line Classification for Pr III")
    print("=" * 60)

    # Step 1: Read energy levels
    levels_dict, levels_list = read_energy_levels()

    # Step 2: Read observed spectral lines
    observed_lines, all_assigned, assigned_index = read_observed_lines(levels_dict)

    # Step 3: Read DREAM calculated transitions
    calc_trans_index = read_dream_transitions(levels_dict, assigned_index)

    # Step 4: Generate all possible transitions
    all_possible = generate_all_possible_transitions(levels_list, calc_trans_index)

    # Step 5: Match & grade
    transition_assignments = match_and_grade(observed_lines, all_possible)

    # Step 6: Resolve conflicts & output
    resolve_conflicts(transition_assignments)
    df = build_output(observed_lines)
    write_output(df)

    print("=" * 60)
    print("Done.")
    print("=" * 60)


if __name__ == '__main__':
    main()
