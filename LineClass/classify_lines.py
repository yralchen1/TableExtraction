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

import argparse
import csv
import os
import math
from typing import List
import numpy as np
import bisect
import pandas as pd
import openpyxl
import config
import gA_imputation
from models import EnergyLevel, SpectralLine, Transition, UNASSIGNED

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Every file name, worksheet name and column name lives in lineclass_config.toml
# (see config.py).  The module-level names below are kept because the companion
# scripts (level_shifts.py, chance_mc.py, decoy_mc.py) refer to them as
# cl.ICALC_FILE, cl.WN_MIN and so on; apply_config() re-derives them whenever a
# different configuration file is given with --config.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

CFG = None          # the config.Config in force
LEVELS_FILE = ''    # workbook of the adopted energy levels
LINES_FILE = ''     # workbook of the observed lines
ICALC_FILE = ''     # workbook of the calculated transitions
OUTPUT_FILE = ''    # output workbook
OUTPUT_CSV = ''     # the same table as csv
LEVEL_OVERRIDES = ''  # csv of revised adopted energies, '' if none
WN_MIN = 0.0        # wavenumber range for possible transitions (cm^-1)
WN_MAX = 0.0
MISSING_POLICY = 'none'  # missing_gA.policy in force: 'none' or 'impute'
_IMPUTED = None          # cached imputation constants; see imputed_values()


def apply_config(cfg, policy: str = None) -> None:
    """Adopt `cfg` (a config.Config) as the configuration of this run.

    `policy` overrides missing_gA.policy for this run (the --missing-gA
    switch); None keeps the value written in the configuration file.
    """
    global CFG, LEVELS_FILE, LINES_FILE, ICALC_FILE, OUTPUT_FILE, OUTPUT_CSV
    global LEVEL_OVERRIDES, LINE_DECISIONS, NEW_LEVELS, ICALC_EXTRA
    global WN_MIN, WN_MAX, MISSING_POLICY, _IMPUTED
    CFG = cfg
    LEVELS_FILE = cfg.levels_file
    LINES_FILE = cfg.lines_file
    ICALC_FILE = cfg.icalc_file
    OUTPUT_FILE = cfg.output_file
    OUTPUT_CSV = cfg.output_csv
    LEVEL_OVERRIDES = cfg.level_overrides
    LINE_DECISIONS = cfg.line_decisions
    NEW_LEVELS = cfg.new_levels
    ICALC_EXTRA = cfg.icalc_extra
    WN_MIN = cfg.wn_min
    WN_MAX = cfg.wn_max
    MISSING_POLICY = policy if policy is not None else cfg.missing_gA['policy']
    if MISSING_POLICY not in ('none', 'impute'):
        raise config.ConfigError(
            f"missing_gA.policy = {MISSING_POLICY!r}; expected 'none' or 'impute'")
    _IMPUTED = None     # recalibrate on demand, against the new configuration


apply_config(config.load())


# ---------------------------------------------------------------------------
# The intensity given to a transition absent from the calculated-transition file
# ---------------------------------------------------------------------------
# Cowan's codes printed a transition only when gA (statistical weight times
# transition probability, s^-1) reached a cutoff, 1e3 s^-1 here.  A pair of
# levels missing from Icalc.xlsx is therefore not a transition of unknown
# strength: it is one known to be weaker than the cutoff.  Under
# missing_gA.policy = "impute" such a pair is given the intensity implied by a
# gA just below that cutoff, so that its predicted weakness counts against the
# identification instead of exempting it from the intensity tests; under
# policy = "none" it keeps no calculated intensity at all, as before.
def _imputation_calibration():
    """The constants used to stand for a censored transition, computed once.

    Returns (gA, u_ln, C, kT): the gA whose upper one-standard-deviation edge
    sits exactly on the printing cutoff, the uncertainty belonging to it (on
    the logarithmic scale), and the two constants of the intensity relation
    Icalc = C*gA*(rwn/1e8)*exp(-Eup/kT) that the file obeys (the predicted
    intensity is an energy flux, hence proportional to the wavenumber).
    gA_imputation.py
    derives all four and explains the choice.
    """
    global _IMPUTED
    if _IMPUTED is None:
        mg = CFG.missing_gA
        df = gA_imputation.load_icalc(CFG)
        gA, u_ln = gA_imputation.estimate_missing_gA(
            df, CFG.gA_cutoff, mg['u_ln_window'], mg['u_ln_estimator'],
            mg['self_consistent'], mg['fit_range_decades'])
        C, kT, _C_fit, _kT_fit, _agrees = gA_imputation.check_intensity_model(df, CFG)
        _IMPUTED = (gA, u_ln, C, kT)
        print(f"  Transitions absent from {os.path.basename(ICALC_FILE)} are imputed "
              f"with gA = {gA:.1f} s^-1 (cutoff {CFG.gA_cutoff:g} s^-1), "
              f"u_ln = {u_ln:.4f},")
        print(f"    through Icalc = {C:g} * gA * (rwn/1e8) * exp(-Eup/{kT:g}).")
    return _IMPUTED


def imputed_values(lower: EnergyLevel, upper: EnergyLevel):
    """(calc_intensity, u_calc) for a level pair absent from that file.

    Returns (None, None) under policy "none", so that the two policies run
    from the same code.  The intensity depends on the pair through the energy
    of the upper level and the wavenumber of the transition.

    The imputed gA is the one whose upper one-standard-deviation bound sits on
    the printing cutoff of the calculation; see the WARNING in
    `gA_imputation.estimate_missing_gA` for the justification that is correct
    and the one that is not.
    """
    if MISSING_POLICY != 'impute':
        return None, None
    gA, u_ln, C, kT = _imputation_calibration()
    rwn = upper.energy - lower.energy
    if rwn <= 0:
        return None, None
    I_calc = gA_imputation.impute_intensity(gA, upper.energy, rwn, C, kT)
    return float(I_calc), u_ln

# Atomic mass and plasma temperature for Doppler width calculation
ATOMIC_MASS = 140.90765  # u, standard mass unit
T = 1.6   # eV, plasma temperature

USE_INTENSITY_ADJUSTMENT = 0

# Weighting of the per-level intensity-factor fits: 'unweighted' or 'mandel_paule'.
# The residuals r = ln(I_obs/I_calc_orig) correlate strongly with u_calc (strong
# lines have small u_calc and systematically negative r), so 1/u_calc^2-weighted
# means are dragged to the strong-line offset and increase the overall intensity
# discrepancy instead of reducing it; unweighted means avoid this bias.
FACTOR_WEIGHTING = 'unweighted'


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


def column_index(ws, layout, source: str) -> dict:
    """0-based positions of the columns of `layout` in worksheet `ws`.

    The names are looked up in the header row (row 1) of the worksheet, so the
    readers do not depend on the order of the columns in the input file.
    Returns {logical name: index}; raises if a configured column is missing.
    """
    header = next(ws.iter_rows(min_row=1, max_row=1))
    return config.resolve_columns(header, layout.columns,
                                  f"{os.path.basename(source)}[{layout.sheet}]")


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
    ws = wb[CFG.levels.sheet]
    col = column_index(ws, CFG.levels, LEVELS_FILE)

    levels_dict = {}
    levels_list = []

    for row in ws.iter_rows(min_row=2):  # skip header
        j_str_val = row[col['J']].value        # adopted J
        # An empty E_ASD means the level is absent from the ASD, i.e. new in
        # Wyart's list.
        e_asd_val = row[col['E_ASD']].value
        energy_val = row[col['E']].value       # adopted energy
        parity_val = row[col['parity']].value
        level_id_val = row[col['id']].value

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
            J_val=j_val,
            is_new=1 if (e_asd_val is None or str(e_asd_val).strip() == '') else 0
        )
        levels_dict[level_id] = lev
        levels_list.append(lev)

    wb.close()
    print(f"  Read {len(levels_list)} energy levels.")
    if NEW_LEVELS:
        add_new_levels(levels_dict, levels_list, NEW_LEVELS)
    if LEVEL_OVERRIDES:
        apply_energy_overrides(levels_dict, LEVEL_OVERRIDES)
    return levels_dict, levels_list


def add_new_levels(levels_dict: dict, levels_list: list, path: str) -> int:
    """Append the levels of `path` to the level list read from the workbook.

    A level identified after the adopted level list was published exists in no
    input file of this pipeline: there is no row for it in the level workbook,
    which is an external published list and is never edited, and there are no
    calculated transitions for it in Icalc.xlsx.  Nothing could therefore ever
    propose an observed line for it, and a line_decisions.csv row accepting one
    of its lines would rule on a candidate that is never generated.  This file
    is how such a level enters.

    Columns: level_id, E (cm^-1), J (as written, "7/2" or "3"), parity ("e" or
    "o"); a comment column, and any other column, are ignored.  The level id
    is the next one free - the last six digits of the largest id in use, plus
    one - so that it can never collide with a published one.

    The level is marked is_new = 1, like a level of the published list that is
    absent from the ASD: it is a level whose reality the work is establishing,
    so the decoy calibration must treat it as one of the levels under test.
    is_added = 1 records that it came from here rather than from the workbook.
    """
    added = 0
    with open(path, newline='', encoding='utf-8-sig') as fh:
        rdr = csv.DictReader(fh)
        missing = [c for c in ('level_id', 'E', 'J', 'parity')
                   if c not in (rdr.fieldnames or [])]
        if missing:
            raise ValueError(f"{os.path.basename(path)}: missing column(s) "
                             f"{', '.join(missing)}")
        for rec in rdr:
            lid = to_str_id(rec['level_id'])
            if lid == '':
                continue
            if lid in levels_dict:
                raise ValueError(
                    f"{os.path.basename(path)}: level id {lid} is already in "
                    f"the level list; a new level must take the next id free")
            j_str = str(rec['J']).strip()
            parity = str(rec['parity']).strip()
            if parity not in ('e', 'o'):
                raise ValueError(f"{os.path.basename(path)}: level {lid} has "
                                 f"parity {parity!r}; expected 'e' or 'o'")
            lev = EnergyLevel(level_id=lid, energy=float(rec['E']),
                              parity=parity, J_str=j_str,
                              J_val=parse_J(j_str), is_new=1, is_added=1)
            levels_dict[lid] = lev
            levels_list.append(lev)
            added += 1
            print(f"    {lid}  E = {lev.energy:.4f} cm^-1, J = {j_str}, "
                  f"parity {parity}")
    if added:
        print(f"  Added {added} level(s) found since the level list was "
              f"published, from {os.path.basename(path)}.")
    return added


def read_energy_overrides(path: str) -> dict[str, float]:
    """The revised adopted energies of a csv, as level_id -> energy (cm^-1).

    The file has a header and the columns level_id and E_input; any further
    column (a comment recording where the value comes from) is ignored.  It is
    how a level that the identification work has MOVED enters the pipeline:
    the candidate transitions, the decoys planted around the level and the
    starting energy the optimizer is measured against must all use the new
    position, and the adopted-level workbook - an external, published list -
    is left untouched.
    """
    emap = {}
    with open(path, newline='', encoding='utf-8-sig') as fh:
        rdr = csv.DictReader(fh)
        missing = [c for c in ('level_id', 'E_input')
                   if c not in (rdr.fieldnames or [])]
        if missing:
            raise ValueError(f"{os.path.basename(path)}: missing column(s) "
                             f"{', '.join(missing)}")
        for rec in rdr:
            lid = to_str_id(rec['level_id'])
            if lid == '':
                continue
            emap[lid] = float(rec['E_input'])
    return emap


def apply_energy_overrides(levels_dict: dict, path: str) -> None:
    """Move the levels listed in `path` to their revised energies, in place.

    Raises if the file names a level that is not in the level list, so that a
    typo cannot pass silently as "no level was moved".
    """
    emap = read_energy_overrides(path)
    unknown = sorted(set(emap) - set(levels_dict))
    if unknown:
        raise ValueError(f"{os.path.basename(path)}: level id(s) not in the "
                         f"level list: {', '.join(unknown)}")
    print(f"  Energies overridden for {len(emap)} level(s) from "
          f"{os.path.basename(path)}:")
    for lid, e in sorted(emap.items()):
        lev = levels_dict[lid]
        print(f"    {lid}  {lev.energy:.4f} -> {e:.4f} "
              f"({e - lev.energy:+.4f} cm^-1)")
        lev.energy = e


# ---------------------------------------------------------------------------
# The decision ledger: the identifications ruled on by hand
# ---------------------------------------------------------------------------
# The classification is re-derived from scratch on every run, so without a
# record of them the verdicts reached by eye - in IDEN2, and in trial LOPT
# runs - are lost, and a candidate the analyst has already refused is proposed
# again next round.  The ledger is that record, and it is what lets a run be
# repeated without the accepted set moving: after all the automatic steps have
# run, every assignment the ledger names is set to the verdict written there,
# so only assignments never ruled on can change between rounds.  It is the
# exact analogue of the revised-energies file: a small csv under the analyst's
# hand, applied to an input that is never edited.
DECISIONS_ACCEPT = 'accept'
DECISIONS_REJECT = 'reject'
DECISIONS_WN_MATCH = 0.01   # cm^-1; how far a ledger wavenumber may miss its line


def read_line_decisions(path: str) -> dict:
    """The manual verdicts of a csv, as {(wn_obs, low_id, upp_id): (decision, reason)}.

    The file has a header and the columns wn_obs (the observed wavenumber of
    the line, cm^-1), low_id and upp_id (the two levels of the assignment) and
    decision ("accept" or "reject").  A `reason` column, and any other column
    (a date, say), are optional; `reason` is carried into the report so that
    the ground for the verdict travels with it.  The wavenumber is matched to
    the nearest observed line within DECISIONS_WN_MATCH, so it may be written
    to fewer decimals than the line list carries; the closest two observed
    lines of this spectrum are 0.10 cm^-1 apart, ten times that window, so the
    match cannot be ambiguous.
    """
    dmap = {}
    with open(path, newline='', encoding='utf-8-sig') as fh:
        rdr = csv.DictReader(fh)
        fields = rdr.fieldnames or []
        missing = [c for c in ('wn_obs', 'low_id', 'upp_id', 'decision')
                   if c not in fields]
        if missing:
            raise ValueError(f"{os.path.basename(path)}: missing column(s) "
                             f"{', '.join(missing)}")
        for n, rec in enumerate(rdr, start=2):
            if not (rec.get('wn_obs') or '').strip():
                continue
            decision = (rec['decision'] or '').strip().lower()
            if decision not in (DECISIONS_ACCEPT, DECISIONS_REJECT):
                raise ValueError(f"{os.path.basename(path)} line {n}: decision "
                                 f"= {rec['decision']!r}; expected "
                                 f"'{DECISIONS_ACCEPT}' or '{DECISIONS_REJECT}'")
            key = (float(rec['wn_obs']),
                   to_str_id(rec['low_id']), to_str_id(rec['upp_id']))
            if key in dmap and dmap[key][0] != decision:
                raise ValueError(f"{os.path.basename(path)} line {n}: "
                                 f"contradicts an earlier row on the same "
                                 f"assignment {key[1]}-{key[2]} at {key[0]}")
            dmap[key] = (decision, (rec.get('reason') or '').strip())
    return dmap


def attach_line_decisions(observed_lines: list, path: str) -> dict:
    """Hang the verdicts of `path` on the observed lines they rule on.

    Raises if a row names a wavenumber that is not an observed line, so that a
    mistyped wavenumber cannot pass silently as "no assignment was ruled on".
    A row whose two levels never come up as a candidate for that line is legal
    - the assignment may have moved out of matching range - and is reported at
    the end of the run by report_unapplied_decisions().

    Returns the verdicts as a list of (key, line, verdict) triples, key being
    (wn_obs, low_id, upp_id) as written in the file.
    """
    dmap = read_line_decisions(path)
    lines_sorted = sorted(observed_lines, key=lambda l: l.wavenumber)
    wns = [l.wavenumber for l in lines_sorted]

    def _line_at(wn):
        j = bisect.bisect_left(wns, wn)
        best, best_d = None, DECISIONS_WN_MATCH
        for k in (j - 1, j):
            if 0 <= k < len(wns) and abs(wns[k] - wn) <= best_d:
                best, best_d = lines_sorted[k], abs(wns[k] - wn)
        return best

    sites, unknown = [], []
    n_acc = n_rej = 0
    for key, verdict in dmap.items():
        wn, low, upp = key
        line = _line_at(wn)
        if line is None:
            unknown.append(wn)
            continue
        line.decisions[(low, upp)] = verdict
        sites.append((key, line, verdict))
        if verdict[0] == DECISIONS_ACCEPT:
            n_acc += 1
        else:
            n_rej += 1
    if unknown:
        raise ValueError(f"{os.path.basename(path)}: no observed line within "
                         f"{DECISIONS_WN_MATCH} cm^-1 of "
                         + ', '.join(f"{w:.4f}" for w in sorted(set(unknown)))
                         + " cm^-1")
    print(f"  Read {len(sites)} manual decision(s) from "
          f"{os.path.basename(path)}: {n_acc} accepted, {n_rej} rejected.")
    return sites


def manual_verdict(transition, line=None):
    """The ledger's verdict on `transition`, or None if it does not rule on it."""
    line = line if line is not None else transition.assigned_to
    if line is None or not line.decisions:
        return None
    if not (transition.lower_level and transition.upper_level):
        return None
    return line.decisions.get((transition.lower_level.level_id,
                               transition.upper_level.level_id))


def apply_line_decisions(line) -> None:
    """Overwrite this line's decisions with the verdicts of the ledger.

    Called at the very end of the weeding of a line, after Steps 1-3 and after
    the oscillation blacklist, so that the ledger is the last word on every
    assignment it names - which is what makes a rerun reproduce the analyst's
    accepted set instead of re-deriving a slightly different one.  The
    automatic verdict it replaces is kept in notes2, so the disagreements can
    be read off the output.
    """
    if not line.decisions:
        return
    for t in line.assigned_transitions:
        if t is UNASSIGNED:
            continue
        verdict = manual_verdict(t, line)
        if verdict is None:
            continue
        decision, reason = verdict
        want = 1 if decision == DECISIONS_ACCEPT else 0
        word = 'Accepted' if want else 'Rejected'
        t.manual = decision
        if t.accepted == want:
            t.notes2 = f"Manual {word} (agrees): {t.notes2}"
        else:
            t.notes2 = f"Manual {word} (overrides): {t.notes2}"
        t.accepted = want


def report_unapplied_decisions(sites: list) -> int:
    """Print the ledger rows that never met the assignment they rule on.

    Since an "accept" row creates its assignment (force_ledger_assignments),
    only "reject" rows should turn up here: a rejected pair can fall outside
    the matching tolerance when its levels move, or lose the transition to
    another line in the conflict resolution, and either way the verdict is
    moot because the assignment is not in the run to begin with.  A stale
    "accept" would mean the analyst's identification is NOT in this run's
    output; it should no longer be possible, and is still reported - and
    counted separately - in case it becomes so again.
    """
    missing = []
    for key, line, verdict in sorted(sites):
        pair = (key[1], key[2])
        if not any(t is not UNASSIGNED
                   and (t.lower_level.level_id, t.upper_level.level_id) == pair
                   for t in line.assigned_transitions):
            missing.append((key, verdict))
    if not missing:
        return 0
    lost = [k for k, v in missing if v[0] == DECISIONS_ACCEPT]
    print(f"  Warning: {len(missing)} manual decision(s) matched no candidate "
          f"of this run ({len(lost)} of them acceptances):")
    for (wn, low, upp), (decision, reason) in missing:
        print(f"    {wn:12.4f}  {low}-{upp}  {decision}"
              + (f"  ({reason})" if reason else ""))
    return len(missing)


def check_forced_decisions(sites: list, levels_dict: dict) -> None:
    """Abort if a ledger row cannot be carried out exactly as it is written.

    An "accept" row is an order: the assignment it names is put into the run
    whether or not the automatic matching ever proposes it (see
    force_ledger_assignments).  An order that cannot be obeyed must therefore
    stop the run rather than be quietly dropped, and two orders that
    contradict each other must stop it as well.  The cases:

      * a level id that is not in the level list - a typo, or a level that has
        been renumbered; this is an error for a "reject" row too, since such a
        row can never rule on anything;
      * an accepted pair that breaks the electric-dipole selection rules used
        to generate candidates (opposite parity, |J(upper) - J(lower)| <= 1,
        not both J = 0), or whose Ritz wavenumber E(upper) - E(lower) falls
        outside the range [WN_MIN, WN_MAX] the run covers - no such transition
        can exist, so the row is a mistake;
      * the same pair of levels accepted on two different observed lines - one
        transition can belong to only one line, so the two orders contradict
        each other.

    (Two rows that accept and reject the same pair on the SAME line are caught
    earlier, in read_line_decisions.  Accepting a pair on one line while
    rejecting it on another is not a contradiction: that is how an assignment
    is moved from one line to another, and both rows are obeyed.)

    All the faults are collected and reported together, so a file with several
    mistakes in it is fixed in one pass instead of one row per run.
    """
    faults = []
    accepted_at = {}
    for (wn, low, upp), line, (decision, _reason) in sorted(sites):
        lev_lo, lev_up = levels_dict.get(low), levels_dict.get(upp)
        for lid, lev in ((low, lev_lo), (upp, lev_up)):
            if lev is None:
                faults.append(f"{wn:12.4f}  {low}-{upp}  {decision}: "
                              f"level {lid} is not in the level list")
        if decision != DECISIONS_ACCEPT or lev_lo is None or lev_up is None:
            continue
        if lev_lo.parity == lev_up.parity:
            faults.append(f"{wn:12.4f}  {low}-{upp}  accept: both levels have "
                          f"parity {lev_lo.parity!r}; an electric-dipole "
                          f"transition connects opposite parities only")
        dj = abs(lev_up.J_val - lev_lo.J_val)
        if dj > 1.5:
            faults.append(f"{wn:12.4f}  {low}-{upp}  accept: "
                          f"|J(upper) - J(lower)| = {dj:g} > 1")
        elif lev_up.J_val == 0.0 and lev_lo.J_val == 0.0:
            faults.append(f"{wn:12.4f}  {low}-{upp}  accept: J = 0 on both "
                          f"levels; 0 - 0 is forbidden")
        ritz = lev_up.energy - lev_lo.energy
        if not (WN_MIN <= ritz <= WN_MAX):
            faults.append(f"{wn:12.4f}  {low}-{upp}  accept: the Ritz "
                          f"wavenumber {ritz:.4f} cm^-1 is outside the range "
                          f"[{WN_MIN:g}, {WN_MAX:g}] cm^-1 of this run")
        seen = accepted_at.get((low, upp))
        if seen is not None:
            faults.append(f"{wn:12.4f}  {low}-{upp}  accept: this pair is "
                          f"already accepted on the line at {seen:.4f} cm^-1; "
                          f"one transition cannot belong to two lines")
        else:
            accepted_at[(low, upp)] = wn
    if faults:
        raise ValueError("the decision ledger cannot be carried out:\n    "
                         + "\n    ".join(faults))


def force_ledger_assignments(obs_line, matches: list, proto: dict,
                             n_forced: list = None) -> list:
    """Put the assignments the ledger orders onto `obs_line`, and return them.

    The automatic matching proposes a pair only when the observed wavenumber
    sits within 5.5 combined standard deviations of the Ritz wavenumber
    E(upper) - E(lower).  An identification made by eye can be further out
    than that - the analyst has evidence the matching does not use (the
    appearance of the line on the plate, the branch structure of the level,
    the behaviour of a trial LOPT fit), and a hand-made identification is
    frequently what pulls a level back to where it belongs.  Such a pair was
    previously unreachable: with no candidate on the line there was nothing
    for the "accept" verdict to be applied to, and the row was reported as
    unapplied at the end of the run.

    So an "accept" row now CREATES the assignment when the matching does not
    propose it.  `proto` holds one Transition per ordered pair, taken from the
    generated candidate list, from which the calculated intensity and its
    uncertainty are copied; the new Transition is appended to the line and
    returned among the matches, so it is graded, takes part in the conflict
    resolution (where a manual acceptance wins the pair outright) and is set
    to accepted by apply_line_decisions at the end of the weeding, exactly
    like any other candidate the ledger rules on.
    """
    n_forced = [0] if n_forced is None else n_forced
    if not obs_line.decisions:
        return matches
    have = {(t.lower_level.level_id, t.upper_level.level_id)
            for t in obs_line.assigned_transitions if t is not UNASSIGNED}
    have |= {(t.lower_level.level_id, t.upper_level.level_id) for t in matches}
    for pair, (decision, _reason) in obs_line.decisions.items():
        if decision != DECISIONS_ACCEPT or pair in have:
            continue
        t = proto.get(pair)
        if t is None:                 # checked by check_forced_decisions
            continue
        forced = Transition(
            lower_level=t.lower_level,
            upper_level=t.upper_level,
            calc_intensity=t.calc_intensity,
            u_calc=t.u_calc,
            assigned_to=obs_line,
            orig_calc_intensity=t.orig_calc_intensity,
            orig_u_calc=t.orig_u_calc,
            is_imputed=t.is_imputed,
        )
        obs_line.assigned_transitions.append(forced)
        matches.append(forced)
        n_forced[0] += 1
    return matches


def forced_pair_prototypes(observed_lines: list,
                           all_possible_transitions: list) -> dict:
    """{(low_id, upp_id): Transition} for every pair the ledger accepts.

    One pass over the generated candidates, so that force_ledger_assignments
    can build its transitions without searching the list again for each line.
    """
    wanted = set()
    for lin in observed_lines:
        for pair, (decision, _reason) in lin.decisions.items():
            if decision == DECISIONS_ACCEPT:
                wanted.add(pair)
    if not wanted:
        return {}
    proto = {}
    for t in all_possible_transitions:
        pair = (t.lower_level.level_id, t.upper_level.level_id)
        if pair in wanted:
            proto[pair] = t
    return proto


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
    ws = wb[CFG.icalc.sheet]
    col = column_index(ws, CFG.icalc, ICALC_FILE)

    calc_trans_index = {}
    n_below = 0

    for row in ws.iter_rows(min_row=2):  # skip header
        id1_val = to_str_id(row[col['id1']].value)   # lower level
        id2_val = to_str_id(row[col['id2']].value)   # upper level

        if not id1_val or not id2_val:
            continue

        if id1_val not in levels_dict or id2_val not in levels_dict:
            continue

        # A few rows fall marginally below the printing cutoff of Cowan's codes
        # after the rescaling to Ritz wavenumbers; the configuration decides
        # whether they are read like any other row.
        if not CFG.allow_below_cutoff:
            gA_val = row[col['gA']].value
            try:
                if gA_val is not None and float(gA_val) < CFG.gA_cutoff:
                    n_below += 1
                    continue
            except (ValueError, TypeError):
                pass

        ua_pcnt_val = row[col['u_pct_gA']].value
        icalc_val = row[col['Icalc']].value

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
    if n_below:
        print(f"  Skipped {n_below} rows with gA < {CFG.gA_cutoff:g} s^-1 "
              f"(icalc.completeness.allow_below_cutoff = false).")
    if ICALC_EXTRA:
        read_extra_transitions(levels_dict, calc_trans_index, ICALC_EXTRA)
    return calc_trans_index


def read_extra_transitions(levels_dict: dict, calc_trans_index: dict,
                           path: str) -> int:
    """Merge the calculated transitions of a supplementary workbook.

    A level added by files.new_levels has no rows in Icalc.xlsx, so its
    transitions are supplied separately, in the same layout (the same
    worksheet name and the same column names).  Only the pairs involving such
    a level need be listed: a pair present in neither file is treated by the
    completeness rule exactly as before, as one whose gA falls below the
    printing cutoff of Cowan's codes.

    THE PREDICTED INTENSITY IS RECOMPUTED HERE, from gA, by the same relation
    the main file obeys,

        Icalc = C * gA * (rwn/1e8) * exp(-Eup/kT) ,

    with C and kT taken from [intensity_model] of the configuration.  It is
    not read from the Icalc column of the supplementary file.  The reason is
    that the main file's Icalc column is rewritten by tools/fit_boltzmann.py
    every time C and kT move, and this file is not: reading its column would
    put a handful of transitions on a stale intensity scale, which is exactly
    the sort of silent inconsistency the intensity tests would then blame on
    the identification.  Recomputing costs nothing and cannot drift.

    Returns the number of rows merged.
    """
    C = float(CFG.intensity_model['C'])
    kT = float(CFG.intensity_model['kT'])
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[CFG.icalc.sheet]
    col = column_index(ws, CFG.icalc, path)

    n, n_skipped, n_replaced = 0, 0, 0
    for row in ws.iter_rows(min_row=2):
        id1_val = to_str_id(row[col['id1']].value)
        id2_val = to_str_id(row[col['id2']].value)
        if not id1_val or not id2_val:
            continue
        if id1_val not in levels_dict or id2_val not in levels_dict:
            n_skipped += 1
            continue

        u_calc = None
        ua_pcnt_val = row[col['u_pct_gA']].value
        if ua_pcnt_val is not None:
            try:
                u_calc = math.log(float(ua_pcnt_val) / 100.0 + 1.0)
            except (ValueError, TypeError):
                pass

        calc_intensity = None
        try:
            gA = float(row[col['gA']].value)
            rwn = float(row[col['rwn']].value)
            Eup = float(row[col['Eup']].value)
            calc_intensity = C * gA * (rwn / 1e8) * math.exp(-Eup / kT)
        except (ValueError, TypeError):
            pass

        if (id1_val, id2_val) in calc_trans_index:
            n_replaced += 1
        calc_trans_index[(id1_val, id2_val)] = {
            'calc_intensity': calc_intensity,
            'u_calc': u_calc,
            'assigned_to': None
        }
        n += 1
    wb.close()
    print(f"  Read {n} further calculated transitions from "
          f"{os.path.basename(path)} (Icalc recomputed from gA with "
          f"C = {C:.6g}, kT = {kT:.6g} cm^-1).")
    if n_replaced:
        print(f"    {n_replaced} of them replaced a row of the main file.")
    if n_skipped:
        print(f"    {n_skipped} row(s) skipped: a level id not in the level list.")
    return n


# ===========================================================================
# Decoy (shadow) levels for in-situ false-positive calibration
# ===========================================================================
DECOY_PREFIX = 'D_'


def decoy_energy(energy: float, parity: str, delta: float) -> float:
    """Energy of the decoy copy of a level.

    The displacement sign alternates with parity so that EVERY transition
    involving a decoy is displaced from its true wavenumber: by delta for
    decoy-real pairs and by 2*delta for decoy-decoy pairs (a transition always
    connects levels of opposite parity, whose displacements have opposite
    signs and therefore never cancel in the energy difference).
    """
    return energy + (-delta if parity == 'o' else delta)


def add_decoy_levels(levels_dict: dict, levels_list: list,
                     calc_trans_index: dict, delta: float,
                     decoy_ids=None) -> int:
    """Insert a decoy (shadow) copy of every level to be validated.

    decoy_ids selects the levels to copy (an iterable of level ids). If it is
    None, the levels absent from ASD (is_new == 1) are copied. decoy_mc.py
    passes the ids of the new* levels — those supported by more new than old
    identifications in the baseline run — so that there is exactly one decoy
    per level whose reality is being tested.

    A decoy duplicates a Wyart new level in everything the pipeline can see —
    J, parity, connectivity to the level system, and calculated transition
    intensities — except that its energy is displaced by +-delta (sign
    alternating with parity, see decoy_energy). The decoys then compete with
    the real levels under completely realistic conditions: they are matched
    against the same observed lines, take part in conflict resolution and
    weeding, and are re-optimized together with the real levels, while the
    real levels stay anchored by the true identifications and most observed
    lines are already consumed by them.

    Because a decoy has, by construction, no true lines in the list, every
    line the pipeline accepts for it is a false positive obtained in situ.
    The per-decoy statistics (number of accepted lines n; wander |dE| of the
    optimized energy from the input value) calibrate the false-positive rate
    of new-level confirmations: see decoy_mc.py and level_shifts.py.

    |delta| must exceed the largest matching tolerance (5.5 * max u_obs) so
    that a decoy transition can never re-match the true line of its original;
    a few cm^-1 does not change the local level and line densities.

    Returns the number of decoy levels added.
    """
    if decoy_ids is None:
        new_levels = [lev for lev in levels_list if lev.is_new == 1 and lev.is_decoy == 0]
    else:
        decoy_ids = set(decoy_ids)
        missing = decoy_ids - set(levels_dict)
        if missing:
            raise ValueError(f"decoy_ids absent from the level list: {sorted(missing)[:5]} ...")
        new_levels = [lev for lev in levels_list
                      if lev.level_id in decoy_ids and lev.is_decoy == 0]
    for lev in new_levels:
        d = EnergyLevel(
            level_id=DECOY_PREFIX + lev.level_id,
            energy=decoy_energy(lev.energy, lev.parity, delta),
            parity=lev.parity,
            J_str=lev.J_str,
            J_val=lev.J_val,
            is_new=1,
            is_decoy=1,
        )
        levels_dict[d.level_id] = d
        levels_list.append(d)

    # Give decoy transitions the same calculated intensities as their
    # originals: for every Icalc entry involving a new level, register the
    # decoy-substituted id pairs ('assigned_to' stays None: decoys have no
    # legacy identifications). Both orientations are stored because the
    # energy ordering of a close pair can flip under the +-delta displacement.
    new_ids = {lev.level_id for lev in new_levels}
    extra = {}
    for (id1, id2), data in calc_trans_index.items():
        variants = set()
        if id1 in new_ids:
            variants.add((DECOY_PREFIX + id1, id2))
        if id2 in new_ids:
            variants.add((id1, DECOY_PREFIX + id2))
        if id1 in new_ids and id2 in new_ids:
            variants.add((DECOY_PREFIX + id1, DECOY_PREFIX + id2))
        for a, b in variants:
            entry = {'calc_intensity': data['calc_intensity'],
                     'u_calc': data['u_calc'],
                     'assigned_to': None}
            extra[(a, b)] = entry
            extra[(b, a)] = entry
    calc_trans_index.update(extra)
    print(f"  Added {len(new_levels)} decoy levels (displaced by +-{abs(delta):g} cm^-1) "
          f"and {len(extra)} decoy Icalc entries.")
    return len(new_levels)


# ===========================================================================
# STEP 3: Read observed spectral lines
# ===========================================================================
# noinspection PyTypeChecker,PyUnresolvedReferences
def read_observed_lines(levels_dict: dict, calc_trans_index: dict, wn_shift: float = 0.0,
                        drop_legacy: bool = False):
    """Read observed spectral lines from Pr3_lines.xlsx.

    Args:
        wn_shift: constant added to every observed wavenumber (cm^-1).
                  Used for Monte-Carlo estimation of chance coincidences:
                  a nonzero shift destroys all true line-transition
                  correspondences while preserving the statistical structure
                  of the line list (density, intensities, uncertainties).
        drop_legacy: if True, do not attach the legacy identifications from
                  the input file, so the pipeline classifies from scratch and
                  every candidate is treated as new. Required for
                  chance-coincidence runs: shifted legacy assignments would
                  otherwise be accepted by the Step-3 old-candidate default
                  despite grossly wrong residuals, and the level optimization
                  would then drag level energies to re-absorb the shift.

    Returns:
        observed_lines: list of SpectralLine objects
    """
    print("Step 3: Reading observed spectral lines...")
    if wn_shift != 0.0:
        print(f"  NOTE: all observed wavenumbers shifted by {wn_shift:+.3f} cm^-1 (chance-coincidence run)")
    wb = openpyxl.load_workbook(LINES_FILE, read_only=True, data_only=True)
    ws = wb[CFG.lines.sheet]
    col = column_index(ws, CFG.lines, LINES_FILE)

    observed_lines = []
    prev_line = None

    for row in ws.iter_rows(min_row=2):  # skip header
        wn_val = row[col['wn']].value            # wavenumber
        unc_val = row[col['u_wn']].value         # its uncertainty
        intens_val = row[col['intensity']].value
        char_val = row[col['character']].value   # line character

        if wn_val is None or unc_val is None or intens_val is None:
            continue

        try:
            wavenumber = float(wn_val) + wn_shift
            uncertainty = float(unc_val)
            intensity = float(intens_val)
        except (ValueError, TypeError):
            continue

        line_char = str(char_val).strip() if char_val is not None else ''

        id1_val = to_str_id(row[col['id1']].value)  # lower level_id
        id2_val = to_str_id(row[col['id2']].value)  # upper level_id

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

        if id1_val and id2_val and not drop_legacy:
            lower = levels_dict.get(id1_val)
            upper = levels_dict.get(id2_val)
            if lower is not None and upper is not None:
                calc_data = calc_trans_index.get((id1_val, id2_val))
                if calc_data is not None:
                    i_calc = calc_data['calc_intensity']
                    u_calc = calc_data['u_calc']
                    imputed = 0
                else:
                    # Absent from the calculated-transition file. A legacy
                    # identification must meet the same intensity test as a new
                    # one, so it is imputed by exactly the same rule.
                    i_calc, u_calc = imputed_values(lower, upper)
                    imputed = 1 if i_calc is not None else 0
                tr = Transition(
                    lower_level=lower,
                    upper_level=upper,
                    calc_intensity=i_calc,
                    u_calc=u_calc,
                    assigned_to=current_line,
                    orig_calc_intensity=i_calc,
                    orig_u_calc=u_calc,
                    is_imputed=imputed,
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
                # The pair is absent from the calculated-transition file, which
                # means gA below the printing cutoff. Imputing here rather than
                # in read_transitions gives decoy pairs the same treatment as
                # real ones, which keeps the false-positive calibration fair.
                i_calc, u_calc = imputed_values(lev_lo, lev_up)
                tr = Transition(
                    lower_level=lev_lo,
                    upper_level=lev_up,
                    calc_intensity=i_calc,
                    u_calc=u_calc,
                    orig_calc_intensity=i_calc,
                    orig_u_calc=u_calc,
                    is_imputed=1 if i_calc is not None else 0,
                )

            all_possible.append(tr)

    # Sort by calculated wavenumber
    all_possible.sort(key=lambda t: t.calculated_wavenumber)
    print(f"  Generated {len(all_possible)} possible transitions.")
    n_imputed = sum(1 for t in all_possible if t.is_imputed)
    if n_imputed:
        print(f"  {n_imputed} of them are absent from "
              f"{os.path.basename(ICALC_FILE)} and were given an imputed intensity.")
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


def match_and_grade(observed_lines: list, all_possible_transitions: list,
                    levels_dict: dict = None, verbose: bool = False):
    """For each observed line, find matching transitions and grade them.

    Modifies observed_lines in place (updates assigned_transitions and grades).
    Also tracks all new assignments for conflict resolution.

    The matching tolerance is 5.5 * combined_sigma, where combined_sigma
    folds in each candidate's level-energy uncertainties (EnergyLevel.u_energy,
    from compute_level_uncertainties(), as estimated from the previous
    cycle's accepted transitions) in quadrature with the line's own
    wavenumber uncertainty, not just the line's uncertainty alone. This lets
    a transition whose Ritz wavenumber depends on a not-yet-precisely-known
    new level be matched even when the observed line itself is very precise.

    Returns:
        transition_assignments: dict keyed by (lower_id, upper_id) ->
                                list of (Transition, SpectralLine, wn_diff_abs)
                                for conflict resolution
    """
    if verbose: print("Step 5.1: Matching observed lines to transitions & grading...")

    # Build a list of wavenumbers for binary search
    all_wn = [t.calculated_wavenumber for t in all_possible_transitions]

    # Widest possible per-level contribution to the tolerance, used only to
    # size the coarse bisect search window; the exact per-candidate tolerance
    # (below) still uses each candidate's own two levels.
    max_u_energy = max((lev.u_energy for lev in levels_dict.values()), default=0.0) if levels_dict else 0.0

    # Track all assignments: (lower_id, upper_id) -> [(transition, line, wn_diff)]
    transition_assignments = {}

    # The pairs the decision ledger orders into the run even when the matching
    # tolerance does not reach them (see force_ledger_assignments).
    proto = forced_pair_prototypes(observed_lines, all_possible_transitions)
    n_forced = [0]

    for line_idx, obs_line in enumerate(observed_lines):
        # if (line_idx + 1) % 1000 == 0:
        #     print(f"  Processing line {line_idx + 1}/{len(observed_lines)}...")

        search_tolerance = 5.5 * math.sqrt(obs_line.wn_uncertainty ** 2 + 2 * max_u_energy ** 2)
        wn_lo = obs_line.wavenumber - search_tolerance
        wn_hi = obs_line.wavenumber + search_tolerance

        idx_lo = bisect.bisect_left(all_wn, wn_lo)
        idx_hi = bisect.bisect_right(all_wn, wn_hi)

        matches = []
        for t in all_possible_transitions[idx_lo:idx_hi]:
            wn_diff_abs = abs(obs_line.wavenumber - t.calculated_wavenumber)
            u_lower = t.lower_level.u_energy if t.lower_level else 0.0
            u_upper = t.upper_level.u_energy if t.upper_level else 0.0
            tolerance = 5.5 * math.sqrt(obs_line.wn_uncertainty ** 2 + u_lower ** 2 + u_upper ** 2)
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
                        is_imputed=t.is_imputed,
                    )
                    obs_line.assigned_transitions.append(target_tr)
                
                matches.append(target_tr)

        # The analyst's own identifications, which the tolerance above may
        # not reach, are put on the line here so that everything downstream
        # treats them as ordinary candidates.
        if proto:
            matches = force_ledger_assignments(obs_line, matches, proto,
                                               n_forced)

        # Grade all candidates on this line: the fresh matches plus any seeded
        # original assignments whose Ritz wavenumber fell outside the matching
        # window this cycle. The latter still take part in the weeding (e.g.
        # the Step-2 blend logic), so they must not keep grade=None/new=None.
        matched_ids = {id(m) for m in matches}
        stale = [otr for otr in obs_line.assigned_transitions
                 if otr is not UNASSIGNED and id(otr) not in matched_ids]
        if matches or stale:
            assign_grades(obs_line, matches + stale)
        if matches:
            for m in matches:
                key = (m.lower_level.level_id, m.upper_level.level_id)
                if key not in transition_assignments:
                    transition_assignments[key] = []
                transition_assignments[key].append((m, obs_line, abs(obs_line.wavenumber - m.calculated_wavenumber)))

    if verbose:
        print(f"  Matching complete.")
        if n_forced[0]:
            print(f"  {n_forced[0]} assignment(s) created from the decision "
                  f"ledger that the matching tolerance does not reach.")
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
            # A pair the decision ledger accepts wins the transition outright,
            # one it rejects gives it up: a manual verdict could not be
            # honoured later by weed_assignments_line() if the transition were
            # taken away from the line here, before the weeding ever sees it.
            verdict = manual_verdict(t, lin)
            manual = 0 if verdict is None else (
                -1 if verdict[0] == DECISIONS_ACCEPT else 1)
            return manual, is_new, wn_diff_abs, tier, z, int_err

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
        n_accepted_line = sum(1 for t in obs_line.assigned_transitions
                              if t is not UNASSIGNED and t.accepted == 1)
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
            'imputed': 0,
            'intens_from_f': np.nan,
            'intens_to_f': np.nan,
            'dif_wn_O-C': np.nan,
            'grade': '',
            'notes1': '',
            'notes2': '',
            'manual': '',
            'new': '',
            'accepted': np.nan,
            'n_accepted': n_accepted_line,
            'low_E': np.nan,
            'upp_E': np.nan,
            'rwn': np.nan,
            'BF': np.nan
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
                # Output the branching fraction (BF) of this component, NOT the LOPT input
                # weight. The LOPT weight is BF**2 (the factor by which LOPT multiplies its
                # own 1/u_own**2); writing BF here avoids accidentally pasting BF**2 into
                # LOPT's "weight" input column. LOPT (centroid model) squares BF internally.
                bf = math.sqrt(weights[id(tr)] * (u_own**2)) if tr.accepted is not None and tr.accepted == 1 else 0.0

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
                    'imputed': tr.is_imputed,
                    'intens_from_f': tr.upper_level.intens_from_factor if tr.upper_level else np.nan,
                    'intens_to_f': tr.lower_level.intens_to_factor if tr.lower_level else np.nan,
                    'dif_wn_O-C': wn_diff,
                    'grade': tr.grade if tr.grade else '',
                    'notes1': tr.notes1,
                    'notes2': tr.notes2 if tr.notes2 else '',
                    'manual': tr.manual,
                    'new': tr.new if tr.new is not None else '',
                    'accepted': tr.accepted if tr.accepted is not None else np.nan,
                    'n_accepted': n_accepted_line,
                    'low_E': tr.lower_level.energy if tr.lower_level else np.nan,
                    'upp_E': tr.upper_level.energy if tr.upper_level else np.nan,
                    'rwn': tr.upper_level.energy-tr.lower_level.energy if tr.upper_level and tr.lower_level else np.nan,
                    'BF': bf
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
        'BF': '0.000'
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


def write_unstable_report(observed_lines: list, blacklist, path: str) -> int:
    """Write the assignments withdrawn for oscillating, one row each.

    These are the only rejections the pipeline makes on no evidence about the
    identification itself: the transition was withdrawn because its acceptance
    would not settle, which keeps the iteration from cycling but says nothing
    about whether the identification is right.  They therefore have to be
    visible, so they can be looked at in IDEN2 and settled for good in the
    decision ledger.  Returns the number of rows written.
    """
    info = getattr(blacklist, 'info', {})
    rows = []
    for line in observed_lines:
        for t in line.assigned_transitions:
            if t is UNASSIGNED:
                continue
            key = trans_key(t)
            if key not in blacklist:
                continue
            states, source = info.get(key, ((None, None, None), ''))
            rows.append({
                'wn_obs': line.wavenumber,
                'low_id': t.lower_level.level_id,
                'upp_id': t.upper_level.level_id,
                'obs_intens': line.intensity,
                'calc_intens': t.calc_intensity,
                'grade': t.grade or '',
                'new': t.new if t.new is not None else '',
                'oscillated_between': source,
                'states': '->'.join('-' if v is None else str(v) for v in states),
                'accepted': t.accepted if t.accepted is not None else '',
                'notes2': t.notes2 or '',
            })
    rows.sort(key=lambda r: -r['wn_obs'])
    print(f"  Writing {len(rows)} unstable candidate(s) to "
          f"{os.path.basename(path)}...")
    pd.DataFrame(rows, columns=['wn_obs', 'low_id', 'upp_id', 'obs_intens',
                                'calc_intens', 'grade', 'new',
                                'oscillated_between', 'states', 'accepted',
                                'notes2']).to_csv(path, index=False)
    return len(rows)


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
            t.manual = ""


def trans_key(t) -> tuple:
    """Stable identity of an assignment: (lower_id, upper_id, line).

    Transition objects are re-created by match_and_grade() in each outer
    Step-5 cycle, so id(t) cannot identify a transition across cycles.
    SpectralLine objects persist for the whole run, hence this key is
    stable and unique (a line holds at most one transition per level pair).
    """
    return t.lower_level.level_id, t.upper_level.level_id, id(t.assigned_to)


def snapshot_accepted(observed_lines: list) -> dict:
    """Return {trans_key(transition): accepted} for all non-UNASSIGNED transitions."""
    snap = {}
    for line in observed_lines:
        for t in line.assigned_transitions:
            if t is UNASSIGNED:
                continue
            snap[trans_key(t)] = t.accepted
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


def _level_energy_determinations(lev) -> tuple:
    """Return the energies that `lev`'s accepted transitions imply for it.

    Each accepted transition of the level is one independent determination of
    the level's energy: the partner level's energy plus (or minus) the
    observed wavenumber, with the uncertainty of that observed wavenumber
    combined in quadrature with the partner level's own energy uncertainty.

    Returns (keys, values, uncertainties): three parallel lists, one entry per
    accepted transition. The keys are trans_key() values, stable across
    Step-5 cycles, so that one determination can later be identified and
    dropped (see u_energy_excluding()).
    """
    keys, values, uncertainties = [], [], []
    for t in lev.from_transitions:  # lev is the upper level of t
        if t.accepted != 1:
            continue
        keys.append(trans_key(t))
        values.append(t.lower_level.energy + t.assigned_to.wavenumber)
        uncertainties.append(math.sqrt(t.assigned_to.wn_uncertainty ** 2 + t.lower_level.u_energy ** 2))
    for t in lev.to_transitions:  # lev is the lower level of t
        if t.accepted != 1:
            continue
        keys.append(trans_key(t))
        values.append(t.upper_level.energy - t.assigned_to.wavenumber)
        uncertainties.append(math.sqrt(t.assigned_to.wn_uncertainty ** 2 + t.upper_level.u_energy ** 2))
    return keys, values, uncertainties


def u_energy_excluding(level, transition) -> float:
    """u_energy of `level` recomputed without `transition`'s own determination.

    Every acceptance test compares a candidate's Ritz mismatch with the
    combined uncertainty sqrt(u_obs^2 + u_energy(lower)^2 + u_energy(upper)^2).
    Taken at face value, u_energy would make that test circular: a transition
    that disagrees with its level raises the scatter of the level's accepted
    determinations, the chi-squared inflation in _wm_red_chi() turns that
    scatter into a larger u_energy, and the larger u_energy softens the very
    test the transition has to pass. The worse the disagreement, the easier
    the transition gets in - which is how a level could end up holding three
    mutually incompatible lines.

    Leaving the candidate's own determination out removes the circularity at
    its source: each transition is judged against the level as defined by the
    *other* accepted transitions, which cannot be influenced by it. A level
    that rests on this transition alone is left with u_energy = 0, so the
    candidate is judged on the observed wavenumber uncertainty alone - the
    honest answer, since without this transition that level has no measured
    energy at all.

    The determinations left after the removal may still disagree among
    themselves, and their u_energy is still inflated accordingly; that is not
    circular, and it resolves itself over the Step-5 cycles: the worst
    offender fails its own test, is dropped, and the level re-forms from what
    remains.
    """
    if level is None:
        return 0.0
    contrib = level.u_contrib
    if not contrib or transition.assigned_to is None:
        return level.u_energy
    rest = [vu for k, vu in contrib.items() if k != trans_key(transition)]
    if len(rest) == len(contrib):     # this transition is not one of them
        return level.u_energy
    if not rest:
        return 0.0
    _, u_wm = _wm_red_chi([v for v, _ in rest], [u for _, u in rest])
    return u_wm if u_wm is not None else 0.0


def compute_level_uncertainties(levels_dict: dict, max_sweeps: int = 50, tol: float = 1e-4) -> None:
    """Estimate each level's u_energy from the scatter of energies implied by
    its currently accepted transitions (level.from_transitions/to_transitions,
    as left by the previous Step-5 cycle's weeding), and write it back onto
    each EnergyLevel. Must be called before matching/weeding starts for a
    cycle, since it consumes the *previous* cycle's accepted transitions.

    For a level L with an accepted transition to partner P (observed
    wavenumber uncertainty u_obs), that transition implies E_L = E_P +- wn_obs
    with combined uncertainty sqrt(u_obs^2 + u_energy(P)^2). u_energy(L) is
    then the chi-squared-inflated weighted-mean uncertainty of all such
    implied energies (_wm_red_chi, the same statistic already used for the
    per-level intensity factors).

    Because u_energy(P) is itself being solved for, this is done by
    fixed-point (Gauss-Seidel-style) iteration over all levels: levels
    connected only through chains of other new levels pick up their partners'
    still-converging uncertainty from the previous sweep, so uncertainty
    propagates outward from the well-determined levels without needing to
    identify a single "defining transition" for each level of a chain.

    Levels with no accepted transitions (including every level on the first
    cycle, before any assignment exists) keep u_energy = 0.0 (unknown); they
    pick up a real estimate once assignments accumulate in later cycles.
    """
    for lev in levels_dict.values():
        lev.u_energy = 0.0
        lev.u_contrib = {}

    for _sweep in range(max_sweeps):
        new_u = {}
        for lid, lev in levels_dict.items():
            _, values, uncertainties = _level_energy_determinations(lev)
            if values:
                _, u_wm = _wm_red_chi(values, uncertainties)
                new_u[lid] = u_wm if u_wm is not None else 0.0
            else:
                new_u[lid] = 0.0

        max_delta = max(abs(new_u[lid] - lev.u_energy) for lid, lev in levels_dict.items())
        for lid, lev in levels_dict.items():
            lev.u_energy = new_u[lid]
        if max_delta < tol:
            break

    # Keep the individual determinations behind each u_energy, so that the
    # acceptance tests can leave a transition's own determination out when
    # they judge that transition (u_energy_excluding()).
    for lev in levels_dict.values():
        keys, values, uncertainties = _level_energy_determinations(lev)
        lev.u_contrib = {k: (v, u) for k, v, u in zip(keys, values, uncertainties)}


def _compute_factor_for_transition_list(transitions: list, min_n: int=5, offset=None) -> tuple:
    """Compute a mean intensity adjustment factor from a list of transitions.

    Qualifying transitions: accepted == 1, orig_calc_intensity and orig_u_calc > 0,
    the parent line has exactly one accepted transition (unblended), and the
    calculated intensity is a real one - an imputed intensity is a bound read
    off the printing cutoff, not a measurement of how well the calculation
    reproduces this level, so it must not enter the fit.

    Args:
        transitions: candidate Transition list (a level's from- or to-list).
        min_n: minimum number of qualifying transitions.
        offset: optional callable t -> float, subtracted from each residual
                r = ln(I_obs/I_calc_orig) before averaging. The backfitting
                loop uses it to fit this factor on the partial residuals left
                by the complementary (from/to) factor.

    Returns (factor, u_factor, n_qualifying).
    If n_qualifying < min_n, returns (0.0, 0.0, n_qualifying).
    """
    CLAMP = 5.0

    qualifying = []
    for t in transitions:
        if t.accepted != 1:
            continue
        if t.is_imputed:
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
        r = math.log(t.assigned_to.intensity / t.orig_calc_intensity)
        if abs(r) > 3.0:
            continue  # Skip candidates with strongly deviating intensities
        qualifying.append((t, r))

    n = len(qualifying)
    if n < min_n or USE_INTENSITY_ADJUSTMENT != 1:
        return 0.0, 0.0, n

    # Residuals to average: r_i = ln(I_obs / I_calc_orig) - offset_i
    r_values = []
    r_uncertainties = []
    for t, r in qualifying:
        r_uncertainties += [t.orig_u_calc]
        r_values += [r - (offset(t) if offset else 0.0)]

    if FACTOR_WEIGHTING == 'mandel_paule':
        factor, u_factor = _wm_mandel_paule(r_values, r_uncertainties)
    else:
        factor = sum(r_values) / n
        u_factor = math.sqrt(sum((r - factor) ** 2 for r in r_values) / (n * (n - 1))) if n > 1 else 0.0
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

    The from- and to-factors are fitted by backfitting: the from-factors are
    computed on the residuals left by the current to-factors and vice versa,
    alternating until the factors stabilize. This solves the joint two-way
    model r = f_from(upper) + f_to(lower) + noise. Fitting each factor
    independently on the full residuals (as done previously) counts shared
    discrepancy structure twice, so the combined adjustment overshoots.

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

    # Backfitting: alternately refit the from-factors on the residuals left by
    # the to-factors and vice versa, until the factors stabilize.
    MAX_BACKFIT_SWEEPS = 50
    BACKFIT_TOL = 1e-4
    raw_from_map, raw_to_map = {}, {}
    u_from_map, u_to_map = {}, {}
    n_from_map, n_to_map = {}, {}
    for _sweep in range(MAX_BACKFIT_SWEEPS):
        max_delta = 0.0
        for lid, level in levels_dict.items():
            f, u_f, n = _compute_factor_for_transition_list(
                level.from_transitions, min_n,
                offset=lambda t: raw_to_map.get(t.lower_level.level_id, 0.0))
            max_delta = max(max_delta, abs(f - raw_from_map.get(lid, 0.0)))
            raw_from_map[lid], u_from_map[lid], n_from_map[lid] = f, u_f, n
        for lid, level in levels_dict.items():
            f, u_f, n = _compute_factor_for_transition_list(
                level.to_transitions, min_n,
                offset=lambda t: raw_from_map.get(t.upper_level.level_id, 0.0))
            max_delta = max(max_delta, abs(f - raw_to_map.get(lid, 0.0)))
            raw_to_map[lid], u_to_map[lid], n_to_map[lid] = f, u_f, n
        if max_delta < BACKFIT_TOL:
            break

    for lid, level in levels_dict.items():
        # Upper-level (emission) correction
        raw_from, u_from, n_from = raw_from_map[lid], u_from_map[lid], n_from_map[lid]
        if prev_factors is not None and lid in prev_factors:
            prev_from = prev_factors[lid][0]
            applied_from = alpha * raw_from + (1.0 - alpha) * prev_from
        else:
            applied_from = raw_from
        level.intens_from_factor = applied_from
        level.u_intens_from_factor = u_from
        if n_from >= min_n:
            n_from_computed += 1
            if applied_from != 0.0:
                from_factors.append(applied_from)
        else:
            n_from_skipped += 1

        # Lower-level correction
        raw_to, u_to, n_to = raw_to_map[lid], u_to_map[lid], n_to_map[lid]
        if prev_factors is not None and lid in prev_factors:
            prev_to = prev_factors[lid][1]
            applied_to = alpha * raw_to + (1.0 - alpha) * prev_to
        else:
            applied_to = raw_to
        level.intens_to_factor = applied_to
        level.u_intens_to_factor = u_to
        if n_to >= min_n:
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


class Blacklist(set):
    """The transitions withdrawn for oscillating, with the evidence for each.

    A plain set of trans_key() tuples - `key in blacklist` works as before -
    that also remembers, in `info`, why each key was put in: the three
    successive acceptance states that made it an oscillator and whether they
    were seen between weeding iterations or between outer cycles.  That is
    what write_unstable_report() prints, so that a withdrawal made on
    stability grounds alone (no evidence about the identification itself) is
    visible and can be overruled through the decision ledger.
    """

    def __init__(self):
        super().__init__()
        self.info = {}

    def record(self, key, states: tuple, source: str) -> None:
        self.add(key)
        self.info[key] = (states, source)


def _find_oscillations(snapshots: list, blacklist: set, source: str) -> int:
    """Blacklist the transitions whose acceptance goes A->B->A over the last
    three of `snapshots`, and return how many were newly blacklisted.

    A -> B -> A is the whole criterion.  A transition that changes its state
    once and keeps the new one has not oscillated: the surrounding loop is a
    damped fixed-point iteration, so a single flip is what convergence looks
    like from close up, and withdrawing on that alone would kill any candidate
    that happens to sit near a decision threshold while the intensity factors
    are still settling - a set that depends on every other transition in the
    run, and therefore changes from run to run.

    Absence from a snapshot means the transition was not a candidate at all
    that time - the levels had moved and it fell outside the matching
    tolerance - which is not a wavering decision but a missing one, so a key
    absent from any of the three is passed over.  Counting absence as a state
    (as this did before 2026-09-03) blacklisted 168 transitions in the third
    cycle of a normal run, while the level energies were still moving by of
    the order of 1 cm^-1 and candidates were coming and going for that reason
    alone.
    """
    if len(snapshots) < 3:
        return 0
    s1, s2, s3 = snapshots[-3], snapshots[-2], snapshots[-1]
    n_new = 0
    for key in set(s1) & set(s2) & set(s3):
        v1, v2, v3 = s1[key], s2[key], s3[key]
        if v1 == v3 and v2 != v1 and key not in blacklist:
            if isinstance(blacklist, Blacklist):
                blacklist.record(key, (v1, v2, v3), source)
            else:
                blacklist.add(key)
            n_new += 1
    return n_new


def detect_oscillations(accepted_history: list, trans_map: dict, verbose: bool=False):
    """Identify transitions that flip acceptance state A→B→A across last 3 snapshots.

    Args:
        accepted_history: list of snapshot dicts [{trans_key(t): accepted}, ...], oldest first.
        trans_map: {trans_key(t): Transition} for reverse lookup.
        verbose: if True, print diagnostics (default False).
    """
    if len(accepted_history) < 3:
        return
    s1, s2, s3 = accepted_history[-3], accepted_history[-2], accepted_history[-1]
    oscillating = []
    for tid in s1:
        if tid not in s2 and tid in s3:
            v1, v3 = s1[tid], s3[tid]
            if v1 != v3:
                oscillating.append((tid, v1, None, v3))
                continue
        if tid in s2 and tid not in s3:
            v1, v2 = s1[tid], s2[tid]
            if v1 != v2:
                oscillating.append((tid, v1, v2, None))
                continue
        if tid not in s2 or tid not in s3:
            continue
        v1, v2, v3 = s1[tid], s2[tid], s3[tid]
        if v1 == v3 and v1 != v2:
            oscillating.append((tid, v1, v2, v3))
    if not verbose:
        oscillating = []
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
        # Each level's uncertainty is taken as its *other* accepted
        # transitions leave it, so that a candidate cannot soften its own
        # test by disagreeing with its level (see u_energy_excluding()).
        u_lower = u_energy_excluding(candidate.lower_level, candidate)
        u_upper = u_energy_excluding(candidate.upper_level, candidate)
        combined_sigma = math.sqrt(line.wn_uncertainty ** 2 + u_lower ** 2 + u_upper ** 2)
        return abs(line.wavenumber - candidate.calculated_wavenumber) / combined_sigma

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

        The observed line profile is a convolution of the Doppler (thermal,
        FWHM-based)
        profile and the instrumental profile.  The effective resolution is
        sqrt(doppler_width² + wn_uncertainty²), so at low wavenumber where
        Doppler widths are small the measurement uncertainty dominates
        (recovering sigma-spread behavior), and at high wavenumber the
        Doppler width dominates.

        Returns dimensionless value:
            effective_spread / sqrt(doppler_width² + wn_uncertainty²)
        """

        if not group: return 999.0
        if temperature_eV <= 0: temperature_eV = 1.0
        if atomic_mass_u <= 0: atomic_mass_u = 100.0

        # --- constants ---
        eV_to_J = 1.602176634e-19
        u_to_kg = 1.66053906660e-27
        c = 2.99792458e8  # m/s (SI, consistent with sqrt(kT/m) below)

        # --- Doppler width (FWHM) in cm^-1 ---
        T_J = temperature_eV * eV_to_J
        m_kg = atomic_mass_u * u_to_kg

        # FWHM factor: sqrt(8*ln2 * kT / m) / c  (dimensionless; 8.22e-6 for
        # Pr at T = 1.6 eV)
        doppler_fwhm_factor = math.sqrt(8.0 * math.log(2.0) * T_J / m_kg) / c

        # Use representative wavenumber (CoG is best)
        cog = _cog(group)
        doppler_width = cog * doppler_fwhm_factor

        # --- combined resolution width ---
        effective_resolution = math.sqrt(doppler_width ** 2 + obs_wn_unc ** 2)

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
                # Do not let intensity alone accept a new line whose Ritz position is mediocre
                # (> 2.5 sigma); keep consistent with the single-candidate new policy.
                if _ritz_sigma(candidate) > 2.5:
                    return False, None, None
                # Compare Icalc with Iobs
                # noinspection PyUnresolvedReferences
                u_sys_tot = np.sqrt(candidate.u_calc ** 2 + u_obs_ln ** 2)
                # Z-test on the TOTAL sum vs the observation (asymmetric: too-weak vs too-strong)
                ln_ratio_tot = np.log(candidate.calc_intensity / line.intensity)
                d_ln_abs = abs(ln_ratio_tot)
                z_score = d_ln_abs / (S * u_sys_tot)
                if ln_ratio_tot > 0:
                    # The candidate's own Ic too strong vs Iobs: relax rejection threshold
                    z_reject = 1.0 * _fudge_factor_for_asym_intensity([candidate],
                                                                      candidate.calc_intensity, line.intensity, 1.0)
                    if z_score <= 1.0:
                        note = f"{step_label} Accepted New: Ic stat significant & compatible with Iobs (asym: Ic too strong)"
                        if _lower_id_resonance(candidate):
                            note += ", resonance lower"
                        return True, 1, note
                    if z_score > z_reject:
                        note = f"{step_label} Rejected new: Ic statistically incompatible with Iobs (asym: Ic too strong, z>{z_reject:.1f})"
                        return True, 0, note
                elif ln_ratio_tot < 0:
                    # The candidate's own Ic too weak vs Iobs: stricter rejection
                    if z_score <= 1.0:
                        return True, 1, f"{step_label} Accepted New: Ic stat significant & compatible with Iobs (asym: Ic too weak)"
                    if z_score > 1.25:
                        return True, 0, f"{step_label} Rejected New: Ic statistically incompatible with Iobs (asym: Ic too weak, z>1.25)"
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
                             if t.accepted != 0 and not (blacklist and trans_key(t) in blacklist)]
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
                   if t.accepted != 0 and not (blacklist and trans_key(t) in blacklist)]
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
                if blacklist and trans_key(trans) in blacklist:
                    continue
                # 1. Initialize cumulative values from already accepted (accepted == 1) transitions
                accepted_already = _accepted_already()

                # Track cumulative intensity and weighted wavenumber sum
                i_cum = _i_cum(accepted_already)

                # Ritz mismatch in sigma units
                ritz_sigma = _ritz_sigma(trans)
                u_sys = _u_sys(trans)

                if trans.new == 1:  # Test ritz_sigma first
                    if ritz_sigma > 3.0:
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
                    if ritz_sigma > 3.0:
                        trans.accepted = 0
                        num_undecided = num_undecided - 1
                        trans.notes2 = "Step1 Rejected New: Strong Ritz mismatch"
                    else:
                        if trans.calc_intensity is not None:
                            d_signed = np.log(trans.calc_intensity/line.intensity)
                            # Intensity can only rescue a new line while its Ritz position is
                            # still good (<= 2.5 sigma). Above that, no intensity-only accept:
                            # the line must be corroborated by a pair/triple (Step 2) or it is
                            # rejected in Step 3. This removes new IDs accepted on weak
                            # positional evidence (the moderate-tail NP-plot humps).
                            z_thresh = 1.0 if ritz_sigma <= 2.5 else 0.0
                            # OPEN QUESTION (2026-08-30): the intensity threshold here
                            # is 0.2, where every other call site passes 1.0. With 1.0
                            # the rejection threshold is relaxed only when the predicted
                            # intensity actually exceeds the observed one; with 0.2 it is
                            # relaxed whenever the prediction reaches a fifth of the
                            # observation, i.e. for most candidates. The commented-out
                            # line below shows the same 1/5 intent, so this looks
                            # deliberate, but the reason for it is no longer on record.
                            # Revisit if related anomalies show up in the output.
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
                if blacklist and trans_key(trans) in blacklist:
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
                          if t.accepted is None and not (blacklist and trans_key(t) in blacklist)]:
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
            # Step 3: Reject all undecided new candidates and accept all undecided old candidates.
            # Blacklisted (oscillating) transitions are always rejected: their instability is
            # evidence of an unreliable assignment, so they must stay dropped even if old.
            for rank, trans in enumerate(undecided, start=1):
                if blacklist and trans_key(trans) in blacklist:
                    trans.accepted = 0
                    trans.notes2 = "Step3 Rejected: decisions oscillate between weeding iterations or cycles"
                elif trans.new == 1:
                    trans.accepted = 0
                    trans.notes2 = "Step3 Rejected: new, no solid evidence for acceptance"
                else:
                    trans.accepted = 1
                    trans.notes2 = "Step3 Accepted: old, no solid evidence for rejection"
                num_undecided = num_undecided - 1

    # The analyst's own verdicts are the last word (see apply_line_decisions).
    apply_line_decisions(line)

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


def weed_assignments(observed_lines: list, levels_dict: dict, blacklist: set,
                     max_iterations: int = 59, alpha: float = 0.5, min_n: int = 10, verbose: bool = False):
    """Iterative weeding with per-level intensity adjustment factors.

    1. Run baseline weeding (no adjustments)
    2. Compute per-level factors from accepted transitions (with damping)
    3. Reset state, apply adjusted intensities, re-run weeding
    4. Repeat until convergence or max_iterations reached

    Args:
        observed_lines: list of SpectralLine objects.
        levels_dict: dictionary of energy levels, keyed by level ID.
        blacklist: set of trans_key() tuples of oscillating transitions.
                   Created once before the outer Step-5 cycle and shared
                   across cycles, so oscillations between cycles persist.
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
        trans_key(t): t
        for line in observed_lines
        for t in line.assigned_transitions
        if t is not UNASSIGNED
    }

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

        # Blacklist the transitions that genuinely oscillate, A->B->A, over the
        # last three iterations.  Before 2026-09-03 a single flip between two
        # consecutive iterations was enough, which withdrew candidates for
        # nothing worse than sitting near a threshold while the damped
        # intensity factors settled; see _find_oscillations().
        n_new_bl = _find_oscillations(accepted_history, blacklist, 'weeding iterations')
        if n_new_bl and verbose:
            print(f"  Blacklisted {n_new_bl} transition(s) oscillating "
                  f"between weeding iterations.")

        if n_changed == 0:
            if verbose: print(f"  Converged after {iteration} iteration(s).")
            break

        prev_snapshot = curr_snapshot

    if blacklist and verbose:
        print(f"  Blacklisted {len(blacklist)} oscillating transition(s) (cumulative).")
    return stats


def blacklist_cycle_oscillations(cycle_snapshots: list, blacklist: set) -> int:
    """Blacklist transitions whose end-of-cycle acceptance oscillates A→B→A
    across the last three outer Step-5 cycles.

    The weeding-internal detection in weed_assignments() cannot see these
    flips: it compares snapshots only between its own inner iterations, so a
    transition that is stable within every cycle but lands on the opposite
    decision in alternating cycles never produces an inner-iteration change.

    Absence from a snapshot (transition not matched in that cycle) counts as
    a distinct state, so present→absent→present is also treated as an
    oscillation: whenever the transition reappears it will stay rejected.

    Args:
        cycle_snapshots: list of end-of-cycle snapshot_accepted() dicts, oldest first.
        blacklist: shared set of trans_key() tuples; oscillators are added to it.

    Returns the number of newly blacklisted transitions.
    """
    return _find_oscillations(cycle_snapshots, blacklist, 'Step-5 cycles')


def assignment_cycle(observed_lines: list, all_possible: list, levels_dict: dict,
                     blacklist: set, verbose: bool = False) -> int:

    # Step 5.1: Match & grade
    transition_assignments = match_and_grade(observed_lines, all_possible, levels_dict, verbose=verbose)

    # Step 5.2: Resolve conflicts & output
    max_num_assignments, num_conflicts = resolve_conflicts(transition_assignments, verbose=verbose)
    if verbose:
        print(f'  Max number of conflicting assignments: {max_num_assignments}')
        print(f'  Total number of conflicting assignments: {num_conflicts}')

    # Step 5.3: Weed assignments (iterative with per-level intensity adjustments)
    stats = weed_assignments(observed_lines, levels_dict, blacklist, max_iterations=100, alpha=0.5, min_n=10, verbose=verbose)
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
            bf = intens / sum_i
            # Weight components of a blend by square of branching fraction to approximate the "centroid"
            # model in LOPT v. >= 5. Note that the "component" model used in earlier versions of LOPT
            # treats each blend component as an independent observation, which is not appropriate for blended lines.
            w = bf**2 / (t.assigned_to.wn_uncertainty ** 2)
            weights[id(t)] = w
    return weights


def optimize_levels(levels_list: list, levels_history: list, weights: dict, verbose: bool = False) -> float:
    """Optimize all level energies at once by weighted least squares.

    Every accepted transition is one observation equation

        E(upper level) - E(lower level) = wn(observed)

    weighted by w = weights[t] (the branching-fraction weight built by
    calc_weights) or, if the transition has no entry there, by 1/u^2, u
    being the uncertainty of the observed wavenumber.  The energies that
    minimize the weighted sum of squared residuals are the solution of the
    normal equations  (A^T W A) x = A^T W y , which are formed here by
    accumulating the 2x2 block each transition contributes and solved
    directly.  This is the exact minimum of the same model the previous
    iterative version approached; the earlier sweep stopped when a single
    pass moved no level by more than 0.001 cm^-1, which bounds the step,
    not the distance still left to the minimum, and so left up to
    ~0.05 cm^-1 of avoidable error in the strongly coupled levels.

    Levels sitting at 0.0 are held fixed (the ground level), as are levels
    that no accepted transition touches - nothing constrains those.  Should
    a group of levels be connected only to each other and not, through any
    chain of transitions, to a fixed level, its energies are determined only
    up to a common offset; one level of such a group is then held at its
    current energy so that the group keeps the position it has.

    Returns the largest energy change made, which the caller uses to decide
    whether the classification cycle has settled.
    """
    levels_history.append([(lev.level_id, lev.energy) for lev in levels_list])

    # ---- the observation equations -------------------------------------
    # A transition is reachable from both of its levels, so it is collected
    # through its identity to be counted once.
    trans = {}
    for lev in levels_list:
        for t in lev.from_transitions:
            if t.accepted == 1: trans[id(t)] = t
        for t in lev.to_transitions:
            if t.accepted == 1: trans[id(t)] = t
    trans = list(trans.values())
    if not trans: return 0.0

    def weight_of(t):
        if id(t) in weights: return weights[id(t)]
        return 1.0 / (t.assigned_to.wn_uncertainty ** 2)

    # ---- which levels are free to move ---------------------------------
    fixed = {id(lev) for lev in levels_list if lev.energy == 0.0}
    touched = set()
    for t in trans:
        touched.add(id(t.lower_level))
        touched.add(id(t.upper_level))

    # A group of levels tied only to each other floats as a whole; find the
    # groups with a simple union-find and hold one level of each floating
    # group in place.  ANCHOR stands for "any fixed level".
    ANCHOR = 0
    parent = {ANCHOR: ANCHOR}
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    def union(a, b):
        parent.setdefault(a, a); parent.setdefault(b, b)
        ra, rb = find(a), find(b)
        if ra == rb: return
        # ANCHOR stays the root of its group, so a group tied to a fixed
        # level is recognised as anchored however the union is written.
        if ra == ANCHOR: parent[rb] = ra
        else: parent[ra] = rb
    for t in trans:
        a, b = id(t.lower_level), id(t.upper_level)
        union(ANCHOR if a in fixed else a, ANCHOR if b in fixed else b)
    pinned = {}
    for lev in levels_list:
        if id(lev) not in touched or id(lev) in fixed: continue
        root = find(id(lev))
        if root == ANCHOR: continue
        if root not in pinned: pinned[root] = id(lev)
    fixed |= set(pinned.values())

    free = [lev for lev in levels_list
            if id(lev) in touched and id(lev) not in fixed]
    if not free: return 0.0
    col = {id(lev): j for j, lev in enumerate(free)}
    n = len(free)

    # ---- normal equations ----------------------------------------------
    energy = {id(lev): lev.energy for lev in levels_list}
    normal = np.zeros((n, n))
    rhs = np.zeros(n)
    for t in trans:
        w = weight_of(t)
        if w <= 0: continue
        y = t.assigned_to.wavenumber
        ju = col.get(id(t.upper_level))
        jl = col.get(id(t.lower_level))
        if ju is None: y -= energy[id(t.upper_level)]     # move to the rhs
        if jl is None: y += energy[id(t.lower_level)]
        if ju is not None:
            normal[ju, ju] += w
            rhs[ju] += w * y
        if jl is not None:
            normal[jl, jl] += w
            rhs[jl] -= w * y
        if ju is not None and jl is not None:
            normal[ju, jl] -= w
            normal[jl, ju] -= w

    try:
        solution = np.linalg.solve(normal, rhs)
    except np.linalg.LinAlgError:
        # Should not happen once every floating group is pinned; the
        # least-squares solution is used rather than failing outright.
        solution = np.linalg.lstsq(normal, rhs, rcond=None)[0]

    max_change = 0.0
    for lev, e_new in zip(free, solution):
        max_change = max(max_change, abs(e_new - lev.energy))
        lev.energy = float(e_new)
    if verbose:
        print(f"Level optimization: {n} levels from {len(trans)} transitions, "
              f"max energy change = {max_change:.6f} cm^-1")
    return max_change


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
def main(max_cycles: int = 20, wn_shift: float = 0.0, write_files: bool = True,
         decoy_shift: float = 0.0, decoy_ids=None) -> pd.DataFrame:
    """Run the full classification pipeline.

    Args:
        max_cycles: maximum number of Step-5 cycles.
        wn_shift: constant added to all observed wavenumbers (cm^-1);
                  nonzero values are used for chance-coincidence Monte Carlo.
        write_files: write line_classifications.xlsx/.csv (disable for
                     Monte-Carlo runs so the real output is not overwritten).
        decoy_shift: if nonzero, add a decoy copy of every level selected by
                  decoy_ids, displaced by +-decoy_shift (see
                  add_decoy_levels), for the in-situ false-positive
                  calibration; the observed lines and legacy identifications
                  are used unchanged.
        decoy_ids: ids of the levels to copy as decoys (defaults to the
                  levels absent from ASD).

    Returns the output DataFrame.
    """
    print("=" * 60)
    print("classify_lines.py — Spectral Line Classification for Pr III")
    if wn_shift != 0.0:
        print(f"CHANCE-COINCIDENCE RUN: wavenumber shift {wn_shift:+.3f} cm^-1")
    if decoy_shift != 0.0:
        print(f"DECOY RUN: shadow copies of the new levels displaced by +-{abs(decoy_shift):g} cm^-1")
    print("=" * 60)

    # Step 1: Read energy levels
    levels_dict, levels_list = read_energy_levels()

    # Step 2: Read calculated transitions
    calc_trans_index = read_transitions(levels_dict)

    # Optional: insert decoy levels for the in-situ false-positive calibration
    if decoy_shift != 0.0:
        add_decoy_levels(levels_dict, levels_list, calc_trans_index, decoy_shift,
                         decoy_ids=decoy_ids)

    # Step 3: Read observed spectral lines
    # Chance-coincidence runs classify from scratch (no legacy identifications)
    observed_lines = read_observed_lines(levels_dict, calc_trans_index, wn_shift=wn_shift,
                                         drop_legacy=(wn_shift != 0.0))

    # The manual verdicts.  A calibration run (a chance-coincidence shift or a
    # decoy planting) measures what the ALGORITHM does with an input it should
    # find nothing in, so the analyst's own decisions must not be laid over it;
    # they are read only for the real classification.
    decisions = {}
    if LINE_DECISIONS and wn_shift == 0.0 and decoy_shift == 0.0:
        decisions = attach_line_decisions(observed_lines, LINE_DECISIONS)

    # Step 4: Generate all possible transitions
    all_possible = generate_all_possible_transitions(levels_list, calc_trans_index)

    # A ledger row that cannot be obeyed stops the run here, before any of the
    # work is done, rather than being reported as unapplied at the end of it.
    if decisions:
        check_forced_decisions(decisions, levels_dict)

    levels_history = []  # List of level snapshots per cycle
    i, na_prev, num_accepted, max_lev_change = 0, 0, 0, 0.0
    weights = {}
    blacklist = Blacklist()  # oscillating transitions, keyed by trans_key(); persists across cycles
    cycle_snapshots = []  # end-of-cycle acceptance snapshots for cross-cycle oscillation detection
    # Step 5 - Main cycle: match & grade, resolve conflicts, weed assignments, optimize levels
    for i in range(max_cycles):
        print(f"Step 5 cycle {i+1}:")
        if i > 0: clear_assignments(all_possible, observed_lines)
        # Level.from_transitions/to_transitions still hold the previous cycle's
        # accepted transitions at this point (clear_assignments doesn't touch
        # them); on cycle 0 they're empty, so every level starts at u_energy=0.
        compute_level_uncertainties(levels_dict)
        num_accepted = assignment_cycle(observed_lines, all_possible, levels_dict, blacklist, verbose=True)
        cycle_snapshots.append(snapshot_accepted(observed_lines))
        n_new_bl = blacklist_cycle_oscillations(cycle_snapshots, blacklist)
        if n_new_bl:
            print(f"  Blacklisted {n_new_bl} transition(s) oscillating between cycles.")
        weights = calc_weights(observed_lines)
        max_lev_change = optimize_levels(levels_list, levels_history, weights, verbose=False)
        print(f"  Accepted assignments: {num_accepted}, max energy change in optimization: {max_lev_change:.6f} cm^-1")
        if num_accepted == na_prev and max_lev_change < 0.001: break
        na_prev = num_accepted
    if i >= max_cycles and (num_accepted != na_prev or max_lev_change >= 0.001):
        print(f"Warning: Step 5 iterations did not converge.")
    levels_history.append([(lev.level_id, lev.energy) for lev in levels_list])
    df = build_output(observed_lines, weights)
    if decisions:
        report_unapplied_decisions(decisions)
    if write_files:
        write_output(df)
        write_unstable_report(observed_lines, blacklist,
                              os.path.join(os.path.dirname(OUTPUT_CSV),
                                           'unstable_candidates.csv'))

    print("=" * 60)
    print("Done.")
    print("=" * 60)
    return df


def cli(argv=None) -> None:
    """Command-line entry point."""
    ap = argparse.ArgumentParser(
        description='Classify the observed lines of Pr III against all '
                    'transitions allowed between the known energy levels.')
    ap.add_argument('--config', default=config.DEFAULT_PATH, metavar='FILE',
                    help='configuration file (default: %(default)s)')
    ap.add_argument('--max-cycles', type=int, default=20, metavar='N',
                    help='maximum number of classify/optimize cycles '
                         '(default: %(default)s)')
    ap.add_argument('--missing-gA', choices=('none', 'impute'), default=None,
                    dest='missing_gA', metavar='{none,impute}',
                    help='how to treat a candidate transition absent from the '
                         'calculated-transition file: "none" leaves it without '
                         'a predicted intensity, "impute" gives it the intensity '
                         'implied by a gA just below the printing cutoff '
                         '(default: missing_gA.policy in the configuration)')
    args = ap.parse_args(argv)
    apply_config(config.load(args.config), policy=args.missing_gA)
    main(max_cycles=args.max_cycles)


if __name__ == '__main__':
    cli()
