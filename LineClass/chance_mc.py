"""Monte-Carlo estimate of chance-coincidence classifications.

Runs the full classification pipeline (classify_lines.main) with all observed
wavenumbers shifted by a constant. A uniform shift cannot be absorbed by the
level optimization (level energies enter only as differences), so every true
line-transition correspondence is destroyed, while the statistical structure
of the problem (line density, intensities, uncertainties, level list) is
preserved. Everything the pipeline ACCEPTS in a shifted run is therefore a
chance coincidence that survived matching, conflict resolution, weeding, and
level optimization — including the self-consolidation of spurious level
support through the optimization cycles.

All |shift| values exceed the largest matching half-window of the pipeline
(5.5 * max u_obs ≈ 6.3 cm^-1), so no true line-transition correspondence can
survive the shift even for the most uncertain lines; and the values are
incommensurate (not multiples of each other or of round numbers) to avoid
aliasing against any periodicity in the spectrum.

Shifted runs drop the legacy identifications from the input file (see
read_observed_lines drop_legacy): retained legacy assignments would be
accepted by the Step-3 old-candidate default despite grossly wrong residuals,
and the level optimization would drag level energies to re-absorb the shift.
Thus every candidate in a shifted run is new, and the accepted count is the
false-positive estimate for from-scratch identifications (the accepted-old
column is a consistency check and must be 0).

Outputs, aggregated over all shifts:
  - expected number of chance-accepted NEW classifications (the headline
    false-positive estimate);
  - per-level distribution of spurious support: how many levels acquire
    exactly n chance-accepted new lines. Comparing this with the same
    distribution from the real (unshifted) output calibrates the
    "questionable level" criterion.

Interpretation caveat: these chance numbers answer the question "how many
lines would the pipeline accept for a level that has NO relation to the
observed line list?". Wyart's new levels are not such levels — he derived
them from real, previously unclassified lines of this very list — so raw
acceptance counts must not be compared with them directly. The meaningful
test compares, per level, the number of supporting lines n together with the
level shift dE = E_final − E_input against the same two quantities of levels
that are false by construction; see level_shifts.py and decoy_mc.py.

Outputs written next to this script:
  chance_mc_summary.csv — per-shift acceptance counts and support histogram;
  chance_mc_levels.csv  — per level and shift: chance support and level shift dE;
  chance_mc_lines.csv   — accepted (chance) lines of every shifted run.

Usage:
    python chance_mc.py                      # full run, all shifts
    python chance_mc.py --smoke              # quick plumbing test: one shift, 2 cycles
    python chance_mc.py --missing-gA impute  # override missing_gA.policy
"""
import os
import sys
import time

import numpy as np
import openpyxl
import pandas as pd

import classify_lines as cl
import config

SHIFTS = [-11.3, -9.7, -8.3, -6.9, 6.9, 8.3, 9.7, 11.3]  # cm^-1, all > 5.5*max(u_obs)
MAX_SUPPORT_BIN = 8  # per-level support histogram: 1..MAX_SUPPORT_BIN+
SUMMARY_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'chance_mc_summary.csv')
LEVELS_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'chance_mc_levels.csv')
LINES_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'chance_mc_lines.csv')
# Accepted-line columns kept from each shifted run (for later residual analyses)
LINE_COLS = ['wn_obs', 'unc_wn_obs', 'obs_intens', 'low_id', 'upp_id',
             'dif_wn_O-C', 'grade', 'u_calc', 'calc_intens', 'imputed']


def accepted_stats(df: pd.DataFrame) -> dict:
    """Extract acceptance counts and per-level support from an output DataFrame."""
    new_flag = pd.to_numeric(df['new'], errors='coerce')
    assigned = df[df['upp_id'].astype(str).str.len() > 0]
    acc = df[df['accepted'] == 1]
    acc_new = acc[new_flag[acc.index] == 1]
    acc_old = acc[new_flag[acc.index] == 0]

    # Per-level support: each accepted new line supports both of its levels
    ids = pd.concat([acc_new['low_id'], acc_new['upp_id']]).astype(str)
    support = ids.value_counts()  # level_id -> number of accepted new lines
    hist = {}
    for n_lines, count in support.value_counts().items():
        hist[min(int(n_lines), MAX_SUPPORT_BIN)] = hist.get(min(int(n_lines), MAX_SUPPORT_BIN), 0) + int(count)

    return {
        'n_candidates': int(len(assigned)),
        'n_accepted': int(len(acc)),
        'n_accepted_new': int(len(acc_new)),
        'n_accepted_old': int(len(acc_old)),
        'support_hist': hist,   # {n_lines (capped): n_levels}
    }


def save_table(df: pd.DataFrame, csv_path: str, decimals: dict = None):
    """Write a table as CSV plus an .xlsx twin, in a form Excel reads cleanly.

    decimals maps column name -> number of decimal places to keep. Rounding
    serves two purposes: the files carry only physically meaningful digits,
    and Excel can read every number as a number. (pandas would otherwise
    print floats with up to 17 significant digits — the number of digits
    needed to reproduce a 64-bit float bit-for-bit — and Excel converts any
    number longer than 15 digits to text.) The .xlsx twin is written because
    Excel mangles text cells of a CSV, e.g. J = "3/2" is silently converted
    to a date; in .xlsx such cells stay text.
    """
    out = df.copy()
    if decimals:
        for col, nd in decimals.items():
            if col in out.columns:
                out[col] = out[col].round(nd)

    def _write(path, writer):
        # A file opened in Excel is locked on Windows; fall back to an
        # alternative name instead of aborting a long run.
        try:
            writer(path)
            return path
        except PermissionError:
            root, ext = os.path.splitext(path)
            alt = root + '_new' + ext
            writer(alt)
            print(f"  NOTE: {os.path.basename(path)} is locked (open in Excel?) "
                  f"- written to {os.path.basename(alt)} instead")
            return alt

    written = _write(csv_path, lambda p: out.to_csv(p, index=False))
    _write(os.path.splitext(csv_path)[0] + '.xlsx',
           lambda p: out.to_excel(p, index=False))
    print(f"Written: {written} (+ .xlsx twin)")


# Rounding used for the per-line tables saved by the MC drivers
LINE_DECIMALS = {'wn_obs': 4, 'unc_wn_obs': 4, 'obs_intens': 2,
                 'dif_wn_O-C': 4, 'u_calc': 4, 'calc_intens': 3}
LEVEL_DECIMALS = {'E_input': 4, 'E_final': 4, 'dE': 4}


def _text(val) -> str:
    """A cell value as a stripped string; an empty cell gives ''."""
    return str(val).strip() if val is not None else ''


def read_input_levels() -> pd.DataFrame:
    """Input (Wyart) level energies and old/new status from the levels file.

    A level is 'new' (found by Wyart, absent from the ASD compilation) iff its
    E_ASD cell is blank; E_adopt is the energy the pipeline starts from.
    note, J and parity are carried along for reporting.
    """
    wb = openpyxl.load_workbook(cl.LEVELS_FILE, read_only=True, data_only=True)
    ws = wb[cl.CFG.levels.sheet]
    col = cl.column_index(ws, cl.CFG.levels, cl.LEVELS_FILE)
    recs = []
    for row in ws.iter_rows(min_row=2):  # skip header
        level_id = cl.to_str_id(row[col['id']].value)
        energy_val = row[col['E']].value
        if not level_id or energy_val is None:
            continue
        e_asd = row[col['E_ASD']].value
        is_new = 1 if (e_asd is None or str(e_asd).strip() == '') else 0
        note = _text(row[col['note']].value)
        j_str = _text(row[col['J']].value)
        parity = _text(row[col['parity']].value)
        recs.append((level_id, float(energy_val), is_new, note, j_str, parity))
    wb.close()
    return pd.DataFrame(recs, columns=['level_id', 'E_input', 'is_new_level',
                                       'note', 'J', 'parity'])


def per_level_table(df: pd.DataFrame, levels: pd.DataFrame) -> pd.DataFrame:
    """Per-level accepted-line support and final energy from one run's output.

    The level-id columns of df must hold strings (when reading the output CSV,
    pass dtype={'low_id': str, 'upp_id': str}). Returns one row per level of
    `levels` with n_old / n_new (accepted lines touching the level), E_final
    (from the run output; E_input where the level never appears in a candidate
    row, i.e. it was not moved) and dE = E_final - E_input.
    """
    d = df.copy()
    for c in ('low_id', 'upp_id'):
        d[c] = d[c].fillna('').astype(str)
    new_flag = pd.to_numeric(d['new'], errors='coerce')
    counts = {}
    for label, flag in (('n_old', 0), ('n_new', 1)):
        acc = d[(d['accepted'] == 1) & (new_flag == flag)]
        counts[label] = pd.concat([acc['low_id'], acc['upp_id']]).value_counts()

    # Final optimized energy: identical on every candidate row of a level
    parts = [d.loc[d[c] != '', [c, e]].rename(columns={c: 'level_id', e: 'E_final'}).dropna()
             for c, e in (('low_id', 'low_E'), ('upp_id', 'upp_E'))]
    ef = pd.concat(parts, ignore_index=True)
    g = ef.groupby('level_id')['E_final']
    if len(ef):
        spread = float((g.max() - g.min()).max())
        if spread > 1e-6:
            raise ValueError(f"inconsistent final energies within a level "
                             f"(max spread {spread:g} cm^-1)")

    out = levels.copy()
    out['n_old'] = out['level_id'].map(counts['n_old']).fillna(0).astype(int)
    out['n_new'] = out['level_id'].map(counts['n_new']).fillna(0).astype(int)
    out['E_final'] = out['level_id'].map(g.mean()).fillna(out['E_input'])
    out['dE'] = out['E_final'] - out['E_input']
    return out


def hist_row(hist: dict) -> list:
    return [hist.get(n, 0) for n in range(1, MAX_SUPPORT_BIN + 1)]


def print_stats(label: str, s: dict):
    bins = ' '.join(f"{n}:{c:4d}" for n, c in zip(range(1, MAX_SUPPORT_BIN + 1),
                                                  hist_row(s['support_hist'])))
    print(f"{label}: candidates={s['n_candidates']}, accepted={s['n_accepted']} "
          f"(new={s['n_accepted_new']}, old={s['n_accepted_old']})")
    print(f"    levels by number of accepted-new lines  {bins}   "
          f"({MAX_SUPPORT_BIN} = {MAX_SUPPORT_BIN} or more)")


def apply_policy_option(argv=None) -> str:
    """Honour `--missing-gA {none,impute}` on the command line.

    A Monte-Carlo run must use the same treatment of the level pairs absent
    from the calculated-transition file as the real run it calibrates, so both
    Monte-Carlo scripts accept the same switch as classify_lines.py.  Without
    the switch the policy of the configuration file is used.  Returns the
    policy in force.
    """
    argv = sys.argv if argv is None else argv
    if '--missing-gA' in argv:
        i = argv.index('--missing-gA')
        if i + 1 >= len(argv) or argv[i + 1] not in ('none', 'impute'):
            raise SystemExit("--missing-gA takes 'none' or 'impute'")
        cl.apply_config(config.load(), policy=argv[i + 1])
    print(f"missing-gA policy: {cl.MISSING_POLICY}")
    return cl.MISSING_POLICY


def main():
    smoke = '--smoke' in sys.argv
    apply_policy_option()
    shifts = SHIFTS[:1] if smoke else SHIFTS
    max_cycles = 2 if smoke else 20
    if smoke:
        print("SMOKE TEST: one shift, max_cycles=2 — numbers are NOT meaningful.\n")

    levels = read_input_levels()
    rows = []
    level_rows = []
    line_rows = []
    for shift in shifts:
        t0 = time.time()
        df = cl.main(max_cycles=max_cycles, wn_shift=shift, write_files=False)
        s = accepted_stats(df)
        if s['n_accepted_old'] != 0:
            raise RuntimeError(f"shift {shift}: {s['n_accepted_old']} accepted OLD "
                               f"assignments in a shifted run (legacy ids not dropped?)")
        s['shift'] = shift
        s['minutes'] = (time.time() - t0) / 60.0
        rows.append(s)

        lvl = per_level_table(df, levels).rename(columns={'n_new': 'n_chance_lines'})
        lvl.insert(0, 'shift', shift)
        level_rows.append(lvl[['shift', 'level_id', 'is_new_level',
                               'n_chance_lines', 'E_input', 'E_final', 'dE']])
        acc_lines = df.loc[df['accepted'] == 1, LINE_COLS].copy()
        acc_lines.insert(0, 'shift', shift)
        line_rows.append(acc_lines)

        print()
        print_stats(f"SHIFT {shift:+.1f} cm^-1 ({s['minutes']:.1f} min)", s)
        print()

    # ------------------------------------------------------------------
    # Aggregate over shifts
    # ------------------------------------------------------------------
    print("=" * 60)
    print("CHANCE-COINCIDENCE SUMMARY")
    print("=" * 60)
    for s in rows:
        print_stats(f"shift {s['shift']:+.1f}", s)

    acc_new = np.array([s['n_accepted_new'] for s in rows], dtype=float)
    print(f"\nChance-accepted NEW classifications per run: "
          f"mean = {acc_new.mean():.1f}, std = {acc_new.std(ddof=1) if len(rows) > 1 else 0.0:.1f}, "
          f"range = {acc_new.min():.0f}..{acc_new.max():.0f}")

    mean_hist = np.array([hist_row(s['support_hist']) for s in rows], dtype=float).mean(axis=0)
    print("\nMean number of levels acquiring n chance-accepted new lines per run:")
    for n in range(1, MAX_SUPPORT_BIN + 1):
        suffix = '+' if n == MAX_SUPPORT_BIN else ' '
        # P(support >= n): mean count of levels with at least n chance lines
        ge = mean_hist[n - 1:].sum()
        print(f"  n = {n}{suffix}: {mean_hist[n - 1]:6.2f} levels   (>= {n}: {ge:6.2f})")

    # ------------------------------------------------------------------
    # Comparison with the real (unshifted) output, if available
    # ------------------------------------------------------------------
    if os.path.exists(cl.OUTPUT_CSV) and not smoke:
        real = pd.read_csv(cl.OUTPUT_CSV, dtype={'low_id': str, 'upp_id': str})
        s = accepted_stats(real)
        print("\nReal (unshifted) run, for comparison "
              f"[{cl.OUTPUT_CSV}]:")
        print_stats("real", s)
        # Support over ALL accepted lines (old + new), matching the
        # "questionable level" criterion based on identified transitions
        acc = real[real['accepted'] == 1]
        ids = pd.concat([acc['low_id'], acc['upp_id']]).astype(str)
        support_all = ids.value_counts().value_counts().sort_index()
        bins_all = ' '.join(
            f"{n}:{int(support_all[support_all.index >= n].sum()) if n == MAX_SUPPORT_BIN else int(support_all.get(n, 0)):4d}"
            for n in range(1, MAX_SUPPORT_BIN + 1))
        print(f"    levels by number of ALL accepted lines   {bins_all}   "
              f"({MAX_SUPPORT_BIN} = {MAX_SUPPORT_BIN} or more)")

    # ------------------------------------------------------------------
    # Per-shift CSV
    # ------------------------------------------------------------------
    out = pd.DataFrame([{
        'shift': s['shift'],
        'n_candidates': s['n_candidates'],
        'n_accepted': s['n_accepted'],
        'n_accepted_new': s['n_accepted_new'],
        'n_accepted_old': s['n_accepted_old'],
        **{f'levels_support_{n}{"plus" if n == MAX_SUPPORT_BIN else ""}': hist_row(s['support_hist'])[n - 1]
           for n in range(1, MAX_SUPPORT_BIN + 1)},
        'minutes': round(s['minutes'], 2),
    } for s in rows])
    print()
    save_table(out, SUMMARY_CSV, decimals={'minutes': 2})
    save_table(pd.concat(level_rows, ignore_index=True), LEVELS_CSV,
               decimals=LEVEL_DECIMALS)
    save_table(pd.concat(line_rows, ignore_index=True), LINES_CSV,
               decimals=LINE_DECIMALS)


if __name__ == '__main__':
    main()
