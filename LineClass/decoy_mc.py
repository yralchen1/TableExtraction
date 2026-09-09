"""In-situ decoy (shadow-level) estimate of false new-level confirmations.

Purpose. The pipeline confirms an energy level when the weeding accepts
observed lines for the level's predicted transitions. For the new levels of
Wyart's list the question is how often such a confirmation could arise by
chance. The decoy experiment measures this directly, under fully realistic
run conditions.

Which levels are tested. The levels whose reality needs validation are the
"new*" levels: those supported in the baseline (no-decoy) run by more new
than old identifications (n_new > n_old). This set is read off the baseline
output line_classifications.csv; it comprises the levels absent from the ASD
compilation plus the levels whose energies, although present in ASD, rest
mainly on the new identifications.

Method. Each run repeats the REAL classification (unshifted lines, legacy
identifications in place) with one addition: for every new* level a decoy
copy is inserted (classify_lines.add_decoy_levels). The decoy inherits
everything the pipeline can see — J, parity, connectivity to the level
system, and the calculated intensities of its transitions — but its energy
is displaced by +-delta, with the sign alternating with parity so that no
transition involving a decoy can coincide with a true wavenumber. |delta|
exceeds the largest matching half-window (5.5 * max u_obs = 6.3 cm^-1), so a
decoy transition can never re-match the true line of its original; several
delta values sample different chance environments.

The decoys compete with the real levels in matching, conflict resolution,
weeding and level optimization, while the real levels stay anchored by the
true identifications and most observed lines are consumed by them — exactly
the environment a spurious level would face in the real run. Every line
accepted for a decoy is therefore a false positive obtained in situ, and,
with one decoy per tested level, the mean number of decoys passing a
criterion IS the expected number of false level confirmations under that
criterion (in the worst case that all tested levels were spurious).

Outputs (written next to this script):
  decoy_mc_levels.csv      — per decoy and run: accepted lines n and energy
                             wander dE = E_final − E_input;
  decoy_mc_real_levels.csv — per real level and run: support and dE, for
                             quantifying the perturbation of the real result
                             by the decoy competition;
  decoy_mc_lines.csv       — the accepted decoy lines themselves;
  decoy_mc_summary.csv     — per-run summary.

The (n, |dE|) criterion tables comparing the real new levels against this
calibration are printed by level_shifts.py.

Usage:
    python decoy_mc.py                      # full run, all deltas
    python decoy_mc.py --smoke              # plumbing test: one delta, 2 cycles
    python decoy_mc.py --missing-gA impute  # override missing_gA.policy
"""
import os
import sys
import time

import numpy as np
import pandas as pd

import classify_lines as cl
import chance_mc as mc
import output_files

DELTAS = [-11.3, -9.7, -8.3, -6.9, 6.9, 8.3, 9.7, 11.3]  # cm^-1
HERE = os.path.dirname(os.path.abspath(__file__))
LEVELS_CSV = os.path.join(HERE, 'decoy_mc_levels.csv')
REAL_LEVELS_CSV = os.path.join(HERE, 'decoy_mc_real_levels.csv')
LINES_CSV = os.path.join(HERE, 'decoy_mc_lines.csv')
SUMMARY_CSV = os.path.join(HERE, 'decoy_mc_summary.csv')


def decoy_level_frame(levels: pd.DataFrame, delta: float, ids: list) -> pd.DataFrame:
    """Input-side table of one run's decoy levels (id, E_input, J, parity)."""
    dec = levels[levels['level_id'].isin(set(ids))].copy()
    dec['E_input'] = [cl.decoy_energy(e, p, delta)
                      for e, p in zip(dec['E_input'], dec['parity'])]
    dec['level_id'] = cl.DECOY_PREFIX + dec['level_id']
    return dec


def main():
    smoke = '--smoke' in sys.argv
    mc.apply_policy_option()
    output_files.require_writable(
        [p for csv in (LEVELS_CSV, REAL_LEVELS_CSV, LINES_CSV, SUMMARY_CSV)
         for p in output_files.with_twin(csv)], 'table')
    deltas = DELTAS[:1] if smoke else DELTAS
    max_cycles = 2 if smoke else 20
    if smoke:
        print("SMOKE TEST: one delta, max_cycles=2 — numbers are NOT meaningful.\n")

    levels = mc.read_input_levels()

    # Baseline (no decoys): defines the tested (new*) levels and serves as
    # the reference for the perturbation check
    base = pd.read_csv(cl.OUTPUT_CSV, dtype={'low_id': str, 'upp_id': str})
    base_real = mc.per_level_table(base, levels)
    base_acc = int((base['accepted'] == 1).sum())
    star_ids = list(base_real.loc[base_real['n_new'] > base_real['n_old'], 'level_id'])
    n_new = len(star_ids)
    print(f"Tested (new*) levels: {n_new} of {len(levels)} "
          f"(supported by more new than old identifications in the baseline run)")

    sum_rows, dec_rows, real_rows, line_rows = [], [], [], []
    for delta in deltas:
        t0 = time.time()
        df = cl.main(max_cycles=max_cycles, wn_shift=0.0, write_files=False,
                     decoy_shift=delta, decoy_ids=star_ids)
        d = df.copy()
        for c in ('low_id', 'upp_id'):
            d[c] = d[c].fillna('').astype(str)
        is_decoy_row = (d['low_id'].str.startswith(cl.DECOY_PREFIX) |
                        d['upp_id'].str.startswith(cl.DECOY_PREFIX))
        acc = d['accepted'] == 1
        n_acc_real = int((acc & ~is_decoy_row).sum())
        n_acc_decoy = int((acc & is_decoy_row).sum())

        # Per-decoy support and energy wander
        dlev = mc.per_level_table(d, decoy_level_frame(levels, delta, star_ids))
        if int(dlev['n_old'].sum()) != 0:
            raise RuntimeError("decoy with legacy support — impossible")
        dlev = dlev.rename(columns={'n_new': 'n_decoy_lines'})
        dlev.insert(0, 'delta', delta)
        dec_rows.append(dlev[['delta', 'level_id', 'J', 'parity',
                              'n_decoy_lines', 'E_input', 'E_final', 'dE']])

        # Real levels in this run: perturbation caused by the decoy competition
        # (support counted over real-real rows only)
        rlev = mc.per_level_table(d[~is_decoy_row], levels)
        rlev.insert(0, 'delta', delta)
        real_rows.append(rlev[['delta', 'level_id', 'is_new_level',
                               'n_old', 'n_new', 'E_final', 'dE']])
        changed = (rlev['n_old'] + rlev['n_new']
                   != base_real['n_old'] + base_real['n_new'])

        lines = d.loc[acc & is_decoy_row, mc.LINE_COLS].copy()
        lines.insert(0, 'delta', delta)
        line_rows.append(lines)

        sup = dlev[dlev['n_decoy_lines'] > 0]
        s = {
            'delta': delta,
            'accepted_real': n_acc_real,
            'd_accepted_real_vs_baseline': n_acc_real - base_acc,
            'accepted_decoy_lines': n_acc_decoy,
            **{f'decoys_ge{n}': int((dlev['n_decoy_lines'] >= n).sum())
               for n in range(1, 9)},
            'median_absdE_supported': (round(float(sup['dE'].abs().median()), 4)
                                       if len(sup) else np.nan),
            'real_levels_support_changed': int(changed.sum()),
            'minutes': round((time.time() - t0) / 60.0, 2),
        }
        sum_rows.append(s)
        print(f"\nDELTA {delta:+.1f} cm^-1 ({s['minutes']:.1f} min): "
              f"real accepted={n_acc_real} ({s['d_accepted_real_vs_baseline']:+d} vs baseline), "
              f"accepted decoy lines={n_acc_decoy}, decoys with >=1 line: {len(sup)} of {n_new}")
        print(f"  real levels with changed support: {int(changed.sum())} of {len(levels)}\n")

    # ------------------------------------------------------------------
    # Aggregate over runs
    # ------------------------------------------------------------------
    print("=" * 60)
    print("IN-SITU DECOY SUMMARY")
    print("=" * 60)
    summ = pd.DataFrame(sum_rows)
    n_runs = len(summ)
    print(f"{n_runs} runs x {n_new} decoys = {n_runs * n_new} decoy trials")
    print(f"\nAccepted decoy lines per run: mean = {summ['accepted_decoy_lines'].mean():.1f}, "
          f"range = {summ['accepted_decoy_lines'].min()}..{summ['accepted_decoy_lines'].max()}")
    print("Mean decoys per run with at least n accepted lines "
          "(= expected false new-level confirmations at support >= n, any dE):")
    for n in range(1, 9):
        print(f"  n >= {n}: {summ[f'decoys_ge{n}'].mean():6.2f} of {n_new}")
    print(f"\nPerturbation of the real result by the decoy competition: "
          f"mean accepted real assignments {summ['accepted_real'].mean():.1f} "
          f"(baseline {base_acc}); real levels with changed support: "
          f"{summ['real_levels_support_changed'].mean():.1f} of {len(levels)} per run")

    print()
    mc.save_table(pd.concat(dec_rows, ignore_index=True), LEVELS_CSV,
                  decimals=mc.LEVEL_DECIMALS)
    mc.save_table(pd.concat(real_rows, ignore_index=True), REAL_LEVELS_CSV,
                  decimals=mc.LEVEL_DECIMALS)
    mc.save_table(pd.concat(line_rows, ignore_index=True), LINES_CSV,
                  decimals=mc.LINE_DECIMALS)
    mc.save_table(summ, SUMMARY_CSV, decimals={'median_absdE_supported': 4})
    print("\nRun level_shifts.py for the criterion tables and the per-level "
          "spurious probabilities.")


if __name__ == '__main__':
    main()
