#!/usr/bin/env python
"""Compare two runs of the classification pipeline.

Usage:
    python tools/compare_runs.py OLD_DIR NEW_DIR

Each directory must contain a `line_classifications.csv`, and optionally a
`level_shift_report.csv`; a path to a single .csv file may be given instead of
a directory, in which case only that pair of files is compared.

An "assignment" is one row of line_classifications.csv, identified by the
triple (wn_obs, low_id, upp_id): the observed wavenumber of the line and the
ids of the lower and upper energy levels of the candidate transition.  The
tool prints only counts and small tables - never a dump of the rows - so that
it stays usable inside a limited context.

Reported for the assignments:
  * rows present only in one of the runs (the candidate list itself changed);
  * acceptance changes among rows present in both, split by `new`
    (new = 1: identification proposed by this code; new = 0: legacy
    identification taken from the input line list);
  * how many of the changed rows are ones the missing-gA policy acts on: those
    whose level pair is absent from the calculated-transition file. They carry
    `imputed` = 1 when the run imputed an intensity for them, and an empty
    `orig_calc_intens` when it did not.

Reported for the level-shift report: the number of levels whose
`question_status` changed, the transition table of those statuses, and the
largest changes of `dE` (final minus input energy, cm^-1) and of `p_spur`
(estimated probability that the level is spurious).

Exit status is 0 when the two runs are identical in every compared column,
1 otherwise, so the tool can be used as a regression check.
"""
import argparse
import os
import sys

import pandas as pd

KEY = ['wn_obs', 'low_id', 'upp_id']
# Columns compared value-by-value for rows present in both runs.
NUM_COLS = ['calc_intens', 'orig_calc_intens', 'u_calc', 'dif_wn_O-C',
            'low_E', 'upp_E', 'rwn', 'BF', 'intens_from_f', 'intens_to_f']
TOL = 1e-9


def _indent(text: str, n: int) -> str:
    pad = ' ' * n
    return chr(10).join(pad + line for line in text.splitlines())


def _resolve(path: str, default_name: str) -> str:
    return path if os.path.isfile(path) and path.lower().endswith('.csv') \
        else os.path.join(path, default_name)


def _load_lines(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={'low_id': str, 'upp_id': str})
    df['wn_obs'] = df['wn_obs'].round(6)
    return df.set_index(KEY)


def compare_lines(old_path: str, new_path: str) -> bool:
    """Compare two line_classifications.csv. Returns True if identical."""
    a, b = _load_lines(old_path), _load_lines(new_path)
    print(f"line_classifications.csv:  {len(a)} rows -> {len(b)} rows")

    only_a = a.index.difference(b.index)
    only_b = b.index.difference(a.index)
    common = a.index.intersection(b.index)
    if len(only_a) or len(only_b):
        print(f"  candidate rows only in OLD: {len(only_a)}   only in NEW: {len(only_b)}")
    ca, cb = a.loc[common], b.loc[common]

    acc_a = ca['accepted'].fillna(-1).astype(int)
    acc_b = cb['accepted'].fillna(-1).astype(int)
    print(f"  accepted, counted on the {len(common)} rows common to both runs: "
          f"{int((acc_a == 1).sum())} -> {int((acc_b == 1).sum())}")

    changed = acc_a != acc_b
    identical = not (len(only_a) or len(only_b) or changed.any())
    if changed.any():
        # -1 stands for "no decision recorded" (empty cell in the csv).
        tab = pd.crosstab(acc_a[changed], acc_b[changed],
                          rownames=['old accepted'], colnames=['new accepted'])
        print(f"  acceptance changed on {int(changed.sum())} rows:")
        print(_indent(tab.to_string(), 4))
        lost = changed & (acc_a == 1) & (acc_b != 1)
        gained = changed & (acc_a != 1) & (acc_b == 1)
        for label, mask in (('lost', lost), ('gained', gained)):
            if not mask.any():
                continue
            m = ca[mask]
            n_new = int((m['new'] == 1).sum())
            n_leg = int((m['new'] == 0).sum())
            # A row is one of the censored pairs if either run treated it as
            # such: `imputed` = 1 in the run that imputed, an empty calculated
            # intensity in the run that did not.
            censored = m['orig_calc_intens'].isna()
            for frame in (ca[mask], cb[mask]):
                if 'imputed' in frame.columns:
                    censored = censored | (frame['imputed'].fillna(0) == 1).to_numpy()
            n_noi = int(censored.sum())
            print(f"    {label}: {int(mask.sum())}  (new={n_new}, legacy={n_leg},"
                  f" absent from the calculated-transition file={n_noi})")

    for col in NUM_COLS:
        if col not in ca.columns or col not in cb.columns:
            continue
        x, y = ca[col], cb[col]
        both_nan = x.isna() & y.isna()
        diff = (~both_nan) & ((x.isna() != y.isna()) | ((x - y).abs() > TOL))
        if diff.any():
            identical = False
            worst = (x - y).abs().max()
            print(f"  {col}: {int(diff.sum())} rows differ (max |diff| = {worst:.6g})")
    return identical


def compare_shifts(old_path: str, new_path: str) -> bool:
    """Compare two level_shift_report.csv. Returns True if identical."""
    a = pd.read_csv(old_path, dtype={'level_id': str}).set_index('level_id')
    b = pd.read_csv(new_path, dtype={'level_id': str}).set_index('level_id')
    print(f"level_shift_report.csv:  {len(a)} levels -> {len(b)} levels")
    common = a.index.intersection(b.index)
    if len(common) != len(a) or len(common) != len(b):
        print(f"  levels only in OLD: {len(a.index.difference(b.index))}"
              f"   only in NEW: {len(b.index.difference(a.index))}")
    a, b = a.loc[common], b.loc[common]
    identical = len(common) == len(a) == len(b)

    st_a, st_b = a['question_status'].fillna(''), b['question_status'].fillna('')
    changed = st_a != st_b
    if changed.any():
        identical = False
        print(f"  question_status changed on {int(changed.sum())} levels:")
        print(_indent(pd.crosstab(st_a[changed], st_b[changed],
                                  rownames=['old'], colnames=['new']).to_string(), 4))
    else:
        print("  question_status: unchanged")

    for col in ('dE', 'p_spur', 'n_tot', 'pattern_V'):
        if col not in a.columns:
            continue
        d = (a[col] - b[col]).abs()
        n = int((d > TOL).sum())
        if n:
            identical = False
            worst = d.idxmax()
            print(f"  {col}: {n} levels differ, largest change {d.max():.6g}"
                  f" at level {worst} ({a.loc[worst, col]:.6g} -> {b.loc[worst, col]:.6g})")
    return identical


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('old', help='directory (or .csv) of the reference run')
    ap.add_argument('new', help='directory (or .csv) of the run to check')
    args = ap.parse_args()

    ok = True
    ol = _resolve(args.old, 'line_classifications.csv')
    nl = _resolve(args.new, 'line_classifications.csv')
    ok &= compare_lines(ol, nl)

    # A .csv given explicitly names the line table alone; the shift report is
    # only looked for when both arguments are directories.
    dirs = os.path.isdir(args.old) and os.path.isdir(args.new)
    os_, ns = _resolve(args.old, 'level_shift_report.csv'), _resolve(args.new, 'level_shift_report.csv')
    if dirs and os.path.isfile(os_) and os.path.isfile(ns):
        print()
        ok &= compare_shifts(os_, ns)
    else:
        print("\nlevel_shift_report.csv: not present in both runs, skipped")

    print("\nIDENTICAL" if ok else "\nDIFFERENT")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
