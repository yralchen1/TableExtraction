#!/usr/bin/env python
"""Iterate the Boltzmann fit and the classification to a fixed point.

Why an iteration is needed
--------------------------
The two constants of the intensity model, C and kT, are fitted on the
identifications the classification accepted; and which identifications the
classification accepts depends on the calculated intensities, hence on C and
kT.  Fitting once therefore leaves the model and the line list disagreeing
with each other.  This script closes the loop:

    fit C, kT on the accepted lines
      -> recompute the whole Icalc column with them
        -> reclassify the observed lines against the new Icalc
          -> fit again on the new set of accepted lines
            -> ...

until the calculated intensities stop moving.

The convergence criterion
-------------------------
The test is on the quantity the rest of the pipeline actually uses, the
calculated intensities themselves, and not on C and kT:

    max over all transitions of | Icalc_now / Icalc_previous - 1 |  <  tol

with tol = 0.01 by default, i.e. no calculated intensity in the file moved by
as much as one per cent in the last round.  (C and kT can shift and still
leave Icalc almost unchanged, because a larger C compensates a smaller kT over
the energy range where most of the lines are; the criterion above cannot be
fooled that way, since it looks at every transition, including those at the
ends of the energy range where such a compensation fails.)

Each round ends with a classification run made with the Icalc that round
produced, so the files left behind are always mutually consistent.

Usage:
    python tools/iterate_boltzmann.py                 # iterate to convergence
    python tools/iterate_boltzmann.py --tol 0.005 --max-rounds 15
    python tools/iterate_boltzmann.py --dry-run       # one fit, change nothing
"""
import argparse
import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
import config          # noqa: E402
import fit_boltzmann   # noqa: E402


def fit_once(cfg, lines_csv, clip):
    """The two-pass Boltzmann fit; returns (C, kT, n_points, rms)."""
    df, n_all = fit_boltzmann.load_points(cfg, lines_csv)
    if len(df) < 10:
        raise SystemExit(f'{lines_csv}: only {len(df)} fittable points')
    C1, kT1, _, _, _ = fit_boltzmann.straight_line(df['Eup'].to_numpy(),
                                                   df['y'].to_numpy())
    pred = fit_boltzmann.model_intensity(C1, kT1, df['gA'], df['Eup'],
                                         df['rwn'])
    keep = np.abs(np.log(pred / df['obs_intens'])) <= clip
    k = df[keep]
    C2, kT2, _, _, rms2 = fit_boltzmann.straight_line(k['Eup'].to_numpy(),
                                                      k['y'].to_numpy())
    return C2, kT2, n_all, len(df), int(keep.sum()), rms2


def set_constants(path, C, kT):
    """Write C and kT into the [intensity_model] table of the config file."""
    with open(path, encoding='utf-8') as fh:
        text = fh.read()
    head, sep, tail = text.partition('[intensity_model]')
    if not sep:
        raise SystemExit(f'{path}: no [intensity_model] table')
    tail = re.sub(r'(?m)^C\s*=.*$', f'C  = {C:.6f}', tail, count=1)
    tail = re.sub(r'(?m)^kT\s*=.*$', f'kT = {kT:.6f}', tail, count=1)
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(head + sep + tail)


def classify(cfg_path, log_path):
    """One classification run; raises if it fails."""
    with open(log_path, 'w', encoding='utf-8') as log:
        p = subprocess.run([sys.executable, 'classify_lines.py',
                            '--config', cfg_path],
                           cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    if p.returncode != 0:
        raise SystemExit(f'classify_lines.py failed (see {log_path})')


def accepted_count(csv_path):
    import pandas as pd
    d = pd.read_csv(csv_path, usecols=['accepted'])
    return int((d['accepted'].fillna(0) == 1).sum())


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default=config.DEFAULT_PATH)
    ap.add_argument('--tol', type=float, default=0.01,
                    help='converged when no Icalc moves by more than this '
                         'fraction in a round (default 0.01)')
    ap.add_argument('--max-rounds', type=int, default=12)
    ap.add_argument('--clip', type=float, default=2.0,
                    help='|ln(Icalc/Iobs)| beyond which a point is dropped '
                         'before the second pass of each fit (default 2)')
    ap.add_argument('--dry-run', action='store_true',
                    help='fit once, print, and change nothing')
    a = ap.parse_args(argv)

    logdir = os.path.join(ROOT, 'boltzmann_iterations')
    os.makedirs(logdir, exist_ok=True)
    history = []
    for rnd in range(1, a.max_rounds + 1):
        cfg = config.load(a.config)
        t0 = time.time()
        C, kT, n_acc, n_fit, n_kept, rms = fit_once(cfg, cfg.output_csv,
                                                    a.clip)
        print(f'round {rnd}: {n_acc} accepted lines, {n_fit} with a gA, '
              f'{n_kept} kept -> C = {C:.6g}, kT = {kT:.6g} cm^-1, '
              f'rms ln I = {rms:.4g}', flush=True)
        if a.dry_run:
            return 0

        _, n, old, new, _ = fit_boltzmann.rewrite_icalc(cfg, C, kT)
        dev = float(np.max(np.abs(new / old - 1.0)))
        med = float(np.median(np.abs(new / old - 1.0)))
        print(f'         Icalc rewritten ({n} values): max |change| '
              f'{dev:.4g}, median {med:.4g}', flush=True)
        set_constants(a.config, C, kT)

        classify(a.config, os.path.join(logdir, f'classify_round{rnd}.log'))
        n_after = accepted_count(cfg.output_csv)
        print(f'         reclassified: {n_after} accepted identifications '
              f'({time.time() - t0:.0f} s)', flush=True)
        history.append((rnd, C, kT, dev, n_after))
        if dev < a.tol:
            print(f'\nconverged after {rnd} rounds: no calculated intensity '
                  f'moved by more than {dev:.3g} (< {a.tol:g}) in the last '
                  f'round.', flush=True)
            break
    else:
        print(f'\nNOT converged in {a.max_rounds} rounds; last max |change| '
              f'{history[-1][3]:.4g}', flush=True)

    print('\nround      C          kT (cm^-1)   max|dIcalc|   accepted')
    for rnd, C, kT, dev, n_after in history:
        print(f'{rnd:5d}  {C:9.4f}   {kT:10.3f}   {dev:11.4g}   {n_after:8d}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
