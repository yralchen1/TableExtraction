"""The systematic uncertainty of a set's LOPT fit, by perturbation runs.

Work_on_hfs_plan.md, Step 9.  Run from LineClass/:

    python lopt_perturb.py --set final                  # 16 LOPT runs at a time
    python lopt_perturb.py --set final --jobs 4 --out DIR

WHY NOT LOPT'S OWN GROUPS
=========================
LOPT carries a systematic error as groups of lines: each group is shifted by
one random number times a fixed function, the groups independent of one
another, and the effect on each level is the derivative of the solution
(Kramida 2011, Eqs. 22-26).  The calibration of Sugar's wavelengths
(wavelength_calibration.py) fits a constant d_lambda in each of 151 groups,
together with the level energies, and its covariance
(wavelength_calibration_cov.csv) says the groups are far from independent:
all of them positively correlated, and a quarter of the variance in one
pattern, d_lambda proportional to lambda - an error of the whole wavelength
scale, which the level system can hardly tell from a common factor on every
energy.  Independent groups would miss most of it.

WHAT THIS PROGRAM DOES
======================
The fit is linear in the measured wavenumbers, so the propagation can be made
exact without a Monte Carlo.  The covariance C of the groups' d_lambda is
split into its eigenvectors, C = sum_k lambda_k v_k v_k^T: 151 independent
patterns, pattern k shifting group g by a_kg = sqrt(lambda_k) v_kg angstrom.
For each pattern, every LOPT record of a line in group g is moved by

    d_wn = -a_kg * wn^2 * 1e-8          (a wavelength too long is a
                                         wavenumber too small)

LOPT is run on the moved records, and the shift of every level is that
pattern's 1-sigma effect.  The calibration part of a level's uncertainty is
the quadrature sum over the patterns, and that of a Ritz wavenumber is the
quadrature sum of the differences of its two levels' shifts - the patterns
are independent, the two levels of one pattern are not.  The calibration's
scale is kept as fitted: the scale pattern is one of the 151, and its
uncertainty enters with the rest.

LOPT prints its input wavenumbers to 3 decimals and a small pattern would
vanish in that rounding.  Each pattern is therefore scaled so that its
largest record moves by SCALE_TO cm^-1, the records are written with every
decimal the 12-character field holds, and the levels' shifts are divided
back; the fit is linear, so this changes nothing but the precision.  One
pattern is also run at half that size, and the report gives how far the two
disagree: a check that the fit really is linear.

Which group a record's line is in comes from the corrections file of the
same calibration run (wavelength_calibration_corrections.csv, by Sugar's
wavenumber, the line's wn_key), through the set's own classification table
and its LOPT_hfs_shifts.txt.  A line without a group, or in a group the
covariance does not hold, is not moved and is counted in the report.  The
corrections and the covariance must come from the same run: each line's
u_d_lambda_A is checked against its group's variance.

WHAT IS WRITTEN (in --out, the set's directory by default)
==========================================================
calibration_sys_levels.csv  per level: the energy of the base run, LOPT's
                            D1, the calibration part u_cal, and the part of
                            it owed to the largest pattern (the scale).
calibration_sys_modes.csv   per level, its shift in each pattern (cm^-1):
                            the matrix whose rows give any Ritz wavenumber's
                            calibration uncertainty (ritz_unc).
calibration_sys_report.txt  what was run and what came out.

The LOPT runs are made in a temporary directory (--work), one subdirectory
per run, removed afterwards unless --keep.  The set's own files are only
read.
"""
import argparse
import concurrent.futures
import csv
import datetime
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

import config
import hfs_correction
import make_LOPT_input

HERE = os.path.dirname(os.path.abspath(__file__))
CORRECTIONS = os.path.join(HERE, 'wavelength_calibration_corrections.csv')
COVARIANCE = os.path.join(HERE, 'wavelength_calibration_cov.csv')
LINES_IN = 'LOPT_input_lines.txt'
FIXLEV = 'LOPT_fixlev.txt'
PAR = 'LOPT.par'
LEVELS_OUT = 'LOPT_output_levels.txt'
LINES_OUT = 'LOPT_output_lines.txt'
OUT_LEVELS = 'calibration_sys_levels.csv'
OUT_MODES = 'calibration_sys_modes.csv'
OUT_REPORT = 'calibration_sys_report.txt'

#: LOPT runs at a time: one per physical core of the 7950X (each run is
#: single-threaded; past 16 the cores' second threads add little)
JOBS = 16
#: the largest record shift of a scaled pattern, cm^-1
SCALE_TO = 1.0
#: the wavenumber field of a LOPT record, 1-based inclusive columns
WN_FIELD = make_LOPT_input.FIELDS['wavenumber']
LOW_FIELD = make_LOPT_input.FIELDS['lower_level']
UPP_FIELD = make_LOPT_input.FIELDS['upper_level']
#: a line's u_d_lambda_A and its group's sqrt(variance) must agree to this,
#: relative, beyond the 5-decimal rounding of the corrections file
COV_CHECK_TOL = 0.02


# ---------------------------------------------------------------------------
# The records
# ---------------------------------------------------------------------------
def field(text, span):
    first, last = span
    return text[first - 1:last].strip()


def read_records(path):
    """The records of a LOPT transitions file: `[(text, wn, low, upp)]`,
    the text without its line ending."""
    out = []
    with open(path, encoding='ascii', newline='') as fh:
        for raw in fh:
            text = raw.rstrip('\r\n')
            if not text.strip():
                continue
            out.append((text, float(field(text, WN_FIELD)),
                        field(text, LOW_FIELD), field(text, UPP_FIELD)))
    return out


def format_wn(wn):
    """`wn` with every decimal the wavenumber field holds."""
    width = WN_FIELD[1] - WN_FIELD[0] + 1
    whole = len('%d' % abs(int(wn))) + (1 if wn < 0 else 0)
    decimals = width - whole - 1
    if decimals < 3:
        raise ValueError('%r does not fit in the %d-character wavenumber '
                         'field with 3 decimals' % (wn, width))
    text = '%.*f' % (decimals, wn)
    return text if len(text) <= width else '%.*f' % (decimals - 1, wn)


def with_wn(text, wn):
    """The record `text` with its wavenumber field replaced by `wn`."""
    first, last = WN_FIELD
    width = last - first + 1
    new = format_wn(wn).ljust(width)
    return text[:first - 1] + new + text[last:]


def write_records(records, shifts, path):
    """Write the transitions file with record i moved by `shifts[i]`."""
    with open(path, 'w', encoding='ascii', newline='') as fh:
        for (text, wn, _low, _upp), d in zip(records, shifts):
            fh.write(with_wn(text, wn + d) + make_LOPT_input.EOL_LINES)


# ---------------------------------------------------------------------------
# Which calibration group each record's line is in
# ---------------------------------------------------------------------------
def read_groups(path):
    """`{'%.4f' % wn_key: (group, u_d_lambda_A or None)}` from the
    corrections file; lines without a group are left out."""
    out = {}
    with open(path, encoding='utf-8-sig', newline='') as fh:
        for r in csv.DictReader(fh):
            g = (r.get('group') or '').strip()
            if not g:
                continue
            u = (r.get('u_d_lambda_A') or '').strip()
            out['%.4f' % float(r['wn_obs'])] = (g, float(u) if u else None)
    return out


def read_covariance(path):
    """`(names, C)`: the group names and their covariance, angstrom^2."""
    with open(path, encoding='utf-8', newline='') as fh:
        rows = list(csv.reader(fh))
    names = rows[0][1:]
    if [r[0] for r in rows[1:]] != names:
        raise SystemExit('%s: the row names are not the column names' % path)
    C = np.array([[float(x) for x in r[1:]] for r in rows[1:]])
    return names, 0.5 * (C + C.T)


def table_keys(path):
    """`{(low_id, upp_id): [(wn_obs, wn_key)]}` of a classification table."""
    out = {}
    with open(path, encoding='utf-8', newline='') as fh:
        for r in csv.DictReader(fh):
            out.setdefault((r['low_id'].strip(), r['upp_id'].strip()),
                           []).append((float(r['wn_obs']),
                                       float(r['wn_key'])))
    return out


#: how far a record's wavenumber may lie from its table row's wn_obs: the
#: 3-decimal rounding of make_LOPT_input.py
MATCH_TOL = 0.0006


def record_groups(records, keys, shifts, groups):
    """The group of every record (None if it has none) and the counts of
    those without one, by reason."""
    out, missing = [], {'not in the table': 0, 'no group': 0}
    for _text, wn, low, upp in records:
        measured = hfs_correction.measured(low, upp, wn, shifts)
        near = [(abs(w - measured), k) for w, k in keys.get((low, upp), ())
                if abs(w - measured) <= MATCH_TOL]
        key = min(near)[1] if near else None
        if key is None:
            out.append(None)
            missing['not in the table'] += 1
            continue
        g = groups.get('%.4f' % key)
        out.append(None if g is None else g[0])
        if g is None:
            missing['no group'] += 1
    return out, missing


def check_covariance(groups, names, C):
    """Raise unless each line's u_d_lambda_A is its group's sqrt(variance):
    the corrections and the covariance come from one calibration run."""
    index = {n: k for k, n in enumerate(names)}
    worst, bad = 0.0, []
    for g, u in groups.values():
        if g not in index or u is None or u <= 0:
            continue
        s = math.sqrt(C[index[g], index[g]])
        worst = max(worst, abs(s - u) / u)
        if abs(s - u) > COV_CHECK_TOL * u + 0.6e-5:
            bad.append(g)
    if bad:
        raise SystemExit('the corrections and the covariance disagree in %d '
                         'group(s), e.g. %s (up to %.1f%%): they are not '
                         'from the same calibration run'
                         % (len(set(bad)), bad[0], 100 * worst))
    return worst


def patterns(C):
    """The independent patterns of the covariance `C`: rows a_k with
    sum_k a_k a_k^T = C, largest first."""
    ev, V = np.linalg.eigh(C)
    order = np.argsort(ev)[::-1]
    ev = np.clip(ev[order], 0.0, None)
    V = V[:, order]
    return (V * np.sqrt(ev)).T, ev


def record_shifts(records, rec_group, index, a):
    """The wavenumber shift of every record in the pattern `a` (angstrom per
    group): -a_g * wn^2 * 1e-8, zero for a record without a group."""
    out = []
    for (_text, wn, _low, _upp), g in zip(records, rec_group):
        k = index.get(g) if g is not None else None
        out.append(0.0 if k is None else -a[k] * wn * wn * 1e-8)
    return out


# ---------------------------------------------------------------------------
# Running LOPT
# ---------------------------------------------------------------------------
def read_levels(path):
    """`{level_id: (energy, D1)}` of a LOPT levels file."""
    out = {}
    with open(path, encoding='utf-8', newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            out[r['Designation'].strip()] = (float(r['Energy']),
                                             float(r['D1']))
    return out


def run_lopt(run_dir):
    """Run `lopt.bat LOPT.par` in `run_dir`; returns the levels."""
    with open(os.path.join(run_dir, 'lopt.log'), 'w') as log:
        subprocess.run(['cmd', '/c', 'lopt.bat', PAR], cwd=run_dir,
                       stdout=log, stderr=subprocess.STDOUT, check=False)
    path = os.path.join(run_dir, LEVELS_OUT)
    if not os.path.exists(path):
        raise RuntimeError('LOPT wrote no levels in %s (see lopt.log)'
                           % run_dir)
    return read_levels(path)


def prepare(work, name, set_dir, records, shifts):
    """A run directory: the set's fixed levels, its parameter file (with
    `stat`, make_LOPT_input.write_par_file), and the moved records."""
    run_dir = os.path.join(work, name)
    os.makedirs(run_dir)
    shutil.copy(os.path.join(set_dir, FIXLEV), run_dir)
    make_LOPT_input.write_par_file(
        os.path.join(set_dir, PAR), os.path.join(run_dir, PAR),
        LINES_IN, FIXLEV, LEVELS_OUT, LINES_OUT)
    write_records(records, shifts, os.path.join(run_dir, LINES_IN))
    return run_dir


def run_all(jobs, runs, runner=run_lopt, progress=None):
    """Run every `{name: run_dir}`, `jobs` at a time; `{name: levels}`."""
    out = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(runner, d): n for n, d in runs.items()}
        for f in concurrent.futures.as_completed(futures):
            out[futures[f]] = f.result()
            if progress:
                progress(len(out), len(runs))
    return out


def responses(base, runs, factors):
    """`{level_id: [shift per pattern]}`: each run's energies less the base
    run's, divided back by the factor its pattern was scaled by."""
    ids = list(base)
    out = {}
    for lid in ids:
        out[lid] = [(runs[k][lid][0] - base[lid][0]) / f
                    for k, f in enumerate(factors)]
    return out


def ritz_unc(modes, low, upp):
    """The calibration part of the uncertainty of the Ritz wavenumber
    E(upp) - E(low), from the patterns' level shifts."""
    a, b = modes[low], modes[upp]
    return math.sqrt(sum((y - x) ** 2 for x, y in zip(a, b)))


# ---------------------------------------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--set', default='final', dest='set_dir',
                   help='the set whose LOPT input is perturbed '
                        '(default: %(default)s)')
    p.add_argument('--corrections', default=CORRECTIONS)
    p.add_argument('--covariance', default=COVARIANCE)
    p.add_argument('--jobs', type=int, default=JOBS,
                   help='LOPT runs at a time (default: %(default)s)')
    p.add_argument('--out', default=None,
                   help="where the results go (default: the set's directory)")
    p.add_argument('--work', default=None,
                   help='where the runs are made (default: a temporary '
                        'directory)')
    p.add_argument('--keep', action='store_true',
                   help='keep the run directories')
    p.add_argument('--limit', type=int, default=None,
                   help='run only the first N patterns (a test)')
    args = p.parse_args(argv)

    set_dir = os.path.join(HERE, args.set_dir) \
        if not os.path.isabs(args.set_dir) else args.set_dir
    cfg = config.load(os.path.join(set_dir, config.CONFIG_NAME))
    out_dir = args.out or set_dir
    out_paths = [os.path.join(out_dir, n)
                 for n in (OUT_LEVELS, OUT_MODES, OUT_REPORT)]
    import output_files
    output_files.require_writable(out_paths)

    records = read_records(os.path.join(set_dir, LINES_IN))
    shifts = hfs_correction.read_shifts(os.path.join(set_dir, LINES_IN))
    groups = read_groups(args.corrections)
    names, C = read_covariance(args.covariance)
    worst = check_covariance(groups, names, C)
    rec_group, missing = record_groups(
        records, table_keys(cfg.output_csv), shifts, groups)
    index = {n: k for k, n in enumerate(names)}
    unknown = sorted({g for g in rec_group if g is not None
                      and g not in index})
    a, ev = patterns(C)
    n_pat = len(a) if args.limit is None else min(args.limit, len(a))

    work = args.work or tempfile.mkdtemp(prefix='lopt_perturb_')
    os.makedirs(work, exist_ok=True)
    t0 = time.time()
    runs, factors = {}, []
    runs['base'] = prepare(work, 'base', set_dir, records,
                           [0.0] * len(records))
    for k in range(n_pat):
        d = record_shifts(records, rec_group, index, a[k])
        big = max(abs(x) for x in d)
        f = SCALE_TO / big if big > 0 else 1.0
        factors.append(f)
        runs[k] = prepare(work, 'p%03d' % (k + 1), set_dir, records,
                          [f * x for x in d])
    # the linearity check: the largest pattern again, at half the size
    d = record_shifts(records, rec_group, index, a[0])
    half = 0.5 * factors[0]
    runs['half'] = prepare(work, 'half', set_dir, records,
                           [half * x for x in d])

    def progress(n, total):
        print('\r  %d of %d LOPT runs done' % (n, total), end='', flush=True)

    print('%d LOPT runs, %d at a time, in %s' % (len(runs), args.jobs, work))
    done = run_all(args.jobs, runs, progress=progress)
    print()
    base = done.pop('base')
    half_run = done.pop('half')
    modes = responses(base, [done[k] for k in range(n_pat)], factors)
    lin = max(abs((half_run[lid][0] - base[lid][0]) / half - modes[lid][0])
              for lid in base)
    lin_ref = max(abs(m[0]) for m in modes.values())

    set_levels = read_levels(os.path.join(set_dir, LEVELS_OUT))
    base_diff = max(abs(base[lid][0] - set_levels[lid][0])
                    for lid in base if lid in set_levels)

    rows = []
    for lid in base:
        u = math.sqrt(sum(x * x for x in modes[lid]))
        rows.append({'level_id': lid, 'energy': '%.4f' % base[lid][0],
                     'D1': '%.4f' % base[lid][1], 'u_cal': '%.4f' % u,
                     'u_cal_scale': '%.4f' % abs(modes[lid][0])})
    with open(out_paths[0], 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    with open(out_paths[1], 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, lineterminator='\n')
        w.writerow(['level_id'] + ['p%03d' % (k + 1) for k in range(n_pat)])
        for lid in base:
            w.writerow([lid] + ['%+.6e' % x for x in modes[lid]])

    u_cal = np.array([float(r['u_cal']) for r in rows])
    d1 = np.array([float(r['D1']) for r in rows])
    E = np.array([float(r['energy']) for r in rows])
    hi = E > 1.0
    text = [
        'lopt_perturb.py, %s' % datetime.datetime.now().strftime(
            '%Y-%m-%d %H:%M'),
        'set: %s (%d records); corrections %s; covariance %s'
        % (set_dir, len(records), args.corrections, args.covariance),
        'groups %d; patterns run %d of %d, the largest carrying %.1f%% of '
        'the variance' % (len(names), n_pat, len(a), 100 * ev[0] / ev.sum()),
        'corrections against covariance: largest relative difference of a '
        'group\'s u %.2f%%' % (100 * worst),
        'records not moved: %d not in the table, %d without a group, %d in '
        'a group the covariance does not hold%s'
        % (missing['not in the table'], missing['no group'],
           sum(1 for g in rec_group if g in unknown),
           (' (' + ', '.join(unknown) + ')') if unknown else ''),
        'base run against the set\'s own LOPT levels: largest |diff| %.2e '
        'cm^-1' % base_diff,
        'linearity: the largest pattern at half size differs by %.2e '
        'cm^-1 at most (its largest level shift %.4f)' % (lin, lin_ref),
        '%d runs, %d at a time, %.0f s' % (len(runs), args.jobs,
                                          time.time() - t0),
        '',
        'u_cal over %d levels (ground excluded): median %.4f, 90%% %.4f, '
        'max %.4f cm^-1' % (hi.sum(), np.median(u_cal[hi]),
                            np.percentile(u_cal[hi], 90), u_cal[hi].max()),
        'u_cal / D1: median %.2f, 90%% %.2f' % (
            np.median(u_cal[hi] / d1[hi]),
            np.percentile(u_cal[hi] / d1[hi], 90)),
        'the largest pattern\'s share of u_cal^2: median %.2f'
        % np.median([(float(r['u_cal_scale']) / float(r['u_cal'])) ** 2
                     for r in rows if float(r['u_cal']) > 0]),
    ]
    by_E = sorted(zip(E[hi], u_cal[hi], d1[hi]))
    for lo, up in ((0, 50000), (50000, 100000), (100000, 1e9)):
        sel = [(u, d) for e, u, d in by_E if lo <= e < up]
        if sel:
            text.append('  E %6d-%-6s: %3d levels, median u_cal %.4f, '
                        'median D1 %.4f' % (
                            lo, '' if up > 1e8 else '%d' % up, len(sel),
                            float(np.median([u for u, _ in sel])),
                            float(np.median([d for _, d in sel]))))
    text.append('written: %s' % ', '.join(out_paths))
    with open(out_paths[2], 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(text) + '\n')
    print('\n'.join(text))
    if not args.keep and not args.work:
        shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
