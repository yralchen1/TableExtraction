"""The systematic uncertainty of a set's LOPT fit, by perturbation runs.

Work_on_hfs_plan.md, Step 9.  Run from LineClass/:

    python lopt_perturb.py --set final                  # the calibration
    python lopt_perturb.py --set final --source hfs     # the A constants
                                                        # and the kappas
    python lopt_perturb.py --set final --jobs 4 --out DIR

16 LOPT runs are made at a time (--jobs).

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

THE HYPERFINE MODEL (--source hfs)
==================================
In Set 3 a line is moved to the center of gravity of its pattern by
cg_shift = -sum BF (kappa D + R) (hfs_correction.Model.cg_shift), so an
error in an A constant or in a kappa moves every line the model moves.  A
pattern here is a change of the model: the lines of the source set
(iter_hfs/ for final/) are given their new cg_shift, and each record of the
set's lines is moved by the difference.  LOPT then gives the effect exactly:
a plain line carries only kappa of an A error and a `c` line almost none,
the ground level stays at its center of gravity, and the levels sharing a
parameter move together.  The patterns (hfs_patterns, kappa_patterns):

* the calculated constants (hfs_A_theory): one per radial parameter, moved by
  its uncertainty, the outer shells copied from it moving with it; one per
  outer shell's copying factor (u_factor); one per untested configuration,
  10% of the part of A that configuration gives, each configuration on its
  own (user, 2026-10-10);
* the semiempirical scale hfs_kappa.U_A_SCALE, on the composition constants;
* the constants fitted to Sugar's resolved components: the eigenvectors of
  their covariance in the anchored fit (hfs_A_fit), which is recomputed and
  must reproduce the table;
* what is left of each level's u_A once those are taken out: its own, one
  pattern per level (none below OWN_MIN of u_A);
* one per kappa class, moved by its uncertainty in [hfs.kappa].

The ground level's A comes from the 4f radial parameter of 4f3, so that
pattern moves nearly every level by the same amount: the energies are
measured from the ground level's center of gravity.  It cancels in a Ritz
wavenumber between excited levels, and the modes file keeps that.

WHAT IS WRITTEN (in --out, the set's directory by default)
==========================================================
calibration_sys_levels.csv  per level: the energy of the base run, LOPT's
                            D1 and D2stat, the calibration part u_cal, and
                            the part of
                            it owed to the largest pattern (the scale).
calibration_sys_modes.csv   per level, its shift in each pattern (cm^-1):
                            the matrix whose rows give any Ritz wavenumber's
                            calibration uncertainty (ritz_unc).
calibration_sys_report.txt  what was run and what came out.
hfs_sys_levels.csv          (--source hfs) per level: the energy, D1 and
                            D2stat, u_hfs
                            from the A constants with its correlated part
                            u_hfs_corr and own part u_hfs_own, and u_kappa.
hfs_sys_modes.csv           per level, its shift in each pattern, the
hfs_sys_report.txt          columns named by the pattern (radial:, factor:,
                            config:, scale, fit:, level:, kappa:).

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

#: LOPT runs at a time: one per physical core of the 7950X (each run is
#: single-threaded; past 16 the cores' second threads add little)
JOBS = 16
#: the largest record shift of a scaled pattern, cm^-1
SCALE_TO = 1.0
#: a level's own part of u_A smaller than this fraction of u_A gets no run
OWN_MIN = 0.01
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
    """`{level_id: (energy, D1, D2stat)}` of a LOPT levels file."""
    out = {}
    with open(path, encoding='utf-8', newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            out[r['Designation'].strip()] = (float(r['Energy']),
                                             float(r['D1']),
                                             float(r['D2stat']))
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
# The hyperfine model's patterns
# ---------------------------------------------------------------------------
def model_lines(source_csv):
    """`{key4: make_set3.Line}` of the source set's table: each line with its
    character, wn_key and accepted transitions, as make_set3 moves it."""
    import make_set3
    return make_set3.read_table(source_csv)[2]


def cg_shifts(model, lines, keys=None):
    """`{key: cg_shift}` of `lines` (of `keys` only, if given) under `model`;
    0 for a line without an accepted transition."""
    out = {}
    for key in (lines if keys is None else keys):
        ln = lines[key]
        out[key] = (model.cg_shift(ln.char, ln.wn_key, ln.components).shift
                    if ln.components else 0.0)
    return out


def lines_of_levels(lines):
    """`{level_id: [key]}`: the lines whose accepted transitions touch each
    level."""
    out = {}
    for key, ln in lines.items():
        for low, upp, _bf in ln.components:
            for lid in (low, upp):
                out.setdefault(lid, set()).add(key)
    return {k: sorted(v) for k, v in out.items()}


def line_moves(model, lines, base, change, by_level):
    """`{key: d cg_shift}` when `change` is applied to `model`, the model
    restored afterwards.  `change` is ('A', {level_id: dA}) or
    ('kappa', class, d kappa)."""
    if change[0] == 'A':
        saved, keys = {}, set()
        for lid, dA in change[1].items():
            if lid not in model.levels:
                continue
            saved[lid] = model.levels[lid]
            J, A, u_A, det = saved[lid]
            model.levels[lid] = (J, A + dA, u_A, det)
            keys.update(by_level.get(lid, ()))
        try:
            new = cg_shifts(model, lines, sorted(keys))
        finally:
            model.levels.update(saved)
    else:
        _, cls, dk = change
        saved = model.kappa[cls]
        model.kappa[cls] = (saved[0] + dk, saved[1])
        try:
            new = cg_shifts(model, lines)
        finally:
            model.kappa[cls] = saved
    return {k: v - base[k] for k, v in new.items() if v != base[k]}


def theory_patterns(table):
    """The correlated patterns of the calculated A constants (the rows of
    `table` whose source is hfs_A_theory's): `[(name, {level_id: dA})]` and
    `{level_id: their variance}`.

    * one per radial parameter, moved by its uncertainty, the outer shells
      copied from it moving with it (Params.apply_scaling);
    * one per outer shell's copying factor, moved by its u_factor;
    * one per untested configuration: `fraction` of the part of A that
      configuration gives (hfs_A_theory.untested_part, taken configuration
      by configuration - user, 2026-10-10).
    """
    import copy
    import hfs_A_theory as T
    params = T.Params.read()
    levels = T.read_levels()
    T.attach_level_ids(levels, log=lambda *a: None)
    sigmas, fraction, min_measured = T.read_uncertainties(T.PARAMS_FILE,
                                                         params)
    known = {lv.level_id: lv for lv in levels if lv.level_id}
    tested = T.configs_tested(list(known.values()), T.read_measured(),
                              min_measured)
    calc = [lid for lid, r in table.items()
            if r['source'].startswith(T.CALC_SOURCE) and lid in known]
    nominal = {lid: known[lid].A(params)[0] for lid in calc}

    def moved(p2):
        out = {}
        for lid in calc:
            a1 = known[lid].A(p2)[0]
            if a1 is not None and nominal[lid] is not None:
                d = a1 - nominal[lid]
                if abs(d) > 1e-12:
                    out[lid] = d
        return out

    pats = []
    for name in params.names():
        key, kind = name.rsplit('.', 1)
        sig = sigmas.get((key, kind), 0.0)
        if sig <= 0:
            continue
        p2 = copy.deepcopy(params)
        p2.set(name, params.table[key][kind] + sig)
        p2.apply_scaling()
        pats.append(('radial:' + name, moved(p2)))
    for nl, sc in params.scaled.items():
        if sc.u_factor <= 0:
            continue
        p2 = copy.deepcopy(params)
        p2.scaled[nl].factor = sc.factor * (1.0 + sc.u_factor)
        p2.apply_scaling()
        pats.append(('factor:' + nl, moved(p2)))
    confs = sorted({conf for lid in calc for (conf, _nl, _k)
                    in known[lid].theta})
    for conf in confs:
        if conf in tested:
            continue
        d = {}
        for lid in calc:
            part = sum(th * (params.value(c, nl, k) or 0.0)
                       for (c, nl, k), th in known[lid].theta.items()
                       if c == conf)
            if abs(part) > 1e-12:
                d[lid] = fraction * part
        pats.append(('config:' + conf, d))
    pats = [(n, d) for n, d in pats if d]
    var = {}
    for _n, d in pats:
        for lid, x in d.items():
            var[lid] = var.get(lid, 0.0) + x * x
    return pats, var


def fit_patterns(table):
    """The patterns of the constants fitted to Sugar's resolved components
    (rows whose source begins 'resolved components'): the eigenvectors of
    their covariance in the anchored fit (hfs_A_fit), recomputed here.
    Returns `[(name, {level_id: dA})]`, `{level_id: variance}` and the
    largest difference between the fit's values and the table's."""
    import hfs_A_fit as F
    import hfs_kappa
    rows = F.read_A_table(os.path.join(HERE, F.A_LEVELS))
    priors = F.read_priors(os.path.join(HERE, F.PRIORS))
    anchors = F.anchors_of(rows, priors=priors)
    pats_, _ = F.read_patterns(os.path.join(HERE, F.LINES), rows,
                               hfs_kappa.read_levels())
    res = F.fit(pats_, anchors)
    M, _b, _n = F._equations(pats_, anchors, res.ids)
    cov = np.linalg.pinv(M.T @ M, rcond=1e-12)
    ids = [lid for lid, r in table.items()
           if r['source'].startswith('resolved components')
           and lid in res.A]
    ix = [res.ids.index(lid) for lid in ids]
    sub = cov[np.ix_(ix, ix)]
    worst = max((abs(res.A[lid] - float(table[lid]['A_cm-1']))
                 for lid in ids), default=0.0)
    a, _ev = patterns(sub)
    out = []
    for k, row in enumerate(a):
        d = {lid: float(x) for lid, x in zip(ids, row) if abs(x) > 1e-12}
        if d:
            out.append(('fit:%02d' % (k + 1), d))
    return out, {lid: float(sub[i, i]) for i, lid in enumerate(ids)}, worst


def hfs_patterns(model, table):
    """Every pattern of the A constants, on the scale the model uses:
    `([(name, {level_id: dA})], notes for the report)`.  The correlated
    parts come first (theory_patterns, the semiempirical scale,
    fit_patterns); what is left of each level's u_A is its own, one pattern
    per level."""
    import hfs_kappa
    pats, var = theory_patterns(table)
    notes = ['calculated constants: %d correlated patterns (%d radial, %d '
             'factors, %d configurations)' % (
                 len(pats), sum(n.startswith('radial:') for n, _ in pats),
                 sum(n.startswith('factor:') for n, _ in pats),
                 sum(n.startswith('config:') for n, _ in pats))]
    scale = {}
    for lid, r in table.items():
        if hfs_kappa.is_semiempirical(r['source']):
            scale[lid] = float(r['A_cm-1']) * hfs_kappa.U_A_SCALE
    if scale:
        pats.append(('scale', scale))
        for lid, x in scale.items():
            var[lid] = var.get(lid, 0.0) + x * x
    fpats, fvar, worst = fit_patterns(table)
    pats += fpats
    for lid, v in fvar.items():
        var[lid] = var.get(lid, 0.0) + v
    notes.append('semiempirical scale: %d levels; fitted constants: %d '
                 'patterns over %d levels, the fit reproducing the table to '
                 '%.5f cm^-1' % (len(scale), len(fpats), len(fvar), worst))
    over = []
    for lid in sorted(model.levels):
        J, A, u_A, det = model.levels[lid]
        if not det or u_A <= 0:
            continue
        rest = u_A * u_A - var.get(lid, 0.0)
        if rest < -1e-10:
            over.append((math.sqrt(var[lid]) / u_A, lid))
            continue
        if rest > (OWN_MIN * u_A) ** 2:
            pats.append(('level:' + lid, {lid: math.sqrt(rest)}))
    notes.append('levels whose correlated parts reach their u_A (no own part '
                 'left): %d%s' % (len(over), (', largest ratio %.2f (%s)'
                                              % max(over)) if over else ''))
    return pats, notes


def kappa_patterns(model):
    """One pattern per fitted kappa class, moved by its uncertainty."""
    return [('kappa:' + cls, ('kappa', cls, u))
            for cls, (v, u) in sorted(model.kappa.items()) if u > 0]


def record_lines(records, keys, shifts):
    """The line (make_set3.key4 of its wn_key) of every record, None if it
    has none."""
    import make_set3
    out = []
    for _text, wn, low, upp in records:
        measured = hfs_correction.measured(low, upp, wn, shifts)
        near = [(abs(w - measured), k) for w, k in keys.get((low, upp), ())
                if abs(w - measured) <= MATCH_TOL]
        out.append(make_set3.key4(min(near)[1]) if near else None)
    return out


# ---------------------------------------------------------------------------
# The runs
# ---------------------------------------------------------------------------
def perturb(set_dir, records, named, jobs, work, check=0):
    """Run LOPT on the set's records moved by each of `named`,
    `[(name, [shift per record])]`, plus the base run, and the pattern
    `check` again at half size.  Returns `(base, modes, info)`: the base
    run's levels, `{level_id: [shift per pattern]}`, and the report lines
    of the linearity check, of the base run against the set's own fit, and
    of the count of runs and the time."""
    t0 = time.time()
    runs, factors = {}, []
    runs['base'] = prepare(work, 'base', set_dir, records,
                           [0.0] * len(records))
    for k, (_name, d) in enumerate(named):
        big = max(abs(x) for x in d)
        f = SCALE_TO / big if big > 0 else 1.0
        factors.append(f)
        runs[k] = prepare(work, 'p%04d' % (k + 1), set_dir, records,
                          [f * x for x in d])
    half = 0.5 * factors[check]
    runs['half'] = prepare(work, 'half', set_dir, records,
                           [half * x for x in named[check][1]])

    def progress(n, total):
        print('\r  %d of %d LOPT runs done' % (n, total), end='', flush=True)

    print('%d LOPT runs, %d at a time, in %s' % (len(runs), jobs, work))
    done = run_all(jobs, runs, progress=progress)
    print()
    base = done.pop('base')
    half_run = done.pop('half')
    modes = responses(base, [done[k] for k in range(len(named))], factors)
    lin = max(abs((half_run[lid][0] - base[lid][0]) / half
                  - modes[lid][check]) for lid in base)
    lin_ref = max(abs(m[check]) for m in modes.values())
    own = read_levels(os.path.join(set_dir, LEVELS_OUT))
    base_diff = max(abs(base[lid][0] - own[lid][0])
                    for lid in base if lid in own)
    info = [
        'base run against the set\'s own LOPT levels: largest |diff| %.2e '
        'cm^-1' % base_diff,
        'linearity: the pattern %s at half size differs by %.2e cm^-1 at '
        'most (its largest level shift %.4f)' % (named[check][0], lin,
                                                 lin_ref),
        '%d runs, %d at a time, %.0f s' % (len(runs), jobs, time.time() - t0)]
    return base, modes, info


def write_rows(path, rows):
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def write_modes(path, base, names, modes):
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, lineterminator='\n')
        w.writerow(['level_id'] + list(names))
        for lid in base:
            w.writerow([lid] + ['%+.6e' % x for x in modes[lid]])


def summary(rows, col):
    """Report lines: the median and spread of `col` against D2stat, LOPT's
    statistical uncertainty, by energy."""
    u = np.array([float(r[col]) for r in rows])
    d2 = np.array([float(r['D2stat']) for r in rows])
    E = np.array([float(r['energy']) for r in rows])
    hi = E > 1.0
    text = ['%s over %d levels (ground excluded): median %.4f, 90%% %.4f, '
            'max %.4f cm^-1; / D2stat median %.2f, 90%% %.2f'
            % (col, hi.sum(), np.median(u[hi]), np.percentile(u[hi], 90),
               u[hi].max(), np.median(u[hi] / d2[hi]),
               np.percentile(u[hi] / d2[hi], 90))]
    for lo, up in ((0, 50000), (50000, 100000), (100000, 1e9)):
        sel = hi & (E >= lo) & (E < up)
        if sel.any():
            text.append('  E %6d-%-6s: %3d levels, median %s %.4f, median '
                        'D2stat %.4f' % (lo, '' if up > 1e8 else '%d' % up,
                                         sel.sum(), col, np.median(u[sel]),
                                         np.median(d2[sel])))
    return text


# ---------------------------------------------------------------------------
def run_calibration(args, set_dir, cfg, records, work):
    """The calibration: one pattern per eigenvector of the groups'
    covariance."""
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
    named = [('p%03d' % (k + 1), record_shifts(records, rec_group, index,
                                               a[k])) for k in range(n_pat)]
    base, modes, info = perturb(set_dir, records, named, args.jobs, work)
    rows = []
    for lid in base:
        u = math.sqrt(sum(x * x for x in modes[lid]))
        rows.append({'level_id': lid, 'energy': '%.4f' % base[lid][0],
                     'D1': '%.4f' % base[lid][1],
                     'D2stat': '%.4f' % base[lid][2], 'u_cal': '%.4f' % u,
                     'u_cal_scale': '%.4f' % abs(modes[lid][0])})
    text = [
        'corrections %s; covariance %s' % (args.corrections, args.covariance),
        'groups %d; patterns run %d of %d, the largest carrying %.1f%% of '
        'the variance' % (len(names), n_pat, len(a), 100 * ev[0] / ev.sum()),
        'corrections against covariance: largest relative difference of a '
        'group\'s u %.2f%%' % (100 * worst),
        'records not moved: %d not in the table, %d without a group, %d in '
        'a group the covariance does not hold%s'
        % (missing['not in the table'], missing['no group'],
           sum(1 for g in rec_group if g in unknown),
           (' (' + ', '.join(unknown) + ')') if unknown else '')] + info
    text += [''] + summary(rows, 'u_cal')
    text.append('the largest pattern\'s share of u_cal^2: median %.2f'
                % np.median([(float(r['u_cal_scale']) / float(r['u_cal'])) ** 2
                             for r in rows if float(r['u_cal']) > 0]))
    return rows, [n for n, _ in named], base, modes, text


def run_hfs(args, set_dir, cfg, records, work):
    """The A constants and the kappa classes: each pattern a change of the
    model, every line of the source set moved to its new center of gravity
    (Model.cg_shift, as make_set3.py moves it), and the records of the set's
    lines moved by the difference."""
    import make_set3
    source = (make_set3.SOURCE_SET if config.derived_from(set_dir)
              else set_dir)
    src_cfg = config.load(os.path.join(source, config.CONFIG_NAME))
    model = hfs_correction.Model(src_cfg.hfs)
    with open(src_cfg.hfs.A_levels, encoding='utf-8', newline='') as fh:
        table = {r['level_id'].strip(): r for r in csv.DictReader(fh)}
    lines = model_lines(src_cfg.output_csv)
    base_shift = cg_shifts(model, lines)
    by_level = lines_of_levels(lines)
    shifts = hfs_correction.read_shifts(os.path.join(set_dir, LINES_IN))
    rec_line = record_lines(records, table_keys(cfg.output_csv), shifts)
    pats, notes = hfs_patterns(model, table)
    pats = [(n, ('A', d)) for n, d in pats] + kappa_patterns(model)
    if args.limit is not None:
        pats = pats[:args.limit]
    named, dropped = [], 0
    for name, change in pats:
        moves = line_moves(model, lines, base_shift, change, by_level)
        d = [moves.get(k, 0.0) if k is not None else 0.0 for k in rec_line]
        if max((abs(x) for x in d), default=0.0) < 1e-9:
            dropped += 1
            continue
        named.append((name, d))
    base, modes, info = perturb(set_dir, records, named, args.jobs, work)
    names = [n for n, _ in named]
    groups = {'corr': [k for k, n in enumerate(names)
                       if not n.startswith(('level:', 'kappa:'))],
              'own': [k for k, n in enumerate(names)
                      if n.startswith('level:')],
              'kappa': [k for k, n in enumerate(names)
                        if n.startswith('kappa:')]}
    rows = []
    for lid in base:
        m = modes[lid]
        part = {g: math.sqrt(sum(m[k] ** 2 for k in ks))
                for g, ks in groups.items()}
        rows.append({'level_id': lid, 'energy': '%.4f' % base[lid][0],
                     'D1': '%.4f' % base[lid][1],
                     'D2stat': '%.4f' % base[lid][2],
                     'u_hfs': '%.4f' % math.hypot(part['corr'], part['own']),
                     'u_hfs_corr': '%.4f' % part['corr'],
                     'u_hfs_own': '%.4f' % part['own'],
                     'u_kappa': '%.4f' % part['kappa']})
    text = ['source set %s; A table %s' % (source, src_cfg.hfs.A_levels)]
    text += notes
    text += ['records without a line: %d; patterns that move no record: %d'
             % (sum(1 for k in rec_line if k is None), dropped),
             'patterns run: %d (%d correlated, %d own, %d kappa)'
             % (len(named), len(groups['corr']), len(groups['own']),
                len(groups['kappa']))] + info
    text += ['']
    for col in ('u_hfs', 'u_hfs_corr', 'u_hfs_own', 'u_kappa'):
        text += summary(rows, col)
    for k in groups['kappa']:
        u = np.array([abs(modes[lid][k]) for lid in base])
        text.append('  %s: median %.4f, max %.4f cm^-1'
                    % (names[k], np.median(u), u.max()))
    return rows, names, base, modes, text


OUTPUTS = {'calibration': ('calibration_sys_levels.csv',
                           'calibration_sys_modes.csv',
                           'calibration_sys_report.txt'),
           'hfs': ('hfs_sys_levels.csv', 'hfs_sys_modes.csv',
                   'hfs_sys_report.txt')}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--source', choices=sorted(OUTPUTS), default='calibration',
                   help='what is perturbed (default: %(default)s)')
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
    out_paths = [os.path.join(out_dir, n) for n in OUTPUTS[args.source]]
    import output_files
    output_files.require_writable(out_paths)

    records = read_records(os.path.join(set_dir, LINES_IN))
    work = args.work or tempfile.mkdtemp(prefix='lopt_perturb_')
    os.makedirs(work, exist_ok=True)
    run = run_calibration if args.source == 'calibration' else run_hfs
    rows, names, base, modes, body = run(args, set_dir, cfg, records, work)

    write_rows(out_paths[0], rows)
    write_modes(out_paths[1], base, names, modes)
    text = ['lopt_perturb.py --source %s, %s' % (
                args.source,
                datetime.datetime.now().strftime('%Y-%m-%d %H:%M')),
            'set: %s (%d records)' % (set_dir, len(records))] + body
    text.append('written: %s' % ', '.join(out_paths))
    with open(out_paths[2], 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(text) + '\n')
    print('\n'.join(text))
    if not args.keep and not args.work:
        shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
