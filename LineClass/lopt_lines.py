"""Read a run from LOPT's line-output file instead of from classify_lines.py.

Why this exists.  The validation report (level_shifts.py) describes one "run":
a set of accepted identifications together with the level energies the
least-squares optimization produced from them.  Normally both come from
classify_lines.py, which does its own weighted least squares internally
(optimize_levels()).  But identifications are also revised by hand - in IDEN2,
or by editing the LOPT input file directly - and LOPT is then re-run on the
revised list.  This module turns LOPT's own line output back into the table
classify_lines.py would have written, so that the same report can be built on
the revised run.

Everything the report needs is in LOPT's line output:

    accepted  - the line entered the fit, i.e. its Weight column is > 0.
                Rows flagged 'P' carry Weight 0: they are predicted (Ritz)
                positions of transitions LOPT was asked to print but not to
                fit, i.e. rejected or merely listed candidates.  Rows flagged
                'c' are the components of a resolved blend and carry their
                branching fraction as the weight; they are accepted.
    low_E, upp_E - the optimized level energies.  By default they are NOT
                LOPT's: they are re-solved from the accepted lines with the
                pipeline's own weighted least squares (refit_energies), so
                that the energy shift dE the report tests is measured with
                the same instrument as the decoy and chance calibrations it
                is tested against.  Pass energies='lopt' to use LOPT's own
                values instead - from its level-output file if one is given
                (four decimals for every level), otherwise from the E1/E2
                columns of the line output, which are rounded to each level's
                own precision.
    new       - not in LOPT's output; recovered from the observed-line
                workbook.  An identification is "old" (new = 0) when the same
                pair of levels stands in the id1/id2 columns of the line list,
                which hold the legacy (published) identifications, against a
                line at the same wavenumber; otherwise it is new.  Level pairs
                are unique across the line list - no pair is the legacy
                identification of two different lines - so the pair alone is a
                safe key, and the wavenumber is checked only as a guard.
    wn_obs, obs_intens - LOPT prints both, but rounded to the precision of the
                line's uncertainty.  The exact values are taken from the line
                workbook instead, matching on the rounded wavenumber.

Columns of the classify_lines.py table that LOPT cannot supply (the graded
candidate's calculated intensity, its grade, the weeding notes) are written as
blanks.  They are used only by the --detail printout, which says so.

Checked against the classify run it was produced from: for
LOPT_output_lines.txt the reconstruction reproduces all 4891 accepted
identifications of line_classifications.csv, all 1476 of them flagged new,
with no disagreement on a single row.

Usage:
    import lopt_lines
    # calibration-consistent energies (what level_shifts.py wants):
    real = lopt_lines.read_run('LOPT_output_lines_revised.txt')
    # LOPT's own optimized energies:
    real = lopt_lines.read_run('LOPT_output_lines_revised.txt',
                               energies='lopt',
                               levels_path='LOPT_output_levels_revised.txt')
"""
import csv
import os

import numpy as np
import openpyxl
import pandas as pd

import classify_lines as cl

# The columns of line_classifications.csv, in order; those LOPT cannot supply
# are filled with NaN so that the table stays interchangeable with the one
# classify_lines.py writes.
CLASSIFY_COLUMNS = [
    'wn_obs', 'unc_wn_obs', 'obs_intens', 'char', 'low_id', 'upp_id',
    'calc_intens', 'orig_calc_intens', 'u_calc', 'imputed', 'intens_from_f',
    'intens_to_f', 'dif_wn_O-C', 'grade', 'notes1', 'notes2', 'new',
    'accepted', 'n_accepted', 'low_E', 'upp_E', 'rwn', 'BF',
]

WN_TOL = 0.05   # cm^-1: LOPT rounds the wavenumber to the precision of its
                # uncertainty, so a printed value may differ from the exact
                # one in the line list by half of the last printed digit.


def _num(text: str) -> float:
    """One numeric cell of a LOPT output file; '_' and blank mean 'not given'."""
    text = (text or '').strip()
    if text in ('', '_'):
        return np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


def read_observed_lines() -> pd.DataFrame:
    """The whole observed line list: wn_obs, unc_wn_obs, obs_intens, char.

    This is every measured line, including the ones no candidate transition
    ever matched.  The masking check needs all of them: a line that hides a
    predicted transition need not itself be identified.
    """
    wb = openpyxl.load_workbook(cl.LINES_FILE, read_only=True, data_only=True)
    ws = wb[cl.CFG.lines.sheet]
    col = cl.column_index(ws, cl.CFG.lines, cl.LINES_FILE)
    recs = []
    for row in ws.iter_rows(min_row=2):
        wn = row[col['wn']].value
        if wn is None:
            continue
        unc = row[col['u_wn']].value
        inten = row[col['intensity']].value
        char = row[col['character']].value
        recs.append((float(wn),
                     float(unc) if unc not in (None, '') else np.nan,
                     float(inten) if inten not in (None, '') else np.nan,
                     '' if char is None else str(char).strip()))
    wb.close()
    return pd.DataFrame(recs, columns=['wn_obs', 'unc_wn_obs',
                                       'obs_intens', 'char'])


def read_legacy_pairs() -> dict:
    """{(lower_id, upper_id) sorted: [wavenumbers]} of legacy identifications.

    These are the id1/id2 columns of the observed-line workbook: the
    identifications published before this work, against which an accepted
    identification counts as "old" rather than "new".
    """
    wb = openpyxl.load_workbook(cl.LINES_FILE, read_only=True, data_only=True)
    ws = wb[cl.CFG.lines.sheet]
    col = cl.column_index(ws, cl.CFG.lines, cl.LINES_FILE)
    legacy = {}
    for row in ws.iter_rows(min_row=2):
        wn = row[col['wn']].value
        if wn is None:
            continue
        i1 = cl.to_str_id(row[col['id1']].value)
        i2 = cl.to_str_id(row[col['id2']].value)
        if i1 and i2:
            legacy.setdefault(tuple(sorted((i1, i2))), []).append(float(wn))
    wb.close()
    return legacy


def read_lopt_lines(path: str) -> pd.DataFrame:
    """LOPT's line-output file as a table, keeping the columns needed here.

    The file is tab-separated with one header row; which columns it holds
    depends on the print options of the LOPT parameter file, so they are looked
    up by name.  Numbers LOPT does not give are written by it as '_'.
    """
    with open(path, newline='', encoding='utf-8', errors='replace') as f:
        rd = csv.reader(f, delimiter='\t')
        header = next(rd)
        ix = {h.strip(): i for i, h in enumerate(header)}
        for need in ('wn_o', 'Iobs', 'L1', 'L2', 'E1', 'E2', 'F', 'Weight'):
            if need not in ix:
                raise ValueError(
                    f"{os.path.basename(path)}: column {need!r} is missing - "
                    f"switch its printing on in the LOPT parameter file")
        recs = []
        for r in rd:
            if not r or len(r) <= ix['Weight']:
                continue
            l1 = r[ix['L1']].strip()
            l2 = r[ix['L2']].strip()
            if not l1 or not l2:
                continue
            recs.append((_num(r[ix['wn_o']]), _num(r[ix['Iobs']]), l1, l2,
                         _num(r[ix['E1']]), _num(r[ix['E2']]),
                         r[ix['F']].strip(), _num(r[ix['Weight']])))
    return pd.DataFrame(recs, columns=['wn_lopt', 'I_lopt', 'low_id', 'upp_id',
                                       'low_E', 'upp_E', 'flag', 'weight'])


def read_run(path: str, energies: str = 'refit', levels_path: str = None,
             e_input: dict = None, verbose: bool = True) -> pd.DataFrame:
    """A LOPT line output as the table classify_lines.py would have written.

    One row per candidate transition of the LOPT run, with accepted = 1 where
    the line entered the fit, new = 1 where the identification is not the
    legacy one, and low_E / upp_E the optimized level energies.

    energies chooses where those level energies come from:
      'refit'  (default) re-solve them from the accepted lines with the
               pipeline's own least squares - refit_energies() explains why
               this is what the validation report needs;
      'lopt'   LOPT's own optimized energies, taken from levels_path when a
               level-output file is given (four decimals for every level) and
               otherwise from the E1/E2 columns of the line output (rounded
               to each level's own precision).
    e_input is passed on to refit_energies() as the anchor energies.
    """
    if energies not in ('refit', 'lopt'):
        raise ValueError(f"energies must be 'refit' or 'lopt', not {energies!r}")
    lo = read_lopt_lines(path)
    obs = read_observed_lines()
    legacy = read_legacy_pairs()

    # LOPT prints a wavenumber rounded to the precision of its uncertainty, so
    # the printed value must be matched to the exact one by nearest neighbour
    # within WN_TOL, not by equality.
    obs = obs.sort_values('wn_obs').reset_index(drop=True)
    wn_sorted = obs['wn_obs'].to_numpy()
    obs_rec = list(obs.itertuples(index=False, name=None))

    def nearest(wn):
        j = int(np.searchsorted(wn_sorted, wn))
        best = None
        for k in (j - 1, j):
            if 0 <= k < len(wn_sorted):
                d = abs(wn_sorted[k] - wn)
                if d <= WN_TOL and (best is None or d < best[0]):
                    best = (d, obs_rec[k])
        return None if best is None else best[1]

    rows = []
    n_unmatched = 0
    for r in lo.itertuples(index=False):
        pair = tuple(sorted((r.low_id, r.upp_id)))
        hit = nearest(r.wn_lopt)
        if hit is None:
            n_unmatched += 1
            wn, unc, inten, char = r.wn_lopt, np.nan, r.I_lopt, ''
        else:
            wn, unc, inten, char = hit
        is_new = 0 if any(abs(x - r.wn_lopt) <= WN_TOL
                          for x in legacy.get(pair, ())) else 1
        rows.append({
            'wn_obs': wn, 'unc_wn_obs': unc, 'obs_intens': inten, 'char': char,
            'low_id': r.low_id, 'upp_id': r.upp_id,
            'new': is_new, 'accepted': 1 if r.weight > 0 else 0,
            'low_E': r.low_E, 'upp_E': r.upp_E, 'BF': r.weight,
            'notes1': r.flag,
        })
    df = pd.DataFrame(rows)
    for c in CLASSIFY_COLUMNS:
        if c not in df.columns:
            df[c] = np.nan
    df['n_accepted'] = df.groupby('wn_obs')['accepted'].transform('sum')
    if verbose:
        n_acc = int(df['accepted'].sum())
        n_new = int(((df['accepted'] == 1) & (df['new'] == 1)).sum())
        print(f"  {os.path.basename(path)}: {len(df)} candidate rows, "
              f"{n_acc} accepted ({n_new} of them new), "
              f"{df['wn_obs'].nunique()} distinct observed lines")
        if n_unmatched:
            print(f"  note: {n_unmatched} LOPT rows had no line of "
                  f"{os.path.basename(cl.LINES_FILE)} at their wavenumber; "
                  f"LOPT's own rounded values were used for them")

    if energies == 'refit':
        df = apply_energies(df, refit_energies(df, e_input=e_input,
                                               verbose=verbose))
    elif levels_path:
        lev = read_levels(levels_path)
        df = apply_energies(df, dict(zip(lev['level_id'], lev['E_lopt'])))
        if verbose:
            print(f"  level energies read from "
                  f"{os.path.basename(levels_path)} ({len(lev)} levels)")
    return df[CLASSIFY_COLUMNS]


# ----------------------------------------------------------------------
# LOPT's level-output file
# ----------------------------------------------------------------------
def read_levels(path: str) -> pd.DataFrame:
    """LOPT's level-output file: level_id, E_lopt, u_lopt, n_lines_lopt.

    Why this file is worth reading even though the line output also carries
    the energies.  In the line output LOPT prints E1 and E2 rounded to the
    precision of the level's own uncertainty - four decimals for a well
    determined level, fewer for a poorly determined one - whereas the level
    output prints four decimals throughout, and the energy shift
    dE = E_final - E_input of a well determined level is itself of order
    0.001 cm^-1.  For LOPT_output_lines_revised.txt the two sources happen to
    agree exactly on all 593 levels, so nothing is lost by omitting this file;
    it simply removes a dependence on how much precision LOPT chose to print.

    Columns of the file: Designation, Energy, D1, D2stat, D2sys, D2tot,
    N_lines, Comments.  D2tot is the total uncertainty of the level relative
    to the ground level.  N_lines is LOPT's own line count, written as
    "total, nL, nD", so only the leading number is kept.
    """
    lv = pd.read_csv(path, sep='\t', dtype=str).rename(columns=str.strip)
    for need in ('Designation', 'Energy'):
        if need not in lv.columns:
            raise ValueError(f"{os.path.basename(path)}: column {need!r} is "
                             f"missing - this is not a LOPT level output file")
    out = pd.DataFrame({
        'level_id': lv['Designation'].str.strip(),
        'E_lopt': lv['Energy'].map(_num),
        'u_lopt': (lv['D2tot'].map(_num) if 'D2tot' in lv.columns
                   else pd.Series(np.nan, index=lv.index)),
        'n_lines_lopt': (lv['N_lines'].str.split(',').str[0].map(_num)
                         if 'N_lines' in lv.columns
                         else pd.Series(np.nan, index=lv.index)),
    })
    return out[out['level_id'].astype(bool) & out['E_lopt'].notna()]


def adopted_energies() -> dict:
    """{level_id: E_adopt} of the adopted-level workbook."""
    wb = openpyxl.load_workbook(cl.LEVELS_FILE, read_only=True, data_only=True)
    ws = wb[cl.CFG.levels.sheet]
    col = cl.column_index(ws, cl.CFG.levels, cl.LEVELS_FILE)
    out = {}
    for row in ws.iter_rows(min_row=2):
        lid = cl.to_str_id(row[col['id']].value)
        e = row[col['E']].value
        if lid and e is not None:
            out[lid] = float(e)
    wb.close()
    return out


# ----------------------------------------------------------------------
# Re-optimization in the pipeline's own metric
# ----------------------------------------------------------------------
def refit_energies(df: pd.DataFrame, e_input: dict = None,
                   verbose: bool = True) -> dict:
    """Re-solve the level energies from an accepted-line set, exactly the way
    classify_lines.optimize_levels() does.

    Why a second optimizer is needed at all.  The validation report compares
    the energy shift dE = E_final - E_input of each tested level against the
    same quantity measured on levels that are false by construction: the
    decoys of decoy_mc.py and the chance levels of chance_mc.py.  Those
    reference populations are produced by re-running the pipeline, so their
    E_final comes from classify_lines.optimize_levels().  If the tested
    levels' E_final came from LOPT instead, the two sides of the comparison
    would be measured with different instruments - LOPT tunes single-line
    levels, treats blends by its centroid model and rounds its output, and
    its energies differ from the pipeline's by up to 0.2 cm^-1, the same
    order as the dE being tested.  Re-solving the REVISED identifications
    with the pipeline's own least squares removes that mismatch, so the only
    difference between a revised report and the baseline one is the thing
    being studied: which lines are assigned to which levels.

    The model, identical to optimize_levels().  Every accepted transition is
    one observation equation E(upper) - E(lower) = wn_obs, weighted by
    w = BF^2 / u_wn^2.  BF is the branching fraction of the line: 1 for an
    unblended line, and for a component of a blend the share of the measured
    intensity it carries.  It enters squared so that the components of one
    blend together count as a single observation, which is LOPT's "centroid"
    model.  u_wn is the measurement uncertainty of the wavenumber.  The
    energies minimizing the weighted sum of squared residuals solve the
    normal equations (A^T W A) x = A^T W y, formed here directly.

    Levels sitting at 0.0 in the adopted-level list are held fixed (the
    ground level).  Levels no accepted transition touches are absent from the
    result and keep their input energy.  A group of levels tied only to each
    other, with no chain of transitions reaching a fixed level, is determined
    only up to a common offset; one level of such a group is then held at its
    input energy so that the group keeps the position it has.

    e_input: {level_id: energy} to anchor on; the adopted-level workbook by
    default.  Pass the revised energies when a level has been moved, so that
    a floating group is anchored where the revision put it.

    Checked: fed the accepted set of LOPT_output_lines.txt - the same
    identifications the pipeline itself made - it reproduces the low_E/upp_E
    columns of line_classifications.csv to a median of 1.3e-6 cm^-1 and a
    maximum of 0.0012 cm^-1, the residue being LOPT's rounding of the
    branching fractions it prints to four decimals.

    Returns {level_id: energy} for the levels the accepted lines determine.
    """
    acc = df[(df['accepted'] == 1) & df['low_id'].notna() &
             df['upp_id'].notna()].copy()
    if not len(acc):
        return {}
    bf = pd.to_numeric(acc['BF'], errors='coerce').fillna(1.0).to_numpy()
    u = pd.to_numeric(acc['unc_wn_obs'], errors='coerce').to_numpy()
    y = pd.to_numeric(acc['wn_obs'], errors='coerce').to_numpy()
    ok = np.isfinite(bf) & np.isfinite(u) & (u > 0) & np.isfinite(y)
    if not ok.all():
        acc, bf, u, y = acc[ok], bf[ok], u[ok], y[ok]
    w = bf ** 2 / u ** 2

    ids = sorted(set(acc['low_id']) | set(acc['upp_id']))
    ix = {k: i for i, k in enumerate(ids)}
    lo = acc['low_id'].map(ix).to_numpy()
    up = acc['upp_id'].map(ix).to_numpy()
    n = len(ids)

    e_in = dict(adopted_energies() if e_input is None else e_input)
    x = np.array([e_in.get(k, 0.0) for k in ids], dtype=float)

    # Anchors: the ground level, plus one level of every group of levels not
    # connected to it through the accepted transitions (union-find).
    anchored = {ix[k] for k in ids if abs(e_in.get(k, np.nan)) < 1e-12}
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for i, j in zip(lo, up):
        ri, rj = find(int(i)), find(int(j))
        if ri != rj:
            parent[ri] = rj
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    for members in groups.values():
        if not any(m in anchored for m in members):
            anchored.add(min(members))

    normal = np.zeros((n, n))
    rhs = np.zeros(n)
    np.add.at(normal, (up, up), w)
    np.add.at(normal, (lo, lo), w)
    np.add.at(normal, (up, lo), -w)
    np.add.at(normal, (lo, up), -w)
    np.add.at(rhs, up, w * y)
    np.add.at(rhs, lo, -w * y)

    anc = sorted(anchored)
    free = [i for i in range(n) if i not in anchored]
    if free:
        # the anchored energies are known, so their terms move to the
        # right-hand side and only the free levels are solved for
        b = rhs[free] - normal[np.ix_(free, anc)] @ x[anc]
        x[free] = np.linalg.solve(normal[np.ix_(free, free)], b)
    if verbose:
        print(f"  re-optimized {len(free)} of {n} level energies from "
              f"{len(acc)} accepted lines with the pipeline's own least "
              f"squares ({len(anc)} level(s) anchored)")
    return dict(zip(ids, x))


def apply_energies(df: pd.DataFrame, energies: dict) -> pd.DataFrame:
    """Overwrite low_E / upp_E with the energies of the given mapping.

    Rows whose level is not in the mapping keep the energy already there.
    """
    out = df.copy()
    for col, ecol in (('low_id', 'low_E'), ('upp_id', 'upp_E')):
        mapped = out[col].map(energies)
        out[ecol] = mapped.where(mapped.notna(), out[ecol])
    return out
