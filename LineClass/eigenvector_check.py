"""Which Cowan eigenvectors the observations confirm: a recursive check by
hyperfine constants and line intensities (the user's idea, 2026-10-07).

WHAT IT IS FOR
==============
A calculated A (``hfs_A_theory.py``) is only as good as the level's
eigenvector, the mixture of basis states the Cowan fit gives it.  Two kinds
of observation test an eigenvector:

a) the level's measured A agrees with the calculated one;
b) the intensities of the level's lines agree with Icalc.

Icalc = C * gA * (rwn/1e8) * exp(-E_up/kT) contains the population of the
upper level through exp(-E_up/kT).  Within one upper level's branch, the set
of its lines, the population is a common factor, so the *ratios* inside a
branch do not depend on the population model.  They depend on the gA alone,
that is on the eigenvectors of the upper level and of each lower one.  If
the lower levels are already known to be good, a branch that keeps its
ratios tests the upper level; if the upper level is good, the line to a
level X, set against the rest of the branch, tests X.

THE RECURSION
=============
Seeds: every level whose measured A agrees with A_calc within Z_SEED, the
`z` of hfs_A_theory.csv (measured uncertainty and the file's precision).

Each round, every level not yet good is tested on two sets of lines:

* down: its own branch to good lower levels, with one free constant (its
  population); needs 2 lines and gives n - 1 degrees of freedom;
* up: its line from each good upper level whose branch has at least 2 other
  lines to good levels.  The branch's constant comes from those other lines,
  and the residual is the line's distance from it.

The residual of a line is ln(I_obs/I_calc), with the uncertainty
sqrt(u_calc^2 + s_obs^2).  u_calc is Cowan's gA uncertainty, as the
classification table carries it.  s_obs is the observational scatter: the
robust standard deviation of ln(I_obs/I_calc) about each upper level's
median, over branches of at least 4 lines, less the median u_calc in
quadrature.  It is measured on every run.  A level with at least 2 degrees
of freedom and chi^2 probability >= P_GOOD becomes good, and the next round
can use it.  The rounds stop when none is added.

Only accepted lines with one classification and both intensities are used.
A blended line's intensity belongs to no single transition.

THE VERDICT
===========
After the last round, every level is tested once more against the final
good set, itself left out, and gets a status:

    good_A       a seed
    good_I       made good in round `round`
    fail         chi^2/n > VETO (the intensities contradict the eigenvector)
    doubtful     not good, chi^2/n <= VETO: typically a level with many
                 lines, where the noise model rather than the eigenvector
                 is the likelier fault
    untested     fewer than 2 degrees of freedom

VETO = 2.4 is the 90th percentile of chi^2/n over the seeds on 2026-10-07;
each run reports the current percentile.  The check is coarse: one line's
intensity scatters by a factor of about 2, so it finds gross errors of a
composition, not amplitude errors of 0.05-0.10.  Passing it vouches for
nothing finer.  ``hfs_A_candidates.py`` uses only the veto: a level that
fails goes to tier 3.

USAGE
=====
    python eigenvector_check.py            writes eigenvector_check.csv/.log
    python eigenvector_check.py --no-write
"""

import argparse
import collections
import csv
import math
import os

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
THEORY_CSV = os.path.join(HERE, 'hfs_A_theory.csv')
LINES = os.path.join(HERE, 'iter_hfs', 'line_classifications.csv')
ID_MAP = os.path.join(HERE, 'iter_hfs', 'IDEN2', 'IDEN_level_ids.txt')
OUT_CSV = os.path.join(HERE, 'eigenvector_check.csv')
OUT_LOG = os.path.join(HERE, 'eigenvector_check.log')

Z_SEED = 2.0
P_GOOD = 0.01
VETO = 2.4
MIN_BRANCH = 4             # lines of a branch used for the scatter

FIELDS = ('level_id', 'iden2_row', 'J', 'leading', 'status', 'round',
          'z_A', 'n_down', 'n_up', 'dof', 'chi2_n', 'p', 'mean_pull_up')


def read_lines(path):
    """The usable lines: accepted, one classification, both intensities;
    each with x = ln(I_obs/I_calc) and its u_calc."""
    out = []
    with open(path, encoding='utf-8', newline='') as fh:
        for r in csv.DictReader(fh):
            if r['accepted'] not in ('1', '1.0') or \
                    r['n_accepted'] not in ('1', '1.0'):
                continue
            io, ic = float(r['obs_intens'] or 0), float(r['calc_intens'] or 0)
            if io <= 0 or ic <= 0:
                continue
            out.append({'low': r['low_id'], 'upp': r['upp_id'],
                        'x': math.log(io / ic),
                        'u_calc': float(r['u_calc'] or 0)})
    return out


def observational_scatter(lines, min_branch=MIN_BRANCH):
    """s_obs: the scatter about each branch's median, less u_calc."""
    by_up = collections.defaultdict(list)
    for ln in lines:
        by_up[ln['upp']].append(ln)
    dev, uc = [], []
    for branch in by_up.values():
        if len(branch) < min_branch:
            continue
        x = np.array([b['x'] for b in branch])
        dev.extend(x - np.median(x))
        uc.extend(b['u_calc'] for b in branch)
    sd = 1.4826 * float(np.median(np.abs(dev)))
    u = float(np.median(uc))
    return math.sqrt(max(sd * sd - u * u, 0.01)), sd, u


class Check:
    """The recursion on one set of lines."""

    def __init__(self, lines, s_obs):
        for ln in lines:
            ln['s'] = math.hypot(ln['u_calc'], s_obs)
        self.lines = lines
        self.by_up = collections.defaultdict(list)
        self.by_low = collections.defaultdict(list)
        for ln in lines:
            self.by_up[ln['upp']].append(ln)
            self.by_low[ln['low']].append(ln)

    def branch_constant(self, upp, good, exclude):
        """(c, u_c) of a good upper level's branch to good levels, the line
        to `exclude` left out; None below 2 lines."""
        sel = [ln for ln in self.by_up[upp]
               if ln['low'] in good and ln['low'] != exclude]
        if len(sel) < 2:
            return None
        w = np.array([1.0 / ln['s'] ** 2 for ln in sel])
        x = np.array([ln['x'] for ln in sel])
        return float(np.sum(w * x) / np.sum(w)), float(1.0 / math.sqrt(
            np.sum(w)))

    def test(self, lid, good):
        """(n_down, n_up, dof, chi2, mean_pull_up) of `lid` against `good`
        (without `lid` itself)."""
        pulls, dof = [], 0
        down = [ln for ln in self.by_up.get(lid, [])
                if ln['low'] in good and ln['low'] != lid]
        if len(down) >= 2:
            w = np.array([1.0 / ln['s'] ** 2 for ln in down])
            x = np.array([ln['x'] for ln in down])
            c = np.sum(w * x) / np.sum(w)
            pulls.extend((x - c) / np.array([ln['s'] for ln in down]))
            dof += len(down) - 1
        up = []
        for ln in self.by_low.get(lid, []):
            if ln['upp'] not in good or ln['upp'] == lid:
                continue
            bc = self.branch_constant(ln['upp'], good, lid)
            if bc is None:
                continue
            up.append((ln['x'] - bc[0]) / math.hypot(ln['s'], bc[1]))
        pulls.extend(up)
        dof += len(up)
        chi2 = float(np.sum(np.square(pulls))) if pulls else 0.0
        return (len(down) if len(down) >= 2 else 0, len(up), dof, chi2,
                float(np.mean(up)) if up else None)

    def run(self, seeds, levels):
        """{level_id: round} of the good levels (0 for the seeds)."""
        good = {lid: 0 for lid in seeds}
        rnd = 0
        while True:
            rnd += 1
            new = []
            for lid in levels:
                if lid in good:
                    continue
                _, _, dof, chi2, _ = self.test(lid, good)
                if dof >= 2 and stats.chi2.sf(chi2, dof) >= P_GOOD:
                    new.append(lid)
            if not new:
                return good
            for lid in new:
                good[lid] = rnd


def status_of(is_good, rnd, dof, chi2_n):
    if is_good:
        return 'good_A' if rnd == 0 else 'good_I'
    if dof < 2:
        return 'untested'
    return 'fail' if chi2_n > VETO else 'doubtful'


def read_status(path=OUT_CSV):
    """{level_id: row} of eigenvector_check.csv; {} if there is none."""
    if not os.path.isfile(path):
        return {}
    with open(path, encoding='utf-8', newline='') as fh:
        return {r['level_id']: r for r in csv.DictReader(fh)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--lines', default=LINES,
                    help='line table of the working set')
    ap.add_argument('--no-write', action='store_true')
    args = ap.parse_args(argv)

    import cowan_gA

    with open(THEORY_CSV, encoding='utf-8', newline='') as fh:
        theory = {r['level_id']: r for r in csv.DictReader(fh)
                  if r['level_id']}
    iden = {lid: row for row, lid in cowan_gA.read_id_map(ID_MAP).items()}
    seeds = {lid for lid, r in theory.items()
             if r['z'] and abs(float(r['z'])) <= Z_SEED}

    lines = read_lines(args.lines)
    s_obs, sd, u_med = observational_scatter(lines)
    chk = Check(lines, s_obs)
    levels = sorted(set(chk.by_up) | set(chk.by_low))
    good = chk.run(seeds, levels)

    rows = []
    for lid in sorted(set(levels) | seeds):
        n_down, n_up, dof, chi2, mean_up = chk.test(lid, good)
        chi2_n = chi2 / dof if dof else None
        t = theory.get(lid, {})
        rows.append({
            'level_id': lid,
            'iden2_row': iden.get(lid, ''),
            'J': t.get('J', ''),
            'leading': t.get('leading', ''),
            'status': status_of(lid in good, good.get(lid), dof,
                                chi2_n or 0.0),
            'round': good.get(lid, ''),
            'z_A': t.get('z', ''),
            'n_down': n_down,
            'n_up': n_up,
            'dof': dof,
            'chi2_n': '' if chi2_n is None else '%.2f' % chi2_n,
            'p': '' if not dof else '%.3g' % stats.chi2.sf(chi2, dof),
            'mean_pull_up': '' if mean_up is None else '%+.2f' % mean_up,
        })

    seed_chi = [float(r['chi2_n']) for r in rows
                if r['status'] == 'good_A' and r['chi2_n']]
    count = collections.Counter(r['status'] for r in rows)
    rounds = collections.Counter(v for v in good.values() if v)
    report = [
        'eigenvector_check.py',
        '',
        'lines used (accepted, single, both intensities): %d' % len(lines),
        'scatter of ln(Iobs/Icalc) about the branch median %.2f, median '
        'u_calc %.2f: s_obs %.2f' % (sd, u_med, s_obs),
        'seeds (|z| of A <= %.0f): %d' % (Z_SEED, len(seeds)),
        'rounds: %s' % ', '.join('%d: +%d' % kv
                                 for kv in sorted(rounds.items())),
        'chi2/n of the seeds: median %.2f, 90th percentile %.2f (VETO %.1f)'
        % (np.median(seed_chi), np.percentile(seed_chi, 90), VETO),
        'status: %s' % ', '.join('%s %d' % kv for kv in sorted(
            count.items())),
        '',
        'failing levels (chi2/n > %.1f):' % VETO,
    ]
    for r in sorted((r for r in rows if r['status'] == 'fail'),
                    key=lambda r: -float(r['chi2_n'])):
        report.append('  %s %5s  %-24s dof %3d  chi2/n %5s  z_A %s'
                      % (r['level_id'], r['iden2_row'], r['leading'],
                         r['dof'], r['chi2_n'], r['z_A'] or '-'))
    print('\n'.join(report))

    if not args.no_write:
        import output_files
        output_files.require_writable([OUT_CSV, OUT_LOG])
        with open(OUT_CSV, 'w', encoding='utf-8', newline='') as fh:
            wr = csv.DictWriter(fh, fieldnames=FIELDS, lineterminator='\n')
            wr.writeheader()
            wr.writerows(rows)
        with open(OUT_LOG, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write('\n'.join(report) + '\n')
        print('wrote %s and %s' % (os.path.basename(OUT_CSV),
                                   os.path.basename(OUT_LOG)))
    return rows


if __name__ == '__main__':
    main()
