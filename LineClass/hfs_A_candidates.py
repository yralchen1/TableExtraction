"""Calculated A constants proposed for A_hfs_levels.csv: the review file for
the user, and the step that writes the adopted rows (2026-10-07).

WHAT IT DOES
============
Over 200 identified levels have no usable A in ``A_hfs_levels.csv``: no row, or a
row "not determined".  The pipeline gives such a level S = I*A*J = 0 at its
head, so the energy LOPT finds for it stands about kappa*S from its center of
gravity, kappa being the fraction of the way from the center of gravity to
the head at which its lines are measured.  ``hfs_A_theory.py`` calculates
an A for most of them from the Cowan eigenvectors, with a total uncertainty
``u_total``.  This program sorts the levels into three tiers, as agreed with
the user on 2026-10-07, and writes the evidence for each one to
``hfs_A_candidates.csv``:

* tier 1: ``cancel`` < CANCEL_MAX, u_total <= REL_MAX * |A|, and the
  leading configuration is tested by measured constants (4f2 6s, 4f5d2,
  4f3);
* tier 2: the same, but the leading configuration is not tested (4f2 5d,
  4f5d6s, 4f5d6p, ...);
* tier 3: A is a cancellation residue or too uncertain, or the eigenvector
  is suspect: the level lies more than E_DEV_MAX from its calculated
  energy, or its line intensities fail ``eigenvector_check.py`` (status
  "fail").

Until 2026-10-08 a tier-3 level got a "not determined" row, with
u_A = sqrt(A^2 + u_total^2 + (0.1 * sum|theta*a|)^2) (``undetermined_u``),
and a calculated row whose level fell to tier 3 was withdrawn to one.  Since
then (user, 2026-10-08) a tier-3 level is adopted with its calculated A as
a *test* row instead - whether it has no row, a "not determined" one, or a
calculated one.  A "not determined" row leaves S = 0 and so tests nothing;
a calculated A moves the level's lines, and IDEN2 and LOPT show whether
they still fit, so that discordant assignments can be dropped.  A test
row's u_A is ``test_u`` = sqrt(u_total^2 + (0.1 * sum|theta*a|)^2), and its
source ends in TEST_TAG.  A measured or semiempirical row is never replaced
by a tier-3 value, and a resolved level never takes one.  A level of
``discarded_levels.csv`` (files.discarded_levels) is left out altogether: its
position has been given up.

Tiers 1 and 2 are adopted (user, 2026-10-07) where the table has no row or
a "not determined" one, and in place of a measured A when the calculated one
is much more accurate and agrees with it (``hfs_A_fit.much_more_accurate``:
an uncertainty at most a third of the measured one, a difference within 3
combined standard uncertainties).  A measured A that disagrees is left to
the user, except where the user has decided (USER_DECISIONS).  The
semiempirical "composition" rows stay: they are the priors of
``hfs_A_fit.py``.  A level of ``[hfs] resolved_levels`` keeps its measured A.

THE COLUMNS
===========
level_id, iden2_row, J, leading     the level (iden2_row: the working set's
                                    IDEN2 enlev.dat row)
table_now           the level's row of A_hfs_levels.csv: "absent", or its
                    source; A_table, u_table its value
tier                1, 2 or 3
A_calc, u_total     from hfs_A_theory.csv, cm^-1; also u_amp, u_par, u_cfg
rel_u               u_total / |A_calc|
cancel              sum |theta*a| / |A|: 1 when nothing cancels
z_table             (A_calc - A_table) / combined uncertainty
eigenvector         eigenvector_check.csv: status (chi^2/n, dof)
E_dev               E_obs - E_calc, cm^-1; E_obs from the working set
S, u_S              I*A*J and I*J*u_total: the head's distance from the
                    center of gravity
n_lines             the level's accepted lines in the working set
kappa_mean          their kappa, weighted by BF / unc^2 as LOPT weighs them
dE_expected         -kappa_mean * S: how far the level's energy should move
                    when the calculated A is adopted and the chain is rerun
n_flagged           Sugar's *r / *v lines of the level, and how many of them
flag_disagree       the calculated A contradicts: *r means D > 0, *v D < 0,
                    D = S_upper - S_lower; the partner's A from the table,
                    else from the calculation
wide_unflagged      unflagged lines inside the range of Sugar's flags, as
                    intense as the weakest tenth of his flagged ones or more,
                    for which D, so calculated, is at least WIDE_D
A_components        A from Sugar's printed component positions of the
                    level's lines (the patterns of hfs_A_fit.py, SET_ASIDE
                    left out), every partner held at its table A, else its
                    calculated one; the partners' uncertainties included
action              what --write does with the row
proposed_*          the row for A_hfs_levels.csv
notes               questionable level, resolved level, HOLD, ...

USAGE
=====
    python hfs_A_candidates.py          writes hfs_A_candidates.csv
    python hfs_A_candidates.py --no-write
    python hfs_A_candidates.py --write  also writes the adopted rows into
                                        A_hfs_levels.csv; a measured row
                                        replaced goes to
                                        A_hfs_levels_superseded.csv
"""

import argparse
import csv
import datetime
import math
import os
import tomllib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
THEORY_CSV = os.path.join(HERE, 'hfs_A_theory.csv')
A_LEVELS = os.path.join(HERE, 'A_hfs_levels.csv')
LINES = os.path.join(HERE, 'iter_hfs', 'line_classifications.csv')
ID_MAP = os.path.join(HERE, 'iter_hfs', 'IDEN2', 'IDEN_level_ids.txt')
QUESTIONABLE = os.path.join(HERE, 'questionable_levels.txt')
CONFIG = os.path.join(HERE, 'lineclass_config.toml')
OUT_CSV = os.path.join(HERE, 'hfs_A_candidates.csv')

I_NUC = 2.5
CANCEL_MAX = 2.0
REL_MAX = 0.25
E_DEV_MAX = 300.0          # cm^-1: a level further from its calculated
                           # energy has a suspect eigenvector (user)
WIDE_D = 0.3               # cm^-1
W_NOTE = 0.01              # weight without parameters worth a note
DATE = datetime.date.today().isoformat()   # stamped on the rows written
DISCARDED = os.path.join(HERE, 'discarded_levels.csv')
TEST_TAG = 'a test of the eigenvector'     # source suffix of tier-3 rows
                                           # adopted as a test (2026-10-08)

#: levels whose measured A the user has dismissed for the calculated one.
USER_DECISIONS = {
    '059003.000151': 'dismissed by the user 2026-10-07: a 4f2 5d level with '
                     'no other sign of a wide pattern; the one printed '
                     'companion (29150.549 of 29151.147 *r) is taken to be '
                     'another line (hfs_A_fit.SET_ASIDE)',
    '059003.000198': 'dismissed by the user 2026-10-07: the center of gravity '
                     'of 21793.442 (1116-1081), a line too weak to be '
                     'trusted',
}

SUPERSEDED_FIELDS = ('level_id', 'cfg', 'J', 'A_cm-1', 'u_A', 'source',
                     'n_flagged', 'superseded_on', 'reason')


def _f(text):
    return float(text) if text not in ('', None) else None


def tier_of(A, u_total, cancel, tested, E_dev=None, eig_fail=False):
    if A is None or cancel is None or cancel >= CANCEL_MAX \
            or u_total > REL_MAX * abs(A) \
            or (E_dev is not None and abs(E_dev) > E_DEV_MAX) or eig_fail:
        return 3
    return 1 if tested else 2


def test_u(A, u_total, cancel):
    """u_A of a tier-3 level adopted as a test: its calculated uncertainty
    and a 10 per cent error in any one of the terms that make up A."""
    terms = (cancel or 1.0) * abs(A or 0.0)
    return math.hypot(u_total, 0.1 * terms)


def is_test_row(source):
    """True for a calculated row adopted from tier 3 as a test."""
    return bool(source) and source.endswith(TEST_TAG)


def read_discarded(path=DISCARDED):
    """The level_ids of discarded_levels.csv (none if there is no file)."""
    import classify_lines
    if not os.path.isfile(path):
        return set()
    return {r['level_id'] for r in classify_lines.read_discarded_records(path)}


def undetermined_u(A, u_total, cancel):
    """u_A of a tier-3 level's "not determined" row."""
    terms = (cancel or 0.0) * abs(A or 0.0)
    return math.sqrt((A or 0.0) ** 2 + u_total ** 2 + (0.1 * terms) ** 2)


def kind_of(source):
    """What a table row is: absent, undetermined, calculated, semiempirical
    or measured."""
    import hfs_A_theory
    import hfs_correction
    import hfs_kappa
    if source is None:
        return 'absent'
    if source.startswith(hfs_A_theory.CALC_SOURCE):
        return 'calculated'
    if not hfs_correction.is_determined(source):
        return 'undetermined'
    if hfs_kappa.is_semiempirical(source):
        return 'semiempirical'
    return 'measured'


def action_of(tier, kind, better, agrees, decided, resolved, test=False):
    """What --write does: 'add', 'replace' (either possibly 'as a test:
    tier 3'), or a reason for leaving the row alone.  `test`: the table's
    row is a calculated one already adopted from tier 3 as a test
    (is_test_row)."""
    if kind == 'calculated':
        return ('in the table' if tier < 3 or test
                else 'replace as a test: tier 3')
    if decided:
        return 'replace (user decision)'
    if resolved:
        return 'keep: resolved level, its A is measured from the sublevels'
    if kind == 'absent':
        return 'add' if tier < 3 else 'add as a test: tier 3'
    if kind == 'undetermined':
        return 'replace' if tier < 3 else 'replace as a test: tier 3'
    if kind == 'semiempirical':
        return 'keep: semiempirical prior of hfs_A_fit'
    if tier == 3:
        return 'keep: tier 3'
    if not agrees:
        return 'user decides: disagrees with the measured A'
    return 'replace' if better else 'keep: not much more accurate'


def components_value(patterns, lid, partner):
    """(A, u, n_positions) of `lid` from the printed component positions of
    its patterns, each partner at `partner(id) -> (A, u)`; None if no
    pattern reaches the level or a partner has no A."""
    import hfs_A_fit
    G, R, W, H = [], [], [], {}
    for p in patterns:
        if lid not in (p.low_id, p.upp_id):
            continue
        other = p.upp_id if p.low_id == lid else p.low_id
        pa = partner(other)
        if pa is None:
            return None
        f1, f2 = I_NUC + p.J1, I_NUC + p.J2
        w = 1.0 / hfs_A_fit.sigma_of(p.wn) ** 2
        for k, d in enumerate(p.dwn, 1):
            c_low = 0.5 * k * (2 * f1 + 1 - k)
            c_upp = -0.5 * k * (2 * f2 + 1 - k)
            g, h = (c_low, c_upp) if p.low_id == lid else (c_upp, c_low)
            G.append(g)
            R.append(d - h * pa[0])
            W.append(w)
            H.setdefault(other, (pa[1], []))[1].append((len(G) - 1, h))
    if not G:
        return None
    G, R, W = np.array(G), np.array(R), np.array(W)
    gg = float(np.sum(W * G * G))
    A = float(np.sum(W * G * R)) / gg
    var = 1.0 / gg
    for u_p, terms in H.values():
        dA = -sum(W[i] * G[i] * h for i, h in terms) / gg
        var += (dA * u_p) ** 2
    return A, math.sqrt(var), len(G)


def read_csv(path, delimiter=','):
    with open(path, encoding='utf-8', newline='') as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def write_table(out, path=A_LEVELS, date=DATE):
    """Write the rows `out` adopts into the A table; a measured row replaced
    is appended to hfs_A_theory.SUPERSEDED_NAME beside it.  Returns
    (added, replaced, superseded) counts."""
    import hfs_A_theory
    import output_files
    old_path = os.path.join(os.path.dirname(os.path.abspath(path)),
                            hfs_A_theory.SUPERSEDED_NAME)
    output_files.require_writable([path, old_path])
    with open(path, encoding='utf-8', newline='') as fh:
        rd = csv.DictReader(fh)
        fields = rd.fieldnames
        table = {r['level_id']: r for r in rd}
    superseded, added, replaced = [], 0, 0
    for x in out:
        act = x['action']
        if not (act.startswith('add') or act.startswith('replace')):
            continue
        lid = x['level_id']
        old = table.get(lid)
        if old is not None and kind_of(old['source']) == 'measured':
            reason = USER_DECISIONS.get(
                lid, 'replaced: the calculated A is much more accurate and '
                     'agrees')
            superseded.append(dict(old, superseded_on=date, reason=reason))
        table[lid] = {
            'level_id': lid,
            'cfg': old['cfg'] if old and old['cfg'] else x['proposed_cfg'],
            'J': old['J'] if old else x['J'],
            'A_cm-1': x['proposed_A'],
            'u_A': x['proposed_u_A'],
            'source': x['proposed_source'],
            'n_flagged': x['n_flagged'],
        }
        if old is None:
            added += 1
        else:
            replaced += 1
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=fields, lineterminator='\n')
        wr.writeheader()
        wr.writerows(table[k] for k in sorted(table))
    if superseded:
        new_file = not os.path.isfile(old_path)
        with open(old_path, 'a', encoding='utf-8', newline='') as fh:
            wr = csv.DictWriter(fh, fieldnames=SUPERSEDED_FIELDS,
                                lineterminator='\n')
            if new_file:
                wr.writeheader()
            wr.writerows(superseded)
    return added, replaced, len(superseded)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--lines', default=LINES,
                    help='line table of the working set')
    ap.add_argument('--no-write', action='store_true')
    ap.add_argument('--write', action='store_true',
                    help='also write the adopted rows into A_hfs_levels.csv')
    args = ap.parse_args(argv)

    import cowan_gA
    import eigenvector_check
    import hfs_A_fit
    import hfs_A_theory
    import hfs_correction
    import hfs_kappa

    theory = {r['level_id']: r for r in read_csv(THEORY_CSV) if r['level_id']}
    eig = eigenvector_check.read_status()
    if not eig:
        print('eigenvector_check.csv not found: no intensity veto')
    table = {r['level_id']: r for r in read_csv(A_LEVELS)}
    iden = {lid: row for row, lid in cowan_gA.read_id_map(ID_MAP).items()}
    questionable = {r['level_id']: r['status']
                    for r in read_csv(QUESTIONABLE, '\t')
                    if r['status'] != 'firm'}
    with open(CONFIG, 'rb') as fh:
        resolved = set(tomllib.load(fh).get('hfs', {})
                       .get('resolved_levels', []))

    def J_of(lid):
        r = theory.get(lid) or table.get(lid)
        return float(r['J']) if r else None

    def known(lid):
        """(A, u) the pipeline has for a level, on the measured scale."""
        r = table.get(lid)
        if r is None or not hfs_correction.is_determined(r['source']):
            return None
        return hfs_kappa.scaled_A(float(r['A_cm-1']), float(r['u_A']),
                                  r['source'])

    def calc(lid):
        r = theory.get(lid)
        if not r or not r['A_calc']:
            return None
        return float(r['A_calc']), float(r['u_total'])

    lines = [r for r in read_csv(args.lines) if r['accepted'] == '1.0']
    E_obs = {}
    for r in lines:
        E_obs[r['low_id']] = float(r['low_E'])
        E_obs[r['upp_id']] = float(r['upp_E'])

    def E_dev(lid, r):
        e = _f(r['E_obs']) if r['E_obs'] else E_obs.get(lid)
        return None if e is None or not r['E_calc'] else e - float(r['E_calc'])

    discarded = read_discarded()
    candidates, n_semi_better, keeps_table = {}, 0, set()
    for lid, r in theory.items():
        if not r['A_calc'] or lid in discarded:
            continue
        t = table.get(lid)
        kind = kind_of(t['source'] if t else None)
        if kind in ('measured', 'semiempirical'):
            A_t, u_t = float(t['A_cm-1']), float(t['u_A'])
            A_c, u_c = float(r['A_calc']), float(r['u_total'])
            better = u_c <= hfs_A_fit.MUCH_MORE_ACCURATE * u_t
            agrees = abs(A_c - A_t) <= hfs_A_fit.AGREE_SIGMA * math.hypot(
                u_c, u_t)
            if kind == 'semiempirical':
                n_semi_better += better
                continue
            if not (better or not agrees or lid in USER_DECISIONS):
                continue
            if not ((better and agrees) or lid in USER_DECISIONS):
                keeps_table.add(lid)
        candidates[lid] = (r, kind)

    def partner(lid):
        """A partner's (A, u): the calculated one where the table has none
        or is to take it, else the table's."""
        if lid in candidates and lid not in keeps_table:
            return calc(lid)
        return known(lid) or calc(lid)

    a_table = hfs_A_fit.read_A_table(A_LEVELS)
    patterns, _ = hfs_A_fit.read_patterns(args.lines, a_table,
                                          hfs_kappa.read_levels())

    flagged = [r for r in lines if r['char'] in ('*r', '*v')]
    wn_lo = min(float(r['wn_obs']) for r in flagged)
    wn_hi = max(float(r['wn_obs']) for r in flagged)
    i_min = float(np.percentile([float(r['obs_intens'] or 0)
                                 for r in flagged], 10))

    by_level = {}
    for r in lines:
        for lid in (r['low_id'], r['upp_id']):
            if lid in candidates:
                by_level.setdefault(lid, []).append(r)

    out = []
    for lid, (r, kind) in candidates.items():
        A = _f(r['A_calc'])
        J = float(r['J'])
        u_total = float(r['u_total'])
        cancel = _f(r['cancel'])
        dev = E_dev(lid, r)
        ev = eig.get(lid, {})
        eig_fail = ev.get('status') == 'fail'
        tier = tier_of(A, u_total, cancel, r['tested'] == 'yes', dev,
                       eig_fail)
        S = I_NUC * A * J
        t = table.get(lid)
        A_t = float(t['A_cm-1']) if t else None
        u_t = float(t['u_A']) if t else None
        z = None
        if kind == 'measured':
            z = (A - A_t) / math.hypot(u_total, u_t)
        better = (kind == 'measured'
                  and u_total <= hfs_A_fit.MUCH_MORE_ACCURATE * u_t)
        agrees = z is None or abs(z) <= hfs_A_fit.AGREE_SIGMA

        mine = by_level.get(lid, [])
        w = [(_f(x['BF']) or 1.0) / (_f(x['unc_wn_obs']) or 0.01) ** 2
             for x in mine]
        k = [float(x['kappa']) for x in mine]
        kappa = (sum(a * b for a, b in zip(w, k)) / sum(w)) if mine else None
        if lid in resolved:
            kappa = None            # its S stays 0: nothing moves

        n_flag, disagree, wide = 0, [], []
        for x in mine:
            lo, up = x['low_id'], x['upp_id']
            if lo in resolved or up in resolved:
                continue
            pa_lo = (A, 0) if lo == lid else partner(lo)
            pa_up = (A, 0) if up == lid else partner(up)
            J_lo, J_up = J_of(lo), J_of(up)
            if None in (pa_lo, pa_up, J_lo, J_up):
                continue
            D = I_NUC * (pa_up[0] * J_up - pa_lo[0] * J_lo)
            wn = float(x['wn_obs'])
            if x['char'] in ('*r', '*v'):
                n_flag += 1
                if D * (1 if x['char'] == '*r' else -1) <= 0:
                    disagree.append('%s%s(D%+.2f)' % (x['wn_obs'], x['char'],
                                                       D))
            elif (x['char'] == '' and abs(D) >= WIDE_D
                  and wn_lo <= wn <= wn_hi
                  and float(x['obs_intens'] or 0) >= i_min):
                wide.append('%s(D%+.2f)' % (x['wn_obs'], D))

        comp = components_value(patterns, lid, partner)

        notes = []
        if lid in questionable:
            notes.append('questionable level: the calculated A will show '
                         'whether its lines line up in IDEN2')
        if lid in resolved:
            notes.append('resolved level: its S stays 0 in the pipeline '
                         'whatever its A')
        if lid in hfs_A_fit.HOLD:
            notes.append('HOLD in hfs_A_fit: ' + hfs_A_fit.HOLD[lid])
        if lid in USER_DECISIONS:
            notes.append(USER_DECISIONS[lid])
        if dev is not None and abs(dev) > E_DEV_MAX:
            notes.append('|E_obs - E_calc| = %.0f cm^-1: eigenvector suspect'
                         % abs(dev))
        if eig_fail:
            notes.append('line intensities contradict the eigenvector '
                         '(eigenvector_check: chi^2/n %s, %s dof)'
                         % (ev['chi2_n'], ev['dof']))
        if (_f(r['w_missing']) or 0.0) >= W_NOTE:
            notes.append('weight %.3f without parameters (%s)'
                         % (float(r['w_missing']), r['missing']))

        decided = lid in USER_DECISIONS
        test = (tier == 3 and not decided and lid not in resolved
                and kind in ('absent', 'undetermined', 'calculated'))
        if test:
            notes.append('tier 3 adopted as a test: its lines show whether '
                         'the calculated A holds')
            prop_A, prop_u = A, test_u(A, u_total, cancel)
            prop_src = '%s, hfs_A_theory %s, tier 3, %s' % (
                hfs_A_theory.CALC_SOURCE, DATE, TEST_TAG)
        elif tier == 3 and not decided:
            prop_A, prop_u = 0.0, undetermined_u(A, u_total, cancel)
            prop_src = 'not determined'
        else:
            prop_A, prop_u = A, u_total
            prop_src = '%s, hfs_A_theory %s, tier %d' % (
                hfs_A_theory.CALC_SOURCE, DATE, tier)
        out.append({
            'level_id': lid,
            'iden2_row': iden.get(lid, ''),
            'J': r['J'],
            'leading': r['leading'],
            'table_now': 'absent' if t is None else t['source'],
            'A_table': '' if t is None else t['A_cm-1'],
            'u_table': '' if t is None else t['u_A'],
            'tier': tier,
            'A_calc': '%+.4f' % A,
            'u_total': r['u_total'],
            'u_amp': r['u_amp'],
            'u_par': r['u_par'],
            'u_cfg': r['u_cfg'],
            'rel_u': '%.2f' % (u_total / abs(A)) if A else '',
            'cancel': r['cancel'],
            'z_table': '' if z is None else '%+.1f' % z,
            'eigenvector': ('%s (%s, %s)' % (ev['status'], ev['chi2_n'],
                                             ev['dof'])
                            if ev.get('chi2_n') else ev.get('status', '')),
            'E_dev': '' if dev is None else '%+.0f' % dev,
            'S': '%+.3f' % S,
            'u_S': '%.3f' % (I_NUC * J * u_total),
            'n_lines': len(mine),
            'kappa_mean': '' if kappa is None else '%.3f' % kappa,
            'dE_expected': '' if kappa is None else '%+.3f' % (-kappa * S),
            'n_flagged': n_flag,
            'flag_disagree': ' '.join(disagree),
            'wide_unflagged': ' '.join(wide),
            'A_components': ('' if comp is None else '%+.4f +- %.4f (%d)'
                             % comp),
            'action': action_of(tier, kind, better, agrees, decided,
                                lid in resolved, test),
            'proposed_cfg': r['leading'].split()[1],
            'proposed_A': '%+.4f' % prop_A,
            'proposed_u_A': '%.4f' % prop_u,
            'proposed_source': prop_src,
            'notes': '; '.join(notes),
        })
    out.sort(key=lambda x: (x['tier'], -abs(float(x['dE_expected'] or 0))))

    print('levels with a calculated A that has a use: %d (discarded levels '
          'left out: %d)' % (len(out), len(discarded & set(theory))))
    for tier in (1, 2, 3):
        sel = [x for x in out if x['tier'] == tier]
        dE = [abs(float(x['dE_expected'])) for x in sel if x['dE_expected']]
        confs = {}
        for x in sel:
            confs[x['proposed_cfg']] = confs.get(x['proposed_cfg'], 0) + 1
        print('  tier %d: %3d levels (%s); |dE_expected| median %.3f, max '
              '%.3f cm^-1' % (tier, len(sel), ', '.join(
                  '%s %d' % kv for kv in sorted(confs.items(),
                                                key=lambda kv: -kv[1])),
                  np.median(dE) if dE else 0, max(dE) if dE else 0))
    acts = {}
    for x in out:
        acts[x['action']] = acts.get(x['action'], 0) + 1
    for a, n in sorted(acts.items(), key=lambda kv: -kv[1]):
        print('  %-48s %3d' % (a, n))
    print('  semiempirical rows the calculation is much more accurate than, '
          'kept: %d' % n_semi_better)
    print('  flagged lines checked: %d, contradicted: %d; unflagged lines '
          'with |D| >= %.1f: %d (range %.0f-%.0f cm^-1, I >= %.0f)'
          % (sum(x['n_flagged'] for x in out),
             sum(len(x['flag_disagree'].split()) for x in out), WIDE_D,
             sum(len(x['wide_unflagged'].split()) for x in out),
             wn_lo, wn_hi, i_min))
    print('  with a component value: %d'
          % sum(bool(x['A_components']) for x in out))

    if args.write:
        added, replaced, old = write_table(out)
        print('A_hfs_levels.csv: %d rows added, %d replaced (%d measured '
              'rows to %s)' % (added, replaced, old,
                               hfs_A_theory.SUPERSEDED_NAME))
    if not args.no_write:
        import output_files
        output_files.require_writable([OUT_CSV])
        with open(OUT_CSV, 'w', encoding='utf-8', newline='') as fh:
            wr = csv.DictWriter(fh, fieldnames=list(out[0].keys()),
                                lineterminator='\n')
            wr.writeheader()
            wr.writerows(out)
        print('wrote %s' % os.path.basename(OUT_CSV))
    return out


if __name__ == '__main__':
    main()
