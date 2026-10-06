"""Tests for `hfs_A_fit.py`, the anchored fit of Step 4.

The arithmetic is tested on synthetic patterns made from known constants:
the fit must give them back, must report a level no prior reaches as not
determined, and must replace only the rows it is allowed to replace.
"""

import csv
import math

import pytest

import hfs_A_fit
import hfs_patterns


def pattern(low, upp, J1, J2, A1, A2, n, wn=30000.0, char='*v'):
    """A pattern whose n listed companions sit exactly on the ladder."""
    dwn = [hfs_patterns.rung(k, J1, A1, J2, A2) for k in range(1, n + 1)]
    return hfs_patterns.Pattern(wn, char, low, upp, J1, J2, A1, A2, 0.0, 0.0,
                                '', '', dwn)


TRUE = {'L1': 0.090, 'L2': -0.040, 'F1': 0.025, 'F2': 0.031, 'F3': 0.012}
J = {'L1': 4.5, 'L2': 2.5, 'F1': 3.5, 'F2': 5.5, 'F3': 2.5}


def patterns_of(pairs):
    return [pattern(a, b, J[a], J[b], TRUE[a], TRUE[b], n, wn=29000.0 + i)
            for i, (a, b, n) in enumerate(pairs)]


def test_the_fit_gives_back_the_constants_it_was_made_from():
    pats = patterns_of([('L1', 'F1', 3), ('L2', 'F1', 2), ('L1', 'F2', 4),
                        ('L2', 'F3', 2), ('F3', 'F2', 1)])
    anchors = {'L1': (TRUE['L1'], 0.003), 'L2': (TRUE['L2'], 0.003)}
    res = hfs_A_fit.fit(pats, anchors)
    assert not res.null
    for x in ('F1', 'F2', 'F3'):
        assert res.A[x] == pytest.approx(TRUE[x], abs=1e-9)
        assert 0 < res.u[x] < 0.01
    assert res.chi2 == pytest.approx(0.0, abs=1e-12)


def test_a_cluster_no_prior_reaches_is_not_determined():
    """Adding one constant to every A of a cluster changes its spacings only
    through J1 - J2; without a prior the cluster's offset is free."""
    pats = patterns_of([('L1', 'F1', 3)])
    pats += [pattern('X1', 'X2', 1.5, 1.5, 0.01, 0.03, 1, wn=11000.0)]
    res = hfs_A_fit.fit(pats, {'L1': (TRUE['L1'], 0.003)})
    assert res.null == {'X1', 'X2'}
    assert 'F1' not in res.null


def test_a_pulled_prior_is_reported():
    pats = patterns_of([('L1', 'F1', 3), ('L2', 'F1', 3)])
    anchors = {'L1': (TRUE['L1'] + 0.02, 0.003), 'L2': (TRUE['L2'], 0.003)}
    res = hfs_A_fit.fit(pats, anchors)
    assert abs(res.pulls['L1']) > 2


def test_each_pattern_alone_gives_the_same_value_on_exact_data():
    pats = patterns_of([('L1', 'F1', 3), ('L2', 'F1', 2)])
    res = hfs_A_fit.fit(pats, {'L1': (TRUE['L1'], 0.003),
                               'L2': (TRUE['L2'], 0.003)})
    single = hfs_A_fit.single_pattern_values(pats, res, 'F1')
    assert len(single) == 2
    for _, a, u in single:
        assert a == pytest.approx(TRUE['F1'], abs=1e-9)
        assert u > 0


def test_the_weight_is_sugars_precision_taken_twice():
    # 0.003 A at 3000 A (33333 cm^-1) is 0.0333 cm^-1; a spacing, sqrt(2) more
    wn = 1e8 / 3000.0
    assert hfs_A_fit.sigma_of(wn) == pytest.approx(
        math.sqrt(2) * 0.003 * 1e-8 * wn * wn)


def flag_row(wn, char, low, upp):
    return {'wn_obs': wn, 'char': char, 'low_id': low, 'upp_id': upp}


def test_flags_bound_A_from_the_right_side():
    """D = I (A2 J2 - A1 J1): an *r line needs D > 0, a *v line D < 0."""
    J_of = {'P': 2.5, 'X': 2.5}
    A_of = {'P': 0.04}.get
    known = lambda x: x == 'P'
    # X upper, *r: A_X J_X > A_P J_P, so A_X > 0.04
    lo, hi, n = hfs_A_fit.flag_bounds([flag_row('1', '*r', 'P', 'X')], 'X',
                                      A_of, J_of, known)
    assert (lo, hi, n) == (pytest.approx(0.04), math.inf, 1)
    # X upper, *v: A_X < 0.04
    lo, hi, n = hfs_A_fit.flag_bounds([flag_row('1', '*v', 'P', 'X')], 'X',
                                      A_of, J_of, known)
    assert (lo, hi) == (-math.inf, pytest.approx(0.04))
    # X lower, *r: A_P J_P > A_X J_X, so A_X < 0.04
    lo, hi, n = hfs_A_fit.flag_bounds([flag_row('1', '*r', 'X', 'P')], 'X',
                                      A_of, J_of, known)
    assert (lo, hi) == (-math.inf, pytest.approx(0.04))
    # a blend bounds nothing
    rows = [flag_row('1', '*r', 'P', 'X'), flag_row('1', '*r', 'Q', 'X')]
    assert hfs_A_fit.flag_bounds(rows, 'X', A_of, J_of, known)[2] == 0


def test_the_lines_route_recovers_A_and_the_level_shift():
    """Lines of a level whose A was taken as 0: a class of kappa sits
    (1 - kappa) I J A off; two classes separate A from the level shift."""
    A, delta, J_x = 0.03, 0.05, 3.5
    rows = []
    for kappa, s in ((0.926, 1), (0.717, 1), (0.276, -1), (0.717, -1),
                     (0.926, -1)):
        f = (1 - kappa) * 2.5 * J_x
        oc = -(s * f * A - s * delta)
        low, upp = ('o', 'X') if s == 1 else ('X', 'o')
        rows.append({'low_id': low, 'upp_id': upp, 'kappa': str(kappa),
                     'BF': '1.0', 'dif_wn_O-C': repr(oc),
                     'unc_wn_obs': '0.05'})
    a, u, n = hfs_A_fit.lines_value(rows, 'X', 0.0, J_x)
    assert n == 5
    assert a == pytest.approx(A, abs=1e-9)


def test_write_levels_replaces_only_what_it_may(tmp_path):
    path = tmp_path / 'A.csv'
    path.write_text(
        'level_id,cfg,J,A_cm-1,u_A,source,n_flagged\n'
        'L1,f26s,4.5,+0.0934,0.0030,Reader & Sugar 1965 (semiempirical),0\n'
        'C1,f26p,2.5,+0.0251,0.0055,composition,5\n'
        'F1,f5d2,3.5,+0.0000,0.0500,not determined,2\n'
        'F2,f5d2,5.5,+0.0300,0.0480,flag interval,1\n'
        'F3,f5d2,2.5,+0.0100,0.0080,"center of gravity of 1, 2026-10-02",0\n'
        'H,f5d6s,4.5,+0.0000,0.0500,not determined,1\n',
        encoding='utf-8', newline='\n')
    table = hfs_A_fit.read_A_table(str(path))
    res = hfs_A_fit.Fit(
        ids=['L1', 'C1', 'F1', 'F2', 'F3', 'H'],
        A={'L1': 0.09, 'C1': 0.024, 'F1': 0.025, 'F2': 0.031, 'F3': 0.02,
           'H': 0.08},
        u={'L1': 0.003, 'C1': 0.004, 'F1': 0.002, 'F2': 0.002, 'F3': 0.002,
           'H': 0.007},
        chi2=0.0, dof=1, n_positions=1, n_priors=1, null=set(), pulls={},
        resid_rms=0.0)
    counts = {x: (2, 3) for x in res.ids}
    saved = hfs_A_fit.HOLD
    hfs_A_fit.HOLD = {'H': 'test'}
    try:
        changed = hfs_A_fit.write_levels(table, res, counts, str(path),
                                         '2026-10-04')
    finally:
        hfs_A_fit.HOLD = saved
    # Reader and Sugar's value is freed and replaced (2026-10-05)
    assert changed == ['L1', 'F1', 'F2']
    raw = path.read_bytes()
    assert b'\r' not in raw
    rows = {r['level_id']: r for r in csv.DictReader(
        raw.decode('utf-8').splitlines())}
    assert rows['F1']['A_cm-1'] == '+0.0250'
    assert rows['F1']['source'].startswith('resolved components, 2 patterns')
    assert rows['F2']['A_cm-1'] == '+0.0310'
    assert rows['L1']['A_cm-1'] == '+0.0900'
    assert rows['L1']['source'].startswith('resolved components')
    # the composition anchor, the hand measurement and the held level stay
    assert rows['C1']['A_cm-1'] == '+0.0251'
    assert rows['F3']['source'].startswith('center of gravity')
    assert rows['H']['source'] == 'not determined'


def test_reader_and_sugar_are_no_priors():
    table = {
        'L1': {'A_cm-1': '+0.0934', 'u_A': '0.0030',
               'source': 'Reader & Sugar 1965 (semiempirical)'},
        'C1': {'A_cm-1': '+0.0251', 'u_A': '0.0055', 'source': 'composition'},
        'F2': {'A_cm-1': '+0.0300', 'u_A': '0.0480',
               'source': 'flag interval'}}
    assert set(hfs_A_fit.anchors_of(table)) == {'C1'}


def test_a_written_source_is_read_as_a_measurement():
    """The new source must not be mistaken for a semiempirical constant (which
    would be scaled again) nor for an undetermined one."""
    import hfs_correction
    import hfs_kappa
    src = hfs_A_fit.source_of(3, 7, '2026-10-04')
    assert not hfs_kappa.is_semiempirical(src)
    assert hfs_correction.is_determined(src)
    assert hfs_kappa.scaled_A(0.02, 0.002, src) == (0.02, 0.002)


def test_a_level_the_fit_never_met_is_left_alone():
    table = {'Z': {'level_id': 'Z', 'source': 'not determined'}}
    res = hfs_A_fit.Fit(ids=[], A={}, u={}, chi2=0.0, dof=1, n_positions=0,
                        n_priors=0, null=set(), pulls={}, resid_rms=0.0)
    assert not hfs_A_fit.adoptable('Z', table, res)


def test_a_second_run_refreshes_its_own_values_only():
    table = {'F1': {'level_id': 'F1', 'source': hfs_A_fit.source_of(
                 2, 3, '2026-10-04')},
             'C': {'level_id': 'C', 'source': 'composition'}}
    res = hfs_A_fit.Fit(ids=['F1', 'C'], A={'F1': 0.02, 'C': 0.03},
                        u={'F1': 0.002, 'C': 0.003}, chi2=0.0, dof=1,
                        n_positions=0, n_priors=0, null=set(), pulls={},
                        resid_rms=0.0)
    assert hfs_A_fit.adoptable('F1', table, res)
    assert not hfs_A_fit.adoptable('C', table, res)
