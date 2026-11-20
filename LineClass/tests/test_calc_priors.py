"""Tests for the calculated A constants as priors of the anchored fit, used
only where no composition prior reaches a pattern (2026-11-20):
hfs_A_fit.full_fit, calc_combinable, write_combined, calc_priors_of;
hfs_kappa.is_calc_combined; hfs_A_theory.read_measured; and
hfs_A_candidates' "combined" kind and refresh_priors.

The first case was 24593.418 *r (000310 - 000655, J 11/2 - 11/2): its two
printed companions measure only A(000655) - A(000310).

Run from the LineClass directory:  python -m pytest tests -q
"""
import pytest

import hfs_A_candidates as C
import hfs_A_fit
import hfs_kappa
import hfs_patterns

CALC = 'calculated, hfs_A_theory 2026-10-08, tier 2'


def pattern(low, upp, J1, J2, A1, A2, n, wn=30000.0, char='*v'):
    """A pattern whose n listed companions sit exactly on the ladder."""
    dwn = [hfs_patterns.rung(k, J1, A1, J2, A2) for k in range(1, n + 1)]
    return hfs_patterns.Pattern(wn, char, low, upp, J1, J2, A1, A2, 0.0, 0.0,
                                '', '', dwn)


def world():
    """An anchored part (L1, with a composition prior, and F1) and an
    isolated pair X1 - X2 of equal J."""
    pats = [pattern('L1', 'F1', 4.5, 3.5, 0.090, 0.025, 3),
            pattern('X1', 'X2', 5.5, 5.5, -0.0107, 0.0318, 2, wn=24593.0,
                    char='*r')]
    anchors = {'L1': (0.090, 0.003)}
    calc = {'X1': (-0.0051, 0.0109), 'X2': (0.0309, 0.0043),
            'F1': (0.02, 0.01)}
    return pats, anchors, calc


def test_calculated_priors_are_used_only_where_nothing_else_reaches():
    pats, anchors, calc = world()
    res0 = hfs_A_fit.fit(pats, anchors)
    assert res0.null == {'X1', 'X2'}
    res, extra, dependent = hfs_A_fit.full_fit(pats, anchors, calc)
    # F1 has a calculated A too, but the composition prior reaches it
    assert set(extra) == {'X1', 'X2'} and dependent == {'X1', 'X2'}
    assert not res.null
    assert res.A['F1'] == pytest.approx(res0.A['F1'], abs=1e-12)
    assert res.u['F1'] == pytest.approx(res0.u['F1'], abs=1e-12)
    # the spacings fix the difference, and the priors share out the rest by
    # their variances: nearly all of it goes to the less certain X1
    assert res.A['X2'] - res.A['X1'] == pytest.approx(0.0425, abs=3e-4)
    assert res.u['X1'] / calc['X1'][1] < 0.5
    assert res.u['X2'] / calc['X2'][1] > 0.85
    # with nothing undetermined nothing is added
    assert hfs_A_fit.full_fit(pats[:1], anchors, calc)[1:] == ({}, set())


def table():
    return {'X1': {'A_cm-1': '-0.0051', 'u_A': '0.0109', 'source': CALC},
            'X2': {'A_cm-1': '+0.0309', 'u_A': '0.0043', 'source': CALC},
            'X3': {'A_cm-1': '+0.0000', 'u_A': '0.0300',
                   'source': 'not determined'}}


def test_a_combination_is_written_only_where_it_improves_the_value():
    pats, anchors, calc = world()
    res, extra, dependent = hfs_A_fit.full_fit(pats, anchors, calc)
    t = table()
    assert hfs_A_fit.calc_combinable('X1', t, res, extra, dependent)
    assert not hfs_A_fit.calc_combinable('X2', t, res, extra, dependent)
    # a prior the spacings pull too far is not combined
    pulled = res._replace(pulls=dict(res.pulls, X1=3.5))
    assert not hfs_A_fit.calc_combinable('X1', t, pulled, extra, dependent)
    # an undetermined row determined through its partners' priors is ...
    res3 = res._replace(A=dict(res.A, X3=0.01), u=dict(res.u, X3=0.004))
    assert hfs_A_fit.calc_combinable('X3', t, res3, extra,
                                     dependent | {'X3'})
    # ... but never by --write-levels, as if it were a measurement
    assert not hfs_A_fit.adoptable('X3', t, res3, dependent | {'X3'})
    assert hfs_A_fit.adoptable('X3', t, res3)


def test_writing_a_calculated_combination_keeps_its_prior(tmp_path):
    pats, anchors, calc = world()
    res, extra, dependent = hfs_A_fit.full_fit(pats, anchors, calc)
    path, priors_path = tmp_path / 'A.csv', tmp_path / 'priors.csv'
    path.write_text('level_id,cfg,J,A_cm-1,u_A,source,n_flagged\n'
                    'X1,f5d6s,5.5,-0.0051,0.0109,"%s",0\n'
                    'X2,fd6p,5.5,+0.0309,0.0043,"%s",1\n' % (CALC, CALC),
                    encoding='utf-8', newline='\n')
    rows = hfs_A_fit.read_A_table(str(path))
    changed = hfs_A_fit.write_combined(
        rows, {}, res, anchors, {'X1': (1, 2), 'X2': (1, 2)}, str(path),
        str(priors_path), '2026-11-20', extra, dependent)
    assert changed == ['X1']
    rows = hfs_A_fit.read_A_table(str(path))
    assert rows['X1']['A_cm-1'] == '%+.4f' % res.A['X1']
    assert rows['X1']['u_A'] == '%.4f' % res.u['X1']     # no scale added
    src = rows['X1']['source']
    assert src.startswith(hfs_A_fit.CALC_COMBINED_PREFIX)
    assert hfs_kappa.is_calc_combined(src)
    assert not hfs_kappa.is_combined(src)
    assert not hfs_kappa.is_semiempirical(src)
    assert rows['X2']['source'] == CALC                   # gains too little
    priors = hfs_A_fit.read_priors(str(priors_path))
    assert list(priors) == ['X1'] and priors['X1']['source'] == CALC
    assert priors['X1']['moved'] == '2026-11-20'
    # the next run takes the prior, never the combined value, and no
    # calculated constant becomes an anchor of the whole fit
    again = hfs_A_fit.calc_priors_of(rows, priors)
    assert again == {'X1': (-0.0051, 0.0109), 'X2': (0.0309, 0.0043)}
    assert hfs_A_fit.anchors_of(rows, priors=priors) == {}
    res2, _, _ = hfs_A_fit.full_fit(pats, anchors, again)
    assert res2.A['X1'] == pytest.approx(res.A['X1'], abs=1e-12)
    assert b'\r' not in path.read_bytes()


def test_a_calculated_combination_without_its_prior_stops_the_fit():
    src = hfs_A_fit.calc_combined_source(1, 2, '2026-11-20')
    t = {'X1': {'A_cm-1': '-0.0107', 'u_A': '0.0042', 'source': src}}
    with pytest.raises(SystemExit, match='no prior'):
        hfs_A_fit.calc_priors_of(t)
    with pytest.raises(SystemExit, match='no prior'):
        hfs_A_fit.anchors_of(t)


def test_a_calculated_combination_is_no_measurement_for_the_theory(tmp_path):
    import hfs_A_theory
    import hfs_correction
    src = hfs_A_fit.calc_combined_source(1, 2, '2026-11-20')
    assert hfs_correction.is_determined(src)
    assert hfs_kappa.scaled_A(-0.0107, 0.0042, src) == (-0.0107, 0.0042)
    path = tmp_path / 'A.csv'
    path.write_text('level_id,cfg,J,A_cm-1,u_A,source,n_flagged\n'
                    'X1,f5d6s,5.5,-0.0107,0.0042,"%s",0\n'
                    'K,f26s,4.5,+0.0900,0.0020,"resolved components, 2 '
                    'patterns, 3 positions (Step 4 anchored fit, '
                    '2026-10-05)",0\n' % src, encoding='utf-8', newline='\n')
    assert set(hfs_A_theory.read_measured(str(path))) == {'K'}


# ---- hfs_A_candidates ----------------------------------------------------------

def test_a_calculated_combination_is_its_own_kind_and_is_kept():
    src = hfs_A_fit.calc_combined_source(1, 2, '2026-11-20')
    assert C.kind_of(src) == 'combined'
    for tier in (1, 2, 3):
        assert C.action_of(tier, 'combined', True, True, False,
                           False).startswith('keep')
        assert C.action_of(tier, 'combined', True, True, False, False,
                           stale_prior=True) == C.PRIOR_REFRESH


def test_refresh_priors_renews_only_the_prior(tmp_path):
    path = tmp_path / 'priors.csv'
    path.write_text('level_id,cfg,J,A_cm-1,u_A,source,n_flagged,moved\n'
                    'X1,f5d6s,5.5,-0.0051,0.0109,"calculated, old",0,'
                    '2026-11-20\n'
                    'C1,f26p,2.5,+0.0251,0.0055,composition,5,2026-10-09\n',
                    encoding='utf-8', newline='\n')
    out = [{'level_id': 'X1', 'J': '5.5', 'proposed_cfg': 'f5d6s',
            'action': C.PRIOR_REFRESH, 'proposed_A': '-0.0060',
            'proposed_u_A': '0.0100', 'proposed_source': 'calculated, new'},
           {'level_id': 'C1', 'action': 'keep: semiempirical prior'}]
    assert C.refresh_priors(out, str(path)) == ['X1']
    p = hfs_A_fit.read_priors(str(path))
    assert (p['X1']['A_cm-1'], p['X1']['u_A'], p['X1']['source'],
            p['X1']['moved']) == ('-0.0060', '0.0100', 'calculated, new',
                                  '2026-11-20')
    assert p['C1']['A_cm-1'] == '+0.0251'
