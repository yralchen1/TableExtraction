"""Tests for hfs_A_candidates.py: the tiers of the calculated A constants
proposed for levels that have none (2026-10-07).

Run from the LineClass directory:  python -m pytest tests -q
"""
import csv
import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import hfs_A_candidates as C            # noqa: E402
import hfs_kappa                        # noqa: E402
import hfs_correction                   # noqa: E402
import hfs_A_theory                     # noqa: E402


def test_tiers():
    assert C.tier_of(0.03, 0.003, 1.1, True) == 1
    assert C.tier_of(0.03, 0.003, 1.1, False) == 2
    assert C.tier_of(0.03, 0.003, 2.0, True) == 3        # cancellation
    assert C.tier_of(0.03, 0.008, 1.1, True) == 3        # over 25 per cent
    assert C.tier_of(0.03, 0.0075, 1.1, True) == 1       # exactly 25
    assert C.tier_of(None, 0.0, None, True) == 3
    assert C.tier_of(0.03, 0.003, 1.1, True, E_dev=-299.0) == 1
    assert C.tier_of(0.03, 0.003, 1.1, True, E_dev=+301.0) == 3  # eigenvector
    assert C.tier_of(0.03, 0.003, 1.1, True, eig_fail=True) == 3  # intensities


def test_undetermined_u_covers_value_and_cancelling_terms():
    # A = 0.004 left from terms summing to 0.04 in magnitude (cancel 10)
    u = C.undetermined_u(0.004, 0.003, 10.0)
    assert u == pytest.approx(math.sqrt(0.004 ** 2 + 0.003 ** 2
                                        + 0.004 ** 2))
    assert C.undetermined_u(None, 0.002, None) == pytest.approx(0.002)


def test_a_calculated_row_is_used_unscaled_and_as_determined():
    """The proposed source must reach the pipeline as a determined A on the
    measured scale, not as a semiempirical one A_SCALE multiplies again."""
    src = '%s, hfs_A_theory %s, tier 1' % (hfs_A_theory.CALC_SOURCE, C.DATE)
    assert not hfs_kappa.is_semiempirical(src)
    assert hfs_kappa.scaled_A(0.0213, 0.0013, src) == (0.0213, 0.0013)
    assert hfs_correction.is_determined(src)


def test_actions():
    assert C.action_of(1, 'calculated', False, True, False, False) \
        == 'in the table'
    # a calculated row whose level has fallen to tier 3 becomes a test row
    # (user, 2026-10-08; it was withdrawn to "not determined" before)
    assert C.action_of(3, 'calculated', False, True, False, False) \
        == 'replace as a test: tier 3'
    assert C.action_of(1, 'absent', False, True, False, False) == 'add'
    # a tier-3 level without a row is adopted as a test (user, 2026-10-08),
    # and its test row then stays
    assert C.action_of(3, 'absent', False, True, False, False) \
        == 'add as a test: tier 3'
    assert C.action_of(3, 'calculated', False, True, False, False, True) \
        == 'in the table'
    assert C.action_of(2, 'undetermined', False, True, False, False) \
        == 'replace'
    assert C.action_of(3, 'undetermined', False, True, False, False) \
        == 'replace as a test: tier 3'
    # a measured A: replaced only by a much more accurate, agreeing value
    assert C.action_of(1, 'measured', True, True, False, False) == 'replace'
    assert C.action_of(1, 'measured', False, True, False, False).startswith(
        'keep')
    assert C.action_of(1, 'measured', True, False, False, False).startswith(
        'user decides')
    assert C.action_of(2, 'measured', False, False, True, False).startswith(
        'replace')
    assert C.action_of(3, 'measured', True, True, False, False).startswith(
        'keep')
    assert C.action_of(1, 'measured', True, True, False, True).startswith(
        'keep')
    # a resolved level never takes a calculated A, even with no row
    assert C.action_of(1, 'absent', False, True, False, True).startswith(
        'keep')
    assert C.action_of(1, 'semiempirical', True, True, False, False) \
        .startswith('keep')


def test_test_rows():
    """A tier-3 test row: determined, unscaled, recognized by its source;
    its u_A adds 10 per cent of the terms of A to u_total."""
    src = '%s, hfs_A_theory %s, tier 3, %s' % (hfs_A_theory.CALC_SOURCE,
                                               C.DATE, C.TEST_TAG)
    assert C.is_test_row(src) and C.kind_of(src) == 'calculated'
    assert hfs_correction.is_determined(src)
    assert hfs_kappa.scaled_A(0.02, 0.004, src) == (0.02, 0.004)
    assert not C.is_test_row('%s, hfs_A_theory %s, tier 2'
                             % (hfs_A_theory.CALC_SOURCE, C.DATE))
    assert C.test_u(0.02, 0.003, 2.5) == pytest.approx(
        math.hypot(0.003, 0.1 * 2.5 * 0.02))


def test_read_discarded(tmp_path):
    path = tmp_path / 'discarded_levels.csv'
    path.write_text('level_id,iden2_row,date,ln_R,reason\n'
                    '059003.000513,436,2026-09-20,-3.887,"no support"\n',
                    encoding='utf-8')
    assert C.read_discarded(str(path)) == {'059003.000513'}
    assert C.read_discarded(str(tmp_path / 'none.csv')) == set()


def test_components_value_gives_back_the_constant():
    import hfs_patterns
    A1, A2, J1, J2 = 0.031, -0.012, 3.5, 2.5
    dwn = [hfs_patterns.rung(k, J1, A1, J2, A2) for k in (1, 2)]
    p = hfs_patterns.Pattern(30000.0, '*v', 'X', 'P', J1, J2, A1, A2, 0, 0,
                             '', '', dwn)
    A, u, n = C.components_value([p], 'X', lambda lid: (A2, 0.0))
    assert A == pytest.approx(A1, abs=1e-9) and n == 2
    _, u_p, _ = C.components_value([p], 'X', lambda lid: (A2, 0.002))
    assert u_p > u
    assert C.components_value([p], 'X', lambda lid: None) is None
    assert C.components_value([p], 'Y', lambda lid: (A2, 0.0)) is None


def test_write_table_adds_replaces_and_keeps_the_measurement(tmp_path):
    path = tmp_path / 'A_hfs_levels.csv'
    path.write_text(
        'level_id,cfg,J,A_cm-1,u_A,source,n_flagged\n'
        '059003.000002,f5d2,3.5,+0.0380,0.0081,"resolved components, 1",1\n'
        '059003.000003,f5d6s,4.5,+0.0000,0.0500,not determined,1\n'
        '059003.000004,f26p,2.5,+0.0250,0.0050,composition,0\n',
        encoding='utf-8', newline='\n')

    def row(lid, action, A, u, src):
        return {'level_id': lid, 'action': action, 'J': '2.5',
                'proposed_cfg': 'f25d', 'proposed_A': A, 'proposed_u_A': u,
                'proposed_source': src, 'n_flagged': 0}
    calc = 'calculated, hfs_A_theory 2026-10-07, tier 1'
    out = [row('059003.000001', 'add', '+0.0300', '0.0020', calc),
           row('059003.000002', 'replace', '+0.0243', '0.0025', calc),
           row('059003.000003', 'replace', '+0.0806', '0.0115', calc),
           row('059003.000004', 'keep: semiempirical prior', '+0.03',
               '0.003', calc),
           row('059003.000005', 'add not determined', '+0.0000', '0.0100',
               'not determined')]
    assert C.write_table(out, str(path)) == (2, 2, 1)
    raw = path.read_bytes()
    assert b'\r' not in raw
    rows = {r['level_id']: r for r in csv.DictReader(
        raw.decode('utf-8').splitlines())}
    assert list(rows) == sorted(rows)
    assert rows['059003.000002']['A_cm-1'] == '+0.0243'
    assert rows['059003.000002']['cfg'] == 'f5d2'          # kept
    assert rows['059003.000004']['source'] == 'composition'
    assert rows['059003.000005']['source'] == 'not determined'
    with open(tmp_path / hfs_A_theory.SUPERSEDED_NAME, encoding='utf-8') as fh:
        old = list(csv.DictReader(fh))
    assert [r['level_id'] for r in old] == ['059003.000002']
    assert old[0]['A_cm-1'] == '+0.0380'
    assert old[0]['reason'].startswith('replaced')
    # the parameters' fit still sees the measurement
    assert hfs_A_theory.read_measured(str(path))['059003.000002'][0] == 0.038
