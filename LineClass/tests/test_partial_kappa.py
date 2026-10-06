"""Tests for the partial classes of the kappa exception registry: a line
measured on a part of its hfs pattern (hfs_patterns.partial_pattern,
hfs_correction.Model.row_kappa) and the search for them in
kappa_candidates.py, which is information only.

Run from the LineClass directory:  python -m pytest tests -q

The pattern of 1116-1021 (J 3.5 -> 3.5, A -0.0377 and +0.0242) is the worked
case of 2026-10-05: 30074.3893 sits at kappa -0.42, where the pattern without
its head (kappa -0.36) would be measured.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import classify_lines as cl                # noqa: E402
import hfs_correction as H                 # noqa: E402
import hfs_patterns as P                   # noqa: E402
import kappa_candidates as K               # noqa: E402
from models import SpectralLine            # noqa: E402
from test_kappa_exceptions import (        # noqa: E402
    LOW, UPP, RESOLVED, OTHER, S_LOW, S_UPP, level, make_model,
    switched_on)

PATTERN = (3.5, -0.03768458, 3.5, 0.02419138)


# ---------------------------------------------------------------------------
# The pattern
# ---------------------------------------------------------------------------
def test_the_pattern_without_its_head():
    D = P.head_displacement(*PATTERN)
    assert D == pytest.approx(0.5414, abs=1e-4)
    part = P.partial_pattern(*PATTERN, 1)
    assert part.kappa == pytest.approx(-0.3622, abs=1e-4)
    assert part.share == pytest.approx(0.7292, abs=1e-4)
    side, where, share = part.missing
    assert side == 'head' and share == pytest.approx(0.2708, abs=1e-4)
    # the kept part and the missing one balance about the cg
    assert part.kappa * D * part.share + where * share == pytest.approx(0.0)


def test_the_head_and_rung_1_only():
    part = P.partial_pattern(*PATTERN, 0, 1)
    assert part.kappa == pytest.approx(0.6687, abs=1e-4)
    assert part.missing[0] == 'tail'
    assert part.share + part.missing[2] == pytest.approx(1.0)


def test_every_component_joins_one_rung():
    comps = P.rung_of(*PATTERN)
    assert len(comps) == len(P.components(*PATTERN))
    assert {k for _, _, k in comps} == set(range(6))
    # the head and its off-ladder neighbor F 6-5 (kappa 0.73) are rung 0
    D = P.head_displacement(*PATTERN)
    assert sorted(round(p / D, 2) for p, _, k in comps if k == 0) \
        == [0.73, 1.0]


@pytest.mark.parametrize('first, last, message', [
    (6, None, 'outside the ladder'),
    (2, 1, 'outside the ladder'),
    (0, 5, 'the whole pattern'),
    (0, None, 'the whole pattern'),
])
def test_a_part_that_is_not_one(first, last, message):
    with pytest.raises(ValueError, match=message):
        P.partial_pattern(*PATTERN, first, last)


def test_a_pattern_with_no_displacement():
    with pytest.raises(ValueError, match='D = 0'):
        P.partial_pattern(1.5, 0.0, 1.5, 0.0, 1)


# ---------------------------------------------------------------------------
# The registry and the model
# ---------------------------------------------------------------------------
def test_the_names():
    assert H.parse_partial('partial:1-') == (1, None)
    assert H.parse_partial('partial:0-2') == (0, 2)
    assert H.parse_partial('cg') is None
    assert H.parse_partial('partial:x') is None
    assert H.partial_name(1) == 'partial:1-'
    assert H.partial_name(0, 2) == 'partial:0-2'


def test_a_partial_row_takes_its_patterns_kappa(tmp_path):
    m = make_model(tmp_path, ('30000.1234', LOW, UPP, 'partial:1-', '',
                              '10/5/2026', 'no head'))
    row = m.exception(30000.1234, LOW, UPP)
    J1, A1, _, J2, A2, _ = m.pattern_of(LOW, UPP)
    want = P.partial_pattern(J1, A1, J2, A2, 1).kappa
    kappa, u = m.row_kappa(row)
    assert kappa == pytest.approx(want)
    assert 0.0 < u < 0.05             # from the two A uncertainties
    with pytest.raises(ValueError, match='depends on the transition'):
        m.kappa_value('partial:1-')
    D = S_UPP - S_LOW
    shift, _, k = m.line_shift('w', 30000.1234, [(LOW, UPP, 1.0)])
    assert k == pytest.approx(want)
    assert shift == pytest.approx((1.0 - want) * D)
    # the calibration holds it: an anchor-class line with D times kappa
    cls, d, fixed = m.fit_terms('w', 30000.1234, LOW, UPP)
    assert (cls, fixed) == ('flag', 'partial:1-')
    assert d == pytest.approx(want * D)
    assert not m.is_head_line(30000.1234)


@pytest.mark.parametrize('row, message', [
    (('30000.5', LOW, UPP, 'partial:6-', '', '', ''), 'outside the ladder'),
    (('30000.5', LOW, UPP, 'partial:0-5', '', '', ''), 'whole pattern'),
    (('30000.5', RESOLVED, UPP, 'partial:1-', '', '', ''),
     'resolved_levels'),
])
def test_a_partial_row_that_cannot_be(tmp_path, row, message):
    with pytest.raises(ValueError, match=message):
        make_model(tmp_path, row)


def test_a_resolved_level_named_unresolved_has_its_pattern(tmp_path):
    m = make_model(tmp_path, ('31000.5', RESOLVED, OTHER, 'partial:1-',
                              RESOLVED, '', ''))
    kappa = m.row_kappa(m.exception(31000.5, RESOLVED, OTHER))[0]
    J1, A1, _, J2, A2, _ = m.pattern_of(RESOLVED, OTHER, {RESOLVED})
    assert kappa == pytest.approx(P.partial_pattern(J1, A1, J2, A2, 1).kappa)


def test_a_part_with_no_strength():
    # J 0.5 -> 2.5 has no allowed component at all
    with pytest.raises(ValueError, match='no strength'):
        P.partial_pattern(0.5, 0.5, 2.5, 0.3, 1)


def test_the_classification_attaches_the_partial_kappa(monkeypatch,
                                                       tmp_path):
    m = make_model(tmp_path, ('30000.1234', LOW, UPP, 'partial:0-1', '',
                              '', ''))
    switched_on(monkeypatch, m, str(tmp_path / 'exc.txt'))
    levels = {x: level(x, 0.0) for x in (LOW, UPP, OTHER, RESOLVED)}
    line = SpectralLine(30000.1234, 0.1, 1.0, '', wn_key=30000.1234)
    assert cl.attach_kappa_exceptions([line], levels) == 1
    kappa, u, unresolved = line.hfs_exceptions[(LOW, UPP)]
    assert (kappa, u) == m.partial_kappa(LOW, UPP, 0, 1)
    assert unresolved == frozenset()


# ---------------------------------------------------------------------------
# The search (kappa_candidates.py)
# ---------------------------------------------------------------------------
def test_the_family_and_the_best_part(tmp_path):
    m = make_model(tmp_path)
    fam = K.partial_family(m, LOW, UPP)
    names = [n for n, _, _ in fam]
    # ladder 0-5: the head side 1- ... 5-, the tail side 0-1 ... 0-4
    assert names == ['partial:%d-' % k for k in range(1, 6)] + \
        ['partial:0-%d' % k for k in range(1, 5)]
    assert K.partial_family(m, RESOLVED, UPP) == []
    # a line measured exactly where 'partial:2-' puts it
    D = S_UPP - S_LOW
    part = dict((n, p) for n, p, _ in fam)['partial:2-']
    r_loo = (part.kappa - 0.944) * D          # fitted at 0.944
    name, best, z = K.best_partial(fam, r_loo, 0.944, D, 0.05)
    assert name == 'partial:2-' and z == pytest.approx(0.0, abs=1e-9)


def test_where_the_missing_part_lies():
    part = P.partial_pattern(*PATTERN, 1)
    D = P.head_displacement(*PATTERN)
    # 30074.3893: measured at kappa -0.42; the head lies 0.758 above it
    assert K.missing_offset(part, -0.42, D) == pytest.approx(
        part.missing[1] + 0.42 * D)
    assert K.missing_offset(part, -0.42, D) == pytest.approx(0.755, abs=0.005)


def test_the_line_itself_is_never_the_missing_part():
    keys = [100.0, 100.02, 101.0]
    assert K.line_near(keys, 100.0, 0.1) == 100.0
    assert K.line_near(keys, 100.0, 0.1, exclude=100.0) == 100.02
    assert K.line_near([100.0, 101.0], 100.0, 0.1, exclude=100.0) is None
