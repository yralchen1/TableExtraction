"""Tests for hfs_A_theory.py: hyperfine A constants from the RCEOUT
eigenvectors (2026-10-06).

Run from the LineClass directory:  python -m pytest tests -q

The tests fix:

  * the one-electron operator: one p, d or f electron gives
    A_j = a l(l+1)/[j(j+1)] with a01 = a12 = a, and one s electron A = a10;
  * the many-electron states: a configuration's J blocks hold every
    determinant once, and a pure LS state of 4f2 6s gives the vector-model
    projections of its 4f orbital and 6s contact terms;
  * RCG's label rule, and that every component RCEOUT prints for the
    configurations computed here names exactly one basis state;
  * the parameter file: ties, configuration-specific entries, missing kinds;
  * the fit: synthetic constants give back the parameters they were made
    with.

Everything that builds states needs the Cowan repository (COWAN_REPO); those
tests are skipped where it is not.
"""
import math
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import hfs_A_theory as H                   # noqa: E402

needs_cowan = pytest.mark.skipif(
    not os.path.isfile(os.path.join(H.COWAN_REPO, 'CODE', 'ING11.CFP')),
    reason='the Cowan repository is not at COWAN_REPO')
needs_rceout = pytest.mark.skipif(not os.path.isfile(H.RCEOUT_FILE),
                                  reason='no RCEOUT')


# ---------------------------------------------------------------------------
# one electron
# ---------------------------------------------------------------------------
@needs_cowan
@pytest.mark.parametrize('l', [1, 2, 3])
def test_one_electron_A_j(l):
    basis = H.ConfigBasis('x', card=[('nl', l, 1)])
    for J2 in (2 * l - 1, 2 * l + 1):
        m = basis.matrices(J2)
        assert m[('nl', 'a01')].shape == (1, 1)
        j = J2 / 2.0
        a = m[('nl', 'a01')][0, 0] + m[('nl', 'a12')][0, 0]
        assert a == pytest.approx(l * (l + 1) / (j * (j + 1)), abs=1e-10)


@needs_cowan
def test_one_s_electron_A_is_a10():
    basis = H.ConfigBasis('x', card=[('6s', 0, 1)])
    m = basis.matrices(1)
    assert set(m) == {('6s', 'a10')}
    assert m[('6s', 'a10')][0, 0] == pytest.approx(1.0, abs=1e-12)


# ---------------------------------------------------------------------------
# many electrons
# ---------------------------------------------------------------------------
@needs_cowan
@pytest.mark.parametrize('name, ndet', [('f26s', 91 * 2), ('f25d', 91 * 10),
                                        ('4f3', 364), ('f5d6s', 14 * 10 * 2)])
def test_blocks_hold_every_determinant_once(name, ndet):
    basis = H.config_basis(name)
    total = 0
    for J2 in range(0, 40):
        chains, labels = basis.block(J2)
        total += (J2 + 1) * len(chains)
    assert total == ndet


def _vector_model(J, L, S, S1):
    """The 6s contact and the 4f orbital projection of (S1 L) 6s; S L J."""
    jj = J * (J + 1)
    s_on_S = (S * (S + 1) + 0.75 - S1 * (S1 + 1)) / (2 * S * (S + 1))
    S_on_J = (jj + S * (S + 1) - L * (L + 1)) / (2 * jj)
    L_on_J = (jj + L * (L + 1) - S * (S + 1)) / (2 * jj)
    return s_on_S * S_on_J, L_on_J


@needs_cowan
def test_pure_LS_states_of_4f2_6s_follow_the_vector_model():
    basis = H.config_basis('f26s')
    checked = 0
    for J2 in range(1, 16, 2):
        chains, labels = basis.block(J2)
        m = basis.matrices(J2)
        for i, ch in enumerate(chains):
            L = ch[-1][1] / 2.0
            S = ch[-1][2] / 2.0
            S1 = ch[0][2] / 2.0
            contact, orbit = _vector_model(J2 / 2.0, L, S, S1)
            assert m[('6s', 'a10')][i, i] == pytest.approx(contact, abs=1e-10)
            assert m[('4f', 'a01')][i, i] == pytest.approx(orbit, abs=1e-10)
            checked += 1
    assert checked > 20


@needs_cowan
def test_the_4P_2P_contact_cross_term_of_J_one_half():
    """The term that decides 000190 and 000187: the 6s contact operator
    between (3P) 6s 4P1/2 and 2P1/2.  Total spin S is diagonal, so the
    spin of the two 4f electrons must cancel the cross term of the 6s
    spin exactly, and their sum on the diagonal is the projection of S on
    J: 5/3 for 4P1/2, -1/3 for 2P1/2.  The sign of the cross term is
    Cowan's convention; the fit of a_6s to the measured constants is what
    tests it."""
    basis = H.config_basis('f26s')
    chains, labels = basis.block(1)
    m = basis.matrices(1)
    s6, s4 = m[('6s', 'a10')], m[('4f', 'a10')]
    i, j = labels.index('3P4P'), labels.index('3P2P')
    assert abs(s6[i, j]) == pytest.approx(4.0 * math.sqrt(2.0) / 9.0,
                                          abs=1e-10)
    assert s4[i, j] == pytest.approx(-s6[i, j], abs=1e-12)
    assert s4[i, i] + s6[i, i] == pytest.approx(5.0 / 3.0, abs=1e-10)
    assert s4[j, j] + s6[j, j] == pytest.approx(-1.0 / 3.0, abs=1e-10)


@needs_cowan
def test_matrices_are_symmetric():
    basis = H.config_basis('f5d2')
    for m in basis.matrices(7).values():
        assert np.allclose(m, m.T, atol=1e-12)


# ---------------------------------------------------------------------------
# labels
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('name, isubj', [
    ('f26s', -1),      # 4f2 is the parent
    ('f5d2', -3),      # 5d2
    ('4f3', -1),
    ('f5d6s', 3),      # the running term after 4f 5d
    ('fd6p', 3),
    ('6s26p', 1),      # nothing marked: subshell 1, 4f0, i.e. 1S
    ('f6s2', 1),
])
def test_rcg_label_position(name, isubj):
    assert H.rcg_label_position(H.CONFIGS[name]) == isubj


@needs_cowan
def test_labels_of_4f2_6s():
    basis = H.config_basis('f26s')
    assert set(basis.labels.values()) == {
        '1S2S', '3P2P', '3P4P', '1D2D', '3F2F', '3F4F', '1G2G', '3H2H',
        '3H4H', '1I2I'}


@needs_cowan
def test_labels_carry_the_alpha_of_4f3():
    labels = set(H.config_basis('4f3').labels.values())
    assert {'2H2H1', '2H2H2', '2D2D1', '2D2D2', '4I4I'} <= labels


@needs_cowan
def test_norm_label():
    assert H.norm_label('(3F) 4D  ') == '3F4D'
    assert H.norm_label('(2H) 2H2') == '2H2H2'


@needs_cowan
@needs_rceout
def test_every_printed_component_names_one_state():
    levels = H.read_levels()
    assert len(levels) > 5000
    unmatched = {u for lv in levels for u in lv.unmatched}
    assert not unmatched
    skipped = [lv.w_skipped for lv in levels if lv.e_obs is not None]
    assert max(skipped) < 0.05


# ---------------------------------------------------------------------------
# parameters
# ---------------------------------------------------------------------------
def _params():
    return H.Params({
        '4f': {'a01': 0.03, 'a12': 'a01', 'a10': 0.0},
        '4f3/4f': {'a01': 0.025, 'a12': 'a01', 'a10': 0.0},
        '6s': {'a10': 0.6},
        '5d': {'a12': 'a01'},
    })


def test_params_ties_and_overrides():
    p = _params()
    assert p.value('f26s', '4f', 'a12') == 0.03
    assert p.value('4f3', '4f', 'a01') == 0.025
    assert p.value('4f3', '4f', 'a12') == 0.025
    assert p.value('f26s', '6s', 'a10') == 0.6
    assert p.value('f26s', '6s', 'a01') is None
    assert p.value('f5d2', '5d', 'a12') is None       # tied to a missing one
    assert p.value('f5d2', '7s', 'a10') is None


def test_params_reject_unknown_kinds():
    with pytest.raises(ValueError):
        H.Params({'4f': {'b01': 1.0}})
    with pytest.raises(ValueError):
        H.Params({'4f': {'a12': 'b01'}})


def test_prepare_free_copies_the_plain_entry():
    p = _params()
    H.prepare_free(p, ['f5d2/4f.a01', '5d.a01'])
    assert p.table['f5d2/4f'] == {'a01': 0.03, 'a12': 'a01', 'a10': 0.0}
    assert p.value('f5d2', '5d', 'a12') == 0.0


# ---------------------------------------------------------------------------
# the fit
# ---------------------------------------------------------------------------
class _Level:
    """Enough of a LevelTheta for `fit`: thetas, weights, an A."""

    def __init__(self, n, theta, conf='f5d2'):
        self.level_id = 'L%03d' % n
        self.J = 2.5
        self.theta = {(conf, nl, kind): v for (nl, kind), v in theta.items()}
        self.weights = {conf: 1.0}

    def A(self, params):
        total = 0.0
        for (conf, nl, kind), th in self.theta.items():
            v = params.value(conf, nl, kind)
            if v is not None:
                total += th * v
        return total, [], 0.0

    def uncertainties(self, params, A, miss=None):
        return 0.0, 0.0, 1.0


def test_fit_recovers_the_parameters_it_was_made_with():
    rng = np.random.default_rng(3)
    truth = {'4f': 0.027, '5d.a01': 0.022, '5d.a10': -0.03}
    levels, measured = [], {}
    for n in range(40):
        th = {('4f', 'a01'): rng.uniform(0.3, 0.9),
              ('4f', 'a12'): rng.uniform(-0.1, 0.1),
              ('5d', 'a01'): rng.uniform(0.0, 0.6),
              ('5d', 'a12'): rng.uniform(-0.2, 0.2),
              ('5d', 'a10'): rng.uniform(-0.3, 0.3)}
        A = (truth['4f'] * (th[('4f', 'a01')] + th[('4f', 'a12')])
             + truth['5d.a01'] * (th[('5d', 'a01')] + th[('5d', 'a12')])
             + truth['5d.a10'] * th[('5d', 'a10')])
        lv = _Level(n, th)
        levels.append(lv)
        measured[lv.level_id] = (A, 0.001, 'synthetic')
    p = H.Params({'4f': {'a01': 0.03, 'a12': 'a01', 'a10': 0.0},
                  '5d': {'a12': 'a01'}})
    sol, sig, chi2, dof, used = H.fit(levels, p, ['4f.a01', '5d.a01',
                                                  '5d.a10'], measured)
    assert len(used) == 40 and dof == 37
    assert sol == pytest.approx([0.027, 0.022, -0.03], abs=1e-9)
    assert chi2 == pytest.approx(0.0, abs=1e-12)
    assert p.value('f5d2', '5d', 'a12') == pytest.approx(0.022, abs=1e-9)


def test_fit_leaves_out_levels_mostly_without_parameters():
    lv_ok = _Level(1, {('4f', 'a01'): 0.8})
    lv_bad = _Level(2, {('4f', 'a01'): 0.8}, conf='f25f')
    lv_bad.theta[('f25f', '5f', 'a01')] = 0.2
    lv_ok2 = _Level(3, {('4f', 'a01'): 0.5})
    measured = {lv.level_id: (0.02, 0.001, 's')
                for lv in (lv_ok, lv_bad, lv_ok2)}
    p = H.Params({'4f': {'a01': 0.03, 'a12': 'a01', 'a10': 0.0}})
    _sol, _sig, _chi2, _dof, used = H.fit([lv_ok, lv_bad, lv_ok2], p,
                                          ['4f.a01'], measured)
    assert [lv.level_id for lv in used] == ['L001', 'L003']
