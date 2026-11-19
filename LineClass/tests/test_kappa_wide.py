"""Tests for the split of the plain 1974 class by the width of the pattern
(2026-11-19): a component of a plain 1974 line whose |D| exceeds [hfs]
kappa_wide_D takes plain_1974_wide; the narrower ones, and every 1969 line,
keep the plain class of their era.

Run from the LineClass directory:  python -m pytest tests -q

The tests fix:

  * the configuration: the class and the threshold go together;
  * the model: the class a component takes, in a single line and in a
    blend, and no split without the class;
  * the classification's factor for a candidate (models.SpectralLine);
  * the calibration: the class a line votes on, and its column;
  * level_positions' kappa of a candidate.
"""
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import config                              # noqa: E402
import hfs_correction as H                 # noqa: E402
import hfs_kappa                           # noqa: E402
import models                              # noqa: E402
import wavelength_calibration as W         # noqa: E402

KAPPA = (('flag', (1.0, 0.0)), ('plain_1974', (0.24, 0.04)),
         ('plain_1969', (0.15, 0.03)), ('c', (0.04, 0.04)),
         ('plain_1974_wide', (0.79, 0.06)))

LOW, WIDE, NARROW = '059003.000100', '059003.000300', '059003.000301'
L74, L69 = 30000.1234, 60000.5678
# D = I * (A_u J_u - A_l J_l), I = 5/2
D_WIDE = 2.5 * (0.10 * 2.5 - -0.04 * 3.5)       # 0.975
D_NARROW = 2.5 * (-0.02 * 1.5 - -0.04 * 3.5)    # 0.275


def a_table(tmp_path):
    path = tmp_path / 'A.csv'
    path.write_text(
        'level_id,cfg,J,A_cm-1,u_A,source,n_flagged\n'
        f'{LOW},f26s,3.5,-0.0400,0.0040,flag interval,0\n'
        f'{WIDE},f26p,2.5,+0.1000,0.0100,flag interval,0\n'
        f'{NARROW},f5d2,1.5,-0.0200,0.0030,flag interval,0\n',
        encoding='utf-8', newline='\n')
    return str(path)


def make_model(tmp_path, kappa=KAPPA, wide_D=0.4):
    return H.Model(config.HfsSettings(apply=True, A_levels=a_table(tmp_path),
                                      kappa=kappa, kappa_wide_D=wide_D))


# ---- the configuration -----------------------------------------------------

def baseline_copy(tmp_path, drop=()):
    """The baseline configuration as text, the lines that begin with one of
    `drop` left out (test_flag_unlisted.baseline_copy)."""
    with open(config.DEFAULT_PATH, encoding='utf-8') as fh:
        text = [ln for ln in fh if not ln.startswith(tuple(drop))]
    path = tmp_path / 'c.toml'
    path.write_text(''.join(text), encoding='utf-8', newline='\n')
    return str(path)


def test_the_baseline_splits_the_plain_1974_class():
    cfg = config.load()
    assert cfg.hfs.kappa_wide_D == pytest.approx(0.4)
    k, u = cfg.hfs.kappa_of('plain_1974_wide')
    assert k > cfg.hfs.kappa_of('plain_1974')[0] and 0.0 < u < 0.2


@pytest.mark.parametrize('drop', [('kappa_wide_D',), ('plain_1974_wide',)])
def test_the_class_and_the_threshold_go_together(tmp_path, drop):
    with pytest.raises(config.ConfigError, match='kappa_wide_D'):
        config.load(baseline_copy(tmp_path, drop=drop))


def test_without_both_there_is_no_split(tmp_path):
    cfg = config.load(baseline_copy(
        tmp_path, drop=('kappa_wide_D', 'plain_1974_wide')))
    assert cfg.hfs.kappa_wide_D == 0.0
    assert 'plain_1974_wide' not in dict(cfg.hfs.kappa)


# ---- the model -------------------------------------------------------------

def test_the_class_follows_the_width_of_the_pattern(tmp_path):
    m = make_model(tmp_path)
    assert m.D(LOW, WIDE)[0] == pytest.approx(D_WIDE)
    assert m.D(LOW, NARROW)[0] == pytest.approx(D_NARROW)
    assert m.plain_class(L74, D_WIDE) == 'plain_1974_wide'
    assert m.plain_class(L74, -D_WIDE) == 'plain_1974_wide'
    assert m.plain_class(L74, D_NARROW) == 'plain_1974'
    assert m.plain_class(L74, 0.4) == 'plain_1974'          # beyond, not at
    assert m.plain_class(L74) == 'plain_1974'               # no D given
    assert m.plain_class(L69, D_WIDE) == 'plain_1969'       # 1969 not split
    assert m.kappa_of('', L74, D_WIDE) == (0.79, 0.06)
    assert m.kappa_of('w', L74, D_WIDE) == (0.79, 0.06)     # any character
    assert m.kappa_of('c', L74, D_WIDE) == (0.04, 0.04)     # c is its own
    assert m.kappa_of('*v', L74, D_WIDE) == (1.0, 0.0)      # the anchor


def test_no_split_without_the_class(tmp_path):
    m = make_model(tmp_path, kappa=KAPPA[:4], wide_D=0.4)
    assert m.wide_D is None
    assert m.plain_class(L74, D_WIDE) == 'plain_1974'
    assert make_model(tmp_path, wide_D=0.0).wide_D is None


def test_a_single_line_is_shifted_by_its_own_class(tmp_path):
    m = make_model(tmp_path)
    shift, _, kappa = m.line_shift('', L74, [(LOW, WIDE, 1.0)])
    assert kappa == pytest.approx(0.79)
    assert shift == pytest.approx((1 - 0.79) * D_WIDE)
    shift, _, _ = m.line_shift('', L74, [(LOW, NARROW, 1.0)])
    assert shift == pytest.approx((1 - 0.24) * D_NARROW)
    cg = m.cg_shift('', L74, [(LOW, WIDE, 1.0)])
    assert cg.shift == pytest.approx(-0.79 * D_WIDE)


def test_each_component_of_a_blend_takes_its_own_class(tmp_path):
    m = make_model(tmp_path)
    parts = m.shift_parts('', L74, [(LOW, WIDE, 0.6), (LOW, NARROW, 0.4)])
    assert [p[1] for p in parts] == [0.79, 0.24]
    # in a flagged blend the components other than the strongest take the
    # plain class their own D chooses
    parts = m.shift_parts('*v', L74, [(LOW, NARROW, 0.6), (LOW, WIDE, 0.4)])
    assert [p[1] for p in parts] == [1.0, 0.79]


# ---- the classification ----------------------------------------------------

def test_the_candidate_factor_follows_the_width_of_the_pattern():
    line = models.SpectralLine(wavenumber=L74, wn_uncertainty=0.01,
                               intensity=1.0, line_character='', wn_key=L74)
    line.hfs_factor, line.hfs_factor_wide, line.hfs_wide_D = 0.76, 0.21, 0.4
    assert line.hfs_factor_for(0.975) == 0.21
    assert line.hfs_factor_for(-0.975) == 0.21
    assert line.hfs_factor_for(0.275) == 0.76
    line.hfs_factor_wide = None
    assert line.hfs_factor_for(0.975) == 0.76


# ---- the calibration -------------------------------------------------------

def line(cls='plain', era=1974, D=0.5, low=LOW, upp=WIDE):
    return hfs_kappa.Line(wn=30000.0 if era == 1974 else 60000.0, char='',
                          era=era, cls=cls, ucls=('', era), low=low,
                          upp=upp, D=D, dJ=None)


def test_the_class_a_line_votes_on(monkeypatch):
    monkeypatch.setattr(W, 'WIDE_D', 0.4)
    assert W.kappa_name(line(D=0.5)) == 'plain_1974_wide'
    assert W.kappa_name(line(D=-0.5)) == 'plain_1974_wide'
    assert W.kappa_name(line(D=0.3)) == 'plain_1974'
    assert W.kappa_name(line(D=0.9, era=1969)) == 'plain_1969'
    assert W.kappa_name(line('c', D=0.9)) == 'c'
    monkeypatch.setattr(W, 'WIDE_D', None)
    assert W.kappa_name(line(D=0.5)) == 'plain_1974'


def build(lines, fixed=None):
    n = len(lines)
    return W.build(lines, [3000.0] * n, [0] * n, [0] * n, [0.01] * n,
                   [True] * n, 0, lambda *a: [], fixed)


def test_a_column_for_the_wide_class_only_where_it_has_lines(monkeypatch):
    monkeypatch.setattr(W, 'WIDE_D', 0.4)
    assert 'plain_1974_wide' not in build([line(D=0.2, upp=NARROW)])[4]
    rows, levels, lidx, off, kidx, A, y, w = build(
        [line(D=0.2, upp=NARROW), line(D=0.7)])
    assert A[1, kidx['plain_1974_wide']] == pytest.approx(0.7)
    assert A[0, kidx['plain_1974_wide']] == 0.0
    assert A[0, kidx['plain_1974']] == pytest.approx(0.2)
    # held, the class takes its kappa * D off the right-hand side
    parts = build([line(D=0.7)], {'plain_1974_wide': 0.79})
    assert 'plain_1974_wide' not in parts[4]
    assert parts[6][0] == pytest.approx(30000.0 - 0.79 * 0.7)


# ---- level_positions -------------------------------------------------------

def test_level_positions_takes_the_wide_kappa_for_a_wide_candidate():
    import level_positions as LP

    class Ctx:
        pass
    ctx = Ctx()
    ctx.hfs = object()
    ctx.unc_o = np.array([0.01])
    ctx.wn_o = np.array([L74])
    ctx.kap_o, ctx.ukap_o = np.array([0.24]), np.array([0.04])
    ctx.kapw_o, ctx.ukapw_o = np.array([0.79]), np.array([0.06])
    ctx.wide_D = 0.4
    ctx.head_o = np.array([False])
    ctx.spread_o = np.array([0.0])
    z = np.zeros(2)
    tt = {'D': np.array([0.975, 0.275]), 'w2a': z, 'w2u': z, 'u2u': z,
          'u2k': z, 'acc_i': np.array([-1, -1]), 'acc_wn': z, 'acc_unc': z}
    nu = np.array([[L74], [L74]])
    var, d = LP.candidate_terms(ctx, tt, np.array([0]), nu)
    assert d[0, 0] == pytest.approx((1 - 0.79) * 0.975)
    assert d[1, 0] == pytest.approx((1 - 0.24) * 0.275)
    assert var[0, 0] == pytest.approx(0.01 ** 2 + (0.06 * 0.975) ** 2)
    assert var[1, 0] == pytest.approx(0.01 ** 2 + (0.04 * 0.275) ** 2)
