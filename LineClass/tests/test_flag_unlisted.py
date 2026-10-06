"""Tests for the split of the flagged lines (2026-10-06): a flagged line whose
resolved components Sugar printed (files.hfs_components) is the anchor,
kappa = 1; one whose components he did not print is of the class
flag_unlisted, which the plate calibration fits.

Run from the LineClass directory:  python -m pytest tests -q

The tests fix:

  * the class of a line, and the join of its wn_key to the list;
  * the configuration: the list needs [hfs.kappa] flag_unlisted, and a
    configuration without the list keeps every flagged line the anchor's;
  * the model: the kappa, the blend rule, the head-frame test and the terms
    the calibration sees;
  * the calibration: a column for flag_unlisted only where it has lines, and
    the writer of [hfs.kappa], which writes it only where it is given.
"""
import os
import sys
import tomllib

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import config                              # noqa: E402
import hfs_correction as H                 # noqa: E402
import hfs_kappa                           # noqa: E402
import wavelength_calibration as W         # noqa: E402

KAPPA = (('flag', (1.0, 0.0)), ('plain_1974', (0.81, 0.013)),
         ('plain_1969', (0.63, 0.02)), ('c', (0.16, 0.05)),
         ('flag_unlisted', (0.90, 0.024)))

LOW, UPP, OTHER = '059003.000100', '059003.000300', '059003.000301'
LISTED = 30000.1234          # a flagged line whose components are printed
UNLISTED = 31000.5678        # one whose components are not


def a_table(tmp_path):
    path = tmp_path / 'A.csv'
    path.write_text(
        'level_id,cfg,J,A_cm-1,u_A,source,n_flagged\n'
        f'{LOW},f26s,3.5,-0.0400,0.0040,flag interval,0\n'
        f'{UPP},f26p,2.5,+0.1000,0.0100,flag interval,0\n'
        f'{OTHER},f5d2,1.5,+0.0200,0.0030,flag interval,0\n',
        encoding='utf-8', newline='\n')
    return str(path)


def components(tmp_path):
    path = tmp_path / 'components.csv'
    path.write_text(
        'wn_line,char,wn_component,dwn\n'
        f'{LISTED + 0.004:.6f},*v,{LISTED + 0.6:.6f},+0.5960\n'
        f'{LISTED + 0.004:.6f},*v,{LISTED + 1.1:.6f},+1.0960\n',
        encoding='utf-8', newline='\n')
    return str(path)


def make_model(tmp_path, with_list=True, kappa=KAPPA):
    return H.Model(config.HfsSettings(
        apply=True, A_levels=a_table(tmp_path), kappa=kappa,
        components=components(tmp_path) if with_list else ''))


D = 2.5 * (0.1 * 2.5 - -0.04 * 3.5)


# ---- the class -------------------------------------------------------------

def test_the_class_of_a_line():
    assert hfs_kappa.kappa_class('*v') == 'flag'
    assert hfs_kappa.kappa_class('*r', listed=False) == 'flag_unlisted'
    # only the flags are split
    assert hfs_kappa.kappa_class('c', listed=False) == 'c'
    assert hfs_kappa.kappa_class('', listed=False) == 'plain'


def test_the_join_to_the_list(tmp_path):
    listed = hfs_kappa.read_listed(components(tmp_path))
    assert len(listed) == 1                # one line, two components
    assert LISTED in listed                # 0.004 away: the same line
    assert LISTED + 0.004 + 0.007 not in listed
    assert UNLISTED not in listed


def test_the_project_list_reads():
    listed = hfs_kappa.read_listed(os.path.join(ROOT, 'hfs_components.csv'))
    assert 200 < len(listed) < 300
    assert 29602.34223667826 in listed     # 1098-985 *v, four components


# ---- the configuration -----------------------------------------------------

def baseline_copy(tmp_path, drop=()):
    """The baseline configuration, copied as text with the lines that begin
    with one of `drop` left out.  Its paths stay relative, which loading does
    not check."""
    with open(config.DEFAULT_PATH, encoding='utf-8') as fh:
        text = [ln for ln in fh if not ln.startswith(tuple(drop))]
    path = tmp_path / 'c.toml'
    path.write_text(''.join(text), encoding='utf-8', newline='\n')
    return str(path)


def test_the_list_needs_the_class(tmp_path):
    with pytest.raises(config.ConfigError, match='flag_unlisted'):
        config.load(baseline_copy(tmp_path, drop=('flag_unlisted',)))


def test_without_the_list_the_class_is_not_needed(tmp_path):
    cfg = config.load(baseline_copy(
        tmp_path, drop=('flag_unlisted', 'hfs_components')))
    assert cfg.hfs.components == ''
    assert 'flag_unlisted' not in dict(cfg.hfs.kappa)


def test_the_baseline_splits_the_flags():
    cfg = config.load()
    assert cfg.hfs.components == os.path.join(ROOT, 'hfs_components.csv')
    k, u = cfg.hfs.kappa_of('flag_unlisted')
    assert 0.0 < k < 1.0 and 0.0 < u < 0.2


# ---- the model -------------------------------------------------------------

def test_the_kappa_of_a_line(tmp_path):
    m = make_model(tmp_path)
    assert m.line_class('*v', LISTED) == 'flag'
    assert m.line_class('*v', UNLISTED) == 'flag_unlisted'
    assert m.kappa_of('*v', LISTED) == (1.0, 0.0)
    assert m.kappa_of('*r', UNLISTED) == (0.90, 0.024)
    assert m.kappa_of('', UNLISTED) == (0.81, 0.013)


def test_without_the_list_every_flag_is_the_anchor(tmp_path):
    m = make_model(tmp_path, with_list=False, kappa=KAPPA[:4])
    assert m.listed is None
    assert m.line_class('*v', UNLISTED) == 'flag'
    assert m.kappa_of('*v', UNLISTED) == (1.0, 0.0)


def test_a_list_without_the_class_stops(tmp_path):
    with pytest.raises(ValueError, match='flag_unlisted'):
        make_model(tmp_path, kappa=KAPPA[:4])


def test_the_shift_of_an_unlisted_flag(tmp_path):
    m = make_model(tmp_path)
    shift, u, kappa = m.line_shift('*v', UNLISTED, [(LOW, UPP, 1.0)])
    assert kappa == pytest.approx(0.90)
    assert shift == pytest.approx(0.10 * D)
    assert m.line_shift('*v', LISTED, [(LOW, UPP, 1.0)])[0] == 0.0


def test_an_unlisted_flagged_blend(tmp_path):
    """The strongest component takes the line's own kappa, the others the
    plain kappa of its era, as in a flagged blend of the anchor class."""
    m = make_model(tmp_path)
    got = m.component_kappas('*v', UNLISTED, [0.3, 0.7])
    assert got == [(0.81, 0.013), (0.90, 0.024)]
    got = m.component_kappas('*v', LISTED, [0.3, 0.7])
    assert got == [(0.81, 0.013), (1.0, 0.0)]


def test_only_a_listed_flag_sits_at_the_head(tmp_path):
    m = make_model(tmp_path)
    assert m.sits_at_head('*v', LISTED)
    assert not m.sits_at_head('*v', UNLISTED)
    assert not m.sits_at_head('', LISTED)


def test_the_terms_the_calibration_sees(tmp_path):
    m = make_model(tmp_path)
    assert m.fit_terms('*v', UNLISTED, LOW, UPP) == ('flag_unlisted',
                                                     pytest.approx(D), '')
    assert m.fit_terms('*r', LISTED, LOW, UPP)[0] == 'flag'


# ---- the calibration -------------------------------------------------------

def line(cls, era=1974, low=LOW, upp=UPP, D=0.5):
    return hfs_kappa.Line(wn=30000.0 if era == 1974 else 50000.0, char='',
                          era=era, cls=cls, ucls=('', era), low=low,
                          upp=upp, D=D, dJ=None)


def test_the_kappa_a_line_votes_on():
    assert W.kappa_name(line('plain')) == 'plain_1974'
    assert W.kappa_name(line('plain', era=1969)) == 'plain_1969'
    assert W.kappa_name(line('c')) == 'c'
    assert W.kappa_name(line('flag_unlisted')) == 'flag_unlisted'


def build(lines):
    n = len(lines)
    return W.build(lines, [3000.0] * n, [0] * n, [0] * n, [0.01] * n,
                   [True] * n, 0, lambda *a: [])


def test_a_column_for_flag_unlisted_only_where_it_has_lines():
    lines = [line('flag'), line('plain', low=OTHER), line('c', upp=OTHER)]
    kidx = build(lines)[4]
    assert list(kidx) == ['plain_1974', 'plain_1969', 'c']
    lines.append(line('flag_unlisted', D=0.7))
    rows, levels, lidx, off, kidx, A, y, w = build(lines)
    assert list(kidx) == ['plain_1974', 'plain_1969', 'c', 'flag_unlisted']
    assert A[3, kidx['flag_unlisted']] == pytest.approx(0.7)
    assert np.all(A[:3, kidx['flag_unlisted']] == 0.0)
    # a line of the anchor class is no vote: its D is on the right-hand side
    assert y[0] == pytest.approx(30000.0 - 0.5)


PARENT = """\
[hfs]
apply = false

[hfs.kappa]
flag          = [1.0, 0.0]
plain_1974    = [0.926, 0.012]
plain_1969    = [0.717, 0.021]
c             = [0.276, 0.044]
flag_unlisted = [0.900, 0.030]
"""

NEW = {'plain_1974': (0.81, 0.013), 'plain_1969': (0.62, 0.02),
       'c': (0.16, 0.05)}


def test_the_writer_writes_flag_unlisted_when_given(tmp_path):
    path = tmp_path / 'k.toml'
    path.write_text(PARENT, encoding='utf-8', newline='\n')
    config.set_hfs_kappa(str(path), NEW)
    with open(path, 'rb') as fh:
        got = tomllib.load(fh)['hfs']['kappa']
    assert got['flag_unlisted'] == [0.9, 0.03]      # not given: left alone
    config.set_hfs_kappa(str(path), dict(NEW, flag_unlisted=(0.895, 0.024)))
    text = path.read_text(encoding='utf-8')
    assert 'flag_unlisted = [0.895, 0.024]' in text
    assert 'plain_1974    = [0.810, 0.013]' in text     # alignment kept


def test_the_writer_needs_the_line_of_a_class_it_is_given(tmp_path):
    path = tmp_path / 'k.toml'
    path.write_text(PARENT.replace('flag_unlisted = [0.900, 0.030]\n', ''),
                    encoding='utf-8', newline='\n')
    with pytest.raises(config.ConfigError, match='flag_unlisted'):
        config.set_hfs_kappa(str(path), dict(NEW, flag_unlisted=(0.9, 0.02)))
    with pytest.raises(config.ConfigError, match='plain_1969'):
        config.set_hfs_kappa(str(path), {'plain_1974': (0.8, 0.01),
                                         'c': (0.1, 0.05)})
