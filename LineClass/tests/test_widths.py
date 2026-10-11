"""Tests for a set's own width model: the [widths] table of its configuration
(config.py) and level_positions.apply_widths.

Run from the LineClass directory:  python -m pytest tests -q

The tests fix:

  * that [widths] is optional and validated: eras, "era/code" keys, no
    negative width;
  * that apply_widths replaces exactly the constants the table names, gives
    each code one form, and turns the per-level hyperfine widths off;
  * the kappa scatter: (s(era) |D|)^2, and that it needs the hfs model;
  * that final/'s configuration carries its table and the others none.
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

import config                              # noqa: E402
import level_positions as lp               # noqa: E402


def write_config(tmp_path, widths_text, parent='lineclass_config.toml'):
    path = tmp_path / 'lineclass_config.toml'
    base = os.path.join(ROOT, parent).replace('\\', '/')
    path.write_text('inherit = "%s"\n\n[files]\noutput = "o.xlsx"\n'
                    'output_csv = "o.csv"\n\n%s' % (base, widths_text),
                    encoding='utf-8')
    return str(path)


@pytest.fixture
def fresh(monkeypatch):
    """level_positions' width constants, restored after the test."""
    for name in ('DLAM_ERA', 'CHAR_DLAM', 'CHAR_DWN', 'KAPPA_SCATTER',
                 '_D_CACHE'):
        monkeypatch.setattr(lp, name, dict(getattr(lp, name)))
    for name in ('UNC_FLOOR', 'LEVEL_HFS', 'HFS_MODEL'):
        monkeypatch.setattr(lp, name, getattr(lp, name))
    return lp


def quiet(*_a):
    pass


# ---------------------------------------------------------------------------
# The configuration
# ---------------------------------------------------------------------------
def test_a_configuration_without_widths_has_none(tmp_path):
    assert config.load(write_config(tmp_path, '')).widths_of == {}


def test_the_widths_table_is_read(tmp_path):
    w = config.load(write_config(tmp_path, (
        '[widths]\nunc_floor = 0.01\nlevel_hfs_widths = false\n'
        '[widths.reading]\n1974 = 0.002\n'
        '[widths.char_dlam]\n"1974/c" = 0.01\n'
        '[widths.kappa_scatter]\n1969 = 0.25\n'))).widths_of
    assert w['unc_floor'] == 0.01 and w['level_hfs_widths'] is False
    assert w['reading'] == {'1974': 0.002}
    assert w['char_dlam'] == {'1974/c': 0.01}
    assert w['kappa_scatter'] == {'1969': 0.25}


@pytest.mark.parametrize('text', [
    '[widths.reading]\n1975 = 0.002\n',               # no such era
    '[widths.char_dlam]\n"c" = 0.01\n',               # no era
    '[widths.char_dlam]\n"1974/" = 0.01\n',           # no code
    '[widths.kappa_scatter]\n"1974/c" = 0.1\n',       # an era table
    '[widths.char_dwn]\n"1974/w" = -0.01\n',          # negative
    '[widths]\nunc_floor = -0.01\n',
    '[widths]\nfloor = 0.01\n',                       # unknown key
    '[widths.reading]\n1974 = "x"\n',                 # not a number
])
def test_a_bad_widths_table_is_refused(tmp_path, text):
    with pytest.raises(config.ConfigError):
        config.load(write_config(tmp_path, text))


def test_final_carries_its_own_widths_and_the_others_none():
    final = config.load(os.path.join(ROOT, 'final', config.CONFIG_NAME))
    assert final.widths_of['level_hfs_widths'] is False
    assert set(final.widths_of['kappa_scatter']) == {'1974', '1969'}
    for d in ('', 'iter', 'iter_hfs'):
        cfg = config.load(os.path.join(ROOT, d, config.CONFIG_NAME))
        assert cfg.widths_of == {}


# ---------------------------------------------------------------------------
# apply_widths
# ---------------------------------------------------------------------------
def test_no_table_changes_nothing(tmp_path, fresh):
    before = (lp.UNC_FLOOR, dict(lp.DLAM_ERA), dict(lp.CHAR_DLAM),
              dict(lp.CHAR_DWN), lp.LEVEL_HFS)
    lp.apply_widths(config.load(write_config(tmp_path, '')), quiet)
    assert (lp.UNC_FLOOR, lp.DLAM_ERA, lp.CHAR_DLAM, lp.CHAR_DWN,
            lp.LEVEL_HFS) == before
    assert lp.KAPPA_SCATTER == {} and lp.HFS_MODEL is None


def test_the_table_replaces_what_it_names(tmp_path, fresh):
    keep = lp.CHAR_DLAM[('1974', 'bl')]
    lp.apply_widths(config.load(write_config(tmp_path, (
        '[widths]\nunc_floor = 0.0125\n'
        '[widths.reading]\n1969 = 0.0024\n'
        '[widths.char_dlam]\n"1974/c" = 0.0104\n'))), quiet)
    assert lp.UNC_FLOOR == 0.0125
    assert lp.DLAM_ERA == {'1974': 0.0030, '1969': 0.0024}
    assert lp.CHAR_DLAM[('1974', 'c')] == 0.0104
    assert lp.CHAR_DLAM[('1974', 'bl')] == keep         # not named: kept
    # the floor reaches a plain line where the reading error is small
    s = float(lp.meas_sigma([9000.0], [''])[0])
    assert s == pytest.approx(math.hypot(0.0030 * 9000.0 ** 2 * 1e-8, 0.0125))


def test_a_code_has_one_form(tmp_path, fresh):
    lp.apply_widths(config.load(write_config(tmp_path, (
        '[widths.char_dwn]\n"1974/c" = 0.02\n'
        '[widths.char_dlam]\n"1974/w" = 0.003\n'))), quiet)
    assert ('1974', 'c') not in lp.CHAR_DLAM
    assert lp.CHAR_DWN[('1974', 'c')] == 0.02
    assert ('1974', 'w') not in lp.CHAR_DWN
    assert lp.CHAR_DLAM[('1974', 'w')] == 0.003
    # a wavenumber-constant excess is the same everywhere
    a = lp.meas_sigma([20000.0, 40000.0], ['c', 'c'])
    b = lp.meas_sigma([20000.0, 40000.0], ['', ''])
    assert np.sqrt(a ** 2 - b ** 2) == pytest.approx([0.02, 0.02])


def test_the_level_widths_can_be_turned_off(tmp_path, fresh):
    lp.apply_widths(config.load(write_config(
        tmp_path, '[widths]\nlevel_hfs_widths = false\n')), quiet)
    assert lp.read_hfs_widths(os.path.join(ROOT, 'level_hfs_widths.csv')) \
        == {}


def test_the_kappa_scatter_needs_the_hfs_model(tmp_path, fresh):
    with pytest.raises(SystemExit, match='hfs'):
        lp.apply_widths(config.load(write_config(
            tmp_path, '[widths.kappa_scatter]\n1969 = 0.25\n')), quiet)


def test_the_kappa_scatter(tmp_path, fresh):
    lp.apply_widths(config.load(write_config(
        tmp_path, '[widths.kappa_scatter]\n1974 = 0.05\n1969 = 0.25\n',
        parent='iter_hfs/lineclass_config.toml')), quiet)
    assert lp.HFS_MODEL is not None
    v = lp.kappa_scatter_var([30000.0, 60000.0], [0.4, 0.4])
    assert v == pytest.approx([(0.05 * 0.4) ** 2, (0.25 * 0.4) ** 2])
    low, upp = '059003.000141', '059003.000452'
    d = lp.abs_D([low], [upp])[0]
    assert d == pytest.approx(abs(lp.HFS_MODEL.D(low, upp)[0]))
    assert (low, upp) in lp._D_CACHE


def test_without_the_model_every_D_is_zero(fresh):
    lp.HFS_MODEL = None
    assert list(lp.abs_D(['a', 'b'], ['c', 'd'])) == [0.0, 0.0]
    assert list(lp.kappa_scatter_var([30000.0], [0.4])) == [0.0]
