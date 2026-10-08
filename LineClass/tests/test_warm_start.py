"""Tests for the warm start of classify_lines.py (files.start_levels,
2026-10-08): the levels start from the set's last LOPT fit instead of their
published energies.

The tests fix:

  * the reader of LOPT_output_levels.txt, and that an empty file - what a
    failed LOPT run leaves - reads as no fit at all;
  * the rule for a level revised by hand: it starts from the fit once the
    fit has seen the revision, and from the revision while it has not;
  * apply_start_levels: levels the fit does not hold keep their energy, and
    a missing or empty file changes nothing;
  * that files.start_levels is read, resolved against the set's directory,
    and left empty where a configuration does not name it.

Run from the LineClass directory:  python -m pytest tests -q
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import classify_lines as cl
import config
from models import EnergyLevel

HEADER = 'Designation\tEnergy\tD1\tD2stat\tD2sys\tD2tot\tN_lines\tComments\n'


def level(lid, energy):
    return EnergyLevel(level_id=lid, energy=energy, parity='e', J_str='2.5',
                       J_val=2.5)


def lopt_file(tmp_path, rows):
    path = tmp_path / 'LOPT_output_levels.txt'
    path.write_text(HEADER + ''.join(
        '%s\t%.7f\t0.01\t0.01\t0\t0.01\t3\t\n' % r for r in rows),
        encoding='utf-8')
    return str(path)


@pytest.fixture(autouse=True)
def no_revisions(monkeypatch):
    monkeypatch.setattr(cl, 'LEVEL_REVISIONS', {})


def test_read_start_energies(tmp_path):
    path = lopt_file(tmp_path, [('059003.000001', 0.0),
                                ('059003.000544', 105190.8373)])
    assert cl.read_start_energies(path) == {
        '059003.000001': 0.0, '059003.000544': pytest.approx(105190.8373)}


def test_an_empty_file_is_no_fit(tmp_path):
    path = tmp_path / 'LOPT_output_levels.txt'
    path.write_text('', encoding='utf-8')
    assert cl.read_start_energies(str(path)) == {}


def test_a_file_without_the_columns_stops_the_run(tmp_path):
    path = tmp_path / 'LOPT_output_levels.txt'
    path.write_text('level_id\tE\n059003.000001\t0\n', encoding='utf-8')
    with pytest.raises(ValueError):
        cl.read_start_energies(str(path))


def test_start_from_fit_rule():
    assert cl.start_from_fit(105190.84)                     # not revised
    # revised from 141328.85 to 141440.84, and fitted there: the fit
    assert cl.start_from_fit(141440.78, (141328.85, 141440.84, ''))
    # revised, but the fit is still at the old position: the revision
    assert not cl.start_from_fit(141328.90, (141328.85, 141440.84, ''))
    # exactly halfway counts as not seen
    assert not cl.start_from_fit(141384.845, (141328.85, 141440.84, ''))


def test_apply_start_levels(tmp_path, monkeypatch):
    levels = {lid: level(lid, e) for lid, e in [
        ('059003.000544', 105191.13),       # published, fitted lower
        ('059003.000622', 141440.8418),     # revised, the fit has seen it
        ('059003.000700', 120000.00),       # revised, the fit has not
        ('059003.000701', 90000.00)]}       # not in the fit
    monkeypatch.setattr(cl, 'LEVEL_REVISIONS', {
        '059003.000622': (141328.85, 141440.8418, ''),
        '059003.000700': (119000.00, 120000.00, '')})
    path = lopt_file(tmp_path, [('059003.000544', 105190.8373),
                                ('059003.000622', 141440.7840),
                                ('059003.000700', 119000.0300)])
    assert cl.apply_start_levels(levels, path) == 2
    assert levels['059003.000544'].energy == pytest.approx(105190.8373)
    assert levels['059003.000622'].energy == pytest.approx(141440.7840)
    assert levels['059003.000700'].energy == 120000.00
    assert levels['059003.000701'].energy == 90000.00


def test_a_failed_fit_changes_nothing(tmp_path):
    levels = {'059003.000544': level('059003.000544', 105191.13)}
    empty = tmp_path / 'LOPT_output_levels.txt'
    empty.write_text('', encoding='utf-8')
    assert cl.apply_start_levels(levels, str(empty)) == 0
    assert cl.apply_start_levels(levels, str(tmp_path / 'none.txt')) == 0
    assert levels['059003.000544'].energy == 105191.13


def test_the_configuration_names_the_file(tmp_path):
    base = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), 'lineclass_config.toml')
    assert config.load(base).start_levels == ''
    (tmp_path / 'lineclass_config.toml').write_text(
        'inherit = %r\n\n[files]\noutput = "t.xlsx"\noutput_csv = "t.csv"\n'
        'start_levels = "LOPT_output_levels.txt"\n'
        % base.replace('\\', '/'), encoding='utf-8')
    cfg = config.load(str(tmp_path / 'lineclass_config.toml'))
    assert cfg.start_levels == os.path.normpath(
        str(tmp_path / 'LOPT_output_levels.txt'))
