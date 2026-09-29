"""The hand-set uncertainties of inflated_unc_lines.txt reach classify_lines.

Run from the LineClass directory:  python -m pytest tests -q

A line whose uncertainty was widened by hand - in IDEN2, or after a trial
LOPT run - must be classified on that uncertainty, or the classification
throws out again what the analyst accepted and LOPT is given a value the
classification never used.
"""
import io
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (ROOT, os.path.dirname(ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

import config                                  # noqa: E402
import classify_lines as C                     # noqa: E402
from models import SpectralLine                # noqa: E402

from test_working_sets import BASE, CHILD, write   # noqa: E402


REGISTRY = ('wn_key\tunc_wn\tdate\treason\n'
            '44652.5408\t0.1948\t9/23/2026\toutlier of the block fit\n'
            '37577.0988\t0.2500\t9/27/2026\tunc inflated (hfs of f26s)\n')


def line(wn, unc, key=None):
    return SpectralLine(wavenumber=wn, wn_uncertainty=unc, intensity=1.0,
                        line_character='', wn_key=wn if key is None else key)


def test_a_listed_line_takes_the_registry_value(tmp_path):
    path = write(str(tmp_path / 'inflated_unc_lines.txt'), REGISTRY)
    a = line(44652.5408, 0.05)
    b = line(37577.0988, 0.11)
    c = line(30000.0, 0.03)
    assert C.apply_inflated_uncertainties([a, b, c], path) == 2
    assert (a.wn_uncertainty, b.wn_uncertainty) == (0.1948, 0.25)
    assert c.wn_uncertainty == 0.03


def test_the_line_is_found_by_its_name_not_its_corrected_wavenumber(tmp_path):
    """In iter/ the wavenumber has moved by the calibration; wn_key has not."""
    path = write(str(tmp_path / 'inflated_unc_lines.txt'), REGISTRY)
    a = line(44652.5612, 0.05, key=44652.5408)
    b = line(37577.3000, 0.11, key=37577.0988)
    C.apply_inflated_uncertainties([a, b], path)
    assert (a.wn_uncertainty, b.wn_uncertainty) == (0.1948, 0.25)


def test_the_registry_never_narrows_a_line(tmp_path):
    """In the baseline, Sugar's stated value can exceed the registry's."""
    path = write(str(tmp_path / 'inflated_unc_lines.txt'), REGISTRY)
    a = line(44652.5408, 0.30)
    b = line(37577.0988, 0.11)
    assert C.apply_inflated_uncertainties([a, b], path) == 1
    assert (a.wn_uncertainty, b.wn_uncertainty) == (0.30, 0.25)


def test_an_entry_naming_no_line_stops_the_run(tmp_path):
    path = write(str(tmp_path / 'inflated_unc_lines.txt'), REGISTRY)
    with pytest.raises(ValueError, match='37577.0988'):
        C.apply_inflated_uncertainties([line(44652.5408, 0.05)], path)


def test_a_missing_registry_changes_nothing(tmp_path):
    a = line(44652.5408, 0.05)
    assert C.apply_inflated_uncertainties(
        [a], str(tmp_path / 'inflated_unc_lines.txt')) == 0
    assert a.wn_uncertainty == 0.05


def test_the_registry_is_shared_by_the_sets(tmp_path):
    """Named in the baseline, it is the baseline's file for iter/ too."""
    root = str(tmp_path / 'LineClass')
    write(os.path.join(root, 'lineclass_config.toml'),
          BASE.replace('line_decisions = "line_decisions.csv"',
                       'line_decisions = "line_decisions.csv"\n'
                       'inflated_unc = "inflated_unc_lines.txt"'))
    write(os.path.join(root, 'iter', 'lineclass_config.toml'), CHILD)
    base = config.load(os.path.join(root, 'lineclass_config.toml'))
    iter_ = config.load(os.path.join(root, 'iter', 'lineclass_config.toml'))
    assert iter_.inflated_unc == base.inflated_unc
    assert iter_.inflated_unc == os.path.join(root, 'inflated_unc_lines.txt')


def test_the_registry_is_optional_in_the_configuration(tmp_path):
    root = str(tmp_path / 'LineClass')
    write(os.path.join(root, 'lineclass_config.toml'), BASE)
    assert config.load(os.path.join(root, 'lineclass_config.toml')) \
        .inflated_unc == ''


def test_a_hand_typed_key_with_fewer_decimals_finds_its_line(tmp_path):
    """53657.57 in the registry names the line 53657.5706, as in the ledger."""
    path = write(str(tmp_path / 'inflated_unc_lines.txt'),
                 'wn_key\tunc_wn\tdate\treason\n53657.57\t0.5\t\tretained\n')
    a = line(53657.5706, 0.1153)
    b = line(53657.6900, 0.1)
    C.apply_inflated_uncertainties([a, b], path)
    assert (a.wn_uncertainty, b.wn_uncertainty) == (0.5, 0.1)


def test_a_corrected_wavenumber_is_not_a_name(tmp_path):
    """100444.8406 is iter's own_corr of the line named 100444.8853."""
    path = write(str(tmp_path / 'inflated_unc_lines.txt'),
                 'wn_key\tunc_wn\tdate\treason\n100444.8406\t0.5\t\tblend\n')
    with pytest.raises(ValueError, match='100444.8406'):
        C.apply_inflated_uncertainties(
            [line(100444.8406, 0.2886, key=100444.8853)], path)


def test_two_entries_for_one_line_must_agree(tmp_path):
    path = write(str(tmp_path / 'inflated_unc_lines.txt'),
                 'wn_key\tunc_wn\tdate\treason\n'
                 '53657.57\t0.5\t\ta\n53657.5706\t0.3\t\tb\n')
    with pytest.raises(ValueError, match='different'):
        C.apply_inflated_uncertainties([line(53657.5706, 0.1153)], path)


# The registry itself, as every reader sees it (hfs_kappa.read_inflated).

MIXED = ('wn_key\tunc_wn\tdate\treason\n'
         '53657.57\t0.5\t\tcopied from LOPT\n'
         '97588.5\t0.4\t\tExcel dropped the trailing zero\n'
         '41473.79906049907\t0.1048\t\tfull precision\n'
         '44652.5408\t0.1948\t\tfour decimals\n')


def test_an_entry_is_matched_at_the_precision_it_was_written_with(tmp_path):
    import hfs_kappa
    reg = hfs_kappa.read_inflated(
        write(str(tmp_path / 'inflated_unc_lines.txt'), MIXED))
    assert sorted(reg) == ['41473.7991', '44652.5408', '53657.57', '97588.5']
    assert reg.lookup(53657.57057117759) == 0.5
    assert reg.lookup(97588.49517756669) == 0.4
    assert reg.lookup(41473.79906049907) == 0.1048
    assert reg.lookup(44652.54081234) == 0.1948
    # within 0.01 of an entry, but not that entry at its own precision
    assert reg.lookup(53657.5649) is None
    assert reg.lookup(44652.5413) is None


def test_an_entry_too_short_to_tell_two_lines_apart_raises(tmp_path):
    import hfs_kappa
    reg = hfs_kappa.read_inflated(
        write(str(tmp_path / 'inflated_unc_lines.txt'), MIXED))
    with pytest.raises(ValueError,
                       match=r'53657\.57 \(53657\.5706, 53657\.5732\)'):
        reg.check([53657.5706, 53657.5732])
    path = str(tmp_path / 'inflated_unc_lines.txt')
    with pytest.raises(ValueError, match='more than one line'):
        C.apply_inflated_uncertainties(
            [line(53657.5706, 0.1), line(53657.5732, 0.1),
             line(97588.4952, 0.1), line(41473.7991, 0.1),
             line(44652.5408, 0.1)], path)


def test_check_returns_the_entries_that_name_no_line(tmp_path):
    import hfs_kappa
    reg = hfs_kappa.read_inflated(
        write(str(tmp_path / 'inflated_unc_lines.txt'), MIXED))
    assert reg.check([53657.5706, 44652.5408]) == ['41473.7991', '97588.5']


def test_one_key_written_twice_with_two_values_raises(tmp_path):
    import hfs_kappa
    path = write(str(tmp_path / 'inflated_unc_lines.txt'),
                 'wn_key\tunc_wn\tdate\treason\n'
                 '53657.5706\t0.5\t\ta\n53657.57057\t0.3\t\tb\n')
    with pytest.raises(ValueError, match='twice'):
        hfs_kappa.read_inflated(path)


def test_review_mismatches_does_not_add_a_line_listed_under_a_short_key(
        tmp_path):
    import review_mismatches as rm
    path = write(str(tmp_path / 'inflated_unc_lines.txt'), MIXED)
    said = []
    added = rm.append_inflated(
        path, [('53657.57057117759', '0.6', 'again'),
               ('12046.21296057019', '0.2', 'new')], said.append)
    assert added == [('12046.2130', '0.2')]
    assert any('53657.5706 is already in the registry as 53657.57' in s
               for s in said)
