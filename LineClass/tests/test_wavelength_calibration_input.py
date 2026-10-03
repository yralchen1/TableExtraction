"""Which wavenumber the calibration fit reads from a classification table.

The baseline's `line_classifications.csv` carries Sugar's published
wavenumber in `wn_obs` and has no `wn_key` column.  A working set's table
(`iter/`, `iter_hfs/`) carries the calibrated value in `wn_obs` and Sugar's in
`wn_key`.  The calibration is a correction to what Sugar published, so it must
read `wn_key` from a set's table; reading `wn_obs` there would fit the
residual of the last calibration and apply it on top of itself.
"""
import csv
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import hfs_kappa                               # noqa: E402
import wavelength_calibration as W            # noqa: E402

BASE = ['wn_obs', 'char', 'low_id', 'upp_id', 'accepted', 'n_accepted']


def _write(path, header, rows):
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)
    return str(path)


@pytest.fixture
def tables(tmp_path):
    """The same three lines as the baseline and as a calibrated set write
    them: one singly accepted, one blend, one unclassified."""
    base = _write(tmp_path / 'base.csv', BASE, [
        ['23083.6428', '', 'a', 'b', '1.0', '1'],
        ['30000.1000', 'c', 'a', 'c', '1.0', '2'],
        ['30000.1000', 'c', 'b', 'd', '1.0', '2'],
        ['40000.5000', 'w', '', '', '', '0'],
    ])
    hdr = ['wn_obs', 'wn_key'] + BASE[1:]
    iset = _write(tmp_path / 'set.csv', hdr, [
        ['23083.6750', '23083.6428', '', 'a', 'b', '1.0', '1'],
        ['30000.0800', '30000.1000', 'c', 'a', 'c', '1.0', '2'],
        ['30000.0800', '30000.1000', 'c', 'b', 'd', '1.0', '2'],
        ['40000.5400', '40000.5000', 'w', '', '', '', '0'],
    ])
    return base, iset


def _read(path, raw):
    return hfs_kappa.read_lines(path, constants={}, J_of={}, raw=raw)


def test_a_sets_table_is_read_on_wn_key(tables):
    base, iset = tables
    got = _read(iset, raw=True)
    assert [ln.wn for ln in got] == [23083.6428]
    assert [ln.key for ln in got] == [23083.6428]
    # the baseline gives the same line, on the same value
    assert [ln.wn for ln in _read(base, raw=True)] == [23083.6428]
    assert [ln.wn for ln in _read(base, raw=False)] == [23083.6428]
    # without raw a set's table gives its calibrated value, as before
    assert [ln.wn for ln in _read(iset, raw=False)] == [23083.6750]


def test_every_observed_line_is_listed_on_wn_key(tables):
    base, iset = tables
    want = [(23083.6428, '', 23083.6428), (30000.1, 'c', 30000.1),
            (40000.5, 'w', 40000.5)]
    assert W.observed_wavenumbers(iset) == want
    assert W.observed_wavenumbers(base) == want


def test_the_table_is_found_or_the_run_stops(tables, tmp_path):
    _, iset = tables
    assert W.classifications_path(iset) == os.path.abspath(iset)
    assert W.classifications_path(None) == os.path.join(W.HERE,
                                                        hfs_kappa.LINES)
    with pytest.raises(SystemExit):
        W.classifications_path(str(tmp_path / 'missing.csv'))
