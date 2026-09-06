"""Tests for the TRANS.DAT half of `tools/estimate_snr.py`.

Run from the LineClass directory:  python -m pytest tests -q

A row of IDEN2's TRANS.DAT carries two intensities and two wavenumbers in
fixed-width fields.  The one the conversion needs is the wavenumber of the
predicted transition in columns 25-38, not the energy of the level at the
other end of it in columns 11-22; the two are different numbers, so reading
the wrong field evaluates the plate calibration and the noise level at a
wavelength belonging to some other part of the spectrum entirely.  These tests
pin the layout down on synthetic rows.
"""
import math
import os
import sys

import numpy as np
import pytest

TOOLS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     'tools')
sys.path.insert(0, TOOLS)

import estimate_snr as snr                      # noqa: E402


# a row of the real file, and the same row with the asterisk that marks a
# level whose energy is known experimentally (it shifts nothing: the fields
# are cut at fixed positions, but a naive split on spaces would move them)
FREE = ('+ 248  -34  139896.000      10099.300     0'
        '       0.000      0.000     0')
ASSIGNED = ('+ 689   73  118610.935 *    22829.907    59'
            '   22829.901     -0.006  5650')


class Opt(object):
    """The handful of attributes the two functions read off the options."""

    def __init__(self, **kw):
        self.iden2_k = 10.0
        self.trans_scale = 1000.0
        self.trans_zero_value = 0.0
        self.scale = 'log'
        self.match_tol = 0.02
        self.__dict__.update(kw)


class Ones(object):
    """A region whose polynomial is 0 in ln, i.e. a factor of 1 everywhere."""

    lo, hi = 1.0, 1.0e9

    def __call__(self, lam):
        return np.zeros(np.shape(lam))


@pytest.fixture
def flat():
    """A calibration, or a noise level, that is 1 at every wavelength, so that
    the only thing left in the written number is the wavelength used."""
    return [Ones()]


def write(tmp_path, rows):
    p = tmp_path / 'trans.dat'
    p.write_text('$   1    J= 1.5   f5d6d~3D4D   149995.300   \n'
                 + '\n'.join(rows) + '\n', encoding='latin-1')
    return str(p)


# ---------------------------------------------------------------------------
# the fields
# ---------------------------------------------------------------------------

def test_the_fields_are_where_the_constants_say():
    assert snr.trans_field(FREE, snr.TRANS_V) == -34.0
    assert snr.trans_field(FREE, snr.TRANS_EPART) == 139896.0
    assert snr.trans_field(FREE, snr.TRANS_WN) == 10099.3
    assert snr.trans_field(FREE, snr.TRANS_LINE) == 0.0

    assert snr.trans_field(ASSIGNED, snr.TRANS_V) == 73.0
    assert snr.trans_field(ASSIGNED, snr.TRANS_EPART) == 118610.935
    assert snr.trans_field(ASSIGNED, snr.TRANS_WN) == 22829.907
    assert snr.trans_field(ASSIGNED, snr.TRANS_VOBS) == 59.0
    assert snr.trans_field(ASSIGNED, snr.TRANS_WNOBS) == 22829.901
    assert snr.trans_field(ASSIGNED, snr.TRANS_LINE) == 5650.0


def test_the_asterisk_does_not_move_the_wavenumber():
    """Splitting on spaces would; cutting at fixed positions does not."""
    assert ASSIGNED[22:24] == ' *'
    assert snr.trans_field(ASSIGNED, snr.TRANS_WN) == 22829.907


# ---------------------------------------------------------------------------
# the conversion
# ---------------------------------------------------------------------------

def test_the_wavelength_used_is_the_transitions_own(tmp_path, flat):
    """With a flat calibration and a flat noise level the written value is
    v - 10*ln(1000) whatever the wavelength, so the test uses a noise level
    that is 1 at the transition's wavelength and e^5 at the wavelength the
    partner energy would give, and looks at which one came out."""
    lam_true = 1e8 / 10099.3          # 9902 A
    lam_wrong = 1e8 / 139896.0        # 715 A

    class Step(object):
        lo, hi = 1.0, 1.0e9

        def __call__(self, lam):
            lam = np.asarray(lam, float)
            return np.where(lam < 5000.0, 5.0, 0.0)

    opt = Opt()
    src = write(tmp_path, [FREE])
    dst = str(tmp_path / 'out.dat')
    snr.rewrite_trans(src, dst, [Step()], flat, opt,
                      snr.Matcher([1.0], [1.0]), lambda m: None)
    got = snr.trans_field(open(dst, encoding='latin-1').read().splitlines()[1],
                          snr.TRANS_V)

    expect = 10 * (math.log(math.exp(-34 / 10.0) / 1000.0) - 0.0)
    assert lam_true > 5000.0 > lam_wrong
    assert got == pytest.approx(round(expect), abs=0.5)


def test_the_observed_column_is_rewritten_from_the_line_list(tmp_path, flat):
    opt = Opt()
    src = write(tmp_path, [ASSIGNED])
    dst = str(tmp_path / 'out.dat')
    # the line list holds the observed line of the row, at SNR = e^2
    matcher = snr.Matcher([22829.901], [math.exp(2.0)])
    snr.rewrite_trans(src, dst, flat, flat, opt, matcher, lambda m: None)
    row = open(dst, encoding='latin-1').read().splitlines()[1]
    assert snr.trans_field(row, snr.TRANS_VOBS) == 20.0
    # and nothing else on the row moved
    assert len(row) == len(ASSIGNED)
    assert row[snr.TRANS_EPART[0]:snr.TRANS_EPART[1]] == \
        ASSIGNED[snr.TRANS_EPART[0]:snr.TRANS_EPART[1]]
    assert row[snr.TRANS_VOBS[1]:] == ASSIGNED[snr.TRANS_VOBS[1]:]


def test_an_unidentified_row_keeps_its_observed_column(tmp_path, flat):
    opt = Opt()
    src = write(tmp_path, [FREE])
    dst = str(tmp_path / 'out.dat')
    snr.rewrite_trans(src, dst, flat, flat, opt,
                      snr.Matcher([10099.3], [math.exp(2.0)]), lambda m: None)
    row = open(dst, encoding='latin-1').read().splitlines()[1]
    assert row[snr.TRANS_VOBS[0]:snr.TRANS_VOBS[1]] == \
        FREE[snr.TRANS_VOBS[0]:snr.TRANS_VOBS[1]]


def test_the_zero_point_uses_the_rows_own_identification(tmp_path, flat):
    """The offset must be measured against the line named in columns 45-57,
    not against whatever line happens to sit at the partner energy."""
    opt = Opt()
    src = write(tmp_path, [ASSIGNED])
    # predicted: exp(73/10)/1000 = 1.4880 ; observed SNR e^1 -> z = 1 - ln(1.488)
    matcher = snr.Matcher([22829.901], [math.exp(1.0)])
    z = snr.trans_zero(src, flat, flat, opt, matcher, lambda m: None)
    assert z == pytest.approx(1.0 - math.log(math.exp(7.3) / 1000.0), abs=1e-9)


def test_the_zero_point_is_zero_when_nothing_matches(tmp_path, flat):
    opt = Opt()
    src = write(tmp_path, [FREE])
    z = snr.trans_zero(src, flat, flat, opt, snr.Matcher([1.0], [1.0]),
                       lambda m: None)
    assert z == 0.0
