"""Tests for the intensity term of the line uncertainty model (2026-11-19):
a weak line is set on less accurately, so its class variance is multiplied by
hfs_kappa.intensity_factor(SNR, k), k per era, fitted in the plate
calibration with the classes and the spread of kappa (hfs_kappa.fit_classes).

Run from the LineClass directory:  python -m pytest tests -q

The tests fix:

  * the factor: 1 at SNR_REF, larger for weaker lines, 1 without an SNR or k;
  * the joint fit recovers a known k from synthetic residuals, and k = 0
    where the residuals do not depend on the SNR or no line has one;
  * a line's sigma and the character floor carry the factor, the spread of
    kappa does not;
  * the SNR file reader.
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

import hfs_kappa                           # noqa: E402


def synthetic(n, a, b, k, seed, s=0.0, era_wn=(20000.0, 47000.0)):
    """Residuals of n lines under the model with intensity coefficient k:
    SNR log-uniform over 0.5-1000, |D| as in test_kappa_spread."""
    rng = np.random.default_rng(seed)
    wn = rng.uniform(*era_wn, n)
    d = np.abs(rng.normal(0.0, 0.2, n))
    h = rng.uniform(0.0, 0.4, n)
    snr = 10.0 ** rng.uniform(math.log10(0.5), 3.0, n)
    f = np.array([hfs_kappa.intensity_factor(x, k) for x in snr])
    u = np.sqrt(f * ((a * 1e-8 * wn * wn) ** 2 + b * b) + (s * d) ** 2)
    r = rng.normal(0.0, 1.0, n) * u * np.sqrt(1.0 - h)
    return [(float(r[i]), float(h[i]), float(wn[i]), float(d[i]),
             float(snr[i])) for i in range(n)]


def test_the_factor_is_one_at_the_reference_and_grows_toward_weak_lines():
    k = 2.0
    assert hfs_kappa.intensity_factor(hfs_kappa.SNR_REF, k) == pytest.approx(1)
    assert hfs_kappa.intensity_factor(1.0, k) > 1.0
    assert hfs_kappa.intensity_factor(1000.0, k) < 1.0
    assert (hfs_kappa.intensity_factor(0.5, k)
            > hfs_kappa.intensity_factor(1.0, k))
    assert hfs_kappa.intensity_factor(1.0, 0.0) == 1.0
    assert hfs_kappa.intensity_factor(None, k) == 1.0
    assert hfs_kappa.intensity_factor(0.0, k) == 1.0


def test_the_joint_fit_recovers_k():
    groups = {('', 1974): synthetic(5000, 0.0022, 0.016, 2.0, 1)}
    fitted, _, k = hfs_kappa.fit_classes(groups, spread=False)
    a, b = fitted[('', 1974)]
    assert k[1974] == pytest.approx(2.0, rel=0.25)
    assert a == pytest.approx(0.0022, rel=0.15)
    assert b == pytest.approx(0.016, rel=0.2)


def test_k_and_the_spread_are_fitted_together():
    groups = {('', 1974): synthetic(6000, 0.0022, 0.016, 2.0, 2, s=0.07)}
    _, s, k = hfs_kappa.fit_classes(groups)
    assert k[1974] == pytest.approx(2.0, rel=0.3)
    assert s[1974] == pytest.approx(0.07, rel=0.2)


def test_no_k_where_the_residuals_do_not_depend_on_the_snr():
    groups = {('', 1974): synthetic(4000, 0.0022, 0.016, 0.0, 3)}
    _, _, k = hfs_kappa.fit_classes(groups, spread=False)
    assert k[1974] < 0.3


def test_no_k_without_an_snr():
    groups = {('', 1974): [v[:4] for v in synthetic(500, 0.0022, 0.016,
                                                    2.0, 4)]}
    _, _, k = hfs_kappa.fit_classes(groups, spread=False)
    assert k[1974] == 0.0
    _, _, k = hfs_kappa.fit_classes(groups, spread=False, intensity=False)
    assert k == {}


def _line(char='', snr=None, D=0.0, cls='plain_1974'):
    return hfs_kappa.Line(wn=30000.0, char=char, era=1974, cls=cls,
                          ucls=(char, 1974), low='a', upp='b', D=D, dJ=0,
                          snr=snr)


def test_a_weak_line_gets_a_larger_sigma_and_the_spread_is_not_scaled():
    ab = (0.0022, 0.0143)
    k = 2.0
    base = hfs_kappa.line_sigma(_line(snr=10.0), ab, ab, 0.0, k)
    weak = hfs_kappa.line_sigma(_line(snr=1.0), ab, ab, 0.0, k)
    assert base == pytest.approx(hfs_kappa.two_term(*ab, 30000.0))
    assert weak == pytest.approx(
        base * math.sqrt(hfs_kappa.intensity_factor(1.0, k)))
    # the spread of kappa adds in quadrature, unscaled
    spread = hfs_kappa.line_sigma(_line(snr=1.0, D=0.3), ab, ab, 0.1, k)
    assert spread == pytest.approx(math.hypot(weak, 0.1 * 0.3))
    # no SNR: the model before the term
    assert (hfs_kappa.line_sigma(_line(), ab, ab, 0.0, k)
            == pytest.approx(base))


def test_the_character_floor_compares_with_a_plain_line_as_strong():
    plain = (0.0022, 0.0143)
    f = hfs_kappa.intensity_factor(1.0, 2.0)
    floor = hfs_kappa.character_floor('w', 1974, 30000.0, plain, f)
    assert floor == pytest.approx(hfs_kappa.CHARACTER_FACTOR * max(
        hfs_kappa.two_term(*plain, 30000.0) * math.sqrt(f), hfs_kappa.FLOOR))
    assert hfs_kappa.character_floor('', 1974, 30000.0, plain, f) == 0.0


def test_read_snr(tmp_path):
    p = tmp_path / 'snr.csv'
    p.write_text('wn,lambda_A,Iorig,noise,snr,iden2\n'
                 '121663.0955,821.9419,10,1.234,8.106,21\n'
                 '121585.40794,822.4671,15,1.228,12.21,25\n', encoding='utf-8')
    got = hfs_kappa.read_snr(str(p))
    assert got == {121663.0955: 8.106, 121585.4079: 12.21}
    assert hfs_kappa.read_snr(str(tmp_path / 'none.csv')) == {}


def test_the_registry_only_widens_in_the_calibration(monkeypatch):
    """A line of inflated_unc_lines.txt gets the larger of its registry
    value and the model's (2026-11-19); before, the registry replaced the
    model even where it was the smaller."""
    rng = np.random.default_rng(7)
    lines, residual, leverage = [], [], []
    for i in range(80):
        wn = 20000.0 + 300.0 * i
        lines.append(hfs_kappa.Line(wn=wn, char='', era=1974,
                                    cls='plain_1974', ucls=('', 1974),
                                    low='a', upp='b%d' % i, D=0.0, dJ=0,
                                    key=wn, snr=1.0 + i))
        residual.append(float(rng.normal(0.0, 0.03)))
        leverage.append(0.1)
    narrow, wide = lines[3].key, lines[5].key
    reg = hfs_kappa.Registry({'%.4f' % narrow: 0.0001, '%.4f' % wide: 5.0})
    monkeypatch.setattr(hfs_kappa, 'read_inflated', lambda *a, **k: reg)
    sigma, table = hfs_kappa.refit_uncertainties(lines, residual, leverage,
                                                 spread=False)
    ab = table[('', 1974)][:2]
    k = table.intensity.get(1974, 0.0)
    model = hfs_kappa.line_sigma(lines[3], ab, ab, 0.0, k)
    assert sigma[3] == pytest.approx(model)        # registry smaller: model
    assert sigma[5] == 5.0                         # registry larger: kept
