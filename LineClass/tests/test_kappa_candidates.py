"""Tests for kappa_candidates.py, the search for lines measured otherwise
than their kappa class says.

Run from the LineClass directory:  python -m pytest tests -q

The arithmetic is fixed on numbers worked by hand: the leave-one-out
residual and the uncertainty of the rest of the network from LOPT's two
uncertainties, the residual under another kappa, the verdicts, the model
uncertainty read from the calibration report, and the rung-1 evidence.
"""
import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import hfs_kappa                           # noqa: E402
import hfs_patterns                        # noqa: E402
import kappa_candidates as K               # noqa: E402

REPORT = """\
some text before
char          era       n       a_A     b_cm1    adopted      after  chi2/dof
(plain)      1974    1647    0.0022    0.0192     0.0383     0.0565      0.96
w            1974     427    0.0023    0.0199     0.0463     0.0821      0.90
other        1974      18       nan       nan        nan     0.0418      0.69
(plain)      1969    1794    0.0026    0.0557     0.1746     0.1894      0.97
other        1969      22    0.0000    0.3211     0.3211     0.3212      1.00

The lines that carry it
"""


def test_the_class_table_is_read(tmp_path):
    path = tmp_path / 'report.txt'
    path.write_text(REPORT, encoding='utf-8')
    t = K.read_unc_table(str(path))
    assert t[('', 1974)] == (0.0022, 0.0192)
    assert t[('w', 1974)] == (0.0023, 0.0199)
    assert ('other', 1974) not in t              # nan: left out
    assert t[('other', 1969)] == (0.0, 0.3211)
    # the model value: the class's own, the pooled one, Sugar's stated one
    wn = 30000.0
    assert K.model_unc(t, 'w', wn) == pytest.approx(
        hfs_kappa.two_term(0.0023, 0.0199, wn))
    assert K.model_unc(t, 'bl', 60000.0) == pytest.approx(0.3211)
    assert K.model_unc(t, 'zz', wn) == pytest.approx(max(
        hfs_kappa.stated_uncertainty('zz', 1974, wn), hfs_kappa.FLOOR))


def test_a_report_without_the_table_stops(tmp_path):
    path = tmp_path / 'report.txt'
    path.write_text('nothing here\n', encoding='utf-8')
    with pytest.raises(SystemExit, match='no table of uncertainty'):
        K.read_unc_table(str(path))


def test_leave_one_out():
    # u_calc = 0.03 against u_obs = 0.05: h = 0.36
    h, r_loo, u_rest = K.leave_one_out(0.064, 0.05, 0.03)
    assert h == pytest.approx(0.36)
    assert r_loo == pytest.approx(0.1)
    assert u_rest == pytest.approx(0.0375)
    # the rest and the line combine back into u_calc
    assert 1 / u_rest ** 2 + 1 / 0.05 ** 2 == pytest.approx(1 / 0.03 ** 2)
    # a line that alone fixes its Ritz value cannot be tested
    assert K.leave_one_out(0.0, 0.05, 0.05)[1:] == (None, None)


def test_the_residual_under_another_kappa():
    # given at kappa 0.8 with D = +0.5: at the cg it would sit 0.4 lower
    z = K.hypotheses(0.0, 0.8, 0.8, 0.5, 0.1)
    assert z['now'] == z['class'] == pytest.approx(0.0)
    assert z['head'] == pytest.approx(-1.0)
    assert z['cg'] == pytest.approx(4.0)
    # a held line (kappa 0 now) seen at its class
    z = K.hypotheses(0.0, 0.0, 0.8, 0.5, 0.1)
    assert z['now'] == pytest.approx(0.0)
    assert z['class'] == pytest.approx(-4.0)
    assert z['cg'] == pytest.approx(0.0)
    assert 'held' not in z
    # a registry row newer than the fit: fitted at the class (0.8), held at
    # the cg, and the cg is where the line fits
    z = K.hypotheses(0.4, 0.8, 0.8, 0.5, 0.1, kappa_held=0.0)
    assert z['now'] == pytest.approx(4.0)
    assert z['held'] == z['cg'] == pytest.approx(8.0)
    z = K.hypotheses(-0.4, 0.8, 0.8, 0.5, 0.1, kappa_held=0.0)
    assert z['held'] == pytest.approx(0.0)
    assert K.verdict(z, 'cg', kappa_class=0.8)[0] == 'audit-ok'


def test_the_registry_advice():
    # no inflation: nothing to say
    assert K.registry_advice(None, {'now': 5.0}, 0.05, 0.01, 0.02) == ''
    # 40253.7799 held at the cg: 1.2 sigma, the inflation goes
    z = {'now': 11.8, 'held': 1.2}
    assert K.registry_advice(0.6, z, 0.0471, 0.0138, 0.0159) == 'delete'
    # 44494.0809: 2.6 sigma of 0.0946 is 0.246; less the hfs shift (0.0771)
    # and the rest (0.0233), the line's own value must be 0.232
    z = {'now': 2.6, 'held': 2.6}
    assert K.registry_advice(0.25, z, 0.0946, 0.0771, 0.0233) == \
        'needs 0.23 (now 0.25)'
    # a line not held is judged where it is fitted
    assert K.registry_advice(0.3, {'now': -1.5}, 0.05, 0.0, 0.0) == 'delete'
    # a candidate is judged where it would be held, and says so
    z = {'now': 10.8, 'cg': 0.9, 'head': 12.0}
    assert K.registry_advice(0.5, z, 0.05, 0.0, 0.0, at='cg') == \
        'at cg: delete'
    assert K.registry_advice(0.5, z, 0.05, 0.0, 0.0, at='head') == \
        'at head: needs 0.60 (now 0.5)'


def test_the_advice_by_verdict():
    # 40790.2000 fails its class (0.817): judged there, then at the cg
    # (0.308 cm-1 off: an own value of 0.31 with the shift and the rest)
    sigma, u_shift, u_rest = 0.0481, 0.0152, 0.0195
    z = {'now': 15.9, 'class': 15.9, 'head': 18.0, 'cg': 0.308 / sigma}
    assert K.advice_of('class-fails', 0.5, z, sigma, u_shift, u_rest,
                       kappa_class=0.817) == \
        'needs 0.76 (now 0.5); at cg: needs 0.31 (now 0.5)'
    # a flagged class has no head alternative, however close the head is
    z = {'now': 5.0, 'class': 5.0, 'head': 5.0, 'cg': 1.0}
    assert K.advice_of('class-fails', 0.3, z, 0.05, 0.0, 0.0,
                       kappa_class=1.0) == \
        'needs 0.25 (now 0.3); at cg: delete'
    # a candidate is judged only where it would be held; no inflation, none
    z = {'now': 6.0, 'class': 6.0, 'head': 7.0, 'cg': 0.5}
    assert K.advice_of('candidate-cg', 0.5, z, 0.05, 0.0, 0.0) == \
        'at cg: delete'
    assert K.advice_of('class-fails', None, z, 0.05, 0.0, 0.0) == ''
    assert K.advice_of('fits-class', 0.5, {'now': 1.0}, 0.05, 0.0, 0.0) == \
        'delete'


def test_an_inflation_that_covers_the_residual_passes_the_audit():
    # 37509.3400 held at the cg, -5.4 sigma of 0.0546: the residual is
    # 0.2923, and less the shift (0.0302) and the rest (0.025) it needs 0.29,
    # which is what the registry has
    sigma, u_shift, u_rest = 0.0546, 0.0302, 0.025
    z = {'now': -0.2923 / sigma, 'held': -0.2923 / sigma, 'class': -11.6,
         'head': -13.1, 'cg': -0.2923 / sigma}
    assert K.registry_advice(0.29, z, sigma, u_shift, u_rest) == \
        'needs 0.29 (now 0.29)'
    assert K.inflation_covers(0.29, z, sigma, u_shift, u_rest)
    assert K.verdict(z, 'cg', kappa_class=0.803, covered=True)[0] == \
        'audit-ok'
    # an inflation short of the need, or none, does not
    assert not K.inflation_covers(0.25, z, sigma, u_shift, u_rest)
    assert not K.inflation_covers(None, z, sigma, u_shift, u_rest)
    assert K.verdict(z, 'cg', kappa_class=0.803)[0] == 'audit-fails'
    # a larger inflation covers too; the advice then asks for less
    assert K.inflation_covers(0.4, z, sigma, u_shift, u_rest)


def test_the_median():
    assert K.median([3.0, 1.0, 2.0]) == 2.0
    assert K.median([4.0, 1.0, 2.0, 3.0]) == 2.5
    assert K.median([0.7]) == 0.7


@pytest.mark.parametrize('z, held, kc, want', [
    ({'class': 5.0, 'head': 6.0, 'cg': 0.5}, '', 0.8, 'candidate-cg'),
    ({'class': -4.0, 'head': 0.3, 'cg': -9.0}, '', 0.8, 'candidate-head'),
    ({'class': 5.0, 'head': 6.0, 'cg': 3.0}, '', 0.8, 'class-fails'),
    ({'class': 1.0, 'head': 2.0, 'cg': 4.0}, '', 0.8, 'fits-class'),
    # 2.5 sigma: not rejected, and the alternative does not gain 9
    ({'class': 2.5, 'head': 3.0, 'cg': 1.0}, '', 0.8, 'fits-class'),
    # a held line is audited at the kappa it is held at, failing as a class
    # does, beyond 3 sigma: 2.3 and 2.6 sigma pass (23792.4971, 44494.0809)
    ({'now': 0.5, 'class': 5.0, 'head': 6.0, 'cg': 0.5}, 'cg', 0.8,
     'audit-ok'),
    ({'now': 2.6, 'class': 5.0, 'head': 6.0, 'cg': 2.6}, 'cg', 0.8,
     'audit-ok'),
    ({'now': 2.3, 'class': -2.3, 'head': 2.3, 'cg': -22.6}, 'head', 0.817,
     'audit-ok'),
    ({'now': 3.5, 'class': 5.0, 'head': 6.0, 'cg': 3.5}, 'cg', 0.8,
     'audit-fails'),
    # a borrowed class is audited at its own kappa, not the line's class
    ({'now': 0.4, 'class': 4.0, 'head': 5.0, 'cg': -3.0}, 'c', 0.8,
     'audit-ok'),
    # a row newer than the fit is audited at its held kappa, not at 'now'
    ({'now': 11.8, 'held': 1.2, 'class': 11.8, 'head': 14.1, 'cg': 1.2},
     'cg', 0.817, 'audit-ok'),
    ({'now': 0.5, 'held': 4.0, 'class': 0.5, 'head': 2.0, 'cg': 4.0},
     'cg', 0.817, 'audit-fails'),
    # a flagged line has no head alternative: its class is the head
    ({'class': -4.0, 'head': -4.0, 'cg': -9.0}, '', 1.0, 'class-fails'),
])
def test_verdicts(z, held, kc, want):
    assert K.verdict(z, held, kappa_class=kc)[0] == want


def test_the_rung_and_the_line_there():
    class M:
        resolved = {'R'}
        levels = {'L': (3.5, -0.03, 0.0, True), 'U': (2.5, 0.05, 0.0, True),
                  'R': (0.5, -0.157, 0.0, True)}
    assert K.rung1(M, 'L', 'U') == pytest.approx(
        hfs_patterns.rung(1, 3.5, -0.03, 2.5, 0.05))
    assert K.rung1(M, 'R', 'U') is None
    keys = [100.0, 100.4, 101.0]
    assert K.line_near(keys, 100.45) == 100.4
    assert K.line_near(keys, 100.7) is None
    assert K.nearest_other(keys, [1, 2, 3], 1) == (pytest.approx(-0.4), 1)
