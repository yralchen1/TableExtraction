"""Tests for `blend_centroid` and the two level fits that use it.

Run from the LineClass directory:  python -m pytest tests -q

A blend is fitted as the centroid of its components, weighted by their
calculated intensities, and those intensities are uncertain.  The tests fix:

  * the formula: two equal components delta apart give u_ln * delta / sqrt(8),
    a lopsided pair less, one component nothing;
  * that classify_lines.calc_weights widens a blend's weight by it, keeps the
    value it used as the line's fit_uncertainty, and still writes the
    branching fraction (not something divided by u_I) in the table;
  * that make_LOPT_input gives a blend's records one uncertainty widened by
    it, and leaves a single line and the weights alone.
"""
import io
import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (ROOT, os.path.dirname(ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

import blend_centroid as B                 # noqa: E402
import classify_lines as cl                # noqa: E402
import config                              # noqa: E402
import make_LOPT_input as M                # noqa: E402
from models import EnergyLevel, SpectralLine, Transition   # noqa: E402

LOW, UP1, UP2 = '059003.000100', '059003.000300', '059003.000301'


# ---------------------------------------------------------------------------
# the formula
# ---------------------------------------------------------------------------
def test_two_equal_components():
    assert B.intensity_unc([1.0, 1.0], [0.0, 0.8], 0.33) == pytest.approx(
        0.33 * 0.8 / math.sqrt(8))


def test_a_lopsided_pair_gives_less():
    # weights 0.9 and 0.1: each term is 0.09 * delta
    assert B.intensity_unc([9.0, 1.0], [0.0, 1.0], 0.5) == pytest.approx(
        0.5 * math.sqrt(2) * 0.09)


def test_the_origin_of_the_positions_does_not_matter():
    assert B.intensity_unc([3.0, 1.0, 2.0], [50000.1, 50000.4, 49999.9],
                           0.33) == pytest.approx(
        B.intensity_unc([3.0, 1.0, 2.0], [0.1, 0.4, -0.1], 0.33))


@pytest.mark.parametrize('intens, pos, u', [
    ([1.0], [0.0], 0.33),              # one component
    ([1.0, 1.0], [0.0, 1.0], 0.0),     # the term switched off
    ([0.0, 0.0], [0.0, 1.0], 0.33),    # intensities that say nothing
])
def test_nothing_to_add(intens, pos, u):
    assert B.intensity_unc(intens, pos, u) == 0.0


def test_the_project_configuration_sets_it():
    assert config.load().blend_u_ln == pytest.approx(0.33)


# ---------------------------------------------------------------------------
# classify_lines
# ---------------------------------------------------------------------------
def a_blend():
    levels = {x: EnergyLevel(level_id=x, energy=E, parity='e', J_str='7/2',
                             J_val=3.5)
              for x, E in ((LOW, 1000.0), (UP1, 31000.0), (UP2, 31000.6))}
    line = SpectralLine(30000.3, 0.02, 100.0, '', wn_key=30000.3)
    t1 = Transition(lower_level=levels[LOW], upper_level=levels[UP1],
                    assigned_to=line, accepted=1, calc_intensity=300.0,
                    grade='2A')
    t2 = Transition(lower_level=levels[LOW], upper_level=levels[UP2],
                    assigned_to=line, accepted=1, calc_intensity=100.0,
                    grade='2A')
    line.assigned_transitions.extend([t1, t2])
    return line, t1, t2


def test_calc_weights_widens_a_blend(monkeypatch):
    monkeypatch.setattr(cl, 'HFS', None)
    monkeypatch.setattr(cl, 'BLEND_U_LN', 0.33)
    line, t1, t2 = a_blend()
    u_I = B.intensity_unc([300.0, 100.0], [30000.0, 30000.6], 0.33)
    assert u_I > 0.02
    w = cl.calc_weights([line])
    assert w[id(t1)] == pytest.approx(0.75 ** 2 / (0.02 ** 2 + u_I ** 2))
    assert w[id(t2)] == pytest.approx(0.25 ** 2 / (0.02 ** 2 + u_I ** 2))
    assert line.fit_uncertainty == pytest.approx(math.hypot(0.02, u_I))
    # the table still says the measurement's uncertainty and the true BF
    df = cl.build_output([line], w).set_index('upp_id')
    assert df.loc[UP1, 'BF'] == pytest.approx(0.75)
    assert df.loc[UP2, 'BF'] == pytest.approx(0.25)
    assert df.loc[UP1, 'unc_wn_obs'] == pytest.approx(0.02)


def test_calc_weights_without_the_term(monkeypatch):
    monkeypatch.setattr(cl, 'HFS', None)
    monkeypatch.setattr(cl, 'BLEND_U_LN', 0.0)
    line, t1, t2 = a_blend()
    w = cl.calc_weights([line])
    assert w[id(t1)] == pytest.approx(0.75 ** 2 / 0.02 ** 2)
    assert line.fit_uncertainty == pytest.approx(0.02)


def test_a_single_line_is_not_widened(monkeypatch):
    monkeypatch.setattr(cl, 'HFS', None)
    monkeypatch.setattr(cl, 'BLEND_U_LN', 0.33)
    line, t1, t2 = a_blend()
    t2.accepted = 0
    w = cl.calc_weights([line])
    assert w[id(t1)] == pytest.approx(1.0 / 0.02 ** 2)
    assert line.fit_uncertainty == pytest.approx(0.02)


# ---------------------------------------------------------------------------
# make_LOPT_input
# ---------------------------------------------------------------------------
def row(wn, unc, low, upp, dif, accepted=1, calc=1.0):
    return {'wn_obs': wn, 'unc_wn_obs': unc, 'obs_intens': 100.0,
            'low_id': low, 'upp_id': upp, 'accepted': accepted,
            'calc_intens': calc, 'dif_wn_O-C': dif}


def read_back(path):
    out = {}
    for rec in io.open(path, encoding='ascii', newline=''):
        rec = rec.rstrip('\r\n')
        if rec.strip():
            got = {k: rec[a - 1:b].strip() for k, (a, b) in M.FIELDS.items()}
            out.setdefault(got['wavenumber'], []).append(
                (float(got['uncertainty']), got['flags'],
                 float(got['weight'] or 0)))
    return out


def test_centroid_unc_reads_the_positions_from_the_residuals():
    rows = [row('30000.3', '0.02', LOW, UP1, '0.3', calc=300.0),
            row('30000.3', '0.02', LOW, UP2, '-0.3', calc=100.0),
            row('30000.3', '0.02', LOW, '059003.000302', '0.9', accepted=0)]
    weights = M.blend_weights(rows[:2])
    assert M.centroid_unc(rows, weights, 0.33) == pytest.approx(
        B.intensity_unc([300.0, 100.0], [0.0, 0.6], 0.33))


def test_the_written_blend_is_widened(tmp_path):
    rows = [row('30000.3', '0.02', LOW, UP1, '0.3', calc=300.0),
            row('30000.3', '0.02', LOW, UP2, '-0.3', calc=100.0),
            row('25000.0', '0.02', LOW, '059003.000302', '0.01')]
    path = str(tmp_path / 'lines.txt')
    stats = {}
    M.write_lines_file(rows, path, {}, u_ln_blend=0.33, blend_stats=stats)
    got = read_back(path)
    u_I = B.intensity_unc([300.0, 100.0], [0.0, 0.6], 0.33)
    want = math.hypot(0.02, u_I)
    uncs = {u for u, _f, _w in got['30000.300']}
    assert len(uncs) == 1
    assert uncs.pop() == pytest.approx(want, abs=6e-4)
    assert sorted(w for _u, _f, w in got['30000.300']) == pytest.approx(
        [0.25, 0.75])
    assert got['25000.000'][0][0] == pytest.approx(0.02)
    assert stats == {'widened': 1, 'largest': pytest.approx(u_I)}


def test_off_by_default(tmp_path):
    rows = [row('30000.3', '0.02', LOW, UP1, '0.3', calc=300.0),
            row('30000.3', '0.02', LOW, UP2, '-0.3', calc=100.0)]
    path = str(tmp_path / 'lines.txt')
    M.write_lines_file(rows, path, {})
    assert {u for u, _f, _w in read_back(path)['30000.300']} == {0.02}
