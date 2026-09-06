"""Tests for `tools/calibrate_intensities.py`.

Run from the LineClass directory:  python -m pytest tests -q

The tool is standalone, so these tests import it directly from the tools
directory and work only on synthetic data: a spectrum is invented whose
intensity scale is deliberately spoiled by a known piecewise function, and the
tool is asked to find that function back.
"""
import math
import os
import sys

import numpy as np
import pytest

TOOLS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     'tools')
sys.path.insert(0, TOOLS)

import calibrate_intensities as cal    # noqa: E402


# ---------------------------------------------------------------------------
# the file format
# ---------------------------------------------------------------------------


def test_correction_file_round_trip(tmp_path):
    """Writing regions and reading them back must give the same function."""
    path = str(tmp_path / 'corr.txt')
    regions = [
        cal.Region(1000.0, 2000.0, np.polynomial.Polynomial([1.0, -1e-3])),
        cal.Region(2000.0, 4000.0, np.polynomial.Polynomial([-2.0, 2e-3, -1e-7])),
    ]
    cal.write_correction(path, regions)
    back = cal.read_correction(path)
    lam = np.linspace(1000.0, 4000.0, 500)
    assert np.allclose(cal.apply_correction(regions, lam),
                       cal.apply_correction(back, lam), atol=1e-9)


def test_outside_the_regions_is_clamped_not_extrapolated():
    """A wavelength beyond the ends takes the value at the nearest end."""
    r = cal.Region(1000.0, 2000.0, np.polynomial.Polynomial([0.0, 1e-3]))
    out = cal.apply_correction([r], np.array([500.0, 1500.0, 5000.0]))
    assert out[0] == pytest.approx(r(1000.0))
    assert out[1] == pytest.approx(1.5)
    assert out[2] == pytest.approx(r(2000.0))


# ---------------------------------------------------------------------------
# the model
# ---------------------------------------------------------------------------


def test_calculated_intensity_is_an_energy_flux():
    """Icalc must be proportional to the wavenumber, not to the wavelength.

    Doubling the wavenumber of a transition, everything else held fixed, must
    double its calculated intensity: the detector records energy, and a photon
    of twice the wavenumber carries twice the energy.
    """
    a = cal.calc_intensity(1.0, 1e4, 1e6, 20000.0, 30000.0)
    b = cal.calc_intensity(1.0, 1e4, 1e6, 40000.0, 30000.0)
    assert b == pytest.approx(2 * a)


def test_joint_fit_recovers_the_constants_with_no_correction():
    rng = np.random.default_rng(0)
    n = 400
    Eup = rng.uniform(2e4, 1.2e5, n)
    rwn = rng.uniform(1e4, 1e5, n)
    gA = 10 ** rng.uniform(4, 8, n)
    C, kT = 130.0, 12000.0
    d = dict(Eup=Eup, rwn=rwn, gA=gA, lam=1e8 / rwn,
             I0=cal.calc_intensity(C, kT, gA, rwn, Eup))
    Cf, kTf, P, regions = cal.joint_fit(d, np.ones(n, bool), [])
    assert Cf == pytest.approx(C, rel=1e-8)
    assert kTf == pytest.approx(kT, rel=1e-8)
    assert regions == []
    assert np.allclose(P, 0.0)


# ---------------------------------------------------------------------------
# finding a known correction back
# ---------------------------------------------------------------------------


def synthetic(n=1500, seed=1, scatter=0.15):
    """A spectrum whose intensity scale is spoiled by a known function.

    The truth is two regions with a jump between them: a straight ramp below
    2000 A and a downward one above it.  The tool has to find both the
    breakpoint and the two polynomials, knowing neither.
    """
    rng = np.random.default_rng(seed)
    lam = np.sort(rng.uniform(1000.0, 4000.0, n))
    rwn = 1e8 / lam
    Eup = rwn + rng.uniform(0.0, 2e4, n)
    gA = 10 ** rng.uniform(4, 8, n)
    C, kT = 100.0, 12000.0
    Icalc = cal.calc_intensity(C, kT, gA, rwn, Eup)
    truth = np.where(lam < 2000.0, 1.0 - 1e-3 * lam, -3.0 + 6e-4 * lam)
    # Iobs is what the plate reported: the true intensity spoiled by the
    # scale error, plus a little scatter.
    I0 = Icalc * np.exp(-truth + rng.normal(0.0, scatter, n))
    d = dict(lam=lam, rwn=rwn, Eup=Eup, gA=gA, I0=I0, Iraw=I0,
             Elow=Eup - rwn, u=np.full(n, 10.0), wn=rwn)
    return d, truth, C, kT


def test_the_known_correction_is_recovered():
    d, truth, C, kT = synthetic()
    opt = cal.parse_args(['--no-gaps'])
    use = np.ones(len(d['lam']), bool)
    fit = cal.calibrate(d, use, opt, opt.tol_final, lambda m: None)
    assert fit.converged
    # the correction is defined up to an overall constant, which the fit
    # levels to zero, so compare the two after removing their means
    got = fit.P - fit.P.mean()
    want = truth - truth.mean()
    assert np.max(np.abs(got - want)) < 0.5
    assert fit.kT == pytest.approx(kT, rel=0.15)
    # the breakpoint at 2000 A must be found to within a few angstroms
    cuts = [r.lo for r in fit.regions[1:]]
    assert any(abs(c - 2000.0) < 20.0 for c in cuts)


def test_a_smooth_scale_needs_only_one_region():
    """With no jump and a gentle slope, no breakpoint should be invented."""
    d, truth, C, kT = synthetic()
    d['I0'] = d['I0'] * np.exp(truth)          # undo the spoiling
    opt = cal.parse_args(['--no-gaps'])
    fit = cal.calibrate(d, np.ones(len(d['lam']), bool), opt, opt.tol_final,
                        lambda m: None)
    assert len(fit.regions) == 1
    assert np.max(np.abs(fit.P - fit.P.mean())) < 0.3


def test_no_polynomial_runs_away_inside_its_region():
    """Every fitted polynomial must stay near the points that constrain it."""
    d, truth, C, kT = synthetic()
    opt = cal.parse_args(['--no-gaps'])
    fit = cal.calibrate(d, np.ones(len(d['lam']), bool), opt, opt.tol_final,
                        lambda m: None)
    target = fit.dln + fit.P
    where = cal.region_of(fit.regions, d['lam'])
    for k, r in enumerate(fit.regions):
        m = where == k
        g = np.linspace(r.lo, r.hi, 300)
        p = r(g)
        assert p.max() - np.percentile(target[m], 99) < opt.max_swing + 1e-9
        assert np.percentile(target[m], 1) - p.min() < opt.max_swing + 1e-9


# ---------------------------------------------------------------------------
# the coverage gaps
# ---------------------------------------------------------------------------


def test_a_gap_forces_a_boundary_even_where_the_scale_is_smooth():
    """A gap is a boundary by decree, not by evidence.

    The scale here is perfectly smooth, so nothing in the data asks for a cut;
    the gap must produce one all the same, with the regions ending exactly at
    its edges and nothing covering the gap itself.
    """
    d, truth, C, kT = synthetic()
    d['I0'] = d['I0'] * np.exp(truth)          # undo the spoiling: no jump
    opt = cal.parse_args(['--gap', '2500', '2600'])
    fit = cal.calibrate(d, np.ones(len(d['lam']), bool), opt, opt.tol_final,
                        lambda m: None)
    assert len(fit.regions) == 2
    assert fit.regions[0].hi == pytest.approx(2500.0)
    assert fit.regions[1].lo == pytest.approx(2600.0)
    # no region covers the gap
    for r in fit.regions:
        assert not (r.lo < 2550.0 < r.hi)


def test_no_region_spans_a_gap():
    d, truth, C, kT = synthetic()
    opt = cal.parse_args([])                   # the two Pr III gaps
    fit = cal.calibrate(d, np.ones(len(d['lam']), bool), opt, opt.tol_final,
                        lambda m: None)
    for lo, hi in opt.gaps:
        for r in fit.regions:
            assert not (r.lo < hi and lo < r.hi), (r.lo, r.hi, lo, hi)


def test_a_jump_at_a_gap_is_reproduced_however_large():
    """Across a gap the two sides are fitted independently, so a step of any
    size is followed exactly - not smeared over the neighbourhood."""
    d, truth, C, kT = synthetic()
    step = 3.0                                 # a factor 20 in intensity
    extra = np.where(d['lam'] > 2600.0, step, 0.0)
    d['I0'] = d['I0'] * np.exp(truth - extra)  # smooth scale, plus a big step
    opt = cal.parse_args(['--gap', '2500', '2600'])
    fit = cal.calibrate(d, np.ones(len(d['lam']), bool), opt, opt.tol_final,
                        lambda m: None)
    # the lines that happen to fall inside the invented gap are not part of
    # any region (a real gap holds no lines at all), so they are left out
    out = (d['lam'] < 2500.0) | (d['lam'] > 2600.0)
    got = (fit.P - extra)[out]
    assert np.max(np.abs(got - got.mean())) < 0.5
    # the correction really does step across the gap
    below = fit.regions[0](fit.regions[0].hi)
    above = fit.regions[1](fit.regions[1].lo)
    assert above - below == pytest.approx(step, abs=0.3)


def test_lines_inside_a_gap_are_left_out_of_the_fits():
    d, truth, C, kT = synthetic()
    said = []
    opt = cal.parse_args(['--gap', '2000', '2600'])   # data does live there
    cal.build_regions(d['lam'], truth, opt, None, said.append)
    assert any('inside the coverage gap' in m for m in said)


def test_gaps_are_validated():
    with pytest.raises(SystemExit):
        cal.parse_args(['--gap', '2600', '2500'])
    with pytest.raises(SystemExit):
        cal.parse_args(['--gap', '2000', '2600', '--gap', '2500', '2700'])


# ---------------------------------------------------------------------------
# undoing an existing correction
# ---------------------------------------------------------------------------


def test_restore_undoes_a_correction(tmp_path):
    """Dividing out the function that was applied gives the originals back."""
    path = str(tmp_path / 'corr.txt')
    regions = [cal.Region(1000.0, 4000.0,
                          np.polynomial.Polynomial([1.0, -5e-4]))]
    cal.write_correction(path, regions)
    lam = np.linspace(1000.0, 4000.0, 300)
    original = np.round(np.linspace(1, 500, 300))
    d = dict(lam=lam, Iraw=1000.0 * original * np.exp(cal.apply_correction(regions, lam)))
    opt = cal.parse_args(['--restore-from', path])
    got = cal.restore_originals(d, opt, lambda m: None)
    assert np.allclose(got, original)


# ---------------------------------------------------------------------------
# the edge of the spectrum: extrapolation
# ---------------------------------------------------------------------------


def test_extrapolation_is_measured_from_the_outermost_fitted_point():
    """_extrapolation reports how far P moves past the last line, not the scatter."""
    x = np.linspace(1000.0, 2000.0, 50)
    poly = np.polynomial.Polynomial.fit(x, 1e-6 * (x - 1500.0) ** 2, 2)
    r = cal.Region(1000.0, 2000.0, poly)
    assert cal._extrapolation(r, x) == pytest.approx(0.0, abs=1e-9)
    r_stretched = cal.Region(500.0, 2000.0, poly)
    moved = cal._extrapolation(r_stretched, x)
    assert moved == pytest.approx(abs(poly(500.0) - poly(1000.0)), rel=1e-6)
    assert moved > 0.5


def test_a_stretched_region_is_flattened_until_it_stops_moving():
    """The degree of a stretched region falls until --max-extrap is satisfied.

    The points say "a parabola"; the region has to cover far more than they
    span, and out there the parabola climbs steeply.  With the guard on, the
    degree must come down; with it off, the parabola survives - which is the
    behaviour that let a cubic dive by three units below 900 A.
    """
    x = np.linspace(1000.0, 2000.0, 60)
    y = 1e-6 * (x - 2000.0) ** 2
    poly, _ = cal._fit(x, y, 2)
    for max_extrap, expect in ((0.1, 1), (float('inf'), 2)):
        r = cal.Region(200.0, 2000.0, poly, len(x))
        cal._tame(r, x, y, max_swing=10.0, max_extrap=max_extrap)
        assert r.degree <= expect
    r = cal.Region(200.0, 2000.0, poly, len(x))
    cal._tame(r, x, y, max_swing=10.0, max_extrap=0.1)
    assert cal._extrapolation(r, x) <= 0.1 + 1e-9


def test_the_outlier_cut_is_a_schedule_by_default():
    """--dln-cut takes a list; one value still means a fixed cut."""
    assert cal.parse_args(['--no-gaps']).dln_cut == [4.0, 3.0, 2.0]
    assert cal.parse_args(['--no-gaps', '--dln-cut', '2']).dln_cut == [2.0]
    assert cal.parse_args(['--no-gaps', '--dln-cut', '5', '2']).dln_cut == [5.0, 2.0]
