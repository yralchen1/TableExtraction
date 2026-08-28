"""Tests for `gA_imputation`.

Run from the LineClass directory:  python -m pytest tests -q

The tests fall into two groups.  Most of them work on a small synthetic table
and check the arithmetic and the failure modes; a few, marked `real_file`,
open the real calculated-transition workbook named in the configuration and
check the numbers the pipeline will actually use.  The second group is
skipped when that workbook is not present, so the first group still runs on a
machine that has only the code.
"""
import math
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config                      # noqa: E402
import gA_imputation as imp        # noqa: E402

CUTOFF = 1.0e3


def synthetic(n_per_decade=500, seed=0):
    """A table shaped like the real one: u_ln falling linearly with log10(gA).

    Every decade from 10^3 to 10^8 gets `n_per_decade` rows whose u_ln is
    -0.25*log10(gA) + 2.5 plus a little symmetric scatter, and whose Icalc
    follows the intensity model exactly with C = 100, kT = 10000 cm^-1.
    """
    rng = np.random.default_rng(seed)
    x = np.concatenate([rng.uniform(d, d + 1, n_per_decade)
                        for d in range(3, 8)])
    gA = 10.0 ** x
    u_ln = -0.25 * x + 2.5 + rng.normal(0.0, 0.05, x.size)
    Eup = rng.uniform(10000.0, 60000.0, x.size)
    rwn = rng.uniform(9000.0, 100000.0, x.size)
    Icalc = 100.0 * gA * np.exp(-Eup / 10000.0) / rwn
    return pd.DataFrame({'gA': gA, 'u_ln': u_ln, 'Icalc': Icalc,
                         'Eup': Eup, 'rwn': rwn})


# --- the defining relation ---------------------------------------------------

@pytest.mark.parametrize('estimator', ['rms', 'mean', 'median'])
def test_plain_estimate_sits_on_the_cutoff(estimator):
    """gA*exp(u_ln) must land exactly on the censoring threshold."""
    df = synthetic()
    gA, u_ln = imp.estimate_missing_gA(df, CUTOFF, 10.0, estimator, False)
    assert abs(gA * math.exp(u_ln) / CUTOFF - 1.0) < 1e-12


@pytest.mark.parametrize('estimator', ['rms', 'mean', 'median'])
def test_self_consistent_estimate_sits_on_the_cutoff(estimator):
    df = synthetic()
    gA, u_ln = imp.estimate_missing_gA(df, CUTOFF, 10.0, estimator, True,
                                       (3.0, 6.0))
    assert abs(gA * math.exp(u_ln) / CUTOFF - 1.0) < 1e-12


def test_plain_estimate_uses_the_window_only():
    """The plain estimate is the chosen summary of the window, nothing else."""
    df = synthetic()
    sel = (df['gA'] >= CUTOFF) & (df['gA'] <= 10.0 * CUTOFF)
    expected = float(np.sqrt(np.mean(df.loc[sel, 'u_ln'].to_numpy() ** 2)))
    _, u_ln = imp.estimate_missing_gA(df, CUTOFF, 10.0, 'rms', False)
    assert u_ln == pytest.approx(expected, rel=0, abs=1e-15)


# --- the self-consistent solution -------------------------------------------

def test_self_consistent_solution_is_a_fixed_point():
    """The returned pair satisfies both defining equations at once."""
    df = synthetic()
    slope, intercept = imp.fit_u_ln(df, 'rms', (3.0, 6.0))
    gA, u_ln = imp.estimate_missing_gA(df, CUTOFF, 10.0, 'rms', True,
                                       (3.0, 6.0))
    assert u_ln == pytest.approx(slope * math.log10(gA) + intercept, abs=1e-12)
    assert gA == pytest.approx(CUTOFF * math.exp(-u_ln), rel=1e-15)


def test_self_consistent_estimate_is_below_the_plain_one():
    """Extrapolating a rising u_ln must push the imputed gA further down."""
    df = synthetic()
    gA_plain, u_plain = imp.estimate_missing_gA(df, CUTOFF, 10.0, 'rms', False)
    gA_sc, u_sc = imp.estimate_missing_gA(df, CUTOFF, 10.0, 'rms', True,
                                          (3.0, 6.0))
    assert u_sc > u_plain
    assert gA_sc < gA_plain


def test_self_consistent_converges_from_any_start():
    """The iteration is a contraction: it does not depend on where it starts.

    Solving the fixed-point equation directly gives the same answer as the
    iteration, so the iteration has genuinely converged rather than stopped.
    """
    df = synthetic()
    slope, intercept = imp.fit_u_ln(df, 'rms', (3.0, 6.0))
    # log10(gA) = log10(cutoff) - (slope*log10(gA) + intercept)/ln(10)
    ln10 = math.log(10.0)
    exact_lg = (math.log10(CUTOFF) - intercept / ln10) / (1.0 + slope / ln10)
    gA, _ = imp.estimate_missing_gA(df, CUTOFF, 10.0, 'rms', True, (3.0, 6.0))
    assert math.log10(gA) == pytest.approx(exact_lg, abs=1e-12)


def test_too_steep_a_slope_is_refused():
    """A slope steeper than ln(10) per decade has no unique solution."""
    x = np.repeat(np.arange(3.0, 8.0) + 0.5, 100)
    df = pd.DataFrame({'gA': 10.0 ** x, 'u_ln': -3.0 * x + 30.0})
    with pytest.raises(imp.ImputationError, match='per decade'):
        imp.estimate_missing_gA(df, CUTOFF, 10.0, 'rms', True, (3.0, 6.0))


# --- what the estimate must ignore ------------------------------------------

def test_rows_below_the_cutoff_are_left_untouched():
    """Rows with gA < cutoff must not move either estimate.

    They are the few rows that slipped below the printing threshold when the
    wavenumbers were rescaled, not a sample of the censored population, so
    adding wildly wrong ones must change nothing.
    """
    df = synthetic()
    intruders = pd.DataFrame({'gA': [1.0, 10.0, 100.0, 999.0],
                              'u_ln': [50.0, 50.0, 50.0, 50.0],
                              'Icalc': [1.0] * 4, 'Eup': [1.0] * 4,
                              'rwn': [1.0] * 4})
    spiked = pd.concat([df, intruders], ignore_index=True)
    for self_consistent in (False, True):
        a = imp.estimate_missing_gA(df, CUTOFF, 10.0, 'rms', self_consistent,
                                    (3.0, 6.0))
        b = imp.estimate_missing_gA(spiked, CUTOFF, 10.0, 'rms',
                                    self_consistent, (3.0, 6.0))
        assert a == pytest.approx(b, rel=1e-15)


def test_an_empty_window_is_an_error_not_a_guess():
    df = synthetic()
    with pytest.raises(imp.ImputationError, match='too few'):
        imp.estimate_missing_gA(df, 1.0e12, 10.0, 'rms', False)


def test_an_unknown_estimator_is_refused():
    with pytest.raises(imp.ImputationError, match='unknown estimator'):
        imp.estimate_missing_gA(synthetic(), CUTOFF, 10.0, 'mode', False)


# --- the intensity model ----------------------------------------------------

def test_intensity_model_is_recovered_from_synthetic_data():
    df = synthetic()
    C, kT = imp.fit_intensity_model(df)
    assert C == pytest.approx(100.0, rel=1e-8)
    assert kT == pytest.approx(10000.0, rel=1e-8)
    assert imp.intensity_model_error(df, C, kT) < 1e-8


def test_impute_intensity_matches_the_definition():
    got = imp.impute_intensity(175.8, 40000.0, 30000.0, 135.8, 12905.0)
    want = 135.8 * 175.8 * math.exp(-40000.0 / 12905.0) / 30000.0
    assert float(got) == pytest.approx(want, rel=1e-14)


def test_impute_intensity_is_vectorised():
    gA = np.array([100.0, 200.0])
    got = imp.impute_intensity(gA, np.array([1e4, 2e4]),
                               np.array([1e4, 2e4]), 135.8, 12905.0)
    assert got.shape == (2,)
    assert got[0] == pytest.approx(
        135.8 * 100.0 * math.exp(-1e4 / 12905.0) / 1e4, rel=1e-14)


# --- the real file ----------------------------------------------------------

@pytest.fixture(scope='module')
def real():
    """(configuration, calculated transitions) or a skip if the file is gone."""
    cfg = config.load()
    if not os.path.isfile(cfg.icalc_file):
        pytest.skip(f"{cfg.icalc_file} is not present")
    return cfg, imp.load_icalc(cfg)


@pytest.mark.real_file
def test_real_file_default_estimate(real):
    """The calibrated numbers the pipeline will use: 176 s^-1 and 1.738."""
    cfg, df = real
    gA, u_ln = imp.estimate_missing_gA(df, cfg.gA_cutoff, 10.0, 'rms', False)
    assert u_ln == pytest.approx(1.738, abs=0.001)
    assert gA == pytest.approx(176.0, abs=0.5)


@pytest.mark.real_file
def test_real_file_self_consistent_estimate(real):
    """The extrapolated variant, quoted in the plan as 114-148 / 1.91-2.17."""
    cfg, df = real
    gA, u_ln = imp.estimate_missing_gA(df, cfg.gA_cutoff, 10.0, 'rms', True,
                                       (3.0, 6.0))
    assert 114.0 <= gA <= 148.0
    assert 1.91 <= u_ln <= 2.17


@pytest.mark.real_file
def test_real_file_follows_the_intensity_model(real):
    """Icalc is reproduced from gA, Eup and rwn to better than 1e-5."""
    _, df = real
    C, kT = imp.fit_intensity_model(df)
    assert imp.intensity_model_error(df, C, kT) < 1e-5


@pytest.mark.real_file
def test_real_file_agrees_with_the_configured_model(real):
    """The C and kT written in the configuration still describe the file."""
    cfg, df = real
    C, kT, C_fit, kT_fit, agrees = imp.check_intensity_model(df, cfg,
                                                             stream=None)
    assert agrees, (f"configured C={C:g}, kT={kT:g}; "
                    f"fitted C={C_fit:g}, kT={kT_fit:g}")


@pytest.mark.real_file
def test_report_runs_and_names_its_numbers(real):
    """The --report output is produced and quotes the calibrated values."""
    import io
    cfg, _ = real
    buf = io.StringIO()
    imp.report(cfg, stream=buf)
    text = buf.getvalue()
    for expected in ('gA_missing', 'u_ln_missing', 'Intensity model',
                     '1.7383', '175.8'):
        assert expected in text
