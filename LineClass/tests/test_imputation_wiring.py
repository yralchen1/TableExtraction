"""Tests for the missing-gA policy as it is wired into the pipeline.

Run from the LineClass directory:  python -m pytest tests -q

`gA_imputation` decides *what value* stands for a transition that Cowan's
codes never printed; these tests are about *where the pipeline uses it*: that
policy "none" leaves the pipeline exactly as it was, that policy "impute"
gives every candidate pair absent from the calculated-transition file an
intensity built from the same model, and that such an imputed value is kept
out of the per-level intensity fits, where it would be mistaken for evidence
about the quality of the calculation.

Most of the tests replace the calibration with fixed numbers instead of
reading the real workbook, so they run on a machine that has only the code;
the two marked `real_file` check that the calibration the pipeline computes
for itself is the one `gA_imputation` reports.
"""
import math
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config                      # noqa: E402
import classify_lines as cl        # noqa: E402
import gA_imputation as imp        # noqa: E402
from models import EnergyLevel, SpectralLine, Transition   # noqa: E402

# The calibration used by the offline tests: (gA, u_ln, C, kT), the shape of
# what _imputation_calibration() returns.
FIXED = (175.8, 1.7383, 135.8, 12905.0)


@pytest.fixture
def policy():
    """Set the policy (and the calibration) for one test, then put it back."""
    saved = (cl.MISSING_POLICY, cl._IMPUTED)

    def _set(name, calibration=FIXED):
        cl.MISSING_POLICY = name
        cl._IMPUTED = calibration
    yield _set
    cl.MISSING_POLICY, cl._IMPUTED = saved


def levels():
    """One even and one odd level, 20000 cm^-1 apart: a candidate transition."""
    lo = EnergyLevel(level_id='lo', energy=0.0, parity='e', J_str='4', J_val=4.0)
    up = EnergyLevel(level_id='up', energy=20000.0, parity='o', J_str='5',
                     J_val=5.0)
    return lo, up


# --- the policy switch ------------------------------------------------------

def test_configuration_default_is_reproduced():
    """The module adopts the policy written in the configuration file."""
    cfg = config.load()
    saved = cl.MISSING_POLICY
    try:
        cl.apply_config(cfg)
        assert cl.MISSING_POLICY == cfg.missing_gA['policy']
        cl.apply_config(cfg, policy='impute')     # the --missing-gA override
        assert cl.MISSING_POLICY == 'impute'
    finally:
        cl.apply_config(cfg, policy=saved)


def test_an_unknown_policy_is_refused():
    cfg = config.load()
    saved = cl.MISSING_POLICY
    try:
        with pytest.raises(config.ConfigError, match='missing_gA.policy'):
            cl.apply_config(cfg, policy='guess')
    finally:
        cl.apply_config(cfg, policy=saved)


# --- the imputed value ------------------------------------------------------

def test_policy_none_imputes_nothing(policy):
    policy('none', None)
    lo, up = levels()
    assert cl.imputed_values(lo, up) == (None, None)


def test_policy_impute_follows_the_intensity_model(policy):
    policy('impute')
    lo, up = levels()
    gA, u_ln, C, kT = FIXED
    i_calc, u_calc = cl.imputed_values(lo, up)
    assert u_calc == u_ln
    assert i_calc == pytest.approx(
        C * gA * math.exp(-up.energy / kT) / (up.energy - lo.energy), rel=1e-14)


def test_a_transition_of_zero_wavenumber_is_not_imputed(policy):
    """The model divides by the wavenumber, so a degenerate pair has no value."""
    policy('impute')
    lo, up = levels()
    up.energy = lo.energy
    assert cl.imputed_values(lo, up) == (None, None)


# --- where the pipeline uses it ---------------------------------------------

def test_generated_transitions_are_untouched_under_none(policy):
    policy('none', None)
    lo, up = levels()
    got = cl.generate_all_possible_transitions([lo, up], {})
    assert len(got) == 1
    assert got[0].is_imputed == 0
    assert got[0].calc_intensity is None
    assert got[0].orig_calc_intensity is None


def test_generated_transitions_are_imputed_under_impute(policy):
    policy('impute')
    lo, up = levels()
    got = cl.generate_all_possible_transitions([lo, up], {})
    assert len(got) == 1
    tr = got[0]
    assert tr.is_imputed == 1
    assert tr.calc_intensity == pytest.approx(cl.imputed_values(lo, up)[0])
    # The imputation must fill the "original" fields too: the intensity
    # adjustment recomputes calc_intensity from them on every cycle.
    assert tr.orig_calc_intensity == tr.calc_intensity
    assert tr.orig_u_calc == tr.u_calc


def test_a_pair_present_in_the_file_is_never_imputed(policy):
    """The imputation applies to absence only, not to a weak calculated line."""
    policy('impute')
    lo, up = levels()
    index = {('lo', 'up'): {'calc_intensity': 0.5, 'u_calc': 0.3,
                            'assigned_to': None}}
    got = cl.generate_all_possible_transitions([lo, up], index)
    assert got[0].is_imputed == 0
    assert got[0].calc_intensity == 0.5


# --- what must not see it ---------------------------------------------------

def _accepted_transition(i_calc, u_calc, is_imputed):
    """An accepted transition on a line of its own, as the factor fit wants."""
    lo, up = levels()
    line = SpectralLine(wavenumber=20000.0, wn_uncertainty=0.01,
                        intensity=i_calc, line_character='')
    tr = Transition(lower_level=lo, upper_level=up, calc_intensity=i_calc,
                    u_calc=u_calc, assigned_to=line, accepted=1,
                    orig_calc_intensity=i_calc, orig_u_calc=u_calc,
                    is_imputed=is_imputed)
    line.assigned_transitions.append(tr)
    return tr


def test_imputed_transitions_stay_out_of_the_intensity_factors():
    """An imputed intensity is a bound, not a measurement: it must not be fitted.

    The third return value is the number of transitions that qualified for the
    fit, which is the number the imputed ones must be missing from.
    """
    real = [_accepted_transition(1.0 + 0.1 * i, 0.3, 0) for i in range(6)]
    fake = [_accepted_transition(1.0 + 0.1 * i, 1.7383, 1) for i in range(6)]
    _, _, n_real = cl._compute_factor_for_transition_list(real)
    _, _, n_both = cl._compute_factor_for_transition_list(real + fake)
    assert n_real == 6
    assert n_both == 6


# --- the real file ----------------------------------------------------------

@pytest.fixture
def real_cfg():
    cfg = config.load()
    if not os.path.isfile(cfg.icalc_file):
        pytest.skip(f"{cfg.icalc_file} is not present")
    return cfg


@pytest.mark.real_file
def test_the_pipeline_calibrates_itself_as_gA_imputation_does(real_cfg, policy):
    """The calibration the pipeline makes is the one the report quotes."""
    policy('impute', None)          # None: force a real calibration
    gA, u_ln, C, kT = cl._imputation_calibration()
    df = imp.load_icalc(real_cfg)
    mg = real_cfg.missing_gA
    want_gA, want_u = imp.estimate_missing_gA(
        df, real_cfg.gA_cutoff, mg['u_ln_window'], mg['u_ln_estimator'],
        mg['self_consistent'], mg['fit_range_decades'])
    assert (gA, u_ln) == pytest.approx((want_gA, want_u), rel=1e-15)
    assert (C, kT) == pytest.approx((real_cfg.intensity_model['C'],
                                     real_cfg.intensity_model['kT']), rel=1e-15)
    assert gA * math.exp(u_ln) == pytest.approx(real_cfg.gA_cutoff, rel=1e-12)


@pytest.mark.real_file
def test_the_calibration_is_computed_only_once(real_cfg, policy):
    """Reading the workbook again for every candidate pair would be ruinous."""
    policy('impute', None)
    first = cl._imputation_calibration()
    assert cl._imputation_calibration() is first
