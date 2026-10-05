"""Regression tests for the convention fit of `hfs_kappa.py`.

Two kinds of test live here.  The first kind is arithmetic: the degeneracy
theorem the whole reading of kappa rests on, the conversion of Sugar's stated
uncertainties, and the solver's ability to recover a kappa that was put into
synthetic data by hand.  Those are exact and never move.

The second kind runs the fit on the project's own line list, and its numbers
are the ones `Work_on_hfs_plan.md` sections 2.2 and 3.3 quote.  That file is
regenerated whenever the assignment work moves, so these assertions are stated
as the ranges the conclusions need rather than as exact values: what must not
change without being noticed is that the unflagged lines are measured near
their strongest component, that the two eras disagree, and that the
displacements carry the signal - not the last digit of any of them.
"""

import math
import os

import pytest

import hfs_correction
import hfs_kappa


@pytest.fixture(scope='module')
def fitted():
    lines = hfs_kappa.read_lines()
    sigma, keep, table, fit = hfs_kappa.adopt_uncertainties(lines)
    return lines, sigma, keep, table, fit


# --- the arithmetic ---------------------------------------------------------
def test_stated_uncertainty_converts_angstrom_to_wavenumber():
    """u_wn = u_lambda * wn^2, with u_lambda in angstrom and wn in cm^-1."""
    # a plain 1974 line at 20000 cm^-1 (5000 angstrom), stated 0.0030 angstrom
    assert hfs_kappa.stated_uncertainty('', 1974, 20000.0) == pytest.approx(
        0.0030 * 1e-8 * 20000.0 ** 2)
    assert hfs_kappa.stated_uncertainty('', 1974, 20000.0) == pytest.approx(
        0.012, abs=1e-6)
    # the same stated wavelength error is worth far more at short wavelength
    assert hfs_kappa.stated_uncertainty('', 1969, 90000.0) == pytest.approx(
        0.324, abs=1e-3)
    # every character but the blank one starts from the stated average
    # deviation from Ritz, 0.0070 angstrom
    assert hfs_kappa.stated_uncertainty('w', 1974, 20000.0) == pytest.approx(
        0.0070 * 1e-8 * 20000.0 ** 2)


def test_only_calculated_constants_are_scaled():
    """Step 2's scale error belongs to the calculated A constants; a constant
    measured from flags or from resolved components is already on the
    measured scale and is read as written."""
    s, u_s = hfs_kappa.A_SCALE, hfs_kappa.U_A_SCALE
    for src in ('composition', 'Reader & Sugar 1965 (calculated)',
                'composition, flags consistent (2026-09-29 review)'):
        assert hfs_kappa.is_calculated(src)
        assert hfs_kappa.scaled_A(0.2, 0.01, src) == pytest.approx(
            (0.2 * s, math.hypot(0.01 * s, 0.2 * u_s)))
    for src in ('flag interval', 'not determined',
                'estimate: unsplit lines (2026-10-01 review)',
                'center of gravity of 35378.125 (1116-991), 2026-10-02',
                'rungs 1 and 2 of 1228-1116 (21962.950 and 21963.279), '
                '2026-10-02'):
        assert not hfs_kappa.is_calculated(src)
        assert hfs_kappa.scaled_A(0.2, 0.01, src) == (0.2, 0.01)


def test_the_project_table_is_read_on_the_measured_scale():
    """Every row of A_hfs_levels.csv comes back scaled by its source."""
    import csv
    path = os.path.join(os.path.dirname(hfs_kappa.__file__),
                        hfs_kappa.A_LEVELS)
    with open(path, encoding='utf-8', newline='') as fh:
        rows = {r['level_id']: r for r in csv.DictReader(fh)}
    read = hfs_kappa.read_A_constants()
    assert set(read) == set(rows)
    for lid, r in rows.items():
        A = read[lid][1]
        if not hfs_correction.is_determined(r['source']):
            assert A == 0.0
        elif hfs_kappa.is_calculated(r['source']):
            assert A == pytest.approx(hfs_kappa.A_SCALE * float(r['A_cm-1']))
        else:
            assert A == float(r['A_cm-1'])


def test_era_and_class_of_a_line():
    assert hfs_kappa.era_of(hfs_kappa.ERA_SPLIT - 1) == 1974
    assert hfs_kappa.era_of(hfs_kappa.ERA_SPLIT) == 1969
    assert hfs_kappa.kappa_class('*r') == 'flag'
    assert hfs_kappa.kappa_class('*v') == 'flag'
    assert hfs_kappa.kappa_class('c') == 'c'
    # `ch` is a character of its own and is not the complex class
    assert hfs_kappa.kappa_class('ch') == 'plain'
    assert hfs_kappa.kappa_class('') == 'plain'


def _toy(kappa_plain=0.5):
    """Four levels, six lines, exact data made with a known kappa."""
    S = {'a': 0.0, 'b': 0.30, 'c': -0.20, 'd': 0.10}
    E = {'a': 0.0, 'b': 1000.0, 'c': 2500.0, 'd': 4000.0}
    pairs = [('a', 'b', '*r'), ('a', 'c', '*r'), ('b', 'd', '*v'),
             ('a', 'd', ''), ('b', 'c', ''), ('c', 'd', '')]
    lines = []
    for low, upp, char in pairs:
        D = S[upp] - S[low]
        k = 1.0 if char else kappa_plain
        wn = E[upp] - E[low] + k * D
        lines.append(hfs_kappa.Line(
            wn=wn, char=char, era=1974, cls=hfs_kappa.kappa_class(char),
            ucls=(char, 1974), low=low, upp=upp, D=D, dJ=1.0))
    return lines, S, E


def test_the_solver_recovers_a_kappa_put_in_by_hand():
    lines, _, _ = _toy(kappa_plain=0.5)
    fit = hfs_kappa.solve(lines, lambda ln: ln.cls, 'flag', 1.0,
                          sigma=[0.01] * len(lines))
    assert fit['kappa']['plain'] == pytest.approx(0.5, abs=1e-9)
    assert fit['rms'] == pytest.approx(0.0, abs=1e-9)


def test_the_degeneracy_theorem_holds_exactly():
    """Adding c to every kappa and c*S to every level changes nothing.

    This is why one class has to be anchored from outside the line list, and
    it is the reason Step 2 - which measures kappa(flag) from Sugar's own
    resolved components - is what the whole correction rests on.
    """
    lines, S, E = _toy(kappa_plain=0.5)
    for c in (0.25, 0.5, 1.0):
        for ln in lines:
            kappa = 1.0 if ln.cls == 'flag' else 0.5
            # move every level by c times its own S and every kappa by -c
            low = E[ln.low] + c * S[ln.low]
            upp = E[ln.upp] + c * S[ln.upp]
            assert ln.wn == pytest.approx(
                (upp - low) + (kappa - c) * ln.D, abs=1e-12)
    # and the solver returns the shifted kappa when the anchor is shifted
    sigma = [0.01] * len(lines)
    for anchor in (1.0, 0.75, 0.5):
        fit = hfs_kappa.solve(lines, lambda ln: ln.cls, 'flag', anchor,
                              sigma=sigma)
        assert fit['kappa']['plain'] == pytest.approx(
            0.5 - (1.0 - anchor), abs=1e-9)
        assert fit['rms'] == pytest.approx(0.0, abs=1e-9)


def test_three_free_kappas_are_rank_deficient_by_one():
    """The same theorem, as the solver sees it.

    With every class free the system has exactly one null direction; anchoring
    one class removes it and nothing else.
    """
    lines, _, _ = _toy()
    sigma = [0.01] * len(lines)
    free_all = hfs_kappa.solve(lines, lambda ln: ln.cls + '_free', 'none',
                               1.0, sigma=sigma)
    anchored = hfs_kappa.solve(lines, lambda ln: ln.cls, 'flag', 1.0,
                               sigma=sigma)
    assert free_all['rank'] == free_all['unknowns'] - 2   # energy zero + kappa
    assert anchored['rank'] == anchored['unknowns'] - 1   # energy zero only


def test_leverage_is_one_where_a_line_stands_alone():
    """A line that alone determines a level cannot measure its own scatter."""
    lines, _, _ = _toy()
    alone = hfs_kappa.Line(wn=9000.0, char='', era=1974, cls='plain',
                           ucls=('', 1974), low='a', upp='z', D=0.0, dJ=0.0)
    fit = hfs_kappa.solve(lines + [alone], lambda ln: ln.cls, 'flag', 1.0,
                          sigma=[0.01] * (len(lines) + 1))
    assert fit['leverage'][-1] == pytest.approx(1.0, abs=1e-9)
    assert fit['residual'][-1] == pytest.approx(0.0, abs=1e-9)


# --- the fit on the project's own line list ---------------------------------
def test_the_unflagged_lines_are_not_measured_at_the_center_of_gravity(fitted):
    """The result of Step 3, and it contradicts the hypothesis of section 1.1.

    kappa(plain) = 0 would mean Sugar measured the center of gravity of an
    unresolved pattern.  He did not: the unflagged lines sit within a tenth of
    the way from the strongest component to the center of gravity.
    """
    _, _, _, _, fit = fitted
    assert 0.88 < fit['kappa']['plain'] < 0.98
    assert fit['u_kappa']['plain'] < 0.03
    assert fit['kappa']['plain'] / fit['u_kappa']['plain'] > 20


def test_the_complex_lines_are_measured_near_the_center_of_gravity(fitted):
    _, _, _, _, fit = fitted
    assert 0.0 < fit['kappa']['c'] < 0.45


def test_the_fit_is_rank_deficient_only_in_the_energy_zero(fitted):
    _, _, _, _, fit = fitted
    assert fit['rank'] == fit['unknowns'] - 1
    assert fit['chi2'] / fit['dof'] == pytest.approx(1.0, abs=0.05)


def test_the_two_eras_disagree(fitted):
    """Section 3.3's first test, and the model does not survive it.

    A single kappa for the unflagged lines is refused by the data: the lines
    Sugar measured in 1969, above the era boundary, sit markedly closer to the
    center of gravity than the 1974 ones.  Step 5 must therefore carry one
    kappa per era, not one per class.
    """
    lines, sigma, keep, _, _ = fitted
    fit = hfs_kappa.era_test(lines, sigma, keep)
    d, u = hfs_kappa.difference(fit, ('plain', 1974), ('plain', 1969))
    assert fit['kappa'][('plain', 1974)] > fit['kappa'][('plain', 1969)]
    assert d / u > 5.0


def test_the_era_difference_is_not_an_artifact_of_the_constants(fitted):
    """The same split, confined to the levels both eras have in common.

    If the two eras differed only because they see different levels with
    differently reliable A constants, this would wash the difference out.  It
    does not.
    """
    lines, sigma, keep, _, _ = fitted
    fit, shared, counts = hfs_kappa.shared_level_test(lines, sigma, keep)
    assert len(shared) > 50
    d, u = hfs_kappa.difference(fit, ('plain', 1974), ('plain', 1969))
    assert d / u > 4.0


def test_the_flag_branches_agree(fitted):
    """Section 3.3's second test, which the model does survive.

    Step 2 showed the extreme component is the strongest one for every J pair
    the line list contains, so no per-branch geometric factor is needed; the
    residuals agree, within about two sigma of the dJ = +1 branch.
    """
    lines, sigma, keep, _, _ = fitted
    fit = hfs_kappa.branch_test(lines, sigma, keep)
    values = [fit['kappa'][k] for k in fit['kappa'] if isinstance(k, tuple)]
    assert len(values) == 3
    assert max(values) - min(values) < 0.15


def test_permuting_the_displacements_destroys_kappa(fitted):
    """The control: kappa(plain) must be produced by hyperfine structure.

    Shuffling the displacements among the plain lines leaves the level values,
    the weights and the class sizes alone and breaks only the pairing of a
    line with its own pattern.  kappa then comes out at zero.
    """
    lines, sigma, keep, _, fit = fitted
    draws = [hfs_kappa.null_test(lines, sigma, keep,
                                 seed=s)['kappa']['plain']
             for s in range(4)]
    mean = sum(draws) / len(draws)
    assert abs(mean) < 0.05
    assert fit['kappa']['plain'] - mean > 0.5


def test_the_adopted_uncertainties_are_sane(fitted):
    """The Step 1 table: one value per character per era, never below the floor.

    The 1969 values are an order of magnitude larger than the 1974 ones for
    the same stated wavelength error, because the same error in angstrom is
    worth wn^2 more in cm^-1 at short wavelength; that is arithmetic, not a
    statement about Sugar's plates.
    """
    _, sigma, _, table, _ = fitted
    assert all(u.a >= 0.0 and u.b >= 0.0 for u in table.values())
    assert all(u.n >= 1 for u in table.values())
    assert all(x >= hfs_kappa.FLOOR for x in sigma)
    plain74 = table[('', 1974)].adopted
    plain69 = table[('', 1969)].adopted
    assert 0.02 < plain74 < 0.08
    assert plain69 > 5 * plain74


def test_the_wavelength_term_is_smaller_than_Sugar_stated(fitted):
    """`a` is the part of the uncertainty that is a distance on the plate.

    Sugar states 0.0030 angstrom for a plain line of 1974 and 0.0040 for one
    of 1969, and those are total uncertainties - the plate-reading error with
    everything else already added in quadrature.  The term fitted here is the
    plate-reading error alone, so it has to come out below what he states.
    """
    _, _, _, table, _ = fitted
    assert 0.0 < table[('', 1974)].a < hfs_kappa.STATED_PLAIN[1974]
    assert 0.0 < table[('', 1969)].a < hfs_kappa.STATED_PLAIN[1969]


def test_a_registry_line_keeps_its_own_uncertainty(fitted):
    """A line listed in `inflated_unc_lines.txt` keeps the value entered
    there, and is not set aside by the outlier filter either."""
    lines, sigma, keep, _, _ = fitted
    fixed = hfs_kappa.read_inflated()
    if not fixed:
        pytest.skip('the registry is empty')
    seen = 0
    for i, ln in enumerate(lines):
        want = fixed.lookup(hfs_kappa.line_wn_key(ln))
        if want is None:
            continue
        seen += 1
        assert sigma[i] == pytest.approx(want)
        assert keep[i]
    assert seen


def test_only_a_handful_of_lines_are_set_aside(fitted):
    """The outlier filter must not be doing the fit's work for it."""
    lines, _, keep, _, _ = fitted
    dropped = sum(1 for k in keep if not k)
    assert dropped < 0.01 * len(lines)
