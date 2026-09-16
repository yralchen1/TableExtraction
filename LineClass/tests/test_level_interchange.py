"""Tests for `level_interchange.py`, the interchanged-identity detector.

Run from the LineClass directory:  python -m pytest tests -q

The tool asks whether two energy levels of the same parity and J are wearing
each other's CALCULATED transition probabilities gA - which would leave their
energies untouched and so escape every other test in the validation.  These
tests exercise the three pieces separately: reading the theoretical labels and
calculated energies out of IDEN2's enlev.dat, the per-configuration energy
window that decides which pairs are even confusable, and the intensity
comparison that decides between the assignment and its interchange.

They use synthetic data throughout, so nothing here depends on the real input
workbooks.
"""
import math
import os
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import level_interchange as li                  # noqa: E402


# Rows of the real enlev.dat.  The first is a level the calculation predicts
# but nobody has observed (uncertainty 5000, no asterisk, E_obs a copy of
# E_calc); the second and third are observed levels, and the third has no
# leading space before its label, which is why the configuration cannot be cut
# at a fixed column.
UNKNOWN = '   1  149995.300  5000.000  149995.300         0.000  1.5 / f5d6d~3D4D/'
KNOWN = '1253      82.000     0.004       0.000 *     -82.000  4.5 / 4f3  ~4I4I/'
TIGHT = '  10  149651.900     0.500  149671.900 *      20.000  5.5 /p5f3d_4I6Ga/'

# A level that HAS been found - it is starred - but whose uncertainty column
# still holds the 5000 placeholder the conversion from Cowan's output writes
# and IDEN2 changes only when the user orders it, level by level.  This is the
# real row 474 as it stood on 2026-09-05.
STARRED_PLACEHOLDER = (
    ' 474  133050.400  5000.000  132958.176 *     -92.224  3.5 / f25f ~3P4D/')


@pytest.fixture
def enlev(tmp_path):
    p = tmp_path / 'enlev.dat'
    p.write_text('\n'.join((UNKNOWN, KNOWN, TIGHT)) + '\n', encoding='latin-1')
    return str(p)


# ---------------------------------------------------------------------------
# The labels
# ---------------------------------------------------------------------------
def test_the_label_splits_at_the_first_separator():
    # "_" is LS coupling, "~" is jj; the configuration is everything before it
    assert li.split_label('f25g _3H4G') == ('f25g', '_3H4G')
    assert li.split_label('f5d6d~3D4D') == ('f5d6d', '~3D4D')
    assert li.split_label('4f3  ~4I4I') == ('4f3', '~4I4I')
    assert li.split_label('p5f3d_4I6Ga') == ('p5f3d', '_4I6Ga')


def test_a_label_without_a_separator_is_all_configuration():
    assert li.split_label(' f25g ') == ('f25g', '')


def test_enlev_fields_and_the_known_flag(enlev):
    en = li.read_enlev(enlev).set_index('idx')
    assert list(en.index) == [1, 1253, 10]

    unknown = en.loc[1]
    assert unknown['E_calc'] == pytest.approx(149995.3)
    assert not unknown['known']          # no star: never observed
    assert unknown['cfg'] == 'f5d6d'

    ground = en.loc[1253]
    assert ground['known']
    assert ground['E_obs'] == pytest.approx(0.0)
    assert ground['E_calc'] == pytest.approx(82.0)
    assert ground['omc'] == pytest.approx(-82.0)   # E_obs - E_calc
    assert ground['J'] == pytest.approx(4.5)
    assert (ground['cfg'], ground['term']) == ('4f3', '~4I4I')


def test_the_star_alone_says_the_level_has_been_found(tmp_path):
    """The uncertainty column is not a second opinion on the star: a level
    IDEN2 has been told about is starred whatever the uncertainty says."""
    p = tmp_path / 'enlev.dat'
    p.write_text(STARRED_PLACEHOLDER + '\n', encoding='latin-1')
    row = li.read_enlev(str(p)).iloc[0]
    assert row['u_obs'] == pytest.approx(5000.0)      # the placeholder
    assert row['known']                               # and yet: found
    assert row['E_obs'] == pytest.approx(132958.176)
    assert (row['cfg'], row['term']) == ('f25f', '~3P4D')


def test_a_found_level_carrying_the_placeholder_is_still_matched(tmp_path):
    """The regression this cost: the level was in the file all along, and the
    run reported it as having no calculated counterpart."""
    p = tmp_path / 'enlev.dat'
    p.write_text('\n'.join((UNKNOWN, STARRED_PLACEHOLDER)) + '\n',
                 encoding='latin-1')
    en = li.read_enlev(str(p))
    levels = pd.DataFrame({'level_id': ['059003.000623'],
                           'E_final': [132958.179]})
    out, unmatched = li.attach_identities(levels, en,
                                          li.configuration_windows(en))
    assert unmatched == []
    assert out.loc[0, 'cfg'] == 'f25f'
    assert out.loc[0, 'E_calc'] == pytest.approx(133050.400)


def test_an_unstarred_level_is_not_matched_however_tight_its_uncertainty():
    """The other direction: E_obs of an unstarred row is a copy of E_calc, a
    position nobody has measured, and must not be matched against."""
    row = '  99  120000.000     0.200  120000.000         0.000  2.5 /f25f _3P4D/'
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.dat', delete=False,
                                     encoding='latin-1') as fh:
        fh.write(KNOWN + '\n' + row + '\n')
        path = fh.name
    en = li.read_enlev(path)
    os.unlink(path)
    assert not bool(en[en['idx'] == 99]['known'].iloc[0])
    levels = pd.DataFrame({'level_id': ['L'], 'E_final': [120000.0]})
    _out, unmatched = li.attach_identities(levels, en,
                                           li.configuration_windows(en))
    assert unmatched == ['L']


def test_the_nearest_row_is_reported_starred_or_not(enlev):
    """What the report prints beside a level it could not match, so that a
    file out of step is told apart from a level the calculation lacks."""
    en = li.read_enlev(enlev)
    row, d = li.nearest_enlev_row(en, 149995.0)
    assert int(row['idx']) == 1 and not row['known']
    assert d == pytest.approx(0.3)
    row, d = li.nearest_enlev_row(en, 0.5)
    assert int(row['idx']) == 1253 and row['known']
    assert d == pytest.approx(0.5)


def test_a_label_flush_against_the_slash_keeps_its_configuration(enlev):
    """The label field is not padded by a fixed amount, so a column cut would
    lose the first character of the longest labels."""
    en = li.read_enlev(enlev).set_index('idx')
    assert en.loc[10, 'cfg'] == 'p5f3d'


# ---------------------------------------------------------------------------
# The energy window
# ---------------------------------------------------------------------------
def test_the_window_is_the_rms_over_the_observed_levels_only():
    en = pd.DataFrame({
        'cfg': ['a', 'a', 'a', 'b'],
        'omc': [3.0, -4.0, 999.0, 10.0],
        # the third level was never observed: its E_obs is a copy of E_calc,
        # so its omc is meaningless and must not enter the rms
        'known': [True, True, False, True],
    })
    w = li.configuration_windows(en)
    assert w['a'] == pytest.approx(math.sqrt((9.0 + 16.0) / 2.0))
    assert w['b'] == pytest.approx(10.0)


def test_the_pair_window_is_the_larger_of_the_two_configurations():
    a = {'W': 20.0}
    b = {'W': 100.0}
    assert li.pair_window(a, b, 1.0) == pytest.approx(100.0)
    assert li.pair_window(a, b, 1.5) == pytest.approx(150.0)


def test_a_level_without_a_label_borrows_its_partner_window():
    a = {'W': float('nan')}
    b = {'W': 100.0}
    assert li.pair_window(a, b, 1.0) == pytest.approx(100.0)
    assert math.isnan(li.pair_window(a, {'W': float('nan')}, 1.0))


def test_attach_identities_matches_by_energy_and_reports_the_misses(enlev):
    en = li.read_enlev(enlev)
    windows = li.configuration_windows(en)
    levels = pd.DataFrame({'level_id': ['L1', 'L2'],
                           'E_final': [0.001, 55555.0]})
    out, unmatched = li.attach_identities(levels, en, windows)
    assert unmatched == ['L2']
    assert out.loc[0, 'cfg'] == '4f3'
    assert out.loc[0, 'E_calc'] == pytest.approx(82.0)
    assert out.loc[1, 'cfg'] == ''
    assert math.isnan(out.loc[1, 'E_calc'])


# ---------------------------------------------------------------------------
# The intensity comparison
# ---------------------------------------------------------------------------
C = 100.0        # the constant of the Boltzmann intensity model
KT = 10000.0     # its temperature parameter, cm^-1
NO_BIAS = {None: (1.0, 1.0)}      # no scale correction, unit spread of ln


def i_of(g, wn, e_up):
    return li.predicted_intensity(g, wn, e_up, C, KT)


def test_the_intensity_model_is_the_pipeline_relation():
    assert i_of(2.0, 5.0e4, 1.0e4) == pytest.approx(
        C * 2.0 * (5.0e4 / 1.0e8) * math.exp(-1.0))


def test_residuals_are_ln_iobs_minus_ln_icalc_with_the_band_weight():
    gA = {('A', 'X'): 7.0}
    lines = [('X', 5.0e4, 3.0, 4.0e4, 1, 0.0)]
    x, y, r, w = li.level_residuals(lines, 'A', gA, 0.0,
                                    {None: (1.0, 2.0)}, C, KT)
    i_calc = i_of(7.0, 5.0e4, 4.0e4)
    assert x[0] == pytest.approx(math.log(i_calc))
    assert y[0] == pytest.approx(math.log(3.0))
    assert r[0] == pytest.approx(math.log(3.0 / i_calc))
    assert w[0] == pytest.approx(0.25)          # 1 / sd^2, sd = 2


def test_the_band_scale_factor_multiplies_the_prediction():
    gA = {('A', 'X'): 7.0}
    lines = [('X', 5.0e4, 3.0, 4.0e4, 1, 0.0)]
    plain = li.level_residuals(lines, 'A', gA, 0.0, NO_BIAS, C, KT)[2][0]
    tenth = li.level_residuals(lines, 'A', gA, 0.0, {None: (0.1, 1.0)},
                               C, KT)[2][0]
    assert tenth - plain == pytest.approx(math.log(10.0))


# ---------------------------------------------------------------------------
# Blended lines
# ---------------------------------------------------------------------------
def test_an_unblended_line_keeps_the_whole_measured_intensity():
    assert li.blend_share(5.0, 1, None) == 1.0
    assert li.blend_share(5.0, 1, 95.0) == 1.0


def test_a_blend_is_split_by_the_calculated_intensities():
    """Two transitions on one measured line, the other calculated three times
    as strong: this one is given a quarter of what was measured."""
    assert li.blend_share(10.0, 2, 30.0) == pytest.approx(0.25)


def test_a_blend_with_an_uncalculated_component_is_split_evenly():
    """`sib` is None when one of the other components has no calculated
    intensity: there is nothing to compare with, so each takes 1/n."""
    assert li.blend_share(10.0, 4, None) == pytest.approx(0.25)
    assert li.blend_share(10.0, 2, 0.0) == pytest.approx(1.0)


def test_the_split_follows_the_hypothesis_being_scored():
    """The same measured line under two identities.  The OTHER component of the
    blend is not touched by the swap, so its calculated intensity stays put; the
    component being scored takes a tenth of the feature under the identity that
    makes it nine times fainter than its neighbour, and half under the identity
    that makes it equal - the share moves with the hypothesis."""
    gA = {('A', 'X'): 1.0, ('B', 'X'): 9.0}
    sib = i_of(9.0, 5.0e4, 4.0e4)           # the other component, unchanged
    lines = [('X', 5.0e4, 100.0, 4.0e4, 2, sib)]
    y_a = li.level_residuals(lines, 'A', gA, 0.0, NO_BIAS, C, KT)[1][0]
    y_b = li.level_residuals(lines, 'B', gA, 0.0, NO_BIAS, C, KT)[1][0]
    assert y_a == pytest.approx(math.log(100.0 * 0.1))
    assert y_b == pytest.approx(math.log(100.0 * 0.5))


def test_a_blend_component_is_no_longer_given_the_whole_feature():
    """The point of the split: an unblended line and a blend component of the
    same measured intensity must not enter with the same observed value."""
    gA = {('A', 'X'): 1.0}
    i_calc = i_of(1.0, 5.0e4, 4.0e4)
    alone = li.level_residuals([('X', 5.0e4, 100.0, 4.0e4, 1, 0.0)],
                               'A', gA, 0.0, NO_BIAS, C, KT)[1][0]
    shared = li.level_residuals([('X', 5.0e4, 100.0, 4.0e4, 2, i_calc)],
                                'A', gA, 0.0, NO_BIAS, C, KT)[1][0]
    assert alone - shared == pytest.approx(math.log(2.0))


def test_accepted_lines_by_level_reads_the_blend_off_the_table():
    """Two accepted rows at one wavenumber are one blend; a row at another
    wavenumber is alone.  Both rows of the blend get the sibling's calculated
    intensity, and the shares they imply are the BF column of the pipeline."""
    real = pd.DataFrame({
        'accepted': [1, 1, 1],
        'low_id': ['L1', 'L2', 'L3'],
        'upp_id': ['U1', 'U1', 'U2'],
        'wn_obs': ['1000.0', '1000.0', '2000.0'],
        'obs_intens': [100.0, 100.0, 50.0],
        'calc_intens': [30.0, 10.0, 7.0],
    })
    e_final = {k: 0.0 for k in ('L1', 'L2', 'L3')}
    e_final.update({'U1': 1000.0, 'U2': 2000.0})
    by = li.accepted_lines_by_level(real, e_final)
    (_p, _wn, _i, _e, n1, sib1), = by['L1']
    (_p, _wn, _i, _e, n2, sib2), = by['L2']
    (_p, _wn, _i, _e, n3, sib3), = by['L3']
    assert (n1, sib1) == (2, 10.0)      # L1's sibling carries 10
    assert (n2, sib2) == (2, 30.0)      # L2's carries 30
    assert (n3, sib3) == (1, 0.0)       # nothing shares that line
    assert li.blend_share(30.0, n1, sib1) == pytest.approx(0.75)
    assert li.blend_share(10.0, n2, sib2) == pytest.approx(0.25)


def test_a_missing_calculated_intensity_makes_the_whole_blend_even():
    real = pd.DataFrame({
        'accepted': [1, 1],
        'low_id': ['L1', 'L2'],
        'upp_id': ['U1', 'U1'],
        'wn_obs': ['1000.0', '1000.0'],
        'obs_intens': [100.0, 100.0],
        'calc_intens': [30.0, np.nan],
    })
    e_final = {'L1': 0.0, 'L2': 0.0, 'U1': 1000.0}
    by = li.accepted_lines_by_level(real, e_final)
    assert all(t[4] == 2 and t[5] is None
               for lid in ('L1', 'L2') for t in by[lid])


def test_the_slope_of_a_perfect_prediction_is_one():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    assert li.slope(x, 2.0 + x) == pytest.approx(1.0)
    assert li.slope(x, np.zeros(4)) == pytest.approx(0.0)


def test_a_slope_needs_three_points_and_some_spread():
    assert math.isnan(li.slope(np.array([1.0, 2.0]), np.array([1.0, 2.0])))
    assert math.isnan(li.slope(np.ones(4), np.arange(4.0)))


def test_common_lines_drops_what_one_hypothesis_cannot_score():
    # under policy "none" (g_imp = 0) a pair absent from the calculated file
    # has no predicted intensity at all, so the line is left out of BOTH sides
    gA = {('A', 'X'): 1.0, ('B', 'X'): 1.0, ('A', 'Y'): 1.0}
    lines = [('X', 1.0e4, 1.0, 1.0e4, 1, 0.0),
             ('Y', 1.0e4, 1.0, 1.0e4, 1, 0.0)]
    assert [t[0] for t in li.common_lines(lines, 'A', 'B', gA, 0.0)] == ['X']
    # with imputation both sides can score Y, so nothing is dropped
    assert [t[0] for t in li.common_lines(lines, 'A', 'B', gA, 1.0)] == \
        ['X', 'Y']


def test_common_lines_drops_a_partner_that_is_the_other_level_of_the_pair():
    """Such a line would have its own gA exchanged by the same swap, so it
    cannot be used to judge it (and cannot exist: the two share a parity)."""
    gA = {('A', 'B'): 1.0, ('A', 'X'): 1.0, ('B', 'X'): 1.0}
    lines = [('B', 1.0e4, 1.0, 1.0e4, 1, 0.0),
             ('X', 1.0e4, 1.0, 1.0e4, 1, 0.0)]
    assert [t[0] for t in li.common_lines(lines, 'A', 'B', gA, 1.0)] == ['X']


# --- a constructed interchange ---------------------------------------------
#
# Two levels A and B, four common partners.  The observed intensities of A
# follow one set of gA values and those of B follow the other, so the data are
# a clean interchange: the SWAPPED hypothesis reproduces them exactly and the
# assignment does not.
PARTNERS = ['P1', 'P2', 'P3', 'P4']
WN = 5.0e4
E_UP = 4.0e4
GA_A = dict(zip(PARTNERS, [1.0e3, 1.0e4, 1.0e5, 1.0e6]))
GA_B = dict(zip(PARTNERS, [1.0e6, 1.0e5, 1.0e4, 1.0e3]))
GA = {tuple(sorted(('A', p))): GA_A[p] for p in PARTNERS}
GA.update({tuple(sorted(('B', p))): GA_B[p] for p in PARTNERS})


def lines_following(ga_of_partner):
    return [(p, WN, i_of(ga_of_partner[p], WN, E_UP), E_UP, 1, 0.0)
            for p in PARTNERS]


def test_an_interchange_is_preferred_and_shows_it_in_every_column():
    by_level = {'A': lines_following(GA_B),    # A emits B's pattern
                'B': lines_following(GA_A)}    # and B emits A's
    ev = li.pair_evidence('A', 'B', by_level, GA, 0.0, NO_BIAS, C, KT)
    assert ev['lnR'] < 0                       # the swap fits better
    assert ev['lnR_A'] < 0 and ev['lnR_B'] < 0   # and does so on both levels
    assert ev['rms_swapped'] == pytest.approx(0.0, abs=1e-9)
    assert ev['rms_assigned'] > 1.0
    assert ev['slope_A_swapped'] == pytest.approx(1.0)
    assert ev['slope_B_swapped'] == pytest.approx(1.0)
    assert li.slopes_prefer_swap(ev)


def test_a_correct_assignment_is_preferred():
    by_level = {'A': lines_following(GA_A),
                'B': lines_following(GA_B)}
    ev = li.pair_evidence('A', 'B', by_level, GA, 0.0, NO_BIAS, C, KT)
    assert ev['lnR'] > 0
    assert ev['rms_assigned'] == pytest.approx(0.0, abs=1e-9)
    assert ev['slope_A_assigned'] == pytest.approx(1.0)
    assert not li.slopes_prefer_swap(ev)


def test_lnR_splits_into_the_two_levels_contributions():
    by_level = {'A': lines_following(GA_B), 'B': lines_following(GA_B)}
    ev = li.pair_evidence('A', 'B', by_level, GA, 0.0, NO_BIAS, C, KT)
    assert ev['lnR'] == pytest.approx(ev['lnR_A'] + ev['lnR_B'])
    # only A's pattern is the other level's, so only A asks for the swap
    assert ev['lnR_A'] < 0 < ev['lnR_B']


def test_a_pair_with_no_scorable_line_returns_nothing():
    assert li.pair_evidence('A', 'B', {'A': [], 'B': []}, GA, 0.0,
                            NO_BIAS, C, KT) is None


def test_slopes_prefer_swap_needs_both_levels_and_a_defined_slope():
    both = {'slope_A_assigned': 0.2, 'slope_A_swapped': 0.9,
            'slope_B_assigned': 0.3, 'slope_B_swapped': 0.8}
    assert li.slopes_prefer_swap(both)
    one = dict(both, slope_B_swapped=0.1)      # B gets worse
    assert not li.slopes_prefer_swap(one)
    unknown = dict(both, slope_A_swapped=float('nan'))
    assert not li.slopes_prefer_swap(unknown)


def test_mark_rows_flags_only_on_all_three_conditions():
    good = {'slope_A_assigned': 0.2, 'slope_A_swapped': 0.9,
            'slope_B_assigned': 0.3, 'slope_B_swapped': 0.8}
    null = np.array([10.0, 20.0, 30.0, 40.0])   # four false swaps
    rows = li.mark_rows([
        dict(good, lnR=-5.0),        # below every false swap, slopes agree
        dict(good, lnR=+25.0),       # the swap does not even fit better
        dict(good, lnR=-5.0, slope_B_swapped=0.1),   # B gets worse
    ], null, p_flag=0.05)
    assert [r['flagged'] for r in rows] == [1, 0, 0]
    assert rows[0]['p_null'] == pytest.approx(0.0)
    assert rows[1]['p_null'] == pytest.approx(0.5)


def test_without_a_calibration_sample_p_null_is_unknown():
    row = li.mark_rows([{'lnR': -5.0, 'slope_A_assigned': 0.2,
                         'slope_A_swapped': 0.9, 'slope_B_assigned': 0.3,
                         'slope_B_swapped': 0.8}], np.array([]), 0.05)[0]
    assert math.isnan(row['p_null'])
    assert row['flagged'] == 1          # decided by lnR and the slopes alone


# ---------------------------------------------------------------------------
# The search
# ---------------------------------------------------------------------------
def levels_frame():
    """Four levels: two close pairs of one J and parity, plus a level of
    another J that must never be paired with them."""
    return pd.DataFrame({
        'level_id': ['A', 'B', 'C', 'D'],
        'J': ['5/2', '5/2', '5/2', '7/2'],
        'parity': ['o', 'o', 'o', 'o'],
        'E_final': [1000.0, 1050.0, 1400.0, 1010.0],
        'n_tot': [5, 5, 5, 5],
        'cfg': ['x', 'x', 'x', 'x'],
        'W': [100.0, 100.0, 100.0, 100.0],
    })


def test_the_search_separates_candidates_from_the_calibration_sample():
    cand, null = li.search_pairs(levels_frame(), min_lines=3, window=1.0,
                                 null_window=10.0)
    got = {tuple(sorted((a['level_id'], b['level_id']))) for a, b, _w in cand}
    assert got == {('A', 'B')}                 # 50 cm^-1 apart, window 100
    null_got = {tuple(sorted((a['level_id'], b['level_id'])))
                for a, b, _w in null}
    assert null_got == {('A', 'C'), ('B', 'C')}   # 400 and 350 cm^-1 apart
    # D has another J: it appears in neither list
    assert all('D' not in p for p in got | null_got)


def test_a_level_with_too_few_lines_is_not_searched():
    lv = levels_frame()
    lv.loc[1, 'n_tot'] = 2
    cand, _null = li.search_pairs(lv, min_lines=3, window=1.0,
                                  null_window=10.0)
    assert cand == []


def test_the_calibration_sample_stops_at_the_outer_edge():
    lv = levels_frame()
    lv.loc[2, 'E_final'] = 9000.0             # 8000 cm^-1 = 80 x the rms
    _cand, null = li.search_pairs(lv, min_lines=3, window=1.0,
                                  null_window=10.0)
    assert null == []


def test_a_level_in_the_id_table_is_tied_to_its_row_not_by_energy():
    """CLAUDE.md: join through IDEN_level_ids.txt, never by energy.

    The level has just been moved in the run, so its energy is nowhere near
    the E_obs enlev.dat still holds, and the nearest starred row is a different
    level.  The table's row is the one taken.
    """
    rows = ['  10  120000.000     0.200  120050.000 *      50.000  2.5 /f25f _3P4D/',
            '  11  121000.000     0.200  121000.000 *       0.000  2.5 /f27p ~3H4H/']
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.dat', delete=False,
                                     encoding='latin-1') as fh:
        fh.write('\n'.join(rows) + '\n')
        path = fh.name
    en = li.read_enlev(path)
    os.unlink(path)
    levels = pd.DataFrame({'level_id': ['059003.000900'],
                           'E_final': [120990.0]})
    win = li.configuration_windows(en)
    by_energy, _ = li.attach_identities(levels, en, win, tol=50.0)
    assert by_energy.loc[0, 'cfg'] == 'f27p'
    out, unmatched = li.attach_identities(levels, en, win, tol=50.0,
                                          row_of_id={'059003.000900': 10})
    assert unmatched == []
    assert (out.loc[0, 'cfg'], out.loc[0, 'term']) == ('f25f', '_3P4D')
    assert out.loc[0, 'E_calc'] == pytest.approx(120000.0)
    assert out.loc[0, 'omc'] == pytest.approx(990.0)
