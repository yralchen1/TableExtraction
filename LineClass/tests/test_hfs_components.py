"""Regression tests for the measured hfs components and what is read from them.

The extraction itself is slow (it opens a 3 MB workbook), so the workbook is
read once per session and the cheap checks run against the written file,
`hfs_components.csv`, which is a tracked input of the project: if it ever
stops matching the workbook these tests say so.
"""
import csv
import os

import pytest

import hfs_components
import hfs_patterns

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# the two levels `hfs_line_corrections.csv` still gives S = 0 (section 3.2)
ZEROED_IN_CORRECTIONS = {'059003.000223', '059003.000248'}


@pytest.fixture(scope='module')
def extracted():
    records, report = hfs_components.build()
    return records, report


@pytest.fixture(scope='module')
def written():
    path = os.path.join(HERE, hfs_components.OUTPUT)
    with open(path, encoding='utf-8', newline='') as fh:
        return list(csv.DictReader(fh))


# --- the extraction ---------------------------------------------------------
def test_counts(extracted):
    """500 rows carry Intens = 0; two of them are the curation's own.

    Row 2053 was a third until the curation of 2026-09-23 restored its
    intensity, so nothing is dropped any more.
    """
    _, report = extracted
    assert report['zero_intensity'] == 500
    assert report['sugar_components'] == 498
    assert report['curated_kept'] == 2
    assert report['curated_dropped'] == []


def test_every_component_has_a_flagged_line(extracted):
    _, report = extracted
    assert report['orphans'] == []


def test_the_side_is_right_in_every_case(extracted):
    """The standing test of both the extraction and the reading of the flags.

    A component of an `*r` line must lie at lower wavenumber than the
    tabulated line and one of a `*v` line at higher.  Nothing in the pairing
    enforces it, so a failure here means one of the two is wrong.
    """
    records, report = extracted
    assert report['side_violations'] == []
    assert len(records) == 500


def test_pairing_does_not_depend_on_the_rule(extracted):
    """Imposing the flagged side first gives the same parent for every row."""
    records, _ = extracted
    rows = hfs_components.read_table1()
    flagged = {r['xl_row']: r for r in rows if r['char'] in hfs_components.FLAGS}
    for rec in records:
        wn = float(rec['wn_component'])
        ok = [r for r in flagged.values()
              if abs(r['wn'] - wn) <= hfs_components.PAIR_WINDOW
              and ((r['char'] == '*r' and wn < r['wn'])
                   or (r['char'] == '*v' and wn > r['wn']))]
        nearest = min(ok, key=lambda r: abs(r['wn'] - wn))
        assert nearest['xl_row'] == int(rec['xl_row_line'])


def test_written_file_matches_the_workbook(extracted, written):
    records, _ = extracted
    assert [dict(r) for r in written] == [
        {k: str(v) for k, v in rec.items()} for rec in records]


# --- the file itself --------------------------------------------------------
def test_flag_populations(written):
    """193 components on `*r` lines and 307 on `*v`, over 248 lines."""
    assert sum(1 for r in written if r['char'] == '*r') == 193
    assert sum(1 for r in written if r['char'] == '*v') == 307
    assert len({r['wn_line'] for r in written}) == 248


def test_displacements_have_the_sign_of_the_flag(written):
    for r in written:
        assert (float(r['dwn']) < 0) == (r['char'] == '*r')


def test_n_components_agrees_with_the_rows(written):
    seen = {}
    for r in written:
        seen.setdefault(r['wn_line'], []).append(r)
    for wn, group in seen.items():
        assert all(int(r['n_components']) == len(group) for r in group)


def test_a_level_pair_is_written_only_for_a_single_classification(written):
    for r in written:
        if r['n_class'] == '1':
            assert r['low_id'] and r['upp_id']
        else:
            assert not r['low_id'] and not r['upp_id']


# --- the physics read out of them -------------------------------------------
def test_sixj_against_known_values():
    """One tabulated value, the orthogonality sum, and a triangle-rule zero."""
    assert hfs_patterns.sixj(1, 1, 1, 1, 1, 1) == pytest.approx(1 / 6, abs=1e-12)
    assert hfs_patterns.sixj(1, 1, 5, 1, 1, 1) == 0.0
    # sum_x (2x+1) {a b x; c d p}^2 = 1/(2p+1) when the symbol is diagonal
    total = sum((2 * x + 1) * hfs_patterns.sixj(1, 1, x, 1, 1, 1) ** 2
                for x in (0, 1, 2))
    assert total == pytest.approx(1 / 3, abs=1e-12)


def test_strongest_component_is_the_extreme_one():
    """The head of the ladder, F = I+J1 -> I+J2, is the strongest component.

    The whole reading of the tabulated wavelength rests on this, so it is
    checked over every J pair the line list can present.
    """
    for j1 in (0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5):
        for j2 in (j1 - 1, j1, j1 + 1):
            if j2 < 0.5:
                continue
            comps = hfs_patterns.components(j1, 0.08, j2, 0.02)
            strongest = max(comps, key=lambda c: c[1])
            head = (hfs_patterns.I_SPIN + j1, hfs_patterns.I_SPIN + j2)
            if j1 == j2 == 0.5:
                # the one exception in the whole table: for 1/2 -> 1/2 the
                # F = 3 -> 2 component beats F = 3 -> 3.  No flagged line of
                # this set is a 1/2 -> 1/2 transition, but the extraction
                # would have to be read differently if one appeared.
                assert (strongest[2], strongest[3]) != head
                continue
            assert (strongest[2], strongest[3]) == head


def test_head_displacement_is_the_D_of_the_plan():
    """D = I (A2 J2 - A1 J1) is where the strongest component sits."""
    j1, a1, j2, a2 = 4.5, 0.0942, 3.5, 0.0173
    comps = hfs_patterns.components(j1, a1, j2, a2)
    strongest = max(comps, key=lambda c: c[1])
    assert strongest[0] == pytest.approx(
        hfs_patterns.head_displacement(j1, a1, j2, a2), abs=1e-12)


def test_rung_reproduces_the_component_positions():
    j1, a1, j2, a2 = 6.5, 0.0758, 5.5, 0.0251
    comps = {(c[2], c[3]): c[0] for c in hfs_patterns.components(j1, a1, j2, a2)}
    head = comps[(hfs_patterns.I_SPIN + j1, hfs_patterns.I_SPIN + j2)]
    for k in range(hfs_patterns.ladder_length(j1, j2) + 1):
        f1 = hfs_patterns.I_SPIN + j1 - k
        f2 = hfs_patterns.I_SPIN + j2 - k
        assert comps[(f1, f2)] - head == pytest.approx(
            hfs_patterns.rung(k, j1, a1, j2, a2), abs=1e-12)


def test_the_calculated_constants_reproduce_the_measurements():
    """The Step 2 test: no free parameter, then one scale factor.

    The numbers are the ones `Work_on_hfs_plan.md` section 3.5 quotes, and the
    test is here so that a change to `A_hfs_levels.csv`, to the component file
    or to the classification of a flagged line cannot pass unnoticed.
    """
    patterns, _ = hfs_patterns.read_patterns()
    tested = hfs_patterns.with_both_A(patterns)
    assert len(tested) == 152
    t = hfs_patterns.scale_tests(tested)
    assert t['n'] == 309
    assert t['rms_none'] == pytest.approx(0.087, abs=0.004)
    assert t['scale'] == pytest.approx(0.958, abs=0.010)
    assert t['rms_scale'] == pytest.approx(0.080, abs=0.004)
    # with the scale free, the tabulated line has no measurable offset from
    # the strongest component: that is what says it IS that component.
    assert abs(t['both'][1]) < 3 * t['u_offset']


def test_kappa_one_is_what_the_pipeline_already_applies():
    """`delta_cm-1` of `hfs_line_corrections.csv` is identically D.

    So the measured patterns confirm the correction the pipeline is built to
    apply, rather than calling for a different one.
    """
    path = os.path.join(HERE, 'hfs_line_corrections.csv')
    with open(path, encoding='utf-8', newline='') as fh:
        corrections = {r['wn_obs']: r for r in csv.DictReader(fh)}
    patterns, _ = hfs_patterns.read_patterns()
    checked = stale = 0
    for p in hfs_patterns.with_both_A(patterns):
        row = corrections.get('%.3f' % p.wn)
        if row is None or row['class'] != 'flag':
            continue
        assert float(row['kappa']) == 1.0
        if {row['low_id'], row['upp_id']} & ZEROED_IN_CORRECTIONS:
            # the correction file was written while section 3.2 labeled these
            # two levels as contradicting their own flags, and gives them
            # S = 0.  The review of 2026-09-29 found both conflicts came from
            # rejected lines and a blend, and the labels were withdrawn from
            # `A_hfs_levels.csv`; the file has not been rewritten since.
            assert float(row['S_low' if row['low_id'] in ZEROED_IN_CORRECTIONS
                             else 'S_upp']) == 0.0
            stale += 1
            continue
        assert float(row['delta_cm-1']) == pytest.approx(
            hfs_patterns.head_displacement(p.J1, p.A1, p.J2, p.A2), abs=1e-4)
        checked += 1
    assert (checked, stale) == (145, 7)


def test_no_level_is_labeled_as_conflicting_any_more():
    """Section 3.2's two conflicts were withdrawn on 2026-09-29.

    `059003.000223`: its one contradicting flag, 34926.912 `*v`, was rejected
    on 9/12.  `059003.000248`: 39917.569 `*r` was rejected on 9/12, and the
    `*v` of 33417.350 belongs to the stronger blend component
    `059003.000133`-`059003.000213`.  Every remaining flag of both levels
    agrees with the calculated A.
    """
    with open(os.path.join(HERE, 'A_hfs_levels.csv'), encoding='utf-8',
              newline='') as fh:
        rows = {r['level_id']: r for r in csv.DictReader(fh)}
    assert not [k for k, r in rows.items() if 'CONFLICT' in r['source']]
    for lid in ZEROED_IN_CORRECTIONS:
        assert rows[lid]['source'].startswith('composition')


def test_no_pattern_is_wider_than_its_predicted_span():
    """A listed component cannot lie outside the pattern it belongs to.

    Two of the 152 do, by a few hundredths; both are two-component patterns
    on levels whose A is itself uncertain.  The test holds the count so that a
    worse calculation shows up.
    """
    patterns, _ = hfs_patterns.read_patterns()
    over = 0
    for p in hfs_patterns.with_both_A(patterns):
        xs = [c[0] for c in hfs_patterns.components(p.J1, p.A1, p.J2, p.A2)]
        if max(abs(d) for d in p.dwn) > max(xs) - min(xs):
            over += 1
    assert over <= 2


def test_the_global_fit_is_reported_with_its_singular_direction():
    """Adding a constant to every A is nearly invisible to these data.

    The fit must come out rank-deficient; if it ever comes out full rank the
    interpretation of its uncertainties changes and the text that explains
    them has to change with it.
    """
    patterns, _ = hfs_patterns.read_patterns()
    _, _, _, _, diag = hfs_patterns.fit_levels(patterns)
    assert diag['rank'] < diag['unknowns']
    assert diag['rms'] == pytest.approx(0.031, abs=0.004)


def test_a_uniform_shift_of_every_A_barely_moves_the_prediction():
    """Why the absolute scale is weakly determined, stated as arithmetic."""
    j1, j2, k = 4.5, 4.5, 2
    a = hfs_patterns.rung(k, j1, 0.05, j2, 0.02)
    b = hfs_patterns.rung(k, j1, 0.05 + 0.03, j2, 0.02 + 0.03)
    assert a == pytest.approx(b, abs=1e-12)          # exactly null at dJ = 0
    c = hfs_patterns.rung(k, j1, 0.05, j2 - 1, 0.02)
    d = hfs_patterns.rung(k, j1, 0.05 + 0.03, j2 - 1, 0.02 + 0.03)
    assert abs(d - c) == pytest.approx(k * 0.03, abs=1e-12)
