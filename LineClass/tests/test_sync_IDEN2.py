"""Tests for `cowan_gA` and `sync_IDEN2`.

Run from the LineClass directory:  python -m pytest tests -q

Two groups.  Most tests build a miniature IDEN2 - four levels, a handful of
transitions - in a temporary directory and check that the rewrite does what it
says: the fixed-width format survives, an identified line is never dropped or
disturbed, the predicted wavenumber is the difference of the two energies
exactly as they are written in the file, and the intensity code is
round(10*ln(Icalc)).  A few tests, marked `real_file`, open the real IDEN2
directory and the real transition-probability workbook and check the numbers
the analysis will actually use; they are skipped when those files are absent,
so the first group still runs on a machine that has only the code.
"""
import math
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import cowan_gA                            # noqa: E402
import sync_IDEN2 as sync                  # noqa: E402
import swap_line_assignments_IDEN as IDEN  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REAL_IDEN2 = os.path.join(HERE, 'IDEN2')
REAL_TP = cowan_gA.TP_FILE

C = 269.979733
KT = 12233.068631

real_file = pytest.mark.skipif(
    not (os.path.exists(os.path.join(REAL_IDEN2, 'trans.dat'))
         and os.path.exists(REAL_TP)),
    reason='the real IDEN2 directory or the transition workbook is absent')


# ---------------------------------------------------------------------------
# A miniature IDEN2
# ---------------------------------------------------------------------------
# Four levels.  1 and 2 are high and not found; 3 and 4 are low and found.
# enlev.dat is in order of decreasing calculated energy, which is the order
# IDEN2 keeps it in and the order the level numbers follow.
MINI_LEVELS = [
    # index, E_calc, unc, E_obs, known, J, label
    (1, 130000.0, 500.0, 130000.0, False, 4.5, ' fd6p _3H4H'),
    (2, 120000.0, 500.0, 120000.0, False, 3.5, ' fd6p _3H4G'),
    (3, 10000.0, 0.005, 10000.500, True, 5.5, ' f25d _3H4K'),
    (4, 100.0, 0.004, 120.250, True, 4.5, '4f3  ~4I4I '),
]
# (lid, index) - the calculation's own numbering, deliberately not IDEN2's
MINI_LIDS = {104: 1, 103: 2, 102: 3, 101: 4}


def enlev_record(index, e_calc, unc, e_obs, known, J, label):
    return ('%4d%12.3f%10.3f%12.3f%s%12.3f%5s /%s/'
            % (index, e_calc, unc, e_obs,
               IDEN.STAR_ON if known else IDEN.STAR_OFF, e_obs - e_calc,
               '%.1f' % J, label))


def trans_header(index, J, label, e_obs, known):
    return ('$%4d    J=%4s  %-11s%13.3f%s '
            % (index, '%.1f' % J, label, e_obs,
               IDEN.STAR_ON if known else IDEN.STAR_OFF))


def mini_iden2(tmp_path, rows):
    """Write a miniature enlev.dat and trans.dat.  `rows` is a list of
    (owner, partner, code, assignment tail)."""
    d = tmp_path / 'IDEN2'
    d.mkdir(parents=True)
    enlev = [enlev_record(*lev) for lev in MINI_LEVELS]
    (d / 'enlev.dat').write_bytes(
        ''.join(r + '\r\n' for r in enlev).encode('latin-1'))

    e_obs = {lev[0]: lev[3] for lev in MINI_LEVELS}
    known = {lev[0]: lev[4] for lev in MINI_LEVELS}
    text = []
    for owner in sorted({r[0] for r in rows}):
        lev = [x for x in MINI_LEVELS if x[0] == owner][0]
        text.append(trans_header(owner, lev[5], lev[6], e_obs[owner],
                                 known[owner]))
        for o, partner, code, tail in sorted(rows):
            if o != owner:
                continue
            text.append(sync.transition_row(
                partner, code, e_obs[partner], known[partner],
                e_obs[owner] - e_obs[partner], tail))
    (d / 'trans.dat').write_bytes(
        ''.join(r + '\r\n' for r in text).encode('latin-1'))

    with open(d / 'IDEN_level_ids.txt', 'w', encoding='utf-8',
              newline='\n') as fh:
        fh.write('level_id\tIDEN_id\n')
        fh.write('059003.000003\t3\n')
        fh.write('059003.000004\t4\n')
    return d


def mini_lopt(tmp_path, e3=10000.500, e4=120.250):
    p = tmp_path / 'LOPT_output_levels.txt'
    with open(p, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('Designation\tEnergy\tD1\tD2stat\tD2sys\tD2tot\tN_lines\t'
                 'Comments\n')
        fh.write('059003.000003\t%.4f\t0.0050\t0.0\t0.0\t0.0030\t5\t\n' % e3)
        fh.write('059003.000004\t%.4f\t0.0040\t0.0\t0.0\t0.0020\t7\t\n' % e4)
    return p


def mini_tp(tmp_path, pairs):
    """A stand-in for the transition workbook.  `pairs` is
    (lid_low, lid_up, gA)."""
    rows = []
    e_of = {lid: [x for x in MINI_LEVELS if x[0] == idx][0][1]
            for lid, idx in MINI_LIDS.items()}
    j_of = {lid: [x for x in MINI_LEVELS if x[0] == idx][0][5]
            for lid, idx in MINI_LIDS.items()}
    for low, up, gA in pairs:
        rows.append(dict(
            lid1=low, conf1='4f3.5p6', term1='4I', parity1='odd',
            J1=j_of[low],
            lid2=up, conf2='4f.5p6.5d.(3H*).6p', term2='4H', parity2='even',
            J2=j_of[up],
            E_low=e_of[low], E_up=e_of[up],
            E_low_exp=np.nan, E_up_exp=np.nan, id1='', id2='',
            A=gA / (2.0 * j_of[up] + 1.0), u_gA_pct=50.0, gA=gA))
    columns = ['lid1', 'conf1', 'term1', 'parity1', 'J1',
               'lid2', 'conf2', 'term2', 'parity2', 'J2',
               'E_low', 'E_up', 'E_low_exp', 'E_up_exp', 'id1', 'id2',
               'A', 'u_gA_pct', 'gA']
    return pd.DataFrame(rows, columns=columns)


def run_mini(tmp_path, rows, pairs, **kw):
    """Sync a miniature IDEN2 and return its rewritten files."""
    d = mini_iden2(tmp_path, rows)
    lopt = mini_lopt(tmp_path, **kw.pop('lopt', {}))
    tp = mini_tp(tmp_path, pairs)

    enlev = IDEN.Enlev(str(d / 'enlev.dat'))
    trans = IDEN.Trans(str(d / 'trans.dat'))
    id_of_row = cowan_gA.read_id_map(str(d / 'IDEN_level_ids.txt'))
    energies, changes = sync.adopted_energies(
        enlev, id_of_row, sync.read_lopt_levels(str(lopt)), lambda *a: None)
    pred = sync.predicted_transitions(tp, MINI_LIDS, energies, C, KT)
    records, report = sync.rebuild_trans(trans, enlev, energies, pred,
                                         kw.get('cutoff', -41),
                                         lambda *a: None)
    for index, (E, unc, is_known) in energies.items():
        enlev.set_measurement(index, unc, E, is_known)
    return enlev, records, report, energies, changes


# ---------------------------------------------------------------------------
# The reader
# ---------------------------------------------------------------------------
def test_level_id_is_the_string_the_pipeline_uses():
    assert cowan_gA.level_id(59003.000042) == '059003.000042'
    assert cowan_gA.level_id(np.nan) == ''
    assert cowan_gA.level_id(None) == ''


def test_alignment_survives_an_order_reversal():
    """Two neighbouring levels that come out in opposite order in the two
    files must still be paired with each other, not left unmatched."""
    left = pd.DataFrame({'idx': [1, 2, 3], 'E_calc': [100.0, 200.5, 201.0],
                         'E_obs': [0.0, 0.0, 0.0], 'J': [0.5, 1.5, 2.5],
                         'label': ['a', 'b', 'c']})
    right = pd.DataFrame({'lid': [11, 12, 13],
                          'conf': ['x'] * 3, 'term': ['y'] * 3,
                          'parity': ['odd'] * 3, 'J': [0.5, 2.5, 1.5],
                          'E_calc': [100.1, 200.4, 201.2],
                          'E_exp': [np.nan] * 3, 'level_id': [''] * 3})
    mapping, report = cowan_gA.match_to_enlev(right, left)
    assert mapping == {11: 1, 12: 3, 13: 2}
    # The two reversed levels cross, so a no-crossing alignment can pair at
    # most one of them and the second pass has to pick the other up.
    assert report['n_reordered'] >= 1
    assert not report['unmatched_lid'] and not report['unmatched_idx']


def test_alignment_never_pairs_different_J():
    left = pd.DataFrame({'idx': [1], 'E_calc': [100.0], 'E_obs': [0.0],
                         'J': [0.5], 'label': ['a']})
    right = pd.DataFrame({'lid': [11], 'conf': ['x'], 'term': ['y'],
                          'parity': ['odd'], 'J': [4.5], 'E_calc': [100.0],
                          'E_exp': [np.nan], 'level_id': ['']})
    mapping, report = cowan_gA.match_to_enlev(right, left)
    assert mapping == {}
    assert report['unmatched_lid'] == [11]
    assert report['unmatched_idx'] == [1]


def test_a_contradicted_level_id_is_an_error():
    """The 594 levels whose IDEN2 row and experimental identifier are both
    known are the check on the whole alignment; one of them coming out wrong
    must stop the run rather than return a mapping."""
    left = pd.DataFrame({'idx': [1], 'E_calc': [100.0], 'E_obs': [0.0],
                         'J': [0.5], 'label': ['a']})
    right = pd.DataFrame({'lid': [11], 'conf': ['x'], 'term': ['y'],
                          'parity': ['odd'], 'J': [0.5], 'E_calc': [100.0],
                          'E_exp': [np.nan], 'level_id': ['059003.000999']})
    with pytest.raises(cowan_gA.MatchError):
        cowan_gA.match_to_enlev(right, left, {1: '059003.000001'})


# ---------------------------------------------------------------------------
# The uncertainty rule and the energies
# ---------------------------------------------------------------------------
def test_uncertainty_is_the_larger_of_D1_and_D2tot(tmp_path):
    lopt = sync.read_lopt_levels(str(mini_lopt(tmp_path)))
    assert lopt['059003.000003'] == (10000.5, 0.005)   # D1 = 0.005 wins
    assert lopt['059003.000004'] == (120.25, 0.004)    # D1 = 0.004 wins


def test_energies_are_rounded_to_the_three_decimals_the_file_holds(tmp_path):
    """IDEN2 reads a wavenumber back as the difference of the two energies as
    they are written.  Rounding has to happen before the subtraction, not
    after, or the file contradicts itself in the third decimal."""
    rows = [(2, 3, 10, IDEN.BLANK_OBS)]
    _, records, _, energies, _ = run_mini(
        tmp_path, rows, [(102, 103, 1.0e6)],
        lopt=dict(e3=10000.5006, e4=120.2504))
    assert energies[3][0] == 10000.501
    assert energies[4][0] == 120.250


# ---------------------------------------------------------------------------
# The rewrite
# ---------------------------------------------------------------------------
def test_the_rewrite_keeps_the_fixed_width_format(tmp_path):
    rows = [(1, 3, 20, IDEN.BLANK_OBS), (1, 4, 15, IDEN.BLANK_OBS),
            (2, 4, 5, IDEN.BLANK_OBS)]
    pairs = [(102, 104, 1.0e7), (101, 104, 5.0e6), (101, 103, 2.0e6)]
    _, records, _, _, _ = run_mini(tmp_path, rows, pairs)
    for rec in records:
        assert len(rec) == (sync.HEADER_WIDTH if rec.startswith('$')
                            else sync.ROW_WIDTH)


def test_the_wavenumber_is_the_difference_of_the_two_written_energies(
        tmp_path):
    rows = [(1, 3, 20, IDEN.BLANK_OBS), (1, 4, 15, IDEN.BLANK_OBS)]
    pairs = [(102, 104, 1.0e7), (101, 104, 5.0e6)]
    enlev, records, _, _, _ = run_mini(tmp_path, rows, pairs)
    owner = None
    seen = 0
    for rec in records:
        if rec.startswith('$'):
            owner = int(IDEN.field(rec, IDEN.TH_INDEX))
            assert IDEN.number(rec, IDEN.TH_EOBS) == enlev.e_obs(owner)
        else:
            partner = int(IDEN.field(rec, IDEN.TR_PARTNER))
            assert IDEN.number(rec, IDEN.TR_EPART) == enlev.e_obs(partner)
            assert (IDEN.number(rec, IDEN.TR_WN)
                    == pytest.approx(enlev.e_obs(owner)
                                     - enlev.e_obs(partner), abs=5e-4))
            seen += 1
    assert seen == 2


def test_the_intensity_code_is_ten_times_the_logarithm(tmp_path):
    gA = 1.0e7
    rows = [(1, 4, 0, IDEN.BLANK_OBS)]
    _, records, _, energies, _ = run_mini(tmp_path, rows, [(101, 104, gA)])
    row = [r for r in records if r.startswith('+')][0]
    rwn = energies[1][0] - energies[4][0]
    Icalc = C * gA * (rwn / 1.0e8) * math.exp(-energies[1][0] / KT)
    assert int(IDEN.field(row, sync.TR_ICALC)) == round(10.0 * math.log(Icalc))


def test_an_identified_line_is_never_dropped(tmp_path):
    """A transition that carries an observed line stays, however weak the new
    intensity code makes it - the identification is the analyst's work and the
    cutoff is a convenience."""
    tail = '   104   12000.500      0.000  6426'
    rows = [(1, 4, 20, tail)]
    _, records, report, _, _ = run_mini(tmp_path, rows, [(101, 104, 1.0e-30)],
                                        cutoff=0)
    kept = [r for r in records if r.startswith('+')]
    assert len(kept) == 1
    assert int(IDEN.field(kept[0], sync.TR_ICALC)) < 0
    assert IDEN.obs_wavenumber(IDEN.assignment(kept[0])) == 12000.5
    assert IDEN.obs_row(IDEN.assignment(kept[0])) == 6426
    assert not report['dropped']


def test_a_transition_with_no_gA_keeps_its_row_and_its_code(tmp_path):
    """The hand-added row for Sugar's identification of a line below the gA
    floor of the Cowan run has no calculated strength to recompute from."""
    tail = '   104   12000.500      0.000  6426'
    rows = [(1, 4, 0, tail)]
    _, records, report, energies, _ = run_mini(tmp_path, rows, [])
    kept = [r for r in records if r.startswith('+')]
    assert len(kept) == 1
    assert int(IDEN.field(kept[0], sync.TR_ICALC)) == 0
    assert report['kept_no_gA'] == [(1, 4)]
    # the geometry is still brought up to date
    assert (IDEN.number(kept[0], IDEN.TR_WN)
            == pytest.approx(energies[1][0] - energies[4][0], abs=5e-4))


def test_a_transition_with_no_gA_and_no_line_is_dropped(tmp_path):
    rows = [(1, 4, -20, IDEN.BLANK_OBS)]
    _, records, report, _, _ = run_mini(tmp_path, rows, [])
    assert not [r for r in records if r.startswith('+')]
    assert len(report['dropped']) == 1
    assert report['dropped'][0][2] != 'below the cutoff'


def test_the_observed_minus_predicted_follows_the_new_prediction(tmp_path):
    """The level moves; the identified line does not.  The column that says
    how far the line is from the prediction has to move with the level."""
    tail = '   104   12000.500      0.000  6426'
    rows = [(1, 4, 20, tail)]
    _, records, _, energies, _ = run_mini(
        tmp_path, rows, [(101, 104, 1.0e7)], lopt=dict(e4=130.250))
    row = [r for r in records if r.startswith('+')][0]
    predicted = energies[1][0] - energies[4][0]
    assert (IDEN.obs_omc(IDEN.assignment(row))
            == pytest.approx(12000.5 - predicted, abs=5e-4))


def test_a_transition_above_the_cutoff_is_added(tmp_path):
    rows = [(1, 4, 20, IDEN.BLANK_OBS)]
    pairs = [(101, 104, 1.0e7), (101, 103, 1.0e7)]
    _, records, report, _, _ = run_mini(tmp_path, rows, pairs)
    assert len(report['added']) == 1
    assert len([r for r in records if r.startswith('+')]) == 2


def test_the_cutoff_is_obeyed(tmp_path):
    rows = [(1, 4, 20, IDEN.BLANK_OBS)]
    pairs = [(101, 104, 1.0e7)]
    for cutoff, expect in ((-1000, 1), (1000, 0)):
        _, records, _, _, _ = run_mini(tmp_path / str(cutoff), rows, pairs,
                                       cutoff=cutoff)
        assert len([r for r in records if r.startswith('+')]) == expect


def test_rows_are_in_the_order_IDEN2_keeps_them(tmp_path):
    """Blocks in order of the upper level's number, rows within a block in
    order of the partner's number - the order the file already has."""
    rows = [(1, 3, 20, IDEN.BLANK_OBS), (1, 4, 15, IDEN.BLANK_OBS),
            (2, 4, 10, IDEN.BLANK_OBS)]
    pairs = [(102, 104, 1.0e7), (101, 104, 1.0e7), (101, 103, 1.0e7)]
    _, records, _, _, _ = run_mini(tmp_path, rows, pairs)
    owners, order = [], []
    for rec in records:
        if rec.startswith('$'):
            owners.append(int(IDEN.field(rec, IDEN.TH_INDEX)))
            order.append([])
        else:
            order[-1].append(int(IDEN.field(rec, IDEN.TR_PARTNER)))
    assert owners == sorted(owners)
    assert all(part == sorted(part) for part in order)


def test_a_level_that_is_not_found_is_left_alone(tmp_path):
    rows = [(1, 4, 20, IDEN.BLANK_OBS)]
    enlev, _, _, energies, _ = run_mini(tmp_path, rows, [(101, 104, 1.0e7)])
    assert energies[1] == (130000.0, 500.0, False)
    assert not enlev.known(1)
    assert enlev.e_obs(1) == 130000.0


def _enlev_of(tmp_path, levels):
    """An `Enlev` over an arbitrary level list, for the uncertainty rules.

    `levels` is a list of (index, E_calc, unc, E_obs, known, J, label); the
    file has to be written and read back because the uncertainty rules work
    from the label, and the label is a fixed-width field of the file.
    """
    p = tmp_path / 'enlev.dat'
    p.write_bytes(''.join(enlev_record(*lev) + '\r\n'
                          for lev in levels).encode('latin-1'))
    return IDEN.Enlev(str(p))


def test_an_unfound_level_takes_the_rms_of_its_configuration(tmp_path):
    """A level nobody has found carries no measurement, so its uncertainty
    says instead how far from its calculated position it is likely to be:
    the rms of E_obs - E_calc over the levels of the same configuration that
    HAVE been found.  Here that is sqrt((30^2 + 40^2 + 50^2)/3) = 40.825."""
    levels = [(1, 5000.0, 5000.0, 5000.0, False, 4.5, ' f25g _3H4H'),
              (2, 1000.0, 0.01, 1030.0, True, 4.5, ' f25g _3H4G'),
              (3, 2000.0, 0.01, 2040.0, True, 3.5, ' f25g _3H4F'),
              (4, 3000.0, 0.01, 2950.0, True, 2.5, ' f25g _3H4D')]
    enlev = _enlev_of(tmp_path, levels)
    energies = {n: (lev[3], lev[2], lev[4])
                for n, lev in zip(range(1, 5), levels)}
    energies, rep = sync.configuration_uncertainties(enlev, energies,
                                                     lambda *a: None)
    expect = round(math.sqrt((30.0 ** 2 + 40.0 ** 2 + 50.0 ** 2) / 3.0), 3)
    assert energies[1] == (5000.0, expect, False)
    assert rep['changes'] == [(1, 'f25g', 5000.0, expect)]
    # The found levels keep the uncertainty of their own measurement.
    assert energies[2] == (1030.0, 0.01, True)


def test_a_configuration_with_no_found_level_is_left_alone(tmp_path):
    """Six configurations of Pr III have not one level found in them.  There
    is nothing to take an rms of, and the rms of the whole list would be a
    claim about them that nothing supports, so those levels keep whatever
    enlev.dat already holds - the placeholder 5000 for a level nobody has
    touched."""
    levels = [(1, 5000.0, 5000.0, 5000.0, False, 4.5, ' f5d6d_3H4H'),
              (2, 1000.0, 0.01, 1030.0, True, 4.5, ' f25g _3H4G'),
              (3, 2000.0, 0.01, 2040.0, True, 3.5, ' f25g _3H4F'),
              (4, 3000.0, 0.01, 2950.0, True, 2.5, ' f25g _3H4D')]
    enlev = _enlev_of(tmp_path, levels)
    energies = {n: (lev[3], lev[2], lev[4])
                for n, lev in zip(range(1, 5), levels)}
    energies, rep = sync.configuration_uncertainties(enlev, energies,
                                                     lambda *a: None)
    assert energies[1] == (5000.0, 5000.0, False)
    assert rep['changes'] == []
    row = [r for r in rep['table'] if r[0] == 'f5d6d'][0]
    assert row[1] == 0 and row[3] == 1 and row[4] is True


def test_too_few_found_levels_is_no_rms_at_all(tmp_path):
    """An rms over one level is not a measure of anything."""
    levels = [(1, 5000.0, 5000.0, 5000.0, False, 4.5, ' f5d6d_3H4H'),
              (2, 1000.0, 0.01, 1030.0, True, 4.5, ' f5d6d_3H4G')]
    enlev = _enlev_of(tmp_path, levels)
    energies = {1: (5000.0, 5000.0, False), 2: (1030.0, 0.01, True)}
    energies, rep = sync.configuration_uncertainties(enlev, energies,
                                                     lambda *a: None)
    assert energies[1][1] == 5000.0
    assert sync.MIN_CFG_LEVELS > 1


def test_a_change_of_uncertainty_is_reported(tmp_path):
    """The uncertainty of a found level follows max(D1, D2tot) of the LOPT
    run, and a change in it has to be reported even when the energy has not
    moved at all."""
    rows = [(1, 4, 20, IDEN.BLANK_OBS)]
    _, _, _, energies, changes = run_mini(
        tmp_path, rows, [(101, 104, 1.0e7)],
        lopt=dict(e3=10000.500, e4=120.250))
    # Level 4 stands where it stood; its uncertainty was 0.004 in the file
    # and max(D1, D2tot) = max(0.0040, 0.0020) = 0.004, so nothing moves.
    assert energies[4] == (120.250, 0.004, True)
    # Level 3: the file said 0.005 and D1 = 0.005 as well, so the energy and
    # the uncertainty both stand.  Move the level and the change is listed.
    _, _, _, energies, changes = run_mini(
        tmp_path / 'moved', rows, [(101, 104, 1.0e7)],
        lopt=dict(e3=10000.900, e4=120.250))
    moved = [c for c in changes if c[0] == 3][0]
    assert moved[2] == 10000.500 and moved[3] == 10000.900
    assert moved[4] == 0.005 and moved[5] == 0.005


def test_a_header_built_from_scratch_keeps_the_label_intact(tmp_path):
    """The spaces inside a label separate the configuration from the parent
    term, so the eleven characters between the slashes have to come through
    unaltered - stripping the leading one shifts the whole label left."""
    d = mini_iden2(tmp_path, [(1, 4, 20, IDEN.BLANK_OBS)])
    enlev = IDEN.Enlev(str(d / 'enlev.dat'))
    for lev in MINI_LEVELS:
        head = sync.build_header(lev[0], enlev)
        assert len(head) == sync.HEADER_WIDTH
        assert IDEN.field(head, sync.TH_LABEL) == lev[6]
        assert int(IDEN.field(head, IDEN.TH_INDEX)) == lev[0]


# ---------------------------------------------------------------------------
# The real files
# ---------------------------------------------------------------------------
@real_file
def test_the_real_levels_line_up_with_enlev():
    """Every level of the calculation finds its row in enlev.dat, and the
    594 that carry an experimental identifier confirm the alignment."""
    trans = cowan_gA.read_transitions()
    enlev = cowan_gA.read_enlev_levels(os.path.join(REAL_IDEN2, 'enlev.dat'))
    id_of_row = cowan_gA.read_id_map(
        os.path.join(REAL_IDEN2, 'IDEN_level_ids.txt'))
    mapping, report = cowan_gA.match_to_enlev(cowan_gA.levels(trans), enlev,
                                              id_of_row)
    assert report['n_matched'] == len(enlev)
    assert report['n_checked'] == len(id_of_row)
    assert not report['unmatched_idx']
    assert report['max_dE'] < 10.0


@real_file
def test_the_real_gA_is_the_gA_of_Icalc():
    """`Icalc.xlsx` is drawn from this workbook: where they hold the same
    transition they must hold the same strength."""
    icalc_path = os.path.join(HERE, 'Icalc.xlsx')
    if not os.path.exists(icalc_path):
        pytest.skip('Icalc.xlsx is absent')
    icalc = pd.read_excel(icalc_path)
    icalc['key'] = [frozenset((cowan_gA.level_id(a), cowan_gA.level_id(b)))
                    for a, b in zip(icalc['id1'], icalc['id2'])]
    trans = cowan_gA.read_transitions()
    both = trans[(trans['id1'] != '') & (trans['id2'] != '')]
    gA_of = dict(zip([frozenset((a, b))
                      for a, b in zip(both['id1'], both['id2'])],
                     both['gA']))
    got = icalc['key'].map(gA_of)
    assert got.notna().all()
    assert np.allclose(got.to_numpy(float), icalc['gA'].to_numpy(float),
                       rtol=1e-4)
