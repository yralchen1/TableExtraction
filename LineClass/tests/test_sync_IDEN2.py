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


# ---------------------------------------------------------------------------
# The observed lines: dlv.dat and the identifications
# ---------------------------------------------------------------------------
def dlv_record(code, wn, lam, character, u_lam, row):
    """One row of dlv.dat at its fixed width of 66 characters."""
    rec = ('%5d%14.3f%14.4f  /%10s/%13.4f%6d'
           % (code, wn, lam, character, u_lam, row))
    assert len(rec) == sync.DLV_WIDTH, len(rec)
    return rec


def test_read_lopt_transitions_lets_an_accepted_record_win():
    """A transition with one accepted record and one flagged P is accepted."""
    path = None
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False,
                                     newline='') as fh:
        path = fh.name
        for flag in ('P', ' '):
            fh.write(make_lopt_record(20000.125, 0.050, '059003.000004',
                                      '059003.000003', flag))
        fh.write(make_lopt_record(30000.500, 0.060, '059003.000004',
                                  '059003.000002', 'P'))
    try:
        got = sync.read_lopt_transitions(path)
    finally:
        os.unlink(path)
    assert got[('059003.000004', '059003.000003')] == (True, 20000.125)
    assert got[('059003.000004', '059003.000002')] == (False, 30000.500)


def make_lopt_record(wn, unc, low, upp, flag):
    """A record of LOPT_input_lines.txt in make_LOPT_input's own columns."""
    buf = [' '] * 93
    def put(span, text):
        buf[span[0]:span[0] + len(text)] = text
    put(sync.FIELD_WN, '%.3f' % wn)
    put(sync.FIELD_UNC, '%.3f' % unc)
    put(sync.FIELD_LOW, low)
    put(sync.FIELD_UPP, upp)
    put(sync.FIELD_FLAGS, flag)
    put((89, 93), 'cm-1')
    return ''.join(buf).rstrip() + '\r\n'


def test_dlv_gets_the_new_wavenumber_and_the_uncertainty_in_angstroms():
    """The three measured fields change and nothing else does.

    The uncertainty in dlv.dat is a wavelength uncertainty, so a line at
    20000 cm^-1 (5000 A) whose wavenumber is known to 0.2 cm^-1 is known to
    0.2 * 5000 / 20000 = 0.05 A.
    """
    old = dlv_record(62, 20000.100, 4999.9750, 'c', 0.0070, 5215)
    lines = [(20000.100, 20000.000, 0.2)]
    new, rep = sync.rewrite_dlv([old], lines, lambda *a: None)
    assert len(new) == 1 and len(new[0]) == sync.DLV_WIDTH
    got = new[0]
    assert float(got[sync.DLV_WN[0]:sync.DLV_WN[1]]) == 20000.000
    assert abs(float(got[sync.DLV_UNC[0]:sync.DLV_UNC[1]]) - 0.0500) < 5e-5
    # the wavelength is the standard air wavelength of the new wavenumber,
    # whatever the row carried before
    assert float(got[sync.DLV_LAMBDA[0]:sync.DLV_LAMBDA[1]]) ==         round(sync.wavelength_of(20000.000), 4)
    # the row number, the intensity and the character are untouched
    assert got[:sync.DLV_WN[0]] == old[:sync.DLV_WN[0]]
    assert got[sync.DLV_LAMBDA[1]:sync.DLV_UNC[0]] == \
        old[sync.DLV_LAMBDA[1]:sync.DLV_UNC[0]]
    assert got[sync.DLV_ROW[0]:] == old[sync.DLV_ROW[0]:]
    assert len(rep['changed']) == 1 and not rep['unmatched']


def test_a_dlv_row_with_no_line_is_left_alone_and_reported():
    old = dlv_record(62, 20000.100, 4999.9750, 'c', 0.0070, 7)
    new, rep = sync.rewrite_dlv([old], [(31000.000, 31000.000, 0.1)],
                                lambda *a: None)
    assert new == [old]
    assert rep['unmatched'] == [(7, 20000.100)]
    assert [round(a[0], 3) for a in rep['absent']] == [31000.000]


def test_a_dlv_row_already_on_the_corrected_scale_is_still_its_line():
    """After the first sync a corrected set's dlv.dat carries own_corr, which
    differs from wn_key by the calibration correction; the second sync must
    still find the row and refresh it."""
    old = dlv_record(62, 20000.100, 4999.9750, 'c', 0.0070, 7)
    lines = [(20000.000, 20000.100, 0.2), (20000.100, 20000.300, 0.1)]
    new, rep = sync.rewrite_dlv([old], lines, lambda *a: None)
    # the row carries the first line's current wavenumber, and the second
    # line's key: the current wavenumber is what it is
    assert float(new[0][sync.DLV_WN[0]:sync.DLV_WN[1]]) == 20000.100
    assert abs(float(new[0][sync.DLV_UNC[0]:sync.DLV_UNC[1]]) - 0.0500) < 5e-5
    assert not rep['unmatched']
    assert rep['absent'] == [(20000.100, 20000.300, 0.1)]


def test_a_row_inserted_in_iden2_under_a_later_line_number_is_matched():
    rows = [dlv_record(62, 30000.000, 3333.3333, 'c', 0.0010, 1),
            dlv_record(62, 25000.000, 4000.0000, 'c', 0.0010, 3),
            dlv_record(62, 20000.000, 5000.0000, 'c', 0.0010, 2)]
    lines = [(20000.0, 20000.0, 0.02), (25000.0, 25000.0, 0.02),
             (30000.0, 30000.0, 0.02)]
    new, rep = sync.rewrite_dlv(rows, lines, lambda *a: None)
    assert not rep['unmatched'] and not rep['absent']
    assert [r[sync.DLV_ROW[0]:] for r in new] == \
        [r[sync.DLV_ROW[0]:] for r in rows]


def test_the_standard_wavelength_is_vacuum_above_50000_and_air_below():
    """Peck and Reeder: n - 1 = 2.79e-4 at 20000 cm^-1 (5000 A), so the air
    wavelength is 1.3944 A shorter than the vacuum one; a row at 62500
    cm^-1 carries its vacuum wavelength 1600 A."""
    assert sync.wavelength_of(62500.0) == 1600.0
    assert abs(sync.air_index(20000.0) - 1.000279) < 1e-6
    assert abs(1e8 / 20000.0 - sync.wavelength_of(20000.0) - 1.3944) < 1e-4
    # row 6578 of the corrected set's dlv.dat as IDEN2 saved it on
    # 2026-09-30, its wavenumber rebuilt from its wavelength
    row = dlv_record(54, 11618.419, 8604.6590, '', 0.0138, 6578)
    assert sync.dispersion_disagreements([row]) == ([], 1)


def test_the_dlv_rewrite_changes_nothing_the_second_time():
    """The wavelength is a function of the line list's wavenumber alone, so
    a sync run on its own output writes exactly the same bytes.  The old
    rescaling moved 8604.6586 A to 8604.6590 and then 8604.6594 on an
    unmoved line, because it rescaled by a wavenumber rounded to 0.001."""
    lines = [(11618.4205, 11618.4205, 0.02), (29710.2085, 29710.2085, 0.01),
             (61680.8807, 61680.8807, 0.4)]
    rows = [dlv_record(62, 61680.881, 1621.2458, 'c', 0.0105, 1),
            dlv_record(62, 29710.209, 3364.8839, 'c', 0.0011, 2),
            dlv_record(62, 11618.421, 8604.6586, 'c', 0.0148, 3)]
    once, _rep = sync.rewrite_dlv(rows, lines, lambda *a: None)
    twice, rep = sync.rewrite_dlv(once, lines, lambda *a: None)
    assert twice == once and rep['changed'] == []


def test_every_rewritten_row_has_its_wavenumbers_wavelength():
    """IDEN2 rebuilds each air row's wavenumber from its wavelength when it
    saves the file, so a row written here must give back its own wavenumber
    to printing precision."""
    lines = [(w, w, 0.01) for w in (61680.8807, 49997.2141, 29710.2085,
                                   11618.4205)]
    rows = [dlv_record(62, w + 0.0017, 1e8 / w, 'c', 0.001, i + 1)
            for i, (_k, w, _u) in enumerate(lines)]
    out, _rep = sync.rewrite_dlv(rows, lines, lambda *a: None)
    assert sync.dispersion_disagreements(out) == ([], 4)
    for rec, (_k, w, _u) in zip(out, lines):
        lam = float(rec[sync.DLV_LAMBDA[0]:sync.DLV_LAMBDA[1]])
        wn = 1e8 / lam
        if w <= sync.VACUUM_ABOVE:
            for _ in range(8):
                wn = 1e8 / (lam * sync.air_index(wn))
        assert abs(wn - w) <= 0.00005 * wn / lam + 1e-9


def test_a_row_drifted_by_0_0017_is_still_its_line_and_is_put_back():
    """The case check_sync.py found on 2026-09-30: 53 rows 0.0015 to 0.0017
    from their line, too far for the old match of 0.0015."""
    old = dlv_record(62, 29710.207, 3364.8841, 'c', 0.0011, 5026)
    new, rep = sync.rewrite_dlv([old], [(29710.2085, 29710.2085, 0.01)],
                                lambda *a: None)
    assert not rep['unmatched'] and not rep['absent']
    assert float(new[0][sync.DLV_WN[0]:sync.DLV_WN[1]]) == 29710.209


def test_a_file_on_another_dispersion_formula_stops_the_sync():
    """If most rows disagree with the formula, the formula is not the one
    the file was written with, and nothing is rewritten with it."""
    rows = [dlv_record(62, w, 1e8 / w, 'c', 0.001, i + 1)
            for i, w in enumerate((30000.0, 25000.0, 20000.0))]
    with pytest.raises(sync.SyncError):
        sync.check_dispersion(rows)
    good = [dlv_record(62, w, sync.wavelength_of(w), 'c', 0.001, i + 1)
            for i, w in enumerate((30000.0, 25000.0, 20000.0))]
    assert sync.check_dispersion(good + rows[:1])[0][0][0] == 1


def test_a_recalibrated_line_is_found_through_the_keys_record():
    """The case check_sync.py found on 2026-10-03: a second calibration moved
    a line by 0.05 since dlv.dat was written, so the row shows neither its
    line's new wavenumber nor its wn_key.  Without the record it is left
    alone; with it, it is the line it was written for."""
    old = dlv_record(62, 20000.100, sync.wavelength_of(20000.100), 'c',
                     0.0050, 7)
    lines = [(20000.000, 20000.150, 0.1)]
    _new, rep = sync.rewrite_dlv([old], lines, lambda *a: None)
    assert rep['unmatched'] == [(7, 20000.100)]
    new, rep = sync.rewrite_dlv([old], lines, lambda *a: None,
                                keys={7: (20000.000, 20000.100)})
    assert not rep['unmatched'] and not rep['absent']
    assert float(new[0][sync.DLV_WN[0]:sync.DLV_WN[1]]) == 20000.150
    assert rep['n_by_key'] == 1
    assert rep['keyed'] == [(7, 20000.000, 20000.150)]


def test_a_row_edited_since_it_was_written_is_not_taken_by_its_record():
    """The record names the line a row was written for only while the row
    still shows what was written; a row moved in IDEN2 is found by what it
    shows now."""
    old = dlv_record(62, 25000.000, sync.wavelength_of(25000.000), 'c',
                     0.0050, 7)
    lines = [(20000.000, 20000.150, 0.1), (25000.000, 25000.000, 0.1)]
    new, rep = sync.rewrite_dlv([old], lines, lambda *a: None,
                                keys={7: (20000.000, 20000.100)})
    assert float(new[0][sync.DLV_WN[0]:sync.DLV_WN[1]]) == 25000.000
    assert rep['keyed'] == [(7, 25000.000, 25000.000)] and rep['n_by_key'] == 0


def test_an_unmatched_row_keeps_its_entry_in_the_keys_record():
    old = dlv_record(62, 20000.100, sync.wavelength_of(20000.100), 'c',
                     0.0050, 7)
    _new, rep = sync.rewrite_dlv([old], [(31000.0, 31000.0, 0.1)],
                                 lambda *a: None,
                                 keys={7: (20000.000, 20000.100)})
    assert rep['unmatched'] == [(7, 20000.100)]
    assert rep['keyed'] == [(7, 20000.000, 20000.100)]


def test_the_keys_record_reads_back_what_was_written(tmp_path):
    key = 15714.85007185795
    sync.write_keys(sync.keys_path(str(tmp_path)),
                    [(6213, key, 15714.848), (1, 30000.0, 30000.1)])
    assert sync.read_keys(str(tmp_path)) == {6213: (key, 15714.848),
                                             1: (30000.0, 30000.1)}
    assert sync.read_keys(str(tmp_path / 'nowhere')) == {}
    lines = (tmp_path / sync.KEYS_FILE).read_text().splitlines()
    assert lines[0] == 'line\twn_key\twn_written'
    assert lines[1].startswith('1\t')


REGISTRY = ('wn_key\tunc_wn\tdate\treason\n'
            '20000.1\t0.3\t\twidened in IDEN2\n'
            '30000.0\t0.01\t\tnarrower than the list\n')


def test_the_registry_widens_a_line_and_never_narrows_one(tmp_path):
    path = tmp_path / 'inflated_unc_lines.txt'
    path.write_text(REGISTRY)
    lines = [(20000.1004, 20000.2, 0.05), (25000.0, 25000.1, 0.05),
             (30000.0, 30000.1, 0.05)]
    got = sync.apply_registry(lines, str(path), lambda *a: None)
    assert got == [(20000.1004, 20000.2, 0.3), (25000.0, 25000.1, 0.05),
                   (30000.0, 30000.1, 0.05)]


def test_a_registry_entry_naming_no_line_stops_the_run(tmp_path):
    path = tmp_path / 'inflated_unc_lines.txt'
    path.write_text(REGISTRY)
    with pytest.raises(sync.SyncError, match='30000.0'):
        sync.apply_registry([(20000.1004, 20000.2, 0.05)], str(path),
                            lambda *a: None)


def test_no_registry_changes_nothing(tmp_path):
    lines = [(20000.1, 20000.2, 0.05)]
    assert sync.apply_registry(lines, '', lambda *a: None) == lines
    assert sync.apply_registry(lines, str(tmp_path / 'none.txt'),
                               lambda *a: None) == lines


def test_the_line_list_collapses_repeated_rows_of_one_line(tmp_path):
    """A blend named once per component is one row of dlv.dat."""
    path = tmp_path / 'lines.xlsx'
    pd.DataFrame({'own': [100.0, 100.0, 200.0],
                  'own_corr': [100.5, 100.5, 200.5],
                  'unc_own_corr': [0.01, 0.01, 0.02]}).to_excel(
                      path, index=False)
    cfg = FakeConfig(str(path), {'wn_key': 'own', 'wn': 'own_corr',
                                 'u_wn': 'unc_own_corr'})
    got = sync.read_line_list(cfg, lambda *a: None)
    assert got == [(100.0, 100.5, 0.01), (200.0, 200.5, 0.02)]


def test_repeated_rows_that_disagree_are_a_stop(tmp_path):
    """dlv.dat has one row per line and no way to hold two wavenumbers."""
    path = tmp_path / 'lines.xlsx'
    pd.DataFrame({'own': [100.0, 100.0],
                  'own_corr': [100.5, 100.9],
                  'unc_own_corr': [0.01, 0.01]}).to_excel(path, index=False)
    cfg = FakeConfig(str(path), {'wn_key': 'own', 'wn': 'own_corr',
                                 'u_wn': 'unc_own_corr'})
    with pytest.raises(sync.SyncError) as exc:
        sync.read_line_list(cfg, lambda *a: None)
    assert '100.000' in str(exc.value)


class FakeLayout(object):
    def __init__(self, columns):
        self.columns = columns
        self.sheet = None


class FakeConfig(object):
    def __init__(self, lines_file, columns):
        self.lines_file = lines_file
        self.lines = FakeLayout(columns)


# --- the identifications ----------------------------------------------------
ID3 = '059003.000003'
ID4 = '059003.000004'


def assignment_sync(tmp_path, rows, lopt_lines, dlv, keep_unlisted=False):
    """sync_assignments over a miniature IDEN2."""
    d = mini_iden2(tmp_path, rows)
    trans = IDEN.Trans(str(d / 'trans.dat'))
    id_of_row = cowan_gA.read_id_map(str(d / 'IDEN_level_ids.txt'))
    wanted, rep = sync.sync_assignments(
        trans, id_of_row, lopt_lines, sync.dlv_rows_by_wavenumber(dlv),
        keep_unlisted, lambda *a: None)
    return trans, wanted, rep


def test_an_excluded_transition_loses_its_identification(tmp_path):
    """Every record flagged P means the identification has been withdrawn."""
    tail = sync.make_assignment(62, 9880.250, 0.000, 11)
    trans, wanted, rep = assignment_sync(
        tmp_path, [(3, 4, 40, tail)],
        {(ID4, ID3): (False, 9880.250)},
        [dlv_record(62, 9880.250, 10121.1, 'c', 0.05, 11)])
    k = trans.row_of[(3, 4)]
    assert not IDEN.has_line(IDEN.assignment(trans.records[k]))
    assert [r[0] for r in rep['removed_p']] == [9880.250]
    assert not rep['removed_unlisted'] and not wanted


def test_a_transition_the_lopt_input_never_mentions_loses_it_too(tmp_path):
    tail = sync.make_assignment(62, 9880.250, 0.000, 11)
    trans, _wanted, rep = assignment_sync(
        tmp_path, [(3, 4, 40, tail)], {},
        [dlv_record(62, 9880.250, 10121.1, 'c', 0.05, 11)])
    k = trans.row_of[(3, 4)]
    assert not IDEN.has_line(IDEN.assignment(trans.records[k]))
    assert [r[0] for r in rep['removed_unlisted']] == [9880.250]


def test_keep_unlisted_leaves_a_hand_marked_line_where_it_is(tmp_path):
    """The switch for lines marked in IDEN2 and not yet classified."""
    tail = sync.make_assignment(62, 9880.250, 0.000, 11)
    trans, _wanted, rep = assignment_sync(
        tmp_path, [(3, 4, 40, tail)], {},
        [dlv_record(62, 9880.250, 10121.1, 'c', 0.05, 11)],
        keep_unlisted=True)
    k = trans.row_of[(3, 4)]
    assert IDEN.assignment(trans.records[k]) == tail
    assert not rep['removed_unlisted']


def test_an_accepted_transition_has_its_wavenumber_refreshed(tmp_path):
    """A corrected set moves the observed wavenumber; the row must follow."""
    tail = sync.make_assignment(62, 9880.250, 0.000, 11)
    trans, _wanted, rep = assignment_sync(
        tmp_path, [(3, 4, 40, tail)],
        {(ID4, ID3): (True, 9880.200)},
        [dlv_record(77, 9880.200, 10121.2, 'c', 0.05, 11)])
    obs = IDEN.assignment(trans.records[trans.row_of[(3, 4)]])
    assert IDEN.obs_wavenumber(obs) == 9880.200
    assert IDEN.obs_row(obs) == 11
    assert [r[1] for r in rep['refreshed']] == [9880.200]


def test_an_accepted_transition_with_no_row_is_added(tmp_path):
    """The transition is in the fit and not on the screen at all."""
    trans, wanted, rep = assignment_sync(
        tmp_path, [(3, 4, 40, IDEN.BLANK_OBS)],
        {(ID4, ID3): (True, 9880.250)},
        [dlv_record(62, 9880.250, 10121.1, 'c', 0.05, 11)])
    obs = IDEN.assignment(trans.records[trans.row_of[(3, 4)]])
    assert IDEN.obs_wavenumber(obs) == 9880.250
    assert not wanted                     # the row was already there
    assert [r[0] for r in rep['added']] == [9880.250]


def test_a_level_the_lookup_table_has_no_row_for_is_a_stop(tmp_path):
    """The rule for IDEN2 lookups: report it, never guess."""
    with pytest.raises(sync.SyncError) as exc:
        assignment_sync(tmp_path, [(3, 4, 40, IDEN.BLANK_OBS)],
                        {('059003.000999', ID3): (True, 9880.250)}, [])
    assert '059003.000999' in str(exc.value)


@real_file
def test_the_real_dlv_uncertainty_is_a_wavelength_uncertainty():
    """u_lambda * wn^2 / 1e8 is the line's uncertainty in cm^-1.

    This is what the rewrite of dlv.dat rests on, and it is stated nowhere in
    the file, so it is checked against the line list the baseline set uses.
    """
    import config
    cfg = config.load(os.path.join(HERE, 'lineclass_config.toml'))
    if not os.path.exists(cfg.lines_file):
        pytest.skip('the line list is absent')
    lines = sync.read_line_list(cfg, lambda *a: None)
    u_of = {round(a, 3): c for a, _b, c in lines}
    records, _ends = IDEN.read_records(os.path.join(REAL_IDEN2, 'dlv.dat'))
    checked, agree = 0, 0
    for rec in records:
        if len(rec) < sync.DLV_WIDTH or not rec.strip():
            continue
        wn = float(rec[sync.DLV_WN[0]:sync.DLV_WN[1]])
        u_lam = float(rec[sync.DLV_UNC[0]:sync.DLV_UNC[1]])
        want = u_of.get(round(wn, 3))
        if want is None:
            continue
        checked += 1
        if abs(u_lam * wn * wn / 1e8 - want) < 0.002 + 0.02 * want:
            agree += 1
    assert checked > 6000, checked
    # All but a handful: the exceptions are lines whose uncertainty was
    # inflated by hand in the line list and never written back here, which is
    # one of the things this rewrite is for.
    assert agree > checked - 20, (agree, checked)
