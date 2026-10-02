"""Tests for the head-frame hyperfine correction ([hfs] apply).

Run from the LineClass directory:  python -m pytest tests -q

The correction (hfs_correction.py, Work_on_hfs_plan.md Step 5) predicts a
line that is candidate transition i at Ritz_i - (1 - kappa) * D_i, and gives
LOPT the line at wn_obs + (1 - kappa) * sum BF_i * D_i.  The tests fix:

  * that nothing changes with the switch off - no S, no factor, no column;
  * the three groups of levels (determined, undetermined, absent) and the
    resolved ones, and the classes of kappa;
  * the arithmetic of a line's shift and its uncertainty, blends included;
  * the widths and registry allowances the correction replaces (D13);
  * the LOPT input: one shifted wavenumber per observed line, the record of
    the shifts, and a table written under the other setting refused;
  * that sync_IDEN2 and check_sync take the measured wavenumber back;
  * the registry of resolved companions: its reader, the kappa = 1 of the
    transition it names on the main line, the companion's row and its
    absence from the LOPT input.
"""
import io
import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import check_sync                          # noqa: E402
import classify_lines as cl                # noqa: E402
import config                              # noqa: E402
import hfs_correction as H                 # noqa: E402
import hfs_patterns                        # noqa: E402
import make_LOPT_input as M                # noqa: E402
import sync_IDEN2 as sync                  # noqa: E402
from models import EnergyLevel, SpectralLine, Transition   # noqa: E402

KAPPA = (('flag', (1.0, 0.0)), ('plain_1974', (0.944, 0.014)),
         ('plain_1969', (0.633, 0.035)), ('c', (0.20, 0.12)))

LOW, UPP, UNDET, ABSENT, RESOLVED = (
    '059003.000100', '059003.000300', '059003.000301', '059003.000302',
    '059003.000303')


def a_table(tmp_path):
    """A table of A constants: two determined levels, one undetermined, one
    determined but resolved; ABSENT is not listed."""
    path = tmp_path / 'A.csv'
    path.write_text(
        'level_id,cfg,J,A_cm-1,u_A,source,n_flagged\n'
        f'{LOW},f26s,3.5,-0.0400,0.0040,composition,0\n'
        f'{UPP},f26p,2.5,+0.1000,0.0100,Reader & Sugar 1965 (calculated),0\n'
        f'{UNDET},f5d2,4.5,+0.0000,0.0500,not determined,1\n'
        f'{RESOLVED},f26s,0.5,+0.5000,0.0100,composition,0\n',
        encoding='utf-8', newline='\n')
    return str(path)


@pytest.fixture
def model(tmp_path):
    return H.Model(config.HfsSettings(apply=True, A_levels=a_table(tmp_path),
                                      kappa=KAPPA,
                                      resolved_levels=(RESOLVED,)))


# ---------------------------------------------------------------------------
# The configuration
# ---------------------------------------------------------------------------
def write_config(tmp_path, body, name='set.toml'):
    base = os.path.join(ROOT, 'lineclass_config.toml').replace('\\', '/')
    path = tmp_path / name
    path.write_text(f'inherit = "{base}"\n' + body, encoding='utf-8',
                    newline='\n')
    return str(path)


def test_the_baseline_has_it_off():
    cfg = config.load()
    assert cfg.hfs.apply is False
    assert dict(cfg.hfs.kappa)['plain_1969'] == (0.633, 0.035)


def test_a_set_turns_it_on_with_one_line(tmp_path):
    cfg = config.load(write_config(tmp_path, '[hfs]\napply = true\n'))
    assert cfg.hfs.apply is True
    # everything else is inherited, the path resolved against the baseline
    assert cfg.hfs.A_levels == os.path.join(ROOT, 'A_hfs_levels.csv')
    assert cfg.hfs.kappa_of('flag') == (1.0, 0.0)
    assert '059003.000642' in cfg.hfs.resolved_levels


def test_on_without_a_table_of_constants_is_refused(tmp_path):
    path = write_config(tmp_path, '[files]\nhfs_A_levels = ""\n'
                                  '[hfs]\napply = true\n')
    with pytest.raises(config.ConfigError, match='hfs_A_levels'):
        config.load(path)


def test_a_kappa_class_must_be_a_pair(tmp_path):
    path = write_config(tmp_path, '[hfs.kappa]\nc = 0.2\n')
    with pytest.raises(config.ConfigError, match='two numbers'):
        config.load(path)


# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------
def test_the_groups_of_levels(model):
    assert model.S(LOW) == pytest.approx((2.5 * -0.04 * 3.5, 2.5 * 3.5 * 0.004))
    # undetermined: no displacement, the unknown A carried as uncertainty
    assert model.S(UNDET) == pytest.approx((0.0, 2.5 * 4.5 * 0.05))
    assert model.S(ABSENT) == (0.0, 0.0)
    # resolved: its own A counts as nothing, whatever the table says
    assert model.S(RESOLVED) == (0.0, 0.0)
    assert [model.is_corrected(x) for x in (LOW, UNDET, ABSENT, RESOLVED)] \
        == [True, False, False, True]
    assert [model.is_undetermined(x)
            for x in (LOW, UNDET, ABSENT, RESOLVED)] \
        == [False, True, False, False]


def test_the_classes_of_kappa(model):
    assert model.kappa_of('*r', 30000.0)[0] == 1.0
    assert model.kappa_of('*v', 90000.0)[0] == 1.0
    assert model.kappa_of('c', 30000.0)[0] == 0.20
    assert model.kappa_of('', 30000.0)[0] == 0.944       # 1974
    assert model.kappa_of('w', 60000.0)[0] == 0.633      # 1969
    assert model.kappa_of(None, 47499.9)[0] == 0.944


def test_a_blend_takes_one_shift(model):
    d1 = model.D(LOW, UPP)[0]                 # 0.625 + 0.35
    d2 = model.D(LOW, ABSENT)[0]              # 0.35
    shift, u, kappa = model.line_shift(
        '', 60000.0, [(LOW, UPP, 0.75), (LOW, ABSENT, 0.25)])
    assert shift == pytest.approx((1 - 0.633) * (0.75 * d1 + 0.25 * d2))
    assert u > 0
    assert kappa == 0.633
    # a flagged line sits on the head-frame Ritz value
    assert model.line_shift('*r', 60000.0, [(LOW, UPP, 1.0)])         == (0.0, 0.0, 1.0)
    # a line with nothing accepted is not corrected
    assert model.line_shift('c', 60000.0, []) == (0.0, 0.0, 0.20)


def test_a_flagged_blend_takes_the_flag_on_its_strongest_component(model):
    d1 = model.D(LOW, UPP)[0]
    d2, u_d2 = model.D(LOW, ABSENT)
    assert model.component_kappas('*v', 30000.0, [0.25, 0.75])         == [(0.944, 0.014), (1.0, 0.0)]
    # a single flagged line, and an unflagged blend, keep their one kappa
    assert model.component_kappas('*v', 30000.0, [1.0]) == [(1.0, 0.0)]
    assert model.component_kappas('', 60000.0, [0.4, 0.6])         == [(0.633, 0.035)] * 2
    shift, u, kappa = model.line_shift(
        '*v', 30000.0, [(LOW, UPP, 0.75), (LOW, ABSENT, 0.25)])
    # only the weaker component, the plain one, is shifted
    assert shift == pytest.approx(0.25 * (1 - 0.944) * d2)
    assert u == pytest.approx(math.hypot(0.25 * (1 - 0.944) * u_d2,
                                         0.25 * 0.014 * d2))
    # and the line's kappa is the BF-weighted mean
    assert kappa == pytest.approx(0.75 * 1.0 + 0.25 * 0.944)
    assert d1 != d2


def test_the_uncertainty_of_a_shift(model):
    d, u_d = model.D(LOW, UNDET)
    shift, u, _ = model.line_shift('c', 30000.0, [(LOW, UNDET, 1.0)])
    assert u == pytest.approx(math.hypot(0.8 * u_d, 0.12 * d))


def test_hfs_reasons():
    assert H.is_hfs_reason('hfs of f26s')
    assert H.is_hfs_reason('unresolved HFS')
    assert not H.is_hfs_reason('Sugar measured the position badly')
    assert not H.is_hfs_reason('hfsx')


# ---------------------------------------------------------------------------
# The objects classify_lines works with
# ---------------------------------------------------------------------------
def level(lid, E, S=0.0):
    return EnergyLevel(level_id=lid, energy=E, parity='e', J_str='1',
                       J_val=1.0, hfs_S=S)


def test_with_the_switch_off_nothing_moves():
    lo, up = level(LOW, 1000.0, S=-0.35), level(UPP, 31000.0, S=0.625)
    line = SpectralLine(wavenumber=30000.3, wn_uncertainty=0.1,
                        intensity=1.0, line_character='')
    t = Transition(lower_level=lo, upper_level=up, assigned_to=line)
    # hfs_factor 0: the plain Ritz value and the plain measurement, exactly
    assert t.predicted_for(line) == t.calculated_wavenumber == 30000.0
    assert t.observed_head == 30000.3


def test_a_candidate_is_predicted_below_the_head_frame_ritz():
    lo, up = level(LOW, 1000.0, S=-0.35), level(UPP, 31000.0, S=0.625)
    line = SpectralLine(wavenumber=29999.7, wn_uncertainty=0.1,
                        intensity=1.0, line_character='', hfs_factor=0.367)
    t = Transition(lower_level=lo, upper_level=up, assigned_to=line)
    assert t.hfs_D == pytest.approx(0.975)
    assert t.predicted_for() == pytest.approx(30000.0 - 0.367 * 0.975)
    assert t.observed_head == pytest.approx(29999.7 + 0.367 * 0.975)


def switched_on(monkeypatch, model):
    monkeypatch.setattr(cl, 'HFS', model)
    monkeypatch.setattr(cl, 'HFS_MAX_D', 0.0)


def test_apply_hfs_model(monkeypatch, model):
    switched_on(monkeypatch, model)
    levels = {x: EnergyLevel(level_id=x, energy=0.0, parity='e',
                             J_str=j, J_val=float(eval(j)))
              for x, j in ((LOW, '7/2'), (UPP, '5/2'), (ABSENT, '1'))}
    decoy = EnergyLevel(level_id=cl.DECOY_PREFIX + UPP, energy=0.0,
                        parity='e', J_str='5/2', J_val=2.5, is_decoy=1)
    levels[decoy.level_id] = decoy
    lines = [SpectralLine(30000.0, 0.1, 1.0, '*r', wn_key=30000.0),
             SpectralLine(60000.0, 0.1, 1.0, '', wn_key=60000.0)]
    cl.apply_hfs_model(levels, lines)
    assert levels[UPP].hfs_S == pytest.approx(0.625)
    assert decoy.hfs_S == levels[UPP].hfs_S      # the decoy copies its level
    assert levels[ABSENT].hfs_S == 0.0
    assert [ln.hfs_factor for ln in lines] == pytest.approx([0.0, 0.367])
    assert cl.HFS_MAX_D == pytest.approx(1.25)


def test_a_table_out_of_step_with_the_levels_stops(monkeypatch, model):
    switched_on(monkeypatch, model)
    levels = {LOW: EnergyLevel(level_id=LOW, energy=0.0, parity='e',
                               J_str='5/2', J_val=2.5)}      # table says 7/2
    with pytest.raises(ValueError, match='out of date'):
        cl.apply_hfs_model(levels, [])


def accepted(line, lo, up, weight_of):
    t = Transition(lower_level=lo, upper_level=up, assigned_to=line,
                   accepted=1, calc_intensity=1.0)
    line.assigned_transitions.append(t)
    weight_of[id(t)] = None
    return t


def test_the_allowance_goes_where_the_correction_accounts(monkeypatch,
                                                          model):
    switched_on(monkeypatch, model)
    lo, up, un = level(LOW, 0.0), level(UPP, 0.0), level(UNDET, 0.0)
    ok = SpectralLine(30000.0, 0.3, 1.0, '', hfs_factor=0.056,
                      unc_before_hfs_allowance=0.02)
    kept = SpectralLine(31000.0, 0.3, 1.0, '', hfs_factor=0.056,
                        unc_before_hfs_allowance=0.02)
    flag = SpectralLine(32000.0, 0.3, 1.0, '*r', hfs_factor=0.0,
                        unc_before_hfs_allowance=0.02)
    w = {}
    accepted(ok, lo, up, w)
    accepted(kept, lo, un, w)       # touches a level without A
    accepted(flag, lo, un, w)       # but a flag sits on the Ritz value
    assert cl.release_hfs_allowances([ok, kept, flag]) == 2
    assert (ok.wn_uncertainty, kept.wn_uncertainty, flag.wn_uncertainty) \
        == (0.02, 0.3, 0.02)


def test_the_shift_written_for_a_line(monkeypatch, model):
    switched_on(monkeypatch, model)
    lo = level(LOW, 0.0, S=model.S(LOW)[0])
    up = level(UPP, 0.0, S=model.S(UPP)[0])
    ab = level(ABSENT, 0.0)
    line = SpectralLine(60000.0, 0.2, 1.0, '', hfs_factor=0.367,
                        hfs_u_kappa=0.035, wn_key=60000.0)
    w = {}
    t1 = accepted(line, lo, up, w)
    t2 = accepted(line, lo, ab, w)
    # calc_weights' own form: w = BF^2 / u^2
    w[id(t1)], w[id(t2)] = 0.75 ** 2 / 0.04, 0.25 ** 2 / 0.04
    shift, _, kappa = cl.hfs_line_shift(line, w)
    assert shift == pytest.approx(0.367 * (0.75 * t1.hfs_D + 0.25 * t2.hfs_D))
    assert kappa == 0.633
    # flagged, the strongest component (t1) takes kappa = 1
    line.line_character = '*v'
    shift, _, kappa = cl.hfs_line_shift(line, w)
    assert shift == pytest.approx(0.25 * 0.367 * t2.hfs_D)
    assert kappa == pytest.approx(0.75 + 0.25 * 0.633)


# ---------------------------------------------------------------------------
# The LOPT input
# ---------------------------------------------------------------------------
def row(wn, low, upp, shift, kappa, accepted=1, calc=1.0, unc='0.100',
        u_shift='0.0300', char='', u_undet='0', allowance='0'):
    return {'wn_obs': wn, 'unc_wn_obs': unc, 'obs_intens': 100.0,
            'char': char, 'low_id': low, 'upp_id': upp, 'accepted': accepted,
            'calc_intens': calc, 'kappa': kappa, 'hfs_D': '0.5',
            'hfs_shift': shift, 'u_hfs_shift': u_shift,
            'u_hfs_undet': u_undet, 'hfs_allowance': allowance}


def read_back(path):
    out = []
    for rec in io.open(path, encoding='ascii', newline=''):
        rec = rec.rstrip('\r\n')
        got = {k: rec[a - 1:b].strip() for k, (a, b) in M.FIELDS.items()}
        out.append(got)
    return out


def test_the_lopt_input_carries_one_shift_per_line(tmp_path, model):
    rows = [row('60000.1234', LOW, UPP, '0.3000', '0.633', calc=3.0),
            row('60000.1234', LOW, ABSENT, '0.3000', '0.633', calc=1.0),
            row('60000.1234', UPP, ABSENT, '0.3000', '0.633', accepted=0),
            row('30000.0000', LOW, UNDET, '0.0', '1.0')]
    path = str(tmp_path / 'lines.txt')
    shifts = []
    M.write_lines_file(rows, path, {}, hfs=model, shifts_out=shifts)
    got = read_back(path)
    assert {g['wavenumber'] for g in got} == {'60000.423', '30000.000'}
    # the shared uncertainty has the shift's in it
    unc = {g['uncertainty'] for g in got if g['wavenumber'] == '60000.423'}
    assert unc == {'%.3f' % math.hypot(0.1, 0.03)}
    # the record: every record that moved, with the table's own text
    assert sorted((s[0], s[1], s[2], s[3]) for s in shifts) == sorted(
        (lo, up, '60000.1234', 60000.423)
        for lo, up in ((LOW, UPP), (LOW, ABSENT), (UPP, ABSENT)))


def test_widths_the_correction_replaces(tmp_path, model):
    w_hfs = {LOW: 0.2, UNDET: 0.3, ABSENT: 0.4}
    rows = [row('60000.0000', LOW, ABSENT, '0.0', '0.633', u_shift='0'),
            row('50000.0000', LOW, UNDET, '0.0', '1.0', u_shift='0',
                char='*r')]
    stats = {}
    path = str(tmp_path / 'lines.txt')
    M.write_lines_file(rows, path, w_hfs, hfs=model, hfs_stats=stats)
    unc = {g['wavenumber']: g['uncertainty'] for g in read_back(path)}
    # LOW's width is replaced; ABSENT keeps its own
    assert unc['60000.000'] == '%.3f' % math.hypot(0.1, 0.4)
    # a flagged line keeps none
    assert unc['50000.000'] == '0.100'
    assert stats['widths_left_out'] == 2


def test_a_level_without_A_is_counted_once(tmp_path, model):
    w_hfs = {LOW: 0.2, UNDET: 0.3, ABSENT: 0.4}
    rows = [
        # the width (0.3) is larger than the unknown A (0.2): the width
        row('60000.0000', LOW, UNDET, '0.0', '0.633', u_shift='0.25',
            u_undet='0.2'),
        # the unknown A (0.45) is larger: u_hfs_shift as it is
        row('50000.0000', UPP, UNDET, '0.0', '0.20', u_shift='0.5',
            u_undet='0.45', char='c'),
        # a registry allowance stands for both; ABSENT keeps its width
        row('40000.0000', UNDET, ABSENT, '0.0', '0.633', unc='0.500',
            u_shift='0.2', u_undet='0.2', allowance='1')]
    stats = {}
    path = str(tmp_path / 'lines.txt')
    M.write_lines_file(rows, path, w_hfs, hfs=model, hfs_stats=stats)
    unc = {g['wavenumber']: g['uncertainty'] for g in read_back(path)}
    assert unc['60000.000'] == '%.3f' % math.sqrt(0.1 ** 2 + 0.15 ** 2
                                                  + 0.3 ** 2)
    assert unc['50000.000'] == '%.3f' % math.hypot(0.1, 0.5)
    assert unc['40000.000'] == '%.3f' % math.hypot(0.5, 0.4)
    assert (stats['undet_counted_once'], stats['undet_by_allowance']) \
        == (2, 1)


def test_the_columns_for_a_level_without_A(monkeypatch, model):
    switched_on(monkeypatch, model)
    lo = level(LOW, 0.0, S=model.S(LOW)[0])
    up = level(UPP, 0.0, S=model.S(UPP)[0])
    un = level(UNDET, 0.0)
    for lev in (lo, up, un):
        lev.u_hfs_S = model.S(lev.level_id)[1]
        lev.hfs_undetermined = model.is_undetermined(lev.level_id)
    line = SpectralLine(30000.0, 0.3, 1.0, 'c', wn_key=30000.0,
                        hfs_factor=0.8,
                        unc_before_hfs_allowance=0.02)
    known = SpectralLine(31000.0, 0.3, 1.0, 'c', wn_key=31000.0,
                         hfs_factor=0.8,
                         unc_before_hfs_allowance=0.02)
    w = {}
    w[id(accepted(line, lo, un, w))] = 1.0 / 0.09
    w[id(accepted(known, lo, up, w))] = 1.0 / 0.09
    # the part owed to UNDET is (1 - kappa) * I*J*u_A of it, nothing else
    assert cl.hfs_undetermined_part(line, w) == pytest.approx(
        0.8 * 2.5 * 4.5 * 0.05)
    assert cl.hfs_undetermined_part(known, w) == 0.0
    # the allowance stays on the line touching UNDET, and only there
    assert cl.release_hfs_allowances([line, known]) == 1
    assert [cl.hfs_allowance_kept(x) for x in (line, known)] == [True, False]


def test_the_record_of_shifts_round_trips(tmp_path):
    lopt = tmp_path / 'LOPT_input_lines.txt'
    H.write_shifts(H.shifts_path(str(lopt)),
                   [(LOW, UPP, '60000.1234', 60000.423, 0.3, 0.03)])
    shifts = H.read_shifts(str(lopt))
    assert H.measured(LOW, UPP, 60000.423, shifts) == 60000.1234
    # LOPT prints its output rounded; the rounding is kept
    assert H.measured(LOW, UPP, 60000.42, shifts, tol=0.05) \
        == pytest.approx(60000.1204)
    # an unlisted record, or one at another wavenumber, is as it is
    assert H.measured(UPP, LOW, 60000.423, shifts) == 60000.423
    assert H.measured(LOW, UPP, 61000.0, shifts) == 61000.0
    assert H.read_shifts(str(tmp_path / 'elsewhere' / 'x.txt')) == {}


def test_a_table_written_under_the_other_setting_is_refused(tmp_path):
    on = write_config(tmp_path, '[hfs]\napply = true\n', 'on.toml')
    off = write_config(tmp_path, '', 'off.toml')
    plain = [{'wn_obs': '1', 'low_id': LOW, 'upp_id': UPP}]
    corrected = [dict(plain[0], kappa='1', hfs_D='0', hfs_shift='0',
                      u_hfs_shift='0', u_hfs_undet='0', hfs_allowance='0')]
    # a table written before the levels without A had their columns
    before = [dict(plain[0], kappa='1', hfs_D='0', hfs_shift='0',
                   u_hfs_shift='0')]
    with pytest.raises(SystemExit, match='counted twice'):
        M.hfs_model(on, before, 'table')
    with pytest.raises(SystemExit, match='written\\s+with it off'):
        M.hfs_model(on, plain, 'table')
    with pytest.raises(SystemExit, match='has it off'):
        M.hfs_model(off, corrected, 'table')
    assert M.hfs_model(off, plain, 'table') is None
    assert M.hfs_model(on, corrected, 'table') is not None


def test_rows_of_one_line_must_agree():
    rows = [row('60000.0', LOW, UPP, '0.3', '0.633'),
            row('60000.0', LOW, ABSENT, '0.2', '0.633')]
    with pytest.raises(SystemExit, match='different hfs_shift'):
        M.hfs_line_shifts(rows)


# ---------------------------------------------------------------------------
# The readers downstream
# ---------------------------------------------------------------------------
def lopt_file(tmp_path, records):
    path = tmp_path / 'LOPT_input_lines.txt'
    with open(path, 'w', encoding='ascii', newline='') as fh:
        for wn, low, upp in records:
            fh.write(M.format_line(wn, 0.1, 100.0, low, upp, '', 1.0)
                     + '\r\n')
    return str(path)


def test_sync_and_check_sync_see_the_measured_wavenumber(tmp_path):
    path = lopt_file(tmp_path, [(60000.423, LOW, UPP),
                                (30000.0, LOW, ABSENT)])
    H.write_shifts(H.shifts_path(path),
                   [(LOW, UPP, '60000.1234', 60000.423, 0.3, 0.03)])
    got = sync.read_lopt_transitions(path)
    assert got[(LOW, UPP)] == (True, 60000.1234)
    assert got[(LOW, ABSENT)] == (True, 30000.0)
    got = check_sync.read_lopt_input(path)
    assert got[check_sync.pair_key(LOW, UPP)][0][0] == 60000.1234
    assert got[check_sync.pair_key(LOW, ABSENT)][0][0] == 30000.0


def test_check_sync_reads_the_measured_value_from_lopt_output(tmp_path):
    lopt = lopt_file(tmp_path, [(60000.423, LOW, UPP)])
    H.write_shifts(H.shifts_path(lopt),
                   [(LOW, UPP, '60000.1234', 60000.423, 0.3, 0.03)])
    out = tmp_path / 'LOPT_output_lines.txt'
    out.write_text('wn_o\tL1\tL2\tE1\tE2\tWeight\n'
                   f'60000.42\t{LOW}\t{UPP}\t0\t60000.4\t1.000\n',
                   encoding='latin-1', newline='')
    got = check_sync.read_lopt_output_lines(str(out))
    assert got[check_sync.pair_key(LOW, UPP)][0][0] == pytest.approx(
        60000.1204)


# ---------------------------------------------------------------------------
# The registry of resolved companions (files.hfs_satellites)
# ---------------------------------------------------------------------------
SAT_HEADER = 'wn_key\tmain_wn_key\tlow_id\tupp_id\trung\tdate\treason\n'


def satellites(tmp_path, *rows, name='sat.txt'):
    path = tmp_path / name
    path.write_text(SAT_HEADER + ''.join('\t'.join(r) + '\n' for r in rows),
                    encoding='utf-8', newline='\n')
    return str(path)


def test_the_registry_of_companions(tmp_path):
    path = satellites(
        tmp_path,
        ('30000.4567', '30000.1234', LOW, UPP, '1', '9/30/2026', 'rung 1'),
        ('30000.8', '30000.1234', LOW, UPP, '2', '9/30/2026', 'rung 2'),
        ('#31000.0', '31000.5', LOW, ABSENT, '1', '', 'a comment'))
    sats = H.read_satellites(path)
    assert len(sats) == 2
    row = sats.companion(30000.45671)
    assert (row.main_key, row.low_id, row.upp_id, row.rung) \
        == ('30000.1234', LOW, UPP, 1)
    # an entry is matched at the precision it was written with
    assert sats.companion(30000.8049).rung == 2
    assert sats.head_pairs(30000.12344) == {(LOW, UPP)}
    assert sats.head_pairs(31000.5) == frozenset()     # the comment row
    sats.check([30000.45671, 30000.8049, 30000.12344, 29000.0], 'sat.txt')
    with pytest.raises(ValueError, match='no observed line'):
        sats.check([30000.45671, 30000.12344], 'sat.txt')
    with pytest.raises(ValueError, match='more than one line'):
        sats.check([30000.45671, 30000.8049, 30000.79, 30000.12344],
                   'sat.txt')
    assert len(H.read_satellites(str(tmp_path / 'none.txt'))) == 0
    assert len(H.read_satellites('')) == 0


@pytest.mark.parametrize('bad, message', [
    (('30000.1', '30000.1', LOW, UPP, '1'), 'its own companion'),
    (('30000.5', '30000.1', LOW, UPP, '0'), 'rung of 1 or more'),
    (('30000.5', '30000.1', '', UPP, '1'), 'rung of 1 or more'),
    (('30000.1234', '29000.0', LOW, UPP, '1'), 'both as a companion'),
    (('30000.4567', '30000.2', LOW, UPP, '2'), 'entered twice'),
])
def test_a_registry_row_that_cannot_mean_what_it_says(tmp_path, bad,
                                                      message):
    path = satellites(tmp_path,
                      ('30000.4567', '30000.1234', LOW, UPP, '1', '', ''),
                      bad + ('', ''))
    with pytest.raises(ValueError, match=message):
        H.read_satellites(path)


def test_a_marked_component_takes_kappa_one(model):
    # a plain line: the marked transition 1, the other the line's own
    assert model.component_kappas('', 30000.0, [0.3, 0.7], [True, False]) \
        == [(1.0, 0.0), (0.944, 0.014)]
    # a flagged blend: the mark, not the larger BF, says which is the head
    assert model.component_kappas('*r', 30000.0, [0.3, 0.7],
                                  [True, False]) \
        == [(1.0, 0.0), (0.944, 0.014)]
    # no mark: as before
    assert model.component_kappas('', 30000.0, [1.0], [False]) \
        == [(0.944, 0.014)]


def test_the_main_line_sits_on_the_head_frame_ritz_value():
    lo, up = level(LOW, 1000.0, S=-0.35), level(UPP, 31000.0, S=0.625)
    ab = level(ABSENT, 1000.0, S=0.2)
    line = SpectralLine(wavenumber=29999.9, wn_uncertainty=0.1,
                        intensity=1.0, line_character='', hfs_factor=0.056,
                        hfs_head_pairs=frozenset({(LOW, UPP)}))
    marked = Transition(lower_level=lo, upper_level=up, assigned_to=line)
    other = Transition(lower_level=ab, upper_level=up, assigned_to=line)
    assert marked.predicted_for() == marked.calculated_wavenumber
    assert marked.observed_head == 29999.9
    assert other.predicted_for() == pytest.approx(
        30000.0 - 0.056 * other.hfs_D)


def test_the_main_line_takes_kappa_one_in_its_shift(monkeypatch, model):
    switched_on(monkeypatch, model)
    lo = level(LOW, 0.0, S=model.S(LOW)[0])
    up = level(UPP, 0.0, S=model.S(UPP)[0])
    ab = level(ABSENT, 0.0)
    line = SpectralLine(30000.0, 0.2, 1.0, '', hfs_factor=0.056,
                        hfs_u_kappa=0.014, wn_key=30000.0,
                        hfs_head_pairs=frozenset({(LOW, UPP)}))
    w = {}
    t1 = accepted(line, lo, up, w)
    t2 = accepted(line, up, ab, w)
    w[id(t1)], w[id(t2)] = 0.25 ** 2 / 0.04, 0.75 ** 2 / 0.04
    shift, _, kappa = cl.hfs_line_shift(line, w)
    assert shift == pytest.approx(0.75 * 0.056 * t2.hfs_D)
    assert kappa == pytest.approx(0.25 + 0.75 * 0.944)
    # and the registry allowance goes, as a flagged line's does
    line.unc_before_hfs_allowance = 0.02
    assert cl.hfs_accounted(line)


def lines_and_levels():
    levels = {x: level(x, E) for x, E in ((LOW, 1000.0), (UPP, 31000.0),
                                           (ABSENT, 2000.0))}
    main = SpectralLine(30000.0, 0.02, 100.0, '', wn_key=29999.9876)
    comp = SpectralLine(30000.45, 0.02, 60.0, '', wn_key=30000.4376)
    other = SpectralLine(29000.0, 0.02, 60.0, '', wn_key=29000.0)
    return levels, [main, comp, other]


def test_attach_hfs_satellites(tmp_path, monkeypatch):
    monkeypatch.setattr(cl, 'HFS', None)
    levels, (main, comp, other) = lines_and_levels()
    path = satellites(tmp_path, ('30000.4376', '29999.9876', LOW, UPP, '1',
                                 '9/30/2026', 'rung 1'))
    assert cl.attach_hfs_satellites([main, comp, other], levels, path) == 1
    assert main.hfs_head_pairs == {(LOW, UPP)}
    assert comp.hfs_companion == (main, levels[LOW], levels[UPP], 1)
    assert other.hfs_companion is None and not other.hfs_head_pairs
    # the row it is written as: its transition, grade hfs, never accepted
    df = cl.build_output([main, comp, other], {})
    r = df[df.wn_obs == 30000.45].iloc[0]
    assert (r.low_id, r.upp_id, r.grade, r.accepted) \
        == (LOW, UPP, H.COMPANION_GRADE, 0)
    assert r['dif_wn_O-C'] == pytest.approx(0.45)
    assert r.notes2 == 'hfs companion, rung 1, of the line 29999.9876'
    assert 'kappa' not in df.columns        # the switch is off


def test_a_companion_cannot_be_an_identification_too(tmp_path, monkeypatch):
    monkeypatch.setattr(cl, 'HFS', None)
    levels, (main, comp, other) = lines_and_levels()
    path = satellites(tmp_path, ('30000.4376', '29999.9876', LOW, UPP, '1',
                                 '', ''))
    comp.decisions[(ABSENT, UPP)] = ('accept', 'by hand')
    with pytest.raises(ValueError, match='identified as'):
        cl.attach_hfs_satellites([main, comp, other], levels, path)
    del levels[LOW]
    comp.decisions.clear()
    with pytest.raises(ValueError, match='not in this run'):
        cl.attach_hfs_satellites([main, comp, other], levels, path)


def test_the_lopt_input_leaves_a_companion_out(tmp_path):
    path = tmp_path / 'lc.csv'
    path.write_text(
        'wn_obs,wn_key,low_id,upp_id,grade,accepted\n'
        f'30000.45,30000.4376,{LOW},{UPP},hfs,0\n'
        f'30000.0,29999.9876,{LOW},{UPP},2A,1\n'
        '29000.0,29000.0,,,,\n', encoding='utf-8', newline='\n')
    rows = M.read_classifications(str(path))
    assert [r['wn_obs'] for r in rows] == ['30000.0']


# --- blended companions (column blend) ----------------------------------------
BLEND_HEADER = ('wn_key\tmain_wn_key\tlow_id\tupp_id\trung\tblend\tdate\t'
                'reason\n')


def blended(tmp_path, *rows):
    path = tmp_path / 'blend.txt'
    path.write_text(BLEND_HEADER + ''.join('\t'.join(r) + '\n' for r in rows),
                    encoding='utf-8', newline='\n')
    return str(path)


def test_the_registry_reads_the_blend_column(tmp_path):
    path = blended(tmp_path,
                   ('30000.4376', '29999.9876', LOW, UPP, '1', '1', '', ''),
                   ('30000.8', '29999.9876', LOW, UPP, '2', '', '', ''),
                   ('30000.9', '29999.9876', LOW, UPP, '3', '0', '', ''))
    sats = H.read_satellites(path)
    assert [r.blend for r in sats.rows] == [True, False, False]
    # a file without the column has no blended companion
    plain = satellites(tmp_path, ('30000.4', '30000.1', LOW, UPP, '1', '',
                                  ''))
    assert not H.read_satellites(plain).rows[0].blend
    bad = blended(tmp_path, ('30000.4', '30000.1', LOW, UPP, '1', 'yes', '',
                             ''))
    with pytest.raises(ValueError, match='blend'):
        H.read_satellites(bad)


def test_the_share_of_a_rung():
    for J1, J2 in ((3.5, 4.5), (4.5, 4.5), (4.5, 3.5), (1.5, 0.5)):
        comps = hfs_patterns.components(J1, 0.0, J2, 0.0)
        total = sum(c[1] for c in comps)
        for k in range(H.hfs_patterns.ladder_length(J1, J2) + 1):
            f1, f2 = H.I_SPIN + J1 - k, H.I_SPIN + J2 - k
            want = sum(s for _, s, a, b in comps if (a, b) == (f1, f2))
            assert H.rung_share(J1, J2, k) == pytest.approx(want / total)
    # the strongest component of 7/2 - 9/2 carries a quarter of the line, its
    # first rung 0.197; past the ladder there is nothing
    assert H.rung_share(3.5, 4.5, 0) == pytest.approx(0.25, abs=5e-4)
    assert H.rung_share(3.5, 4.5, 1) == pytest.approx(0.197, abs=5e-4)
    assert H.rung_share(1.5, 0.5, 2) == 0.0


def blend_lines_and_levels():
    levels = {x: EnergyLevel(level_id=x, energy=E, parity='e', J_str=j,
                             J_val=float(eval(j)))
              for x, E, j in ((LOW, 1000.0, '7/2'), (UPP, 31000.0, '9/2'),
                              (ABSENT, 1000.45, '7/2'))}
    main = SpectralLine(30000.0, 0.02, 100.0, '', wn_key=29999.9876)
    comp = SpectralLine(30000.45, 0.02, 60.0, '', wn_key=30000.4376)
    return levels, main, comp


def test_a_blended_companion_keeps_its_own_identification(tmp_path,
                                                          monkeypatch):
    monkeypatch.setattr(cl, 'HFS', None)
    levels, main, comp = blend_lines_and_levels()
    path = blended(tmp_path, ('30000.4376', '29999.9876', LOW, UPP, '1', '1',
                              '', ''))
    # Sugar's identification of the companion, and a ledger accept: both
    # stand beside the hfs component
    own = Transition(lower_level=levels[ABSENT], upper_level=levels[UPP],
                     assigned_to=comp)
    comp.original_assignments.append(own)
    comp.decisions[(ABSENT, UPP)] = ('accept', 'by hand')
    assert cl.attach_hfs_satellites([main, comp], levels, path) == 1
    assert main.hfs_head_pairs == {(LOW, UPP)}
    assert comp.hfs_companion is None
    assert comp.hfs_blend_companion == (main, levels[LOW], levels[UPP], 1)

    # the classification: the main transition on the main line, and the
    # companion's own accepted on it
    t_main = Transition(lower_level=levels[LOW], upper_level=levels[UPP],
                        assigned_to=main, accepted=1, calc_intensity=400.0,
                        grade='2A')
    main.assigned_transitions.append(t_main)
    t_own = Transition(lower_level=levels[ABSENT], upper_level=levels[UPP],
                       assigned_to=comp, accepted=1, calc_intensity=20.0,
                       grade='2B')
    comp.assigned_transitions.append(t_own)
    rung_i = 400.0 * H.rung_share(3.5, 4.5, 1)
    assert cl.hfs_rung_intensity(comp) == pytest.approx(rung_i)
    share = 20.0 / (20.0 + rung_i)
    w = cl.calc_weights([main, comp])
    assert w[id(t_main)] == pytest.approx(1.0 / 0.02 ** 2)
    assert w[id(t_own)] == pytest.approx(share ** 2 / 0.02 ** 2)

    # the table: the companion's own row, and the component beside it
    df = cl.build_output([main, comp], w)
    rows = df[df.wn_obs == 30000.45].set_index('grade')
    assert set(rows.index) == {'2B', H.COMPANION_GRADE}
    assert rows.loc['2B', 'BF'] == pytest.approx(share)
    assert rows.loc['2B', 'accepted'] == 1
    hfs_row = rows.loc[H.COMPANION_GRADE]
    assert (hfs_row.low_id, hfs_row.upp_id, hfs_row.accepted) \
        == (LOW, UPP, 0)
    assert hfs_row.calc_intens == pytest.approx(rung_i)
    assert hfs_row.BF == pytest.approx(1.0 - share)
    assert hfs_row.notes2.endswith("blended with the line's own")


def test_a_blended_companion_without_a_transition_of_its_own(tmp_path,
                                                             monkeypatch):
    monkeypatch.setattr(cl, 'HFS', None)
    levels, main, comp = blend_lines_and_levels()
    path = blended(tmp_path, ('30000.4376', '29999.9876', LOW, UPP, '1', '1',
                              '', ''))
    cl.attach_hfs_satellites([main, comp], levels, path)
    df = cl.build_output([main, comp], {})
    rows = df[df.wn_obs == 30000.45]
    # one row, the component's: no empty row for the line beside it
    assert list(rows.grade) == [H.COMPANION_GRADE]
    assert rows.iloc[0].BF == 0.0


def test_the_shift_of_a_blended_companion_ignores_the_component(
        monkeypatch, model):
    switched_on(monkeypatch, model)
    lo = level(LOW, 0.0, S=model.S(LOW)[0])
    up = level(UPP, 0.0, S=model.S(UPP)[0])
    line = SpectralLine(30000.0, 0.2, 1.0, '', hfs_factor=0.056,
                        hfs_u_kappa=0.014, wn_key=30000.0,
                        hfs_blend_companion=(None, lo, up, 1))
    w = {}
    t = accepted(line, lo, up, w)
    w[id(t)] = 0.4 ** 2 / 0.2 ** 2          # the component took 0.6
    shift, _, _ = cl.hfs_line_shift(line, w)
    assert shift == pytest.approx(0.056 * t.hfs_D)


def test_the_lopt_input_widens_a_blended_companion(tmp_path):
    path = tmp_path / 'lc.csv'
    path.write_text(
        'wn_obs,wn_key,low_id,upp_id,grade,accepted,calc_intens\n'
        f'30000.45,30000.4376,{LOW},{UPP},hfs,0,60.0\n'
        f'30000.45,30000.4376,{ABSENT},{UPP},2B,1,20.0\n'
        f'30000.0,29999.9876,{LOW},{UPP},2A,1,400.0\n'
        f'29000.0,29000.0,{LOW},{ABSENT},hfs,0,\n',
        encoding='utf-8', newline='\n')
    comps = M.read_hfs_components(str(path))
    assert comps == {'30000.4376': 60.0}      # a pure companion has none
    rows = [dict(row('30000.45', ABSENT, UPP, '0', '0.944', calc=20.0),
                 wn_key='30000.4376'),
            dict(row('30000.0', LOW, UPP, '0', '0.944', calc=400.0),
                 wn_key='29999.9876')]
    out = str(tmp_path / 'lines.txt')
    M.write_lines_file(rows, out, components=comps)
    got = {g['wavenumber']: g for g in read_back(out)}
    # the line's own transition keeps weight 1, LOPT normalizing it anyway,
    # and its uncertainty is divided by its share, 20 / (20 + 60)
    assert float(got['30000.450']['weight']) == 1.0
    assert got['30000.450']['uncertainty'] == '0.400'
    assert got['30000.000']['uncertainty'] == '0.100'
    assert M.accepted_share(rows[:1], 0.0) == 1.0


def test_review_mismatches_leaves_a_companion_out_of_the_pairs(tmp_path):
    import review_mismatches as rm
    path = tmp_path / 'lc.csv'
    # the companion lies below its main line, so it comes after it
    path.write_text(
        'wn_obs,wn_key,low_id,upp_id,grade,accepted,calc_intens\n'
        f'30000.0,29999.9876,{LOW},{UPP},2A,1,400.0\n'
        f'29999.55,29999.5376,{LOW},{UPP},hfs,0,60.0\n',
        encoding='utf-8', newline='\n')
    _, by_line, by_pair = rm.read_classifications(str(path))
    assert by_pair[(LOW, UPP)]['wn_obs'] == '30000.0'
    assert len(by_line['29999.5376']) == 1


def test_the_main_line_gives_up_its_widths(tmp_path):
    sat = satellites(tmp_path, ('60000.4', '60000.0', LOW, ABSENT, '1', '',
                                ''))
    m = H.Model(config.HfsSettings(apply=True, A_levels=a_table(tmp_path),
                                   kappa=KAPPA, satellites=sat))
    rows = [dict(row('60000.0000', LOW, ABSENT, '0.0', '1.0', u_shift='0'),
                 wn_key='60000.0'),
            dict(row('50000.0000', UPP, ABSENT, '0.0', '0.944', u_shift='0'),
                 wn_key='50000.0')]
    path = str(tmp_path / 'lines.txt')
    M.write_lines_file(rows, path, {ABSENT: 0.4}, hfs=m)
    unc = {g['wavenumber']: g['uncertainty'] for g in read_back(path)}
    assert unc['60000.000'] == '0.100'
    assert unc['50000.000'] == '%.3f' % math.hypot(0.1, 0.4)


def test_the_baseline_names_the_registry():
    cfg = config.load()
    assert os.path.basename(cfg.hfs.satellites) == 'hfs_satellites.txt'
    sats = H.read_satellites(cfg.hfs.satellites)
    assert len(sats) > 0
