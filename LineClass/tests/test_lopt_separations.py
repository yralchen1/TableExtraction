"""Tests of lopt_separations.py: LOPT 6's separations file of blend components."""

import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import lopt_separations as S             # noqa: E402
import make_LOPT_input as M              # noqa: E402


def rec(wn, low, upp, flag='', weight=1.0):
    return M.format_line(wn, 0.05, 10.0, low, upp, flag, weight)


def write_lines(path, records):
    path.write_text(''.join(r + '\r\n' for r in records), encoding='ascii')


LOW = ['059003.0000%02d' % k for k in range(1, 10)]


def flat(wn):
    return 1.0                            # an FWHM of 1 cm-1 everywhere


def test_the_blends_are_lopts_own(tmp_path):
    # same wavenumber text = one blend; P records and zero weights are not
    # members; a lone record is no blend
    path = tmp_path / 'lines.txt'
    write_lines(path, [rec(30000.0, LOW[0], LOW[1], weight=0.6),
                       rec(30000.0, LOW[2], LOW[3], weight=0.4),
                       rec(30000.0, LOW[4], LOW[5], flag='P', weight=0.0),
                       rec(29000.0, LOW[0], LOW[3]),
                       rec(28000.0, LOW[1], LOW[2], weight=0.5),
                       rec(28000.0, LOW[4], LOW[5], weight=0.0)])
    got = S.blends(S.read_records(str(path)))
    assert [(t, [(r['low'], r['upp']) for r in rs]) for t, _, rs in got] == \
        [('30000.000', [(LOW[0], LOW[1]), (LOW[2], LOW[3])])]


def test_pair_uncertainties_are_lopts_weights():
    # two equal components: u = tau; three: tau * sqrt(3/2) per pair
    assert S.pair_uncertainties([0.2, 0.2]) == [(0, 1, pytest.approx(0.2))]
    three = S.pair_uncertainties([0.2] * 3)
    assert [p[:2] for p in three] == [(0, 1), (0, 2), (1, 2)]
    assert all(p[2] == pytest.approx(0.2 * math.sqrt(1.5)) for p in three)
    # unequal: sqrt(sum w / (w_i w_j)), w = 2 / tau^2
    w = [2 / 0.1 ** 2, 2 / 0.3 ** 2]
    assert S.pair_uncertainties([0.1, 0.3])[0][2] == \
        pytest.approx(math.sqrt(sum(w) / (w[0] * w[1])))


def test_a_wide_character_gets_the_wide_factor(tmp_path):
    path = tmp_path / 'lines.txt'
    write_lines(path, [rec(30000.0, LOW[0], LOW[1]), rec(30000.0, LOW[2], LOW[3]),
                       rec(20000.0, LOW[0], LOW[2]), rec(20000.0, LOW[1], LOW[3]),
                       rec(10000.0, LOW[4], LOW[5]), rec(10000.0, LOW[6], LOW[7])])
    chars = {30000.0: 'w**', 20000.0: 'cl', 10000.0: ''}
    text, stats = S.build(S.read_records(str(path)),
                          lambda t, wn: chars[round(wn, 3)], 0.5, 1.0,
                          ('w', 'c', 'd'), fwhm=flat)
    records = [l for l in text.splitlines() if l and not l.startswith(';')]
    # w** is w (the multiple-classification mark stripped); cl is not c
    assert [float(l.split()[-1]) for l in records] == [1.0, 0.5, 0.5]
    assert stats['blends'] == 3 and stats['wide'] == 1
    assert stats['records'] == 3 and stats['extra_dof'] == 0
    # four labels of the label field's width, then the value and u
    first, last = M.FIELDS['lower_level']
    w = last - first + 1
    assert records[0][:4 * w] == LOW[0] + LOW[1] + LOW[2] + LOW[3]
    assert records[0][4 * w:].split() == ['0', '1.0000']


def test_a_blend_of_three_adds_one_degree_of_freedom(tmp_path):
    path = tmp_path / 'lines.txt'
    write_lines(path, [rec(30000.0, LOW[k], LOW[k + 3]) for k in range(3)])
    _, stats = S.build(S.read_records(str(path)), lambda t, wn: '', 0.5, 1.0,
                       ('w',), fwhm=flat)
    assert stats['records'] == 3 and stats['extra_dof'] == 1


def test_the_effective_width_is_instrument_and_doppler():
    import level_shifts
    wn = 20000.0                                  # 5000 A
    want = math.hypot(level_shifts.INSTR_FWHM_A,
                      level_shifts.doppler_fwhm_factor() * 5000.0) * 4.0
    assert S.fwhm_cm(wn) == pytest.approx(want)


def test_refresh_follows_the_transitions_file(tmp_path):
    # a tool exchanges two levels: the rebuilt file names the new pairs, keeps
    # the rules and the character it recorded, and takes a new blend's
    # character from the table
    lines = tmp_path / 'LOPT_input_lines.txt'
    sep = tmp_path / S.DEF_NAME
    write_lines(lines, [rec(30000.0, LOW[0], LOW[1]),
                        rec(30000.0, LOW[2], LOW[3])])
    text, _ = S.build(S.read_records(str(lines)), lambda t, wn: 'd', 0.5, 1.0,
                      ('d',), fwhm=flat)
    S.write(str(sep), text)
    write_lines(lines, [rec(30000.0, LOW[0], LOW[4]),
                        rec(30000.0, LOW[2], LOW[3]),
                        rec(25000.0, LOW[5], LOW[6]),
                        rec(25000.0, LOW[7], LOW[8])])
    stats = S.refresh(str(lines), str(sep), {25000.0: 'w'}, fwhm=flat)
    body = [l for l in sep.read_text().splitlines()
            if l and not l.startswith(';')]
    assert body[0].startswith(LOW[0] + LOW[4] + LOW[2] + LOW[3])
    assert [float(l.split()[-1]) for l in body] == [1.0, 0.5]   # d wide, w not
    assert stats['blends'] == 2 and stats['unknown_char'] == 0


def test_refresh_without_a_file_does_nothing(tmp_path):
    assert S.refresh(str(tmp_path / 'x.txt'), str(tmp_path / 'none.txt')) \
        is None


def test_the_parameter_file_names_the_file_and_turns_lopts_own_off():
    par = ['a.txt ; TRANSITIONS input file name',
           'centroid; BLEND TREATMENT [component/centroid]']
    out = S.set_par_lines(par, 'LOPT_separations.txt')
    assert out[:2] == par
    assert out[2].split(';')[0].strip() == 'LOPT_separations.txt'
    assert S.PAR_SEP_FILE_KEY in out[2].upper()
    assert out[3].split(';')[0].strip() == '0'
    assert S.PAR_SEP_FACTOR_KEY in out[3].upper()
    # rewritten in place, not appended twice; None takes both keys out
    again = S.set_par_lines(out, 'other.txt')
    assert len(again) == 4 and again[2].startswith('other.txt')
    assert S.set_par_lines(out, None) == par


def test_separations_beside_reads_the_par_file(tmp_path):
    (tmp_path / 'LOPT.par').write_text(
        'x ; TRANSITIONS input file name\n'
        'LOPT_separations.txt ; SEPARATIONS INPUT FILE NAME\n')
    got = S.separations_beside(str(tmp_path / 'LOPT_input_lines.txt'))
    assert got == os.path.join(str(tmp_path), 'LOPT_separations.txt')
    (tmp_path / 'LOPT.par').write_text('x ; TRANSITIONS input file name\n')
    assert S.separations_beside(str(tmp_path / 'LOPT_input_lines.txt')) is None


def test_write_par_file_keeps_names_or_drops_the_separations(tmp_path):
    sample = tmp_path / 'sample.par'
    sample.write_text('a ; TRANSITIONS input file name\n'
                      'b ; FIXED LEVELS input file name\n'
                      'c ; Levels OUTPUT file name\n'
                      'd ; Transitions OUTPUT file name\n'
                      'centroid; BLEND TREATMENT\n')
    named, kept, dropped = (tmp_path / n for n in ('n.par', 'k.par', 'd.par'))
    M.write_par_file(str(sample), str(named), 'L', 'F', 'O1', 'O2',
                     separations=S.DEF_NAME)
    assert S.par_value(str(named), S.PAR_SEP_FILE_KEY) == S.DEF_NAME
    assert S.par_value(str(named), S.PAR_SEP_FACTOR_KEY) == '0'
    # the default copies the sample's keys as they are (lopt_perturb.py)
    M.write_par_file(str(named), str(kept), 'L', 'F', 'O1', 'O2')
    assert S.par_value(str(kept), S.PAR_SEP_FILE_KEY) == S.DEF_NAME
    M.write_par_file(str(named), str(dropped), 'L', 'F', 'O1', 'O2',
                     separations=None)
    assert S.par_value(str(dropped), S.PAR_SEP_FILE_KEY) is None
    assert S.par_value(str(dropped), S.PAR_SEP_FACTOR_KEY) is None


def test_a_registry_factor_overrides_the_character_and_survives_refresh(
        tmp_path):
    lines = tmp_path / 'LOPT_input_lines.txt'
    sep = tmp_path / S.DEF_NAME
    write_lines(lines, [rec(30000.0, LOW[0], LOW[1]),
                        rec(30000.0, LOW[2], LOW[3]),
                        rec(20000.0, LOW[4], LOW[5]),
                        rec(20000.0, LOW[6], LOW[7])])
    own = {30000.0: 2.5}
    text, stats = S.build(S.read_records(str(lines)), lambda t, wn: 'w', 0.5,
                          1.0, ('w',), fwhm=flat, own_of=lambda t, wn:
                          own.get(round(wn, 3)))
    assert stats['own'] == 1
    S.write(str(sep), text)

    def body():
        return [float(l.split()[-1]) for l in sep.read_text().splitlines()
                if l and not l.startswith(';')]
    assert body() == [2.5, 1.0]
    # a tool rewrites the file: the blend keeps its factor of its own
    S.refresh(str(lines), str(sep), fwhm=flat)
    assert body() == [2.5, 1.0]


def test_the_registry_is_read_and_checked(tmp_path):
    reg = tmp_path / 'blend_separations.txt'
    reg.write_text('# a comment\nwn_key\tfactor\tdate\treason\n'
                   '34861.88980\t3.0\t2026-11-18\tflagged\n', encoding='utf-8')
    assert S.read_registry(str(reg)) == {34861.8898: 3.0}
    assert S.read_registry('') == {}
    reg.write_text('wn_key\tfactor\tdate\treason\n1.0\t1\tx\ty\n1.0\t2\tx\ty\n',
                   encoding='utf-8')
    with pytest.raises(ValueError, match='twice'):
        S.read_registry(str(reg))


def test_the_project_registry_parses():
    path = os.path.join(os.path.dirname(HERE), 'blend_separations.txt')
    reg = S.read_registry(path)
    assert reg and all(f > 0 for f in reg.values())
