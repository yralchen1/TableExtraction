"""Tests for the registry of lines measured otherwise than their class says
(files.kappa_exceptions; hfs_correction.py).

Run from the LineClass directory:  python -m pytest tests -q

A row names one transition of one line and the class it was measured in -
head (kappa = 1), cg (kappa = 0) or a class of [hfs.kappa] - and may name
levels of [hfs] resolved_levels the line does not resolve, which then count
with their whole pattern.  The tests fix:

  * the reader, and the rows it refuses;
  * the model: the kappa, the D and the shift of the named transition, the
    other transitions of the line left alone, a head row giving up widths;
  * the classification: the candidate's prediction, the line's shift, the
    lines the registry is attached to and the entries it cannot place;
  * the calibration: the line is an anchor-class line with D multiplied by
    its kappa, resolved levels count S = 0 there as in the classification,
    and the main line of resolved companions is held at kappa = 1.
"""
import csv
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import classify_lines as cl                # noqa: E402
import config                              # noqa: E402
import hfs_correction as H                 # noqa: E402
import hfs_kappa                           # noqa: E402
import wavelength_calibration as W         # noqa: E402
from models import EnergyLevel, SpectralLine, Transition   # noqa: E402

KAPPA = (('flag', (1.0, 0.0)), ('plain_1974', (0.944, 0.014)),
         ('plain_1969', (0.633, 0.035)), ('c', (0.20, 0.12)))

LOW, UPP, OTHER, RESOLVED = ('059003.000100', '059003.000300',
                             '059003.000301', '059003.000303')

S_LOW = 2.5 * -0.04 * 3.5
S_UPP = 2.5 * 0.1 * hfs_kappa.A_SCALE * 2.5
S_OTHER = 2.5 * 0.02 * 1.5
S_RES = 2.5 * 0.5 * hfs_kappa.A_SCALE * 0.5

HEADER = '\t'.join(H.EXCEPTION_COLUMNS) + '\n'


def a_table(tmp_path):
    path = tmp_path / 'A.csv'
    path.write_text(
        'level_id,cfg,J,A_cm-1,u_A,source,n_flagged\n'
        f'{LOW},f26s,3.5,-0.0400,0.0040,flag interval,0\n'
        f'{UPP},f26p,2.5,+0.1000,0.0100,Reader & Sugar 1965 (semiempirical),0\n'
        f'{OTHER},f5d2,1.5,+0.0200,0.0030,flag interval,0\n'
        f'{RESOLVED},f26s,0.5,+0.5000,0.0100,composition,0\n',
        encoding='utf-8', newline='\n')
    return str(path)


def registry(tmp_path, *rows, name='exc.txt'):
    path = tmp_path / name
    path.write_text(HEADER + ''.join('\t'.join(r) + '\n' for r in rows),
                    encoding='utf-8', newline='\n')
    return str(path)


def make_model(tmp_path, *rows, satellites=''):
    return H.Model(config.HfsSettings(
        apply=True, A_levels=a_table(tmp_path), kappa=KAPPA,
        resolved_levels=(RESOLVED,), satellites=satellites,
        kappa_exceptions=registry(tmp_path, *rows)))


ROWS = (('30000.1234', LOW, UPP, 'cg', '', '10/5/2026', 'char w'),
        ('30000.1234', OTHER, UPP, 'head', '', '10/5/2026', 'its head'),
        ('31000.5', RESOLVED, UPP, 'cg', RESOLVED, '10/5/2026', 'between'),
        ('32000.0', LOW, UPP, 'c', '', '10/5/2026', 'complex in effect'),
        ('#33000.0', LOW, UPP, 'cg', '', '', 'a comment'))


# ---------------------------------------------------------------------------
# The reader
# ---------------------------------------------------------------------------
def test_the_registry_is_read(tmp_path):
    reg = H.read_kappa_exceptions(registry(tmp_path, *ROWS), (RESOLVED,),
                                  dict(KAPPA))
    assert len(reg) == 4
    assert reg.find(30000.1234, LOW, UPP).cls == 'cg'
    assert reg.find(30000.1234, OTHER, UPP).cls == 'head'
    assert reg.find(30000.1234, UPP, LOW) is None
    assert reg.find(31000.5, RESOLVED, UPP).unresolved == {RESOLVED}
    # a row written with fewer decimals names the line by them
    assert reg.find(31000.4987, RESOLVED, UPP) is not None
    assert reg.find(33000.0, LOW, UPP) is None            # the comment
    assert len(reg.rows_of(30000.1234)) == 2


def test_no_registry_is_an_empty_one(tmp_path):
    assert len(H.read_kappa_exceptions('')) == 0
    assert len(H.read_kappa_exceptions(str(tmp_path / 'none.txt'))) == 0


@pytest.mark.parametrize('bad, message', [
    (('30000.5', LOW, UPP, 'wide', '', '', ''), 'has the class'),
    (('30000.5', LOW, UPP, 'cg', LOW, '', ''), 'not in .hfs. resolved'),
    (('30000.5', LOW, UPP, 'cg', RESOLVED, '', ''), 'not one of its levels'),
    (('30000.5', '', UPP, 'cg', '', '', ''), 'needs low_id and upp_id'),
    (('30000.1234', LOW, UPP, 'head', '', '', ''), 'entered twice'),
])
def test_a_row_that_cannot_mean_what_it_says(tmp_path, bad, message):
    path = registry(tmp_path, ROWS[0], bad)
    with pytest.raises(ValueError, match=message):
        H.read_kappa_exceptions(path, (RESOLVED,), dict(KAPPA))


# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------
def test_the_named_transition_takes_the_written_kappa(tmp_path):
    m = make_model(tmp_path, *ROWS)
    assert m.kappa_value('cg') == (0.0, 0.0)
    assert m.kappa_value('head') == (1.0, 0.0)
    assert m.kappa_value('c') == (0.20, 0.12)
    # a plain 1974 line: the cg row moves its own transition only
    assert m.component_kappas('w', 30000.0, [0.6, 0.4], [False, False],
                              [(0.0, 0.0), None]) \
        == [(0.0, 0.0), (0.944, 0.014)]
    shift, _, kappa = m.line_shift('w', 30000.1234, [(LOW, UPP, 1.0)])
    assert shift == pytest.approx(S_UPP - S_LOW)          # (1 - 0) * D
    assert kappa == 0.0
    # a transition the registry does not name keeps the line's class
    shift, _, kappa = m.line_shift('w', 30000.1234, [(LOW, OTHER, 1.0)])
    assert shift == pytest.approx(0.056 * (S_OTHER - S_LOW))
    # a borrowed class takes that class's value
    shift, _, kappa = m.line_shift('', 32000.0, [(LOW, UPP, 1.0)])
    assert shift == pytest.approx(0.8 * (S_UPP - S_LOW))


def test_an_unresolved_level_counts_with_its_pattern(tmp_path):
    m = make_model(tmp_path, *ROWS)
    assert m.S(RESOLVED) == (0.0, 0.0)
    assert m.S(RESOLVED, pattern=True)[0] == pytest.approx(S_RES)
    assert m.D(RESOLVED, UPP)[0] == pytest.approx(S_UPP)
    assert m.D(RESOLVED, UPP, {RESOLVED})[0] == pytest.approx(S_UPP - S_RES)
    shift, u, _ = m.line_shift('d', 31000.5, [(RESOLVED, UPP, 1.0)])
    assert shift == pytest.approx(S_UPP - S_RES)
    assert u > 0
    # another line of the level ends on its sublevel, as before
    shift, _, _ = m.line_shift('d', 31500.0, [(RESOLVED, UPP, 1.0)])
    assert shift == pytest.approx(0.056 * S_UPP)


def test_a_head_row_makes_a_head_line(tmp_path):
    m = make_model(tmp_path, *ROWS)
    assert m.is_head_line(30000.1234)
    assert not m.is_head_line(31000.5)                    # cg only
    assert not m.is_head_line(29000.0)


def test_fit_terms(tmp_path):
    sat = tmp_path / 'sat.txt'
    sat.write_text('\t'.join(H.SATELLITE_COLUMNS) + '\n'
                   + f'40000.4\t40000.0\t{LOW}\t{UPP}\t1\t\t\t\n',
                   encoding='utf-8', newline='\n')
    m = make_model(tmp_path, *ROWS, satellites=str(sat))
    D = S_UPP - S_LOW
    assert m.fit_terms('w', 30000.1234, LOW, UPP) == ('flag', 0.0, 'cg')
    cls, d, fixed = m.fit_terms('', 32000.0, LOW, UPP)
    assert (cls, fixed) == ('flag', 'c') and d == pytest.approx(0.2 * D)
    cls, d, fixed = m.fit_terms('d', 31000.5, RESOLVED, UPP)
    assert (cls, d, fixed) == ('flag', 0.0, 'cg')
    # the main line of resolved companions: held at kappa = 1
    cls, d, fixed = m.fit_terms('', 40000.0, LOW, UPP)
    assert (cls, fixed) == ('flag', 'head') and d == pytest.approx(D)
    # an ordinary line: its class, and a resolved level counting S = 0
    cls, d, fixed = m.fit_terms('w', 35000.0, RESOLVED, UPP)
    assert (cls, fixed) == ('plain', '') and d == pytest.approx(S_UPP)
    assert m.fit_terms('c', 35000.0, LOW, UPP)[0] == 'c'


# ---------------------------------------------------------------------------
# The classification
# ---------------------------------------------------------------------------
def level(lid, E, S=0.0, S_pattern=None, u=0.0):
    return EnergyLevel(level_id=lid, energy=E, parity='e', J_str='1',
                       J_val=1.0, hfs_S=S, u_hfs_S=u,
                       hfs_S_pattern=S if S_pattern is None else S_pattern,
                       u_hfs_S_pattern=u)


def test_the_candidate_is_predicted_with_the_written_kappa():
    lo, up = level(LOW, 1000.0, S=-0.35), level(UPP, 31000.0, S=0.6)
    res = level(RESOLVED, 1000.0, S=0.0, S_pattern=0.3)
    line = SpectralLine(wavenumber=29999.2, wn_uncertainty=0.1,
                        intensity=1.0, line_character='w', hfs_factor=0.056,
                        hfs_exceptions={(LOW, UPP): (0.0, 0.0, frozenset()),
                                        (RESOLVED, UPP): (
                                            0.0, 0.0,
                                            frozenset({RESOLVED}))})
    named = Transition(lower_level=lo, upper_level=up, assigned_to=line)
    assert named.predicted_for() == pytest.approx(30000.0 - 0.95)
    assert named.observed_head == pytest.approx(29999.2 + 0.95)
    unres = Transition(lower_level=res, upper_level=up, assigned_to=line)
    assert unres.hfs_D == pytest.approx(0.6)
    assert unres.hfs_D_on(line) == pytest.approx(0.3)
    assert unres.predicted_for() == pytest.approx(30000.0 - 0.3)
    # not named: the line's factor and the plain D
    other = Transition(lower_level=level(OTHER, 1000.0, S=0.1),
                       upper_level=up, assigned_to=line)
    assert other.predicted_for() == pytest.approx(30000.0 - 0.056 * 0.5)
    # a flagged line (factor 0) still takes its exception
    line.hfs_factor = 0.0
    assert named.predicted_for() == pytest.approx(30000.0 - 0.95)


def switched_on(monkeypatch, m, path):
    monkeypatch.setattr(cl, 'HFS', m)
    monkeypatch.setattr(cl, 'HFS_MAX_D', 0.0)
    monkeypatch.setattr(cl, 'CFG', type('C', (), {
        'hfs': config.HfsSettings(kappa_exceptions=path)})())


def test_the_registry_is_attached_and_used(monkeypatch, tmp_path):
    m = make_model(tmp_path, *ROWS)
    switched_on(monkeypatch, m, str(tmp_path / 'exc.txt'))
    levels = {x: level(x, 0.0) for x in (LOW, UPP, OTHER, RESOLVED)}
    lines = [SpectralLine(30000.1234, 0.1, 1.0, 'w', wn_key=30000.1234),
             SpectralLine(31000.4999, 0.1, 1.0, 'd', wn_key=31000.4999),
             SpectralLine(32000.0, 0.1, 1.0, '', wn_key=32000.0),
             SpectralLine(39000.0, 0.1, 1.0, '', wn_key=39000.0)]
    assert cl.attach_kappa_exceptions(lines, levels) == 4
    assert lines[0].hfs_exceptions[(LOW, UPP)] == (0.0, 0.0, frozenset())
    assert lines[1].hfs_exceptions[(RESOLVED, UPP)][2] == {RESOLVED}
    assert lines[2].hfs_exceptions[(LOW, UPP)][:2] == (0.20, 0.12)
    assert lines[3].hfs_exceptions == {}
    # the line's shift: the cg transition and the head one
    lo, up = level(LOW, 0.0, S=-0.35, u=0.01), level(UPP, 0.0, S=0.6, u=0.02)
    lines[0].hfs_factor = 0.056
    w = {}
    for low in (lo, level(OTHER, 0.0, S=0.1)):
        t = Transition(lower_level=low, upper_level=up,
                       assigned_to=lines[0], accepted=1, calc_intensity=1.0)
        lines[0].assigned_transitions.append(t)
        w[id(t)] = 0.5 ** 2 / 0.1 ** 2      # BF 0.5 each: sqrt(w) * u = BF
    shift, _, kappa = cl.hfs_line_shift(lines[0], w)
    assert shift == pytest.approx(0.5 * 0.95 + 0.5 * 0.0)
    assert kappa == pytest.approx(0.5 * 0.0 + 0.5 * 1.0)
    # an entry whose transition is not accepted on its line is reported
    assert sorted(cl.unmet_kappa_exceptions(lines)) == sorted(
        [(32000.0, LOW, UPP), (31000.4999, RESOLVED, UPP)])


def test_an_entry_that_names_no_line_or_level_stops(monkeypatch, tmp_path):
    m = make_model(tmp_path, *ROWS)
    switched_on(monkeypatch, m, str(tmp_path / 'exc.txt'))
    levels = {x: level(x, 0.0) for x in (LOW, UPP, OTHER, RESOLVED)}
    lines = [SpectralLine(30000.1234, 0.1, 1.0, 'w', wn_key=30000.1234)]
    with pytest.raises(ValueError, match='no observed line'):
        cl.attach_kappa_exceptions(lines, levels)
    lines += [SpectralLine(31000.5, 0.1, 1.0, '', wn_key=31000.5),
              SpectralLine(32000.0, 0.1, 1.0, '', wn_key=32000.0)]
    del levels[OTHER]
    with pytest.raises(ValueError, match='not in this run'):
        cl.attach_kappa_exceptions(lines, levels)


# ---------------------------------------------------------------------------
# The calibration
# ---------------------------------------------------------------------------
def table(tmp_path, rows):
    path = tmp_path / 'lines.csv'
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, lineterminator='\n')
        w.writerow(['wn_obs', 'wn_key', 'char', 'low_id', 'upp_id',
                    'accepted', 'n_accepted'])
        w.writerows(rows)
    return str(path)


def test_the_calibration_sees_the_lines_as_the_classification(tmp_path):
    m = make_model(tmp_path, *ROWS)
    path = table(tmp_path, [
        ['30000.2', '30000.1234', 'w', LOW, UPP, '1.0', '1'],
        ['31000.6', '31000.5', 'd', RESOLVED, UPP, '1.0', '1'],
        ['35000.0', '35000.0', 'w', RESOLVED, UPP, '1.0', '1'],
        ['36000.0', '36000.0', 'c', LOW, OTHER, '1.0', '1'],
    ])
    got = {ln.key: ln for ln in hfs_kappa.read_lines(
        path, J_of={}, raw=True, model=m)}
    assert (got[30000.1234].cls, got[30000.1234].D,
            got[30000.1234].fixed) == ('flag', 0.0, 'cg')
    assert got[31000.5].fixed == 'cg'
    assert got[35000.0].D == pytest.approx(S_UPP)        # resolved: S = 0
    assert (got[36000.0].cls, got[36000.0].fixed) == ('c', '')
    # without a model, as before
    plain = {ln.key: ln for ln in hfs_kappa.read_lines(
        path, constants={}, J_of={}, raw=True)}
    assert plain[30000.1234].cls == 'plain' and plain[30000.1234].fixed == ''


def test_the_configuration_beside_the_table(tmp_path):
    t = tmp_path / 'lines.csv'
    t.write_text('', encoding='utf-8')
    assert W.config_for(str(t)) == config.DEFAULT_PATH
    (tmp_path / 'lineclass_config.toml').write_text('', encoding='utf-8')
    assert W.config_for(str(t)) == str(tmp_path / 'lineclass_config.toml')
    assert W.config_for(str(t), 'x.toml') == 'x.toml'


def test_the_baseline_names_the_registry_and_it_reads():
    cfg = config.load()
    assert os.path.basename(cfg.hfs.kappa_exceptions) == \
        'kappa_exceptions.txt'
    m = W.hfs_model_of(config.DEFAULT_PATH)
    for row in m.exceptions.rows:
        assert row.low_id in m.levels and row.upp_id in m.levels
