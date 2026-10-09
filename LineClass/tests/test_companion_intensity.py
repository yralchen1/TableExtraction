"""Tests for the light of resolved hfs companions added to their main lines.

Run from the LineClass directory:  python -m pytest tests -q

Where Sugar gave a resolved part of a pattern an intensity of its own, the
main line's printed intensity is only part of the transition, while the
calculated intensity it is tested against is all of it.  The tests fix:

  * that an unblended companion's whole intensity is added to its main line,
    and the companion keeps its own value;
  * that a blended companion gives its main line the hfs component's share
    of its light, by the calculated intensities of the component and of the
    line's own identifications (published unless rejected, ledger accepts);
  * that nothing is given when the head is not observed, or when the share
    cannot be computed;
  * the table column obs_intens_hfs;
  * that wavelength_calibration corrects a line of the line list that no
    classification table holds yet.
"""
import os
import sys

import openpyxl
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (ROOT, os.path.dirname(ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

import classify_lines as cl                 # noqa: E402
import hfs_correction as H                  # noqa: E402
import wavelength_calibration as W          # noqa: E402
from models import EnergyLevel, SpectralLine, Transition   # noqa: E402

LOW, UPP, OWN = '059003.000100', '059003.000300', '059003.000101'
HEADER = 'wn_key\tmain_wn_key\tlow_id\tupp_id\trung\tblend\tdate\treason\n'


def registry(tmp_path, *rows):
    path = tmp_path / 'sat.txt'
    path.write_text(HEADER + ''.join('\t'.join(r) + '\n' for r in rows),
                    encoding='utf-8', newline='\n')
    return str(path)


def levels():
    return {x: EnergyLevel(level_id=x, energy=E, parity='e', J_str=j,
                           J_val=float(eval(j)))
            for x, E, j in ((LOW, 1000.0, '7/2'), (UPP, 31000.0, '9/2'),
                            (OWN, 1000.45, '7/2'))}


def lines():
    main = SpectralLine(30000.0, 0.02, 100.0, '*r', wn_key=29999.9876)
    comp = SpectralLine(29999.8, 0.02, 40.0, '', wn_key=29999.7876)
    other = SpectralLine(29000.0, 0.02, 60.0, '', wn_key=29000.0)
    return main, comp, other


def test_an_unblended_companion_gives_all_its_light(tmp_path, monkeypatch):
    monkeypatch.setattr(cl, 'HFS', None)
    lev = levels()
    main, comp, other = lines()
    path = registry(tmp_path, ('29999.7876', '29999.9876', LOW, UPP, '1', '',
                               '', ''))
    cl.attach_hfs_satellites([main, comp, other], lev, path)
    assert main.intensity == pytest.approx(140.0)
    assert main.intensity_from_companions == pytest.approx(40.0)
    assert comp.intensity == pytest.approx(40.0)       # kept as measured
    assert other.intensity == pytest.approx(60.0)
    df = cl.build_output([main, comp, other], {})
    got = df.drop_duplicates('wn_key').set_index('wn_key')
    assert got.loc[29999.9876, 'obs_intens'] == pytest.approx(140.0)
    assert got.loc[29999.9876, 'obs_intens_hfs'] == pytest.approx(40.0)
    assert got.loc[29999.7876, 'obs_intens_hfs'] == 0.0
    assert got.loc[29000.0, 'obs_intens_hfs'] == 0.0


def test_two_companions_of_one_line(tmp_path, monkeypatch):
    monkeypatch.setattr(cl, 'HFS', None)
    lev = levels()
    main, comp, other = lines()
    comp2 = SpectralLine(29999.6, 0.02, 10.0, '', wn_key=29999.5876)
    path = registry(tmp_path,
                    ('29999.7876', '29999.9876', LOW, UPP, '1', '', '', ''),
                    ('29999.5876', '29999.9876', LOW, UPP, '2', '', '', ''))
    cl.attach_hfs_satellites([main, comp, comp2, other], lev, path)
    assert main.intensity == pytest.approx(150.0)
    assert main.intensity_from_companions == pytest.approx(50.0)


def test_a_headless_companion_gives_nothing(tmp_path, monkeypatch):
    monkeypatch.setattr(cl, 'HFS', None)
    lev = levels()
    main, comp, other = lines()
    path = registry(tmp_path, ('29999.7876', '', LOW, UPP, '1', '', '', ''))
    cl.attach_hfs_satellites([main, comp, other], lev, path)
    assert main.intensity == pytest.approx(100.0)
    assert main.intensity_from_companions == 0.0


def blended_setup(tmp_path):
    lev = levels()
    main, comp, other = lines()
    path = registry(tmp_path, ('29999.7876', '29999.9876', LOW, UPP, '1', '1',
                               '', ''))
    index = {(LOW, UPP): {'calc_intensity': 400.0},
             (OWN, UPP): {'calc_intensity': 20.0}}
    sugar = Transition(lower_level=lev[OWN], upper_level=lev[UPP],
                       assigned_to=comp, calc_intensity=20.0)
    comp.original_assignments.append(sugar)
    return lev, main, comp, other, path, index


def test_a_blended_companion_gives_the_component_s_share(tmp_path,
                                                         monkeypatch):
    monkeypatch.setattr(cl, 'HFS', None)
    lev, main, comp, other, path, index = blended_setup(tmp_path)
    cl.attach_hfs_satellites([main, comp, other], lev, path, index)
    rung_i = 400.0 * H.rung_share(3.5, 4.5, 1)
    part = 40.0 * rung_i / (rung_i + 20.0)
    assert main.intensity_from_companions == pytest.approx(part)
    assert main.intensity == pytest.approx(100.0 + part)
    assert comp.intensity == pytest.approx(40.0)


def test_a_rejected_identification_does_not_share(tmp_path, monkeypatch):
    monkeypatch.setattr(cl, 'HFS', None)
    lev, main, comp, other, path, index = blended_setup(tmp_path)
    comp.decisions[(OWN, UPP)] = ('reject', 'not this one')
    comp.decisions[(LOW, OWN)] = ('accept', 'by hand')
    index[(LOW, OWN)] = {'calc_intensity': 60.0}
    cl.attach_hfs_satellites([main, comp, other], lev, path, index)
    rung_i = 400.0 * H.rung_share(3.5, 4.5, 1)
    assert main.intensity_from_companions == pytest.approx(
        40.0 * rung_i / (rung_i + 60.0))


def test_no_share_without_calculated_intensities(tmp_path, monkeypatch):
    monkeypatch.setattr(cl, 'HFS', None)
    lev, main, comp, other, path, index = blended_setup(tmp_path)
    # without the index the main transition's intensity is unknown
    cl.attach_hfs_satellites([main, comp, other], lev, path)
    assert main.intensity_from_companions == 0.0
    # an own identification without one leaves the split undetermined
    lev, main, comp, other, path, index = blended_setup(tmp_path)
    del index[(OWN, UPP)]
    comp.original_assignments[0].calc_intensity = None
    cl.attach_hfs_satellites([main, comp, other], lev, path, index)
    assert main.intensity_from_companions == 0.0


# ---------------------------------------------------------------------------
# wavelength_calibration: a line the table does not hold yet
# ---------------------------------------------------------------------------
def line_list(tmp_path):
    path = str(tmp_path / 'lines.xlsx')
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(['own', 'unc_own', 'Icor', 'Ch.'])
    ws.append([30000.1234, 0.01, 100.0, '*r'])
    ws.append([30000.5, 0.3, 50.0, 'c'])
    ws.append([31000.0, 0.01, 10.0, None])
    ws.append([31000.0, 0.01, 10.0, '**'])     # its second identification
    wb.save(path)
    return path


def test_the_line_list_is_read_once_per_line(tmp_path):
    listed = W.read_line_list(line_list(tmp_path))
    # the first row of a line with two, as classify_lines.py reads it
    assert listed == [(30000.1234, '*r'), (30000.5, 'c'), (31000.0, '')]


def test_the_calibration_takes_lines_no_table_holds(tmp_path):
    listed = W.read_line_list(line_list(tmp_path))
    observed = [(30000.1234, '*r', 30000.1234)]
    got = W.with_unclassified(observed, listed)
    assert got == [(30000.1234, '*r', 30000.1234), (30000.5, 'c', 30000.5),
                   (31000.0, '', 31000.0)]


TABLE = ('wn_obs,wn_key,char,low_id,upp_id,accepted,n_accepted\n'
         '30000.14,30000.1234,bl,059003.000100,059003.000300,1.0,1\n'
         '30000.52,30000.5,c,,,,0\n')


def test_the_calibration_takes_the_character_from_the_line_list(tmp_path):
    # the table was classified before the list's 'bl' was deleted
    table = tmp_path / 'lc.csv'
    table.write_text(TABLE, encoding='utf-8', newline='\n')
    chars = W.list_characters([(30000.1234, ''), (30000.5, 'c')])
    assert W.observed_wavenumbers(str(table), chars) == [
        (30000.1234, '', 30000.1234), (30000.5, 'c', 30000.5)]
    assert W.observed_wavenumbers(str(table))[0][1] == 'bl'
    import hfs_kappa
    got = hfs_kappa.read_lines(str(table), constants={}, J_of={}, raw=True,
                               chars=chars)
    assert [ln.char for ln in got] == ['']
    assert [ln.char for ln in hfs_kappa.read_lines(
        str(table), constants={}, J_of={}, raw=True)] == ['bl']
