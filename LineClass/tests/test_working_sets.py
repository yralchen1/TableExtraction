"""Tests for the three working sets: baseline, `iter/`, `final/`.

Run from the LineClass directory:  python -m pytest tests -q

A working set is a directory holding its own configuration file, its own
classification table and its own LOPT files, and sharing with every other set
the files the analyst keeps by hand.  Three things make that work, and they
are what is fixed here:

  * `inherit` in a configuration file, so that a set says only what differs
    from the baseline and the baseline's own relative paths still mean what
    they meant where they were written;
  * `wn_key`, the immutable name of an observed line, so that one decision
    ledger rules on a set whose wavenumbers have been moved by the wavelength
    calibration - the corrections reach 0.66 cm^-1 against the 0.01 cm^-1
    within which a ledger row is matched to its line;
  * `check_sync.py --set`, which looks in the set first and falls back to the
    project directory, which is exactly where the shared files are.
"""
import io
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (ROOT, os.path.dirname(ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

import config                                  # noqa: E402
import classify_lines as C                     # noqa: E402
from models import SpectralLine                # noqa: E402
from swap_paths import working_path            # noqa: E402


BASE = """
[files]
levels     = "../shared/levels.xlsm"
lines      = "lines.xlsx"
icalc      = "Icalc.xlsx"
output     = "out.xlsx"
output_csv = "out.csv"
line_decisions = "line_decisions.csv"

[range]
wn_min = 9327.0
wn_max = 121665.0

[levels.layout]
sheet = "S"
[levels.layout.columns]
id = "ASD_id"

[lines.layout]
sheet = "Sheet1"
[lines.layout.columns]
wn = "own"
u_wn = "unc_own"

[icalc.layout]
sheet = "Icalc"
[icalc.layout.columns]
id1 = "id1"

[icalc.completeness]
gA_cutoff = 1.0e3
allow_below_cutoff = true

[missing_gA]
policy = "impute"
u_ln_window = 10.0
u_ln_estimator = "rms"
self_consistent = false
fit_range_decades = [3.0, 6.0]

[intensity_model]
C = 269.979733
kT = 12233.068631
verify_tolerance = 0.01
"""

CHILD = """
inherit = "../lineclass_config.toml"

[files]
lines      = "lines_corrected.xlsx"
output     = "out.xlsx"
output_csv = "out.csv"

[lines.layout.columns]
wn = "own_corr"
u_wn = "unc_own_corr"
wn_key = "own"
"""


def write(path, text):
    d = os.path.dirname(path)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with io.open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(text)
    return path


@pytest.fixture
def two_sets(tmp_path):
    """A baseline configuration and a set that inherits it."""
    root = str(tmp_path / 'LineClass')
    write(os.path.join(root, 'lineclass_config.toml'), BASE)
    write(os.path.join(root, 'iter', 'lineclass_config.toml'), CHILD)
    return root


# ---------------------------------------------------------------------------
# inherit
# ---------------------------------------------------------------------------
def test_the_set_takes_what_it_does_not_override(two_sets):
    """The intensity model is written once and both sets use it."""
    base = config.load(os.path.join(two_sets, 'lineclass_config.toml'))
    iter_ = config.load(os.path.join(two_sets, 'iter',
                                     'lineclass_config.toml'))
    assert iter_.intensity_model == base.intensity_model
    assert (iter_.wn_min, iter_.wn_max) == (base.wn_min, base.wn_max)


def test_an_inherited_path_keeps_its_own_directory(two_sets):
    """`../shared/levels.xlsm` names the same workbook for both sets.

    Resolving it against the inheriting file's directory instead would put it
    beside the set, which is the failure this arrangement exists to avoid.
    """
    base = config.load(os.path.join(two_sets, 'lineclass_config.toml'))
    iter_ = config.load(os.path.join(two_sets, 'iter',
                                     'lineclass_config.toml'))
    assert iter_.levels_file == base.levels_file
    assert os.path.basename(os.path.dirname(iter_.levels_file)) == 'shared'


def test_the_hand_kept_files_are_shared_and_the_outputs_are_not(two_sets):
    base = config.load(os.path.join(two_sets, 'lineclass_config.toml'))
    iter_ = config.load(os.path.join(two_sets, 'iter',
                                     'lineclass_config.toml'))
    # one ledger, ruling on both sets
    assert iter_.line_decisions == base.line_decisions
    # its own table, so a run here cannot desync the baseline
    assert iter_.output_csv != base.output_csv
    assert os.path.basename(os.path.dirname(iter_.output_csv)) == 'iter'
    # its own line list
    assert os.path.basename(iter_.lines_file) == 'lines_corrected.xlsx'


def test_a_table_is_merged_key_by_key(two_sets):
    """The set names two columns; the other four come from the baseline."""
    iter_ = config.load(os.path.join(two_sets, 'iter',
                                     'lineclass_config.toml'))
    assert iter_.lines.columns == {'wn': 'own_corr', 'u_wn': 'unc_own_corr',
                                   'wn_key': 'own'}
    assert iter_.lines.sheet == 'Sheet1'


def test_a_configuration_that_inherits_itself_stops_the_run(tmp_path):
    p = write(str(tmp_path / 'a.toml'), 'inherit = "a.toml"\n')
    with pytest.raises(config.ConfigError) as err:
        config.load(p)
    assert 'inherits itself' in str(err.value)


def test_the_environment_variable_chooses_the_set(two_sets, monkeypatch):
    """A dozen programs import classify_lines, which loads the configuration
    while it is being imported; only something read before the import can
    choose the set for all of them."""
    import importlib
    monkeypatch.setenv('LINECLASS_CONFIG',
                       os.path.join(two_sets, 'iter', 'lineclass_config.toml'))
    reloaded = importlib.reload(config)
    try:
        assert os.path.basename(
            os.path.dirname(reloaded.DEFAULT_PATH)) == 'iter'
        assert os.path.basename(reloaded.load().lines_file) == \
            'lines_corrected.xlsx'
    finally:
        monkeypatch.delenv('LINECLASS_CONFIG')
        importlib.reload(config)


# ---------------------------------------------------------------------------
# wn_key
# ---------------------------------------------------------------------------
def line(wn, key=None):
    return SpectralLine(wavenumber=wn, wn_uncertainty=0.05, intensity=10.0,
                        line_character='', wn_key=key if key is not None else wn)


LEDGER = ('wn_obs,low_id,upp_id,decision,reason\n'
          '45994.3206,059003.000141,059003.000218,accept,off the block trend\n')


def test_the_ledger_finds_its_line_on_the_published_scale(tmp_path):
    """The baseline: name and wavenumber are the same column."""
    path = write(str(tmp_path / 'line_decisions.csv'), LEDGER)
    lines = [line(45994.3206), line(46172.4502)]
    sites = C.attach_line_decisions(lines, path)
    assert len(sites) == 1
    assert lines[0].decisions[('059003.000141', '059003.000218')][0] == 'accept'


def test_the_ledger_finds_its_line_after_a_calibration_shift(tmp_path):
    """0.66 cm^-1 is sixty-six times DECISIONS_WN_MATCH, so a ledger keyed on
    the wavenumbers this set works with would miss every line it names."""
    path = write(str(tmp_path / 'line_decisions.csv'), LEDGER)
    shifted = line(45994.3206 + 0.66, key=45994.3206)
    sites = C.attach_line_decisions([shifted, line(46172.4502)], path)
    assert len(sites) == 1
    assert shifted.decisions[('059003.000141', '059003.000218')][0] == 'accept'
    # the key written in the ledger is what is reported back, not the shifted
    # wavenumber the set happens to be working with
    assert sites[0][0][0] == pytest.approx(45994.3206)


def test_a_wavenumber_that_is_no_line_still_stops_the_run(tmp_path):
    path = write(str(tmp_path / 'line_decisions.csv'), LEDGER)
    with pytest.raises(ValueError) as err:
        C.attach_line_decisions([line(46172.4502)], path)
    assert '45994.3206' in str(err.value)


def test_a_line_names_itself_where_the_set_does_not_say_otherwise():
    """`wn_key` defaults to the `wn` column, so nothing about a baseline run
    changes: `col.get('wn_key', col['wn'])` in read_observed_lines()."""
    cols = {'wn': 3, 'u_wn': 4}
    assert cols.get('wn_key', cols['wn']) == 3


# ---------------------------------------------------------------------------
# check_sync --set
# ---------------------------------------------------------------------------
def test_the_set_gets_its_own_files_and_shares_the_rest(tmp_path,
                                                         monkeypatch):
    """`working_path(name, cwd=set)` is the whole of `--set`: the set's own
    LOPT files and classification table, the project's ledger."""
    import swap_paths
    project = str(tmp_path / 'LineClass')
    set_dir = os.path.join(project, 'iter')
    write(os.path.join(project, 'line_decisions.csv'), 'wn_obs\n')
    write(os.path.join(project, 'line_classifications.csv'), 'wn_obs\n')
    write(os.path.join(set_dir, 'line_classifications.csv'), 'wn_obs\n')
    monkeypatch.setattr(swap_paths, 'PROJECT', project)

    assert working_path('line_classifications.csv', cwd=set_dir) == \
        os.path.join(set_dir, 'line_classifications.csv')
    assert working_path('line_decisions.csv', cwd=set_dir) == \
        os.path.join(project, 'line_decisions.csv')


# ---------------------------------------------------------------------------
# the corrections are applied once
# ---------------------------------------------------------------------------
def test_correcting_an_already_corrected_table_stops_the_run():
    """`--corrections` on the baseline table and classifying the corrected set
    are alternatives; doing both would shift every wavenumber twice."""
    import make_LOPT_input as M
    rows = [{'wn_obs': '45994.9806', 'wn_key': '45994.3206',
             'unc_wn_obs': '0.12'}]
    with pytest.raises(SystemExit) as err:
        M.apply_corrections(rows, {'45994.9806': (0.66, 0.02)})
    assert 'twice' in str(err.value)


def test_the_baseline_table_is_corrected_as_before():
    """Its wn_key equals its wn_obs, so nothing is in the way."""
    import make_LOPT_input as M
    rows = [{'wn_obs': '45994.3206', 'wn_key': '45994.3206',
             'unc_wn_obs': '0.12'}]
    shifted, largest, unknown = M.apply_corrections(
        rows, {'45994.3206': (0.66, 0.02)})
    assert (shifted, unknown) == (1, [])
    assert largest == pytest.approx(0.66)
    assert float(rows[0]['wn_obs']) == pytest.approx(45994.9806)


def test_a_table_without_the_column_is_corrected_as_before():
    """Tables written before wn_key existed still work."""
    import make_LOPT_input as M
    rows = [{'wn_obs': '45994.3206', 'unc_wn_obs': '0.12'}]
    shifted, _, _ = M.apply_corrections(rows, {'45994.3206': (0.66, 0.02)})
    assert shifted == 1
