"""Tests for the writing of the fitted kappas into the configuration
(config.kappa_file, config.set_hfs_kappa), which wavelength_calibration.py
does at the end of every run that writes its files.

Run from the LineClass directory:  python -m pytest tests -q

The kappas go into the file of the configuration chain that holds the
[hfs.kappa] table - the baseline's, which a working set inherits - and only
the lines of the three fitted classes and the program's own note change.
"""
import os
import sys
import tomllib

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import config                              # noqa: E402
import wavelength_calibration as W         # noqa: E402

PARENT = """\
# the baseline
[hfs]
apply = false

[hfs.kappa]
# a hand comment that must survive
#   flag   the anchor
flag       = [1.0, 0.0]
plain_1974 = [0.926, 0.012]
plain_1969 = [0.717, 0.021]
c          = [0.276, 0.044]

[other]
x = 1
"""

CHILD = """\
inherit = "../parent.toml"

[hfs]
apply = true
"""

NEW = {'plain_1974': (0.85, 0.0131), 'plain_1969': (0.6514, 0.021),
       'c': (0.172, 0.051)}


@pytest.fixture
def chain(tmp_path):
    parent = tmp_path / 'parent.toml'
    parent.write_text(PARENT, encoding='utf-8', newline='\n')
    (tmp_path / 'set').mkdir()
    child = tmp_path / 'set' / 'child.toml'
    child.write_text(CHILD, encoding='utf-8', newline='\n')
    return str(parent), str(child)


def kappa_of(path):
    with open(path, 'rb') as fh:
        return tomllib.load(fh)['hfs']['kappa']


def test_the_file_that_holds_the_table(chain, tmp_path):
    parent, child = chain
    assert config.kappa_file(child) == os.path.abspath(parent)
    assert config.kappa_file(parent) == os.path.abspath(parent)
    # a set that writes its own table is its own file
    own = tmp_path / 'own.toml'
    own.write_text(CHILD.replace('apply = true', 'apply = true\n\n'
                                 '[hfs.kappa]\nflag = [1.0, 0.0]\n'),
                   encoding='utf-8', newline='\n')
    assert config.kappa_file(str(own)) == str(own)
    lone = tmp_path / 'lone.toml'
    lone.write_text('[hfs]\napply = false\n', encoding='utf-8')
    with pytest.raises(config.ConfigError, match='no file of the chain'):
        config.kappa_file(str(lone))


def test_the_values_are_written_and_nothing_else(chain):
    parent, child = chain
    written, before = config.set_hfs_kappa(child, NEW, 'line one\nline two')
    assert written == os.path.abspath(parent)
    assert before['plain_1974'] == (0.926, 0.012)
    got = kappa_of(parent)
    assert got['plain_1974'] == [0.85, 0.013]
    assert got['plain_1969'] == [0.651, 0.021]
    assert got['c'] == [0.172, 0.051]
    assert got['flag'] == [1.0, 0.0]
    text = open(parent, encoding='utf-8', newline='').read()
    assert '\r' not in text
    assert '# a hand comment that must survive' in text
    assert 'plain_1974 = [0.850, 0.013]' in text
    assert 'c          = [0.172, 0.051]' in text
    assert '[other]\nx = 1' in text
    # the note sits just above flag
    lines = text.split('\n')
    i = lines.index('flag       = [1.0, 0.0]')
    assert lines[i - 2:i] == [config.KAPPA_NOTE + 'line one',
                              config.KAPPA_NOTE + 'line two']
    # the child is not touched
    assert open(child, encoding='utf-8').read() == CHILD


def test_a_second_run_replaces_the_note(chain):
    parent, child = chain
    config.set_hfs_kappa(child, NEW, 'first')
    config.set_hfs_kappa(child, {'plain_1974': (0.9, 0.01),
                                 'plain_1969': (0.7, 0.02),
                                 'c': (0.2, 0.05)}, 'second')
    text = open(parent, encoding='utf-8').read()
    assert config.KAPPA_NOTE + 'first' not in text
    assert text.count(config.KAPPA_NOTE) == 1
    assert kappa_of(parent)['plain_1974'] == [0.9, 0.01]


def test_crlf_is_kept(chain):
    parent, child = chain
    with open(parent, 'w', encoding='utf-8', newline='') as fh:
        fh.write(PARENT.replace('\n', '\r\n'))
    config.set_hfs_kappa(child, NEW, 'note')
    text = open(parent, encoding='utf-8', newline='').read()
    assert '\r\n' in text and '\n' not in text.replace('\r\n', '')
    assert kappa_of(parent)['c'] == [0.172, 0.051]


def test_a_table_it_cannot_write_is_left_as_it_was(chain):
    parent, child = chain
    broken = PARENT.replace('c          = [0.276, 0.044]\n', '')
    with open(parent, 'w', encoding='utf-8', newline='') as fh:
        fh.write(broken)
    with pytest.raises(config.ConfigError, match='lacks the line of c'):
        config.set_hfs_kappa(child, NEW, 'note')
    assert open(parent, encoding='utf-8', newline='').read() == broken
    with pytest.raises(config.ConfigError, match='no value for c'):
        config.set_hfs_kappa(child, {'plain_1974': (0.9, 0.01),
                                     'plain_1969': (0.7, 0.02)})


def test_the_note_of_the_calibration(chain):
    parent, _ = chain
    note = W.kappa_note(NEW, W.before_of(parent),
                        os.path.join(W.HERE, 'iter_hfs',
                                     'line_classifications.csv'), 4494)
    first, second = note.split('\n')
    assert first.endswith('from iter_hfs/line_classifications.csv '
                          '(4494 lines);')
    assert second == ('the values before: plain_1974 0.926(12), '
                      'plain_1969 0.717(21), c 0.276(44)')


def test_the_baseline_holds_the_table_iter_hfs_uses():
    path = os.path.join(ROOT, 'iter_hfs', 'lineclass_config.toml')
    assert config.kappa_file(path) == os.path.join(ROOT,
                                                   'lineclass_config.toml')
