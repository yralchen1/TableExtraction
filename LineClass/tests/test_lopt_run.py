import os
import shutil
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import lopt_run as R  # noqa: E402


def test_the_python_path_is_dropped_and_the_rest_kept():
    base = {'PYTHONPATH': 'F:/UMD/AI/as-ai', 'PythonHome': 'x',
            'PYTHONSTARTUP': 'y', 'PATH': 'C:/scripts', 'OTHER': '1'}
    got = R.env(base, OMP_NUM_THREADS='1')
    assert got == {'PATH': 'C:/scripts', 'OTHER': '1', 'OMP_NUM_THREADS': '1'}
    assert 'PYTHONPATH' in base          # the caller's dict is left alone


def test_the_command_names_lopt_bat():
    cmd = R.command('LOPT.par')
    assert cmd[-2:] == ['lopt.bat', 'LOPT.par']


@pytest.mark.skipif(os.name != 'nt' or not shutil.which('lopt.bat'),
                    reason='needs lopt.bat (LOPT 6) on the PATH')
def test_lopt_starts_with_the_project_on_pythonpath(tmp_path):
    # the project root holds a statistics.py that hides the standard
    # library's; with env() LOPT 6 still imports and prints its version
    root = os.path.dirname(os.path.dirname(HERE))
    caller = dict(os.environ, PYTHONPATH=root)
    out = subprocess.run(['cmd', '/c', 'lopt.bat', '--version'], cwd=tmp_path,
                         env=R.env(caller), capture_output=True, text=True)
    assert out.returncode == 0, out.stdout + out.stderr
    assert 'lopt' in out.stdout.lower()
