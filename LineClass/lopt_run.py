"""
lopt_run.py
===========
How the project starts LOPT: the command and the environment it runs in.

LOPT 6 is a Python program with its own interpreter and packages
(F:/LOPT/v6/.venv), started by C:/scripts/lopt.bat, which is on the PATH and,
since lopt.exe was uninstalled (2026-11-18), the only way to start it.  It
must not see the caller's Python path.  A PYTHONPATH that names the project
directory F:/UMD/AI/as-ai - set by an IDE that adds its content roots, and by
insert_new_level.run for its own Python child programs - puts the project's
statistics.py (the Mandel-Paule mean) in front of the standard library's
statistics, and LOPT 6 stops at once on ``from statistics import NormalDist``
with exit status 1 (lopt_perturb, 2026-11-18).  So every LOPT run takes its
environment from env(), which drops the variables that point an interpreter
at another interpreter's modules.
"""

import os

#: the variables that make a Python import modules of somebody else's
PYTHON_PATH_VARS = ('PYTHONPATH', 'PYTHONHOME', 'PYTHONSTARTUP',
                    'PYTHONUSERBASE', 'PYTHONEXECUTABLE')


def env(base=None, **extra):
    """The environment for a LOPT run: `base` (default os.environ) without
    PYTHON_PATH_VARS, with `extra` added."""
    out = dict(os.environ if base is None else base)
    for name in list(out):
        if name.upper() in PYTHON_PATH_VARS:
            del out[name]
    out.update(extra)
    return out


def command(par):
    """The command that runs LOPT on the parameter file `par` (in the
    directory the process is started in).  A .bat file is not an executable
    image, so on Windows it is handed to the command interpreter."""
    return (['cmd', '/c', 'lopt.bat', par] if os.name == 'nt'
            else ['lopt.bat', par])
