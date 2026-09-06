"""Where the swap tools look for the files they rewrite.

The four ``swap_line_assignments*`` scripts live in the project directory
next to the working copies of the LOPT files, of ``IDEN2/`` and of the
pipeline's overlay files.  That is the ordinary place to work, but not the
only one: a trial exchange is usually made in an iteration folder
(``iter22/`` and the like) holding its own copy of the LOPT parameter file,
its own transitions file and its own output, so that the run can be looked
at and thrown away without disturbing the working copies.

So a default file name is resolved twice:

  * first against the directory the command was run from - if the file is
    there, that is the one meant, because standing in ``iter22`` and asking
    for the transitions file can only mean ``iter22``'s;
  * then against the directory holding the scripts, which is the project
    directory.  This is what makes the tools work from anywhere for the
    files that exist in only one place - the decision ledger, the level
    overrides, the published level list - and it is what happens when the
    command is run from the project directory itself, where the two answers
    coincide.

Nothing is resolved at import time: the working directory may change between
importing a module and parsing its arguments, and the tests rely on that.
Each script therefore calls ``working_path`` while building its argument
parser, and prints every path it settled on before it writes anything, so
that a file picked up from the wrong one of the two directories is visible
in the report rather than only in the result.
"""

import os

PROJECT = os.path.dirname(os.path.abspath(__file__))


def working_path(name, cwd=None):
    """``name`` in the current directory if it is there, else in PROJECT.

    Works for a directory (``IDEN2``) exactly as for a file.  A name that
    exists in neither place comes back as the project's, so that the error
    message a caller then produces names the file it actually expected.
    """
    here = os.path.abspath(cwd if cwd is not None else os.getcwd())
    candidate = os.path.join(here, name)
    if os.path.exists(candidate):
        return candidate
    return os.path.join(PROJECT, name)


def describe(name, path):
    """One line of the report: which of the two directories a path came from."""
    where = os.path.dirname(os.path.abspath(path))
    if where == os.path.abspath(os.getcwd()):
        return '%s (the current directory)' % path
    if where == PROJECT:
        return '%s (the project directory)' % path
    return path
