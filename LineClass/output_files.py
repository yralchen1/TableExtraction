"""Can the files a run will write be written, before the run does the work?

Excel holds the workbook or comma-separated file it has open locked against
writing.  A run whose output file is open in Excel therefore fails at the very
last moment - after every level has been read, every candidate weighed and
every cycle turned - with

    PermissionError: [Errno 13] Permission denied: 'line_classifications.xlsx'

and the whole run is lost.  The check here is made at the start instead: it
names the files that are held open and stops before any of the work is done.

Typical use, first thing in a main():

    import output_files
    output_files.require_writable([OUTPUT_FILE, OUTPUT_CSV])

`require_writable` raises SystemExit with a message naming the files when any
of them cannot be written; it returns quietly when they all can.  Paths that
are None or empty are ignored, so an output that this run is not going to
write can be passed unconditionally.

The test is a real one: the existing file is opened for update, which takes
the same write access the eventual write will need and changes nothing.  It
is not a guess from the file name, and it costs nothing measurable.

There remains a gap of the length of the run between the check and the write:
nothing stops the analyst from opening the file in Excel while the run is
going.  The check turns the common case - the file was already open when the
run started - into an immediate, explanatory stop, which is what it is for.

READING is affected by the same lock.  A program that opens a file with the
sharing mode Excel and several editors use denies everyone else even read
access while it holds it, and the hold can last a fraction of a second - long
enough to kill one subprocess of a walk over hundreds of levels and no longer.
:func:`read_retry` is for that: it repeats a read that failed with
PermissionError a few times, a moment apart, and only then gives up.  It is
the read-side companion of `require_writable`, and it is used for the files
that every run reads and some runs write - IDEN2/enlev.dat,
IDEN2/IDEN_level_ids.txt, the Cowan transition list and its cache.
"""
import os
import time

__all__ = ['why_unwritable', 'unwritable', 'require_writable', 'with_twin',
           'read_retry']


def read_retry(read, path, tries: int = 5, pause: float = 0.4, log=None):
    """`read(path)`, repeated while the file is momentarily locked.

    `read` is any callable taking the path - `pd.read_csv`, `open`, a lambda
    that wraps either.  A PermissionError is taken to mean that some other
    program holds the file open right now, and the read is tried again after
    `pause` seconds, then after twice that, and so on, up to `tries` attempts
    in all - about 6 seconds by default.  Any other exception is raised at
    once: a missing file, a bad format or a genuine read-only file is not
    something waiting can cure.

    When the last attempt fails the original PermissionError is raised with
    the file named and the waiting stated, so that the traceback says what
    happened rather than only `[Errno 13] Permission denied`.  `log` is an
    optional printer for the note written when a read has to be repeated;
    the default prints nothing, because most retries succeed on the second
    attempt and are of no interest.
    """
    say = log or (lambda *_: None)
    wait = pause
    for attempt in range(1, tries + 1):
        try:
            return read(path)
        except PermissionError as exc:
            if attempt == tries:
                raise PermissionError(
                    '%s is held open by another program (Excel locks the file '
                    'it has open); %d attempts over %.1f s all failed'
                    % (os.path.basename(path), tries,
                       pause * (2 ** (tries - 1) - 1))) from exc
            say('  %s is locked (open in Excel?) - waiting %.1f s and trying '
                'again (%d of %d)'
                % (os.path.basename(path), wait, attempt, tries - 1))
            time.sleep(wait)
            wait *= 2


def _excel_owner_file(path: str) -> str:
    """The path of Excel's owner file for `path`, whether or not it exists.

    Excel writes a hidden companion named `~$` plus the file name beside a
    workbook it has open, to record who has it.  Its presence is a useful
    detail to print, but not a verdict on its own: Excel leaves it behind when
    it is killed, and a stale one would otherwise block every later run.
    """
    return os.path.join(os.path.dirname(path) or '.',
                        '~$' + os.path.basename(path))


def why_unwritable(path: str):
    """Why `path` cannot be written now, as a phrase, or None when it can.

    An existing file is opened in 'r+b' mode - read and write, no truncation -
    which asks the operating system for exactly the access the write will need
    without touching a single byte.  A file that does not exist yet needs only
    a folder to be created in; a missing folder is reported here rather than
    at the end of the run.
    """
    if not os.path.exists(path):
        folder = os.path.dirname(os.path.abspath(path))
        if not os.path.isdir(folder):
            return 'the folder %s does not exist' % folder
        return None
    if os.path.isdir(path):
        return 'it is a folder, not a file'
    try:
        with open(path, 'r+b'):
            pass
    except PermissionError:
        reason = ('it is open in another program (Excel locks the file it has '
                  'open) or is read-only')
        if os.path.exists(_excel_owner_file(path)):
            reason = ('it is open in Excel (%s, Excel\'s owner file, is beside '
                      'it)' % os.path.basename(_excel_owner_file(path)))
        return reason
    except OSError as exc:
        return str(exc)
    return None


def unwritable(paths) -> list:
    """The (path, reason) pairs of those `paths` that cannot be written now.

    Empty and None paths are skipped, and a path named twice is tested once.
    """
    seen = set()
    bad = []
    for path in paths:
        if not path:
            continue
        key = os.path.normcase(os.path.abspath(path))
        if key in seen:
            continue
        seen.add(key)
        reason = why_unwritable(path)
        if reason is not None:
            bad.append((path, reason))
    return bad


def require_writable(paths, what: str = 'output file') -> None:
    """Stop the run now if any of `paths` cannot be written.

    `what` names the kind of file in the message ('output file', 'ledger',
    ...).  Raises SystemExit, which prints the message and returns a failing
    exit status without a traceback; returns None when all is well.
    """
    bad = unwritable(paths)
    if not bad:
        return
    plural = '' if len(bad) == 1 else 's'
    lines = ['', 'Cannot write the %s%s of this run:' % (what, plural)]
    for path, reason in bad:
        lines.append('    %s - %s' % (os.path.basename(path), reason))
    lines.append('')
    lines.append('Close the file%s and run again.  The run stops here, before '
                 'any of the work is done, rather than after all of it.' % plural)
    raise SystemExit('\n'.join(lines))


def with_twin(csv_path: str) -> list:
    """`csv_path` and the .xlsx twin written beside it.

    The tables saved by chance_mc.save_table come in these pairs, and both
    have to be free for the save to succeed.
    """
    if not csv_path:
        return []
    return [csv_path, os.path.splitext(csv_path)[0] + '.xlsx']
