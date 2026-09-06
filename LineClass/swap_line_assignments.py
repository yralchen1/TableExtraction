"""
swap_line_assignments.py
========================
Exchange two levels' measurements everywhere at once: in the LOPT working
files, in the IDEN2 working files, and in the overlay the classification
pipeline reads.

What an exchange is
-------------------
Two levels of the same parity and the same J, close enough in energy that
the calculation cannot tell them apart, can have had their **theoretical
identities interchanged**: both measured energies are right, but each was
attached to the other's calculated level - the other configuration, the
other term, and with them the other whole set of calculated transition
probabilities.  ``level_interchange.py`` is what looks for this.

The repair is not to rename the two levels.  A level identifier such as
``059003.000483`` appears in files in several folders, and IDEN2's own level
numbering is fixed by the *calculated* energy and cannot move at all.  So
the identifiers stay exactly where they are and everything the
**measurement** gave them is exchanged instead: their observed lines, and
hence their energies.  Each identifier keeps its own calculated identity and
takes the other's measured position.

Why it takes three scripts
--------------------------
Because the exchange has to be made in three unrelated places, and none of
them can see the other two.

``swap_line_assignments_LOPT.py``
    The LOPT transitions file, where each observed line names its two
    levels, and the fixed-levels file.  The next LOPT run then returns the
    two energies exchanged.

``swap_line_assignments_IDEN.py``
    IDEN2's ``enlev.dat`` and ``trans.dat``, so that IDEN2 shows each
    calculated level against the observed lines that now belong to it and
    the branching patterns can be judged by eye.

``swap_line_assignments_pipeline.py``
    ``revised_level_energies.csv`` and ``line_decisions.csv`` - the two
    small files under the analyst's hand that the classification reads.
    Without them the exchange is only in working files: ``classify_lines.py``
    would put both levels back at their published energies on the next run,
    and ``make_LOPT_input.py`` would then rebuild the LOPT transitions file
    from that.

Nothing in ``Icalc.xlsx``, in ``icalc_new.xlsx`` or in the published line
workbook changes.  The calculated transitions are keyed by the pair of level
identifiers, and an identifier keeps its calculated identity through the
exchange; it is the measurement that moves.

The order of work
-----------------
This script runs the three in the order above and stops at the first one
that refuses, so that a pair the LOPT step rejects - a typo in an
identifier, a line joining the two levels, a fixed energy on one of them
only - never reaches the other two.  Nothing is written until every check
of that step has passed, and ``--dry-run`` runs all three as reports.

Between the first step and the third it also tries to run LOPT itself: the
command ``lopt`` (``--lopt-command``), in the directory the LOPT files it has
just rewritten were found in.  That refit is not needed to record the
exchange - an exchange writes each identifier the energy the other already
had, and the third step can do that arithmetic on its own - but it is what
turns those two numbers into least-squares energies refitted with every other
level, which differ from the plain exchange by a few thousandths of a
wavenumber, and it is also the check that the rewritten transitions file is
one LOPT will still read.

If the command is not there, or returns an error, that is reported as a
warning and the run continues.  The third step then sees a level table that
does not yet know about the exchange, says so, and exchanges the two energies
itself.  Running LOPT by hand afterwards and then

    python swap_line_assignments.py id1 id2 --after-lopt

replaces those two energies with the fitted ones and leaves everything else
it wrote alone; the same command is worth running after any later LOPT run
that has moved the two levels.

Usage
-----
    python swap_line_assignments.py 059003.000483 059003.000398 --dry-run
    python swap_line_assignments.py 059003.000483 059003.000398
    python swap_line_assignments.py id1 id2 --index 721 723
    python swap_line_assignments.py id1 id2 --only pipeline
    python swap_line_assignments.py id1 id2 --no-lopt-run
    python swap_line_assignments.py id1 id2 --after-lopt

Both arguments are experimental level identifiers: the LOPT files and the
pipeline know levels by nothing else.  ``--index`` gives IDEN2's two row
numbers when they are known, which saves that step a lookup.

Every file name each step needs is looked for first in the directory the
command was run from and then in the directory holding these scripts, so
that working in an iteration folder acts on that folder's copies while the
files that exist in only one place - the decision ledger, the level
overrides - are still found.
"""

import argparse
import os
import subprocess
import sys

import swap_line_assignments_IDEN as iden
import swap_line_assignments_LOPT as lopt
import swap_line_assignments_pipeline as pipeline
import swap_paths

STEPS = ('lopt', 'iden', 'pipeline')

TITLE = {
    'lopt': 'The LOPT working files',
    'iden': 'The IDEN2 working files',
    'pipeline': 'The files the classification pipeline reads',
}


def parse_args(argv):
    p = argparse.ArgumentParser(
        description='Exchange two levels\' observed lines and observed '
                    'energies in the LOPT files, the IDEN2 files and the '
                    'pipeline overlay, keeping both identifiers where they '
                    'are.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('id1', help='identifier of the first level')
    p.add_argument('id2', help='identifier of the second level')
    p.add_argument('--only', default=','.join(STEPS),
                   help='which steps to run, comma-separated: '
                        + ', '.join(STEPS))
    p.add_argument('--after-lopt', action='store_true',
                   help='the exchange has been made and LOPT has been run '
                        'again: run only the pipeline step, and take the two '
                        'energies from the level table as they stand instead '
                        'of exchanging them')
    p.add_argument('--lopt-command', default='lopt',
                   help='the command that runs LOPT, tried after the LOPT '
                        'files have been rewritten and before the pipeline '
                        'step reads the new level table')
    p.add_argument('--no-lopt-run', action='store_true',
                   help='do not try to run LOPT; the pipeline step then '
                        'exchanges the two energies itself')
    p.add_argument('--index', nargs=2, type=int, metavar=('N1', 'N2'),
                   help='the two levels\' row numbers in enlev.dat, passed '
                        'to the IDEN2 step so that it need not look them up')
    p.add_argument('--rejects-only', action='store_true',
                   help='passed to the pipeline step: for an identification '
                        'of the published line list write only the reject at '
                        'the level it names')
    p.add_argument('--date', help='passed to the pipeline step: the date '
                                  'written into the rows it adds')
    p.add_argument('--force', action='store_true',
                   help='passed to the IDEN2 step: exchange two levels of '
                        'different J, which is normally refused')
    p.add_argument('--no-backup', action='store_true',
                   help='do not keep the previous contents of any file in '
                        '<file>.bak')
    p.add_argument('--dry-run', action='store_true',
                   help='run every step as a report and write nothing')
    return p.parse_args(argv)


def wanted(args):
    """The steps to run, in the order they have to be run in."""
    if args.after_lopt:
        return ['pipeline']
    names = [s.strip().lower() for s in args.only.split(',') if s.strip()]
    unknown = [n for n in names if n not in STEPS]
    if unknown:
        raise ValueError('--only: no step called %s; the steps are %s'
                         % (', '.join(unknown), ', '.join(STEPS)))
    if not names:
        raise ValueError('--only: no step named')
    return [s for s in STEPS if s in names]


def lopt_directory():
    """Where to run LOPT: the directory its parameter file was found in.

    LOPT reads ``LOPT.par`` from the directory it is started in, and that
    file names the transitions file the first step has just rewritten, so
    running it anywhere else would refit a different set of lines.  The
    parameter file is looked for by the same rule as everything else - the
    current directory first, then the directory holding these scripts.
    """
    return os.path.dirname(swap_paths.working_path('LOPT.par'))


def run_lopt(command, directory):
    """Run LOPT.  Returns None if it worked, or a sentence saying why not.

    Nothing here is fatal.  The refit only sharpens the two energies by a few
    thousandths of a wavenumber; without it the pipeline step exchanges the
    two energies of the old level table itself, which is the same exchange.
    So a missing command or a non-zero exit is reported and the run goes on.
    """
    try:
        # shell=True so that a .bat or a shell function is found the way it
        # would be if the command were typed.
        code = subprocess.call(command, shell=True, cwd=directory)
    except OSError as exc:
        return 'it could not be started (%s)' % exc
    if code != 0:
        return 'it stopped with status %d' % code
    return None


def step_argv(step, args):
    """The command line one step is given."""
    argv = [args.id1, args.id2]
    if args.dry_run:
        argv.append('--dry-run')
    if args.no_backup:
        argv.append('--no-backup')
    if step == 'iden':
        if args.index:
            argv += ['--index', str(args.index[0]), str(args.index[1])]
        if args.force:
            argv.append('--force')
    elif step == 'pipeline':
        if args.rejects_only:
            argv.append('--rejects-only')
        if args.date:
            argv += ['--date', args.date]
    return argv


def lopt_step(args):
    """Refit the levels between the first step and the third, if we can."""
    print('')
    print('=' * 74)
    print('The refit  (%s)' % args.lopt_command)
    print('=' * 74)
    if args.dry_run:
        print('--dry-run: LOPT is not run.')
        return
    if args.no_lopt_run:
        print('--no-lopt-run: LOPT is not run.  The step below will exchange '
              'the two energies of the level table as it stands.')
        return
    where = lopt_directory()
    print('Running %s in %s' % (args.lopt_command, where))
    trouble = run_lopt(args.lopt_command, where)
    if trouble is None:
        print('')
        print('LOPT finished.  The level table below is the refitted one.')
        return
    print('')
    print('WARNING: %s was not run - %s.' % (args.lopt_command, trouble))
    print('LOPT_output_levels.txt therefore still holds the energies as they')
    print('were before the exchange, and the step below will exchange those')
    print('two energies itself.  That is the same exchange to within the few')
    print('thousandths of a wavenumber the refit would move them.  Run LOPT')
    print('by hand and then:')
    print('    python swap_line_assignments.py %s %s --after-lopt'
          % (args.id1, args.id2))


def main(argv=None):
    args = parse_args(argv)
    steps = wanted(args)
    modules = {'lopt': lopt, 'iden': iden, 'pipeline': pipeline}

    print('swap_line_assignments.py')
    print('  exchanging       : %s  <->  %s' % (args.id1, args.id2))
    print('  steps            : %s' % ', '.join(steps))
    if args.dry_run:
        print('  --dry-run        : every step reports and writes nothing')

    for step in steps:
        if step == 'pipeline' and 'lopt' in steps:
            lopt_step(args)
        print('')
        print('=' * 74)
        print('%s  (%s)' % (TITLE[step], modules[step].__name__ + '.py'))
        print('=' * 74)
        name = modules[step].__name__ + '.py'
        try:
            code = modules[step].main(step_argv(step, args))
        except ValueError as exc:
            # Each step raises rather than returning when a check of its own
            # fails.  Say which step it was: the message names files, and
            # three scripts read overlapping sets of them.
            sys.stderr.write('\n%s: %s\n' % (name, exc))
            sys.stderr.write('The steps after it were not run.\n')
            return 1
        if code:
            print('')
            print('%s stopped with status %d; the steps after it were not '
                  'run.' % (name, code))
            return code

    print('')
    print('=' * 74)
    if args.dry_run:
        print('Nothing was written.  Run again without --dry-run to make the '
              'exchange.')
        return 0
    if 'lopt' in steps:
        print('If LOPT was not run above, run it (lopt.bat) and then')
        print('    python swap_line_assignments.py %s %s --after-lopt'
              % (args.id1, args.id2))
        print('to replace the exchanged energies in the overrides with the '
              'fitted ones.  Everything else already written is recognised '
              'and left alone.')
    if 'pipeline' in steps:
        print('Then run classify_lines.py, and check in '
              'lopt_vs_classify_levels.csv that each of the two levels came '
              'out at the energy that was written for it.  A level fitted '
              'between two groups of lines has kept an identification at its '
              'old position.')
    if 'iden' in steps:
        print('Open both levels in IDEN2 and look at the branching patterns '
              'before accepting the exchange.')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except ValueError as exc:
        sys.stderr.write('\nswap_line_assignments.py: %s\n' % exc)
        sys.exit(1)
