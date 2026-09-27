"""Bring IDEN2's ``enlev.dat``, ``trans.dat`` and ``dlv.dat`` up to date.

WHAT IS OUT OF DATE, AND WHY
============================
IDEN2 keeps three files of its own.  ``enlev.dat`` is the level list it shows
on screen: one row per calculated level, carrying the level's calculated
energy, the measured energy adopted for it, the uncertainty of that measurement
and an asterisk when the level has been found at all.  ``trans.dat`` is the
predicted transition list: one block per upper level, one row per transition,
carrying a copy of the partner's energy, the predicted wavenumber, a code for
the predicted intensity, and - where the transition has been identified with an
observed line - that line and its departure from the prediction.  ``dlv.dat``
is the observed line list: one row per measured line, carrying its intensity,
its wavenumber, its standard wavelength, its character and the uncertainty of
the measurement.

All three drift out of step with the pipeline, in five separate ways.

*The level energies.*  Every LOPT run moves the measured levels a little.
407 of the 594 levels differ from ``LOPT_output_levels.txt`` as things stand,
by up to 0.19 cm^-1.  Every predicted wavenumber in ``trans.dat`` is a
difference of two of those energies, so all of them inherit the drift.

*The predicted intensities.*  ``trans.dat`` stores the intensity as

    Icalc_IDEN2 = round(10 * ln(Icalc)) ,

a small signed integer - the file's whole range is -37 to 146.  Those codes
were computed years ago, from an excitation temperature of about 13100 cm^-1
and a normalisation that suited the intensity scale of the time.  The
pipeline now works at kT = 12233 cm^-1 with a different C, and the two scales
differ by a factor 0.70 and a Boltzmann tilt.  The ranking the old codes give
is still very nearly right - the rank correlation with the current Icalc is
0.9994 - but the numbers themselves are not the numbers the rest of the
pipeline uses.

*The transition list itself.*  ``trans.dat`` holds 104272 transitions.  The
Cowan table it was built from, ``tp_E1_no_trials.xlsx``, holds 120273 whose
two levels both appear in ``enlev.dat``; the rest were cut off at the bottom
of the old intensity scale, at code -37.  Which transitions that cut keeps is
not the same on the new scale as on the old one.

*The observed lines.*  Every wavenumber in ``dlv.dat`` and every uncertainty
beside it is a measurement, and both are revised: a corrected set applies the
wavelength calibration correction to every line, and the uncertainty model is
refitted whenever the evidence for it changes.  A line whose uncertainty is
wrong here is a line the eye judges against the wrong window, which is the
whole purpose the file serves.

*Which lines are identified.*  The pipeline's statement of an identification is
a record of ``LOPT_input_lines.txt``: a record not flagged ``P`` is one the fit
uses, and a transition all of whose records are flagged ``P`` is one that was
considered and rejected.  Every classification run changes that set, and
``trans.dat`` goes on showing the identifications of the run before it.  A
withdrawn identification left on the screen is worse than no identification at
all, because it also consumes a line that some other transition may want.

WHAT THIS PROGRAM DOES
======================
It rewrites both files from the current fit:

 1. ``LOPT_output_levels.txt`` gives the measured energy of every found level
    and its uncertainty, which is the larger of LOPT's ``D1`` and ``D2tot``
    rounded to a thousandth of a wavenumber - the rule the file already
    follows.  These go into ``enlev.dat``, with the observed-minus-calculated
    column recomputed from each level's own calculated energy.  A level that
    has not been found keeps its calculated energy as its adopted energy -
    nothing in this run knows better - but its uncertainty is rewritten too,
    and it means something quite different there: how far from its calculated
    position such a level is likely to be found.  That is measured for each
    configuration separately, as the rms of E_obs - E_calc over the levels of
    that configuration that HAVE been found, because the calculation is wrong
    by different amounts in different parts of the spectrum.  Both kinds of
    change are listed in the report.

 2. ``tp_E1_no_trials.xlsx`` gives gA for every calculated transition.  With
    the C and kT of ``[intensity_model]`` in ``lineclass_config.toml`` - the
    same constants ``Icalc.xlsx`` obeys, checked against it before anything is
    written - each becomes

        Icalc = C * gA * (rwn/1e8) * exp(-Eup/kT)

    with rwn the predicted wavenumber and Eup the upper level's adopted
    energy, and then the code round(10*ln(Icalc)).  For a transition between
    two found levels this reproduces the ``Icalc`` column of ``Icalc.xlsx``
    to the last digit, because that column is this same formula applied to
    this same gA - except where the two disagree about where a level is.
    ``Icalc.xlsx`` carries the wavenumbers and upper-level energies of the
    level list it was written against, and some levels have been revised by
    several hundred wavenumbers since; on 108 transitions as things stand the
    code comes out one unit, a tenth of a natural logarithm, away from what
    that file would give.  The energies used here are the current ones.

 3. Every row of ``trans.dat`` is rewritten with the partner's new energy, the
    new predicted wavenumber and the new intensity code, and every
    identification is recomputed against the new prediction.

 4. The set's own line list - ``[files] lines`` of its configuration, read on
    the ``wn`` and ``u_wn`` columns its ``[lines.layout.columns]`` names -
    gives the wavenumber and the uncertainty of every observed line, and those
    go into ``dlv.dat``.  Its row numbers are how ``trans.dat`` names a line,
    so no row is added, removed or reordered: only the three measured fields
    of each row change.  The uncertainty there is a WAVELENGTH uncertainty in
    angstroms, as ``numset.dat`` states it too, and it is converted.

 5. ``LOPT_input_lines.txt`` says which transitions the fit is given, and the
    identifications in ``trans.dat`` are made to match it: an accepted
    transition keeps or gets its line, a transition every record of which is
    flagged ``P`` loses it, and so does one the file does not mention at all.
    Each removal and each addition is named in the report.

    ``--keep-unlisted`` suspends the third of those.  It is for the one case
    where the pipeline is the party out of date rather than IDEN2: lines just
    marked by hand on the screen, which ``classify_lines.py`` has not yet been
    run over.  ``--no-lines`` suspends steps 4 and 5 altogether.

WHAT IS NEVER LOST
==================
**An identification the pipeline still holds is never dropped.**  Such a
transition stays in ``trans.dat`` however weak its new intensity code makes it,
and whether or not the Cowan table still has a gA for it.  One transition has
no gA at all: it is Sugar's identification of a line that fell below the
gA = 1e3 s^-1 floor of the Cowan run and so was never calculated, adopted
because the observed line's ``*v`` character matches the large hyperfine
structure of the level, as the ``*r``, ``cl``, ``c`` and ``w`` characters of
the transitions upward from it do.  It was added to ``trans.dat`` by hand with
an intensity code of 0, it is the one imputed ``Icalc`` in
``line_classifications.csv``, and this program carries its row through
untouched but for the geometry.

What is dropped is an identification the pipeline has withdrawn, and the report
lists every one of them with its wavenumber, its two levels and its row in
``dlv.dat``, so that a removal can be looked at and, if it is wrong, put back
by hand.  The previous files are kept under ``--backup-suffix`` in any case.

WHICH TRANSITIONS ARE LISTED
============================
A transition is written when it carries an identified line, or when its new
intensity code reaches ``--cutoff``.  The default cutoff is -41, which is the
value that leaves the present list most nearly alone: it adds 420 transitions
that the old scale had cut off and drops 47 - 38 that no longer reach the
floor and 9 that have gone from the Cowan table - out of 104272.  A larger
cutoff makes a shorter file - each step of one removes about 700 transitions,
being a factor 1.105 in intensity - and
``--cutoff -1000`` writes every calculated transition there is.

USAGE
=====
    python sync_IDEN2.py                 # rewrite both files, keeping backups
    python sync_IDEN2.py --dry-run       # say what would change, write nothing
    python sync_IDEN2.py --cutoff -45    # a longer list, down to weaker lines
    python sync_IDEN2.py --report sync_IDEN2_report.txt
    python ../sync_IDEN2.py              # run from iter/: syncs iter/IDEN2
    python sync_IDEN2.py --set iter      # the same, named from anywhere
    python sync_IDEN2.py --no-lines      # levels and predictions only
    python sync_IDEN2.py --keep-unlisted # keep identifications made by hand

THE WORKING SET
===============
The files rewritten are the ones of the set the command is run from, not the
ones beside the script.  ``IDEN2``, ``LOPT_output_levels.txt``,
``LOPT_input_lines.txt`` and ``lineclass_config.toml`` are looked for in that
directory first and in the project directory only if the set has not got them,
which is what ``swap_paths`` describes and what ``check_sync.py --set`` already
does.  An IDEN2 from one set and a level or transitions
table from another are refused: the corrected set's levels are on the
corrected wavenumber scale and the baseline's are on Sugar's, so writing one
into the other would move every Ritz wavenumber silently.  Both resolved
paths are printed before anything is written.

Close IDEN2 first.  It holds its files open and rewrites them from memory when
it exits, which would undo everything this program has done.
"""

import argparse
import bisect
import io
import math
import os
import shutil
import sys

import numpy as np
import pandas as pd

import config
import cowan_gA
import gA_imputation
import level_interchange
import output_files
import swap_line_assignments_IDEN as IDEN
from swap_paths import working_path

HERE = os.path.dirname(os.path.abspath(__file__))
IDEN2_DIR = os.path.join(HERE, 'IDEN2')
LOPT_LEVELS = os.path.join(HERE, 'LOPT_output_levels.txt')
ICALC_FILE = os.path.join(HERE, 'Icalc.xlsx')
DEF_CUTOFF = -41
# An rms taken over one or two levels is not a measure of anything, so a
# configuration with fewer found levels than this is left alone.
MIN_CFG_LEVELS = 3
BACKUP_SUFFIX = '.presync'
# dlv.dat, one row per observed line, fixed columns.  The uncertainty is a
# WAVELENGTH uncertainty in angstroms - which is how numset.dat states it too -
# and the wavelength is the standard one, vacuum below 2000 A and air above, so
# it is not 1e8 divided by the wavenumber and cannot be rebuilt from it without
# a dispersion formula.  It is rescaled instead; see rewrite_dlv.
DLV_INTENS = (0, 5)
DLV_WN = (5, 19)
DLV_LAMBDA = (19, 33)
DLV_UNC = (47, 60)
DLV_ROW = (60, 66)
DLV_WIDTH = 66
# A row of dlv.dat and a row of the line list are the same observed line when
# their wavenumbers agree to this much.  dlv.dat carries three decimals and the
# line list carries more, so the difference is rounding and nothing else.
DLV_MATCH = 0.001
# The intensity code is round(10*ln(Icalc)); an Icalc of zero or less has no
# logarithm, and a transition of exactly zero wavenumber none either.
MIN_INTENSITY = 1e-300


class SyncError(Exception):
    """The files cannot be brought into step, and nothing has been written."""


# ---------------------------------------------------------------------------
# The pieces that are read
# ---------------------------------------------------------------------------
def read_lopt_levels(path):
    """``{level_id: (energy, uncertainty)}`` from ``LOPT_output_levels.txt``.

    The uncertainty IDEN2 carries is the larger of LOPT's ``D1``, the
    statistical uncertainty of the level relative to its neighbours, and
    ``D2tot``, its uncertainty relative to the ground level, rounded to a
    thousandth of a wavenumber.  That is the rule the present ``enlev.dat``
    already follows, to the last digit, for all 636 levels.
    """
    t = pd.read_csv(path, sep='\t', dtype={'Designation': str})
    for column in ('Designation', 'Energy', 'D1', 'D2tot'):
        if column not in t.columns:
            raise SyncError('%s has no %s column' % (path, column))
    out = {}
    for row in t.itertuples():
        unc = round(max(float(row.D1), float(row.D2tot)), 3)
        out[row.Designation] = (float(row.Energy), unc)
    return out


def read_lopt_transitions(path):
    """``{(low_id, upp_id): (accepted, wn)}`` from ``LOPT_input_lines.txt``.

    A record flagged ``P`` is one LOPT is told to exclude from the level
    optimization.  A transition can have more than one record - the hyperfine
    components of one observed line are written separately - so it counts as
    accepted if any of its records is not flagged.
    """
    out = {}
    with io.open(path, encoding='latin-1', newline='') as fh:
        for rec in fh:
            text = rec.rstrip('\r\n')
            if not text.strip():
                continue
            wn = float(text[FIELD_WN[0]:FIELD_WN[1]])
            low = text[FIELD_LOW[0]:FIELD_LOW[1]].strip()
            upp = text[FIELD_UPP[0]:FIELD_UPP[1]].strip()
            accepted = 'P' not in text[FIELD_FLAGS[0]:FIELD_FLAGS[1]]
            was = out.get((low, upp))
            if was is None or (accepted and not was[0]):
                out[(low, upp)] = (accepted, wn)
    return out


# The columns of LOPT_input_lines.txt, as make_LOPT_input.FIELDS writes them,
# in the half-open form the rest of this program uses.
FIELD_WN = (0, 12)
FIELD_UNC = (13, 19)
FIELD_LOW = (42, 55)
FIELD_UPP = (58, 71)
FIELD_FLAGS = (72, 77)


def read_line_list(cfg, log):
    """``(wn_key, wn, u_wn)`` for every observed line of the set, sorted.

    The set's own workbook and its own column names, so that the baseline set
    is read on Sugar's wavenumbers and a corrected set on the corrected ones.
    ``wn_key`` is the immutable name of the line - Sugar's wavenumber, which
    never moves - and it is what a row of ``dlv.dat`` is matched on, because
    that file was built on that scale and its row numbers are quoted all
    through ``trans.dat``.
    """
    columns = cfg.lines.columns
    key_name = columns.get('wn_key') or columns['wn']
    frame = pd.read_excel(cfg.lines_file, sheet_name=cfg.lines.sheet or 0)
    for name in (key_name, columns['wn'], columns['u_wn']):
        if name not in frame.columns:
            raise SyncError('%s has no %s column' % (cfg.lines_file, name))
    table = frame[[key_name, columns['wn'], columns['u_wn']]].dropna()
    log('observed lines: %s' % cfg.lines_file)
    log('  %d lines; the wavenumber is its %s column and the uncertainty its '
        '%s' % (len(table), columns['wn'], columns['u_wn']))
    if key_name != columns['wn']:
        log('  matched to dlv.dat on %s, the wavenumber that never moves'
            % key_name)
    # The line list names a blended feature once per component, so the same
    # observed line can appear on several rows.  dlv.dat has one row per
    # observed line and no way to hold two wavenumbers for it, so the
    # repetitions are collapsed - and they must agree, because if two rows
    # sharing a wavenumber carried different corrected wavenumbers there would
    # be no answer to which of them IDEN2 should show.
    seen, clash = {}, []
    for a, b, c in table.to_numpy():
        key = round(float(a), 3)
        value = (float(a), float(b), float(c))
        if key in seen:
            if (abs(seen[key][1] - value[1]) > 5e-4
                    or abs(seen[key][2] - value[2]) > 5e-5):
                clash.append((key, seen[key], value))
            continue
        seen[key] = value
    if clash:
        raise SyncError(
            '%d observed wavenumber(s) appear more than once in %s with '
            'different corrected values, and dlv.dat has one row for each: '
            '%s.  Nothing has been written.'
            % (len(clash), os.path.basename(cfg.lines_file),
               ', '.join('%.3f' % c[0] for c in clash[:10])))
    if len(seen) != len(table):
        log('  %d of them are repeated rows of the same observed line, which '
            'is one row of dlv.dat; %d distinct lines remain'
            % (len(table) - len(seen), len(seen)))
    out = sorted(seen.values())
    return out


def intensity_model(icalc_path, log):
    """(C, kT) for Icalc = C*gA*(rwn/1e8)*exp(-Eup/kT), verified on the file.

    The constants come from ``[intensity_model]`` in the configuration, which
    is where the rest of the pipeline takes them from.  ``Icalc.xlsx`` is then
    re-fitted and the two compared: if the file no longer obeys the configured
    constants, the codes written here would not be the codes the pipeline
    works in, and that is a stop rather than a warning.
    """
    model = config.load().intensity_model
    C = float(model['C'])
    kT = float(model['kT'])
    tol = float(model.get('verify_tolerance', 0.01))
    log(f"intensity model  Icalc = C * gA * (rwn/1e8) * exp(-Eup/kT)")
    log(f"  C = {C:g}, kT = {kT:g} cm^-1, from the configuration")
    if not os.path.exists(icalc_path):
        log(f"  {os.path.basename(icalc_path)} not found; the constants were "
            f"not checked against it")
        return C, kT
    df = pd.read_excel(icalc_path)
    fit_C, fit_kT = gA_imputation.fit_intensity_model(df)
    err = gA_imputation.intensity_model_error(df, C, kT)
    log(f"  a fresh fit of {os.path.basename(icalc_path)} gives "
        f"C = {fit_C:g}, kT = {fit_kT:g}")
    log(f"  the configured constants reproduce its Icalc column to "
        f"{err:.2e}")
    if abs(fit_C / C - 1.0) > tol or abs(fit_kT / kT - 1.0) > tol:
        raise SyncError(
            'the configured intensity model (C = %g, kT = %g) disagrees with '
            'a fit of %s (C = %g, kT = %g) by more than %g.  The intensity '
            'codes written into trans.dat would not be the ones the rest of '
            'the pipeline works in.  Re-fit the file or correct '
            '[intensity_model] in the configuration.'
            % (C, kT, os.path.basename(icalc_path), fit_C, fit_kT, tol))
    return C, kT


# ---------------------------------------------------------------------------
# The new numbers
# ---------------------------------------------------------------------------
def rewrite_dlv(records, lines, log):
    """``dlv.dat`` on the set's own wavenumbers, and a report of the change.

    Every row keeps its number, its place in the file, its intensity code and
    its character: the row number is how ``trans.dat`` names the line, so
    nothing may be inserted, removed or reordered here.  What is rewritten is
    the wavenumber, the wavelength and the uncertainty.

    The wavelength is the standard one - vacuum in the ultraviolet, air above
    2000 A - so it is not a function of the wavenumber alone, and the
    dispersion formula the file was built with is not recorded anywhere.  It
    is therefore rescaled rather than recomputed,

        lambda_new = lambda_old * wn_old / wn_new ,

    which leaves the refractive index the row already implies exactly where it
    is: the corrections are a few hundredths of a wavenumber in ten thousand,
    and the index changes by nothing measurable over that.  The uncertainty is
    converted to the wavelength scale the file states it on by
    ``u_lambda = u_wn * lambda / wn``.

    A row whose wavenumber matches no line of the list is left as it stands
    and reported; so is a line of the list that has no row here, which cannot
    be given one without renumbering the file.
    """
    keys = [k for k, _wn, _u in lines]
    out, changed, unmatched = [], [], []
    for rec in records:
        if len(rec) < DLV_WIDTH or not rec.strip():
            out.append(rec)
            continue
        wn_old = float(rec[DLV_WN[0]:DLV_WN[1]])
        lam_old = float(rec[DLV_LAMBDA[0]:DLV_LAMBDA[1]])
        u_old = float(rec[DLV_UNC[0]:DLV_UNC[1]])
        row = int(rec[DLV_ROW[0]:DLV_ROW[1]])
        k = bisect.bisect_left(keys, wn_old)
        best = None
        for m in (k - 1, k, k + 1):
            if 0 <= m < len(keys) and abs(keys[m] - wn_old) <= DLV_MATCH:
                if best is None or abs(keys[m] - wn_old) < abs(keys[best]
                                                               - wn_old):
                    best = m
        if best is None:
            unmatched.append((row, wn_old))
            out.append(rec)
            continue
        _key, wn_new, u_wn = lines[best]
        lam_new = lam_old * wn_old / wn_new
        u_lam = u_wn * lam_new / wn_new
        rec = IDEN.put(rec, DLV_WN, '%14.3f' % wn_new)
        rec = IDEN.put(rec, DLV_LAMBDA, '%14.4f' % lam_new)
        rec = IDEN.put(rec, DLV_UNC, '%13.4f' % u_lam)
        out.append(rec)
        if (abs(wn_new - wn_old) > 0.0005 or abs(u_lam - u_old) > 0.00005):
            changed.append((row, wn_old, wn_new, u_old, u_lam, u_wn))
    absent = []
    dlv_keys = sorted(float(rec[DLV_WN[0]:DLV_WN[1]]) for rec in records
                      if len(rec) >= DLV_WIDTH and rec.strip())
    for key, wn_new, u_wn in lines:
        k = bisect.bisect_left(dlv_keys, key)
        near = min((abs(dlv_keys[m] - key) for m in (k - 1, k)
                    if 0 <= m < len(dlv_keys)), default=1e9)
        if near > DLV_MATCH:
            absent.append((key, wn_new, u_wn))
    report = {'n_rows': len(records), 'changed': changed,
              'unmatched': unmatched, 'absent': absent}
    return out, report


def dlv_rows_by_wavenumber(records):
    """``{wavenumber rounded to 3 decimals: row number}`` for dlv.dat.

    How a transition's assignment names its observed line: the last field of
    the assignment is the row number, and this is the way back from a
    wavenumber the pipeline knows to the row IDEN2 knows.
    """
    out = {}
    for rec in records:
        if len(rec) < DLV_WIDTH or not rec.strip():
            continue
        out[round(float(rec[DLV_WN[0]:DLV_WN[1]]), 3)] = (
            int(rec[DLV_ROW[0]:DLV_ROW[1]]),
            int(float(rec[DLV_INTENS[0]:DLV_INTENS[1]])))
    return out


def adopted_energies(enlev, id_of_row, lopt, log):
    """The energy and uncertainty every level should carry, and what changed.

    Returns ``(energies, changes)``.  ``energies`` is ``{index: (E, unc,
    known)}`` for every row of ``enlev.dat``: the LOPT energy for a level that
    has been found, and the level's own calculated energy, untouched, for one
    that has not.  ``changes`` lists the found levels whose energy or
    uncertainty moves.

    The energies are rounded here, to the three decimals ``enlev.dat`` and
    ``trans.dat`` hold, and everything downstream works in the rounded values.
    That is what makes the two files agree with each other: IDEN2 takes a
    predicted wavenumber to be the difference of the two energies exactly as
    they are written, and a difference of unrounded energies would miss it by
    a thousandth wherever the two roundings went opposite ways.
    """
    energies, changes = {}, []
    missing = []
    for index in sorted(enlev.row_of):
        lid = id_of_row.get(index)
        old_E, old_u = enlev.e_obs(index), IDEN.number(enlev.record(index),
                                                       IDEN.EN_UNC)
        if lid is None:
            if enlev.known(index):
                missing.append(index)
            energies[index] = (old_E, old_u, enlev.known(index))
            continue
        if lid not in lopt:
            raise SyncError(
                'enlev.dat row %d is mapped to %s by IDEN_level_ids.txt, and '
                'that level is not in the LOPT output.  The two files '
                'describe different sets of levels; sort that out first.'
                % (index, lid))
        E, unc = lopt[lid]
        E = round(E, 3)
        energies[index] = (E, unc, True)
        if abs(E - old_E) > 0.0005 or abs(unc - old_u) > 0.0005:
            changes.append((index, lid, old_E, E, old_u, unc))
    if missing:
        log(f"  warning: {len(missing)} levels are starred in enlev.dat and "
            f"are not in IDEN_level_ids.txt, so this run has no measured "
            f"energy for them and leaves them as they are: "
            + ', '.join(str(v) for v in missing[:10]))
    return energies, changes


def configuration_uncertainties(enlev, energies, log):
    """Give every level that has not been found the uncertainty of its
    configuration, and say what moved.

    The uncertainty column of ``enlev.dat`` carries two different quantities,
    depending on the row.  On a found level it is the uncertainty of a
    measurement: how well the fit knows where that level is.  On a level that
    has not been found there is no measurement, and the only honest reading of
    the column is a prediction: how far from its CALCULATED position the level
    is likely to turn out to be, once somebody finds it.  IDEN2 uses the
    column that way - it is the half-width of the window it will search - and
    the placeholder 5000 the conversion from Cowan's output writes there says
    only that nobody has ever filled it in.

    The scale of the calculation's error is not one number for the whole
    spectrum: it depends on the configuration, because the calculation
    describes some configurations much better than others.  So it is measured
    per configuration, as the rms of E_obs - E_calc over the levels of that
    configuration that have been found, using the energies this run is about
    to write.  The configuration is read from the label of ``enlev.dat``
    itself (``f26p``, ``fd6p``, ...), which is where IDEN2 keeps it.

    A configuration with fewer than ``MIN_CFG_LEVELS`` found levels has
    nothing to average - six of them in Pr III have none at all - and those
    levels are left exactly as they are, carrying whatever ``enlev.dat``
    already holds, which for a level nobody has ever touched is the
    placeholder 5000.  Borrowing the rms of the level list as a whole would
    be worse than saying nothing: it would claim a window of about 130 cm^-1
    for a level whose configuration may well be one of the badly calculated
    ones, where the true scatter is three times that.

    Returns ``(energies, report)``.  ``energies`` is the input with the
    uncertainty of every unfound level of a measured configuration replaced;
    ``report`` carries ``table`` - one ``(cfg, n_known, rms, n_unfound,
    left_alone)`` row per configuration - ``global_rms``, which is quoted for
    scale only, and ``changes``.
    """
    cfg_of, omc = {}, {}
    for index in sorted(enlev.row_of):
        cfg = level_interchange.split_label(enlev.label(index))[0]
        cfg_of[index] = cfg
        E, unc, known = energies[index]
        if known:
            omc.setdefault(cfg, []).append(E - round(enlev.e_calc(index), 3))

    def rms(values):
        return round(float(np.sqrt(np.mean(np.square(values)))), 3)

    every = [v for values in omc.values() for v in values]
    global_rms = rms(every) if every else 0.0
    per_cfg = {cfg: rms(v) for cfg, v in omc.items()
               if len(v) >= MIN_CFG_LEVELS}

    changes, counted = [], {}
    for index in sorted(enlev.row_of):
        E, unc, known = energies[index]
        if known:
            continue
        cfg = cfg_of[index]
        counted[cfg] = counted.get(cfg, 0) + 1
        if cfg not in per_cfg:
            continue
        new = per_cfg[cfg]
        energies[index] = (E, new, known)
        if abs(new - unc) > 0.0005:
            changes.append((index, cfg, unc, new))

    table = []
    for cfg in sorted(set(cfg_of.values())):
        table.append((cfg, len(omc.get(cfg, [])),
                      per_cfg.get(cfg, float('nan')), counted.get(cfg, 0),
                      cfg not in per_cfg))
    return energies, {'table': table, 'global_rms': global_rms,
                      'changes': changes}


def predicted_transitions(trans_table, mapping, energies, C, kT):
    """Every calculated transition between two levels of ``enlev.dat``.

    Returns a DataFrame with ``upper`` and ``lower`` (IDEN2 indices, the upper
    being the level of higher adopted energy), ``rwn`` the predicted
    wavenumber, ``Eup`` the upper level's adopted energy, ``gA``, ``Icalc``
    and ``code``, the integer round(10*ln(Icalc)) that ``trans.dat`` carries.
    """
    t = trans_table
    i1 = t['lid1'].map(mapping)
    i2 = t['lid2'].map(mapping)
    keep = i1.notna() & i2.notna()
    i1 = i1[keep].to_numpy(int)
    i2 = i2[keep].to_numpy(int)
    gA = t['gA'].to_numpy(float)[keep.to_numpy()]

    E1 = np.array([energies[k][0] for k in i1])
    E2 = np.array([energies[k][0] for k in i2])
    upper = np.where(E2 >= E1, i2, i1)
    lower = np.where(E2 >= E1, i1, i2)
    Eup = np.maximum(E1, E2)
    rwn = np.abs(E2 - E1)

    Icalc = gA_imputation.impute_intensity(gA, Eup, rwn, C, kT)
    Icalc = np.where(np.isfinite(Icalc) & (Icalc > MIN_INTENSITY),
                     Icalc, MIN_INTENSITY)
    return pd.DataFrame({'upper': upper, 'lower': lower, 'rwn': rwn,
                         'Eup': Eup, 'gA': gA, 'Icalc': Icalc,
                         'code': np.round(10.0 * np.log(Icalc)).astype(int)})


# ---------------------------------------------------------------------------
# Writing trans.dat
# ---------------------------------------------------------------------------
def sync_assignments(trans, id_of_row, lopt_lines, dlv_by_wn,
                     keep_unlisted, log):
    """Make the identified lines of ``trans.dat`` the ones LOPT is given.

    IDEN2 shows an identification as a line assigned to a predicted
    transition.  The pipeline's statement of the same thing is a record of
    ``LOPT_input_lines.txt``: a transition with a record that is not flagged
    ``P`` is one the fit uses, and a transition all of whose records are
    flagged ``P`` is one that was considered and excluded.  The two drift
    apart at every classification run, and an assignment IDEN2 still shows for
    an excluded transition is an assignment that has been withdrawn.

    So each of the four cases gets its own treatment:

    * the transition is accepted - the assignment stays, and its observed
      wavenumber and its row in ``dlv.dat`` are refreshed, because a corrected
      set moves every observed wavenumber;
    * every record of it is flagged ``P`` - the assignment is removed;
    * the LOPT input does not mention it at all - the classification run no
      longer proposes it, so the assignment is removed as well, unless
      ``keep_unlisted`` says to leave it.  That switch is there for the one
      case where removing it would be wrong: lines just marked by hand in
      IDEN2, which the pipeline has not been told about yet;
    * it is accepted and carries no assignment - one is made, when the
      observed line has a row in ``dlv.dat`` to point at.

    ``trans.records`` is rewritten in place for the transitions that have a
    row.  The ones that have none are returned in ``wanted``, for
    ``rebuild_trans`` to write when it builds the file.
    """
    row_of_id = {}
    for row, level_id in id_of_row.items():
        row_of_id[level_id] = row
    state, unknown_ids = {}, set()
    for (low, upp), (accepted, wn) in lopt_lines.items():
        if low not in row_of_id or upp not in row_of_id:
            unknown_ids.add(low if low not in row_of_id else upp)
            continue
        key = tuple(sorted((row_of_id[low], row_of_id[upp])))
        was = state.get(key)
        if was is None or (accepted and not was[0]):
            state[key] = (accepted, wn)
    if unknown_ids:
        raise SyncError(
            'the LOPT input names %d level(s) that IDEN_level_ids.txt has no '
            'row for: %s.  Every level the pipeline uses is in enlev.dat, so '
            'the lookup table is out of date rather than the level being '
            'absent; nothing has been written.'
            % (len(unknown_ids), ', '.join(sorted(unknown_ids)[:10])))

    removed_p, removed_unlisted, refreshed, no_dlv_row = [], [], [], []
    added, not_addable = [], []
    for (owner, partner), k in sorted(trans.row_of.items()):
        rec = trans.records[k]
        tail = IDEN.assignment(rec)
        if not IDEN.has_line(tail):
            continue
        key = tuple(sorted((owner, partner)))
        known = state.get(key)
        if known is None:
            if keep_unlisted:
                continue
            removed_unlisted.append((IDEN.obs_wavenumber(tail), key,
                                     IDEN.obs_row(tail)))
            trans.records[k] = rec[:IDEN.TR_OBS] + IDEN.BLANK_OBS
            continue
        accepted, wn = known
        if not accepted:
            removed_p.append((IDEN.obs_wavenumber(tail), key,
                              IDEN.obs_row(tail)))
            trans.records[k] = rec[:IDEN.TR_OBS] + IDEN.BLANK_OBS
            continue
        hit = dlv_by_wn.get(round(wn, 3))
        if hit is None:
            no_dlv_row.append((wn, key))
            continue
        row, code = hit
        old_wn = IDEN.obs_wavenumber(tail)
        new_tail = make_assignment(code, wn, IDEN.obs_omc(tail), row)
        if new_tail != tail:
            refreshed.append((old_wn, wn, key))
            trans.records[k] = rec[:IDEN.TR_OBS] + new_tail

    marked = set()
    for (owner, partner), k in trans.row_of.items():
        if IDEN.has_line(IDEN.assignment(trans.records[k])):
            marked.add(tuple(sorted((owner, partner))))
    wanted = {}
    for key, (accepted, wn) in sorted(state.items()):
        if not accepted or key in marked:
            continue
        hit = dlv_by_wn.get(round(wn, 3))
        if hit is None:
            not_addable.append((wn, key))
            continue
        row, code = hit
        # The departure from the prediction is not known until the row's
        # predicted wavenumber is in hand, so it is left at zero here and
        # recomputed by rebuild_trans along with every other row's.
        tail = make_assignment(code, wn, 0.0, row)
        k = trans.row(*key)
        if k is None:
            wanted[key] = tail
        else:
            trans.records[k] = trans.records[k][:IDEN.TR_OBS] + tail
        added.append((wn, key, k is None))

    report = {'removed_p': removed_p, 'removed_unlisted': removed_unlisted,
              'refreshed': refreshed, 'no_dlv_row': no_dlv_row,
              'added': added, 'not_addable': not_addable,
              'n_lopt': len(state),
              'n_accepted': sum(1 for v in state.values() if v[0])}
    return wanted, report


def make_assignment(code, wn, omc, row):
    """The assignment tail of a ``trans.dat`` row, at its fixed widths."""
    text = '%6d%12.3f%11.3f%6d' % (code, wn, omc, row)
    if len(text) != len(IDEN.BLANK_OBS):
        raise SyncError('built an assignment of %d characters, not %d: %r'
                        % (len(text), len(IDEN.BLANK_OBS), text))
    return text


def existing_rows(trans):
    """``{(upper, lower): (code, assignment)}`` for the file as it stands.

    The assignment is the tail of the row from column 38 on - the identified
    line, its departure from the prediction and its row number in IDEN2's own
    line list - kept verbatim so that nothing the analyst has done is lost.
    """
    out = {}
    for (owner, partner), k in trans.row_of.items():
        rec = trans.records[k]
        key = (min(owner, partner), max(owner, partner))
        out[key] = (int(IDEN.field(rec, TR_ICALC)),
                    IDEN.assignment(rec))
    return out


TR_ICALC = (5, 10)       # the intensity code, round(10*ln(Icalc))
TH_J = (11, 15)          # the block header: J of the level
TH_LABEL = (17, 28)      # the block header: the level's label
ROW_WIDTH = 72
HEADER_WIDTH = 44


def transition_row(partner, code, e_partner, partner_known, rwn, assignment):
    """One ``+`` row of ``trans.dat``, at its fixed width of 72 characters."""
    rec = '+%4d%5d%12.3f%s%13.3f' % (partner, code, e_partner,
                                     IDEN.STAR_ON if partner_known
                                     else IDEN.STAR_OFF, rwn)
    if len(rec) != IDEN.TR_OBS:
        raise SyncError('built a transition row of %d characters, not the %d '
                        'the format allows: %r'
                        % (len(rec), IDEN.TR_OBS, rec))
    return rec + assignment


def header_row(old_header, index, e_obs, known):
    """The ``$`` row that opens a level's block, with its energy refreshed."""
    rec = IDEN.put(old_header, IDEN.TH_EOBS, '%13.3f' % e_obs)
    rec = IDEN.put(rec, IDEN.TH_STAR,
                   IDEN.STAR_ON if known else IDEN.STAR_OFF)
    if int(IDEN.field(rec, IDEN.TH_INDEX)) != index:
        raise SyncError('header for level %d carries the number %s'
                        % (index, IDEN.field(rec, IDEN.TH_INDEX)))
    return rec


def build_header(index, enlev):
    """A block header for a level that has none, built from ``enlev.dat``.

    The label is taken between the two slashes of the ``enlev.dat`` row
    exactly as it stands, spaces and all.  It is eleven characters wide and
    the spaces are part of it - they separate the configuration from the
    parent term - so it can be neither stripped nor rebuilt.
    """
    label = enlev.record(index).partition('/')[2].partition('/')[0]
    rec = ('$%4d    J=%4s  %-11.11s%13.3f%s '
           % (index, enlev.j(index), label, enlev.e_obs(index),
              IDEN.STAR_ON if enlev.known(index) else IDEN.STAR_OFF))
    if len(rec) != HEADER_WIDTH:
        raise SyncError('built a block header of %d characters, not %d: %r'
                        % (len(rec), HEADER_WIDTH, rec))
    return rec


def rebuild_trans(trans, enlev, energies, predictions, cutoff, log,
                  wanted=None):
    """The new ``trans.dat`` as a list of records, and a report of the change.

    Blocks come in order of the upper level's IDEN2 index and rows within a
    block in order of the partner's index, which is the order the file already
    has.  A transition is written when it carries an identified line or when
    its intensity code reaches ``cutoff``.
    """
    old = existing_rows(trans)
    pred = {(int(u), int(l)): (int(c), float(w), float(i))
            for u, l, c, w, i in zip(predictions['upper'],
                                     predictions['lower'],
                                     predictions['code'],
                                     predictions['rwn'],
                                     predictions['Icalc'])}

    wanted = dict(wanted or {})
    rows = {}                       # (upper, lower) -> (code, assignment)
    kept_no_gA, added, dropped, code_change = [], [], [], []
    no_row = []                     # wanted assignments with nowhere to go
    for key, (old_code, tail) in old.items():
        upper, lower = key
        if energies[upper][0] < energies[lower][0]:
            upper, lower = lower, upper
        assigned = IDEN.has_line(tail)
        got = pred.get((upper, lower)) or pred.get((lower, upper))
        if got is None:
            # No gA: the hand-added row, or a transition the Cowan table no
            # longer carries.  Keep it when it holds an identification and the
            # code it already has; there is nothing to recompute it from.
            if assigned:
                kept_no_gA.append(key)
                rows[(upper, lower)] = (old_code, tail)
            else:
                dropped.append((key, old_code, 'the Cowan table no longer '
                                               'has this transition'))
            continue
        code = got[0]
        if code < cutoff and not assigned:
            dropped.append((key, old_code, 'below the cutoff'))
            continue
        rows[(upper, lower)] = (code, tail)
        if code != old_code:
            code_change.append(code - old_code)

    for (upper, lower), (code, _wn, _I) in pred.items():
        if (upper, lower) in rows or (lower, upper) in rows:
            continue
        tail = wanted.pop(tuple(sorted((upper, lower))), None)
        # A transition the pipeline has accepted is written whatever its
        # intensity code, exactly as one that already carries a line is: the
        # cutoff decides what is worth looking at, not what has been found.
        if tail is not None or code >= cutoff:
            rows[(upper, lower)] = (code, tail or IDEN.BLANK_OBS)
            added.append((upper, lower, code))
    # An accepted transition with neither a row in the old file nor a
    # calculated one to hang a row on cannot be shown at all.
    no_row = sorted(wanted)

    # ---- assemble ----------------------------------------------------------
    by_owner = {}
    for (upper, lower), value in rows.items():
        by_owner.setdefault(upper, []).append((lower, value))
    records, n_reassigned = [], 0
    for owner in sorted(by_owner):
        e_owner, _, known = energies[owner]
        if owner in trans.header_of:
            head = header_row(trans.records[trans.header_of[owner]],
                              owner, e_owner, known)
        else:
            head = build_header(owner, enlev)
        records.append(head)
        for lower, (code, tail) in sorted(by_owner[owner]):
            e_partner, _, partner_known = energies[lower]
            rwn = e_owner - e_partner
            if IDEN.has_line(tail):
                tail = IDEN.set_omc(tail, IDEN.obs_wavenumber(tail) - rwn)
                n_reassigned += 1
            records.append(transition_row(lower, code, e_partner,
                                          partner_known, rwn, tail))

    empty = sorted(set(trans.header_of) - set(by_owner))
    # How much of what could be seen is on the list: a transition between two
    # levels that have both been found is one whose wavenumber is known and
    # which the analyst can therefore go and look for, so the fraction of
    # those the cutoff keeps is the one that matters.
    known_pairs = [(u, l) for (u, l) in pred
                   if energies[u][2] and energies[l][2]]
    report = {
        'n_before': len(old), 'n_after': len(rows),
        'added': added, 'dropped': dropped, 'kept_no_gA': kept_no_gA,
        'code_change': np.array(code_change, dtype=int),
        'n_assignments': n_reassigned,
        'blocks_before': len(trans.header_of), 'blocks_after': len(by_owner),
        'blocks_emptied': empty,
        'n_known_pairs': len(known_pairs),
        'n_known_kept': sum(1 for k in known_pairs
                            if k in rows or (k[1], k[0]) in rows),
        'no_row': no_row,
    }
    return records, report


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------
def endings_for(records, template):
    """End-of-line strings for ``records``, in the style the file already has.

    IDEN2 writes its files with CRLF and git stores them with LF; whichever
    the file on disk has, the rewritten file keeps.
    """
    if template:
        style = max(set(template), key=template.count)
    else:
        style = '\r\n'
    return [style] * len(records)


def backup(path, suffix, log):
    if not suffix:
        return
    dest = path + suffix
    shutil.copy2(path, dest)
    log(f"  the previous {os.path.basename(path)} is kept as "
        f"{os.path.basename(dest)}")


def report_dlv(rep, log):
    """What the rewrite of ``dlv.dat`` changed."""
    changed = rep['changed']
    log(f"dlv.dat: {rep['n_rows']} rows, {len(changed)} rewritten")
    if changed:
        dwn = np.array([new - old for _r, old, new, _uo, _un, _u in changed])
        moved = np.abs(dwn[dwn != 0])
        if len(moved):
            log(f"  {len(moved)} wavenumbers move, by up to "
                f"{moved.max():.3f} cm^-1 (median {np.median(moved):.4f})")
        big = sorted(changed, key=lambda c: -abs(c[2] - c[1]))[:10]
        log(f"  {'row':>5} {'was':>12} {'now':>12} {'move':>8} "
            f"{'u_A was':>9} {'u_A now':>9} {'u_cm-1':>9}")
        for row, old, new, u_old, u_new, u_wn in big:
            log(f"  {row:5d} {old:12.3f} {new:12.3f} {new - old:+8.3f} "
                f"{u_old:9.4f} {u_new:9.4f} {u_wn:9.4f}")
    if rep['unmatched']:
        log(f"  {len(rep['unmatched'])} rows match no line of the list and "
            f"are left exactly as they were:")
        for row, wn in rep['unmatched'][:20]:
            log(f"    row {row:5d}  {wn:12.3f}")
    if rep['absent']:
        log(f"  {len(rep['absent'])} lines of the list have no row here.  A "
            f"row cannot be added without renumbering the file, which every "
            f"assignment in trans.dat refers to, so they are only listed:")
        for key, wn, _u in rep['absent'][:20]:
            log(f"    {wn:12.3f}" + ('' if abs(key - wn) < 5e-4
                                     else f"  (was {key:.3f})"))


def report_assignments(rep, log):
    """What the sync of the identifications changed."""
    log(f"  {rep['n_lopt']} transitions, {rep['n_accepted']} of them "
        f"accepted; the rest are flagged P, which is to say considered and "
        f"excluded")
    log(f"identifications in trans.dat:")
    log(f"  {len(rep['refreshed'])} kept, with the observed wavenumber and "
        f"the dlv.dat row refreshed")
    log(f"  {len(rep['removed_p'])} removed because every record of the "
        f"transition is flagged P")
    for wn, (a, b), row in rep['removed_p'][:40]:
        log(f"    {wn:12.3f}  levels {a:4d} - {b:4d}  dlv row {row:5d}")
    log(f"  {len(rep['removed_unlisted'])} removed because the LOPT "
        f"transitions file does not mention the transition at all")
    for wn, (a, b), row in rep['removed_unlisted'][:40]:
        log(f"    {wn:12.3f}  levels {a:4d} - {b:4d}  dlv row {row:5d}")
    log(f"  {len(rep['added'])} added for accepted transitions that carried "
        f"none")
    for wn, (a, b), new_row in rep['added'][:40]:
        log(f"    {wn:12.3f}  levels {a:4d} - {b:4d}"
            + ('  (a new row)' if new_row else ''))
    if rep['no_dlv_row'] or rep['not_addable']:
        stuck = sorted(rep['no_dlv_row'] + rep['not_addable'])
        log(f"  {len(stuck)} accepted transitions name an observed "
            f"wavenumber that has no row in dlv.dat; they are left as they "
            f"are, and the line has to be put into dlv.dat before IDEN2 can "
            f"show it:")
        for wn, (a, b) in stuck[:20]:
            log(f"    {wn:12.3f}  levels {a:4d} - {b:4d}")


def parse_args(argv):
    p = argparse.ArgumentParser(
        description='Rewrite IDEN2 enlev.dat and trans.dat from the current '
                    'LOPT levels and the current calculated intensities.')
    p.add_argument('--iden2', default=None, metavar='DIR',
                   help="the IDEN2 directory (default: the working set's)")
    p.add_argument('--lopt-levels', default=None, metavar='PATH',
                   help="LOPT output level table (default: the working set's)")
    p.add_argument('--lopt-lines', default=None, metavar='PATH',
                   help="the LOPT transitions file, which says which "
                        "identifications the fit is given (default: the "
                        "working set's LOPT_input_lines.txt)")
    p.add_argument('--config', default=None, metavar='PATH',
                   help="the configuration naming the set's line list "
                        "(default: the working set's lineclass_config.toml)")
    p.add_argument('--no-lines', action='store_true',
                   help='leave dlv.dat and the identifications alone and '
                        'rewrite only the levels and the predictions, as this '
                        'program did before it read the observed lines')
    p.add_argument('--keep-unlisted', action='store_true',
                   help='keep an identification that the LOPT transitions '
                        'file does not mention at all, instead of removing '
                        'it.  Use it when lines have just been marked by '
                        'hand in IDEN2 and not yet put through '
                        'classify_lines.py')
    p.add_argument('--set', metavar='DIR', default=None, dest='set_dir',
                   help='the working set to sync: its IDEN2 and its LOPT '
                        'levels, falling back to the project directory for '
                        'the files it has not got (default: the current '
                        'directory)')
    p.add_argument('--allow-mixed', action='store_true',
                   help='sync an IDEN2 directory with a LOPT level table '
                        'from a different set.  Refused by default, because '
                        'the two sets are on different wavenumber scales')
    p.add_argument('--tp', default=cowan_gA.TP_FILE, metavar='PATH',
                   help='the Cowan transition probabilities '
                        '(default %(default)s)')
    p.add_argument('--icalc', default=ICALC_FILE, metavar='PATH',
                   help='the calculated intensities the configured C and kT '
                        'are checked against (default %(default)s)')
    p.add_argument('--cutoff', type=int, default=DEF_CUTOFF, metavar='N',
                   help='write a transition when its intensity code '
                        'round(10*ln(Icalc)) reaches N.  A transition '
                        'carrying an identified line is written whatever its '
                        'code (default %(default)s)')
    p.add_argument('--no-cache', action='store_true',
                   help='re-read the transition-probability workbook instead '
                        'of the comma-separated cache beside it')
    p.add_argument('--backup-suffix', default=BACKUP_SUFFIX, metavar='EXT',
                   help='keep the previous files under this suffix '
                        '(default %(default)s)')
    p.add_argument('--no-backup', dest='backup_suffix', action='store_const',
                   const='', help='overwrite without keeping a copy')
    p.add_argument('--dry-run', action='store_true',
                   help='report what would change and write nothing')
    p.add_argument('--report', metavar='PATH',
                   help='write the report to this file as well as the screen')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    lines = []

    def log(text=''):
        print(text)
        lines.append(text)

    # The script lives in the project directory but is run from whichever
    # working set is being iterated, so IDEN2 and the LOPT level table are
    # resolved against that set first and against the project only for what
    # the set has not got - the arrangement swap_paths describes and the one
    # check_sync.py --set already follows.  Resolving them against the
    # script's own directory instead would quietly rewrite the baseline
    # IDEN2 from a run made somewhere else.
    set_dir = os.path.abspath(args.set_dir or os.getcwd())
    if not os.path.isdir(set_dir):
        raise SystemExit('sync_IDEN2.py: --set %s is not a directory'
                         % args.set_dir)
    if args.iden2 is None:
        args.iden2 = working_path('IDEN2', cwd=set_dir)
    if args.lopt_levels is None:
        args.lopt_levels = working_path('LOPT_output_levels.txt', cwd=set_dir)
    if args.lopt_lines is None:
        args.lopt_lines = working_path('LOPT_input_lines.txt', cwd=set_dir)
    if args.config is None:
        args.config = working_path('lineclass_config.toml', cwd=set_dir)
    if args.report and not os.path.isabs(args.report):
        args.report = os.path.join(set_dir, args.report)

    enlev_path = os.path.join(args.iden2, 'enlev.dat')
    trans_path = os.path.join(args.iden2, 'trans.dat')
    dlv_path = os.path.join(args.iden2, 'dlv.dat')
    map_path = os.path.join(args.iden2, 'IDEN_level_ids.txt')
    needed = [enlev_path, trans_path, map_path, args.lopt_levels, args.tp]
    written = [enlev_path, trans_path]
    if not args.no_lines:
        needed += [dlv_path, args.lopt_lines, args.config]
        written.append(dlv_path)
    for path in needed:
        if not os.path.exists(path):
            raise SystemExit('%s does not exist' % path)
    if not args.dry_run:
        output_files.require_writable(written, 'IDEN2 file')
    if args.report:
        output_files.require_writable([args.report], 'report file')

    # A level table from one set and an IDEN2 from another describe two
    # different fits: the corrected set's levels are on the corrected
    # wavenumber scale and the baseline's are on Sugar's, so writing one
    # into the other moves every Ritz wavenumber by the calibration
    # correction without saying so.
    iden2_home = os.path.dirname(os.path.abspath(args.iden2))
    homes = [('LOPT levels', args.lopt_levels)]
    if not args.no_lines:
        homes.append(('LOPT lines', args.lopt_lines))
    mixed = [(what, path) for what, path in homes
             if os.path.dirname(os.path.abspath(path)) != iden2_home]
    if mixed and not args.allow_mixed:
        raise SystemExit('\n'.join([
            '',
            'sync_IDEN2.py: these files belong to different working sets:',
            '    IDEN2        %s' % args.iden2,
            ] + ['    %-12s %s' % (what, path) for what, path in mixed] + [
            'Their wavenumbers are on different scales.  Run the set through '
            'the pipeline',
            'and LOPT first, or name the files explicitly, or pass '
            '--allow-mixed if this is meant.',
            '']))

    log('working set: %s' % set_dir)
    log('  IDEN2        %s' % args.iden2)
    log('  LOPT levels  %s' % args.lopt_levels)
    if args.no_lines:
        log('  --no-lines: dlv.dat and the identifications are left alone')
    else:
        log('  LOPT lines   %s' % args.lopt_lines)
        log('  config       %s' % args.config)

    C, kT = intensity_model(args.icalc, log)

    log()
    trans_table = cowan_gA.read_transitions(args.tp, cache=not args.no_cache,
                                            log=log)
    calc_levels = cowan_gA.levels(trans_table)
    enlev_levels = cowan_gA.read_enlev_levels(enlev_path)
    id_of_row = cowan_gA.read_id_map(map_path)
    mapping, match = cowan_gA.match_to_enlev(calc_levels, enlev_levels,
                                             id_of_row)
    log(f"{len(trans_table)} calculated transitions among "
        f"{len(calc_levels)} calculated levels")
    log(f"  {match['n_matched']} of them matched to rows of enlev.dat "
        f"({match['n_reordered']} recovered from an order reversal); "
        f"the calculated energies agree to {match['max_dE']:.1f} cm^-1")
    log(f"  the match was checked against all {match['n_checked']} levels "
        f"whose IDEN2 row and experimental identifier are both known")
    if match['unmatched_lid']:
        log(f"  {len(match['unmatched_lid'])} calculated levels have no row "
            f"in enlev.dat and are ignored: "
            + ', '.join(str(v) for v in match['unmatched_lid'][:10]))
    if match['unmatched_idx']:
        log(f"  {len(match['unmatched_idx'])} rows of enlev.dat have no "
            f"calculated level: "
            + ', '.join(str(v) for v in match['unmatched_idx'][:10]))

    enlev = IDEN.Enlev(enlev_path)
    trans = IDEN.Trans(trans_path)
    lopt = read_lopt_levels(args.lopt_levels)
    row_id = {n: lid for n, lid in id_of_row.items() if n in enlev.row_of}
    energies, changes = adopted_energies(enlev, row_id, lopt, log)

    log()
    log(f"enlev.dat: {len(lopt)} measured levels, {len(enlev.row_of)} rows")
    moves = [c for c in changes if abs(c[3] - c[2]) > 0.0005]
    if moves:
        moved = np.array([abs(new - old) for _, _, old, new, _, _ in moves])
        log(f"  {len(moves)} levels move; the largest move is "
            f"{moved.max():.3f} cm^-1, the median {np.median(moved):.4f}")
        big = sorted(moves, key=lambda c: -abs(c[3] - c[2]))[:10]
        log(f"  {'row':>5}  {'level_id':<14} {'was':>12} {'now':>12} "
            f"{'move':>8}")
        for index, lid, old, new, _, _ in big:
            log(f"  {index:5d}  {lid:<14} {old:12.3f} {new:12.3f} "
                f"{new - old:+8.3f}")
    else:
        log("  every level energy already agrees with the LOPT output")

    # The uncertainty of a found level: how well the fit knows where it is.
    u_moves = [c for c in changes if abs(c[5] - c[4]) > 0.0005]
    if u_moves:
        du = np.array([new - old for _, _, _, _, old, new in u_moves])
        log(f"  {len(u_moves)} found levels change uncertainty to "
            f"max(D1, D2tot) of the LOPT output; {int((du > 0).sum())} grow, "
            f"{int((du < 0).sum())} shrink, the largest change is "
            f"{np.abs(du).max():.3f} cm^-1 and the median "
            f"{np.median(np.abs(du)):.4f}")
        big = sorted(u_moves, key=lambda c: -abs(c[5] - c[4]))[:10]
        log(f"  {'row':>5}  {'level_id':<14} {'was':>10} {'now':>10} "
            f"{'change':>8}")
        for index, lid, _, _, old, new in big:
            log(f"  {index:5d}  {lid:<14} {old:10.3f} {new:10.3f} "
                f"{new - old:+8.3f}")
    else:
        log("  every found level already carries max(D1, D2tot)")

    # The uncertainty of a level that has NOT been found: how far from its
    # calculated position it is likely to be, which is a property of its
    # configuration rather than of the level.
    energies, cfg_rep = configuration_uncertainties(enlev, energies, log)
    log()
    log(f"the levels that have not been found take the uncertainty of their "
        f"configuration - the rms of E_obs - E_calc over the levels of that "
        f"configuration that HAVE been found.  Over the level list as a whole "
        f"that rms is {cfg_rep['global_rms']:.3f} cm^-1, but it varies by a "
        f"factor of nearly twenty between configurations, which is the reason "
        f"for doing it this way:")
    log(f"  {'config':<10} {'found':>6} {'rms':>9} {'unfound':>8}")
    for cfg, n_known, rms_v, n_unfound, alone in cfg_rep['table']:
        if not n_unfound and not n_known:
            continue
        log(f"  {cfg:<10} {n_known:6d} "
            + (f"{'-':>9} " if alone else f"{rms_v:9.3f} ")
            + f"{n_unfound:8d}"
            + ('   no found level of its own; these are left as they are'
               if alone and n_unfound else ''))
    log(f"  {len(cfg_rep['changes'])} unfound levels change uncertainty")

    dlv_records, dlv_rep, assign_rep = None, None, None
    wanted = {}
    if not args.no_lines:
        log()
        lines = read_line_list(config.load(args.config), log)
        dlv_records, dlv_endings = IDEN.read_records(dlv_path)
        dlv_records, dlv_rep = rewrite_dlv(dlv_records, lines, log)
        report_dlv(dlv_rep, log)
        lopt_lines = read_lopt_transitions(args.lopt_lines)
        log()
        log('%s: %d records' % (args.lopt_lines, len(lopt_lines)))
        wanted, assign_rep = sync_assignments(
            trans, row_id, lopt_lines, dlv_rows_by_wavenumber(dlv_records),
            args.keep_unlisted, log)
        report_assignments(assign_rep, log)

    predictions = predicted_transitions(trans_table, mapping, energies, C, kT)
    records, rep = rebuild_trans(trans, enlev, energies, predictions,
                                 args.cutoff, log, wanted=wanted)

    log()
    log(f"trans.dat: {rep['n_before']} transitions before, "
        f"{rep['n_after']} after (cutoff {args.cutoff})")
    log(f"  {len(rep['added'])} added, {len(rep['dropped'])} dropped, "
        f"{rep['n_assignments']} identified lines carried through")
    if rep['kept_no_gA']:
        log(f"  {len(rep['kept_no_gA'])} transitions carry an identified line "
            f"and have no gA in the Cowan table; their rows and their "
            f"intensity codes are kept as they are:")
        for upper, lower in sorted(rep['kept_no_gA']):
            log(f"    levels {upper} - {lower}")
    gone = [d for d in rep['dropped'] if d[2] != 'below the cutoff']
    if gone:
        log(f"  {len(gone)} transitions are no longer in the Cowan table and "
            f"carry no identification; they are dropped:")
        for (upper, lower), code, _ in sorted(gone)[:20]:
            log(f"    levels {upper} - {lower}, old code {code}")
    delta = rep['code_change']
    if len(delta):
        log(f"  {len(delta)} intensity codes change, by "
            f"{np.percentile(delta, 5):+.0f} to "
            f"{np.percentile(delta, 95):+.0f} at the 5th and 95th percentile "
            f"(median {np.median(delta):+.0f}); one code unit is a factor "
            f"{math.exp(0.1):.3f} in intensity")
    if rep['blocks_emptied']:
        log(f"  {len(rep['blocks_emptied'])} levels lose every DOWNWARD "
            f"transition they had - they head no block now, and IDEN2 will "
            f"show nothing below them: "
            + ', '.join(str(v) for v in rep['blocks_emptied'][:10]))
    log(f"  {rep['blocks_after']} levels head a block, "
        f"{rep['blocks_before']} did before")
    log(f"  of the {rep['n_known_pairs']} calculated transitions between two "
        f"levels that have both been found - the ones whose wavenumber is "
        f"known and which can therefore be looked for - "
        f"{rep['n_known_kept']} are on the list; the rest are below the "
        f"cutoff")
    if rep['no_row']:
        log(f"  {len(rep['no_row'])} accepted transitions have neither a row "
            f"in the old file nor a calculated one, so they cannot be shown; "
            f"their identification is not written:")
        for upper, lower in rep['no_row'][:20]:
            log(f"    levels {upper} - {lower}")

    log()
    if args.dry_run:
        log("--dry-run: nothing was written")
    else:
        backup(enlev_path, args.backup_suffix, log)
        backup(trans_path, args.backup_suffix, log)
        if dlv_records is not None:
            backup(dlv_path, args.backup_suffix, log)
            IDEN.write_records(dlv_path, dlv_records, dlv_endings)
            log(f"  {os.path.basename(dlv_path)} rewritten, "
                f"{len(dlv_records)} rows")
        for index, (E, unc, known) in energies.items():
            enlev.set_measurement(index, unc, E, known)
        IDEN.write_records(enlev_path, enlev.records, enlev.endings)
        log(f"  {os.path.basename(enlev_path)} rewritten, "
            f"{len(enlev.records)} rows")
        IDEN.write_records(trans_path, records,
                           endings_for(records, trans.endings))
        log(f"  {os.path.basename(trans_path)} rewritten, "
            f"{len(records)} rows")
        log()
        log("Reopen IDEN2 to pick the new files up.")

    if args.report:
        with open(args.report, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write('\n'.join(lines) + '\n')
        print(f"report written to {args.report}")
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (SyncError, cowan_gA.MatchError) as exc:
        sys.exit(str(exc))
