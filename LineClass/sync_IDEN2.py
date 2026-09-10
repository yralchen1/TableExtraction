"""Bring IDEN2's ``enlev.dat`` and ``trans.dat`` up to date with the fit.

WHAT IS OUT OF DATE, AND WHY
============================
IDEN2 keeps two files of its own.  ``enlev.dat`` is the level list it shows on
screen: one row per calculated level, carrying the level's calculated energy,
the measured energy adopted for it, the uncertainty of that measurement and an
asterisk when the level has been found at all.  ``trans.dat`` is the predicted
transition list: one block per upper level, one row per transition, carrying a
copy of the partner's energy, the predicted wavenumber, a code for the
predicted intensity, and - where the analyst has identified the transition with
an observed line - that line and its departure from the prediction.

Both drift out of step with the pipeline, in three separate ways.

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
    new predicted wavenumber and the new intensity code, and a transition that
    carries an identified line keeps it, with its observed-minus-predicted
    recomputed against the new prediction.

WHAT IS NEVER LOST
==================
**An identification is never dropped.**  A transition that carries an observed
line stays in ``trans.dat`` however weak its new intensity code makes it, and
whether or not the Cowan table still has a gA for it.  One transition has no
gA at all: it is Sugar's identification of a line that fell below the gA = 1e3
s^-1 floor of the Cowan run and so was never calculated, adopted because the
observed line's ``*v`` character matches the large hyperfine structure of the
level, as the ``*r``, ``cl``, ``c`` and ``w`` characters of the transitions
upward from it do.  It was added to ``trans.dat`` by hand with an intensity
code of 0, it is the one imputed ``Icalc`` in ``line_classifications.csv``,
and this program carries its row through untouched but for the geometry.

WHICH TRANSITIONS ARE LISTED
============================
A transition is written when it carries an identified line, or when its new
intensity code reaches ``--cutoff``.  The default cutoff is -41, which is the
value that leaves the present list most nearly alone: it adds 420 transitions
that the old scale had cut off and drops 47 - 38 that no longer reach the
floor and 9 that have gone from the Cowan table - out of 104272.  A larger cutoff makes a shorter file - each step of one
removes about 700 transitions, being a factor 1.105 in intensity - and
``--cutoff -1000`` writes every calculated transition there is.

USAGE
=====
    python sync_IDEN2.py                 # rewrite both files, keeping backups
    python sync_IDEN2.py --dry-run       # say what would change, write nothing
    python sync_IDEN2.py --cutoff -45    # a longer list, down to weaker lines
    python sync_IDEN2.py --report sync_IDEN2_report.txt

Close IDEN2 first.  It holds its files open and rewrites them from memory when
it exits, which would undo everything this program has done.
"""

import argparse
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

HERE = os.path.dirname(os.path.abspath(__file__))
IDEN2_DIR = os.path.join(HERE, 'IDEN2')
LOPT_LEVELS = os.path.join(HERE, 'LOPT_output_levels.txt')
ICALC_FILE = os.path.join(HERE, 'Icalc.xlsx')
DEF_CUTOFF = -41
# An rms taken over one or two levels is not a measure of anything, so a
# configuration with fewer found levels than this is left alone.
MIN_CFG_LEVELS = 3
BACKUP_SUFFIX = '.presync'
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
    already follows, to the last digit, for all 594 levels.
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


def rebuild_trans(trans, enlev, energies, predictions, cutoff, log):
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

    rows = {}                       # (upper, lower) -> (code, assignment)
    kept_no_gA, added, dropped, code_change = [], [], [], []
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
        if code >= cutoff:
            rows[(upper, lower)] = (code, IDEN.BLANK_OBS)
            added.append((upper, lower, code))

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


def parse_args(argv):
    p = argparse.ArgumentParser(
        description='Rewrite IDEN2 enlev.dat and trans.dat from the current '
                    'LOPT levels and the current calculated intensities.')
    p.add_argument('--iden2', default=IDEN2_DIR, metavar='DIR',
                   help='the IDEN2 directory (default %(default)s)')
    p.add_argument('--lopt-levels', default=LOPT_LEVELS, metavar='PATH',
                   help='LOPT output level table (default %(default)s)')
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

    enlev_path = os.path.join(args.iden2, 'enlev.dat')
    trans_path = os.path.join(args.iden2, 'trans.dat')
    map_path = os.path.join(args.iden2, 'IDEN_level_ids.txt')
    for path in (enlev_path, trans_path, map_path, args.lopt_levels, args.tp):
        if not os.path.exists(path):
            raise SystemExit('%s does not exist' % path)
    if not args.dry_run:
        output_files.require_writable([enlev_path, trans_path], 'IDEN2 file')
    if args.report:
        output_files.require_writable([args.report], 'report file')

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

    predictions = predicted_transitions(trans_table, mapping, energies, C, kT)
    records, rep = rebuild_trans(trans, enlev, energies, predictions,
                                 args.cutoff, log)

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

    log()
    if args.dry_run:
        log("--dry-run: nothing was written")
    else:
        backup(enlev_path, args.backup_suffix, log)
        backup(trans_path, args.backup_suffix, log)
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
