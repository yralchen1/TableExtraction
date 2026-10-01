"""
make_LOPT_input.py
==================
Build the three input files of LOPT (the level-optimization code of
A. Kramida) from the classification table produced by classify_lines.py.

LOPT needs
  * a TRANSITIONS file  - one fixed-column record per classified line,
  * a FIXED LEVELS file - the levels whose energies must not be moved,
  * a PARAMETER file    - the file names, the options, and the column
                          positions of every field of the transitions file.

This script reads line_classifications.csv and writes

  LOPT_input_lines.txt  the transitions file, in exactly the column layout
                        of the sample Pr3_line_class13.prn,
  LOPT_fixlev.txt       the fixed levels, copied verbatim from the sample
                        fixed-levels file,
  LOPT.par              the parameter file, copied from the sample parameter
                        file with only the four file names replaced.

Rules applied to the transitions file
  * only rows that carry both a lower and an upper level identifier are
    written (an unclassified observed line has no transition to optimise);
  * a row with accepted = 0 gets the flag "P" in the flags column.  "P"
    means "predicted": LOPT prints the line but gives it no weight, so the
    rejected classifications do not influence the optimised levels, and its
    weight is written as 0.0000;
  * an observed line with a single accepted classification (n_accepted = 1)
    gets the weight 1.0000: the whole measured intensity belongs to it;
  * an observed line with several accepted classifications - a blend - has
    its weight split between them in proportion to their calculated
    intensities (column calc_intens), because that is the best estimate of
    how much of the blend each component contributes.  The proportions are
    normalised so that the components of one blend sum to 1.0000.

Why the weights of a blend are normalised here.  LOPT normalizes them
itself, so only their ratios matter; what normalizing buys is that every
weight then fits the six characters of the weight column at four decimals,
and the transitions file keeps exactly the column layout of the sample - no
column position in the parameter file has to be touched.  Raw calc_intens
values reach ~9e5 and would not fit.  The precision cost is nil for the
present table: the smallest normalised weight is 0.0204, still three
significant figures, and no component rounds to 0.0000.  If a future table
ever holds a blend so lopsided that a real component would round away, the
script says so and stops rather than writing a silent zero.

The uncertainty written for a line is not the one quoted in the line list.
The quoted value, unc_wn_obs, is the uncertainty of the *measurement* - the
reading error of its era, the precision floor, and whatever the character flag
of the line says about its width.  A line also carries the hyperfine structure
of the two levels it joins: a level whose nuclear spin splits it into unresolved
components displaces every one of its lines, and level_positions.py fits one
width per level from the residuals of the 1974 lines and writes them to
level_hfs_widths.csv, in its w_applied column - which for a level with too few
lines of its own is the typical width of the level's configuration, that being
what governs a hyperfine splitting.  That part cannot be quoted per line,
because it depends
on which two levels the line was assigned to, and the two components of a blend
may sit on levels of different widths.  So it is added here, where the
assignment is known:

    unc = sqrt(unc_wn_obs^2 + w_hfs(lower)^2 + w_hfs(upper)^2)

That is a property of the *transition*, and one observed line may be assigned
to several of them, whose levels carry different widths.  LOPT, though, reads
one record per component and takes the uncertainty of each as the uncertainty
of the measurement it constrains, so two components of one blend written with
two uncertainties are two different weights on the same measured wavenumber -
which is a statement the measurement cannot make.  An observed line therefore
gets one uncertainty for all of its records: the mean of the transitions'
values weighted by the very weights LOPT is given, the calculated intensity
fractions, so that the component which carries the line governs its width.
Records flagged P take the same value although they are outside the fit, since
they describe the same measurement.  A line whose records are all flagged has
no weights to average with, and takes the plain mean.

One transition, one record.  The same pair of levels can be reached from
several observed lines - a line of the list is assigned to it and rejected,
another is assigned to it and accepted, or two rejected candidates name it at
two wavenumbers.  A transition, however, has one energy difference, so a second
record for it tells LOPT nothing it does not already have, and only one is
written.  The accepted record wins, there being at most one: a transition
accepted at two wavenumbers would be one energy difference measured twice by
one line list, and the script stops rather than choose.  Where no record is
accepted, all of them are flagged P and carry no weight, so LOPT prints the
first and which wavenumber that is does not matter.  The rest are dropped
before the weights and the uncertainty of their observed lines are worked out,
so a record that is not written cannot affect one that is.

The blending factor k(n) of level_positions.py is deliberately NOT included.
LOPT.par sets BLEND TREATMENT = centroid, so LOPT already compares the observed
wavenumber against the intensity-weighted centroid of a blend's components and
reports the correction as dEcent; k(n) measures the spread of the components
about that centroid, which LOPT has thus already removed, and adding it would
charge the same effect twice.  --no-hfs turns the hyperfine term off.

The hyperfine correction.  With `apply = true` in the `[hfs]` section of the
set's configuration (the `lineclass_config.toml` beside the classification
table, or the project's), the levels are fitted in the head frame of
hfs_correction.py: every observed line with an accepted classification is
written at wn_obs + hfs_shift, the column classify_lines.py computed with the
same switch on, (1 - kappa) * sum BF_i * D_i over the line's accepted
transitions, and its uncertainty u_hfs_shift is added to the line's in
quadrature.  All records of one observed line take the same shift, so that
LOPT still sees the components of a blend at one wavenumber.  The hyperfine
width of a level whose hyperfine structure the correction accounts for (a
determined A constant, or resolved sublevels) is no longer added, nor any
width on a line Sugar flagged, which sits on the head-frame Ritz value; a
level without a determined A keeps its width (Work_on_hfs_plan.md, D13).
What was added to each record is written to LOPT_hfs_shifts.txt beside the
transitions file, so that sync_IDEN2.py and check_sync.py, which find the
observed line of a LOPT record by its wavenumber, can take the measured value
back.  A table written with the switch in the other position stops the run:
a correction applied to one file and not the other would fit a mixture of
frames.

The rows are written in order of decreasing wavenumber, as in the sample.

Usage
    python make_LOPT_input.py
    python make_LOPT_input.py --classifications my_lines.csv --par-out run7.par
    cd iter && python ../make_LOPT_input.py     # builds iter's three files

The working set
    The classification table says which set is being built, and the three
    output files are written beside it.  An input named by its bare name -
    the sample parameter file, the sample fixed-levels file,
    level_hfs_widths.csv - is looked for beside that table first and in the
    project directory only if the set has not got it, so running the script
    from an iteration folder finds the project's copies instead of failing
    on them or silently doing without them.  All three resolved paths are
    printed.
"""

import argparse
import csv
import math
import os

import config
import hfs_correction
import hfs_kappa
import output_files
from swap_paths import working_path

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEF_CLASSIFICATIONS = 'line_classifications.csv'
DEF_SAMPLE_PAR = 'Pr3_line_class13.par'
DEF_SAMPLE_FIXLEV = 'Pr3_line_class5_fixlev.txt'
DEF_HFS_WIDTHS = 'level_hfs_widths.csv'
DEF_CONFIG = 'lineclass_config.toml'
# The columns classify_lines.py writes with [hfs] apply on.
HFS_COLUMNS = ('kappa', 'hfs_D', 'hfs_shift', 'u_hfs_shift')

# The admission rule for a fitted hyperfine width, as in level_positions.py:
# a width below HFS_APPLY changes no sigma measurably and is consistent with
# the noise of a fit that cannot return a negative width, and a level with
# fewer than HFS_MIN_LINES lines of 1974 has too few residuals to fit one.
HFS_APPLY = 0.02        # cm^-1
HFS_MIN_LINES = 4
DEF_LINES_OUT = 'LOPT_input_lines.txt'
DEF_FIXLEV_OUT = 'LOPT_fixlev.txt'
DEF_PAR_OUT = 'LOPT.par'
DEF_LEV_OUT = 'LOPT_output_levels.txt'
DEF_LIN_OUT = 'LOPT_output_lines.txt'

# Column layout of the transitions file, taken from the sample parameter
# file Pr3_line_class13.par.  Each entry is (first column, last column) with
# 1-based, inclusive columns, i.e. exactly the numbers LOPT is given.  The
# fields are left-justified inside their span and padded with blanks.
FIELDS = {
    'wavenumber':  (1, 12),
    'uncertainty': (14, 19),
    'intensity':   (25, 42),
    'lower_level': (43, 55),
    'upper_level': (59, 71),
    'flags':       (73, 77),
    'weight':      (82, 87),
    'units':       (90, 93),
}

# LOPT reads the wavelength/wavenumber unit of each line from the units
# field; "cm-1" says the first field is a wavenumber in cm^-1.
UNITS = 'cm-1'

# The sample transitions file uses DOS line endings; LOPT accepts either,
# but the output is kept byte-compatible with the sample.
EOL_LINES = '\r\n'
EOL_PAR = '\n'


# ---------------------------------------------------------------------------
# The transitions file
# ---------------------------------------------------------------------------
def format_line(wn, unc, intens, low_id, upp_id, flag, weight):
    """Return one record of the transitions file, without the line ending.

    Every value is placed at the column LOPT is told to read it from; the
    record is built as a list of blanks that the fields are written into.
    """
    width = max(last for _, last in FIELDS.values())
    buf = [' '] * width

    def put(field, text):
        first, last = FIELDS[field]
        span = last - first + 1
        if len(text) > span:
            raise ValueError(
                f'value {text!r} does not fit in the {span}-character '
                f'{field} field (columns {first}-{last})')
        buf[first - 1:first - 1 + len(text)] = text

    put('wavenumber', f'{wn:.3f}')
    put('uncertainty', f'{unc:.3f}')
    put('intensity', f'{intens:.3f}')
    put('lower_level', low_id)
    put('upper_level', upp_id)
    put('flags', flag)
    put('weight', f'{weight:.4f}')
    put('units', UNITS)
    return ''.join(buf).rstrip()


def read_classifications(path):
    """Return the classified rows of the classification table.

    A row is classified when it names both an upper and a lower level;
    the unclassified observed lines carry no transition and are dropped.  So
    is the row of a resolved hfs companion (grade hfs, hfs_correction.py):
    it names its transition, but the transition's energy difference is
    measured by the main line, and the companion lies a rung of the pattern
    away from it.
    """
    with open(path, newline='', encoding='utf-8-sig') as fh:
        rows = [r for r in csv.DictReader(fh)
                if (r.get('low_id') or '').strip()
                and (r.get('upp_id') or '').strip()
                and (r.get('grade') or '').strip()
                != hfs_correction.COMPANION_GRADE]
    if not rows:
        raise SystemExit(f'{path}: no classified lines found')
    return rows


def read_corrections(path):
    """`{wavenumber: (correction, statistical uncertainty)}` in cm^-1.

    `wavelength_calibration_corrections.csv` is written by
    `wavelength_calibration.py`.  `own_correction` is what has to be added to
    Sugar's measured wavenumber to put the line on the corrected scale, and
    `u_stat_cm1` is what the line's uncertainty class says its measuring
    scatter is once that correction has been made - which is the uncertainty
    a corrected wavenumber is worth, not the one Sugar stated and not the one
    adopted before the curve was known.  A line the calibration could not
    reach carries a correction of zero and its class uncertainty all the
    same.

    The key is `'%.4f' % wn`, the precision the wavenumbers are quoted to.
    """
    out = {}
    with open(path, newline='', encoding='utf-8-sig') as fh:
        for r in csv.DictReader(fh):
            wn = (r.get('wn_obs') or '').strip()
            if not wn:
                continue
            d = (r.get('own_correction') or '').strip()
            u = (r.get('u_stat_cm1') or '').strip()
            out['%.4f' % float(wn)] = (float(d) if d else 0.0,
                                       float(u) if u else None)
    return out


def apply_corrections(rows, corrections):
    """Put the classified rows on the corrected wavenumber scale.

    Returns `(shifted, largest, unknown)`: how many rows moved, the largest
    shift in cm^-1, and the rows whose wavenumber the corrections file does
    not carry.  The rows are changed in place, so everything downstream -
    the repeat filter, the blend weights, the written file - works on the
    corrected values without knowing that they are corrected.
    """
    # A classification table written by a corrected set already carries
    # corrected wavenumbers, and says so: its wn_key column is the name of the
    # line, which is the wavenumber on the published scale, and the two differ
    # exactly where a correction has been applied.  Correcting again would
    # double every shift, silently, so it stops the run instead.
    already = [r for r in rows
               if (r.get('wn_key') or '').strip()
               and abs(float(r['wn_key']) - float(r['wn_obs'])) > 1e-4]
    if already:
        raise SystemExit(
            'make_LOPT_input.py: --corrections was given a classification '
            'table whose wavenumbers are already corrected (%d of %d rows '
            'differ from their wn_key, the first at %.4f cm-1).  Applying '
            'the corrections to it would count them twice.  Either build '
            'LOPT input from the baseline table with --corrections, or '
            'classify the corrected set and build from its own table without '
            'them - they are alternatives, not steps.'
            % (len(already), len(rows), float(already[0]['wn_key'])))

    shifted, largest, unknown = 0, 0.0, []
    for row in rows:
        wn = float(row['wn_obs'])
        rec = corrections.get('%.4f' % wn)
        if rec is None:
            unknown.append(wn)
            continue
        d, u = rec
        if d:
            row['wn_obs'] = repr(wn + d)
            shifted += 1
            largest = max(largest, abs(d))
        if u is not None:
            row['unc_wn_obs'] = repr(u)
    return shifted, largest, unknown


def is_accepted(row):
    return float(row['accepted'] or 0) == 1


def read_hfs_widths(path):
    """{level_id: hyperfine width} from level_hfs_widths.csv, in cm^-1.

    The file is written by level_positions.py.  Its w_applied column already
    holds that module's whole admission rule - a level the fit can measure
    carries its own fitted width, a level it cannot carries the typical width
    of its configuration - so reading that column is what keeps the two in
    step.  A file written before the column existed is read the old way,
    above HFS_APPLY on at least HFS_MIN_LINES lines, which was the rule then.
    A missing file means no hyperfine term is added at all.
    """
    if not path or not os.path.exists(path):
        return {}
    out = {}
    with open(path, newline='', encoding='utf-8') as fh:
        for rec in csv.DictReader(fh):
            if 'w_applied' in rec:
                try:
                    w = float(rec['w_applied'])
                except (TypeError, ValueError):
                    continue
                if w > 0:
                    out[rec['level_id'].strip()] = w
                continue
            try:
                w = float(rec['w_hfs'])
                n = float(rec['n_1974'])
            except (KeyError, TypeError, ValueError):
                continue
            if w > HFS_APPLY and n >= HFS_MIN_LINES:
                out[rec['level_id'].strip()] = w
    return out


def total_unc(unc, low_id, upp_id, w_hfs):
    """The uncertainty LOPT is given: the measurement and the two level widths."""
    a = w_hfs.get(low_id, 0.0)
    b = w_hfs.get(upp_id, 0.0)
    if not a and not b:
        return unc
    return (unc ** 2 + a ** 2 + b ** 2) ** 0.5


def one_record_per_transition(rows):
    """`rows` reduced to one record per pair of levels, and the ones dropped.

    Returns `(kept, dropped)`, `kept` in the order it was given.  A transition
    is one energy difference, so a second record for it adds nothing to the
    fit; the accepted record is the one kept, or the first where none is
    accepted - those are all flagged P, weightless, and LOPT prints whichever
    wavenumber it is given.  Two accepted records of one transition are a
    contradiction in the classification, not something to choose between, and
    stop the run.
    """
    seen = {}
    for r in rows:
        key = (r['low_id'].strip(), r['upp_id'].strip())
        chosen = seen.get(key)
        if chosen is None:
            seen[key] = r
        elif is_accepted(r):
            if is_accepted(chosen):
                raise SystemExit(
                    f'the transition {key[0]} - {key[1]} is accepted at two '
                    f'observed wavenumbers, {float(chosen["wn_obs"]):.3f} and '
                    f'{float(r["wn_obs"]):.3f} cm-1; one transition is one '
                    f'energy difference, so at most one of them can be its '
                    f'measurement - settle it in the classification')
            seen[key] = r
    kept_ids = {id(r) for r in seen.values()}
    return ([r for r in rows if id(r) in kept_ids],
            [r for r in rows if id(r) not in kept_ids])


def blend_uncertainty(uncs, weights):
    """The one uncertainty the records of a single observed line share.

    `uncs` are the per-transition uncertainties of those records and `weights`
    the weights LOPT is given, zero for a record it does not fit.  The result
    is their weighted mean, or the plain mean where nothing is fitted.
    """
    total = sum(weights)
    if total <= 0:
        return sum(uncs) / float(len(uncs))
    return sum(w * u for w, u in zip(weights, uncs)) / total


def blend_weights(rows):
    """Return {id(row): weight} for the accepted rows of one observed line.

    A single accepted classification takes the whole line, weight 1.  The
    components of a blend divide it in proportion to their calculated
    intensities.  Should those intensities be missing or add up to nothing,
    the components share the line equally - that is the only assumption
    left once the intensities say nothing.
    """
    if len(rows) == 1:
        return {id(rows[0]): 1.0}
    calc = [max(float(r['calc_intens'] or 0.0), 0.0) for r in rows]
    total = sum(calc)
    if total <= 0:
        calc, total = [1.0] * len(rows), float(len(rows))
    return {id(r): c / total for r, c in zip(rows, calc)}


def hfs_line_shifts(rows):
    """`{wn_obs: (hfs_shift, u_hfs_shift, kappa)}`, one entry per observed
    line, from the columns classify_lines.py writes with [hfs] apply on.

    The shift belongs to the line and is written on each of its rows; rows
    of one line that disagree say the table has been edited by hand, and
    stop the run rather than have one of them chosen.
    """
    out = {}
    for r in rows:
        rec = (float(r['hfs_shift'] or 0.0), float(r['u_hfs_shift'] or 0.0),
               float(r['kappa']))
        was = out.setdefault(r['wn_obs'], rec)
        if was != rec:
            raise SystemExit(
                f'the rows of the observed line {r["wn_obs"]} cm-1 carry '
                f'different hfs_shift, u_hfs_shift or kappa values '
                f'({was} and {rec}); they are one line\'s, so the table has '
                f'been changed since classify_lines.py wrote it')
    return out


def write_lines_file(rows, path, w_hfs=None, hfs=None, shifts_out=None,
                     hfs_stats=None):
    """Write the LOPT transitions file.

    Returns (written, accepted, flagged, widened, dropped).  `w_hfs` is
    {level_id: hyperfine width}; the width of the two levels a line joins is
    added to its quoted uncertainty in quadrature.  `widened` counts the
    records that got one, `dropped` the repeated records of a transition that
    is already written.

    `hfs`, an hfs_correction.Model, turns the hyperfine correction on (see
    the module docstring): each record is written at wn_obs + hfs_shift of
    its line, with u_hfs_shift added in quadrature, and the width of a level
    the correction accounts for is left out.  Every record whose wavenumber
    is changed is appended to `shifts_out` as (low_id, upp_id, wn_obs as
    the table gives it, the wavenumber written, shift, its uncertainty),
    and `hfs_stats['widths_left_out']` counts the records that lost a
    level width to it.
    """
    w_hfs = w_hfs or {}
    line_shift = hfs_line_shifts(rows) if hfs is not None else {}
    # A transition is one energy difference, however many observed lines have
    # been assigned to it.  The repeats go before anything else is computed,
    # so that a record which is not written cannot weigh on one that is.
    rows, dropped = one_record_per_transition(rows)
    # The accepted classifications of one observed line share its weight, so
    # they have to be weighed together; wn_obs identifies the observed line.
    groups = {}
    for r in rows:
        if is_accepted(r):
            groups.setdefault(r['wn_obs'], []).append(r)
    weights = {}
    for group in groups.values():
        weights.update(blend_weights(group))

    # And they share its uncertainty, which the hyperfine term would otherwise
    # make differ from one component to the next.
    per_row = {}
    all_rows = {}
    for r in rows:
        low, upp = r['low_id'].strip(), r['upp_id'].strip()
        widths = w_hfs
        if hfs is not None and w_hfs:
            # The widths the correction replaces (D13): none at all on a
            # line Sugar flagged, which sits on the head-frame Ritz value
            # (its blends' other components included), nor on the main line
            # of resolved companions, which does too, and none for a level
            # whose hyperfine structure is computed.
            flagged = (hfs_kappa.kappa_class(
                (r.get('char') or '').strip()) == 'flag'
                or hfs.is_head_line(float(r.get('wn_key') or r['wn_obs'])))
            widths = {lid: w for lid, w in w_hfs.items()
                      if lid in (low, upp) and not flagged
                      and not hfs.is_corrected(lid)}
            if hfs_stats is not None and len(widths) < sum(
                    1 for lid in (low, upp) if lid in w_hfs):
                hfs_stats['widths_left_out'] = (
                    hfs_stats.get('widths_left_out', 0) + 1)
        unc = total_unc(float(r['unc_wn_obs']), low, upp, widths)
        if hfs is not None:
            unc = math.hypot(unc, line_shift[r['wn_obs']][1])
        per_row[id(r)] = unc
        all_rows.setdefault(r['wn_obs'], []).append(r)
    shared = {}
    for key, group in all_rows.items():
        if len(group) > 1:
            shared[key] = blend_uncertainty(
                [per_row[id(r)] for r in group],
                [weights.get(id(r), 0.0) for r in group])

    records = []
    n_accepted_rows = 0
    n_widened = 0
    for r in rows:
        wn = float(r['wn_obs'])
        if hfs is not None:
            shift, u_shift = line_shift[r['wn_obs']][:2]
            if shift:
                wn = float('%.3f' % (wn + shift))
                if shifts_out is not None:
                    shifts_out.append((r['low_id'].strip(),
                                       r['upp_id'].strip(), r['wn_obs'],
                                       wn, shift, u_shift))
        quoted = float(r['unc_wn_obs'])
        unc = shared.get(r['wn_obs'], per_row[id(r)])
        if unc > quoted:
            n_widened += 1
        intens = float(r['obs_intens'])
        if is_accepted(r):
            flag, weight = '', weights[id(r)]
            if weight > 0 and round(weight, 4) == 0:
                raise SystemExit(
                    f'the component {r["low_id"]} - {r["upp_id"]} of the '
                    f'blend at {wn:.3f} cm-1 carries the weight {weight:.3g}, '
                    f'which rounds to zero in the {FIELDS["weight"][1] - FIELDS["weight"][0] + 1}'
                    f'-character weight column; widen that column in '
                    f'FIELDS and in the parameter file')
            n_accepted_rows += 1
        else:
            flag, weight = 'P', 0.0
        records.append((wn, format_line(wn, unc, intens,
                                        r['low_id'].strip(),
                                        r['upp_id'].strip(),
                                        flag, weight)))

    records.sort(key=lambda t: -t[0])
    with open(path, 'w', newline='', encoding='ascii') as fh:
        for _, text in records:
            fh.write(text + EOL_LINES)
    return (len(records), n_accepted_rows,
            len(records) - n_accepted_rows, n_widened, len(dropped))


def hfs_model(cfg_path, rows, source):
    """The hfs_correction.Model if the configuration at `cfg_path` switches
    the hyperfine correction on, else None.

    The classification table must have been written with the switch in the
    same position: its hfs columns are what the correction adds, and a
    table without them has had its candidates judged against the plain Ritz
    values, one with them against the corrected ones.
    """
    settings = config.load(cfg_path).hfs
    has = all(c in rows[0] for c in HFS_COLUMNS)
    if settings.apply and not has:
        raise SystemExit(
            f'make_LOPT_input.py: {cfg_path} switches the hyperfine '
            f'correction on ([hfs] apply = true), but {source} was written '
            f'with it off - it has no hfs_shift column.  Run '
            f'classify_lines.py on this set again first.')
    if has and not settings.apply:
        raise SystemExit(
            f'make_LOPT_input.py: {source} was written with the hyperfine '
            f'correction on, but {cfg_path} has it off ([hfs] apply = '
            f'false).  Run classify_lines.py on this set again, or turn the '
            f'switch back on.')
    return hfs_correction.Model(settings) if settings.apply else None


# ---------------------------------------------------------------------------
# The fixed-levels file
# ---------------------------------------------------------------------------
def write_fixlev_file(sample, path):
    """Copy the sample fixed-levels file byte for byte."""
    if os.path.abspath(sample) == os.path.abspath(path):
        return 0
    with open(sample, 'rb') as src:
        data = src.read()
    with open(path, 'wb') as dst:
        dst.write(data)
    return len(data)


# ---------------------------------------------------------------------------
# The parameter file
# ---------------------------------------------------------------------------
def replace_par_name(line, new_name):
    """Return `line` of the parameter file with its file name replaced.

    A parameter line is "<value><blanks>; <comment>".  Only the value is
    replaced; the comment is kept where it was whenever the new name is
    short enough, so the file stays as readable as the sample.
    """
    semi = line.find(';')
    if semi < 0:                       # no comment: the whole line is the name
        return new_name
    comment = line[semi:]
    pad = max(semi - len(new_name), 1)
    return new_name + ' ' * pad + comment


def write_par_file(sample, path, lines_name, fixlev_name,
                   lev_out_name, lin_out_name):
    """Copy the sample parameter file, replacing only the four file names.

    The first four lines of a LOPT parameter file are, in order, the
    transitions input, the fixed levels input, the levels output and the
    transitions output.  Everything below them - the options and the column
    positions - is copied unchanged.
    """
    with open(sample, encoding='utf-8') as fh:
        par = fh.read().splitlines()
    if len(par) < 4:
        raise SystemExit(f'{sample}: not a LOPT parameter file '
                         f'(only {len(par)} lines)')

    for i, name in enumerate((lines_name, fixlev_name,
                              lev_out_name, lin_out_name)):
        par[i] = replace_par_name(par[i], name)

    with open(path, 'w', newline='', encoding='utf-8') as fh:
        for text in par:
            fh.write(text + EOL_PAR)
    return len(par)


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description='Build the LOPT input files from line_classifications.csv',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--classifications', default=DEF_CLASSIFICATIONS,
                   help='classification table written by classify_lines.py')
    p.add_argument('--sample-par', default=DEF_SAMPLE_PAR,
                   help='parameter file whose options and column layout are '
                        'to be reused')
    p.add_argument('--sample-fixlev', default=DEF_SAMPLE_FIXLEV,
                   help='fixed-levels file to copy')
    p.add_argument('--hfs-widths', default=DEF_HFS_WIDTHS,
                   help='per-level hyperfine widths (level_hfs_widths.csv); '
                        'added to each line uncertainty in quadrature')
    p.add_argument('--no-hfs', action='store_true',
                   help='do not add the hyperfine widths to the uncertainties')
    p.add_argument('--config', default=DEF_CONFIG,
                   help='the configuration whose [hfs] section says whether '
                        'the hyperfine correction is applied; a bare name is '
                        "looked for beside the classification table first")
    p.add_argument('--lines-out', default=DEF_LINES_OUT,
                   help='transitions file to write')
    p.add_argument('--fixlev-out', default=DEF_FIXLEV_OUT,
                   help='fixed-levels file to write')
    p.add_argument('--par-out', default=DEF_PAR_OUT,
                   help='parameter file to write')
    p.add_argument('--levels-output', default=DEF_LEV_OUT,
                   help='name of the levels file LOPT itself will write')
    p.add_argument('--lines-output', default=DEF_LIN_OUT,
                   help='name of the transitions file LOPT itself will write')
    p.add_argument('--corrections', default='',
                   help='wavelength_calibration_corrections.csv: put the '
                        'lines on the corrected wavenumber scale and give '
                        'them the statistical uncertainty that goes with it '
                        '(this is what makes the "iter" and "final" sets)')
    p.add_argument('--outdir', default='',
                   help='directory for the three written files '
                        '(default: alongside the classification table)')
    p.add_argument('--unlock', action='store_true',
                   help='write a locked set (one whose own '
                        'lineclass_config.toml says locked = true, as the '
                        'baseline\'s does)')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    def in_dir(name):
        if os.path.isabs(name) or os.path.dirname(name):
            return name
        base = args.outdir or os.path.dirname(
            os.path.abspath(args.classifications))
        return os.path.join(base, name)

    def in_src(name):
        """An input file: the working set's copy, else the project's.

        The classification table says which set is being built, so an input
        named by its bare name is looked for beside that table first and in
        the project directory only if the set has not got it - the rule
        swap_paths describes.  The set's own copy therefore wins where it
        exists, and the files that exist in one place only - the sample
        parameter file, the sample fixed-levels file, the hyperfine widths -
        are found from anywhere instead of being silently skipped or
        crashing the run.  `--outdir` says where the run's own output goes
        and never where its inputs are looked for.
        """
        if os.path.isabs(name) or os.path.dirname(name):
            return name
        return working_path(
            name, cwd=os.path.dirname(os.path.abspath(args.classifications)))

    lines_out = in_dir(args.lines_out)
    fixlev_out = in_dir(args.fixlev_out)
    par_out = in_dir(args.par_out)

    # Resolved before anything is read or written, and reported below, so
    # that an input picked up from the project directory rather than from
    # the set is visible in the run's own output.
    sample_fixlev = in_src(args.sample_fixlev)
    sample_par = in_src(args.sample_par)
    hfs_widths = in_src(args.hfs_widths)
    cfg_path = in_src(args.config)
    shifts_file = hfs_correction.shifts_path(lines_out)
    for what, path in (('sample fixed levels', sample_fixlev),
                       ('sample parameter file', sample_par)):
        if not os.path.exists(path):
            raise SystemExit('make_LOPT_input.py: %s not found: %s'
                             % (what, path))

    # LOPT's input files are tab-delimited text, and an analyst who has one
    # of them open in Excel would otherwise learn of it only after the whole
    # classification had been re-read and re-written.
    config.require_unlocked([lines_out, fixlev_out, par_out, shifts_file],
                            'make_LOPT_input.py', args.unlock)
    output_files.require_writable([lines_out, fixlev_out, par_out,
                                   shifts_file])

    rows = read_classifications(args.classifications)
    hfs = hfs_model(cfg_path, rows, args.classifications)
    if args.corrections:
        shifted, largest, unknown = apply_corrections(
            rows, read_corrections(args.corrections))
        print(f'{args.corrections}: {shifted} of {len(rows)} records put on '
              f'the corrected scale, largest shift {largest:.4f} cm-1')
        if unknown:
            print(f'  {len(unknown)} record(s) are not in the corrections '
                  f"file and keep Sugar's wavenumber and uncertainty; "
                  f'the first is {unknown[0]:.4f}')
    w_hfs = {} if args.no_hfs else read_hfs_widths(hfs_widths)
    shifts, hfs_stats = [], {}
    written, accepted, flagged, widened, dropped = write_lines_file(
        rows, lines_out, w_hfs, hfs=hfs, shifts_out=shifts,
        hfs_stats=hfs_stats)
    # The record of the shifts describes the transitions file just written,
    # so one left from an earlier run with the correction on is removed
    # when the correction is off: it would describe a file that is gone.
    if hfs is not None:
        hfs_correction.write_shifts(shifts_file, shifts)
    elif os.path.exists(shifts_file):
        os.remove(shifts_file)
        print(f'{shifts_file}: removed; it described a transitions file '
              f'written with the hyperfine correction on')
    write_fixlev_file(sample_fixlev, fixlev_out)
    # The parameter file must name the input files as LOPT will look for
    # them; LOPT resolves them next to itself, so bare names are written.
    write_par_file(sample_par, par_out,
                   os.path.basename(lines_out), os.path.basename(fixlev_out),
                   args.levels_output, args.lines_output)

    print(f'classifications: {args.classifications}')
    print(f'configuration: {cfg_path}')
    print(f'sample fixed levels: {sample_fixlev}')
    print(f'sample parameter file: {sample_par}')
    print(f'{lines_out}: {written} transitions '
          f'({accepted} weighted, {flagged} flagged "P")')
    if dropped:
        print(f'  {dropped} repeated record(s) of a transition already '
              f'written were dropped')
    if w_hfs:
        print(f'  hyperfine widths from {hfs_widths}: {len(w_hfs)} levels '
              f'carry one; {widened} transitions had their uncertainty '
              f'widened by it'
              + (' or by the uncertainty of the hyperfine correction'
                 if hfs is not None else ''))
        if hfs is not None:
            print(f"  {hfs_stats.get('widths_left_out', 0)} transitions "
                  f'lost a level width that the hyperfine correction '
                  f'replaces (a determined A, or a flagged line)')
    elif args.no_hfs:
        print('  no hyperfine widths applied (--no-hfs)')
    else:
        print(f'  no hyperfine widths applied: {hfs_widths} does not exist')
    if hfs is not None:
        lines = {}
        for low, upp, wn_obs, _, d, u in shifts:
            lines[wn_obs] = (d, u)
        big = sorted((abs(d) for d, _ in lines.values()), reverse=True)
        print(f'  hyperfine correction on ([hfs] apply): {len(lines)} '
              f'observed lines ({len(shifts)} records) moved into the head '
              f'frame' + (f', largest |shift| {big[0]:.4f} cm-1, '
                          f'{sum(1 for b in big if b > 0.1)} above 0.1'
                          if big else ''))
        print(f'  {shifts_file}: what was added to each record')
    print(f'{fixlev_out}: copied from {sample_fixlev}')
    print(f'{par_out}: copied from {sample_par} with new file names')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
