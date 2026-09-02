"""
make_LOPT_input.py
==================
Build the three input files of LOPT (the level-optimisation code of
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

Why the weights of a blend are normalised here.  LOPT normalises them
itself, so only their ratios matter; what normalising buys is that every
weight then fits the six characters of the weight column at four decimals,
and the transitions file keeps exactly the column layout of the sample - no
column position in the parameter file has to be touched.  Raw calc_intens
values reach ~9e5 and would not fit.  The precision cost is nil for the
present table: the smallest normalised weight is 0.0204, still three
significant figures, and no component rounds to 0.0000.  If a future table
ever holds a blend so lopsided that a real component would round away, the
script says so and stops rather than writing a silent zero.

The rows are written in order of decreasing wavenumber, as in the sample.

Usage
    python make_LOPT_input.py
    python make_LOPT_input.py --classifications my_lines.csv --par-out run7.par
"""

import argparse
import csv
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEF_CLASSIFICATIONS = 'line_classifications.csv'
DEF_SAMPLE_PAR = 'Pr3_line_class13.par'
DEF_SAMPLE_FIXLEV = 'Pr3_line_class5_fixlev.txt'
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
    the unclassified observed lines carry no transition and are dropped.
    """
    with open(path, newline='', encoding='utf-8-sig') as fh:
        rows = [r for r in csv.DictReader(fh)
                if (r.get('low_id') or '').strip()
                and (r.get('upp_id') or '').strip()]
    if not rows:
        raise SystemExit(f'{path}: no classified lines found')
    return rows


def is_accepted(row):
    return float(row['accepted'] or 0) == 1


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


def write_lines_file(rows, path):
    """Write the LOPT transitions file; return (written, accepted, flagged)."""
    # The accepted classifications of one observed line share its weight, so
    # they have to be weighed together; wn_obs identifies the observed line.
    groups = {}
    for r in rows:
        if is_accepted(r):
            groups.setdefault(r['wn_obs'], []).append(r)
    weights = {}
    for group in groups.values():
        weights.update(blend_weights(group))

    records = []
    n_accepted_rows = 0
    for r in rows:
        wn = float(r['wn_obs'])
        unc = float(r['unc_wn_obs'])
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
    return len(records), n_accepted_rows, len(records) - n_accepted_rows


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
    p.add_argument('--outdir', default='',
                   help='directory for the three written files '
                        '(default: alongside the classification table)')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    def in_dir(name):
        if os.path.isabs(name) or os.path.dirname(name):
            return name
        base = args.outdir or os.path.dirname(
            os.path.abspath(args.classifications))
        return os.path.join(base, name)

    lines_out = in_dir(args.lines_out)
    fixlev_out = in_dir(args.fixlev_out)
    par_out = in_dir(args.par_out)

    rows = read_classifications(args.classifications)
    written, accepted, flagged = write_lines_file(rows, lines_out)
    write_fixlev_file(args.sample_fixlev, fixlev_out)
    # The parameter file must name the input files as LOPT will look for
    # them; LOPT resolves them next to itself, so bare names are written.
    write_par_file(args.sample_par, par_out,
                   os.path.basename(lines_out), os.path.basename(fixlev_out),
                   args.levels_output, args.lines_output)

    print(f'{lines_out}: {written} transitions '
          f'({accepted} weighted, {flagged} flagged "P")')
    print(f'{fixlev_out}: copied from {args.sample_fixlev}')
    print(f'{par_out}: copied from {args.sample_par} with new file names')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
