"""
lopt_separations.py
===================
The separations file of LOPT 6: how far apart the components of a blend may
lie.

LOPT 6 fits a blend (one observed line assigned to several transitions) as
one observation: the intensity-weighted mean of its components' Ritz
wavenumbers has to reproduce the observed wavenumber.  One number cannot place
the components, so LOPT adds "separation equations" that say the components
coincide within tau, tau = BLEND SEPARATION FACTOR x the width (FWHM) of the
line.  Without a resolution in the parameter file it takes the FWHM to be 10
times the line's uncertainty, which follows the measurement uncertainty, not
the spectrograph (user's guide, "Blends").

This module writes those equations itself, into LOPT's separations input file,
and make_LOPT_input.py then sets BLEND SEPARATION FACTOR = 0 so that LOPT's
own equations are off (they are not replaced by a file record, they are added
to it; lopt/observations.py, build_group_observations).  What changes:

* **The width is the line's physical one.**  The effective FWHM is the
  instrumental width level_shifts.INSTR_FWHM_A (0.035 A) combined in
  quadrature with the Doppler width level_shifts.doppler_fwhm_factor() x
  lambda, converted to cm^-1 at the line's wavenumber (as in make_set3.fwhm_cm).
* **Wide blends get a wider allowance.**  A line whose Sugar character, with
  the '**' of a multiple classification stripped, is one of
  [blends] wide_characters (w wide, c complex, d double; user, 2026-11-18)
  has tau = wide_separation_factor x FWHM; the others separation_factor x FWHM.
  On iter_hfs's LOPT 6 fit of 2026-11-18 the fitted component separations of
  those blends were 0.81 FWHM rms against 0.41 for the others, hence 1.0 and
  0.5.
* **A blend can be given a factor of its own** in the hand-kept registry
  files.blend_separations (blend_separations.txt), keyed on Sugar's
  wavenumber of the line (wn_key), as inflated_unc_lines.txt is.  It
  overrides both factors; it is for the blends whose components LOPT puts
  farther apart than their character allows (its warning "deviation > 2
  sigma").

Each pair (i, j) of the n components of a blend gets one record
``(E_up_i - E_low_i) - (E_up_j - E_low_j) = 0 +- u_ij`` with
u_ij = sqrt(sum_k w_k / (w_i w_j)), w_k = 2 / tau_k^2 - the weights LOPT
gives its own equations (create_blend_sep_obs), so a fit with this file and
the same tau is LOPT's own.  One difference: LOPT counts each record of the
file as one degree of freedom, its own equations as n - 1 per blend, so a blend
of 3 components adds 1 to the degrees of freedom and one of 4 adds 3.

The blends are found as LOPT finds them, from the transitions file itself:
records with the same text in the wavenumber field and the same unit, flag
neither P nor M, weight above zero.  So the file can be rebuilt from any
transitions file (refresh), which every tool that rewrites one must do -
a stale record holds together two transitions that are no longer one blend.

The file opens with comments that record the factors and the width model,
and each blend's records follow a comment naming the line and its character,
so that refresh() can rebuild the file with the same rules.
"""

import csv
import io
import math
import os
import re

import make_LOPT_input

#: the default name of the file, beside LOPT.par
DEF_NAME = 'LOPT_separations.txt'
#: LOPT's parameter keys (matched in upper case anywhere in a line)
PAR_SEP_FILE_KEY = 'SEPARATIONS INPUT FILE NAME'
PAR_SEP_FACTOR_KEY = 'BLEND SEPARATION FACTOR'
#: the '**' Sugar's list puts after a multiply classified line
MULTIPLE_MARK = '**'

_HEAD_RE = re.compile(r'^;\s*rules:\s*factor\s+(\S+)\s+wide_factor\s+(\S+)'
                      r'\s+wide_characters\s+(\S*)\s*$')
_BLEND_RE = re.compile(r'^;\s*blend\s+(\S+)\s+char\s+(\S*)\s')
_OWN_RE = re.compile(r'registry factor\s+([0-9.eE+-]+)')


def read_registry(path):
    """`{wn_key rounded to 4 decimals: factor}` of the hand-kept registry
    of blends with an allowance of their own: tab-delimited, header
    `wn_key  factor  date  reason`, `#` starts a comment line.  {} if `path`
    is '' or absent."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with io.open(path, encoding='utf-8-sig', newline='') as fh:
        rows = [l for l in fh.read().splitlines()
                if l.strip() and not l.lstrip().startswith('#')]
    for row in csv.DictReader(rows, delimiter='\t'):
        key = round(float(row['wn_key']), 4)
        factor = float(row['factor'])
        if not factor > 0:
            raise ValueError('%s: the factor of %s must be positive'
                             % (path, row['wn_key']))
        if key in out:
            raise ValueError('%s: %s is listed twice' % (path, row['wn_key']))
        out[key] = factor
    return out


def base_char(char):
    """Sugar's character with the multiple-classification mark stripped."""
    c = (char or '').strip()
    if c.endswith(MULTIPLE_MARK):
        c = c[:-len(MULTIPLE_MARK)]
    return c.strip()


def is_wide(char, wide_chars):
    """Whether a line of character `char` is a wide one."""
    return base_char(char) in set(wide_chars)


def fwhm_cm(wn):
    """The effective line width (FWHM) at `wn`, cm^-1: the instrumental width
    and the Doppler width, both in angstrom, in quadrature."""
    import level_shifts            # heavy; only needed when the file is built
    lam = 1.0e8 / wn
    width_A = math.hypot(level_shifts.INSTR_FWHM_A,
                         level_shifts.doppler_fwhm_factor() * lam)
    return width_A * wn * wn * 1.0e-8


# ---------------------------------------------------------------------------
# The transitions file, read the way LOPT reads it
# ---------------------------------------------------------------------------
def _field(text, name):
    first, last = make_LOPT_input.FIELDS[name]
    return text[first - 1:last]


def read_records(path):
    """The records of a LOPT transitions file as dicts: `wn_text` (the
    wavenumber field as written, LOPT's blend key), `wn`, `low`, `upp`,
    `flag`, `weight`, `unit`."""
    out = []
    with io.open(path, encoding='latin-1', newline='') as fh:
        for raw in fh.read().splitlines():
            if not raw.strip():
                continue
            wtext = _field(raw, 'weight').strip()
            # LOPT: an empty weight, or the exact text 0, means 1
            weight = 1.0 if wtext in ('', '0') else float(wtext)
            out.append({'wn_text': _field(raw, 'wavenumber').strip(),
                        'wn': float(_field(raw, 'wavenumber')),
                        'low': _field(raw, 'lower_level').strip(),
                        'upp': _field(raw, 'upper_level').strip(),
                        'flag': _field(raw, 'flags').upper(),
                        'weight': weight,
                        'unit': _field(raw, 'units').strip()})
    return out


def blends(records):
    """`[(wn_text, unit, [record, ...])]` of the blends, in file order: the
    records LOPT puts in one blend (same wavenumber text and unit, not P or
    M, a weight above zero), where there are at least two."""
    groups = {}
    order = []
    for r in records:
        if 'P' in r['flag'] or 'M' in r['flag'] or not r['weight'] > 0:
            continue
        key = (r['wn_text'], r['unit'])
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(r)
    return [(k[0], k[1], groups[k]) for k in order if len(groups[k]) > 1]


# ---------------------------------------------------------------------------
# The file
# ---------------------------------------------------------------------------
def pair_uncertainties(taus):
    """`[(i, j, u_ij)]` for every pair of components with the per-component
    allowances `taus`: LOPT's own weights, w = 2 / tau^2."""
    w = [2.0 / (t * t) for t in taus]
    sw = sum(w)
    return [(i, j, math.sqrt(sw / (w[i] * w[j])))
            for i in range(len(w)) for j in range(i + 1, len(w))]


def label(text, width):
    """A level label padded to the width of LOPT's label field."""
    if len(text) > width:
        raise ValueError('level label %r is wider than the %d-character '
                         'label field' % (text, width))
    return text.ljust(width)


def build(records, char_of, factor, wide_factor, wide_chars, fwhm=fwhm_cm,
          own_of=None):
    """The text of the separations file and its statistics.

    `char_of(wn_text, wn)` is the Sugar character of the observed line of a
    blend, or None when it is not known (treated as not wide, and counted).
    `own_of(wn_text, wn)`, if given, is the blend's factor of its own (the
    registry), or None.
    """
    first, last = make_LOPT_input.FIELDS['lower_level']
    width = last - first + 1
    out = ['; LOPT separations input file, written by lopt_separations.py.',
           '; One record per pair of components of every blend:',
           ';   (E_up1 - E_low1) - (E_up2 - E_low2) = 0 +- u',
           '; with u from the allowance tau of each component, tau = factor x',
           '; the effective line width (0.035 A instrumental with Doppler); a',
           '; line of a wide character gets wide_factor instead.',
           '; rules: factor %g wide_factor %g wide_characters %s'
           % (factor, wide_factor, ','.join(wide_chars)),
           '']
    stats = {'blends': 0, 'wide': 0, 'own': 0, 'records': 0,
             'unknown_char': 0, 'extra_dof': 0}
    for wn_text, unit, comps in blends(records):
        if unit and unit.lower() != 'cm-1':
            raise ValueError('blend at %s %s: only wavenumbers in cm-1 are '
                             'handled' % (wn_text, unit))
        wn = comps[0]['wn']
        char = char_of(wn_text, wn)
        if char is None:
            stats['unknown_char'] += 1
            char = ''
        wide = is_wide(char, wide_chars)
        own = own_of(wn_text, wn) if own_of is not None else None
        fw = fwhm(wn)
        tau = (own if own is not None else
               wide_factor if wide else factor) * fw
        out.append('; blend %s char %s %d components, FWHM %.4f cm-1, tau '
                   '%.4f%s%s' % (wn_text, char or '-', len(comps), fw, tau,
                                 ', wide' if wide else '',
                                 ', registry factor %g' % own
                                 if own is not None else ''))
        for i, j, u in pair_uncertainties([tau] * len(comps)):
            a, b = comps[i], comps[j]
            out.append('%s%s%s%s  0  %.4f'
                       % (label(a['low'], width), label(a['upp'], width),
                          label(b['low'], width), label(b['upp'], width), u))
            stats['records'] += 1
        stats['blends'] += 1
        stats['wide'] += 1 if wide else 0
        stats['own'] += 1 if own is not None else 0
        n = len(comps)
        stats['extra_dof'] += n * (n - 1) // 2 - (n - 1)
    return '\n'.join(out) + '\n', stats


def write(path, text):
    with io.open(path, 'w', encoding='ascii', newline='') as fh:
        fh.write(text)


def read_rules(path):
    """`(factor, wide_factor, wide_chars, {wn_text: char}, {wn_text: own
    factor})` from the comments of an existing file, or None if it has no
    rules line."""
    rules, chars, owns = None, {}, {}
    with io.open(path, encoding='latin-1') as fh:
        for line in fh:
            m = _HEAD_RE.match(line.strip())
            if m:
                rules = (float(m.group(1)), float(m.group(2)),
                         tuple(c for c in m.group(3).split(',') if c))
                continue
            m = _BLEND_RE.match(line)
            if m:
                chars[m.group(1)] = '' if m.group(2) == '-' else m.group(2)
                own = _OWN_RE.search(line)
                if own:
                    owns[m.group(1)] = float(own.group(1))
    if rules is None:
        return None
    return rules + (chars, owns)


def chars_from_table(path):
    """`{round(wn_obs, 3): char}` of the accepted rows of a classification
    table, for blends the separations file has not seen yet."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with io.open(path, encoding='utf-8', newline='') as fh:
        for r in csv.DictReader(fh):
            if (r.get('accepted') or '') in ('1.0', '1', 'True', 'true'):
                out.setdefault(round(float(r['wn_obs']), 3),
                               (r.get('char') or '').strip())
    return out


def refresh(lines_path, sep_path, table_chars=None, fwhm=fwhm_cm):
    """Rebuild the separations file `sep_path` for the transitions file
    `lines_path`, with the rules, characters and registry factors its own
    comments record; a blend it has no comment for takes its character from
    `table_chars` (chars_from_table) by its wavenumber, and no factor of its
    own (a registry row for it reaches the file with the next
    make_LOPT_input.py).  Returns the statistics, or None when `sep_path`
    does not exist or records no rules."""
    if not os.path.exists(sep_path):
        return None
    got = read_rules(sep_path)
    if got is None:
        return None
    factor, wide_factor, wide_chars, chars, owns = got
    table_chars = table_chars or {}

    def char_of(wn_text, wn):
        if wn_text in chars:
            return chars[wn_text]
        return table_chars.get(round(wn, 3))

    text, stats = build(read_records(lines_path), char_of, factor,
                        wide_factor, wide_chars, fwhm,
                        lambda wn_text, wn: owns.get(wn_text))
    write(sep_path, text)
    return stats


# ---------------------------------------------------------------------------
# The parameter file
# ---------------------------------------------------------------------------
def par_value(par_path, key):
    """The value of the parameter `key` in a LOPT parameter file, or None."""
    if not os.path.exists(par_path):
        return None
    with io.open(par_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            semi = line.find(';')
            if semi >= 0 and key in line[semi:].upper():
                return line[:semi].strip()
    return None


def separations_beside(lines_path, par_name='LOPT.par'):
    """The separations file the LOPT.par beside `lines_path` names, or None."""
    d = os.path.dirname(os.path.abspath(lines_path))
    name = par_value(os.path.join(d, par_name), PAR_SEP_FILE_KEY)
    if not name or name.lower() in ('none', "''"):
        return None
    return name if os.path.isabs(name) else os.path.join(d, name)


def refresh_beside(lines_path, log=None):
    """Rebuild the separations file of the set `lines_path` is in, after a
    tool has rewritten that transitions file.  The characters of new blends
    come from the set's line_classifications.csv.  Does nothing (and returns
    None) when the set's LOPT.par names no separations file."""
    sep = separations_beside(lines_path)
    if sep is None:
        return None
    table = os.path.join(os.path.dirname(os.path.abspath(lines_path)),
                         make_LOPT_input.DEF_CLASSIFICATIONS)
    stats = refresh(lines_path, sep, chars_from_table(table))
    if stats is not None and log is not None:
        log('   %s rebuilt: %d blends (%d wide), %d records%s'
            % (os.path.basename(sep), stats['blends'], stats['wide'],
               stats['records'],
               ', %d blend(s) of unknown character taken as not wide'
               % stats['unknown_char'] if stats['unknown_char'] else ''))
    return stats


def set_par_lines(par, sep_name):
    """`par` (a list of parameter-file lines) with the separations file named
    and LOPT's own separation equations switched off; the keys are replaced
    where they are and appended where they are not.  `sep_name` None removes
    both, which leaves LOPT's own equations at its default factor."""
    keys = {PAR_SEP_FILE_KEY: sep_name, PAR_SEP_FACTOR_KEY: '0'}
    out = []
    seen = set()
    for text in par:
        semi = text.find(';')
        hit = [k for k in keys if semi >= 0 and k in text[semi:].upper()]
        if hit:
            seen.add(hit[0])
            if sep_name is not None:
                out.append('%-8s%s' % (keys[hit[0]] + ' ', text[semi:]))
            continue
        out.append(text)
    if sep_name is not None:
        if PAR_SEP_FILE_KEY not in seen:
            out.append('%s ; %s' % (sep_name, PAR_SEP_FILE_KEY))
        if PAR_SEP_FACTOR_KEY not in seen:
            out.append('0       ; %s (0: the separations file holds every '
                       'blend)' % PAR_SEP_FACTOR_KEY)
    return out
