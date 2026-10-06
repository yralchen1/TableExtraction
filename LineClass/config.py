"""Loader for `lineclass_config.toml`, the configuration of the pipeline.

The configuration holds every constant that describes an input file (its name,
the worksheet to read, the names of the columns) and every constant of the
missing-gA policy, so that none of them is written in the code.  "gA" is the
statistical weight of the upper level times the transition probability, in
s^-1; it is the quantity Cowan's atomic-structure codes print, and they print
it only above a cutoff, which is what the `[icalc.completeness]` and
`[missing_gA]` sections are about.

Typical use:

    import config
    cfg = config.load()                 # the file next to this module
    cfg = config.load('other.toml')     # an explicit one

    idx = config.resolve_columns(header_row, cfg.icalc.columns, 'Icalc.xlsx')
    gA  = row[idx['gA']].value

`load` validates the file against the schema below: an unknown key, a missing
key or a value of the wrong type raises `ConfigError` and stops the run, rather
than being silently ignored.  Paths are returned already resolved against the
directory that holds the configuration file.
"""
import os
import tomllib
from dataclasses import dataclass

#: The configuration `load()` reads when it is given no path.  The
#: environment variable comes first so that a whole working set - the
#: corrected one in `iter/`, the merged one in `final/` - can be selected for
#: every program at once:
#:
#:     LINECLASS_CONFIG=iter/lineclass_config.toml python unfound_levels.py
#:
#: It selects only what the configuration names.  A program that also reads
#: IDEN2 or the LOPT output takes `--set DIR` instead (check_sync.py,
#: sync_IDEN2.py, level_positions.py, find_unknown_levels.py).
#:
#: A per-program `--config` cannot do that job.  `classify_lines.py` loads the
#: configuration when it is imported, and a dozen other programs import it, so
#: the configuration is already chosen by the time any of them parses its own
#: arguments.  The environment variable is read before the import, which is
#: why it reaches all of them.
DEFAULT_PATH = os.environ.get('LINECLASS_CONFIG') or os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'lineclass_config.toml')

#: The name of a working set's configuration file.  A directory that holds
#: one is a set, and every file below it, down to the next such directory,
#: belongs to that set (see set_of).
CONFIG_NAME = 'lineclass_config.toml'

#: The override of a lock, for a program run by another: insert_new_level.py,
#: move_level.py and discard_level.py, given --unlock, set it to 1 for the
#: programs of the chain they run.  See require_unlocked.
UNLOCK_ENV = 'LINECLASS_UNLOCK'


class ConfigError(Exception):
    """Raised on any problem with the configuration file."""


# --- schema -----------------------------------------------------------------
# Leaves of the schema are either a python type, or one of the two sentinels.
class _StrMap:
    """A table whose keys and values are all strings (a column mapping)."""


class _FloatPair:
    """A list of exactly two numbers."""


class _StrList:
    """A list of strings, possibly empty."""


_ENUMS = {
    'missing_gA.policy': ('none', 'impute'),
    'missing_gA.u_ln_estimator': ('rms', 'mean', 'median'),
    'hfs.iden2_display': ('measured', 'lopt'),
}

#: The classes of the hfs convention model, each given as
#: [kappa, its uncertainty] under [hfs.kappa]; see hfs_correction.py.
#: `flag_unlisted` is needed only by a configuration that names the list of
#: Sugar's resolved components (files.hfs_components): it is the class of the
#: flagged lines whose components he did not list (2026-10-06).
HFS_KAPPA_CLASSES = ('flag', 'plain_1974', 'plain_1969', 'c', 'flag_unlisted')

# Keys the configuration file may leave out; they take the default written
# into load() below.
_OPTIONAL = {'inherit', 'locked',
             'files.level_overrides', 'files.line_decisions',
             'files.new_levels', 'files.icalc_extra',
             'files.discarded_levels', 'files.inflated_unc',
             'files.hfs_A_levels', 'files.hfs_satellites',
             'files.kappa_exceptions', 'files.hfs_components',
             'hfs.kappa.flag_unlisted',
             'decisions', 'decisions.max_forced_offset',
             'hfs', 'hfs.resolved_levels', 'hfs.iden2_display'}

_SCHEMA = {
    'inherit': str,
    'locked': bool,
    'files': {'levels': str, 'lines': str, 'icalc': str,
              'output': str, 'output_csv': str, 'level_overrides': str,
              'line_decisions': str, 'new_levels': str,
              'icalc_extra': str, 'discarded_levels': str,
              'inflated_unc': str, 'hfs_A_levels': str,
              'hfs_satellites': str, 'kappa_exceptions': str,
              'hfs_components': str},
    'range': {'wn_min': float, 'wn_max': float},
    'levels': {'layout': {'sheet': str, 'columns': _StrMap}},
    'lines': {'layout': {'sheet': str, 'columns': _StrMap}},
    'icalc': {'layout': {'sheet': str, 'columns': _StrMap},
              'completeness': {'gA_cutoff': float, 'allow_below_cutoff': bool}},
    'missing_gA': {'policy': str, 'u_ln_window': float, 'u_ln_estimator': str,
                   'self_consistent': bool, 'fit_range_decades': _FloatPair},
    'intensity_model': {'C': float, 'kT': float, 'verify_tolerance': float},
    'decisions': {'max_forced_offset': float},
    'hfs': {'apply': bool, 'resolved_levels': _StrList, 'iden2_display': str,
            'kappa': {k: _FloatPair for k in HFS_KAPPA_CLASSES}},
}


def _validate(node, schema, path=''):
    if not isinstance(node, dict):
        raise ConfigError(f"'{path.rstrip('.')}' must be a table")
    for key in node:
        if key not in schema:
            raise ConfigError(f"unknown key '{path}{key}'; allowed here: "
                              + ', '.join(sorted(schema)))
    for key, expected in schema.items():
        full = path + key
        if key not in node:
            if full in _OPTIONAL:
                continue
            raise ConfigError(f"missing key '{full}'")
        val = node[key]
        if isinstance(expected, dict):
            _validate(val, expected, full + '.')
        elif expected is _StrMap:
            if not isinstance(val, dict) or not all(
                    isinstance(k, str) and isinstance(v, str)
                    for k, v in val.items()):
                raise ConfigError(
                    f"'{full}' must be a table of  name = \"column\"  pairs")
        elif expected is _FloatPair:
            if (not isinstance(val, list) or len(val) != 2
                    or not all(isinstance(v, (int, float))
                               and not isinstance(v, bool) for v in val)):
                raise ConfigError(f"'{full}' must be a list of two numbers")
            node[key] = [float(v) for v in val]
        elif expected is _StrList:
            if (not isinstance(val, list)
                    or not all(isinstance(v, str) for v in val)):
                raise ConfigError(f"'{full}' must be a list of strings")
        elif expected is bool:
            if not isinstance(val, bool):
                raise ConfigError(f"'{full}' must be true or false")
        elif expected is float:
            if isinstance(val, bool) or not isinstance(val, (int, float)):
                raise ConfigError(f"'{full}' must be a number")
            node[key] = float(val)
        elif expected is str:
            if not isinstance(val, str):
                raise ConfigError(f"'{full}' must be a string")
            allowed = _ENUMS.get(full)
            if allowed is not None and val not in allowed:
                raise ConfigError(f"'{full}' must be one of "
                                  + ', '.join(repr(a) for a in allowed)
                                  + f", not {val!r}")


# --- typed view -------------------------------------------------------------
@dataclass(frozen=True)
class Layout:
    """Where to find a table: its worksheet and its columns.

    `sheet` is the worksheet name; `columns` maps the logical name used in the
    code to the column name written in the header row of that worksheet.
    """
    sheet: str
    columns: dict


@dataclass(frozen=True)
class HfsSettings:
    """The [hfs] section: the head-frame hyperfine correction.

    `apply` switches it on; with it off (the default, and the state of a
    configuration without the section) no program changes anything.
    `A_levels` is the table of hyperfine A constants (files.hfs_A_levels,
    '' if none), `kappa` maps each class of HFS_KAPPA_CLASSES the file
    gives to (kappa, its uncertainty), and `resolved_levels` names the levels whose
    sublevels are resolved and whose own A therefore counts as zero in a
    line's displacement.  `satellites` is the registry of resolved hfs
    companions (files.hfs_satellites, '' if none), which is read whether
    `apply` is on or not.  `kappa_exceptions` is the registry of lines
    measured otherwise than their class says (files.kappa_exceptions, ''
    if none), read with `apply` on and by the plate calibration.
    `components` is the list of Sugar's resolved hfs components
    (files.hfs_components, '' if none): with it, a flagged line whose
    components he did not list is of the class `flag_unlisted`, which
    [hfs.kappa] must then give.  `iden2_display` says where sync_IDEN2.py puts the
    lines LOPT uses in IDEN2's dlv.dat: 'measured' (the default) or 'lopt',
    at the wavenumber and uncertainty LOPT is given, which needs `apply`.
    hfs_correction.py says what each of them means.
    """
    apply: bool = False
    A_levels: str = ''
    satellites: str = ''
    kappa_exceptions: str = ''
    components: str = ''
    kappa: tuple = ()           # ((class, (kappa, u_kappa)), ...)
    resolved_levels: tuple = ()
    iden2_display: str = 'measured'

    def kappa_of(self, cls: str) -> tuple:
        """(kappa, u_kappa) of the class `cls`."""
        return dict(self.kappa)[cls]


@dataclass(frozen=True)
class Config:
    path: str                 # the configuration file this was read from
    levels_file: str          # workbook of the adopted energy levels
    lines_file: str           # workbook of the observed lines
    icalc_file: str           # workbook of the calculated transitions
    output_file: str          # output workbook
    output_csv: str           # the same table as csv
    level_overrides: str      # csv of revised adopted energies, '' if none
    line_decisions: str       # csv of manual accept/reject verdicts, '' if none
    new_levels: str           # csv of levels found since the level list, '' if none
    icalc_extra: str          # workbook of their calculated transitions, '' if none
    discarded_levels: str     # csv of levels whose position was given up, '' if none
    inflated_unc: str         # tab-delimited registry of hand-set line uncertainties, '' if none
    wn_min: float             # lower end of the observed range, cm^-1
    wn_max: float             # upper end of the observed range, cm^-1
    levels: Layout
    lines: Layout
    icalc: Layout
    gA_cutoff: float          # gA below which Cowan's codes printed nothing, s^-1
    allow_below_cutoff: bool  # read rows that fall marginally below it
    missing_gA: dict          # the [missing_gA] table, keys as in the file
    intensity_model: dict     # the [intensity_model] table
    max_forced_offset: float = 5.0
    # How far, in cm^-1, the Ritz wavenumber of a pair the decision ledger
    # accepts may sit from the wavenumber of the line it is accepted on before
    # the run stops.  See classify_lines.check_forced_decisions().
    hfs: HfsSettings = HfsSettings()


def _merge(base_tbl: dict, over: dict) -> dict:
    """`over` laid on top of `base_tbl`, table by table.

    A table is merged key by key, so a working set can name the two columns it
    reads its wavenumbers from without repeating the other four; anything that
    is not a table is replaced outright.
    """
    out = dict(base_tbl)
    for key, val in over.items():
        if isinstance(val, dict) and isinstance(out.get(key), dict):
            out[key] = _merge(out[key], val)
        else:
            out[key] = val
    return out


def _read(path: str, seen=()) -> dict:
    """The raw table of `path`, with anything it inherits already under it.

    `inherit` names another configuration file, relative to this one.  The
    working sets - the calibrated one in `iter/`, the merged one in `final/` -
    are the whole point of it: each differs from the baseline in a handful of
    settings, and a copy of the file would have to be kept in step with the
    baseline by hand, which is how the intensity model of one set comes to be
    a year older than the other's.

    The inherited file's own paths are resolved to absolute names before the
    merge, against the directory holding *it*; the inheriting file's are
    resolved later against the directory holding itself.  So `levels =
    "../TableExtraction/..."` in the baseline goes on naming the same workbook
    when `iter/lineclass_config.toml` inherits it, and `lines =
    "Pr3_lines_corrected.xlsx"` written in the set names the set's own.
    """
    if path in seen:
        raise ConfigError("configuration inherits itself: "
                          + ' -> '.join(seen + (path,)))
    if not os.path.isfile(path):
        raise ConfigError(f"configuration file not found: {path}")
    with open(path, 'rb') as fh:
        try:
            raw = tomllib.load(fh)
        except tomllib.TOMLDecodeError as exc:
            raise ConfigError(f"{path}: {exc}") from None
    parent_name = raw.pop('inherit', None)
    if parent_name is None:
        return raw
    if not isinstance(parent_name, str):
        raise ConfigError(f"{path}: 'inherit' must be a string")
    here = os.path.dirname(path)
    parent = _read(os.path.abspath(os.path.join(here, parent_name)),
                   seen + (path,))
    p_base = os.path.dirname(os.path.abspath(os.path.join(here, parent_name)))
    # A lock belongs to the directory whose file says it, not to the sets
    # that inherit from that file: iter/ inherits the locked baseline and is
    # where the work is done.
    parent.pop('locked', None)
    for name, value in (parent.get('files') or {}).items():
        if isinstance(value, str) and value:
            parent['files'][name] = os.path.normpath(
                os.path.join(p_base, value))
    return _merge(parent, raw)


#: the classes of [hfs.kappa] the plate calibration fits; `flag` is its anchor
#: and is never rewritten.  `flag_unlisted` is fitted only where the
#: configuration splits the flags (files.hfs_components); the other three
#: always are.
FITTED_KAPPA_CLASSES = ('plain_1974', 'plain_1969', 'c', 'flag_unlisted')
REQUIRED_KAPPA_CLASSES = ('plain_1974', 'plain_1969', 'c')

#: what begins each comment line `set_hfs_kappa` writes, so that the next run
#: can find and replace them; nothing else in the file is touched.
KAPPA_NOTE = '# [wavelength_calibration.py] '


def kappa_file(path: str = None) -> str:
    """The configuration file that holds the `[hfs.kappa]` a run of `path`
    uses: `path` itself if it writes the table, else the nearest file it
    inherits from that does (today the baseline's, which iter_hfs/
    inherits).  Raises if none of them does."""
    path = os.path.abspath(path or DEFAULT_PATH)
    seen = ()
    while True:
        if path in seen:
            raise ConfigError("configuration inherits itself: "
                              + ' -> '.join(seen + (path,)))
        if not os.path.isfile(path):
            raise ConfigError(f"configuration file not found: {path}")
        with open(path, 'rb') as fh:
            raw = tomllib.load(fh)
        if 'kappa' in (raw.get('hfs') or {}):
            return path
        parent = raw.get('inherit')
        if not parent:
            first = seen[0] if seen else path
            raise ConfigError(f"no file of the chain of {first} has an "
                              f"[hfs.kappa] table")
        seen += (path,)
        path = os.path.abspath(os.path.join(os.path.dirname(path), parent))


def set_hfs_kappa(path: str, values: dict, note: str = '') -> tuple:
    """Write the fitted kappas `values` ({class: (kappa, u_kappa)}, the
    classes of FITTED_KAPPA_CLASSES) into the `[hfs.kappa]` the configuration
    `path` uses (`kappa_file`), and return `(file written, the values it
    held before)`.

    The classes of REQUIRED_KAPPA_CLASSES must be in `values`; those of
    FITTED_KAPPA_CLASSES beyond them are written when `values` gives them.
    The file is edited as text, so that every comment and every other line
    stays as it is: only the line of each class is rewritten, in the form
    `plain_1974 = [0.850, 0.013]`, and the comment lines beginning with
    KAPPA_NOTE are replaced by `note`, placed just above the `flag` line.
    The file is read back afterwards; if it does not give exactly the new
    values, the old text is put back and ConfigError raised.
    """
    missing = [c for c in REQUIRED_KAPPA_CLASSES if c not in values]
    wanted = [c for c in FITTED_KAPPA_CLASSES if c in values]
    if missing:
        raise ConfigError(f"set_hfs_kappa: no value for "
                          f"{', '.join(missing)}")
    target = kappa_file(path)
    with open(target, 'rb') as fh:
        before = {k: tuple(v) for k, v in
                  tomllib.load(fh)['hfs']['kappa'].items()}
    with open(target, encoding='utf-8', newline='') as fh:
        text = fh.read()
    eol = '\r\n' if '\r\n' in text else '\n'
    lines = text.split(eol)
    try:
        start = next(i for i, ln in enumerate(lines)
                     if ln.strip() == '[hfs.kappa]')
    except StopIteration:
        raise ConfigError(f"{target}: [hfs.kappa] is not written as a "
                          f"table header of its own") from None
    end = next((i for i in range(start + 1, len(lines))
                if lines[i].lstrip().startswith('[')), len(lines))
    # the names are padded as the table's own keys are, so that a file keeps
    # its alignment whichever classes it holds
    keys = [ln.split('=', 1)[0].strip() for ln in lines[start + 1:end]
            if '=' in ln and not ln.lstrip().startswith('#')]
    width = max([len(k) for k in keys] + [len(c) for c in wanted])
    body, done = [], set()
    for ln in lines[start + 1:end]:
        if ln.startswith(KAPPA_NOTE):
            continue
        key = (ln.split('=', 1)[0].strip()
               if '=' in ln and not ln.lstrip().startswith('#') else None)
        if key == 'flag' and note:
            body += [KAPPA_NOTE + part for part in note.splitlines()]
        if key in wanted:
            k, u = values[key]
            ln = '%s = [%.3f, %.3f]' % (key.ljust(width), k, u)
            done.add(key)
        body.append(ln)
    if done != set(wanted):
        raise ConfigError(f"{target}: [hfs.kappa] lacks the line of "
                          f"{', '.join(sorted(set(wanted) - done))}")
    new_text = eol.join(lines[:start + 1] + body + lines[end:])
    with open(target, 'w', encoding='utf-8', newline='') as fh:
        fh.write(new_text)
    try:
        with open(target, 'rb') as fh:
            after = tomllib.load(fh)['hfs']['kappa']
        for c in wanted:
            want = tuple(float('%.3f' % x) for x in values[c])
            if tuple(after[c]) != want:
                raise ConfigError(f"{target}: [hfs.kappa] {c} reads back as "
                                  f"{after[c]}, not {list(want)}")
    except (tomllib.TOMLDecodeError, ConfigError, KeyError) as exc:
        with open(target, 'w', encoding='utf-8', newline='') as fh:
            fh.write(text)
        raise ConfigError(f"set_hfs_kappa: {exc}; {target} "
                          f"restored") from None
    return target, before


def load(path: str = None) -> Config:
    """Read, validate and return the configuration."""
    path = os.path.abspath(path or DEFAULT_PATH)
    raw = _read(path)
    _validate(raw, _SCHEMA)

    base = os.path.dirname(path)

    def _p(name):
        return os.path.normpath(os.path.join(base, raw['files'][name]))

    satellites = (_p('hfs_satellites')
                  if raw['files'].get('hfs_satellites') else '')
    exceptions = (_p('kappa_exceptions')
                  if raw['files'].get('kappa_exceptions') else '')
    components = (_p('hfs_components')
                  if raw['files'].get('hfs_components') else '')
    hfs = HfsSettings(satellites=satellites, kappa_exceptions=exceptions,
                      components=components)
    if 'hfs' in raw:
        h = raw['hfs']
        a_levels = (_p('hfs_A_levels')
                    if raw['files'].get('hfs_A_levels') else '')
        if h['apply'] and not a_levels:
            raise ConfigError("[hfs] apply = true needs files.hfs_A_levels, "
                              "the table of hyperfine A constants")
        display = h.get('iden2_display', 'measured')
        if display == 'lopt' and not h['apply']:
            raise ConfigError("[hfs] iden2_display = 'lopt' needs apply = "
                              "true: with the correction off LOPT is given "
                              "the measured lines, so there is nothing else "
                              "to show")
        if components and 'flag_unlisted' not in h['kappa']:
            raise ConfigError("files.hfs_components splits the flagged lines, "
                              "so [hfs.kappa] needs flag_unlisted, the class "
                              "of those whose components Sugar did not list")
        hfs = HfsSettings(
            apply=h['apply'], A_levels=a_levels, satellites=satellites,
            kappa_exceptions=exceptions, components=components,
            kappa=tuple((k, tuple(h['kappa'][k]))
                        for k in HFS_KAPPA_CLASSES if k in h['kappa']),
            resolved_levels=tuple(h.get('resolved_levels', ())),
            iden2_display=display)

    return Config(
        path=path,
        levels_file=_p('levels'),
        lines_file=_p('lines'),
        icalc_file=_p('icalc'),
        output_file=_p('output'),
        output_csv=_p('output_csv'),
        level_overrides=(_p('level_overrides')
                         if raw['files'].get('level_overrides') else ''),
        line_decisions=(_p('line_decisions')
                        if raw['files'].get('line_decisions') else ''),
        new_levels=(_p('new_levels')
                    if raw['files'].get('new_levels') else ''),
        icalc_extra=(_p('icalc_extra')
                     if raw['files'].get('icalc_extra') else ''),
        discarded_levels=(_p('discarded_levels')
                          if raw['files'].get('discarded_levels') else ''),
        inflated_unc=(_p('inflated_unc')
                      if raw['files'].get('inflated_unc') else ''),
        wn_min=raw['range']['wn_min'],
        wn_max=raw['range']['wn_max'],
        levels=Layout(**raw['levels']['layout']),
        lines=Layout(**raw['lines']['layout']),
        icalc=Layout(**raw['icalc']['layout']),
        gA_cutoff=raw['icalc']['completeness']['gA_cutoff'],
        allow_below_cutoff=raw['icalc']['completeness']['allow_below_cutoff'],
        missing_gA=raw['missing_gA'],
        intensity_model=raw['intensity_model'],
        max_forced_offset=float(
            raw.get('decisions', {}).get('max_forced_offset', 5.0)),
        hfs=hfs,
    )


def resolve_columns(header_row, columns: dict, source: str) -> dict:
    """Map the logical column names to 0-based positions in a worksheet.

    `header_row` is the first row of the worksheet as openpyxl returns it (a
    tuple of cells); `columns` is a Layout.columns mapping, logical name ->
    column name as written in that header; `source` names the file, for the
    error message.  Returns {logical name: index}.  A column named in the
    configuration but absent from the header stops the run: silently reading
    the wrong column would be far worse than failing here.
    """
    header = {}
    for i, cell in enumerate(header_row):
        name = cell.value
        if name is None:
            continue
        name = str(name).strip()
        header.setdefault(name, i)   # on duplicate names the first one wins
    idx = {}
    for logical, physical in columns.items():
        if physical not in header:
            raise ConfigError(
                f"{source}: no column named {physical!r} (wanted for "
                f"{logical!r}); the header row has: "
                + ', '.join(repr(h) for h in header))
        idx[logical] = header[physical]
    return idx


# --- the lock on a working set ------------------------------------------------
# A set is locked by `locked = true` at the top of its own configuration file.
# The baseline (LineClass/ itself) is locked so that no run made by habit,
# without --set or --config, overwrites the files its LOPT fit and its IDEN2
# were last brought into step with.  The lock is about the files a run
# writes, not about the configuration it reads: a set whose configuration
# forgot to name its own output, or which has no IDEN2 of its own and so
# falls back on the project's, would write the baseline's files, and it is
# the files that are checked.  The hand-kept files - line_decisions.csv,
# new_levels.txt and the rest - live in the baseline's directory and are
# shared by every set; the programs that write them do not ask this question
# of them.
def set_of(path: str):
    """The working set `path` belongs to: the nearest directory, from the
    file's own upward, that holds a lineclass_config.toml.  None if no
    directory above it does (a scratch copy, a test's temporary folder)."""
    d = os.path.dirname(os.path.abspath(path))
    while True:
        if os.path.isfile(os.path.join(d, CONFIG_NAME)):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            return None
        d = parent


def is_locked(set_dir: str) -> bool:
    """True if the set's own configuration file says `locked = true`.

    Only that file is read.  A configuration that inherits a locked one is
    not locked by it, which is what lets iter/ inherit the baseline.
    """
    path = os.path.join(set_dir, CONFIG_NAME)
    try:
        with open(path, 'rb') as fh:
            raw = tomllib.load(fh)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise ConfigError(f"{path}: {exc}") from None
    value = raw.get('locked', False)
    if not isinstance(value, bool):
        raise ConfigError(f"{path}: 'locked' must be true or false")
    return value


def require_unlocked(paths, program: str, unlock: bool = False) -> list:
    """Stop the run if any of `paths` belongs to a locked set.

    `unlock` (the program's --unlock), or the environment variable
    UNLOCK_ENV set to 1, lets it through; the locked sets it then writes are
    printed and returned, so that the run says what it was allowed to do.
    Called before any work is done, as output_files.require_writable is.
    """
    locked = {}
    for path in paths:
        home = set_of(path)
        if home is not None and is_locked(home):
            locked.setdefault(home, []).append(os.path.abspath(path))
    if not locked:
        return []
    if unlock or os.environ.get(UNLOCK_ENV) == '1':
        for home in sorted(locked):
            print(f'{program}: writing the locked set {home} (--unlock)')
        return sorted(locked)
    lines = ['', f'{program}: this run would write a locked set.']
    for home, files in sorted(locked.items()):
        lines.append(f'  {home}   ({CONFIG_NAME} there has locked = true)')
        lines += [f'    {f}' for f in files]
    lines += ['Run it on a working set instead, or pass --unlock to write '
              'this one anyway.', '']
    raise SystemExit('\n'.join(lines))


def require_own_output(config_path: str, paths, program: str) -> None:
    """Stop the run if any of `paths` belongs to another set than the
    configuration `config_path` does.

    A set's configuration inherits its parent's [files] table resolved
    against the parent's directory, so a set that does not name its own
    output inherits its parent's table and would overwrite it: an iter_hfs/
    without `output` and `output_csv` would rewrite iter's table in the head
    frame.  The lock does not catch that, since iter/ is not locked.  A file
    belongs to the set of the nearest directory above it that holds a
    lineclass_config.toml (set_of); the configuration file to the set of its
    own directory, or the nearest above it.  There is no override: a table
    meant for another set is written by that set's configuration.
    """
    def key(d):
        return None if d is None else os.path.normcase(d)

    own = set_of(config_path)
    stray = [os.path.abspath(f) for f in paths if key(set_of(f)) != key(own)]
    if not stray:
        return
    lines = ['', f'{program}: this configuration would write another set\'s '
             'files.', f'  configuration: {os.path.abspath(config_path)}',
             f'  its set:       {own}']
    for f in stray:
        lines.append(f'    {f}   (belongs to {set_of(f)})')
    lines += ['Name this set\'s own files in its [files] table, e.g.',
              '  output     = "line_classifications.xlsx"',
              '  output_csv = "line_classifications.csv"', '']
    raise SystemExit('\n'.join(lines))
