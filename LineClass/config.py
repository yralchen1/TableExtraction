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
#:     LINECLASS_CONFIG=iter/lineclass_config.toml python level_positions.py
#:
#: A per-program `--config` cannot do that job.  `classify_lines.py` loads the
#: configuration when it is imported, and a dozen other programs import it, so
#: the configuration is already chosen by the time any of them parses its own
#: arguments.  The environment variable is read before the import, which is
#: why it reaches all of them.
DEFAULT_PATH = os.environ.get('LINECLASS_CONFIG') or os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'lineclass_config.toml')


class ConfigError(Exception):
    """Raised on any problem with the configuration file."""


# --- schema -----------------------------------------------------------------
# Leaves of the schema are either a python type, or one of the two sentinels.
class _StrMap:
    """A table whose keys and values are all strings (a column mapping)."""


class _FloatPair:
    """A list of exactly two numbers."""


_ENUMS = {
    'missing_gA.policy': ('none', 'impute'),
    'missing_gA.u_ln_estimator': ('rms', 'mean', 'median'),
}

# Keys the configuration file may leave out; they take the default written
# into load() below.
_OPTIONAL = {'inherit',
             'files.level_overrides', 'files.line_decisions',
             'files.new_levels', 'files.icalc_extra',
             'files.discarded_levels',
             'decisions', 'decisions.max_forced_offset'}

_SCHEMA = {
    'inherit': str,
    'files': {'levels': str, 'lines': str, 'icalc': str,
              'output': str, 'output_csv': str, 'level_overrides': str,
              'line_decisions': str, 'new_levels': str,
              'icalc_extra': str, 'discarded_levels': str},
    'range': {'wn_min': float, 'wn_max': float},
    'levels': {'layout': {'sheet': str, 'columns': _StrMap}},
    'lines': {'layout': {'sheet': str, 'columns': _StrMap}},
    'icalc': {'layout': {'sheet': str, 'columns': _StrMap},
              'completeness': {'gA_cutoff': float, 'allow_below_cutoff': bool}},
    'missing_gA': {'policy': str, 'u_ln_window': float, 'u_ln_estimator': str,
                   'self_consistent': bool, 'fit_range_decades': _FloatPair},
    'intensity_model': {'C': float, 'kT': float, 'verify_tolerance': float},
    'decisions': {'max_forced_offset': float},
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
    for name, value in (parent.get('files') or {}).items():
        if isinstance(value, str) and value:
            parent['files'][name] = os.path.normpath(
                os.path.join(p_base, value))
    return _merge(parent, raw)


def load(path: str = None) -> Config:
    """Read, validate and return the configuration."""
    path = os.path.abspath(path or DEFAULT_PATH)
    raw = _read(path)
    _validate(raw, _SCHEMA)

    base = os.path.dirname(path)

    def _p(name):
        return os.path.normpath(os.path.join(base, raw['files'][name]))

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
