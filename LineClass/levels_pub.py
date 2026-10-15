"""Extend levels_pub.xlsx with Set 3's level energies and uncertainties.

levels_pub.xlsx is the table of the article: one row per calculated level of
the Cowan calculation, with its composition, its calculated energy Ec and,
where it is known, its published energy (Sugar's or Wyart's) and level_id.
This tool completes its level_id column and writes seven columns after the
existing ones:

  lid       the level's number in the Cowan calculation (tp_E1_no_trials.xlsx)
  Eopt      the optimized energy, cm^-1, from final/LOPT_output_levels.txt
  unc_stat  max(D1, D2stat), 0 for the fixed ground level
  unc_sys   sqrt(u_cal^2 + u_hfs^2 + u_kappa^2), the systematic parts that
            lopt_perturb.py measures (calibration, A constants, kappa)
  unc_tot   sqrt(unc_stat^2 + unc_sys^2)
  N_lines   the number of lines that determine the level: LOPT's count
            without its "L" and "D" parts
  Notes     N new (new_levels.txt), R revised (revised_level_energies.csv),
            D discarded (discarded_levels.csv); letters combine, e.g. "NR"

The uncertainties come from final/levels_set3.csv, which make_set3.py
--levels writes from the set's LOPT fit and the two lopt_perturb.py passes
(user, 2026-10-11).  Its unc_tot is that table's D2tot*.  So the order after
any change to the calibration, the uncertainty model or kappa is: the chains
of iter_hfs and final, lopt_perturb.py twice, make_set3.py --levels, and then
this tool.  It stops if levels_set3.csv does not belong to the fit in
final/LOPT_output_levels.txt (another energy or line count), or has no
systematic parts.

Each row is joined to the calculated level list by its level_id, through
classify_lines.cowan_lid_ids (the IDEN2 table, never by energy).  A row
without a level_id is joined by parity, J and Ec, the calculated energy,
which never changes; Ec is printed to 0.1 cm^-1 in the workbook.  A level of
the fit that ends up in no row is reported and the run stops: the table must
hold every level.

The workbook is rewritten in place unless --out names another file.  The
seven columns are found by their headers and replaced, so a rerun after the
next calibration overwrites them.  Nothing else in the workbook is changed
except the empty level_id cells that the join fills.

The report also lists, for the information of the analyst, the levels whose
Eopt departs most from the published energy without being marked R.  "R" is
taken from revised_level_energies.csv only: a revision means a different set
of lines, which a large shift by itself does not prove.

Usage (from LineClass/):
    python levels_pub.py                  # update levels_pub.xlsx
    python levels_pub.py --check          # report only, write nothing
    python levels_pub.py --out levels_pub_new.xlsx
"""
import argparse
import csv
import math
import os
import sys

import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import config                              # noqa: E402
import output_files                        # noqa: E402

PUB = os.path.join(HERE, 'levels_pub.xlsx')
SET_DIR = os.path.join(HERE, 'final')
LOPT_LEVELS = 'LOPT_output_levels.txt'
SET3_LEVELS = 'levels_set3.csv'
TOTAL = 'D2tot*'

NEW_COLUMNS = ['lid', 'Eopt', 'unc_stat', 'unc_sys', 'unc_tot', 'N_lines',
               'Notes']
#: the workbook's Ec is not the E_calc of tp_E1_no_trials.xlsx to the
#: digit: the median difference is 0.1-0.2 cm^-1 and a few reach 5 (measured
#: 2026-10-15), so a row without a level_id is joined to the nearest level of
#: its parity and J within this window, and only if no other comes close
EC_TOL = 6.0
#: how many unmarked levels with the largest |Eopt - energy| the report lists
N_SHIFTS = 15
PARITY = {'o': 'odd', 'e': 'even', 'odd': 'odd', 'even': 'even'}


def norm_id(value):
    """'059003.000624' from a cell, '' for an empty one."""
    return '' if value is None else str(value).strip()


# ---------------------------------------------------------------------------
# The calculated levels
# ---------------------------------------------------------------------------
def calculated_levels():
    """`(calc, lid_of_id)`: the Cowan level list as records (lid, parity, J,
    E_calc) and {level_id: lid} for every level that has an identifier."""
    import classify_lines
    import cowan_gA
    trans = cowan_gA.read_transitions(log=lambda m: None)
    table = cowan_gA.levels(trans)
    calc = [dict(lid=int(r['lid']), parity=str(r['parity']),
                 J=float(r['J']), E_calc=float(r['E_calc']),
                 conf=str(r['conf']), term=str(r['term']))
            for r in table.to_dict('records')]
    lid_of_id = {v: k for k, v in classify_lines.cowan_lid_ids(trans).items()}
    return calc, lid_of_id


def match_rows(rows, calc, lid_of_id, tol=EC_TOL):
    """The Cowan lid of each workbook row.

    `rows` are dicts with parity ('o'/'e'), J, Ec and level_id.  Returns
    `(lids, problems)`: one lid or None per row, and the report lines of the
    rows that could not be joined or whose join looks wrong.
    """
    by_lid = {c['lid']: c for c in calc}
    lids, problems, taken = [None] * len(rows), [], set()
    # first the rows that name their level
    for i, r in enumerate(rows):
        lid = lid_of_id.get(r['level_id']) if r['level_id'] else None
        if r['level_id'] and lid is None:
            problems.append('row %d: level_id %s is not in the calculated '
                            'level list' % (r['row'], r['level_id']))
            continue
        if lid is None:
            continue
        c = by_lid[lid]
        if (c['parity'] != PARITY.get(r['parity'])
                or abs(c['J'] - r['J']) > 1e-6):
            problems.append('row %d: %s is lid %d (%s J=%g), the row says %s '
                            'J=%g' % (r['row'], r['level_id'], lid,
                                      c['parity'], c['J'], r['parity'],
                                      r['J']))
            continue
        if abs(c['E_calc'] - r['Ec']) > tol:
            problems.append('row %d: %s is lid %d with E_calc %.1f, the '
                            'row has Ec %.1f (kept: the identifier decides)'
                            % (r['row'], r['level_id'], lid, c['E_calc'],
                               r['Ec']))
        lids[i] = lid
        taken.add(lid)
    # then the rest, by parity, J and the calculated energy
    pending = []
    for i, r in enumerate(rows):
        if lids[i] is not None or r['level_id']:
            continue
        cand = sorted((c for c in calc if c['lid'] not in taken
                       and c['parity'] == PARITY.get(r['parity'])
                       and abs(c['J'] - r['J']) < 1e-6
                       and abs(c['E_calc'] - r['Ec']) <= tol),
                      key=lambda c: abs(c['E_calc'] - r['Ec']))
        # the nearest must be clearly the nearest, or the row is reported
        clear = len(cand) == 1 or (
            len(cand) > 1 and abs(cand[1]['E_calc'] - r['Ec'])
            >= 2 * abs(cand[0]['E_calc'] - r['Ec']) + 0.2)
        if not clear:
            pending.append((i, cand))
            continue
        lids[i] = cand[0]['lid']
        taken.add(cand[0]['lid'])
    # a row near two levels: the table is in lid order within each parity,
    # so the one between the lids of its neighbors is the one
    for i, cand in pending:
        r = rows[i]
        same = [(rows[k]['row'], lids[k]) for k in range(len(rows))
                if lids[k] is not None and rows[k]['parity'] == r['parity']]
        lo = max((l for row, l in same if row < r['row']), default=-1)
        hi = min((l for row, l in same if row > r['row']), default=10 ** 9)
        fit = [c for c in cand if lo < c['lid'] < hi and c['lid'] not in taken]
        if len(fit) != 1:
            problems.append('row %d: %s J=%g Ec %.1f matches %d calculated '
                            'levels' % (r['row'], r['parity'], r['J'],
                                        r['Ec'], len(cand)))
            continue
        lids[i] = fit[0]['lid']
        taken.add(fit[0]['lid'])
    return lids, problems


# ---------------------------------------------------------------------------
# The set's levels
# ---------------------------------------------------------------------------
def read_tsv(path):
    with open(path, encoding='utf-8', newline='') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    return [{k.strip(): (v or '').strip() for k, v in r.items()}
            for r in rows]


def read_csv(path):
    with open(path, encoding='utf-8', newline='') as fh:
        return [{k.strip(): (v or '').strip() for k, v in r.items()}
                for r in csv.DictReader(fh)]


def level_values(lopt_rows, set3_rows):
    """{level_id: (Eopt, unc_stat, unc_sys, unc_tot, N_lines, S)}, S the
    distance of the level's hfs head from its center of gravity.

    `lopt_rows` are LOPT_output_levels.txt's records, `set3_rows`
    levels_set3.csv's.  Raises SystemExit when the two are not of one fit or
    the systematic parts are missing.
    """
    set3 = {r['level_id']: r for r in set3_rows}
    out, stale = {}, []
    for r in lopt_rows:
        lid = r['Designation']
        s = set3.get(lid)
        n_lines = r['N_lines'].split(',')[0].strip()
        if (s is None or abs(float(s['E_cg']) - float(r['Energy'])) > 1e-6
                or s['N_lines'] != n_lines):
            stale.append(lid)
            continue
        if not s.get(TOTAL) or not s.get('u_cal') or not s.get('u_hfs'):
            raise SystemExit('%s has no systematic parts for %s: run '
                             'lopt_perturb.py twice and make_set3.py --levels'
                             % (SET3_LEVELS, lid))
        u_stat = float(s['u_stat'])
        u_sys = math.sqrt(sum(float(s[k] or 0.0) ** 2
                              for k in ('u_cal', 'u_hfs', 'u_kappa')))
        u_tot = math.hypot(u_stat, u_sys)
        if abs(u_tot - float(s[TOTAL])) > 2e-7:
            raise SystemExit('%s: unc_tot %.7f differs from %s %s'
                             % (lid, u_tot, TOTAL, s[TOTAL]))
        out[lid] = (float(r['Energy']), u_stat, u_sys, u_tot, int(n_lines),
                    float(s.get('S') or 0.0))
    extra = sorted(set(set3) - {r['Designation'] for r in lopt_rows})
    if stale or extra:
        raise SystemExit(
            '%s does not belong to the fit in %s (%d levels differ in energy '
            'or line count, e.g. %s; %d not in the fit): run lopt_perturb.py '
            'twice and make_set3.py --levels on the present fit first'
            % (SET3_LEVELS, LOPT_LEVELS, len(stale),
               ', '.join(stale[:3]) or '-', len(extra)))
    return out


def notes_of(new_ids, revised_ids, discarded_ids):
    """{level_id: 'N'/'R'/'D' letters, in that order}."""
    out = {}
    for letter, ids in (('N', new_ids), ('R', revised_ids),
                        ('D', discarded_ids)):
        for i in ids:
            out[i] = out.get(i, '') + letter
    return out


def id_column(path, name='level_id', sep=','):
    if not path or not os.path.exists(path):
        return set()
    with open(path, encoding='utf-8', newline='') as fh:
        return {(r.get(name) or '').strip()
                for r in csv.DictReader(fh, delimiter=sep)} - {''}


# ---------------------------------------------------------------------------
# The workbook
# ---------------------------------------------------------------------------
def header_map(ws):
    return {str(c.value).strip(): c.column for c in ws[1]
            if c.value is not None}


def read_rows(ws, cols):
    rows = []
    for r in range(2, ws.max_row + 1):
        if ws.cell(r, cols['Ec']).value is None:
            continue
        rows.append(dict(row=r,
                         parity=str(ws.cell(r, cols['parity']).value).strip(),
                         J=float(ws.cell(r, cols['J']).value),
                         Ec=float(ws.cell(r, cols['Ec']).value),
                         energy=ws.cell(r, cols['energy']).value,
                         level_id=norm_id(ws.cell(r, cols['level_id']).value)))
    return rows


def insert_missing(ws, cols, rows, lids, calc, lid_of_id, id_of_lid, wanted):
    """Insert a row for each level of `wanted` that no row holds or can hold.

    The workbook stops at a calculated energy (about 141400 cm^-1 on
    2026-10-15), and a level found above it has no row.  Its row goes where
    its lid puts it within its parity, with what the calculation knows of
    it.  Returns [(level_id, lid, row number)] in the order inserted.
    """
    have = {r['level_id'] for r in rows if r['level_id']}
    have |= {id_of_lid[lid] for lid in lids if lid in id_of_lid}
    by_lid = {c['lid']: c for c in calc}
    todo = sorted((lid_of_id[i], i) for i in wanted
                  if i not in have and i in lid_of_id)
    out = []
    for lid, level_id in sorted(todo, reverse=True):
        c = by_lid[lid]
        par = c['parity'][0]
        same = [(l, r['row']) for r, l in zip(rows, lids)
                if l is not None and r['parity'] == par]
        before = [row for l, row in same if l < lid]
        at = (max(before) + 1 if before
              else min((row for _, row in same), default=ws.max_row + 1))
        ws.insert_rows(at)
        values = {'parity': par, 'J': c['J'], 'Ec': c['E_calc'],
                  'level_id': level_id, 'conf': c['conf'],
                  'term': c['term'] + ('*' if par == 'o' else '')}
        for name, value in values.items():
            if name in cols:
                ws.cell(at, cols[name]).value = value
        # the rows below have moved down by one
        for r in rows:
            if r['row'] >= at:
                r['row'] += 1
        out.append((level_id, lid, at))
    return out


def build(pub_path, set_dir, cfg, calc=None, lid_of_id=None):
    """Read everything and work the new cells out.

    Returns `(wb, ws, cells, report)`: the workbook, its sheet, the cells to
    write as {(row, column name): value}, and the report's lines.  Raises
    SystemExit when the table cannot be completed.
    """
    wb = openpyxl.load_workbook(pub_path)
    ws = wb.active
    cols = header_map(ws)
    missing = [c for c in ('parity', 'J', 'Ec', 'energy', 'level_id')
               if c not in cols]
    if missing:
        raise SystemExit('%s has no column %s' % (pub_path,
                                                  ', '.join(missing)))
    if calc is None:
        calc, lid_of_id = calculated_levels()
    id_of_lid = {v: k for k, v in lid_of_id.items()}
    values = level_values(read_tsv(os.path.join(set_dir, LOPT_LEVELS)),
                          read_csv(os.path.join(set_dir, SET3_LEVELS)))
    discarded = id_column(cfg.discarded_levels)

    rows = read_rows(ws, cols)
    lids, _ = match_rows(rows, calc, lid_of_id)
    inserted = insert_missing(ws, cols, rows, lids, calc, lid_of_id,
                              id_of_lid, set(values) | discarded)
    if inserted:
        rows = read_rows(ws, cols)
    lids, problems = match_rows(rows, calc, lid_of_id)
    notes = notes_of(id_column(cfg.new_levels, sep='\t'),
                     id_column(cfg.level_overrides), discarded)

    cells, filled, placed = {}, [], set()
    for r, lid in zip(rows, lids):
        level_id = r['level_id']
        if not level_id and lid is not None and lid in id_of_lid:
            level_id = id_of_lid[lid]
            cells[(r['row'], 'level_id')] = level_id
            filled.append(level_id)
        cells[(r['row'], 'lid')] = lid
        v = values.get(level_id)
        if v is not None:
            for name, x in zip(NEW_COLUMNS[1:6], v[:5]):
                cells[(r['row'], name)] = x
            placed.add(level_id)
        else:
            for name in NEW_COLUMNS[1:6]:
                cells[(r['row'], name)] = None
            if level_id and level_id not in discarded:
                problems.append('row %d: %s has no energy in the fit and is '
                                'not discarded' % (r['row'], level_id))
        cells[(r['row'], 'Notes')] = notes.get(level_id) or None
        r['level_id'] = level_id
    unplaced = sorted(set(values) - placed)
    if unplaced:
        raise SystemExit('%d levels of the fit are in no row of %s: %s'
                         % (len(unplaced), pub_path, ', '.join(unplaced)))
    unknown_notes = sorted(set(notes) - {r['level_id'] for r in rows})
    if unknown_notes:
        problems.append('levels named in the note files but in no row: %s'
                        % ', '.join(unknown_notes))

    report = ['%s: %d rows, %d with an energy of the fit; %d level_id cells '
              'filled' % (os.path.basename(pub_path), len(rows), len(placed),
                          len(filled))]
    for level_id, lid, row in inserted:
        report.append('  row %d inserted for %s (lid %d, beyond the table): '
                      'parity, conf, term, J and Ec from the calculation; '
                      'composition and Lande g left empty'
                      % (row, level_id, lid))
    if filled:
        report.append('  filled: ' + ' '.join(sorted(filled)))
    count = {k: sum(1 for r in rows if k in (notes.get(r['level_id']) or ''))
             for k in 'NRD'}
    report.append('Notes: %d N, %d R, %d D' % (count['N'], count['R'],
                                               count['D']))
    report += ['WARNING: ' + p for p in problems]
    report += shift_report(rows, values, notes)
    return wb, ws, cells, report


def shift_report(rows, values, notes):
    """Report lines: Eopt - energy for the revised levels, and the largest
    for the levels not marked R or N, after the change of frame.

    Eopt is a center of gravity measured from the ground level's center of
    gravity; a published energy is near the hfs head, measured from the
    ground level's head.  The change of frame alone moves a level by about
    -(S - S_ground), S the head's distance from the center of gravity
    (levels_set3.csv), so the list is ranked by what is left after it."""
    out, plain = [], []
    rev = []
    s0 = next((v[5] for i, v in values.items()
               if v[1] == 0.0 and v[0] == 0.0), 0.0)
    for r in rows:
        v = values.get(r['level_id'])
        if v is None or not isinstance(r['energy'], (int, float)):
            continue
        d = v[0] - float(r['energy'])
        mark = notes.get(r['level_id']) or ''
        if 'R' in mark:
            rev.append((abs(d), r['level_id'], d, v[3], v[5]))
        elif 'N' not in mark:
            left = d + (v[5] - s0)
            plain.append((abs(left), r['level_id'], d, v[3], left))
    rev.sort(reverse=True)
    plain.sort(reverse=True)
    if rev:
        out.append('revised levels (R), Eopt - energy (unc_tot), cm^-1:')
        out += ['  %s %+9.3f (%.3f)' % (i, d, u) for _, i, d, u, _ in rev]
    if plain:
        out.append('largest Eopt - energy of the levels not marked R or N, '
                   'for information (R comes from the revision list only), '
                   'ranked by what is left after the change of frame from '
                   'the hfs head to the center of gravity, -(S - S_ground), '
                   'S_ground = %.3f:' % s0)
        out.append('  %-13s %9s %7s %9s' % ('level_id', 'Eopt-E', 'unc_tot',
                                          'left'))
        out += ['  %s %+9.3f %7.3f %+9.3f' % (i, d, u, left)
                for _, i, d, u, left in plain[:N_SHIFTS]]
    return out


def write(wb, ws, cells, path):
    cols = header_map(ws)
    nxt = ws.max_column + 1
    for name in NEW_COLUMNS:
        if name not in cols:
            ws.cell(1, nxt, name)
            cols[name] = nxt
            nxt += 1
    # a rerun replaces the seven columns whole
    for name in NEW_COLUMNS:
        for r in range(2, ws.max_row + 1):
            ws.cell(r, cols[name]).value = None
    for (r, name), value in cells.items():
        ws.cell(r, cols[name]).value = value
    wb.save(path)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--pub', default=PUB, help='the workbook (default: '
                   '%(default)s)')
    p.add_argument('--set', dest='set_dir', default=SET_DIR,
                   help='the set whose fit is published (default: '
                        '%(default)s)')
    p.add_argument('--out', default=None,
                   help='write here instead of updating --pub in place')
    p.add_argument('--check', action='store_true',
                   help='report only, write nothing')
    args = p.parse_args(argv)
    out = args.out or args.pub
    if not args.check:
        output_files.require_writable([out])
    cfg = config.load()
    wb, ws, cells, report = build(args.pub, args.set_dir, cfg)
    if not args.check:
        write(wb, ws, cells, out)
        report.append('written: %s' % out)
    print('\n'.join(report))
    return 0


if __name__ == '__main__':
    sys.exit(main())
