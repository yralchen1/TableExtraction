"""Prepare a new level audit for review, beside the verdicts already given.

Run from LineClass/ after `level_positions.py --audit`:

    python prepare_audit.py                    # the iter_hfs set
    python prepare_audit.py --set iter         # another set
    python prepare_audit.py --out FILE.xlsx    # write elsewhere (a trial)

WHAT IT DOES
============
The audit (`level_positions.py --audit`) writes `<set>/level_positions.xlsx`,
one row per level; a questionable level has `question` = '?', and an
alternate position, if the audit found one, `dE_alt` away from its energy E.
The user's decisions on questionable levels are kept in one place only, the
verdicts file (level_verdicts.txt, files.level_verdicts; level_verdicts.py):
`firm` or `?`, the alternate positions looked at and denied (E_alt), and a
note.  This program writes the review copy,
`<set>/level_positions_audit_<set>.xlsx`: level_positions.xlsx with these
columns added after its last one (action), as values:

  IDEN_id   the level's row number in IDEN2's enlev.dat, through
            `<set>/IDEN2/IDEN_level_ids.txt` (never by energy).  A level
            missing from that table stops the run: the table is broken.
  E_new     E + dE_alt, the alternate position the audit suggests now
  verdict   the level's verdict in the verdicts file
  E_alt     its denied alternate nearest E_new (the first, without E_new)
  d         E_new - E_alt
  note      the verdict's note
  open      why a questionable level needs looking at: `no verdict`, or
            `new alternate` - one no E_alt lies within --tol (0.5 cm^-1) of.
            Blank where a verdict settles it.

The open levels are the ones to look at in IDEN2; a decision on one goes
into the verdicts file, and a rerun of this program shows it.  The review
copy holds nothing of its own: it is replaced on every run (git keeps the
previous one), and level_positions.xlsx is only read, so it keeps its time
stamp.
"""
import argparse
import os
import sys
import tempfile

import openpyxl

import config
import cowan_gA
import level_verdicts
import output_files

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SET = 'iter_hfs'
AUDIT_NAME = 'level_positions.xlsx'
REVIEW_NAME = 'level_positions_audit_%s.xlsx'
ID_MAP = os.path.join('IDEN2', 'IDEN_level_ids.txt')
#: the columns this program adds, in order
NEW_COLUMNS = ('IDEN_id', 'E_new', 'verdict', 'E_alt', 'd', 'note', 'open')


def number(v):
    """`v` as a float, or None if the cell holds no number."""
    if v is None or isinstance(v, bool):
        return None
    if isinstance(v, (int, float)):
        return float(v)
    try:
        return float(str(v).strip())
    except ValueError:
        return None


def header_index(ws, path, names):
    """`{name: 1-based column}` of `names` in the first row of `ws`; stops if
    one is missing."""
    head = [c.value for c in ws[1]]
    out = {}
    for n in names:
        if n not in head:
            raise SystemExit('prepare_audit.py: %s has no column %r'
                             % (path, n))
        out[n] = head.index(n) + 1
    return out


def open_reason(questioned, e_new, rec, tol):
    """Why a level needs looking at, or '' if it does not: only a
    questioned level does, when no verdict is given, or when the audit
    suggests an alternate no E_alt of the verdict lies within `tol` of."""
    if not questioned:
        return ''
    if rec is None:
        return 'no verdict'
    if e_new is not None and not any(abs(e_new - e) < tol
                                     for e in rec.e_alt):
        return 'new alternate'
    return ''


def prepare(audit_path, id_map_path, verdicts_path, out_path,
            tol=level_verdicts.VERDICT_TOL):
    """Write the review copy `out_path` from `audit_path`; returns the
    summary's lines."""
    id_of = {lid: iden for iden, lid in
             cowan_gA.read_id_map(id_map_path).items()}
    verdicts = level_verdicts.read(verdicts_path)

    wb = openpyxl.load_workbook(audit_path)
    ws = wb.active
    col = header_index(ws, audit_path, ('level_id', 'E', 'dE_alt',
                                        'question'))
    head = [c.value for c in ws[1]]
    clash = [n for n in NEW_COLUMNS if n in head]
    if clash:
        raise SystemExit('prepare_audit.py: %s already has the column(s) %s;'
                         ' it is a review copy, not the audit'
                         % (audit_path, ', '.join(clash)))
    last = len(head)
    while last and head[last - 1] in (None, ''):
        last -= 1
    c = {n: last + 1 + k for k, n in enumerate(NEW_COLUMNS)}
    for n, k in c.items():
        ws.cell(row=1, column=k, value=n)

    missing, seen, opened = [], set(), []
    n_q = n_alt = n_q_alt = n_settled = 0
    for row in range(2, ws.max_row + 1):
        lid = ws.cell(row=row, column=col['level_id']).value
        if lid is None or str(lid).strip() == '':
            continue
        lid = str(lid).strip()
        seen.add(lid)
        questioned = bool(ws.cell(row=row, column=col['question']).value)
        iden = id_of.get(lid)
        if iden is None:
            missing.append(lid)
        ws.cell(row=row, column=c['IDEN_id'], value=iden)
        E = number(ws.cell(row=row, column=col['E']).value)
        dE = number(ws.cell(row=row, column=col['dE_alt']).value)
        e_new = round(E + dE, 4) if E is not None and dE is not None else None
        ws.cell(row=row, column=c['E_new'], value=e_new)
        n_q += questioned
        n_alt += e_new is not None
        n_q_alt += questioned and e_new is not None
        rec = verdicts.get(lid)
        e_alt = None
        if rec is not None:
            ws.cell(row=row, column=c['verdict'], value=rec.verdict)
            ws.cell(row=row, column=c['note'], value=rec.note or None)
            e_alt = level_verdicts.nearest(rec, e_new)
            if e_alt is not None:
                ws.cell(row=row, column=c['E_alt'], value=round(e_alt, 4))
                if e_new is not None:
                    ws.cell(row=row, column=c['d'],
                            value=round(e_new - e_alt, 4))
        why = open_reason(questioned, e_new, rec, tol)
        if why:
            ws.cell(row=row, column=c['open'], value=why)
            opened.append((row, lid, why, e_new, e_alt))
        elif questioned:
            n_settled += 1
    if missing:
        raise SystemExit(
            'prepare_audit.py: %d level(s) of %s are not in %s: %s\n'
            'Every level of the set is in enlev.dat, so the lookup table is '
            'broken; nothing was written.'
            % (len(missing), audit_path, id_map_path, ', '.join(missing[:10])))

    ws.freeze_panes = 'A2'
    ws.auto_filter.ref = 'A1:%s%d' % (
        openpyxl.utils.get_column_letter(c[NEW_COLUMNS[-1]]), ws.max_row)
    # written beside the target and moved over it, so that a failed save
    # leaves the previous copy whole
    fd, tmp = tempfile.mkstemp(suffix='.xlsx',
                               dir=os.path.dirname(os.path.abspath(out_path)))
    os.close(fd)
    try:
        wb.save(tmp)
        os.replace(tmp, out_path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)

    stray = sorted(set(verdicts) - seen)
    text = ['%s: %d levels, %d with an alternate position; %d questionable, '
            '%d of them with an alternate position and %d without (the '
            '[audit] policy marks, in q_why)'
            % (audit_path, ws.max_row - 1, n_alt, n_q, n_q_alt,
               n_q - n_q_alt),
            'verdicts: %s (%d levels: %d firm, %d ?)'
            % (verdicts_path or 'none', len(verdicts),
               sum(v.verdict == 'firm' for v in verdicts.values()),
               sum(v.verdict == '?' for v in verdicts.values())),
            'questionable levels settled by a verdict: %d' % n_settled,
            'questionable levels open (to look at in IDEN2): %d' % len(opened)]
    for row, lid, why, e_new, e_alt in opened:
        extra = ''
        if why == 'new alternate':
            extra = ' at %.4f' % e_new + (
                ', the nearest denied %.4f' % e_alt if e_alt is not None
                else ', none denied')
        text.append('  row %d  %s  %s%s' % (row, lid, why, extra))
    if stray:
        text.append('WARNING: %d level(s) of the verdicts file are not in the '
                    'audit: %s' % (len(stray), ', '.join(stray)))
    text.append('written: %s' % out_path)
    return text


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--set', default=DEFAULT_SET, dest='set_dir',
                   help='the set whose audit is prepared (default: '
                        '%(default)s)')
    p.add_argument('--out', default=None,
                   help='write here instead of over the review copy')
    p.add_argument('--verdicts', default=None,
                   help="the verdicts file (default: the set's "
                        'files.level_verdicts)')
    p.add_argument('--tol', type=float, default=level_verdicts.VERDICT_TOL,
                   help='|E_new - E_alt| within which an alternate is the one '
                        'denied, cm^-1 (default: %(default)s)')
    args = p.parse_args(argv)

    set_dir = (args.set_dir if os.path.isabs(args.set_dir)
               else os.path.join(HERE, args.set_dir))
    name = os.path.basename(os.path.normpath(set_dir))
    audit = os.path.join(set_dir, AUDIT_NAME)
    out = args.out or os.path.join(set_dir, REVIEW_NAME % name)
    if not os.path.isfile(audit):
        raise SystemExit('prepare_audit.py: no audit %s; run '
                         'level_positions.py --audit first' % audit)
    verdicts = args.verdicts
    if verdicts is None:
        verdicts = config.load(os.path.join(
            set_dir, config.CONFIG_NAME)).level_verdicts
    output_files.require_writable([out])
    print('\n'.join(prepare(audit, os.path.join(set_dir, ID_MAP), verdicts,
                            out, args.tol)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
