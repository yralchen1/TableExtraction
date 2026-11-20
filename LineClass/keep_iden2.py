"""classify_lines.py --keep-IDEN2-ids: hold the classification to IDEN2.

Once the set of identifications is stable, IDEN2 is where it is kept: the
analyst judges every assignment there.  This mode makes a classification run
reproduce exactly what IDEN2/trans.dat of the working set identifies, instead
of the analyst copying each disagreement into line_decisions.csv by hand:

  * every transition IDEN2 identifies with a line is accepted on that line,
    whether or not the matching proposes it (it is then forced onto the line,
    as an "accept" row of the ledger would force it) and whatever the
    weeding says of it;
  * every other candidate is rejected.

Each enforced verdict that disagrees with what the run would otherwise have
done is written into the decision ledger, so that later runs, with or
without the option, keep it:

  * an IDEN2 assignment the matching did not propose, or the weeding
    rejected: an "accept" row, reason "kept IDEN2 assignment";
  * an assignment the weeding accepted that IDEN2 does not show: a "reject"
    row, reason "kept IDEN2 rejection";
  * a ledger row that contradicts IDEN2 is reversed in place, its old
    verdict and reason kept in the new reason.

An assignment the run accepts anyway, or rejects anyway, needs no row.  An
accepted row whose Ritz wavenumber misses its line by more than
decisions.max_forced_offset would stop every later run
(classify_lines.check_forced_decisions); since IDEN2 holds it deliberately,
its reason carries "offset-ok" as well.

Every enforced decision is listed in `LOG_NAME` in the set's directory,
rewritten by each run, and counted in a warning on the console.

The join to IDEN2 follows the project's rules: levels through
IDEN2/IDEN_level_ids.txt, never by energy; lines through IDEN2/dlv_keys.txt
(the wn_key each row of dlv.dat was written for), or, for a row edited in
IDEN2 since, by its wavenumber within classify_lines.DECISIONS_WN_MATCH.
Anything that does not join stops the run: a broken join is reported, not
guessed around.
"""
import bisect
import csv
import datetime
import io
import os

ACCEPT_REASON = 'kept IDEN2 assignment'
REJECT_REASON = 'kept IDEN2 rejection'
LOG_NAME = 'keep_IDEN2_log.txt'

# the categories of the log, in the order they are listed
MISSING = 'IDEN2 assignment not proposed by the matching: forced, accepted'
REJECTED = 'IDEN2 assignment rejected by the weeding: accepted'
EXTRA = 'accepted by the weeding, not in IDEN2: rejected'
REVERSED = 'ledger row contradicting IDEN2: reversed'
OFFSET = 'IDEN2 assignment beyond decisions.max_forced_offset: offset-ok'
CATEGORIES = (MISSING, REJECTED, EXTRA, REVERSED, OFFSET)


class _Collect(object):
    """A check_sync.Report stand-in that keeps the errors and warnings."""

    def __init__(self):
        self.errors, self.warnings = [], []

    def error(self, section, msg, *a, **k):
        self.errors.append(msg)

    def warn(self, section, msg, *a, **k):
        self.warnings.append(msg)

    def ok(self, *a, **k):
        pass

    def act(self, *a, **k):
        pass


def read_assignments(iden2_dir):
    """[(dlv line number, wn_key or None, wavenumber, level_id, level_id)]
    for every transition trans.dat identifies with an observed line.

    Raises ValueError if the IDEN2 files are missing or unreadable, or if an
    identified transition joins a level IDEN_level_ids.txt does not list.
    """
    import check_sync
    import swap_line_assignments_IDEN as IDEN
    rep = _Collect()
    iden = check_sync.load_iden2(iden2_dir, rep)
    if iden is None or rep.errors:
        raise ValueError('--keep-IDEN2-ids: the IDEN2 files of %s cannot be '
                         'read: %s' % (iden2_dir,
                                       '; '.join(rep.errors + rep.warnings)))
    out, unmapped = [], []
    for (owner, partner), k in iden.trans.row_of.items():
        obs = IDEN.assignment(iden.trans.records[k])
        if not IDEN.has_line(obs):
            continue
        row = IDEN.obs_row(obs)
        wn = iden.wn_of_line.get(row)
        if wn is None:
            wn = IDEN.obs_wavenumber(obs)
        a, b = iden.id_of_row.get(owner), iden.id_of_row.get(partner)
        if a is None or b is None:
            unmapped.append('rows %d - %d at %.3f' % (owner, partner, wn))
            continue
        out.append((row, iden.key_of_line.get(row), wn, a, b))
    if unmapped:
        raise ValueError('--keep-IDEN2-ids: %d identified transitions of '
                         'trans.dat join levels that IDEN_level_ids.txt does '
                         'not list: %s' % (len(unmapped), ', '.join(unmapped)))
    return out


def join(assignments, observed_lines, levels_dict, wn_match):
    """{id(line): (line, {(lower_id, upper_id), ...})} for the assignments.

    A line is found by the wn_key its dlv.dat row was written for, else by
    the row's wavenumber against the run's own, within `wn_match`.  The order
    of the two levels is the order of their energies in this run.  Raises
    ValueError, listing every failure, if a line or a level is not found.
    """
    by_key = {round(l.wn_key, 6): l for l in observed_lines}
    srt = sorted(observed_lines, key=lambda l: l.wavenumber)
    wns = [l.wavenumber for l in srt]

    def nearest(wn):
        j = bisect.bisect_left(wns, wn)
        best, best_d = None, wn_match
        for k in (j - 1, j):
            if 0 <= k < len(wns) and abs(wns[k] - wn) <= best_d:
                best, best_d = srt[k], abs(wns[k] - wn)
        return best

    out, faults = {}, []
    for row, key, wn, a, b in assignments:
        line = by_key.get(round(key, 6)) if key is not None else None
        if line is None:
            line = nearest(wn)
        if line is None:
            faults.append('dlv.dat line %d (%.4f): no observed line%s'
                          % (row, wn, '' if key is None
                             else ' with wn_key %.4f' % key))
            continue
        la, lb = levels_dict.get(a), levels_dict.get(b)
        missing = [x for x, lev in ((a, la), (b, lb)) if lev is None]
        if missing:
            faults.append('%.4f  %s - %s: level %s not in the level list'
                          % (line.wn_key, a, b, ', '.join(missing)))
            continue
        pair = (a, b) if la.energy <= lb.energy else (b, a)
        out.setdefault(id(line), (line, set()))[1].add(pair)
    if faults:
        raise ValueError('--keep-IDEN2-ids: %d identifications of IDEN2 do '
                         'not join the run:\n    %s'
                         % (len(faults), '\n    '.join(faults)))
    return out


class KeepState(object):
    """What the mode did before the run, for finish() to complete.

    `provisional` - {(id(line), pair)}: accepted verdicts laid on the lines
    for IDEN2's assignments that the ledger does not already accept;
    `reversed` - [(ledger key, line, new verdict, old verdict, old reason)];
    `offset` - [(line, pair, offset)] of IDEN2's assignments beyond the
    forced-offset limit.
    """

    def __init__(self):
        self.provisional = set()
        self.reversed = []
        self.offset = []
        self.offset_keys = set()
        self.n_iden2 = 0


def lay_on(observed_lines, levels_dict, sites, iden2_dir, accept, reject,
           wn_match):
    """Read IDEN2, give every line its `iden2_pairs`, and lay IDEN2's
    verdicts over the ledger's.

    `sites` is attach_line_decisions' list of (key, line, verdict); it is
    edited in place, so that check_forced_decisions and
    report_unapplied_decisions see the verdicts actually obeyed.  Returns a
    KeepState.
    """
    joined = join(read_assignments(iden2_dir), observed_lines, levels_dict,
                  wn_match)
    st = KeepState()
    for line in observed_lines:
        got = joined.get(id(line))
        line.iden2_pairs = set(got[1]) if got else set()
        st.n_iden2 += len(line.iden2_pairs)
    for i, (key, line, verdict) in enumerate(sites):
        pair = key[1:]
        decision, reason = verdict
        kept = pair in line.iden2_pairs
        if decision == accept and not kept:
            new = (reject, REJECT_REASON)
        elif decision == reject and kept:
            new = (accept, ACCEPT_REASON)
        else:
            continue
        line.decisions[pair] = new
        sites[i] = (key, line, new)
        st.reversed.append((key, line, new[0], decision, reason))
    for line in observed_lines:
        for pair in sorted(line.iden2_pairs):
            if pair not in line.decisions:
                line.decisions[pair] = (accept, ACCEPT_REASON)
                sites.append(((line.wn_key,) + pair, line,
                              (accept, ACCEPT_REASON)))
                st.provisional.add((id(line), pair))
    return st


def lenient_keys(st, sites, accept):
    """The ledger keys whose forced offset is to be warned of, not stopped
    on: every accepted verdict that IDEN2 holds."""
    return {key for key, line, (decision, _r) in sites
            if decision == accept and key[1:] in (line.iden2_pairs or ())}


def finish(observed_lines, st, accept, reject, offset_ok, unassigned):
    """Check that the run obeyed IDEN2 and collect the ledger rows it needs.

    Returns (new_rows, log_rows): new_rows as (line, low, upp, decision,
    reason) for the rows to append; log_rows as (category, line, low, upp,
    note).  Raises ValueError if an IDEN2 assignment is not accepted at the
    end of the run, or anything else is - which would mean the enforcement
    failed, and nothing must be written.
    """
    faults, new_rows, log_rows = [], [], []
    offset_pairs = {(id(l), p) for l, p, _o in st.offset}
    for line in observed_lines:
        if line.iden2_pairs is None:
            continue
        on = {}
        for t in line.assigned_transitions:
            if t is unassigned or not (t.lower_level and t.upper_level):
                continue
            on[(t.lower_level.level_id, t.upper_level.level_id)] = t
        for pair in sorted(line.iden2_pairs):
            t = on.get(pair)
            if t is None or t.accepted != 1:
                faults.append('%.4f  %s - %s: identified in IDEN2, %s at the '
                              'end of the run'
                              % (line.wn_key, pair[0], pair[1],
                                 'not a candidate' if t is None
                                 else 'not accepted'))
                continue
            if (id(line), pair) not in st.provisional:
                continue
            if t.forced:
                cat = MISSING
            elif t.algo_accepted != 1:
                cat = REJECTED
            else:
                continue
            reason = ACCEPT_REASON
            if (id(line), pair) in offset_pairs:
                reason += ' (%s)' % offset_ok
            new_rows.append((line, pair[0], pair[1], accept, reason))
            log_rows.append((cat, line, pair[0], pair[1], t.algo_note))
        for pair, t in sorted(on.items()):
            if pair in line.iden2_pairs:
                continue
            if t.accepted == 1:
                faults.append('%.4f  %s - %s: accepted at the end of the run '
                              'and not identified in IDEN2'
                              % (line.wn_key, pair[0], pair[1]))
            elif t.algo_accepted == 1 and pair not in line.decisions:
                new_rows.append((line, pair[0], pair[1], reject,
                                 REJECT_REASON))
                log_rows.append((EXTRA, line, pair[0], pair[1], t.algo_note))
    for key, line, new, old, reason in st.reversed:
        log_rows.append((REVERSED, line, key[1], key[2],
                         'was %s: %s' % (old, reason)))
    for line, pair, off in st.offset:
        log_rows.append((OFFSET, line, pair[0], pair[1],
                         'Ritz - observed %+.4f cm^-1' % off))
    if faults:
        raise ValueError('--keep-IDEN2-ids could not hold the run to IDEN2 '
                         '(nothing written):\n    ' + '\n    '.join(faults))
    return new_rows, log_rows


def today():
    """The ledger's date format, m/d/yyyy as Excel writes it."""
    d = datetime.date.today()
    return '%d/%d/%d' % (d.month, d.day, d.year)


def update_ledger(path, new_rows, st, accept, offset_ok):
    """Append `new_rows` to the ledger at `path` and reverse the rows of
    `st.reversed` in place, keeping the file's line ending and columns; a
    reversed row that now accepts beyond the forced-offset limit carries
    `offset_ok` too.  Returns (rows appended, rows reversed)."""
    if not (new_rows or st.reversed):
        return 0, 0
    with open(path, 'rb') as fh:
        raw = fh.read()
    bom = raw.startswith(b'\xef\xbb\xbf')
    text = raw.decode('utf-8-sig')
    nl = '\r\n' if '\r\n' in text else '\n'
    rows = list(csv.reader(io.StringIO(text, newline='')))
    head = rows[0]
    col = {c: i for i, c in enumerate(head)}
    date = today()

    def blank():
        return [''] * len(head)

    def put(r, name, value):
        if name in col:
            while len(r) <= col[name]:
                r.append('')
            r[col[name]] = value

    flip = {}
    for key, line, new, old, reason in st.reversed:
        flip[key] = (new, old, reason)
    n_rev = 0
    for r in rows[1:]:
        try:
            key = (float(r[col['wn_key']]), r[col['low_id']].strip(),
                   r[col['upp_id']].strip())
        except (ValueError, IndexError):
            continue
        hit = flip.get(key)
        if hit is None:
            continue
        new, old, reason = hit
        put(r, 'decision', new)
        put(r, 'date', date)
        put(r, 'reason', '%s (%swas %s: %s)' % (
            ACCEPT_REASON if new == accept else REJECT_REASON,
            offset_ok + '; ' if key in st.offset_keys else '', old, reason))
        n_rev += 1
    for line, low, upp, decision, reason in new_rows:
        r = blank()
        put(r, 'wn_key', '%.4f' % line.wn_key)
        put(r, 'low_id', low)
        put(r, 'upp_id', upp)
        put(r, 'decision', decision)
        put(r, 'date', date)
        put(r, 'reason', reason)
        rows.append(r)
    buf = io.StringIO(newline='')
    w = csv.writer(buf, lineterminator=nl)
    w.writerows(rows)
    data = buf.getvalue().encode('utf-8')
    with open(path, 'wb') as fh:
        fh.write((b'\xef\xbb\xbf' if bom else b'') + data)
    return len(new_rows), n_rev


def write_log(path, log_rows, st, ledger, n_new, n_rev):
    """The record of every enforced decision of the run."""
    count = {c: 0 for c in CATEGORIES}
    for row in log_rows:
        count[row[0]] += 1
    out = ['classify_lines.py --keep-IDEN2-ids, %s'
           % datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
           '%d identifications read from IDEN2' % st.n_iden2,
           '%s: %d rows appended, %d reversed' % (ledger, n_new, n_rev), '']
    for c in CATEGORIES:
        out.append('%4d  %s' % (count[c], c))
    for c in CATEGORIES:
        rows = [r for r in log_rows if r[0] == c]
        if not rows:
            continue
        out += ['', c, '-' * len(c),
                '%12s %12s  %-28s %s' % ('wn_key', 'wn_obs', 'lower - upper',
                                         'note')]
        for _c, line, low, upp, note in sorted(
                rows, key=lambda r: (r[1].wn_key, r[2], r[3])):
            out.append('%12.4f %12.4f  %s - %s %s'
                       % (line.wn_key, line.wavenumber, low, upp, note or ''))
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(out) + '\n')
    return count
