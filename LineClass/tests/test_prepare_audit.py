"""Tests for `prepare_audit`: a new audit made into the review copy, beside
the user's verdicts (level_verdicts.txt).

Run from the LineClass directory:  python -m pytest tests -q
"""
import os
import sys

import openpyxl
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import prepare_audit as P                  # noqa: E402

A, B, C, D, E = ('059003.000551', '059003.000678', '059003.000587',
                 '059003.000484', '059003.000001')
AUDIT_HEAD = ['level_id', 'E', 'dE_alt', 'question', 'q_why', 'action']
IDS = [(A, 849), (B, 690), (C, 414), (D, 700), (E, 1253)]


def workbook(path, head, rows):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(head)
    for r in rows:
        ws.append(r)
    wb.save(path)
    return str(path)


def id_map(d, pairs=IDS):
    path = d / 'IDEN_level_ids.txt'
    path.write_text('level_id\tIDEN_id\n' + ''.join(
        '%s\t%d\n' % p for p in pairs), encoding='utf-8')
    return str(path)


def verdicts(d, rows):
    path = d / 'level_verdicts.txt'
    path.write_text('level_id\tverdict\tE_alt\tdate\tnote\n' + ''.join(
        '\t'.join(r) + '\n' for r in rows), encoding='utf-8')
    return str(path)


def read(path):
    ws = openpyxl.load_workbook(path).active
    it = ws.iter_rows(values_only=True)
    head = list(next(it))
    return head, {r[0]: dict(zip(head, r)) for r in it}, ws


def audit(d):
    return workbook(d / 'level_positions.xlsx', AUDIT_HEAD, [
        [A, 106384.4862, -18.7666, '?', 'alternate', 'weak'],
        [B, 118670.0995, -66.0999, '?', 'alternate', None],
        [C, 134480.0, 47.1124, '?', 'alternate', None],
        [D, 116000.0, None, '?', 'S -4.05 < 3', None],
        [E, 0.0, None, None, None, None]])


def test_the_columns_come_from_the_verdicts(tmp_path):
    v = verdicts(tmp_path, [
        (A, '?', '106365.7196', '2026-10-16', ''),
        (B, 'firm', '118605.0', '2026-10-16', 'JK pair'),   # 1 cm^-1 off
        (C, 'firm', '134000.0', '2026-10-16', ''),
        (C, 'firm', '134527.0', '2026-11-19', 'second look')])
    out = str(tmp_path / 'review.xlsx')
    text = P.prepare(audit(tmp_path), id_map(tmp_path), v, out)
    head, rows, ws = read(out)
    assert head == AUDIT_HEAD + list(P.NEW_COLUMNS)
    a, b, c, d, e = (rows[k] for k in (A, B, C, D, E))
    assert a['IDEN_id'] == 849 and e['IDEN_id'] == 1253
    assert a['E_new'] == pytest.approx(106365.7196)
    assert (a['verdict'], a['E_alt'], a['d'], a['open']) == \
        ('?', pytest.approx(106365.7196), 0, None)
    # firm, but the audit's alternate is not the one denied: open
    assert b['verdict'] == 'firm' and b['note'] == 'JK pair'
    assert b['d'] == pytest.approx(-1.0004) and b['open'] == 'new alternate'
    # of two alternates denied, the one nearest the suggestion is shown
    assert c['E_alt'] == pytest.approx(134527.0)
    assert c['d'] == pytest.approx(0.1124) and c['open'] is None
    assert c['note'] == 'second look'
    # questioned by the policy and without a verdict: open; not questioned:
    # nothing to do
    assert d['open'] == 'no verdict' and d['E_new'] is None
    assert e['open'] is None and e['verdict'] is None
    assert ws.freeze_panes == 'A2'
    assert 'questionable levels settled by a verdict: 2' in text
    assert 'questionable levels open (to look at in IDEN2): 2' in text
    assert ('  row 3  %s  new alternate at 118603.9996, the nearest denied '
            '118605.0000' % B) in text
    assert '  row 5  %s  no verdict' % D in text


def test_a_level_of_the_verdicts_missing_from_the_audit_is_reported(tmp_path):
    v = verdicts(tmp_path, [('059003.000999', 'firm', '', '', '')])
    text = P.prepare(audit(tmp_path), id_map(tmp_path), v,
                     str(tmp_path / 'review.xlsx'))
    assert any(t.startswith('WARNING') and '059003.000999' in t
               for t in text)


def test_without_verdicts_every_questionable_level_is_open(tmp_path):
    out = str(tmp_path / 'review.xlsx')
    P.prepare(audit(tmp_path), id_map(tmp_path), '', out)
    _, rows, _ = read(out)
    assert [rows[k]['open'] for k in (A, B, C, D, E)] == \
        ['no verdict'] * 4 + [None]


def test_a_level_missing_from_the_lookup_table_stops_the_run(tmp_path):
    out = tmp_path / 'review.xlsx'
    with pytest.raises(SystemExit, match='lookup table is broken'):
        P.prepare(audit(tmp_path), id_map(tmp_path, IDS[1:]), '', str(out))
    assert not out.exists()


def test_a_review_copy_is_not_taken_for_the_audit(tmp_path):
    path = workbook(tmp_path / 'level_positions.xlsx',
                    AUDIT_HEAD + ['IDEN_id'], [[A, 1.0, None, None, None,
                                                None, 1]])
    with pytest.raises(SystemExit, match='review copy'):
        P.prepare(path, id_map(tmp_path), '', str(tmp_path / 'out.xlsx'))


def test_the_program_writes_the_review_copy_and_reads_the_audit_only(
        tmp_path):
    set_dir = tmp_path / 'iter_hfs'
    (set_dir / 'IDEN2').mkdir(parents=True)
    id_map(set_dir / 'IDEN2')
    v = verdicts(tmp_path, [(A, 'firm', '106365.72', '', '')])
    path = audit(set_dir)
    stamp = os.path.getmtime(path)
    P.main(['--set', str(set_dir), '--verdicts', v])
    _, rows, _ = read(str(set_dir / (P.REVIEW_NAME % 'iter_hfs')))
    assert rows[A]['verdict'] == 'firm' and rows[A]['open'] is None
    assert os.path.getmtime(path) == stamp
    # the copy holds nothing of its own: a second run simply replaces it
    P.main(['--set', str(set_dir), '--verdicts', v])
