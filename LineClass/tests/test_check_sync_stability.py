"""Tests for the ledger half of `check_sync.check_stability`.

Run from the LineClass directory:  python -m pytest tests -q

`line_decisions.csv` (the ledger) holds the analyst's verdicts, one row per
transition on one observed line, named by Sugar's wavenumber (`wn_key`).  An
assignment is moved from one line to another by two rows: a `reject` on the
old line and an `accept` on the new one.  The check must read such a pair of
rows as one ruling obeyed, and warn only when a rejected transition has
settled on a line that no ledger row names.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import check_sync                          # noqa: E402

LOW, UPP = '000260', '000642'


def classification(rows):
    """line_classifications.csv rows: (wn_key, low_id, upp_id, accepted)."""
    return pd.DataFrame([dict(wn_key=w, wn_obs=w, low_id=lo, upp_id=up,
                              accepted=acc) for w, lo, up, acc in rows])


def run(tmp_path, cls_rows, dec_rows):
    path = tmp_path / check_sync.DECISIONS
    pd.DataFrame([dict(wn_key=w, low_id=lo, upp_id=up, decision=d)
                  for w, lo, up, d in dec_rows]).to_csv(path, index=False)
    cls = classification(cls_rows)
    rep = check_sync.Report()
    check_sync.check_stability({check_sync.DECISIONS: str(path)}, cls,
                               check_sync.Context(cls=cls), rep, 10)
    return [f for f in rep.findings if f['check'] == 'Stability'
            and 'ledger' in f['message']]


def test_moved_assignment_is_not_a_warning(tmp_path):
    """The 000642 case: rejected on the F=2 line, accepted on F=3."""
    found = run(tmp_path,
                [(21847.197, LOW, UPP, 1)],
                [(21848.700, LOW, UPP, 'reject'),
                 (21847.197, LOW, UPP, 'accept')])
    assert [f['severity'] for f in found] == [check_sync.OK]


def test_unruled_new_line_is_still_a_warning(tmp_path):
    """Rejected on one line, which is no longer a candidate for the pair, and
    settled on another with no ledger row there."""
    found = run(tmp_path,
                [(21847.197, LOW, UPP, 1)],
                [(21848.700, LOW, UPP, 'reject')])
    assert [f['severity'] for f in found] == [check_sync.WARN]
    assert len(found[0]['items']) == 1


def test_accept_on_a_third_line_does_not_excuse(tmp_path):
    """An accept row for the pair on some other line is not a ruling on the
    line the transition settled on."""
    found = run(tmp_path,
                [(21847.197, LOW, UPP, 1)],
                [(21848.700, LOW, UPP, 'reject'),
                 (21900.000, LOW, UPP, 'accept')])
    assert check_sync.WARN in [f['severity'] for f in found]
