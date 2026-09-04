#!/usr/bin/env python
"""Stage 5: what the missing-gA policy actually changed, assignment by assignment.

Compares two archived runs of ``classify_lines.py`` (typically
``baseline/policy_none`` and ``baseline/policy_impute``) and answers the three
questions of Stage 5 of PLAN_missing_gA.md:

1. Which accepted identifications were lost (accepted in the first run, rejected
   in the second) and which were gained, split by ``new`` (0 = legacy, taken from
   the input line list; 1 = proposed by this code) and by whether the observed
   line carried other accepted identifications - "blended" - or only this one -
   "sole".
2. Which decision channel rejected each lost identification.  The channel is read
   from the ``notes2`` text written by ``weed_assignments_line``:
     rel-intensity          = the relative-intensity filter (a candidate
                              contributing less than 10% - new - or 2% - legacy -
                              of the summed predicted intensity of the line);
     z-test                 = the comparison of predicted with observed intensity;
     insufficient-evidence  = a new candidate on a line with no predicted
                              intensity at all;
     undecided/other        = no rejecting note, i.e. the identification simply
                              never reached acceptance in the second run.
3. How the levels touched by those identifications moved: the level shift dE, the
   spurious-level probability p_spur, the number of identifications n_tot, and
   the question_status verdict.

Usage:  python tools/analyze_policy_diff.py OLD_DIR NEW_DIR [--list]
"""
import argparse
import os
import sys

import pandas as pd

KEY = ['wn_obs', 'low_id', 'upp_id']


def _indent(text, n):
    pad = ' ' * n
    return chr(10).join(pad + ln for ln in text.splitlines())


def load_lines(d):
    p = d if d.lower().endswith('.csv') else os.path.join(d, 'line_classifications.csv')
    df = pd.read_csv(p, dtype={'low_id': str, 'upp_id': str})
    df['wn_obs'] = df['wn_obs'].round(6)
    return df.dropna(subset=['low_id', 'upp_id']).set_index(KEY)


def load_levels(d):
    p = os.path.join(d, 'level_shift_report.csv')
    if not os.path.isfile(p):
        return None
    return pd.read_csv(p, dtype={'level_id': str}).set_index('level_id')


def channel(note):
    """Name the decision channel from the note left by the classifier."""
    t = '' if not isinstance(note, str) else note
    if 'of I_cum' in t:
        return 'rel-intensity'
    if 'incompatible with Iobs' in t:
        return 'z-test'
    if 'Insufficient evidence' in t:
        return 'insufficient-evidence'
    if t.strip() == '':
        return 'undecided'
    return 'other'


def accepted_per_line(df):
    """Number of accepted identifications of each observed line."""
    s = df['accepted'].fillna(0)
    return s.groupby(level='wn_obs').sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('old')
    ap.add_argument('new')
    ap.add_argument('--list', action='store_true',
                    help='print one line per lost/gained identification')
    a = ap.parse_args()

    A, B = load_lines(a.old), load_lines(a.new)
    nA, nB = accepted_per_line(A), accepted_per_line(B)
    common = A.index.intersection(B.index)
    accA = A.loc[common, 'accepted'].fillna(0)
    accB = B.loc[common, 'accepted'].fillna(0)

    lost = common[(accA == 1) & (accB != 1)]
    gained = common[(accA != 1) & (accB == 1)]
    print(f"accepted identifications: {int(A['accepted'].fillna(0).sum())} -> "
          f"{int(B['accepted'].fillna(0).sum())}   (on the {len(common)} rows "
          f"common to both runs: {int(accA.sum())} -> {int(accB.sum())})")
    print(f"lost: {len(lost)}    gained: {len(gained)}\n")

    for label, idx in (('LOST', lost), ('GAINED', gained)):
        if len(idx) == 0:
            continue
        rows = []
        for k in idx:
            r_old, r_new = A.loc[k], B.loc[k]
            # the note that explains the outcome: for a lost identification the
            # rejection written by the new run, for a gained one the acceptance
            # written by the new run as well
            deciding = r_new
            rows.append({
                'wn_obs': k[0], 'low_id': k[1], 'upp_id': k[2],
                'origin': 'new' if r_old['new'] == 1 else 'legacy',
                'imputed': int(r_new.get('imputed', 0) or 0),
                # context in the run where the identification WAS accepted
                'context': 'sole' if (nA if label == 'LOST' else nB)[k[0]] <= 1
                           else 'blended',
                'channel': channel(deciding['notes2']),
                'obs_intens': r_old['obs_intens'],
                'Ic_old': r_old['calc_intens'], 'Ic_new': r_new['calc_intens'],
                'note': deciding['notes2'],
            })
        t = pd.DataFrame(rows)
        print(f"=== {label} ({len(t)}) ===")
        print("  by origin:      ", dict(t['origin'].value_counts()))
        print(f"  imputed itself:  {int(t['imputed'].sum())} of {len(t)}")
        print("  line context where accepted:", dict(t['context'].value_counts()))
        print("  deciding channel x origin:")
        print(_indent(pd.crosstab(t['channel'], t['origin']).to_string(), 4))
        print("  deciding channel x imputed:")
        print(_indent(pd.crosstab(t['channel'], t['imputed']).to_string(), 4))
        print("  deciding channel x line context:")
        print(_indent(pd.crosstab(t['channel'], t['context']).to_string(), 4))
        if a.list:
            print("  detail:")
            for _, r in t.sort_values('wn_obs').iterrows():
                ic = r['Ic_new'] if pd.notna(r['Ic_new']) else float('nan')
                print(f"    {r['wn_obs']:12.4f} {r['low_id']}-{r['upp_id']} "
                      f"{r['origin']:6s} imp={r['imputed']} {r['context']:7s} "
                      f"Iobs={r['obs_intens']:.3g} Ic={ic:.3g} "
                      f"| {str(r['note'])[:95]}")
        print()

    LA, LB = load_levels(a.old), load_levels(a.new)
    if LA is None or LB is None:
        return
    lev = sorted(set(lost.get_level_values('low_id')) |
                 set(lost.get_level_values('upp_id')))
    lev = [x for x in lev if x in LA.index and x in LB.index]
    print(f"=== levels touched by the lost identifications: {len(lev)} ===")
    cmp = pd.DataFrame({
        'is_new_level': LA.loc[lev, 'is_new_level'],
        'dE_old': LA.loc[lev, 'dE'], 'dE_new': LB.loc[lev, 'dE'],
        'p_old': LA.loc[lev, 'p_spur'], 'p_new': LB.loc[lev, 'p_spur'],
        'n_old': LA.loc[lev, 'n_tot'], 'n_new': LB.loc[lev, 'n_tot'],
        'q_old': LA.loc[lev, 'question_status'].fillna('-'),
        'q_new': LB.loc[lev, 'question_status'].fillna('-'),
    })
    cmp['|dE|_change'] = cmp['dE_new'].abs() - cmp['dE_old'].abs()
    for tag, sub in (('all touched levels', cmp),
                     ('new levels only', cmp[cmp['is_new_level'] == 1])):
        if len(sub) == 0:
            continue
        print(f"  {tag} ({len(sub)}):")
        print(f"    rms |dE|:    {(sub['dE_old'] ** 2).mean() ** .5:.4f} -> "
              f"{(sub['dE_new'] ** 2).mean() ** .5:.4f} cm^-1")
        print(f"    mean p_spur: {sub['p_old'].mean():.4f} -> {sub['p_new'].mean():.4f}")
        print(f"    |dE| smaller on {(sub['|dE|_change'] < -1e-9).sum()} levels, "
              f"larger on {(sub['|dE|_change'] > 1e-9).sum()}, "
              f"unchanged on {(sub['|dE|_change'].abs() <= 1e-9).sum()}")
        ch = sub[sub['q_old'] != sub['q_new']]
        if len(ch):
            print(f"    question_status changed on {len(ch)}:")
            print(_indent(pd.crosstab(ch['q_old'], ch['q_new']).to_string(), 6))

    print("\n=== the whole level list ===")
    both = LA.index.intersection(LB.index)
    for tag, sub in (('all levels', both),
                     ('new levels', [x for x in both if LA.loc[x, 'is_new_level'] == 1])):
        d_old, d_new = LA.loc[sub, 'dE'], LB.loc[sub, 'dE']
        print(f"  {tag} ({len(sub)}): rms |dE| {(d_old ** 2).mean() ** .5:.4f} -> "
              f"{(d_new ** 2).mean() ** .5:.4f} cm^-1;  mean p_spur "
              f"{LA.loc[sub, 'p_spur'].mean():.4f} -> {LB.loc[sub, 'p_spur'].mean():.4f}")
    ch = pd.crosstab(LA.loc[both, 'question_status'].fillna('-'),
                     LB.loc[both, 'question_status'].fillna('-'))
    print("  question_status (old down, new across):")
    print(_indent(ch.to_string(), 4))


if __name__ == '__main__':
    sys.exit(main())
