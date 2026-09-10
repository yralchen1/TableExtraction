#!/usr/bin/env python
"""Which levels nobody has found yet are the ones most likely to be found.

Run from inside LineClass/:

    python unfound_levels.py                  # the ranked table
    python unfound_levels.py --detail 742     # the transitions of one level
    python unfound_levels.py --min-prom 3     # a shorter list


1.  The question
----------------
Cowan's calculation gives Pr III 1253 levels.  594 of them have been found in
the spectrum; the other 659 have never been placed, and each of them is an
energy the calculation predicts and the observed line list has never been
searched for.  Searching for one is work - it means putting a trial energy into
IDEN2 and looking at what lines fall on the transitions the level would then
have - so the question this module answers is which of the 659 are worth the
work.

A level can only be found through its lines, so the answer is a count: how many
of the transitions the calculation gives it would have been recorded on the
plates, if the level is where the calculation puts it.  A level with one such
transition cannot be found at all - one line can be made to fit any energy, so
a single coincidence is no evidence - and a level with ten is a level whose
position, if it is there at all, is heavily overdetermined.

That count is ``n_prom`` below, and the table is sorted by it.


2.  What counts as a transition that would have been recorded
-------------------------------------------------------------
For an unfound level L at a trial energy E and a partner level M that HAS been
found, at its measured energy E_M, the transition L-M would sit at the
wavenumber

    nu = |E - E_M|                       cm^-1

and would have the predicted intensity, on the same scale as every other
predicted intensity of this pipeline,

    I = C gA (nu/1e8) exp(-E_up/kT)      E_up = max(E, E_M)

with gA (the statistical weight of the upper level times the transition
probability, s^-1) read from Cowan's own transition list through
``cowan_gA.py``, and C and kT the constants of ``[intensity_model]`` in
``lineclass_config.toml``, fitted to the accepted observed intensities.  This
is exactly the number ``Icalc.xlsx`` carries for the transitions between two
found levels - recomputing those from gA here reproduces the file to an rms of
0.2 % over all 29260 of them - so an unfound level's transitions arrive on the
intensity scale the rest of the pipeline already works in, which is the whole
point of computing them from gA rather than from anywhere else.  ``Icalc.xlsx``
itself cannot be used: it holds only transitions whose two levels are both
known, and by definition none of these are.

The probability that such a line would have appeared in Sugar's list at all is
then the same one every other module uses,

    P_obs = c(lambda) D(z)

- ``c`` the coverage function measured by ``tools/coverage_map.py`` (what
fraction of the plates were looking at that wavelength at all) and ``D`` the
detection curve measured by ``tools/obscuration_rate.py`` (how often a line of
that brightness relative to the local noise was in fact recorded), reached
through ``level_positions.p_obs_vec``.  A transition is **promising** when
P_obs is at least ``P_PROMISING`` (0.5): better than an even chance of being
on the plates.  ``n_obs`` counts the weaker bar ``P_SEEN`` (0.2) that
``level_positions.py`` uses, and is reported beside it for comparison.


3.  Where the level might be, and why that has to be averaged over
------------------------------------------------------------------
E is not known.  All that is known is that the level is somewhere near its
calculated energy, and how near is a property of the CONFIGURATION rather than
of the level, because the calculation describes some configurations far better
than others: the rms of E_obs - E_calc over the found levels of the same
configuration - ``level_interchange.configuration_windows``, the same W the
alternate-position scan uses - runs from 24 cm^-1 for 4f^2.5g to 404 cm^-1 for
4f.5d.6p.

So nu is uncertain by that much, and with it lambda, the local noise level and
therefore P_obs.  Rather than pretend E = E_calc, every probability here is
averaged over a Gaussian of width W centred on E_calc: P_obs is evaluated on
the grid E_calc + z W for the z of ``Z_GRID`` and combined with the weights of
the normal density.  What comes out is then a probability in the ordinary
sense - the chance that the line is on the plates, given both the physics and
our ignorance of where the level sits - and summing it over a level's
transitions gives ``sum_P``, the expected NUMBER of its lines that were
recorded.

Six configurations - 5d3, f26g, f5d6d, f5d7s, f6s2, p5f3d - have not one level
found in them, so they have no W of their own.  Their levels are still listed,
using the rms over the level list as a whole (132 cm^-1) and marked ``list`` in
the ``W_src`` column, because a ranking is a suggestion about where to spend an
afternoon and not a number written into a file.  Read those rows knowing that
the width used for them is a guess.


4.  What the count does not say
-------------------------------
``n_prom`` says the evidence would exist, not that it can be recognised.  Three
things it does not know:

- **Whether the lines are still free.**  A promising transition may fall on a
  feature some other level has already been given.  Where the level might be
  anywhere in a window hundreds of wavenumbers wide, there is no way to say in
  advance which lines it would reach; ``rho``, the density of recorded lines in
  the neighbourhood of its transitions, is the closest thing to a warning the
  table can carry.  A level whose promising transitions all lie in a crowded
  region will find something wherever it is put, and finding something is not
  the same as finding it.
- **Whether the partners are sound.**  Every partner here is a found level, but
  found levels differ in how well they are established; a transition to a level
  that is itself in doubt is weaker evidence than the count suggests.
  ``level_positions.py`` is what says which is which.
- **Whether the calculated strengths are right.**  gA is a calculated quantity
  and the observed intensities scatter against the predicted ones by a factor
  of about three (s = 1.2 in the log).  P_obs already carries that spread
  through the detection curve, but a transition predicted at ten times the
  noise can still be absent.

The tool that answers "and is it there?" is the scan of ``level_positions.py``
at a trial energy.  This one only says where to point it.
"""
import argparse
import os

import numpy as np
import pandas as pd

import chance_mc as mc
import classify_lines as cl
import cowan_gA
import level_interchange as li
import level_positions as lp
import output_files

HERE = os.path.dirname(os.path.abspath(__file__))
ENLEV = os.path.join(HERE, 'IDEN2', 'enlev.dat')
IDS = os.path.join(HERE, 'IDEN2', 'IDEN_level_ids.txt')
OUT_CSV = os.path.join(HERE, 'unfound_levels.csv')

P_PROMISING = 0.5    # P_obs at which a transition is worth searching on
MIN_PROM = 2         # promising transitions a level needs to be searchable at
                     # all: one line fits any energy, so one is no evidence
# The trial positions, in units of the configuration's own W, and their
# weights: a normal density truncated at two sigma and renormalised.  Nine
# points rather than three because the noise level and the coverage vary
# sharply with wavelength, and a level of 4f.5d.6p is being moved over
# +/- 800 cm^-1.
Z_GRID = np.linspace(-2.0, 2.0, 9)
Z_WEIGHT = np.exp(-0.5 * Z_GRID ** 2)
Z_WEIGHT = Z_WEIGHT / Z_WEIGHT.sum()


# ---------------------------------------------------------------------------
# The theoretical level list, and its transitions to the found levels
# ---------------------------------------------------------------------------
def read_theory(enlev=ENLEV, ids=IDS, log=print):
    """The calculated levels of IDEN2, joined to Cowan's transition list.

    Returns ``(en, trans, mapping)``: ``en`` the frame of
    ``level_interchange.read_enlev`` indexed by the IDEN2 row number, ``trans``
    Cowan's E1 transition list from ``cowan_gA.read_transitions``, and
    ``mapping`` the ``{Cowan level number: IDEN2 row}`` correspondence, which
    ``cowan_gA.match_to_enlev`` checks against every level whose identifier
    both files know before it returns.
    """
    en = li.read_enlev(enlev)
    trans = cowan_gA.read_transitions(log=None)
    calc = cowan_gA.levels(trans)
    mapping, rep = cowan_gA.match_to_enlev(
        calc, cowan_gA.read_enlev_levels(enlev), cowan_gA.read_id_map(ids))
    log(f"calculated levels: {len(en)}, of which {int(en['known'].sum())} "
        f"found and {int((~en['known']).sum())} not")
    log(f"  matched to Cowan's transition list: {rep['n_matched']} "
        f"({rep['n_checked']} checked against their level_id, largest energy "
        f"disagreement {rep['max_dE']:.1f} cm^-1)")
    log(f"calculated E1 transitions: {len(trans)}")
    return en.set_index('idx'), trans, mapping


def windows(en, log=print):
    """``({configuration: W}, W of the whole list)``.

    W is the rms of E_obs - E_calc over the FOUND levels of a configuration:
    how far from its calculated position a level of that configuration turns
    out to be.  A configuration with no found level has no W and takes the
    list-wide value, flagged in the table.
    """
    known = en[en['known']]
    per = {c: float(np.sqrt(float((g['omc'] ** 2).mean())))
           for c, g in known.groupby('cfg')}
    whole = float(np.sqrt(float((known['omc'] ** 2).mean())))
    missing = sorted(set(en['cfg']) - set(per))
    if missing:
        log(f"  no found level in {len(missing)} configurations "
            f"({', '.join(missing)}); their unfound levels are placed with the "
            f"list-wide W = {whole:.1f} cm^-1 and marked 'list'")
    return per, whole


def unfound_transitions(en, trans, mapping, e_meas, log=print):
    """One row per transition between an unfound level and a found one.

    Columns ``idx`` (the unfound level, as an IDEN2 row number), ``partner``
    (the found one), ``gA``, ``E_calc`` (of the unfound level), ``E_partner``
    (measured) and ``upper`` (``'L'`` when the unfound level is the upper of
    the two at its calculated energy, ``'M'`` when the partner is).

    Transitions between two unfound levels are dropped: neither end is known,
    so no wavenumber can be predicted for them and they cannot be searched for.
    Transitions between two found levels are not this module's business.
    """
    a = np.array([mapping.get(int(v), -1) for v in trans['lid1']])
    b = np.array([mapping.get(int(v), -1) for v in trans['lid2']])
    gA = trans['gA'].to_numpy(dtype=float)
    ok = (a > 0) & (b > 0) & (gA > 0)
    a, b, gA = a[ok], b[ok], gA[ok]

    known = en['known'].to_dict()
    ka = np.array([known[i] for i in a], dtype=bool)
    kb = np.array([known[i] for i in b], dtype=bool)
    one = ka ^ kb                      # exactly one end found
    a, b, gA = a[one], b[one], gA[one]
    unk = np.where(ka[one], b, a)
    prt = np.where(ka[one], a, b)

    e_calc = en['E_calc'].to_dict()
    E_u = np.array([e_calc[i] for i in unk])
    E_p = np.array([e_meas[i] for i in prt])
    log(f"  {len(unk)} of them join an unfound level to a found one, over "
        f"{len(set(unk.tolist()))} unfound levels")
    return pd.DataFrame(dict(idx=unk, partner=prt, gA=gA,
                             E_calc=E_u, E_partner=E_p,
                             upper=np.where(E_u > E_p, 'L', 'M')))


def found_transitions(en, trans, mapping, e_meas):
    """The same frame as :func:`unfound_transitions`, for the FOUND levels.

    Every found level in turn plays the part of the unfound one - it is put at
    its CALCULATED energy, as though nobody knew where it really was - and its
    transitions to the other found levels are collected.  This is what
    ``--validate`` scores the count on: for these levels the answer is known,
    because the run says how many accepted lines each of them actually has.
    """
    a = np.array([mapping.get(int(v), -1) for v in trans['lid1']])
    b = np.array([mapping.get(int(v), -1) for v in trans['lid2']])
    gA = trans['gA'].to_numpy(dtype=float)
    ok = (a > 0) & (b > 0) & (gA > 0)
    a, b, gA = a[ok], b[ok], gA[ok]
    known = en['known'].to_dict()
    both = (np.array([known[i] for i in a], dtype=bool)
            & np.array([known[i] for i in b], dtype=bool))
    a, b, gA = a[both], b[both], gA[both]
    idx = np.concatenate([a, b])
    prt = np.concatenate([b, a])
    e_calc = en['E_calc'].to_dict()
    E_u = np.array([e_calc[i] for i in idx])
    E_p = np.array([e_meas[i] for i in prt])
    return pd.DataFrame(dict(idx=idx, partner=prt,
                             gA=np.concatenate([gA, gA]),
                             E_calc=E_u, E_partner=E_p,
                             upper=np.where(E_u > E_p, 'L', 'M')))


def validate(ctx, en, trans, mapping, per_cfg, whole, e_meas,
             p_prom=P_PROMISING, ids=IDS, log=print):
    """Does the count predict anything?  Score it on the levels already found.

    A level's ``n_prom`` is computed from the calculation alone - calculated
    energy, calculated gA, the measured plate coverage and noise - and never
    looks at whether any line was in fact assigned to the level.  So it can be
    scored on the 594 levels that HAVE been found by computing it for them as
    though they had not, and comparing it with the number of accepted lines
    each of them really carries in the run.

    Returns the frame it scored.
    """
    t = found_transitions(en, trans, mapping, e_meas)
    t = attach_windows(t, en, per_cfg, whole)
    t = observation_probabilities(t, ctx)
    tab = level_table(t, en, p_prom)
    id_of = cowan_gA.read_id_map(ids)
    tab['level_id'] = [id_of.get(int(i), '') for i in tab['idx']]
    tab = tab[tab['level_id'] != ''].copy()
    tab['n_acc'] = [ctx.n_acc_level.get(l, 0) for l in tab['level_id']]

    log(f"\nthe count scored on the {len(tab)} levels that have been found, "
        f"each of them put back at its\ncalculated energy and counted as "
        f"though nobody knew where it was:")
    log(f"\n  {'n_prom':>8} {'levels':>7} {'median lines':>13} "
        f"{'mean lines':>11}")
    for lo, hi, name in [(0, 0, '0'), (1, 1, '1'), (2, 4, '2-4'),
                         (5, 9, '5-9'), (10, 10 ** 6, '10+')]:
        m = tab[(tab['n_prom'] >= lo) & (tab['n_prom'] <= hi)]
        if not len(m):
            continue
        log(f"  {name:>8} {len(m):7d} {m['n_acc'].median():13.0f} "
            f"{m['n_acc'].mean():11.1f}")
    rho_p = tab[['n_prom', 'n_acc']].corr(method='spearman').iloc[0, 1]
    rho_s = tab[['sum_P', 'n_acc']].corr(method='spearman').iloc[0, 1]
    log(f"\n  rank correlation with the number of accepted lines: "
        f"n_prom {rho_p:.3f}, sum_P {rho_s:.3f}")
    n_low = int((tab['n_prom'] < MIN_PROM).sum())
    log(f"  {len(tab) - n_low} of the {len(tab)} found levels reach "
        f"n_prom >= {MIN_PROM}, and {n_low} do not - which is the bar itself "
        f"being checked:\n  a level that could not be found this way is "
        f"hardly ever a level that was found")
    log("\n  The found levels are not a fair sample - they were found because "
        "they had lines - so\n  this says the count ranks levels correctly, "
        "not that a level with n_prom = 10 will\n  certainly turn up.  What "
        "it rules out is the opposite error: a level the count calls\n  "
        "hopeless is not one that a search would have found anyway.")
    return tab


def measured_energies(en, ids=IDS, e_final=None, log=print):
    """``{IDEN2 row: measured energy}`` for the found levels.

    The energies of the run are preferred where the level's identifier is
    known, so that a partner sits exactly where the rest of the pipeline puts
    it; ``enlev.dat``'s own observed energy is the fallback, and it is the same
    number to a thousandth unless the run has moved the level since IDEN2 was
    last brought up to date.
    """
    id_of = cowan_gA.read_id_map(ids)
    out, n_run, moved = {}, 0, 0.0
    for i, k in en['known'].items():
        if not k:
            continue
        e_file = float(en['E_obs'][i])
        e_run = (e_final or {}).get(id_of.get(i))
        if e_run is None:
            out[i] = e_file
        else:
            out[i] = float(e_run)
            n_run += 1
            moved = max(moved, abs(e_run - e_file))
    log(f"  partner energies: {n_run} taken from the run, "
        f"{len(out) - n_run} from enlev.dat; the largest disagreement between "
        f"the two is {moved:.3f} cm^-1")
    return out


# ---------------------------------------------------------------------------
# The probability that each transition was recorded
# ---------------------------------------------------------------------------
def observation_probabilities(t, ctx, z_grid=Z_GRID, z_weight=Z_WEIGHT):
    """Add ``W``, ``P_obs`` and the quantities at E_calc to the frame ``t``.

    ``P_obs`` is averaged over the trial positions E_calc + z W with the
    weights of the normal density (section 3 of the module docstring), so it is
    the probability that the line is on the plates given both the physics and
    the ignorance about E.  ``wn``, ``lam``, ``I_pred`` and ``P_at_calc`` are
    the same quantities evaluated at E_calc alone, for the detail table.
    """
    C = float(cl.CFG.intensity_model['C'])
    kT = float(cl.CFG.intensity_model['kT'])
    W = t['W'].to_numpy(dtype=float)
    E_c = t['E_calc'].to_numpy(dtype=float)
    E_p = t['E_partner'].to_numpy(dtype=float)
    gA = t['gA'].to_numpy(dtype=float)

    def at(E):
        nu = np.abs(E - E_p)
        i_pred = C * gA * (nu / 1.0e8) * np.exp(-np.maximum(E, E_p) / kT)
        inside = (nu >= cl.WN_MIN) & (nu <= cl.WN_MAX)
        p = np.zeros(len(nu))
        if inside.any():
            p[inside] = lp.p_obs_vec(nu[inside], i_pred[inside],
                                     ctx.calib, ctx.bias)
        return nu, i_pred, p

    total = np.zeros(len(t))
    for z, wt in zip(z_grid, z_weight):
        total += wt * at(E_c + z * W)[2]

    nu, i_pred, p_calc = at(E_c)
    out = t.copy()
    out['wn'] = nu
    out['lam'] = np.where(nu > 0, 1.0e8 / np.where(nu > 0, nu, 1.0), np.nan)
    out['I_pred'] = i_pred
    out['P_at_calc'] = p_calc
    out['P_obs'] = total
    # how crowded the neighbourhood is: recorded lines per cm^-1 there
    out['rho'] = np.interp(nu, ctx.wn_o, ctx.rho_o)
    return out


# ---------------------------------------------------------------------------
# The table
# ---------------------------------------------------------------------------
def level_table(t, en, p_prom=P_PROMISING):
    """One row per unfound level, best first.

    ``n_prom`` - promising transitions, P_obs >= ``p_prom`` - is the verdict;
    everything else is context for reading it.
    """
    rows = []
    for idx, g in t.groupby('idx'):
        p = g['P_obs'].to_numpy(dtype=float)
        prom = g[p >= p_prom]
        strongest = g.loc[g['I_pred'].idxmax()] if len(g) else None
        rows.append(dict(
            idx=int(idx),
            label=en['label'][idx],
            J=en['J'][idx],
            parity=li.configuration_parity(en['cfg'][idx]),
            cfg=en['cfg'][idx],
            E_calc=float(en['E_calc'][idx]),
            W=float(g['W'].iloc[0]),
            W_src=g['W_src'].iloc[0],
            n_pred=len(g),
            n_obs=int((p >= lp.P_SEEN).sum()),
            n_prom=len(prom),
            sum_P=float(p.sum()),
            I_max=float(strongest['I_pred']),
            wn_I_max=float(strongest['wn']),
            rho_med=float(np.median(prom['rho'])) if len(prom)
            else float(np.median(g['rho'])),
        ))
    tab = pd.DataFrame(rows)
    if tab.empty:
        return tab
    # rounded here rather than on the way out, so that what is printed and
    # what is written are the same numbers
    for col, nd in dict(E_calc=1, W=1, sum_P=2, I_max=1, wn_I_max=2,
                        rho_med=3).items():
        tab[col] = tab[col].round(nd)
    return tab.sort_values(['n_prom', 'sum_P'], ascending=False)\
              .reset_index(drop=True)


def attach_partner_names(t, en, ids=IDS):
    """Give every row the partner's IDEN2 label and its pipeline level_id.

    Both are wanted for reading a detail table beside IDEN2: the label is what
    IDEN2 shows on the screen, the level_id what every other file of the
    pipeline calls the same level.
    """
    id_of = cowan_gA.read_id_map(ids)
    out = t.copy()
    out['partner_label'] = [en['label'][i] for i in t['partner']]
    out['partner_id'] = [id_of.get(int(i), '') for i in t['partner']]
    return out


def attach_windows(t, en, per_cfg, whole):
    """Give every row the W of its unfound level's configuration."""
    cfg = np.array([en['cfg'][i] for i in t['idx']])
    out = t.copy()
    out['W'] = [per_cfg.get(c, whole) for c in cfg]
    out['W_src'] = ['cfg' if c in per_cfg else 'list' for c in cfg]
    return out


def summarise(tab, log=print):
    """What the table says, in words."""
    if tab.empty:
        log("no unfound level has a transition to a found one")
        return
    n = len(tab)
    searchable = tab[tab['n_prom'] >= MIN_PROM]
    log(f"\n{n} unfound levels have at least one calculated transition to a "
        f"found level.")
    log(f"  {int((tab['n_prom'] == 0).sum())} have no transition that would "
        f"probably have been recorded at all - nothing to search on")
    log(f"  {int((tab['n_prom'] == 1).sum())} have exactly one, which is not "
        f"enough: one line can be made to fit any energy")
    log(f"  {len(searchable)} have {MIN_PROM} or more, and are the levels "
        f"worth a search; {int((tab['n_prom'] >= 5).sum())} have five or more")
    if len(searchable):
        by_cfg = searchable.groupby('cfg').size().sort_values(ascending=False)
        log("  by configuration: "
            + ', '.join(f"{c} {v}" for c, v in by_cfg.items()))
        guessed = searchable[searchable['W_src'] == 'list']
        if len(guessed):
            log(f"  {len(guessed)} of them are in a configuration with no "
                f"found level, so the width their probabilities were averaged "
                f"over is the list-wide guess")


def print_detail(t, tab, idx, p_prom=P_PROMISING, show_all=False, log=print):
    """The transition table of one unfound level."""
    g = t[t['idx'] == idx]
    if g.empty:
        log(f"IDEN2 row {idx} has no calculated transition to a found level")
        return
    head = tab[tab['idx'] == idx].iloc[0]
    log(f"\nIDEN2 row {idx}: {head['label']}   J = {head['J']:g}   "
        f"{'odd' if head['parity'] == 'o' else 'even'} parity")
    log(f"  calculated energy {head['E_calc']:.1f} cm^-1, and levels of "
        f"{head['cfg']} are found within {head['W']:.1f} cm^-1 of their "
        f"calculated position"
        + ('' if head['W_src'] == 'cfg' else
           ' (no found level of this configuration: the list-wide width)'))
    log(f"  {head['n_pred']} transitions to found levels, {head['n_prom']} of "
        f"them promising (P_obs >= {p_prom:g}), expected number recorded "
        f"{head['sum_P']:.1f}")
    show = g if show_all else g[g['P_obs'] >= lp.P_SEEN]
    show = show.sort_values('P_obs', ascending=False)
    log(f"\n  {'partner':>7} {'label':<11} {'level_id':<13} {'E_partner':>11} "
        f"{'L is':>5} {'wn':>10} {'lambda':>9} {'I_pred':>10} {'P_obs':>6} "
        f"{'rho':>6}")
    for r in show.itertuples(index=False):
        log(f"  {r.partner:7d} {r.partner_label:<11} {r.partner_id:<13} "
            f"{r.E_partner:11.3f} "
            f"{('upper' if r.upper == 'L' else 'lower'):>5} {r.wn:10.2f} "
            f"{r.lam:9.1f} {r.I_pred:10.2f} {r.P_obs:6.2f} {r.rho:6.3f}")
    if not show_all and len(show) < len(g):
        log(f"  ... and {len(g) - len(show)} transitions too faint to have "
            f"been recorded (P_obs below {lp.P_SEEN:g}); --detail-all lists "
            f"them")
    log("\n  'L is' says whether the unfound level is the upper or the lower "
        "of the two.\n  wn and lambda are where the transition would fall if "
        "the level sat exactly at its\n  calculated energy; P_obs is averaged "
        f"over the {head['W']:.0f} cm^-1 the level may be out by."
        "\n  rho is the density "
        "of recorded lines there, per cm^-1 - the higher it is, the more\n  "
        "easily a line will be found near any trial energy whether the level "
        "is there or not.")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def build(args, log=print):
    """Everything the report needs: the per-transition frame and the table."""
    ctx = lp.build(args, log=(lambda *a: None) if args.quiet_run else log)
    en, trans, mapping = read_theory(args.enlev, args.ids, log)
    per_cfg, whole = windows(en, log)
    e_meas = measured_energies(en, args.ids, ctx.e_final, log)
    t = unfound_transitions(en, trans, mapping, e_meas, log)
    t = attach_partner_names(t, en, args.ids)
    t = attach_windows(t, en, per_cfg, whole)
    t = observation_probabilities(t, ctx)
    tab = level_table(t, en, args.p_prom)
    return ctx, en, t, tab, dict(trans=trans, mapping=mapping,
                                 per_cfg=per_cfg, whole=whole, e_meas=e_meas)


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--detail', type=int, metavar='IDEN2_ROW',
                   help='the per-transition table of one unfound level, by '
                        'its row number in IDEN2/enlev.dat')
    p.add_argument('--detail-all', action='store_true',
                   help='--detail lists every transition, including the ones '
                        'too faint to have been recorded')
    p.add_argument('--validate', action='store_true',
                   help='score the count on the levels that have been found, '
                        'by computing it for them as though they had not')
    p.add_argument('--min-prom', type=int, default=MIN_PROM,
                   help='print only the levels with at least this many '
                        'promising transitions (default %(default)s); the csv '
                        'keeps every level whatever this says')
    p.add_argument('--p-prom', type=float, default=P_PROMISING,
                   help='P_obs at which a transition counts as promising '
                        '(default %(default)s)')
    p.add_argument('--top', type=int, default=40,
                   help='rows to print (default %(default)s); 0 for all')
    p.add_argument('--enlev', default=ENLEV)
    p.add_argument('--ids', default=IDS)
    p.add_argument('--out', default=OUT_CSV,
                   help='output csv (an .xlsx twin is written beside it)')
    p.add_argument('--quiet-run', action='store_true',
                   help='do not print what the run measures')
    # the arguments level_positions.build reads off the namespace
    p.add_argument('--lopt', default=None,
                   help='build the run from a LOPT line-output file')
    p.add_argument('--lopt-levels', default=None)
    p.add_argument('--energies', default='refit', choices=['refit', 'lopt'])
    p.add_argument('--energies-csv', default=None)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if not args.detail and not args.validate:
        output_files.require_writable(output_files.with_twin(args.out),
                                      'report file')
    ctx, en, t, tab, extra = build(args)

    if args.detail:
        print_detail(t, tab, args.detail, args.p_prom, args.detail_all)
        return 0

    if args.validate:
        validate(ctx, en, extra['trans'], extra['mapping'], extra['per_cfg'],
                 extra['whole'], extra['e_meas'], args.p_prom, args.ids)
        return 0

    summarise(tab)
    mc.save_table(tab, args.out)
    print(f"\nwritten to {os.path.basename(args.out)} and its .xlsx twin")

    show = tab[tab['n_prom'] >= args.min_prom]
    if args.top:
        show = show.head(args.top)
    with pd.option_context('display.width', 200, 'display.max_columns', 24):
        print(f"\nthe {len(show)} unfound levels most likely to be found:")
        print(show.to_string(index=False, na_rep='-'))
    print("\nn_prom is the count that matters: transitions to a found level "
          "that would probably have\nbeen recorded (P_obs >= "
          f"{args.p_prom:g}).  sum_P is the expected number of the level's "
          "lines actually\non the plates.  n_obs uses the weaker bar of "
          f"level_positions.py (P_obs >= {lp.P_SEEN:g}).\nrho_med is how "
          "crowded the promising transitions' neighbourhoods are, in recorded "
          "lines\nper cm^-1: where it is high, a trial energy will find lines "
          "whether the level is there or not.")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
