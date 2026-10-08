"""Magnetic-dipole hyperfine constants A from the eigenvectors of the Cowan fit.

WHAT IT COMPUTES
================
The hyperfine constant A of a level is the expectation value of the
magnetic-dipole hyperfine operator, divided by J.  In the effective-operator
form of Sandars and Beck that operator is, for every open shell nl,

    T(1) = sum_i [ a01 l_i  -  sqrt(10) a12 (s_i C2_i)(1)  +  a10 s_i ]

summed over the electrons i of that shell.  The three numbers a01 (orbital),
a12 (spin-dipolar) and a10 (contact, or core polarization for l > 0) are the
*radial parameters* of the shell, in cm^-1.  Nonrelativistically a01 = a12 and
a10 = 0 for l > 0; a single s electron has only a10.  The checks: one p, d or f
electron gives A_j = a l(l+1)/[j(j+1)] with a01 = a12 = a, and one s electron
gives A = a10.

So A of every level is linear in the radial parameters,

    A = sum over shells and kinds of  theta(shell, kind) * a(shell, kind),

and the angular coefficients theta are what the eigenvector decides.  This
program computes the thetas of every level, evaluates A with a set of radial
parameters (``hfs_radial_params.toml``), and can fit chosen parameters to the
measured constants of ``A_hfs_levels.csv``.

WHERE THE EIGENVECTORS COME FROM
================================
``RCEOUT``, the output of Cowan's RCE fit, lists for every level up to ten
leading components: the signed amplitude times 100, rounded to an integer,
the configuration, and a 12-character label of the LS basis state.  The
signs matter: A is a quadratic form in the amplitudes, and its cross terms
between basis states of different S or different parent terms are as large as
its diagonal ones.

The basis states are rebuilt here exactly as RCG builds them:

* the states of one shell l^w come from Cowan's own coefficients of
  fractional parentage (``CowanLSShell`` of the Cowan repository's
  ``CODE/lsjj/cowan_cfp.py``, reading ``CODE/ING11.CFP``), so they carry
  RCG's alpha labels and RCG's signs;
* shells are coupled in the order of the configuration card, the running L
  and S first and the new shell second, as RCG's CALCFC recouples them; L and
  S are coupled to J with L first;
* the label of each state is formed by RCG's own rule (CALCFC, rcg11kd.f,
  the loop that sets ISUBJ and FORMAT 41), so an RCEOUT label names exactly
  one state.

The operator is applied to these states in the m-scheme - Slater
determinants, no Racah algebra - so no phase convention has to be
reconciled.  Matrix elements between configurations are neglected, as is
usual for effective radial parameters.

Three configurations are left out: ``f4p5``, ``p5f3d`` and ``p5f3s``, the
5p^5 core-excited ones.  Their labels do not identify a state (RCG
renames duplicates), there are no radial parameters for a 5p hole, and
together they hold well under 1 per cent of the observed levels'
compositions.  Their components are counted in ``w_skipped``.

HOW GOOD THE TEN COMPONENTS ARE
================================
Two errors come from RCEOUT itself, and both are computed per level.

* Rounding: an amplitude printed as an integer per cent is uncertain by
  +-0.005 (a uniform distribution, standard deviation 0.0029).  ``u_round``
  propagates that through the quadratic form.
* Truncation: the components after the tenth are missing.  Their total
  weight is ``1 - norm``, where ``norm`` is the sum of the squares of the
  printed amplitudes.  ``u_trunc`` is the largest change in A that weight
  could make (Cauchy-Schwarz): twice the missing amplitude times the part of
  T|psi> lying outside the printed states, plus the missing weight times the
  largest eigenvalue of T in the block.  It is a bound, not a standard
  deviation.

The printed amplitudes are normalized to unit length before A is formed.

HOW GOOD THE EIGENVECTOR IS
===========================
Both errors above measure only how well the file reproduces Cowan's own
eigenvector, not how close that eigenvector is to the true one.  The
amplitudes of a level-structure fit are themselves uncertain, by about
0.05-0.10 in the user's experience (2026-10-07).  ``u_amp`` gives each
printed amplitude an independent error of standard deviation ``--amp-sigma``
(default ``AMP_SIGMA`` = 0.10) and propagates it to A to first order, exactly
as ``u_round`` propagates the rounding: u_amp = u_round * amp_sigma /
U_ROUND.  A Monte Carlo of the same model (random amplitudes, renormalized,
A recomputed) agrees with the first-order value to about 10 per cent.  It is
a lower bound on the real error: mixing with a state that is not printed is
not in it.  u_amp is not part of ``z`` or of the fit weights.

Two more terms come from the radial parameters, both set in
``hfs_radial_params.toml``:

* ``u_par``: the parameters' own uncertainties (``[uncertainty]``), taken
  as independent, times dA/da;
* ``u_cfg``: a configuration is *tested* when at least ``min_measured``
  levels with a measured constant have it as their leading configuration
  (column ``tested``).  Nothing checks that a parameter fitted in one
  configuration holds in another, so the part of A the untested
  configurations give is uncertain by ``fraction`` of itself.

OUTER SHELLS NO MEASUREMENT REACHES
===================================
The measured constants test the 4f, 5d, 6s (and, by estimate, 6p)
parameters only.  The outer electrons of the high configurations - 4f2 7s,
8s, 6d, 7d, 7p, 5f, 6f, 5g, 6g - get theirs from the table ``[scaled]`` of
the parameter file (2026-10-08): each is a copy of a fitted shell times a
factor read from the Cowan fit itself, RCEOUT's parameter listing
(``read_rceout_parameters``).  For an s shell the contact parameter goes as
1/n*^3 (Fermi-Segre), with n* from the shell's E_av and the limit of the
4f2 ns series (``s_series``: limit 186714 cm^-1, one quantum defect 3.405,
fitted to 4f2 6s, 7s, 8s).  For l > 0 the factor is the ratio of the
spin-orbit parameters zeta, which goes as <r^-3> just as a does.  The
factor's relative uncertainty (``fraction``) enters u_par.  Without these
the 271 observed levels of those configurations had no A at all.

``u_total`` combines u_round, u_trunc, u_amp, u_par and u_cfg in
quadrature.  It is the uncertainty to give a calculated A that is used
(``hfs_A_candidates.py``); such rows of ``A_hfs_levels.csv`` carry a source
beginning with CALC_SOURCE, and the fit here never uses them.

LEVEL IDENTIFIERS
=================
RCEOUT carries no level identifier.  Its levels are paired with the levels of
``tp_E1_no_trials.xlsx`` - the same calculation - by aligning the two
calculated-energy lists within each parity and J (``cowan_gA._align``), and
the Cowan level number of each pair is turned into a level_id through
``classify_lines.cowan_lid_ids``, which goes through
``IDEN2/IDEN_level_ids.txt``.  Where both files carry an observed energy the
two are compared, and the pairs that disagree by more than ``OBS_TOL`` are
reported.

USAGE
=====
    python hfs_A_theory.py                 thetas and A for every identified
                                           level, written to hfs_A_theory.csv
    python hfs_A_theory.py --fit 5d.a01,5d.a10
                                           also fit those radial parameters
                                           to the measured A constants
    python hfs_A_theory.py --amp-sigma 0.05
                                           u_amp for amplitude errors of 0.05
    python hfs_A_theory.py --no-write      print the report only

Only hfs_A_candidates.py reads the output; no A constant used by the
pipeline changes.
"""

import argparse
import copy
import csv
import math
import os
import re
import sys
import tomllib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
COWAN_REPO = os.environ.get('COWAN_REPO', r'F:\COWAN\repo')
RCEOUT_FILE = os.path.join(HERE, 'RCEOUT')
PARAMS_FILE = os.path.join(HERE, 'hfs_radial_params.toml')
A_LEVELS_FILE = os.path.join(HERE, 'A_hfs_levels.csv')
OUT_CSV = os.path.join(HERE, 'hfs_A_theory.csv')
OUT_LOG = os.path.join(HERE, 'hfs_A_theory.log')

KINDS = ('a01', 'a12', 'a10')
LSYM = 'SPDFGHIKLMNOQRTUVWXYZABC'      # LSYMB of rcg11kd.f
LETTER_L = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5}
U_ROUND = 0.005 / math.sqrt(3.0)        # sd of an amplitude printed to 0.01
AMP_SIGMA = 0.10                        # sd of an amplitude of the fit itself
CALC_SOURCE = 'calculated'              # source prefix of A_hfs_levels.csv
                                        # rows taken from this program
OBS_TOL = 0.5                           # cm^-1, observed-energy agreement
SKIP = ('f4p5', 'p5f3d', 'p5f3s')
MISSING_MAX = 0.05                      # weight without parameters that
                                        # still allows an A
TAIL_FROM = 7                           # components kept when the tail
                                        # estimate is calibrated

# The configuration cards of the Pr III deck (IN36 of the Cowan repository's
# WORK/Pr3, columns 33 on), in the order of the deck.  RCG numbers the
# subshells of a configuration by their position on this card, and both the
# coupling order and the label rule go by that position.
CARDS = """
4f3    4f3  5p6  5d0
f26p   4f2  5p6  5d0  6s0  6p
f27p   4f2  5p6  5d0  6s0  6p0  6d0  7s0  7p
f25f   4f2  5p6  5d0  6s0  6p0  6d0  7s0  5f
f26f   4f2  5p6  5d0  6s0  6p0  6d0  7s0  6f
f5d2   4f   5p6  5d2
f5d6s  4f   5p6  5d   6s
f5d7s  4f   5p6  5d   6s0  6p0  6d0  7s0  7s
f5d6d  4f   5p6  5d   6s0  6p0  6d
f5d7d  4f   5p6  5d   6s0  6p0  6d0  7s0  7d
f6s2   4f   5p6  5d0  6s2
f6d2   4f   5p6  5d0  6s0  6p0  6d2
f6s7s  4f   5p6  5d0  6s   6p0  6d0  7s0  7s
f7s2   4f   5p6  5d0  6s0  6p0  6d0  7s2
f6p2   4f   5p6  5d0  6s0  6p2
f6p7p  4f   5p6  5d0  6s0  6p   6d0  7s0  7p
d6s6p  4f0  5p6  5d   6s   6p
d26p   4f0  5p6  5d2  6s0  6p
6s26p  4f0  5p6  5d0  6s2  6p
6p3    4f0  5p6  5d0  6s0  6p3
f4p5   4f4  5p5
f26s   4f2  5p6  5d0  6s   6p0
f27s   4f2  5p6  5d0  6s0  6p0  6d0  7s
f28s   4f2  5p6  5d0  6s0  6p0  6d0  7s0  8s
f25d   4f2  5p6  5d
f26d   4f2  5p6  5d0  6s0  6p0  6d
f27d   4f2  5p6  5d0  6s0  6p0  6d0  7s0  7d
f25g   4f2  5p6  5d0  6s0  6p0  6d0  7s0  5g
f26g   4f2  5p6  5d0  6s0  6p0  6d0  7s0  6g
fd6p   4f   5p6  5d   6s0  6p
fd7p   4f   5p6  5d   6s0  6p0  6d0  7s0  7p
5d3    4f0  5p6  5d3
5d26s  4f0  5p6  5d2  6s
5d26d  4f0  5p6  5d2  6s0  6p0  6d
5d6s2  4f0  5p6  5d   6s2
5d6p2  4f0  5p6  5d   6s0  6p2
6s6p2  4f0  5p6  5d0  6s   6p2
p5f3d  4f3  5p5  5d
p5f3s  4f3  5p5  5d0  6s
"""


def parse_card(text):
    """`'4f2 5p6 6p'` -> [('4f', 3, 2), ('5p', 1, 6), ('6p', 1, 1)]."""
    out = []
    for tok in text.split():
        m = re.fullmatch(r'(\d+)([spdfgh])(\d*)', tok)
        if m is None:
            raise ValueError('cannot read subshell %r' % tok)
        out.append((m.group(1) + m.group(2), LETTER_L[m.group(2)],
                    int(m.group(3)) if m.group(3) else 1))
    return out


def read_cards(text=CARDS):
    out = {}
    for line in text.strip().splitlines():
        name, rest = line.split(None, 1)
        out[name] = parse_card(rest)
    return out


CONFIGS = read_cards()


# ---------------------------------------------------------------------------
# the Cowan repository: RCEOUT reader and the CFP-built shell states
# ---------------------------------------------------------------------------
_COWAN = {}


def cowan():
    """`(rce, cowan_cfp, lsjj, decks)` from the Cowan repository.

    The lsjj directory goes on the path before CODE/, because CODE/lsjj is
    also a package of that name and cowan_cfp imports the module inside it.
    Appended, not inserted, so nothing here shadows a module of LineClass.
    """
    if not _COWAN:
        code = os.path.join(COWAN_REPO, 'CODE')
        if not os.path.isdir(code):
            raise SystemExit('The Cowan repository is not at %s; set '
                             'COWAN_REPO to where it is.' % COWAN_REPO)
        for p in (os.path.join(code, 'lsjj'), code):
            if p not in sys.path:
                sys.path.append(p)
        from toolbox import rce
        import cowan_cfp
        import lsjj
        _COWAN.update(rce=rce, cowan_cfp=cowan_cfp, lsjj=lsjj,
                      decks=cowan_cfp.read_deck_file(
                          os.path.join(code, 'ING11.CFP')))
    c = _COWAN
    return c['rce'], c['cowan_cfp'], c['lsjj'], c['decks']


def cg(j1, m1, j2, m2, j, m):
    """Clebsch-Gordan coefficient, every argument doubled."""
    return cowan()[2].cg(j1, m1, j2, m2, j, m)


def threej(j1, m1, j2, m2, j3, m3):
    return cowan()[2].threej(j1, m1, j2, m2, j3, m3)


# ---------------------------------------------------------------------------
# the one-electron operators
# ---------------------------------------------------------------------------
def orbitals(l):
    """The one-electron states (2ml, 2ms) in the order LSShell uses."""
    return [(ml2, ms2) for ml2 in range(2 * l, -2 * l - 1, -2)
            for ms2 in (1, -1)]


def _spin(q, ms2p, ms2):
    """< ms' | s(1)_q | ms >."""
    if q == 0:
        return ms2 / 2.0 if ms2p == ms2 else 0.0
    if q == 1:
        return -1.0 / math.sqrt(2.0) if (ms2p, ms2) == (1, -1) else 0.0
    return 1.0 / math.sqrt(2.0) if (ms2p, ms2) == (-1, 1) else 0.0


def _ck(l, k, q, ml2p, ml2):
    """< l ml' | C(k)_q | l ml >."""
    return ((-1) ** (ml2p // 2) * (2 * l + 1)
            * threej(2 * l, 0, 2 * k, 0, 2 * l, 0)
            * threej(2 * l, -ml2p, 2 * k, 2 * q, 2 * l, ml2))


def one_electron(l, kind):
    """The q = 0 component of one kind of the operator, over `orbitals(l)`.

    a01: l_z.  a10: s_z.  a12: -sqrt(10) (s(1) C(2))(1)_0.
    """
    orbs = orbitals(l)
    m = np.zeros((len(orbs), len(orbs)))
    for i, (lp, sp) in enumerate(orbs):
        for j, (lq, sq) in enumerate(orbs):
            if kind == 'a01':
                v = lq / 2.0 if i == j else 0.0
            elif kind == 'a10':
                v = sq / 2.0 if i == j else 0.0
            else:
                v = 0.0
                for q1 in (-1, 0, 1):
                    c = cg(2, 2 * q1, 4, -2 * q1, 2, 0)
                    if c:
                        v += c * _spin(q1, sp, sq) * _ck(l, 2, -q1, lp, lq)
                v *= -math.sqrt(10.0)
            m[i, j] = v
    return m


def _sign(mask, p):
    """(-1) to the number of occupied orbitals below p."""
    return -1.0 if bin(mask & ((1 << p) - 1)).count('1') % 2 else 1.0


class Shell:
    """The LS states of l^w, Cowan's, and the operators acting on them."""

    _cache = {}

    def __new__(cls, l, w):
        key = (l, w)
        if key not in cls._cache:
            self = super().__new__(cls)
            self._build(l, w)
            cls._cache[key] = self
        return cls._cache[key]

    def _build(self, l, w):
        _rce, cowan_cfp, _lsjj, decks = cowan()
        self.l, self.w = l, w
        built = cowan_cfp.CowanLSShell(l, w, deck_table=decks)
        # terms: (2L, 2S, alpha, {(2ML, 2MS): {mask: amplitude}})
        self.terms = []
        for (L2, S2), states in sorted(built.terms.items()):
            for alpha, _suffix, _degen, cols in states:
                self.terms.append((L2, S2, alpha, cols))
        self.ops = {k: one_electron(l, k) for k in KINDS}
        self._act = {k: {} for k in KINDS}

    def act(self, kind, mask):
        """{mask': coefficient} of one operator kind on one determinant."""
        cache = self._act[kind]
        if mask in cache:
            return cache[mask]
        m = self.ops[kind]
        out = {}
        for q in range(m.shape[0]):
            if not mask >> q & 1:
                continue
            s_q = _sign(mask, q)
            rest = mask & ~(1 << q)
            for p in np.nonzero(m[:, q])[0]:
                p = int(p)
                if rest >> p & 1:
                    continue
                new = rest | (1 << p)
                out[new] = out.get(new, 0.0) + s_q * _sign(rest, p) * m[p, q]
        cache[mask] = out
        return out


# ---------------------------------------------------------------------------
# the basis of one configuration
# ---------------------------------------------------------------------------
def rcg_label_position(card):
    """RCG's ISUBJ for a configuration card (rcg11kd.f, CALCFC, loop 380).

    Negative: the label's parent is the term of that subshell, the lowest
    one occupied by more than one electron (or hole).  Positive: the label's
    parent is the running term after that subshell, the last singly occupied
    one before the very last.  Positions count from 1.
    """
    isubj = isub = 0
    for i1 in range(len(card), 0, -1):
        _nl, l, w = card[i1 - 1]
        full = 4 * l + 2
        if w == 0 or w > full - 1:
            continue
        if w == 1 or w == full - 1:
            if isubj == 0 and i1 < isub:
                isubj = i1
            isub = i1
        else:
            isubj = -i1
    return isubj or 1


def norm_label(text):
    """A basis-state label with the parentheses and blanks removed."""
    return re.sub(r'[()\s]', '', text)


class ConfigBasis:
    """The LS basis states of one configuration, and the operator in it.

    `shells` holds the open subshells in card order as (position, nl, Shell).
    A basis state is a chain, one (term index, 2L, 2S) per open subshell,
    the last two being the running L and S after that subshell.
    """

    def __init__(self, name, card=None):
        self.name = name
        self.card = CONFIGS[name] if card is None else card
        self.shells = [(i + 1, nl, Shell(l, w))
                       for i, (nl, l, w) in enumerate(self.card)
                       if 0 < w < 4 * l + 2]
        self.parity = sum(l * w for _nl, l, w in self.card) % 2
        self.chains = self._chains()
        self.labels = {ch: self._label(ch) for ch in self.chains}
        self._memo = {}
        self._blocks = {}

    def _chains(self):
        sh = self.shells[0][2]
        chains = [((t, L2, S2),) for t, (L2, S2, _a, _c) in enumerate(sh.terms)]
        for _pos, _nl, sh in self.shells[1:]:
            new = []
            for ch in chains:
                Lr, Sr = ch[-1][1], ch[-1][2]
                for t, (L2, S2, _a, _c) in enumerate(sh.terms):
                    for L in range(abs(Lr - L2), Lr + L2 + 1, 2):
                        for S in range(abs(Sr - S2), Sr + S2 + 1, 2):
                            new.append(ch + ((t, L, S),))
            chains = new
        return chains

    def _label(self, chain):
        """RCG's label of a chain, with parentheses and blanks removed."""
        isubj = rcg_label_position(self.card)
        pos = abs(isubj)
        term = (0, 0, '')               # (2L, 2S, alpha) at that position
        run = (0, 0)
        for k, (p, _nl, sh) in enumerate(self.shells):
            if p > pos:
                break
            L2, S2, alpha, _c = sh.terms[chain[k][0]]
            run = (chain[k][1], chain[k][2])
            if p == pos:
                term = (L2, S2, alpha)
        if isubj < 0:
            mult, L = term[1] + 1, term[0]
        else:
            mult, L = run[1] + 1, run[0]
        Lt, St = chain[-1][1], chain[-1][2]
        return '%d%s%d%s%s' % (mult, LSYM[L // 2], St + 1, LSYM[Lt // 2],
                               term[2])

    def block(self, J2):
        """`(chains, labels)` of the basis states of total 2J."""
        out = [ch for ch in self.chains
               if abs(ch[-1][1] - ch[-1][2]) <= J2 <= ch[-1][1] + ch[-1][2]
               and (ch[-1][1] + ch[-1][2] - J2) % 2 == 0]
        return out, [self.labels[ch] for ch in out]

    # -- m-scheme states ----------------------------------------------------
    def _comp(self, chain, ML2, MS2):
        """{(mask, ...): amplitude} of a chain prefix at (2ML, 2MS)."""
        key = (chain, ML2, MS2)
        if key in self._memo:
            return self._memo[key]
        k = len(chain) - 1
        t, L, S = chain[-1]
        L2t, S2t, _a, cols = self.shells[k][2].terms[t]
        out = {}
        if k == 0:
            for mask, v in cols.get((ML2, MS2), {}).items():
                out[(mask,)] = v
        else:
            Lr, Sr = chain[-2][1], chain[-2][2]
            for mlb in range(-L2t, L2t + 1, 2):
                mla = ML2 - mlb
                if abs(mla) > Lr:
                    continue
                a = cg(Lr, mla, L2t, mlb, L, ML2)
                if not a:
                    continue
                for msb in range(-S2t, S2t + 1, 2):
                    msa = MS2 - msb
                    if abs(msa) > Sr:
                        continue
                    b = cg(Sr, msa, S2t, msb, S, MS2)
                    if not b:
                        continue
                    right = cols.get((mlb, msb), {})
                    for lm, lv in self._comp(chain[:-1], mla, msa).items():
                        for rm, rv in right.items():
                            kk = lm + (rm,)
                            out[kk] = out.get(kk, 0.0) + a * b * lv * rv
        self._memo[key] = out
        return out

    def state(self, chain, J2):
        """{(mask, ...): amplitude} of a basis state at M = J."""
        L, S = chain[-1][1], chain[-1][2]
        out = {}
        for ML2 in range(-L, L + 1, 2):
            MS2 = J2 - ML2
            if abs(MS2) > S:
                continue
            c = cg(L, ML2, S, MS2, J2, J2)
            if not c:
                continue
            for kk, v in self._comp(chain, ML2, MS2).items():
                out[kk] = out.get(kk, 0.0) + c * v
        return {k: v for k, v in out.items() if abs(v) > 1e-12}

    def _apply(self, state, k, kind):
        sh = self.shells[k][2]
        out = {}
        for kk, v in state.items():
            for nm, c in sh.act(kind, kk[k]).items():
                key = kk[:k] + (nm,) + kk[k + 1:]
                out[key] = out.get(key, 0.0) + c * v
        return out

    def matrices(self, J2):
        """{(nl, kind): matrix of <i|T_z|j>/J} over `block(J2)` states.

        Divided by J, so that c.M.c is directly the coefficient theta of the
        radial parameter in A.  Empty for J = 0, which has no A.
        """
        if J2 in self._blocks:
            return self._blocks[J2]
        chains, _labels = self.block(J2)
        out = {}
        if J2 > 0 and chains:
            states = [self.state(ch, J2) for ch in chains]
            for k, (_pos, nl, sh) in enumerate(self.shells):
                for kind in KINDS:
                    if sh.l == 0 and kind != 'a10':
                        continue
                    m = np.zeros((len(chains), len(chains)))
                    for j, sj in enumerate(states):
                        img = self._apply(sj, k, kind)
                        for i, si in enumerate(states):
                            m[i, j] = _dot(si, img)
                    out[(nl, kind)] = 0.5 * (m + m.T) / (J2 / 2.0)
        self._blocks[J2] = out
        return out


def _dot(a, b):
    if len(a) > len(b):
        a, b = b, a
    return sum(v * b.get(k, 0.0) for k, v in a.items())


_BASES = {}


def config_basis(name):
    if name not in _BASES:
        _BASES[name] = ConfigBasis(name)
    return _BASES[name]


# ---------------------------------------------------------------------------
# radial parameters
# ---------------------------------------------------------------------------
class Params:
    """Radial parameters: `value(config, nl, kind)` in cm^-1, or None.

    Read from a TOML table `[radial]` whose keys are a shell (`"5d"`) or a
    configuration and a shell (`"f5d2/5d"`, which overrides the plain one in
    that configuration), each holding up to three kinds.  A kind may be a
    number, or the name of another kind of the same entry (`a12 = "a01"`),
    which ties the two.  A kind not given has no value, and a level that
    needs it gets no A.

    The table `[scaled]` adds outer shells no measured constant reaches,
    each a copy of a shell of `[radial]` times a factor taken from the Cowan
    fit (`scale_outer_shells`).  They are made again by `apply_scaling`
    whenever a parameter they copy changes.
    """

    def __init__(self, table, scalings=None):
        self.table = {k: dict(v) for k, v in table.items()}
        for key, entry in self.table.items():
            for kind, v in entry.items():
                if kind not in KINDS:
                    raise ValueError('%s: unknown kind %r' % (key, kind))
                if isinstance(v, str) and v not in KINDS:
                    raise ValueError('%s.%s: %r is not a kind' % (key, kind, v))
        self.scaled = {}            # {nl: Scaling}
        for sc in (scalings or []):
            if sc.nl in self.table:
                raise ValueError('[scaled] %s is also in [radial]' % sc.nl)
            if sc.ref not in self.table:
                raise ValueError('[scaled] %s: %s is not in [radial]'
                                 % (sc.nl, sc.ref))
            self.scaled[sc.nl] = sc
        self.apply_scaling()

    @classmethod
    def read(cls, path=PARAMS_FILE, rceout=RCEOUT_FILE):
        with open(path, 'rb') as fh:
            doc = tomllib.load(fh)
        scalings = None
        if doc.get('scaled'):
            scalings = scale_outer_shells(doc['scaled'],
                                          read_rceout_parameters(rceout))
        return cls(doc.get('radial', {}), scalings)

    def apply_scaling(self):
        """(Re)make every `[scaled]` entry from the shell it copies: numbers
        times the factor, ties kept."""
        for nl, sc in self.scaled.items():
            self.table[nl] = {
                kind: v if isinstance(v, str) else sc.factor * float(v)
                for kind, v in self.table[sc.ref].items()}

    def entry(self, config, nl):
        key = '%s/%s' % (config, nl)
        if key in self.table:
            return key
        return nl if nl in self.table else None

    def resolve(self, key, kind):
        """The (entry, kind) a value is actually stored under."""
        v = self.table[key].get(kind)
        if isinstance(v, str):
            return key, v
        return key, kind

    def value(self, config, nl, kind):
        key = self.entry(config, nl)
        if key is None:
            return None
        key, kind = self.resolve(key, kind)
        v = self.table[key].get(kind)
        return None if v is None or isinstance(v, str) else float(v)

    def stored(self, config, nl, kind):
        """The (entry, kind) the value of `kind` of `nl` in `config` comes
        from, or None where there is no value."""
        if self.value(config, nl, kind) is None:
            return None
        return self.resolve(self.entry(config, nl), kind)

    def set(self, name, value):
        key, kind = name.rsplit('.', 1)
        self.table.setdefault(key, {})[kind] = float(value)

    def names(self):
        """Every stored (not tied, not scaled) parameter as
        `"entry.kind"`."""
        return ['%s.%s' % (k, kind) for k, e in self.table.items()
                if k not in self.scaled
                for kind, v in e.items() if not isinstance(v, str)]


# ---------------------------------------------------------------------------
# outer shells scaled from the Cowan fit
# ---------------------------------------------------------------------------
Z_CORE = 3                  # charge of the Pr IV core an outer electron sees
RYDBERG = 109736.9          # cm^-1, for the mass of 141Pr


class Scaling:
    """An outer shell `nl` whose radial parameters are those of `ref` times
    `factor`; `fraction` is the factor's relative uncertainty, `how` says
    where the factor came from and `alt` is the n*^-3 estimate beside it
    (None where there is none)."""

    def __init__(self, nl, ref, factor, fraction, how, alt=None):
        self.nl, self.ref, self.factor = nl, ref, factor
        self.fraction, self.how, self.alt = fraction, how, alt


def read_rceout_parameters(path=RCEOUT_FILE):
    """{configuration: (E_av, {shell position: zeta})} from the parameter
    listing of RCEOUT, in cm^-1 (RCEOUT prints 1000 cm^-1).  The position is
    RCG's: the 1-based place of the shell on the configuration card
    (`CONFIGS`), which ZETA's index gives.  A later listing replaces an
    earlier one."""
    out, conf = {}, None
    eav = re.compile(r'EAV (\S+)\s+-?\d+\s+(-?\d+\.\d+)')
    zeta = re.compile(r'ZETA\s*(\d+)\s+-?\d+\s+(-?\d+\.\d+)')
    with open(path, encoding='ascii', errors='replace') as fh:
        for line in fh:
            m = eav.match(line)
            if m:
                conf = m.group(1)
                out[conf] = (1000.0 * float(m.group(2)), {})
                continue
            m = zeta.match(line)
            if m and conf is not None:
                out[conf][1][int(m.group(1))] = 1000.0 * float(m.group(2))
    return out


def outer_config(nl):
    """The configuration 4f2 nl of the deck: 4f2, a closed 5p6, nl once,
    every other shell empty; None if there is none."""
    for name, shells in CONFIGS.items():
        occ = {s: w for s, _l, w in shells if w}
        if occ == {'4f': 2, '5p': 6, nl: 1}:
            return name
    return None


def _zeta(rce, conf, nl):
    """zeta of shell `nl` in `conf`, cm^-1, or None."""
    if conf not in rce:
        return None
    for pos, (s, _l, _w) in enumerate(CONFIGS[conf], start=1):
        if s == nl:
            return rce[conf][1].get(pos)
    return None


def s_series(rce):
    """The 4f2 ns series of the deck: `(limit, delta, {nl: n*}, {nl:
    residual})`, every E_av fitted by E_av = limit - Z_CORE^2 R/(n -
    delta)^2 with one quantum defect delta (least squares; exact for two
    members)."""
    from scipy.optimize import minimize_scalar
    E = {}
    for nl in ('%ds' % n for n in range(5, 12)):
        conf = outer_config(nl)
        if conf in rce:
            E[nl] = (int(nl[:-1]), rce[conf][0])
    if len(E) < 2:
        raise ValueError('RCEOUT has fewer than two 4f2 ns configurations')
    zr = Z_CORE ** 2 * RYDBERG

    def limit(d):
        return float(np.mean([e + zr / (n - d) ** 2 for n, e in E.values()]))

    def ss(d):
        L = limit(d)
        return sum((e - (L - zr / (n - d) ** 2)) ** 2 for n, e in E.values())

    nmin = min(n for n, _e in E.values())
    d = minimize_scalar(ss, bounds=(0.0, nmin - 0.5), method='bounded',
                        options={'xatol': 1e-10}).x
    L = limit(d)
    return (L, d, {nl: math.sqrt(zr / (L - e)) for nl, (n, e) in E.items()},
            {nl: e - (L - zr / (n - d) ** 2) for nl, (n, e) in E.items()})


def nstar(rce, nl, series):
    """n* of the outer electron of 4f2 nl, from its E_av and the limit of
    the ns series, or None."""
    conf = outer_config(nl)
    if conf not in rce or series[0] <= rce[conf][0]:
        return None
    return math.sqrt(Z_CORE ** 2 * RYDBERG / (series[0] - rce[conf][0]))


def scale_outer_shells(table, rce):
    """`Scaling`s for the `[scaled]` table of the parameter file, each entry
    `nl = {from = ref, fraction = f}`.

    * An s shell: the contact parameter goes as 1/n*^3 (Fermi-Segre), so
      factor = (n*_ref / n*_nl)^3, n* from the E_av of 4f2 ref and 4f2 nl
      and the limit of the ns series (`s_series`).
    * Any other shell: a01 and a10 go as <r^-3>, and so, nearly, does the
      spin-orbit parameter zeta, so factor = zeta(nl) / zeta(ref), each in
      its own 4f2 nl configuration; a ref of 4f takes the zeta of the 4f
      core of 4f2 nl.  The n*^-3 ratio is kept in `alt` as a check where
      ref is an outer shell too.
    """
    series = s_series(rce)
    out = []
    for nl, spec in table.items():
        ref = spec['from']
        fraction = float(spec.get('fraction', 0.0))
        l_nl, l_ref = LETTER_L[nl[-1]], LETTER_L[ref[-1]]
        if (l_nl == 0) != (l_ref == 0):
            raise ValueError('[scaled] %s: an s shell scales only from an '
                             's shell' % nl)
        conf = outer_config(nl)
        if conf is None or conf not in rce:
            raise ValueError('[scaled] %s: RCEOUT has no 4f2 %s' % (nl, nl))
        n_nl, n_ref = nstar(rce, nl, series), nstar(rce, ref, series)
        alt = (n_ref / n_nl) ** 3 if n_nl and n_ref else None
        if l_nl == 0:
            if alt is None:
                raise ValueError('[scaled] %s: no n* for %s or %s'
                                 % (nl, nl, ref))
            out.append(Scaling(nl, ref, alt, fraction,
                               '(n* %s / n* %s)^3 = (%.3f / %.3f)^3'
                               % (ref, nl, n_ref, n_nl)))
            continue
        conf_ref = conf if ref == '4f' else outer_config(ref)
        z_nl, z_ref = _zeta(rce, conf, nl), _zeta(rce, conf_ref, ref)
        if not z_nl or not z_ref:
            raise ValueError('[scaled] %s: RCEOUT gives no zeta for %s in %s '
                             'or %s in %s' % (nl, nl, conf, ref, conf_ref))
        out.append(Scaling(nl, ref, z_nl / z_ref, fraction,
                           'zeta %s (%s) / zeta %s (%s) = %.1f / %.1f'
                           % (nl, conf, ref, conf_ref, z_nl, z_ref),
                           None if ref == '4f' else alt))
    return out


def read_uncertainties(path=PARAMS_FILE, params=None):
    """`(sigmas, fraction, min_measured)` from the parameter file.

    `sigmas` is {(entry, kind): standard uncertainty} of the table
    `[uncertainty]`; `fraction` and `min_measured` are those of `[untested]`
    (0 and 1 when it is absent).  With `params`, every scaled entry gets an
    uncertainty too: the copied parameter's, times the factor, and the
    factor's own `fraction` of the value, in quadrature (the two are taken
    as independent, though the first is shared with the shell copied).
    """
    with open(path, 'rb') as fh:
        doc = tomllib.load(fh)
    sigmas = {}
    for key, entry in doc.get('uncertainty', {}).items():
        for kind, v in entry.items():
            if kind not in KINDS:
                raise ValueError('[uncertainty] %s: unknown kind %r'
                                 % (key, kind))
            sigmas[(key, kind)] = float(v)
    for nl, sc in (params.scaled if params is not None else {}).items():
        for kind, v in params.table[nl].items():
            if isinstance(v, str):
                continue
            sigmas[(nl, kind)] = math.hypot(
                sc.factor * sigmas.get((sc.ref, kind), 0.0),
                sc.fraction * float(v))
    untested = doc.get('untested', {})
    return (sigmas, float(untested.get('fraction', 0.0)),
            int(untested.get('min_measured', 1)))


# ---------------------------------------------------------------------------
# levels
# ---------------------------------------------------------------------------
class LevelTheta:
    """The thetas of one RCEOUT level, and what limits them."""

    def __init__(self, level, kset, parity):
        self.level = level
        self.kset = kset
        self.parity = parity
        self.J = level.j
        self.e_calc = level.e_calc * 1000.0
        self.e_obs = level.e_obs * 1000.0 if level.observed else None
        self.level_id = ''
        self.lid = None
        self.norm = sum((c.percent / 100.0) ** 2 for c in level.components)
        self.w_skipped = 0.0
        self.unmatched = []
        # {config: (indices into the block, amplitudes)}, the amplitudes
        # divided by sqrt(norm), and {config: their weight}
        self.parts = {}
        self.weights = {}
        self.theta = {}         # {(config, nl, kind): theta}
        self._compute()

    def _compute(self):
        J2 = int(round(2 * self.J))
        by_conf = {}
        for c in self.level.components:
            by_conf.setdefault(c.configuration, []).append(c)
        n = math.sqrt(self.norm) if self.norm > 0 else 1.0
        for conf, comps in by_conf.items():
            if conf in SKIP or conf not in CONFIGS:
                self.w_skipped += sum((c.percent / 100.0) ** 2 for c in comps)
                continue
            basis = config_basis(conf)
            _chains, labels = basis.block(J2)
            index = {}
            for i, lab in enumerate(labels):
                index.setdefault(lab, []).append(i)
            idx, amp = [], []
            for c in comps:
                hit = index.get(norm_label(c.label), [])
                if len(hit) != 1:
                    self.unmatched.append('%s %s' % (conf, c.label.strip()))
                    self.w_skipped += (c.percent / 100.0) ** 2
                    continue
                idx.append(hit[0])
                amp.append(c.percent / 100.0 / n)
            if not idx:
                continue
            self.parts[conf] = (np.array(idx), np.array(amp))
            self.weights[conf] = float(np.sum(np.array(amp) ** 2))
            if J2 == 0:
                continue
            for (nl, kind), m in basis.matrices(J2).items():
                ii, a = self.parts[conf]
                self.theta[(conf, nl, kind)] = float(a @ m[np.ix_(ii, ii)] @ a)

    def A(self, params):
        """`(A, missing, w_missing)` with `params`.

        A shell without a value of a kind the level needs contributes
        nothing.  `missing` names those parameters and `w_missing` is the
        weight of the configurations they belong to; A is None when that
        weight exceeds MISSING_MAX.
        """
        total, missing, confs = 0.0, set(), set()
        for (conf, nl, kind), th in self.theta.items():
            if abs(th) < 1e-9:
                continue
            v = params.value(conf, nl, kind)
            if v is None:
                missing.add('%s.%s' % (nl, kind))
                confs.add(conf)
                continue
            total += th * v
        w_missing = sum(self.weights[c] for c in confs)
        return ((None if w_missing > MISSING_MAX else total), sorted(missing),
                w_missing)

    def _operator(self, conf, params):
        """The block matrix of T/J with `params`, or None."""
        J2 = int(round(2 * self.J))
        mats = config_basis(conf).matrices(J2)
        if not mats:
            return None
        out = None
        for (nl, kind), m in mats.items():
            v = params.value(conf, nl, kind)
            if not v:
                continue
            out = v * m if out is None else out + v * m
        return out

    def uncertainties(self, params, A, miss=None):
        """`(u_round, bound, cancel)` with `params`.

        u_round is what the printed precision of the amplitudes allows.
        bound is the largest change in A that a missing part of the
        eigenvector of weight `miss` (by default all of 1 - norm) could make:
        2 sqrt(miss) |P T psi| + miss |T|, with P the projection on the basis
        states not printed and |T| the largest eigenvalue of T in the block.
        cancel is the sum of the absolute contributions divided by |A|: 1
        where every term has the same sign, large where A is the residue of
        a cancellation.
        """
        images, total = [], 0.0
        for conf, (ii, a) in self.parts.items():
            T = self._operator(conf, params)
            if T is None:
                continue
            Ta = T[:, ii] @ a                       # T|psi> over the block
            images.append((T, ii, a, Ta))
            total += float(a @ Ta[ii])
        var, img_out, lam = 0.0, 0.0, 0.0
        for T, ii, a, Ta in images:
            # with a normalized, d(a.T.a / a.a)/da_i = 2[(Ta)_i - A a_i]
            var += float(np.sum((2.0 * (Ta[ii] - total * a)) ** 2))
            outside = np.delete(Ta, ii)
            img_out += float(outside @ outside)
            lam = max(lam, float(np.max(np.abs(np.linalg.eigvalsh(T)))))
        u_round = U_ROUND * math.sqrt(var) / math.sqrt(max(self.norm, 1e-12))
        if miss is None:
            miss = max(1.0 - self.norm, 0.0)
        bound = 2.0 * math.sqrt(miss) * math.sqrt(img_out) + miss * lam
        terms = [abs(th * (params.value(c, nl, k) or 0.0))
                 for (c, nl, k), th in self.theta.items()]
        cancel = sum(terms) / abs(A) if A else float('inf')
        return u_round, bound, cancel

    def u_params(self, params, sigmas):
        """The uncertainty of A that the radial parameters' own give, to
        first order, the parameters taken as independent: dA/da is the sum
        of the thetas whose values a supplies, ties included."""
        deriv = {}
        for (conf, nl, kind), th in self.theta.items():
            src = params.stored(conf, nl, kind)
            if src is not None:
                deriv[src] = deriv.get(src, 0.0) + th
        return math.sqrt(sum((d * sigmas.get(src, 0.0)) ** 2
                             for src, d in deriv.items()))

    def untested_part(self, params, tested):
        """The part of A that configurations outside `tested` give."""
        total = 0.0
        for (conf, nl, kind), th in self.theta.items():
            v = params.value(conf, nl, kind)
            if conf not in tested and v is not None:
                total += th * v
        return total

    @property
    def leading_config(self):
        comps = self.level.components
        return comps[0].configuration if comps else ''

    def truncated(self, k):
        """The same level computed from its first k components only."""
        level = copy.copy(self.level)
        level.components = self.level.components[:k]
        out = LevelTheta(level, self.kset, self.parity)
        out.level_id, out.lid = self.level_id, self.lid
        return out

    def theta_by_shell(self):
        """{(nl, kind): theta} summed over configurations."""
        out = {}
        for (_c, nl, kind), th in self.theta.items():
            out[(nl, kind)] = out.get((nl, kind), 0.0) + th
        return out

    @property
    def leading(self):
        c = self.level.components[0] if self.level.components else None
        return '' if c is None else '%d %s %s' % (c.percent, c.configuration,
                                                   c.label.strip())


def read_levels(path=RCEOUT_FILE):
    """Every RCEOUT level as a LevelTheta."""
    rce = cowan()[0]
    out = []
    for cset in rce.read_rceout(path):
        parity = None
        for name in cset.configurations:
            if name in CONFIGS:
                parity = 'odd' if config_basis(name).parity else 'even'
                break
        for level in cset.levels:
            out.append(LevelTheta(level, cset.index, parity))
    return out


def attach_level_ids(levels, log=print):
    """Give each level its Cowan number and level_id, where it has them.

    Returns the list of (level, E_exp) pairs whose observed energies
    disagree by more than OBS_TOL.
    """
    import cowan_gA
    import classify_lines
    trans = cowan_gA.read_transitions(log=lambda m: None)
    calc = cowan_gA.levels(trans)
    ids = classify_lines.cowan_lid_ids(trans)
    bad = []
    blocks = {}
    for lv in levels:
        blocks.setdefault((lv.parity, lv.J), []).append(lv)
    for (parity, J), group in blocks.items():
        mine = sorted(group, key=lambda lv: lv.e_calc)
        theirs = calc[(calc.parity == parity) & (calc.J == J)]
        theirs = theirs.sort_values('E_calc').reset_index(drop=True)
        if theirs.empty:
            continue
        # One J at a time: two levels of different J that swap order
        # between the two runs would otherwise cross, and the alignment,
        # which cannot cross, would leave one of them unpaired.
        e_mine = np.array([lv.e_calc for lv in mine])
        pairs, spare_l, spare_r = cowan_gA._align(
            e_mine, np.zeros(len(mine)), theirs.E_calc.to_numpy(),
            np.zeros(len(theirs)))
        # An order reversal inside one J leaves two unpaired levels
        # side by side; pair them by energy, as cowan_gA does.
        for j in sorted(spare_r):
            free = [i for i in spare_l
                    if abs(e_mine[i] - theirs.E_calc[j]) <= cowan_gA.CROSS_TOL]
            if free:
                i = min(free, key=lambda i: abs(e_mine[i] - theirs.E_calc[j]))
                spare_l.remove(i)
                pairs.append((i, j))
        for i, j in pairs:
            row = theirs.iloc[j]
            lv = mine[i]
            lv.lid = int(row.lid)
            lv.level_id = ids.get(lv.lid, '') or ''
            if (lv.e_obs is not None and row.E_exp == row.E_exp
                    and abs(lv.e_obs - row.E_exp) > OBS_TOL):
                bad.append((lv, float(row.E_exp)))
    return bad


# ---------------------------------------------------------------------------
# measured constants, and the fit
# ---------------------------------------------------------------------------
SUPERSEDED_NAME = 'A_hfs_levels_superseded.csv'


def read_measured(path=A_LEVELS_FILE):
    """{level_id: (A, u_A, source)} of the measured constants.

    Semiempirical rows (composition, Reader & Sugar), undetermined ones,
    hand estimates, rows marked CONFLICT and rows taken from this program
    (CALC_SOURCE) are left out.  A measured row that a calculated one has
    replaced (`hfs_A_candidates.py --write`) is read back from
    SUPERSEDED_NAME beside `path`, so that the parameters stay fitted to
    every measurement; one the user dismissed (reason "dismissed ...") is
    not.
    """
    import hfs_kappa

    def usable(src):
        return not (hfs_kappa.is_semiempirical(src)
                    or src == 'not determined' or src.startswith('estimate')
                    or 'CONFLICT' in src or src.startswith(CALC_SOURCE))

    out = {}
    with open(path, encoding='utf-8', newline='') as fh:
        for row in csv.DictReader(fh):
            if usable(row['source']):
                out[row['level_id']] = (float(row['A_cm-1']),
                                        float(row['u_A']), row['source'])
    old = os.path.join(os.path.dirname(os.path.abspath(path)),
                       SUPERSEDED_NAME)
    if os.path.isfile(old):
        with open(old, encoding='utf-8', newline='') as fh:
            for row in csv.DictReader(fh):
                if (row['level_id'] not in out and usable(row['source'])
                        and not row['reason'].startswith('dismissed')):
                    out[row['level_id']] = (float(row['A_cm-1']),
                                            float(row['u_A']), row['source'])
    return out


def tail_calibration(levels, params, first=TAIL_FROM):
    """How the truncation bound compares with what truncation actually does.

    Every level that prints all ten components is recomputed from its first
    `first` only.  The change that brings back components first+1 to 10 is
    known exactly, and so is the bound for a missing part of their weight;
    their ratio, over all such levels, measures how pessimistic the bound
    is.  Returns `(median, 90th percentile, number of levels)`.
    """
    ratios = []
    for lv in levels:
        if len(lv.level.components) < 10 or lv.J <= 0:
            continue
        A10 = lv.A(params)[0]
        if A10 is None:
            continue
        short = lv.truncated(first)
        A_short = short.A(params)[0]
        if A_short is None:
            continue
        w = sum((c.percent / 100.0) ** 2
                for c in lv.level.components[first:10]) / short.norm
        _u, bound, _c = short.uncertainties(params, A_short, miss=w)
        if bound > 1e-6:
            ratios.append(abs(A10 - A_short) / bound)
    if not ratios:
        return 1.0, 1.0, 0
    ratios = np.array(ratios)
    return (float(np.median(ratios)), float(np.percentile(ratios, 90)),
            len(ratios))


def _design(lv, params, free):
    """`(x, fixed, w_missing)`: the row of the fit for one level.

    x holds the thetas of the free parameters, fixed the contribution of
    every other parameter, and w_missing the weight of the configurations
    that need a parameter with no value.
    """
    x = np.zeros(len(free))
    fixed, confs = 0.0, set()
    for (conf, nl, kind), th in lv.theta.items():
        if abs(th) < 1e-9:
            continue
        key = params.entry(conf, nl)
        if key is None:
            confs.add(conf)
            continue
        key, kind = params.resolve(key, kind)
        name = '%s.%s' % (key, kind)
        if name in free:
            x[free.index(name)] += th
            continue
        v = params.table[key].get(kind)
        if v is None or isinstance(v, str):
            confs.add(conf)
            continue
        fixed += th * float(v)
    return x, fixed, sum(lv.weights[c] for c in confs)


def prepare_free(params, free):
    """Make an entry for every free parameter the file does not have.

    A configuration-specific one (`4f3/4f.a01`) starts as a copy of the
    plain entry of its shell, ties included.
    """
    for name in free:
        key, kind = name.rsplit('.', 1)
        if kind not in KINDS:
            raise SystemExit('--fit %s: %r is not one of %s'
                             % (name, kind, ', '.join(KINDS)))
        if key.split('/')[-1] in params.scaled:
            raise SystemExit('--fit %s: %s is scaled from %s; fit that'
                             % (name, key, params.scaled[
                                 key.split('/')[-1]].ref))
        if key not in params.table:
            plain = key.split('/')[-1]
            params.table[key] = dict(params.table.get(plain, {}))
        if params.table[key].get(kind) is None or isinstance(
                params.table[key].get(kind), str):
            params.table[key][kind] = 0.0


def fit(levels, params, free, measured, tail=0.0, rounds=3):
    """Least-squares values of the `free` parameters.

    Each measured constant is weighted by 1/(u_A^2 + u_round^2 + u_tail^2),
    u_tail being `tail` times the truncation bound; the calculated
    uncertainties depend on the parameters, so the weights are refreshed for
    a few rounds.  Returns `(values, sigmas, chi2, dof, levels used)`;
    `params` is left holding the fit.
    """
    prepare_free(params, free)
    data = [lv for lv in levels if lv.level_id in measured and lv.J > 0]
    for _round in range(rounds):
        used, X, y, w = [], [], [], []
        for lv in data:
            A_meas, u_meas, _src = measured[lv.level_id]
            x, fixed, w_missing = _design(lv, params, free)
            if w_missing > MISSING_MAX:
                continue
            A_now = lv.A(params)[0] or 0.0
            u_round, bound, _c = lv.uncertainties(params, A_now)
            used.append(lv)
            X.append(x)
            y.append(A_meas - fixed)
            w.append(1.0 / (u_meas ** 2 + u_round ** 2 + (tail * bound) ** 2))
        X, y, w = np.array(X), np.array(y), np.array(w)
        cov = np.linalg.inv((X * w[:, None]).T @ X)
        sol = cov @ (X * w[:, None]).T @ y
        for name, v in zip(free, sol):
            params.set(name, v)
        params.apply_scaling()      # the scaled shells follow their refs
    resid = y - X @ sol
    chi2 = float(np.sum(w * resid ** 2))
    return sol, np.sqrt(np.diag(cov)), chi2, len(y) - len(free), used


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------
def shell_columns(levels):
    seen = set()
    for lv in levels:
        seen.update(k for k, v in lv.theta_by_shell().items() if abs(v) > 1e-9)
    order = {nl: i for i, nl in enumerate(
        ['4f', '5d', '6s', '6p', '7s', '8s', '6d', '7d', '5f', '6f', '5g',
         '6g', '7p', '5p'])}
    return sorted(seen, key=lambda k: (order.get(k[0], 99), KINDS.index(k[1])))


def tested_configs(levels, measured, min_measured):
    """{configuration: n} of the leading configurations of at least
    `min_measured` levels with a measured constant."""
    count = {}
    for lv in levels:
        if lv.level_id in measured and lv.J > 0:
            c = lv.leading_config
            count[c] = count.get(c, 0) + 1
    return {c: n for c, n in count.items() if n >= min_measured}


def rows_of(levels, params, measured, tail, amp_sigma=AMP_SIGMA,
            sigmas=None, tested=None, untested_fraction=0.0):
    """One output row per level.  `tested` None counts every configuration
    as tested (no u_cfg)."""
    cols = shell_columns(levels)
    out = []
    for lv in levels:
        A, missing, w_missing = lv.A(params)
        u_round, bound, cancel = lv.uncertainties(params, A or 0.0)
        u_calc = math.hypot(u_round, tail * bound)
        # the same first-order propagation as u_round, with sd amp_sigma
        u_amp = u_round * amp_sigma / U_ROUND
        u_par = lv.u_params(params, sigmas or {})
        u_cfg = (0.0 if tested is None else
                 untested_fraction * abs(lv.untested_part(params, tested)))
        u_total = math.sqrt(u_calc ** 2 + u_amp ** 2 + u_par ** 2
                            + u_cfg ** 2)
        is_tested = tested is None or lv.leading_config in tested
        meas = measured.get(lv.level_id)
        th = lv.theta_by_shell()
        row = {
            'level_id': lv.level_id, 'parity': lv.parity,
            'J': '%g' % lv.J,
            'E_obs': '' if lv.e_obs is None else '%.3f' % lv.e_obs,
            'E_calc': '%.1f' % lv.e_calc,
            'leading': lv.leading,
            'norm': '%.4f' % lv.norm,
            'w_skipped': '%.4f' % lv.w_skipped,
        }
        for nl, kind in cols:
            row['th_%s_%s' % (nl, kind)] = '%+.5f' % th.get((nl, kind), 0.0)
        z = ''
        if meas is not None and A is not None:
            z = '%+.1f' % ((meas[0] - A) / math.hypot(meas[1], u_calc))
        row.update({
            'A_calc': '' if A is None else '%+.4f' % A,
            'u_round': '%.4f' % u_round,
            'u_trunc': '%.4f' % (tail * bound),
            'u_amp': '%.4f' % u_amp,
            'u_par': '%.4f' % u_par,
            'u_cfg': '%.4f' % u_cfg,
            'u_total': '%.4f' % u_total,
            'tested': 'yes' if is_tested else 'no',
            'trunc_bound': '%.4f' % bound,
            'cancel': '' if A is None or not math.isfinite(cancel)
                      else '%.1f' % cancel,
            'missing': ' '.join(missing),
            'w_missing': '%.4f' % w_missing,
            'A_meas': '' if meas is None else '%+.4f' % meas[0],
            'u_meas': '' if meas is None else '%.4f' % meas[1],
            'resid': '' if meas is None or A is None
                     else '%+.4f' % (meas[0] - A),
            'z': z,
        })
        out.append(row)
    return out


def say_scaling(say, rce, params):
    """Report how the `[scaled]` shells were made."""
    L, delta, ns, resid = s_series(rce)
    say('')
    say('scaled shells ([scaled] of the parameter file), factors from the '
        'Cowan fit (RCEOUT):')
    say('  4f2 ns series: limit %.0f cm^-1, quantum defect %.3f; n* %s; '
        'E_av - fit %s cm^-1'
        % (L, delta, ', '.join('%s %.3f' % kv for kv in ns.items()),
           ', '.join('%s %+.0f' % kv for kv in resid.items())))
    for nl, sc in params.scaled.items():
        say('  %-3s from %-3s factor %.4f +- %.0f%%: %s%s'
            % (nl, sc.ref, sc.factor, 100 * sc.fraction, sc.how,
               '' if sc.alt is None else '; n*^-3 would give %.4f' % sc.alt))


def _level_id(text):
    """`'151'` or `'000151'` -> `'059003.000151'`."""
    text = text.strip()
    return text if '.' in text else '059003.%06d' % int(text)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='hyperfine A constants from the RCEOUT eigenvectors')
    ap.add_argument('--rceout', default=RCEOUT_FILE)
    ap.add_argument('--params', default=PARAMS_FILE)
    ap.add_argument('--fit', default='',
                    help='comma-separated parameters to fit to the measured '
                         'A constants, e.g. 5d.a01,5d.a10 or 4f3/4f.a01')
    ap.add_argument('--exclude', default='',
                    help='comma-separated level_ids (zero-padded or bare '
                         'numbers) whose measured A the fit leaves out')
    ap.add_argument('--all', action='store_true',
                    help='write every RCEOUT level, not only those with a '
                         'level_id')
    ap.add_argument('--amp-sigma', type=float, default=AMP_SIGMA,
                    help='standard deviation of an eigenvector amplitude, '
                         'propagated to u_amp (default %(default)s)')
    ap.add_argument('--no-write', action='store_true')
    args = ap.parse_args(argv)

    report = ['hfs_A_theory.py %s' % ' '.join(sys.argv[1:] if argv is None
                                             else argv), '']

    def say(msg=''):
        print(msg)
        report.append(msg)

    params = Params.read(args.params, args.rceout)
    levels = read_levels(args.rceout)
    bad = attach_level_ids(levels)
    known = [lv for lv in levels if lv.level_id]
    say('RCEOUT: %d levels; %d paired with a level_id'
        % (len(levels), len(known)))
    for lv, e_exp in bad:
        say('  observed energies disagree: %s J=%g RCEOUT %.3f, Cowan table '
            '%.3f' % (lv.level_id or 'lid %s' % lv.lid, lv.J, lv.e_obs, e_exp))
    unmatched = sorted({u for lv in known for u in lv.unmatched})
    if unmatched:
        say('  component labels with no single basis state: %s'
            % ', '.join(unmatched[:20]))
    measured = read_measured()

    if args.fit:
        free = [s.strip() for s in args.fit.split(',') if s.strip()]
        out = {_level_id(x) for x in args.exclude.split(',') if x.strip()}
        data = {k: v for k, v in measured.items() if k not in out}
        # The truncation estimate needs a complete parameter set, which
        # the fit may be what supplies: fit, calibrate, fit again.
        fit(known, params, free, data)
        tail = tail_calibration(known, params)[0]
        sol, sig, chi2, dof, used = fit(known, params, free, data, tail=tail)
        scale = math.sqrt(chi2 / dof) if dof > 0 else float('nan')
        say('')
        say('fit to %d measured constants, %d free; chi2/dof = %.2f'
            % (len(used), len(free), chi2 / dof if dof else float('nan')))
        if out:
            say('  left out: %s' % ', '.join(sorted(out)))
        for name, v, s in zip(free, sol, sig):
            say('  %-12s %+.5f +- %.5f  (+- %.5f scaled by sqrt(chi2/dof))'
                % (name, v, s, s * max(scale, 1.0)))

    say('')
    say('radial parameters used (cm^-1):')
    for key, entry in params.table.items():
        say('  %-8s %s%s' % (key, ', '.join(
            '%s = %s' % (k, v if isinstance(v, str) else '%+.5f' % v)
            for k, v in entry.items()),
            '   (scaled from %s)' % params.scaled[key].ref
            if key in params.scaled else ''))
    if params.scaled:
        say_scaling(say, read_rceout_parameters(args.rceout), params)

    median, p90, n = tail_calibration(known, params)
    say('')
    say('truncation: over %d levels, the change brought by components %d-10 '
        'is a median %.2f (90th percentile %.2f) of its bound; u_trunc is '
        'the bound times %.2f' % (n, TAIL_FROM + 1, median, p90, median))

    sigmas, fraction, min_measured = read_uncertainties(args.params, params)
    tested = tested_configs(known, measured, min_measured)
    say('tested configurations (leading configuration of >= %d measured '
        'constants): %s; the part of A from the others is uncertain by %.0f%%'
        % (min_measured, ', '.join('%s %d' % kv for kv in sorted(
            tested.items(), key=lambda kv: -kv[1])), 100 * fraction))
    rows = rows_of(levels if args.all else known, params, measured, median,
                   args.amp_sigma, sigmas, tested, fraction)
    have = [r for r in rows if r['z']]
    if have:
        res = np.array([float(r['resid']) for r in have])
        z = np.array([float(r['z']) for r in have])
        say('measured constants compared: %d; rms residual %.4f cm^-1; '
            'median |z| %.2f' % (len(have), math.sqrt(np.mean(res ** 2)),
                                 float(np.median(np.abs(z)))))
        far = sorted((r for r in have if abs(float(r['z'])) > 3),
                     key=lambda r: -abs(float(r['z'])))
        for r in far:
            say('  %s J=%-4s %-24s A_meas %s +- %s  A_calc %s  z %s  cancel %s'
                % (r['level_id'], r['J'], r['leading'], r['A_meas'],
                   r['u_meas'], r['A_calc'], r['z'], r['cancel']))
    obs = [r for r in rows if r['E_obs']]
    if obs:
        norms = np.array([float(r['norm']) for r in obs])
        u_tr = np.array([float(r['u_trunc']) for r in obs])
        u_rd = np.array([float(r['u_round']) for r in obs])
        say('observed levels: printed weight median %.3f, 5th percentile '
            '%.3f, lowest %.3f; u_round median %.4f; u_trunc median %.4f, '
            '90th percentile %.4f cm^-1'
            % (np.median(norms), np.percentile(norms, 5), norms.min(),
               np.median(u_rd), np.median(u_tr), np.percentile(u_tr, 90)))
        with_A = [r for r in obs if r['A_calc']]
        u_am = np.array([float(r['u_amp']) for r in with_A])
        n_sign = sum(abs(float(r['A_calc'])) < float(r['u_amp'])
                     for r in with_A)
        say('amplitude error %.2f: u_amp median %.4f, 90th percentile %.4f, '
            'largest %.4f cm^-1; |A_calc| < u_amp for %d of %d observed '
            'levels' % (args.amp_sigma, np.median(u_am),
                        np.percentile(u_am, 90), u_am.max(), n_sign,
                        len(with_A)))
        u_tot = np.array([float(r['u_total']) for r in with_A])
        say('u_total (u_round, u_trunc, u_amp, u_par, u_cfg in quadrature): '
            'median %.4f, 90th percentile %.4f cm^-1'
            % (np.median(u_tot), np.percentile(u_tot, 90)))

    if not args.no_write:
        import output_files
        output_files.require_writable([OUT_CSV, OUT_LOG])
        with open(OUT_CSV, 'w', encoding='utf-8', newline='') as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()),
                                lineterminator='\n')
            wr.writeheader()
            wr.writerows(rows)
        with open(OUT_LOG, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write('\n'.join(report) + '\n')
        print('wrote %s and %s' % (os.path.basename(OUT_CSV),
                                   os.path.basename(OUT_LOG)))


if __name__ == '__main__':
    main()
