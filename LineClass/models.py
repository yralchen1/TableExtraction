from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class EnergyLevel:
    level_id: str  # Unique identifier from Excel
    energy: float
    parity: str
    J_str: str
    J_val: float
    is_new: int = 0                        # 1 = level absent from ASD (found by Wyart), 0 = previously known
    is_decoy: int = 0                      # 1 = decoy (shadow) copy used for false-positive calibration
    is_added: int = 0                      # 1 = level found since the adopted level list was published,
                                           #     read from files.new_levels rather than from the workbook
    iden2_row: int = 0                     # row of IDEN2/enlev.dat this level is, 0 = unknown; set for added levels
    cowan_lid: int = 0                     # level number of the Cowan calculation (tp_E1_no_trials.xlsx),
                                           #     0 = unknown; how an added level's calculated transitions are found
    intens_from_factor: float = 0.0        # weighted mean of ln(I_obs/I_calc) for transitions FROM this level (upper)
    u_intens_from_factor: float = 0.0      # uncertainty (chi²-inflated)
    intens_to_factor: float = 0.0          # same, for transitions TO this level (lower)
    u_intens_to_factor: float = 0.0
    u_energy: float = 0.0                  # estimated uncertainty of the adopted energy (cm^-1); 0.0 = not yet estimated
    hfs_S: float = 0.0                     # I*A*J, the offset of the F = I+J sublevel from the center of
                                           #     gravity (cm^-1); 0.0 unless [hfs] apply is on
                                           #     (see hfs_correction.py)
    u_hfs_S: float = 0.0                   # its uncertainty (cm^-1)
    hfs_undetermined: bool = False         # listed in the table of A constants without a usable
                                           #     A: u_hfs_S is the unknown A, S is 0
    hfs_S_pattern: float = 0.0             # the S of a level of [hfs] resolved_levels in a line that
                                           #     does not separate its sublevels (hfs_S is 0 for it);
                                           #     hfs_S for every other level
    u_hfs_S_pattern: float = 0.0           # its uncertainty (cm^-1)
    from_transitions: List['Transition'] = field(default_factory=list)  # transitions where this is upper_level
    to_transitions: List['Transition'] = field(default_factory=list)    # transitions where this is lower_level
    u_contrib: Dict[tuple, tuple] = field(default_factory=dict)
    # The individual energy determinations behind u_energy, keyed by
    # trans_key(transition) -> (implied energy, its uncertainty); filled by
    # compute_level_uncertainties() so that a single transition's own
    # contribution can be removed again (see u_energy_excluding()).


@dataclass
class SpectralLine:
    wavenumber: float
    wn_uncertainty: float   # Uncertainty of wavenumber
    intensity: float        # Observed intensity reduced to a uniform scale (>0)
    line_character: str
    wn_key: float = 0.0
    # The immutable name of this observed line: the wavenumber it is filed
    # under in the hand-kept ledgers, whatever scale `wavenumber` is on.  It
    # comes from the column named by `[lines.layout.columns] wn_key`, which
    # defaults to the `wn` column, so a baseline run has wn_key == wavenumber.
    # A corrected set reads its wavenumbers from another column and leaves
    # this one on Sugar's original value, and then `line_decisions.csv`,
    # `new_levels.txt`, `revised_level_energies.csv` and
    # `inflated_unc_lines.txt` are the same files for both sets: a line's
    # identity does not change when the scale it is measured on does.
    assigned_transitions: List['Transition'] = field(default_factory=list)
    original_assignments: List['Transition'] = field(default_factory=list)
    decisions: Dict[tuple, tuple] = field(default_factory=dict)
    # The manual verdicts on this line, (lower_id, upper_id) ->
    # ('accept'|'reject', reason), read from the decision ledger by
    # classify_lines.attach_line_decisions().  Empty unless the ledger names
    # this line.
    hfs_factor: float = 0.0
    # 1 - kappa, the fraction of a transition's hyperfine displacement D by
    # which this line was measured below the head-frame Ritz wavenumber;
    # kappa is the measurement convention of the line's class (flag, era,
    # character).  0.0 unless [hfs] apply is on, and then 0.0 for a line
    # Sugar flagged, which sits on the Ritz value (see hfs_correction.py).
    hfs_u_kappa: float = 0.0
    # the uncertainty of kappa.
    hfs_head_pairs: frozenset = frozenset()
    # The (lower_id, upper_id) pairs whose resolved hfs companions the
    # registry files.hfs_satellites lists: for those this line is the
    # strongest component of the pattern and sits on the head-frame Ritz
    # value, kappa = 1, whatever its class (hfs_correction.py).  Set by
    # classify_lines.attach_hfs_satellites().
    hfs_exceptions: Dict[tuple, tuple] = field(default_factory=dict)
    # (lower_id, upper_id) -> (kappa, u_kappa, unresolved level ids): the
    # transitions of this line the registry files.kappa_exceptions says were
    # measured otherwise than the line's class (hfs_correction.py).  Set by
    # classify_lines.attach_kappa_exceptions().
    hfs_companion: Optional[tuple] = None
    # (main line, lower level, upper level, rung) if the registry names this
    # line as a resolved hfs companion of that transition on the main line
    # (main line None if the head of the pattern is not observed): it is
    # classified as that, never accepted, and no candidate is sought for it.
    # None for every other line.
    hfs_blend_companion: Optional[tuple] = None
    # The same tuple if the registry names this line as a companion blended
    # with transitions of its own (column blend): it is classified as any
    # other line, and the hfs component takes its share of the line's
    # calculated intensity (classify_lines.calc_weights).  None otherwise.
    unc_before_hfs_allowance: float = 0.0
    # The line list's own uncertainty, kept when a registry row tagged hfs
    # (inflated_unc_lines.txt) widened it and [hfs] apply is on, so that the
    # widening can be withdrawn once the correction accounts for the line's
    # hfs (classify_lines.release_hfs_allowances).  0.0 = not widened so.
    fit_uncertainty: float = 0.0
    # The uncertainty the last level fit weighted this line with: its own,
    # and for a blend the uncertainty of the centroid owed to the calculated
    # intensities in quadrature (classify_lines.calc_weights,
    # blend_centroid.py).  0.0 = not weighted (no accepted transition).
    legacy_keys: Optional[set] = None
    # The (lower_id, upper_id) pairs this line holds a published
    # identification for, as classify_lines.retag_legacy_identifications()
    # works them out: the identifications of the line workbook, with a level
    # whose measured position was exchanged with another's replaced by that
    # other and an identification naming a re-positioned level dropped.  These
    # are the pairs the classification treats as old.  None means the step was
    # never run, and then the transitions seeded from the workbook are taken
    # as the old ones.


@dataclass
class Transition:
    lower_level: Optional[EnergyLevel] = None
    upper_level: Optional[EnergyLevel] = None
    calc_intensity: Optional[float] = None    # On the same scale as observed intensities in lines (>=0); None for transitions without known theoretical intensity
    u_calc: Optional[float] = None            # Uncertainty of calculated intensity on log scale
    CF: Optional[float] = None                # Retained for future use
    assigned_to: Optional['SpectralLine'] = None  # Can be None
    grade: Optional[str] = None               # Can be None
    notes1: str = ""                          # Metadata: 'F' (Conflicting), 'R' (Revised)
    notes2: str = ""                          # Reason for accepting or rejecting
    new: Optional[int] = None                 # 1 = new classification, 0 = original, None = unclassified
    accepted: Optional[int] = None            # 1 = accepted, 0 = rejected, None = undecided
    orig_calc_intensity: Optional[float] = None   # original I_calc from input, never modified
    orig_u_calc: Optional[float] = None           # original u_calc from input, never modified
    manual: str = ""                          # 'accept' / 'reject' if the decision ledger rules on
                                              #     this assignment, '' if it does not
    is_imputed: int = 0                       # 1 = the pair is absent from the calculated-transition
                                              #     file, so its intensity was imputed from the
                                              #     printing cutoff of gA (see gA_imputation.py);
                                              #     0 = the intensity comes from the file, or none

    @property
    def calculated_wavenumber(self) -> float:
        if self.lower_level and self.upper_level:
            return self.upper_level.energy - self.lower_level.energy
        return 0.0

    @property
    def hfs_D(self) -> float:
        """D = S(upper) - S(lower): how far the strongest hyperfine component
        lies from the center of gravity of the line (cm^-1)."""
        if self.lower_level and self.upper_level:
            return self.upper_level.hfs_S - self.lower_level.hfs_S
        return 0.0

    def hfs_D_on(self, line: Optional['SpectralLine'] = None) -> float:
        """D of this transition as `line` (default: the line it is assigned
        to) measures it: hfs_D, unless the registry files.kappa_exceptions
        names levels of [hfs] resolved_levels this line does not resolve,
        which then count with their whole pattern (hfs_S_pattern)."""
        line = self.assigned_to if line is None else line
        exc = line.hfs_exceptions.get(
            (self.lower_level.level_id, self.upper_level.level_id)) \
            if line is not None and line.hfs_exceptions else None
        if not exc or not exc[2]:
            return self.hfs_D
        low, upp = self.lower_level, self.upper_level
        s_u = upp.hfs_S_pattern if upp.level_id in exc[2] else upp.hfs_S
        s_l = low.hfs_S_pattern if low.level_id in exc[2] else low.hfs_S
        return s_u - s_l

    def hfs_offset(self, line: Optional['SpectralLine'] = None) -> float:
        """How far below the Ritz wavenumber `line` (default: the line this
        transition is assigned to) is expected to have been measured, if it
        is this transition: (1 - kappa) * D.  0.0 with [hfs] apply off, and
        0.0 for a transition whose resolved companions the registry of
        files.hfs_satellites lists on this line.  A transition the registry
        files.kappa_exceptions names on this line takes the kappa written
        there, and the D of hfs_D_on."""
        line = self.assigned_to if line is None else line
        if line is None:
            return 0.0
        if line.hfs_exceptions:
            exc = line.hfs_exceptions.get((self.lower_level.level_id,
                                           self.upper_level.level_id))
            if exc is not None:
                return (1.0 - exc[0]) * self.hfs_D_on(line)
        if not line.hfs_factor:
            return 0.0
        if line.hfs_head_pairs and (self.lower_level.level_id,
                                    self.upper_level.level_id) \
                in line.hfs_head_pairs:
            return 0.0
        return line.hfs_factor * self.hfs_D

    def predicted_for(self, line: Optional['SpectralLine'] = None) -> float:
        """The wavenumber at which `line` (default: the line this transition
        is assigned to) would have been measured if it is this transition:
        the Ritz wavenumber in the head frame, less (1 - kappa) * D.  With
        [hfs] apply off it is the Ritz wavenumber itself."""
        return self.calculated_wavenumber - self.hfs_offset(line)

    @property
    def observed_head(self) -> float:
        """The observed wavenumber of the line this transition is assigned
        to, carried into the head frame for this transition: the measured
        value plus (1 - kappa) * D.  It is what the level energies are fitted
        to.  With [hfs] apply off it is the measured value itself."""
        return self.assigned_to.wavenumber + self.hfs_offset()

# Special instance for rejected classifications
UNASSIGNED = Transition(notes1='R')
