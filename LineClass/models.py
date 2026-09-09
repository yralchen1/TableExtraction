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
    intens_from_factor: float = 0.0        # weighted mean of ln(I_obs/I_calc) for transitions FROM this level (upper)
    u_intens_from_factor: float = 0.0      # uncertainty (chi²-inflated)
    intens_to_factor: float = 0.0          # same, for transitions TO this level (lower)
    u_intens_to_factor: float = 0.0
    u_energy: float = 0.0                  # estimated uncertainty of the adopted energy (cm^-1); 0.0 = not yet estimated
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
    assigned_transitions: List['Transition'] = field(default_factory=list)
    original_assignments: List['Transition'] = field(default_factory=list)
    decisions: Dict[tuple, tuple] = field(default_factory=dict)
    # The manual verdicts on this line, (lower_id, upper_id) ->
    # ('accept'|'reject', reason), read from the decision ledger by
    # classify_lines.attach_line_decisions().  Empty unless the ledger names
    # this line.
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

# Special instance for rejected classifications
UNASSIGNED = Transition(notes1='R')
