from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class EnergyLevel:
    level_id: str  # Unique identifier from Excel
    energy: float
    parity: str
    J_str: str
    J_val: float


@dataclass
class SpectralLine:
    wavenumber: float
    wn_uncertainty: float   # Uncertainty of wavenumber
    intensity: float        # Observed intensity reduced to a uniform scale (>0)
    line_character: str
    assigned_transitions: List['Transition'] = field(default_factory=list)
    original_assignments: List['Transition'] = field(default_factory=list)


@dataclass
class Transition:
    lower_level: Optional[EnergyLevel] = None
    upper_level: Optional[EnergyLevel] = None
    calc_intensity: Optional[float] = None    # On the same scale as observed intensities in lines (>=0); None for transitions without known theoretical intensity
    u_calc: Optional[float] = None            # Uncertainty of calculated intensity on log scale
    CF: Optional[float] = None                # Retained for future use
    assigned_to: Optional['SpectralLine'] = None  # Can be None
    grade: Optional[str] = None               # Can be None
    notes: str = ""                           # Metadata: 'N' (New), 'F' (Conflicting), 'R' (Revised)

    @property
    def calculated_wavenumber(self) -> float:
        if self.lower_level and self.upper_level:
            return self.upper_level.energy - self.lower_level.energy
        return 0.0

# Special instance for rejected classifications
UNASSIGNED = Transition(notes='R')
