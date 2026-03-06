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


@dataclass
class Transition:
    lower_level: EnergyLevel
    upper_level: EnergyLevel
    calc_intensity: Optional[float] = None    # On the same scale as observed intensities in lines (>=0); None for transitions without known theoretical intensity
    CF: Optional[float] = None                # Cancellation factor in theoretical calculation of intensities (>=0); None for transitions without known theoretical intensity
    assigned_to: Optional['SpectralLine'] = None  # Can be None
    grade: Optional[str] = None               # Can be None

    @property
    def calculated_wavenumber(self) -> float:
        return self.upper_level.energy - self.lower_level.energy