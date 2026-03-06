from dataclasses import dataclass
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
    assigned_transitions: List['Transition'] = None


@dataclass
class Transition:
    lower_level: EnergyLevel
    upper_level: EnergyLevel
    calc_intensity: float     # On the same scale as observed intensities in lines (>=0); Should be None for transitions without known theoretical intensity
    CF: float                 # Cancellation factor in theoretical calculation of intensities (>=0); Should be None for transitions without known theoretical intensity
    assigned_to: SpectralLine # Can be None, default None
    grade: str            # Can be None, default None

    @property
    def calculated_wavenumber(self) -> float:
        return self.upper_level.energy - self.lower_level.energy