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
    intensity: float
    line_character: str
    possible_transitions: List['Transition'] = None


@dataclass
class Transition:
    upper_level: EnergyLevel
    lower_level: EnergyLevel

    @property
    def calculated_wavenumber(self) -> float:
        return self.upper_level.energy - self.lower_level.energy