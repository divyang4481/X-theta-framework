"""
X-Theta experiment metadata schema.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Optional

@dataclass
class SpacetimeMetadata:
    """
    Metadata required for physical X-Theta predictions.
    """
    emitter_coords: Optional[Dict[str, float]] = None
    detector_coords: Optional[Dict[str, float]] = None
    altitude: Optional[float] = None
    path_length: Optional[float] = None
    gravitational_potential: Optional[float] = None
    satellite_params: Optional[Dict[str, float]] = None

    def has_sufficient_metadata(self) -> bool:
        """Checks if enough metadata exists for a physical prediction."""
        # Simple heuristic for now
        return all([self.emitter_coords, self.detector_coords, self.gravitational_potential])
