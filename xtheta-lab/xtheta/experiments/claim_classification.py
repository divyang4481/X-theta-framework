"""
X-Theta claim classification system.
"""
from __future__ import annotations

class ClaimLevel:
    THEOREM = "theorem"
    SIMULATION = "simulation"
    PHENOMENOLOGICAL_FIT = "phenomenological_fit"
    PHYSICAL_PREDICTION = "physical_prediction"

def classify_claim(has_spacetime_metadata: bool, is_synthetic: bool = False) -> str:
    """
    Classifies a scientific claim based on available metadata.
    """
    if is_synthetic:
        return ClaimLevel.SIMULATION
    if has_spacetime_metadata:
        return ClaimLevel.PHYSICAL_PREDICTION
    return ClaimLevel.PHENOMENOLOGICAL_FIT
