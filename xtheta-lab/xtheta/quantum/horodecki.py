"""
X-Theta Horodecki and state diagnostics.
"""
from __future__ import annotations
import numpy as np
from xtheta.quantum.chsh import s_max_horodecki

def get_concurrence(Phi: float) -> float:
    """Returns the concurrence C(Phi) = |cos(2Phi)|."""
    return float(abs(np.cos(2 * Phi)))

def get_horodecki_smax(Phi: float) -> float:
    """Returns S_max = 2 * sqrt(1 + cos^2(2Phi))."""
    c = get_concurrence(Phi)
    return 2.0 * np.sqrt(1.0 + c**2)
