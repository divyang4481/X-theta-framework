"""
X-Theta Horodecki and state diagnostics.
"""
from __future__ import annotations
import numpy as np
from xtheta.quantum.chsh import s_max_horodecki

def get_concurrence(phi: float) -> float:
    """Returns the concurrence C(phi) = |cos(2phi)|."""
    return float(abs(np.cos(2 * phi)))

def get_horodecki_smax(phi: float) -> float:
    """Returns S_max = 2 * sqrt(1 + cos^2(2phi))."""
    c = get_concurrence(phi)
    return 2.0 * np.sqrt(1.0 + c**2)
