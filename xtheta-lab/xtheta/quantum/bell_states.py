"""
X-Theta Bell state definitions.
"""
from __future__ import annotations
import numpy as np

def bell_singlet():
    """Returns the Bell singlet state |psi-> = (|01> - |10>)/sqrt(2)."""
    return np.array([0, 1, -1, 0]) / np.sqrt(2)

def bell_triplet_phi_plus():
    """Returns |phi+> = (|00> + |11>)/sqrt(2)."""
    return np.array([1, 0, 0, 1]) / np.sqrt(2)
