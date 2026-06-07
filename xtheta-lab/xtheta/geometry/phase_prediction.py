"""
X-Theta relational phase prediction from geometry.
"""
from __future__ import annotations
import numpy as np

def predict_phase_schwarzschild(rs, r_emit, r_det, theta):
    """
    Predicts the relational phase in a Schwarzschild geometry.
    Phi_pred = rs * theta * (1/r_emit - 1/r_det)

    Scientific status:
    Kinematic prediction based on Schwarzschild spacetime geometry.
    """
    if r_emit <= 0 or r_det <= 0:
        return 0.0
    return float(rs * theta * (1.0/r_emit - 1.0/r_det))
