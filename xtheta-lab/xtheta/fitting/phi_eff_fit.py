"""
X-Theta effective phase fitting.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize_scalar

def fit_phi_eff(S_observed: float, geometry: str = 'smax-envelope') -> dict:
    """
    Fits an effective phi value from an observed CHSH S value.
    """
    def objective(phi):
        if geometry == 'xy':
            S_theory = 2 * np.sqrt(2) * abs(np.cos(2 * phi))
        elif geometry == 'xz':
            S_theory = 2 * np.sqrt(2) * (np.cos(phi)**2)
        elif geometry == 'smax-envelope':
            S_theory = 2 * np.sqrt(1 + np.cos(2 * phi)**2)
        else:
            raise ValueError(f"Unknown geometry: {geometry}")
        return (abs(S_observed) - S_theory)**2

    res = minimize_scalar(objective, bounds=(0, np.pi/4), method='bounded')
    phi_eff = float(res.x)
    r_theta_eff = float(2 * (np.sin(2 * phi_eff)**2))

    return {
        "phi_eff": phi_eff,
        "R_theta_eff": r_theta_eff,
        "fit_status": "Success" if res.success else "Failed"
    }
