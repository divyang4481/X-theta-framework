"""
Fit effective X-Theta phase from observed CHSH values.
"""
import numpy as np

def fit_phi_eff_from_smax(S_observed: float) -> dict:
    """
    Uses: S_max = 2 * sqrt(1 + cos^2(2phi))
    Inverse: cos^2(2phi) = (S/2)^2 - 1
    """
    S = abs(S_observed)

    if S > 2.0 * np.sqrt(2) + 1e-9:
        return {
            "fit_status": "warning",
            "warning": "S_observed exceeds Tsirelson bound (2sqrt(2))",
            "phi_eff": 0.0,
            "R_theta_eff": 0.0
        }

    if S < 2.0:
         return {
            "fit_status": "below_bell",
            "warning": "S_observed does not violate Bell inequality (<2)",
            "phi_eff": np.pi/4,
            "R_theta_eff": 2.0
        }

    val = (S / 2.0)**2 - 1.0
    val = max(0.0, min(1.0, val)) # Clamp for safety

    cos_2phi = np.sqrt(val)
    phi_eff = 0.5 * np.arccos(cos_2phi)

    return {
        "fit_status": "success",
        "phi_eff": float(phi_eff),
        "R_theta_eff": float(anisotropy_from_phi(phi_eff))
    }

def anisotropy_from_phi(phi: float) -> float:
    """R_theta = 2sin^2(2phi)"""
    return float(2.0 * np.sin(2.0 * phi)**2)

def fit_phi_eff_from_xy(S_observed: float) -> dict:
    """Uses: S_XY = 2sqrt(2) * |cos(2phi)|"""
    S = abs(S_observed)
    limit = 2.0 * np.sqrt(2)

    if S > limit + 1e-9:
        return {
            "fit_status": "warning",
            "warning": "S_observed exceeds Tsirelson bound (2sqrt(2))",
            "phi_eff": 0.0,
            "R_theta_eff": 0.0
        }

    val = S / limit
    val = min(1.0, val)
    phi_eff = 0.5 * np.arccos(val)

    return {
        "fit_status": "success",
        "phi_eff": float(phi_eff),
        "R_theta_eff": float(anisotropy_from_phi(phi_eff))
    }

def fit_phi_eff_from_xz(S_observed: float) -> dict:
    """Uses: S_XZ = 2sqrt(2) * cos^2(phi)"""
    S = abs(S_observed)
    limit = 2.0 * np.sqrt(2)

    if S > limit + 1e-9:
        return {
            "fit_status": "warning",
            "warning": "S_observed exceeds Tsirelson bound (2sqrt(2))",
            "phi_eff": 0.0,
            "R_theta_eff": 0.0
        }

    val = S / limit
    val = min(1.0, val)
    phi_eff = np.arccos(np.sqrt(val))

    return {
        "fit_status": "success",
        "phi_eff": float(phi_eff),
        "R_theta_eff": float(anisotropy_from_phi(phi_eff))
    }
