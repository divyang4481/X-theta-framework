import numpy as np
from scipy.constants import G, c

def get_schwarzschild_radius(mass):
    """Returns r_s = 2GM/c^2."""
    return 2 * G * mass / c**2

def compute_phi_rel(mass, theta, r_emit, r_det):
    """
    Computes the relational phase:
    Phi_rel = r_s * theta * (1/r_emit - 1/r_det)

    mass: Mass of the central body (kg)
    theta: Angular separation (radians)
    r_emit: Emission radius (meters)
    r_det: Detection radius (meters)
    """
    rs = get_schwarzschild_radius(mass)
    phi = rs * theta * (1.0/r_emit - 1.0/r_det)
    return phi

def compute_phi_rel_simplified(rs, theta, r_emit, r_det):
    """Computes Phi_rel given r_s directly."""
    return rs * theta * (1.0/r_emit - 1.0/r_det)
