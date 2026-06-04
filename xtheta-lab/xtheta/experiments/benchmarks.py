import numpy as np
from xtheta.geometry.schwarzschild import compute_phi_rel
from scipy.constants import g

# Constants
M_EARTH = 5.972e24
R_EARTH = 6.371e6

def get_micius_params():
    """
    Returns benchmark parameters for Micius.
    Altitude ~ 500km
    Distance ~ 1200km
    """
    altitude = 500e3
    r_emit = R_EARTH + altitude
    r_det = R_EARTH
    # theta approx distance / R_EARTH
    theta = 1200e3 / R_EARTH
    return {
        "name": "Micius",
        "mass": M_EARTH,
        "theta": theta,
        "r_emit": r_emit,
        "r_det": r_det
    }

def get_gps_params():
    """
    Returns benchmark parameters for GPS.
    Altitude ~ 20200km
    """
    altitude = 20200e3
    r_emit = R_EARTH + altitude
    r_det = R_EARTH
    theta = np.pi / 2 # Typical separation
    return {
        "name": "GPS",
        "mass": M_EARTH,
        "theta": theta,
        "r_emit": r_emit,
        "r_det": r_det
    }

def get_neutron_star_params():
    """
    Returns benchmark parameters for a Neutron Star scenario.
    M = 1.4 M_sun
    R = 10km
    """
    M_SUN = 1.989e30
    r_emit = 12e3 # Just above surface
    r_det = 1e6   # Distant observer
    theta = 1.0   # 1 radian
    return {
        "name": "Neutron Star",
        "mass": 1.4 * M_SUN,
        "theta": theta,
        "r_emit": r_emit,
        "r_det": r_det
    }

def get_black_hole_params():
    """
    Returns benchmark parameters for a Stellar Black Hole scenario.
    M = 10 M_sun
    """
    M_SUN = 1.989e30
    rs = 2 * 6.674e-11 * 10 * M_SUN / (3e8**2) # ~30km
    r_emit = 3 * rs # ISCO-ish
    r_det = 1e9    # Far away
    theta = np.pi / 2
    return {
        "name": "Black Hole",
        "mass": 10 * M_SUN,
        "theta": theta,
        "r_emit": r_emit,
        "r_det": r_det
    }

def run_benchmark_scenarios():
    scenarios = [
        get_micius_params(),
        get_gps_params(),
        get_neutron_star_params(),
        get_black_hole_params()
    ]

    results = []
    for s in scenarios:
        phi = compute_phi_rel(s['mass'], s['theta'], s['r_emit'], s['r_det'])
        results.append({
            "name": s['name'],
            "phi_rel": phi
        })
    return results
