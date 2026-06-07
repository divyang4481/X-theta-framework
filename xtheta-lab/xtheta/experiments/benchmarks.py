import numpy as np
from xtheta.geometry.schwarzschild import compute_Phi_rel
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

def compute_scenario_results(name, mass, theta, r_emit, r_det):
    from xtheta.quantum.engine import (
        evolve_state, compute_correlation_tensor, compute_chsh_max,
        compute_invariants, compute_concurrence_from_Phi, compute_purity,
        compute_chsh_xy, compute_chsh_xz
    )

    Phi = compute_Phi_rel(mass, theta, r_emit, r_det)
    psi = evolve_state(Phi)
    T = compute_correlation_tensor(psi)
    i_theta, r_theta = compute_invariants(T)
    s_max = compute_chsh_max(T)
    conc = compute_concurrence_from_Phi(Phi)
    purity = compute_purity(psi)
    s_xy = compute_chsh_xy(Phi)
    s_xz = compute_chsh_xz(Phi)

    return {
        "scenario": name,
        "mass_kg": mass,
        "theta_rad": theta,
        "r_emit_m": r_emit,
        "r_det_m": r_det,
        "Phi_rel": Phi,
        "abs_Phi_rel": abs(Phi),
        "Txx": T[0,0],
        "Tyy": T[1,1],
        "Tzz": T[2,2],
        "I_theta": i_theta,
        "R_theta": r_theta,
        "concurrence": conc,
        "purity": purity,
        "S_xz": s_xz,
        "S_xy": s_xy,
        "S_max": s_max,
        "delta_S_from_tsirelson": abs(s_max - 2*np.sqrt(2)),
        # Scientific notation fields
        "Phi_rel_scientific": f"{Phi:.4e}",
        "R_theta_scientific": f"{r_theta:.4e}",
        "delta_S_scientific": f"{abs(s_max - 2*np.sqrt(2)):.4e}"
    }

def run_benchmark_scenarios(custom_scenarios=None):
    scenarios = [
        get_micius_params(),
        get_gps_params(),
        get_neutron_star_params(),
        get_black_hole_params()
    ]
    if custom_scenarios:
        scenarios.extend(custom_scenarios)

    results = []
    for s in scenarios:
        results.append(compute_scenario_results(s['name'], s['mass'], s['theta'], s['r_emit'], s['r_det']))
    return results
