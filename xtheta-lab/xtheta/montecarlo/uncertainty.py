import numpy as np
from xtheta.geometry.schwarzschild import compute_Phi_rel
from xtheta.quantum.engine import (
    evolve_state,
    compute_correlation_tensor,
    compute_chsh_max,
    compute_chsh_xy,
    compute_chsh_xz,
    compute_invariants,
    compute_concurrence_from_Phi
)

def monte_carlo_sensitivity(mass, theta, r_emit, r_det,
                          mass_std=0.0, theta_std=0.0, r_emit_std=0.0, r_det_std=0.0,
                          n_samples=1000, return_samples=False):
    """
    Runs a Monte Carlo simulation to estimate the uncertainty in V2 observables.

    Protects against unphysical sampled values (mass, r_emit, r_det <= 0, theta < 0)
    by resampling.
    """
    results = {
        "Phi": [],
        "s_max": [],
        "s_xy": [],
        "s_xz": [],
        "r_theta": [],
        "concurrence": [],
        "txx": [],
        "tyy": [],
        "tzz": []
    }

    count = 0
    while count < n_samples:
        m_s = np.random.normal(mass, mass_std) if mass_std > 0 else mass
        t_s = np.random.normal(theta, theta_std) if theta_std > 0 else theta
        re_s = np.random.normal(r_emit, r_emit_std) if r_emit_std > 0 else r_emit
        rd_s = np.random.normal(r_det, r_det_std) if r_det_std > 0 else r_det

        # Protection against unphysical values
        if m_s <= 0 or re_s <= 0 or rd_s <= 0 or t_s < 0:
            continue

        Phi = compute_Phi_rel(m_s, t_s, re_s, rd_s)
        psi = evolve_state(Phi)
        T = compute_correlation_tensor(psi)

        s_max = compute_chsh_max(T)
        s_xy = compute_chsh_xy(Phi)
        s_xz = compute_chsh_xz(Phi)
        _, r_theta = compute_invariants(T)
        concurrence = compute_concurrence_from_Phi(Phi)

        results["Phi"].append(Phi)
        results["s_max"].append(s_max)
        results["s_xy"].append(s_xy)
        results["s_xz"].append(s_xz)
        results["r_theta"].append(r_theta)
        results["concurrence"].append(concurrence)
        results["txx"].append(T[0, 0])
        results["tyy"].append(T[1, 1])
        results["tzz"].append(T[2, 2])

        count += 1

    output = {
        "samples": n_samples
    }

    for key, values in results.items():
        output[f"{key}_mean"] = np.mean(values)
        output[f"{key}_std"] = np.std(values)
        if return_samples:
            output[f"{key}_samples"] = np.array(values)

    return output
