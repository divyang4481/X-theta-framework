import numpy as np
from xtheta.geometry.schwarzschild import compute_phi_rel
from xtheta.quantum.engine import evolve_state, compute_correlation_tensor, compute_chsh_max

def monte_carlo_sensitivity(mass, theta, r_emit, r_det,
                          mass_std=0.0, theta_std=0.0, r_emit_std=0.0, r_det_std=0.0,
                          n_samples=1000):
    """
    Runs a Monte Carlo simulation to estimate the uncertainty in Phi_rel and S_max.
    """
    phis = []
    s_maxs = []

    for _ in range(n_samples):
        m_s = np.random.normal(mass, mass_std) if mass_std > 0 else mass
        t_s = np.random.normal(theta, theta_std) if theta_std > 0 else theta
        re_s = np.random.normal(r_emit, r_emit_std) if r_emit_std > 0 else r_emit
        rd_s = np.random.normal(r_det, r_det_std) if r_det_std > 0 else r_det

        phi = compute_phi_rel(m_s, t_s, re_s, rd_s)
        psi = evolve_state(phi)
        T = compute_correlation_tensor(psi)
        s_max = compute_chsh_max(T)

        phis.append(phi)
        s_maxs.append(s_max)

    return {
        "phi_mean": np.mean(phis),
        "phi_std": np.std(phis),
        "s_max_mean": np.mean(s_maxs),
        "s_max_std": np.std(s_maxs),
        "samples": n_samples
    }
