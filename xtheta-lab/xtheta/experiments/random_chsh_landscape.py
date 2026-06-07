import numpy as np
import pandas as pd
from xtheta.quantum.engine import evolve_state, compute_correlation_tensor, compute_chsh_from_tensor, compute_chsh_max

def random_unit_vector(rng: np.random.Generator) -> np.ndarray:
    """
    Generate a random unit vector uniformly distributed on the Bloch sphere
    using the Gaussian-normalization method.
    """
    v = rng.normal(size=3)
    norm = np.linalg.norm(v)
    if norm < 1e-12:
        return random_unit_vector(rng)
    return v / norm

def random_detector_quadruple(rng: np.random.Generator):
    """
    Return a0, a1, b0, b1 as random unit vectors.
    """
    return (
        random_unit_vector(rng),
        random_unit_vector(rng),
        random_unit_vector(rng),
        random_unit_vector(rng)
    )

def simulate_random_chsh_landscape(
    Phi_values: np.ndarray = None,
    samples_per_Phi: int = 1000,
    seed: int = 42
) -> pd.DataFrame:
    """
    Simulate random CHSH landscape for a range of Phi values.
    Returns a DataFrame with columns:
    Phi, sample_id, S_random, S_abs, S_max, bell_limit, tsirelson_limit
    """
    if Phi_values is None:
        Phi_values = np.linspace(0.0, np.pi / 2, 101)

    rng = np.random.default_rng(seed)
    data = []

    bell_limit = 2.0
    tsirelson_limit = 2 * np.sqrt(2)

    for Phi in Phi_values:
        psi = evolve_state(Phi)
        T = compute_correlation_tensor(psi)
        s_max = compute_chsh_max(T)

        for i in range(samples_per_Phi):
            a0, a1, b0, b1 = random_detector_quadruple(rng)
            s_random = compute_chsh_from_tensor(T, a0, a1, b0, b1)

            data.append({
                "Phi": Phi,
                "sample_id": i,
                "S_random": s_random,
                "S_abs": abs(s_random),
                "S_max": s_max,
                "bell_limit": bell_limit,
                "tsirelson_limit": tsirelson_limit
            })

    return pd.DataFrame(data)
