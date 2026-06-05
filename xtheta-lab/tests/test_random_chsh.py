import numpy as np
import pytest
from xtheta.experiments.random_chsh_landscape import random_unit_vector, simulate_random_chsh_landscape

def test_random_unit_vector():
    rng = np.random.default_rng(42)
    for _ in range(100):
        v = random_unit_vector(rng)
        assert abs(np.linalg.norm(v) - 1.0) < 1e-12

def test_random_chsh_never_exceeds_horodecki_envelope():
    # Use small samples for speed in test
    phi_values = np.linspace(0.0, np.pi/2, 11)
    df = simulate_random_chsh_landscape(phi_values=phi_values, samples_per_phi=50, seed=42)

    # abs(S_random) <= S_max + 1e-9
    # S_max <= 2sqrt(2) + 1e-9
    assert all(df['S_abs'] <= df['S_max'] + 1e-9)
    assert all(df['S_max'] <= 2 * np.sqrt(2) + 1e-9)
