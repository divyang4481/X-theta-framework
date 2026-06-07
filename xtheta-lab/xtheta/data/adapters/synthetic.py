"""
Synthetic X-Theta Bell data generator.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from typing import Iterator
from xtheta.quantum.correlation_tensor import get_correlation_tensor

def generate_synthetic_xtheta_data(
    Phi: float,
    n_trials: int = 1000,
    geometry: str = 'smax-envelope',
    seed: int = 42,
    chunksize: int = 1000
) -> Iterator[pd.DataFrame]:
    """
    Generate synthetic Bell-test data following X-Theta anisotropy.
    """
    rng = np.random.default_rng(seed)

    # Correlation tensor T(Phi) = diag(-cos(2Phi), -cos(2Phi), -1)
    T = get_correlation_tensor(Phi)

    # Standard CHSH settings for maximal violation (for Phi=0)
    # A1 = Z, A2 = X
    # B1 = (Z+X)/sqrt(2), B2 = (Z-X)/sqrt(2)

    A_vecs = [
        np.array([0, 0, 1]), # a1: Z
        np.array([1, 0, 0])  # a2: X
    ]
    B_vecs = [
        np.array([1, 0, 1]) / np.sqrt(2),  # b1
        np.array([-1, 0, 1]) / np.sqrt(2)  # b2
    ]

    for i in range(0, n_trials, chunksize):
        this_chunk = min(chunksize, n_trials - i)

        alice_settings = rng.choice([0, 1], size=this_chunk)
        bob_settings = rng.choice([0, 1], size=this_chunk)

        alice_outcomes = []
        bob_outcomes = []

        for s_a, s_b in zip(alice_settings, bob_settings):
            A = A_vecs[s_a]
            B = B_vecs[s_b]

            # Correlation E = A^T T B
            # Flip sign to match +++- convention for singlet state
            E = -(A @ T @ B)

            # Joint probabilities P(a,b) = 1/4 (1 + a*b*E)
            p_pp = 0.25 * (1 + E)
            p_mm = 0.25 * (1 + E)
            p_pm = 0.25 * (1 - E)
            p_mp = 0.25 * (1 - E)

            probs = np.clip([p_pp, p_pm, p_mp, p_mm], 0, 1)
            probs /= probs.sum()

            choice = rng.choice([0, 1, 2, 3], p=probs)
            a_out = 1 if choice in [0, 1] else -1
            b_out = 1 if choice in [0, 2] else -1

            alice_outcomes.append(a_out)
            bob_outcomes.append(b_out)

        yield pd.DataFrame({
            "alice_setting": alice_settings,
            "bob_setting": bob_settings,
            "alice_outcome": alice_outcomes,
            "bob_outcome": bob_outcomes,
            "Phi_true": Phi
        })
