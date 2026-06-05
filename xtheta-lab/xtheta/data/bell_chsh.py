"""
Correlation and CHSH calculations for Bell-test data.
"""
import numpy as np
import math
from typing import Tuple, List

class RunningAB:
    """Aggregates for Alice-Bob outcomes by (a,b) setting pairs."""
    def __init__(self):
        self.count = np.zeros(4, dtype=np.int64)
        self.sum_ab = np.zeros(4, dtype=np.int64)

    def update(self, a: int, b: int, ab: int):
        """
        a, b in {0, 1}
        ab in {-1, 1}
        """
        idx = a * 2 + b
        self.count[idx] += 1
        self.sum_ab[idx] += ab

    def expectation(self) -> np.ndarray:
        E = np.zeros(4, dtype=float)
        nonzero = self.count > 0
        E[nonzero] = self.sum_ab[nonzero] / self.count[nonzero]
        return E

    def chsh(self) -> float:
        E = self.expectation()
        # S = E(0,0) + E(0,1) + E(1,0) - E(1,1)
        return float(E[0] + E[1] + E[2] - E[3])

    def chsh_se(self) -> float:
        """Standard error using Var(AB)=1-E^2."""
        E = self.expectation()
        se_terms = []
        for i in range(4):
            n = self.count[i]
            if n <= 0:
                se_terms.append(0.0)
                continue
            var = max(0.0, 1.0 - float(E[i] ** 2))
            se_terms.append(var / n)
        return float(math.sqrt(sum(se_terms)))

def bootstrap_chsh(alice_out: np.ndarray, bob_out: np.ndarray,
                   alice_set: np.ndarray, bob_set: np.ndarray,
                   samples: int = 1000, seed: int = 42) -> dict:
    """
    Perform percentile bootstrap for CHSH S-statistic.
    Expects arrays of same length.
    """
    np.random.seed(seed)
    n = len(alice_out)
    s_values = []

    # Pre-calculate setting masks to speed up
    masks = []
    for a in [0, 1]:
        for b in [0, 1]:
            masks.append((alice_set == a) & (bob_set == b))

    ab = alice_out * bob_out

    for _ in range(samples):
        idx = np.random.choice(n, n, replace=True)
        # In a real bootstrap for CHSH, we should resample the whole event set
        # and recompute expectations for the 4 settings.

        # Simplified for efficiency:
        E = []
        for m in masks:
            m_resampled = m[idx]
            if np.any(m_resampled):
                E.append(np.mean(ab[idx][m_resampled]))
            else:
                E.append(0.0)

        if len(E) == 4:
            s_values.append(E[0] + E[1] + E[2] - E[3])

    if not s_values:
        return {}

    s_values = np.sort(s_values)
    return {
        "S_ci_low_95": float(np.percentile(s_values, 2.5)),
        "S_ci_high_95": float(np.percentile(s_values, 97.5)),
        "S_bootstrap_mean": float(np.mean(s_values)),
        "S_bootstrap_std": float(np.std(s_values))
    }
