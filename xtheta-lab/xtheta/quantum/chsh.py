"""
X-Theta CHSH and Bell inequality utilities.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

def compute_chsh_variants(E00: float, E01: float, E10: float, E11: float) -> dict:
    """
    Compute all 4 sign variants of the CHSH S-statistic and find the maximum absolute value.
    Returns a dict with the variants, max_abs, and max_abs_convention.
    """
    variants = {
        "+++-": E00 + E01 + E10 - E11,
        "++-+": E00 + E01 - E10 + E11,
        "+-++": E00 - E01 + E10 + E11,
        "-+++": -E00 + E01 + E10 + E11,
    }

    max_abs_val = -1.0
    max_abs_conv = ""

    for conv, val in variants.items():
        if abs(val) > max_abs_val:
            max_abs_val = abs(val)
            max_abs_conv = conv

    result = dict(variants)
    result["max_abs"] = float(max_abs_val)
    result["max_abs_convention"] = max_abs_conv

    return result

def chsh_s_from_tensor(T: np.ndarray, a_vecs: list, b_vecs: list) -> float:
    """
    Computes the CHSH S-statistic for a given correlation tensor and measurement directions.
    S = E(a1,b1) + E(a1,b2) + E(a2,b1) - E(a2,b2)
    where E(a,b) = a^T T b.
    """
    E00 = a_vecs[0] @ T @ b_vecs[0]
    E01 = a_vecs[0] @ T @ b_vecs[1]
    E10 = a_vecs[1] @ T @ b_vecs[0]
    E11 = a_vecs[1] @ T @ b_vecs[1]

    variants = compute_chsh_variants(E00, E01, E10, E11)
    return variants["max_abs"]

def s_max_horodecki(T: np.ndarray) -> float:
    """
    Computes the Horodecki maximum CHSH S-statistic.
    S_max = 2 * sqrt(M) where M is the sum of the two largest squared singular values of T.
    """
    sing = np.linalg.svd(T, compute_uv=False)
    m = np.sort(sing**2)
    return 2.0 * np.sqrt(m[-1] + m[-2])

def bootstrap_chsh(alice_out: np.ndarray, bob_out: np.ndarray,
                   alice_set: np.ndarray, bob_set: np.ndarray,
                   samples: int = 1000, seed: int = 42) -> dict:
    """
    Perform percentile bootstrap for CHSH S-statistic.
    """
    np.random.seed(seed)
    n = len(alice_out)
    s_values = []

    masks = []
    for a in [0, 1]:
        for b in [0, 1]:
            masks.append((alice_set == a) & (bob_set == b))

    ab = alice_out * bob_out

    for _ in range(samples):
        idx = np.random.choice(n, n, replace=True)
        E = []
        for m in masks:
            m_resampled = m[idx]
            if np.any(m_resampled):
                E.append(np.mean(ab[idx][m_resampled]))
            else:
                E.append(0.0)

        if len(E) == 4:
            vars_res = compute_chsh_variants(E[0], E[1], E[2], E[3])
            s_values.append(vars_res["max_abs"])

    if not s_values:
        return {}

    s_values = np.sort(s_values)
    return {
        "S_ci_low_95": float(np.percentile(s_values, 2.5)),
        "S_ci_high_95": float(np.percentile(s_values, 97.5)),
        "S_bootstrap_mean": float(np.mean(s_values)),
        "S_bootstrap_std": float(np.std(s_values))
    }
