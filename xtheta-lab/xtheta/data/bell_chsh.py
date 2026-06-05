"""
Correlation and CHSH calculations for Bell-test data.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from xtheta.data.schema import BellEventSchema


@dataclass(init=False)
class RunningAB:
    count: np.ndarray  # shape (4,)
    sum_ab: np.ndarray 
      
    """Aggregates for Alice-Bob outcomes by (a,b) setting pairs."""
    def __init__(self, count=None, sum_ab=None):
        self.count = count if count is not None else np.zeros(4, dtype=np.int64)
        self.sum_ab = sum_ab if sum_ab is not None else np.zeros(4, dtype=np.int64)

    def update(self, a: int, b: int, ab: int, weight: int = 1):
        """
        a, b in {0, 1}
        ab in {-1, 1}
        """
        idx = a * 2 + b
        self.count[idx] += weight
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
        return float(np.sqrt(np.sum(se_terms)))

@dataclass(init=False)
class ThetaBinnedAB:
    bins: int
    count: np.ndarray  # shape (bins,4)
    sum_ab: np.ndarray  # shape (bins,4)

    def __init__(self, bins, count=None, sum_ab=None):
        self.bins = bins
        self.count = count if count is not None else np.zeros((bins, 4), dtype=np.int64)
        self.sum_ab = sum_ab if sum_ab is not None else np.zeros((bins, 4), dtype=np.int64)

    def update(self, bin_id: int, a: int, b: int, ab: int, weight: int = 1):
        idx = a * 2 + b
        self.count[bin_id, idx] += weight
        self.sum_ab[bin_id, idx] += ab

    def chsh_by_bin(self) -> Tuple[np.ndarray, np.ndarray]:
        S = np.zeros(self.bins, dtype=float)
        Nmin = np.zeros(self.bins, dtype=int)
        for k in range(self.bins):
            cnt = self.count[k]
            sab = self.sum_ab[k]
            E = np.zeros(4, dtype=float)
            nz = cnt > 0
            E[nz] = sab[nz] / cnt[nz]
            S[k] = E[0] + E[1] + E[2] - E[3]
            Nmin[k] = int(cnt.min())
        return S, Nmin

def compute_correlations_per_setting(df: pd.DataFrame) -> pd.DataFrame:
    schema = BellEventSchema()
    valid_df = df[(df[schema.alice_outcome].isin([1, -1])) & (df[schema.bob_outcome].isin([1, -1]))].copy()
    valid_df['ab'] = valid_df[schema.alice_outcome] * valid_df[schema.bob_outcome]
    grouped = valid_df.groupby([schema.alice_setting, schema.bob_setting])
    correlations = grouped['ab'].agg(['mean', 'count', 'std']).reset_index()
    correlations = correlations.rename(columns={'mean': 'E', 'std': 'E_std'})
    correlations['E_sem'] = correlations['E_std'] / np.sqrt(correlations['count'])
    return correlations

def calculate_chsh_from_correlations(correlations: pd.DataFrame, a0, a1, b0, b1) -> dict:
    def get_e(a, b):
        row = correlations[(correlations[correlations.columns[0]] == a) & (correlations[correlations.columns[1]] == b)]
        if row.empty: return 0.0, 0
        return row['E'].values[0], row['count'].values[0]

    E_a0b0, n00 = get_e(a0, b0)
    E_a0b1, n01 = get_e(a0, b1)
    E_a1b0, n10 = get_e(a1, b0)
    E_a1b1, n11 = get_e(a1, b1)
    S = E_a0b0 + E_a0b1 + E_a1b0 - E_a1b1
    return {
        "S": S,
        "counts": {"n00": n00, "n01": n01, "n10": n10, "n11": n11},
        "correlations": {"E_a0b0": E_a0b0, "E_a0b1": E_a0b1, "E_a1b0": E_a1b0, "E_a1b1": E_a1b1}
    }

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
