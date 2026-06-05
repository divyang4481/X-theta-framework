"""
Logic to convert outcome values to ±1.
"""
import numpy as np
from collections import Counter
from dataclasses import dataclass

@dataclass
class ValueEncoder:
    """Stream-safe factor encoder: maps arbitrary values to small ints as they appear."""
    mapping: dict
    next_id: int = 0

    def encode(self, v) -> int:
        if v in self.mapping:
            return self.mapping[v]
        self.mapping[v] = self.next_id
        self.next_id += 1
        return self.mapping[v]

class OutcomeCoder:
    """
    Convert outcomes to ±1.
    Handles:
      - already ±1
      - 0/1 -> -1/+1
      - True/False -> -1/+1
      - arbitrary categorical -> factorize then map first->-1 second->+1 (only sensible for 2-category)
    """
    def __init__(self):
        self.seen = Counter()
        self.categorical_encoder = ValueEncoder(mapping={})

    def to_pm1(self, v) -> int:
        if v is True:
            return 1
        if v is False:
            return -1

        if isinstance(v, (np.generic,)):
            v = v.item()

        if isinstance(v, (int, np.integer)):
            if v == 1:
                return 1
            if v == 0:
                return -1
            if v == -1:
                return -1

        if isinstance(v, (float, np.floating)):
            if abs(v - 1.0) < 1e-12:
                return 1
            if abs(v) < 1e-12:
                return -1
            if abs(v + 1.0) < 1e-12:
                return -1

        # categorical fallback
        self.seen[v] += 1
        cid = self.categorical_encoder.encode(v)
        return 1 if (cid % 2 == 1) else -1

    def array_to_pm1(self, arr: np.ndarray) -> np.ndarray:
        out = np.empty(len(arr), dtype=np.int8)
        for i, v in enumerate(arr):
            out[i] = self.to_pm1(v)
        return out
