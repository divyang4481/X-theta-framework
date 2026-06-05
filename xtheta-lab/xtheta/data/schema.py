from dataclasses import dataclass
import pandas as pd
import numpy as np
from collections import Counter
from typing import Dict

REQUIRED_COLUMNS = [
    "trial_id",
    "timestamp",
    "alice_setting",
    "bob_setting",
    "alice_outcome",
    "bob_outcome",
    "source_file"
]

@dataclass
class BellEventSchema:
    trial_id: str = "trial_id"
    timestamp: str = "timestamp"
    alice_setting: str = "alice_setting"
    bob_setting: str = "bob_setting"
    alice_outcome: str = "alice_outcome"
    bob_outcome: str = "bob_outcome"
    source_file: str = "source_file"

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

    def encode_array(self, arr: np.ndarray) -> np.ndarray:
        out = np.empty(len(arr), dtype=np.int16)
        for i, v in enumerate(arr):
            out[i] = self.encode(v)
        return out

    def top_k_by_frequency(self, counts: Counter, k: int) -> list:
        items = counts.most_common(k)
        return [self.mapping[val] for val, _ in items]

@dataclass
class OutcomeCoder:
    """Convert outcomes to +/- 1."""
    seen: Counter
    categorical_encoder: ValueEncoder

    def to_pm1(self, v) -> int:
        if v is True: return 1
        if v is False: return -1
        if isinstance(v, (np.generic,)): v = v.item()
        if isinstance(v, (int, np.integer)):
            if v == 1: return 1
            if v == 0 or v == -1: return -1
        if isinstance(v, (float, np.floating)):
            if abs(v - 1.0) < 1e-12: return 1
            if abs(v) < 1e-12 or abs(v + 1.0) < 1e-12: return -1

        self.seen[v] += 1
        cid = self.categorical_encoder.encode(v)
        return 1 if (cid % 2 == 1) else -1

    def array_to_pm1(self, arr: np.ndarray) -> np.ndarray:
        out = np.empty(len(arr), dtype=np.int8)
        for i, v in enumerate(arr):
            out[i] = self.to_pm1(v)
        return out

def normalize_outcomes(val):
    """Simple normalization utility."""
    try:
        fval = float(val)
        if fval > 0: return 1
        if fval <= 0: return -1
        return -1
    except:
        return -1

def validate_bell_schema(df: pd.DataFrame) -> dict:
    """
    Return validation summary for a Bell-event DataFrame:
    - missing required columns
    - row count
    - null counts
    - unique settings
    - outcome value summary
    """
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]

    summary = {
        "row_count": len(df),
        "missing_required_columns": missing,
        "null_counts": df.isnull().sum().to_dict(),
    }

    if "alice_setting" in df.columns:
        summary["unique_alice_settings"] = df["alice_setting"].unique().tolist()
    if "bob_setting" in df.columns:
        summary["unique_bob_settings"] = df["bob_setting"].unique().tolist()

    if "alice_outcome" in df.columns:
        summary["alice_outcome_counts"] = df["alice_outcome"].value_counts().to_dict()
    if "bob_outcome" in df.columns:
        summary["bob_outcome_counts"] = df["bob_outcome"].value_counts().to_dict()

    return summary
