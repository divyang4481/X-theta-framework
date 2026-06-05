"""
X-Theta Adapter for BIG Bell Test 2018 public datasets.
"""
import pandas as pd
from pathlib import Path
from typing import Iterator

def load_big_bell_test_dataset(path: str, chunksize: int = 200_000) -> Iterator[pd.DataFrame]:
    """
    Load BIG Bell Test dataset.
    Normalizes human-setting Bell tests to canonical schema.
    """
    # Placeholder for BIG Bell Test schema normalization
    # These often vary by lab (NIST repo has many labs)
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"BIG Bell Test data not found at: {path}")

    # Heuristic: assume CSV with standard names if generic
    for chunk in pd.read_csv(p, chunksize=chunksize):
        # Normalization logic goes here
        if "setting_a" in chunk.columns:
            chunk = chunk.rename(columns={"setting_a": "alice_setting", "setting_b": "bob_setting"})
        if "outcome_a" in chunk.columns:
            chunk = chunk.rename(columns={"outcome_a": "alice_outcome", "outcome_b": "bob_outcome"})

        chunk["source_file"] = p.name
        yield chunk
