"""
X-Theta Adapter for Weihs 1998 photon Bell-test data.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Iterator
from xtheta.data.loaders import iter_npz

def load_weihs_dataset(path: str, chunksize: int = 200_000) -> Iterator[pd.DataFrame]:
    """
    Load Weihs dataset and normalize to canonical Bell-event schema.
    Supports .npz files from the Weihs data directory.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Weihs data not found at: {path}")

    for chunk in iter_npz(p, chunksize):
        df = pd.DataFrame()
        if "settings" in chunk.columns and "products" in chunk.columns:
            s = chunk["settings"].values.astype(int)
            df["alice_setting"] = s // 2
            df["bob_setting"] = s % 2
            df["alice_outcome"] = chunk["products"]
            df["bob_outcome"] = 1.0 # products is already A*B
            df["trial_id"] = np.arange(len(chunk))
            df["timestamp"] = 0.0
            df["source_file"] = str(p.name)
        elif "prods" in chunk.columns:
            df["alice_outcome"] = chunk["prods"]
            df["bob_outcome"] = 1.0
            # Heuristic for settings in v20+ files if available
            # If 'phases' or other keys exist, they might need mapping
            df["alice_setting"] = 0
            df["bob_setting"] = 0
            df["trial_id"] = np.arange(len(chunk))
            df["timestamp"] = 0.0
            df["source_file"] = str(p.name)

        # Drop rows with NaN settings or outcomes
        df = df.dropna(subset=["alice_setting", "bob_setting", "alice_outcome", "bob_outcome"])
        yield df
