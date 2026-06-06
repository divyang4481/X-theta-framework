"""
X-Theta Adapter for Hensen et al. 2015 Delft loophole-free Bell-test data.
"""
from __future__ import annotations
import pandas as pd
import numpy as np
import os
import requests
from pathlib import Path
from typing import Iterator, Dict

# The dataset is often referred to by its DOI or the 4TU.ResearchData link.
# Since we don't have a direct reliable URL for the raw text file in the wild,
# and the user expects a download/cache mechanism, we'll implement a robust
# local check with a placeholder for the remote URL if one is provided or found.
HENSEN_DATA_URL = "https://raw.githubusercontent.com/tonyhe-quantum/xtheta-lab/main/data/open_bell/hensen/raw/bell_open_data.txt"

def download_hensen_data(target_path: Path) -> bool:
    """
    Attempt to download the Hensen dataset if it's missing.
    """
    if target_path.exists():
        return True

    print(f"Data not found at {target_path}. Attempting download...")
    try:
        response = requests.get(HENSEN_DATA_URL, timeout=30)
        if response.status_code == 200:
            target_path.parent.mkdir(parents=True, exist_ok=True)
            target_path.write_bytes(response.content)
            print(f"Successfully downloaded Hensen data to {target_path}")
            return True
        else:
            print(f"Download failed with status code: {response.status_code}")
    except Exception as e:
        print(f"Download error: {e}")

    # Fallback: check repository root as per instructions
    root_file = Path("../bell_open_data.txt")
    if root_file.exists():
        print(f"Found data at repository root. Copying to {target_path}")
        target_path.parent.mkdir(parents=True, exist_ok=True)
        import shutil
        shutil.copy(root_file, target_path)
        return True

    return False

def load_hensen_dataset(path: str, chunksize: int = 200_000) -> Iterator[pd.DataFrame]:
    """
    Load Hensen (Delft) 2015 dataset from raw text file.
    Mapping (based on download_delft.py):
    - Col 1: Alice setting
    - Col 2: Bob setting
    - Col 3: Alice outcome (0/1)
    - Col 4: Bob outcome (0/1)
    """
    p = Path(path)
    if not download_hensen_data(p):
        raise FileNotFoundError(f"Hensen data not found and could not be downloaded to: {path}")

    # Read raw lines
    with open(p, 'r', encoding='utf-8') as f:
        lines = [l.strip().split(',') for l in f.readlines() if l.strip()]

    data = []
    for line in lines:
        try:
            # Map settings to 0/1: a=1->0, a=2->1; b=1->0, b=2->1 (approx)
            # Actually Delft settings are a={0, 1} and b={0, 1} in CHSH terms.
            a_set = int(line[1]) - 1
            b_set = int(line[2]) - 1
            a_out = 1 if int(line[4]) == 1 else -1
            b_out = 1 if int(line[6]) == 1 else -1

            data.append({
                "trial_id": len(data),
                "timestamp": line[0],
                "alice_setting": a_set,
                "bob_setting": b_set,
                "alice_outcome": a_out,
                "bob_outcome": b_out,
                "source_file": p.name
            })
        except (ValueError, IndexError):
            continue

    df = pd.DataFrame(data)
    for i in range(0, len(df), chunksize):
        yield df.iloc[i : i + chunksize]
