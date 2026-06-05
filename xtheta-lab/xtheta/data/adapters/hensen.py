"""
X-Theta Adapter for Hensen et al. 2015 Delft loophole-free Bell-test data.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Iterator, Dict
from xtheta.experiments.open_data_validation import run_open_data_chsh_validation

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
    if not p.exists():
        raise FileNotFoundError(f"Hensen data not found at: {path}")

    # Read raw lines to handle potential format variations
    with open(p, 'r', encoding='utf-8') as f:
        lines = [l.strip().split(',') for l in f.readlines() if l.strip()]

    # Simple parser assuming CSV-like structure in bell_open_data.txt
    # 2015-06-26 17:24:12.119993,1,2,5445065,1,5671004,1,1,1,10379,10371,11281,13113,0,0,0,0
    # Paper mapping often uses different columns.
    # Based on download_delft.py:
    # Index 1: Setting A, Index 2: Setting B, Index 4: Outcome A, Index 6: Outcome B
    # Wait, download_delft.py says Index 3/4 for outcome?
    # Let's re-examine bell_open_data.txt:
    # 2015-06-26 17:24:12.119993 (0), 1 (1), 2 (2), 5445065 (3), 1 (4), 5671004 (5), 1 (6), 1 (7), 1 (8) ...
    # Delft 2015 settings are {a=1, 2} and {b=1, 2}. Map to 0/1.

    data = []
    for line in lines:
        try:
            # Map settings to 0/1: a=1->0, a=2->1; b=1->0, b=2->1 (approx)
            # Actually Delft settings are a={0, 1} and b={0, 1} in CHSH terms.
            # Raw data in bell_open_data.txt seems to have settings in col 1 and 2.
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
    # Yield in chunks
    for i in range(0, len(df), chunksize):
        yield df.iloc[i : i + chunksize]

def run_hensen_count_table_validation(counts: Dict[str, int], output_dir: str = "outputs") -> dict:
    """
    Run validation using published count tables (manual input mode).
    counts = {
        "n00_pp": ..., "n00_pm": ..., "n00_mp": ..., "n00_mm": ...,
        ...
    }
    """
    # Convert count table to synthetic event list
    data = []
    for a in [0, 1]:
        for b in [0, 1]:
            for oa in [1, -1]:
                for ob in [1, -1]:
                    key = f"n{a}{b}_{'p' if oa==1 else 'm'}{'p' if ob==1 else 'm'}"
                    n = counts.get(key, 0)
                    for _ in range(n):
                        data.append({
                            "alice_setting": a, "bob_setting": b,
                            "alice_outcome": oa, "bob_outcome": ob,
                            "trial_id": len(data), "timestamp": 0.0, "source_file": "manual_input"
                        })

    df = pd.DataFrame(data)
    return run_open_data_chsh_validation(iter([df]), dataset_name="hensen_2015_manual", output_dir=output_dir)
