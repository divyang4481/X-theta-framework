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
    Implements the official filtering and mapping logic from the 2015 Nature paper.
    """
    p = Path(path)
    if not download_hensen_data(p):
        raise FileNotFoundError(f"Hensen data not found and could not be downloaded to: {path}")

    # Read raw data
    # The file has no header. Col 0 is timestamp, followed by 16 data columns.
    df_raw = pd.read_csv(p, header=None)

    # Official constants for event-ready (heralding) and readout windows
    EVENT_READY_WINDOW_START_CH0 = 5426350
    EVENT_READY_WINDOW_START_CH1 = 5425700
    EVENT_READY_WINDOW_LENGTH = 52450
    EVENT_READY_WINDOW_SEPARATION = 250000
    READOUT_WINDOW_START = 10620
    READOUT_WINDOW_LENGTH = 3700
    CHECK_FOR_INVALID_MARKER_IN_PAST = 250

    # Column mapping (0-indexed based on raw file)
    ER_CLICK1_TIME = 3
    ER_CLICK1_CH = 4
    ER_CLICK2_TIME = 5
    ER_CLICK2_CH = 6
    RN_A = 7
    RN_B = 8
    RO_CLICK_A_TIME = 11
    RO_CLICK_B_TIME = 12
    CLICK_AFTER_EXCITE_A = 13
    CLICK_AFTER_EXCITE_B = 14
    INVALID_MARKER_A = 15
    INVALID_MARKER_B = 16

    # 1. Heralding Filters (Event Ready)
    t1 = df_raw[ER_CLICK1_TIME]
    ch1 = df_raw[ER_CLICK1_CH]
    t2 = df_raw[ER_CLICK2_TIME]
    ch2 = df_raw[ER_CLICK2_CH]

    filter_w1_ch0 = (EVENT_READY_WINDOW_START_CH0 <= t1) & (t1 < EVENT_READY_WINDOW_START_CH0 + EVENT_READY_WINDOW_LENGTH) & (ch1 == 0)
    filter_w1_ch1 = (EVENT_READY_WINDOW_START_CH1 <= t1) & (t1 < EVENT_READY_WINDOW_START_CH1 + EVENT_READY_WINDOW_LENGTH) & (ch1 == 1)
    w1_filter = filter_w1_ch0 | filter_w1_ch1

    filter_w2_ch0 = (EVENT_READY_WINDOW_START_CH0 + EVENT_READY_WINDOW_SEPARATION <= t2) & (t2 < EVENT_READY_WINDOW_START_CH0 + EVENT_READY_WINDOW_SEPARATION + EVENT_READY_WINDOW_LENGTH) & (ch2 == 0)
    filter_w2_ch1 = (EVENT_READY_WINDOW_START_CH1 + EVENT_READY_WINDOW_SEPARATION <= t2) & (t2 < EVENT_READY_WINDOW_START_CH1 + EVENT_READY_WINDOW_SEPARATION + EVENT_READY_WINDOW_LENGTH) & (ch2 == 1)
    w2_filter = filter_w2_ch0 | filter_w2_ch1

    psi_min_filter = (ch1 != ch2)
    ready_filter = w1_filter & w2_filter & psi_min_filter

    # 2. Signal Integrity Filters
    inv_a = df_raw[INVALID_MARKER_A]
    inv_b = df_raw[INVALID_MARKER_B]
    no_invalid_marker = ((inv_a == 0) | (inv_a > CHECK_FOR_INVALID_MARKER_IN_PAST)) & \
                        ((inv_b == 0) | (inv_b > CHECK_FOR_INVALID_MARKER_IN_PAST))

    exc_a = df_raw[CLICK_AFTER_EXCITE_A]
    exc_b = df_raw[CLICK_AFTER_EXCITE_B]
    no_excitation = (exc_a == 0) & (exc_b == 0)

    # Final Bell Trial Filter
    bell_trial_filter = ready_filter & no_invalid_marker & no_excitation
    df_filtered = df_raw[bell_trial_filter].copy()

    if df_filtered.empty:
        return iter([])

    # 3. Derive Outcomes from Readout Windows
    ro_a = df_filtered[RO_CLICK_A_TIME]
    ro_b = df_filtered[RO_CLICK_B_TIME]
    det_a = (ro_a > READOUT_WINDOW_START) & (ro_a <= READOUT_WINDOW_START + READOUT_WINDOW_LENGTH)
    det_b = (ro_b > READOUT_WINDOW_START) & (ro_b <= READOUT_WINDOW_START + READOUT_WINDOW_LENGTH)

    # 4. Map to Canonical Schema
    df_final = pd.DataFrame({
        "trial_id": df_filtered.index,
        "timestamp": df_filtered[0],
        "alice_setting": df_filtered[RN_A].astype(int),
        "bob_setting": df_filtered[RN_B].astype(int),
        "alice_outcome": np.where(det_a, 1, -1),
        "bob_outcome": np.where(det_b, 1, -1),
        "source_file": p.name
    })

    for i in range(0, len(df_final), chunksize):
        yield df_final.iloc[i : i + chunksize]
