from __future__ import annotations
import subprocess
import sys
import os
from pathlib import Path

# This test should be runnable from both the project root and xtheta-lab/
# We determine the path to xtheta-lab/ based on this file's location.
XTHETA_LAB_DIR = Path(__file__).parent.parent

def test_batch_runner_graceful_skip():
    # Force skip by using a non-existent path in a temporary environment
    # or just renaming the existing data file if it exists.

    data_file = XTHETA_LAB_DIR / "data/open_bell/hensen/raw/bell_open_data.txt"
    backup = XTHETA_LAB_DIR / "data/open_bell/hensen/raw/bell_open_data.txt.bak"

    moved = False
    if data_file.exists():
        data_file.rename(backup)
        moved = True

    try:
        cmd = [sys.executable, "scripts/run_all_open_data_validation.py"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=XTHETA_LAB_DIR)

        # Check that it didn't crash
        assert result.returncode == 0

        # Check that it reported a skip for hensen
        assert "[SKIP] hensen" in result.stdout
        assert "missing data file" in result.stdout
        assert "python scripts/download_open_data.py --dataset hensen" in result.stdout
    finally:
        if moved:
            backup.rename(data_file)

def test_batch_runner_comparison_csv_generation_skipped():
    data_file = XTHETA_LAB_DIR / "data/open_bell/hensen/raw/bell_open_data.txt"
    backup = XTHETA_LAB_DIR / "data/open_bell/hensen/raw/bell_open_data.txt.bak"

    moved = False
    if data_file.exists():
        data_file.rename(backup)
        moved = True

    try:
        cmd = [sys.executable, "scripts/run_all_open_data_validation.py"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=XTHETA_LAB_DIR)
        assert "No datasets were successfully processed" in result.stdout
    finally:
        if moved:
            backup.rename(data_file)
