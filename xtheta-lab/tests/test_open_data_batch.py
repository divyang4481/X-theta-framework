from __future__ import annotations
import subprocess
import sys
import os
from pathlib import Path

def test_batch_runner_graceful_skip():
    # Force skip by using a non-existent path in a temporary environment
    # or just renaming the existing data file if it exists.

    data_file = Path("data/open_bell/hensen/raw/bell_open_data.txt")
    backup = Path("data/open_bell/hensen/raw/bell_open_data.txt.bak")

    moved = False
    if data_file.exists():
        data_file.rename(backup)
        moved = True

    try:
        cmd = [sys.executable, "scripts/run_all_open_data_validation.py"]
        result = subprocess.run(cmd, capture_output=True, text=True)

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
    data_file = Path("data/open_bell/hensen/raw/bell_open_data.txt")
    backup = Path("data/open_bell/hensen/raw/bell_open_data.txt.bak")

    moved = False
    if data_file.exists():
        data_file.rename(backup)
        moved = True

    try:
        cmd = [sys.executable, "scripts/run_all_open_data_validation.py"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert "No datasets were successfully processed" in result.stdout
    finally:
        if moved:
            backup.rename(data_file)
