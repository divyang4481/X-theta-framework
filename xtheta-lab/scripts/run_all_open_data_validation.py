#!/usr/bin/env python3
"""
Batch runner for all available open Bell/CHSH datasets.
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

# Important: All paths are relative to xtheta-lab/
DATASETS = [
    {
        "name": "hensen",
        "dataset": "hensen",
        "data": Path("data/open_bell/hensen/raw/bell_open_data.txt"),
        "output": Path("outputs/open_data/hensen"),
    },
]

def run_dataset(item: dict) -> dict | None:
    data_path = item["data"]

    if not data_path.exists():
        print(f"[SKIP] {item['name']}: missing data file: {data_path}")
        print(f"Run: python scripts/download_open_data.py --dataset {item['dataset']}")
        return None

    cmd = [
        sys.executable,
        "scripts/run_open_data_validation.py",
        "--dataset",
        item["dataset"],
        "--data",
        str(data_path),
        "--output",
        str(item["output"]),
        "--bootstrap-samples",
        "1000",
    ]

    print(f"\n[RUN] {item['name']}")
    print(" ".join(cmd))

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Validation failed for {item['name']}: {e}")
        return None

    summary_file = item["output"] / "data" / f"{item['dataset']}_chsh_summary.csv"
    if not summary_file.exists():
        print(f"[WARN] Summary not found: {summary_file}")
        return None

    with summary_file.open("r", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
        return rows[0] if rows else None

def main() -> None:
    results = []

    for item in DATASETS:
        result = run_dataset(item)
        if result:
            results.append(result)

    comparison_dir = Path("outputs/open_data/comparison")
    comparison_dir.mkdir(parents=True, exist_ok=True)
    comparison_file = comparison_dir / "open_data_comparison.csv"

    if results:
        # Final columns as required
        fields = [
            "dataset_name",
            "row_count",
            "CHSH_S",
            "CHSH_S_se",
            "S_ci_low_95",
            "S_ci_high_95",
            "S_bootstrap_mean",
            "S_bootstrap_std",
            "phi_eff",
            "R_theta_eff",
            "fit_status",
            "interpretation_warning",
        ]

        with comparison_file.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(results)

        print(f"\n[DONE] Comparison saved to: {comparison_file}")
    else:
        print("\n[WARN] No datasets were successfully processed.")

if __name__ == "__main__":
    main()
