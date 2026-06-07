#!/usr/bin/env python3
"""
Synthetic validation script for X-Theta Phase recovery.
"""
from __future__ import annotations
import os
import sys
import pandas as pd
import numpy as np
import json
from pathlib import Path

# Ensure xtheta-lab is in path
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from xtheta.data.adapters.synthetic import generate_synthetic_xtheta_data
from xtheta.data.validation import run_open_data_chsh_validation

def run_synthetic_validation():
    phi_values = [0.0, 0.01, 0.05, 0.10, 0.20, 0.30]
    n_trials = 10000
    output_dir = Path("xtheta-lab/reports/synthetic_validation")
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = []

    print("Starting Synthetic X-Theta Validation...")

    for phi in phi_values:
        print(f"\nProcessing Phi = {phi}")
        dataset_name = f"synthetic_phi_{phi:.2f}"

        data_iterator = generate_synthetic_xtheta_data(phi, n_trials=n_trials, seed=42)

        res = run_open_data_chsh_validation(
            data_iterator,
            dataset_name=dataset_name,
            output_dir=str(output_dir / dataset_name),
            bootstrap_samples=100, # Lower for speed in synthetic
            claim_level="simulation"
        )

        res["phi_true"] = phi
        res["phi_error"] = res["phi_eff"] - phi
        all_results.append(res)

    summary_df = pd.DataFrame(all_results)
    summary_df.to_csv(output_dir / "synthetic_phi_recovery.csv", index=False)

    # Generate Markdown Report
    report_path = output_dir / "synthetic_phi_recovery.md"
    with open(report_path, "w") as f:
        f.write("# Synthetic X-Theta Phase Recovery Report\n\n")
        f.write("This report validates the ability of the X-Theta pipeline to recover a known phase from synthetic Bell-test data.\n\n")
        f.write("## Summary Table\n\n")
        f.write(summary_df[["phi_true", "phi_eff", "phi_error", "CHSH_S", "R_theta_eff"]].to_markdown(index=False))
        f.write("\n\n## Conclusion\n\n")

        max_error = summary_df["phi_error"].abs().max()
        if max_error < 0.05:
            f.write(f"SUCCESS: Maximum phase recovery error is {max_error:.4f}, within tolerance.\n")
        else:
            f.write(f"WARNING: Maximum phase recovery error is {max_error:.4f}, exceeding optimal tolerance.\n")

    print(f"\nSynthetic validation complete. Reports saved to {output_dir}")

if __name__ == "__main__":
    run_synthetic_validation()
