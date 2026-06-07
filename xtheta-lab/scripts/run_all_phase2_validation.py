#!/usr/bin/env python3
"""
Master Phase 2 Research Pipeline Runner.
"""
from __future__ import annotations
import os
import sys
import subprocess
import pandas as pd
from pathlib import Path

# Ensure xtheta-lab is in path
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

def run_phase2_pipeline():
    print("=== X-Theta Phase 2 Research Pipeline ===")

    root = Path("xtheta-lab")
    reports_dir = root / "reports"
    reports_dir.mkdir(exist_ok=True)

    # 1. Hensen Audit
    print("\n[1/4] Running Hensen Open-Data Audit...")
    hensen_data = root / "data/open_bell/hensen/raw/bell_open_data.txt"
    if hensen_data.exists():
        cmd_hensen = [
            sys.executable, str(root / "scripts/run_open_data_validation.py"),
            "--dataset", "hensen",
            "--data", str(hensen_data),
            "--output", str(reports_dir / "hensen_audit")
        ]
        subprocess.run(cmd_hensen, check=True)
    else:
        print("Warning: Hensen data not found. Skipping audit.")

    # 2. Synthetic Validation
    print("\n[2/4] Running Synthetic Validation...")
    cmd_synth = [sys.executable, str(root / "scripts/run_synthetic_validation.py")]
    subprocess.run(cmd_synth, check=True)

    # 3. Summary Report Generation
    print("\n[3/4] Generating Phase 2 Summary...")
    generate_summary(reports_dir)

    print("\n[4/4] Pipeline Complete.")

def generate_summary(reports_dir):
    summary_path = reports_dir / "phase2_summary.md"

    with open(summary_path, "w") as f:
        f.write("# X-Theta Phase 2 Research Summary\n\n")
        f.write("## Milestone Status\n\n")
        f.write("| Milestone | Status | Result |\n")
        f.write("| --- | --- | --- |\n")

        # Check Hensen
        hensen_audit = reports_dir / "hensen_audit/hensen_audit.csv"
        if hensen_audit.exists():
            df = pd.read_csv(hensen_audit)
            s = df['CHSH_S'].iloc[0]
            f.write(f"| Hensen Audit | PASS | S = {s:.3f} |\n")
        else:
            f.write("| Hensen Audit | SKIPPED | - |\n")

        # Check Synthetic
        synth_csv = reports_dir / "synthetic_validation/synthetic_phi_recovery.csv"
        if synth_csv.exists():
            df = pd.read_csv(synth_csv)
            max_err = df['phi_error'].abs().max()
            status = "PASS" if max_err < 0.05 else "WARNING"
            f.write(f"| Synthetic Validation | {status} | Max Error = {max_err:.4f} |\n")
        else:
            f.write("| Synthetic Validation | SKIPPED | - |\n")

        f.write("\n## Next Steps\n")
        f.write("- Refine Pillar B theoretical derivations.\n")
        f.write("- Design metadata schema for Micius/Satellite experiments.\n")

    print(f"Summary generated at {summary_path}")

if __name__ == "__main__":
    run_phase2_pipeline()
