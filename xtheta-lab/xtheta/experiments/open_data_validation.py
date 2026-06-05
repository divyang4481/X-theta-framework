"""
Generic Bell/CHSH data validation experiment runner.
"""
import pandas as pd
import numpy as np
import os
from pathlib import Path
from xtheta.data.bell_chsh import RunningAB, bootstrap_chsh
from xtheta.data.effective_fit import fit_phi_eff_from_smax, anisotropy_from_phi
from xtheta.data.schema import validate_bell_schema

SCIENTIFIC_WARNING = (
    "Phi_eff is an effective phenomenological parameter only. "
    "Without gravitational path, altitude, curvature, or spacetime-baseline metadata, "
    "this is not evidence of spacetime-induced X-Theta holonomy."
)

def run_open_data_chsh_validation(
    data_iterator,
    dataset_name: str = "unknown",
    output_dir: str = "outputs",
    bootstrap_samples: int = 1000,
    seed: int = 42
) -> dict:
    """
    Run CHSH validation on provided data and compute effective X-Theta fit.
    """
    print(f"\n--- Running Open Data CHSH Validation: {dataset_name} ---")
    print(SCIENTIFIC_WARNING)

    rab = RunningAB()
    all_events = []
    total_rows = 0
    schema_summaries = []

    for chunk in data_iterator:
        if chunk.empty:
            continue

        if total_rows < 1:
            schema_summaries.append(validate_bell_schema(chunk))

        total_rows += len(chunk)

        for a in [0, 1]:
            for b in [0, 1]:
                mask = (chunk["alice_setting"] == a) & (chunk["bob_setting"] == b)
                count = mask.sum()
                if count > 0:
                    prod = (chunk.loc[mask, "alice_outcome"] * chunk.loc[mask, "bob_outcome"]).sum()
                    idx = a * 2 + b
                    rab.count[idx] += count
                    rab.sum_ab[idx] += int(prod)

        if bootstrap_samples > 0:
            all_events.append(chunk[["alice_setting", "bob_setting", "alice_outcome", "bob_outcome"]].copy())

    if total_rows == 0:
        print("Error: No valid events processed.")
        return {}

    S = rab.chsh()
    S_se = rab.chsh_se()

    results = {
        "dataset_name": dataset_name,
        "row_count": total_rows,
        "CHSH_S": float(S),
        "CHSH_S_se": float(S_se),
        "interpretation_warning": SCIENTIFIC_WARNING
    }

    E = rab.expectation()
    for i in range(4):
        a_s, b_s = i // 2, i % 2
        results[f"count_{a_s}{b_s}"] = int(rab.count[i])
        results[f"E_{a_s}{b_s}"] = float(E[i])

    if bootstrap_samples > 0 and len(all_events) > 0:
        print(f"Calculating bootstrap (n={bootstrap_samples})...")
        full_df = pd.concat(all_events)
        boot = bootstrap_chsh(full_df["alice_outcome"].values, full_df["bob_outcome"].values,
                              full_df["alice_setting"].values, full_df["bob_setting"].values,
                              samples=bootstrap_samples, seed=seed)
        results.update(boot)
        results["bootstrap_samples"] = bootstrap_samples

    fit = fit_phi_eff_from_smax(S)
    results["phi_eff"] = fit["phi_eff"]
    results["R_theta_eff"] = fit["R_theta_eff"]
    results["fit_status"] = fit["fit_status"]
    if "warning" in fit:
        results["fit_warning"] = fit["warning"]

    os.makedirs(os.path.join(output_dir, "data"), exist_ok=True)
    summary_path = os.path.join(output_dir, "data", f"{dataset_name}_chsh_summary.csv")
    pd.DataFrame([results]).to_csv(summary_path, index=False)

    counts_df = pd.DataFrame({
        "setting_pair": ["00", "01", "10", "11"],
        "count": [results[f"count_{s}"] for s in ["00", "01", "10", "11"]],
        "expectation": [results[f"E_{s}"] for s in ["00", "01", "10", "11"]]
    })
    counts_path = os.path.join(output_dir, "data", f"{dataset_name}_setting_counts.csv")
    counts_df.to_csv(counts_path, index=False)

    os.makedirs(os.path.join(output_dir, "reports"), exist_ok=True)
    report_path = os.path.join(output_dir, "reports", f"{dataset_name}_validation_report.md")
    with open(report_path, "w") as f:
        f.write(f"# CHSH Validation Report: {dataset_name}\n\n")
        f.write(f"**Scientific Warning:** {SCIENTIFIC_WARNING}\n\n")
        f.write(f"## Summary Results\n\n")
        f.write(f"- **Total Row Count:** {total_rows}\n")
        f.write(f"- **CHSH S-statistic:** {S:.6f} ± {S_se:.6f} (Standard Error)\n")
        if "S_ci_low_95" in results:
            f.write(f"- **95% Bootstrap Confidence Interval:** [{results['S_ci_low_95']:.6f}, {results['S_ci_high_95']:.6f}]\n")
            f.write(f"- **Bootstrap Samples:** {bootstrap_samples}\n")

        f.write(f"\n## Effective X-Theta Fit\n\n")
        f.write(r"- **Effective Phase ($\Phi_{eff}$):** " + f"{results['phi_eff']:.6f} rad\n")
        f.write(r"- **Effective Anisotropy ($R_{\Theta, eff}$):** " + f"{results['R_theta_eff']:.6f}\n")
        f.write(f"- **Fit Status:** {results['fit_status']}\n")
        if "fit_warning" in results:
            f.write(f"- **Fit Warning:** {results['fit_warning']}\n")

        f.write(f"\n## Setting Expectations and Counts\n\n")
        f.write(counts_df.to_markdown(index=False))
        f.write("\n\n")

        if schema_summaries:
            f.write(f"## Schema Validation (First Chunk)\n\n")
            s = schema_summaries[0]
            f.write(f"- Missing columns: {s['missing_required_columns']}\n")
            if "unique_alice_settings" in s:
                f.write(f"- Alice unique settings: {s['unique_alice_settings']}\n")
            if "unique_bob_settings" in s:
                f.write(f"- Bob unique settings: {s['unique_bob_settings']}\n")

    print(f"Results saved to {output_dir}")
    print(f"S = {S:.4f} ± {S_se:.4f}")
    print(f"Phi_eff = {results['phi_eff']:.4f}")

    return results
