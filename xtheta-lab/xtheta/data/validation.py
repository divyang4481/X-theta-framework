import pandas as pd
import numpy as np
import os
from xtheta.data.bell_chsh import RunningAB, bootstrap_chsh, compute_correlations_per_setting, calculate_chsh_from_correlations, compute_chsh_variants
from xtheta.data.schema import BellEventSchema, validate_bell_schema
from xtheta.data.loaders import load_bell_data
from xtheta.fitting.phi_eff_fit import fit_phi_eff

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
    seed: int = 42,
    claim_level: str = "phenomenological_fit",
    raw_row_count: int | None = None
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
                    rab.update(a, b, int(prod), weight=int(count))

        if bootstrap_samples > 0:
            all_events.append(chunk[["alice_setting", "bob_setting", "alice_outcome", "bob_outcome"]].copy())

    if total_rows == 0:
        print("Error: No valid events processed.")
        return {}

    # Validation: fail if any setting pair has near-zero count
    for i in range(4):
        if rab.count[i] < 5:  # Arbitrary threshold for "near-zero"
            a_s, b_s = i // 2, i % 2
            print(f"Error: Setting pair {a_s}{b_s} has insufficient counts ({rab.count[i]}).")
            return {"status": "failed", "error": f"Insufficient counts for setting {a_s}{b_s}"}

    S = rab.chsh()
    S_se = rab.chsh_se()

    results = {
        "dataset_name": dataset_name,
        "row_count": total_rows,
        "CHSH_S": float(S),
        "CHSH_S_se": float(S_se),
        "interpretation_warning": SCIENTIFIC_WARNING,
        "claim_level": claim_level
    }
    if raw_row_count is not None:
        results["raw_row_count"] = raw_row_count
        results["rejected_row_count"] = raw_row_count - total_rows

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

    # Fit phi_eff
    # Use max_abs variant for fitting to handle sign conventions
    variants = compute_chsh_variants(E[0], E[1], E[2], E[3])
    S_max_abs = variants["max_abs"]

    fit = fit_phi_eff(S_max_abs, 'smax-envelope')
    results["phi_eff"] = fit["phi_eff"]
    results["R_theta_eff"] = fit["R_theta_eff"]
    results["fit_status"] = fit["fit_status"]
    results["fit_warning"] = SCIENTIFIC_WARNING

    os.makedirs(output_dir, exist_ok=True)

    # Generic audit CSV
    audit_df = pd.DataFrame([results])
    audit_df.to_csv(os.path.join(output_dir, f"{dataset_name}_audit.csv"), index=False)

    counts_df = pd.DataFrame({
        "setting_pair": ["00", "01", "10", "11"],
        "count": [results[f"count_{s}"] for s in ["00", "01", "10", "11"]],
        "expectation": [results[f"E_{s}"] for s in ["00", "01", "10", "11"]]
    })

    report_path = os.path.join(output_dir, f"{dataset_name}_audit.md")
    with open(report_path, "w") as f:
        f.write(f"# X-Theta Audit Report: {dataset_name}\n\n")
        f.write(f"**Scientific Warning:** {SCIENTIFIC_WARNING}\n\n")
        f.write(f"## Claim Level\n\n")
        f.write(f"- **Level:** {claim_level}\n\n")
        f.write(f"## Summary Results\n\n")
        if raw_row_count is not None:
            f.write(f"- **Raw Row Count:** {raw_row_count}\n")
            f.write(f"- **Rejected Row Count:** {results['rejected_row_count']}\n")
        f.write(f"- **Valid Bell Trial Count:** {total_rows}\n")
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

