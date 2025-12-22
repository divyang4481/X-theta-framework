import json
import math
import os
from dataclasses import asdict, dataclass
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


# Constants
HBAR = 1.054_571_817e-34  # J*s


@dataclass
class FitResult:
    params: Dict[str, float]
    r2: float


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


# -------------------------
# Exp 1: θ-AB fringe vs phi
# -------------------------
def ab_fringe_model(phi: np.ndarray, A: float, B: float, phi0: float, n_harm: float) -> np.ndarray:
    """
    Flexible fringe model: I(φ) = A + B * cos(n * φ + φ0)
    - For simple AB, n ≈ 1; some setups may look like cos^2 → effectively n ≈ 2 with offset.
    We'll fit n as a real to capture periodicity in the data (bounded to [0.5, 2.5]).
    """
    return A + B * np.cos(n_harm * phi + phi0)


def fit_exp1_theta_ab_fringe(csv_path: str, out_dir: str) -> FitResult:
    df = pd.read_csv(csv_path)
    phi = df["phi_theta"].to_numpy()
    I = df["intensity"].to_numpy()

    # Initial guesses
    A0 = float(np.mean(I))
    B0 = float((np.max(I) - np.min(I)) / 2)
    phi0_0 = 0.0
    n0 = 1.0

    # bounds: B >= 0, n in [0.5, 2.5] to cover cos vs cos^2-like
    bounds = ([-np.inf, 0.0, -2 * np.pi, 0.5], [np.inf, np.inf, 2 * np.pi, 2.5])
    popt, pcov = curve_fit(ab_fringe_model, phi, I, p0=[A0, B0, phi0_0, n0], bounds=bounds, maxfev=50_000)
    A, B, phi0, n_harm = [float(x) for x in popt]

    I_fit = ab_fringe_model(phi, *popt)
    ss_res = float(np.sum((I - I_fit) ** 2))
    ss_tot = float(np.sum((I - np.mean(I)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    # Plot
    ensure_dir(out_dir)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.scatter(phi, I, s=10, alpha=0.6, label="data")

    phi_dense = np.linspace(phi.min(), phi.max(), 1000)
    ax.plot(phi_dense, ab_fringe_model(phi_dense, *popt), color="C1", label=f"fit (n={n_harm:.2f})")
    ax.set_xlabel("holonomy phase φ_θ [rad]")
    ax.set_ylabel("normalized intensity")
    ax.set_title("Exp1: θ-AB fringe vs φ_θ")
    ax.legend()
    ax.grid(ls=":", alpha=0.4)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "exp1_theta_ab_fringe_fit.png"), dpi=200)
    plt.close(fig)

    return FitResult(params={"A": A, "B": B, "phi0": phi0, "n_harm": n_harm}, r2=r2)


# -----------------------------------
# Exp 2: Cross-Hall drift ~ k * T^2
# -----------------------------------
def linear_model(x: np.ndarray, m: float, b: float) -> np.ndarray:
    return m * x + b


def fit_exp2_cross_hall(csv_path: str, out_dir: str) -> FitResult:
    df = pd.read_csv(csv_path)
    # Prefer using precomputed T2 column if present
    if "T2" in df.columns:
        x = df["T2"].to_numpy(dtype=float)
    else:
        T = df["T"].to_numpy(dtype=float)
        x = T ** 2
    y = df["delta_y"].to_numpy(dtype=float)

    # Robust initial guesses
    m0 = 0.0 if np.allclose(np.var(x), 0) else (np.cov(x, y)[0, 1] / np.var(x))
    b0 = float(np.mean(y) - m0 * np.mean(x))
    popt, pcov = curve_fit(linear_model, x, y, p0=[m0, b0], maxfev=50_000)
    m, b = [float(x) for x in popt]

    y_fit = linear_model(x, *popt)
    ss_res = float(np.sum((y - y_fit) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    # Plot
    ensure_dir(out_dir)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.scatter(x, y, s=20, alpha=0.8, label="data")
    x_dense = np.linspace(x.min(), x.max(), 200)
    ax.plot(x_dense, linear_model(x_dense, *popt), color="C1", label=f"fit Δy = k T² + c")
    ax.set_xlabel("T² [s²]")
    ax.set_ylabel("Δy [m]")
    ax.set_title("Exp2: Cross-Hall drift vs T²")
    ax.legend()
    ax.grid(ls=":", alpha=0.4)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "exp2_cross_hall_T2_fit.png"), dpi=200)
    plt.close(fig)

    return FitResult(params={"kappa": m, "offset": b}, r2=r2)


# -------------------------------------------------
# Exp 3: Rotor spectrum E(ell) = ħ²/(2I) (ell - φ/2π)²
# -------------------------------------------------
def rotor_energy_model(ell: np.ndarray, I: float, phi: float) -> np.ndarray:
    # Ensure physical parameterization (I>0). curve_fit can handle this with bounds too.
    return (HBAR ** 2 / (2.0 * I)) * (ell - phi / (2.0 * np.pi)) ** 2


def fit_exp3_rotor_levels(csv_path: str, out_dir: str) -> FitResult:
    df = pd.read_csv(csv_path)
    ell = df["ell"].to_numpy(dtype=float)
    E = df["Energy_J"].to_numpy(dtype=float)

    # Initial guesses: estimate curvature to get I, and center to get phi
    # Fit a simple quadratic y = a (ell - c)^2 → a = ħ²/(2I), c = φ/2π
    # Use numpy polyfit on (ell, sqrt(E)) is tricky; do nonlinear directly with guesses.
    a0 = float(np.max(E) / max((np.max(ell) - np.min(ell)) ** 2, 1e-16))
    I0 = HBAR ** 2 / (2 * max(a0, 1e-80))
    phi0 = 0.0

    # Widen bounds and clamp initial guess into bounds to avoid ValueError
    I_min, I_max = 1e-52, 1e-38
    bounds = ([I_min, -2 * np.pi], [I_max, 2 * np.pi])
    I0 = float(min(max(I0, I_min * 1.001), I_max * 0.999))
    try:
        popt, pcov = curve_fit(rotor_energy_model, ell, E, p0=[I0, phi0], bounds=bounds, maxfev=200_000)
    except ValueError:
        # Fallback: reparameterize with K = ħ²/(2I) to improve conditioning
        def rotor_energy_model_K(l, K, phi):
            return K * (l - (phi / (2.0 * np.pi))) ** 2

        K0 = (HBAR ** 2) / (2.0 * I0)
        K_min, K_max = (HBAR ** 2) / (2.0 * I_max), (HBAR ** 2) / (2.0 * I_min)
        K0 = float(min(max(K0, K_min * 1.001), K_max * 0.999))
        popt, pcov = curve_fit(rotor_energy_model_K, ell, E, p0=[K0, phi0], bounds=([K_min, -2*np.pi], [K_max, 2*np.pi]), maxfev=200_000)
        K_fit, phi_fit = [float(x) for x in popt]
        I_fit = (HBAR ** 2) / (2.0 * K_fit)
        E_fit = rotor_energy_model_K(ell, K_fit, phi_fit)
        ss_res = float(np.sum((E - E_fit) ** 2))
        ss_tot = float(np.sum((E - np.mean(E)) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
        # Plot
        ensure_dir(out_dir)
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.scatter(ell, E, s=25, alpha=0.9, label="data")
        ell_dense = np.linspace(np.min(ell), np.max(ell), 400)
        ax.plot(ell_dense, rotor_energy_model_K(ell_dense, K_fit, phi_fit), color="C1", label="fit (K)")
        ax.set_xlabel("rotor quantum number ℓ")
        ax.set_ylabel("Energy [J]")
        ax.set_title("Exp3: Rotor spectrum vs ℓ")
        ax.legend()
        ax.grid(ls=":", alpha=0.4)
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "exp3_rotor_levels_fit.png"), dpi=200)
        plt.close(fig)
        return FitResult(params={"I": I_fit, "phi": float(phi_fit)}, r2=r2)
    I_fit, phi_fit = [float(x) for x in popt]

    E_fit = rotor_energy_model(ell, *popt)
    ss_res = float(np.sum((E - E_fit) ** 2))
    ss_tot = float(np.sum((E - np.mean(E)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    # Plot
    ensure_dir(out_dir)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.scatter(ell, E, s=25, alpha=0.9, label="data")
    ell_dense = np.linspace(np.min(ell), np.max(ell), 400)
    ax.plot(ell_dense, rotor_energy_model(ell_dense, *popt), color="C1", label="fit")
    ax.set_xlabel("rotor quantum number ℓ")
    ax.set_ylabel("Energy [J]")
    ax.set_title("Exp3: Rotor spectrum vs ℓ")
    ax.legend()
    ax.grid(ls=":", alpha=0.4)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "exp3_rotor_levels_fit.png"), dpi=200)
    plt.close(fig)

    return FitResult(params={"I": I_fit, "phi": phi_fit}, r2=r2)


def main():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    paper_dir = os.path.join(repo_root, "paper")
    data_dir = os.path.join(paper_dir, "data")
    build_dir = os.path.join(paper_dir, "build")
    figs_out = os.path.join(build_dir, "figs", "regression")
    analysis_out = os.path.join(build_dir, "analysis")
    ensure_dir(figs_out)
    ensure_dir(analysis_out)

    results: Dict[str, Dict[str, float]] = {}

    # Exp 1
    exp1_csv = os.path.join(data_dir, "exp1_theta_ab_fringe.csv")
    if os.path.exists(exp1_csv):
        fit1 = fit_exp1_theta_ab_fringe(exp1_csv, figs_out)
        results["exp1_theta_ab_fringe"] = {**fit1.params, "r2": fit1.r2}
    else:
        print(f"Warning: missing {exp1_csv}")

    # Exp 2
    exp2_csv = os.path.join(data_dir, "exp2_drift_T2.csv")
    if os.path.exists(exp2_csv):
        fit2 = fit_exp2_cross_hall(exp2_csv, figs_out)
        results["exp2_cross_hall_T2"] = {**fit2.params, "r2": fit2.r2}
    else:
        print(f"Warning: missing {exp2_csv}")

    # Exp 3
    exp3_csv = os.path.join(data_dir, "exp3_rotor_levels.csv")
    if os.path.exists(exp3_csv):
        fit3 = fit_exp3_rotor_levels(exp3_csv, figs_out)
        results["exp3_rotor_levels"] = {**fit3.params, "r2": fit3.r2}
    else:
        print(f"Warning: missing {exp3_csv}")

    # Write JSON summary
    summary_path = os.path.join(analysis_out, "regression_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"Wrote regression summary → {summary_path}")
    print(f"Figures → {figs_out}")


if __name__ == "__main__":
    main()
