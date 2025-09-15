# Add fitting & figure-factory to the X–θ simulation toolkit
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, Tuple
import json
import pandas as pd
from caas_jupyter_tools import display_dataframe_to_user

# Reuse/define physical constants
hbar = 1.054_571_817e-34

# -----------------------------
# Utilities
# -----------------------------

def linear_fit(x, y):
    X = np.vstack([np.ones_like(x), x]).T
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    resid = y - yhat
    dof = max(1, len(y) - 2)
    sigma2 = (resid @ resid) / dof
    cov = sigma2 * np.linalg.inv(X.T @ X)
    return {"intercept": beta[0], "slope": beta[1], "sigma2": sigma2, "cov": cov, "yhat": yhat}

def cosine_fit(phi, I):
    # Fit I = a + b cos(phi) + c sin(phi)
    X = np.vstack([np.ones_like(phi), np.cos(phi), np.sin(phi)]).T
    beta, *_ = np.linalg.lstsq(X, I, rcond=None)
    a, b, c = beta
    I0 = a
    V = np.sqrt(b**2 + c**2) / max(1e-12, np.abs(I0))
    phi0 = np.arctan2(-c, b)  # so that I ≈ I0 [1 + V cos(phi - phi0)]
    I_fit = I0 * (1.0 + V * np.cos(phi - phi0))
    resid = I - I_fit
    dof = max(1, len(I) - 3)
    sigma2 = (resid @ resid) / dof
    return {"I0": I0, "V": V, "phi0": phi0, "sigma2": sigma2, "I_fit": I_fit}

def double_integral(dt, signal):
    # Return S(t) = ∫0^t dt1 ∫0^t1 signal(t2) dt2 (two cumulative integrals)
    v = np.cumsum(signal) * dt
    x = np.cumsum(v) * dt
    return x

# -----------------------------
# 1) θ–AB: phase vs Φθ slope from intensity-only data
# -----------------------------

def fit_theta_ab_from_intensity(Phi_theta, I_signal):
    # Bin or smooth if needed; here we directly extract phase at each setting via cosine fit on small windows.
    # For synthetic scan with monotonically increasing Φθ, we can do a sliding window to estimate local phase.
    # Simpler: treat each point as I(Φθ) ≈ I0 [1 + V cos(α Φθ + φ0)] and fit α,φ0,I0,V globally.
    # Build design for non-linear α: do grid search for α, φ0; solve linear for I0,V given cos term.
    Phi = Phi_theta
    alphas = np.linspace(0.8, 1.2, 161)  # ±20% search around 1
    phi0s = np.linspace(-np.pi, np.pi, 181)
    best = None
    for a in alphas:
        cos_term = np.cos(a*Phi)
        sin_term = np.sin(a*Phi)
        X = np.vstack([np.ones_like(Phi), cos_term, sin_term]).T
        beta, *_ = np.linalg.lstsq(X, I_signal, rcond=None)
        I_fit = X @ beta
        resid = I_signal - I_fit
        rss = resid @ resid
        # Recover V and φ0 from b,c
        I0, b, c = beta
        V = np.sqrt(b**2 + c**2)/max(1e-12, np.abs(I0))
        phi0 = np.arctan2(-c, b)
        cand = {"alpha": a, "phi0": phi0, "I0": I0, "V": V, "rss": rss, "I_fit": I_fit}
        if best is None or rss < best["rss"]:
            best = cand
    return best

# -----------------------------
# 2) Cross-Hall drift: estimate F0/m from +Ω and −Ω runs
# -----------------------------

def fit_cross_hall(t, x_plus, x_minus, theta_dot_plus, theta_dot_minus, dt):
    # Model: x(t) ≈ a * S_theta(t) + b + c t, where S_theta is double integral of θdot
    S_plus  = double_integral(dt, theta_dot_plus)
    S_minus = double_integral(dt, theta_dot_minus)
    # Build antisymmetric (odd) combo to cancel b,c: x_odd = (x+ - x-)/2 ≈ a * (S+ - S-)/2
    x_odd = 0.5*(x_plus - x_minus)
    S_odd = 0.5*(S_plus - S_minus)
    X = S_odd.reshape(-1,1)
    # Fit a single-parameter slope (force/mass scale)
    a_hat = np.linalg.lstsq(X, x_odd, rcond=None)[0][0]
    # Uncertainty estimate
    yhat = X[:,0]*a_hat
    resid = x_odd - yhat
    dof = max(1, len(x_odd)-1)
    sigma2 = (resid @ resid) / dof
    var_a = sigma2 / (X[:,0]@X[:,0])
    return {"a_hat": a_hat, "sigma_a": np.sqrt(var_a), "t": t, "x_odd": x_odd, "S_odd": S_odd, "yhat": yhat}

# -----------------------------
# 3) Sidebands: fit I (inertia) and Φθ shift from peaks
# -----------------------------

def fit_sidebands(n, E_meas, w=1.0):
    # Grid over φ in [0,1); linear least-squares for c = ħ^2/(2I) at each φ
    phis = np.linspace(0.0, 1.0, 1001, endpoint=False)
    best = None
    W = np.eye(len(n)) * w if np.isscalar(w) else np.diag(w)
    for phi in phis:
        f = (n - phi)**2
        # Solve E ≈ c * f  (weighted)
        A = f.reshape(-1,1)
        Aw = W @ A
        Ew = W @ E_meas
        c_hat = np.linalg.lstsq(Aw, Ew, rcond=None)[0][0]
        E_fit = c_hat * f
        resid = E_meas - E_fit
        rss = float(resid @ (W @ resid))
        if (best is None) or (rss < best["rss"]):
            best = {"phi": phi, "c": c_hat, "rss": rss, "E_fit": E_fit}
    # Convert c to inertia I
    I_hat = hbar**2/(2.0*best["c"]) if best["c"]>0 else np.inf
    return {"phi_hat": best["phi"], "I_hat": I_hat, "c_hat": best["c"], "rss": best["rss"], "E_fit": best["E_fit"]}

# -----------------------------
# Generate synthetic datasets (reusing earlier simulator logic)
# -----------------------------

# θ–AB synthetic scan
Phi = np.linspace(-2*np.pi, 2*np.pi, 800)
I0_true, V_true, alpha_true, phi0_true, sigma_phase = 1.0, 0.94, 1.03, 0.2, 0.01
# Build intensity directly from a cosine model (lab observable)
I_signal = I0_true * (1.0 + V_true * np.cos(alpha_true*Phi - phi0_true)) + np.random.normal(0, 0.005, size=Phi.size)

# Cross-Hall synthetic (+Ω and −Ω)
T, dt = 0.03, 2e-5
t = np.arange(0, T, dt)
theta_dot_plus  = np.full_like(t, 1500.0)    # rad/s
theta_dot_minus = np.full_like(t, -1500.0)   # rad/s
a_true = 2.5e-21   # dimensionless scale a = F0/m  (in our toy units)
# x(t) = a * S_theta(t) + b + c t + noise
S_plus  = double_integral(dt, theta_dot_plus)
S_minus = double_integral(dt, theta_dot_minus)
b_true, c_true = 2e-9, -5e-9
x_plus  = a_true*S_plus  + b_true + c_true*t + np.random.normal(0, 2e-10, size=t.size)
x_minus = a_true*S_minus + b_true + c_true*t + np.random.normal(0, 2e-10, size=t.size)

# Sidebands synthetic
n = np.arange(-5, 6)
I_inertia_true = 1.2e-46
phi_true = 0.28  # Φθ/(2π)
c_true = hbar**2/(2*I_inertia_true)
E_true = c_true * (n - phi_true)**2
E_meas = E_true + np.random.normal(0, 0.03*np.max(E_true), size=E_true.size)

# -----------------------------
# Fits
# -----------------------------

theta_ab_fit = fit_theta_ab_from_intensity(Phi, I_signal)
cross_hall_fit = fit_cross_hall(t, x_plus, x_minus, theta_dot_plus, theta_dot_minus, dt)
sideband_fit = fit_sidebands(n, E_meas)

# -----------------------------
# Plots: figure factory templates
# -----------------------------

# Fig A: θ–AB fit
plt.figure()
plt.plot(Phi, I_signal, ".", ms=2, label="data")
a = theta_ab_fit["alpha"]; phi0 = theta_ab_fit["phi0"]; I0 = theta_ab_fit["I0"]; V = theta_ab_fit["V"]
I_fit = I0 * (1.0 + V * np.cos(a*Phi - phi0))
plt.plot(Phi, I_fit, "-", label=f"fit: α={a:.4f}, φ0={phi0:.3f}, V={V:.3f}")
plt.xlabel("Φ_θ  [rad]")
plt.ylabel("Intensity (a.u.)")
plt.title("θ–AB: Intensity vs Φ_θ with global cosine fit")
plt.legend()
plt.tight_layout()
plt.savefig("/mnt/data/fig_theta_ab_fit.png", dpi=160)

# Fig B: Cross-Hall odd component & fit
plt.figure()
plt.plot(t*1e3, cross_hall_fit["x_odd"], label="x_odd data")
plt.plot(t*1e3, cross_hall_fit["yhat"], label=f"fit: a={cross_hall_fit['a_hat']:.2e} ± {cross_hall_fit['sigma_a']:.1e}")
plt.xlabel("time [ms]")
plt.ylabel("x_odd [m]")
plt.title("Cross-Hall drift: odd-in-θ̇ component")
plt.legend()
plt.tight_layout()
plt.savefig("/mnt/data/fig_cross_hall_fit.png", dpi=160)

# Fig C: Sidebands (measured vs fit)
plt.figure()
plt.stem(n, E_meas/np.max(E_meas), linefmt='C0-', markerfmt='C0o', basefmt=" ", label="meas (norm)")
plt.stem(n+0.05, sideband_fit["E_fit"]/np.max(E_meas), linefmt='C1-', markerfmt='C1s', basefmt=" ", label="fit (shifted)")
plt.xlabel("Sideband index n")
plt.ylabel("Relative energy")
plt.title(f"Universal sidebands fit: φ̂={sideband_fit['phi_hat']:.3f}, Î={sideband_fit['I_hat']:.2e} kg·m²")
plt.legend()
plt.tight_layout()
plt.savefig("/mnt/data/fig_sidebands_fit.png", dpi=160)

# -----------------------------
# Summaries to CSV/Markdown
# -----------------------------

theta_ab_summary = {
    "alpha_hat": theta_ab_fit["alpha"],
    "phi0_hat":  theta_ab_fit["phi0"],
    "I0_hat":    theta_ab_fit["I0"],
    "V_hat":     theta_ab_fit["V"],
    "rss":       theta_ab_fit["rss"],
}

cross_hall_summary = {
    "a_hat": float(cross_hall_fit["a_hat"]),
    "sigma_a": float(cross_hall_fit["sigma_a"]),
}

sideband_summary = {
    "phi_hat": float(sideband_fit["phi_hat"]),
    "I_hat":   float(sideband_fit["I_hat"]),
    "c_hat":   float(sideband_fit["c_hat"]),
    "rss":     float(sideband_fit["rss"]),
}

summary_df = pd.DataFrame([
    {"section": "theta_ab", **theta_ab_summary},
    {"section": "cross_hall", **cross_hall_summary},
    {"section": "sidebands", **sideband_summary},
])
summary_path = "/mnt/data/fit_summary.csv"
summary_df.to_csv(summary_path, index=False)

display_dataframe_to_user("Fit summary (θ–AB, cross-Hall, sidebands)", summary_df)

# Markdown report
report = f"""# X–θ Simulation Fit Report (auto)

## θ–AB (intensity-only global fit)
- α̂ (slope vs Φ_θ): **{theta_ab_summary['alpha_hat']:.5f}**
- φ̂₀: **{theta_ab_summary['phi0_hat']:.4f} rad**
- Visibility V̂: **{theta_ab_summary['V_hat']:.4f}**
- RSS: {theta_ab_summary['rss']:.3e}

**Figure:** fig_theta_ab_fit.png

## Cross-Hall drift (odd-in-chirp estimator)
- â (proportional to F₀/m): **{cross_hall_summary['a_hat']:.3e} ± {cross_hall_summary['sigma_a']:.1e}** (1σ)

**Figure:** fig_cross_hall_fit.png

## Universal sidebands (rotor)
- φ̂ (Φ_θ / 2π): **{sideband_summary['phi_hat']:.4f}**
- Î (internal inertia): **{sideband_summary['I_hat']:.3e} kg·m²**
- ĉ = ħ²/(2Î): **{sideband_summary['c_hat']:.3e} J**
- RSS: {sideband_summary['rss']:.3e}

**Figure:** fig_sidebands_fit.png

---

### Notes
- θ–AB fit uses a grid over α and analytic linear solve for I₀, V, φ₀-equivalents; robust to small nonidealities.
- Cross-Hall fit cancels even-in-chirp backgrounds by working with x_odd = (x₊−x₋)/2 and regressing on the double integral of θ̇.
- Sideband fit scans φ ∈ [0,1) and solves for c at each; converts to Î = ħ²/(2c).

Replace the synthetic arrays with your lab CSVs (same shapes) and re-run to get publication-ready numbers and plots.
"""
with open("/mnt/data/Fit_Report.md", "w") as f:
    f.write(report)

"/mnt/data/fig_theta_ab_fit.png, fig_cross_hall_fit.png, fig_sidebands_fit.png, fit_summary.csv, Fit_Report.md are ready."

