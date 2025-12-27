#!/usr/bin/env python3
"""
weihs_paper_plot_v20.py

THE RIGOROUS PLOTTER + SINE FITTING
================================================================================
VERSION HISTORY:
  v20: Added WEIGHTED SINE FIT. 
       Now calculates goodness-of-fit (Chi-Squared, R^2) to test if the 
       oscillation is statistically significant compared to a flat line.
  v19: Adapted thresholds for low count rates (Set00).
  v18: Added Bin Population Diagnostics.

SCIENTIFIC CONTEXT:
  This script analyzes Bell Test data using the "X-Theta" framework.
  It overlays a sinusoidal fit S(theta) = c0 + A*sin(theta + phi)
  to quantify the phase dependence of the Bell parameter.

USAGE:
   python weihs_paper_plot_v20.py --alice-zip "A.zip" --bob-zip "B.zip" --reverse-bits
================================================================================
"""

import argparse
import zipfile
import numpy as np
import math
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from collections import deque
from typing import Dict, List, Tuple

# === CONFIGURATION ===
RESONANCE_KAPPA = 0.8       # The coupling constant for the geometric phase
BOOTSTRAP_ROUNDS = 2000     # Number of resamples for error bars
MIN_COUNTS_PER_SETTING = 5  # Low threshold for sparse data
BIN_COUNT = 12              # Number of angular bins for the plot

# ==========================================
# PART 1: CORE ENGINE & DATA LOADING
# ==========================================

def _norm_member(name: str) -> str:
    """Normalizes zip file paths."""
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"): return "/".join(parts[1:])
    return n

def load_raw_data(z, stem, reverse_bits):
    """Extracts Time (t), Setting (s), and Outcome (o)."""
    nm = { _norm_member(n).lower(): n for n in z.namelist() }
    vk, ck = (stem+"_V.dat").lower(), (stem+"_C.dat").lower()
    
    if vk not in nm:
        vk = next((k for k in nm if k.endswith(stem.lower() + "_v.dat")), None)
    if ck not in nm:
        ck = next((k for k in nm if k.endswith(stem.lower() + "_c.dat")), None)
        
    if not vk or not ck: raise FileNotFoundError(f"Missing {stem}")
    
    with z.open(nm[vk], "r") as f: vb = f.read()
    with z.open(nm[ck], "r") as f: cb = f.read()
    
    t = np.frombuffer(vb, dtype=">f8")
    c = np.frombuffer(cb, dtype=">u2")
    n = min(len(t), len(c))
    t = t[:n].astype(np.float64, copy=False)
    c = c[:n]
    
    if reverse_bits:
        det = ((c >> 1) & 1).astype(np.uint8)
        s   = (c & 1).astype(np.uint8)
    else:
        det = (c & 1).astype(np.uint8)
        s   = ((c >> 1) & 1).astype(np.uint8)
        
    o = (2 * det - 1).astype(np.int8)
    return t, s, o

State = Tuple[int, int]
def get_orientation(states: List[State]) -> int:
    """Calculates winding direction."""
    area = 0.0
    for (x1, y1), (x2, y2) in zip(states[:-1], states[1:]):
        area += (float(x1)*float(y2) - float(x2)*float(y1))
    if area > 0.1: return 1
    if area < -0.1: return -1
    return 0

@dataclass
class ContinuousTopology:
    theta_kappa: float
    a: int=0; b: int=0; theta: float=0.0
    history: deque = field(default_factory=lambda: deque(maxlen=5))
    bins: int = BIN_COUNT
    raw_outcomes: List[List[List[int]]] = field(default_factory=list)

    def __post_init__(self):
        self.raw_outcomes = [[[] for _ in range(4)] for _ in range(self.bins)]

    def update(self, actor, val):
        if actor == 'A': self.a = int(val)
        else: self.b = int(val)
        self.history.append((self.a, self.b))
        
        if len(self.history) == 5:
            h = list(self.history)
            if h[0] == h[-1] and len(set(h[:-1])) == 4:
                self.theta += self.theta_kappa * get_orientation(h)

    def capture(self, a, b, prod):
        ph = self.theta % (2 * math.pi)
        idx = int((ph / (2 * math.pi)) * self.bins) % self.bins
        k = a * 2 + b
        self.raw_outcomes[idx][k].append(prod)

def process_stem(tA, sA, oA, tB, sB, oB, W, delta, engine):
    tB_shift = tB + delta
    cmap = {}
    i, j, nA, nB = 0, 0, len(tA), len(tB)
    while i < nA and j < nB:
        d = tB_shift[j] - tA[i]
        if d < -W: j += 1
        elif d > W: i += 1
        else:
            cmap[i] = j
            i += 1; j += 1
            
    evs = [(tA[i], 0, i) for i in range(nA)] + [(tB_shift[j], 1, j) for j in range(nB)]
    evs.sort(key=lambda x: x[0])
    
    count = 0
    for _, src, idx in evs:
        if src == 0:
            engine.update('A', sA[idx])
            if idx in cmap:
                bi = cmap[idx]
                engine.capture(sA[idx], sB[bi], oA[idx] * oB[bi])
                count += 1
        else:
            engine.update('B', sB[idx])
    return count

# ==========================================
# PART 2: RIGOROUS STATISTICS & FITTING
# ==========================================

def calculate_S_robust(bin_outcomes, rounds=2000):
    counts = [len(bin_outcomes[k]) for k in range(4)]
    if any(c < MIN_COUNTS_PER_SETTING for c in counts):
        return None, None

    means = []
    for k in range(4):
        means.append(np.mean(bin_outcomes[k]))
    S_mean = means[0] + means[1] + means[2] - means[3]
    
    s_boots = []
    for _ in range(rounds):
        b_means = []
        for k in range(4):
            dat = np.array(bin_outcomes[k])
            resamp = np.random.choice(dat, size=len(dat), replace=True)
            b_means.append(np.mean(resamp))
        s_boots.append(b_means[0] + b_means[1] + b_means[2] - b_means[3])
    
    return S_mean, np.std(s_boots)

def weighted_linear_fit(X, y, sigma):
    """Weighted least squares for y ~ X @ beta, with weights 1/sigma^2."""
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    sigma = np.maximum(sigma, 1e-12)

    w = 1.0 / (sigma ** 2)
    sw = np.sqrt(w)
    Xw = X * sw[:, None]
    yw = y * sw

    XtX = Xw.T @ Xw
    Xty = Xw.T @ yw
    beta = np.linalg.solve(XtX, Xty)

    cov_beta = np.linalg.inv(XtX)
    yhat = X @ beta
    chi2 = np.sum(((y - yhat) / sigma) ** 2)
    dof = max(1, len(y) - X.shape[1])
    return beta, cov_beta, chi2, dof, yhat

def fit_sine_model(phases, s_means, s_errs):
    """Fits S(theta) = c0 + a sin(theta) + b cos(theta)."""
    th = np.asarray(phases, dtype=float)
    y  = np.asarray(s_means, dtype=float)
    se = np.asarray(s_errs, dtype=float)

    if len(y) < 4: return None

    # Sine model
    X = np.column_stack([np.ones_like(th), np.sin(th), np.cos(th)])
    beta, cov, chi2, dof, yhat = weighted_linear_fit(X, y, se)
    c0, a, b = beta

    A = float(np.hypot(a, b))
    phi = float(np.arctan2(b, a))

    # Approx amplitude uncertainty
    if A > 1e-12:
        da = a / A
        db = b / A
        varA = (da*da)*cov[1,1] + (db*db)*cov[2,2] + 2*(da*db)*cov[1,2]
        sigA = float(np.sqrt(max(0.0, varA)))
    else:
        sigA = float("nan")

    # Weighted R^2
    w = 1.0 / (np.maximum(se, 1e-12) ** 2)
    ybar = np.sum(w * y) / np.sum(w)
    sst = np.sum(w * (y - ybar) ** 2)
    sse = np.sum(w * (y - yhat) ** 2)
    r2w = 1.0 - (sse / sst if sst > 0 else np.nan)

    # Flat model baseline
    X0 = np.ones((len(th), 1))
    beta0, cov0, chi2_0, dof0, yhat0 = weighted_linear_fit(X0, y, se)

    return {
        "c0": float(c0), "a": float(a), "b": float(b),
        "A": A, "phi": phi, "sigA": sigA,
        "chi2": float(chi2), "dof": int(dof), "red_chi2": float(chi2 / dof),
        "r2w": float(r2w),
        "chi2_flat": float(chi2_0), "dof_flat": int(dof0), "red_chi2_flat": float(chi2_0 / dof0),
        "delta_chi2": float(chi2_0 - chi2)
    }

# ==========================================
# PART 3: MAIN EXECUTION
# ==========================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alice-zip", required=True)
    parser.add_argument("--bob-zip", required=True)
    parser.add_argument("--reverse-bits", action="store_true")
    args = parser.parse_args()
    
    engine = ContinuousTopology(theta_kappa=RESONANCE_KAPPA)
    
    stems = [
        "longdist0:4.5", "longdist1:-10.0", "longdist23:-2.0",
        "longdist3:-14.0", "longdist31:29.0", "longdist32:-5.0",
        "longdist35:16.0", "longdist36:-9.0", "longdist37:-6.0",
        "longdist5:-22.0"
    ]
    
    print("==========================================")
    print("    WEIH Paper Plotter v20 (Fitting)      ")
    print("==========================================")
    
    with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
        full_paths = []
        for n in za.namelist():
             if n.lower().endswith("_v.dat"):
                 cl = n.replace("\\", "/").lstrip("/")
                 if cl.lower().startswith("alice/"): cl = cl[6:]
                 full_paths.append(cl[:-6])
        
        for item in stems:
            short, d_str = item.split(":")
            match = next((p for p in full_paths if p.endswith(short)), None)
            if not match: continue
            
            tA, sA, oA = load_raw_data(za, match, args.reverse_bits)
            tB, sB, oB = load_raw_data(zb, match, args.reverse_bits)
            process_stem(tA, sA, oA, tB, sB, oB, 1.3e-9, float(d_str)*1e-9, engine)

    # Compile Data
    phases, s_means, s_errs = [], [], []
    
    print("\n[ANALYSIS] Calculating S parameters...")
    print(f"   Phase | S_mean | S_err (BS={BOOTSTRAP_ROUNDS})")
    print("   --------------------------------------")
    
    for i in range(engine.bins):
        S, err = calculate_S_robust(engine.raw_outcomes[i], BOOTSTRAP_ROUNDS)
        if S is None: continue
        ph = (i + 0.5) * (2*math.pi/engine.bins)
        phases.append(ph)
        s_means.append(S)
        s_errs.append(err)
        print(f"   {ph:4.2f}  | {S:5.2f}  | {err:5.2f}")

    # --- FIT SINE WAVE ---
    fit = None
    if len(phases) >= 4:
        fit = fit_sine_model(phases, s_means, s_errs)

    if fit:
        print("\n[SINE FIT] Weighted fit: S(theta)=c0 + A sin(theta + phi)")
        print(f"  Amplitude A  = {fit['A']:.4f} (± {fit['sigA']:.4f})")
        print(f"  Offset c0    = {fit['c0']:.4f}")
        print(f"  Phase phi    = {fit['phi']:.4f} rad")
        print(f"  Goodness: Red. Chi2 = {fit['red_chi2']:.3f}, Weighted R2 = {fit['r2w']:.3f}")
        print(f"  Comparison: Delta Chi2 (Flat - Sine) = {fit['delta_chi2']:.2f}")

    # PLOT
    plt.figure(figsize=(10, 6))
    
    plt.axhline(2.0, color='red', linestyle='--', linewidth=1.5, label='Classical Bound (S=2)')
    plt.axhline(2.828, color='green', linestyle=':', linewidth=1.5, label='Quantum Limit')
    
    # Plot Fit Curve
    if fit:
        th_dense = np.linspace(0, 2*math.pi, 400)
        y_dense = fit["c0"] + fit["a"]*np.sin(th_dense) + fit["b"]*np.cos(th_dense)
        label_fit = f"Sine Fit ($R_w^2={fit['r2w']:.2f}, \Delta\chi^2={fit['delta_chi2']:.1f}$)"
        plt.plot(th_dense, y_dense, color='orange', linewidth=2, label=label_fit)

    # Plot Data
    plt.errorbar(phases, s_means, yerr=s_errs, fmt='o', color='blue', 
                 ecolor='black', capsize=5, label='X-Theta Data')
    
    plt.xlabel('Geometric Phase $\\theta$ (radians)')
    plt.ylabel('Bell Parameter $S$')
    plt.title(f'Phase-Dependent Bell Violation ($\kappa={RESONANCE_KAPPA}$)')
    plt.legend(loc='lower right')
    plt.grid(True, alpha=0.3)
    plt.ylim(-0.5, 3.5)
    
    plt.savefig("weihs_paper_plot_v20.png")
    print("\nSaved plot to: weihs_paper_plot_v20.png")

if __name__ == "__main__":
    main()