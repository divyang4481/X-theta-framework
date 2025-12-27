#!/usr/bin/env python3
"""
weihs_paper_plot_v21.py

THE RIGOROUS PLOTTER + PERMUTATION TESTING
================================================================================
VERSION HISTORY:
  v21: Added PERMUTATION TEST (p-value).
       - Decouples Data Collection from Analysis.
       - Runs 500 Monte Carlo permutations (shuffling Theta vs Outcome).
       - Calculates empirical p-value for the Sine Fit improvement.
       - Updated Tsirelson label.
  v20: Added Weighted Sine Fit.

SCIENTIFIC CONTEXT:
  This script tests the "X-Theta" hypothesis by permuting the geometric phases.
  Null Hypothesis (H0): The Bell parameter S is independent of the geometric phase.
  
  If H0 is true, shuffling the theta values among events should yield a fit 
  quality (Delta Chi2) similar to the original data.

USAGE:
   python weihs_paper_plot_v21.py --alice-zip "A.zip" --bob-zip "B.zip" --reverse-bits
=========================== =====================================================
"""

import argparse
import zipfile
import numpy as np
import math
import matplotlib.pyplot as plt
import random
from dataclasses import dataclass, field
from collections import deque
from typing import Dict, List, Tuple

# === CONFIGURATION ===
RESONANCE_KAPPA = 0.8       
BOOTSTRAP_ROUNDS = 2000     
MIN_COUNTS_PER_SETTING = 5  
BIN_COUNT = 6              
PERMUTATION_ROUNDS = 500    # Number of null-hypothesis shuffles

# ==========================================
# PART 1: CORE ENGINE & DATA EXTRACTION
# ==========================================

def _norm_member(name: str) -> str:
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"): return "/".join(parts[1:])
    return n

def load_raw_data(z, stem, reverse_bits):
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
    
    # Store RAW events instead of binning immediately
    # List of (Theta, SettingIndex, Product)
    collected_events: List[Tuple[float, int, int]] = field(default_factory=list)

    def update(self, actor, val):
        if actor == 'A': self.a = int(val)
        else: self.b = int(val)
        self.history.append((self.a, self.b))
        
        if len(self.history) == 5:
            h = list(self.history)
            if h[0] == h[-1] and len(set(h[:-1])) == 4:
                self.theta += self.theta_kappa * get_orientation(h)

    def capture(self, a, b, prod):
        # Just store the raw data point
        k = a * 2 + b
        self.collected_events.append((self.theta, k, prod))

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
    
    for _, src, idx in evs:
        if src == 0:
            engine.update('A', sA[idx])
            if idx in cmap:
                bi = cmap[idx]
                engine.capture(sA[idx], sB[bi], oA[idx] * oB[bi])
        else:
            engine.update('B', sB[idx])

# ==========================================
# PART 2: ANALYSIS PIPELINE (Binning -> S -> Fit)
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
    try:
        beta = np.linalg.solve(XtX, Xty)
        cov_beta = np.linalg.inv(XtX)
    except np.linalg.LinAlgError:
        return None, None, float('inf'), 0, np.zeros_like(y)

    yhat = X @ beta
    chi2 = np.sum(((y - yhat) / sigma) ** 2)
    dof = max(1, len(y) - X.shape[1])
    return beta, cov_beta, chi2, dof, yhat

def run_full_analysis(events: List[Tuple[float, int, int]]):
    """
    Takes a list of (Theta, Setting, Product), bins them, calculates S, fits Sine.
    Returns dictionary of results.
    """
    # 1. Binning
    raw_outcomes = [[[] for _ in range(4)] for _ in range(BIN_COUNT)]
    for theta, k, prod in events:
        ph = theta % (2 * math.pi)
        idx = int((ph / (2 * math.pi)) * BIN_COUNT) % BIN_COUNT
        raw_outcomes[idx][k].append(prod)

    # 2. Calculate S
    phases, s_means, s_errs = [], [], []
    for i in range(BIN_COUNT):
        S, err = calculate_S_robust(raw_outcomes[i], BOOTSTRAP_ROUNDS)
        if S is not None:
            phases.append((i + 0.5) * (2*math.pi/BIN_COUNT))
            s_means.append(S)
            s_errs.append(err)

    if len(phases) < 4:
        return None

    # 3. Fit Sine
    th = np.asarray(phases)
    y  = np.asarray(s_means)
    se = np.asarray(s_errs)

    # Sine Model
    X = np.column_stack([np.ones_like(th), np.sin(th), np.cos(th)])
    beta, cov, chi2, dof, yhat = weighted_linear_fit(X, y, se)
    if beta is None: return None
    
    c0, a, b = beta
    A = float(np.hypot(a, b))
    phi = float(np.arctan2(b, a))
    
    # Flat Model
    X0 = np.ones((len(th), 1))
    beta0, cov0, chi2_0, dof0, yhat0 = weighted_linear_fit(X0, y, se)

    # R2 weighted
    w = 1.0 / (np.maximum(se, 1e-12) ** 2)
    ybar = np.sum(w * y) / np.sum(w)
    sst = np.sum(w * (y - ybar) ** 2)
    sse = np.sum(w * (y - yhat) ** 2)
    r2w = 1.0 - (sse / sst if sst > 0 else np.nan)

    return {
        "phases": phases, "s_means": s_means, "s_errs": s_errs,
        "c0": c0, "a": a, "b": b, "A": A, "phi": phi,
        "chi2": chi2, "red_chi2": chi2/dof,
        "delta_chi2": chi2_0 - chi2,
        "r2w": r2w
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
    print("    WEIH Paper Plotter v21 (Permutation)  ")
    print("==========================================")
    
    # --- LOAD DATA ---
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

    events = engine.collected_events
    print(f"Total Events Collected: {len(events)}")

    # --- REAL ANALYSIS ---
    print("\n[REAL DATA] Running Analysis...")
    real_res = run_full_analysis(events)
    if not real_res:
        print("Error: Insufficient data for real analysis.")
        return

    obs_delta_chi2 = real_res["delta_chi2"]
    print(f"Observed Delta Chi2: {obs_delta_chi2:.4f}")
    print(f"Observed Weighted R2: {real_res['r2w']:.4f}")

    # --- PERMUTATION TEST ---
    print(f"\n[PERMUTATION] Running {PERMUTATION_ROUNDS} shuffles (Scrambling Theta)...")
    
    better_count = 0
    # Extract columns for shuffling
    thetas = [e[0] for e in events]
    fixed_data = [(e[1], e[2]) for e in events] # (Setting, Product) pairs kept together
    
    for i in range(PERMUTATION_ROUNDS):
        # Shuffle Thetas ONLY. This destroys the link between history(theta) and outcome(S)
        # while preserving the global Bell violation and count statistics.
        shuffled_thetas = list(thetas)
        random.shuffle(shuffled_thetas)
        
        # Reconstruct event list
        perm_events = []
        for j in range(len(events)):
            perm_events.append((shuffled_thetas[j], fixed_data[j][0], fixed_data[j][1]))
            
        perm_res = run_full_analysis(perm_events)
        
        if perm_res:
            null_delta = perm_res["delta_chi2"]
            if null_delta >= obs_delta_chi2:
                better_count += 1
        
        if (i+1) % 50 == 0:
            print(f"  .. {i+1}/{PERMUTATION_ROUNDS} done (Better so far: {better_count})")

    p_value = (better_count + 1) / (PERMUTATION_ROUNDS + 1)
    print(f"\n[RESULT] Empirical p-value = {p_value:.4f}")
    
    # --- PLOT ---
    plt.figure(figsize=(10, 6))
    
    plt.axhline(2.0, color='red', linestyle='--', linewidth=1.5, label='Classical Bound (S=2)')
    plt.axhline(2.828, color='green', linestyle=':', linewidth=1.5, label='Tsirelson Bound ($2\\sqrt{2}$)')
    
    # Plot Fit
    th_dense = np.linspace(0, 2*math.pi, 400)
    fit = real_res
    y_dense = fit["c0"] + fit["a"]*np.sin(th_dense) + fit["b"]*np.cos(th_dense)
    
    label_fit = f"Sine Fit ($p={p_value:.3f}, R_w^2={fit['r2w']:.2f}$)"
    plt.plot(th_dense, y_dense, color='orange', linewidth=2, label=label_fit)

    # Plot Data
    plt.errorbar(fit["phases"], fit["s_means"], yerr=fit["s_errs"], fmt='o', color='blue', 
                 ecolor='black', capsize=5, label='X-Theta Data')
    
    plt.xlabel('Geometric Phase $\\theta$ (radians)')
    plt.ylabel('Bell Parameter $S$')
    plt.title(f'Phase-Dependent Bell Violation ($\kappa={RESONANCE_KAPPA}$)\nPermutation Test: p={p_value:.3f} (N={PERMUTATION_ROUNDS})')
    plt.legend(loc='lower right')
    plt.grid(True, alpha=0.3)
    plt.ylim(-0.5, 3.5)
    
    plt.savefig("weihs_paper_plot_v21.png")
    print("\nSaved plot to: weihs_paper_plot_v21.png")

if __name__ == "__main__":
    main()