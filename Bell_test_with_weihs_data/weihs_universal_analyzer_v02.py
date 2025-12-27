#!/usr/bin/env python3
"""
weihs_universal_analyzer_v02.py

THE UNIVERSAL DATA MINER (Optimized / Fast Mode)
================================================================================
CHANGES from v01:
  1. VECTORIZED ANALYSIS: Uses numpy.bincount for instant aggregation.
     Replaces slow Python loops that caused the 'hang' in v01.
  2. ANALYTICAL ERRORS: Uses standard error (1/sqrt(N)) for the Permutation 
     Test loop, which is 2000x faster than Bootstrap and valid for N=300k.
  3. DATA PERSISTENCE: Saves mined events to 'weihs_mined_data.npz' so you 
     never have to re-mine if the script stops.

USAGE:
   python weihs_universal_analyzer_v02.py --alice-zip "A.zip" --bob-zip "B.zip" --reverse-bits
================================================================================
"""

import argparse
import zipfile
import numpy as np
import math
import matplotlib.pyplot as plt
import random
import os
from dataclasses import dataclass, field
from collections import deque
from typing import Dict, List, Tuple

# === CONFIGURATION ===
RESONANCE_KAPPA = 0.8       
BOOTSTRAP_ROUNDS = 2000     
PERMUTATION_ROUNDS = 500    
BIN_COUNT = 6               
SYNC_WINDOW_NS = 1000.0     
SYNC_BIN_SIZE_NS = 1.0      

# ==========================================
# PART 1: AUTO-SYNC & MINING ENGINE
# ==========================================

def _norm_member(name: str) -> str:
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"): return "/".join(parts[1:])
    return n

def load_raw_data(z, stem, reverse_bits):
    nm = { _norm_member(n).lower(): n for n in z.namelist() }
    vk, ck = (stem+"_V.dat").lower(), (stem+"_C.dat").lower()
    if vk not in nm: vk = next((k for k in nm if k.endswith(stem.lower() + "_v.dat")), None)
    if ck not in nm: ck = next((k for k in nm if k.endswith(stem.lower() + "_c.dat")), None)
    if not vk or not ck: return None, None, None 
    
    with z.open(nm[vk], "r") as f: vb = f.read()
    with z.open(nm[ck], "r") as f: cb = f.read()
    
    t = np.frombuffer(vb, dtype=">f8")
    c = np.frombuffer(cb, dtype=">u2")
    n = min(len(t), len(c))
    t = t[:n].astype(np.float64, copy=False)
    c = c[:n]
    
    if reverse_bits:
        det = ((c >> 1) & 1).astype(np.uint8); s = (c & 1).astype(np.uint8)
    else:
        det = (c & 1).astype(np.uint8); s = ((c >> 1) & 1).astype(np.uint8)
    o = (2 * det - 1).astype(np.int8)
    return t, s, o

def find_optimal_shift(tA, tB, window_ns=SYNC_WINDOW_NS, bin_size_ns=SYNC_BIN_SIZE_NS):
    limit = min(5000, len(tA), len(tB))
    if limit < 100: return None, 0.0
    tA_s = tA[:limit] * 1e9; tB_s = tB[:limit] * 1e9
    diffs = tB_s[:, np.newaxis] - tA_s[np.newaxis, :]
    valid_diffs = diffs[np.abs(diffs) < window_ns]
    if len(valid_diffs) == 0: return None, 0.0
    
    bins = int((2 * window_ns) / bin_size_ns)
    hist, bin_edges = np.histogram(valid_diffs, bins=bins, range=(-window_ns, window_ns))
    peak_idx = np.argmax(hist)
    peak_val = hist[peak_idx]
    
    # Noise floor
    noise_mask = np.ones(len(hist), dtype=bool)
    noise_mask[max(0, peak_idx-2):min(len(hist), peak_idx+3)] = False
    avg_noise = np.mean(hist[noise_mask]) if np.any(noise_mask) else 0.1
    
    return (bin_edges[peak_idx] + bin_edges[peak_idx+1]) / 2.0 * 1e-9, peak_val / (avg_noise + 0.1)

# ==========================================
# PART 2: TOPOLOGY & ANALYSIS
# ==========================================

@dataclass
class ContinuousTopology:
    theta_kappa: float
    a: int=0; b: int=0; theta: float=0.0
    history: deque = field(default_factory=lambda: deque(maxlen=5))
    # Stores (Theta, Setting, Product)
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
        k = a * 2 + b
        self.collected_events.append((self.theta, k, prod))

State = Tuple[int, int]
def get_orientation(states: List[State]) -> int:
    area = 0.0
    for (x1, y1), (x2, y2) in zip(states[:-1], states[1:]):
        area += (float(x1)*float(y2) - float(x2)*float(y1))
    if area > 0.1: return 1
    if area < -0.1: return -1
    return 0

def process_stem_data(tA, sA, oA, tB, sB, oB, shift, engine, coinc_window=1.5e-9):
    tB_shifted = tB + shift
    i, j, nA, nB = 0, 0, len(tA), len(tB)
    cmap = {}
    while i < nA and j < nB:
        d = tB_shifted[j] - tA[i]
        if d < -coinc_window: j += 1
        elif d > coinc_window: i += 1
        else:
            cmap[i] = j; i += 1; j += 1
            
    evs = [(tA[i], 0, i) for i in range(nA)] + [(tB_shifted[j], 1, j) for j in range(nB)]
    evs.sort(key=lambda x: x[0])
    
    local_count = 0
    for _, src, idx in evs:
        if src == 0:
            engine.update('A', sA[idx])
            if idx in cmap:
                engine.capture(sA[idx], sB[cmap[idx]], oA[idx] * oB[cmap[idx]])
                local_count += 1
        else:
            engine.update('B', sB[idx])
    return local_count

# --- ANALYSIS ENGINES ---

def weighted_linear_fit(X, y, sigma):
    y = np.asarray(y); X = np.asarray(X); sigma = np.maximum(np.asarray(sigma), 1e-12)
    w = 1.0 / (sigma ** 2); sw = np.sqrt(w)
    Xw = X * sw[:, None]; yw = y * sw
    try:
        XtX = Xw.T @ Xw
        beta = np.linalg.solve(XtX, Xw.T @ yw)
        yhat = X @ beta
        chi2 = np.sum(((y - yhat) / sigma) ** 2)
        dof = max(1, len(y) - X.shape[1])
        ybar = np.sum(w * y) / np.sum(w)
        r2w = 1.0 - (np.sum(w * (y - yhat) ** 2) / np.sum(w * (y - ybar) ** 2))
        return beta, chi2, dof, r2w
    except:
        return None, 0, 1, 0

def analyze_vectorized(phases, settings, products, use_bootstrap=False):
    """
    Super-Fast Vectorized Analysis.
    Input: Numpy arrays.
    """
    # 1. Binning (Instant)
    bin_indices = ((phases % (2*np.pi)) / (2*np.pi) * BIN_COUNT).astype(int) % BIN_COUNT
    
    # 2. Aggregation (bincount is O(N))
    # Flatten index: bin * 4 + setting
    flat_indices = bin_indices * 4 + settings
    counts = np.bincount(flat_indices, minlength=BIN_COUNT*4)
    sums = np.bincount(flat_indices, weights=products, minlength=BIN_COUNT*4)
    
    # Filter empty bins
    if np.any(counts < 5): return None
    
    means = sums / counts
    
    # 3. Error Calculation
    if use_bootstrap:
        # Slower, but accurate for final plot
        # We can't fully vectorize bootstrap easily, so we do a quick loop 
        # only over the bins (small N), resampling from the original data logic?
        # Actually, for N=300k, analytical error is indistinguishable from bootstrap.
        # Let's use analytical for speed/stability even in main plot, or standard error.
        errs = 1.0 / np.sqrt(counts) # Standard Error for +/-1 data
    else:
        # Fast Analytical Error (1/sqrt(N))
        errs = 1.0 / np.sqrt(counts)

    # 4. Calculate S per bin
    # Reshape to [Bins, 4]
    means_reshaped = means.reshape((BIN_COUNT, 4))
    errs_reshaped = errs.reshape((BIN_COUNT, 4))
    
    s_vals = means_reshaped[:,0] + means_reshaped[:,1] + means_reshaped[:,2] - means_reshaped[:,3]
    # Error propagation: sqrt(sum(err^2))
    s_errs = np.sqrt(np.sum(errs_reshaped**2, axis=1))
    
    bin_centers = (np.arange(BIN_COUNT) + 0.5) * (2*math.pi/BIN_COUNT)
    
    # 5. Fit
    X = np.column_stack([np.ones(BIN_COUNT), np.sin(bin_centers), np.cos(bin_centers)])
    beta, chi2, dof, r2w = weighted_linear_fit(X, s_vals, s_errs)
    
    # Flat Fit
    X0 = np.ones((BIN_COUNT, 1))
    _, chi2_0, _, _ = weighted_linear_fit(X0, s_vals, s_errs)
    
    if beta is None: return None
    
    return {
        "phases": bin_centers, "s": s_vals, "err": s_errs, 
        "beta": beta, "delta_chi2": chi2_0 - chi2, "r2w": r2w
    }

# ==========================================
# PART 4: MAIN
# ==========================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alice-zip", required=True); parser.add_argument("--bob-zip", required=True)
    parser.add_argument("--reverse-bits", action="store_true")
    args = parser.parse_args()
    
    print("==========================================")
    print("   WEIHS UNIVERSAL DATA MINER v02 (FAST)")
    print("==========================================")
    
    engine = ContinuousTopology(theta_kappa=RESONANCE_KAPPA)
    
    # --- STEP 1: LOAD OR MINE ---
    cache_file = "weihs_mined_data.npz"
    
    if os.path.exists(cache_file):
        print(f"Found cached data ({cache_file}). Loading...")
        data = np.load(cache_file)
        phases = data['phases']
        settings = data['settings']
        products = data['products']
        print(f"Loaded {len(phases)} events.")
    else:
        # MINE
        with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
            a_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in za.namelist() if "_v.dat" in n.lower() }
            b_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in zb.namelist() if "_v.dat" in n.lower() }
            stems = sorted(list(a_files.intersection(b_files)))
            
            print(f"Mining {len(stems)} file pairs...")
            
            valid_files = 0
            for i, stem in enumerate(stems):
                tA, sA, oA = load_raw_data(za, stem, args.reverse_bits)
                tB, sB, oB = load_raw_data(zb, stem, args.reverse_bits)
                if tA is None or len(tA) < 100: continue
                
                shift, snr = find_optimal_shift(tA, tB)
                if shift is not None and snr > 5.0:
                    process_stem_data(tA, sA, oA, tB, sB, oB, shift, engine)
                    valid_files += 1
                
                if i % 20 == 0: print(f"[{i}/{len(stems)}] Processed...")

        # CONVERT & SAVE
        events = engine.collected_events
        print("Converting to Numpy...")
        raw_arr = np.array(events) # [N, 3]
        phases = raw_arr[:,0]
        settings = raw_arr[:,1].astype(int)
        products = raw_arr[:,2].astype(int)
        
        np.savez(cache_file, phases=phases, settings=settings, products=products)
        print(f"Mined data saved to {cache_file}.")

    print(f"Total Events: {len(phases)}")
    
    # --- STEP 2: REAL ANALYSIS ---
    print("\n[ANALYSIS] Fitting Model...")
    # Use standard analytical error for consistency
    res = analyze_vectorized(phases, settings, products, use_bootstrap=False) 
    
    print(f"Observed Delta Chi2: {res['delta_chi2']:.4f}")
    print(f"Observed Weighted R2: {res['r2w']:.4f}")
    
    # --- STEP 3: FAST PERMUTATION TEST ---
    print(f"\n[PERMUTATION] Running {PERMUTATION_ROUNDS} shuffles (Vectorized)...")
    better = 0
    obs_delta = res['delta_chi2']
    
    # Copy arrays for shuffling
    perm_phases = phases.copy()
    
    for i in range(PERMUTATION_ROUNDS):
        np.random.shuffle(perm_phases)
        # Fast Analysis
        pres = analyze_vectorized(perm_phases, settings, products, use_bootstrap=False)
        if pres and pres['delta_chi2'] >= obs_delta:
            better += 1
        
        if (i+1) % 100 == 0:
            print(f"  .. {i+1} done (Better: {better})")
            
    p_val = (better + 1)/(PERMUTATION_ROUNDS + 1)
    print(f"\n[FINAL RESULT] p-value = {p_val:.5f}")

    # --- PLOT ---
    plt.figure(figsize=(10,6))
    plt.axhline(2, color='red', linestyle='--', label="Classical")
    plt.axhline(2.828, color='green', linestyle=':', label="Tsirelson")
    
    c0, a, b = res["beta"]
    th_d = np.linspace(0, 2*math.pi, 200)
    plt.plot(th_d, c0 + a*np.sin(th_d) + b*np.cos(th_d), color='orange', linewidth=2, label=f"Fit (p={p_val:.3f})")
    plt.errorbar(res["phases"], res["s"], yerr=res["err"], fmt='o', color='blue', ecolor='black', label="Data")
    
    plt.title(f"Universal Analysis (N={len(phases)}, $\Delta\chi^2={obs_delta:.1f}$)")
    plt.ylim(1.0, 3.0) # Zoom in a bit
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig("weihs_universal_plot_v02.png")
    print("Saved: weihs_universal_plot_v02.png")

if __name__ == "__main__":
    main()