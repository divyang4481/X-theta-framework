#!/usr/bin/env python3
"""
weihs_universal_analyzer_v01.py

THE UNIVERSAL DATA MINER & ANALYZER
================================================================================
DESCRIPTION:
  This script moves beyond manual "stems" and attempts to mine the ENTIRE 
  Weihs dataset (~400+ files).
  
  1. AUTO-DISCOVERY: Scans Zip files for matching Alice/Bob pairs.
  2. AUTO-SYNC: Uses cross-correlation to mathematically find the time shift 
     (delta) for each file pair automatically.
  3. UNIVERSAL ANALYSIS: Feeds all valid, synchronized data into the 
     X-Theta Continuous Topology engine.
  4. PERMUTATION TEST: Runs the rigorous shuffling test on the massive dataset.

ALGORITHM (Time Sync):
  For each file pair, we calculate the time differences (t_Bob - t_Alice) 
  for the first 5000 events. We histogram these differences. 
  The "Peak" of the histogram corresponds to the correct optical delay. 
  If the peak is sharp (high Signal-to-Noise), we keep the file.

USAGE:
   python weihs_universal_analyzer_v01.py --alice-zip "A.zip" --bob-zip "B.zip" --reverse-bits
================================================================================
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
PERMUTATION_ROUNDS = 500    
BIN_COUNT = 6               # Kept at 6 for robust statistics
MIN_COUNTS_PER_SETTING = 5
SYNC_WINDOW_NS = 1000.0     # Search for time shifts within +/- 1000 ns
SYNC_BIN_SIZE_NS = 1.0      # Resolution of time shift search

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
    
    # Robust lookup
    if vk not in nm: vk = next((k for k in nm if k.endswith(stem.lower() + "_v.dat")), None)
    if ck not in nm: ck = next((k for k in nm if k.endswith(stem.lower() + "_c.dat")), None)
    if not vk or not ck: return None, None, None # Fail gracefully
    
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
    """
    Automatic Time-Synchronization.
    Returns: (best_shift_seconds, peak_quality_score)
    """
    # Use a slice of data for speed (first 5000 events)
    limit = min(5000, len(tA), len(tB))
    if limit < 100: return None, 0.0
    
    tA_s = tA[:limit] * 1e9 # Convert to ns for histogramming
    tB_s = tB[:limit] * 1e9
    
    # Calculate all differences within window
    # We use broadcasting for a brute-force check on the slice (efficient enough for N=5000)
    # diffs = tB - tA
    # Matrix subtraction:
    diffs = tB_s[:, np.newaxis] - tA_s[np.newaxis, :]
    
    # Flatten and filter by window
    valid_mask = np.abs(diffs) < window_ns
    valid_diffs = diffs[valid_mask]
    
    if len(valid_diffs) == 0: return None, 0.0
    
    # Histogram
    bins = int((2 * window_ns) / bin_size_ns)
    hist, bin_edges = np.histogram(valid_diffs, bins=bins, range=(-window_ns, window_ns))
    
    peak_idx = np.argmax(hist)
    peak_val = hist[peak_idx]
    
    # Calculate noise floor (average of non-peak bins)
    # Exclude the peak region
    noise_mask = np.ones(len(hist), dtype=bool)
    noise_mask[max(0, peak_idx-2):min(len(hist), peak_idx+3)] = False
    avg_noise = np.mean(hist[noise_mask]) if np.any(noise_mask) else 0.1
    
    signal_to_noise = peak_val / (avg_noise + 0.1)
    
    # Refine peak center
    center_ns = (bin_edges[peak_idx] + bin_edges[peak_idx+1]) / 2.0
    
    return center_ns * 1e-9, signal_to_noise

# ==========================================
# PART 2: CONTINUOUS TOPOLOGY ENGINE
# ==========================================

@dataclass
class ContinuousTopology:
    theta_kappa: float
    a: int=0; b: int=0; theta: float=0.0
    history: deque = field(default_factory=lambda: deque(maxlen=5))
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
    
    # Linear Scan for Coincidences
    i, j = 0, 0
    nA, nB = len(tA), len(tB)
    cmap = {}
    
    while i < nA and j < nB:
        d = tB_shifted[j] - tA[i]
        if d < -coinc_window: j += 1
        elif d > coinc_window: i += 1
        else:
            cmap[i] = j
            i += 1; j += 1
            
    # Chronological Replay
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

# ==========================================
# PART 3: ANALYSIS & STATISTICS
# ==========================================

def calculate_S_robust(bin_outcomes, rounds=2000):
    counts = [len(bin_outcomes[k]) for k in range(4)]
    if any(c < MIN_COUNTS_PER_SETTING for c in counts): return None, None
    means = [np.mean(bin_outcomes[k]) for k in range(4)]
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
    y = np.asarray(y); X = np.asarray(X); sigma = np.maximum(np.asarray(sigma), 1e-12)
    w = 1.0 / (sigma ** 2); sw = np.sqrt(w)
    Xw = X * sw[:, None]; yw = y * sw
    try:
        XtX = Xw.T @ Xw
        beta = np.linalg.solve(XtX, Xw.T @ yw)
        yhat = X @ beta
        chi2 = np.sum(((y - yhat) / sigma) ** 2)
        dof = max(1, len(y) - X.shape[1])
        # R2
        ybar = np.sum(w * y) / np.sum(w)
        sst = np.sum(w * (y - ybar) ** 2)
        sse = np.sum(w * (y - yhat) ** 2)
        r2w = 1.0 - (sse / sst if sst > 0 else 0)
        return beta, chi2, dof, r2w
    except:
        return None, 0, 1, 0

def analyze_events(events):
    # Binning
    raw = [[[] for _ in range(4)] for _ in range(BIN_COUNT)]
    for th, k, p in events:
        idx = int((th % (2*math.pi)) / (2*math.pi) * BIN_COUNT) % BIN_COUNT
        raw[idx][k].append(p)
        
    phases, s_means, s_errs = [], [], []
    for i in range(BIN_COUNT):
        S, err = calculate_S_robust(raw[i], BOOTSTRAP_ROUNDS)
        if S:
            phases.append((i+0.5)*(2*math.pi/BIN_COUNT))
            s_means.append(S); s_errs.append(err)
            
    if len(phases) < 4: return None
    
    # Fit
    th = np.array(phases)
    X = np.column_stack([np.ones_like(th), np.sin(th), np.cos(th)])
    beta, chi2, dof, r2w = weighted_linear_fit(X, s_means, s_errs)
    
    # Flat
    X0 = np.ones((len(th), 1))
    _, chi2_0, _, _ = weighted_linear_fit(X0, s_means, s_errs)
    
    if beta is None: return None
    return {"r2w": r2w, "delta_chi2": chi2_0 - chi2, "phases": phases, "s": s_means, "err": s_errs, "beta": beta}

# ==========================================
# PART 4: MAIN
# ==========================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alice-zip", required=True); parser.add_argument("--bob-zip", required=True)
    parser.add_argument("--reverse-bits", action="store_true")
    args = parser.parse_args()
    
    print("==========================================")
    print("   WEIHS UNIVERSAL DATA MINER v01")
    print("==========================================")
    
    engine = ContinuousTopology(theta_kappa=RESONANCE_KAPPA)
    
    # 1. Discover Pairs
    with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
        a_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in za.namelist() if "_v.dat" in n.lower() }
        b_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in zb.namelist() if "_v.dat" in n.lower() }
        stems = sorted(list(a_files.intersection(b_files)))
        
        print(f"Found {len(stems)} potential file pairs. Starting Mining Operation...")
        print(f"Sync Logic: Searching +/- {SYNC_WINDOW_NS}ns for coincidence peaks.")
        
        valid_files = 0
        total_coinc = 0
        
        for i, stem in enumerate(stems):
            tA, sA, oA = load_raw_data(za, stem, args.reverse_bits)
            tB, sB, oB = load_raw_data(zb, stem, args.reverse_bits)
            
            if tA is None or len(tA) < 100: continue
            
            # Auto-Sync
            shift, snr = find_optimal_shift(tA, tB)
            
            status = "SKIP"
            if shift is not None and snr > 5.0: # Threshold for valid sync
                c = process_stem_data(tA, sA, oA, tB, sB, oB, shift, engine)
                total_coinc += c
                valid_files += 1
                status = f"LOCK ({c} ev)"
            
            if i % 10 == 0:
                print(f"[{i+1}/{len(stems)}] {stem[:15]}... \tSNR={snr:4.1f} -> {status}")

    print(f"\nMINING COMPLETE.")
    print(f"Valid Files: {valid_files} / {len(stems)}")
    print(f"Total Events: {total_coinc}")
    
    if total_coinc < 100:
        print("Error: Not enough data mined. Check Sync Parameters.")
        return

    # 2. Analyze
    print("\n[ANALYSIS] Fitting Model...")
    res = analyze_events(engine.collected_events)
    if not res: return
    
    print(f"Observed Delta Chi2: {res['delta_chi2']:.4f}")
    print(f"Observed Weighted R2: {res['r2w']:.4f}")
    
    # 3. Permutation Test
    print(f"\n[PERMUTATION] Running {PERMUTATION_ROUNDS} shuffles...")
    better = 0
    events = engine.collected_events
    thetas = [e[0] for e in events]
    fixed = [(e[1], e[2]) for e in events]
    
    for i in range(PERMUTATION_ROUNDS):
        random.shuffle(thetas)
        perm_ev = [(thetas[j], fixed[j][0], fixed[j][1]) for j in range(len(events))]
        pres = analyze_events(perm_ev)
        if pres and pres["delta_chi2"] >= res["delta_chi2"]: better += 1
        if i % 50 == 0: print(f"..{i} (Better: {better})")
            
    p_val = (better + 1)/(PERMUTATION_ROUNDS + 1)
    print(f"\n[FINAL RESULT] p-value = {p_val:.5f}")
    
    # 4. Plot
    plt.figure(figsize=(10,6))
    plt.axhline(2, color='red', linestyle='--'); plt.axhline(2.828, color='green', linestyle=':')
    
    th_d = np.linspace(0, 2*math.pi, 200)
    c0, a, b = res["beta"]
    plt.plot(th_d, c0 + a*np.sin(th_d) + b*np.cos(th_d), color='orange', linewidth=2, label=f"Fit (p={p_val:.3f})")
    plt.errorbar(res["phases"], res["s"], yerr=res["err"], fmt='o', color='blue', ecolor='black', label="Data")
    
    plt.title(f"Universal Analysis (N={total_coinc}, Files={valid_files})")
    plt.ylim(-0.5, 3.5); plt.legend(); plt.grid(True, alpha=0.3)
    plt.savefig("weihs_universal_plot_v01.png")
    print("Saved: weihs_universal_plot_v01.png")

if __name__ == "__main__":
    main()