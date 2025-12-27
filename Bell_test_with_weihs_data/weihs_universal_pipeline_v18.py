#!/usr/bin/env python3
"""
weihs_universal_pipeline_v18.py

THE UNIFIED DISCOVERY PIPELINE ("Strict Quality Control")
================================================================================
RESPONSES TO REVIEWER:
  1. THETA RESET: The topology engine now resets theta=0 at the start of 
     every file. This prevents artificial correlations from file ordering.
  2. PER-FILE LOCK: We calculate S for EACH file individually. We strictly 
     DISCARD any file where S < 1.5. This ensures we are only analyzing 
     valid Bell tests, not noise.
  3. PLUMBING FIXES: 
     - Sync sign corrected (tB - shift).
     - Flip applied only to final product.
     - Nearest-neighbor pairing used to handle jitter.

USAGE:
   python weihs_universal_pipeline_v18.py --alice-zip "A.zip" --bob-zip "B.zip"
================================================================================
"""

import argparse
import zipfile
import numpy as np
import matplotlib.pyplot as plt
import os
from dataclasses import dataclass, field
from collections import deque, Counter
from typing import List, Tuple

# === CONFIGURATION ===
RESONANCE_KAPPA = 0.8       
PERMUTATION_ROUNDS = 10000  
BIN_COUNT = 6               
WIDE_WINDOW_NS = 100000.0   # +/- 100us for coarse sync
COINCIDENCE_WINDOW = 3.5e-9 # 3.5ns window
S_THRESHOLD = 1.5           # Discard files with S below this (Noise filter)
CACHE_FILE = "weihs_mined_data_v18.npz"

# ==========================================
# PART 1: CORE UTILITIES
# ==========================================

def _norm_member(name: str) -> str:
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"): return "/".join(parts[1:])
    return n

def get_raw_buffers(z, stem):
    nm = { _norm_member(n).lower(): n for n in z.namelist() }
    vk, ck = (stem+"_V.dat").lower(), (stem+"_C.dat").lower()
    if vk not in nm: vk = next((k for k in nm if k.endswith(stem.lower() + "_v.dat")), None)
    if ck not in nm: ck = next((k for k in nm if k.endswith(stem.lower() + "_c.dat")), None)
    if not vk or not ck: return None, None
    with z.open(nm[vk], "r") as f: vb = f.read()
    with z.open(nm[ck], "r") as f: cb = f.read()
    return np.frombuffer(vb, dtype=">f8"), np.frombuffer(cb, dtype=">u2")

def match_nearest_neighbor(tA, tB, window):
    """
    Robust pairing: finds nearest photon within window.
    """
    if len(tB) == 0: return np.array([]), np.array([])
    idx = np.searchsorted(tB, tA)
    idx = np.clip(idx, 0, len(tB)-1)
    
    dt_right = tB[idx] - tA
    idx_left = np.maximum(idx - 1, 0)
    dt_left = tB[idx_left] - tA
    
    use_left = np.abs(dt_left) < np.abs(dt_right)
    best_idx = np.where(use_left, idx_left, idx)
    best_dt = tB[best_idx] - tA
    
    mask = np.abs(best_dt) < window
    return np.where(mask)[0], best_idx[mask]

def find_shift_and_sync(tA, tB):
    # Use larger slice for robust sync finding
    limit = min(len(tA), len(tB), 20000)
    if limit < 100: return None
    
    tA_s = tA[:limit] * 1e9
    tB_s = tB[:limit] * 1e9
    diffs = tB_s[:, None] - tA_s[None, :]
    
    # Wide Search
    valid = diffs[np.abs(diffs) < WIDE_WINDOW_NS] 
    if len(valid) == 0: return None
    
    # Coarse Hist
    hist, edges = np.histogram(valid, bins=5000)
    peak_bin = np.argmax(hist)
    
    if hist[peak_bin] < 10: return None 
    
    # Fine Search
    peak_coarse = (edges[peak_bin] + edges[peak_bin+1])/2.0
    valid_fine = valid[np.abs(valid - peak_coarse) < 100] 
    if len(valid_fine) == 0: return None
    
    hist_f, edges_f = np.histogram(valid_fine, bins=200) # 1ns bins
    peak_f = np.argmax(hist_f)
    shift = (edges_f[peak_f] + edges_f[peak_f+1])/2.0 * 1e-9
    
    # SNR Metric
    bg_mean = np.mean(hist)
    snr = hist[peak_bin] / (bg_mean + 0.1)
    
    return shift, snr

# ==========================================
# PART 2: CALIBRATION
# ==========================================

def calibrate_best_file(za, zb, stems):
    """
    Scans for the single best file to determine bitmask.
    """
    print(f"  [CALIBRATION] Scanning {len(stems)} files for Golden Key...")
    
    global_max_S = 0.0
    global_best_cfg = (0, 1, 1) # Default
    
    # Only test the two standard configs
    configs = [(0, 1), (1, 0)]
    
    for stem in stems:
        tA, cA = get_raw_buffers(za, stem)
        tB, cB = get_raw_buffers(zb, stem)
        if tA is None or len(tA) < 5000: continue
        
        res = find_shift_and_sync(tA, tB)
        if res is None: continue
        shift, snr = res
        
        if snr < 5.0: continue # Skip noisy files for calibration
        
        tB_shifted = tB - shift
        idxA, idxB = match_nearest_neighbor(tA, tB_shifted, COINCIDENCE_WINDOW)
        if len(idxA) < 100: continue
        
        for (s, d) in configs:
            sA = ((cA[idxA] >> s) & 1).astype(int)
            oA = (2 * ((cA[idxA] >> d) & 1).astype(int)) - 1
            sB = ((cB[idxB] >> s) & 1).astype(int)
            oB = (2 * ((cB[idxB] >> d) & 1).astype(int)) - 1
            
            if np.std(oA) < 0.1: continue
            
            # Local S check
            prod = oA * oB
            k = sA * 2 + sB
            counts = np.bincount(k, minlength=4)
            sums = np.bincount(k, weights=prod, minlength=4)
            if np.min(counts) < 10: continue
            
            S = (sums/counts)[0] + (sums/counts)[1] + (sums/counts)[2] - (sums/counts)[3]
            
            if abs(S) > global_max_S:
                global_max_S = abs(S)
                flip = -1 if S < 0 else 1
                global_best_cfg = (s, d, flip)
                print(f"    New Best: {stem} -> S={global_max_S:.4f} (Cfg={global_best_cfg})")

    print(f"  WINNER: S={global_max_S:.4f} with Config {global_best_cfg}")
    return global_best_cfg

# ==========================================
# PART 3: MINING & TOPOLOGY
# ==========================================

@dataclass
class ContinuousTopology:
    theta_kappa: float
    a: int=0; b: int=0; theta: float=0.0
    history: deque = field(default_factory=lambda: deque(maxlen=5))
    collected_events: List[Tuple[float, int, int]] = field(default_factory=list)

    def reset(self):
        """CRITICAL: Resets theta for new run to prevent cross-file artifacts."""
        self.a = 0
        self.b = 0
        self.theta = 0.0
        self.history.clear()

    def update(self, actor, val):
        if actor == 'A': self.a = int(val)
        else: self.b = int(val)
        self.history.append((self.a, self.b))
        if len(self.history) == 5:
            h = list(self.history)
            if h[0] == h[-1] and len(set(h[:-1])) == 4:
                self.theta += self.theta_kappa * self._get_orientation(h)

    def capture(self, a, b, prod):
        k = a * 2 + b
        self.collected_events.append((self.theta, k, prod))

    def _get_orientation(self, states):
        area = 0.0
        for (x1, y1), (x2, y2) in zip(states[:-1], states[1:]):
            area += (float(x1)*float(y2) - float(x2)*float(y1))
        return 1 if area > 0.1 else (-1 if area < -0.1 else 0)

def mine_data(alice_zip, bob_zip):
    print("\n[MINING] Starting Universal Miner v18...")
    engine = ContinuousTopology(theta_kappa=RESONANCE_KAPPA)
    
    with zipfile.ZipFile(alice_zip, "r") as za, zipfile.ZipFile(bob_zip, "r") as zb:
        a_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in za.namelist() if "_v.dat" in n.lower() }
        b_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in zb.namelist() if "_v.dat" in n.lower() }
        all_stems = sorted(list(a_files.intersection(b_files)))
        
        stems = [s for s in all_stems if "longdist" in s]
        print(f"  Filtered to {len(stems)} 'longdist' files.")
        
        # CALIBRATE
        s_bit, d_bit, flip = calibrate_best_file(za, zb, stems)
        
        kept_files = 0
        total_events = 0
        
        for i, stem in enumerate(stems):
            tA, cA = get_raw_buffers(za, stem)
            tB, cB = get_raw_buffers(zb, stem)
            if tA is None or len(tA) < 500: continue
            
            # Sync
            res = find_shift_and_sync(tA, tB)
            if res is None: continue
            shift, snr = res
            
            # Extract
            sA = ((cA >> s_bit) & 1).astype(int)
            oA = (2 * ((cA >> d_bit) & 1).astype(int)) - 1
            sB = ((cB >> s_bit) & 1).astype(int)
            oB = (2 * ((cB >> d_bit) & 1).astype(int)) - 1
            
            # Match
            tB_shifted = tB - shift
            idxA, idxB = match_nearest_neighbor(tA, tB_shifted, COINCIDENCE_WINDOW)
            
            if len(idxA) < 50: continue
            
            # PER-FILE VALIDATION (The Reviewer's Demand)
            # Calculate S for THIS file only
            cur_sA = sA[idxA]; cur_sB = sB[idxB]
            cur_prod = oA[idxA] * oB[idxB] * flip
            
            k = cur_sA * 2 + cur_sB
            counts = np.bincount(k, minlength=4)
            sums = np.bincount(k, weights=cur_prod, minlength=4)
            if np.min(counts) < 5: continue
            
            file_S = (sums/counts)[0] + (sums/counts)[1] + (sums/counts)[2] - (sums/counts)[3]
            
            # REJECT JUNK
            if file_S < S_THRESHOLD:
                # print(f"    [Reject] {stem}: S={file_S:.2f} (Too low)")
                continue
            
            # ACCEPT & ACCUMULATE
            kept_files += 1
            engine.reset() # CRITICAL: Reset theta for new run
            
            # We iterate matched pairs to simulate topological evolution
            # Note: Ideally we'd run the whole stream, but for now we run the matched stream
            for k in range(len(idxA)):
                engine.update('A', cur_sA[k])
                engine.update('B', cur_sB[k])
                engine.capture(cur_sA[k], cur_sB[k], cur_prod[k])
                
            total_events += len(idxA)
            if kept_files % 5 == 0:
                print(f"  Accepted {kept_files} files (Last S={file_S:.2f}). Total Events: {total_events}")

    if kept_files == 0:
        print("CRITICAL: No files passed the S > 1.5 check. Pipeline halted.")
        return

    arr = np.array(engine.collected_events)
    np.savez(CACHE_FILE, phases=arr[:,0], settings=arr[:,1], products=arr[:,2])
    print(f"  Mining Complete. {len(arr)} events saved from {kept_files} valid files.")

# ==========================================
# PART 4: ANALYSIS
# ==========================================

def weighted_fit(X, y, sigma):
    w = 1.0 / (np.maximum(sigma, 1e-12) ** 2)
    sw = np.sqrt(w)
    Xw = X * sw[:, None]; yw = y * sw
    try:
        XtX = Xw.T @ Xw
        beta = np.linalg.solve(XtX, Xw.T @ yw)
        yhat = X @ beta
        chi2 = np.sum(((y - yhat) / sigma) ** 2)
        return beta, chi2
    except: return None, 0

def analyze_binned(phases, settings, products):
    bin_idx = ((phases % (2*np.pi)) / (2*np.pi) * BIN_COUNT).astype(int) % BIN_COUNT
    flat_idx = bin_idx * 4 + settings.astype(int)
    counts = np.bincount(flat_idx, minlength=BIN_COUNT*4)
    sums = np.bincount(flat_idx, weights=products, minlength=BIN_COUNT*4)
    
    if np.any(counts < 5): return None
    E = sums / counts; Err = 1.0 / np.sqrt(counts)
    
    E_r = E.reshape((BIN_COUNT, 4)); Err_r = Err.reshape((BIN_COUNT, 4))
    S = E_r[:,0] + E_r[:,1] + E_r[:,2] - E_r[:,3]
    S_err = np.sqrt(np.sum(Err_r**2, axis=1))
    
    centers = (np.arange(BIN_COUNT) + 0.5) * (2*math.pi/BIN_COUNT)
    X = np.column_stack([np.ones(BIN_COUNT), np.sin(centers), np.cos(centers)])
    beta, chi2 = weighted_fit(X, S, S_err)
    X0 = np.ones((BIN_COUNT, 1))
    _, chi2_0 = weighted_fit(X0, S, S_err)
    
    return {"s": S, "err": S_err, "phases": centers, "beta": beta, "delta": chi2_0 - chi2}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alice-zip", required=True); parser.add_argument("--bob-zip", required=True)
    args = parser.parse_args()

    print("==========================================")
    print("   WEIHS UNIVERSAL PIPELINE v18")
    print("==========================================")
    
    if not os.path.exists(CACHE_FILE):
        mine_data(args.alice_zip, args.bob_zip)
    else:
        print(f"[CACHE] Loaded {CACHE_FILE}")
        
    data = np.load(CACHE_FILE)
    phases, settings, products = data['phases'], data['settings'].astype(int), data['products']
    print(f"\n[ANALYSIS] Analyzing {len(phases)} events...")
    
    g_counts = np.bincount(settings, minlength=4)
    g_sums = np.bincount(settings, weights=products, minlength=4)
    g_E = g_sums / np.maximum(g_counts, 1)
    g_S = g_E[0] + g_E[1] + g_E[2] - g_E[3]
    print(f"  Global S (Raw): {g_S:.4f}")
    
    if g_S < 1.5:
        print("\nWARNING: Global S is below 2.0. This may indicate poor quality data despite filtering.")
    
    res = analyze_binned(phases, settings, products)
    print(f"  Observed Delta Chi2: {res['delta']:.2f}")
    
    print(f"\n[PERMUTATION] Running {PERMUTATION_ROUNDS} stratified shuffles...")
    better = 0; perm_prod = products.copy(); indices = [np.where(settings == k)[0] for k in range(4)]
    for i in range(PERMUTATION_ROUNDS):
        for k in range(4):
            idx = indices[k]; subset = perm_prod[idx]
            np.random.shuffle(subset); perm_prod[idx] = subset
        pres = analyze_binned(phases, settings, perm_prod)
        if pres and pres['delta'] >= res['delta']: better += 1
        if (i+1) % 1000 == 0: print(f"  ..{i+1}")
        
    p_val = (better + 1) / (PERMUTATION_ROUNDS + 1)
    print(f"  Final p-value: {p_val:.6f}")
    
    plt.figure(figsize=(10, 6))
    plt.errorbar(res['phases'], res['s'], yerr=res['err'], fmt='o', color='blue', ecolor='black', capsize=3, label='Data')
    c0, a, b = res['beta']
    th = np.linspace(0, 2*np.pi, 200)
    plt.plot(th, c0 + a*np.sin(th) + b*np.cos(th), color='orange', linewidth=2, label=f'Fit (p={p_val:.5f})')
    plt.axhline(2.0, color='red', linestyle='--', label='Classical')
    plt.axhline(2.828, color='green', linestyle=':', label='Tsirelson')
    y_min = np.min(res['s'] - res['err']) - 0.1; y_max = np.max(res['s'] + res['err']) + 0.1
    plt.ylim(y_min, y_max)
    plt.title(f"Phase-Dependent Violation\n$N={len(phases)}, S_{{global}}={g_S:.2f}, p={p_val:.5f}$")
    plt.legend(loc='lower right'); plt.grid(True, alpha=0.3)
    plt.savefig("weihs_final_result.png")
    print("Saved: weihs_final_result.png")

if __name__ == "__main__":
    main()