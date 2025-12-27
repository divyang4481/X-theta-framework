#!/usr/bin/env python3
"""
weihs_universal_pipeline_v07.py

THE UNIFIED DISCOVERY PIPELINE (Wide-Net Sync Fix)
================================================================================
CHANGES FROM v06:
  1. WIDE-NET SYNC: Calibration now searches a massive +/- 100,000 ns window 
     (100 microseconds) to find the optical delay, fixing the "S=0" issue caused 
     by delays larger than the previous 2us limit.
  2. VISUAL FEEDBACK: Prints "Peak Height" during scanning so you can see 
     if it finds real coincidences (e.g. Peak > 50).
  3. ADAPTIVE WINDOW: Once calibrated, it narrows the window for speed.

USAGE:
   python weihs_universal_pipeline_v07.py --alice-zip "A.zip" --bob-zip "B.zip"
================================================================================
"""

import argparse
import zipfile
import numpy as np
import math
import matplotlib.pyplot as plt
import os
import sys
from dataclasses import dataclass, field
from collections import deque
from typing import List, Tuple

# === CONFIGURATION ===
RESONANCE_KAPPA = 0.8       
PERMUTATION_ROUNDS = 10000  
BIN_COUNT = 6               
CALIB_WINDOW_NS = 100000.0  # Huge window for first lock (+/- 100us)
MINING_WINDOW_NS = 2000.0   # Narrow window for speed after lock
CACHE_FILE = "weihs_mined_data_v7.npz"

# ==========================================
# PART 1: WIDE-NET CALIBRATION
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

def check_file_for_S(za, zb, target):
    """
    Helper: Tries to calculate max S using WIDE search.
    """
    tA, cA = get_raw_buffers(za, target)
    tB, cB = get_raw_buffers(zb, target)
    if tA is None or len(tA) < 2000: return 0, None, 0
    
    # Use decent chunk for sync
    limit = min(len(tA), len(tB), 20000)
    tA, cA, tB, cB = tA[:limit], cA[:limit], tB[:limit], cB[:limit]
    
    # --- WIDE SYNC SEARCH ---
    # We only need to find the shift ONCE per file, independent of bitmask
    diffs = (tB[:5000] * 1e9)[:,None] - (tA[:5000] * 1e9)[None,:]
    
    # Coarse Histogram (Large bins)
    hist, edges = np.histogram(diffs, bins=5000, range=(-CALIB_WINDOW_NS, CALIB_WINDOW_NS))
    peak_bin = np.argmax(hist)
    peak_val_ns = (edges[peak_bin] + edges[peak_bin+1])/2.0
    peak_height = hist[peak_bin]
    
    if peak_height < 10: 
        return 0, None, peak_height # No signal found
        
    # Refine Peak (Zoom in +/- 100ns around coarse peak)
    valid_diffs = diffs[np.abs(diffs - peak_val_ns) < 100]
    if len(valid_diffs) == 0: return 0, None, peak_height
    
    # Fine Histogram
    hist_fine, edges_fine = np.histogram(valid_diffs, bins=200, range=(peak_val_ns-50, peak_val_ns+50))
    peak_idx_fine = np.argmax(hist_fine)
    fine_shift = (edges_fine[peak_idx_fine] + edges_fine[peak_idx_fine+1])/2.0 * 1e-9 # Seconds

    local_best_S = 0.0
    local_best_cfg = None
    
    # Pre-calculate shifted tB for speed
    dt_matrix = diffs - (fine_shift * 1e9)
    pairs = np.argwhere(np.abs(dt_matrix) < 2.5) # 2.5ns coincidence
    
    if len(pairs) < 50: return 0, None, peak_height
    idxA, idxB = pairs[:, 1], pairs[:, 0]

    for s_bit in [0, 1, 2]:
        for d_bit in [0, 1, 2]:
            if s_bit == d_bit: continue
            
            sA = ((cA >> s_bit) & 1).astype(int)
            oA = (2 * ((cA >> d_bit) & 1).astype(int)) - 1
            sB = ((cB >> s_bit) & 1).astype(int)
            oB = (2 * ((cB >> d_bit) & 1).astype(int)) - 1
            
            # Entropy Check
            if np.std(oA) < 0.2 or np.std(oB) < 0.2: continue 
            
            # Calculate S on matched pairs
            v_sA = sA[idxA]; v_sB = sB[idxB]
            v_prod = oA[idxA] * oB[idxB]
            
            counts = np.bincount(v_sA*2 + v_sB, minlength=4)
            sums = np.bincount(v_sA*2 + v_sB, weights=v_prod, minlength=4)
            
            if np.min(counts) < 5: continue
            
            E = sums / counts
            S = E[0] + E[1] + E[2] - E[3]
            
            if abs(S) > local_best_S:
                local_best_S = abs(S)
                flip = -1 if S < 0 else 1
                local_best_cfg = (s_bit, d_bit, flip)
                
    return local_best_S, local_best_cfg, peak_height

def calibrate_bitmask(za, zb):
    print(f"[CALIBRATION] Wide-Net Mode (+/- {int(CALIB_WINDOW_NS)}ns)...")
    
    a_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in za.namelist() if "_v.dat" in n.lower() }
    b_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in zb.namelist() if "_v.dat" in n.lower() }
    common = sorted(list(a_files.intersection(b_files)))
    
    priority = [c for c in common if "sine" in c] + [c for c in common if "longdist" in c]
    candidates = (priority + common)[:50] 
    
    global_best_S = 0.0
    global_best_cfg = (1, 0, 1) 
    
    for target in candidates:
        S, cfg, height = check_file_for_S(za, zb, target)
        
        # Log progress for user sanity
        status = f"S={S:.4f}" if S > 0 else "No S"
        if height > 50:
            print(f"  {target[:20]}... Peak={height} -> {status}")
        
        if S > global_best_S:
            global_best_S = S
            global_best_cfg = cfg
            if S > 2.0: 
                print(f"  VIOLATION FOUND! Locking configuration.")
                return cfg
            
    print(f"  Tournament Finished. Winner S={global_best_S:.4f}. Config: {global_best_cfg}")
    return global_best_cfg

# ==========================================
# PART 2: MINING ENGINE
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
                self.theta += self.theta_kappa * self._get_orientation(h)

    def capture(self, a, b, prod):
        k = a * 2 + b
        self.collected_events.append((self.theta, k, prod))

    def _get_orientation(self, states):
        area = 0.0
        for (x1, y1), (x2, y2) in zip(states[:-1], states[1:]):
            area += (float(x1)*float(y2) - float(x2)*float(y1))
        return 1 if area > 0.1 else (-1 if area < -0.1 else 0)

def find_shift_and_sync(tA, tB):
    limit = min(5000, len(tA), len(tB))
    if limit < 100: return None, 0
    tA_s = tA[:limit] * 1e9; tB_s = tB[:limit] * 1e9
    diffs = tB_s[:, None] - tA_s[None, :]
    
    # Phase 1: Wide Search (reused logic for robustness)
    # We use a narrower window for mining than calibration to save time, 
    # assuming calibration found typical delay range? 
    # Actually, delay might vary per file. Let's use Medium Window for mining.
    
    valid = diffs[np.abs(diffs) < 20000.0] # +/- 20us for mining
    if len(valid) == 0: return None, 0
    
    hist, edges = np.histogram(valid, bins=2000)
    peak_idx = np.argmax(hist)
    peak_val = hist[peak_idx]
    
    if peak_val < np.mean(hist) * 3: return None, 0 
    
    shift = (edges[peak_idx] + edges[peak_idx+1])/2.0 * 1e-9
    return shift, peak_val

def mine_data(alice_zip, bob_zip):
    print("\n[MINING] Starting Universal Miner v07...")
    engine = ContinuousTopology(theta_kappa=RESONANCE_KAPPA)
    
    with zipfile.ZipFile(alice_zip, "r") as za, zipfile.ZipFile(bob_zip, "r") as zb:
        s_bit, d_bit, flip = calibrate_bitmask(za, zb)
        
        a_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in za.namelist() if "_v.dat" in n.lower() }
        b_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in zb.namelist() if "_v.dat" in n.lower() }
        stems = sorted(list(a_files.intersection(b_files)))
        print(f"  Found {len(stems)} file pairs.")
        
        processed_count = 0
        for i, stem in enumerate(stems):
            tA, cA = get_raw_buffers(za, stem)
            tB, cB = get_raw_buffers(zb, stem)
            if tA is None or len(tA) < 500: continue
            
            sA = ((cA >> s_bit) & 1).astype(int)
            oA = ((2 * ((cA >> d_bit) & 1).astype(int)) - 1) * flip
            sB = ((cB >> s_bit) & 1).astype(int)
            oB = ((2 * ((cB >> d_bit) & 1).astype(int)) - 1) * flip 
            
            shift, snr = find_shift_and_sync(tA, tB)
            if shift is None: continue
            
            tB_s = tB + shift
            i_a, i_b = 0, 0
            nA, nB = len(tA), len(tB)
            cmap = {}
            
            while i_a < nA and i_b < nB:
                d = tB_s[i_b] - tA[i_a]
                if d < -1.5e-9: i_b += 1
                elif d > 1.5e-9: i_a += 1
                else:
                    cmap[i_a] = i_b
                    i_a += 1; i_b += 1
            
            evs = [(tA[k], 0, k) for k in range(nA)] + [(tB_s[k], 1, k) for k in range(nB)]
            evs.sort(key=lambda x: x[0])
            
            for _, src, idx in evs:
                if src == 0:
                    engine.update('A', sA[idx])
                    if idx in cmap:
                        prod = int(oA[idx]) * int(oB[cmap[idx]])
                        engine.capture(sA[idx], sB[cmap[idx]], prod)
                else:
                    engine.update('B', sB[idx])
            
            processed_count += 1
            if processed_count % 20 == 0:
                print(f"  Processed {processed_count} files... (N={len(engine.collected_events)})")

    arr = np.array(engine.collected_events)
    np.savez(CACHE_FILE, phases=arr[:,0], settings=arr[:,1], products=arr[:,2])
    print(f"  Mining Complete. Saved {len(arr)} events to {CACHE_FILE}")

# ==========================================
# PART 3: ANALYSIS & PLOTTING
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
    print("   WEIHS UNIVERSAL PIPELINE v07")
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
    print(f"  Global S Check: {g_S:.4f}")
    
    if abs(g_S) < 1.0:
        print("CRITICAL: S is low. Aborting.")
        return

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
    plt.title(f"Phase-Dependent Violation\n$N={len(phases)}, \Delta\chi^2={res['delta']:.1f}, p={p_val:.5f}$")
    plt.legend(loc='lower right'); plt.grid(True, alpha=0.3)
    plt.savefig("weihs_final_result.png")
    print("Saved: weihs_final_result.png")

if __name__ == "__main__":
    main()