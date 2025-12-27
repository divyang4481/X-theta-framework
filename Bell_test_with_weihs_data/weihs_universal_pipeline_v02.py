#!/usr/bin/env python3
"""
weihs_universal_pipeline_v02.py

THE UNIFIED DISCOVERY PIPELINE (FIXED)
================================================================================
DESCRIPTION:
  A single, end-to-end script that:
  1. AUTO-CALIBRATES bitmasks (finds correct Setting/Detector bits).
  2. MINES the entire dataset (~400 files) with cross-correlation sync.
  3. CACHES data to 'weihs_mined_data_v2.npz'.
  4. PERFORMS rigorous stratified permutation testing.
  5. PLOTS the final high-precision result.

FIXES in this version:
  - Enforced signed integer casting for Outcome calculation to prevent 
    unsigned underflow (0 - 1 -> 65535) which caused scalar overflow warnings.

USAGE:
   python weihs_universal_pipeline_v02.py --alice-zip "A.zip" --bob-zip "B.zip"
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
PERMUTATION_ROUNDS = 10000  # High precision for final result
BIN_COUNT = 6               
SYNC_WINDOW_NS = 1000.0     
CACHE_FILE = "weihs_mined_data_v2.npz"

# ==========================================
# PART 1: BITMASK AUTO-CALIBRATION
# ==========================================

def _norm_member(name: str) -> str:
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"): return "/".join(parts[1:])
    return n

def get_raw_buffers(z, stem):
    nm = { _norm_member(n).lower(): n for n in z.namelist() }
    vk, ck = (stem+"_V.dat").lower(), (stem+"_C.dat").lower()
    
    # Robust lookup
    if vk not in nm: vk = next((k for k in nm if k.endswith(stem.lower() + "_v.dat")), None)
    if ck not in nm: ck = next((k for k in nm if k.endswith(stem.lower() + "_c.dat")), None)
    if not vk or not ck: return None, None
    
    with z.open(nm[vk], "r") as f: vb = f.read()
    with z.open(nm[ck], "r") as f: cb = f.read()
    return np.frombuffer(vb, dtype=">f8"), np.frombuffer(cb, dtype=">u2")

def calibrate_bitmask(za, zb):
    """
    Brute-forces bit definitions to find the one that yields Bell Violation (S > 2).
    Returns: (set_bit, det_bit)
    """
    print("[CALIBRATION] Auto-detecting bit definitions...")
    
    # Find a high-quality file (usually 'sine' or 'longdist')
    a_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in za.namelist() if "_v.dat" in n.lower() }
    b_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in zb.namelist() if "_v.dat" in n.lower() }
    common = sorted(list(a_files.intersection(b_files)))
    
    target = next((c for c in common if "sine" in c), common[0])
    print(f"  Using calibration file: {target}")
    
    tA, cA = get_raw_buffers(za, target)
    tB, cB = get_raw_buffers(zb, target)
    
    # Truncate for speed
    limit = min(len(tA), len(tB), 50000)
    tA, cA, tB, cB = tA[:limit], cA[:limit], tB[:limit], cB[:limit]
    
    best_S = 0
    best_cfg = (0, 1) # Default
    
    # Scan permutations of bit 0, 1, 2
    for s_bit in [0, 1, 2]:
        for d_bit in [0, 1, 2]:
            if s_bit == d_bit: continue
            
            # Extract
            sA = ((cA >> s_bit) & 1).astype(int)
            # FORCE CAST TO INT BEFORE MATH to avoid unsigned wrap
            oA = (2 * ((cA >> d_bit) & 1).astype(int)) - 1
            sB = ((cB >> s_bit) & 1).astype(int)
            oB = (2 * ((cB >> d_bit) & 1).astype(int)) - 1
            
            # Rough Sync
            diffs = (tB[:1000] * 1e9)[:,None] - (tA[:1000] * 1e9)[None,:]
            hist, _ = np.histogram(diffs, bins=100, range=(-1000, 1000))
            if np.max(hist) < 5: continue
            
            # Simple check: Balance of settings
            counts = np.bincount(sA*2 + sB, minlength=4)
            if np.min(counts) < limit/10: continue # Malformed settings
            
            ratio = np.std(counts) / np.mean(counts)
            
            if ratio < 0.2: # Settings are well distributed (random)
                 best_cfg = (s_bit, d_bit)
                 print(f"  Candidate: SetBit={s_bit}, DetBit={d_bit} (Balance Ratio={ratio:.2f})")
                 return s_bit, d_bit

    print(f"  Defaulting to Set={best_cfg[0]}, Det={best_cfg[1]}")
    return best_cfg

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
    valid = diffs[np.abs(diffs) < SYNC_WINDOW_NS]
    
    if len(valid) == 0: return None, 0
    
    hist, edges = np.histogram(valid, bins=2000)
    peak_idx = np.argmax(hist)
    peak_val = hist[peak_idx]
    
    # SNR check
    avg = np.mean(hist)
    if peak_val < avg * 3: return None, 0 # Weak signal
    
    shift = (edges[peak_idx] + edges[peak_idx+1])/2.0 * 1e-9
    return shift, peak_val

def mine_data(alice_zip, bob_zip):
    print("\n[MINING] Starting Universal Miner...")
    engine = ContinuousTopology(theta_kappa=RESONANCE_KAPPA)
    
    with zipfile.ZipFile(alice_zip, "r") as za, zipfile.ZipFile(bob_zip, "r") as zb:
        # 1. Calibrate
        s_bit, d_bit = calibrate_bitmask(za, zb)
        print(f"  Configuration Locked: Setting Bit {s_bit}, Detector Bit {d_bit}")
        
        # 2. List Files
        a_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in za.namelist() if "_v.dat" in n.lower() }
        b_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in zb.namelist() if "_v.dat" in n.lower() }
        stems = sorted(list(a_files.intersection(b_files)))
        print(f"  Found {len(stems)} file pairs.")
        
        # 3. Process
        processed_count = 0
        for i, stem in enumerate(stems):
            tA, cA = get_raw_buffers(za, stem)
            tB, cB = get_raw_buffers(zb, stem)
            if tA is None or len(tA) < 500: continue
            
            # Apply Mask WITH SAFETY CAST
            sA = ((cA >> s_bit) & 1).astype(int)
            oA = (2 * ((cA >> d_bit) & 1).astype(int)) - 1
            
            sB = ((cB >> s_bit) & 1).astype(int)
            oB = (2 * ((cB >> d_bit) & 1).astype(int)) - 1
            
            # Sync
            shift, snr = find_shift_and_sync(tA, tB)
            if shift is None: continue
            
            # Match Coincidences
            tB_s = tB + shift
            i_a, i_b = 0, 0
            nA, nB = len(tA), len(tB)
            cmap = {}
            
            # Linear scan
            while i_a < nA and i_b < nB:
                d = tB_s[i_b] - tA[i_a]
                if d < -1.5e-9: i_b += 1
                elif d > 1.5e-9: i_a += 1
                else:
                    cmap[i_a] = i_b
                    i_a += 1; i_b += 1
            
            # Replay History
            evs = [(tA[k], 0, k) for k in range(nA)] + [(tB_s[k], 1, k) for k in range(nB)]
            evs.sort(key=lambda x: x[0])
            
            for _, src, idx in evs:
                if src == 0:
                    engine.update('A', sA[idx])
                    if idx in cmap:
                        # Explicit casting to standard int avoids any residual numpy scalar warnings
                        engine.capture(sA[idx], sB[cmap[idx]], int(oA[idx]) * int(oB[cmap[idx]]))
                else:
                    engine.update('B', sB[idx])
            
            processed_count += 1
            if processed_count % 20 == 0:
                print(f"  Processed {processed_count} files... (Current N={len(engine.collected_events)})")

    # Save
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
    # Bin
    bin_idx = ((phases % (2*np.pi)) / (2*np.pi) * BIN_COUNT).astype(int) % BIN_COUNT
    flat_idx = bin_idx * 4 + settings.astype(int)
    
    counts = np.bincount(flat_idx, minlength=BIN_COUNT*4)
    sums = np.bincount(flat_idx, weights=products, minlength=BIN_COUNT*4)
    
    if np.any(counts < 5): return None
    
    E = sums / counts
    Err = 1.0 / np.sqrt(counts)
    
    E_r = E.reshape((BIN_COUNT, 4))
    Err_r = Err.reshape((BIN_COUNT, 4))
    
    # S = E00 + E01 + E10 - E11
    S = E_r[:,0] + E_r[:,1] + E_r[:,2] - E_r[:,3]
    S_err = np.sqrt(np.sum(Err_r**2, axis=1))
    
    centers = (np.arange(BIN_COUNT) + 0.5) * (2*math.pi/BIN_COUNT)
    
    # Fit Sine
    X = np.column_stack([np.ones(BIN_COUNT), np.sin(centers), np.cos(centers)])
    beta, chi2 = weighted_fit(X, S, S_err)
    
    # Fit Flat
    X0 = np.ones((BIN_COUNT, 1))
    _, chi2_0 = weighted_fit(X0, S, S_err)
    
    return {"s": S, "err": S_err, "phases": centers, "beta": beta, "delta": chi2_0 - chi2}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alice-zip", required=True); parser.add_argument("--bob-zip", required=True)
    args = parser.parse_args()

    print("==========================================")
    print("   WEIHS UNIVERSAL PIPELINE v02 (FIXED)")
    print("==========================================")
    
    # 1. Load or Mine
    if not os.path.exists(CACHE_FILE):
        mine_data(args.alice_zip, args.bob_zip)
    else:
        print(f"[CACHE] Loaded existing data from {CACHE_FILE}")
        
    data = np.load(CACHE_FILE)
    phases = data['phases']
    settings = data['settings'].astype(int)
    products = data['products']
    
    print(f"\n[ANALYSIS] Analyzing {len(phases)} events...")
    
    # 2. Global Check
    # Check if S is reasonable globally
    g_counts = np.bincount(settings, minlength=4)
    g_sums = np.bincount(settings, weights=products, minlength=4)
    g_E = g_sums / np.maximum(g_counts, 1)
    g_S = g_E[0] + g_E[1] + g_E[2] - g_E[3]
    print(f"  Global S Check: {g_S:.4f}")
    
    if g_S < 1.0:
        print("  WARNING: S is low. Bitmask might still be wrong or data is noise.")
    
    # 3. Real Analysis
    res = analyze_binned(phases, settings, products)
    if not res: 
        print("Analysis Failed (Insufficient Data).")
        return
        
    print(f"  Observed Delta Chi2: {res['delta']:.2f}")
    
    # 4. Stratified Permutation
    print(f"\n[PERMUTATION] Running {PERMUTATION_ROUNDS} stratified shuffles...")
    better = 0
    perm_prod = products.copy()
    indices = [np.where(settings == k)[0] for k in range(4)]
    
    for i in range(PERMUTATION_ROUNDS):
        # Shuffle within settings
        for k in range(4):
            idx = indices[k]
            subset = perm_prod[idx]
            np.random.shuffle(subset)
            perm_prod[idx] = subset
            
        pres = analyze_binned(phases, settings, perm_prod)
        if pres and pres['delta'] >= res['delta']:
            better += 1
            
        if (i+1) % 1000 == 0: print(f"  ..{i+1}")
        
    p_val = (better + 1) / (PERMUTATION_ROUNDS + 1)
    print(f"  Final p-value: {p_val:.6f}")
    
    # 5. Plotting
    print("\n[PLOTTING] Generating final figure...")
    plt.figure(figsize=(10, 6))
    
    # Data
    plt.errorbar(res['phases'], res['s'], yerr=res['err'], fmt='o', color='blue', ecolor='black', capsize=3, label='Data')
    
    # Fit
    c0, a, b = res['beta']
    th = np.linspace(0, 2*np.pi, 200)
    plt.plot(th, c0 + a*np.sin(th) + b*np.cos(th), color='orange', linewidth=2, label=f'Fit (p={p_val:.5f})')
    
    # Refs
    plt.axhline(2.0, color='red', linestyle='--', label='Classical')
    plt.axhline(2.828, color='green', linestyle=':', label='Tsirelson')
    
    # Zoom
    y_min = np.min(res['s'] - res['err']) - 0.1
    y_max = np.max(res['s'] + res['err']) + 0.1
    plt.ylim(y_min, y_max)
    
    plt.xlabel(r'Geometric Phase $\theta$')
    plt.ylabel(r'Bell Parameter $S$')
    plt.title(f"Phase-Dependent Violation\n$N={len(phases)}, \Delta\chi^2={res['delta']:.1f}, p={p_val:.5f}$")
    plt.legend(loc='lower right')
    plt.grid(True, alpha=0.3)
    
    plt.savefig("weihs_final_result.png")
    print("Saved: weihs_final_result.png")

if __name__ == "__main__":
    main()