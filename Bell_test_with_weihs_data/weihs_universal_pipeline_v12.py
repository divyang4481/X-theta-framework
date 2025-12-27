#!/usr/bin/env python3
"""
weihs_universal_pipeline_v12.py

THE UNIFIED DISCOVERY PIPELINE (Ensemble Calibration)
================================================================================
STRATEGY CHANGE:
  - Previous versions failed because single files ('longdist0') have too few 
    coincidences (~30) to statistically prove S > 2.0 individually.
  - v12 assumes the config from your successful v19 run (Set=0, Det=1).
  - It mines ALL data first, then checks S on the AGGREGATE dataset.
  - This relies on the "Law of Large Numbers" to reveal the Bell violation 
    that is hidden in the noise of single files.

USAGE:
   python weihs_universal_pipeline_v12.py --alice-zip "A.zip" --bob-zip "B.zip"
================================================================================
"""

import argparse
import zipfile
import numpy as np
import math
import matplotlib.pyplot as plt
import os
from dataclasses import dataclass, field
from collections import deque
from typing import List, Tuple

# === CONFIGURATION ===
RESONANCE_KAPPA = 0.8       
PERMUTATION_ROUNDS = 10000  
BIN_COUNT = 6               
WIDE_WINDOW_NS = 100000.0   # +/- 100us for robust sync
CACHE_FILE = "weihs_mined_data_v12.npz"

# === CORE ENGINE ===

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
    # Load enough data to find a peak
    limit = min(len(tA), len(tB), 50000)
    if limit < 100: return None
    
    tA_s = tA[:limit] * 1e9; tB_s = tB[:limit] * 1e9
    diffs = tB_s[:, None] - tA_s[None, :]
    
    # Wide Search
    valid = diffs[np.abs(diffs) < WIDE_WINDOW_NS] 
    if len(valid) == 0: return None
    
    # Coarse Hist
    hist, edges = np.histogram(valid, bins=5000)
    peak_bin = np.argmax(hist)
    peak_coarse = (edges[peak_bin] + edges[peak_bin+1])/2.0
    
    if hist[peak_bin] < 10: return None # No sync
    
    # Fine Search
    valid_fine = valid[np.abs(valid - peak_coarse) < 50]
    if len(valid_fine) == 0: return None
    
    hist_f, edges_f = np.histogram(valid_fine, bins=100)
    peak_f = np.argmax(hist_f)
    shift = (edges_f[peak_f] + edges_f[peak_f+1])/2.0 * 1e-9
    
    return shift

def mine_data(alice_zip, bob_zip):
    print("\n[MINING] Starting Ensemble Miner v12...")
    engine = ContinuousTopology(theta_kappa=RESONANCE_KAPPA)
    
    # === HYPOTHESIS: Set=0, Det=1 (From v19) ===
    s_bit, d_bit = 0, 1
    print(f"  Using Config: Set={s_bit}, Det={d_bit} (Ensemble Verification)")
    
    with zipfile.ZipFile(alice_zip, "r") as za, zipfile.ZipFile(bob_zip, "r") as zb:
        a_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in za.namelist() if "_v.dat" in n.lower() }
        b_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in zb.namelist() if "_v.dat" in n.lower() }
        stems = sorted(list(a_files.intersection(b_files)))
        
        priority = [c for c in stems if "longdist" in c]
        stems = priority + [c for c in stems if c not in priority]
        
        print(f"  Found {len(stems)} file pairs. Mining ALL...")
        
        processed = 0
        for i, stem in enumerate(stems):
            tA, cA = get_raw_buffers(za, stem)
            tB, cB = get_raw_buffers(zb, stem)
            if tA is None or len(tA) < 500: continue
            
            # Apply Mask
            sA = ((cA >> s_bit) & 1).astype(int)
            oA = (2 * ((cA >> d_bit) & 1).astype(int)) - 1
            sB = ((cB >> s_bit) & 1).astype(int)
            oB = (2 * ((cB >> d_bit) & 1).astype(int)) - 1
            
            # Sync
            shift = find_shift_and_sync(tA, tB)
            if shift is None: continue
            
            tB_s = tB + shift
            i_a, i_b = 0, 0
            nA, nB = len(tA), len(tB)
            cmap = {}
            
            # 2.5ns window
            while i_a < nA and i_b < nB:
                d = tB_s[i_b] - tA[i_a]
                if d < -2.5e-9: i_b += 1
                elif d > 2.5e-9: i_a += 1
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
            
            processed += 1
            if processed % 20 == 0:
                print(f"  Processed {processed} files... (Events: {len(engine.collected_events)})")

    arr = np.array(engine.collected_events)
    np.savez(CACHE_FILE, phases=arr[:,0], settings=arr[:,1], products=arr[:,2])
    print(f"  Mining Complete. {len(arr)} total events saved.")

# === ANALYSIS ===

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
    print("   WEIHS UNIVERSAL PIPELINE v12")
    print("==========================================")
    
    if not os.path.exists(CACHE_FILE):
        mine_data(args.alice_zip, args.bob_zip)
    else:
        print(f"[CACHE] Loaded {CACHE_FILE}")
        
    data = np.load(CACHE_FILE)
    phases, settings, products = data['phases'], data['settings'].astype(int), data['products']
    print(f"\n[ANALYSIS] Analyzing {len(phases)} events...")
    
    # GLOBAL CHECKS
    g_counts = np.bincount(settings, minlength=4)
    g_sums = np.bincount(settings, weights=products, minlength=4)
    g_E = g_sums / np.maximum(g_counts, 1)
    g_S = g_E[0] + g_E[1] + g_E[2] - g_E[3]
    print(f"  Global S (Raw): {g_S:.4f}")
    
    # AUTO-FLIP IF NEEDED
    if g_S < -1.8:
        print("  NOTICE: Global S is negative. Flipping outcomes...")
        products = products * -1
        g_S = -g_S
        
    print(f"  Global S (Final): {g_S:.4f}")
    
    if g_S < 1.0:
        print("\nWARNING: Global S is low. The bitmask (Set=0, Det=1) might be wrong,")
        print("or the data is dominated by noise/accidental coincidences.")
        # Proceeding anyway to show the plot, but with a warning.
    
    res = analyze_binned(phases, settings, products)
    if not res: return
    
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