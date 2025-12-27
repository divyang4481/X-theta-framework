#!/usr/bin/env python3
"""
weihs_drift_diagnostic_v27.py

THE "KILLER" CONTROL EXPERIMENT
================================================================================
PURPOSE:
  To determine if the observed signal at low Kappa is real Geometric Phase
  or just Instrumental Drift (Time Non-stationarity).

TESTS:
  1. Phase Coverage: Does theta actually wrap 0->2pi, or is it just a linear ramp?
  2. Correlation: Is theta just a proxy for Time? (Spearman r)
  3. Time-Bin Control: What is S(time)? Does it match S(theta)?
  4. Block-Shuffle: If we preserve time-structure but scramble topology,
     does the signal vanish? (If YES -> Physics is real. If NO -> It's Drift.)

USAGE:
   python weihs_drift_diagnostic_v27.py --alice-zip "..\Bell_data\7185335\Alice.zip" --bob-zip "..\Bell_data\7185335\Bob.zip" --path-contains "timetags/longdist" --kappa 0.01
"""

import argparse
import zipfile
import numpy as np
import math
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from collections import deque
from scipy.stats import spearmanr

# === CONFIG ===
# We default to the "suspicious" kappa found in the scan
DEFAULT_KAPPA = 0.01 
BIN_COUNT = 8
BLOCK_SIZE = 5000 # Events for block shuffle

# ============================================================
# CORE UTILS
# ============================================================

def _norm_member(name): return name.replace("\\", "/").lstrip("/").split("/")[-1]

def load_raw_data(z, stem):
    nm = { _norm_member(n).lower(): n for n in z.namelist() }
    vk = next((k for k in nm if k.endswith("_v.dat") and stem.lower() in k), None)
    ck = next((k for k in nm if k.endswith("_c.dat") and stem.lower() in k), None)
    if not vk or not ck: return None, None
    with z.open(nm[vk]) as f: t = np.frombuffer(f.read(), ">f8")
    with z.open(nm[ck]) as f: c = np.frombuffer(f.read(), ">u2")
    n = min(len(t), len(c))
    return t[:n].copy(), c[:n].copy()

def decode_bits(c, sb, db, invS, invO, flip):
    s = ((c >> sb) & 1).astype(np.int32)
    if invS: s = 1 - s
    o = (2*((c >> db) & 1).astype(np.int32) - 1)
    if invO: o = -o
    return s, o * flip

def quick_sync_match(tA, tB):
    # Fast sync for diagnostics
    limit = min(20000, len(tA))
    diffs = tB[:limit]*1e9 - tA[:limit, None]*1e9
    valid = diffs[np.abs(diffs) < 200000]
    if len(valid) < 100: return None
    hist, edges = np.histogram(valid, bins=5000)
    peak = (edges[np.argmax(hist)] + edges[np.argmax(hist)+1])/2.0
    
    # Fine
    vf = valid[np.abs(valid - peak) < 100]
    if len(vf) < 10: return None
    hf, ef = np.histogram(vf, bins=100)
    shift = (ef[np.argmax(hf)] + ef[np.argmax(hf)+1])/2.0 * 1e-9
    
    # Greedy Match (Fast)
    tB_s = tB - shift
    idx = np.searchsorted(tB_s, tA)
    idx = np.clip(idx, 0, len(tB)-1)
    dt = np.abs(tB_s[idx] - tA)
    mask = dt < 1.5e-9
    return np.where(mask)[0], idx[mask], shift

@dataclass
class ThetaEngine:
    kappa: float
    history: deque = field(default_factory=lambda: deque(maxlen=5))
    theta: float = 0.0
    a: int = 0; b: int = 0
    
    def reset(self):
        self.history.clear()
        self.theta = 0.0
        self.a = 0
        self.b = 0
    
    def update(self, who, val):
        if who==0: self.a=val
        else: self.b=val
        self.history.append((self.a, self.b))
        if len(self.history)==5:
            h = list(self.history)
            area = 0.0
            for i in range(4):
                area += float(h[i][0]*h[i+1][1] - h[i+1][0]*h[i][1])
            ori = 1 if area > 0.1 else (-1 if area < -0.1 else 0)
            self.theta += self.kappa * ori

# ============================================================
# ANALYSIS ENGINES
# ============================================================

def get_stats(x_vals, ks, prods, bins=8):
    # Generic Binner (works for Theta OR Time)
    idx = (x_vals * bins).astype(int) % bins
    flat = idx * 4 + ks
    c = np.bincount(flat, minlength=bins*4).reshape(bins, 4)
    s = np.bincount(flat, weights=prods, minlength=bins*4).reshape(bins, 4)
    
    valid = np.all(c > 10, axis=1)
    if np.sum(valid) < 4: return None
    
    m = s[valid]/c[valid]
    S = m[:,0]+m[:,1]+m[:,2]-m[:,3]
    err = np.sqrt(np.sum((1-m**2)/c[valid], axis=1))
    centers = (np.where(valid)[0]+0.5)/bins
    
    # Fit Sine
    w = 1.0/np.maximum(err, 1e-12)**2; sw = np.sqrt(w)
    th = centers * 2 * np.pi
    X = np.column_stack([np.ones_like(th), np.sin(th), np.cos(th)])
    beta = np.linalg.lstsq(X*sw[:,None], S*sw, rcond=None)[0]
    yhat = X @ beta
    chi2 = np.sum(((S-yhat)/err)**2)
    
    X0 = np.ones((len(th), 1))
    beta0 = np.linalg.lstsq(X0*sw[:,None], S*sw, rcond=None)[0]
    chi2_0 = np.sum(((S - X0@beta0)/err)**2)
    
    return chi2_0 - chi2, S, centers

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alice-zip", required=True); parser.add_argument("--bob-zip", required=True)
    parser.add_argument("--path-contains", default="timetags/longdist")
    parser.add_argument("--kappa", type=float, default=DEFAULT_KAPPA)
    args = parser.parse_args()

    print(f"=== WEIHS DRIFT DIAGNOSTIC (Kappa={args.kappa}) ===")
    
    # 1. MINE DATA (Quickly)
    with zipfile.ZipFile(args.alice_zip,"r") as za, zipfile.ZipFile(args.bob_zip,"r") as zb:
        # Calibration Hardcoded (from v26)
        sb, db, invS, invO, flip = 0, 1, True, False, -1
        
        # Proper list logic
        all_files = sorted([_norm_member(n) for n in za.namelist()])
        stems = [n.replace("_V.dat", "").replace("_v.dat","") for n in all_files if "longdist" in n and "_v.dat" in n.lower()]
        stems = sorted(list(set(stems)))[:10] # LIMIT TO 10 FILES FOR SPEED
        
        print(f"Scanning first {len(stems)} files...")
        
        phases_all, time_norm_all, ks_all, prods_all = [], [], [], []
        
        total_wraps = 0
        total_corr = 0
        n_files = 0
        
        engine = ThetaEngine(kappa=args.kappa)
        
        for stem in stems:
            tA, cA = load_raw_data(za, stem); tB, cB = load_raw_data(zb, stem)
            if tA is None: continue
            
            sA, oA = decode_bits(cA, sb, db, invS, invO, 1)
            sB, oB = decode_bits(cB, sb, db, False, False, 1)
            
            match = quick_sync_match(tA, tB)
            if not match: continue
            ai, bi, shift = match
            if len(ai) < 500: continue
            
            # Reconstruct Phase
            t_merge = np.concatenate([tA, tB])
            who = np.concatenate([np.zeros(len(tA), int), np.ones(len(tB), int)])
            ord_ = np.argsort(t_merge)
            
            engine.reset()
            map_b = -np.ones(len(tA), int); map_b[ai] = bi
            
            ph_stem, k_stem, p_stem, t_stem = [], [], [], []
            
            # Run stream
            t0 = tA[0]
            
            for w, i in zip(who[ord_], np.concatenate([np.arange(len(tA)), np.arange(len(tB))])[ord_]):
                if w==0:
                    engine.update(0, sA[i])
                    if map_b[i] >= 0:
                        ph_stem.append(engine.theta)
                        k_stem.append(sA[i]*2 + sB[map_b[i]])
                        p_stem.append(oA[i] * oB[map_b[i]] * flip)
                        t_stem.append(tA[i] - t0)
                else:
                    engine.update(1, sB[i])
            
            if not ph_stem: continue
            
            # DIAGNOSTIC 1: Phase Coverage
            ph_arr = np.array(ph_stem)
            wraps = (ph_arr.max() - ph_arr.min()) / (2*np.pi)
            total_wraps += wraps
            
            # DIAGNOSTIC 2: Correlation
            t_arr = np.array(t_stem)
            corr, _ = spearmanr(ph_arr, t_arr)
            total_corr += abs(corr)
            
            phases_all.append(ph_arr)
            time_norm_all.append(t_arr / t_arr.max()) # Normalize time 0..1 per file
            ks_all.append(k_stem)
            prods_all.append(p_stem)
            n_files += 1
            
            print(f"  {stem}: Wraps={wraps:.1f}, Corr(th, t)={corr:.2f}")

    # CONSOLIDATE
    phases = np.concatenate(phases_all)
    times = np.concatenate(time_norm_all)
    ks = np.concatenate(ks_all)
    prods = np.concatenate(prods_all)
    
    print("\n=== DIAGNOSTIC RESULTS ===")
    print(f"Avg Wraps per File: {total_wraps/n_files:.2f}")
    print(f"Avg |Corr(theta, time)|: {total_corr/n_files:.2f}")
    
    if total_wraps/n_files < 1.5:
        print("WARNING: Phase wraps < 1.5. Theta is effectively a linear ramp (Time Proxy).")
    
    # TEST 3: Time Binning
    # We treat 'Normalized Time' (0..1) as a phase (0..2pi) for binning
    res_time = get_stats(times, ks, prods, BIN_COUNT)
    dchi2_time = res_time[0] if res_time else 0
    print(f"\n[CONTROL 1] Time-Binning (S vs Time): DeltaChi2 = {dchi2_time:.2f}")
    
    # TEST 4: Original Theta
    # Wrap theta to 0..1 for generic binner
    th_wrapped = (phases % (2*np.pi)) / (2*np.pi)
    res_theta = get_stats(th_wrapped, ks, prods, BIN_COUNT)
    dchi2_theta = res_theta[0] if res_theta else 0
    print(f"[SIGNAL]   Theta-Binning (S vs Theta): DeltaChi2 = {dchi2_theta:.2f}")
    
    # TEST 5: Block Shuffle
    # Shuffle theta within blocks of 5000 events
    # This preserves drift but destroys geometry
    print("\n[CONTROL 2] Block-Shuffle Test...")
    ph_shuffled = phases.copy()
    n_blocks = len(ph_shuffled) // BLOCK_SIZE
    rng = np.random.default_rng(42)
    
    for i in range(n_blocks + 1):
        start = i * BLOCK_SIZE
        end = min((i+1)*BLOCK_SIZE, len(ph_shuffled))
        if end > start:
            chunk = ph_shuffled[start:end]
            rng.shuffle(chunk)
            ph_shuffled[start:end] = chunk
            
    th_shuff_wrapped = (ph_shuffled % (2*np.pi)) / (2*np.pi)
    res_shuff = get_stats(th_shuff_wrapped, ks, prods, BIN_COUNT)
    dchi2_shuff = res_shuff[0] if res_shuff else 0
    
    print(f"  Block-Shuffled DeltaChi2: {dchi2_shuff:.2f}")
    
    # CONCLUSION
    print("\n=== VERDICT ===")
    if dchi2_shuff > 0.5 * dchi2_theta:
        print("FAIL: Signal survives Block-Shuffle. It is DRIFT, not Geometry.")
    elif dchi2_time > dchi2_theta:
        print("FAIL: Time-binning explains the data better than Theta. It is DRIFT.")
    else:
        print("PASS: Signal vanishes under Block-Shuffle and beats Time-Binning.")
        print("      This supports a genuine Geometric origin.")

    # PLOT COMPARISON
    if res_time and res_theta:
        plt.figure(figsize=(10,6))
        
        # Plot Time Control
        st = res_time[1]
        xt = np.linspace(0, 2*np.pi, len(st))
        plt.plot(xt, st, 'r--', label=f'Time Control (dX2={dchi2_time:.1f})')
        
        # Plot Theta Signal
        sth = res_theta[1]
        xth = res_theta[2] * 2*np.pi
        plt.plot(xth, sth, 'b-o', label=f'Theta Signal (dX2={dchi2_theta:.1f})')
        
        plt.xlabel("Phase / Normalized Time")
        plt.ylabel("Bell Parameter S")
        plt.title(f"Drift Diagnostic (Kappa={args.kappa})")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig("weihs_drift_diagnostic.png")
        print("Saved plot: weihs_drift_diagnostic.png")

if __name__ == "__main__":
    main()