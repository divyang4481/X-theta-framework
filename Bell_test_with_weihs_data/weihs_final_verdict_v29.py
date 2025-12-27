#!/usr/bin/env python3
"""
weihs_final_verdict_v29.py

THE FINAL VERDICT (Integrated v26 Engine + Verdict Logic)
================================================================================
PURPOSE:
  1. Mines the data using the robust v26 engine (Strict 1-to-1 matching).
  2. Determines if the signal at Kappa=0.1 is GEOMETRY or DRIFT.

FIXES:
  - Restored robust, case-insensitive file matching (fixes "No files matched").
  - Hardcoded the known-good bitmask configuration (S=2.35) to ensure stability.
  - Generates 'weihs_verdict_v29.npz' automatically.

VERDICT LOGIC:
  - We calculate DeltaChi2 for the Real Signal.
  - We Block-Shuffle the phases (preserving time-drift) and re-calculate.
  - If Signal vanishes -> Discovery (Geometry).
  - If Signal survives -> Artifact (Drift).

USAGE:
  python weihs_final_verdict_v29.py --alice-zip "..\Bell_data\7185335\Alice.zip" --bob-zip "..\Bell_data\7185335\Bob.zip" --path-contains "timetags/longdist"
"""

import argparse
import zipfile
import numpy as np
import math
import os
import sys
from dataclasses import dataclass, field
from collections import deque

# === CONFIGURATION ===
KAPPA = 0.1
BIN_COUNT = 8
BLOCK_SIZE = 5000 

# Known "Gold Standard" Config from previous successful runs
# (set=0, det=1, invS_A=True, flip=-1)
CFG_SB = 0
CFG_DB = 1
CFG_INVS_A = True
CFG_INVO_A = False
CFG_INVS_B = False
CFG_INVO_B = False
CFG_FLIP = -1

# ============================================================
# 1. ROBUST UTILITIES (From v20/v26)
# ============================================================

def _norm_member(name: str) -> str:
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"):
        return "/".join(parts[1:])
    return n

def list_stems(z: zipfile.ZipFile, path_contains: str = None):
    """Robust case-insensitive stem finder."""
    names = [_norm_member(n) for n in z.namelist()]
    
    # Filter by path string (case insensitive)
    if path_contains:
        pc = path_contains.lower().replace("\\", "/")
        names = [n for n in names if pc in n.lower()]
        
    stems = set()
    for n in names:
        if n.lower().endswith("_v.dat"):
            stems.add(n[:-6]) # Remove _V.dat
            
    # Verify pairs exist
    nm_lower = {n.lower(): n for n in names}
    valid_stems = []
    for s in stems:
        v_key = (s + "_v.dat").lower()
        c_key = (s + "_c.dat").lower()
        if v_key in nm_lower and c_key in nm_lower:
            valid_stems.append(s)
            
    return sorted(valid_stems)

def load_raw_data(z, stem):
    nm = { _norm_member(n).lower(): n for n in z.namelist() }
    # Fuzzy match
    vk = next((k for k in nm if k.endswith("_v.dat") and stem.lower() in k), None)
    ck = next((k for k in nm if k.endswith("_c.dat") and stem.lower() in k), None)
    
    if not vk or not ck: return None, None
    
    with z.open(nm[vk]) as f: t = np.frombuffer(f.read(), ">f8")
    with z.open(nm[ck]) as f: c = np.frombuffer(f.read(), ">u2")
    n = min(len(t), len(c))
    return t[:n].copy(), c[:n].copy()

def decode_bits(c, sb, db, invS, invO):
    s = ((c >> sb) & 1).astype(np.int32)
    det = ((c >> db) & 1).astype(np.int32)
    if invS: s = 1 - s
    o = (2 * det - 1)
    if invO: o = -o
    return s.astype(np.uint8), o.astype(np.int8)

# ============================================================
# 2. V26 MINING ENGINE
# ============================================================

def estimate_shift_ns(tA, tB):
    # Robust coarse-to-fine sync
    rng = np.random.default_rng(1)
    nA = len(tA)
    m = min(20000, nA)
    idxA = rng.choice(nA, size=m, replace=False)
    a = tA[idxA]
    
    j = np.searchsorted(tB, a)
    j = np.clip(j, 0, len(tB)-1)
    dt = (tB[j] - a) * 1e9
    
    # Filter crazy values
    mask = np.abs(dt) < 200000.0 
    if np.sum(mask) < 100: return 0.0
    
    hist, edges = np.histogram(dt[mask], bins=5000)
    peak_coarse = (edges[np.argmax(hist)] + edges[np.argmax(hist)+1])/2.0
    
    # Fine tune
    mask_fine = np.abs(dt - peak_coarse) < 100
    if np.sum(mask_fine) < 10: return peak_coarse
    
    hist_f, edges_f = np.histogram(dt[mask_fine], bins=500)
    peak_fine = (edges_f[np.argmax(hist_f)] + edges_f[np.argmax(hist_f)+1])/2.0
    return peak_fine

def match_strict(tA, tB, shift):
    # Strict 1-to-1 Matching (Reviewer Requirement)
    tB_s = tB - shift*1e-9
    
    j = np.searchsorted(tB_s, tA)
    j0 = np.clip(j-1, 0, len(tB)-1)
    j1 = np.clip(j, 0, len(tB)-1)
    
    dt0 = (tB_s[j0] - tA)*1e9
    dt1 = (tB_s[j1] - tA)*1e9
    
    cand_A, cand_B, cand_dt = [], [], []
    
    # Collect candidates in window 1.5ns
    m0 = np.abs(dt0) < 1.5; i0 = np.where(m0)[0]
    if len(i0)>0: cand_A.append(i0); cand_B.append(j0[i0]); cand_dt.append(dt0[i0])
        
    m1 = np.abs(dt1) < 1.5; i1 = np.where(m1)[0]
    if len(i1)>0: cand_A.append(i1); cand_B.append(j1[i1]); cand_dt.append(dt1[i1])
    
    if not cand_A: return [], [], []
    
    all_A = np.concatenate(cand_A)
    all_B = np.concatenate(cand_B)
    all_dt = np.concatenate(cand_dt)
    
    # Sort by time difference (best match first)
    order = np.argsort(np.abs(all_dt))
    
    usedA = np.zeros(len(tA), bool)
    usedB = np.zeros(len(tB), bool)
    fA, fB, fdt = [], [], []
    
    for k in order:
        ia, ib = all_A[k], all_B[k]
        if not usedA[ia] and not usedB[ib]:
            usedA[ia] = True
            usedB[ib] = True
            fA.append(ia); fB.append(ib); fdt.append(all_dt[k])
            
    return np.array(fA), np.array(fB), np.array(fdt)

@dataclass
class ThetaEngine:
    kappa: float
    a: int=0; b: int=0; theta: float=0.0
    history: deque = field(default_factory=lambda: deque(maxlen=5))
    
    def reset(self):
        self.a=0; self.b=0; self.theta=0.0; self.history.clear()
        
    def update(self, w, v):
        if w==0: self.a=v
        else: self.b=v
        self.history.append((self.a, self.b))
        if len(self.history)==5:
            h = list(self.history)
            area = 0.0
            for i in range(4):
                area += float(h[i][0]*h[i+1][1] - h[i+1][0]*h[i][1])
            ori = 1 if area > 0.1 else (-1 if area < -0.1 else 0)
            self.theta += self.kappa * ori

def build_events(tA, sA, tB, sB, ai, bi, prod, shift, eng):
    # Reconstruct timeline for phase calculation
    tB_s = tB - shift*1e-9
    times = np.concatenate([tA, tB_s])
    who = np.concatenate([np.zeros(len(tA),int), np.ones(len(tB),int)])
    idx = np.concatenate([np.arange(len(tA)), np.arange(len(tB))])
    
    # Sort by time
    order = np.argsort(times)
    
    mapB = -np.ones(len(tA), int); mapB[ai] = bi
    mapP = np.zeros(len(tA), int); mapP[ai] = prod
    
    ph, ks, pr = [], [], []
    eng.reset()
    
    for w, i in zip(who[order], idx[order]):
        if w==0:
            eng.update(0, int(sA[i]))
            if mapB[i] >= 0:
                ph.append(eng.theta)
                ks.append(int(sA[i]*2 + sB[mapB[i]]))
                pr.append(mapP[i])
        else:
            eng.update(1, int(sB[i]))
            
    return np.array(ph), np.array(ks), np.array(pr)

# ============================================================
# 3. ANALYSIS LOGIC
# ============================================================

def bin_and_score(phases, ks, prods, bins=8):
    # Standard Binning
    b = ((phases % (2*np.pi))/(2*np.pi)*bins).astype(int) % bins
    flat = b*4 + ks
    c = np.bincount(flat, minlength=bins*4).reshape(bins, 4)
    s = np.bincount(flat, weights=prods, minlength=bins*4).reshape(bins, 4)
    
    valid = np.all(c > 5, axis=1)
    if np.sum(valid) < 4: return 0.0, None
    
    m = s[valid]/c[valid]
    S = m[:,0]+m[:,1]+m[:,2]-m[:,3]
    err = np.sqrt(np.sum((1-m**2)/c[valid], axis=1))
    th = (np.where(valid)[0]+0.5)*(2*np.pi/bins)
    
    # Weighted Fit
    w = 1.0/np.maximum(err, 1e-9)**2; sw = np.sqrt(w)
    X = np.column_stack([np.ones_like(th), np.sin(th), np.cos(th)])
    beta = np.linalg.lstsq(X*sw[:,None], S*sw, rcond=None)[0]
    yhat = X @ beta
    chi2 = np.sum(((S-yhat)/err)**2)
    
    X0 = np.ones((len(th), 1))
    beta0 = np.linalg.lstsq(X0*sw[:,None], S*sw, rcond=None)[0]
    chi2_0 = np.sum(((S - X0@beta0)/err)**2)
    
    return chi2_0 - chi2, (th, S, err, beta)

def cv_align_score(phases, ks, prods, offsets, bins):
    # Cross-Validated Alignment
    ph_aligned = []
    ks_aligned, pr_aligned = [], []
    
    for i in range(len(offsets)-1):
        lo, hi = offsets[i], offsets[i+1]
        mid = lo + (hi-lo)//2
        if mid <= lo: continue
        
        # Train on first half
        score, fit = bin_and_score(phases[lo:mid], ks[lo:mid], prods[lo:mid], bins)
        if not fit: continue
        
        # Align second half
        beta = fit[3]
        phi = math.atan2(beta[2], beta[1])
        
        ph_aligned.append(phases[mid:hi] + phi)
        ks_aligned.append(ks[mid:hi])
        pr_aligned.append(prods[mid:hi])
        
    if not ph_aligned: return 0.0
    
    # Score Combined
    ph_f = np.concatenate(ph_aligned)
    ks_f = np.concatenate(ks_aligned)
    pr_f = np.concatenate(pr_aligned)
    
    score, _ = bin_and_score(ph_f, ks_f, pr_f, bins)
    return score

# ============================================================
# MAIN
# ============================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alice-zip", required=True); ap.add_argument("--bob-zip", required=True)
    ap.add_argument("--path-contains", default="")
    ap.add_argument("--cache", default="weihs_verdict_v29.npz")
    ap.add_argument("--force-remine", action="store_true")
    args = ap.parse_args()

    print("=== FINAL VERDICT DIAGNOSTIC (v29) ===")
    
    # 1. MINING (If needed)
    if args.force_remine or not os.path.exists(args.cache):
        print(f"[MINING] Scanning for '{args.path_contains or 'all'}'...")
        with zipfile.ZipFile(args.alice_zip,"r") as za, zipfile.ZipFile(args.bob_zip,"r") as zb:
            stems = list_stems(za, args.path_contains)
            print(f"  Found {len(stems)} matching file pairs.")
            if not stems:
                print("CRITICAL: No files found. Check --path-contains.")
                return

            # Apply Hardcoded Gold Config
            print("[CONFIG] Using calibrated config (S=2.35 settings)")
            
            eng = ThetaEngine(kappa=KAPPA)
            ph_l, k_l, p_l, off_l = [], [], [], [0]
            accepted = 0
            
            for s in stems:
                tA, cA = load_raw_data(za, s); tB, cB = load_raw_data(zb, s)
                if tA is None: continue
                
                sA, oA = decode_bits(cA, CFG_SB, CFG_DB, CFG_INVS_A, CFG_INVO_A)
                sB, oB = decode_bits(cB, CFG_SB, CFG_DB, CFG_INVS_B, CFG_INVO_B)
                
                shift = estimate_shift_ns(tA, tB)
                if shift == 0.0: continue
                
                # Strict Matching
                ai, bi, _ = match_strict(tA, tB, shift)
                if len(ai) < 500: continue
                
                prod = (oA[ai].astype(np.int16) * oB[bi].astype(np.int16)) * CFG_FLIP
                ph, kk, pp = build_events(tA, sA, tB, sB, ai, bi, prod, shift, eng)
                
                ph_l.append(ph); k_l.append(kk); p_l.append(pp)
                off_l.append(off_l[-1]+len(ph))
                accepted += 1
                if accepted % 5 == 0: print(f"  Processed {accepted} files...")
            
            if accepted == 0:
                print("CRITICAL: No files passed sync/matching.")
                return
                
            np.savez(args.cache, phases=np.concatenate(ph_l), ks=np.concatenate(k_l), 
                     prods=np.concatenate(p_l), offsets=np.array(off_l))
            print(f"[CACHE] Saved {accepted} stems to {args.cache}")

    # 2. VERDICT TESTS
    d = np.load(args.cache)
    phases = d["phases"]; ks = d["ks"]; prods = d["prods"]; offsets = d["offsets"]
    print(f"\n[DATA] Loaded {len(phases)} matched events.")
    
    # Test A: Real Signal
    real_score = cv_align_score(phases, ks, prods, offsets, BIN_COUNT)
    print(f"\n[TEST 1] Real Signal (Aligned): DeltaChi2 = {real_score:.2f}")
    
    # Test B: Block Shuffle
    print("\n[TEST 2] Block-Shuffle Control (Drift Preservation)")
    ph_shuff = phases.copy()
    n_blocks = len(ph_shuff)//BLOCK_SIZE
    rng = np.random.default_rng(999)
    
    for i in range(n_blocks+1):
        s = i*BLOCK_SIZE; e = min((i+1)*BLOCK_SIZE, len(ph_shuff))
        if e > s:
            chunk = ph_shuff[s:e]
            rng.shuffle(chunk)
            ph_shuff[s:e] = chunk
            
    # RE-ALIGN NOISE
    block_score = cv_align_score(ph_shuff, ks, prods, offsets, BIN_COUNT)
    print(f"  DeltaChi2 = {block_score:.2f}")
    
    # RESULT
    print("\n=== FINAL VERDICT ===")
    if real_score < 5.0:
        print("INCONCLUSIVE: Real signal too weak to test.")
    elif block_score < 0.5 * real_score:
        print("PASS: Signal is GEOMETRY (Discovery).")
        print(f"      Shuffle destroyed the signal ({block_score:.2f} << {real_score:.2f}).")
    else:
        print("FAIL: Signal is DRIFT (Artifact).")
        print(f"      Signal survived time-preserving shuffle ({block_score:.2f} approx {real_score:.2f}).")

if __name__ == "__main__":
    main()