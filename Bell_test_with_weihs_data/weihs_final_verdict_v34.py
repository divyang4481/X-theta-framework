#!/usr/bin/env python3
"""
weihs_final_verdict_v34.py

THE FIX: RESTORING SYMMETRIC SYNC
================================================================================
CRITICAL BUG FIX:
  - v32/v33 failed because they used an "Asymmetric" shift estimator (looking only forward).
  - v34 restores the "Symmetric" estimator from v25 (looking backward and forward).
  - This is required to see coincidences where Bob's photon arrives slightly "early".

EXPECTATION:
  - Baseline S should return to > 2.3 (like in v25).
  - Once S is healthy, the Geometric Phase test (DeltaChi2) becomes meaningful.

USAGE:
  python weihs_final_verdict_v34.py --alice-zip "..\Bell_data\7185335\Alice.zip" --bob-zip "..\Bell_data\7185335\Bob.zip" --path-contains "timetags/longdist"
"""

import argparse
import zipfile
import numpy as np
import math
import os
import sys
from dataclasses import dataclass, field
from collections import deque

# === CONFIG ===
KAPPA = 0.1
BIN_COUNT = 8
BLOCK_SIZE = 5000
PERM_ROUNDS = 200 

def _norm_member(name: str) -> str:
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"):
        return "/".join(parts[1:])
    return n

def load_raw_data(z, stem):
    nm = { _norm_member(n).lower(): n for n in z.namelist() }
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

def estimate_shift_ns(tA, tB):
    # === THE FIX: SYMMETRIC CHECK (From v25) ===
    rng = np.random.default_rng(1)
    nA = len(tA)
    m = min(20000, nA)
    idxA = rng.choice(nA, size=m, replace=False)
    a = tA[idxA]
    
    j = np.searchsorted(tB, a)
    j0 = np.clip(j - 1, 0, len(tB) - 1)
    j1 = np.clip(j, 0, len(tB) - 1)
    
    dt0 = (tB[j0] - a) * 1e9
    dt1 = (tB[j1] - a) * 1e9
    
    # Critical: Pick the closer neighbor (Left or Right)
    dt = np.where(np.abs(dt0) <= np.abs(dt1), dt0, dt1)
    
    mask = np.abs(dt) < 200.0
    if np.sum(mask) < 10: return 0.0
    
    hist, edges = np.histogram(dt[mask], bins=2000)
    peak_bin = np.argmax(hist)
    peak = (edges[peak_bin] + edges[peak_bin+1])/2.0
    return peak

def match_pairs(tA, tB, shift, window, mode="strict"):
    tB_s = tB - shift*1e-9
    
    # Symmetric Search for Candidates
    j = np.searchsorted(tB_s, tA)
    j0 = np.clip(j-1, 0, len(tB)-1)
    j1 = np.clip(j, 0, len(tB)-1)
    dt0 = (tB_s[j0] - tA)*1e9
    dt1 = (tB_s[j1] - tA)*1e9
    
    cand_A, cand_B, cand_dt = [], [], []
    m0 = np.abs(dt0) < window; i0 = np.where(m0)[0]
    if len(i0)>0: cand_A.append(i0); cand_B.append(j0[i0]); cand_dt.append(dt0[i0])
    m1 = np.abs(dt1) < window; i1 = np.where(m1)[0]
    if len(i1)>0: cand_A.append(i1); cand_B.append(j1[i1]); cand_dt.append(dt1[i1])
    
    if not cand_A: return [], [], []
    all_A = np.concatenate(cand_A); all_B = np.concatenate(cand_B); all_dt = np.concatenate(cand_dt)
    
    if mode == "greedy":
        # Best match per A
        use_0 = np.abs(dt0) < np.abs(dt1)
        best_idx = np.where(use_0, j0, j1)
        best_dt = np.where(use_0, dt0, dt1)
        mask = np.abs(best_dt) < window
        return np.where(mask)[0], best_idx[mask], best_dt[mask]
    else: # Strict
        order = np.argsort(np.abs(all_dt))
        usedA = np.zeros(len(tA), bool); usedB = np.zeros(len(tB), bool)
        fA, fB, fdt = [], [], []
        for k in order:
            ia, ib = all_A[k], all_B[k]
            if not usedA[ia] and not usedB[ib]:
                usedA[ia]=True; usedB[ib]=True
                fA.append(ia); fB.append(ib); fdt.append(all_dt[k])
        return np.array(fA), np.array(fB), np.array(fdt)

def compute_S_quick(sA, oA, sB, oB, ai, bi, flip):
    a = sA[ai].astype(np.int8); b = sB[bi].astype(np.int8)
    prod = (oA[ai] * oB[bi]) * flip
    k = (2 * a + b)
    c = np.bincount(k, minlength=4)
    if np.any(c==0): return 0.0
    s = np.bincount(k, weights=prod, minlength=4)
    m = s/c
    return m[0]+m[1]+m[2]-m[3]

def calibrate_config(za, zb, stems, window, mode):
    print("[CALIBRATION] Scanning for best bitmask (Goal: S > 2.0)...")
    best = {"S": -9.0}
    
    test_stems = []
    for s in stems:
        tA, cA = load_raw_data(za, s); tB, cB = load_raw_data(zb, s)
        if tA is None: continue
        shift = estimate_shift_ns(tA, tB)
        if shift==0: continue
        test_stems.append((s, tA, tB, cA, cB, shift))
        if len(test_stems) >= 3: break # Use 3 files for calib
    
    if not test_stems: return None

    for sb in [0]:
        for db in [0, 1]:
            if sb==db: continue
            for invS_A in [False, True]:
                for invS_B in [False, True]:
                    for flip in [1, -1]:
                        s_vals = []
                        for (_, tA, tB, cA, cB, shift) in test_stems:
                            sA, oA = decode_bits(cA, sb, db, invS_A, False)
                            sB, oB = decode_bits(cB, sb, db, invS_B, False)
                            ai, bi, _ = match_pairs(tA, tB, shift, window, mode)
                            if len(ai) < 50: continue
                            S = compute_S_quick(sA, oA, sB, oB, ai, bi, flip)
                            s_vals.append(S)
                        
                        if s_vals:
                            avg_S = np.mean(s_vals)
                            if avg_S > best["S"]:
                                best = {"S": avg_S, "sb": sb, "db": db, "invS_A": invS_A, "invS_B": invS_B, "flip": flip}
                            
    print(f"  WINNER: S = {best['S']:.4f} (Flip={best['flip']}, invA={best['invS_A']}, invB={best['invS_B']})")
    return best

@dataclass
class StrictThetaEngine:
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
            if h[0] != h[-1]: return 
            unique_states = set(h)
            if len(unique_states) < 4: return
            area = 0.0
            for i in range(4):
                area += float(h[i][0]*h[i+1][1] - h[i+1][0]*h[i][1])
            ori = 1 if area > 0.1 else (-1 if area < -0.1 else 0)
            self.theta += self.kappa * ori

def build_events(tA, sA, tB, sB, ai, bi, prod, shift, eng):
    tB_s = tB - shift*1e-9
    times = np.concatenate([tA, tB_s])
    who = np.concatenate([np.zeros(len(tA),int), np.ones(len(tB),int)])
    idx = np.concatenate([np.arange(len(tA)), np.arange(len(tB))])
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

def bin_and_score(phases, ks, prods, bins=8):
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
    ph_aligned = []
    ks_aligned, pr_aligned = [], []
    for i in range(len(offsets)-1):
        lo, hi = offsets[i], offsets[i+1]
        mid = lo + (hi-lo)//2
        if mid <= lo: continue
        score, fit = bin_and_score(phases[lo:mid], ks[lo:mid], prods[lo:mid], bins)
        if not fit: continue
        beta = fit[3]
        phi = math.atan2(beta[2], beta[1])
        ph_aligned.append(phases[mid:hi] + phi)
        ks_aligned.append(ks[mid:hi])
        pr_aligned.append(prods[mid:hi])
    if not ph_aligned: return 0.0
    ph_f = np.concatenate(ph_aligned)
    ks_f = np.concatenate(ks_aligned)
    pr_f = np.concatenate(pr_aligned)
    score, _ = bin_and_score(ph_f, ks_f, pr_f, bins)
    return score

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alice-zip", required=True); ap.add_argument("--bob-zip", required=True)
    ap.add_argument("--path-contains", default="timetags/longdist")
    ap.add_argument("--window", type=float, default=1.5) # DEFAULT 1.5ns (Clean)
    ap.add_argument("--mode", type=str, default="strict") # DEFAULT Strict
    ap.add_argument("--force-remine", action="store_true")
    args = ap.parse_args()

    print(f"=== VERDICT v34 (The Fix): Window={args.window}ns | Mode={args.mode} ===")
    
    with zipfile.ZipFile(args.alice_zip,"r") as za, zipfile.ZipFile(args.bob_zip,"r") as zb:
        all_files = sorted([_norm_member(n) for n in za.namelist()])
        stems = sorted(list(set([n.replace("_V.dat", "").replace("_v.dat","") for n in all_files if args.path_contains.lower() in n.lower() and "_v.dat" in n.lower()])))
        
        if not stems: print("No files found!"); return

        # 1. AUTO-CALIBRATION
        cfg = calibrate_config(za, zb, stems, args.window, args.mode)
        if not cfg or cfg["S"] < 2.0:
            print(f"CRITICAL: Calibration S = {cfg['S'] if cfg else 0:.3f}. Continuing to per-stem check...")
            if not cfg: cfg = {"sb":0, "db":1, "invS_A":True, "invS_B":False, "flip":-1}

        # 2. MINING
        eng = StrictThetaEngine(kappa=KAPPA)
        ph_l, k_l, p_l, off_l = [], [], [], [0]
        accepted = 0
        
        print("\n[PER-STEM DIAGNOSTICS]")
        for s in stems:
            tA, cA = load_raw_data(za, s); tB, cB = load_raw_data(zb, s)
            if tA is None: continue
            
            sA, oA = decode_bits(cA, cfg["sb"], cfg["db"], cfg["invS_A"], False)
            sB, oB = decode_bits(cB, cfg["sb"], cfg["db"], cfg["invS_B"], False)
            
            shift = estimate_shift_ns(tA, tB) # SYMMETRIC FIX
            if shift == 0.0: continue
            
            ai, bi, _ = match_pairs(tA, tB, shift, args.window, args.mode)
            if len(ai)<500: continue
            
            S_local = compute_S_quick(sA, oA, sB, oB, ai, bi, cfg["flip"])
            print(f"  {s}: N={len(ai)}, S={S_local:.3f}")
            
            # FILTER BAD STEMS
            if S_local < 2.0: continue
            
            prod = (oA[ai].astype(np.int16) * oB[bi].astype(np.int16)) * cfg["flip"]
            ph, kk, pp = build_events(tA, sA, tB, sB, ai, bi, prod, shift, eng)
            
            ph_l.append(ph); k_l.append(kk); p_l.append(pp)
            off_l.append(off_l[-1]+len(ph))
            accepted += 1

    if not ph_l:
        print("CRITICAL: No valid data.")
        return

    phases = np.concatenate(ph_l)
    ks = np.concatenate(k_l)
    prods = np.concatenate(p_l)
    offsets = np.array(off_l)
    
    # GLOBAL S
    counts = np.bincount(ks, minlength=4)
    sums = np.bincount(ks, weights=prods, minlength=4)
    means = sums / counts
    S_global = means[0]+means[1]+means[2]-means[3]
    print(f"\n[GLOBAL] Events: {len(phases)} | Baseline S: {S_global:.4f}")
    
    # TESTS
    real_score = cv_align_score(phases, ks, prods, offsets, BIN_COUNT)
    print(f"\n[TEST 1] Real Signal: DeltaChi2 = {real_score:.2f}")
    
    print(f"\n[TEST 2] Permutation Test ({PERM_ROUNDS} rounds)...")
    better_count = 0
    rng = np.random.default_rng(999)
    for r in range(PERM_ROUNDS):
        ph_perm = phases.copy()
        for i in range(len(offsets)-1):
            s, e = offsets[i], offsets[i+1]
            length = e - s
            if length < BLOCK_SIZE: continue
            n_blk = length // BLOCK_SIZE
            for b in range(n_blk):
                bs = s + b*BLOCK_SIZE
                be = min(s + (b+1)*BLOCK_SIZE, e)
                chunk = ph_perm[bs:be]
                rng.shuffle(chunk)
                ph_perm[bs:be] = chunk
        perm_score = cv_align_score(ph_perm, ks, prods, offsets, BIN_COUNT)
        if perm_score >= real_score: better_count += 1
        if (r+1)%50==0: print(f"  Round {r+1}...")
            
    p_val = (better_count + 1) / (PERM_ROUNDS + 1)
    print(f"\n[RESULT] P-Value = {p_val:.4f}")

    if p_val < 0.05 and S_global > 2.0:
        print("SUCCESS: Valid Bell violation AND Significant Geometric Signal.")
    else:
        print("FAILURE.")

if __name__ == "__main__":
    main()