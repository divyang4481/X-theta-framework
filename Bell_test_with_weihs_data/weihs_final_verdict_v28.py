#!/usr/bin/env python3
"""
weihs_final_verdict_v28.py

THE FINAL VERDICT: DRIFT vs GEOMETRY
================================================================================
PURPOSE:
  To determine if the signal found in v26 (DeltaChi2 ~ 29) is real Geometric Phase
  or just the massive Time Drift (DeltaChi2 ~ 216) leaking in.

METHOD:
  We use the 'v26' pipeline (Strict Matching + CV Alignment) which we KNOW 
  can see the signal.
  
  We compare the Real Signal against a "Block Shuffled" signal.
  - Block Shuffle preserves Time Drift (locally).
  - But it destroys Geometric Topology.

VERDICT LOGIC:
  - If Block Score is HIGH (close to Real Score) -> It's DRIFT (Artifact).
  - If Block Score is LOW (close to 0) -> It's GEOMETRY (Discovery).

USAGE:
  python weihs_final_verdict_v28.py --alice-zip "..\Bell_data\7185335\Alice.zip" --bob-zip "..\Bell_data\7185335\Bob.zip" --path-contains "timetags/longdist"
"""

import argparse
import zipfile
import numpy as np
import math
import os
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from collections import deque

# === CONFIG ===
KAPPA = 0.1
BIN_COUNT = 8
BLOCK_SIZE = 5000  # Number of events to shuffle within (preserves slow drift)

# ============================================================
# CORE PIPELINE (v26 ENGINE)
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

def decode_from_c(c, sb, db, invS, invO):
    s = ((c >> sb) & 1).astype(np.int32)
    det = ((c >> db) & 1).astype(np.int32)
    if invS: s = 1 - s
    o = (2 * det - 1)
    if invO: o = -o
    return s.astype(np.uint8), o.astype(np.int8)

def estimate_shift_ns(tA, tB):
    m = min(20000, len(tA))
    dt0 = tB[:m] - tA[:m] # Rough check
    # Assume pre-sync or quick scan (simplified for brevity as v26 works)
    # Actually, let's copy the robust v26 sync to be safe
    rng = np.random.default_rng(1)
    idxA = rng.choice(len(tA), size=m, replace=False)
    a = tA[idxA]
    j = np.searchsorted(tB, a)
    j = np.clip(j, 0, len(tB)-1)
    dt = (tB[j] - a) * 1e9
    mask = np.abs(dt) < 200.0
    if np.sum(mask) < 100: return 0.0
    hist, edges = np.histogram(dt[mask], bins=400)
    peak = (edges[np.argmax(hist)] + edges[np.argmax(hist)+1])/2.0
    return peak

def match_strict(tA, tB, shift):
    tB_s = tB - shift*1e-9
    j = np.searchsorted(tB_s, tA)
    j0 = np.clip(j-1, 0, len(tB)-1)
    j1 = np.clip(j, 0, len(tB)-1)
    dt0 = (tB_s[j0] - tA)*1e9
    dt1 = (tB_s[j1] - tA)*1e9
    
    cand_A, cand_B, cand_dt = [], [], []
    
    m0 = np.abs(dt0) < 1.5; i0 = np.where(m0)[0]
    if len(i0)>0: cand_A.append(i0); cand_B.append(j0[i0]); cand_dt.append(dt0[i0])
        
    m1 = np.abs(dt1) < 1.5; i1 = np.where(m1)[0]
    if len(i1)>0: cand_A.append(i1); cand_B.append(j1[i1]); cand_dt.append(dt1[i1])
        
    if not cand_A: return [], [], []
    all_A = np.concatenate(cand_A); all_B = np.concatenate(cand_B); all_dt = np.concatenate(cand_dt)
    
    # Sort by quality
    order = np.argsort(np.abs(all_dt))
    usedA = np.zeros(len(tA), bool); usedB = np.zeros(len(tB), bool)
    fA, fB, fdt = [], [], []
    
    for k in order:
        ia, ib = all_A[k], all_B[k]
        if not usedA[ia] and not usedB[ib]:
            usedA[ia]=True; usedB[ib]=True
            fA.append(ia); fB.append(ib); fdt.append(all_dt[k])
            
    return np.array(fA), np.array(fB), np.array(fdt)

@dataclass
class ThetaEngine:
    kappa: float
    a: int=0; b: int=0; theta: float=0.0
    history: deque = field(default_factory=lambda: deque(maxlen=5))
    def reset(self): self.a=0; self.b=0; self.theta=0.0; self.history.clear()
    def update(self, w, v):
        if w==0: self.a=v
        else: self.b=v
        self.history.append((self.a, self.b))
        if len(self.history)==5:
            h=list(self.history)
            area = 0.0
            for i in range(4): area += float(h[i][0]*h[i+1][1] - h[i+1][0]*h[i][1])
            ori = 1 if area>0.1 else (-1 if area<-0.1 else 0)
            self.theta += self.kappa*ori

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

# ============================================================
# ANALYSIS
# ============================================================

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
        
        # Train
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
    ap.add_argument("--cache", default="weihs_verdict.npz")
    ap.add_argument("--force-remine", action="store_true")
    args = ap.parse_args()

    print("=== FINAL VERDICT DIAGNOSTIC (v28) ===")
    
    # 1. LOAD/MINE DATA
    if args.force_remine or not os.path.exists(args.cache):
        print("[MINING] Extracting data using v26 standards...")
        with zipfile.ZipFile(args.alice_zip,"r") as za, zipfile.ZipFile(args.bob_zip,"r") as zb:
            # Stems
            all_files = sorted([_norm_member(n) for n in za.namelist()])
            stems = sorted(list(set([n.replace("_V.dat", "").replace("_v.dat","") for n in all_files if args.path_contains in n and "_v.dat" in n.lower()])))
            
            # Calibration Winner Config
            sb, db, invS, invO, flip = 0, 1, True, False, -1
            
            eng = ThetaEngine(kappa=KAPPA)
            ph_l, k_l, p_l, off_l = [], [], [], [0]
            
            for s in stems:
                tA, cA = load_raw_data(za, s); tB, cB = load_raw_data(zb, s)
                if tA is None: continue
                
                sA, oA = decode_from_c(cA, sb, db, invS, invO)
                sB, oB = decode_from_c(cB, sb, db, False, False)
                shift = estimate_shift_ns(tA, tB)
                ai, bi, _ = match_strict(tA, tB, shift)
                if len(ai)<500: continue
                
                prod = (oA[ai].astype(np.int16) * oB[bi].astype(np.int16)) * flip
                ph, kk, pp = build_events(tA, sA, tB, sB, ai, bi, prod, shift, eng)
                
                ph_l.append(ph); k_l.append(kk); p_l.append(pp)
                off_l.append(off_l[-1]+len(ph))
                
            np.savez(args.cache, phases=np.concatenate(ph_l), ks=np.concatenate(k_l), 
                     prods=np.concatenate(p_l), offsets=np.array(off_l))
            print(f"[CACHE] Saved {len(ph_l)} stems.")

    # 2. RUN TEST
    d = np.load(args.cache)
    phases = d["phases"]; ks = d["ks"]; prods = d["prods"]; offsets = d["offsets"]
    
    print("\n[TEST 1] Real Signal (Aligned)")
    real_score = cv_align_score(phases, ks, prods, offsets, BIN_COUNT)
    print(f"  DeltaChi2 = {real_score:.2f}")
    
    print("\n[TEST 2] Block-Shuffle Control (Preserves Drift, Kills Geometry)")
    # Shuffle theta within 5000-event blocks
    # If the signal is Drift, shuffling theta locally won't kill it (because drift is slow/local)
    # If the signal is Geometry, shuffling theta locally will kill the fine-structure
    
    ph_shuff = phases.copy()
    n_blocks = len(ph_shuff)//BLOCK_SIZE
    rng = np.random.default_rng(999)
    
    for i in range(n_blocks+1):
        s = i*BLOCK_SIZE; e = min((i+1)*BLOCK_SIZE, len(ph_shuff))
        if e > s:
            chunk = ph_shuff[s:e]
            rng.shuffle(chunk)
            ph_shuff[s:e] = chunk
            
    # CRITICAL: We must re-run alignment on the shuffled data!
    # Because if drift is present, the alignment algo might latch onto it.
    block_score = cv_align_score(ph_shuff, ks, prods, offsets, BIN_COUNT)
    print(f"  DeltaChi2 = {block_score:.2f}")
    
    # VERDICT
    print("\n=== FINAL VERDICT ===")
    ratio = block_score / (real_score + 0.01)
    
    if real_score < 5.0:
        print("INCONCLUSIVE: Real signal is too weak.")
    elif ratio > 0.5:
        print("FAIL: Signal is DRIFT. (Block shuffle didn't kill it)")
        print(f"      The signal is largely time-dependent artifacts ({ratio*100:.1f}% survival).")
    else:
        print("PASS: Signal is GEOMETRY. (Block shuffle killed it)")
        print(f"      The signal depends on the exact topological order ({100-ratio*100:.1f}% drop).")

if __name__ == "__main__":
    main()