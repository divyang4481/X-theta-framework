#!/usr/bin/env python3
"""
weihs_final_verdict_v36.py

THE RIGOROUS VERDICT: TEMPORAL ALIGNMENT & DRIFT-PRESERVING NULL
================================================================================
FIXES APPLIED:
  1. TEMPORAL ALIGNMENT: Products computed inside build_events() loop.
  2. EXACT STEM RESOLUTION: Uses dictionary nm.get() for exact file matching.
  3. RESIDUAL SHUFFLE: Shuffles remainder chunks in the block-permutation test.
  4. SYNC RELIABILITY: Only skips stems if histogram matches are too low.
  5. S-GUARD: Prevents divide-by-zero in global CHSH calculation.

USAGE:
  python weihs_final_verdict_v36.py --alice-zip "..\Bell_data\7185335\Alice.zip" --bob-zip "..\Bell_data\7185335\Bob.zip" --window 1.0 --mode strict
"""

import argparse
import zipfile
import numpy as np
import math
import os
from dataclasses import dataclass, field
from collections import deque

# === CONFIG ===
KAPPA = 0.1
BIN_COUNT = 8
BLOCK_SIZE = 5000
PERM_ROUNDS = 2000 

# Calibrated Config from previous successful runs
CFG = {"sb": 0, "db": 1, "invS_A": True, "invS_B": False, "flip": -1}

def _norm_member(name: str) -> str:
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"):
        return "/".join(parts[1:])
    return n

def load_raw_data(z, stem):
    """Exact lookup to prevent stem collisions (e.g., longdist1 vs 10)."""
    nm = { _norm_member(n).lower(): n for n in z.namelist() }
    stem0 = stem.lower()
    vk = nm.get(f"{stem0}_v.dat")
    ck = nm.get(f"{stem0}_c.dat")
    
    if not vk or not ck: return None, None
    
    with z.open(vk) as f: t = np.frombuffer(f.read(), ">f8")
    with z.open(ck) as f: c = np.frombuffer(f.read(), ">u2")
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
    """Symmetric shift estimator with reliability check."""
    rng = np.random.default_rng(1)
    nA = len(tA)
    m = min(20000, nA)
    idxA = rng.choice(nA, size=m, replace=False)
    a = tA[idxA]
    j = np.searchsorted(tB, a)
    j0, j1 = np.clip(j-1, 0, len(tB)-1), np.clip(j, 0, len(tB)-1)
    dt0, dt1 = (tB[j0]-a)*1e9, (tB[j1]-a)*1e9
    dt = np.where(np.abs(dt0) <= np.abs(dt1), dt0, dt1)
    mask = np.abs(dt) < 200.0
    if np.sum(mask) < 20: return None # Reliability check
    hist, edges = np.histogram(dt[mask], bins=2000)
    return (edges[np.argmax(hist)] + edges[np.argmax(hist)+1])/2.0

def match_pairs(tA, tB, shift, window, mode="strict"):
    tB_s = tB - shift*1e-9
    j = np.searchsorted(tB_s, tA)
    j0, j1 = np.clip(j-1, 0, len(tB)-1), np.clip(j, 0, len(tB)-1)
    dt0, dt1 = (tB_s[j0]-tA)*1e9, (tB_s[j1]-tA)*1e9
    
    if mode == "greedy":
        best_dt = np.where(np.abs(dt0) < np.abs(dt1), dt0, dt1)
        best_idx = np.where(np.abs(dt0) < np.abs(dt1), j0, j1)
        mask = np.abs(best_dt) < window
        return np.where(mask)[0], best_idx[mask]
    else: # Strict 1-to-1 matching
        cand_A, cand_B, cand_dt = [], [], []
        m0, m1 = np.abs(dt0) < window, np.abs(dt1) < window
        if np.any(m0): cand_A.append(np.where(m0)[0]); cand_B.append(j0[m0]); cand_dt.append(dt0[m0])
        if np.any(m1): cand_A.append(np.where(m1)[0]); cand_B.append(j1[m1]); cand_dt.append(dt1[m1])
        if not cand_A: return [], []
        a_arr, b_arr, d_arr = np.concatenate(cand_A), np.concatenate(cand_B), np.concatenate(cand_dt)
        order = np.argsort(np.abs(d_arr))
        uA, uB = np.zeros(len(tA), bool), np.zeros(len(tB), bool)
        fA, fB = [], []
        for k in order:
            ia, ib = a_arr[k], b_arr[k]
            if not uA[ia] and not uB[ib]:
                uA[ia], uB[ib] = True, True
                fA.append(ia); fB.append(ib)
        return np.array(fA), np.array(fB)

@dataclass
class StrictThetaEngine:
    kappa: float
    a: int=0; b: int=0; theta: float=0.0
    history: deque = field(default_factory=lambda: deque(maxlen=5))
    def reset(self): self.a=0; self.b=0; self.theta=0.0; self.history.clear()
    def update(self, w, v):
        if w==0: self.a=v
        else: self.b=v
        self.history.append((self.a, self.b))
        if len(self.history)==5:
            h = list(self.history)
            if h[0] == h[-1] and len(set(h)) == 4: # Strict Plaquette check
                area = 0.0
                for i in range(4): area += float(h[i][0]*h[i+1][1] - h[i+1][0]*h[i][1])
                self.theta += self.kappa * (1 if area > 0.1 else (-1 if area < -0.1 else 0))

def build_events(tA, sA, oA, tB, sB, oB, ai, bi, shift, eng, flip):
    """Restores temporal alignment of products (prods)."""
    tB_s = tB - shift*1e-9
    times = np.concatenate([tA, tB_s])
    who = np.concatenate([np.zeros(len(tA),int), np.ones(len(tB),int)])
    idx = np.concatenate([np.arange(len(tA)), np.arange(len(tB))])
    order = np.argsort(times)
    mapB = -np.ones(len(tA), int); mapB[ai] = bi
    ph, ks, pr = [], [], []
    eng.reset()
    for w, i in zip(who[order], idx[order]):
        if w == 0:
            eng.update(0, int(sA[i]))
            j = mapB[i]
            if j >= 0:
                ph.append(eng.theta); ks.append(int(sA[i]*2 + sB[j]))
                pr.append(int(oA[i]) * int(oB[j]) * flip) # Correct Temporal Prod
        else: eng.update(1, int(sB[i]))
    return np.array(ph), np.array(ks), np.array(pr, dtype=np.int16)

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
    chi2 = np.sum(((S-X@beta)/err)**2)
    chi2_0 = np.sum(((S - np.mean(S))/err)**2)
    return chi2_0 - chi2, beta

def cv_align_score(phases, ks, prods, offsets, bins):
    ph_a, ks_a, pr_a = [], [], []
    for i in range(len(offsets)-1):
        s, e = offsets[i], offsets[i+1]
        mid = s + (e-s)//2
        _, beta = bin_and_score(phases[s:mid], ks[s:mid], prods[s:mid], bins)
        if beta is not None:
            phi = math.atan2(beta[2], beta[1])
            ph_a.append(phases[mid:e] + phi); ks_a.append(ks[mid:e]); pr_a.append(prods[mid:e])
    if not ph_a: return 0.0
    sc, _ = bin_and_score(np.concatenate(ph_a), np.concatenate(ks_a), np.concatenate(pr_a), bins)
    return sc

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alice-zip", required=True); ap.add_argument("--bob-zip", required=True)
    ap.add_argument("--path-contains", default="timetags/longdist")
    ap.add_argument("--window", type=float, default=1.0)
    ap.add_argument("--mode", type=str, default="strict")
    args = ap.parse_args()

    print(f"=== VERDICT v36_fixed: Window={args.window}ns | Mode={args.mode} ===")
    
    with zipfile.ZipFile(args.alice_zip,"r") as za, zipfile.ZipFile(args.bob_zip,"r") as zb:
        raw_names = [_norm_member(n) for n in za.namelist()]
        stems = sorted(list(set([n[:-6] for n in raw_names if args.path_contains.lower() in n.lower() and n.lower().endswith("_v.dat")])))
        
        eng = StrictThetaEngine(kappa=KAPPA)
        ph_l, k_l, p_l, off_l, proc_stems = [], [], [], [0], []
        
        for s in stems:
            tA, cA = load_raw_data(za, s); tB, cB = load_raw_data(zb, s)
            if tA is None: continue
            shift = estimate_shift_ns(tA, tB)
            if shift is None: continue # Skip unreliable sync
            ai, bi = match_pairs(tA, tB, shift, args.window, args.mode)
            if len(ai) < 500: continue
            
            sA, oA = decode_bits(cA, CFG["sb"], CFG["db"], CFG["invS_A"], False)
            sB, oB = decode_bits(cB, CFG["sb"], CFG["db"], CFG["invS_B"], False)
            ph, kk, pr = build_events(tA, sA, oA, tB, sB, oB, ai, bi, shift, eng, CFG["flip"])
            
            ph_l.append(ph); k_l.append(kk); p_l.append(pr); off_l.append(off_l[-1]+len(ph))
            proc_stems.append(s)

    if not ph_l: return print("No valid matches found.")
    phases, ks, prods, offsets = np.concatenate(ph_l), np.concatenate(k_l), np.concatenate(p_l), np.array(off_l)
    
    # GLOBAL S WITH GUARD
    c = np.bincount(ks, minlength=4)
    if np.any(c == 0): print("Warning: Zero counts in setting pairs."); c = np.maximum(c, 1)
    m = np.bincount(ks, weights=prods, minlength=4)/c
    print(f"\n[GLOBAL] S = {m[0]+m[1]+m[2]-m[3]:.4f} | Events: {len(phases)} | Stems: {len(proc_stems)}")
    
    real_score = cv_align_score(phases, ks, prods, offsets, BIN_COUNT)
    print(f"[TEST 1] Real Signal DeltaChi2 = {real_score:.2f}")
    
    print(f"[TEST 2] Permutation Test ({PERM_ROUNDS} rounds)...")
    better, rng = 0, np.random.default_rng(999)
    for r in range(PERM_ROUNDS):
        ph_p = phases.copy()
        for i in range(len(offsets)-1):
            s, e = offsets[i], offsets[i+1]
            # SHUFFLE INCLUDING REMAINDER CHUNKS
            for bs in range(s, e, BLOCK_SIZE):
                be = min(bs + BLOCK_SIZE, e)
                chunk = ph_p[bs:be]; rng.shuffle(chunk); ph_p[bs:be] = chunk
        if cv_align_score(ph_p, ks, prods, offsets, BIN_COUNT) >= real_score: better += 1
        if (r+1)%400==0: print(f"  Round {r+1}...")
            
    print(f"\n[RESULT] P-Value = {(better+1)/(PERM_ROUNDS+1):.4f}")

if __name__ == "__main__": main()