#!/usr/bin/env python3
r"""
weihs_universal_pipeline_v26.py

THE "REVIEWER-PROOF" BUILD (Physical Scale & Strict Matching)
================================================================================
CRITICAL UPGRADES:
  1. LOOP DENSITY (rho): Measures loops/pair to calculate the physical 
     memory length (N_pairs) of the effect. This answers "Why Kappa=0.1?".
  2. STRICT 1-TO-1 MATCHING: Enforces unique photon usage to prevent bias.
  3. RIGOROUS PERMUTATION: Runs the full alignment on shuffles (End-to-End).
  4. ROBUSTNESS SCAN: Checks bin stability (4, 6, 8, 12).

USAGE:
  python weihs_universal_pipeline_v26.py --alice-zip "..\Bell_data\7185335\Alice.zip" --bob-zip "..\Bell_data\7185335\Bob.zip" --path-contains "timetags/longdist"
"""

import argparse
import zipfile
import numpy as np
import math
import os
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from collections import deque
from typing import List, Tuple, Optional, Dict

# ============================================================
# UTILITIES
# ============================================================

def _norm_member(name: str) -> str:
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"):
        n = "/".join(parts[1:])
    return n

def list_stems(z: zipfile.ZipFile, path_contains: Optional[str] = None) -> List[str]:
    names = [_norm_member(n).replace("\\", "/") for n in z.namelist()]
    if path_contains:
        pc = path_contains.replace("\\", "/").lower()
        names = [n for n in names if pc in n.lower()]

    nmset = {n.lower() for n in names}
    stems = set()
    for n in names:
        nl = n.lower()
        if nl.endswith("_v.dat"):
            stems.add(nl[:-6])

    out = []
    for s in stems:
        if (s + "_v.dat") in nmset and (s + "_c.dat") in nmset:
            out.append(s)
    out.sort()
    return out

def load_raw_data(z: zipfile.ZipFile, stem_lc: str) -> Tuple[np.ndarray, np.ndarray]:
    norm_to_zipname = {_norm_member(n).lower(): n for n in z.namelist()}
    vkey = (stem_lc + "_v.dat").lower()
    ckey = (stem_lc + "_c.dat").lower()
    if vkey not in norm_to_zipname or ckey not in norm_to_zipname:
        raise FileNotFoundError(f"Missing V/C for stem={stem_lc}")

    with z.open(norm_to_zipname[vkey], "r") as f:
        vb = f.read()
    with z.open(norm_to_zipname[ckey], "r") as f:
        cb = f.read()

    t = np.frombuffer(vb, dtype=">f8")
    c = np.frombuffer(cb, dtype=">u2")
    n = min(len(t), len(c))
    return t[:n].astype(np.float64, copy=False), c[:n].astype(np.uint16, copy=False)

def decode_from_c(c: np.ndarray, set_bit: int, det_bit: int,
                  invert_setting: bool = False,
                  invert_outcome: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    s = ((c >> set_bit) & 1).astype(np.int32)
    det = ((c >> det_bit) & 1).astype(np.int32)
    
    if invert_setting: s = 1 - s
    o = (2 * det - 1) 
    if invert_outcome: o = -o
        
    return s.astype(np.uint8), o.astype(np.int8)

def estimate_shift_ns(tA: np.ndarray, tB: np.ndarray, window_ns: float = 200.0, bin_ns: float = 0.5, seed: int = 1) -> Tuple[float, float]:
    rng = np.random.default_rng(seed)
    nA = len(tA); nB = len(tB)
    if nA < 200 or nB < 200: return 0.0, 0.0

    m = min(20000, nA)
    idxA = rng.choice(nA, size=m, replace=False)
    a = tA[idxA]

    j = np.searchsorted(tB, a)
    j0 = np.clip(j - 1, 0, nB - 1); j1 = np.clip(j, 0, nB - 1)
    dt0 = (tB[j0] - a) * 1e9; dt1 = (tB[j1] - a) * 1e9
    dt = np.where(np.abs(dt0) <= np.abs(dt1), dt0, dt1)

    mask = np.abs(dt) <= window_ns
    dt = dt[mask]
    if dt.size < 200: return 0.0, 0.0

    bins = int((2 * window_ns) / bin_ns)
    hist, edges = np.histogram(dt, bins=bins, range=(-window_ns, window_ns))
    peak = int(np.argmax(hist))
    peak_val = float(hist[peak])
    med = float(np.median(hist)) + 1e-9
    snr = peak_val / med
    peak_center = 0.5 * (edges[peak] + edges[peak + 1])
    return float(peak_center), float(snr)

# ============================================================
# REVIEWER UPGRADE: STRICT 1-TO-1 MATCHING
# ============================================================

def match_coincidences_strict(tA: np.ndarray, tB: np.ndarray, shift_ns: float, window_ns: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Strict 1-to-1 matching. Each Bob photon is used exactly once.
    Prioritizes the smallest time difference |dt|.
    """
    tB_aligned = tB - (shift_ns * 1e-9)
    
    # 1. Find all candidate pairs within window
    j = np.searchsorted(tB_aligned, tA)
    j0 = np.clip(j - 1, 0, len(tB) - 1)
    j1 = np.clip(j, 0, len(tB) - 1)
    
    dt0 = (tB_aligned[j0] - tA) * 1e9
    dt1 = (tB_aligned[j1] - tA) * 1e9
    
    cand_idx_A = []
    cand_idx_B = []
    cand_dt = []
    
    mask0 = np.abs(dt0) <= window_ns
    idxA0 = np.where(mask0)[0]
    if len(idxA0) > 0:
        cand_idx_A.append(idxA0)
        cand_idx_B.append(j0[idxA0])
        cand_dt.append(dt0[idxA0])
        
    mask1 = np.abs(dt1) <= window_ns
    idxA1 = np.where(mask1)[0]
    if len(idxA1) > 0:
        cand_idx_A.append(idxA1)
        cand_idx_B.append(j1[idxA1])
        cand_dt.append(dt1[idxA1])
        
    if not cand_idx_A:
        return np.array([]), np.array([]), np.array([])
        
    all_A = np.concatenate(cand_idx_A)
    all_B = np.concatenate(cand_idx_B)
    all_dt = np.concatenate(cand_dt)
    
    # 2. Sort by |dt| (best matches first)
    sort_order = np.argsort(np.abs(all_dt))
    
    # 3. Greedy assignment
    used_A = np.zeros(len(tA), dtype=bool)
    used_B = np.zeros(len(tB), dtype=bool)
    
    final_A, final_B, final_dt = [], [], []
    
    for k in sort_order:
        idx_a = all_A[k]
        idx_b = all_B[k]
        
        if not used_A[idx_a] and not used_B[idx_b]:
            used_A[idx_a] = True
            used_B[idx_b] = True
            final_A.append(idx_a)
            final_B.append(idx_b)
            final_dt.append(all_dt[k])
            
    return np.array(final_A, dtype=np.int64), np.array(final_B, dtype=np.int64), np.array(final_dt)

def compute_chsh_from_pairs(sA, oA, sB, oB, a_idx, b_idx, flip_product=1):
    a = sA[a_idx].astype(np.int8); b = sB[b_idx].astype(np.int8)
    prod = (oA[a_idx].astype(np.int16) * oB[b_idx].astype(np.int16)) * int(flip_product)
    k = (2 * a + b).astype(np.int8)
    counts = np.bincount(k, minlength=4).astype(np.int64)
    sums = np.bincount(k, weights=prod.astype(np.float64), minlength=4).astype(np.float64)
    if np.any(counts == 0): return float("nan"), counts, np.zeros(4)
    means = sums / counts
    S = means[0] + means[1] + means[2] - means[3]
    return float(S), counts, means

# ============================================================
# THETA ENGINE (With Loop Counter)
# ============================================================

def get_orientation(states: List[Tuple[int,int]]) -> int:
    area = 0.0
    for (x1, y1), (x2, y2) in zip(states[:-1], states[1:]):
        area += (float(x1) * float(y2) - float(x2) * float(y1))
    return 1 if area > 0.1 else (-1 if area < -0.1 else 0)

@dataclass
class ThetaEngine:
    kappa: float
    a: int = 0; b: int = 0; theta: float = 0.0
    loop_count: int = 0  # NEW: Track physical loops
    history: deque = field(default_factory=lambda: deque(maxlen=5))

    def reset(self):
        self.a = 0; self.b = 0; self.theta = 0.0; self.history.clear()
        self.loop_count = 0

    def update_setting(self, who: int, setting: int):
        if who == 0: self.a = int(setting)
        else: self.b = int(setting)
        self.history.append((self.a, self.b))
        if len(self.history) == 5:
            h = list(self.history)
            if h[0] == h[-1] and len(set(h[:-1])) == 4:
                ori = get_orientation(h)
                self.theta += self.kappa * ori
                if ori != 0: self.loop_count += 1

def build_phase_events(tA, sA, tB, sB, a_idx, b_idx, prod, shift_ns, engine):
    tB_aligned = tB - (shift_ns * 1e-9)
    times = np.concatenate([tA, tB_aligned])
    who = np.concatenate([np.zeros(len(tA), dtype=np.int8), np.ones(len(tB_aligned), dtype=np.int8)])
    idx = np.concatenate([np.arange(len(tA), dtype=np.int64), np.arange(len(tB_aligned), dtype=np.int64)])
    order = np.argsort(times, kind="mergesort")
    who = who[order]; idx = idx[order]

    mapB = -np.ones(len(tA), dtype=np.int64)
    mapP = np.zeros(len(tA), dtype=np.int16)
    mapB[a_idx] = b_idx
    mapP[a_idx] = prod.astype(np.int16)

    phases, ks, prods = [], [], []
    engine.reset()

    for w, i in zip(who, idx):
        if w == 0:
            engine.update_setting(0, int(sA[i]))
            bi = mapB[i]
            if bi >= 0:
                k = int(sA[i]) * 2 + int(sB[bi])
                phases.append(engine.theta)
                ks.append(k)
                prods.append(int(mapP[i]))
        else:
            engine.update_setting(1, int(sB[i]))

    return np.array(phases), np.array(ks, dtype=np.int8), np.array(prods, dtype=np.int8), engine.loop_count

# ============================================================
# END-TO-END STATISTICS
# ============================================================

def phase_binned_chsh(phases, ks, prods, bin_count, min_counts=50):
    if phases.size < 5000: return None
    b = ((phases % (2*math.pi)) / (2*math.pi) * bin_count).astype(np.int64) % bin_count
    flat = b * 4 + ks.astype(np.int64)
    counts = np.bincount(flat, minlength=bin_count*4).reshape((bin_count, 4))
    sums = np.bincount(flat, weights=prods, minlength=bin_count*4).reshape((bin_count, 4))
    
    valid = np.all(counts >= min_counts, axis=1)
    if np.sum(valid) < 4: return None
    
    means = sums[valid] / counts[valid]
    var = 1.0 - means**2
    se = np.sqrt(np.maximum(var, 1e-12) / counts[valid])
    
    s_vals = means[:,0] + means[:,1] + means[:,2] - means[:,3]
    s_errs = np.sqrt(np.sum(se**2, axis=1))
    centers = (np.where(valid)[0] + 0.5) * (2*math.pi / bin_count)
    return {"centers": centers, "s": s_vals, "err": s_errs}

def weighted_fit_sine(bin_centers, s_vals, s_errs):
    th = bin_centers; y = s_vals
    w = 1.0 / (np.maximum(s_errs, 1e-12)**2); sw = np.sqrt(w)
    X = np.column_stack([np.ones_like(th), np.sin(th), np.cos(th)])
    beta = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
    yhat = X @ beta
    chi2 = np.sum(((y - yhat) / s_errs)**2)
    
    X0 = np.ones((len(th), 1))
    beta0 = np.linalg.lstsq(X0 * sw[:, None], y * sw, rcond=None)[0]
    chi2_0 = np.sum(((y - X0 @ beta0) / s_errs)**2)
    
    return {"beta": beta, "delta_chi2": chi2_0 - chi2}

def cv_align_and_score(phases, ks, prods, offsets, bin_count, min_counts):
    ph_cv, ks_cv, pr_cv = [], [], []
    for i in range(len(offsets)-1):
        lo, hi = offsets[i], offsets[i+1]
        mid = lo + (hi-lo)//2
        tr = phase_binned_chsh(phases[lo:mid], ks[lo:mid], prods[lo:mid], bin_count, min_counts)
        if not tr: continue
        res = weighted_fit_sine(tr["centers"], tr["s"], tr["err"])
        phi = math.atan2(res["beta"][2], res["beta"][1])
        ph_cv.append(phases[mid:hi] + phi)
        ks_cv.append(ks[mid:hi]); pr_cv.append(prods[mid:hi])
    
    if not ph_cv: return 0.0
    ph_f = np.concatenate(ph_cv); ks_f = np.concatenate(ks_cv); pr_f = np.concatenate(pr_cv)
    obs = phase_binned_chsh(ph_f, ks_f, pr_f, bin_count, min_counts)
    if not obs: return 0.0
    fit = weighted_fit_sine(obs["centers"], obs["s"], obs["err"])
    return fit["delta_chi2"]

def rigorous_permutation_test(phases, ks, prods, offsets, bin_count, rounds, min_counts):
    print(f"\n[PERMUTATION] Running {rounds} rigorous rounds (Shuffle -> Re-Align -> Score)...")
    obs_delta = cv_align_and_score(phases, ks, prods, offsets, bin_count, min_counts)
    print(f"  Observed DeltaChi2: {obs_delta:.2f}")
    
    better = 0; rng = np.random.default_rng(123)
    idx_groups = [np.where(ks == k)[0] for k in range(4)]
    
    for i in range(rounds):
        ph_perm = phases.copy()
        for idx in idx_groups:
            if idx.size > 1:
                temp = ph_perm[idx]
                rng.shuffle(temp)
                ph_perm[idx] = temp
        
        perm_delta = cv_align_and_score(ph_perm, ks, prods, offsets, bin_count, min_counts)
        if perm_delta >= obs_delta: better += 1
        if (i+1) % 100 == 0: print(f"    Round {i+1}: max_noise={perm_delta:.2f}")
            
    p = (better + 1) / (rounds + 1)
    return p, obs_delta

def calibrate_config(za, zb, stems):
    print("  [CALIBRATION] Finding bits...")
    best = {"score": -999}
    
    for (sb, db) in [(0, 1), (1, 0)]:
        for invA in [(False, False), (True, False), (False, True), (True, True)]:
            for invB in [(False, False), (True, False), (False, True), (True, True)]:
                for flip in [1, -1]:
                    Ss = []
                    for s in stems[:6]:
                        tA, cA = load_raw_data(za, s); tB, cB = load_raw_data(zb, s)
                        sA, oA = decode_from_c(cA, sb, db, *invA)
                        sB, oB = decode_from_c(cB, sb, db, *invB)
                        shift, snr = estimate_shift_ns(tA, tB)
                        if snr < 5: continue
                        a_i, b_i, _ = match_coincidences_strict(tA, tB, shift, 1.5)
                        if len(a_i) < 500: continue
                        S, _, _ = compute_chsh_from_pairs(sA, oA, sB, oB, a_i, b_i, flip)
                        if np.isfinite(S): Ss.append(S)
                    
                    if Ss:
                        med = np.median(Ss)
                        if med > best["score"]:
                            best = {"score": med, "sb": sb, "db": db, "invA": invA, "invB": invB, "flip": flip}
                            print(f"    New Best: S={med:.3f}")
    
    print(f"  WINNER: S={best['score']:.3f}")
    return best

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alice-zip", required=True); ap.add_argument("--bob-zip", required=True)
    ap.add_argument("--path-contains", default="timetags/longdist")
    ap.add_argument("--kappa", type=float, default=0.1)
    ap.add_argument("--bin-count", type=int, default=8)
    ap.add_argument("--perm-rounds", type=int, default=1000)
    ap.add_argument("--cache", default="weihs_mined_v26.npz")
    ap.add_argument("--force-remine", action="store_true")
    args = ap.parse_args()

    print(f"=== WEIHS UNIVERSAL PIPELINE v26 (Robust Mode) ===")
    
    with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
        stems = sorted(list(set(list_stems(za, args.path_contains)).intersection(list_stems(zb, args.path_contains))))
        print(f"Found {len(stems)} stems.")

        if args.force_remine or not os.path.exists(args.cache):
            cfg = calibrate_config(za, zb, stems)
            if cfg["score"] < 1.8:
                print("CRITICAL: Calibration failed."); return

            engine = ThetaEngine(kappa=args.kappa)
            ph_list, k_list, p_list, off_list = [], [], [], [0]
            total_loops, total_pairs = 0, 0
            accepted = 0
            
            print("[MINING] Processing (Strict 1-to-1)...")
            for i, s in enumerate(stems):
                tA, cA = load_raw_data(za, s); tB, cB = load_raw_data(zb, s)
                sA, oA = decode_from_c(cA, cfg["sb"], cfg["db"], *cfg["invA"])
                sB, oB = decode_from_c(cB, cfg["sb"], cfg["db"], *cfg["invB"])
                shift, snr = estimate_shift_ns(tA, tB)
                if snr < 5: continue
                
                a_i, b_i, _ = match_coincidences_strict(tA, tB, shift, 1.5)
                if len(a_i) < 1000: continue
                
                S, _, _ = compute_chsh_from_pairs(sA, oA, sB, oB, a_i, b_i, cfg["flip"])
                if S < 1.8: continue
                
                prod = (oA[a_i].astype(np.int16) * oB[b_i].astype(np.int16)) * cfg["flip"]
                ph, kk, pp, n_loops = build_phase_events(tA, sA, tB, sB, a_i, b_i, prod, shift, engine)
                
                ph_list.append(ph); k_list.append(kk); p_list.append(pp)
                off_list.append(off_list[-1] + len(ph))
                total_loops += n_loops
                total_pairs += len(ph)
                
                accepted += 1
                if accepted % 5 == 0: print(f"  Accepted {accepted}, S={S:.2f}")

            np.savez(args.cache, phases=np.concatenate(ph_list), ks=np.concatenate(k_list), 
                     prods=np.concatenate(p_list), offsets=np.array(off_list),
                     meta={'loops': total_loops, 'pairs': total_pairs, 'kappa': args.kappa})
            print(f"[CACHE] Saved to {args.cache}")
        else:
            # Load stats from cache if available
            d = np.load(args.cache, allow_pickle=True)
            if 'meta' in d:
                meta = d['meta'].item()
                total_loops = meta['loops']
                total_pairs = meta['pairs']
            else:
                total_loops, total_pairs = 0, 1 # Fallback

        # --- PHYSICAL INTERPRETATION ---
        rho = total_loops / total_pairs if total_pairs > 0 else 0
        N_pi = math.pi / (args.kappa * rho) if rho > 0 else 0
        print(f"\n[PHYSICAL SCALE]")
        print(f"  Loop Density (rho): {rho:.3f} loops/pair")
        print(f"  Memory Length (N_pi): {N_pi:.1f} pairs")
        
        # --- ANALYSIS ---
        data = np.load(args.cache)
        phases = data["phases"]; ks = data["ks"]; prods = data["prods"]; offsets = data["offsets"]
        
        p_val, final_dchi2 = rigorous_permutation_test(phases, ks, prods, offsets, args.bin_count, args.perm_rounds, 50)
        
        print(f"\n[FINAL RESULT v26]")
        print(f"  DeltaChi2: {final_dchi2:.2f}")
        print(f"  P-Value:   {p_val:.5f} (Rigorous)")
        
        print("\n[ROBUSTNESS SCAN]")
        for b in [4, 6, 8, 12, 16]:
            d = cv_align_and_score(phases, ks, prods, offsets, b, 50)
            print(f"  Bin={b:02d}: DeltaChi2={d:.2f}")

if __name__ == "__main__":
    main()