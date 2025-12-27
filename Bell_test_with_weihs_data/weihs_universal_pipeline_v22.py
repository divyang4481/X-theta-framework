#!/usr/bin/env python3
r"""
weihs_universal_pipeline_v22.py  (CV Align + cache backward-compat)

Fixes included:
- Backward compatible cache load:
    If cache lacks 'offsets', we auto-ignore it and re-mine.
- Case-insensitive stem/path filter fix (Timetags vs timetags).
- Cross-validated per-stem phase offset alignment (held-out halves).
- Stratified permutation test (shuffle theta within each setting-pair k).
- Filenames now include kappa and bin_count parameters.

Usage:
  python weihs_universal_pipeline_v22.py --alice-zip "..\Bell_data\7185335\Alice.zip" --bob-zip "..\Bell_data\7185335\Bob.zip" --path-contains "timetags/longdist" --inspect-c

  python weihs_universal_pipeline_v22.py --alice-zip "..\Bell_data\7185335\Alice.zip" --bob-zip "..\Bell_data\7185335\Bob.zip" --path-contains "timetags/longdist" --kappa 0.8 --bin-count 12 --perm-rounds 2000


  python weihs_universal_pipeline_v22.py   --alice-zip "..\Bell_data\7185335\Alice.zip"   --bob-zip "..\Bell_data\7185335\Bob.zip"   --path-contains "timetags/longdist"   --kappa 0.8 --bin-count 12 --perm-rounds 2000   --stratified-perm --save-plots   --cache weihs_mined_v22.npz


  python weihs_universal_pipeline_v22.py --alice-zip "..\Bell_data\7185335\Alice.zip" --bob-zip "..\Bell_data\7185335\Bob.zip" --path-contains "timetags/longdist" --force-remine --inspect-c


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
# Zip path normalization
# ============================================================

def _norm_member(name: str) -> str:
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"):
        n = "/".join(parts[1:])
    return n

def list_stems(z: zipfile.ZipFile, path_contains: Optional[str] = None) -> List[str]:
    """
    Return LOWERCASED stems where both *_V.dat and *_C.dat exist.
    path_contains is case-insensitive.
    """
    names = [_norm_member(n).replace("\\", "/") for n in z.namelist()]
    if path_contains:
        pc = path_contains.replace("\\", "/").lower()
        names = [n for n in names if pc in n.lower()]

    nmset = {n.lower() for n in names}

    stems = set()
    for n in names:
        nl = n.lower()
        if nl.endswith("_v.dat"):
            stems.add(nl[:-6])  # strip "_v.dat"

    out = []
    for s in stems:
        if (s + "_v.dat") in nmset and (s + "_c.dat") in nmset:
            out.append(s)
    out.sort()
    return out

# ============================================================
# Load Weihs raw files
# ============================================================

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

    t = np.frombuffer(vb, dtype=">f8")  # big-endian float64 seconds
    c = np.frombuffer(cb, dtype=">u2")  # big-endian uint16
    n = min(len(t), len(c))
    return t[:n].astype(np.float64, copy=False), c[:n].astype(np.uint16, copy=False)

# ============================================================
# C.dat inspection
# ============================================================

def inspect_c_bits(c: np.ndarray, top_bits: int = 10) -> Dict:
    uniq = np.unique(c)
    bit_p = []
    for b in range(16):
        bit_p.append(float(np.mean((c >> b) & 1)))

    informative = [b for b, p in enumerate(bit_p) if 0.05 < p < 0.95]
    informative = sorted(informative, key=lambda b: abs(bit_p[b] - 0.5))[:top_bits]
    return {
        "unique_count": int(len(uniq)),
        "unique_min": int(uniq.min()) if len(uniq) else 0,
        "unique_max": int(uniq.max()) if len(uniq) else 0,
        "bit_p": bit_p,
        "candidate_bits": informative,
    }

# ============================================================
# Decode (setting, outcome) from c
# ============================================================

def decode_from_c(c: np.ndarray, set_bit: int, det_bit: int,
                  invert_setting: bool = False,
                  invert_outcome: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    s = ((c >> set_bit) & 1).astype(np.uint8)
    det = ((c >> det_bit) & 1).astype(np.uint8)
    if invert_setting:
        s = (1 - s).astype(np.uint8)

    o = (2 * det - 1).astype(np.int8)  # 0->-1, 1->+1
    if invert_outcome:
        o = (-o).astype(np.int8)
    return s, o

# ============================================================
# Shift estimation (nearest-neighbor dt sampling)
# ============================================================

def estimate_shift_ns(tA: np.ndarray, tB: np.ndarray,
                      window_ns: float = 200.0,
                      bin_ns: float = 0.5,
                      sample: int = 20000,
                      seed: int = 1) -> Tuple[float, float]:
    rng = np.random.default_rng(seed)
    nA = len(tA); nB = len(tB)
    if nA < 200 or nB < 200:
        return 0.0, 0.0

    m = min(sample, nA)
    idxA = rng.choice(nA, size=m, replace=False)
    a = tA[idxA]

    j = np.searchsorted(tB, a)
    j0 = np.clip(j - 1, 0, nB - 1)
    j1 = np.clip(j,     0, nB - 1)

    dt0 = (tB[j0] - a) * 1e9
    dt1 = (tB[j1] - a) * 1e9
    dt = np.where(np.abs(dt0) <= np.abs(dt1), dt0, dt1)

    mask = np.abs(dt) <= window_ns
    dt = dt[mask]
    if dt.size < 200:
        return 0.0, 0.0

    bins = int((2 * window_ns) / bin_ns)
    hist, edges = np.histogram(dt, bins=bins, range=(-window_ns, window_ns))
    peak = int(np.argmax(hist))
    peak_val = float(hist[peak])
    med = float(np.median(hist)) + 1e-9
    snr = peak_val / med
    peak_center = 0.5 * (edges[peak] + edges[peak + 1])
    return float(peak_center), float(snr)

# ============================================================
# One-to-one coincidence matching (greedy)
# ============================================================

def match_coincidences(tA: np.ndarray, tB: np.ndarray,
                       shift_ns: float,
                       window_ns: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Align Bob by subtracting shift (dt = tB - tA)
    tB_aligned = tB - (shift_ns * 1e-9)

    nB = len(tB_aligned)
    j = np.searchsorted(tB_aligned, tA)

    j0 = np.clip(j - 1, 0, nB - 1)
    j1 = np.clip(j,     0, nB - 1)

    dt0 = (tB_aligned[j0] - tA) * 1e9
    dt1 = (tB_aligned[j1] - tA) * 1e9
    choose0 = np.abs(dt0) <= np.abs(dt1)

    bj = np.where(choose0, j0, j1).astype(np.int64)
    dt = np.where(choose0, dt0, dt1).astype(np.float64)

    ok = np.abs(dt) <= window_ns
    a_idx = np.nonzero(ok)[0].astype(np.int64)
    b_idx = bj[ok]
    dt_ns = dt[ok]

    if a_idx.size == 0:
        return a_idx, b_idx, dt_ns

    order = np.argsort(np.abs(dt_ns))
    usedB = np.zeros(nB, dtype=bool)
    keep_a, keep_b, keep_dt = [], [], []
    for q in order:
        bi = int(b_idx[q])
        if not usedB[bi]:
            usedB[bi] = True
            keep_a.append(int(a_idx[q]))
            keep_b.append(bi)
            keep_dt.append(float(dt_ns[q]))

    return (np.array(keep_a, dtype=np.int64),
            np.array(keep_b, dtype=np.int64),
            np.array(keep_dt, dtype=np.float64))

# ============================================================
# CHSH computation
# ============================================================

def compute_chsh_from_pairs(sA: np.ndarray, oA: np.ndarray,
                            sB: np.ndarray, oB: np.ndarray,
                            a_idx: np.ndarray, b_idx: np.ndarray,
                            flip_product: int = 1) -> Tuple[float, np.ndarray, np.ndarray]:
    a = sA[a_idx].astype(np.int8)
    b = sB[b_idx].astype(np.int8)
    prod = (oA[a_idx].astype(np.int16) * oB[b_idx].astype(np.int16)) * int(flip_product)
    k = (2 * a + b).astype(np.int8)  # 0..3

    counts = np.bincount(k, minlength=4).astype(np.int64)
    sums = np.bincount(k, weights=prod.astype(np.float64), minlength=4).astype(np.float64)

    means = np.zeros(4, dtype=np.float64)
    ok = counts > 0
    means[ok] = sums[ok] / counts[ok]

    if np.any(counts == 0):
        return float("nan"), counts, means

    S = means[0] + means[1] + means[2] - means[3]
    return float(S), counts, means

# ============================================================
# Theta engine (holonomy rule)
# ============================================================

State = Tuple[int, int]

def get_orientation(states: List[State]) -> int:
    area = 0.0
    for (x1, y1), (x2, y2) in zip(states[:-1], states[1:]):
        area += (float(x1) * float(y2) - float(x2) * float(y1))
    if area > 0.1:
        return 1
    if area < -0.1:
        return -1
    return 0

@dataclass
class ThetaEngine:
    kappa: float
    a: int = 0
    b: int = 0
    theta: float = 0.0
    history: deque = field(default_factory=lambda: deque(maxlen=5))

    def reset(self):
        self.a = 0
        self.b = 0
        self.theta = 0.0
        self.history.clear()

    def update_setting(self, who: int, setting: int):
        if who == 0:
            self.a = int(setting)
        else:
            self.b = int(setting)
        self.history.append((self.a, self.b))
        if len(self.history) == 5:
            h = list(self.history)
            if h[0] == h[-1] and len(set(h[:-1])) == 4:
                self.theta += self.kappa * get_orientation(h)

def build_phase_events_per_stem(tA: np.ndarray, sA: np.ndarray,
                                tB: np.ndarray, sB: np.ndarray,
                                a_idx: np.ndarray, b_idx: np.ndarray,
                                prod: np.ndarray,
                                shift_ns: float,
                                engine: ThetaEngine) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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

    return (np.array(phases, dtype=np.float64),
            np.array(ks, dtype=np.int8),
            np.array(prods, dtype=np.int8))

# ============================================================
# Phase-binned CHSH + sine fit
# ============================================================

def phase_binned_chsh(phases: np.ndarray, ks: np.ndarray, prods: np.ndarray,
                      bin_count: int,
                      min_counts_per_k: int = 50) -> Optional[Dict]:
    if phases.size < 5000:
        return None

    two_pi = 2 * math.pi
    b = ((phases % two_pi) / two_pi * bin_count).astype(np.int64) % bin_count
    flat = b * 4 + ks.astype(np.int64)

    counts = np.bincount(flat, minlength=bin_count * 4).astype(np.int64).reshape((bin_count, 4))
    sums = np.bincount(flat, weights=prods.astype(np.float64), minlength=bin_count * 4).astype(np.float64).reshape((bin_count, 4))

    valid_bins = np.all(counts >= min_counts_per_k, axis=1)
    if np.sum(valid_bins) < max(4, bin_count // 2):
        return None

    means = np.zeros_like(sums, dtype=np.float64)
    means[valid_bins] = sums[valid_bins] / counts[valid_bins]

    var = 1.0 - np.clip(means, -1.0, 1.0) ** 2
    se = np.sqrt(np.maximum(var, 1e-12) / np.maximum(counts, 1))

    s_vals = means[:, 0] + means[:, 1] + means[:, 2] - means[:, 3]
    s_errs = np.sqrt(se[:, 0]**2 + se[:, 1]**2 + se[:, 2]**2 + se[:, 3]**2)

    idx = np.where(valid_bins)[0]
    centers = (idx + 0.5) * (two_pi / bin_count)

    return {
        "bin_centers": centers.astype(np.float64),
        "s_vals": s_vals[valid_bins].astype(np.float64),
        "s_errs": s_errs[valid_bins].astype(np.float64),
        "valid_bin_idx": idx.astype(np.int64),
    }

def weighted_fit_sine(bin_centers: np.ndarray, s_vals: np.ndarray, s_errs: np.ndarray) -> Dict:
    th = bin_centers
    y = s_vals
    sig = np.maximum(s_errs, 1e-12)
    w = 1.0 / (sig * sig)
    sw = np.sqrt(w)

    X = np.column_stack([np.ones_like(th), np.sin(th), np.cos(th)])
    Xw = X * sw[:, None]
    yw = y * sw

    beta = np.linalg.lstsq(Xw, yw, rcond=None)[0]
    yhat = X @ beta
    chi2 = float(np.sum(((y - yhat) / sig) ** 2))

    X0 = np.ones((len(th), 1))
    X0w = X0 * sw[:, None]
    beta0 = np.linalg.lstsq(X0w, yw, rcond=None)[0]
    yhat0 = X0 @ beta0
    chi2_0 = float(np.sum(((y - yhat0) / sig) ** 2))

    ybar = float(np.sum(w * y) / np.sum(w))
    sst = float(np.sum(w * (y - ybar) ** 2))
    sse = float(np.sum(w * (y - yhat) ** 2))
    r2w = 1.0 - (sse / sst if sst > 0 else 0.0)

    return {
        "beta": beta,
        "chi2": chi2,
        "chi2_flat": chi2_0,
        "delta_chi2": chi2_0 - chi2,
        "r2w": r2w,
    }

# ============================================================
# Stratified permutation test (correct null)
# ============================================================

def permutation_test_stratified(phases: np.ndarray, ks: np.ndarray, prods: np.ndarray,
                                bin_count: int, rounds: int,
                                min_counts_per_k: int,
                                seed: int = 123) -> Dict:
    rng = np.random.default_rng(seed)

    obs = phase_binned_chsh(phases, ks, prods, bin_count, min_counts_per_k=min_counts_per_k)
    if obs is None:
        return {"ok": False, "reason": "insufficient bins"}

    fit = weighted_fit_sine(obs["bin_centers"], obs["s_vals"], obs["s_errs"])
    obs_delta = float(fit["delta_chi2"])

    better = 0
    used = 0

    idx_by_k = [np.where(ks == k)[0] for k in range(4)]

    for _ in range(rounds):
        ph_perm = phases.copy()
        for k in range(4):
            idx = idx_by_k[k]
            if idx.size > 1:
                rng.shuffle(ph_perm[idx])

        pres = phase_binned_chsh(ph_perm, ks, prods, bin_count, min_counts_per_k=min_counts_per_k)
        if pres is None:
            continue

        fitp = weighted_fit_sine(pres["bin_centers"], pres["s_vals"], pres["s_errs"])
        d = float(fitp["delta_chi2"])

        used += 1
        if d >= obs_delta:
            better += 1

    used = max(1, used)
    p = (better + 1) / (used + 1)
    return {"ok": True, "obs_delta_chi2": obs_delta, "p_value": float(p), "rounds_used": int(used), "better": int(better)}

# ============================================================
# Cross-validated per-stem phase alignment (θ0 nuisance)
# ============================================================

def crossval_align_heldout(phases: np.ndarray, ks: np.ndarray, prods: np.ndarray,
                           offsets: np.ndarray,
                           bin_count: int,
                           min_counts_per_k: int) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    phases_out = []
    ks_out = []
    prods_out = []

    for i in range(len(offsets) - 1):
        lo = int(offsets[i])
        hi = int(offsets[i + 1])
        n = hi - lo
        if n < 8000:
            continue

        mid = lo + n // 2
        ph_tr = phases[lo:mid]
        ks_tr = ks[lo:mid]
        pr_tr = prods[lo:mid]

        obs_tr = phase_binned_chsh(ph_tr, ks_tr, pr_tr, bin_count, min_counts_per_k=min_counts_per_k)
        if obs_tr is None:
            continue

        fit_tr = weighted_fit_sine(obs_tr["bin_centers"], obs_tr["s_vals"], obs_tr["s_errs"])
        _, a, b = fit_tr["beta"]
        amp = float(math.sqrt(a*a + b*b))
        if amp < 1e-6:
            continue

        phi = float(math.atan2(b, a))

        ph_te = phases[mid:hi] + phi
        ks_te = ks[mid:hi]
        pr_te = prods[mid:hi]

        phases_out.append(ph_te.astype(np.float64, copy=False))
        ks_out.append(ks_te.astype(np.int8, copy=False))
        prods_out.append(pr_te.astype(np.int8, copy=False))

    if not phases_out:
        return None

    return (np.concatenate(phases_out),
            np.concatenate(ks_out),
            np.concatenate(prods_out))

# ============================================================
# Calibration search
# ============================================================

def calibrate_config_on_files(
    za: zipfile.ZipFile, zb: zipfile.ZipFile,
    stems_lc: List[str],
    max_files: int,
    sync_window_ns: float,
    sync_bin_ns: float,
    match_window_ns: float,
    required_s: float,
    verbose: bool = True
) -> Dict:
    tested = stems_lc[:max_files]
    if len(tested) == 0:
        raise RuntimeError("No stems to calibrate.")

    _, cA0 = load_raw_data(za, tested[0])
    infoA = inspect_c_bits(cA0, top_bits=8)
    cand_bits = infoA["candidate_bits"]
    if len(cand_bits) < 2:
        cand_bits = [0, 1]

    bit_pairs = []
    for sb in cand_bits:
        for db in cand_bits:
            if db != sb:
                bit_pairs.append((sb, db))
    if not bit_pairs:
        bit_pairs = [(0, 1), (1, 0)]

    toggles = [(False, False), (True, False), (False, True), (True, True)]
    best = {"score": -1e18}

    for (set_bit, det_bit) in bit_pairs:
        for invS_A, invO_A in toggles:
            for invS_B, invO_B in toggles:
                for flip_prod in (1, -1):
                    Ss = []
                    ok_files = 0
                    for j, stem in enumerate(tested):
                        tA, cA = load_raw_data(za, stem)
                        tB, cB = load_raw_data(zb, stem)
                        sA, oA = decode_from_c(cA, set_bit, det_bit, invS_A, invO_A)
                        sB, oB = decode_from_c(cB, set_bit, det_bit, invS_B, invO_B)

                        shift_ns, snr = estimate_shift_ns(tA, tB, window_ns=sync_window_ns, bin_ns=sync_bin_ns, seed=11 + j)
                        if snr < 5.0:
                            continue

                        a_idx, b_idx, _ = match_coincidences(tA, tB, shift_ns=shift_ns, window_ns=match_window_ns)
                        if a_idx.size < 500:
                            continue

                        S, counts, _ = compute_chsh_from_pairs(sA, oA, sB, oB, a_idx, b_idx, flip_product=flip_prod)
                        if not np.isfinite(S):
                            continue
                        if np.any(counts < 50):
                            continue

                        Ss.append(S)
                        ok_files += 1

                    if ok_files == 0:
                        continue

                    medS = float(np.median(Ss))
                    score = medS + 0.05 * ok_files
                    if score > best["score"]:
                        best = {
                            "score": score,
                            "medS": medS,
                            "ok_files": ok_files,
                            "set_bit": set_bit,
                            "det_bit": det_bit,
                            "invS_A": invS_A,
                            "invO_A": invO_A,
                            "invS_B": invS_B,
                            "invO_B": invO_B,
                            "flip_prod": flip_prod,
                            "candidate_bits": cand_bits,
                        }
                        if verbose:
                            print(f"  New best: medS={medS:.3f} ok_files={ok_files} cfg={(set_bit, det_bit, invS_A, invO_A, invS_B, invO_B, flip_prod)}")

    if "set_bit" not in best or best.get("medS", -999) < required_s:
        print("\nCRITICAL: Calibration could not find a config with strong Bell violation.")
    else:
        print(f"\nCALIBRATION WINNER: medS={best['medS']:.3f} using cfg={(best['set_bit'], best['det_bit'], best['invS_A'], best['invO_A'], best['invS_B'], best['invO_B'], best['flip_prod'])}")
    return best

# ============================================================
# Plot helper
# ============================================================

def save_fit_plot(outfile: str, centers: np.ndarray, svals: np.ndarray, serrs: np.ndarray, beta: np.ndarray, title: str):
    c0, a, b = beta
    th = np.linspace(0, 2 * math.pi, 400)
    y = c0 + a * np.sin(th) + b * np.cos(th)

    plt.figure(figsize=(10, 6))
    plt.errorbar(centers, svals, yerr=serrs, fmt='o')
    plt.plot(th, y)
    plt.axhline(2.0, linestyle="--")
    plt.axhline(2.828, linestyle=":")
    plt.title(title)
    plt.xlabel("theta (rad)")
    plt.ylabel("S(theta)")
    plt.grid(True, alpha=0.3)
    plt.savefig(outfile, dpi=160)
    plt.close()

# ============================================================
# MAIN
# ============================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alice-zip", required=True)
    ap.add_argument("--bob-zip", required=True)
    ap.add_argument("--path-contains", default="", help="Case-insensitive substring filter on normalized zip paths.")
    ap.add_argument("--inspect-c", action="store_true")

    ap.add_argument("--kappa", type=float, default=0.8)
    ap.add_argument("--bin-count", type=int, default=12)
    ap.add_argument("--perm-rounds", type=int, default=2000)
    ap.add_argument("--min-counts-per-k", type=int, default=50)


    # Compatibility toggles (older commands may pass these flags)
    ap.add_argument("--save-plots", action="store_true",
                    help="(Compatibility) Save plots (default behavior is to save unless --no-plots).")
    ap.add_argument("--no-plots", action="store_true",
                    help="Disable saving plots to disk.")
    ap.add_argument("--stratified-perm", action="store_true",
                    help="(Compatibility) Use stratified permutation within setting-pair k (default).")
    ap.add_argument("--no-stratified-perm", action="store_true",
                    help="Disable stratified permutation (NOT recommended; mostly for debugging).")
    ap.add_argument("--cache", default="weihs_mined_v22.npz")
    ap.add_argument("--calib-files", type=int, default=6)
    ap.add_argument("--required-s", type=float, default=2.0)

    ap.add_argument("--sync-window-ns", type=float, default=200.0)
    ap.add_argument("--sync-bin-ns", type=float, default=0.5)
    ap.add_argument("--coinc-window-ns", type=float, default=1.5)

    ap.add_argument("--max-files", type=int, default=999999)
    ap.add_argument("--disable-cv-align", action="store_true")
    ap.add_argument("--force-remine", action="store_true", help="Ignore cache and remine.")
    args = ap.parse_args()
    # Compatibility flags:
    if getattr(args, "save_plots", False):
        args.no_plots = False
    if getattr(args, "stratified_perm", False):
        args.no_stratified_perm = False


    path_contains = args.path_contains.strip() or None

    print("==========================================")
    print("   WEIHS UNIVERSAL PIPELINE v22 (CV Align)")
    print("==========================================")

    with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
        stemsA = set(list_stems(za, path_contains=path_contains))
        stemsB = set(list_stems(zb, path_contains=path_contains))
        stems = sorted(list(stemsA.intersection(stemsB)))
        if args.max_files < len(stems):
            stems = stems[:args.max_files]

        if not stems:
            raise RuntimeError("No common stems found. Check --path-contains filter.")

        print(f"Found {len(stems)} common stems. (filter={path_contains})")
        print("First 5 stems:", stems[:5])

        if args.inspect_c:
            _, cA = load_raw_data(za, stems[0])
            info = inspect_c_bits(cA, top_bits=12)
            print("\n[C INSPECT]")
            print(f"  stem={stems[0]}")
            print(f"  unique_count={info['unique_count']} min={info['unique_min']} max={info['unique_max']}")
            print("  candidate_bits:", info["candidate_bits"])
            print("  per-bit p(1):")
            for b, p in enumerate(info["bit_p"]):
                print(f"    bit{b:02d}: {p:.4f}")
            return

        # -------------- CACHE LOAD (backward compatible) --------------
        phases_all = ks_all = prods_all = offsets = None
        cfg = {}
        cache_ok = False

        if (not args.force_remine) and os.path.exists(args.cache):
            print(f"[CACHE] Loading {args.cache} ...")
            try:
                data = np.load(args.cache, allow_pickle=True)
                files = set(data.files)
                # New format requires offsets
                if {"phases", "ks", "prods"}.issubset(files) and ("offsets" in files):
                    phases_all = data["phases"].astype(np.float64)
                    ks_all = data["ks"].astype(np.int8)
                    prods_all = data["prods"].astype(np.int8)
                    offsets = data["offsets"].astype(np.int64)
                    if "cfg" in files:
                        try:
                            cfg = data["cfg"].item()
                        except Exception:
                            cfg = {}
                    cache_ok = True
                    print(f"[CACHE] Loaded N={len(phases_all)} phase-events from {len(offsets)-1} stems.")
                else:
                    print("[CACHE] Cache is old/missing 'offsets'. Ignoring cache and re-mining to enable CV-align.")
                    cache_ok = False
            except Exception as e:
                print(f"[CACHE] Failed to load cache ({e}). Will re-mine.")
                cache_ok = False

        # -------------- MINE if no usable cache --------------
        if not cache_ok:
            print("[CALIBRATION] Finding decode config that yields Bell violation...")
            cfg = calibrate_config_on_files(
                za, zb, stems,
                max_files=min(args.calib_files, len(stems)),
                sync_window_ns=args.sync_window_ns,
                sync_bin_ns=args.sync_bin_ns,
                match_window_ns=args.coinc_window_ns,
                required_s=args.required_s,
                verbose=True
            )

            set_bit = int(cfg.get("set_bit", 0))
            det_bit = int(cfg.get("det_bit", 1))
            invS_A = bool(cfg.get("invS_A", False))
            invO_A = bool(cfg.get("invO_A", False))
            invS_B = bool(cfg.get("invS_B", False))
            invO_B = bool(cfg.get("invO_B", False))
            flip_prod = int(cfg.get("flip_prod", 1))

            engine = ThetaEngine(kappa=float(args.kappa))

            phases_chunks = []
            ks_chunks = []
            prods_chunks = []
            offsets_list = [0]

            accepted = 0
            print("\n[MINING] Mining stems with per-stem CHSH sanity gate...")

            for i, stem in enumerate(stems):
                tA, cA = load_raw_data(za, stem)
                tB, cB = load_raw_data(zb, stem)
                sA, oA = decode_from_c(cA, set_bit, det_bit, invS_A, invO_A)
                sB, oB = decode_from_c(cB, set_bit, det_bit, invS_B, invO_B)

                shift_ns, snr = estimate_shift_ns(tA, tB, window_ns=args.sync_window_ns, bin_ns=args.sync_bin_ns, seed=1000 + i)
                if snr < 5.0:
                    continue

                a_idx, b_idx, _ = match_coincidences(tA, tB, shift_ns=shift_ns, window_ns=args.coinc_window_ns)
                if a_idx.size < 1000:
                    continue

                S, counts, _ = compute_chsh_from_pairs(sA, oA, sB, oB, a_idx, b_idx, flip_product=flip_prod)
                if not np.isfinite(S):
                    continue
                if np.any(counts < 100):
                    continue
                if S < (args.required_s - 0.1):
                    continue

                prod = (oA[a_idx].astype(np.int16) * oB[b_idx].astype(np.int16)) * flip_prod
                prod = prod.astype(np.int8)

                ph, kk, pp = build_phase_events_per_stem(
                    tA, sA, tB, sB,
                    a_idx, b_idx,
                    prod,
                    shift_ns=shift_ns,
                    engine=engine
                )

                if ph.size < 5000:
                    continue

                phases_chunks.append(ph)
                ks_chunks.append(kk)
                prods_chunks.append(pp)

                offsets_list.append(offsets_list[-1] + ph.size)
                accepted += 1

                if accepted % 5 == 0:
                    print(f"  accepted={accepted}  lastS={S:.3f}  stem={stem}")

            if accepted == 0:
                raise RuntimeError("No stems passed CHSH sanity. Still not matching photons correctly (or wrong decode bits).")

            phases_all = np.concatenate(phases_chunks).astype(np.float64)
            ks_all = np.concatenate(ks_chunks).astype(np.int8)
            prods_all = np.concatenate(prods_chunks).astype(np.int8)
            offsets = np.array(offsets_list, dtype=np.int64)

            np.savez(args.cache, phases=phases_all, ks=ks_all, prods=prods_all, offsets=offsets, cfg=cfg)
            print(f"\n[MINING] Done. accepted_files={accepted}, N_phase_events={len(phases_all)}")
            print(f"[CACHE] Saved: {args.cache}")

        # ---------------------------
        # Unaligned analysis
        # ---------------------------
        print("\n[ANALYSIS-A] Unaligned global phase-binned CHSH...")
        obsA = phase_binned_chsh(phases_all, ks_all, prods_all, bin_count=args.bin_count, min_counts_per_k=args.min_counts_per_k)
        if obsA is None:
            raise RuntimeError("Unaligned: insufficient bins. Try fewer bins or lower min-count thresholds.")

        fitA = weighted_fit_sine(obsA["bin_centers"], obsA["s_vals"], obsA["s_errs"])
        ampA = float(math.sqrt(fitA["beta"][1]**2 + fitA["beta"][2]**2))
        print(f"  delta_chi2={fitA['delta_chi2']:.3f}  R2w={fitA['r2w']:.3f}  amp={ampA:.6f}")
        
        # FIX: Include kappa and bin_count in filename
        out_png_unaligned = f"weihs_v22_unaligned_k{args.kappa}_b{args.bin_count}.png"
    if not args.no_plots:
        save_fit_plot(out_png_unaligned, obsA["thetas"], obsA["s_vals"], obsA["s_errs"], fitA["beta"],
                      title=f"Unaligned (bins={args.bin_count})  Δχ²={fitA['delta_chi2']:.2f}",
                      show_refs=True)
        ks_perm = ks_all if not args.no_stratified_perm else np.zeros_like(ks_all)
        permA = permutation_test_stratified(phases_all, ks_perm, prods_all, args.bin_count, args.perm_rounds,
                                            min_counts_per_k=args.min_counts_per_k, seed=777)
        if permA["ok"]:
            print(f"  perm p={permA['p_value']:.6f}  (better={permA['better']}/{permA['rounds_used']})")
        else:
            print("  permutation failed:", permA.get("reason", "?"))

        # ---------------------------
        # CV-aligned held-out analysis
        # ---------------------------
        if args.disable_cv_align:
            print("\n[ANALYSIS-B] CV alignment disabled (--disable-cv-align). Done.")
            return

        print("\n[ANALYSIS-B] Cross-validated per-stem phase alignment (held-out halves)...")
        cv = crossval_align_heldout(phases_all, ks_all, prods_all, offsets,
                                    bin_count=args.bin_count, min_counts_per_k=args.min_counts_per_k)
        if cv is None:
            raise RuntimeError("CV alignment produced no usable stems (too small / insufficient bins).")

        ph_cv, ks_cv, pr_cv = cv
        print(f"  CV pooled held-out N={len(ph_cv)}")

        obsB = phase_binned_chsh(ph_cv, ks_cv, pr_cv, bin_count=args.bin_count, min_counts_per_k=args.min_counts_per_k)
        if obsB is None:
            raise RuntimeError("CV-aligned: insufficient bins. Try fewer bins or lower min-count thresholds.")

        fitB = weighted_fit_sine(obsB["bin_centers"], obsB["s_vals"], obsB["s_errs"])
        ampB = float(math.sqrt(fitB["beta"][1]**2 + fitB["beta"][2]**2))
        print(f"  delta_chi2={fitB['delta_chi2']:.3f}  R2w={fitB['r2w']:.3f}  amp={ampB:.6f}")
        
        # FIX: Include kappa and bin_count in filename
        out_png_cvaligned = f"weihs_v22_cvaligned_k{args.kappa}_b{args.bin_count}.png"
    if not args.no_plots:
        save_fit_plot(out_png_cvaligned, obsB["thetas"], obsB["s_vals"], obsB["s_errs"], fitB["beta"],
                      title=f"CV-aligned held-out (bins={args.bin_count})  Δχ²={fitB['delta_chi2']:.2f}",
                      show_refs=True)
        ks_permB = ks_cv if not args.no_stratified_perm else np.zeros_like(ks_cv)
        permB = permutation_test_stratified(ph_cv, ks_permB, pr_cv, args.bin_count, args.perm_rounds,
                                            min_counts_per_k=args.min_counts_per_k, seed=999)
        if permB["ok"]:
            print(f"  perm p={permB['p_value']:.6f}  (better={permB['better']}/{permB['rounds_used']})")
        else:
            print("  permutation failed:", permB.get("reason", "?"))
    if not args.no_plots:
        print(f"\nSaved plots: {out_png_unaligned}, {out_png_cvaligned}")
    else:
        print("\nPlots not saved (--no-plots).")
    print("DONE.")

if __name__ == "__main__":
    main()
