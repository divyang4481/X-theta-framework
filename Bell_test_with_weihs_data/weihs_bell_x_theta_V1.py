    #!/usr/bin/env python3

"""
weihs_bell_x_theta.py

Stream one "stem" at a time from the Weihs et al. Bell timetag dataset (Zenodo 7185335),
stored as Alice.zip and Bob.zip.

Files per stem:
- *_V.dat : big-endian float64 times (seconds)
- *_C.dat : big-endian uint16 codes (two low bits encode detector + setting)

This tool:
1) Builds one-to-one coincidences using:
     | tA - (tB + delta) | <= W
2) Computes:
   - E(a,b) and CHSH S
   - A simple X–θ "loop orientation" probe in setting-space (CW vs CCW square cycles)
   - Optional θ-conditioned CHSH (toy θ updated on loop closures)
3) Aggregates POOLED statistics across all stems for robust hypothesis testing.

PowerShell example:
  python .\weihs_bell_x_theta_V1.py `
    --alice-zip "..\Bell_data\7185335\Alice.zip" `
    --bob-zip   "..\Bell_data\7185335\Bob.zip" `
    --folder longdist `
    --window-ns 6 `
    --delta-auto `
    --reverse-bits `
    --max-stems 20 `
    --out "weihs_results_revbits_v1.json"
"""

from __future__ import annotations

import argparse
import json
import math
import zipfile
import random
from dataclasses import dataclass, field
from collections import deque
from typing import Dict, List, Tuple, Optional

import numpy as np


# -----------------------------
# Setting-space loop helpers
# -----------------------------
State = Tuple[int, int]  # (a,b), each in {0,1}

def hamming01(u: State, v: State) -> int:
    return (u[0] != v[0]) + (u[1] != v[1])

def signed_area(poly: List[State]) -> float:
    # poly is closed (last == first)
    area = 0.0
    for (x1, y1), (x2, y2) in zip(poly[:-1], poly[1:]):
        area += (x1 * y2 - x2 * y1)
    return 0.5 * area

def is_square_cycle(states: List[State]) -> bool:
    """
    A "square cycle" in the 2x2 setting lattice:
      - length 5 (closed)
      - 4 unique corners visited
      - each step flips exactly one of (a,b)
    """
    if len(states) != 5:
        return False
    if states[0] != states[-1]:
        return False
    if len(set(states[:-1])) != 4:
        return False
    for u, v in zip(states[:-1], states[1:]):
        if hamming01(u, v) != 1:
            return False
    return True


# -----------------------------
# Stats accumulators
# -----------------------------
@dataclass
class RunningAB:
    count: np.ndarray  # (4,)
    sum_ab: np.ndarray # (4,)

    def update(self, a: int, b: int, ab: int) -> None:
        idx = a * 2 + b
        self.count[idx] += 1
        self.sum_ab[idx] += ab

    def E(self) -> np.ndarray:
        E = np.zeros(4, dtype=float)
        nz = self.count > 0
        E[nz] = self.sum_ab[nz] / self.count[nz]
        return E

    def S(self) -> float:
        E = self.E()
        # Convention: S = E00 + E01 + E10 - E11
        return float(E[0] + E[1] + E[2] - E[3])


@dataclass
class ThetaBinned:
    bins: int
    count: np.ndarray  # (bins, 4)
    sum_ab: np.ndarray # (bins, 4)

    def update(self, bin_id: int, a: int, b: int, ab: int) -> None:
        idx = a * 2 + b
        self.count[bin_id, idx] += 1
        self.sum_ab[bin_id, idx] += ab

    def S_by_bin(self) -> Tuple[np.ndarray, np.ndarray]:
        S = np.zeros(self.bins, dtype=float)
        Nmin = np.zeros(self.bins, dtype=int)
        for k in range(self.bins):
            cnt = self.count[k]
            sab = self.sum_ab[k]
            E = np.zeros(4, dtype=float)
            nz = cnt > 0
            E[nz] = sab[nz] / cnt[nz]
            S[k] = E[0] + E[1] + E[2] - E[3]
            Nmin[k] = int(cnt.min())
        return S, Nmin

@dataclass
class PooledAnalysis:
    """Aggregates statistics across multiple stems for global hypothesis testing."""
    global_cw_vals: List[float] = field(default_factory=list)
    global_ccw_vals: List[float] = field(default_factory=list)
    
    # Global CHSH accumulators
    global_counts: np.ndarray = field(default_factory=lambda: np.zeros(4, dtype=np.int64))
    global_sums: np.ndarray = field(default_factory=lambda: np.zeros(4, dtype=np.int64))
    
    total_coincidences: int = 0

    def add_stem_result(self, raw_data: Dict[str, object]):
        # Loop test pooling
        self.global_cw_vals.extend(raw_data['cw_vals'])
        self.global_ccw_vals.extend(raw_data['ccw_vals'])
        
        # CHSH pooling
        self.global_counts += np.array(raw_data['counts'], dtype=np.int64)
        self.global_sums   += np.array(raw_data['sums'], dtype=np.int64)
        
        self.total_coincidences += raw_data.get('coincidences', 0)

    def report(self) -> Dict[str, object]:
        # 1. Global CHSH
        E = np.zeros(4, dtype=float)
        nz = self.global_counts > 0
        E[nz] = self.global_sums[nz] / self.global_counts[nz]
        S_global = E[0] + E[1] + E[2] - E[3]

        # 2. Global Loop Z-Score
        n_cw = len(self.global_cw_vals)
        n_ccw = len(self.global_ccw_vals)
        
        cw_mean = np.mean(self.global_cw_vals) if n_cw > 0 else 0.0
        ccw_mean = np.mean(self.global_ccw_vals) if n_ccw > 0 else 0.0
        
        # Pooled Standard Error
        # Var(diff) = Var(CW)/N_cw + Var(CCW)/N_ccw
        cw_var = np.var(self.global_cw_vals, ddof=1) if n_cw > 1 else 0.0
        ccw_var = np.var(self.global_ccw_vals, ddof=1) if n_ccw > 1 else 0.0
        
        se_diff = math.sqrt((cw_var / n_cw) + (ccw_var / n_ccw)) if (n_cw > 0 and n_ccw > 0) else 0.0
        
        diff = ccw_mean - cw_mean
        z_score = (diff / se_diff) if se_diff > 0 else 0.0

        print(f"\n====== POOLED ANALYSIS (All Stems) ======")
        print(f"Total Coincidences: {self.total_coincidences}")
        print(f"Global CHSH S:      {S_global:.5f}")
        print(f"-----------------------------------------")
        print(f"Loop Orientation Test (Pooled):")
        print(f"  CW cycles:  {n_cw} (mean {cw_mean:.4f})")
        print(f"  CCW cycles: {n_ccw} (mean {ccw_mean:.4f})")
        print(f"  Diff (CCW-CW): {diff:.5f}")
        print(f"  Z-Score:       {z_score:.5f}")
        print(f"=========================================\n")

        return {
            "global_S": S_global,
            "global_loop_z": z_score,
            "total_coincidences": self.total_coincidences,
            "n_cw": n_cw,
            "n_ccw": n_ccw,
            "diff": diff
        }


# -----------------------------
# Zip dataset parsing (Zenodo 7185335 format)
# -----------------------------
def _norm_member(name: str) -> str:
    """Normalize zip member names: use forward slashes, strip leading /, and strip leading 'Alice/' or 'Bob/'."""
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"):
        n = "/".join(parts[1:])
    return n

def _zip_name_map(z: zipfile.ZipFile) -> Dict[str, str]:
    """Map normalized lowercase -> original member name in zip."""
    m: Dict[str, str] = {}
    for name in z.namelist():
        key = _norm_member(name).lower()
        m[key] = name
    return m

def list_stems(z: zipfile.ZipFile, folder_filter: Optional[str]) -> List[str]:
    names = [_norm_member(n) for n in z.namelist()]
    if folder_filter:
        ff = folder_filter.replace("\\", "/").strip("/").lower()
        def ok(path: str) -> bool:
            p = path.lower()
            return p.startswith(ff + "/") or ("/" + ff + "/") in p
        names = [n for n in names if ok(n)]

    have_v, have_c = set(), set()
    for n in names:
        nl = n.lower()
        if nl.endswith("_v.dat"):
            have_v.add(n[:-6])
        elif nl.endswith("_c.dat"):
            have_c.add(n[:-6])

    return sorted(have_v.intersection(have_c))

def load_events(z: zipfile.ZipFile, stem: str, reverse_bits: bool) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    nm = _zip_name_map(z)
    v_key = (stem + "_V.dat").lower()
    c_key = (stem + "_C.dat").lower()
    if v_key not in nm or c_key not in nm:
        raise FileNotFoundError(f"Missing V/C files for stem: {stem}")

    with z.open(nm[v_key], "r") as f:
        vb = f.read()
    with z.open(nm[c_key], "r") as f:
        cb = f.read()

    t = np.frombuffer(vb, dtype=">f8")
    c = np.frombuffer(cb, dtype=">u2")

    n = min(len(t), len(c))
    t = t[:n].astype(np.float64, copy=False)
    c = c[:n]

    if not reverse_bits:
        det = (c & 1).astype(np.uint8, copy=False)
        s   = ((c >> 1) & 1).astype(np.uint8, copy=False)
    else:
        det = ((c >> 1) & 1).astype(np.uint8, copy=False)
        s   = (c & 1).astype(np.uint8, copy=False)

    o = (2 * det - 1).astype(np.int8, copy=False)
    return t, s, o


# -----------------------------
# Coincidence matching + X–θ tests
# -----------------------------
def count_coincidences(tA: np.ndarray, tB: np.ndarray, W: float, delta: float) -> int:
    tB2 = tB + delta
    i = j = 0
    nA = len(tA)
    nB = len(tB2)
    c = 0
    while i < nA and j < nB:
        ta = tA[i]
        tb = tB2[j]
        if tb < ta - W:
            j += 1
        elif tb > ta + W:
            i += 1
        else:
            c += 1
            i += 1
            j += 1
    return c

def auto_delta_ns(
    tA: np.ndarray,
    tB: np.ndarray,
    W: float,
    scan_range_ns: float = 50.0,
    coarse_step_ns: float = 1.0,
    fine_step_ns: float = 0.05,
    coarse_take: int = 200_000,
    fine_take: int = 400_000,
) -> float:
    tA1 = tA[:min(len(tA), coarse_take)]
    tB1 = tB[:min(len(tB), coarse_take)]

    deltas = np.arange(-scan_range_ns, scan_range_ns + 1e-12, coarse_step_ns) * 1e-9
    best_d = 0.0
    best_c = -1
    for d in deltas:
        c = count_coincidences(tA1, tB1, W, d)
        if c > best_c:
            best_c, best_d = c, d

    tA2 = tA[:min(len(tA), fine_take)]
    tB2 = tB[:min(len(tB), fine_take)]
    center_ns = best_d * 1e9
    fine_deltas = np.arange(center_ns - 2.0, center_ns + 2.0 + 1e-12, fine_step_ns) * 1e-9
    best_d2 = best_d
    best_c2 = -1
    for d in fine_deltas:
        c = count_coincidences(tA2, tB2, W, d)
        if c > best_c2:
            best_c2, best_d2 = c, d

    return float(best_d2 * 1e9)

def process_stem(
    tA: np.ndarray, sA: np.ndarray, oA: np.ndarray,
    tB: np.ndarray, sB: np.ndarray, oB: np.ndarray,
    W: float,
    delta_ns: float,
    theta_kappa: float,
    theta_bins: int
) -> Dict[str, object]:
    
    delta = delta_ns * 1e-9
    tB2 = tB + delta

    ab = RunningAB(count=np.zeros(4, dtype=np.int64), sum_ab=np.zeros(4, dtype=np.int64))
    tb = ThetaBinned(
        bins=theta_bins,
        count=np.zeros((theta_bins, 4), dtype=np.int64),
        sum_ab=np.zeros((theta_bins, 4), dtype=np.int64),
    )

    last_states: deque[State] = deque(maxlen=5)
    last_ab: deque[int] = deque(maxlen=5)
    cw_vals: List[float] = []
    ccw_vals: List[float] = []
    theta = 0.0

    i = j = 0
    nA = len(tA)
    nB = len(tB2)
    coincidences = 0

    while i < nA and j < nB:
        ta = tA[i]
        tbj = tB2[j]

        if tbj < ta - W:
            j += 1
            continue
        if tbj > ta + W:
            i += 1
            continue

        a = int(sA[i])
        b = int(sB[j])
        prod = int(oA[i] * oB[j])

        ab.update(a, b, prod)

        th = theta % (2.0 * math.pi)
        bin_id = int((th / (2.0 * math.pi)) * theta_bins)
        if bin_id >= theta_bins:
            bin_id = theta_bins - 1
        tb.update(bin_id, a, b, prod)

        last_states.append((a, b))
        last_ab.append(prod)
        if len(last_states) == 5 and is_square_cycle(list(last_states)):
            area = signed_area(list(last_states))
            orient = 1 if area > 0 else (-1 if area < 0 else 0)

            cycle_val = float(np.mean(list(last_ab)[1:]))
            if orient > 0:
                ccw_vals.append(cycle_val)
            elif orient < 0:
                cw_vals.append(cycle_val)

            if orient != 0 and theta_kappa != 0.0:
                theta += theta_kappa * orient

        coincidences += 1
        i += 1
        j += 1

    E = ab.E()
    S = ab.S()

    def mean_se(x: List[float]) -> Tuple[float, float, int]:
        if not x:
            return (float("nan"), float("nan"), 0)
        arr = np.asarray(x, dtype=float)
        m = float(arr.mean())
        if len(arr) < 2:
            return (m, float("nan"), len(arr))
        se = float(arr.std(ddof=1) / math.sqrt(len(arr)))
        return (m, se, len(arr))

    cw_m, cw_se, cw_n = mean_se(cw_vals)
    ccw_m, ccw_se, ccw_n = mean_se(ccw_vals)

    diff = (ccw_m - cw_m) if (cw_n and ccw_n) else float("nan")
    diff_se = (
        math.sqrt((0 if math.isnan(cw_se) else cw_se**2) + (0 if math.isnan(ccw_se) else ccw_se**2))
        if (cw_n > 1 and ccw_n > 1)
        else float("nan")
    )
    z = (diff / diff_se) if (diff_se and not math.isnan(diff_se) and diff_se > 0) else float("nan")

    S_bins, Nmin = tb.S_by_bin()

    return {
        "coincidences": int(coincidences),
        "E00_E01_E10_E11": [float(E[0]), float(E[1]), float(E[2]), float(E[3])],
        "CHSH_S": float(S),
        "delta_ns_used": float(delta_ns),
        "x_theta_loop_test": {
            "cw":  {"mean": cw_m,  "se": cw_se,  "n_cycles": cw_n},
            "ccw": {"mean": ccw_m, "se": ccw_se, "n_cycles": ccw_n},
            "diff_ccw_minus_cw": float(diff),
            "diff_se": float(diff_se),
            "z_score_approx": float(z),
        },
        "theta_conditioned": {
            "theta_kappa": float(theta_kappa),
            "theta_bins": int(theta_bins),
            "S_by_bin": [float(x) for x in S_bins.tolist()],
            "min_counts_by_bin": [int(x) for x in Nmin.tolist()],
            "S_bin_min": float(np.nanmin(S_bins)) if len(S_bins) else float("nan"),
            "S_bin_max": float(np.nanmax(S_bins)) if len(S_bins) else float("nan"),
        },
        # RAW DATA for Pooling
        "raw_data": {
            "cw_vals": cw_vals,
            "ccw_vals": ccw_vals,
            "counts": ab.count.tolist(),
            "sums": ab.sum_ab.tolist(),
            "coincidences": coincidences
        }
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--alice-zip", required=True)
    ap.add_argument("--bob-zip", required=True)
    ap.add_argument("--folder", default=None, help="Filter stems to those under this folder, e.g. longdist")
    ap.add_argument("--max-stems", type=int, default=20, help="Process at most this many matching stems (for speed)")

    ap.add_argument("--window-ns", type=float, default=6.0, help="Coincidence half-window W in ns")

    ap.add_argument("--delta-ns", type=float, default=None, help="Fixed delta shift applied to Bob times in ns")
    ap.add_argument("--delta-auto", action="store_true", help="Auto-estimate delta per stem (two-stage scan)")
    ap.add_argument("--delta-scan-range-ns", type=float, default=50.0)

    ap.add_argument("--reverse-bits", action="store_true", help="Swap low two bits in _C.dat decoding")
    ap.add_argument("--theta-kappa", type=float, default=0.25, help="Theta increment per detected loop (radians)")
    ap.add_argument("--theta-bins", type=int, default=16)

    # NEW: Null test
    ap.add_argument("--null-shuffle", action="store_true", help="Randomly shuffle Bob's outcomes to destroy correlations (Null Test)")

    ap.add_argument("--out", default="weihs_x_theta_results.json")
    ap.add_argument("--list-only", action="store_true", help="Only list matching stems and exit")
    args = ap.parse_args()

    W = args.window_ns * 1e-9

    with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
        stems_a = set(list_stems(za, args.folder))
        stems_b = set(list_stems(zb, args.folder))
        stems = sorted(list(stems_a.intersection(stems_b)))

        if not stems:
            raise SystemExit("No matching stems found between Alice.zip and Bob.zip (check --folder).")

        print(f"Found {len(stems)} matching stems.")
        print("First 20 stems:")
        for s in stems[:20]:
            print("  ", s)

        if args.list_only:
            return

        stems = stems[:args.max_stems]
        results = {
            "alice_zip": args.alice_zip,
            "bob_zip": args.bob_zip,
            "folder_filter": args.folder,
            "window_ns": args.window_ns,
            "reverse_bits": bool(args.reverse_bits),
            "null_shuffle": bool(args.null_shuffle),
            "processed_stems": [],
            "global_stats": None 
        }

        # Initialize Pooled Analysis
        pooled_stats = PooledAnalysis()

        for stem in stems:
            print(f"\n--- Processing stem: {stem}")
            tA, sA, oA = load_events(za, stem, reverse_bits=args.reverse_bits)
            tB, sB, oB = load_events(zb, stem, reverse_bits=args.reverse_bits)

            # Optional Null Test: Shuffle Bob's outcomes
            if args.null_shuffle:
                # We shuffle oB in place to destroy correlations with A
                np.random.shuffle(oB)
                # Note: We do NOT shuffle timestamps, so coincidence counts remain valid, 
                # but the spin correlations should vanish.

            if args.delta_ns is not None:
                delta_ns = float(args.delta_ns)
            elif args.delta_auto:
                delta_ns = auto_delta_ns(
                    tA, tB, W,
                    scan_range_ns=float(args.delta_scan_range_ns)
                )
            else:
                delta_ns = 0.0

            r = process_stem(
                tA, sA, oA,
                tB, sB, oB,
                W=W,
                delta_ns=delta_ns,
                theta_kappa=float(args.theta_kappa),
                theta_bins=int(args.theta_bins),
            )
            r["stem"] = stem
            
            # Feed data to pooled analysis
            pooled_stats.add_stem_result(r["raw_data"])
            
            # Remove raw data from per-stem JSON output to keep file size reasonable
            del r["raw_data"]
            
            results["processed_stems"].append(r)

            E = r["E00_E01_E10_E11"]
            print(f"coincidences={r['coincidences']}")
            print("E00,E01,E10,E11 =", E)
            print(f"CHSH S = {r['CHSH_S']:.6f}")
            lt = r["x_theta_loop_test"]
            print("loop diff (CCW-CW) =", lt["diff_ccw_minus_cw"], "z≈", lt["z_score_approx"])
            print("delta_ns_used =", r["delta_ns_used"])

        # Final Report
        global_report = pooled_stats.report()
        results["global_stats"] = global_report

        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)

        print("\nWrote:", args.out)


if __name__ == "__main__":
    main()