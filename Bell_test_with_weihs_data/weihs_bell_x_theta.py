#!/usr/bin/env python3
"""
weihs_bell_x_theta.py

Weihs et al. Bell dataset (Zenodo 7185335): Alice.zip + Bob.zip. :contentReference[oaicite:3]{index=3}

Reads per-stem triplets:
  *_V.dat : big-endian float64 times (seconds)
  *_C.dat : big-endian uint16 code bits:
            LSB = detector (0/1), next bit = switch setting (0/1). :contentReference[oaicite:4]{index=4}

Coincidences: |tA - (tB + delta)| <= W

Computes:
- E(a,b) 2x2
- CHSH S (Weihs form): |E(a0,b0)-E(a1,b0)| + |E(a0,b1)+E(a1,b1)| :contentReference[oaicite:5]{index=5}
- S_max over allowed relabelings (row/col sign flips + swaps), so conventions can’t hide violation
- X–θ loop-orientation toy test in setting-space (CW vs CCW square cycles)
- CHSH conditioned on a toy θ updated on loop closure

PowerShell example:
  python .\weihs_bell_x_theta.py `
    --alice-zip "..\Bell_data\7185335\Alice.zip" `
    --bob-zip   "..\Bell_data\7185335\Bob.zip" `
    --folder longdist `
    --window-ns 6 `
    --delta-auto `
    --max-stems 10 `
    --out weihs_results.json
"""

from __future__ import annotations

import argparse
import json
import math
import zipfile
from dataclasses import dataclass
from collections import deque
from typing import Dict, List, Tuple, Optional

import numpy as np


# ------------------------- JSON helpers -------------------------
def jfloat(x: float) -> Optional[float]:
    if x is None:
        return None
    if isinstance(x, (np.floating,)):
        x = float(x)
    if not isinstance(x, (float, int)):
        return None
    if math.isnan(x) or math.isinf(x):
        return None
    return float(x)

def jint(x: int) -> int:
    return int(x)


# ------------------------- Setting-space loop helpers -------------------------
def hamming01(u: Tuple[int, int], v: Tuple[int, int]) -> int:
    return (u[0] != v[0]) + (u[1] != v[1])

def signed_area(poly: List[Tuple[int, int]]) -> float:
    # poly is closed (last == first)
    area = 0.0
    for (x1, y1), (x2, y2) in zip(poly[:-1], poly[1:]):
        area += (x1 * y2 - x2 * y1)
    return 0.5 * area

def is_square_cycle(states: List[Tuple[int, int]]) -> bool:
    # Expect 5 states: 4 unique corners + return to start
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


# ------------------------- Stats accumulators -------------------------
@dataclass
class RunningAB:
    # index = a*2 + b
    count: np.ndarray  # (4,)
    sum_ab: np.ndarray # (4,)

    def update(self, a: int, b: int, ab: int) -> None:
        idx = a * 2 + b
        self.count[idx] += 1
        self.sum_ab[idx] += ab

    def E4(self) -> np.ndarray:
        E = np.full(4, np.nan, dtype=float)
        nz = self.count > 0
        E[nz] = self.sum_ab[nz] / self.count[nz]
        return E

    def E2x2(self) -> np.ndarray:
        E4 = self.E4()
        return E4.reshape(2, 2)

@dataclass
class ThetaBinned:
    bins: int
    count: np.ndarray  # (bins,4)
    sum_ab: np.ndarray # (bins,4)

    def update(self, bin_id: int, a: int, b: int, ab: int) -> None:
        idx = a * 2 + b
        self.count[bin_id, idx] += 1
        self.sum_ab[bin_id, idx] += ab

    def S_by_bin(self, min_pair_count: int = 0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Returns:
          S[bins] with NaN for bins that fail the min_pair_count threshold
          min_counts[bins] = min over 4 setting-pairs in that bin
        """
        S = np.full(self.bins, np.nan, dtype=float)
        Nmin = np.zeros(self.bins, dtype=int)

        for k in range(self.bins):
            cnt = self.count[k]
            Nmin[k] = int(cnt.min())
            if Nmin[k] < min_pair_count:
                continue
            sab = self.sum_ab[k]
            E = np.full(4, np.nan, dtype=float)
            nz = cnt > 0
            E[nz] = sab[nz] / cnt[nz]
            # Weihs CHSH form needs E00,E01,E10,E11; per-bin we store the 4 E's
            E00, E01, E10, E11 = E[0], E[1], E[2], E[3]
            # Use Weihs Eq.(1) structure with default ordering
            if not (np.isnan(E00) or np.isnan(E10) or np.isnan(E01) or np.isnan(E11)):
                S[k] = abs(E00 - E10) + abs(E01 + E11)

        return S, Nmin


# ------------------------- Zip index + parsing -------------------------
def _norm_member(name: str) -> str:
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    # Some zips contain top-level "Alice/" or "Bob/" folder; normalize away.
    if parts and parts[0].lower() in ("alice", "bob"):
        n = "/".join(parts[1:])
    return n

def build_zip_index(z: zipfile.ZipFile) -> Dict[str, str]:
    """
    Map normalized lowercase path -> original zip member name.
    """
    m: Dict[str, str] = {}
    for name in z.namelist():
        key = _norm_member(name).lower()
        m[key] = name
    return m

def list_stems(z: zipfile.ZipFile, zidx: Dict[str, str], folder_filter: Optional[str]) -> List[str]:
    keys = list(zidx.keys())  # normalized lowercase
    if folder_filter:
        ff = folder_filter.replace("\\", "/").strip("/").lower()
        # match either ".../<ff>/..." or prefix "<ff>/..."
        keys = [k for k in keys if (f"/{ff}/" in k) or k.startswith(ff + "/")]

    have_v, have_c = set(), set()
    for k in keys:
        if k.endswith("_v.dat"):
            have_v.add(k[:-6])  # drop "_v.dat"
        elif k.endswith("_c.dat"):
            have_c.add(k[:-6])  # drop "_c.dat"

    # stems are stored as normalized strings (lowercase paths)
    return sorted(have_v.intersection(have_c))

def decode_C(c: np.ndarray, reverse_bits: bool) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Input:
      c: uint16 array
    Output:
      det: uint8 0/1
      s:   uint8 0/1 (switch)
      o:   int8  ±1 (outcome mapped from det)
    """
    if reverse_bits:
        # swap the two lowest bits (LSB <-> next LSB)
        det = ((c >> 1) & 1).astype(np.uint8, copy=False)
        s   = (c & 1).astype(np.uint8, copy=False)
    else:
        det = (c & 1).astype(np.uint8, copy=False)
        s   = ((c >> 1) & 1).astype(np.uint8, copy=False)

    o = (2 * det - 1).astype(np.int8, copy=False)  # 0->-1, 1->+1
    return det, s, o

def load_events(
    z: zipfile.ZipFile,
    zidx: Dict[str, str],
    stem_lower: str,
    reverse_bits: bool
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    stem_lower: normalized lowercase stem without _V/_C suffix, e.g. "timetags/longdist/longdist0"
    Returns:
      t: float64 seconds (ascending)
      s: uint8 (0/1 setting)
      o: int8 outcome ±1
    """
    v_key = (stem_lower + "_v.dat")
    c_key = (stem_lower + "_c.dat")
    if v_key not in zidx or c_key not in zidx:
        raise FileNotFoundError(f"Missing V/C files for stem: {stem_lower}")

    with z.open(zidx[v_key], "r") as f:
        vb = f.read()
    with z.open(zidx[c_key], "r") as f:
        cb = f.read()

    t = np.frombuffer(vb, dtype=">f8")  # big-endian float64
    c = np.frombuffer(cb, dtype=">u2")  # big-endian uint16

    n = min(len(t), len(c))
    t = t[:n].astype(np.float64, copy=False)
    c = c[:n]

    _, s, o = decode_C(c, reverse_bits=reverse_bits)
    return t, s, o


# ------------------------- Delta estimation (scan for coincidence peak) -------------------------
def count_coincidences_fast(tA: np.ndarray, tB: np.ndarray, W: float, delta: float) -> int:
    """
    Fast approximate count (not enforcing one-to-one):
      For each tA, checks nearest neighbor in shifted tB (searchsorted).
    Used only for picking delta.
    """
    b = tB + delta
    idx = np.searchsorted(b, tA, side="left")

    # check idx and idx-1
    nB = len(b)
    best = np.full(len(tA), np.inf, dtype=float)

    m = idx < nB
    if np.any(m):
        best[m] = np.minimum(best[m], np.abs(b[idx[m]] - tA[m]))

    m2 = idx > 0
    if np.any(m2):
        best[m2] = np.minimum(best[m2], np.abs(b[idx[m2] - 1] - tA[m2]))

    return int(np.sum(best <= W))

def estimate_delta_by_scan(
    tA: np.ndarray,
    tB: np.ndarray,
    W: float,
    scan_range_ns: float,
    scan_step_ns: float,
    max_events: int,
    refine_rounds: int = 2
) -> float:
    """
    Returns delta (seconds) that maximizes coincidence peak height on a subset.
    """
    # Subsample to limit cost
    tA0 = tA[:max_events]
    tB0 = tB[:max_events]

    # Light thinning to reduce searchsorted cost if arrays huge
    if len(tA0) > 200_000:
        tA0 = tA0[::2]
    if len(tB0) > 200_000:
        tB0 = tB0[::2]

    center = 0.0
    step = scan_step_ns * 1e-9
    rng = scan_range_ns * 1e-9

    best_delta = 0.0
    best_cnt = -1

    for _ in range(max(1, refine_rounds)):
        deltas = np.arange(center - rng, center + rng + 0.5 * step, step, dtype=float)
        for d in deltas:
            cnt = count_coincidences_fast(tA0, tB0, W, d)
            if cnt > best_cnt:
                best_cnt = cnt
                best_delta = float(d)

        # refine around best
        center = best_delta
        rng = max(rng / 5.0, 5e-9)   # don't shrink below 5 ns span
        step = max(step / 5.0, 0.05e-9)  # don't go below 0.05 ns step

    return best_delta


# ------------------------- CHSH (Weihs form + maxima over relabelings) -------------------------
def S_weihs(E: np.ndarray) -> float:
    """
    Weihs Eq.(1) with default mapping:
      S = |E00 - E10| + |E01 + E11|
    where E[a,b] with a in {0,1}, b in {0,1}
    """
    E00, E01 = E[0, 0], E[0, 1]
    E10, E11 = E[1, 0], E[1, 1]
    if np.isnan(E00) or np.isnan(E01) or np.isnan(E10) or np.isnan(E11):
        return float("nan")
    return float(abs(E00 - E10) + abs(E01 + E11))

def S_weihs_max_relabel(E: np.ndarray) -> Tuple[float, Dict[str, object]]:
    """
    Max over:
      - swapping Alice settings (a0<->a1)
      - swapping Bob settings (b0<->b1)
      - flipping outcome labels per setting (row/col sign flips)
    This prevents “wrong convention” from hiding violation.
    """
    best = -1.0
    best_cfg: Dict[str, object] = {}

    # swaps
    for swap_a in (False, True):
        for swap_b in (False, True):
            E2 = E.copy()
            if swap_a:
                E2 = E2[[1, 0], :]
            if swap_b:
                E2 = E2[:, [1, 0]]

            # row/col sign flips: outcome relabeling can be per setting
            for ra0 in (1.0, -1.0):
                for ra1 in (1.0, -1.0):
                    for cb0 in (1.0, -1.0):
                        for cb1 in (1.0, -1.0):
                            Em = E2.copy()
                            Em[0, :] *= ra0
                            Em[1, :] *= ra1
                            Em[:, 0] *= cb0
                            Em[:, 1] *= cb1
                            s = S_weihs(Em)
                            if not np.isnan(s) and s > best:
                                best = float(s)
                                best_cfg = {
                                    "swap_a": swap_a,
                                    "swap_b": swap_b,
                                    "row_flips": [ra0, ra1],
                                    "col_flips": [cb0, cb1],
                                }

    return best, best_cfg


# ------------------------- Main coincidence processing + X–θ toys -------------------------
def process_stem(
    tA: np.ndarray, sA: np.ndarray, oA: np.ndarray,
    tB: np.ndarray, sB: np.ndarray, oB: np.ndarray,
    W: float,
    delta: float,
    theta_kappa: float,
    theta_bins: int,
    theta_min_pair_count: int,
    keep_cycles: bool
) -> Dict[str, object]:
    """
    Greedy one-to-one coincidence pairing.
    """
    tB2 = tB + delta

    ab = RunningAB(count=np.zeros(4, dtype=np.int64), sum_ab=np.zeros(4, dtype=np.int64))
    tb = ThetaBinned(
        bins=theta_bins,
        count=np.zeros((theta_bins, 4), dtype=np.int64),
        sum_ab=np.zeros((theta_bins, 4), dtype=np.int64),
    )

    last_states = deque(maxlen=5)
    last_ab = deque(maxlen=5)
    cw_vals: List[float] = []
    ccw_vals: List[float] = []
    theta = 0.0

    i = 0
    j = 0
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

        # match i with j
        a = int(sA[i])
        b = int(sB[j])
        prod = int(oA[i] * oB[j])

        ab.update(a, b, prod)

        # theta-binned update (theta BEFORE loop closure update)
        th = theta % (2.0 * math.pi)
        bin_id = int((th / (2.0 * math.pi)) * theta_bins)
        if bin_id >= theta_bins:
            bin_id = theta_bins - 1
        tb.update(bin_id, a, b, prod)

        # loop detection in setting-space (sequence of (a,b) at coincident events)
        last_states.append((a, b))
        last_ab.append(prod)
        if len(last_states) == 5 and is_square_cycle(list(last_states)):
            area = signed_area(list(last_states))
            orient = 1 if area > 0 else (-1 if area < 0 else 0)

            cycle_val = float(np.mean(list(last_ab)[1:]))  # mean over 4 steps
            if orient > 0:
                ccw_vals.append(cycle_val)
            elif orient < 0:
                cw_vals.append(cycle_val)

            if orient != 0 and theta_kappa != 0.0:
                theta += theta_kappa * orient

        coincidences += 1
        i += 1
        j += 1

    E2 = ab.E2x2()
    E4 = ab.E4()

    # Weihs CHSH + max over relabelings
    S0 = S_weihs(E2)
    Smax, cfg = S_weihs_max_relabel(E2)

    # loop stats
    def mean_se(x: List[float]) -> Tuple[Optional[float], Optional[float], int]:
        if not x:
            return (None, None, 0)
        arr = np.asarray(x, dtype=float)
        m = float(arr.mean())
        if len(arr) < 2:
            return (jfloat(m), None, len(arr))
        se = float(arr.std(ddof=1) / math.sqrt(len(arr)))
        return (jfloat(m), jfloat(se), len(arr))

    cw_m, cw_se, cw_n = mean_se(cw_vals)
    ccw_m, ccw_se, ccw_n = mean_se(ccw_vals)

    if cw_m is not None and ccw_m is not None:
        diff = ccw_m - cw_m
    else:
        diff = None

    if cw_se is not None and ccw_se is not None:
        diff_se = math.sqrt(cw_se**2 + ccw_se**2)
        z = (diff / diff_se) if (diff_se > 0 and diff is not None) else None
    else:
        diff_se = None
        z = None

    S_bins, Nmin = tb.S_by_bin(min_pair_count=theta_min_pair_count)

    out: Dict[str, object] = {
        "coincidences": jint(coincidences),
        "E00_E01_E10_E11": [jfloat(E4[0]), jfloat(E4[1]), jfloat(E4[2]), jfloat(E4[3])],
        "CHSH_weihs_S": jfloat(S0),
        "CHSH_weihs_Smax_relabel": jfloat(Smax),
        "CHSH_weihs_best_relabeling": cfg,
        "x_theta_loop_test": {
            "cw":  {"mean": cw_m,  "se": cw_se,  "n_cycles": jint(cw_n)},
            "ccw": {"mean": ccw_m, "se": ccw_se, "n_cycles": jint(ccw_n)},
            "diff_ccw_minus_cw": jfloat(diff) if diff is not None else None,
            "diff_se": jfloat(diff_se) if diff_se is not None else None,
            "z_score_approx": jfloat(z) if z is not None else None,
        },
        "theta_conditioned": {
            "theta_kappa": jfloat(theta_kappa),
            "theta_bins": jint(theta_bins),
            "theta_min_pair_count": jint(theta_min_pair_count),
            "S_by_bin": [jfloat(x) for x in S_bins.tolist()],
            "min_counts_by_bin": [jint(x) for x in Nmin.tolist()],
            "S_bin_min": jfloat(np.nanmin(S_bins)) if np.any(~np.isnan(S_bins)) else None,
            "S_bin_max": jfloat(np.nanmax(S_bins)) if np.any(~np.isnan(S_bins)) else None,
        },
    }

    if keep_cycles:
        out["x_theta_loop_cycles"] = {
            "cw_vals": [jfloat(v) for v in cw_vals],
            "ccw_vals": [jfloat(v) for v in ccw_vals],
        }

    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--alice-zip", required=True)
    ap.add_argument("--bob-zip", required=True)
    ap.add_argument("--folder", default=None, help="Filter to folder name, e.g. longdist")
    ap.add_argument("--max-stems", type=int, default=10, help="Process up to this many stems (0 = all)")
    ap.add_argument("--window-ns", type=float, default=6.0, help="Coincidence half-window W in ns")
    ap.add_argument("--delta-ns", type=float, default=None, help="Fixed delta (ns) added to Bob times")
    ap.add_argument("--delta-auto", action="store_true", help="Estimate delta by scanning coincidence peak")
    ap.add_argument("--delta-scan-range-ns", type=float, default=100.0)
    ap.add_argument("--delta-scan-step-ns", type=float, default=1.0)
    ap.add_argument("--delta-scan-max-events", type=int, default=200000)

    ap.add_argument("--reverse-bits", action="store_true",
                    help="Swap the two lowest bits when decoding C.dat (try if results look wrong)")

    ap.add_argument("--theta-kappa", type=float, default=0.25, help="Theta increment per detected loop (radians)")
    ap.add_argument("--theta-bins", type=int, default=16)
    ap.add_argument("--theta-min-pair-count", type=int, default=20,
                    help="Bins with min(setting-pair counts) below this become null in S_by_bin")
    ap.add_argument("--keep-cycles", action="store_true", help="Store raw CW/CCW cycle values in JSON")

    ap.add_argument("--out", default="weihs_x_theta_results.json")
    ap.add_argument("--list-only", action="store_true", help="Only list matching stems and exit")
    args = ap.parse_args()

    W = args.window_ns * 1e-9

    with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
        ia = build_zip_index(za)
        ib = build_zip_index(zb)

        stems_a = set(list_stems(za, ia, args.folder))
        stems_b = set(list_stems(zb, ib, args.folder))
        stems = sorted(stems_a.intersection(stems_b))

        if not stems:
            raise SystemExit("No matching stems found (check --folder and zip contents).")

        print(f"Found {len(stems)} matching stems.")
        print("First 20 stems:")
        for s in stems[:20]:
            print("  ", s)

        if args.list_only:
            return

        if args.max_stems == 0:
            sel = stems
        else:
            sel = stems[:args.max_stems]

        results: Dict[str, object] = {
            "alice_zip": args.alice_zip,
            "bob_zip": args.bob_zip,
            "folder_filter": args.folder,
            "window_ns": args.window_ns,
            "delta_ns_input": args.delta_ns,
            "delta_auto": bool(args.delta_auto),
            "reverse_bits": bool(args.reverse_bits),
            "processed_stems": [],
        }

        for stem in sel:
            print("\n--- Processing stem:", stem)

            tA, sA, oA = load_events(za, ia, stem, reverse_bits=args.reverse_bits)
            tB, sB, oB = load_events(zb, ib, stem, reverse_bits=args.reverse_bits)

            if args.delta_ns is not None:
                delta = args.delta_ns * 1e-9
                delta_used_ns = args.delta_ns
            elif args.delta_auto:
                delta = estimate_delta_by_scan(
                    tA, tB, W,
                    scan_range_ns=args.delta_scan_range_ns,
                    scan_step_ns=args.delta_scan_step_ns,
                    max_events=args.delta_scan_max_events,
                    refine_rounds=2
                )
                delta_used_ns = delta * 1e9
            else:
                delta = 0.0
                delta_used_ns = 0.0

            r = process_stem(
                tA, sA, oA,
                tB, sB, oB,
                W=W,
                delta=delta,
                theta_kappa=args.theta_kappa,
                theta_bins=args.theta_bins,
                theta_min_pair_count=args.theta_min_pair_count,
                keep_cycles=args.keep_cycles
            )
            r["stem"] = stem
            r["delta_ns_used"] = jfloat(delta_used_ns)

            results["processed_stems"].append(r)

            print(f"coincidences={r['coincidences']}")
            print("E00,E01,E10,E11 =", r["E00_E01_E10_E11"])
            print("CHSH (Weihs) S =", r["CHSH_weihs_S"])
            print("CHSH Smax (relabel) =", r["CHSH_weihs_Smax_relabel"])
            lt = r["x_theta_loop_test"]
            print("loop diff (CCW-CW) =", lt["diff_ccw_minus_cw"], "z≈", lt["z_score_approx"])
            print("delta_ns_used =", r["delta_ns_used"])

        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2, allow_nan=False)

        print("\nWrote:", args.out)


if __name__ == "__main__":
    main()
