#!/usr/bin/env python3

"""
weihs_topology_lab_v3.py

A rigorous tool for the X-Theta Bell Test.
Separates "Time Filtering" (Finding the events) from "Topological Analysis" (The Theta Loop).

MODES:
  1) calibration: Scans time-shifts (delta) to find the "Physical Signal" (High CHSH S).
  2) analysis:    Runs the full Topological Test using a FIXED delta.

USAGE:
  python src/weihs_topology_lab_v3.py --alice-zip "..\Alice.zip" --bob-zip "..\Bob.zip" --folder longdist calibration --stem longdist20
  python src/weihs_topology_lab_v3.py --alice-zip "..\Alice.zip" --bob-zip "..\Bob.zip" --folder longdist analysis --delta-ns 1.5
"""

import argparse
import json
import math
import zipfile
import numpy as np
from dataclasses import dataclass, field
from collections import deque
from typing import Dict, List, Tuple, Optional

# ==========================================
# PART 1: HELPERS (Embedded to avoid import errors)
# ==========================================

def _norm_member(name: str) -> str:
    """Normalize zip member names."""
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
    """Return stems that have both _V.dat and _C.dat."""
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


# ==========================================
# PART 2: THE TOPOLOGICAL MACHINE
# ==========================================

State = Tuple[int, int]  # (a,b)

def is_square_cycle(states: List[State]) -> bool:
    if len(states) != 5: return False
    if states[0] != states[-1]: return False
    if len(set(states[:-1])) != 4: return False
    for u, v in zip(states[:-1], states[1:]):
        dist = (u[0] != v[0]) + (u[1] != v[1])
        if dist != 1: return False
    return True

def get_orientation(states: List[State]) -> int:
    area = 0.0
    for (x1, y1), (x2, y2) in zip(states[:-1], states[1:]):
        area += (x1 * y2 - x2 * y1)
    return 1 if area > 0 else -1

@dataclass
class TopologyEngine:
    theta_kappa: float
    theta: float = 0.0
    history: deque = field(default_factory=lambda: deque(maxlen=5))
    history_ab: deque = field(default_factory=lambda: deque(maxlen=5))
    cw_outcomes: List[float] = field(default_factory=list)
    ccw_outcomes: List[float] = field(default_factory=list)

    def process_event(self, a: int, b: int, prod: int):
        self.history.append((a, b))
        self.history_ab.append(prod)
        if is_square_cycle(list(self.history)):
            orient = get_orientation(list(self.history))
            loop_val = float(np.mean(list(self.history_ab)[1:]))
            if orient == 1:
                self.ccw_outcomes.append(loop_val)
            elif orient == -1:
                self.cw_outcomes.append(loop_val)
            self.theta += self.theta_kappa * orient

# ==========================================
# PART 3: ORCHESTRATION
# ==========================================

def run_coincidence_stream(tA, sA, oA, tB, sB, oB, W, delta, callback):
    tB_shifted = tB + delta
    i = j = 0
    nA, nB = len(tA), len(tB)
    while i < nA and j < nB:
        diff = tB_shifted[j] - tA[i]
        if diff < -W:
            j += 1
        elif diff > W:
            i += 1
        else:
            callback(sA[i], sB[j], oA[i] * oB[j])
            i += 1
            j += 1

def run_calibration(args, stems, za, zb):
    print(f"\n=== CALIBRATION MODE ===")
    target_stem = args.stem if args.stem else stems[0]
    print(f"Scanning Stem: {target_stem}")
    print(f"Search Range: {args.cal_min} ns to {args.cal_max} ns")
    
    tA, sA, oA = load_events(za, target_stem, args.reverse_bits)
    tB, sB, oB = load_events(zb, target_stem, args.reverse_bits)
    W = args.window_ns * 1e-9
    
    print("\nDelta(ns) | Coinc |   S_bell   | Status")
    print("-" * 45)
    
    best_S = 0.0
    best_delta = 0.0
    
    for d_ns in np.arange(args.cal_min, args.cal_max, args.cal_step):
        counts = np.zeros(4, dtype=int)
        sums = np.zeros(4, dtype=int)
        
        def count_accumulator(a, b, prod):
            idx = a * 2 + b
            counts[idx] += 1
            sums[idx] += prod
            
        run_coincidence_stream(tA, sA, oA, tB, sB, oB, W, d_ns * 1e-9, count_accumulator)
        
        E = np.zeros(4)
        nz = counts > 0
        E[nz] = sums[nz] / counts[nz]
        S = E[0] + E[1] + E[2] - E[3]
        total = counts.sum()
        
        status = ""
        if abs(S) > 2.0: status = "** VIOLATION **"
        elif abs(S) > 1.5: status = "*"
        
        print(f"{d_ns:8.2f} | {total:5d} | {S:10.4f} | {status}")
        
        if abs(S) > abs(best_S):
            best_S = S
            best_delta = d_ns

    print("-" * 45)
    print(f"RECOMMENDATION: Use --delta-ns {best_delta:.2f}")

def run_analysis(args, stems, za, zb):
    if args.delta_ns is None:
        print("ERROR: Analysis mode requires --delta-ns (use calibration first!)")
        return

    print(f"\n=== TOPOLOGICAL ANALYSIS MODE ===")
    print(f"Fixed Delta: {args.delta_ns} ns")
    print(f"Pooling {len(stems)} stems...")
    
    topo_engine = TopologyEngine(theta_kappa=args.theta_kappa)
    global_counts = np.zeros(4, dtype=int)
    global_sums = np.zeros(4, dtype=int)
    total_coinc = 0
    W = args.window_ns * 1e-9
    delta = args.delta_ns * 1e-9
    
    for stem in stems:
        tA, sA, oA = load_events(za, stem, args.reverse_bits)
        tB, sB, oB = load_events(zb, stem, args.reverse_bits)
        
        if args.null_shuffle:
            np.random.shuffle(oB)

        def pipeline(a, b, prod):
            idx = a * 2 + b
            global_counts[idx] += 1
            global_sums[idx] += prod
            topo_engine.process_event(a, b, prod)
            
        prev_coinc = total_coinc
        run_coincidence_stream(tA, sA, oA, tB, sB, oB, W, delta, pipeline)
        
        new_hits = global_counts.sum() - prev_coinc
        print(f"  > {stem}: +{new_hits} events")
        total_coinc = global_counts.sum()

    print("\n=== FINAL RESULTS ===")
    E = np.zeros(4)
    nz = global_counts > 0
    E[nz] = global_sums[nz] / global_counts[nz]
    S_global = E[0] + E[1] + E[2] - E[3]
    print(f"Global CHSH S: {S_global:.5f}")

    cw = topo_engine.cw_outcomes
    ccw = topo_engine.ccw_outcomes
    n_cw, n_ccw = len(cw), len(ccw)
    
    if n_cw > 1 and n_ccw > 1:
        m_cw = np.mean(cw)
        m_ccw = np.mean(ccw)
        se_pool = math.sqrt(np.var(cw, ddof=1)/n_cw + np.var(ccw, ddof=1)/n_ccw)
        z = (m_ccw - m_cw) / se_pool
        print(f"\nTopological Loop Test:")
        print(f"  CW Loops:  {n_cw} (mean {m_cw:.4f})")
        print(f"  CCW Loops: {n_ccw} (mean {m_ccw:.4f})")
        print(f"  Z-Score:   {z:.5f}")
    else:
        print("\nTopological Loop Test: Insufficient data")

def main():
    parser = argparse.ArgumentParser()
    # Global Args (MUST COME FIRST)
    parser.add_argument("--alice-zip", required=True)
    parser.add_argument("--bob-zip", required=True)
    parser.add_argument("--folder", default="longdist")
    parser.add_argument("--window-ns", type=float, default=6.0)
    parser.add_argument("--reverse-bits", action="store_true", default=True)

    subparsers = parser.add_subparsers(dest="mode", required=True)
    
    # Calibration Args
    p_cal = subparsers.add_parser("calibration")
    p_cal.add_argument("--stem", help="Specific stem")
    p_cal.add_argument("--cal-min", type=float, default=-5.0)
    p_cal.add_argument("--cal-max", type=float, default=20.0)
    p_cal.add_argument("--cal-step", type=float, default=0.5)
    
    # Analysis Args
    p_ana = subparsers.add_parser("analysis")
    p_ana.add_argument("--delta-ns", type=float)
    p_ana.add_argument("--theta-kappa", type=float, default=0.25)
    p_ana.add_argument("--null-shuffle", action="store_true")
    
    args = parser.parse_args()
    
    with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
        # Correctly find stems present in BOTH files
        stems_a = set(list_stems(za, args.folder))
        stems_b = set(list_stems(zb, args.folder))
        stems = sorted(list(stems_a.intersection(stems_b)))
        
        if not stems:
            print("No stems found.")
            return

        if args.mode == "calibration":
            run_calibration(args, stems, za, zb)
        elif args.mode == "analysis":
            run_analysis(args, stems, za, zb)

if __name__ == "__main__":
    main()