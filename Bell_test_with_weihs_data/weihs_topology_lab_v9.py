#!/usr/bin/env python3

"""
weihs_topology_lab_v9.py

The "Deep Harvest" Edition.
- Scans a wider time range (+/- 30ns) to catch drifting stems.
- Prints a full diagnostic table of ALL stems (Accepted & Rejected).
- Aggregates maximum possible events for the Topological Loop Test.

USAGE:
  python src/weihs_topology_lab_v9.py --alice-zip "..\Alice.zip" --bob-zip "..\Bob.zip" --folder longdist

  python .\weihs_topology_lab_v9.py   --alice-zip "..\Bell_data\7185335\Alice.zip"   --bob-zip   "..\Bell_data\7185335\Bob.zip"   --folder longdist   --reverse-bits   --scan-range 30.0
"""

import argparse
import zipfile
import numpy as np
import math
from dataclasses import dataclass, field
from collections import deque
from typing import Dict, List, Tuple, Optional

# ==========================================
# PART 1: HELPERS
# ==========================================

def _norm_member(name: str) -> str:
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"):
        n = "/".join(parts[1:])
    return n

def _zip_name_map(z: zipfile.ZipFile) -> Dict[str, str]:
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

def calculate_S(counts, sums):
    E = np.zeros(4)
    nz = counts > 0
    E[nz] = sums[nz] / counts[nz]
    return E[0] + E[1] + E[2] - E[3]

# ==========================================
# PART 2: TOPOLOGY ENGINE
# ==========================================

State = Tuple[int, int]

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
        x1f, y1f = float(x1), float(y1)
        x2f, y2f = float(x2), float(y2)
        area += (x1f * y2f - x2f * y1f)
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alice-zip", required=True)
    parser.add_argument("--bob-zip", required=True)
    parser.add_argument("--folder", default="longdist")
    
    # Aggressive Harvest Defaults
    parser.add_argument("--window-ns", type=float, default=1.3)
    parser.add_argument("--min-s", type=float, default=1.8)
    parser.add_argument("--scan-range", type=float, default=30.0, help="Wide scan for drifting clocks")
    
    parser.add_argument("--reverse-bits", action="store_true", help="Use bit-reversed decoding")
    parser.add_argument("--null-shuffle", action="store_true", help="Run as Null Test")
    
    args = parser.parse_args()

    print(f"\n=== X-THETA DEEP HARVESTER v9 ===")
    print(f"Window: {args.window_ns} ns | Min S: {args.min_s} | Scan Range: +/- {args.scan_range} ns")
    
    with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
        stems = list_stems(za, args.folder)
        print(f"Found {len(stems)} stems. Starting Deep Scan...\n")
        
        topo_engine = TopologyEngine(theta_kappa=0.25)
        
        total_stems_used = 0
        global_counts = np.zeros(4, dtype=int)
        global_sums = np.zeros(4, dtype=int)
        
        print(f"{'STEM':<20} | {'DELTA (ns)':<10} | {'EVENTS':<8} | {'MAX S':<8} | {'STATUS'}")
        print("-" * 70)

        for stem in stems:
            try:
                tA, sA, oA = load_events(za, stem, args.reverse_bits)
                tB, sB, oB = load_events(zb, stem, args.reverse_bits)
            except:
                continue

            W = args.window_ns * 1e-9
            best_S = 0.0
            best_delta = 0.0
            best_count = 0
            
            # Coarse Scan first to save time
            coarse_points = np.arange(-args.scan_range, args.scan_range, 1.0)
            for d in coarse_points:
                d_ns = d * 1e-9
                c_temp = np.zeros(4, dtype=int)
                s_temp = np.zeros(4, dtype=int)
                def quick_cb(a, b, prod):
                    idx = a * 2 + b
                    c_temp[idx]+=1
                    s_temp[idx]+=prod
                run_coincidence_stream(tA, sA, oA, tB, sB, oB, W, d_ns, quick_cb)
                
                count = c_temp.sum()
                if count > 10:
                    S_val = calculate_S(c_temp, s_temp)
                    if abs(S_val) > abs(best_S):
                        best_S = S_val
                        best_delta = d_ns
                        best_count = count

            stem_short = stem.split("/")[-1]
            
            if abs(best_S) >= args.min_s:
                status = "ACCEPT"
                total_stems_used += 1
                
                if args.null_shuffle:
                    np.random.shuffle(oB)
                
                def pipeline(a, b, prod):
                    idx = a * 2 + b
                    global_counts[idx] += 1
                    global_sums[idx] += prod
                    topo_engine.process_event(a, b, prod)
                
                # Re-run strict stream at best delta
                run_coincidence_stream(tA, sA, oA, tB, sB, oB, W, best_delta, pipeline)
            else:
                status = "REJECT"

            print(f"{stem_short:<20} | {best_delta*1e9:10.2f} | {best_count:<8} | {best_S:8.3f} | {status}")

        total_events = global_counts.sum()
        print("-" * 70)
        print(f"Stems Used:   {total_stems_used} / {len(stems)}")
        print(f"Total Events: {total_events}")
        
        if total_events == 0:
            return

        S_global = calculate_S(global_counts, global_sums)
        print(f"Global CHSH S: {S_global:.5f}")

        cw = topo_engine.cw_outcomes
        ccw = topo_engine.ccw_outcomes
        n_cw, n_ccw = len(cw), len(ccw)
        
        if n_cw + n_ccw > 0:
            m_cw = np.mean(cw) if n_cw else 0
            m_ccw = np.mean(ccw) if n_ccw else 0
            var_cw = np.var(cw, ddof=1) if n_cw > 1 else 0
            var_ccw = np.var(ccw, ddof=1) if n_ccw > 1 else 0
            se_pool = math.sqrt(var_cw/max(1,n_cw) + var_ccw/max(1,n_ccw))
            
            z = (m_ccw - m_cw) / se_pool if se_pool > 0 else 0
            
            print(f"\nTopological Loop Test:")
            print(f"  CW Loops:  {n_cw} (mean {m_cw:.4f})")
            print(f"  CCW Loops: {n_ccw} (mean {m_ccw:.4f})")
            print(f"  Z-Score:   {z:.5f}")
        else:
            print(f"\nTopological Loop Test: No loops found.")

if __name__ == "__main__":
    main()