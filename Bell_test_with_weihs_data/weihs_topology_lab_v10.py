#!/usr/bin/env python3

"""
weihs_topology_lab_v10.py

THE "S-MODULATION" ENGINE.
Updates v10 to calculate the full CHSH S-parameter per phase bin.
This prevents the "Signal Dilution" where E00, E01, etc. cancel each other out.

USAGE:
  python src/weihs_topology_lab_v11.py --alice-zip "..\Alice.zip" --bob-zip "..\Bob.zip" --folder longdist --reverse-bits
"""

import argparse
import zipfile
import numpy as np
import math
from dataclasses import dataclass, field
from collections import deque
from typing import Dict, List, Tuple

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

def list_stems(z: zipfile.ZipFile, folder_filter: str) -> List[str]:
    names = [_norm_member(n) for n in z.namelist()]
    ff = folder_filter.replace("\\", "/").strip("/").lower()
    have_v, have_c = set(), set()
    for n in names:
        nl = n.lower()
        if (ff in nl) and nl.endswith("_v.dat"):
            have_v.add(n[:-6])
        elif (ff in nl) and nl.endswith("_c.dat"):
            have_c.add(n[:-6])
    return sorted(have_v.intersection(have_c))

def load_raw_data(z: zipfile.ZipFile, stem: str, reverse_bits: bool):
    nm = _zip_name_map(z)
    v_key = (stem + "_V.dat").lower()
    c_key = (stem + "_C.dat").lower()
    if v_key not in nm or c_key not in nm: raise FileNotFoundError(stem)
    
    with z.open(nm[v_key], "r") as f: vb = f.read()
    with z.open(nm[c_key], "r") as f: cb = f.read()
    
    t = np.frombuffer(vb, dtype=">f8")
    c = np.frombuffer(cb, dtype=">u2")
    n = min(len(t), len(c))
    t = t[:n].astype(np.float64, copy=False)
    c = c[:n]
    
    if reverse_bits:
        det = ((c >> 1) & 1).astype(np.uint8, copy=False)
        s   = (c & 1).astype(np.uint8, copy=False)
    else:
        det = (c & 1).astype(np.uint8, copy=False)
        s   = ((c >> 1) & 1).astype(np.uint8, copy=False)
        
    o = (2 * det - 1).astype(np.int8, copy=False)
    return t, s, o

# ==========================================
# PART 2: TOPOLOGY STATE MACHINE
# ==========================================

State = Tuple[int, int]

def get_orientation(states: List[State]) -> int:
    area = 0.0
    for (x1, y1), (x2, y2) in zip(states[:-1], states[1:]):
        x1f, y1f = float(x1), float(y1)
        x2f, y2f = float(x2), float(y2)
        area += (x1f * y2f - x2f * y1f)
    if area > 0.1: return 1   # CCW
    if area < -0.1: return -1 # CW
    return 0

@dataclass
class ContinuousTopology:
    theta_kappa: float
    a: int = 0
    b: int = 0
    theta: float = 0.0
    history: deque = field(default_factory=lambda: deque(maxlen=5))
    
    # Binned Results: [Bin][Setting_Index]
    # Setting Index = a*2 + b (00, 01, 10, 11)
    bins: int = 8
    phase_sums: np.ndarray = field(default_factory=lambda: np.zeros((8, 4)))
    phase_counts: np.ndarray = field(default_factory=lambda: np.zeros((8, 4)))

    def update_setting(self, actor: str, val: int):
        if actor == 'A': self.a = int(val)
        if actor == 'B': self.b = int(val)
        
        self.history.append((self.a, self.b))
        
        if len(self.history) == 5:
            h = list(self.history)
            if h[0] == h[-1] and len(set(h[:-1])) == 4:
                orient = get_orientation(h)
                self.theta += self.theta_kappa * orient

    def capture_coincidence(self, a: int, b: int, prod: int):
        phase = self.theta % (2 * math.pi)
        bin_idx = int((phase / (2 * math.pi)) * self.bins) % self.bins
        
        idx = a * 2 + b
        self.phase_sums[bin_idx, idx] += prod
        self.phase_counts[bin_idx, idx] += 1

# ==========================================
# PART 3: MAIN
# ==========================================

def process_continuous_stem(tA, sA, oA, tB, sB, oB, W, delta, engine: ContinuousTopology):
    tB_shifted = tB + delta
    
    coinc_map = {}
    i = j = 0
    nA, nB = len(tA), len(tB)
    while i < nA and j < nB:
        diff = tB_shifted[j] - tA[i]
        if diff < -W: j += 1
        elif diff > W: i += 1
        else:
            coinc_map[i] = j
            i += 1
            j += 1
            
    events = []
    for i in range(nA): events.append((tA[i], 0, i))
    for j in range(nB): events.append((tB_shifted[j], 1, j))
    events.sort(key=lambda x: x[0])
    
    for _, src, idx in events:
        if src == 0: # Alice
            engine.update_setting('A', sA[idx])
            if idx in coinc_map:
                b_idx = coinc_map[idx]
                prod = oA[idx] * oB[b_idx]
                # Pass current settings to capture
                engine.capture_coincidence(sA[idx], sB[b_idx], prod)
        else: # Bob
            engine.update_setting('B', sB[idx])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alice-zip", required=True)
    parser.add_argument("--bob-zip", required=True)
    parser.add_argument("--folder", default="longdist")
    parser.add_argument("--reverse-bits", action="store_true")
    parser.add_argument("--window-ns", type=float, default=1.3)
    
    # GOLDEN STEMS from Harvest v9
    default_stems = [
        "longdist0:4.5", "longdist1:-10.0", "longdist23:-2.0",
        "longdist3:-14.0", "longdist31:29.0", "longdist32:-5.0",
        "longdist35:16.0", "longdist36:-9.0", "longdist37:-6.0",
        "longdist5:-22.0"
    ]
    
    args = parser.parse_args()
    
    print(f"\n=== X-THETA LAB v11 (S-Modulation) ===")
    print("Bins: 8 | Kappa: 0.5 | Metric: CHSH S(theta)")
    
    engine = ContinuousTopology(theta_kappa=0.5)
    
    with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
        all_stems = list_stems(za, args.folder)
        for item in default_stems:
            name, d_str = item.split(":")
            # Match
            matches = [s for s in all_stems if s.endswith(name)]
            if not matches: continue
            full_name = matches[0]
            
            print(f"Processing {name:<12} (Delta {d_str} ns)...")
            tA, sA, oA = load_raw_data(za, full_name, args.reverse_bits)
            tB, sB, oB = load_raw_data(zb, full_name, args.reverse_bits)
            
            process_continuous_stem(
                tA, sA, oA, tB, sB, oB, 
                args.window_ns * 1e-9, 
                float(d_str) * 1e-9, 
                engine
            )

    print("\n=== S-VALUE MODULATION RESULTS ===")
    print(f"{'Phase':<8} | {'N_events':<8} | {'S-Value':<10} | {'Status'}")
    print("-" * 50)
    
    s_values = []
    
    for i in range(engine.bins):
        # Calculate S for this bin
        # S = E00 + E01 + E10 - E11
        E = np.zeros(4)
        for k in range(4):
            cnt = engine.phase_counts[i, k]
            sm  = engine.phase_sums[i, k]
            E[k] = sm / cnt if cnt > 10 else 0.0 # Min stats protection
            
        S = E[0] + E[1] + E[2] - E[3]
        total_ev = engine.phase_counts[i].sum()
        s_values.append(S)
        
        center = (i + 0.5) * (2 * math.pi / engine.bins)
        status = "** VIOLATION" if abs(S) > 2.0 else "Classical"
        print(f"{center:4.2f} rad | {int(total_ev):<8} | {S:6.3f}     | {status}")

    print("-" * 50)
    
    # Check Modulation
    valid_s = [s for s, c in zip(s_values, engine.phase_counts.sum(axis=1)) if c > 50]
    if valid_s:
        mod = max(valid_s) - min(valid_s)
        print(f"S-Modulation Depth: {mod:.3f}")
        if mod > 0.4:
            print(">> SIGNIFICANT MODULATION: Phase affects Quantum Strength!")
    else:
        print("Not enough data per bin to determine modulation.")

if __name__ == "__main__":
    main()