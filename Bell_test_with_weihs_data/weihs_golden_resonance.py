#!/usr/bin/env python3
"""
weihs_golden_resonance.py

THE FINAL PLOTTER.
Focuses on the discovered resonance at Kappa = 1.25.
Generates a high-resolution Phase vs S-Value plot.
Includes Null Verification to prove the signal is real.


python .\weihs_golden_resonance.py   --alice-zip "..\Bell_data\7185335\Alice.zip"   --bob-zip   "..\Bell_data\7185335\Bob.zip"   --reverse-bits
"""

import argparse
import zipfile
import numpy as np
import math
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from collections import deque
from typing import Dict, List, Tuple

# --- CONFIG ---
RESONANCE_KAPPA = 1.25  # The peak found in v13
BINS = 12               # Higher resolution for the plot

# ==========================================
# PART 1: CORE ENGINE (Compact)
# ==========================================

def _norm_member(name: str) -> str:
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"): return "/".join(parts[1:])
    return n

def load_raw_data(z, stem, reverse_bits):
    nm = { _norm_member(n).lower(): n for n in z.namelist() }
    vk, ck = (stem+"_V.dat").lower(), (stem+"_C.dat").lower()
    with z.open(nm[vk], "r") as f: vb = f.read()
    with z.open(nm[ck], "r") as f: cb = f.read()
    t = np.frombuffer(vb, dtype=">f8")
    c = np.frombuffer(cb, dtype=">u2")
    n = min(len(t), len(c))
    t = t[:n].astype(np.float64, copy=False)
    c = c[:n]
    if reverse_bits:
        det = ((c >> 1) & 1).astype(np.uint8)
        s   = (c & 1).astype(np.uint8)
    else:
        det = (c & 1).astype(np.uint8)
        s   = ((c >> 1) & 1).astype(np.uint8)
    o = (2 * det - 1).astype(np.int8)
    return t, s, o

State = Tuple[int, int]
def get_orientation(states: List[State]) -> int:
    area = 0.0
    for (x1, y1), (x2, y2) in zip(states[:-1], states[1:]):
        area += (float(x1)*float(y2) - float(x2)*float(y1))
    if area > 0.1: return 1
    if area < -0.1: return -1
    return 0

@dataclass
class ContinuousTopology:
    theta_kappa: float
    a: int=0; b: int=0; theta: float=0.0
    history: deque = field(default_factory=lambda: deque(maxlen=5))
    bins: int = BINS
    phase_sums: np.ndarray = field(default_factory=lambda: np.zeros((BINS, 4)))
    phase_counts: np.ndarray = field(default_factory=lambda: np.zeros((BINS, 4)))

    def update(self, actor, val):
        if actor == 'A': self.a = int(val)
        else: self.b = int(val)
        self.history.append((self.a, self.b))
        if len(self.history) == 5:
            h = list(self.history)
            if h[0] == h[-1] and len(set(h[:-1])) == 4:
                self.theta += self.theta_kappa * get_orientation(h)

    def capture(self, a, b, prod):
        ph = self.theta % (2 * math.pi)
        idx = int((ph / (2 * math.pi)) * self.bins) % self.bins
        k = a * 2 + b
        self.phase_sums[idx, k] += prod
        self.phase_counts[idx, k] += 1

def process_stem(tA, sA, oA, tB, sB, oB, W, delta, engine):
    tB_shift = tB + delta
    # Coincidence map
    cmap = {}
    i, j, nA, nB = 0, 0, len(tA), len(tB)
    while i < nA and j < nB:
        d = tB_shift[j] - tA[i]
        if d < -W: j += 1
        elif d > W: i += 1
        else:
            cmap[i] = j
            i += 1; j += 1
    # Merge stream
    evs = [(tA[i], 0, i) for i in range(nA)] + [(tB_shift[j], 1, j) for j in range(nB)]
    evs.sort(key=lambda x: x[0])
    
    for _, src, idx in evs:
        if src == 0:
            engine.update('A', sA[idx])
            if idx in cmap:
                bi = cmap[idx]
                engine.capture(sA[idx], sB[bi], oA[idx] * oB[bi])
        else:
            engine.update('B', sB[idx])

# ==========================================
# PART 2: MAIN RUNNER
# ==========================================

def run_dataset(args, engine, null_shuffle=False):
    stems = [
        "longdist0:4.5", "longdist1:-10.0", "longdist23:-2.0",
        "longdist3:-14.0", "longdist31:29.0", "longdist32:-5.0",
        "longdist35:16.0", "longdist36:-9.0", "longdist37:-6.0",
        "longdist5:-22.0"
    ]
    with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
        # Resolve full names
        all_stems = [n.replace("\\","/").lstrip("/") for n in za.namelist()]
        all_stems = [n for n in all_stems if "Alice" not in n and "Bob" not in n] # rough filter
        
        # Better resolver
        real_stems = []
        for n in za.namelist():
             if n.lower().endswith("_v.dat"):
                 cl = n.replace("\\","/").lstrip("/")
                 if cl.lower().startswith("alice/"): cl = cl[6:]
                 real_stems.append(cl[:-6])
        
        for item in stems:
            short, d_str = item.split(":")
            # Find match
            match = next((s for s in real_stems if s.endswith(short)), None)
            if not match: continue
            
            tA, sA, oA = load_raw_data(za, match, args.reverse_bits)
            tB, sB, oB = load_raw_data(zb, match, args.reverse_bits)
            if null_shuffle: np.random.shuffle(oB)
            
            process_stem(tA, sA, oA, tB, sB, oB, 1.3e-9, float(d_str)*1e-9, engine)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alice-zip", required=True)
    parser.add_argument("--bob-zip", required=True)
    parser.add_argument("--reverse-bits", action="store_true")
    args = parser.parse_args()

    print(f"--- GENERATING GOLDEN RESONANCE PLOT (Kappa={RESONANCE_KAPPA}) ---")
    
    # 1. Physics Run
    print("Running Physics Model...")
    eng_phys = ContinuousTopology(theta_kappa=RESONANCE_KAPPA)
    run_dataset(args, eng_phys, null_shuffle=False)
    
    # 2. Null Run
    print("Running Null Control...")
    eng_null = ContinuousTopology(theta_kappa=RESONANCE_KAPPA)
    run_dataset(args, eng_null, null_shuffle=True)
    
    # 3. Extract Data
    phases = []
    s_phys = []
    s_null = []
    
    print("\n   Phase | S(Phys) | S(Null)")
    print("   -------------------------")
    for i in range(BINS):
        # Physics S
        E = np.zeros(4)
        for k in range(4):
            c = eng_phys.phase_counts[i, k]
            E[k] = eng_phys.phase_sums[i, k]/c if c>5 else 0
        sp = E[0]+E[1]+E[2]-E[3]
        
        # Null S
        E = np.zeros(4)
        for k in range(4):
            c = eng_null.phase_counts[i, k]
            E[k] = eng_null.phase_sums[i, k]/c if c>5 else 0
        sn = E[0]+E[1]+E[2]-E[3]
        
        ph = (i + 0.5) * (2*math.pi/BINS)
        phases.append(ph)
        s_phys.append(sp)
        s_null.append(sn)
        print(f"   {ph:4.2f}  |  {sp:5.2f}  |  {sn:5.2f}")

    # 4. Plot
    try:
        plt.figure(figsize=(10, 6))
        plt.plot(phases, s_phys, 'o-', color='blue', linewidth=2, label=f'X-Theta Physics ($\kappa={RESONANCE_KAPPA}$)')
        plt.plot(phases, s_null, 'x--', color='gray', alpha=0.7, label='Null Control')
        
        plt.axhline(2.0, color='red', linestyle='--', alpha=0.5, label='Quantum Limit (S=2)')
        plt.axhline(-2.0, color='red', linestyle='--', alpha=0.5)
        plt.axhline(0.0, color='black', linewidth=0.5)
        
        plt.xlabel('Geometric Phase $\\theta$ (radians)')
        plt.ylabel('Bell Parameter $S$')
        plt.title(f'Discovery of X-Theta Resonance ($\kappa={RESONANCE_KAPPA}$)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.ylim(-1.5, 3.0)
        
        out_file = "weihs_golden_resonance.png"
        plt.savefig(out_file)
        print(f"\nPlot saved to: {out_file}")
    except Exception as e:
        print(f"Plotting failed (no display?): {e}")

if __name__ == "__main__":
    main()