#!/usr/bin/env python3
"""
weihs_validator_v14.py

THE HOLDOUT VALIDATOR.
1. Splits valid stems into Set A (Discovery) and Set B (Validation).
2. Scans Kappa on Set A to find the "Natural Frequency".
3. Tests that specific Kappa on Set B.
4. If Set B modulates, the effect is PREDICTIVE (Real Physics).
"""

import argparse
import zipfile
import numpy as np
import math
import random
from dataclasses import dataclass, field
from collections import deque
from typing import Dict, List, Tuple

# ==========================================
# PART 1: CORE ENGINE
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
    bins: int = 8
    phase_sums: np.ndarray = field(default_factory=lambda: np.zeros((8, 4)))
    phase_counts: np.ndarray = field(default_factory=lambda: np.zeros((8, 4)))

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
    cmap = {}
    i, j, nA, nB = 0, 0, len(tA), len(tB)
    while i < nA and j < nB:
        d = tB_shift[j] - tA[i]
        if d < -W: j += 1
        elif d > W: i += 1
        else:
            cmap[i] = j
            i += 1; j += 1
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
# PART 2: LOGIC
# ==========================================

def get_modulation(data_tuples, kappa):
    engine = ContinuousTopology(theta_kappa=kappa)
    for (tA, sA, oA, tB, sB, oB, d) in data_tuples:
        process_stem(tA, sA, oA, tB, sB, oB, 1.3e-9, d, engine)
    
    s_vals = []
    for i in range(engine.bins):
        if engine.phase_counts[i].sum() < 20: continue
        E = np.zeros(4)
        for k in range(4):
            c = engine.phase_counts[i, k]
            E[k] = engine.phase_sums[i, k]/c if c>5 else 0
        s_vals.append(E[0]+E[1]+E[2]-E[3])
    
    if len(s_vals) < 4: return 0.0
    return max(s_vals) - min(s_vals)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alice-zip", required=True)
    parser.add_argument("--bob-zip", required=True)
    parser.add_argument("--reverse-bits", action="store_true")
    args = parser.parse_args()
    
    # 1. Load Data into RAM
    stems = [
        "longdist0:4.5", "longdist1:-10.0", "longdist23:-2.0",
        "longdist3:-14.0", "longdist31:29.0", "longdist32:-5.0",
        "longdist35:16.0", "longdist36:-9.0", "longdist37:-6.0",
        "longdist5:-22.0"
    ]
    
    print("Loading data...")
    loaded = []
    with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
        # Resolve names
        real_map = {}
        for n in za.namelist():
            if n.lower().endswith("_v.dat"):
                 cl = n.replace("\\","/").lstrip("/")
                 if cl.lower().startswith("alice/"): cl = cl[6:]
                 real_map[cl[:-6]] = n[:-6] # map short -> full

        for item in stems:
            short, d_str = item.split(":")
            # find match
            match_key = next((k for k in real_map if k.endswith(short)), None)
            if not match_key: continue
            
            full = real_map[match_key]
            # Strip Alice/Bob prefix from full name if present for loading
            if full.lower().startswith("alice/") or full.lower().startswith("bob/"):
                full = full.split("/", 1)[1]

            tA, sA, oA = load_raw_data(za, full, args.reverse_bits)
            tB, sB, oB = load_raw_data(zb, full, args.reverse_bits)
            loaded.append((tA, sA, oA, tB, sB, oB, float(d_str)*1e-9))

    # 2. Split Data
    random.seed(42) # Fixed seed for reproducibility
    random.shuffle(loaded)
    mid = len(loaded) // 2
    set_A = loaded[:mid]
    set_B = loaded[mid:]
    
    print(f"\nData Split: {len(set_A)} stems in Discovery (A), {len(set_B)} stems in Validation (B)")
    
    # 3. Discover on A
    print("\n--- PHASE 1: DISCOVERY (Scanning Set A) ---")
    kappas = [0.1, 0.2, 0.25, 0.5, 0.8, 1.0, 1.25, 1.5]
    best_k = 0
    best_mod = 0
    
    for k in kappas:
        mod = get_modulation(set_A, k)
        print(f"Kappa {k:<4.2f} -> Modulation {mod:.3f}")
        if mod > best_mod:
            best_mod = mod
            best_k = k
            
    print(f"Result: Best Kappa on Set A is {best_k} (Mod: {best_mod:.3f})")
    
    # 4. Validate on B
    print(f"\n--- PHASE 2: VALIDATION (Testing Set B at K={best_k}) ---")
    val_mod = get_modulation(set_B, best_k)
    print(f"Modulation on Set B: {val_mod:.3f}")
    
    print("-" * 40)
    if val_mod > 0.4: # Arbitrary "Strong Signal" threshold
        print("SUCCESS: The resonance predicts the holdout set!")
    else:
        print("WARNING: The signal did not generalize. It might be overfitting.")

if __name__ == "__main__":
    main()