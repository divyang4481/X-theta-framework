#!/usr/bin/env python3
"""
weihs_paper_plot_v18.py

THE RIGOROUS PLOTTER (Debug Edition)
================================================================================
DESCRIPTION:
This script implements the "X-Theta" framework analysis on experimental Bell Test 
data (specifically the Weihs et al. dataset). It treats the geometric phase 
(Theta) as a loop holonomy that accumulates based on the history of setting 
choices (00 -> 01 -> 11 -> 10 -> 00).

KEY ALGORITHMS:
1. HOLONOMY TRACKING:
   - The class `ContinuousTopology` tracks the history of Alice and Bob's settings.
   - When a closed loop (cycle) in setting space is detected, the geometric phase 
     `theta` is incremented by `theta_kappa` * orientation (+1/-1).

2. STRICT BINNING (The "Rigor" Filter):
   - Data is bucketed into `bins` based on the current value of `theta`.
   - To calculate the Bell Parameter S for a bin, we need data for ALL 4 setting 
     permutations (00, 01, 10, 11).
   - If ANY setting permutation in a bin has fewer than `MIN_COUNTS_PER_SETTING` 
     events, that entire bin is discarded. This prevents "divide by zero" errors 
     and prevents artificially high S values caused by missing subtraction terms.

3. ROBUST BOOTSTRAP:
   - Error bars are calculated by resampling the data 2000 times (with replacement)
     to generate a standard deviation (sigma), providing statistical confidence.

USAGE:
   python weihs_paper_plot_v18.py --alice-zip "A.zip" --bob-zip "B.zip" --reverse-bits
================================================================================
"""

import argparse
import zipfile
import numpy as np
import math
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from collections import deque
from typing import Dict, List, Tuple

# === CONFIGURATION ===
RESONANCE_KAPPA = 0.8       # The coupling constant for the geometric phase
BOOTSTRAP_ROUNDS = 2000     # Number of resamples for error bars
MIN_COUNTS_PER_SETTING = 20 # Minimum counts required per setting (++/+-/-+/--) per bin
BIN_COUNT = 12              # Number of angular bins for the plot

# ==========================================
# PART 1: CORE ENGINE & DATA LOADING
# ==========================================

def _norm_member(name: str) -> str:
    """Normalizes zip file paths to handle Windows/Linux slashes."""
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"): return "/".join(parts[1:])
    return n

def load_raw_data(z, stem, reverse_bits):
    """Extracts Time (t), Setting (s), and Outcome (o) from binary files."""
    nm = { _norm_member(n).lower(): n for n in z.namelist() }
    vk, ck = (stem+"_V.dat").lower(), (stem+"_C.dat").lower()
    
    # Robust file lookup
    if vk not in nm:
        vk = next((k for k in nm if k.endswith(stem.lower() + "_v.dat")), None)
    if ck not in nm:
        ck = next((k for k in nm if k.endswith(stem.lower() + "_c.dat")), None)
        
    if not vk or not ck: raise FileNotFoundError(f"Missing {stem}")
    
    with z.open(nm[vk], "r") as f: vb = f.read()
    with z.open(nm[ck], "r") as f: cb = f.read()
    
    t = np.frombuffer(vb, dtype=">f8") # Big-endian double (Time)
    c = np.frombuffer(cb, dtype=">u2") # Big-endian unsigned short (Choices)
    n = min(len(t), len(c))
    t = t[:n].astype(np.float64, copy=False)
    c = c[:n]
    
    # Bit parsing based on hardware setup
    if reverse_bits:
        det = ((c >> 1) & 1).astype(np.uint8)
        s   = (c & 1).astype(np.uint8)
    else:
        det = (c & 1).astype(np.uint8)
        s   = ((c >> 1) & 1).astype(np.uint8)
        
    o = (2 * det - 1).astype(np.int8) # Convert 0/1 to -1/+1
    return t, s, o

State = Tuple[int, int]
def get_orientation(states: List[State]) -> int:
    """Calculates the winding direction of a setting loop using Shoelace formula."""
    area = 0.0
    for (x1, y1), (x2, y2) in zip(states[:-1], states[1:]):
        area += (float(x1)*float(y2) - float(x2)*float(y1))
    if area > 0.1: return 1   # Clockwise
    if area < -0.1: return -1 # Counter-Clockwise
    return 0

@dataclass
class ContinuousTopology:
    theta_kappa: float
    a: int=0; b: int=0; theta: float=0.0
    history: deque = field(default_factory=lambda: deque(maxlen=5))
    bins: int = BIN_COUNT
    # Structure: [Bin Index][Setting (0-3)] -> List of outcomes (+1 or -1)
    raw_outcomes: List[List[List[int]]] = field(default_factory=list)

    def __post_init__(self):
        self.raw_outcomes = [[[] for _ in range(4)] for _ in range(self.bins)]

    def update(self, actor, val):
        """Updates the state history. If a loop closes, updates Theta."""
        if actor == 'A': self.a = int(val)
        else: self.b = int(val)
        
        self.history.append((self.a, self.b))
        
        # LOGIC: Update ONLY on closed loops (Holonomy)
        # We look for a pattern like A->B->C->D->A (Length 5, start==end)
        if len(self.history) == 5:
            h = list(self.history)
            if h[0] == h[-1] and len(set(h[:-1])) == 4:
                self.theta += self.theta_kappa * get_orientation(h)

    def capture(self, a, b, prod):
        """bins the coincidence product based on current Theta."""
        # Wrap theta to [0, 2pi)
        ph = self.theta % (2 * math.pi)
        
        # Calculate Bin Index
        idx = int((ph / (2 * math.pi)) * self.bins) % self.bins
        
        # Calculate Setting Index (00=0, 01=1, 10=2, 11=3)
        k = a * 2 + b
        
        self.raw_outcomes[idx][k].append(prod)

def process_stem(tA, sA, oA, tB, sB, oB, W, delta, engine):
    """
    Matches Alice and Bob timestamps.
    Returns: number of coincidences found.
    """
    tB_shift = tB + delta
    cmap = {}
    i, j, nA, nB = 0, 0, len(tA), len(tB)
    
    # Linear scan for coincidences window W
    while i < nA and j < nB:
        d = tB_shift[j] - tA[i]
        if d < -W: j += 1
        elif d > W: i += 1
        else:
            cmap[i] = j
            i += 1; j += 1
            
    # Chronological Replay
    evs = [(tA[i], 0, i) for i in range(nA)] + [(tB_shift[j], 1, j) for j in range(nB)]
    evs.sort(key=lambda x: x[0])
    
    count = 0
    for _, src, idx in evs:
        if src == 0: # Alice Event
            engine.update('A', sA[idx])
            if idx in cmap:
                bi = cmap[idx]
                # We have a coincidence!
                engine.capture(sA[idx], sB[bi], oA[idx] * oB[bi])
                count += 1
        else: # Bob Event
            engine.update('B', sB[idx])
            
    return count

# ==========================================
# PART 2: RIGOROUS STATISTICS
# ==========================================

def calculate_S_robust(bin_outcomes, rounds=2000):
    """
    Returns (S_mean, S_std) OR (None, None) if data is insufficient.
    """
    # 1. STRICT CHECK: Do we have enough data for every setting?
    counts = [len(bin_outcomes[k]) for k in range(4)]
    if any(c < MIN_COUNTS_PER_SETTING for c in counts):
        return None, None # Reject this bin

    # 2. Calculate Mean S
    # E(a,b) = Mean of products in that setting
    means = []
    for k in range(4):
        means.append(np.mean(bin_outcomes[k]))
    
    # CHSH Inequality: S = E(0,0) + E(0,1) + E(1,0) - E(1,1)
    # Mapping: 0->(0,0), 1->(0,1), 2->(1,0), 3->(1,1)
    S_mean = means[0] + means[1] + means[2] - means[3]
    
    # 3. Bootstrap Error
    s_boots = []
    for _ in range(rounds):
        b_means = []
        for k in range(4):
            dat = np.array(bin_outcomes[k])
            # Resample with replacement
            resamp = np.random.choice(dat, size=len(dat), replace=True)
            b_means.append(np.mean(resamp))
        s_boots.append(b_means[0] + b_means[1] + b_means[2] - b_means[3])
    
    return S_mean, np.std(s_boots)

# ==========================================
# PART 3: MAIN EXECUTION
# ==========================================

def main():
    parser = argparse.ArgumentParser(description="X-Theta Bell Test Plotter")
    parser.add_argument("--alice-zip", required=True, help="Path to Alice's data zip")
    parser.add_argument("--bob-zip", required=True, help="Path to Bob's data zip")
    parser.add_argument("--reverse-bits", action="store_true", help="Swap detector/setting bits")
    args = parser.parse_args()
    
    engine = ContinuousTopology(theta_kappa=RESONANCE_KAPPA)
    
    # Validated time-shift stems
    stems = [
        "longdist0:4.5", "longdist1:-10.0", "longdist23:-2.0",
        "longdist3:-14.0", "longdist31:29.0", "longdist32:-5.0",
        "longdist35:16.0", "longdist36:-9.0", "longdist37:-6.0",
        "longdist5:-22.0"
    ]
    
    print("==========================================")
    print("    WEIH Paper Plotter v18 (Debug)        ")
    print("==========================================")
    
    total_coincidences = 0
    
    with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
        full_paths = []
        # Index files
        for n in za.namelist():
             if n.lower().endswith("_v.dat"):
                 cl = n.replace("\\", "/").lstrip("/")
                 if cl.lower().startswith("alice/"): cl = cl[6:]
                 full_paths.append(cl[:-6])
        
        print(f"Found {len(full_paths)} data files in zip. Matching stems...")
        
        for item in stems:
            short, d_str = item.split(":")
            match = next((p for p in full_paths if p.endswith(short)), None)
            
            if not match: 
                print(f"Skipping {short} (Not found in zip)")
                continue
            
            # Load Data
            tA, sA, oA = load_raw_data(za, match, args.reverse_bits)
            tB, sB, oB = load_raw_data(zb, match, args.reverse_bits)
            
            # Process
            c_count = process_stem(tA, sA, oA, tB, sB, oB, 1.3e-9, float(d_str)*1e-9, engine)
            total_coincidences += c_count
            print(f" -> {short}: {c_count} coincidences processed.")

    print(f"\nTotal Coincidences Accumulated: {total_coincidences}")

    # --- DEBUG: PRINT BIN STATISTICS ---
    print("\n[DIAGNOSTICS] Bin Population Report:")
    print("Bin | Set00 | Set01 | Set10 | Set11 | Status")
    print("-------------------------------------------------")
    has_valid_data = False
    for i in range(engine.bins):
        counts = [len(engine.raw_outcomes[i][k]) for k in range(4)]
        valid = all(c >= MIN_COUNTS_PER_SETTING for c in counts)
        status = "READY" if valid else f"SKIP (<{MIN_COUNTS_PER_SETTING})"
        print(f"{i:3d} | {counts[0]:5d} | {counts[1]:5d} | {counts[2]:5d} | {counts[3]:5d} | {status}")
        if valid: has_valid_data = True

    if not has_valid_data:
        print("\nCRITICAL WARNING: No bins met the strict data requirements.")
        print("Try reducing MIN_COUNTS_PER_SETTING or checking if theta is evolving.")
        # We proceed to plot anyway, but it will be empty.

    # Compile Data for Plotting
    phases = []
    s_means = []
    s_errs = []
    
    print("\n[ANALYSIS] Calculating S parameters...")
    print(f"   Phase | S_mean | S_err (BS={BOOTSTRAP_ROUNDS})")
    print("   --------------------------------------")
    
    for i in range(engine.bins):
        S, err = calculate_S_robust(engine.raw_outcomes[i], BOOTSTRAP_ROUNDS)
        
        if S is None: continue
            
        ph = (i + 0.5) * (2*math.pi/engine.bins)
        phases.append(ph)
        s_means.append(S)
        s_errs.append(err)
        print(f"   {ph:4.2f}  | {S:5.2f}  | {err:5.2f}")

    # PLOT
    plt.figure(figsize=(10, 6))
    
    # Reference Lines
    plt.axhline(2.0, color='red', linestyle='--', linewidth=1.5, label='Classical Bound (S=2)')
    plt.axhline(2.828, color='green', linestyle=':', linewidth=1.5, label='Quantum Limit (Tsirelson)')
    plt.axhline(0.0, color='gray', linewidth=0.5)
    
    label_str = r'X-Theta Data ($\kappa=' + str(RESONANCE_KAPPA) + r'$)'
    
    if len(phases) > 0:
        plt.errorbar(phases, s_means, yerr=s_errs, fmt='o-', color='blue', 
                     ecolor='black', capsize=5, elinewidth=1, markeredgewidth=1,
                     label=label_str)
    else:
        plt.text(3, 1.5, "NO VALID DATA POINTS", color='red', fontsize=14, ha='center')
    
    plt.xlabel('Geometric Phase $\\theta$ (radians)')
    plt.ylabel('Bell Parameter $S$')
    plt.title('Phase-Dependent Bell Violation (Strict Binning & Bootstrap)')
    plt.legend(loc='lower right')
    plt.grid(True, alpha=0.3)
    plt.ylim(-0.5, 3.5) # Standard Bell Test range
    
    # Save
    out_name = "weihs_paper_plot_v18.png"
    plt.savefig(out_name)
    print(f"\nSaved plot to: {out_name}")

if __name__ == "__main__":
    main()