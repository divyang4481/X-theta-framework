#!/usr/bin/env python3

"""
weihs_topology_lab.py

A rigorous tool for the X–Theta Bell Test.
Separates "Time Filtering" (Finding the events) from "Topological Analysis" (The Theta Loop).

MODES:
  1) --mode calibration
     Scans a range of time-shifts (delta) on a SINGLE stem to find the "Physical Signal".
     Maximizes CHSH S, not just coincidence counts.
     
  2) --mode analysis
     Runs the full X-Theta Topological Test using a FIXED delta (found in calibration).
     Pools statistics across all stems to remove noise.
     
USAGE:
  1. CALIBRATE:
     python src/weihs_topology_lab.py --mode calibration --stem longdist20 ...
     (Look for Delta where |S| > 2.0)
     
  2. ANALYZE:
     python src/weihs_topology_lab.py --mode analysis --delta-ns 1.5 ...

EXAMPLE:

Step 1: Calibrate (Find the Clock) You need to find the specific delta where the Bell violation appears.

  python .\weihs_topology_lab_v2.py calibration `
  --alice-zip "..\Bell_data\7185335\Alice.zip" `
  --bob-zip "..\Bell_data\7185335\Bob.zip" `
  --folder longdist `
  --stem longdist20 `
  --cal-min -5 --cal-max 20

  python .\weihs_topology_lab_v2.py --alice-zip "..\Bell_data\7185335\Alice.zip" --bob-zip "..\Bell_data\7185335\Bob.zip" --folder longdist calibration --stem longdist20 --cal-min -5 --cal-max 20

Step 2: Analyze (Test the Topology) Now run the real test using that fixed value. The code will pool all stems together.

    python .\weihs_topology_lab_v2.py analysis `
    --alice-zip "..\Bell_data\7185335\Alice.zip" `
    --bob-zip "..\Bell_data\7185335\Bob.zip" `
    --folder longdist `
    --delta-ns 1.50 `
    --reverse-bits

Step 3: Verification (Null Test) Run the analysis again with --null-shuffle.

    python .\weihs_topology_lab_v2.py analysis ... --delta-ns 1.50 --null-shuffle

"""

import argparse
import json
import math
import zipfile
import sys
import numpy as np
from dataclasses import dataclass, field
from collections import deque
from typing import Dict, List, Tuple, Optional

# ==========================================
# PART 1: THE TOPOLOGICAL MACHINE (Time-Independent)
# ==========================================
# This section contains NO time variables. It operates purely on 
# sequences of states (a,b).

State = Tuple[int, int]  # (a,b)

def is_square_cycle(states: List[State]) -> bool:
    """
    Checks if the last 5 states form a closed square loop in configuration space.
    Purely geometric check.
    """
    if len(states) != 5: return False
    if states[0] != states[-1]: return False # Not closed
    
    # Must visit 4 unique corners
    if len(set(states[:-1])) != 4: return False
    
    # Check Hamming distance for each step (must be 1)
    for u, v in zip(states[:-1], states[1:]):
        dist = (u[0] != v[0]) + (u[1] != v[1])
        if dist != 1: return False
    return True

def get_orientation(states: List[State]) -> int:
    """Returns +1 (CCW) or -1 (CW) for a valid square cycle."""
    # Signed area calculation
    area = 0.0
    for (x1, y1), (x2, y2) in zip(states[:-1], states[1:]):
        area += (x1 * y2 - x2 * y1)
    return 1 if area > 0 else -1

@dataclass
class TopologyEngine:
    """
    The X-Theta State Machine.
    Ingests discrete events (a, b, outcome_prod).
    Maintains internal angle Theta.
    """
    theta_kappa: float
    
    # Internal state
    theta: float = 0.0
    history: deque = field(default_factory=lambda: deque(maxlen=5))
    history_ab: deque = field(default_factory=lambda: deque(maxlen=5))
    
    # Output statistics
    cw_outcomes: List[float] = field(default_factory=list)
    ccw_outcomes: List[float] = field(default_factory=list)

    def process_event(self, a: int, b: int, prod: int):
        """
        Step the machine forward by one vertex in configuration space.
        """
        self.history.append((a, b))
        self.history_ab.append(prod)
        
        # Check for loop closure
        if is_square_cycle(list(self.history)):
            orient = get_orientation(list(self.history))
            
            # Record the mean outcome of this loop
            # (We use the mean of the 4 steps involved in the loop)
            loop_val = float(np.mean(list(self.history_ab)[1:]))
            
            if orient == 1:
                self.ccw_outcomes.append(loop_val)
            elif orient == -1:
                self.cw_outcomes.append(loop_val)
                
            # HOLONOMY: Update internal Theta based on loop
            # This is the "Memory" effect
            self.theta += self.theta_kappa * orient


# ==========================================
# PART 2: THE TIME FILTER (Data Construction)
# ==========================================
# This section handles the dirty reality of experimental timing.
# Its ONLY job is to output a clean stream of (a,b) pairs.

def load_stem_data(za, zb, stem, reverse_bits):
    """Loads raw timestamps and codes from Zip files."""
    def get_arrays(z, s):
        name = [n for n in z.namelist() if s in n and n.endswith("_V.dat")][0]
        # Hacky reliable loader from before
        # ... (omitted for brevity, assuming standard structure)
        return np.array([]), np.array([]), np.array([]) 

    # Re-using your robust loader logic here
    # (Compressed for brevity - reusing logic from previous script)
    from weihs_bell_x_theta import load_events 
    tA, sA, oA = load_events(za, stem, reverse_bits)
    tB, sB, oB = load_events(zb, stem, reverse_bits)
    return tA, sA, oA, tB, sB, oB

def run_coincidence_stream(tA, sA, oA, tB, sB, oB, W, delta, callback):
    """
    Streaming Coincidence Logic.
    For every valid pair found, calls 'callback(a, b, prod)'.
    """
    tB_shifted = tB + delta
    i = j = 0
    nA, nB = len(tA), len(tB)
    
    # Standard linear merge
    while i < nA and j < nB:
        diff = tB_shifted[j] - tA[i]
        if diff < -W:
            j += 1
        elif diff > W:
            i += 1
        else:
            # HIT! Valid geometric vertex found.
            # Pass to the topological layer.
            callback(sA[i], sB[j], oA[i] * oB[j])
            i += 1
            j += 1

# ==========================================
# PART 3: ORCHESTRATION (Main Script)
# ==========================================

def run_calibration(args, stems, za, zb):
    """Finds the 'Golden Delta' where physics exists."""
    print(f"\n=== CALIBRATION MODE ===")
    target_stem = args.stem if args.stem else stems[0]
    print(f"Scanning Stem: {target_stem}")
    print(f"Search Range: {args.cal_min} ns to {args.cal_max} ns")
    
    tA, sA, oA, tB, sB, oB = load_stem_data(za, zb, target_stem, args.reverse_bits)
    W = args.window_ns * 1e-9
    
    print("\nDelta(ns) | Coinc |   S_bell   | Status")
    print("-" * 45)
    
    best_S = 0.0
    best_delta = 0.0
    
    # Coarse Scan
    for d_ns in np.arange(args.cal_min, args.cal_max, args.cal_step):
        # Accumulate stats locally
        counts = np.zeros(4, dtype=int)
        sums = np.zeros(4, dtype=int)
        
        def count_accumulator(a, b, prod):
            idx = a * 2 + b
            counts[idx] += 1
            sums[idx] += prod
            
        run_coincidence_stream(tA, sA, oA, tB, sB, oB, W, d_ns * 1e-9, count_accumulator)
        
        # Calc S
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
    """Runs the full Topological Test with Fixed Delta."""
    if args.delta_ns is None:
        print("ERROR: Analysis mode requires --delta-ns (use calibration first!)")
        return

    print(f"\n=== TOPOLOGICAL ANALYSIS MODE ===")
    print(f"Fixed Delta: {args.delta_ns} ns")
    print(f"Pooling {len(stems)} stems...")
    
    # Global Pool
    topo_engine = TopologyEngine(theta_kappa=args.theta_kappa)
    global_counts = np.zeros(4, dtype=int)
    global_sums = np.zeros(4, dtype=int)
    total_coinc = 0
    
    W = args.window_ns * 1e-9
    delta = args.delta_ns * 1e-9
    
    for stem in stems:
        # Load Data
        tA, sA, oA, tB, sB, oB = load_stem_data(za, zb, stem, args.reverse_bits)
        
        # Null Shuffle Check
        if args.null_shuffle:
            np.random.shuffle(oB) # Destroys physics, keeps timing

        # Define Callback to feed the Engine
        def pipeline(a, b, prod):
            # 1. Update CHSH Stats
            idx = a * 2 + b
            global_counts[idx] += 1
            global_sums[idx] += prod
            
            # 2. Update Topological Engine
            topo_engine.process_event(a, b, prod)
            
        # Run Stream
        prev_coinc = total_coinc
        run_coincidence_stream(tA, sA, oA, tB, sB, oB, W, delta, pipeline)
        
        # Progress
        new_hits = global_counts.sum() - prev_coinc
        print(f"  > {stem}: +{new_hits} events")
        total_coinc = global_counts.sum()

    # Final Report
    print("\n=== FINAL RESULTS ===")
    
    # 1. CHSH Baseline
    E = np.zeros(4)
    nz = global_counts > 0
    E[nz] = global_sums[nz] / global_counts[nz]
    S_global = E[0] + E[1] + E[2] - E[3]
    print(f"Global CHSH S: {S_global:.5f}")
    if abs(S_global) < 2.0 and not args.null_shuffle:
        print("WARNING: No Bell violation. Check calibration!")

    # 2. Topological Test
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
        print("\nTopological Loop Test: Insufficient data (not enough loops found)")

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode", required=True)
    
    # Common Args
    parser.add_argument("--alice-zip", required=True)
    parser.add_argument("--bob-zip", required=True)
    parser.add_argument("--folder", default="longdist")
    parser.add_argument("--window-ns", type=float, default=6.0)
    parser.add_argument("--reverse-bits", action="store_true", default=True) # Defaulting to True as per your findings
    
    # Calibration Args
    p_cal = subparsers.add_parser("calibration")
    p_cal.add_argument("--stem", help="Specific stem to scan (default: first found)")
    p_cal.add_argument("--cal-min", type=float, default=-5.0)
    p_cal.add_argument("--cal-max", type=float, default=20.0)
    p_cal.add_argument("--cal-step", type=float, default=0.5)
    
    # Analysis Args
    p_ana = subparsers.add_parser("analysis")
    p_ana.add_argument("--delta-ns", type=float, help="Fixed delta from calibration")
    p_ana.add_argument("--theta-kappa", type=float, default=0.25)
    p_ana.add_argument("--null-shuffle", action="store_true")
    
    args = parser.parse_args()
    
    # Load Zips
    from weihs_bell_x_theta import list_stems, load_events # Helper reuse
    # (Note: In standalone, you'd paste those helpers here. 
    #  For now, assume they exist or imports work.)
    
    with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
        stems = list_stems(za, args.folder)
        if not stems:
            print("No stems found.")
            return

        if args.mode == "calibration":
            run_calibration(args, stems, za, zb)
        elif args.mode == "analysis":
            run_analysis(args, stems, za, zb)

if __name__ == "__main__":
    main()