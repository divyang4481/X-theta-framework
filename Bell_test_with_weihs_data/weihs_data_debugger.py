#!/usr/bin/env python3
"""
weihs_data_debugger_v2.py

THE "X-RAY" TOOL (Crash-Proof Edition)
================================================================================
CHANGES:
  - Added strict checks to ensure the target file exists in BOTH Alice and Bob 
    zip archives to prevent 'NoneType' crashes.
  - If the default target (longdist33) is missing/broken, it automatically 
    falls back to the first valid 'longdist' file it finds.

USAGE:
   python weihs_data_debugger_v2.py --alice-zip "A.zip" --bob-zip "B.zip"
================================================================================
"""

import argparse
import zipfile
import numpy as np
import sys

# === CONFIG ===
PREFERRED_TARGET = "longdist33" 
COINCIDENCE_WINDOW = 3.5e-9

def _norm_member(name: str) -> str:
    n = name.replace("\\", "/").lstrip("/")
    parts = n.split("/")
    if parts and parts[0].lower() in ("alice", "bob"): return "/".join(parts[1:])
    return n

def get_raw_buffers(z, stem):
    nm = { _norm_member(n).lower(): n for n in z.namelist() }
    vk, ck = (stem+"_V.dat").lower(), (stem+"_C.dat").lower()
    
    # Try exact match first
    if vk not in nm: 
        # Fallback to suffix match
        vk = next((k for k in nm if k.endswith(stem.lower() + "_v.dat")), None)
    if ck not in nm: 
        ck = next((k for k in nm if k.endswith(stem.lower() + "_c.dat")), None)
        
    if not vk or not ck: return None, None
    
    with z.open(nm[vk], "r") as f: vb = f.read()
    with z.open(nm[ck], "r") as f: cb = f.read()
    return np.frombuffer(vb, dtype=">f8"), np.frombuffer(cb, dtype=">u2")

def match_nearest_neighbor(tA, tB, window):
    if len(tB) == 0: return np.array([]), np.array([])
    idx = np.searchsorted(tB, tA)
    idx = np.clip(idx, 0, len(tB)-1)
    dt_right = tB[idx] - tA
    idx_left = np.maximum(idx - 1, 0)
    dt_left = tB[idx_left] - tA
    use_left = np.abs(dt_left) < np.abs(dt_right)
    best_idx = np.where(use_left, idx_left, idx)
    best_dt = tB[best_idx] - tA
    mask = np.abs(best_dt) < window
    return np.where(mask)[0], best_idx[mask]

def find_shift(tA, tB):
    # Quick Sync
    limit = min(len(tA), len(tB), 20000)
    tA_s = tA[:limit] * 1e9; tB_s = tB[:limit] * 1e9
    diffs = tB_s[:, None] - tA_s[None, :]
    valid = diffs[np.abs(diffs) < 100000.0]
    if len(valid) == 0: return 0.0
    
    hist, edges = np.histogram(valid, bins=5000)
    peak_bin = np.argmax(hist)
    peak_val = (edges[peak_bin] + edges[peak_bin+1])/2.0
    
    # Refine
    valid_fine = valid[np.abs(valid - peak_val) < 100]
    if len(valid_fine) == 0: return 0.0
    
    hist_f, edges_f = np.histogram(valid_fine, bins=200)
    peak_f = np.argmax(hist_f)
    return (edges_f[peak_f] + edges_f[peak_f+1])/2.0 * 1e-9

def analyze_bits(vals):
    """Returns which bits are actually toggling in the dataset."""
    active = []
    for i in range(16):
        bit_col = (vals >> i) & 1
        if np.std(bit_col) > 0.01: # It varies
            active.append(i)
    return active

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alice-zip", required=True); parser.add_argument("--bob-zip", required=True)
    args = parser.parse_args()

    print("==========================================")
    print("   WEIHS DATA X-RAY v2 (DEBUGGER)")
    print("==========================================")
    
    with zipfile.ZipFile(args.alice_zip, "r") as za, zipfile.ZipFile(args.bob_zip, "r") as zb:
        # 1. Compute Safe Intersection of Files
        a_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in za.namelist() if "_v.dat" in n.lower() }
        b_files = { _norm_member(n).lower().replace("_v.dat","").replace("_c.dat","") for n in zb.namelist() if "_v.dat" in n.lower() }
        all_stems = sorted(list(a_files.intersection(b_files)))
        
        # 2. Select Target
        target = next((s for s in all_stems if PREFERRED_TARGET in s), None)
        if not target:
            print(f"Warning: '{PREFERRED_TARGET}' not found in both archives.")
            # Fallback to any longdist file
            target = next((s for s in all_stems if "longdist" in s), None)
            
        if not target:
            print("CRITICAL: No 'longdist' files found in both archives.")
            return

        print(f"Inspecting file: {target}")
        
        # 3. Load Safely
        tA, cA = get_raw_buffers(za, target)
        tB, cB = get_raw_buffers(zb, target)
        
        if tA is None or tB is None:
            print("Error: Failed to load binary buffers.")
            return

        # 4. Sync
        shift = find_shift(tA, tB)
        print(f"Sync Shift: {shift*1e9:.1f} ns")
        
        # 5. Match
        tB_shifted = tB - shift
        idxA, idxB = match_nearest_neighbor(tA, tB_shifted, COINCIDENCE_WINDOW)
        print(f"Matched Pairs: {len(idxA)}")
        
        if len(idxA) == 0: 
            print("No matches found. Check sync window.")
            return

        # === RAW DUMP ===
        print("\n[RAW DATA DUMP - First 20 Matches]")
        print("Idx | Alice Raw (Int) | Alice Binary     | Bob Raw (Int) | Bob Binary       ")
        print("----|-----------------|------------------|---------------|------------------")
        
        for k in range(min(20, len(idxA))):
            a_val = cA[idxA[k]]
            b_val = cB[idxB[k]]
            a_bin = f"{a_val:016b}"
            b_bin = f"{b_val:016b}"
            print(f"{k:3d} | {a_val:15d} | {a_bin} | {b_val:13d} | {b_bin}")

        # === BIT STATISTICS ===
        print("\n[BIT STATISTICS]")
        a_active = analyze_bits(cA[idxA])
        b_active = analyze_bits(cB[idxB])
        print(f"Alice Active Bits: {a_active}")
        print(f"Bob   Active Bits: {b_active}")
        
        # Check One-Hot Hypothesis for lowest bits
        a_low = cA[idxA] & 3 # Bits 0 and 1
        b_low = cB[idxB] & 3
        print(f"\n[Low-Bit Histogram (0-3)]")
        valA, countA = np.unique(a_low, return_counts=True)
        valB, countB = np.unique(b_low, return_counts=True)
        
        print(f"Alice Values: {dict(zip(valA, countA))}")
        print(f"Bob   Values: {dict(zip(valB, countB))}")
        
        print("\nINTERPRETATION GUIDE:")
        print(" - If values are {0, 1}, Bit 0 is the data.")
        print(" - If values are {1, 2}, it is ONE-HOT (Bit 0 = '+', Bit 1 = '-').")
        print(" - If active bits are [4, 5], your mining scan (0-3) missed them.")

if __name__ == "__main__":
    main()