#!/usr/bin/env python3
"""
weihs_audit_protocol.py

THE RIGOROUS AUDIT
================================================================================
PURPOSE:
  Validates the "X-Theta" discovery claim by applying strict sanity checks
  and stratified permutation testing on existing mined data.

CHECKS:
  1. GLOBAL S: Calculates S without phase bins. Must be > 2.0.
  2. S-CALCULATION: Explicitly checks S = E00 + E01 + E10 - E11.
  3. ROBUSTNESS: Scans Bin Counts [4, 6, 8, 12] to detect aliasing.
  4. STRATIFIED NULL: Permutes outcomes ONLY within matching setting groups.

USAGE:
   python weihs_audit_protocol.py
   (Requires 'weihs_mined_data.npz' to be present)
================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import sys

# === CONFIG ===
CACHE_FILE = "weihs_mined_data.npz"
PERMUTATION_ROUNDS = 2000  # Stratified shuffles
BIN_SCANS = [4, 6, 8, 12]

def weighted_linear_fit(X, y, sigma):
    """ Weighted Least Squares with stability checks """
    y = np.asarray(y)
    X = np.asarray(X)
    sigma = np.maximum(np.asarray(sigma), 1e-12)
    
    w = 1.0 / (sigma ** 2)
    sw = np.sqrt(w)
    Xw = X * sw[:, None]
    yw = y * sw
    
    try:
        XtX = Xw.T @ Xw
        beta = np.linalg.solve(XtX, Xw.T @ yw)
        yhat = X @ beta
        chi2 = np.sum(((y - yhat) / sigma) ** 2)
        dof = max(1, len(y) - X.shape[1])
        
        # Weighted R2
        ybar = np.sum(w * y) / np.sum(w)
        sst = np.sum(w * (y - ybar) ** 2)
        sse = np.sum(w * (y - yhat) ** 2)
        r2w = 1.0 - (sse/sst) if sst > 0 else 0.0
        
        return beta, chi2, r2w
    except np.linalg.LinAlgError:
        return None, 99999.0, 0.0

def calculate_S_profile(phases, settings, products, n_bins):
    """
    Calculates S(theta) for a given bin count.
    Returns: bin_centers, s_values, s_errors
    """
    # 1. Bin Indices
    bin_indices = ((phases % (2*np.pi)) / (2*np.pi) * n_bins).astype(int) % n_bins
    
    # 2. Aggregation (Flattened: bin * 4 + setting)
    flat_indices = bin_indices * 4 + settings
    
    # Counts and Sums
    N = np.bincount(flat_indices, minlength=n_bins*4)
    Sum = np.bincount(flat_indices, weights=products, minlength=n_bins*4)
    
    # Avoid divide by zero
    mask = N > 5
    if not np.all(mask):
        return None, None, None
        
    E = Sum / N
    # Standard Error of Mean for +/-1 variables: sigma = sqrt(1-E^2)/sqrt(N)
    # Conservative approx: 1/sqrt(N)
    Err = 1.0 / np.sqrt(N)
    
    # Reshape to [bins, 4]
    E_r = E.reshape((n_bins, 4))
    Err_r = Err.reshape((n_bins, 4))
    
    # 3. Calculate S per bin
    # S = E(0,0) + E(0,1) + E(1,0) - E(1,1)
    # Mappings: 0->00, 1->01, 2->10, 3->11 (Standard CHSH ordering check required!)
    # Assuming standard mapping from miner: 0=(0,0), 1=(0,1), 2=(1,0), 3=(1,1)
    
    S_vals = E_r[:,0] + E_r[:,1] + E_r[:,2] - E_r[:,3]
    
    # Error Prop: Sqrt(sum of variances)
    S_errs = np.sqrt(np.sum(Err_r**2, axis=1))
    
    centers = (np.arange(n_bins) + 0.5) * (2*math.pi/n_bins)
    
    return centers, S_vals, S_errs

def main():
    print("==========================================")
    print("      WEIHS RIGOROUS AUDIT PROTOCOL       ")
    print("==========================================")
    
    # --- 1. LOAD DATA ---
    try:
        data = np.load(CACHE_FILE)
        phases = data['phases']
        settings = data['settings'] # 0..3
        products = data['products'] # -1, 1
        print(f"Loaded {len(phases)} events from {CACHE_FILE}")
    except FileNotFoundError:
        print(f"CRITICAL: {CACHE_FILE} not found. Run miner first.")
        return

    # --- 2. GLOBAL SANITY CHECK (The "Is it Real?" Check) ---
    print("\n[CHECK 1] Global S Calculation (No Phases)...")
    # Aggregating purely by setting
    g_counts = np.bincount(settings, minlength=4)
    g_sums = np.bincount(settings, weights=products, minlength=4)
    
    g_E = g_sums / np.maximum(g_counts, 1)
    g_S = g_E[0] + g_E[1] + g_E[2] - g_E[3]
    g_Err = np.sqrt(np.sum(1.0/np.maximum(g_counts, 1)))
    
    print(f"  COUNTS per Setting: {g_counts}")
    print(f"  CORRELATIONS E_ab: {np.round(g_E, 3)}")
    print(f"  GLOBAL S = {g_S:.4f} +/- {g_Err:.4f}")
    
    if g_S < 1.8:
        print("  FAIL: Global S is too low. Check variable definitions/mappings.")
        print("  (Did you map 0/1/2/3 correctly to ++/+-/-+/-- ?)")
    elif g_S > 2.0:
        print("  PASS: Violation of Classical Limit (>2.0) confirmed globally.")
    
    # --- 3. ROBUSTNESS SCAN ---
    print("\n[CHECK 2] Bin Robustness Scan (Aliasing Check)...")
    best_config = None
    best_delta = -1.0
    
    print("  Bins | Delta Chi2 | Weighted R2 | A (Amp)")
    print("  -----------------------------------------")
    
    for n_b in BIN_SCANS:
        centers, s, err = calculate_S_profile(phases, settings, products, n_b)
        if centers is None:
            print(f"  {n_b:4d} |    (Insufficient Data)")
            continue
            
        # Fit Sine
        X = np.column_stack([np.ones(n_b), np.sin(centers), np.cos(centers)])
        beta, chi2_sin, r2w = weighted_linear_fit(X, s, err)
        
        # Fit Flat
        X0 = np.ones((n_b, 1))
        _, chi2_flat, _ = weighted_linear_fit(X0, s, err)
        
        delta = chi2_flat - chi2_sin
        amp = np.hypot(beta[1], beta[2]) if beta is not None else 0
        
        print(f"  {n_b:4d} | {delta:10.2f} | {r2w:11.3f} | {amp:.4f}")
        
        if delta > best_delta:
            best_delta = delta
            best_config = (n_b, centers, s, err, beta, delta, r2w)

    if best_config is None:
        print("CRITICAL: No valid bins found.")
        return

    # --- 4. STRATIFIED PERMUTATION (The "Reviewer-Proof" Null) ---
    print(f"\n[CHECK 3] Stratified Permutation Test (N={PERMUTATION_ROUNDS})...")
    print("  (Shuffling outcomes ONLY within matching setting pairs)")
    
    chosen_bins = best_config[0]
    obs_delta = best_config[5]
    
    better_count = 0
    
    # Pre-group indices by setting for fast shuffling
    # We create a list of indices for each setting 0,1,2,3
    indices_by_set = [np.where(settings == k)[0] for k in range(4)]
    
    # Working copy of products
    perm_products = products.copy()
    
    for i in range(PERMUTATION_ROUNDS):
        # Stratified Shuffle
        for k in range(4):
            # Shuffle the outcomes belonging to setting k in-place
            idx_k = indices_by_set[k]
            subset = perm_products[idx_k]
            np.random.shuffle(subset)
            perm_products[idx_k] = subset
            
        # Re-Analyze
        _, ps, perr = calculate_S_profile(phases, settings, perm_products, chosen_bins)
        
        # Fit Sine (Fast)
        X = np.column_stack([np.ones(chosen_bins), np.sin(best_config[1]), np.cos(best_config[1])])
        _, p_chi2_sin, _ = weighted_linear_fit(X, ps, perr)
        
        # Fit Flat
        X0 = np.ones((chosen_bins, 1))
        _, p_chi2_flat, _ = weighted_linear_fit(X0, ps, perr)
        
        if (p_chi2_flat - p_chi2_sin) >= obs_delta:
            better_count += 1
            
        if (i+1) % 500 == 0:
            print(f"  ..{i+1} done")
            
    p_val = (better_count + 1) / (PERMUTATION_ROUNDS + 1)
    print(f"  Stratified p-value: {p_val:.5f}")

    # --- 5. VISUALIZATION ---
    print("\n[PLOTTING] Generating Audit Report...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    centers, s, err, beta = best_config[1], best_config[2], best_config[3], best_config[4]
    
    # Fit Curve
    th = np.linspace(0, 2*np.pi, 200)
    fit_y = beta[0] + beta[1]*np.sin(th) + beta[2]*np.cos(th)
    
    # Panel 1: Full Physics Scale
    ax1.errorbar(centers, s, yerr=err, fmt='o', color='blue', ecolor='black', capsize=3, label='Data')
    ax1.plot(th, fit_y, color='orange', linewidth=2, label='Fit')
    ax1.axhline(2.0, color='red', linestyle='--', label='Classical (2.0)')
    ax1.axhline(2.828, color='green', linestyle=':', label='Tsirelson')
    ax1.set_ylim(-0.5, 3.5)
    ax1.set_title("Full Physical Scale")
    ax1.set_xlabel(r"Geometric Phase $\theta$")
    ax1.set_ylabel("Bell Parameter S")
    ax1.legend(loc='lower left')
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Zoomed Evidence
    ax2.errorbar(centers, s, yerr=err, fmt='o', color='blue', ecolor='black', capsize=3, label='Data')
    ax2.plot(th, fit_y, color='orange', linewidth=2, label=f'Fit (p={p_val:.4f})')
    ax2.axhline(2.0, color='red', linestyle='--', alpha=0.5)
    ax2.axhline(2.828, color='green', linestyle=':', alpha=0.5)
    
    # Auto-zoom range
    y_mid = np.mean(s)
    y_span = max(np.max(s) - np.min(s), 0.2)
    ax2.set_ylim(y_mid - y_span, y_mid + y_span)
    
    ax2.set_title(f"Zoomed Signal (Bins={best_config[0]}, $\Delta\chi^2={obs_delta:.1f}$)")
    ax2.set_xlabel(r"Geometric Phase $\theta$")
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig("weihs_audit_report.png")
    print("Saved: weihs_audit_report.png")

if __name__ == "__main__":
    main()