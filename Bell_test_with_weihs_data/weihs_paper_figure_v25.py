#!/usr/bin/env python3
"""
weihs_paper_figure_v25.py

THE FINAL PLOT GENERATOR
================================================================================
Reads 'weihs_mined_v25.npz' and generates the publication-quality figure.
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import os

CACHE_FILE = "weihs_mined_v25.npz"
BIN_COUNT = 8

def weighted_fit_sine(th, y, sig):
    w = 1.0/sig**2; sw = np.sqrt(w)
    X = np.column_stack([np.ones_like(th), np.sin(th), np.cos(th)])
    beta = np.linalg.lstsq(X*sw[:,None], y*sw, rcond=None)[0]
    yhat = X @ beta
    chi2 = np.sum(((y-yhat)/sig)**2)
    X0 = np.ones((len(th), 1))
    beta0 = np.linalg.lstsq(X0*sw[:,None], y*sw, rcond=None)[0]
    chi2_0 = np.sum(((y-X0 @ beta0)/sig)**2)
    return beta, chi2_0 - chi2

def bin_data(ph, ks, pr, bins):
    b = ((ph % (2*np.pi))/(2*np.pi)*bins).astype(int) % bins
    flat = b*4 + ks
    c = np.bincount(flat, minlength=bins*4).reshape(bins,4)
    s = np.bincount(flat, weights=pr, minlength=bins*4).reshape(bins,4)
    valid = np.all(c > 10, axis=1)
    if np.sum(valid) < 4: return None
    m = s[valid]/c[valid]
    S = m[:,0]+m[:,1]+m[:,2]-m[:,3]
    err = np.sqrt(np.sum((1-m**2)/c[valid], axis=1))
    return (np.where(valid)[0]+0.5)*(2*np.pi/bins), S, err

def main():
    if not os.path.exists(CACHE_FILE):
        print(f"Error: {CACHE_FILE} not found. Run pipeline v25 first.")
        return
    
    d = np.load(CACHE_FILE)
    phases = d["phases"]; ks = d["ks"]; prods = d["prods"]; offsets = d["offsets"]
    print(f"Loaded {len(phases)} events.")

    # CV Align
    ph_cv, ks_cv, pr_cv = [], [], []
    for i in range(len(offsets)-1):
        lo, hi = offsets[i], offsets[i+1]
        mid = lo + (hi-lo)//2
        res = bin_data(phases[lo:mid], ks[lo:mid], prods[lo:mid], BIN_COUNT)
        if not res: continue
        beta, _ = weighted_fit_sine(*res)
        phi = math.atan2(beta[2], beta[1])
        ph_cv.append(phases[mid:hi] + phi)
        ks_cv.append(ks[mid:hi]); pr_cv.append(prods[mid:hi])
            
    ph_f = np.concatenate(ph_cv); ks_f = np.concatenate(ks_cv); pr_f = np.concatenate(pr_cv)
    res = bin_data(ph_f, ks_f, pr_f, BIN_COUNT)
    x, y, err = res
    beta, dchi2 = weighted_fit_sine(x, y, err)
    
    print(f"Final DeltaChi2: {dchi2:.2f}")
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.errorbar(x, y, yerr=err, fmt='o', color='blue', capsize=4, label='Experimental Data')
    th = np.linspace(0, 2*np.pi, 200)
    plt.plot(th, beta[0] + beta[1]*np.sin(th) + beta[2]*np.cos(th), 'orange', lw=2, label=f'Model Fit ($\\Delta\\chi^2={dchi2:.1f}$)')
    plt.axhline(2.0, color='red', ls='--', label='Classical Limit')
    plt.axhline(2.828, color='green', ls=':', label='Tsirelson Bound')
    plt.xlabel('Geometric Phase $\\theta$ (rad)')
    plt.ylabel('Bell Parameter $S$')
    plt.title(f"Geometric Phase Modulation (Weihs Data, $\\kappa=0.1, N={len(ph_f)}$)")
    plt.legend(loc='lower right')
    plt.grid(True, alpha=0.3)
    plt.ylim(min(y)-0.1, max(y)+0.1)
    plt.savefig("weihs_discovery_plot_v25.png", dpi=300)
    print("Saved weihs_discovery_plot_v25.png")

if __name__ == "__main__":
    main()