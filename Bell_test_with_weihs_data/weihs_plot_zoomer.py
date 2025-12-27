import numpy as np
import matplotlib.pyplot as plt
import math

# === CONFIG ===
CACHE_FILE = "weihs_mined_data.npz"
BIN_COUNT = 6
PERMUTATION_ROUNDS = 10000  # High precision run

def weighted_linear_fit(X, y, sigma):
    w = 1.0 / (np.maximum(sigma, 1e-12) ** 2)
    sw = np.sqrt(w)
    Xw = X * sw[:, None]
    yw = y * sw
    try:
        XtX = Xw.T @ Xw
        beta = np.linalg.solve(XtX, Xw.T @ yw)
        yhat = X @ beta
        chi2 = np.sum(((y - yhat) / sigma) ** 2)
        return beta, chi2
    except:
        return None, 0

def analyze_vectorized(phases, settings, products):
    # Binning
    bin_indices = ((phases % (2*np.pi)) / (2*np.pi) * BIN_COUNT).astype(int) % BIN_COUNT
    flat_indices = bin_indices * 4 + settings
    
    counts = np.bincount(flat_indices, minlength=BIN_COUNT*4)
    sums = np.bincount(flat_indices, weights=products, minlength=BIN_COUNT*4)
    
    if np.any(counts < 5): return None
    
    means = sums / counts
    errs = 1.0 / np.sqrt(counts) # Analytical Error
    
    means_r = means.reshape((BIN_COUNT, 4))
    errs_r = errs.reshape((BIN_COUNT, 4))
    
    s_vals = means_r[:,0] + means_r[:,1] + means_r[:,2] - means_r[:,3]
    s_errs = np.sqrt(np.sum(errs_r**2, axis=1))
    bin_centers = (np.arange(BIN_COUNT) + 0.5) * (2*math.pi/BIN_COUNT)
    
    # Fit Sine
    X = np.column_stack([np.ones(BIN_COUNT), np.sin(bin_centers), np.cos(bin_centers)])
    beta, chi2 = weighted_linear_fit(X, s_vals, s_errs)
    
    # Fit Flat
    X0 = np.ones((BIN_COUNT, 1))
    _, chi2_0 = weighted_linear_fit(X0, s_vals, s_errs)
    
    return {"s": s_vals, "err": s_errs, "phases": bin_centers, "beta": beta, "delta_chi2": chi2_0 - chi2}

def main():
    print(f"Loading {CACHE_FILE}...")
    try:
        data = np.load(CACHE_FILE)
    except FileNotFoundError:
        print("Error: Run the miner (v02) first to generate the .npz file!")
        return

    phases = data['phases']
    settings = data['settings']
    products = data['products']
    print(f"Loaded {len(phases)} events.")
    
    # 1. Real Analysis
    res = analyze_vectorized(phases, settings, products)
    obs_delta = res['delta_chi2']
    print(f"Observed Delta Chi2: {obs_delta:.2f}")

    # 2. High-Precision Permutation
    print(f"Running {PERMUTATION_ROUNDS} permutations...")
    better = 0
    perm_phases = phases.copy()
    
    for i in range(PERMUTATION_ROUNDS):
        np.random.shuffle(perm_phases)
        pres = analyze_vectorized(perm_phases, settings, products)
        if pres and pres['delta_chi2'] >= obs_delta:
            better += 1
        if (i+1) % 1000 == 0: print(f"..{i+1}")
            
    p_val = (better + 1) / (PERMUTATION_ROUNDS + 1)
    print(f"Final p-value: {p_val:.6f}")

    # 3. Zoomed Plot
    plt.figure(figsize=(10, 6))
    
    # Plot Data
    plt.errorbar(res['phases'], res['s'], yerr=res['err'], fmt='o', color='blue', 
                 ecolor='black', capsize=3, label='Data')
    
    # Plot Fit
    c0, a, b = res['beta']
    th = np.linspace(0, 2*np.pi, 200)
    plt.plot(th, c0 + a*np.sin(th) + b*np.cos(th), color='orange', linewidth=2, 
             label=f'Fit (p={p_val:.5f})')
    
    # Reference Lines
    plt.axhline(2.0, color='red', linestyle='--', alpha=0.5, label='Classical')
    plt.axhline(2.828, color='green', linestyle=':', alpha=0.5, label='Tsirelson')
    
    # Auto-Zoom
    y_min = np.min(res['s'] - res['err']) - 0.05
    y_max = np.max(res['s'] + res['err']) + 0.05
    plt.ylim(y_min, y_max)
    
    plt.xlabel(r'Geometric Phase $\theta$')
    plt.ylabel(r'Bell Parameter $S$')
    plt.title(f"Phase-Dependent Violation (Zoomed)\n$N={len(phases)}, \Delta\chi^2={obs_delta:.1f}$")
    plt.legend(loc='lower right')
    plt.grid(True, alpha=0.3)
    
    plt.savefig("weihs_result_zoomed.png")
    print("Saved: weihs_result_zoomed.png")

if __name__ == "__main__":
    main()