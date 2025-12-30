#!/usr/bin/env python3
r"""
run in order :
 python x_theta_report_v1.py --files   out\x_theta_v3_none_N2000000_T64.npz   out\x_theta_v3_hard_N2000000_T64.npz

"""

import argparse
import numpy as np
import matplotlib.pyplot as plt


def load_npz(path):
    z = np.load(path, allow_pickle=True)
    return {
        "path": path,
        "theta": z["theta_grid"],
        "deltaS": z["deltaS"],
        "se_boot": z["se_deltaS_boot"],
        "p_perm": float(z["p_perm"]),
        "spec": z["spec_obs"],
        "kept": int(z["kept"]),
        "cfg": z["cfg"].item() if z["cfg"].dtype == object else z["cfg"],
    }


def dominant_harmonic(spec):
    if len(spec) <= 1:
        return 0
    return int(np.argmax(spec[1:]) + 1)


def summary_row(d):
    cfg = d["cfg"]
    systematics = cfg["systematics"]
    hol = float(cfg["holonomy_amp"])
    kdom = dominant_harmonic(d["spec"])
    # weighted deltaS (approx): average over theta
    dS_w = float(np.average(d["deltaS"]))
    stat = float(np.max(np.abs(d["deltaS"])))
    return systematics, hol, d["kept"], dS_w, stat, d["p_perm"], kdom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--files", nargs="+", required=True, help="npz files to compare")
    ap.add_argument("--title", default="X–Θ Report")
    args = ap.parse_args()

    runs = [load_npz(p) for p in args.files]

    # ---- Summary table ----
    print("\n=== SUMMARY ===")
    print(
        "systematics  holonomy  kept        mean(ΔS)    max|ΔS|     p_perm       k_dom"
    )
    for r in runs:
        systematics, hol, kept, dS_w, stat, p, kdom = summary_row(r)
        print(
            f"{systematics:10s} {hol:7.3f}  {kept:9d}  {dS_w:10.5f}  {stat:10.5f}  {p:11.6g}  {kdom:5d}"
        )

    # ---- Plot ΔS(θ) with error bands ----
    plt.figure()
    for r in runs:
        cfg = r["cfg"]
        label = f"{cfg['systematics']}, hol={cfg['holonomy_amp']}"
        th = r["theta"]
        y = r["deltaS"]
        se = r["se_boot"]
        plt.plot(th, y, label=label)
        plt.fill_between(th, y - se, y + se, alpha=0.2)
    plt.xlabel("θ")
    plt.ylabel("ΔS(θ) = S_CW(θ) - S_CCW(θ)")
    plt.title(args.title + " — ΔS(θ)")
    plt.legend()
    plt.grid(True, alpha=0.2)

    # ---- Plot FFT magnitude ----
    plt.figure()
    for r in runs:
        cfg = r["cfg"]
        label = f"{cfg['systematics']}, hol={cfg['holonomy_amp']}"
        spec = r["spec"]
        k = np.arange(len(spec))
        plt.plot(k[1:], spec[1:], label=label)  # skip DC
    plt.xlabel("harmonic k")
    plt.ylabel("|FFT(ΔS)|")
    plt.title(args.title + " — FFT spectrum (ΔS)")
    plt.legend()
    plt.grid(True, alpha=0.2)

    plt.show()


if __name__ == "__main__":
    main()
