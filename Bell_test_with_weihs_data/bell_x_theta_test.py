#!/usr/bin/env python3
"""
bell_x_theta_test.py - Restored original functionality as a wrapper.
"""
import argparse
import json
import math
import numpy as np
import pandas as pd
from pathlib import Path
from collections import deque, Counter
from xtheta.data.schema import ValueEncoder, OutcomeCoder
from xtheta.data.loaders import infer_columns, pick_data_files, iter_dataframes
from xtheta.data.loops import signed_area, is_square_cycle
from xtheta.data.bell_chsh import RunningAB, ThetaBinnedAB

def analyze(data_files, cols_override, chunksize, theta_kappa, theta_bins, max_files):
    alice_set_counts = Counter()
    bob_set_counts = Counter()
    a_out_coder = OutcomeCoder(seen=Counter(), categorical_encoder=ValueEncoder(mapping={}))
    b_out_coder = OutcomeCoder(seen=Counter(), categorical_encoder=ValueEncoder(mapping={}))
    preview_rows_target = 250_000
    inferred_cols = None
    preview_seen = 0

    for fpath in data_files[:max_files]:
        for df in iter_dataframes(fpath, chunksize=min(chunksize, 100_000)):
            if cols_override:
                parts = [p.strip() for p in cols_override.split(",")]
                inferred_cols = {"alice_setting": parts[0], "bob_setting": parts[1], "alice_outcome": parts[2], "bob_outcome": parts[3]}
            else:
                inferred_cols = infer_columns(df)
            a_vals = df[inferred_cols["alice_setting"]].to_numpy()
            b_vals = df[inferred_cols["bob_setting"]].to_numpy()
            for v in a_vals: alice_set_counts[v] += 1
            for v in b_vals: bob_set_counts[v] += 1
            preview_seen += len(df)
            if preview_seen >= preview_rows_target: break
        if preview_seen >= preview_rows_target: break

    top2_a = [v for v, _ in alice_set_counts.most_common(2)]
    top2_b = [v for v, _ in bob_set_counts.most_common(2)]
    a_map = {top2_a[0]: 0, top2_a[1]: 1}
    b_map = {top2_b[0]: 0, top2_b[1]: 1}

    global_ab = RunningAB(count=np.zeros(4, dtype=np.int64), sum_ab=np.zeros(4, dtype=np.int64))
    theta_ab = ThetaBinnedAB(bins=theta_bins, count=np.zeros((theta_bins, 4), dtype=np.int64), sum_ab=np.zeros((theta_bins, 4), dtype=np.int64))
    last_states = deque(maxlen=5)
    last_ab = deque(maxlen=5)
    cw_stats, ccw_stats = [], []
    theta = 0.0
    processed, filtered_out = 0, 0

    for fpath in data_files[:max_files]:
        for df in iter_dataframes(fpath, chunksize=chunksize):
            A = a_out_coder.array_to_pm1(df[inferred_cols["alice_outcome"]].to_numpy())
            B = b_out_coder.array_to_pm1(df[inferred_cols["bob_outcome"]].to_numpy())
            AB = (A.astype(np.int16) * B.astype(np.int16)).astype(np.int8)
            a_set_arr = df[inferred_cols["alice_setting"]].to_numpy()
            b_set_arr = df[inferred_cols["bob_setting"]].to_numpy()

            for i in range(len(df)):
                av, bv = a_set_arr[i], b_set_arr[i]
                if av not in a_map or bv not in b_map:
                    filtered_out += 1; continue
                a, b, ab = a_map[av], b_map[bv], int(AB[i])
                global_ab.update(a, b, ab)
                theta_mod = theta % (2.0 * math.pi)
                bin_id = min(int((theta_mod / (2.0 * math.pi)) * theta_bins), theta_bins - 1)
                theta_ab.update(bin_id, a, b, ab)
                last_states.append((a, b)); last_ab.append(ab)
                if len(last_states) == 5:
                    states = list(last_states)
                    if is_square_cycle(states):
                        area = signed_area(states)
                        orient = 1 if area > 0 else (-1 if area < 0 else 0)
                        cycle_val = float(np.mean(list(last_ab)[1:]))
                        if orient > 0: ccw_stats.append(cycle_val)
                        elif orient < 0: cw_stats.append(cycle_val)
                        if orient != 0 and theta_kappa != 0.0: theta += theta_kappa * orient
                processed += 1

    E = global_ab.expectation()
    def mean_se(x):
        if not x: return float("nan"), float("nan"), 0
        arr = np.asarray(x, dtype=float)
        m = float(arr.mean())
        if len(arr) < 2: return m, float("nan"), len(arr)
        return m, float(arr.std(ddof=1) / math.sqrt(len(arr))), len(arr)

    cw_m, cw_se, cw_n = mean_se(cw_stats)
    ccw_m, ccw_se, ccw_n = mean_se(ccw_stats)
    S_bins, Nmin = theta_ab.chsh_by_bin()

    return {
        "global": {"CHSH_S": float(global_ab.chsh()), "CHSH_S_se": float(global_ab.chsh_se())},
        "x_theta_loop_test": {"cw": {"mean": cw_m, "se": cw_se, "n_cycles": cw_n}, "ccw": {"mean": ccw_m, "se": ccw_se, "n_cycles": ccw_n}},
        "x_theta_conditioned_CHSH": {"S_by_bin": S_bins.tolist(), "min_counts_by_bin": Nmin.tolist()}
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default="bell_x_theta_results.json")
    ap.add_argument("--chunksize", type=int, default=200_000)
    ap.add_argument("--cols", default=None)
    ap.add_argument("--theta-kappa", type=float, default=0.25)
    ap.add_argument("--theta-bins", type=int, default=16)
    ap.add_argument("--max-files", type=int, default=3)
    args = ap.parse_args()
    files = pick_data_files(Path(args.data))
    results = analyze(files, args.cols, args.chunksize, args.theta_kappa, args.theta_bins, args.max_files)
    print(f"\nS = {results['global']['CHSH_S']:.6f} ± {results['global']['CHSH_S_se']:.6f}")
    Path(args.out).write_text(json.dumps(results, indent=2))

if __name__ == "__main__": main()
