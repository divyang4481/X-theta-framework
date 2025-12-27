#!/usr/bin/env python3
"""
bell_x_theta_test.py

Streaming Bell/CHSH analysis + X–θ path-memory tests on large (~GB) datasets.

Works with:
- Single file OR a directory containing one or more data files.
- CSV (chunked), Parquet (via pyarrow), NPZ (loads arrays), JSONL (basic).

Core outputs:
1) Standard CHSH S-statistic + per-setting E(a,b)
2) Loop orientation effect (CW vs CCW) in setting-space square loops
3) CHSH conditioned on a causal toy holonomy theta (binned)

Usage examples (Windows):
  python bell_x_theta_test.py --data "Physics\\Bell_data"
  python bell_x_theta_test.py --data "Physics\\Bell_data\\events.csv" --chunksize 200000
  python bell_x_theta_test.py --data "Physics\\Bell_data" --cols "alice_setting,bob_setting,alice_outcome,bob_outcome,timestamp"

If your columns aren't auto-detected, pass --cols (see below).
"""

from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from collections import deque, Counter
from typing import Dict, Iterable, Iterator, List, Optional, Tuple, Union

import numpy as np

try:
    import pandas as pd
except ImportError:
    raise SystemExit("This script requires pandas. Install with: pip install pandas")


# -----------------------------
# Helpers: schema inference
# -----------------------------

CAND_SETTING = ["setting", "basis", "angle", "a_setting", "b_setting", "alice_setting", "bob_setting"]
CAND_OUTCOME = ["outcome", "result", "value", "bit", "x", "y", "alice_outcome", "bob_outcome", "a_outcome", "b_outcome"]
CAND_TIME = ["time", "timestamp", "t", "ts"]
CAND_TRIAL = ["trial", "trial_id", "index", "event", "n"]


def _find_col(cols: List[str], candidates: List[str]) -> Optional[str]:
    low = {c.lower(): c for c in cols}
    for cand in candidates:
        for c in cols:
            if cand == c.lower():
                return c
        for c in cols:
            if cand in c.lower():
                return c
    # last attempt: exact candidate in low mapping
    for cand in candidates:
        if cand in low:
            return low[cand]
    return None


def infer_columns(df: "pd.DataFrame") -> Dict[str, Optional[str]]:
    cols = list(df.columns)
    # Try to infer Alice/Bob settings/outcomes
    a_set = None
    b_set = None
    a_out = None
    b_out = None

    # Strong matches first
    for c in cols:
        cl = c.lower()
        if "alice" in cl and any(k in cl for k in ["setting", "basis", "angle"]):
            a_set = c
        if "bob" in cl and any(k in cl for k in ["setting", "basis", "angle"]):
            b_set = c
        if "alice" in cl and any(k in cl for k in ["outcome", "result", "value", "bit", "x"]):
            a_out = c
        if "bob" in cl and any(k in cl for k in ["outcome", "result", "value", "bit", "y"]):
            b_out = c

    # Fallback: guess pairs by generic names
    if a_set is None:
        a_set = _find_col(cols, ["alice_setting", "a_setting", "setting_a", "a_basis"] + CAND_SETTING)
    if b_set is None:
        b_set = _find_col(cols, ["bob_setting", "b_setting", "setting_b", "b_basis"] + CAND_SETTING)
    if a_out is None:
        a_out = _find_col(cols, ["alice_outcome", "a_outcome", "outcome_a", "result_a", "x"] + CAND_OUTCOME)
    if b_out is None:
        b_out = _find_col(cols, ["bob_outcome", "b_outcome", "outcome_b", "result_b", "y"] + CAND_OUTCOME)

    ts = _find_col(cols, CAND_TIME)
    trial = _find_col(cols, CAND_TRIAL)

    return {
        "alice_setting": a_set,
        "bob_setting": b_set,
        "alice_outcome": a_out,
        "bob_outcome": b_out,
        "timestamp": ts,
        "trial_id": trial,
    }


# -----------------------------
# Encoders (stream-safe)
# -----------------------------

@dataclass
class ValueEncoder:
    """Stream-safe factor encoder: maps arbitrary values to small ints as they appear."""
    mapping: Dict[object, int]
    next_id: int = 0

    def encode(self, v) -> int:
        if v in self.mapping:
            return self.mapping[v]
        self.mapping[v] = self.next_id
        self.next_id += 1
        return self.mapping[v]

    def encode_array(self, arr: np.ndarray) -> np.ndarray:
        out = np.empty(len(arr), dtype=np.int16)
        for i, v in enumerate(arr):
            out[i] = self.encode(v)
        return out

    def top_k_by_frequency(self, counts: Counter, k: int) -> List[int]:
        # returns encoder ids for top-k values
        items = counts.most_common(k)
        return [self.mapping[val] for val, _ in items]


@dataclass
class OutcomeCoder:
    """
    Convert outcomes to ±1.
    Handles:
      - already ±1
      - 0/1 -> -1/+1
      - True/False -> -1/+1
      - arbitrary categorical -> factorize then map first->-1 second->+1 (only sensible for 2-category)
    """
    seen: Counter
    categorical_encoder: ValueEncoder

    def to_pm1(self, v) -> int:
        # normalize common numeric / boolean forms
        if v is True:
            return 1
        if v is False:
            return -1

        # numpy scalars -> python types
        if isinstance(v, (np.generic,)):
            v = v.item()

        if isinstance(v, (int, np.integer)):
            if v == 1:
                return 1
            if v == 0:
                return -1
            if v == -1:
                return -1

        # float-ish forms
        if isinstance(v, (float, np.floating)):
            if abs(v - 1.0) < 1e-12:
                return 1
            if abs(v) < 1e-12:
                return -1
            if abs(v + 1.0) < 1e-12:
                return -1

        # categorical fallback
        self.seen[v] += 1
        cid = self.categorical_encoder.encode(v)
        # map first category -> -1, second -> +1, others alternate
        return 1 if (cid % 2 == 1) else -1

    def array_to_pm1(self, arr: np.ndarray) -> np.ndarray:
        out = np.empty(len(arr), dtype=np.int8)
        for i, v in enumerate(arr):
            out[i] = self.to_pm1(v)
        return out


# -----------------------------
# Setting-space loop detector
# -----------------------------

def hamming01(u: Tuple[int, int], v: Tuple[int, int]) -> int:
    return (u[0] != v[0]) + (u[1] != v[1])


def signed_area(vertices: List[Tuple[int, int]]) -> float:
    # vertices assumed closed: last == first
    area = 0.0
    for (x1, y1), (x2, y2) in zip(vertices[:-1], vertices[1:]):
        area += (x1 * y2 - x2 * y1)
    return 0.5 * area


def is_square_cycle(states: List[Tuple[int, int]]) -> bool:
    """
    states length 5: v0,v1,v2,v3,v4 with v4 == v0.
    Must visit all 4 vertices exactly once (except closure), and move along edges.
    """
    if len(states) != 5:
        return False
    if states[0] != states[-1]:
        return False
    uniq = states[:-1]
    if len(set(uniq)) != 4:
        return False
    # must be edges (hamming distance 1) between consecutive vertices
    for u, v in zip(states[:-1], states[1:]):
        if hamming01(u, v) != 1:
            return False
    return True


@dataclass
class RunningAB:
    # aggregates for AB by (a,b)
    count: np.ndarray  # shape (4,)
    sum_ab: np.ndarray  # shape (4,)

    def update(self, a: int, b: int, ab: int):
        idx = a * 2 + b
        self.count[idx] += 1
        self.sum_ab[idx] += ab

    def expectation(self) -> np.ndarray:
        E = np.zeros(4, dtype=float)
        nonzero = self.count > 0
        E[nonzero] = self.sum_ab[nonzero] / self.count[nonzero]
        return E

    def chsh(self) -> float:
        E = self.expectation()
        # E00 + E01 + E10 - E11
        return float(E[0] + E[1] + E[2] - E[3])

    def chsh_se(self) -> float:
        """
        Standard error using Var(AB)=1-E^2 for ±1 AB.
        """
        E = self.expectation()
        se_terms = []
        for i in range(4):
            n = self.count[i]
            if n <= 0:
                se_terms.append(0.0)
                continue
            var = max(0.0, 1.0 - float(E[i] ** 2))
            se_terms.append(var / n)
        # S = E00 + E01 + E10 - E11 => variances add
        return float(math.sqrt(se_terms[0] + se_terms[1] + se_terms[2] + se_terms[3]))


@dataclass
class ThetaBinnedAB:
    # aggregates for AB by theta-bin and (a,b)
    bins: int
    count: np.ndarray  # shape (bins,4)
    sum_ab: np.ndarray  # shape (bins,4)

    def update(self, bin_id: int, a: int, b: int, ab: int):
        idx = a * 2 + b
        self.count[bin_id, idx] += 1
        self.sum_ab[bin_id, idx] += ab

    def chsh_by_bin(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Returns:
          S_bins: (bins,)
          N_min: (bins,) min count across the 4 setting pairs (useful for sanity)
        """
        S = np.zeros(self.bins, dtype=float)
        Nmin = np.zeros(self.bins, dtype=int)
        for k in range(self.bins):
            cnt = self.count[k]
            sab = self.sum_ab[k]
            E = np.zeros(4, dtype=float)
            nz = cnt > 0
            E[nz] = sab[nz] / cnt[nz]
            S[k] = E[0] + E[1] + E[2] - E[3]
            Nmin[k] = int(cnt.min())
        return S, Nmin


# -----------------------------
# Data readers (streamed)
# -----------------------------

def iter_csv(path: Path, chunksize: int) -> Iterator["pd.DataFrame"]:
    for chunk in pd.read_csv(path, chunksize=chunksize):
        yield chunk


def iter_jsonl(path: Path, chunksize: int) -> Iterator["pd.DataFrame"]:
    # simple JSON Lines streaming
    buf = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            buf.append(json.loads(line))
            if len(buf) >= chunksize:
                yield pd.DataFrame(buf)
                buf.clear()
    if buf:
        yield pd.DataFrame(buf)


def iter_parquet(path: Path, chunksize: int) -> Iterator["pd.DataFrame"]:
    try:
        import pyarrow.parquet as pq
    except ImportError:
        raise SystemExit("Parquet file detected but pyarrow is not installed. Install with: pip install pyarrow")

    pf = pq.ParquetFile(path)
    for batch in pf.iter_batches(batch_size=chunksize):
        yield batch.to_pandas()


def iter_npz(path: Path, chunksize: int) -> Iterator["pd.DataFrame"]:
    data = np.load(path, allow_pickle=True)
    # Heuristic: expect arrays named like alice_setting, bob_setting, alice_outcome, bob_outcome
    keys = list(data.keys())
    # build DataFrame from all arrays of same length
    lens = [len(data[k]) for k in keys if hasattr(data[k], "__len__")]
    if not lens:
        raise SystemExit(f"NPZ {path} has no array-like content.")
    n = min(lens)
    cols = {k: data[k][:n] for k in keys if hasattr(data[k], "__len__") and len(data[k]) >= n}
    df = pd.DataFrame(cols)
    # yield in chunks
    for start in range(0, n, chunksize):
        yield df.iloc[start:start + chunksize].copy()


def pick_data_files(data_path: Path) -> List[Path]:
    if data_path.is_file():
        return [data_path]

    exts = {".csv", ".parquet", ".pq", ".jsonl", ".ndjson", ".npz"}
    files = [p for p in data_path.rglob("*") if p.is_file() and p.suffix.lower() in exts]
    if not files:
        raise SystemExit(f"No supported data files found under: {data_path}")

    # process largest first (often the main dataset)
    files.sort(key=lambda p: p.stat().st_size, reverse=True)
    return files


def iter_dataframes(path: Path, chunksize: int) -> Iterator["pd.DataFrame"]:
    ext = path.suffix.lower()
    if ext == ".csv":
        return iter_csv(path, chunksize)
    if ext in (".jsonl", ".ndjson"):
        return iter_jsonl(path, chunksize)
    if ext in (".parquet", ".pq"):
        return iter_parquet(path, chunksize)
    if ext == ".npz":
        return iter_npz(path, chunksize)
    raise SystemExit(f"Unsupported file extension: {ext}")


# -----------------------------
# Main analysis
# -----------------------------

def analyze(
    data_files: List[Path],
    cols_override: Optional[str],
    chunksize: int,
    theta_kappa: float,
    theta_bins: int,
    max_files: int,
) -> Dict[str, object]:

    # Encoders and counters (stream stable)
    alice_set_enc = ValueEncoder(mapping={})
    bob_set_enc = ValueEncoder(mapping={})
    alice_set_counts = Counter()
    bob_set_counts = Counter()

    a_out_coder = OutcomeCoder(seen=Counter(), categorical_encoder=ValueEncoder(mapping={}))
    b_out_coder = OutcomeCoder(seen=Counter(), categorical_encoder=ValueEncoder(mapping={}))

    # We will first collect enough to decide which two settings are "0/1" for each side.
    # Strategy: one light pass over the beginning, but without reading all data into memory.
    # Then do the full streaming pass with filtering to those top-2 settings.
    preview_rows_target = 250_000

    inferred_cols: Optional[Dict[str, Optional[str]]] = None
    preview_seen = 0

    for fpath in data_files[:max_files]:
        for df in iter_dataframes(fpath, chunksize=min(chunksize, 100_000)):
            if cols_override:
                # user provided explicit mapping
                parts = [p.strip() for p in cols_override.split(",")]
                # allowed forms:
                # 4 cols: a_set,b_set,a_out,b_out
                # 5 cols: +timestamp
                # 6 cols: +timestamp,+trial_id
                if len(parts) < 4:
                    raise SystemExit("--cols must include at least 4 comma-separated names.")
                inferred_cols = {
                    "alice_setting": parts[0],
                    "bob_setting": parts[1],
                    "alice_outcome": parts[2],
                    "bob_outcome": parts[3],
                    "timestamp": parts[4] if len(parts) >= 5 else None,
                    "trial_id": parts[5] if len(parts) >= 6 else None,
                }
            else:
                inferred_cols = infer_columns(df)

            a_set = inferred_cols["alice_setting"]
            b_set = inferred_cols["bob_setting"]
            a_out = inferred_cols["alice_outcome"]
            b_out = inferred_cols["bob_outcome"]

            if not all([a_set, b_set, a_out, b_out]):
                raise SystemExit(
                    "Could not auto-detect required columns. "
                    "Provide them via --cols 'alice_setting,bob_setting,alice_outcome,bob_outcome[,timestamp,trial_id]'."
                )

            # update setting frequency counts using encoded raw values
            a_vals = df[a_set].to_numpy()
            b_vals = df[b_set].to_numpy()
            for v in a_vals:
                alice_set_counts[v] += 1
            for v in b_vals:
                bob_set_counts[v] += 1

            preview_seen += len(df)
            if preview_seen >= preview_rows_target:
                break
        if preview_seen >= preview_rows_target:
            break

    # Pick top-2 settings for each side
    top2_a_vals = [v for v, _ in alice_set_counts.most_common(2)]
    top2_b_vals = [v for v, _ in bob_set_counts.most_common(2)]
    if len(top2_a_vals) < 2 or len(top2_b_vals) < 2:
        raise SystemExit("Need at least 2 distinct settings for Alice and Bob to compute CHSH.")

    # fixed mapping: most common -> 0, second -> 1
    a_map = {top2_a_vals[0]: 0, top2_a_vals[1]: 1}
    b_map = {top2_b_vals[0]: 0, top2_b_vals[1]: 1}

    # Aggregators
    global_ab = RunningAB(count=np.zeros(4, dtype=np.int64), sum_ab=np.zeros(4, dtype=np.int64))
    theta_ab = ThetaBinnedAB(
        bins=theta_bins,
        count=np.zeros((theta_bins, 4), dtype=np.int64),
        sum_ab=np.zeros((theta_bins, 4), dtype=np.int64),
    )

    # Loop detection state (streaming)
    last_states: deque[Tuple[int, int]] = deque(maxlen=5)
    last_ab: deque[int] = deque(maxlen=5)

    # Loop stats (small arrays kept in RAM)
    cw_stats: List[float] = []
    ccw_stats: List[float] = []

    # theta state (causal: theta affects next trials after loop completes)
    theta = 0.0

    processed = 0
    filtered_out = 0

    for fpath in data_files[:max_files]:
        for df in iter_dataframes(fpath, chunksize=chunksize):
            # columns (same inference as preview)
            if not inferred_cols:
                inferred_cols = infer_columns(df)
            a_set = inferred_cols["alice_setting"]
            b_set = inferred_cols["bob_setting"]
            a_out = inferred_cols["alice_outcome"]
            b_out = inferred_cols["bob_outcome"]

            # Extract arrays
            a_set_arr = df[a_set].to_numpy()
            b_set_arr = df[b_set].to_numpy()
            a_out_arr = df[a_out].to_numpy()
            b_out_arr = df[b_out].to_numpy()

            # Convert outcomes to ±1
            A = a_out_coder.array_to_pm1(a_out_arr)
            B = b_out_coder.array_to_pm1(b_out_arr)
            AB = (A.astype(np.int16) * B.astype(np.int16)).astype(np.int8)

            # Process row-wise (streaming logic + loop detection)
            for i in range(len(df)):
                av = a_set_arr[i]
                bv = b_set_arr[i]

                # Filter to top-2 settings each side (CHSH requires 2x2)
                if av not in a_map or bv not in b_map:
                    filtered_out += 1
                    continue

                a = a_map[av]
                b = b_map[bv]
                ab = int(AB[i])

                # Update global CHSH
                global_ab.update(a, b, ab)

                # Update theta-binned CHSH (theta BEFORE any loop-closure update on this trial)
                theta_mod = theta % (2.0 * math.pi)
                bin_id = int((theta_mod / (2.0 * math.pi)) * theta_bins)
                if bin_id == theta_bins:
                    bin_id = theta_bins - 1
                theta_ab.update(bin_id, a, b, ab)

                # Loop detection update
                last_states.append((a, b))
                last_ab.append(ab)

                if len(last_states) == 5:
                    states = list(last_states)
                    if is_square_cycle(states):
                        # orientation via signed area of the polygon
                        # treat vertices as points (a,b) in plane
                        area = signed_area(states)
                        orient = 1 if area > 0 else (-1 if area < 0 else 0)

                        # cycle statistic: mean AB on the 4 "step" trials after the first vertex
                        # (uses indices 1..4 so the start isn't double-counted)
                        cycle_val = float(np.mean(list(last_ab)[1:]))

                        if orient > 0:
                            ccw_stats.append(cycle_val)
                        elif orient < 0:
                            cw_stats.append(cycle_val)

                        # causal holonomy proxy: update theta AFTER detecting completed loop
                        if orient != 0 and theta_kappa != 0.0:
                            theta += theta_kappa * orient

                processed += 1

    # Compute global CHSH results
    E = global_ab.expectation()
    S = global_ab.chsh()
    S_se = global_ab.chsh_se()

    # Loop effect summary
    def mean_se(x: List[float]) -> Tuple[float, float, int]:
        if not x:
            return (float("nan"), float("nan"), 0)
        arr = np.asarray(x, dtype=float)
        m = float(arr.mean())
        if len(arr) < 2:
            return (m, float("nan"), len(arr))
        se = float(arr.std(ddof=1) / math.sqrt(len(arr)))
        return (m, se, len(arr))

    cw_m, cw_se, cw_n = mean_se(cw_stats)
    ccw_m, ccw_se, ccw_n = mean_se(ccw_stats)

    diff = float(ccw_m - cw_m) if (cw_n > 0 and ccw_n > 0) else float("nan")
    diff_se = float(math.sqrt((cw_se if not math.isnan(cw_se) else 0.0) ** 2 +
                              (ccw_se if not math.isnan(ccw_se) else 0.0) ** 2)) if (cw_n > 1 and ccw_n > 1) else float("nan")
    z = float(diff / diff_se) if (diff_se and not math.isnan(diff_se) and diff_se > 0) else float("nan")

    # Theta-binned CHSH
    S_bins, Nmin = theta_ab.chsh_by_bin()
    S_bins_max = float(np.nanmax(S_bins)) if len(S_bins) else float("nan")
    S_bins_min = float(np.nanmin(S_bins)) if len(S_bins) else float("nan")

    results = {
        "data_files_used": [str(p) for p in data_files[:max_files]],
        "rows_processed_filtered_space": int(processed),
        "rows_filtered_out_non_2x2_settings": int(filtered_out),
        "settings_top2_alice": [repr(top2_a_vals[0]), repr(top2_a_vals[1])],
        "settings_top2_bob": [repr(top2_b_vals[0]), repr(top2_b_vals[1])],
        "global": {
            "E00_E01_E10_E11": [float(E[0]), float(E[1]), float(E[2]), float(E[3])],
            "CHSH_S": float(S),
            "CHSH_S_se": float(S_se),
        },
        "x_theta_loop_test": {
            "definition": "Detect 5-trial square cycles in setting-space; compare mean(AB) for CCW vs CW loops.",
            "cw": {"mean": cw_m, "se": cw_se, "n_cycles": cw_n},
            "ccw": {"mean": ccw_m, "se": ccw_se, "n_cycles": ccw_n},
            "diff_ccw_minus_cw": diff,
            "diff_se": diff_se,
            "z_score_approx": z,
        },
        "x_theta_conditioned_CHSH": {
            "theta_kappa": float(theta_kappa),
            "theta_bins": int(theta_bins),
            "S_by_bin": [float(x) for x in S_bins.tolist()],
            "min_counts_by_bin": [int(x) for x in Nmin.tolist()],
            "S_bin_max": S_bins_max,
            "S_bin_min": S_bins_min,
            "note": "If 'conditioning on theta' is meaningful, you may see structured variation across bins "
                    "(but beware of small N per bin).",
        },
        "notes": [
            "This script uses a *toy* causal theta updated only when a full square loop closes.",
            "If your dataset has explicit timing/space channels, we can upgrade theta to a physically-defined loop functional.",
        ],
    }

    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="Path to Physics/Bell_data directory or a single file.")
    ap.add_argument("--out", default="bell_x_theta_results.json", help="Output JSON filename.")
    ap.add_argument("--chunksize", type=int, default=200_000, help="Chunk size for streaming CSV/Parquet/JSONL.")
    ap.add_argument("--cols", default=None,
                    help="Comma-separated column names: a_setting,b_setting,a_outcome,b_outcome[,timestamp,trial_id]. "
                         "Use this if auto-detection fails.")
    ap.add_argument("--theta-kappa", type=float, default=0.25,
                    help="Theta increment per detected loop (radians). 0 disables theta update.")
    ap.add_argument("--theta-bins", type=int, default=16, help="Bins for theta-conditioned CHSH.")
    ap.add_argument("--max-files", type=int, default=3, help="Max files to process (largest first).")
    args = ap.parse_args()

    data_path = Path(args.data)
    files = pick_data_files(data_path)

    results = analyze(
        data_files=files,
        cols_override=args.cols,
        chunksize=args.chunksize,
        theta_kappa=args.theta_kappa,
        theta_bins=args.theta_bins,
        max_files=args.max_files,
    )

    # Print a human-readable summary
    g = results["global"]
    lt = results["x_theta_loop_test"]
    ct = results["x_theta_conditioned_CHSH"]

    print("\n=== Global CHSH ===")
    print("E00,E01,E10,E11 =", g["E00_E01_E10_E11"])
    print(f"S = {g['CHSH_S']:.6f}  ± {g['CHSH_S_se']:.6f} (SE)")

    print("\n=== X–θ Loop Test (Setting-space square loops) ===")
    print("CW : mean =", lt["cw"]["mean"], "SE =", lt["cw"]["se"], "n =", lt["cw"]["n_cycles"])
    print("CCW: mean =", lt["ccw"]["mean"], "SE =", lt["ccw"]["se"], "n =", lt["ccw"]["n_cycles"])
    print("diff (CCW-CW) =", lt["diff_ccw_minus_cw"], "z ≈", lt["z_score_approx"])

    print("\n=== X–θ Conditioned CHSH (theta bins) ===")
    print("theta_kappa =", ct["theta_kappa"], "bins =", ct["theta_bins"])
    print("S_bin_min =", ct["S_bin_min"], "S_bin_max =", ct["S_bin_max"])
    print("min count per setting-pair per bin (sanity) =", min(ct["min_counts_by_bin"]))

    out_path = Path(args.out)
    out_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\nWrote: {out_path.resolve()}\n")


if __name__ == "__main__":
    main()
