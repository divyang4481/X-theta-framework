import pandas as pd
import numpy as np
import os
import json
from pathlib import Path
from typing import Dict, Optional, List, Iterator
from xtheta.data.schema import BellEventSchema, normalize_outcomes

CAND_SETTING = ["setting", "basis", "angle", "a_setting", "b_setting", "alice_setting", "bob_setting"]
CAND_OUTCOME = ["outcome", "result", "value", "bit", "x", "y", "alice_outcome", "bob_outcome", "a_outcome", "b_outcome"]
CAND_TIME = ["time", "timestamp", "t", "ts"]
CAND_TRIAL = ["trial", "trial_id", "index", "event", "n"]

def _find_col(cols: List[str], candidates: List[str]) -> Optional[str]:
    low = {c.lower(): c for c in cols}
    for cand in candidates:
        for c in cols:
            if cand == c.lower(): return c
        for c in cols:
            if cand in c.lower(): return c
    for cand in candidates:
        if cand in low: return low[cand]
    return None

def infer_columns(df: pd.DataFrame) -> Dict[str, Optional[str]]:
    cols = list(df.columns)
    a_set = None
    b_set = None
    a_out = None
    b_out = None

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

    if a_set is None: a_set = _find_col(cols, ["alice_setting", "a_setting", "setting_a", "a_basis"] + CAND_SETTING)
    if b_set is None: b_set = _find_col(cols, ["bob_setting", "b_setting", "setting_b", "b_basis"] + CAND_SETTING)
    if a_out is None: a_out = _find_col(cols, ["alice_outcome", "a_outcome", "outcome_a", "result_a", "x"] + CAND_OUTCOME)
    if b_out is None: b_out = _find_col(cols, ["bob_outcome", "b_outcome", "outcome_b", "result_b", "y"] + CAND_OUTCOME)

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

def pick_data_files(data_path: Path) -> List[Path]:
    if data_path.is_file(): return [data_path]
    exts = {".csv", ".parquet", ".pq", ".jsonl", ".ndjson", ".npz"}
    files = [p for p in data_path.rglob("*") if p.is_file() and p.suffix.lower() in exts]
    files.sort(key=lambda p: p.stat().st_size, reverse=True)
    return files

def iter_csv(path: Path, chunksize: int) -> Iterator[pd.DataFrame]:
    for chunk in pd.read_csv(path, chunksize=chunksize):
        yield chunk

def iter_jsonl(path: Path, chunksize: int) -> Iterator[pd.DataFrame]:
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

def iter_parquet(path: Path, chunksize: int) -> Iterator[pd.DataFrame]:
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(path)
    for batch in pf.iter_batches(batch_size=chunksize):
        yield batch.to_pandas()

def iter_npz(path: Path, chunksize: int) -> Iterator[pd.DataFrame]:
    data = np.load(path, allow_pickle=True)
    keys = list(data.keys())
    lens = [len(data[k]) for k in keys if hasattr(data[k], "__len__")]
    if not lens:
        raise ValueError(f"NPZ {path} has no array-like content.")
    n = min(lens)
    cols = {k: data[k][:n] for k in keys if hasattr(data[k], "__len__") and len(data[k]) >= n}
    df = pd.DataFrame(cols)
    for start in range(0, n, chunksize):
        yield df.iloc[start:start + chunksize].copy()

def iter_dataframes(path: Path, chunksize: int) -> Iterator[pd.DataFrame]:
    ext = path.suffix.lower()
    if ext == ".csv":
        return iter_csv(path, chunksize)
    if ext in (".jsonl", ".ndjson"):
        return iter_jsonl(path, chunksize)
    if ext in (".parquet", ".pq"):
        return iter_parquet(path, chunksize)
    if ext == ".npz":
        return iter_npz(path, chunksize)
    raise ValueError(f"Unsupported file extension: {ext}")

def get_loader(data_path):
    """
    Returns a loader function for the given data path.
    If data_path is a directory, the loader will iterate through all supported files.
    """
    path = Path(data_path)

    def loader(p, chunksize=200_000):
        p = Path(p)
        if p.is_dir():
            files = pick_data_files(p)
            for file_path in files:
                yield from iter_dataframes(file_path, chunksize=chunksize)
        else:
            yield from iter_dataframes(p, chunksize=chunksize)

    return loader

def load_bell_data(filepath, schema_map=None, **kwargs):
    """
    Load Bell test data from various formats.
    """
    path = Path(filepath)
    df = next(iter_dataframes(path, chunksize=10000000)) # Load first chunk
    if schema_map:
        df = df.rename(columns=schema_map)

    schema = BellEventSchema()
    for attr in vars(schema):
        col = getattr(schema, attr)
        if col not in df.columns:
            df[col] = np.nan

    # Normalize outcomes
    df[schema.alice_outcome] = df[schema.alice_outcome].apply(normalize_outcomes)
    df[schema.bob_outcome] = df[schema.bob_outcome].apply(normalize_outcomes)

    df[schema.source_file] = os.path.basename(filepath)

    return df
