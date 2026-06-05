"""
CSV / Parquet / JSONL / NPZ loaders for Bell-test data.
"""
import pandas as pd
import numpy as np
import json
from pathlib import Path
from typing import Iterator, List

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

def get_loader(path: Path):
    ext = path.suffix.lower()
    if ext == ".csv":
        return iter_csv
    if ext in (".jsonl", ".ndjson"):
        return iter_jsonl
    if ext in (".parquet", ".pq"):
        return iter_parquet
    if ext == ".npz":
        return iter_npz
    raise ValueError(f"Unsupported file extension: {ext}")
