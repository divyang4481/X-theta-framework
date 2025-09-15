great call — let’s lock down a **repo-wide policy** for *where we generate* data/figures and *where we show* them, then give you drop-in cells to enforce it.

## The policy (simple + future-proof)

* **Development (quick previews):** save to repo-root `figs/` (and display inline in the notebook).
* **Paper artifacts (versioned & citable):** save to `paper/figs`, `paper/tables`, `paper/data`, and optional text to `paper/analysis/`.
* **LaTeX includes:** `\graphicspath{{paper/figs/}}` and `\input{paper/tables/...}` from `main.tex` or `sections/*.tex`.

You’ll choose the target with a single switch: `SAVE_TARGET = "paper"` or `"dev"`.

---

### Cell A — Code (Repo-aware paths + switches)

```python
# PATH POLICY: one switch to pick where artifacts go, dev vs paper
from pathlib import Path
import os

SAVE_TARGET = "paper"   # "paper" or "dev"
SHOW_IN_NOTEBOOK = True # always True for now; set False for headless runs

# Detect repo root (folder name anchor)
PROJECT_ROOT_NAME = "X-theta-framework"
repo_root = Path.cwd().resolve()
for parent in [repo_root, *repo_root.parents]:
    if parent.name == PROJECT_ROOT_NAME:
        repo_root = parent
        break

# Targets
DEV_FIG_DIR = repo_root / "figs"
PAPER_FIG_DIR = repo_root / "paper" / "figs"
PAPER_TAB_DIR = repo_root / "paper" / "tables"
PAPER_DAT_DIR = repo_root / "paper" / "data"
PAPER_TXT_DIR = repo_root / "paper" / "analysis"

# Select active dirs
if SAVE_TARGET == "paper":
    FIG_DIR = PAPER_FIG_DIR
    TAB_DIR = PAPER_TAB_DIR
    DAT_DIR = PAPER_DAT_DIR
    TXT_DIR = PAPER_TXT_DIR
else:
    FIG_DIR = DEV_FIG_DIR
    TAB_DIR = repo_root / "tables"      # local dev tables (optional)
    DAT_DIR = repo_root / "data"        # local dev data (optional)
    TXT_DIR = repo_root / "others"      # local dev notes (optional)

# Ensure directories exist
for d in [FIG_DIR, TAB_DIR, DAT_DIR, TXT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

print("Repo root         :", repo_root)
print("Saving target     :", SAVE_TARGET)
print("FIG_DIR (active)  :", FIG_DIR)
print("TAB_DIR (active)  :", TAB_DIR)
print("DAT_DIR (active)  :", DAT_DIR)
print("TXT_DIR (active)  :", TXT_DIR)
```

---

### Cell B — Code (Helpers: figure/table/data paths + exporters)

```python
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import json

def fig_base(name: str):
    """Base path for a figure without extension."""
    return (FIG_DIR / name).with_suffix("")

def export_fig(name: str, png=True, svg=True, dpi=180, transparent=False):
    """
    Save current figure to FIG_DIR with both PNG & SVG (optional).
    Returns dict {ext: fullpath}.
    """
    base = fig_base(name)
    base.parent.mkdir(parents=True, exist_ok=True)
    out = {}
    if png:
        p = base.with_suffix(".png")
        plt.savefig(p, dpi=dpi, transparent=transparent)
        out["png"] = str(p.resolve())
    if svg:
        p = base.with_suffix(".svg")
        plt.savefig(p, transparent=transparent)
        out["svg"] = str(p.resolve())
    # echo
    for ext, fp in out.items():
        p = Path(fp)
        print(f"Saved {ext.upper():>3}: {p}  size={p.stat().st_size if p.exists() else 0}")
    return out

def table_path(name: str, ext="tex"):
    """Path for LaTeX/CSV tables in TAB_DIR."""
    p = (TAB_DIR / name).with_suffix(f".{ext}")
    p.parent.mkdir(parents=True, exist_ok=True)
    return p

def data_path(name: str, ext="csv"):
    """Path for data artifacts in DAT_DIR."""
    p = (DAT_DIR / name).with_suffix(f".{ext}")
    p.parent.mkdir(parents=True, exist_ok=True)
    return p

def write_table(df: pd.DataFrame, name: str, tex=True, csv=True, floatfmt=6):
    """Write a DataFrame to CSV and/or LaTeX in TAB_DIR and return paths."""
    paths = {}
    if csv:
        p_csv = table_path(name, "csv")
        df.to_csv(p_csv, index=False)
        paths["csv"] = str(p_csv.resolve())
        print("Saved CSV:", p_csv)
    if tex:
        p_tex = table_path(name, "tex")
        df.round(floatfmt).to_latex(p_tex, index=False)
        paths["tex"] = str(p_tex.resolve())
        print("Saved TEX:", p_tex)
    return paths

def write_json(obj, name: str):
    p = data_path(name, "json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)
    print("Saved JSON:", p)
    return str(p.resolve())
```

---

### Cell C — Code (Where to **generate** vs **show**)

```python
# HOW WE USE IT IN PRACTICE

# Generate: numbers/tables/data always written to DAT_DIR/TAB_DIR.
# Show: figures always displayed in the notebook if SHOW_IN_NOTEBOOK=True,
#       AND saved to FIG_DIR (dev or paper depending on SAVE_TARGET).
#
# Example stubs (replace with your actual variables from earlier cells):

# 1) AB baseline figure
import numpy as np
phi_fracs = np.linspace(-1.0, 1.0, 13)
# assume numerical_em_phase_from_fraction exists from earlier
phase_num = np.array([numerical_em_phase_from_fraction(f) for f in phi_fracs])
phase_the = 2*np.pi*phi_fracs

plt.figure(figsize=(6,4))
plt.plot(phi_fracs, phase_num, 'o-', label='Numerical (EM only)')
plt.plot(phi_fracs, phase_the, '--', label='Theory 2π f')
plt.xlabel(r'$\Phi/\Phi_0$'); plt.ylabel('phase (rad)')
plt.title('AB phase vs flux fraction (θ off)')
plt.legend()
export_fig("ab_baseline")
if SHOW_IN_NOTEBOOK:
    plt.show()
else:
    plt.close()

# 2) AB baseline table
import pandas as pd
df_ab = pd.DataFrame({
    "phi_frac_EM": phi_fracs,
    "phase_numeric_rad": phase_num,
    "phase_theory_rad": phase_the,
    "residual_rad": phase_num - phase_the,
})
write_table(df_ab, "ab_baseline_table")  # saves to TAB_DIR; CSV + TEX
write_json({
    "slope_target": float(2*np.pi),
    "residual_mean": float((phase_num - phase_the).mean()),
    "residual_std": float((phase_num - phase_the).std()),
}, "ab_baseline_metrics")  # DAT_DIR
```

---

### Cell D — Markdown (Where to **show** in the paper)

````markdown
## Where results are shown (paper integration)

- **Figures** live in `paper/figs` (or `figs/` in dev mode).
  In LaTeX, add once:
  ```tex
  \graphicspath{{paper/figs/}}
````

and include with:

```tex
\begin{figure}[t]
  \centering
  \includegraphics[width=0.7\linewidth]{ab_baseline}
  \caption{AB phase vs flux fraction with $\theta$ off. Slope $\approx 2\pi$.}
  \label{fig:ab_baseline}
\end{figure}
```

* **Tables** live in `paper/tables`. Include with:

  ```tex
  \input{paper/tables/ab_baseline_table}
  ```

* **Data** live in `paper/data` for archiving / reproducibility (CSV/JSON). Cite in text as:
  “See data in `paper/data/ab_baseline_metrics.json`.”

````

---

### Cell E — Code (Quick sanity: list active outputs)
```python
from pathlib import Path

def _lsdir(path: Path, exts=(".png",".svg",".tex",".csv",".json")):
    path = Path(path)
    found = []
    for ext in exts:
        found.extend(sorted(p for p in path.glob(f"*{ext}")))
    return found

print("FIG_DIR contents:")
for p in _lsdir(FIG_DIR, (".png",".svg")):
    print(" ", p.name, p.stat().st_size, "bytes")

print("\nTAB_DIR contents:")
for p in _lsdir(TAB_DIR, (".tex",".csv")):
    print(" ", p.name, p.stat().st_size, "bytes")

print("\nDAT_DIR contents:")
for p in _lsdir(DAT_DIR, (".json",".csv")):
    print(" ", p.name, p.stat().st_size, "bytes")
````

---

### What this gives you

* A **single switch** (`SAVE_TARGET`) to decide **where to generate** artifacts.
* All **publication** assets land in `paper/{figs,tables,data,analysis}`.
* All **dev** assets land in `{figs,tables,data,others}` at repo root.
* Clear LaTeX include patterns so you know exactly **where to show** them in the paper.

Flip `SAVE_TARGET` to `"paper"` when you’re producing camera-ready outputs; keep `"dev"` during exploratory runs.
