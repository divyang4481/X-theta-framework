# Regression simulations

This folder contains a driver that fits simple physics-motivated models to the three experiments described in the paper and produces figures plus a JSON summary under `paper/build`.

Outputs:
- `paper/build/figs/regression/*.png`
- `paper/build/analysis/regression_summary.json`

Prereqs: create the conda env from `env/environment.yml`.

On Windows PowerShell:

```powershell
conda env create -f env/environment.yml ; conda activate xtheta
python scripts/run_regression_simulations.py
```

If `conda` isn't on PATH in PowerShell, launch "Anaconda Prompt" first or use the Miniforge/Mamba equivalent.
