# notebook2

This folder is for the **latest** simulation notebooks used to generate:
- figures saved into `paper2/figs/`
- tables/CSVs saved into `paper2/data/`
- quick analysis plots and sanity checks

## Notebooks
- `00_environment_check.ipynb`: verifies imports and paths
- `01_xtheta_cross_hall_drift.ipynb`: runs the baseline X–θ drift sim and writes a CSV + figure
- `02_analysis_and_citations.ipynb`: loads the CSV, reproduces plots, and checks BibTeX keys exist in `paper2/refs.bib`

## Conventions
- Keep outputs under `paper2/` so the paper stays reproducible.
- Use BibTeX keys from `paper2/refs.bib` inside notebook markdown (e.g. `Aharonov1959`).
