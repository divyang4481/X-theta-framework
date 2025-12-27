# X-Theta Manim Animations

This module holds reusable Manim scenes that visualize the core simulations in `src/sim`. Each scene produces a short explanatory animation that can be rendered at multiple quality levels (low for previews, high for publication cuts).

## Folder layout

- `scenes.py` – Manim scene definitions. Each scene imports the corresponding simulation code and converts its output into animated primitives.
- `__init__.py` – Declares this directory as a Python package so Manim can import the scenes when launched from the project root.
- `assets/` *(optional)* – Place custom SVGs, transparency plates, or audio beds here if needed for future scenes. Leave empty for now.

## Installation

1. Ensure the base environment is active (`conda activate xtheta` if you use the provided `env/environment.yml`).
2. Install Manim's dependencies (already captured in `requirements.txt` / `environment.yml`). On Windows you also need FFmpeg in your `PATH` (Manim will prompt if it is missing).

```
pip install -r requirements.txt  # or conda env update -f env/environment.yml
```

## Rendering a scene

From the repository root, run Manim and point it at the new module:

```
manim -pql animations/manim/scenes.py DoubleSlitThetaScene
```

Replace `DoubleSlitThetaScene` with `XThetaDriftScene` (centroid evolution) or `XThetaMonteCarloScene` (Monte Carlo scatter) as needed. Use quality flags (`-pql` preview, `-pqh` high quality) per the Manim docs.

Each command writes MP4s to `media/` (Manim default). Copy the resulting clips into `paper/figs/` or any presentation deck as required.
