# X-Theta Lab V2: Advanced Computational Research Stack

X-Theta Lab is a research framework for simulating the effects of curved spacetime on quantum entanglement, specifically focusing on the **relational entanglement anisotropy** and **entanglement lensing** phenomenological model.

## Core Theory (X-Theta V2)

The framework implements a kinematic phenomenological model where the presence of a gravitational potential induces a relational phase $\Phi_{\rm rel}$ between entangled qubits.

### Key Equations

1.  **Relational Phase**:
    $$\Phi_{\rm rel} = r_s \theta \left( \frac{1}{r_{\rm emit}} - \frac{1}{r_{\rm det}} \right)$$
    where $r_s = 2GM/c^2$ is the Schwarzschild radius.

2.  **Relational Generator and Unitary**:
    $$G_{\rm rel} = \frac{1}{2} (X \otimes Y - Y \otimes X)$$
    $$U_{\rm rel}(\Phi) = \exp(i \Phi G_{\rm rel})$$

3.  **Correlation Tensor Anisotropy**:
    The evolved Bell state $|\Psi^- \rangle \rightarrow U_{\rm rel}(\Phi)|\Psi^- \rangle$ results in a correlation tensor $T_{ij} = \text{Tr}(\rho (\sigma_i \otimes \sigma_j))$:
    $$T(\Phi) = \text{diag}[-\cos(2\Phi), -\cos(2\Phi), -1]$$

4.  **Entanglement Lensing**:
    The deformation of the correlation sphere into an ellipsoid. We define two surfaces:
    - **Direct Correlation-Strength Ellipsoid**: Shows actual observable correlation magnitudes.
      Radii: $r_x = r_y = |\cos(2\Phi)|, r_z = 1$
      Equation: $\frac{x^2}{\cos^2(2\Phi)} + \frac{y^2}{\cos^2(2\Phi)} + z^2 = 1$
    - **Dual Response Ellipsoid**: The inverse surface associated with $v^T(T^T T)v = 1$.
      Radii: $r_x = r_y = 1/|\cos(2\Phi)|, r_z = 1$
      Equation: $\cos^2(2\Phi)x^2 + \cos^2(2\Phi)y^2 + z^2 = 1$

5.  **Invariants**:
    $$R_{\Theta} = 3 - \text{Tr}(T^T T) = 2 \sin^2(2\Phi)$$
    $$C(\Phi) = |\cos(2\Phi)| \quad \text{(Concurrence)}$$
    $$S_{\max} = 2\sqrt{1 + C^2}$$

## Scientific Framing

**Note**: The current implementation is a **kinematic phenomenological model**.

1. It does not yet derive the relational generator $G_{\rm rel}$ from a fundamental action $S[g, \Theta]$.
2. It does not yet solve the full covariant surface-selection problem for the relational surface $\Sigma$.
3. It serves to validate the computational signature of curvature-driven entanglement anisotropy.

## Project Structure

- `xtheta/`: Core Python package.
  - `quantum/`: Quantum state evolution, correlation tensor, and CHSH projections.
  - `geometry/`: Schwarzschild phase calculations.
  - `experiments/`: Benchmark scenarios (Micius, GPS, Neutron Star, etc.).
  - `montecarlo/`: Uncertainty propagation for all V2 observables.
  - `visualization/`: Ellipsoid (Strength/Dual) and anisotropy plotting.
- `notebooks/`: Research notebooks for analysis.
  - `01_internal_consistency.ipynb`: Verifies the mathematical heart of the theory.
  - `02_benchmark_scenarios.ipynb`: Computes predictions for real-world and extreme astrophysical cases.
  - `03_entanglement_lensing.ipynb`: Visualizes Direct vs Dual correlation surfaces.
  - `04_monte_carlo_uncertainty.ipynb`: Analyzes sensitivity to experimental uncertainties.
  - `05_concurrence_chsh_geometry.ipynb`: Explores the geometry of CHSH projections.
- `tests/`: Unit tests for all modules.

## Installation

```bash
pip install -r requirements.txt
pip install -e .
```

## Running the Research Stack

You can explore the framework through the provided Jupyter notebooks in the `notebooks/` directory.

### Running notebooks

From the `xtheta-lab` directory, install the package in editable mode:

```bash
pip install -e .
```

Then run notebooks from `xtheta-lab/notebooks`.

If running directly from an IDE, each notebook also includes a small fallback cell that adds the project root to `sys.path`.

## Testing

Run unit tests using `pytest`:

```bash
export PYTHONPATH=$PYTHONPATH:$(pwd)/xtheta-lab
python3 -m pytest xtheta-lab/tests/
```

```powershell
$env:PYTHONPATH = "$env:PYTHONPATH;$PWD\xtheta-lab"
python -m pytest xtheta-lab\tests\
```

## Open Bell/CHSH Data Validation

The framework includes a validation pipeline for real Bell-test datasets. This pipeline computes the CHSH S-statistic and fits an effective phenomenological X-Theta phase ($\Phi_{eff}$) and anisotropy ($R_{\Theta, eff}$).

### Scientific Warning
**Phi_eff is an effective phenomenological parameter only.** Without gravitational path, altitude, curvature, or spacetime-baseline metadata, this is not evidence of spacetime-induced X-Theta holonomy.

### Running Validation

You can run validation on generic CSV data or supported specific datasets (Weihs, Hensen, BIG Bell Test):

```bash
python scripts/run_open_data_validation.py \
  --dataset hensen \
  --data ../bell_open_data.txt \
  --output outputs/open_data/hensen
```

### Interpreting Phi_eff
- $\Phi_{eff} \approx 0$ indicates maximal Bell violation ($S \approx 2\sqrt{2}$) and minimal anisotropy.
- $\Phi_{eff} \approx \pi/4$ indicates a result at the classical Bell limit ($S \approx 2$).
- Structured variation in $\Phi_{eff}$ across datasets with different spacetime baselines would be required for physical X-Theta claims.
