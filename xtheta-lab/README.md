# X-Theta Lab V2: Advanced Computational Research Stack

X-Theta Lab is a research framework for simulating the effects of curved spacetime on quantum entanglement, specifically focusing on the **relational entanglement anisotropy** and **entanglement lensing** kinematic phenomenological model.

## Core Theory (X-Theta V2)

The framework implements a kinematic phenomenological model where the presence of a gravitational potential induces a relational phase $\Phi_{\rm rel}$ between entangled qubits.

Unitary relational evolution preserves state purity and norm, while entanglement concurrence follows the curve $C(\Phi) = |\cos(2\Phi)|$. The resulting CHSH variation is understood as a detector-geometry projection of the underlying correlation tensor anisotropy.

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

4.  **Invariants**:
    $$R_{\Theta} = 3 - \text{Tr}(T^T T) = 2 \sin^2(2\Phi)$$
    $$C(\Phi) = |\cos(2\Phi)| \quad \text{(Concurrence)}$$
    $$S_{\max} = 2\sqrt{1 + C^2}$$

## Scientific Framing

**Note**: The current implementation is a **kinematic phenomenological model**.

1. It serves to validate the computational signature of curvature-driven entanglement anisotropy.
2. It does not yet derive the relational generator $G_{\rm rel}$ from a fundamental action.
3. Earth-orbit predictions are many orders below current practical sensitivity. Compact-object scenarios are theoretical stress tests.

## Project Structure

- `xtheta/`: Core Python package.
  - `quantum/`: Quantum state evolution, correlation tensor, spectrum, and purity diagnostics.
  - `geometry/`: Schwarzschild phase calculations.
  - `experiments/`: Benchmark scenarios, random CHSH landscape, and open data validation.
  - `data/`: Open data adapter layer (CSV, Parquet, NPZ loaders).
  - `montecarlo/`: Uncertainty propagation for all V2 observables.
  - `visualization/`: Ellipsoid (Strength/Dual), anisotropy, and landscape plotting.
- `scripts/`:
  - `generate_paper1_artifacts.py`: One-command generation of all research artifacts.
  - `run_open_data_validation.py`: Validate pipeline against real datasets.
- `notebooks/`: Research notebooks for interactive analysis.
- `tests/`: Extensive unit test suite.

## Installation

```bash
pip install -r requirements.txt
pip install -e .
```

## Running the Research Stack

### 1. Generate Paper 1 Artifacts
To regenerate all figures, data files, and the validation summary:
```bash
python scripts/generate_paper1_artifacts.py
```
Outputs will be saved in the `outputs/` directory.

### 2. Run Tests
```bash
pytest tests/ -v
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
  --data path/to/dataset.txt \
  --output outputs/open_data/hensen
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
