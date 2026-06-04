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
    *   **Direct Correlation-Strength Ellipsoid**: Shows actual observable correlation magnitudes.
        Radii: $r_x = r_y = |\cos(2\Phi)|, r_z = 1$
        Equation: $\frac{x^2}{\cos^2(2\Phi)} + \frac{y^2}{\cos^2(2\Phi)} + z^2 = 1$
    *   **Dual Response Ellipsoid**: The inverse surface associated with $v^T(T^T T)v = 1$.
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
```

## Running the Research Stack

You can explore the framework through the provided Jupyter notebooks in the `notebooks/` directory.

## Testing

Run unit tests using `pytest`:

```bash
export PYTHONPATH=$PYTHONPATH:$(pwd)/xtheta-lab
python3 -m pytest xtheta-lab/tests/
```
