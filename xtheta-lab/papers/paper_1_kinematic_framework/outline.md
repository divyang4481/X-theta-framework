# Paper 1: Predictive Relational Geometry (X-Theta Kinematic Framework)

## Structure

1. **Motivation**
   - The need for a predictive relational geometry in quantum correlation.
   - Transition from phenomenological fitting to predictive validation.

2. **X-Theta relational phase model**
   - Derivation of the relational phase $\Phi$ from Schwarzschild geometry.
   - Equations for $\Phi$ in terms of $r_s, r_{emit}, r_{det}, \theta$.
   - *Supporting Notebook:* `02_benchmark_scenarios` (Regime definitions).

3. **Correlation tensor deformation**
   - UNITARY evolution of the singlet state under $G_{rel}$.
   - Resulting anisotropic correlation tensor $T(\Phi) = \mathrm{diag}[-\cos(2\Phi), -\cos(2\Phi), -1]$.
   - *Supporting Notebook:* `01_internal_consistency` (Mathematical theorem).

4. **CHSH, concurrence, and anisotropy invariants**
   - Derivation of $C(\Phi) = |\cos(2\Phi)|$ and $S_{\max}(\Phi) = 2\sqrt{1+\cos^2(2\Phi)}$.
   - The anisotropy invariant $R_\Theta = 2\sin^2(2\Phi)$.
   - *Supporting Notebook:* `05_concurrence_chsh_geometry` (Geometric relations).

5. **Benchmark regimes: satellite, neutron star, black hole**
   - Numerical analysis of Earth-orbit (null regime) vs extreme gravity.
   - Uncertainty propagation using Monte Carlo.
   - *Supporting Notebooks:* `02_benchmark_scenarios`, `04_monte_carlo_uncertainty`.

6. **Open-data validation using Hensen**
   - Pipeline validation using Delft loophole-free data.
   - Calculation of $S \approx 2.42$ and effective $\Phi_{eff} \approx 0.41$.
   - Explicit disclaimer on the lack of gravitational metadata.
   - *Supporting Notebook:* `06_hensen_open_data_audit`.

7. **Synthetic recovery and model-comparison methodology**
   - Distinguishing detector phase shift $\delta$ from tensor anisotropy $\Phi$.
   - Residual signature test: $\Delta E \approx \delta \sin(\theta_a - \theta_b)$ vs tensor signature.
   - *Supporting Notebook:* `07_synthetic_phase_recovery`.

8. **Claim classification and falsifiability**
   - 4-level system: Theorem, Simulation, Phenomenological Fit, Physical Prediction.
   - Falsification rules: Rule 1 (Zero anisotropy) and Rule 2 (Inconsistency).
   - *Supporting Notebook:* `08_claim_classification_and_falsification`.

9. **Limitations and future physical experiment**
   - Requirements for metadata-complete Bell tests.
   - Correlation-space lensing and visual diagnostics.
   - *Supporting Notebook:* `03_entanglement_lensing`.

## Notebook Traceability Table

| Section | Title | Primary Notebook |
| :--- | :--- | :--- |
| 1 | Motivation | N/A |
| 2 | X-Theta relational phase model | `02_benchmark_scenarios` |
| 3 | Correlation tensor deformation | `01_internal_consistency` |
| 4 | CHSH, concurrence, and invariants | `05_concurrence_chsh_geometry` |
| 5 | Benchmark regimes | `02_benchmark_scenarios`, `04_monte_carlo_uncertainty` |
| 6 | Open-data validation | `06_hensen_open_data_audit` |
| 7 | Synthetic recovery methodology | `07_synthetic_phase_recovery` |
| 8 | Claim classification and falsifiability | `08_claim_classification_and_falsification` |
| 9 | Limitations and future experiment | `03_entanglement_lensing` |
