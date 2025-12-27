# Mercury Benchmark Note — X–θ Simulation & Math

## Why Mercury?
Mercury’s perihelion precession is the cleanest “no-miracles” filter for any modification of central-force gravity.

- Newton + known planetary perturbations leave an anomalous residual ≈ **43 arcseconds/century**.
- GR reproduces this via a distinctive weak-field correction that behaves like an effective **\(\propto r^{-3}\)** term in the orbit dynamics.

**Goal (not to “beat GR”):**
Use Mercury as a **constraint + calibration** test so X–θ remains falsifiable and numerically disciplined.

---

## Baseline: GR target formula
A standard expression for GR perihelion advance per orbit is

\[
\Delta\phi_{\mathrm{GR}} \approx \frac{6\pi GM}{c^2\,a\,(1-e^2)}
\]

where
- \(G\) gravitational constant, \(M\) Solar mass,
- \(a\) semi-major axis, \(e\) eccentricity,
- \(c\) speed of light.

Useful conversion:

- arcsec per orbit: \(\Delta\phi_{\text{arcsec/orbit}} = \Delta\phi\;\times\;\frac{180}{\pi}\times 3600\)
- arcsec per century: multiply by (orbits/century)

Orbits/century approx:
\[
N_{\text{orbits/century}} \approx \frac{100\;\text{years}}{T_{\text{Mercury}}}
\]

---

## X–θ hook: how the fiber can generate an orbit correction
In X–θ, the “fiber energy” appears as

\[
H_{\text{fiber}}=\frac{1}{2I}\,(p_\theta-q_\theta A_\theta)^2
\]

with
\[
p_\theta = I\dot{\theta}+q_\theta A_\theta
\]

If \(A_\theta\) depends on position \(r\), then \(H_{\text{fiber}}\) contributes an **effective potential** term to base motion.

Expanding:
\[
H_{\text{fiber}}=\frac{p_\theta^2}{2I} - \frac{q_\theta p_\theta}{I}A_\theta(r) + \frac{q_\theta^2}{2I}A_\theta(r)^2
\]

So the induced correction is generically
\[
\delta V_{X\theta}(r)\;\sim\;-\frac{q_\theta p_\theta}{I}\,A_\theta(r)\; +\; \frac{q_\theta^2}{2I}\,A_\theta(r)^2
\]

### Matching the GR “shape”
Mercury’s precession is sensitive to corrections that act like
\[
\delta V(r) \propto -\frac{\kappa}{r^3}
\]

Therefore, a practical (phenomenological) constraint is:

> Choose/derive \(A_\theta(r)\) (and possibly \(I(r)\)) so that the **leading** X–θ correction behaves like \(\delta V_{X\theta}(r)\approx -\kappa/r^3\) in the weak-field regime.

Then \(\kappa\) is either:
- fixed by matching \(\Delta\phi\) to the GR benchmark (calibration), or
- bounded tightly by ephemeris residuals (exclusion).

---

## Two-test protocol (Mercury hammer)

### Test A — Short-horizon ephemeris residual (days → months)
**Input:** JPL/Horizons initial state \((\mathbf r_0,\mathbf v_0)\) at epoch \(t_0\)

1) Integrate the orbit under:
- Newtonian gravity baseline
- Newton + X–θ correction (your \(A_\theta(r)\) model)

2) Compare to Horizons position at \(t_1\):
\[
\varepsilon(t_1)=\|\mathbf r_{\text{sim}}(t_1)-\mathbf r_{\text{Horizons}}(t_1)\|
\]

**Rule:** if \(\varepsilon\) blows up beyond numerical/ephemeris tolerance, the chosen coupling is ruled out.


### Test B — Long-horizon perihelion drift (many orbits)
1) Detect perihelion events via local minima of \(r(t)=\|\mathbf r(t)\|\)
2) Record perihelion angle \(\phi_n\) at each event.
3) Measure per-orbit advance:
\[
\Delta\phi_{\text{sim}} = \phi_{n+1}-\phi_n
\]
4) Convert to arcsec/century and compare to the target residual (~43″/century).

**Rule:** X–θ must match the benchmark within uncertainty, or be excluded.

---

## Practical simulation details

### Numerical integrator
- Symplectic (preferred for long runs): Velocity Verlet / Leapfrog
- Or standard high-accuracy: RK45 with tight tolerances (good for quick validation)

### Force form (if you encode \(\delta V = -\kappa/r^3\))
If
\[
\delta V(r)=-\frac{\kappa}{r^3}
\]
then
\[
\delta F_r(r) = -\frac{d}{dr}\delta V = -\frac{d}{dr}\left(-\kappa r^{-3}\right)= -\left(3\kappa r^{-4}\right)= -\frac{3\kappa}{r^4}
\]
Direction: radial (along \(-\hat{r}\) if attractive).

This gives you a clean “knob” \(\kappa\) to sweep and fit.

---

## How this ties back to our X–θ intuition
We already use mixed curvature as a physical effect (e.g., a cross-Hall drift from gradients in \(A_\theta\)).

Mercury is the same idea, but in the **radial channel**:
- If \(A_\theta\) varies with \(r\), the fiber sector contributes a tiny, structured correction to central motion.
- Perihelion precession is exquisitely sensitive to exactly that kind of correction.

So Mercury acts as a reality-anchor: **either X–θ reduces to GR-like weak-field behavior (within bounds), or it gets crushed by data.**

---

## Implementation hooks (your notebooks)
- `Mercury_vs_NASA_Challenge.ipynb` → use for Horizons comparison + perihelion extraction
- `Cross_Hall_Drift.ipynb` → reuse the “Aθ gradient → force” scaffolding
- `Draft.md` → source of the fiber-Hamiltonian definitions and narrative

---

## Minimal deliverables for a paper section
1) Define \(H_{\text{fiber}}\) and the induced \(\delta V_{X\theta}(r)\)
2) State the GR benchmark \(\Delta\phi_{\mathrm{GR}}\)
3) Describe the two tests (ephemeris residual + perihelion drift)
4) Report a sweep over \(\kappa\) (or over a parameterized \(A_\theta(r)\))
5) Conclude with: **allowed parameter band** or **exclusion**

