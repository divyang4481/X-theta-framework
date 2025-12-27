# X–θ + Topological Transformer Research Canvas (Chronological Roadmap)

This canvas consolidates **our working theory**, the **simulation plan**, and **experiment / validation ideas**, in the same order we iteratively discovered and improved them.

---

## 0) The Core Hypothesis (the seed idea)

**Working claim:** Reality (or at least our effective model of it) can be represented as a trivial bundle

\[
Q = X \times S^1
\]

- \(X\): ordinary space(-time) coordinates.
- \(\theta \in S^1\): a compact “fiber” phase.

**Key move:** introduce a **connection** (gauge-like field) \(A_\mu\) that couples \(\theta\) to motion in \(X\). The connection can have nonzero curvature (“Berry curvature-like”):

\[
F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu
\]

**Consequences (predictions):** not just local forces, but **path-memory** effects (holonomy / geometric phase).

---

## 1) First falsifiable effect: Cross-Hall Drift (Simulation Module #1)

### 1.1 Minimal prediction
A gradient in the fiber connection can induce a transverse drift in \(X\):

\[
\Delta x_\perp \propto (\nabla A_\theta)\, \dot{\theta}
\]

Interpretation: even if the main flow is chaotic or symmetric, a “hidden” fiber drive can bias transport.

### 1.2 First simulation (synthetic)
- Particles move in a 2D flow field \(u(x,t)\) (turbulent-ish or random).
- Each particle carries a phase \(\theta\) evolving with \(\dot{\theta} = \Omega + \text{noise}\).
- Add coupling via \(A_\theta(x)\) with a controlled gradient.

**Success criterion:** the ensemble COM shifts in \(x\) (or \(y\)) compared to control.

### 1.3 PINN verification attempt (Physics-Informed Neural Net)
We tried a PINN-style network to validate that the drift is consistent with the PDE/ODE constraint (rather than a numerical artifact).

**Improvement:** use the PINN as a *consistency check* (residual minimization), not as the main simulator.

---

## 2) Upgrade: Macroscopic Coherence / “Levitation” Thought Experiment (Module #2)

### 2.1 Scaling argument
If many microscopic phases align (cohere), small drifts add:

\[
F_{\text{total}} \approx \sum_{i=1}^N f_i \sim N\, f_{\text{atom}}
\]

**Working claim:** coherence is the lever; thermal randomness cancels the effect.

### 2.2 What this really means (as a testable idea)
- Not “magic propulsion.”
- A claim about **collective phase locking** generating a measurable bias in transport or acceleration.

**Experiment direction (safer, nearer-term):** look for drift/transport bias in systems where coherence is already engineered:
- cold atoms / BEC-like setups,
- driven oscillators / metamaterials,
- coupled resonator lattices,
- photonic waveguides with synthetic dimensions.

---

## 3) Real-data grounding: Turbulent Radiative Layer / THE WELL-type data (Module #3)

### 3.1 Motivation
Synthetic turbulence is too forgiving. Real datasets have:
- multi-scale structure,
- nontrivial vorticity,
- genuine transport/mixing complexity.

### 3.2 “Chiral Topological Resonance” anomaly
We saw a surprising behavior:
- A “Hidden” (unreversed) mode outperformed a “Perfect” mode at a specific winding \(N\) (notably \(N=-10\) in our notes).

**Interpretation (working):** the hidden coupling aligned with an intrinsic chirality/vorticity in the data, producing a resonant transport bias.

### 3.3 Concrete validation checklist
To treat this like science (not vibes), we need:
- fixed seeds,
- ablation of \(A_\theta\) terms,
- sweep \(N\) across a wide range,
- multiple flow regimes / datasets,
- report uncertainty (mean ± std over runs).

**Metrics to report:**
- drift magnitude vs time,
- mixing/dispersion measures,
- correlation of drift with vorticity sign,
- robustness to noise in \(\theta\).

---

## 4) Holonomy as the “Memory”: steering-angle analogy (Conceptual Bridge)

We reframed holonomy using a human analogy:
- walk/drive a loop in \(X\), return to same location,
- but \(\theta\) changes by a net amount.

\[
\Delta \theta_{\text{loop}} = \oint A \cdot d\ell \quad (\text{mod } 2\pi)
\]

Repeat the loop \(N\) times:

\[
\Delta \theta_N \approx N\, \Delta \theta_1 \; (\text{mod } 2\pi)
\]

**Why this matters:** it gives a clean experimental signature: path-dependence independent of speed.

---

## 5) Bell/CHSH angle: “Geometry generates probabilities” (Module #4)

### 5.1 The discipline rule
No cheating: simulated outcomes cannot use real outcomes.

### 5.2 Case-B protocol (probabilistic generation)
For each trial with detector settings \(a, b\):
1) compute a geometric phase/holonomy \(\Delta\phi\) (the *only* allowed dependency is settings + hidden \(\theta\)-structure).
2) set probabilities from geometry:

\[
P(++) = \cos^2(\Delta\phi)
\]

3) sample \(A_{sim}, B_{sim}\) from those probabilities.
4) compute CHSH \(S\) on simulated data and compare.

### 5.3 What we’re actually testing
Whether a geometric bundle model can reproduce correlations while respecting the no-cheating constraint.

**Key outputs:**
- \(S\) value,
- per-setting correlation table \(E(a,b)\),
- calibration curve between \(\Delta\phi\) and observed correlations.

---

## 6) “A_MODE” / data-field hygiene (Operational Improvement)

When moving to real datasets, we must document every field we use.

**Example:** if a dataset includes something like `A_MODE` (common in scientific sims), we must define it precisely:
- what physical mode/variable it encodes,
- whether it is a control flag, a boundary condition, or a label,
- whether it leaks information (data leakage risk).

**Rule:** any feature that can encode the answer indirectly must be treated like a forbidden shortcut.

---

## 7) Parallel track: Topological Transformer for 16K context on 6GB VRAM (Module #5)

This is the *engineering sibling* of the physics story: use topology/geometry to compress parameters and extend context.

### 7.1 Constraint
- VRAM: ~6GB.
- Target: **16K context length**.
- Need small params + low KV-cache cost.
- Use KD (knowledge distillation) from a strong teacher.

### 7.2 Key idea A: Fourier / winding-mode heads
Heads lie on a circle. Instead of separate matrices per head, define head weights as a low-mode Fourier series:

\[
W_Q(h) = \sum_{m=0}^{M-1}\big(A_m\cos(m\theta_h) + B_m\sin(m\theta_h)\big),\quad M \ll H
\]

Same for \(W_K(h), W_V(h)\). This gives:
- parameter sharing across heads,
- smooth head-to-head structure,
- controllable capacity via \(M\).

**Drop-in module goal:** replace per-head projection matrices with a mode generator.

### 7.3 Key idea B: Q/K/V as derivatives of a shared base
Instead of storing three independent projections, define a base map \(W\) and derive others:

- **Shared-base:** \(W_Q = W\)
- **Derivative-tied:** \(W_K = \partial_{\theta} W\), \(W_V = \partial_{\theta\theta} W\) (discrete analogs acceptable)
- or **low-rank tying:** \(W_Q = U\Lambda_Q V^T\), \(W_K = U\Lambda_K V^T\), \(W_V = U\Lambda_V V^T\)

This can cut parameters while preserving expressive structure.

### 7.4 Key idea C: Complex embeddings (phase-rich) without param explosion
Use complex-valued representations (or real+imag packing) to encode phase-like features.

We must keep it practical:
- represent complex vector \(z = x + i y\) as concatenated \([x;y]\),
- use complex RoPE-like rotations as structured transforms.

### 7.5 Practical recipe to hit 16K in 6GB
- Use FlashAttention or memory-efficient attention.
- Reduce KV-cache via grouped-query attention (GQA) or multi-query attention (MQA).
- Prefer smaller \(d\_model\) + more layers only if stable.
- Distill with a schedule: start CE-heavy, then ramp KD.
- Initialize embeddings from teacher (projection) to stabilize early steps.

**Success metrics:** perplexity vs baseline, long-context pass rate, inference throughput, VRAM peak.

---

## 8) New areas to check (where the theory might “bite”)

These are domains where a **geometric phase / synthetic dimension** story is already plausible and measurable.

### 8.1 Physics / systems
- Photonics: waveguides, ring resonators, synthetic dimensions.
- Cold atoms: shaken lattices, artificial gauge fields.
- Metamaterials: chiral structures, topological edge modes.
- Fluid mixing: vorticity-coupled transport biases (our current strongest lead).

### 8.2 ML / computation
- Long-context LLMs: topology-inspired weight tying + KV reduction.
- Robustness: holonomy as an inductive bias for “path-memory” tasks.
- Compression: mode-based parameterization (Fourier heads) as a general technique.

### 8.3 “Bridge” experiments
- Any lab setup where you can:
  1) drive a cyclic path in control space,
  2) measure a net phase-like shift independent of speed,
  3) test additivity over repeated cycles.

---

## 9) Roadmap: improvement order (what to do next, in sequence)

### Phase 1 — Make the evidence clean
1) Re-run Cross-Hall drift with strict ablations and fixed seeds.
2) Port to real turbulence datasets (The Well / radiative layer) with documented fields.
3) Sweep \(N\) and check for “resonant” windows; report uncertainty.

### Phase 2 — Turn anomaly into mechanism
4) Quantify chirality coupling: correlate drift sign with vorticity sign.
5) Test alternative definitions of \(\Delta\phi\) (holonomy functions) and see which generalize.
6) Create a minimal “toy theorem” explaining why a particular \(N\) would resonate.

### Phase 3 — Falsification-grade protocols
7) Publish a strict falsification checklist: what would conclusively kill the claim.
8) Run the Bell/CHSH Case-B simulation with a locked protocol, publish code + seeds.

### Phase 4 — Engineering spinout (Topological Transformer)
9) Implement Fourier-head module + train tiny student with KD.
10) Add Q/K/V tying; measure quality vs parameter drop.
11) Target 16K context with MQA/GQA + efficient attention; publish memory/throughput plots.

---

## 10) The “One-page pitch” (for collaborators)

We are exploring a unifying concept:
- **Geometric phase (holonomy) = path memory**
- In physics: it can bias transport/mixing in turbulent flows.
- In ML: it can compress attention parameters and extend context.

The program is intentionally falsifiable:
- ablations,
- real-data validation,
- strict no-leak protocols (especially for CHSH-style tests).

---

## Appendix: Definitions we must keep consistent

- **Fiber phase:** \(\theta \in S^1\)
- **Connection:** \(A\) (a 1-form) that defines how \(\theta\) changes along paths in \(X\)
- **Curvature:** \(F = dA\)
- **Holonomy:** \(\Delta\theta = \oint A\cdot d\ell\) (mod \(2\pi\))
- **Winding:** repeated loops add phase: \(\Delta\theta_N \approx N\Delta\theta_1\) (mod \(2\pi\))

---

### Notes (living document)
As we run new validations, we will add:
- exact dataset names + versions,
- exact hyperparameters + seeds,
- links to code folders / notebooks,
- tables of results for each \(N\) sweep.

