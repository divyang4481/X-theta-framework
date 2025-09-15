# FRW Bounce vs. GR Singularity — Pre/Post Analysis Notes

**Run timestamp:** 2025-09-12T19:05:43  


This document summarizes the setup, assumptions, equations, and results for the
side‑by‑side simulations comparing a regulated X–θ cosmology to unregulated GR.

---

## 1. Objective (plain language)
Show that adding a compact internal phase θ induces an effective **regulator**
so that the total energy density saturates at a critical scale \( \rho_c \),
turning the would‑be singular collapse into a **non‑singular bounce**.
Compare with plain GR, which has no such cutoff.

---

## 2. Equations (minimal set)
- **Components:** matter \( \propto a^{-3} \), radiation \( \propto a^{-4} \),
  θ‑stiff \( \propto a^{-6} \).  
  Total density:  
  \[
  \rho(a)=\rho_m a^{-3}+\rho_r a^{-4}+\rho_s a^{-6}.
  \]
- **Hubble law (unregulated GR):**
  \[
  H^2 = \tfrac{8\pi G}{3}\,\rho(a).
  \]
- **Hubble law (X–θ regulated):**
  \[
  H^2 = \tfrac{8\pi G}{3}\,\rho(a)\Big(1-\tfrac{\rho(a)}{\rho_c}\Big).
  \]
  The bracketed term is the **regulator**; as \( \rho\to\rho_c \), \( H\to 0 \) and the
  universe bounces instead of hitting a singularity.
- **Dynamics:** \( \dot a = a H \).
- **Curvature proxy (flat FRW):**  
  \[
  K(t) = 12\left[\left(\tfrac{\ddot a}{a}\right)^2 + \left(\tfrac{\dot a}{a}\right)^4\right].
  \]
  This grows without bound in the GR collapse; it stays finite through the
  bounce in X–θ.

---

## 3. Parameters (dimensionless units)

- Absorb \( \tfrac{8\pi G}{3}=1 \) into units.  
- Critical density (regulator): \( \rho_c = 1.0 \).  
- Present‑day densities at \( a=1 \):  
  \( \rho_m=0.3 \), \( \rho_r=10^{-4} \), \( \rho_s=10^{-5} \).  
- Integration: \( a_0=1 \), \( dt=2\times10^{-3} \), \( t_{\max}=30 \).  
- Floors/caps to keep numerics stable: \( a_{\min}=10^{-6} \), \( a_{\max}=10 \).

These values are arbitrary but convenient; the qualitative result is robust.

---

## 4. Methods (repro recipe)
1. Define \( \rho(a) \) and the two Hubble functions (with/without regulator).
2. Integrate \( \dot a = a H \) forward and backward in time (Euler step is fine
   for this demo; use RK4 if you want cleaner derivatives).
3. Compute \( K(t) \) from finite differences for \( \dot a \) and \( \ddot a \).
4. Plot **a(t)**, **H(t)**, and **\(K(t)\)** for both models on the same axes.

---

## 5. Results (headline numbers)

- **Regulated bounce:** \( a_{\min} \approx 0.669567 \) at \( t \approx -1.012 \).  
- **Unregulated GR:** reaches \( a \approx 10^{-6} \) (numerical floor).

**Figures (from the notebook):**
- *Scale factor a(t):* X–θ (solid) shows a clean bounce; GR (dashed) heads for zero.  
- *Hubble rate:* X–θ pins \( H\to 0 \) at the bounce; GR diverges during collapse.  
- *Curvature proxy K(t):* finite across the X–θ bounce; skyrockets for GR as \( a\to 0 \).

---

## 6. Diagnostics & sanity checks
- Changing \( \rho_c \) shifts \( a_{\min} \) but preserves nonsingularity.  
- Removing the \( a^{-6} \) θ‑stiff term delays when the regulator becomes relevant,
  but the \( 1-\rho/\rho_c \) factor alone still enforces a bounce once \( \rho\to\rho_c \).
- Using a higher‑order integrator improves the smoothness of \( K(t) \) near the bounce
  (finite‑difference second derivatives are numerically noisy).

---

## 7. Interpretation (for undergrads)
- **Why GR fails:** with only \( H^2\propto\rho \), nothing stops \( \rho \) from
  blowing up as \( a\to 0 \).  
- **Why X–θ helps:** compact \( \theta \) means the energy density can’t exceed
  a maximum \( \rho_c \). The factor \( 1-\rho/\rho_c \) is like an “elastic floor”:
  when you hit it, expansion rate goes to zero and reverses—**a bounce**.

---

## 8. Next steps
- Map \( \rho_c \) to microscopic θ parameters \( (R_\theta, I, g_\theta) \) via your
  chosen derivation (Born–Infeld‑style saturation or holonomy).  
- Port to **black‑hole interiors** (Kantowski–Sachs ODEs) to show an
  \( r\to r_{\min}>0 \) bounce instead of \( r=0 \).  
- Export plots as PDF/PGF for LaTeX and include this note as a short
  **Methods** subsection.

---

## 9. Repro cell (pseudocode)
```python
# define rho(a), H(a) with/without (1 - rho/rho_c)
# integrate dot{a} = a H
# compute K(t) from finite differences
# plot a(t), H(t), K(t) for both models
```
