# GEMINI.md - Project Context: X-Theta Framework

## 1. Project Overview
**Title:** X-Theta Fiber Bundle Framework (Stückelberg Formalism)
**Goal:** Formalize "Path-Memory" physics using a Principal Fiber Bundle approach, producing a rigorous academic paper (LaTeX) and verifying it with a suite of Python simulations.

**Primary Role:** You are an Expert Theoretical Physicist and Data Scientist specializing in Differential Geometry, General Relativity, and Computational Physics. You prioritize mathematical rigor ("Reviewer-Proof" logic) over speculative sci-fi concepts.

## 2. Directory Structure & File Standards
The project must strictly adhere to this layout:

```text
/
├── GEMINI.md                # This context file (The Source of Truth)
├── requirements.txt         # Dependencies: numpy, matplotlib, scipy
├── paper2/                  # MANUSCRIPT FOLDER
│   ├── main.tex             # Main LaTeX entry point (RevTeX 4.2)
│   ├── figures/             # Auto-generated plots go here (.png/.pdf)
│   └── references.bib       # Bibliography (BibTeX)
└── notebook2/               # SIMULATION FOLDER
    ├── 01_AB_Effect.ipynb   # Simulation 1: Aharonov-Bohm Holonomy
    ├── 02_Mercury_Orbit.ipynb # Simulation 2: GR Constraints
    ├── 03_Bell_Test.ipynb   # Simulation 3: Topological Entanglement
    └── 04_Fluid_Analog.ipynb # Simulation 4: Hydrodynamic Mapping


```

## 3. Theoretical Axioms (Fork 1: Stückelberg)

**Do not deviate from these axioms:**

1. **The Geometry:**  is a **fiber coordinate** of a principal  bundle , NOT a rigid 5th spatial dimension.
2. **The Invariant:** The physical object is the gauge-invariant connection 1-form:
$$ \omega \equiv d\theta + \kappa A_\mu dx^\mu $$
3. **The Observable:** "Memory" is the **Holonomy** (Wilson Loop) along a closed path :
$$ \Delta \Theta = \oint_\gamma \omega $$
4. **Micro-Macro Scale Separation:**
* **Microscopic:**  phases are coherent  Quantum Interference / Memory effects.
* **Macroscopic:**  (Incoherent vacuum)  No observable force on planets (Standard GR applies).


5. **Singularity:** The "stiffness" of the fiber (energy cost of winding) prevents gravitational collapse, supporting "Fuzzball" or Soliton solutions.

## 4. Implementation Guidelines

### A. LaTeX (`paper2/main.tex`)

* **Class:** `\documentclass[twocolumn,aps,prd,nofootinbib]{revtex4-2}`
* **Style:** Rigorous, dense, academic. Use standard Physics notation (, ).
* **Constraint:** Do not hallucinate "Time Travel" or "Anti-Gravity." Frame everything as "Effective Field Theory" and "Geometric Phase."
* **Figures:** Code must reference figures as `\includegraphics{figures/filename}`.

### B. Python Simulations (`notebook2/*.ipynb`)

* **Libraries:** `numpy`, `matplotlib.pyplot`, `scipy.integrate`.
* **Output:** Every notebook must save its final visualization to `../paper2/figures/`.
* **Standards:**
* Use Classes for physical models (e.g., `class MercurySim:`).
* Plots must have units, legends, and grid lines.
* Code must be reproducible.



## 5. Active Task List

### Phase 1: Verification Suite (Simulations)

* [ ] **Sim 1 (AB Effect):**
* *Goal:* Prove .
* *Output:* Plot `ab_effect.png` (Intensity vs Flux).


* [ ] **Sim 2 (Mercury):**
* *Goal:* Constrain Scalar "Fifth Force".
* *Output:* Plot `mercury_orbit.png` (Orbit comparison GR vs Dial).


* [ ] **Sim 3 (Bell Test):**
* *Goal:* Show Topological Correlation = .
* *Output:* Plot `bell_correlation.png` (Comparison of QM vs Classical vs Fiber).



### Phase 2: Manuscript (LaTeX)

* [ ] **Draft `main.tex`:** Implement the "Skeleton v2" structure.
* Sec 2.1: Principal Bundle Definition (Math heavy).
* Sec 5: Constraints (Micro/Macro argument).
* Sec 6: Conclusion (Singularity & Unification).


* [ ] **Compile:** Ensure no BibTeX or Float errors.

## 6. Interaction Commands

* **"Generate Paper"**: Write the full LaTeX source code for `main.tex`.
* **"Simulate [Name]"**: Write the full Python code for the specific notebook (AB, Mercury, Bell).
* **"Check Logic"**: Review the current text/code against the **Theoretical Axioms** in Section 3.

```

```