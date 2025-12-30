import math
import numpy as np

try:
    from scipy.special import ellipe  # parameter m
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False
    import mpmath as mp

def _ellipe_param_m(m: np.ndarray) -> np.ndarray:
    """
    Complete elliptic integral of the second kind E(m) with parameter m:
    E(m) = ∫_0^{π/2} sqrt(1 - m sin^2 t) dt
    """
    if _HAS_SCIPY:
        return ellipe(m)
    # mpmath fallback (slow)
    out = np.empty_like(m, dtype=np.float64)
    for i, mi in enumerate(m):
        out[i] = float(mp.ellipe(mi))
    return out

def G_and_H(lam: np.ndarray, zeta: np.ndarray, tau: float):
    """
    Implements paper's G(λ,ζ;τ) and H(λ,ζ;τ) as shown around Eqs (17)-(18).
    Notes:
      - The paper uses λ, ζ as the top two strain eigenvalues.
      - All arrays are elementwise.
    """
    lam = np.asarray(lam, dtype=np.float64)
    zeta = np.asarray(zeta, dtype=np.float64)

    # Common exponent combos
    a = (2*lam + zeta) * tau
    b = (zeta - lam) * tau

    # Parameter for elliptic E(·) as rendered in the PDF text layer
    # If you later verify from the PDF equation image that the parameter differs,
    # you only need to change this one line.
    m = 1.0 - np.exp(2.0*(2*zeta + lam)*tau)  # as in PDF text layer

    E = _ellipe_param_m(m)

    # Eq (17): <ln f>_{θ,φ} = G(λ,ζ;τ)
    G = (
        lam*tau
        + (2.0/math.pi) * np.exp(-a) * E
        - 0.25 * (np.exp(-2.0*a) + np.exp(2.0*b))
        - math.log(2.0)
    )

    # Eq (18): σ^2 = ln H(λ,ζ;τ)
    # The PDF text layer expresses H via the angle-averaged moments; we implement it directly.
    # <f>:
    g = (2.0/math.pi) * np.exp(-a) * E

    # <f^2>:
    h = (3.0/2.0) * np.exp(-2.0*a) * E - 0.125 * (np.exp(-2.0*a) + np.exp(2.0*b))

    # H = <f^2>/<f>^2
    H = h / (g*g + 1e-300)

    return G, H

def topological_entropy_Eulerian(lam: np.ndarray, zeta: np.ndarray, tau: float) -> float:
    """
    Eq (19): S = (1/τ)[ (1/2)<G> + (1/4) ln <H> ]
    """
    G, H = G_and_H(lam, zeta, tau)
    meanG = float(np.mean(G))
    meanH = float(np.mean(H))
    return (1.0/tau) * (0.5*meanG + 0.25*math.log(meanH))

if __name__ == "__main__":
    # quick self-test with dummy samples
    rng = np.random.default_rng(0)
    lam = rng.normal(0.6, 0.1, size=5000)
    zeta = rng.normal(0.1, 0.1, size=5000)
    tau = 0.2
    print("S(Eulerian) =", topological_entropy_Eulerian(lam, zeta, tau))
    