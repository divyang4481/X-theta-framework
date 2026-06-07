"""
X-Theta correlation tensor and relational evolution.
"""
from __future__ import annotations
import numpy as np

def get_g_rel():
    """
    Returns the relational generator: G_rel = 1/2 (X⊗Y - Y⊗X).

    Scientific status:
    Kinematic minimal generator ansatz.
    Not yet derived from a variational action principle.
    """
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])

    # Kronecker products
    XY = np.kron(X, Y)
    YX = np.kron(Y, X)

    return 0.5 * (XY - YX)

def get_correlation_tensor(phi: float) -> np.ndarray:
    """
    Returns the X-Theta correlation tensor T(phi) = diag(-cos(2phi), -cos(2phi), -1).

    Scientific status:
    Mathematical theorem derived from unitary evolution of Bell singlet
    under U_rel(phi) = exp(i phi G_rel).
    """
    cos2phi = np.cos(2 * phi)
    return np.diag([-cos2phi, -cos2phi, -1.0])

def get_anisotropy_invariant(phi: float) -> float:
    """
    Returns the anisotropy invariant R_theta = 2 * sin^2(2phi).
    """
    return 2.0 * (np.sin(2 * phi)**2)
