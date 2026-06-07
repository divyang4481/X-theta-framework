"""
X-Theta relational connection and generator definitions.
"""
from __future__ import annotations
import numpy as np
from xtheta.quantum.correlation_tensor import get_g_rel

def relational_generator_minimal():
    """
    Returns the currently postulated minimal non-factorizable generator:
        G_rel = 1/2 (X⊗Y - Y⊗X)

    Scientific status:
        Kinematic ansatz / minimal generator.
        Not yet derived from a variational action.
    """
    return get_g_rel()
