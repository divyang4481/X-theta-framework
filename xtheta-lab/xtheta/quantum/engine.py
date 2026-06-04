from qutip import sigmax, sigmay, sigmaz, tensor, Qobj
import numpy as np

def get_paulis():
    """Returns the Pauli matrices X, Y, Z."""
    return sigmax(), sigmay(), sigmaz()

def get_bell_states():
    """Returns the Bell states |Psi-> and |Psi+>."""
    # |01> - |10> / sqrt(2)
    psi_minus = (tensor(Qobj([[1],[0]]), Qobj([[0],[1]])) -
                 tensor(Qobj([[0],[1]]), Qobj([[1],[0]]))).unit()

    # |01> + |10> / sqrt(2)
    psi_plus = (tensor(Qobj([[1],[0]]), Qobj([[0],[1]])) +
                tensor(Qobj([[0],[1]]), Qobj([[1],[0]]))).unit()

    return psi_minus, psi_plus

def get_relational_generator():
    """Returns the relational generator G_rel = 1/2 * (X⊗Y - Y⊗X)."""
    X, Y, Z = get_paulis()
    return 0.5 * (tensor(X, Y) - tensor(Y, X))

def get_unitary_evolution(phi):
    """Returns the unitary operator U_rel(phi) = exp(i * phi * G_rel)."""
    G_rel = get_relational_generator()
    return (1j * phi * G_rel).expm()

def evolve_state(phi):
    """
    Evolves the |Psi-> state by phi.
    Returns |Psi(phi)> = cos(phi)|Psi-> + sin(phi)|Psi+>.
    Also verifies it matches U_rel(phi) * |Psi->.
    """
    psi_minus, psi_plus = get_bell_states()
    # Analytic formula
    psi_phi_analytic = np.cos(phi) * psi_minus + np.sin(phi) * psi_plus

    # Numerical evolution
    U = get_unitary_evolution(phi)
    psi_phi_numeric = U * psi_minus

    return psi_phi_numeric

def compute_correlation_tensor(state):
    """
    Computes the correlation tensor T_ij = Tr(rho * (sigma_i ⊗ sigma_j)).
    Returns a 3x3 numpy array.
    """
    if state.type == 'ket':
        rho = state * state.dag()
    else:
        rho = state

    paulis = get_paulis()
    T = np.zeros((3, 3))

    for i in range(3):
        for j in range(3):
            op = tensor(paulis[i], paulis[j])
            T[i, j] = (rho * op).tr().real

    return T

def compute_chsh_max(T):
    """
    Computes S_max = 2 * sqrt(1 + cos^2(2*phi)) using the singular values of T.
    For X-Theta T = diag[-cos(2phi), -cos(2phi), -1].
    The two largest singular values squared are used.
    """
    # Singular values of T
    # T.T @ T = diag[cos^2(2phi), cos^2(2phi), 1]
    # S_max = 2 * sqrt(u1^2 + u2^2) where u1, u2 are the two largest singular values
    u = np.linalg.svd(T, compute_uv=False)
    u_sorted = np.sort(u)[::-1]
    s_max = 2 * np.sqrt(u_sorted[0]**2 + u_sorted[1]**2)
    return s_max

def compute_invariants(T):
    """
    Computes:
    I_Theta = Tr(T.T @ T)
    R_Theta = 3 - I_Theta
    """
    TT = T.T @ T
    I_theta = np.trace(TT)
    R_theta = 3 - I_theta
    return I_theta, R_theta
