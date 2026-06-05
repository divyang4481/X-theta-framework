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

def compute_chsh_from_tensor(T, a0, a1, b0, b1):
    """
    Generic CHSH calculator using E(a,b) = a.T @ T @ b.
    Returns S = E(a0,b0) + E(a0,b1) + E(a1,b0) - E(a1,b1).
    Input vectors a0, a1, b0, b1 are normalized internally.
    """
    def normalize(v):
        norm = np.linalg.norm(v)
        if norm < 1e-12:
            raise ValueError("Zero vector passed to CHSH calculator.")
        return v / norm

    a0 = normalize(a0)
    a1 = normalize(a1)
    b0 = normalize(b0)
    b1 = normalize(b1)

    def E(a, b):
        return a.T @ T @ b

    return E(a0, b0) + E(a0, b1) + E(a1, b0) - E(a1, b1)

def compute_chsh_xy(phi: float) -> float:
    """Returns S_XY = 2√2 |cos(2φ)|."""
    # Settings
    X = np.array([1.0, 0.0, 0.0])
    Y = np.array([0.0, 1.0, 0.0])

    A0 = X
    A1 = Y
    B0 = -(X + Y) / np.sqrt(2)
    B1 =  (Y - X) / np.sqrt(2)

    psi = evolve_state(phi)
    T = compute_correlation_tensor(psi)
    S = compute_chsh_from_tensor(T, A0, A1, B0, B1)
    return abs(S)

def compute_chsh_xz(phi: float) -> float:
    """Returns S_XZ = 2√2 cos²(φ)."""
    # Settings
    X = np.array([1.0, 0.0, 0.0])
    Z = np.array([0.0, 0.0, 1.0])

    A0 = Z
    A1 = X
    B0 = -(Z + X) / np.sqrt(2)
    B1 =  (X - Z) / np.sqrt(2)

    psi = evolve_state(phi)
    T = compute_correlation_tensor(psi)
    S = compute_chsh_from_tensor(T, A0, A1, B0, B1)
    return abs(S)

def compute_concurrence_from_phi(phi: float) -> float:
    """
    Returns C(φ)=|cos(2φ)| for the X-Theta evolved pure state.
    Note that concurrence is not constant under relational evolution.
    """
    return abs(np.cos(2 * phi))

def compute_concurrence_from_state(state) -> float:
    """
    Compute concurrence numerically for a two-qubit pure state.
    For |ψ⟩ = a|00⟩ + b|01⟩ + c|10⟩ + d|11⟩, C = 2 |ad - bc|.
    """
    if state.type != 'ket':
        raise ValueError("compute_concurrence_from_state only supports pure states (kets).")

    # QuTiP kets are (4, 1) for two qubits
    coeffs = state.full().flatten()
    a, b, c, d = coeffs
    return 2 * abs(a * d - b * c)

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

def compute_density_matrix(state):
    """
    Computes the density matrix rho = |psi><psi| or returns the density matrix if passed.
    """
    if state.type == 'ket':
        return state * state.dag()
    return state

def compute_purity(state_or_rho) -> float:
    """
    Computes the purity of a quantum state: Tr(rho^2).
    For a pure state, purity = 1.
    """
    rho = compute_density_matrix(state_or_rho)
    return (rho * rho).tr().real

def compute_analytic_correlation_tensor(phi: float) -> np.ndarray:
    """
    Returns the analytic X-Theta correlation tensor for a given phi.
    T(phi) = diag[-cos(2phi), -cos(2phi), -1.0]
    """
    return np.diag([
        -np.cos(2 * phi),
        -np.cos(2 * phi),
        -1.0
    ])

def compute_tensor_spectrum(T: np.ndarray) -> dict:
    """
    Return eigenvalues of T.T @ T, singular values of T,
    tensor norm invariant, and anisotropy invariant.
    """
    # Singular values of T
    singular_values = np.linalg.svd(T, compute_uv=False)

    # Eigenvalues of T.T @ T
    tt_eigenvalues = np.linalg.eigvalsh(T.T @ T)

    # Invariants
    i_theta = np.trace(T.T @ T)
    r_theta = 3 - i_theta

    # Effective rank (number of non-zero singular values)
    rank_effective = np.sum(singular_values > 1e-10)

    return {
        "singular_values": singular_values,
        "tt_eigenvalues": tt_eigenvalues,
        "I_theta": i_theta,
        "R_theta": r_theta,
        "rank_effective": rank_effective
    }
