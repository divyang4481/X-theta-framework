from xtheta.quantum.engine import (
    get_paulis, get_bell_states, get_relational_generator,
    get_unitary_evolution, evolve_state, compute_correlation_tensor,
    compute_chsh_max, compute_invariants, compute_chsh_from_tensor,
    compute_chsh_xy, compute_chsh_xz, compute_concurrence_from_phi,
    compute_concurrence_from_state
)
import numpy as np
import pytest

def test_paulis():
    X, Y, Z = get_paulis()
    assert X.shape == (2, 2)
    assert Y.shape == (2, 2)
    assert Z.shape == (2, 2)
    assert (X*X - 1).norm() < 1e-10

def test_bell_states():
    psi_minus, psi_plus = get_bell_states()
    assert (psi_minus.norm() - 1) < 1e-10
    assert (psi_plus.norm() - 1) < 1e-10
    assert abs(psi_minus.overlap(psi_plus)) < 1e-10

def test_relational_generator():
    G_rel = get_relational_generator()
    assert G_rel.isherm
    psi_minus, psi_plus = get_bell_states()
    res = G_rel * psi_minus
    assert (res + 1j * psi_plus).norm() < 1e-10

def test_unitary_evolution():
    phi = np.pi / 4
    U = get_unitary_evolution(phi)
    assert U.isunitary

    psi_minus, psi_plus = get_bell_states()
    psi_phi = U * psi_minus

    expected = np.cos(phi) * psi_minus + np.sin(phi) * psi_plus
    assert (psi_phi - expected).norm() < 1e-10

def test_correlation_tensor():
    phi = 0.1
    psi_phi = evolve_state(phi)
    T = compute_correlation_tensor(psi_phi)

    expected_T = np.diag([-np.cos(2*phi), -np.cos(2*phi), -1])
    np.testing.assert_allclose(T, expected_T, atol=1e-10)

def test_chsh_max():
    phi = 0.2
    psi_phi = evolve_state(phi)
    T = compute_correlation_tensor(psi_phi)
    s_max = compute_chsh_max(T)

    expected_s_max = 2 * np.sqrt(1 + np.cos(2*phi)**2)
    assert abs(s_max - expected_s_max) < 1e-10

def test_invariants():
    phi = 0.3
    T = compute_correlation_tensor(evolve_state(phi))
    I_theta, R_theta = compute_invariants(T)

    expected_I = 2 * np.cos(2*phi)**2 + 1
    expected_R = 2 * np.sin(2*phi)**2

    assert abs(I_theta - expected_I) < 1e-10
    assert abs(R_theta - expected_R) < 1e-10

def test_chsh_projections():
    for phi in [0.0, 0.1, 0.3, 0.4]:
        s_xy = compute_chsh_xy(phi)
        s_xz = compute_chsh_xz(phi)

        expected_xy = 2 * np.sqrt(2) * abs(np.cos(2*phi))
        expected_xz = 2 * np.sqrt(2) * (np.cos(phi)**2)

        assert abs(s_xy - expected_xy) < 1e-10
        assert abs(s_xz - expected_xz) < 1e-10

def test_chsh_generic():
    phi = 0.2
    psi = evolve_state(phi)
    T = compute_correlation_tensor(psi)

    X = np.array([1.0, 0.0, 0.0])
    Y = np.array([0.0, 1.0, 0.0])
    Z = np.array([0.0, 0.0, 1.0])

    # Test XY geometry
    A0 = X
    A1 = Y
    B0 = -(X + Y) / np.sqrt(2)
    B1 =  (Y - X) / np.sqrt(2)
    s_xy = compute_chsh_from_tensor(T, A0, A1, B0, B1)
    assert abs(abs(s_xy) - compute_chsh_xy(phi)) < 1e-10

    # Test normalization and error
    with pytest.raises(ValueError):
        compute_chsh_from_tensor(T, np.array([0,0,0]), A1, B0, B1)

def test_concurrence():
    for phi in [0.0, 0.1, 0.4, np.pi/4]:
        c_phi = compute_concurrence_from_phi(phi)
        psi = evolve_state(phi)
        c_state = compute_concurrence_from_state(psi)

        expected = abs(np.cos(2*phi))

        assert abs(c_phi - expected) < 1e-10
        assert abs(c_state - expected) < 1e-10

        # Verify S_max relation
        T = compute_correlation_tensor(psi)
        s_max = compute_chsh_max(T)
        assert abs(s_max - 2 * np.sqrt(1 + c_phi**2)) < 1e-10
