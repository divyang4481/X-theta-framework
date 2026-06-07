from xtheta.quantum.engine import (
    get_paulis, get_bell_states, get_relational_generator,
    get_unitary_evolution, evolve_state, compute_correlation_tensor,
    compute_chsh_max, compute_invariants, compute_chsh_from_tensor,
    compute_chsh_xy, compute_chsh_xz, compute_concurrence_from_Phi,
    compute_concurrence_from_state, compute_density_matrix, compute_purity,
    compute_analytic_correlation_tensor, compute_tensor_spectrum
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
    Phi = np.pi / 4
    U = get_unitary_evolution(Phi)
    assert U.isunitary

    psi_minus, psi_plus = get_bell_states()
    psi_Phi = U * psi_minus

    expected = np.cos(Phi) * psi_minus + np.sin(Phi) * psi_plus
    assert (psi_Phi - expected).norm() < 1e-10

def test_correlation_tensor():
    Phi = 0.1
    psi_Phi = evolve_state(Phi)
    T = compute_correlation_tensor(psi_Phi)

    expected_T = np.diag([-np.cos(2*Phi), -np.cos(2*Phi), -1])
    np.testing.assert_allclose(T, expected_T, atol=1e-10)

def test_chsh_max():
    Phi = 0.2
    psi_Phi = evolve_state(Phi)
    T = compute_correlation_tensor(psi_Phi)
    s_max = compute_chsh_max(T)

    expected_s_max = 2 * np.sqrt(1 + np.cos(2*Phi)**2)
    assert abs(s_max - expected_s_max) < 1e-10

def test_invariants():
    Phi = 0.3
    T = compute_correlation_tensor(evolve_state(Phi))
    I_theta, R_theta = compute_invariants(T)

    expected_I = 2 * np.cos(2*Phi)**2 + 1
    expected_R = 2 * np.sin(2*Phi)**2

    assert abs(I_theta - expected_I) < 1e-10
    assert abs(R_theta - expected_R) < 1e-10

def test_chsh_projections():
    for Phi in [0.0, 0.1, 0.3, 0.4]:
        s_xy = compute_chsh_xy(Phi)
        s_xz = compute_chsh_xz(Phi)

        expected_xy = 2 * np.sqrt(2) * abs(np.cos(2*Phi))
        expected_xz = 2 * np.sqrt(2) * (np.cos(Phi)**2)

        assert abs(s_xy - expected_xy) < 1e-10
        assert abs(s_xz - expected_xz) < 1e-10

def test_chsh_generic():
    Phi = 0.2
    psi = evolve_state(Phi)
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
    assert abs(abs(s_xy) - compute_chsh_xy(Phi)) < 1e-10

    # Test normalization and error
    with pytest.raises(ValueError):
        compute_chsh_from_tensor(T, np.array([0,0,0]), A1, B0, B1)

def test_concurrence():
    for Phi in [0.0, 0.1, 0.4, np.pi/4]:
        c_Phi = compute_concurrence_from_Phi(Phi)
        psi = evolve_state(Phi)
        c_state = compute_concurrence_from_state(psi)

        expected = abs(np.cos(2*Phi))

        assert abs(c_Phi - expected) < 1e-10
        assert abs(c_state - expected) < 1e-10

        # Verify S_max relation
        T = compute_correlation_tensor(psi)
        s_max = compute_chsh_max(T)
        assert abs(s_max - 2 * np.sqrt(1 + c_Phi**2)) < 1e-10

def test_density_matrix_trace():
    for Phi in [0.0, 0.1, 0.5]:
        psi = evolve_state(Phi)
        rho = compute_density_matrix(psi)
        assert abs(rho.tr() - 1.0) < 1e-10

def test_purity_preserved_under_unitary():
    for Phi in np.linspace(0, np.pi, 20):
        psi = evolve_state(Phi)
        purity = compute_purity(psi)
        # For a pure state evolved unitarily, purity should remain 1.0
        assert abs(purity - 1.0) < 1e-10

def test_analytic_tensor_vs_numerical():
    Phis = np.linspace(0, np.pi / 2, 100)
    for Phi in Phis:
        psi = evolve_state(Phi)
        T_num = compute_correlation_tensor(psi)
        T_ana = compute_analytic_correlation_tensor(Phi)
        np.testing.assert_allclose(T_num, T_ana, atol=1e-10)

def test_tensor_spectrum():
    Phis = [0.0, 0.1, 0.3, 0.4]
    for Phi in Phis:
        psi = evolve_state(Phi)
        T = compute_correlation_tensor(psi)
        spec = compute_tensor_spectrum(T)

        # singular values = [1, |cos(2Phi)|, |cos(2Phi)|] (unsorted in SVD)
        expected_sv = sorted([1.0, abs(np.cos(2*Phi)), abs(np.cos(2*Phi))], reverse=True)
        np.testing.assert_allclose(sorted(spec["singular_values"], reverse=True), expected_sv, atol=1e-10)

        # I_theta = 1 + 2cos^2(2Phi)
        expected_I = 1 + 2 * np.cos(2*Phi)**2
        assert abs(spec["I_theta"] - expected_I) < 1e-10

        # R_theta = 2sin^2(2Phi)
        expected_R = 2 * np.sin(2*Phi)**2
        assert abs(spec["R_theta"] - expected_R) < 1e-10
