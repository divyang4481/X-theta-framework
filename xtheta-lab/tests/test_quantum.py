from xtheta.quantum.engine import get_paulis, get_bell_states, get_relational_generator, get_unitary_evolution, evolve_state, compute_correlation_tensor, compute_chsh_max, compute_invariants
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
    # G_rel should act on |Psi-> to give i|Psi+>
    # G_rel = 1/2(XY - YX)
    # G_rel |Psi-> = i |Psi+> ?
    # Let's check numerically
    psi_minus, psi_plus = get_bell_states()
    res = G_rel * psi_minus
    # G_rel |Psi-> = -i |Psi+>
    assert (res + 1j * psi_plus).norm() < 1e-10

def test_unitary_evolution():
    phi = np.pi / 4
    U = get_unitary_evolution(phi)
    assert U.isunitary

    psi_minus, psi_plus = get_bell_states()
    psi_phi = U * psi_minus

    # |Psi(phi)> = cos(phi)|Psi-> + i*sin(phi)|Psi+> ?
    # Wait, my engine says psi_phi_analytic = np.cos(phi) * psi_minus + np.sin(phi) * psi_plus
    # If G_rel |Psi-> = i|Psi+>, then exp(i phi G_rel) |Psi-> = (cos(phi) + i sin(phi) (i|Psi+><Psi-| + ...)) |Psi->
    # exp(i phi G_rel) |Psi-> = cos(phi)|Psi-> + i sin(phi) (G_rel/ (norm of G_rel effect)) ...
    # G_rel^2 |Psi-> = G_rel (i|Psi+>)
    # X|0> = |1>, X|1> = |0>
    # Y|0> = i|1>, Y|1> = -i|0>
    # |Psi-> = 1/sqrt(2) (|01> - |10>)
    # X⊗Y |01> = |1> ⊗ (-i|0>) = -i|10>
    # X⊗Y |10> = |0> ⊗ (i|1>) = i|01>
    # X⊗Y |Psi-> = 1/sqrt(2) (-i|10> - i|01>) = -i |Psi+>
    # Y⊗X |01> = (i|1>) ⊗ |0> = i|10>
    # Y⊗X |10> = (-i|0>) ⊗ |1> = -i|01>
    # Y⊗X |Psi-> = 1/sqrt(2) (i|10> + i|01>) = i |Psi+>
    # G_rel |Psi-> = 1/2 (-i|Psi+> - i|Psi+>) = -i |Psi+>
    # So U_rel |Psi-> = exp(i phi G_rel) |Psi-> = cos(phi)|Psi-> + i sin(phi) (-i|Psi+>) = cos(phi)|Psi-> + sin(phi)|Psi+>
    # Correct!

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
