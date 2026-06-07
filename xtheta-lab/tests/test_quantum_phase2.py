from __future__ import annotations
import pytest
import numpy as np
from xtheta.quantum.correlation_tensor import get_correlation_tensor, get_g_rel, get_anisotropy_invariant
from xtheta.quantum.chsh import chsh_s_from_tensor, s_max_horodecki
from xtheta.quantum.horodecki import get_horodecki_smax

def test_correlation_tensor_structure():
    T = get_correlation_tensor(0.0)
    assert np.allclose(T, np.diag([-1, -1, -1]))

    T_Phi = get_correlation_tensor(np.pi/4)
    # cos(2*pi/4) = cos(pi/2) = 0
    assert np.allclose(T_Phi, np.diag([0, 0, -1]))

def test_g_rel_hermitian():
    G = get_g_rel()
    # For G = 1/2(XY - YX):
    # (XY)^dagger = Y^dagger X^dagger = YX (since X, Y are hermitian)
    # (YX)^dagger = XY
    # So G^dagger = 1/2(YX - XY) = -G if it was anti-hermitian.
    # WAIT, (XY)^dagger = YX. So G^dagger = 1/2(YX - XY) = -G.
    # Let me re-check my manual calculation vs code.
    # Actually, in the code output G was [[0,0,0,0],[0,0,i,0],[0,-i,0,0],[0,0,0,0]]
    # G^dagger = [[0,0,0,0],[0,0,i,0],[0,-i,0,0],[0,0,0,0]]
    # It IS hermitian.
    assert np.allclose(G, G.conj().T)

def test_horodecki_smax():
    # Phi = 0 should give 2*sqrt(2)
    s_max = get_horodecki_smax(0.0)
    assert pytest.approx(s_max) == 2 * np.sqrt(2)

    # Phi = pi/4 should give 2
    s_max_45 = get_horodecki_smax(np.pi/4)
    assert pytest.approx(s_max_45) == 2.0

def test_chsh_from_tensor():
    T = get_correlation_tensor(0.0)
    # Standard settings for max violation
    A_vecs = [np.array([0,0,1]), np.array([1,0,0])]
    B_vecs = [np.array([1,0,1])/np.sqrt(2), np.array([-1,0,1])/np.sqrt(2)]

    S = chsh_s_from_tensor(T, A_vecs, B_vecs)
    assert pytest.approx(S) == 2 * np.sqrt(2)

def test_s_max_horodecki():
    T = get_correlation_tensor(0.1)
    s_max = s_max_horodecki(T)
    s_horo = get_horodecki_smax(0.1)
    assert pytest.approx(s_max) == s_horo
