from __future__ import annotations
import pytest
from xtheta.data.bell_chsh import compute_chsh_variants

def test_compute_chsh_variants_basic():
    # Standard CHSH case
    E00, E01, E10, E11 = 0.5, 0.5, 0.5, -0.5
    # +++- : 0.5 + 0.5 + 0.5 - (-0.5) = 2.0
    # ++-+ : 0.5 + 0.5 - 0.5 + (-0.5) = 0.0
    # +-++ : 0.5 - 0.5 + 0.5 + (-0.5) = 0.0
    # -+++ : -0.5 + 0.5 + 0.5 + (-0.5) = 0.0
    res = compute_chsh_variants(E00, E01, E10, E11)

    assert res["+++-"] == 2.0
    assert res["++-+"] == 0.0
    assert res["+-++"] == 0.0
    assert res["-+++"] == 0.0
    assert res["max_abs"] == 2.0
    assert res["max_abs_convention"] == "+++-"

def test_compute_chsh_variants_negative_max():
    # Case where max absolute is negative variant
    E00, E01, E10, E11 = -0.7, -0.7, -0.7, 0.7
    res = compute_chsh_variants(E00, E01, E10, E11)

    assert res["+++-"] == -2.8
    assert res["max_abs"] == 2.8
    assert res["max_abs_convention"] == "+++-"

def test_compute_chsh_variants_mixed():
    E00, E01, E10, E11 = 0.1, -0.8, 0.2, 0.5
    # +++- : 0.1 - 0.8 + 0.2 - 0.5 = -1.0
    # ++-+ : 0.1 - 0.8 - 0.2 + 0.5 = -0.4
    # +-++ : 0.1 + 0.8 + 0.2 + 0.5 = 1.6
    # -+++ : -0.1 - 0.8 + 0.2 + 0.5 = -0.2

    res = compute_chsh_variants(E00, E01, E10, E11)
    assert pytest.approx(res["+-++"]) == 1.6
    assert pytest.approx(res["max_abs"]) == 1.6
    assert res["max_abs_convention"] == "+-++"
