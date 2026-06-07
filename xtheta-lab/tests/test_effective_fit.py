import numpy as np
import pytest
from xtheta.data.effective_fit import fit_Phi_eff_from_smax, anisotropy_from_Phi

def test_fit_Phi_eff_smax_limit():
    # S = 2sqrt(2) => Phi = 0
    res = fit_Phi_eff_from_smax(2.0 * np.sqrt(2))
    assert abs(res["Phi_eff"]) < 1e-6
    assert abs(res["R_theta_eff"]) < 1e-6

def test_fit_Phi_eff_bell_limit():
    # S = 2 => Phi = pi/4
    res = fit_Phi_eff_from_smax(2.0)
    assert abs(res["Phi_eff"] - np.pi/4) < 1e-6
    assert abs(res["R_theta_eff"] - 2.0) < 1e-6

def test_fit_Phi_eff_warning():
    res = fit_Phi_eff_from_smax(3.0)
    assert res["fit_status"] == "warning"
    assert "exceeds Tsirelson" in res["warning"]

def test_anisotropy_from_Phi():
    assert abs(anisotropy_from_Phi(0)) < 1e-6
    assert abs(anisotropy_from_Phi(np.pi/4) - 2.0) < 1e-6
