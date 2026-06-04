import numpy as np
from xtheta.montecarlo.uncertainty import monte_carlo_sensitivity
import pytest

def test_monte_carlo_no_noise():
    # With no noise, result should be deterministic (within precision)
    res = monte_carlo_sensitivity(5.972e24, 0.1, 6.371e6+500e3, 6.371e6, n_samples=10)
    assert res['phi_std'] < 1e-18
    assert res['s_max_std'] < 1e-15

def test_monte_carlo_with_noise():
    res = monte_carlo_sensitivity(5.972e24, 0.1, 6.371e6+500e3, 6.371e6,
                                theta_std=0.01, n_samples=100)
    assert res['phi_std'] > 0
    assert res['samples'] == 100
