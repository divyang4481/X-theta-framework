import numpy as np
from xtheta.montecarlo.uncertainty import monte_carlo_sensitivity
import pytest

def test_monte_carlo_no_noise():
    # With no noise, result should be deterministic (within precision)
    res = monte_carlo_sensitivity(5.972e24, 0.1, 6.371e6+500e3, 6.371e6, n_samples=10)
    assert res['phi_std'] < 1e-18
    assert res['s_max_std'] < 1e-15
    assert 'concurrence_mean' in res
    assert 'r_theta_mean' in res

def test_monte_carlo_with_noise():
    res = monte_carlo_sensitivity(5.972e24, 0.1, 6.371e6+500e3, 6.371e6,
                                theta_std=0.01, n_samples=100, return_samples=True)
    assert res['phi_std'] > 0
    assert res['samples'] == 100
    assert 'phi_samples' in res
    assert len(res['phi_samples']) == 100
    assert 'txx_samples' in res

def test_monte_carlo_unphysical():
    # Should resample and still return n_samples
    # We pass a large std and a small mean to trigger unphysical values
    res = monte_carlo_sensitivity(1.0, 0.1, 1.0, 1.0,
                                mass_std=2.0, n_samples=10)
    assert res['samples'] == 10
