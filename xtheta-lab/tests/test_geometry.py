from xtheta.geometry.schwarzschild import get_schwarzschild_radius, compute_phi_rel
from xtheta.experiments.benchmarks import run_benchmark_scenarios
import numpy as np
import pytest

def test_schwarzschild_radius():
    # Sun: rs ~ 3km
    M_SUN = 1.989e30
    rs = get_schwarzschild_radius(M_SUN)
    assert 2950 < rs < 3000

def test_phi_rel_scaling():
    # If r_emit == r_det, phi should be 0
    phi = compute_phi_rel(1e30, 1.0, 1e6, 1e6)
    assert phi == 0

    # If theta == 0, phi should be 0
    phi = compute_phi_rel(1e30, 0.0, 1e5, 1e6)
    assert phi == 0

def test_benchmarks():
    results = run_benchmark_scenarios()
    assert len(results) == 4

    micius = next(r for r in results if r['scenario'] == "Micius")
    # Micius phi should be very small
    assert abs(micius['phi_rel']) < 1e-10

    ns = next(r for r in results if r['scenario'] == "Neutron Star")
    # NS phi should be significant
    assert abs(ns['phi_rel']) > 0.1
