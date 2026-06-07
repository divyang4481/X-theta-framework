from xtheta.geometry.schwarzschild import get_schwarzschild_radius, compute_Phi_rel
from xtheta.experiments.benchmarks import run_benchmark_scenarios
import numpy as np
import pytest

def test_schwarzschild_radius():
    # Sun: rs ~ 3km
    M_SUN = 1.989e30
    rs = get_schwarzschild_radius(M_SUN)
    assert 2950 < rs < 3000

def test_Phi_rel_scaling():
    # If r_emit == r_det, Phi should be 0
    Phi = compute_Phi_rel(1e30, 1.0, 1e6, 1e6)
    assert Phi == 0

    # If theta == 0, Phi should be 0
    Phi = compute_Phi_rel(1e30, 0.0, 1e5, 1e6)
    assert Phi == 0

def test_benchmarks():
    results = run_benchmark_scenarios()
    assert len(results) == 4

    micius = next(r for r in results if r['scenario'] == "Micius")
    # Micius Phi should be very small
    assert abs(micius['Phi_rel']) < 1e-10

    ns = next(r for r in results if r['scenario'] == "Neutron Star")
    # NS Phi should be significant
    assert abs(ns['Phi_rel']) > 0.1
