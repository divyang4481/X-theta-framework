import numpy as np
import pytest
from xtheta.data.bell_chsh import RunningAB, bootstrap_chsh

def test_running_ab_perfect_correlation():
    rab = RunningAB()
    # E(0,0)=1, E(0,1)=1, E(1,0)=1, E(1,1)=-1 => S=4 (Perfect Bell)
    for _ in range(10):
        rab.update(0, 0, 1)
        rab.update(0, 1, 1)
        rab.update(1, 0, 1)
        rab.update(1, 1, -1)

    assert rab.chsh() == 4.0
    assert rab.chsh_se() == 0.0

def test_running_ab_no_correlation():
    rab = RunningAB()
    # Random ±1 outcomes should give S ≈ 0
    np.random.seed(42)
    for _ in range(1000):
        rab.update(np.random.randint(2), np.random.randint(2), np.random.choice([1, -1]))

    assert abs(rab.chsh()) < 0.5

def test_bootstrap_chsh():
    n = 100
    alice_set = np.random.randint(2, size=n)
    bob_set = np.random.randint(2, size=n)
    alice_out = np.random.choice([1, -1], size=n)
    bob_out = np.random.choice([1, -1], size=n)

    boot = bootstrap_chsh(alice_out, bob_out, alice_set, bob_set, samples=100)
    assert "S_ci_low_95" in boot
    assert "S_ci_high_95" in boot
    assert boot["S_ci_low_95"] <= boot["S_ci_high_95"]
