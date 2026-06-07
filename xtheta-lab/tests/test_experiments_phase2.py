from __future__ import annotations
import pytest
from xtheta.experiments.claim_classification import classify_claim, ClaimLevel
from xtheta.experiments.falsification_tests import check_falsification

def test_classify_claim():
    assert classify_claim(has_spacetime_metadata=False, is_synthetic=False) == ClaimLevel.PHENOMENOLOGICAL_FIT
    assert classify_claim(has_spacetime_metadata=False, is_synthetic=True) == ClaimLevel.SIMULATION
    assert classify_claim(has_spacetime_metadata=True, is_synthetic=False) == ClaimLevel.PHYSICAL_PREDICTION

def test_falsification_rule_1():
    # Zero observed anisotropy when non-zero predicted
    phi_pred = 0.05
    phi_eff = 0.001
    phi_eff_se = 0.01
    report = check_falsification(phi_pred, phi_eff, phi_eff_se)
    assert report['constrained'] is True
    assert any(r['id'] == 'RULE_1' for r in report['rules'])

def test_falsification_rule_2():
    # Inconsistency
    phi_pred = 0.10
    phi_eff = 0.20
    phi_eff_se = 0.01
    report = check_falsification(phi_pred, phi_eff, phi_eff_se)
    assert any(r['id'] == 'RULE_2' for r in report['rules'])
