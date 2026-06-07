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
    Phi_pred = 0.05
    Phi_eff = 0.001
    Phi_eff_se = 0.01
    report = check_falsification(Phi_pred, Phi_eff, Phi_eff_se)
    assert report['constrained'] is True
    assert any(r['id'] == 'RULE_1' for r in report['rules'])

def test_falsification_rule_2():
    # Inconsistency
    Phi_pred = 0.10
    Phi_eff = 0.20
    Phi_eff_se = 0.01
    report = check_falsification(Phi_pred, Phi_eff, Phi_eff_se)
    assert any(r['id'] == 'RULE_2' for r in report['rules'])
