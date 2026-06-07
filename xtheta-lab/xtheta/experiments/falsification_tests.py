"""
X-Theta falsification rules and scientific checks.
"""
from __future__ import annotations
import numpy as np

def check_falsification(phi_pred: float, phi_eff: float, phi_eff_se: float) -> dict:
    """
    Applies falsification rules to X-Theta results.
    """
    rules = []

    # Rule 1: Zero observed anisotropy when non-zero predicted
    if abs(phi_pred) > 1e-6 and abs(phi_eff) < 2 * phi_eff_se:
        rules.append({
            "id": "RULE_1",
            "status": "CONSTRAINED",
            "message": "Predicted phase is non-zero, but measured anisotropy is zero within sensitivity."
        })

    # Rule 2: Inconsistency between predicted and effective phase
    if abs(phi_pred - phi_eff) > 3 * phi_eff_se and phi_eff_se > 0:
        rules.append({
            "id": "RULE_2",
            "status": "INCONSISTENT",
            "message": "Significant deviation between predicted phase and effective fitted phase."
        })

    return {
        "falsified": any(r["status"] == "FALSIFIED" for r in rules),
        "constrained": any(r["status"] == "CONSTRAINED" for r in rules),
        "rules": rules
    }
