from __future__ import annotations
import pytest
import numpy as np
from xtheta.experiments.benchmarks import run_benchmark_scenarios

def test_run_benchmark_scenarios_schema():
    results = run_benchmark_scenarios()

    assert len(results) >= 4

    required_keys = [
        "scenario",
        "phi_rel",
        "S_max",
        "delta_S_from_tsirelson",
        "R_theta",
        "concurrence",
        "purity"
    ]

    for row in results:
        # Check presence
        for key in required_keys:
            assert key in row, f"Missing key '{key}' in benchmark result row: {row.keys()}"

        # Check types
        assert isinstance(row["scenario"], str)
        assert isinstance(row["phi_rel"], (int, float, np.number))
        assert isinstance(row["S_max"], (int, float, np.number))
        assert isinstance(row["delta_S_from_tsirelson"], (int, float, np.number))
        assert isinstance(row["R_theta"], (int, float, np.number))
        assert isinstance(row["concurrence"], (int, float, np.number))
        assert isinstance(row["purity"], (int, float, np.number))

        # Sanity checks
        assert row["S_max"] >= 0
        assert -0.0000001 <= row["purity"] <= 1.0000001
