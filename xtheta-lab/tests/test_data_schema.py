import pandas as pd
import pytest
from xtheta.data.schema import validate_bell_schema, REQUIRED_COLUMNS

def test_validate_bell_schema_valid():
    df = pd.DataFrame({
        "trial_id": [1, 2],
        "timestamp": [0.1, 0.2],
        "alice_setting": [0, 1],
        "bob_setting": [0, 1],
        "alice_outcome": [1, -1],
        "bob_outcome": [1, 1],
        "source_file": ["test.csv", "test.csv"]
    })
    summary = validate_bell_schema(df)
    assert summary["row_count"] == 2
    assert len(summary["missing_required_columns"]) == 0
    assert summary["unique_alice_settings"] == [0, 1]

def test_validate_bell_schema_missing():
    df = pd.DataFrame({
        "alice_setting": [0],
        "bob_setting": [0]
    })
    summary = validate_bell_schema(df)
    assert "trial_id" in summary["missing_required_columns"]
    assert "alice_outcome" in summary["missing_required_columns"]
