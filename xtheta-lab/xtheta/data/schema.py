"""
Canonical Bell-event schema definitions and validation logic.
"""
import pandas as pd
from typing import Dict

REQUIRED_COLUMNS = [
    "trial_id",
    "timestamp",
    "alice_setting",
    "bob_setting",
    "alice_outcome",
    "bob_outcome",
    "source_file"
]

def validate_bell_schema(df: pd.DataFrame) -> dict:
    """
    Return validation summary for a Bell-event DataFrame:
    - missing required columns
    - row count
    - null counts
    - unique settings
    - outcome value summary
    """
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]

    summary = {
        "row_count": len(df),
        "missing_required_columns": missing,
        "null_counts": df.isnull().sum().to_dict(),
    }

    if "alice_setting" in df.columns:
        summary["unique_alice_settings"] = df["alice_setting"].unique().tolist()
    if "bob_setting" in df.columns:
        summary["unique_bob_settings"] = df["bob_setting"].unique().tolist()

    if "alice_outcome" in df.columns:
        summary["alice_outcome_counts"] = df["alice_outcome"].value_counts().to_dict()
    if "bob_outcome" in df.columns:
        summary["bob_outcome_counts"] = df["bob_outcome"].value_counts().to_dict()

    return summary
