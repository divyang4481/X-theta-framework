import pandas as pd
import numpy as np
import os

def load_delft_local_only():
    # Force search in the script's own folder
    target_path = "bell_open_data.txt"
    
    if not os.path.exists(target_path):
        print(f"FAILED: Please move 'bell_open_data.txt' into: {os.getcwd()}")
        return None

    try:
        # Read raw lines to bypass stream locks
        with open(target_path, 'r', encoding='utf-8') as f:
            lines = [l.strip().split() for l in f.readlines() if l.strip()]
        
        # Convert to DataFrame
        df_raw = pd.DataFrame(lines).apply(pd.to_numeric, errors='coerce')
        
        # Drop headers/NaNs
        df_raw = df_raw.dropna().reset_index(drop=True)
        
        # Mapping for Hensen (Delft) 2015:
        # Index 1: Setting A, Index 2: Setting B, Index 3: Outcome A, Index 4: Outcome B
        df_mapped = pd.DataFrame({
            'a': df_raw[1].astype(int),
            'b': df_raw[2].astype(int),
            'prod': (df_raw[3] * df_raw[4]).astype(int)
        })
        
        print(f"--- SUCCESS ---")
        print(f"Loaded {len(df_mapped)} events from Delft 2015.")
        return df_mapped
        
    except Exception as e:
        print(f"Local Load Failure: {e}")
        return None

# Execute
df_real = load_delft_local_only()
if df_real is not None:
    print(df_real.head())