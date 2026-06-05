import os
import numpy as np
import pandas as pd
from xtheta.data.loaders import load_bell_data
from xtheta.data.bell_chsh import compute_correlations_per_setting, calculate_chsh_from_correlations
from xtheta.experiments.open_data_validation import fit_effective_phi_from_chsh

def test_data_pipeline_with_dummy():
    # Create dummy data
    dummy_path = "dummy_bell_data.npz"
    data = {
        'alice_setting': [0,0,1,1]*10,
        'bob_setting': [0,1,0,1]*10,
        'alice_outcome': [1,1,1,1]*10,
        'bob_outcome': [1,1,1,1]*10
    }
    np.savez(dummy_path, **data)

    df = load_bell_data(dummy_path)
    correlations = compute_correlations_per_setting(df)
    res = calculate_chsh_from_correlations(correlations, 0, 1, 0, 1)

    # E should be 1.0 for all settings
    # S = 1 + 1 + 1 - 1 = 2.0
    assert abs(res['S'] - 2.0) < 1e-10

    os.remove(dummy_path)

def test_phi_fitting():
    # S_max = 2*sqrt(2) approx 2.828 -> phi = 0
    res = fit_effective_phi_from_chsh(2.8284271247, 'smax-envelope')
    assert abs(res['phi_eff']) < 1e-5
