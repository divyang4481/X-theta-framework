import matplotlib
matplotlib.use('Agg') # Set backend to Agg for headless environment
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft
from tqdm import tqdm

# --- CORE PHYSICS ENGINE (Optimized SSFM from previous step) ---
class XThetaPhysicsEngine:
    def __init__(self, N_y=256, N_theta=64, L_y=40.0, mass=1.0, Inertia=0.5, q_theta=1.0):
        self.Ny, self.Nt = N_y, N_theta
        self.Ly = L_y
        self.m, self.I, self.q = mass, Inertia, q_theta

        self.y = np.linspace(-L_y/2, L_y/2, N_y)
        self.dy = self.y[1] - self.y[0]
        self.theta = np.linspace(0, 2*np.pi, N_theta, endpoint=False)

        # Momentum Grids
        self.ky = 2 * np.pi * np.fft.fftfreq(N_y, d=self.dy)
        # Fix: Correct wavenumber calculation for integer l modes
        self.kt = 2 * np.pi * np.fft.fftfreq(N_theta, d=2*np.pi/N_theta)
        self.Y, self.Theta = np.meshgrid(self.y, self.theta, indexing='ij')
        self.KY, self.KT = np.meshgrid(self.ky, self.kt, indexing='ij')

    def run_experiment(self, A_grad, Duration, temp=1.0):
        """Returns exact center-of-mass drift for a single particle realization."""
        # 1. Thermal State Generation
        hbar = 1.0
        energies = (hbar * self.kt)**2 / (2 * self.I)
        probs = np.exp(-energies / temp)
        probs /= np.sum(probs)
        l_mode = np.random.choice(self.kt, p=probs)

        # 2. SSFM Evolution (Simplified for speed)
        # Theoretical drift derived from Ehrenfest: Dy = 0.5 * (q/m) * (l/I) * A_grad * T^2
        # We use the analytical approximation here to generate 1000s of points quickly,
        # but add 'quantum jitter' based on the wavepacket width.

        drift_theoretical = 0.5 * (self.q / self.m) * (l_mode / self.I) * A_grad * (Duration**2)

        # Add Quantum Delocalization Noise (Wavepacket spreading)
        # Sigma(t) ~ Sigma_0 * sqrt(1 + t^2)
        spread = 0.5 * np.sqrt(1 + Duration**2)
        random_spread = np.random.normal(0, spread)

        return drift_theoretical + random_spread, l_mode

# --- VIRTUAL LAB EQUIPMENT (Adds Noise & Artifacts) ---
def virtual_interferometer(
    n_samples=500,
    grad_range=(0.0, 1.0),
    time_range=(1.0, 5.0),
    noise_level=0.2,
    efficiency=0.9
):
    engine = XThetaPhysicsEngine()
    data_records = []

    print(f"Generating {n_samples} experimental shots...")

    for _ in range(n_samples):
        # 1. Set Experimental Controls
        # We vary Gradient randomly, but Time in discrete steps (like a real scan)
        A_grad = np.random.uniform(*grad_range)
        duration = np.random.choice(np.linspace(*time_range, 5))

        # 2. Run Physics Engine
        true_drift, l_state = engine.run_experiment(A_grad, duration, temp=2.0)

        # 3. Apply Detector Physics
        # Efficiency: 10% chance the detector misses the particle entirely (NaN or outlier)
        if np.random.random() > efficiency:
            measured_y = np.nan # Dropped packet
        else:
            # Shot Noise + Readout Error
            # The detector has finite pixel size (simulated by rounding/noise)
            readout_noise = np.random.normal(0, noise_level)
            measured_y = true_drift + readout_noise

        # 4. Save Record
        record = {
            "run_id": np.random.randint(10000, 99999),
            "gradient_strength_T_m": round(A_grad, 4), # Tesla/meter equivalent
            "pulse_duration_s": round(duration, 2),    # Seconds
            "detected_y_mm": measured_y,               # Millimeters
            "temperature_k": 2.0 + np.random.normal(0, 0.05), # Slight environmental drift
            "system_status": "OK" if np.random.random() > 0.02 else "ERROR_SYNC" # Occasional glitches
        }
        data_records.append(record)

    return pd.DataFrame(data_records)

if __name__ == "__main__":
    # --- EXECUTE & SAVE ---
    df_synthetic = virtual_interferometer(n_samples=1000, noise_level=0.5)

    # Add some "Outliers" (e.g. cosmic rays hitting the detector)
    n_outliers = 20
    outlier_indices = np.random.choice(df_synthetic.index, n_outliers)
    df_synthetic.loc[outlier_indices, "detected_y_mm"] += np.random.uniform(10, 50, n_outliers) * np.random.choice([-1, 1], n_outliers)

    # Save to CSV
    filename = "paper/data/synthetic_cross_hall_exp.csv"
    df_synthetic.to_csv(filename, index=False)

    print(f"\nDataset generated: {filename}")
    print(df_synthetic.head())

    # --- QUICK PREVIEW PLOT ---
    plt.figure(figsize=(10, 6))
    clean_data = df_synthetic.dropna()
    clean_data = clean_data[clean_data["system_status"] == "OK"]

    plt.scatter(
        clean_data["gradient_strength_T_m"] * clean_data["pulse_duration_s"]**2,
        clean_data["detected_y_mm"],
        alpha=0.5, c='b', edgecolors='none', label='Noisy Data'
    )
    plt.title("Synthetic Experiment: Cross-Hall Drift\n(x-axis: $\\nabla A \\cdot T^2$, y-axis: Drift)")
    plt.xlabel("Control Parameter $\\xi = \\nabla A \\cdot T^2$")
    plt.ylabel("Measured Displacement $y$ (mm)")
    plt.grid(True, alpha=0.3)
    output_path = "paper/figs/synthetic_data_preview.png"
    plt.savefig(output_path)
    print(f"Plot saved to {output_path}")
