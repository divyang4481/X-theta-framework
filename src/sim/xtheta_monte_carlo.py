import sys
import matplotlib

# Only force a headless backend when not running in Jupyter.
if 'ipykernel' not in sys.modules:
    matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftshift
from tqdm import tqdm

class XThetaMonteCarlo:
    def __init__(self, N_y=256, N_theta=64, L_y=40.0, mass=1.0, Inertia=0.5, q_theta=1.0):
        # --- Config Space ---
        self.Ny = N_y
        self.Nt = N_theta
        self.Ly = L_y
        self.Lt = 2 * np.pi
        self.m = mass
        self.I = Inertia
        self.q = q_theta

        # --- Grids ---
        self.y = np.linspace(-L_y/2, L_y/2, N_y)
        self.dy = self.y[1] - self.y[0]
        self.theta = np.linspace(0, self.Lt, N_theta, endpoint=False)

        # --- Momentum Space (k_y and k_theta) ---
        self.ky = 2 * np.pi * np.fft.fftfreq(N_y, d=self.dy)
        # Fix: The user provided formula was * N_theta which yields non-integers ~10.
        # Correct physics for S^1 domain of length 2pi requires integer modes.
        # fftfreq returns cycles/unit_length. Multiplied by 2pi gives angular wavenumber.
        self.kt = 2 * np.pi * np.fft.fftfreq(N_theta, d=self.Lt/N_theta) # Integer modes l

        # Meshgrids for vectorized ops
        self.Y, self.Theta = np.meshgrid(self.y, self.theta, indexing='ij')
        self.KY, self.KT = np.meshgrid(self.ky, self.kt, indexing='ij')

    def generate_thermal_state(self, temp=1.0, y0=0.0, sigma_y=1.0):
        """
        Creates a mixed state: Spatial Gaussian * Thermal Rotor Mixture
        Returns: A single realization (pure state) drawn from the thermal ensemble.
        """
        # 1. Spatial Wavefunction (Coherent Gaussian)
        psi_space = np.exp(-(self.y - y0)**2 / (4 * sigma_y**2))
        psi_space /= np.sqrt(np.sum(np.abs(psi_space)**2) * self.dy)

        # 2. Sample an internal angular momentum 'l' from Boltzmann dist
        energies = (hbar * self.kt)**2 / (2 * self.I) # Rotor energies
        probs = np.exp(-energies / temp)
        probs /= np.sum(probs)

        # Monte Carlo Step: Pick one 'l' mode based on probability
        l_selected = np.random.choice(self.kt, p=probs)

        # 3. Combine
        psi_theta = np.exp(1j * l_selected * self.theta)
        # Full 2D wavefunction
        psi = np.outer(psi_space, psi_theta)
        return psi, l_selected

    def evolve(self, psi_in, steps=200, dt=0.05, A_grad=0.2):
        """
        Split-Step Fourier Evolution on the Cylinder
        Hamiltonian: H = Py^2/2m + (P_theta - q*A(y))^2 / 2I
        Expanded:    H = Py^2/2m + P_theta^2/2I - (q/I)A(y)P_theta + q^2 A(y)^2/2I
        """
        psi = psi_in.copy()

        # Pre-compute Operators
        # 1. Kinetic Y (Diagonal in Ky)
        U_kin_y = np.exp(-1j * dt * (hbar*self.KY)**2 / (2*self.m))

        # 2. Kinetic Theta (Diagonal in Kt)
        U_kin_t = np.exp(-1j * dt * (hbar*self.KT)**2 / (2*self.I))

        # 3. Gauge Magnitude (Diagonal in Y)
        A_field = A_grad * self.Y
        V_gauge_mag = (self.q**2 * A_field**2) / (2*self.I)
        U_pot = np.exp(-1j * dt * V_gauge_mag)

        # 4. Mixed Interaction (Diagonal in Y AND Kt) -- The Drift Term
        # This term couples Y-position to Theta-momentum
        # Term: - (q / I) * A(y) * P_theta
        # Since we operate in (Y, Kt) basis during the split, we can apply this exactly.
        drift_phase = dt * (self.q / self.I) * A_field * (hbar * self.KT)
        U_drift = np.exp(1j * drift_phase)

        trajectory = []

        for _ in range(steps):
            # Record Center of Mass
            prob = np.sum(np.abs(psi)**2, axis=1) # Integrate out theta
            prob /= np.sum(prob)
            y_cm = np.sum(self.y * prob)
            trajectory.append(y_cm)

            # --- STEP 1: Y-Kinetic (FFT Y) ---
            psi = np.fft.fft(psi, axis=0) # To Ky
            psi *= U_kin_y
            psi = np.fft.ifft(psi, axis=0) # Back to Y

            # --- STEP 2: Theta-Kinetic (FFT Theta) ---
            psi = np.fft.fft(psi, axis=1) # To Kt
            psi *= U_kin_t

            # --- STEP 3: The "Theta Kick" (Mixed Term) ---
            # We are currently in (Y, Kt) basis.
            # Ideally suited for the mixed operator!
            psi *= U_drift

            psi = np.fft.ifft(psi, axis=1) # Back to Theta

            # --- STEP 4: Potential (Real Space) ---
            psi *= U_pot

        return np.array(trajectory)

if __name__ == "__main__":
    # --- RUN SIMULATION ---
    hbar = 1.0
    sim = XThetaMonteCarlo(N_y=256, N_theta=64, L_y=30.0, Inertia=0.5, q_theta=1.0)

    # Statistical Parameters
    N_trials = 50  # Number of Monte Carlo runs
    temp = 2.0     # Temperature (Controls spread of l-modes)
    A_gradient = 0.4 # Strength of the gauge gradient

    print(f"Running {N_trials} Monte Carlo trials at T={temp} with Gradient={A_gradient}...")

    drift_results = []
    l_modes = []

    plt.figure(figsize=(12, 5))

    # Monte Carlo Loop
    for i in tqdm(range(N_trials)):
        psi_init, l_val = sim.generate_thermal_state(temp=temp, y0=0.0)
        traj = sim.evolve(psi_init, steps=150, dt=0.05, A_grad=A_gradient)

        drift_results.append(traj[-1] - traj[0]) # Final displacement
        l_modes.append(l_val)

        # Plot individual faint lines
        color = 'red' if l_val > 0 else 'blue'
        if l_val == 0: color = 'gray'
        plt.subplot(1, 2, 1)
        plt.plot(traj, color=color, alpha=0.3, linewidth=1)

    # --- Analysis & Plotting ---
    plt.subplot(1, 2, 1)
    plt.title(f"Trajectory Ensemble (Red: $l>0$, Blue: $l<0$)")
    plt.xlabel("Time Steps")
    plt.ylabel("Center of Mass $y$")
    plt.grid(True, alpha=0.3)

    plt.subplot(1, 2, 2)
    # Correlation Plot: Drift vs L-mode
    plt.scatter(l_modes, drift_results, c=drift_results, cmap='coolwarm', edgecolors='k')
    plt.axhline(0, color='k', linestyle='--')
    plt.axvline(0, color='k', linestyle='--')
    plt.title("Statistical Validation: Drift vs Internal Momentum ($l$)")
    plt.xlabel("Internal Quantum Number $l$")
    plt.ylabel("Induced Spatial Drift $\Delta y$")
    plt.grid(True)

    plt.tight_layout()
    output_path = "paper/figs/monte_carlo_drift.png"
    plt.savefig(output_path)
    print(f"Plot saved to {output_path}")

    # --- STATISTICAL OUTPUT ---
    drifts = np.array(drift_results)
    ls = np.array(l_modes)
    correlation = np.corrcoef(ls, drifts)[0,1]

    print("\n--- STATISTICAL ANALYSIS ---")
    print(f"Correlation (L-mode vs Drift): {correlation:.4f}")
    print(f"Mean Drift (L > 0): {np.mean(drifts[ls > 0]):.4f}")
    print(f"Mean Drift (L < 0): {np.mean(drifts[ls < 0]):.4f}")
    print(f"Mean Drift (L = 0): {np.mean(drifts[ls == 0]):.4f} (Control)")
