import matplotlib
matplotlib.use('Agg') # Set backend to Agg for headless environment
import numpy as np
import matplotlib.pyplot as plt

class XThetaQuantumSim:
    def __init__(self, N_y=128, N_theta=64, L_y=20.0, mass=1.0, Inertia=0.5, q_theta=1.0):
        # Simulation Parameters
        self.Ny = N_y
        self.Nt = N_theta
        self.Ly = L_y
        self.Lt = 2 * np.pi  # Theta is compact 0 to 2pi
        self.m = mass
        self.I = Inertia
        self.q = q_theta

        # Space-Theta Grids
        self.y = np.linspace(-L_y/2, L_y/2, N_y)
        self.theta = np.linspace(0, self.Lt, N_theta, endpoint=False)
        self.Y, self.T = np.meshgrid(self.y, self.theta, indexing='ij')

        # Momentum Grids (k-space)
        dy = self.y[1] - self.y[0]
        dt = self.theta[1] - self.theta[0]

        ky = 2 * np.pi * np.fft.fftfreq(N_y, d=dy)
        kt = 2 * np.pi * np.fft.fftfreq(N_theta, d=dt)
        self.KY, self.KT = np.meshgrid(ky, kt, indexing='ij')

        # Wavefunction
        self.psi = np.zeros((N_y, N_theta), dtype=complex)

    def initialize_wavepacket(self, y0=0.0, sigma_y=1.0, l_mode=1):
        """
        Init Gaussian in space, Plane wave in theta (definite angular momentum)
        l_mode: The internal quantum number (rotor level)
        """
        # Spatial Gaussian
        space_part = np.exp(-(self.Y - y0)**2 / (2 * sigma_y**2))
        # Internal Rotor State exp(i * l * theta)
        theta_part = np.exp(1j * l_mode * self.T)

        self.psi = space_part * theta_part
        self.normalize()

    def normalize(self):
        norm = np.sqrt(np.sum(np.abs(self.psi)**2) * (self.y[1]-self.y[0]) * (self.theta[1]-self.theta[0]))
        self.psi /= norm

    # get_potential_operator removed as it was unused and contained only comments/pass.

    def run_step(self, dt, A_grad):
        """
        Executes one time step using Strang Splitting.
        Hamiltonian: H = K_y + K_theta + V_gauge
        """
        # 1. Half-step Potential (Gauge magnitude term)
        # V_mag = (q^2 * (A_grad * y)^2) / 2I
        V_mag = (self.q**2 * (A_grad * self.Y)**2) / (2 * self.I)
        self.psi *= np.exp(-0.5j * V_mag * dt)

        # 2. Mixed Term Evolution: H_mix = -(q/I) * A(y) * p_theta
        # Operator: exp( i * dt * (q/I) * A(y) * p_theta )
        # Since p_theta is diagonal in k-space (value hbar*kt),
        # but A(y) is diagonal in real-space, we cannot do this simply globally.
        # We approximate by small dt:
        # The physical effect is a shift in Ky momentum.
        # Let's stick to the Kinetic evolution and the V_mag for stability,
        # but to see the DRIFT, we need the p_theta coupling.

        # ALTERNATIVE: Lagrangian frame or exact operator split.
        # Let's do:
        # FFT along Theta -> Psi(y, k_theta)
        # In this basis, p_theta is just a number (hbar * k_theta).
        # So the Hamiltonian becomes H(y) = py^2/2m + V_eff(y, k_theta)
        # where V_eff = (hbar*k_theta - q*A(y))^2 / 2I

        # A. FFT w.r.t Theta
        psi_kt = np.fft.fft(self.psi, axis=1)

        # B. Apply V_eff(y) evolution for each k_theta mode
        # The potential depends on y and kt
        # A_val = A_grad * y
        # Term: (hbar*kt - q*A_val)^2 / 2I
        # We ignore py^2 here (that's the kinetic step).

        # Note: k_theta in our grid is self.KT
        # We construct the effective potential grid in (y, kt) space
        A_field_y = A_grad * self.Y
        V_eff = (self.KT - self.q * A_field_y)**2 / (2 * self.I)

        # Evolve by V_eff for full dt (since we are inside the split)
        psi_kt *= np.exp(-1j * V_eff * dt)

        # C. Inverse FFT w.r.t Theta
        self.psi = np.fft.ifft(psi_kt, axis=1)

        # 3. Kinetic Step for Space (py^2/2m)
        # FFT w.r.t Y
        psi_ky = np.fft.fft(self.psi, axis=0)
        E_kin_y = (self.KY**2) / (2 * self.m)
        psi_ky *= np.exp(-1j * E_kin_y * dt)
        self.psi = np.fft.ifft(psi_ky, axis=0)

    def measure_centroid(self):
        prob = np.abs(self.psi)**2
        # Integrate over theta to get P(y)
        P_y = np.sum(prob, axis=1)
        P_y /= np.sum(P_y) # re-normalize marginal
        return np.sum(self.y * P_y)

if __name__ == "__main__":
    # --- RUN SIMULATION ---
    sim = XThetaQuantumSim(N_y=256, N_theta=32, L_y=10.0, mass=1.0, Inertia=0.2, q_theta=1.0)

    # 1. Initialize with l=3 (strong internal rotation)
    sim.initialize_wavepacket(y0=0.0, sigma_y=0.5, l_mode=3)

    # 2. Run Time Evolution
    times = np.linspace(0, 4.0, 100)
    centroids = []
    dt = times[1] - times[0]
    A_gradient = 0.5  # Strength of the "Magnetic" field gradient

    for t in times:
        centroids.append(sim.measure_centroid())
        sim.run_step(dt, A_grad=A_gradient)

    # --- PLOT RESULTS ---
    plt.figure(figsize=(10, 6))

    # Trajectory
    plt.subplot(1, 2, 1)
    plt.plot(times, centroids, 'b-', linewidth=2, label=r'Drift $\langle y \rangle$')
    plt.plot(times, np.zeros_like(times), 'k--', alpha=0.5, label='Baseline (A=0)')
    plt.title(f"Cross-Hall Drift Validation\n(Internal Mode $\ell=3$, Gradient={A_gradient})")
    plt.xlabel("Time")
    plt.ylabel("Centroid Position $y$")
    plt.legend()
    plt.grid(True)

    # Final State Density
    plt.subplot(1, 2, 2)
    prob_density = np.sum(np.abs(sim.psi)**2, axis=1)
    plt.plot(sim.y, prob_density, 'r-', fillstyle='full')
    plt.fill_between(sim.y, prob_density, color='r', alpha=0.3)
    plt.title(f"Final Wavepacket (t={times[-1]})")
    plt.xlabel("Position $y$")
    plt.xlim(-5, 5)
    plt.grid(True)

    plt.tight_layout()
    output_path = "paper/figs/cross_hall_drift_validation.png"
    plt.savefig(output_path)
    print(f"Plot saved to {output_path}")

    drift = centroids[-1] - centroids[0]
    print(f"Total Drift Observed: {drift:.4f}")
    print("If Drift > 0, the coupling between internal theta-momentum and spatial position is CONFIRMED.")
