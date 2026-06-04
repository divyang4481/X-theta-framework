import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_correlation_ellipsoid(phi, ax=None):
    """
    Plots the correlation ellipsoid cos^2(2phi)x^2 + cos^2(2phi)y^2 + z^2 = 1.
    """
    if ax is None:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')

    # Radii
    # a*x^2 + b*y^2 + c*z^2 = 1 => radii are 1/sqrt(a), 1/sqrt(b), 1/sqrt(c)
    # Here a = b = cos^2(2phi), c = 1
    # Radii: rx = ry = 1/|cos(2phi)|, rz = 1

    cos2phi = np.cos(2*phi)
    rx = 1.0 / np.abs(cos2phi) if np.abs(cos2phi) > 1e-10 else 10.0 # Cap for visualization
    ry = rx
    rz = 1.0

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = rx * np.outer(np.cos(u), np.sin(v))
    y = ry * np.outer(np.sin(u), np.sin(v))
    z = rz * np.outer(np.ones(np.size(u)), np.cos(v))

    ax.plot_surface(x, y, z, color='b', alpha=0.3)

    # Plot axes
    max_r = max(rx, ry, rz)
    ax.set_xlim([-max_r, max_r])
    ax.set_ylim([-max_r, max_r])
    ax.set_zlim([-max_r, max_r])

    ax.set_xlabel('X (Correlation)')
    ax.set_ylabel('Y (Correlation)')
    ax.set_zlabel('Z (Correlation)')
    ax.set_title(f'Correlation Ellipsoid ($\phi$ = {phi:.4f})')

    return ax

def plot_anisotropy_curve(phi_range):
    """Plots R_theta and CHSH max vs phi."""
    from xtheta.quantum.engine import compute_chsh_max, compute_correlation_tensor, evolve_state, compute_invariants

    rs = []
    chshs = []

    for phi in phi_range:
        psi = evolve_state(phi)
        T = compute_correlation_tensor(psi)
        _, r_theta = compute_invariants(T)
        s_max = compute_chsh_max(T)

        rs.append(r_theta)
        chshs.append(s_max)

    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(phi_range, rs, label='$R_{\Theta}$')
    plt.xlabel('$\phi$')
    plt.ylabel('$R_{\Theta}$')
    plt.title('Entanglement Anisotropy Invariant')
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(phi_range, chshs, label='$S_{max}$')
    plt.axhline(y=2*np.sqrt(2), color='r', linestyle='--', label='$2\sqrt{2}$')
    plt.xlabel('$\phi$')
    plt.ylabel('$S_{max}$')
    plt.title('Maximum Bell Violation')
    plt.legend()

    plt.tight_layout()
    plt.show()
