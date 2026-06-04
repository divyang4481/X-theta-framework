import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def get_correlation_ellipsoid_radii(phi, mode="strength", cap=10.0):
    """
    Returns the radii for the correlation ellipsoid.

    mode="strength":
        rx = ry = |cos(2phi)|
        rz = 1

    mode="dual":
        rx = ry = 1/|cos(2phi)|
        rz = 1
        (capped by 'cap' parameter)
    """
    cos2phi = np.cos(2*phi)

    if mode == "strength":
        rx = np.abs(cos2phi)
        ry = rx
        rz = 1.0
    elif mode == "dual":
        if np.abs(cos2phi) < 1.0/cap:
            rx = cap
        else:
            rx = 1.0 / np.abs(cos2phi)
        ry = rx
        rz = 1.0
    else:
        raise ValueError(f"Unknown mode: {mode}")

    return rx, ry, rz

def plot_correlation_ellipsoid(phi, ax=None, mode="strength", cap=10.0):
    """
    Plots the correlation ellipsoid.

    mode="strength" -> direct correlation-strength ellipsoid (rx=ry=|cos(2phi)|, rz=1)
    mode="dual"     -> dual response ellipsoid (rx=ry=1/|cos(2phi)|, rz=1)
    """
    if ax is None:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')

    rx, ry, rz = get_correlation_ellipsoid_radii(phi, mode=mode, cap=cap)

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = rx * np.outer(np.cos(u), np.sin(v))
    y = ry * np.outer(np.sin(u), np.sin(v))
    z = rz * np.outer(np.ones(np.size(u)), np.cos(v))

    color = 'b' if mode == "strength" else 'r'
    ax.plot_surface(x, y, z, color=color, alpha=0.3)

    # Plot axes
    max_r = max(rx, ry, rz)
    ax.set_xlim([-max_r, max_r])
    ax.set_ylim([-max_r, max_r])
    ax.set_zlim([-max_r, max_r])

    ax.set_xlabel('X (Correlation)')
    ax.set_ylabel('Y (Correlation)')
    ax.set_zlabel('Z (Correlation)')

    title_suffix = "Strength" if mode == "strength" else "Dual Response"
    ax.set_title(f'Correlation Ellipsoid ({title_suffix}) [$\phi$ = {phi:.4f}]')

    return ax

def plot_correlation_strength_ellipsoid(phi, ax=None):
    """Wrapper for strength mode."""
    return plot_correlation_ellipsoid(phi, ax=ax, mode="strength")

def plot_dual_response_ellipsoid(phi, ax=None, cap=10.0):
    """Wrapper for dual mode."""
    return plot_correlation_ellipsoid(phi, ax=ax, mode="dual", cap=cap)

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
