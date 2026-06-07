import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def get_correlation_ellipsoid_radii(Phi, mode="strength", cap=10.0):
    """
    Returns the radii for the correlation ellipsoid.

    mode="strength":
        rx = ry = |cos(2Phi)|
        rz = 1

    mode="dual":
        rx = ry = 1/|cos(2Phi)|
        rz = 1
        (capped by 'cap' parameter)
    """
    cos2Phi = np.cos(2 * Phi)

    if mode == "strength":
        rx = np.abs(cos2Phi)
        ry = rx
        rz = 1.0
    elif mode == "dual":
        if np.abs(cos2Phi) < 1.0/cap:
            rx = cap
        else:
            rx = 1.0 / np.abs(cos2Phi)
        ry = rx
        rz = 1.0
    else:
        raise ValueError(f"Unknown mode: {mode}")

    return rx, ry, rz

def plot_correlation_ellipsoid(Phi, ax=None, mode="strength", cap=10.0):
    """
    Plots the correlation ellipsoid.

    mode="strength" -> direct correlation-strength ellipsoid (rx=ry=|cos(2Phi)|, rz=1)
    mode="dual"     -> dual response ellipsoid (rx=ry=1/|cos(2Phi)|, rz=1)
    """
    if ax is None:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')

    rx, ry, rz = get_correlation_ellipsoid_radii(Phi, mode=mode, cap=cap)

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
    ax.set_title(rf'Correlation Ellipsoid ({title_suffix}) [$\Phi$ = {Phi:.4f}]')

    return ax

def plot_correlation_strength_ellipsoid(Phi, ax=None):
    """Wrapper for strength mode."""
    return plot_correlation_ellipsoid(Phi, ax=ax, mode="strength")

def plot_dual_response_ellipsoid(Phi, ax=None, cap=10.0):
    """Wrapper for dual mode."""
    return plot_correlation_ellipsoid(Phi, ax=ax, mode="dual", cap=cap)

def plot_anisotropy_curve(Phi_range):
    """Plots R_theta and CHSH max vs Phi."""
    from xtheta.quantum.engine import compute_chsh_max, compute_correlation_tensor, evolve_state, compute_invariants

    rs = []
    chshs = []

    for Phi in Phi_range:
        psi = evolve_state(Phi)
        T = compute_correlation_tensor(psi)
        _, r_theta = compute_invariants(T)
        s_max = compute_chsh_max(T)

        rs.append(r_theta)
        chshs.append(s_max)

    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(Phi_range, rs, label=r'$R_{\Theta}$')
    plt.xlabel(r'$\Phi$')
    plt.ylabel(r'$R_{\Theta}$')
    plt.title('Entanglement Anisotropy Invariant')
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(Phi_range, chshs, label=r'$S_{max}$')
    plt.axhline(y=2*np.sqrt(2), color='r', linestyle='--', label=r'$2\sqrt{2}$')
    plt.xlabel(r'$\Phi$')
    plt.ylabel(r'$S_{max}$')
    plt.title('Maximum Bell Violation')
    plt.legend()

    plt.tight_layout()
    plt.show()

def plot_random_chsh_landscape(df, output_path=None):
    """
    Plots the random CHSH landscape.
    df: DataFrame from simulate_random_chsh_landscape
    """
    import matplotlib.pyplot as plt
    import numpy as np

    plt.figure(figsize=(10, 6))

    # Plot random CHSH cloud
    plt.scatter(df['Phi'], df['S_abs'], alpha=0.1, s=1, color='gray', label='Random CHSH')

    # Plot S_max envelope (it's the same for all samples at a given Phi)
    Phi_unique = df['Phi'].unique()
    s_max_unique = df.groupby('Phi')['S_max'].first()
    plt.plot(Phi_unique, s_max_unique, 'r-', linewidth=2, label=r'Horodecki $S_{max}$')

    # Limits
    plt.axhline(y=2.0, color='blue', linestyle='--', label='Bell Limit (2.0)')
    plt.axhline(y=2*np.sqrt(2), color='green', linestyle=':', label=r'Tsirelson Limit ($2\sqrt{2}$)')

    plt.xlabel(r'Relational Phase $\Phi$')
    plt.ylabel('CHSH $|S|$')
    plt.title('Random CHSH Landscape and Horodecki Envelope')
    plt.legend(loc='upper right')
    plt.grid(True, alpha=0.3)

    if output_path:
        plt.savefig(output_path)
        print(f"Plot saved to {output_path}")

    plt.show()
