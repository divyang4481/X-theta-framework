import numpy as np

def simulate_interference(x, slit_sep, wavelength, L, theta_amp=0.0, theta_freq=0.0, t=0.0):
    """
    Minimal scalar-field toy model: two coherent sources + internal phase modulation.
    I(x) ~ | e^{i(k r1 + phi_theta)} + e^{i(k r2 - phi_theta)} |^2
    """
    k = 2*np.pi / wavelength
    r1 = np.sqrt(L**2 + (x - slit_sep/2)**2)
    r2 = np.sqrt(L**2 + (x + slit_sep/2)**2)

    # internal theta-phase (drive-locked)
    phi_theta = theta_amp * np.sin(2*np.pi*theta_freq * t)

    psi = np.exp(1j*(k*r1 + phi_theta)) + np.exp(1j*(k*r2 - phi_theta))
    I = (psi * np.conj(psi)).real
    return I / I.max()
