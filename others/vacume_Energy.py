# X–θ vacuum energy: vanilla QFT vs GR, massless + massive KK on R^3 × S^1,
# stabilizer solve, joint (α, L) search, FRW relaxation with w(t), and PROGRESS/ETA.
#
# Requirements: numpy, mpmath
# Optional: set MPWORKPREC to raise precision for extreme scans

import os, sys, time
import numpy as np
import mpmath as mp

# ==============================
#    CONFIG / VERBOSITY
# ==============================
VERBOSE = True          # global on/off for progress prints
FRW_LOG_EVERY = 500     # print FRW status every N steps
FRW_CHI_MIN = 1e-13     # lower bound on L in meters (as χ=ln L reflecting wall)
FRW_CHI_MAX = 1e-6      # upper bound on L in meters (as χ=ln L reflecting wall)

# ---------- Precision ----------
mp.mp.dps = int(os.environ.get("MPWORKPREC", "50"))  # default 50 decimal digits

# ---------- Physical constants (SI) ----------
hbar = 1.054_571_817e-34       # J·s
c    = 2.997_924_58e8          # m/s
pi   = mp.pi
eV   = 1.602_176_634e-19       # J
GeV  = 1e9 * eV

# Observed vacuum energy (Planck-era ΛCDM)—our GR target
rho_obs = 6.0e-10  # J/m^3  ~ (2.24 meV)^4

# Reduced Planck mass (for FRW)
Mpl = 2.435e18 * eV  # J

# ---------- Safe scientific formatter for float/mpf ----------
def fe(x, digits=3):
    try:
        return f"{float(x):.{digits}e}"
    except Exception:
        s = mp.nstr(x, n=digits+6)
        try:
            return f"{float(s):.{digits}e}"
        except Exception:
            return s

# ---------- micro progress helper ----------
class Progress:
    def __init__(self, total, prefix=""):
        self.total = max(int(total), 1)
        self.prefix = prefix
        self.start = time.perf_counter()
        self.last_print = self.start

    def update(self, i, extra=""):
        if not VERBOSE:
            return
        now = time.perf_counter()
        # throttle to ~15 Hz max
        if now - self.last_print < 1/15:
            return
        self.last_print = now
        i = int(i)
        frac = min(max(i / self.total, 0.0), 1.0)
        elapsed = now - self.start
        eta = (elapsed / frac - elapsed) if frac > 0 else float("inf")
        bar_len = 24
        filled = int(frac * bar_len)
        bar = "█" * filled + "░" * (bar_len - filled)
        msg = (f"\r{self.prefix} [{bar}] {i}/{self.total} "
               f"({100*frac:5.1f}%) | elapsed {elapsed:6.1f}s | eta {eta:6.1f}s")
        if extra:
            msg += f" | {extra}"
        sys.stdout.write(msg)
        sys.stdout.flush()

    def done(self, extra=""):
        if not VERBOSE:
            return
        now = time.perf_counter()
        elapsed = now - self.start
        msg = f"\r{self.prefix} [████████████████████████] {self.total}/{self.total} (100.0%)"
        msg += f" | elapsed {elapsed:6.1f}s"
        if extra:
            msg += f" | {extra}"
        sys.stdout.write(msg + "\n")
        sys.stdout.flush()

# ---------- section timer ----------
class SectionTimer:
    def __init__(self, name):
        self.name = name
    def __enter__(self):
        self.t0 = time.perf_counter()
        if VERBOSE:
            print(f"\n>> {self.name} ...")
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        t1 = time.perf_counter()
        if VERBOSE:
            print(f"<< {self.name} done in {t1 - self.t0:0.3f}s")

# ---------- Vanilla QFT hard-cutoff zero-point (for context) ----------
def rho_zpe_cutoff(Emax_J, Nb=1, Nf=0):
    kmax = Emax_J / (hbar * c)
    return (hbar * c / (16 * pi**2)) * kmax**4 * (Nb - Nf)

# ---------- Massless KK on R^3 × S^1 via Bernoulli (exact) ----------
def B4(alpha):
    a = mp.fmod(alpha, 1.0)
    if a < 0: a += 1.0
    return a**4 - 2*a**3 + a**2 - mp.mpf(1)/mp.mpf(30)

def B4p(alpha):
    a = mp.fmod(alpha, 1.0)
    if a < 0: a += 1.0
    return 4*a**3 - 6*a**2 + 2*a

def V_KK_massless(alpha, L, sigma=+1):
    return -sigma * (hbar*c) * (pi**2/90) * B4(alpha) / (L**4)

def dVdalpha_massless(alpha, L, sigma=+1):
    return -sigma * (hbar*c) * (pi**2/90) * B4p(alpha) / (L**4)

def dVdL_massless(alpha, L, sigma=+1):
    Vkk = V_KK_massless(alpha, L, sigma=sigma)
    return -4*Vkk / L  # Vkk ∝ L^{-4}

# ---------- Massive KK on R^3 × S^1 (Bessel-K sums) ----------
def f_massive_term(n, x):
    if x == 0:
        return mp.mpf(1) / (n**4)
    z = 2*pi*n*x
    return (x**2 / (n**2)) * mp.besselk(2, z)

def f_massive_term_dx(n, x):
    if x == 0:
        return mp.mpf(0)
    n = mp.mpf(n)
    z = 2*pi*n*x
    K2 = mp.besselk(2, z)
    K1 = mp.besselk(1, z)
    K3 = mp.besselk(3, z)
    dK2_dz = -mp.mpf(1)/2 * (K1 + K3)
    return (2*x/(n**2)) * K2 + (x**2/(n**2)) * (2*pi*n) * dK2_dz

def V_KK_massive(alpha, L, m_J, sigma=+1, Nmax=2000, tol=mp.mpf('1e-30')):
    x = (m_J * L) / (hbar / c)  # = m c L / ħ
    pref = sigma * (hbar*c) / (pi**2 * L**4)
    total = mp.mpf('0')
    for n in range(1, Nmax+1):
        termf = f_massive_term(n, x)
        term  = termf * mp.cos(2*pi*n*alpha)
        total += term
        if abs(term) < tol and n > 10:
            break
    return pref * total

def dVdalpha_massive(alpha, L, m_J, sigma=+1, Nmax=2000, tol=mp.mpf('1e-30')):
    x = (m_J * L) / (hbar / c)
    pref = sigma * (hbar*c) / (pi**2 * L**4)
    total = mp.mpf('0')
    for n in range(1, Nmax+1):
        termf = f_massive_term(n, x)
        term  = (-2*pi*n) * termf * mp.sin(2*pi*n*alpha)
        total += term
        if abs(term) < tol and n > 10:
            break
    return pref * total

def dVdL_massive(alpha, L, m_J, sigma=+1, Nmax=2000, tol=mp.mpf('1e-30')):
    x = (m_J * L) / (hbar / c)
    pref = sigma * (hbar*c) / (pi**2)
    total = mp.mpf('0')
    for n in range(1, Nmax+1):
        f    = f_massive_term(n, x)
        dfdx = f_massive_term_dx(n, x)
        term = mp.cos(2*pi*n*alpha) * (-4*f / (L**5) + (dfdx * (x/L)) / (L**4))
        total += term
        if abs(term) < tol and n > 10:
            break
    return pref * total

# ---------- Field content ----------
class Field:
    def __init__(self, mass_eV, sigma=+1, bc_shift=0.0, name=""):
        self.mass_J = mass_eV * eV
        self.sigma = mp.mpf(sigma)
        self.shift = mp.mpf(bc_shift)  # 0 periodic (boson-like), 1/2 antiperiodic (fermion-like)
        self.name = name

def V_KK_total(alpha, L, fields, massless_use_closed=True, Nmax=2000):
    tot = mp.mpf('0')
    for f in fields:
        a = alpha + f.shift
        if f.mass_J == 0 and massless_use_closed:
            tot += V_KK_massless(a, L, sigma=f.sigma)
        else:
            tot += V_KK_massive(a, L, f.mass_J, sigma=f.sigma, Nmax=Nmax)
    return tot

def dVdalpha_total(alpha, L, fields, massless_use_closed=True, Nmax=2000):
    tot = mp.mpf('0')
    for f in fields:
        a = alpha + f.shift
        if f.mass_J == 0 and massless_use_closed:
            tot += dVdalpha_massless(a, L, sigma=f.sigma)
        else:
            tot += dVdalpha_massive(a, L, f.mass_J, sigma=f.sigma, Nmax=Nmax)
    return tot

def dVdL_total(alpha, L, fields, massless_use_closed=True, Nmax=2000):
    tot = mp.mpf('0')
    for f in fields:
        a = alpha + f.shift
        if f.mass_J == 0 and massless_use_closed:
            tot += dVdL_massless(a, L, sigma=f.sigma)
        else:
            tot += dVdL_massive(a, L, f.mass_J, sigma=f.sigma, Nmax=Nmax)
    return tot

# ---------- Two-term stabilizer: solve A,B at (α*, L*) ----------
def solve_stabilizer_AB(alpha_star, L_star, fields, p=2.0, q=6.0, rho_star=rho_obs):
    Vkk    = V_KK_total(alpha_star, L_star, fields)
    dVdLkk = dVdL_total(alpha_star, L_star, fields)
    a1, b1 = 1/(L_star**p), 1/(L_star**q)
    c1     = rho_star - Vkk
    a2, b2 = -p/(L_star**(p+1)), -q/(L_star**(q+1))
    c2     = -dVdLkk
    det = a1*b2 - a2*b1
    if abs(det) < mp.mpf('1e-80'):
        raise RuntimeError("Stabilizer equations nearly singular; choose different (p,q).")
    A = (c1*b2 - c2*b1)/det
    B = (a1*c2 - a2*c1)/det
    return A, B

def V_total(alpha, L, fields, A, B, p, q):
    return V_KK_total(alpha, L, fields) + A/(L**p) + B/(L**q)

def dVdalpha_total_with_stab(alpha, L, fields, A, B, p, q):
    return dVdalpha_total(alpha, L, fields)

def dVdL_total_with_stab(alpha, L, fields, A, B, p, q):
    return dVdL_total(alpha, L, fields) - (p*A)/(L**(p+1)) - (q*B)/(L**(q+1))

# ---------- Joint (α, L) search (coarse-to-fine) with progress ----------
def minimize_alpha_L(alpha_guess, L_guess, fields, p=2.0, q=6.0, rho_star=rho_obs):
    alpha = mp.mpf(alpha_guess)
    L     = mp.mpf(L_guess)
    A, B  = solve_stabilizer_AB(alpha, L, fields, p, q, rho_star)

    zooms = [mp.mpf('3.0'), mp.mpf('1.6'), mp.mpf('1.25'), mp.mpf('1.10')]
    for zi, span in enumerate(zooms, 1):
        with SectionTimer(f"zoom {zi}/{len(zooms)} span×{float(span):.2f}"):
            Ls = [L*(span**(mp.mpf(i)/60 - 0.5)) for i in range(121)]  # 121 L-points
            alphas = [mp.mpf(i)/1200 for i in range(1201)]             # 1201 α-points
            prog = Progress(total=len(Ls), prefix=f"  scanning L (zoom {zi})")
            best = (mp.mpf('Inf'), None, None)
            for i, Lc in enumerate(Ls, 1):
                vals = [abs(V_total(a, Lc, fields, A, B, p, q) - rho_star) for a in alphas]
                j = int(np.argmin([float(v) for v in vals]))
                s = vals[j]
                if s < best[0]:
                    best = (s, alphas[j], Lc)
                # live update with residual at current best
                prog.update(i, extra=f"best |V-ρ|={fe(best[0],3)} at α≈{fe(best[1],3)} L≈{fe(best[2],3)}")
            prog.done(extra=f"final best |V-ρ|={fe(best[0],3)}")
            _, alpha, L = best
            A, B = solve_stabilizer_AB(alpha, L, fields, p, q, rho_star)

    resid = V_total(alpha, L, fields, A, B, p, q) - rho_star
    return dict(alpha=float(alpha), L=float(L), A=float(A), B=float(B), resid=float(resid))

# =========================================================
#                    FRW (χ = ln L) — STABLE + PROGRESS
# =========================================================
f_alpha = 1.0 * Mpl
M_chi   = 1.0 * Mpl

def Hubble_chi(alpha, chi, adot, chidot, fields, A, B, p, q):
    L = mp.e**chi
    K = 0.5*(f_alpha**2)*(adot**2) + 0.5*(M_chi**2)*(chidot**2)
    V = V_total(alpha, L, fields, A, B, p, q)
    rho = K + V
    return mp.sqrt(max(rho, mp.mpf('1e-80'))/(3*Mpl**2))

def dV_dalpha_exact(alpha, chi, fields, A, B, p, q):
    L = mp.e**chi
    return dVdalpha_total_with_stab(alpha, L, fields, A, B, p, q)

def dV_dchi_exact(alpha, chi, fields, A, B, p, q):
    L = mp.e**chi
    dVdL = dVdL_total_with_stab(alpha, L, fields, A, B, p, q)
    return L * dVdL

def rhs_chi(alpha, chi, adot, chidot, fields, A, B, p, q):
    H = Hubble_chi(alpha, chi, adot, chidot, fields, A, B, p, q)
    a_dd   = -3*H*adot   - (1.0/f_alpha**2) * dV_dalpha_exact(alpha, chi, fields, A, B, p, q)
    chi_dd = -3*H*chidot - (1.0/M_chi**2)   * dV_dchi_exact(alpha, chi, fields, A, B, p, q)
    return a_dd, chi_dd

def rk4_step_chi(state, dt, fields, A, B, p, q,
                 chi_min=mp.log(FRW_CHI_MIN), chi_max=mp.log(FRW_CHI_MAX)):
    alpha, chi, adot, chidot = state
    def f(S):
        a, x, ad, xd = S
        a_dd, x_dd = rhs_chi(a, x, ad, xd, fields, A, B, p, q)
        return mp.matrix([ad, xd, a_dd, x_dd])

    S = mp.matrix(state)
    k1 = f(S)
    k2 = f(S + 0.5*dt*k1)
    k3 = f(S + 0.5*dt*k2)
    k4 = f(S + dt*k3)
    S_next = S + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

    # reflect χ softly at bounds, damp velocity on bounce
    if not mp.isfinite(S_next[1]) or S_next[1] < chi_min:
        S_next[1]  = chi_min
        S_next[3] *= -0.5
    elif S_next[1] > chi_max:
        S_next[1]  = chi_max
        S_next[3] *= -0.5

    # guard α (periodic physics)
    if not mp.isfinite(S_next[0]):
        S_next[0] = mp.mpf('0.5')
        S_next[2] *= -0.5

    return [float(S_next[i]) for i in range(4)]

# ---------- Demo / Harness ----------
if __name__ == "__main__":
    with SectionTimer("Vanilla QFT vs GR (sanity)"):
        for E in [100*GeV, 1e3*GeV, 1e6*GeV, 1.22e19*GeV]:
            rho = rho_zpe_cutoff(E, Nb=1, Nf=0)
            print(f"Emax={E/eV:8.2e} eV  ->  rho_zpe/rho_obs = {float(rho/rho_obs): .3e}")

    # Field content examples
    fields_massless = [Field(0.0, sigma=+1, bc_shift=0.0, name="scalar0")]
    fields_mixed = [
        Field(0.0,        sigma=+1,  bc_shift=0.0, name="boson0"),
        Field(100*GeV/eV, sigma=-1,  bc_shift=0.5, name="fermion100GeV"),
    ]

    with SectionTimer("Massless closed-form check at α=1/2"):
        for L in [1e-6, 1e-9, 1e-12]:
            Vmin = V_KK_massless(0.5, L, sigma=+1)
            print(f"L={L:.0e} m -> V_KK(α=1/2) = {fe(Vmin,3)} J/m^3 ; ratio={fe(Vmin/rho_obs,3)}")

    # Stabilizer (massless)
    with SectionTimer("Solve stabilizer at (α*=1/2, L*=1 nm) MASSLESS"):
        alpha_star, L_star = 0.5, 1e-9
        p, q = 2.0, 6.0
        A, B = solve_stabilizer_AB(alpha_star, L_star, fields_massless, p, q, rho_star=rho_obs)
        Vcheck   = V_total(alpha_star, L_star, fields_massless, A, B, p, q)
        dLcheck  = dVdL_total_with_stab(alpha_star, L_star, fields_massless, A, B, p, q)
        dAlphach = dVdalpha_total_with_stab(alpha_star, L_star, fields_massless, A, B, p, q)
        print(f"A={fe(A,3)}, B={fe(B,3)}, V*={fe(Vcheck,3)}, dV/dL|*={fe(dLcheck,3)}")
        print(f"sanity: ΔV=V*-ρ_obs = {fe(Vcheck - rho_obs,3)} ; dV/dα|* = {fe(dAlphach,3)}")

    with SectionTimer("Joint (α,L) search — MASSLESS"):
        res = minimize_alpha_L(0.49, 1.0e-9, fields_massless, p, q, rho_star=rho_obs)
        print(res)

    # Stabilizer (mixed)
    with SectionTimer("Solve stabilizer at (α*=1/2, L*=1 nm) MIXED"):
        A2, B2 = solve_stabilizer_AB(alpha_star, L_star, fields_mixed, p, q, rho_star=rho_obs)
        V2   = V_total(alpha_star, L_star, fields_mixed, A2, B2, p, q)
        dL2  = dVdL_total_with_stab(alpha_star, L_star, fields_mixed, A2, B2, p, q)
        dA2  = dVdalpha_total_with_stab(alpha_star, L_star, fields_mixed, A2, B2, p, q)
        print(f"A={fe(A2,3)}, B={fe(B2,3)}, V*={fe(V2,3)}, dV/dL|*={fe(dL2,3)}")
        print(f"sanity: ΔV=V*-ρ_obs = {fe(V2 - rho_obs,3)} ; dV/dα|* = {fe(dA2,3)}")

    with SectionTimer("Joint (α,L) search — MIXED"):
        res2 = minimize_alpha_L(0.48, 1.0e-9, fields_mixed, p, q, rho_star=rho_obs)
        print(res2)

    # FRW (χ = ln L) with periodic progress
    print("\n== FRW relaxation near the minimum (massless demo; χ = ln L) ==")
    alpha0 = 0.49
    chi0   = mp.log(0.97e-9)
    adot0, chidot0 = 0.0, 0.0
    state = [alpha0, chi0, adot0, chidot0]

    dt = 1e8         # seconds (~3.17 years)
    steps = 12000    # total integration steps
    prog = Progress(total=steps, prefix="  FRW")

    t_phys = 0.0
    for n in range(1, steps+1):
        state = rk4_step_chi(state, dt, fields_massless, A, B, p, q)
        t_phys += dt
        if (n % FRW_LOG_EVERY) == 0 or n == steps:
            alpha_t, chi_t, adot_t, chidot_t = state
            L_t = float(mp.e**chi_t)
            V_t = V_total(alpha_t, L_t, fields_massless, A, B, p, q)
            K_t = 0.5*(f_alpha**2)*(adot_t**2) + 0.5*(M_chi**2)*(chidot_t**2)
            w_t = (K_t - V_t)/(K_t + V_t)
            extra = f"t≈{fe(t_phys,3)} s, α≈{alpha_t:.4f}, L≈{fe(L_t,3)} m, w≈{fe(w_t,3)}"
            prog.update(n, extra=extra)
    prog.done()

    # final FRW summary print
    alpha_t, chi_t, adot_t, chidot_t = state
    L_t = float(mp.e**chi_t)
    V_t = V_total(alpha_t, L_t, fields_massless, A, B, p, q)
    K_t = 0.5*(f_alpha**2)*(adot_t**2) + 0.5*(M_chi**2)*(chidot_t**2)
    w_eff = (K_t - V_t)/(K_t + V_t)
    print(f"final: alpha≈{alpha_t:.6f}, L≈{L_t:.3e} m, V≈{fe(V_t,3)} J/m^3, w_eff≈{fe(w_eff,6)}")
