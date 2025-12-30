#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim

# ============================================================
# Setup
# ============================================================


def set_seed(seed=123):
    np.random.seed(seed)
    torch.manual_seed(seed)


def pick_device():
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


TS_LIMIT = 2.0 * math.sqrt(2.0)

# Canonical CHSH angles (radians)
# Alice: a0=0°, a1=45°
# Bob:   b0=22.5°, b1=-22.5°
A_ANGLES = torch.tensor([0.0, math.pi / 4], dtype=torch.float32)
B_ANGLES = torch.tensor([math.pi / 8, -math.pi / 8], dtype=torch.float32)

# ============================================================
# Model: continuous-angle Bell correlator (REALISM SAFE)
# ============================================================


class AngleBellModel(nn.Module):
    """
    Variants:
      - classical: c = mA*mB
      - theta_sig: mixes settings (signaling)
      - theta_hol: holonomy correlator, no-signaling marginals
    """

    def __init__(self, hidden_dim=64, variant="theta_hol", theta_max=1.2):
        super().__init__()
        self.variant = variant
        self.theta_max = float(theta_max)

        self.source = nn.Sequential(
            nn.Linear(16, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.Tanh(),
        )

        # Local marginals depend ONLY on (lambda, local angle)
        self.margA = nn.Sequential(
            nn.Linear(hidden_dim + 1, 64), nn.ReLU(), nn.Linear(64, 1), nn.Tanh()
        )
        self.margB = nn.Sequential(
            nn.Linear(hidden_dim + 1, 64), nn.ReLU(), nn.Linear(64, 1), nn.Tanh()
        )

        # Reparameterized theta: bounded
        self.raw_theta = nn.Parameter(torch.tensor(0.2))

        # Connection δ = tanh(net(λ, Δ)) ∈ [-1,1], scaled by bounded theta
        self.conn = nn.Sequential(
            nn.Linear(hidden_dim + 1, 64), nn.Tanh(), nn.Linear(64, 1), nn.Tanh()
        )

    @property
    def device(self):
        return next(self.parameters()).device

    def theta(self):
        # bounded theta in [-theta_max, +theta_max]
        return self.theta_max * torch.tanh(self.raw_theta)

    def forward_fixed(self, alpha, beta, noise=None, lam=None):
        """
        alpha,beta: (N,1)
        returns: c, mA, mB (each (N,1))
        """
        if lam is None:
            if noise is None:
                noise = torch.randn(alpha.shape[0], 16, device=alpha.device)
            lam = self.source(noise)

        mA = self.margA(torch.cat([lam, alpha], dim=1))
        mB = self.margB(torch.cat([lam, beta], dim=1))

        if self.variant == "classical":
            return (mA * mB), mA, mB

        if self.variant == "theta_sig":
            th = self.theta()
            alpha_eff = alpha + 0.5 * th * beta
            beta_eff = beta + 0.5 * th * alpha
            mA2 = self.margA(torch.cat([lam, alpha_eff], dim=1))
            mB2 = self.margB(torch.cat([lam, beta_eff], dim=1))
            return (mA2 * mB2), mA2, mB2

        if self.variant == "theta_hol":
            d = alpha - beta
            th = self.theta()
            delta = th * self.conn(torch.cat([lam, d], dim=1))  # bounded
            c = -torch.cos(2.0 * d + delta)
            # keep tiny classical mix for early stability
            c_eff = 0.98 * c + 0.02 * (mA * mB)
            return c_eff, mA, mB

        raise ValueError(f"Unknown variant: {self.variant}")


# ============================================================
# Helpers
# ============================================================


def sample_angles(batch_size, device):
    alpha = (torch.rand(batch_size, 1, device=device) - 0.5) * math.pi
    beta = (torch.rand(batch_size, 1, device=device) - 0.5) * math.pi
    return alpha, beta


def physicality_penalty(mA, mB, c):
    p_pp = (1 + mA + mB + c) / 4
    p_pm = (1 + mA - mB - c) / 4
    p_mp = (1 - mA + mB - c) / 4
    p_mm = (1 - mA - mB + c) / 4
    return (
        torch.relu(-p_pp).mean()
        + torch.relu(-p_pm).mean()
        + torch.relu(-p_mp).mean()
        + torch.relu(-p_mm).mean()
    )


def chsh_in_graph(model: AngleBellModel, n_per=4096):
    """
    Differentiable CHSH estimate: computes S and |S| as tensors.
    """
    pairs = [(0, 0, +1), (0, 1, +1), (1, 0, +1), (1, 1, -1)]
    S = 0.0
    for ai, bi, sgn in pairs:
        alpha = A_ANGLES[ai].to(model.device).view(1, 1).repeat(n_per, 1)
        beta = B_ANGLES[bi].to(model.device).view(1, 1).repeat(n_per, 1)
        c, _, _ = model.forward_fixed(alpha, beta)
        E = c.mean()
        S = S + float(sgn) * E
    return S, torch.abs(S)


@torch.no_grad()
def chsh_report(model: AngleBellModel, n_per=200000):
    pairs = [(0, 0, +1), (0, 1, +1), (1, 0, +1), (1, 1, -1)]
    S = 0.0
    for ai, bi, sgn in pairs:
        alpha = A_ANGLES[ai].to(model.device).view(1, 1).repeat(n_per, 1)
        beta = B_ANGLES[bi].to(model.device).view(1, 1).repeat(n_per, 1)
        c, _, _ = model.forward_fixed(alpha, beta)
        S += float(sgn) * c.mean().item()
    return S, abs(S)


# @torch.no_grad()
# def nosig_grid_audit(model: AngleBellModel, grid_n=9, n_per=40000, seed=0):
#     rng = np.random.default_rng(seed)
#     grid = np.linspace(-math.pi / 2, math.pi / 2, grid_n)
#     max_devA = 0.0
#     max_devB = 0.0

#     for alpha_val in grid:
#         betas = rng.choice(grid, size=3, replace=False)
#         vals = []
#         for b in betas:
#             alpha = torch.full((n_per, 1), float(alpha_val), device=model.device)
#             beta = torch.full((n_per, 1), float(b), device=model.device)
#             _, mA, _ = model.forward_fixed(alpha, beta)
#             vals.append(mA.mean().item())
#         max_devA = max(max_devA, max(vals) - min(vals))

#     for beta_val in grid:
#         alphas = rng.choice(grid, size=3, replace=False)
#         vals = []
#         for a in alphas:
#             alpha = torch.full((n_per, 1), float(a), device=model.device)
#             beta = torch.full((n_per, 1), float(beta_val), device=model.device)
#             _, _, mB = model.forward_fixed(alpha, beta)
#             vals.append(mB.mean().item())
#         max_devB = max(max_devB, max(vals) - min(vals))

#     return max_devA, max_devB


@torch.no_grad()
def nosig_full_audit(model: AngleBellModel, grid_n=9, n_per=40000, seed=0):
    rng = np.random.default_rng(seed)
    grid = np.linspace(-math.pi / 2, math.pi / 2, grid_n)
    max_A = 0.0
    max_B = 0.0

    for alpha_val in grid:
        betas = rng.choice(grid, size=3, replace=False)
        vals = []
        for b in betas:
            alpha = torch.full((n_per, 1), float(alpha_val), device=model.device)
            beta = torch.full((n_per, 1), float(b), device=model.device)
            _, mA, _ = model.forward_fixed(alpha, beta)
            pA_plus = ((1.0 + mA) / 2.0).mean().item()
            vals.append(pA_plus)
        max_A = max(max_A, max(vals) - min(vals))

    for beta_val in grid:
        alphas = rng.choice(grid, size=3, replace=False)
        vals = []
        for a in alphas:
            alpha = torch.full((n_per, 1), float(a), device=model.device)
            beta = torch.full((n_per, 1), float(beta_val), device=model.device)
            _, _, mB = model.forward_fixed(alpha, beta)
            pB_plus = ((1.0 + mB) / 2.0).mean().item()
            vals.append(pB_plus)
        max_B = max(max_B, max(vals) - min(vals))

    return max_A, max_B


# ============================================================
# Training
# ============================================================


def train_realism(
    model: AngleBellModel, epochs=6000, batch_size=2048, lr=0.003, log_every=500
):

    opt = optim.Adam(model.parameters(), lr=lr)

    # weights
    w_neutral = 30.0  # mean marginals ~0
    w_phys = 60.0  # non-negative probs
    w_nosig = 80.0  # mean drift test (marginals)
    w_S = 28.0  # CHSH drive (|S|)
    w_wall = 900.0  # Tsirelson wall (holonomy only)
    w_thetaL2 = 0.3  # gently discourages maxing theta

    # NOTE: sharpness term removed for holonomy realism (it can fight neutrality)
    # If you want it, add it only for classical, not for theta_hol.

    for ep in range(1, epochs + 1):
        opt.zero_grad()

        alpha, beta = sample_angles(batch_size, model.device)
        c, mA, mB = model.forward_fixed(alpha, beta)

        loss = 0.0

        # neutrality
        loss = loss + w_neutral * (mA.mean() ** 2 + mB.mean() ** 2)

        # no-signaling drift (mA shouldn't change when swapping beta, etc.)
        beta2 = (torch.rand(batch_size, 1, device=model.device) - 0.5) * math.pi
        _, mA2, _ = model.forward_fixed(alpha, beta2)
        loss = loss + w_nosig * ((mA - mA2).mean() ** 2)

        alpha2 = (torch.rand(batch_size, 1, device=model.device) - 0.5) * math.pi
        _, _, mB2 = model.forward_fixed(alpha2, beta)
        loss = loss + w_nosig * ((mB - mB2).mean() ** 2)

        # physicality
        loss = loss + w_phys * physicality_penalty(mA, mB, c)

        # maximize |S| (differentiable)
        S_t, Sabs_t = chsh_in_graph(model, n_per=4096)

        if model.variant == "theta_hol":
            excess = torch.relu(Sabs_t - TS_LIMIT)
            loss = (
                loss
                - w_S
                * torch.minimum(Sabs_t, torch.tensor(TS_LIMIT, device=model.device))
                + w_wall * excess**2
            )
        else:
            loss = loss - w_S * Sabs_t

        # keep theta from saturating (prevents phase-scramble collapse)
        th = model.theta()
        loss = loss + w_thetaL2 * (th**2)

        loss.backward()
        opt.step()

        if ep % log_every == 0:
            S_now, Sabs_now = chsh_report(model, n_per=120000)
            th_now = float(model.theta().detach().cpu().item())
            print(
                f"[ep {ep:5d}] loss={loss.item(): .4f}  |S|≈{Sabs_now: .4f}  S≈{S_now: .4f}  theta={th_now: .4f}"
            )


# ============================================================
# Main
# ============================================================


def main():
    set_seed(123)
    dev = pick_device()
    print("Using device:", dev)

    print("\nTraining CLASSICAL (angle-realism)...")
    classical = AngleBellModel(variant="classical", theta_max=1.2).to(dev)
    classical.raw_theta.data.zero_()
    train_realism(classical, epochs=10000, log_every=500)
    S, Sabs = chsh_report(classical, n_per=250000)

    dA, dB = nosig_full_audit(classical, grid_n=9, n_per=40000, seed=0)
    print(f"Classical CHSH: S={S:.4f}  |S|={Sabs:.4f}")
    print("\nNo-signaling audit (marginal probability drift):")
    print(f"  max Δ P(A=+1|α) across β choices: {dA:.4e}")
    print(f"  max Δ P(B=+1|β) across α choices: {dB:.4e}")

    print("\nTraining THETA_SIG (should signal)...")
    sig = AngleBellModel(variant="theta_sig", theta_max=1.2).to(dev)
    train_realism(sig, epochs=10000, log_every=500)
    S, Sabs = chsh_report(sig, n_per=250000)

    dA, dB = nosig_full_audit(sig, grid_n=9, n_per=40000, seed=1)
    print(f"Theta_SIG CHSH: S={S:.4f}  |S|={Sabs:.4f}")
    print("\nNo-signaling audit (marginal probability drift):")
    print(f"  max Δ P(A=+1|α) across β choices: {dA:.4e}")
    print(f"  max Δ P(B=+1|β) across α choices: {dB:.4e}")

    print("\nTraining THETA_HOL (Tsirelson wall)...")
    hol = AngleBellModel(variant="theta_hol", theta_max=1.2).to(dev)
    train_realism(hol, epochs=10000, log_every=500)
    S, Sabs = chsh_report(hol, n_per=300000)

    dA, dB = nosig_full_audit(hol, grid_n=9, n_per=50000, seed=2)
    print(f"Theta_HOL CHSH: S={S:.4f}  |S|={Sabs:.4f}")
    print("\nNo-signaling audit (marginal probability drift):")
    print(f"  max Δ P(A=+1|α) across β choices: {dA:.4e}")
    print(f"  max Δ P(B=+1|β) across α choices: {dB:.4e}")

    print("\nDone.")


if __name__ == "__main__":
    main()
