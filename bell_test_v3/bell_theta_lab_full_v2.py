#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim

# ============================================================
# Utilities
# ============================================================


def set_seed(seed: int = 123):
    np.random.seed(seed)
    torch.manual_seed(seed)


def device_name():
    return "cuda" if torch.cuda.is_available() else "cpu"


TS_LIMIT = 2.0 * math.sqrt(2.0)

# ============================================================
# Model
# ============================================================


class NeuralBellModel(nn.Module):
    def __init__(self, hidden_dim=32, variant="classical", learn_setting_offsets=True):
        super().__init__()
        self.variant = variant
        self.learn_setting_offsets = learn_setting_offsets

        self.source = nn.Sequential(
            nn.Linear(16, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.Tanh(),
        )
        self.alice = nn.Sequential(
            nn.Linear(hidden_dim + 1, 64), nn.ReLU(), nn.Linear(64, 1), nn.Tanh()
        )
        self.bob = nn.Sequential(
            nn.Linear(hidden_dim + 1, 64), nn.ReLU(), nn.Linear(64, 1), nn.Tanh()
        )

        self.theta_param = nn.Parameter(torch.tensor(0.5))

        if learn_setting_offsets:
            self.a_off = nn.Parameter(torch.zeros(2, 1))
            self.b_off = nn.Parameter(torch.zeros(2, 1))

        self._device_cached = None

    @property
    def device(self):
        if self._device_cached is None:
            self._device_cached = next(self.parameters()).device
        return self._device_cached

    def forward(self, batch_size: int):
        n = batch_size // 4
        a = torch.cat(
            [torch.full((n, 1), float(i), device=self.device) for i in [0, 0, 1, 1]],
            dim=0,
        )
        b = torch.cat(
            [torch.full((n, 1), float(i), device=self.device) for i in [0, 1, 0, 1]],
            dim=0,
        )
        c, mA, mB = self.forward_fixed(a, b, lam=None)
        return a, b, c, mA, mB

    def forward_fixed(self, a_t, b_t, lam=None, noise=None):
        if lam is None:
            if noise is None:
                noise = torch.randn(a_t.shape[0], 16, device=a_t.device)
            lam = self.source(noise)

        mA = self.alice(torch.cat([lam, a_t], dim=1))
        mB = self.bob(torch.cat([lam, b_t], dim=1))

        if self.variant == "classical":
            return (mA * mB), mA, mB

        if self.variant == "theta_v1":
            a_eff = a_t + (b_t * self.theta_param)
            b_eff = b_t + (a_t * self.theta_param)
            mA2 = self.alice(torch.cat([lam, a_eff], dim=1))
            mB2 = self.bob(torch.cat([lam, b_eff], dim=1))
            return (mA2 * mB2), mA2, mB2

        if self.variant == "theta_v2":
            if self.learn_setting_offsets:
                ai = a_t.long().squeeze(-1)
                bi = b_t.long().squeeze(-1)
                uA = torch.tanh(lam[:, 0:1] + self.a_off[ai])
                uB = torch.tanh(lam[:, 1:2] + self.b_off[bi])
            else:
                uA = torch.tanh(lam[:, 0:1] + 1.0 * (2 * a_t - 1))
                uB = torch.tanh(lam[:, 1:2] + 1.0 * (2 * b_t - 1))

            alpha = math.pi * uA
            beta = math.pi * uB

            g = torch.tanh(lam[:, 2:3])
            phi = (2 * a_t - 1) * (2 * b_t - 1)  # ±1
            delta = self.theta_param * phi * g

            c = torch.cos(alpha - beta + delta)
            c_eff = 0.95 * c + 0.05 * (mA * mB)
            return c_eff, mA, mB

        raise ValueError(f"Unknown variant: {self.variant}")


# ============================================================
# Training (Constitutional)
# ============================================================


def train_constrained(
    model: NeuralBellModel,
    targets=None,
    maximize_S=False,
    epochs=3000,
    lr=0.005,
    log_every=500,
):
    opt = optim.Adam(model.parameters(), lr=lr)

    w_target = 20.0
    w_S = 15.0
    w_marg = 40.0
    w_nosig = 120.0
    w_phys = 80.0

    for ep in range(1, epochs + 1):
        opt.zero_grad()
        a, b, c, mA, mB = model(1024)

        loss = 0.0

        # neutrality + no-signaling (mean-based)
        for val in [0, 1]:
            loss = loss + w_marg * (mA[a == val].mean() ** 2 + mB[b == val].mean() ** 2)
            loss = loss + w_nosig * (
                (mA[(a == val) & (b == 0)].mean() - mA[(a == val) & (b == 1)].mean())
                ** 2
            )
            loss = loss + w_nosig * (
                (mB[(b == val) & (a == 0)].mean() - mB[(b == val) & (a == 1)].mean())
                ** 2
            )

        # moment-positivity / physicality
        p_pp = (1 + mA + mB + c) / 4
        p_pm = (1 + mA - mB - c) / 4
        p_mp = (1 - mA + mB - c) / 4
        p_mm = (1 - mA - mB + c) / 4

        loss = loss + w_phys * (
            torch.relu(-p_pp).mean()
            + torch.relu(-p_pm).mean()
            + torch.relu(-p_mp).mean()
            + torch.relu(-p_mm).mean()
        )

        if targets is not None:
            for (A_s, B_s), t in targets.items():
                loss = loss + w_target * ((c[(a == A_s) & (b == B_s)].mean() - t) ** 2)

        if maximize_S:
            E00 = c[(a == 0) & (b == 0)].mean()
            E01 = c[(a == 0) & (b == 1)].mean()
            E10 = c[(a == 1) & (b == 0)].mean()
            E11 = c[(a == 1) & (b == 1)].mean()
            S = E00 + E01 + E10 - E11

            if model.variant == "theta_v2":
                excess = torch.relu(S - TS_LIMIT)
                loss = (
                    loss
                    - w_S
                    * torch.minimum(S, torch.tensor(TS_LIMIT, device=model.device))
                    + 500.0 * excess * excess
                )
            else:
                loss = loss - w_S * S

        loss.backward()
        opt.step()

        if ep % log_every == 0:
            with torch.no_grad():
                E00 = c[(a == 0) & (b == 0)].mean().item()
                E01 = c[(a == 0) & (b == 1)].mean().item()
                E10 = c[(a == 1) & (b == 0)].mean().item()
                E11 = c[(a == 1) & (b == 1)].mean().item()
                S_now = E00 + E01 + E10 - E11
                th = float(model.theta_param.detach().cpu().item())
            print(
                f"[ep {ep:5d}] loss={loss.item(): .4f}  S≈{S_now: .4f}  theta={th: .4f}"
            )


# ============================================================
# Sampling from moments -> outcomes
# ============================================================


def sample_joint_from_moments(rng: np.random.Generator, mA, mB, c):
    p_pp = (1 + mA + mB + c) / 4
    p_pm = (1 + mA - mB - c) / 4
    p_mp = (1 - mA + mB - c) / 4
    p_mm = (1 - mA - mB + c) / 4

    P = np.stack([p_pp, p_mm, p_pm, p_mp], axis=1)  # ++, --, +-, -+
    P = np.clip(P, 1e-12, 1.0)
    P = P / P.sum(axis=1, keepdims=True)

    outs = np.array([[1, 1], [-1, -1], [1, -1], [-1, 1]], dtype=int)
    idx = np.array([rng.choice(4, p=P[i]) for i in range(P.shape[0])], dtype=int)
    xy = outs[idx]
    return xy[:, 0], xy[:, 1]


# ============================================================
# Mixed-batch generator + stress transforms
# ============================================================


def make_mixed_settings(rng, n_trials):
    per = n_trials // 4
    a = np.concatenate(
        [np.zeros(per), np.zeros(per), np.ones(per), np.ones(per)]
    ).astype(np.int64)
    b = np.concatenate(
        [np.zeros(per), np.ones(per), np.zeros(per), np.ones(per)]
    ).astype(np.int64)
    idx = rng.permutation(len(a))
    return a[idx], b[idx]


def apply_stress(a, b, lam, rng, mode):
    """
    Stress modes:
      - None
      - permute_settings: permute (a,b) jointly (breaks ab structure vs lam)
      - permute_lambda: permute lam only
      - permute_both: permute settings and lam independently
      - swap_b_within_a: for each a separately, swap b labels among samples (keeps P(a), kills correlations that require b tied to lam within a)
      - swap_a_within_b: analogous
      - flip_phi_half: flip phi(a,b) on half samples while preserving marginals by swapping b within fixed a (surgical parity attack)
    """
    n = len(a)

    if mode is None:
        return a, b, lam

    if mode == "permute_settings":
        p = rng.permutation(n)
        return a[p], b[p], lam

    if mode == "permute_lambda":
        p = rng.permutation(n)
        return a, b, lam[p]

    if mode == "permute_both":
        p1 = rng.permutation(n)
        p2 = rng.permutation(n)
        return a[p1], b[p1], lam[p2]

    if mode == "swap_b_within_a":
        b2 = b.copy()
        for aval in [0, 1]:
            idx = np.where(a == aval)[0]
            b2[idx] = b2[rng.permutation(idx)]
        return a, b2, lam

    if mode == "swap_a_within_b":
        a2 = a.copy()
        for bval in [0, 1]:
            idx = np.where(b == bval)[0]
            a2[idx] = a2[rng.permutation(idx)]
        return a2, b, lam

    if mode == "flip_phi_half":
        # Flip parity on half the samples by swapping b within each a block
        # This preserves P(a), P(b|a) and keeps marginals stable, but scrambles phi-sign structure.
        b2 = b.copy()
        half = rng.choice(n, size=n // 2, replace=False)
        # In the selected half, flip b -> 1-b (changes phi sign for fixed a)
        b2[half] = 1 - b2[half]
        return a, b2, lam

    raise ValueError(f"Unknown stress mode: {mode}")


# ============================================================
# Audit (counts + strict no-signaling)
# ============================================================


def counts_audit_from_samples(a_np, b_np, x, y):
    counts = {
        (a, b): {(xx, yy): 0 for xx in [1, -1] for yy in [1, -1]}
        for a in [0, 1]
        for b in [0, 1]
    }
    for ai, bi, xi, yi in zip(a_np, b_np, x, y):
        counts[(int(ai), int(bi))][(int(xi), int(yi))] += 1
    return counts


def compute_E_and_marginals(counts):
    Eab, EA, EB, Ns, PAplus, PBplus = {}, {}, {}, {}, {}, {}
    for a in [0, 1]:
        for b in [0, 1]:
            d = counts[(a, b)]
            n = sum(d.values())
            Ns[(a, b)] = n

            eab = (d[(1, 1)] + d[(-1, -1)] - d[(1, -1)] - d[(-1, 1)]) / n
            ea = ((d[(1, 1)] + d[(1, -1)]) - (d[(-1, 1)] + d[(-1, -1)])) / n
            eb = ((d[(1, 1)] + d[(-1, 1)]) - (d[(1, -1)] + d[(-1, -1)])) / n

            Eab[(a, b)] = eab
            EA[(a, b)] = ea
            EB[(a, b)] = eb

            PAplus[(a, b)] = (d[(1, 1)] + d[(1, -1)]) / n
            PBplus[(a, b)] = (d[(1, 1)] + d[(-1, 1)]) / n

    return Eab, EA, EB, Ns, PAplus, PBplus


def sigma_diff_p(p1, n1, p2, n2):
    # Binomial SE for difference
    v1 = p1 * (1 - p1) / n1
    v2 = p2 * (1 - p2) / n2
    return abs(p1 - p2) / math.sqrt(v1 + v2 + 1e-18)


def audit_mixed(
    title, model: NeuralBellModel, n_trials=120_000, stress_mode=None, seed=1234
):
    rng = np.random.default_rng(seed)
    a_np, b_np = make_mixed_settings(rng, n_trials)

    a_t = torch.from_numpy(a_np.astype(np.float32)).to(model.device).view(-1, 1)
    b_t = torch.from_numpy(b_np.astype(np.float32)).to(model.device).view(-1, 1)

    with torch.no_grad():
        noise = torch.randn(len(a_np), 16, device=model.device)
        lam = model.source(noise)

        a2, b2, lam2 = apply_stress(a_np, b_np, lam, rng, stress_mode)
        a_t2 = torch.from_numpy(a2.astype(np.float32)).to(model.device).view(-1, 1)
        b_t2 = torch.from_numpy(b2.astype(np.float32)).to(model.device).view(-1, 1)

        c_t, mA_t, mB_t = model.forward_fixed(a_t2, b_t2, lam=lam2)

        mA = mA_t.squeeze(-1).cpu().numpy()
        mB = mB_t.squeeze(-1).cpu().numpy()
        c = c_t.squeeze(-1).cpu().numpy()

    x, y = sample_joint_from_moments(rng, mA, mB, c)
    counts = counts_audit_from_samples(a2, b2, x, y)
    Eab, EA, EB, Ns, PAplus, PBplus = compute_E_and_marginals(counts)

    # CHSH S
    S = Eab[(0, 0)] + Eab[(0, 1)] + Eab[(1, 0)] - Eab[(1, 1)]
    S_var = 0.0
    for a in [0, 1]:
        for b in [0, 1]:
            e = Eab[(a, b)]
            n = Ns[(a, b)]
            S_var += (1 - e * e) / n
    S_se = math.sqrt(S_var)

    # STRICT no-signaling on probabilities P(A=+|a,b) and P(B=+|a,b)
    max_z = 0.0
    for a in [0, 1]:
        zA = sigma_diff_p(PAplus[(a, 0)], Ns[(a, 0)], PAplus[(a, 1)], Ns[(a, 1)])
        max_z = max(max_z, zA)
    for b in [0, 1]:
        zB = sigma_diff_p(PBplus[(0, b)], Ns[(0, b)], PBplus[(1, b)], Ns[(1, b)])
        max_z = max(max_z, zB)

    hdr = title + (f"  [STRESS={stress_mode}]" if stress_mode else "")
    print("\n" + "=" * 68)
    print(f"AUDIT: {hdr}")
    print("-" * 68)
    for a in [0, 1]:
        for b in [0, 1]:
            print(
                f"({a},{b})  Eab={Eab[(a,b)]:+0.4f}  EA={EA[(a,b)]:+0.4f}  EB={EB[(a,b)]:+0.4f}   "
                f"P(A+)= {PAplus[(a,b)]:0.2f}  P(B+)= {PBplus[(a,b)]:0.2f}"
            )
    print("-" * 68)
    z_lhv = (abs(S) - 2.0) / S_se if S_se > 0 else float("inf")
    print(f"S = {S:0.4f} ± {S_se:0.4f}   |  (|S|-2)/SE = {z_lhv:+0.2f}σ")
    print(
        f"No-signaling (prob-based): max deviation = {max_z:0.2f}σ  -> {'FAIL' if max_z > 4.0 else 'PASS'}"
    )
    print("=" * 68)


# ============================================================
# Main
# ============================================================


def main():
    set_seed(123)
    dev = torch.device(device_name())
    print("Using device:", dev)

    # Classical: fit interior point targets
    print("\nTraining Classical (baseline to interior targets)...")
    classical = NeuralBellModel(variant="classical").to(dev)
    classical.theta_param.data.zero_()
    targets = {(0, 0): 0.65, (0, 1): 0.65, (1, 0): 0.65, (1, 1): 0.35}
    train_constrained(
        classical, targets=targets, maximize_S=False, epochs=4000, log_every=500
    )

    # Theta V1: signaling
    print("\nTraining Theta-V1 (signaling cheat, maximize S)...")
    theta_v1 = NeuralBellModel(variant="theta_v1").to(dev)
    train_constrained(theta_v1, maximize_S=True, epochs=6000, log_every=500)

    # Theta V2: holonomy
    print("\nTraining Theta-V2 (holonomy, maximize S with Tsirelson wall)...")
    theta_v2 = NeuralBellModel(variant="theta_v2").to(dev)
    train_constrained(theta_v2, maximize_S=True, epochs=6000, log_every=500)

    # Audits
    print("\n" + "=" * 68)
    print("MIXED AUDITS + STRESS BATTERY")
    print("=" * 68)

    for m, mdl in [
        ("Classical", classical),
        ("Theta-V1 (Sig)", theta_v1),
        ("Theta-V2 (Holonomy)", theta_v2),
    ]:
        audit_mixed(m, mdl, stress_mode=None, seed=100 + hash(m) % 1000)
        audit_mixed(m, mdl, stress_mode="permute_settings", seed=200 + hash(m) % 1000)
        audit_mixed(m, mdl, stress_mode="permute_lambda", seed=300 + hash(m) % 1000)
        audit_mixed(m, mdl, stress_mode="permute_both", seed=400 + hash(m) % 1000)
        audit_mixed(m, mdl, stress_mode="swap_b_within_a", seed=500 + hash(m) % 1000)
        audit_mixed(m, mdl, stress_mode="swap_a_within_b", seed=600 + hash(m) % 1000)
        audit_mixed(m, mdl, stress_mode="flip_phi_half", seed=700 + hash(m) % 1000)


if __name__ == "__main__":
    main()
