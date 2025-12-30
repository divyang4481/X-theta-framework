#!/usr/bin/env python3
"""
x_theta_simulation.py
=========================================================
Full simulation + analysis pipeline for X–Θ holonomy tests.

What you get:
1) Synthetic Bell-CHSH trials with controllable Θ modulation
2) Two-loop protocol (CW vs CCW) that shares endpoints but differs by path
3) Adversarial “systematics generator” (timing drift, basis-dependent efficiency,
   coincidence-window artifact, memory in RNG, EOM-induced loss, etc.)
4) GPU-optional (PyTorch). Works on CPU too.
5) Analysis:
   - E_ab(θ), S(θ) with error bars
   - ΔS(θ) = S_CW(θ) - S_CCW(θ)
   - Block bootstrap SE
   - Permutation test for ΔS
   - FFT spectrum of E_ab(θ) with θ-shuffle / phase-scramble surrogates

Run examples:
  python x_theta_simulation.py --N 2000000 --theta-bins 64 --use-gpu 1
  python x_theta_simulation.py --N 1000000 --systematics hard --use-gpu 1


Suggested runs (do these in order)

Clean holonomy present
python x_theta_simulation.py --N 2000000 --theta-bins 64 --systematics none --holonomy-amp 0.2 --use-gpu 1

Holonomy off (should vanish)
python x_theta_simulation.py --N 2000000 --theta-bins 64 --systematics none --holonomy-amp 0.0 --use-gpu 1

Hard systematics, holonomy off (pipeline should not hallucinate)
python x_theta_simulation.py --N 2000000 --theta-bins 64 --systematics hard --holonomy-amp 0.0 --use-gpu 1

Hard systematics + holonomy on (can you still detect it?)
python x_theta_simulation.py --N 2000000 --theta-bins 64 --systematics hard --holonomy-amp 0.2 --use-gpu 1

It writes an .npz in out/ containing everything you need for plots.
"""

import argparse
import math
import os
from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np
import torch

# -----------------------------
# Utilities
# -----------------------------


def pick_device(use_gpu: int) -> torch.device:
    if use_gpu and torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def set_seed(seed: int):
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def wrap_pi(x: torch.Tensor) -> torch.Tensor:
    # wrap into (-pi, pi]
    return (x + math.pi) % (2 * math.pi) - math.pi


def sigmoid(x: torch.Tensor) -> torch.Tensor:
    return 1 / (1 + torch.exp(-x))


# -----------------------------
# Configuration
# -----------------------------


@dataclass
class SimConfig:
    N: int = 1_000_000
    theta_bins: int = 64
    block_size: int = 20_000  # for block bootstrap + permutation
    bootstrap_rounds: int = 300
    perm_rounds: int = 400
    seed: int = 123
    use_gpu: int = 1

    # Physics-ish parameters
    # Base correlation strength; ideal singlet gives S up to 2√2 with right angles.
    vis: float = 0.98

    # CHSH measurement angles (radians) for maximal QM violation
    # Alice: a0=0, a1=pi/4; Bob: b0=pi/8, b1=-pi/8
    a0: float = 0.0
    a1: float = math.pi / 4
    b0: float = math.pi / 8
    b1: float = -math.pi / 8

    # Θ modulation
    theta_amp: float = 0.7  # amplitude of Θ contribution to effective phase
    theta_freq: int = 1  # harmonic in θ
    holonomy_amp: float = 0.20  # CW/CCW holonomy offset amplitude
    holonomy_harm: int = 1  # harmonic in θ for holonomy

    # Systematics toggle
    systematics: str = "none"  # "none" | "mild" | "hard"

    # Coincidence window artifact model
    # Think of it as: detection time jitter depends on setting + Θ,
    # and coincidence selection biases correlations.
    coincidence_window: float = 3.0e-9  # seconds (toy)
    time_jitter_sigma: float = 0.8e-9  # seconds (toy)

    # Memory in settings RNG
    rng_memory: float = 0.0  # 0..1 ; higher => settings correlate with previous

    # EOM-induced loss correlated with Θ and maybe with basis
    eom_loss_amp: float = 0.0

    # Basis-dependent detection efficiency (loophole-ish)
    basis_eff_amp: float = 0.0

    # Slow drift in timing offset (can fake effects)
    timing_drift_amp: float = 0.0

    # Add local “electronics bias” in outcomes (dangerous)
    local_bias_amp: float = 0.0


def apply_systematics_preset(cfg: SimConfig) -> SimConfig:
    if cfg.systematics == "none":
        return cfg
    if cfg.systematics == "mild":
        cfg.rng_memory = 0.05
        cfg.eom_loss_amp = 0.03
        cfg.basis_eff_amp = 0.03
        cfg.timing_drift_amp = 0.4e-9
        cfg.local_bias_amp = 0.01
        return cfg
    if cfg.systematics == "hard":
        cfg.rng_memory = 0.15
        cfg.eom_loss_amp = 0.08
        cfg.basis_eff_amp = 0.08
        cfg.timing_drift_amp = 1.2e-9
        cfg.local_bias_amp = 0.03
        return cfg
    raise ValueError(f"Unknown systematics preset: {cfg.systematics}")


# -----------------------------
# X–Θ model: generate trials
# -----------------------------


def generate_theta_grid(theta_bins: int, device: torch.device) -> torch.Tensor:
    # center bins for analysis; for trials we also assign each event to a bin
    return torch.linspace(
        0, 2 * math.pi, theta_bins, device=device, dtype=torch.float64
    )


def generate_settings_with_memory(
    N: int, memory: float, device: torch.device
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Generate a,b ∈ {0,1} with optional Markov-ish memory.
    memory=0 => iid fair.
    memory>0 => tendency to repeat previous setting.
    """
    # Start with iid
    a = torch.randint(0, 2, (N,), device=device, dtype=torch.int64)
    b = torch.randint(0, 2, (N,), device=device, dtype=torch.int64)
    if memory <= 0:
        return a, b

    # With probability=memory, repeat previous
    u = torch.rand((N,), device=device, dtype=torch.float64)
    # shift previous; first stays iid
    a_prev = torch.cat([a[:1], a[:-1]], dim=0)
    b_prev = torch.cat([b[:1], b[:-1]], dim=0)
    a = torch.where(u < memory, a_prev, a)
    b = torch.where(u < memory, b_prev, b)
    return a, b


def measurement_angles(
    a: torch.Tensor, b: torch.Tensor, cfg: SimConfig, device: torch.device
) -> Tuple[torch.Tensor, torch.Tensor]:
    alpha = torch.where(
        a == 0,
        torch.tensor(cfg.a0, device=device, dtype=torch.float64),
        torch.tensor(cfg.a1, device=device, dtype=torch.float64),
    )
    beta = torch.where(
        b == 0,
        torch.tensor(cfg.b0, device=device, dtype=torch.float64),
        torch.tensor(cfg.b1, device=device, dtype=torch.float64),
    )
    return alpha, beta


def loop_label(N: int, device: torch.device) -> torch.Tensor:
    """
    Loop label L ∈ {0,1}: 0=CCW, 1=CW.
    Here we just randomize half/half per trial.
    """
    return torch.randint(0, 2, (N,), device=device, dtype=torch.int64)


def theta_assignments(N: int, theta_bins: int, device: torch.device) -> torch.Tensor:
    # assign each trial to a theta bin index
    return torch.randint(0, theta_bins, (N,), device=device, dtype=torch.int64)


def x_theta_effective_phase(
    theta: torch.Tensor, L: torch.Tensor, cfg: SimConfig, device: torch.device
) -> torch.Tensor:
    """
    This is where your X–Θ lives:
    - Standard “endpoint phase” contribution: theta_amp * sin(kθ)
    - Holonomy contribution: direction-dependent offset: ± holonomy_amp * sin(mθ)
      L=1 (CW): +, L=0 (CCW): -
    """
    base = cfg.theta_amp * torch.sin(cfg.theta_freq * theta)
    hol = cfg.holonomy_amp * torch.sin(cfg.holonomy_harm * theta)
    sign = torch.where(
        L == 1,
        torch.tensor(1.0, device=device, dtype=torch.float64),
        torch.tensor(-1.0, device=device, dtype=torch.float64),
    )
    return base + sign * hol


def simulate_outcomes(
    N: int,
    theta_bin: torch.Tensor,
    L: torch.Tensor,
    a: torch.Tensor,
    b: torch.Tensor,
    cfg: SimConfig,
    theta_grid: torch.Tensor,
    device: torch.device,
) -> Dict[str, torch.Tensor]:
    """
    Produces:
      x,y ∈ {+1,-1}
      detected flags (to model detection inefficiency / selection)
      timetags tA,tB (toy)
    """

    # Gather theta values
    theta = theta_grid[theta_bin]

    # Angles
    alpha, beta = measurement_angles(a, b, cfg, device)

    # Effective extra phase from X–Θ
    phi_theta = x_theta_effective_phase(theta, L, cfg, device)

    # Ideal singlet correlation: E = -cos(2(α-β + phi_theta))
    # The 2 factor is typical for polarization angles; adjust if you choose other conventions.
    E = -cfg.vis * torch.cos(2.0 * (alpha - beta + phi_theta))

    # Convert correlation into joint distribution with unbiased marginals:
    # P(x=y)= (1+E)/2 for x*y = +1, P(x≠y)=(1-E)/2 for x*y=-1
    # We sample x ~ ±1 uniform; then y = x with prob p_same else -x
    x = torch.where(
        torch.rand((N,), device=device, dtype=torch.float64) < 0.5,
        torch.tensor(1, device=device, dtype=torch.int8),
        torch.tensor(-1, device=device, dtype=torch.int8),
    )

    p_same = (1.0 + E) / 2.0
    same = torch.rand((N,), device=device, dtype=torch.float64) < p_same
    y = torch.where(same, x, -x).to(torch.int8)

    # Add local electronics bias (can mimic weird effects)
    if cfg.local_bias_amp > 0:
        # flip probability depends weakly on theta and setting
        bias = cfg.local_bias_amp * torch.sin(
            theta + 0.3 * a.to(torch.float64) - 0.2 * b.to(torch.float64)
        )
        p_flip_x = torch.clamp(0.5 * (bias + cfg.local_bias_amp), 0.0, 1.0)
        p_flip_y = torch.clamp(0.5 * (-bias + cfg.local_bias_amp), 0.0, 1.0)
        flip_x = torch.rand((N,), device=device, dtype=torch.float64) < p_flip_x
        flip_y = torch.rand((N,), device=device, dtype=torch.float64) < p_flip_y
        x = torch.where(flip_x, -x, x)
        y = torch.where(flip_y, -y, y)

    # Detection efficiency model: depends on basis and theta (danger!)
    detected = torch.ones((N,), device=device, dtype=torch.bool)
    if cfg.basis_eff_amp > 0 or cfg.eom_loss_amp > 0:
        # baseline efficiency
        eta0 = 0.92
        basis_term = cfg.basis_eff_amp * (
            (a.to(torch.float64) - 0.5) + (b.to(torch.float64) - 0.5)
        )
        eom_term = cfg.eom_loss_amp * torch.sin(theta)
        eta = torch.clamp(eta0 + basis_term - eom_term, 0.02, 0.999)
        detected = torch.rand((N,), device=device, dtype=torch.float64) < eta

    # Timing model with coincidence-window selection bias
    # Create timetags with jitter that depends on setting and theta; then select coincidences
    tA = torch.randn((N,), device=device, dtype=torch.float64) * cfg.time_jitter_sigma
    tB = torch.randn((N,), device=device, dtype=torch.float64) * cfg.time_jitter_sigma

    if cfg.timing_drift_amp > 0:
        # slow drift over trial index
        idx = torch.arange(N, device=device, dtype=torch.float64)
        drift = cfg.timing_drift_amp * torch.sin(2 * math.pi * idx / max(N, 1))
        tB = tB + drift

    # jitter bias depends on theta and settings (classic gotcha)
    jitter_bias = 0.4e-9 * torch.sin(
        theta + 0.8 * a.to(torch.float64) - 0.6 * b.to(torch.float64)
    )
    tB = tB + jitter_bias

    # Coincidence condition (toy): |tA - tB| < window
    coinc = torch.abs(tA - tB) < cfg.coincidence_window

    # Final keep mask
    keep = detected & coinc

    return {
        "theta_bin": theta_bin,
        "loop": L,
        "a": a,
        "b": b,
        "x": x,
        "y": y,
        "keep": keep,
        "E_true": E,  # useful for sanity
    }


# -----------------------------
# Statistics: E_ab and S per theta-bin and loop
# -----------------------------


def compute_counts(
    data: Dict[str, torch.Tensor], theta_bins: int, device: torch.device
) -> torch.Tensor:
    """
    Returns N_abxy[loop, theta, a, b, x, y] where x,y index: 0->+1, 1->-1.
    """
    keep = data["keep"]
    loop = data["loop"][keep].to(torch.int64)
    th = data["theta_bin"][keep].to(torch.int64)
    a = data["a"][keep].to(torch.int64)
    b = data["b"][keep].to(torch.int64)
    x = data["x"][keep].to(torch.int64)
    y = data["y"][keep].to(torch.int64)

    # map ±1 to 0/1
    xi = torch.where(
        x == 1, torch.tensor(0, device=device), torch.tensor(1, device=device)
    )
    yi = torch.where(
        y == 1, torch.tensor(0, device=device), torch.tensor(1, device=device)
    )

    # Flatten index
    # dims: loop(2), theta(T), a(2), b(2), x(2), y(2)
    T = theta_bins
    idx = ((((loop * T + th) * 2 + a) * 2 + b) * 2 + xi) * 2 + yi

    out = torch.zeros((2 * T * 2 * 2 * 2 * 2,), device=device, dtype=torch.int64)
    ones = torch.ones_like(idx, dtype=torch.int64)
    out.scatter_add_(0, idx, ones)
    out = out.view(2, T, 2, 2, 2, 2)
    return out


def counts_to_E_S(
    N_abxy: torch.Tensor, eps: float = 1e-12
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Input: N_abxy[loop, theta, a, b, x, y]
    Output:
      E[loop, theta, a, b]
      S[loop, theta]
      N_ab[loop, theta, a, b]
    """
    # N_ab = sum over x,y
    N_ab = N_abxy.sum(dim=(-1, -2)).to(torch.float64)  # [2, T, 2, 2]

    # Convert x,y indices back to ±1 for expectation:
    # x index 0 => +1, 1 => -1
    # y index 0 => +1, 1 => -1
    # E = (N++ + N-- - N+- - N-+) / N
    Npp = N_abxy[..., 0, 0].to(torch.float64)
    Npm = N_abxy[..., 0, 1].to(torch.float64)
    Nmp = N_abxy[..., 1, 0].to(torch.float64)
    Nmm = N_abxy[..., 1, 1].to(torch.float64)

    E = (Npp + Nmm - Npm - Nmp) / torch.clamp(N_ab, min=eps)

    # CHSH S = E00 + E01 + E10 - E11
    S = E[..., 0, 0] + E[..., 0, 1] + E[..., 1, 0] - E[..., 1, 1]  # [2, T]

    return E, S, N_ab


def approx_se_E(E: torch.Tensor, N: torch.Tensor, eps: float = 1.0) -> torch.Tensor:
    """
    Approximate standard error for correlation estimator:
      SE ≈ sqrt((1 - E^2) / N)
    """
    N_eff = torch.clamp(N, min=eps)
    return torch.sqrt(torch.clamp(1.0 - E**2, min=0.0) / N_eff)


def approx_se_S(E: torch.Tensor, N_ab: torch.Tensor) -> torch.Tensor:
    """
    Propagate SE across the four E_ab terms in CHSH.
    """
    seE = approx_se_E(E, N_ab)
    # sum variances of four terms
    v = (
        seE[..., 0, 0] ** 2
        + seE[..., 0, 1] ** 2
        + seE[..., 1, 0] ** 2
        + seE[..., 1, 1] ** 2
    )
    return torch.sqrt(v)


# -----------------------------
# Block bootstrap + permutation test for ΔS
# -----------------------------


def block_indices(N: int, block_size: int, device: torch.device) -> torch.Tensor:
    n_blocks = N // block_size
    idx = torch.arange(n_blocks * block_size, device=device)
    return idx.view(n_blocks, block_size)


def compute_deltaS_from_trials(
    data: Dict[str, torch.Tensor], theta_bins: int, device: torch.device
) -> torch.Tensor:
    N_abxy = compute_counts(data, theta_bins, device)
    _, S, _ = counts_to_E_S(N_abxy)
    deltaS = S[1] - S[0]  # CW - CCW
    return deltaS


def block_bootstrap_deltaS(
    data: Dict[str, torch.Tensor], cfg: SimConfig, device: torch.device
) -> torch.Tensor:
    """
    Bootstrap over blocks of kept trials (not perfect, but robust vs mild dependence).
    Returns bootstrap samples of deltaS[theta].
    """
    keep = data["keep"]
    # pull kept trials
    kept_idx = torch.nonzero(keep, as_tuple=False).flatten()
    if kept_idx.numel() < cfg.block_size * 3:
        # too few; fallback: no bootstrap
        deltaS = compute_deltaS_from_trials(data, cfg.theta_bins, device)
        return deltaS.unsqueeze(0)

    # Build blocks in terms of kept indices (preserve order)
    M = kept_idx.numel()
    n_blocks = M // cfg.block_size
    kept_idx = kept_idx[: n_blocks * cfg.block_size]
    blocks = kept_idx.view(n_blocks, cfg.block_size)

    samples = []
    for _ in range(cfg.bootstrap_rounds):
        pick = torch.randint(0, n_blocks, (n_blocks,), device=device)
        sample_idx = blocks[pick].reshape(-1)

        # Build resampled dataset (same keys)
        d = {}
        for k, v in data.items():
            if k in ("E_true",):  # not needed
                continue
            d[k] = v[sample_idx]
        # keep mask is all True now (since sampled from kept)
        d["keep"] = torch.ones_like(d["a"], dtype=torch.bool, device=device)

        samples.append(
            compute_deltaS_from_trials(d, cfg.theta_bins, device).unsqueeze(0)
        )

    return torch.cat(samples, dim=0)  # [B, theta]


def permutation_test_deltaS(
    data: Dict[str, torch.Tensor], cfg: SimConfig, device: torch.device
) -> Tuple[torch.Tensor, float]:
    """
    Permutation: shuffle loop labels within blocks (keeps time structure).
    Test statistic: max |ΔS(θ)| over θ-bins.
    """
    keep = data["keep"]
    kept_idx = torch.nonzero(keep, as_tuple=False).flatten()
    if kept_idx.numel() < cfg.block_size * 3:
        obs = compute_deltaS_from_trials(data, cfg.theta_bins, device)
        stat_obs = torch.max(torch.abs(obs)).item()
        return obs, float("nan")

    M = kept_idx.numel()
    n_blocks = M // cfg.block_size
    kept_idx = kept_idx[: n_blocks * cfg.block_size]
    blocks = kept_idx.view(n_blocks, cfg.block_size)

    # Observed
    obs = compute_deltaS_from_trials(data, cfg.theta_bins, device)
    stat_obs = torch.max(torch.abs(obs))

    # Permute within blocks: shuffle loop labels
    stats = []
    loop_full = data["loop"].clone()
    for _ in range(cfg.perm_rounds):
        loop_perm = loop_full.clone()
        for bi in range(n_blocks):
            ix = blocks[bi]
            # random permutation of loop labels inside block
            perm = ix[torch.randperm(ix.numel(), device=device)]
            loop_perm[ix] = loop_full[perm]
        d = dict(data)
        d["loop"] = loop_perm
        dstat = torch.max(
            torch.abs(compute_deltaS_from_trials(d, cfg.theta_bins, device))
        )
        stats.append(dstat)

    stats = torch.stack(stats)
    # p-value: fraction >= observed
    p = (torch.sum(stats >= stat_obs) + 1).item() / (cfg.perm_rounds + 1)
    return obs, p


# -----------------------------
# FFT + surrogates
# -----------------------------


def fft_spectrum(y: np.ndarray) -> np.ndarray:
    # y is real, length T
    Y = np.fft.rfft(y - np.mean(y))
    return np.abs(Y)


def surrogate_theta_shuffle(y: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    ys = y.copy()
    rng.shuffle(ys)
    return ys


def surrogate_phase_scramble(y: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    # Keep magnitudes, randomize phases (for rfft components excluding DC/Nyquist)
    Y = np.fft.rfft(y - np.mean(y))
    mag = np.abs(Y)
    phase = np.angle(Y)
    # random phases for components 1..end-1
    if len(Y) > 2:
        new_phase = phase.copy()
        new_phase[1:-1] = rng.uniform(0, 2 * np.pi, size=len(Y) - 2)
        Y2 = mag * np.exp(1j * new_phase)
    else:
        Y2 = Y
    y2 = np.fft.irfft(Y2, n=len(y))
    return y2.real


# -----------------------------
# Main runner
# -----------------------------


def run(cfg: SimConfig):
    cfg = apply_systematics_preset(cfg)
    set_seed(cfg.seed)
    device = pick_device(cfg.use_gpu)

    print(f"[device] {device}")
    print(
        f"[cfg] N={cfg.N:,} theta_bins={cfg.theta_bins} systematics={cfg.systematics}"
    )
    print(
        f"[X–Θ] theta_amp={cfg.theta_amp}, holonomy_amp={cfg.holonomy_amp}, vis={cfg.vis}"
    )

    # Build theta grid
    theta_grid = generate_theta_grid(cfg.theta_bins, device)

    # Generate trials
    a, b = generate_settings_with_memory(cfg.N, cfg.rng_memory, device)
    L = loop_label(cfg.N, device)
    th_bin = theta_assignments(cfg.N, cfg.theta_bins, device)

    data = simulate_outcomes(cfg.N, th_bin, L, a, b, cfg, theta_grid, device)

    # Counts -> E,S
    N_abxy = compute_counts(data, cfg.theta_bins, device)
    E, S, N_ab = counts_to_E_S(N_abxy)
    seS = approx_se_S(E, N_ab)

    # ΔS = CW - CCW
    deltaS = S[1] - S[0]
    se_deltaS = torch.sqrt(seS[1] ** 2 + seS[0] ** 2)

    # Print headline results
    # Aggregate over theta (weighted)
    N_theta = N_ab.sum(dim=(-1, -2))  # [2, T]
    w = torch.clamp(N_theta, min=1.0)
    S_weighted = (S * w).sum(dim=1) / w.sum(dim=1)
    dS_weighted = S_weighted[1] - S_weighted[0]

    print("\n[headline]")
    print(f"  kept trials: {int(data['keep'].sum().item()):,} / {cfg.N:,}")
    print(f"  S_weighted CCW: {S_weighted[0].item(): .4f}")
    print(f"  S_weighted  CW: {S_weighted[1].item(): .4f}")
    print(f"  ΔS_weighted (CW-CCW): {dS_weighted.item(): .4f}")

    # Permutation test
    obs_deltaS, p_perm = permutation_test_deltaS(data, cfg, device)
    max_abs = torch.max(torch.abs(obs_deltaS)).item()
    print("\n[permutation test]")
    print(f"  stat = max_theta |ΔS(θ)| = {max_abs:.5f}")
    print(f"  p-value (block-permute loop labels) = {p_perm}")

    # Bootstrap SE
    boot = block_bootstrap_deltaS(data, cfg, device)  # [B, T] or [1,T]
    boot_np = boot.detach().cpu().numpy()
    deltaS_np = deltaS.detach().cpu().numpy()
    se_boot = (
        np.std(boot_np, axis=0)
        if boot_np.shape[0] > 1
        else se_deltaS.detach().cpu().numpy()
    )

    # FFT on one representative channel: use ΔS(θ) itself (most direct holonomy signal)
    rng = np.random.default_rng(cfg.seed + 999)
    spec_obs = fft_spectrum(deltaS_np)

    # Surrogates
    K = 200
    spec_shuffle = np.zeros((K, spec_obs.size))
    spec_phase = np.zeros((K, spec_obs.size))
    for k in range(K):
        ys = surrogate_theta_shuffle(deltaS_np, rng)
        yp = surrogate_phase_scramble(deltaS_np, rng)
        spec_shuffle[k] = fft_spectrum(ys)
        spec_phase[k] = fft_spectrum(yp)

    # Simple “ladder” heuristic: peaks beyond DC
    # Report top few frequencies where observed exceeds 99th percentile of shuffle
    thr = np.quantile(spec_shuffle, 0.99, axis=0)
    sig_bins = np.where(spec_obs > thr)[0]
    # ignore DC
    sig_bins = sig_bins[sig_bins > 0]

    print("\n[FFT sanity]")
    if sig_bins.size == 0:
        print(
            "  No FFT bins exceed 99% θ-shuffle surrogate threshold (good: no fake ladder)."
        )
    else:
        top = sig_bins[:10]
        print("  FFT bins exceeding 99% θ-shuffle threshold (first 10):", top.tolist())
        print("  (bin index k corresponds to harmonic k over the θ grid)")

    # Save results to NPZ (easy to plot elsewhere)
    out = {
        "theta_grid": theta_grid.detach().cpu().numpy(),
        "E": E.detach().cpu().numpy(),  # [loop, theta, a, b]
        "S": S.detach().cpu().numpy(),  # [loop, theta]
        "seS": seS.detach().cpu().numpy(),
        "deltaS": deltaS_np,  # [theta]
        "se_deltaS_approx": se_deltaS.detach().cpu().numpy(),
        "se_deltaS_boot": se_boot,
        "p_perm": p_perm,
        "spec_obs": spec_obs,
        "spec_shuffle_q99": thr,
        "cfg": vars(cfg),
    }
    os.makedirs("out", exist_ok=True)
    path = os.path.join(
        "out", f"x_theta_sim_{cfg.systematics}_N{cfg.N}_T{cfg.theta_bins}.npz"
    )
    np.savez(path, **out)
    print(f"\n[saved] {path}")

    # Print a few θ-bin lines (so you can eyeball quickly without plotting)
    print("\n[sample θ bins]")
    for k in np.linspace(0, cfg.theta_bins - 1, num=min(8, cfg.theta_bins), dtype=int):
        print(
            f"  k={k:3d} θ={out['theta_grid'][k]:.3f}  ΔS={out['deltaS'][k]: .5f}  SE≈{out['se_deltaS_boot'][k]:.5f}"
        )


def parse_args() -> SimConfig:
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=1_000_000)
    p.add_argument("--theta-bins", type=int, default=64)
    p.add_argument("--block-size", type=int, default=20_000)
    p.add_argument("--bootstrap-rounds", type=int, default=300)
    p.add_argument("--perm-rounds", type=int, default=400)
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--use-gpu", type=int, default=1)

    p.add_argument(
        "--systematics", type=str, default="none", choices=["none", "mild", "hard"]
    )

    # X–Θ knobs
    p.add_argument("--vis", type=float, default=0.98)
    p.add_argument("--theta-amp", type=float, default=0.7)
    p.add_argument("--theta-freq", type=int, default=1)
    p.add_argument("--holonomy-amp", type=float, default=0.20)
    p.add_argument("--holonomy-harm", type=int, default=1)

    args = p.parse_args()
    cfg = SimConfig(
        N=args.N,
        theta_bins=args.theta_bins,
        block_size=args.block_size,
        bootstrap_rounds=args.bootstrap_rounds,
        perm_rounds=args.perm_rounds,
        seed=args.seed,
        use_gpu=args.use_gpu,
        systematics=args.systematics,
        vis=args.vis,
        theta_amp=args.theta_amp,
        theta_freq=args.theta_freq,
        holonomy_amp=args.holonomy_amp,
        holonomy_harm=args.holonomy_harm,
    )
    return cfg


if __name__ == "__main__":
    cfg = parse_args()
    run(cfg)
