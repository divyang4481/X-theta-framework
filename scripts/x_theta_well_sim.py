#!/usr/bin/env python3
"""
X-θ Drift Simulation on The Well (single-file script)

What it does:
- Loads a Well dataset window (streaming from HF by default, or local HDF5)
- Extracts a 2D velocity field (u,v) from the dataset channels
- Simulates tracer particles with an added X-θ drift term:
      x_dot = u(x,t) + kappa * ∇A_theta(x) * theta_dot
- Runs winding test N in [0,1,2,4,8], plots Δx vs N, and a scatter “money shot”.

Requirements:
  pip install the-well torch numpy matplotlib
Optional for local validation/repair:
  pip install h5py

Examples:
  # Streaming mode (recommended)
  python x_theta_well_sim.py --mode stream

  # Local mode: download + validate + repair if needed, then run
  python x_theta_well_sim.py --mode local --base-path C:\\workspace\\Physics\\X-theta-framework\\well_data --download --repair

  # Use a different dataset window
  python x_theta_well_sim.py --mode stream --window-index 5

Notes:
- If your dataset doesn’t expose "velocity_x"/"velocity_y", the script will try fallbacks.
- Some Well datasets may not have velocities at all. Then pick another dataset.
"""

import argparse
import glob
import math
import os
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np

# torch + matplotlib
import torch
import matplotlib.pyplot as plt


def require_pkg(name: str, pip_name: Optional[str] = None) -> None:
    """Friendly error message if import isn't available."""
    try:
        __import__(name)
    except ImportError as e:
        pip_name = pip_name or name
        raise SystemExit(
            f"Missing package '{name}'. Install with:\n"
            f"  python -m pip install {pip_name}\n"
        ) from e


# the-well (PyPI name "the-well", import is "the_well")
require_pkg("the_well", "the-well")
from the_well.data import WellDataset


def resolve_store_path(base_path: Path) -> Path:
    """
    The Well downloader sometimes creates base_path/datasets/<dataset>/...
    This resolves the correct root for WellDataset.
    """
    candidate = base_path / "datasets"
    return candidate if candidate.exists() else base_path


def list_hdf5_files(train_dir: Path) -> List[str]:
    return sorted(
        glob.glob(str(train_dir / "*.h5")) + glob.glob(str(train_dir / "*.hdf5"))
    )


def validate_hdf5_files(files: List[str]) -> List[Tuple[str, str]]:
    """
    Returns list of (bad_file, error_msg).
    Requires h5py.
    """
    require_pkg("h5py")
    import h5py  # noqa: F401

    bad = []
    for f in files:
        try:
            import h5py
            with h5py.File(f, "r") as hf:
                _ = list(hf.keys())[:5]  # touch metadata/keys
        except Exception as e:
            bad.append((f, str(e)))
    return bad


def download_dataset(base_path: Path, dataset: str, split: str) -> None:
    require_pkg("the_well", "the-well")
    from the_well.utils.download import well_download
    well_download(base_path=str(base_path), dataset=dataset, split=split)


def flatten_field_names(field_names_dict) -> List[str]:
    out = []
    for _, group in field_names_dict.items():
        out.extend(group)
    return out


def find_velocity_indices(names: List[str]) -> Tuple[int, int]:
    """
    Try common channel naming patterns.
    """
    # common exact pairs
    pairs = [
        ("velocity_x", "velocity_y"),
        ("u", "v"),
        ("vel_x", "vel_y"),
        ("vx", "vy"),
    ]
    for a, b in pairs:
        if a in names and b in names:
            return names.index(a), names.index(b)

    # fallback: first two containing "velocity"
    vel = [i for i, n in enumerate(names) if "velocity" in n.lower()]
    if len(vel) >= 2:
        return vel[0], vel[1]

    raise RuntimeError(
        "No velocity channels found. Field names were:\n" + "\n".join(names)
    )


def finite_diff_grad(A: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    A: [H,W], periodic central differences
    Returns dA_dx, dA_dy each [H,W].
    """
    A_left = torch.roll(A, +1, 1)
    A_right = torch.roll(A, -1, 1)
    dA_dx = 0.5 * (A_right - A_left)

    A_up = torch.roll(A, +1, 0)
    A_down = torch.roll(A, -1, 0)
    dA_dy = 0.5 * (A_down - A_up)
    return dA_dx, dA_dy


def bilinear_sample(field_hw: torch.Tensor, x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    """
    field_hw: [H,W] or [H,W,C]
    x,y: [P] float in index space
    periodic boundaries
    """
    H, W = field_hw.shape[0], field_hw.shape[1]
    x = x % W
    y = y % H

    x0 = torch.floor(x).long()
    y0 = torch.floor(y).long()
    x1 = (x0 + 1) % W
    y1 = (y0 + 1) % H

    wx = (x - x0.float())
    wy = (y - y0.float())

    def g(ix, iy):
        return field_hw[iy, ix] if field_hw.ndim == 2 else field_hw[iy, ix, :]

    f00 = g(x0, y0)
    f10 = g(x1, y0)
    f01 = g(x0, y1)
    f11 = g(x1, y1)

    if field_hw.ndim == 3:
        f0 = f00 * (1 - wx).unsqueeze(-1) + f10 * wx.unsqueeze(-1)
        f1 = f01 * (1 - wx).unsqueeze(-1) + f11 * wx.unsqueeze(-1)
        return f0 * (1 - wy).unsqueeze(-1) + f1 * wy.unsqueeze(-1)
    else:
        f0 = f00 * (1 - wx) + f10 * wx
        f1 = f01 * (1 - wx) + f11 * wx
        return f0 * (1 - wy) + f1 * wy


def simulate_tracers(
    vel_seq: torch.Tensor,
    dA_dx: torch.Tensor,
    dA_dy: torch.Tensor,
    N_wind: int,
    kappa: float,
    dt: float,
    n_particles: int = 4000,
    noise_D: float = 0.0,
    seed: int = 0,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    vel_seq: [T,H,W,2]
    X-θ drift: kappa * ∇Aθ(x) * θ̇
    """
    torch.manual_seed(seed)
    T, H, W, _ = vel_seq.shape

    x = torch.rand(n_particles, device=vel_seq.device) * (W - 1)
    y = torch.rand(n_particles, device=vel_seq.device) * (H - 1)
    x0, y0 = x.clone(), y.clone()

    theta_dot = (2.0 * math.pi * N_wind) / (max(1, T - 1) * dt)

    gradA = torch.stack([dA_dx, dA_dy], dim=-1)  # [H,W,2]

    for t in range(T):
        uv_p = bilinear_sample(vel_seq[t], x, y)       # [P,2]
        gradA_p = bilinear_sample(gradA, x, y)         # [P,2]

        drift = kappa * gradA_p * theta_dot
        step = (uv_p + drift) * dt

        if noise_D > 0:
            step = step + torch.randn_like(step) * math.sqrt(2 * noise_D * dt)

        x = x + step[:, 0]
        y = y + step[:, 1]

    dx = x - x0
    dy = y - y0
    return dx, dy, x, y


def load_dataset(args) -> Tuple[WellDataset, dict]:
    """
    Returns (ds, item) for the chosen mode.
    """
    if args.mode == "stream":
        ds = WellDataset(
            well_base_path="hf://datasets/polymathic-ai/",
            well_dataset_name=args.dataset,
            well_split_name=args.split,
            n_steps_input=args.steps_in,
            n_steps_output=args.steps_out,
            use_normalization=False,
        )
        item = ds[args.window_index]
        return ds, item

    # local mode
    base_path = Path(args.base_path).resolve()
    base_path.mkdir(parents=True, exist_ok=True)

    if args.download:
        download_dataset(base_path, args.dataset, args.split)

    store_path = resolve_store_path(base_path)
    train_dir = store_path / args.dataset / "data" / args.split

    if not train_dir.exists():
        raise SystemExit(
            f"Local train_dir does not exist:\n  {train_dir}\n"
            f"Try running with --download, or check your base path."
        )

    files = list_hdf5_files(train_dir)
    if len(files) == 0:
        raise SystemExit(
            f"No .h5/.hdf5 files found in:\n  {train_dir}\n"
            f"Try --download or check network/proxy."
        )

    if args.repair:
        bad = validate_hdf5_files(files)
        if bad:
            print(f"[repair] Found {len(bad)} corrupt HDF5 files. Deleting them...")
            for f, msg in bad[:10]:
                print(" -", os.path.basename(f), "|", msg.splitlines()[0][:120])
            for f, _ in bad:
                try:
                    os.remove(f)
                except Exception as e:
                    print("Could not delete:", f, "error:", e)
            print("[repair] Re-downloading to replace deleted files...")
            download_dataset(base_path, args.dataset, args.split)

    # Re-resolve in case download created /datasets
    store_path = resolve_store_path(base_path)

    ds = WellDataset(
        well_base_path=str(store_path),
        well_dataset_name=args.dataset,
        well_split_name=args.split,
        n_steps_input=args.steps_in,
        n_steps_output=args.steps_out,
        use_normalization=False,
    )
    item = ds[args.window_index]
    return ds, item


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["stream", "local"], default="stream",
                   help="stream = HF streaming; local = local HDF5 download")
    p.add_argument("--base-path", default=str((Path.cwd() / "well_data").resolve()),
                   help="Local base path (only used in --mode local)")
    p.add_argument("--download", action="store_true",
                   help="In local mode: download dataset split into base-path")
    p.add_argument("--repair", action="store_true",
                   help="In local mode: validate & delete corrupt HDF5 then re-download (needs h5py)")

    p.add_argument("--dataset", default="turbulent_radiative_layer_2D")
    p.add_argument("--split", default="train")
    p.add_argument("--window-index", type=int, default=0)
    p.add_argument("--steps-in", type=int, default=8)
    p.add_argument("--steps-out", type=int, default=8)

    p.add_argument("--kappa", type=float, default=1.0, help="X-θ coupling strength")
    p.add_argument("--alpha", type=float, default=1.0, help="Aθ ramp strength")
    p.add_argument("--noise-d", type=float, default=0.0, help="Diffusion coefficient")
    p.add_argument("--particles", type=int, default=4000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--N", default="0,1,2,4,8", help="Comma-separated winding list")
    p.add_argument("--no-scatter", action="store_true", help="Skip endpoint scatter plot")

    args = p.parse_args()

    print("[env] torch:", torch.__version__, "| cuda:", torch.cuda.is_available())
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print("[env] device:", device)
    print("[run] mode:", args.mode, "| dataset:", args.dataset, "| split:", args.split, "| window:", args.window_index)

    # Load dataset window
    ds, item = load_dataset(args)

    # Build velocity sequence
    field_names = flatten_field_names(ds.metadata.field_names)
    ix, iy = find_velocity_indices(field_names)
    print("[data] velocity channels:", field_names[ix], field_names[iy])

    fields = torch.cat([item["input_fields"], item["output_fields"]], dim=0).to(torch.float32)  # [T,H,W,F]
    vel = fields[..., [ix, iy]].to(device)  # [T,H,W,2]
    T, H, W, _ = vel.shape
    print("[data] vel shape:", tuple(vel.shape))

    # dt from time grid if present
    dt = 1.0
    t_in = item.get("input_time_grid", None)
    t_out = item.get("output_time_grid", None)
    if t_in is not None and t_out is not None:
        t_full = torch.cat([t_in, t_out], dim=0).to(torch.float32)
        if len(t_full) >= 2:
            dt = float((t_full[1] - t_full[0]).abs().cpu().item()) or 1.0
    print("[data] dt:", dt)

    # Aθ = alpha * x ramp along width direction
    x_coords = torch.linspace(0, 1, W, device=device)
    A = args.alpha * x_coords.unsqueeze(0).repeat(H, 1)  # [H,W]
    dA_dx, dA_dy = finite_diff_grad(A)

    # Parse N list
    N_list = [int(s.strip()) for s in args.N.split(",") if s.strip()]
    print("[run] winding N list:", N_list)

    results = []
    for N in N_list:
        dx_c, dy_c, _, _ = simulate_tracers(
            vel, dA_dx, dA_dy, N_wind=N, kappa=0.0, dt=dt,
            n_particles=args.particles, noise_D=args.noise_d, seed=args.seed
        )
        dx_x, dy_x, _, _ = simulate_tracers(
            vel, dA_dx, dA_dy, N_wind=N, kappa=args.kappa, dt=dt,
            n_particles=args.particles, noise_D=args.noise_d, seed=args.seed
        )

        mean_dx_control = dx_c.mean().item()
        mean_dx_xtheta = dx_x.mean().item()
        delta = mean_dx_xtheta - mean_dx_control
        results.append((N, mean_dx_control, mean_dx_xtheta, delta))
        print(f"[N={N:>2}] mean_dx(control)={mean_dx_control:+.6f} | mean_dx(Xθ)={mean_dx_xtheta:+.6f} | Δ={delta:+.6f}")

    Ns = np.array([r[0] for r in results], dtype=float)
    deltas = np.array([r[3] for r in results], dtype=float)

    # Plot Δx vs N
    plt.figure()
    plt.plot(Ns, deltas, marker="o")
    plt.xlabel("Winding number N")
    plt.ylabel("Mean Δx (Xθ - control)")
    plt.title(f"X-θ Winding Test on {args.dataset} ({args.mode})")
    plt.grid(True)
    plt.show()

    # Fit linearity
    Afit = np.vstack([Ns, np.ones_like(Ns)]).T
    m, b = np.linalg.lstsq(Afit, deltas, rcond=None)[0]
    pred = m * Ns + b
    ss_res = np.sum((deltas - pred) ** 2)
    ss_tot = np.sum((deltas - deltas.mean()) ** 2) + 1e-12
    r2 = 1 - ss_res / ss_tot
    print(f"[fit] Δx ≈ {m:.6f} * N + {b:.6f} | R^2={r2:.6f}")

    # Optional scatter
    if not args.no_scatter and (len(N_list) > 0):
        N_show = max(N_list)
        _, _, x_c, y_c = simulate_tracers(
            vel, dA_dx, dA_dy, N_wind=N_show, kappa=0.0, dt=dt,
            n_particles=min(3000, args.particles), noise_D=args.noise_d, seed=args.seed
        )
        _, _, x_x, y_x = simulate_tracers(
            vel, dA_dx, dA_dy, N_wind=N_show, kappa=args.kappa, dt=dt,
            n_particles=min(3000, args.particles), noise_D=args.noise_d, seed=args.seed
        )

        plt.figure(figsize=(6, 6))
        plt.scatter(x_c.detach().cpu(), y_c.detach().cpu(), s=1, alpha=0.35, label="Control")
        plt.scatter(x_x.detach().cpu(), y_x.detach().cpu(), s=1, alpha=0.35, label="X-θ")
        plt.gca().invert_yaxis()
        plt.title(f"Tracer endpoints (N={N_show})")
        plt.legend()
        plt.show()


if __name__ == "__main__":
    main()
