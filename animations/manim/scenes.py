"""Reusable Manim scenes that visualize the X-Theta simulations."""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Tuple

import numpy as np
from manim import (
    BLUE,
    BLUE_E,
    DOWN,
    GREEN,
    LEFT,
    ORANGE,
    RED,
    RIGHT,
    UP,
    WHITE,
    Axes,
    DecimalNumber,
    Dot,
    FadeIn,
    LaggedStart,
    Line,
    MathTex,
    Scene,
    VGroup,
    ValueTracker,
    always_redraw,
)

# Ensure Manim can import the project modules when run from repo root
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

from src.sim.double_slit import simulate_interference  # noqa: E402
from src.sim.xtheta_sim import XThetaQuantumSim  # noqa: E402
import src.sim.xtheta_monte_carlo as xtheta_mc_module  # noqa: E402
from src.sim.xtheta_monte_carlo import XThetaMonteCarlo  # noqa: E402

HBAR = 1.0
xtheta_mc_module.hbar = HBAR


def _double_slit_profile(t_value: float) -> Tuple[np.ndarray, np.ndarray]:
    """Helper that returns (x, I(x, t)) for the toy double slit model."""
    x_vals = np.linspace(-5.0, 5.0, 400)
    params = dict(slit_sep=1.6, wavelength=0.55, L=12.0, theta_amp=0.9, theta_freq=0.6)
    intensities = simulate_interference(x_vals, t=t_value, **params)
    return x_vals, intensities


def _run_xtheta_trajectory(total_time: float = 4.0, steps: int = 120) -> Tuple[np.ndarray, np.ndarray]:
    """Generates a centroid trajectory using the coarse quantum simulator."""
    sim = XThetaQuantumSim(N_y=160, N_theta=48, L_y=12.0, mass=1.0, Inertia=0.25, q_theta=1.0)
    sim.initialize_wavepacket(y0=0.0, sigma_y=0.6, l_mode=3)

    dt = total_time / (steps - 1)
    times = [0.0]
    centroids = [sim.measure_centroid()]
    for idx in range(1, steps):
        sim.run_step(dt, A_grad=0.45)
        times.append(idx * dt)
        centroids.append(sim.measure_centroid())
    return np.array(times), np.array(centroids)


def _collect_monte_carlo_samples(trials: int = 24) -> Tuple[np.ndarray, np.ndarray]:
    """Runs a handful of Monte Carlo draws for the drift statistics plot."""
    rng = np.random.default_rng(seed=23)
    sim = XThetaMonteCarlo(N_y=160, N_theta=64, L_y=24.0, Inertia=0.4, q_theta=1.0)

    l_modes = []
    drifts = []
    for _ in range(trials):
        psi_init, l_val = sim.generate_thermal_state(temp=2.0, y0=0.0, sigma_y=0.8)
        traj = sim.evolve(psi_init, steps=90, dt=0.04, A_grad=0.4)
        # Mild randomization keeps the scatter visually interesting
        jitter = 0.02 * rng.normal()
        l_modes.append(float(l_val))
        drifts.append(float(traj[-1] - traj[0] + jitter))
    return np.array(l_modes), np.array(drifts)


class DoubleSlitThetaScene(Scene):
    """Animated double-slit intensity with a driven theta-phase."""

    def construct(self) -> None:
        axes = Axes(x_range=[-5, 5, 1], y_range=[0, 1.2, 0.2], x_length=11, y_length=4)
        axes.to_edge(UP)
        labels = axes.get_axis_labels("x", "I(x)")

        time_tracker = ValueTracker(0.0)

        # Redraw the intensity envelope every frame as the theta-phase evolves
        def intensity_curve():
            x_vals, intensities = _double_slit_profile(time_tracker.get_value())
            return axes.plot_line_graph(x_vals, intensities, line_color=ORANGE, add_vertex_dots=False)

        curve = always_redraw(intensity_curve)

        phi_text = MathTex(r"\phi_\theta(t) = A_\theta \sin(2\pi f_\theta t)")
        phi_text.next_to(axes, DOWN, buff=0.8)
        time_display = VGroup(
            MathTex("t ="),
            always_redraw(lambda: DecimalNumber(time_tracker.get_value(), num_decimal_places=2)),
        )
        time_display.arrange(RIGHT)
        time_display.next_to(phi_text, DOWN)

        self.play(FadeIn(axes), FadeIn(labels))
        self.play(FadeIn(curve), FadeIn(phi_text), FadeIn(time_display))
        self.play(time_tracker.animate.set_value(8.0), run_time=8.0)
        self.play(time_tracker.animate.set_value(16.0), run_time=8.0)
        self.wait(1.5)


class XThetaDriftScene(Scene):
    """Shows the centroid drift that arises from the cross Hall coupling."""

    def construct(self) -> None:
        times, centroids = _run_xtheta_trajectory()
        t_range = [0, float(times[-1]), 1.0]
        y_margin = 0.1
        y_min = float(np.min(centroids)) - y_margin
        y_max = float(np.max(centroids)) + y_margin
        axes = Axes(x_range=t_range, y_range=[y_min, y_max, 0.2], x_length=11, y_length=5)
        axes.to_edge(UP)
        axes_labels = axes.get_axis_labels("t", r"\langle y \rangle")

        graph = axes.plot_line_graph(times, centroids, add_vertex_dots=False, line_color=BLUE)
        tracker = ValueTracker(0.0)

        # Interpolate so the tracer moves smoothly, not stepwise between samples
        def moving_dot():
            t_val = tracker.get_value()
            y_val = float(np.interp(t_val, times, centroids))
            return Dot(color=WHITE).move_to(axes.coords_to_point(t_val, y_val))

        dot = always_redraw(moving_dot)

        legend = VGroup(
            MathTex(r"A_y \propto \theta\text{-momentum}", color=BLUE).scale(0.8),
            MathTex(r"\Delta y > 0 \Rightarrow\ \text{drift}", color=WHITE).scale(0.8),
        ).arrange(DOWN, aligned_edge=LEFT)
        legend.next_to(axes, DOWN, buff=0.8)

        self.play(FadeIn(axes), FadeIn(axes_labels))
        self.play(FadeIn(graph), FadeIn(dot))
        self.play(tracker.animate.set_value(float(times[-1])), run_time=6.0)
        self.play(FadeIn(legend))
        self.wait(2.0)


class XThetaMonteCarloScene(Scene):
    """Scatter plot of Monte Carlo drifts vs internal rotor mode."""

    def construct(self) -> None:
        l_modes, drifts = _collect_monte_carlo_samples()
        x_range = [min(l_modes) - 0.5, max(l_modes) + 0.5, 1]
        y_pad = 0.2
        y_range = [float(np.min(drifts)) - y_pad, float(np.max(drifts)) + y_pad, 0.2]

        axes = Axes(x_range=x_range, y_range=y_range, x_length=10, y_length=5)
        axes.to_edge(UP)
        axes_labels = axes.get_axis_labels("l", r"\Delta y")

        dots = VGroup()
        for l_val, drift in zip(l_modes, drifts):
            color = RED if drift > 0 else BLUE_E
            dots.add(Dot(axes.coords_to_point(l_val, drift), color=color, radius=0.06))

        corr_text = MathTex(r"\rho(l, \Delta y) = %.2f" % np.corrcoef(l_modes, drifts)[0, 1])
        corr_text.next_to(axes, DOWN, buff=0.6)

        guide_line = Line(axes.coords_to_point(x_range[0], 0), axes.coords_to_point(x_range[1], 0), color=GREEN, stroke_width=2)

        self.play(FadeIn(axes), FadeIn(axes_labels))
        self.play(FadeIn(guide_line))
        self.play(LaggedStart(*(FadeIn(dot) for dot in dots), lag_ratio=0.05))
        self.play(FadeIn(corr_text))
        self.wait(2.0)
