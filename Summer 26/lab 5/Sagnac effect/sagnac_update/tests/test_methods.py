"""Synthetic checks of the phase methods in methods.py: profiles and pixels with a known phase sweep."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import methods  # noqa: E402

S = np.arange(1000, 260000, 2000.0)      # the r^2 bins of fringes.py
P = 107500.0


def profiles(delta, gain=1.0):
    """Two-beam ring profiles in r^2 with phase delta per frame (rings move inwards as delta grows)."""
    g = np.broadcast_to(gain, np.shape(delta))
    return np.array([gk * (80 * np.exp(-S / 2.6e5) + 85 * np.exp(-S / 1.4e5) * np.cos(2 * np.pi * S / P + d) + 7)
                     for d, gk in zip(delta, g)])


def ring_intensity(ys, xs, delta, centre=(450.0, 380.0), dark=2.0):
    """Two-beam ring pattern at pixels (ys, xs) for each phase in delta: photographs x pixels."""
    s = (xs - centre[0]) ** 2 + (ys - centre[1]) ** 2
    return np.array([dark + 80 * np.exp(-s / 2.6e5) + 70 * np.exp(-s / 1.6e5) * (1 + np.cos(2 * np.pi * s / P + d))
                     for d in delta])


def test_points_drawn_reproducibly_inside_region():
    ys, xs = methods.fringe_region((700, 900), 450.0, 380.0)
    a, b = methods.draw_points(ys, xs, 7), methods.draw_points(ys, xs, 7)
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])
    y, x = a
    assert len(set(zip(y, x))) == methods.POINTS_N
    th = np.degrees(np.arctan2(-(y - 380.0), x - 450.0))
    assert np.all((th > methods.THETA_DEG[0]) & (th < methods.THETA_DEG[1]))
    assert np.all((x - 450.0) ** 2 + (y - 380.0) ** 2 < methods.S2_FIT_MAX)


def test_offsets_amplitudes_is_the_bounded_least_squares():
    rng = np.random.default_rng(3)
    K, J, dark = 40, 6, 5.0
    X = np.cos(rng.uniform(-np.pi, np.pi, (K, J)) * rng.uniform(0.1, 1))
    M = rng.uniform(0, 200, J) + rng.uniform(-80, 120, J) * X + rng.normal(0, 20, (K, J))
    a, b = methods.offsets_amplitudes(M, X, dark)
    assert np.all(b >= 0) and np.all(a - b >= dark - 1e-9)
    for j in range(J):
        ref = optimize.lsq_linear(np.c_[np.ones(K), 1 + X[:, j]], M[:, j],
                                  bounds=([dark, 0], [np.inf, np.inf]), tol=1e-12).x
        sse_ref = np.sum((M[:, j] - ref[0] - ref[1] * (1 + X[:, j])) ** 2)
        assert np.sum((M[:, j] - a[j] - b[j] * X[:, j]) ** 2) <= sse_ref * (1 + 1e-9)


def test_points_recover_half_fringe_sweep():
    rng = np.random.default_rng(4)
    delta = np.r_[np.full(10, -0.75), np.linspace(-0.4, 2.5, 60)]          # the measured range, 190 deg
    fr = pd.DataFrame(dict(S=np.r_[np.zeros(10), np.repeat(np.arange(1.0, 7.0), 10)], time_s=np.arange(70)))
    ys, xs = methods.fringe_region((700, 900), 450.0, 380.0)
    for seed in range(5):
        y, x = methods.draw_points(ys, xs, seed)
        V = ring_intensity(y, x, delta) + rng.normal(0, 2.0, (len(delta), len(x)))
        psi = 2 * np.pi * ((x - 450.0) ** 2 + (y - 380.0) ** 2) / P
        phi = methods.point_phases(fr, V, psi, 2.0)[0]
        truth = delta - delta[:10].mean()
        assert np.abs(phi - truth).max() < 0.15
        assert abs(np.polyfit(truth, phi, 1)[0] - 1) < 0.03


def test_template_recovers_shift():
    rng = np.random.default_rng(1)
    delta = np.r_[np.zeros(10), np.linspace(0.3, 3.3, 11)]
    M = profiles(delta) + rng.normal(0, 1.5, (len(delta), len(S)))
    _, P_T, _, phi = methods.template_phases(S, M, delta == 0)
    assert abs(P_T / P - 1) < 0.02
    assert np.abs(phi - delta).max() < 0.1
