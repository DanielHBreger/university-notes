"""Synthetic checks of the phase methods in methods.py: profiles with a known phase sweep."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import methods  # noqa: E402

S = np.arange(1000, 260000, 2000.0)      # the r^2 bins of fringes.py
P = 107500.0


def profiles(delta, gain=1.0):
    """Two-beam ring profiles in r^2 with phase delta per frame (rings move inwards as delta grows)."""
    g = np.broadcast_to(gain, np.shape(delta))
    return np.array([gk * (80 * np.exp(-S / 2.6e5) + 85 * np.exp(-S / 1.4e5) * np.cos(2 * np.pi * S / P + d) + 7)
                     for d, gk in zip(delta, g)])


def test_aia_recovers_half_fringe_sweep():
    rng = np.random.default_rng(0)
    delta = np.sort(rng.uniform(-0.75, 2.55, 161))                # the measured range, 190 deg
    M = profiles(delta) + rng.normal(0, 1.5, (161, len(S)))
    d = np.unwrap(methods.aia(M)[0])
    d = d if np.polyfit(delta, d, 1)[0] > 0 else -d
    assert np.abs((d - d.mean()) - (delta - delta.mean())).max() < 0.03


def test_fourier_recovers_sweep():
    rng = np.random.default_rng(2)
    delta = np.r_[np.zeros(10), np.linspace(0.3, 3.3, 11)]
    s = np.linspace(0, methods.S2_FIT_MAX, methods.FOURIER_GRID[1], endpoint=False)
    prof = [80 * np.exp(-s / 2.6e5) + 85 * np.exp(-s / 1.4e5) * np.cos(2 * np.pi * s / P + d) + 7 for d in delta]
    maps = np.repeat(np.array(prof)[:, None, :], methods.FOURIER_GRID[0], axis=1)    # straight in (theta, s)
    maps += rng.normal(0, 1.5, maps.shape)
    phi = np.unwrap(methods.fourier_phases(maps, delta == 0)[0])
    assert np.abs(phi - delta).max() < 0.1


def test_template_recovers_shift():
    rng = np.random.default_rng(1)
    delta = np.r_[np.zeros(10), np.linspace(0.3, 3.3, 11)]
    M = profiles(delta) + rng.normal(0, 1.5, (len(delta), len(S)))
    _, P_T, _, phi = methods.template_phases(S, M, delta == 0)
    assert abs(P_T / P - 1) < 0.02
    assert np.abs(phi - delta).max() < 0.1
