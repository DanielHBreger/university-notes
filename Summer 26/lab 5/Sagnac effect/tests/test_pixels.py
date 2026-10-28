"""Synthetic check of the fixed-pixel readout (pixels.py): intensities with a known phase sweep."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pixels  # noqa: E402


def test_fit_recovers_phases():
    rng = np.random.default_rng(3)
    phases = np.r_[np.zeros(10), np.linspace(0.3, 3.3, 11)]           # rest, then a sweep of half a fringe
    S = np.r_[np.zeros(10), np.arange(1, 12)]
    psi = rng.uniform(0, 2 * np.pi, 10)                                  # ten pixels at random fringe phases
    a, b = rng.uniform(60, 120, 10), rng.uniform(20, 50, 10)
    dark = 5.0
    I = a + b * np.cos(psi + phases[:, None]) + rng.normal(0, 1.0, (len(phases), 10))
    delta, a_fit, b_fit, rms, rounds = pixels.fit_phases(I, psi, dark)
    found = pixels.relative_to_rest(S, np.arange(len(S)), delta)
    assert rounds < pixels.MAX_ROUNDS
    assert np.abs(found - phases).max() < 0.05
    assert np.abs(a_fit - a).max() < 2 and np.abs(b_fit - b).max() < 2
