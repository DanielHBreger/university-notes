"""Synthetic check of template matching (methods.py): profiles with a known phase sweep."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import methods  # noqa: E402

S_BINS = np.arange(1000, 260000, 2000.0)      # the r^2 bins of fringes.py
PERIOD = 107500.0


def make_profiles(phases):
    """Two-beam ring profiles in r^2, one per phase (the rings move inwards as the phase grows)."""
    return np.array([80 * np.exp(-S_BINS / 2.6e5) + 85 * np.exp(-S_BINS / 1.4e5) * np.cos(2 * np.pi * S_BINS / PERIOD + phi) + 7
                     for phi in phases])


def test_template_recovers_shift():
    rng = np.random.default_rng(1)
    phases = np.r_[np.zeros(10), np.linspace(0.3, 3.3, 11)]
    profiles = make_profiles(phases) + rng.normal(0, 1.5, (len(phases), len(S_BINS)))
    _, period, _, found = methods.template_phases(S_BINS, profiles, phases == 0)
    assert abs(period / PERIOD - 1) < 0.02
    assert np.abs(found - phases).max() < 0.1
