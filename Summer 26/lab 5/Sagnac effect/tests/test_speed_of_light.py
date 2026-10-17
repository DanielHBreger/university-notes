"""Speed of light from the slope: Eq. (5) solved for c must invert the theory slope."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sagnac import C_LIGHT, N_TURNS, RADIUS, U_N_TURNS, U_RADIUS, speed_of_light, theory  # noqa: E402


def test_theory_slope_gives_defined_c():
    a_th, _ = theory()
    c, comps = speed_of_light(a_th, {'fit': 0.01 * a_th})
    assert np.isclose(c, C_LIGHT, rtol=1e-12)
    assert np.isclose(comps['fit'], 0.01 * C_LIGHT)                   # c ~ 1/a: same relative size
    assert np.isclose(comps['N'], C_LIGHT * U_N_TURNS / N_TURNS)
    assert np.isclose(comps['R'], 2 * C_LIGHT * U_RADIUS / RADIUS)    # A = pi R^2


def test_c_inverse_in_slope():
    c1, _ = speed_of_light(1.0, {})
    c2, _ = speed_of_light(2.0, {})
    assert np.isclose(c1, 2 * c2)
