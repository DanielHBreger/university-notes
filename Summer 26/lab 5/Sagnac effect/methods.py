"""Template matching: a third way to read the fringe phase, with no fringe model.

fringes.py reads the phase of every photograph from the dark-ring radii and from a cosine fit.
Here the at-rest profile (intensity against s = r^2, averaged over the rest photographs) is the
template T(s). Each photograph's profile is matched to g T(s + D) + c: the template is slid
along r^2 by D, with a free gain g and offset c that absorb brightness changes (in s, an
exponential envelope shifted by D only changes by a constant factor). Rings that move inwards
give D > 0, and the phase relative to rest is 2 pi D / P_T, where P_T is the template's own
period: the shift that maps it onto itself one fringe further out. In words: slide the at-rest
profile along r^2 until it lies on top of the new one.

Reads fringes_frames.csv and fringes_profiles.npz (fringes.py) and writes methods_frames.csv,
the phase of every photograph, which sagnac.py fits like the other two methods.
Run from this folder after fringes.py:   python methods.py   (about a minute)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import interpolate, optimize

from fringes import P_GUESS

HERE = Path(__file__).resolve().parent
MATCH_RANGE = (2e4, 1.8e5)      # px^2 of s over which a profile is compared with the shifted template
SHIFT_RANGE = (-2e4, 7.8e4)     # px^2, allowed shifts D: s + D stays inside the profile
SELF_RANGE = (1e3, 1.1e5)       # px^2, comparison range for the template's own period


def misfit(D, s, profile, template):
    """Sum of squared differences between profile(s) and g template(s + D) + c, with the gain g and
    offset c solved by linear least squares for this shift D. Returns the sum and (g, c)."""
    columns = np.c_[template(s + D), np.ones_like(s)]
    gain_offset, *_ = np.linalg.lstsq(columns, profile, rcond=None)
    return float(np.sum((profile - columns @ gain_offset) ** 2)), gain_offset


def best_shift(s, profile, template, lo, hi):
    """The shift D in [lo, hi] (px^2) with the smallest misfit: first a grid every 1000 px^2, then
    a 1D minimiser between the grid points either side of the best one."""
    grid = np.arange(lo, hi + 1, 1000.0)
    best = int(np.argmin([misfit(D, s, profile, template)[0] for D in grid]))
    bracket = (grid[max(best - 1, 0)], grid[min(best + 1, len(grid) - 1)])
    return optimize.minimize_scalar(lambda D: misfit(D, s, profile, template)[0],
                                    bounds=bracket, method='bounded').x


def template_phases(s, profiles, rest, match_range=MATCH_RANGE):
    """Template matching of every profile (rows of `profiles`; `rest` marks the rest photographs).
    The template is the mean rest profile as a cubic spline, so it can be shifted by any D.
    Returns the template, its period P_T, the shifts D and the phases 2 pi D / P_T."""
    template = interpolate.CubicSpline(s, profiles[rest].mean(axis=0))
    inside = (s >= SELF_RANGE[0]) & (s <= SELF_RANGE[1])
    P_T = best_shift(s[inside], template(s[inside]), template, 0.8 * P_GUESS, 1.25 * P_GUESS)
    inside = (s >= match_range[0]) & (s <= match_range[1])
    shifts = np.array([best_shift(s[inside], profile[inside], template, *SHIFT_RANGE) for profile in profiles])
    return template, P_T, shifts, 2 * np.pi * shifts / P_T


def main():
    photos = pd.read_csv(HERE / 'fringes_frames.csv')
    saved = np.load(HERE / 'fringes_profiles.npz')
    s, profiles = saved['s'], saved['profiles']            # s: bin centres (px^2); profiles: photographs x bins
    assert list(saved['files']) == list(photos.file), 'profiles and frame table differ: rerun fringes.py'
    rest = (photos.S == 0).to_numpy()

    template, P_T, shifts, phases = template_phases(s, profiles, rest)
    pd.DataFrame(dict(file=photos.file, phi_template_rad=np.round(phases, 5),
                      shift_template_px2=np.round(shifts, 1))).to_csv(HERE / 'methods_frames.csv', index=False)
    print(f'template: mean of {rest.sum()} rest profiles, own period P_T = {P_T / 1e3:.2f}k px^2 '
          f'(ring minima: median P = {photos.P_px2.median() / 1e3:.2f}k); shifts D = {shifts.min() / 1e3:.1f} to '
          f'{shifts.max() / 1e3:.1f}k px^2; {len(photos)} phases -> methods_frames.csv')


if __name__ == '__main__':
    main()
