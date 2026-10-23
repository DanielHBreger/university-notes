# `methods.py`

Template matching: the third way the report reads the phase of every photograph, as a check of the ring minima and the cosine fit of `fringes.py`. The per-photo phases go to `methods_frames.csv`, where `sagnac.py` fits them like the other two ([sagnac.md](sagnac.md)). About a minute.

Inputs, from the script's own folder: `fringes_frames.csv` and `fringes_profiles.npz` (`fringes.py`: per photo the ring center, P, and the intensity averaged over the fringe sector in 130 bins of s = r², 2000 px² wide, up to 2.6 × 10⁵ px²).

## The idea

All methods read the same photos, so they are checks, not independent measurements; they make different assumptions, and a real Sagnac slope must come out of all of them. The minima and the cosine fit read the phase at the ring center. Template matching reads how far the whole profile moved along r², a phase averaged over its comparison window, with no fringe model at all. A change of the fringe period P between photos therefore enters it, but not the minima.

## `template_phases(s, profiles, rest, match_range=MATCH_RANGE)`

```python
template = interpolate.CubicSpline(s, profiles[rest].mean(axis=0))
inside = (s >= SELF_RANGE[0]) & (s <= SELF_RANGE[1])
P_T = best_shift(s[inside], template(s[inside]), template, 0.8 * P_GUESS, 1.25 * P_GUESS)
```

The template is the mean of the 10 rest profiles, made a smooth function of s by a cubic spline so that it can be shifted by any amount. Each profile is matched to g T(s + D) + c over `match_range` (2 × 10⁴ to 1.8 × 10⁵ px²). Rings that moved inward make the profile look like the rest profile further out, so D > 0, and the phase relative to rest is

```
Δφ = 2π D / P_T
```

P_T is the template's own period, the shift that maps T onto itself one fringe further out, searched between 0.8 and 1.25 `P_GUESS` over `SELF_RANGE`: 108.02k px², against 107.37k for the median of the minima fits.

## `misfit(D, s, profile, template)` and `best_shift(s, profile, template, lo, hi)`

```python
columns = np.c_[template(s + D), np.ones_like(s)]
gain_offset, *_ = np.linalg.lstsq(columns, profile, rcond=None)
```

For each trial D the gain g and offset c are solved exactly by linear least squares, so they need no start values. They absorb changes of brightness; in s, a single exponential envelope shifted by D only changes by a constant factor, which g absorbs too. `best_shift` scans D every 1000 px² over [lo, hi] and refines with `minimize_scalar` between the grid points either side of the best one. `SHIFT_RANGE` allows −2 × 10⁴ to 7.8 × 10⁴ px², so s + D stays inside the profile; the measured D run from −1.6k to 55.4k px².

The phases are relative to the rest template but not set to zero at rest (rest mean −0.0008 rad, `sagnac_speeds.csv`); the free intercept of the line in `sagnac.py` makes that irrelevant.

## `main()`

```python
assert list(saved['files']) == list(photos.file), 'profiles and frame table differ: rerun fringes.py'
```

The profiles and the frame table must hold the same photos in the same order. `rest` marks S = 0. The output `methods_frames.csv` has `file`, `phi_template_rad` (rad, 5 decimals) and `shift_template_px2` (D in px², 1 decimal); one line is printed with P_T and the range of D.

## Constants

| name | value |
|---|---|
| `MATCH_RANGE` | (2e4, 1.8e5) px², template comparison range |
| `SHIFT_RANGE` | (-2e4, 7.8e4) px², allowed D |
| `SELF_RANGE` | (1e3, 1.1e5) px², range for P_T |
| from `fringes.py` | `P_GUESS` 107.5e3 px² |

## Caveats

- Template matching does not read the phase at the ring center: its slope moves with the comparison range (1.834 for 5 × 10⁴ to 1.8 × 10⁵ px² against 1.768 for the main range, in the scans of 22 October). The main range is the one that uses the three dark rings.
- `tests/test_methods.py` checks it on synthetic profiles: P_T within 2 % and every phase within 0.1 rad of a known 3.3 rad sweep.
- Earlier versions (`previous_reports/methods_2026-10-13.py`, `sagnac_update/methods.py`) added further methods that the final analysis never used: ten fixed random pixels (the presentation stage; its 1000-draw scan is still read as stored by `sagnac.py`, prep Q5), phase shifting with unknown steps (AIA) and a 2D Fourier sideband, each with robustness scans and a six-panel illustration (`figures/methods_illustration`). They were removed on 23 October to keep the pipeline to what the report uses.
