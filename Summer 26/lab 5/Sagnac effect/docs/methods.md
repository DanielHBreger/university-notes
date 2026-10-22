# `sagnac_update/methods.py`

> **Note (22 October).** This describes `sagnac_update/methods.py`, the presentation stage. The report uses the top-level `methods.py`, whose template matching is identical; its robustness scans now use `sagnac.regress` (equal-weight regression, third value s instead of σ_x). `sagnac_update/methods.py` imports the old `sagnac_update/sagnac.py`, which needs the earlier `calibration.py` (`previous_reports/calibration_eps_2026-10-13.py`) to run.

Two more ways to read the phase of every photo, as checks of the ring minima: template matching of the r² profile, and ten random pixels fixed in the image. The per-photo phases go to `methods_frames.csv`, where `sagnac.py` fits them like the minima ([sagnac.md](sagnac.md)); the robustness scans go to `methods.json` and `methods.txt`. Numbers below are from `sagnac_update/methods.txt` unless another file is named.

Inputs, from the script's own folder: `fringes_frames.csv` and `fringes_profiles.npz` (`fringes.py`: per photo the ring center, P, and the intensity averaged over the fringe sector in 130 bins of s = r², 2000 px² wide, up to 2.6 × 10⁵ px²), `calibration.json` (scan slopes only) and, for the pixels, the 161 photos.

## The idea

All methods read the same photos, so they are checks, not independent measurements; they make different assumptions, and a real Sagnac slope must come out of all of them. The minima and the cosine fit read the phase at the ring center. Template matching reads how far the whole profile moved, a phase averaged over its comparison window. A fixed pixel reads the local phase at its own radius. A change of the fringe period P between photos therefore enters both new methods, and a move of the ring center enters the pixels; neither enters the minima.

## `main()`, step 1

```python
assert list(z['files']) == list(fr.file), 'profiles and frame table differ: rerun fringes.py'
```

The profiles and the frame table must hold the same photos in the same order. `rest` marks S = 0; `model` is `calibration.json['conversion']`.

## `template_phases(s, M, rest, match=MATCH_S)`

```python
T = interpolate.CubicSpline(s, M[rest].mean(axis=0))
sel = (s >= SELF_S[0]) & (s <= SELF_S[1])
P_T = best_shift(s[sel], T(s[sel]), T, 0.8 * P_GUESS, 1.25 * P_GUESS)
```

The template T is the mean of the 10 rest profiles, made a smooth function of s by a cubic spline so that it can be shifted by any amount. Each profile is matched to g T(s + D) + c over `match` (2 × 10⁴ to 1.8 × 10⁵ px²). Rings that moved inward make the profile look like the rest profile further out, so D > 0, and the phase relative to rest is

```
Δφ = 2π D / P_T
```

P_T is the template's own period, the shift that maps T onto itself one fringe further out, searched between 0.8 and 1.25 `P_GUESS` over `SELF_S`: 108.02k px², against 107.37k for the median of the minima fits. No fringe shape is assumed.

## `match_residual(D, s, prof, T)` and `best_shift(s, prof, T, lo, hi)`

```python
A = np.c_[T(s + D), np.ones_like(s)]               # columns: shifted template, constant
gc, *_ = np.linalg.lstsq(A, prof, rcond=None)
```

For each trial D the gain g and offset c are solved exactly by linear least squares, so they need no start values. They absorb changes of brightness; in s, a single exponential envelope shifted by D only changes by a constant factor, which g absorbs too. `best_shift` scans D every 1000 px² over [lo, hi] and refines with `minimize_scalar` between the grid points either side of the best one. `SHIFT_RANGE` allows −2 × 10⁴ to 7.8 × 10⁴ px², so s + D stays inside the profile; the measured D run from −1.6k to 55.4k px².

The phases are relative to the rest template but not set to zero at rest (rest mean −0.0008 rad, `sagnac_speeds.csv`); the free intercept of the line in `sagnac.py` makes that irrelevant.

## Choosing the pixels: `fringe_region`, `draw_points`, `dark_level`

The rest geometry is the median of the rest photos' fits: center (897.8, 507.7) px, P0 = 107.55k px². `rest_img` is the mean of the smoothed rest images.

- `fringe_region`: every pixel between −75° and 180° (`THETA_DEG`) with r² < 2.6 × 10⁵ px² (`S2_FIT_MAX`) about that center.
- `draw_points`: `POINTS_N` = 10 different pixels, each equally likely (`default_rng(seed).choice(..., replace=False)`). Inside the sector equal steps of r² enclose equal areas, so the draw is even in r², that is, in fringe phase. The seed, 20261011 (the measurement date), was fixed in advance.
- `dark_level`: the median of the four 100 × 100 px corners of `rest_img`, far outside the beam: 2.0 counts.

All 31 draws are read in one pass, from the same Gaussian-smoothed red channel (σ = 4 px) as the other methods, at the same image positions in every photo.

## `point_fit(M, psi, dark)`

```
I_kj = a_j + b_j cos(ψ_j + δ_k),     ψ_j = 2π s_j / P0
```

Pixel j, photo k. ψ_j is the pixel's fringe phase at rest, from its r² about the rest center; δ_k is the photo's phase (the convention of the center phase in `fringes.py`); a_j and b_j are the pixel's offset and amplitude. Alternating least squares fits the 161 phases and 20 pixel parameters, starting from a_j = the pixel's mean and b_j = half its range:

1. a, b fixed: every δ_k (`best_phases`);
2. δ fixed: every a_j, b_j (`offsets_amplitudes`);

until no δ_k moves by more than 10⁻⁶ rad (`POINTS_TOL`), at most 500 rounds (`POINTS_MAX_ITER`; the round count is returned, so a fit that hit the limit can be spotted). `best_phases` computes each photo's misfit on 721 phases from −π to π (half a degree apart), all photos at once, and refines the best with a parabola through it and its two neighbors. The main draw: 116 rounds, rms residual 10.19 counts, amplitudes 21.0–139.6 counts, phases spanning 3.23 rad.

## `offsets_amplitudes(M, X, dark)`: the two limits

```
b_j ≥ 0,     a_j − b_j ≥ dark
```

a_j − b_j is the darkest the pixel's fringe gets, and it cannot be darker than the camera's black. A negative b_j would be the same as moving the pixel half a fringe (ψ_j + π), undoing the geometry. The limits are needed because the largest shift is only about half a fringe (3.3 rad, RESULTS.md): each pixel follows a piece of a cosine, and ten such pieces cannot tell a large phase change at low contrast from a small one at high contrast. The dark-level limit caps the contrast, and with it how small the phase change can be made.

With c = a − b the model is M = c + b(1 + X), X = cos(ψ_j + δ_k): a two-parameter least squares with c ≥ dark and b ≥ 0. Its best point is the unconstrained one if that is allowed, otherwise the best point on one of the two edges:

```python
cands = [(c_u, b_u, (c_u >= dark) & (b_u >= 0)),                                  # interior
         (np.full_like(Sy, dark), np.maximum(0.0, (Szy - dark * Sz) / Szz), True),   # on c = dark
         (np.maximum(dark, Sy / n), np.zeros_like(Sy), True)]                     # on b = 0
```

Per pixel the allowed candidate with the smallest sum of squares wins. This is the exact bounded solution; `sagnac_update/tests/test_methods.py` compares it with `scipy.optimize.lsq_linear`.

## `point_phases` and `relative_to_rest(fr, phase)`

```python
order = np.lexsort((fr.time_s, fr.S))              # sort by S, ties by time
du = np.empty_like(phase); du[order] = np.unwrap(phase[order])
return du - du[(fr.S == 0).to_numpy()].mean()
```

The photos are ordered by S, then time; `np.unwrap` removes 2π jumps between neighbors (it assumes they differ by less than π), and the rest mean is subtracted. No sign needs fixing, since ψ ties δ to the ring geometry.

## Step 4: `methods_frames.csv` and `illustration(...)`

`methods_frames.csv`: `file`, `phi_template_rad` (rad, 5 decimals), `shift_template_px2` (D in px², 1 decimal), `phi_points_rad` (rad, 5 decimals).

`figures/methods_illustration` (.pdf and .png, full text width), 2 × 2: (a) the cosine fit and (b) the template match of one S = 90 photo, the one with the median D there, with the comparison range shaded; (c) the ten pixels, numbered, on the mean rest image (±540 px about the rest center, marked +); (d) their intensities against the fitted phase relative to rest, each with its fitted cosine over the range of the data. This version is in `sagnac_update/figures/`; the copy in the top-level `figures/` is the earlier six-panel figure.

## Step 5: robustness scans, with `slope(fr, phi, model)`

`slope` puts any per-photo phases through `sagnac.py`'s `setting_means`, `calibration_part` and `fit_line`, and returns a, its fit error (no calibration-curve part) and σ_x.

- Template, the four ranges of `MATCH_SCAN`: 1.768 ± 0.017 (the main range), 1.770 ± 0.015, 1.764 ± 0.015, and 1.834 ± 0.035 for 5 × 10⁴ to 1.8 × 10⁵ px².
- Pixels, the main draw and seeds 1–30, ψ from the rest geometry: main draw 1.786 ± 0.037 (σ_x = 0.068 rad); all 31: median 1.764, 16–84 % 1.506–1.898, full range 1.414–2.119; one fit stopped at 500 rounds.
- Diagnostic: the same pixels with ψ from each photo's own fitted center and period (`psi_frame`), which removes P changes and center motion: median 1.779, 16–84 % 1.681–1.978.

`methods.json`: `calibration_model`, `calibration_degree`; `template_scan` (`range`, `slope`, `u_fit`, `sigma_x`); `points` (`seed`, `n`, `dark`, `centre`, `P0`, `pixels_xy`, `offsets`, `amplitudes`, `rms`, `iterations`); `points_scan` and `points_frame_geometry_scan`, 31 rows each with the main draw first (`seed`, `slope`, `u_fit`, `sigma_x`, `rms`, `iterations`). `methods.txt` says the same in words, with the pixel coordinates.

## Constants

| name | value |
|---|---|
| `MATCH_S` | (2e4, 1.8e5) px², template comparison range |
| `SHIFT_RANGE` | (-2e4, 7.8e4) px², allowed D |
| `SELF_S` | (1e3, 1.1e5) px², range for P_T |
| `MATCH_SCAN` | (2e4, 1.8e5), (1e4, 1.2e5), (2e4, 1.4e5), (5e4, 1.8e5) px² |
| `POINTS_N`, `POINTS_SEED`, `POINTS_SCAN_SEEDS` | 10, 20261011, range(1, 31) |
| `POINTS_TOL`, `POINTS_MAX_ITER` | 1e-6 rad, 500 |
| `PHASE_GRID` | np.linspace(-np.pi, np.pi, 721) |
| from `fringes.py` | `P_GUESS` 107.5e3 px², `S2_FIT_MAX` 2.6e5 px², `SMOOTH_PX` 4.0 px, `THETA_DEG` (-75.0, 180.0) deg |

## How the top-level version differs

Template matching is the same code with the same constants. Instead of the pixels there are two other methods:

- **AIA** (phase shifting with unknown steps, Wang & Han 2004): alternating least squares on the r² profile bins, each bin with its own sinusoid a_j + c_j cos δ_k + d_j sin δ_k, so no geometry is imposed. It starts from the first two principal components, picks the sign so that δ grows with S, stops at 10⁻¹⁰ rad or 3000 rounds, and is scanned over the largest s used (`AIA_SCAN`), with and without a brightness factor per photo.
- **2D Fourier sideband** (Takeda): each image resampled on a 256 × 256 (θ, r²) grid about its own ring center, where the rings are straight stripes; background removed, a Gaussian filter of width 0.26 × the carrier frequency (`FOURIER_SCAN` 0.15, 0.26, 0.40), phase arg Σ Z_k conj(Z_rest).

`methods_frames.csv` has `phi_aia_rad` and `phi_fourier_rad` instead of `phi_points_rad`; `methods.json` has `template_scan`, `aia_scan`, `fourier_scan`; the figure has six panels. It ran with the quadratic conversion, so its template slopes differ (1.798 for the main range, top-level `methods.txt`), while the template phases are identical. AIA's slope ranged over 1.61–2.11 with window and variant, which is why the earlier `sagnac.py` left it out of the method spread.

## Caveats

- `sagnac_update/` holds none of these inputs (they are in the parent folder), so `main()` stops at once where it sits. Template matching needs only the two `fringes.py` outputs: rerun on the top-level copies it reproduces `phi_template_rad`, `shift_template_px2` and, with the linear conversion, the template scan exactly. The pixel part needs all 161 photos and was not rerun for this check; `random_points_scan.py`, which calls the same functions, reproduces its 31 slopes to 10⁻⁷ ([random_points_scan.md](random_points_scan.md)).
- The calibration enters only the scan slopes, not the per-photo phases.
- Neither method reads the phase at the ring center: the template slope moves with the comparison range (1.834 against 1.768), and the pixel slope with the geometry used for ψ (diagnostic median 1.779 against 1.764).
- The pixel slope depends on the draw far more than one draw's fit error says (16–84 % range 1.506–1.898, fit error 0.037). The presentation's 1000-draw scan with these functions gives SD 0.196, 3.4 times the typical fit error (RESULTS.md, [random_points_scan.md](random_points_scan.md)). That is why `sagnac.py` keeps the pixels out of `u_method`.
