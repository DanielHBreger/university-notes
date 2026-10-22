# `fringes.py`

Reduces each of the 161 photographs to its ring center, the radii of the innermost dark rings, the ring spacing P and, above all, φ_c, the fringe phase at the ring center, read two ways: from the dark-ring radii (main method) and from a cosine fit of the intensity profile (second method, prep Q5). Every later phase analysis starts from its two outputs. About a second per photograph.

## The idea

**Why rings, evenly spaced in r².** The two beams leave the fiber as Gaussian beams, which carry the transverse phase k r²/(2R(z)). If they reach the camera with different wavefront radii R1 ≠ R2, their phase difference at distance r from the ring center is

```
phi(r) = (pi/lambda) (1/R1 - 1/R2) r^2 + phi_c  =  2 pi r^2 / P + phi_c,      P = 2 lambda / |1/R1 - 1/R2|
```

Going out by P in s = r² adds one fringe, the same P for every ring, so in r the rings crowd together outward. In the S = 90 photograph `20261011163338` the dark rings are at 102.3, 343.5 and 475.3 px (`fringes_frames.csv`): 241 and then 132 px apart in r, but 1.08 × 10⁵ px² both times in r². Rotation adds the same phase everywhere, so it changes φ_c and the rings slide in or out; a tilt between the beams mainly moves the center, refocusing mainly changes P. Hence the center and P are fitted in every photograph, and φ_c is tracked.

**From ring radii to φ_c.** A dark ring is where φ = π + 2πm, so

```
R_m^2 = P (m + eta),      phi_c = pi - 2 pi eta
```

and a straight line through the points (m, R_m²) gives P (slope) and Pη (intercept). This is the manual's Eq. (7), the shift of a minimum over the local fringe period, done in r², where every ring has the same period. η is a ratio of two numbers in px², so no pixel-to-length calibration is needed.

**Sign convention.** The model takes the phase to grow outward (+2πr²/P, beams labeled so that 1/R1 − 1/R2 > 0). The intensity depends on cos φ, which cannot tell φ from −φ, so this is a choice; the opposite one would flip the sign of every φ_c. With it, contracting rings (η falling) mean φ_c rising: η is 0.62 at rest and 0.10 at the top speed (RESULTS.md), and the setting mean of φ_c goes from −0.75 to 2.54 rad (`phi` in `sagnac_speeds.csv`). The rotation direction was not recorded, so only the size of the shift is compared.

**Why minima.** A minimum does not move under any monotonic camera response (gamma), so the main method needs no linear camera. The docstring says the same of the brightness envelope; strictly, an envelope that changes across a fringe pulls a dip slightly toward its darker side. The cosine fit uses every bin, but depends on its envelope model and a linear camera.

## `analyse_frame(path, start=START)`

The driver for one photograph. Step 1 smooths the red channel once (`ndimage.gaussian_filter`, σ = 4 px) against speckle; all later steps use the smoothed image. Steps 2–4, `ring_guesses(sector_profile(...))`, `ray_minima` and `fit_concentric`, run three times (`N_ITER`), each round centered on the previous round's fitted center (the first on `START`); the saved center, radii and errors are the third round's. The three innermost rings are kept, then `ring_phase` (step 5) and `cosine_fit` (step 6) run about the final center.

## `load_red(path)` and `frame_time(path)`

`load_red` returns channel 0 (red) of the RGBA PNG as float counts; the other channels are dark, and red never saturates (maximum 239 of 255, RESULTS.md). `frame_time` turns characters 8–13 of the file name `YYYYMMDDhhmmss` into seconds since midnight.

## `sector_profile(img, cx, cy)`

Mean intensity in 1 px wide annuli about (cx, cy), r = 0 … 640 px, over the pixels in the fringe sector only:

```python
th = np.degrees(np.arctan2(-(y - cy), x - cx))      # minus: image y grows downwards
```

Angles run anticlockwise from +x with "up" positive, hence the minus sign. The sector −75° to 180° leaves out the lower left, where the photographs show no fringes.

## `ring_guesses(profile)`

The profile is smoothed lightly (σ = 3 px along r) and its dips are found as peaks of its negative (`find_peaks`), keeping those whose prominence exceeds 3 % of the profile's maximum. Dips within 40 px of the center are ignored, since an annulus there holds few pixels. Each dip is a first guess of a dark-ring radius, in whole pixels.

## `ray_minima(img, cx, cy, guesses)`

Along 170 rays, every 1.5° from −75° to 178.5°, the image is sampled at r = 0, 1, …, 640 px by bilinear interpolation (`map_coordinates`, `order=1`; past the image edge `mode='nearest'` repeats the edge pixel). The rings are evenly spaced in r², so the local period in r at radius g is about P/(2g), from d(r²) = 2r dr. For each guess g the window is ±0.3 of that period (`lam = P_GUESS / (2 * g)`), at most ±0.8 g so it never reaches the center. The lowest sample q in the window is kept only if it is a real dip:

```python
if q <= lo or q >= hi or min(seg[0], seg[-1]) - p[q] < MIN_DEPTH:
    continue
```

A minimum on the window's edge means the intensity is still falling there (a flat stretch past the image edge fails here too, since its lowest sample is the first); a dip less than 3 counts below the lower window end is noise. The kept minimum is refined below a pixel by a parabola through the samples within h = max(2, int(lam/10)) px, a tenth of a period on each side:

```python
r = -c1 / (2 * c2) if c2 > 0 else _RR[q]    # vertex of c2 r^2 + c1 r + c0
if abs(r - _RR[q]) > h:                     # vertex outside the fitted points: keep the sample
    r = _RR[q]
```

Only an upward-opening parabola has a minimum, and a vertex outside the fitted points would be an extrapolation; in both cases the sample is used. Each point is stored as (ring index, x, y), at most one per ray and ring.

## `fit_concentric(pts, cx, cy, guesses)`

One common center and one radius per ring, fitted to all ray minima at once. A point's residual is its distance from the center minus its ring's radius:

```python
sol = optimize.least_squares(res, p0, loss='soft_l1', f_scale=4.0)
```

Rings with fewer than 30 minima (`MIN_POINTS`) are dropped first. The `soft_l1` loss acts like least squares for residuals well under `f_scale` = 4 px and grows only linearly beyond, so a few stray minima (speckle) barely pull the circles. The covariance is `pinv(JᵀJ)` times the sum of squared residuals per degree of freedom; scipy's J is already down-weighted for stray points, while the sum uses the raw residuals. Its square root is also returned as the rms distance of the minima from their circles.

## `ring_phase(R, P_fallback=P_GUESS)`

```python
m = np.r_[0, np.cumsum(np.rint(np.diff(s) / P_GUESS))]     # ring order, gaps allowed
P, s0 = np.polyfit(m, s, 1)                                 # s0: s of ring order 0
return P, float(np.angle(np.exp(1j * (np.pi - 2 * np.pi * s0 / P))))
```

With s the sorted R², the innermost ring found gets m = 0 and the others are numbered by their gap in r² over `P_GUESS`, rounded, so a missing ring does not shift the numbering. Then η = s0/P, and φ_c = π − 2πη is wrapped into (−π, π] by `angle(exp(i x))`. Two rings fix the line exactly and three leave one degree of freedom, so the line has no fit statistic; its check is the circle fit. A single ring would get P = `P_GUESS`; every photograph had two or three.

## `sector_profile_s2`, `cos_model` and `cosine_fit(img, cx, cy)`

The second method. `sector_profile_s2(img, cx, cy, edges)` is the same sector average in bins of s = r²: equal steps in r² enclose equal areas, so every bin holds about the same number of pixels, and the fringes are evenly spaced along this axis. `cosine_fit` uses 130 bins of 2000 px² up to `S2_FIT_MAX` = 2.6 × 10⁵ px² (r up to 510 px) and fits `cos_model`:

```
I(s) = C exp(-s/sC) + B exp(-s/sB) cos(2 pi s/P + phi) + D
```

What this assumes: two-beam interference, I1 + I2 + 2√(I1 I2) cos φ, with the phase linear in s and one P over the whole range; Gaussian beams, so the fringe amplitude 2√(I1 I2) is exactly one exponential in s, while the background I1 + I2 (two exponentials) is approximated by one; a constant floor D; a camera linear in intensity; and all bins equally good (plain least squares). B ≥ 0, so φ = 0 means a bright center, the same convention as `ring_phase`.

P may move ±30 % from `P_GUESS` and D may go down to −50 counts (the other bounds are in the code). From a poor start a fit can settle on the wrong fringe, so it is run from four starting phases (−π, −π/2, 0, π/2); a start that does not converge is skipped, and the smallest sum of squares wins. φ is wrapped afterwards, and the rms is √(sum of squares / (130 − 7)), in counts.

## `frames(data=DATA)` and `main()`

`frames` yields (S, path) for every `Measurements/<S>/*.png`, settings in increasing S and photographs in time order (the names are time stamps). `main` runs `analyse_frame` on each, prints one line per photograph and writes the outputs.

## Outputs

`fringes_frames.csv`, one row per photograph (161):

| column | content |
|---|---|
| `S` | setting (folder name; 0 = rest) |
| `file` | path from the script folder, `Measurements/<S>/<YYYYMMDDhhmmss>.png` |
| `time_s` | seconds since midnight, from the file name |
| `cx`, `cy` | fitted ring center (px; image coordinates, y downward) |
| `u_cx`, `u_cy` | their standard errors from the circle fit (px) |
| `n_rings` | dark rings recorded (2 or 3 in this run) |
| `R1`, `R2`, `R3` | radii of the innermost dark rings (px); empty if fewer rings |
| `uR1`, `uR2`, `uR3` | their standard errors (px) |
| `n1`, `n2`, `n3` | ray minima per ring (at most 170) |
| `rms_px` | rms distance of the ray minima from their circles (px), over all rings in the fit |
| `P_px2`, `phi_rad` | ring spacing P in r² (px²) and φ_c (rad, in (−π, π]), ring minima |
| `P_cos_px2`, `phi_cos_rad` | the same from the cosine fit |
| `rms_cos` | rms difference between the cosine fit and the profile (counts) |

`fringes_profiles.npz`: `s`, the 130 bin centers (1000 … 259 000 px²); `profiles`, 161 × 130, the smoothed red intensity averaged in each bin about each photograph's fitted center (the profile the cosine fit used); `files`, the `file` column in the same order. `methods.py` and the presentation code read both outputs.

## Constants

| name | value | meaning |
|---|---|---|
| `SMOOTH_PX` | 4.0 | px, Gaussian σ of the image smoothing |
| `THETA_DEG` | (−75.0, 180.0) | deg, fringe sector (anticlockwise from +x, up positive) |
| `THETA_STEP` | 1.5 | deg, ray spacing: 170 rays |
| `R_MAX` | 640 | px, reach of the rays and the annulus profile |
| `MIN_DEPTH` | 3.0 | counts, shallowest dip accepted on a ray |
| `P_GUESS` | 107.5e3 | px², typical ring spacing (windows, numbering, cosine-fit start and bounds) |
| `MIN_POINTS` | 30 | fewest ray minima for a ring to be kept |
| `START` | (900.0, 505.0) | px, starting center (x, y) |
| `N_ITER` | 3 | rounds of profile, minima and circle fit |
| `MAX_RINGS` | 3 | innermost rings recorded |
| `S2_BIN` | 2000.0 | px², bin width of the r² profile |
| `S2_FIT_MAX` | 2.6e5 | px², upper end of the cosine fit |

`START`, `P_GUESS` and `THETA_DEG` are set for the camera view of 11 October.

## Caveats

- The circle-fit errors (`u_cx`, `u_cy`, `uR1` …) treat every ray minimum as independent, but neighboring rays are only 2.6 px apart at r = 100 px (14 px at 530 px), close to the smoothing and speckle scale. RESULTS.md takes the phase uncertainties from the spread of the photographs instead.
- The ring numbering starts at the innermost ring found. A new ring appearing at the center would change φ_c by 2π, which the wrap absorbs; a total shift beyond ±π would need unwrapping downstream. Here `phi_rad` stays between −0.84 and 2.57 rad (`fringes_frames.csv`).
- The cosine fit reads 0.10–0.19 rad higher than the minima at every moving setting, because at rest the bright flat center is not a cosine (rms 9.5 counts at rest, 1.4–2.3 above S = 44; RESULTS.md). Its D is not a measured dark level: rerunning `cosine_fit` on one photograph per setting gives D below zero every time, at the −50 bound for the four lowest settings. The bins share speckle, so the fit reports an rms and no χ².
