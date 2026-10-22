# `presentation_slides/source/deck.py`, `slides.py` (analysis part) and `collect_numbers.py`

> **Presentation stage (19 October).** This describes the code that made the slides. The report's numbers now come from `sagnac.py` ([sagnac.md](sagnac.md)), which gives the same slope and c with a simpler uncertainty treatment (no speed-reproducibility term ε, no redrawn calibration lines, no wander estimate). `numbers.json` is superseded by `sagnac.json`. To rerun this code, copy `previous_reports/calibration_eps_2026-10-13.json` to `calibration.json` first.

These three files make the final numbers: a = 1.84 ± 0.07 rad/(rad/s), its comparison with the prediction, c = (2.8 ± 0.2) × 10⁸ m/s, and the checks behind them. `deck.py` loads the data, the analysis functions of `slides.py` do the statistics, and `collect_numbers.py` writes every number to `numbers.json`, the source of `RESULTS.md`. Numbers quoted here are from `numbers.json` unless a file is named.

## The idea

- **One point per setting, a plain regression.** Each setting gives the mean phase of its photos against ω from the calibration line. Weighted by the photos' standard errors (SEM) and the table speed, the line fails the χ² test: the phase wanders slowly, so the photos of one setting are not independent. The line is therefore fitted with equal weights, and the typical distance of the points from the line, s = √(Σr²/(n − 2)), is used as every point's uncertainty; the slope error is s/√(Σ(ω − ω̄)²). χ²/dof = 1 then holds by construction, so R² and an added ω² term are the checks.
- **The calibration line moves all points together**, which s cannot see. It is propagated by redrawing its two coefficients and refitting.
- **Three parts of the slope uncertainty**, added in quadrature: regression and calibration line, phase-reading method, calibration model. For c the relative uncertainties of N, R (twice) and a add in quadrature.

## `deck.py`: data

`DATA` is `$SAGNAC_DATA` or two folders up (the 'Sagnac effect' folder), `UPD = DATA / 'sagnac_update'`, and `OUT` is `$SAGNAC_OUT` or `presentation_slides/`. Every file but the scan is read at import, so a missing one stops the import even if it is unused.

| object | file | use |
|---|---|---|
| `RES` | `sagnac_update/sagnac.json` | a_th = 1.738 ± 0.137 (`theory`, `u_theory`); `slope` = 1.8394 for the speed term of `fit8('photo')` |
| `CAL` | `calibration.json` | the calibration fits |
| `SPEEDS` | `sagnac_update/sagnac_speeds.csv` | per setting: `phi`, `sem`, `omega`, `dphi`, `P`, `drift`, `p_drift`, the other methods' phases and SEMs |
| `FRAMES` | `fringes_frames.csv` | per photo: time, center, `R1`–`R3`, `phi_rad`, `P_px2`, circle-fit rms |
| `SCAN` | `random_points_scan.json` (or `$SAGNAC_SCAN`) | the draws, `None` if missing |

`RES_Q`, `METHODS`, `MFRAMES` and `PROF` (`sagnac_update/optional_quadratic/sagnac.json`, `sagnac_update/methods.json`, `sagnac_update/methods_frames.csv`, `fringes_profiles.npz`) are loaded but used by no script in `source/`.

```python
LIN = CAL['fits_with_eps']['1']          # the linear conversion used by sagnac_update
C0, C1 = LIN['coef']
EPS = CAL['conversion']['eps']
```

`coef` starts with the constant: `C0` = ω₀ = −0.0212 rad/s, `C1` = k = 0.020386 rad/s per unit S (`lin`). ε = 0.48 % (`eps`), the speed reproducibility, is set on the quadratic fit. `omega_lin(S)` is `C0 + C1 * S` for S > 0 and exactly 0 at rest (see `fit8`). It gives `FRAMES` the columns `omega`, `dphi` (`phi_rad` minus the mean of the ten rest photos) and `minutes`; asserts check that `SPEEDS` used the same line.

`calibration_table()` turns the laps read by `load_turns()` into T, u(T) = SD/√(n(n + 1)) (split timing, see `calibration.md`), ω = 2π/T, u_sw = ω u(T)/T and u = √(u_sw² + (εω)²): slide 4 and `cal_table`. Constants: `N_TURNS, R_LOOP, LAM, C_DEF = 13, .57, 640e-9, 299792458.0`, `AREA = π R_LOOP²`.

## `deck.py`: helpers

- `representative(S)`: the photo of setting S whose `phi_rad` is nearest the setting's median (every "typical photo" on the slides).
- `rings(row)`: the finite radii among `R1`–`R3` (two at rest, three from S = 44).
- `crop(img, cx, cy, half=470)`: a 2·half square about (cx, cy), zero-padded past the image edge.
- `sector_profile_r(file, cx, cy, r_max=560, smooth=4.0)`: the red channel smoothed with σ = 4 px and averaged in 1 px annuli over the sector `THETA = (-75.0, 180.0)` deg, as `fringes.sector_profile` (slide 6).
- `ray_minima_for(file, start=(900.0, 505.0))`: the three rounds of the `fringes.py` ring finder on one photo, with what slide 6 draws (center, ray minima, all radii).
- `paren(x, u)`: the notation 0.02039(4). u goes to one significant digit, x to the same place, and a 9.6 that rounds to 10 moves one place up.

## `ring_line(row)` (slides.py)

The line R_m² = P(m + η) of one photo, drawn on slide 7, by the steps of `fringes.ring_phase`, which made the phases in `FRAMES`. Ring numbers m come from the gaps in s = R² over the typical period `107.5e3` px², so a missing ring does not shift the numbering; `np.polyfit(m, s, 1)` gives P and s0, and the function returns m, s, P and η = s0/P. The center phase is φ_c = π − 2πη.

## `fit8(mode, col='phi')` (slides.py)

Δφ = aω + b through the ten setting means of `SPEEDS[col]` (`phi`: ring minima; `phi_cos`, `phi_template`: the other methods), cached per `(mode, col)`. y is the absolute phase, so `b` (−0.797) is not the Δφ intercept.

```python
def omega_of(k0, k1):
    w = k0 + k1 * S
    w[S == 0] = 0.0
    return w
```

**The rest point is at ω = 0 exactly.** The line's intercept is ω₀ = −0.0212 rad/s (tested only for S = 5.78–90), so at S = 0 it would give a negative rotation that did not happen. The table at rest does not turn: ω = 0 has no calibration error and stays 0 when the line is redrawn, which moves only the nine moving points.

**`'regression'`**, the result:

```python
q, *_ = np.linalg.lstsq(A, y, rcond=None)
r = y - A @ q
sd = np.sqrt(r @ r / (n - 2))
return q, sd ** 2 * np.linalg.inv(A.T @ A), np.full(n, sd)
```

Ordinary least squares, A = [ω, 1], equal weights. `sd` is s (n − 2: two parameters came from the same ten points) and becomes every point's uncertainty. s²(AᵀA)⁻¹ is the covariance of (a, b), and its first element is u_stat². In `reg`: a = 1.8390, s = 0.0447, u_stat = 0.0239.

**`'photo'`**, the check that rejects the photos' error bars: weighted least squares with u = √(SEM² + (aεω)²), the photos' SEM and the table speed as a phase, a held at `RES['slope']`. `photo_weighted`: slope 1.843, χ² = 56.6 for 8 dof, p = 2 × 10⁻⁹; the points lie √(56.6/8) = 2.7 times farther from the line than these error bars allow. A third mode, `'wander'`, adds σ_x of `wander_sigma` to u; it only draws the spare slide `08_option_C_measured_wander`.

**Calibration line**, in every mode:

```python
k = np.random.default_rng(1).multivariate_normal([c0, c1], np.array(LIN['cov']), 4000)
u_cal = float(np.std([line(omega_of(k0, k1))[0][0] for k0, k1 in k]))
```

4000 pairs (ω₀, k) from a normal distribution with the line's covariance, seed 1, so every call sees the same lines. For each, ω is recomputed and the line refitted; the SD of the 4000 slopes is u_cal = 0.0037, and u_a = √(u_stat² + u_cal²) = 0.0242.

**Checks and comparison** (in `reg`):

- `chi` = Σ(r/u)² and its `p` for 8 dof. In `'regression'` it is exactly 8 (p = 0.43), since s came from the same residuals: it tests nothing.
- `r2` = 1 − Σr²/Σ(y − ȳ)² = 0.9986.
- `t_q`, `p_q`, the added ω² term: y = qω² + a′ω + b′, unweighted in every mode, covariance from Σr²/(n − 3), t = q/u_q in a two-sided t test with 7 dof. t = 0.48, p = 0.65.
- `z` = (a − a_th)/√(u_a² + u_th²) = 0.72, `p_z` = 2(1 − Φ(|z|)) = 0.47, `ratio` a/a_th = 1.058 ± 0.085. This uses u_a only; the full-uncertainty comparison is in `collect_numbers.py`.

## `wander_sigma(n_boot=20000)` (slides.py)

How much a setting mean wanders beyond photo noise, from the photos alone. Each setting's photos, in time order, are split into the first `n // 2` and the rest:

```python
rows.append((S, a2.mean() - a1.mean(), np.sqrt(a1.var(ddof=1) / len(a1) + a2.var(ddof=1) / len(a2))))
```

Δ is the difference of the half means and u_Δ its photo-noise error. If each half wanders independently by σ_h, Δ has variance 2σ_h² + u_Δ² and the setting mean wanders by σ_h/√2, so σ_x = √⟨Δ² − u_Δ²⟩/2, averaged over the ten settings (0 if negative). Its error: SD over 20000 bootstrap resamplings of the ten settings (seed 0). `wander` = 0.039 ± 0.007 rad, close to s = 0.045.

## `cal_model_slopes()` and `final_result()` (slides.py)

`cal_model_slopes()` repeats the equal-weight regression with ω from the line and from the quadratic fit with ε (`CAL['fits_with_eps']['2']`), rest at 0 in both: a_lin = 1.8390, a_quad = 1.8706 (`final`).

```python
u_meth = (max(sl) - min(sl)) / 2
u_model = abs(a_quad - a_lin)
u = float(np.sqrt(u_reg ** 2 + u_meth ** 2 + u_model ** 2))
```

The parts: u_a of the ring minima (regression and calibration line, 0.024); half the range of the three full-image slopes 1.8390, 1.8761, 1.7676 (method, 0.054, one standard uncertainty by convention; the random pixels are not among them); the full linear-to-quadratic difference (calibration model, 0.032). In quadrature 0.067, so a = 1.84 ± 0.07. Then c = 8πN·AREA/(λa) = 2.834 × 10⁸ m/s. Its relative uncertainty `rel` has five parts: N, `1 / N_TURNS` = 7.7 %; method, u_meth/a = 3.0 %; R, `2 * .005 / R_LOOP` = 1.8 % (R enters squared); calibration model 1.7 %; regression and calibration line 1.3 %. In quadrature 8.7 %, u_c = 0.246 × 10⁸ m/s; λ has no uncertainty given. `z` = (c − c₀)/u_c = −0.67, `p_z` = 0.50.

## `scan_stats(n_boot=2000)` (slides.py)

The statistics of `random_points_scan.json` (see `random_points_scan.md`). The draw with `ORIGINAL_SEED = 20261011` is the single draw of `sagnac_update` (slope 1.786, `orig`) and is kept out; `ten` are the other 1000 ten-pixel draws. For `'fixed'` (the method) and `'frame'` (the diagnostic), `summary` drops slopes ≤ 0 as failed fits (2) and returns the mean with SD/√N, the SD with an error from 2000 bootstrap resamplings (one generator, seed 0, for all summaries in turn), the median of the draws' own fit errors (all 1000), and the count of fits stopped at `max_iter` = 500 rounds (18 fixed-position fits, kept). `by_n` does the same for the fixed-position slopes of 10, 20, 40 and 80 pixels. `scan`: mean 1.702(6), SD 0.196(5), median fit error 0.057.

## `drift_fit()` (slides.py)

Δφ = aω + Dt + b by equal-weight regression, t the mean time of each setting's photos after the first photo, covariance from Σr²/7, and a two-sided t test of D with 7 dof. R² is kept as `drift_fit.r2`. ω and t rose together, so a and D trade off. `drift`: a = 2.14 ± 0.17, D = −0.85 ± 0.47 mrad/s, p = 0.12, correlation −0.99.

## `collect_numbers.py`

Sets `SAGNAC_OUT` to a new temporary folder (no image is written), imports `slides`, calls the functions above and writes `numbers.json` beside itself:

- **Calibration**: `cal_keys`, `cal_fits_keys` (the structure of `calibration.json`); `lin` (k, ω₀, errors, χ² with and without ε); `fit_eps_1`, `_2`, `_3` (the linear, quadratic, cubic fits with ε); `eps`; `cal_table` (slide 4's points); `stopwatch_rel_pct` (range of u_sw/ω, %); `omega90_lin`, `omega90_quad` (ω at S = 90 from each model); `S90` (n, T, SD, ω of that timing).
- **Settings**: `settings` (S, n, first photo's clock time, ω, u_ω, Δφ, SEM, P); `spcols` (columns of `sagnac_speeds.csv`); `theory` [a_th, u_th].
- **Line**: `reg` (every scalar of `fit8('regression')`); `reg_resid` (its residuals); `intercept` [b, u_b, p] (the line refitted to Δφ, t test with 8 dof); `origin_fit` [a, u] (the line through the origin, 9 dof); `photo_weighted` (`fit8('photo')`); `all_photos` [a, u] (a regression through all 161 photos, which treats correlated photos as independent); `wander` [σ_x, u]; `parabola` (= `reg.t_q`, `reg.p_q`).
- **Final**: `final` (`final_result` without `fits`, `rel`); `final_rel` (the five relative parts of c); `methods` (a, u_a, s, R², z, p_z per method); `methods_photo_chi` ([χ², 8, p] with each method's own SEM); `z_full`, `p_full`, `ratio_full` (slope against a_th with the full u: z = 0.66, p = 0.51); `c_over_c0`; `rel_meas` (the three measurement parts, 3.7 %); `N_implied` (13 a/a_th and its errors without and with R).
- **Largest shift**: `max_shift` [Δφ at S = 90, √2 s, Δφ/2π]; `dt_meas_fs` (Δφλ/(2πc)); `dt_pred_fs` (4NAω/c²); `period_fs` (λ/c); `phi_pred_top` (a_th ω).
- **Random pixels**: `scan` (the scalars of `scan_stats`); `scan_robust` (half the 16–84 % width of the good slopes); `scan_pct` (those two percentiles).
- **Drift**: `drift` (`drift_fit`, D in mrad/s, `rho` the a–D correlation); `drift_within` (range, plain mean and SD/√10 of the within-setting drifts; `n_sig` with p < 0.05).
- **Rings**: `rings` (circle-fit rms, radius errors and minima per ring, center errors and SD, median and extreme P).

Nested integers come out as floats (`"n_failed": 2.0`). Rerun on a copy on 22 October 2026, it reproduced the stored `numbers.json` byte for byte.

## Slide drawing

| functions | files in `presentation_slides/` |
|---|---|
| `s01_title`, `s03_system`, `s04_calibration`, `s06_rings`, `s07_extraction` | `01_title`, `03_system`, `04_calibration`, `06_rings_in_r2`, `07_phase_extraction` |
| `s08_result` = `slide8_option('regression', ...)` | `08_phase_vs_omega`; spare modes `08_option_B_...`, `08_option_C_...` |
| `s09_methods`, `s10_random_pixels`, `s11_speed_of_light`, `s12_conclusions`, `s13_thank_you` | `09_methods`, `10_random_pixels`, `11_speed_of_light`, `12_conclusions`, `13_thank_you` |
| `b1_quadratic`, `b2_error_budget`, `b3_data_table`, `b4_drift` | `B1_quadratic_calibration`, `B2_error_budget`, `B3_data_table`, `B4_drift` |

`make_all.py` runs every function in `slides.REGISTRY` (main slides, then backups) and both animations. `anim_concept.py` writes slide 02 and its GIF (two pulses on a turning loop). `anim_measurements.py` writes slide 05 and its GIF (all 161 photos in recording order).

## Caveats

- `fit8('photo', col)` takes the ring-minima SEM `sem` whatever `col` is (line 436). Only `col='phi'` is used that way; `methods_photo_chi` is computed separately.
- `'photo'` holds a at `RES['slope']` where `sagnac_update/sagnac.py` iterates it: χ² = 56.607 here, 56.612 in `sagnac_update/sagnac.json`.
- `scan_stats` shares one generator between summaries, so the same 998 slopes get two bootstrap errors: `scan.fixed.u_sd` = 0.0055 and `scan.by_n['10'].u_sd` = 0.0056, printed 0.196(5) and 0.196(6). Only the SD's error is a bootstrap, though the docstring says both.
- Some slide numbers are typed text, not computed: in `s12_conclusions`, in `b2_error_budget` and in the titles of slides 9, 11 and B1. They agree with `numbers.json` now but will not follow the data.
- `RESULTS.md` also quotes numbers that `numbers.json` lacks: the halves test (8 of 10 settings), r = +0.81 and r = 0.89, D with the first photo's time, and the stopwatch-only quadratic χ² and lag-1 correlations (in `calibration.json`). Recomputed from the data `deck.py` loads, the first four agree.
