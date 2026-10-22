# `sagnac.py`

Tasks 3 and 5 and the speed of light: the slope a of the phase against ω with its uncertainty, the comparison with the manual's Eq. (5), the same fit for the other phase-reading methods, Eq. (5) solved for c, and the drift check. Every number in `RESULTS.md` comes from here. It reads `calibration.json`, `fringes_frames.csv` and `methods_frames.csv` (if present) from its own folder and writes `sagnac_speeds.csv`, `sagnac.txt` and `sagnac.json`. A second.

## The idea

Eq. (5) predicts Δφ = a_th ω with a_th = 8π²NR²/(λc). The ten setting means are fitted with the line Δφ = aω + b by ordinary least squares (Schwartz 2017, Eq. 7): every point's uncertainty is s, the rms distance of the points from the line, and u(a) = s/√Σ(ω − ω̄)². Weighted by the photos' SEM instead, the line has χ²/dof = 197/8: the settings differ five times more than the photo noise allows, because the phase wanders over minutes, so a setting mean is less well known than its SEM and the regression's own s is the right uncertainty.

The slope's uncertainty has three parts, added in quadrature: the regression (0.024), the phase-reading method (half the range of the three methods' slopes, 0.054) and the calibration (the change of a when ω comes from the quadratic calibration instead of the line, 0.031). a = 1.8398 ± 0.0670 → 1.84 ± 0.07 rad/(rad/s).

## `round1(x, u)`

The report form: u to one significant figure, x to the same digit, half up on the printed digits (`Decimal`, so no binary rounding surprises; `round(2.25, 1)` would give 2.2). `tests/test_rounding.py` pins it down.

## `theory()` and `speed_of_light(a, u_a_parts)`

a_th = 8π²NR²/(λc) with its components from N (u_N/N = 7.7 %) and R (2u_R/R = 1.8 %); c = 8π²NR²/(λa) with the slope's parts carried over at the same relative size (c ∝ 1/a) plus N and R. `U_WAVELENGTH = None`: λ's uncertainty was not given and is left out, stated in the output. `tests/test_speed_of_light.py` checks that a_th gives back the defined c.

## `regress(x, y)`

```python
q, *_ = np.linalg.lstsq(A, y, rcond=None)          # A = [x, 1]
r = y - A @ q
s2 = r @ r / (len(y) - 2)
cov = s2 * np.linalg.inv(A.T @ A)
```

Ordinary least squares with the covariance from the residuals (n − 2: two parameters came from the same points). Returns a, b, u_a, u_b, s, R² and the residuals. `methods.py` imports it for its robustness scans.

## `setting_means(fr, col)`

Per setting: n, the mean phase, SEM = SD/√n, the time of the first photo, the median ring spacing P, and dphi = phase minus the rest mean. Two safeguards for angles known modulo 2π: each photo is brought within π of the setting's circular mean before averaging, and each setting mean within π of the previous one. Here neither changes a number (the whole run spans 3.3 rad, the settings differ by under 0.6 rad).

## `scan_summary(path=SCAN)`

Optional (prep Q5): reads the presentation stage's `random_points_scan.json` if present and summarizes, per number of pixels per draw, the draws, the mean and SD of their slopes (negative slopes = failed fits, left out) and the median fit error of one draw. Those slopes were fitted by the presentation's weighted line, not by `regress`.

## `main()`, step by step

1. `setting_means(fr, 'phi_rad')`; ω from `calibration.json['line']['coef']` through `calibration.omega_of`, exactly 0 at rest; u_omega = `rel_rms` × ω (the line's 1.2 % scatter); ω from the quadratic as well, for step 6.
2. `regress(omega, dphi)`: a, b, u_stat, s, R².
3. Checks: the line weighted by the SEMs (χ², p); an added ω² term by least squares, with a t test on its coefficient (7 dof); the line through the origin; a t test on the intercept.
4. The other methods (`OTHER_METHODS`: cosine fit from `fringes_frames.csv`, template matching from `methods_frames.csv`) through the same `setting_means` and `regress`; u_meth = half the range of the three slopes.
5. `regress(omega_quadratic, dphi)`; u_model = |a_quad − a|.
6. u_a in quadrature; a_th and the comparison (ratio, z, p, both with the full u_a and with u_stat alone); N implied = 13 a/a_th; the largest shift (top setting minus rest, u = √2 s), in fringes and as Δt = Δφ λ/(2πc).
7. `speed_of_light` with the three parts; u_c in quadrature; z and p against c₀.
8. Drift: Δφ = aω + Dt + b by least squares (t the mean photo time of each setting), covariance from the residuals, a t test on D and the correlation of a and D.
9. Outputs.

## Outputs

`sagnac_speeds.csv`, one row per setting: `S`, `n`, `phi` (mean center phase, rad), `sem`, `t_first` (s after midnight), `P` (median, px²), `dphi` (phi minus the rest mean), `omega`, `u_omega` (rad/s), `resid` (dphi minus the line), `dphi_cos`, `sem_cos`, `dphi_template`, `sem_template`.

`sagnac.txt`: the per-setting table, the line and its checks, the other methods, the result with its parts, the theory and comparison, the largest shift, c with its budget, the drift check, the random-pixel summary, and a report-form block (one significant figure).

`sagnac.json`: `slope` (a, u_stat, u_method, u_calibration, u, b, u_b, p_b, s, r2, dof, chi2_sem_weighted, p_sem_weighted, omega2_term, u_omega2_term, t_omega2, p_omega2, a_through_origin, u_a_through_origin, a_quadratic_calibration), `methods` (a, u_a, s, r2 per method), `theory` (a_th, u_a_th, components, per_nm), `comparison` (ratio, u_ratio, z, p, z_regression_only, N_implied, u_N_implied), `largest_shift`, `speed_of_light` (c, u_c, components, rel, ratio_to_defined, z, p), `drift` (a, u_a, D, u_D, p_D, correlation), `random_pixels`, `report_form`.

## Constants

`N_TURNS, U_N_TURNS = 13, 1`; `RADIUS, U_RADIUS = 0.570, 0.005` (m); `WAVELENGTH = 640e-9` (m); `U_WAVELENGTH = None`; `C_LIGHT = 299_792_458.0` (m/s, exact); `OTHER_METHODS`: cosine fit (`phi_cos_rad`), template matching (`phi_template_rad`); `SCAN`: `presentation_slides/source/random_points_scan.json`.

## Against the presentation (19 October)

The presentation's numbers (`presentation_slides/source/numbers.json`, `previous_reports/`) were a = 1.8390 ± 0.0242 (regression 0.0239 plus a 0.0037 calibration-line part from 4000 redrawn lines), method 0.054, calibration model 0.032, total 0.067; here a = 1.8398 ± 0.0240, 0.054, 0.031, total 0.067. The line is the same to four digits (relative weights there came from the ε term, here from 1/ω²), so every reported value agrees after rounding; the one visible change is the largest shift, 3.28 ± 0.06 rad (from the unrounded 3.2850, where the presentation rounded a six-digit 3.285 up).

## Caveats

- The three methods read the same photographs; half their range is a convention for one standard uncertainty.
- A steady drift in time is in no uncertainty, and the sign of the shift is not compared (rotation direction not recorded).
- The earlier versions (`previous_reports/sagnac_quadratic_2026-10-13.py`, `sagnac_update/sagnac.py`) fitted a weighted line with an extra per-setting term set by χ² = dof and propagated the calibration curve by derivatives; they need the earlier `calibration.py`.
