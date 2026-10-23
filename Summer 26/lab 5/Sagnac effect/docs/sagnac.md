# `sagnac.py`

Tasks 3 and 5 and the speed of light: the slope a of the phase against ω with its uncertainty, the comparison with the manual's Eq. (5), the same fit for the other phase-reading methods, Eq. (5) solved for c, and the drift check. Every number in `RESULTS.md` comes from here. It reads `calibration.json`, `fringes_frames.csv` and `methods_frames.csv` (if present) from its own folder and writes `sagnac_speeds.csv`, `sagnac.txt` and `sagnac.json`. A second.

## The idea

Eq. (5) predicts Δφ = a_th ω with a_th = 8π²NR²/(λc). The ten setting means are fitted with the line Δφ = aω + b by `scipy.stats.linregress` (Schwartz 2017, Eq. 7): the slope's uncertainty comes from how far the points lie from the line. Measured against the photos' SEM, the line has χ²/dof = 205/8: the settings differ five times more than the photo noise allows, because the phase wanders over minutes, so a setting mean is less well known than its SEM and the regression's own scatter is the right uncertainty.

The slope's uncertainty has three parts, added in quadrature: the regression (0.024), the phase-reading method (half the range of the three methods' slopes, 0.054) and the calibration (every ω is predicted by the same line from the same table, so the relative prediction uncertainty √((u_k/k)² + rel_rms²) = 1.3 % is a common scale error of ω and hence of a: 0.024). a = 1.8261 ± 0.0634 → 1.83 ± 0.06 rad/(rad/s).

## `round1(x, u)`

The report form: u to one significant figure, x to the same digit, half up on the printed digits (`Decimal`, so no binary rounding surprises; `round(2.25, 1)` would give 2.2). `tests/test_rounding.py` pins it down.

## `theory()` and `speed_of_light(a, u_a_parts)`

a_th = 8π²NR²/(λc) with its components from N (u_N/N = 7.7 %) and R (2u_R/R = 1.8 %); c = 8π²NR²/(λa) with the slope's parts carried over at the same relative size (c ∝ 1/a) plus N and R. `U_WAVELENGTH = None`: λ's uncertainty was not given and is left out, stated in the output. `tests/test_speed_of_light.py` checks that a_th gives back the defined c.

## `setting_means(photos, phase_col)`

Per setting, by a pandas groupby: n, the mean phase, SEM (`sem()`, SD/√n), the time of the first photo, the median ring spacing P, and dphi = phase minus the rest mean. The setting means are passed through `np.unwrap` in order of S, a safeguard for angles (here no 2π jump occurs: the whole run spans 3.3 rad and the settings differ by under 0.6 rad).

## `rms_from_line(x, y, fit)`

The rms distance of the points from a `linregress` line with n − 2 degrees of freedom: s, which sets the uncertainty of the largest shift (√2 s) and is reported per method.

## `drift_model(omega_and_t, a, D, b)`

Δφ = aω + Dt + b, the model of the drift check for `curve_fit`, with (ω, t) stacked as the independent variable.

## `random_pixel_summary(path=SCAN)`

Optional (prep Q5): reads the presentation stage's `random_points_scan.json` if present and summarizes, per number of pixels per draw, the draws, the mean and SD of their slopes (negative slopes = failed fits, left out) and the median fit error of one draw. Those slopes were fitted by the presentation's weighted line, not by `linregress`.

## `main()`, step by step

1. `setting_means(photos, 'phi_rad')`; ω from the `line` entry of `calibration.json` through `calibration.omega_of`, exactly 0 at rest; u_omega through `calibration.u_omega_pred` (1.3 %, common to all settings).
2. `stats.linregress(omega, dphi)`: a, b, their standard errors, R²; s from `rms_from_line`; a t test on the intercept.
3. Checks: χ² of the fitted line against the SEMs; an added ω² term by `np.polyfit(..., 2, cov=True)` with a t test on its coefficient (7 dof).
4. The other methods (`OTHER_METHODS`: cosine fit from `fringes_frames.csv`, template matching from `methods_frames.csv`) through the same `setting_means` and `linregress`; u_method = half the range of the three slopes.
5. u_cal = a √((u_k/k)² + rel_rms²) from the `line` entry (ω₀ shifts the moving points together and is absorbed by the intercept, apart from the fixed rest point).
6. u_a in quadrature; a_th and the comparison (ratio, z, p, both with the full u_a and with the regression part alone); N implied = 13 a/a_th; the largest shift (top setting minus rest, u = √2 s), in fringes and as Δt = Δφ λ/(2πc).
7. `speed_of_light` with the three parts; u_c in quadrature; z and p against c₀.
8. Drift: `drift_model` by `curve_fit` with (ω, t) as the two independent variables (t the mean photo time of each setting); a t test on D and the correlation of a and D from the covariance.
9. Outputs.

## Outputs

`sagnac_speeds.csv`, one row per setting: `S`, `n`, `phi` (mean center phase, rad), `sem`, `t_first` (s after midnight), `P` (median, px²), `dphi` (phi minus the rest mean), `omega`, `u_omega` (rad/s), `resid` (dphi minus the line), `dphi_cos`, `sem_cos`, `dphi_template`, `sem_template`.

`sagnac.txt`: the per-setting table, the line and its checks, the other methods, the result with its parts, the theory and comparison, the largest shift, c with its budget, the drift check, the random-pixel summary, and a report-form block (one significant figure).

`sagnac.json`: `slope` (a, u_stat, u_method, u_calibration, u, b, u_b, p_b, s, r2, dof, chi2_vs_sem, p_vs_sem, omega2_term, u_omega2_term, t_omega2, p_omega2, u_k_rel, rel_calibration), `methods` (a, u_a, s, r2 per method), `theory` (a_th, u_a_th, components, per_nm), `comparison` (ratio, u_ratio, z, p, z_regression_only, N_implied, u_N_implied), `largest_shift`, `speed_of_light` (c, u_c, components, rel, ratio_to_defined, z, p), `drift` (a, u_a, D, u_D, p_D, correlation), `random_pixels`, `report_form`.

## Constants

`N_TURNS, U_N_TURNS = 13, 1`; `RADIUS, U_RADIUS = 0.570, 0.005` (m); `WAVELENGTH = 640e-9` (m); `U_WAVELENGTH = None`; `C_LIGHT = 299_792_458.0` (m/s, exact); `OTHER_METHODS`: cosine fit (`phi_cos_rad`), template matching (`phi_template_rad`); `SCAN`: `presentation_slides/source/random_points_scan.json`.

## Against the presentation (19 October)

The presentation's numbers (`presentation_slides/source/numbers.json`, `previous_reports/`) were a = 1.8390 ± 0.0242 (regression plus a calibration-line part from 4000 redrawn lines), method 0.054, calibration model (line against quadratic) 0.032, total 0.067 → 1.84 ± 0.07; here a = 1.8261 ± 0.0236, 0.054, 0.024, total 0.063 → 1.83 ± 0.06. The slope moved by 0.7 % because the calibration line is now fitted with the stopwatch uncertainties (lap scatter and reaction time) as weights, which favor the slow settings, so ω at the top settings is about 1 % higher than before, and because the calibration file was completed on 23 October (10 turns instead of 8 at the four fastest settings). The calibration part (1.3 %) is now the line's prediction uncertainty (k's 0.24 % and the 1.27 % scatter of the timed settings about the line) instead of the quadratic difference; the 3.0 % bend at the top setting is a stated limitation.

## Caveats

- The three methods read the same photographs; half their range is a convention for one standard uncertainty.
- The calibration line sits 3.0 % above the timed ω at S = 90 (the table bends at high settings); a curve through the timed points would lower a by a few percent. The rms of the deviations is in u(a) through the 1.27 % scatter term; the 3.0 % at the top setting is not. Stated in `RESULTS.md`, Task 5.
- A steady drift in time is in no uncertainty, and the sign of the shift is not compared (rotation direction not recorded).
- The earlier versions (`previous_reports/sagnac_quadratic_2026-10-13.py`, `sagnac_update/sagnac.py`) fitted a weighted line with an extra per-setting term set by χ² = dof and propagated the calibration curve by derivatives; they need the earlier `calibration.py`.
