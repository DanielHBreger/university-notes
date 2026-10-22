# `calibration.py`

Turns the stopwatch timing in `Rotation calibration.csv` into the table's angular velocity at each controller setting, ω(S), and fits the line ω = kS + ω₀ that `sagnac.py` uses to convert the photo settings (task 1). Writes `calibration.json` and `calibration.txt`. A fraction of a second.

## The idea

**Timing.** At each setting n = 5 or 8 consecutive turns were timed. T is the mean turn time, u(T) = SD/√n its standard error, ω = 2π/T and u(ω)/ω = u(T)/T: 0.1–0.5 % of ω. (The turns were timed as consecutive laps, so neighboring turn times share a button press and SD/√n is slightly conservative; nothing downstream depends on it.)

**The line, with relative weights.** The nine points scatter about any smooth curve by about 1 % of ω, ten times the stopwatch uncertainty, and that scatter is a roughly constant fraction of ω. The line is therefore fitted by least squares with weights 1/ω² (minimizing the relative residuals); an unweighted fit would be set by the fast settings and miss the slowest one by 13 %. The coefficient covariance is scaled by the residual scatter (χ²_w/dof), so u_k says how well the points fix the line, not how good the stopwatch is.

**What comes out.** k = 0.02038 ± 0.00011 rad/s per unit S, ω₀ = −0.0212 ± 0.0015 rad/s; relative residuals −0.6, +1.3, +0.8, +0.7, +1.1, +0.8, −0.4, −1.2, −2.3 % (rms 1.15 %); against the stopwatch errors χ²/dof = 184/7. The table turns slower than linear at high settings: 2.3 % below the line at S = 90. A quadratic fits to 0.40 %. The photos use the line, as task 1 asks, with the 1.2 % scatter as the uncertainty of ω at a setting; the quadratic is the check that sets the calibration part of the Sagnac slope's uncertainty (`sagnac.py`).

## `load_turns(path=CSV)`

`Speed` is the setting S, `turn 1` … `turn 10` are single-turn times in s: 5 turns at the five lowest settings, 8 at the four highest; empty cells (NaN) are dropped row by row. `encoding='utf-8-sig'` because the file starts with a byte-order mark.

## `per_setting(S, turns)`

One row per setting: n, T, the sample SD of the turns (`ddof=1`), u_T = SD/√n, ω = 2π/T and u_omega = ω u_T/T.

## `fit_poly(x, y, deg, weights)`

Weighted least squares written out with the design matrix A (columns 1, x, x², …):

```python
cov0 = np.linalg.inv(A.T @ (A * weights[:, None]))
c = cov0 @ (A.T @ (weights * y))
r = y - A @ c
cov = cov0 * float(np.sum(weights * r ** 2)) / dof        # scaled by the residual scatter
```

Returns the coefficients (constant first), their covariance, the residuals and the rms of the relative residuals (y − curve)/y. `main()` calls it with `weights = 1 / omega**2` for the line (deg 1) and the quadratic (deg 2).

## `omega_of(S, coef)`

The fitted polynomial at the settings S; `sagnac.py`, `methods.py` and `report_figures.py` import it with the coefficients from `calibration.json`.

## `main()`

Steps 1–3 as above, then `chi2_stopwatch` = Σ((ω − line)/u_omega)² for the line, and the outputs.

## Outputs

`calibration.json` (ω in rad/s; coefficients constant first):

| key | content |
|---|---|
| `settings` | one record per timed setting: `S`, `n`, `T`, `sd_T`, `u_T`, `omega`, `u_omega` |
| `line` | `coef`, `cov`, `k`, `u_k`, `omega0`, `u_omega0`, `rel_rms` (rms of the relative residuals, 0.0115), `resid_rel`, `chi2_stopwatch`, `dof`, `p_stopwatch`, `omega90` |
| `quadratic` | `coef`, `rel_rms` (0.0040), `resid_rel`, `omega90` |
| `omega90_timed` | 1.7724 |

`calibration.txt`: the per-setting table, both fits with their residuals, the comparison at S = 90, and ω at the nine photo settings from the line (± the rms scatter) and from the quadratic.

## Constants

| name | value | meaning |
|---|---|---|
| `CSV` | `Rotation calibration.csv` | input, beside the script |
| `PHOTO_SETTINGS` | 12.59 … 90.0 | the nine photo settings listed in `calibration.txt` |

## Caveats

- The calibration settings (5.78 … 90) are not the photographed ones (12.59 … 90); only S = 90 is common, so ω at the photographs is read off the line and was not timed during them. The 1.2 % scatter stands for how well the table repeats a setting.
- The line's intercept ω₀ is not physical (the line was tested for S = 5.78–90 only); at rest the later scripts set ω = 0.
- The earlier version of this script (13 October, `previous_reports/calibration_eps_2026-10-13.py`) used a split-timing formula for u(T), a speed-reproducibility term ε set by χ² = dof of the quadratic, and a cubic; the presentation code (`presentation_slides/source/`, `sagnac_update/`) reads its `calibration.json` (`previous_reports/calibration_eps_2026-10-13.json`). The line itself is the same to four digits.
