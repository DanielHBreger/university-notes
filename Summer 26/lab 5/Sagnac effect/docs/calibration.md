# `calibration.py`

Turns the stopwatch timing in `Rotation calibration.csv` into the table's angular velocity at each controller setting, ω(S), and fits the line ω = kS + ω₀ that `sagnac.py` uses to convert the photo settings (task 1), with the line's uncertainty. Writes `calibration.json` and `calibration.txt`. A fraction of a second.

## The idea

**Timing.** At each setting n = 5 or 8 consecutive turns were timed, one press per turn. T is the mean turn time and ω = 2π/T. Two things limit T: the scatter of the single laps, SD/√n, and the reaction time of the presses. The laps are consecutive, so the mean is (t_n − t_0)/n and only the first and last press enter: with σ_r = `REACTION_TIME_S` = 0.2 s per press that part is √2 σ_r/n. In quadrature, u(T)² = (SD/√n)² + (√2 σ_r/n)², and u(ω)/ω = u(T)/T: 0.14–1.1 % of ω, the reaction time dominating at the fast settings. (0.2 s is a conventional value; the lap scatter implies a smaller press jitter, so this is conservative.)

**The line.** Weighted least squares with weights 1/u(ω)². The points scatter more than their uncertainties allow (χ²/dof = 22.5/7: the table is not exactly linear in S, it turns slower than the line at high settings, 3.2 % below it at S = 90), so the coefficient covariance is multiplied by χ²/dof, the usual remedy. k = 0.02058 ± 0.00010 rad/s per unit S (0.48 %), ω₀ = −0.0229 ± 0.0006 rad/s.

**ω at a photo setting** (not timed, so predicted by the line): kS + ω₀ with the prediction uncertainty u(ω)² = (1, S) Cov (1, S)ᵀ + (rel_rms ω)², the line's own uncertainty (0.3–0.45 %) and the rms scatter of the timed settings about the line (1.38 %), 1.45 % in all. Both parts are common to every photo setting (the same line, the same table), so `sagnac.py` takes the calibration part of the Sagnac slope's uncertainty as the same relative size, √((u_k/k)² + rel_rms²) = 1.46 %.

## `load_turns(path=CSV)`

`Speed` is the setting S, `turn 1` … `turn 10` are single-turn times in s: 5 turns at the five lowest settings, 8 at the four highest; empty cells (NaN) are dropped row by row. `encoding='utf-8-sig'` because the file starts with a byte-order mark.

## `per_setting(S, turns, sigma_r=REACTION_TIME_S)`

One row per setting: n, T, the sample SD of the turns (`ddof=1`), `u_T_stat` = SD/√n, `u_T_reaction` = √2 σ_r/n, their quadrature sum `u_T`, ω = 2π/T and `u_omega` = ω u_T/T.

## `fit_line(x, y, u)`

Weighted least squares written out with the design matrix A = [1, x]:

```python
cov = np.linalg.inv(A.T @ (A * W[:, None]))          # W = 1/u^2
c = cov @ (A.T @ (W * y))
chi2 = float(np.sum(((y - A @ c) / u) ** 2))
scale = max(1.0, chi2 / dof)
return c, cov * scale, chi2, dof, p, scale
```

The covariance is enlarged by χ²/dof when the points scatter more than their error bars (never shrunk when they scatter less). Returns the coefficients (constant first), the scaled covariance, χ², dof, p and the scale.

## `omega_of(S, coef)`, `u_omega_of(S, cov)` and `u_omega_pred(S, line)`

The line at the settings S; its own uncertainty √((1, S) Cov (1, S)ᵀ); and the prediction uncertainty, that in quadrature with `rel_rms` × ω (`line` is the `line` entry of `calibration.json`). `sagnac.py` imports `omega_of` and `u_omega_pred`; `methods.py` and `report_figures.py` import `omega_of`.

## `main()`

Steps 1–3 as above, then the outputs.

## Outputs

`calibration.json` (ω in rad/s; coefficients constant first):

| key | content |
|---|---|
| `reaction_time_s` | 0.2 |
| `settings` | one record per timed setting: `S`, `n`, `T`, `sd_T`, `u_T_stat`, `u_T_reaction`, `u_T`, `omega`, `u_omega` |
| `line` | `coef`, `cov` (scaled), `k`, `u_k`, `omega0`, `u_omega0`, `chi2`, `dof`, `p`, `cov_scale`, `rel_rms` (rms of the relative deviations, 0.014), `resid_rel`, `omega90`, `u_omega90` |
| `omega90_timed` | 1.7724 |

`calibration.txt`: the per-setting table with both parts of u(T), the line with χ²/dof and its deviations, the comparison at S = 90, and ω at the nine photo settings with the prediction uncertainty (and the line's alone, for comparison).

## Constants

| name | value | meaning |
|---|---|---|
| `CSV` | `Rotation calibration.csv` | input, beside the script |
| `REACTION_TIME_S` | 0.2 | s, human reaction time per press |
| `PHOTO_SETTINGS` | 12.59 … 90.0 | the nine photo settings listed in `calibration.txt` |

## Caveats

- The calibration settings (5.78 … 90) are not the photographed ones (12.59 … 90); only S = 90 is common, so ω at the photographs is read off the line and was not timed during them.
- The weights favor the slow settings (their relative uncertainty is smallest), so the line follows them and sits 3.2 % above the timed ω at S = 90. That bend is a systematic deviation from the line, not scatter; its rms over the nine settings is in the 1.38 % term, its full size at the top setting is not (`RESULTS.md`, Task 5).
- The line's intercept ω₀ is not physical (the line was tested for S = 5.78–90 only); at rest the later scripts set ω = 0.
- The earlier version of this script (13 October, `previous_reports/calibration_eps_2026-10-13.py`) used a split-timing formula for u(T), a speed-reproducibility term ε set by χ² = dof of a quadratic, and a cubic; the presentation code (`presentation_slides/source/`, `sagnac_update/`) reads its `calibration.json` (`previous_reports/calibration_eps_2026-10-13.json`).
