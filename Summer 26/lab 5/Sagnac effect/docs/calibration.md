# `calibration.py`

Turns the stopwatch timing in `Rotation calibration.csv` into the table's angular velocity at each controller setting, ω(S), and fits the line ω = kS + ω₀ that `sagnac.py` uses to convert the photo settings (task 1), with the uncertainty of ω at those settings. Writes `calibration.json` and `calibration.txt`. A fraction of a second.

## The idea

**Timing.** At each setting n = 5 (the five lowest) or 10 (the four highest) consecutive turns were timed, one press per turn. T is the mean turn time and ω = 2π/T. Two things limit T: the scatter of the single laps, SD/√n, and the reaction time of the presses. The laps are consecutive, so the mean is (t_n − t_0)/n and only the first and last press enter: with σ_r = `REACTION_TIME_S` = 0.25 s per press that part is √2 σ_r/n. In quadrature, u(T)² = (SD/√n)² + (√2 σ_r/n)², and u(ω)/ω = u(T)/T: 0.15–1.1 % of ω, the reaction time dominating at the fast settings. (0.25 s is our own reaction time, measured with an online reaction-time test, used as the uncertainty of each press; the lap scatter implies a smaller press jitter, so this is conservative.)

**The line.** `scipy.optimize.curve_fit(line, S, omega, sigma=u_omega, absolute_sigma=True)`: each point weighted by its timing uncertainty, the uncertainties of k and ω₀ propagated from those, nothing rescaled. k = 0.02054 ± 0.00006 rad/s per unit S (0.29 %), ω₀ = −0.0227 ± 0.0004 rad/s. The points lie farther from the line than their uncertainties allow (χ²/dof = 20.2/7, p = 0.005): the table is not exactly linear in S, it turns slower than the line at high settings, 3.0 % below it at S = 90. That misfit is not hidden in the coefficients; it is carried as the rms scatter of the timed settings about the line, 1.28 % of ω.

**ω at a photo setting** (not timed, so predicted by the line): kS + ω₀ with the prediction uncertainty u(ω)² = (line uncertainty)² + (rel_rms ω)², the line's own uncertainty (0.2–0.3 %) and the 1.28 % scatter, 1.3 % in all. Both parts are common to every photo setting (the same line, the same table), so `sagnac.py` takes the calibration part of the Sagnac slope's uncertainty as the same relative size, √((u_k/k)² + rel_rms²) = 1.3 %.

## `load_turns(path=CSV)`

`Speed` is the setting S, `turn 1` … `turn 10` are single-turn times in s: 5 turns at the five lowest settings, 10 at the four highest; empty cells (NaN) are dropped row by row. `encoding='utf-8-sig'` because the file starts with a byte-order mark.

## `per_setting(S, turns, sigma_r=REACTION_TIME_S)`

One row per setting: n, T, the sample SD of the turns (`ddof=1`), `u_T_stat` = SD/√n, `u_T_reaction` = √2 σ_r/n, their quadrature sum `u_T`, ω = 2π/T and `u_omega` = ω u_T/T.

## `line(S, omega0, k)`, `omega_of(S, fit)` and `u_omega_pred(S, fit)`

The model for `curve_fit`; the line at the settings S (`fit` is the `line` entry of `calibration.json`); and the prediction uncertainty, the line's own uncertainty from the covariance of (ω₀, k) in quadrature with `rel_rms` × ω. `sagnac.py` imports `omega_of` and `u_omega_pred`; `methods.py` imports `omega_of`; `report_figures.py` imports both (the band).

## `main()`

Steps 1–3 as above, then the outputs.

## Outputs

`calibration.json` (ω in rad/s):

| key | content |
|---|---|
| `reaction_time_s` | 0.2 |
| `settings` | one record per timed setting: `S`, `n`, `T`, `sd_T`, `u_T_stat`, `u_T_reaction`, `u_T`, `omega`, `u_omega` |
| `line` | `k`, `u_k`, `omega0`, `u_omega0`, `cov` (of ω₀, k), `chi2`, `dof`, `p`, `rel_rms` (0.0128), `resid_rel`, `omega90`, `u_omega90` |
| `omega90_timed` | 1.7724 |

`calibration.txt`: the per-setting table with both parts of u(T), the line with χ²/dof and its deviations, the comparison at S = 90, and ω at the nine photo settings with the prediction uncertainty.

## Constants

| name | value | meaning |
|---|---|---|
| `CSV` | `Rotation calibration.csv` | input, beside the script |
| `REACTION_TIME_S` | 0.2 | s, human reaction time per press |
| `PHOTO_SETTINGS` | 12.59 … 90.0 | the nine photo settings listed in `calibration.txt` |

## Caveats

- The calibration settings (5.78 … 90) are not the photographed ones (12.59 … 90); only S = 90 is common, so ω at the photographs is read off the line and was not timed during them.
- The weights favor the slow settings (their uncertainty is smallest), so the line follows them and sits 3.0 % above the timed ω at S = 90. That bend is a systematic deviation from the line, not scatter; its rms over the nine settings is in the 1.27 % term, its full size at the top setting is not (`RESULTS.md`, Task 5).
- The line's intercept ω₀ is not physical (the line was tested for S = 5.78–90 only); at rest the later scripts set ω = 0.
- The CSV was completed on 23 October: the four fastest settings now hold 10 turns (two were missing before).
- The earlier version of this script (13 October, `previous_reports/calibration_eps_2026-10-13.py`) used a split-timing formula for u(T), a speed-reproducibility term ε set by χ² = dof of a quadratic, and a cubic; the presentation code (`presentation_slides/source/`, `sagnac_update/`) reads its `calibration.json` (`previous_reports/calibration_eps_2026-10-13.json`).
