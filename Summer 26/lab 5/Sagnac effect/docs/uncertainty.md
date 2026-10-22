# `uncertainty/calculate.py`

> **Superseded (22 October).** `uncertainty/calculate.py` and its outputs were moved to `previous_reports/uncertainty_audit_2026-10-22/`; they import the presentation code and the earlier `calibration.json`. The uncertainty treatment of the report is in `uncertainty/REPORT.md` and `sagnac.py` ([sagnac.md](sagnac.md)).

Recomputes every uncertainty component of the final result from the inputs, states the model behind each, and checks the ones that rest on an assumption. It changes nothing: the final numbers stay those of `RESULTS.md`, and the script checks itself against `presentation_slides/source/numbers.json`. `uncertainty/REPORT.md` is the written version, section by section.

## The idea

Like `report_figures.py`, the script imports the final analysis (`deck.py` for the data, `slides.py` for `fit8`, `final_result`, `wander_sigma`) instead of repeating it, and adds only what the analysis does not already print: the derivation behind each component, the checks, and the alternatives. The given inputs and the unknown ones come from `uncertainty/instrument_inputs.json`, where unknown means `null`, never zero.

## Sections, in the order of the output

1. **Stopwatch.** `deck.load_turns()` gives the laps of each timed setting. Per setting:

   ```python
   u_T = sd / math.sqrt(n * (n + 1))                  # lap splits share their presses
   s_press = sd * math.sqrt(n / (2 * (n + 1)))        # timing error of one press
   u_res = INP['stopwatch_display_resolution_s'] / math.sqrt(12) / math.sqrt(n)
   ```

   With an independent error σ per press, a lap has variance 2σ² and neighbors covariance −σ², so E[SD²] = 2σ²(n+1)/n and Var(T) = 2σ²/n²; the two lines invert that. The display resolution (0.01 s, read from the data) adds a uniform rounding error per displayed lap; it is printed to show it is negligible. The lag-1 check is read from `calibration.json`.

2. **Speed reproducibility ε.** `chi2_poly(deg, eps)` is a weighted polynomial fit of ω(S) with errors √(u_sw² + (εω)²). `brentq` solves χ²(quadratic) = dof for ε (0.483 %, the value in `calibration.json`) and, for comparison, χ²(line) = dof (1.30 %). The table of χ²/dof and p for degree 1–3 with and without ε is the one of `RESULTS.md`.

3. **Calibration line.** The line's coefficients and covariance from `calibration.json` `fits_with_eps['1']`; for each photo setting, u_curve = √((1, S) Cov (1, S)ᵀ) and u_repro = εω, checked against `sagnac_update/sagnac_speeds.csv`.

4. **Phase of one photo.** `phi_c(R)` repeats `fringes.ring_phase` (unwrapped): squared radii, ring numbers from the gaps divided by `P_GUESS`, a straight line, φ_c = π − 2π s₀/P. For every photo, the derivative of φ_c with respect to each radius (finite differences, 10⁻³ px) times the radius uncertainty from the circle fit (`uR1`–`uR3`), added in quadrature as if the radii were independent. Per setting, the median of this precision is set against the SD of the photos' phases and the within-setting drift of `sagnac_speeds.csv`. A check only; it does not enter the result.

5. **Setting means.** SEMs, the rest SEM, and `slides.wander_sigma()`: the difference between the second and first half of each setting's photos, with its photo-noise error, gives the wander of a setting mean (bootstrap over the ten settings); the number of settings whose halves differ by more than 2 errors, and by how much.

6. **Regression.** `slides.fit8('regression')` (equal weights; s; the calibration-line part from 4000 redrawn lines) and, for contrast, `fit8('photo')` (weights from the SEM and ε only), with the comparison against a_th using the regression error only. Then the Pearson correlation of the points' distances from the line (`fit8`'s `resid`) with each setting's fringe period, mean and median of the photos' `P_px2`, and of the ring-fit with the cosine-fit mean periods: the r = +0.81 and r = 0.89 that `RESULTS.md` quotes and `numbers.json` lacks.

7. **Final slope.** `slides.final_result()`: the three parts and their quadrature sum; then the same comparison with the parts combined in other reasonable ways (straight sum, method part as range/√12, as the full range, and the calibration-line part enlarged by √(χ²/dof) of the line, since `calibration.py` does not scale the line's covariance).

8. **Prediction.** a_th from the inputs file, its uncertainty from N and R, the sensitivity to λ per nm, the λ error that would equal our own 3.7 %, and the comparison if ±1 turn is a bound (u = 1/√3 turn).

9. **Speed of light.** The first-order relative budget of `final_result`, a Monte Carlo of c = 8π²NR²/(λa) with N, R and a drawn from normal distributions (400 000 draws, fixed seed), and the budget with N as a bound.

10. **Largest shift.** Top setting minus rest, u = √2 s, in fringes and as Δt = Δφλ/(2πc).

11. **Missing inputs.** Every `null` in `instrument_inputs.json`, plus what the data cannot test (sign, steady drift, ω during the photos).

12. **Agreement.** Ten values (slope, its regression uncertainty, s, the final uncertainty, c and u_c, a_th and its uncertainty, ε, the wander) compared with `numbers.json` to a relative 10⁻⁶; the count of differences is printed.

## Outputs

- `uncertainty/calculate.txt`: the printout above.
- `uncertainty/budget.json`: every number, unrounded, under keys `1_stopwatch` … `12_numbers_json_agreement`.

```
python uncertainty/calculate.py     # from the 'Sagnac effect' folder, a few seconds
```

## Caveats

- The radius uncertainties of section 4 come from the circle fit, which treats neighboring ray minima as independent although they share speckle; the photo precision there is probably somewhat optimistic. The section's conclusion rests on the settings where the photos differ several times more than that precision.
- The Monte Carlo uses normal distributions for N and R as given; with N as a bound, section 8's numbers apply instead.
