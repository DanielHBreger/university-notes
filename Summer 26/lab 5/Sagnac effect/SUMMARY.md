# Sagnac effect: summary (11 October 2026 data)

The final analysis (`calibration.py`, `sagnac.py`; figures in `report_figures/`). Details: `RESULTS.md`. Uncertainties are one standard deviation, rounded to one significant digit with the value to the same digit. Comparisons use unrounded values.

## Result

| | Value |
|---|---|
| Sagnac slope $a = d\Delta\phi/d\omega$ | **1.84 ± 0.07 rad/(rad/s)** |
| Prediction $a_{\mathrm{th}} = 8\pi^2 N R^2/(\lambda c)$ | 1.7 ± 0.1 rad/(rad/s) |
| Agreement | 0.7σ, p = 0.51 (a / a_th = 1.06 ± 0.09) |
| Largest shift | 3.28 ± 0.06 rad = 0.52 fringe at ω = 1.81 rad/s (Δt = 1.1 fs) |
| Speed of light, $c = 8\pi^2 N R^2/(\lambda a)$ | **(2.8 ± 0.2) × 10⁸ m/s**, 0.7σ from c₀, p = 0.50 |

The phase grows linearly with ω and agrees with the Sagnac formula. The comparison is limited by the given turn count N = 13 ± 1 (7.7 %). A steady phase drift in time is not covered by the uncertainty, because the speeds were taken in increasing order.

## Setup

| N (turns) | R | λ | Fiber |
|---|---|---|---|
| 13 ± 1 (given) | 57.0 ± 0.5 cm (given) | 640 nm (no uncertainty given) | SMF-28e+ |

- Loop centered on the rotation axis. Everything on the table turns with it.
- Calibration: 5–8 timed turns at each of 9 controller settings, S = 5.78–90.
- Fringes: 161 photos, 10–22 per setting, rest and 9 speeds (S = 12.59–90), in increasing order, 16:21–16:33.
- Manual tasks done: 1, 2, 3, 5. Task 4 (off-center small loop) skipped on the instructor's instruction.

## Task 1: rotation calibration

- $T$ = mean of 5 or 8 turn times, $u(T) = \mathrm{SD}/\sqrt{n}$, $\omega = 2\pi/T$. Stopwatch: 0.1–0.5 % of ω.
- Line $\omega = kS + \omega_0$ (least squares with equal relative weights): k = 0.0204 ± 0.0001 rad/s per unit S, ω₀ = −0.021 ± 0.002 rad/s. The timed points scatter about it by 1.2 % rms, far more than the stopwatch allows (χ²/dof = 184/7): the table turns slower than linear at high settings, 2.3 % below the line at S = 90.
- A quadratic fits to 0.4 %. The line is used, as the task asks; the quadratic raises the Sagnac slope by 1.7 %, and that is the calibration part of the uncertainty.

## Task 2: fringe patterns

- Rings, evenly spaced in $r^2$: $\phi(r) = 2\pi r^2/P + \phi_c$ for two Gaussian beams of different wavefront radii, $P = 2\lambda/|1/R_1 - 1/R_2|$. Median P = 1.07 × 10⁵ px².
- Rotation changes only $\phi_c$: the rings contract, the first dark ring moves from 258 px to 101 px, the center goes from bright to nearly dark.
- Ring finder: minima along 170 rays (−75° to 180°), concentric circles fitted. Minima lie 3–8 px (rms) from their circles, ring radii known to about ±1 px.

## Task 3: phase shift against ω

1. Per photo: $R_m^2 = P(m + \eta)$ fitted against ring number m, $\phi_c = \pi - 2\pi\eta$ (manual Eq. 7, Schwartz Eq. 5, in $r^2$). No pixel calibration needed.
2. $\Delta\phi = \phi_c - \phi_c(\text{rest})$, mean per setting.
3. Linear regression $\Delta\phi = a\omega + b$ through the ten setting means, equal weights (Schwartz Eq. 7).

| S | Photos | ω (rad/s) | Δφ (rad) |
|---|---|---|---|
| rest | 10 | 0 | 0 |
| 12.59 | 15 | 0.235(3) | 0.39(3) |
| 20.74 | 15 | 0.402(5) | 0.64(2) |
| 31.11 | 19 | 0.613(7) | 1.02(2) |
| 44.07 | 22 | 0.88(1) | 1.59(2) |
| 55.83 | 18 | 1.12(1) | 1.984(8) |
| 64.00 | 16 | 1.28(1) | 2.375(8) |
| 73.65 | 16 | 1.48(2) | 2.700(4) |
| 82.18 | 15 | 1.65(2) | 2.958(4) |
| 90.00 | 15 | 1.81(2) | 3.285(5) |

ω: the line's 1.2 % scatter. Δφ: standard error of the mean of the photos (SEM); the rest SEM (0.02 rad) is common to all rows.

| Fit result | Value |
|---|---|
| a (regression) | 1.84 ± 0.02 rad/(rad/s) |
| rms distance of the points from the line, s | 0.045 rad |
| R² | 0.9986 |
| Added ω² term | p = 0.65, no bend |
| Intercept | −0.05 ± 0.03 rad, p = 0.11 |
| Weighted by the photos' SEM instead | χ²/dof = 197/8: the settings vary five times more than the photo noise, so the regression's own s is used |

## Phase-reading methods

| Method | a (rad/(rad/s)) | R² |
|---|---|---|
| Ring minima (main) | 1.84 ± 0.02 | 0.9986 |
| Cosine fit of the $r^2$ profile | 1.88 ± 0.03 | 0.9982 |
| Template matching to the rest profile | 1.77 ± 0.02 | 0.9993 |
| 10 fixed pixels, 1000 draws (not used) | 1.70, spread 0.20 between draws | |

- Method uncertainty: half the range of the three full-image methods, 0.054 (3.0 %).
- Fixed pixels: one draw of ten pixels varies by 0.20 between draws, 3.4 times its own fit error. The rings move only half a fringe, and fixed pixels read center and spacing changes as phase.

**Final slope:** a = 1.84 ± 0.07 rad/(rad/s), from 1.3 % regression, 3.0 % method and 1.7 % calibration, added in quadrature (3.6 %).

## Speed of light

| Source | Relative |
|---|---|
| Turn count N | 7.7 % |
| Phase-reading method | 3.0 % |
| Loop radius R (squared) | 1.8 % |
| Calibration | 1.7 % |
| Regression | 1.3 % |
| **Total** | **8.7 %** |
| λ (not included) | 0.16 % per nm |

- **c = (2.8 ± 0.2) × 10⁸ m/s**, c / c₀ = 0.94 ± 0.08, 0.7σ, p = 0.50. The same comparison as a / a_th, not an independent test.
- Read the other way, the slope fixes N = 13.8 ± 0.6 turns.

## Task 5: uncertainties

| Source | Kind | Size | Included |
|---|---|---|---|
| Photo-to-photo variation | statistical | SEM 0.004–0.028 rad | yes, in s |
| Table speed | statistical | 1.2 % of ω | yes, in s |
| Wander between settings | statistical | about 0.04 rad | yes, in s |
| Calibration curve (line vs quadratic) | systematic | 1.7 % of a | yes |
| Phase-reading method | systematic | 3.0 % of a | yes |
| N, R | given | 7.7 %, 1.8 % | yes (prediction, c) |
| λ | not given | 0.16 % per nm | no |
| Steady drift in time | systematic | not separable | no |

- Drift: joint regression $\Delta\phi = a\omega + Dt + b$ gives a = 2.1 ± 0.2, D = −0.8 ± 0.5 mrad/s, p = 0.12 for D, correlation −0.99: the data cannot tell a steady drift from the slope.
- Plausible source of the wander (not tested): the fiber guides several modes at 640 nm (V = 5.6).
- Not checked: the sign of the shift (rotation direction not recorded).
- Next time: rest photos between speeds, reversed rotation, random order, count the turns, single-mode coupling.

## Figures

`report_figures/`: setup, calibration, patterns, rings, extraction, phase_vs_omega, methods, drift, random_pixels (optional). The slides of 19 October are in `presentation_slides/`; the earlier versions of this file and `RESULTS.md` in `previous_reports/`.
