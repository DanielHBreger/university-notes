# Sagnac effect: summary (11 October 2026 data)

The final analysis (`calibration.py`, `sagnac.py`; figures in `report_figures/`). Details: `RESULTS.md`. Uncertainties are one standard deviation, rounded to one significant digit with the value to the same digit. Comparisons use unrounded values.

## Result

| | Value |
|---|---|
| Sagnac slope $a = d\Delta\phi/d\omega$ | **1.83 ± 0.06 rad/(rad/s)** |
| Prediction $a_{\mathrm{th}} = 8\pi^2 N R^2/(\lambda c)$ | 1.7 ± 0.1 rad/(rad/s) |
| Agreement | 0.6σ, p = 0.56 (a / a_th = 1.05 ± 0.09) |
| Largest shift | 3.28 ± 0.06 rad = 0.52 fringe at ω = 1.83 rad/s (Δt = 1.1 fs) |
| Speed of light, $c = 8\pi^2 N R^2/(\lambda a)$ | **(2.9 ± 0.2) × 10⁸ m/s**, 0.6σ from c₀, p = 0.56 |

The phase grows linearly with ω and agrees with the Sagnac formula. The comparison is limited by the given turn count N = 13 ± 1 (7.7 %). A steady phase drift in time is not covered by the uncertainty, because the speeds were taken in increasing order; nor is the bend of the calibration at high settings beyond its rms.

## Setup

| N (turns) | R | λ | Fiber |
|---|---|---|---|
| 13 ± 1 (given) | 57.0 ± 0.5 cm (given) | 640 nm (no uncertainty given) | SMF-28e+ |

- Loop centered on the rotation axis. Everything on the table turns with it.
- Calibration: 5 or 10 timed turns at each of 9 controller settings, S = 5.78–90.
- Fringes: 161 photos, 10–22 per setting, rest and 9 speeds (S = 12.59–90), in increasing order, 16:21–16:33.
- Manual tasks done: 1, 2, 3, 5. Task 4 (off-center small loop) skipped on the instructor's instruction.

## Task 1: rotation calibration

- $T$ = mean of n = 5 or 10 consecutive turn times, $\omega = 2\pi/T$. $u(T)^2 = (\mathrm{SD}/\sqrt n)^2 + (\sqrt2\,\sigma_r/n)^2$ with the reaction time σ_r = 0.25 s per press (only the first and last press set the mean of consecutive laps): 0.15–1.1 % of ω.
- Line $\omega = kS + \omega_0$, weighted by the timing uncertainties: k = 0.02054 ± 0.00006 rad/s per unit S, ω₀ = −0.0227 ± 0.0004 rad/s. The points lie farther from the line than their error bars allow (χ²/dof = 20.2/7, p = 0.005): the table turns slower than linear at high settings, 3.0 % below the line at S = 90.
- ω at the photo settings is predicted by the line: its uncertainty combines the line's (0.2–0.3 %) with the 1.3 % scatter of the timed settings about the line, 1.3 % in all; both are common to every setting, so the calibration part of the slope's uncertainty is 1.3 %.

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
| 12.59 | 15 | 0.236(3) | 0.39(3) |
| 20.74 | 15 | 0.403(5) | 0.64(2) |
| 31.11 | 19 | 0.616(8) | 1.02(2) |
| 44.07 | 22 | 0.88(1) | 1.59(2) |
| 55.83 | 18 | 1.12(1) | 1.984(8) |
| 64.00 | 16 | 1.29(2) | 2.375(8) |
| 73.65 | 16 | 1.49(2) | 2.700(4) |
| 82.18 | 15 | 1.67(2) | 2.958(4) |
| 90.00 | 15 | 1.83(2) | 3.285(5) |

ω: the calibration line's prediction uncertainty (1.3 %). Δφ: the table shows the SEM of the photos, the photo noise; the uncertainty of a setting mean in the report is s = 0.04 rad for every setting (the scatter about the line), decision of 24 October.

| Fit result | Value |
|---|---|
| a (regression) | 1.83 ± 0.02 rad/(rad/s) |
| rms distance of the points from the line, s | 0.044 rad |
| R² | 0.9987 |
| Added ω² term | p = 0.67, no bend |
| Intercept | −0.05 ± 0.03 rad, p = 0.12 |
| The line against the photos' SEM | χ²/dof = 205/8: the settings vary five times more than the photo noise, so the regression's own s is used |

## Phase-reading methods

| Method | a (rad/(rad/s)) | R² |
|---|---|---|
| Ring minima (main) | 1.83 ± 0.02 | 0.9987 |
| Cosine fit of the $r^2$ profile | 1.86 ± 0.03 | 0.9982 |
| Template matching to the rest profile | 1.75 ± 0.02 | 0.9993 |
| Fixed pixels, main draw of ten (`pixels.py`; in Fig. 7, not in the method range) | 1.8 ± 0.2 (spread between 1000 draws) | 0.9966 |

- Method uncertainty: half the range of the three full-image methods, 0.054 (2.9 %).
- Fixed pixels: the slope of a draw of ten pixels varies by 0.19 (11 %) between random draws, 3.4 times its own fit error; with 80 pixels still 3 %. The rings move only half a fringe, and fixed pixels read center and spacing changes as phase.

**Final slope:** a = 1.83 ± 0.06 rad/(rad/s), from 1.3 % regression, 2.9 % method and 1.3 % calibration, added in quadrature (3.5 %).

## Speed of light

| Source | Relative |
|---|---|
| Turn count N | 7.7 % |
| Phase-reading method | 2.9 % |
| Loop radius R (squared) | 1.8 % |
| Calibration | 1.3 % |
| Regression | 1.3 % |
| **Total** | **8.6 %** |
| λ (not included) | 0.16 % per nm |

- **c = (2.9 ± 0.2) × 10⁸ m/s**, c / c₀ = 0.95 ± 0.08, 0.6σ, p = 0.56. The same comparison as a / a_th, not an independent test.
- Read the other way, the slope fixes N = 13.7 ± 0.5 turns.

## Task 5: uncertainties

| Source | Kind | Size | Included |
|---|---|---|---|
| Photo-to-photo variation | statistical | SEM 0.004–0.028 rad | yes, in s |
| Wander between settings | statistical | about 0.04 rad | yes, in s |
| Stopwatch (lap scatter, reaction time 0.25 s) | statistical | 0.15–1.1 % of a timed ω | yes, through the line |
| Calibration line and the table's scatter about it | statistical, correlated | 1.3 % of a | yes |
| Phase-reading method | systematic | 2.9 % of a | yes |
| N, R | given | 7.7 %, 1.8 % | yes (prediction, c) |
| λ | not given | 0.16 % per nm | no |
| Bend of the calibration at high settings | systematic | up to 3 % of ω | partly (its rms is in the scatter) |
| Steady drift in time | systematic | not separable | no |

- Drift: joint regression $\Delta\phi = a\omega + Dt + b$ gives a = 2.1 ± 0.2, D = −0.8 ± 0.5 mrad/s, p = 0.12 for D, correlation −0.99: the data cannot tell a steady drift from the slope.
- Plausible source of the wander (not tested): the fiber guides several modes at 640 nm (V = 5.6).
- Not checked: the sign of the shift (rotation direction not recorded).
- Next time: rest photos between speeds, reversed rotation, random order, count the turns, time the turns at the photo settings, single-mode coupling.

## Figures

`report_figures/`: setup, calibration, conversion, patterns, rings, extraction, phase_vs_omega, methods, drift, random_pixels (optional); data and fits only, no residual panels. The slides of 19 October are in `presentation_slides/`; the earlier versions of this file and `RESULTS.md` in `previous_reports/`.
