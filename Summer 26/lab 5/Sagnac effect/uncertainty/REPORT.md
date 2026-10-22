# Measurement uncertainty: every component of the result (22 October 2026)

The uncertainty side of `RESULTS.md`, written so that the report's task 5 and an appendix can be built from it. Every number is from `calibration.txt` and `sagnac.txt` (`calibration.py`, `sagnac.py`). All uncertainties are standard uncertainties (one standard deviation). Kinds: **statistical** (shrinks with more data), **systematic** (a spread between reasonable alternatives), **given** (N and R as stated in the manual). λ has no given uncertainty and is reported as missing, never as zero. The given inputs and the unknown ones are listed in `instrument_inputs.json` (`null` = unknown).

## Overview

| Quantity | How it is obtained | Uncertainty | Kind | Size |
|---|---|---|---|---|
| T of a timed setting | mean of n = 5 or 8 consecutive turns | lap scatter SD/√n and reaction time √2 σ_r/n in quadrature | statistical | 0.14–1.1 % of ω |
| Calibration line ω = kS + ω₀ | weighted least squares, weights 1/u(ω)² | covariance × χ²/dof (3.2) | statistical | k to 0.48 % |
| ω of a photo setting | predicted by the line | the line's uncertainty and the rms scatter of the timed settings about it, in quadrature | statistical, common to all settings | 1.45 % of ω |
| φ_c of one photo | ring radii → R_m² = P(m + η) → φ_c = π − 2πη | circle-fit precision (check only) | statistical | about 0.02 rad |
| Δφ of a setting | mean of its photos minus the rest mean | SEM = SD/√n (assumes independent photos) | statistical | 0.004–0.028 rad |
| Slope a, regression | equal-weight line through the ten setting means | s, the rms distance of the points from the line | statistical | 0.024 (1.3 %) |
| Slope a, method | three ways of reading the phase | half their range | systematic | 0.054 (2.9 %) |
| Slope a, calibration | every ω predicted by the same line | a √((u_k/k)² + rel_rms²) | statistical, correlated | 0.027 (1.5 %) |
| **Slope a** | | three parts in quadrature | | **1.82 ± 0.06 (3.5 %)** |
| Prediction a_th | 8π²NR²/(λc) | N (7.7 %), R squared (1.8 %); λ not given | given | 1.7 ± 0.1 (7.9 %) |
| **Speed of light c** | 8π²NR²/(λa) | relative parts of a, N, R in quadrature | | **(2.9 ± 0.2) × 10⁸ m/s (8.6 %)** |

## 1. Timing of one setting

At each of nine settings, n = 5 or 8 consecutive turns were timed with one press per turn. T is their mean and SD their standard deviation. The turns are consecutive laps, so the mean turn time is (t_n − t_0)/n and only the first and the last press enter it; with a human reaction time σ_r = 0.2 s per press that part is √2 σ_r/n. The lap scatter gives the statistical part SD/√n. Together

```
u(T)^2 = (SD/sqrt(n))^2 + (sqrt(2) sigma_r / n)^2,     omega = 2 pi / T,     u(omega)/omega = u(T)/T
```

| S | n | ⟨T⟩ (s) | SD (s) | SD/√n (s) | √2σ_r/n (s) | u(T) (s) | ω (rad/s) | u(ω)/ω |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5.78 | 5 | 65.446 | 0.156 | 0.070 | 0.057 | 0.090 | 0.0960 | 0.14 % |
| 14.67 | 5 | 22.334 | 0.091 | 0.041 | 0.057 | 0.070 | 0.2813 | 0.31 % |
| 25.66 | 5 | 12.424 | 0.039 | 0.018 | 0.057 | 0.059 | 0.5057 | 0.48 % |
| 36.63 | 5 | 8.602 | 0.041 | 0.018 | 0.057 | 0.059 | 0.7304 | 0.69 % |
| 45.09 | 5 | 6.920 | 0.039 | 0.018 | 0.057 | 0.059 | 0.9080 | 0.86 % |
| 56.03 | 8 | 5.565 | 0.045 | 0.016 | 0.035 | 0.039 | 1.1291 | 0.70 % |
| 68.85 | 8 | 4.565 | 0.053 | 0.019 | 0.035 | 0.040 | 1.3764 | 0.88 % |
| 79.21 | 8 | 3.994 | 0.042 | 0.015 | 0.035 | 0.038 | 1.5733 | 0.96 % |
| 90.00 | 8 | 3.545 | 0.054 | 0.019 | 0.035 | 0.040 | 1.7724 | 1.13 % |

The reaction time dominates at the fast settings, where a turn takes 3–7 s. σ_r = 0.2 s is a conventional value; the measured lap scatter (0.04–0.16 s) implies a press-to-press jitter well below it, so the uncertainties are on the safe side. SD/√n treats the laps as independent, which is also conservative for consecutive laps.

## 2. The calibration line and ω at the photo settings

The line ω = kS + ω₀ is fitted by weighted least squares with weights 1/u(ω)². Its χ²/dof = 22.5/7 = 3.2 (p = 0.002): the points scatter more than their uncertainties allow, because the table is not exactly linear in S (it turns slower than the line at high settings: the deviations are −0.1, +0.8, +0.1, −0.1, +0.3, −0.1, −1.3, −2.1, −3.2 % in order of S). The coefficient covariance is therefore multiplied by χ²/dof, the usual remedy, so that the uncertainties reflect the actual scatter:

```
k = 0.02058 +- 0.00010 rad/s per unit S (0.48 %),   omega0 = -0.0229 +- 0.0006 rad/s
```

A photo setting was not timed: the line predicts its ω, and a prediction carries two things, how well the line is known and how far a single setting lies from it. The second is the rms scatter of the timed settings about the line, 1.38 % of ω (the table does not repeat a setting better than that). So

```
u(omega)^2 = (1, S) Cov (1, S)^T + (0.0138 omega)^2        1.45 % of omega
```

(the line's part is 0.3–0.45 %). Both parts are common to all photo settings (the same line, the same table), so they do not average out between settings; their effect on the slope is a separate part (section 6).

What is only partly in these numbers: the bend of the table's speed at high settings is a systematic deviation from the line (3.2 % at S = 90), not scatter. Its rms over the nine settings is in the 1.38 %; its full size at the top setting is not. A curve through the timed points instead of the line would lower the Sagnac slope by a few percent; this is stated in `RESULTS.md` (Task 5).

## 3. Phase of one photo and of one setting

Per photo, the dark-ring radii come from concentric circles fitted to the minima found along 170 rays, and the line R_m² = P(m + η) through the two or three rings gives φ_c = π − 2πη. The line has at most one degree of freedom, so no fit statistic of its own; the circle fit is its check: the single minima lie 2.8–8.2 px rms from their circles (speckle), and with 54–170 minima per ring a radius is known to 0.3–2.3 px (median 0.8 px), which propagates to about 0.02 rad on φ_c. This precision is a check only; the setting means take their uncertainty from the photos.

Per setting, Δφ is the mean of the photos minus the mean of the 10 rest photos, with u = SEM = SD/√n (0.004–0.028 rad). The rest mean (SEM 0.020 rad) is a common offset of every setting, which a line with a free intercept absorbs. The SEM assumes that the photos of a setting are independent; they are not, because the phase wanders over minutes (section 4 and the drift section of `RESULTS.md`), so a setting mean is less well known than its SEM.

## 4. The regression

```
Delta phi = a omega + b,   ordinary least squares, ten setting means, equal weights (Schwartz 2017, Eq. 7)
s = sqrt( sum (Delta phi_i - a omega_i - b)^2 / (10 - 2) ) = 0.0444 rad
u(a) = s / sqrt( sum (omega_i - mean omega)^2 ) = 0.0235
```

s is the rms distance of the points from the line and is used as every point's uncertainty: it contains everything that moves a setting mean (photo noise, the table speed, the wander between settings). Weighted by the photos' SEM instead, the line gives χ²/dof = 197/8 (p = 3 × 10⁻³⁸): the points lie five times farther from the line than the photo noise allows, which is why the regression is the right tool here. a = 1.8228 ± 0.0235, reported as 1.82 ± 0.02 rad/(rad/s); R² = 0.9987; an added ω² term is 0.5 standard errors from zero (p = 0.67).

With the regression uncertainty only: z = (1.8228 − 1.7381)/√(0.0235² + 0.1371²) = 0.61, two-sided p = 0.54.

## 5. The method part

The cosine fit and template matching, through the same averaging and regression, give 1.8595 and 1.7520 against the minima's 1.8228. Half the range of the three, 0.0538 (2.9 %), is taken as one standard uncertainty. The three read the same photographs, so they are correlated checks, not independent measurements; half the range is a convention, and it is the largest part of u(a).

## 6. The calibration part

Every moving ω is predicted by the same line from the same table, so the relative prediction uncertainty of section 2 is a common scale error of ω, and a scale error of ω is the same scale error of a: u_cal(a)/a = √((u_k/k)² + 0.0138²) = √(0.48² + 1.38²) % = 1.46 %, u_cal(a) = 0.0266. (ω₀ shifts the nine moving points together and is absorbed by the line's intercept, apart from the fixed rest point; its effect is far below the rounding.) This is the only place where the calibration enters u(a); the per-setting u(ω) of section 2 are the same information point by point. Treating the scatter as common is the conservative reading: if it were independent between settings it would partly average out in the slope, but the deviations are a bend, not noise.

## 7. The final slope

| Part | On a | Relative | Kind |
|---|---:|---:|---|
| Regression | 0.0235 | 1.3 % | statistical |
| Phase-reading method | 0.0538 | 2.9 % | systematic |
| Calibration | 0.0266 | 1.5 % | statistical (correlated across settings) |
| **Quadrature** | **0.0645** | **3.5 %** | a = 1.8228 ± 0.0645 → **1.82 ± 0.06 rad/(rad/s)** |

The three parts come from different things and are added as independent. Against a_th: a/a_th = 1.049 ± 0.091, z = 0.56, p = 0.58. Read the other way, with R, λ and c known, N = 13 a/a_th = 13.63 ± 0.54 turns.

## 8. The prediction and the given inputs

```
a_th = 8 pi^2 N R^2 / (lambda c) = 1.7381 s
u(a_th)/a_th = sqrt( (u_N/N)^2 + (2 u_R/R)^2 ) = sqrt( 7.69^2 + 1.75^2 ) % = 7.89 %   ->   1.7 +- 0.1 rad/(rad/s)
```

c is exact (defined). λ has no given uncertainty; a_th and c change by 0.16 % per nm, so λ would have to be off by about 22 nm to matter as much as all our own uncertainties together. Whether the given ±1 turn and ±0.5 cm are standard uncertainties or bounds is not stated; they are treated as standard uncertainties.

## 9. The speed of light

```
c = 8 pi^2 N R^2 / (lambda a),     (u_c/c)^2 = (u_a/a)^2 + (u_N/N)^2 + (2 u_R/R)^2   (first order; c ~ N R^2 / a)
```

| Source | Relative uncertainty of c |
|---|---:|
| Turn count N = 13 ± 1 | 7.7 % |
| Phase-reading method | 2.9 % |
| Loop radius R (squared) | 1.8 % |
| Calibration | 1.5 % |
| Regression | 1.3 % |
| **Total** | **8.6 %** |

c = 2.859 ± 0.247 × 10⁸ m/s → (2.9 ± 0.2) × 10⁸ m/s; c/c₀ = 0.95 ± 0.08, z = −0.56, p = 0.57. This is a/a_th read the other way (c/c₀ = a_th/a), not an independent test.

## 10. The largest shift

```
Delta phi(S = 90) - Delta phi(rest) = 3.285 rad,   u = sqrt(2) s = 0.063 rad (two setting means, each known to s)
= 0.523 +- 0.010 fringe;   Delta t = Delta phi lambda / (2 pi c) = 1.12 +- 0.02 fs
```

Predicted at ω = 1.83 rad/s: 3.18 rad, 1.08 fs, half an optical period (λ/c = 2.13 fs).

## 11. What the quoted uncertainty does not cover

- **The bend of the calibration beyond its rms.** The table turns slower than the line at high settings (3.2 % at S = 90); a curve through the timed points would lower a by a few percent. The rms of the deviations (1.4 %) is in the calibration part; the full 3.2 % at the top setting is not.
- **A steady drift of the phase in time.** The settings were photographed in increasing order with rest photos only at the start, so ω and time rise together. The joint regression Δφ = aω + Dt + b gives a = 2.1 ± 0.2, D = −0.8 ± 0.5 mrad/s (p = 0.13) with a correlation of −0.99 between a and D: the data cannot separate a steady drift from the slope. The wander between settings is in s; a steady drift is not.
- **The wavelength** (no tolerance given; section 8).
- **The sign of the shift.** The rotation direction was not recorded; only magnitudes are compared.
- **ω during the photographs.** Not measured; it comes from the calibration at the same setting, done at other times.
- **Camera geometry.** The phase uses only ratios of squared radii (η, P), so no pixel-to-length calibration enters. A pixel aspect ratio different from 1 would turn the circles into ellipses; the minima lie only 3–8 px rms from the fitted circles, so any such effect is small, but the aspect ratio itself is not known.
- **A plausible source of the wander, untested:** at 640 nm the fiber guides several modes (V = 5.6 with a = 4.1 µm, NA = 0.14; single-mode needs V < 2.405), and light in several modes is sensitive to stress and temperature; the camera records the whole field with no single-mode filter at the common port.

## Files

- `instrument_inputs.json`: the given values and the explicit missing inputs (`null`).
- `calibration.py`, `sagnac.py` (the folder above): every number here; `calibration.txt`, `sagnac.txt`, `sagnac.json` hold them unrounded. The reaction time is `REACTION_TIME_S` in `calibration.py`.
- The earlier, longer audit (22 October, presentation stage) is in `previous_reports/uncertainty_audit_2026-10-22/`.
