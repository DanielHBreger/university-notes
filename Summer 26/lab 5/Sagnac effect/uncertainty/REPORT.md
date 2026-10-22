# Measurement uncertainty: every component of the result (22 October 2026)

The uncertainty side of `RESULTS.md`, written so that the report's task 5 and an appendix can be built from it. Every number is from `calibration.txt` and `sagnac.txt` (`calibration.py`, `sagnac.py`). All uncertainties are standard uncertainties (one standard deviation). Kinds: **statistical** (shrinks with more data), **systematic** (a spread between reasonable alternatives), **given** (N and R as stated in the manual). λ has no given uncertainty and is reported as missing, never as zero. The given inputs and the unknown ones are listed in `instrument_inputs.json` (`null` = unknown).

## Overview

| Quantity | How it is obtained | Uncertainty | Kind | Size |
|---|---|---|---|---|
| ω of a timed setting | 5 or 8 consecutive turns, stopwatch | u(T) = SD/√n, u(ω)/ω = u(T)/T | statistical | 0.1–0.5 % of ω |
| Calibration line ω = kS + ω₀ | least squares, weights 1/ω² | from the scatter of the points about the line | statistical | k to 0.55 % |
| ω of a photo setting | the line | rms scatter of the timed points about the line | statistical | 1.2 % of ω |
| φ_c of one photo | ring radii → R_m² = P(m + η) → φ_c = π − 2πη | circle-fit precision (check only) | statistical | about 0.02 rad |
| Δφ of a setting | mean of its photos minus the rest mean | SEM = SD/√n (assumes independent photos) | statistical | 0.004–0.028 rad |
| Slope a, regression | equal-weight line through the ten setting means | s, the rms distance of the points from the line | statistical | 0.024 (1.3 %) |
| Slope a, method | three ways of reading the phase | half their range | systematic | 0.054 (3.0 %) |
| Slope a, calibration | line against quadratic calibration | the difference | systematic | 0.031 (1.7 %) |
| **Slope a** | | three parts in quadrature | | **1.84 ± 0.07 (3.6 %)** |
| Prediction a_th | 8π²NR²/(λc) | N (7.7 %), R squared (1.8 %); λ not given | given | 1.7 ± 0.1 (7.9 %) |
| **Speed of light c** | 8π²NR²/(λa) | relative parts of a, N, R in quadrature | | **(2.8 ± 0.2) × 10⁸ m/s (8.7 %)** |

## 1. Timing of one setting

At each of nine settings, n = 5 or 8 consecutive turns were timed; T is their mean, SD their standard deviation, and

```
u(T) = SD / sqrt(n),     omega = 2 pi / T,     u(omega) / omega = u(T) / T
```

| S | n | ⟨T⟩ (s) | SD (s) | u(T) (s) | ω (rad/s) | u(ω)/ω |
|---:|---:|---:|---:|---:|---:|---:|
| 5.78 | 5 | 65.446 | 0.156 | 0.070 | 0.0960 | 0.11 % |
| 14.67 | 5 | 22.334 | 0.091 | 0.041 | 0.2813 | 0.18 % |
| 25.66 | 5 | 12.424 | 0.039 | 0.018 | 0.5057 | 0.14 % |
| 36.63 | 5 | 8.602 | 0.041 | 0.018 | 0.7304 | 0.21 % |
| 45.09 | 5 | 6.920 | 0.039 | 0.018 | 0.9080 | 0.25 % |
| 56.03 | 8 | 5.565 | 0.045 | 0.016 | 1.1291 | 0.29 % |
| 68.85 | 8 | 4.565 | 0.053 | 0.019 | 1.3764 | 0.41 % |
| 79.21 | 8 | 3.994 | 0.042 | 0.015 | 1.5733 | 0.37 % |
| 90.00 | 8 | 3.545 | 0.054 | 0.019 | 1.7724 | 0.54 % |

The turns were timed as consecutive laps, so neighboring turn times share a button press and SD/√n is slightly conservative. Nothing depends on it: the scatter about the calibration line (section 2) is ten times larger.

## 2. The calibration line and ω at the photo settings

The line ω = kS + ω₀ is fitted by least squares with each point weighted by 1/ω², i.e. by minimizing Σ((ω_i − kS_i − ω₀)/ω_i)². The weights are relative because the points scatter about any smooth curve by a roughly constant fraction of ω; an unweighted fit would be set by the fast settings and miss the slowest one by 13 %. The coefficient covariance is scaled by the residual scatter (χ²_w/dof), so it says how well the points fix the line:

```
k = 0.02038 +- 0.00011 rad/s per unit S (0.55 %),   omega0 = -0.0212 +- 0.0015 rad/s
(omega - line)/omega: -0.6, +1.3, +0.8, +0.7, +1.1, +0.8, -0.4, -1.2, -2.3 %   (rms 1.15 %)
```

Against the stopwatch uncertainties alone the line has χ²/dof = 184/7 (p = 3 × 10⁻³⁶), so the 1.2 % scatter is the table, not the stopwatch: the speed is not exactly reproducible at a setting and not exactly linear in S (a quadratic fits to 0.40 %, and the line is 2.3 % high at S = 90). The photos take ω from the line, as task 1 asks, and

```
u(omega at a photo setting) = 1.2 % of omega
```

covers the two together. In the regression of Δφ on ω (section 4) this uncertainty is one of the things that move the points off the line and is contained in s; its systematic part, the bend, is measured separately by the quadratic (section 6). The 0.55 % uncertainty of k describes the same misfit and is not propagated a second time.

## 3. Phase of one photo and of one setting

Per photo, the dark-ring radii come from concentric circles fitted to the minima found along 170 rays, and the line R_m² = P(m + η) through the two or three rings gives φ_c = π − 2πη. The line has at most one degree of freedom, so no fit statistic of its own; the circle fit is its check: the single minima lie 2.8–8.2 px rms from their circles (speckle), and with 54–170 minima per ring a radius is known to 0.3–2.3 px (median 0.8 px), which propagates to about 0.02 rad on φ_c. This precision is a check only; the setting means take their uncertainty from the photos.

Per setting, Δφ is the mean of the photos minus the mean of the 10 rest photos, with u = SEM = SD/√n (0.004–0.028 rad). The rest mean (SEM 0.020 rad) is a common offset of every setting, which a line with a free intercept absorbs. The SEM assumes that the photos of a setting are independent; they are not, because the phase wanders over minutes (section 4 and the drift section of `RESULTS.md`), so a setting mean is less well known than its SEM.

## 4. The regression

```
Delta phi = a omega + b,   ordinary least squares, ten setting means, equal weights (Schwartz 2017, Eq. 7)
s = sqrt( sum (Delta phi_i - a omega_i - b)^2 / (10 - 2) ) = 0.0448 rad
u(a) = s / sqrt( sum (omega_i - mean omega)^2 ) = 0.0240
```

s is the rms distance of the points from the line and is used as every point's uncertainty: it contains everything that moves a setting mean (photo noise, the table speed, the wander between settings). Weighted by the photos' SEM instead, the line gives χ²/dof = 197/8 (p = 2 × 10⁻³⁸): the points lie five times farther from the line than the photo noise allows, which is why the regression is the right tool here. a = 1.8398 ± 0.0240, reported as 1.84 ± 0.02 rad/(rad/s); R² = 0.9986; an added ω² term is 0.5 standard errors from zero (p = 0.65).

With the regression uncertainty only: z = (1.8398 − 1.7381)/√(0.0240² + 0.1371²) = 0.73, two-sided p = 0.47.

## 5. The method part

The cosine fit and template matching, through the same averaging and regression, give 1.8769 and 1.7684 against the minima's 1.8398. Half the range of the three, 0.0543 (3.0 %), is taken as one standard uncertainty. The three read the same photographs, so they are correlated checks, not independent measurements; half the range is a convention, and it is the largest part of u(a).

## 6. The calibration part

The regression repeated with ω from the quadratic calibration gives a = 1.8709, 1.7 % above the line's 1.8398; the difference, 0.0311, is the calibration uncertainty of a. It is the only place where the line's misfit enters the result.

## 7. The final slope

| Part | On a | Relative | Kind |
|---|---:|---:|---|
| Regression | 0.0240 | 1.3 % | statistical |
| Phase-reading method | 0.0543 | 3.0 % | systematic |
| Calibration | 0.0311 | 1.7 % | systematic |
| **Quadrature** | **0.0670** | **3.6 %** | a = 1.8398 ± 0.0670 → **1.84 ± 0.07 rad/(rad/s)** |

The three parts come from different things and are added as independent. Against a_th: a/a_th = 1.058 ± 0.092, z = 0.67, p = 0.51. Read the other way, with R, λ and c known, N = 13 a/a_th = 13.76 ± 0.56 turns.

## 8. The prediction and the given inputs

```
a_th = 8 pi^2 N R^2 / (lambda c) = 1.7381 s
u(a_th)/a_th = sqrt( (u_N/N)^2 + (2 u_R/R)^2 ) = sqrt( 7.69^2 + 1.75^2 ) % = 7.89 %   ->   1.7 +- 0.1 rad/(rad/s)
```

c is exact (defined). λ has no given uncertainty; a_th and c change by 0.16 % per nm, so λ would have to be off by about 23 nm to matter as much as all our own uncertainties together. Whether the given ±1 turn and ±0.5 cm are standard uncertainties or bounds is not stated; they are treated as standard uncertainties (as a bound, ±1 turn would be u = 0.58 turn and the comparison z = 0.9).

## 9. The speed of light

```
c = 8 pi^2 N R^2 / (lambda a),     (u_c/c)^2 = (u_a/a)^2 + (u_N/N)^2 + (2 u_R/R)^2   (first order; c ~ N R^2 / a)
```

| Source | Relative uncertainty of c |
|---|---:|
| Turn count N = 13 ± 1 | 7.7 % |
| Phase-reading method | 3.0 % |
| Loop radius R (squared) | 1.8 % |
| Calibration | 1.7 % |
| Regression | 1.3 % |
| **Total** | **8.7 %** |

c = 2.832 ± 0.246 × 10⁸ m/s → (2.8 ± 0.2) × 10⁸ m/s; c/c₀ = 0.94 ± 0.08, z = −0.67, p = 0.50. This is a/a_th read the other way (c/c₀ = a_th/a), not an independent test.

## 10. The largest shift

```
Delta phi(S = 90) - Delta phi(rest) = 3.285 rad,   u = sqrt(2) s = 0.063 rad (two setting means, each known to s)
= 0.523 +- 0.010 fringe;   Delta t = Delta phi lambda / (2 pi c) = 1.12 +- 0.02 fs
```

Predicted at ω = 1.81 rad/s: 3.15 rad, 1.07 fs, half an optical period (λ/c = 2.13 fs).

## 11. What the quoted uncertainty does not cover

- **A steady drift of the phase in time.** The settings were photographed in increasing order with rest photos only at the start, so ω and time rise together. The joint regression Δφ = aω + Dt + b gives a = 2.1 ± 0.2, D = −0.8 ± 0.5 mrad/s (p = 0.12) with a correlation of −0.99 between a and D: the data cannot separate a steady drift from the slope. The wander between settings is in s; a steady drift is not.
- **The wavelength** (no tolerance given; section 8).
- **The sign of the shift.** The rotation direction was not recorded; only magnitudes are compared.
- **ω during the photographs.** Not measured; it comes from the calibration at the same setting, done at other times. The 1.2 % scatter stands for how well the table repeats a setting.
- **Camera geometry.** The phase uses only ratios of squared radii (η, P), so no pixel-to-length calibration enters. A pixel aspect ratio different from 1 would turn the circles into ellipses; the minima lie only 3–8 px rms from the fitted circles, so any such effect is small, but the aspect ratio itself is not known.
- **A plausible source of the wander, untested:** at 640 nm the fiber guides several modes (V = 5.6 with a = 4.1 µm, NA = 0.14; single-mode needs V < 2.405), and light in several modes is sensitive to stress and temperature; the camera records the whole field with no single-mode filter at the common port.

## Files

- `instrument_inputs.json`: the given values and the explicit missing inputs (`null`).
- `calibration.py`, `sagnac.py` (the folder above): every number here; `calibration.txt`, `sagnac.txt`, `sagnac.json` hold them unrounded.
- The earlier, longer audit (22 October, presentation stage) is in `previous_reports/uncertainty_audit_2026-10-22/`.
