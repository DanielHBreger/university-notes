# Sagnac effect: current linear-calibration results

Measurements: 11 October 2026. Standard uncertainties, rounded to one significant figure in report values. Comparisons use unrounded values. Full numerical output: `calibration.txt`, `sagnac.txt`, `sagnac_speeds.csv`.

## Setup and analysis

N = 13 ± 1 turns, R = 57.0 ± 0.5 cm, λ = 640 nm. The large loop is centred on the rotation axis. The data contain 161 photographs at rest and nine moving settings, plus rotation timing at nine settings. The off-centre-loop task was not performed.

Dark minima along rays over −75° to 180° define concentric circles. A fit of squared ring radius against ring order gives the fringe period P and centre phase φc. Differences from rest give Δφ. The rings contract with increasing rotation and are equally spaced in r², with P ≈ 1.074 × 10⁵ px².

## Default rotation calibration

ω(S) = ω₀ + kS, with k = 0.02039 ± 0.00004 rad/s per unit S and ω₀ = -0.0212 ± 0.0006 rad/s. This line supplies every moving-setting angular velocity below. The rest setting has ω = 0 exactly.

Stopwatch splits give u(T) = SD/√[n(n+1)]. The original speed reproducibility ε = 0.48% is retained, estimated from the quadratic residuals. The selected linear coefficients and their covariance are unchanged from the saved straight-line fit.

## Phase shift and rotation

| S | Frames | ω (rad/s) | Δφ (rad) |
|---|---:|---:|---:|
| 0 | 10 | 0 | 0.00 ± 0.05 |
| 12.59 | 15 | 0.235 ± 0.001 | 0.39 ± 0.05 |
| 20.74 | 15 | 0.402 ± 0.002 | 0.64 ± 0.04 |
| 31.11 | 19 | 0.613 ± 0.003 | 1.02 ± 0.05 |
| 44.07 | 22 | 0.877 ± 0.005 | 1.59 ± 0.04 |
| 55.83 | 18 | 1.117 ± 0.006 | 1.98 ± 0.04 |
| 64 | 16 | 1.283 ± 0.007 | 2.38 ± 0.04 |
| 73.65 | 16 | 1.480 ± 0.008 | 2.70 ± 0.04 |
| 82.18 | 15 | 1.654 ± 0.009 | 2.96 ± 0.04 |
| 90 | 15 | 1.814 ± 0.009 | 3.29 ± 0.04 |

Phase entries carry the per-setting fit weight. They are not independent errors on differences from the shared rest reference. The fit uses absolute centre phases with a free intercept.

| Quantity | Result |
|---|---|
| Slope, fit + calibration | 1.84 ± 0.02 rad/(rad/s) |
| Slope, including method contribution | 1.84 ± 0.06 rad/(rad/s) |
| Theory 8π²NR²/(λc) | 1.7 ± 0.1 rad/(rad/s) |
| Difference from theory, including method | 0.68 combined σ, p = 0.50 |
| Maximum measured shift | 3.29 ± 0.06 rad = 0.52 ± 0.01 fringe |
| Quadratic phase term | 0.02 ± 0.05 rad/(rad/s)² |

## Phase-extraction methods

| Method | Slope a, rad/(rad/s) | Phase observable |
|---|---:|---|
| Minima (main) | 1.84 ± 0.02 | Ring centre |
| cosine fit | 1.88 ± 0.03 | Ring centre |
| template matching | 1.77 ± 0.02 | Middle of comparison window |
| random points | 1.79 ± 0.04 | Local phase at ten fixed pixels; excluded from method spread |

Errors in this table combine phase-fit and calibration covariance. Half the range of the three retained slopes gives a method contribution of 0.05468 rad/(rad/s), or 3.0%. The procedures use the same photographs and are not independent measurements. The minima method does not assume an intensity envelope. The random-points method reads ten pixels at fixed image positions; it is excluded because its slope depends on which pixels are drawn (central 68 % of 31 draws: 1.51–1.90 rad/(rad/s)). The refreshed robustness scans are in `methods.txt` and `methods.json`.

## Inferred speed of light and uncertainty

c = 8π²NR²/(λa). This is the same theory comparison expressed as an inferred c, not an independent test.

**c = (2.8 ± 0.2) × 10⁸ m/s**, c/c₀ = 0.94 ± 0.08. Difference from the defined value: -0.68 σ, p = 0.49.

| Contribution | Relative standard uncertainty in c |
|---|---:|
| fit | 1.31% |
| calibration curve | 0.20% |
| method | 2.97% |
| N | 7.69% |
| R | 1.75% |
| Combined | 8.53% |

The measurement alone contributes 3.26%. The per-setting phase error combines SEM, speed reproducibility and σx = 0.0408 rad. Without σx, χ² = 56.6 for 8 dof (frame SEMs alone: 197.4). σx is set by χ² = dof. Calibration-curve covariance propagates as a correlated contribution to the slope.

Steady time drift and wavelength uncertainty are outside the quoted budget. Settings were recorded in increasing order, with rest frames only at the start. The joint rotation/time fit has parameter correlation -0.992, so it cannot separate a rotation slope from a steady drift. The rotation direction was not recorded.

## Optional quadratic calibration check

The default line retains residual curvature: χ²/dof = 47.2/7, p = 5.1e-08. The quadratic has χ²/dof = 6.0/6 by construction because its residuals set ε. The line is the requested baseline approximation, not the statistically preferred calibration model. Residual model curvature is separate from the propagated coefficient covariance.

| Downstream quantity | Linear baseline | Quadratic check |
|---|---:|---:|
| a, rad/(rad/s) | 1.839449 | 1.871349 |
| c, 10⁸ m/s | 2.832794 | 2.784505 |

The quadratic changes the slope by +1.73%. Its complete output is isolated in `optional_quadratic/`, and does not feed the main figures or numbers. Calibration-model sensitivity is reported separately from the quoted standard uncertainty.

## Presentation assets

`presentation_figures/` contains eight main PNGs and one optional model-comparison PNG. The source is `presentation_source/make_figures.py`. No animation or presentation file is delivered.

Reproduce the complete pipeline using `README.md`. Refresh this document and `SUMMARY.md` with `python update_results.py` after both calibration runs.
