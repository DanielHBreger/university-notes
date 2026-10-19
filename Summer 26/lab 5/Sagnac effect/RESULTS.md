# Sagnac effect: results (11 October 2026 data)

This is the final analysis, the one presented in `presentation_slides/`. It replaces the earlier version of this document, which used the quadratic calibration (kept in `previous_reports/`). The short version is `SUMMARY.md`.

**Conventions.** Statistical terms (SEM, regression, R², χ²/dof, two-sided p, t test, bootstrap, robust spread) are explained in `presentation_slides/presenter_guide.md`, Part 1.5. All uncertainties are standard uncertainties (one standard deviation). Values are given rounded, with the uncertainty to its first significant digit and the value to the same digit. A number like 0.02039(4) means 0.02039 ± 0.00004. Comparisons (σ, p) always use the unrounded values, which are listed where they are needed. Every fit reports χ²/dof and p. A plain linear regression, whose points have no independent error bars, reports R² instead (Task 3 explains why).

**Where the numbers come from.** `presentation_slides/source/collect_numbers.py` computes every number below with the presentation's own code (`slides.py`, `deck.py`) and writes them to `numbers.json`. Its inputs are `Rotation calibration.csv`, `calibration.json`, `fringes_frames.csv`, the linear-calibration outputs in `sagnac_update/` and the random-pixel scan `presentation_slides/source/random_points_scan.json`. The top-level scripts `calibration.py` and `sagnac.py` still default to the quadratic calibration, so the top-level `sagnac.txt` and `sagnac.json` belong to the earlier analysis.

**What changed from the earlier version.**

| | Earlier version | This version |
|---|---|---|
| Calibration for the photos | quadratic | linear, as task 1 asks. The linear-to-quadratic difference is an uncertainty |
| Line through the ten settings | weighted fit with an extra per-setting term σ_x set by χ²/dof = 1 | ordinary linear regression with equal weights (Schwartz Eq. 7). On these data it gives the same slope and error as the σ_x-weighted fit |
| Phase-reading methods for the method uncertainty | minima, cosine fit, template matching, 2D Fourier (AIA shown, not used) | minima, cosine fit, template matching |
| Random-pixel method | not used | 1000 random draws (the intermediate `sagnac_update` analysis used one draw) |
| Slope a | 1.87 ± 0.06 | 1.84 ± 0.07 |
| c | (2.8 ± 0.2) × 10⁸ m/s, 0.9σ from c₀ | (2.8 ± 0.2) × 10⁸ m/s, 0.7σ from c₀ |

---

## Setup and given parameters

| Quantity | Value | Source |
|---|---|---|
| Fiber turns $N$ | 13 ± 1 | given (manual, Fig. 1 caption) |
| Loop radius $R$ | 57.0 ± 0.5 cm | given (manual, Fig. 1 caption) |
| Wavelength $\lambda$ | 640 nm, no uncertainty given | given |
| Fiber | Corning SMF-28e+ (the manual's text says SMF-28e, its printed datasheet and the preparation sheet say SMF-28e+). Core radius 4.1 µm, NA 0.14 | datasheet in the manual |
| Enclosed area of one turn | $A = \pi R^2 = 1.02\ \mathrm{m^2}$ | |

- The large loop is centered on the rotation axis. Everything on the table (laser, polarizers, beam splitter, coupling lenses, fiber, camera) turns with it.
- Data: 161 photographs at rest and at nine motor settings, 10–22 photos per setting, taken in increasing order of speed between 16:21:40 and 16:33:46. Rotation timing at nine settings.
- Manual tasks done: 1, 2, 3 and 5. Task 4 (the off-center small loop, R = 15.5 ± 0.1 cm, 46 ± 1 turns) was skipped on the instructor's instruction.

**Prediction.** The Sagnac phase (manual Eq. 5) for $N$ turns enclosing $A$ each is

$$\Delta\phi = \frac{8\pi N A}{\lambda c}\,\omega \equiv a_{\mathrm{th}}\,\omega, \qquad a_{\mathrm{th}} = \frac{8\pi^2 N R^2}{\lambda c} = 1.7 \pm 0.1\ \mathrm{rad/(rad/s)}$$

Unrounded $a_{\mathrm{th}} = 1.738 \pm 0.137$. The uncertainty is from $N$ (7.7 %) and $R$ (squared, 2 × 0.88 % = 1.8 %), $\sqrt{7.7^2 + 1.8^2} = 7.9\ \%$. The unit rad/(rad/s) equals seconds.

---

## Task 1: rotation calibration

Figure: `presentation_slides/04_calibration.png`, backup `B1_quadratic_calibration.png`.

**Timing.** At each of nine controller settings $S$ (a number on the motor controller, no unit, which sets the motor voltage) 5 or 8 consecutive turns were timed as stopwatch splits. The mean turn time is $T = (t_n - t_0)/n$, and $\omega = 2\pi/T$.

- Splits share their button presses, so only the first and last press set $T$. If each press has an independent error $\sigma_r$, then $\mathrm{Var}(T) = 2\sigma_r^2/n^2$ while the turn times have $E[\mathrm{SD}^2] = 2\sigma_r^2 (n+1)/n$. So $u(T) = \mathrm{SD}/\sqrt{n(n+1)}$.
- Check that the turns are splits: neighboring turn times are anticorrelated, lag-1 correlation −0.34, against −0.43 ± 0.09 expected for splits (simulated).
- Stopwatch uncertainty: 0.04–0.18 % of ω. Example, S = 90: 8 turns, ⟨T⟩ = 3.545 s, SD = 0.054 s, u(T) = 0.006 s, ω = 1.772 rad/s ± 0.18 %.
- Speed reproducibility ε = 0.48 % of ω: the spread of the settings about a smooth curve that the stopwatch cannot explain. It is set by χ²/dof = 1 for the quadratic and added in quadrature to every ω.

| Model | Errors | χ²/dof | p |
|---|---|---|---|
| Linear | stopwatch only | 807/7 | 6 × 10⁻¹⁷⁰ |
| **Linear (used)** | + ε | 47.2/7 | 5 × 10⁻⁸ |
| Quadratic | stopwatch only | 154.9/6 | 7 × 10⁻³¹ |
| Quadratic (check) | + ε | 6.0/6 (by construction) | 0.42 |
| Cubic (check) | + ε | 2.6/5 | 0.75 |

**Linear calibration (used).** $\omega = kS + \omega_0$ with $k = 0.02039(4)$ rad/s per unit S and $\omega_0 = -0.0212(6)$ rad/s. The line is rejected (p = 5 × 10⁻⁸): the table turns slower than linear at high settings. At S = 90 the line gives 1.814 rad/s, 2.3 % above the timed 1.772 rad/s. The quadratic gives 1.782 rad/s there. The intercept is not physical: the line was only tested for S = 5.78–90.

The photos take ω from the line, as the task asks. The quadratic is used only as a check, and the 1.7 % it changes the Sagnac slope is part of the final uncertainty (calibration model, see Task 3).

| Downstream quantity | Linear (used) | Quadratic (check) | Change |
|---|---|---|---|
| ω at S = 90 (rad/s) | 1.814 | 1.782 | −1.8 % |
| Sagnac slope a (rad/(rad/s)) | 1.839 | 1.871 | +1.7 % |
| a / a_th | 1.058 | 1.076 | |
| c (10⁸ m/s) | 2.834 | 2.786 | −1.7 % |

---

## Task 2: fringe patterns

Figures: `presentation_slides/05_measurements.png` with `05_measurements_ANIMATION.gif`, `06_rings_in_r2.png`.

- The pattern is a set of rings. As ω rises, the rings contract and the first dark ring at rest (radius 258 px) moves in to 101 px at the top speed. The center goes from bright to nearly dark: about half a fringe.
- The rings are evenly spaced in $r^2$, not in $r$. In the fastest photo the dark rings are 241 px and then 132 px apart in $r$, but 1.08 × 10⁵ px² both times in $r^2$.
- **Why.** The two beams leave the fiber as Gaussian beams with different wavefront radii $R_1$, $R_2$ at the camera. A Gaussian beam carries the transverse phase $k\rho^2/(2R(z))$ (QE3, PDF p. 10). The phase difference across the camera is therefore
  $$\phi(r) = 2\pi\frac{r^2}{P} + \phi_c, \qquad P = \frac{2\lambda}{|1/R_1 - 1/R_2|}$$
  Going out by $P$ in $r^2$ adds one fringe. Rotation adds the same phase across the whole beam, so it changes only $\phi_c$ and the rings slide in or out. A tilt between the beams mainly moves the ring center. Refocusing mainly changes $P$.
- Ring spacing over all photos: median P = 1.074 × 10⁵ px². Per-setting means 1.058–1.094 × 10⁵ px² (±1.7 %).
- The camera never saturates (maximum 239 of 255 in the red channel, the only channel with light).

**Ring finder** (`fringes.py`, steps 1–4).

1. Red channel smoothed with a Gaussian of σ = 4 px against speckle.
2. 170 rays every 1.5° from −75° to 180° (the lower left has no fringes).
3. Along each ray, each dark ring's lowest point is searched within ±0.3 of the local fringe period around a first guess and refined to below a pixel with a parabola through neighboring samples. Shallow minima are dropped.
4. Concentric circles (one center, one radius per ring) are fitted to all minima with a robust loss. The rays are re-centered and the search repeated, three rounds in all.

**Circle-fit quality** (all 161 photos). The single minima lie 2.8–8.2 px (rms, median 4.7 px) from their circles. This is speckle on single points, small next to the 130–240 px between rings. Each ring has 54–170 minima (at most one per ray, median 136), so its fitted radius is known to 0.3–2.3 px (median 0.8 px), and the center to about ±0.6–0.7 px. The fitted center wanders by 3 px (SD) over the run. There is no χ² for these fits because the single minima have no independent error bars. The rms is their measured noise.

---

## Task 3: phase shift against ω

Figures: `presentation_slides/07_phase_extraction.png`, `08_phase_vs_omega.png`, backup `B3_data_table.png`.

**Phase of one photo.** A dark ring is where $\phi = \pi + 2\pi m$ (manual Eq. 6 with cos = −1). With the ring model above, the dark rings satisfy

$$R_m^2 = P\,(m + \eta), \qquad \phi_c = \pi - 2\pi\eta$$

A straight-line fit of $R_m^2$ against the ring number $m$ gives the slope $P$ and the intercept $P\eta$. With two or three rings per photo this line has at most one degree of freedom, so it carries no fit statistic. Its quality check is the circle fit above. This is the manual's Eq. 7 (shift of each minimum over the local fringe period) and Schwartz's $\Delta\phi = 2\pi\,\Delta x/\delta$ (Eq. 5), done in $r^2$, where every ring has the same period. $\eta$ is a ratio of two numbers in px², so no pixel-to-length calibration is needed. $P$ is fitted in every photo, so a change of the ring spacing does not move $\phi_c$. The phase change is $\Delta\phi = \phi_c - \phi_c(\text{rest})$, with the rest value the mean of the 10 rest photos. $\eta$ stays between 0 and 1 for the whole run (0.62 at rest, 0.10 at the top speed), so the ring numbering never shifts. At rest only two rings are in the fringe region, from S = 44 on every photo has three.

**Per setting** (mean of the photos, Δφ error the standard error of the mean, SEM):

| Setting S | Photos | First photo | ω (rad/s) | Δφ (rad) | P (10⁵ px²) |
|---|---|---|---|---|---|
| rest | 10 | 16:21:40 | 0 | 0 (reference) | 1.075 |
| 12.59 | 15 | 16:24:10 | 0.235(1) | 0.39(3) | 1.071 |
| 20.74 | 15 | 16:25:44 | 0.402(2) | 0.64(2) | 1.058 |
| 31.11 | 19 | 16:26:38 | 0.613(3) | 1.02(2) | 1.065 |
| 44.07 | 22 | 16:27:54 | 0.877(5) | 1.59(2) | 1.093 |
| 55.83 | 18 | 16:29:08 | 1.117(6) | 1.984(8) | 1.058 |
| 64.00 | 16 | 16:30:16 | 1.283(7) | 2.375(8) | 1.094 |
| 73.65 | 16 | 16:31:28 | 1.480(8) | 2.700(4) | 1.087 |
| 82.18 | 15 | 16:32:32 | 1.654(9) | 2.958(4) | 1.068 |
| 90.00 | 15 | 16:33:24 | 1.814(9) | 3.285(5) | 1.073 |

ω errors combine the calibration line and ε. The rest SEM (0.02 rad) is common to every row.

### The line through the ten settings

$\Delta\phi = a\omega + b$, fitted by ordinary least squares with equal weights, as Schwartz obtains his slopes (Eq. 7). The regression takes the uncertainty of each point from how far the points lie from the line, $s = \sqrt{\sum(\Delta\phi - \text{line})^2/(10-2)}$, and the slope error is $u_a = s/\sqrt{\sum(\omega - \bar\omega)^2}$. The calibration line's own uncertainty is propagated by redrawing its two coefficients (4000 draws) and added in quadrature.

| Quantity | Value | Unrounded |
|---|---|---|
| Slope a (regression and calibration line) | 1.84 ± 0.02 rad/(rad/s) | 1.8390 ± 0.0242 (regression 0.0239, calibration line 0.0037) |
| Distance of the points from the line, s (per-point uncertainty) | 0.04 rad | 0.0447 |
| R² | 0.9986 | |
| Intercept b | −0.05 ± 0.03 rad, p = 0.11 (t test) | −0.048 ± 0.027 |
| Added ω² term | 0.5 standard errors from zero, p = 0.65 | no sign of a bend |
| Line forced through the origin (check) | 1.80 ± 0.01 | 1.802 ± 0.014 |

**Why a regression and not the photos' error bars.** Weighting each setting by its SEM (0.004–0.028 rad) and the table-speed term gives χ²/dof = 56.6/8, p = 2 × 10⁻⁹ (slope 1.843). The setting means lie 2.7 times farther from the line than those error bars allow. The SEM assumes the photos of one setting are independent, and they are not, because the phase wanders slowly (see the drift section). The regression's $s$ contains every effect that moves a setting mean (photo noise, table speed, wander). A weighted fit with an extra per-setting term σ_x = 0.041 rad, set by χ²/dof = 1, gives the same 1.839 ± 0.024 (`sagnac_update/sagnac.txt`). A regression through all 161 photos would treat correlated photos as independent and give a misleadingly small 1.85 ± 0.01. In a regression χ²/dof = 1 by construction, so R² is the fit-quality number.

Two checks on that extra variation (computed for the presentation):
- Splitting each setting's photos into first and second halves, the halves differ by more than their photo noise allows in 8 of 10 settings (2.5–6 standard errors). The implied wander of a setting mean is 0.039(7) rad, close to $s$ (unrounded 0.045 rad).
- The points' distances from the line follow the mean fringe period P of each setting, r = +0.81, p = 0.005 (the ring-fit and cosine-fit periods agree, r = 0.89). Something changes the beams slightly between settings, and the phase moves with it.

**Comparison with the prediction** (regression error only): a / a_th = 1.06 ± 0.08. Difference 0.101 against a combined uncertainty $\sqrt{0.024^2 + 0.137^2} = 0.139$, so 0.7σ, two-sided p = 0.47. For scale, Schwartz's large ring gave 1.304 ± 0.038 against a prediction of 1.344 ± 0.006 (agreement within 3 %).

**Largest shift.** At ω = 1.81 rad/s, Δφ = 3.29 ± 0.06 rad = 0.52 ± 0.01 fringe. The ± 0.06 is $\sqrt{2}\,s$, top setting minus rest. As a travel-time difference $\Delta t = \Delta\phi\,\lambda/(2\pi c) = 1.12 \pm 0.02$ fs. The prediction there is 3.2 ± 0.2 rad (unrounded 3.15), $\Delta t = 1.07$ fs, half an optical period (2.13 fs).

### Phase-reading methods

Figure: `presentation_slides/09_methods.png`. Each method's per-setting phases go through the same regression.

| Method | What it reads | a (rad/(rad/s)) | s (rad) | R² | Photo-error fit χ²/dof, p |
|---|---|---|---|---|---|
| Ring minima (main) | phase at the ring center, from the ring radii | 1.84 ± 0.02 | 0.04 | 0.9986 | 56.6/8, 2 × 10⁻⁹ |
| Cosine fit of the $r^2$ profile (prep Q5) | phase at the center from a fitted fringe model | 1.88 ± 0.03 | 0.05 | 0.9982 | 131.6/8, 1 × 10⁻²⁴ |
| Template matching to the rest profile | shift of the whole profile, read mid-window | 1.77 ± 0.02 | 0.03 | 0.9993 | 42.7/8, 1 × 10⁻⁶ |

Unrounded slopes 1.8390, 1.8761, 1.7676. All three see the same extra variation between settings (last column).

- **Cosine fit.** Intensity against $s = r^2$: $C e^{-s/s_C} + B e^{-s/s_B}\cos(2\pi s/P + \phi) + D$. The exponentials are the envelopes of two Gaussian beams of different size (QE3, PDF p. 7), D the dark level. It uses every pixel but relies on the envelope model and a linear camera response. Its fits report an rms and no χ², because neighboring pixels share speckle, so the profile points are not independent. It reads 0.10–0.19 rad higher than the minima at every moving setting, because at rest the bright flat center is not a cosine (fit rms 9.5 counts at rest, 1.4–2.3 above S = 44).
- **Template matching.** Each $r^2$ profile is shifted (and scaled) to match the averaged rest profile, and the shift $D$ gives $2\pi D/P$. It needs no fringe model but reads the phase in the middle of the comparison window, where a change of P adds to it.
- **Method uncertainty.** The methods read the same photographs, so they are correlated checks, not independent measurements. Half the range of the three slopes, (1.8761 − 1.7676)/2 = 0.054 rad/(rad/s) (3.0 % of a), is taken as one standard uncertainty. This is a convention. Range/√12 would give 0.031 and the full range 0.108.

### Random-pixel method (1000 draws)

Figure: `presentation_slides/10_random_pixels.png`. Script: `presentation_slides/source/random_points_scan.py` (the method of `sagnac_update/methods.py`, unchanged).

**Model.** Ten pixels at fixed image positions in the fringe sector are read in all 161 photos: $I_{jk} = a_j + b_j\cos(\psi_j + \delta_k)$, with $\psi_j = 2\pi s_j/P_0$ the pixel's phase at rest and $\delta_k$ the photo's phase change. An alternating least-squares fit (pixel parameters with the phases fixed, then the phases with the pixel parameters fixed, repeated until it settles) finds all $a_j$, $b_j$, $\delta_k$, with $b_j \ge 0$ and $a_j - b_j \ge$ dark level. Each draw's phases go through the same averaging per setting and the same regression.

| 1000 draws of ten pixels | Value |
|---|---|
| Mean slope (998 good draws) | 1.702(6) |
| Spread between draws (SD) | 0.196(5). Middle 68 % of the draws: 1.52–1.89 (half-width 0.187) |
| Typical fit error of one draw (median) | 0.057 |
| Spread ÷ fit error | 3.4 |
| Failed fits (negative slope, left out) | 2 |
| Stopped at the 500-round limit (kept) | 18 |
| The single draw of `sagnac_update` (seed 20261011) | 1.79 |

| Pixels per draw | Draws | Mean | Spread (SD) | Fit error of one draw |
|---|---|---|---|---|
| 10 | 1000 | 1.70 | 0.196(5) | 0.057 |
| 20 | 60 | 1.71 | 0.091(8) | 0.041 |
| 40 | 60 | 1.73 | 0.073(6) | 0.040 |
| 80 | 60 | 1.73 | 0.056(5) | 0.037 |

- One draw of ten pixels gives the slope only to about ±0.2 (11 %), 3.4 times worse than its own error bar says. More pixels reduce the spread but it stays above the fit error, and the mean settles near 1.73.
- **Why.** Our largest shift is only half a fringe, so each pixel passes through less than one bright-dark cycle, and the fit has little to separate a phase change from a change of brightness. Prep Q5 asks which method is more accurate when Δφ_max is large: there, fixed points pass through many full cycles. Our run is the small-shift case.
- A fixed pixel also reads a moving ring center or a change of P as a phase change. Recomputing each pixel's $\psi_j$ from the photo's own fitted center and period (diagnostic) moves the mean from 1.70 to 1.78 and narrows the spread to 0.17. So about half of the 0.14 difference from the ring minima comes from the fixed positions.
- The method is not used for the result or the method uncertainty. On slide 9 it is shown as the mean of the draws with one draw's typical fit error (0.06).

### Final slope

$$a = 1.84 \pm 0.07\ \mathrm{rad/(rad/s)}$$

| Part | On a (unrounded) | Relative |
|---|---|---|
| Regression and calibration line | 0.024 | 1.3 % |
| Phase-reading method (half range of three) | 0.054 | 3.0 % |
| Calibration model (linear vs quadratic) | 0.032 | 1.7 % |
| **Added in quadrature** | **0.067** | **3.7 %** |

The parts come from different things (noise between settings, choice of phase reading, choice of calibration curve) and are added in quadrature as independent. A straight sum (0.11, the worst case) or the full method range (total 0.12) would still leave the result 0.6σ from the prediction. `sagnac_update/RESULTS.md` quotes ± 0.06 because it lists the calibration-model change separately.

**Against the prediction** (full uncertainty): a / a_th = 1.06 ± 0.09, difference 0.7σ (z = 0.66), p = 0.51.

---

## Speed of light

Figure: `presentation_slides/11_speed_of_light.png`.

$$c = \frac{8\pi^2 N R^2}{\lambda a} = (2.8 \pm 0.2) \times 10^8\ \mathrm{m/s}$$

Unrounded 2.834 ± 0.246 (8.7 %). c / c₀ = 0.95 ± 0.08, with c₀ = 299 792 458 m/s (defined). Difference −0.7σ (z = −0.67), two-sided p = 0.50. This is the same comparison as a / a_th ($c/c_0 = a_{\mathrm{th}}/a$), shown as a speed of light, not an independent test. The p differs from the slope's 0.51 in the second digit only, because for c the uncertainties are propagated as relative uncertainties (c ∝ 1/a, to first order). It is inferred through the Sagnac formula, not timed.

| Source | Relative uncertainty of c |
|---|---|
| Turn count N = 13 ± 1 | 7.7 % |
| Phase-reading method | 3.0 % |
| Loop radius R = 57.0 ± 0.5 cm (enters squared) | 1.8 % |
| Calibration model | 1.7 % |
| Regression and calibration line | 1.3 % |
| **Total of the listed sources** | **8.7 %** |
| Wavelength (not given, not included) | 0.16 % per nm |

- Our measurement (regression, calibration line, calibration model, method) is 3.7 %. The given geometry is the rest.
- Read the other way, with c, R and λ taken as known, the slope fixes the number of turns: N = 13.8 ± 0.5 (± 0.6 with R's uncertainty), about twice as precise as the given 13 ± 1.
- If the given ±1 turn is a bound rather than one standard deviation (flat distribution, ± 0.58 turns), the uncertainty of c drops to about 6 %.
- With the laser frequency ν known instead of λ, $c = \sqrt{8\pi^2 N R^2 \nu / a}$ would have half the relative uncertainty.

---

## Task 5: uncertainties

Figure: `presentation_slides/B2_error_budget.png`.

| Source | Kind | Size | Included |
|---|---|---|---|
| Photo-to-photo variation | statistical | SEM 0.004–0.028 rad per setting | yes, through the regression |
| Table speed at a setting | statistical (calibration) | ε = 0.48 % of ω | yes, through the regression |
| Wander between settings | statistical | about 0.04 rad per setting | yes, through the regression |
| Calibration line | calibration, correlated across settings | 0.2 % of a | yes |
| Phase-reading method | systematic | 3.0 % of a | yes |
| Linear vs quadratic calibration | model choice | 1.7 % of a | yes |
| Turn count N = 13 ± 1 | given | 7.7 % | yes (prediction, c) |
| Loop radius R = 57.0 ± 0.5 cm | given | 1.8 % | yes (prediction, c) |
| Wavelength λ = 640 nm | not given | 0.16 % per nm | no |
| Steady drift in time | systematic | not separable from a | no |

Not checked: the sign of the shift (the rotation direction was not recorded). ω comes from the calibration at the same setting, not from timing during the photographs (ε covers that).

### Drift in time

Figure: `presentation_slides/B4_drift.png`.

The settings were photographed in increasing order, with rest photos only at the start, so ω and time rise together. A phase that drifts steadily in time would change the slope and cannot be told apart from it.

- **Joint regression** $\Delta\phi = a\omega + Dt + b$ (t the mean time of each setting's photos): a = 2.1 ± 0.2 rad/(rad/s), D = −0.8 ± 0.5 mrad/s, R² = 0.9991. The drift term is not significant (t test, p = 0.12), and the correlation between a and D is −0.99. With the first photo's time instead, D = −1.1 ± 0.5, p = 0.07. A negative D, if real, would mean the slope without it is too small.
- **Within each setting** the phase drifts by −6.9 to +8.2 mrad/s, significant (p < 0.05) in 6 of 10 settings, and changing sign between settings. Their plain mean, −1 ± 1 mrad/s, cannot rule a steady drift in or out. (A mean weighted by their own error bars would be misleading, because the ten drifts differ far more than those error bars allow.)
- The quoted uncertainties do not include a steady drift.

**A plausible source of the wander (not tested).** At 640 nm SMF-28e+ guides several modes: the manual gets V = 5.6 (core radius 4.1 µm, NA 0.14), above the single-mode limit 2.405 (prep Q12). Light in several modes is sensitive to fiber stress and temperature, and the camera records the whole multimode field with no single-mode filter at the common port, so the two directions need not stay exactly reciprocal. How much light enters the fundamental mode depends on the coupling: the manual computes up to 93 % with a matched 15.5 mm lens and at most 47 % with the setup's 8 mm lens. Coupling into one mode is possible with this fiber but was very hard to align.

### Next time

- Take rest photos between the speeds and reverse the rotation. The Sagnac phase changes sign with the direction of rotation, a drift does not (Sagnac measured this way, Pascoli Sec. 2). Take the speeds in random order.
- Count the fiber turns. N = 13 ± 1 is the largest uncertainty.
- Couple the light into a single fiber mode, with a longer-focal-length coupling lens.
- Record the rotation direction.

---

## Files

| File | Content |
|---|---|
| `SUMMARY.md` | short version of this document |
| `presentation_slides/` | the slides (PNG and GIF), `presenter_notes.md`, `presenter_guide.md`, `README.md` |
| `presentation_slides/source/collect_numbers.py`, `numbers.json` | every number in this document |
| `presentation_slides/source/slides.py`, `deck.py` | the presentation analysis and figures |
| `sagnac_update/` | linear-calibration outputs the presentation builds on |
| `fringes.py`, `fringes_frames.csv` | ring finder and per-photo results |
| `calibration.py`, `calibration.json`, `calibration.txt` | rotation calibration (linear and quadratic fits) |
| `previous_reports/` | the earlier `RESULTS.md` and `SUMMARY.md` (quadratic calibration) |
