# Sagnac effect: results (11 October 2026 data)

The final analysis, ready for the report. Every number is computed by `calibration.py` and `sagnac.py` (`calibration.txt`, `sagnac.txt`, `sagnac.json`, `sagnac_speeds.csv`); the figures are in `report_figures/`. The short version is `SUMMARY.md`; the uncertainty derivations are in `uncertainty/REPORT.md`.

**Conventions.** All uncertainties are standard uncertainties (one standard deviation). Values are given rounded, with the uncertainty to its first significant digit and the value to the same digit; comparisons (σ, p) use the unrounded values, which are listed where needed. A regression reports R² and the rms distance of the points from the line; a fit against given error bars reports χ²/dof and p. Figures show data and fits only, no residual panels.

---

## Setup and given parameters

| Quantity | Value | Source |
|---|---|---|
| Fiber turns $N$ | 13 ± 1 | given (manual, Fig. 1 caption) |
| Loop radius $R$ | 57.0 ± 0.5 cm | given (manual, Fig. 1 caption) |
| Wavelength $\lambda$ | 640 nm, no uncertainty given | given |
| Fiber | Corning SMF-28e+ (the manual's text says SMF-28e, its datasheet and the preparation sheet say SMF-28e+). Core radius 4.1 µm, NA 0.14 | datasheet in the manual |
| Enclosed area of one turn | $A = \pi R^2 = 1.02\ \mathrm{m^2}$ | |

- The loop is centered on the rotation axis. Everything on the table (laser, polarizers, beam splitter, coupling lenses, fiber, camera) turns with it.
- Data: 161 photographs at rest and at nine motor settings, 10–22 per setting, taken in increasing order of speed between 16:21:40 and 16:33:46. Rotation timing at nine settings.
- Manual tasks done: 1, 2, 3 and 5. Task 4 (the off-center small loop, R = 15.5 ± 0.1 cm, 46 ± 1 turns) was skipped on the instructor's instruction.

**Prediction.** The Sagnac phase (manual Eq. 5) for $N$ turns enclosing $A$ each is

$$\Delta\phi = \frac{8\pi N A}{\lambda c}\,\omega \equiv a_{\mathrm{th}}\,\omega, \qquad a_{\mathrm{th}} = \frac{8\pi^2 N R^2}{\lambda c} = 1.7 \pm 0.1\ \mathrm{rad/(rad/s)}$$

Unrounded $a_{\mathrm{th}} = 1.738 \pm 0.137$: $N$ gives 7.7 %, $R$ (squared) 2 × 0.88 % = 1.8 %, together 7.9 %. The unit rad/(rad/s) is a second. The prediction does not depend on where the rotation axis lies (Schwartz 2017, Eq. 4; prep Q8) nor on the fiber's refractive index (Schwartz, Sec. 2; prep Q10).

---

## Task 1: rotation calibration

Figures: `report_figures/calibration` (the timed points, the line and its ±1σ prediction band) and `report_figures/conversion` (the nine photo settings read off that line, with the prediction uncertainty at the line: the ω column of Table 2).

**Timing.** At each of nine controller settings $S$ (a number on the motor controller, no unit, which sets the motor voltage) n = 5 (the five lowest settings) or 10 (the four highest) consecutive turns were timed with a stopwatch, one press per turn. $T$ is the mean turn time and $\omega = 2\pi/T$. Two things limit $T$: the scatter of the single turns, SD/√n, and the human reaction time of the presses. The n turns are consecutive laps, so the mean turn time is $(t_n - t_0)/n$ and only the first and the last press enter it; with a reaction time $\sigma_r = 0.25$ s per press that part is $\sqrt{2}\sigma_r/n$. Together

$$u(T)^2 = \left(\frac{\mathrm{SD}}{\sqrt n}\right)^2 + \left(\frac{\sqrt2\,\sigma_r}{n}\right)^2, \qquad \frac{u(\omega)}{\omega} = \frac{u(T)}{T}$$

which gives 0.15–1.1 % of ω (the reaction time dominates at the fast settings). Example, S = 90: 10 turns, ⟨T⟩ = 3.545 s, SD = 0.048 s, SD/√n = 0.015 s, √2σ_r/n = 0.035 s, u(T) = 0.038 s, ω = 1.772 ± 0.019 rad/s. The 0.25 s is our own reaction time, measured with an online reaction-time test, used as the uncertainty of each press; the lap scatter itself implies a smaller press-to-press jitter, so it is on the safe side.

**The line.** $\omega = kS + \omega_0$, fitted with each timed point weighted by its uncertainty; the uncertainties of k and ω₀ follow from the timing uncertainties:

| | Value |
|---|---|
| $k$ | 0.02054 ± 0.00006 rad/s per unit S (0.29 %) |
| $\omega_0$ | −0.0227 ± 0.0004 rad/s |
| χ²/dof, p | 20.2/7, 0.005 |

The points lie farther from the line than their uncertainties allow (χ²/dof = 2.9). The deviations from the line are −0.1, +0.9, +0.3, +0.1, +0.5, 0.0, −0.9, −1.9, −3.0 % in order of S (rms 1.3 %): the table turns slower than linear at high settings. At S = 90 the line gives 1.826 rad/s, 3.0 % above the timed 1.772 rad/s. The intercept is not physical (the line was only tested for S = 5.78–90); at rest ω = 0 is used.

**ω at the photo settings.** A photo setting was not timed: the line predicts its ω. The prediction carries the line's own uncertainty, $(1, S)\,\mathrm{Cov}\,(1, S)^T$ (0.2–0.3 % of ω), and the scatter of the timed settings about the line, 1.3 % of ω (the table does not repeat a setting better than that); in quadrature 1.3 % of ω (Table 2). Both parts are common to every photo setting (the same line, the same table), so the calibration part of the Sagnac slope's uncertainty has the same relative size, $\sqrt{0.29^2 + 1.28^2}\,\% = 1.3\ \%$ (Task 3).

---

## Task 2: fringe patterns

Figures: `report_figures/patterns` (one photo per setting), `report_figures/rings`.

- The pattern is a set of rings. As ω rises, the rings contract: the first dark ring at rest (radius 258 px) moves in to 101 px at the top speed, and the center goes from bright to nearly dark, about half a fringe.
- The rings are evenly spaced in $r^2$, not in $r$. In the fastest photo the dark rings are 241 px and then 132 px apart in $r$, but 1.08 × 10⁵ px² both times in $r^2$.
- **Checked over the run** (23 October, from the radii R1–R3 in `fringes_frames.csv`): in the 126 photos with three dark rings the two gaps in $r^2$, $R_2^2 - R_1^2$ and $R_3^2 - R_2^2$, agree to 0.1 % on average, with a spread of 2 % of a fringe (SD of their difference over P; largest 7.5 %), while in $r$ the ratio of the two gaps is 1.35–1.86 (mean 1.61).
- **Why the manual's Eq. 7 is carried out in $r^2$.** Eq. 7 reads a phase as the shift of a minimum over its local period, Δr_i/Λ_i, and Schwartz (p. 6, Eq. 5) reads his circular fringes the same way, with an estimated width. That is exact only while the width is the same over the shift. Ours reach half a fringe and the width changes under the moving ring: the first dark ring moves from 258 to 101 px (median radii of the rest and top-speed photos), 157 px, while its gap to the next ring is 159 px at rest and 242 px at the top speed, so Eq. 7 in $r$ gives 6.2, 4.1 or 4.9 rad with the rest width, the top-speed width or their mean, against 3.30 rad in $r^2$, (258² − 101²)/P with P = 1.074 × 10⁵ px², and the report's 3.28 rad. In $r^2$ the period is one number for every ring and every photo, so Eq. 7 is exact there.
- **Why.** The two beams leave the fiber as Gaussian beams with different wavefront radii $R_1$, $R_2$ at the camera. A Gaussian beam carries the transverse phase $k\rho^2/(2R(z))$ (QE3, PDF p. 10), so the phase difference across the camera is
  $$\phi(r) = 2\pi\frac{r^2}{P} + \phi_c, \qquad P = \frac{2\lambda}{|1/R_1 - 1/R_2|}$$
  Going out by $P$ in $r^2$ adds one fringe. Rotation adds the same phase across the whole beam, so it changes only $\phi_c$: the rings slide in or out. A tilt between the beams mainly moves the ring center; refocusing mainly changes $P$ (prep Q7).
- Ring spacing over all photos: median P = 1.074 × 10⁵ px²; per-setting medians 1.058–1.094 × 10⁵ px² (±1.7 %).
- The camera never saturates (maximum 239 of 255 in the red channel, the only channel with light).

**Ring finder** (`fringes.py`): the red channel smoothed with a Gaussian of σ = 4 px against speckle; 170 rays every 1.5° from −75° to 180° (the lower left has no fringes); along each ray the lowest point of each dark ring, refined to below a pixel with a parabola; concentric circles (one center, one radius per ring) fitted to all minima with a robust loss; three rounds, each re-centered. The single minima lie 2.8–8.2 px (rms, median 4.7 px) from their circles, which is speckle on single points, small next to the 130–240 px between rings; with dozens to 170 minima per ring, the radii and the center are fixed to about a pixel. The center wanders by 3 px (SD) over the run.

---

## Task 3: phase shift against ω

Figures: `report_figures/extraction`, `report_figures/phase_vs_omega`.

**Phase of one photo.** A dark ring is where $\phi = \pi + 2\pi m$ (manual Eq. 6 with cos = −1), so with the ring model above

$$R_m^2 = P\,(m + \eta), \qquad \phi_c = \pi - 2\pi\eta$$

A straight line through $R_m^2$ against the ring number $m$ gives the slope $P$ and the intercept $P\eta$. This is the manual's Eq. 7 (shift of each minimum over the local fringe period, averaged over the minima) and Schwartz's $\Delta\phi = 2\pi\,\Delta x/\delta$ (Eq. 5), done in $r^2$, where every ring has the same period. $\eta$ is a ratio of two numbers in px², so no pixel-to-length calibration is needed, and $P$ is fitted in every photo, so a change of the ring spacing does not move $\phi_c$. With two or three rings per photo the line has at most one degree of freedom and no fit statistic; its quality check is the circle fit above. $\Delta\phi = \phi_c - \phi_c(\text{rest})$, with the rest value the mean of the 10 rest photos. $\eta$ stays between 0 and 1 for the whole run (0.62 at rest, 0.10 at the top speed), so the ring numbering never shifts. At rest two rings are in the fringe region; from S = 44 on every photo has three.

**Per setting** (Table 2; mean of the photos. The SEM column is the photo noise; the uncertainty of a setting mean in the report is s = 0.044 rad, the scatter about the line, the same for every setting: decision of 24 October):

| Setting S | Photos | First photo | ω (rad/s) | Δφ (rad) | P (10⁵ px²) |
|---|---|---|---|---|---|
| rest | 10 | 16:21:40 | 0 | 0 (reference) | 1.075 |
| 12.59 | 15 | 16:24:10 | 0.236 ± 0.003 | 0.39 ± 0.03 | 1.071 |
| 20.74 | 15 | 16:25:44 | 0.403 ± 0.005 | 0.64 ± 0.02 | 1.058 |
| 31.11 | 19 | 16:26:38 | 0.616 ± 0.008 | 1.02 ± 0.02 | 1.065 |
| 44.07 | 22 | 16:27:54 | 0.88 ± 0.01 | 1.59 ± 0.02 | 1.093 |
| 55.83 | 18 | 16:29:08 | 1.12 ± 0.01 | 1.984 ± 0.008 | 1.058 |
| 64.00 | 16 | 16:30:16 | 1.29 ± 0.02 | 2.375 ± 0.008 | 1.094 |
| 73.65 | 16 | 16:31:28 | 1.49 ± 0.02 | 2.700 ± 0.004 | 1.087 |
| 82.18 | 15 | 16:32:32 | 1.67 ± 0.02 | 2.958 ± 0.004 | 1.068 |
| 90.00 | 15 | 16:33:24 | 1.83 ± 0.02 | 3.285 ± 0.005 | 1.073 |

ω from the calibration line with its prediction uncertainty (1.3 %, common to all settings). The rest SEM (0.02 rad) is the photo noise of the reference; in the report every setting mean, the rest one included, carries ± 0.04 rad (s). P is the median of the setting's photos.

### The line through the ten settings

$\Delta\phi = a\omega + b$ by ordinary least squares, ten points, equal weights, as Schwartz obtains his slopes (Eq. 7). The regression takes the uncertainty of every point from the points themselves: $s = \sqrt{\sum(\Delta\phi - \text{line})^2/(10-2)}$, the rms distance of the points from the line, and $u(a) = s/\sqrt{\sum(\omega - \bar\omega)^2}$.

| Quantity | Value | Unrounded |
|---|---|---|
| Slope a (regression) | 1.83 ± 0.02 rad/(rad/s) | 1.8258 ± 0.0236 |
| s, the rms distance of the points from the line | 0.04 rad | 0.0444 |
| R² | 0.9987 | |
| Intercept b | −0.05 ± 0.03 rad, p = 0.12 (t test) | −0.047 ± 0.027 |
| Added ω² term | 0.5 standard errors from zero, p = 0.67 | no bend visible; the calibration's bend alone would give −0.04, inside the ± 0.05, so the check cannot tell the bend from none |
| Line forced through the origin (check) | 1.79 ± 0.01 | 1.791 ± 0.014 |

**Why the uncertainty of a setting mean is s = 0.044 rad, not the photos' standard error** (presented so in the report from Table 2 on, 24 October). The fitted line misses the setting means by five times their SEM (0.004–0.028 rad; χ²/dof = 205/8, p = 6 × 10⁻⁴⁰). The SEM assumes the photos of one setting are independent, and they are not: the phase wanders slowly (see the drift section), so a setting mean is less well known than its SEM says. The regression's $s$ contains everything that moves a setting mean (photo noise, the table speed, the wander), which is why it is used. A regression through all 161 photos would treat correlated photos as independent and give a misleadingly small uncertainty.

**Comparison with the prediction** (regression uncertainty only): a / a_th = 1.05 ± 0.08; the difference 0.088 against $\sqrt{0.024^2 + 0.137^2} = 0.139$ is 0.6σ, two-sided p = 0.53. For scale, Schwartz's large ring gave 1.304 ± 0.038 against a prediction of 1.344 ± 0.006 (agreement within 3 %).

**Largest shift.** At ω = 1.83 rad/s, Δφ = 3.28 ± 0.06 rad = 0.52 ± 0.01 fringe (the ± 0.06 is $\sqrt{2}\,s$: top setting minus rest, each known to s). As a travel-time difference $\Delta t = \Delta\phi\,\lambda/(2\pi c) = 1.12 \pm 0.02$ fs. The prediction there is 3.2 ± 0.3 rad (unrounded 3.17), $\Delta t = 1.08$ fs, half an optical period (2.13 fs).

### Phase-reading methods (prep Q5)

Figure: `report_figures/methods`. Each method's per-setting phases go through the same regression.

| Method | What it reads | a (rad/(rad/s)) | s (rad) | R² |
|---|---|---|---|---|
| Ring minima (main) | phase at the ring center, from the ring radii | 1.83 ± 0.02 | 0.044 | 0.9987 |
| Cosine fit of the $r^2$ profile | phase at the center from a fitted fringe model | 1.86 ± 0.03 | 0.053 | 0.9982 |
| Template matching to the rest profile | shift of the whole profile, read mid-window | 1.75 ± 0.02 | 0.031 | 0.9993 |
| Fixed pixels (main draw of ten; in Fig. 7, not in the method range) | the brightness of ten fixed pixels, one phase per photo | 1.8 ± 0.2 (spread between draws; the fit alone 1.77 ± 0.04) | 0.069 | 0.9966 |

Unrounded slopes 1.8258, 1.8626, 1.7549; fixed pixels 1.7734 (not in the range).

- **Cosine fit.** Intensity against $s = r^2$: $C e^{-s/s_C} + B e^{-s/s_B}\cos(2\pi s/P + \phi) + D$. The exponentials are the envelopes of two Gaussian beams of different size (QE3, PDF p. 7), D a floor. It uses every pixel but relies on the envelope model and a linear camera response. Its fits report an rms and no χ², because neighboring pixels share speckle. It reads 0.10–0.19 rad higher than the minima at every moving setting, because at rest the bright flat center is not a cosine (fit rms 9.5 counts at rest, 1.4–2.3 above S = 44).
- **Template matching** (our addition). Each $r^2$ profile is shifted and scaled to match the averaged rest profile; the shift $D$ gives $2\pi D/P$. It needs no fringe model but reads the phase in the middle of the comparison window, where a change of P adds to it.
- **Method uncertainty.** The methods read the same photographs, so they are correlated checks, not independent measurements. Half the range of the three slopes, (1.8626 − 1.7549)/2 = 0.054 rad/(rad/s) (2.9 % of a), is taken as one standard uncertainty. This is a convention; it is the largest part of the slope's uncertainty.
- **Fixed pixels** (the second half of prep Q5; `pixels.py`, 24 October; in Fig. 7 beside the other readouts, outside the method range; the instructor asked for this readout). Ten pixels at fixed image positions are read in all 161 photos and fitted with $I = a_j + b_j\cos(\psi_j + \delta_k)$ for every pixel j and photo k, ψ_j the pixel's fringe phase at rest, under two physical limits: b_j ≥ 0, and a_j − b_j ≥ the camera's dark level (2.0 counts). Without them, ten pixels that see only half a fringe cannot tell a large phase change at low contrast from a small one at high contrast. The main draw (seed fixed in advance) gives, through the same averaging and regression, a = 1.7734 ± 0.0367 (fit), s = 0.069 rad, R² = 0.9966, and lies 0 to 0.12 rad below the minima per setting. Over 1001 random draws of ten pixels (6 failed with a negative slope) the slope has mean 1.689 and spreads by 0.193 (SD, 11 %) between draws, 3.4 times the median fit error of one draw (0.057); with 20, 40 and 80 pixels the spread falls to 0.090, 0.072 and 0.055 (5.3, 4.2, 3.2 %) but stays above the fit error. Our largest shift is only half a fringe, so each pixel passes through less than one bright–dark cycle, and a fixed pixel also reads a moving ring center or a change of P as a phase change. For a large Δφ_max, where fixed points pass through many cycles, this readout would use the whole pattern and could be the most precise; our run is the small-shift case, so it is shown with the spread between draws as its uncertainty and left out of the method range. (The presentation stage's scan with the earlier calibration line gave the same: mean 1.70, SD 0.20.)

### Final slope

$$a = 1.83 \pm 0.06\ \mathrm{rad/(rad/s)}$$

| Part | On a (unrounded) | Relative | How |
|---|---|---|---|
| Regression | 0.024 | 1.3 % | s and the spread of the ten ω |
| Phase-reading method | 0.054 | 2.9 % | half the range of the three methods |
| Calibration | 0.024 | 1.3 % | the line's prediction uncertainty, common to every ω: $\sqrt{(u_k/k)^2 + (1.28\,\%)^2}$ |
| **Added in quadrature** | **0.064** | **3.5 %** | |

The three parts come from different things (noise between settings, how the phase is read, the calibration line) and are added as independent.

**Against the prediction** (full uncertainty): a / a_th = 1.05 ± 0.09, difference 0.6σ (z = 0.58), p = 0.56. Read the other way, with c, R and λ taken as known, the slope fixes the number of turns: N = 13.7 ± 0.5, against the given 13 ± 1.

---

## Speed of light

$$c = \frac{8\pi^2 N R^2}{\lambda a} = (2.9 \pm 0.2) \times 10^8\ \mathrm{m/s}$$

Unrounded 2.854 ± 0.246 (8.6 %). c / c₀ = 0.95 ± 0.08, with c₀ = 299 792 458 m/s (defined); −0.6σ (z = −0.59), two-sided p = 0.56. This is the same comparison as a / a_th ($c/c_0 = a_{\mathrm{th}}/a$), shown as a speed of light, not an independent test: the slope is inferred through the Sagnac formula with the given geometry.

| Source | Relative uncertainty of c |
|---|---|
| Turn count N = 13 ± 1 | 7.7 % |
| Phase-reading method | 2.9 % |
| Loop radius R = 57.0 ± 0.5 cm (enters squared) | 1.8 % |
| Calibration | 1.3 % |
| Regression | 1.3 % |
| **Total (in quadrature)** | **8.6 %** |
| Wavelength (not given, not included) | 0.16 % per nm |

Our measurement (regression, method, calibration) is 3.5 %; the given geometry is the rest. The phase fixes the product λc, so c needs the wavelength; with the laser frequency ν known instead, $c = \sqrt{8\pi^2 N R^2 \nu/a}$ would have half the relative uncertainty.

---

## Task 5: uncertainties

Figure: `report_figures/drift`. Derivations: `uncertainty/REPORT.md`.

| Source | Kind | Size | Included |
|---|---|---|---|
| Photo-to-photo variation | statistical | SEM 0.004–0.028 rad per setting | yes, in s |
| Wander of the phase between settings | statistical | about 0.04 rad per setting | yes, in s (it dominates s) |
| Stopwatch: lap scatter and reaction time (0.25 s per press) | statistical | 0.15–1.1 % of a timed ω | yes, through the line |
| Calibration line and the table's scatter about it | statistical (correlated across settings) | 1.3 % of a | yes |
| Phase-reading method | systematic | 2.9 % of a | yes |
| Turn count N = 13 ± 1 | given | 7.7 % | yes (prediction, c) |
| Loop radius R = 57.0 ± 0.5 cm | given | 1.8 % | yes (prediction, c) |
| Wavelength λ = 640 nm | not given | 0.16 % per nm | no |
| Bend of the calibration at high settings (line 3.0 % above the timed ω at S = 90) | systematic | up to 3 % of ω at the top settings; +2.2 to +2.5 % on a if a curve replaced the line | partly: its rms is in the 1.3 % scatter, the 3.0 % at the top setting is not |
| Steady drift of the phase in time | systematic | not separable from a | no |

Not checked: the sign of the shift (the rotation direction was not recorded). ω comes from the calibration at the same setting, not from timing during the photographs. The calibration's bend at high settings would raise the slope if the line were replaced by a curve through the timed points: by 2.2 % with a weighted quadratic through them and 2.5 % with straight segments between them, 0.04 rad/(rad/s), because the line's ω is too high at the fast settings (`bend` block of `sagnac.txt`, 24 October; "lower" here before was a sign error). Its average size is in the scatter term, its full size is not. An ω² term in Δφ(ω) cannot tell this bend from none: the bend alone would give −0.04 rad/(rad/s)² against the fitted +0.02 ± 0.05.

### Drift in time

The settings were photographed in increasing order, with rest photos only at the start, so ω and time rise together. A phase that drifts steadily in time would change the slope and cannot be told apart from it. The joint regression $\Delta\phi = a\omega + Dt + b$ (t the mean time of each setting's photos) gives a = 2.1 ± 0.2 rad/(rad/s), D = −0.8 ± 0.5 mrad/s, p = 0.12 for D, and a correlation of −0.99 between a and D. Within the settings the phase drifts by up to ±8 mrad/s, changing sign (Fig. 8). The quoted uncertainties do not include a steady drift.

**A plausible source of the wander (not tested).** At 640 nm SMF-28e+ guides several modes: V = 5.6 (core radius 4.1 µm, NA 0.14), above the single-mode limit 2.405 (prep Q12). Light in several modes is sensitive to fiber stress and temperature, and the camera records the whole multimode field with no single-mode filter at the common port, so the two directions need not stay exactly reciprocal. Single-mode coupling is possible with this fiber (the manual's 15.5 mm lens gives up to 93 %, the setup's 8 mm lens at most 47 %) but was very hard to align.

### Next time

- Take rest photos between the speeds and reverse the rotation: the Sagnac phase changes sign with the direction of rotation, a drift does not. Take the speeds in random order.
- Count the fiber turns. N = 13 ± 1 is the largest uncertainty.
- Time the turns at the photo settings themselves, with more turns at the fast settings, where the reaction time dominates.
- Couple the light into a single fiber mode, with a longer-focal-length coupling lens.
- Record the rotation direction.

---

## Files

| File | Content |
|---|---|
| `SUMMARY.md` | short version of this document |
| `calibration.py` → `calibration.json`, `calibration.txt` | task 1 |
| `fringes.py` → `fringes_frames.csv`, `fringes_profiles.npz` | ring finder and per-photo phases (minima, cosine fit) |
| `methods.py` → `methods_frames.csv` | template matching per photo |
| `pixels.py` → `pixels_frames.csv`, `pixels.json` | the fixed-pixel readout: the main draw per photo, and the spread of the slope over random draws |
| `sagnac.py` → `sagnac_speeds.csv`, `sagnac.txt`, `sagnac.json` | tasks 3 and 5, c; every number above |
| `report_figures.py` → `report_figures/` | the figures |
| `REPORT_figures_tables.tex`, `REPORT_TODO.md`, `PLAN_CHECKLIST.md` | report blocks, plan, manual checklist |
| `uncertainty/REPORT.md` | the uncertainty derivations |
| `presentation_slides/`, `sagnac_update/`, `previous_reports/` | the presentation of 19 October and the earlier stages (see `README.md`) |
