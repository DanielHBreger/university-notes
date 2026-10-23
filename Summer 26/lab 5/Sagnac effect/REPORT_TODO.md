# Lab 5 report (Sagnac): plan and open items

Started 22 October 2026, before the first draft; revised after the review of the analysis and the feedback on the presentation (no residual panels in figures; reaction time in the stopwatch uncertainty; no quadratic check for now; the calibration file completed on 23 October with 10 turns at the four fastest settings). As in `chaotic circuit/REPORT_TODO.md`: tick an item (`[x]`) or delete it when a draft covers it, and put each draft review at the top of the file, newest first. Numbers here are the final analysis and match `RESULTS.md`, `sagnac.txt` and `calibration.txt`.

What goes with this file:

| File | What it is |
|---|---|
| `report_figures/` (`report_figures.py`) | the report's figures, PDF + 600 dpi PNG, in the report style, data and fits only |
| `REPORT_figures_tables.tex` | the figure and table blocks (numbers and `\includegraphics` only, captions left to you) |
| `PLAN_CHECKLIST.md` | the manual's tasks and the preparation questions against what was done |
| `uncertainty/REPORT.md` | every uncertainty component, its model and kind; material for task 5 and an appendix |
| `docs/` | the code explained script by script (for writing 1.2 and for questions) |

## Rules for the report

- Every number with its uncertainty and what kind it is (statistical, systematic, given). Uncertainty rounded to its first significant figure and the value to the same digit (1.7 ± 0.1, not 1.74 ± 0.14), but every comparison (σ, p) computed from the unrounded values, given below where needed.
- A regression reports R² and s (the rms distance of the points from the line); a fit against given error bars reports χ²/dof and p. Say once (2.3) why the main line is a regression.
- Figures from `report_figures/` (column width unless figure*): data and fits only, no residual panels, no fit results inside the figures; the numbers go in captions and text.
- **The prose is yours.** This file gives content (facts, numbers, order), not sentences; the figure and table blocks hold no text. (The Chua 2.4/2.5 lesson: suggested text pasted almost verbatim had to be rewritten before submission.)
- Methods and terms from the lab documents: the manual, Schwartz 2017, the QE3 Gaussian-beam slides, the preparation sheet. Template matching and the fixed-pixel method are our own additions: say so where they appear. Use only terms you can explain.
- Do not call fit statistics "scatter" or "curvature": "the points lie farther from the line than their error bars allow", "s, the rms distance of the points from the line", "an added ω² term".
- American spelling (center, fiber), as in `RESULTS.md`. Siunitx as in the Chua report (`\qty`, `\num` with `\pm`, compact output).
- One unit for the slope everywhere: rad/(rad s⁻¹) (= s). The figures use rad/(rad s⁻¹); `REPORT_figures_tables.tex` defines `\slopeunit` once.
- ω always from the linear calibration (the top speed is 1.83 rad/s). The old `figures/` folder, `previous_reports/` and `presentation_slides/source/numbers.json` belong to earlier stages: do not take numbers from them.

## BLOCKING before submission

- [ ] All prose written in your own words (abstract, every section, captions).
- [ ] Every number checked against `RESULTS.md` (or `sagnac.txt`); the slope is 1.83 ± 0.06 with the full uncertainty and 1.83 ± 0.02 with the regression only: say which one each comparison uses.
- [ ] Acknowledgments with the AI disclosure (section "Acknowledgments and AI disclosure" below), and the in-text mentions you decide on.
- [ ] References and citations: the manual, the preparation sheet, Schwartz, QE3 at least (list below).
- [ ] Task 4 (off-center small loop) stated as not done, on the instructor's instruction.
- [ ] Every figure and table cited in the text; each figure* defined before the text that cites it.

## Structure (proposed; the Chua report's layout)

| Section | Content | Figures, tables |
|---|---|---|
| Abstract | one paragraph, 150–200 words, no references to figures or sections | |
| 1.1 Theory | 1.1.1 the Sagnac phase; 1.1.2 two-beam interference and the ring pattern; 1.1.3 fiber modes (short) | Fig. 4 can be cited from 1.1.2 |
| 1.2 Setup and procedure | apparatus, given parameters, timing, photographs, the analysis pipeline | Fig. 1 |
| 2.1 Rotation calibration (task 1) | timing with its uncertainty, the line, the conversion | Fig. 2, Fig. 2b (optional), Table 1 (optional) |
| 2.2 Interference patterns (task 2) | the rings and what rotation does to them; equal spacing in r² | Fig. 3 (figure*), Fig. 4 |
| 2.3 Phase shift and comparison with theory (task 3) | phase of one photo; setting means; the regression; comparison; largest shift | Fig. 5, Fig. 6, Table 2 |
| 2.4 Phase-reading methods (prep Q5) | three full-image methods; the method uncertainty; fixed pixels (optional) | Fig. 7, Table 3 (optional), random pixels (optional) |
| 2.5 Speed of light (extension) | c from the slope; budget; N read the other way | Table 4 (optional) |
| 2.6 Uncertainties and drift (task 5) | budget; the final slope; drift; what is not covered; next time | Table 5 (figure-width), Fig. 8 |
| 3 Conclusions | two short paragraphs | |
| Acknowledgments, References | | |
| Appendix A (optional) | uncertainty propagation, from `uncertainty/REPORT.md` | |

Results and Discussion merged, as in the Chua report (the discussion items below say where each one goes).

## Figures and tables

Numbering as in the structure above, assuming every optional table is used (LaTeX renumbers by order of appearance, so refer to them by label). Each caption must say the items listed; write them in your words. The same lists are in `REPORT_figures_tables.tex`.

| | File | Width | Section | The caption must say |
|---|---|---|---|---|
| Fig. 1 | `setup.pdf` | column | 1.2 | top view, after Fig. 1 of the manual (cite); P₁, P₂ polarizers, L₁, L₂ coupling lenses, L₃ camera lens; beams 1 and 2 go around the same loop (N turns, radius R) in opposite directions and meet at the splitter; everything turns with the table about the loop's center; the ω arrow is schematic (direction not recorded); the small loop of task 4 is not drawn |
| Fig. 2 | `calibration.pdf` | column | 2.1 | ω = 2π/⟨T⟩ from 5 or 10 consecutive timed turns against the controller setting S; error bars: lap scatter and the 0.2 s reaction time of the presses, mostly smaller than the markers; line: the weighted least-squares fit used for the photos; band: its ±1σ prediction uncertainty (the line's own and the 1.3 % rms scatter of the timed settings about it), the uncertainty of ω at any photo setting |
| Fig. 2b (optional) | `conversion.pdf` | column | 2.1 or 2.3 | how the nine photo settings become ω: dotted guides from each setting up to the line of Fig. 2 (same line and band) and across to the ω axis, the bar at the line the prediction uncertainty; nothing was measured at these settings; these are the ω values of Table 2, so use the figure or the table |
| Fig. 3 | `patterns.pdf` | figure* | 2.2 | one photo per setting in the order taken, the one with the median phase of its setting; red channel in gray, cropped around the fitted center (+); dashed: fitted dark rings; ω from the linear calibration; no fringes in the lower left |
| Fig. 4 | `rings.pdf` | column | 2.2 | intensity averaged over the fringe sector (−75° to 180°) in 1 px rings about each photo's fitted center, rest and top speed; (a) against r, (b) against r²; dashed: fitted dark-ring radii; arrows: the top-speed spacings, unequal in r, equal in r² |
| Fig. 5 | `extraction.pdf` | column | 2.3 | (a) top-speed photo: darkest point of each ring on each ray (rays every 1.5°, one color per ring) and the concentric circles fitted to them; (b) R_m²/P against m at rest and at top speed, the line R_m² = P(m + η) of each; intercepts 0.52 fringe apart (3.3 rad) |
| Fig. 6 | `phase_vs_omega.pdf` | column | 2.3 | setting means of Δφ (ring minima) against ω (error bars: the calibration line's prediction uncertainty in ω, 1.35 %, and the SEM in Δφ); solid: equal-weight regression; dashed with band: a_th ω ± 1σ (N, R) |
| Fig. 7 | `methods.pdf` | column | 2.4 | (a) Δφ by the cosine fit and by template matching minus the ring-minima value, each relative to its own rest mean; (b) each method's slope by the same regression, error bars the regression; band: prediction ± 1σ |
| Fig. 8 | `drift.pdf` | column | 2.6 | (a) ω of each photo against time: settings in increasing order; (b) each photo's Δφ minus the line of Fig. 6, with a straight line per setting against time: the phase wanders within and between settings |
| (optional) | `random_pixels.pdf` | column | 2.4 | (a) slopes of 1000 draws of ten fixed pixels (one alternating least-squares fit per draw), outline with each photo's own center and P, lines at the minima slope and the prediction, 2 failed draws off the plot; (b) spread between draws and one draw's fit error against pixels per draw |
| Table 1 (optional) | | column | 2.1 | the nine timed settings: S, n, ⟨T⟩ with u(T), ω |
| Table 2 | | column | 2.3 | per setting: photos, time of the first, ω (line, with its prediction uncertainty), Δφ (SEM; the rest mean's 0.02 rad is common to all rows), P (median of the setting's photos) |
| Table 3 (optional) | | column | 2.4 | slope, s, R² of each method; method uncertainty = half the range |
| Table 4 (optional) | | column | 2.5 | relative budget of c, λ not included |
| Table 5 | | figure-width | 2.6 | every source, kind, size, whether included; how the slope's three parts combine; N and R enter the prediction and c; the bend of the calibration and the sign not covered |

`setup.pdf` replaces the manual's Fig. 1 (an image of a printed figure, with the small loop); if you prefer the manual's own figure, cite it as "taken from [manual]".

## Section plans

### 1.1 Theory

**1.1.1 The Sagnac phase.** Purpose: the formula tested, derived (prep Q4).
- Lab frame, loop of radius R turning at ω, the splitter moves at v = ωR, light at c both ways: t± = l/(c ∓ v), l = 2πR, so Δt = 2lv/(c² − v²) ≈ 4Aω/c² with A = πR² (manual Eq. 1–2; Pascoli Eq. 1; Tartaglia Eq. 1). N turns: Δt = 4NAω/c².
- A delay is a phase 2π(c/λ)Δt: **Δφ = 8πNAω/(λc)** (manual Eq. 3–5; Schwartz Eq. 1). Label `eq:sagnac`. So the slope is a_th = 8π²NR²/(λc).
- Independent of where the axis is (Schwartz Eq. 4, prep Q8) and of the fiber's refractive index (Schwartz Sec. 2, prep Q10): one sentence each, cited. These two explain why only N, R, λ enter.
- Optional: first order in v/c; at the rim v ≈ 1 m/s, v/c ≈ 3 × 10⁻⁹, and relativity changes Δt only at order (v/c)² (Pascoli Sec. 3; Tartaglia). The effect is measurable because interference resolves a fraction of an optical period (2.13 fs).
- The sign of Δφ follows the sense of rotation (needed for 2.6: not checked).

**1.1.2 Interference and the ring pattern.** Purpose: why the fringes are rings, evenly spaced in r², and what the analysis reads.
- Two coherent beams: I = I₁ + I₂ + 2√(I₁I₂) cos φ (manual Eq. 6); dark where φ = π + 2πm. Plane waves at an angle give straight fringes with period λ/(2 sin θ) (manual Eq. 6 and Fig. 2; prep Q1–Q2); ours are rings (Schwartz Sec. 3: his pattern is not linear either).
- A Gaussian beam carries the transverse phase kρ²/(2R(z)) (QE3, PDF p. 10; prep Q3). Two beams on a common axis with wavefront radii R₁ ≠ R₂ at the camera differ by **φ(r) = 2πr²/P + φ_c, P = 2λ/|1/R₁ − 1/R₂|**. Label `eq:rings`. Equal spacing in r² (like Newton's rings).
- Dark rings at **R_m² = P(m + η), φ_c = π − 2πη**. Label `eq:dark_rings`. Rotation adds the same phase everywhere, so it changes only φ_c: the rings slide in or out. A tilt moves the center; a lens–fiber distance changes P (prep Q7); neither moves φ_c (to first order).
- This is the manual's Eq. 7 (shift of each minimum over the local fringe period) and Schwartz's Eq. 5 (Δφ = 2πΔx/δ), done in r², where every ring has the same period P. η is a ratio of two lengths squared, so no pixel-to-length calibration is needed.
- Prep Q9 (wavelength): P ∝ λ for given wavefronts, and the Sagnac phase ∝ 1/λ. Optional.

**1.1.3 Fiber modes (short; needed for the discussion in 2.6).**
- V = (2π/λ) a NA (manual Eq. 9); one mode for V < 2.405. SMF-28e+, a = 4.1 µm, NA = 0.14: V = 5.6 at 640 nm (manual p. 10), so the fiber guides several modes (prep Q12: no, not single-mode at 640 nm).
- Coupling into the fundamental mode: at most 47 % with the setup's 8 mm lens, 93 % with a matched 15.5 mm lens (manual p. 10). Optional detail.

### 1.2 Setup and procedure

Purpose: what was done, so that the results can be followed. Fig. 1.
- Apparatus (manual Fig. 1, p. 12): 640 nm laser → P₁ → beam splitter → L₁ and L₂ into the two ends of the fiber loop → back through the splitter → P₂ → L₃ → camera. Everything on the table turns with it; the loop is centered on the rotation axis. The polarizers suppress light reflected from the lenses, holders and splitter that did not go around the loop (prep Q6; manual preparation question 3).
- Given (manual, Fig. 1 caption): N = 13 ± 1 turns, R = 57.0 ± 0.5 cm, λ = 640 nm without an uncertainty. Fiber: Corning SMF-28e+ (the manual's text says SMF-28e, its datasheet and the preparation sheet say SMF-28e+; say it once). A = πR² = 1.02 m² per turn.
- Rotation: the motor controller's setting S (a number, no unit, which sets the motor voltage; task 1 calls it the voltage). Calibration: at nine settings, S = 5.78–90, 5 consecutive turns (lowest five) or 10 (highest four) timed with a stopwatch, one press per turn.
- Photographs: 161, 1920 × 1080, at rest and at nine settings S = 12.59–90, 10–22 per setting, in increasing order of S, 16:21:40–16:33:46 on 11 October 2026, rest photos only at the start. Only the red channel carries light, maximum 239 of 255 (nothing saturated).
- Task 4 (small off-center loop, R = 15.5 ± 0.1 cm, 46 ± 1 turns) not done, on the instructor's instruction; Schwartz (2017) did exactly this test.
- Analysis, one paragraph: per photo, the ring finder (smoothing σ = 4 px against speckle; darkest point of each dark ring along 170 rays every 1.5° from −75° to 180°; concentric circles fitted to them; three rounds) → radii R_m → line R_m² against m → φ_c (`eq:dark_rings`); Δφ = φ_c − mean φ_c at rest; mean per setting; ω from the calibration line; regression Δφ = aω + b. Detail in 2.3; `docs/fringes.md` explains every step.
- AI note (your decision, see Acknowledgments): the analysis code was written with an AI assistant.

### 2.1 Rotation calibration (task 1)

Purpose: the conversion from S to ω that the photos need; task 1 asks for the linear ratio. Fig. 2, Fig. 2b (optional), Table 1 (optional).
- T = mean of n consecutive turn times, ω = 2π/T. Two parts of u(T): the lap scatter SD/√n, and the reaction time of the presses. The laps are consecutive, so the mean is (t_n − t_0)/n and only the first and last press enter it: with σ_r = 0.2 s per press that part is √2σ_r/n. In quadrature: 0.14–0.9 % of ω, the reaction time dominating at the fast settings. Example S = 90: 10 turns, ⟨T⟩ = 3.545 s, SD = 0.048 s, SD/√n = 0.015 s, √2σ_r/n = 0.028 s, u(T) = 0.032 s, ω = 1.772 ± 0.016 rad/s. Say that 0.2 s is a conventional value and that the lap scatter implies a smaller press jitter (safe side).
- The line ω = kS + ω₀ by weighted least squares (weights 1/u(ω)²): k = 0.0205 ± 0.0001 rad s⁻¹ per unit S (unrounded 0.02054 ± 0.00010, 0.50 %), ω₀ = −0.0226 ± 0.0007 rad/s. χ²/dof = 28.7/7 (p = 2 × 10⁻⁴): the points scatter more than their uncertainties allow, so the coefficient covariance is multiplied by χ²/dof (4.1), the usual remedy; the uncertainties above include that.
- The deviations from the line are small at low S and grow negative at high S (−0.9, −1.8, −3.0 % at the three fastest settings): the table turns slower than linear at high settings, 3.0 % below the line at S = 90 (line 1.826, timed 1.772 rad/s). The intercept has no physical meaning (the line was tested for S = 5.78–90 only).
- Used as task 1 asks. A photo setting was not timed, so the line predicts its ω; the prediction carries the line's own uncertainty (0.3–0.5 %) and the 1.3 % rms scatter of the timed settings about the line, 1.35 % in all (Table 2, Fig. 2b). Both are common to every setting, so the calibration part of the Sagnac slope's uncertainty is √(0.50² + 1.27²) % = 1.4 % (2.6). The 3.0 % at the top setting goes beyond that rms: state it in 2.6 as a limitation.

### 2.2 Interference patterns (task 2)

Purpose: the pattern at each speed (task 2) and what rotation does to it. Fig. 3, Fig. 4.
- Rings; as ω rises they contract: the first dark ring moves from 258 px at rest to 101 px at the top speed, and the center goes from bright to nearly dark, about half a fringe (Fig. 3).
- Evenly spaced in r², not in r (Fig. 4): in the fastest photo the dark rings are 241 px and then 132 px apart in r, but 1.08 × 10⁵ px² both times in r² (`eq:rings`). Ring spacing P: median 1.074 × 10⁵ px² over all photos; per-setting medians (Table 2) 1.058–1.094 × 10⁵ px² (±1.7 %).
- Rotation changes only φ_c (1.1.2): the pattern slides, P stays; the center wanders by about 3 px (SD) over the run.
- The camera never saturates (maximum 239 of 255); the minima positions do not depend on the camera's intensity response.

### 2.3 Phase shift and comparison with theory (task 3)

Purpose: deduce the phase shift and compare with the prediction. Fig. 5, Fig. 6, Table 2.

*Phase of one photo (Fig. 5).*
- Minima along the rays, concentric circles fitted: the single minima lie 2.8–8.2 px rms (median 4.7) from their circles (speckle, small against the 130–240 px between rings); each ring has 54–170 minima (median 136), so its radius is known to 0.3–2.3 px (median 0.8 px), the center to about ±0.6–0.7 px.
- Line R_m² = P(m + η) through the two or three rings (`eq:dark_rings`): at most one degree of freedom, so no fit statistic; its check is the circle fit. P is fitted in every photo, so a change of the ring spacing does not move φ_c.
- η stays between 0 and 1 for the whole run (0.62 at rest, 0.10 at the top speed), so the ring numbering never shifts and there is no 2π ambiguity. At rest two rings are in the fringe region; from S = 44.07 on every photo has three.
- Δφ = φ_c − mean of the 10 rest photos.

*Setting means (Table 2).* Mean of the photos, uncertainty the SEM (0.004–0.028 rad). The rest SEM (0.02 rad) is common to all rows.

*The line (Fig. 6).*
- Δφ = aω + b by ordinary least squares, ten points, equal weights, as Schwartz obtains his slopes (Eq. 7). The regression takes every point's uncertainty from the points: s = 0.044 rad, u(a) = s/√Σ(ω − ω̄)². Name each piece in words you can defend.
- **a = 1.83 ± 0.02 rad/(rad s⁻¹)** (regression; unrounded 1.8261 ± 0.0236), R² = 0.9987.
- Linear: an added ω² term is 0.5 standard errors from zero (p = 0.67); intercept b = −0.05 ± 0.03 rad (t test, p = 0.12), consistent with zero; forcing the line through the origin gives 1.79 ± 0.01 (check only).
- Why not the photos' error bars (one paragraph): weighted by the SEM, χ²/dof = 197/8 (p = 3 × 10⁻³⁸, slope 1.809): the setting means lie five times farther from the line than those error bars allow, because the photos of a setting are not independent (the phase wanders slowly, 2.6, Fig. 8). The regression's s contains every effect that moves a setting mean (photo noise, table speed, wander). All 161 photos as independent points would give a misleadingly small uncertainty.

*Comparison.*
- Prediction a_th = 1.7 ± 0.1 rad/(rad s⁻¹) (unrounded 1.738 ± 0.137, from N: 7.7 % and R: 2 × 0.88 % = 1.8 %).
- With the regression error only: a/a_th = 1.05 ± 0.08, 0.6σ (z = 0.63), p = 0.53. With the full uncertainty of 2.6: 1.05 ± 0.09, 0.6σ (z = 0.58), p = 0.56. Quote the second as the result; the first can go in parentheses.
- For scale: Schwartz's large ring 1.304 ± 0.038 against 1.344 ± 0.006 (agreement within 3 %).

*Largest shift.* Δφ = 3.28 ± 0.06 rad = 0.52 ± 0.01 fringe at ω = 1.83 rad/s (u = √2 s, top setting minus rest); Δt = Δφλ/(2πc) = 1.12 ± 0.02 fs. Predicted 3.2 ± 0.3 rad (3.17), 1.08 fs: half an optical period (2.13 fs).

### 2.4 Phase-reading methods (prep Q5)

Purpose: prep Q5 asks for two ways to extract Δφ and which is better when Δφ_max is large; and the slope should not depend much on how the phase is read. Fig. 7, Table 3 (optional).
- (i) Ring minima (main, 2.3): positions of minima, independent of the camera's intensity response (a monotonic response moves no minimum) and nearly independent of the beams' brightness envelopes (an envelope that changes across a fringe pulls a dip slightly toward its darker side; `docs/fringes.md`).
- (ii) Cosine fit of the sector-averaged intensity against s = r²: C e^(−s/s_C) + B e^(−s/s_B) cos(2πs/P + φ) + D (Gaussian envelopes, QE3 PDF p. 7; D a constant floor), over the three dark rings (s < 2.6 × 10⁵ px²). Uses every pixel but relies on the envelope model and a linear camera response. Its fits report an rms, not a χ², because neighboring pixels share speckle. Do not call D the camera's dark level: refitted on one photo per setting it comes out negative every time (`docs/fringes.md`), so it absorbs what the envelope model misses.
- (iii) Template matching (our addition): each r² profile is shifted and scaled onto the averaged rest profile, phase 2πD/P_T. No fringe model, but it reads the phase in the middle of its window, so a change of P adds to it.
- Slopes by the same regression: 1.83 ± 0.02 (R² 0.9987), 1.86 ± 0.03 (0.9982), 1.76 ± 0.02 (0.9993). Unrounded 1.8261, 1.8630, 1.7552.
- The cosine fit reads 0.10–0.19 rad higher than the minima at every moving setting (Fig. 7a), because at rest the bright flat center is not a cosine: fit rms 9.5 counts at rest, 1.4–2.3 above S = 44.07.
- Method uncertainty: half the range of the three slopes, (1.8630 − 1.7552)/2 = 0.054 (2.9 %), taken as one standard uncertainty. Say it is a convention and that the methods read the same photographs, so they are correlated checks, not independent measurements.
- Optional, answering the second half of prep Q5 (random-pixel figure): ten pixels at fixed image positions, I = a_j + b_j cos(ψ_j + δ_k), one alternating least-squares fit for all photos (the idea of Wang and Han's phase-shifting algorithm; cite if used). Over 1000 random draws: mean slope 1.70, spread between draws 0.20 (SD), 3.4 times one draw's own fit error (0.057). More pixels per draw reduce the spread (20: 0.09, 40: 0.07, 80: 0.06) but not below the fit error. Why: our largest shift is only half a fringe, so each pixel passes through less than one bright–dark cycle; a fixed pixel also reads a moving center or a change of P as phase (with each photo's own center and P the mean moves to 1.78). For a large Δφ_max, where fixed points pass through many cycles, it would do better; our run is the small-shift case. Not used for the result. These slopes are the presentation stage's scan (its weighted line and earlier calibration, not the regression; `docs/random_points_scan.md`).

### 2.5 Speed of light (extension, beyond the manual)

Purpose: the slope read through `eq:sagnac` gives c; it is the same comparison as a/a_th, not an independent test. Could be one paragraph at the end of 2.3. Table 4 (optional).
- **c = 8π²NR²/(λa) = (2.9 ± 0.2) × 10⁸ m/s** (unrounded 2.854 ± 0.246, 8.6 %); c/c₀ = 0.95 ± 0.08; −0.6σ (z = −0.59), p = 0.56; c₀ = 299 792 458 m/s (defined). Inferred through the formula, not timed.
- Budget: N 7.7 %, method 2.9 %, R 1.8 %, calibration 1.4 %, regression 1.3 %; λ not included (0.16 % per nm). Our own measurement is 3.5 %; the given geometry is the rest.
- Read the other way (c, R, λ known): N = 13.7 ± 0.5 turns, against the given 13 ± 1.
- Optional: with the laser frequency ν instead of λ, c = √(8π²NR²ν/a) with half the relative uncertainty.

### 2.6 Uncertainties and drift (task 5)

Purpose: every source, its kind and size, what is in the quoted uncertainty and what is not. Table 5, Fig. 8. Material: `uncertainty/REPORT.md`.
- Table 5 (sources). The slope's final uncertainty: regression 1.3 %, method 2.9 %, calibration 1.4 %, in quadrature 3.5 %: **a = 1.83 ± 0.06 rad/(rad s⁻¹)**. They come from different things (noise between settings, choice of phase reading, the calibration line and the table's scatter about it) and are added as independent.
- Not covered, say it plainly: the bend of the calibration beyond its rms (the line is 3.0 % above the timed ω at S = 90; the rms of the deviations, 1.3 %, is in the calibration part, the 3.0 % at the top is not; a curve through the timed points would lower a by a few percent); a steady drift in time; λ; the sign.
- Drift (Fig. 8): settings in increasing order, rest photos only at the start, so ω and time rise together. Joint regression Δφ = aω + Dt + b: a = 2.1 ± 0.2, D = −0.8 ± 0.5 mrad/s, p = 0.12 for D, correlation between a and D −0.99: the data cannot tell a steady drift from the slope. Within each setting the phase drifts by up to ±8 mrad/s, changing sign. A steady drift is not in the quoted uncertainty; the wander between settings is (through s).
- Plausible source of the wander, not tested: the fiber guides several modes at 640 nm (1.1.3); light in several modes is sensitive to stress and temperature, and the camera records the whole field with no single-mode filter at the common port, so the two directions need not stay exactly reciprocal.
- Not checked: the sign of the shift (the rotation direction was not recorded). λ: no uncertainty given. ω during the photos: from the calibration at the same setting.
- Next time: rest photos between the speeds and the rotation reversed (the Sagnac phase changes sign with the direction, a drift does not; Pascoli Sec. 2 describes Sagnac doing this); speeds in random order; count the fiber turns (N is the largest uncertainty); time the turns at the photo settings themselves, with more turns at the fast settings where the reaction time dominates; couple into one fiber mode with a longer-focal-length lens; record the direction.

## Discussion items (where they go)

- [ ] What limits the comparison: the given N (7.7 %) against our 3.5 % (2.3 or 2.5; conclusions).
- [ ] The linear calibration, used as asked, and the bend of the table's speed at high settings (2.1, 2.6).
- [ ] Why a regression: the photos of a setting are not independent; the wander (2.3, 2.6).
- [ ] The cosine fit's offset at rest; template matching and P; fixed pixels for half a fringe (2.4).
- [ ] The steady drift that the order of the measurements cannot separate; how to test it (2.6).
- [ ] The multimode fiber as a plausible source of the wander, untested (2.6, with 1.1.3).
- [ ] Linearity: no ω² term (p = 0.67), intercept consistent with zero (2.3).
- [ ] Comparison with Schwartz's agreement (3 % on his large ring) (2.3, optional).

## Abstract plan

One paragraph, 150–200 words, no figure or section references, only numbers in the report. (1) What: the Sagnac effect in a fiber loop of 13 turns rotating with the whole interferometer. (2) How: the table speed calibrated against the controller setting by timing turns; 161 photographs at rest and at nine speeds; the fringes are rings evenly spaced in r², and the phase at their center follows from a straight-line fit to the squared ring radii. (3) Results: Δφ linear in ω, slope 1.83 ± 0.06 rad/(rad s⁻¹) against the predicted 1.7 ± 0.1 (0.6σ); largest shift 3.28 ± 0.06 rad (half a fringe, 1.1 fs); c = (2.9 ± 0.2) × 10⁸ m/s. (4) Limits: the given turn count; a slow drift of the phase that the order of the measurements cannot separate. (5) Takeaway, one clause. Recheck sentence (1) once Section 1 exists.

## Conclusions plan

Two short paragraphs, about 200 words, only numbers already in Section 2. P1, what was measured: calibration (linear, as asked; the table bends by 3 % at the top); rings evenly spaced in r², phase at the center from the ring radii; Δφ linear in ω (R² = 0.9987, no ω² term); a = 1.83 ± 0.06 against 1.7 ± 0.1 (0.6σ, p = 0.56); largest shift 3.28 ± 0.06 rad; c = (2.9 ± 0.2) × 10⁸ m/s, or N = 13.7 ± 0.5 read the other way. P2, limits and next steps: N dominates; the method choice 3 %; the settings wander by about 0.04 rad and a steady drift is untested because of the increasing order; rest photos between speeds, reversed rotation, random order, counted turns, single-mode coupling.

## Acknowledgments and AI disclosure

Facts to cover (the Chua acknowledgments were checked against APS, AIP, IOP and Optica policies on 3 October: developers and versions, what the AI covered including figures, how its output was directed and verified, and that the authors take responsibility):
- [ ] Which assistants were used, with versions (fill in from your own records; this session is Claude, Anthropic). If `sagnac_update/` came from another assistant, name it too.
- [ ] What they covered: all analysis code (`calibration.py`, `fringes.py`, `methods.py`, `sagnac.py`, `report_figures.py`, the earlier stages), every figure, the presentation, the results and uncertainty documents, and the derivations and analysis choices suggested along the way (the r² ring model and the ring finder, the regression treatment, the method uncertainty, the drift checks).
- [ ] How you directed and verified it: the analysis choices you made (linear calibration as task 1 asks, the regression as in Schwartz, which methods count, the review of 22 October that simplified the uncertainty treatment, the feedback after the presentation) and the checks you ran (against the manual and the papers, the tests, the presentation review). Only what you actually did.
- [ ] In-text mentions where a method is described (your decision; the Chua report put one in 2.3 and one in 2.5): e.g. 1.2 or 2.3 for the ring finder and the phase extraction.
- [ ] Your partner's and the instructor's contributions as you see fit. The friend's `sagnac_2d.py` (2D Fourier method) is only a check in `methods.py`; credit it only if the report mentions that method.

## References (checked against the PDFs in `Documents/`)

1. Sagnac experiment, laboratory manual, Physics Lab 5, Technion (with the Corning SMF-28e+ product sheet). `Documents/Sagnac_manual.pdf`
2. Sagnac effect: preparing for the first meeting (preparation sheet), Physics Lab 5, Technion. `Documents/Entry - Sagnac effect experiment.pdf`
3. E. Schwartz, "Sagnac effect in an off-center rotating ring frame of reference," Eur. J. Phys. **38**, 015301 (2017), doi:10.1088/0143-0807/38/1/015301.
4. G. Pascoli, "The Sagnac effect and its interpretation by Paul Langevin," C. R. Physique **18**, 563–569 (2017).
5. A. Tartaglia and M. L. Ruggiero, "The Sagnac effect and pure geometry," Am. J. Phys. **83**, 427–432 (2015), doi:10.1119/1.4904319.
6. J. T. Costello, "The Gaussian beam," lecture slides (QE3). `Documents/QE3.pdf` (the phase of a Gaussian beam: PDF p. 10)
7. Optional: G. Sagnac, C. R. Acad. Sci. **157**, 708 and 1410 (1913) (as cited by Schwartz); K. U. Schreiber and J.-P. R. Wells, Rev. Sci. Instrum. **84**, 041101 (2013) (ring lasers, the manual's ref. 2); Z. Wang and B. Han, Opt. Lett. **29**, 1671 (2004), doi:10.1364/OL.29.001671 (only with the fixed-pixel method); JCGM 100:2008 (GUM) (only with an uncertainty appendix).

## Known traps

- ω from the earlier calibrations (1.78 rad/s at the top in the old `figures/patterns.png`, 1.81 in the presentation): the report uses the current line (1.83).
- Earlier documents quoted a = 1.84 ± 0.07 (presentation and first review), 1.84 ± 0.06 (`sagnac_update/RESULTS.md`) or 1.82 ± 0.06 (22 October, before the calibration file was completed), and a largest shift of 3.29 ± 0.06: the report's values are 1.83 ± 0.06 and 3.28 ± 0.06.
- 1.83 and 1.7 look 0.13 apart; the unrounded values are 0.088 apart. Compute σ and p from the unrounded numbers.
- c/c₀ and a/a_th are the same comparison; do not present c as an independent confirmation.
- R² = 0.9987 does not by itself show the model is right; the ω² term and the intercept do (there is no residual panel to point at).
- The SEM error bars of Table 2 are not the uncertainties the line uses (s is); say so once.
- "Δφ_max large" is not our case: half a fringe.
- The lower left of every photo has no fringes; the sector −75° to 180° is a choice, say it.
- Task 1 says "voltage on the motor": we have the controller setting S, not a measured voltage.
- Table 2's P column and the 1.058–1.094 range are the medians of each setting's photos.
- The random-pixel slopes are the presentation stage's (weighted line, earlier calibration), not the regression; say so if the figure is used.
- The 0.2 s reaction time is an assumed value (`REACTION_TIME_S` in `calibration.py`); if the course prescribes another, change it there and rerun.

## Coming next

- Drafts: review them here against this plan (numbers, figure and table citations, rules), newest first.
