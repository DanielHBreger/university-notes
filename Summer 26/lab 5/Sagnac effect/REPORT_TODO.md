# Lab 5 report (Sagnac): plan and open items

Started 22 October 2026, before the first draft; revised 23 October after the presentation feedback (no residual panels; physics and findings, not analysis mechanics) and the completed calibration file. As in `chaotic circuit/REPORT_TODO.md`: tick an item (`[x]`) or delete it when a draft covers it, and put each draft review at the top of the file, newest first. Numbers here match `RESULTS.md`, `sagnac.txt` and `calibration.txt`; `RESULTS.md` holds the detail behind every number, `docs/` the code.

What goes with this file:

| File | What it is |
|---|---|
| `report_figures/` (`report_figures.py`) | the report's figures, PDF + 600 dpi PNG, in the report style, data and fits only |
| `REPORT_figures_tables.tex` | the figure and table blocks (numbers and `\includegraphics` only, captions left to you) |
| `PLAN_CHECKLIST.md` | the manual's tasks and the preparation questions against what was done |
| `uncertainty/REPORT.md` | every uncertainty component and its derivation; material for task 5 and an appendix |
| `docs/` | the code explained script by script, for questions, not for the report |

## Rules for the report

- The reader wants the physics and how it was found. Say what was measured, what was fitted, what came out and what it means. Leave out how the fit was computed, which routine or convention produced an uncertainty, pixel counts, smoothing widths, ray spacings and other mechanics that change nothing. If a detail would not change a number or a conclusion, it does not go in.
- Every number with its uncertainty and what kind it is (statistical, systematic, given). Uncertainty rounded to its first significant figure and the value to the same digit (1.7 ± 0.1, not 1.74 ± 0.14); every comparison (σ, p) computed from the unrounded values, given below where needed.
- Figures from `report_figures/` (column width unless figure*): data and fits only, no residual panels, no fit results inside the figures; the numbers go in captions and text.
- **The prose is yours.** This file gives content (facts, numbers, order), not sentences; the figure and table blocks hold no text. (The Chua 2.4/2.5 lesson: suggested text pasted almost verbatim had to be rewritten before submission.)
- Methods and terms from the lab documents: the manual, Schwartz 2017, the QE3 Gaussian-beam slides, the preparation sheet. Template matching and the fixed-pixel method are our own additions: say so where they appear.
- Plain words for fit quality: "the points lie farther from the line than their error bars allow", "the rms distance of the points from the line", "an added ω² term". Not "scatter" or "curvature" for a statistic.
- American spelling (center, fiber), as in `RESULTS.md`. Siunitx as in the Chua report (`\qty`, `\num` with `\pm`, compact output).
- One unit for the slope everywhere: rad/(rad s⁻¹) (= s). `REPORT_figures_tables.tex` defines `\slopeunit` once.
- ω always from the current calibration line (top speed 1.83 rad/s). `previous_reports/`, the old `figures/` folder and `presentation_slides/source/numbers.json` are earlier stages: take no numbers from them.

## BLOCKING before submission

- [ ] All prose written in your own words (abstract, every section, captions).
- [ ] Every number checked against `RESULTS.md` (or `sagnac.txt`); the slope is 1.83 ± 0.06 with the full uncertainty and 1.83 ± 0.02 from the fit alone: say which one each comparison uses.
- [ ] Acknowledgments with the AI disclosure (section below), and the in-text mentions you decide on.
- [ ] References and citations: the manual, the preparation sheet, Schwartz, QE3 at least (list below).
- [ ] Task 4 (off-center small loop) stated as not done, on the instructor's instruction.
- [ ] Every figure and table cited in the text; each figure* defined before the text that cites it.

## Structure (proposed; the Chua report's layout)

| Section | Content | Figures, tables |
|---|---|---|
| Abstract | one paragraph, 150–200 words, no references to figures or sections | |
| 1.1 Theory | 1.1.1 the Sagnac phase; 1.1.2 two-beam interference and the ring pattern; 1.1.3 fiber modes (short) | Fig. 4 can be cited from 1.1.2 |
| 1.2 Setup and procedure | apparatus, given parameters, timing, photographs, the analysis in one paragraph | Fig. 1 |
| 2.1 Rotation calibration (task 1) | the timing and its uncertainty, the line, the conversion | Fig. 2, Table 1 |
| 2.2 Interference patterns (task 2) | the rings and what rotation does to them; equal spacing in r² | Fig. 3 (figure*), Fig. 4 |
| 2.3 Phase shift and comparison with theory (task 3) | phase of one photo; setting means; the line; comparison; largest shift | Fig. 5, Fig. 6, Table 2 |
| 2.4 Phase-reading methods (prep Q5) | three ways of reading the phase; the method uncertainty; fixed pixels (optional) | Fig. 7, Table 3 (optional) |
| 2.5 Speed of light (extension) | c from the slope; budget; N read the other way | Table 4 (optional) |
| 2.6 Uncertainties and drift (task 5) | budget; the final slope; drift; what is not covered; next time | Table 5 (figure-width), Fig. 8 |
| 3 Conclusions | two short paragraphs | |
| Acknowledgments, References | | |
| Appendix A (optional) | uncertainty propagation, from `uncertainty/REPORT.md` | |

Results and Discussion merged, as in the Chua report. The conversion figure (`conversion.pdf`) is optional: it draws the ω column of Table 2 on the line of Fig. 2; use the table unless you want the picture.

## Figures and tables

Numbering as in the structure above, assuming every optional table is used (LaTeX renumbers by order of appearance, so refer to them by label). Each caption must say the items listed; write them in your words. The same lists are in `REPORT_figures_tables.tex`.

| | File | Width | Section | The caption must say |
|---|---|---|---|---|
| Fig. 1 | `setup.pdf` | column | 1.2 | top view, after Fig. 1 of the manual (cite); P₁, P₂ polarizers, L₁, L₂ coupling lenses, L₃ camera lens; beams 1 and 2 go around the same loop (N turns, radius R) in opposite directions and meet at the splitter; everything turns with the table about the loop's center; the ω arrow is schematic (direction not recorded); the small loop of task 4 is not drawn |
| Fig. 2 | `calibration.pdf` | column | 2.1 | ω of the table at each controller setting S, from the timed turns (error bars: the timing uncertainty, mostly smaller than the markers); line: the linear calibration used for the photos; band: its ±1σ uncertainty, which is the uncertainty of ω at any photo setting |
| (optional) | `conversion.pdf` | column | 2.1 | the nine photo settings read off the line of Fig. 2: dotted guides from each setting to the line and across to ω; nothing was measured at these settings; the values are those of Table 2 |
| Fig. 3 | `patterns.pdf` | figure* | 2.2 | one photo per setting in the order taken; the fitted dark rings (dashed) and ring center (+); ω from the calibration; the rings move inward as ω rises; no fringes in the lower left |
| Fig. 4 | `rings.pdf` | column | 2.2 | intensity averaged around the ring center for the rest and top-speed photos, (a) against r, (b) against r²; dashed: the dark-ring radii; arrows: the top-speed spacings, unequal in r, equal in r² |
| Fig. 5 | `extraction.pdf` | column | 2.3 | (a) top-speed photo with the located dark-ring minima (one color per ring) and the circles fitted to them; (b) R_m²/P against the ring number m at rest and at top speed, with the fitted lines; the intercepts differ by 0.52 fringe (3.3 rad) |
| Fig. 6 | `phase_vs_omega.pdf` | column | 2.3 | phase shift of each setting against ω (error bars: calibration in ω, standard error of the photos in Δφ); solid: the fitted line; dashed with band: the Sagnac prediction a_th ω ± 1σ from N and R |
| Fig. 7 | `methods.pdf` | column | 2.4 | (a) the phase shift by the cosine fit and by template matching minus the ring-minima value; (b) the slope from each method with its fit uncertainty; band: the prediction ± 1σ |
| Fig. 8 | `drift.pdf` | column | 2.6 | (a) ω against time: the settings were taken in increasing order; (b) each photo's phase shift minus the line of Fig. 6 against time, with a trend line per setting: the phase wanders within and between settings |
| (optional) | `random_pixels.pdf` | column | 2.4 | the slope from ten fixed pixels over 1000 random draws, and the spread between draws against the number of pixels per draw |
| Table 1 | | column | 2.1 | the nine timed settings: S, n, ⟨T⟩ with its uncertainty, ω |
| Table 2 | | column | 2.3 | per setting: photos, time of the first, ω (from the calibration, with its uncertainty), Δφ (standard error of the photos; the rest value's 0.02 rad is common to all rows), ring spacing P |
| Table 3 (optional) | | column | 2.4 | the slope from each method; method uncertainty = half the range |
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
- Photographs: 161, at rest and at nine settings S = 12.59–90, 10–22 per setting, in increasing order of S, 16:21–16:34 on 11 October 2026, rest photos only at the start. Only the red channel carries light; nothing saturated.
- Task 4 (small off-center loop, R = 15.5 ± 0.1 cm, 46 ± 1 turns) not done, on the instructor's instruction; Schwartz (2017) did exactly this test.
- Analysis, one paragraph: in each photo the dark rings are located and concentric circles fitted to them; the squared radii against the ring number give the phase at the ring center (`eq:dark_rings`); Δφ is that phase minus its value at rest, averaged over the photos of a setting; ω comes from the calibration; a straight line through the ten settings gives the slope. Detail in 2.3.
- AI note (your decision, see Acknowledgments): the analysis code was written with an AI assistant.

### 2.1 Rotation calibration (task 1)

Purpose: the conversion from S to ω that the photos need; task 1 asks for the linear ratio. Fig. 2, Table 1.
- At each setting the mean turn time T gives ω = 2π/T. Its uncertainty has two parts: the spread of the single turns, and the reaction time of the stopwatch presses (0.2 s per press). The turns are consecutive, so the mean is set by the first and last press only, and the reaction-time part is √2 × 0.2 s divided by the number of turns. Together 0.14–0.9 % of ω, the reaction time dominating at the fast settings. One worked example, S = 90: 10 turns, ⟨T⟩ = 3.545 s, u(T) = 0.032 s, ω = 1.772 ± 0.016 rad/s. Table 1.
- The line ω = kS + ω₀: k = 0.02054 ± 0.00005 rad s⁻¹ per unit S, ω₀ = −0.0226 ± 0.0003 rad/s. The points lie farther from the line than their error bars allow (χ²/dof = 4.1): the table turns slower than linear at the fast settings, 3.0 % below the line at S = 90 (line 1.826, timed 1.772 rad/s). The intercept has no physical meaning (the line was tested for S = 5.78–90 only). Fig. 2.
- Used as task 1 asks. A photo setting was not timed, so ω is read from the line; its uncertainty combines how well the line is known with how far the timed settings lie from it, 1.3 % of ω, the same for every setting (Table 2). The bend at the fast settings goes beyond that: state it in 2.6 as a limitation.

### 2.2 Interference patterns (task 2)

Purpose: the pattern at each speed (task 2) and what rotation does to it. Fig. 3, Fig. 4.
- Rings; as ω rises they contract: the first dark ring moves from 258 px at rest to 101 px at the top speed, and the center goes from bright to nearly dark, about half a fringe (Fig. 3).
- Evenly spaced in r², not in r (Fig. 4): in the fastest photo the dark rings are 241 px and then 132 px apart in r, but 1.08 × 10⁵ px² both times in r² (`eq:rings`). Ring spacing P: median 1.07 × 10⁵ px²; per-setting values within ±1.7 % (Table 2).
- Rotation changes only φ_c (1.1.2): the pattern slides, P stays.
- The camera never saturates; the positions of minima do not depend on the camera's intensity response.

### 2.3 Phase shift and comparison with theory (task 3)

Purpose: deduce the phase shift and compare with the prediction. Fig. 5, Fig. 6, Table 2.

*Phase of one photo (Fig. 5).*
- The dark rings are located and concentric circles fitted; a ring radius is known to about 1 px, the center to about 1 px.
- The line R_m² = P(m + η) through the two or three rings (`eq:dark_rings`) gives φ_c. P is fitted in every photo, so a change of the ring spacing does not move φ_c.
- η stays between 0 and 1 for the whole run (0.62 at rest, 0.10 at the top speed), so the ring numbering never shifts and there is no 2π ambiguity.
- Δφ = φ_c − its mean over the 10 rest photos.

*Setting means (Table 2).* Mean of the photos, uncertainty the standard error of the mean (0.004–0.028 rad). The rest value's 0.02 rad is common to all rows.

*The line (Fig. 6).*
- A straight line Δφ = aω + b through the ten setting means, as Schwartz obtains his slopes (Eq. 7). **a = 1.83 ± 0.02 rad/(rad s⁻¹)** (unrounded 1.8261 ± 0.0236), R² = 0.9987.
- The uncertainty of a comes from how far the points lie from the line (0.044 rad rms), not from the photos' standard errors: the setting means lie five times farther from the line than those allow, because the phase wanders slowly between settings (2.6, Fig. 8). One or two sentences.
- Linear: an added ω² term is 0.5 standard errors from zero (p = 0.67); the intercept b = −0.05 ± 0.03 rad is consistent with zero (p = 0.12).

*Comparison.*
- Prediction a_th = 1.7 ± 0.1 rad/(rad s⁻¹) (unrounded 1.738 ± 0.137, from N: 7.7 % and R: 2 × 0.88 % = 1.8 %).
- With the fit uncertainty only: a/a_th = 1.05 ± 0.08, 0.6σ (z = 0.63), p = 0.53. With the full uncertainty of 2.6: 1.05 ± 0.09, 0.6σ (z = 0.58), p = 0.56. Quote the second as the result; the first can go in parentheses.
- For scale: Schwartz's large ring 1.304 ± 0.038 against 1.344 ± 0.006 (agreement within 3 %).

*Largest shift.* Δφ = 3.28 ± 0.06 rad = 0.52 ± 0.01 fringe at ω = 1.83 rad/s; Δt = Δφλ/(2πc) = 1.12 ± 0.02 fs. Predicted 3.2 ± 0.3 rad (3.17), 1.08 fs: half an optical period (2.13 fs).

### 2.4 Phase-reading methods (prep Q5)

Purpose: prep Q5 asks for two ways to extract Δφ and which is better when Δφ_max is large; and the slope should not depend much on how the phase is read. Fig. 7, Table 3 (optional).
- (i) Ring minima (main, 2.3): the positions of the dark rings. Independent of the camera's intensity response.
- (ii) Cosine fit: a cosine in r² with the envelopes of two Gaussian beams, fitted to the intensity averaged around the ring center; the same φ_c from every pixel, at the price of trusting the envelope model and a linear camera.
- (iii) Template matching (our addition): each photo's intensity profile is shifted onto the rest profile; the shift is the phase change. No fringe model, but it reads an average phase over the profile, so a change of P enters it.
- Slopes: 1.83 ± 0.02, 1.86 ± 0.03, 1.76 ± 0.02 (unrounded 1.8261, 1.8630, 1.7552). The cosine fit reads 0.10–0.19 rad higher than the minima at every moving setting, because at rest the bright flat center is not a cosine (Fig. 7a).
- Method uncertainty: half the range of the three slopes, 0.054 (2.9 %), taken as one standard uncertainty. Say it is a convention and that the methods read the same photographs, so they are checks of each other, not independent measurements.
- Optional, the second half of prep Q5: ten pixels at fixed image positions, the phase of every photo from their brightness (Wang and Han's phase-shifting idea; cite if used). Over 1000 random draws of ten pixels the slope has mean 1.70 and spreads by 0.20 between draws, 3.4 times the fit error of one draw; more pixels narrow the spread but not to the fit error. Why: our largest shift is only half a fringe, so each pixel passes through less than one bright–dark cycle, and a fixed pixel also reads a moving center or a change of P as phase. For a large Δφ_max this method would do well; our run is the small-shift case. Not used for the result.

### 2.5 Speed of light (extension, beyond the manual)

Purpose: the slope read through `eq:sagnac` gives c; it is the same comparison as a/a_th, not an independent test. Could be one paragraph at the end of 2.3. Table 4 (optional).
- **c = 8π²NR²/(λa) = (2.9 ± 0.2) × 10⁸ m/s** (unrounded 2.854 ± 0.246, 8.6 %); c/c₀ = 0.95 ± 0.08; −0.6σ (z = −0.59), p = 0.56; c₀ = 299 792 458 m/s (defined). Inferred through the formula, not timed.
- Budget: N 7.7 %, method 2.9 %, R 1.8 %, calibration 1.3 %, fit 1.3 %; λ not included (0.16 % per nm). Our own measurement is 3.5 %; the given geometry is the rest.
- Read the other way (c, R, λ known): N = 13.7 ± 0.5 turns, against the given 13 ± 1.
- Optional: with the laser frequency ν instead of λ, c = √(8π²NR²ν/a) with half the relative uncertainty.

### 2.6 Uncertainties and drift (task 5)

Purpose: every source, its kind and size, what is in the quoted uncertainty and what is not. Table 5, Fig. 8. Material: `uncertainty/REPORT.md`.
- Table 5 (sources). The slope's uncertainty has three parts: the fit (1.3 %, how far the settings lie from the line), the choice of phase-reading method (2.9 %) and the calibration (1.3 %, how well the line gives ω). In quadrature 3.5 %: **a = 1.83 ± 0.06 rad/(rad s⁻¹)**. They come from different things and are added as independent.
- Not covered, say it plainly: the bend of the calibration at the fast settings beyond its average (the line is 3.0 % above the timed ω at S = 90; a curve through the timed points would lower a by a few percent); a steady drift in time; λ; the sign.
- Drift (Fig. 8): settings in increasing order with rest photos only at the start, so ω and time rise together. Fitting a drift term together with the slope gives a = 2.1 ± 0.2, D = −0.8 ± 0.5 mrad/s (p = 0.12), and the two cannot be told apart. Within each setting the phase drifts by up to ±8 mrad/s, changing sign. A steady drift is not in the quoted uncertainty; the wander between settings is.
- Plausible source of the wander, not tested: the fiber guides several modes at 640 nm (1.1.3); light in several modes is sensitive to stress and temperature, and the camera records the whole field with no single-mode filter at the common port, so the two directions need not stay exactly reciprocal.
- Not checked: the sign of the shift (the rotation direction was not recorded). λ: no uncertainty given. ω during the photos: from the calibration at the same setting.
- Next time: rest photos between the speeds and the rotation reversed (the Sagnac phase changes sign with the direction, a drift does not; Pascoli Sec. 2 describes Sagnac doing this); speeds in random order; count the fiber turns (N is the largest uncertainty); time the turns at the photo settings themselves, with more turns at the fast settings; couple into one fiber mode with a longer-focal-length lens; record the direction.

## Discussion items (where they go)

- [ ] What limits the comparison: the given N (7.7 %) against our 3.5 % (2.3 or 2.5; conclusions).
- [ ] The linear calibration, used as asked, and the bend of the table's speed at the fast settings (2.1, 2.6).
- [ ] Why the slope's uncertainty is larger than the photos' standard errors: the wander (2.3, 2.6).
- [ ] The cosine fit's offset at rest; template matching and P; fixed pixels for half a fringe (2.4).
- [ ] The steady drift that the order of the measurements cannot separate; how to test it (2.6).
- [ ] The multimode fiber as a plausible source of the wander, untested (2.6, with 1.1.3).
- [ ] Linearity: no ω² term, intercept consistent with zero (2.3).
- [ ] Comparison with Schwartz's agreement (3 % on his large ring) (2.3, optional).

## Abstract plan

One paragraph, 150–200 words, no figure or section references, only numbers in the report. (1) What: the Sagnac effect in a fiber loop of 13 turns rotating with the whole interferometer. (2) How: the table speed calibrated against the controller setting by timing turns; 161 photographs at rest and at nine speeds; the fringes are rings evenly spaced in r², and the phase at their center follows from the squared ring radii. (3) Results: Δφ linear in ω, slope 1.83 ± 0.06 rad/(rad s⁻¹) against the predicted 1.7 ± 0.1 (0.6σ); largest shift 3.28 ± 0.06 rad (half a fringe, 1.1 fs); c = (2.9 ± 0.2) × 10⁸ m/s. (4) Limits: the given turn count; a slow drift of the phase that the order of the measurements cannot separate. (5) Takeaway, one clause. Recheck sentence (1) once Section 1 exists.

## Conclusions plan

Two short paragraphs, about 200 words, only numbers already in Section 2. P1, what was measured: calibration (linear, as asked; the table bends by 3 % at the top); rings evenly spaced in r², phase at the center from the ring radii; Δφ linear in ω (no ω² term); a = 1.83 ± 0.06 against 1.7 ± 0.1 (0.6σ, p = 0.56); largest shift 3.28 ± 0.06 rad; c = (2.9 ± 0.2) × 10⁸ m/s, or N = 13.7 ± 0.5 read the other way. P2, limits and next steps: N dominates; the method choice 3 %; the settings wander by about 0.04 rad and a steady drift is untested because of the increasing order; rest photos between speeds, reversed rotation, random order, counted turns, single-mode coupling.

## Acknowledgments and AI disclosure

Facts to cover (the Chua acknowledgments were checked against APS, AIP, IOP and Optica policies on 3 October: developers and versions, what the AI covered including figures, how its output was directed and verified, and that the authors take responsibility):
- [ ] Which assistants were used, with versions (fill in from your own records; this session is Claude, Anthropic). If `sagnac_update/` came from another assistant, name it too.
- [ ] What they covered: all analysis code, every figure, the presentation, the results and uncertainty documents, and the analysis choices suggested along the way (the ring model in r², the method uncertainty, the drift checks).
- [ ] How you directed and verified it: the choices you made (linear calibration as task 1 asks, which methods count, the review and simplification of 22 October, the feedback after the presentation) and the checks you ran (against the manual and the papers, the tests, the presentation review). Only what you actually did.
- [ ] In-text mentions where a method is described (your decision; the Chua report put one in 2.3 and one in 2.5).
- [ ] Your partner's and the instructor's contributions as you see fit.

## References (checked against the PDFs in `Documents/`)

1. Sagnac experiment, laboratory manual, Physics Lab 5, Technion (with the Corning SMF-28e+ product sheet). `Documents/Sagnac_manual.pdf`
2. Sagnac effect: preparing for the first meeting (preparation sheet), Physics Lab 5, Technion. `Documents/Entry - Sagnac effect experiment.pdf`
3. E. Schwartz, "Sagnac effect in an off-center rotating ring frame of reference," Eur. J. Phys. **38**, 015301 (2017), doi:10.1088/0143-0807/38/1/015301.
4. G. Pascoli, "The Sagnac effect and its interpretation by Paul Langevin," C. R. Physique **18**, 563–569 (2017).
5. A. Tartaglia and M. L. Ruggiero, "The Sagnac effect and pure geometry," Am. J. Phys. **83**, 427–432 (2015), doi:10.1119/1.4904319.
6. J. T. Costello, "The Gaussian beam," lecture slides (QE3). `Documents/QE3.pdf` (the phase of a Gaussian beam: PDF p. 10)
7. Optional: G. Sagnac, C. R. Acad. Sci. **157**, 708 and 1410 (1913) (as cited by Schwartz); K. U. Schreiber and J.-P. R. Wells, Rev. Sci. Instrum. **84**, 041101 (2013) (ring lasers, the manual's ref. 2); Z. Wang and B. Han, Opt. Lett. **29**, 1671 (2004), doi:10.1364/OL.29.001671 (only with the fixed-pixel method); JCGM 100:2008 (GUM) (only with an uncertainty appendix).

## Known traps

- Earlier documents quoted a = 1.84 ± 0.07, 1.84 ± 0.06 or 1.82 ± 0.06 and a top speed of 1.78 or 1.81 rad/s: the report's values are 1.83 ± 0.06 and 1.83 rad/s.
- 1.83 and 1.7 look 0.13 apart; the unrounded values are 0.088 apart. Compute σ and p from the unrounded numbers.
- c/c₀ and a/a_th are the same comparison; do not present c as an independent confirmation.
- The standard errors in Table 2 are not the uncertainty the slope uses; say so once (2.3).
- "Δφ_max large" is not our case: half a fringe.
- Task 1 says "voltage on the motor": we have the controller setting S, not a measured voltage.
- The 0.2 s reaction time is an assumed value; if the course prescribes another, say so and rerun (`REACTION_TIME_S` in `calibration.py`).

## Coming next

- Drafts: review them here against this plan (numbers, figure and table citations, rules), newest first.
