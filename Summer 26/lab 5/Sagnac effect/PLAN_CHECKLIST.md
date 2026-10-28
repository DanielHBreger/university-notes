# Sagnac manual and preparation sheet against what was done

Status key: **done** (complete, and its output is in the folder), **done, caveat** (complete, with something the reader must know), **partly**, **not verifiable** (an instruction that leaves no trace in the data), **not done**. "Report" points to the section of the proposed structure in `REPORT_TODO.md`.

Sources: `Documents/Sagnac_manual.pdf` (minimum tasks and Fig. 1 on p. 12, fiber coupling on p. 10, the SMF-28e+ datasheet on p. 11, preparation task on pp. 11–12), `Documents/Entry - Sagnac effect experiment.pdf` (preparation sheet, Q1–Q12), Schwartz 2017; results in `RESULTS.md`, pipeline in `README.md`, code explanations in `docs/`, uncertainties in `uncertainty/REPORT.md`.

## The setup (manual, Fig. 1 and its caption)

| Manual | Status | Description | Notes |
|---|---|---|---|
| Laser, polarizer P₁, beam splitter, lenses L₁ and L₂ into the two fiber ends, P₂, L₃, camera; everything on the rotating table | done | As used on 11 October; `report_figures/setup.pdf` redraws the manual's Fig. 1 (without the small loop). | The polarizer orientations are not in the data files. |
| Large loop: R = 57.0 ± 0.5 cm, 13 ± 1 turns, centered on the axis | done, caveat | The values the analysis uses (`sagnac.py` constants, `uncertainty/instrument_inputs.json`). | Given, not measured by us. Whether ±1 turn is a standard uncertainty or a bound is not stated; N is the largest uncertainty of the comparison (7.7 %), and the slope fixes N = 13.7 ± 0.5. |
| Small loop: R = 15.5 ± 0.1 cm, 46 ± 1 turns (task 4) | not done | | Skipped on the instructor's instruction. |
| λ = 640 nm | done, caveat | Used as given. | No tolerance given and the laser model was not recorded: λ is left out of every uncertainty (0.16 % of a_th and c per nm). |
| Fiber SMF-28e (manual text) / SMF-28e+ (datasheet, preparation sheet) | done, caveat | Datasheet values: core 8.2 µm (a = 4.1 µm), NA 0.14, cable cutoff ≤ 1260 nm. | The two names differ; the report should say SMF-28e+ once with the note. At 640 nm the fiber is not single-mode (V = 5.6). |
| Coupling lens f = 8 mm (manual p. 10: at most 47 % into the fundamental mode) | not verifiable | | Coupled power and mode content were not measured. Single-mode coupling was attempted and found very hard to align (RESULTS, drift section). |

## Minimum experimental tasks (manual p. 12)

| Task | Status | Description | Notes |
|---|---|---|---|
| **1** Measure the (linear) ratio between the voltage on the motor and its angular velocity | done, caveat | `calibration.py` → `calibration.json/.txt`: 5 or 10 timed turns at nine settings S = 5.78–90; u(T) from the lap scatter and the 0.25 s reaction time of the presses; the line ω = kS + ω₀ weighted by the timing uncertainties: k = 0.02054 ± 0.00006 rad s⁻¹ per unit S, ω₀ = −0.0227 ± 0.0004 rad/s. Report 2.1, Fig. 2, Table 1. | The variable is the controller setting S, which sets the motor voltage; no voltage was read. The table turns slower than linear at high settings (the line is 3.0 % above the timed ω at S = 90); the line is used, as the task asks, and the bend beyond its rms is stated as a limitation. Calibration settings differ from the photo settings (only S = 90 is common). |
| **2** Spin the system at different speeds and picture the interference pattern at each | done | 161 photographs: rest (10) and nine settings S = 12.59–90 (15–22 each), 16:21–16:33. Report 2.2, Fig. 3 (one per setting), Fig. 4. | Taken in increasing order of S with rest photos only at the start, which makes a steady drift inseparable from the slope (task 5). |
| **3** Deduce the phase shift and compare with theory | done | `fringes.py` (ring finder, phase at the ring center), `sagnac.py` (setting means, regression): a = 1.83 ± 0.06 rad/(rad s⁻¹) against a_th = 1.7 ± 0.1, 0.6σ, p = 0.56. Report 2.2 (the phase of each photo, Fig. 5, Table 2) and 2.3 (the line and the comparison, Figs. 6–7). | The manual's Eq. 7 (minimum shift over local fringe period, averaged) is done in r², where every ring has the same period; this is also Schwartz's Eq. 5. The line is Schwartz's Eq. 7 (equal-weight regression). |
| **4** Repeat for the off-centered small circle | not done | | Instructor's instruction. The axis-independence it would test is Schwartz's result (his Fig. 4); the report can say so in one sentence. |
| **5** Evaluate the measurement errors | done | RESULTS Task 5 (budget, drift), `uncertainty/REPORT.md` (every component with its model and kind). Report 2.3 (the slope's three parts with their kinds) and 2.4 (what is not covered; Fig. 8), restated in brief in the conclusions. | Three parts of the slope's uncertainty (regression 1.3 %, method 2.9 %, calibration 1.3 %: the line's prediction uncertainty). Not covered by the quoted uncertainty: the 3.0 % bend at the top calibration setting beyond its rms, a steady drift in time, λ, the sign of the shift (rotation direction not recorded). |

## The manual's analysis instructions

| Manual | Status | Description | Notes |
|---|---|---|---|
| Eq. 5: Δφ = N (8π/(λc)) A·ω | done | a_th = 8π²NR²/(λc) = 1.738 ± 0.137 s (N, R). | λ and c without uncertainty (c exact; λ not given). |
| Eq. 6–7: fringe minima at y_i = Δφ + Λ i; Δφ from ⟨Δr_i/Λ_i⟩ over the minima | done, caveat | The fringes are rings evenly spaced in r², not straight fringes; the dark-ring radii give R_m² = P(m + η), φ_c = π − 2πη, the same average done in r² with P measured in every photo. | Explain in 1.1.2 why the local period differs between rings in r but not in r² (two Gaussian beams of different wavefront radius, QE3). |
| Eq. 9: V number | done | V = 5.6 at 640 nm (manual p. 10): several modes. | Used for the discussion of the wander (untested cause). |

## Preparation task in the manual (pp. 11–12)

These were answered before the experiment (pre-report) and are not part of the analysis.

| Question | Status | Where the report uses it |
|---|---|---|
| 1. Optimum NA and focal length of the coupling lens (D = 4 mm) | pre-report | One clause in 2.4's next-time list: the manual's matched 15.5 mm lens (93 % into the fundamental mode against 47 % with the 8 mm lens). |
| 2. Coupling efficiency of the RMS40X objective; the best of the four | pre-report | Not needed. |
| 3. Polarizers (and Berry phase) against unwanted reflections | pre-report | 1.2: why there are two polarizers. |

## Preparation sheet questions Q1–Q12

| Question | Status | Where it is answered or used | Notes |
|---|---|---|---|
| Q1 Interference of two plane waves at a relative angle | pre-report | 1.1.2 (manual Eq. 6) | |
| Q2 Distance between adjacent maxima | pre-report | 1.1.2: Λ = λ/(2 sin θ) | Our fringes are rings; Λ is only the plane-wave picture. |
| Q3 The Gaussian beam | pre-report | 1.1.2: transverse phase kρ²/(2R(z)) (QE3 PDF p. 10); envelopes in the cosine fit (QE3 PDF p. 7) | The basis of the ring model. |
| Q4 Derive Δφ = N(8π/(λc)) A·ω | pre-report | 1.1.1 (manual Eqs. 1–5; Pascoli; Tartaglia) | |
| Q5 Two ways to extract Δφ; which is more accurate when Δφ_max is large | done | 2.3: ring minima, cosine fit, template matching (method spread 2.9 %) and the fixed-pixel readout, all four in Fig. 7; the fixed pixels (the instructor asked for them) in their own paragraph, with the spread of the slope between 1000 random draws of ten pixels (11 %) as their uncertainty and the reason it is large at half a fringe, outside the method range | Our largest shift is half a fringe: the small-shift case. |
| Q6 Polarizers (Malus) and Berry phase | pre-report | 1.2 | |
| Q7 Lens–fiber distance, light from the fiber and from its surface | pre-report | 1.1.2: changes the wavefront radius, so P, not φ_c (to first order) | P varies by ±1.7 % between settings. |
| Q8 Influence of the center of rotation | done (theory) | 1.1.1: none (Schwartz Eq. 4) | Not tested here (task 4 skipped). |
| Q9 Influence of the wavelength on the pattern | pre-report | Not in the report (REPORT_TODO.md, "Left out on purpose") | |
| Q10 Influence of the refractive index (Schwartz) | done (theory) | 1.1.1: none (Schwartz Sec. 2) | |
| Q11 SMF-28e+ core diameter, NA, effective index at 1310 and 1550 nm | pre-report | 1.2 if wanted: core 8.2 µm, NA 0.14, effective group index 1.4674 (1310 nm) and 1.4679 (1550 nm) (datasheet in the manual, p. 11) | |
| Q12 Is SMF-28e+ single-mode at 640 nm? | done | 1.1.3: no, V = 5.6 > 2.405 (cutoff ≤ 1260 nm) | The plausible, untested cause of the wander (2.4). |

## Schwartz 2017 (assigned paper)

| Paper | Status | Description | Notes |
|---|---|---|---|
| Phase from the dark-fringe shift scaled by the fringe width, Eq. 5 | done | In r², with the width (P) measured in every photo | Schwartz notes that ring fringes have different widths in r; in r² they do not. |
| Slopes by linear regression with standard errors, Eq. 7 | done | Equal-weight regression of the ten setting means; slope error from the points | The line misses the setting means by five times their SEM (χ²/dof = 205/8), which is why the regression's own scatter is the right uncertainty. Figures show the data and the fits only, no residual panels. |
| Comparison with Eq. 1 | done | a/a_th = 1.05 ± 0.09 | Schwartz: 1.304 ± 0.038 against 1.344 ± 0.006 on his large ring. |

## Beyond the manual

| Item | Where | Why it is there |
|---|---|---|
| Speed of light from the slope, c = (2.9 ± 0.2) × 10⁸ m/s; N read the other way, 13.7 ± 0.5 | `sagnac.py`, RESULTS | The same comparison as a/a_th, shown as a speed; not an independent test. |
| Reaction time in the stopwatch uncertainty; the calibration line's prediction uncertainty for ω at the photo settings | `calibration.py`, `uncertainty/REPORT.md` 1–2 | The standard treatment of a timed calibration. |
| Ring model in r² and the ring finder | `fringes.py`, `docs/fringes.md` | The pattern is rings, evenly spaced in r². |
| Cosine fit, template matching, fixed pixels (1000 draws) | `fringes.py`, `methods.py`, `pixels.py` (24 October, with the current calibration) | Method uncertainty; prep Q5. |
| Drift check (joint fit with time; trends per setting in Fig. 8) | `sagnac.py`, `report_figures.py` | The settings were taken in increasing order. |
| Report figures and table blocks | `report_figures.py`, `REPORT_figures_tables.tex` | The report. |
| Tests: 5 at the top level | `tests/`, `docs/tests.md` | Synthetic checks of template matching and of the fixed-pixel fit, of the report-form rounding and of Eq. (5) for c. |

## Open points

- λ's tolerance and the laser model: ask the instructor or read the laser's label, if the report should include λ in the budget.
- N: if the instructor knows how 13 ± 1 was obtained (counted? estimated from length?), say so; it is the dominant uncertainty.
- A steady drift and the sign of the shift cannot be checked with these data; only new measurements (rest photos between speeds, reversed rotation, random order) can.
- The presentation code (`presentation_slides/source/`, `sagnac_update/`) is kept as the record of 19 October and needs the earlier calibration files from `previous_reports/` to rerun (`README.md`, "Earlier stages"). Its random-pixel scan was superseded on 24 October by `pixels.py`, which repeats it with the current calibration (same numbers to within the rounding).
