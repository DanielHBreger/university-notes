# Lab 5 report (Sagnac): what is left (25 October 2026, after draft (17))

Numbers from `sagnac.txt` and `calibration.txt`; `RESULTS.md` has the detail behind each. The full plan with its history is archived in `previous_reports/REPORT_TODO_full_2026-10-25.md`; wording-level items (typos, phrasing, names, caption nits) were dropped from this file on 25 October and are in that archive's draft (12) lists if ever needed. Reviews of new drafts go at the top of this file, short, newest first.

Each open item has a code block: where it goes, then a wording suggestion. The suggestions fall under the own-words item below: reword them, do not paste them.

Figure numbers in this file are draft (17)'s: Fig. 1 calibration, Fig. 2 patterns, Fig. 3 extraction, Fig. 4 drift, Fig. 5 phase against ω, Fig. 6 readouts. The setup figure (`report_figures/setup.pdf`, for 1.2; its block with the label `fig:setup` is at lines 25–37 of `REPORT_figures_tables.tex`) is not in the report yet; once it is, it becomes Fig. 1 and the others shift, so cite every figure by its label.

Section 1, Hila's Word draft (27 October): checked against the manual (Eqs. 1–7, the p. 12 figure and its values), Schwartz (Eqs. 4–5) and `fringes_frames.csv`; no physics errors. Converted to `REPORT_section1.tex` (her wording; in the equations λC → λc twice and the 2π of the manual's Eq. 7 added, which the manual drops: without it the average is in fringes). It covers the Task 4 blocking item. Open: "Section 1: open items" below.

Draft (20) (25 October): the abstract written (purpose, method, result, verdict, limits; the drift left out, your call); the acknowledgment thanks Peck and Yakobovich. Open: one note on the abstract, the three on the conclusions, the one-liners, Section 1 and the references.

Draft (19) (25 October): the conclusions written, numbers right; the acknowledgment rewritten (done); "read". Open: the one-liners below, three notes on the conclusions, and the abstract, Section 1 and the references.

Draft (18) (25 October): 2.1 ¶1 rewritten (consecutive turns, a press per turn, first and last press; Table 1 and Fig. 1 cited); the multimode sentence right; "fixed pixels" in 2.3 ¶3 and the method in ¶6. Still open: the items below.

Draft (17) (25 October): 2.4 dissolved in the text (heading, roadmap clause and three forward references gone; the wander with its ±8 mrad/s, the multimode cause and the drift figure in 2.2 ¶4; the two not-covered items named in 2.3 ¶7); the S = 90 line stated in 2.1 ¶3; "linearly with r²"; the 0.52 → 3.3 rad link and the ring numbering in 2.2 ¶3; 2.3 ¶1's intercept sentence right; ¶2's regression number 0.024; ¶5's "left out of the range"; P10's a_th/a inline. Still open: everything below.

## Rules that still bind
- The prose is yours; every number with its uncertainty and kind, one significant figure, comparisons from the unrounded values; a number gets its ± once, complete.
- Methods presented as the lab documents' methods and cited where used; template matching and the drift check marked as ours; the half-range readout uncertainty "taken as".
- One unit for the slope everywhere; figures: data and fits only; each number in one form.

## Blocking before submission
- [ ] All prose in your own words. About twenty sentences of 2.1–2.3 and the captions are light edits of the archived plan's suggested sentences (2.1's misfit sentence and last paragraph; 2.2 ¶2 from the 241/132 px sentence on, ¶3's slope–intercept sentence, ¶4's uncertainty sentences; 2.3 ¶2, ¶5, ¶7; the Fig. 1, 3, 5 and 6 captions). Hila's Section 1 follows the Section 1 blocks below sentence by sentence in the same way.
- [ ] Every number checked against `sagnac.txt` and `calibration.txt`; the slope's ± once, 1.83 ± 0.06.
- [x] Acknowledgments: done in draft (19); confirm the models from your records.
- [ ] References and in-text citations (below); Section 2 has none yet, and the placeholders "[cref …]", "[section 1.2]", "[subsection 1.2]", "[cite entry report]", "[cite]" (Eq. 7) and "[cref to manual eq 5]" are still in the text.
- [ ] Task 4 stated as not done, on the instructor's instruction (1.2).
- [ ] Every figure and table cited in the text: Fig. 4, the drift figure, is cited nowhere (cite it in 2.2 ¶4 at the ±8 mrad/s).
- [ ] The slope's unit: rad/(rad/s) in the text against rad/(rad s⁻¹) in the figure axes; one form (change the text, or the labels in `report_figures.py`).

## Section 2: open items (draft (17))

Dropped on 25 October at your call: the bend's size in 2.1 ¶2; the drift check's numbers in 2.3 ¶7 (your sentence stays); the fixed pixels' attribution and low mean; the reasons for the readout offsets in 2.3 ¶4; the η range in 2.2 ¶3; why the intercept is left free in 2.3 ¶1.

- [ ] What is left, one line each:
  - 2.1 ¶1: "the reaction time's effect is reduced" → "enters as √2 × 0.25 s divided by the number of turns". Why: without the factor the 1.1 % at S = 90 cannot be reproduced.
  - 2.1 ¶3: at rest ω = 0 by definition. Why: Table 2 shows 0 while the line gives −0.023.
  - 2.2 ¶2: "holds true to within 2 %" → "agree to 2 % of a fringe (rms)". Why: "to within" claims a bound the data break at 7.5 %.
  - 2.2 ¶4: say that Table 2 lists the mean over a setting's photos. Why: "each setting's mean" is never defined.
  - 2.3 ¶2: "contributes 0.024 rad/(rad/s) from each setting". Why: as written, each setting contributes 0.024.
  - Table 2 and Fig. 5 captions: say what the uncertainties and error bars are (ω 1.3 %, common to all settings; Δφ 0.04 rad for every setting, the scatter about the line of 2.3; the rest row the reference; the date). Why: the two uncertainties are of different kinds, and that difference is what 2.3 builds on.

## Section 1: open items (Hila's draft, 27 October)

- [ ] 1.2 ¶3: ten turns were timed only at the four fastest settings, five at the five slowest. Why: Table 1 shows n = 5 at five settings.
- [ ] 1.2 ¶2: the fiber guides several modes at 640 nm (V = 5.6, manual Eq. 9 and p. 10). Why: 2.2 ¶4 cites Section 1.2 for it.
- [ ] Citations throughout, and "after the manual's figure on p. 12" in the setup caption. Why: the figure is redrawn from it; the blocks below mark each spot.
- [ ] 1.1.2 ¶1: "consistent beams" → coherent.
- [ ] 1.1.2 ¶3: rotation adds the same phase across the whole beam. Why: that is why only φ_c changes.
- [ ] 1.1.2 ¶3: "none of these effects changes φ_c" → "mainly", as RESULTS.md says. Why: a tilt also changes the phase at the moved center, and our ring center moves 9 px by the top speed.
- [ ] 1.1.2 ¶5: 157 ± 1, 159.2 ± 0.3 and 242.0 ± 0.2 px (standard errors of the setting means).
- [ ] 1.1.1 ¶2: "after beam 1" → after the counter-rotating beam. Why: beam 1 is defined only in the setup caption.
- [ ] 1.1.2: R₁, R₂ (beam curvatures) clash with the ring radii R_m at m = 1, 2.

## Section 1: content to write

### 1.1 Theoretical background (two parts; no fiber-modes subsection)

1.1.1 The Sagnac phase, ending in the coefficient the experiment measures.

```
Where: 1.1.1, first paragraph.

Consider a loop of radius R that turns at angular velocity ω, with light sent around it in
both directions from a beam splitter that moves with it at v = ωR. In the laboratory frame
the light travels at c both ways, so the beam running with the rotation needs
t₊ = l/(c − v) to return to the splitter, which has moved on, and the counter-running beam
t₋ = l/(c + v), with l = 2πR [manual, Eqs. 1–2; entry report, Eqs. 12–16]. Their
difference is Δt = 2lv/(c² − v²) ≈ 4Aω/c² for v ≪ c, where A = πR² is the area enclosed by
the loop; a fiber wound N times gives N times that.
```

```
Where: 1.1.1, second paragraph.

A delay Δt between two beams of wavelength λ is a phase difference 2πcΔt/λ [manual,
Eq. 3], so the interferometer measures
    Δφ = 8πNAω/(λc),                                                    (eq:sagnac)
the Sagnac formula for N turns [manual, Eq. 5; entry report, Eq. 19; Schwartz, Eq. 1].
```

```
Where: 1.1.1, third paragraph.

The shift does not depend on where the rotation axis lies, only on the enclosed area
[Schwartz, Eq. 4], nor on the refractive index of the fiber, which rotates with the
interferometer [Schwartz, Sec. 2]; its sign follows the sense of rotation.
```

```
Where: 1.1.1, last paragraph.

The phase shift is therefore a straight line through the origin against ω, with slope
a_th = 8π²NR²/(λc). This coefficient, fixed by the loop and the wavelength, is what the
experiment tests.

Optional: At our top speed it corresponds to a delay of about 1 fs, half an optical
period, which interference resolves.
```

1.1.2 Interference and the ring pattern: why rings, why evenly spaced in r², how a photo becomes a phase.

```
Where: 1.1.2, first paragraph.

Two coherent beams of intensities I₁ and I₂ that meet with a phase difference φ give
I = I₁ + I₂ + 2√(I₁I₂) cos φ [manual, Eq. 6], dark where φ = π + 2πm. Two plane waves
crossing at an angle give straight fringes of a single period [manual, Eq. 6 and Fig. 2];
our pattern is a set of rings, as in Schwartz's apparatus [Schwartz, p. 6].
```

```
Where: 1.1.2, second paragraph.

The rings come from the beams' curvature. A Gaussian beam carries a transverse phase
kr²/(2R(z)), with R(z) the radius of curvature of its wavefront [QE3; entry report,
Eq. 11], so two beams on a common axis with different curvatures at the camera differ by
    φ(r) = φ_c + 2πr²/P,   P = 2λ/|1/R₁ − 1/R₂|,                         (eq:rings)
a phase that grows linearly with r²: the rings are evenly spaced in r², like Newton's
rings.
```

```
Where: 1.1.2, third paragraph.

Dark rings lie where φ = π + 2πm, that is at
    R_m² = P(m + η),   φ_c = π − 2πη,                                    (eq:dark_rings)
so the squared radii of the dark rings lie on a straight line against the ring number,
with slope P and an intercept that gives the phase at the center. Rotation adds the same
phase across the whole beam, so it changes only φ_c: the rings slide inward or outward
while P stays. A tilt between the beams moves the ring center and a change of the
lens-to-fiber distance changes P; neither moves φ_c.
```

```
Where: 1.1.2, last paragraph.

Reading the phase from the dark-ring radii in this way is the manual's Eq. 7, the shift of
each minimum over its own period averaged over the minima [manual, Eq. 7], carried out in
r², where every ring has the same period P; Schwartz reads his circular fringes the same
way, with the shift and the width measured in r [Schwartz, Eq. 5]. A reading over a local
width is exact only while the width stays the same over the shift, which ours does not at
half a fringe: the first dark ring moves by 157 px between rest and the top speed while
its gap to the next ring changes from 159 to 242 px. In r² the period is one number for
every ring and every photo (checked in Section 2.2), and because η is a ratio of two areas
in pixels, no pixel-to-length calibration is needed.
```

### 1.2 Experimental setup and procedure (with the setup figure, not yet in the report: `report_figures/setup.pdf`, column width, label `fig:setup`; it becomes Fig. 1 and renumbers the others)

```
Where: 1.2, first paragraph, with the setup figure (\cref{fig:setup}).

\cref{fig:setup} shows the setup, after the manual's Fig. 1 [manual]. A 640 nm laser passes
a polarizer P₁ and a beam splitter, and the two beams are coupled by lenses L₁ and L₂ into
the two ends of a fiber wound N times in a loop of radius R, so that each travels the loop
in the sense opposite to the other. They return through the splitter, a second polarizer
P₂ and the camera lens L₃ to a camera. The polarizers suppress reflections from the
lenses, holders and splitter that did not travel the loop [preparation sheet, Q6]. The
whole interferometer sits on a turntable, with the loop centered on the rotation axis,
and turns with it.

Caption (fig:setup): Top view of the interferometer, after the manual's Fig. 1 [manual].
P₁, P₂ polarizers; L₁, L₂ coupling lenses; L₃ camera lens. Beams 1 and 2 travel the same
fiber loop (N turns, radius R) in opposite directions and meet at the splitter. Everything
turns with the table about the loop's center; the direction of ω is schematic, as it was
not recorded. The small loop of task 4 is not drawn.
```

```
Where: 1.2, second paragraph: the given values, the fiber (now required, since 2.2 ¶4
names the modes), then a_th from them (2.3 cites this value).

The loop has N = 13 ± 1 turns of radius R = 57.0 ± 0.5 cm [manual, Fig. 1 caption]; the
wavelength is given as 640 nm without a tolerance. The fiber is Corning SMF-28e+, which at
640 nm guides several modes (V = 5.6 [manual, Eq. 9; entry report, Eq. 30]); the 8 mm
coupling lenses couple at most 47 % of the light into its fundamental mode, against 93 %
for the matched 15.5 mm lens [manual, p. 10]. With these values, Eq. (sagnac) predicts
a_th = 8π²NR²/(λc) = 1.7(1) rad/(rad/s); the uncertainty is that of the given N (7.7 %)
and R (1.8 %, as it enters squared), and λ, having no stated tolerance, contributes none
(a_th moves by 0.16 % per nm).
```

```
Where: 1.2, third paragraph: the procedure.

The table is driven by a motor whose controller setting S, a number without unit that
sets the motor voltage, was our only handle on the speed. To calibrate it we timed ten
consecutive turns (five at the five slowest settings) at nine settings between S = 5.78
and 90 with a stopwatch, one press per turn. We then photographed the fringe pattern 161
times: ten at rest and 15 to 22 at each of nine settings between S = 12.59 and 90, taken
in increasing order of S between 16:21 and 16:34 on 11 October 2026, with the rest photos
only at the start. The direction of rotation was not recorded.
```

```
Where: 1.2, last paragraph.

The manual's task 4, repeating the measurement with a small off-center loop
(R = 15.5 ± 0.1 cm, 46 ± 1 turns), was not done, on the instructor's instruction;
Schwartz (2017) reports that test.
```

No analysis-chain paragraph: the Section 2 roadmap is that chain.

## Abstract

Written in draft (20), 115 words; it states the purpose, the method, the result with its numbers, the verdict and the limits, with no references. One note:
- "Two laser beams travel through the loop in opposite directions causing a phase difference … to grow proportionally to the angular velocity" makes the beams the cause; the rotation is. Say that the rotation gives the two beams a phase difference proportional to ω, and that the loop rotates with the whole interferometer (the point of the test). Why: the one sentence that states the effect should state it right.

## Conclusions (the template's Summary)

Written in draft (19); every number matches Section 2. Three notes, one line each:
- "We were able to convert this to a measurement of the speed of light" claims more than P10, which says c/c₀ is the slope comparison read the other way, not an independent test; say "read as a speed of light" or "which, read the other way, gives". Why: the two statements contradict each other as written.
- The reversed rotation needs its reason: the Sagnac shift changes sign with the direction while a drift does not. Why: without it the reader cannot see why it tests the drift, and the rubric scores improvements on being specific.
- "13 ± 1" against the "(1)" notation used everywhere else, in the same paragraph as 1.83(6). Why: one notation per report.

## Acknowledgments

Done in draft (19): developer, models and tool; what the AI helped with (the code, the figures, the tests, some analysis choices, the outline, the drafts); your checks; the text and the responsibility yours. Left: confirm the models from your records (by the commit attributions: Fable 5.1 on 22–23 October, Opus 5.5 on 24 October), and name ChatGPT if it touched the Sagnac work.

## References (checked against `Documents/`)
1. Sagnac experiment, laboratory manual, Physics Lab 5, Technion (with the Corning SMF-28e+ sheet). `Documents/Sagnac_manual.pdf`
2. Sagnac effect: preparation sheet, Physics Lab 5, Technion. `Documents/Entry - Sagnac effect experiment.pdf`
3. E. Schwartz, "Sagnac effect in an off-center rotating ring frame of reference," Eur. J. Phys. 38, 015301 (2017), doi:10.1088/0143-0807/38/1/015301.
4. G. Pascoli, "The Sagnac effect and its interpretation by Paul Langevin," C. R. Physique 18, 563–569 (2017): only if cited in 1.1.
5. A. Tartaglia and M. L. Ruggiero, "The Sagnac effect and pure geometry," Am. J. Phys. 83, 427–432 (2015), doi:10.1119/1.4904319: only if cited.
6. J. T. Costello, "The Gaussian beam," lecture slides (QE3). `Documents/QE3.pdf` (the phase of a Gaussian beam: PDF p. 10)
7. D. H. Breger and H. Hemo, Sagnac effect entry report, 24 October 2026: our answers to the preparation sheet; cite for Eq. 11 and Eqs. 26–27.

```
Where: the citations in Section 2 (Section 1's are in its blocks).

2.1 ¶2  "the linear ratio task 1 asks for"                      → \cite{manual} (p. 12)
2.2 ¶2  "Gaussian beams with different wavefront curvatures"    → \cite{qe3, entry} (Eq. 11)
2.2 ¶3  "Equation 7 from the experiment manual[cite]"           → \cite{manual}; add
        "Schwartz reads his circular fringes the same way"      → \cite{schwartz} (Eq. 5)
2.2 ¶4  "the fiber guides several modes at 640 nm"              → Section 1.2 (manual Eq. 9)
2.3 ¶1  the equal-weight regression                              → \cite{schwartz} (p. 7)
2.3 ¶3  "Our entry report proposed two"                          → \cite{entry} (Eqs. 26–27)
2.3 ¶3  "intensity envelopes of the two Gaussian beams"          → \cite{qe3}
2.3 P10 "[cref to manual eq 5]"                                  → \cite{manual}
```

## Numbers (unrounded, for comparisons)
Slope 1.8258 ± 0.0236 (regression) ± 0.0539 (readout) ± 0.0239 (calibration) = ± 0.0635 (3.5 %); b = −0.047 ± 0.027 (1.7σ); s = 0.0444 rad; R² = 0.9987; cosine 1.8626, template 1.7549; fixed pixels: main draw 1.7734 ± 0.0367, draws mean 1.689, SD 0.193 (11.4 %), median fit error 0.057; 20/40/80 pixels SD 5.3/4.2/3.2 %; a_th = 1.7381 ± 0.1371; a/a_th = 1.050 ± 0.091, z = 0.58, p = 0.56; c = 2.854 ± 0.246 × 10⁸ m/s (parts 1.3, 2.9, 1.3, 7.7, 1.8 %); c/c₀ = 0.952 ± 0.082. Calibration: k = 0.02054 ± 0.00006, ω₀ = −0.0227 ± 0.0004, χ² = 20.2/7, p = 0.005, the line 3.0 % above the timed ω at S = 90, scatter 1.28 % rms; bend: slope 1.8658 (quadratic, +2.2 %) or 1.8723 (segments, +2.5 %); ω² term +0.023 ± 0.051; drift fit a = 2.12 ± 0.17, D = −0.82 ± 0.47 mrad/s, correlation −0.992; η 0.62 at rest, 0.10 at the top speed; the two ring gaps in r² agree to 2 % of a fringe (rms) over the 126 photos with three rings, the largest 7.5 %.
