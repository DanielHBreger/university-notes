# Sagnac effect: presenter notes

The main talk is slides 01–14 and runs about 10 minutes. Slides 10 (random pixels) and 12 (drift) are marked as **optional**. Skip them if time is short. B1–B3 are backup slides for questions.

**Where the explanations come from.** Every formula and claim below is tied to one of the project documents, cited in brackets:

- **[Manual]** Sagnac_manual.pdf (lab manual): equation numbers as printed, the V-number section, the minimum tasks.
- **[Prep]** "Entry - Sagnac effect experiment.pdf": preparation questions Q1–Q12.
- **[QE3]** QE3.pdf: Gaussian beam lecture slides.
- **[Schwartz]** E. Schwartz, Eur. J. Phys. 38, 015301 (2017): equation numbers as printed.
- **[Tartaglia]** A. Tartaglia and M. L. Ruggiero, Am. J. Phys. 83, 427 (2015).
- **[Pascoli]** G. Pascoli, C. R. Physique 18, 563 (2017).
- **[Analysis]** Our analysis documentation: `sagnac_update/README.md`, `sagnac_update/RESULTS.md` and the text outputs `sagnac.txt`, `calibration.txt`, `methods.txt`.

**Uncertainty notation.** A number like 0.02039(4) means 0.02039 ± 0.00004. The digit in parentheses is the uncertainty in the last digit shown. All uncertainties are one standard deviation.

**χ²/dof and p (used on several slides).** For a fit, χ² = Σ (data − fit)² / u², summed over the points, where u is each point's uncertainty. If the uncertainties are right and the model fits, χ² is about equal to the number of degrees of freedom, dof = number of points − number of fitted parameters, so χ²/dof ≈ 1. The p value is the probability of getting a χ² at least this large if the model and the uncertainties were correct. A tiny p (for example 10⁻⁸) means the model, or the error bars, are not enough to describe the data.

---

## 01 Title (0:20)

"We measured the Sagnac effect: light sent both ways round a fibre loop on a rotating table. From the size of the effect we get the speed of light."

The five circles are real photographs from our run, each centred on its fitted ring centre. From left to right the table turns faster, and the dark ring closes in on the centre.

---

## 02 Rotation delays the beam that travels with it (1:00). Let the animation run once.

**Say:** "A beam splitter on the table sends light both ways around a closed loop. While the light travels, the splitter moves with the table. The beam going with the rotation (blue) has to chase the splitter, so it needs a little more than one full turn. The beam going against the rotation (orange) meets the splitter early. The bar at the bottom shows the two return times and the gap Δt between them."

**The formulas, step by step:**

1. **Travel times.** In the lab frame both beams move at c. The splitter moves at v = ωR. The co-rotating beam must cover the loop length l plus the distance the splitter moved, so t₊ = l/(c − v). The counter-rotating beam needs t₋ = l/(c + v). [Pascoli Eq. 1, Manual "Classical derivation"]
2. **Delay.** Δt = t₊ − t₋ = 2lv/(c² − v²). For a circle, l = 2πR and v = ωR. Since v ≪ c this becomes Δt ≈ 4Aω/c² with A = πR². With N turns the delay is N times larger: Δt = 4NAω/c². [Manual Eq. 1–2, Tartaglia Eq. 1–2, Pascoli Eq. 2, Schwartz Eq. 2–3]
3. **Phase.** The camera sees a phase difference equal to the light's angular frequency 2πc/λ times the delay: Δφ = (2πc/λ)Δt = 8πNAω/(λc). This is the Sagnac formula, the equation we test. [Manual Eq. 3–5, Schwartz Eq. 1, Pascoli Eq. 3]

**Grey box (predicted for our loop at the top speed, ω = 1.81 rad/s).** Δt = 4·13·(1.02 m²)·(1.81 rad/s)/c² = 1.07 fs. One optical period at 640 nm is λ/c = 2.13 fs, so the delay is half a period and the phase is about π ≈ 3.2 rad, half a fringe.

**Not to scale.** The animation draws the rim moving at 0.25c so the effect is visible. With that speed the return times are 2π/(1 − 0.25) and 2π/(1 + 0.25) in units of R/c. In the experiment ωR is only about 1 m/s, about 3×10⁻⁹ c.

**Footer claims.** The phase does not depend on where the rotation axis is, only on the enclosed area. Schwartz shows this by splitting the velocity into a rotation about the loop's centre plus a constant vector, whose closed-loop integral vanishes [Schwartz Eq. 4, Prep Q8]. The phase is also independent of the fibre's refractive index [Schwartz Sec. 2, refs 19–20, Prep Q10].

---

## 03 Both beams travel the same 13-turn fibre loop (0:45)

**Say:** "The laser goes through a polarizer to the beam splitter. One beam is coupled into the fibre at L₁, the other at L₂. They travel the same 13 turns in opposite directions, recombine at the splitter and go to the camera. Everything you see turns with the table. We measure two things independently: ω with a stopwatch, and Δφ from photos of the fringes."

- The drawing follows the apparatus figure of the lab manual. The polarizers reduce unwanted back-reflections from the fibre ends [Prep Q6, Manual prep question 3].
- System parameters N = 13 ± 1, R = 57.0 ± 0.5 cm and λ = 640 nm were given. The area is A = πR² = 1.02 m².
- **Prediction.** Equation 5 is linear in ω, Δφ = a_th·ω, with slope a_th = 8πNA/(λc) = 8π²NR²/(λc) = 1.74 ± 0.14 rad/(rad/s) [Manual Eq. 5 with A = πR²]. The unit rad/(rad/s) is radians of phase per rad/s of rotation (it equals seconds).
- The asterisk: "Most of this ±0.14 comes from the turn count. We come back to that with the speed of light."

---

## 04 The rings are evenly spaced in r², not in r (0:50)

**Say:** "The camera does not see straight fringes. It sees rings, and they are not evenly spaced in radius: 241 pixels between the first two dark rings, 132 between the next two. Plotted against r², the spacing is the same, 1.08 × 10⁵ px² both times."

**Why rings, and why r².**
- Two plane waves meeting at an angle give straight, evenly spaced fringes with period Λ = λ/(2 sin θ) [Manual Eq. 6, Prep Q1–Q2]. Our beams are not plane waves. They leave the fibre as Gaussian beams.
- A Gaussian beam has the transverse phase kρ²/(2R(z)), where R(z) is the wavefront radius of curvature [QE3, "The phase of the Gaussian Beam" slide, PDF p. 10, Prep Q3].
- The two output beams reach the camera with different curvature radii R₁ and R₂. Their phase difference is then (k r²/2)(1/R₁ − 1/R₂) plus a constant. Writing k = 2π/λ gives φ(r) = 2π r²/P + φ_c with P = 2λ/|1/R₁ − 1/R₂| [Analysis README, "Fringe phase"].
- So the phase grows linearly in r², and the dark rings are equally spaced in r² with period P. Schwartz makes the same observation: the pattern is not linear like Young's fringes, and fringe widths differ between the centre and the edge [Schwartz Sec. 3].
- **What moves the rings.** Rotation adds the same phase across the whole beam. That changes only φ_c, so the rings slide inward or outward. A relative tilt of the beams adds a linear phase, which only moves the ring centre. A change of a lens–fibre distance changes the curvature, so it changes P but not φ_c [Analysis README, Prep Q7]. This is why we track φ_c.

**Figure.** Left: the fastest-rotation photo with the fitted dark rings (dashed). Right: the intensity averaged over the fringe sector against r (top) and against r² (bottom). The lower-left part of the image has no fringes, so the averaging uses −75° to 180° only.

---

## 05 A stopwatch calibration converts motor setting to ω (0:50)

**Say:** "The motor controller gives a setting S, not an angular velocity. At nine settings we timed five to eight full turns with a stopwatch, converted the mean turn time to ω = 2π/T and fitted a straight line, as the manual asks." [Manual minimum task 1]

- **Fit result.** k = 0.02039(4) rad/s per unit of S, that is 0.02039 ± 0.00004. The intercept ω₀ is part of the fit but we do not quote it. Below the lowest timed setting the table does not turn at all (friction stops it), so the line means nothing there.
- **χ²/dof = 47.2/7, p = 5×10⁻⁸.** The line does not describe the points within their error bars. The residual panel shows why: the points curve, and at the top speed the table turns 2.3 % slower than the line predicts. The manual asks for a linear relation, so the linear fit is what we use in the rest of the talk.
- **Error bars, relative.** Both numbers are relative uncertainties, u(ω)/ω.
  - *Stopwatch timing.* The turns were taken as stopwatch splits, so neighbouring laps share an endpoint. The mean turn time is (t_n − t₀)/n with uncertainty u(T) = SD/√(n(n+1)) [Analysis README, "Calibration"]. Since ω = 2π/T, u(ω)/ω = u(T)/T. Example at S = 90: T = 3.545 s, u(T) = 0.006 s, so u(ω)/ω = 0.18 % and u(ω) = 0.003 rad/s. Over the nine settings this is 0.04–0.18 %.
  - *Table speed scatter ε = 0.48 %.* If the table is set to the same S again, its real speed differs by about 0.48 % of ω. We do not measure this directly. It is the extra scatter of the settings about a smooth curve that the stopwatch alone cannot explain [Analysis README, "Calibration"]. It is added in quadrature to every error bar and later to every ω.
- The dashed vertical lines are the nine settings used for the photographs.

---

## 06 Faster rotation pulls the rings inward (0:40). Let the animation play once.

**Say:** "These are all 161 photographs in the order we took them, 0.1 s each. Watch the inner dark ring close in on the centre. On the right each photo adds a point: its ω and its phase. Large dots are the mean at each setting. At the end the fitted line and the prediction appear."

- The bar at the bottom is the recording clock (16:21 to 16:33, 12 minutes). The settings were taken in increasing speed. This matters on slide 12.
- The phase is measured at the ring centre, relative to the mean at rest. The next slide shows how.
- The photos in the animation are lightly smoothed for display only. The analysis uses the original photos.
- **Right axis.** "fringes" is Δφ/2π. At the top speed the pattern has moved by about half a fringe, as predicted on slide 2.

---

## 07 Squared ring radii give the phase at the centre (0:50)

**Say:** "For every photo we find the dark rings, fit concentric circles to get the centre and the radii, and plot R² against the ring number. That gives a straight line. Its intercept tells us how far the innermost ring is from the centre, in units of one fringe."

**Formulas.**
- From slide 4, φ(r) = 2π r²/P + φ_c. A dark ring is where the two beams cancel, φ = π + 2πm [Manual Eq. 6 with cos = −1, Prep Q1]. Solving for r² gives R_m² = P(m + η) with η = ½ − φ_c/2π [Analysis README, `fringes.py`].
- So a straight line of R_m² against m has slope P and intercept Pη. Rearranging, φ_c = π − 2πη.
- **This pair.** η goes from 0.62 at rest to 0.10 at the top speed. The phase change is Δφ = 2π × 0.52 ≈ 3.3 rad.
- This is the manual's method (shift of each minimum divided by the local fringe period, averaged over the minima) [Manual Eq. 7] and Schwartz's Δφ = 2π Δx/δ [Schwartz Eq. 5], done in r². In r² every ring has the same period P, so one straight-line fit uses all the rings at once.
- **Why it is robust.** P is fitted for each photo, so a change of the ring spacing does not move φ_c. Only R²/P enters, so no pixel-to-length calibration is needed. The positions of minima also do not depend on the camera's brightness response [Analysis README].
- The circle fits have rms 3–8 px over all photos.

---

## 08 The phase grows linearly with ω, as Sagnac predicts (1:00)

**Say:** "Each point is the mean phase of one setting, against its calibrated ω. A weighted straight line gives a slope of 1.84 ± 0.02 rad per rad/s. The prediction from the loop geometry is 1.74 ± 0.14. They agree: the difference is 0.7 standard deviations."

- **Fit.** Δφ = aω + b, weighted least squares. Each setting's uncertainty is u = √(SEM² + (aεω)² + σ_x²) [Analysis README, "Phase against omega"]:
  - SEM is the standard error of the mean of that setting's photos (photo-to-photo scatter, 0.004–0.028 rad).
  - aεω converts the 0.48 % table-speed scatter of slide 5 into phase.
  - σ_x is explained next.
- **Why σ_x.** With only SEM and the table-speed term, χ²/dof = 56.6/8 and p = 2×10⁻⁹. The settings scatter about the line more than these errors allow. Something adds scatter between settings that the photos within one setting do not see. Slide 12 shows the phase slowly wandering in time, which is such a source. We model it as one extra uncertainty σ_x, the same for every setting, added in quadrature. Its size is chosen so that χ²/dof = 1, which gives σ_x = 0.041 rad. This means the slope's uncertainty reflects the actual scatter of the points about the line [Analysis README, `sagnac.txt`]. With the photo-to-photo scatter alone the slope would be 1.822(6), an uncertainty that is too small.
- **Slope uncertainty.** 0.024 from the fit plus 0.004 from the calibration line, combined: a = 1.84 ± 0.02 rad/(rad/s).
- **Comparison.** a/a_th = 1.06 ± 0.08. The difference in units of the combined uncertainty is z = (a − a_th)/√(u_a² + u_th²) = 0.7, two-sided p = 0.47 [`sagnac.txt`]. For scale, Schwartz's large ring gave 1.304 ± 0.038 against a prediction of 1.344 ± 0.006, a 3–8 % agreement [Schwartz Eq. 6–7].
- **Figure.** In the upper panel all error bars are smaller than the markers. The lower panel shows the residuals with their error bars. The rest phase is set to zero for display.

---

## 09 Four ways of reading the phase agree within ±3 % (0:50)

**Say:** "Does the result depend on how we read the phase from the photos? We applied three more methods to the same photographs."

- **Cosine fit** [Prep Q5 asks for two ways to extract Δφ]. The intensity against s = r² is fitted with C e^(−s/s_C) + B e^(−s/s_B) cos(2πs/P + φ) + D. The exponentials are the envelopes of two Gaussian beams of different size [QE3, intensity slide, PDF p. 7], D is the dark level, φ is the phase at the centre [Analysis README]. It uses every pixel, but relies on the envelope model and on a linear camera response.
- **Template matching.** Each r² profile is shifted and scaled to match the averaged rest profile. The shift D in r² gives the phase 2πD/P [Analysis README]. It uses no fringe model, but reads the phase in the middle of the comparison window rather than at the centre.
- **10 random pixels** (instructor's method). Ten pixels at fixed image positions are followed through all 161 photos (next slide).
- **Left panel.** Each method's phase change minus the ring-minima phase change, per setting. The cosine fit reads 0.1–0.19 rad higher at every moving setting. The difference comes from the rest photos, where the bright flat centre does not look like a cosine (fit rms 9.5 counts at rest, against 1.4–2.3 at high speed) [`sagnac.txt`].
- **Right panel.** The fitted slopes 1.84, 1.88, 1.77 and 1.79, with fit and calibration uncertainty. The grey band is the prediction ± 1σ (1.60–1.88). All four slopes lie inside it.
- **Method uncertainty (explain here, not on the slide).** All four methods read the same photographs, so their differences are not random noise. They show how much the answer depends on the choice of method. We take half the range of the slopes, (1.877 − 1.768)/2 = 0.055 rad/(rad/s), about 3 %, as a systematic uncertainty [Analysis README]. It is combined in quadrature with the fit and calibration uncertainty: √(0.024² + 0.055²) = 0.06. This gives the bold line, a = 1.84 ± 0.06 rad/(rad/s). The random-pixel slope lies inside the range of the other three, so including it does not change this number.
- **Footer.** Each method's slope comes from the same weighted line as slide 8, each with its own σ_x. Without σ_x every method's line is rejected (p < 10⁻⁶), so every method sees the same extra scatter between settings.

---

## 10 The slope from ten random pixels depends on the draw (0:40, optional)

**Say:** "In the random-pixel method we pick ten pixels in the fringe region and follow their brightness through all photos. As the rings move, each pixel gets brighter and darker. One fit finds the phase of every photo from all ten pixels together."

- **Model.** I_jk = a_j + b_j cos(ψ_j + δ_k): pixel j in photo k has an offset a_j and an amplitude b_j, ψ_j = 2π s_j/P₀ is its fringe phase at rest and δ_k is the photo's phase change. An alternating least-squares fit finds all a_j, b_j and δ_k [Analysis README, `methods.txt`]. This is the same cosine law of two-beam interference as Manual Eq. 6.
- **The histogram.** We repeated the method for 31 random draws of ten pixels. The slope ranges from 1.41 to 2.12, with the middle 68 % between 1.51 and 1.90. The draw used on the previous slide gave 1.79.
- **Shortcomings.**
  1. A fixed pixel reads the phase at its own radius. If the ring spacing P changes (the wavefront curvature changes, slide 4) or the ring centre moves, the local phase changes too, and the pixel cannot tell this from a Sagnac shift. The purple outline repeats the fit with each photo's own fitted centre and P. The spread narrows (68 %: 1.68–1.98) but does not disappear [`methods.txt`].
  2. Over the whole run the pattern moves only half a fringe (the fitted phases span 3.2 rad). Ten pixels that each see half a cycle cannot cleanly separate a large phase change at low contrast from a small one at high contrast. The fit needs limits on a_j and b_j to stay physical.
- The ring-minima method avoids both problems, because it uses all the rings and fits the centre and P for every photo.

---

## 11 The slope gives c = (2.8 ± 0.2) × 10⁸ m/s (0:45)

**Say:** "The prediction for the slope contains the speed of light. If we trust N, R and λ, we can turn the measured slope around and solve for c."

- **Formula.** a = 8π²NR²/(λc) [Manual Eq. 5 with A = πR², differentiated with respect to ω] gives c = 8π²NR²/(λa). With N = 13, R = 0.570 m, λ = 640 nm and a = 1.839 s, c = 2.83 × 10⁸ m/s. Note that rad/(rad/s) is a time, so a is in seconds.
- **Comparison.** c/c₀ = 0.94 ± 0.08, 0.7σ from the defined value c₀ = 299 792 458 m/s, two-sided p = 0.49 [`sagnac.txt`].
- **Uncertainty budget.** Relative uncertainties add in quadrature: (u_c/c)² = (u_a/a)² + (u_N/N)² + (2u_R/R)². R enters squared, so its 0.9 % uncertainty counts twice (1.8 %). The turn count, 13 ± 1, alone gives 7.7 %. Everything from our own measurement (fit, calibration, method) is 3.3 %. The total is 8.5 %. λ was given without an uncertainty. Each nanometre would change c by 0.16 % [`sagnac.txt`].
- **Turn count.** Read the other way, if c, R and λ are taken as known, the measured slope fixes N = 13.8 ± 0.3 turns, three times more precise than the given 13 ± 1.
- This is the same comparison as slide 8 (c/c₀ = a_th/a), shown as a speed of light. It is not an independent second test.

---

## 12 A steady phase drift would mimic part of the slope (0:45, optional)

**Say:** "One systematic effect we cannot rule out. We photographed the speeds in increasing order, with rest photos only at the start. So ω and time increase together. A phase that drifts steadily in time would look exactly like extra slope."

- **Top.** ω against time during the run: a staircase that only goes up.
- **Bottom.** Each photo's phase minus the fitted line, against time, with a straight line through each setting's photos. Within single settings the phase drifts by −6.9 to +8.2 mrad/s, and 6 of the 10 drifts are significant (p < 0.05). Between settings the means wander by about σ_x = 0.04 rad. This is the extra scatter of slide 8.
- **Joint fit.** Δφ = aω + Dt + b, with t the time of each setting [Analysis RESULTS.md, `sagnac.txt`]. It gives a = 2.15 ± 0.20 rad/(rad/s) and D = −0.9 ± 0.5 mrad/s, χ²/dof = 5.5/7, p = 0.60. The correlation between a and D is −0.99. The data cannot tell a slope in ω from a drift in time, so we report the slope without D (1.84) and state that a steady drift is not covered by the quoted uncertainty.
- **Likely source.** SMF-28e+ is not single-mode at 640 nm. Its V number there is about 5.6, above the single-mode limit 2.405, so it guides several modes [Manual V-number section, Prep Q12]. Coupling between these modes changes with fibre stress and temperature, which can add a slowly varying phase that is different for the two directions.
- **How to separate it next time.** Take rest photos between the speeds, and reverse the rotation. The Sagnac phase changes sign with the direction of rotation, a drift does not. Sagnac himself measured this way, rotating in one direction and then the other [Pascoli Sec. 2].

---

## 13 The Sagnac phase is linear in ω and matches theory (0:30)

**Say the three numbers:**
1. The slope of Δφ against ω is 1.84 ± 0.06 rad/(rad/s) with the method spread included. The prediction 8π²NR²/(λc) = 1.74 ± 0.14 is 0.7σ away.
2. The largest shift is 3.29 ± 0.06 rad, half a fringe, at ω = 1.81 rad/s. It comes from a travel-time difference of only about 1.1 fs.
3. Read as a speed of light, c = (2.8 ± 0.2) × 10⁸ m/s, 0.7σ from the defined value, limited mainly by the turn count.

**Next steps:** interleave rest photos and reverse the rotation to separate a drift. Count the fibre turns. Use a fibre that is single-mode at 640 nm (or filter the output to one mode).

---

## 14 Thank you

Open for questions. The backup slides cover the most likely ones.

---

## Backup slides and likely questions

**B1 A quadratic calibration would raise the slope by 1.7 %.** "Why a straight line when the residuals curve?" The manual asks for the linear relation [Manual minimum task 1]. A quadratic describes the points (χ²/dof = 6.0/6, by construction because ε is set from it) and changes the slope to 1.871 and c to 2.785 × 10⁸ m/s. Both shifts are 1.7 %, inside the 3.3 % measurement uncertainty. A cubic adds nothing (χ²/dof = 2.6/5, p = 0.75) [`calibration.txt`, `sagnac_update/RESULTS.md`].

**B2 Every error source.** "Which uncertainties are included?" The table is manual task 5. Statistical parts enter the fit. Calibration, method, N and R are added in quadrature afterwards. Not included: the wavelength (not given), a steady drift (not separable) and the calibration model choice (reported separately). Not checked: the sign of the shift, because the rotation direction was not recorded [Analysis RESULTS.md].

**B3 Measured values at each setting.** The numbers behind slide 8, in the notation explained above (0.235(1) = 0.235 ± 0.001). Δφ uncertainties are the per-setting fit weights including σ_x.

**Other questions.**
- *Does it matter that the loop is centred on the axis?* No. Only the enclosed area enters [Schwartz Eq. 4, Prep Q8]. The off-centre loop was not measured.
- *Does the refractive index of the fibre matter?* No [Schwartz Sec. 2, Prep Q10].
- *Why not measure the laser frequency instead of using λ?* The phase fixes the product λc. With a known frequency ν instead of λ, c = √(8π²NR²ν/a), with half the relative uncertainty [Analysis README].
- *Is the effect relativistic?* The classical derivation gives the same first-order result. Sagnac took it as evidence for an ether, but Langevin showed it is fully consistent with relativity [Pascoli, Tartaglia Sec. II].
