# Sagnac presenter notes

## 01 — Measuring rotation with the Sagnac effect

We use interference to measure the optical phase caused by rotation. The experiment has two measured quantities: the angular velocity of the table and the displacement of the camera's interference rings. We obtain them independently and then test the relation between them.

Transition: first show how rotation enters the optical phase.

## 02 — The Sagnac interferometer

Trace the light from the laser to the beam splitter, through the large fibre loop in opposite directions, and back to the camera. During a rotation, the two beams have different travel times. The first equation is the Sagnac travel-time difference for N turns, each enclosing area A. Multiplying that difference by the light's angular frequency gives the second equation and hence the displayed phase–rotation relation. The sign depends on the chosen direction and phase convention.

The measured configuration is the centred large loop: 13 ± 1 turns, radius 57.0 ± 0.5 cm, laser wavelength 640 nm. The smaller loop also appears in the manual's schematic but is not the configuration analysed here. The apparatus predicts a linear relation. We will evaluate its slope when we have measured both axes.

Transition: the motor's setting is not itself an angular velocity, so the horizontal axis needs a calibration.

Source: the experiment manual, apparatus schematic on page 12; setup parameters in the same caption. Original schematic: `presentation_source/assets/optical_setup_manual.png`. The equations agree with the Sagnac relation used in `sagnac.py`.

## 03 — Rotation calibration

Time several complete rotations at each motor setting. The inset shows the actual five turn times at S = 45.09 and their mean, 6.92 s. Convert the mean period to angular velocity with 2π/⟨T⟩. Repeat at nine settings to obtain the main graph, then fit ω = kS + ω₀.

The saved linear calibration is the conversion used for every subsequent angular velocity and result. k is per unit of the dimensionless motor setting. The intercept applies to the rotating motor data; the stationary photograph series is independently set to ω = 0. The plotted timing errors include the retained 0.48% speed reproducibility term. The optional image compares calibration models.

Transition: with the rotation axis calibrated, inspect what the camera actually sees.

Source: the timing CSV loaded by `calibration.load_turns()`, and `calibration.json`. The plotted calibration is the saved degree-one conversion, without refitting it for this presentation.

## 04 — Rotation shifts the interference rings

The three stills run from rest through 0.88 rad/s to 1.81 rad/s. They are original photographs, using identical image dimensions. Follow a dark ring as the rotation increases: it moves inward. The central region changes from bright to dark.

The visual observation alone does not yet give a phase in radians. The next image measures the ring displacement relative to a fringe period.

Sources: `Measurements/0/20261011162151.png`, `Measurements/44.07/20261011162806.png`, `Measurements/90/20261011163338.png`. Rotation rates are from `sagnac_speeds.csv` using the linear calibration.

## 05 — The ring displacement gives the phase shift

The circles overlay the fitted dark-ring radii on the rest and fastest-rotation photographs. A common fitted centre removes the effect of a simple image translation. Dark rings are equally spaced in squared radius, rather than in radius.

For each frame, fit Rₘ² = Pm + s₀, where P is that frame's fringe period in squared pixels. Dividing by P gives the plotted straight line Rₘ²/P = m + η, with η = s₀/P. The intercept is the fraction of a fringe between the centre and the innermost dark ring. At a dark ring the optical phase is π modulo 2π, so the centre phase is φc = π − 2πη. Its change is therefore Δφ = −2πΔη.

The intercept falls by about 0.52 fringe between the two photographs. That corresponds to approximately 3.3 rad. These are representative frames illustrating the operation. The next graph uses the averages of all 161 frames at the ten settings, not only these two photographs. The fringe period is fitted separately for each frame, so this procedure also allows the ring spacing to change.

Source: the two photographs listed above, with their actual centres and radii from `fringes_frames.csv`. The normalized fit and intercept displacement are recomputed directly from those saved radii. The representative-frame shift is illustrative; the aggregate maximum is 3.29 ± 0.06 rad.

## 06 — Phase shift versus angular velocity

Apply the same extraction to all photographs and plot the mean phase shift at each setting against its calibrated angular velocity. The orange line fits the data with a free intercept. Its slope, including fit and calibration-curve uncertainty, is 1.84 ± 0.02 rad/(rad/s).

Now return to the phase relation on image 02. Its coefficient of angular velocity predicts a slope of 1.738 ± 0.137 rad/(rad/s) from the loop dimensions, turn count and wavelength. The grey line and band show that prediction and its standard uncertainty. The measured trend lies within that band. The largest mean shift is 3.29 ± 0.06 rad, slightly over half a fringe.

The line is fitted to absolute phases with a free intercept; subtracting the measured rest phase only changes the vertical origin for display. Vertical error bars are the per-setting weights of the fit, including the fitted excess-scatter term. They are not independent errors of the plotted phase differences, since every difference shares the rest reference. Horizontal errors combine calibration-curve and reproducibility terms. The quoted endpoint uncertainty separately combines the rest and highest-speed fit errors.

Transition: before using the slope to infer c, check how much the result depends on the phase extraction.

Sources: `sagnac_speeds.csv`, `sagnac.json`, `calibration.json`. The main fit is linear in angular velocity and uses the linear motor calibration.

## 07 — Sensitivity to the phase-extraction method

Apply the alternative methods to the same photographs. The left graph subtracts the ring-minimum phase change from each alternative phase change, so small differences that would be hidden on the previous full-range plot become visible. Every method uses its own rest reference. The right graph summarizes the fitted slopes.

The cosine fit estimates the phase at the centre, as the minima method does. Template matching estimates phase over an image window, so changes in the fringe period can affect its comparison with centre phase. The slope range of the minima, cosine fit and template matching is 1.768–1.877 rad/(rad/s). Use half that accepted range, 0.0547 rad/(rad/s), as the method contribution. Combining it in quadrature with the fit and calibration uncertainty gives the displayed main slope, 1.84 ± 0.06 rad/(rad/s). The shaded range is the spread of the central method estimates, not a confidence interval.

The random-points method follows the brightness of ten pixels at fixed image positions and finds the one phase change per photograph that explains all ten. It is drawn with an open marker and is not part of the accepted-method range: a fixed pixel also responds to changes in the fringe spacing and to small moves of the ring centre, and its slope depends on which ten pixels are drawn (central 68 % of 31 random draws: 1.51–1.90 rad/(rad/s)). It remains documented in the analysis files.

Transition: solve the same Sagnac scale relation for the speed of light and propagate the uncertainty.

Sources: `sagnac_speeds.csv`, `sagnac.json`, `methods.py` and `methods.json`. Error bars on individual slopes include fit and calibration covariance. They do not yet include the shared method-spread contribution.

## 08 — The speed of light from the Sagnac slope

Invert the phase–rotation relation with the measured slope. The result is c = (2.8 ± 0.2) × 10⁸ m/s. The dashed line marks the exact defined value, 299,792,458 m/s.

Read the uncertainty graph from top to bottom. The central estimate stays the same. Each row adds the named independent standard-uncertainty component in quadrature, widening the interval around it. Fit and rotation calibration give about 1.3%. Adding the extraction-method spread gives 3.3%. Adding loop-radius uncertainty gives 3.7%, and turn-count uncertainty brings the reported total to 8.5%. The assumed 13 ± 1 turns is the largest contribution. The full interval includes the defined value. This is the slope comparison expressed as an inferred c, rather than an independent second measurement.

Discuss the limits here once: the quoted budget covers the components drawn in the graph. Wavelength uncertainty was not supplied. Rotation increased with acquisition time, so steady optical drift and the Sagnac slope are difficult to separate in this run. The joint drift fit has correlation −0.992 between their fitted coefficients. The 0.041 rad excess-scatter term in the main phase fit accounts for observed scatter but does not establish the absence of a drift bias. The linear motor calibration has residual curvature; the optional image quantifies that model sensitivity. A future run that interleaves rest measurements and reverses the rotation would help separate rotation from drift. Only the centred large-loop configuration was measured.

Sources: `sagnac.json`, including `speed_of_light` and `drift_check`; the reported independent uncertainty components are propagated by `sagnac.py`. The defined c is exact in SI. No wavelength or unresolved-drift component has been added to the displayed uncertainty.

## Optional A1 — Calibration-model sensitivity

The calibration residuals reveal curvature about the linear baseline. The quadratic check reduces the residual χ² from 47.2/7 to 6.0/6 with the same retained reproducibility model. Its speed-scatter estimate is set by χ² = dof, so that second value is not an independent validation of the noise model.

Rerunning the downstream analysis with the quadratic calibration changes the slope from 1.839 to 1.871 rad/(rad/s), about +1.7%, and c from 2.833 to 2.785 × 10⁸ m/s, about −1.7%. This is an optional sensitivity check. Every main image and the main reported result continue to use the linear calibration. The displayed bottom-row numbers retain extra digits only to expose this small model difference.

Sources: `calibration.json`, `sagnac.json`, `optional_quadratic/sagnac.json`.
