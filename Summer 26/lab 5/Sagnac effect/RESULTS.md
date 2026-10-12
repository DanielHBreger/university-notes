# Sagnac effect: results (11 October 2026 data)

Produced by the scripts in `README.md`; full tables in `calibration.txt` and `sagnac.txt`.
Uncertainties are standard uncertainties; their kind is given each time. In "report form" the
uncertainty is rounded to one significant figure and the value to the same digit. Comparisons
use unrounded values.

## Task 1: rotation calibration (`figures/calibration`)

| quantity | value |
|---|---|
| stopwatch (splits) | u(T) = SD/sqrt(n(n+1)); omega to 0.04-0.18 % per setting; one press 0.025-0.10 s |
| split model check | lap lag-1 correlation -0.34; splits predict -0.43 +- 0.09 (simulated); independent laps give this or lower with p = 0.05 |
| speed reproducibility | eps = 0.48 % of omega (scatter about the quadratic, chi2 = dof) |
| linear fit (manual's model), with eps | k = d omega/dS = 0.02039 +- 0.00004 rad s^-1 per unit S (fit), omega0 = -0.0212 +- 0.0006 rad/s; chi2 = 47.2, dof = 7, p = 5e-8 |
| linear, stopwatch only | chi2 = 807, dof = 7, p = 6e-170 |
| quadratic (used) | stopwatch only: chi2 = 154.9, dof = 6, p = 7e-31; with eps: chi2 = dof by construction |
| cubic, with eps (check) | chi2 = 2.6, dof = 5, p = 0.75: no need for a cubic term |
| linear-model error | up to 0.041 rad/s (2.3 % of omega) at S = 90 |

The linear relation the manual asks for is curved: the table turns slower than linear at high
settings. The stopwatch is far more precise than the table is reproducible, so the 0.48 %
scatter between settings dominates. omega at the photographed settings (quadratic; curve,
correlated, then reproducibility, independent; rad/s):
0.2383(6)(12), 0.4073(11)(20), 0.6201(15)(30), 0.8825(18)(43), 1.1173(21)(54), 1.2786(25)(62),
1.4671(35)(71), 1.6319(47)(79), 1.7816(61)(86). Report form, combined: 0.238(1), 0.407(2),
0.620(3), 0.883(5), 1.117(6), 1.279(7), 1.467(8), 1.632(9), 1.78(1).

## Task 2: patterns (`figures/patterns`, `figures/profiles`)

`patterns` shows one photograph per setting, the frame closest to its setting's median phase,
with the fitted dark rings and centre. The centre goes from bright (at rest) to nearly dark (top
speed) as the rings contract. `profiles` shows the sector-averaged intensity against r^2: the
fringes are equally spaced in r^2 (period 107.4k px^2), and the whole pattern moves inwards with omega.

## Task 3: phase shift against omega (`figures/phase_vs_omega`)

| quantity | value | report form |
|---|---|---|
| largest shift | 3.285 +- 0.069 rad (fit errors of both settings) = 0.52 +- 0.01 fringe, at omega = 1.782 +- 0.011 rad/s | 3.29(7) rad |
| slope a | 1.871 +- 0.027 (fit) +- 0.006 (calibration curve) = 1.871 +- 0.028 rad/(rad/s) | 1.87(3) |
| theory a_th = 8 pi^2 N R^2/(lambda c) | 1.738 +- 0.137 (N: 0.134, R: 0.031; lambda not given) | 1.7(1) |
| a / a_th | 1.077 +- 0.086 | 1.08(9) |
| a - a_th | +0.13 +- 0.14 rad/(rad/s) = 0.95 combined sigma, p = 0.34 | |
| N implied by a (given R, lambda) | 14.00 +- 0.32 turns | 14.0(3) |

Second method (prep Q5), the cosine fit of each frame's profile (shown as open squares in
`figures/phase_vs_omega`):

| quantity | value | report form |
|---|---|---|
| slope a_cos | 1.910 +- 0.025 (fit) +- 0.007 (calibration curve) = 1.910 +- 0.025 rad/(rad/s) | 1.91(3) |
| a_cos / a_th | 1.099 +- 0.088; 1.2 combined sigma from theory, p = 0.22 | 1.10(9) |
| largest shift | 3.395 +- 0.062 rad (fit errors of both settings) = 0.54 +- 0.01 fringe | 3.39(6) rad |
| between-setting scatter | sigma_x = 0.041 rad (minima 0.045); line without it: chi2 = 83, dof = 8, p = 1e-14 | |
| fit rms, median per setting | 9.5 counts at rest, 5.5-6.7 at the three lowest speeds, 1.4-2.3 above | |
| cosine minus minima, frame by frame | slope +0.038 +- 0.026 rad/(rad/s), 1.5 sigma | |

Both methods agree with Eq. (5). They differ in slope by 1.5 sigma, and the per-setting
differences show where: relative to rest, every moving setting reads +0.10 to +0.19 rad higher
with the cosine fit. So the two methods disagree mainly about the rest frames. That is where the
cosine model fits worst: the bright, flat-topped centre at rest is not a cosine (rms 9.5 counts).
The difference is a method systematic of about 0.04 in the slope, larger than the cosine fit's
statistical advantage. On precision, the prep's expectation (Q5) that the fit uses all the points
and is therefore more precise does not show here: the per-setting SEMs are alike (0.003-0.035
against 0.004-0.028 rad), and both slope errors are set by the slow phase wander (0.025 against
0.027). The minima method stays the main result, because it needs no model of the envelope or of
the camera response.

### Phase-extraction methods compared (`figures/methods_comparison`, `figures/methods_illustration`)

Three more methods (`methods.py`; scans in `methods.txt`):

| method | what it reads | slope a, rad/(rad/s) | a / a_th | minus minima, frame by frame | robustness scan |
|---|---|---|---|---|---|
| minima (main) | phase at the ring centre | 1.871 +- 0.028 | 1.08(9) | | |
| cosine fit | phase at the ring centre | 1.910 +- 0.025 | 1.10(9) | +0.038 +- 0.026 (1.5 sigma) | fit range: 1.896 at s < 3.0e5 (one-off check, not in the pipeline) |
| template matching | mean shift over s = 0.2-1.8e5 px^2 | 1.798 +- 0.022 | 1.03(8) | -0.073 +- 0.034 (2.2 sigma) | 1.794-1.801 over three ranges; 1.864 +- 0.041 without the centre |
| 2D Fourier sideband on (theta, r^2) | fringe-strength-weighted mean over the map | 1.896 +- 0.039 | 1.09(9) | +0.025 +- 0.038 (0.7 sigma) | 1.895-1.897 over sideband widths 0.15-0.40 x carrier |
| AIA (plain) | weighted mean over the window | 1.609 +- 0.027 | 0.93(7) | -0.26 +- 0.04 (6.8 sigma) | 1.61-1.88 over windows; 2.05-2.11 with a brightness factor |

(Slope errors are fit and calibration curve combined. sigma_x: minima 0.045, cosine 0.041,
template 0.034, Fourier 0.070, AIA 0.047 rad. Report form: 1.87(3), 1.91(3), 1.80(2), 1.90(4),
1.61(3).)

- **Template matching** gives the smallest scatter between settings (sigma_x = 0.034 rad) and is
  stable against its comparison range. Its per-setting differences from the minima follow the
  fringe period of each setting (correlation -0.92, p < 0.001, n = 9). The fitted slope,
  -0.050 +- 0.008 rad per 1k px^2 of period change, matches the -0.049 expected for a phase read
  at r^2 = 0.9e5 px^2, the middle of the matching range. So template matching measures the phase
  mid-window, where a change in wavefront curvature adds to it. The minima and the cosine fit
  measure it at the ring centre, where a curvature change adds nothing, so the centre phase is
  the right observable for a Sagnac phase that is the same across the beam. The template slope
  sits below because the three settings with the largest period (44.07, 64.00, 73.65) happen to
  be at higher speeds.
- **2D Fourier sideband,** adapted from a friend's `sagnac_2d.py` by first resampling each image on
  (theta, r^2), where the rings are straight. Its slope agrees with the minima (0.7 sigma frame by
  frame) and does not move with the sideband width. Its per-setting differences from the minima
  also follow the fringe period, but less than template matching's (correlation -0.74, p = 0.02;
  -0.038 rad per 1k px^2). That fits a phase averaged mostly over the bright inner fringes,
  around r^2 = 0.7e5 px^2. It has the largest scatter between settings (sigma_x = 0.070 rad), so
  it is the least precise. Panel (f) of `methods_illustration` shows that the local phase shift
  is uniform over the map to within about 0.05 rad.
- **AIA** works on synthetic profiles: it recovers a 3.3 rad sweep to 0.01 rad
  (`tests/test_methods.py`). On the photographs its slope moves from 1.61 to 1.88 with the
  window, and to 2.05-2.11 with a per-frame brightness factor. Its rms residual is 5.1 counts
  (full window), against 1.4-2.3 counts for the single-frame cosine fits above S = 44. Panel (c) of
  `methods_illustration` shows why: the r^2 bins
  do not follow single sinusoids of the frame phase, most clearly at rest. AIA assumes that only
  the phase changes between frames. Here the fringe period (1.7 %), the brightness (template gain
  0.92-1.04) and the bright flat centre all change too, and a 1.7 % period jitter alone costs
  about 0.2 rad per frame on synthetic data. AIA is shown for completeness but is not used.
  Its bin phases psi(s) still come out nearly linear in r^2 (P = 105.7k px^2; deviations up to
  0.56 rad), which is a model-free sign that the rings are equally spaced in r^2.
- **Method systematic.** The four usable methods give 1.798-1.910. Half that range, 0.056,
  covers the choice of phase extraction:
  a = 1.871 +- 0.028 (fit and calibration) +- 0.056 (method) rad/(rad/s). With this term,
  a - a_th = 0.88 combined sigma, p = 0.38. Report form: a = 1.87(6) rad/(rad/s).

The measurement agrees with Eq. (5). The theory uncertainty is dominated by N = 13 +- 1
(7.7 %). The measured slope is five times more precise than the prediction; read the other way,
it fixes the number of turns to 14.0 +- 0.3 if R and lambda are right.

Fit details: with the frame SEMs (0.004-0.028 rad) and the speed reproducibility (0.002-0.016
rad in phase) alone, the line gives chi2 = 62, dof = 8, p = 2e-10 (frame SEMs alone: chi2 = 183).
The settings scatter about the line by more than either explains. sigma_x = 0.045 rad per setting
is added in quadrature (fixed by chi2 = dof). With the frame SEMs alone the slope is
1.875 +- 0.007 (fit): about the same value, with an error that is too small. A quadratic term is
0.07 +- 0.06 rad/(rad/s)^2, 1.3 sigma from zero: no sign of nonlinearity.

Checks (not needed for the report):
- Eq. (7) ring by ring with one global period instead of the centre phase: slopes 1.864 (inner
  ring) and 1.836 (second ring) rad/(rad/s). The two rings differ by up to 0.14 rad at single
  settings, because the fringe period varies over 105.8-109.4k px^2 between settings (1.7 %);
  the centre phase accounts for this. By the beam picture in `README.md` (prep Q3, Q7) a P change
  is a change of the wavefront curvature difference, e.g. a lens-fibre distance moving as the
  table turns.
- Dropping the first two frames at S = 12.59 (taken 26 s before the rest, ring-fit rms 8 px
  against 3-6): slope 1.876, a change of 0.004 (0.2 sigma).

## Task 4: small loop off the axis

No photographs yet. The code is ready (`README.md`). Prediction from Eq. (5) with N = 46 +- 1,
R = 15.5 +- 0.1 cm: a_th = 0.455 +- 0.011 rad/(rad/s) (N: 0.010, R: 0.006), i.e. 0.81 rad =
0.13 fringe at omega = 1.78 rad/s, a quarter of the large loop's shift. If the phase wanders as
much as on 11 October (sigma_x = 0.045 rad per setting, same settings), the slope would be
measured to about 0.027 (6 % of a_th), and a - a_th would carry about 0.029 with the theory's
0.011. The small loop's test is therefore limited by the phase wander, not by N or R.

## Task 5: errors

| source | kind | size | in the result |
|---|---|---|---|
| frame-to-frame phase scatter | statistical | SD 0.016-0.11 rad per frame, SEM 0.004-0.028 rad per setting | yes (fit) |
| speed reproducibility at a setting | statistical, from the calibration | 0.48 % of omega; 0.002-0.016 rad in phase | yes (fit, as error in omega) |
| slow phase wander between settings | statistical, from the fit | sigma_x = 0.045 rad per setting | yes (fit) |
| calibration curve | calibration, correlated across settings | 0.006 on a (0.3 %) | yes, separately |
| phase-extraction method | systematic, half the range of minima, cosine fit, template matching, 2D Fourier | 0.056 on a (3 %) | yes, separately |
| N = 13 +- 1, R = 57.0 +- 0.5 cm | given | 7.7 % and 1.8 % on a_th | yes (theory) |
| wavelength | not given | 0.16 % on a_th per nm | **no**, stated |
| drift linear in time | systematic, not separable | see below | **no** |

The settings were photographed in increasing order with the only rest frames at the start. A
phase drift linear in time therefore has the same signature as the Sagnac slope. Fitting both
gives a = 2.27 +- 0.23 rad/(rad/s) and D = -1.1 +- 0.6 mrad/s with correlation -0.993, so the data
cannot separate them. Within each setting the phase drifts by -6.9 to +8.2 mrad/s, changing sign
between settings. The residuals about the line are only 0.045 rad over 12 minutes, so the phase
wanders without random-walking. That does not exclude a steady drift. Rest frames at the end, or
reversing the rotation, would.

A likely cause of the wander (an interpretation, not tested): SMF-28e+ carries a few modes at
640 nm (V = 5.6, prep Q12), and the camera records the whole output with no single-mode filter at
the common port. A fibre gyroscope needs such a filter for the two directions to be exactly
reciprocal, so mode coupling that changes with fibre stress and temperature leaves a slowly
varying non-reciprocal phase. The speckle in the photographs fits the same picture.

Also: omega is the calibration at the same setting, not measured during the photographs (eps
covers the scatter seen in the calibration). The sign of the shift is not compared because the
rotation direction is not recorded.
