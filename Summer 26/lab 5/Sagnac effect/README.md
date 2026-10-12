# Sagnac effect, lab 5

Analysis of the 11 October 2026 measurements against the Sagnac manual (minimum tasks 1-5) and
Schwartz, Eur. J. Phys. 38, 015301 (2017). Numbers and their meaning are in `RESULTS.md`.

## Data

- `Rotation calibration.csv`: times of consecutive single turns of the table, taken as
  stopwatch splits: 5 turns at each of the five lowest speed settings, 8 at the four highest
  (columns `turn 1` ... `turn 10`, empty cells skipped). `Speed` is the motor controller setting S.
- `Measurements/<S>/<YYYYMMDDhhmmss>.png`: photographs of the fringe pattern with the large loop
  centred on the axis, 10-22 per setting, S = 0 (at rest), 12.59, 20.74, 31.11, 44.07, 55.83, 64.00,
  73.65, 82.18, 90. 1920x1080 RGBA. Only the red channel carries light (max 239 counts, nothing
  saturated). The settings were taken in increasing order between 16:21 and 16:33; the time in the
  file name is used.
- `_transfer/meas.tgz`: an archive of the same photographs and the calibration CSV (not read by
  the code).
- Fibre loops (constants in `sagnac.py`, `LOOPS`): large loop N = 13 +- 1 turns, radius
  R = 57.0 +- 0.5 cm; small loop (task 4, off the axis) N = 46 +- 1 turns, R = 15.5 +- 0.1 cm.
  lambda = 640 nm, no uncertainty given.

## Pipeline

Run from this folder with the lab's Python (`..\chaotic circuit\.venv\Scripts\python.exe`, or any
Python 3.11+ with numpy, scipy, pandas, matplotlib and Pillow; LaTeX is used for figure text
when it is installed):

```
python calibration.py   # task 1 -> calibration.json, calibration.txt, figures/calibration
python fringes.py       # dark rings, minima and cosine-fit phases of the 161 photographs -> fringes_frames.csv, fringes_profiles.npz (a few minutes)
python methods.py       # template matching, AIA, 2D Fourier -> methods_frames.csv, methods.txt/.json, figures/methods_illustration (figure*) (a few minutes)
python sagnac.py        # task 3, all methods -> sagnac_speeds.csv, sagnac.txt, sagnac.json, figures/phase_vs_omega, figures/methods_comparison
python patterns.py      # task 2 -> figures/patterns (figure*), figures/profiles
python -m pytest -q -p no:cacheprovider tests   # synthetic checks of the three methods, report-form rounding
```

Task 4 (small loop off the axis), once its photographs are in a folder laid out like
`Measurements/` (one subfolder per setting, the same calibration):

```
python fringes.py --data "<folder>" --out fringes_frames_small.csv
python sagnac.py --loop small    # -> sagnac_small.txt/.json, sagnac_speeds_small.csv, figures/phase_vs_omega_small (methods.py covers the large loop only)
python patterns.py --loop small  # -> figures/patterns_small, figures/profiles_small
```

`fringes.py`'s starting centre, fringe-period guess and angular sector (`START`, `P_GUESS`,
`THETA_DEG`) are set for the 11 October camera view. Check them, and the ring overlays in
`patterns_small`, on new photographs.

Figures follow the lab report style (`report_style.py`, copied unchanged from `chaotic circuit/`):
column width unless noted, vector PDF plus 600 dpi PNG, no titles, no fit statistics in the image.

## Method

**Calibration** (`calibration.py`). The turns are stopwatch splits, so neighbouring laps share
an endpoint. The mean lap is (t_n - t_0)/n with stopwatch uncertainty SD/sqrt(n(n+1)). The laps'
lag-1 correlation (-0.34) agrees with this model (simulated -0.43 +- 0.09); independent laps
would give -0.16. omega = 2 pi/T. The settings scatter about any smooth curve by more than the
stopwatch allows. This is modelled as a reproducibility eps of the table speed at a fixed setting
(a fraction of omega, added in quadrature), with eps set by chi2 = dof for the quadratic fit. The
manual's linear model is reported but rejected even with eps. The quadratic converts the
measurement settings to omega with two uncertainty parts: the curve's covariance (correlated
across settings) and eps omega (independent per setting).

**Fringe phase** (`fringes.py`). A Gaussian beam carries the phase k r^2/(2R(z)). The two output
beams have different wavefront radii R1, R2 at the camera, so the fringes are rings with phase
2 pi r^2/P + phi_c, where P = 2 lambda/|1/R1 - 1/R2|. They are equally spaced in s = r^2, and the
Sagnac phase shifts phi_c. A relative tilt only moves the ring centre; a change of a lens-fibre
distance changes P, not phi_c.

In each frame, the minimum of every dark ring is found along rays over the sector with fringes
(-75 to 180 deg; the lower left has none). One centre and the ring radii come from a
concentric-circle fit to these minima. s_m = R_m^2 regressed on ring order m gives P and phi_c,
the phase at the ring centre. This is the manual's Eq. (7) (minimum shift over local fringe
period, averaged over minima) done in r^2, where the period is the same for every ring. It is
the r^2 version of Schwartz's Eq. (5). Because Delta phi = 2 pi Delta s / P is a ratio, no
pixel-to-length calibration is needed, and minimum positions do not depend on the camera's
intensity response.

Second method (prep Q5, the cosine fit): the intensity averaged over the same sector, in bins of
s about the same centre (s < 2.6e5 px^2, the three dark rings), is fitted with
C e^(-s/sC) + B e^(-s/sB) cos(2 pi s/P + phi) + D. The two exponentials are the envelopes of two
Gaussian beams of different size, and D is the dark floor. phi is the same phase at the ring
centre. This method uses every bin, but it relies on the envelope model and on a linear camera
response; the minima rely on neither.

Three more methods (`methods.py`):
- **Template matching:** each r^2 profile is matched to the averaged rest profile, shifted by D
  and scaled (gain and offset free). The phase is 2 pi D / P_T, where P_T is the template's own
  period. It uses no fringe model.
- **Phase shifting with unknown steps (AIA, Wang & Han 2004):** every r^2 bin brightens and dims
  sinusoidally with the frame's phase. Alternating least squares finds each bin's sinusoid and
  each frame's phase together, so no period, envelope or cosine shape in r is assumed. It
  assumes that only the phase changes between frames.
- **2D Fourier sideband (Takeda),** adapted from a friend's `sagnac_2d.py`. A single carrier
  does not fit rings, so each image is first resampled on (theta, r^2) about its ring centre,
  where the rings become straight fringes. Then the slow background is removed, the positive
  sideband is kept with a Gaussian filter, and the phase relative to rest is
  arg sum Z_k conj(Z_rest) over the map.

These three read an average phase over their window, not the phase at the ring centre, so a
change of the fringe period between frames enters them. Every method goes through the same
averaging and fit in `sagnac.py` and is compared with the minima frame by frame, since all of
them read the same photographs. The robustness scans (other ranges and windows, AIA with a
per-frame brightness factor, other sideband widths) are in `methods.txt`.
`tests/test_methods.py` checks the three new methods on synthetic data.

**Phase against omega** (`sagnac.py`). Per setting: mean of the frame phases, u = SD/sqrt(n).
Weighted line phi = a omega + b with per-setting error sqrt(SEM^2 + (a eps omega)^2 + sigma_x^2).
The speed reproducibility enters as an error in omega (effective variance, iterated on a). With
SEM and reproducibility alone the line is rejected (chi2 = 62, dof = 8, p = 2e-10): the phase
wanders on time scales longer than one setting's photographs. sigma_x = 0.045 rad per setting is
set by chi2 = dof. The slope's fit and calibration-curve uncertainties are given separately; the
curve covariance is propagated through the fit. Theory: a_th = 8 pi^2 N R^2/(lambda c), Eq. (5),
independent of where the axis lies (Schwartz; prep Q8) and of the fibre's index (prep Q10).

## Conventions

- Every number carries its uncertainty and its kind: fit (statistical), SEM of frames,
  calibration curve, speed reproducibility, or given (N, R). The wavelength uncertainty was not
  given; it is reported as missing, not set to zero.
- Every chi2 comes with dof and p.
- `*.txt` keep more digits. The "report form" lines (`sagnac.round1`) round the uncertainty to
  one significant figure and the value to the same digit, half up on the decimal digits (2.25
  becomes 2.3, not the 2.2 of float round()). Comparisons are computed from unrounded values.
- Not covered by the data so far: task 4 (no photographs of the small loop yet) and the rotation
  direction (sign).
