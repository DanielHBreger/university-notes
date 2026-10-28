# Sagnac effect, lab 5

Analysis of the 11 October 2026 measurements against the Sagnac manual (minimum tasks 1-3 and 5) and
Schwartz, Eur. J. Phys. 38, 015301 (2017). Numbers and their meaning are in `RESULTS.md` (short
version `SUMMARY.md`); the report material is `report_figures/`, `REPORT_figures_tables.tex`,
`REPORT_TODO.md` and `PLAN_CHECKLIST.md`.

## Data

- `Rotation calibration.csv`: times of consecutive single turns of the table: 5 turns at each of
  the five lowest speed settings, 10 at the four highest (columns `turn 1` ... `turn 10`, empty cells
  skipped). `Speed` is the motor controller setting S.
- `Measurements/<S>/<YYYYMMDDhhmmss>.png`: photographs of the fringe pattern with the fibre loop
  centred on the axis, 10-22 per setting, S = 0 (at rest), 12.59, 20.74, 31.11, 44.07, 55.83, 64.00,
  73.65, 82.18, 90. 1920x1080 RGBA. Only the red channel carries light (max 239 counts, nothing
  saturated). The settings were taken in increasing order between 16:21 and 16:33; the time in the
  file name is used.
- Fibre loop, centred on the rotation axis (constants in `sagnac.py`): N = 13 +- 1 turns, radius
  R = 57.0 +- 0.5 cm. lambda = 640 nm, no uncertainty given.

## Pipeline

Run from this folder with the lab's Python (`..\chaotic circuit\.venv\Scripts\python.exe`, or any
Python 3.11+ with numpy, scipy, pandas, matplotlib and Pillow; LaTeX is used for figure text
when it is installed):

```
python calibration.py     # task 1 -> calibration.json, calibration.txt
python fringes.py         # dark rings, minima and cosine-fit phases of the 161 photographs -> fringes_frames.csv, fringes_profiles.npz (a few minutes)
python methods.py         # template matching -> methods_frames.csv (about a minute)
python pixels.py          # fixed pixels, 1000 random draws -> pixels_frames.csv, pixels.json (a few minutes)
python sagnac.py          # tasks 3 and 5, the speed of light -> sagnac_speeds.csv, sagnac.txt, sagnac.json
python report_figures.py  # report_figures/*.pdf, *.png (about a minute)
python -m pytest -q -p no:cacheprovider tests   # synthetic checks of template matching and the fixed-pixel fit, report-form rounding, Eq. (5) for c
```

`fringes.py`'s starting centre, fringe-period guess and angular sector (`START_CENTRE`, `P_GUESS`,
`SECTOR_DEG`) are set for the 11 October camera view. Figures follow the lab report style
(`report_style.py`, copied unchanged from `chaotic circuit/`): column width unless noted, vector PDF
plus 600 dpi PNG, no titles, no fit statistics in the image, no residual panels.

## Method

**Calibration** (`calibration.py`). Per setting, T is the mean of n = 5 or 10 consecutive timed
turns and omega = 2 pi/T. u(T) combines the lap scatter SD/sqrt(n) with the reaction time of the
presses: the laps are consecutive, so only the first and last press set the mean, and with
sigma_r = 0.25 s per press that part is sqrt(2) sigma_r/n. Together 0.15-1.1 % of omega. The line
omega = k S + omega0 is fitted with scipy's curve_fit, weighted by the timing uncertainties. The
points lie farther from the line than those allow (chi2/dof = 2.9, p = 0.005): the table turns slower than
linear at high settings, 3.0 % below the line at S = 90. The photos use the line, as task 1 asks.
A photo setting was not timed, so its omega carries the line's prediction uncertainty: the line's
own (0.2-0.3 %) and the 1.3 % rms scatter of the timed settings about it, 1.3 % in all, common to
every setting; that relative size is the calibration part of the Sagnac slope's uncertainty. The
3.0 % bend at the top setting is stated as a limitation.

**Fringe phase** (`fringes.py`). A Gaussian beam carries the phase k r^2/(2R(z)). The two output
beams have different wavefront radii R1, R2 at the camera, so the fringes are rings with phase
2 pi r^2/P + phi_c, where P = 2 lambda/|1/R1 - 1/R2|. They are equally spaced in s = r^2, and the
Sagnac phase shifts phi_c. A relative tilt only moves the ring centre; a change of a lens-fibre
distance changes P, not phi_c.

In each frame, the minimum of every dark ring is found along rays over the sector with fringes
(-75 to 180 deg; the lower left has none). One centre and the ring radii come from a
concentric-circle fit to these minima. s_m = R_m^2 regressed on ring order m gives P and phi_c,
the phase at the ring centre. This is the manual's Eq. (7) (minimum shift over local fringe
period, averaged over minima) done in r^2, where the period is the same for every ring, and the
r^2 version of Schwartz's Eq. (5). Because Delta phi = 2 pi Delta s / P is a ratio, no
pixel-to-length calibration is needed, and minimum positions do not depend on the camera's
intensity response.

Second method (prep Q5, the cosine fit): the intensity averaged over the same sector, in bins of
s about the same centre (s < 2.6e5 px^2, the three dark rings), is fitted with
C e^(-s/sC) + B e^(-s/sB) cos(2 pi s/P + phi) + D. The two exponentials are the envelopes of two
Gaussian beams of different size, and D is a floor. phi is the same phase at the ring centre. This
method uses every bin, but it relies on the envelope model and on a linear camera response; the
minima rely on neither.

Third method (`methods.py`, template matching): each r^2 profile is matched to the averaged rest
profile, shifted by D and scaled (gain and offset free). The phase is 2 pi D / P_T, where P_T is the
template's own period. It uses no fringe model but reads an average phase over its window, so a
change of the fringe period between frames enters it.

Fourth method (`pixels.py`, fixed pixels; prep Q5, asked for by the instructor): ten pixels at
fixed image positions are read in every frame and fitted with I_kj = a_j + b_j cos(psi_j + delta_k),
one phase delta_k per frame, under two physical limits (b_j >= 0; a_j - b_j >= the camera's dark
level), which the fit needs because no pixel sweeps a whole fringe over this run. A fixed pixel
reads the local phase at its own radius, so a moving centre or spacing reads as phase, and the
slope depends on which pixels are drawn: over 1000 random draws of ten it spreads by 0.19 (11 %).
It is shown beside the other readouts with that spread as its uncertainty and is not in the
method range.

**Phase against omega** (`sagnac.py`). Per setting: mean of the frame phases; its uncertainty in the
report is s, the scatter of the settings about the line (the SEM is the photo noise only).
A straight line Delta phi = a omega + b through the ten setting means (scipy's linregress; Schwartz
Eq. 7): the slope's uncertainty comes from how far the points lie from the line (s = 0.044 rad
rms). Measured against the photos' SEM, the line has chi2/dof = 205/8: the settings differ five
times more than the photo noise allows, because the phase wanders slowly, so the regression's own
scatter is the right uncertainty. The slope's uncertainty has three parts added in quadrature: the
regression (1.3 %), the phase-reading method (half the range of the three methods, 2.9 %) and the
calibration (the line's prediction uncertainty, 1.3 %).
Theory: a_th = 8 pi^2 N R^2/(lambda c), Eq. (5), independent of where the axis lies (Schwartz;
prep Q8) and of the fibre's index (prep Q10).

**Speed of light** (`sagnac.py`). Eq. (5) solved for c with the measured slope:
c = 8 pi^2 N R^2/(lambda a). The slope's parts carry over with the same relative size (c ~ 1/a);
N and R enter as given; lambda has no uncertainty given and is left out. Since c/c_def = a_th/a,
this is the theory comparison read the other way.

## Conventions

- Every number carries its uncertainty and its kind: statistical (regression), systematic
  (method, calibration) or given (N, R). The wavelength uncertainty was not given; it is reported
  as missing, not set to zero.
- Report form (`sagnac.round1`): the uncertainty to one significant figure and the value to the
  same digit, half up on the decimal digits. Comparisons are computed from unrounded values.
- Not covered by the data: the rotation direction (sign), a steady drift in time.

## Earlier stages, kept for the record

The analysis was presented on 19 October (`presentation_slides/`, PNG slides and presenter notes)
and revised for the report on 22 October (this pipeline: simpler uncertainty treatment, same
numbers). The presentation's own code (`presentation_slides/source/`, which reads
`sagnac_update/`) expects the calibration file of that stage; to rerun it, copy
`previous_reports/calibration_eps_2026-10-13.json` to `calibration.json` first. The folders
`figures/` (quadratic-calibration figures of 13 October) and `sagnac_update/` (the linear-calibration
stage of 18 October), the earlier results documents in `previous_reports/` and
`presentation_slides/source/numbers.json` are superseded by `sagnac.txt` and `RESULTS.md`.
`docs/` explains each script.
