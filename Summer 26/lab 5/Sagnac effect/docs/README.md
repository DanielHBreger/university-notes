# Code explanations

One file per script, as in `chaotic circuit/docs/`. Each file walks the script function by function, in the order the data flows through it, and stops on the lines where the code does something that is not obvious. Read them in this order and the whole pipeline follows.

| Order | File | Script | What it settles |
|---|---|---|---|
| 1 | [calibration.md](calibration.md) | `calibration.py` | ω of the timed settings with the reaction-time uncertainty, the weighted line and its uncertainty |
| 2 | [fringes.md](fringes.md) | `fringes.py` | the ring finder, the phase at the ring center of every photo, the cosine fit |
| 3 | [methods.md](methods.md) | `methods.py` | template matching |
| 4 | [pixels.md](pixels.md) | `pixels.py` | the fixed-pixel readout and the spread of its slope between random draws |
| 5 | [sagnac.md](sagnac.md) | `sagnac.py` | setting means, the regression, the three-part uncertainty, c, the drift check, the bend: every number of `RESULTS.md` |
| 6 | [report_figures.md](report_figures.md) | `report_figures.py` | the report's figures |
| 7 | [tests.md](tests.md) | `tests/` | what each test pins down |
| — | [final_numbers.md](final_numbers.md), [random_points_scan.md](random_points_scan.md), [uncertainty.md](uncertainty.md), [patterns.md](patterns.md) | presentation-stage and earlier scripts | kept for the record; each starts with a note on its status |

## Which code made which numbers

The report's numbers are made by the top-level pipeline, `calibration.py` → `fringes.py` → `methods.py` → `pixels.py` → `sagnac.py` → `report_figures.py` (22 October; simplified to library fits and plain variable names on 23 October, same numbers; `pixels.py` added on 24 October, a port of the presentation stage's fixed-pixel readout with the current calibration, because the instructor asked for that readout). Two earlier stages are kept:

- **Presentation (19 October)**: `presentation_slides/source/` reads `sagnac_update/` (the linear-calibration stage of 18 October) and the earlier `calibration.json`; its `numbers.json` was the source of the slides and of the earlier `RESULTS.md` (`previous_reports/`). It used a speed-reproducibility term ε, a split-timing formula for u(T), 4000 redrawn calibration lines and a wander estimate; the report pipeline gives the same slope and c without them. To rerun it, copy `previous_reports/calibration_eps_2026-10-13.json` to `calibration.json` (and the `.py` of the same name for `sagnac_update/`).
- **Earlier analysis (13 October)**: the quadratic calibration (`previous_reports/*_2026-10-13.*`, `figures/`).

## The physics, once

The prediction and the quantity measured (manual Eq. 5; Schwartz Eq. 1):

```
Δφ = 8π N A ω / (λ c) = a_th ω,     A = π R²,     a_th = 8π² N R² / (λ c)
```

The fringes are rings: two Gaussian beams on a common axis with wavefront radii R₁ ≠ R₂ at the camera have the phase difference

```
φ(r) = 2π r²/P + φ_c,     P = 2λ / |1/R₁ − 1/R₂|
```

so the dark rings sit at R_m² = P(m + η), with φ_c = π − 2πη. Rotation changes only φ_c. Every phase method reads φ_c or a phase averaged over part of the pattern; Δφ is that phase minus its mean at rest.

The rotation calibration: T the mean of the timed turns with u(T) from the lap scatter and the 0.25 s reaction time of the presses, ω = 2π/T, and the line ω = kS + ω₀ (weighted by the timing uncertainties) that converts the photo settings, with the line's prediction uncertainty as the uncertainty of ω.

## Naming conventions

- **Photos**: `Measurements/<S>/<YYYYMMDDhhmmss>.png`, S the motor-controller setting (`0` at rest); the time in the file name is the photo's time (seconds since midnight in `time_s`).
- **Per-photo tables**: `fringes_frames.csv` (center, ring radii, rms of the circle fit, P, φ_c by the minima and by the cosine fit), `methods_frames.csv` (template matching per photo), `pixels_frames.csv` (the fixed pixels' main draw per photo).
- **Per-setting table**: `sagnac_speeds.csv` (mean phase, SEM, first photo time, median P, Δφ, ω with its uncertainty, the residual, and the other methods' Δφ and SEM).
- **Results**: `calibration.json`/`.txt`, `sagnac.json`/`.txt` beside the scripts that wrote them.
- **Figures**: `report_figures/` (report); `figures/` (13 October), `sagnac_update/figures/` (18 October), `presentation_slides/` (slides).

## How the pieces connect

```
Rotation calibration.csv ──calibration.py──> calibration.json, calibration.txt ──────────┐
                                                                                         │
Measurements/<S>/*.png ──fringes.py──> fringes_frames.csv, fringes_profiles.npz ──┐       │
                                                                                  ▼       │
                                       methods.py ──> methods_frames.csv
                                       pixels.py ───> pixels_frames.csv, pixels.json
                                                                                  │       │
                                                                                  ▼       ▼
                                       sagnac.py ───> sagnac_speeds.csv, sagnac.txt, sagnac.json
                                                                                  │
                                       report_figures.py ───> report_figures/*.pdf, *.png
```
