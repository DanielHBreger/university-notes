# `report_figures.py`

Draws the report's figures in the report's style from the outputs of `calibration.py` and `sagnac.py`, so a number read off a figure is the number in `sagnac.txt` and `RESULTS.md`. It only draws; it fits nothing.

## Inputs

`calibration.json` (the timed settings and the line), `sagnac.json` (the regression, the methods, a_th), `sagnac_speeds.csv` (the ten setting means), `fringes_frames.csv` (every photo's center, radii and phase; ω per photo is recomputed from the line, 0 at rest) and, for `patterns`, `rings` and `extraction`, the representative photographs themselves. `methods` also reads the fixed pixels' columns of `sagnac_speeds.csv` and `fixed_pixels` of `sagnac.json`. No figure has a residual panel (points minus fit).

## Helpers

- `RAD_S`, `SLOPE_UNIT`: the units rad s⁻¹ and rad/(rad s⁻¹), with a thin space when LaTeX typesets the text.
- `panel_letter(ax, text, ...)`: the panel letter, inside the axes.
- `omega_text(omega)`: "ω = 0" or "ω = 1.83 rad s⁻¹".
- `legend_in_order(ax, labels)`: a legend in the given order, whatever order the artists were drawn in.
- `representative(S)`: the photo of setting S whose phase is closest to the setting's median (the same choice as the slides).
- `ring_radii(row)`: the finite radii among `R1`–`R3`.
- `crop(row, half)`: the red channel of a photo, a square of 2·half px around its fitted center, zero-padded past the image edge.
- `profile_against_r(row, r_max)`: `fringes.sector_profile` of the smoothed red channel about the photo's center, in 1 px annuli.
- `ring_line(row)`: the line R_m² = P(m + η) of one photo by the steps of `fringes.ring_phase` (`fringes.ring_numbers` for m).
- `save(fig, name)`: `report_style.save` into `report_figures/`, then closes the figure.

## The figures

- `setup()`: the top view of the apparatus, after Fig. 1 of the manual: turntable, the fiber loop as a bundle of circles, laser, P₁, beam splitter (diagonal as in the manual, so each beam is reflected once and transmitted once on its way to the camera), L₁ and L₂ to the two fiber ends, P₂, L₃, camera, the two beam directions, the sense of rotation, the radius R. Drawing only.
- `calibration()`: ω of the nine timed settings against S with their error bars (lap scatter and reaction time) on `line_axes()`: the fitted line with its ±1σ prediction band (`calibration.u_omega_pred`: the line's uncertainty and the rms scatter of the timed settings about it). One panel.
- `conversion()`: on the same `line_axes()`, dotted guides from each photo setting up to the line and across to the ω axis, with a bare error bar (no marker) at the line for the prediction uncertainty: the ω column of Table 2, drawn as a read-off rather than as data. One panel; it carries the same numbers as Table 2, so the report can use either.
- `patterns()`: one square panel per setting (the figure*, 512 pt wide), the representative photo cropped 500 px around its center, the fitted dark rings dashed, the center marked, the ω label from the line.
- `rings_figure()`: the rest and top-speed photos' intensity against r and against r² (1 px annuli, fringe sector), dashed lines at the fitted radii, arrows with the top-speed spacings (241 and 132 px; 1.08 and 1.08 × 10⁵ px²).
- `extraction()`: (a) the ring finder's minima on the top-speed photo (`fringes.analyse_photo` rerun on that photo), one color per ring, with the fitted circles; (b) R_m²/P against m at rest and at the top speed with the lines m + η, and the difference of the intercepts (0.52 fringe).
- `phase_vs_omega()`: the ten setting means against ω (error bars: the calibration line's uncertainty in ω, the per-setting uncertainty s in Δφ), the regression line `slope.a`, `slope.b`, and the prediction `theory.a_th` ± 1σ. One panel.
- `methods()`: (a) `dphi_cos − dphi`, `dphi_template − dphi` and `dphi_pixels − dphi` per setting; (b) the four slopes over the prediction's band: the first three with their regression errors, the fixed pixels with the SD of the slope between random draws of ten pixels (`fixed_pixels.u_draws`), which is why they are outside the method range.
- `drift()`: (a) ω of every photo against minutes since the first photo; (b) each photo's Δφ minus the line, with a straight line fitted to each setting's photos.

## Outputs and usage

`report_figures/<name>.pdf` and `.png` for `setup`, `calibration`, `conversion`, `patterns`, `rings`, `extraction`, `phase_vs_omega`, `methods`, `drift` (`random_pixels` was dropped on 24 October: the fixed pixels are in `methods`). Column figures are about 245 pt wide; `patterns` is 512 pt.

```
python report_figures.py                    # all, about a minute
python report_figures.py rings extraction   # only these
```

## Caveats

- With LaTeX, the text needs the `tgpagella` and `newpxmath` packages (and `dvipng` for the PNG); MiKTeX installs them on demand. Without LaTeX, matplotlib's STIX fonts are used and the figures look slightly different.
- If the analysis changes, rerun `calibration.py`, `sagnac.py` and this script in that order.
