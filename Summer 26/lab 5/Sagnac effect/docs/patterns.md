# `patterns.py`

> **Earlier stage.** `patterns.py` still runs (it reads `sagnac_speeds.csv`, whose `omega` column is now the line's value), but the report uses `report_figures.py` for the pattern figure.

The two task-2 figures of the earlier analysis: one photograph per setting with its fitted dark rings (`figures/patterns`), and each setting's mean intensity profile against r² (`figures/profiles`). The report now uses `report_figures.py` for its figures instead.

Inputs: the photographs, `fringes_frames.csv` (centers, radii and phases from `fringes.py`) and `sagnac_speeds.csv` (ω per setting, from the top-level `sagnac.py`). Run after both. Both figures are written as PDF and 600 dpi PNG through `report_style.py`.

## `representative(fr, S)`

Picks the photograph shown for setting S: the one whose phase is closest to the setting's median phase, so the panel shows a typical photograph rather than the first.

```python
ph = np.angle(np.exp(1j * (g.phi_rad - np.angle(np.mean(np.exp(1j * g.phi_rad))))))
return g.iloc[int(np.argmin(np.abs(ph - np.median(ph))))]
```

The phases are first taken relative to the setting's circular mean (the direction of the average of the unit vectors e^(iφ)) and wrapped, so that a setting whose phases straddle ±π is not split in two before the median is taken. In this run no phase is near ±π (`phi_rad` between −0.84 and 2.57 rad in `fringes_frames.csv`), so this step changes nothing. With an even number of photographs the median falls halfway between the two middle ones, both are equally close, and floating-point rounding decides which is shown. That is the case at rest (10 photographs) and at S = 44.07, 55.83, 64.00 and 73.65 (22, 18, 16, 16).

## `unit(tex)`

The unit rad s⁻¹ for labels, written for LaTeX or for matplotlib's own math text, depending on whether `report_style` found a LaTeX installation.

## `gallery(fr, sp, path)`

`figures/patterns`, full text width (figure*): 2 rows of 5 panels, one per row of `sagnac_speeds.csv` (rest and the nine settings, in increasing S). Each panel shows

- the red channel of the representative photograph, unsmoothed, in gray from 0 to 220 counts (brighter pixels show white; the brightest pixel in the data is 239);
- a 940 px square crop around the fitted center (`HALF` = 470 px each way), cut at the rounded center but drawn with `extent` in image pixel coordinates, so the overlays can use the unrounded fitted values:

  ```python
  ax.imshow(crop, cmap='gray', vmin=0, vmax=220, interpolation='antialiased',
            extent=(x0, x0 + crop.shape[1], y0 + crop.shape[0], y0))
  ```

- the fitted dark rings `R1`–`R3` of that photograph as dashed circles and the center as +, both in `rs.COLORS[4]` (orange);
- the label ω = … rad s⁻¹ with two decimals from the `omega` column of `sagnac_speeds.csv`, or ω = 0 at rest.

The y axis runs from cy + HALF at the bottom to cy − HALF at the top, so the photograph stays upright (image rows count downward). Panels beyond the number of settings would be switched off; with ten settings there are none.

## `profiles(fr, sp, path)`

`figures/profiles`, column width: for each setting, the intensity averaged over the fringe sector against r² about each photograph's own fitted center (`sector_profile_s2` from `fringes.py`), then averaged over all photographs of the setting. Evenly spaced fringes along r² that slide to smaller r² as ω grows are the picture of the Sagnac shift.

```python
edges = np.arange(0, (R_MAX - 60) ** 2, S2_BIN)
```

Bins of 2000 px² out to 336 000 px² (r about 580 px), further out than the cosine fit's 2.6 × 10⁵ px². Unlike `fringes.py`, it reads the unsmoothed red channel, and it opens every photograph (161 files). The curves are colored by ω on the `Blues` map, normalized from −0.45 ω_max to ω_max so that the rest curve is not the palest blue; the color bar shows only 0 to ω_max, with ticks every 0.5 rad/s. x is in units of 10⁵ px², intensity in counts.

## `main()`

Reads `fringes_frames.csv` and `sagnac_speeds.csv` from the script's folder and writes `figures/patterns` and `figures/profiles`.

## Constants

| name | value | meaning |
|---|---|---|
| `HALF` | 470 | px, half-width of the square crop in `gallery` |
| `R_MAX` | 640 | px, imported from `fringes.py`; `profiles` goes out to `R_MAX` − 60 |
| `S2_BIN` | 2000.0 | px², imported from `fringes.py`, bin width in r² |

## Caveats

- **The ω values are from the earlier, quadratic analysis.** `main` reads the top-level `sagnac_speeds.csv`, which the top-level `sagnac.py` writes with the quadratic calibration: its `omega` column (0.2383 … 1.7816 rad/s) equals the quadratic conversion printed in `calibration.txt`. The final analysis uses the linear calibration, `sagnac_update/sagnac_speeds.csv` (0.2354 … 1.8135 rad/s, as in RESULTS.md, Task 3). At two decimals five of the nine labels differ: 0.41, 0.62, 1.47, 1.63, 1.78 here against 0.40, 0.61, 1.48, 1.65, 1.81. The same column sets the colors and the top of the color bar in `profiles`. The existing `figures/patterns.png` and `figures/profiles.png` carry these values.
- For three settings the photograph shown is not the one in the slides and the report. Those use `representative` in `presentation_slides/source/deck.py`, which compares the plain phases with the pandas median and breaks the tie the other way at rest, S = 55.83 and S = 73.65.
- The dashed circles are the fitted rings of the one photograph shown, not a setting average.
