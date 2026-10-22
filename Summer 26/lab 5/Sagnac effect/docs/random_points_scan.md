# `presentation_slides/source/random_points_scan.py`

> **Presentation stage.** The stored `random_points_scan.json` is summarized by `sagnac.py` (`scan_summary`) and drawn by `report_figures.py`; its slopes were fitted with the presentation's weighted line, not the report's regression (the difference is small, see below). The script itself needs the earlier `calibration.json` and `sagnac_update/` to rerun.

Repeats the random-pixel method of `sagnac_update/methods.py` for 1000 draws of ten pixels, plus 60 draws each of 20, 40 and 80 pixels, and writes every draw's slope to `random_points_scan.json`. It settles how much the slope from one draw of ten pixels depends on which pixels were drawn (slide 10, `RESULTS.md`).

## The idea

The method (see `methods.md`): ten pixels at fixed image positions are read in all 161 photos, and pixel j in photo k is modeled as I = a_j + b_j cos(ψ_j + δ_k), with ψ_j = 2πs_j/P₀ the pixel's phase at rest. An alternating least-squares fit finds every a_j, b_j and photo phase δ_k, and the phases go through the per-setting averaging and the line. One draw reports its slope with a fit error, but that error covers only the noise of that draw's fit, not the luck of which pixels were drawn. Repeating the whole method with many seeds measures the second directly: the SD of the slopes over the draws, set against the typical fit error of one draw. Draws with more pixels show whether the spread shrinks toward the fit error.

The script changes nothing in the method: it imports `sagnac_update/methods.py` and calls its functions. The 31 draws that `sagnac_update/methods.json` also holds (seed 20261011 and seeds 1–30) have the same slopes there to within 10⁻⁷.

## Paths and constants

```python
DATA = Path(os.environ.get('SAGNAC_DATA', HERE.parents[1]))          # the 'Sagnac effect' folder
sys.path[:0] = [str(DATA / 'sagnac_update'), str(DATA)]
import methods as M
```

`sagnac_update` comes first, so `methods` and the `sagnac` it imports are the linear-calibration versions; `fringes.py` and `calibration.py` come from the main folder. `N_DRAWS` is the first command-line argument (default 1000), `EXTRA_N = [(20, 60), (40, 60), (80, 60)]` (pixels per draw, number of draws) and `WORKERS` the number of CPUs. `G` is a module-level dict that `main()` fills and the worker processes inherit.

## `linear_model()`

The linear calibration with ε from `calibration.json` (`fits_with_eps['1']`, ε from `conversion`) as the dict `{degree, coef, cov, eps}` that `methods.slope` expects. It is the same line as in `deck.py`.

## `main()`

1. Reads `fringes_frames.csv` and makes the file paths absolute.
2. Rest geometry: `cx0`, `cy0` and `P0` are the medians of the ten rest photos' fitted center and period. `rest_img` is the mean of their red channels, smoothed with σ = `SMOOTH_PX` = 4 px. `dark = M.dark_level(rest_img)` is the median of the four 100 × 100 px image corners. `M.fringe_region` gives every pixel in the sector −75° to 180° with r² < 2.6 × 10⁵ px² about the rest center.
3. The plan of draws:

```python
plan = [(M.POINTS_SEED, 10)] + [(sd, 10) for sd in range(1, N_DRAWS + 1)]
for n, k in EXTRA_N:
    plan += [(100000 * n + sd, n) for sd in range(1, k + 1)]
```

The original draw (seed 20261011) first, then seeds 1–1000, then seeds 2000001–2000060, 4000001–4000060 and 8000001–8000060 for 20, 40 and 80 pixels, so that no larger draw reuses a ten-pixel seed. `M.draw_points` picks n different pixels of the region, each equally likely, with `np.random.default_rng(seed)`.

4. One pass over the photographs: each is smoothed once and sampled at every drawn pixel into `V_all` (photos × all pixels, `float32`). `offsets` marks where each draw's pixels start.
5. A pool of `WORKERS` processes started by `fork` runs `one_draw` for every draw; `imap` keeps the plan order. Every 100 rows, the rows so far are written to `random_points_scan.partial.json` as a bare list, which `deck.py` also accepts (through `SAGNAC_SCAN`).
6. Writes `random_points_scan.json`.

## `one_draw(i)`

```python
psi_rest = 2 * np.pi * ((d[1] - G['cx0']) ** 2 + (d[0] - G['cy0']) ** 2) / G['P0']
psi_frame = 2 * np.pi * ((d[1][None] - G['CX']) ** 2 + (d[0][None] - G['CY']) ** 2) / G['PP']
```

`d` holds the (y, x) of the draw's pixels. `'fixed'` takes each pixel's phase at rest from the rest center and P₀: the method as used. `'frame'` is a diagnostic, the same pixels with ψ from each photo's own fitted center and period. It removes the two things a fixed pixel cannot tell from a phase change, a moving center and a changing P.

For each, `M.point_phases` runs the alternating fit and makes the phases continuous and relative to the rest mean. `M.slope` then makes a slope as `sagnac_update/sagnac.py` does: per-setting means (`setting_means`), ω from the line with 0 at rest, and the weighted line of `fit_line` (per-point error from SEM, aεω and σ_x set by χ²/dof = 1). The row keeps the slope, its fit error `u_fit`, `sigma_x`, the rms of the pixel fit and the number of rounds.

## `random_points_scan.json`

```
{"note": "...", "dark": 2.0000000000000004, "centre": [897.81, 507.685], "P0": 107546.85, "max_iter": 500,
 "draws": [{"seed": 20261011, "n": 10,
            "fixed": {"slope": ..., "u_fit": ..., "sigma_x": ..., "rms": ..., "iterations": ...},
            "frame": {the same keys}},
           ...]}
```

1181 draws in plan order (1 + 1000 + 3 × 60). `slope` and `u_fit` are in rad/(rad/s), `sigma_x` in rad, `rms` in counts; `iterations` = `max_iter` means the fit stopped at the limit. `dark` (counts), `centre` (px) and `P0` (px²) are the rest geometry above, the same values as in `sagnac_update/methods.json`.

`slides.scan_stats` reads the file (see `final_numbers.md`). Results (`numbers.json` key `scan`, `RESULTS.md`): over the 998 good ten-pixel draws, mean slope 1.702(6), SD 0.196(5) and median fit error 0.057, so the spread is 3.4 times one draw's fit error; 2 fits failed (negative slope) and 18 stopped at 500 rounds. With 20, 40 and 80 pixels the SD falls to 0.091, 0.073 and 0.056 while the fit error stays at 0.041, 0.040 and 0.037. The `'frame'` diagnostic moves the mean to 1.78, with SD 0.17.

## Caveats

- The docstring says the script "Runs in a few minutes"; `make_all.py` and `presentation_slides/README.md` say about an hour on two cores.
- `mp.get_context('fork')` exists only on Linux and macOS. On Windows the script stops at that line.
- Each draw's slope comes from the σ_x-weighted line of `sagnac_update/sagnac.py`, not from the equal-weight regression of `fit8`, although `RESULTS.md` and the presenter guide say "the same regression". The two nearly coincide when σ_x dominates every point's error (seed 20261011: σ_x = 0.068 rad, `sem_points` at most 0.016 rad in `sagnac_speeds.csv`).
- `u_fit` has no calibration-line part: `methods.slope` asks `calibration_part` only for ω.
- Two things differ from `main()` of `methods.py`: the pixel values are stored as `float32` (`methods.py` keeps `float64`), and the calibration model is built from `fits_with_eps['1']` of the top-level `calibration.json`, where `methods.py` reads the `conversion` entry of a `calibration.json` in its own folder. The slopes still agree with `methods.json` to 10⁻⁷.
- A fit that converges in exactly 500 rounds is also counted as stopped at the limit.
