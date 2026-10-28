# `pixels.py`

The fixed-pixel readout of preparation question 5, as a fourth way to read the fringe phase (24 October 2026; a port of the presentation stage's `sagnac_update/methods.py` and `presentation_slides/source/random_points_scan.py` into the current pipeline, with the current calibration). Run after `fringes.py`; `sagnac.py` fits its output like the other readouts and the methods figure shows it with the spread between draws as its uncertainty.

**The idea.** Ten pixels at fixed positions in the image are read in every photograph. As the rings move, each pixel brightens and dims: I_kj = a_j + b_j cos(ψ_j + δ_k) for pixel j and photograph k, with ψ_j = 2π s_j/P₀ the pixel's fringe phase at rest (s_j its r² from the rest ring center, P₀ the rest spacing), δ_k the photograph's phase (the unknown; the same convention as the ring readout, the phase at the ring center) and a_j, b_j the pixel's offset and amplitude. One alternating least-squares fit finds all of them.

**Why the fit needs limits.** Over this run the pattern moves by half a fringe, so no pixel passes through a whole bright–dark cycle. Without limits, ten such pixels cannot tell a large phase change at low contrast from a small one at high contrast. Two physical limits make the fit well posed: b_j ≥ 0, and a_j − b_j ≥ the camera's dark level (the darkest point of a pixel's fringe cannot be darker than the camera's black, the median of the four 100 × 100 px image corners of the mean rest image, 2.0 counts).

**Why many draws.** A fixed pixel reads the local phase at its own radius, so a move of the ring center or a change of the spacing between photographs reads as a phase change; which pixels are drawn therefore matters. The fit is repeated for 1000 random draws of ten pixels (seeds 1–1000) besides the main draw (seed 20261011, the measurement date, fixed before any result), and for 60 draws each of 20, 40 and 80 pixels. The spread (SD) of the slope between draws is the method's uncertainty; fits ending with a negative slope are counted as failed.

## Functions

- `fringe_region(shape, cx, cy)`: the (y, x) of every pixel inside the fringe sector `SECTOR_DEG` with r² < `S_MAX` about the rest center (both from `fringes.py`).
- `draw_pixels(ys, xs, seed, n)`: n different pixels of the region, every pixel equally likely.
- `dark_level(img)`: the median of the four image corners.
- `best_phases(I, a, b, psi)`: with a, b held, the phase of every photograph: the best point of a 1° grid over (−π, π], refined by a parabola through it and its neighbours.
- `offsets_amplitudes(I, X, dark)`: with the phases held (X = cos(ψ_j + δ_k)), a_j and b_j by least squares under the two limits. With c = a − b the model is I = c + b(1 + X), a two-parameter least squares with bounds: the best allowed point is the free solution if allowed, else the best point on one of the two limits, compared per pixel.
- `fit_phases(I, psi, dark)`: the alternation, until no phase moves by more than `TOL` = 10⁻⁶ rad (at most `MAX_ROUNDS` = 500). Returns the wrapped phases, a, b, the rms residual in counts and the number of rounds.
- `relative_to_rest(S, time_s, delta)`: the wrapped phases unwrapped in order of S, then time, and made relative to the rest mean.
- `slope_of(S, omega, phases)`: the setting means against ω by an equal-weight line, as `sagnac.py` fits the other readouts (the per-draw slopes).
- `one_draw(task)`: fit and regress one draw; the draws run in a process pool.

## Outputs

- `pixels_frames.csv`: `file`, `phi_pixels_rad` (the main draw's phase of every photograph, relative to rest).
- `pixels.json`: `main` (seed, n, the pixels' (x, y), dark level, center, P₀, slope, fit error, rms, rounds), `summary` (per number of pixels: draws, failed, mean, SD, median of the slope, median fit error) and `draws` (every draw's seed, n, slope, fit error, rms, rounds).

## Numbers of 24 October

Main draw: slope 1.7734 ± 0.0367 (fit), rms residual 10.2 counts, 116 rounds. Ten pixels, 1001 draws (6 failed): mean 1.689, SD 0.193 (11.4 %), median fit error 0.057. With 20, 40 and 80 pixels the SD falls to 0.090, 0.072 and 0.055 (5.3, 4.2, 3.2 %), still above the fit error of one draw. These agree with the presentation stage's scan, which used the earlier calibration line.

`sagnac.py` reports the main draw through the same averaging and regression (1.7734 ± 0.0367, s = 0.069 rad, R² = 0.9966; per setting 0 to −0.12 rad below the minima) and leaves the readout out of the method range: its slope depends on which pixels are drawn, so it cannot constrain a at the 3 % level. Test: `tests/test_pixels.py` (synthetic intensities with a known phase sweep).
