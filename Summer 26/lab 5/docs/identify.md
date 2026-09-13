# `chua/identify.py`

Recovers the circuit that was actually on the bench from the oscillator records alone: C1, C2, the inductor's inductance and loss as functions of current, and the nonlinear element as the circuit sees it, mapped onto Kennedy's two-op-amp realisation. Nothing here uses a nominal component value. The output `identified.json` is what `simulate.py --model static` and `--model bench` read.

Other outputs, beside the script: `identified.txt` (the readable report with the per-record table), `identified_inductor.png`, `identified_diode.png`, `identified_c1.png`.

## The idea

The periodic records (the period-1 cycle below the Hopf point, the small orbit near the origin, the large outer cycle) are averaged over their hundreds of cycles. Averaging removes the scope's quantisation and leaves waveforms clean enough to differentiate. The circuit equations then give each element directly:

- **node 1**, `C1 v1' = (v2 − v1)/Rt − g(v1)`. For any static element the loop integral of g dv1 over a cycle is zero, so `C1 = ∮ (v2 − v1)/Rt · v1' dt / ∮ v1'² dt`, whatever g is. With C1 known, `i_NR = (v2 − v1)/Rt − C1 v1'` against v1 is the element's curve as the circuit sees it, and any loop in it is a dynamic effect.
- **node 2**, `(v1 − v2)/Rt = C2 v2' + iL` and `L iL' = v2 − r iL`. C2 comes from the tank admittance `(v1 − v2)/Rt ÷ v2` harmonic by harmonic. With C2 known, `iL = (v1 − v2)/Rt − C2 v2'`, and the flux `∫ v2 dt` against iL is the inductor's own loop: secant slope = effective inductance, enclosed area = loss per cycle.
- **the element**: the five-segment law fitted to all records by integrated KCL with a spline basis, saturation segments from the plateaus of the large cycle, and the eight segment numbers mapped onto Kennedy's circuit.

## Cycle averaging

### `cycle_markers(t, v1, T_est)`

Times at which v1 crosses the middle of its range upward, at least 0.6·T_est apart, with linear interpolation between samples. The crossing is where v1 moves fastest, so it marks a cycle far more sharply than a maximum on a flat top; on the large cycle, whose maxima sit on a plateau, peak-based markers produced spurious jitter. A test checks the markers of a sinusoid are one period apart.

### `averaged_cycle(t, v, t_peaks, nph=2048, tol=0.05)`

For every cycle whose length is within 5 % of the median, resamples each column of `v` onto 2048 phase points by interpolation and accumulates. Returns the mean cycle, the mean period and the number of cycles used, or `None` if fewer than 20 qualify.

### `harmonics`, `derivative`, `lowpass`

Fourier tools on the averaged cycle. `derivative(x, T, kmax=80)` multiplies each harmonic by `i·2πk/T` and zeroes everything above harmonic 80, so the derivative is exact for the periodic signal and does not amplify the residual noise. `lowpass` keeps the same 80 harmonics.

## Element identification per record

### `loop_c1(v1, v2, dv1, Rt)`

The loop-integral formula above, as two sums over the phase points. Exact for a static element of any shape; a test recovers 11 nF from a synthetic cycle to 1e-9.

### `inductor_loop(v1, v2, dv2, Rt, C2, T)`

```python
iL = (v1 - v2) / Rt - C2 * dv2
iL = iL - iL.mean()             # the channel offsets, not the tiny DC drop
phi = np.cumsum(v2c) * dt; phi -= phi.mean()
I = 0.5 * np.ptp(iL)
L_eff = np.ptp(phi) / np.ptp(iL)
r_eff = np.sum(v2c * iL) / np.sum(iL * iL)
```

Returns the current, the flux, the current amplitude, the secant inductance of the loop and the loss resistance from the in-phase power.

### `admittance_rows(v1, v2, Rt, T, kmax=24, vmin=2e-4)`

`(ω, Y, |V2|)` per harmonic up to 24, keeping only harmonics where v2 has at least 0.2 mV; the third element is the weight.

### `fit_tank(rows)`

Least squares of `Y = jωC2 + 1/(r + jωL)` over all rows, real and imaginary parts stacked, in units of nF, mH and ohm with bounds. Returns `(C2, L, r)` and the relative residual. A test recovers a synthetic tank to 2 %.

## The element from integrated KCL

### `spline_design(v)`, `KNOTS`, `NB`

A cubic B-spline basis on −8 to 7.4 V with 0.2 V knot spacing. g(v1) is represented as `B(v1) · c`.

### `kcl_normal_equations(path, rpot, window_s=30e-6, stride=25, nrows=400000)`

Integrating node 1 over a window w of 30 µs:

```
∫_w (v2 − v1)/Rt dt = C1 [v1]_w + ∫_w B(v1) dt · c
```

Integration is what makes raw, quantised records usable: it averages the quantisation and turns the derivative into a difference. The design row for each window is `[v1(end) − v1(start), IB(end) − IB(start)]` with `IB` the cumulative integral of the basis columns, and the target is the integrated coupling current. Only the normal equations `X'X, X'y, y'y` are returned (scaled by the window count), so hundreds of records accumulate into one small system.

### `solve_kcl(XtX, Xty, yty, lam=1e-9)`

Solves the accumulated system with a second-difference penalty on the spline coefficients (not on C1). Returns the apparent C1, the coefficients and the relative residual.

### `pwl_from_spline(cs)`

Reads the five-segment numbers off the fitted spline: the slopes Gb_left, Ga, Gb_right, Gc_left, Gc_right as straight-line fits over fixed voltage windows, the inner breakpoints where the spline's derivative crosses the midpoint between Ga and the mean Gb, the outer ones where it crosses 1.5 mS. Also g(0).

### `saturation_from_plateau(v1, v2, dv1, Rt, C1, side, frac=0.03)`

The saturation segments are known better from the large cycle's plateaus than from the spline's edges. While v1 sits on the steep segment v1' is tiny, so `i_NR = (v2 − v1)/Rt − C1 v1'` is known without any model of the element, and it moves along the segment as v2 swings. A line through the plateau samples (|dv1| below 3 % of its maximum, |v1| > 4 V) gives the slope and intercept.

### `kennedy_from_pwl(Ga, Gb, Gc, bp_in_left, bp_in_right, bp_out_left, bp_out_right, RB=22e3)`

Inverts the formulas in `integration.static_pwl`. With RB fixed at its nominal 22 kΩ the three slopes give RA, AA, AB exactly, and the breakpoints times the gains give the rails of each stage. A test checks the round trip.

### `rayleigh_fit(I, L_eff, r_eff, bin_mA=1.0)`

For `L(i) = L0 + ν|i|` and `r(i) = r0 + ρ|i|` driven by a sinusoid of amplitude I, the secant inductance is `L0 + νI/2` and the loss resistance is `r0 + (8/3π)ρI`. Two weighted straight-line fits recover the four numbers. The weights give every 1 mA amplitude bin the same total weight, so the forty-odd large-cycle records at one amplitude do not outvote the few small cycles that fix the small-signal end. A test recovers known values.

## Small-signal checks

### `small_signal_period(g, Rt, c1, c2, l0, r0)`

Finds the negative equilibrium (root of `g(v) + v/(Rt + r0)`), the local slope G of g there by central difference, the 3×3 Jacobian, and returns 2π over the imaginary part of the complex pair.

### `small_signal_hopf(g, c1, c2, l0, r0)`

Scans Rt from 1200 to 3200 ohm for the sign change of the complex pair's real part and refines it with `brentq`.

## Per record: `analyse_record(path, rpot)`

1. Reads all three channels, finds the maxima with `maxima_from_arrays`, needs at least 60.
2. Cycle markers from the mid-level crossings; period and jitter from them.
3. Divider identity: fits `CH3 = A·CH1 + B·CH2` and stores `A + B`, which equals 1 when the CH1 and CH2 gains match.
4. `periodic = jitter < 1.5 % and spread of maxima < 0.12 V`. Non-periodic records return here with only the summary numbers.
5. Averaged cycle, low-passed, differentiated; the loop C1; the cycle arrays and the admittance rows stored for later.

## `main()`

**1. Records.** Every record with a clean divider fit in `forward_rpot.csv` and `back_rpot.csv` (every fourth with `--quick`, but always all of `forward/trace402` upward, where the small cycles are). The periodic ones are sorted into `small` (period-1 cycles: v1 max < 5.5, min < −1.5), `tiny` (the near-origin orbit: v1 max < 5.5, min ≥ −1.5) and `large` (the outer cycle: v1 max > 5.5).

**2. C1 and a first C2.** A joint admittance fit over the period-1 cycles with at least five harmonics gives a first C2. Three loop-integral C1 values are averaged: at small amplitude (records whose v1 stays below −2 V, the shoulders only), on the full period-1 cycle, and on the large cycle. They differ (11.1, 11.5, 12.6 nF) because the element is not static; the growth is the op-amp lag.

**3. The element.** Two KCL fits: the `core` records (forward sweep, 330 to 900 ohm, every other record) for the inner three lines, and the `large_cycle` records (forward below 315 ohm, back below 680) as a check. The saturation lines come from the plateaus of every large-cycle record, median over records. The outer breakpoints are where the shoulder lines through the inner breakpoints meet the saturation lines:

```python
cL = Ga * bpL - GbL * bpL                    # intercept of the left shoulder
bp_out_left = (c_left - cL) / (GbL - Gc_left)
```

Then `kennedy_from_pwl` with the averaged Gb and Gc.

**4. C2, the inductor, the closure.** The period-1 fit pins L·C2 but splits it poorly: at 3 kHz and 2 mA the inductor is already nonlinear over the harmonics it uses. The near-origin orbit (1.4 kHz, 0.7 mA, ten harmonics) is where the inductor is most nearly linear and the harmonics reach highest, so C2 is taken from it when it exists (90.5 nF). With that C2 every periodic record gives a flux loop and `rayleigh_fit` gives L0, ν, r0, ρ over records with I > 0.6 mA (below that the loop is a few quantisation steps wide).

The small-signal end of the inductance law is then closed on the measured onset period: the flux loops of the smallest cycles are the least certain numbers here, while the period at onset is measured to 0.3 %. L0 is solved so that the linearised circuit (small-signal C1, C2, r0, the fitted element) oscillates at the onset period at the onset resistance, and ν is refitted to the loops with that L0 (same bin weights). The loop-only values are kept in the JSON as `L0_loops_H`, `nu_loops_H_per_A` for comparison.

Two checks follow: the model's period and Hopf point from the linearisation, and the measured Hopf point by extrapolating the squared amplitude of the small cycles (0.15 to 1 V) linearly in R to zero.

**5. Dynamic deviation.** On the highest-R large-cycle record, `i_NR` from the measured cycle minus the static law, inside |v1| < 3 V: the peak deviation in mA and the v1 slew rate at which it occurs. This is the number that says the diode's op-amps cannot follow the fast crossings.

**Report and figures.** The JSON (fields listed in the summary at the top of `identified.txt`), the text report with the per-record inductor table, and three figures: flux loops at four amplitudes with L_eff and r_eff against amplitude and the Rayleigh lines; the element as seen from a period-1 record and a large-cycle record over the V–I trace, with both spline fits and the breakpoints; and the loop C1 against cycle amplitude.
