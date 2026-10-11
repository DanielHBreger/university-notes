# `chua/identify.py`

Recovers the circuit that was actually on the bench from the oscillator records alone: node 1's capacitance and how it grows with the swing of v1 (with stage B's effective lag), C2, the inductor's inductance and loss as functions of current, and the nonlinear element as the circuit sees it, mapped onto Kennedy's two-op-amp realisation. Nothing here uses a nominal component value. The output `identified.json` is what `simulate.py --model static` and `--model bench` read.

Other outputs, beside the script: `identified.txt` (the readable report with the per-record table), `identified_inductor.png`, `identified_diode.png`, `identified_c1.png`.

## The idea

The periodic records (the period-1 cycle below the Hopf point, the small orbit near the origin, the large outer cycle) are averaged over their hundreds of cycles. Averaging reduces quantization scatter and leaves waveforms suitable for differentiation; it does not remove calibration uncertainty. The circuit equations then give each element directly:

- **node 1**, `C1 v1' = (v2 − v1)/Rt − g(v1)`. For any static element the loop integral of g dv1 over a cycle is zero, so `C1 = ∮ (v2 − v1)/Rt · v1' dt / ∮ v1'² dt`, whatever g is. With C1 known, `i_NR = (v2 − v1)/Rt − C1 v1'` against v1 is the element's curve as the circuit sees it, and any loop in it is a dynamic effect. If node 1's capacitance varies, the same ratio is its v1'²-weighted mean along the record, and over a long record the element's term is a vanishing boundary term even without a cycle; so it is taken over every record, and `node1_law` fits how it varies.
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

## Node 1's capacitance in every record

Two estimators of the loop integral that keep the scope's noise out of the denominator ∫ v1'² dt. Both were checked on the bench model run with constant C1 and sampled like the scope (400 ns or 1 µs, 5 mV of noise, the 12 to 80 mV steps of the ranges used): within 0.4 % from 0.15 V peak-to-peak up, in every regime. The earlier estimator (cycle markers on the raw v1, one averaged cycle) read the smallest cycles 2.8 % low at 0.15 V peak-to-peak, because markers on the noisy v1 select its noise and sharpen the average. A test recovers 11 nF through a cubic element from a noisy, quantised synthetic record with both.

### `cycle_node1(t, v1, v2, Rt)`

Periodic records. The cycles are marked on v1 low-passed at 20 kHz (`MARKER_FC`); the numerator comes from the average of all cycles, the denominator from the product of the derivatives averaged over the first half of the cycles and over the second half, whose noise does not correlate. Returns a `Node1(c, weights, v1, swing, method)` with the per-phase weights, v1 and swing, or `None` if the record cannot be averaged.

### `window_node1(t, v1, v2, Rt)`

Any record. Over windows of 20 µs (`WINDOW_S`), `∫_w (v2 − v1)/Rt dt = ΔQ + ∫_w g dt`, and `Σ y_w Δv_w / Σ Δv_w²` is the loop integral with every derivative a finite difference over the window: the noise enters Δv only through two end values (each averaged over three samples), and the element's term is a boundary term up to corrections at its kinks. The record is cut where v1 passes its starting value in the starting direction, which removes the boundary term the DC current through R would otherwise leave (2 % on a small cycle over 200 ms). On the large cycle the 20 µs windows straddle the kinks during the fast crossing and read 3.7 % high, so periodic records use `cycle_node1`.

### `node1_features(est, knots, bp)`

The estimator-weighted means of the quantities a law may depend on: `min(x, k)` and `max(x − k, 0)` for each candidate knot k (x the swing, `integration.swing`), the fraction of weight where stage B is linear (inside the inner breakpoints) and where stage A is saturated (beyond the outer ones).

### `node1_law(rows, small_labels, pwl, kennedy)`

```
C_rec = C_s + a0 <min(x, k)> + a1 <max(x − k, 0)> + s_B <B linear> + s_A <A linear>
```

fitted over every record. `s_A = tau_A·AA/RA` is stage A's lag at the datasheet speed (it is linear almost everywhere, so the records cannot tell it from `C_s`); `C_s` is the capacitor's share on the onset cycles, the smallest swing recorded; `a0`, `a1`, `s_B` are least squares, and the knot is the candidate (`LAW_KNOTS`, 0.3 to 1.2 V) with the least residual. A lag tau_B of stage B draws a current A_B·tau_B·v1'/R_B, a capacitance, but only while the stage is linear, so `s_B` gives stage B's effective lag `tau_B = s_B·R_B/A_B`. The fit is repeated on the periodic records alone, predicting the chaotic ones, as a check. Returns the law, its standard errors, the residual by group and the per-record values for the figure.

Why this law. The measured capacitance rises with the size of the swing on every kind of record (onset cycles, period-1 cycles, chaotic records, large cycle), and the rise is best described by the swing since the last turning point: the Rayleigh law of a ferroelectric dielectric (a class-2 ceramic capacitor's permittivity grows linearly with the AC field amplitude, by domain-wall motion), the same law the ferrite inductor follows. Laws in the voltage itself (a static C(v), the DC-bias effect) or in the excursion from a running mean fit worse, especially out of sample. The extra capacitance while stage B is linear is a step exactly at the element's breakpoints (a scan of the step's width puts it there), which is how a lag behaves and a smooth C(v) would not.

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

### `rayleigh_fit(I, L_eff, r_eff, bin_mA=1.0)`, `bin_weights(I, bin_mA=1.0)`

For `L = L0 + ν|d|` with the loss voltage `r0 i + ρ|d|d`, d the excursion of the current from its mean, driven by a sinusoid of amplitude I about that mean, the secant inductance is `L0 + νI/2` and the loss resistance is `r0 + (8/3π)ρI`. Two weighted straight-line fits recover the four numbers. The weights give every 1 mA amplitude bin the same total weight, so the forty-odd large-cycle records at one amplitude do not outvote the few small cycles that fix the small-signal end. A test recovers known values.

## The linearised circuit

### `negative_equilibrium(g, Rt, r0)`, `complex_pair(g, Rt, c1, c2, l0, r0)`

The negative equilibrium is the last root of `g(v) + v/(Rt + r0)` below zero; `complex_pair` takes the local slope G of g there by central difference, builds the 3×3 Jacobian and returns its eigenvalues with an imaginary part (None without an equilibrium).

### `small_signal_period(g, Rt, c1, c2, l0, r0)`

2π over the imaginary part of the complex pair.

### `small_signal_hopf(g, c1, c2, l0, r0)`

Scans Rt from 1200 to 3200 ohm for the sign change of the complex pair's real part and refines it with `brentq`.

## Per record: `analyse_record(path, rpot)`

1. Reads all three channels, finds the maxima with `maxima_from_arrays`, needs at least 60.
2. Cycle markers from the mid-level crossings; period and jitter from them.
3. Divider identity: fits `CH3 = A·CH1 + B·CH2` and stores `A + B`, which equals 1 when the CH1 and CH2 gains match.
4. `periodic = jitter < 1.5 % and spread of maxima < 0.12 V`. Non-periodic records (and periodic ones `cycle_node1` cannot average) return here with the summary numbers and `window_node1`'s estimate as `node1`.
5. Averaged cycle, low-passed, differentiated; `cycle_node1`'s estimate as `node1` and its capacitance as `C1_loop`; the cycle arrays and the admittance rows stored for later.

## `main()`

Each step below is a function of its own; `main` runs them in order and writes the outputs.

**1. Records** (`analyse_sweeps`). Every record with a clean divider fit in `forward_rpot.csv` and `back_rpot.csv` (every fourth with `--quick`, but always all of `forward/trace402` upward, where the small cycles are). The periodic ones are sorted into `small` (period-1 cycles: v1 max < 5.5, min < −1.5), `tiny` (the near-origin orbit: v1 max < 5.5, min ≥ −1.5) and `large` (the outer cycle: v1 max > 5.5).

**2. C1** (`loop_c1_values`). Loop-integral C1 values averaged over groups of periodic records: at small amplitude (the cycles `onset.py`'s A1² fit uses, A1 < 200 mV, the smallest swing recorded; this is the small-signal value the linearisation and L0 use), on the full period-1 cycle, and on the large cycle. They differ (10.8, 11.5, 12.6 nF): node 1's capacitance grows with the swing, which step 3b describes. The mean over every cycle on the shoulder (the "small-signal" value before 27 September) is kept as `shoulder_mean_loop`.

**3. The element** (`kcl_records`, `kcl_fit`, `five_segment_law`). Two KCL fits: the `core` records (forward sweep, 330 to 900 ohm, every other record) for the inner three lines, and the `large_cycle` records (forward below 315 ohm, back below 680) as a check. The saturation lines come from the plateaus of every large-cycle record, median over records. The outer breakpoints are where the shoulder lines through the inner breakpoints meet the saturation lines:

```python
cL = Ga * bpL - GbL * bpL                    # intercept of the left shoulder
bp_out_left = (c_left - cL) / (GbL - Gc_left)
```

Then `kennedy_from_pwl` with the averaged Gb and Gc.

**3b. Node 1's law** (`node1_law`), over every record, with the breakpoints and gains from step 3. Stored in the JSON as `C1_nF.law`, and what `simulate.py --model bench` integrates.

**4. C2, the inductor, the closure** (`tank_fits`, `inductor_loops`, `rayleigh_fit`, `calibrate_l0`). A joint admittance fit over the period-1 cycles with at least five harmonics gives a first C2. That fit pins L·C2 but splits it poorly: at 3 kHz and 2 mA the inductor is already nonlinear over the harmonics it uses. The near-origin orbit (1.4 kHz, 0.7 mA, ten harmonics) is where the inductor is most nearly linear and the harmonics reach highest, so C2 is taken from it when it exists (90.5 nF). With that C2 every periodic record gives a flux loop and `rayleigh_fit` gives L0, ν, r0, ρ over records with I > 0.6 mA (below that the loop is a few quantisation steps wide).

The small-signal end of the inductance law is then closed on the measured onset frequency: the flux loops of the smallest cycles are the least certain numbers here. L0 is solved so that the linearised circuit (small-signal C1, C2, r0, the fitted element) oscillates at the frequency of the first limit cycle of the forward sweep at its resistance (`onset_records.csv`: forward trace14, 3045 Hz at 904.9 ohm), and ν is refitted to the loops with that L0 (same bin weights). (Until 27 September the mean period of the 18 smallest cycles was matched at their mean rpot.py resistance, a finite-amplitude period at a resistance rpot.py reads low.) The records near the onset enter every step at `onset.py`'s resistances (`sweep_records`). The loop-only values are kept in the JSON as `L0_loops_H`, `nu_loops_H_per_A` for comparison.

A check follows: the model's period at the onset resistance and its Hopf point from the linearisation (the measured onset is `onset.py`'s).

**5. Dynamic deviation** (`dynamic_deviation`). On the highest-R large-cycle record, `i_NR` from the measured cycle (node 1 at the record's own loop-integral capacitance) minus the static law, inside |v1| < 3 V: the peak deviation in mA and the v1 slew rate at which it occurs. This is the number that says the diode's op-amps cannot follow the fast crossings. (Before the node-1 law the small-signal C1 was used here, which booked the large cycle's extra capacitance as element current.)

**Report and figures** (`report_text`, `plot_inductor`, `plot_diode`, `plot_c1`). The JSON (fields listed in the summary at the top of `identified.txt`), the text report with the per-record inductor table, and three figures: flux loops at four amplitudes with L_eff and r_eff against amplitude and the Rayleigh lines; the element as seen from a period-1 record and a large-cycle record (each at its own C1) over the V–I trace, with both spline fits and the breakpoints; and node 1's capacitance per record against Rpot with the law, and the capacitor's share against the swing.
