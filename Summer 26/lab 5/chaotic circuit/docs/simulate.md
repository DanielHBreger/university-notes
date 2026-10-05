# `chua/simulate.py`

Simulates the oscillator with three models of the same circuit and puts each beside the measurements: phase portraits next to the nearest record, a bifurcation sweep over the measured points, and the transition resistances.

Outputs beside the script, `--tag` adding a suffix: `simulated_nr<tag>.png` (the element and load lines), `simulated_portraits<tag>.png`, `simulated_timeseries<tag>.png`, `simulated_bifurcation<tag>.png`, `simulated_sweep<tag>.npz` (every maximum of both sweeps, for the report figure), `simulated_transitions<tag>.json` (the Table 2 rows and the regime runs), optionally `simulated<tag>/` with scope-format CSVs. The printout is kept as `simulate<tag>.txt`.

## The three models

- `--model ideal`: the plan's model. Constant nominal C1 = 10 nF, C2 = 100 nF, L = 18 mH; the element is Table 1's: the five slopes and four breakpoints fitted to the V–I record (`diode_fit.json`), joined into the continuous law through the origin (`table1_g`). `--element intersections` joins the fitted lines where they intersect instead (inner breakpoints −1.26/+1.20 V rather than −1.09/+0.58 V); with C1 = 11.5 nF that is the model of the presentation's animations.
- `--model static`: the identified circuit (`identified.json`) with everything held constant: small-signal C1, C2, L0, r0 and the identified five-segment law, integrated like the ideal model. It isolates what the right numbers do without the right physics.
- `--model bench`: the identified circuit in full: node 1's capacitor growing with the swing of v1 since its last turning point (the Rayleigh law of a ferroelectric ceramic, `identify.py`'s node-1 law), Rayleigh inductor acting on the current's excursion from its running mean, Kennedy diode with op-amp lag and slew limit, integrated with the bench kernel. Stage B's lag is the one the node-1 law measures (a lag acts as a capacitance A_B·tau_B/R_B at node 1 while the stage is linear, which is how the loop integrals see it); stage A's lag and the slew rate are the TL082 datasheet values (`TAU_A = 0.059 µs` at 3 MHz gain-bandwidth, `SLEW = 13 V/µs`). `TAU_A` adds tau_A·AA/RA = 0.26 nF to node 1, which is subtracted from the measured C1 (see `bench_params`). `--constant-c1` runs the model as it was before the node-1 law: constant small-signal C1 and the datasheet stage-B lag `TAU_B = 0.40 µs`. `--datasheet-lag` keeps the capacitor law but puts stage B at `TAU_B`: the report's model (b) (`--tag _bench_datasheet`), which leaves the measured lag out; `report_figures.py` draws it. `--tau-b` and `--slew` override the op-amps.

## The circuit: `Circuit`

A frozen dataclass of the constant components, SI units: `c1`, `c2`, `l`, the inductor's series resistance `rl`, and the fixed resistor `r0` in series with the potentiometer; `circuit.rt(rpot)` is the total coupling resistance R0 + Rpot. The defaults are the plan's nominal values, and `NOMINAL = Circuit()`. Every function that needs component values takes a `Circuit` explicitly; nothing reads or changes module-level component values. The bench model's linearisation is `bench_small_signal(P)`, also a `Circuit`.

## The ideal element

- `load_vi_fit`: reads the segments from `diode_fit.json` and refuses a fit not made against the element voltage. `vi_segments()` reads it once and caches it.
- `fitted_offset`: g(0) of the raw fit, the intercept of the segment containing v = 0. One scope code of channel offset.
- `pwl_from_segments(segments, i0)`: subtracts the offset from every intercept and puts each breakpoint where adjacent lines intersect, `(b2 − b1)/(m1 − m2)`. Raises if the intersections are out of order.
- `scale_inner_breakpoints`, `symmetrize_inner`: the `--bp-scale` and `--symmetric` experiments.
- `make_g(...)`: the line-intersection element as a callable (`--element intersections`); `table1_g()`: Table 1's slopes and breakpoints as the continuous law through the origin (the default). `pwl_g(bp, G)` builds the callable from knots (`integration.pwl_knots`) and attaches `.knots`, `.segments`, `.bp`, `.G`. Anything with `.knots` runs in the compiled kernel.

## Equilibria and stability

- `equilibria(g, Rt, rL)`: at a fixed point iL = (v1 − v2)/Rt and v2 = rL·iL, so `g(v1) = −v1/(Rt + rL)`, the load line through the origin. Roots by sign change on a fine grid plus `brentq`, with exact zeros caught separately (an exact zero gives no sign change) and near-duplicates merged. A test checks the three roots satisfy all three equations.
- `stability(g, v1, Rt, circuit)`: the Jacobian at the fixed point with the local slope G of g by central difference, and its eigenvalues.
- `hopf_point(g, circuit, side=-1)`: scans Rt from 1200 to 3200 ohm for the sign change of the complex pair's real part at the outer equilibrium on `side`, refines it, and returns (Rt, v1*, period). This is the DC to limit-cycle transition of the forward sweep. For the bench model it is called with `bench_static_g(P)` and `bench_small_signal(P)`.

## Integration, ideal model

- `rhs`, `rk4_step`: a vectorised numpy RK4 over a batch of runs, used only when the element has no `.knots` (an arbitrary callable); `collect_maxima` advances such a batch and keeps the maxima of v1.
- `integrate(g, circuit, rpot, y0, t_end, dt, keep=True, t_skip=0)`: a batch of runs at the given resistances from the states `y0` of shape (3, n). With `.knots` it dispatches to `kern.trajectory` (keep) or `kern.hold` (maxima only). A test compares it with the matrix exponential of a linear circuit, including the time grid after `t_skip`.
- `midpoint(v1, v2, rpot, r0)`: what CH3 would show, for the exported CSVs.

## The sweep protocol

- `sweep_starts(g, circuit, rpot)`: the two starting states. `forward` is the negative outer equilibrium nudged by 0.05 V in v1 (the origin where that side has no equilibrium); `back` is (7 V, 6 V, 0), inside the basin of the large outer cycle wherever it exists.
- `hold_batch`: settle then collect at fixed Rt for a batch, via `kern.hold` or the numpy loop.
- `sweep_order(rsw, direction)`: the order a sweep visits the resistances, descending for `down`, ascending for `up`.
- `sweep_continuation(g, circuit, rsw, dt, t_settle, t_collect, direction)`: visits the resistances in that order, starting on the corresponding start state at the first one (held `FIRST_SETTLE` = 0.2 s longer) and carrying the actual state from each resistance to the next, as the knob does. Returns the maxima per resistance in the order of `rsw`. A test with a mocked `hold_batch` checks the state really is carried and the output re-ordered. Defaults: the whole dial, 0–1000 ohm in 1 ohm steps, 20 ms settle and 100 ms of maxima per step.

## Bench model wrappers

- `bench_params(path, tau_a, tau_b=None, slew, tau_m, c1_law=True)`: the kernel's parameter vector from `identified.json`. The small-signal capacitance `identify.py` measures on node 1 is C1 + tau_A·AA/RA (the lag of stage A looks like a capacitor), so the kernel's C1 is that value minus the lag's share. With `c1_law` (and a law in the file) the capacitor follows the node-1 law (`c1_a0`, `c1_a1`, `c1_xk`) and stage B's lag is the law's `stage_B_lag_us`; without, C1 is constant and the lag is `TAU_B`. An explicit `tau_b` overrides either.
- `bench_static_g(P)`: the diode with infinitely fast op-amps, for equilibria and starting points.
- `bench_small_signal(P, r0)`: the linearisation as a `Circuit` (apparent C1, C2, L0, winding resistance r0). At rest there is no swing, so the node-1 law adds nothing here. With `bench_static_g` it is also the static model.
- `bench_sweep(P, rsw, dt, t_settle, t_collect, direction, r0)`: `sweep_continuation` for the bench kernel, from the same starting states (`kern.bench_state`), carrying the whole state (turning-point tracker included) from each resistance to the next.
- `bench_runs(P, rpots, states, t_end, dt, t_skip, r0)`: trajectories from given three-state points at each resistance.

## Measured records for comparison

- `find_records(folder)`: `{Rpot: (path, label)}` from legacy names containing "ohm" and from every `<sweep>_rpot.csv` sidecar whose record exists and whose fit is clean (status ok, residual ≤ 5 %, 0 to 1000 ohm). Sidecars are visited with `forward` first, and `setdefault` keeps the first entry, so the forward sweep wins ties.
- `nearest_record(records, rpot, tol=5, prefer='forward/')`: the record closest to the requested value within 5 ohm, preferring the forward sweep.
- `load_measured_bifurcation(folder)`: the `(rpot, max_v, sweep)` rows of every `*_bifurcation_points.csv` whose sweep folder exists, keeping only rows whose `sweep` column matches the file's own sweep and whose record still exists (checked once per record). A test checks that a points file without a sweep column and a foreign sweep row are ignored.

## Regime labels and transitions

- `maxima_period(m)`: the smallest lag p ≤ 64 at which the maxima repeat (mean |m_i − m_(i+p)| below 2 mV); 0 if none (aperiodic). The 2 mV floor makes the period-8 label 1–3 ohm late.
- `label(m, eqs)`: the regime at one resistance from its maxima and the equilibria there: rest, large cycle (a maximum above 5 V), double scroll (maxima beyond half of each outer equilibrium), origin orbit (a small orbit when the origin is the only equilibrium), `P<p>` or chaos.
- `transitions(rsw, max_down, max_up, g, circuit)`: the rows of the report's Table 2, the first resistance of each regime going down, the ends of the double scroll, and the last resistance of the unbroken large cycle going up; `regime_runs` lists the runs of equal labels. A test feeds a synthetic sweep and checks every reported value. `stacked_maxima` flattens a sweep for the `.npz`; `mean_period(t, v1)` is the mean spacing of maxima for the run table.

## Figures

- `plot_nr`: the model's element over the V–I trace (offset removed), breakpoints, and the load line −v/Rt for each `--r`.
- `plot_portraits`: one row per `--r`: the run from the negative equilibrium, the run from the large cycle, and the nearest measured record's v2 against v1.
- `plot_timeseries`: 10 ms of v1(t) for both starts per `--r`.
- `plot_bifurcation`: the measured points in grey underneath, the simulated down-sweep in blue and up-sweep in red on top.

## `main()`

1. `parse_args`, then `build_model`: the element, the `Circuit`, the bench parameters when needed and dt (0.5 µs ideal and static, 0.1 µs bench), with the components printed.
2. `report_equilibria`: the element's segments, the equilibria and their stability at each `--r`, and the Hopf point.
3. `portrait_runs`: `sweep_starts` gives both starts per `--r`; the states are interleaved (`y0[:, 0::2]`, `y0[:, 1::2]`) so one batch integrates everything; the run table lists ranges, maxima, distinct levels and period per start.
4. Names the measured record found for each `--r`, draws the three figures, and with `--export` writes each run as a scope-format CSV (`sim <R> ohm plus/minus.csv`, `export_runs`).
5. Unless `--no-sweep`, `bifurcation_sweep`: the continuation sweeps down and up (or the `seeded` protocol for the ideal model, which restarts every resistance from its seed), the bifurcation figure over the measured points, the regime runs, and the transitions, written to the `.npz` and `.json`.
