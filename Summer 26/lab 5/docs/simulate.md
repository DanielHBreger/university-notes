# `chua/simulate.py`

Simulates the oscillator with three models of the same circuit and puts each beside the measurements: phase portraits next to the nearest record, a bifurcation sweep over the measured points, and the transition resistances.

Outputs beside the script, `--tag` adding a suffix: `simulated_nr<tag>.png` (the element and load lines), `simulated_portraits<tag>.png`, `simulated_timeseries<tag>.png`, `simulated_bifurcation<tag>.png`, optionally `simulated<tag>/` with scope-format CSVs. The printout is kept as `simulate<tag>.txt`.

## The three models

- `--model ideal`: the plan's model. Constant C1, C2, L; the element is the five segments fitted to the V–I record (`diode_fit.json`), with the fitted current offset removed and adjacent lines joined at their intersections. Defaults C2 = 100 nF and L = 18 mH nominal, C1 = 11.5 nF, the value that puts the Hopf point where it was measured (`--c1 10` for the nominal value, which puts it near 1000 ohm, beyond the dial).
- `--model static`: the identified circuit (`identified.json`) with everything held constant: small-signal C1, C2, L0, r0 and the identified five-segment law, integrated like the ideal model. It isolates what the right numbers do without the right physics.
- `--model bench`: the identified circuit in full: Rayleigh inductor, Kennedy diode with op-amp lag and slew limit, integrated with the five-state kernel. Two constants are not measured and are set here: `TAU_B = 2.4 µs` (stage B, gain 7.6, a 0.5 MHz gain-bandwidth op-amp) and `SLEW = 0.5 V/µs` (a 741-class part). `TAU_A = 0.05 µs` adds tau_A·AA/RA = 0.2 nF to node 1, which is subtracted from the measured C1 (see `bench_params`).

## Module constants and `set_components`

`C1, C2, L, R0` are module globals because `rhs`, `stability`, `integrate` and `_hold` read them. `set_components` overrides them from the command line (or from `feigenbaum.py` for the static model). Functions that must use other values (the bench model's linearisation) pass them explicitly through `c1=, c2=, l=`.

## The ideal element

- `load_vi_fit`: reads the segments from `diode_fit.json` and refuses a fit not made against the element voltage.
- `fitted_offset`: g(0) of the raw fit, the intercept of the segment containing v = 0. One scope code of channel offset.
- `pwl_from_segments(segments, i0)`: subtracts the offset from every intercept and puts each breakpoint where adjacent lines intersect, `(b2 − b1)/(m1 − m2)`. Raises if the intersections are out of order.
- `scale_inner_breakpoints`, `symmetrize_inner`: the `--bp-scale` and `--symmetric` experiments.
- `make_g(...)`: the ideal element as a callable; `pwl_g(bp, G)` builds the callable from knots (`integration.pwl_knots`) and attaches `.knots`, `.segments`, `.bp`, `.G`. Anything with `.knots` runs in the compiled kernel.

## Equilibria and stability

- `equilibria(g, Rt, rL)`: at a fixed point iL = (v1 − v2)/Rt and v2 = rL·iL, so `g(v1) = −v1/(Rt + rL)`, the load line through the origin. Roots by sign change on a fine grid plus `brentq`, with exact zeros caught separately (an exact zero gives no sign change) and near-duplicates merged. A test checks the three roots satisfy all three equations.
- `stability(g, v1, Rt, rL, ...)`: the Jacobian at the fixed point with the local slope G of g by central difference, and its eigenvalues.
- `hopf_point(g, rL, side, ...)`: scans Rt from 1200 to 3200 ohm for the sign change of the complex pair's real part at the outer equilibrium on `side`, refines it, and returns (Rt, v1*, period). This is the DC to limit-cycle transition of the forward sweep.

## Integration, ideal model

- `rhs`, `rk4_step`: a vectorised numpy RK4 over a batch of runs, used only when the element has no `.knots` (an arbitrary callable).
- `integrate(rpot, v1_0, t_end, dt, rL, g, keep=True, t_skip=0, y0=None)`: a batch of runs at the given resistances, from (v1_0, 0, 0) or from full states `y0` of shape (3, n). With `.knots` it dispatches to `kern.trajectory` (keep) or `kern.hold` (maxima only). A test compares it with the matrix exponential of a linear circuit, including the time grid after `t_skip`.
- `midpoint(v1, v2, rpot)`: what CH3 would show, for the exported CSVs.

## The sweep protocol

- `sweep_starts(g, rpot, rL)`: the two starting states. `forward` is the negative outer equilibrium nudged by 0.05 V in v1 (the origin where that side has no equilibrium); `back` is (7 V, 6 V, 0), inside the basin of the large outer cycle wherever it exists.
- `_hold`: settle then collect at fixed Rt for a batch, via `kern.hold` or the numpy loop.
- `sweep_continuation(g, rsw, rL, dt, t_settle, t_collect, direction)`: visits the resistances in order (descending for `down`, ascending for `up`), starting on the corresponding start state at the first one and carrying the actual state from each resistance to the next, as the knob does. Returns the maxima per resistance in the order of `rsw`. A test with mocked `_hold` checks the state really is carried and the output re-ordered.

## Bench model wrappers

- `bench_params(path, tau_a, tau_b, slew)`: the kernel's parameter vector from `identified.json`. The small-signal capacitance `identify.py` measures on node 1 is C1 + tau_A·AA/RA (the lag of stage A looks like a capacitor), so the kernel's C1 is that value minus the lag's share.
- `bench_static_g(P)`: the diode with infinitely fast op-amps, for equilibria and starting points.
- `bench_small_signal(P)`: (apparent C1, C2, L0, r0) of the linearisation. `bench_hopf` uses them.
- `bench_sweep`: `sweep_continuation` for the five-state kernel, starting from `kern.bench_state`.
- `bench_runs(P, rpots, states, t_end, dt, t_skip)`: trajectories from given three-state points at each resistance.

## Measured records for comparison

- `find_records(folder)`: `{Rpot: (path, label)}` from legacy names containing "ohm" and from every `<sweep>_rpot.csv` sidecar whose record exists and whose fit is clean (status ok, residual ≤ 5 %, 0 to 1000 ohm). Sidecars are visited with `forward` first, and `setdefault` keeps the first entry, so the forward sweep wins ties.
- `nearest_record(records, rpot, tol=5, prefer='forward/')`: the record closest to the requested value within 5 ohm, preferring the forward sweep.
- `load_measured_bifurcation(folder)`: the `(rpot, max_v, sweep)` rows of every `*_bifurcation_points.csv` whose sweep folder exists, keeping only rows whose `sweep` column matches the file's own sweep and whose record still exists. A test checks that a points file without a sweep column and a foreign sweep row are ignored.

## Regime metrics

- `describe(m)`: a one-line label for a list of maxima: large cycle (max > 5 V), double scroll (a maximum above 0.9 V and one below −0.3 V), or the number of distinct levels and the range.
- `metrics(rsw, max_down, max_up, persist=4)`: the transition resistances. Each regime predicate (`is_lc`, `is_ds`, `is_p1`, `is_split`, `is_chaos`) is applied per resistance, and a transition is reported at the first resistance where the predicate holds for `persist` consecutive steps in sweep order, so the slow transients just below the Hopf point do not count. `r1` requires a period-1 stretch to have been seen first. A test feeds a synthetic sweep and checks every reported value. `period_of(t, v1)` is the mean spacing of maxima for the run table.

## Figures

- `plot_nr`: the model's element over the V–I trace (offset removed), breakpoints, and the load line −v/Rt for each `--r`.
- `plot_portraits`: one row per `--r`: the run from the negative equilibrium, the run from the large cycle, and the nearest measured record's v2 against v1.
- `plot_timeseries`: 10 ms of v1(t) for both starts per `--r`.
- `plot_bifurcation`: the measured points in grey underneath, the simulated down-sweep in blue and up-sweep in red on top.

## `main()`

1. Chooses dt (0.5 µs ideal and static, 0.1 µs bench), sets the components, builds the element and prints its segments.
2. Prints the equilibria and their stability at each `--r`, and the Hopf point (from the bench linearisation for the bench model).
3. Portrait runs: `sweep_starts` gives both starts per `--r`; the states are interleaved (`y0[:, 0::2]`, `y0[:, 1::2]`) so one batch integrates everything; the run table lists ranges, maxima, distinct levels and period per start.
4. Names the measured record found for each `--r`, draws the three figures, and with `--export` writes each run as a scope-format CSV (`sim <R> ohm plus/minus.csv`).
5. Unless `--no-sweep`: the continuation sweeps down and up (or the `seeded` protocol for the ideal model, which restarts every resistance from its seed), the bifurcation figure over the measured points, a regime table every 25 ohm, and the transitions from `metrics`.
