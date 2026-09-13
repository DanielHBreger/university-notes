# `tests/test_analysis.py`

Forty-one regression tests. Each checks one function against an independent equation or against deliberately bad data, so a change that breaks a result fails here before it reaches a figure. Run from the project root:

```bash
.venv\Scripts\python.exe -m pytest tests -q -p no:cacheprovider
```

`write_scope(path, t, *channels)` is the helper that writes a scope-format CSV from arrays.

## Reader and divider (`scope_data.py`, `rpot.py`, `sweeplib.py`)

- `test_divider_recovers_resistance_despite_offset_and_midpoint_gain`: a midpoint built by `simulate.midpoint` at 673 ohm, then scaled by 1.07 and offset by 0.09 V, is recovered to 1e-10. The midpoint channel's gain and offset cancel, as `rpot.md` claims.
- `test_divider_gain_ratio`: end channels with gains 1.02 and 0.97 are recovered when `g21 = 0.97/1.02` is given.
- `test_divider_rejects_unidentifiable_and_nonphysical_data` (flat, collinear, negative coefficient): each raises.
- `test_reader_rejects_broken_time_series` (NaN, duplicate timestamp, gap, reversed time): each raises.
- `test_reader_keeps_named_channel_order_and_float64_time`: asking for `('CH2','CH1')` returns them in that order, and 10 ns steps at t = 1 s survive in float64.
- `test_csv_discovery_and_sidecar_validation`: natural ordering (`trace2` before `trace10`), a summary CSV is excluded, `load_rpot` skips rows with NaN, `sibling` ignores a trailing slash, `folder_labels` disambiguates two `forward` folders.

## Maxima (`lorenz_map.py`)

- `test_peak_amplitude_period_and_return_time`: a 100 Hz sinusoid gives 20 maxima of amplitude 2, f0 within 2 %, and return times of 10 ms.
- `test_peaks_do_not_crash_for_short_flat_and_large_window`: three samples, a flat record, and a window larger than the record all return gracefully.

## Segment fit (`diode_analysis.py`)

- `test_piecewise_fit_recovers_known_slopes_without_duplicate_samples`: five known lines are recovered to 1e-8, every sample is assigned exactly once, R² = 1.
- `test_ranges_do_not_hide_invalid_slope_order`: an impossible slope order is reported `valid = False`; a valid side gives the right bounds.

## Lyapunov (`lyapunov.py`, `rosenstein.py`)

- `test_local_regression_handles_large_coordinate_offset`: data at x ≈ 10⁶ fits to 1e-7 (the centring in `local_fit`).
- `test_logistic_map_exponent_against_ln2`: the logistic map at r = 4 has lambda = ln 2; the return-map estimator gets within 0.08. A shuffled version is refused.
- `test_return_mean_not_median_and_nonmonotonic_times_rejected`: return times 1, 1, 4 repeated give mean 2 and median 1; a repeated peak time is refused.
- `test_exact_zero_derivative_not_silently_discarded`: a zero slope makes the mean −∞ rather than being dropped.
- `test_width_uncertainty_changes_window_even_when_floor_dominates`: the window-doubling systematic is nonzero.
- `test_nearest_recurrent_matches_exhaustive_search_even_k1`: the batched neighbour search equals brute force with the Theiler exclusion.
- `test_rosenstein_invalid_options`: six impossible settings raise.

## Simulation (`simulate.py`, `integration.py`, `benchmark_lyapunov.py`)

- `test_equilibria_include_exact_origin_and_satisfy_three_equations`: with the raw fitted offset the origin is an exact root; all three roots satisfy the circuit equations.
- `test_rk4_matches_exact_linear_circuit_and_time_grid`: `integrate` with a linear element against the matrix exponential, and the time grid starts at `t_skip`.
- `test_outer_segments_extrapolate_linearly_and_offset_is_removed`: the outer segments continue linearly to ±30 V and g(0) = 0.
- `test_simulation_skips_sidecars_without_a_sweep_column`: `load_measured_bifurcation` keeps only rows from the right sweep whose record exists.
- `test_continuation_carries_the_actual_previous_state`: with `_hold` mocked to add Rt to the state, the sweep visits the resistances in sorted order, carries the state, and returns results in input order.
- `test_compiled_kernel_matches_matrix_exponential`: `integration.trajectory` against the exact linear solution.
- `test_variational_exponent_agrees_with_linear_eigenvalues`: the tangent-vector exponent of `benchmark_lyapunov.simulate` at a stable fixed point equals the largest eigenvalue real part.
- `test_pwl_knots_are_continuous_with_the_given_slopes`: the knot slopes equal G and g(0) = 0.
- `test_kennedy_mapping_reproduces_the_segments`: `kennedy_from_pwl` then `static_pwl` returns the input breakpoints and slopes to 1e-9.
- `test_bench_kernel_matches_matrix_exponential_in_the_linear_limit`: the five-state kernel with rails far away and no Rayleigh terms against the exact linear solution.
- `test_bench_op_amps_saturate_and_hold_at_the_rails`: a railed output pushed outward has zero rate; pushed inward it moves at the slew limit.
- `test_regime_metrics_read_a_synthetic_sweep`: every transition `metrics` reports on a hand-built sweep.

## Identification (`identify.py`)

- `test_loop_integral_c1_and_flux_loop_recover_a_linear_circuit`: `loop_c1` recovers 11 nF to 1e-9 from a synthetic cycle; `inductor_loop` recovers the loss resistance to 0.2 %; `fit_tank` recovers C2, L, r to 2 % from three harmonics.
- `test_rayleigh_fit_and_cycle_markers`: `rayleigh_fit` inverts the secant relations exactly; `cycle_markers` spaces a sinusoid's crossings one period apart; `averaged_cycle` reproduces the waveform.
