# Chua oscillator, lab 5

Analysis of the bench measurements against `Lab5_Chua_experiment_plan` (Summer 2026, v5),
branch `lab-5-review`. Results and the model/measurement comparison are in `RESULTS.md`.

## Data

- `chua/forward/`: the sweep with the knob turned down, 418 records, 905 -> 4 ohm analysed
  (the 13 flagged records at the top of the dial are excluded).
  `chua/back/`: the sweep back up, 80 records, 3 -> 900 ohm analysed (6 flagged records
  excluded), dense between 685 and 780 ohm.
  These two folders and the V-I record are all the report uses; earlier sessions remain in
  the history of branch `lab-5`.
- Oscillator CSVs: `Time(s),CH1(V),CH2(V),CH3(V)` with CH1 = V1 (across the element),
  CH2 = V2, CH3 = V3 (node between R0 and the potentiometer). Records are 200 or 500 ms.
- `trace1.csv` at the root is the V-I record: CH1 = VV (source side of the 216 ohm shunt),
  CH2 = VI (across the element); current = (CH1 - CH2) / 216.

## Pipeline

Run with `.venv\Scripts\python.exe` (or any Python 3.11+ with `requirements.txt`; numba is
needed for the simulations). The first script runs from this folder, the rest from `chua/`:

```
python find_breakpoints.py                       # M1/M2 -> diode_fit.json, diode_segments.csv, current_vs_voltage.png, optimized_breakpoints.png
cd chua
python batch_rpot.py forward                     # M3 -> forward_rpot.csv        (same for back)
python lorenz_map.py forward --pooled            # M4 -> forward_lorenz/          (same for back)
python bifurcation.py forward                    # N3 -> forward_bifurcation.png + _points.csv   (same for back)
python bifurcation.py forward back -o forward_vs_back_hysteresis.png   # N4 overlay
python lyapunov.py forward --rosenstein --each   # M5 -> forward_lyapunov.png/.csv + forward_lyapunov_each/   (same for back)
python cascade_periods.py forward                # N5 from the bench -> cascade_periods_forward.txt   (same for back)
python delta_bench.py                            # N5: R1-R3 (brackets) and Feigenbaum delta_1 -> delta_bench.txt/.json
python period8_check.py                          # which cascade records are a steady period 8 -> period8_check.txt/.png
python onset.py                                  # N1 from the bench -> onset.txt, onset_records.csv, onset.pdf/.png
python report_figures.py                         # report figures for 2.2-2.5 -> ../report_figures/bifurcation, gallery, cascade, return_maps, lyapunov, simulation (.pdf); run after simulate.py
python gallery.py                                # N2 -> gallery.png (limit cycle, period 2, single and double scroll, large cycle)
python shilnikov.py                              # N6 -> shilnikov.txt (eigenvalues at the three equilibria, |sigma|/|gamma|)
python identify.py                               # the circuit from the records -> identified.json/.txt + identified_*.png (10 min)
python simulate.py --tag _nominal                # the plan's model: nominal parts, Table 1 element -> simulated_*_nominal.*, simulate_nominal.txt
python simulate.py --model static --tag _static  # identified circuit, constant components -> simulated_*_static.png, simulate_static.txt
python simulate.py --model bench --tag _bench    # identified circuit, datasheet op-amps -> simulated_*_bench.*, simulate_bench.txt
python feigenbaum.py                             # N5 from the plan's model -> feigenbaum_ideal.txt   (same with --model static, --model bench; not rerun since 27 Sep)
python verify_models.py                          # checks behind the simulation (kernels, Hopf, node-1 C1) -> verify_models.txt
python benchmark_lyapunov.py                     # estimator calibration -> benchmark_lyapunov.txt
cd .. ; python -m pytest -q -p no:cacheprovider --basetemp=<a writable temp dir> tests
```

## What the scripts do

Function-by-function explanations of every script, with the naming conventions and the
flow between them, are in `docs/` (start at `docs/README.md`).

- `scope_data.py`: strict reader of the scope CSVs; `R0 = 992 ohm`, the value measured on this
  bench (the plan quotes 990; the difference is a 0.2 % scale on every Rpot).
- `rpot.py`, `batch_rpot.py`: the divider fit V3 = A V1 + B V2 + C, Rpot = R0 B/A. A record is
  dropped when the fit residual exceeds 5 % or Rpot falls outside 0-1000 ohm.
- `lorenz_map.py`, `bifurcation.py`, `lyapunov.py`, `rosenstein.py`: maxima of V1
  (Savitzky-Golay at period/20, 2 % prominence, quarter-period spacing), return maps, the
  bifurcation diagram, the return-map and Rosenstein exponents. `cascade_periods.py` reads the
  period of each record from lag distances and locates R1, R2, R3 by the sqrt law.
- `gallery.py`: the measured phase portraits of the regimes the plan names (N2), nearest clean
  forward record to each resistance. `shilnikov.py`: eigenvalues at the three equilibria from
  the measured slopes and the Shilnikov ratio |sigma|/|gamma| (N6), for the plan's model and
  the identified circuit.
- `diode_analysis.py` (`find_breakpoints.py`, `plot_current.py`): M1, the five segments fitted
  against the element voltage, with the sweep-direction check and the M2 range.
- `identify.py`: recovers the circuit from the oscillator records themselves. Periodic records
  are averaged over their cycles; the loop integral of node 1 gives C1 for any static element,
  the tank admittance gives C2, the flux-current loop gives the inductor's inductance and loss
  against amplitude (a Rayleigh law), and integrated KCL over every record gives the five-segment
  law as the circuit sees it, which maps onto Kennedy's two-op-amp diode.
- `integration.py`: numba RK4 kernels for both models. `simulate.py`: `--model ideal` is the
  plan's circuit (nominal components, Table 1's element as a continuous law through the origin);
  `--model static` is the identified circuit with everything held constant; `--model bench` is
  the identified circuit in full (Rayleigh inductor acting on the current's excursion from its
  mean, two-op-amp diode with TL082 datasheet speeds). All three sweep the whole dial as a
  continuation, down from the negative equilibrium and up from the large cycle, as the knob
  does, label the regime at each resistance and write the Table 2 transitions
  (`simulated_transitions_<tag>.json`); they draw the nearest forward-sweep record beside each
  simulated portrait. The checks behind this (27 September) are in RESULTS.md.
  `simulate.txt`, `simulated_*.png` without a tag and `*_sym*` are from the pre-27-September
  models (C1 = 11.5 nF, line-intersection element) and are superseded.
- `onset.py`: N1, the Hopf onset, on rpot.py's resistance scale. rpot.py's divider fit reads
  small cycles low (noise on V1, V2 attenuates the fitted coefficients), so the accepted records
  nearest the onset are placed by their mean voltages,
  calibrated against rpot.py on large cycles of the same V/div range; V1's fundamental amplitude
  comes from a sine fit, and A1^2 is extrapolated to zero. rpot.py itself is left as handed out.
- `feigenbaum.py`: the period-doubling points of either model by bisection on the orbit period
  and the sqrt law of the splitting. `benchmark_lyapunov.py`: the two exponent estimators run
  on a simulated double scroll with known exponent.
- `chua/exploration/`: scratch from the 11 September search for the mismatch, superseded by
  `identify.py`; not part of the pipeline.

## Phase-space animations

`python chua/animate_resistance.py` creates `measurements_resistance_sweep.gif`
and `simulation_resistance_sweep.gif` in `chua/animations/`. Each frame shows a
20 ms portrait at a different Rpot, decreasing across the clean forward sweep.
Both use the same measured resistance values, fixed full-view and zoom axes,
and an Rpot indicator. The ideal simulation carries its final state to the next
resistance, settling for 100 ms at each value. Its model is the presentation's:
the fitted V-I lines joined at their intersections with C1 = 11.5 nF
(`animate_resistance.ANIMATION_CIRCUIT`), not the report's `simulate.py --model ideal`. Measured portraits are selected
near an evenly spaced resistance grid with extra frames around 700–810 ohm;
they are not interpolated. Parameters and source records are saved in
`resistance_sweep_metadata.json`.

The same script also creates `comparison_equal_r.gif`, placing measurements and
the ideal simulation side by side at identical Rpot, with matched full-view and
zoom axes. Use `python chua/animate_resistance.py --comparison-only` to render
only this combined animation.

`python chua/animate_resistance.py --model bench --comparison-only` creates
`comparison_equal_r_bench.gif` with the identified circuit described in
`RESULTS.md`: measured component values, the Rayleigh inductor, and Kennedy diode
with amplifier lag and slew. It uses `identified.json` and the TL082 datasheet
amplifier parameters from `simulate.py`, a 0.1 µs integration step, and carries all
six states between resistances. Exact parameters are saved in
`resistance_sweep_bench_metadata.json`.

`python chua/animate_phase.py --r 550 --window-ms 8` creates two looping GIFs in
`chua/animations/`: the nearest clean forward measurement and the animations' ideal
model at that record's fitted resistance. Both show V2 against V1 with common axis limits,
a growing path, a 0.25 ms highlighted trail, and elapsed time. The displayed 8 ms
plays over 12 seconds; the simulation discards 100 ms of settling first. The two
windows are not phase aligned. `metadata.json` records the source and parameters.

## Conventions and data-quality notes

- Sidecars are per sweep folder (`forward_*.csv`, `back_*.csv`); the hysteresis overlay is
  named `forward_vs_back_hysteresis*` so that no script mistakes it for a sweep.
- The vertical range was changed several times within both sweeps, always on all three
  channels together, so the divider ratio is unaffected (the divider identity A + B = 1 holds
  with mean 1.0029 and maximum deviation 1.25%, `identified.txt`; this does not calibrate the channel gain ratio). Between forward trace329 and trace330 the fitted
  Rpot steps by 11 ohm in a sweep that otherwise falls 1-2 ohm per record, so forward values
  below about 467 ohm carry a +2 % offset relative to those above.
- The back sweep's resistance scale reads about 6 ohm lower than the forward sweep's near
  the cascade (R1 768 against 775), a 0.8 % difference between the two sweeps (both were
  taken in one sitting; the cause is not established).
- Records flagged by the divider (13 forward, 6 back) are all at the top of the dial: records
  where the circuit rests at its equilibrium (forward trace1-11, back trace76-80), which give
  the divider fit nothing to use, and the three smallest limit cycles (forward trace12-13, back
  trace75), whose fitted resistance is biased low by 70-140 ohm. No record was taken with the
  knob moving. They are excluded everywhere; `onset.py` places the slightly biased accepted
  cycles below them by their mean voltages.
- The resistance scale shifts by 0.5-2 % between vertical ranges (RESULTS data-quality notes).

## History

`lab-5` (commit c75206a) is the analysis as it stood after the third meeting; `lab-5-chatgpt-raw`
preserves a ChatGPT revision of it, of which the strict CSV reader, the M1 fit against the
element voltage, the input validation and the tests were kept and the rest rejected
(its `simulate.py` never oscillated above 850 ohm because it kept the fitted current offset,
and its conclusions were drawn from data sets no longer in use).

`chua/rpot.py` was handed out with the lab (commit f5ad118) and must be used as is. The copy
here is the ChatGPT rewrite of it (strict CSV reader, input validation, errors instead of
numbers for non-physical fits). Checked on 25 September: on every one of the 498 records where
the current file returns a resistance, the handed-out file returns the same value to 0.005 ohm;
on the other 8 (the DC records forward trace1-3, back trace76-80) it returns -2 to -20 kohm at
11-14 % residual. Do not modify it; corrections belong in separate scripts (`onset.py`).

## Measurement uncertainty

The instrument audit and propagation equations are in [uncertainty/REPORT.md](uncertainty/REPORT.md).
Run `python uncertainty/calculate.py` to rebuild the voltage-resolution, divider-fit, diode-fit,
cascade and power-supply calculations, then regenerate the presentation figures.
Settings and calibration inputs belong in `uncertainty/instrument_inputs.json`; missing values
are explicit and are never treated as zero error. `python -m pytest tests/test_uncertainty.py`
checks covariance propagation and paired block resampling.
