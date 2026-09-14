# Experiment plan v5 against what was done

Status key: **done** (the item is complete and its output is in the folder), **done, caveat** (complete, with something the reader must know), **partly** (some of the item is missing), **not verifiable** (a bench-time instruction that leaves no trace in the records), **not done**.

Sources: `Lab5_Chua_experiment_plan.pdf` (v5), `Fig1._NR_scheme.pdf`, `Fig2_Chua_Oscillator.png`; results in `RESULTS.md`, pipeline in `README.md`, code explanations in `docs/`.

## The setup

| plan | status | description | notes |
|---|---|---|---|
| Component values: L 18 mH, C2 100 nF, C1 10 nF, R0 990 ohm measured, R 0-1 kohm ten-turn, Rs 216 ohm | done, caveat | The nominal values are the ideal model's defaults (`simulate.py`) and Rs = 216 is in `diode_analysis.py`. R0 is 992 ohm in the code (`scope_data.R0`), the value measured on this bench; the plan's 990 would scale every Rpot by 0.2 %. | The oscillator records say the components are not nominal: C1 11.1 nF, C2 90.5 nF, L 19.1 mH rising to 24 mH with current (`identify.py`, RESULTS "The circuit as measured"). Every model comparison in RESULTS is done both ways. |
| Element as in Fig. 1 (Kennedy two-op-amp diode, TL082, R1 = R2 = 220, R3 = 2.2 k, R4 = R5 = 22 k, R6 = 3.3 k) | done | The identified element maps onto exactly this circuit: AA 1.115 (schematic 1.100), AB 7.62 (7.667), Ga -0.762 mS (-0.758), Gb -0.415 (-0.409), rails 6.7-7.7 V on 9 V. Table in RESULTS. | Gc is 3.9-4.2 mS against the schematic's 4.6 (RA 249 against 220 ohm): the saturated op-amp output adds some tens of ohms. The op-amp constants the bench model uses (2.4 us lag, 0.5 V/us) are not a TL082's small-signal figures; they stand in for its recovery from saturation. See RESULTS. |
| Probe points 1 = V1 (CH1), 2 = V2 (CH2), 3 = V3 (CH3) | done, caveat | Every script uses this mapping, and the records confirm it (CH1 swings several volts, CH2 under 2 V, divider inside 0-1000 ohm). | `Fig2_Chua_Oscillator.png` labels the channels the other way round (node 2 -> CH1, node 1 -> CH2). The plan's own Figure 1 and the data agree with the code; the drawing is wrong. Noted in RESULTS. |
| Oscilloscope DSO-X 3014A | done | `s3014a.py` is the capture program; 1 Mpt records. | |

## Reading R without opening the circuit

| plan | status | description | notes |
|---|---|---|---|
| Fit V3 = A V1 + B V2 + C, R = R0 B/A; residual is the health check (1 % normal, 5 % bad) | done | `rpot.py`, `batch_rpot.py`; `forward_rpot.csv` (405 of 418 clean), `back_rpot.csv` (74 of 80). The 5 % cut is `sweeplib.CLEAN_DIVIDER_MAX_PCT`. | The fit is done on scaled, centred regressors with a rank test, so it refuses records that cannot identify the ratio instead of returning a number (`docs/rpot.md`). |
| Absolute accuracy a few percent from the relative channel gain; relative accuracy a fraction of an ohm | done | Checked directly: the identity A + B = 1 holds to 1.003 +- 0.004 over all records (`identify.py`), so the two end channels have equal gain to 0.5 %. | The back sweep reads about 6 ohm (0.8 %) lower than the forward sweep; forward values below 467 ohm carry a +2 % step (RESULTS data-quality notes). |

## Must have

| plan | status | description | notes |
|---|---|---|---|
| **M1** V-I curve: K1 open, K2 = VI, triangular sweep 20-100 Hz past +-7 V; v = VI, i = (VV - VI)/Rs | done, caveat | `trace1.csv`, `find_breakpoints.py` -> `diode_fit.json`, `diode_segments.csv`, `current_vs_voltage.png`, `optimized_breakpoints.png`. Both saturation segments captured (element reaches -8.6 and +7.7 V). | The drive was 9.2 Hz, below the plan's 20-100 Hz. Harmless: slower is what the plan asks for if a loop appears. |
| M1: check that forward and backward sweeps overlap | done | Rising and falling sweeps agree to 0.027 mA rms, far below the 0.6 mA current code; no loop. Lower panel of `current_vs_voltage.png`. | |
| M1: fit Ga, Gb left and right separately, Gc, breakpoints | done | Five segments by a global dynamic-programming fit over the voltage codes (`fit_segments`); Gb left -0.4176, right -0.4337 mS, reported separately throughout. | The inner segment spans only three current codes (R2 0.59): Ga = -0.717 mS is the least certain number of the trace. The oscillator itself gives -0.762 mS, and the schematic predicts -0.758. |
| M1 trap: boundaries at the found breakpoints, not by eye | done | Breakpoints are the optimiser's, placed between adjacent codes. | |
| M1 trap: intercept of the inner segment | done | -0.413 mA = 0.7 of one code of channel offset; removed in the models (`fitted_offset`), kept in the JSON. | |
| **M2** 1/abs(Ga) < Rtotal < 1/abs(Gb), converted to dial readings, written before looking | done, caveat | `slope_ranges` in `diode_analysis.py`: Rpot > 403 ohm, Rpot < 1314 (right) / 1403 (left) ohm. The observed double scroll (668-328) lies inside. | The lower bound from the trace (403) is above the observed lower end (328); with the inner slope the circuit sees (-0.762 mS) the bound is 320. Whether the prediction was written down before the sweep is not recorded. |
| **M3** find the double scroll: K2 = CHAOS, knob turned slowly, record on every change of character, especially at the two edges of the double scroll; 200 ms records | done | 418 forward and 80 back records, 200 or 500 ms, 570-1400 windings each. Double scroll 668-328 ohm (forward); edges resolved to 1-2 ohm. RESULTS M3 table. | Records were taken every 1-2 ohm rather than only on changes, which is what makes N3-N5 possible. |
| M3: run rpot.py, sort by R, compare with M2 | done | `batch_rpot.py`, `bifurcation.py`; the range lies inside M2. | |
| M3: at least three long clean records well inside the range | done | forward/trace342 (450 ohm), trace270 (549), trace229 (650): 1170-1300 maxima each, two-branch maps with within-branch R2 above 0.9. RESULTS M5 table. | |
| **M4** Lorenz map: Savitzky-Golay at 1/20 period, 2 % prominence, M(n+1) against M(n), diagonal drawn | done | `lorenz_map.py` -> `forward_lorenz/` (grid pages, pooled map, pairs and summary CSVs), `back_lorenz/`; the three picks have their maps with the local fit in `forward_lyapunov_each/`. | The period is measured per record from the peak spacing (the time base changed mid-sweep), the prominence floor is two quantisation steps, and duplicate flat-top peaks are merged (`docs/lorenz_map.md`). |
| M4: curve rather than cloud; repeat at two or three R and show how the map changes | done | The double-scroll maps are curves (R2 0.9 within branches); the sequence period 1 -> 2 -> 4 -> 8 -> single scroll -> double scroll is the grid pages. | |
| **M5** lambda = <ln abs(f')> / <T> from a sliding linear fit | done | `lyapunov.py` -> `forward_lyapunov.png/.csv`, `back_*`; 286 of 405 forward records get a map exponent; the three picks 2360-3080 /s. | The window is a fixed width in M(n), the fit is judged within branches, and periodic records are refused with a reason (`docs/lyapunov.md`). |
| M5 trap: mean return time, not median (difference above 10 %) | done | Both reported per record. In the double scroll the mean exceeds the median by 8 % (median over records), by more than 10 % in 40 % of records; in the single-scroll band they agree to 1 %. | The plan's "about twice as long" per lobe switch is consistent: the tail is there, the average excess is smaller than the plan's example. |
| M5 trap: average over the points visited, including the turning point | done | Every visited point enters; slopes below their standard error are floored at se/e and the unresolved fraction is reported. | |
| M5: quote as a return-map estimate; Rosenstein typically 10-20 % lower | done, caveat | `rosenstein.py` runs on every admitted record: direct exponents 1500-2000 /s on the three picks, map/direct median ratio 1.32 in the double scroll. `benchmark_lyapunov.py` calibrates both on a simulated record: map 1.05x, direct 0.92x the truth at this bench's parameters; 1.42x and 1.09x on Matsumoto's scroll. | The measured gap (map above direct by 30 %) is larger than the plan's 10-20 %; the benchmark shows why (the map folds at the lobe switch). RESULTS quotes the direct value as the measurement. |

## Nice to have

| plan | status | description | notes |
|---|---|---|---|
| **N1** Hopf point R_Hopf = -Gb(1 + C1/C2)/(Gb^2 + C1^2/(C2 L)); measure the onset and the frequency there | done | Plan formula with nominal C1 = 10 nF: 966-1006 ohm, beyond the dial. Measured onset 893 ohm (forward scale; squared amplitude extrapolated to zero), period 330.3 us. Identified circuit: 904 ohm, 330 us. | The plan's formula and the exact eigenvalue calculation (`simulate.hopf_point`) agree to the ohm; both are used. The mismatch with nominal values is C1 (11.1 measured against 10) and C2 (90.5 against 100), not the formula. The models' onset is abrupt where the bench's amplitude grows over 50 ohm (RESULTS N1). |
| **N2** gallery: limit cycle, period 2, single scroll, double scroll | done | `gallery.py` -> `chua/gallery.png`: trace54 (804 ohm), trace82 (760), trace164 (700), trace270 (549), plus the large cycle trace411 (191). Same records beside the simulations in `simulated_portraits_*.png`. | Added on 13 September; before that the four records existed but only inside the simulation figures. |
| **N3** bifurcation diagram, maxima of V1 against R, as a scatter | done | `bifurcation.py` -> `forward_bifurcation.png` + `_points.csv`, `back_*`. 405 columns forward, 74 back. Clipped records dropped (none over 2 %). | The sampling is 1-2 ohm, fine enough that the cascade is resolved to period 8. |
| **N4** sweep both ways without opening the circuit, overlay; coexisting attractors; be honest about sparse sampling | done | `forward_vs_back_hysteresis.png` + `_points.csv`. Double scroll survives down to 328 ohm; the large cycle survives up to 675 ohm on the way back. The bench model reproduces both (344, 683). | The back sweep is dense only between 685 and 780 ohm; below that it holds the large cycle, so the double-scroll edge on the way up is not sampled and is stated as such. |
| **N5** period doubling and Feigenbaum; needs a fine even step | done, caveat | `cascade_periods.py` (lag distances, sqrt law): R1 774.6 +- 1.5, R2 755.6 +- 0.8, R3 750.5 +- 0.4 ohm, delta_1 = 3.7 +- 0.7 (universal 4.669). `feigenbaum.py` for the models: bench model R1 776.8, R2 754.9; plan's model 771.9, 732.5. | delta_1 from the first three doublings is not expected to be 4.669 (convergence is geometric); one ratio with 20 % error is what the data supports. No model resolves period 8: in all three the period-4 orbit gives way to chaos within 0.05 ohm. |
| **N6** eigenvalues at all three equilibria from the measured slopes; check abs(sigma) < abs(gamma) | done | `shilnikov.py` -> `chua/shilnikov.txt`, table in RESULTS N6. Outer equilibria: ratio 0.02-0.04 everywhere. Origin: 0.25 at 806 ohm rising to 0.69 at 400 and 1.41 at 328, where the condition fails; the double scroll ends at 328. Pair frequency 3.0-3.2 kHz = the winding frequency. | Added on 13 September. |

## Oscilloscope settings (check at the start of every session)

| plan | status | description | notes |
|---|---|---|---|
| 1. Channel mapping CH1 = V1, CH2 = V2, CH3 = V3 | done | Verified from the records by the swing test the plan gives and by the divider result. | The supplied Fig. 2 drawing has CH1 and CH2 swapped; see above. |
| 2. DC coupling on all three channels | done | The single-scroll records sit on one side of zero (V1 -3.4 .. +0.2 V at 700 ohm), which AC coupling would have erased. | |
| 3. 20 MHz bandwidth limit on | not verifiable | Leaves no trace in a 2.85 kHz record. | |
| 4. Fix V/div once and leave it | not done, caveat | The vertical range was changed several times within both sweeps, always on all three channels together. The divider identity A + B = 1 holds on every range, so the gain ratio survived; a +2 % step in Rpot appears between forward trace329 and trace330 (RESULTS data-quality notes). | The cascade region is unaffected. |
| 5. Probe attenuation matches the probe switch | not verifiable | | A mismatch would be a factor 10, which the amplitudes rule out. |
| 6. x10 probes, compensated | not verifiable | | A x1 probe's 90 pF is 0.9 % of C1; the measured C1 excess is 11 %, so probes are not the explanation. |
| 7. Acquisition Normal | not verifiable | | The quantisation steps (43 / 201 mV) are those of normal acquisition; averaging would have produced finer codes. |
| 8. 20 ms per division, 200 ms record | done | Records are 200 or 500 ms. | |

## Time budget and the two things that can ruin it

| plan | status | description | notes |
|---|---|---|---|
| Meeting A: setup, M1, fit, M2 prediction, first hunt | done | M1/M2 outputs exist. | Which meeting produced what is in the git history of branch `lab-5` (commits "meeting 3 set 3+4", "set 5"), not in this checklist. |
| Meeting B: both boundaries, long records, start M4 | done | | |
| Meeting C: M4, M5, then N1, N2 | done | All of N1-N6 done. | |
| Wrong mapping / AC coupling / changed V/div | done, caveat | Mapping and coupling are right; V/div did change, with the consequence quantified above. | |
| Records too short | done | 570-1400 windings per record. | |

## Beyond the plan

| item | where | why it is there |
|---|---|---|
| Identification of the actual circuit from the records (C1, C2, Rayleigh inductor, Kennedy diode with rails) | `identify.py`, `identified.*`, RESULTS "The circuit as measured" | The plan's model does not reproduce the sweeps with any constant components; the records show why. |
| Bench model with the identified non-idealities, reproducing every transition to 1-4 % | `simulate.py --model bench`, `simulated_*_bench.png`, RESULTS "Models against the bench" | Closes the gap between simulation and measurement. Two op-amp constants are set, not measured (see the Fig. 1 row above). |
| Calibration of both Lyapunov estimators on a simulated record | `benchmark_lyapunov.py`, `benchmark_lyapunov.txt` | Puts a number on the plan's "10-20 % lower". |
| 42 regression tests | `tests/test_analysis.py` | Every estimator against an independent equation or bad data. |
| Function-by-function explanations | `docs/` | |

## Open points

- The bench model's op-amp constants (2.4 us, 0.5 V/us) reproduce the data but are not a TL082's small-signal figures (0.4 us at gain 7.6, 13 V/us). The mechanism is recovery from output saturation; modelling it as a delay rather than a slew limit is the next refinement and would likely also fix the 7 % short large-cycle period. Not done.
- `Fig2_Chua_Oscillator.png` should have its CH1/CH2 labels corrected before it goes into a report.
- The period-8 window is resolved in the bench records but in none of the models at the default settle time; a longer `--settle` in `feigenbaum.py` might resolve it. Not tried beyond the defaults.
