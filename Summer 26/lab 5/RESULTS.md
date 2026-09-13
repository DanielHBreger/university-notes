# Results (rerun of 12 September 2026)

Inputs: `chua/forward/` (418 records, knob turned down), `chua/back/` (80 records, knob turned
back up) and the V-I record `trace1.csv`. Divider fits use R0 = 992 ohm (`chua/scope_data.py`).
Every number below is reproduced by the pipeline in `README.md`.

## Summary

- The measurements themselves are sound: the divider identity A + B = 1 holds to 0.5 % on every
  record, so the resistance axis is right to that level, and the whole regime sequence
  (Hopf point, period doubling, single scroll, double scroll, large outer cycle, hysteresis) is
  in the data with 570 to 1400 windings per record.
- The plan's model (constant components, five straight segments from the V-I trace) does not
  reproduce the measurement, and no choice of its constants makes it: the components are not
  what the model assumes. Identified from the oscillator records themselves (`identify.py`):
  the inductor is a ferromagnetic-core choke whose inductance grows from 19 to 24 mH and whose
  loss grows from 2 to 57 ohm between the smallest cycle and the large cycle (a textbook
  Rayleigh law, seen directly as lens-shaped flux-current loops), C2 is 90 nF rather than 100,
  and the element is Kennedy's two-op-amp diode whose op-amps are slow enough to be seen
  (0.7 mA of dynamic deviation in the large cycle).
- With those measured non-idealities put into the simulation (`simulate.py --model bench`)
  every transition of the bench is reproduced to within 1 to 4 % in resistance: Hopf point
  905 against 893 to 897 ohm, first doubling 776 against 775, chaos 750 against 750, double
  scroll 697 to 344 against 668 to 328, large cycle surviving to 683 against 675 ohm on the
  way up, with the periods, amplitudes and the tank current of every regime within 2 to 7 %.
  What remains is the large cycle's period (7 % short) and the abrupt onset of the small
  cycle, both from the piecewise-linear treatment of the op-amp knees.

## M1  V-I characteristic  (`find_breakpoints.py` -> `diode_fit.json`)

| Segment | VI range (V) | slope (mS) | intercept (mA) |
|---|---|---|---|
| Gc left | -8.62 .. -6.75 | +3.623 | +27.14 |
| Gb left | -6.75 .. -1.09 | -0.4176 | -0.035 |
| Ga inner | -1.09 .. +0.58 | -0.7171 | -0.413 |
| Gb right | +0.58 .. +5.98 | -0.4337 | -0.752 |
| Gc right | +5.98 .. +7.72 | +3.786 | -26.06 |

Combined R2 0.964. Rising and falling sweeps agree to 0.027 mA rms, far below the 0.6 mA current
code, so there is no loop (drive 9.2 Hz). The inner intercept of -0.413 mA equals 89 mV across the
two channels, 0.7 of one 129 mV code: a channel offset, which the models remove. The inner
segment spans only three current codes (R2 0.59), so its slope and breakpoints are the least
certain numbers of the trace; the circuit itself gives them ten times more precisely (below).

## M2  Predicted range for three equilibria

1/|Ga| < R0 + Rpot < 1/|Gb| gives Rpot > 403 ohm and Rpot < 1314 (right shoulder) / 1403 (left)
ohm, i.e. the upper bound lies beyond the dial. With the inner slope the circuit sees
(-0.762 mS) the lower bound is 320 ohm, and the double scroll is observed down to 328 ohm.

## M3  Sweeps  (`batch_rpot.py`, `bifurcation.py`)

| | forward (down) | back (up) |
|---|---|---|
| usable records (residual < 5 %, 0-1000 ohm) | 405 of 418 | 74 of 80 |
| oscillation seen up to | 891 ohm, DC above | 873 ohm, DC at 901 |
| period-1 | 891 .. 775 | 873 .. 768 |
| period-2 from | 774.6 | 768.1 |
| period-4 from | 755.6 | 744.7 (one record) |
| period-8 from | 750.5 | - |
| single-scroll chaos | 750 .. 722 (periodic window 719-713) | 737 .. 717 |
| single scroll, no lobe change | 721 .. 670 | 713 .. 685 |
| double scroll | 668 .. 328 | (not reached from below: the large cycle holds) |
| small orbit around the origin | 328 .. 315 (period 660-710 us) | - |
| large outer cycle | 314 .. 4 | 3 .. 675 |

The observed double-scroll range sits inside the M2 range. The back sweep's resistance scale
reads about 6 ohm lower than the forward sweep's (R1 768 against 775, onset 874 against 880),
a 0.8 % gain difference between the two sessions.

## N1  Hopf point

| | value |
|---|---|
| measured onset | squared amplitude of the small cycles extrapolates to zero at 893 ohm (forward scale) |
| measured period at onset | 330.3 us (18 records between 875 and 891 ohm); 327-329 us at 891-901 |
| plan formula, nominal C1 = 10 nF, C2 = 100 nF, L = 18 mH, trace slopes | 968-1008 ohm: beyond the dial |
| plan formula, C1 = 11.5 nF | 859-889 ohm |
| identified circuit (C1 11.1 nF, C2 90.5 nF, L0 19.1 mH, r0 2.3 ohm) | 904 ohm, 330 us |

The amplitude grows over 50 ohm below the onset (peak-to-peak V1 0.8 V at 890, 1.5 V at 871,
3.8 V at 845 ohm); the piecewise-linear models reach full size within a few ohm.

## N4  Hysteresis  (`forward_vs_back_hysteresis.png`)

The large outer cycle and the double scroll coexist: turning down, the double scroll survives to
328 ohm; turning up, the large cycle survives to 675 ohm. The bench model reproduces both
(344 and 683 ohm); the constant-component models keep the large cycle to 770-805 ohm.

## N5  Period doubling from the bench  (`cascade_periods.py`)

| | forward | back |
|---|---|---|
| R1 (1 -> 2) | 774.6 +- 1.5 ohm | 768.1 +- 0.2 ohm |
| R2 (2 -> 4) | 755.6 +- 0.8 ohm | not resolved (one period-4 record) |
| R3 (4 -> 8) | 750.5 +- 0.4 ohm | - |
| delta_1 = (R1-R2)/(R2-R3) | 3.7 +- 0.7 | - |

Universal value 4.669. The models' values are in the Feigenbaum section below.

## M5  Lyapunov exponents  (`lyapunov.py --rosenstein --each`)

Forward sweep: a return-map exponent for 286 of the 405 records (the periodic and large-cycle
records have no curve to fit) and a direct (Rosenstein) exponent for all 405. In the double
scroll, 340 to 700 ohm, the map gives a median of 1996 /s and the direct method 1513 /s, a
median ratio of 1.32; the direct exponent exceeds 500 /s from 328 up to 726 ohm and is near zero
or negative on the cycles and the large cycle. Both of the plan's traps show in the data: in the
double scroll the mean return time exceeds the median by 8 % (median over records) and by more
than 10 % in 40 % of the records, whereas in the single-scroll band (700 to 760 ohm) the two agree
to 1 %; and the map exponent reads above the direct one throughout, as the benchmark below
predicts. Per-record diagnostics: `chua/forward_lyapunov_each/`, summary `chua/forward_lyapunov.csv`.

Three long, clean double-scroll records for M4/M5 (two-branch maps with within-branch R2 above
0.9, chosen nearest 450, 550 and 650 ohm):

| record | Rpot (ohm) | maxima | mean / median T (us) | map lambda (/s) | direct lambda (/s) |
|---|---|---|---|---|---|
| forward/trace342 | 450.2 | 1170 | 427 / 420 (+1.8 %) | 2362 +- 429 | 1977 +- 27 |
| forward/trace270 | 549.3 | 1300 | 385 / 330 (+14.2 %) | 2365 +- 175 | 1495 +- 43 |
| forward/trace229 | 650.3 | 1292 | 387 / 350 (+9.5 %) | 3078 +- 102 | 1793 +- 43 |

Back sweep: a map exponent for the 17 chaotic records between 685 and 737 ohm, 60 to 1627 /s,
and a direct exponent for all 74 records, -239 to 1055 /s; the direct value exceeds 500 /s only
between 685 and 713 ohm, and the mean and median return times differ by under 1 % there because
these records rarely switch lobes.

Estimator calibration (`benchmark_lyapunov.txt`, simulated double scrolls quantised like the
bench): with this bench's alpha and beta the map estimate reads 1.05 and the direct estimate
0.92 times the true exponent; on Matsumoto's double scroll 1.42 and 1.09. Quote the map value as
an estimate from the return map and the direct value as the measurement, as the plan says.

## The circuit as measured  (`identify.py` -> `identified.json`, `identified_*.png`)

479 records; 109 are periodic (58 period-1 cycles, 3 orbits around the origin, 48 large cycles)
and were averaged over their cycles, which removes the scope's quantisation. Channel gains:
A + B = 1.003 +- 0.004 over all records.

| element | how | value |
|---|---|---|
| C1 | loop integral of node 1 (exact for any static element) | 11.10 nF at small amplitude; 11.50 nF on the period-1 cycle; 12.56 nF on the large cycle |
| C2 | tank admittance, harmonics of the orbit around the origin (1.4 kHz, 26 harmonics, residual 0.9 %) | 90.5 nF (the period-1 cycles alone give 85.2 nF with L 20.9 mH: they cannot split L*C2) |
| inductor | flux-current loop of every periodic record | L_eff 19.0-19.6 mH and r_eff 3-12 ohm below 2.5 mA; 22.4-23.8 mH and 30-57 ohm at 8-14 mA; lens-shaped Rayleigh loops (`identified_inductor.png`) |
| inductor law | L(i) = L0 + nu abs(i), r(i) = r0 + rho abs(i) | L0 = 19.14 mH (from the onset period; the loops extrapolate to 18.2), nu = 0.64 H/A, r0 = 2.30 ohm, rho = 4244 ohm/A |
| element | integrated KCL over 194 double-scroll and cascade records (residual 1.0 %) and the plateaus of 48 large cycles | breakpoints -6.78, -1.015, +0.875, +6.03 V; slopes +3.89, -0.416, -0.762, -0.415, +4.22 mS |
| element as Kennedy's diode | the eight segment numbers map exactly onto two op-amp stages | RA = 249 ohm, AA = 1.115 (stage A), RB = 22 k, AB = 7.62 (stage B), rails +6.73/-7.56 V (A) and +6.67/-7.73 V (B) |
| element dynamics | i_NR - static law inside |v1| < 3 V of the large cycle | up to 0.71 mA at 0.14 V/us: stage B stuck at its rail while v1 crosses the inner region (op-amp slew rate about 0.5 V/us); 3.6 % residual of a static fit to the large cycle against 1.0 % elsewhere |

The apparent C1 grows with amplitude (11.1 -> 11.5 -> 12.6 nF) because the loop integral
attributes the op-amps' lag to a capacitance; the bench model carries C1 = 10.9 nF plus the
stage-A lag (0.2 nF) and reproduces the small-signal value. Two op-amp constants are not
measurable from two node voltages and are set in `simulate.py`: the stage-B lag (2.4 us, an
op-amp of about 0.5 MHz gain-bandwidth at gain 7.6; 1.2 us moves the low end of the double
scroll from 344 to 390 ohm and nothing else) and the slew rate (0.5 V/us; 0.3 V/us kills the
large cycle 35 ohm early, above 1 V/us the measured deviation loop is not reproduced).

## Models against the bench  (`simulate.py`, transitions from the continuation sweeps)

| | bench | plan's model, nominal C1 = 10 nF | plan's model, C1 = 11.5 nF | identified circuit, constant components (`--model static`) | identified circuit (`--model bench`) |
|---|---|---|---|---|---|
| Hopf point, period | 893-897 ohm, 330 us | 1006 ohm, 338 us | 887 ohm, 330 us | 905 ohm, 331 us | 905 ohm, 331 us |
| first doubling R1 | 774.6 | none: chaos from 900 down | 772 | 845 | 776 |
| chaos onset | 750 | 900 | 720 | 823 | 750 |
| double scroll from (down) | 668 | 866 | 579 | 778 | 697 |
| double scroll ends (down) | 328 | 489 | 519 | 447 | 344 |
| large cycle survives to (up) | 675 | 900 | 769 | 805 | 683 |
| period-1 cycle at 806 ohm: period, V1 range | 346 us, -4.33..-0.67 V | chaotic | 340 us, -4.34..-0.61 | chaotic | 348 us, -4.24..-0.62 |
| large cycle at 320 ohm: period, V2 amplitude, iL amplitude | 356 us, 4.95 V, 11.7 mA | 290 us, 7.6 V | 293 us, 7.2 V | 290 us, 7.1 V | 331 us, 5.08 V, 10.6 mA |
| large cycle at 628 ohm: period | 374 us | - | - (dead) | - | 355 us |

The three constant-component columns show what the plan's model cannot do whatever its numbers:
C1 = 11.5 nF puts the Hopf point in place but leaves the double scroll a 60 ohm band around
550 ohm and the large cycle alive to 770 ohm; the identified constants with a static element
(fourth column) are worse still, because the loss the inductor has at 2 mA and above is what
ends the large cycle and lowers the cascade. Turning on the inductor's amplitude dependence and
the op-amp dynamics (last column) puts every transition within 1-4 % of the bench.

Figures: `chua/simulated_bifurcation_bench.png` (bench model over the measured points),
`chua/simulated_bifurcation.png` and `_nominal.png`, `_static.png` (the others), and the
portrait galleries `chua/simulated_portraits_*.png`, each panel beside the nearest
forward-sweep record.

What the bench model still gets wrong, and why:

- The large cycle's period is 5-7 % short (331 against 356 us at 320 ohm, 355 against 374 at
  628). The whole deficit is in the dwell on the positive plateau (V1 near +6.3 V, stage A at its
  rail), which lasts 25 us longer on the bench; the model's saturated stage is a resistor to a
  fixed rail, the real op-amp's overload recovery is not modelled.
- The small cycle appears within a few ohm at 850 ohm in the model, whereas on the bench it grows
  over 50 ohm below 893. A piecewise-linear element has no curvature between its corners to limit
  a growing oscillation gently.
- The double-scroll onset (697 against 668 ohm) and its lower end (344 against 328) are 2-4 %
  off; the onset is the transition most sensitive to the inner-region slope and breakpoints,
  which are known to about 1 %.

## N5  Period doubling from the models  (`feigenbaum.py`)

| | bench (records) | plan's model, C1 = 11.5 nF | identified circuit, constant components | identified circuit (bench model) |
|---|---|---|---|---|
| R1 (1 -> 2) | 774.6 +- 1.5 | 771.9 | 845.0 | 776.8 |
| R2 (2 -> 4) | 755.6 +- 0.8 | 732.5 | 826.8 | 754.9 |
| g1 = R1 - R2 | 19.0 | 39.4 | 18.3 | 21.9 |
| R3 (4 -> 8) | 750.5 +- 0.4 | not resolved | not resolved | not resolved |
| delta_1 | 3.7 +- 0.7 | - | - | - |

The doubling points are found by bisection on the orbit's period and the sqrt law of the newest
splitting (`feigenbaum_<model>.txt`). In all three models the period-4 orbit gives way to chaos
within the 0.05 ohm resolution of the search (the bench model at 747.4 ohm, with a period-3
window at 746 ohm), so no period-8 window and no delta could be taken from them; the bench
records resolve period 8 over 5 ohm. The bench model reproduces R1 and R2 to 2 ohm; the plan's
model with C1 = 11.5 nF has the first doubling in place but the second 23 ohm too low.

## Data-quality notes

- Records are 200 or 500 ms (570 to 1400 windings), DC coupled, CH1 swings 3.6-5.8 V against
  0.7-1.8 V on CH2 in the oscillating region, and no record has more than 2 % of samples at a rail.
- The vertical range was changed several times within both sweeps, always on all three channels
  together (the divider identity A + B = 1 holds on every range). Between forward trace329 and
  trace330 the divider result steps from 455.6 to 466.9 ohm in a sweep that otherwise falls
  1-2 ohm per record, so forward values below about 467 ohm carry a +2 % offset relative to those
  above; the cascade region is unaffected.
- The 13 forward and 6 back records flagged by the divider are the DC records at the top of the
  dial and the records taken while the knob was still moving. Forward trace14-18 pass the divider
  cut but are near-DC records whose resistance is unreliable (their maxima sit at the equilibrium
  while neighbours at the same nominal resistance oscillate); they are excluded from the onset
  numbers above by the amplitude criterion.
