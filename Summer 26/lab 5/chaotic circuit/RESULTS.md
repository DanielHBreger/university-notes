# Results (rerun of 12 September 2026)

Inputs: `chua/forward/` (418 records, knob turned down), `chua/back/` (80 records, knob turned
back up) and the V-I record `trace1.csv`. Divider fits use R0 = 992 ohm (`chua/scope_data.py`).
Nominal point estimates below come from the pipeline in `README.md`.

**Uncertainty update (15 September):** see [the error analysis](uncertainty/REPORT.md) and
`uncertainty/*_uncertainty.csv`. Fit errors, acquisition brackets and instrument limits are
different quantities. The regenerated presentation figures label the calculated fit
components. The instrument models alone do not recover the missing voltage scales,
calibration history or resistor calibration. Tables below retain nominal values unless noted.

**Onset update (25 September):** the Hopf onset moves from 893 to 907 ohm. rpot.py's divider
fit reads low on small oscillations, which is exactly where the onset is; `chua/onset.py`
places the accepted records near the onset by their mean voltages on rpot.py's scale (N1,
data-quality notes). The records rpot.py flags stay excluded; rpot.py and the sidecars are
unchanged.

## Summary

- The measurements show the regime sequence, including period doubling and hysteresis.
  The divider identity is a consistency check, not a calibration of the resistance axis.
  The identified data give mean A+B = 1.0029 and a maximum deviation of 1.25%; neither
  number establishes relative channel gain accuracy.
- The plan's model (nominal components, the element as Table 1's five segments) runs through the
  bench's whole sequence of regimes in the same order (Hopf point, period 2, 4, 8, chaos, double
  scroll, large cycle, hysteresis), but everything sits 180-260 ohm higher on the dial. Identified
  from the oscillator records themselves (`identify.py`): the inductor is a ferromagnetic-core
  choke whose inductance grows from 19 to 24 mH and whose loss grows from 2 to 57 ohm between the
  smallest cycle and the large cycle (a Rayleigh law, seen directly as lens-shaped flux-current
  loops), C2 is 90 nF rather than 100, and the element is Kennedy's two-op-amp diode with
  G_a = -0.762 mS.
- Node 1's capacitance grows with the swing of V1: the loop integral, taken over every record,
  gives 10.78 +- 0.03 nF on the smallest cycles, 11.4-11.7 in the cascade and the double scroll and
  12.4-12.6 on the large cycle, and a two-parameter Rayleigh law in the swing since V1's last
  turning point (the behaviour of a ferroelectric ceramic capacitor, as the inductor's core
  follows its own Rayleigh law) plus stage B's lag (2.66 +- 0.03 us, measured the same way)
  describes all 479 records to 0.4 %, the chaotic ones out of sample.
- With everything measured from the records (`simulate.py --model bench`, nothing adjusted to the
  transitions) the model puts the onset at 925 +- 3 ohm (bench 907 +- 5) and every later
  transition within 26 ohm of the bench (-20 at the first doubling, -2 at the double scroll's
  start, +25 at its end), except the large cycle's end going up (-59); it also has the bench's
  orbit about the origin between the double scroll and the large cycle. With C1 constant the gap
  grew with the size of the oscillation, to +122 ohm. The report's model (b) leaves stage B's
  measured lag out (datasheet op-amps): the same 25 ohm mean gap, the cascade within 11 ohm, the
  double scroll +25 to +59 ohm late, and no orbit about the origin. See "Node 1's capacitance law" and
  "Simulation: verification of 27 September" (the earlier 0-5 % agreement rested on an error).

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
| records with a divider-fit resistance (residual < 5 %, 0-1000 ohm) | 405 of 418 | 74 of 80 |
| highest accepted record (resistance from `onset.py`; rpot.py reads it low) | trace14, 904.9 ohm (rpot.py 870.5) | trace74, 899.5 ohm (rpot.py 900.6) |
| period-1 | 905 .. 775 | 900 .. 768 |
| period-2 from | 775.7 | 768.1 |
| period-4 from | 755.6 | 744.7 (one record) |
| period-8 from | 750.3 | - |
| single-scroll chaos | 750 .. 722 (periodic window 719-713) | 737 .. 717 |
| single scroll, no lobe change | 721 .. 670 | 713 .. 685 |
| double scroll | 668 .. 328 | (not reached from below: the large cycle holds) |
| small orbit around the origin | 328 .. 315 (period 660-710 us) | - |
| large outer cycle | 314 .. 4 | 3 .. 675 |

The observed double-scroll range sits inside the M2 range. The back sweep's resistance scale
reads about 6 ohm lower than the forward sweep's in the cascade (R1 768 against 775). A relative
gain difference is one possible explanation; it is not independently calibrated. (The earlier
comparison of the onsets, 874 against 880, used the biased small-cycle values; see N1.)

## N1  Hopf point

| | value |
|---|---|
| measured onset (`onset.py`) | 907.4 +- 0.8 ohm: A1^2 of the four smallest accepted limit cycles (897.8-904.9 ohm) extrapolates to zero, chi2 = 0.2, dof = 2, p = 0.91; systematic about +-5 ohm. With the flagged records included: 907.8 +- 0.5 ohm and a direct bracket 907.3-907.8 ohm |
| measured period at onset | 328.40 us (f = 3045.06 +- 0.03 Hz, trace14, the smallest accepted cycle; SE over 10 blocks of the record, timebase accuracy not included); 330.3 us averaged over the 18 smallest cycles (`identify.py`) |
| fit range | the n smallest cycles give R_H = 907.7 +- 1.1 (n = 3), 907.4 +- 0.8 (4, used), 907.1 +- 0.7 (5), 906.5 +- 0.6 ohm (6), all with p > 0.3; at n = 7 the sqrt law is rejected (p = 0.001). k = 3.2 +- 0.4 x 10^3 mV^2/ohm for n = 4 |
| plan formula, nominal C1 = 10 nF, C2 = 100 nF, L = 18 mH, trace slopes, R0 = 992 ohm | 1006 ohm and 2963 Hz (unrounded G_b = -0.41756 mS; the rounded -0.418 gives 1005, 2961) with the left shoulder (the one that applies: the circuit starts oscillating about the negative equilibrium, V1 near -3 V); 965 ohm and 2908 Hz with the right. +-5 ohm from G_b; component tolerances (not measured) dominate: C1 +10 % moves it -79 ohm, L +10 % +44, C2 -10 % -33 |
| plan formula, C1 = 11.5 nF | 886 ohm, 3028 Hz (left shoulder); 856 ohm, 2983 Hz (right) |
| slope at the equilibrium | the circuit starts oscillating about V1 = -3.30 V (mean of V1, forward trace14), inside the left shoulder, where the slope is curved: -0.390 +- 0.010 mS over the middle third, -0.33 to -0.35 mS from local lines over +-0.4 to +-0.8 V (6-12 codes of 129 mV), -0.429 +- 0.004 mS from a quadratic over the shoulder. The prediction moves 25 ohm and 30 Hz per 0.01 mS (nominal parts): -0.39 mS gives 1074 ohm, 3048 Hz; -0.429 gives 977 ohm, 2924 Hz. The C1 that reproduces the onset is 10.9-11.8 nF over these slopes, so the onset does not isolate C1. The frequency is robust (+-2 % over the slopes), the onset is not (+-10 %) |
| identified circuit (C1 10.80 nF, C2 90.5 nF, L0 18.87 mH, r0 2.3 ohm; L0 set by the onset frequency) | 923 +- 1 ohm (C1 and C2 fit errors), 3043 Hz |

`onset.py` (figure `chua/onset.pdf`): on small cycles noise on V1 and V2 attenuates rpot.py's
fitted coefficients (errors in variables), B more than A since V2 swings less, so R0 B/A reads
low: by up to 13 ohm on the accepted records nearest the onset (forward trace14-23, back
trace73-74, V1 half-range under 0.4 V), under 1 ohm above 0.6 V. The divider also holds for the mean voltages,
m3 = a (R0 m1 + R m2)/(R0 + R) + C, with a = A + B from the fit and the offset C calibrated to
reproduce rpot.py on large-cycle records of the same V/div range (forward trace29-36, 0.7 ohm
rms; the back sweep's top range has none and borrows trace69-71, 1.4 ohm rms). The amplitude
A1 of V1's fundamental comes from a least-squares sine fit. The previous value, 893 ohm, was
the same extrapolation on rpot.py's biased values.

The amplitude then grows over about 60 ohm below the onset (peak-to-peak V1 0.8 V at 890,
1.5 V at 871, 3.8 V at 845 ohm); the piecewise-linear models reach full size within a few ohm.

## N2  Phase-portrait gallery  (`gallery.py` -> `chua/gallery.png`)

Four clean forward-sweep records, one per regime the plan names, plus the large outer cycle,
each as V2 against V1 with 4 ms of V1(t) underneath:

| regime | record | Rpot (ohm) | V1 range (V) |
|---|---|---|---|
| limit cycle (period 1) | forward/trace54 | 804.4 | -4.36 .. -0.62 |
| period 2 | forward/trace82 | 760.2 | -3.92 .. -0.20 |
| single scroll | forward/trace164 | 699.7 | -3.38 .. +0.23 |
| double scroll | forward/trace270 | 549.3 | -2.40 .. +1.98 |
| large outer cycle | forward/trace411 | 191.0 | -7.33 .. +6.58 |

The same records appear beside the simulations in `chua/simulated_portraits_*.png`.

## N4  Hysteresis  (`forward_vs_back_hysteresis.png`)

The large outer cycle and the double scroll coexist: turning down, the double scroll survives to
328 ohm; turning up, the large cycle survives to 675 ohm. The loop is wider than the double
scroll's range: turning down, the circuit passes through the small orbit and reaches the large
cycle only below 314.4 +- 0.3 ohm, while turning up it is on the large cycle at 319.5 ohm (back
trace11) and leaves it at 680 +- 5 ohm (half-bracket, trace35/36, 675.0 and 685.4 ohm). The two
sweeps therefore disagree from 314.4 to 680 ohm, a loop 366 +- 5 ohm wide (plus the 0.5-2 %
range-to-range scale; the two ends are on different vertical ranges). The bench model reproduces both
(344 and 683 ohm); the constant-component models keep the large cycle to 770-805 ohm.

## N5  Period doubling from the bench  (`delta_bench.py`; periods from `cascade_periods.py`)

| | forward | back |
|---|---|---|
| R1 (1 -> 2) | 775.7 +- 1.9 ohm | 768.1 +- 0.2 ohm |
| R2 (2 -> 4) | 755.6 +- 0.8 ohm | not resolved (one period-4 record) |
| R3 (4 -> 8) | 750.3 +- 0.8 ohm | - |
| delta_1 = (R1-R2)/(R2-R3) | 3.8 +- 1.0 | - |
| accumulation point, universal delta | 748.9 +- 1.1 ohm | - |
| start of chaos | 748.3 +- 0.7 ohm | - |

Every value is the midpoint of a bracket, with half the bracket as its uncertainty, as for the
other transitions; delta_1 by linear propagation with R2 shared. Brackets between records whose
pattern is steady (the newest splitting keeps its sign through the record): R1 between trace63
(777.62, period 1) and trace62 (773.79, period 2); R2 between trace85 (756.35) and trace87
(754.84); R3 between trace91 (751.11) and trace92 (749.48). trace67 (774.55) and trace94
(750.29) flip the phase of their splitting: too close to a doubling to say on which side, left
out. The whole cascade (trace61-100) is on one oscilloscope range, so the range-to-range scale
cancels in delta_1. The universal 4.669 is 0.9 of the uncertainty away; delta_1 is only the
first ratio.

Earlier values: 3.7 +- 0.8 (`cascade_periods.py`) placed R1 and R3 by the square-root growth of
the splitting from the three nearest records, which included trace67 and trace94; the
presentation's 4 +- 1 used trace67 and trace94 as bracket edges. A square-root-fit analysis
tried on 26 September gave 3.6 +- 0.5, consistent, and was dropped for the simpler brackets.

Near R1 the knob was turned back up for four records (forward trace63-66, 777.6 up to 783.0 ohm,
all on the 26.5 mV-step range): the circuit returned to period 1 and was in a steady period 2
again at trace68 (773.58 ohm). The first doubling shows no hysteresis.

Which records are period 8 (`period8_check.py` -> `period8_check.txt`, `period8_check.png`):
cascade_periods.py labels five records period 8 (trace92-96) from their average lag distances,
but only three hold a steady period 8: trace92, trace93 and trace95 (749.48, 748.96, 749.25 ohm).
In each, maxima four apart differ by 13-75 mV, period 4 is rejected (chi2/dof = 725-1186, 4 dof,
p < 1e-300), nothing repeats only every 16 (p = 0.59-0.99), the differences agree between the two
halves of the record to within 5 mV, and the level deviations are uncorrelated (autocorrelation
-0.20 to +0.05). The period-4 records trace84-91 give p = 0.10-0.95 on the same test. trace94
(750.29 ohm) has the period-8 pattern, but it collapses and re-forms with the opposite phase
during the record. trace96 (748.57 ohm) reverses phase repeatedly and fluctuates by tens of mV:
not a steady orbit. The splittings are 0.5-3 ADC steps (26.5 mV), which is why period 8 was hard
to see on the scope. Chaos therefore starts between the last steady period-8 record (trace93,
748.96 ohm) and the first chaotic one (trace97, 747.66 ohm): 748.3 +- 0.7 ohm, 0.5 +- 1.2 ohm
below the accumulation point.

Below it the maxima still nearly repeat every four cycles (trace97-98, lag-4 distance
35-61 mV against about 3 mV on the periodic records), then every two (trace99 on, lag-2 distance
about 200 mV against 500-570 mV at lag 1): the chaotic bands merge in the reverse order of the
doublings. The model values are below.

## M5  Lyapunov exponents  (`lyapunov.py --rosenstein --each`)

Forward sweep: a return-map exponent for 286 of the 405 records (the periodic and large-cycle
records have no curve to fit) and a direct (Rosenstein) exponent for all 405. In the double
scroll, 340 to 700 ohm, the map gives a median of 1996 /s and the direct method 1513 /s, a
median ratio of 1.32; the direct exponent exceeds 500 /s from 328 up to 726 ohm and is near zero
or negative on the cycles and the large cycle. Both of the plan's traps show in the data: in the
double scroll the mean return time exceeds the median by 10.4 % (median over records) and by more
than 10 % in 51 % of the records (relative to the median; the csv column mean_vs_median_pct is
relative to the mean), whereas in the single-scroll band (700 to 760 ohm) the two agree
to 1 %; and the map exponent reads above the direct one throughout, as the benchmark below
predicts. Per-record diagnostics: `chua/forward_lyapunov_each/`, summary `chua/forward_lyapunov.csv`.

Three long, clean double-scroll records for M4/M5 (two-branch maps with within-branch R2 above
0.9, chosen nearest 450, 550 and 650 ohm):

| record | Rpot (ohm) | maxima | mean / median T (us) | map lambda (/s) | direct lambda (/s) |
|---|---|---|---|---|---|
| forward/trace342 | 450.2 | 1170 | 427.5(32) / 420(4) (+1.8(7) %) | 2362 +- 429 | 1977 +- 27 |
| forward/trace270 | 549.3 | 1300 | 384.7(30) / 330.0(12) (+16.6(7) %) | 2365 +- 175 | 1495 +- 43 |
| forward/trace229 | 650.3 | 1292 | 386.6(27) / 350.0(19) (+10.5(7) %) | 3078 +- 102 | 1793 +- 43 |

How the map changes across the three records (M4): positive-lobe maxima 1.84-2.50 V at 650 ohm,
1.37-1.96 V at 549 ohm, 1.15-1.51 V at 450 ohm; lobe switches 14.1 +- 1.0, 18.0 +- 1.1 and
31.0 +- 1.4 % of turns.

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
and were averaged over their cycles, reducing quantization scatter without removing
calibration error. Divider consistency:
A + B = 1.003 +- 0.004 over all records.

| element | how | value |
|---|---|---|
| C1 | loop integral of node 1 (exact for any static element), every record | 10.78 +- 0.03 nF on the four cycles of the onset fit (A1 < 200 mV; 10.80 +- 0.02 with the earlier, biased estimator), 11.50 nF on the period-1 cycle, 11.47 and 11.70 nF on the chaotic records above and below 669 ohm, 12.56 nF on the large cycle; described by node 1's law (next section) |
| node 1 law | capacitor grows with the swing x of v1 since its last turning point; stage B adds its lag's share while linear | C_s = 10.515 +- 0.033 nF + 0.843 +- 0.009 nF/V min(x, 0.5 V) + 0.1775 +- 0.0012 nF/V max(x - 0.5 V, 0); stage B 0.921 +- 0.010 nF, a lag of 2.66 +- 0.03 us |
| C2 | tank admittance, harmonics of the orbit around the origin (1.4 kHz, 26 harmonics, residual 0.9 %) | 90.5 +- 1.3 nF (fit error) (the period-1 cycles alone give 85.2 nF with L 20.9 mH: they cannot split L*C2) |
| inductor | flux-current loop of every periodic record | L_eff 19.0-19.6 mH and r_eff 3-12 ohm below 2.5 mA; 22.4-23.8 mH and 30-57 ohm at 8-14 mA; lens-shaped Rayleigh loops (`identified_inductor.png`) |
| inductor law | L = L0 + nu abs(d), loss voltage r0 i + rho abs(d) d, with d the excursion of iL from its mean | L0 = 18.87 mH (from the onset frequency, 3045 Hz at 904.9 ohm; the loops extrapolate to 18.2), nu = 0.69 H/A, r0 = 2.30 ohm, rho = 4244 ohm/A |
| element | integrated KCL over 194 double-scroll and cascade records (residual 1.0 %) and the plateaus of 48 large cycles | breakpoints -6.78, -1.015, +0.875, +6.03 V; slopes +3.89, -0.416, -0.762, -0.415, +4.22 mS |
| element as Kennedy's diode | the eight segment numbers map exactly onto two op-amp stages | RA = 249 ohm, AA = 1.115 (stage A), RB = 22 k, AB = 7.62 (stage B), rails +6.73/-7.56 V (A) and +6.67/-7.73 V (B) |
| element dynamics | i_NR - static law inside |v1| < 3 V of the large cycle, node 1 at the record's own capacitance | up to 0.61 mA at 0.14 V/us (0.71 with the small-signal C1, which booked the large cycle's extra capacitance as element current): stage B still at its old rail while v1 crosses the inner region (a spike a few us wide; not reproduced by a lag and slew model, see below); 3.6 % residual of a static fit to the large cycle against 1.0 % elsewhere |

Against the schematic of the element (`Fig1._NR_scheme.pdf`: stage A with R1 = R2 = 220 ohm,
R3 = 2.2 k; stage B with R4 = R5 = 22 k, R6 = 3.3 k; TL082 on +-9 V):

| | schematic | identified from the oscillator | V-I trace fit (M1) |
|---|---|---|---|
| stage A gain AA = 1 + R2/R3, resistor RA = R1 | 1.100, 220 ohm | 1.115, 249 ohm | - |
| stage B gain AB = 1 + R5/R6, resistor RB = R4 | 7.667, 22 k | 7.62, 22 k (fixed) | - |
| Ga = -(1/R3 + 1/R6) | -0.758 mS | -0.762 mS | -0.717 mS |
| Gb = -1/R3 + 1/R4 | -0.409 mS | -0.416 / -0.414 mS | -0.418 / -0.434 mS |
| Gc = 1/R1 + 1/R4 | 4.59 mS | 3.89 / 4.22 mS | 3.62 / 3.79 mS |
| inner breakpoints Vsat/AB | 0.98 V at +-7.5 V rails | -1.015 / +0.875 V | -1.09 / +0.58 V |
| outer breakpoints Vsat/AA | 6.8 V at +-7.5 V rails | -6.78 / +6.03 V | -6.75 / +5.98 V |

The identified element is the schematic to within the resistor tolerances (Ga to 0.5 %, Gb to
1.5 %, AB to 0.6 %), which settles the inner slope the V-I trace resolved poorly (-0.717 against
-0.758 expected). The saturation slope is lower than 1/R1 + 1/R4 on both sides (3.9-4.2 against
4.6 mS, i.e. RA 249 against 220 ohm): a saturated TL082 output is not a stiff rail but has some
tens of ohms of output resistance, which adds to R1. The rails, 6.7 to 7.7 V on a 9 V supply, are
the TL082's typical swing.

The op-amps' speed shows in two ways. The large cycle shows a dynamic deviation of i_NR from the
static law at each crossing of the inner region: a spike of 0.6-1.1 mA, a few us wide, of
opposite sign on the two crossings (stage B still at its old rail). A first-order lag with a slew
limit does not produce that shape for any constants: it gives a broad bump of 0.3-0.4 mA at best
(0.5 V/us). And a lag acts on node 1 as a capacitance, A tau / R, while its stage is linear: node
1's law (next section) measures stage B's that way, 0.921 +- 0.010 nF inside the inner
breakpoints, a lag of 2.66 +- 0.03 us (6.6 times the TL082 datasheet's 0.40 us at gain 7.6). Stage
A is linear almost everywhere, so its share cannot be told from the capacitor; the model keeps the
datasheet value there (0.059 us, 0.26 nF of node 1's small-signal capacitance) and the 13 V/us
slew rate.

## Node 1's capacitance law  (`identify.py` node1_law, 27 September)

Node 1's loop integral, int (v2 - v1)/Rt v1' dt / int v1'^2 dt, removes any static element. If
node 1's capacitance varies, it is the v1'^2-weighted mean of the incremental capacitance along
the record; over a long record the element's term is a vanishing boundary term even without a
cycle. So it is measured on all 479 records, not only the periodic ones:
- estimators, checked on the bench model run with constant C1 and sampled like the scope (400 ns
  or 1 us, 5 mV noise, the 12-80 mV steps of the ranges used), within 0.4 % from 0.15 V
  peak-to-peak up in every regime (`identify.cycle_node1`, `window_node1`; a test):
  - periodic records: cycles marked on v1 low-passed at 20 kHz, the denominator from two
    independent halves of the cycles. The earlier estimator (markers on the raw v1, one average)
    read the smallest cycles low, 2.8 % at 0.15 V peak-to-peak and 0.6 % at 0.3 V (markers on the
    noisy v1 select its noise and sharpen the average). The small-signal value 10.80 +- 0.02 nF
    was biased by that; it is 10.78 +- 0.03 nF.
  - the other records: 20 us windows, sum y_w dv_w / sum dv_w^2, the record cut where v1 passes
    its starting value in the starting direction (removes the boundary term of the DC current).
- what the records show (mean apparent capacitance per group, nF): the four onset cycles 10.78,
  all cycles below 0.7 V peak-to-peak 10.93, the larger shoulder and period-1 cycles 11.38,
  chaotic records above 669 ohm 11.47, double scroll 11.70 (rising to 12.05 as it shrinks
  towards the origin), orbits around the origin 12.01, large cycle 12.56.
  A change of scope range moves it by about 1 % (records 59-60 against 61-67), the floor of any
  fit here.
- the law, fitted over all 479 records (`identified.json` C1_nF.law; statistical errors):

      C_rec = C_s + a0 <min(x, k)> + a1 <max(x - k, 0)> + s_B <stage B linear> + s_A <stage A linear>

  x = |v1 - v_rev|, the swing since v1's last turning point (a turning point counts once v1 is
  back 20 mV from its extreme; the memory fades over 10 ms, so that the equilibrium's slow drift
  while the knob turns is not counted; `integration.swing`, the model's own tracker). C_s =
  10.515 +- 0.033 nF (the capacitor on the onset cycles), a0 = 0.843 +- 0.009 nF/V up to the knot
  k = 0.5 V (0.4-0.6 V fit equally), a1 = 0.1775 +- 0.0012 nF/V beyond, s_B = 0.921 +- 0.010 nF
  (stage B's lag, 2.66 +- 0.03 us), s_A = 0.264 nF (datasheet). Residual 0.046 nF rms (0.4 %); by
  group -0.11 (onset cycles), +0.01 (other periodic below 5.5 V), 0.00 (large cycle), -0.02 and
  +0.01 nF (chaotic records above and below 669 ohm). Fitted on the 109 periodic records alone,
  a1 = 0.176 nF/V and s_B = 0.875 nF, and it predicts the 370 chaotic records to 0.045 nF rms (a
  constant fitted to the periodic records predicts them to 0.35 nF).
- why this law (rms of the fit / of the out-of-sample prediction, nF):
  - the capacitor grows with the swing, not with the voltage: constant 0.35 / 0.35; C(v) even in
    v (a ceramic's DC-bias effect) 0.33 / 0.75; the swing from the 3 ms running mean (the form of
    the inductor's law) with stage B 0.086 / 0.14; the swing since the last turning point with
    stage B 0.046 / 0.045. A capacitance growing with the swing since the last reversal is the
    Rayleigh law of a ferroelectric dielectric: in a class-2 ceramic (BaTiO3) the domain-wall
    part of the permittivity grows linearly with the AC field amplitude, with hysteresis, the
    dielectric twin of the ferrite inductor's law. The steeper first 0.5 V is empirical.
  - the extra 0.92 nF sits exactly between the element's inner breakpoints: scanning the width of
    that step, the fit is best at the breakpoints (scale 1.0-1.1; 0.5, 1.5 or 2 times them predict
    the chaotic records 2.3-3.5 times worse), where a lag of stage B acts as A_B tau_B / R_B and
    nowhere else; a smooth bump of C(v) fits only with its width as a free parameter.
  - the capacitor's type is not recorded (the plan says only "C1 = 10 nF"); 10.5-10.8 nF against
    10 nF, and C2 90.5 against 100 nF, are within class-2 ceramic tolerance. Looking at the part
    would settle it.
- the element does not depend on it: the core KCL fit repeated with the law's charge subtracted
  moves the slopes by at most 0.3 % (Gb +0.15/+0.09 %, Ga -0.28 %) and the inner breakpoints by
  5 mV, the spline's resolution (residual 0.0098 -> 0.0094).
- in the model (`simulate.py --model bench`): C1(x) in the kernel with the same tracker (a test
  checks the kernel's x equals `swing` on the same signal), stage B's lag 2.66 us. The model's own
  loop integrals, sampled like the scope, reproduce the records' at the same resistance
  (verify_models.txt 4.): 12.500 / 12.499, 12.552 / 12.552, 12.575 / 12.605 nF on the large cycle,
  11.809 / 11.815 and 11.621 / 11.581 in the double scroll, 11.50 / 11.44 in the chaotic single
  scroll and cascade, 11.37-11.45 / 11.47-11.54 on period-1 cycles (model / measured).

## Simulation: verification of 27 September

Everything the simulation section depends on was rechecked before writing 2.5.

What holds:
- Circuit equations and topology are those of the plan's Figure 1 (R0 on the C2 side, R between
  node 3 and node 1, L and C2 at node 2, C1 and N_R at node 1; CH1 = V1, CH2 = V2).
- Both compiled RK4 kernels solve their equations (`verify_models.py`): against scipy's DOP853 at
  rtol 1e-10 the difference in V1 over 1 ms is at most 7e-4 V (plan's model, large cycle, a small
  phase drift at the kinks, the same at a quarter of the step) and 9e-5 V (bench model with node
  1's law, the turning-point tracker updated after every DOP853 step as the kernel does); halving
  the step moves no transition by more than one 1 ohm step.
- `pwl_knots` builds the continuous five-segment law exactly (1e-18 A against an independent
  integral of the slopes); `pwl_eval` equals np.interp to 3e-17 A and is 7 times faster.
- The linear-stability code reproduces the plan's Hopf formula to 0.01 ohm (1005.79 ohm,
  2962.8 Hz with the unrounded G_b = -0.41756 mS and the nominal parts). (The report's 2.3 quotes
  1005 +- 5 ohm and 2961 +- 7 Hz, computed with the rounded -0.418 mS: they should read
  1006 +- 5 ohm and 2963 +- 7 Hz; 98 +- 7 ohm and 2.8 +- 0.2 % are unchanged.)
- The averaging time of the inductor's mean current (3 ms) moves no transition by more than
  5 ohm between 1 and 10 ms; 500 ms of maxima per step instead of 100 ms moves neither the start of
  the double scroll nor the end of the large cycle in either model; a 0.2 mV period threshold
  instead of 2 mV moves only the first doubling, by one step.

What was wrong, and is corrected (each checked against the records):
1. Bench model, inductor law. The simulation applied the Rayleigh law to |iL|, while
   `identify.py` measures it on the excursion from the record's mean current. Near the negative
   equilibrium iL carries -1.7 mA of DC, so the simulated circuit had 16.7 ohm of incremental loss
   and 20.1 mH there, where the smallest cycles show 2.3-4.7 ohm and 16-18 mH, and its real Hopf
   point was 846 ohm and 2923 Hz, not the 905 ohm quoted (which came from a linearisation with L0
   and r0). The law now acts on the excursion d = iL - im, with im the mean current over 3 ms
   (`integration.py`); a DC current sees only L0 and the winding resistance. The full Jacobian of
   the corrected model and the Hopf point of its linearisation (`simulate.hopf_point`) agree to
   0.01 ohm.
2. C1 at small amplitude. `identify.py` averaged the loop integral over every cycle on the
   shoulder (0.18-1.89 V peak-to-peak), over which it rises from 10.77 to 11.33 nF, and called the
   mean (11.08 nF) small-signal. The linearisation needs the zero-amplitude value: the four cycles
   of the onset fit (A1 < 200 mV) sit on a plateau at 10.80 +- 0.02 nF. The near-onset records also
   enter `identify.py` at onset.py's resistances; rpot.py reads them up to 34 ohm low, which had
   raised their loop C1 by up to 1.8 %. (Later the same day: that plateau was partly the
   estimator's noise bias; with the corrected estimator the four cycles give 10.78 +- 0.03 nF, and
   the capacitance keeps rising with the swing from there; see node 1's law.)
3. L0 was matched to the mean period of the 18 smallest cycles (330.3 us, finite amplitude) at
   their mean rpot.py resistance (882 ohm). It is now matched to the onset frequency, 3045 Hz at
   904.9 ohm (forward trace14, `onset.py`): L0 = 18.87 mH.
4. The two op-amp constants of the bench model (stage B 2.4 us, 0.5 V/us) were set against the
   ends of the double scroll and of the large cycle, i.e. against transitions the model is
   compared with, and do not reproduce the measured deviation's shape (above). Replaced by the
   TL082 datasheet values: lags gain/(2 pi 3 MHz), 0.059 us (stage A) and 0.40 us (stage B), and
   13 V/us. (Later the same day stage B's lag was measured from node 1's loop integrals, 2.66 us,
   independently of any transition; see node 1's law.)
5. Plan's model, element. The segments were joined where the fitted lines intersect, which puts
   the inner breakpoints at -1.26 and +1.20 V instead of Table 1's -1.09 and +0.58 V (the M1 fit
   is discontinuous there: at +0.58 V the fitted lines differ by 0.17 mA). A continuous law
   through the origin can carry all five slopes and all four breakpoints of Table 1 exactly; the
   plan's model now uses that, so that it is the element the report tabulates.
6. Protocol. Sweeps covered 300-900 ohm with 30 ms of maxima per step; the dial is 0-1000 ohm (the
   plan's model's cascade is at 954-933 ohm) and the records are 200-500 ms. Sweeps now cover
   0-1000 ohm in 1 ohm steps, 20 ms settle and 100 ms of maxima per step, with a 200 ms settle at
   the first step. Regimes are labelled from the maxima: period p if they repeat at lag p <= 64 to
   2 mV, chaos otherwise; double scroll when maxima lie beyond half of each outer equilibrium;
   large cycle above 5 V; a regime must hold for two steps.
7. The model table compared with superseded bench values (775, 750, 668, 675); Table 2's are
   776(2), 748.3(7), 669.3(9), 680(5).

Results (ohm; `simulate_nominal.txt`, `simulate_static.txt`, `simulate_bench.txt`, transitions
in `simulated_transitions_*.json`; each model transition lies between two 1 ohm steps and is
given as the step's midpoint, +-0.5 ohm, +-1 ohm for the first doubling; the model onsets are the
linear-stability values):

| | bench (Table 2) | plan's model: nominal parts, Table 1 element | identified circuit, constant components | identified circuit, constant C1, datasheet op-amps | node 1's capacitor law, datasheet op-amps (the report's model (b)) | node 1's law with stage B's measured lag (the bench model) |
|---|---|---|---|---|---|---|
| onset | 907(5), 3045 Hz | 1006(5) (from G_b), 2963 Hz | 925(3), 3043 Hz | 925(3), 3043 Hz | 925(3), 3043 Hz | 925(3), 3043 Hz |
| period 1 -> 2 | 776(2) | 954(1) | 870(1) | 832(1) | 766(1) | 756(1) |
| period 2 -> 4 | 755.6(8) | 938.5 | 853.5 | 817.5 | 744.5 | 729.5 |
| period 4 -> 8 | 750.3(8) | 934.5 | 849.5 | 813.5 | 739.5 | 723.5 |
| chaos begins | 748.3(7) | 933.5 | 849.5 | 812.5 | 738.5 | 722.5 |
| period-3 window | 718.9(3) - 713.1(5) | none at 1 ohm steps | 838.5 - 837.5 | 802.5 - 801.5 | 725.5 - 723.5 | 705.5 - 703.5 |
| double scroll | 669.3(9) -> 327.9(2) | 886.5 -> 583.5 | 808.5 -> 473.5 | 776.5 -> 449.5 | 694.5 -> 386.5 | 667.5 -> 352.5 |
| orbit about the origin | 327.9(2) - 314.4(3) | none | none | none | none | 352.5 - 298.5 (smaller orbits about the origin, a short rest at 317-313) |
| large cycle, up | to 680(5) | to 927.5 | to 823.5 | to 724.5 | to 632.5 | to 621.5 |

(`simulate_bench_datasheet.txt`, `simulated_*_bench_datasheet.*`: `simulate.py --model bench
--datasheet-lag`. The report leaves stage B's measured lag out and shows this column as model (b);
`report_figures.py` draws it.)

(The constant-C1 column is `--constant-c1`, with the corrected small-signal C1 of 10.78 nF; before
the correction it read 832(1), 815.5, 811.5, 810.5, 800.5 - 798.5, 773.5 -> 446.5, 723.5. Model
onsets are the linearisation's, 924.6 ohm, +-3 from C1's and C2's errors with L0 re-matched to the
onset frequency.)

Reading across: the plan's model has the bench's whole sequence, shifted up by 180-260 ohm. The
measured small-amplitude components (C1 10.78 nF, C2 90.5 nF, L0 18.86 mH, r0 2.3 ohm, the element
as the oscillator shows it) bring the onset to 925 ohm and everything else 80-110 ohm down; the
inductor's loss growing with amplitude (Rayleigh law) moves the cascade and the double scroll
another 24-38 ohm down and ends the large cycle 99 ohm earlier. With C1 constant, what is left
grows with the size of the oscillation: +57 ohm at the first doubling, +64 at chaos, +107 and +122
at the double scroll's ends. Node 1's law removes that pattern. In 2 ohm sweeps, the capacitor
law alone moves the cascade 64-72 ohm and the double scroll's ends 78 and 60 ohm down; adding
stage B's measured lag moves the cascade a further 12-16 ohm and the double scroll's ends 28 and
36 ohm (the effects do not add: stage B's lag alone, with C1 constant, lowers the double scroll's
end by 75 ohm). With both, model minus bench: -20(2) at the first doubling, -26(1) at 2 -> 4, 4 -> 8 and
chaos, -13 and -10 at the period-3 window's edges, -2(1) at the double scroll's start, +25 at its
end, -16 where the large cycle takes over going down, -59(5) at the large cycle's end going up;
the onset stays +18(6). The mean size of the gap over the seven transitions after the onset
falls from 74 to 26 ohm, the gaps now have both signs and no longer grow with the amplitude, and
the model now has an orbit about the origin between the double scroll and the large cycle, as
the bench does (wider: 54 against 13.5 ohm). What is left is not identified; the model's node 1
reproduces the records' loop integrals (node 1's law), so it is not C1.

Without stage B's measured lag (datasheet op-amps, the report's model (b)) the mean gap is the
same, 25 ohm, with another pattern: the cascade within 11 ohm (-10(2) at the first doubling,
-11.1(9), -10.8(9), -9.8(9)), the period-3 window's edges +6.6(6) and +10.4(7), the double scroll
from +25(1) to +58.6(5), then straight to the large cycle (+72.1(6), no orbit about the origin),
the large cycle's end going up -48(5); onset +17(6). This model's node 1 falls short of the
records wherever the orbit spends time inside the inner breakpoints (the lag's share is missing
there): by 0.04-0.08 nF in the chaotic single scroll, 0.12 nF on the large cycle and 0.24-0.41 nF
(2-3.5 %) in the double scroll (549 and 450 ohm), while matching on the period-1 cycles.

The plan's model with the line-intersection element (the old default) goes through the cascade
at 976-933 ohm and then jumps to the large cycle at 895 ohm before any double scroll, a sequence
the bench never shows.

## N6  Shilnikov ratios  (`shilnikov.py` -> `chua/shilnikov.txt`)

Eigenvalues of the Jacobian at the three equilibria, gamma the real one and sigma +- j omega
the pair, from the measured slopes: Table 1's element with the plan's nominal components
and the identified circuit. The double scroll needs the origin to be a saddle-focus with a
one-dimensional unstable manifold (gamma > 0, sigma < 0) and the outer equilibria the other kind
(gamma < 0, sigma > 0); Shilnikov's condition is |sigma| / |gamma| < 1.

| Rpot (ohm) | origin: gamma, sigma (1/s), ratio | outer: gamma, sigma (1/s), ratio | verdict |
|---|---|---|---|
| 806 | +24 200, -5 970, 0.25 | -19 860, +430, 0.02 | geometry right, condition holds |
| 700 | +21 720, -6 490, 0.30 | -24 190, +830, 0.03 | holds |
| 600 | +19 000, -7 010, 0.37 | -28 560, +1 140, 0.04 | holds |
| 500 | +15 700, -7 500, 0.48 | -33 320, +1 390, 0.04 | holds |
| 400 | +11 150, -7 660, 0.69 | -38 600, +1 580, 0.04 | holds |
| 328 | +4 350, -6 240, 1.43 | -42 800, +1 700, 0.04 | fails at the origin |

(identified circuit, 27 September; the plan's model with nominal parts gives the same picture, ratios 0.26 to 0.59 at the
origin between 806 and 500 ohm, and has a single equilibrium below 403 ohm.) The condition
holds throughout the observed double-scroll range 668 .. 328 ohm at the outer equilibria, and at
the origin down to about 350 ohm, where the ratio passes 1 as gamma of the origin shrinks; the
measured double scroll ends at 328 ohm. The frequency of the pair, 3.0 to 3.2 kHz at the outer
equilibria, is the winding frequency seen in every record (return times 330 to 430 us).

Report values (26 September): the plan's recipe as written, the Table 1 slopes with the nominal
components (C1 = 10 nF, C2 = 100 nF, L = 18 mH, R0 = 992 ohm). Same geometry: the origin a
saddle-focus with a one-dimensional unstable manifold, the outer equilibria saddle-foci with
two-dimensional ones. At the outer equilibria |sigma|/|gamma| = 0.045-0.054. At the origin it is
0.34 +- 0.04 at 669 ohm and reaches 1 at 420 +- 40 ohm (G_a = -0.72 +- 0.02 mS; 465 and 378 ohm at
the ends of that range). Three equilibria exist above 397 ohm, so both limits lie above the
observed end of the double scroll, 327.9 ohm: the M2 tension again (with the oscillator's
G_a = -0.762 mS both move down).

## Data-quality notes

- The oscillator drawing supplied with the plan (`Fig2_Chua_Oscillator.png`) labels node 2 as
  CH1 and node 1 as CH2; the plan's own Figure 1, its oscilloscope-settings list and the
  records themselves (CH1 swings several volts, CH2 under 2 V; the divider gives 0-1000 ohm only
  with CH1 = V1) have CH1 = V1, CH2 = V2, which is what every script assumes. The drawing's
  labels are the ones that are wrong.
- The plan quotes R0 = 990 ohm measured; the sidecars use 992 ohm (`scope_data.R0`), the value
  measured on this bench. The difference is a 0.2 % scale on every Rpot.

- Records are 200 or 500 ms (570 to 1400 windings), DC coupled, CH1 swings 3.6-5.8 V against
  0.7-1.8 V on CH2 in the oscillating region, and no record has more than 2 % of samples at a rail.
- The vertical range was changed several times within both sweeps, always on all three channels
  together (the divider identity A + B = 1 holds on every range). Between forward trace329 and
  trace330 the divider result steps from 455.6 to 466.9 ohm in a sweep that otherwise falls
  1-2 ohm per record, so forward values below about 467 ohm carry a +2 % offset relative to those
  above; the cascade region is unaffected.
- The 13 forward and 6 back records flagged by the divider are all at the top of the dial: the
  records where the circuit rests at its equilibrium (forward trace1-11, back trace76-80; the
  divider fit has nothing to use) and the three smallest limit cycles (forward trace12-13, back
  trace75; fitted resistance 70-140 ohm too low). They are excluded from every result. No record
  was taken with the knob moving. The same bias leaves up to 13 ohm on the accepted records up to
  a V1 half-range of 0.4 V (forward trace14-23, back trace73-74); `onset.py` places those from
  their mean voltages. No other result uses records in that range.
- The resistance scale changes with the vertical range. In acquisition order the V1 maxima rise
  monotonically (-2.36 V at trace31 to -0.94 V at trace44), but Rpot jumps at the range changes,
  in either direction (+3.8 ohm at trace36/37, +9.3 at trace42/43, -17.6 at trace43/44 in a sweep
  that otherwise falls 1-6 ohm per record): trace37-39 (the first 14.1 mV-step range) read at
  least 4 ohm high, and trace43, alone on its range, reads 862.1 ohm where its maxima place it
  between 852.9 and 844.5 ohm. The step between trace329 and trace330 is
  the same effect. Estimate: 0.5-2 % between ranges; not correctable without a per-range calibration.
