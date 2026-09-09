# Chua lab verification — 9 September 2026

The analysis contained substantive errors, which have been corrected. All 525 original recordings passed structural checks, all original measurement hashes were preserved, and 34 regression tests passed. All 40 original figures were regenerated and visually reviewed; the final gallery contains 55 figures, including new requirement-specific comparisons and diagnostics.

The available data do **not** establish full agreement with every experimental requirement. In particular, the two double-scroll boundaries are not both verified, and the measured nonlinear-element model does not reproduce all recorded regimes. The results retain those discrepancies rather than treating acquisition filenames or an assumed offset correction as ground truth.

## The resistor split, checked against page 1

Source: [the supplied experiment plan](<C:/Users/danib/Downloads/Lab5_Chua_experiment_plan (1).pdf>), Figure 1 and the divider derivation on page 2. A local image of the schematic is [available here](verification/plan_diagram_0.png).

In CHAOS mode, the series coupling path contains the fixed **990 Ω resistor and the 0–1000 Ω potentiometer**:

`Rtotal = 990 Ω + Rpot`.

The **216 Ω shunt belongs to the separate V–I measurement** and is switched out in CHAOS mode. It must not be added to the oscillator resistance. The nonlinear element is represented by its measured terminal current–voltage characteristic; its internal resistors must not be added separately to the coupling resistance either. Inductor winding resistance is a separate optional model parameter, currently zero because no measurement was supplied.

The node-3 divider is `V3 = A V1 + B V2 + C`, with `Rpot = 990 B/A` for equal channel gains. This follows the schematic's resistor positions. The previously used fixed value of 992 Ω changed the result by only about 0.2%; that was corrected throughout.

| Quantity | Potentiometer | Total coupling resistance |
|---|---:|---:|
| Lower bound from the fitted inner slope | 404.58 Ω | 1394.58 Ω |
| Upper bound from the left shoulder | 1404.87 Ω | 2394.87 Ω |
| Upper bound from the right shoulder | 1315.91 Ω | 2305.91 Ω |
| Record labelled “double-scroll start” | 654.77 Ω | 1644.77 Ω |
| Record labelled “double-scroll end” | 315.63 Ω | 1305.63 Ω |

The physical potentiometer caps the accessible predicted upper value at 1000 Ω. Adding 990 Ω to **both** the prediction and measurements leaves the lower-end shortfall at **88.94 Ω**. The series split therefore does not explain it. Adding the shunt would be inconsistent with the switch position in the schematic.

The “end” record is also a large, nearly periodic orbit: its maxima span only about 0.072 V, and its direct divergence estimate is approximately zero. Its filename alone cannot identify the disappearance boundary. Moreover, the slope-only bound assumes the usual approximately centred characteristic. The measured inner intercept is −0.413 mA; retaining it substantially changes the actual equilibrium count and positions. These are material limitations on the comparison, not arithmetic corrections to the series resistance. See [both resistance scales](verification/M2_M3_resistance_comparison.png) and [equilibrium calculations](verification/resistance_check.json).

## Requirements and evidence

| Plan item | Result of verification |
|---|---|
| M1: V–I curve, separate slopes, fitted breaks, intercept, sweep overlap | Implemented and regenerated. Corrected the horizontal axis to element voltage. Five disjoint segments are optimized over the recorded voltage codes; both shoulders and both saturation branches are reported. Rising/falling comparison is included. |
| M2: predicted three-equilibrium range | Both asymmetric shoulder bounds and conversion to potentiometer values are provided above. The data cannot verify whether the prediction was written down before observing the oscillator. The large intercept limits the slope-only interpretation. |
| M3: two boundaries and three long, clean records | Three 500 ms records with divider residuals below 1% were selected and checked visually. Both boundary locations remain unverified; the labelled lower endpoint is nearly periodic. |
| M4: three Lorenz maps with diagonal | Produced from the same peaks used by M5, with light smoothing around period/20 and a 2% prominence threshold subject to a noise floor. The maps show structure but also multiple branches at the same input; a strictly one-dimensional horseshoe is not established by these plots alone. |
| M5: visited local slopes / mean return time | Implemented with per-record resolution and fit-quality gates, coverage reporting, return-time histograms, window sensitivity, and direct time-series comparisons. Unresolved maps are omitted with reasons. Estimates are not forced to match an expected percentage difference. |
| N1: Hopf comparison | Slope-only predictions are calculated: total 1997.79 Ω / 1958.08 Ω, corresponding to pot 1007.79 Ω / 968.08 Ω and about 2.963 / 2.909 kHz. A matched experimental onset measurement is not established. |
| N2: phase gallery | Added with raw voltage codes and explicitly labelled light smoothing. Acquisition regime names are retained as labels, not automatically certified classifications. |
| N3: bifurcation scatter | Rebuilt for all four sweeps at their individually measured, uneven resistance values. |
| N4: forward/backward comparison | Added overlays for both measurement sets. Differences are visible; sparse sampling and changes in channel gain prevent assigning precise coexistence boundaries from an overlay alone. |
| N5: Feigenbaum ratio | Not claimed. Fine, even steps and adequately resolved successive period-doubling thresholds are not established by these records. |
| N6: eigenvalues and Shilnikov ratios | Calculated for every equilibrium at four representative resistances under both offset assumptions. The magnitude test alone is not a proof of a homoclinic orbit or chaos. Roots inside an artificial joining interval are model-sensitive. |

## Corrected V–I results

Current is `(VV − VI)/216 Ω`, plotted against **VI**, not VV. For the root V–I file, CH1 is the triangular source and CH2 is the element voltage; this differs from the oscillator channel meanings. The drive shape supports this mapping, but the saved CSV cannot independently verify the historical probe wiring.

| Segment | Slope, mS | Fitted element-voltage interval, V |
|---|---:|---:|
| Left saturation | +3.62319 | −8.619 to −6.754 |
| Left shoulder | −0.417559 | −6.754 to −1.093 |
| Inner | −0.717064 | −1.093 to +0.579 |
| Right shoulder | −0.433668 | +0.579 to +5.982 |
| Right saturation | +3.78643 | +5.982 to +7.719 |

The inner intercept is −0.412971 mA. The combined R² is 0.9638, but the inner segment alone has R² ≈ 0.588: coarse voltage/current codes limit the fit. Reported ordinary least-squares standard errors omit gain calibration, correlated errors, drift and quantisation effects. The original outer source trimming range (−9.33 to +8.28 V) was retained to exclude drive-edge regions; the four interior breaks are optimized, not hand-set.

The rising/falling binned current difference is about 0.0265 mA RMS. This comparison does not reveal a large dynamic loop relative to the milliamperes-scale curve. The measured drive frequency is about 9.17 Hz, below the plan's suggested 20–100 Hz; that acquisition difference is recorded rather than hidden. Results: [fit figure](optimized_breakpoints.png), [curve and direction check](current_vs_voltage.png), [full fitted values](diode_fit.json).

## Three records for M3–M5

All three come from the first forward sweep and last 500 ms. Their divider residuals are 0.91%, 0.72%, and 0.60%, and their V1 code spacing is about 43 mV. See [phase portraits and full time records](verification/M3_selected_records.png) and [the three return maps](verification/M4_return_maps.png).

| Record | Rpot, Ω | Maxima | Mean / median return, µs | Map estimate, s⁻¹ | Direct estimate, s⁻¹ |
|---|---:|---:|---:|---:|---:|
| trace105 | 442.04 | 1171 | 427.11 / 415.00 | 3147 ± 205 | 1841 ± 31 |
| trace101 | 546.49 | 1284 | 389.50 / 334.00 | 2342 ± 270 | 1521 ± 49 |
| trace94 | 650.58 | 1320 | 378.70 / 353.00 | 3218 ± 95 | 1917 ± 43 |

The mean/median differences are 2.84%, 14.25%, and 6.79%, respectively. The code uses the mean in every case; the plan's “exceeds 10%” is supported by the middle example but is not universal in these records.

Map error bars combine nominal sampling error and window sensitivity. Direct error bars are block standard errors. Neither includes all systematic uncertainty. Direct fit-window alternatives span approximately 1730–1876, 1482–1547, and 1912–1917 s⁻¹. High straight-line fit R² does not by itself establish an unbiased exponent. The direct estimates are roughly 35–42% lower than the maps here, exceeding the plan's illustrative 10–20%; the values have not been adjusted to meet that expectation.

The neighbouring-trajectory calculation follows the approach of [Rosenstein, Collins and De Luca (1993)](https://physionet.org/files/lyapunov/1.0.0/RosensteinM93.pdf): temporally separated neighbours in reconstructed state space, followed through time, with a slope fitted to mean logarithmic separation. The chosen reconstruction and fit interval remain analysis choices. [Full results and sensitivities](verification/M4_M5_results.json) and the three `trace*_M5.png` figures preserve those checks.

## Code fixes and numerical checks

- Unified strict, named-column, float64 scope reading; invalid times or missing/non-finite samples are rejected instead of silently joined. Timestamp precision is preserved in peak finding and return times.
- Stabilized divider fitting and rejected degenerate/nonphysical fits. Applied the plan's 5% residual threshold and physical potentiometer range consistently. Named records now report fitted resistance alongside the filename value.
- Excluded generated CSVs from raw-record discovery, fixed output paths for trailing separators and other working directories, and distinguished identical sweep folder names from different sets.
- Corrected short/flat-record handling, local regression numerical stability, unresolved/zero derivatives, actual window doubling when a minimum width dominates, and single-neighbour search shape handling. Missing estimates no longer imply a continuous trend through rejected regimes.
- Removed copied nonlinear-fit constants from the simulation. Corrected outer extrapolation, exact-zero equilibria and integration timing. Simulation figures retain both offset assumptions and use comparable voltage axes. Measurement overlays no longer ingest simulation exports or duplicate records.
- Corrected acquisition display timebase scaling and simultaneous-channel decimation; protected mode restoration after failed transfers and closed unused instrument resources. Hardware changes were checked with mocks, not a connected scope.

The 34 tests include independent references: the logistic-map exponent ln(2), exact matrix-exponential evolution of a linear circuit, exhaustive nearest-neighbour comparison, and the eigenvalue growth rate of a linear variational system. These verify behavior rather than merely restating implementation formulas. Ten analysis CLIs also passed help checks from a different working directory.

For the measured circuit model, 0.5 µs RK4 trajectories agreed with an independent DOP853 reference within **0.324 mV** over a 0.5 ms interval at three resistances under both offset assumptions. Halving the step to 0.25 µs changed post-transient V1 standard deviations by less than 0.5% in the six tested cases. Chaotic trajectories need not remain pointwise identical at longer times, and tail quantiles still show finite-record variability. This supports the tested integrator settings; it does not validate the physical calibration. See [simulation validation](verification/simulation_validation.json).

An independent double-scroll benchmark gave a variational exponent of about 1984 s⁻¹, versus 3128 s⁻¹ from the return map and 2339 s⁻¹ from direct divergence. This confirms that finite-record/projection bias can be substantial even when the implementation runs correctly. See [benchmark output](verification/benchmark.log).

Every exported measured bifurcation maximum and Lorenz pair was compared with the audited peaks: 411,846 maxima and 411,413 return pairs across 433 admitted sweep records passed. All 525 original file hashes remained unchanged. See [data audit](verification/data_audit.csv), [output checks](verification/output_checks.json), and [test results](verification/tests.log).

## Figure review and remaining limits

All 55 final PNGs passed image integrity/non-blank checks, and all 40 original figure paths are represented. Visual review covered every original and final grid page, axes, units, diagonal lines, uncertainty displays and simulation comparisons. Improvements include readable labels, resolution and quality annotations, explicit rejected-record status, consistent comparison axes, and new requirement-specific figures. The [gallery](verification/figures.html) links to full-resolution images; the [manifest](verification/figure_manifest.json) records their dimensions and hashes. Original figure copies remain in the ignored `verification/before` folder.

The unresolved experimental points are the historical probe settings and gain calibration, the nonlinear-element current offset, measured inductor losses/component tolerances, and precise boundary/onset locations. The simulation joins independent fitted segments over 0.05 V; those bridges can introduce steep local slopes and extra roots, so equilibria near a bridge are not robust physical predictions. Removing the inner offset improves some qualitative behavior while failing to reproduce others.

The supplied code and available records have been checked and corrected within that evidence. They support a reproducible analysis with explicit limitations; they do not justify a claim of zero possible bugs, verified hardware acquisition, or complete experimental agreement.
