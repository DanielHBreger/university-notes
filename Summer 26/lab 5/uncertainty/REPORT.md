# Measurement uncertainty audit — 15 September 2026

## What changed

The numerical analysis now separates sampling resolution, within-record fit stability, instrument specifications and calibration. Presentation figures display the calculated fit components in simple language. Detailed scope and resistor limitations belong here, rather than in slide footers.

- Calculated paired time-block bootstrap uncertainty for all **405 forward and 74 backward** clean divider fits. Each record uses 400 draws with both 5 ms and 10 ms blocks; the larger standard deviation is reported. The simultaneous channels stay together. These estimates quantify stability within each recording, not absolute resistance accuracy.
- Calculated voltage-code spacing for every channel, per-sample quantization components, and the covariance between element voltage and derived current.
- Recalculated the diode fit's slope/intercept uncertainty with paired time blocks, preserving their covariance and the nominal segment boundaries.
- Corrected the period-doubling uncertainty propagation: R2 occurs in both intervals, so their covariance cannot be discarded. Existing half-brackets were also incorrectly presented as if they were standard errors.
- Added the confirmed 9 V / 0.03 A settings and the IPS 3303S specifications.
- Removed the unsupported claim that the divider identity establishes 0.5% resistance accuracy, and softened quantitative model-agreement claims accordingly.
- Regenerated the presentation figures and the main equal-R bench-model animation. Nominal readings and model parameters have not been shifted: no calibration correction is available to justify doing so.

## Instrument information

Confirmed by the experimenter: Agilent InfiniiVision DSO-X-3014A, two probes at 1:1, ISO-TECH IPS 3303S, both supply channels at 9 V and 0.03 A. The current entry is recorded as a setting, not assumed to be the load current. The third oscilloscope channel's probe setting was not confirmed.

The CSV writer saved only time and voltage values. It queried scales, offsets and waveform preambles but did not save them. No calibration date or uncertainty for the measured 992 Ω and 216 Ω resistors is available. The code therefore stores these quantities as `null`, not zero, in `instrument_inputs.json`. Probe attenuation by itself does not establish probe accuracy or loading.

## Oscilloscope specification and voltage propagation

The applicable manufacturer sheet specifies gain accuracy as **±2% of full scale**, not ±2% of the voltage reading. For the typical single DC-reading envelope, combine this with the offset term and an additional 0.25% of full scale. At 1:1, with S in V/div and O the offset setting:

```
FS = 8 S
B_DC = 0.0225 FS + 0.1 S + 0.002 V + 0.01 |O|
```

The magnified 1/2 mV/div settings use the 4 mV/div input-equivalent scale. Instrument conditions include warm-up and temperature requirements. This is a DC/cursor accuracy envelope, not a distribution of independent noise on every sample. It must not be divided by sqrt(number of samples). Unknown S and O prevent assigning numerical full-accuracy bars retrospectively; ADC code spacing alone does not uniquely recover those settings. [Keysight 5990-6619EN, pp. 18–19, 22](https://www.farnell.com/datasheets/1847045.pdf).

For the observed voltage-code spacing q, the uniform-within-bin model gives:

```
u_quantization(V) = q / sqrt(12)
```

This is a standard uncertainty component, not an estimate of all analogue noise. The audit checks that occupied code gaps lie on the inferred grid. Repeated quantized points are not independent precision measurements, and a smooth interpolated or averaged curve is not evidence of a better absolute calibration.

## Current and the I–V fit

For I = (VV − VI)/Rs, with independent resistor calibration:

```
u(I)^2 = [u(VV)^2 + u(VI)^2 − 2 cov(VV,VI) + I^2 u(Rs)^2] / Rs^2
cov(VI,I) = [cov(VI,VV) − u(VI)^2] / Rs
```

VI appears in both plotted coordinates, giving negative covariance even when the two channels' quantization errors are independent. Omitting that covariance produces an incorrect uncertainty ellipse and an incorrect errors-in-variables model.

For the V–I recording, q = 128.643 mV on each channel. The independent-bin components are **u(VI) = 37.14 mV**, **u(I) = 0.2431 mA**, and cov(VI,I) = −6.385 × 10⁻⁶ V·A. Resistor and channel-calibration terms must be added when known.

Slope results (mS; uncertainty is the block-bootstrap standard deviation conditional on fixed segment boundaries):

| Segment | Slope | Fit component u |
|---|---:|---:|
| Left outer | 3.623 | 0.043 |
| Left shoulder | −0.4176 | 0.0054 |
| Inner | −0.717 | 0.035 |
| Right shoulder | −0.4337 | 0.0036 |
| Right outer | 3.786 | 0.120 |

These replace interpreting the original OLS standard errors as a complete error budget. They do not include breakpoint selection or errors-in-variables bias. The V–I trace has limited independent slow-sweep coverage, so these conditional fit estimates should not be promoted to inter-session repeatability. The script resamples paired time blocks, not independently shuffled voltage and current samples. Figure bands use the full slope/intercept covariance.

## Resistance

The measurement equation is

```
Rpot = R0 (B/A) (g2/g1)
```

For z = (A,B), the fit contribution is J Cov(z) Jᵀ, where

```
J = (−R0 g21 B/A², R0 g21/A)
```

For independently characterized R0 and gain ratio g21:

```
u(Rpot)^2 = u_fit^2 + (Rpot/R0)^2 u(R0)^2 + (Rpot/g21)^2 u(g21)^2
```

Retain covariance if these quantities come from a joint calibration. CH3's multiplicative gain and constant channel offsets cancel in this regression ratio, but CH1/CH2 relative gain does not. A+B≈1 alone cannot determine the three channel gains or the ratio. Same V/div does not prove same gain.

Examples of the calculated within-record fit component:

| Record | Rpot (Ω) | u_fit (Ω) | Voltage quantization u (mV, CH1/CH2) |
|---|---:|---:|---:|
| forward/trace54 | 804.36 | 0.12 | 6.96 / 6.96 |
| forward/trace82 | 760.17 | 0.15 | 7.66 / 7.66 |
| forward/trace164 | 699.74 | 0.24 | 7.66 / 7.66 |
| forward/trace270 | 549.34 | 0.26 | 10.44 / 10.44 |
| back/trace21 | 546.45 | 0.053 | 23.21 / 23.21 |

Forward median u_fit is 0.220 Ω; backward median is 0.109 Ω. These small values do not establish absolute resistance accuracy. Regression predictor quantization, probe loading, calibration, between-record drift and the documented range changes can matter beyond the calculated resampling component. Do not RSS-add quantization again to the resampling scatter without a model separating overlapping contributions.

## Period doubling and Feigenbaum ratio

The nominal estimates remain 774.612, 755.595 and 750.448 Ω. Their acquisition brackets are [774.55,777.62], [754.84,756.35] and [750.29,751.11] Ω. The half-widths 1.535, 0.755 and 0.410 Ω describe resistance sampling, not measured Gaussian standard errors.

With D = R2−R3 and δ = (R1−R2)/D:

```
∂δ/∂(R1,R2,R3) = (1, −(1+δ), δ)/D
u(δ)^2 = gradient · Cov(R1,R2,R3) · gradientᵀ
```

The previous calculation treated the two gaps as independent, dropping cov(R1−R2,R2−R3)=−var(R2) when the transition errors are independent. Substituting the existing half-widths as sensitivity scales gives **δ = 3.695 with a 0.806 propagated resolution sensitivity**, not the old ±0.70. It is deliberately not labeled 1σ.

As a separate explicit interval model, independently uniform transition locations within the three actual acquisition brackets give median δ = **4.18** and a central 95% interval **[3.34,5.38]** (200,000 draws, fixed seed). This does not replace the nominal fitted δ or constitute a full experimental confidence interval. It expresses what the record spacing alone permits. The universal value 4.669 lies inside that conditional interval.

A common multiplicative resistance calibration cancels exactly in δ. Independently assigning the same systematic calibration term to every transition would destroy this cancellation and exaggerate uncertainty. Changes of gain between ranges or sessions need a different covariance model.

## Timing, frequency and Lyapunov estimates

The time-base specification is 25 ppm plus 5 ppm per year of aging. Because calibration age is unknown, 25 ppm alone is not a complete timing limit. For a time-scale fractional standard uncertainty uε, u(T)=|T|uε and u(f)=|f|uε. An inverse-time Lyapunov exponent similarly has u_clock(λ)=|λ|uε. Cursor-specific screen-width terms are not automatically applicable to digital peak fitting.

The 25 ppm component alone, modeled as rectangular, contributes 0.00477 µs at T=330.3 µs and about 0.0216 s⁻¹ at λ=1495 s⁻¹. For a single period determined by two independently rounded sample indices, a sampling component is dt/sqrt(6), but averaging many intervals requires their shared-endpoint covariance. These terms do not describe peak-picking noise, filtering sensitivity or estimator bias.

The displayed Lyapunov error bars retain the existing fit uncertainty, and now have horizontal within-record R errors. At 549.34 Ω the direct result remains **1495 ± 43 s⁻¹ (fit component)**. Rescaling one measured voltage channel by a constant does not by itself change the ideal log-divergence slope, so a simple percentage voltage error should not be added mechanically to λ. Noise, neighborhood selection and model bias require separate robustness analysis.

## Power supply

At the confirmed 9 V setting, the IPS 3303S CH1/CH2 voltage programming **or** readback limit is **±12.7 mV**. A rectangular model gives a 7.33 mV standard component. Do not add programming and readback specifications as independent errors for the same setting. The line/load regulation limits at 9 V are each 3.9 mV; ripple is specified at ≤1 mV RMS over its stated band. These are separate operating characteristics, not independent calibration terms to sum indiscriminately. [ISO-TECH IPS X303 manual, pp. 49–50](https://docs.rs-online.com/a534/0900766b81409dea.pdf).

The 0.03 A current setting has a programming limit of ±10.09 mA. This is relevant to the current-limit setting, not evidence that this current flowed. The voltage specifications require operation in constant-voltage mode and the specified warm-up/temperature conditions. Tracking mode is unknown; older datasheets and the manual also give different tracking figures, so no tracking uncertainty has been selected without knowing the mode/revision.

Supply error changes the circuit operating conditions; it is not an extra independent measurement error to add to every oscilloscope voltage. The nonlinear-element rails in the identified model were inferred from the circuit data and are not equal to the nominal supply voltage. No unit transfer from supply voltage to effective diode rails has been established, so the code does not invent a simulation uncertainty envelope by shifting those rails ±12.7 mV.

## Model comparison and conclusions

Period doubling, differing sweep histories and the broad chaotic regimes remain visible in the measured data. The nominal identified model still improves several transition locations. The existing numerical deviations are useful descriptive comparisons, not proof of agreement within experimental uncertainty. The largest of the tabulated double-scroll boundary discrepancies is about 4.9%, so the earlier blanket “every transition within 1–4%” statement was also too strong.

A defensible full prediction band would require a joint parameter covariance from the identification procedure, retaining the correlations with the same data used to evaluate the model. That covariance is not supplied by the existing point-estimate JSON. The uncertainty code preserves nominal simulations and reports measurement components rather than constructing independent fake error bars on fitted circuit parameters.

## Files and reproduction

- `instrument_inputs.json`: user-confirmed values, settings and explicit missing inputs.
- `forward_uncertainty.csv`, `back_uncertainty.csv`: every clean record's resolution and divider-fit components.
- `diode_uncertainty.json`: paired I–V covariance and fit covariance.
- `cascade_uncertainty.json`: corrected covariance propagation and the separate interval-only assessment.
- `*_lyapunov_uncertainty.csv`: fit results with R components and the identifiable clock contribution.
- `supply_uncertainty.json`: calculations for the confirmed power-supply settings.
- `chua/uncertainty.py`: reusable propagation functions.

Run `python uncertainty/calculate.py`, then `python presentation_figures/make_figures.py --animation` from the repository root. For the original cascade log, run `python cascade_periods.py forward` from `chua/`.

The propagation conventions follow the [JCGM Guide to the Expression of Uncertainty in Measurement](https://www.bipm.org/en/web/guest/publications/guides): retain covariance, distinguish limits from standard uncertainties, and state the assumed distributions. Unknown calibration inputs must never be encoded as zero uncertainty.
