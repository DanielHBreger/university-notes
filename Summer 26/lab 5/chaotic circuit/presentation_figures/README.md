# Figures for the 10-minute presentation

Use PNG for straightforward insertion into PowerPoint or Keynote. SVG copies retain sharp text and lines when resized (dense scatter layers are rasterized). All generated static figures are exactly 1920 × 1080 pixels (16:9), with slide-ready titles. The model-comparison animation is 1280 × 720 pixels (16:9); the measured resistance-sweep animation is 1920 × 1080 for use as an entire slide. All figures use the actual project data, with no smoothing of the displayed measured phase portraits.

## Suggested selection and order

| Figure | Placement | What to say |
|---|---|---|
| 01_title_double_scroll | Title | A measured double-scroll attractor introduces the phenomenon. |
| 02_nonlinear_element | Methods | The nonlinear element has negative differential resistance; the segment fit supplies the initial model. |
| 03_measured_regimes | Main results | The measured orbit changes as Rpot decreases. All four portraits share axes. |
| measurements_resistance_sweep.gif | Full-slide animation | The measured downward sweep, with fixed full and zoom views and a resistance indicator. |
| comparison_equal_r.gif | Full-slide model comparison | Measurements and the ideal simulation at the same resistance, with matched full and zoom views. |
| comparison_equal_r_bench.gif | Full-slide model comparison | Measurements and the model with measured components at the same resistance, with matched full and zoom views. |
| 04_bifurcation | Main results | Each vertical slice contains the maxima from one recording. The right panel resolves period doubling. |
| 05_hysteresis | Main results | Different sweep directions give different attractors at nearby resistances. |
| 06_lyapunov | Main results | Positive direct Lyapunov estimates support local trajectory divergence in the chaotic regimes. |
| 07_identified_inductor | Supporting figure or backup | Inductance and loss vary with oscillation amplitude, motivating the bench model. |
| 08_model_transitions | Discussion | The identified model improves several transition thresholds beyond the ideal model. |
| 12_model_transitions_deviation | Discussion, alternative to 08 | Shorter bars mean the predicted transition is closer to the measurement. |
| 09_equal_r_bench.gif | Model comparison | Play once, then use the corresponding still to discuss the agreement. |
| 09_equal_r_bench_still | Model comparison | An equal-resistance comparison in the double-scroll regime. |
| 10_model_onset_difference | Optional backup | Another equal-resistance comparison near 700 Ω. |
| 11_period_doubling_orbits | Main results | Period-1, 2, 4 and 8 measured orbits, with 16 successive peaks below each portrait. |

For a ten-minute talk, use the phase portraits, bifurcation, Lyapunov figure, and model animation as the central evidence. Show the hysteresis and transition comparison briefly. Keep the inductor detail and second comparison still available for questions.

## Important interpretation notes

- R means Rpot in the plots. Rtotal = 992 Ω + Rpot.
- The hysteresis still compares 549.34 Ω on the downward sweep with 546.45 Ω on the upward sweep. It does not claim exact resistance matching.
- The model animation matches the exact measured Rpot at every frame. It uses 130 selected values, each showing 20 ms; the bench simulation settles for 100 ms and carries all five states to the next resistance. No phase synchronization of trajectories is assumed.
- Both animation panels share the same scales. Below 360 Ω the animation explicitly switches to a full-range view to contain the large orbit. All portraits are projections onto V1–V2, not the full state space.
- The benchmark model uses identified.json and the effective amplifier settings in simulate.py (0.05/2.4 µs lags and 0.5 V/µs slew). These settings are part of the existing model, not independent measurements or new fits.
- Figure 08 reports the rounded continuation-sweep thresholds from RESULTS.md, rather than mixing those values with the separately refined bisection estimates.
- The Lyapunov error bars show the reported fitting uncertainty, not a total error budget. The analysis has additional estimator bias and measurement limitations.
- Near-DC forward traces 14–18 are excluded from the bifurcation and Lyapunov plots because RESULTS.md identifies their resistance estimates as unreliable.
- The raw I–V plot retains the channel offset to match its segment fit. The simulation removes the fitted current offset.
- The identified model uses these characterization data; agreement is not an independent validation experiment.

## Reproduction

From the repository root, with the project requirements and Pillow installed:

```sh
python uncertainty/calculate.py
python presentation_figures/make_figures.py --animation
```

Omit `--animation` to render the static figures 01–08, 11 and 12. `manifest.json` identifies the data sources. The original project figures and earlier animations remain available.

To rebuild the full-slide measured sweep separately, run `python presentation_figures/make_measurements_sweep.py`. It uses the 130 records from `chua/animations/resistance_sweep_metadata.json`, preserving the original 120 ms frame timing (15.6 seconds, looping). It writes a GIF and matching PNG here and updates `chua/animations/measurements_resistance_sweep.gif` and its PNG. The original GIF is kept as `chua/animations/measurements_resistance_sweep_original.gif`. The fixed full view is ±9 V on both axes; the orange box marks the ±4 V / ±1 V zoom. This preserves the original selection, including the near-DC first trace whose resistance estimate is less reliable.

To rebuild both full-slide comparison animations, run `python presentation_figures/make_comparison_sweeps.py` (or choose `--model ideal` / `--model bench`). Both are 1920 × 1080, with 130 frames and the same 15.6-second loop as their originals. The script reads the original sweep metadata, reruns the existing continuation simulation, and preserves the full/zoom views and matched scales. GIFs and matching PNGs are saved here and replace their counterparts in `chua/animations`; the previous GIFs are retained with an `_original.gif` suffix. The earlier two-panel `09_equal_r_bench.gif` remains a separate presentation option.

## Error bars and captions

The figures use short labels: resistance fit errors, segment-fit bands, and Lyapunov fit errors are one-standard-deviation components. Acquisition intervals in the bifurcation figure describe resistance sampling. Detailed error calculations and instrument conditions are in [the uncertainty report](../uncertainty/REPORT.md).

The 1/2/4/8 comparison uses forward traces 58, 82, 90 and 92. Their periods are checked by the existing lag-distance analysis before rendering. The upper panels are raw voltage portraits; the lower panels use peaks from the same smoothed analysis as the bifurcation diagram.

## Slide readability updates

- The I–V title states how the measurement defines the nonlinear model.
- The bifurcation overview outlines the exact zoom region with an orange dashed box; the detail panel has a matching border.
- Increasing-resistance data use larger purple triangles and a solid, enlarged legend marker.
- The Lyapunov plot has no point annotation arrow. Its second panel expands 540–560 Ω to display the original fit uncertainties with capped error bars.
- Figure 12 compares absolute model deviations in ohms using the same rounded thresholds as figure 08, which is retained. These bars are differences, not measurement errors.
