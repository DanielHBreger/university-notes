# `chua/lorenz_map.py`

Finds the successive local maxima of one channel of a record, draws the first-return map M(n+1) against M(n), and defines which records of a sweep count. The last part is the important one: `bifurcation.py` and `lyapunov.py` both build on `maxima` and `collect` from this file, so the three figures agree by construction about what a record's maxima are.

Outputs, in `<folder>_lorenz/`: `lorenz_grid_NN.png` (pages of small panels), `lorenz_pairs.csv` (every (M(n), M(n+1)) pair), `lorenz_summary.csv` (one row per record), optionally `each/` and `lorenz_pooled.png`.

## Why the return map

A periodic orbit of period P visits P maxima in turn, so its map is P points. Chaos on a single scroll gives the one-humped curve; the double scroll gives two clusters. It is the simplest picture of the dynamics that a single channel supports.

## `quantum(x)`

The smallest gap between the distinct values present in the channel: the scope's quantisation step, 43 mV or 201 mV depending on the range in use. It sets the floor of everything downstream: how finely a maximum can be known, and what prominence a peak must have to be real.

## `peak_prominence(x_range, prominence, q)`

```python
return max(prominence * x_range, 2.0 * q)
```

The nominal threshold is 2 % of the signal range, but never less than two quantisation steps. On the coarse-range records 2 % of the range is under one step, and then every wiggle of the smoothed staircase would count as a maximum.

## `period_samples(x, limit=200000, max_lag=50000)`

The winding period in samples, used only to size the smoothing window and the minimum separation of peaks. Two stages:

1. **Does it oscillate at all?** The autocorrelation is computed through the FFT (`np.fft.rfft` on the zero-padded signal, multiplied by its conjugate, inverse transformed). A direct `np.correlate` would cost seven seconds per record. If there is no autocorrelation peak above 0.1 within `max_lag`, or the first one is under 20 samples, the record is not oscillating and the function returns that lag (or 0).
2. **How long is one winding?** Not from the autocorrelation. On the low-resistance double-scroll records the lobe switching decorrelates successive windings, the first autocorrelation peak lands at six to eight windings, and a smoothing window sized from it would flatten the very maxima the chain is built on. Instead, the signal is lightly smoothed (window 9), its peaks found, and the gaps between peaks examined:

```python
g = np.diff(pk)
g = g[g > 0.5 * float(np.percentile(g, 75))]
return int(np.median(g))
```

Quantised flat tops produce duplicate peaks a few samples apart. Duplicates can only shorten gaps, never lengthen them, so the upper quartile of the gaps is a safe stand-in for the period; gaps under half of it are dropped and the median of the rest is the period. If fewer than five gaps exist the autocorrelation lag is used.

## `_blank(dt, per, status)`

The `info` dict with every field present, whichever way `maxima` exits, so callers can index it without checking.

## `maxima(path, ch='CH1', prominence=0.02, period=None)`

Validates the arguments, reads only the time column and the one channel (a fifth of the parsing time of reading all four), and calls `maxima_from_arrays`.

## `maxima_from_arrays(t, x, prominence=0.02, period=None)`

Returns `(M, info)`.

- `dt = validate_time(t)` and `per = period_samples(x)` unless a fixed period was given.
- `info['quantum']`, `info['clip_pct']`: the quantisation step and the percentage of samples pinned at the extreme values. A record that runs past the scope's input range sits at the rail for a large part of every cycle and its "maxima" are all the clip level; `collect` drops such records.
- `if per < 20`: a lag that short is the noise floor, not an orbit. Status "no periodicity found", no maxima.
- Smoothing: Savitzky–Golay of order 3 with window `max(5, (per // 20) | 1)`, about a twentieth of a period, forced odd. The window is per record because the sweep folders change time base partway (400 ns then 1 µs per sample), so any fixed window would be wrong for half of them.
- Peak finding:

```python
p, _ = find_peaks(xs, prominence=peak_prominence(rng, prominence, info['quantum']),
                  distance=max(per // 4, 1))
```

`distance` is the other half of the flat-top problem. Two equal plateaus a few samples apart both pass the prominence test (each looks down to the deep valley on its far side); without a minimum separation each such top gives two maxima, a spurious point on the diagonal of the return map and a return time of a few microseconds. A quarter of a winding keeps every real maximum (the shortest real return is about 0.75 of a winding) and merges the duplicates.

- `info['t_peaks'] = t[p]`: when each maximum happened, so return times are measured from these same peaks rather than re-detected elsewhere.

## `collect(folder, ch, prominence, period, max_residual=5.0, max_clip=2.0, max_files=None)`

The single definition of "a record that counts". Returns `(records, dropped)`, where each `Record` is `(name, rpot, residual, M, info)` and `dropped` pairs each rejected name with the reason. The gates, in order:

1. no entry in `<folder>_rpot.csv` (the sidecar is required here);
2. Rpot non-finite or outside 0 to 1000 ohm;
3. divider residual above `max_residual` (5 %, the plan's cut);
4. `maxima` raised;
5. fewer than two maxima (not oscillating);
6. more than `max_clip` percent of samples at the rail.

The bifurcation diagram and the Lyapunov plot both call this, so they are always built from the same set of records.

## Plotting

- `frame(ax, lim)` draws the identity diagonal and forces square, equal axes; a period-1 orbit lands on the diagonal.
- `draw(ax, M, title, lim)` is one panel: a scatter of `M[:-1]` against `M[1:]`.
- `write_csvs` writes the pairs (only for records with status `ok`) and the summary (every record, with its status).
- `plot_grid` lays the panels out in pages of `rows × cols` on common axes, so panels are comparable across the sweep. Excluded records are drawn on a grey background with the reason written in.
- `plot_each` writes one full-size map per record into `each/`.
- `plot_pooled` puts every record's pairs on one map, coloured by Rpot.

## `main()`

Resolves the folder, lists its records, loads the Rpot sidecar if present (only to label panels; this script does not require it), runs `maxima` on every record, and assigns a status: clipping, no valid resistance fit, resistance out of range, or residual too high. Common axis limits come from the records with status `ok`. Then the CSVs and the figures.
