# `chua/rosenstein.py`

The largest Lyapunov exponent straight from the time series (Rosenstein, Collins and De Luca 1993). It needs no map and no assumption that the flow has collapsed to a one-dimensional function: it reconstructs the attractor from the two measured channels and watches how fast nearby points separate.

`lyapunov.py --rosenstein` calls `rosenstein()` for every admitted record; the script also runs on its own for one file and writes `<record>_rosenstein.png`.

## `load_channels(path)`

`(t, V1, V2)` through `read_scope`, at least 100 samples.

## `embed(v1, v2, tau, smooth_window=0)`

The state is three-dimensional but only v1 and v2 are measured. The inductor current is stood in for by a delayed copy of v1:

```python
X = np.column_stack([v1[tau:], v2[tau:], v1[:n]])
```

with `tau` a quarter period. Each coordinate is then centred and scaled to unit standard deviation, so v2, which swings a tenth of v1, counts equally in distances. With `smooth_window` both channels are Savitzky–Golay smoothed first; that is what takes the quantisation out of the distances, though it cannot recover what the recording does not contain.

## `nearest_recurrent(X, last, theiler, stride, tree_stride, k_query)`

For every reference point (every `stride`-th sample below `last`) find its nearest neighbour that is at least `theiler` samples away in time. The temporal exclusion is essential: without it the nearest neighbour is the next sample on the same orbit, which never diverges.

Implementation details that keep this fast on a million-sample record:

- The k-d tree holds only every `tree_stride`-th sample; a neighbour need not be sample-exact.
- Samples inside the Theiler window are nearly always the closest, so `k_query` (256) neighbours are asked for at once and the first admissible one taken. References whose whole batch was excluded are retried with the list doubled, in batches sized to bound memory (`batch = max(1, 500000 // k)`), until every reference has a partner or the tree is exhausted.

A test compares the result against an exhaustive search, including the k = 1 case.

## `divergence(X, theiler, k_max, stride, tree_stride, k_query, blocks, k_step)`

Follows every pair forward k samples and averages ln of the distance:

```python
d = np.linalg.norm(X[i + k] - X[j + k], axis=1)
sums[:, c] = np.bincount(block[good], weights=np.log(d[good]), minlength=blocks)
```

The averages are kept per time block (four blocks of the record) as well as overall; the spread of the block slopes is the quoted uncertainty. A zero distance (identical smoothed states, rare) is left out of that step's average rather than contributing −∞.

## `rosenstein(v1, v2, dt, period, tau=None, theiler=None, fit=(0.5, 2.5), follow=3.5, stride=None, smooth=20)`

Defaults derived from the winding period (from `lorenz_map.period_samples`): delay a quarter period, Theiler window one period, smoothing window period/20, pairs followed 3.5 periods, reference points thinned to about 50 000, curve sampled every period/200 samples.

The mean log-distance rises linearly while divergence is exponential and bends over when the pairs are as far apart as the attractor is wide. The slope of the straight part, in 1/s, is lambda:

```python
sel = (ks >= fit[0] * period) & (ks <= fit[1] * period) & np.isfinite(curve)
slope, icpt = np.polyfit(t[sel], curve[sel], 1)
```

The fit window is a fixed range in periods, 0.5 to 2.5 by default. Starting after half a period skips the initial stretch where the nearest pairs are atypically close and relax outward faster than lambda; stopping before the bend avoids saturation. Both the curve and the window are drawn on every diagnostic so the choice can be checked by eye. The returned dict carries the slope, the block-SEM error, the fit R² and everything `draw_curve` needs.

A parametrised test checks that impossible settings (zero dt, zero period, a reversed or oversized fit window, zero stride or delay) raise.

## `draw_curve(ax, res, dt_period, title)`

The curve in units of periods, the fitted window shaded, the fitted line with lambda and R² in the legend.

## `main()`

For each file: load, validate time, find the period, run `rosenstein`, print period, pair count, lambda and the block slopes, and save the figure.
