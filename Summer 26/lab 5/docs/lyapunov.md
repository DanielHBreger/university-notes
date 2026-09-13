# `chua/lyapunov.py`

M5: the largest Lyapunov exponent from the first-return map, with the direct time-series estimate (`rosenstein.py`) run alongside when asked. The plan's formula:

```
lambda ≈ < ln |f'(M_n)| > / <T>
```

`f'` is the local slope of the return map at each maximum actually visited, `<T>` the mean return time between successive maxima. Both come from the same peaks `lorenz_map.maxima` found.

Outputs: `<first input>_lyapunov.png`, `.csv` (one row per admitted record, columns listed in `CSV_FIELDS`), and with `--each` a diagnostic PNG per record in `<first input>_lyapunov_each/`.

## The two traps the plan names

- **Mean, not median return time.** A lobe switch takes longer than an ordinary winding, so the return-time distribution has a tail and the median sits below the mean. The formula needs the mean (the total time divided by the number of returns). Both are reported and `regime_report` summarises the gap by regime, so the choice is defended by the data.
- **Average over the points visited, including the turning point.** Near the top of the hump the slope goes through zero and ln|f'| plunges. Leaving those points out would bias lambda upward. They are kept, with the regularisation described under `mean_ln_slope`.

## What it takes to have a slope at all

The formula assumes M(n+1) is a function of M(n). A periodic orbit is a few tight clusters of noise: there is no curve, and a fit through the clusters measures the distance between them, not a derivative. The gates below refuse such records with a reason instead of printing a number. Periodic windows carry no lambda by design.

## `split_branches(m_n, gap)`

Sorts the points by M(n) and cuts the sorted index at any gap larger than `gap` (5 %) of the span. A double-scroll map lives on two separated clusters; a fit crossing the void between them would measure the void. Each branch comes back in ascending M(n) order.

## `local_fit(x, y, width, min_points)`

A straight line fitted at every point over the window |x' − x| ≤ width/2, returning slope, fitted value and slope standard error per point, with NaN where the window holds fewer than `min_points` points or has no spread.

Two details matter:

```python
origin_x, origin_y = x.mean(), y.mean()
x, y = x - origin_x, y - origin_y
```

Centring before forming prefix sums. The sums of x² over a window are differences of large cumulative sums, and for data at an offset the cancellation destroys the precision (a test puts the data at 10⁶ and checks the slope to 1e-7).

```python
lo = np.searchsorted(x, x - h, 'left')
hi = np.searchsorted(x, x + h, 'right')
cx = np.concatenate([z, np.cumsum(x)])   # and cy, cxx, cxy, cyy
```

The window at every point is found by binary search on the sorted x, and every sum over a window is a difference of prefix sums. Each fit is O(1) and the whole pass O(n) whatever the width. The slope, prediction and standard error are the usual closed forms, evaluated with `np.errstate` suppressed and then masked.

Why a window of fixed **width** in M(n) and not a fixed number of neighbours: a fixed count squeezes the x-range of the window to a sliver on a dense cluster while M(n+1) keeps its full noise, and the slope of pure noise then comes out above 1.

## `fit_map(m_n, m_next, width, min_width, gap, min_spread, min_points)`

Runs `local_fit` on each branch. A branch narrower than `min_spread` (0.1 V) is a cluster and is skipped; a branch with fewer than `min_points` is skipped. The window is `max(width · branch span, min_width)`. Returns the per-point slopes, fitted values and standard errors, the within-branch R², and the number of branches and points fitted.

The R² is computed **within branches**: `ss_tot` sums each branch's own variance about its own mean. Against the pooled variance a fit that explains nothing inside either cluster of a period-2 orbit would still score 0.99, because the pooled variance is the distance between the clusters.

## `mean_ln_slope(slopes, se)`

```python
floor = se[ok] / np.e
unresolved = a <= se[ok]
a = np.maximum(a, floor)
```

Points whose |slope| is below its own standard error are "unresolved": the data cannot tell the slope from zero. Their |f'| is floored at se/e, so ln|f'| is at most one unit below ln(se). This is a heuristic motivated by averaging log|s| across a zero crossing, not a correction with guaranteed accuracy; the unresolved fraction is reported alongside so the reader can see how much of the average rests on it.

An exactly zero slope (possible on quantised maxima) makes the mean −∞ and the function returns that, rather than deleting the point, which would bias the invariant average upward. A test pins this behaviour.

## `lyapunov(M, t_peaks, ...)`

One record. Fills a result dict with every field present, then applies the gates in order, returning early with a status string:

1. invalid inputs or fewer than two maxima;
2. return times not strictly increasing;
3. fewer than `2·min_points + 2` maxima;
4. map spread below `min_spread` (period-1: one point, nothing to fit);
5. map spanning fewer than `min_steps` (30) quantisation steps: with 201 mV steps a 3 V map is 15 codes and no local slope is reliable;
6. no branch wide enough to fit;
7. fewer than `min_curve_frac` (50 %) of the points on fitted branches;
8. within-branch R² below `min_r2` (0.8).

Only then is lambda computed. The uncertainty combines three parts in quadrature:

```python
out['lam_stat'] = sem / mean_T                 # standard error of <ln|f'|>
out['lam_sys']  = sys_ln / mean_T              # half the change when the window is doubled
out['lam_err']  = hypot(hypot(lam_stat, lam_sys), |lam| · rel_T)   # plus the SE of <T>
```

The window sensitivity uses a **doubled** window, not a halved one: on quantised maxima a half window holds one or two distinct M(n) values. A test checks that the doubled window changes the estimate even when the floor dominates. These errors describe precision; serial correlation can make them optimistic, and the map projection can bias the value itself (see [benchmark_lyapunov.md](benchmark_lyapunov.md)).

## `diagnostic(...)`

The per-record figure: the return map coloured by ln|f'| with the local fit in black, unresolved points in grey and unfitted points always shown (dropping them makes a cloud look like a curve); the return-time histogram with mean and median marked; and the Rosenstein divergence curve when computed.

## `regime_report(label, rows)`

Groups the admitted records by number of fitted branches (one, or two or more) and prints the average and maximum |mean − median| return-time gap in each group. Descriptive only: the branch count does not establish the number of physical scrolls.

## `rpot_from_name`, `single_records`

Support for records named with their resistance (`chaos - 716.9 ohm.csv`). The value from the name is the nominal one; if the file has a CH3 the divider fit is run and its value replaces it, with the source recorded in the CSV.

## `main()`

- Splits the inputs into single CSVs and folders; folders go through `collect`, files through `single_records`.
- For each record: `lyapunov(...)` with the record's quantum; if `--rosenstein`, the direct estimate is run for every admitted record, whether or not the map was accepted, and its result stored alongside.
- Prints one line per record (status, or the full set of numbers), plots map estimates as filled circles and direct estimates as open squares with error bars against Rpot, and writes the CSV.
- If nothing could be plotted the figure says so and points at the CSV for the reasons.
