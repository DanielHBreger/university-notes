# `diode_analysis.py` (and the wrappers `find_breakpoints.py`, `plot_current.py`)

M1 and M2 of the plan: fit five straight segments to the measured V–I characteristic of the nonlinear element, report the slopes and breakpoints, and derive from the slopes the range of Rpot for which the circuit has three equilibria.

The two wrappers only import `main` from this file and call it. They exist because the plan names the steps separately; `find_breakpoints.py` is the one to run.

## The V–I record

`trace1.csv` was taken with the element driven through a 216 ohm shunt (`RS`). CH1 is the source side of the shunt, CH2 the element side. So

```
i = (CH1 - CH2) / 216
```

is the current into the element and CH2 is the voltage across it. The fit is made against CH2, the element voltage (`voltage_basis = 'element'` in the JSON). `simulate.py` refuses a JSON fitted against anything else.

Outputs, beside the script: `diode_fit.json`, `diode_segments.csv`, `current_vs_voltage.png`, `optimized_breakpoints.png`.

## `block_average(arr, window)`

Averages consecutive blocks of `window` samples. Kept for the higher-resolution use case where the segment fit would otherwise see more than 2500 distinct voltage codes; not used in the default run.

## `fit_segments(x, y, n_segments=5, min_points=50, min_width=0.3)`

The heart of M1. Finds the partition of the voltage axis into five contiguous segments that minimises the total sum of squared residuals of five independent straight-line fits. It is a global optimum, found by dynamic programming, not a search from a starting guess.

The trick that makes it feasible: the scope quantises voltage, so the record has only a few hundred distinct voltage values ("codes"). Every sample with the same code must fall in the same segment, and a breakpoint can only sit between two adjacent codes. So the problem is: choose four cut positions among m codes.

Step by step:

```python
order = np.argsort(x, kind='stable')
x, y = x[order], y[order]
codes, starts = np.unique(x, return_index=True)
cuts = np.r_[starts, len(x)]
```

Sort by voltage; `codes` are the distinct voltages and `cuts[a]:cuts[b]` is the slice of samples covering codes a to b−1.

```python
prefixes = [np.r_[0., np.cumsum(v)] for v in (x, y, x*x, x*y, y*y)]
```

Prefix sums of x, y, x², xy, y². With these, the sums over any slice are a subtraction, so the SSE of a straight-line fit to any code range is O(1).

```python
for a in range(m):
    b = np.arange(a + 2, m + 1)
    n = cuts[b] - cuts[a]
    sx, sy, sxx, sxy, syy = [v[cuts[b]] - v[cuts[a]] for v in prefixes]
    vx = sxx - sx*sx/n
    valid = (n >= min_points) & (codes[b-1] - codes[a] >= min_width) & (vx > 0)
    sse = syy - sy*sy/n - (sxy - sx*sy/n)**2 / np.maximum(vx, 1e-300)
    costs[a, b[valid]] = np.maximum(sse[valid], 0)
```

`costs[a, b]` is the SSE of the best line through codes a..b−1. The closed form is the usual one: total variance of y minus the part explained by the regression, `Sxy²/Sxx`. A range is only allowed if it holds at least 50 samples, spans at least 0.3 V and has spread in x. Everything else stays at infinity.

```python
dp = np.full((n_segments+1, m+1), np.inf)
dp[0, 0] = 0
for k in range(1, n_segments+1):
    for b in range(1, m+1):
        vals = dp[k-1, :b] + costs[:b, b]
        a = int(np.argmin(vals))
        dp[k, b], prev[k, b] = vals[a], a
```

`dp[k, b]` is the least total SSE for covering codes 0..b−1 with k segments. It is the best k−1 segment cover of some prefix plus the cost of one final segment. `prev` remembers the cut. The walk back from `dp[5, m]` recovers the four cuts.

```python
bounds = [float(x[0])] + [float((codes[e-1]+codes[e])/2) for e in ends[1:-1]] + [float(x[-1])]
```

Each breakpoint is placed halfway between the last code of one segment and the first of the next.

Then each segment is refitted with `np.polyfit(..., cov=True)` to get slope, intercept, their standard errors and an R². The function returns the list of segment dicts and the combined R² of the whole partition. Labels are `Gc left, Gb left, Ga inner, Gb right, Gc right` in Chua's notation: Ga the inner negative slope, Gb the shoulders, Gc the saturation segments.

A test builds a synthetic five-line signal with known slopes and checks the fit recovers them to 1e-8.

## `slope_ranges(fits, r0=R0)`

M2. The equilibria are where the element curve meets the load line i = −v/Rt. There are three of them exactly when the load line is steeper than the shoulders and shallower than the inner segment, so

```
-1/Ga < Rt < -1/Gb
```

with both slopes negative and Ga < Gb. The function reports this for each shoulder separately, as total resistance and as Rpot (minus R0), and sets `valid = False` if the slope order is wrong. A test checks that an invalid order is not hidden.

## `analyze(path, source, element, rs, source_range)`

The full M1 computation.

- Reads the CSV, forms the current, and keeps only samples where the source voltage is within `source_range = (-9.33, 8.28)`. Outside that range the drive amplifier saturates and the samples are not the element's characteristic.
- Calls `fit_segments` on (element voltage, current).
- Sweep-direction check: smooths the source voltage with a Savitzky–Golay filter, differentiates it, and marks samples as rising or falling. The current is then averaged in 80 voltage bins separately for the two directions; the RMS of the difference (`directional_rms_difference_A`) tells whether the element shows hysteresis or a dynamic loop at the drive frequency. `current_vs_voltage.png` plots the two directions and their difference.
- Drive frequency from the peaks and troughs of the smoothed source voltage.
- Assembles the result dict that becomes `diode_fit.json`: schema version, column names, the voltage basis, the shunt, R0, the fit range, the segments, the M2 ranges, and an honest `uncertainty_note`.

## `plot_analysis(result, arrays, out)`

Two figures. `current_vs_voltage.png` is the characteristic as two direction-averaged curves with the five fitted lines on top and a lower panel of rising minus falling. `optimized_breakpoints.png` is the raw scatter with each segment's line in its own colour, dotted breakpoint lines with their voltages, and a residual panel per segment.

## `main()`

Parses the arguments (defaults reproduce the shipped outputs), runs `analyze`, writes the JSON and the CSV of segments, draws the figures and prints the slopes, the M2 ranges and the sweep-direction check. The printed reminder that "M2 is a slope-only necessary estimate" is there because the three-equilibria condition from slopes alone assumes the segments meet at the origin; the actual intersections with the fitted offset are what `simulate.equilibria` computes.
