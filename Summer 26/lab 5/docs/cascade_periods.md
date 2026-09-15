# `chua/cascade_periods.py`

The period of each swept record, read from the maxima that `bifurcation.py` already extracted, and the doubling points R1, R2, R3 that follow, including the 4 → 8 doubling. The kept outputs are `cascade_periods_forward.txt` and `cascade_periods_back.txt`.

## Why not count levels

Grouping a record's maxima into levels with a tolerance cannot see the deep cascade. By period 8 the newest splittings are 10 to 45 mV while a tolerance loose enough to survive the scatter is about 30 mV, so levels merge and the count comes out as 3, or 4, or anything. The splittings are in the data; the counting method cannot resolve them.

## The lag-distance test

A settled period-P orbit has a maxima sequence that repeats exactly every P. So

```
D(p) = < |m_i - m_(i+p)| >
```

collapses to the record's noise floor at p = P and at every multiple, and is large at every other lag. The period is the smallest lag at which D reaches the floor, with no tolerance anywhere. Two by-products:

- D(P/2) is the splitting created by the most recent doubling. It obeys the square-root law of a period-doubling normal form, so D(P/2)² is linear in R and its zero is R_n.
- D(64)/floor is a stationarity check: about 1 on a clean record, growing when the pot was still moving during acquisition.

## `lag_distances(M, drop=0.25)`

D(p) for p in `LAGS = (1, 2, 4, 8, 16, 32, 64)`, over the last three quarters of the sequence (the first quarter may hold a transient).

## `period_of(D, factor=2.0)`

```python
floor = min(D.values())
P = next((p for p in LAGS if p in D and D[p] < factor * floor), 0)
return P, (D[P // 2] if P >= 2 else np.nan), floor, D[max(D)] / floor
```

The floor is the smallest D over all lags. The period is the first lag within a factor 2 of it. Returns `(period, splitting, floor, drift)`.

## `load_sidecar(tag, base='.')`

Reads `<tag>_bifurcation_points.csv` into `{filename: (rpot, maxima in acquisition order)}`. Because `bifurcation.py` wrote the maxima in order, the sequence is intact.

## `bifurcation_points(rows, nfit=3)`

For each n = 1, 2, 3 (P = 2, 4, 8):

- the records labelled period P, sorted by descending R; the highest one is `lo`, the lower edge of the bracket;
- the lowest R among the records labelled P/2 above it is `hi`, the upper edge;
- the `nfit` records nearest the transition are fitted: splitting² against R, straight line, zero crossing = R_n;
- if the root lands inside the bracket it is used ("sqrt fit"); otherwise the bracket midpoint, with half the bracket width reported as a resolution descriptor, not a standard error.

## `main()`

- `--record path`: one raw record through `lorenz_map.maxima`, then its D(p) table, period and, for period ≥ 8, the P levels as mean ± sd over the repeats. This is the diagnostic for one file.
- Otherwise: the sidecar of the named sweep, restricted to `--range` (default 744 to 790 ohm, the cascade), one line per record with period, floor, splitting and drift, records with drift above `--drift-max` flagged as "pot moving" and excluded from the fits. Then R1, R2, R3 with their brackets, the gaps g1 = R1 − R2 and g2 = R2 − R3, and delta_1 = g1/g2 against the universal 4.669. The model's values come from `feigenbaum.py`.

## Uncertainty update

`delta_covariance` propagates R1, R2 and R3 jointly, retaining the shared R2 contribution.
The historical half-brackets are not 1σ fit errors; the log now labels the resulting value
as a resolution sensitivity. A separate uniform-within-bracket calculation is saved in
`uncertainty/cascade_uncertainty.json` and explained in `uncertainty/REPORT.md`.
