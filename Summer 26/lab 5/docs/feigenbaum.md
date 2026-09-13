# `chua/feigenbaum.py`

N5 for the models: the period-doubling points R1, R2, ... of the negative single-scroll branch, for the ideal model, the static identified model or the bench model, and the Feigenbaum ratios between them. The kept outputs are `feigenbaum_ideal.txt`, `feigenbaum_static.txt`, `feigenbaum_bench.txt`.

## Method

Doubling points are not found by counting attractor levels; near a bifurcation the new levels are closer than any clustering tolerance. Each R_k comes from the normal form of the doubling. For a settled 2^k orbit the maxima m_i repeat every 2^k, and the splitting created by the most recent doubling,

```
s_k = < |m_i - m_(i + 2^(k-1))| >
```

is zero on the 2^(k−1) orbit above R_k and grows as sqrt(R_k − R) below it, so s_k² is linear in R with its zero at R_k. The steps are:

1. locate the 2^k window by bisection on the period read from lag distances (the same test `cascade_periods.py` applies to the bench records);
2. sample s_k at a dozen resistances inside the window;
3. fit a quadratic to s_k² and take its root nearest the top of the window.

The orbit at every resistance is settled from the continuation state of the previous one, as in the sweeps.

## `class Model`

Wraps one model so the rest of the script does not care which.

- `__init__(kind, dt, tau_b, slew)`: for `bench`, the parameter vector from `simulate.bench_params`, the static five-segment element (used only to find equilibria for starting points), the small-signal inductor loss `rL = P[4]` and dt = 0.1 µs. For `static`, the identified small-signal components are pushed into `simulate`'s globals with `set_components`, and the ideal three-state kernel is used with dt = 0.5 µs. For `ideal`, the M1 element and the nominal-plus-calibrated components.
- `start(rpot)`: places the state on the negative outer equilibrium, nudged by +0.05 V in v1, with the inductor current consistent with it (`iL = v/(R0 + rpot + rL)`), as the forward sweep does.
- `hold(rpot, t_settle, t_collect)`: settle then collect through `kern.hold` or `kern.bench_hold`, carrying the state along. Returns the maxima.

## `period_of(m, factor=2.0, drop=0.25)`

The lag-distance test adapted to a noise-free simulation:

```python
floor = max(min(D.values()), 2e-4)   # a settled orbit repeats to well under a millivolt
if D[1] < 0.02:  return 1            # period-1: no structure at all
if floor > 0.25 * D[1]:  return 0    # no lag brings the sequence back onto itself: aperiodic
```

Both guards were added after failures. A deterministic period-1 orbit has D = 0 at every lag, so the raw floor is 0 and every ratio test is garbage; the 2e-4 V floor fixes that. A chaotic orbit has similar D at all lags, so the smallest of them is within a factor 2 of the "floor" and the orbit was being labelled periodic; requiring the floor to be well below D(1) fixes that. Sequences shorter than 70 maxima return 0.

## `splitting(m, k, ntail=256)`

s_k over the last 256 maxima: the mean |m_i − m_(i+2^(k−1))|.

## `find_window(model, k, r_hi, r_lo, t_settle, t_collect, tol=0.05)`

Two bisections. The first, between `r_lo` and `r_hi`, finds the top of the 2^k window: the resistance where period 2^(k−1) turns into 2^k. The second, between `r_lo` and that top, finds the bottom: where 2^k turns into 2^(k+1) or worse. Each probe restarts on the equilibrium, settles, then settles again and collects (the extra short hold lets the state pass the starting nudge). `r_lo` stays at the range's low end for every k, so a chaotic or higher-period probe below the window always counts as "not this period".

## `doubling_point(model, k, top, bottom, t_settle, t_collect, npts=12)`

Twelve resistances from just below the top to just above the bottom, visited from the top down with the state carried along. s_k at each, then

```python
coef = np.polyfit(rs[ok], sp[ok] ** 2, 2)
roots = ... real roots above the highest sampled R ...
return float(roots.min())
```

A quadratic rather than a line, because far from R_k the square-root law picks up its next term. The root is required to lie above the sampled range, on the side where s_k vanishes.

## `main()`

For k = 1..`--kmax`: find the window (with 30 % of the settle time, since a period label needs less settling than a splitting), stop if the window is narrower than 0.2 ohm (unresolved), else fit R_k and move `r_hi` down to the top of that window. Then the gaps g_n = R_n − R_(n+1), the ratios delta_n = g_n/g_(n+1) against 4.669, and the accumulation point estimate R_inf ≈ R_k − g_(k−1)/(delta − 1). Written to `feigenbaum_<model>.txt`.

In all three models the period-8 window was not resolved at the default settle time, so each model gives R1, R2 and one delta. The bench model also shows a period-3 window at 746 ohm.
