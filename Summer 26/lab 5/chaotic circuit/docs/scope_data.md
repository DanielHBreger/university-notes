# `chua/scope_data.py`

The one place a scope CSV is turned into arrays. Every other script goes through `read_scope`, so a bad record is refused the same way everywhere.

## Why it is strict

A scope record is a time series. If a row with a bad value were dropped, the samples on either side of it would be joined and every quantity that depends on the spacing of samples (return times, derivatives, periods) would be wrong at that point without any warning. The reader therefore rejects the whole file instead of cleaning it.

## Constants

```python
R0 = 992.0
RPOT_MAX = 1000.0
```

`R0` is the fixed resistor in series with the potentiometer, as measured on the bench. The plan says 990 ohm. Every committed `_rpot.csv` was computed with 992, and every Rpot in the project scales as R0 if this is changed (Rpot = R0 · B/A, see [rpot.md](rpot.md)). `RPOT_MAX` is the dial's end stop, used as a sanity limit on fitted values.

## `read_scope(path, channels=('CH1', 'CH2'), min_samples=3)`

Returns `(t, data)` where `t` is the time column as float64 and `data` has one column per requested channel, in the order requested.

Line by line:

- `head = pd.read_csv(path, nrows=0)` reads only the header. The column names look like `Time(s)`, `CH1(V)`; the unit in brackets is stripped so the caller can ask for `'CH1'`.
- The first column must be time; channel names must be unique; every requested channel must exist and must not be the time column.
- `cols = [head.columns[0]] + [...]` builds the list of real column names in the requested order, and `pd.read_csv(path, usecols=cols, dtype=np.float64)[cols]` reads only those columns. Reading two columns instead of four is a fifth of the parsing time on a 28 MB record, which is why `lorenz_map.maxima` asks for one channel only.
- `len(data) < min_samples` and `not np.isfinite(data).all()` reject short and broken files.
- `validate_time(data[:, 0])` is the timing check below.

## `validate_time(t)`

Returns the sample interval `dt` after checking that time is finite, strictly increasing and uniform.

The tolerance line is the subtle one:

```python
tolerance = max(0.005 * dt, 2e-6 * float(np.max(np.abs(t))), 1e-15)
```

The scope writes timestamps as `%.6e`, six significant digits. A record that starts at t = 1 s with 1 µs steps has steps that are exact to only a few parts in 10⁶ of the absolute time, so a fixed relative tolerance on `dt` would reject perfectly good files. The tolerance is therefore the larger of half a percent of the step and two parts per million of the largest timestamp.
