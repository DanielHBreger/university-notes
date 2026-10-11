# `chua/bifurcation.py`

The measured bifurcation diagram: every local maximum of v1 in every admitted record, plotted at that record's Rpot. A period-1 cycle is one point on its vertical line, period-2 a pair, chaos a filled band. Given two sweep folders it overlays them in different colours, which is how the hysteresis figure `forward_vs_back_hysteresis.png` was made.

Outputs: `<first folder>_bifurcation.png` and the same name with `_points.csv`, one row per maximum with `sweep, filename, rpot_ohm, rpot_residual_pct, max_v`. That CSV is what `cascade_periods.py` and `simulate.py` read back, so the diagram, the measured periods and the model comparison all rest on the same maxima.

## `add_record_args(p)`

The command-line options that decide which records count (`--ch`, `--prominence`, `--period-samples`, `--max-residual`, `--max-clip`, `--max-files`). Shared with `lyapunov.py` by import, so the two scripts admit the same records for the same flags.

## `report(label, records, dropped)`

Prints what was admitted, the range of Rpot, the records admitted with a residual over the 5 % convention (possible only if `--max-residual` was raised), and every dropped record with its reason. Nothing is silently discarded.

## `main()`

- `folders = resolve_folders(a.folders)`: named folders or the repeated picker.
- Output name from the first folder via `sibling`; the points CSV is the PNG name with `_points.csv`.
- For each folder: `collect(...)` from `lorenz_map.py`, then one scatter for the whole sweep:

```python
xs = np.concatenate([np.full(len(r.M), r.rpot) for r in records])
ys = np.concatenate([r.M for r in records])
ax.scatter(xs, ys, s=a.size, alpha=a.alpha, lw=0, color=colors[k], ..., rasterized=True)
```

One scatter call per sweep instead of one per record, because matplotlib's per-call overhead would dominate with 400 records. `rasterized=True` keeps the PNG small at 450 dpi.

- The rows for the CSV are collected in the same loop.
- `--descending` inverts the x axis so the forward sweep (knob turned down) reads left to right; `--xlim`, `--ylim` crop.
- The legend handles are set to full alpha and a visible size, since the plotted markers are tiny and translucent.

The measured diagram carries no Rpot uncertainty bars; the residual column in the CSV is the per-record quality of the Rpot fit.
