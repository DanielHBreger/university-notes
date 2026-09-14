# `chua/gallery.py`

N2 of the plan: four clean measured records, one per regime (limit cycle, period 2, single scroll, double scroll), for the report. A fifth panel shows the large outer cycle because it is the attractor the hysteresis item is about. Output: `gallery.png` beside the script.

## `REGIMES`

The default `(label, Rpot)` pairs. The resistances sit inside each regime as read off the bifurcation diagram and `cascade_periods_forward.txt`: 806 ohm for period 1, 760 for period 2 (well below R1 = 774.6 so the two loops are clearly apart), 700 for the single scroll, 550 for the double scroll, 190 for the large cycle.

## `nearest(folder, rpot, max_residual=5.0)`

Reads `<folder>_rpot.csv` through `sweeplib.load_rpot`, keeps records whose divider residual is under 5 %, and returns the name and value nearest the requested resistance. The forward sweep has a record every 1 to 2 ohm, so the nearest one is within about 1 ohm of the request.

## `main()`

For each resistance: load the record with `read_scope`, draw the first `--xy-samples` (150 000) samples as V2 against V1 in the top row, V2 vertical and V1 horizontal as the plan's XY display, and the first `--window` (4 ms, about twelve windings) of V1(t) in the bottom row. The title of each panel names the record and its fitted Rpot, and the same information is printed with the V1 range.

`--r` and `--labels` change the set; `--sweep back` takes the records from the back sweep.
