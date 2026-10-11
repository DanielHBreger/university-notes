# `chua/sweeplib.py`

File plumbing shared by the sweep scripts. Nothing here touches the physics. It is kept to the standard library so `batch_rpot.py` can run without numpy or matplotlib installed.

## The one convention

```python
def sibling(folder, suffix):
    folder = os.path.normpath(folder)
    return os.path.join(os.path.dirname(folder), os.path.basename(folder) + suffix)
```

Every artefact derived from a sweep folder is written beside the folder and named after it: `forward/` gives `forward_rpot.csv`, `forward_lorenz/`, `forward_bifurcation.png`. `normpath` strips a trailing slash so `"forward/"` and `"forward"` give the same name (a test checks this). Because the producer (`batch_rpot.py`) and the consumers (`load_rpot`, and through it every plotting script) call the same function, they cannot disagree about where the sidecar is.

## `CLEAN_DIVIDER_MAX_PCT = 5.0`

The residual, in percent of the midpoint channel's span, above which a divider fit is called "not a clean divider". `rpot.py` prints the warning, `batch_rpot.py` writes the verdict into the `status` column, and the plotting scripts re-derive it from the numeric residual. One constant, three users.

## Finding records

- `natural_key(name)` splits a name at digit runs and turns the digits into integers, so `trace2` sorts before `trace10`. Acquisition order is the natural order.
- `is_scope_csv(path)` reads one line and checks that the first header cell is `Time(s)`. This keeps generated summaries (`lorenz_summary.csv`, an `_rpot.csv` reached by a recursive walk) out of the input list.
- `list_csvs(folder, recursive=False, full=False)` lists the scope CSVs in natural order, as names or full paths. Non-recursive listing does not need to exclude sidecars because they live outside the folder.

## Choosing folders

- `script_dir()` is where the running script lives; the sweep folders sit beside it, so dialogs and prompts open there.
- `pick_folder` opens a Tk directory chooser and falls back to `input()` when there is no display. In the prompt fallback, a bare name is resolved beside the script, matching what the dialog would have done.
- `pick_folders` keeps asking for another folder until Cancel, starting each dialog beside the folder just chosen because sweeps of one run live together.
- `resolve_folders(named, multi=True)` takes named folders when given, otherwise asks, and exits if none survive. `multi=False` asks once.

## Reading the sidecar back

```python
def load_rpot(folder, required=False):
```

Returns `{filename: (Rpot, residual %)}` from `<folder>_rpot.csv`. Rows whose value or residual is blank or non-finite are skipped, so records that failed the fit drop out here and do not have to be filtered by every caller. With `required=True` the function exits with the exact `batch_rpot.py` command to run when the sidecar is missing. The scripts that put Rpot on an axis (bifurcation, Lyapunov) require it; `lorenz_map.py` only labels panels with it and carries on without.

## Plot helpers

- `folder_labels(folders)` returns basenames, prefixed with the parent folder only when two inputs share a basename (two `forward` folders from different sets).
- `sweep_colors(n)` gives one colour per sweep from `tab10`, or `tab20` above ten. The matplotlib import is inside the function so the module stays import-free for `batch_rpot.py`.
