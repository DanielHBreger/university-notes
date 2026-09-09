# Chua oscillator analysis

Start with [the verification report](VERIFICATION_REPORT.md) and [the figure gallery](verification/figures.html). The report distinguishes corrected calculations from experimental conclusions that the existing records cannot establish.

## Reproduce the analysis

Run from this folder using Python 3.11 or newer. A local `.venv` has already been prepared. For a fresh environment:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

Core commands:

```powershell
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe find_breakpoints.py
.\.venv\Scripts\python.exe chua/rpot.py "chua/doublescroll 2 - 656.7ohm.csv"
.\.venv\Scripts\python.exe chua/batch_rpot.py "chua/sweep forward"
.\.venv\Scripts\python.exe chua/lorenz_map.py "chua/sweep forward"
.\.venv\Scripts\python.exe chua/bifurcation.py "chua/sweep forward" "chua/sweep back" -o chua/hysteresis_bifurcation.png
.\.venv\Scripts\python.exe chua/lyapunov.py "chua/sweep forward" --rosenstein
.\.venv\Scripts\python.exe chua/simulate.py
.\.venv\Scripts\python.exe chua/simulate.py --i0 auto --tag _corrected
```

Repeat the sweep commands for `chua/sweep back`, `chua/set 2/sweep forward`, and `chua/set 2/sweep back`. Direct Lyapunov estimates are optional and substantially slower; the supplied set 2 summaries contain map estimates only. `plot_current.py` and `find_breakpoints.py` both regenerate the fitted V–I results through the same analysis implementation.

The verification helpers reproduce the full audit and figure set without repeatedly parsing every large recording:

```powershell
.\.venv\Scripts\python.exe verification/audit_data.py
.\.venv\Scripts\python.exe verification/rebuild_figures.py
.\.venv\Scripts\python.exe verification/required_results.py
.\.venv\Scripts\python.exe verification/selected_portraits.py
.\.venv\Scripts\python.exe verification/validate_simulation.py
.\.venv\Scripts\python.exe verification/finish_results.py
.\.venv\Scripts\python.exe verification/redraw_simulation.py
.\.venv\Scripts\python.exe verification/check_outputs.py
.\.venv\Scripts\python.exe verification/make_gallery.py
```

Run the two simulation commands before `finish_results.py` and `redraw_simulation.py`; they use the saved simulation trajectories and maxima. The audit caches are ignored by Git and can be regenerated from the original scope CSVs. Rebuild the audit after changing measurements or peak-finding settings. The audit records input hashes, and `check_outputs.py` detects changes since that audit.

## Conventions

- Oscillator recordings: CH1 = V1, CH2 = V2, CH3 = divider midpoint V3.
- Root `trace1.csv` is the separate V–I acquisition: CH1 = source voltage VV and CH2 = element voltage VI. Current is `(CH1 - CH2)/216`, in amperes. The fitting CLI permits explicit column overrides.
- `Rpot` is the potentiometer alone. The oscillator coupling resistance is `990 Ω + Rpot`. The 216 Ω V–I shunt is switched out in CHAOS mode.
- Generated `_rpot.csv`, `_bifurcation_points.csv`, Lorenz and Lyapunov summaries are derived results. Scope recordings are not overwritten.
- Simulation uses the fitted element voltage curve in `diode_fit.json`. `--i0 auto` removes its inner current intercept as an alternative calibration assumption. The historical `_corrected` suffix does **not** mean that calibration has been verified.
- `s3014a.py` controls real hardware when launched. Its transfer/decimation logic was reviewed and tested with mocks; no connected oscilloscope was available for an end-to-end acquisition test.
