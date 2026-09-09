# Chua oscillator, lab 5

Branch `lab-5-review`: the analysis code of `lab-5` (commit c75206a) plus the parts of the
ChatGPT review that were worth keeping, after checking everything against
`Lab5_Chua_experiment_plan` (Summer 2026, v5). The complete ChatGPT revision is preserved on
`lab-5-chatgpt-raw`; its `simulate.py`, its scope-driver edits and its `VERIFICATION_REPORT.md`
were rejected (reasons at the end).

## Data

- `chua/set 4/`: down-sweep, 418 records, 891 -> 4 ohm. `chua/set 5/`: up-sweep, 80 records,
  3 -> 901 ohm, dense between 685 and 780 ohm. These two sets and the V-I record are what the
  report uses; `sweep forward/back`, `set 2/` and `set 3/` are earlier sessions.
- Oscillator CSVs: `Time(s),CH1(V),CH2(V),CH3(V)` with CH1 = V1, CH2 = V2, CH3 = V3 (node
  between R0 and the potentiometer). Records are 200 or 500 ms.
- `trace1.csv` at the root is the V-I record: CH1 = VV (source side of the 216 ohm shunt),
  CH2 = VI (across the element); current = (CH1 - CH2) / 216.

## Pipeline

Run from this folder with `.venv\Scripts\python.exe` (or any Python 3.11+ with `requirements.txt`).

```
python find_breakpoints.py                                  # M1/M2: five-segment fit against VI -> diode_fit.json, figures
python chua/batch_rpot.py "chua/set 4"                      # M3: Rpot per record -> chua/set 4_rpot.csv
python chua/lorenz_map.py "chua/set 4"                      # M4: return maps -> chua/set 4_lorenz/
python chua/bifurcation.py "chua/set 4" "chua/set 5"        # N3/N4: overlay -> chua/set 4_bifurcation.png + _points.csv
python chua/lyapunov.py "chua/set 4" --rosenstein --each    # M5: map and direct exponents
python chua/simulate.py                                     # calibrated model: C1 = 11.5 nF, offset removed, shunt-converted fit
python chua/feigenbaum.py                                   # N5 from the model
python chua/cascade_periods.py "set 4"                      # N5 from the bench (run inside chua/)
python -m pytest -q                                         # regression tests
```

## Conventions and open points

- R0: every committed `*_rpot.csv` sidecar was computed with R0 = 992 ohm, which is the value in
  `chua/scope_data.py`. The plan quotes 990 ohm measured, and `simulate.py`/`diode_analysis.py`
  use 990 for the dial conversion. Decide on one value; switching `scope_data.py` to 990 rescales
  every Rpot by -0.2 % and needs `batch_rpot.py` rerun for sets 3-5.
- Divider residual above 5 % drops a record (plan), as does a fit outside 0-1000 ohm.
- The V-I fit is made against the element voltage VI (plan M1). `simulate.py` carries its own
  copy of the earlier VV-basis fit and converts it with `--shunt 216`; the two agree to within the
  inner-segment uncertainty (Ga between -0.72 and -0.78 mS depending on the inner window, so the
  M2 lower bound for Rpot is 295-405 ohm, not one number).
- The vertical range was changed several times within sets 4 and 5, always on all three channels
  together. One visible effect: between set 4 trace329 and trace330 Rpot steps from 455.6 to
  466.9 ohm, so set 4 values below about 467 ohm carry a +2 % offset relative to those above.

## What came from the ChatGPT review

Kept: `chua/scope_data.py` (strict CSV reader), the divider, peak, bifurcation and Lyapunov
scripts (same numbers on sets 4-5 as before; 5 % residual cut, 0-1000 ohm filter, input
validation, Rosenstein fit R2), `diode_analysis.py` with thin `find_breakpoints.py` /
`plot_current.py` wrappers (fit against VI, rising/falling overlap check, dial conversion),
`requirements.txt`, `.gitignore`, `tests/`.

Rejected: its `simulate.py` (nominal C1 puts the Hopf point at about 1000 ohm, outside the dial,
against 890-905 ohm observed; the retained -0.41 mA intercept leaves a single stable equilibrium
at +6 V so the model never oscillates above 850 ohm; `sweep_starts`/`hopf_point` needed by
`feigenbaum.py` were deleted; transients cut to 20 ms), its `s3014a.py` edits (untested on the
instrument; the `lab-5` driver produced sets 3-5), its regenerated outputs for sets 1-2, and
`verification/` + `VERIFICATION_REPORT.md` (conclusions drawn from sets 1-2 and from its own
simulation defaults).

Small fixes applied on this branch to `simulate.py`: the time axis of stored states is now the
step they were taken at (was off by up to one step when `--skip` is not a multiple of `--dt`);
`equilibria()` no longer misses an exact zero on its grid (the origin once the offset is removed);
`load_measured_bifurcation()` skips sidecars without a `sweep` column instead of crashing.
