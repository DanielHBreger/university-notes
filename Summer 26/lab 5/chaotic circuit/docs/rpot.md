# `chua/rpot.py` and `chua/batch_rpot.py`

The potentiometer was never read off the dial. Every record carries its own measurement of Rpot, taken from the third scope channel without opening the circuit.

## The idea

R0 and the potentiometer are in series between node 1 (v1, CH1) and node 2 (v2, CH2), with CH3 on the junction between them. No current flows into the probe, so the junction sits on a plain voltage divider:

```
v_m = A v1 + B v2 + C
```

where in the ideal case A = Rpot/(R0 + Rpot) and B = R0/(R0 + Rpot), so B/A = R0/Rpot and

```
Rpot = R0 · (B/A) · (g2/g1)
```

`g2/g1` is the ratio of the two end channels' gains. Using the ratio B/A instead of A alone has two consequences that matter:

- the gain of the midpoint channel cancels completely (it scales A and B together);
- only the relative gain of CH1 and CH2 survives. Equal V/div settings do not guarantee equal gains. The identity A + B ≈ 1 is a consistency check, not a calibration of their ratio. See `uncertainty/REPORT.md` for the error budget and unresolved calibration terms.

The offset C absorbs any DC offset in the channels. Fitting over the whole record, with v1 and v2 swinging independently, is what makes A and B separable.

## `rpot(path, r0, g21, c1='CH1', c2='CH2', mid='CH3')`

Validates the arguments, reads the three channels through `read_scope` (at least ten samples) and hands the array to `fit_divider`. Returns `(Rpot, residual fraction)`.

## `fit_divider(d, r0=R0, g21=1.0)`

The least-squares fit, written to fail loudly on data that cannot identify the ratio.

```python
v1, v2, vm = d.T
M = np.column_stack([v1 - v1.mean(), v2 - v2.mean()])
```

Centring the regressors and the target removes the offset C from the problem without fitting it.

```python
scales = M.std(axis=0)
if np.any(scales <= 1e-12) or np.ptp(vm) <= 1e-12:
    raise ValueError('flat channels cannot identify the divider')
Z = M / scales
```

Each regressor is scaled to unit standard deviation before the solve. v2 swings a tenth of v1, so without this the normal equations would be badly conditioned and the rank test below would be meaningless.

```python
coef, _, rank, singular = np.linalg.lstsq(Z, vm - vm.mean(), rcond=None)
if rank < 2 or singular[-1] / singular[0] < 1e-4:
    raise ValueError('end channels are collinear; ...')
```

If v1 and v2 are proportional (a DC record, or a pure sinusoid seen on both channels in phase) there is no way to split the midpoint between them. The ratio of singular values is the check.

```python
A, B = coef / scales
if A <= 1e-10 or B < 0:
    raise ValueError('non-physical divider coefficients; check channels/clipping')
```

Undoing the scaling gives the real coefficients. A must be positive (Rpot > 0), B must not be negative. A clipped channel typically produces a negative coefficient here.

```python
resid = float(np.sqrt(np.mean((vm - vm.mean() - Z @ coef) ** 2)))
span = float(vm.max() - vm.min())
return r0 * (B / A) * g21, resid / max(span, 1e-12)
```

The residual is the RMS misfit of the midpoint channel as a fraction of its span. Above 5 % (`CLEAN_DIVIDER_MAX_PCT`) the record is flagged; the plotting scripts drop such records because a wrong Rpot places a good slice at the wrong place on the axis.

## `main()` of `rpot.py`

Accepts file names or glob patterns, prints one line per record with the value, the residual and the `CHECK` warnings for a dirty divider or a value outside 0 to 1000 ohm.

## `batch_rpot.py`

Runs `rpot` over every scope CSV in a folder and writes `<folder>_rpot.csv` with the columns `filename, rpot_ohm, residual_pct, status`. Points to note:

- `files = [f for f in files if os.path.abspath(f) != os.path.abspath(out)]` makes sure the output of a previous run is never read as an input, in case someone points `-o` inside the folder.
- A record that raises is kept in the table with blank numbers and `status = error: ...`, so the row count always equals the record count and the failure is visible.
- `status` is `ok`, `CHECK: not a clean divider` (residual above 5 %), or `CHECK: outside 0-1000 ohm potentiometer range`. `load_rpot` in `sweeplib.py` keeps any row with finite numbers; the consumers then apply their own residual and range cuts, so the status column is documentation rather than a filter.
- The two sweep sidecars that ship with the project (`forward_rpot.csv`, `back_rpot.csv`) were made with R0 = 992 and g21 = 1.
