# Code explanations

One file per script. Each file walks the script function by function, in the order the data flows through it, and stops on individual lines where the code is doing something that is not obvious. Read them in the order below and the whole pipeline follows.

| order | file | script | what it settles |
|---|---|---|---|
| 1 | [scope_data.md](scope_data.md) | `chua/scope_data.py` | how a record is read and rejected, and the value of R0 |
| 2 | [sweeplib.md](sweeplib.md) | `chua/sweeplib.py` | where files live and how artefacts are named |
| 3 | [rpot.md](rpot.md) | `chua/rpot.py`, `chua/batch_rpot.py` | the potentiometer value of every record |
| 4 | [diode_analysis.md](diode_analysis.md) | `diode_analysis.py`, `find_breakpoints.py`, `plot_current.py` | M1 and M2: the five-segment fit of the V–I record |
| 5 | [lorenz_map.md](lorenz_map.md) | `chua/lorenz_map.py` | the maxima of a record, the return map, and which records count |
| 6 | [bifurcation.md](bifurcation.md) | `chua/bifurcation.py` | the measured bifurcation diagram and hysteresis overlay |
| 7 | [lyapunov.md](lyapunov.md) | `chua/lyapunov.py` | the return-map Lyapunov exponent and its gates |
| 8 | [rosenstein.md](rosenstein.md) | `chua/rosenstein.py` | the direct time-series exponent |
| 9 | [benchmark_lyapunov.md](benchmark_lyapunov.md) | `chua/benchmark_lyapunov.py` | how far either exponent can be trusted |
| 10 | [cascade_periods.md](cascade_periods.md) | `chua/cascade_periods.py` | the measured period of each record and R1, R2, R3 |
| 11 | [identify.md](identify.md) | `chua/identify.py` | the circuit as it actually was: C1, C2, inductor law, diode |
| 12 | [integration.md](integration.md) | `chua/integration.py` | the compiled integrators of both models |
| 13 | [simulate.md](simulate.md) | `chua/simulate.py` | the three models, their sweeps, portraits and transitions |
| 14 | [feigenbaum.md](feigenbaum.md) | `chua/feigenbaum.py` | the model's doubling points and Feigenbaum delta |
| 15 | [tests.md](tests.md) | `tests/test_analysis.py` | what each test pins down |
| 16 | [s3014a.md](s3014a.md) | `chua/s3014a.py` | the scope program used at the bench |

## The circuit, once

Every script uses the same three equations. Node 1 is the capacitor C1 with the nonlinear element across it, node 2 is C2 with the inductor across it, and the two nodes are joined by the resistance Rt = R0 + Rpot.

```
C1 dv1/dt = (v2 - v1)/Rt - i_NR(v1)
C2 dv2/dt = (v1 - v2)/Rt - iL
L  diL/dt = v2 - rL iL
```

Channel convention in the oscillator records: CH1 = v1, CH2 = v2, CH3 = the junction between R0 and the potentiometer. In the V–I record (`trace1.csv`) the convention is different: CH1 is the source side of a 216 ohm shunt and CH2 is the voltage across the element.

## Naming conventions

- **Records**: `forward/traceN.csv` and `back/traceN.csv`. N counts up in acquisition order. In `forward/` the knob was turned down as N rose, in `back/` it was turned up.
- **Sidecars**: anything derived from a sweep folder is written beside it and named `<folder><suffix>`, by `sweeplib.sibling`. So `forward/` produces `forward_rpot.csv`, `forward_bifurcation.png`, `forward_bifurcation_points.csv`, `forward_lorenz/`, `forward_lyapunov.csv`, `forward_lyapunov.png`, `forward_lyapunov_each/`. The overlay of both sweeps is `forward_vs_back_hysteresis.png` with its `_points.csv`.
- **`_points.csv`**: the numbers behind the figure of the same name.
- **Simulation outputs**: `simulated_<figure><tag>.png` and the log `simulate<tag>.txt`. The tag names the model run: empty for the plan's model with C1 = 11.5 nF, `_nominal` for C1 = 10 nF, `_sym` for the symmetrised element, `_static` for the identified circuit with constant components, `_bench` for the identified circuit in full.
- **`identified.*`**: outputs of `identify.py`. The JSON is what `simulate.py` reads.
- **`feigenbaum_<model>.txt`**, **`cascade_periods_<sweep>.txt`**, **`benchmark_lyapunov.txt`**: the kept printout of the script of that name.
- **Legacy record names** such as `chaos - 716.9 ohm.csv` are still recognised by `lyapunov.py` and `simulate.py`, which read the resistance out of the name.

## How the pieces connect

```
trace1.csv ──diode_analysis──> diode_fit.json ──────────────────────┐
                                                                    │
forward/, back/ ──batch_rpot──> <sweep>_rpot.csv                    │
      │                              │                              │
      │      ┌───────────────────────┤                              │
      ▼      ▼                       ▼                              ▼
   lorenz_map.collect ──> bifurcation ──> <sweep>_bifurcation_points.csv    simulate --model ideal
      │                        │                     │
      ▼                        ▼                     ▼
   lyapunov (+rosenstein)   hysteresis overlay   cascade_periods
      │
      ▼
   benchmark_lyapunov (calibration on a simulated record)

forward/, back/ + <sweep>_rpot.csv ──identify──> identified.json ──> simulate --model static / bench
                                                                 └──> feigenbaum --model ...
```
