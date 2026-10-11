# `chua/benchmark_lyapunov.py`

Calibrates the two Lyapunov estimators on a simulated record whose true exponent is known. It backs the statements in `lyapunov.py` and `rosenstein.py` about how far each estimate can be trusted. The kept printout is `benchmark_lyapunov.txt` (the default run, then Matsumoto's parameters appended).

## The simulated circuit

Chua's dimensionless form:

```
x' = alpha (y - x - f(x)),   y' = x - y + z,   z' = -beta y
f(x) = m1 x + (m0 - m1)(|x + 1| - |x - 1|)/2
```

with alpha = C2/C1 and beta = C2 R²/L from the nominal components (`C1, C2, L = 10e-9, 100e-9, 18e-3`) and R = R0 + Rpot. x is in units of the breakpoint voltage (about 1 V), so the simulated v1 spans a few volts like the real one; time is scaled by C2·R so one sample is 1 µs as on the bench. The canonical slopes m0 = −8/7, m1 = −5/7 are the defaults; `--alpha 9 --beta 14.286` is Matsumoto's double scroll.

## `simulate(alpha, beta, m0, m1, n, dtau, seed=0, transient=50000)`

Fixed-step RK4 on a six-dimensional system: the three states plus one tangent vector. The tangent vector obeys the variational equation (the Jacobian applied to it, with `fp = m1 if |x| > 1 else m0` as the local slope of f):

```python
return np.array([alpha * (y - x - fx(x)), x - y + z, -beta * y,
                 -alpha * (1 + fp) * vx + alpha * vy, vx - vy + vz, -beta * vy])
```

After each step the tangent vector is renormalised and the log of its growth accumulated:

```python
nv = np.linalg.norm(s[3:]); lnsum += np.log(nv); s[3:] /= nv
```

`lnsum / (n · dtau)` is the largest Lyapunov exponent per unit dimensionless time, the standard calculation. A test runs this on a parameter set where the attractor is a fixed point and checks the exponent equals the largest real part of the Jacobian's eigenvalues.

## `main()`

- Builds alpha, beta and the time unit, simulates, and converts the true exponent to 1/s.
- Quantises v1 and v2 the way the scope does (`--lsb 0.04 0.004` V) after adding a little Gaussian noise, and writes a scope-format CSV with a fake CH3 into a temporary folder.
- Pushes that file through exactly the measurement pipeline: `lorenz_map.maxima` → `lyapunov.lyapunov` for the map estimate, `rosenstein.rosenstein` for the direct one.
- Prints the true value, each estimate, and each as a multiple of the truth.

What the table shows, and why it matters: the map is only approximately a one-dimensional function of M(n), and where it folds (the lobe switch) the projection has a slope the flow does not; the direct estimate needs the channels smoothed before the quantisation lets it see the divergence. The comparison is case-specific. It does not establish a universal correction factor, and neither script claims one.
