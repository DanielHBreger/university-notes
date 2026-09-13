# `chua/integration.py`

The compiled integrators. Everything numerically heavy in `simulate.py` and `feigenbaum.py` runs here under numba's `@njit`, which turns the Python into machine code on first call (cached on disk). Every circuit value is an argument, never a frozen global, so the same kernel serves any parameter set.

Two models share the file: the **ideal** model with three states and a static piecewise-linear element, and the **bench** model with five states, where the element is Kennedy's two-op-amp diode integrated with its own dynamics.

## `pwl_knots(bp, G, extend=40.0)`

Turns four breakpoints and five slopes into interpolation knots `(v, i)` of the continuous five-segment law through the origin. The curve is built outward from v = 0 in both directions:

```python
j0 = int(np.searchsorted(vs, 0.0)) - 1     # the segment containing v = 0
cur = 0.0
for j in range(j0, len(vs) - 1):           # rightward: i accumulates G_j * width
    lo = 0.0 if j == j0 else vs[j]
    cur += G[j] * (vs[j + 1] - lo)
    ik[j + 1] = cur
```

and the mirror loop leftward. Starting from the origin guarantees g(0) = 0 exactly and continuity at every breakpoint. The knots are extended to ±40 V so `np.interp`, which clamps outside its range, never has to extrapolate inside any simulation. A test checks the slopes between knots equal G and that g(0) = 0.

## Ideal model

```python
@njit(cache=True)
def derivative(v1, v2, i, rt, rl, c1, c2, ind, vk, ik):
    g = np.interp(v1, vk, ik)
    ir = (v2 - v1) / rt
    return (ir - g) / c1, (-ir - i) / c2, (v2 - rl * i) / ind
```

The three circuit equations, with the element as a table lookup.

- `step`: one classical RK4 step, written out in scalars (numba is fastest on scalars).
- `trajectory(y, rt, rl, c1, c2, ind, vk, ik, dt, nsteps, nskip)`: the states at steps `nskip..nsteps`, shape `(nsteps − nskip + 1, 3)`. Used for portraits and time series.
- `hold(y, ..., dt, nsettle, ncollect)`: integrates `nsettle` steps without recording, then `ncollect` steps recording every local maximum of v1, found by the three-point test `p1 > p2 and p1 >= v1` on consecutive samples. Returns the final state and the maxima. This is the building block of a continuation sweep: the returned state is the start of the next resistance.

A test integrates a linear element and compares the result with the matrix exponential of the exact linear system to 1e-8.

## Bench model

State `(v1, v2, iL, voA, voB)` where `voA`, `voB` are the outputs of the two op-amp stages. Parameters travel as one vector `P` in the order `BENCH_FIELDS = ('c1','c2','l0','nu','r0','rho','RA','AA','RB','AB','vpA','vnA','vpB','vnB','sr','tauA','tauB')`; `bench_params(...)` builds it.

The equations:

```
i_A = (v1 - voA)/RA,  i_B = (v1 - voB)/RB,  i_NR = i_A + i_B
voX' = clip((AX v1 - voX)/tauX, -SR, +SR),  voX held inside [VnX, VpX]
C1 v1' = (v2 - v1)/Rt - i_NR
C2 v2' = (v1 - v2)/Rt - iL
(L0 + nu |iL|) iL' = v2 - (r0 + rho |iL|) iL
```

Each op-amp stage is a non-inverting amplifier of gain AX feeding node 1 through RX. Its output follows AX·v1 with a first-order lag tauX and a slew-rate limit SR, and stops at its rails. The inductor follows the Rayleigh law of a ferromagnetic core: inductance and loss both grow linearly with |current|.

### `static_pwl(P)`

With the op-amps infinitely fast the element is a five-segment law. Where neither stage is saturated the current is (v1 − AA v1)/RA + (v1 − AB v1)/RB, so

```
Ga = (1 - AA)/RA + (1 - AB)/RB           inner segment
Gb = (1 - AA)/RA + 1/RB                  stage B saturated: its current is (v1 - VpB)/RB, slope 1/RB
Gc = 1/RA + 1/RB                         both saturated
```

Stage B saturates first, at v1 = VpB/AB and VnB/AB (the inner breakpoints); stage A at VpA/AA and VnA/AA (the outer ones). A test checks that `kennedy_from_pwl` in `identify.py` and this function are inverses.

### `bench_state(v1, v2, iL, P)`

The full five-state vector with each op-amp output placed where it would sit statically, clamped to its rails. Used to start a run from a three-state point.

### `_opamp_rate(target, vo, tau, sr, vn, vp)`

```python
d = (target - vo) / tau                 # first-order lag
d = clip(d, -sr, sr)                    # slew limit
if (vo >= vp and d > 0) or (vo <= vn and d < 0): d = 0.0   # railed and pushed outward: held
```

The output can leave a rail only inward. A test checks the held and the slew-limited cases.

### `bench_derivative`, `bench_step`, `bench_trajectory`, `bench_hold`

The five-state equivalents of the ideal kernels. In `bench_step` the op-amp outputs are clamped to their rails at every intermediate RK4 stage and again after the update, so no stage of the step ever evaluates an output beyond a rail.

### The step size

An explicit RK4 step on `vo' = −vo/tau` is stable only for dt < 2.78·tau. Stage A has tau = 0.05 µs, so the bench model is integrated with dt = 0.1 µs (`BENCH_DT` in `simulate.py`). A 0.2 µs step was tried and produced spurious 22 µs "periods" from the unstable lag equation.

A test integrates the bench kernel with the rails far away and the Rayleigh terms zero, where it is a five-state linear system, and checks it against the matrix exponential to 1e-6.
