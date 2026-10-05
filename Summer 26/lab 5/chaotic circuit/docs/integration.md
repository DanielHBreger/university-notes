# `chua/integration.py`

The compiled integrators. Everything numerically heavy in `simulate.py` and `feigenbaum.py` runs here under numba's `@njit`, which turns the Python into machine code on first call (cached on disk). Every circuit value is an argument, never a frozen global, so the same kernel serves any parameter set.

Two models share the file: the **ideal** model with three states and a static piecewise-linear element, and the **bench** model with six continuous states and a three-number turning-point tracker, where the element is Kennedy's two-op-amp diode integrated with its own dynamics, the inductor's core law acts on the current's excursion from its running mean, and node 1's capacitance grows with the swing of v1 since its last turning point.

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
    g = pwl_eval(v1, vk, ik)
    ir = (v2 - v1) / rt
    return (ir - g) / c1, (-ir - i) / c2, (v2 - rl * i) / ind
```

The three circuit equations, with the element as a table lookup. `pwl_eval` is `np.interp` written out for a handful of knots (equal to 3e-17 A, seven times faster).

- `step`: one classical RK4 step, written out in scalars (numba is fastest on scalars).
- `trajectory(y, rt, rl, c1, c2, ind, vk, ik, dt, nsteps, nskip)`: the states at steps `nskip..nsteps`, shape `(nsteps − nskip + 1, 3)`. Used for portraits and time series.
- `hold(y, ..., dt, nsettle, ncollect)`: integrates `nsettle` steps without recording, then `ncollect` steps recording every local maximum of v1, found by the three-point test `p1 > p2 and p1 >= v1` on consecutive samples. Returns the final state and the maxima. This is the building block of a continuation sweep: the returned state is the start of the next resistance.

A test integrates a linear element and compares the result with the matrix exponential of the exact linear system to 1e-8.

## Bench model

State `(v1, v2, iL, voA, voB, im, v_rev, v_ext, rising)` where `voA`, `voB` are the outputs of the two op-amp stages, `im` the running mean of iL, and the last three the turning-point tracker of v1 (`N_BENCH` = 9). Parameters travel as one vector `P` in the order `BENCH_FIELDS = ('c1','c2','l0','nu','r0','rho','RA','AA','RB','AB','vpA','vnA','vpB','vnB','sr','tauA','tauB','tauM','c1_a0','c1_a1','c1_xk')`; `bench_params(...)` builds it (`tauM` defaults to `TAU_M` = 3 ms, the node-1 law to off: `c1_a0 = c1_a1 = 0`).

The equations:

```
i_A = (v1 - voA)/RA,  i_B = (v1 - voB)/RB,  i_NR = i_A + i_B
voX' = clip((AX v1 - voX)/tauX, -SR, +SR),  voX held inside [VnX, VpX]
C1(x) v1' = (v2 - v1)/Rt - i_NR,   C1(x) = c1 + a0 min(x, xk) + a1 max(x - xk, 0),   x = |v1 - v_rev|
C2 v2' = (v1 - v2)/Rt - iL
(L0 + nu |d|) iL' = v2 - r0 iL - rho |d| d,   d = iL - im
im' = (iL - im) / tauM
```

Each op-amp stage is a non-inverting amplifier of gain AX feeding node 1 through RX. Its output follows AX·v1 with a first-order lag tauX and a slew-rate limit SR, and stops at its rails. The inductor follows the Rayleigh law of a ferromagnetic core: inductance and core loss both grow linearly with the current's excursion from its mean over a few periods, which is how `identify.py` measures them (each record's flux loop about its mean current). A steady current sees L0 and the winding resistance r0 only. Before 27 September the law acted on |iL| itself, which charged the -1.7 mA of DC at the negative equilibrium with 14 ohm of core loss the small cycles do not show and put the model's real Hopf point at 846 ohm; a test now checks the steady-current case.

The capacitor at node 1 follows the Rayleigh law of a ferroelectric (class-2 ceramic) dielectric: its incremental capacitance grows with the swing of v1 since v1 last turned, steeply over the first `xk` and by `a1` per volt beyond. `identify.py` measures all four numbers (`node1_law`). The tracker is updated after every RK4 step by `_turn`: a turning point counts once v1 has come back `TURN_V` = 20 mV from its extreme, and the remembered point fades towards v1 over `TAU_R` = 10 ms. The fading is there so that the slow drift of the equilibrium while the knob turns (no turning points at all above the onset) is not counted as swing; it is far longer than an oscillation period (0.3–0.7 ms), so it changes the swing within a cycle by about 1 %. `swing(v, dt)` computes the same x along a sampled signal, which is how `identify.py` measures it; a test checks the two agree. With the law off (`a0 = a1 = 0`) the other six states evolve bit for bit as before the tracker was added.

### `c1_of_swing(x, P)`

Node 1's capacitance at a swing x, as the kernel uses it.

`GBW` = 3 MHz, the TL082's datasheet gain-bandwidth, lives here too, so that `identify.py` (stage A's lag share of node 1) and `simulate.py` (the lags) use one value: a stage of gain A lags by A/(2π·GBW).

### `static_pwl(P)`

With the op-amps infinitely fast the element is a five-segment law. Where neither stage is saturated the current is (v1 − AA v1)/RA + (v1 − AB v1)/RB, so

```
Ga = (1 - AA)/RA + (1 - AB)/RB           inner segment
Gb = (1 - AA)/RA + 1/RB                  stage B saturated: its current is (v1 - VpB)/RB, slope 1/RB
Gc = 1/RA + 1/RB                         both saturated
```

Stage B saturates first, at v1 = VpB/AB and VnB/AB (the inner breakpoints); stage A at VpA/AA and VnA/AA (the outer ones). A test checks that `kennedy_from_pwl` in `identify.py` and this function are inverses.

### `bench_state(v1, v2, iL, P)`

The full state vector with each op-amp output placed where it would sit statically, clamped to its rails, the running mean of the current equal to the current, and the tracker at v1 (no swing yet). Used to start a run from a three-state point.

### `_opamp_rate(target, vo, tau, sr, vn, vp)`

```python
d = (target - vo) / tau                 # first-order lag
d = clip(d, -sr, sr)                    # slew limit
if (vo >= vp and d > 0) or (vo <= vn and d < 0): d = 0.0   # railed and pushed outward: held
```

The output can leave a rail only inward. A test checks the held and the slew-limited cases.

### `bench_derivative`, `bench_step`, `bench_trajectory`, `bench_hold`

The bench equivalents of the ideal kernels. `bench_derivative(v1, v2, iL, voA, voB, im, v_rev, rt, P)` returns the six continuous rates, with v_rev held through the step. In `bench_step` the op-amp outputs are clamped to their rails at every intermediate RK4 stage and again after the update, so no stage of the step ever evaluates an output beyond a rail, and the tracker is updated from the new v1.

### The step size

An explicit RK4 step on `vo' = −vo/tau` is stable only for dt < 2.78·tau. Stage A has tau = 0.059 µs, so the bench model is integrated with dt = 0.1 µs (`BENCH_DT` in `simulate.py`). A 0.2 µs step was tried and produced spurious 22 µs "periods" from the unstable lag equation.

A test integrates the bench kernel with the rails far away and the Rayleigh terms zero, where it is a six-state linear system, and checks it against the matrix exponential to 1e-6. `verify_models.py` checks both kernels against scipy's DOP853 at rtol 1e-10 (for the bench model with the tracker updated after every accepted DOP853 step, as the kernel does).
