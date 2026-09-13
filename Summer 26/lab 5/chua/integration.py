"""
Compiled RK4 kernels for simulate.py. Every circuit value is an argument, never
a frozen global, so the same code serves any parameter set.

Ideal model (three states v1, v2, iL; the nonlinear element is a static
piecewise-linear law given as interpolation knots):

    C1 v1' = (v2 - v1)/Rt - g(v1)
    C2 v2' = (v1 - v2)/Rt - iL
    L  iL' = v2 - rL iL

Bench model (five states v1, v2, iL, voA, voB), identify.py's circuit:

    i_A = (v1 - voA)/RA,   i_B = (v1 - voB)/RB,   i_NR = i_A + i_B
    voX' = clip((AX v1 - voX)/tauX, -SR, +SR), voX held inside [VnX, VpX]
    C1 v1' = (v2 - v1)/Rt - i_NR
    C2 v2' = (v1 - v2)/Rt - iL
    (L0 + nu |iL|) iL' = v2 - (r0 + rho |iL|) iL

The two op-amp stages of Kennedy's diode saturate at their output rails
(VpA, VnA) and (VpB, VnB); with SR -> inf and tau -> 0 the element is the
five-segment law
    Ga = (1-AA)/RA + (1-AB)/RB,  Gb = (1-AA)/RA + 1/RB,  Gc = 1/RA + 1/RB,
    inner breakpoints VnB/AB, VpB/AB,  outer breakpoints VnA/AA, VpA/AA.
The inductor follows the Rayleigh law of a ferromagnetic core at low field:
inductance and loss both grow linearly with the current.

The explicit RK4 step must satisfy dt < 2.78 tau for the op-amp equations,
so the bench model is integrated with dt = 0.1 us (tauA = 0.05 us).
"""
import numpy as np
from numba import njit

# index of every bench-model parameter in the parameter vector P
BENCH_FIELDS = ('c1', 'c2', 'l0', 'nu', 'r0', 'rho', 'RA', 'AA', 'RB', 'AB',
                'vpA', 'vnA', 'vpB', 'vnB', 'sr', 'tauA', 'tauB')


def pwl_knots(bp, G, extend=40.0):
    """
    Knots (v, i) of the continuous five-segment law through the origin with
    breakpoints bp (4 values, V) and slopes G (5 values, S), extended to
    +-extend V so that np.interp never extrapolates.
    """
    vs = np.r_[-extend, np.asarray(bp, float), extend]
    ik = np.zeros(len(vs))
    j0 = int(np.searchsorted(vs, 0.0)) - 1
    cur = 0.0
    for j in range(j0, len(vs) - 1):
        lo = 0.0 if j == j0 else vs[j]
        cur += G[j] * (vs[j + 1] - lo)
        ik[j + 1] = cur
    cur = 0.0
    for j in range(j0, -1, -1):
        hi = 0.0 if j == j0 else vs[j + 1]
        cur -= G[j] * (hi - vs[j])
        ik[j] = cur
    return vs, ik


# ---- ideal model -------------------------------------------------------------
@njit(cache=True)
def derivative(v1, v2, i, rt, rl, c1, c2, ind, vk, ik):
    g = np.interp(v1, vk, ik)
    ir = (v2 - v1) / rt
    return (ir - g) / c1, (-ir - i) / c2, (v2 - rl * i) / ind


@njit(cache=True)
def step(v1, v2, i, dt, rt, rl, c1, c2, ind, vk, ik):
    a, b, c = derivative(v1, v2, i, rt, rl, c1, c2, ind, vk, ik)
    d, e, f = derivative(v1 + dt * .5 * a, v2 + dt * .5 * b, i + dt * .5 * c, rt, rl, c1, c2, ind, vk, ik)
    h, j, k = derivative(v1 + dt * .5 * d, v2 + dt * .5 * e, i + dt * .5 * f, rt, rl, c1, c2, ind, vk, ik)
    m, n, p = derivative(v1 + dt * h, v2 + dt * j, i + dt * k, rt, rl, c1, c2, ind, vk, ik)
    return v1 + dt / 6 * (a + 2 * d + 2 * h + m), v2 + dt / 6 * (b + 2 * e + 2 * j + n), i + dt / 6 * (c + 2 * f + 2 * k + p)


@njit(cache=True)
def trajectory(y, rt, rl, c1, c2, ind, vk, ik, dt, nsteps, nskip):
    """States (v1, v2, iL) at steps nskip..nsteps from the state y."""
    out = np.empty((nsteps - nskip + 1, 3))
    v1, v2, i = y
    for k in range(nsteps + 1):
        if k >= nskip:
            out[k - nskip, 0] = v1; out[k - nskip, 1] = v2; out[k - nskip, 2] = i
        if k < nsteps:
            v1, v2, i = step(v1, v2, i, dt, rt, rl, c1, c2, ind, vk, ik)
    return out


@njit(cache=True)
def hold(y, rt, rl, c1, c2, ind, vk, ik, dt, nsettle, ncollect):
    """Integrate nsettle steps, then ncollect more collecting the maxima of v1."""
    v1, v2, i = y
    for k in range(nsettle):
        v1, v2, i = step(v1, v2, i, dt, rt, rl, c1, c2, ind, vk, ik)
    maxima = np.empty(ncollect // 2 + 1)
    count = 0
    p2 = v1; p1 = v1
    for k in range(ncollect):
        v1, v2, i = step(v1, v2, i, dt, rt, rl, c1, c2, ind, vk, ik)
        if p1 > p2 and p1 >= v1:
            maxima[count] = p1; count += 1
        p2 = p1; p1 = v1
    return np.array([v1, v2, i]), maxima[:count]


# ---- bench model -------------------------------------------------------------
def bench_params(c1, c2, l0, nu, r0, rho, RA, AA, RB, AB, vpA, vnA, vpB, vnB, sr, tauA, tauB):
    """The parameter vector P of the bench kernels (SI units), in BENCH_FIELDS order."""
    return np.array([c1, c2, l0, nu, r0, rho, RA, AA, RB, AB, vpA, vnA, vpB, vnB, sr, tauA, tauB], float)


def static_pwl(P):
    """(breakpoints, slopes) of the bench diode with SR -> inf, tau -> 0."""
    c1, c2, l0, nu, r0, rho, RA, AA, RB, AB, vpA, vnA, vpB, vnB, sr, tauA, tauB = P
    Ga = (1 - AA) / RA + (1 - AB) / RB
    Gb = (1 - AA) / RA + 1 / RB
    Gc = 1 / RA + 1 / RB
    return [vnA / AA, vnB / AB, vpB / AB, vpA / AA], [Gc, Gb, Ga, Gb, Gc]


def bench_state(v1, v2, iL, P):
    """Full state with the op-amp outputs where they would sit statically."""
    AA, AB, vpA, vnA, vpB, vnB = P[7], P[9], P[10], P[11], P[12], P[13]
    return np.array([v1, v2, iL, min(max(AA * v1, vnA), vpA), min(max(AB * v1, vnB), vpB)])


@njit(cache=True)
def _clamp(x, lo, hi):
    return lo if x < lo else (hi if x > hi else x)


@njit(cache=True)
def _opamp_rate(target, vo, tau, sr, vn, vp):
    d = (target - vo) / tau
    if d > sr:
        d = sr
    elif d < -sr:
        d = -sr
    if (vo >= vp and d > 0.0) or (vo <= vn and d < 0.0):
        d = 0.0
    return d


@njit(cache=True)
def bench_derivative(v1, v2, i, va, vb, rt, P):
    c1, c2, l0, nu, r0, rho, RA, AA, RB, AB, vpA, vnA, vpB, vnB, sr, tauA, tauB = P
    inr = (v1 - va) / RA + (v1 - vb) / RB
    ir = (v2 - v1) / rt
    ai = abs(i)
    return ((ir - inr) / c1, (-ir - i) / c2, (v2 - (r0 + rho * ai) * i) / (l0 + nu * ai),
            _opamp_rate(AA * v1, va, tauA, sr, vnA, vpA), _opamp_rate(AB * v1, vb, tauB, sr, vnB, vpB))


@njit(cache=True)
def bench_step(y, dt, rt, P):
    vpA, vnA, vpB, vnB = P[10], P[11], P[12], P[13]
    k1 = bench_derivative(y[0], y[1], y[2], y[3], y[4], rt, P)
    k2 = bench_derivative(y[0] + .5 * dt * k1[0], y[1] + .5 * dt * k1[1], y[2] + .5 * dt * k1[2],
                          _clamp(y[3] + .5 * dt * k1[3], vnA, vpA), _clamp(y[4] + .5 * dt * k1[4], vnB, vpB), rt, P)
    k3 = bench_derivative(y[0] + .5 * dt * k2[0], y[1] + .5 * dt * k2[1], y[2] + .5 * dt * k2[2],
                          _clamp(y[3] + .5 * dt * k2[3], vnA, vpA), _clamp(y[4] + .5 * dt * k2[4], vnB, vpB), rt, P)
    k4 = bench_derivative(y[0] + dt * k3[0], y[1] + dt * k3[1], y[2] + dt * k3[2],
                          _clamp(y[3] + dt * k3[3], vnA, vpA), _clamp(y[4] + dt * k3[4], vnB, vpB), rt, P)
    out = np.empty(5)
    for j in range(5):
        out[j] = y[j] + dt / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])
    out[3] = _clamp(out[3], vnA, vpA)
    out[4] = _clamp(out[4], vnB, vpB)
    return out


@njit(cache=True)
def bench_trajectory(y, rt, dt, nsteps, nskip, P):
    """States (v1, v2, iL, voA, voB) at steps nskip..nsteps from the state y."""
    out = np.empty((nsteps - nskip + 1, 5))
    y = y.copy()
    for k in range(nsteps + 1):
        if k >= nskip:
            out[k - nskip] = y
        if k < nsteps:
            y = bench_step(y, dt, rt, P)
    return out


@njit(cache=True)
def bench_hold(y, rt, dt, nsettle, ncollect, P):
    """Integrate nsettle steps, then ncollect more collecting the maxima of v1."""
    y = y.copy()
    for k in range(nsettle):
        y = bench_step(y, dt, rt, P)
    maxima = np.empty(ncollect // 2 + 1)
    count = 0
    p2 = y[0]; p1 = y[0]
    for k in range(ncollect):
        y = bench_step(y, dt, rt, P)
        if p1 > p2 and p1 >= y[0]:
            maxima[count] = p1; count += 1
        p2 = p1; p1 = y[0]
    return y, maxima[:count]
