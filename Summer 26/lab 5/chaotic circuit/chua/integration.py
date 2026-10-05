"""
Compiled RK4 kernels for simulate.py. Every circuit value is an argument, never
a frozen global, so the same code serves any parameter set.

Ideal model (three states v1, v2, iL; the nonlinear element is a static
piecewise-linear law given as interpolation knots):

    C1 v1' = (v2 - v1)/Rt - g(v1)
    C2 v2' = (v1 - v2)/Rt - iL
    L  iL' = v2 - rL iL

Bench model (six states v1, v2, iL, voA, voB, im, and three that track the
turning points of v1: v_rev, v_ext, rising), identify.py's circuit:

    i_A = (v1 - voA)/RA,   i_B = (v1 - voB)/RB,   i_NR = i_A + i_B
    voX' = clip((AX v1 - voX)/tauX, -SR, +SR), voX held inside [VnX, VpX]
    C1(x) v1' = (v2 - v1)/Rt - i_NR,               x = |v1 - v_rev|
    C2 v2' = (v1 - v2)/Rt - iL
    (L0 + nu |d|) iL' = v2 - r0 iL - rho |d| d,     d = iL - im
    im' = (iL - im)/tauM

The two op-amp stages of Kennedy's diode saturate at their output rails
(VpA, VnA) and (VpB, VnB); with SR -> inf and tau -> 0 the element is the
five-segment law
    Ga = (1-AA)/RA + (1-AB)/RB,  Gb = (1-AA)/RA + 1/RB,  Gc = 1/RA + 1/RB,
    inner breakpoints VnB/AB, VpB/AB,  outer breakpoints VnA/AA, VpA/AA.
The inductor follows the Rayleigh law of a ferromagnetic core at low field:
inductance and core loss both grow linearly with the current's excursion d
from its running mean im (time constant tauM, several periods), which is how
identify.py measures them: L_eff and r_eff of each periodic record are taken
about the record's mean current. A steady current therefore sees only L0 and
the winding resistance r0.
The capacitor at node 1 follows the Rayleigh law of a ferroelectric dielectric:
its incremental capacitance grows with the swing x of v1 since v1 last turned,
    C1(x) = c1 + a0 min(x, xk) + a1 max(x - xk, 0),
steeper over the first xk of swing than beyond (identify.py measures all four
numbers). A turning point counts once v1 has come back TURN_V from its extreme,
and the memory of it fades towards v1 over TAU_R, so that the slow drift of
the equilibrium during a sweep is not taken for swing. The tracker is updated
after every step. With a0 = a1 = 0, C1 is the constant c1 and the other six
states evolve exactly as without it.

The explicit RK4 step must satisfy dt < 2.78 tau for the op-amp equations,
so the bench model is integrated with dt = 0.1 us (tauA = 0.059 us).
"""
import numpy as np
from numba import njit

# index of every bench-model parameter in the parameter vector P
BENCH_FIELDS = ('c1', 'c2', 'l0', 'nu', 'r0', 'rho', 'RA', 'AA', 'RB', 'AB',
                'vpA', 'vnA', 'vpB', 'vnB', 'sr', 'tauA', 'tauB', 'tauM', 'c1_a0', 'c1_a1', 'c1_xk')
GBW = 3e6           # Hz, TL082 gain-bandwidth (datasheet, typical): a stage of gain A lags by A/(2 pi GBW)
TAU_M = 3e-3        # s, running mean of iL: about nine periods of the ~330 us oscillation
TURN_V = 0.02       # V, how far v1 must come back from an extreme for it to count as a turning point
TAU_R = 10e-3       # s, fading of the turning-point memory: many periods, well inside a sweep step
N_BENCH = 9         # states of the bench model


def pwl_knots(bp, G, extend=40.0):
    """
    Knots (v, i) of the continuous five-segment law through the origin with
    breakpoints bp (4 values, V) and slopes G (5 values, S), extended to
    +-extend V so that np.interp never extrapolates.
    """
    vs = np.r_[-extend, np.asarray(bp, float), extend]
    ik = np.zeros(len(vs))
    j0 = int(np.searchsorted(vs, 0.0)) - 1       # the segment containing v = 0
    current = 0.0
    for j in range(j0, len(vs) - 1):             # integrate the slopes outward from the origin
        lo = 0.0 if j == j0 else vs[j]
        current += G[j] * (vs[j + 1] - lo)
        ik[j + 1] = current
    current = 0.0
    for j in range(j0, -1, -1):
        hi = 0.0 if j == j0 else vs[j + 1]
        current -= G[j] * (hi - vs[j])
        ik[j] = current
    return vs, ik


# ---- ideal model -------------------------------------------------------------
@njit(cache=True)
def pwl_eval(v, vk, ik):
    """np.interp(v, vk, ik) for the few knots of a piecewise-linear law, without its overhead."""
    n = len(vk)
    if v <= vk[0]:
        return ik[0]
    if v >= vk[n - 1]:
        return ik[n - 1]
    j = 0
    while v > vk[j + 1]:
        j += 1
    return ik[j] + (ik[j + 1] - ik[j]) * (v - vk[j]) / (vk[j + 1] - vk[j])


@njit(cache=True)
def derivative(v1, v2, i, rt, rl, c1, c2, ind, vk, ik):
    """(v1', v2', iL') of the ideal model."""
    g = pwl_eval(v1, vk, ik)
    ir = (v2 - v1) / rt
    return (ir - g) / c1, (-ir - i) / c2, (v2 - rl * i) / ind


@njit(cache=True)
def step(v1, v2, i, dt, rt, rl, c1, c2, ind, vk, ik):
    """One RK4 step of the ideal model."""
    k1 = derivative(v1, v2, i, rt, rl, c1, c2, ind, vk, ik)
    k2 = derivative(v1 + dt * .5 * k1[0], v2 + dt * .5 * k1[1], i + dt * .5 * k1[2], rt, rl, c1, c2, ind, vk, ik)
    k3 = derivative(v1 + dt * .5 * k2[0], v2 + dt * .5 * k2[1], i + dt * .5 * k2[2], rt, rl, c1, c2, ind, vk, ik)
    k4 = derivative(v1 + dt * k3[0], v2 + dt * k3[1], i + dt * k3[2], rt, rl, c1, c2, ind, vk, ik)
    return (v1 + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
            v2 + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]),
            i + dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]))


@njit(cache=True)
def trajectory(y, rt, rl, c1, c2, ind, vk, ik, dt, nsteps, nskip):
    """States (v1, v2, iL) at steps nskip..nsteps from the state y."""
    out = np.empty((nsteps - nskip + 1, 3))
    v1, v2, i = y
    for k in range(nsteps + 1):
        if k >= nskip:
            out[k - nskip, 0] = v1
            out[k - nskip, 1] = v2
            out[k - nskip, 2] = i
        if k < nsteps:
            v1, v2, i = step(v1, v2, i, dt, rt, rl, c1, c2, ind, vk, ik)
    return out


@njit(cache=True)
def hold(y, rt, rl, c1, c2, ind, vk, ik, dt, nsettle, ncollect):
    """Integrate nsettle steps, then ncollect more collecting the maxima of v1. Returns (state, maxima)."""
    v1, v2, i = y
    for k in range(nsettle):
        v1, v2, i = step(v1, v2, i, dt, rt, rl, c1, c2, ind, vk, ik)
    maxima = np.empty(ncollect // 2 + 1)
    count = 0
    before = v1          # v1 two steps back
    last = v1            # v1 one step back
    for k in range(ncollect):
        v1, v2, i = step(v1, v2, i, dt, rt, rl, c1, c2, ind, vk, ik)
        if last > before and last >= v1:
            maxima[count] = last
            count += 1
        before = last
        last = v1
    return np.array([v1, v2, i]), maxima[:count]


# ---- bench model -------------------------------------------------------------
def bench_params(c1, c2, l0, nu, r0, rho, RA, AA, RB, AB, vpA, vnA, vpB, vnB, sr, tauA, tauB, tauM=TAU_M,
                 c1_a0=0.0, c1_a1=0.0, c1_xk=0.5):
    """The parameter vector P of the bench kernels (SI units), in BENCH_FIELDS order."""
    return np.array([c1, c2, l0, nu, r0, rho, RA, AA, RB, AB, vpA, vnA, vpB, vnB, sr, tauA, tauB, tauM,
                     c1_a0, c1_a1, c1_xk], float)


def c1_of_swing(x, P):
    """Node 1's capacitance at a swing x of v1 since its last turning point."""
    return P[0] + P[18] * np.minimum(x, P[20]) + P[19] * np.maximum(x - P[20], 0.0)


@njit(cache=True)
def swing(v, dt, thr=TURN_V, tau=TAU_R):
    """
    |v - v_rev| along a sampled signal, with v_rev tracked exactly as the bench
    kernel tracks it (turning points with hysteresis thr, memory fading over tau).
    """
    out = np.empty_like(v)
    rev, ext, rising = v[0], v[0], 1.0
    for i in range(len(v)):
        rev, ext, rising = _turn(v[i], rev, ext, rising, dt, thr, tau)
        out[i] = abs(v[i] - rev)
    return out


@njit(cache=True)
def _turn(x, rev, ext, rising, dt, thr, tau):
    """The turning-point tracker after v1 has moved to x: (v_rev, v_ext, rising)."""
    if rising > 0:
        if x > ext:
            ext = x
        elif ext - x > thr:
            rev, ext, rising = ext, x, -1.0
    else:
        if x < ext:
            ext = x
        elif x - ext > thr:
            rev, ext, rising = ext, x, 1.0
    rev += dt / tau * (x - rev)
    return rev, ext, rising


def static_pwl(P):
    """(breakpoints, slopes) of the bench diode with SR -> inf, tau -> 0."""
    RA, AA, RB, AB, vpA, vnA, vpB, vnB = P[6:14]
    Ga = (1 - AA) / RA + (1 - AB) / RB
    Gb = (1 - AA) / RA + 1 / RB
    Gc = 1 / RA + 1 / RB
    return [vnA / AA, vnB / AB, vpB / AB, vpA / AA], [Gc, Gb, Ga, Gb, Gc]


def bench_state(v1, v2, iL, P):
    """Full state: the op-amp outputs where they would sit statically, the mean current at iL, no swing yet."""
    AA, AB, vpA, vnA, vpB, vnB = P[7], P[9], P[10], P[11], P[12], P[13]
    return np.array([v1, v2, iL, min(max(AA * v1, vnA), vpA), min(max(AB * v1, vnB), vpB), iL, v1, v1, 1.0])


@njit(cache=True)
def _clamp(x, lo, hi):
    return lo if x < lo else (hi if x > hi else x)


@njit(cache=True)
def _opamp_rate(target, vo, tau, sr, vn, vp):
    """Output rate of an op-amp stage: first-order lag, slew limit, and no motion past a rail."""
    d = (target - vo) / tau
    if d > sr:
        d = sr
    elif d < -sr:
        d = -sr
    if (vo >= vp and d > 0.0) or (vo <= vn and d < 0.0):
        d = 0.0
    return d


@njit(cache=True)
def bench_derivative(v1, v2, i, va, vb, im, vrev, rt, P):
    """Time derivative of the six continuous bench-model states (v_rev held)."""
    c1, c2, l0, nu, r0, rho, RA, AA, RB, AB, vpA, vnA, vpB, vnB, sr, tauA, tauB, tauM, a0, a1, xk = P
    inr = (v1 - va) / RA + (v1 - vb) / RB
    ir = (v2 - v1) / rt
    d = i - im
    ad = abs(d)
    x = abs(v1 - vrev)
    cap = c1 + a0 * min(x, xk) + a1 * max(x - xk, 0.0)
    return ((ir - inr) / cap, (-ir - i) / c2, (v2 - r0 * i - rho * ad * d) / (l0 + nu * ad),
            _opamp_rate(AA * v1, va, tauA, sr, vnA, vpA), _opamp_rate(AB * v1, vb, tauB, sr, vnB, vpB),
            d / tauM)


@njit(cache=True)
def _stage(y, k, h, vnA, vpA, vnB, vpB):
    """The state y + h k with the op-amp outputs held inside their rails."""
    return (y[0] + h * k[0], y[1] + h * k[1], y[2] + h * k[2],
            _clamp(y[3] + h * k[3], vnA, vpA), _clamp(y[4] + h * k[4], vnB, vpB), y[5] + h * k[5])


@njit(cache=True)
def bench_step(y, dt, rt, P):
    """One RK4 step of the bench model, with the op-amp outputs clamped at every stage, then the tracker."""
    vpA, vnA, vpB, vnB = P[10], P[11], P[12], P[13]
    k1 = bench_derivative(y[0], y[1], y[2], y[3], y[4], y[5], y[6], rt, P)
    a = _stage(y, k1, .5 * dt, vnA, vpA, vnB, vpB)
    k2 = bench_derivative(a[0], a[1], a[2], a[3], a[4], a[5], y[6], rt, P)
    b = _stage(y, k2, .5 * dt, vnA, vpA, vnB, vpB)
    k3 = bench_derivative(b[0], b[1], b[2], b[3], b[4], b[5], y[6], rt, P)
    c = _stage(y, k3, dt, vnA, vpA, vnB, vpB)
    k4 = bench_derivative(c[0], c[1], c[2], c[3], c[4], c[5], y[6], rt, P)
    out = np.empty(N_BENCH)
    for j in range(6):
        out[j] = y[j] + dt / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])
    out[3] = _clamp(out[3], vnA, vpA)
    out[4] = _clamp(out[4], vnB, vpB)
    out[6], out[7], out[8] = _turn(out[0], y[6], y[7], y[8], dt, TURN_V, TAU_R)
    return out


@njit(cache=True)
def bench_trajectory(y, rt, dt, nsteps, nskip, P):
    """States (v1, v2, iL, voA, voB, im, v_rev, v_ext, rising) at steps nskip..nsteps from the state y."""
    out = np.empty((nsteps - nskip + 1, N_BENCH))
    y = y.copy()
    for k in range(nsteps + 1):
        if k >= nskip:
            out[k - nskip] = y
        if k < nsteps:
            y = bench_step(y, dt, rt, P)
    return out


@njit(cache=True)
def bench_hold(y, rt, dt, nsettle, ncollect, P):
    """Integrate nsettle steps, then ncollect more collecting the maxima of v1. Returns (state, maxima)."""
    y = y.copy()
    for k in range(nsettle):
        y = bench_step(y, dt, rt, P)
    maxima = np.empty(ncollect // 2 + 1)
    count = 0
    before = y[0]
    last = y[0]
    for k in range(ncollect):
        y = bench_step(y, dt, rt, P)
        if last > before and last >= y[0]:
            maxima[count] = last
            count += 1
        before = last
        last = y[0]
    return y, maxima[:count]
