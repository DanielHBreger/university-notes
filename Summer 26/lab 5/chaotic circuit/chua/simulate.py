"""
Simulate the Chua oscillator of the lab with three models of the same circuit.

    C1 dv1/dt = (v2 - v1)/Rt - i_NR(v1)        Rt = R0 + Rpot
    C2 dv2/dt = (v1 - v2)/Rt - iL
    L  diL/dt = v2 - rL iL

--model ideal   The plan's model: constant nominal components (C1 = 10 nF,
                C2 = 100 nF, L = 18 mH) and the nonlinear element of Table 1:
                the five slopes and four breakpoints fitted to the V-I record
                (M1, diode_fit.json), joined into the continuous law through
                the origin (the fit's current offset, one scope code of
                channel offset, removed). --element intersections instead
                joins the fitted lines where they intersect (inner breakpoints
                -1.26/+1.20 V instead of -1.09/+0.58 V; the element of the
                presentation's animations). --symmetric averages the two
                shoulders and inner breakpoints; --bp-scale moves the inner
                breakpoints.

--model static  The identified circuit (identified.json) with everything held
                constant: the small-signal C1, C2, L0 and r0 and the five-segment
                law the records give, integrated like the ideal model. It isolates
                what the constant-component model can and cannot do with the
                right numbers.

--model bench   The circuit identify.py recovers from the sweep records
                (identified.json): C2 as measured, a capacitor at node 1 whose
                capacitance grows with the swing of v1 since its last turning
                point (Rayleigh law of a ferroelectric ceramic, node 1's law),
                an inductor whose inductance and loss grow with current
                (Rayleigh law of its core), and the element as Kennedy's
                two-op-amp diode with the rails and gains the fitted segments
                imply, integrated with a first-order lag and a slew-rate limit
                per op-amp. Stage B's lag is the one node 1's law measures (a
                lag adds A_B tau_B / R_B of capacitance while the stage is
                linear); stage A's and the slew rate are the TL082 datasheet
                values (3 MHz gain-bandwidth: 0.059 us at gain 1.1; 13 V/us).
                --constant-c1 gives the model without the node-1 law (constant
                C1, datasheet stage-B lag TAU_B = 0.40 us); --datasheet-lag
                keeps the capacitor law with stage B at TAU_B (the report's
                model, which leaves the measured lag out); --tau-b, --slew
                change the op-amps.

Sweep protocol. The bifurcation sweep carries the state from one resistance to
the next, as the knob does: down from the negative outer equilibrium (the
bench's forward sweep) and up from the large outer cycle (the back sweep), over
the whole dial by default. Each resistance is held --settle + --collect seconds
(the first one FIRST_SETTLE longer) and the maxima of v1 in the collect part are
kept; `label` names the regime of each resistance from its maxima and
`transitions` reads off the rows of the report's Table 2. --protocol seeded
(ideal model only) restarts every resistance from the seed instead. Integration
is fixed-step RK4 compiled with numba (integration.py); the bench model needs
dt = 0.1 us for its op-amp lags.

The --r portrait runs start on the negative outer equilibrium (nudged, as the
forward sweep does) and on the large outer cycle, so that both attractors are
shown where they coexist.

Measured records for the portrait panels come from the sweep folders: every
"<sweep>_rpot.csv" beside this script maps a record to its fitted Rpot, and the
record nearest each --r value (within 5 ohm, forward sweep first) is drawn.
Channel convention: CH1 = v1 (across the element), CH2 = v2.

Outputs beside this script (--tag adds a suffix):
    simulated_nr.png            i_NR(v) with the load lines of the chosen R
    simulated_portraits.png     v2 against v1 per --r: the two starts and the record
    simulated_timeseries.png    v1(t) per --r
    simulated_bifurcation.png   maxima of v1 against Rpot, measured points underneath
    simulated_sweep.npz, simulated_transitions.json   the sweep's maxima and Table 2 rows
    simulated/                  with --export: one scope-format CSV per run

Usage:
    python simulate.py --tag _nominal               # the plan's model, nominal parts, Table 1 element
    python simulate.py --model bench --tag _bench   # identified circuit
    python simulate.py --model bench --datasheet-lag --tag _bench_datasheet   # the report's model (b)
    python simulate.py --model bench --r 806 628 320 --export --t 0.2
"""
import argparse
import csv
import glob
import json
import math
import os
import re
from dataclasses import dataclass
from functools import lru_cache
from itertools import dropwhile, takewhile
from typing import Callable, NamedTuple, Optional

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scope_data import R0, read_scope
import integration as kern

HERE = os.path.dirname(os.path.abspath(__file__))
V_EXTEND = 40.0     # V, how far the outer segments are extrapolated


@dataclass(frozen=True)
class Circuit:
    """
    The constant components of the three-state model, SI units: the two
    capacitors, the inductor and its series resistance, and the fixed
    resistor R0 in series with the potentiometer (Rt = R0 + Rpot).
    The defaults are the plan's nominal values.
    """
    c1: float = 10e-9
    c2: float = 100e-9
    l: float = 18e-3
    rl: float = 0.0
    r0: float = R0

    def rt(self, rpot):
        """Total coupling resistance at the potentiometer setting `rpot`."""
        return self.r0 + rpot


NOMINAL = Circuit()

# ---- bench model: the op-amp constants identify.py cannot measure ------------------
GBW = kern.GBW      # Hz, TL082 gain-bandwidth (datasheet, typical)
TAU_A = 1.115 / (2 * np.pi * GBW)   # s, stage A at its identified gain 1.115: 0.059 us; adds tau_A AA/RA = 0.26 nF to node 1
TAU_B = 7.62 / (2 * np.pi * GBW)    # s, stage B at its identified gain 7.62: 0.40 us
SLEW = 13e6         # V/s, both stages (TL082 datasheet)
BENCH_DT = 0.1e-6   # s, RK4 step the op-amp lags need
IDEAL_DT = 0.5e-6   # s, RK4 step of the three-state model
FIRST_SETTLE = 0.2  # s, extra settling at the first resistance of a sweep


# ---- the nonlinear element of the ideal model --------------------------------
def load_vi_fit(path=os.path.join(HERE, '..', 'diode_fit.json')):
    """The M1 segments [(v_lo, v_hi, slope S, intercept A), ...] against the element voltage."""
    with open(path) as fh:
        cfg = json.load(fh)
    if cfg.get('voltage_basis') != 'element':
        raise ValueError(f'{path} must be fitted against the element voltage (run find_breakpoints.py)')
    segs = [(s['v_lo'], s['v_hi'], s['slope_S'], s['intercept_A']) for s in cfg['segments']]
    if len(segs) != 5:
        raise ValueError('five segments expected')
    return segs


@lru_cache(maxsize=None)
def vi_segments():
    """The M1 segments of diode_fit.json, read once."""
    return tuple(load_vi_fit())


def fitted_offset(segments=None):
    """g(0) of the raw fit: the intercept of the segment containing v = 0."""
    for lo, hi, m, b in (segments or vi_segments()):
        if lo <= 0 <= hi:
            return b
    return 0.0


def pwl_from_segments(segments, i0=0.0):
    """(breakpoints, slopes) with the offset i0 removed and adjacent lines intersected."""
    lines = [(m, b - i0) for _, _, m, b in segments]
    bp = [(b2 - b1) / (m1 - m2) for (m1, b1), (m2, b2) in zip(lines[:-1], lines[1:])]
    if np.any(np.diff(bp) <= 0):
        raise ValueError('fitted segments do not intersect in order')
    return bp, [m for m, _ in lines]


def scale_inner_breakpoints(bp, G, scale):
    """Move the two inner breakpoints inward by `scale`, keeping every slope."""
    return [bp[0], bp[1] * scale, bp[2] * scale, bp[3]], list(G)


def symmetrize_inner(bp, G):
    """
    Odd-symmetric shoulders: the two Gb slopes are averaged and the inner
    breakpoints placed at +-(mean of their magnitudes). The saturation
    segments keep their fitted slopes and breakpoints.
    """
    b = 0.5 * (abs(bp[1]) + abs(bp[2]))
    Gb = 0.5 * (G[1] + G[3])
    return [bp[0], -b, b, bp[3]], [G[0], Gb, G[2], Gb, G[4]]


def make_g(i0=None, bp_scale=1.0, symmetric=False, segments=None):
    """
    The element with the fitted lines joined where they intersect, as a
    callable g(v) [A] (see pwl_g). Defaults: the M1 fit with its current
    offset removed, breakpoints as fitted.
    """
    segments = segments or vi_segments()
    if i0 is None:
        i0 = fitted_offset(segments)
    bp, G = pwl_from_segments(segments, i0)
    bp, G = scale_inner_breakpoints(bp, G, bp_scale)
    if symmetric:
        bp, G = symmetrize_inner(bp, G)
    return pwl_g(bp, G)


def table1_g(segments=None):
    """
    The element of Table 1: the fitted slopes and breakpoints (diode_fit.json)
    as the continuous five-segment law through the origin. The M1 lines do not
    meet at the fitted breakpoints (the inner segment spans three current
    codes), so a continuous law cannot keep their intercepts; it keeps the
    nine numbers the report tabulates.
    """
    segments = segments or vi_segments()
    return pwl_g([seg[1] for seg in segments[:4]], [seg[2] for seg in segments])


def pwl_g(bp, G):
    """
    Callable g(v) of the continuous five-segment law (bp: 4 V, G: 5 S), with
    attributes .knots (for the compiled kernels), .segments, .bp and .G.
    """
    vk, ik = kern.pwl_knots(bp, G, V_EXTEND)

    def g(v):
        return np.interp(v, vk, ik)
    g.knots = (vk, ik)
    g.segments = [(vk[j], vk[j + 1], G[j], ik[j] - G[j] * vk[j]) for j in range(5)]
    g.bp, g.G = list(bp), list(G)
    return g


# ---- equilibria and their stability -------------------------------------------------
def equilibria(g, Rt, rL):
    """
    Fixed points v1 of the circuit: iL = (v1 - v2)/Rt and v2 = rL*iL, so
    g(v1) = -v1/(Rt + rL), the load line through the origin.
    """
    from scipy.optimize import brentq

    def h(v):
        return g(v) + v / (Rt + rL)
    vs = np.linspace(-15, 15, 30001)
    hv = h(vs)
    roots = [float(v) for v in vs[hv == 0]]       # an exact zero gives no sign change
    idx = np.flatnonzero(hv[:-1] * hv[1:] < 0)
    roots += [brentq(h, vs[k], vs[k + 1]) for k in idx]
    roots.sort()
    return [r for i, r in enumerate(roots) if i == 0 or r - roots[i - 1] > 1e-9]


def stability(g, v1, Rt, circuit, dv=1e-4):
    """(Jacobian eigenvalues, local slope G of the element) at the fixed point v1."""
    c1, c2, l = circuit.c1, circuit.c2, circuit.l
    G = (g(v1 + dv) - g(v1 - dv)) / (2 * dv)
    J = np.array([[(-1 / Rt - G) / c1, 1 / (Rt * c1), 0.0],
                  [1 / (Rt * c2), -1 / (Rt * c2), -1 / c2],
                  [0.0, 1 / l, -circuit.rl / l]])
    return np.linalg.eigvals(J), G


def hopf_point(g, circuit, side=-1, rt_range=(1200.0, 3200.0)):
    """
    (Rt, v1*, period) where the outer equilibrium on `side` loses stability
    through its complex pair as Rt decreases; None if it never does. This is
    the DC -> limit cycle transition of the measured forward sweep.
    """
    from scipy.optimize import brentq

    def outer(Rt):
        return [v for v in equilibria(g, Rt, circuit.rl) if np.sign(v) == side]

    def growth_rate(Rt):
        eqs = outer(Rt)
        if not eqs:
            return np.nan
        ev, _ = stability(g, eqs[0], Rt, circuit)
        ev = ev[np.abs(ev.imag) > 1.0]
        return ev.real.max() if len(ev) else np.nan

    Rs = np.arange(rt_range[0], rt_range[1], 10.0)
    vals = np.array([growth_rate(R) for R in Rs])
    idx = [k for k in np.flatnonzero(np.sign(vals[:-1]) != np.sign(vals[1:]))
           if np.isfinite(vals[k]) and np.isfinite(vals[k + 1])]
    if not idx:
        return None
    Rh = brentq(growth_rate, Rs[idx[-1]], Rs[idx[-1] + 1])
    v = outer(Rh)[0]
    ev, _ = stability(g, v, Rh, circuit)
    return Rh, v, 2 * np.pi / np.abs(ev.imag).max()


# ---- ideal model: integrator ---------------------------------------------------
def rhs(y, Rt, g, circuit):
    """Time derivative of a batch: y has shape (3, n) with rows v1, v2, iL, and Rt shape (n,)."""
    v1, v2, iL = y
    ir = (v2 - v1) / Rt
    return np.stack([(ir - g(v1)) / circuit.c1, (-ir - iL) / circuit.c2, (v2 - circuit.rl * iL) / circuit.l])


def rk4_step(y, dt, Rt, g, circuit):
    k1 = rhs(y, Rt, g, circuit)
    k2 = rhs(y + 0.5 * dt * k1, Rt, g, circuit)
    k3 = rhs(y + 0.5 * dt * k2, Rt, g, circuit)
    k4 = rhs(y + dt * k3, Rt, g, circuit)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def collect_maxima(y, dt, Rt, g, circuit, n_steps, skip=0):
    """
    Advance a batch by n_steps of the numpy RK4 and collect the maxima of v1
    per column after the first `skip` steps. Returns (state, maxima lists).
    """
    maxima = [[] for _ in range(y.shape[1])]
    before = y[0].copy()
    last = y[0].copy()
    for k in range(1, n_steps + 1):
        y = rk4_step(y, dt, Rt, g, circuit)
        current = y[0]
        if k > skip:
            for i in np.flatnonzero((last > before) & (last >= current)):
                maxima[i].append(last[i])
        before, last = last, current.copy()
    return y, maxima


def integrate(g, circuit, rpot, y0, t_end, dt, keep=True, t_skip=0.0):
    """
    Integrate a batch of runs at the resistances `rpot` from the states y0 of
    shape (3, n). An element with .knots runs in the compiled kernel, any
    other callable in the vectorised numpy RK4.

    keep=True  -> (t, Y) with Y of shape (nt, 3, n) for t >= t_skip.
    keep=False -> the list of v1-maxima per column found after t_skip.
    """
    if (not np.isfinite([t_end, dt, t_skip, circuit.rl]).all() or dt <= 0 or not 0 <= t_skip < t_end
            or circuit.rl < 0):
        raise ValueError('require dt > 0, 0 <= t_skip < t_end and rL >= 0')
    rpot = np.atleast_1d(np.asarray(rpot, float))
    Rt = circuit.rt(rpot)
    n = len(rpot)
    y = np.array(y0, float).reshape(3, n)
    nsteps = int(round(t_end / dt))
    n_skip = int(round(t_skip / dt))
    c1, c2, l, rl = circuit.c1, circuit.c2, circuit.l, circuit.rl

    if hasattr(g, 'knots'):
        vk, ik = g.knots
        if keep:
            out = np.stack([kern.trajectory(y[:, j], Rt[j], rl, c1, c2, l, vk, ik, dt, nsteps, n_skip)
                            for j in range(n)], axis=2)
            return dt * np.arange(n_skip, nsteps + 1), out
        return [kern.hold(y[:, j], Rt[j], rl, c1, c2, l, vk, ik, dt, n_skip, nsteps - n_skip)[1]
                for j in range(n)]

    if not keep:
        return collect_maxima(y, dt, Rt, g, circuit, nsteps, skip=n_skip)[1]
    out = np.empty((nsteps - n_skip + 1, 3, n))
    out[0] = y
    j = 1
    for k in range(1, nsteps + 1):
        y = rk4_step(y, dt, Rt, g, circuit)
        if k == n_skip:
            out[0] = y
            j = 1
        elif k > n_skip:
            out[j] = y
            j += 1
    return dt * np.arange(n_skip, nsteps + 1), out


def midpoint(v1, v2, rpot, r0=R0):
    """Node 3 (CH3): between R0 on the C2 side and Rpot on the C1 side."""
    return v2 + (v1 - v2) * r0 / (r0 + rpot)


def sweep_starts(g, circuit, rpot, side=-1, eps=0.05):
    """
    Initial states (3, n) for the two bifurcation sweeps.

    forward: on the outer equilibrium of the chosen side, nudged by `eps` in
             v1; where that side has no equilibrium the origin is used.
    back:    on the large outer limit cycle, (v1, v2, iL) = (7 V, 6 V, 0),
             which is inside its basin wherever it exists.
    """
    rpot = np.atleast_1d(np.asarray(rpot, float))
    fwd = np.zeros((3, len(rpot)))
    for k, r in enumerate(rpot):
        Rt = circuit.rt(r)
        eqs = [v for v in equilibria(g, Rt, circuit.rl) if np.sign(v) == side]
        v = eqs[0] if eqs else 0.0
        iL = v / (Rt + circuit.rl)
        fwd[:, k] = (v + eps, circuit.rl * iL, iL)
    back = np.zeros((3, len(rpot)))
    back[0], back[1] = 7.0, 6.0
    return fwd, back


def hold_batch(y, Rt, g, circuit, dt, n_settle, n_collect):
    """Integrate a batch (3, n) at fixed Rt: settle, then collect the v1 maxima. Returns (state, maxima)."""
    if hasattr(g, 'knots'):
        vk, ik = g.knots
        res = [kern.hold(y[:, j], Rt[j], circuit.rl, circuit.c1, circuit.c2, circuit.l, vk, ik,
                         dt, n_settle, n_collect)
               for j in range(y.shape[1])]
        return np.column_stack([r[0] for r in res]), [list(r[1]) for r in res]
    for _ in range(n_settle):
        y = rk4_step(y, dt, Rt, g, circuit)
    return collect_maxima(y, dt, Rt, g, circuit, n_collect)


def sweep_order(rsw, direction):
    """Indices of `rsw` in the order the sweep visits them: 'down' from the top, 'up' from the bottom."""
    if direction not in ('down', 'up'):
        raise ValueError('direction must be down or up')
    return np.argsort(rsw)[::-1] if direction == 'down' else np.argsort(rsw)


def sweep_continuation(g, circuit, rsw, dt, t_settle, t_collect, direction):
    """
    Maxima of v1 along a sweep that carries the state from one resistance to
    the next. 'down' starts on the negative outer equilibrium at the highest
    R (the bench's forward sweep); 'up' starts on the large outer cycle at the
    lowest R (the back sweep). Returns the list of maxima per R, in the order
    of `rsw`.
    """
    if dt <= 0 or t_settle < 0 or t_collect <= 0:
        raise ValueError('invalid sweep integration times')
    rsw = np.asarray(rsw, float)
    if rsw.ndim != 1 or not len(rsw) or not np.isfinite(rsw).all():
        raise ValueError('need a nonempty finite resistance vector')
    order = sweep_order(rsw, direction)
    fwd, back = sweep_starts(g, circuit, rsw[order[:1]])
    y = fwd if direction == 'down' else back
    n_settle, n_collect = int(round(t_settle / dt)), int(round(t_collect / dt))
    y, _ = hold_batch(y, circuit.rt(rsw[order[:1]]), g, circuit, dt, int(round(FIRST_SETTLE / dt)), 2)
    out = [None] * len(rsw)
    for k in order:
        y, mx = hold_batch(y, circuit.rt(rsw[k:k + 1]), g, circuit, dt, n_settle, n_collect)
        out[k] = mx[0]
    return out


# ---- bench model ---------------------------------------------------------------
def bench_params(path=os.path.join(HERE, 'identified.json'), tau_a=TAU_A, tau_b=None, slew=SLEW, tau_m=kern.TAU_M,
                 c1_law=True):
    """
    The bench kernel's parameter vector from identified.json. The small-signal
    capacitance identify.py measures on node 1 is C1 + tau_A AA/RA, so C1 is
    that value minus the lag's share. With c1_law (and a law in the file) the
    capacitor grows with the swing as identify.py's node-1 law says, and stage
    B's lag is the one that law measures; otherwise C1 is constant and the lag
    is the datasheet's. An explicit tau_b overrides either.
    """
    with open(path) as fh:
        d = json.load(fh)
    k, ray = d['diode']['kennedy'], d['inductor']['rayleigh']
    law = d['C1_nF'].get('law') if c1_law else None
    c1 = d['C1_nF']['small_signal_loop'] * 1e-9 - tau_a * k['AA'] / k['RA_ohm']
    if tau_b is None:
        tau_b = law['stage_B_lag_us'] * 1e-6 if law else TAU_B
    shape = (law['a0_nF_per_V'] * 1e-9, law['a1_nF_per_V'] * 1e-9, law['knot_V']) if law else (0.0, 0.0, 0.5)
    return kern.bench_params(c1, d['C2_nF'] * 1e-9, ray['L0_H'], ray['nu_H_per_A'], ray['r0_ohm'],
                             ray['rho_ohm_per_A'], k['RA_ohm'], k['AA'], k['RB_ohm'], k['AB'],
                             k['VpA_V'], k['VnA_V'], k['VpB_V'], k['VnB_V'], slew, tau_a, tau_b, tau_m, *shape)


def bench_static_g(P):
    """The bench diode with the op-amps infinitely fast: a five-segment law."""
    bp, G = kern.static_pwl(P)
    return pwl_g(bp, G)


def bench_small_signal(P, r0=R0):
    """
    The bench model linearised: C1 plus the stage-A lag's share (the apparent
    C1), C2, L0 and the winding resistance r0 as a constant-component Circuit.
    With bench_static_g this is the static model.
    """
    return Circuit(c1=P[0] + P[15] * P[7] / P[6], c2=P[1], l=P[2], rl=P[4], r0=r0)


def bench_sweep(P, rsw, dt, t_settle, t_collect, direction, r0=R0):
    """Continuation sweep of the bench model (see sweep_continuation)."""
    rsw = np.asarray(rsw, float)
    order = sweep_order(rsw, direction)
    r_first = rsw[order[0]]
    fwd, back = sweep_starts(bench_static_g(P), bench_small_signal(P, r0), [r_first])
    y = kern.bench_state(*(fwd if direction == 'down' else back)[:, 0], P)
    n_settle, n_collect = int(round(t_settle / dt)), int(round(t_collect / dt))
    y, _ = kern.bench_hold(y, r0 + r_first, dt, int(round(FIRST_SETTLE / dt)), 2, P)
    out = [None] * len(rsw)
    for k in order:
        y, mx = kern.bench_hold(y, r0 + rsw[k], dt, n_settle, n_collect, P)
        out[k] = list(mx)
    return out


def bench_runs(P, rpots, states, t_end, dt, t_skip, r0=R0):
    """(t, Y) with Y (nt, 6, n) for runs from the (3, n) states (v1, v2, iL) at each rpot."""
    nsteps, n_skip = int(round(t_end / dt)), int(round(t_skip / dt))
    Y = np.stack([kern.bench_trajectory(kern.bench_state(*states[:, j], P), r0 + r, dt, nsteps, n_skip, P)
                  for j, r in enumerate(rpots)], axis=2)
    return dt * np.arange(n_skip, nsteps + 1), Y


# ---- measured records ------------------------------------------------------------
def find_records(folder):
    """
    {Rpot: (path, label)} of the measured records available for comparison:
    files named like 'chaos - 716.9 ohm.csv' in `folder`, plus every record
    listed in a '<sweep>_rpot.csv' sidecar beside it (batch_rpot.py output)
    whose divider fit is clean, at the Rpot that fit gave it. The forward
    sweep takes precedence where two sweeps land on the same value.
    """
    out = {}
    for f in glob.glob(os.path.join(folder, '*.csv')):
        m = re.search(r'(\d+(?:\.\d+)?)\s*ohm', os.path.basename(f), re.I)
        if m:
            out[float(m.group(1))] = (f, os.path.basename(f))
    sidecars = sorted(glob.glob(os.path.join(folder, '*_rpot.csv')),
                      key=lambda f: not os.path.basename(f).startswith('forward'))
    for side in sidecars:
        sweep = os.path.basename(side)[:-len('_rpot.csv')]
        sweep_dir = os.path.join(folder, sweep)
        if not os.path.isdir(sweep_dir):
            continue
        with open(side, newline='') as fh:
            for row in csv.DictReader(fh):
                try:
                    r = float(row['rpot_ohm'])
                    residual = float(row['residual_pct'])
                except (KeyError, ValueError, TypeError):
                    continue
                record_path = os.path.join(sweep_dir, row['filename'])
                if (row.get('status', 'ok') == 'ok' and np.isfinite(r) and 0 <= r <= 1000
                        and np.isfinite(residual) and 0 <= residual <= 5 and os.path.isfile(record_path)):
                    out.setdefault(r, (record_path, f'{sweep}/{row["filename"]}'))
    return out


def nearest_record(records, rpot, tol=5.0, prefer='forward/'):
    """
    (Rpot, path, label) of the record closest to `rpot` within `tol` ohm.
    Records of the `prefer` sweep win when they have one within tolerance.
    """
    for pool in ({r: v for r, v in records.items() if v[1].startswith(prefer)}, records):
        if pool:
            r = min(pool, key=lambda x: abs(x - rpot))
            if abs(r - rpot) <= tol:
                return (r,) + tuple(pool[r])
    return None


def load_record(path, n_max=None):
    """(t, v1, v2) of a scope record, the first n_max samples."""
    t, d = read_scope(path)
    return t[:n_max], d[:n_max, 0], d[:n_max, 1]


def load_measured_bifurcation(folder):
    """(Rpot, maximum, sweep) of every measured maximum in the *_bifurcation_points.csv sidecars in `folder`."""
    rows = []
    for path in glob.glob(os.path.join(folder, '*_bifurcation_points.csv')):
        sweep = os.path.basename(path)[:-len('_bifurcation_points.csv')]
        sweep_dir = os.path.join(folder, sweep)
        if not os.path.isdir(sweep_dir):
            continue
        with open(path, newline='') as fh:
            reader = csv.DictReader(fh)
            if not {'sweep', 'filename', 'rpot_ohm', 'max_v'} <= set(reader.fieldnames or ()):
                continue        # not a bifurcation.py sidecar
            exists = {}         # one file-system check per record, not per maximum
            for row in reader:
                try:
                    r, m, name = float(row['rpot_ohm']), float(row['max_v']), row['filename']
                    if name not in exists:
                        exists[name] = os.path.isfile(os.path.join(sweep_dir, name))
                except (TypeError, ValueError):
                    continue        # a short or malformed row
                if (row['sweep'] == sweep and math.isfinite(r) and math.isfinite(m) and 0 <= r <= 1000
                        and exists[name]):
                    rows.append((r, m, sweep))
    return rows


# ---- regime labels and transitions ---------------------------------------------------
PERIOD_FLOOR_V = 2e-3     # maxima that repeat at lag p to this (mean |m_i - m_(i+p)|) have period p
MAX_PERIOD = 64


def maxima_period(m, floor=PERIOD_FLOOR_V):
    """Smallest lag p <= MAX_PERIOD at which the maxima repeat; 0 if none (aperiodic), -1 if too few."""
    m = np.asarray(m, float)
    if len(m) < 2 * MAX_PERIOD + 10:
        return -1
    for p in range(1, MAX_PERIOD + 1):
        if np.mean(np.abs(m[:-p] - m[p:])) < floor:
            return p
    return 0


def label(m, eqs):
    """
    Regime at one resistance from the maxima m of v1 and the equilibria eqs:
    'rest' (no oscillation), 'large cycle' (maxima above 5 V), 'double scroll'
    (maxima beyond half of each outer equilibrium), 'origin orbit' (a small
    orbit when the origin is the only equilibrium), 'P<p>' or 'chaos'.
    """
    m = np.asarray(m, float)
    if len(m) == 0:
        return 'rest'
    vplus = max(eqs) if max(eqs) > 0.1 else None
    vminus = min(eqs) if min(eqs) < -0.1 else None
    if m.max() > 5.0:
        return 'large cycle'
    if len(m) >= 50 and np.ptp(m[-50:]) < 1e-3 and vminus is not None and abs(m[-1] - vminus) < 0.02:
        return 'rest'                   # a decaying transient about the equilibrium
    if vplus is not None and vminus is not None and m.max() > 0.5 * vplus and m.min() < 0.5 * vminus:
        return 'double scroll'
    if vplus is None:
        return 'origin orbit' if np.ptp(m) < 0.5 else 'other'
    p = maxima_period(m)
    if p > 0:
        return f'P{p}'
    if p == 0:
        return 'chaos'
    return 'short'


def regime_runs(labels):
    """[(label, first R, last R), ...] of consecutive equal labels."""
    out = []
    for r, lab in labels:
        if out and out[-1][0] == lab:
            out[-1][2] = r
        else:
            out.append([lab, r, r])
    return [tuple(x) for x in out]


def transitions(rsw, max_down, max_up, g, circuit):
    """
    The rows of the report's Table 2 from the two sweeps: the first resistance
    (going down) of each regime, the ends of the double scroll, and the last
    resistance (going up) of the unbroken large cycle. Also returns the
    labelled sweeps, (R, label) in sweep order.
    """
    rsw = np.asarray(rsw, float)

    def labelled(maxima, order):
        return [(float(rsw[k]), label(maxima[k], equilibria(g, circuit.rt(rsw[k]), circuit.rl))) for k in order]
    down = labelled(max_down, sweep_order(rsw, 'down'))
    up = labelled(max_up, sweep_order(rsw, 'up'))
    moving = [x for x in down if x[1] != 'rest']

    def first(pred, seq, hold=1):
        """First R whose label, and the next hold - 1 labels, satisfy pred (a one-step transient does not count)."""
        for k in range(len(seq) - hold + 1):
            if all(pred(lab) for _, lab in seq[k:k + hold]):
                return seq[k][0]
        return None
    res = {'oscillation begins (sweep)': first(lambda lab: True, moving, 2)}
    for name, done in [('period 1 -> 2', ('P1',)), ('period 2 -> 4', ('P1', 'P2')), ('period 4 -> 8', ('P1', 'P2', 'P4'))]:
        res[name] = first(lambda lab, done=done: lab not in done, moving, 2)
    # chaos begins where the orbit leaves the 1-2-4-8 cascade for good (the 2^n windows close within 1 ohm)
    res['chaos begins'] = first(lambda lab: lab not in ('P1', 'P2', 'P4', 'P8'), moving, 2)
    # the period-3 window: the first unbroken run of P3, P6, P12, P24 below the start of chaos
    window = ('P3', 'P6', 'P12', 'P24')
    chaos = res['chaos begins']
    below = [x for x in down if x[0] < chaos] if chaos is not None else []
    run = list(takewhile(lambda x: x[1] in window, dropwhile(lambda x: x[1] not in window, below)))
    res['period-3 window'] = (run[0][0], run[-1][0]) if run else None
    ds = [r for r, lab in down if lab == 'double scroll']
    res['double scroll from (down)'] = max(ds) if ds else None
    res['double scroll ends (down)'] = min(ds) if ds else None
    res['large cycle from (down)'] = first(lambda lab: lab == 'large cycle', down)
    # going up, the large cycle holds from the bottom of the dial until it first breaks
    held = list(takewhile(lambda x: x[1] == 'large cycle', up))
    res['large cycle survives to (up)'] = held[-1][0] if held else None
    return res, down, up


def stacked_maxima(rsw, maxima, order):
    """(R, M) of one sweep's maxima in sweep order, flattened: one R per maximum."""
    R = np.concatenate([np.full(len(maxima[k]), rsw[k]) for k in order])
    M = np.concatenate([np.asarray(maxima[k], float) for k in order])
    return R, M


def mean_period(t, v1):
    """Mean spacing of the maxima of v1(t), in seconds (nan if fewer than three)."""
    pk = np.flatnonzero((v1[1:-1] > v1[:-2]) & (v1[1:-1] >= v1[2:])) + 1
    return float(np.mean(np.diff(t[pk]))) if len(pk) > 2 else np.nan


# ---- figures ---------------------------------------------------------------------
def plot_nr(g, rpots, path, title, r0=R0):
    """The element's law over the V-I record, with the load line of each resistance."""
    v = np.linspace(-10, 10, 2001)
    fig, ax = plt.subplots(figsize=(8, 5))
    raw = os.path.join(HERE, '..', 'trace1.csv')
    if os.path.exists(raw):
        d = np.genfromtxt(raw, delimiter=',', skip_header=1)
        i_m = (d[:, 1] - d[:, 2]) / 216.0
        i0 = np.nanmedian(i_m[np.abs(d[:, 2]) < 0.15])
        ax.plot(d[:, 2], (i_m - i0) * 1e3, '.', ms=1, alpha=0.15, color='steelblue',
                label='measured (trace1.csv, offset removed)')
    ax.plot(v, g(v) * 1e3, 'crimson', lw=2, label='i_NR(v) of the model')
    for b in g.bp:
        ax.axvline(b, color='grey', ls='--', lw=0.7)
    for r in rpots:
        ax.plot(v, -v / (r0 + r) * 1e3, lw=0.9, label=f'load line, Rpot = {r:g} ohm')
    ax.set_xlabel('v across N_R  (V)')
    ax.set_ylabel('current into N_R  (mA)')
    ax.set_title(title)
    ax.set_ylim(-6, 6)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_portraits(runs, records, path, n_meas):
    """runs: list of (rpot, t, Y_from_equilibrium, Y_from_large_cycle) with Y of shape (nt, >=2)."""
    nrow = len(runs)
    fig, axes = plt.subplots(nrow, 3, figsize=(12, 3.2 * nrow), squeeze=False)
    for row, (r, t, yp, ym) in enumerate(runs):
        for col, (y, tag) in enumerate([(yp, 'from the negative equilibrium'), (ym, 'from the large cycle')]):
            ax = axes[row, col]
            ax.plot(y[:, 0], y[:, 1], lw=0.3, color='C0')
            ax.set_title(f'simulated, Rpot = {r:g} ohm, {tag}', fontsize=9)
        ax = axes[row, 2]
        rec = nearest_record(records, r)
        if rec:
            r_meas, path_meas, rec_label = rec
            _, v1, v2 = load_record(path_meas, n_meas)
            ax.plot(v1, v2, lw=0.3, color='C3')
            ax.set_title(f'measured: {rec_label}, Rpot = {r_meas:.1f} ohm', fontsize=9)
        else:
            ax.set_title('no measured record at this Rpot', fontsize=9)
            ax.set_axis_off()
        for ax in axes[row]:
            ax.grid(alpha=0.3)
            ax.set_xlabel('v1 = v_C1 (V)')
            ax.set_ylabel('v2 = v_C2 (V)')
    fig.suptitle('Chua oscillator: phase portraits (v_C2 against v_C1)')
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_timeseries(runs, path, window=0.01):
    """v1(t) of both runs at each resistance, the first `window` seconds after the transient."""
    nrow = len(runs)
    fig, axes = plt.subplots(nrow, 1, figsize=(12, 1.9 * nrow), squeeze=False, sharex=True)
    for ax, (r, t, yp, ym) in zip(axes[:, 0], runs):
        m = t <= t[0] + window
        ax.plot((t[m] - t[0]) * 1e3, yp[m, 0], lw=0.6, color='C0', label='from the negative equilibrium')
        ax.plot((t[m] - t[0]) * 1e3, ym[m, 0], lw=0.6, color='C1', label='from the large cycle')
        ax.set_ylabel('v1 (V)')
        ax.set_title(f'Rpot = {r:g} ohm', fontsize=9, loc='left')
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=8, loc='upper right')
    axes[-1, 0].set_xlabel('time after transient (ms)')
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_bifurcation(rpot, max_fwd, max_back, measured, path, subtitle='',
                     tags=('sim, forward: R turned down from the negative equilibrium',
                           'sim, back: R turned up from the large cycle')):
    """Simulated maxima of both sweeps over the measured ones."""
    fig, ax = plt.subplots(figsize=(14, 7.5))
    if measured:
        sweeps = sorted({s for _, _, s in measured})
        for s, c in zip(sweeps, ['0.65', '0.8']):
            pts = np.array([(r, m) for r, m, sw in measured if sw == s])
            ax.plot(pts[:, 0], pts[:, 1], '.', ms=0.7, color=c, alpha=0.6, label=f'measured, {s}', zorder=1)
    for maxima, c, tag in [(max_fwd, 'C0', tags[0]), (max_back, 'C3', tags[1])]:
        xs = np.concatenate([np.full(len(m), r) for r, m in zip(rpot, maxima)])
        ys = np.concatenate([np.asarray(m, float) for m in maxima])
        ax.plot(xs, ys, '.', ms=0.45, color=c, alpha=0.8, label=tag, zorder=2)
    ax.set_xlabel('Rpot (ohm)')
    ax.set_ylabel('local maxima of v1 = v_C1 (V)')
    ax.set_title('Bifurcation diagram: simulated maxima of v1 against Rpot' + (f'  ({subtitle})' if subtitle else ''))
    ax.grid(alpha=0.3)
    ax.legend(markerscale=15, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=450)
    plt.close(fig)


# ---- command line ------------------------------------------------------------------
class Model(NamedTuple):
    """The model chosen on the command line."""
    kind: str               # 'ideal', 'static' or 'bench'
    g: Callable             # the static element (bench: op-amps infinitely fast)
    circuit: Circuit        # its constant components (bench: the linearisation)
    P: Optional[np.ndarray]  # bench kernel parameters (bench only)
    dt: float
    subtitle: str


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--model', choices=['ideal', 'static', 'bench'], default='ideal')
    p.add_argument('--r', type=float, nargs='+', default=[806.0, 770.0, 745.0, 700.0, 650.0, 450.0, 320.0],
                   help='potentiometer values for portraits/time series (ohm)')
    p.add_argument('--element', choices=['table1', 'intersections'], default='table1',
                   help="ideal: Table 1's slopes and breakpoints (default) or the fitted lines joined at their intersections")
    p.add_argument('--i0', default='auto',
                   help='ideal: current offset (mA) subtracted from the V-I fit; "auto" removes the fitted g(0)')
    p.add_argument('--bp-scale', type=float, default=1.0, help='ideal: scale of the two inner breakpoints')
    p.add_argument('--symmetric', action='store_true', help='ideal: average the shoulders and inner breakpoints')
    p.add_argument('--c1', type=float, default=NOMINAL.c1 * 1e9,
                   help=f'ideal: C1 in nF (default {NOMINAL.c1 * 1e9:g}, nominal)')
    p.add_argument('--c2', type=float, default=NOMINAL.c2 * 1e9, help='ideal: C2 in nF')
    p.add_argument('--l', type=float, default=NOMINAL.l * 1e3, help='ideal: L in mH')
    p.add_argument('--rl', type=float, default=0.0, help='ideal: inductor series resistance (ohm)')
    p.add_argument('--r0', type=float, default=R0, help='fixed part of the coupling resistance (ohm)')
    p.add_argument('--identified', default=os.path.join(HERE, 'identified.json'), help='bench: identify.py output')
    p.add_argument('--tau-b', type=float, default=None,
                   help='bench: stage-B op-amp lag (us); default the node-1 law\'s, or the datasheet\'s with --constant-c1')
    p.add_argument('--constant-c1', action='store_true',
                   help='bench: node 1 at its constant small-signal capacitance and the datasheet lag (the model '
                        'before the node-1 law)')
    p.add_argument('--datasheet-lag', action='store_true',
                   help='bench: stage B at the TL082 datasheet lag instead of the one node 1\'s law measures, with '
                        'the capacitor law kept (the report\'s model)')
    p.add_argument('--slew', type=float, default=SLEW * 1e-6, help='bench: op-amp slew rate (V/us)')
    p.add_argument('--dt', type=float, default=None, help='RK4 step (s); default 0.5 us ideal, 0.1 us bench')
    p.add_argument('--t', type=float, default=0.1, help='length of each portrait run (s)')
    p.add_argument('--skip', type=float, default=0.06, help='transient discarded before plotting (s)')
    p.add_argument('--sweep', type=float, nargs=3, default=[0, 1000, 1], metavar=('RMIN', 'RMAX', 'STEP'))
    p.add_argument('--protocol', choices=['continuation', 'seeded'], default='continuation',
                   help='continuation: carry the state along (default); seeded: restart every R (ideal only)')
    p.add_argument('--sweep-t', type=float, default=0.1, help='seeded: length of each run (s)')
    p.add_argument('--settle', type=float, default=0.02, help='continuation: seconds discarded after each R step')
    p.add_argument('--collect', type=float, default=0.1, help='continuation: seconds of maxima kept at each R')
    p.add_argument('--no-sweep', action='store_true')
    p.add_argument('--records', default=HERE, help='folder with the sweep folders and their _rpot.csv sidecars')
    p.add_argument('--n-meas', type=int, default=60000, help='samples of each measured record to draw')
    p.add_argument('--export', action='store_true', help='write scope-format CSVs of the --r runs')
    p.add_argument('--tag', default='', help='suffix for the output names')
    p.add_argument('--out', default=HERE, help='where the figures go')
    a = p.parse_args()
    if a.protocol == 'seeded' and a.model != 'ideal':
        p.error('--protocol seeded is only available for the ideal model')
    return a


def ideal_element(a):
    """(g, description) of the ideal model's element as the command line asks for it."""
    if a.element == 'intersections':
        i0 = fitted_offset() if a.i0 == 'auto' else float(a.i0) * 1e-3
        return (make_g(i0, a.bp_scale, a.symmetric),
                f'fitted lines joined at their intersections, offset removed = {i0 * 1e3:.3f} mA')
    g = table1_g()
    bp, G = scale_inner_breakpoints(g.bp, g.G, a.bp_scale)
    if a.symmetric:
        bp, G = symmetrize_inner(bp, G)
    return pwl_g(bp, G), "Table 1's slopes and breakpoints, continuous through the origin"


def build_model(a):
    """The model of the command line, with its components printed."""
    if a.model == 'ideal':
        circuit = Circuit(c1=a.c1 * 1e-9, c2=a.c2 * 1e-9, l=a.l * 1e-3, rl=a.rl, r0=a.r0)
        g, how = ideal_element(a)
        print(f'ideal model: C1 = {circuit.c1 * 1e9:g} nF, C2 = {circuit.c2 * 1e9:g} nF, L = {circuit.l * 1e3:g} mH, '
              f'R0 = {circuit.r0:g} ohm, rL = {circuit.rl:g} ohm')
        print(f'N_R: M1 fit (diode_fit.json), {how}, breakpoint scale = {a.bp_scale:g}'
              + (', shoulders symmetrised' if a.symmetric else ''))
        subtitle = (f'ideal ({a.element}), C1 = {circuit.c1 * 1e9:g} nF, rL = {circuit.rl:g} ohm'
                    + (', symmetric' if a.symmetric else '')
                    + (f', bp scale {a.bp_scale:g}' if a.bp_scale != 1 else '') + f', {a.protocol}')
        return Model('ideal', g, circuit, None, a.dt or IDEAL_DT, subtitle)

    tau_b = TAU_B if a.datasheet_lag else (None if a.tau_b is None else a.tau_b * 1e-6)
    P = bench_params(a.identified, TAU_A, tau_b, a.slew * 1e6, c1_law=not a.constant_c1)
    g = bench_static_g(P)
    circuit = bench_small_signal(P, a.r0)
    if a.model == 'static':
        print(f'static model ({os.path.basename(a.identified)}): C1 = {circuit.c1 * 1e9:.2f} nF, '
              f'C2 = {circuit.c2 * 1e9:.1f} nF, L = {circuit.l * 1e3:.2f} mH, rL = {circuit.rl:.2f} ohm, '
              f'R0 = {circuit.r0:g} ohm; the identified five-segment element')
        return Model('static', g, circuit, None, a.dt or IDEAL_DT,
                     'static model: identified circuit with constant components and a static element')
    print(f'bench model ({os.path.basename(a.identified)}): C1 = {P[0] * 1e9:.2f} nF (+ {(circuit.c1 - P[0]) * 1e9:.2f} nF '
          f'from the stage-A lag), C2 = {P[1] * 1e9:.1f} nF, L = {P[2] * 1e3:.2f} mH + {P[3]:.2f} H/A |d|, '
          f'r = {P[4]:.1f} ohm, core loss {P[5]:.0f} ohm/A |d| d, d = iL - its {P[17] * 1e3:g} ms mean, R0 = {circuit.r0:g} ohm')
    print(f'N_R: Kennedy diode RA = {P[6]:.0f} ohm, AA = {P[7]:.3f}, RB = {P[8]:.0f} ohm, AB = {P[9]:.2f}, '
          f'rails A {P[10]:+.2f}/{P[11]:+.2f} V, B {P[12]:+.2f}/{P[13]:+.2f} V, slew {P[14] * 1e-6:g} V/us, '
          f'lags {P[15] * 1e6:g}/{P[16] * 1e6:g} us')
    law = P[18] or P[19]
    print('node 1: ' + (f'C1 grows with the swing x since the last turning point, + {P[18] * 1e9:.3f} nF/V up to '
                        f'{P[20]:g} V, + {P[19] * 1e9:.4f} nF/V beyond' if law else 'constant C1'))
    return Model('bench', g, circuit, P, a.dt or BENCH_DT,
                 f'bench model, {"node-1 law" if law else "constant C1"}, tau_B = {P[16] * 1e6:.3g} us, '
                 f'slew {P[14] * 1e-6:g} V/us')


def report_equilibria(model, rpots):
    """Print the element's segments and the equilibria at each resistance; return the Hopf point."""
    print(f'{"v_lo":>8} {"v_hi":>8} {"G (mS)":>9} {"c (mA)":>9}')
    for lo, hi, G, c in model.g.segments:
        print(f'{lo:8.3f} {hi:8.3f} {G * 1e3:9.4f} {c * 1e3:9.4f}')
    print('\nequilibria, g(v1) = -v1/(Rt + rL):')
    for r in rpots:
        Rt = model.circuit.rt(r)
        parts = []
        for v in equilibria(model.g, Rt, model.circuit.rl):
            ev, G = stability(model.g, v, Rt, model.circuit)
            parts.append(f'v1 = {v:6.2f} V, G = {G * 1e3:6.3f} mS, {"UNSTABLE" if ev.real.max() > 0 else "stable"}')
        print(f'   Rpot = {r:6.1f}:  ' + '  |  '.join(parts))
    hopf = hopf_point(model.g, model.circuit)
    if hopf:
        Rh, vh, Th = hopf
        print(f'\nHopf point of the negative outer equilibrium: Rpot = {Rh - model.circuit.r0:.1f} ohm, '
              f'v1* = {vh:.2f} V, period {Th * 1e6:.1f} us ({1 / Th:.0f} Hz)')
    return hopf


def portrait_runs(model, rpots, t_end, t_skip):
    """[(rpot, t, Y from the negative equilibrium, Y from the large cycle), ...] with a summary printed."""
    rp = np.repeat(rpots, 2)
    y_fwd, y_back = sweep_starts(model.g, model.circuit, rpots)
    y0 = np.empty((3, 2 * len(rpots)))
    y0[:, 0::2], y0[:, 1::2] = y_fwd, y_back
    print(f'\nintegrating {len(rp)} runs of {t_end} s at dt = {model.dt} ...')
    if model.kind == 'bench':
        t, Y = bench_runs(model.P, rp, y0, t_end, model.dt, t_skip, model.circuit.r0)
    else:
        t, Y = integrate(model.g, model.circuit, rp, y0, t_end, model.dt, keep=True, t_skip=t_skip)
    runs = [(r, t, Y[:, :, 2 * k], Y[:, :, 2 * k + 1]) for k, r in enumerate(rpots)]

    print(f'{"Rpot":>7}  {"start":>6}  {"v1 min":>7}  {"v1 max":>7}  {"v2 min":>7}  {"v2 max":>7}  '
          f'{"n maxima":>8}  {"distinct":>8}  {"period":>8}')
    for r, _, yp, ym in runs:
        for tag, y in [('eq', yp), ('cycle', ym)]:
            v1 = y[:, 0]
            pk = v1[1:-1][(v1[1:-1] > v1[:-2]) & (v1[1:-1] >= v1[2:])]
            T = mean_period(t, v1)
            print(f'{r:7.1f}  {tag:>6}  {v1.min():7.3f}  {v1.max():7.3f}  {y[:, 1].min():7.3f}  {y[:, 1].max():7.3f}  '
                  f'{len(pk):8d}  {len(np.unique(np.round(pk, 2))):8d}  {T * 1e6:6.1f} us')
    return runs


def export_runs(runs, folder, r0):
    """One scope-format CSV per run, with CH3 at the divider's midpoint."""
    os.makedirs(folder, exist_ok=True)
    for r, t, yp, ym in runs:
        for tag, y in [('plus', yp), ('minus', ym)]:
            v3 = midpoint(y[:, 0], y[:, 1], r, r0)
            f = os.path.join(folder, f'sim {r:.1f} ohm {tag}.csv')
            np.savetxt(f, np.column_stack([t - t[0], y[:, 0], y[:, 1], v3]), delimiter=',', fmt='%.6e',
                       header='Time(s),CH1(V),CH2(V),CH3(V)', comments='')
    print(f'\nwrote {2 * len(runs)} records to {folder}')


def bifurcation_sweep(model, a, hopf):
    """Both sweeps, their figure, and the Table 2 transitions (printed and written)."""
    rmin, rmax, step = a.sweep
    rsw = np.arange(rmin, rmax + step / 2, step)
    if a.protocol == 'continuation':
        print(f'\nsweeping {len(rsw)} values of Rpot down from the negative equilibrium and up from the '
              f'large cycle, carrying the state along ({a.settle * 1e3:.0f} ms settle + '
              f'{a.collect * 1e3:.0f} ms collect per step) ...')
        if model.kind == 'bench':
            mf = bench_sweep(model.P, rsw, model.dt, a.settle, a.collect, 'down', model.circuit.r0)
            mb = bench_sweep(model.P, rsw, model.dt, a.settle, a.collect, 'up', model.circuit.r0)
        else:
            mf = sweep_continuation(model.g, model.circuit, rsw, model.dt, a.settle, a.collect, 'down')
            mb = sweep_continuation(model.g, model.circuit, rsw, model.dt, a.settle, a.collect, 'up')
        tags = ('sim, forward: R turned down from the negative equilibrium',
                'sim, back: R turned up from the large cycle')
    else:
        print(f'\nsweeping {len(rsw)} values of Rpot x 2 starts, {a.sweep_t} s each ...')
        y_fwd, y_back = sweep_starts(model.g, model.circuit, rsw)
        mf = integrate(model.g, model.circuit, rsw, y_fwd, a.sweep_t, model.dt, keep=False, t_skip=a.skip)
        mb = integrate(model.g, model.circuit, rsw, y_back, a.sweep_t, model.dt, keep=False, t_skip=a.skip)
        tags = ('sim, forward (each R restarted on the negative equilibrium)',
                'sim, back (each R restarted on the large cycle)')
    measured = load_measured_bifurcation(a.records)
    plot_bifurcation(rsw, mf, mb, measured, os.path.join(a.out, f'simulated_bifurcation{a.tag}.png'),
                     subtitle=model.subtitle, tags=tags)

    res, down, up = transitions(rsw, mf, mb, model.g, model.circuit)
    res['Hopf point (linear)'] = None if not hopf else float(hopf[0] - model.circuit.r0)
    res['Hopf frequency (Hz)'] = None if not hopf else float(1 / hopf[2])
    print('\nregimes, down (forward): ' + ', '.join(f'{lab} {x:g}-{y:g}' for lab, x, y in regime_runs(down)))
    print('regimes, up (back):      ' + ', '.join(f'{lab} {x:g}-{y:g}' for lab, x, y in regime_runs(up)))
    print(f'\ntransitions (ohm; sweep step {step:g} ohm):')
    for k, v in res.items():
        print(f'   {k:32s} {v}')
    R_down, M_down = stacked_maxima(rsw, mf, sweep_order(rsw, 'down'))
    R_up, M_up = stacked_maxima(rsw, mb, sweep_order(rsw, 'up'))
    np.savez_compressed(os.path.join(a.out, f'simulated_sweep{a.tag}.npz'),
                        R_down=R_down, M_down=M_down, R_up=R_up, M_up=M_up)
    with open(os.path.join(a.out, f'simulated_transitions{a.tag}.json'), 'w') as fh:
        json.dump(dict(model=a.model, subtitle=model.subtitle, sweep=[rmin, rmax, step], settle_s=a.settle,
                       collect_s=a.collect, transitions=res,
                       regimes_down=regime_runs(down), regimes_up=regime_runs(up)), fh, indent=1)


def main():
    a = parse_args()
    model = build_model(a)
    hopf = report_equilibria(model, a.r)
    runs = portrait_runs(model, a.r, a.t, a.skip)

    records = find_records(a.records)
    for r in a.r:
        rec = nearest_record(records, r)
        print(f'measured record for Rpot = {r:g} ohm: ' + (f'{rec[2]} at {rec[0]:.1f} ohm' if rec else 'none within 5 ohm'))
    element = 'bench model (op-amps infinitely fast)' if a.model == 'bench' else f'ideal model (M1 fit, {a.element})'
    plot_nr(model.g, a.r, os.path.join(a.out, f'simulated_nr{a.tag}.png'), f'Nonlinear element of the {element}',
            model.circuit.r0)
    plot_portraits(runs, records, os.path.join(a.out, f'simulated_portraits{a.tag}.png'), a.n_meas)
    plot_timeseries(runs, os.path.join(a.out, f'simulated_timeseries{a.tag}.png'))
    if a.export:
        export_runs(runs, os.path.join(a.out, f'simulated{a.tag}'), model.circuit.r0)
    if not a.no_sweep:
        bifurcation_sweep(model, a, hopf)
    print('\nfigures written to', a.out)


if __name__ == '__main__':
    main()
