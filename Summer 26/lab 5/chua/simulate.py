"""
Simulate the Chua oscillator of the lab with two models of the same circuit.

    C1 dv1/dt = (v2 - v1)/Rt - i_NR(v1)        Rt = R0 + Rpot
    C2 dv2/dt = (v1 - v2)/Rt - iL
    L  diL/dt = v2 - rL iL

--model ideal   The plan's model: constant components and the nonlinear
                element as the five straight segments fitted to the V-I record
                (M1, diode_fit.json, element-voltage basis). The fitted
                current offset (the inner segment's intercept, one scope code
                of channel offset) is removed so that the element passes no
                current at v = 0; --i0 keeps or sets it. Adjacent fitted lines
                are joined where they intersect. Defaults: C2 = 100 nF and
                L = 18 mH nominal, C1 = 11.5 nF, the value that puts the
                measured Hopf point in place (--c1 10 for the nominal value,
                which puts it at about 1000 ohm, beyond the dial).
                --symmetric averages the two shoulders and inner breakpoints;
                --bp-scale moves the inner breakpoints.

--model static  The identified circuit (identified.json) with everything held
                constant: the small-signal C1, C2, L0 and r0 and the five-segment
                law the records give, integrated like the ideal model. It isolates
                what the constant-component model can and cannot do with the
                right numbers.

--model bench   The circuit identify.py recovers from the sweep records
                (identified.json): C1 and C2 as measured, an inductor whose
                inductance and loss grow with current (Rayleigh law of its
                core), and the element as Kennedy's two-op-amp diode with the
                rails and gains the fitted segments imply, integrated with a
                first-order lag and a slew-rate limit per op-amp. Two constants
                are not measured and are set here: the stage-B lag TAU_B and
                the slew rate SLEW (see RESULTS.md); --tau-b, --slew change them.

Sweep protocol. The bifurcation sweep carries the state from one resistance to
the next, as the knob does: down from the negative outer equilibrium (the
bench's forward sweep) and up from the large outer cycle (the back sweep). Each
resistance is held --settle + --collect seconds and the maxima of v1 in the
collect part are kept. --protocol seeded (ideal model only) restarts every
resistance from the seed instead. Integration is fixed-step RK4 compiled with
numba (integration.py); the bench model needs dt = 0.1 us for its op-amp lags.

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
    simulated/                  with --export: one scope-format CSV per run

Usage:
    python simulate.py                              # ideal model, C1 = 11.5 nF
    python simulate.py --c1 10 --tag _nominal       # ideal model, nominal C1
    python simulate.py --model bench --tag _bench   # identified circuit
    python simulate.py --model bench --r 806 628 320 --export --t 0.2
"""
import argparse
import csv
import glob
import json
import os
import re
import warnings

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scope_data import R0
import integration as kern

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- ideal model: circuit values ----------------------------------------------
# C2 and L nominal; C1 is the calibrated value (see the docstring). main()
# overrides them from the command line, so read them through the module globals.
C1_NOMINAL = 10e-9
C1 = 11.5e-9
C2 = 100e-9
L = 18e-3
V_EXTEND = 40.0     # V, how far the outer segments are extrapolated

# ---- bench model: the two op-amp constants that identify.py does not measure ----
TAU_A = 0.05e-6     # s, stage A (gain 1.1); adds tau_A AA/RA = 0.2 nF to node 1
TAU_B = 2.4e-6      # s, stage B (gain 7.6); an op-amp of ~0.5 MHz gain-bandwidth
SLEW = 0.5e6        # V/s, both stages (a 741-class op-amp)
BENCH_DT = 0.1e-6   # s, RK4 step the op-amp lags need


def set_components(c1=None, c2=None, l=None, r0=None):
    """Override the ideal model's component values used by every function here."""
    global C1, C2, L, R0
    if c1 is not None: C1 = c1
    if c2 is not None: C2 = c2
    if l is not None: L = l
    if r0 is not None: R0 = r0


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


try:
    VI_SEGMENTS = load_vi_fit()
except (OSError, ValueError, KeyError) as exc:      # pragma: no cover - depends on the working tree
    VI_SEGMENTS = None
    _VI_ERROR = exc


def fitted_offset(segments=None):
    """g(0) of the raw fit: the intercept of the segment containing v = 0."""
    for lo, hi, m, b in (segments or VI_SEGMENTS):
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


def make_g(i0=None, bp_scale=1.0, symmetric=False, segments=None, blend=None, shunt=None):
    """
    The ideal model's element as a callable g(v) [A] with .knots and .segments.
    Defaults: the M1 fit with its current offset removed, breakpoints as fitted.
    (`blend` and `shunt` are accepted for compatibility and ignored: the lines
    are joined at their intersections and the fit is already against the
    element voltage.)
    """
    segments = segments or VI_SEGMENTS
    if segments is None:
        raise RuntimeError(f'no V-I fit available: {_VI_ERROR}')
    if i0 is None:
        i0 = fitted_offset(segments)
    bp, G = pwl_from_segments(segments, i0)
    bp, G = scale_inner_breakpoints(bp, G, bp_scale)
    if symmetric:
        bp, G = symmetrize_inner(bp, G)
    return pwl_g(bp, G)


def pwl_g(bp, G):
    """Callable g(v) of the continuous five-segment law (bp: 4 V, G: 5 S)."""
    vk, ik = kern.pwl_knots(bp, G, V_EXTEND)

    def g(v):
        return np.interp(v, vk, ik)
    g.knots = (vk, ik)
    g.segments = [(vk[j], vk[j + 1], G[j], ik[j] - G[j] * vk[j]) for j in range(5)]
    g.bp, g.G = list(bp), list(G)
    return g


# ---- equilibria ----------------------------------------------------------------
def equilibria(g, Rt, rL):
    """
    Fixed points v1 of the circuit: iL = (v1 - v2)/Rt and v2 = rL*iL, so
    g(v1) = -v1/(Rt + rL), the load line through the origin.
    """
    from scipy.optimize import brentq
    h = lambda v: g(v) + v / (Rt + rL)
    vs = np.linspace(-15, 15, 30001)
    hv = h(vs)
    roots = [float(v) for v in vs[hv == 0]]       # an exact zero gives no sign change
    idx = np.flatnonzero(hv[:-1] * hv[1:] < 0)
    roots += [brentq(h, vs[k], vs[k + 1]) for k in idx]
    roots.sort()
    return [r for i, r in enumerate(roots) if i == 0 or r - roots[i - 1] > 1e-9]


def stability(g, v1, Rt, rL, dv=1e-4, c1=None, c2=None, l=None):
    """(Jacobian eigenvalues, local slope G of the element) at the fixed point v1."""
    c1, c2, l = c1 or C1, c2 or C2, l or L
    G = (g(v1 + dv) - g(v1 - dv)) / (2 * dv)
    J = np.array([[(-1 / Rt - G) / c1, 1 / (Rt * c1), 0.0],
                  [1 / (Rt * c2), -1 / (Rt * c2), -1 / c2],
                  [0.0, 1 / l, -rL / l]])
    return np.linalg.eigvals(J), G


def hopf_point(g, rL, side=-1, rt_range=(1200.0, 3200.0), **comps):
    """
    (Rt, v1*, period) where the outer equilibrium on `side` loses stability
    through its complex pair as Rt decreases; None if it never does. This is
    the DC -> limit cycle transition of the measured forward sweep.
    """
    from scipy.optimize import brentq

    def re_part(Rt):
        eqs = [v for v in equilibria(g, Rt, rL) if np.sign(v) == side]
        if not eqs:
            return np.nan
        ev, _ = stability(g, eqs[0], Rt, rL, **comps)
        ev = ev[np.abs(ev.imag) > 1.0]
        return ev.real.max() if len(ev) else np.nan

    Rs = np.arange(rt_range[0], rt_range[1], 10.0)
    vals = np.array([re_part(R) for R in Rs])
    idx = [k for k in np.flatnonzero(np.sign(vals[:-1]) != np.sign(vals[1:]))
           if np.isfinite(vals[k]) and np.isfinite(vals[k + 1])]
    if not idx:
        return None
    Rh = brentq(re_part, Rs[idx[-1]], Rs[idx[-1] + 1])
    v = [v for v in equilibria(g, Rh, rL) if np.sign(v) == side][0]
    ev, _ = stability(g, v, Rh, rL, **comps)
    return Rh, v, 2 * np.pi / np.abs(ev.imag).max()


# ---- ideal model: integrator ---------------------------------------------------
def rhs(y, Rt, rL, g):
    """y has shape (3, n): rows v1, v2, iL; Rt shape (n,)."""
    v1, v2, iL = y
    ir = (v2 - v1) / Rt
    return np.stack([(ir - g(v1)) / C1, (-ir - iL) / C2, (v2 - rL * iL) / L])


def rk4_step(y, dt, Rt, rL, g):
    k1 = rhs(y, Rt, rL, g)
    k2 = rhs(y + 0.5 * dt * k1, Rt, rL, g)
    k3 = rhs(y + 0.5 * dt * k2, Rt, rL, g)
    k4 = rhs(y + dt * k3, Rt, rL, g)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def integrate(rpot, v1_0, t_end, dt, rL, g, keep=True, t_skip=0.0, y0=None):
    """
    Integrate a batch from (v1, v2, iL) = (v1_0, 0, 0), or from the full
    initial state y0 of shape (3, n) when given. rpot and v1_0 broadcast to
    the same length. An element with .knots runs in the compiled kernel, any
    other callable in the vectorised numpy RK4.

    keep=True  -> (t, Y) with Y of shape (nt, 3, n) for t >= t_skip.
    keep=False -> the list of v1-maxima per column found after t_skip.
    """
    if not np.isfinite([t_end, dt, t_skip, rL]).all() or dt <= 0 or not 0 <= t_skip < t_end or rL < 0:
        raise ValueError('require dt > 0, 0 <= t_skip < t_end and rL >= 0')
    rpot = np.atleast_1d(np.asarray(rpot, float))
    Rt = R0 + rpot
    n = len(rpot)
    if y0 is None:
        y = np.zeros((3, n))
        y[0] = np.broadcast_to(np.asarray(v1_0, float), rpot.shape)
    else:
        y = np.array(y0, float).reshape(3, n)
    nsteps = int(round(t_end / dt))
    n_skip = int(round(t_skip / dt))

    if hasattr(g, 'knots'):
        vk, ik = g.knots
        if keep:
            out = np.stack([kern.trajectory(y[:, j], Rt[j], rL, C1, C2, L, vk, ik, dt, nsteps, n_skip)
                            for j in range(n)], axis=2)
            return dt * np.arange(n_skip, nsteps + 1), out
        return [kern.hold(y[:, j], Rt[j], rL, C1, C2, L, vk, ik, dt, n_skip, nsteps - n_skip)[1]
                for j in range(n)]

    if keep:
        out = np.empty((nsteps - n_skip + 1, 3, n))
        out[0] = y
        j = 1
        for k in range(1, nsteps + 1):
            y = rk4_step(y, dt, Rt, rL, g)
            if k == n_skip:
                out[0] = y
                j = 1
            elif k > n_skip:
                out[j] = y
                j += 1
        return dt * np.arange(n_skip, nsteps + 1), out

    maxima = [[] for _ in range(n)]
    prev2 = y[0].copy()
    prev1 = y[0].copy()
    for k in range(1, nsteps + 1):
        y = rk4_step(y, dt, Rt, rL, g)
        cur = y[0]
        if k > n_skip:
            for i in np.flatnonzero((prev1 > prev2) & (prev1 >= cur)):
                maxima[i].append(prev1[i])
        prev2, prev1 = prev1, cur.copy()
    return maxima


def midpoint(v1, v2, rpot):
    """Node 3 (CH3): between R0 on the C2 side and Rpot on the C1 side."""
    return v2 + (v1 - v2) * R0 / (R0 + rpot)


def sweep_starts(g, rpot, rL, side=-1, eps=0.05):
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
        Rt = R0 + r
        eqs = [v for v in equilibria(g, Rt, rL) if np.sign(v) == side]
        v = eqs[0] if eqs else 0.0
        iL = v / (Rt + rL)
        fwd[:, k] = (v + eps, rL * iL, iL)
    back = np.zeros((3, len(rpot)))
    back[0], back[1] = 7.0, 6.0
    return fwd, back


def _hold(y, Rt, rL, g, dt, n_settle, n_collect):
    """Integrate a batch (3, n) at fixed Rt: settle, then collect the v1 maxima."""
    if hasattr(g, 'knots'):
        vk, ik = g.knots
        res = [kern.hold(y[:, j], Rt[j], rL, C1, C2, L, vk, ik, dt, n_settle, n_collect)
               for j in range(y.shape[1])]
        return np.column_stack([r[0] for r in res]), [list(r[1]) for r in res]
    for _ in range(n_settle):
        y = rk4_step(y, dt, Rt, rL, g)
    mx = [[] for _ in range(y.shape[1])]
    prev2 = y[0].copy()
    prev1 = y[0].copy()
    for _ in range(n_collect):
        y = rk4_step(y, dt, Rt, rL, g)
        cur = y[0]
        for i in np.flatnonzero((prev1 > prev2) & (prev1 >= cur)):
            mx[i].append(prev1[i])
        prev2, prev1 = prev1, cur.copy()
    return y, mx


def sweep_continuation(g, rsw, rL, dt, t_settle, t_collect, direction):
    """
    Maxima of v1 along a sweep that carries the state from one resistance to
    the next. 'down' starts on the negative outer equilibrium at the highest
    R (the bench's forward sweep); 'up' starts on the large outer cycle at the
    lowest R (the back sweep). Returns the list of maxima per R, in the order
    of `rsw`.
    """
    if direction not in ('down', 'up'):
        raise ValueError('direction must be down or up')
    if dt <= 0 or t_settle < 0 or t_collect <= 0:
        raise ValueError('invalid sweep integration times')
    rsw = np.asarray(rsw, float)
    if rsw.ndim != 1 or not len(rsw) or not np.isfinite(rsw).all():
        raise ValueError('need a nonempty finite resistance vector')
    order = np.argsort(rsw)[::-1] if direction == 'down' else np.argsort(rsw)
    seq = rsw[order]
    fwd, back = sweep_starts(g, seq[:1], rL)
    y = fwd if direction == 'down' else back
    n_settle, n_collect = int(round(t_settle / dt)), int(round(t_collect / dt))
    out = [None] * len(seq)
    for index, r in zip(order, seq):
        y, mx = _hold(y, np.array([R0 + r]), rL, g, dt, n_settle, n_collect)
        out[index] = mx[0]
    return out


# ---- bench model ---------------------------------------------------------------
def bench_params(path=os.path.join(HERE, 'identified.json'), tau_a=TAU_A, tau_b=TAU_B, slew=SLEW):
    """
    The bench kernel's parameter vector from identified.json. The small-signal
    capacitance identify.py measures on node 1 is C1 + tau_A AA/RA, so C1 is
    that value minus the lag's share.
    """
    with open(path) as fh:
        d = json.load(fh)
    k, ray = d['diode']['kennedy'], d['inductor']['rayleigh']
    c1_app = d['C1_nF']['small_signal_loop'] * 1e-9
    c1 = c1_app - tau_a * k['AA'] / k['RA_ohm']
    return kern.bench_params(c1, d['C2_nF'] * 1e-9, ray['L0_H'], ray['nu_H_per_A'], ray['r0_ohm'],
                             ray['rho_ohm_per_A'], k['RA_ohm'], k['AA'], k['RB_ohm'], k['AB'],
                             k['VpA_V'], k['VnA_V'], k['VpB_V'], k['VnB_V'], slew, tau_a, tau_b)


def bench_static_g(P):
    """The bench diode with the op-amps infinitely fast: a five-segment law."""
    bp, G = kern.static_pwl(P)
    return pwl_g(bp, G)


def bench_small_signal(P):
    """(c1_apparent, c2, l0, r0) of the bench model's linearisation."""
    return P[0] + P[15] * P[7] / P[6], P[1], P[2], P[4]


def bench_hopf(P, side=-1):
    c1, c2, l0, r0 = bench_small_signal(P)
    return hopf_point(bench_static_g(P), r0, side, c1=c1, c2=c2, l=l0)


def bench_sweep(P, rsw, dt, t_settle, t_collect, direction):
    """Continuation sweep of the bench model (see sweep_continuation)."""
    if direction not in ('down', 'up'):
        raise ValueError('direction must be down or up')
    rsw = np.asarray(rsw, float)
    order = np.argsort(rsw)[::-1] if direction == 'down' else np.argsort(rsw)
    g = bench_static_g(P)
    r0 = P[4]
    r_first = rsw[order[0]]
    if direction == 'down':
        eqs = [v for v in equilibria(g, R0 + r_first, r0) if v < 0]
        v = eqs[0] if eqs else 0.0
        iL = v / (R0 + r_first + r0)
        y = kern.bench_state(v + 0.05, r0 * iL, iL, P)
    else:
        y = kern.bench_state(7.0, 6.0, 0.0, P)
    n_settle, n_collect = int(round(t_settle / dt)), int(round(t_collect / dt))
    out = [None] * len(rsw)
    for k in order:
        y, mx = kern.bench_hold(y, R0 + rsw[k], dt, n_settle, n_collect, P)
        out[k] = list(mx)
    return out


def bench_runs(P, rpots, states, t_end, dt, t_skip):
    """(t, Y) with Y (nt, 5, n) for runs from the (3, n) states (v1, v2, iL) at each rpot."""
    nsteps, n_skip = int(round(t_end / dt)), int(round(t_skip / dt))
    Y = np.stack([kern.bench_trajectory(kern.bench_state(*states[:, j], P), R0 + r, dt, nsteps, n_skip, P)
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
    from scope_data import read_scope
    t, d = read_scope(path)
    return t[:n_max], d[:n_max, 0], d[:n_max, 1]


def load_measured_bifurcation(folder):
    """(rpot, max_v, sweep) from the *_bifurcation_points.csv sidecars."""
    rows = []
    for f in glob.glob(os.path.join(folder, '*_bifurcation_points.csv')):
        source_sweep = os.path.basename(f)[:-len('_bifurcation_points.csv')]
        source_dir = os.path.join(folder, source_sweep)
        if not os.path.isdir(source_dir):
            continue
        with open(f) as fh:
            head = fh.readline().strip().split(',')
            try:
                ir, im, isw, ifile = (head.index('rpot_ohm'), head.index('max_v'),
                                      head.index('sweep'), head.index('filename'))
            except ValueError:
                continue        # not a bifurcation.py sidecar
            for line in fh:
                p = line.strip().split(',')
                try:
                    r, m = float(p[ir]), float(p[im])
                    if (p[isw] == source_sweep and np.isfinite([r, m]).all() and 0 <= r <= 1000
                            and os.path.isfile(os.path.join(source_dir, p[ifile]))):
                        rows.append((r, m, p[isw]))
                except (ValueError, IndexError):
                    pass
    return rows


# ---- regime metrics ----------------------------------------------------------------
def describe(m):
    """One-line regime label for a list of v1 maxima."""
    m = np.asarray(m, float)
    if len(m) == 0:
        return 'no maxima (settled)'
    k = len(np.unique(np.round(m, 2)))
    lo, hi = m.min(), m.max()
    if hi > 5:
        return f'large cycle, max {hi:4.2f}'
    if hi > 0.9 > -0.3 > lo:
        return f'double scroll, {k:3d} lvls'
    return f'{k:3d} lvls  [{lo:5.2f},{hi:5.2f}]'


def metrics(rsw, max_down, max_up, persist=4):
    """
    The transition resistances of the two sweeps, each required to hold for
    `persist` consecutive steps so that the slow transients just below the
    Hopf point do not count.
    """
    rsw = np.asarray(rsw, float)
    lvls = lambda m: len(np.unique(np.round(m, 2))) if len(m) else 0
    is_lc = lambda m: len(m) > 0 and min(m) > 5
    is_ds = lambda m: len(m) > 0 and max(m) > 0.9 and min(m) < -0.3 and not is_lc(m)
    is_p1 = lambda m: len(m) > 0 and max(m) < 0.9 and np.ptp(m) < 0.02
    is_split = lambda m: len(m) > 0 and max(m) < 0.9 and lvls(m) >= 2 and np.ptp(m) >= 0.03
    is_chaos = lambda m: len(m) > 0 and max(m) < 0.9 and lvls(m) > 8

    def first_persisting(cond, seq, direction):
        idx = np.argsort(rsw)[::-1] if direction == 'down' else np.argsort(rsw)
        flags = [cond(seq[k]) for k in idx]
        for j in range(len(idx) - persist + 1):
            if all(flags[j:j + persist]):
                return float(rsw[idx[j]])
        return None

    def r1(seq):
        idx = np.argsort(rsw)[::-1]
        seen = False
        flags = [is_split(seq[k]) or is_chaos(seq[k]) or is_ds(seq[k]) for k in idx]
        for j, k in enumerate(idx):
            if is_p1(seq[k]):
                seen = True
            elif seen and all(flags[j:j + persist]):
                return float(rsw[k])
        return None
    ds_down = [rsw[k] for k in range(len(rsw)) if is_ds(max_down[k])]
    lc_up = [rsw[k] for k in range(len(rsw)) if is_lc(max_up[k])]
    return {
        'first doubling R1 (down)': r1(max_down),
        'chaos onset (down)': first_persisting(is_chaos, max_down, 'down'),
        'double scroll from (down)': first_persisting(is_ds, max_down, 'down'),
        'double scroll ends (down)': float(min(ds_down)) if ds_down else None,
        'large cycle from (down)': first_persisting(is_lc, max_down, 'down'),
        'large cycle survives to (up)': float(max(lc_up)) if lc_up else None,
        'double scroll from (up)': first_persisting(is_ds, max_up, 'up'),
    }


def period_of(t, v1):
    """Mean spacing of the maxima of v1(t), in seconds (nan if fewer than three)."""
    pk = np.flatnonzero((v1[1:-1] > v1[:-2]) & (v1[1:-1] >= v1[2:])) + 1
    return float(np.mean(np.diff(t[pk]))) if len(pk) > 2 else np.nan


# ---- figures ---------------------------------------------------------------------
def plot_nr(g, rpots, path, title):
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
        ax.plot(v, -v / (R0 + r) * 1e3, lw=0.9, label=f'load line, Rpot = {r:g} ohm')
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
            r_meas, path_meas, label = rec
            _, v1, v2 = load_record(path_meas, n_meas)
            ax.plot(v1, v2, lw=0.3, color='C3')
            ax.set_title(f'measured: {label}, Rpot = {r_meas:.1f} ohm', fontsize=9)
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


# ---- main ------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--model', choices=['ideal', 'static', 'bench'], default='ideal')
    p.add_argument('--r', type=float, nargs='+', default=[806.0, 770.0, 745.0, 700.0, 650.0, 450.0, 320.0],
                   help='potentiometer values for portraits/time series (ohm)')
    p.add_argument('--i0', default='auto',
                   help='ideal: current offset (mA) subtracted from the V-I fit; "auto" removes the fitted g(0)')
    p.add_argument('--bp-scale', type=float, default=1.0, help='ideal: scale of the two inner breakpoints')
    p.add_argument('--symmetric', action='store_true', help='ideal: average the shoulders and inner breakpoints')
    p.add_argument('--c1', type=float, default=C1 * 1e9,
                   help=f'ideal: C1 in nF (default {C1 * 1e9:g}, calibrated; nominal {C1_NOMINAL * 1e9:g})')
    p.add_argument('--c2', type=float, default=C2 * 1e9, help='ideal: C2 in nF')
    p.add_argument('--l', type=float, default=L * 1e3, help='ideal: L in mH')
    p.add_argument('--rl', type=float, default=0.0, help='ideal: inductor series resistance (ohm)')
    p.add_argument('--r0', type=float, default=R0, help='fixed part of the coupling resistance (ohm)')
    p.add_argument('--identified', default=os.path.join(HERE, 'identified.json'), help='bench: identify.py output')
    p.add_argument('--tau-b', type=float, default=TAU_B * 1e6, help='bench: stage-B op-amp lag (us)')
    p.add_argument('--slew', type=float, default=SLEW * 1e-6, help='bench: op-amp slew rate (V/us)')
    p.add_argument('--dt', type=float, default=None, help='RK4 step (s); default 0.5 us ideal, 0.1 us bench')
    p.add_argument('--t', type=float, default=0.1, help='length of each portrait run (s)')
    p.add_argument('--skip', type=float, default=0.06, help='transient discarded before plotting (s)')
    p.add_argument('--sweep', type=float, nargs=3, default=[300, 900, 1], metavar=('RMIN', 'RMAX', 'STEP'))
    p.add_argument('--protocol', choices=['continuation', 'seeded'], default='continuation',
                   help='continuation: carry the state along (default); seeded: restart every R (ideal only)')
    p.add_argument('--sweep-t', type=float, default=0.1, help='seeded: length of each run (s)')
    p.add_argument('--settle', type=float, default=0.02, help='continuation: seconds discarded after each R step')
    p.add_argument('--collect', type=float, default=0.03, help='continuation: seconds of maxima kept at each R')
    p.add_argument('--no-sweep', action='store_true')
    p.add_argument('--records', default=HERE, help='folder with the sweep folders and their _rpot.csv sidecars')
    p.add_argument('--n-meas', type=int, default=60000, help='samples of each measured record to draw')
    p.add_argument('--export', action='store_true', help='write scope-format CSVs of the --r runs')
    p.add_argument('--tag', default='', help='suffix for the output names')
    p.add_argument('--out', default=HERE, help='where the figures go')
    a = p.parse_args()
    if a.protocol == 'seeded' and a.model != 'ideal':
        p.error('--protocol seeded is only available for the ideal model')
    dt = a.dt or (BENCH_DT if a.model == 'bench' else 0.5e-6)
    set_components(a.c1 * 1e-9, a.c2 * 1e-9, a.l * 1e-3, a.r0)

    if a.model == 'ideal':
        i0 = fitted_offset() if a.i0 == 'auto' else float(a.i0) * 1e-3
        g = make_g(i0, a.bp_scale, a.symmetric)
        rL = a.rl
        comps = {}
        print(f'ideal model: C1 = {C1 * 1e9:g} nF, C2 = {C2 * 1e9:g} nF, L = {L * 1e3:g} mH, '
              f'R0 = {R0:g} ohm, rL = {rL:g} ohm')
        print(f'N_R: M1 fit (diode_fit.json), offset removed = {i0 * 1e3:.3f} mA, breakpoint scale = {a.bp_scale:g}'
              + (', shoulders symmetrised' if a.symmetric else ''))
        subtitle = f'ideal, C1 = {C1 * 1e9:g} nF, rL = {rL:g} ohm' + (', symmetric' if a.symmetric else '') \
            + (f', bp scale {a.bp_scale:g}' if a.bp_scale != 1 else '') + f', {a.protocol}'
    elif a.model == 'static':
        P = bench_params(a.identified, TAU_A, a.tau_b * 1e-6, a.slew * 1e6)
        c1a, c2b, l0, r0b = bench_small_signal(P)
        set_components(c1a, c2b, l0)
        g = bench_static_g(P)
        rL = r0b
        comps = {}
        print(f'static model ({os.path.basename(a.identified)}): C1 = {C1 * 1e9:.2f} nF, C2 = {C2 * 1e9:.1f} nF, '
              f'L = {L * 1e3:.2f} mH, rL = {rL:.2f} ohm, R0 = {R0:g} ohm; the identified five-segment element')
        subtitle = 'static model: identified circuit with constant components and a static element'
    else:
        P = bench_params(a.identified, TAU_A, a.tau_b * 1e-6, a.slew * 1e6)
        g = bench_static_g(P)
        c1a, c2b, l0, r0b = bench_small_signal(P)
        rL = r0b
        comps = dict(c1=c1a, c2=c2b, l=l0)
        print(f'bench model ({os.path.basename(a.identified)}): C1 = {P[0] * 1e9:.2f} nF (+ {(c1a - P[0]) * 1e9:.2f} nF '
              f'from the stage-A lag), C2 = {P[1] * 1e9:.1f} nF, L = {P[2] * 1e3:.1f} mH + {P[3]:.2f} H/A |iL|, '
              f'r = {P[4]:.1f} ohm + {P[5]:.0f} ohm/A |iL|, R0 = {R0:g} ohm')
        print(f'N_R: Kennedy diode RA = {P[6]:.0f} ohm, AA = {P[7]:.3f}, RB = {P[8]:.0f} ohm, AB = {P[9]:.2f}, '
              f'rails A {P[10]:+.2f}/{P[11]:+.2f} V, B {P[12]:+.2f}/{P[13]:+.2f} V, slew {P[14] * 1e-6:g} V/us, '
              f'lags {P[15] * 1e6:g}/{P[16] * 1e6:g} us')
        subtitle = f'bench model, tau_B = {P[16] * 1e6:g} us, slew {P[14] * 1e-6:g} V/us'
    print(f'{"v_lo":>8} {"v_hi":>8} {"G (mS)":>9} {"c (mA)":>9}')
    for lo, hi, G, c in g.segments:
        print(f'{lo:8.3f} {hi:8.3f} {G * 1e3:9.4f} {c * 1e3:9.4f}')

    print('\nequilibria, g(v1) = -v1/(Rt + rL):')
    for r in a.r:
        Rt = R0 + r
        parts = []
        for v in equilibria(g, Rt, rL):
            ev, G = stability(g, v, Rt, rL, **comps)
            parts.append(f'v1 = {v:6.2f} V, G = {G * 1e3:6.3f} mS, {"UNSTABLE" if ev.real.max() > 0 else "stable"}')
        print(f'   Rpot = {r:6.1f}:  ' + '  |  '.join(parts))
    hopf = bench_hopf(P) if a.model == 'bench' else hopf_point(g, rL)
    if a.model == 'static':
        dt = a.dt or 0.5e-6
    if hopf:
        Rh, vh, Th = hopf
        print(f'\nHopf point of the negative outer equilibrium: Rpot = {Rh - R0:.0f} ohm, '
              f'v1* = {vh:.2f} V, period {Th * 1e6:.0f} us')

    # ---- single runs: from the negative outer equilibrium and from the large cycle ----
    rp = np.repeat(a.r, 2)
    y_fwd, y_back = sweep_starts(g, a.r, rL)
    y0 = np.empty((3, 2 * len(a.r)))
    y0[:, 0::2], y0[:, 1::2] = y_fwd, y_back
    print(f'\nintegrating {len(rp)} runs of {a.t} s at dt = {dt} ...')
    if a.model == 'bench':
        t, Y = bench_runs(P, rp, y0, a.t, dt, a.skip)
    else:
        t, Y = integrate(rp, None, a.t, dt, rL, g, keep=True, t_skip=a.skip, y0=y0)
    runs = [(r, t, Y[:, :, 2 * k], Y[:, :, 2 * k + 1]) for k, r in enumerate(a.r)]

    print(f'{"Rpot":>7}  {"start":>6}  {"v1 min":>7}  {"v1 max":>7}  {"v2 min":>7}  {"v2 max":>7}  '
          f'{"n maxima":>8}  {"distinct":>8}  {"period":>8}')
    for r, _, yp, ym in runs:
        for tag, y in [('eq', yp), ('cycle', ym)]:
            v1 = y[:, 0]
            pk = v1[1:-1][(v1[1:-1] > v1[:-2]) & (v1[1:-1] >= v1[2:])]
            T = period_of(t, v1)
            print(f'{r:7.1f}  {tag:>6}  {v1.min():7.3f}  {v1.max():7.3f}  {y[:, 1].min():7.3f}  {y[:, 1].max():7.3f}  '
                  f'{len(pk):8d}  {len(np.unique(np.round(pk, 2))):8d}  {T * 1e6:6.1f} us')

    records = find_records(a.records)
    for r in a.r:
        rec = nearest_record(records, r)
        print(f'measured record for Rpot = {r:g} ohm: ' + (f'{rec[2]} at {rec[0]:.1f} ohm' if rec else 'none within 5 ohm'))
    plot_nr(g, a.r, os.path.join(a.out, f'simulated_nr{a.tag}.png'),
            'Nonlinear element of the ' + ('bench model (op-amps infinitely fast)' if a.model == 'bench' else 'ideal model (M1 fit)'))
    plot_portraits(runs, records, os.path.join(a.out, f'simulated_portraits{a.tag}.png'), a.n_meas)
    plot_timeseries(runs, os.path.join(a.out, f'simulated_timeseries{a.tag}.png'))

    if a.export:
        d = os.path.join(a.out, f'simulated{a.tag}')
        os.makedirs(d, exist_ok=True)
        for r, t, yp, ym in runs:
            for tag, y in [('plus', yp), ('minus', ym)]:
                v3 = midpoint(y[:, 0], y[:, 1], r)
                f = os.path.join(d, f'sim {r:.1f} ohm {tag}.csv')
                np.savetxt(f, np.column_stack([t - t[0], y[:, 0], y[:, 1], v3]), delimiter=',', fmt='%.6e',
                           header='Time(s),CH1(V),CH2(V),CH3(V)', comments='')
        print(f'\nwrote {2 * len(runs)} records to {d}')

    # ---- bifurcation sweep ----
    if not a.no_sweep:
        rmin, rmax, step = a.sweep
        rsw = np.arange(rmin, rmax + step / 2, step)
        if a.protocol == 'continuation':
            print(f'\nsweeping {len(rsw)} values of Rpot down from the negative equilibrium and up from the '
                  f'large cycle, carrying the state along ({a.settle * 1e3:.0f} ms settle + '
                  f'{a.collect * 1e3:.0f} ms collect per step) ...')
            if a.model == 'bench':
                mf = bench_sweep(P, rsw, dt, a.settle, a.collect, 'down')
                mb = bench_sweep(P, rsw, dt, a.settle, a.collect, 'up')
            else:
                mf = sweep_continuation(g, rsw, rL, dt, a.settle, a.collect, 'down')
                mb = sweep_continuation(g, rsw, rL, dt, a.settle, a.collect, 'up')
            tags = ('sim, forward: R turned down from the negative equilibrium',
                    'sim, back: R turned up from the large cycle')
        else:
            print(f'\nsweeping {len(rsw)} values of Rpot x 2 starts, {a.sweep_t} s each ...')
            y_fwd, y_back = sweep_starts(g, rsw, rL)
            mf = integrate(rsw, None, a.sweep_t, dt, rL, g, keep=False, t_skip=a.skip, y0=y_fwd)
            mb = integrate(rsw, None, a.sweep_t, dt, rL, g, keep=False, t_skip=a.skip, y0=y_back)
            tags = ('sim, forward (each R restarted on the negative equilibrium)',
                    'sim, back (each R restarted on the large cycle)')
        measured = load_measured_bifurcation(a.records)
        plot_bifurcation(rsw, mf, mb, measured, os.path.join(a.out, f'simulated_bifurcation{a.tag}.png'),
                         subtitle=subtitle, tags=tags)
        print(f'\n{"Rpot":>7}  {"down (forward)":>24}  {"up (back)":>24}')
        for r, m1, m2 in zip(rsw, mf, mb):
            if abs(r - round(r / 25) * 25) > step / 2:
                continue
            print(f'{r:7.1f}  {describe(m1):>24}  {describe(m2):>24}')
        print('\ntransitions (ohm):')
        for k, v in metrics(rsw, mf, mb).items():
            print(f'   {k:32s} {v if v is None else f"{v:.0f}"}')

    print('\nfigures written to', a.out)


if __name__ == '__main__':
    main()
