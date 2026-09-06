"""
Simulate the Chua oscillator of Fig. 2 with the MEASURED nonlinear element.

The circuit (K1 closed, K2 in CHAOS): C1 in parallel with N_R at node 1,
L and C2 at node 2, and R0 + Rpot between the two nodes. With v1 = v_C1,
v2 = v_C2 and iL the inductor current flowing from node 2 to ground:

    C1 dv1/dt = (v2 - v1)/Rt - g(v1)
    C2 dv2/dt = (v1 - v2)/Rt - iL           Rt = R0 + Rpot
    L  diL/dt = v2 - rL iL                   rL = inductor series resistance

g(v) is the current into N_R. It is NOT the textbook three-segment odd
function: it is the five-segment piecewise-linear fit that find_breakpoints.py
found on trace1.csv (V-I measurement, K2 in VI), copied verbatim into
NR_SEGMENTS below. Two things about that fit matter for the simulation and
are exposed as options:

  --shunt   find_breakpoints.py plots the current against CH1 = V_V, the
            SOURCE side of the 216 ohm shunt. The voltage across N_R is
            V_I = V_V - Rs*i (Fig. 2 caption). --shunt 216 re-expresses the
            fitted lines against V_I; --shunt 0 (default) takes the fit
            exactly as printed.
  --i0      the fitted curve passes through g(0) = -0.51 mA, i.e. the
            whole curve sits 0.5 mA low (0.11 V of CH1-CH2 offset, about one
            scope quantum). A Chua diode passes no current at v = 0, and the
            measured oscillator lives further on the NEGATIVE side, which a
            downward offset cannot produce. --i0 auto removes it; the
            default keeps it.

With the fit taken literally the load line -v/Rt misses the negative-slope
segments for every Rpot above ~600 ohm: the only equilibrium is in the
saturated (positive-slope) branch near +5.4 V and it is stable, so the
simulated circuit sits at DC where the bench shows limit cycles and chaos.
The fitted segments also do not quite meet at the breakpoints (steps of up
to 0.2 mA); they are bridged linearly over --blend volts so the integrator
sees a continuous function. Beyond the fitted range (-9.33 V .. 8.28 V,
where the trace saturates) the outer segments are extended, which is the
physically right thing: there the op-amps are railed and N_R is a plain
resistor.

Integration is a fixed-step RK4 vectorised over a batch of resistances, so a
551-point bifurcation sweep costs about as much as a handful of single runs.
Every R is run twice, from v1(0) = +1 V and from v1(0) = -1 V, because the
single-scroll attractors around the two outer equilibria coexist; that is the
simulated stand-in for the forward/back hysteresis of the measured sweeps.

Channel convention: in the records CH1 = v_C1 (the large swing, across N_R)
and CH2 = v_C2, as rpot.py documents; the labels in the Fig. 2 caption are
the other way round.

Outputs, written beside this script (--tag adds a suffix to the names):

    simulated_nr.png            g(v) with the load lines -v/Rt of the chosen R
    simulated_portraits.png     v2 vs v1 for each --r: both starts and, when a
                                "<...> <R> ohm.csv" record exists, the scope
    simulated_timeseries.png    v1(t) for each --r
    simulated_bifurcation.png   maxima of v1 against Rpot, measured sweep
                                points (from bifurcation.py) underneath
    simulated/                  with --export: one scope-format CSV per --r
                                (Time, CH1 = v1, CH2 = v2, CH3 = midpoint) so
                                batch_rpot / lorenz_map / lyapunov run on it

Usage:
    python simulate.py                                  # fit taken literally
    python simulate.py --shunt 216 --i0 auto --tag _corrected
    python simulate.py --shunt 216 --i0 auto --rl 30 --sweep 300 900 0.5
    python simulate.py --r 773 716.9 656.7 316.1 --export --t 0.2
"""
import argparse
import glob
import os
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- circuit values from Fig. 2 -------------------------------------------
C1 = 10e-9        # F
C2 = 100e-9       # F
L = 18e-3         # H
R0 = 990.0        # ohm, fixed part of the coupling resistance
RS_BENCH = 216.0  # ohm, the shunt of the V-I measurement

# ---- N_R: output of find_breakpoints.py on ../trace1.csv (R^2 = 0.9757) ----
# (V_lo, V_hi, slope [S], intercept [A]); i = slope * V + intercept, V = CH1
NR_SEGMENTS = [
    (-9.330, -6.164,  2.0483e-03,  1.5342e-02),
    (-6.164, -1.142, -4.6218e-04, -4.6379e-05),
    (-1.142,  0.856, -8.9042e-04, -5.1227e-04),
    ( 0.856,  5.371, -4.7317e-04, -8.4943e-04),
    ( 5.371,  8.280,  2.1063e-03, -1.4491e-02),
]
V_EXTEND = 40.0   # V, how far the outer segments are extrapolated


# ---- the nonlinear element --------------------------------------------------
def nr_segments_v(segments=NR_SEGMENTS, shunt=0.0, i0=0.0):
    """
    The fitted segments re-expressed against the voltage ACROSS N_R.

    Each fitted line i = m*V_V + b becomes, with v = V_V - Rs*i,

        i = G*v + c,   G = m / (1 - Rs*m),   c = (b - i0) / (1 - Rs*m)

    and its end points move by Rs*i. `i0` is subtracted from the fit first.
    shunt = 0 and i0 = 0 reproduce the fit exactly as printed.

    Returns [(v_lo, v_hi, G, c), ...]. Because the fitted segments do not
    meet at the breakpoints, neighbouring ranges may overlap or leave a gap
    by a few tens of mV once the shunt is accounted for.
    """
    out = []
    for lo, hi, m, b in segments:
        k = 1.0 - shunt * m
        G, c = m / k, (b - i0) / k
        i_lo, i_hi = m * lo + b - i0, m * hi + b - i0
        out.append((lo - shunt * i_lo, hi - shunt * i_hi, G, c))
    return out


def nr_knots(segments_v, blend=0.05, extend=V_EXTEND):
    """
    Knots (v, i) of the continuous piecewise-linear g(v).

    Between two segments the breakpoint is the mean of the two fitted end
    points (identical to the fit when shunt = 0), and the interval of width
    `blend` around it is a straight bridge from one line to the next, so the
    ~0.1 mA steps of the raw fit become steep but finite. The outer segments
    are extended to +-`extend` volts.
    """
    n = len(segments_v)
    vs, cs = [], []
    for k, (lo, hi, G, c) in enumerate(segments_v):
        a = -extend if k == 0 else 0.5 * (segments_v[k - 1][1] + lo) + blend / 2
        z = extend if k == n - 1 else 0.5 * (hi + segments_v[k + 1][0]) - blend / 2
        vs += [a, z]
        cs += [G * a + c, G * z + c]
    vs, cs = np.array(vs), np.array(cs)
    if np.any(np.diff(vs) <= 0):
        raise ValueError('N_R knots not monotonic; increase --blend')
    return vs, cs


def make_g(blend=0.05, shunt=0.0, i0=0.0):
    segs = nr_segments_v(NR_SEGMENTS, shunt, i0)
    vk, ik = nr_knots(segs, blend=blend)

    def g(v):
        return np.interp(v, vk, ik)
    g.knots = (vk, ik)
    g.segments = segs
    return g


def fitted_offset():
    """g(0) of the raw fit: the intercept of the segment containing V = 0."""
    for lo, hi, m, b in NR_SEGMENTS:
        if lo <= 0 <= hi:
            return b
    return 0.0


# ---- equilibria ------------------------------------------------------------
def equilibria(g, Rt, rL):
    """
    Fixed points v1 of the circuit: iL = (v1 - v2)/Rt and v2 = rL*iL, so
    g(v1) = -v1/(Rt + rL), the load line through the origin.
    """
    from scipy.optimize import brentq
    h = lambda v: g(v) + v / (Rt + rL)
    vs = np.linspace(-15, 15, 30001)
    hv = h(vs)
    idx = np.flatnonzero(hv[:-1] * hv[1:] < 0)
    return [brentq(h, vs[k], vs[k + 1]) for k in idx]


def stability(g, v1, Rt, rL, dv=1e-4):
    """(Jacobian eigenvalues, local slope G of N_R) at the fixed point v1."""
    G = (g(v1 + dv) - g(v1 - dv)) / (2 * dv)
    J = np.array([[(-1 / Rt - G) / C1, 1 / (Rt * C1), 0.0],
                  [1 / (Rt * C2), -1 / (Rt * C2), -1 / C2],
                  [0.0, 1 / L, -rL / L]])
    return np.linalg.eigvals(J), G


# ---- integrator ------------------------------------------------------------
def rhs(y, Rt, rL, g):
    """y has shape (3, n): rows v1, v2, iL; Rt shape (n,)."""
    v1, v2, iL = y
    ir = (v2 - v1) / Rt
    return np.stack([(ir - g(v1)) / C1,
                     (-ir - iL) / C2,
                     (v2 - rL * iL) / L])


def rk4_step(y, dt, Rt, rL, g):
    k1 = rhs(y, Rt, rL, g)
    k2 = rhs(y + 0.5 * dt * k1, Rt, rL, g)
    k3 = rhs(y + 0.5 * dt * k2, Rt, rL, g)
    k4 = rhs(y + dt * k3, Rt, rL, g)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def integrate(rpot, v1_0, t_end, dt, rL, g, keep=True, t_skip=0.0):
    """
    Integrate a batch from (v1, v2, iL) = (v1_0, 0, 0). rpot and v1_0
    broadcast to the same length.

    keep=True  -> returns (t, Y) with Y of shape (nt, 3, n) for t >= t_skip.
    keep=False -> returns the list of v1-maxima per column found after
                  t_skip (for the bifurcation diagram); nothing stored.
    """
    rpot = np.atleast_1d(np.asarray(rpot, float))
    v1_0 = np.broadcast_to(np.asarray(v1_0, float), rpot.shape)
    Rt = R0 + rpot
    n = len(rpot)
    y = np.zeros((3, n))
    y[0] = v1_0
    nsteps = int(round(t_end / dt))
    n_skip = int(round(t_skip / dt))

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
        t = t_skip + dt * np.arange(out.shape[0])
        return t, out

    maxima = [[] for _ in range(n)]
    prev2 = y[0].copy()
    prev1 = y[0].copy()
    for k in range(1, nsteps + 1):
        y = rk4_step(y, dt, Rt, rL, g)
        cur = y[0]
        if k > n_skip:
            hit = (prev1 > prev2) & (prev1 >= cur)
            for i in np.flatnonzero(hit):
                maxima[i].append(prev1[i])
        prev2, prev1 = prev1, cur.copy()
    return maxima


def midpoint(v1, v2, rpot):
    """Node 3 (CH3): between R0 on the C2 side and Rpot on the C1 side."""
    return v2 + (v1 - v2) * R0 / (R0 + rpot)


# ---- measured records ------------------------------------------------------
def find_records(folder):
    """{Rpot: path} for files named like 'chaos - 716.9 ohm.csv'."""
    out = {}
    for f in glob.glob(os.path.join(folder, '*.csv')):
        m = re.search(r'(\d+(?:\.\d+)?)\s*ohm', os.path.basename(f), re.I)
        if m:
            out[float(m.group(1))] = f
    return out


def load_record(path, n_max=None):
    d = np.genfromtxt(path, delimiter=',', skip_header=1, max_rows=n_max)
    return d[:, 0], d[:, 1], d[:, 2]     # t, CH1 = v_C1, CH2 = v_C2


def load_measured_bifurcation(folder):
    """(rpot, max_v, sweep) from the *_bifurcation_points.csv sidecars."""
    rows = []
    for f in glob.glob(os.path.join(folder, '*_bifurcation_points.csv')):
        with open(f) as fh:
            head = fh.readline().strip().split(',')
            ir, im, isw = (head.index('rpot_ohm'), head.index('max_v'),
                           head.index('sweep'))
            for line in fh:
                p = line.strip().split(',')
                try:
                    rows.append((float(p[ir]), float(p[im]), p[isw]))
                except (ValueError, IndexError):
                    pass
    return rows


# ---- figures ---------------------------------------------------------------
def plot_nr(g, rpots, path, shunt=0.0, i0=0.0):
    v = np.linspace(-10, 10, 2001)
    fig, ax = plt.subplots(figsize=(8, 5))
    raw = os.path.join(HERE, '..', 'trace1.csv')
    if os.path.exists(raw):
        d = np.genfromtxt(raw, delimiter=',', skip_header=1)
        i_m = (d[:, 1] - d[:, 2]) / RS_BENCH
        ax.plot(d[:, 1] - shunt * i_m, (i_m - i0) * 1e3, '.', ms=1,
                alpha=0.15, color='steelblue',
                label='measured (trace1.csv)'
                      + (', same corrections' if shunt or i0 else ''))
    ax.plot(v, g(v) * 1e3, 'crimson', lw=2, label='g(v): fitted segments')
    for (_, hi, _, _), (lo, _, _, _) in zip(g.segments[:-1], g.segments[1:]):
        ax.axvline(0.5 * (hi + lo), color='grey', ls='--', lw=0.7)
    for r in rpots:
        ax.plot(v, -v / (R0 + r) * 1e3, lw=0.9,
                label=f'load line, Rpot = {r:g} ohm')
    ax.set_xlabel('v across N_R  (V)')
    ax.set_ylabel('current into N_R  (mA)')
    what = (f'shunt {shunt:g} ohm, {i0 * 1e3:.2f} mA offset removed'
            if shunt or i0 else 'fit taken literally')
    ax.set_title(f'Nonlinear element used in the simulation ({what})')
    ax.set_ylim(-6, 6)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_portraits(runs, records, path, n_meas):
    """
    runs: list of (rpot, t, Yplus, Yminus) with Y of shape (nt, 3).
    One row per R: sim from +1 V, sim from -1 V, measured (if any).
    """
    nrow = len(runs)
    fig, axes = plt.subplots(nrow, 3, figsize=(12, 3.2 * nrow), squeeze=False)
    for row, (r, t, yp, ym) in enumerate(runs):
        for col, (y, tag) in enumerate([(yp, 'v1(0) = +1 V'),
                                        (ym, 'v1(0) = -1 V')]):
            ax = axes[row, col]
            ax.plot(y[:, 0], y[:, 1], lw=0.3, color='C0')
            ax.set_title(f'simulated, Rpot = {r:g} ohm, {tag}', fontsize=9)
        ax = axes[row, 2]
        rec = records.get(r)
        if rec:
            _, v1, v2 = load_record(rec, n_meas)
            ax.plot(v1, v2, lw=0.3, color='C3')
            ax.set_title(f'measured: {os.path.basename(rec)}', fontsize=9)
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
    fig, axes = plt.subplots(nrow, 1, figsize=(12, 1.9 * nrow), squeeze=False,
                             sharex=True)
    for ax, (r, t, yp, ym) in zip(axes[:, 0], runs):
        m = t <= t[0] + window
        ax.plot((t[m] - t[0]) * 1e3, yp[m, 0], lw=0.6, color='C0',
                label='v1(0) = +1 V')
        ax.plot((t[m] - t[0]) * 1e3, ym[m, 0], lw=0.6, color='C1',
                label='v1(0) = -1 V')
        ax.set_ylabel('v1 (V)')
        ax.set_title(f'Rpot = {r:g} ohm', fontsize=9, loc='left')
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=8, loc='upper right')
    axes[-1, 0].set_xlabel('time after transient (ms)')
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_bifurcation(rpot, max_plus, max_minus, measured, path):
    fig, ax = plt.subplots(figsize=(14, 7.5))
    if measured:
        sweeps = sorted({s for _, _, s in measured})
        for s, c in zip(sweeps, ['0.65', '0.8']):
            pts = np.array([(r, m) for r, m, sw in measured if sw == s])
            ax.plot(pts[:, 0], pts[:, 1], '.', ms=0.7, color=c, alpha=0.6,
                    label=f'measured, {s}', zorder=1)
    for maxima, c, tag in [(max_plus, 'C0', 'sim, v1(0) = +1 V'),
                           (max_minus, 'C3', 'sim, v1(0) = -1 V')]:
        xs = np.concatenate([np.full(len(m), r) for r, m in zip(rpot, maxima)])
        ys = np.concatenate([np.asarray(m, float) for m in maxima])
        ax.plot(xs, ys, '.', ms=0.45, color=c, alpha=0.8, label=tag, zorder=2)
    ax.set_xlabel('Rpot (ohm)')
    ax.set_ylabel('local maxima of v1 = v_C1 (V)')
    ax.set_title('Bifurcation diagram: simulated maxima of v1 against Rpot')
    ax.grid(alpha=0.3)
    ax.legend(markerscale=15, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=450)
    plt.close(fig)


def describe(m):
    """One-line regime label for a list of v1 maxima."""
    m = np.asarray(m, float)
    if len(m) == 0:
        return 'no maxima (settled)'
    k = len(np.unique(np.round(m, 2)))
    lo, hi = m.min(), m.max()
    if hi > 0 > lo:
        return f'double scroll, {k:3d} lvls'
    return f'{k:3d} lvls  [{lo:5.2f},{hi:5.2f}]'


# ---- main ------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--r', type=float, nargs='+',
                   default=[773.0, 764.6, 716.9, 666.0, 656.7, 316.1],
                   help='potentiometer values for portraits/time series (ohm); '
                        'defaults are the six recorded regimes')
    p.add_argument('--shunt', type=float, default=0.0,
                   help='shunt Rs (ohm) used to convert the fit from CH1 = V_V '
                        'to the voltage across N_R; 0 = fit taken literally, '
                        '216 = the bench value')
    p.add_argument('--i0', default='0',
                   help='current offset (mA) subtracted from the fit; '
                        '"auto" removes the fitted g(0)')
    p.add_argument('--rl', type=float, default=0.0,
                   help='inductor series resistance (ohm)')
    p.add_argument('--blend', type=float, default=0.05,
                   help='width (V) over which the fitted segments are joined')
    p.add_argument('--dt', type=float, default=0.5e-6, help='RK4 step (s)')
    p.add_argument('--t', type=float, default=0.06,
                   help='length of each portrait run (s)')
    p.add_argument('--skip', type=float, default=0.02,
                   help='transient discarded before plotting/maxima (s)')
    p.add_argument('--sweep', type=float, nargs=3, default=[300, 900, 1],
                   metavar=('RMIN', 'RMAX', 'STEP'),
                   help='bifurcation sweep of Rpot (ohm)')
    p.add_argument('--sweep-t', type=float, default=0.05,
                   help='length of each sweep run (s); --skip is discarded')
    p.add_argument('--no-sweep', action='store_true')
    p.add_argument('--records', default=HERE,
                   help='folder with "<name> - <R> ohm.csv" scope records')
    p.add_argument('--n-meas', type=int, default=60000,
                   help='samples of each measured record to draw')
    p.add_argument('--export', action='store_true',
                   help='write scope-format CSVs of the --r runs to simulated/')
    p.add_argument('--tag', default='', help='suffix for the output names')
    p.add_argument('--out', default=HERE, help='where the figures go')
    a = p.parse_args()

    i0 = fitted_offset() if a.i0 == 'auto' else float(a.i0) * 1e-3
    g = make_g(a.blend, a.shunt, i0)
    print(f'N_R model: shunt = {a.shunt:g} ohm, offset removed = '
          f'{i0 * 1e3:.3f} mA, rL = {a.rl:g} ohm')
    print(f'{"v_lo":>8} {"v_hi":>8} {"G (mS)":>9} {"c (mA)":>9}')
    for lo, hi, G, c in g.segments:
        print(f'{lo:8.3f} {hi:8.3f} {G * 1e3:9.4f} {c * 1e3:9.4f}')

    print('\nequilibria, g(v1) = -v1/(Rt + rL):')
    for r in a.r:
        Rt = R0 + r
        parts = []
        for v in equilibria(g, Rt, a.rl):
            ev, G = stability(g, v, Rt, a.rl)
            re = ev.real.max()
            parts.append(f'v1 = {v:6.2f} V, G = {G * 1e3:6.3f} mS, '
                         f'{"UNSTABLE" if re > 0 else "stable"}')
        print(f'   Rpot = {r:6.1f}:  ' + '  |  '.join(parts))

    # ---- single runs (both starts at once, one batch) ----
    rp = np.repeat(a.r, 2)
    v0 = np.tile([1.0, -1.0], len(a.r))
    print(f'\nintegrating {len(rp)} runs of {a.t} s at dt = {a.dt} ...')
    t, Y = integrate(rp, v0, a.t, a.dt, a.rl, g, keep=True, t_skip=a.skip)
    runs = [(r, t, Y[:, :, 2 * k], Y[:, :, 2 * k + 1])
            for k, r in enumerate(a.r)]

    print(f'{"Rpot":>7}  {"start":>6}  {"v1 min":>7}  {"v1 max":>7}  '
          f'{"v2 min":>7}  {"v2 max":>7}  {"n maxima":>8}  {"distinct":>8}')
    for r, _, yp, ym in runs:
        for tag, y in [('+1 V', yp), ('-1 V', ym)]:
            v1 = y[:, 0]
            pk = v1[1:-1][(v1[1:-1] > v1[:-2]) & (v1[1:-1] >= v1[2:])]
            distinct = len(np.unique(np.round(pk, 2)))
            print(f'{r:7.1f}  {tag:>6}  {v1.min():7.3f}  {v1.max():7.3f}  '
                  f'{y[:, 1].min():7.3f}  {y[:, 1].max():7.3f}  '
                  f'{len(pk):8d}  {distinct:8d}')

    records = find_records(a.records)
    plot_nr(g, a.r, os.path.join(a.out, f'simulated_nr{a.tag}.png'),
            a.shunt, i0)
    plot_portraits(runs, records,
                   os.path.join(a.out, f'simulated_portraits{a.tag}.png'),
                   a.n_meas)
    plot_timeseries(runs,
                    os.path.join(a.out, f'simulated_timeseries{a.tag}.png'))

    if a.export:
        d = os.path.join(a.out, f'simulated{a.tag}')
        os.makedirs(d, exist_ok=True)
        for r, t, yp, ym in runs:
            for tag, y in [('plus', yp), ('minus', ym)]:
                v3 = midpoint(y[:, 0], y[:, 1], r)
                f = os.path.join(d, f'sim {r:.1f} ohm {tag}.csv')
                np.savetxt(f, np.column_stack([t - t[0], y[:, 0], y[:, 1], v3]),
                           delimiter=',', fmt='%.6e',
                           header='Time(s),CH1(V),CH2(V),CH3(V)', comments='')
        print(f'\nwrote {2 * len(runs)} records to {d}')

    # ---- bifurcation sweep ----
    if not a.no_sweep:
        rmin, rmax, step = a.sweep
        rsw = np.arange(rmin, rmax + step / 2, step)
        print(f'\nsweeping {len(rsw)} values of Rpot x 2 starts, '
              f'{a.sweep_t} s each ...')
        mp = integrate(rsw, 1.0, a.sweep_t, a.dt, a.rl, g, keep=False,
                       t_skip=a.skip)
        mm = integrate(rsw, -1.0, a.sweep_t, a.dt, a.rl, g, keep=False,
                       t_skip=a.skip)
        measured = load_measured_bifurcation(a.records)
        plot_bifurcation(rsw, mp, mm, measured,
                         os.path.join(a.out, f'simulated_bifurcation{a.tag}.png'))
        print(f'\n{"Rpot":>7}  {"from +1 V":>24}  {"from -1 V":>24}')
        for r, a1, a2 in zip(rsw, mp, mm):
            if abs(r - round(r / 25) * 25) > step / 2:
                continue
            print(f'{r:7.1f}  {describe(a1):>24}  {describe(a2):>24}')

    print('\nfigures written to', a.out)


if __name__ == '__main__':
    main()
