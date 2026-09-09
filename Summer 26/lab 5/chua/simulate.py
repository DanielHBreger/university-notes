"""Simulate the measured Chua circuit with fixed-step RK4.

C1 dv1/dt = (v2-v1)/Rt - g(v1)
C2 dv2/dt = (v1-v2)/Rt - iL
L diL/dt = v2 - rL*iL; Rt = 990 ohm + Rpot.

The five measured segments are read from ../diode_fit.json, produced by
find_breakpoints.py using voltage ACROSS the nonlinear element. No shunt
correction is needed. The default retains the measured current offset.
--i0 auto is an explicit alternative calibration assumption, not a measured
correction: it sets g(0)=0. Its effect on the equilibria must be reported.
Independent fits are joined over --blend volts. Both starts (+1 and -1 V)
are simulated at every resistance; this tests initial-condition dependence,
not experimental up/down continuation hysteresis. Model/measurement disagreement
must remain visible, and a time-step convergence check is required.

Run: python simulate.py
     python simulate.py --i0 auto --tag _corrected
"""
import argparse
import glob
import os
import re
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scope_data import read_scope, R0

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- circuit values from page 1 of the experiment plan --------------------
C1 = 10e-9        # F
C2 = 100e-9       # F
L = 18e-3         # H
RS_BENCH = 216.0  # ohm, the shunt of the V-I measurement

# ---- Reproducible fitted data: element voltage in V, slope in S, current in A ----
def load_nr_segments(path=None):
    path = path or os.path.join(HERE, '..', 'diode_fit.json')
    with open(path, encoding='utf-8') as fh:
        fit = json.load(fh)
    if fit.get('voltage_basis') != 'element':
        raise ValueError('diode_fit.json must use voltage across the element')
    return [(f['v_lo'], f['v_hi'], f['slope_S'], f['intercept_A'])
            for f in fit['segments']]


NR_SEGMENTS = load_nr_segments()
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
        if not np.isfinite(k) or k <= 0:
            raise ValueError('shunt transformation is not single-valued')
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
    if blend <= 0 or not np.isfinite(blend):
        raise ValueError('blend must be positive and finite')
    vs, cs = [], []
    for k, (lo, hi, G, c) in enumerate(segments_v):
        a = -extend if k == 0 else 0.5 * (segments_v[k - 1][1] + lo) + blend / 2
        z = extend if k == n - 1 else 0.5 * (hi + segments_v[k + 1][0]) - blend / 2
        vs += [a, z]
        cs += [G * a + c, G * z + c]
    vs, cs = np.array(vs), np.array(cs)
    if np.any(np.diff(vs) <= 0):
        raise ValueError('N_R knots not monotonic; use a smaller positive --blend')
    return vs, cs


def make_g(blend=0.05, shunt=0.0, i0=0.0, segments=None):
    if shunt != 0 and segments is None:
        raise ValueError('the saved fit already uses element voltage; use --shunt 0')
    segs = nr_segments_v(NR_SEGMENTS if segments is None else segments, shunt, i0)
    vk, ik = nr_knots(segs, blend=blend)

    def g(v):
        v = np.asarray(v)
        value = np.interp(v, vk, ik)
        value = np.where(v < vk[0], segs[0][2]*v + segs[0][3], value)
        return np.where(v > vk[-1], segs[-1][2]*v + segs[-1][3], value)
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
    if Rt <= 0 or rL < 0:
        raise ValueError('Rt must be positive and rL nonnegative')
    vs = g.knots[0] if hasattr(g, 'knots') else np.linspace(-40, 40, 80001)
    hv = h(vs)
    idx = np.flatnonzero(hv[:-1] * hv[1:] < 0)
    roots = list(vs[np.abs(hv) < 1e-12])
    roots.extend(brentq(h, vs[k], vs[k + 1]) for k in idx)
    if hasattr(g, 'segments'):
        for side,k in [(-1,0),(1,-1)]:
            _,_,G,c = g.segments[k]
            slope = G + 1/(Rt+rL)
            if abs(slope)>1e-15:
                root = -c/slope
                if side*(root-vs[k]) > 0:
                    roots.append(root)
    roots = sorted(roots)
    return [r for i,r in enumerate(roots) if i==0 or abs(r-roots[i-1])>1e-7]


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
    if rpot.ndim != 1 or not np.isfinite(rpot).all() or np.any(rpot < 0):
        raise ValueError('potentiometer resistances must be finite and nonnegative')
    if not np.isfinite([t_end, dt, rL, t_skip]).all() or dt <= 0 or rL < 0 or not 0 <= t_skip < t_end:
        raise ValueError('require dt > 0, rL >= 0 and 0 <= skip < duration')
    v1_0 = np.broadcast_to(np.asarray(v1_0, float), rpot.shape)
    Rt = R0 + rpot
    n = len(rpot)
    y = np.zeros((3, n))
    y[0] = v1_0
    if not np.isfinite(v1_0).all():
        raise ValueError('initial voltage must be finite')
    nsteps = int(np.floor(t_end / dt + 1e-10))
    n_skip = int(np.ceil(t_skip / dt - 1e-10))
    if nsteps <= n_skip:
        raise ValueError('no integration steps remain after the transient')

    if keep:
        out = np.empty((nsteps - n_skip + 1, 3, n))
        out[0] = y
        j = 1
        for k in range(1, nsteps + 1):
            y = rk4_step(y, dt, Rt, rL, g)
            if k % 1000 == 0 and (not np.isfinite(y).all() or np.max(np.abs(y))>1e6):
                raise ValueError('integration diverged; reduce --dt or check the model')
            if k == n_skip:
                out[0] = y
                j = 1
            elif k > n_skip:
                out[j] = y
                j += 1
        if not np.isfinite(out).all():
            raise ValueError('integration produced non-finite states')
        t = dt * np.arange(n_skip, nsteps + 1)
        return t, out

    maxima = [[] for _ in range(n)]
    prev2 = y[0].copy()
    prev1 = y[0].copy()
    for k in range(1, nsteps + 1):
        y = rk4_step(y, dt, Rt, rL, g)
        if k % 1000 == 0 and (not np.isfinite(y).all() or np.max(np.abs(y))>1e6):
            raise ValueError('integration diverged; reduce --dt or check the model')
        cur = y[0]
        if k > n_skip:
            hit = (prev1 - prev2 > 1e-10) & (prev1 >= cur)
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
    t, channels = read_scope(path)
    if n_max is not None:
        if n_max < 3:
            raise ValueError('n_max must allow at least three samples')
        t, channels = t[:n_max], channels[:n_max]
    return t, channels[:, 0], channels[:, 1]


def load_measured_bifurcation(folder):
    """Read measurement sidecars once per record, excluding simulation exports."""
    import csv
    records = {}
    for path in sorted(glob.glob(os.path.join(folder, '*_bifurcation_points.csv'))):
        with open(path, newline='') as fh:
            reader = csv.DictReader(fh)
            if not {'sweep','filename','rpot_ohm','max_v'}.issubset(reader.fieldnames or []):
                continue
            current = {}
            for row in reader:
                try:
                    value = (float(row['rpot_ohm']), float(row['max_v']), row['sweep'])
                except (ValueError, TypeError):
                    continue
                if not np.isfinite(value[:2]).all():
                    continue
                current.setdefault((row['sweep'],row['filename']),[]).append(value)
            for key,values in current.items():
                records.setdefault(key,values)
    return [value for values in records.values() for value in values]


# ---- figures ---------------------------------------------------------------
def plot_nr(g, rpots, path, shunt=0.0, i0=0.0, rL=0.0):
    v = np.linspace(-10, 10, 2001)
    fig, ax = plt.subplots(figsize=(8, 5))
    raw = os.path.join(HERE, '..', 'trace1.csv')
    if os.path.exists(raw):
        d = np.genfromtxt(raw, delimiter=',', skip_header=1)
        i_m = (d[:, 1] - d[:, 2]) / RS_BENCH
        ax.plot(d[:, 2], (i_m - i0) * 1e3, '.', ms=1,
                alpha=0.15, color='steelblue',
                label='measured (trace1.csv)'
                      + (', same corrections' if shunt or i0 else ''))
    ax.plot(v, g(v) * 1e3, 'crimson', lw=2, label='g(v): fitted segments')
    for (_, hi, _, _), (lo, _, _, _) in zip(g.segments[:-1], g.segments[1:]):
        ax.axvline(0.5 * (hi + lo), color='grey', ls='--', lw=0.7)
    for r in rpots:
        ax.plot(v, -v / (R0 + r + rL) * 1e3, lw=0.9,
                label=f'load line, Rpot = {r:g} ohm')
    ax.set_xlabel('v across N_R  (V)')
    ax.set_ylabel('current into N_R  (mA)')
    what = (f'{i0 * 1e3:.3f} mA offset removed (assumption)'
            if i0 else 'Measured current offset retained')
    ax.set_title(f'Nonlinear element used in the simulation\n{what}')
    ax.set_ylim(-6, 6)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_portraits(runs, records, path, n_meas, model_note='Measured current offset retained'):
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
            if np.ptp(y[:,0]) < 1e-7:
                ax.plot(y[-1,0],y[-1,1],'o',color='C0',ms=4)
                ax.text(.03,.92,'Settled equilibrium',transform=ax.transAxes,fontsize=8)
            else:
                ax.plot(y[:, 0], y[:, 1], lw=0.35, alpha=.7, color='C0')
            ax.set_title(f'simulated, Rpot = {r:g} ohm, {tag}', fontsize=9)
        ax = axes[row, 2]
        rec = records.get(r)
        if rec:
            from lorenz_map import period_samples
            from scipy.signal import savgol_filter
            from rpot import rpot as measure_rpot
            _, v1, v2 = load_record(rec, n_meas)
            window = max(5, (period_samples(v1)//20)|1)
            window = min(window, len(v1) if len(v1)%2 else len(v1)-1)
            step = max(1, len(v1)//60000)
            ax.scatter(v1[::step], v2[::step], s=2, alpha=.08, color='0.5', lw=0,
                       rasterized=True)
            if window > 3:
                a,b = savgol_filter(v1,window,3),savgol_filter(v2,window,3)
                ax.scatter(a[::step],b[::step],s=1,alpha=.15,color='C3',lw=0,
                           rasterized=True)
            measured_r = measure_rpot(rec,R0,1.0)[0]
            ax.set_title(f'Measured Rpot = {measured_r:.1f} Ω (nominal {r:g} Ω)\n'
                         f'Grey: raw codes; red: SG {window} samples', fontsize=9)
        else:
            ax.set_title('no measured record at this Rpot', fontsize=9)
            ax.set_axis_off()
        for ax in axes[row]:
            ax.grid(alpha=0.3)
            ax.set_xlabel('v1 = v_C1 (V)')
            ax.set_ylabel('v2 = v_C2 (V)')
        available = [yp[:,[0,1]], ym[:,[0,1]]]
        if rec:
            available.append(np.column_stack([v1,v2]))
        points = np.concatenate(available)
        for dim in (0,1):
            lo,hi = points[:,dim].min(),points[:,dim].max()
            pad = max(.08*(hi-lo),.05)
            for ax in axes[row]:
                (ax.set_xlim if dim==0 else ax.set_ylim)(lo-pad,hi+pad)
    fig.suptitle(f'Chua phase portraits · {model_note}\nShared voltage limits within each row')
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_timeseries(runs, path, window=0.01, model_note='Measured current offset retained'):
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
    fig.suptitle(f'Simulated V1 · {model_note}')
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_bifurcation(rpot, max_plus, max_minus, measured, path,
                     model_note='Measured current offset retained'):
    fig, ax = plt.subplots(figsize=(14, 7.5))
    if measured:
        sweeps = sorted({s for _, _, s in measured})
        for k,s in enumerate(sweeps):
            c = plt.get_cmap('Greys')(.45+.4*k/max(len(sweeps)-1,1))
            pts = np.array([(r, m) for r, m, sw in measured if sw == s])
            ax.plot(pts[:, 0], pts[:, 1], '.', ms=0.7, color=c, alpha=0.6,
                    label=f'measured, {s}', zorder=1)
    for maxima, c, tag in [(max_plus, 'C0', 'sim, v1(0) = +1 V'),
                           (max_minus, 'C3', 'sim, v1(0) = -1 V')]:
        xs = np.concatenate([np.full(len(m), r) for r, m in zip(rpot, maxima)])
        ys = np.concatenate([np.asarray(m, float) for m in maxima])
        ax.plot(xs, ys, '.', ms=1.0, color=c, alpha=0.5, label=tag, zorder=2)
    ax.set_xlabel('Rpot (ohm)')
    ax.set_ylabel('local maxima of v1 = v_C1 (V)')
    ax.set_title(f'Simulation and measurement · {model_note}\n'
                 'Two initial conditions at each resistance; this is not a continuation sweep')
    ax.grid(alpha=0.3)
    ax.legend(markerscale=4, fontsize=9)
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
                        'to the voltage across N_R; saved fits require 0, '
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
    os.makedirs(a.out, exist_ok=True)
    if a.n_meas is not None and a.n_meas < 3:
        p.error('--n-meas must be at least 3')

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
    model_note = 'Offset removed (assumption)' if i0 else 'Measured current offset retained'
    plot_nr(g, a.r, os.path.join(a.out, f'simulated_nr{a.tag}.png'),
            a.shunt, i0, a.rl)
    plot_portraits(runs, records,
                   os.path.join(a.out, f'simulated_portraits{a.tag}.png'),
                   a.n_meas, model_note)
    plot_timeseries(runs,
                    os.path.join(a.out, f'simulated_timeseries{a.tag}.png'), model_note=model_note)
    np.savez_compressed(os.path.join(a.out, f'simulated_runs{a.tag}.npz'),
                        t=t, states=Y, rpot=rp, initial_v1=v0,
                        dt=a.dt, rL=a.rl, removed_offset_A=i0)

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
        if not np.isfinite([rmin,rmax,step]).all() or step <= 0 or rmin < 0 or rmax < rmin:
            p.error('sweep requires 0 <= RMIN <= RMAX and STEP > 0')
        rsw = np.arange(rmin, rmax + step / 2, step)
        print(f'\nsweeping {len(rsw)} values of Rpot x 2 starts, '
              f'{a.sweep_t} s each ...')
        mp = integrate(rsw, 1.0, a.sweep_t, a.dt, a.rl, g, keep=False,
                       t_skip=a.skip)
        mm = integrate(rsw, -1.0, a.sweep_t, a.dt, a.rl, g, keep=False,
                       t_skip=a.skip)
        measured = load_measured_bifurcation(a.records)
        plot_bifurcation(rsw, mp, mm, measured,
                         os.path.join(a.out, f'simulated_bifurcation{a.tag}.png'), model_note)
        import csv
        with open(os.path.join(a.out,f'simulated_bifurcation{a.tag}_points.csv'),'w',newline='') as fh:
            writer=csv.writer(fh); writer.writerow(['rpot_ohm','initial_v1_V','max_v'])
            for start,maps in [(1,mp),(-1,mm)]:
                for r,values in zip(rsw,maps):
                    writer.writerows((r,start,float(value)) for value in values)
        print(f'\n{"Rpot":>7}  {"from +1 V":>24}  {"from -1 V":>24}')
        for r, a1, a2 in zip(rsw, mp, mm):
            if abs(r - round(r / 25) * 25) > step / 2:
                continue
            print(f'{r:7.1f}  {describe(a1):>24}  {describe(a2):>24}')

    print('\nfigures written to', a.out)


if __name__ == '__main__':
    main()
