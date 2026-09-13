"""
Identify the circuit that was actually on the bench from the oscillator records.

Every quantity here comes from the sweep records themselves (V1, V2 and the
divider-fitted Rpot), not from the nominal component values. The periodic
records (the period-1 cycle below the Hopf point, the large outer cycle and the
small orbit near the origin) are averaged over their hundreds of cycles, which
removes the scope's quantisation and leaves waveforms clean enough to
differentiate. The circuit equations then give each element directly:

  node 1   C1 v1' = (v2 - v1)/Rt - g(v1)
           For a static element the loop integral of g dv1 over a cycle is
           zero, so   C1 = oint (v2 - v1)/Rt v1' dt / oint v1'^2 dt   for ANY g.
           i_NR = (v2 - v1)/Rt - C1 v1' against v1 is then the element's curve
           as the circuit sees it, and any loop in it is a dynamic effect.
  node 2   (v1 - v2)/Rt = C2 v2' + iL,   L iL' = v2 - r iL
           C2 comes from the admittance (v1 - v2)/Rt / v2 harmonic by harmonic.
           With C2 known, iL = (v1 - v2)/Rt - C2 v2' and the flux
           phi = int v2 dt against iL is the inductor's own loop: its secant
           slope is the effective inductance, its area the loss per cycle.
  diode    the five-segment law is fitted to ALL records (the double scroll
           covers -5..+5 V continuously, the large cycle reaches the
           saturation segments) by integrated KCL with a spline basis, and
           the eight numbers (four breakpoints, three distinct slopes) are
           mapped onto Kennedy's two-op-amp realisation, which they fit.

Outputs, beside this script:
    identified.json        every number below, used by simulate.py --model bench
    identified.txt         the same as a readable report
    identified_inductor.png   flux-current loops and L_eff, r_eff against amplitude
    identified_diode.png      i_NR(v1) from the cycles, the fitted PWL and the V-I trace
    identified_c1.png         loop-integral C1 against the record's amplitude

Usage:
    python identify.py               # forward/ and back/ beside this script
    python identify.py --quick       # every 4th record (about 2 min)
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import BSpline
from scipy.optimize import least_squares

from scope_data import R0, read_scope
from lorenz_map import maxima_from_arrays
from integration import pwl_knots

HERE = os.path.dirname(os.path.abspath(__file__))
NPH = 2048          # phase points of an averaged cycle
KMAX = 80           # harmonics kept when differentiating the averaged cycle
RB_NOMINAL = 22e3   # ohm, the 22 k feedback resistor of Kennedy's second stage


# ---- cycle averaging -------------------------------------------------------
def averaged_cycle(t, v, t_peaks, nph=NPH, tol=0.05):
    """
    Mean of the columns of v over the cycles between successive maxima, on
    `nph` phase points. Cycles whose length differs from the median by more
    than `tol` are left out. Returns (cycle, mean period, cycles used).
    """
    T = np.diff(t_peaks)
    Tm = np.median(T)
    keep = np.flatnonzero(np.abs(T - Tm) < tol * Tm)
    if len(keep) < 20:
        return None
    ph = np.arange(nph) / nph
    acc = np.zeros((nph, v.shape[1]))
    for i in keep:
        tt = t_peaks[i] + ph * T[i]
        for j in range(v.shape[1]):
            acc[:, j] += np.interp(tt, t, v[:, j])
    return acc / len(keep), float(T[keep].mean()), len(keep)


def harmonics(x):
    return np.fft.rfft(x) / len(x)


def derivative(x, T, kmax=KMAX):
    X = np.fft.rfft(x)
    k = np.arange(len(X))
    X = X * (1j * 2 * np.pi * k / T)
    X[k > kmax] = 0
    return np.fft.irfft(X, len(x))


def lowpass(x, kmax=KMAX):
    X = np.fft.rfft(x)
    X[kmax + 1:] = 0
    return np.fft.irfft(X, len(x))


def loop_c1(v1, v2, dv1, Rt):
    """C1 from the loop integral: exact for any static element g(v1)."""
    ir = (v2 - v1) / Rt
    return float(np.sum(ir * dv1) / np.sum(dv1 * dv1))


def inductor_loop(v1, v2, dv2, Rt, C2, T):
    """
    (iL, flux, I_amp, L_eff, r_eff): the inductor current from node 2, the
    flux int v2 dt, the current amplitude, the secant inductance of the
    flux loop and the loss resistance from the in-phase power.
    """
    iL = (v1 - v2) / Rt - C2 * dv2
    iL = iL - iL.mean()             # the channel offsets, not the tiny DC drop
    v2c = v2 - v2.mean()
    dt = T / len(v1)
    phi = np.cumsum(v2c) * dt
    phi -= phi.mean()
    I = 0.5 * np.ptp(iL)
    L_eff = np.ptp(phi) / np.ptp(iL)
    r_eff = float(np.sum(v2c * iL) / np.sum(iL * iL))
    return iL, phi, I, float(L_eff), r_eff


def admittance_rows(v1, v2, Rt, T, kmax=24, vmin=2e-4):
    """(omega, Y) per harmonic of the tank admittance (v1 - v2)/Rt / v2."""
    IR = harmonics((v1 - v2) / Rt)
    V2 = harmonics(v2 - v2.mean())
    w = 2 * np.pi * np.arange(len(V2)) / T
    ks = np.arange(1, kmax + 1)
    good = ks[np.abs(V2[ks]) > vmin]
    return [(float(w[k]), complex(IR[k] / V2[k]), float(abs(V2[k]))) for k in good]


def fit_tank(rows):
    """Least squares of Y = j w C2 + 1/(r + j w L) over (w, Y, weight) rows."""
    w = np.array([r[0] for r in rows]); Y = np.array([r[1] for r in rows])
    wt = np.array([r[2] for r in rows])

    def fun(p):
        C2, L, r = p[0] * 1e-9, p[1] * 1e-3, p[2]
        res = (1j * w * C2 + 1 / (r + 1j * w * L) - Y) * wt
        return np.r_[res.real, res.imag]
    f = least_squares(fun, [100, 18, 5], bounds=([10, 1, 0], [500, 100, 500]),
                      xtol=1e-12, ftol=1e-12, gtol=1e-12)
    rel = np.linalg.norm(f.fun) / np.linalg.norm(np.r_[(Y * wt).real, (Y * wt).imag])
    return f.x, float(rel)


# ---- the diode from integrated KCL -----------------------------------------
VLO, VHI, DV = -8.0, 7.4, 0.2
KNOTS = np.r_[[VLO] * 4, np.arange(VLO + DV, VHI - DV / 2, DV), [VHI] * 4]
NB = len(KNOTS) - 4


def spline_design(v):
    return BSpline.design_matrix(np.clip(v, VLO, VHI), KNOTS, 3, extrapolate=False).toarray()


def kcl_normal_equations(path, rpot, window_s=30e-6, stride=25, nrows=400000):
    """
    Accumulate X'X, X'y of   int_w (v2 - v1)/Rt dt = C1 dv1|_w + int_w B(v1) dt c
    over windows w of one record, with a spline basis B for g(v1).
    """
    t, d = read_scope(path, ('CH1', 'CH2'))
    t, d = t[:nrows], d[:nrows]
    dt = float(np.median(np.diff(t)))
    n = max(2, int(round(window_s / dt)))
    v1, v2 = d[:, 0], d[:, 1]
    IB = cumulative_trapezoid(spline_design(v1), dx=dt, axis=0, initial=0)
    ir = cumulative_trapezoid((v2 - v1) / (R0 + rpot), dx=dt, initial=0)
    ix = np.arange(0, len(t) - n, stride)
    X = np.column_stack([v1[ix + n] - v1[ix], IB[ix + n] - IB[ix]])
    y = ir[ix + n] - ir[ix]
    X, y = X / np.sqrt(len(ix)), y / np.sqrt(len(ix))
    return X.T @ X, X.T @ y, float(y @ y)


def solve_kcl(XtX, Xty, yty, lam=1e-9):
    D = np.diff(np.eye(NB), n=2, axis=0)
    P = np.zeros_like(XtX)
    P[1:, 1:] = D.T @ D * lam
    coef = np.linalg.solve(XtX + P, Xty)
    resid2 = yty - 2 * coef @ Xty + coef @ XtX @ coef
    return coef[0], coef[1:], float(np.sqrt(max(resid2, 0) / yty))


def pwl_from_spline(cs, inner=True):
    """Slopes and breakpoints of the five-segment law read off the spline."""
    v = np.linspace(VLO, VHI, 3081)
    s = BSpline(KNOTS, cs, 3)
    g, dg = s(v), s.derivative()(v)
    out = {}
    for lo, hi, name in [(-6.2, -1.6, 'Gb_left'), (-0.5, 0.4, 'Ga'), (1.4, 5.2, 'Gb_right'),
                         (-7.6, -7.0, 'Gc_left'), (6.3, 7.0, 'Gc_right')]:
        m = (v > lo) & (v < hi)
        out[name] = float(np.polyfit(v[m], g[m], 1)[0])
    mid_in = 0.5 * (out['Ga'] + 0.5 * (out['Gb_left'] + out['Gb_right']))
    for side, name in [(-1, 'bp_in_left'), (1, 'bp_in_right')]:
        m = (v * side > 0.2) & (v * side < 3)
        idx = np.flatnonzero(np.diff(np.sign(dg[m] - mid_in)))
        out[name] = float(v[m][idx[0]]) if len(idx) else np.nan
    for side, name in [(-1, 'bp_out_left'), (1, 'bp_out_right')]:
        m = (v * side > 4) & (v * side < 7.4)
        idx = np.flatnonzero(np.diff(np.sign(dg[m] - 1.5e-3)))
        out[name] = float(v[m][idx[0]]) if len(idx) else np.nan
    out['g0_A'] = float(s(0.0))
    return out, v, g


def saturation_from_plateau(v1, v2, dv1, Rt, C1, side, frac=0.03):
    """
    (slope, intercept) of the saturation segment from the plateau of a large
    cycle: while v1 sits on the steep segment v1' is tiny, so
    i_NR = (v2 - v1)/Rt - C1 v1' is known without any model of the element,
    and it moves along the segment as v2 swings.
    """
    i = (v2 - v1) / Rt - C1 * dv1
    m = (np.abs(dv1) < frac * np.abs(dv1).max()) & (v1 * side > 4.0)
    if m.sum() < 20 or np.ptp(v1[m]) < 0.1:
        return None
    G, c = np.polyfit(v1[m], i[m], 1)
    return float(G), float(c)


def kennedy_from_pwl(Ga, Gb, Gc, bp_in_left, bp_in_right, bp_out_left, bp_out_right, RB=RB_NOMINAL):
    """
    Kennedy's diode: stage A (gain AA, resistor RA) saturates at the outer
    breakpoints, stage B (gain AB, resistor RB) at the inner ones, each at
    its own output rails:
        Ga = (1 - AA)/RA + (1 - AB)/RB,  Gb = (1 - AA)/RA + 1/RB,  Gc = 1/RA + 1/RB
        inner breakpoints VnB/AB, VpB/AB;  outer VnA/AA, VpA/AA.
    With RB fixed at its nominal 22 k the eight segment numbers give the
    six parameters exactly.
    """
    RA = 1.0 / (Gc - 1.0 / RB)
    AA = 1.0 - (Gb - 1.0 / RB) * RA
    AB = (Gb - Ga) * RB
    return dict(RA_ohm=float(RA), AA=float(AA), RB_ohm=float(RB), AB=float(AB),
                VpA_V=float(AA * bp_out_right), VnA_V=float(AA * bp_out_left),
                VpB_V=float(AB * bp_in_right), VnB_V=float(AB * bp_in_left))


def rayleigh_fit(I, L_eff, r_eff, bin_mA=1.0):
    """
    L_eff = L0 + nu I / 2  and  r_eff = r0 + (8/3pi) rho I: the secant
    inductance and the loss resistance of L(i) = L0 + nu|i|, r(i) = r0 + rho|i|
    for a sinusoid of amplitude I. Every amplitude bin of `bin_mA` carries the
    same total weight, so the forty-odd large-cycle records do not outvote the
    small cycles that fix the small-signal end of the law.
    """
    I, L_eff, r_eff = np.asarray(I, float), np.asarray(L_eff, float), np.asarray(r_eff, float)
    bins = np.floor(I * 1e3 / bin_mA).astype(int)
    counts = {b: int(np.sum(bins == b)) for b in np.unique(bins)}
    w = np.array([1.0 / counts[b] for b in bins])
    a = np.polyfit(I, L_eff, 1, w=np.sqrt(w))
    b = np.polyfit(I, r_eff, 1, w=np.sqrt(w))
    return dict(L0_H=float(a[1]), nu_H_per_A=float(2 * a[0]),
                r0_ohm=float(b[1]), rho_ohm_per_A=float(b[0] / (8 / (3 * np.pi))))


# ---- records ---------------------------------------------------------------
def small_signal_period(g, Rt, c1, c2, l0, r0):
    """Period of the linearisation about the negative equilibrium of the circuit with a static element g."""
    from scipy.optimize import brentq
    h = lambda v: g(v) + v / (Rt + r0)
    vs = np.linspace(-15, -1e-3, 20001)
    hv = h(vs)
    idx = np.flatnonzero(hv[:-1] * hv[1:] < 0)
    if not len(idx):
        return np.nan
    v = brentq(h, vs[idx[-1]], vs[idx[-1] + 1])
    dv = 1e-4
    G = (g(v + dv) - g(v - dv)) / (2 * dv)
    J = np.array([[(-1 / Rt - G) / c1, 1 / (Rt * c1), 0.0],
                  [1 / (Rt * c2), -1 / (Rt * c2), -1 / c2],
                  [0.0, 1 / l0, -r0 / l0]])
    ev = np.linalg.eigvals(J)
    ev = ev[np.abs(ev.imag) > 1.0]
    return 2 * np.pi / np.abs(ev.imag).max() if len(ev) else np.nan


def small_signal_hopf(g, c1, c2, l0, r0, rt_range=(1200.0, 3200.0)):
    """Rt at which the negative equilibrium's complex pair crosses the axis (None if never)."""
    from scipy.optimize import brentq

    def re_part(Rt):
        h = lambda v: g(v) + v / (Rt + r0)
        vs = np.linspace(-15, -1e-3, 20001)
        hv = h(vs)
        idx = np.flatnonzero(hv[:-1] * hv[1:] < 0)
        if not len(idx):
            return np.nan
        v = brentq(h, vs[idx[-1]], vs[idx[-1] + 1])
        dv = 1e-4
        G = (g(v + dv) - g(v - dv)) / (2 * dv)
        J = np.array([[(-1 / Rt - G) / c1, 1 / (Rt * c1), 0.0],
                      [1 / (Rt * c2), -1 / (Rt * c2), -1 / c2],
                      [0.0, 1 / l0, -r0 / l0]])
        ev = np.linalg.eigvals(J)
        ev = ev[np.abs(ev.imag) > 1.0]
        return ev.real.max() if len(ev) else np.nan
    Rs = np.arange(rt_range[0], rt_range[1], 10.0)
    vals = np.array([re_part(R) for R in Rs])
    idx = [k for k in np.flatnonzero(np.sign(vals[:-1]) != np.sign(vals[1:]))
           if np.isfinite(vals[k]) and np.isfinite(vals[k + 1])]
    return brentq(re_part, Rs[idx[-1]], Rs[idx[-1] + 1]) if idx else None


def sweep_records(sweep):
    tab = pd.read_csv(os.path.join(HERE, f'{sweep}_rpot.csv'))
    tab = tab[tab.status.eq('ok')]
    return [(os.path.join(HERE, sweep, r.filename), float(r.rpot_ohm), f'{sweep}/{r.filename}')
            for r in tab.itertuples()]


def cycle_markers(t, v1, T_est):
    """
    Times at which v1 crosses the middle of its range upwards, at least
    0.6 T_est apart. The crossing is where v1 moves fastest, so it marks a
    cycle far more sharply than a maximum on a flat top (the large cycle).
    """
    mid = 0.5 * (v1.max() + v1.min())
    x = v1 - mid
    idx = np.flatnonzero((x[:-1] < 0) & (x[1:] >= 0))
    if not len(idx):
        return np.empty(0)
    tc = t[idx] + (t[idx + 1] - t[idx]) * (-x[idx]) / (x[idx + 1] - x[idx])
    keep = [tc[0]]
    for c in tc[1:]:
        if c - keep[-1] > 0.6 * T_est:
            keep.append(c)
    return np.array(keep)


def analyse_record(path, rpot):
    t, d = read_scope(path, ('CH1', 'CH2', 'CH3'))
    M, info = maxima_from_arrays(t, d[:, 0])
    if len(M) < 60:
        return None
    T_est = np.mean(np.diff(info['t_peaks']))
    tp = cycle_markers(t, d[:, 0], T_est)
    if len(tp) < 60:
        return None
    Tall = np.diff(tp)
    out = dict(rpot=rpot, Rt=R0 + rpot, n_max=len(M), T_us=Tall.mean() * 1e6,
               jitter=float(Tall.std() / Tall.mean()), max_spread=float(np.ptp(M)))
    # the divider identity A + B = 1 tests the relative gain of CH1 and CH2
    X = np.column_stack([d[:, 0] - d[:, 0].mean(), d[:, 1] - d[:, 1].mean()])
    A, B = np.linalg.lstsq(X, d[:, 2] - d[:, 2].mean(), rcond=None)[0]
    out['divider_A_plus_B'] = float(A + B)
    out['periodic'] = out['jitter'] < 0.015 and out['max_spread'] < 0.12
    if not out['periodic']:
        return out
    res = averaged_cycle(t, d[:, :2], tp)
    if res is None:
        out['periodic'] = False
        return out
    cyc, T, ncyc = res
    v1, v2 = lowpass(cyc[:, 0]), lowpass(cyc[:, 1])
    dv1, dv2 = derivative(cyc[:, 0], T), derivative(cyc[:, 1], T)
    out.update(T_cyc_us=T * 1e6, n_cycles=ncyc, v1_min=float(v1.min()), v1_max=float(v1.max()),
               v2_amp=float(0.5 * np.ptp(v2)), C1_loop=loop_c1(v1, v2, dv1, R0 + rpot))
    out['cycle'] = (v1, v2, dv1, dv2, T)
    out['harm'] = admittance_rows(v1, v2, R0 + rpot, T)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--quick', action='store_true', help='every 4th record')
    p.add_argument('--out', default=HERE)
    a = p.parse_args()
    step = 4 if a.quick else 1

    # ---- 1. cycle-averaged periodic records ----
    rows = []
    for sweep in ['forward', 'back']:
        recs = sweep_records(sweep)
        for k, (path, rpot, label) in enumerate(recs):
            n = int(os.path.basename(path)[5:-4])
            if k % step and not (sweep == 'forward' and n >= 402):
                continue
            r = analyse_record(path, rpot)
            if r is None:
                continue
            r['label'] = label
            rows.append(r)
            print(f'{label:22s} R={rpot:7.2f}  T={r["T_us"]:6.1f} us  '
                  + (f'C1 {r["C1_loop"]*1e9:6.3f} nF  v1 {r["v1_min"]:6.2f}..{r["v1_max"]:6.2f}'
                     if r['periodic'] else 'not periodic'), flush=True)
    per = [r for r in rows if r['periodic']]
    small = [r for r in per if r['v1_max'] < 5.5 and r['v1_min'] < -1.5]      # period-1 cycles
    tiny = [r for r in per if r['v1_max'] < 5.5 and r['v1_min'] >= -1.5]     # the orbit near the origin
    large = [r for r in per if r['v1_max'] > 5.5]                              # the large outer cycle
    gain = np.array([r['divider_A_plus_B'] for r in rows])

    # ---- 2. C2 from the admittance fits (>= 5 harmonics), then closed on the onset period ----
    hr = [h for r in small if len(r['harm']) >= 5 for h in r['harm']]
    (c2_adm_nF, L_fit_mH, r_fit), tank_resid = fit_tank(hr)
    C2 = c2_adm_nF * 1e-9

    # C1 at small amplitude (shoulders only), for the linearisation below
    onset = [r for r in small if r['v1_max'] < -2.0]
    c1_small = float(np.mean([r['C1_loop'] for r in onset])) if onset else np.nan      # F
    c1_p1 = float(np.mean([r['C1_loop'] for r in small if r['v1_max'] > -1.0]))         # F
    c1_lc = float(np.mean([r['C1_loop'] for r in large]))                                 # F

    # ---- 3. the diode: integrated KCL, inner law from the core records, saturation from the large cycle ----
    core, lc = [], []
    for sweep in ['forward', 'back']:
        for k, (path, rpot, label) in enumerate(sweep_records(sweep)):
            n = int(os.path.basename(path)[5:-4])
            if sweep == 'forward' and 330 < rpot < 900 and n % (2 * step) == 0:
                core.append((path, rpot))
            if ((sweep == 'forward' and rpot < 315) or (sweep == 'back' and rpot < 680)) and n % step == 0:
                lc.append((path, rpot))
    fits = {}
    for name, sel in [('core', core), ('large_cycle', lc)]:
        XtX = np.zeros((NB + 1, NB + 1)); Xty = np.zeros(NB + 1); yty = 0.0
        for path, rpot in sel:
            A, b, yy = kcl_normal_equations(path, rpot)
            XtX += A; Xty += b; yty += yy
        c1, cs, resid = solve_kcl(XtX, Xty, yty)
        pw, vgrid, ggrid = pwl_from_spline(cs)
        fits[name] = dict(records=len(sel), C1_apparent_nF=c1 * 1e9, relative_residual=resid, pwl=pw,
                          curve_v=vgrid, curve_i=ggrid)
        print(f'KCL fit {name}: {len(sel)} records, apparent C1 {c1*1e9:.3f} nF, residual {resid:.4f}', flush=True)
    inner = fits['core']['pwl']
    sat = {-1: [], 1: []}
    for r in large:
        v1, v2, dv1, dv2, T = r['cycle']
        for side in (-1, 1):
            f = saturation_from_plateau(v1, v2, dv1, r['Rt'], fits['core']['C1_apparent_nF'] * 1e-9, side)
            if f:
                sat[side].append(f)
    if not sat[-1] or not sat[1]:
        raise SystemExit('no large-cycle plateau to fit the saturation segments from')
    Gc_left, c_left = np.median(sat[-1], axis=0)
    Gc_right, c_right = np.median(sat[1], axis=0)
    # the shoulder lines through the inner breakpoints, then their intersections with the saturation lines
    Ga, GbL, GbR = inner['Ga'], inner['Gb_left'], inner['Gb_right']
    bpL, bpR = inner['bp_in_left'], inner['bp_in_right']
    cL = Ga * bpL - GbL * bpL                    # i = GbL v + cL on the left shoulder
    cR = Ga * bpR - GbR * bpR
    bp_out_left = (c_left - cL) / (GbL - Gc_left)
    bp_out_right = (c_right - cR) / (GbR - Gc_right)
    Gb = 0.5 * (GbL + GbR)
    Gc = 0.5 * (Gc_left + Gc_right)
    pwl = dict(bp_V=[float(bp_out_left), bpL, bpR, float(bp_out_right)],
               G_S=[float(Gc_left), GbL, Ga, GbR, float(Gc_right)],
               saturation_records=[len(sat[-1]), len(sat[1])])
    kennedy = kennedy_from_pwl(Ga, Gb, Gc, bpL, bpR, bp_out_left, bp_out_right)

    # ---- 4. C2 from the near-origin orbit, then the inductor loops and the Rayleigh law ----
    # The period-1 admittance fit pins L*C2 but splits it poorly between C2 and L:
    # at 3 kHz and 2 mA the inductor is already nonlinear over the harmonics the
    # fit uses. The small orbit around the origin (about 1.4 kHz, 0.7 mA, ten
    # harmonics) is the record where the inductor is most nearly linear and the
    # harmonics reach highest, so C2 is taken from it when it exists.
    hr_tiny = [h for r in tiny if len(r['harm']) >= 6 for h in r['harm']]
    if hr_tiny:
        (c2_nF, L_tiny_mH, r_tiny), tiny_resid = fit_tank(hr_tiny)
        c2_source = f'near-origin orbit ({len(tiny)} records, {len(hr_tiny)} harmonics)'
    else:
        c2_nF, L_tiny_mH, r_tiny, tiny_resid = c2_adm_nF, L_fit_mH, r_fit, tank_resid
        c2_source = 'period-1 admittance fit (no near-origin record)'
    C2 = c2_nF * 1e-9
    tank = []
    for r in per:
        v1, v2, dv1, dv2, T = r['cycle']
        iL, phi, I, L_eff, r_eff = inductor_loop(v1, v2, dv2, r['Rt'], C2, T)
        r['iL'], r['phi'] = iL, phi
        r.update(I_amp=I, L_eff=L_eff, r_eff=r_eff)
        tank.append((r['label'], r['rpot'], r['T_cyc_us'], I, L_eff, r_eff))
    I = np.array([t[3] for t in tank]); Le = np.array([t[4] for t in tank]); re = np.array([t[5] for t in tank])
    ok = I > 0.6e-3            # below this the flux loop is a few quantisation steps wide
    rayleigh = rayleigh_fit(I[ok], Le[ok], re[ok])
    # The small-signal end of the inductance law comes from the period of the
    # oscillation at onset, which is measured to 0.3 %: the flux loops of the
    # smallest cycles are only a few quantisation steps wide, so their
    # extrapolation to zero amplitude is the least certain number here. L0 is
    # the inductance at which the linearised circuit (C1 small-signal, C2, r0,
    # the fitted element) oscillates at the onset period; nu is then refitted
    # to the loops with that L0.
    from scipy.optimize import brentq
    onset_T = [r for r in small if r['v1_max'] - r['v1_min'] < 1.5]
    T_onset = float(np.mean([r['T_cyc_us'] for r in onset_T])) * 1e-6
    R_onset = float(np.mean([r['rpot'] for r in onset_T]))
    vk_s, ik_s = pwl_knots(pwl['bp_V'], pwl['G_S'])
    g_static = lambda v: np.interp(v, vk_s, ik_s)
    rayleigh['L0_loops_H'] = rayleigh['L0_H']
    rayleigh['nu_loops_H_per_A'] = rayleigh['nu_H_per_A']
    f = lambda l0: small_signal_period(g_static, R0 + R_onset, c1_small, C2, l0, rayleigh['r0_ohm']) - T_onset
    if np.isfinite(f(10e-3)) and np.isfinite(f(40e-3)) and f(10e-3) * f(40e-3) < 0:
        rayleigh['L0_H'] = float(brentq(f, 10e-3, 40e-3, xtol=1e-7))
        bins = np.floor(I[ok] * 1e3).astype(int)
        counts = {b: int(np.sum(bins == b)) for b in np.unique(bins)}
        w = np.array([1.0 / counts[b] for b in bins])
        rayleigh['nu_H_per_A'] = float(2 * np.sum(w * I[ok] * (Le[ok] - rayleigh['L0_H'])) / np.sum(w * I[ok] ** 2))
    T_model = small_signal_period(g_static, R0 + R_onset, c1_small, C2, rayleigh['L0_H'], rayleigh['r0_ohm'])
    R_hopf_model = small_signal_hopf(g_static, c1_small, C2, rayleigh['L0_H'], rayleigh['r0_ohm'])
    # the measured onset: the squared amplitude of the small cycles is linear in R and vanishes at the Hopf point
    amp = np.array([0.5 * (r['v1_max'] - r['v1_min']) for r in small])
    rr = np.array([r['rpot'] for r in small])
    m_on = (amp > 0.15) & (amp < 1.0)
    if m_on.sum() >= 4:
        k_on, c_on = np.polyfit(rr[m_on], amp[m_on] ** 2, 1)
        R_hopf_meas = float(-c_on / k_on)
    else:
        R_hopf_meas = float('nan')

    # ---- 5. the dynamic deviation of the diode in the large cycle ----
    dev = {}
    if large:
        r = max(large, key=lambda r: r['rpot'])
        v1, v2, dv1, dv2, T = r['cycle']
        vk, ik = pwl_knots(pwl['bp_V'], pwl['G_S'])
        inr = (v2 - v1) / r['Rt'] - c1_small * dv1
        d = inr - np.interp(v1, vk, ik)
        m = np.abs(v1) < 3
        dev = dict(record=r['label'], rpot=r['rpot'], peak_mA=float(np.max(np.abs(d[m])) * 1e3),
                   max_dv1dt_V_per_us=float(np.max(np.abs(dv1)) * 1e-6))

    # ---- report ----
    result = dict(
        records_analysed=len(rows), periodic_records=len(per),
        divider_gain_check=dict(mean_A_plus_B=float(gain.mean()), max_deviation=float(np.max(np.abs(gain - 1)))),
        C1_nF=dict(small_signal_loop=c1_small * 1e9, period1_loop=c1_p1 * 1e9, large_cycle_loop=c1_lc * 1e9,
                   kcl_core=fits['core']['C1_apparent_nF'], kcl_large_cycle=fits['large_cycle']['C1_apparent_nF']),
        C2_nF=float(c2_nF),
        tank_admittance_fit=dict(C2_nF=float(c2_adm_nF), L_mH=float(L_fit_mH), r_ohm=float(r_fit),
                                 relative_residual=tank_resid, harmonics=len(hr)),
        C2_source=c2_source,
        tank_admittance_near_origin=dict(C2_nF=float(c2_nF), L_mH=float(L_tiny_mH), r_ohm=float(r_tiny),
                                         relative_residual=float(tiny_resid), harmonics=len(hr_tiny)),
        small_signal_check=dict(onset_period_us=T_onset * 1e6, onset_records=len(onset_T), onset_rpot_mean=R_onset,
                                model_period_us=T_model * 1e6,
                                hopf_rpot_model=None if R_hopf_model is None else R_hopf_model - R0,
                                hopf_rpot_measured=R_hopf_meas),
        inductor=dict(rayleigh=rayleigh, per_record=[dict(label=t[0], rpot=t[1], T_us=t[2], I_amp_mA=t[3] * 1e3,
                                                          L_eff_mH=t[4] * 1e3, r_eff_ohm=t[5]) for t in tank]),
        diode=dict(pwl=pwl, kennedy=kennedy,
                   kcl_core=dict(records=fits['core']['records'], relative_residual=fits['core']['relative_residual'],
                                 **fits['core']['pwl']),
                   kcl_large_cycle=dict(records=fits['large_cycle']['records'],
                                        relative_residual=fits['large_cycle']['relative_residual'],
                                        **fits['large_cycle']['pwl']),
                   large_cycle_dynamic_deviation=dev),
        notes=[
            'C1 (small_signal_loop) is the capacitance node 1 shows at small amplitude; the growth with amplitude and the large-cycle deviation are the diode op-amps (see RESULTS.md).',
            'inductor: L_eff and r_eff are the secant inductance and the loss resistance of the measured flux-current loop; the Rayleigh law L(i) = L0 + nu|i|, r(i) = r0 + rho|i| reproduces their growth with amplitude.',
            'diode: pwl is the five-segment law as the circuit sees it (inner three lines from the core records, saturation lines from the large cycle); kennedy maps it onto the two-op-amp realisation.'])
    with open(os.path.join(a.out, 'identified.json'), 'w') as fh:
        json.dump(result, fh, indent=2)

    lines = [f'identify.py: {len(rows)} records, {len(per)} periodic ({len(small)} period-1, {len(tiny)} near-origin, {len(large)} large cycle)',
             f'channel gains: divider A + B = {gain.mean():.4f} (max deviation {np.max(np.abs(gain-1)):.4f}) -> CH1 and CH2 gains equal',
             f'C1 from the loop integral: {c1_small*1e9:.2f} nF at small amplitude, {c1_p1*1e9:.2f} nF on the period-1 cycle, {c1_lc*1e9:.2f} nF on the large cycle',
             f'C1 from the KCL fits: {fits["core"]["C1_apparent_nF"]:.2f} nF (core records), {fits["large_cycle"]["C1_apparent_nF"]:.2f} nF (large cycle)',
             f'C2 = {c2_nF:.1f} nF from the {c2_source} (L = {L_tiny_mH:.1f} mH, r = {r_tiny:.1f} ohm there, residual {tiny_resid:.3f});'
             f' the period-1 cycles alone give C2 = {c2_adm_nF:.1f} nF, L = {L_fit_mH:.1f} mH, r = {r_fit:.1f} ohm ({len(hr)} harmonics, residual {tank_resid:.3f})',
             f'inductor, Rayleigh law: L0 = {rayleigh["L0_H"]*1e3:.2f} mH (from the onset period {T_onset*1e6:.1f} us at {R_onset:.0f} ohm, {len(onset_T)} records;'
             f' the flux loops extrapolate to {rayleigh["L0_loops_H"]*1e3:.2f} mH), nu = {rayleigh["nu_H_per_A"]:.2f} H/A'
             f' (loops alone {rayleigh["nu_loops_H_per_A"]:.2f}), r0 = {rayleigh["r0_ohm"]:.2f} ohm, rho = {rayleigh["rho_ohm_per_A"]:.0f} ohm/A',
             f'small-signal check: the linearised circuit oscillates at {T_model*1e6:.1f} us at {R_onset:.0f} ohm and its Hopf point is at'
             f' {R_hopf_model - R0 if R_hopf_model else float("nan"):.0f} ohm; the squared amplitude of the small cycles extrapolates to zero at {R_hopf_meas:.0f} ohm',
             f'{"record":22s} {"Rpot":>7} {"T us":>6} {"I mA":>6} {"L_eff mH":>8} {"r_eff":>6}']
    for t in sorted(tank, key=lambda t: t[3]):
        lines.append(f'{t[0]:22s} {t[1]:7.1f} {t[2]:6.1f} {t[3]*1e3:6.2f} {t[4]*1e3:8.2f} {t[5]:6.1f}')
    lines += [f'diode PWL: breakpoints {np.round(pwl["bp_V"], 3).tolist()} V, slopes {np.round(np.array(pwl["G_S"])*1e3, 4).tolist()} mS',
              f'  inner three lines: KCL fit over {fits["core"]["records"]} core records (residual {fits["core"]["relative_residual"]:.4f}); '
              f'saturation lines: plateaus of {pwl["saturation_records"]} large-cycle records; '
              f'KCL fit over the large cycle alone: residual {fits["large_cycle"]["relative_residual"]:.4f}, apparent C1 {fits["large_cycle"]["C1_apparent_nF"]:.2f} nF',
              f'Kennedy diode: RA = {kennedy["RA_ohm"]:.0f} ohm, AA = {kennedy["AA"]:.3f}, RB = {kennedy["RB_ohm"]:.0f} ohm, AB = {kennedy["AB"]:.2f}, '
              f'rails stage A {kennedy["VpA_V"]:+.2f}/{kennedy["VnA_V"]:+.2f} V, stage B {kennedy["VpB_V"]:+.2f}/{kennedy["VnB_V"]:+.2f} V']
    if dev:
        lines.append(f'large-cycle diode deviation ({dev["record"]}): up to {dev["peak_mA"]:.2f} mA inside |v1| < 3 V at {dev["max_dv1dt_V_per_us"]:.2f} V/us')
    txt = '\n'.join(lines)
    print(txt)
    with open(os.path.join(a.out, 'identified.txt'), 'w') as fh:
        fh.write(txt + '\n')

    # ---- figures ----
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    ax = axes[0]
    show = sorted(per, key=lambda r: r['I_amp'])
    pick = [show[0], show[len(show) // 3], show[2 * len(show) // 3], show[-1]] if len(show) >= 4 else show
    for r in pick:
        ax.plot(r['iL'] * 1e3, r['phi'] * 1e3, lw=1, label=f'{r["label"]} ({r["rpot"]:.0f} ohm), I = {r["I_amp"]*1e3:.1f} mA')
    ax.set_xlabel('i_L (mA)'); ax.set_ylabel('flux = int v2 dt (mWb)'); ax.set_title('inductor: flux against current, one averaged cycle'); ax.legend(fontsize=7); ax.grid(alpha=.3)
    ax = axes[1]
    ax.plot(I * 1e3, Le * 1e3, 'o', ms=4, label='secant inductance of the loop')
    xx = np.linspace(0, I.max() * 1e3, 50)
    ax.plot(xx, (rayleigh['L0_H'] + 0.5 * rayleigh['nu_H_per_A'] * xx * 1e-3) * 1e3, '-', label=f'L0 + nu I/2, L0 = {rayleigh["L0_H"]*1e3:.1f} mH')
    ax.set_xlabel('current amplitude (mA)'); ax.set_ylabel('L_eff (mH)'); ax.legend(fontsize=8); ax.grid(alpha=.3); ax.set_title('effective inductance against amplitude')
    ax = axes[2]
    ax.plot(I * 1e3, re, 'o', ms=4, label='loss resistance of the loop')
    ax.plot(xx, rayleigh['r0_ohm'] + (8 / (3 * np.pi)) * rayleigh['rho_ohm_per_A'] * xx * 1e-3, '-', label=f'r0 + 0.85 rho I, r0 = {rayleigh["r0_ohm"]:.1f} ohm')
    ax.set_xlabel('current amplitude (mA)'); ax.set_ylabel('r_eff (ohm)'); ax.legend(fontsize=8); ax.grid(alpha=.3); ax.set_title('loss resistance against amplitude')
    fig.tight_layout(); fig.savefig(os.path.join(a.out, 'identified_inductor.png'), dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 6))
    vi = os.path.join(HERE, '..', 'trace1.csv')
    if os.path.exists(vi):
        d = np.genfromtxt(vi, delimiter=',', skip_header=1)
        i_m = (d[:, 1] - d[:, 2]) / 216.0
        off = np.nanmedian(i_m[np.abs(d[:, 2]) < 0.15])
        ax.plot(d[:, 2], (i_m - off) * 1e3, '.', ms=1, alpha=0.1, color='steelblue', label='V-I trace (trace1.csv, offset removed)')
    for r, c in [(small[len(small) // 2] if small else None, 'C1'), (max(large, key=lambda r: r['rpot']) if large else None, 'C3')]:
        if r is None:
            continue
        v1, v2, dv1, dv2, T = r['cycle']
        inr = (v2 - v1) / r['Rt'] - c1_small * dv1
        ax.plot(v1, inr * 1e3, '.', ms=1.5, color=c, label=f'i_NR from {r["label"]} ({r["rpot"]:.0f} ohm)')
    vgrid = fits['core']['curve_v']
    ax.plot(vgrid, fits['core']['curve_i'] * 1e3, 'k-', lw=1.2, label='KCL spline fit, core records')
    ax.plot(fits['large_cycle']['curve_v'], fits['large_cycle']['curve_i'] * 1e3, 'k--', lw=1, label='KCL spline fit, large cycle')
    for b in pwl['bp_V']:
        ax.axvline(b, color='0.6', lw=0.6, ls=':')
    ax.set_xlim(-8, 7.4); ax.set_ylim(-4, 4); ax.set_xlabel('v1 (V)'); ax.set_ylabel('current into N_R (mA)')
    ax.set_title('the nonlinear element as the circuit sees it'); ax.legend(fontsize=8); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, 'identified_diode.png'), dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    xs = [np.ptp([r['v1_min'], r['v1_max']]) for r in per]
    ax.plot(xs, [r['C1_loop'] * 1e9 for r in per], 'o', ms=4)
    ax.set_xlabel('peak-to-peak v1 of the averaged cycle (V)'); ax.set_ylabel('loop-integral C1 (nF)')
    ax.set_title('apparent C1 against amplitude (static element assumed)'); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, 'identified_c1.png'), dpi=150); plt.close(fig)
    print('written identified.json, identified.txt and the three figures to', a.out)


if __name__ == '__main__':
    main()
