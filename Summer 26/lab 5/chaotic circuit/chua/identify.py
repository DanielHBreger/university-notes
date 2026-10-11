"""
Identify the circuit that was actually on the bench from the oscillator records.

Every quantity here comes from the sweep records themselves (V1, V2 and the
divider-fitted Rpot), not from the nominal component values. The periodic
records (the period-1 cycle below the Hopf point, the large outer cycle and the
small orbit near the origin) are averaged over their hundreds of cycles, which
reduces the scope's quantisation scatter and leaves waveforms clean enough to
differentiate. The circuit equations then give each element directly:

  node 1   C1 v1' = (v2 - v1)/Rt - g(v1)
           For a static element the loop integral of g dv1 over a cycle is
           zero, so   C1 = oint (v2 - v1)/Rt v1' dt / oint v1'^2 dt   for ANY g.
           i_NR = (v2 - v1)/Rt - C1 v1' against v1 is then the element's curve
           as the circuit sees it, and any loop in it is a dynamic effect.
           In general the loop integral is the v1'^2-weighted mean of node 1's
           incremental capacitance along the record, so it is taken over EVERY
           record (averaged cycles, or finite-difference windows for the
           chaotic ones), and node1_law fits how it varies: the capacitor
           grows with the swing of v1 since its last turning point (Rayleigh
           law of a ferroelectric ceramic), and stage B of the diode adds its
           lag's share A_B tau_B / R_B while it is linear.
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
    identified_c1.png         node 1's capacitance per record with its law

Usage:
    python identify.py               # forward/ and back/ beside this script
    python identify.py --quick       # every 4th record (about 2 min)
"""
import argparse
import json
import os
from collections import namedtuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import BSpline
from scipy.optimize import brentq, least_squares

from scope_data import R0, read_scope
from lorenz_map import maxima_from_arrays
from integration import GBW, TAU_R, pwl_knots, swing

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
    """Complex Fourier amplitudes of one period of x."""
    return np.fft.rfft(x) / len(x)


def derivative(x, T, kmax=KMAX):
    """Time derivative of one period of x (period T) from its first kmax harmonics."""
    X = np.fft.rfft(x)
    k = np.arange(len(X))
    X = X * (1j * 2 * np.pi * k / T)
    X[k > kmax] = 0
    return np.fft.irfft(X, len(x))


def lowpass(x, kmax=KMAX):
    """One period of x with only its first kmax harmonics."""
    X = np.fft.rfft(x)
    X[kmax + 1:] = 0
    return np.fft.irfft(X, len(x))


def loop_c1(v1, v2, dv1, Rt):
    """C1 from the loop integral: exact for any static element g(v1)."""
    ir = (v2 - v1) / Rt
    return float(np.sum(ir * dv1) / np.sum(dv1 * dv1))


# ---- node 1's capacitance in every record ---------------------------------------
# Multiplying node 1's current balance by v1' and integrating, the static element
# drops out: int (v2 - v1)/Rt v1' dt / int v1'^2 dt is the v1'^2-weighted mean of
# node 1's incremental capacitance along the record. Two estimators keep the scope's
# noise out of the denominator; both were checked on scope-like simulated records.
MARKER_FC = 20e3    # Hz, low-pass of v1 for the cycle markers
WINDOW_S = 20e-6    # s, window of the finite-difference estimator
Node1 = namedtuple('Node1', 'c weights v1 swing method')


def lowpassed(x, dt, fc):
    """x with every Fourier component above fc removed."""
    X = np.fft.rfft(x - x.mean())
    X[np.fft.rfftfreq(len(x), dt) > fc] = 0
    return np.fft.irfft(X, len(x)) + x.mean()


def cycle_node1(t, v1, v2, Rt):
    """
    Periodic records: cycles marked on v1 low-passed at MARKER_FC (markers on the
    raw v1 pick up its noise and sharpen the average, inflating int v1'^2); the
    numerator from the average of all cycles, the denominator from the product of
    the derivatives averaged over the first and the second half of the cycles,
    whose noise does not correlate. None if the record cannot be averaged.
    """
    ref = lowpassed(v1, t[1] - t[0], MARKER_FC)
    _, info = maxima_from_arrays(t, ref)
    if len(info['t_peaks']) < 60:
        return None
    tp = cycle_markers(t, ref, np.mean(np.diff(info['t_peaks'])))
    data = np.column_stack([v1, v2])
    m = len(tp) // 2
    parts = [averaged_cycle(t, data, tp), averaged_cycle(t, data, tp[:m + 1]), averaged_cycle(t, data, tp[m:])]
    if len(tp) < 60 or any(p is None for p in parts):
        return None
    (cyc, T, _), (h1, T1, _), (h2, T2, _) = parts
    a, b = lowpass(cyc[:, 0]), lowpass(cyc[:, 1])
    w = derivative(h1[:, 0], T1) * derivative(h2[:, 0], T2)
    c = float(np.sum((b - a) / Rt * derivative(cyc[:, 0], T)) / np.sum(w))
    reps = int(np.ceil(10 * TAU_R / T)) + 2                  # the swing memory settles over repeated cycles
    x = swing(np.tile(a, reps), T / len(a))[-len(a):]
    return Node1(c, w, a, x, 'cycle')


def window_node1(t, v1, v2, Rt):
    """
    Any record: over windows of WINDOW_S, int_w (v2 - v1)/Rt dt = Delta Q + int_w g dt,
    and sum y_w Delta v_w / sum Delta v_w^2 is the loop integral with finite
    differences, so noise enters Delta v only through two end values. The record
    is cut where v1 passes its starting value in the starting direction, which
    removes the boundary term (the DC current through R times the net change of v1).
    """
    dt = float(np.median(np.diff(t)))
    n = max(2, int(round(WINDOW_S / dt)))
    q = cumulative_trapezoid((v2 - v1) / Rt, dx=dt, initial=0)
    v_end = np.convolve(v1, np.ones(3) / 3, mode='same')
    vf = lowpassed(v1, dt, 100e3)
    rising = vf[n + 1] > vf[n]
    x = vf - vf[n]
    cross = np.flatnonzero((x[:-1] < 0) & (x[1:] >= 0) if rising else (x[:-1] > 0) & (x[1:] <= 0))
    cross = cross[cross > 2 * n]
    last = n + ((cross[-1] if len(cross) else len(t) - 2) - n) // n * n
    i0 = np.arange(n, last - n + 1, n)
    dv = v_end[i0 + n] - v_end[i0]
    c = float(np.sum((q[i0 + n] - q[i0]) * dv) / np.sum(dv * dv))
    mid = i0 + n // 2
    return Node1(c, dv * dv, vf[mid], swing(vf, dt)[mid], 'window')


def node1_features(est, knots, bp):
    """Weighted means over a record: min(x, k) and max(x - k, 0) per knot, stage B linear, stage A saturated."""
    w = est.weights / est.weights.sum()
    v, x = est.v1, est.swing
    out = dict(swing=float(w @ x), B_lin=float(w @ ((v > bp[1]) & (v < bp[2]))),
               A_sat=float(w @ ((v < bp[0]) | (v > bp[3]))))
    out.update({('lo', k): float(w @ np.minimum(x, k)) for k in knots})
    out.update({('hi', k): float(w @ np.maximum(x - k, 0.0)) for k in knots})
    return out


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
    w = np.array([r[0] for r in rows])
    Y = np.array([r[1] for r in rows])
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
    """The cubic B-spline basis for g(v1), evaluated at every sample."""
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
    """(apparent C1, spline coefficients, relative residual), with a second-difference penalty on the spline."""
    D = np.diff(np.eye(NB), n=2, axis=0)
    P = np.zeros_like(XtX)
    P[1:, 1:] = D.T @ D * lam
    coef = np.linalg.solve(XtX + P, Xty)
    resid2 = yty - 2 * coef @ Xty + coef @ XtX @ coef
    return coef[0], coef[1:], float(np.sqrt(max(resid2, 0) / yty))


def pwl_from_spline(cs):
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


def bin_weights(I, bin_mA=1.0):
    """Weights that give every current-amplitude bin of `bin_mA` the same total weight."""
    bins = np.floor(I * 1e3 / bin_mA).astype(int)
    counts = {b: int(np.sum(bins == b)) for b in np.unique(bins)}
    return np.array([1.0 / counts[b] for b in bins])


def rayleigh_fit(I, L_eff, r_eff, bin_mA=1.0):
    """
    L_eff = L0 + nu I / 2  and  r_eff = r0 + (8/3pi) rho I: the secant
    inductance and the loss resistance of L = L0 + nu|d| with the loss voltage
    r0 i + rho|d|d, for a current excursion d = I sin(wt) about the record's
    mean (integration.py). Every amplitude bin of `bin_mA` carries the
    same total weight, so the forty-odd large-cycle records do not outvote the
    small cycles that fix the small-signal end of the law.
    """
    I, L_eff, r_eff = np.asarray(I, float), np.asarray(L_eff, float), np.asarray(r_eff, float)
    w = bin_weights(I, bin_mA)
    a = np.polyfit(I, L_eff, 1, w=np.sqrt(w))
    b = np.polyfit(I, r_eff, 1, w=np.sqrt(w))
    return dict(L0_H=float(a[1]), nu_H_per_A=float(2 * a[0]),
                r0_ohm=float(b[1]), rho_ohm_per_A=float(b[0] / (8 / (3 * np.pi))))


# ---- the linearised circuit ------------------------------------------------
def negative_equilibrium(g, Rt, r0):
    """v1 of the circuit's negative equilibrium with a static element g, or None."""
    def h(v):
        return g(v) + v / (Rt + r0)
    vs = np.linspace(-15, -1e-3, 20001)
    hv = h(vs)
    idx = np.flatnonzero(hv[:-1] * hv[1:] < 0)
    if not len(idx):
        return None
    return brentq(h, vs[idx[-1]], vs[idx[-1] + 1])


def complex_pair(g, Rt, c1, c2, l0, r0):
    """Eigenvalues with an imaginary part of the linearisation about the negative equilibrium (None without one)."""
    v = negative_equilibrium(g, Rt, r0)
    if v is None:
        return None
    dv = 1e-4
    G = (g(v + dv) - g(v - dv)) / (2 * dv)
    J = np.array([[(-1 / Rt - G) / c1, 1 / (Rt * c1), 0.0],
                  [1 / (Rt * c2), -1 / (Rt * c2), -1 / c2],
                  [0.0, 1 / l0, -r0 / l0]])
    ev = np.linalg.eigvals(J)
    return ev[np.abs(ev.imag) > 1.0]


def small_signal_period(g, Rt, c1, c2, l0, r0):
    """Period of the linearisation about the negative equilibrium of the circuit with a static element g."""
    ev = complex_pair(g, Rt, c1, c2, l0, r0)
    return 2 * np.pi / np.abs(ev.imag).max() if ev is not None and len(ev) else np.nan


def small_signal_hopf(g, c1, c2, l0, r0, rt_range=(1200.0, 3200.0)):
    """Rt at which the negative equilibrium's complex pair crosses the axis (None if never)."""
    def growth_rate(Rt):
        ev = complex_pair(g, Rt, c1, c2, l0, r0)
        return ev.real.max() if ev is not None and len(ev) else np.nan
    Rs = np.arange(rt_range[0], rt_range[1], 10.0)
    vals = np.array([growth_rate(R) for R in Rs])
    idx = [k for k in np.flatnonzero(np.sign(vals[:-1]) != np.sign(vals[1:]))
           if np.isfinite(vals[k]) and np.isfinite(vals[k + 1])]
    return brentq(growth_rate, Rs[idx[-1]], Rs[idx[-1] + 1]) if idx else None


# ---- records ---------------------------------------------------------------
def onset_table():
    """onset.py's records at the top of the dial (onset_records.csv) with its two thresholds."""
    from onset import CAL_AMP_V, OSC_MIN_V
    top = pd.read_csv(os.path.join(HERE, 'onset_records.csv'))
    return top[top.sidecar_status.eq('ok')], CAL_AMP_V, OSC_MIN_V


def sweep_records(sweep):
    """
    (path, Rpot, label) of every record rpot.py accepts. rpot.py reads the small
    cycles at the top of the dial low (N1, onset.py), so those records take
    onset.py's resistance, as in the report's figures.
    """
    tab = pd.read_csv(os.path.join(HERE, f'{sweep}_rpot.csv'))
    tab = tab[tab.status.eq('ok')]
    top, cal_amp, osc_min = onset_table()
    top = top[top.sweep.eq(sweep) & (top.v1_half_range_V < cal_amp) & (top.A1_mV >= 1e3 * osc_min)]
    fix = dict(zip(top.filename, top.rpot_means_ohm.astype(float)))
    return [(os.path.join(HERE, sweep, r.filename), fix.get(r.filename, float(r.rpot_ohm)), f'{sweep}/{r.filename}')
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
    # Divider consistency check; A+B alone cannot identify relative channel gains.
    X = np.column_stack([d[:, 0] - d[:, 0].mean(), d[:, 1] - d[:, 1].mean()])
    A, B = np.linalg.lstsq(X, d[:, 2] - d[:, 2].mean(), rcond=None)[0]
    out['divider_A_plus_B'] = float(A + B)
    out['periodic'] = out['jitter'] < 0.015 and out['max_spread'] < 0.12
    res = averaged_cycle(t, d[:, :2], tp) if out['periodic'] else None
    node1 = cycle_node1(t, d[:, 0], d[:, 1], out['Rt']) if res is not None else None
    if node1 is None:
        out['periodic'] = False
        out['node1'] = window_node1(t, d[:, 0], d[:, 1], out['Rt'])
        return out
    cyc, T, ncyc = res
    v1, v2 = lowpass(cyc[:, 0]), lowpass(cyc[:, 1])
    dv1, dv2 = derivative(cyc[:, 0], T), derivative(cyc[:, 1], T)
    out.update(T_cyc_us=T * 1e6, n_cycles=ncyc, v1_min=float(v1.min()), v1_max=float(v1.max()),
               v2_amp=float(0.5 * np.ptp(v2)), C1_loop=node1.c, node1=node1)
    out['cycle'] = (v1, v2, dv1, dv2, T)
    out['harm'] = admittance_rows(v1, v2, R0 + rpot, T)
    return out


def record_number(path):
    """N of a record named traceN.csv."""
    return int(os.path.basename(path)[5:-4])


# ---- the steps of main ---------------------------------------------------------
def analyse_sweeps(step):
    """analyse_record on every accepted record (every step-th, but all of forward from trace402 on)."""
    rows = []
    for sweep in ['forward', 'back']:
        for k, (path, rpot, label) in enumerate(sweep_records(sweep)):
            if k % step and not (sweep == 'forward' and record_number(path) >= 402):
                continue
            r = analyse_record(path, rpot)
            if r is None:
                continue
            r['label'] = label
            rows.append(r)
            print(f'{label:22s} R={rpot:7.2f}  T={r["T_us"]:6.1f} us  '
                  + (f'C1 {r["C1_loop"]*1e9:6.3f} nF  v1 {r["v1_min"]:6.2f}..{r["v1_max"]:6.2f}'
                     if r['periodic'] else 'not periodic'), flush=True)
    return rows


def loop_c1_values(small, large, top):
    """
    The loop-integral C1 (F) at small amplitude and on larger cycles. Node 1's
    apparent capacitance grows with the amplitude of the cycle (node1_law), so
    the small-signal value is taken from the cycles onset.py's A1^2 fit uses
    (A1 < FIT_MAX_V), the smallest swing recorded.
    """
    from onset import FIT_MAX_V
    near = {f'{s}/{f}' for s, f, a1 in zip(top.sweep, top.filename, top.A1_mV) if a1 < 1e3 * FIT_MAX_V}
    onset = [r for r in small if r['label'] in near]
    c1_onset = [r['C1_loop'] for r in onset]
    return dict(small=float(np.mean(c1_onset)) if onset else np.nan,
                small_se=float(np.std(c1_onset, ddof=1) / np.sqrt(len(onset))) if len(onset) > 1 else np.nan,
                small_records=[r['label'] for r in onset], fit_max_V=FIT_MAX_V,
                shoulder=float(np.mean([r['C1_loop'] for r in small if r['v1_max'] < -2.0])),
                period1=float(np.mean([r['C1_loop'] for r in small if r['v1_max'] > -1.0])),
                large_cycle=float(np.mean([r['C1_loop'] for r in large])))


LAW_KNOTS = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2)     # V, candidate knots of the node-1 law


def node1_law(rows, small_labels, pwl, kennedy):
    """
    Node 1's capacitance as a law in the swing x of v1 since its last turning
    point (integration.swing), fitted to the loop-integral capacitance of every
    record:
        C_rec = C_s + a0 <min(x, k)> + a1 <max(x - k, 0)> + s_B <B linear> + s_A <A linear>,
    <> the estimator's weighted mean over the record. s_A = tau_A AA/RA is stage
    A's lag at the datasheet speed (the records cannot tell it from C_s); C_s is
    the capacitor's share on the onset cycles, the smallest swing recorded; a0,
    a1 and s_B are least squares, and the knot k is the one with the least
    residual. A lag tau_B of stage B acts as a capacitance s_B = tau_B AB/RB only
    while stage B is linear, inside the inner breakpoints. The fit is repeated on
    the periodic records alone, predicting the chaotic ones, as a check.
    """
    feats = [node1_features(r['node1'], LAW_KNOTS, pwl['bp_V']) for r in rows]
    C = np.array([r['node1'].c for r in rows])
    periodic = np.array([r['node1'].method == 'cycle' for r in rows])
    s_a = kennedy['AA'] / (2 * np.pi * GBW) * kennedy['AA'] / kennedy['RA_ohm']
    lin_a = 1.0 - np.array([f['A_sat'] for f in feats])
    b_lin = np.array([f['B_lin'] for f in feats])
    onset = np.array([r['label'] in small_labels for r in rows])
    c_s_records = C[onset] - s_a * lin_a[onset]
    c_s = float(c_s_records.mean())
    y = C - s_a * lin_a - c_s

    def solve(k, idx):
        A = np.column_stack([[f[('lo', k)] for f in feats], [f[('hi', k)] for f in feats], b_lin])
        p = np.linalg.lstsq(A[idx], y[idx], rcond=None)[0]
        r = y - A @ p
        cov = np.linalg.inv(A[idx].T @ A[idx]) * np.sum(r[idx] ** 2) / (idx.sum() - A.shape[1])
        return p, np.sqrt(np.diag(cov)), r

    everything = np.ones(len(rows), bool)
    fits = {k: solve(k, everything) for k in LAW_KNOTS}
    k = min(fits, key=lambda knot: np.sum(fits[knot][2] ** 2))
    (a0, a1, s_b), (e0, e1, eb), r = fits[k]
    (_, a1_p, s_b_p), _, r_p = solve(k, periodic)
    rpot = np.array([r_['rpot'] for r_ in rows])
    groups = {'onset cycles': onset,
              'other periodic, V1 < 5.5 V': periodic & ~onset & np.array([r_.get('v1_max', 0) < 5.5 for r_ in rows]),
              'large cycle': periodic & np.array([r_.get('v1_max', 0) > 5.5 for r_ in rows]),
              'chaotic, above 669 ohm': ~periodic & (rpot > 669),
              'chaotic, 669 ohm and below': ~periodic & (rpot <= 669)}
    return dict(
        C_s_nF=c_s * 1e9, C_s_se_nF=float(c_s_records.std(ddof=1) / np.sqrt(len(c_s_records))) * 1e9,
        knot_V=k, a0_nF_per_V=a0 * 1e9, a0_se_nF_per_V=e0 * 1e9, a1_nF_per_V=a1 * 1e9, a1_se_nF_per_V=e1 * 1e9,
        stage_A_share_nF=s_a * 1e9, stage_B_share_nF=s_b * 1e9, stage_B_share_se_nF=eb * 1e9,
        stage_B_lag_us=s_b * kennedy['RB_ohm'] / kennedy['AB'] * 1e6,
        stage_B_lag_se_us=eb * kennedy['RB_ohm'] / kennedy['AB'] * 1e6,
        residual_rms_nF=float(np.sqrt(np.mean(r ** 2))) * 1e9, records=len(rows), periodic_records=int(periodic.sum()),
        residual_by_group_nF={g: float(np.mean(r[m])) * 1e9 for g, m in groups.items() if m.any()},
        periodic_only=dict(a1_nF_per_V=a1_p * 1e9, stage_B_share_nF=s_b_p * 1e9,
                           chaotic_rms_nF=float(np.sqrt(np.mean(r_p[~periodic] ** 2))) * 1e9),
        knot_scan_rms_nF={str(kk): float(np.sqrt(np.mean(f[2] ** 2))) * 1e9 for kk, f in fits.items()},
        per_record=dict(rpot=rpot, C=C, capacitor=C - s_a * lin_a - s_b * b_lin, swing=np.array([f['swing'] for f in feats]),
                        law=C - r, periodic=periodic))


def kcl_records(step):
    """The records of the two KCL fits: 'core' (the inner three lines) and 'large_cycle' (a check)."""
    sets = {'core': [], 'large_cycle': []}
    for sweep in ['forward', 'back']:
        for path, rpot, label in sweep_records(sweep):
            n = record_number(path)
            if sweep == 'forward' and 330 < rpot < 900 and n % (2 * step) == 0:
                sets['core'].append((path, rpot))
            if ((sweep == 'forward' and rpot < 315) or (sweep == 'back' and rpot < 680)) and n % step == 0:
                sets['large_cycle'].append((path, rpot))
    return sets


def kcl_fit(records):
    """Integrated-KCL fit of g(v1) over a set of records, and the five-segment numbers read off it."""
    XtX, Xty, yty = np.zeros((NB + 1, NB + 1)), np.zeros(NB + 1), 0.0
    for path, rpot in records:
        A, b, yy = kcl_normal_equations(path, rpot)
        XtX += A
        Xty += b
        yty += yy
    c1, cs, resid = solve_kcl(XtX, Xty, yty)
    pw, vgrid, ggrid = pwl_from_spline(cs)
    return dict(records=len(records), C1_apparent_nF=c1 * 1e9, relative_residual=resid, pwl=pw,
                curve_v=vgrid, curve_i=ggrid)


def five_segment_law(inner, large, c1_kcl):
    """
    (pwl, kennedy): the inner three lines of the core KCL fit, the saturation
    lines from the plateaus of every large-cycle record (median over records),
    and the outer breakpoints where the shoulder lines through the inner
    breakpoints meet the saturation lines.
    """
    sat = {-1: [], 1: []}
    for r in large:
        v1, v2, dv1, dv2, T = r['cycle']
        for side in (-1, 1):
            f = saturation_from_plateau(v1, v2, dv1, r['Rt'], c1_kcl, side)
            if f:
                sat[side].append(f)
    if not sat[-1] or not sat[1]:
        raise SystemExit('no large-cycle plateau to fit the saturation segments from')
    Gc_left, c_left = np.median(sat[-1], axis=0)
    Gc_right, c_right = np.median(sat[1], axis=0)
    Ga, GbL, GbR = inner['Ga'], inner['Gb_left'], inner['Gb_right']
    bpL, bpR = inner['bp_in_left'], inner['bp_in_right']
    cL = Ga * bpL - GbL * bpL                    # i = GbL v + cL on the left shoulder
    cR = Ga * bpR - GbR * bpR
    bp_out_left = (c_left - cL) / (GbL - Gc_left)
    bp_out_right = (c_right - cR) / (GbR - Gc_right)
    pwl = dict(bp_V=[float(bp_out_left), bpL, bpR, float(bp_out_right)],
               G_S=[float(Gc_left), GbL, Ga, GbR, float(Gc_right)],
               saturation_records=[len(sat[-1]), len(sat[1])])
    Gb = 0.5 * (GbL + GbR)
    Gc = 0.5 * (Gc_left + Gc_right)
    kennedy = kennedy_from_pwl(Ga, Gb, Gc, bpL, bpR, bp_out_left, bp_out_right)
    return pwl, kennedy


def tank_fits(small, tiny):
    """
    (period-1 fit, near-origin fit, source of C2) of the tank admittance. The
    period-1 fit pins L*C2 but splits it poorly between C2 and L: at 3 kHz and
    2 mA the inductor is already nonlinear over the harmonics the fit uses. The
    small orbit around the origin (about 1.4 kHz, 0.7 mA, ten harmonics) is the
    record where the inductor is most nearly linear and the harmonics reach
    highest, so C2 is taken from it when it exists.
    """
    hr = [h for r in small if len(r['harm']) >= 5 for h in r['harm']]
    (c2_nF, L_mH, r_ohm), resid = fit_tank(hr)
    period1 = dict(C2_nF=float(c2_nF), L_mH=float(L_mH), r_ohm=float(r_ohm), relative_residual=resid,
                   harmonics=len(hr))
    hr_tiny = [h for r in tiny if len(r['harm']) >= 6 for h in r['harm']]
    if hr_tiny:
        (c2_nF, L_mH, r_ohm), resid = fit_tank(hr_tiny)
        source = f'near-origin orbit ({len(tiny)} records, {len(hr_tiny)} harmonics)'
    else:
        source = 'period-1 admittance fit (no near-origin record)'
    near_origin = dict(C2_nF=float(c2_nF), L_mH=float(L_mH), r_ohm=float(r_ohm), relative_residual=float(resid),
                       harmonics=len(hr_tiny))
    return period1, near_origin, source


Loop = namedtuple('Loop', 'label rpot T_us I L_eff r_eff')


def inductor_loops(per, C2):
    """The flux-current loop of every periodic record (also kept on the record, for the figure)."""
    loops = []
    for r in per:
        v1, v2, dv1, dv2, T = r['cycle']
        iL, phi, I, L_eff, r_eff = inductor_loop(v1, v2, dv2, r['Rt'], C2, T)
        r['iL'], r['phi'] = iL, phi
        r.update(I_amp=I, L_eff=L_eff, r_eff=r_eff)
        loops.append(Loop(r['label'], r['rpot'], r['T_cyc_us'], I, L_eff, r_eff))
    return loops


def onset_record(top, osc_min):
    """The first limit cycle of the forward sweep (onset.py): its label, period and resistance."""
    first = top[top.sweep.eq('forward') & (top.A1_mV >= 1e3 * osc_min)].sort_values('rpot_means_ohm').iloc[-1]
    return dict(label=f'forward/{first.filename}', T=1.0 / float(first.f_Hz), R=float(first.rpot_means_ohm))


def calibrate_l0(rayleigh, I, L_eff, g_static, c1, C2, onset):
    """
    Close the small-signal end of the inductance law on the onset frequency.
    The flux loops of the smallest cycles are only a few quantisation steps
    wide, so their extrapolation to zero amplitude is the least certain number
    here. L0 is set so that the linearised circuit (small-signal C1, C2, r0,
    the fitted element) oscillates at the frequency of the first limit cycle
    at its resistance; nu is then refitted to the loops with that L0. The
    loop-only values stay in the dict as L0_loops_H and nu_loops_H_per_A.
    """
    rayleigh['L0_loops_H'] = rayleigh['L0_H']
    rayleigh['nu_loops_H_per_A'] = rayleigh['nu_H_per_A']

    def period_mismatch(l0):
        return small_signal_period(g_static, R0 + onset['R'], c1, C2, l0, rayleigh['r0_ohm']) - onset['T']
    lo, hi = period_mismatch(10e-3), period_mismatch(40e-3)
    if np.isfinite(lo) and np.isfinite(hi) and lo * hi < 0:
        rayleigh['L0_H'] = float(brentq(period_mismatch, 10e-3, 40e-3, xtol=1e-7))
        w = bin_weights(I)
        rayleigh['nu_H_per_A'] = float(2 * np.sum(w * I * (L_eff - rayleigh['L0_H'])) / np.sum(w * I ** 2))


def dynamic_deviation(large, g_static):
    """
    On the highest-R large cycle: the peak of i_NR minus the static law inside
    |v1| < 3 V, with node 1 at the record's own loop-integral capacitance, and the peak slew.
    """
    if not large:
        return {}
    r = max(large, key=lambda r: r['rpot'])
    v1, v2, dv1, dv2, T = r['cycle']
    inr = (v2 - v1) / r['Rt'] - r['C1_loop'] * dv1
    d = inr - g_static(v1)
    m = np.abs(v1) < 3
    return dict(record=r['label'], rpot=r['rpot'], peak_mA=float(np.max(np.abs(d[m])) * 1e3),
                max_dv1dt_V_per_us=float(np.max(np.abs(dv1)) * 1e-6))


def report_text(n_rows, per, small, tiny, large, gain, c1, law, fits, c2_fit, c2_origin, c2_source, rayleigh,
                onset, T_model, R_hopf_model, loops, pwl, kennedy, dev):
    """The text of identified.txt."""
    core, lcf = fits['core'], fits['large_cycle']
    lines = [f'identify.py: {n_rows} records, {len(per)} periodic ({len(small)} period-1, '
             f'{len(tiny)} near-origin, {len(large)} large cycle)',
             f'divider consistency: A + B = {gain.mean():.4f} (max deviation {np.max(np.abs(gain - 1)):.4f}); '
             f'does not establish relative gain accuracy',
             f'C1 from the loop integral: {c1["small"] * 1e9:.3f} +- {c1["small_se"] * 1e9:.3f} nF at small amplitude '
             f'({len(c1["small_records"])} cycles with A1 < {c1["fit_max_V"] * 1e3:.0f} mV), '
             f'{c1["shoulder"] * 1e9:.2f} nF over every cycle on the shoulder, '
             f'{c1["period1"] * 1e9:.2f} nF on the period-1 cycle, '
             f'{c1["large_cycle"] * 1e9:.2f} nF on the large cycle',
             f'C1 from the KCL fits: {core["C1_apparent_nF"]:.2f} nF (core records), '
             f'{lcf["C1_apparent_nF"]:.2f} nF (large cycle)',
             f'node 1 law over {law["records"]} records ({law["periodic_records"]} by the cycle estimator): capacitor '
             f'C_s = {law["C_s_nF"]:.3f} +- {law["C_s_se_nF"]:.3f} nF + {law["a0_nF_per_V"]:.3f} +- '
             f'{law["a0_se_nF_per_V"]:.3f} nF/V min(x, {law["knot_V"]} V) + {law["a1_nF_per_V"]:.4f} +- '
             f'{law["a1_se_nF_per_V"]:.4f} nF/V max(x - {law["knot_V"]} V, 0), x the swing since the last turning '
             f'point; stage B linear adds {law["stage_B_share_nF"]:.3f} +- {law["stage_B_share_se_nF"]:.3f} nF '
             f'(lag {law["stage_B_lag_us"]:.2f} +- {law["stage_B_lag_se_us"]:.2f} us), stage A '
             f'{law["stage_A_share_nF"]:.3f} nF (datasheet); residual {law["residual_rms_nF"]:.3f} nF rms, by group '
             + ', '.join(f'{g} {v:+.3f}' for g, v in law['residual_by_group_nF'].items())
             + f'; fitted on the periodic records alone: a1 {law["periodic_only"]["a1_nF_per_V"]:.4f} nF/V, stage B '
             f'{law["periodic_only"]["stage_B_share_nF"]:.3f} nF, predicting the others to '
             f'{law["periodic_only"]["chaotic_rms_nF"]:.3f} nF rms',
             f'C2 = {c2_origin["C2_nF"]:.1f} nF from the {c2_source} (L = {c2_origin["L_mH"]:.1f} mH, '
             f'r = {c2_origin["r_ohm"]:.1f} ohm there, residual {c2_origin["relative_residual"]:.3f}); '
             f'the period-1 cycles alone give C2 = {c2_fit["C2_nF"]:.1f} nF, L = {c2_fit["L_mH"]:.1f} mH, '
             f'r = {c2_fit["r_ohm"]:.1f} ohm ({c2_fit["harmonics"]} harmonics, '
             f'residual {c2_fit["relative_residual"]:.3f})',
             f'inductor, Rayleigh law: L0 = {rayleigh["L0_H"] * 1e3:.2f} mH (from the onset period '
             f'{onset["T"] * 1e6:.2f} us at {onset["R"]:.1f} ohm, {onset["label"]};'
             f' the flux loops extrapolate to {rayleigh["L0_loops_H"] * 1e3:.2f} mH), '
             f'nu = {rayleigh["nu_H_per_A"]:.2f} H/A'
             f' (loops alone {rayleigh["nu_loops_H_per_A"]:.2f}), r0 = {rayleigh["r0_ohm"]:.2f} ohm, '
             f'rho = {rayleigh["rho_ohm_per_A"]:.0f} ohm/A',
             f'small-signal check: the linearised circuit oscillates at {T_model * 1e6:.2f} us at '
             f'{onset["R"]:.1f} ohm and its Hopf point is at'
             f' {R_hopf_model - R0 if R_hopf_model else float("nan"):.1f} ohm (measured onset: onset.py)',
             f'{"record":22s} {"Rpot":>7} {"T us":>6} {"I mA":>6} {"L_eff mH":>8} {"r_eff":>6}']
    for t in sorted(loops, key=lambda t: t.I):
        lines.append(f'{t.label:22s} {t.rpot:7.1f} {t.T_us:6.1f} {t.I * 1e3:6.2f} {t.L_eff * 1e3:8.2f} {t.r_eff:6.1f}')
    lines += [f'diode PWL: breakpoints {np.round(pwl["bp_V"], 3).tolist()} V, '
              f'slopes {np.round(np.array(pwl["G_S"]) * 1e3, 4).tolist()} mS',
              f'  inner three lines: KCL fit over {core["records"]} core records '
              f'(residual {core["relative_residual"]:.4f}); '
              f'saturation lines: plateaus of {pwl["saturation_records"]} large-cycle records; '
              f'KCL fit over the large cycle alone: residual {lcf["relative_residual"]:.4f}, '
              f'apparent C1 {lcf["C1_apparent_nF"]:.2f} nF',
              f'Kennedy diode: RA = {kennedy["RA_ohm"]:.0f} ohm, AA = {kennedy["AA"]:.3f}, '
              f'RB = {kennedy["RB_ohm"]:.0f} ohm, AB = {kennedy["AB"]:.2f}, '
              f'rails stage A {kennedy["VpA_V"]:+.2f}/{kennedy["VnA_V"]:+.2f} V, '
              f'stage B {kennedy["VpB_V"]:+.2f}/{kennedy["VnB_V"]:+.2f} V']
    if dev:
        lines.append(f'large-cycle diode deviation ({dev["record"]}): up to {dev["peak_mA"]:.2f} mA inside '
                     f'|v1| < 3 V at {dev["max_dv1dt_V_per_us"]:.2f} V/us')
    return '\n'.join(lines)


def plot_inductor(per, loops, rayleigh, path):
    """Flux loops at four amplitudes, and L_eff and r_eff against amplitude with the Rayleigh lines."""
    I = np.array([t.I for t in loops])
    Le = np.array([t.L_eff for t in loops])
    re = np.array([t.r_eff for t in loops])
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    ax = axes[0]
    show = sorted(per, key=lambda r: r['I_amp'])
    pick = [show[0], show[len(show) // 3], show[2 * len(show) // 3], show[-1]] if len(show) >= 4 else show
    for r in pick:
        ax.plot(r['iL'] * 1e3, r['phi'] * 1e3, lw=1, label=f'{r["label"]} ({r["rpot"]:.0f} ohm), I = {r["I_amp"]*1e3:.1f} mA')
    ax.set_xlabel('i_L (mA)')
    ax.set_ylabel('flux = int v2 dt (mWb)')
    ax.set_title('inductor: flux against current, one averaged cycle')
    ax.legend(fontsize=7)
    ax.grid(alpha=.3)
    ax = axes[1]
    ax.plot(I * 1e3, Le * 1e3, 'o', ms=4, label='secant inductance of the loop')
    xx = np.linspace(0, I.max() * 1e3, 50)
    ax.plot(xx, (rayleigh['L0_H'] + 0.5 * rayleigh['nu_H_per_A'] * xx * 1e-3) * 1e3, '-',
            label=f'L0 + nu I/2, L0 = {rayleigh["L0_H"]*1e3:.1f} mH')
    ax.set_xlabel('current amplitude (mA)')
    ax.set_ylabel('L_eff (mH)')
    ax.legend(fontsize=8)
    ax.grid(alpha=.3)
    ax.set_title('effective inductance against amplitude')
    ax = axes[2]
    ax.plot(I * 1e3, re, 'o', ms=4, label='loss resistance of the loop')
    ax.plot(xx, rayleigh['r0_ohm'] + (8 / (3 * np.pi)) * rayleigh['rho_ohm_per_A'] * xx * 1e-3, '-',
            label=f'r0 + 0.85 rho I, r0 = {rayleigh["r0_ohm"]:.1f} ohm')
    ax.set_xlabel('current amplitude (mA)')
    ax.set_ylabel('r_eff (ohm)')
    ax.legend(fontsize=8)
    ax.grid(alpha=.3)
    ax.set_title('loss resistance against amplitude')
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_diode(small, large, fits, pwl, path):
    """
    i_NR from a period-1 and a large-cycle record (each at its own loop-integral
    C1) over the V-I trace, with both spline fits and the breakpoints.
    """
    fig, ax = plt.subplots(figsize=(9, 6))
    vi = os.path.join(HERE, '..', 'trace1.csv')
    if os.path.exists(vi):
        d = np.genfromtxt(vi, delimiter=',', skip_header=1)
        i_m = (d[:, 1] - d[:, 2]) / 216.0
        off = np.nanmedian(i_m[np.abs(d[:, 2]) < 0.15])
        ax.plot(d[:, 2], (i_m - off) * 1e3, '.', ms=1, alpha=0.1, color='steelblue',
                label='V-I trace (trace1.csv, offset removed)')
    examples = [(small[len(small) // 2] if small else None, 'C1'),
                (max(large, key=lambda r: r['rpot']) if large else None, 'C3')]
    for r, c in examples:
        if r is None:
            continue
        v1, v2, dv1, dv2, T = r['cycle']
        inr = (v2 - v1) / r['Rt'] - r['C1_loop'] * dv1
        ax.plot(v1, inr * 1e3, '.', ms=1.5, color=c, label=f'i_NR from {r["label"]} ({r["rpot"]:.0f} ohm)')
    ax.plot(fits['core']['curve_v'], fits['core']['curve_i'] * 1e3, 'k-', lw=1.2, label='KCL spline fit, core records')
    ax.plot(fits['large_cycle']['curve_v'], fits['large_cycle']['curve_i'] * 1e3, 'k--', lw=1,
            label='KCL spline fit, large cycle')
    for b in pwl['bp_V']:
        ax.axvline(b, color='0.6', lw=0.6, ls=':')
    ax.set_xlim(-8, 7.4)
    ax.set_ylim(-4, 4)
    ax.set_xlabel('v1 (V)')
    ax.set_ylabel('current into N_R (mA)')
    ax.set_title('the nonlinear element as the circuit sees it')
    ax.legend(fontsize=8)
    ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_c1(law, path):
    """Node 1's loop-integral capacitance per record against Rpot with the law, and the capacitor's share against the swing."""
    rec = law['per_record']
    per = rec['periodic']
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    ax = axes[0]
    ax.plot(rec['rpot'][per], rec['C'][per] * 1e9, 'o', ms=3, label='periodic records (cycle estimator)')
    ax.plot(rec['rpot'][~per], rec['C'][~per] * 1e9, 'o', ms=3, label='other records (window estimator)')
    ax.plot(rec['rpot'], rec['law'] * 1e9, 'k.', ms=1.5, label='the law')
    ax.set_xlabel('Rpot (ohm)')
    ax.set_ylabel('loop-integral capacitance of node 1 (nF)')
    ax.legend(fontsize=8)
    ax = axes[1]
    ax.plot(rec['swing'][per], rec['capacitor'][per] * 1e9, 'o', ms=3)
    ax.plot(rec['swing'][~per], rec['capacitor'][~per] * 1e9, 'o', ms=3)
    x = np.linspace(0, rec['swing'].max(), 200)
    k = law['knot_V']
    ax.plot(x, law['C_s_nF'] + law['a0_nF_per_V'] * np.minimum(x, k) + law['a1_nF_per_V'] * np.maximum(x - k, 0),
            'k-', lw=1, label='C_s + a0 min(x, k) + a1 max(x - k, 0)')
    ax.set_xlabel('weighted swing of v1 since its last turning point (V)')
    ax.set_ylabel('capacitor share (op-amp lags removed, nF)')
    ax.legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--quick', action='store_true', help='every 4th record')
    p.add_argument('--out', default=HERE)
    a = p.parse_args()
    step = 4 if a.quick else 1

    # 1. every accepted record, the periodic ones averaged over their cycles
    rows = analyse_sweeps(step)
    per = [r for r in rows if r['periodic']]
    small = [r for r in per if r['v1_max'] < 5.5 and r['v1_min'] < -1.5]      # period-1 cycles
    tiny = [r for r in per if r['v1_max'] < 5.5 and r['v1_min'] >= -1.5]     # the orbit near the origin
    large = [r for r in per if r['v1_max'] > 5.5]                              # the large outer cycle
    gain = np.array([r['divider_A_plus_B'] for r in rows])

    # 2. C1 from the loop integral
    top, _, osc_min = onset_table()
    c1 = loop_c1_values(small, large, top)

    # 3. the element: integrated KCL for the inner lines, the large cycle's plateaus for saturation
    fits = {}
    for name, records in kcl_records(step).items():
        fits[name] = kcl_fit(records)
        print(f'KCL fit {name}: {len(records)} records, apparent C1 {fits[name]["C1_apparent_nF"]:.3f} nF, '
              f'residual {fits[name]["relative_residual"]:.4f}', flush=True)
    pwl, kennedy = five_segment_law(fits['core']['pwl'], large, fits['core']['C1_apparent_nF'] * 1e-9)
    vk, ik = pwl_knots(pwl['bp_V'], pwl['G_S'])

    def g_static(v):
        return np.interp(v, vk, ik)

    # 3b. node 1's capacitance law over every record, and stage B's effective lag
    law = node1_law(rows, c1['small_records'], pwl, kennedy)

    # 4. C2 from the tank admittance, the inductor loops, and the Rayleigh law closed on the onset
    c2_fit, c2_origin, c2_source = tank_fits(small, tiny)
    C2 = c2_origin['C2_nF'] * 1e-9
    loops = inductor_loops(per, C2)
    I = np.array([t.I for t in loops])
    Le = np.array([t.L_eff for t in loops])
    re = np.array([t.r_eff for t in loops])
    ok = I > 0.6e-3            # below this the flux loop is a few quantisation steps wide
    rayleigh = rayleigh_fit(I[ok], Le[ok], re[ok])
    onset = onset_record(top, osc_min)
    calibrate_l0(rayleigh, I[ok], Le[ok], g_static, c1['small'], C2, onset)
    T_model = small_signal_period(g_static, R0 + onset['R'], c1['small'], C2, rayleigh['L0_H'], rayleigh['r0_ohm'])
    R_hopf_model = small_signal_hopf(g_static, c1['small'], C2, rayleigh['L0_H'], rayleigh['r0_ohm'])

    # 5. the dynamic deviation of the element in the large cycle
    dev = dynamic_deviation(large, g_static)

    result = dict(
        records_analysed=len(rows), periodic_records=len(per),
        divider_gain_check=dict(mean_A_plus_B=float(gain.mean()), max_deviation=float(np.max(np.abs(gain - 1)))),
        C1_nF=dict(small_signal_loop=c1['small'] * 1e9, small_signal_se=c1['small_se'] * 1e9,
                   small_signal_records=c1['small_records'], shoulder_mean_loop=c1['shoulder'] * 1e9,
                   period1_loop=c1['period1'] * 1e9, large_cycle_loop=c1['large_cycle'] * 1e9,
                   kcl_core=fits['core']['C1_apparent_nF'], kcl_large_cycle=fits['large_cycle']['C1_apparent_nF'],
                   law={k: v for k, v in law.items() if k != 'per_record'}),
        C2_nF=c2_origin['C2_nF'],
        tank_admittance_fit=c2_fit,
        C2_source=c2_source,
        tank_admittance_near_origin=c2_origin,
        small_signal_check=dict(onset_record=onset['label'], onset_period_us=onset['T'] * 1e6, onset_rpot=onset['R'],
                                model_period_us=T_model * 1e6,
                                hopf_rpot_model=None if R_hopf_model is None else R_hopf_model - R0),
        inductor=dict(rayleigh=rayleigh, per_record=[dict(label=t.label, rpot=t.rpot, T_us=t.T_us, I_amp_mA=t.I * 1e3,
                                                          L_eff_mH=t.L_eff * 1e3, r_eff_ohm=t.r_eff) for t in loops]),
        diode=dict(pwl=pwl, kennedy=kennedy,
                   kcl_core=dict(records=fits['core']['records'], relative_residual=fits['core']['relative_residual'],
                                 **fits['core']['pwl']),
                   kcl_large_cycle=dict(records=fits['large_cycle']['records'],
                                        relative_residual=fits['large_cycle']['relative_residual'],
                                        **fits['large_cycle']['pwl']),
                   large_cycle_dynamic_deviation=dev),
        notes=[
            'C1 (small_signal_loop) is the capacitance node 1 shows on the onset cycles. C1.law describes '
            'its growth over every record: the capacitor rises with the swing of v1 since its last turning '
            'point (the Rayleigh behaviour of a ferroelectric ceramic dielectric), and stage B of the diode '
            'adds stage_B_share_nF while it is linear (inside the inner breakpoints), an effective lag of '
            'stage_B_lag_us; stage A adds stage_A_share_nF at its datasheet speed (see RESULTS.md).',
            'inductor: L_eff and r_eff are the secant inductance and the loss resistance of the measured '
            'flux-current loop, about the mean current of the record; the Rayleigh law L = L0 + nu|d|, core '
            'loss rho|d|d on top of the winding resistance r0, with d the excursion of iL from its mean, '
            'reproduces their growth with amplitude.',
            'diode: pwl is the five-segment law as the circuit sees it (inner three lines from the core '
            'records, saturation lines from the large cycle); kennedy maps it onto the two-op-amp realisation.'])
    with open(os.path.join(a.out, 'identified.json'), 'w') as fh:
        json.dump(result, fh, indent=2)
    txt = report_text(len(rows), per, small, tiny, large, gain, c1, law, fits, c2_fit, c2_origin, c2_source,
                      rayleigh, onset, T_model, R_hopf_model, loops, pwl, kennedy, dev)
    print(txt)
    with open(os.path.join(a.out, 'identified.txt'), 'w') as fh:
        fh.write(txt + '\n')

    plot_inductor(per, loops, rayleigh, os.path.join(a.out, 'identified_inductor.png'))
    plot_diode(small, large, fits, pwl, os.path.join(a.out, 'identified_diode.png'))
    plot_c1(law, os.path.join(a.out, 'identified_c1.png'))
    print('written identified.json, identified.txt and the three figures to', a.out)


if __name__ == '__main__':
    main()
