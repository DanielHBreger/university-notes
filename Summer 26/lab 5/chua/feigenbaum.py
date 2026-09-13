"""
Feigenbaum delta of the simulated Chua oscillator: the period-doubling
points R_1, R_2, ... of the negative single-scroll branch, for the ideal
model (simulate.py --model ideal) or the identified bench model (--model bench).

Method
------
The doubling points are NOT found by counting attractor levels: near a
bifurcation the new levels are closer together than any clustering
tolerance. Each R_k is found from the normal form of the doubling instead.
For a settled 2^k orbit the sequence of v1 maxima m_i repeats every 2^k, and
the splitting created by the MOST RECENT doubling,

    s_k = < |m_i - m_(i + 2^(k-1))| >_i ,

is zero on the 2^(k-1) orbit above R_k and grows as sqrt(R_k - R) below it,
so s_k^2 is linear in R with its zero at R_k. The 2^k window is located by
bisection on the period read from the lag distances D(p) = <|m_i - m_(i+p)|>
(the same test cascade_periods.py applies to the bench records), s_k is then
sampled at a dozen resistances inside the window and the root of a quadratic
fit of s_k^2 is R_k. The orbit at every resistance is settled from the
continuation state of the previous one, as in the sweeps.

Results are printed and written to feigenbaum_<model>.txt.

Usage:
    python feigenbaum.py                    # ideal model (C1 = 11.5 nF), k = 1..4
    python feigenbaum.py --model static     # identified circuit, constant components
    python feigenbaum.py --model bench      # identified circuit
    python feigenbaum.py --kmax 5 --settle 3
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import simulate as s
import integration as kern

HERE = os.path.dirname(os.path.abspath(__file__))
LAGS = (1, 2, 4, 8, 16, 32, 64)


class Model:
    """Settle-and-collect at one resistance, carrying the state along."""

    def __init__(self, kind, dt=None, tau_b=None, slew=None):
        self.kind = kind
        if kind == 'bench':
            self.P = s.bench_params(tau_b=tau_b or s.TAU_B, slew=slew or s.SLEW)
            self.g = s.bench_static_g(self.P)
            self.rL = self.P[4]
            self.dt = dt or s.BENCH_DT
        elif kind == 'static':
            P = s.bench_params(tau_b=tau_b or s.TAU_B, slew=slew or s.SLEW)
            c1a, c2b, l0, r0b = s.bench_small_signal(P)
            s.set_components(c1a, c2b, l0)
            self.g = s.bench_static_g(P)
            self.rL = r0b
            self.dt = dt or 0.5e-6
        else:
            self.g = s.make_g()
            self.rL = 0.0
            self.dt = dt or 0.5e-6
        self.y = None

    def start(self, rpot):
        g = self.g
        eqs = [v for v in s.equilibria(g, s.R0 + rpot, self.rL) if v < 0]
        v = eqs[0] if eqs else 0.0
        iL = v / (s.R0 + rpot + self.rL)
        if self.kind == 'bench':
            self.y = kern.bench_state(v + 0.05, self.rL * iL, iL, self.P)
        else:
            self.y = np.array([v + 0.05, self.rL * iL, iL])

    def hold(self, rpot, t_settle, t_collect):
        n1, n2 = int(round(t_settle / self.dt)), int(round(t_collect / self.dt))
        if self.kind == 'bench':
            self.y, mx = kern.bench_hold(self.y, s.R0 + rpot, self.dt, n1, n2, self.P)
        else:
            vk, ik = self.g.knots
            self.y, mx = kern.hold(self.y, s.R0 + rpot, self.rL, s.C1, s.C2, s.L, vk, ik, self.dt, n1, n2)
        return np.asarray(mx, float)


def period_of(m, factor=2.0, drop=0.25):
    """Period 2^j of a settled maxima sequence from its lag distances (0 if none)."""
    t = np.asarray(m, float)[int(len(m) * drop):]
    if len(t) < 70:
        return 0
    D = {p: float(np.mean(np.abs(t[:-p] - t[p:]))) for p in LAGS if len(t) > p}
    floor = max(min(D.values()), 2e-4)     # a settled orbit repeats to well under a millivolt
    if D[1] < 0.02:                 # a period-1 orbit: no structure at all
        return 1
    if floor > 0.25 * D[1]:         # no lag brings the sequence back onto itself: aperiodic
        return 0
    for p in LAGS:
        if p in D and D[p] < factor * floor:
            return p
    return 0


def splitting(m, k, ntail=256):
    m = np.asarray(m, float)[-ntail:]
    off = 1 << (k - 1)
    return float(np.mean(np.abs(m[:-off] - m[off:]))) if len(m) > off else np.nan


def find_window(model, k, r_hi, r_lo, t_settle, t_collect, tol=0.05):
    """
    (R_top, R_bottom) of the 2^k regime by bisection on the period label
    between r_hi (period 2^(k-1)) and r_lo (period >= 2^(k+1) or chaotic).
    """
    P = 1 << k
    a, b = r_hi, r_lo
    # top of the window: where 2^(k-1) turns into 2^k
    lo, hi = b, a
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        model.start(mid)
        model.hold(mid, t_settle, 0.001)
        p = period_of(model.hold(mid, t_settle, t_collect))
        if p and p < P:
            hi = mid
        else:
            lo = mid
    top = lo
    # bottom of the window: where 2^k turns into 2^(k+1) or worse
    lo, hi = b, top
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        model.start(mid)
        model.hold(mid, t_settle, 0.001)
        p = period_of(model.hold(mid, t_settle, t_collect))
        if p == P:
            hi = mid
        else:
            lo = mid
    return top, hi


def doubling_point(model, k, top, bottom, t_settle, t_collect, npts=12):
    """R_k from the sqrt law of the splitting sampled inside the window."""
    rs = np.linspace(bottom + 0.05 * (top - bottom), top - 0.1 * (top - bottom), npts)[::-1]
    model.start(rs[0] + 0.02 * (top - bottom))
    model.hold(rs[0], t_settle, 0.001)
    sp = []
    for r in rs:
        m = model.hold(r, t_settle, t_collect)
        sp.append(splitting(m, k))
    sp = np.array(sp)
    ok = np.isfinite(sp) & (sp > 0)
    coef = np.polyfit(rs[ok], sp[ok] ** 2, 2)
    roots = np.roots(coef)
    roots = np.array([r.real for r in roots if abs(r.imag) < 1e-9 and r.real > rs.max() - 1e-9])
    return (float(roots.min()) if len(roots) else np.nan), rs, sp


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--model', choices=['ideal', 'static', 'bench'], default='ideal')
    p.add_argument('--kmax', type=int, default=4)
    p.add_argument('--settle', type=float, default=2.0, help='seconds settled at each resistance in the fits')
    p.add_argument('--collect', type=float, default=0.3)
    p.add_argument('--range', type=float, nargs=2, default=[880.0, 700.0], metavar=('RHI', 'RLO'),
                   help='period-1 resistance above the cascade and a chaotic one below it')
    p.add_argument('--tau-b', type=float, default=None, help='bench: stage-B lag (us)')
    p.add_argument('--slew', type=float, default=None, help='bench: slew rate (V/us)')
    a = p.parse_args()
    model = Model(a.model, tau_b=a.tau_b and a.tau_b * 1e-6, slew=a.slew and a.slew * 1e6)
    lines = [f'feigenbaum.py --model {a.model}' + (f' (C1 = {s.C1 * 1e9:g} nF)' if a.model == 'ideal' else '')]
    R = {}
    r_hi, r_lo = a.range
    for k in range(1, a.kmax + 1):
        top, bottom = find_window(model, k, r_hi, r_lo, 0.3 * a.settle, a.collect)
        if bottom >= top - 0.2:
            lines.append(f'k={k}: window of period {1 << k} not resolved ({top:.2f} .. {bottom:.2f} ohm)')
            print(lines[-1], flush=True)
            break
        Rk, rs, sp = doubling_point(model, k, top, bottom, a.settle, a.collect)
        R[k] = Rk
        lines.append(f'k={k}  period {1 << (k - 1):2d} -> {1 << k:<2d}  R_{k} = {Rk:9.3f} ohm   '
                     f'(window {top:.2f} .. {bottom:.2f}, s from {sp[0]:.4f} to {sp[-1]:.4f} V)')
        print(lines[-1], flush=True)
        r_hi = top
    lines.append(f'\n{"n":>2} {"g_n = R_n - R_(n+1)":>21} {"delta_n":>10}')
    for n in sorted(R):
        if n + 1 not in R:
            continue
        g = R[n] - R[n + 1]
        d = f'{g / (R[n + 1] - R[n + 2]):10.3f}' if n + 2 in R else ''
        lines.append(f'{n:2d} {g:21.4f} {d}')
    lines.append('universal delta = 4.669202')
    if len(R) >= 2:
        k = max(R)
        lines.append(f'accumulation point R_inf ~ {R[k] - (R[k - 1] - R[k]) / (4.669202 - 1):.3f} ohm')
    txt = '\n'.join(lines)
    print(txt)
    with open(os.path.join(HERE, f'feigenbaum_{a.model}.txt'), 'w') as fh:
        fh.write(txt + '\n')


if __name__ == '__main__':
    main()
