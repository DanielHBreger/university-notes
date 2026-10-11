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
    python feigenbaum.py                    # the plan's model (nominal parts, Table 1 element), k = 1..4
    python feigenbaum.py --model static     # identified circuit, constant components
    python feigenbaum.py --model bench      # identified circuit
    python feigenbaum.py --kmax 5 --settle 3
"""
import argparse
import os

import numpy as np

import simulate as sim
import integration as kern
from cascade_periods import LAGS

HERE = os.path.dirname(os.path.abspath(__file__))
DELTA = 4.669202


class Model:
    """One of the three models, settled and sampled at one resistance at a time, carrying its state along."""

    def __init__(self, kind, dt=None, tau_b=None, slew=None):
        self.kind = kind
        self.P = None                                   # bench kernel parameters
        if kind == 'ideal':
            self.g, self.circuit = sim.table1_g(), sim.NOMINAL
        else:
            P = sim.bench_params(tau_b=tau_b, slew=slew or sim.SLEW)
            self.g, self.circuit = sim.bench_static_g(P), sim.bench_small_signal(P)
            if kind == 'bench':
                self.P = P
        self.dt = dt or (sim.BENCH_DT if kind == 'bench' else sim.IDEAL_DT)
        self.y = None

    def start(self, rpot):
        """Restart on the negative outer equilibrium at rpot, nudged as in the forward sweep."""
        fwd, _ = sim.sweep_starts(self.g, self.circuit, [rpot])
        self.y = fwd[:, 0] if self.P is None else kern.bench_state(*fwd[:, 0], self.P)

    def hold(self, rpot, t_settle, t_collect):
        """Settle at rpot, then return the maxima of v1 over t_collect."""
        n_settle, n_collect = int(round(t_settle / self.dt)), int(round(t_collect / self.dt))
        Rt = self.circuit.rt(rpot)
        if self.P is not None:
            self.y, mx = kern.bench_hold(self.y, Rt, self.dt, n_settle, n_collect, self.P)
        else:
            c = self.circuit
            vk, ik = self.g.knots
            self.y, mx = kern.hold(self.y, Rt, c.rl, c.c1, c.c2, c.l, vk, ik, self.dt, n_settle, n_collect)
        return np.asarray(mx, float)

    def period_at(self, rpot, t_settle, t_collect):
        """Period label of the orbit settled at rpot from a fresh start."""
        self.start(rpot)
        self.hold(rpot, t_settle, 0.001)
        return settled_period(self.hold(rpot, t_settle, t_collect))


def settled_period(m, factor=2.0, drop=0.25):
    """Period 2^j of a settled maxima sequence from its lag distances (0 if none)."""
    t = np.asarray(m, float)[int(len(m) * drop):]
    if len(t) < 70:
        return 0
    D = {p: float(np.mean(np.abs(t[:-p] - t[p:]))) for p in LAGS if len(t) > p}
    floor = max(min(D.values()), 2e-4)     # a settled orbit repeats to well under a millivolt
    if D[1] < 0.02:                         # a period-1 orbit: no structure at all
        return 1
    if floor > 0.25 * D[1]:                 # no lag brings the sequence back onto itself: aperiodic
        return 0
    for p in LAGS:
        if p in D and D[p] < factor * floor:
            return p
    return 0


def splitting(m, k, ntail=256):
    """s_k: the mean distance between maxima 2^(k-1) apart over the last ntail maxima."""
    m = np.asarray(m, float)[-ntail:]
    off = 1 << (k - 1)
    return float(np.mean(np.abs(m[:-off] - m[off:]))) if len(m) > off else np.nan


def bisect(lo, hi, in_upper_part, tol):
    """Narrow [lo, hi] to within tol around the resistance above which in_upper_part(r) holds."""
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if in_upper_part(mid):
            hi = mid
        else:
            lo = mid
    return lo, hi


def find_window(model, k, r_hi, r_lo, t_settle, t_collect, tol=0.05):
    """
    (R_top, R_bottom) of the 2^k regime by bisection on the period label
    between r_hi (period 2^(k-1)) and r_lo (period >= 2^(k+1) or chaotic).
    """
    P = 1 << k

    def period(r):
        return model.period_at(r, t_settle, t_collect)
    top, _ = bisect(r_lo, r_hi, lambda r: 0 < period(r) < P, tol)      # where 2^(k-1) turns into 2^k
    _, bottom = bisect(r_lo, top, lambda r: period(r) == P, tol)       # where 2^k turns into 2^(k+1) or worse
    return top, bottom


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
    lines = [f'feigenbaum.py --model {a.model}' + (f' (C1 = {model.circuit.c1 * 1e9:g} nF)' if a.model == 'ideal' else '')]
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
    lines.append(f'universal delta = {DELTA}')
    if len(R) >= 2:
        k = max(R)
        lines.append(f'accumulation point R_inf ~ {R[k] - (R[k - 1] - R[k]) / (DELTA - 1):.3f} ohm')
    txt = '\n'.join(lines)
    print(txt)
    with open(os.path.join(HERE, f'feigenbaum_{a.model}.txt'), 'w') as fh:
        fh.write(txt + '\n')


if __name__ == '__main__':
    main()
