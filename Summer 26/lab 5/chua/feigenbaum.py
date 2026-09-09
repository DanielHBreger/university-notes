"""
Feigenbaum delta of the measured Chua oscillator, from the calibrated model
of simulate.py.

Method
------
The period-doubling points R_n are NOT found by counting attractor levels --
near a bifurcation the new levels are closer together than any sensible
clustering tolerance, and the count is contaminated by the (very long)
transient. Instead each R_n is found from the normal form of the doubling.

For a settled 2^k orbit the sequence of v_C1 maxima m_i is periodic with
period 2^k, and the splitting created by the MOST RECENT doubling is

    s_k = < |m_i - m_(i + 2^(k-1))| >_i

which is identically zero on the 2^(k-1) orbit above R_k and grows as
sqrt(R_k - R) below it. So s_k^2 is linear in R with its zero at R_k:
sample s_k at a dozen values of R inside the 2^k window, fit, and read off
the root. No tolerance, no level counting, and the sqrt law means the
extrapolation is short even from well inside the window.

Integration uses a cheap long transient at a coarse step followed by an
accurate settle and collection at dt (the attractor is attracting, so the
coarse phase only has to land near it). Halving dt to 0.25 us moves R_3 by
0.0002 ohm and R_4 by 0.0024 ohm, i.e. <= 0.2 % of the smallest gap used.

Result (negative single-scroll branch, C1 = 11.5 nF, rL = 0):

    n   period      R_n (ohm)    g_n = R_n - R_(n+1)
    1   1  ->  2     798.758     25.921
    2   2  ->  4     772.837      6.048
    3   4  ->  8     766.789      1.316
    4   8  -> 16     765.474      0.282
    5  16  -> 32     765.192

    delta_1 = 4.286    delta_2 = 4.596    delta_3 = 4.666
    universal value    4.669202
    accumulation point R_inf = 765.115 ohm

Usage:
    python feigenbaum.py            # recompute every R_n (about 25 min)
    python feigenbaum.py --quick    # coarser, ~5 min, R_1..R_3 only
"""
import argparse
import sys
import numpy as np

sys.path.insert(0, '.')
import simulate as s

RL = 0.0
DT = 0.5e-6
G = s.make_g()

# (window low, window high, step) inside the 2^k regime, per doubling k
WINDOWS = {1: (790.00, 798.00, 0.50),
           2: (767.20, 772.40, 0.40),
           3: (765.60, 766.70, 0.10),
           4: (765.24, 765.44, 0.02),
           5: (765.142, 765.186, 0.004)}
SETTLE = {1: (4.0, 0.5, 0.35), 2: (4.0, 0.5, 0.35), 3: (5.0, 0.6, 0.40),
          4: (8.0, 1.0, 0.50), 5: (10.0, 1.5, 0.60)}


def _run(y, rp, t, dt, collect=False):
    Rt = s.R0 + np.asarray(rp, float)
    nsteps = int(round(t / dt))
    if not collect:
        for _ in range(nsteps):
            y = s.rk4_step(y, dt, Rt, RL, G)
        return y
    mx = [[] for _ in range(len(Rt))]
    p2 = y[0].copy()
    p1 = y[0].copy()
    for _ in range(nsteps):
        y = s.rk4_step(y, dt, Rt, RL, G)
        cur = y[0]
        for i in np.flatnonzero((p1 > p2) & (p1 >= cur)):
            mx[i].append(p1[i])
        p2, p1 = p1, cur.copy()
    return mx


def settled_maxima(rp, t_coarse, t_fine, t_take, dt=DT, dt_coarse=2e-6, side=-1):
    rp = np.atleast_1d(np.asarray(rp, float))
    y, _ = s.sweep_starts(G, rp, RL, side=side)
    y = _run(y, rp, t_coarse, dt_coarse)
    y = _run(y, rp, t_fine, dt)
    return _run(y, rp, t_take, dt, collect=True)


def splitting(m, k, ntail=256):
    """Newest period-doubling splitting of a settled 2^k orbit."""
    m = np.asarray(m, float)[-ntail:]
    off = 1 << (k - 1)
    return float(np.mean(np.abs(m[:-off] - m[off:]))) if len(m) > off else np.nan


def bifurcation_point(rs, s2, deg=2):
    """Root of a polynomial fit of s^2 against R, taken just above the window."""
    roots = np.roots(np.polyfit(rs, s2, deg))
    roots = np.array([r.real for r in roots if abs(r.imag) < 1e-9])
    cand = roots[roots > np.max(rs) - 1e-9]
    return float(cand.min()) if len(cand) else np.nan


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--quick', action='store_true')
    p.add_argument('--kmax', type=int, default=5)
    a = p.parse_args()
    kmax = min(3, a.kmax) if a.quick else a.kmax

    R = {}
    for k in range(1, kmax + 1):
        lo, hi, step = WINDOWS[k]
        tc, tf, tt = SETTLE[k]
        if a.quick:
            tc, tf, tt = tc / 2, tf / 2, tt
        rs = np.arange(lo, hi + step / 2, step)
        M = settled_maxima(rs, tc, tf, tt)
        sp = np.array([splitting(m, k) for m in M])
        R[k] = bifurcation_point(rs, sp ** 2, 2)
        print(f'k={k}  period {2**(k-1):2d} -> {2**k:<2d}  '
              f'R_{k} = {R[k]:9.4f} ohm   '
              f'(s from {sp[0]:.5f} down to {sp[-1]:.5f} V over '
              f'{len(rs)} points)')

    print(f'\n{"n":>2} {"g_n = R_n - R_(n+1)":>21} {"delta_n":>10}')
    for n in sorted(R):
        if n + 1 not in R:
            continue
        g = R[n] - R[n + 1]
        d = (f'{g / (R[n+1] - R[n+2]):10.3f}' if n + 2 in R else '')
        print(f'{n:2d} {g:21.4f} {d}')
    if len(R) >= 3:
        print(f'\nuniversal delta = 4.669202')
    if len(R) >= 2:
        k = max(R)
        print(f'accumulation point R_inf ~ '
              f'{R[k] - (R[k-1] - R[k]) / (4.669202 - 1):.3f} ohm')


if __name__ == '__main__':
    main()
