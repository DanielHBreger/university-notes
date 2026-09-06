"""
Calibrate the two Lyapunov estimators on a simulated Chua circuit.

The circuit is integrated in Chua's dimensionless form,

    x' = alpha (y - x - f(x)),   y' = x - y + z,   z' = -beta y,
    f(x) = m1 x + (m0 - m1)(|x + 1| - |x - 1|) / 2,

with alpha = C2/C1 and beta = C2 R^2 / L from this bench's components and
R_total = R0 + Rpot, and the true largest exponent taken from the
variational equation integrated alongside (one tangent vector, renormalised
each step). The simulated V1 and V2 are then quantised the way the scope
does it (40 mV and 4 mV steps), given a little noise, written as a scope
CSV and pushed through EXACTLY the pipeline the measurements use:
lorenz_map.maxima -> lyapunov.lyapunov for the return-map estimate and
rosenstein.rosenstein for the direct one.

This is what backs the statements in lyapunov.py and rosenstein.py about
how far each estimate can be trusted: run it and read the table. The map is
only approximately a one-dimensional function of M_n, and where it folds
(the lobe switch) the projection has a slope the flow does not; the direct
estimate needs the channels smoothed before the quantisation lets it see
the divergence. Try --alpha 9 --beta 14.286 (Matsumoto's double scroll) and
--alpha 8.6 --beta 14.286 (a single scroll) as well as the bench values.

x is in units of the breakpoint voltage (about 1 V here), so the simulated
V1 spans a few volts like the real one. Time is scaled by C2 R so that one
sample is 1 us, as on the bench.

Usage:
    python benchmark_lyapunov.py                     # Rpot 660 ohm, canonical slopes
    python benchmark_lyapunov.py --rpot 720 --m0 -1.14 --m1 -0.71
    python benchmark_lyapunov.py --alpha 9 --beta 14.286   # Matsumoto's double scroll
"""
import argparse
import os
import tempfile

import numpy as np

from lorenz_map import maxima
from lyapunov import lyapunov
from rosenstein import rosenstein

C1, C2, L, R0 = 10e-9, 100e-9, 18e-3, 990.0


def simulate(alpha, beta, m0, m1, n, dtau, seed=0, transient=50000):
    """
    (states, true lambda per unit dimensionless time).

    Fixed-step RK4 on the six-dimensional system (state + tangent vector);
    the tangent vector is renormalised every step and the log of its growth
    accumulated, which is the standard largest-exponent calculation.
    """
    def fx(x):
        return m1 * x + 0.5 * (m0 - m1) * (abs(x + 1) - abs(x - 1))

    def rhs(s):
        x, y, z, vx, vy, vz = s
        fp = m1 if abs(x) > 1 else m0
        return np.array([alpha * (y - x - fx(x)), x - y + z, -beta * y,
                         -alpha * (1 + fp) * vx + alpha * vy, vx - vy + vz,
                         -beta * vy])

    def step(s):
        k1 = rhs(s)
        k2 = rhs(s + 0.5 * dtau * k1)
        k3 = rhs(s + 0.5 * dtau * k2)
        k4 = rhs(s + dtau * k3)
        return s + dtau / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)

    rng = np.random.default_rng(seed)
    s = np.array([0.1 + 0.01 * rng.standard_normal(), 0.0, 0.0, 1.0, 0.0, 0.0])
    for _ in range(transient):
        s = step(s)
        s[3:] /= np.linalg.norm(s[3:])
    out = np.empty((n, 3))
    lnsum = 0.0
    for i in range(n):
        s = step(s)
        nv = np.linalg.norm(s[3:])
        lnsum += np.log(nv)
        s[3:] /= nv
        out[i] = s[:3]
    return out, lnsum / (n * dtau)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--rpot', type=float, default=660.0,
                   help='potentiometer setting, ohm (sets beta)')
    p.add_argument('--alpha', type=float, help='override alpha = C2/C1')
    p.add_argument('--beta', type=float, help='override beta = C2 R^2 / L')
    p.add_argument('--m0', type=float, default=-8 / 7,
                   help='inner slope Ga R (dimensionless)')
    p.add_argument('--m1', type=float, default=-5 / 7,
                   help='outer slope Gb R (dimensionless)')
    p.add_argument('--samples', type=int, default=200000,
                   help='record length in 1 us samples')
    p.add_argument('--lsb', type=float, nargs=2, default=(0.04, 0.004),
                   metavar=('V1', 'V2'), help='quantisation steps, V')
    p.add_argument('--noise', type=float, default=0.005,
                   help='rms noise added before quantisation, V')
    p.add_argument('--seed', type=int, default=0)
    a = p.parse_args()

    R = R0 + a.rpot
    alpha = a.alpha or C2 / C1
    beta = a.beta or C2 * R * R / L
    unit = C2 * R                      # seconds per unit dimensionless time
    dt = 1e-6
    dtau = dt / unit
    print(f'R_total = {R:.0f} ohm  alpha = {alpha:.2f}  beta = {beta:.2f}  '
          f'm0 = {a.m0:.3f}  m1 = {a.m1:.3f}  time unit {unit * 1e6:.0f} us')

    raw, lam_dimless = simulate(alpha, beta, a.m0, a.m1, a.samples, dtau,
                                a.seed)
    lam_true = lam_dimless / unit
    rng = np.random.default_rng(a.seed)
    n = len(raw)
    t = np.arange(n) * dt
    v1 = np.round((raw[:, 0] + rng.normal(0, a.noise, n)) / a.lsb[0]) * a.lsb[0]
    v2 = np.round((raw[:, 1] + rng.normal(0, a.noise / 5, n)) / a.lsb[1]) * a.lsb[1]
    print(f'x in [{raw[:, 0].min():.2f}, {raw[:, 0].max():.2f}]  '
          f'({"double" if raw[:, 0].min() < -1 < 1 < raw[:, 0].max() else "single"} scroll)')

    tmp = tempfile.mkdtemp()
    path = os.path.join(tmp, 'sim.csv')
    np.savetxt(path, np.column_stack([t, v1, v2, 0.5 * (v1 + v2)]),
               delimiter=',', header='Time(s),CH1(V),CH2(V),CH3(V)',
               comments='', fmt='%.6e')
    M, info = maxima(path)
    per = info['period_samples']
    T = np.diff(info['t_peaks'])
    res = lyapunov(M, info['t_peaks'], quantum=info['quantum'])
    print(f'period {per} samples, {len(M)} maxima, <T> = {T.mean() * 1e6:.0f} us')
    print(f'\n  true lambda_1            {lam_true:8.0f} /s   '
          f'(lambda <T> = {lam_true * T.mean():.3f} per return)')
    if res['status'] == 'ok':
        print(f'  return map, as measured  {res["lam"]:8.0f} +/- {res["lam_err"]:.0f} /s   '
              f'({res["lam"] / lam_true:.2f} x true;  R2 = {res["r2"]:.3f}, '
              f'{res["unresolved_pct"]:.0f} % below noise)')
    else:
        print(f'  return map, as measured  no lambda: {res["status"]}')
    if per >= 20:
        ros = rosenstein(v1, v2, dt, per)
        print(f'  Rosenstein, as measured  {ros["lam"]:8.0f} +/- {ros["lam_err"]:.0f} /s   '
              f'({ros["lam"] / lam_true:.2f} x true;  {ros["n_pairs"]} pairs)')
    os.remove(path)
    os.rmdir(tmp)


if __name__ == '__main__':
    main()
