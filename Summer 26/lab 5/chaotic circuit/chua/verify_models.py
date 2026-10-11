"""
Checks behind the simulation section (RESULTS.md, "Simulation: verification of 27 September").

1. Both compiled RK4 kernels against scipy's DOP853 (rtol 1e-10) over 1 ms, on a
   periodic, a chaotic and a large-cycle state, and against themselves at a
   quarter of the step.
2. The linear-stability code against the plan's Hopf formula (nominal parts,
   Table 1's left shoulder).
3. The bench model's Hopf point from a numerical Jacobian of the full kernel
   equations at its equilibrium, against hopf_point on the linearised circuit
   (constant L0, r0, stage A as a capacitance).
4. Node 1's capacitance (identify.py's loop-integral estimators) on records of
   every regime, against the bench model with its node-1 law at the same
   resistance.

Output: verify_models.txt beside this script (about two minutes).

Usage:
    python verify_models.py
"""
import os
import sys

import numpy as np
from scipy.integrate import DOP853, solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import simulate as S
import integration as K
import identify as I

HERE = os.path.dirname(os.path.abspath(__file__))
R0 = S.R0


def ideal_rhs(g, Rt, c1, c2, l):
    """Right-hand side of the ideal model for solve_ivp, independent of the compiled kernel."""
    vk, ik = g.knots

    def f(t, y):
        v1, v2, i = y
        ir = (v2 - v1) / Rt
        return [(ir - np.interp(v1, vk, ik)) / c1, (-ir - i) / c2, v2 / l]
    return f


def bench_rhs(P, Rt, tracker):
    """
    Right-hand side of the bench model's six continuous states for scipy, written
    out again from its equations; tracker[0] is v1's last turning point, which
    the caller updates between steps as the kernel does.
    """
    c1, c2, l0, nu, r0, rho, RA, AA, RB, AB, vpA, vnA, vpB, vnB, sr, tauA, tauB, tauM, a0, a1, xk = P

    def rate(target, vo, tau, vn, vp):
        r = max(-sr, min(sr, (target - vo) / tau))
        return 0.0 if (vo >= vp and r > 0) or (vo <= vn and r < 0) else r

    def f(t, y):
        v1, v2, i, va, vb, im = y
        va, vb = min(max(va, vnA), vpA), min(max(vb, vnB), vpB)
        ir = (v2 - v1) / Rt
        d = i - im
        x = abs(v1 - tracker[0])
        cap = c1 + a0 * min(x, xk) + a1 * max(x - xk, 0.0)
        return [(ir - (v1 - va) / RA - (v1 - vb) / RB) / cap, (-ir - i) / c2,
                (v2 - r0 * i - rho * abs(d) * d) / (l0 + nu * abs(d)),
                rate(AA * v1, va, tauA, vnA, vpA), rate(AB * v1, vb, tauB, vnB, vpB), d / tauM]
    return f


def bench_dop853(P, Rt, y, n, dt):
    """v1 at n + 1 points dt apart from DOP853, the turning-point tracker updated after every accepted step."""
    tracker = [y[6], y[7], y[8]]
    solver = DOP853(bench_rhs(P, Rt, tracker), 0.0, y[:6], n * dt, rtol=1e-10, atol=1e-12, max_step=dt / 2)
    ts, v1 = [0.0], [y[0]]
    while solver.status == 'running':
        t_old = solver.t
        solver.step()
        tracker[:] = K._turn(solver.y[0], *tracker, solver.t - t_old, K.TURN_V, K.TAU_R)
        ts.append(solver.t)
        v1.append(solver.y[0])
    return np.interp(np.arange(n + 1) * dt, ts, v1)


def negative_equilibrium(g, Rt, r0):
    """v1 of the circuit's negative outer equilibrium."""
    return [e for e in S.equilibria(g, Rt, r0) if e < 0][0]


def check_kernels(lines):
    """1. Both kernels against DOP853 and against themselves at a quarter of the step."""
    lines.append('1. kernels against DOP853 (rtol 1e-10) over 1 ms; max |dV1|')
    c1, c2, l = S.NOMINAL.c1, S.NOMINAL.c2, S.NOMINAL.l
    g = S.table1_g()
    vk, ik = g.knots
    for rpot, what in [(950.0, 'period 2'), (700.0, 'double scroll'), (300.0, 'large cycle')]:
        Rt = R0 + rpot
        y0 = np.array([7.0, 6.0, 0.0]) if rpot < 400 else np.array([-2.0, 0.1, -0.5e-3])
        y = K.trajectory(y0, Rt, 0.0, c1, c2, l, vk, ik, 0.5e-6, 60000, 60000)[-1]
        n = 2000
        ker = K.trajectory(y.copy(), Rt, 0.0, c1, c2, l, vk, ik, 0.5e-6, n, 0)
        fine = K.trajectory(y.copy(), Rt, 0.0, c1, c2, l, vk, ik, 0.125e-6, 4 * n, 0)[::4]
        sol = solve_ivp(ideal_rhs(g, Rt, c1, c2, l), (0, n * 0.5e-6), y, method='DOP853', rtol=1e-10,
                        atol=1e-12, t_eval=np.arange(n + 1) * 0.5e-6)
        lines.append(f'   plan\'s model, {what:13s} ({rpot:.0f} ohm): kernel dt 0.5 us vs DOP853 '
                     f'{np.abs(ker[:, 0] - sol.y[0]).max():.1e} V, vs dt 0.125 us {np.abs(ker[:, 0] - fine[:, 0]).max():.1e} V')
    P = S.bench_params()
    for rpot, what in [(806.0, 'period 1'), (550.0, 'double scroll'), (320.0, 'large cycle')]:
        Rt = R0 + rpot
        y0 = K.bench_state(7.0, 6.0, 0.0, P) if rpot < 400 else K.bench_state(-2.0, 0.1, -0.5e-3, P)
        y = K.bench_trajectory(y0, Rt, 0.1e-6, 300000, 300000, P)[-1]
        n = 10000
        ker = K.bench_trajectory(y.copy(), Rt, 0.1e-6, n, 0, P)
        fine = K.bench_trajectory(y.copy(), Rt, 0.025e-6, 4 * n, 0, P)[::4]
        ref = bench_dop853(P, Rt, y, n, 0.1e-6)
        lines.append(f'   bench model,  {what:13s} ({rpot:.0f} ohm): kernel dt 0.1 us vs DOP853 '
                     f'{np.abs(ker[:, 0] - ref).max():.1e} V, vs dt 0.025 us {np.abs(ker[:, 0] - fine[:, 0]).max():.1e} V')


def check_hopf_formula(lines):
    """2. hopf_point against the plan's closed-form Hopf point."""
    g = S.table1_g()
    Gb, C1, C2, L = g.G[1], S.NOMINAL.c1, S.NOMINAL.c2, S.NOMINAL.l
    Rt = -Gb * (1 + C1 / C2) / (Gb ** 2 + C1 ** 2 / (C2 * L))
    f = np.sqrt(Gb / (Rt * C1 * C2) + 1 / (L * C2)) / (2 * np.pi)
    h = S.hopf_point(g, S.NOMINAL)
    lines.append(f'2. Hopf point, nominal parts, G_b = {Gb * 1e3:.5f} mS: plan formula {Rt - R0:.2f} ohm, {f:.1f} Hz; '
                 f'hopf_point {h[0] - R0:.2f} ohm, {1 / h[2]:.1f} Hz')


def check_bench_linearisation(lines):
    """3. The bench model's Hopf point from a numerical Jacobian of the full equations."""
    P = S.bench_params()
    g = S.bench_static_g(P)
    r0 = P[4]

    def pair(rpot):
        Rt = R0 + rpot
        v = negative_equilibrium(g, Rt, r0)
        i = v / (Rt + r0)
        y = K.bench_state(v, r0 * i, i, P)[:6]
        f = bench_rhs(P, Rt, [v])                  # at rest: no swing
        f0 = np.array(f(0, y))
        idx = [0, 1, 2, 3, 5]                      # stage B sits at its rail
        J = np.zeros((5, 5))
        for a, k in enumerate(idx):
            h = 1e-7 * (1e-3 if k in (2, 5) else 1.0)
            yp = y.copy()
            yp[k] += h
            J[:, a] = (np.array(f(0, yp))[idx] - f0[idx]) / h
        ev = np.linalg.eigvals(J)
        ev = ev[np.abs(ev.imag) > 1]
        return ev.real.max(), np.abs(ev.imag).max()
    Rh = brentq(lambda r: pair(r)[0], 700, 1000)
    bh = S.hopf_point(g, S.bench_small_signal(P))
    lines.append(f'3. bench model Hopf point: full Jacobian {Rh:.2f} ohm, {pair(Rh)[1] / 2 / np.pi:.1f} Hz; '
                 f'hopf_point {bh[0] - R0:.2f} ohm, {1 / bh[2]:.1f} Hz')


def model_node1(P, g, rpot, large):
    """
    identify.py's node-1 capacitance of the bench model at rpot: 100 ms settled, 100 ms
    sampled every 400 ns like the scope, cycle estimator if the run is periodic,
    window estimator otherwise. Starts where the sweeps start: on the large cycle, or
    at the negative equilibrium (the origin below the three-equilibria range).
    """
    fwd, back = S.sweep_starts(g, S.bench_small_signal(P), [rpot])
    y0 = K.bench_state(*(back if large else fwd)[:, 0], P)
    tr = K.bench_trajectory(y0, R0 + rpot, 0.1e-6, 2000000, 1000000, P)[::4]
    if np.ptp(tr[:, 0]) < 0.05:
        return np.nan, 'at rest', tr[:, 0].min(), tr[:, 0].max()
    t = np.arange(len(tr)) * 0.4e-6
    M, _ = I.maxima_from_arrays(t, tr[:, 0])
    est = I.cycle_node1(t, tr[:, 0], tr[:, 1], R0 + rpot) if len(M) and np.ptp(M) < 0.12 else None
    est = est or I.window_node1(t, tr[:, 0], tr[:, 1], R0 + rpot)
    return est.c, est.method, tr[:, 0].min(), tr[:, 0].max()


def check_loop_c1(lines):
    """4. Node 1's loop-integral capacitance on records of every regime and on the model at the same resistance."""
    lines.append('4. node 1 capacitance (loop integral, nF; cycle or window estimator): records against the bench model')
    picks = {'forward': (14, 25, 44, 60, 92, 200, 270, 342, 404, 411), 'back': (25, 15)}
    meas = {}
    for sweep, numbers in picks.items():
        for path, rp, lab in I.sweep_records(sweep):
            if int(os.path.basename(path)[5:-4]) in numbers:
                r = I.analyse_record(path, rp)
                if r:
                    d = r['node1']
                    meas[lab] = (rp, d.c * 1e9, d.method, d.v1.min(), d.v1.max())
    P = S.bench_params()
    g = S.bench_static_g(P)
    for lab, (rp, c, how, lo, hi) in sorted(meas.items(), key=lambda x: -x[1][0]):
        cm, how_m, mlo, mhi = model_node1(P, g, rp, hi > 5)
        lines.append(f'   {lab:22s} {rp:6.1f} ohm, V1 {lo:+.2f}..{hi:+.2f} V: measured {c:6.3f} ({how}), '
                     f'model {cm * 1e9:6.3f} ({how_m}, model V1 {mlo:+.2f}..{mhi:+.2f})')


def main():
    lines = ['verify_models.py: checks behind the simulation section']
    for check in (check_kernels, check_hopf_formula, check_bench_linearisation, check_loop_c1):
        check(lines)
        print('\n'.join(lines[-4:]), flush=True)
    txt = '\n'.join(lines)
    with open(os.path.join(HERE, 'verify_models.txt'), 'w') as fh:
        fh.write(txt + '\n')
    print('\n' + txt)


if __name__ == '__main__':
    main()
