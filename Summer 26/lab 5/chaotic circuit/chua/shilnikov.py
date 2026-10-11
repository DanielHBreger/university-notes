"""
N6. Eigenvalues at the three equilibria and the Shilnikov ratio.

At each fixed point the Jacobian of

    C1 v1' = (v2 - v1)/Rt - g(v1),  C2 v2' = (v1 - v2)/Rt - iL,  L iL' = v2 - rL iL

has one real eigenvalue gamma and a complex pair sigma +- j omega. The double
scroll needs the origin to be a saddle-focus with a one-dimensional unstable
manifold (gamma > 0, sigma < 0) and the outer equilibria the opposite kind
(gamma < 0, sigma > 0). Shilnikov's theorem then gives horseshoes near a
homoclinic orbit when

    |sigma| < |gamma|        (ratio |sigma| / |gamma| below 1)

at the equilibrium the orbit returns to. Both ratios are printed for every
resistance, for the plan's model (Table 1's element, nominal C1, C2
and L) and for the identified circuit (identified.json, small-signal values).

Output: shilnikov.txt beside this script.

Usage:
    python shilnikov.py
    python shilnikov.py --r 806 745 668 500 328
"""
import argparse
import os

import numpy as np

import simulate as sim

HERE = os.path.dirname(os.path.abspath(__file__))


def eigen_rows(g, Rt, circuit):
    """One dict per equilibrium: v1, local slope G, gamma, sigma, omega and the ratio |sigma|/|gamma|."""
    rows = []
    for v in sim.equilibria(g, Rt, circuit.rl):
        ev, G = sim.stability(g, v, Rt, circuit)
        real = ev[np.abs(ev.imag) < 1.0].real
        cplx = ev[np.abs(ev.imag) >= 1.0]
        gamma = float(real[np.argmax(np.abs(real))]) if len(real) else np.nan
        sigma = float(cplx.real.max()) if len(cplx) else np.nan
        omega = float(np.abs(cplx.imag).max()) if len(cplx) else np.nan
        rows.append(dict(v1=float(v), G=float(G), gamma=gamma, sigma=sigma, omega=omega,
                         ratio=abs(sigma) / abs(gamma) if np.isfinite(sigma) and gamma else np.nan))
    return rows


def kind(row):
    """Saddle-focus type from the signs, or 'stable' / 'other'."""
    if not np.isfinite(row['sigma']):
        return 'real eigenvalues only'
    if row['gamma'] > 0 > row['sigma']:
        return 'saddle-focus, 1-D unstable'
    if row['gamma'] < 0 < row['sigma']:
        return 'saddle-focus, 2-D unstable'
    if row['gamma'] < 0 and row['sigma'] < 0:
        return 'stable focus'
    return 'unstable in every direction'


def models(identified):
    """{name: (g, circuit)} of the plan's model and, when identified.json exists, the identified circuit."""
    plan = sim.NOMINAL
    out = {f'plan model (Table 1 element, C1 {plan.c1 * 1e9:g} nF, C2 {plan.c2 * 1e9:g} nF, L {plan.l * 1e3:g} mH)':
           (sim.table1_g(), plan)}
    if os.path.exists(identified):
        P = sim.bench_params(identified)
        c = sim.bench_small_signal(P)
        name = f'identified circuit (C1 {c.c1 * 1e9:.1f} nF, C2 {c.c2 * 1e9:.1f} nF, L0 {c.l * 1e3:.1f} mH, r0 {c.rl:.1f} ohm)'
        out[name] = (sim.bench_static_g(P), c)
    return out


def verdict(rows):
    """The double-scroll geometry and the Shilnikov condition at the three equilibria, as one line."""
    outer_left, origin, outer_right = rows
    geometry = origin['gamma'] > 0 > origin['sigma'] and outer_left['gamma'] < 0 < outer_left['sigma']
    at_origin = 'holds' if origin['ratio'] < 1 else 'fails'
    at_outer = 'holds' if max(outer_left['ratio'], outer_right['ratio']) < 1 else 'fails'
    return (f'{"":>6} {"":>8} double-scroll geometry {"yes" if geometry else "no"}; '
            f'Shilnikov at the origin {at_origin} ({origin["ratio"]:.2f}), at the outer equilibria '
            f'{at_outer} ({outer_left["ratio"]:.2f}, {outer_right["ratio"]:.2f})')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--r', type=float, nargs='+', default=[806, 770, 745, 700, 668, 600, 500, 400, 328],
                   help='potentiometer values (ohm)')
    p.add_argument('--identified', default=os.path.join(HERE, 'identified.json'))
    a = p.parse_args()

    lines = ['shilnikov.py: eigenvalues at the three equilibria, gamma real, sigma +- j omega the pair; '
             'ratio = |sigma| / |gamma| (Shilnikov needs < 1)']
    for name, (g, circuit) in models(a.identified).items():
        lines += ['', name,
                  f'{"Rpot":>6} {"v1* (V)":>8} {"G (mS)":>7} {"gamma (1/s)":>12} {"sigma (1/s)":>12} '
                  f'{"f (kHz)":>8} {"ratio":>6}  type']
        for r in a.r:
            rows = eigen_rows(g, circuit.rt(r), circuit)
            for k, row in enumerate(rows):
                lines.append(f'{r if k == 0 else "":>6} {row["v1"]:8.3f} {row["G"] * 1e3:7.3f} {row["gamma"]:12.0f} '
                             f'{row["sigma"]:12.0f} {row["omega"] / 2 / np.pi * 1e-3:8.2f} {row["ratio"]:6.3f}  {kind(row)}')
            if len(rows) == 3:
                lines.append(verdict(rows))
            elif len(rows) == 1:
                lines.append(f'{"":>6} {"":>8} one equilibrium: no double scroll possible')
    txt = '\n'.join(lines)
    print(txt)
    with open(os.path.join(HERE, 'shilnikov.txt'), 'w') as fh:
        fh.write(txt + '\n')


if __name__ == '__main__':
    main()
