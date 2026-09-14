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
resistance, for the plan's model (the V-I segments, C1 = 11.5 nF, nominal C2
and L) and for the identified circuit (identified.json, small-signal values).

Output: shilnikov.txt beside this script.

Usage:
    python shilnikov.py
    python shilnikov.py --r 806 745 668 500 328
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import simulate as s

HERE = os.path.dirname(os.path.abspath(__file__))


def eigen_rows(g, Rt, rL, c1, c2, l):
    """One dict per equilibrium: v1, local slope, gamma, sigma, omega, ratio."""
    rows = []
    for v in s.equilibria(g, Rt, rL):
        ev, G = s.stability(g, v, Rt, rL, c1=c1, c2=c2, l=l)
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
    """{name: (g, rL, c1, c2, l)} of the two element / component sets."""
    out = {'plan model (V-I segments, C1 11.5 nF, C2 100 nF, L 18 mH)':
           (s.make_g(), 0.0, s.C1, s.C2, s.L)}
    if os.path.exists(identified):
        P = s.bench_params(identified)
        c1a, c2b, l0, r0 = s.bench_small_signal(P)
        out[f'identified circuit (C1 {c1a*1e9:.1f} nF, C2 {c2b*1e9:.1f} nF, L0 {l0*1e3:.1f} mH, r0 {r0:.1f} ohm)'] = \
            (s.bench_static_g(P), r0, c1a, c2b, l0)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--r', type=float, nargs='+', default=[806, 770, 745, 700, 668, 600, 500, 400, 328],
                   help='potentiometer values (ohm)')
    p.add_argument('--identified', default=os.path.join(HERE, 'identified.json'))
    a = p.parse_args()

    lines = ['shilnikov.py: eigenvalues at the three equilibria, gamma real, sigma +- j omega the pair; '
             'ratio = |sigma| / |gamma| (Shilnikov needs < 1)']
    for name, (g, rL, c1, c2, l) in models(a.identified).items():
        lines += ['', name,
                  f'{"Rpot":>6} {"v1* (V)":>8} {"G (mS)":>7} {"gamma (1/s)":>12} {"sigma (1/s)":>12} '
                  f'{"f (kHz)":>8} {"ratio":>6}  type']
        for r in a.r:
            rows = eigen_rows(g, s.R0 + r, rL, c1, c2, l)
            for k, row in enumerate(rows):
                lines.append(f'{r if k == 0 else "":>6} {row["v1"]:8.3f} {row["G"]*1e3:7.3f} {row["gamma"]:12.0f} '
                             f'{row["sigma"]:12.0f} {row["omega"]/2/np.pi*1e-3:8.2f} {row["ratio"]:6.3f}  {kind(row)}')
            if len(rows) == 3:
                ok = (rows[1]['gamma'] > 0 > rows[1]['sigma'] and rows[0]['gamma'] < 0 < rows[0]['sigma'])
                lines.append(f'{"":>6} {"":>8} double-scroll geometry {"yes" if ok else "no"}; '
                             f'Shilnikov at the origin {"holds" if rows[1]["ratio"] < 1 else "fails"} '
                             f'({rows[1]["ratio"]:.2f}), at the outer equilibria '
                             f'{"holds" if max(rows[0]["ratio"], rows[2]["ratio"]) < 1 else "fails"} '
                             f'({rows[0]["ratio"]:.2f}, {rows[2]["ratio"]:.2f})')
            elif len(rows) == 1:
                lines.append(f'{"":>6} {"":>8} one equilibrium: no double scroll possible')
    txt = '\n'.join(lines)
    print(txt)
    with open(os.path.join(HERE, 'shilnikov.txt'), 'w') as fh:
        fh.write(txt + '\n')


if __name__ == '__main__':
    main()
