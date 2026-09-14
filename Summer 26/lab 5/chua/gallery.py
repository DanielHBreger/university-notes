"""
N2. Phase-portrait gallery of measured records: limit cycle, period-2,
single scroll, double scroll (and the large outer cycle), each as V2 against
V1 with a short stretch of V1(t) underneath.

The records are the clean forward-sweep records nearest the requested
resistances (from forward_rpot.csv, divider residual under 5 %). Defaults are
inside each regime as read off the bifurcation diagram and cascade_periods.py.

Output: gallery.png beside this script.

Usage:
    python gallery.py
    python gallery.py --r 806 760 700 550 190 --sweep forward
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scope_data import read_scope
from sweeplib import load_rpot

HERE = os.path.dirname(os.path.abspath(__file__))
REGIMES = [('limit cycle (period 1)', 806.0), ('period 2', 760.0), ('single scroll', 700.0),
           ('double scroll', 550.0), ('large outer cycle', 190.0)]


def nearest(folder, rpot, max_residual=5.0):
    """(name, Rpot) of the clean record in `folder` nearest `rpot`."""
    table = {n: r for n, (r, res) in load_rpot(folder, required=True).items() if res <= max_residual}
    name = min(table, key=lambda n: abs(table[n] - rpot))
    return name, table[name]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sweep', default='forward')
    p.add_argument('--r', type=float, nargs='+', default=[r for _, r in REGIMES])
    p.add_argument('--labels', nargs='+', default=None)
    p.add_argument('--xy-samples', type=int, default=150000, help='samples drawn in the XY panel')
    p.add_argument('--window', type=float, default=4e-3, help='seconds of V1(t) shown')
    p.add_argument('--out', default=os.path.join(HERE, 'gallery.png'))
    a = p.parse_args()
    labels = a.labels or ([lab for lab, _ in REGIMES] if len(a.r) == len(REGIMES) else [''] * len(a.r))
    folder = os.path.join(HERE, a.sweep)

    n = len(a.r)
    fig, axes = plt.subplots(2, n, figsize=(3.4 * n, 6.4), squeeze=False, height_ratios=[3, 1.3])
    for k, (r, lab) in enumerate(zip(a.r, labels)):
        name, rp = nearest(folder, r)
        t, d = read_scope(os.path.join(folder, name))
        v1, v2 = d[:, 0], d[:, 1]
        ax = axes[0, k]
        ax.plot(v1[:a.xy_samples], v2[:a.xy_samples], lw=0.25, color='C0')
        ax.set_title(f'{lab}\n{a.sweep}/{name}, Rpot = {rp:.1f} ohm', fontsize=9)
        ax.set_xlabel('V1 (V)'); ax.set_ylabel('V2 (V)'); ax.grid(alpha=0.3)
        ax = axes[1, k]
        m = t <= t[0] + a.window
        ax.plot((t[m] - t[0]) * 1e3, v1[m], lw=0.6, color='C0')
        ax.set_xlabel('t (ms)'); ax.set_ylabel('V1 (V)'); ax.grid(alpha=0.3)
        print(f'{lab:24s} {a.sweep}/{name}  Rpot = {rp:.1f} ohm  V1 {v1.min():.2f}..{v1.max():.2f} V')
    fig.suptitle('Measured phase portraits (V2 against V1) and V1(t)')
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print('wrote', a.out)


if __name__ == '__main__':
    main()
