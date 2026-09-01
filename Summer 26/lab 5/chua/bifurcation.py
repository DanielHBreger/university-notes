"""
Bifurcation diagram: the maxima of V1 against the potentiometer resistance.

Each record contributes one vertical slice - all the local maxima M_n of the
chosen channel, plotted at that record's Rpot. A period-1 limit cycle gives a
single point, period-2 a pair, chaos a filled band, so the period-doubling
cascade and the onset of the double scroll can be read straight off the plot.

The maxima come from lorenz_map.maxima (Savitzky-Golay at period/20, then a
2 % prominence threshold) and the records are chosen by lorenz_map.collect,
so this diagram and the Lorenz maps always agree by construction.

Rpot comes from the "<folder>_rpot.csv" written by batch_rpot.py, which must
exist - it is the x axis. Records whose divider fit is poor are dropped, since
a bad Rpot puts a good slice at the wrong place on the axis; the cut is on the
fit residual (--max-residual, default 10 %). rpot.py's own "not a clean
divider" warning fires at 5 %, so the default admits a borderline band; the
run prints exactly which records that is, and --max-residual 5 excludes them.

Records that clip are dropped too (--max-clip). At the low-resistance end of
these sweeps the waveform runs past the scope's input range and sits pinned at
the rail for a third of every record; its "maxima" are all exactly the clip
level, which would draw a flat line across the diagram that looks like a
period-1 branch but is only the edge of the screen.

Give any number of folders to overlay sweeps, each drawn in its own colour and
named in the legend - a forward and a backward sweep on the same axes show
hysteresis, where the transitions sit at different resistances depending on
which way the pot was turned. With no folders named, the picker opens beside
the script and keeps asking for another until you cancel.

Usage:
    python bifurcation.py                                   # pick, repeatedly
    python bifurcation.py "sweep forward"
    python bifurcation.py "sweep forward" "sweep back"      # hysteresis
    python bifurcation.py "sweep forward" "sweep back" comparisons
    python bifurcation.py "sweep forward" --max-residual 5
    python bifurcation.py "sweep forward" --descending      # R falls to the right
"""
import argparse
import csv
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from lorenz_map import collect
from sweeplib import (CLEAN_DIVIDER_MAX_PCT, resolve_folders, sibling,
                      sweep_colors)


def add_record_args(p):
    """The options that decide which records count, shared with lyapunov.py."""
    p.add_argument('--ch', default='CH1', help='channel to take maxima of')
    p.add_argument('--prominence', type=float, default=0.02,
                   help='peak prominence, fraction of the range')
    p.add_argument('--period-samples', type=int,
                   help='fixed period in samples (default: per-record autocorr)')
    p.add_argument('--max-residual', type=float, default=10.0,
                   help='drop records whose Rpot fit residual exceeds this %%')
    p.add_argument('--max-clip', type=float, default=2.0,
                   help='drop records with more than this %% of samples at the rail')
    p.add_argument('--max-files', type=int, help='stop after this many records')
    return p


def report(label, records, dropped):
    """What was admitted, what was borderline, and what was thrown away."""
    if not records:
        print(f'{label}: nothing usable')
        return
    rs = [rec.rpot for rec in records]
    n = sum(len(rec.M) for rec in records)
    print(f'{label}: {len(records)} records, {n} maxima, '
          f'R = {min(rs):.1f} .. {max(rs):.1f} ohm')
    borderline = [rec for rec in records if rec.residual > CLEAN_DIVIDER_MAX_PCT]
    if borderline:
        print(f'  admitted with residual over {CLEAN_DIVIDER_MAX_PCT:.0f} % '
              f'(--max-residual {CLEAN_DIVIDER_MAX_PCT:.0f} drops these):')
        for rec in borderline:
            print(f'    {rec.name:14s} R={rec.rpot:7.1f} ohm  '
                  f'residual {rec.residual:.1f} %')
    if dropped:
        print(f'  dropped {len(dropped)}:')
        for name, why in dropped:
            print(f'    {name:14s} {why}')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('folders', nargs='*',
                   help='one or more sweep folders (omit to pick them)')
    p.add_argument('-o', '--out',
                   help='output PNG (default: <first folder>_bifurcation.png)')
    add_record_args(p)
    p.add_argument('--descending', action='store_true',
                   help='plot R decreasing to the right')
    p.add_argument('--xlim', nargs=2, type=float, metavar=('LO', 'HI'),
                   help='resistance range to show, ohm')
    p.add_argument('--ylim', nargs=2, type=float, metavar=('LO', 'HI'),
                   help='voltage range to show, V')
    p.add_argument('--size', type=float, default=0.6, help='marker size')
    p.add_argument('--alpha', type=float, default=0.18, help='marker alpha')
    a = p.parse_args()

    folders = resolve_folders(a.folders)
    out = a.out or sibling(folders[0], '_bifurcation.png')
    csv_out = os.path.splitext(out)[0] + '_points.csv'
    colors = sweep_colors(len(folders))

    fig, ax = plt.subplots(figsize=(10, 5.5))
    rows = []

    for k, folder in enumerate(folders):
        label = os.path.basename(folder)
        print(f'\n--- {label} ---')
        records, dropped = collect(folder, a.ch, a.prominence, a.period_samples,
                                   a.max_residual, a.max_clip, a.max_files)
        report(label, records, dropped)
        if not records:
            continue
        # One scatter for the whole sweep; per-record calls would be far slower.
        xs = np.concatenate([np.full(len(r.M), r.rpot) for r in records])
        ys = np.concatenate([r.M for r in records])
        ax.scatter(xs, ys, s=a.size, alpha=a.alpha, lw=0, color=colors[k],
                   label=f'{label}  ({len(records)} records)', rasterized=True)
        for rec in records:
            rows.extend((label, rec.name, round(rec.rpot, 2),
                         round(rec.residual, 3), round(float(m), 6))
                        for m in rec.M)

    if not rows:
        raise SystemExit('no usable records in any folder')

    ax.set_xlabel('$R_{pot}$  ($\\Omega$)')
    ax.set_ylabel(f'maxima of $V_{{C1}}$ ({a.ch})  (V)')
    ax.set_title('Bifurcation diagram - Chua circuit')
    if a.xlim:
        ax.set_xlim(*a.xlim)
    if a.ylim:
        ax.set_ylim(*a.ylim)
    if a.descending:
        ax.invert_xaxis()
    leg = ax.legend(loc='best', fontsize=8, framealpha=0.9)
    for h in leg.legend_handles:
        h.set_alpha(1.0)
        h.set_sizes([12])
    ax.grid(alpha=0.15, lw=0.5)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)

    with open(csv_out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['sweep', 'filename', 'rpot_ohm', 'rpot_residual_pct', 'max_v'])
        w.writerows(rows)

    print(f'\n{len(rows)} maxima plotted\nplot   -> {out}\npoints -> {csv_out}')


if __name__ == '__main__':
    main()
