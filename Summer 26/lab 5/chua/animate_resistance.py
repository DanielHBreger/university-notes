"""Animate complete C1–C2 portraits as Rpot decreases.

Run from the repository root: python chua/animate_resistance.py
Frames use measured resistances; no interpolation between measured portraits.
"""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Rectangle
import numpy as np

from scope_data import read_scope
import simulate

HERE = Path(__file__).resolve().parent
OUT = HERE / 'animations'
WINDOW = .02
SETTLE = .1
DT = .5e-6
FPS = 8


def records():
    with (HERE / 'forward_rpot.csv').open() as source:
        clean = [(row['filename'], float(row['rpot_ohm']))
                 for row in csv.DictReader(source)
                 if row['status'] == 'ok' and float(row['residual_pct']) <= 5
                 and 0 <= float(row['rpot_ohm']) <= 1000]
    targets = np.r_[np.linspace(min(r for _, r in clean), max(r for _, r in clean), 150),
                    np.arange(700, 811, 3)]
    chosen = {min(clean, key=lambda item: abs(item[1]-target)) for target in targets}
    return sorted(chosen, key=lambda item: item[1], reverse=True)


def render(portraits, rows, kind, color, full_limits):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5.6), dpi=100)
    fig.subplots_adjust(left=.08, right=.97, bottom=.25, top=.76, wspace=.27)
    fig.suptitle(f'{kind} · phase portraits as resistance decreases', fontsize=17, y=.97)
    label = fig.text(.5, .86, '', ha='center', fontsize=14)
    source = fig.text(.5, .80, '', ha='center', fontsize=10, color='#555555')
    lines = []
    for ax, limit, title in zip(axes, [full_limits, ((-4, 4), (-1, 1))],
                               ['Full view', 'Zoom · small attractors']):
        ax.set(xlim=limit[0], ylim=limit[1], xlabel='C1 voltage, V1 (V)',
               ylabel='C2 voltage, V2 (V)', title=title)
        ax.grid(alpha=.2)
        line, = ax.plot([], [], color=color, lw=.45, alpha=.8)
        lines.append(line)
    axes[0].add_patch(Rectangle((-4, -1), 8, 2, fill=False, edgecolor='#888888',
                                linewidth=.8, linestyle='--'))
    slider = fig.add_axes([.15, .11, .70, .025])
    slider.set(xlim=(rows[0][1], rows[-1][1]), ylim=(-1, 1))
    slider.set_yticks([])
    slider.set_xlabel('Rpot (Ω) · downward sweep →', fontsize=10)
    slider.axhline(0, color='#aaaaaa', lw=2)
    marker, = slider.plot([], [], 'o', color=color, ms=7)
    fig.text(.5, .015, 'Each frame: one 20 ms phase portrait at fixed R. Playback advances R, not trajectory time.',
             ha='center', fontsize=9, color='#555555')

    def update(index):
        name, r = rows[index]
        xy = portraits[index]
        for line in lines:
            line.set_data(xy[:, 0], xy[:, 1])
        label.set_text(f'Rpot = {r:.2f} Ω     |     Rtotal = {r+simulate.R0:.2f} Ω')
        description = 'Identified bench model' if kind == 'Bench simulation' else 'Ideal model'
        source.set_text(f'forward/{name}' if kind == 'Measurements' else
                        f'{description} · continuation sweep · 100 ms settling per resistance')
        marker.set_data([r], [0])
        return *lines, label, source, marker

    anim = FuncAnimation(fig, update, frames=len(rows), blit=False)
    path = OUT / f'{kind.lower()}_resistance_sweep.gif'
    anim.save(path, writer=PillowWriter(fps=FPS))
    update(min(range(len(rows)), key=lambda i: abs(rows[i][1]-550)))
    fig.savefig(path.with_suffix('.png'))
    plt.close(fig)
    print(f'Wrote {path}', flush=True)


def render_comparison(measured, modeled, rows, full_limits, model='ideal'):
    """Synchronized columns, identical axis limits in each corresponding row."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 8.3), dpi=100)
    fig.subplots_adjust(left=.09, right=.97, bottom=.19, top=.82,
                        hspace=.40, wspace=.27)
    fig.suptitle('Measurements and simulation at equal resistance', fontsize=17, y=.97)
    label = fig.text(.5, .91, '', ha='center', fontsize=14)
    source = fig.text(.5, .865, '', ha='center', fontsize=10, color='#555555')
    lines = []
    model_title = 'Identified bench simulation' if model == 'bench' else 'Ideal simulation'
    for col, (title, color) in enumerate([('Measurements', '#167aab'),
                                         (model_title, '#cc6530')]):
        column = []
        for row, (limits, view) in enumerate([(full_limits, 'full view'),
                                               (((-4, 4), (-1, 1)), 'zoom')]):
            ax = axes[row, col]
            ax.set(xlim=limits[0], ylim=limits[1], xlabel='C1 voltage, V1 (V)',
                   ylabel='C2 voltage, V2 (V)', title=f'{title} · {view}')
            ax.grid(alpha=.2)
            line, = ax.plot([], [], color=color, lw=.45, alpha=.8)
            column.append(line)
        axes[0, col].add_patch(Rectangle((-4, -1), 8, 2, fill=False,
                                          edgecolor='#888888', lw=.8, linestyle='--'))
        lines.append(column)
    slider = fig.add_axes([.18, .09, .64, .02])
    slider.set(xlim=(rows[0][1], rows[-1][1]), ylim=(-1, 1))
    slider.set_yticks([])
    slider.set_xlabel('Rpot (Ω) · downward sweep →', fontsize=10)
    slider.axhline(0, color='#aaaaaa', lw=2)
    marker, = slider.plot([], [], 'o', color='#444444', ms=7)
    fig.text(.5, .022,
             '20 ms portraits · matched voltage scales · simulation carries state and settles 100 ms at each R',
             ha='center', fontsize=9, color='#555555')

    def update(index):
        name, r = rows[index]
        for column, portraits in zip(lines, (measured, modeled)):
            xy = portraits[index]
            for line in column:
                line.set_data(xy[:, 0], xy[:, 1])
        label.set_text(f'Both panels: Rpot = {r:.2f} Ω  |  Rtotal = {r+simulate.R0:.2f} Ω')
        description = ('Identified circuit · Rayleigh inductor · diode lag and slew'
                       if model == 'bench' else 'Ideal model: C1 = 11.5 nF, C2 = 100 nF, L = 18 mH')
        source.set_text(f'forward/{name}  |  {description}')
        marker.set_data([r], [0])

    path = OUT / ('comparison_equal_r_bench.gif' if model == 'bench' else 'comparison_equal_r.gif')
    anim = FuncAnimation(fig, update, frames=len(rows), blit=False)
    anim.save(path, writer=PillowWriter(fps=FPS))
    update(min(range(len(rows)), key=lambda i: abs(rows[i][1]-550)))
    fig.savefig(path.with_suffix('.png'))
    plt.close(fig)
    print(f'Wrote {path}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', choices=['ideal', 'bench'], default='ideal')
    parser.add_argument('--comparison-only', action='store_true',
                        help='render only the synchronized measurement/simulation comparison')
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    rows = records()
    measured, modeled = [], []
    if args.model == 'bench':
        parameters = simulate.bench_params()
        g = simulate.bench_static_g(parameters)
        start, _ = simulate.sweep_starts(g, [rows[0][1]], parameters[4])
        y = simulate.kern.bench_state(*start[:, 0], parameters)
        dt = simulate.BENCH_DT
    else:
        g = simulate.make_g()
        y, _ = simulate.sweep_starts(g, [rows[0][1]], 0.)
        dt = DT
    for index, (name, r) in enumerate(rows):
        t, xy = read_scope(HERE / 'forward' / name)
        xy = xy[t-t[0] <= WINDOW]
        measured.append(xy[::max(1, len(xy)//12000)])
        if args.model == 'bench':
            states = simulate.kern.bench_trajectory(y, simulate.R0+r, dt,
                         round((SETTLE+WINDOW)/dt), round(SETTLE/dt), parameters)
            y = states[-1].copy()  # Carry all five states, including both op-amp outputs.
            portrait = states[::20, :2]
        else:
            st, states = simulate.integrate([r], None, SETTLE+WINDOW, dt, 0., g,
                                            t_skip=SETTLE, y0=y)
            y = states[-1].copy()
            portrait = states[::4, :2, 0]
        assert np.isfinite(states).all()
        modeled.append(portrait.copy())
        if index % 20 == 0:
            print(f'Prepared {index+1}/{len(rows)} portraits; Rpot={r:.2f} Ω', flush=True)
    bounds = [np.ceil(max(np.max(np.abs(xy[:, axis])) for xy in measured+modeled)*1.08)
              for axis in range(2)]
    full_limits = [(-bound, bound) for bound in bounds]
    print(f'Shared full-view limits: {full_limits}', flush=True)
    if not args.comparison_only:
        render(measured, rows, 'Measurements', '#167aab', full_limits)
        render(modeled, rows, 'Simulation' if args.model == 'ideal' else 'Bench simulation', '#cc6530', full_limits)
    render_comparison(measured, modeled, rows, full_limits, args.model)
    metadata = dict(direction='decreasing Rpot', frames=len(rows), fps=FPS,
                    portrait_window_s=WINDOW, model=args.model, settling_s=SETTLE,
                    dt_s=dt, continuation=True, r0_ohm=simulate.R0,
                    records=[dict(filename=n, rpot_ohm=r) for n, r in rows])
    if args.model == 'bench':
        metadata.update(parameter_source='chua/identified.json',
                        parameters_SI=dict(zip(simulate.kern.BENCH_FIELDS, parameters.tolist())))
    else:
        metadata.update(c1_nF=11.5, c2_nF=100, l_mH=18, rL_ohm=0)
    metadata_name = 'resistance_sweep_bench_metadata.json' if args.model == 'bench' else 'resistance_sweep_metadata.json'
    (OUT / metadata_name).write_text(json.dumps(metadata, indent=2)+'\n')


if __name__ == '__main__':
    main()
