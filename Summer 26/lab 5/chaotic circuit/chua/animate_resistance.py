"""
Animate complete C1-C2 portraits as Rpot decreases.

Each frame is a 20 ms portrait at one measured resistance of the clean forward
sweep; nothing is interpolated between records. The simulation carries its
state from one resistance to the next and settles 100 ms at each.

Run from the repository root:
    python chua/animate_resistance.py                                   # measurements, ideal model, comparison
    python chua/animate_resistance.py --comparison-only                 # only the side-by-side comparison
    python chua/animate_resistance.py --model bench --comparison-only   # the identified circuit instead
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
WINDOW = .02        # s, one portrait
SETTLE = .1         # s, settling at each resistance
DT = .5e-6          # s, RK4 step of the ideal model
FPS = 8
ZOOM = ((-4, 4), (-1, 1))
# The ideal model these animations show: the fitted V-I lines joined at their intersections
# (simulate.make_g) with C1 = 11.5 nF, the value that put its Hopf point at the onset rpot.py
# reported. It is the presentation's model, not the report's (simulate.py --model ideal).
ANIMATION_CIRCUIT = simulate.Circuit(c1=11.5e-9)


def records():
    """(filename, Rpot) of the clean forward records nearest an even Rpot grid, densest at 700-810 ohm, by falling Rpot."""
    with (HERE / 'forward_rpot.csv').open() as source:
        clean = [(row['filename'], float(row['rpot_ohm']))
                 for row in csv.DictReader(source)
                 if row['status'] == 'ok' and float(row['residual_pct']) <= 5
                 and 0 <= float(row['rpot_ohm']) <= 1000]
    targets = np.r_[np.linspace(min(r for _, r in clean), max(r for _, r in clean), 150),
                    np.arange(700, 811, 3)]
    chosen = {min(clean, key=lambda item: abs(item[1] - target)) for target in targets}
    return sorted(chosen, key=lambda item: item[1], reverse=True)


def nearest_index(rows, rpot):
    """Index of the row whose Rpot is nearest `rpot` (the still frame)."""
    return min(range(len(rows)), key=lambda i: abs(rows[i][1] - rpot))


def add_slider(fig, rect, rows, color):
    """The Rpot indicator under the portraits; returns its marker."""
    slider = fig.add_axes(rect)
    slider.set(xlim=(rows[0][1], rows[-1][1]), ylim=(-1, 1))
    slider.set_yticks([])
    slider.set_xlabel('Rpot (Ω) · downward sweep →', fontsize=10)
    slider.axhline(0, color='#aaaaaa', lw=2)
    marker, = slider.plot([], [], 'o', color=color, ms=7)
    return marker


def save_animation(fig, update, rows, path):
    """Write the GIF and a PNG still near 550 ohm."""
    anim = FuncAnimation(fig, update, frames=len(rows), blit=False)
    anim.save(path, writer=PillowWriter(fps=FPS))
    update(nearest_index(rows, 550))
    fig.savefig(path.with_suffix('.png'))
    plt.close(fig)
    print(f'Wrote {path}', flush=True)


def render(portraits, rows, kind, color, full_limits):
    """One source (measurements or a model) in a full view and a zoom."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 5.6), dpi=100)
    fig.subplots_adjust(left=.08, right=.97, bottom=.25, top=.76, wspace=.27)
    fig.suptitle(f'{kind} · phase portraits as resistance decreases', fontsize=17, y=.97)
    label = fig.text(.5, .86, '', ha='center', fontsize=14)
    source = fig.text(.5, .80, '', ha='center', fontsize=10, color='#555555')
    lines = []
    for ax, limit, title in zip(axes, [full_limits, ZOOM], ['Full view', 'Zoom · small attractors']):
        ax.set(xlim=limit[0], ylim=limit[1], xlabel='C1 voltage, V1 (V)',
               ylabel='C2 voltage, V2 (V)', title=title)
        ax.grid(alpha=.2)
        line, = ax.plot([], [], color=color, lw=.45, alpha=.8)
        lines.append(line)
    axes[0].add_patch(Rectangle((-4, -1), 8, 2, fill=False, edgecolor='#888888',
                                linewidth=.8, linestyle='--'))
    marker = add_slider(fig, [.15, .11, .70, .025], rows, color)
    fig.text(.5, .015, 'Each frame: one 20 ms phase portrait at fixed R. Playback advances R, not trajectory time.',
             ha='center', fontsize=9, color='#555555')
    description = 'Identified bench model' if kind == 'Bench simulation' else 'Ideal model'

    def update(index):
        name, r = rows[index]
        xy = portraits[index]
        for line in lines:
            line.set_data(xy[:, 0], xy[:, 1])
        label.set_text(f'Rpot = {r:.2f} Ω     |     Rtotal = {r + simulate.R0:.2f} Ω')
        source.set_text(f'forward/{name}' if kind == 'Measurements' else
                        f'{description} · continuation sweep · 100 ms settling per resistance')
        marker.set_data([r], [0])
        return *lines, label, source, marker

    save_animation(fig, update, rows, OUT / f'{kind.lower()}_resistance_sweep.gif')


def render_comparison(measured, modeled, rows, full_limits, model='ideal'):
    """Measurements and the model side by side, with identical axis limits in each row."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 8.3), dpi=100)
    fig.subplots_adjust(left=.09, right=.97, bottom=.19, top=.82,
                        hspace=.40, wspace=.27)
    fig.suptitle('Measurements and simulation at equal resistance', fontsize=17, y=.97)
    label = fig.text(.5, .91, '', ha='center', fontsize=14)
    source = fig.text(.5, .865, '', ha='center', fontsize=10, color='#555555')
    lines = []
    model_title = 'Identified bench simulation' if model == 'bench' else 'Ideal simulation'
    for col, (title, color) in enumerate([('Measurements', '#167aab'), (model_title, '#cc6530')]):
        column = []
        for row, (limits, view) in enumerate([(full_limits, 'full view'), (ZOOM, 'zoom')]):
            ax = axes[row, col]
            ax.set(xlim=limits[0], ylim=limits[1], xlabel='C1 voltage, V1 (V)',
                   ylabel='C2 voltage, V2 (V)', title=f'{title} · {view}')
            ax.grid(alpha=.2)
            line, = ax.plot([], [], color=color, lw=.45, alpha=.8)
            column.append(line)
        axes[0, col].add_patch(Rectangle((-4, -1), 8, 2, fill=False,
                                         edgecolor='#888888', lw=.8, linestyle='--'))
        lines.append(column)
    marker = add_slider(fig, [.18, .09, .64, .02], rows, '#444444')
    fig.text(.5, .022,
             '20 ms portraits · matched voltage scales · simulation carries state and settles 100 ms at each R',
             ha='center', fontsize=9, color='#555555')
    c = ANIMATION_CIRCUIT
    description = ('Identified circuit · Rayleigh inductor · diode lag and slew' if model == 'bench' else
                   f'Ideal model: C1 = {c.c1 * 1e9:g} nF, C2 = {c.c2 * 1e9:g} nF, L = {c.l * 1e3:g} mH')

    def update(index):
        name, r = rows[index]
        for column, portraits in zip(lines, (measured, modeled)):
            xy = portraits[index]
            for line in column:
                line.set_data(xy[:, 0], xy[:, 1])
        label.set_text(f'Both panels: Rpot = {r:.2f} Ω  |  Rtotal = {r + simulate.R0:.2f} Ω')
        source.set_text(f'forward/{name}  |  {description}')
        marker.set_data([r], [0])

    name = 'comparison_equal_r_bench.gif' if model == 'bench' else 'comparison_equal_r.gif'
    save_animation(fig, update, rows, OUT / name)


def portraits(rows, model):
    """(measured, modelled) portraits at every row's Rpot; the model carries its state along. Also returns dt."""
    if model == 'bench':
        parameters = simulate.bench_params()
        start, _ = simulate.sweep_starts(simulate.bench_static_g(parameters),
                                         simulate.bench_small_signal(parameters), [rows[0][1]])
        y = simulate.kern.bench_state(*start[:, 0], parameters)
        dt = simulate.BENCH_DT
    else:
        g = simulate.make_g()
        y, _ = simulate.sweep_starts(g, ANIMATION_CIRCUIT, [rows[0][1]])
        dt = DT
    measured, modeled = [], []
    for index, (name, r) in enumerate(rows):
        t, xy = read_scope(HERE / 'forward' / name)
        xy = xy[t - t[0] <= WINDOW]
        measured.append(xy[::max(1, len(xy) // 12000)])
        if model == 'bench':
            states = simulate.kern.bench_trajectory(y, simulate.R0 + r, dt, round((SETTLE + WINDOW) / dt),
                                                    round(SETTLE / dt), parameters)
            portrait = states[::20, :2]
        else:
            _, states = simulate.integrate(g, ANIMATION_CIRCUIT, [r], y, SETTLE + WINDOW, dt, t_skip=SETTLE)
            portrait = states[::4, :2, 0]
        y = states[-1].copy()       # carry every state, including the op-amp outputs and the mean current
        assert np.isfinite(states).all()
        modeled.append(portrait.copy())
        if index % 20 == 0:
            print(f'Prepared {index + 1}/{len(rows)} portraits; Rpot={r:.2f} Ω', flush=True)
    return measured, modeled, dt


def shared_limits(measured, modeled):
    """Symmetric full-view limits that contain every portrait of both sources."""
    bounds = [np.ceil(max(np.max(np.abs(xy[:, axis])) for xy in measured + modeled) * 1.08)
              for axis in range(2)]
    full_limits = [(-bound, bound) for bound in bounds]
    print(f'Shared full-view limits: {full_limits}', flush=True)
    return full_limits


def write_metadata(rows, model, dt):
    """The frames' records and the model's parameters, beside the animations."""
    metadata = dict(direction='decreasing Rpot', frames=len(rows), fps=FPS,
                    portrait_window_s=WINDOW, model=model, settling_s=SETTLE,
                    dt_s=dt, continuation=True, r0_ohm=simulate.R0,
                    records=[dict(filename=n, rpot_ohm=r) for n, r in rows])
    if model == 'bench':
        parameters = simulate.bench_params()
        metadata.update(parameter_source='chua/identified.json',
                        parameters_SI=dict(zip(simulate.kern.BENCH_FIELDS, parameters.tolist())))
        name = 'resistance_sweep_bench_metadata.json'
    else:
        c = ANIMATION_CIRCUIT
        metadata.update(c1_nF=c.c1 * 1e9, c2_nF=c.c2 * 1e9, l_mH=c.l * 1e3, rL_ohm=c.rl)
        name = 'resistance_sweep_metadata.json'
    (OUT / name).write_text(json.dumps(metadata, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--model', choices=['ideal', 'bench'], default='ideal')
    parser.add_argument('--comparison-only', action='store_true',
                        help='render only the synchronized measurement/simulation comparison')
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    rows = records()
    measured, modeled, dt = portraits(rows, args.model)
    full_limits = shared_limits(measured, modeled)
    if not args.comparison_only:
        render(measured, rows, 'Measurements', '#167aab', full_limits)
        render(modeled, rows, 'Simulation' if args.model == 'ideal' else 'Bench simulation', '#cc6530', full_limits)
    render_comparison(measured, modeled, rows, full_limits, args.model)
    write_metadata(rows, args.model, dt)


if __name__ == '__main__':
    main()
