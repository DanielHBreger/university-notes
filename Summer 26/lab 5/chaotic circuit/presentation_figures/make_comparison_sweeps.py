"""
Rebuild both comparison GIFs (measurements beside the ideal model, and beside the
identified circuit) as 1920 x 1080 presentation slides.

Uses the original sweep metadata in chua/animations, the continuation simulation
of animate_resistance.py, and the original frame order and timing. The GIFs and
stills are written here and replace their counterparts in chua/animations; the
first version of each GIF there is kept as *_original.gif.

Run from the repository root:
    python presentation_figures/make_comparison_sweeps.py [--model ideal|bench|both]
"""
import argparse
import json

import numpy as np
from matplotlib.patches import Rectangle

from make_measurements_sweep import (BLUE, ORANGE, INK, GREY, ROOT, OUT, ANIMATIONS, plt, add_slider,
                                     frame_image, original_timing, replace_original, still_index,
                                     write_gif)
from animate_resistance import ANIMATION_CIRCUIT
from scope_data import read_scope
import simulate


def prepare(metadata):
    """(measured, modelled, full-view limits): the portraits of every frame, the model carrying its state along."""
    rows = metadata['records']
    bench = metadata['model'] == 'bench'
    dt = metadata['dt_s']
    settle, window = metadata['settling_s'], metadata['portrait_window_s']
    first = rows[0]['rpot_ohm']
    if bench:
        # fields added after the metadata was written: tauM (27 September), the node-1 law (off: constant C1)
        stored = {'tauM': simulate.kern.TAU_M, 'c1_a0': 0.0, 'c1_a1': 0.0, 'c1_xk': 0.5, **metadata['parameters_SI']}
        parameters = np.array([stored[key] for key in simulate.kern.BENCH_FIELDS])
        start, _ = simulate.sweep_starts(simulate.bench_static_g(parameters),
                                         simulate.bench_small_signal(parameters), [first])
        y = simulate.kern.bench_state(*start[:, 0], parameters)
    else:
        g = simulate.make_g()
        y, _ = simulate.sweep_starts(g, ANIMATION_CIRCUIT, [first])
    measured, modeled = [], []
    for index, row in enumerate(rows):
        t, xy = read_scope(ROOT / 'chua/forward' / row['filename'])
        xy = xy[t - t[0] <= window]
        measured.append(xy[::max(1, len(xy) // 12000)])
        r = row['rpot_ohm']
        if bench:
            states = simulate.kern.bench_trajectory(y, metadata['r0_ohm'] + r, dt, round((settle + window) / dt),
                                                    round(settle / dt), parameters)
            portrait = states[::20, :2]
        else:
            _, states = simulate.integrate(g, ANIMATION_CIRCUIT, [r], y, settle + window, dt, t_skip=settle)
            portrait = states[::4, :2, 0]
        y = states[-1].copy()
        assert np.isfinite(states).all()
        modeled.append(portrait.copy())
        if index % 20 == 0:
            print(f'{metadata["model"]}: prepared {index + 1}/{len(rows)} portraits', flush=True)
    bounds = [np.ceil(max(np.max(np.abs(xy[:, axis])) for xy in measured + modeled) * 1.08)
              for axis in range(2)]
    return measured, modeled, [(-bound, bound) for bound in bounds]


def portrait_axes(fig, col, row, limits, color):
    """One portrait panel: the full view (row 0, with the zoom box) or the zoom (row 1); returns its line."""
    ax = fig.add_axes([.10 + col * .50, (.535, .27)[row], .35, .19])
    ax.set(xlim=limits[0], ylim=limits[1], ylabel='V2 (V)')
    ax.set_xlabel('V1 (V)' if row else '', fontsize=16)
    ax.tick_params(labelsize=14)
    ax.yaxis.label.set_size(16)
    ax.grid(alpha=.15)
    ax.set_axisbelow(True)
    ax.set_xticks([-4, -2, 0, 2, 4] if row else [-8, -4, 0, 4, 8])
    ax.set_yticks([-1, 0, 1] if row else [-4, 0, 4])
    if row:
        ax.set_title('Zoom · small attractors', fontsize=15, pad=8)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('#85919b')
            spine.set_linestyle('--')
    else:
        ax.text(.5, 1.02, 'Full view', transform=ax.transAxes, ha='center', va='bottom',
                fontsize=10, color=GREY)
        ax.add_patch(Rectangle((-4, -1), 8, 2, fill=False,
                               edgecolor='#85919b', linewidth=1.3, linestyle='--'))
    line, = ax.plot([], [], color=color, lw=.55, alpha=.75)
    return line


def render(model):
    """Rebuild one comparison animation ('ideal' or 'bench')."""
    name = 'comparison_equal_r_bench' if model == 'bench' else 'comparison_equal_r'
    original = ANIMATIONS / (name + '.gif')
    metadata_name = 'resistance_sweep_bench_metadata.json' if model == 'bench' else 'resistance_sweep_metadata.json'
    metadata = json.loads((ANIMATIONS / metadata_name).read_text())
    rows = metadata['records']
    durations, loop = original_timing(original, len(rows))
    measured, modeled, full_limits = prepare(metadata)

    fig = plt.figure(figsize=(40 / 3, 7.5), dpi=144)
    fig.text(.075, .93, 'Measurements and simulation at equal resistance', fontsize=27, weight='bold')
    model_title = 'Model with measured components' if model == 'bench' else 'Ideal circuit model'
    fig.text(.075, .875, f'{model_title} · 20 ms portraits · matched voltage scales in each row',
             fontsize=15, color=GREY)
    label = fig.text(.075, .813, '', fontsize=19, weight='bold')
    source = fig.text(.95, .813, '', fontsize=10, color=GREY, ha='right')
    simulation_title = 'Simulation · measured components' if model == 'bench' else 'Simulation · ideal model'
    lines = []
    for col, (title, color) in enumerate([('Measurements', BLUE), (simulation_title, ORANGE)]):
        fig.text(.275 + col * .50, .753, title, fontsize=18, ha='center', color=color)
        lines.append([portrait_axes(fig, col, row, limits, color)
                      for row, limits in enumerate([full_limits, ((-4, 4), (-1, 1))])])
    marker = add_slider(fig, rows, INK)
    fig.text(.075, .035, 'Fixed scales · Dashed boxes mark the zoom · '
                         'Model carries state and settles for 100 ms at each resistance',
             fontsize=10, color=GREY)

    still = OUT / (name + '.png')
    still_at = still_index(rows)
    frames = []
    for i, row in enumerate(rows):
        for column, collection in zip(lines, [measured, modeled]):
            xy = collection[i]
            for line in column:
                line.set_data(xy[:, 0], xy[:, 1])
        r = row['rpot_ohm']
        label.set_text(f'Rpot = {r:.2f} Ω     |     Rtotal = {r + metadata["r0_ohm"]:.2f} Ω')
        source.set_text(f'forward/{row["filename"]}')
        marker.set_data([r], [0])
        frames.append(frame_image(fig))
        if i == still_at:
            fig.savefig(still, dpi=144)
        if i % 20 == 0:
            print(f'{model}: rendered {i + 1}/{len(rows)} frames', flush=True)
    path = OUT / (name + '.gif')
    write_gif(frames, path, durations, loop)
    replace_original(path, still, original)
    plt.close(fig)
    print(f'Wrote {path}: {len(rows)} frames, {sum(durations) / 1000:.2f} seconds', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--model', choices=['ideal', 'bench', 'both'], default='both')
    args = parser.parse_args()
    for model in (['ideal', 'bench'] if args.model == 'both' else [args.model]):
        render(model)


if __name__ == '__main__':
    main()
