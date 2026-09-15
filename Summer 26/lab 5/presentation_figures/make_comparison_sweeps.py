"""Rebuild both comparison GIFs as 1920 x 1080 presentation slides.

python presentation_figures/make_comparison_sweeps.py
Uses the original metadata, simulation pipeline, frame order, and timing.
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
from PIL import Image

from make_measurements_sweep import BLUE, ORANGE, INK, ROOT, OUT, plt
from matplotlib.patches import Rectangle
import simulate
from scope_data import read_scope


def prepare(metadata):
    rows = metadata['records']
    bench = metadata['model'] == 'bench'
    dt = metadata['dt_s']
    settle, window = metadata['settling_s'], metadata['portrait_window_s']
    if bench:
        parameters = np.array([metadata['parameters_SI'][key] for key in simulate.kern.BENCH_FIELDS])
        g = simulate.bench_static_g(parameters)
        start, _ = simulate.sweep_starts(g, [rows[0]['rpot_ohm']], parameters[4])
        y = simulate.kern.bench_state(*start[:, 0], parameters)
    else:
        g = simulate.make_g()
        y, _ = simulate.sweep_starts(g, [rows[0]['rpot_ohm']], 0.)
    measured, modeled = [], []
    for index, row in enumerate(rows):
        t, xy = read_scope(ROOT / 'chua/forward' / row['filename'])
        xy = xy[t-t[0] <= window]
        measured.append(xy[::max(1, len(xy)//12000)])
        r = row['rpot_ohm']
        if bench:
            states = simulate.kern.bench_trajectory(y, metadata['r0_ohm'] + r, dt,
                round((settle + window)/dt), round(settle/dt), parameters)
            y = states[-1].copy()
            portrait = states[::20, :2]
        else:
            _, states = simulate.integrate([r], None, settle + window, dt, 0., g,
                                           t_skip=settle, y0=y)
            y = states[-1].copy()
            portrait = states[::4, :2, 0]
        assert np.isfinite(states).all()
        modeled.append(portrait.copy())
        if index % 20 == 0:
            print(f'{metadata["model"]}: prepared {index+1}/{len(rows)} portraits', flush=True)
    bounds = [np.ceil(max(np.max(np.abs(xy[:, axis])) for xy in measured+modeled)*1.08)
              for axis in range(2)]
    return measured, modeled, [(-bound, bound) for bound in bounds]


def render(model):
    name = 'comparison_equal_r_bench' if model == 'bench' else 'comparison_equal_r'
    original = ROOT / 'chua/animations' / (name + '.gif')
    backup = original.with_name(name + '_original.gif')
    metadata_name = 'resistance_sweep_bench_metadata.json' if model == 'bench' else 'resistance_sweep_metadata.json'
    metadata = json.loads((original.parent / metadata_name).read_text())
    rows = metadata['records']
    with Image.open(backup if backup.exists() else original) as gif:
        assert gif.n_frames == len(rows)
        durations = []
        for i in range(gif.n_frames):
            gif.seek(i)
            durations.append(gif.info['duration'])
        loop = gif.info.get('loop', 0)
    measured, modeled, full_limits = prepare(metadata)
    fig = plt.figure(figsize=(40/3, 7.5), dpi=144)
    fig.text(.075, .93, 'Measurements and simulation at equal resistance', fontsize=27, weight='bold')
    model_title = 'Model with measured components' if model == 'bench' else 'Ideal circuit model'
    fig.text(.075, .875, f'{model_title} · 20 ms portraits · matched voltage scales in each row',
             fontsize=15, color='#586673')
    label = fig.text(.075, .813, '', fontsize=19, weight='bold')
    source = fig.text(.95, .813, '', fontsize=10, color='#586673', ha='right')
    lines = []
    for col, (title, color) in enumerate([('Measurements', BLUE),
                                         ('Simulation · measured components' if model == 'bench' else 'Simulation · ideal model', ORANGE)]):
        column = []
        fig.text(.275 + col*.50, .753, title, fontsize=18, ha='center', color=color)
        for row, (limits, bottom) in enumerate([(full_limits, .535), (((-4, 4), (-1, 1)), .27)]):
            ax = fig.add_axes([.10+col*.50, bottom, .35, .19])
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
                        fontsize=10, color='#586673')
                ax.add_patch(Rectangle((-4, -1), 8, 2, fill=False,
                                      edgecolor='#85919b', linewidth=1.3, linestyle='--'))
            line, = ax.plot([], [], color=color, lw=.55, alpha=.75)
            column.append(line)
        lines.append(column)
    slider = fig.add_axes([.10, .143, .85, .02])
    slider.set(xlim=(rows[0]['rpot_ohm'], rows[-1]['rpot_ohm']), ylim=(-1, 1))
    slider.set_yticks([])
    slider.set_xticks([800, 700, 600, 500, 400, 300, 200, 100])
    slider.set_xlabel('Rpot (Ω) · decreasing resistance →', fontsize=15, labelpad=6)
    for spine in slider.spines.values():
        spine.set_visible(False)
    slider.axhline(0, color='#c8cfd5', lw=3)
    marker, = slider.plot([], [], 'o', color=INK, ms=10, clip_on=False)
    fig.text(.075, .035, 'Fixed scales · Dashed boxes mark the zoom · Model carries state and settles for 100 ms at each resistance',
             fontsize=10, color='#586673')
    frames = []
    still_index = min(range(len(rows)), key=lambda i: abs(rows[i]['rpot_ohm'] - 550))
    for i, row in enumerate(rows):
        for column, collection in zip(lines, [measured, modeled]):
            xy = collection[i]
            for line in column:
                line.set_data(xy[:, 0], xy[:, 1])
        r = row['rpot_ohm']
        label.set_text(f'Rpot = {r:.2f} Ω     |     Rtotal = {r+metadata["r0_ohm"]:.2f} Ω')
        source.set_text(f'forward/{row["filename"]}')
        marker.set_data([r], [0])
        fig.canvas.draw()
        frame = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert('RGB')
        frames.append(frame.quantize(colors=256, method=Image.Quantize.MEDIANCUT))
        if i == still_index:
            fig.savefig(OUT / (name + '.png'), dpi=144)
        if i % 20 == 0:
            print(f'{model}: rendered {i+1}/{len(rows)} frames', flush=True)
    path = OUT / (name + '.gif')
    frames[0].save(path, save_all=True, append_images=frames[1:], duration=durations,
                   loop=loop, disposal=2, optimize=False)
    with Image.open(path) as gif:
        assert gif.size == (1920, 1080) and gif.n_frames == len(rows)
        for i, duration in enumerate(durations):
            gif.seek(i)
            assert gif.info['duration'] == duration
    if not backup.exists():
        shutil.copy2(original, backup)
    shutil.copy2(path, original)
    shutil.copy2(OUT / (name + '.png'), original.with_suffix('.png'))
    plt.close(fig)
    print(f'Wrote {path}: {len(rows)} frames, {sum(durations)/1000:.2f} seconds', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', choices=['ideal', 'bench', 'both'], default='both')
    args = parser.parse_args()
    for model in (['ideal', 'bench'] if args.model == 'both' else [args.model]):
        render(model)
