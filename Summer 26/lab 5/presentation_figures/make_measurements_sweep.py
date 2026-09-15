"""Render the original measured resistance sweep as a full 1920 x 1080 slide.

Run from the repository root: python presentation_figures/make_measurements_sweep.py
Uses the original selection, raw 20 ms portraits, and GIF frame timing.
"""
import json
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'chua'))
from scope_data import read_scope

OUT = ROOT / 'presentation_figures'
ORIGINAL = ROOT / 'chua/animations/measurements_resistance_sweep.gif'
BACKUP = ORIGINAL.with_name('measurements_resistance_sweep_original.gif')
BLUE, ORANGE, INK = '#147da3', '#d87530', '#243444'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 17,
    'axes.labelsize': 18, 'axes.titlesize': 20, 'xtick.labelsize': 15,
    'ytick.labelsize': 15, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.edgecolor': '#85919b', 'text.color': INK, 'axes.labelcolor': INK,
    'xtick.color': INK, 'ytick.color': INK, 'savefig.facecolor': 'white'})


def main():
    metadata = json.loads((ROOT / 'chua/animations/resistance_sweep_metadata.json').read_text())
    rows = metadata['records']
    with Image.open(BACKUP if BACKUP.exists() else ORIGINAL) as original:
        assert original.n_frames == len(rows)
        durations = []
        for i in range(original.n_frames):
            original.seek(i)
            durations.append(original.info['duration'])
        loop = original.info.get('loop', 0)

    fig = plt.figure(figsize=(40 / 3, 7.5), dpi=144)
    fig.text(.075, .93, 'Measured regimes as resistance decreases', fontsize=27, weight='bold')
    fig.text(.075, .875, 'V1 and V2 are the voltages across C1 and C2; each frame shows 20 ms',
             fontsize=15, color='#586673')
    label = fig.text(.075, .813, '', fontsize=19, weight='bold')
    axes = [fig.add_axes([.10, .285, .35, .435]),
            fig.add_axes([.60, .285, .35, .435])]
    lines = []
    for ax, limits, title in zip(axes, [((-9, 9), (-9, 9)), ((-4, 4), (-1, 1))],
                                  ['Full view', 'Zoom · small attractors']):
        ax.set(xlim=limits[0], ylim=limits[1], xlabel='V1 (V)', ylabel='V2 (V)', title=title)
        ax.grid(alpha=.15)
        ax.set_axisbelow(True)
        line, = ax.plot([], [], color=BLUE, lw=.55, alpha=.75)
        lines.append(line)
    axes[0].set_xticks([-8, -4, 0, 4, 8])
    axes[0].set_yticks([-8, -4, 0, 4, 8])
    axes[1].set_xticks([-4, -2, 0, 2, 4])
    axes[1].set_yticks([-1, -.5, 0, .5, 1])
    axes[0].add_patch(Rectangle((-4, -1), 8, 2, fill=False,
                               edgecolor=ORANGE, linewidth=1.8, linestyle='--'))
    for spine in axes[1].spines.values():
        spine.set_visible(True)
        spine.set_color(ORANGE)
        spine.set_linestyle('--')
        spine.set_linewidth(1.2)
    slider = fig.add_axes([.10, .143, .85, .02])
    slider.set(xlim=(rows[0]['rpot_ohm'], rows[-1]['rpot_ohm']), ylim=(-1, 1))
    slider.set_yticks([])
    slider.set_xticks([800, 700, 600, 500, 400, 300, 200, 100])
    slider.set_xlabel('Rpot (Ω) · decreasing resistance →', fontsize=15, labelpad=6)
    for spine in slider.spines.values():
        spine.set_visible(False)
    slider.axhline(0, color='#c8cfd5', lw=3)
    marker, = slider.plot([], [], 'o', color=BLUE, ms=10, clip_on=False)
    fig.text(.075, .035, 'Fixed voltage scales · Orange box marks the zoom · Playback advances resistance, not trajectory time',
             fontsize=10, color='#586673')
    source = fig.text(.95, .813, '', fontsize=10, color='#586673', ha='right')
    frames = []
    still_index = min(range(len(rows)), key=lambda i: abs(rows[i]['rpot_ohm'] - 550))
    for i, row in enumerate(rows):
        t, xy = read_scope(ROOT / 'chua/forward' / row['filename'])
        xy = xy[t-t[0] <= metadata['portrait_window_s']]
        xy = xy[::max(1, len(xy)//12000)]
        assert np.isfinite(xy).all() and np.max(np.abs(xy)) < 9, 'Full view clips data'
        for line in lines:
            line.set_data(xy[:, 0], xy[:, 1])
        r = row['rpot_ohm']
        label.set_text(f'Rpot = {r:.2f} Ω     |     Rtotal = {r + metadata["r0_ohm"]:.2f} Ω')
        source.set_text(f'forward/{row["filename"]}')
        marker.set_data([r], [0])
        fig.canvas.draw()
        frame = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert('RGB')
        frames.append(frame.quantize(colors=256, method=Image.Quantize.MEDIANCUT))
        if i == still_index:
            fig.savefig(OUT / 'measurements_resistance_sweep.png', dpi=144)
        if i % 20 == 0:
            print(f'Rendered {i+1}/{len(rows)} frames', flush=True)
    path = OUT / 'measurements_resistance_sweep.gif'
    frames[0].save(path, save_all=True, append_images=frames[1:], duration=durations,
                   loop=loop, disposal=2, optimize=False)
    with Image.open(path) as result:
        assert result.size == (1920, 1080) and result.n_frames == len(rows)
        for i, duration in enumerate(durations):
            result.seek(i)
            assert result.info['duration'] == duration
    if not BACKUP.exists():
        shutil.copy2(ORIGINAL, BACKUP)
    shutil.copy2(path, ORIGINAL)
    shutil.copy2(OUT / 'measurements_resistance_sweep.png', ORIGINAL.with_suffix('.png'))
    plt.close(fig)
    print(f'Wrote {path}; {len(rows)} frames; {sum(durations)/1000:.2f} seconds', flush=True)


if __name__ == '__main__':
    main()
