"""
Render the original measured resistance sweep as a full 1920 x 1080 slide.

Uses the original selection of records (chua/animations/resistance_sweep_metadata.json),
raw 20 ms portraits and the original GIF's frame timing. The GIF and its still are
written here and replace their counterparts in chua/animations; the first version of
the GIF there is kept as *_original.gif.

Run from the repository root: python presentation_figures/make_measurements_sweep.py
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
from scope_data import read_scope  # noqa: E402

OUT = ROOT / 'presentation_figures'
ANIMATIONS = ROOT / 'chua/animations'
BLUE, ORANGE, INK = '#147da3', '#d87530', '#243444'
GREY = '#586673'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 17,
                     'axes.labelsize': 18, 'axes.titlesize': 20, 'xtick.labelsize': 15,
                     'ytick.labelsize': 15, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.edgecolor': '#85919b', 'text.color': INK, 'axes.labelcolor': INK,
                     'xtick.color': INK, 'ytick.color': INK, 'savefig.facecolor': 'white'})


# ---- shared with make_comparison_sweeps.py -------------------------------------------
def backup_of(original):
    """Where the first version of a rebuilt animation is kept."""
    return original.with_name(original.stem + '_original.gif')


def original_timing(original, n_frames):
    """(frame durations in ms, loop count) of the animation being rebuilt, from its first version."""
    backup = backup_of(original)
    with Image.open(backup if backup.exists() else original) as gif:
        assert gif.n_frames == n_frames
        durations = []
        for i in range(gif.n_frames):
            gif.seek(i)
            durations.append(gif.info['duration'])
        return durations, gif.info.get('loop', 0)


def add_slider(fig, rows, color):
    """The Rpot indicator along the bottom of a slide; returns its marker."""
    slider = fig.add_axes([.10, .143, .85, .02])
    slider.set(xlim=(rows[0]['rpot_ohm'], rows[-1]['rpot_ohm']), ylim=(-1, 1))
    slider.set_yticks([])
    slider.set_xticks([800, 700, 600, 500, 400, 300, 200, 100])
    slider.set_xlabel('Rpot (Ω) · decreasing resistance →', fontsize=15, labelpad=6)
    for spine in slider.spines.values():
        spine.set_visible(False)
    slider.axhline(0, color='#c8cfd5', lw=3)
    marker, = slider.plot([], [], 'o', color=color, ms=10, clip_on=False)
    return marker


def frame_image(fig):
    """The figure as a 256-colour GIF frame."""
    fig.canvas.draw()
    frame = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert('RGB')
    return frame.quantize(colors=256, method=Image.Quantize.MEDIANCUT)


def write_gif(frames, path, durations, loop):
    """Save the frames with the original timing and check the result."""
    frames[0].save(path, save_all=True, append_images=frames[1:], duration=durations,
                   loop=loop, disposal=2, optimize=False)
    with Image.open(path) as result:
        assert result.size == (1920, 1080) and result.n_frames == len(frames)
        for i, duration in enumerate(durations):
            result.seek(i)
            assert result.info['duration'] == duration


def replace_original(gif, still, original):
    """Copy the new GIF and still over the animation in chua/animations, keeping its first version."""
    backup = backup_of(original)
    if not backup.exists():
        shutil.copy2(original, backup)
    shutil.copy2(gif, original)
    shutil.copy2(still, original.with_suffix('.png'))


def still_index(rows, rpot=550):
    """The frame nearest `rpot`, saved as the still."""
    return min(range(len(rows)), key=lambda i: abs(rows[i]['rpot_ohm'] - rpot))


# ---- the measured sweep ------------------------------------------------------------
def main():
    metadata = json.loads((ANIMATIONS / 'resistance_sweep_metadata.json').read_text())
    rows = metadata['records']
    original = ANIMATIONS / 'measurements_resistance_sweep.gif'
    durations, loop = original_timing(original, len(rows))

    fig = plt.figure(figsize=(40 / 3, 7.5), dpi=144)
    fig.text(.075, .93, 'Measured regimes as resistance decreases', fontsize=27, weight='bold')
    fig.text(.075, .875, 'V1 and V2 are the voltages across C1 and C2; each frame shows 20 ms',
             fontsize=15, color=GREY)
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
    marker = add_slider(fig, rows, BLUE)
    fig.text(.075, .035, 'Fixed voltage scales · Orange box marks the zoom · Playback advances resistance, not trajectory time',
             fontsize=10, color=GREY)
    source = fig.text(.95, .813, '', fontsize=10, color=GREY, ha='right')

    still = OUT / 'measurements_resistance_sweep.png'
    still_at = still_index(rows)
    frames = []
    for i, row in enumerate(rows):
        t, xy = read_scope(ROOT / 'chua/forward' / row['filename'])
        xy = xy[t - t[0] <= metadata['portrait_window_s']]
        xy = xy[::max(1, len(xy) // 12000)]
        assert np.isfinite(xy).all() and np.max(np.abs(xy)) < 9, 'Full view clips data'
        for line in lines:
            line.set_data(xy[:, 0], xy[:, 1])
        r = row['rpot_ohm']
        label.set_text(f'Rpot = {r:.2f} Ω     |     Rtotal = {r + metadata["r0_ohm"]:.2f} Ω')
        source.set_text(f'forward/{row["filename"]}')
        marker.set_data([r], [0])
        frames.append(frame_image(fig))
        if i == still_at:
            fig.savefig(still, dpi=144)
        if i % 20 == 0:
            print(f'Rendered {i + 1}/{len(rows)} frames', flush=True)
    path = OUT / 'measurements_resistance_sweep.gif'
    write_gif(frames, path, durations, loop)
    replace_original(path, still, original)
    plt.close(fig)
    print(f'Wrote {path}; {len(rows)} frames; {sum(durations) / 1000:.2f} seconds', flush=True)


if __name__ == '__main__':
    main()
