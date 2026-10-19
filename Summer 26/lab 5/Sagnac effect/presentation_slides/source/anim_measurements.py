"""Slide 05: the measurements. All 161 photographs in recording order, raw, with the rate of rotation.

Writes out/05_measurements.png (the slide, with a guide box where the animation goes) and
out/05_measurements_ANIMATION.gif (the animated region only). No analysis is shown here: the
photographs are cropped around a fixed point (the ring centre at rest) and the only overlay is a
reference circle, the first dark ring at rest.
"""
import numpy as np
from matplotlib.ticker import MultipleLocator
from PIL import Image

from deck import *

REGION = (96, 192, 1728, 768)          # 0.5 in, 1.0 in, 9.0 in x 4.0 in
HALF = 500                             # photo crop half-width (px of the camera image)
FR = FRAMES.sort_values('time_s').reset_index(drop=True)
REST = FR[FR.S == 0]
CX0, CY0 = REST.cx.median(), REST.cy.median()       # fixed crop centre: the ring centre at rest
R1_REST = REST.R1.median()                          # first dark ring at rest (px)
TH = np.linspace(0, 2 * np.pi, 300)
T0 = FR.time_s.min()
MIN = (FR.time_s - T0) / 60                         # minutes since the first photograph


def display_photo(row, size=600, sigma=1.8, levels=48):
    """Fixed crop resampled to the display size, lightly smoothed and posterised (keeps the GIF small)."""
    from scipy import ndimage
    img, _ = crop(photo(row.file), CX0, CY0, HALF)
    red = np.asarray(Image.fromarray(img[..., 0]).resize((size, size), Image.LANCZOS), float)
    red = ndimage.gaussian_filter(red, sigma)
    step = 256 / levels
    red = np.floor(red / step) * step + step / 2
    rgb = np.zeros((size, size, 3), np.uint8)
    rgb[..., 0] = np.clip(red, 0, 255)
    return rgb


def clock(t):
    h, m, s = int(t // 3600), int(t % 3600 // 60), int(t % 60)
    return f'{h:02d}:{m:02d}:{s:02d}'


def staircase(ax, upto=None, **kw):
    """omega against minutes as a step line, up to photograph index upto (inclusive)."""
    sel = FR.index if upto is None else FR.index[:upto + 1]
    x, y = MIN[sel].to_numpy(), FR.omega[sel].to_numpy()
    ax.step(x, y, where='post', **kw)


def frame(k, path):
    w, h = REGION[2], REGION[3]
    fig = plt.figure(figsize=(w / DPI, h / DPI), dpi=DPI)
    row = FR.iloc[k]
    n_in = int((FR.S == row.S).sum())
    j_in = int((FR.iloc[:k + 1].S == row.S).sum())
    # ---- counters
    fig.text(48 / w, 1 - 42 / h, f'ω = {row.omega:.3f} rad/s', fontsize=22, weight='bold')
    fig.text(500 / w, 1 - 40 / h, ('at rest' if row.S == 0 else f'motor setting S = {row.S:g}') +
             f'   ·   photo {j_in} of {n_in}', fontsize=15, color=GREY)
    fig.text(1 - 20 / w, 1 - 42 / h, clock(row.time_s), fontsize=20, weight='bold', ha='right')
    # ---- photo, fixed crop, with the first dark ring at rest as a reference
    ax = fig.add_axes([48 / w, 1 - 680 / h, 600 / w, 600 / h])
    ax.imshow(display_photo(row), extent=(-HALF, HALF, HALF, -HALF), interpolation='none')
    ax.plot(R1_REST * np.cos(TH), R1_REST * np.sin(TH), color='white', lw=1.6, ls=(0, (5, 3)), alpha=.9)
    ax.text(0, -R1_REST - 14, 'first dark ring at rest', color='white', fontsize=12, ha='center', va='bottom',
            bbox=dict(boxstyle='round,pad=0.2', fc='black', ec='none', alpha=.55))
    ax.set_xlim(-HALF, HALF); ax.set_ylim(HALF, -HALF); ax.axis('off')
    # ---- rate of rotation during the run
    pa = fig.add_axes([815 / w, 1 - 650 / h, 860 / w, 520 / h])
    staircase(pa, color='#c9d1d8', lw=2.2)                                    # the whole run, light
    staircase(pa, upto=k, color=INK, lw=2.6)                                  # recorded so far
    pa.plot(MIN[:k + 1], FR.omega[:k + 1], 'o', color=BLUE, ms=4, zorder=3)
    pa.plot([MIN[k]], [row.omega], 'o', ms=15, mfc='none', mec=BLUE, mew=2.4, zorder=4)
    pa.set(xlim=(-.3, 12.5), ylim=(-.1, 2.0), xlabel='Minutes since the first photograph',
           ylabel=r'Angular velocity $\omega$ (rad/s)')
    pa.yaxis.set_major_locator(MultipleLocator(.5))
    style(pa)
    pa.set_title('Rate of rotation during the run', loc='left', fontsize=16, color=GREY, pad=10)
    fig.savefig(path, dpi=DPI, facecolor='white')
    plt.close(fig)


def base_slide():
    fig = canvas('The rings move inward as the speed is stepped up',
                 'All 161 photographs in recording order: at rest, then nine speeds in increasing order')
    placeholder(fig, REGION, '05_measurements_ANIMATION.gif')
    footer(fig, r'$\omega$ of each photo from the calibration on the previous slide. Dashed circle: the first dark ring '
                'at rest. Photos lightly smoothed for display.')
    return save(fig, '05_measurements')


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        paths, durations = [], []
        for k in range(len(FR)):
            p = tmp / f'f{k:04d}.png'
            frame(k, p)
            paths.append(p)
            durations.append(1000 if k == 0 else (2500 if k == len(FR) - 1 else 100))
        assert Image.open(paths[-1]).size == REGION[2:]
        base_slide()
        write_gif(paths, durations, OUT / '05_measurements_ANIMATION.gif')


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:          # quick preview of single frames
        for k in sys.argv[1:]:
            frame(int(k), OUT / f'_preview_meas_{k}.png')
    else:
        main()
