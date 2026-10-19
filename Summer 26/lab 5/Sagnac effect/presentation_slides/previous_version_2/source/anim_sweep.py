"""Slide 06: all 161 photographs in recording order, with the fitted rings and the phase building up.

Writes out/06_rotation_sweep.png (base slide showing the final frame) and
out/06_rotation_sweep_ANIMATION.gif (the animated region only).
"""
import numpy as np
from matplotlib.patches import Circle, FancyBboxPatch
from PIL import Image

from deck import *
from matplotlib.ticker import MultipleLocator

REGION = (96, 192, 1728, 768)          # 0.5 in, 1.0 in, 9.0 in x 4.0 in
HALF = 500                             # photo crop half-width (px of the camera image)
SETTINGS = SPEEDS.S.to_numpy()
FR = FRAMES.sort_values('time_s').reset_index(drop=True)
FIT_A, FIT_B = RES['slope'], RES['intercept'] - SPEEDS.phi.iloc[0]   # line in the plotted (rest-subtracted) frame
TH = np.linspace(0, 2 * np.pi, 300)


def display_photo(row, size=600, sigma=1.8, levels=48):
    """Crop resampled to the display size, lightly smoothed and posterised (keeps the GIF small)."""
    from scipy import ndimage
    img, _ = crop(photo(row.file), row.cx, row.cy, HALF)
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


def frame(k, path, show_fit=False):
    """Frame k (0..160) of the sweep; show_fit adds the fitted line and the theory band."""
    w, h = REGION[2], REGION[3]
    fig = plt.figure(figsize=(w / DPI, h / DPI), dpi=DPI)
    row = FR.iloc[k]
    n_in_setting = (FR.S == row.S).sum()
    j_in_setting = int((FR.iloc[:k + 1].S == row.S).sum())
    # ---- counters
    fig.text(48 / w, 1 - 42 / h, f'ω = {row.omega:.3f} rad/s', fontsize=22, weight='bold')
    fig.text(462 / w, 1 - 42 / h, '|', fontsize=22, color='#b0b9c1')
    fig.text(495 / w, 1 - 42 / h, f'Δφ = {row.dphi:+.2f} rad', fontsize=22, weight='bold', color=BLUE)
    fig.text(897 / w, 1 - 42 / h, f'({row.dphi / (2 * np.pi):+.2f} fringe)', fontsize=17, color=BLUE)
    fig.text(1 - 20 / w, 1 - 28 / h, (f'motor setting S = {row.S:g}' if row.S else 'at rest (S = 0)') + f'  ·  photo {j_in_setting} of {n_in_setting}',
             fontsize=11, color=GREY, ha='right')
    fig.text(1 - 20 / w, 1 - 50 / h, row.file, fontsize=10, color=GREY, ha='right')
    # ---- photo with the fitted rings
    ax = fig.add_axes([48 / w, 1 - 680 / h, 600 / w, 600 / h])
    ax.imshow(display_photo(row), extent=(-HALF, HALF, HALF, -HALF), interpolation='none')
    for R in rings(row):
        ax.plot(R * np.cos(TH), R * np.sin(TH), color='white', lw=1.5, ls=(0, (5, 3)), alpha=.85)
    ax.plot(0, 0, '+', color='white', ms=13, mew=1.6)
    ax.set_xlim(-HALF, HALF); ax.set_ylim(HALF, -HALF); ax.axis('off')
    # ---- phase against omega
    pa = fig.add_axes([815 / w, 1 - 590 / h, 805 / w, 500 / h])
    done = FR.iloc[:k + 1]
    pa.scatter(done.omega, done.dphi, s=16, color=BLUE, alpha=.28, lw=0, zorder=2)
    finished = [S for S in SETTINGS if (done.S == S).sum() == (FR.S == S).sum()]
    sp = SPEEDS[SPEEDS.S.isin(finished)]
    pa.errorbar(sp.omega, sp.dphi, sp.u_fit, fmt='o', color=BLUE, ms=9, mec='white', mew=1.2, capsize=0,
                elinewidth=2, zorder=4)
    pa.plot([row.omega], [row.dphi], 'o', ms=15, mfc='none', mec=INK, mew=2.2, zorder=5)
    xs = np.linspace(0, 1.9, 50)
    if show_fit:
        th = RES['theory'] * xs
        pa.fill_between(xs, (RES['theory'] - RES['u_theory']) * xs, (RES['theory'] + RES['u_theory']) * xs,
                        color=THEORY, alpha=.16, lw=0, zorder=1)
        pa.plot(xs, th, color=THEORY, lw=2, ls=(0, (6, 3)), zorder=1)
        pa.plot(xs, FIT_A * xs + FIT_B, color=ORANGE, lw=2.6, zorder=3)
        pa.text(.03, .93, f'fitted slope: {RES["slope"]:.2f} ± {RES["u_slope"]:.2f} rad/(rad/s)', color=ORANGE,
                fontsize=16, transform=pa.transAxes, weight='bold')
        pa.text(.03, .84, f'Sagnac prediction: {RES["theory"]:.2f} ± {RES["u_theory"]:.2f}', color=THEORY,
                fontsize=16, transform=pa.transAxes, weight='bold')
    pa.set(xlim=(-.05, 1.92), ylim=(-.35, 3.65), xlabel=r'Angular velocity $\omega$ (rad/s)', ylabel=r'Phase shift $\Delta\phi$ (rad)')
    style(pa)
    pr = pa.twinx()
    pr.set_ylim(np.array(pa.get_ylim()) / (2 * np.pi))
    pr.set_ylabel('fringes', color=GREY)
    pr.tick_params(colors=GREY)
    pr.spines['right'].set_visible(True); pr.spines['right'].set_color('#85919b')
    pr.yaxis.set_major_locator(MultipleLocator(.1))
    # ---- recording clock
    t0, t1 = FR.time_s.min(), FR.time_s.max()
    span = (t1 - t0) / 60
    ta = fig.add_axes([62 / w, 1 - 722 / h, 1420 / w, 22 / h])
    ta.set_xlim(-.05, span + .05)
    ta.set_ylim(-1, 1)
    ta.plot([0, span], [0, 0], color='#d3d9de', lw=6, solid_capstyle='round')
    ta.plot([0, (row.time_s - t0) / 60], [0, 0], color='#9fb7c4', lw=6, solid_capstyle='round')
    ta.plot([(row.time_s - t0) / 60], [0], 'o', ms=13, color=INK, clip_on=False)
    ta.set_yticks([])
    for s in ['left', 'bottom']:
        ta.spines[s].set_visible(False)
    ta.set_xticks(np.arange(0, 13, 2))
    ta.set_xticklabels([f'{m} min' for m in range(0, 13, 2)], fontsize=12, color=GREY)
    ta.tick_params(length=0, pad=6)
    fig.text(1 - 20 / w, 1 - 716 / h, clock(row.time_s), fontsize=17, color=INK, ha='right', va='bottom',
             weight='bold')
    fig.text(1 - 20 / w, 1 - 742 / h, 'recording clock', fontsize=11, color=GREY, ha='right', va='bottom')
    fig.savefig(path, dpi=DPI, facecolor='white')
    plt.close(fig)


def base_slide():
    fig = canvas('Faster rotation pulls the rings inward',
                 'All 161 photographs in recording order, with the fitted dark rings dashed')
    footer(fig, r'$\omega$: linear calibration at each photo’s motor setting. $\Delta\phi$: phase at the ring centre '
                'relative to the mean at rest (method on the next slide).\nSmall dots are single photos, large dots the '
                'mean per setting. Photos lightly smoothed for display. The bar shows the recording clock.')
    return save(fig, '06_rotation_sweep')


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        paths, durations = [], []
        for k in range(len(FR)):
            p = tmp / f'f{k:04d}.png'
            frame(k, p)
            paths.append(p)
            durations.append(900 if k == 0 else 100)
        p = tmp / 'final.png'
        frame(len(FR) - 1, p, show_fit=True)
        paths.append(p); durations.append(3500)
        base = base_slide()
        im = Image.open(base).convert('RGB')
        last = Image.open(paths[-1]).convert('RGB')
        assert last.size == REGION[2:], last.size
        im.paste(last, REGION[:2])
        im.save(base)
        write_gif(paths, durations, OUT / '06_rotation_sweep_ANIMATION.gif')


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:          # quick preview of single frames
        for k in sys.argv[1:]:
            frame(int(k), OUT / f'_preview_sweep_{k}.png', show_fit=int(k) == 160)
    else:
        main()
