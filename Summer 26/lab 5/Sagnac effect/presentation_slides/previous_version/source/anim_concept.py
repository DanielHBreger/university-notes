"""Slide 02: why rotation delays the co-rotating beam. Base slide PNG + region GIF.

Units: loop radius 1, light speed 1 (so light moves 1 rad of loop angle per time unit).
The table turns at OMEGA (rim speed OMEGA c), exaggerated for visibility.
"""
import numpy as np
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle
from PIL import Image

from deck import *

REGION = (96, 192, 864, 768)          # x, y, w, h in slide pixels: 0.5 in, 1.0 in, 4.5 in, 4.0 in
OMEGA = .25
T_MINUS = 2 * np.pi / (1 + OMEGA)     # counter-rotating pulse meets the splitter again
T_PLUS = 2 * np.pi / (1 - OMEGA)      # co-rotating pulse
START = -np.pi / 2                    # the splitter starts at the bottom of the loop
CY = 0.0                              # loop centre (data units) inside the frame


def arc(ax, th0, th1, r=1.0, **kw):
    th = np.linspace(th0, th1, max(2, int(abs(th1 - th0) * 40)))
    ax.plot(r * np.cos(th), CY + r * np.sin(th), **kw)


def frame(t, path):
    w, h = REGION[2], REGION[3]
    fig = plt.figure(figsize=(w / DPI, h / DPI), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    half_w = 1.77
    ax.set_xlim(-half_w, half_w)
    ax.set_ylim(1.30 - 2 * half_w * h / w, 1.30)
    ax.set_aspect('equal'); ax.axis('off')
    rot = OMEGA * t                                        # table angle turned so far
    # table (turns): disk and spokes
    ax.add_patch(Circle((0, CY), 1.12, fc='#f4f6f8', ec='#c9d1d8', lw=1.2, zorder=0))
    for k in range(6):
        a = rot + k * np.pi / 3 + np.pi / 6
        ax.plot([.18 * np.cos(a), 1.08 * np.cos(a)], [CY + .18 * np.sin(a), CY + 1.08 * np.sin(a)],
                color='#dfe5ea', lw=1.5, zorder=0)
    ax.add_patch(Circle((0, CY), 1.0, fill=False, ec='#55616c', lw=5, alpha=.35, zorder=1))   # fibre
    # rotation arrow (static, outside)
    th = np.linspace(np.radians(15), np.radians(60), 40)
    ax.plot(1.27 * np.cos(th), CY + 1.27 * np.sin(th), color=INK, lw=2)
    ax.annotate('', xy=(1.27 * np.cos(th[-1]), CY + 1.27 * np.sin(th[-1])),
                xytext=(1.27 * np.cos(th[-4]), CY + 1.27 * np.sin(th[-4])),
                arrowprops=dict(arrowstyle='-|>,head_length=.8,head_width=.4', color=INK, lw=2))
    ax.text(1.18, CY + 1.12, 'ω', fontsize=22, style='italic', ha='left')
    # where the splitter was at t = 0 (fixed in the lab)
    ax.plot([0, 0], [CY - 1.17, CY - .83], color=GREY, lw=1.3, ls=(0, (3, 2)), zorder=2)
    ax.text(0, CY - 1.22, 'start', fontsize=12, color=GREY, ha='center', va='top')
    # pulses
    s_ang = START + rot
    t_co = min(t, T_PLUS)
    t_ct = min(t, T_MINUS)
    a_co = START + t_co                     # anticlockwise, with the rotation
    a_ct = START - t_ct                     # clockwise, against it
    tail = .9
    if t < T_PLUS:
        arc(ax, a_co - min(tail, t_co), a_co, color=BLUE, lw=7, alpha=.35, solid_capstyle='round', zorder=3)
    if t < T_MINUS:
        arc(ax, a_ct, a_ct + min(tail, t_ct), color=ORANGE, lw=7, alpha=.35, solid_capstyle='round', zorder=3)
    for a, col, t_arr in [(a_co, BLUE, T_PLUS), (a_ct, ORANGE, T_MINUS)]:
        if t < t_arr:
            ax.add_patch(Circle((np.cos(a), CY + np.sin(a)), .065, fc=col, ec='white', lw=1.5, zorder=6))
        else:                                   # arrived: a fading ring where it met the splitter
            age = t - t_arr
            if age < 1.2 or t_arr == T_PLUS:
                r = .08 + .12 * min(age, 1.2) / 1.2
                ax.add_patch(Circle((np.cos(a), CY + np.sin(a)), r, fill=False, ec=col, lw=3, zorder=6,
                                    alpha=max(.15, 1 - age / 1.2) if t_arr == T_MINUS else 1))
    # splitter on the loop (turns with the table)
    sx, sy = np.cos(s_ang), CY + np.sin(s_ang)
    box = Rectangle((sx - .07, sy - .07), .14, .14, fc='white', ec=INK, lw=1.8, zorder=5,
                    transform=matplotlib.transforms.Affine2D().rotate_around(sx, sy, s_ang + np.pi / 2) + ax.transData)
    ax.add_patch(box)
    ax.text(sx + .2 * np.cos(s_ang), sy + .2 * np.sin(s_ang), 'splitter', fontsize=12, color=INK,
            ha='center', va='center', rotation=0)
    # labels in the middle
    ax.text(0, CY + .23, 'with rotation', color=BLUE, fontsize=16, ha='center', weight='bold')
    ax.text(0, CY + .02, 'against rotation', color=ORANGE, fontsize=16, ha='center', weight='bold')
    # timeline
    y0, x0, x1 = -1.50, -1.25, 1.25
    scale = (x1 - x0) / (T_PLUS * 1.05)
    ax.plot([x0, x1], [y0, y0], color='#b8c1c8', lw=2, solid_capstyle='round')
    ax.text(x0, y0 - .1, 'time of return\nto the splitter', fontsize=12, color=GREY, ha='left', va='top', linespacing=1.1)
    ax.plot([x0 + t * scale], [y0], 'o', ms=6, color=INK, zorder=4)
    if t >= T_MINUS:
        xm = x0 + T_MINUS * scale
        ax.plot([xm, xm], [y0 - .07, y0 + .07], color=ORANGE, lw=4)
    if t >= T_PLUS:
        xp = x0 + T_PLUS * scale
        ax.plot([xp, xp], [y0 - .07, y0 + .07], color=BLUE, lw=4)
        ax.annotate('', xy=(xp, y0 - .16), xytext=(x0 + T_MINUS * scale, y0 - .16),
                    arrowprops=dict(arrowstyle='<->', color=INK, lw=1.6, shrinkA=0, shrinkB=0))
        ax.text(.5 * (xp + x0 + T_MINUS * scale), y0 - .19, 'Δt', fontsize=18, ha='center', va='top',
                style='italic')
    fig.savefig(path, dpi=DPI, facecolor='white')
    plt.close(fig)


def base_slide():
    fig = canvas('Rotation delays the beam that travels with it',
                 'The principle of fibre-optic gyroscopes: the delay grows in proportion to the rotation rate')
    x = 1030 / 1920

    def at(y):
        return 1 - y / 1080
    fig.text(x, at(232), 'The co-rotating beam must chase the splitter;', fontsize=17, color=GREY)
    fig.text(x, at(265), 'the counter-rotating beam meets it early.', fontsize=17, color=GREY)
    fig.text(x, at(298), 'For N turns enclosing an area A:', fontsize=17, color=GREY)
    fig.text(x, at(405), r'$\Delta t = \dfrac{4NA}{c^2}\,\omega$', fontsize=30)
    fig.text(x, at(505), 'and the camera sees the phase difference', fontsize=17, color=GREY)
    fig.text(x, at(615), r'$\Delta\phi = \dfrac{2\pi c}{\lambda}\,\Delta t = \dfrac{8\pi NA}{\lambda c}\,\omega$',
             fontsize=30, color=BLUE)
    # scale box
    ax = px_axes(fig, 1020, 705, 800, 215)
    ax.axis('off')
    ax.add_patch(FancyBboxPatch((0.005, 0.02), .99, .96, boxstyle='round,pad=0,rounding_size=0.04',
                                transform=ax.transAxes, fc=LIGHT, ec='none'))
    dt = 4 * N_TURNS * AREA * SPEEDS.omega.iloc[-1] / C_DEF ** 2
    period = LAM / C_DEF
    ax.text(.04, .80, f'Our loop at ω = {SPEEDS.omega.iloc[-1]:.2f} rad/s', fontsize=16, color=GREY,
            transform=ax.transAxes, weight='bold')
    ax.text(.04, .50, f'Δt = {dt * 1e15:.2f} fs  ≈  ½ optical period ({period * 1e15:.2f} fs)',
            fontsize=18, transform=ax.transAxes)
    ax.text(.04, .20, f'→  Δφ ≈ {2 * np.pi * dt / period:.1f} rad, about half a fringe', fontsize=18,
            transform=ax.transAxes)
    footer(fig, 'Animation not to scale: the rim is drawn moving at 0.25 c and the loop as a single turn; in the '
                'experiment ωR ≈ 1 m/s ≈ 3 × 10⁻⁹ c.\nThe result does not depend on where the rotation axis is '
                'or on the fibre\'s refractive index (Schwartz, Eur. J. Phys. 38, 015301, 2017).')
    return save(fig, '02_sagnac_principle')


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        times = np.r_[np.zeros(1), np.linspace(0, T_PLUS, 100)]
        durations = [720] + [40] * 99 + [2600]
        paths = []
        for j, t in enumerate(times):
            p = tmp / f'f{j:04d}.png'
            frame(t, p)
            paths.append(p)
        # base slide with the final frame pasted into the animation region
        base = base_slide()
        im = Image.open(base).convert('RGB')
        last = Image.open(paths[-1]).convert('RGB')
        assert last.size == REGION[2:], last.size
        im.paste(last, REGION[:2])
        im.save(base)
        write_gif(paths, durations, OUT / '02_sagnac_principle_ANIMATION.gif', colors=128)


if __name__ == '__main__':
    main()
