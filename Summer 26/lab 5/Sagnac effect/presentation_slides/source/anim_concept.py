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
    ax.add_patch(Circle((0, CY), 1.0, fill=False, ec='#55616c', lw=5, alpha=.35, zorder=1))   # fiber
    # rotation arrow (static, outside)
    th = np.linspace(np.radians(15), np.radians(60), 40)
    ax.plot(1.27 * np.cos(th), CY + 1.27 * np.sin(th), color=INK, lw=2)
    ax.annotate('', xy=(1.27 * np.cos(th[-1]), CY + 1.27 * np.sin(th[-1])),
                xytext=(1.27 * np.cos(th[-4]), CY + 1.27 * np.sin(th[-4])),
                arrowprops=dict(arrowstyle='-|>,head_length=.8,head_width=.4', color=INK, lw=2))
    ax.text(1.18, CY + 1.12, r'$\omega$', fontsize=24, ha='left')
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
    cs, sn = np.cos(s_ang), np.sin(s_ang)              # label just outside the box, never on top of it
    if sn < -.8:                                       # near the bottom: beside the box, clear of 'start'
        ax.text(sx + .12, sy + .02, 'splitter', fontsize=12, color=INK, zorder=7, ha='left', va='bottom')
    else:
        ax.text(sx + .13 * cs, sy + .13 * sn, 'splitter', fontsize=12, color=INK, rotation=0, zorder=7,
                ha='left' if cs > .35 else ('right' if cs < -.35 else 'center'),
                va='bottom' if sn > .35 else ('top' if sn < -.35 else 'center'))
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
        ax.text(.5 * (xp + x0 + T_MINUS * scale), y0 - .19, r'$\Delta t$', fontsize=18, ha='center', va='top')
    fig.savefig(path, dpi=DPI, facecolor='white')
    plt.close(fig)


def base_slide():
    fig = canvas('Rotation delays the beam that travels with it',
                 'The basis of optical gyroscopes: the delay grows in proportion to the rate of rotation')
    x = 1030 / 1920

    def at(y):
        return 1 - y / 1080
    fig.text(x, at(298), r'Delay after $N$ turns, each enclosing an area $A$:', fontsize=17, color=GREY)
    fig.text(x, at(405), r'$\Delta t = \dfrac{4NA}{c^2}\,\omega$', fontsize=30)
    # why a delay is a phase: the same laser light in both beams, the co-rotating one delayed by dt
    fig.text(x, at(545), 'Why the camera sees it as a phase shift:', fontsize=17, color=GREY)
    ax = px_axes(fig, 1030, 565, 560, 200)
    tt = np.linspace(0, 3.0, 600)
    lag = .3                                            # illustrative delay, in periods
    ax.plot(tt, np.sin(2 * np.pi * tt), color=ORANGE, lw=2.4)
    ax.plot(tt, np.sin(2 * np.pi * (tt - lag)), color=BLUE, lw=2.4)
    ax.annotate('', xy=(1.25 + lag, 1.25), xytext=(1.25, 1.25),
                arrowprops=dict(arrowstyle='<->', color=INK, lw=1.6, shrinkA=0, shrinkB=0))
    ax.text(1.25 + lag + .06, 1.25, r'$\Delta t$', color=INK, fontsize=16, ha='left', va='center')
    ax.annotate('', xy=(2.25, -1.4), xytext=(1.25, -1.4),
                arrowprops=dict(arrowstyle='<->', color=GREY, lw=1.6, shrinkA=0, shrinkB=0))
    ax.text(2.31, -1.4, r'one period $T = \lambda/c$', color=GREY, fontsize=14, ha='left', va='center')
    ax.set(xlim=(0, 3.0), ylim=(-1.75, 1.6)); ax.axis('off')
    fig.text(1605 / 1920, at(640), 'against rotation', color=ORANGE, fontsize=14)
    fig.text(1605 / 1920, at(676), 'with rotation', color=BLUE, fontsize=14)
    fig.text(x, at(905), r'$\Delta\phi = \dfrac{2\pi}{T}\,\Delta t = \dfrac{8\pi NA}{\lambda c}\,\omega$',
             fontsize=28, color=BLUE)
    footer(fig, 'Animation not to scale. The rim is drawn moving at 0.25c and the loop as a single turn.')
    placeholder(fig, REGION, '02_sagnac_principle_ANIMATION.gif')
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
        assert Image.open(paths[-1]).size == REGION[2:]
        base_slide()                                    # the slide, with a guide box where the GIF goes
        write_gif(paths, durations, OUT / '02_sagnac_principle_ANIMATION.gif', colors=128)


if __name__ == '__main__':
    main()
