"""Static Sagnac slide images (1920 x 1080). Run: python slides.py [names...]

Slides 02 and 05 (the animated ones) are made by anim_concept.py and anim_measurements.py.
"""
import sys
from math import floor, log10

import numpy as np
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch
from matplotlib.ticker import MultipleLocator
from matplotlib import patheffects
from scipy import stats

from deck import *


# ---- helpers ----------------------------------------------------------------------------
def at(y):
    """Figure y coordinate of a slide pixel row (top = 0)."""
    return 1 - y / 1080


def fx(x):
    return x / 1920


def fmt_p(p, bold=False):
    """p value: two decimals, or a power of ten in math text when small."""
    if p >= 0.01:
        return f'{p:.2f}'
    e = floor(log10(p))
    m = p / 10 ** e
    if bold:
        return r'$\mathbf{%d\times10^{%d}}$' % (round(m), e)
    return r'$%d\times10^{%d}$' % (round(m), e)


def pm(v, u):
    """'value ± uncertainty' with the uncertainty rounded to one significant digit and the value to the same
    decimal place (1.738, 0.137 -> '1.7 ± 0.1')."""
    d = -int(floor(log10(abs(u))))
    if round(u, d) >= 10 ** (1 - d):                  # 0.096 rounds up to 0.1: one place fewer
        d -= 1
    d = max(d, 0)
    return f'{round(v, d):.{d}f} ± {round(u, d):.{d}f}'.replace('-', '−')


def chi_text(chi, dof, p, decimals=1):
    return (r'$\chi^2/\mathrm{dof}$ = ' + f'{chi:.{decimals}f}/{dof},  ' + r'$p$ = ' + fmt_p(p))


def head(fig, x, y, txt):
    fig.text(fx(x), at(y), txt, fontsize=14, color=GREY, weight='bold')


def draw_table(fig, x0, y0, cols, rows, row_h=46, size=17, head_size=14, width=None, bold_first=False,
               colors=None, header=True):
    """Simple striped table. cols: x offsets (px) of the columns, rows: list of tuples (the first is the
    header when header=True). y0 is the baseline of the first row."""
    width = width or (cols[-1] + 260)
    for i, r in enumerate(rows):
        y = y0 + i * row_h
        is_head = header and i == 0
        if (i % 2 == 1) if header else (i % 2 == 0):
            a = px_axes(fig, x0 - 14, y - row_h * .70, width + 28, row_h); a.axis('off')
            a.add_patch(Rectangle((0, 0), 1, 1, transform=a.transAxes, fc=LIGHT, ec='none'))
        for j, (txt, cx) in enumerate(zip(r, cols)):
            c = GREY if is_head else (colors[i][j] if colors and colors[i] and colors[i][j] else INK)
            w = 'bold' if is_head or (bold_first and j == 0) else 'normal'
            fig.text(fx(x0 + cx), at(y), txt, fontsize=head_size if is_head else size, color=c, weight=w)


def backup_tag(fig, label=''):
    fig.text(1 - 96 / 1920, .878, f'BACKUP {label}'.strip(), fontsize=13, color=GREY, weight='bold', ha='right',
             bbox=dict(boxstyle='round,pad=0.35', fc=LIGHT, ec='none'))


def show_photo(ax, row, half=470, ring_color='white', label_rings=False, cross=True):
    img, (cx, cy) = crop(photo(row.file), row.cx, row.cy, half)
    ax.imshow(img, extent=(-half, half, half, -half))
    th = np.linspace(0, 2 * np.pi, 400)
    for k, R in enumerate(rings(row)):
        ax.plot(R * np.cos(th), R * np.sin(th), color=ring_color, lw=1.6, ls=(0, (5, 3)), alpha=.9)
        if label_rings:
            a = np.radians(50)
            ax.text(R * np.cos(a) + 10, -R * np.sin(a) - 10, f'$m$ = {k}', color='white', fontsize=13,
                    ha='left', va='bottom', weight='bold',
                    bbox=dict(boxstyle='round,pad=0.15', fc='black', ec='none', alpha=.55))
    if cross:
        ax.plot(0, 0, '+', color='white', ms=14, mew=1.8)
    ax.set_xlim(-half, half); ax.set_ylim(half, -half)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def circle_photo(fig, row, x, y, d, half=420, edge='#3b3b3b'):
    img, _ = crop(photo(row.file), row.cx, row.cy, half)
    yy, xx = np.indices(img.shape[:2])
    img = img.copy(); img[np.hypot(xx - half, yy - half) > half] = 255
    ax = px_axes(fig, x, y, d, d)
    ax.imshow(img, extent=(-half, half, half, -half))
    ax.add_patch(Circle((0, 0), half - 1, fill=False, ec=edge, lw=1))
    ax.set_xlim(-half - 4, half + 4); ax.set_ylim(half + 4, -half - 4); ax.axis('off')
    return ax


# ---- 01 ---------------------------------------------------------------------------------
def s01_title():
    fig = canvas('')
    fig.text(.075, .80, 'The Sagnac effect', fontsize=48, weight='bold')
    fig.text(.075, .725, 'The Sagnac phase in a rotating fiber loop, and the speed of light', fontsize=23, va='top')
    fig.text(.075, .585, 'Daniel Haim Breger  &  Hila Hemo', fontsize=20)
    fig.text(.075, .535, 'Instructor: Ron Nordman', fontsize=17, color=GREY)
    picks = [0.0, 20.74, 44.07, 64.00, 90.0]
    d, gap, y0 = 300, 47, 640
    x0 = (1920 - (5 * d + 4 * gap)) / 2                 # centred on the slide
    for j, S in enumerate(picks):
        row = representative(S)
        circle_photo(fig, row, x0 + j * (d + gap), y0, d)
        lab = 'at rest' if S == 0 else f'{row.omega:.2f} rad/s'
        fig.text(fx(x0 + j * (d + gap) + d / 2), at(y0 + d + 40), lab, ha='center', fontsize=17,
                 color=BLUE if S == 0 else INK, weight='bold' if S == 0 else 'normal')
    return save(fig, '01_title')


# ---- 03 ---------------------------------------------------------------------------------
def draw_setup(ax):
    """Top view of the apparatus, after Fig. 1 of the lab manual."""
    ax.set_aspect('equal'); ax.axis('off')
    ax.add_patch(Circle((0, 0), 1.0, fc='#f4f6f8', ec='#c3ccd3', lw=1.5, zorder=0))
    for r in np.linspace(.925, .975, 6):                       # fiber loop, 13 turns drawn as a bundle
        ax.add_patch(Circle((0, 0), r, fill=False, ec='#3a4652', lw=.9, alpha=.85, zorder=1))
    lw = 3.0
    ax.plot([0, 0], [-.70, 0], color=BEAM, lw=lw, zorder=2, solid_capstyle='butt')
    ax.plot([0, 0], [0, .80], color=BEAM, lw=lw, zorder=2)
    ax.plot([-.80, 0], [0, 0], color=BEAM, lw=lw, zorder=2)
    ax.plot([0, .66], [0, 0], color=BEAM, lw=lw, zorder=2, ls=(0, (1.6, 1.0)))

    def lens(x, y, vertical=True, label='', dx=0, dy=0):
        w, h = (.035, .13) if vertical else (.13, .035)
        ax.add_patch(matplotlib.patches.Ellipse((x, y), w, h, fc='white', ec=INK, lw=1.4, zorder=4))
        if label:
            ax.text(x + dx, y + dy, label, fontsize=13, ha='center', va='center', color=INK)

    def polarizer(x, y, vertical=True, label='', dx=0, dy=0):
        w, h = (.03, .14) if vertical else (.14, .03)
        ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, fc='#b9c3cb', ec=INK, lw=1.2, zorder=4))
        if label:
            ax.text(x + dx, y + dy, label, fontsize=13, ha='center', va='center', color=INK)

    s = .085
    ax.add_patch(Rectangle((-s, -s), 2 * s, 2 * s, fc='white', ec=INK, lw=1.6, zorder=5))
    ax.plot([-s, s], [-s, s], color=INK, lw=1.4, zorder=6)
    ax.text(-.13, -.17, 'beam\nsplitter', fontsize=13, ha='right', va='top', color=INK, linespacing=1.1)
    ax.add_patch(FancyBboxPatch((-.05, -.86), .10, .17, boxstyle='round,pad=0.005', fc='#e2a400', ec=INK, lw=1.2,
                                zorder=4))
    ax.text(.08, -.775, 'laser', fontsize=13, ha='left', va='center', color=INK)
    polarizer(0, -.42, vertical=False, label='P₁', dx=.12)
    lens(-.62, 0, True, 'L₁', dy=-.12)
    lens(0, .62, False, 'L₂', dx=.12)
    for x, y in [(-.86, 0), (0, .86)]:
        ax.add_patch(FancyBboxPatch((x - .045, y - .045), .09, .09, boxstyle='round,pad=0.005', fc='#c9d1d8',
                                    ec=INK, lw=1.2, zorder=4))
    polarizer(.38, 0, True, 'P₂', dy=-.13)
    lens(.55, 0, True, 'L₃', dy=-.12)
    ax.add_patch(FancyBboxPatch((.66, -.065), .16, .13, boxstyle='round,pad=0.01', fc=INK, ec=INK, zorder=4))
    ax.add_patch(Rectangle((.62, -.035), .05, .07, fc=INK, ec=INK, zorder=4))
    ax.text(.74, .12, 'camera', fontsize=13, ha='center', color=INK)
    for ang0, ang1, col, r in [(115, 165, BLUE, 1.09), (165, 115, ORANGE, 1.18)]:
        th = np.radians(np.linspace(ang0, ang1, 60))
        ax.plot(r * np.cos(th), r * np.sin(th), color=col, lw=3, solid_capstyle='round')
        ax.annotate('', xy=(r * np.cos(th[-1]), r * np.sin(th[-1])),
                    xytext=(r * np.cos(th[-6]), r * np.sin(th[-6])),
                    arrowprops=dict(arrowstyle='-|>,head_length=.7,head_width=.4', color=col, lw=3))
    ax.text(-1.17, .12, 'beam 1', color=BLUE, fontsize=15, ha='center', va='center', weight='bold')
    ax.text(-.62, 1.17, 'beam 2', color=ORANGE, fontsize=15, ha='right', va='center', weight='bold')
    th = np.radians(np.linspace(20, 70, 60))
    ax.plot(1.12 * np.cos(th), 1.12 * np.sin(th), color=INK, lw=2.2)
    ax.annotate('', xy=(1.12 * np.cos(th[-1]), 1.12 * np.sin(th[-1])),
                xytext=(1.12 * np.cos(th[-5]), 1.12 * np.sin(th[-5])),
                arrowprops=dict(arrowstyle='-|>,head_length=.8,head_width=.45', color=INK, lw=2.2))
    ax.text(1.03, .96, r'$\omega$', fontsize=26, color=INK, ha='left', va='bottom')
    ax.annotate('fiber loop, 13 turns', xy=(.95 * np.cos(np.radians(222)), .95 * np.sin(np.radians(222))),
                xytext=(-1.52, -1.13), fontsize=14, color='#3a4652', ha='left', va='center',
                arrowprops=dict(arrowstyle='-', color='#3a4652', lw=1, shrinkA=2, shrinkB=0))
    ax.set_xlim(-1.62, 1.62); ax.set_ylim(-1.32, 1.32)


def s03_system():
    fig = canvas('Both beams travel the same 13-turn fiber loop',
                 'Top view of the system. Everything drawn here turns with the table.')
    ax = px_axes(fig, 110, 165, 1000, 840)
    draw_setup(ax)
    x0 = 1100
    head(fig, x0, 205, 'SYSTEM PARAMETERS  (given)')
    draw_table(fig, x0, 255, [0, 410], [
        ('Fiber turns  $N$', '13 ± 1'),
        ('Loop radius  $R$', '57.0 ± 0.5 cm'),
        (r'Wavelength  $\lambda$', '640 nm'),
        ('Fiber', 'SMF-28e+'),
    ], row_h=42, size=17, width=690, header=False)
    head(fig, x0, 460, 'WHAT WE MEASURE')
    draw_table(fig, x0, 510, [0, 410], [
        (r'Angular velocity  $\omega$', 'stopwatch'),
        (r'Phase shift  $\Delta\phi$', 'camera photos'),
    ], row_h=42, size=17, width=690, header=False)
    head(fig, x0, 670, 'PREDICTION')
    fig.text(fx(x0), at(770), r'$\Delta\phi = \dfrac{8\pi N A}{\lambda c}\,\omega \equiv a_{\rm th}\,\omega$',
             fontsize=29, color=BLUE)
    fig.text(fx(x0), at(868), r'$a_{\rm th}$ = ' + pm(RES['theory'], RES['u_theory']) + ' rad/(rad/s)', fontsize=20)
    fig.text(fx(x0), at(912), r'with $A = \pi R^2$, the area of one turn', fontsize=15, color=GREY)
    footer(fig, 'Schematic of the setup. P = polarizer, L₁ and L₂ = fiber coupling lenses, '
                'L₃ = camera lens.\nThe polarizers reduce reflections. ' r'$\lambda$ was given without an uncertainty.')
    return save(fig, '03_system')


# ---- 04 ---------------------------------------------------------------------------------
def s04_calibration():
    fig = canvas('A stopwatch calibration converts motor setting to ω',
                 r'Mean turn time $T$ at nine settings $S$ of the motor controller ($S$ sets the motor voltage), $\omega = 2\pi/T$')
    tab = calibration_table()
    lin = CAL['linear']
    ax = px_axes(fig, 175, 180, 1010, 700)
    for j, S in enumerate(SPEEDS.S[1:]):
        ax.axvline(S, color='#9aa5ae', lw=1.3, ls=(0, (4, 3)), zorder=1,
                   label='Motor settings used for the photos' if j == 0 else None)
    xs = np.linspace(tab.S.min() - 1, tab.S.max() + 2, 100)
    ax.plot(xs, C0 + C1 * xs, color=ORANGE, lw=2.6, label='Linear fit', zorder=2)
    ax.errorbar(tab.S, tab.omega, tab.u, fmt='o', color=BLUE, ms=10, capsize=0, zorder=3, label='Timed turns')
    ax.set(xlim=(0, 95), ylim=(0, 2.0), xlabel=r'Motor setting $S$', ylabel=r'Angular velocity $\omega$ (rad/s)')
    ax.yaxis.set_major_locator(MultipleLocator(.5))
    style(ax)
    h, l = ax.get_legend_handles_labels()
    ax.legend([h[2], h[1], h[0]], [l[2], l[1], l[0]], loc='upper left', frameon=True, facecolor='white',
              edgecolor='none', framealpha=1)
    # right column: the fit, the error bars, and how one calibration point is made
    S_all, laps = load_turns()
    j = int(np.argmax(S_all))
    lap = laps[j]
    rel = tab.u_sw / tab.omega * 100
    x = 1270
    head(fig, x, 225, 'LINEAR FIT')
    fig.text(fx(x), at(280), r'$\omega = kS + \omega_0$', fontsize=22)
    fig.text(fx(x), at(333), r'$k$ = ' + paren(lin['k'], lin['u_k']) + ' rad/s per unit S', fontsize=17)
    fig.text(fx(x), at(375), chi_text(lin['chi2'], lin['dof'], lin['p']), fontsize=17)
    dev = np.max(np.abs(tab.omega - (C0 + C1 * tab.S)) / tab.omega) * 100
    fig.text(fx(x), at(412), f'points within {dev:.1f} % of the line,', fontsize=15, color=GREY)
    fig.text(fx(x), at(443), f'error bars {np.max(tab.u / tab.omega) * 100:.1f} % of ω', fontsize=15, color=GREY)
    head(fig, x, 500, f'ONE POINT: S = {S_all[j]:g}, {len(lap)} TURNS TIMED')
    fig.text(fx(x), at(534), 'Each dot is one turn (stopwatch splits).', fontsize=14, color=GREY)
    pa = px_axes(fig, 1365, 570, 465, 305)
    pa.plot(np.arange(1, len(lap) + 1), lap, 'o', color=BLUE, ms=8, zorder=3)
    pa.axhline(lap.mean(), color=ORANGE, lw=2)
    pa.set(xticks=np.arange(1, len(lap) + 1), ylim=(3.40, 3.80), xlim=(.4, len(lap) + .6), yticks=[3.45, 3.55, 3.65],
           xlabel='Turn', ylabel='Turn time (s)')
    style(pa)
    pa.xaxis.label.set_size(15); pa.yaxis.label.set_size(15)
    pa.tick_params(labelsize=14)
    pa.text(.03, .95, r'$\langle T\rangle$' + f' = {lap.mean():.3f} s', transform=pa.transAxes, fontsize=15,
            color=ORANGE, va='top')
    pa.text(.03, .85, r'$\omega = 2\pi/\langle T\rangle$' + f' = {2 * np.pi / lap.mean():.3f} rad/s',
            transform=pa.transAxes, fontsize=15, va='top')
    w_line, w_timed = C0 + C1 * S_all[j], 2 * np.pi / lap.mean()
    footer(fig, f'The photos take ω from this line: {w_line:.2f} rad/s at S = {S_all[j]:g}, '
                f'{(w_line / w_timed - 1) * 100:.1f} % above the timed {w_timed:.3f} rad/s.\n'
                'This difference is part of the final uncertainty (calibration model).')
    return save(fig, '04_calibration')


# ---- 06 ---------------------------------------------------------------------------------
RING_COLOURS = ['#7fdcff', '#ffe066', '#9df08c']


def s06_rings():
    fig = canvas('The rings are evenly spaced in r², not in r',
                 'Fastest-rotation photo: dark-ring minima found along 170 rays, concentric circles fitted to them')
    row = representative(90.0)
    rm = ray_minima_for(row.file)
    half = 560
    ax = px_axes(fig, 150, 185, 520, 520)
    img, _ = crop(photo(row.file), row.cx, row.cy, half)
    ax.imshow(img, extent=(-half, half, half, -half))
    rx, ry = rm['ray_centre'][0] - row.cx, rm['ray_centre'][1] - row.cy
    # rays: faded in from the centre, where 170 of them would otherwise merge into a white disc
    from matplotlib.collections import LineCollection
    rr = np.linspace(50, rm['r_max'], 30)
    fade = .32 * np.clip((.5 * (rr[:-1] + rr[1:]) - 50) / 220, 0, 1)
    segs, cols = [], []
    for t in rm['theta']:
        xs_, ys_ = rx + rr * np.cos(t), ry - rr * np.sin(t)
        for q in range(len(rr) - 1):
            segs.append([(xs_[q], ys_[q]), (xs_[q + 1], ys_[q + 1])])
            cols.append((1, 1, 1, fade[q]))
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=.45, zorder=2, capstyle='butt'))
    pts = rm['pts']
    order = np.argsort(np.argsort([np.median(np.hypot(pts[pts[:, 0] == k, 1] - row.cx, pts[pts[:, 0] == k, 2] - row.cy))
                                   for k in range(int(pts[:, 0].max()) + 1)]))
    for k in range(int(pts[:, 0].max()) + 1):
        p = pts[pts[:, 0] == k]
        ax.plot(p[:, 1] - row.cx, p[:, 2] - row.cy, 'o', ms=2.6, color=RING_COLOURS[order[k] % 3], mew=0, zorder=4)
    th = np.linspace(0, 2 * np.pi, 400)
    for k, R in enumerate(rings(row)):
        ax.plot(R * np.cos(th), R * np.sin(th), color='white', lw=1.5, ls=(0, (5, 3)), alpha=.95, zorder=3)
        # labels straight below the centre, in the sector without rays
        ax.text(0, R + 16, f'$m$ = {k}', color='white', fontsize=13, ha='center', va='top', weight='bold', zorder=5,
                bbox=dict(boxstyle='round,pad=0.15', fc='black', ec='none', alpha=.6))
    ax.plot(0, 0, '+', color='white', ms=14, mew=1.8, zorder=5)
    # quality of the circle fit, written into the empty lower-left sector
    fit_txt = f'circle fit\nrms {rm["rms"]:.0f} px'
    ax.text(-half + 20, half - 20, fit_txt, color='white', fontsize=13, ha='left', va='bottom', zorder=6,
            linespacing=1.3, bbox=dict(boxstyle='round,pad=0.3', fc='black', ec='none', alpha=.6))
    ax.set_xlim(-half, half); ax.set_ylim(half, -half)
    ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    R = np.array(rings(row))
    r, prof = sector_profile_r(row.file, row.cx, row.cy, r_max=560)
    a1 = px_axes(fig, 990, 215, 810, 255)
    a1.plot(r, prof, color=BLUE, lw=2.2)
    for x in R:
        a1.axvline(x, color=INK, lw=1, ls=(0, (4, 3)), alpha=.7)
    for x0, x1 in zip(R[:-1], R[1:]):
        a1.annotate('', xy=(x1, 190), xytext=(x0, 190), arrowprops=dict(arrowstyle='<->', color=ORANGE, lw=1.8,
                                                                        shrinkA=0, shrinkB=0))
        a1.text(.5 * (x0 + x1), 197, f'{x1 - x0:.0f} px', color=ORANGE, ha='center', va='bottom', fontsize=15,
                weight='bold')
    a1.set(xlim=(0, 560), ylim=(0, 225), xlabel=r'Radius $r$ (px)', ylabel='Intensity (counts)')
    style(a1)
    a1.set_title(r'Against $r$: the spacing shrinks outward', loc='left', fontsize=16, color=GREY, pad=14)
    a2 = px_axes(fig, 990, 625, 810, 255)
    S2 = R ** 2 / 1e5
    a2.plot(r ** 2 / 1e5, prof, color=BLUE, lw=2.2)
    for x in S2:
        a2.axvline(x, color=INK, lw=1, ls=(0, (4, 3)), alpha=.7)
    for x0, x1 in zip(S2[:-1], S2[1:]):
        a2.annotate('', xy=(x1, 190), xytext=(x0, 190), arrowprops=dict(arrowstyle='<->', color=ORANGE, lw=1.8,
                                                                        shrinkA=0, shrinkB=0))
        a2.text(.5 * (x0 + x1), 197, f'{x1 - x0:.2f}', color=ORANGE, ha='center', va='bottom', fontsize=15,
                weight='bold')
    a2.set(xlim=(0, 3.14), ylim=(0, 225), xlabel=r'Squared radius $r^2$ ($10^5$ px$^2$)',
           ylabel='Intensity (counts)')
    style(a2)
    a2.set_title(r'Against $r^2$: one constant spacing $P$', loc='left', fontsize=16, color=GREY, pad=14)
    x = 150
    fig.text(fx(x), at(752), 'Two beams with different wavefront radii:', fontsize=16, color=GREY)
    fig.text(fx(x), at(848), r'$\phi(r) = 2\pi\,\dfrac{r^2}{P} + \phi_c$', fontsize=26)
    fig.text(fx(x), at(928), r'Rotation shifts $\phi_c$, so the rings move.', fontsize=16)
    fig.text(fx(x), at(968), r'$P$ is the ring spacing in $r^2$, the same for every ring.', fontsize=16)
    footer(fig, 'Rays every 1.5° from −75° to 180° (no fringes at lower left). Dots: darkest point of each ring on each ray. '
                'All photos: rms 3–8 px, radii typically ±1 px.\n'
                r'$P = 2\lambda\,/\,|1/R_1 - 1/R_2|$ for wavefront radii $R_1$, $R_2$. Median over all photographs: '
                r'$P$' + f' = {np.median(FRAMES.P_px2) / 1e5:.2f}' + r'$\times10^5$ px$^2$.')
    return save(fig, '06_rings_in_r2')


# ---- 07 ---------------------------------------------------------------------------------
def ring_line(row):
    R = np.sort(np.array(rings(row)))
    s = R ** 2
    m = np.r_[0, np.cumsum(np.rint(np.diff(s) / 107.5e3))]
    P, s0 = np.polyfit(m, s, 1)
    return m, s, P, s0 / P


def s07_extraction():
    fig = canvas('Squared ring radii give the phase at the center',
                 'One photograph at rest and one at the highest speed. The same fit is applied to all 161.')
    rest, fast = representative(0.0), representative(90.0)
    for row, y, col, lab in [(rest, 215, BLUE, 'at rest'), (fast, 610, ORANGE, f'ω = {fast.omega:.2f} rad/s')]:
        ax = px_axes(fig, 150, y, 320, 320)
        show_photo(ax, row, half=520, label_rings=False)
        for s in ax.spines.values():
            s.set_visible(True); s.set_color(col); s.set_linewidth(4)
        fig.text(fx(150), at(y - 14), lab, color=col, fontsize=17, weight='bold')
    ax = px_axes(fig, 650, 205, 720, 610)
    etas = []
    for row, col in [(rest, BLUE), (fast, ORANGE)]:
        m, s, P, eta = ring_line(row)
        etas.append(eta)
        xs = np.linspace(0, 2.3, 50)
        ax.plot(xs, xs + eta, color=col, lw=2.6)
        ax.plot(m, s / P, 'o', color=col, ms=13, mec='white', mew=1.5, zorder=4)
    ax.axvline(0, color='#85919b', lw=1)
    ax.annotate('', xy=(-.12, etas[1]), xytext=(-.12, etas[0]),
                arrowprops=dict(arrowstyle='<->', color=INK, lw=2, shrinkA=0, shrinkB=0))
    d_eta = etas[0] - etas[1]
    ax.text(-.19, .5 * (etas[0] + etas[1]), f'{d_eta:.2f}\nfringe', ha='right', va='center', fontsize=17,
            weight='bold', linespacing=1.1)
    ax.set(xlim=(-.8, 2.3), ylim=(-.25, 3.05), xlabel=r'Dark-ring number $m$',
           ylabel=r'$R_m^2\,/\,P$   (fringe periods)')
    ax.set_xticks([0, 1, 2]); ax.yaxis.set_major_locator(MultipleLocator(1))
    style(ax)
    x = 1450
    head(fig, x, 240, 'FIT EACH PHOTO')
    fig.text(fx(x), at(305), r'$R_m^2 = P\,(m + \eta)$', fontsize=23)
    fig.text(fx(x), at(348), r'slope $P$, intercept $P\eta$', fontsize=15, color=GREY)
    head(fig, x, 445, 'PHASE AT THE CENTER')
    fig.text(fx(x), at(510), r'$\phi_c = \pi - 2\pi\eta$', fontsize=23)
    fig.text(fx(x), at(553), r'a dark ring has phase $\pi$', fontsize=15, color=GREY)
    head(fig, x, 650, 'THIS PAIR')
    fig.text(fx(x), at(695), r'$\eta$:', fontsize=19)
    fig.text(fx(x + 45), at(695), f'{etas[0]:.2f}', fontsize=19, color=BLUE, weight='bold')
    fig.text(fx(x + 150), at(695), '→', fontsize=19)
    fig.text(fx(x + 192), at(695), f'{etas[1]:.2f}', fontsize=19, color=ORANGE, weight='bold')
    fig.text(fx(x), at(758), r'$\Delta\phi = 2\pi \times %.2f$' % d_eta, fontsize=23)
    fig.text(fx(x), at(810), r'$\approx$' + f' {2 * np.pi * d_eta:.1f} rad', fontsize=23)
    footer(fig, r'Per photo: center and radii $R_m$ from the rays and circles of the previous slide '
                r'(rms %.0f–%.0f px over all photos), then a straight line in $R_m^2$.' % (
                    FRAMES.rms_px.min(), FRAMES.rms_px.max()) + '\n'
                r'$P$ is fitted for every photo, so a change of ring spacing does not move $\phi_c$.')
    return save(fig, '07_phase_extraction')


# ---- 08 ---------------------------------------------------------------------------------
def s08_result():
    """Linear regression of the ten setting means, the slope error from the scatter about the line
    (Schwartz, Eq. 7)."""
    return slide8_option('regression', '08_phase_vs_omega')


# ---- 08, the two options for the uncertainty of the line ---------------------------------
_FIT8 = {}


def fit8(mode, col='phi'):
    """The line of slide 8 (and, with col, of the other phase methods of slide 9).
    'regression': ordinary least squares on the ten setting means (equal weights), the slope error from
                  the scatter of the points about the line (as Schwartz, Eq. 7).
    'photo':      weights from the photo scatter and the table speed only, u = sqrt(SEM^2 + (a eps w)^2),
                  no sigma_x.
    The calibration-line part of the slope error is propagated by redrawing the line's two coefficients."""
    if (mode, col) in _FIT8:
        return _FIT8[mode, col]
    sp = SPEEDS
    S, y, sem = sp.S.to_numpy(), sp[col].to_numpy(), sp['sem'].to_numpy()
    n = len(y)
    c0, c1 = LIN['coef']

    def omega_of(k0, k1):
        w = k0 + k1 * S
        w[S == 0] = 0.0                                  # at rest the table does not turn
        return w

    def line(w):
        A = np.c_[w, np.ones(n)]
        if mode == 'regression':
            q, *_ = np.linalg.lstsq(A, y, rcond=None)
            r = y - A @ q
            sd = np.sqrt(r @ r / (n - 2))
            return q, sd ** 2 * np.linalg.inv(A.T @ A), np.full(n, sd)
        u = np.sqrt(sem ** 2 + (RES['slope'] * EPS * w) ** 2 + (wander_sigma()['sigma'] ** 2 if mode == 'wander' else 0))
        W = 1 / u ** 2
        C = np.linalg.inv(A.T @ (A * W[:, None]))
        return C @ (A.T @ (W * y)), C, u

    w = omega_of(c0, c1)
    q, C, u = line(w)
    k = np.random.default_rng(1).multivariate_normal([c0, c1], np.array(LIN['cov']), 4000)
    u_cal = float(np.std([line(omega_of(k0, k1))[0][0] for k0, k1 in k]))
    a, b = q
    u_stat = float(np.sqrt(C[0, 0]))
    u_a = float(np.hypot(u_stat, u_cal))
    resid = y - (a * w + b)
    chi = float(np.sum((resid / u) ** 2))
    th, uth = RES['theory'], RES['u_theory']
    ratio = a / th
    z = (a - th) / np.hypot(u_a, uth)
    A2 = np.c_[w ** 2, w, np.ones(n)]                  # linearity check: a quadratic term
    q2, *_ = np.linalg.lstsq(A2, y, rcond=None)
    r2 = y - A2 @ q2
    C2 = (r2 @ r2 / (n - 3)) * np.linalg.inv(A2.T @ A2)
    t_q = q2[0] / np.sqrt(C2[0, 0])
    r2 = 1 - np.sum(resid ** 2) / np.sum((y - y.mean()) ** 2)          # coefficient of determination
    _FIT8[mode, col] = dict(a=a, b=b, w=w, u=u, u_stat=u_stat, u_cal=u_cal, u_a=u_a, resid=resid, chi=chi, dof=n - 2,
                p=stats.chi2.sf(chi, n - 2), s=float(u[0]) if mode == 'regression' else None, ratio=ratio,
                u_ratio=ratio * np.hypot(u_a / a, uth / th), z=z, p_z=2 * stats.norm.sf(abs(z)), t_q=t_q,
                p_q=2 * stats.t.sf(abs(t_q), n - 3), sem=sem, r2=r2)
    return _FIT8[mode, col]


_WANDER = {}


def wander_sigma(n_boot=20000):
    """Phase wander within a setting, measured from the photos themselves: for every setting, the mean phase of
    the second half of its photos minus that of the first half (Delta), with its photo-noise error u_Delta.
    If the two halves wander independently by sigma_h each, Delta has variance 2 sigma_h^2 + u_Delta^2 and the
    setting mean wanders by sigma_h / sqrt(2), so sigma_x = sqrt(<Delta^2 - u_Delta^2>) / 2.
    Uncertainty: bootstrap over the ten settings."""
    if _WANDER:
        return _WANDER
    rows = []
    for S, g in FRAMES.sort_values('time_s').groupby('S'):
        p = g.phi_rad.to_numpy()
        h = len(p) // 2
        a1, a2 = p[:h], p[h:]
        rows.append((S, a2.mean() - a1.mean(), np.sqrt(a1.var(ddof=1) / len(a1) + a2.var(ddof=1) / len(a2))))
    S_, d, ud = map(np.array, zip(*rows))
    est = lambda k: np.sqrt(max(np.mean(d[k] ** 2 - ud[k] ** 2), 0)) / 2
    boot = [est(k) for k in np.random.default_rng(0).integers(0, len(d), (n_boot, len(d)))]
    omega = SPEEDS.set_index('S').loc[S_, 'omega'].to_numpy()
    _WANDER.update(S=S_, omega=omega, d=d, ud=ud, sigma=est(np.arange(len(d))), u_sigma=float(np.std(boot)))
    return _WANDER


def cal_model_slopes():
    """The Sagnac slope (linear regression of slide 8) with the linear and with the quadratic calibration."""
    qc = CAL['fits_with_eps']['2']['coef']
    S_ = SPEEDS.S.to_numpy()
    w_lin = np.where(S_ == 0, 0, C0 + C1 * S_)
    w_quad = np.where(S_ == 0, 0, qc[0] + qc[1] * S_ + qc[2] * S_ ** 2)
    return tuple(np.polyfit(w_, SPEEDS.phi.to_numpy(), 1)[0] for w_ in (w_lin, w_quad))


def final_result():
    """Slope a and c with all included uncertainties: regression (with the calibration line), method spread
    (half range of the three full-image methods) and calibration model (linear vs quadratic)."""
    fits = [fit8('regression', col) for col in ('phi', 'phi_cos', 'phi_template')]
    a, u_reg = fits[0]['a'], fits[0]['u_a']
    sl = [f['a'] for f in fits]
    u_meth = (max(sl) - min(sl)) / 2
    a_lin, a_quad = cal_model_slopes()
    u_model = abs(a_quad - a_lin)
    u = float(np.sqrt(u_reg ** 2 + u_meth ** 2 + u_model ** 2))
    cv = 8 * np.pi * N_TURNS * AREA / (LAM * a)
    rel = {r'Turn count $N$ = 13 ± 1': 1 / N_TURNS, 'Phase-reading method': u_meth / a,
           r'Loop radius $R$ = 57.0 ± 0.5 cm': 2 * .005 / R_LOOP, 'Calibration model': u_model / a,
           'Regression and calibration line': u_reg / a}
    rel_tot = float(np.sqrt(sum(v ** 2 for v in rel.values())))
    uc = cv * rel_tot
    z = (cv - C_DEF) / uc
    return dict(a=a, u_reg=u_reg, u_meth=u_meth, u_model=u_model, u=u, slopes=sl, fits=fits, c=cv, u_c=uc,
                rel=rel, rel_tot=rel_tot, z=z, p_z=2 * stats.norm.sf(abs(z)), a_lin=a_lin, a_quad=a_quad)


def slide8_option(mode, out_name=None):
    """Slide 8 for one treatment of the scatter between settings. The deck uses 'regression' (s08_result)."""
    f = fit8(mode)
    reg = mode == 'regression'
    fig = canvas('The phase grows linearly with ω, as Sagnac predicts',
                 'Mean phase per motor setting against the calibrated angular velocity')
    sp = SPEEDS
    a, b = f['a'], f['b'] - sp.phi.iloc[0]
    u_om = np.hypot(sp.u_omega_curve, sp.u_omega_repro)
    xs = np.linspace(0, 1.95, 100)
    ax = px_axes(fig, 175, 175, 1020, 500)
    th, uth = RES['theory'], RES['u_theory']
    ax.fill_between(xs, (th - uth) * xs, (th + uth) * xs, color=THEORY, alpha=.16, lw=0,
                    label=r'Prediction $\pm1\sigma$')
    ax.plot(xs, th * xs, color=THEORY, lw=2.2, ls=(0, (6, 3)))
    ax.plot(xs, a * xs + b, color=ORANGE, lw=2.8, label='Linear regression' if reg else 'Weighted linear fit')
    ax.errorbar(sp.omega, sp.dphi, None if reg else f['u'], u_om, fmt='o', color=BLUE, ms=10, mec='white', mew=1.3,
                elinewidth=2, capsize=0, zorder=4, label='Measured')
    ax.set(xlim=(-.05, 1.95), ylim=(-.25, 3.6), ylabel=r'Phase shift $\Delta\phi$ (rad)')
    ax.tick_params(labelbottom=False)
    style(ax)
    h, l = ax.get_legend_handles_labels()
    ax.legend([h[2], h[1], h[0]], [l[2], l[1], l[0]], loc='upper left', frameon=True, facecolor='white',
              edgecolor='none', framealpha=1)
    pr = ax.twinx()
    pr.set_ylim(np.array(ax.get_ylim()) / (2 * np.pi)); pr.set_ylabel('fringes', color=GREY)
    pr.tick_params(colors=GREY); pr.spines['right'].set_visible(True); pr.spines['right'].set_color('#85919b')
    pr.yaxis.set_major_locator(MultipleLocator(.1))
    ar = px_axes(fig, 175, 700, 1020, 160)
    ar.axhline(0, color=ORANGE, lw=2)
    if reg:
        ar.plot(sp.omega, f['resid'], 'o', color=BLUE, ms=8, zorder=3)
    else:
        ar.errorbar(sp.omega, f['resid'], f['u'], fmt='o', color=BLUE, ms=8, elinewidth=2, capsize=4, zorder=3)
    ar.set(xlim=(-.05, 1.95), ylim=(-.15, .15), xlabel=r'Angular velocity $\omega$ (rad/s)', ylabel=r'$\Delta\phi$ − line (rad)')
    ar.yaxis.set_major_locator(MultipleLocator(.1))
    style(ar)
    x = 1345
    head(fig, x, 210, 'MEASURED SLOPE')
    fig.text(fx(x), at(262), 'a = ' + pm(a, f['u_a']), fontsize=25, color=ORANGE, weight='bold')
    fig.text(fx(x), at(294), 'rad/(rad/s)', fontsize=14, color=GREY)
    fig.text(fx(x), at(324), 'regression and calibration line', fontsize=14, color=GREY)
    head(fig, x, 375, r'PREDICTED  $8\pi^2NR^2/(\lambda c)$')
    fig.text(fx(x), at(427), r'$\mathbf{a_{th}}$' + ' = ' + pm(th, uth), fontsize=25, color=THEORY,
             weight='bold')
    if reg:
        head(fig, x, 505, 'AGREEMENT')
        gx = px_axes(fig, x + 95, 522, 395, 88)                     # the two slopes on one axis
        gx.errorbar([a], [1], xerr=[f['u_a']], fmt='o', color=ORANGE, ms=10, elinewidth=3, capsize=6, capthick=2)
        gx.errorbar([th], [0], xerr=[uth], fmt='o', color=THEORY, ms=10, elinewidth=3, capsize=6, capthick=2)
        gx.set(xlim=(1.55, 1.95), ylim=(-.6, 1.6), xticks=[1.6, 1.7, 1.8, 1.9], yticks=[0, 1])
        gx.set_yticklabels([r'$a_{\rm th}$', r'$a$'])
        gx.tick_params(labelsize=14)
        style(gx)
        gx.grid(axis='x', alpha=.15)
        fig.text(fx(x), at(698), f'{f["z"]:.1f}σ apart,  ' + r'$p$' + f' = {f["p_z"]:.2f}', fontsize=17)
        fig.text(fx(x), at(730), 'with the regression error only', fontsize=14, color=GREY)
        head(fig, x, 790, 'LINEAR REGRESSION')
        fig.text(fx(x), at(825), 'ten setting means, equal weights', fontsize=14, color=GREY)
        fig.text(fx(x), at(865), r'$R^2$' + f' = {f["r2"]:.4f}', fontsize=17)
        fig.text(fx(x), at(907), r'$\omega^2$ term: $p$' + f' = {f["p_q"]:.2f}, no bend', fontsize=17)
        name = '08_option_A_linear_regression'
    else:
        head(fig, x, 510, 'AGREEMENT')
        fig.text(fx(x), at(556), r'$a\,/\,a_{\rm th}$' + ' = ' + pm(f['ratio'], f['u_ratio']), fontsize=20)
        fig.text(fx(x), at(592), f'{f["z"]:.1f}σ apart,  ' + r'$p$' + f' = {f["p_z"]:.2f}', fontsize=17)
    if reg:
        pass
    elif mode == 'wander':
        wd = wander_sigma()
        head(fig, x, 660, 'WANDER WITHIN A SETTING')
        wx = px_axes(fig, 1440, 685, 390, 110)
        wx.axhline(0, color=INK, lw=1)
        wx.errorbar(wd['omega'], wd['d'], wd['ud'], fmt='o', color=PURPLE, ms=6, elinewidth=1.8, capsize=0, zorder=3)
        wx.set(xlim=(-.1, 1.95), ylim=(-.2, .2), xlabel=r'$\omega$ (rad/s)')
        wx.set_ylabel(r'$\Delta$ (rad)', fontsize=14)
        wx.xaxis.label.set_size(14); wx.tick_params(labelsize=12)
        wx.yaxis.set_major_locator(MultipleLocator(.1))
        style(wx)
        fig.text(fx(x), at(912), r'$\sigma_x = \frac{1}{2}\sqrt{\langle\Delta^2 - u_\Delta^2\rangle}$ = '
                 + paren(wd['sigma'], wd['u_sigma']) + ' rad', fontsize=15)
        fig.text(fx(x), at(955), 'fit: ' + chi_text(f['chi'], f['dof'], f['p']), fontsize=15)
        footer(fig, r'Error bars: $u = \sqrt{\mathrm{SEM}^2 + (a\,\varepsilon\,\omega)^2 + \sigma_x^2}$, all smaller than the '
                    r'markers in the upper panel. $\Delta$: mean phase of the second half of a setting’s photos' + '\n'
                    r'minus the first half, $u_\Delta$ its photo-noise error. $\sigma_x$ is measured from this '
                    'wander, not fitted, and the line then passes the χ² test.')
        name = '08_option_C_measured_wander'
    else:
        head(fig, x, 675, 'FIT STATISTICS')
        fig.text(fx(x), at(710), 'photo scatter and table speed', fontsize=14, color=GREY)
        fig.text(fx(x), at(748), chi_text(f['chi'], f['dof'], f['p']), fontsize=16)
        fig.text(fx(x), at(790), f'The points scatter {np.sqrt(f["chi"] / f["dof"]):.1f} times more', fontsize=14,
                 color=GREY)
        fig.text(fx(x), at(820), 'than their error bars allow.', fontsize=14, color=GREY)
        footer(fig, r'Error bars: $u = \sqrt{\mathrm{SEM}^2 + (a\,\varepsilon\,\omega)^2}$. '
                    r'In the upper panel all error bars are smaller than the markers, the residual panel shows them.')
        name = '08_option_B_without_extra_scatter'
    print(f'{mode}: a = {a:.4f} +- {f["u_a"]:.4f} (stat {f["u_stat"]:.4f}, cal {f["u_cal"]:.4f}), '
          f'chi2 {f["chi"]:.2f}/{f["dof"]}, ratio {f["ratio"]:.3f} +- {f["u_ratio"]:.3f}, z {f["z"]:.2f}, '
          f'p {f["p_z"]:.3f}, quad t {f["t_q"]:.2f} p {f["p_q"]:.2f}')
    return save(fig, out_name or name)


# ---- 09 ---------------------------------------------------------------------------------
METHOD_ROWS = [  # label, results key, column in sagnac_speeds.csv, colour, marker, face
    ('Ring minima\n(main)', None, None, BLUE, 'o', BLUE),
    ('Cosine fit\nof the profile', 'cosine fit', 'phi_cos', ORANGE, 'o', ORANGE),
    ('Template match\nto rest profile', 'template matching', 'phi_template', PURPLE, 's', PURPLE),
    ('10 random pixels\n(mean of 1000)', 'random points', 'phi_points', GREY, 'D', 'white'),
]
SHORT = ['Ring minima (main)', 'Cosine fit', 'Template match', '10 random pixels'
]


def s09_methods():
    fig = canvas('Three full-image methods differ by at most 6 %',
                 'The same 161 photographs analyzed with four ways of reading the phase')
    sp = SPEEDS
    ax = px_axes(fig, 175, 215, 760, 540)
    ax.axhline(0, color=BLUE, lw=2.4, label='Ring minima (main)')
    for (lab, key, col, c, mk, mfc), short in zip(METHOD_ROWS[1:3], SHORT[1:3]):      # the full-image methods
        d = (sp[col] - sp[col].iloc[0]) - sp.dphi
        ax.plot(sp.omega, d, color=c, lw=1.8, marker=mk, ms=9, mfc=mfc, mec=c, mew=1.8, label=short)
    ax.set(xlim=(-.05, 1.9), ylim=(-.22, .28), xlabel=r'Angular velocity $\omega$ (rad/s)',
           ylabel='Difference from minima (rad)')
    ax.yaxis.set_major_locator(MultipleLocator(.1))
    style(ax)
    ax.legend(loc='upper left', ncol=2, fontsize=14)
    ax.set_title('Phase change relative to the ring-minima result', loc='left', fontsize=16, color=GREY, pad=12)
    bx = px_axes(fig, 1255, 215, 540, 540)
    fr_ = final_result()
    st = scan_stats()
    slopes = [(f['a'], f['u_a']) for f in fr_['fits']] + [(st['fixed']['mean'], st['fixed']['u_fit'])]
    th, uth = RES['theory'], RES['u_theory']
    bx.axvspan(th - uth, th + uth, color=THEORY, alpha=.15, lw=0)
    bx.axvline(th, color=THEORY, lw=2, ls=(0, (6, 3)))
    for j, ((lab, key, col, c, mk, mfc), (v, u)) in enumerate(zip(METHOD_ROWS, slopes)):
        y = len(METHOD_ROWS) - 1 - j
        bx.errorbar(v, y, xerr=u, fmt=mk, color=c, ms=11, mfc=mfc, mec=c, mew=2, elinewidth=2.4, capsize=0)
        bx.text(v, y + .27, f'{v:.2f}', ha='center', fontsize=15, color=c, weight='bold', zorder=5,
                path_effects=[patheffects.withStroke(linewidth=4, foreground='white')])
    lo, hi = min(v for v, _ in slopes[:3]), max(v for v, _ in slopes[:3])
    bx.axhline(.5, color='#c9d1d8', lw=1.2, ls=(0, (4, 3)))
    bx.annotate('', xy=(lo, -.5), xytext=(hi, -.5), arrowprops=dict(arrowstyle='|-|', color=INK, lw=1.6,
                                                                     shrinkA=0, shrinkB=0, mutation_scale=6))
    bx.text(.5 * (lo + hi) + .01, -.72, 'range', fontsize=13, ha='center', va='top', color=INK)
    bx.set_yticks(range(len(METHOD_ROWS)))
    bx.set_yticklabels([r[0] for r in METHOD_ROWS[::-1]], fontsize=14)
    bx.set(xlim=(1.45, 2.0), ylim=(-1.05, 3.85), xlabel=r'Fitted slope $a$ (rad/(rad/s))')
    bx.grid(axis='x', alpha=.15)
    bx.set_title(r'Fitted slope  (gray: prediction $\pm1\sigma$)', loc='left', fontsize=15, color=GREY, pad=12)
    fig.text(fx(175), at(893), 'a = ' + pm(fr_['a'], fr_['u']) + ' rad/(rad/s)', fontsize=21, weight='bold')
    pc = {k: fr_[k] / fr_['a'] * 100 for k in ('u_reg', 'u_meth', 'u_model', 'u')}
    fig.text(fx(175), at(935), f'{pc["u_reg"]:.1f} % regression and calibration line, {pc["u_meth"]:.1f} % method, '
             f'{pc["u_model"]:.1f} % calibration model, in quadrature {pc["u"]:.1f} %', fontsize=15, color=GREY)
    sc = [f['s'] for f in fr_['fits']]
    footer(fig, 'Each slope comes from the same linear regression as the previous slide, ' + r'$R^2$ = '
                + ', '.join(f'{f["r2"]:.4f}' for f in fr_['fits']) + ' from top. Error bars: regression and calibration line.')
    return save(fig, '09_methods')


# ---- 10 ---------------------------------------------------------------------------------
ORIGINAL_SEED = 20261011                     # the draw of sagnac_update/methods.py (slide 9)


def scan_stats(n_boot=2000):
    """Slopes of the random-pixel scan: 1000 draws of ten pixels, and the extra draws with more pixels.
    Draws with a negative slope are failed fits and are left out. Centre: mean, spread: standard deviation.
    Their uncertainties come from resampling the draws (bootstrap)."""
    d = SCAN['draws']
    rng = np.random.default_rng(0)

    def summary(x_all):
        x = x_all[x_all > 0]
        k = rng.integers(0, len(x), (n_boot, len(x)))
        return dict(s=x_all, n_failed=int(np.sum(x_all <= 0)), mean=x.mean(), u_mean=x.std(ddof=1) / np.sqrt(len(x)),
                    sd=x.std(ddof=1), u_sd=np.std(x[k].std(axis=1, ddof=1)))

    orig = [r for r in d if r['seed'] == ORIGINAL_SEED][0]
    ten = [r for r in d if r['n'] == 10 and r['seed'] != ORIGINAL_SEED]
    out = dict(orig=orig['fixed']['slope'], n_draws=len(ten))
    for key in ('fixed', 'frame'):
        out[key] = summary(np.array([r[key]['slope'] for r in ten]))
        out[key]['u_fit'] = np.median([r[key]['u_fit'] for r in ten])
        out[key]['not_converged'] = int(sum(r[key]['iterations'] >= SCAN.get('max_iter', 500) for r in ten))
    by_n = {}
    for n in sorted({r['n'] for r in d}):
        rows = [r for r in d if r['n'] == n and r['seed'] != ORIGINAL_SEED]
        sm = summary(np.array([r['fixed']['slope'] for r in rows]))
        by_n[n] = dict(count=len(rows), mean=sm['mean'], sd=sm['sd'], u_sd=sm['u_sd'],
                       u_fit=np.median([r['fixed']['u_fit'] for r in rows]))
    out['by_n'] = by_n
    return out


def s10_random_pixels():
    st = scan_stats()
    nd = st['n_draws']
    fig = canvas('The slope from ten random pixels depends on the draw',
                 f'The pixel method repeated for {nd} random draws of ten pixels')
    fx_, fr_ = st['fixed'], st['frame']
    ax = px_axes(fig, 175, 205, 905, 610)
    lo, hi = 0.8, 2.72
    bins = np.arange(lo, hi + 1e-9, .04)
    n_out = int(np.sum((fx_['s'] < lo) | (fx_['s'] > hi)))
    cnt, _, _ = ax.hist(fx_['s'], bins=bins, color='#8a96a0', edgecolor='white', lw=.6,
                        label='fixed pixel positions (the method)')
    cnt2, _, _ = ax.hist(fr_['s'], bins=bins, histtype='step', color=PURPLE, lw=2.4,
                         label='same pixels, each photo’s own center and $P$')
    top = max(cnt.max(), cnt2.max())
    ymax = top * 1.75
    hl = top * 1.08 / ymax
    ax.axvline(RES['slope'], ymax=hl, color=BLUE, lw=3, label=f'ring minima: {RES["slope"]:.2f}')
    ax.axvline(RES['theory'], ymax=hl, color=THEORY, lw=2.2, ls=(0, (6, 3)), label='prediction: ' + pm(RES['theory'], RES['u_theory']))
    ax.set(xlabel=r'Fitted slope $a$ (rad/(rad/s))', ylabel='Number of draws', ylim=(0, ymax), xlim=(lo, hi))
    style(ax)
    ax.legend(loc='upper left', fontsize=14)
    x = 1200
    head(fig, x, 232, f'SPREAD OVER {nd} DRAWS')
    draw_table(fig, x, 282, [0, 420], [
        ('Mean slope', paren(fx_['mean'], fx_['u_mean'])),
        ('Spread (SD)', paren(fx_['sd'], fx_['u_sd'])),
        ('Typical fit error', f'{fx_["u_fit"]:.3f}'),
        ('Spread ÷ fit error', f'{fx_["sd"] / fx_["u_fit"]:.1f}'),
    ], row_h=44, size=17, width=625, header=False)
    by_n = st['by_n']
    if len(by_n) > 1:
        head(fig, x, 515, 'MORE PIXELS PER DRAW')
        bx = px_axes(fig, x + 95, 545, 520, 170)
        ns = np.array(sorted(by_n))
        sds = np.array([by_n[n]['sd'] for n in ns]); usd = np.array([by_n[n]['u_sd'] for n in ns])
        ufit = np.array([by_n[n]['u_fit'] for n in ns])
        bx.errorbar(ns, sds, usd, fmt='o-', color=INK, ms=8, capsize=0, elinewidth=2, lw=1.4, zorder=3,
                    label='spread between draws')
        bx.plot(ns, ufit, 'o--', color=GREY, mfc='white', mew=2, ms=8, lw=1.4, zorder=2, label='fit error of one draw')
        bx.set_xscale('log')
        bx.set_xticks(ns); bx.set_xticklabels([str(n) for n in ns])
        bx.minorticks_off()
        bx.set(xlabel='Pixels per draw', ylabel='rad/(rad/s)', ylim=(0, .25))
        bx.xaxis.label.set_size(14); bx.yaxis.label.set_size(13); bx.tick_params(labelsize=13)
        bx.yaxis.set_major_locator(MultipleLocator(.1))
        style(bx)
        hh, ll = bx.get_legend_handles_labels()
        bx.legend(hh[::-1], ll[::-1], loc='upper right', fontsize=12, handlelength=1.6)
    head(fig, x, 820, 'LIMITATIONS IN OUR RUN')
    for k, t in enumerate(['The rings move by only half a fringe, so', 'each pixel sees less than one bright-dark',
                           'cycle. A fixed pixel also reads a moving', r'center or a change of $P$ as a phase shift.']):
        fig.text(fx(x), at(856 + 32 * k), t, fontsize=15)
    footer(fig, 'Each draw: 10 pixels in the fringe sector, read at the same image position in all 161 photos. One '
                'least-squares fit gives every photo’s phase.\n'
                + (f'Spread: standard deviation of the slopes. {fx_["n_failed"]} draws with a negative slope (failed fits) '
                   'are left out and lie off the plot.' if fx_['n_failed'] else 'Spread: standard deviation of the slopes.'))
    return save(fig, '10_random_pixels')


# ---- 11 ---------------------------------------------------------------------------------
def s11_speed_of_light():
    from matplotlib.transforms import blended_transform_factory
    fig = canvas('The slope gives the speed of light, 0.7σ from c₀',
                 r'The Sagnac slope solved for $c$, with the measured $a$ and the given $N$, $R$ and $\lambda$')
    fr_ = final_result()
    cv, uc = fr_['c'] / 1e8, fr_['u_c'] / 1e8
    xl, xr, top = 175, 1000, 520                     # left column, right column, panel headings
    # result line
    f1 = fig.text(fx(xl), at(315), r'$c = \dfrac{8\pi^2 N R^2}{\lambda\,a}\; =$', fontsize=30, va='center')
    fig.canvas.draw()
    x2 = f1.get_window_extent().x1 / fig.dpi / fig.get_size_inches()[0] + 16 / 1920
    v, u = pm(cv, uc).split(' ± ')
    fig.text(x2, at(315), f'({v} ± {u}) × 10⁸ m/s', fontsize=40, weight='bold', color=BLUE, va='center')
    # left: comparison with the defined value
    head(fig, xl, top, 'COMPARED WITH THE DEFINED VALUE')
    ax = px_axes(fig, xl, 565, 680, 120)
    ax.axvline(C_DEF / 1e8, color=INK, lw=1.8, ls=(0, (5, 3)))
    ax.text(C_DEF / 1e8 + .012, 1.6, r'$c_0$', fontsize=18, color=INK, va='center', ha='left')
    ax.errorbar(cv, 0.7, xerr=uc, fmt='o', color=BLUE, ms=14, elinewidth=3, capsize=8, capthick=2.5, zorder=3)
    ax.set(xlim=(2.45, 3.25), ylim=(0, 2.0), xlabel=r'Speed of light $c$ ($10^8$ m/s)')
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.grid(axis='x', alpha=.15)
    fig.text(fx(xl), at(827), f'{abs(fr_["z"]):.1f}σ from ' + r'$c_0$,  $p$' + f' = {fr_["p_z"]:.2f}', fontsize=21,
             va='center')
    # right: uncertainty budget as a list, labels and bars in fixed columns
    head(fig, xr, top, 'UNCERTAINTY BUDGET (ADDED IN QUADRATURE)')
    bx = px_axes(fig, xr, 565, 1840 - xr, 280)
    bx.set_xlim(xr, 1840); bx.set_ylim(845, 565); bx.axis('off')         # axes in slide pixels
    x_bar, full = 1505, 245                                              # bar start, length of 10 %
    items = sorted(fr_['rel'].items(), key=lambda kv: -kv[1])
    rows = [(n, v_ * 100, y) for (n, v_), y in zip(items, (590, 632, 674, 716, 758))]
    rows.append(('Total of the listed sources', fr_['rel_tot'] * 100, 827))
    for n, pct, y in rows:
        tot = n.startswith('Total')
        bx.text(xr, y, n, va='center', fontsize=16, weight='bold' if tot else 'normal')
        bx.add_patch(Rectangle((x_bar, y - 15), pct / 10 * full, 30, color=BLUE if tot else '#8a96a0', lw=0))
        bx.text(x_bar + pct / 10 * full + 10, y, f'{pct:.1f} %', va='center', fontsize=16,
                weight='bold' if tot else 'normal')
    bx.plot([xr, 1840], [792, 792], color='#c9d1d8', lw=1.2)
    fig.text(fx(xl), at(915), r'Not included: a steady drift of the phase in time, and $\lambda$ (given without an '
             'uncertainty).', fontsize=17)
    footer(fig, r'$c_0$ = 299 792 458 m/s (defined).')
    return save(fig, '11_speed_of_light')


# ---- 12, 13 -------------------------------------------------------------------------------
def s12_conclusions():
    fig = canvas('The Sagnac phase is linear in ω and matches theory', '')
    fr_ = final_result()
    v, u = pm(fr_['c'] / 1e8, fr_['u_c'] / 1e8).split(' ± ')
    items = [(pm(fr_['a'], fr_['u']), 'rad/(rad/s)',
              r'Slope of $\Delta\phi$ against $\omega$, full uncertainty.',
              r'Prediction $8\pi^2NR^2/(\lambda c)$ = ' + pm(RES['theory'], RES['u_theory']) + ', 0.7σ apart.'),
             ('3.29 ± 0.06 rad', '0.52 fringe',
              r'Largest shift, at $\omega$ = 1.81 rad/s.',
              'That is a travel-time difference of only 1.1 fs.'),
             (f'({v} ± {u}) × 10⁸ m/s', 'speed of light',
              r'Inferred $c$: ' + f'{abs(fr_["z"]):.1f}σ from the defined value.',
              r'Limited by the turn count $N$ = 13 ± 1.')]
    y = 195
    for big, unit, l1, l2 in items:
        fig.text(fx(144), at(y + 45), big, fontsize=29, weight='bold', color=BLUE)
        fig.text(fx(144), at(y + 88), unit, fontsize=16, color=GREY)
        fig.text(fx(900), at(y + 28), l1, fontsize=19)
        fig.text(fx(900), at(y + 72), l2, fontsize=17, color=GREY)
        y += 165
    ax = px_axes(fig, 144, 712, 1680, 2)
    ax.axhline(0, color='#c9d1d8', lw=1.5); ax.axis('off')
    head(fig, 144, 765, 'LIMITATION')
    fig.text(fx(144), at(805), 'A slow drift of the phase in time would change the slope, since the speeds only went up.',
             fontsize=17)
    head(fig, 144, 862, 'NEXT')
    for k, t in enumerate(['•  Take rest photos between the speeds and reverse the rotation, to test for a slow drift',
                           '•  Count the fiber turns (now 13 ± 1)',
                           '•  Couple the light into a single fiber mode (possible with this fiber, but hard to align)']):
        fig.text(fx(144), at(900 + 37 * k), t, fontsize=17)
    return save(fig, '12_conclusions')


def s13_thank_you():
    fig = canvas('')
    fig.text(.075, .62, 'Thank you', fontsize=60, weight='bold')
    fig.text(.075, .53, 'Questions?', fontsize=28, color=GREY)
    fig.text(.075, .26, 'Daniel Haim Breger  &  Hila Hemo', fontsize=20)
    fig.text(.075, .21, 'Instructor: Ron Nordman', fontsize=17, color=GREY)
    row = representative(44.07)
    circle_photo(fig, row, 1080, 160, 760, half=460)
    fig.text(fx(1460), at(965), f'ω = {row.omega:.2f} rad/s', fontsize=16, color=GREY, ha='center')
    return save(fig, '13_thank_you')


# ---- backups ----------------------------------------------------------------------------
def b1_quadratic():
    fig = canvas('A quadratic calibration would raise the slope by 1.7 %',
                 'Calibration residuals for the two models, and what each does to the result')
    backup_tag(fig, 'B1')
    tab = calibration_table()
    f1, f2 = CAL['fits_with_eps']['1'], CAL['fits_with_eps']['2']
    qc = f2['coef']
    lin = C0 + C1 * tab.S
    quad = qc[0] + qc[1] * tab.S + qc[2] * tab.S ** 2
    ax = px_axes(fig, 175, 205, 800, 560)
    ax.axhline(0, color=INK, lw=1)
    off = .7
    ax.errorbar(tab.S - off, (tab.omega - lin) / tab.omega * 100, tab.u / tab.omega * 100, fmt='o', color=ORANGE,
                ms=9, elinewidth=2, capsize=0, label='Linear (used)')
    ax.errorbar(tab.S + off, (tab.omega - quad) / tab.omega * 100, tab.u / tab.omega * 100, fmt='s', color=PURPLE,
                ms=9, elinewidth=2, capsize=0, label='Quadratic')
    ax.set(xlim=(0, 95), ylim=(-3.3, 2.3), xlabel=r'Motor setting $S$', ylabel=r'Measured $\omega$ − model (% of $\omega$)')
    style(ax)
    ax.legend(loc='lower left')
    p2 = stats.chi2.sf(f2['chi2'], f2['dof'])
    # the Sagnac slope with each calibration, by the linear regression of slide 8
    a_lin, a_quad = cal_model_slopes()
    c_of = lambda a_: 8 * np.pi * N_TURNS * AREA / (LAM * a_) / 1e8
    rows = [('', 'Linear (used)', 'Quadratic'),
            (r'$\chi^2/\mathrm{dof}$', f'{f1["chi2"]:.1f}/{f1["dof"]}', f'{f2["chi2"]:.1f}/{f2["dof"]}'),
            (r'$p$', fmt_p(f1['p'], bold=True), f'{p2:.2f}'),
            (r'Slope $a$ (rad/(rad/s))', f'{a_lin:.3f}', f'{a_quad:.3f}'),
            (r'$a\,/\,a_{\rm th}$', f'{a_lin / RES["theory"]:.3f}', f'{a_quad / RES["theory"]:.3f}'),
            (r'$c$ ($10^8$ m/s)', f'{c_of(a_lin):.3f}', f'{c_of(a_quad):.3f}'),
            (r'$\omega$ at $S$ = 90 (rad/s)', f'{C0 + C1 * 90:.3f}', f'{qc[0] + qc[1] * 90 + qc[2] * 8100:.3f}')]
    x0, y0 = 1075, 250
    for i, r in enumerate(rows):
        y = y0 + i * 58
        for j, (txt, xx) in enumerate(zip(r, [x0, 1540, 1755])):
            fig.text(fx(xx), at(y), txt, fontsize=17 if i else 15,
                     color=(GREY if i == 0 else INK) if j == 0 else (ORANGE if j == 1 else PURPLE),
                     weight='bold' if i == 0 or j > 0 else 'normal', ha='left' if j == 0 else 'center')
        if i == 0:
            a = px_axes(fig, x0, y + 18, 750, 2); a.axhline(0, color='#c9d1d8'); a.axis('off')
    for k, t in enumerate(['The main analysis uses the linear model,', 'but the table turns slower than linear',
                           r'at high settings. The choice moves $a$ and $c$', 'by 1.7 %. This difference is included in',
                           'the result as the calibration-model', 'uncertainty.']):
        fig.text(fx(x0), at(700 + 32 * k), t, fontsize=16)
    f3 = CAL['fits_with_eps']['3']
    footer(fig, r'Both fits include the speed reproducibility $\varepsilon$ = 0.48 %, set by $\chi^2/\mathrm{dof}$ = 1 for '
                r'the quadratic, so its $p$ is not an independent test.' '\n'
                'A cubic adds nothing (' + chi_text(f3['chi2'], f3['dof'], stats.chi2.sf(f3['chi2'], f3['dof'])) +
                '). The 1.7 % shift is part of the quoted uncertainty, as the calibration-model term.')
    return save(fig, 'B1_quadratic_calibration')


def b2_error_budget():
    fig = canvas('Uncertainty sources considered',
                 r'$N$, $R$ and $\lambda$ enter the prediction and $c$. All other sources enter the measured slope $a$.')
    backup_tag(fig, 'B2')
    rows = [('Source', 'Kind', 'Size', 'Included'),
            ('Photo-to-photo variation', 'statistical', 'SEM 0.004–0.028 rad *', 'yes, regression'),
            ('Table speed at a setting', 'statistical (calibration)', r'$\varepsilon$ = 0.48 % of $\omega$', 'yes, regression'),
            ('Wander between settings', 'statistical', 'about 0.04 rad', 'yes, regression'),
            ('Calibration line', 'calibration, correlated', r'0.2 % of $a$', 'yes'),
            ('Phase-reading method', 'systematic', r'3.0 % of $a$', 'yes'),
            (r'Turn count $N$ = 13 ± 1', 'given', '7.7 %', r'yes (theory, $c$)'),
            (r'Loop radius $R$ = 57.0 ± 0.5 cm', 'given', '1.8 %', r'yes (theory, $c$)'),
            (r'Wavelength $\lambda$ = 640 nm', 'not given', '0.16 % per nm', 'no'),
            ('Steady drift in time', 'systematic', r'not separable from $a$', 'no'),
            ('Linear vs quadratic calibration', 'model choice', r'1.7 % of $a$', 'yes')]
    colors = [None] + [[None, None, None, '#2e7d4f' if r[3].startswith('yes') else '#b5452b'] for r in rows[1:]]
    xs = [150, 720, 1120, 1515]
    y0 = 215
    for i, r in enumerate(rows):
        y = y0 + i * 66
        if i % 2 == 1:
            a = px_axes(fig, 130, y - 44, 1670, 62); a.axis('off')
            a.add_patch(Rectangle((0, 0), 1, 1, transform=a.transAxes, fc=LIGHT, ec='none'))
        for j, (txt, x) in enumerate(zip(r, xs)):
            col = GREY if i == 0 else (colors[i][j] or INK)
            fig.text(fx(x), at(y), txt, fontsize=14 if i == 0 else 16, color=col,
                     weight='bold' if i == 0 or j == 3 else 'normal')
    footer(fig, '* SEM: standard error of the mean of the photos at one setting. Statistical parts enter through the '
                r'regression. Others in quadrature: $a$ = 1.84 ± 0.07.' '\n'
                r'Not checked: the sign of the shift (rotation direction not recorded). $\omega$ comes from the calibration, '
                'not measured during the photos.')
    return save(fig, 'B2_error_budget')


def b3_data_table():
    fig = canvas('Measured values at each setting', r'Ring-minima phase, $\omega$ from the linear calibration')
    backup_tag(fig, 'B3')
    sp = SPEEDS
    first = FRAMES.groupby('S').time_s.min()
    head_row = ('Motor setting S', 'Photos', 'First photo', 'ω (rad/s)', 'Δφ (rad)', r'Spacing $P$ ($10^5$ px$^2$)')
    xs = [150, 450, 640, 860, 1150, 1440]
    for t, x in zip(head_row, xs):
        fig.text(fx(x), at(215), t, fontsize=14, color=GREY, weight='bold')
    for i, r in sp.iterrows():
        y = 270 + i * 62
        if i % 2 == 0:
            a = px_axes(fig, 130, y - 42, 1670, 60); a.axis('off')
            a.add_patch(Rectangle((0, 0), 1, 1, transform=a.transAxes, fc=LIGHT, ec='none'))
        t = first[r.S]
        u_om = np.hypot(r.u_omega_curve, r.u_omega_repro)
        vals = ['at rest' if r.S == 0 else f'{r.S:.2f}', f'{int(r.n)}',
                f'{int(t // 3600):02d}:{int(t % 3600 // 60):02d}:{int(t % 60):02d}',
                '0' if r.S == 0 else paren(r.omega, u_om), '0 (reference)' if r.S == 0 else paren(r.dphi, r['sem']),
                f'{r.P / 1e5:.3f}']
        for x, v in zip(xs, vals):
            fig.text(fx(x), at(y), v, fontsize=17)
    footer(fig, r'$\Delta\phi$ errors: standard error of the mean of the photos (SEM). All rows share the SEM of the rest '
                'mean (0.02 rad).\n' r'$\omega$ errors combine the calibration line and the 0.48 % '
                'speed reproducibility. Settings were recorded in this order on 11 October 2026.')
    return save(fig, 'B3_data_table')


# ---- backup: drift -----------------------------------------------------------------------
def drift_fit():
    """The joint fit phase = a omega + D t + b (sagnac.py step 8), done as a linear regression like slide 8.
    Returns the parameters, their covariance and the two-sided p value of the drift term D (t test)."""
    sp = SPEEDS
    tm = (FRAMES.groupby('S').time_s.mean() - FRAMES.time_s.min()).to_numpy()
    A = np.c_[sp.omega, tm, np.ones(len(sp))]
    y = sp.phi.to_numpy()
    p, *_ = np.linalg.lstsq(A, y, rcond=None)
    dof = len(y) - 3
    r = y - A @ p
    cov = (r @ r / dof) * np.linalg.inv(A.T @ A)
    t_D = p[1] / np.sqrt(cov[1, 1])
    drift_fit.r2 = 1 - (r @ r) / np.sum((y - y.mean()) ** 2)
    return p, cov, 2 * stats.t.sf(abs(t_D), dof)


def b4_drift():
    fig = canvas('A steady phase drift would mimic part of the slope',
                 'Speeds were photographed in increasing order, with rest photos only at the start')
    backup_tag(fig, 'B4')
    fr = FRAMES.sort_values('time_s')
    a, b = RES['slope'], RES['intercept']
    fr = fr.assign(res=fr.phi_rad - (a * fr.omega + b))
    ax = px_axes(fig, 175, 200, 1020, 270)
    ax.step(fr.minutes, fr.omega, where='post', color=INK, lw=2)
    ax.plot(fr.minutes, fr.omega, 'o', color=INK, ms=3)
    ax.set(ylabel=r'$\omega$ (rad/s)', xlim=(-.3, 12.4), ylim=(-.1, 2.0))
    ax.tick_params(labelbottom=False)
    style(ax)
    ax.set_title('Rate of rotation during the run', loc='left', fontsize=16, color=GREY, pad=10)
    bx = px_axes(fig, 175, 545, 1020, 320)
    bx.axhline(0, color=ORANGE, lw=2)
    for S, g in fr.groupby('S'):
        bx.plot(g.minutes, g.res, 'o', color=BLUE, ms=5, alpha=.7)
        pf = np.polyfit(g.minutes, g.res, 1)
        xx = np.array([g.minutes.min(), g.minutes.max()])
        bx.plot(xx, np.polyval(pf, xx), color=INK, lw=2.2)
    bx.set(xlabel='Minutes since the first photograph', ylabel='Phase − fit (rad)', xlim=(-.3, 12.4),
           ylim=(-.32, .32))
    bx.yaxis.set_major_locator(MultipleLocator(.1))
    style(bx)
    bx.set_title('Each photo’s phase minus the fitted line, with the trend per setting', loc='left', fontsize=16, color=GREY,
                 pad=10)
    p, cov, p_D = drift_fit()
    rho = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    sp = SPEEDS
    x = 1290
    head(fig, x, 215, 'JOINT FIT')
    fig.text(fx(x), at(265), r'$\Delta\phi = a\,\omega + D\,t + b$', fontsize=21)
    fig.text(fx(x), at(320), r'$a$ = ' + pm(p[0], np.sqrt(cov[0, 0])) + ' rad/(rad/s)', fontsize=18)
    fig.text(fx(x), at(358), r'$D$ = ' + pm(p[1] * 1e3, np.sqrt(cov[1, 1]) * 1e3) + ' mrad/s', fontsize=18)
    fig.text(fx(x), at(396), r'correlation$(a, D)$' + f' = {rho:.2f}'.replace('-', '−'), fontsize=18)
    fig.text(fx(x), at(438), r'drift term:  $p$' + f' = {p_D:.2f}', fontsize=17)
    fig.text(fx(x), at(480), r'$R^2$' + f' = {drift_fit.r2:.4f}', fontsize=17)
    fig.text(fx(x), at(522), r'without $D$:  $a$ = ' + pm(RES['slope'], RES['u_slope']), fontsize=16,
             color=GREY)
    head(fig, x, 575, 'WITHIN EACH SETTING')
    fig.text(fx(x), at(620), f'trend {sp.drift.min() * 1e3:.1f} to +{sp.drift.max() * 1e3:.1f} mrad/s'
             .replace('-', '−'), fontsize=18)
    fig.text(fx(x), at(658), f'{int((sp.p_drift < .05).sum())} of {len(sp)} settings: ' + r'$p$ < 0.05', fontsize=18)
    footer(fig, 'At 640 nm SMF-28e+ guides several modes (V number 5.6, above the single-mode limit 2.405). Unless the light '
                'is coupled into one mode,\nwhich is hard to align, it travels in several. Light in several modes is sensitive '
                'to stress and temperature, a plausible source\nof the wander. The quoted uncertainties do not include a steady drift.')
    return save(fig, 'B4_drift')


# ---------------------------------------------------------------------------------------
REGISTRY = {}
for _name, _fn in list(globals().items()):
    if callable(_fn) and (_name[:1] in 'sb') and _name[1:2].isdigit():
        REGISTRY[_name] = _fn

if __name__ == '__main__':
    names = sys.argv[1:] or sorted(REGISTRY)
    for n in names:
        REGISTRY[n]()
