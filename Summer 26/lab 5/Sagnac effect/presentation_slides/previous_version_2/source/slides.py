"""Static Sagnac slide images (1920 x 1080). Run: python slides.py [names...]

Slides 02 and 06 (the animated ones) are made by anim_concept.py and anim_sweep.py.
"""
import sys
from math import floor, log10

import numpy as np
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch
from matplotlib.ticker import MultipleLocator
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


def backup_tag(fig):
    fig.text(1 - 96 / 1920, .878, 'BACKUP', fontsize=13, color=GREY, weight='bold', ha='right',
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
    fig.text(.075, .725, 'Measuring the speed of light with a rotating laser interferometer', fontsize=23, va='top')
    fig.text(.075, .585, 'Daniel Haim Breger  &  Hila Hemo', fontsize=20)
    fig.text(.075, .535, 'Instructor: Ron Nordman', fontsize=17, color=GREY)
    picks = [0.0, 20.74, 44.07, 64.00, 90.0]
    d, gap, x0, y0 = 300, 47, 144, 640
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
    for r in np.linspace(.925, .975, 6):                       # fibre loop, 13 turns drawn as a bundle
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
    ax.annotate('fibre loop, 13 turns', xy=(.95 * np.cos(np.radians(222)), .95 * np.sin(np.radians(222))),
                xytext=(-1.52, -1.13), fontsize=14, color='#3a4652', ha='left', va='center',
                arrowprops=dict(arrowstyle='-', color='#3a4652', lw=1, shrinkA=2, shrinkB=0))
    ax.set_xlim(-1.62, 1.62); ax.set_ylim(-1.32, 1.32)


def s03_setup():
    fig = canvas('Both beams travel the same 13-turn fibre loop',
                 'Top view of the apparatus. Everything drawn here turns with the table.')
    ax = px_axes(fig, 110, 165, 1000, 840)
    draw_setup(ax)
    x0 = 1100
    head(fig, x0, 205, 'SYSTEM PARAMETERS  (given)')
    draw_table(fig, x0, 255, [0, 410], [
        ('Fibre turns  $N$', '13 ± 1'),
        ('Loop radius  $R$', '57.0 ± 0.5 cm'),
        (r'Enclosed area  $A = \pi R^2$', '1.02 m²'),
        (r'Wavelength  $\lambda$', '640 nm'),
        ('Fibre', 'SMF-28e+'),
    ], row_h=42, size=17, width=690, header=False)
    head(fig, x0, 500, 'WHAT WE MEASURE')
    draw_table(fig, x0, 550, [0, 410], [
        (r'Angular velocity  $\omega$', 'stopwatch'),
        (r'Phase shift  $\Delta\phi$', 'camera photos'),
    ], row_h=42, size=17, width=690, header=False)
    head(fig, x0, 690, 'PREDICTION')
    fig.text(fx(x0), at(785), r'$\Delta\phi = \dfrac{8\pi N A}{\lambda c}\,\omega \equiv a_{\rm th}\,\omega$',
             fontsize=29, color=BLUE)
    fig.text(fx(x0), at(880), r'$a_{\rm th}$ = 1.74 ± 0.14 rad/(rad/s) *', fontsize=20)
    footer(fig, 'Schematic after Fig. 1 of the lab manual. P = polarizer, L = coupling lens. The polarizers suppress '
                'back-reflections. ' r'$\lambda$ was given without an uncertainty.' '\n'
                '* Where this uncertainty comes from is discussed with the speed of light.')
    return save(fig, '03_setup')


# ---- 04 ---------------------------------------------------------------------------------
def s04_rings():
    fig = canvas('The rings are evenly spaced in r², not in r',
                 'Fastest-rotation photograph, intensity averaged over the fringe sector around the fitted centre')
    row = representative(90.0)
    ax = px_axes(fig, 150, 185, 500, 500)
    show_photo(ax, row, half=520, label_rings=True)
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
    fig.text(fx(x), at(738), 'Two beams of different wavefront curvature:', fontsize=16, color=GREY)
    fig.text(fx(x), at(838), r'$\phi(r) = 2\pi\,\dfrac{r^2}{P} + \phi_c$', fontsize=26)
    fig.text(fx(x), at(930), r'Rotation shifts $\phi_c$, so the rings move.', fontsize=16)
    fig.text(fx(x), at(963), r'A tilt moves the centre, a refocus changes $P$.', fontsize=16)
    footer(fig, 'Dashed circles: concentric-circle fit to dark-ring minima along 170 rays (−75° to 180°, the lower left '
                'has no fringes). Red channel smoothed with σ = 4 px.\n'
                r'$P = 2\lambda\,/\,|1/R_1 - 1/R_2|$ for wavefront radii $R_1$, $R_2$. Median over all photographs: '
                r'$P$' + f' = {np.median(FRAMES.P_px2) / 1e5:.2f}' + r'$\times10^5$ px$^2$.')
    return save(fig, '04_rings_in_r2')


# ---- 05 ---------------------------------------------------------------------------------
def s05_calibration():
    fig = canvas('A stopwatch calibration converts motor setting to ω',
                 r'Mean turn time $T$ at nine settings, $\omega = 2\pi/T$, with a straight-line fit')
    tab = calibration_table()
    lin = CAL['linear']
    ax = px_axes(fig, 175, 180, 1010, 480)
    for j, S in enumerate(SPEEDS.S[1:]):
        ax.axvline(S, color='#9aa5ae', lw=1.3, ls=(0, (4, 3)), zorder=1,
                   label='Motor settings used for the photos' if j == 0 else None)
    xs = np.linspace(tab.S.min() - 1, tab.S.max() + 2, 100)
    ax.plot(xs, C0 + C1 * xs, color=ORANGE, lw=2.6, label='Linear fit', zorder=2)
    ax.errorbar(tab.S, tab.omega, tab.u, fmt='o', color=BLUE, ms=9, capsize=0, zorder=3, label='Timed turns')
    ax.set(xlim=(0, 95), ylim=(0, 2.0), ylabel=r'Angular velocity $\omega$ (rad/s)')
    ax.tick_params(labelbottom=False)
    style(ax)
    h, l = ax.get_legend_handles_labels()
    ax.legend([h[2], h[1], h[0]], [l[2], l[1], l[0]], loc='upper left', frameon=True, facecolor='white',
              edgecolor='none', framealpha=1)
    res = (tab.omega - (C0 + C1 * tab.S)) / tab.omega * 100
    ures = tab.u / tab.omega * 100
    ar = px_axes(fig, 175, 695, 1010, 190)
    for S in SPEEDS.S[1:]:
        ar.axvline(S, color='#9aa5ae', lw=1.3, ls=(0, (4, 3)), zorder=1)
    ar.axhline(0, color=ORANGE, lw=2)
    ar.errorbar(tab.S, res, ures, fmt='o', color=BLUE, ms=8, capsize=4, zorder=3)
    ar.set(xlim=(0, 95), ylim=(-3.3, 2.2), xlabel=r'Motor setting $S$', ylabel='Residual (%)')
    ar.yaxis.set_major_locator(MultipleLocator(2))
    style(ar)
    rel = tab.u_sw / tab.omega * 100
    x = 1270
    head(fig, x, 225, 'LINEAR FIT')
    fig.text(fx(x), at(280), r'$\omega = kS + \omega_0$', fontsize=22)
    fig.text(fx(x), at(333), r'$k$ = ' + paren(lin['k'], lin['u_k']) + ' rad/s per unit S', fontsize=17)
    fig.text(fx(x), at(375), chi_text(lin['chi2'], lin['dof'], lin['p']), fontsize=17)
    head(fig, x, 470, 'ERROR BARS, RELATIVE')
    fig.text(fx(x), at(515), 'stopwatch timing', fontsize=17)
    fig.text(fx(x + 530), at(515), f'{rel.min():.2f}–{rel.max():.2f} %', fontsize=17, ha='right')
    fig.text(fx(x), at(555), r'table speed scatter $\varepsilon$', fontsize=17)
    fig.text(fx(x + 530), at(555), f'{EPS * 100:.2f} %', fontsize=17, ha='right')
    fig.text(fx(x), at(595), r'(both as $u(\omega)/\omega$)', fontsize=14, color=GREY)
    head(fig, x, 715, 'RESIDUALS')
    fig.text(fx(x), at(755), 'They curve away from the line.', fontsize=16)
    fig.text(fx(x), at(788), 'At the top speed: 2.3 % below.', fontsize=16)
    footer(fig, r'Turn times are stopwatch splits: $u(T) = \mathrm{SD}/\sqrt{n(n+1)}$ for $n$ = 5–8 turns, and '
                r'$u(\omega)/\omega = u(T)/T$. Error bars add $\varepsilon\,\omega$ in quadrature.' + '\n'
                'The linear fit converts every motor setting to ω in this talk.')
    return save(fig, '05_calibration')


# ---- 07 ---------------------------------------------------------------------------------
def ring_line(row):
    R = np.sort(np.array(rings(row)))
    s = R ** 2
    m = np.r_[0, np.cumsum(np.rint(np.diff(s) / 107.5e3))]
    P, s0 = np.polyfit(m, s, 1)
    return m, s, P, s0 / P


def s07_extraction():
    fig = canvas('Squared ring radii give the phase at the centre',
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
    fig.text(fx(x), at(348), r'slope $P$, intercept $\eta$', fontsize=15, color=GREY)
    head(fig, x, 445, 'PHASE AT THE CENTRE')
    fig.text(fx(x), at(510), r'$\phi_c = \pi - 2\pi\eta$', fontsize=23)
    fig.text(fx(x), at(553), r'a dark ring has phase $\pi$', fontsize=15, color=GREY)
    head(fig, x, 650, 'THIS PAIR')
    fig.text(fx(x), at(695), r'$\eta$:', fontsize=19)
    fig.text(fx(x + 45), at(695), f'{etas[0]:.2f}', fontsize=19, color=BLUE, weight='bold')
    fig.text(fx(x + 150), at(695), '→', fontsize=19)
    fig.text(fx(x + 192), at(695), f'{etas[1]:.2f}', fontsize=19, color=ORANGE, weight='bold')
    fig.text(fx(x), at(758), r'$\Delta\phi = 2\pi \times %.2f$' % d_eta, fontsize=23)
    fig.text(fx(x), at(810), r'$\approx$' + f' {2 * np.pi * d_eta:.1f} rad', fontsize=23)
    footer(fig, r'Per photo: dark-ring minima along 170 rays, a concentric-circle fit for the centre and the radii '
                r'$R_m$ (rms %.0f–%.0f px over all photos), then a straight line in $R_m^2$.' % (
                    FRAMES.rms_px.min(), FRAMES.rms_px.max()) + '\n'
                r'$P$ is fitted for every photo, so a change of ring spacing does not move $\phi_c$. Only $R^2/P$ '
                'enters, so no pixel-to-length calibration is needed.')
    return save(fig, '07_phase_extraction')


# ---- 08 ---------------------------------------------------------------------------------
def s08_result():
    fig = canvas('The phase grows linearly with ω, as Sagnac predicts',
                 'Mean phase per motor setting against the calibrated angular velocity')
    sp = SPEEDS
    a, b = RES['slope'], RES['intercept'] - sp.phi.iloc[0]
    u_om = np.hypot(sp.u_omega_curve, sp.u_omega_repro)
    xs = np.linspace(0, 1.95, 100)
    ax = px_axes(fig, 175, 175, 1020, 500)
    th, uth = RES['theory'], RES['u_theory']
    ax.fill_between(xs, (th - uth) * xs, (th + uth) * xs, color=THEORY, alpha=.16, lw=0,
                    label=r'Prediction $\pm1\sigma$')
    ax.plot(xs, th * xs, color=THEORY, lw=2.2, ls=(0, (6, 3)))
    ax.plot(xs, a * xs + b, color=ORANGE, lw=2.8, label='Linear fit')
    ax.errorbar(sp.omega, sp.dphi, sp.u_fit, u_om, fmt='o', color=BLUE, ms=10, mec='white', mew=1.3,
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
    ar.errorbar(sp.omega, sp.dphi - (a * sp.omega + b), sp.u_fit, fmt='o', color=BLUE, ms=8, elinewidth=2,
                capsize=4, zorder=3)
    ar.set(xlim=(-.05, 1.95), ylim=(-.15, .15), xlabel=r'Angular velocity $\omega$ (rad/s)', ylabel='Resid. (rad)')
    ar.yaxis.set_major_locator(MultipleLocator(.1))
    style(ar)
    x = 1345
    head(fig, x, 210, 'MEASURED SLOPE')
    fig.text(fx(x), at(262), f'a = {a:.2f} ± {RES["u_slope"]:.2f}', fontsize=25, color=ORANGE, weight='bold')
    fig.text(fx(x), at(294), 'rad/(rad/s), fit and calibration', fontsize=14, color=GREY)
    head(fig, x, 375, r'PREDICTED  $8\pi^2NR^2/\lambda c$')
    fig.text(fx(x), at(427), r'$\mathbf{a_{th}}$' + f' = {th:.2f} ± {uth:.2f}', fontsize=25, color=THEORY,
             weight='bold')
    head(fig, x, 510, 'AGREEMENT')
    fig.text(fx(x), at(556), r'$a\,/\,a_{\rm th}$' + f' = {RES["ratio"]:.2f} ± {RES["u_ratio"]:.2f}', fontsize=20)
    fig.text(fx(x), at(592), f'{RES["z"]:.1f}σ apart,  ' + r'$p$' + f' = {RES["p_z"]:.2f}', fontsize=17)
    head(fig, x, 675, 'FIT STATISTICS')
    fig.text(fx(x), at(710), 'photo scatter and table speed only', fontsize=14, color=GREY)
    fig.text(fx(x), at(748), chi_text(RES['chi2_without_sigma_x'], RES['dof'], RES['p_without_sigma_x']),
             fontsize=16)
    fig.text(fx(x), at(800), 'so we add a scatter between settings', fontsize=13, color=GREY)
    fig.text(fx(x), at(838), r'$\sigma_x$' + f' = {RES["sigma_x"]:.3f} rad  ' + r'$\rightarrow$  $\chi^2/\mathrm{dof}$ = 1',
             fontsize=16)
    footer(fig, r'Error bars: $u = \sqrt{\mathrm{SEM}^2 + (a\,\varepsilon\,\omega)^2 + \sigma_x^2}$. '
                r'In the upper panel all error bars are smaller than the markers, the residual panel shows them.')
    return save(fig, '08_phase_vs_omega')


# ---- 09 ---------------------------------------------------------------------------------
METHOD_ROWS = [  # label, results key, column in sagnac_speeds.csv, colour, marker, face
    ('Ring minima\n(main)', None, None, BLUE, 'o', BLUE),
    ('Cosine fit\nof the profile', 'cosine fit', 'phi_cos', ORANGE, 'o', ORANGE),
    ('Template match\nto rest profile', 'template matching', 'phi_template', PURPLE, 's', PURPLE),
    ('10 random\npixels', 'random points', 'phi_points', GREY, 'D', 'white'),
]
SHORT = ['Ring minima (main)', 'Cosine fit', 'Template match', '10 random pixels'
]


def s09_methods():
    fig = canvas('Four ways of reading the phase agree within ±3 %',
                 'The same 161 photographs analysed with four phase estimators')
    sp = SPEEDS
    mm = RES['methods']
    ax = px_axes(fig, 175, 215, 760, 540)
    ax.axhline(0, color=BLUE, lw=2.4, label='Ring minima (main)')
    for (lab, key, col, c, mk, mfc), short in zip(METHOD_ROWS[1:], SHORT[1:]):
        d = (sp[col] - sp[col].iloc[0]) - sp.dphi
        ax.plot(sp.omega, d, color=c, lw=1.8, marker=mk, ms=9, mfc=mfc, mec=c, mew=1.8, label=short)
    ax.set(xlim=(-.05, 1.9), ylim=(-.22, .28), xlabel=r'Angular velocity $\omega$ (rad/s)',
           ylabel='Difference from minima (rad)')
    ax.yaxis.set_major_locator(MultipleLocator(.1))
    style(ax)
    ax.legend(loc='upper left', ncol=2, fontsize=14)
    ax.set_title('Phase change relative to the ring-minima result', loc='left', fontsize=16, color=GREY, pad=12)
    bx = px_axes(fig, 1255, 215, 540, 540)
    slopes = [(RES['slope'], RES['u_slope'])] + [(mm[k]['slope'], mm[k]['u_slope']) for _, k, *_ in METHOD_ROWS[1:]]
    th, uth = RES['theory'], RES['u_theory']
    bx.axvspan(th - uth, th + uth, color=THEORY, alpha=.15, lw=0)
    bx.axvline(th, color=THEORY, lw=2, ls=(0, (6, 3)))
    for j, ((lab, key, col, c, mk, mfc), (v, u)) in enumerate(zip(METHOD_ROWS, slopes)):
        y = len(METHOD_ROWS) - 1 - j
        bx.errorbar(v, y, xerr=u, fmt=mk, color=c, ms=11, mfc=mfc, mec=c, mew=2, elinewidth=2.4, capsize=0)
        bx.text(v, y + .27, f'{v:.2f}', ha='center', fontsize=15, color=c, weight='bold')
    lo, hi = min(v for v, _ in slopes), max(v for v, _ in slopes)
    bx.annotate('', xy=(lo, -.5), xytext=(hi, -.5), arrowprops=dict(arrowstyle='|-|', color=INK, lw=1.6,
                                                                     shrinkA=0, shrinkB=0, mutation_scale=6))
    bx.text(.5 * (lo + hi) + .01, -.72, 'range', fontsize=13, ha='center', va='top', color=INK)
    bx.set_yticks(range(len(METHOD_ROWS)))
    bx.set_yticklabels([r[0] for r in METHOD_ROWS[::-1]], fontsize=14)
    bx.set(xlim=(1.55, 1.99), ylim=(-1.05, 3.85), xlabel=r'Fitted slope $a$ (rad/(rad/s))')
    bx.grid(axis='x', alpha=.15)
    bx.set_title(r'Fitted slope  (grey: prediction $\pm1\sigma$)', loc='left', fontsize=15, color=GREY, pad=12)
    fig.text(fx(175), at(925), f'a = {RES["slope"]:.2f} ± {RES["u_slope"]:.2f} (fit + calibration) ± {RES["u_method"]:.2f} '
             f'(method) = 1.84 ± 0.06 rad/(rad/s)', fontsize=19, weight='bold')
    sx = [RES['sigma_x']] + [mm[k]['sigma_x'] for _, k, *_ in METHOD_ROWS[1:]]
    base = [RES['chi2_without_sigma_x']] + method_chi_base()
    footer(fig, 'Each slope comes from the same weighted line as the previous slide, with its own scatter between settings '
                r'$\sigma_x$ = ' + ', '.join(f'{s:.3f}' for s in sx) + ' rad,\nso '
                r'$\chi^2/\mathrm{dof}$ = 1 for each. Without $\sigma_x$: $\chi^2/\mathrm{dof}$ = '
                + ', '.join(f'{c:.0f}/8' for c in base) + r', all $p < 10^{-6}$. Slope errors: fit and calibration.')
    return save(fig, '09_methods')


# ---- 10 ---------------------------------------------------------------------------------
def s10_random_pixels():
    fig = canvas('The slope from ten random pixels depends on the draw',
                 'The pixel method repeated for 31 random draws of ten pixels')
    s1 = np.array([d['slope'] for d in METHODS['points_scan']])
    s2 = np.array([d['slope'] for d in METHODS['points_frame_geometry_scan']])
    ax = px_axes(fig, 175, 210, 1000, 600)
    bins = np.arange(1.35, 2.2, .05)
    ax.hist(s1, bins=bins, color=GREY, alpha=.75, label='fixed pixel positions', edgecolor='white')
    ax.hist(s2, bins=bins, histtype='step', color=PURPLE, lw=2.5,
            label='same pixels, each photo’s own centre and $P$')
    top = .6
    ax.axvline(RES['methods']['random points']['slope'], ymax=top, color=GREY, lw=2.5, ls=(0, (2, 2)),
               label=f'draw used in the comparison: {RES["methods"]["random points"]["slope"]:.2f}')
    ax.axvline(RES['slope'], ymax=top, color=BLUE, lw=3, label=f'ring minima: {RES["slope"]:.2f}')
    ax.axvline(RES['theory'], ymax=top, color=THEORY, lw=2.2, ls=(0, (6, 3)), label=f'prediction: {RES["theory"]:.2f}')
    ax.set(xlabel=r'Fitted slope $a$ (rad/(rad/s))', ylabel='Number of random draws', ylim=(0, 14.5))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    style(ax)
    ax.legend(loc='upper left', fontsize=14)
    x = 1260
    lo, med, hi = np.percentile(s1, [16, 50, 84])
    head(fig, x, 250, 'SPREAD OVER 31 DRAWS')
    fig.text(fx(x), at(297), f'median {med:.2f}', fontsize=19)
    fig.text(fx(x), at(338), f'16–84 %: {lo:.2f} – {hi:.2f}', fontsize=19)
    fig.text(fx(x), at(379), f'full range {s1.min():.2f} – {s2.max() if False else s1.max():.2f}', fontsize=19)
    head(fig, x, 470, 'SHORTCOMINGS')
    for k, t in enumerate(['A fixed pixel reads the phase at', 'its own radius, so a change of', 'ring spacing or a moving',
                           'centre looks like a Sagnac shift.', 'Ten pixels see only half a', 'fringe of motion.']):
        fig.text(fx(x), at(512 + 34 * k), t, fontsize=16)
    footer(fig, 'Each draw: 10 pixels in the fringe sector, read at the same image position in all 161 photos. One '
                'alternating least-squares fit gives every photo’s phase.')
    return save(fig, '10_random_pixels')


# ---- 11 ---------------------------------------------------------------------------------
def s11_speed_of_light():
    fig = canvas('The slope gives c = (2.8 ± 0.2) × 10⁸ m/s',
                 r'The predicted slope solved for $c$, with the measured $a$ and the given $N$, $R$ and $\lambda$')
    c = RES['speed_of_light']
    cv, uc = c['c'] / 1e8, c['u_c'] / 1e8
    um = c['u_measurement'] / 1e8
    x = 175
    fig.text(fx(x), at(262), r'$a = \dfrac{8\pi^2 N R^2}{\lambda\,c}\quad\Rightarrow\quad'
             r'c = \dfrac{8\pi^2 N R^2}{\lambda\,a}$', fontsize=23)
    fig.text(fx(x), at(392), r'$c = \dfrac{8\pi^2\cdot 13\cdot(0.570\,\mathrm{m})^2}{640\,\mathrm{nm}'
             r'\cdot %.3f\,\mathrm{s}} = %.2f\times10^{8}$ m/s' % (RES['slope'], cv), fontsize=19)
    ax = px_axes(fig, 175, 470, 760, 180)
    ax.axvline(C_DEF / 1e8, color=INK, lw=1.8, ls=(0, (5, 3)))
    ax.text(C_DEF / 1e8 + .012, 1.62, r'$c_0$', fontsize=17, color=INK, va='center', ha='left')
    ax.errorbar(cv, 0.55, xerr=uc, fmt='o', color=BLUE, ms=14, elinewidth=3, capsize=8, capthick=2.5, zorder=3)
    ax.errorbar(cv, 0.55, xerr=um, fmt='none', ecolor=ORANGE, elinewidth=9, alpha=.95, zorder=2)
    ax.text(cv, 1.05, f'{cv:.2f} ± {uc:.2f}', fontsize=19, color=BLUE, ha='center', weight='bold')
    ax.set(xlim=(2.45, 3.25), ylim=(0, 2.0), xlabel=r'Speed of light $c$ ($10^8$ m/s)')
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.grid(axis='x', alpha=.15)
    fig.text(fx(x), at(778), '▬ measurement alone (fit, calibration, method): 3.3 %', fontsize=15, color=ORANGE)
    fig.text(fx(x), at(813), r'▬ with the loop geometry $N$ and $R$: 8.5 %', fontsize=15, color=BLUE)
    fig.text(fx(x), at(895), r'$c\,/\,c_0$' + f' = {c["ratio_to_defined"]:.2f} ± {c["u_ratio_to_defined"]:.2f},  '
             f'{abs(c["z"]):.1f}σ,  ' + r'$p$' + f' = {c["p_z"]:.2f}', fontsize=19)
    comps = c['components']
    items = [(r'Turn count $N$ = 13 ± 1', comps['N']), ('Phase-reading method', comps['method']),
             (r'Loop radius $R$ ± 0.5 cm', comps['R']), ('Phase fit', comps['fit']),
             ('Rotation calibration', comps['calibration curve'])]
    bx = px_axes(fig, 1330, 230, 470, 480)
    vals = [v / c['c'] * 100 for _, v in items]
    ys = np.arange(len(items))[::-1] + 1.4
    bx.barh(ys, vals, color=[GREY, ORANGE, GREY, ORANGE, ORANGE], height=.62)
    bx.barh([0], [uc / cv * 100], color=BLUE, height=.62)
    for y, v in zip(list(ys) + [0], vals + [uc / cv * 100]):
        bx.text(v + .15, y, f'{v:.1f} %', va='center', fontsize=15)
    bx.set_yticks(list(ys) + [0])
    bx.set_yticklabels([n for n, _ in items] + ['Total (quadrature)'], fontsize=15)
    bx.get_yticklabels()[-1].set_weight('bold')
    bx.set(xlim=(0, 10), xlabel=r'Relative uncertainty in $c$ (%)')
    bx.grid(axis='x', alpha=.15)
    bx.set_title('Uncertainty budget', loc='left', fontsize=16, color=GREY, pad=12)
    fig.text(fx(1080), at(850), 'The slope also pins down the turn count:', fontsize=16, color=GREY)
    fig.text(fx(1080), at(890), f'N = {RES["N_implied"]:.1f} ± {RES["u_N_implied"]:.1f}  (given 13 ± 1)',
             fontsize=19, weight='bold')
    footer(fig, r'$c_0$ = 299 792 458 m/s (defined). Orange bars come from our measurement, grey bars from the given '
                'loop geometry.\n' r'$\lambda$ = 640 nm was given without an uncertainty. '
                r'This is the slope comparison read as $c$, not an independent measurement of $c$.')
    return save(fig, '11_speed_of_light')


# ---- 12 ---------------------------------------------------------------------------------
def drift_fit():
    """The joint fit phase = a omega + D t + b of sagnac.py (step 8), with its chi2."""
    sp = SPEEDS
    tm = (FRAMES.groupby('S').time_s.mean() - FRAMES.time_s.min()).to_numpy()
    A = np.c_[sp.omega, tm, np.ones(len(sp))]
    u = sp.u_fit.to_numpy(); y = sp.phi.to_numpy()
    W = 1 / u ** 2
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    p = cov @ (A.T @ (W * y))
    chi = float(np.sum(((y - A @ p) / u) ** 2))
    dof = len(y) - 3
    return p, cov, chi, dof, stats.chi2.sf(chi, dof)


def s12_drift():
    fig = canvas('A steady phase drift would mimic part of the slope',
                 'Speeds were photographed in increasing order, with rest photos only at the start')
    fr = FRAMES.sort_values('time_s')
    a, b = RES['slope'], RES['intercept']
    fr = fr.assign(res=fr.phi_rad - (a * fr.omega + b))
    ax = px_axes(fig, 175, 200, 1020, 270)
    ax.step(fr.minutes, fr.omega, where='post', color=INK, lw=2)
    ax.plot(fr.minutes, fr.omega, 'o', color=INK, ms=3)
    ax.set(ylabel=r'$\omega$ (rad/s)', xlim=(-.3, 12.4), ylim=(-.1, 2.0))
    ax.tick_params(labelbottom=False)
    style(ax)
    ax.set_title('Rotation rate during the run', loc='left', fontsize=16, color=GREY, pad=10)
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
    bx.set_title('Residuals from the fitted line, with the trend per setting', loc='left', fontsize=16, color=GREY,
                 pad=10)
    p, cov, chi, dof, pv = drift_fit()
    rho = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    sp = SPEEDS
    x = 1290
    head(fig, x, 215, 'JOINT FIT')
    fig.text(fx(x), at(265), r'$\Delta\phi = a\,\omega + D\,t + b$', fontsize=21)
    fig.text(fx(x), at(320), r'$a$' + f' = {p[0]:.2f} ± {np.sqrt(cov[0, 0]):.2f} rad/(rad/s)', fontsize=18)
    fig.text(fx(x), at(358), r'$D$' + f' = {p[1] * 1e3:.1f} ± {np.sqrt(cov[1, 1]) * 1e3:.1f} mrad/s'.replace('-', '−'),
             fontsize=18)
    fig.text(fx(x), at(396), r'correlation$(a, D)$' + f' = {rho:.2f}'.replace('-', '−'), fontsize=18)
    fig.text(fx(x), at(442), chi_text(chi, dof, pv), fontsize=17)
    fig.text(fx(x), at(490), r'without $D$:  $a$' + f' = {RES["slope"]:.2f} ± {RES["u_slope"]:.2f}', fontsize=16,
             color=GREY)
    head(fig, x, 575, 'WITHIN EACH SETTING')
    fig.text(fx(x), at(620), f'trend {sp.drift.min() * 1e3:.1f} to +{sp.drift.max() * 1e3:.1f} mrad/s'
             .replace('-', '−'), fontsize=18)
    fig.text(fx(x), at(658), f'{int((sp.p_drift < .05).sum())} of {len(sp)} settings: ' + r'$p$ < 0.05', fontsize=18)
    footer(fig, r'SMF-28e+ is not single-mode at 640 nm ($V \approx 5.6 > 2.405$). Coupling between its modes changes '
                'with stress and temperature, a likely source of the wander.\nThe quoted uncertainties do not include a '
                'steady drift.')
    return save(fig, '12_drift')


# ---- 13, 14 -------------------------------------------------------------------------------
def s13_conclusions():
    fig = canvas('The Sagnac phase is linear in ω and matches theory', '')
    items = [('1.84 ± 0.06', 'rad/(rad/s)',
              r'Slope of $\Delta\phi$ against $\omega$, with the method spread.',
              r'Prediction $8\pi^2NR^2/\lambda c$ = 1.74 ± 0.14, 0.7σ apart.'),
             ('3.29 ± 0.06 rad', '0.52 fringe',
              r'Largest shift, at $\omega$ = 1.81 rad/s.',
              'The beams’ travel times differ by only 1.1 fs.'),
             ('(2.8 ± 0.2) × 10⁸ m/s', 'speed of light',
              r'Inferred $c$: 0.7σ from the defined value.',
              r'Limited by the turn count $N$ = 13 ± 1.')]
    y = 215
    for big, unit, l1, l2 in items:
        fig.text(fx(144), at(y + 45), big, fontsize=29, weight='bold', color=BLUE)
        fig.text(fx(144), at(y + 88), unit, fontsize=16, color=GREY)
        fig.text(fx(900), at(y + 28), l1, fontsize=19)
        fig.text(fx(900), at(y + 72), l2, fontsize=17, color=GREY)
        y += 185
    ax = px_axes(fig, 144, 790, 1680, 2)
    ax.axhline(0, color='#c9d1d8', lw=1.5); ax.axis('off')
    head(fig, 144, 850, 'NEXT')
    fig.text(fx(144), at(895), '•  Interleave rest photos and reverse the rotation, to separate a drift from the Sagnac '
             'phase', fontsize=17)
    fig.text(fx(144), at(935), '•  Count the fibre turns (now 13 ± 1)        •  Use a fibre that is single-mode at '
             '640 nm', fontsize=17)
    return save(fig, '13_conclusions')


def s14_thank_you():
    fig = canvas('')
    fig.text(.075, .62, 'Thank you', fontsize=60, weight='bold')
    fig.text(.075, .53, 'Questions?', fontsize=28, color=GREY)
    fig.text(.075, .26, 'Daniel Haim Breger  &  Hila Hemo', fontsize=20)
    fig.text(.075, .21, 'Instructor: Ron Nordman', fontsize=17, color=GREY)
    row = representative(44.07)
    circle_photo(fig, row, 1080, 160, 760, half=460)
    fig.text(fx(1460), at(965), f'ω = {row.omega:.2f} rad/s', fontsize=16, color=GREY, ha='center')
    return save(fig, '14_thank_you')


# ---- backups ----------------------------------------------------------------------------
def b1_quadratic():
    fig = canvas('A quadratic calibration would raise the slope by 1.7 %',
                 'Calibration residuals for the two models, and what each does to the result')
    backup_tag(fig)
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
    rows = [('', 'Linear (used)', 'Quadratic'),
            (r'$\chi^2/\mathrm{dof}$', f'{f1["chi2"]:.1f}/{f1["dof"]}', f'{f2["chi2"]:.1f}/{f2["dof"]}'),
            (r'$p$', fmt_p(f1['p'], bold=True), f'{p2:.2f}'),
            (r'Slope $a$ (rad/(rad/s))', f'{RES["slope"]:.3f}', f'{RES_Q["slope"]:.3f}'),
            (r'$a\,/\,a_{\rm th}$', f'{RES["ratio"]:.3f}', f'{RES_Q["ratio"]:.3f}'),
            (r'$c$ ($10^8$ m/s)', f'{RES["speed_of_light"]["c"] / 1e8:.3f}', f'{RES_Q["speed_of_light"]["c"] / 1e8:.3f}'),
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
    for k, t in enumerate(['The line is the model the manual asks for,', 'but the table turns slower than linear',
                           r'at high settings. The choice moves $a$ and $c$', 'by 1.7 %, inside the 3.3 % uncertainty',
                           'of the measurement.']):
        fig.text(fx(x0), at(700 + 32 * k), t, fontsize=16)
    f3 = CAL['fits_with_eps']['3']
    footer(fig, r'Both fits include the speed reproducibility $\varepsilon$ = 0.48 %, set by $\chi^2/\mathrm{dof}$ = 1 for '
                r'the quadratic, so its $p$ is not an independent test.' '\n'
                'A cubic adds nothing (' + chi_text(f3['chi2'], f3['dof'], stats.chi2.sf(f3['chi2'], f3['dof'])) +
                '). The model-choice shift is reported separately and is not added to the quoted uncertainty.')
    return save(fig, 'B1_quadratic_calibration')


def b2_error_budget():
    fig = canvas('Every error source, and whether it is in the result',
                 r'Sizes refer to the slope $a$. The same relative sizes apply to $c$.')
    backup_tag(fig)
    rows = [('Source', 'Kind', 'Size', 'Included'),
            ('Photo-to-photo scatter', 'statistical', 'SEM 0.004–0.028 rad', 'yes, fit'),
            ('Table speed at a setting', 'statistical (calibration)', r'$\varepsilon$ = 0.48 % of $\omega$', 'yes, fit'),
            ('Wander between settings', 'statistical (fit)', r'$\sigma_x$ = 0.041 rad', 'yes, fit'),
            ('Calibration line', 'calibration, correlated', r'0.2 % of $a$', 'yes'),
            ('Phase-reading method', 'systematic', r'3.0 % of $a$', 'yes'),
            (r'Turn count $N$ = 13 ± 1', 'given', '7.7 %', r'yes (theory, $c$)'),
            (r'Loop radius $R$ = 57.0 ± 0.5 cm', 'given', '1.8 %', r'yes (theory, $c$)'),
            (r'Wavelength $\lambda$ = 640 nm', 'not given', '0.16 % per nm', 'no'),
            ('Steady drift in time', 'systematic', r'not separable from $a$', 'no'),
            ('Linear vs quadratic calibration', 'model choice', r'1.7 % of $a$', 'no, reported')]
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
    footer(fig, r'Statistical parts enter the fit of $\Delta\phi$ against $\omega$. The others are added in quadrature '
                r'afterwards: $a$ = 1.84 ± 0.06, $c$ = (2.8 ± 0.2) × 10⁸ m/s.' '\n'
                r'Not checked: the sign of the shift (rotation direction not recorded). $\omega$ comes from the calibration, '
                'not measured during the photos.')
    return save(fig, 'B2_error_budget')


def b3_data_table():
    fig = canvas('Measured values at each setting', r'Ring-minima phase, $\omega$ from the linear calibration')
    backup_tag(fig)
    sp = SPEEDS
    first = FRAMES.groupby('S').time_s.min()
    head_row = ('Motor setting S', 'Photos', 'Clock', 'ω (rad/s)', 'Δφ (rad)', r'Period $P$ ($10^5$ px$^2$)')
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
        vals = ['at rest' if r.S == 0 else f'{r.S:g}', f'{int(r.n)}',
                f'{int(t // 3600):02d}:{int(t % 3600 // 60):02d}:{int(t % 60):02d}',
                '0' if r.S == 0 else paren(r.omega, u_om), paren(r.dphi, r.u_fit), f'{r.P / 1e5:.3f}']
        for x, v in zip(xs, vals):
            fig.text(fx(x), at(y), v, fontsize=17)
    footer(fig, r'The $\Delta\phi$ errors are the per-setting fit weights. They share the rest reference, so they are not '
                'independent errors of the differences.\n' r'$\omega$ errors combine the calibration line and the 0.48 % '
                'speed reproducibility. Settings were recorded in this order on 11 October 2026.')
    return save(fig, 'B3_data_table')


# ---------------------------------------------------------------------------------------
REGISTRY = {}
for _name, _fn in list(globals().items()):
    if callable(_fn) and (_name[:1] in 'sb') and _name[1:2].isdigit():
        REGISTRY[_name] = _fn

if __name__ == '__main__':
    names = sys.argv[1:] or sorted(REGISTRY)
    for n in names:
        REGISTRY[n]()
