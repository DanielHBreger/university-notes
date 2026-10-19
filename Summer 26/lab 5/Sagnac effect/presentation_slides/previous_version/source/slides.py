"""Static Sagnac slide images (1920 x 1080). Run: python slides.py [names...]"""
import sys

import numpy as np
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch, Polygon, FancyArrowPatch, Arc, Wedge
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from scipy import stats

from deck import *


# ---------------------------------------------------------------------------------------
def s01_title():
    fig = canvas('')
    fig.text(.075, .80, 'The Sagnac effect', fontsize=48, weight='bold')
    fig.text(.075, .725, 'Measuring rotation with counter-propagating light in a 13-turn fibre loop',
             fontsize=23, color=INK, va='top')
    fig.text(.075, .585, 'Daniel Haim Breger  &  Hila Hemo', fontsize=20)
    fig.text(.075, .535, 'Physics laboratory  ·  Technion  ·  measured 11 October 2026', fontsize=15, color=GREY)
    # filmstrip of measured patterns, each centred on its fitted ring centre
    picks = [0.0, 20.74, 44.07, 64.00, 90.0]
    d, gap, x0, y0 = 300, 47, 144, 600
    half = 420
    for j, S in enumerate(picks):
        row = representative(S)
        img, (cx, cy) = crop(photo(row.file), row.cx, row.cy, half)
        yy, xx = np.indices(img.shape[:2])
        img = img.copy(); img[np.hypot(xx - half, yy - half) > half] = 255
        ax = px_axes(fig, x0 + j * (d + gap), y0, d, d)
        ax.imshow(img, extent=(-half, half, half, -half))
        ax.add_patch(Circle((0, 0), half - 1, fill=False, ec='#3b3b3b', lw=1))
        ax.set_xlim(-half - 4, half + 4); ax.set_ylim(half + 4, -half - 4); ax.axis('off')
        lab = 'at rest' if S == 0 else f'{row.omega:.2f} rad/s'
        fig.text((x0 + j * (d + gap) + d / 2) / 1920, 1 - (y0 + d + 38) / 1080, lab, ha='center', fontsize=17,
                 color=BLUE if S == 0 else INK, weight='bold' if S == 0 else 'normal')
    # arrow under the strip
    ax = px_axes(fig, x0, y0 + d + 62, 4 * (d + gap) + d, 30)
    ax.annotate('', xy=(1, .5), xytext=(0, .5), xycoords='axes fraction',
                arrowprops=dict(arrowstyle='-|>,head_length=.9,head_width=.35', color='#9aa5ae', lw=1.6))
    ax.axis('off')
    fig.text(.075, 1 - 1040 / 1080, 'Measured fringe patterns: the rings move inward as the table turns faster',
             fontsize=13, color=GREY)
    return save(fig, '01_title')


# ---------------------------------------------------------------------------------------
def draw_setup(ax):
    """Top view of the apparatus, after Fig. 1 of the lab manual."""
    ax.set_aspect('equal'); ax.axis('off')
    ax.add_patch(Circle((0, 0), 1.0, fc='#f4f6f8', ec='#c3ccd3', lw=1.5, zorder=0))
    # fibre loop: 13 turns drawn as a bundle
    for r in np.linspace(.925, .975, 6):
        ax.add_patch(Circle((0, 0), r, fill=False, ec='#3a4652', lw=.9, alpha=.85, zorder=1))
    # beams
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

    # beam splitter
    s = .085
    ax.add_patch(Rectangle((-s, -s), 2 * s, 2 * s, fc='white', ec=INK, lw=1.6, zorder=5))
    ax.plot([-s, s], [-s, s], color=INK, lw=1.4, zorder=6)
    ax.text(-.13, -.17, 'beam\nsplitter', fontsize=13, ha='right', va='top', color=INK, linespacing=1.1)
    # laser and P1
    ax.add_patch(FancyBboxPatch((-.05, -.86), .10, .17, boxstyle='round,pad=0.005', fc='#e2a400', ec=INK, lw=1.2, zorder=4))
    ax.text(-.08, -.78, 'laser\n640 nm', fontsize=13, ha='right', va='center', color=INK, linespacing=1.1)
    polarizer(0, -.42, vertical=False, label='P₁', dx=.12)
    # fibre couplers (L1 left, L2 top) and the fibre ends at the rim
    lens(-.62, 0, True, 'L₁', dy=-.12)
    lens(0, .62, False, 'L₂', dx=.12)
    for x, y in [(-.86, 0), (0, .86)]:
        ax.add_patch(FancyBboxPatch((x - .045, y - .045), .09, .09, boxstyle='round,pad=0.005', fc='#c9d1d8',
                                    ec=INK, lw=1.2, zorder=4))
    # output: P2, L3, camera
    polarizer(.38, 0, True, 'P₂', dy=-.13)
    lens(.55, 0, True, 'L₃', dy=-.12)
    ax.add_patch(FancyBboxPatch((.66, -.065), .16, .13, boxstyle='round,pad=0.01', fc=INK, ec=INK, zorder=4))
    ax.add_patch(Rectangle((.62, -.035), .05, .07, fc=INK, ec=INK, zorder=4))
    ax.text(.74, .12, 'camera', fontsize=13, ha='center', color=INK)
    # direction of the two beams in the loop (outside the rim, upper left)
    for ang0, ang1, col, r in [(115, 165, BLUE, 1.09), (165, 115, ORANGE, 1.18)]:
        th = np.radians(np.linspace(ang0, ang1, 60))
        ax.plot(r * np.cos(th), r * np.sin(th), color=col, lw=3, solid_capstyle='round')
        ax.annotate('', xy=(r * np.cos(th[-1]), r * np.sin(th[-1])),
                    xytext=(r * np.cos(th[-6]), r * np.sin(th[-6])),
                    arrowprops=dict(arrowstyle='-|>,head_length=.7,head_width=.4', color=col, lw=3))
    ax.text(-1.17, .12, 'beam 1', color=BLUE, fontsize=15, ha='center', va='center', weight='bold')
    ax.text(-.62, 1.17, 'beam 2', color=ORANGE, fontsize=15, ha='right', va='center', weight='bold')
    # rotation
    th = np.radians(np.linspace(20, 70, 60))
    ax.plot(1.12 * np.cos(th), 1.12 * np.sin(th), color=INK, lw=2.2)
    ax.annotate('', xy=(1.12 * np.cos(th[-1]), 1.12 * np.sin(th[-1])),
                xytext=(1.12 * np.cos(th[-5]), 1.12 * np.sin(th[-5])),
                arrowprops=dict(arrowstyle='-|>,head_length=.8,head_width=.45', color=INK, lw=2.2))
    ax.text(1.03, .96, 'ω', fontsize=24, color=INK, ha='left', va='bottom', style='italic')
    ax.annotate('fibre loop, 13 turns', xy=(.95 * np.cos(np.radians(222)), .95 * np.sin(np.radians(222))),
                xytext=(-1.52, -1.13), fontsize=14, color='#3a4652', ha='left', va='center',
                arrowprops=dict(arrowstyle='-', color='#3a4652', lw=1, shrinkA=2, shrinkB=0))
    ax.set_xlim(-1.62, 1.62); ax.set_ylim(-1.32, 1.32)


def s03_setup():
    fig = canvas('Both beams travel the same 13-turn fibre loop',
                 'Top view of the apparatus; everything drawn here turns with the table')
    ax = px_axes(fig, 110, 165, 1000, 840)
    draw_setup(ax)
    x = 1100 / 1920
    def head(y, txt):
        fig.text(x, 1 - y / 1080, txt, fontsize=14, color=GREY, weight='bold')
    def line(y, txt, size=19, color=INK):
        fig.text(x, 1 - y / 1080, txt, fontsize=size, color=color, va='baseline')
    head(215, 'LOOP AND LIGHT  (given)')
    line(270, 'N = 13 ± 1 turns,  R = 57.0 ± 0.5 cm')
    line(312, 'A = πR² = 1.02 m²,  λ = 640 nm')
    head(395, 'WHAT WE MEASURE')
    line(450, 'ω   stopwatch times of table turns')
    line(492, 'Δφ  photos of the output fringes')
    line(530, '      161 photos: at rest and at 9 speeds', 16, GREY)
    head(615, 'PREDICTION')
    fig.text(x, 1 - 722 / 1080, r'$\Delta\phi = \dfrac{8\pi N A}{\lambda c}\,\omega \equiv a_{\rm th}\,\omega$',
             fontsize=29, color=BLUE)
    line(832, r'$a_{\rm th}$ = 1.74 ± 0.14 rad per (rad/s)', 20)
    line(870, 'the ± comes almost entirely from N = 13 ± 1', 15, GREY)
    footer(fig, 'Schematic after Fig. 1 of the lab manual. P = polarizer, L = coupling lens. The polarizers suppress '
                'back-reflections from the fibre ends.\nThe small off-centre loop (46 turns, R = 15.5 cm) was not '
                'measured. λ was given without an uncertainty.')
    return save(fig, '03_setup')


# ---------------------------------------------------------------------------------------
def show_photo(ax, row, half=470, ring_color='white', label_rings=False, cross=True):
    img, (cx, cy) = crop(photo(row.file), row.cx, row.cy, half)
    ax.imshow(img, extent=(-half, half, half, -half))
    th = np.linspace(0, 2 * np.pi, 400)
    for k, R in enumerate(rings(row)):
        ax.plot(R * np.cos(th), R * np.sin(th), color=ring_color, lw=1.6, ls=(0, (5, 3)), alpha=.9)
        if label_rings:
            a = np.radians(50)
            ax.text(R * np.cos(a) + 6, -R * np.sin(a) - 6, f'm = {k}', color='white', fontsize=13,
                    ha='left', va='bottom', weight='bold')
    if cross:
        ax.plot(0, 0, '+', color='white', ms=14, mew=1.8)
    ax.set_xlim(-half, half); ax.set_ylim(half, -half)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def s04_rings():
    fig = canvas('The rings are evenly spaced in r², not in r',
                 'Fastest-rotation photograph; intensity averaged over the fringe sector around the fitted centre')
    row = representative(90.0)
    ax = px_axes(fig, 150, 185, 500, 500)
    show_photo(ax, row, half=520, label_rings=True)
    R = np.array(rings(row))
    r, prof = sector_profile_r(row.file, row.cx, row.cy, r_max=560)
    # r
    a1 = px_axes(fig, 990, 215, 810, 255)
    a1.plot(r, prof, color=BLUE, lw=2.2)
    for x in R:
        a1.axvline(x, color=INK, lw=1, ls=(0, (4, 3)), alpha=.7)
    for x0, x1 in zip(R[:-1], R[1:]):
        a1.annotate('', xy=(x1, 190), xytext=(x0, 190), arrowprops=dict(arrowstyle='<->', color=ORANGE, lw=1.8,
                                                                        shrinkA=0, shrinkB=0))
        a1.text(.5 * (x0 + x1), 197, f'{x1 - x0:.0f} px', color=ORANGE, ha='center', va='bottom', fontsize=15,
                weight='bold')
    a1.set(xlim=(0, 560), ylim=(0, 225), xlabel='Radius r (px)', ylabel='Intensity (counts)')
    style(a1)
    a1.set_title('Against r: the spacing shrinks outward', loc='left', fontsize=16, color=GREY, pad=14)
    # r^2
    a2 = px_axes(fig, 990, 625, 810, 255)
    s = r ** 2 / 1e5
    a2.plot(s, prof, color=BLUE, lw=2.2)
    for x in R ** 2 / 1e5:
        a2.axvline(x, color=INK, lw=1, ls=(0, (4, 3)), alpha=.7)
    S2 = R ** 2 / 1e5
    for x0, x1 in zip(S2[:-1], S2[1:]):
        a2.annotate('', xy=(x1, 190), xytext=(x0, 190), arrowprops=dict(arrowstyle='<->', color=ORANGE, lw=1.8,
                                                                        shrinkA=0, shrinkB=0))
        a2.text(.5 * (x0 + x1), 197, f'{x1 - x0:.2f}', color=ORANGE, ha='center', va='bottom', fontsize=15,
                weight='bold')
    a2.set(xlim=(0, 3.14), ylim=(0, 225), xlabel=r'Squared radius r² (10⁵ px²)', ylabel='Intensity (counts)')
    style(a2)
    a2.set_title('Against r²: one constant spacing P', loc='left', fontsize=16, color=GREY, pad=14)
    # physics under the photo
    x = 150 / 1920
    fig.text(x, 1 - 738 / 1080, 'Two beams of different wavefront curvature:', fontsize=16, color=GREY)
    fig.text(x, 1 - 838 / 1080, r'$\phi(r) = 2\pi\,\dfrac{r^2}{P} + \phi_c$', fontsize=26)
    fig.text(x, 1 - 930 / 1080, 'Rotation shifts φc, so the rings move.', fontsize=16, color=INK)
    fig.text(x, 1 - 963 / 1080, 'Tilt moves the centre; focus changes P.', fontsize=16, color=INK)
    footer(fig, f'Dashed circles: concentric-circle fit to dark-ring minima along 170 rays (−75° to 180°; the lower '
                f'left has no fringes). Red channel smoothed with σ = 4 px.\nP = 2λ / |1/R₁ − 1/R₂| for wavefront radii '
                f'R₁, R₂; median over all photographs P = {np.median(FRAMES.P_px2) / 1e5:.2f} × 10⁵ px².')
    return save(fig, '04_rings_in_r2')


def s05_calibration():
    fig = canvas('A stopwatch calibration converts motor setting to ω',
                 'Mean turn time at nine settings, ω = 2π / T; straight-line fit as in the lab manual')
    tab = calibration_table()
    xs = np.linspace(0, 95, 200)
    ax = px_axes(fig, 175, 180, 1010, 480)
    ax.errorbar(tab.S, tab.omega, tab.u, fmt='o', color=BLUE, ms=9, capsize=0, zorder=3, label='Timed turns')
    ax.plot(xs, C0 + C1 * xs, color=ORANGE, lw=2.6, label='Linear fit', zorder=2)
    ax.set(xlim=(0, 95), ylim=(-.05, 2.0), ylabel='Angular velocity ω (rad/s)')
    ax.tick_params(labelbottom=False)
    style(ax)
    ax.legend(loc='upper left')
    # photographed settings
    ax.plot(SPEEDS.S[1:], np.full(9, -.05), ls='none', marker='^', color=GREY, ms=9, clip_on=False, zorder=4,
            label='Settings used for photos')
    ax.legend(loc='upper left')
    # residuals
    res = (tab.omega - (C0 + C1 * tab.S)) / tab.omega * 100
    ures = tab.u / tab.omega * 100
    ar = px_axes(fig, 175, 695, 1010, 190)
    ar.axhline(0, color=ORANGE, lw=2)
    ar.errorbar(tab.S, res, ures, fmt='o', color=BLUE, ms=8, capsize=4, zorder=3)
    ar.set(xlim=(0, 95), ylim=(-3.3, 2.2), xlabel='Motor setting S', ylabel='Residual (%)')
    ar.yaxis.set_major_locator(MultipleLocator(2))
    style(ar)
    # numbers on the right
    lin = CAL['linear']
    x = 1290 / 1920
    def at(y):
        return 1 - y / 1080
    fig.text(x, at(225), 'FIT   ω = kS + ω₀', fontsize=14, color=GREY, weight='bold')
    fig.text(x, at(275), f'k = {lin["k"]:.5f} ± {lin["u_k"]:.5f}', fontsize=19)
    fig.text(x, at(306), 'rad/s per unit of S', fontsize=14, color=GREY)
    fig.text(x, at(355), f'ω₀ = {lin["omega0"]:.4f} ± {lin["u_omega0"]:.4f} rad/s'.replace('-', '−'), fontsize=18)
    fig.text(x, at(450), 'PRECISION', fontsize=14, color=GREY, weight='bold')
    fig.text(x, at(497), 'Stopwatch: 0.04–0.18 % of ω', fontsize=18)
    fig.text(x, at(532), f'Speed reproducibility: {EPS * 100:.2f} %', fontsize=18)
    fig.text(x, at(630), 'CURVATURE', fontsize=14, color=GREY, weight='bold')
    fig.text(x, at(677), 'At the top speed the table', fontsize=18)
    fig.text(x, at(712), 'runs 2.3 % below the line', fontsize=18)
    fig.text(x, at(747), f'(χ² = {lin["chi2"]:.0f} for {lin["dof"]} dof).', fontsize=18)
    fig.text(x, at(800), 'Quadratic calibration: backup slide', fontsize=15, color=GREY)
    footer(fig, 'Turn times are stopwatch splits (shared endpoints): u(T) = SD / √(n(n+1)), n = 5–8 turns. Error bars add '
                'the 0.48 % speed reproducibility in quadrature\n(set by the scatter about a smooth curve). The linear fit '
                'converts every motor setting to ω in this talk.')
    return save(fig, '05_calibration')


# ---------------------------------------------------------------------------------------
def ring_line(row):
    R = np.sort(np.array(rings(row)))
    s = R ** 2
    m = np.r_[0, np.cumsum(np.rint(np.diff(s) / 107.5e3))]
    P, s0 = np.polyfit(m, s, 1)
    return m, s, P, s0 / P


def s07_extraction():
    fig = canvas('Squared ring radii give the phase at the centre',
                 'One photograph at rest and one at the highest speed; the same fit is applied to all 161')
    rest, fast = representative(0.0), representative(90.0)
    for row, y, col, lab in [(rest, 215, BLUE, 'at rest'), (fast, 610, ORANGE, f'ω = {fast.omega:.2f} rad/s')]:
        ax = px_axes(fig, 150, y, 320, 320)
        show_photo(ax, row, half=520, label_rings=False)
        for s in ax.spines.values():
            s.set_visible(True); s.set_color(col); s.set_linewidth(4)
        fig.text(150 / 1920, 1 - (y - 14) / 1080, lab, color=col, fontsize=17, weight='bold')
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
    ax.set(xlim=(-.8, 2.3), ylim=(-.25, 3.05), xlabel='Dark-ring number m', ylabel=r'$R_m^2\,/\,P$   (fringe periods)')
    ax.set_xticks([0, 1, 2]); ax.yaxis.set_major_locator(MultipleLocator(1))
    style(ax)
    # recipe on the right
    x = 1450 / 1920
    def at(y):
        return 1 - y / 1080
    fig.text(x, at(240), 'FIT EACH PHOTO', fontsize=14, color=GREY, weight='bold')
    fig.text(x, at(305), r'$R_m^2 = P\,(m + \eta)$', fontsize=23)
    fig.text(x, at(348), 'slope P, intercept η', fontsize=15, color=GREY)
    fig.text(x, at(445), 'PHASE AT THE CENTRE', fontsize=14, color=GREY, weight='bold')
    fig.text(x, at(510), r'$\phi_c = \pi - 2\pi\eta$', fontsize=23)
    fig.text(x, at(553), 'a dark ring has phase π', fontsize=15, color=GREY)
    fig.text(x, at(650), 'THIS PAIR', fontsize=14, color=GREY, weight='bold')
    fig.text(x, at(695), 'η:', fontsize=19, color=INK)
    fig.text(x + 45 / 1920, at(695), f'{etas[0]:.2f}', fontsize=19, color=BLUE, weight='bold')
    fig.text(x + 150 / 1920, at(695), '→', fontsize=19, color=INK)
    fig.text(x + 192 / 1920, at(695), f'{etas[1]:.2f}', fontsize=19, color=ORANGE, weight='bold')
    fig.text(x, at(758), r'$\Delta\phi = 2\pi \times %.2f$' % d_eta, fontsize=23, color=INK)
    fig.text(x, at(810), f'= {2 * np.pi * d_eta:.2f} rad', fontsize=23, color=INK)
    footer(fig, 'Per photo: dark-ring minima along 170 rays → concentric-circle fit (common centre, radii Rₘ) → straight '
                'line in Rₘ². P is fitted per photo, so a change in ring\nspacing does not move φc, and only R²/P enters: '
                'no pixel-to-length calibration. This is the manual’s Eq. (7) (minimum shift over the local period) in r².')
    return save(fig, '07_phase_extraction')


def s08_result():
    fig = canvas('The phase grows linearly with ω, as Sagnac predicts',
                 'Mean phase per setting (10–22 photos each) against the calibrated angular velocity')
    sp = SPEEDS
    a, b = RES['slope'], RES['intercept'] - sp.phi.iloc[0]
    u_om = np.hypot(sp.u_omega_curve, sp.u_omega_repro)
    xs = np.linspace(0, 1.95, 100)
    ax = px_axes(fig, 175, 175, 1020, 500)
    th, uth = RES['theory'], RES['u_theory']
    ax.fill_between(xs, (th - uth) * xs, (th + uth) * xs, color=THEORY, alpha=.16, lw=0,
                    label='Prediction ± 1σ')
    ax.plot(xs, th * xs, color=THEORY, lw=2.2, ls=(0, (6, 3)))
    ax.plot(xs, a * xs + b, color=ORANGE, lw=2.8, label='Linear fit')
    ax.errorbar(sp.omega, sp.dphi, sp.u_fit, u_om, fmt='o', color=BLUE, ms=10, mec='white', mew=1.3,
                elinewidth=2, capsize=0, zorder=4, label='Measured')
    ax.set(xlim=(-.05, 1.95), ylim=(-.25, 3.6), ylabel='Phase shift Δφ (rad)')
    ax.tick_params(labelbottom=False)
    style(ax)
    h, l = ax.get_legend_handles_labels()
    ax.legend([h[2], h[1], h[0]], [l[2], l[1], l[0]], loc='upper left')
    pr = ax.twinx()
    pr.set_ylim(np.array(ax.get_ylim()) / (2 * np.pi)); pr.set_ylabel('fringes', color=GREY)
    pr.tick_params(colors=GREY); pr.spines['right'].set_visible(True); pr.spines['right'].set_color('#85919b')
    pr.yaxis.set_major_locator(MultipleLocator(.1))
    ar = px_axes(fig, 175, 700, 1020, 160)
    ar.axhline(0, color=ORANGE, lw=2)
    ar.errorbar(sp.omega, sp.dphi - (a * sp.omega + b), sp.u_fit, fmt='o', color=BLUE, ms=8, elinewidth=2,
                capsize=0, zorder=3)
    ar.set(xlim=(-.05, 1.95), ylim=(-.15, .15), xlabel='Angular velocity ω (rad/s)', ylabel='Resid. (rad)')
    ar.yaxis.set_major_locator(MultipleLocator(.1))
    style(ar)
    x = 1355 / 1920
    def at(y):
        return 1 - y / 1080
    fig.text(x, at(210), 'MEASURED SLOPE', fontsize=14, color=GREY, weight='bold')
    fig.text(x, at(262), f'a = {a:.2f} ± {RES["u_slope"]:.2f}', fontsize=25, color=ORANGE, weight='bold')
    fig.text(x, at(294), 'rad per rad/s  (fit and calibration)', fontsize=14, color=GREY)
    fig.text(x, at(375), 'PREDICTED  8π²NR² / λc', fontsize=14, color=GREY, weight='bold')
    fig.text(x, at(427), r'$\mathbf{a_{th}}$' + f' = {th:.2f} ± {uth:.2f}', fontsize=25, color=THEORY, weight='bold')
    fig.text(x, at(510), 'AGREEMENT', fontsize=14, color=GREY, weight='bold')
    fig.text(x, at(556), r'$a\,/\,a_{\rm th}$' + f' = {RES["ratio"]:.2f} ± {RES["u_ratio"]:.2f}', fontsize=20)
    fig.text(x, at(592), f'{RES["z"]:.1f} σ apart  (p = {RES["p_z"]:.2f})', fontsize=17)
    fig.text(x, at(675), 'LINEARITY', fontsize=14, color=GREY, weight='bold')
    fig.text(x, at(720), f'quadratic term {RES["curvature_q"]:.2f} ± {RES["u_q"]:.2f}', fontsize=17)
    fig.text(x, at(752), f'{RES["curvature_q"] / RES["u_q"]:.1f} σ from zero', fontsize=17)
    fig.text(x, at(835), 'LARGEST SHIFT', fontsize=14, color=GREY, weight='bold')
    fig.text(x, at(878), '3.29 ± 0.06 rad = 0.52 fringe', fontsize=16)
    footer(fig, f'Error bars: √(SEM² + (a ε ω)² + σx²), with σx = {RES["sigma_x"]:.3f} rad of extra scatter between '
                f'settings set by χ²/dof = 1 (without it χ² = {RES["chi2_without_sigma_x"]:.0f} for 8 dof). Horizontal\n'
                'errors (≤ 0.009 rad/s) are smaller than the markers. The fit uses absolute phases with a free intercept; '
                'the rest phase is subtracted for display.')
    return save(fig, '08_phase_vs_omega')


def s09_methods():
    fig = canvas('Three ways of reading the phase agree to ±3 %',
                 'The same 161 photographs analysed with different phase estimators')
    sp = SPEEDS
    mm = RES['methods']
    ax = px_axes(fig, 175, 215, 760, 560)
    ax.axhline(0, color=BLUE, lw=2.4, label='Ring minima (main)')
    for col, lab, mk, mfc, c in [('phi_cos', 'Cosine fit', 'o', ORANGE, ORANGE),
                                 ('phi_template', 'Template match', 's', PURPLE, PURPLE),
                                 ('phi_points', '10 fixed pixels', 'D', 'white', GREY)]:
        d = (sp[col] - sp[col].iloc[0]) - sp.dphi
        ax.plot(sp.omega, d, color=c, lw=1.8, marker=mk, ms=9, mfc=mfc, mec=c, mew=1.8, label=lab,
                alpha=1 if col != 'phi_points' else .8)
    ax.set(xlim=(-.05, 1.9), ylim=(-.22, .28), xlabel='Angular velocity ω (rad/s)',
           ylabel='Difference from minima (rad)')
    ax.yaxis.set_major_locator(MultipleLocator(.1))
    style(ax)
    ax.legend(loc='upper left', ncol=2, fontsize=14)
    ax.set_title('Phase change relative to the ring-minima result', loc='left', fontsize=16, color=GREY, pad=12)
    bx = px_axes(fig, 1255, 215, 540, 560)
    rows = [('Ring minima\n(main)', RES['slope'], RES['u_slope'], BLUE, 'o', BLUE),
            ('Cosine fit\nof the profile', mm['cosine fit']['slope'], mm['cosine fit']['u_slope'], ORANGE, 'o', ORANGE),
            ('Template match\nto rest profile', mm['template matching']['slope'], mm['template matching']['u_slope'],
             PURPLE, 's', PURPLE),
            ('10 fixed pixels\n(not used)', mm['random points']['slope'], mm['random points']['u_slope'], GREY, 'D',
             'white')]
    lo, hi = min(r[1] for r in rows[:3]), max(r[1] for r in rows[:3])
    bx.axvspan(lo, hi, color=BLUE, alpha=.10, lw=0)
    bx.axvline(RES['theory'], color=THEORY, lw=2, ls=(0, (6, 3)))
    for j, (lab, v, u, c, mk, mfc) in enumerate(rows):
        y = len(rows) - 1 - j
        bx.errorbar(v, y, xerr=u, fmt=mk, color=c, ms=11, mfc=mfc, mec=c, mew=2, elinewidth=2.4, capsize=0)
        bx.text(v, y + .27, f'{v:.2f}', ha='center', fontsize=15, color=c, weight='bold')
    bx.set_yticks(range(len(rows)))
    bx.set_yticklabels([r[0] for r in rows[::-1]], fontsize=14)
    bx.set(xlim=(1.55, 1.99), ylim=(-.6, 3.7), xlabel='Fitted slope a (rad per rad/s)')
    bx.text(RES['theory'] - .01, 3.55, 'prediction', color=THEORY, fontsize=13, ha='right', va='center')
    bx.text(.5 * (lo + hi), -.45, 'spread', color=BLUE, fontsize=13, ha='center', va='center')
    bx.grid(axis='x', alpha=.15)
    bx.set_title('Fitted slope', loc='left', fontsize=16, color=GREY, pad=12)
    fig.text(175 / 1920, 1 - 905 / 1080, f'Method uncertainty = half the spread of the three used slopes = '
             f'{RES["u_method"]:.3f}', fontsize=17, color=INK)
    fig.text(175 / 1920, 1 - 945 / 1080, f'→  a = {RES["slope"]:.2f} ± {RES["u_slope"]:.2f} (fit) ± '
             f'{RES["u_method"]:.2f} (method) = 1.84 ± 0.06 rad per rad/s', fontsize=17, color=INK, weight='bold')
    footer(fig, 'Cosine fit: envelope × cosine model of each r² profile (fits worst at rest). Template: shift of each profile '
                'against the mean rest profile. All read the same photos.\n10 fixed pixels: left out, because the slope '
                'depends on which pixels are drawn (16–84 % over 31 draws: 1.51–1.90). Slope errors: fit and calibration.')
    return save(fig, '09_methods')


def s10_speed_of_light():
    fig = canvas('The slope gives c = (2.8 ± 0.2) × 10⁸ m/s',
                 'The Sagnac relation solved for c, with the given N, R and λ')
    c = RES['speed_of_light']
    cv, uc = c['c'] / 1e8, c['u_c'] / 1e8
    um = c['u_measurement'] / 1e8
    ax = px_axes(fig, 175, 230, 760, 270)
    ax.axvline(C_DEF / 1e8, color=INK, lw=1.8, ls=(0, (5, 3)))
    ax.text(C_DEF / 1e8 + .012, 1.6, 'defined c', fontsize=14, color=INK, va='center', ha='left')
    ax.errorbar(cv, 0.6, xerr=uc, fmt='o', color=BLUE, ms=14, elinewidth=3, capsize=8, capthick=2.5, zorder=3)
    ax.errorbar(cv, 0.6, xerr=um, fmt='none', ecolor=ORANGE, elinewidth=9, alpha=.95, zorder=2)
    ax.text(cv, 1.0, f'{cv:.2f} ± {uc:.2f}', fontsize=20, color=BLUE, ha='center', weight='bold')
    ax.set(xlim=(2.45, 3.25), ylim=(0, 2.0), xlabel='Speed of light c (10⁸ m/s)')
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.grid(axis='x', alpha=.15)
    x = 175 / 1920
    fig.text(x, 1 - 615 / 1080, '▬ measurement alone (fit, calibration, method): 3.3 %', fontsize=15, color=ORANGE)
    fig.text(x, 1 - 650 / 1080, '▬ with the loop geometry N and R: 8.5 %', fontsize=15, color=BLUE)
    fig.text(x, 1 - 775 / 1080, r'$c = \dfrac{8\pi^2 N R^2}{\lambda\,a}$', fontsize=28)
    fig.text(x, 1 - 900 / 1080, f'c / c_defined = {c["ratio_to_defined"]:.2f} ± {c["u_ratio_to_defined"]:.2f}'
             f'  ({abs(c["z"]):.1f} σ, p = {c["p_z"]:.2f})', fontsize=18)
    comps = c['components']
    items = [('Turn count N = 13 ± 1', comps['N']), ('Phase-reading method', comps['method']),
             ('Loop radius R ± 0.5 cm', comps['R']), ('Phase fit', comps['fit']),
             ('Rotation calibration', comps['calibration curve'])]
    bx = px_axes(fig, 1330, 230, 470, 480)
    vals = [v / c['c'] * 100 for _, v in items]
    ys = np.arange(len(items))[::-1] + 1.4
    cols = [GREY, ORANGE, GREY, ORANGE, ORANGE]
    bx.barh(ys, vals, color=cols, height=.62)
    bx.barh([0], [uc / cv * 100], color=BLUE, height=.62)
    for y, v in zip(list(ys) + [0], vals + [uc / cv * 100]):
        bx.text(v + .15, y, f'{v:.1f} %', va='center', fontsize=15)
    bx.set_yticks(list(ys) + [0])
    bx.set_yticklabels([n for n, _ in items] + ['Total (quadrature)'], fontsize=15)
    bx.get_yticklabels()[-1].set_weight('bold')
    bx.set(xlim=(0, 10), xlabel='Relative uncertainty in c (%)')
    bx.grid(axis='x', alpha=.15)
    bx.set_title('Uncertainty budget', loc='left', fontsize=16, color=GREY, pad=12)
    fig.text(1080 / 1920, 1 - 850 / 1080, 'The slope also pins down the turn count:', fontsize=16, color=GREY)
    fig.text(1080 / 1920, 1 - 890 / 1080, f'N = {RES["N_implied"]:.1f} ± {RES["u_N_implied"]:.1f}  (given 13 ± 1)',
             fontsize=19, color=INK, weight='bold')
    footer(fig, 'Orange bars: from our measurement; grey bars: from the given loop geometry. λ = 640 nm was given without '
                'an uncertainty (1 nm changes c by 0.16 %).\nThis is the slope comparison read as c, not an '
                'independent measurement of c.')
    return save(fig, '10_speed_of_light')


def s11_drift():
    fig = canvas('A steady phase drift would mimic part of the slope',
                 'Speeds were photographed in increasing order, with rest photos only at the start')
    fr = FRAMES.sort_values('time_s')
    a, b = RES['slope'], RES['intercept']
    fr = fr.assign(res=fr.phi_rad - (a * fr.omega + b))
    ax = px_axes(fig, 175, 200, 1020, 270)
    ax.step(fr.minutes, fr.omega, where='post', color=INK, lw=2)
    ax.plot(fr.minutes, fr.omega, 'o', color=INK, ms=3)
    ax.set(ylabel='ω (rad/s)', xlim=(-.3, 12.4), ylim=(-.1, 2.0))
    ax.tick_params(labelbottom=False)
    style(ax)
    ax.set_title('Rotation rate during the run', loc='left', fontsize=16, color=GREY, pad=10)
    bx = px_axes(fig, 175, 545, 1020, 320)
    bx.axhline(0, color=ORANGE, lw=2)
    for S, g in fr.groupby('S'):
        bx.plot(g.minutes, g.res, 'o', color=BLUE, ms=5, alpha=.7)
        p = np.polyfit(g.minutes, g.res, 1)
        xx = np.array([g.minutes.min(), g.minutes.max()])
        bx.plot(xx, np.polyval(p, xx), color=INK, lw=2.2)
    bx.set(xlabel='Minutes since the first photograph', ylabel='Phase − fit (rad)', xlim=(-.3, 12.4),
           ylim=(-.32, .32))
    bx.yaxis.set_major_locator(MultipleLocator(.1))
    style(bx)
    bx.set_title('Residuals from the fitted line, with the trend per setting',
                 loc='left', fontsize=16, color=GREY, pad=10)
    dc = RES['drift_check']
    x = 1280 / 1920
    def at(y):
        return 1 - y / 1080
    lines = [(215, 'WHAT WE SEE', 'h'),
             (257, 'Within a setting the phase drifts', ''),
             (289, 'by −7 to +8 mrad/s; between', ''),
             (321, f'settings it wanders by {RES["sigma_x"]:.2f} rad.', ''),
             (408, 'WHY IT MATTERS', 'h'),
             (450, 'A drift steady in time looks just', ''),
             (482, 'like a slope in ω. Fit both:', ''),
             (526, f'a = {dc["slope"]:.2f} ± {dc["u_slope"]:.2f} rad per rad/s', 'b'),
             (560, f'D = {dc["drift"] * 1e3:.1f} ± {dc["u_drift"] * 1e3:.1f} mrad/s'.replace('-', '−'), 'b'),
             (596, f'correlation {dc["correlation"]:.2f}: not separable'.replace('-', '−'), ''),
             (683, 'NEXT TIME', 'h'),
             (725, '• rest photos between speeds', ''),
             (759, '• reverse the rotation: Sagnac', ''),
             (791, '   flips sign, a drift does not', '')]
    for y, t, k in lines:
        if k == 'h':
            fig.text(x, at(y), t, fontsize=14, color=GREY, weight='bold')
        else:
            fig.text(x, at(y), t, fontsize=16.5, weight='bold' if k == 'b' else 'normal')
    footer(fig, 'Possible cause (not tested): SMF-28e+ carries a few modes at 640 nm (V ≈ 5.6), so stress- and '
                'temperature-dependent mode coupling can add\na slowly varying phase. The quoted uncertainties do not '
                'include a steady drift.')
    return save(fig, '11_drift')


def s12_conclusions():
    fig = canvas('The Sagnac phase is linear in ω and matches theory', '')
    items = [('1.84 ± 0.06', 'rad per rad/s',
              'Slope of Δφ against ω, with the method spread.',
              'Theory 8π²NR²/(λc) = 1.74 ± 0.14: 0.7 σ apart.'),
             ('3.29 ± 0.06 rad', '0.52 fringe',
              'Largest shift, at ω = 1.81 rad/s.',
              'The beams’ travel times differ by only 1.1 fs.'),
             ('(2.8 ± 0.2) × 10⁸ m/s', 'speed of light',
              'Inferred c: 0.7 σ from the defined value.',
              'Limited by the turn count N = 13 ± 1.')]
    y = 215
    for big, unit, l1, l2 in items:
        fig.text(144 / 1920, 1 - (y + 45) / 1080, big, fontsize=29, weight='bold', color=BLUE)
        fig.text(144 / 1920, 1 - (y + 88) / 1080, unit, fontsize=16, color=GREY)
        fig.text(900 / 1920, 1 - (y + 30) / 1080, l1, fontsize=19, color=INK)
        fig.text(900 / 1920, 1 - (y + 68) / 1080, l2, fontsize=17, color=GREY)
        y += 185
    ax = px_axes(fig, 144, 790, 1680, 2)
    ax.axhline(0, color='#c9d1d8', lw=1.5); ax.axis('off')
    fig.text(144 / 1920, 1 - 850 / 1080, 'NEXT', fontsize=14, color=GREY, weight='bold')
    fig.text(144 / 1920, 1 - 895 / 1080, '•  Interleave rest photos and reverse the rotation, to separate a drift from '
             'the Sagnac phase', fontsize=17)
    fig.text(144 / 1920, 1 - 935 / 1080, '•  Count the fibre turns (now 13 ± 1)        •  Measure the off-centre loop '
             '(manual task 4)', fontsize=17)
    return save(fig, '12_conclusions')


# ---------------------------------------------------------------------------------------
def backup_tag(fig):
    fig.text(1 - 96 / 1920, .878, 'BACKUP', fontsize=13, color=GREY, weight='bold', ha='right',
             bbox=dict(boxstyle='round,pad=0.35', fc=LIGHT, ec='none'))


def b1_quadratic():
    fig = canvas('A quadratic calibration would raise the slope by 1.7 %',
                 'Calibration residuals for the two models, and what each does to the result')
    backup_tag(fig)
    tab = calibration_table()
    qc = CAL['fits_with_eps']['2']['coef']
    lin = C0 + C1 * tab.S
    quad = qc[0] + qc[1] * tab.S + qc[2] * tab.S ** 2
    ax = px_axes(fig, 175, 205, 800, 560)
    ax.axhline(0, color=INK, lw=1)
    off = .7
    ax.errorbar(tab.S - off, (tab.omega - lin) / tab.omega * 100, tab.u / tab.omega * 100, fmt='o', color=ORANGE,
                ms=9, elinewidth=2, capsize=0, label=f'Linear (used): χ² = {CAL["fits_with_eps"]["1"]["chi2"]:.0f} for 7 dof')
    ax.errorbar(tab.S + off, (tab.omega - quad) / tab.omega * 100, tab.u / tab.omega * 100, fmt='s', color=PURPLE,
                ms=9, elinewidth=2, capsize=0, label=f'Quadratic: χ² = {CAL["fits_with_eps"]["2"]["chi2"]:.0f} for 6 dof')
    ax.set(xlim=(0, 95), ylim=(-3.3, 2.3), xlabel='Motor setting S', ylabel='Measured ω − model (% of ω)')
    style(ax)
    ax.legend(loc='lower left')
    # table of consequences
    rows = [('', 'Linear (used)', 'Quadratic'),
            ('Slope a (rad per rad/s)', f'{RES["slope"]:.3f}', f'{RES_Q["slope"]:.3f}'),
            ('a / a_th', f'{RES["ratio"]:.3f}', f'{RES_Q["ratio"]:.3f}'),
            ('c (10⁸ m/s)', f'{RES["speed_of_light"]["c"] / 1e8:.3f}', f'{RES_Q["speed_of_light"]["c"] / 1e8:.3f}'),
            ('ω at S = 90 (rad/s)', f'{C0 + C1 * 90:.3f}', f'{qc[0] + qc[1] * 90 + qc[2] * 8100:.3f}')]
    x0, y0 = 1075, 260
    for i, r in enumerate(rows):
        y = y0 + i * 62
        for j, (txt, xx) in enumerate(zip(r, [x0, 1540, 1755])):
            fig.text(xx / 1920, 1 - y / 1080, txt, fontsize=17 if i else 15,
                     color=(GREY if i == 0 else INK) if j == 0 else (ORANGE if j == 1 else PURPLE),
                     weight='bold' if i == 0 or j > 0 else 'normal', ha='left' if j == 0 else 'center')
        if i == 0:
            a = px_axes(fig, x0, y + 18, 750, 2); a.axhline(0, color='#c9d1d8'); a.axis('off')
    for y, t in [(640, 'The line is the model the manual asks for,'), (672, 'but the table turns slower than linear'),
                 (704, 'at high settings. The choice moves a and c'), (736, 'by 1.7 %, inside the 3.3 % uncertainty'),
                 (768, 'of the measurement.')]:
        fig.text(x0 / 1920, 1 - y / 1080, t, fontsize=16, color=INK)
    footer(fig, 'Both fits include the 0.48 % speed reproducibility, which is set by χ² = dof for the quadratic (so its χ² '
                'is not an independent test). A cubic term is not needed\n(χ² = 2.6 for 5 dof). The model-choice shift is '
                'reported separately and is not added to the quoted uncertainty.')
    return save(fig, 'B1_quadratic_calibration')


def b2_error_budget():
    fig = canvas('Every error source, and whether it is in the result',
                 'Manual task 5. Sizes refer to the slope a; the same relative sizes apply to c')
    backup_tag(fig)
    rows = [('Source', 'Kind', 'Size', 'Included'),
            ('Photo-to-photo scatter', 'statistical', 'SEM 0.004–0.028 rad', 'yes, fit'),
            ('Table speed at a setting', 'statistical (calibration)', '0.48 % of ω', 'yes, fit'),
            ('Wander between settings', 'statistical (fit)', 'σx = 0.041 rad', 'yes, fit'),
            ('Calibration line', 'calibration, correlated', '0.2 % of a', 'yes'),
            ('Phase-reading method', 'systematic', '3.0 % of a', 'yes'),
            ('Turn count N = 13 ± 1', 'given', '7.7 %', 'yes (theory, c)'),
            ('Loop radius R = 57.0 ± 0.5 cm', 'given', '1.8 %', 'yes (theory, c)'),
            ('Wavelength 640 nm', 'not given', '0.16 % per nm', 'no'),
            ('Steady drift in time', 'systematic', 'not separable from a', 'no'),
            ('Linear vs quadratic calibration', 'model choice', '1.7 % of a', 'no, reported')]
    xs = [150, 690, 1110, 1540]
    y0 = 215
    for i, r in enumerate(rows):
        y = y0 + i * 66
        if i % 2 == 1:
            a = px_axes(fig, 130, y - 44, 1670, 62); a.axis('off')
            a.add_patch(Rectangle((0, 0), 1, 1, transform=a.transAxes, fc=LIGHT, ec='none'))
        for j, (txt, x) in enumerate(zip(r, xs)):
            col = GREY if i == 0 else INK
            if j == 3 and i:
                col = '#2e7d4f' if txt.startswith('yes') else '#b5452b'
            fig.text(x / 1920, 1 - y / 1080, txt, fontsize=14 if i == 0 else 16.5, color=col,
                     weight='bold' if i == 0 or j == 3 else 'normal')
    footer(fig, 'Statistical parts enter the fit of Δφ against ω; the others are added in quadrature afterwards: '
                'a = 1.84 ± 0.06, c = (2.8 ± 0.2) × 10⁸ m/s.\nNot checked: the sign of the shift (rotation direction '
                'not recorded). ω comes from the calibration, not measured during the photos.')
    return save(fig, 'B2_error_budget')


def b3_random_pixels():
    fig = canvas('Ten fixed pixels give a slope that depends on the draw',
                 'Why the random-points method is left out of the method spread')
    backup_tag(fig)
    s1 = np.array([d['slope'] for d in METHODS['points_scan']])
    s2 = np.array([d['slope'] for d in METHODS['points_frame_geometry_scan']])
    ax = px_axes(fig, 175, 210, 1000, 600)
    bins = np.arange(1.35, 2.2, .05)
    ax.hist(s1, bins=bins, color=GREY, alpha=.75, label='fixed positions (31 draws)', edgecolor='white')
    ax.hist(s2, bins=bins, histtype='step', color=PURPLE, lw=2.5,
            label='same pixels, each photo’s own centre and P')
    ax.axvline(RES['slope'], color=BLUE, lw=3, label=f'ring minima: {RES["slope"]:.2f}')
    ax.axvline(RES['theory'], color=THEORY, lw=2.2, ls=(0, (6, 3)), label=f'prediction: {RES["theory"]:.2f}')
    ax.set(xlabel='Fitted slope a (rad per rad/s)', ylabel='Number of random draws', ylim=(0, 13.5))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    style(ax)
    ax.legend(loc='upper left', fontsize=14)
    x = 1260 / 1920
    lo, med, hi = np.percentile(s1, [16, 50, 84])
    for y, t, size, col in [(250, 'SPREAD OVER DRAWS', 14, GREY), (297, f'median {med:.2f}', 19, INK),
                            (338, f'16–84 %: {lo:.2f} – {hi:.2f}', 19, INK),
                            (379, f'full range {s1.min():.2f} – {s1.max():.2f}', 19, INK),
                            (470, 'WHY', 14, GREY),
                            (510, 'A fixed pixel reads the local', 16, INK),
                            (544, 'phase at its own radius, so a', 16, INK),
                            (578, 'change in ring spacing or a', 16, INK),
                            (612, 'moving centre looks like a', 16, INK),
                            (646, 'Sagnac shift. Ten pixels see', 16, INK),
                            (680, 'only half a fringe of motion.', 16, INK)]:
        fig.text(x, 1 - y / 1080, t, fontsize=size, color=col, weight='bold' if size == 14 else 'normal')
    footer(fig, 'Each draw: 10 pixels in the fringe sector, read at the same image position in all 161 photos; one '
                'alternating least-squares fit gives every photo’s phase.\nThe draw shown in the method comparison '
                f'(seed 20261011) gives {RES["methods"]["random points"]["slope"]:.2f}.')
    return save(fig, 'B3_random_pixels')


def b4_data_table():
    fig = canvas('Measured values at each setting', 'Ring-minima phase; ω from the linear calibration')
    backup_tag(fig)
    sp = SPEEDS
    first = FRAMES.groupby('S').time_s.min()
    head = ('Motor setting S', 'Photos', 'Clock', 'ω (rad/s)', 'Δφ (rad)', 'Fringe period P (10⁵ px²)')
    xs = [150, 450, 640, 860, 1150, 1440]
    for j, (t, x) in enumerate(zip(head, xs)):
        fig.text(x / 1920, 1 - 215 / 1080, t, fontsize=14, color=GREY, weight='bold')
    for i, r in sp.iterrows():
        y = 270 + i * 62
        if i % 2 == 0:
            a = px_axes(fig, 130, y - 42, 1670, 60); a.axis('off')
            a.add_patch(Rectangle((0, 0), 1, 1, transform=a.transAxes, fc=LIGHT, ec='none'))
        t = first[r.S]
        u_om = np.hypot(r.u_omega_curve, r.u_omega_repro)
        vals = ['at rest' if r.S == 0 else f'{r.S:g}', f'{int(r.n)}',
                f'{int(t // 3600):02d}:{int(t % 3600 // 60):02d}:{int(t % 60):02d}',
                '0' if r.S == 0 else f'{r.omega:.3f} ± {u_om:.3f}', f'{r.dphi:.2f} ± {r.u_fit:.2f}',
                f'{r.P / 1e5:.3f}']
        for x, v in zip(xs, vals):
            fig.text(x / 1920, 1 - y / 1080, v, fontsize=17, color=INK)
    footer(fig, 'Δφ errors are the per-setting fit weights; they share the rest reference, so they are not independent '
                'errors of the differences.\nω errors combine the calibration line and the 0.48 % speed reproducibility. '
                'Settings were recorded in this order on 11 October 2026.')
    return save(fig, 'B4_data_table')


# ---------------------------------------------------------------------------------------
REGISTRY = {}
for _name, _fn in list(globals().items()):
    if callable(_fn) and (_name[:1] in 'sb') and _name[1:2].isdigit():
        REGISTRY[_name] = _fn

if __name__ == '__main__':
    names = sys.argv[1:] or sorted(REGISTRY)
    for n in names:
        REGISTRY[n]()
