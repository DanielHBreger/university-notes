"""
Figures for the Sagnac report, in the report's style (report_style.py): drawn at their printed
size, column width unless noted (figure* = full text width), vector PDF plus a 600 dpi PNG, no
titles and no fit results written into the image (those go in the caption and the text).

Every number behind these figures comes from the outputs of calibration.py and sagnac.py, so a
value read off a figure is the value in sagnac.txt and RESULTS.md. The script only draws.

Reads:
    calibration.json                                the timed settings and the line
    sagnac.json, sagnac_speeds.csv                  the regression, the methods, the ten setting means
    fringes_frames.csv                              ring centre, radii and phase of every photo (fringes.py)
    presentation_slides/source/random_points_scan.json   the random-pixel method, 1000 draws (optional figure)
    Measurements/<S>/*.png                          per setting, the photo whose phase is the setting's median

Writes to report_figures/ (name: what it shows; width):
    setup           top view of the apparatus, after Fig. 1 of the manual; column
    calibration     omega of the timed settings against the motor setting S with the linear fit and its
                    prediction band; column
    conversion      the nine photo settings converted to omega by that line, with their uncertainties; column
    patterns        one photograph per setting with the fitted dark rings and centre; figure*
    rings           intensity around the ring centre against r and against r^2, at rest and at
                    the top speed: equal spacing in r^2; column
    extraction      the ring finder on the top-speed photo (dark-ring minima, concentric circles)
                    and R_m^2/P against the ring number m at rest and at the top speed; column
    phase_vs_omega  setting means of the phase against omega, the regression line and the
                    prediction; column
    methods         the three phase-reading methods: difference from the ring minima per
                    setting, and the slope each one gives; column
    drift           omega and each photo's phase minus the fitted line against time; column
    random_pixels   (optional) slopes of 1000 draws of ten fixed pixels, and the spread against
                    the number of pixels per draw; column

Usage, from this folder (about a minute; LaTeX is used for the text when it is installed):
    python report_figures.py              all figures
    python report_figures.py rings drift  only these
"""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import ndimage

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Rectangle  # noqa: E402

import fringes  # noqa: E402
import report_style as rs  # noqa: E402
from calibration import omega_of, u_omega_pred  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / 'report_figures'

rs.apply()
TEX = plt.rcParams['text.usetex']
BLUE, VERM, GREEN, PINK, ORANGE, SKY = rs.COLORS[:6]      # Okabe-Ito
GREY = '0.45'
RAD_S = r'rad\,s$^{-1}$' if TEX else r'rad s$^{-1}$'      # the unit rad/s in axis labels
SLOPE_UNIT = rf'rad/({RAD_S})'                            # the unit of the Sagnac slope

CAL = json.loads((HERE / 'calibration.json').read_text())
RES = json.loads((HERE / 'sagnac.json').read_text())
SETTINGS = pd.read_csv(HERE / 'sagnac_speeds.csv')        # ten setting means, rest first
PHOTOS = pd.read_csv(HERE / 'fringes_frames.csv')         # one row per photo
PHOTOS['omega'] = np.where(PHOTOS.S > 0, omega_of(PHOTOS.S, CAL['line']), 0.0)
PHOTOS['minutes'] = (PHOTOS.time_s - PHOTOS.time_s.min()) / 60
REST_PHI = PHOTOS.loc[PHOTOS.S == 0, 'phi_rad'].mean()
A_TH, U_A_TH = RES['theory']['a_th'], RES['theory']['u_a_th']    # a_th = 8 pi^2 N R^2 / (lambda c)
SCAN = HERE / 'presentation_slides' / 'source' / 'random_points_scan.json'


# ---- helpers ------------------------------------------------------------------------------
def panel_letter(ax, text, x=0.03, y=0.96, ha='left', va='top'):
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va)


def omega_text(omega):
    return r'$\omega = 0$' if omega == 0 else rf'$\omega = {omega:.2f}$ {RAD_S}'


def legend_in_order(ax, labels, **kw):
    """A legend whose entries follow `labels`, whatever order they were drawn in."""
    handles, drawn = ax.get_legend_handles_labels()
    ax.legend([handles[drawn.index(label)] for label in labels], labels, **kw)


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    rs.save(fig, OUT / name)
    plt.close(fig)
    print('wrote', (OUT / name).with_suffix('.pdf').relative_to(HERE), flush=True)


def representative(S):
    """The photo of setting S whose phase is closest to the setting's median phase."""
    group = PHOTOS[PHOTOS.S == S]
    return group.iloc[int(np.argmin(np.abs(group.phi_rad - group.phi_rad.median())))]


def ring_radii(row):
    """The fitted dark-ring radii of a photo (two at rest, three from S = 44)."""
    return np.array([row[f'R{k}'] for k in (1, 2, 3) if np.isfinite(row[f'R{k}'])])


def crop(row, half):
    """Red channel of a photo, a square of side 2 half around its fitted ring centre, zero-padded
    past the image edge."""
    image = fringes.load_red(HERE / row.file)
    x0, y0 = int(round(row.cx)) - half, int(round(row.cy)) - half
    square = np.zeros((2 * half, 2 * half))
    xa, ya = max(x0, 0), max(y0, 0)
    xb, yb = min(x0 + 2 * half, image.shape[1]), min(y0 + 2 * half, image.shape[0])
    square[ya - y0:yb - y0, xa - x0:xb - x0] = image[ya:yb, xa:xb]
    return square


def profile_against_r(row, r_max=560):
    """Mean of the smoothed red channel in 1 px annuli about the photo's fitted centre, fringe
    sector only. Returns the radii (bin centres) and the profile."""
    image = ndimage.gaussian_filter(fringes.load_red(HERE / row.file), fringes.SMOOTH_PX)
    edges = np.arange(0, r_max + 2, 1.0) ** 2
    return np.arange(r_max + 1) + 0.5, fringes.sector_profile(image, row.cx, row.cy, edges)


def ring_line(row):
    """The line R_m^2 = P (m + eta) of one photo, as fringes.ring_phase fits it: ring numbers m,
    squared radii, P and eta."""
    s = np.sort(ring_radii(row)) ** 2
    m = fringes.ring_numbers(s)
    P, s0 = np.polyfit(m, s, 1)
    return m, s, P, s0 / P


# ---- 1. setup -------------------------------------------------------------------------------
def setup():
    """Top view of the table, after Fig. 1 of the manual: laser, polarizer P1, beam splitter, the two
    coupling lenses L1 and L2 into the two ends of the fiber loop, polarizer P2, lens L3 and the camera.
    The beam splitter is drawn with its diagonal as in the manual, so that each beam is reflected once
    and transmitted once on its way to the camera. The small off-centre loop (task 4) is left out."""
    ink = '0.15'
    fig, ax = plt.subplots(figsize=rs.figsize(aspect=2.40 / 2.60), layout='constrained')   # the aspect of the limits below
    ax.set_aspect('equal')
    ax.axis('off')
    ax.add_patch(Circle((0, 0), 1.0, fc='0.95', ec='0.70', lw=0.6, zorder=0))       # the turntable
    for r in np.linspace(.925, .975, 5):                                            # fiber loop, a bundle of turns
        ax.add_patch(Circle((0, 0), r, fill=False, ec='0.30', lw=0.35, zorder=1))
    beam = dict(color=VERM, lw=1.1, zorder=2, solid_capstyle='butt')
    ax.plot([0, 0], [-.72, .80], **beam)                    # laser -> P1 -> splitter -> L2 -> fiber end
    ax.plot([-.80, 0], [0, 0], **beam)                      # splitter -> L1 -> fiber end
    ax.plot([0, .64], [0, 0], color=VERM, lw=1.1, ls=(0, (1.6, 1.2)), zorder=2)    # both beams back -> camera

    def lens(x, y, vertical, text, dx=0, dy=0):
        w, h = (.035, .14) if vertical else (.14, .035)
        ax.add_patch(Ellipse((x, y), w, h, fc='white', ec=ink, lw=0.6, zorder=4))
        ax.text(x + dx, y + dy, text, ha='center', va='center', fontsize=8)

    def polarizer(x, y, vertical, text, dx=0, dy=0):
        w, h = (.03, .15) if vertical else (.15, .03)
        ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, fc='0.75', ec=ink, lw=0.5, zorder=4))
        ax.text(x + dx, y + dy, text, ha='center', va='center', fontsize=8)

    half_cube = .085                                           # beam-splitter cube, diagonal as in the manual
    ax.add_patch(Rectangle((-half_cube, -half_cube), 2 * half_cube, 2 * half_cube, fc='white', ec=ink, lw=0.7, zorder=5))
    ax.plot([-half_cube, half_cube], [half_cube, -half_cube], color=ink, lw=0.6, zorder=6)
    ax.text(-.11, -.13, 'beam\nsplitter', ha='right', va='top', fontsize=8, linespacing=1.0)
    ax.add_patch(FancyBboxPatch((-.045, -.88), .09, .16, boxstyle='round,pad=0.005', fc=ORANGE, ec=ink,
                                lw=0.5, zorder=4))
    ax.text(.08, -.80, 'laser', ha='left', va='center', fontsize=8)
    polarizer(0, -.42, False, r'$P_1$', dx=.15)
    lens(-.60, 0, True, r'$L_1$', dy=-.14)
    lens(0, .60, False, r'$L_2$', dx=.15)
    for x, y in [(-.86, 0), (0, .86)]:                         # fiber holders: the two ends of the loop
        ax.add_patch(FancyBboxPatch((x - .045, y - .045), .09, .09, boxstyle='round,pad=0.005', fc='0.82',
                                    ec=ink, lw=0.5, zorder=4))
    polarizer(.34, 0, True, r'$P_2$', dy=-.15)
    lens(.50, 0, True, r'$L_3$', dy=-.14)
    ax.add_patch(FancyBboxPatch((.64, -.06), .15, .12, boxstyle='round,pad=0.008', fc='0.25', ec=ink,
                                lw=0.5, zorder=4))
    ax.text(.715, .12, 'camera', ha='center', va='bottom', fontsize=8)
    arrow_head = '-|>,head_length=0.45,head_width=0.22'
    # the two directions around the loop: one arc each, just outside the loop, labelled beside it
    for start_deg, end_deg, color, text, tx, ty in [(100, 140, BLUE, 'beam 1', -.52, 1.17),
                                                    (200, 160, GREEN, 'beam 2', -1.17, -.47)]:
        angles = np.radians(np.linspace(start_deg, end_deg, 60))
        ax.plot(1.08 * np.cos(angles), 1.08 * np.sin(angles), color=color, lw=1.0)
        ax.annotate('', xy=(1.08 * np.cos(angles[-1]), 1.08 * np.sin(angles[-1])),
                    xytext=(1.08 * np.cos(angles[-4]), 1.08 * np.sin(angles[-4])),
                    arrowprops=dict(arrowstyle=arrow_head, color=color, lw=1.0))
        ax.text(tx, ty, text, color=color, ha='center', va='center', fontsize=8)
    angles = np.radians(np.linspace(18, 62, 60))                # sense of rotation
    ax.plot(1.10 * np.cos(angles), 1.10 * np.sin(angles), color=ink, lw=0.9)
    ax.annotate('', xy=(1.10 * np.cos(angles[-1]), 1.10 * np.sin(angles[-1])),
                xytext=(1.10 * np.cos(angles[-4]), 1.10 * np.sin(angles[-4])),
                arrowprops=dict(arrowstyle=arrow_head, color=ink, lw=0.9))
    ax.text(.86, 1.00, r'$\omega$', ha='right', va='center', fontsize=9)
    angle = np.radians(-52)                                     # loop radius, drawn in the empty lower right
    ax.annotate('', xy=(.95 * np.cos(angle), .95 * np.sin(angle)), xytext=(0.13 * np.cos(angle), 0.13 * np.sin(angle)),
                arrowprops=dict(arrowstyle='-|>,head_length=0.4,head_width=0.18', color='0.35', lw=0.6))
    ax.text(.62 * np.cos(angle) - .09, .62 * np.sin(angle) - .03, r'$R$', ha='center', va='center', fontsize=8,
            color='0.25')
    ax.annotate(r'fiber loop, $N$ turns', xy=(.95 * np.cos(np.radians(222)), .95 * np.sin(np.radians(222))),
                xytext=(-1.40, -1.08), ha='left', va='center', fontsize=8,
                arrowprops=dict(arrowstyle='-', color='0.35', lw=0.5, shrinkA=1, shrinkB=0))
    ax.set_xlim(-1.45, 1.15)
    ax.set_ylim(-1.15, 1.25)
    save(fig, 'setup')


# ---- 2. calibration and conversion ----------------------------------------------------------
def line_axes():
    """A column-width axes of omega against S with the calibration line and its +-1 sigma prediction band,
    shared by the calibration and conversion figures."""
    fig, ax = plt.subplots(figsize=rs.figsize(aspect=0.72), layout='constrained')
    S = np.linspace(0, 95, 200)
    omega, u_omega = omega_of(S, CAL['line']), u_omega_pred(S, CAL['line'])
    ax.fill_between(S, omega - u_omega, omega + u_omega, color=BLUE, alpha=0.18, lw=0, label=r'Fit $\pm1\sigma$')
    ax.plot(S, omega, color=BLUE, label='Linear fit')
    ax.set_xlim(0, 95)
    ax.set_ylim(0, 2.0)
    ax.set_yticks(np.arange(0, 2.01, 0.5))
    ax.set_xlabel(r'Motor setting $S$')
    ax.set_ylabel(rf'$\omega$ ({RAD_S})')
    return fig, ax


def calibration():
    """omega = 2 pi / <T> of the timed turns against the motor setting S, with the linear fit that
    converts the photo settings and its prediction band (the line's uncertainty and the rms scatter of
    the timed settings about it, +-1 sigma). Error bars: the uncertainty of each timed setting (lap
    scatter and reaction time), mostly smaller than the markers."""
    timed = pd.DataFrame(CAL['settings'])
    fig, ax = line_axes()
    ax.errorbar(timed.S, timed.omega, timed.u_omega, fmt='o', color='k', ms=2.8, elinewidth=0.6, capsize=0, zorder=3,
                label='Timed turns')
    legend_in_order(ax, ['Timed turns', 'Linear fit', r'Fit $\pm1\sigma$'], loc='upper left')
    save(fig, 'calibration')


def conversion():
    """How the nine motor settings of the photographs become omega: dotted guides from each setting up to
    the calibration line and across to the omega axis, and at the line the prediction uncertainty (the
    line's uncertainty and the rms scatter of the timed settings about it) as a bare error bar. No
    measurement was made at these settings; the values are those of Table 2. The timed settings are in
    the calibration figure."""
    fig, ax = line_axes()
    moving = SETTINGS[SETTINGS.S > 0]
    ax.vlines(moving.S, 0, moving.omega, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=1)    # setting up to the line
    ax.hlines(moving.omega, 0, moving.S, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=1)    # across to the omega axis
    read_off = r'$\omega$ read off the line, $\pm1\sigma$'
    ax.errorbar(moving.S, moving.omega, moving.u_omega, fmt='none', ecolor=VERM, elinewidth=1.1, capsize=2.5,
                capthick=0.8, zorder=3, label=read_off)
    legend_in_order(ax, [read_off, 'Linear fit', r'Fit $\pm1\sigma$'], loc='upper left')
    save(fig, 'conversion')


# ---- 3. patterns ----------------------------------------------------------------------------
def patterns():
    """One photograph per setting (the one whose phase is the median of its setting), red channel on
    a grey scale, cropped around its fitted ring centre (+), with the fitted dark rings dashed.
    omega from the linear calibration. Full text width."""
    half = 500                                     # crop half-width (px); the centre is 505 px below the top edge
    n_cols = 5
    n_rows = int(np.ceil(len(SETTINGS) / n_cols))
    width = rs.TEXT_WIDTH
    gap = 0.03                                     # space between panels, as a fraction of a panel
    panel = width / (n_cols + (n_cols - 1) * gap)  # square panels that fill the full text width
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(width, panel * (n_rows + (n_rows - 1) * gap)),
                             gridspec_kw=dict(left=0, right=1, bottom=0, top=1, wspace=gap, hspace=gap))
    circle = np.linspace(0, 2 * np.pi, 400)
    for ax, (_, setting) in zip(axes.flat, SETTINGS.iterrows()):
        row = representative(setting.S)
        ax.imshow(crop(row, half), cmap='gray', vmin=0, vmax=225, extent=(-half, half, half, -half),
                  interpolation='antialiased')
        for R in ring_radii(row):
            ax.plot(R * np.cos(circle), R * np.sin(circle), ls=(0, (4, 3)), color=ORANGE, lw=0.6)
        ax.plot(0, 0, '+', color=ORANGE, ms=4, mew=0.6)
        ax.set_xlim(-half, half)
        ax.set_ylim(half, -half)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.minorticks_off()
        ax.text(0.04, 0.96, omega_text(setting.omega), transform=ax.transAxes, ha='left', va='top', color='white',
                fontsize=7.5, bbox=dict(boxstyle='round,pad=0.15', fc='black', ec='none', alpha=0.55))
    for ax in list(axes.flat)[len(SETTINGS):]:
        ax.axis('off')
    save(fig, 'patterns')


# ---- 4. rings -------------------------------------------------------------------------------
def rings_figure():
    """Intensity of the red channel averaged over the fringe sector in 1 px annuli around each photo's
    fitted ring centre (smoothed as in the ring finder), at rest and at the top speed, (a) against r and
    (b) against r^2. Dashed: the fitted dark-ring radii of each photo. The arrows give the top-speed
    photo's ring spacings: unequal in r, equal in r^2."""
    rest, fast = representative(0.0), representative(SETTINGS.S.max())
    fig, (ax_r, ax_r2) = plt.subplots(2, 1, figsize=rs.figsize(aspect=0.92), layout='constrained')
    fig.get_layout_engine().set(hspace=0.06, h_pad=0.03)
    for row, color, text in [(rest, BLUE, 'At rest'), (fast, VERM, omega_text(fast.omega))]:
        r, profile = profile_against_r(row, r_max=560)
        radii = np.sort(ring_radii(row))
        for ax, x, x_rings in [(ax_r, r, radii), (ax_r2, r ** 2 / 1e5, radii ** 2 / 1e5)]:
            ax.plot(x, profile, color=color, lw=0.9, label=text)
            for v in x_rings:
                ax.axvline(v, color=color, lw=0.5, ls=(0, (3, 2)))
    radii = np.sort(ring_radii(fast))
    for ax, x_rings, fmt in [(ax_r, radii, lambda d: f'{d:.0f} px'), (ax_r2, radii ** 2 / 1e5, lambda d: f'{d:.2f}')]:
        for x0, x1 in zip(x_rings[:-1], x_rings[1:]):
            ax.annotate('', xy=(x1, 222), xytext=(x0, 222),
                        arrowprops=dict(arrowstyle='<->', color=VERM, lw=0.6, shrinkA=0, shrinkB=0,
                                        mutation_scale=6))
            ax.text(0.5 * (x0 + x1), 228, fmt(x1 - x0), color=VERM, ha='center', va='bottom', fontsize=7.5)
    ax_r.set_xlim(0, 560)
    ax_r2.set_xlim(0, 560 ** 2 / 1e5)
    for ax in (ax_r, ax_r2):
        ax.set_ylim(0, 262)
        ax.set_ylabel('Intensity (counts)')
    ax_r.set_xlabel(r'Distance from the ring center $r$ (px)')
    ax_r2.set_xlabel(r'$r^2$ ($10^5$ px$^2$)')
    ax_r2.legend(loc='upper right', bbox_to_anchor=(1.0, 0.80), frameon=True, facecolor='white', edgecolor='none',
                 framealpha=1)
    panel_letter(ax_r, '(a)', x=0.97, y=0.95, ha='right')
    panel_letter(ax_r2, '(b)', x=0.97, y=0.95, ha='right')
    fig.align_ylabels([ax_r, ax_r2])
    save(fig, 'rings')


# ---- 5. extraction --------------------------------------------------------------------------
RING_COLORS = [SKY, ORANGE, GREEN]


def extraction():
    """(a) The ring finder of fringes.py on the top-speed photo: the darkest point of each dark ring on
    each ray (dots, one colour per ring; rays every 1.5 deg over the sector with fringes) and the
    concentric circles fitted to them (dashed), centre +. (b) For the rest and the top-speed photo,
    R_m^2 / P against the dark-ring number m, with the line R_m^2 = P (m + eta) fitted to each photo;
    the intercepts eta differ by the shift of the pattern, in fringes."""
    rest, fast = representative(0.0), representative(SETTINGS.S.max())
    fig, (ax_photo, ax_line) = plt.subplots(1, 2, figsize=rs.figsize(aspect=0.50), width_ratios=[1, 1.15],
                                            layout='constrained')
    fig.get_layout_engine().set(wspace=0.04, w_pad=0.03)
    half = 500
    minima = fringes.analyse_photo(HERE / fast.file)['minima']      # (ring index, x, y), innermost ring first
    ax_photo.imshow(crop(fast, half), cmap='gray', vmin=0, vmax=225, extent=(-half, half, half, -half),
                    interpolation='antialiased')
    for ring in range(int(minima[:, 0].max()) + 1):
        points = minima[minima[:, 0] == ring]
        ax_photo.plot(points[:, 1] - fast.cx, points[:, 2] - fast.cy, 'o', ms=0.9, mew=0, color=RING_COLORS[ring % 3])
    circle = np.linspace(0, 2 * np.pi, 400)
    for R in ring_radii(fast):
        ax_photo.plot(R * np.cos(circle), R * np.sin(circle), color='white', lw=0.5, ls=(0, (3, 2)))
    ax_photo.plot(0, 0, '+', color='white', ms=4, mew=0.6)
    ax_photo.set_xlim(-half, half)
    ax_photo.set_ylim(half, -half)
    ax_photo.set_xticks([])
    ax_photo.set_yticks([])
    ax_photo.minorticks_off()
    ax_photo.text(0.04, 0.96, '(a)', transform=ax_photo.transAxes, ha='left', va='top', color='white')
    etas = []
    for row, color, marker, text in [(rest, BLUE, 'o', 'At rest'), (fast, VERM, 's', omega_text(fast.omega))]:
        m, s, P, eta = ring_line(row)
        etas.append(eta)
        x = np.linspace(-0.08, 2.25, 20)
        ax_line.plot(x, x + eta, color=color, lw=0.8)
        ax_line.plot(m, s / P, marker, color=color, ms=3.2, label=text, zorder=3)
    ax_line.axvline(0, color=GREY, lw=0.5, zorder=0)
    ax_line.annotate('', xy=(-0.22, etas[1]), xytext=(-0.22, etas[0]),
                     arrowprops=dict(arrowstyle='<->', color='k', lw=0.6, shrinkA=0, shrinkB=0, mutation_scale=6))
    ax_line.text(-0.28, 0.5 * (etas[0] + etas[1]), rf'${etas[0] - etas[1]:.2f}$', ha='right', va='center', fontsize=8)
    ax_line.set_xlim(-0.85, 2.25)
    ax_line.set_ylim(-0.2, 3.0)
    ax_line.set_xticks([0, 1, 2])
    ax_line.set_xlabel(r'Dark ring $m$')
    ax_line.set_ylabel(r'$R_m^2/P$')
    ax_line.legend(loc='upper left', handletextpad=0.2)
    panel_letter(ax_line, '(b)', x=0.96, y=0.04, ha='right', va='bottom')
    save(fig, 'extraction')


# ---- 6. phase against omega -----------------------------------------------------------------
def phase_vs_omega():
    """Mean phase change of each setting's photos (ring minima) against omega, the regression line
    Delta phi = a omega + b fitted to the ten points with equal weights (sagnac.py), and the prediction
    a_th omega with its 1 sigma band (N, R). Error bars: the calibration line's uncertainty in omega and
    the standard error of the mean of the setting's photos in Delta phi, both mostly smaller than the
    markers."""
    a, b = RES['slope']['a'], RES['slope']['b']
    fig, ax = plt.subplots(figsize=rs.figsize(aspect=0.72), layout='constrained')
    omega = np.linspace(0, 1.95, 100)
    prediction = r'Prediction $a_\mathrm{th}\,\omega$, $\pm1\sigma$'
    ax.fill_between(omega, (A_TH - U_A_TH) * omega, (A_TH + U_A_TH) * omega, color=VERM, alpha=0.15, lw=0)
    ax.plot(omega, A_TH * omega, color=VERM, ls='--', label=prediction)
    ax.plot(omega, a * omega + b, color=BLUE, label='Linear regression')
    ax.errorbar(SETTINGS.omega, SETTINGS.dphi, xerr=SETTINGS.u_omega, yerr=SETTINGS['sem'], fmt='o', color=BLUE, ms=3,
                elinewidth=0.6, capsize=0, zorder=3, label='Setting means')
    legend_in_order(ax, ['Setting means', 'Linear regression', prediction], loc='upper left')
    ax.set_xlabel(rf'$\omega$ ({RAD_S})')
    ax.set_ylabel(r'$\Delta\phi$ (rad)')
    ax.set_xlim(-0.05, 1.92)
    ax.set_ylim(-0.25, 3.6)
    save(fig, 'phase_vs_omega')


# ---- 7. methods -----------------------------------------------------------------------------
METHODS = [('Ring minima', None, BLUE, 'o', BLUE),                 # (label, column key, colour, marker, marker face)
           ('Cosine fit', 'cos', GREEN, 's', GREEN),
           ('Template matching', 'template', PINK, 'D', 'white')]


def methods():
    """(a) Phase change of each setting by the cosine fit of the r^2 profile and by template matching to
    the rest profile, minus the ring-minima value; each method's phase is relative to its own rest mean.
    (b) The slope of each method's ten points by the same regression (sagnac.json), error bars from the
    regression; shaded: the prediction a_th +- 1 sigma."""
    fig, (ax_diff, ax_slopes) = plt.subplots(2, 1, figsize=rs.figsize(aspect=1.0), height_ratios=[1.35, 1],
                                             layout='constrained')
    fig.get_layout_engine().set(hspace=0.08, h_pad=0.04)
    ax_diff.axhline(0, color=BLUE, lw=0.9, label='Ring minima')
    for text, key, color, marker, face in METHODS[1:]:
        difference = SETTINGS[f'dphi_{key}'] - SETTINGS.dphi
        ax_diff.plot(SETTINGS.omega, difference, color=color, lw=0.7, marker=marker, ms=3, mfc=face, mec=color, mew=0.7,
                     label=text)
    ax_diff.set_xlim(-0.05, 1.9)
    ax_diff.set_ylim(-0.21, 0.27)
    ax_diff.set_xlabel(rf'$\omega$ ({RAD_S})')
    ax_diff.set_ylabel(r'Difference from minima (rad)')
    ax_diff.legend(loc='upper left', ncol=2, columnspacing=1.0)
    panel_letter(ax_diff, '(a)', x=0.97, y=0.04, ha='right', va='bottom')
    ax_slopes.axvspan(A_TH - U_A_TH, A_TH + U_A_TH, color=VERM, alpha=0.15, lw=0)
    ax_slopes.axvline(A_TH, color=VERM, ls='--', lw=0.9)
    n = len(METHODS)
    slopes = [(RES['slope']['a'], RES['slope']['u_stat'])] + \
             [(RES['methods'][k]['a'], RES['methods'][k]['u_a']) for k in ('cosine fit', 'template matching')]
    for j, ((text, key, color, marker, face), (value, u)) in enumerate(zip(METHODS, slopes)):
        ax_slopes.errorbar(value, n - 1 - j, xerr=u, fmt=marker, color=color, mfc=face, mec=color, ms=3.4, mew=0.8,
                           elinewidth=0.8, capsize=0)
    ax_slopes.text(A_TH + 0.006, -0.62, r'$a_\mathrm{th}$', color=VERM, ha='left', va='center', fontsize=8)
    ax_slopes.set_yticks(range(n))
    ax_slopes.set_yticklabels([method[0] for method in METHODS[::-1]])
    ax_slopes.tick_params(axis='y', which='both', length=0)
    ax_slopes.set_ylim(-0.95, n - 0.4)
    ax_slopes.set_xlim(1.55, 1.95)
    ax_slopes.set_xlabel(rf'Slope $a$ ({SLOPE_UNIT})')
    panel_letter(ax_slopes, '(b)', x=0.97, y=0.93, ha='right')
    save(fig, 'methods')


# ---- 8. drift -------------------------------------------------------------------------------
def drift():
    """(a) omega of each photo (linear calibration) against the time of the photo: the settings were
    taken in increasing order, so omega and time rise together. (b) Each photo's phase change minus the
    fitted line a omega + b, with a straight line fitted to each setting's photos against time."""
    a, b = RES['slope']['a'], RES['slope']['b']
    photos = PHOTOS.sort_values('time_s')
    photos = photos.assign(residual=(photos.phi_rad - REST_PHI) - (a * photos.omega + b))
    fig, (ax_omega, ax_resid) = plt.subplots(2, 1, sharex=True, figsize=rs.figsize(aspect=0.90), height_ratios=[1, 1.6],
                                             layout='constrained')
    fig.get_layout_engine().set(hspace=0.02, h_pad=0.02)
    ax_omega.step(photos.minutes, photos.omega, where='post', color='k', lw=0.7)
    ax_omega.plot(photos.minutes, photos.omega, 'o', color='k', ms=1.4, mew=0)
    ax_omega.set_ylabel(rf'$\omega$ ({RAD_S})')
    ax_omega.set_ylim(-0.1, 2.0)
    panel_letter(ax_omega, '(a)', x=0.97, y=0.08, ha='right', va='bottom')
    ax_resid.axhline(0, color=GREY, lw=0.5, zorder=0)
    ax_resid.plot(photos.minutes, photos.residual, 'o', color=BLUE, ms=1.8, mew=0, alpha=0.8)
    for _, group in photos.groupby('S'):
        trend = np.polyfit(group.minutes, group.residual, 1)
        ends = np.array([group.minutes.min(), group.minutes.max()])
        ax_resid.plot(ends, np.polyval(trend, ends), color='k', lw=0.9)
    ax_resid.set_xlim(-0.3, 12.4)
    ax_resid.set_ylim(-0.3, 0.3)
    ax_resid.set_xlabel('Time since the first photo (min)')
    ax_resid.set_ylabel(r'$\Delta\phi - \mathrm{line}$ (rad)')
    panel_letter(ax_resid, '(b)', x=0.97, y=0.95, ha='right')
    fig.align_ylabels([ax_omega, ax_resid])
    save(fig, 'drift')


# ---- 9. random pixels (optional) ------------------------------------------------------------
def random_pixels():
    """(a) The slope from ten fixed random pixels (the method of sagnac_update/methods.py, every photo's
    phase from one alternating least-squares fit), for 1000 draws: filled, the method as it is (fixed image
    positions); outline, the same pixels with each photo's own fitted centre and P. Lines: the ring-minima
    slope and the prediction. Failed draws (negative slope) lie off the plot. (b) Spread (SD) of the slope
    between draws against the number of pixels per draw, and the median fit error of a single draw
    (sagnac.json, random_pixels)."""
    draws = [d for d in json.loads(SCAN.read_text())['draws'] if d['n'] == 10 and d['seed'] != 20261011]
    fixed = np.array([d['fixed']['slope'] for d in draws])
    own_frame = np.array([d['frame']['slope'] for d in draws])
    fig, (ax_hist, ax_spread) = plt.subplots(2, 1, figsize=rs.figsize(aspect=1.05), height_ratios=[1.6, 1],
                                             layout='constrained')
    fig.get_layout_engine().set(hspace=0.08, h_pad=0.04)
    bins = np.arange(0.8, 2.72 + 1e-9, 0.04)
    counts_fixed, _, _ = ax_hist.hist(fixed, bins=bins, color='0.70', edgecolor='white', lw=0.3,
                                      label='Fixed pixel positions')
    counts_own, _, _ = ax_hist.hist(own_frame, bins=bins, histtype='step', color=PINK, lw=0.9,
                                    label="Each photo's own center and $P$")
    top = max(counts_fixed.max(), counts_own.max())
    ax_hist.set_ylim(0, 1.75 * top)                          # room for the legend above the bars
    line_height = 1.08 / 1.75                                # the two lines stop just above the bars
    ax_hist.axvline(RES['slope']['a'], ymax=line_height, color=BLUE, lw=1.0, label='Ring minima')
    ax_hist.axvline(A_TH, ymax=line_height, color=VERM, ls='--', lw=0.9, label=r'Prediction $a_\mathrm{th}$')
    ax_hist.set_xlim(0.8, 2.72)
    ax_hist.set_xlabel(rf'Slope $a$ ({SLOPE_UNIT})')
    ax_hist.set_ylabel('Draws')
    ax_hist.legend(loc='upper right')
    panel_letter(ax_hist, '(a)', x=0.03, y=0.95)
    by_n = {int(k): v for k, v in RES['random_pixels'].items()}
    pixels = np.array(sorted(by_n))
    ax_spread.plot(pixels, [by_n[k]['sd'] for k in pixels], 'o-', color='k', ms=3, lw=0.7,
                   label='Spread between draws (SD)')
    ax_spread.plot(pixels, [by_n[k]['u_fit'] for k in pixels], 's--', color=GREY, mfc='white', ms=3, lw=0.7,
                   label='Fit error of one draw')
    ax_spread.set_xscale('log')
    ax_spread.set_xlim(7, 115)
    ax_spread.set_xticks(pixels)
    ax_spread.set_xticklabels([str(k) for k in pixels])
    ax_spread.minorticks_off()
    ax_spread.set_ylim(0, 0.25)
    ax_spread.set_xlabel('Pixels per draw')
    ax_spread.set_ylabel(rf'Uncertainty ({SLOPE_UNIT})')
    ax_spread.legend(loc='upper right')
    panel_letter(ax_spread, '(b)', x=0.03, y=0.95)
    save(fig, 'random_pixels')


FIGURES = {'setup': setup, 'calibration': calibration, 'conversion': conversion, 'patterns': patterns,
           'rings': rings_figure, 'extraction': extraction, 'phase_vs_omega': phase_vs_omega, 'methods': methods,
           'drift': drift, 'random_pixels': random_pixels}

if __name__ == '__main__':
    names = sys.argv[1:] or list(FIGURES)
    unknown = [name for name in names if name not in FIGURES]
    if unknown:
        sys.exit(f'unknown figure(s) {unknown}; choose from {list(FIGURES)}')
    for name in names:
        FIGURES[name]()
