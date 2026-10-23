"""Figures of the interference patterns (manual task 2) and of the fringe shift.

  figures/patterns.pdf/.png   one photograph per speed setting (the frame whose phase is the
                              median of its setting), red channel, with the fitted dark rings;
                              full text width (figure*).
  figures/profiles.pdf/.png   intensity averaged over the fringe sector against r^2 about each
                              frame's fitted centre, mean over the frames of each setting: the
                              rings are equally spaced in r^2 and move inwards with omega.

Inputs: the photographs, the frame table (fringes.py) and sagnac_speeds.csv (sagnac.py).
Run from this folder after sagnac.py:  python patterns.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors as mcolors

import report_style as rs
from fringes import R_MAX, S2_BIN, load_red, sector_profile_s2

HERE = Path(__file__).resolve().parent
HALF = 470                     # half-width of the square crop around the ring centre (px)


def representative(fr, S):
    """Frame of setting S whose (unwrapped) phase is closest to the setting's median, so the
    photograph shown for each setting is a typical one rather than the first."""
    g = fr[fr.S == S]
    ph = np.angle(np.exp(1j * (g.phi_rad - np.angle(np.mean(np.exp(1j * g.phi_rad))))))
    return g.iloc[int(np.argmin(np.abs(ph - np.median(ph))))]


def unit(tex):
    """The unit rad/s for axis labels (LaTeX or matplotlib's own math text)."""
    return r'rad\,s$^{-1}$' if tex else r'rad s$^{-1}$'


def gallery(fr, sp, path):
    """figures/patterns: a grid of 5 columns, one panel per setting. Each panel is a square crop
    (2 HALF px wide) of the representative photograph around its fitted ring centre, with the
    fitted dark rings (dashed circles, radii R1-R3 from fringes.py) and the centre (+), labelled
    with omega from sagnac_speeds.csv."""
    rs.apply()
    tex = plt.rcParams['text.usetex']
    n = len(sp)
    ncol = 5
    nrow = int(np.ceil(n / ncol))
    w = rs.TEXT_WIDTH
    fig, axs = plt.subplots(nrow, ncol, figsize=(w, w / ncol * nrow * 1.0),
                            gridspec_kw=dict(wspace=0.04, hspace=0.04))
    th = np.linspace(0, 2 * np.pi, 400)
    for ax, (_, s) in zip(axs.flat, sp.iterrows()):
        row = representative(fr, s.S)
        img = load_red(HERE / row.file)
        cx, cy = int(round(row.cx)), int(round(row.cy))
        y0, x0 = max(cy - HALF, 0), cx - HALF
        crop = img[y0:cy + HALF, x0:cx + HALF]
        ax.imshow(crop, cmap='gray', vmin=0, vmax=220, interpolation='antialiased',
                  extent=(x0, x0 + crop.shape[1], y0 + crop.shape[0], y0))
        for k in (1, 2, 3):
            R = row.get(f'R{k}')
            if pd.notna(R) and R != '':
                R = float(R)
                ax.plot(row.cx + R * np.cos(th), row.cy + R * np.sin(th), '--',
                        color=rs.COLORS[4], lw=0.5)
        ax.plot(row.cx, row.cy, '+', color=rs.COLORS[4], ms=4, mew=0.6)
        ax.set_xlim(cx - HALF, cx + HALF); ax.set_ylim(cy + HALF, cy - HALF)
        ax.set_xticks([]); ax.set_yticks([])
        ax.minorticks_off()
        lab = (rf'$\omega = {s.omega:.2f}$ {unit(tex)}' if s.S > 0 else r'$\omega = 0$')
        ax.text(0.04, 0.96, lab, transform=ax.transAxes, ha='left', va='top', color='white',
                fontsize=7)
    for ax in list(axs.flat)[n:]:
        ax.axis('off')
    rs.save(fig, path)
    plt.close(fig)


def profiles(fr, sp, path):
    """figures/profiles: for every setting, the sector-averaged intensity against r^2 about each
    photograph's own ring centre, averaged over the setting's photographs; colour from light
    (rest) to dark blue (fastest). Evenly spaced fringes along r^2 that slide to smaller r^2
    with omega are the picture of the Sagnac shift."""
    rs.apply()
    tex = plt.rcParams['text.usetex']
    edges = np.arange(0, (R_MAX - 60) ** 2, S2_BIN)
    mid = 0.5 * (edges[1:] + edges[:-1])
    fig, ax = plt.subplots(figsize=rs.figsize(aspect=0.62))
    wmax = sp.omega.max()
    cmap = plt.get_cmap('Blues')
    norm = mcolors.Normalize(vmin=-0.45 * wmax, vmax=wmax)       # skip the palest blues
    for _, s in sp.iterrows():
        g = fr[fr.S == s.S]
        prof = np.mean([sector_profile_s2(load_red(HERE / r.file), r.cx, r.cy, edges)
                        for r in g.itertuples()], axis=0)
        ax.plot(mid / 1e5, prof, color=cmap(norm(s.omega)), lw=0.8)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.05, boundaries=np.linspace(0, wmax, 200))
    cb.set_label(rf'$\omega$ ({unit(tex)})')
    cb.set_ticks(np.arange(0, wmax + 1e-9, 0.5))
    cb.ax.minorticks_off()
    ax.set_xlabel(r'$r^2$ ($10^5$ px$^2$)')
    ax.set_ylabel('intensity (counts)')
    ax.set_xlim(0, mid[-1] / 1e5); ax.set_ylim(0, None)
    rs.save(fig, path)
    plt.close(fig)


def main():
    fr = pd.read_csv(HERE / 'fringes_frames.csv')
    sp = pd.read_csv(HERE / 'sagnac_speeds.csv')
    gallery(fr, sp, HERE / 'figures' / 'patterns')
    profiles(fr, sp, HERE / 'figures' / 'profiles')


if __name__ == '__main__':
    main()
