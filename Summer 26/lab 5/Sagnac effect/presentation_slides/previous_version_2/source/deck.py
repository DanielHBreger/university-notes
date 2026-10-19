"""Shared style, data and helpers for the Sagnac slide images.

Style follows the Chua presentation figures (chaotic circuit/presentation_figures/make_figures.py):
1920 x 1080 canvas, DejaVu Sans, sentence headlines, grey subtitle and footer, teal = measurement,
orange = fit/model, grey = theory.
"""
from pathlib import Path
import json
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

_HERE = Path(__file__).resolve().parent
# Default layout: <Sagnac effect>/presentation_slides/source/deck.py; override with environment variables.
DATA = Path(os.environ.get('SAGNAC_DATA', _HERE.parents[1]))      # the 'Sagnac effect' folder
UPD = DATA / 'sagnac_update'                                      # linear-calibration results (used throughout)
OUT = Path(os.environ.get('SAGNAC_OUT', _HERE.parent))            # where the slide images are written
OUT.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, PURPLE, INK = '#147da3', '#d87530', '#8257a6', '#243444'
GREY = '#586673'
THEORY = '#7d8a95'
LIGHT = '#eef3f6'
BEAM = '#d0312d'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 17, 'mathtext.fontset': 'dejavusans',
                     'axes.labelsize': 18, 'axes.titlesize': 20, 'xtick.labelsize': 15, 'ytick.labelsize': 15,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#85919b',
                     'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
                     'legend.frameon': False, 'legend.fontsize': 15, 'savefig.facecolor': 'white',
                     'figure.facecolor': 'white'})
W, H = 40 / 3, 7.5          # inches; 1920 x 1080 at 144 dpi
DPI = 144


# ---- slide furniture -------------------------------------------------------------------
def canvas(title, subtitle=''):
    fig = plt.figure(figsize=(W, H), dpi=DPI)
    if title:
        fig.text(.075, .93, title, fontsize=27, weight='bold')
    if subtitle:
        fig.text(.075, .875, subtitle, fontsize=15, color=GREY)
    return fig


def footer(fig, text, y=.035):
    fig.text(.075, y, text, fontsize=10.5, color=GREY, linespacing=1.5)


def style(ax):
    ax.grid(alpha=.15)
    ax.set_axisbelow(True)


def px_axes(fig, x, y, w, h, **kw):
    """Axes placed in slide pixels (x, y from the top-left corner)."""
    return fig.add_axes([x / 1920, 1 - (y + h) / 1080, w / 1920, h / 1080], **kw)


def check_layout(fig, name, margin=(90, 1840), ignore_overlap=()):
    """Warn about text outside the slide margins and about overlapping text boxes."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    boxes = []
    def visit(t, kind):
        if not t.get_visible() or not t.get_text().strip():
            return
        import re
        plain = re.sub(r'\$[^$]*\$', '', t.get_text())
        if ';' in plain:
            print(f'  [{name}] SEMICOLON: {t.get_text()[:70]!r}')
        if re.search(r'±\s*0\.00', t.get_text()):
            print(f'  [{name}] ZERO-HEAVY ±: {t.get_text()[:70]!r}')
        bb = t.get_window_extent(r)
        x0, x1 = bb.x0, bb.x1
        y0, y1 = 1080 - bb.y1, 1080 - bb.y0          # top-left origin
        left = 30 if kind in ('label', 'tick') else margin[0]
        if x0 < left - 1 or x1 > margin[1] + 1 or y0 < 0 or y1 > 1078:
            print(f'  [{name}] OUT OF BOUNDS ({x0:.0f},{y0:.0f})-({x1:.0f},{y1:.0f}): {t.get_text()[:60]!r}')
        boxes.append((x0, y0, x1, y1, t.get_text()[:40], kind))
    for t in fig.texts:
        visit(t, 'fig')
    for ax in fig.axes:
        if not ax.get_visible() or not ax.axison and not ax.texts:
            pass
        for t in ax.texts:
            visit(t, 'ax')
        for t in [ax.title, ax._left_title, ax._right_title, ax.xaxis.label, ax.yaxis.label]:
            visit(t, 'label')
        if ax.axison:
            for axis in (ax.xaxis, ax.yaxis):
                lo, hi = sorted(axis.get_view_interval())
                for tick in axis.get_major_ticks():
                    if lo - 1e-9 <= tick.get_loc() <= hi + 1e-9:
                        for lab in (tick.label1, tick.label2):
                            if lab.get_visible():
                                visit(lab, 'tick')
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i], boxes[j]
            if a[5] == 'tick' and b[5] == 'tick':
                continue
            ox = min(a[2], b[2]) - max(a[0], b[0])
            oy = min(a[3], b[3]) - max(a[1], b[1])
            if ox > 2 and oy > 2 and not any(k in a[4] or k in b[4] for k in ignore_overlap):
                print(f'  [{name}] OVERLAP {a[4]!r} <-> {b[4]!r} ({ox:.0f}x{oy:.0f}px)')


def save(fig, name, check=True, **kw):
    path = OUT / f'{name}.png'
    if check:
        check_layout(fig, name, **kw)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    im = Image.open(path)
    assert im.size == (1920, 1080), im.size
    print('Wrote', path.name, flush=True)
    return path


# ---- data ----------------------------------------------------------------------------
RES = json.loads((UPD / 'sagnac.json').read_text())
RES_Q = json.loads((UPD / 'optional_quadratic' / 'sagnac.json').read_text())
CAL = json.loads((DATA / 'calibration.json').read_text())
METHODS = json.loads((UPD / 'methods.json').read_text())
SPEEDS = pd.read_csv(UPD / 'sagnac_speeds.csv')
FRAMES = pd.read_csv(DATA / 'fringes_frames.csv')
MFRAMES = pd.read_csv(UPD / 'methods_frames.csv')
PROF = np.load(DATA / 'fringes_profiles.npz')

LIN = CAL['fits_with_eps']['1']          # the linear conversion used by sagnac_update
C0, C1 = LIN['coef']
EPS = CAL['conversion']['eps']
assert abs(C1 - CAL['linear']['k']) < 1e-12


def omega_lin(S):
    S = np.asarray(S, float)
    return np.where(S > 0, C0 + C1 * S, 0.0)


# per-frame omega (linear calibration) and phase relative to the mean rest phase
FRAMES['omega'] = omega_lin(FRAMES.S)
REST_PHI = FRAMES.loc[FRAMES.S == 0, 'phi_rad'].mean()
FRAMES['dphi'] = FRAMES.phi_rad - REST_PHI
FRAMES['minutes'] = (FRAMES.time_s - FRAMES.time_s.min()) / 60
assert np.allclose(SPEEDS.omega, omega_lin(SPEEDS.S), atol=1e-5)
assert np.allclose(SPEEDS.dphi, SPEEDS.phi - SPEEDS.phi[0], atol=1e-5)

N_TURNS, R_LOOP, LAM, C_DEF = 13, .57, 640e-9, 299792458.0
AREA = np.pi * R_LOOP ** 2


def load_turns():
    """Stopwatch laps per setting from 'Rotation calibration.csv' (empty cells skipped)."""
    df = pd.read_csv(DATA / 'Rotation calibration.csv', encoding='utf-8-sig')
    S = df['Speed'].to_numpy(float)
    laps = [row.dropna().to_numpy(float) for _, row in df.drop(columns='Speed').iterrows()]
    return S, laps


def calibration_table():
    S, laps = load_turns()
    rows = []
    for s, lap in zip(S, laps):
        n = len(lap)
        T = lap.mean()
        u_T = lap.std(ddof=1) / np.sqrt(n * (n + 1))
        w = 2 * np.pi / T
        u_sw = w * u_T / T
        rows.append(dict(S=s, n=n, T=T, omega=w, u_sw=u_sw, u=np.hypot(u_sw, EPS * w)))
    return pd.DataFrame(rows)


def photo(file):
    """RGB array of a fringe photograph (only the red channel carries light)."""
    return np.asarray(Image.open(DATA / file).convert('RGB'))


def crop(img, cx, cy, half=470):
    x0, y0 = int(round(cx)) - half, int(round(cy)) - half
    pad = np.zeros((2 * half, 2 * half, 3), np.uint8)
    xs0, ys0 = max(x0, 0), max(y0, 0)
    xs1, ys1 = min(x0 + 2 * half, img.shape[1]), min(y0 + 2 * half, img.shape[0])
    pad[ys0 - y0:ys1 - y0, xs0 - x0:xs1 - x0] = img[ys0:ys1, xs0:xs1]
    return pad, (cx - x0, cy - y0)


def representative(S):
    """The frame of setting S whose phase is closest to the setting's median phase."""
    g = FRAMES[FRAMES.S == S]
    return g.iloc[int(np.argmin(np.abs(g.phi_rad - g.phi_rad.median())))]


def rings(row):
    return [row[f'R{k}'] for k in (1, 2, 3) if np.isfinite(row[f'R{k}'])]


def fmt_pm(x, u, digits):
    return f'{x:.{digits}f} ± {u:.{digits}f}'


def paren(x, u):
    """Concise uncertainty notation: 0.020386 +- 0.000043 -> '0.02039(4)' (one significant figure)."""
    from math import floor, log10
    e = floor(log10(u))
    d = int(round(u / 10 ** e))
    if d == 10:
        d, e = 1, e + 1
    dec = max(0, -e)
    return f'{x:.{dec}f}({d})'.replace('-', '−')


THETA = (-75.0, 180.0)          # sector with fringes (deg), as in fringes.py


def sector_profile_r(file, cx, cy, r_max=560, smooth=4.0):
    """Mean of the smoothed red channel in 1 px annuli about (cx, cy), fringe sector only."""
    from scipy import ndimage
    img = ndimage.gaussian_filter(photo(file)[..., 0].astype(float), smooth)
    y, x = np.indices(img.shape)
    r = np.hypot(x - cx, y - cy)
    th = np.degrees(np.arctan2(-(y - cy), x - cx))
    sel = (th > THETA[0]) & (th < THETA[1]) & (r <= r_max)
    idx = r[sel].astype(int)
    n = np.bincount(idx, minlength=r_max + 1)[:r_max + 1]
    s = np.bincount(idx, img[sel], minlength=r_max + 1)[:r_max + 1]
    return np.arange(r_max + 1) + .5, s / np.maximum(n, 1)


def write_gif(paths, durations_ms, gif, colors=256):
    """GIF from PNG frames with per-frame durations (multiples of a common tick). Uses ffmpeg (one shared
    palette, only changed rectangles stored) when available, otherwise Pillow."""
    from math import gcd
    from functools import reduce
    if shutil.which('ffmpeg'):
        tick = reduce(gcd, [int(d) for d in durations_ms])
        with tempfile.TemporaryDirectory() as tmp:
            n = 0
            for p, d in zip(paths, durations_ms):          # hold a frame by repeating it at the tick rate
                for _ in range(int(round(d / tick))):
                    shutil.copyfile(p, Path(tmp) / f'seq_{n:05d}.png')
                    n += 1
            subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', f'{1000 / tick:g}',
                            '-i', str(Path(tmp) / 'seq_%05d.png'),
                            '-vf', f'split[a][b];[a]palettegen=max_colors={colors}:stats_mode=full[p];'
                                   '[b][p]paletteuse=dither=none:diff_mode=rectangle',
                            '-loop', '0', str(gif)], check=True)
    else:
        frames = [Image.open(p).convert('RGB').quantize(colors=colors, dither=Image.Dither.NONE) for p in paths]
        frames[0].save(gif, save_all=True, append_images=frames[1:], duration=list(durations_ms), loop=0)
    print('Wrote', Path(gif).name, f'{Path(gif).stat().st_size / 1e6:.1f} MB', flush=True)


def method_chi_base():
    """chi2 of each other method's line with SEM and speed reproducibility only (from sagnac.txt), in the
    order cosine fit, template matching, random points."""
    import re
    txt = (UPD / 'sagnac.txt').read_text()
    vals = [float(m) for m in re.findall(r'line with SEM and speed reproducibility: chi2 = ([0-9.]+)', txt)]
    assert len(vals) == 3, vals
    return vals
