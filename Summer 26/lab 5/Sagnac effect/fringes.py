"""Dark-fringe radii and fringe phase of every photograph (manual Eq. 7).

Input: <data>/<S>/<YYYYMMDDhhmmss>.png, one folder per speed setting S (0 = at rest).
  python fringes.py                                    Measurements/ -> fringes_frames.csv
  python fringes.py --data <folder> --out fringes_frames_small.csv     (the small loop)
The two beams leave the interferometer with different wavefront curvature. A Gaussian beam
carries the phase k r^2 / (2 R(z)), so two beams with wavefront radii R1 != R2 at the camera
differ by (pi / lambda)(1/R1 - 1/R2) r^2: the fringes are rings, and their phase at radius r
from the ring centre is
    phi(r) = 2 pi r^2 / P + phi_c,     P = 2 lambda / |1/R1 - 1/R2|  (in camera px^2),
so the rings are equally spaced in s = r^2. A relative tilt only moves the centre. A change of
a lens-fibre distance changes R1 or R2 and so P, not phi_c. Rotation adds the Sagnac phase to
phi_c uniformly; dark ring m sits at s_m = P (m + 1/2) - P phi_c / (2 pi).

Per frame (red channel only; the camera never saturates, max 239/255):
  1. Gaussian smoothing, sigma = SMOOTH_PX, against the speckle.
  2. Ring guesses: minima of the intensity averaged over the sector THETA_DEG around the
     current centre (the lower-left of the pattern carries no fringes).
  3. Along rays every THETA_STEP deg, the minimum of each dark ring is located inside a window
     of +-0.3 local fringe periods around its guess and refined with a parabola; minima shallower
     than MIN_DEPTH counts are dropped.
  4. One centre and one radius per ring are fitted to all the minima (concentric circles,
     soft-L1 loss). Steps 2-4 are repeated three times from the updated centre.
  5. s_m = R_m^2 is regressed on the ring order m: slope P, intercept P (1/2 - phi_c/2pi).
     phi_c (wrapped to (-pi, pi]) is the fringe phase at the ring centre. This is the manual's
     Eq. (7) (shift of each minimum divided by the local fringe period, averaged over minima)
     carried out in s, where the period is the same for every ring, with the period measured
     in each frame.
  6. Second method (prep Q5, the cosine fit): the sector profile in bins of s about the same
     centre, up to S2_FIT_MAX, is fitted with the two-beam model of cos_model; its phase is the
     same phi_c. It uses every bin, but depends on how well the envelope and the camera's
     intensity response follow the model; the minima do not.
The minimum positions do not depend on the camera's intensity response (gamma), since a
monotonic response moves no extremum.

START, P_GUESS and THETA_DEG are set for the camera view of 11 October; check them on new data.

Output: one row per frame, and the r^2 profiles of step 6 (fringes_profiles*.npz) for methods.py.
Run from this folder.
"""
import argparse
import csv
import os
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage, optimize, signal

HERE = Path(__file__).resolve().parent
DATA = HERE / 'Measurements'
OUT = HERE / 'fringes_frames.csv'

SMOOTH_PX = 4.0                     # Gaussian sigma (px) against speckle
THETA_DEG = (-75.0, 180.0)          # angular sector with fringes (deg, anticlockwise from +x, image up)
THETA_STEP = 1.5                    # ray spacing (deg)
R_MAX = 640                         # px
MIN_DEPTH = 3.0                     # counts; shallower minima on a ray are not used
P_GUESS = 107.5e3                   # px^2, typical fringe period in r^2 (sets the search windows)
MIN_POINTS = 30                     # minimum ray minima for a ring to be kept
START = (900.0, 505.0)              # starting centre (px)
N_ITER = 3
MAX_RINGS = 3
S2_BIN = 2000.0                     # px^2, bin width of the r^2 profile (equal pixel counts per bin)
S2_FIT_MAX = 2.6e5                  # px^2, cosine-fit range: the three dark rings

_RR = np.arange(0, R_MAX + 1, 1.0)
_TH = np.radians(np.arange(THETA_DEG[0], THETA_DEG[1], THETA_STEP))


def load_red(path):
    """Red channel as float (counts 0-255)."""
    return np.asarray(Image.open(path))[..., 0].astype(float)


def frame_time(path):
    """Seconds of the day from the YYYYMMDDhhmmss file name."""
    t = Path(path).stem
    return int(t[8:10]) * 3600 + int(t[10:12]) * 60 + int(t[12:14])


def sector_profile(img, cx, cy):
    """Mean intensity in 1 px annuli over the fringe sector, r = 0 ... R_MAX."""
    y, x = np.indices(img.shape)
    r = np.hypot(x - cx, y - cy)
    th = np.degrees(np.arctan2(-(y - cy), x - cx))
    sel = (th > THETA_DEG[0]) & (th < THETA_DEG[1]) & (r <= R_MAX)
    idx = r[sel].astype(int)
    n = np.bincount(idx, minlength=R_MAX + 1)[:R_MAX + 1]
    s = np.bincount(idx, img[sel], minlength=R_MAX + 1)[:R_MAX + 1]
    return s / np.maximum(n, 1)


def sector_profile_s2(img, cx, cy, edges):
    """Mean intensity in bins of s = r^2 over the fringe sector (equal pixel counts per bin)."""
    y, x = np.indices(img.shape)
    s = (x - cx) ** 2 + (y - cy) ** 2
    th = np.degrees(np.arctan2(-(y - cy), x - cx))
    sel = (th > THETA_DEG[0]) & (th < THETA_DEG[1]) & (s < edges[-1])
    idx = np.digitize(s[sel], edges) - 1
    num = np.bincount(idx, img[sel], minlength=len(edges) - 1)[:len(edges) - 1]
    den = np.bincount(idx, minlength=len(edges) - 1)[:len(edges) - 1]
    return num / np.maximum(den, 1)


def ring_guesses(profile):
    """Radii of the dark rings in the sector profile (prominence > 3 % of its maximum)."""
    p = ndimage.gaussian_filter1d(profile, 3)
    pk, _ = signal.find_peaks(-p, prominence=0.03 * p.max())
    return _RR[pk[_RR[pk] > 40]]


def ray_minima(img, cx, cy, guesses):
    """(ring index, x, y) of the minimum of each guessed ring along every ray."""
    pts = []
    for t in _TH:
        xs, ys = cx + _RR * np.cos(t), cy - _RR * np.sin(t)
        p = ndimage.map_coordinates(img, [ys, xs], order=1, mode='nearest')
        for k, g in enumerate(guesses):
            lam = P_GUESS / (2 * g)                     # local fringe period in r
            w = min(0.3 * lam, 0.8 * g)
            lo, hi = int(max(g - w, 1)), int(min(g + w, R_MAX - 1))
            seg = p[lo:hi + 1]
            q = lo + int(np.argmin(seg))
            if q <= lo or q >= hi or min(seg[0], seg[-1]) - p[q] < MIN_DEPTH:
                continue
            h = max(2, int(lam / 10))
            a, b = max(q - h, 0), min(q + h, R_MAX)
            c2, c1, _ = np.polyfit(_RR[a:b + 1], p[a:b + 1], 2)
            r = -c1 / (2 * c2) if c2 > 0 else _RR[q]
            if abs(r - _RR[q]) > h:
                r = _RR[q]
            pts.append((k, cx + r * np.cos(t), cy - r * np.sin(t)))
    return np.array(pts).reshape(-1, 3)


def fit_concentric(pts, cx, cy, guesses):
    """Common centre and one radius per ring. Returns params, covariance, points per ring, rms."""
    lab = pts[:, 0].astype(int)
    n = np.bincount(lab, minlength=len(guesses))
    keep_rings = np.flatnonzero(n >= MIN_POINTS)
    remap = -np.ones(len(guesses), int); remap[keep_rings] = np.arange(len(keep_rings))
    m = remap[lab] >= 0
    L, X, Y = remap[lab[m]], pts[m, 1], pts[m, 2]

    def res(p):
        return np.hypot(X - p[0], Y - p[1]) - p[2:][L]

    p0 = np.r_[cx, cy, np.asarray(guesses, float)[keep_rings]]
    sol = optimize.least_squares(res, p0, loss='soft_l1', f_scale=4.0)
    dof = max(len(sol.fun) - len(sol.x), 1)
    s2 = np.sum(sol.fun ** 2) / dof
    cov = np.linalg.pinv(sol.jac.T @ sol.jac) * s2
    return sol.x, cov, n[keep_rings], np.sqrt(s2)


def ring_phase(R, P_fallback=P_GUESS):
    """Fringe period P (px^2) and phase at the ring centre phi_c (rad, wrapped) from radii R."""
    s = np.sort(np.asarray(R, float)) ** 2
    if len(s) == 1:
        P, s0 = P_fallback, s[0]
    else:
        m = np.r_[0, np.cumsum(np.rint(np.diff(s) / P_GUESS))]     # ring order, gaps allowed
        P, s0 = np.polyfit(m, s, 1)                                 # s0: s of ring order 0
    return P, float(np.angle(np.exp(1j * (np.pi - 2 * np.pi * s0 / P))))


def cos_model(s, sC, C, sB, B, P, phi, D):
    """Two Gaussian beams in s = r^2: background C e^(-s/sC), fringes B e^(-s/sB) cos(2 pi s/P + phi)
    with B >= 0 (bright centre at phi = 0, as in ring_phase), dark floor D."""
    return C * np.exp(-s / sC) + B * np.exp(-s / sB) * np.cos(2 * np.pi * s / P + phi) + D


def cosine_fit(img, cx, cy):
    """Second method (prep Q5, I = C + B cos(2 pi y / Lambda + phi), here in s = r^2): fit cos_model
    to the sector profile up to S2_FIT_MAX, best of four starting phases.
    Returns the parameters of cos_model, the rms residual (counts) and the profile."""
    edges = np.arange(0, S2_FIT_MAX + 1, S2_BIN)
    s = 0.5 * (edges[1:] + edges[:-1])
    p = sector_profile_s2(img, cx, cy, edges)
    bounds = ([1e4, 0, 1e4, 0, 0.7 * P_GUESS, -10, -50], [1e7, 1e3, 1e7, 1e3, 1.3 * P_GUESS, 10, 255])
    best = None
    for phi0 in np.arange(-np.pi, np.pi, np.pi / 2):
        p0 = [1.5e5, p.max() / 2, 1.5e5, p.max() / 2, P_GUESS, phi0, p.min()]
        try:
            popt, _ = optimize.curve_fit(cos_model, s, p, p0=p0, bounds=bounds, maxfev=20000)
        except RuntimeError:                                       # no convergence from this start
            continue
        ssr = np.sum((p - cos_model(s, *popt)) ** 2)
        if best is None or ssr < best[0]:
            best = (ssr, popt)
    ssr, popt = best
    return popt, float(np.sqrt(ssr / (len(s) - len(popt)))), p


def analyse_frame(path, start=START):
    img = ndimage.gaussian_filter(load_red(path), SMOOTH_PX)
    cx, cy = start
    for _ in range(N_ITER):
        g = ring_guesses(sector_profile(img, cx, cy))
        pts = ray_minima(img, cx, cy, g)
        x, cov, npts, rms = fit_concentric(pts, cx, cy, g)
        cx, cy = x[0], x[1]
    R, uR = x[2:], np.sqrt(np.diag(cov))[2:]
    order = np.argsort(R)
    R, uR, npts = R[order][:MAX_RINGS], uR[order][:MAX_RINGS], npts[order][:MAX_RINGS]
    P, phi = ring_phase(R)
    popt, rms_cos, prof = cosine_fit(img, cx, cy)
    P_cos, phi_cos = popt[4], float(np.angle(np.exp(1j * popt[5])))
    return dict(cx=cx, cy=cy, u_cx=np.sqrt(cov[0, 0]), u_cy=np.sqrt(cov[1, 1]), R=R, uR=uR,
                npts=npts, rms_px=rms, P=P, phi=phi, P_cos=P_cos, phi_cos=phi_cos, rms_cos=rms_cos,
                profile=prof)


def frames(data=DATA):
    """(S, path) for every photograph, settings in increasing order, frames in time order."""
    for d in sorted((p for p in Path(data).iterdir() if p.is_dir()), key=lambda p: float(p.name)):
        for f in sorted(d.glob('*.png')):
            yield float(d.name), f


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--data', type=Path, default=DATA, help='folder with one subfolder per setting')
    ap.add_argument('--out', type=Path, default=OUT, help='output CSV (relative to this folder)')
    args = ap.parse_args()
    out = args.out if args.out.is_absolute() else HERE / args.out
    cols = ['S', 'file', 'time_s', 'cx', 'cy', 'u_cx', 'u_cy', 'n_rings',
            'R1', 'R2', 'R3', 'uR1', 'uR2', 'uR3', 'n1', 'n2', 'n3', 'rms_px', 'P_px2', 'phi_rad',
            'P_cos_px2', 'phi_cos_rad', 'rms_cos']
    rows, profiles = [], []
    for S, f in frames(args.data):
        r = analyse_frame(f)
        profiles.append(r['profile'])
        k = len(r['R'])
        pad = lambda v: list(np.round(v, 3)) + [''] * (MAX_RINGS - k)
        rel = Path(os.path.relpath(f.resolve(), HERE)).as_posix()     # path from this folder
        rows.append([S, rel, frame_time(f), round(r['cx'], 2),
                     round(r['cy'], 2), round(r['u_cx'], 3), round(r['u_cy'], 3), k,
                     *pad(r['R']), *pad(r['uR']), *(list(r['npts']) + [''] * (MAX_RINGS - k)),
                     round(r['rms_px'], 2), round(r['P'], 1), round(r['phi'], 5),
                     round(r['P_cos'], 1), round(r['phi_cos'], 5), round(r['rms_cos'], 2)])
        print(f"S={S:6.2f} {f.name}  c=({r['cx']:.1f}, {r['cy']:.1f})  R={np.round(r['R'], 1)} px  "
              f"P={r['P'] / 1e3:.1f}k px^2  phi_c={r['phi']:+.3f} rad  cosine fit {r['phi_cos']:+.3f} rad",
              flush=True)
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        w.writerows(rows)
    edges = np.arange(0, S2_FIT_MAX + 1, S2_BIN)
    prof_out = out.with_name(out.stem.replace('_frames', '_profiles') + '.npz')
    np.savez(prof_out, s=0.5 * (edges[1:] + edges[:-1]), profiles=np.array(profiles),
             files=np.array([r[1] for r in rows]))
    print(f'{len(rows)} frames -> {out.name}, profiles -> {prof_out.name}')


if __name__ == '__main__':
    main()
