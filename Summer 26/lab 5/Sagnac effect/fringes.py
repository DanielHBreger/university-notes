"""Fringe phase of every photograph: where the dark rings are, and what phase that means.

WHAT THIS SCRIPT IS FOR
Each photograph shows bright and dark rings. Rotation shifts all of them together, and the
size of that shift is the Sagnac phase. This script reduces each of the 161 photographs to a
few numbers, above all phi_c, the fringe phase at the centre of the rings. sagnac.py then
compares phi_c between speeds. Input: Measurements/<S>/<YYYYMMDDhhmmss>.png, one folder per
speed setting S (0 = at rest).   python fringes.py   (a few minutes)

HOW IT RUNS
main() loops over the photographs (frames(): settings in increasing order, each in time order)
and calls analyse_frame() on each. Per photograph:
  1. load_red + Gaussian smoothing (SMOOTH_PX) against speckle. Only the red channel carries
     light, and the camera never saturates (maximum 239 of 255).
  2. ring_guesses(sector_profile(...)): average the intensity in thin rings (annuli) around the
     current centre guess, over the angular sector THETA_DEG where fringes are visible (the
     lower left of the pattern has none); the dips of that profile are the dark-ring radii.
  3. ray_minima: along rays from the centre every THETA_STEP degrees, find each dark ring's
     minimum inside a window of +-0.3 local fringe periods around its guess, and refine it with
     a parabola through the neighbouring pixels. Minima shallower than MIN_DEPTH are dropped.
  4. fit_concentric: fit one common centre and one radius per ring to all those minima
     (concentric circles; a robust loss, so stray points have little pull).
     Steps 2-4 are repeated N_ITER = 3 times, each time from the improved centre.
  5. ring_phase: from the ring radii R_m, the fringe period P and the centre phase phi_c
     (see "WHY RINGS" below). This is the main method ("minima").
  6. cosine_fit: second method (prep Q5). The intensity profile against s = r^2 around the same
     centre is fitted with the two-beam model cos_model; its phase is the same phi_c.
main() writes one row per photograph to fringes_frames.csv and the step-6 profiles to
fringes_profiles.npz (methods.py reads both).

WHY RINGS, AND WHY THEY ARE EVENLY SPACED IN r^2
The two beams leave the interferometer with different wavefront curvature. A Gaussian beam
carries the phase k r^2 / (2 R(z)), so two beams with wavefront radii R1 != R2 at the camera
differ by (pi / lambda)(1/R1 - 1/R2) r^2: the fringes are rings, and their phase at radius r
from the ring centre is
    phi(r) = 2 pi r^2 / P + phi_c,     P = 2 lambda / |1/R1 - 1/R2|  (in camera px^2),
so the rings are equally spaced in s = r^2. A relative tilt only moves the centre. A change of
a lens-fibre distance changes R1 or R2 and so P, not phi_c. Rotation adds the Sagnac phase to
phi_c uniformly; dark ring m sits at s_m = P (m + 1/2) - P phi_c / (2 pi).
So a straight-line fit of s_m against the ring number m gives P (slope) and phi_c (from the
intercept). This is the manual's Eq. (7) (shift of each minimum divided by the local fringe
period, averaged over minima) carried out in s, where the period is the same for every ring,
with the period measured in each frame.

WHY MINIMA RATHER THAN BRIGHTNESS
The position of a minimum does not depend on the camera's intensity response (gamma) or on
the beam's brightness envelope, since a monotonic response moves no extremum. The cosine fit
uses every bin of the profile, but it depends on how well the envelope and the camera's
response follow cos_model; the minima do not.

START, P_GUESS and THETA_DEG are set for the camera view of 11 October.
Run from this folder.
"""
import csv
import os
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage, optimize, signal

HERE = Path(__file__).resolve().parent
DATA = HERE / 'Measurements'
OUT = HERE / 'fringes_frames.csv'

# Settings of the ring finder, chosen for the camera view of 11 October
SMOOTH_PX = 4.0                     # Gaussian sigma (px) against speckle
THETA_DEG = (-75.0, 180.0)          # angular sector with fringes (deg, anticlockwise from +x, image up)
THETA_STEP = 1.5                    # ray spacing (deg): 170 rays across the sector
R_MAX = 640                         # px, how far out from the centre the rays go
MIN_DEPTH = 3.0                     # counts; shallower minima on a ray are not used
P_GUESS = 107.5e3                   # px^2, typical fringe period in r^2 (sets the search windows)
MIN_POINTS = 30                     # minimum ray minima for a ring to be kept
START = (900.0, 505.0)              # starting centre (px, x and y)
N_ITER = 3                          # rounds of "profile -> minima -> circle fit" per photograph
MAX_RINGS = 3                       # at most the three innermost dark rings are recorded
S2_BIN = 2000.0                     # px^2, bin width of the r^2 profile (equal pixel counts per bin)
S2_FIT_MAX = 2.6e5                  # px^2, cosine-fit range: the three dark rings

_RR = np.arange(0, R_MAX + 1, 1.0)                                  # radii sampled along a ray (px)
_TH = np.radians(np.arange(THETA_DEG[0], THETA_DEG[1], THETA_STEP))  # ray angles (rad)


def load_red(path):
    """Red channel of the PNG as float (counts 0-255); the green and blue channels are dark."""
    return np.asarray(Image.open(path))[..., 0].astype(float)


def frame_time(path):
    """Seconds since midnight from the YYYYMMDDhhmmss file name (characters 8-13 = hhmmss)."""
    t = Path(path).stem
    return int(t[8:10]) * 3600 + int(t[10:12]) * 60 + int(t[12:14])


def sector_profile(img, cx, cy):
    """Step 2. Mean intensity in 1 px wide annuli around (cx, cy), r = 0 ... R_MAX, using only the
    pixels inside the fringe sector THETA_DEG. Returns one value per integer radius."""
    y, x = np.indices(img.shape)                       # pixel coordinates of every pixel
    r = np.hypot(x - cx, y - cy)
    th = np.degrees(np.arctan2(-(y - cy), x - cx))      # minus: image y grows downwards
    sel = (th > THETA_DEG[0]) & (th < THETA_DEG[1]) & (r <= R_MAX)
    idx = r[sel].astype(int)                           # annulus number of each selected pixel
    n = np.bincount(idx, minlength=R_MAX + 1)[:R_MAX + 1]               # pixels per annulus
    s = np.bincount(idx, img[sel], minlength=R_MAX + 1)[:R_MAX + 1]     # summed intensity
    return s / np.maximum(n, 1)


def sector_profile_s2(img, cx, cy, edges):
    """The same sector average, but in bins of s = r^2 (bin edges `edges`, px^2). Equal steps in
    r^2 enclose equal areas, so every bin has about the same number of pixels, and the fringes
    are evenly spaced along this axis. Used by the cosine fit, methods.py and patterns.py."""
    y, x = np.indices(img.shape)
    s = (x - cx) ** 2 + (y - cy) ** 2
    th = np.degrees(np.arctan2(-(y - cy), x - cx))
    sel = (th > THETA_DEG[0]) & (th < THETA_DEG[1]) & (s < edges[-1])
    idx = np.digitize(s[sel], edges) - 1
    num = np.bincount(idx, img[sel], minlength=len(edges) - 1)[:len(edges) - 1]
    den = np.bincount(idx, minlength=len(edges) - 1)[:len(edges) - 1]
    return num / np.maximum(den, 1)


def ring_guesses(profile):
    """Step 2. Rough dark-ring radii: the dips of the (lightly smoothed) sector profile that are
    deeper than 3 % of its maximum. Dips within 40 px of the centre are ignored (too few pixels
    there to be reliable)."""
    p = ndimage.gaussian_filter1d(profile, 3)
    pk, _ = signal.find_peaks(-p, prominence=0.03 * p.max())   # minima of p = peaks of -p
    return _RR[pk[_RR[pk] > 40]]


def ray_minima(img, cx, cy, guesses):
    """Step 3. For each ray from (cx, cy) and each guessed ring, the point where the intensity
    along the ray is lowest near that ring. Returns an array of (ring index, x, y).

    The search window is +-0.3 of the local fringe period around the guess. Because the rings
    are evenly spaced in r^2 with period P, the period in r at radius g is about P / (2 g)
    (from d(r^2) = 2 r dr). The lowest sample is refined by fitting a parabola through the
    samples around it; its vertex gives a sub-pixel radius."""
    pts = []
    for t in _TH:
        xs, ys = cx + _RR * np.cos(t), cy - _RR * np.sin(t)
        p = ndimage.map_coordinates(img, [ys, xs], order=1, mode='nearest')   # intensity along the ray
        for k, g in enumerate(guesses):
            lam = P_GUESS / (2 * g)                     # local fringe period in r
            w = min(0.3 * lam, 0.8 * g)
            lo, hi = int(max(g - w, 1)), int(min(g + w, R_MAX - 1))
            seg = p[lo:hi + 1]
            q = lo + int(np.argmin(seg))
            # reject: minimum at the window edge (not a real dip) or too shallow (noise)
            if q <= lo or q >= hi or min(seg[0], seg[-1]) - p[q] < MIN_DEPTH:
                continue
            h = max(2, int(lam / 10))                   # parabola half-width: a tenth of a period
            a, b = max(q - h, 0), min(q + h, R_MAX)
            c2, c1, _ = np.polyfit(_RR[a:b + 1], p[a:b + 1], 2)
            r = -c1 / (2 * c2) if c2 > 0 else _RR[q]    # vertex of c2 r^2 + c1 r + c0
            if abs(r - _RR[q]) > h:                     # vertex outside the fitted points: keep the sample
                r = _RR[q]
            pts.append((k, cx + r * np.cos(t), cy - r * np.sin(t)))
    return np.array(pts).reshape(-1, 3)


def fit_concentric(pts, cx, cy, guesses):
    """Step 4. Fit concentric circles to the ray minima: parameters (cx, cy, R_1, R_2, ...), one
    common centre and one radius per ring, minimising each point's distance from its circle.
    Rings with fewer than MIN_POINTS minima are dropped. The 'soft_l1' loss behaves like least
    squares for small residuals but grows only linearly for large ones, so a few stray minima
    (speckle) do not pull the circles. Returns the parameters, their covariance (from the
    Jacobian, scaled by the residual variance), points per ring, and the rms distance (px)."""
    lab = pts[:, 0].astype(int)                        # which ring each point belongs to
    n = np.bincount(lab, minlength=len(guesses))
    keep_rings = np.flatnonzero(n >= MIN_POINTS)
    remap = -np.ones(len(guesses), int); remap[keep_rings] = np.arange(len(keep_rings))
    m = remap[lab] >= 0                                # points of the rings that are kept
    L, X, Y = remap[lab[m]], pts[m, 1], pts[m, 2]

    def res(p):                                        # distance from the centre minus the ring's radius
        return np.hypot(X - p[0], Y - p[1]) - p[2:][L]

    p0 = np.r_[cx, cy, np.asarray(guesses, float)[keep_rings]]
    sol = optimize.least_squares(res, p0, loss='soft_l1', f_scale=4.0)
    dof = max(len(sol.fun) - len(sol.x), 1)
    s2 = np.sum(sol.fun ** 2) / dof
    cov = np.linalg.pinv(sol.jac.T @ sol.jac) * s2
    return sol.x, cov, n[keep_rings], np.sqrt(s2)


def ring_phase(R, P_fallback=P_GUESS):
    """Step 5, the main method. Fringe period P (px^2) and the phase at the ring centre phi_c
    (rad, wrapped to (-pi, pi]) from the dark-ring radii R.

    Dark ring m sits at s_m = P (m + 1/2) - P phi_c / (2 pi). So a straight-line fit of s_m = R_m^2
    against the ring number m has slope P and intercept s0 = P/2 - P phi_c/(2 pi), which gives
    phi_c = pi - 2 pi s0 / P. Ring numbers come from the gaps in r^2 divided by P_GUESS, so a ring
    missed in between does not shift the numbering. With a single ring, P_GUESS is used for P."""
    s = np.sort(np.asarray(R, float)) ** 2
    if len(s) == 1:
        P, s0 = P_fallback, s[0]
    else:
        m = np.r_[0, np.cumsum(np.rint(np.diff(s) / P_GUESS))]     # ring order, gaps allowed
        P, s0 = np.polyfit(m, s, 1)                                 # s0: s of ring order 0
    # angle(exp(i x)) wraps x into (-pi, pi]
    return P, float(np.angle(np.exp(1j * (np.pi - 2 * np.pi * s0 / P))))


def cos_model(s, sC, C, sB, B, P, phi, D):
    """Two-beam intensity against s = r^2. Two Gaussian beams of widths w1, w2 interfere as
    I1 + I2 + 2 sqrt(I1 I2) cos(phase). A Gaussian intensity exp(-2 r^2 / w^2) is an exponential
    in s, so: background C e^(-s/sC) (the sum of the two beams, approximated by one exponential),
    fringe amplitude B e^(-s/sB), fringes cos(2 pi s/P + phi), dark floor D. B >= 0, so phi = 0
    means a bright centre, the same convention as ring_phase."""
    return C * np.exp(-s / sC) + B * np.exp(-s / sB) * np.cos(2 * np.pi * s / P + phi) + D


def cosine_fit(img, cx, cy):
    """Step 6, the second method (prep Q5: I = C + B cos(2 pi y / Lambda + phi), here in s = r^2).
    Build the sector profile in bins of s up to S2_FIT_MAX (the three dark rings), fit cos_model
    to it from four starting phases (a fit can get stuck in a wrong fringe from a poor start),
    and keep the best. Returns the parameters of cos_model, the rms residual (counts) and the
    profile, which is saved for methods.py."""
    edges = np.arange(0, S2_FIT_MAX + 1, S2_BIN)
    s = 0.5 * (edges[1:] + edges[:-1])                 # bin centres
    p = sector_profile_s2(img, cx, cy, edges)
    # bounds in the order of cos_model's parameters: sC, C, sB, B, P (+-30 % of P_GUESS), phi, D
    bounds = ([1e4, 0, 1e4, 0, 0.7 * P_GUESS, -10, -50], [1e7, 1e3, 1e7, 1e3, 1.3 * P_GUESS, 10, 255])
    best = None
    for phi0 in np.arange(-np.pi, np.pi, np.pi / 2):    # starting phases -pi, -pi/2, 0, pi/2
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
    """Steps 1-6 for one photograph (see the module docstring). Returns the ring centre and its
    uncertainty, the (up to three) dark-ring radii, the period P and phase phi of the minima
    method, the same from the cosine fit, the fit quality of both, and the r^2 profile."""
    img = ndimage.gaussian_filter(load_red(path), SMOOTH_PX)              # step 1
    cx, cy = start
    for _ in range(N_ITER):                                               # steps 2-4, three rounds
        g = ring_guesses(sector_profile(img, cx, cy))
        pts = ray_minima(img, cx, cy, g)
        x, cov, npts, rms = fit_concentric(pts, cx, cy, g)
        cx, cy = x[0], x[1]                            # x = (cx, cy, R_1, R_2, ...)
    R, uR = x[2:], np.sqrt(np.diag(cov))[2:]
    order = np.argsort(R)                              # innermost ring first, keep MAX_RINGS
    R, uR, npts = R[order][:MAX_RINGS], uR[order][:MAX_RINGS], npts[order][:MAX_RINGS]
    P, phi = ring_phase(R)                                                # step 5
    popt, rms_cos, prof = cosine_fit(img, cx, cy)                         # step 6
    P_cos, phi_cos = popt[4], float(np.angle(np.exp(1j * popt[5])))       # popt[5] = phi, wrapped
    return dict(cx=cx, cy=cy, u_cx=np.sqrt(cov[0, 0]), u_cy=np.sqrt(cov[1, 1]), R=R, uR=uR,
                npts=npts, rms_px=rms, P=P, phi=phi, P_cos=P_cos, phi_cos=phi_cos, rms_cos=rms_cos,
                profile=prof)


def frames(data=DATA):
    """(S, path) for every photograph, settings in increasing order, frames in time order."""
    for d in sorted((p for p in Path(data).iterdir() if p.is_dir()), key=lambda p: float(p.name)):
        for f in sorted(d.glob('*.png')):
            yield float(d.name), f


def main():
    """Analyse every photograph and write fringes_frames.csv (one row per photograph; columns:
    setting S, file, time of day, ring centre and its errors, number of rings, radii R1-R3 and
    their errors, minima per ring, circle-fit rms, P and phi of the minima method, P and phi of
    the cosine fit, its rms) and fringes_profiles.npz (the r^2 profiles, for methods.py)."""
    out = OUT
    cols = ['S', 'file', 'time_s', 'cx', 'cy', 'u_cx', 'u_cy', 'n_rings',
            'R1', 'R2', 'R3', 'uR1', 'uR2', 'uR3', 'n1', 'n2', 'n3', 'rms_px', 'P_px2', 'phi_rad',
            'P_cos_px2', 'phi_cos_rad', 'rms_cos']
    rows, profiles = [], []
    for S, f in frames(DATA):
        r = analyse_frame(f)
        profiles.append(r['profile'])
        k = len(r['R'])                                # rings found (usually 2 or 3)
        pad = lambda v: list(np.round(v, 3)) + [''] * (MAX_RINGS - k)   # empty cells for missing rings
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
    prof_out = HERE / 'fringes_profiles.npz'
    np.savez(prof_out, s=0.5 * (edges[1:] + edges[:-1]), profiles=np.array(profiles),
             files=np.array([r[1] for r in rows]))
    print(f'{len(rows)} frames -> {out.name}, profiles -> {prof_out.name}')


if __name__ == '__main__':
    main()
