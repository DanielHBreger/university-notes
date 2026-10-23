"""Fringe phase of every photograph: where the dark rings are, and what phase that means.

Each photograph shows bright and dark rings. Rotation shifts all of them together, and the size
of that shift is the Sagnac phase. This script reduces each of the 161 photographs to a few
numbers, above all phi_c, the fringe phase at the centre of the rings, which sagnac.py compares
between speeds. Input: Measurements/<S>/<YYYYMMDDhhmmss>.png, one folder per speed setting S
(0 = at rest).   python fringes.py   (a few minutes)

HOW ONE PHOTOGRAPH IS ANALYSED (analyse_photo)
  1. Red channel, blurred by SMOOTH_PX against speckle (only the red channel carries light; the
     camera never saturates, maximum 239 of 255).
  2. ring_guesses: the mean intensity in thin annuli around the current centre, over the angular
     sector SECTOR_DEG where fringes are visible (the lower left has none); its dips are rough
     dark-ring radii.
  3. ray_minima: along rays from the centre every RAY_STEP_DEG, the darkest point near each
     guessed ring, refined with a parabola through the neighbouring pixels.
  4. fit_circles: one common centre and one radius per ring through all those minima.
     Steps 2-4 run N_ROUNDS times, each round from the improved centre.
  5. ring_phase: from the ring radii, the fringe period P and the centre phase phi_c (below).
     This is the main method ("ring minima").
  6. cosine_fit: the second method (prep Q5), a two-beam model fitted to the intensity profile
     against s = r^2 around the same centre; its phase is the same phi_c.
main() writes one row per photograph to fringes_frames.csv and the step-6 profiles to
fringes_profiles.npz (methods.py reads both).

WHY RINGS, AND WHY THEY ARE EVENLY SPACED IN r^2
The two beams leave the interferometer with different wavefront curvature. A Gaussian beam
carries the phase k r^2 / (2 R(z)), so two beams with wavefront radii R1 != R2 at the camera
differ by (pi / lambda)(1/R1 - 1/R2) r^2: the fringes are rings with the phase
    phi(r) = 2 pi r^2 / P + phi_c,     P = 2 lambda / |1/R1 - 1/R2|  (in camera px^2),
equally spaced in s = r^2. A relative tilt only moves the centre; a change of a lens-fibre
distance changes P, not phi_c. Rotation adds the Sagnac phase to phi_c; dark ring m sits at
s_m = P (m + 1/2) - P phi_c / (2 pi), so a straight line of s_m against m gives P (slope) and
phi_c (intercept). This is the manual's Eq. (7) (shift of each minimum divided by the local
fringe period, averaged over minima) done in s, where the period is the same for every ring.

WHY MINIMA RATHER THAN BRIGHTNESS
The position of a minimum does not depend on the camera's intensity response or on the beam's
brightness envelope: a monotonic response moves no extremum. The cosine fit uses every bin of
the profile but depends on how well the envelope and the camera follow its model.

START_CENTRE, P_GUESS and SECTOR_DEG are set for the camera view of 11 October. Run from this folder.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage, optimize, signal

HERE = Path(__file__).resolve().parent
PHOTOS = HERE / 'Measurements'

# Settings of the ring finder, chosen for the camera view of 11 October
SMOOTH_PX = 4.0                   # Gaussian blur of the image (px) against speckle
SECTOR_DEG = (-75.0, 180.0)       # angular sector with fringes (deg, anticlockwise from +x, image up)
RAY_STEP_DEG = 1.5                # one ray every 1.5 deg: 170 rays across the sector
R_MAX = 640                       # px, how far from the centre the rays go
MIN_DEPTH = 3.0                   # counts; shallower dips along a ray are not used
P_GUESS = 107.5e3                 # px^2, typical fringe period in r^2 (sets the search windows)
MIN_POINTS = 30                   # a ring needs at least this many ray minima to be kept
START_CENTRE = (900.0, 505.0)     # px (x, y), centre for the first round
N_ROUNDS = 3                      # rounds of "guesses -> minima -> circle fit" per photograph
MAX_RINGS = 3                     # at most the three innermost dark rings are recorded
S_BIN = 2000.0                    # px^2, bin width of the r^2 profile
S_MAX = 2.6e5                     # px^2, end of the r^2 profile: the three dark rings

RAY_RADII = np.arange(0, R_MAX + 1, 1.0)                                       # px, samples along a ray
RAY_ANGLES = np.radians(np.arange(SECTOR_DEG[0], SECTOR_DEG[1], RAY_STEP_DEG))
S_EDGES = np.arange(0, S_MAX + 1, S_BIN)                                       # bins of the r^2 profile
S_CENTRES = 0.5 * (S_EDGES[1:] + S_EDGES[:-1])


def load_red(path):
    """Red channel of the PNG as floats (counts 0-255); the green and blue channels are dark."""
    return np.asarray(Image.open(path))[..., 0].astype(float)


def photo_time(path):
    """Seconds since midnight from the YYYYMMDDhhmmss file name."""
    name = Path(path).stem
    return int(name[8:10]) * 3600 + int(name[10:12]) * 60 + int(name[12:14])


def wrap_angle(phi):
    """phi wrapped into (-pi, pi]."""
    return np.angle(np.exp(1j * phi))


def sector_profile(image, cx, cy, s_edges):
    """Mean intensity in rings around (cx, cy), binned in s = r^2 by s_edges (px^2), over the pixels
    inside the fringe sector only. Equal steps in r^2 enclose equal areas, so every bin holds about
    the same number of pixels; edges (0, 1, 4, 9, ...) give 1 px wide annuli in r."""
    y, x = np.indices(image.shape)
    s = (x - cx) ** 2 + (y - cy) ** 2
    angle = np.degrees(np.arctan2(-(y - cy), x - cx))       # minus: image y grows downwards
    inside = (angle > SECTOR_DEG[0]) & (angle < SECTOR_DEG[1]) & (s < s_edges[-1])
    which_bin = np.digitize(s[inside], s_edges) - 1
    n_bins = len(s_edges) - 1
    total = np.bincount(which_bin, image[inside], minlength=n_bins)
    count = np.bincount(which_bin, minlength=n_bins)
    return total / np.maximum(count, 1)


def ring_guesses(image, cx, cy):
    """Step 2. Rough dark-ring radii (px): the dips of the lightly smoothed intensity profile
    against r that are deeper than 3 % of its maximum, beyond 40 px from the centre (too few
    pixels in an annulus closer in)."""
    profile = ndimage.gaussian_filter1d(sector_profile(image, cx, cy, RAY_RADII ** 2), 3)
    dips, _ = signal.find_peaks(-profile, prominence=0.03 * profile.max())     # minima = peaks of -profile
    return dips[dips > 40].astype(float)


def ray_minima(image, cx, cy, guesses):
    """Step 3. For every ray from the centre and every guessed ring, the darkest sample near that
    ring, refined to a fraction of a pixel by a parabola through the neighbouring samples.
    Returns an array of (ring index, x, y)."""
    minima = []
    for angle in RAY_ANGLES:
        xs, ys = cx + RAY_RADII * np.cos(angle), cy - RAY_RADII * np.sin(angle)
        intensity = ndimage.map_coordinates(image, [ys, xs], order=1, mode='nearest')
        for ring, guess in enumerate(guesses):
            period_r = P_GUESS / (2 * guess)                  # local fringe period in r: d(r^2) = 2 r dr
            half_window = min(0.3 * period_r, 0.8 * guess)
            lo, hi = int(max(guess - half_window, 1)), int(min(guess + half_window, R_MAX - 1))
            window = intensity[lo:hi + 1]
            i_min = lo + int(np.argmin(window))
            # a minimum at the window edge is not a dip; a dip shallower than MIN_DEPTH is noise
            if i_min <= lo or i_min >= hi or min(window[0], window[-1]) - intensity[i_min] < MIN_DEPTH:
                continue
            half_fit = max(2, int(period_r / 10))             # parabola over a tenth of a period each side
            fit_lo, fit_hi = max(i_min - half_fit, 0), min(i_min + half_fit, R_MAX)
            c2, c1, _ = np.polyfit(RAY_RADII[fit_lo:fit_hi + 1], intensity[fit_lo:fit_hi + 1], 2)
            r = -c1 / (2 * c2) if c2 > 0 else RAY_RADII[i_min]            # vertex of c2 r^2 + c1 r + c0
            if abs(r - RAY_RADII[i_min]) > half_fit:          # vertex outside the fitted samples
                r = RAY_RADII[i_min]
            minima.append((ring, cx + r * np.cos(angle), cy - r * np.sin(angle)))
    return np.array(minima).reshape(-1, 3)


def fit_circles(minima, cx, cy, guesses):
    """Step 4. Concentric circles through the ray minima: one common centre and one radius per
    ring, by least squares on each point's distance from its circle. Rings with fewer than
    MIN_POINTS minima are dropped; loss='soft_l1' lets a few stray minima (speckle) pull little.
    Returns the centre (x, y), the radii and the rms distance of the minima from the circles."""
    ring = minima[:, 0].astype(int)
    kept_rings = np.flatnonzero(np.bincount(ring, minlength=len(guesses)) >= MIN_POINTS)
    keep = np.isin(ring, kept_rings)
    ring_index = np.searchsorted(kept_rings, ring[keep])     # 0, 1, ... among the kept rings
    x, y = minima[keep, 1], minima[keep, 2]

    def distance_from_circle(params):                        # params = (cx, cy, R_1, R_2, ...)
        return np.hypot(x - params[0], y - params[1]) - params[2:][ring_index]

    start = np.r_[cx, cy, guesses[kept_rings]]
    fit = optimize.least_squares(distance_from_circle, start, loss='soft_l1', f_scale=4.0)
    rms = np.sqrt(np.sum(fit.fun ** 2) / max(len(fit.fun) - len(fit.x), 1))
    return fit.x[:2], fit.x[2:], rms


def ring_numbers(squared_radii):
    """Dark-ring order m of sorted squared radii: the innermost ring found is m = 0 and the gaps in
    r^2 are counted in units of P_GUESS, so a ring missed in between does not shift the numbering."""
    return np.r_[0, np.cumsum(np.rint(np.diff(squared_radii) / P_GUESS))]


def ring_phase(radii):
    """Step 5, the main method. Fringe period P (px^2) and the phase at the ring centre phi_c (rad,
    wrapped) from the dark-ring radii. Dark ring m sits at R_m^2 = P (m + 1/2) - P phi_c / (2 pi):
    a straight line of R_m^2 against m has slope P and intercept s0 = P/2 - P phi_c / (2 pi), so
    phi_c = pi - 2 pi s0 / P. With a single ring P_GUESS stands in for P."""
    s = np.sort(radii) ** 2
    if len(s) == 1:
        P, s0 = P_GUESS, s[0]
    else:
        P, s0 = np.polyfit(ring_numbers(s), s, 1)
    return P, float(wrap_angle(np.pi - 2 * np.pi * s0 / P))


def cosine_model(s, sC, C, sB, B, P, phi, D):
    """Two-beam intensity against s = r^2: background C e^(-s/sC) (two Gaussian beams are
    exponentials in s, approximated by one), fringe amplitude B e^(-s/sB), fringes
    cos(2 pi s/P + phi), dark floor D. B >= 0, so phi = 0 is a bright centre, as in ring_phase."""
    return C * np.exp(-s / sC) + B * np.exp(-s / sB) * np.cos(2 * np.pi * s / P + phi) + D


def cosine_fit(profile):
    """Step 6, the second method (prep Q5: I = C + B cos(2 pi y / Lambda + phi), here in s = r^2).
    Fit cosine_model to the r^2 profile from four starting phases (a fit can lock onto the wrong
    fringe from a poor start) and keep the best. Returns P, the phase (wrapped) and the rms
    residual in counts."""
    #        sC   C    sB   B    P               phi  D
    lower = [1e4, 0, 1e4, 0, 0.7 * P_GUESS, -10, -50]
    upper = [1e7, 1e3, 1e7, 1e3, 1.3 * P_GUESS, 10, 255]
    best = None
    for phi0 in (-np.pi, -np.pi / 2, 0.0, np.pi / 2):
        start = [1.5e5, profile.max() / 2, 1.5e5, profile.max() / 2, P_GUESS, phi0, profile.min()]
        try:
            params, _ = optimize.curve_fit(cosine_model, S_CENTRES, profile, p0=start,
                                           bounds=(lower, upper), maxfev=20000)
        except RuntimeError:                                 # no convergence from this start
            continue
        ssr = np.sum((profile - cosine_model(S_CENTRES, *params)) ** 2)
        if best is None or ssr < best[0]:
            best = (ssr, params)
    ssr, params = best
    P, phi = params[4], params[5]
    return P, float(wrap_angle(phi)), float(np.sqrt(ssr / (len(profile) - len(params))))


def analyse_photo(path):
    """Steps 1-6 for one photograph. Returns the ring centre, the (up to MAX_RINGS) dark-ring
    radii, the circle-fit rms, P and phi_c by the ring minima, P, phi_c and the rms by the cosine
    fit, the r^2 profile, and the ray minima of the last round (for the extraction figure)."""
    image = ndimage.gaussian_filter(load_red(path), SMOOTH_PX)                       # step 1
    cx, cy = START_CENTRE
    for _ in range(N_ROUNDS):                                                         # steps 2-4
        guesses = ring_guesses(image, cx, cy)
        minima = ray_minima(image, cx, cy, guesses)
        (cx, cy), radii, rms = fit_circles(minima, cx, cy, guesses)
    radii = np.sort(radii)[:MAX_RINGS]
    P, phi = ring_phase(radii)                                                        # step 5
    profile = sector_profile(image, cx, cy, S_EDGES)
    P_cos, phi_cos, rms_cos = cosine_fit(profile)                                     # step 6
    return dict(cx=cx, cy=cy, radii=radii, rms_px=rms, P=P, phi=phi, P_cos=P_cos, phi_cos=phi_cos,
                rms_cos=rms_cos, profile=profile, minima=minima)


def photos():
    """(S, path) of every photograph: settings in increasing order, photographs in time order."""
    for folder in sorted((p for p in PHOTOS.iterdir() if p.is_dir()), key=lambda p: float(p.name)):
        for path in sorted(folder.glob('*.png')):
            yield float(folder.name), path


def main():
    rows, profiles = [], []
    for S, path in photos():
        result = analyse_photo(path)
        radii = list(np.round(result['radii'], 3)) + [np.nan] * (MAX_RINGS - len(result['radii']))
        rows.append(dict(S=S, file=path.relative_to(HERE).as_posix(), time_s=photo_time(path),
                         cx=round(result['cx'], 2), cy=round(result['cy'], 2), n_rings=len(result['radii']),
                         R1=radii[0], R2=radii[1], R3=radii[2], rms_px=round(result['rms_px'], 2),
                         P_px2=round(result['P'], 1), phi_rad=round(result['phi'], 5),
                         P_cos_px2=round(result['P_cos'], 1), phi_cos_rad=round(result['phi_cos'], 5),
                         rms_cos=round(result['rms_cos'], 2)))
        profiles.append(result['profile'])
        print(f"S={S:6.2f} {path.name}  centre=({result['cx']:.1f}, {result['cy']:.1f})  "
              f"R={np.round(result['radii'], 1)} px  P={result['P'] / 1e3:.1f}k px^2  "
              f"phi_c={result['phi']:+.3f} rad  cosine fit {result['phi_cos']:+.3f} rad", flush=True)
    pd.DataFrame(rows).to_csv(HERE / 'fringes_frames.csv', index=False)
    np.savez(HERE / 'fringes_profiles.npz', s=S_CENTRES, profiles=np.array(profiles),
             files=np.array([row['file'] for row in rows]))
    print(f'{len(rows)} photographs -> fringes_frames.csv, profiles -> fringes_profiles.npz')


if __name__ == '__main__':
    main()
