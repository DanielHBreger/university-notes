"""Three more ways to read the fringe phase, as cross-checks of fringes.py (minima, cosine fit).

WHAT THIS SCRIPT IS FOR
fringes.py reads the phase from the positions of the dark rings (and from a cosine fit). If the
Sagnac slope is real, other ways of reading the same photographs must give about the same
slope. This script implements three independent ones and writes their phase for every
photograph to methods_frames.csv; sagnac.py fits all methods in the same way and compares them.
Input: fringes_frames.csv and fringes_profiles.npz (fringes.py): for every photograph, the ring
centre and the intensity averaged over the fringe sector in bins of s = r^2 about that centre.

HOW IT RUNS (main(), top to bottom)
  1. Load the frame table, the r^2 profiles M (one row per photograph) and the calibration.
  2. template_phases   template matching (below)               -> phi_t
  3. aia_phases        phase shifting with unknown steps (AIA)  -> phi_aia
  4. unwrap_image + fourier_phases   2D Fourier sideband        -> phi_f
  5. Save the three phases per photograph (methods_frames.csv) and draw the six-panel
     figure methods_illustration (one example of each method on real data).
  6. Robustness scans: rerun each method with other choices (comparison range, window,
     brightness factor, filter width) and record the resulting slope (methods.json, methods.txt).

THE THREE METHODS
  Template matching. The rest frames (S = 0) are averaged into a template T(s). Each profile
  is matched to g T(s + D) + c by least squares over the shift D, with gain g and offset c free.
  In s, an exponential envelope shifted by D only changes by a constant factor, which g absorbs.
  Rings that move inwards give D > 0, and the phase relative to rest is 2 pi D / P_T, where P_T
  is the shift that maps the template onto itself one fringe further out. No fringe model.
  In words: slide the at-rest profile along r^2 until it lies on top of the new one.

  Phase shifting with unknown steps (AIA, Wang & Han, Opt. Lett. 29, 1671 (2004)). Bin j of
  frame k is I_kj = a_j + b_j cos(delta_k - psi_j): every bin brightens and dims with the frame's
  phase delta_k. Alternating least squares: given delta, each bin's a, b cos psi, b sin psi;
  given those, each frame's cos delta, sin delta. Started from the first two principal
  components of the profiles. No fringe period, envelope or cosine shape in r is assumed;
  psi(s) comes out of the fit and should be linear in s (rings equally spaced in r^2). delta is
  defined up to an offset and a sign; the sign is taken so that delta grows with S.
  In words: every point of the profile goes up and down as the fringes pass; find the one
  "clock" per photograph (delta) that explains all those oscillations at once.

  2D Fourier sideband (Takeda), after the friend's sagnac_2d.py. A single carrier does not fit
  rings, so each image is first resampled on (theta, s = r^2) about its ring centre, where the
  rings are straight fringes. Slow background removed along s, the positive sideband around
  the carrier of the mean rest map kept with a Gaussian filter, and the phase relative to rest
  taken as arg sum Z_k conj(Z_rest) over the map, which lines up the local phase before summing.
  In words: unroll the rings into straight stripes, keep only the stripe frequency, and read
  how far the stripes moved from the phase of that frequency component.

All three give an average phase over the window rather than the phase at the ring centre, so a
change of the fringe period between frames (fringes.py) enters them, but not the minima or the
cosine fit. Each is also run with other comparison ranges (template), windows and a per-frame
brightness factor (AIA), or sideband widths (Fourier), to show how much its slope depends on
those choices.

Outputs: methods_frames.csv (phi_template_rad, phi_aia_rad, phi_fourier_rad per frame), methods.txt,
methods.json (the robustness scans), figures/methods_illustration.pdf/.png (figure*).
Run from this folder after calibration.py and fringes.py, before sagnac.py (a few minutes).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import interpolate, ndimage, optimize

from fringes import P_GUESS, S2_FIT_MAX, SMOOTH_PX, THETA_DEG, cos_model, cosine_fit, load_red
from calibration import omega_of
from sagnac import regress, setting_means

HERE = Path(__file__).resolve().parent
MATCH_S = (2e4, 1.8e5)          # px^2, template-matching comparison range
SHIFT_RANGE = (-2e4, 7.8e4)     # px^2, allowed shifts: s + D stays inside the profile
SELF_S = (1e3, 1.1e5)           # px^2, comparison range for the template's own period
AIA_TOL = 1e-10                 # rad, convergence of the frame phases
MATCH_SCAN = [(2e4, 1.8e5), (1e4, 1.2e5), (2e4, 1.4e5), (5e4, 1.8e5)]   # px^2, template robustness
AIA_SCAN = [2.6e5, 1.8e5, 1.3e5, 0.9e5]                                # px^2, largest s used by AIA
FOURIER_GRID = (256, 256)       # (theta, s) samples of the unwrapped image
FOURIER_WIDTH = 0.26            # sideband sigma / carrier frequency (the friend's 0.030 / 0.114 rad/px)
FOURIER_SCAN = [0.15, 0.26, 0.40]


def wrap(x):
    """Angle(s) wrapped into (-pi, pi]."""
    return np.angle(np.exp(1j * np.asarray(x)))


# ---------------------------------------------------------------- 1. template matching

def match_residual(D, s, prof, T):
    """How badly the template, shifted by D, fits the profile: the sum of squared differences
    between prof(s) and g T(s + D) + c, where the gain g and offset c are solved exactly by
    linear least squares for this D. Returns that sum and (g, c)."""
    A = np.c_[T(s + D), np.ones_like(s)]               # columns: shifted template, constant
    gc, *_ = np.linalg.lstsq(A, prof, rcond=None)
    return float(np.sum((prof - A @ gc) ** 2)), gc


def best_shift(s, prof, T, lo, hi):
    """The shift D in [lo, hi] (px^2) that fits best: first a coarse grid every 1000 px^2, then a
    1D minimiser between the grid points either side of the best one (sub-grid precision)."""
    grid = np.arange(lo, hi + 1, 1000.0)
    i = int(np.argmin([match_residual(D, s, prof, T)[0] for D in grid]))
    a, b = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    return optimize.minimize_scalar(lambda D: match_residual(D, s, prof, T)[0],
                                    bounds=(a, b), method='bounded').x


def template_phases(s, M, rest, match=MATCH_S):
    """Template matching for every profile. The template T is the mean of the rest profiles,
    made a smooth function of s by a cubic spline (so it can be shifted by any D). Its own period
    P_T is found by matching T against itself one fringe further out (the shift near P_GUESS that
    maps T onto T). Each profile is then compared with the shifted template over `match` (px^2).
    Returns T, P_T, the shifts D and the phases 2 pi D / P_T (relative to rest)."""
    T = interpolate.CubicSpline(s, M[rest].mean(axis=0))
    sel = (s >= SELF_S[0]) & (s <= SELF_S[1])
    P_T = best_shift(s[sel], T(s[sel]), T, 0.8 * P_GUESS, 1.25 * P_GUESS)
    sel = (s >= match[0]) & (s <= match[1])
    D = np.array([best_shift(s[sel], p[sel], T, *SHIFT_RANGE) for p in M])
    return T, P_T, D, 2 * np.pi * D / P_T


# ---------------------------------------------------------------- 2. AIA (phase shifting)

def aia(M, gain=False):
    """Advanced iterative algorithm on the profile matrix M (rows = photographs k, columns = r^2
    bins j). Model: M_kj = a_j + c_j cos(delta_k) + d_j sin(delta_k), i.e. each bin oscillates as
    a_j + b_j cos(delta_k - psi_j) with c = b cos psi, d = b sin psi (psi = atan2(d, c)).
    Alternate two linear least-squares problems until the frame phases delta stop changing:
      (i)  delta fixed  -> a, c, d for every bin (one regression per bin, all bins at once);
      (ii) a, c, d fixed -> cos delta_k, sin delta_k for every frame -> delta_k = atan2(sin, cos).
    Start: the first two principal components of M (via SVD); for a pure phase sweep they trace
    a circle whose angle is delta. gain=True also fits a brightness factor g_k per frame:
    M_kj = g_k (a_j + c_j cos delta_k + d_j sin delta_k).
    Returns delta, a, c, d and the number of iterations."""
    U, S, _ = np.linalg.svd(M - M.mean(axis=0), full_matrices=False)
    delta = np.arctan2(U[:, 1] * S[1], U[:, 0] * S[0])                 # PCA start
    g = np.ones(len(M))
    for it in range(1, 3001):                                           # at most 3000 rounds
        # (i) delta known: design matrix with columns 1, cos delta, sin delta (times the gain);
        #     one least-squares solve gives a, c, d for all bins at once (columns of M)
        D = g[:, None] * np.c_[np.ones_like(delta), np.cos(delta), np.sin(delta)]
        a, c, d = np.linalg.lstsq(D, M, rcond=None)[0]                  # per bin
        # (ii) a, c, d known: each frame's profile is a combination of the vectors c and d (plus a);
        #      its coefficients u = cos delta_k, v = sin delta_k (times g_k with gain=True)
        if gain:                                                        # per frame: g, g cos, g sin
            gk, u, v = np.linalg.lstsq(np.c_[a, c, d], M.T, rcond=None)[0]
            g = gk / gk.mean()
        else:                                                           # per frame: cos, sin
            u, v = np.linalg.lstsq(np.c_[c, d], (M - a).T, rcond=None)[0]
        new = np.arctan2(v, u)
        done = np.max(np.abs(wrap(new - delta))) < AIA_TOL             # no frame phase moved
        delta = new
        if done:
            break
    return delta, a, c, d, it


def aia_phases(fr, M, gain=False):
    """AIA frame phases relative to the rest frames, unwrapped along increasing S, sign so that
    they grow with S; also delta, a, c, d (same sign convention) and the iteration count.
    (AIA cannot tell delta from -delta, so the sign is fixed by requiring a positive trend.)"""
    delta, a, c, d, n_it = aia(M, gain)
    du = relative_to_rest(fr, delta)
    if np.polyfit(fr.S, du, 1)[0] < 0:
        delta, d, du = -delta, -d, -du
    return du, delta, a, c, d, n_it


def relative_to_rest(fr, phase):
    """Make wrapped frame phases continuous and relative to rest. The frames are put in order of
    S (then time) and np.unwrap removes the 2 pi jumps between neighbours (it assumes successive
    frames differ by less than pi); then the mean of the rest frames is subtracted."""
    order = np.lexsort((fr.time_s, fr.S))              # sort by S, ties by time
    du = np.empty_like(phase); du[order] = np.unwrap(phase[order])
    return du - du[(fr.S == 0).to_numpy()].mean()


# ---------------------------------------------------------------- 3. 2D Fourier sideband

def unwrap_image(img, cx, cy):
    """The image resampled on a (theta, s = r^2) grid around the ring centre, over the fringe
    sector: FOURIER_GRID[0] angles by FOURIER_GRID[1] values of s. Rings (constant s) become
    straight vertical stripes, evenly spaced, which is what a Fourier method needs."""
    th = np.radians(np.linspace(*THETA_DEG, FOURIER_GRID[0], endpoint=False))
    r = np.sqrt(np.linspace(0, S2_FIT_MAX, FOURIER_GRID[1], endpoint=False))
    return ndimage.map_coordinates(img, [cy - np.outer(np.sin(th), r), cx + np.outer(np.cos(th), r)],
                                   order=1, mode='nearest')


def analytic(m, G):
    """The complex fringe signal of a (theta, s) map: subtract the slowly varying background
    (a wide Gaussian blur along s, sigma = 0.44 fringe periods, as the friend's 24 px for 55 px
    fringes), Fourier transform, keep only the region around the positive fringe frequency
    (multiply by the filter G), transform back. Its angle at each point is the local fringe phase."""
    ds = S2_FIT_MAX / m.shape[-1]                      # px^2 per sample along s
    hp = m - ndimage.gaussian_filter1d(m, 0.44 * P_GUESS / ds, axis=-1, mode='nearest')
    return np.fft.ifft2(np.fft.fft2(hp) * G)


def fourier_phases(maps, rest, width=FOURIER_WIDTH):
    """2D Fourier-sideband phase (Takeda) of every (theta, s) map relative to the mean rest map.
    The carrier is the strongest peak of the rest map's spectrum on the positive-s side near the
    expected fringe frequency; the sideband filter is a Gaussian of sigma = width x carrier. The
    phase difference is arg sum Z_k conj(Z_rest), which lines up the local phase before summing.
    Returns the phases (wrapped), the filter G and Z_rest.

    Why the product with conj(Z_rest): Z_k conj(Z_rest) has, at every point, the angle
    (local phase now) - (local phase at rest), i.e. the local shift. Summing these complex
    numbers averages the shift, weighted by fringe strength, without the 2 pi ambiguity that
    averaging raw angles would have."""
    ds = S2_FIT_MAX / maps.shape[-1]
    rest_map = maps[rest].mean(axis=0)
    ref = rest_map - ndimage.gaussian_filter1d(rest_map, 0.44 * P_GUESS / ds, axis=-1, mode='nearest')
    ft, fs = np.meshgrid(np.fft.fftfreq(maps.shape[1]), np.fft.fftfreq(maps.shape[2]), indexing='ij')
    # spectrum of the rest map, with a Hann window against edge artefacts
    spec = np.abs(np.fft.fft2(ref * np.outer(np.hanning(ref.shape[0]), np.hanning(ref.shape[1]))))
    f_guess = ds / P_GUESS                                       # cycles per sample along s
    # the carrier: strongest spectral peak with positive s-frequency within a factor 2 of f_guess
    cand = (fs > 0) & (np.hypot(ft, fs) > 0.5 * f_guess) & (np.hypot(ft, fs) < 2 * f_guess)
    i = np.unravel_index(np.argmax(np.where(cand, spec, 0)), spec.shape)
    # Gaussian sideband filter centred on the carrier
    G = np.exp(-((ft - ft[i]) ** 2 + (fs - fs[i]) ** 2) / (2 * (width * np.hypot(ft[i], fs[i])) ** 2))
    Z_ref = analytic(rest_map, G)
    phases = np.array([np.angle(np.sum(analytic(m, G) * np.conj(Z_ref))) for m in maps])
    return phases, G, Z_ref


# ---------------------------------------------------------------- robustness scans

def slope(fr, phi, model):
    """Sagnac slope, its regression error and s (the rms distance of the setting means from the
    line) for per-frame phases phi, as sagnac.py fits them: the same per-setting averaging and the
    same equal-weight regression, omega from the calibration line (model: its coefficients). Used
    for the robustness scans."""
    tab = setting_means(fr.assign(x=phi), 'x')
    S = tab.S.to_numpy()
    r = regress(np.where(S > 0, omega_of(S, model), 0.0), tab.dphi.to_numpy())
    return r['a'], r['u_a'], r['s']


def illustration(fr, s, M, T, D, phi_aia, delta, a, c, d, maps, G, Z_ref, path):
    """figure*: (a) cosine fit, (b) template match, (c) AIA bins against the frame phase, (d) psi(s),
    (e) the (theta, r^2) map of the same frame, (f) its local Fourier phase about the mean."""
    import matplotlib.pyplot as plt
    import report_style as rs
    rs.apply()
    tex = plt.rcParams['text.usetex']
    w = rs.TEXT_WIDTH
    fig, axs = plt.subplots(3, 2, figsize=(w, 0.95 * w), gridspec_kw=dict(hspace=0.42, wspace=0.22))
    top = np.flatnonzero(fr.S == fr.S.max())
    k = top[np.argsort(D[top])[len(top) // 2]]             # median-shift frame at the top speed
    x = s / 1e5
    xlab = r'$r^2$ ($10^5$ px$^2$)'

    # (a) cosine fit of that frame
    ax = axs[0, 0]
    row = fr.iloc[k]
    img = ndimage.gaussian_filter(load_red(HERE / row.file), SMOOTH_PX)
    popt, _, prof = cosine_fit(img, row.cx, row.cy)
    ax.plot(x, prof, 'o', color='k', ms=1.6, mew=0, label='profile')
    ax.plot(x, cos_model(s, *popt), '-', color=rs.COLORS[2], label='cosine fit')
    ax.set_xlabel(xlab); ax.set_ylabel('intensity (counts)'); ax.legend(loc='upper right')

    # (b) template matching: rest template, the frame, the shifted template
    ax = axs[0, 1]
    sel = (s >= MATCH_S[0]) & (s <= MATCH_S[1])
    _, (g, c0) = match_residual(D[k], s[sel], M[k][sel], T)
    ax.plot(x, T(s), '--', color='0.55', label='rest template')
    ax.plot(x, M[k], 'o', color='k', ms=1.6, mew=0, label='profile')
    ax.plot(x[sel], g * T(s[sel] + D[k]) + c0, '-', color=rs.COLORS[3],
            label=rf'template shifted by {D[k] / 1e3:.0f}k px$^2$' if tex
            else f'template shifted by {D[k] / 1e3:.0f}k px²')
    ax.axvspan(*(np.array(MATCH_S) / 1e5), color='0.92', zorder=0, lw=0)
    ax.set_ylim(0, 1.2 * T(s).max())
    ax.set_xlabel(xlab); ax.legend(loc='upper right')

    # (c) AIA: three bins oscillate with the frame phase (relative to rest)
    ax = axs[1, 0]
    off = np.angle(np.mean(np.exp(1j * (delta - phi_aia))))
    xx = np.linspace(phi_aia.min() - 0.2, phi_aia.max() + 0.2, 200)
    js = np.searchsorted(s, [3e4, 9e4, 1.5e5])
    for j, col, mk in zip(js, (rs.COLORS[0], rs.COLORS[4], rs.COLORS[5]), ('o', 's', '^')):
        ax.plot(phi_aia, M[:, j], mk, color=col, ms=2, mew=0, alpha=0.8,
                label=rf'$r^2 = {s[j] / 1e5:.1f}\times10^5$ px$^2$' if tex else f'r² = {s[j] / 1e5:.1f}e5 px²')
        ax.plot(xx, a[j] + c[j] * np.cos(xx + off) + d[j] * np.sin(xx + off), '-', color=col, lw=0.7)
    ax.set_xlabel('AIA frame phase, relative to rest (rad)'); ax.set_ylabel('intensity (counts)')
    ax.set_ylim(0, 1.1 * M[:, js].max())
    ax.legend(loc='upper right', fontsize=7)

    # (d) AIA: the fitted phase of every bin against r^2
    ax = axs[1, 1]
    psi = np.unwrap(np.arctan2(d, c))
    b = np.hypot(c, d)
    good = b > 0.2 * b.max()
    k1, k0 = np.polyfit(s[good], psi[good], 1, w=b[good])
    ax.plot(x[good], psi[good], 'o', color=rs.COLORS[3], ms=2, mew=0,
            label=r'$\psi$ from AIA' if tex else 'psi from AIA')
    ax.plot(x, k1 * s + k0, '-', color='0.4', lw=0.7,
            label=rf'linear, $P = {2 * np.pi / abs(k1) / 1e3:.1f}$k px$^2$' if tex
            else f'linear, P = {2 * np.pi / abs(k1) / 1e3:.1f}k px²')
    ax.set_xlabel(xlab); ax.set_ylabel(r'bin phase $\psi$ (rad)'); ax.legend(loc='upper right')

    # (e) the frame on (theta, r^2): the rings become straight fringes
    ext = (0, S2_FIT_MAX / 1e5, *THETA_DEG)
    ax = axs[2, 0]
    ax.imshow(maps[k], cmap='gray', vmin=0, vmax=220, origin='lower', aspect='auto', extent=ext)
    ax.set_xlabel(xlab); ax.set_ylabel(r'$\theta$ (deg)')
    ax.minorticks_off()

    # (f) local phase of Z_k conj(Z_rest) about its mean, where the fringes are strong
    ax = axs[2, 1]
    prod = analytic(maps[k], G) * np.conj(Z_ref)
    local = np.angle(prod * np.exp(-1j * np.angle(prod.sum())))
    local[np.abs(prod) < 0.1 * np.abs(prod).max()] = np.nan
    v = float(f'{np.nanpercentile(np.abs(local), 99):.1g}')            # symmetric, one figure
    im = ax.imshow(local, cmap='RdBu_r', vmin=-v, vmax=v, origin='lower', aspect='auto', extent=ext)
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.06)
    cb.set_label('local phase $-$ mean (rad)' if tex else 'local phase - mean (rad)')
    cb.ax.minorticks_off()
    ax.set_xlabel(xlab); ax.set_ylabel(r'$\theta$ (deg)')
    ax.minorticks_off()
    for ax, lab in zip(axs.flat, 'abcdef'):
        ax.text(-0.1, 1.02, f'({lab})', transform=ax.transAxes, va='bottom', ha='left')
    rs.save(fig, path)
    plt.close(fig)
    return 2 * np.pi / abs(k1), np.max(np.abs(psi[good] - (k1 * s[good] + k0)))


def main():
    # 1. Inputs: frame table and r^2 profiles from fringes.py (same order of photographs), the
    #    rest frames (S = 0) and the calibration model (for the slopes in the scans).
    fr = pd.read_csv(HERE / 'fringes_frames.csv')
    z = np.load(HERE / 'fringes_profiles.npz')
    s, M = z['s'], z['profiles']                       # s: bin centres (px^2); M: photographs x bins
    assert list(z['files']) == list(fr.file), 'profiles and frame table differ: rerun fringes.py'
    rest = (fr.S == 0).to_numpy()
    model = json.loads((HERE / 'calibration.json').read_text())['line']['coef']

    # 2-4. The three methods. The Fourier method needs the images themselves, resampled on
    #      (theta, r^2) about each photograph's own ring centre from fringes.py.
    T, P_T, D, phi_t = template_phases(s, M, rest)
    phi_aia, delta, a, c, d, n_it = aia_phases(fr, M)
    maps = np.array([unwrap_image(ndimage.gaussian_filter(load_red(HERE / r.file), SMOOTH_PX), r.cx, r.cy)
                     for r in fr.itertuples()], dtype=np.float32)
    ph_f, G, Z_ref = fourier_phases(maps, rest)
    phi_f = relative_to_rest(fr, ph_f)

    # 5. Per-photograph phases for sagnac.py, and the illustration figure.
    pd.DataFrame(dict(file=fr.file, phi_template_rad=np.round(phi_t, 5),
                      shift_template_px2=np.round(D, 1), phi_aia_rad=np.round(phi_aia, 5),
                      phi_fourier_rad=np.round(phi_f, 5))
                 ).to_csv(HERE / 'methods_frames.csv', index=False)
    P_aia, dev = illustration(fr, s, M, T, D, phi_aia, delta, a, c, d, maps, G, Z_ref,
                              HERE / 'figures' / 'methods_illustration')

    # 6. Robustness: the same slope from other comparison ranges / windows / an AIA brightness
    #    factor / other sideband widths. A method whose slope moves a lot with such choices is
    #    not reliable here (AIA); one that does not (template, Fourier) can be trusted.
    t_scan = [(m, *slope(fr, template_phases(s, M, rest, m)[3], model)) for m in MATCH_SCAN]
    a_scan = [(smax, gain, *slope(fr, aia_phases(fr, M[:, s < smax], gain)[0], model))
              for gain in (False, True) for smax in AIA_SCAN]
    f_scan = [(wd, *slope(fr, relative_to_rest(fr, fourier_phases(maps, rest, wd)[0]), model))
              for wd in FOURIER_SCAN]
    (HERE / 'methods.json').write_text(json.dumps(dict(
        template_scan=[dict(range=m, slope=k, u_fit=u, s=x) for m, k, u, x in t_scan],
        aia_scan=[dict(s_max=m, gain=g, slope=k, u_fit=u, s=x) for m, g, k, u, x in a_scan],
        fourier_scan=[dict(width=wd, slope=k, u_fit=u, s=x) for wd, k, u, x in f_scan]),
        indent=1, default=float))

    # Text report (methods.txt); resid = how far the AIA model is from the profiles, in counts
    resid = M - (a + np.outer(np.cos(delta), c) + np.outer(np.sin(delta), d))
    L = ['Three more phase methods (methods.py)\n',
         f'Template matching: template = mean of the {rest.sum()} rest profiles; its own period '
         f'P_T = {P_T / 1e3:.2f}k px^2 (minima median {fr.P_px2.median() / 1e3:.2f}k); shifts D = '
         f'{D.min() / 1e3:.1f} to {D.max() / 1e3:.1f}k px^2 over {MATCH_S[0]:.0f}-{MATCH_S[1]:.0f} px^2.',
         '  slope against the comparison range (regression error, s):']
    L += [f'    {m[0]:7.0f}-{m[1]:6.0f} px^2: {k:.3f} +- {u:.3f} rad/(rad/s), s = {x:.3f} rad'
          for m, k, u, x in t_scan]
    L += [f'AIA: {len(fr)} frames x {len(s)} bins, converged in {n_it} iterations; rms residual '
          f'{np.sqrt(np.mean(resid ** 2)):.2f} counts; frame phases span {np.ptp(phi_aia):.2f} rad '
          f'({np.degrees(np.ptp(phi_aia)):.0f} deg).',
          f'  bin phases psi(s): linear fit P = {P_aia / 1e3:.1f}k px^2, largest deviation {dev:.2f} rad '
          '(bins with > 20 % of the peak modulation).',
          '  slope against the largest s used, without and with a brightness factor per frame:']
    L += [f'    s < {m:7.0f} px^2, {"gain " if g else "plain"}: {k:.3f} +- {u:.3f} rad/(rad/s), '
          f's = {x:.3f} rad' for m, g, k, u, x in a_scan]
    L += ['  AIA assumes that only the phase changes between frames. On synthetic profiles '
          '(tests/test_methods.py) it recovers a 3.3 rad sweep to 0.01 rad; a 5 % brightness drift '
          'leaves the slope right but puts single frames off by up to 0.1 rad, and a 1.7 % '
          'frame-to-frame change of the fringe period, as measured, gives errors of order 0.2 rad. '
          'Its spread over windows and variants here measures how far the photographs break that '
          'assumption.',
          f'2D Fourier sideband on (theta, r^2) maps of {FOURIER_GRID[0]} x {FOURIER_GRID[1]} samples '
          f'({THETA_DEG[0]:.0f} to {THETA_DEG[1]:.0f} deg, r^2 < {S2_FIT_MAX:.0f} px^2) about each '
          f'frame\'s ring centre; sideband sigma = {FOURIER_WIDTH} x carrier; rest frames scatter by '
          f'{phi_f[rest].std(ddof=1):.3f} rad.',
          '  slope against the sideband width (regression error, s):']
    L += [f'    width {wd:.2f} x carrier: {k:.3f} +- {u:.3f} rad/(rad/s), s = {x:.3f} rad'
          for wd, k, u, x in f_scan]
    L += ['Per-frame phases in methods_frames.csv; statistics and the comparison with the other '
          'methods in sagnac.txt.']
    (HERE / 'methods.txt').write_text('\n'.join(L) + '\n', encoding='utf-8')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
