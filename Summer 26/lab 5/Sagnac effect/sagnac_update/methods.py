"""Two more ways to read the fringe phase, as cross-checks of fringes.py (minima, cosine fit).

WHAT THIS SCRIPT IS FOR
fringes.py reads the phase from the positions of the dark rings (and from a cosine fit). If the
Sagnac slope is real, other ways of reading the same photographs must give about the same
slope. This script implements two more and writes their phase for every photograph to
methods_frames.csv; sagnac.py fits all methods in the same way and compares them.
Input: fringes_frames.csv and fringes_profiles.npz (fringes.py): for every photograph, the ring
centre and the intensity averaged over the fringe sector in bins of s = r^2 about that centre;
the random-points method also reads the photographs themselves.

HOW IT RUNS (main(), top to bottom)
  1. Load the frame table, the r^2 profiles M (one row per photograph) and the calibration.
  2. template_phases   template matching (below)                          -> phi_t
  3. point_phases      ten random pixels, fixed in the image (below)      -> phi_p
  4. Save the phases per photograph (methods_frames.csv) and draw the four-panel figure
     methods_illustration (one example of each method on real data).
  5. Robustness scans: rerun each method with other choices (template: comparison range;
     random points: other random draws, and a diagnostic with per-photograph geometry) and
     record the resulting slopes (methods.json, methods.txt).

THE TWO METHODS
  Template matching. The rest frames (S = 0) are averaged into a template T(s). Each profile
  is matched to g T(s + D) + c by least squares over the shift D, with gain g and offset c free.
  In s, an exponential envelope shifted by D only changes by a constant factor, which g absorbs.
  Rings that move inwards give D > 0, and the phase relative to rest is 2 pi D / P_T, where P_T
  is the shift that maps the template onto itself one fringe further out. No fringe model.
  In words: slide the at-rest profile along r^2 until it lies on top of the new one.

  Random fixed points. POINTS_N = 10 pixels are drawn at random (seed POINTS_SEED) inside the
  fringe region and kept at the same image positions in every photograph. As the rings move,
  each pixel brightens and dims:
      I_kj = a_j + b_j cos(psi_j + delta_k),
  pixel j, photograph k. psi_j = 2 pi s_j / P0 is the pixel's local fringe phase at rest (s_j =
  its r^2 from the rest ring centre, P0 the rest fringe period); delta_k is the photograph's
  phase, the unknown; a_j and b_j are the pixel's unknown offset and amplitude. One alternating
  least-squares fit finds all of them (point_fit). Two physical limits keep it well posed:
  b_j >= 0 and a_j - b_j >= the camera's dark level (a pixel's darkest fringe cannot be darker
  than the camera's black). Without them, ten pixels that see only half a fringe cannot tell a
  large phase change at low contrast from a small one at high contrast.
  In words: watch ten fixed pixels get brighter and darker, and find the one phase change per
  photograph that explains all ten at once.

Template matching reads an average phase over its comparison window, and a fixed pixel reads
the local phase at its own radius, not the phase at the ring centre. A change of the fringe
period between photographs (fringes.py), and for the pixels also a move of the ring centre,
therefore enters both, but not the minima or the cosine fit.

Outputs: methods_frames.csv (phi_template_rad, phi_points_rad per frame), methods.txt,
methods.json (the chosen pixels and the robustness scans), figures/methods_illustration.pdf/.png
(figure*). Run from this folder after calibration.py and fringes.py, before sagnac.py
(a few minutes).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import interpolate, ndimage, optimize

from fringes import P_GUESS, S2_FIT_MAX, SMOOTH_PX, THETA_DEG, cos_model, cosine_fit, load_red
from sagnac import calibration_part, fit_line, setting_means

HERE = Path(__file__).resolve().parent
MATCH_S = (2e4, 1.8e5)          # px^2, template-matching comparison range
SHIFT_RANGE = (-2e4, 7.8e4)     # px^2, allowed shifts: s + D stays inside the profile
SELF_S = (1e3, 1.1e5)           # px^2, comparison range for the template's own period
MATCH_SCAN = [(2e4, 1.8e5), (1e4, 1.2e5), (2e4, 1.4e5), (5e4, 1.8e5)]   # px^2, template robustness
POINTS_N = 10                   # number of fixed pixels
POINTS_SEED = 20261011          # fixed in advance (the measurement date), not chosen by its result
POINTS_SCAN_SEEDS = range(1, 31)   # 30 more draws for the robustness scan
POINTS_TOL = 1e-6               # rad, convergence of the photograph phases
POINTS_MAX_ITER = 500
PHASE_GRID = np.linspace(-np.pi, np.pi, 721)    # half-degree grid for the phase search (refined)


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

# ---------------------------------------------------------------- 2. random fixed points

def fringe_region(shape, cx, cy):
    """(y, x) of every pixel inside the fringe sector THETA_DEG with r^2 < S2_FIT_MAX about (cx, cy):
    the part of the image where rings are visible (the lower left of the pattern has none)."""
    y, x = np.indices(shape)
    s = (x - cx) ** 2 + (y - cy) ** 2
    th = np.degrees(np.arctan2(-(y - cy), x - cx))      # minus: image y grows downwards
    sel = (th > THETA_DEG[0]) & (th < THETA_DEG[1]) & (s < S2_FIT_MAX)
    return y[sel], x[sel]


def draw_points(ys, xs, seed, n=POINTS_N):
    """n different pixels of the region, drawn uniformly (every pixel equally likely, so equal areas
    get equal numbers of points). The same seed always gives the same pixels."""
    i = np.random.default_rng(seed).choice(len(ys), n, replace=False)
    return ys[i], xs[i]


def dark_level(img):
    """The camera's dark level (counts): the median of the four 100 x 100 px image corners, which
    lie far outside the beam. No pixel of the fringe pattern can be darker than this."""
    c = np.concatenate([img[:100, :100].ravel(), img[:100, -100:].ravel(),
                        img[-100:, :100].ravel(), img[-100:, -100:].ravel()])
    return float(np.median(c))


def best_phases(M, a, b, psi):
    """For every photograph k, the phase delta_k that best explains its pixel intensities M_kj with
    the pixels' offsets a_j and amplitudes b_j known: minimise sum_j (M_kj - a_j - b_j cos(psi_j +
    delta))^2. psi is (pixels,) or (photographs, pixels). A grid search over PHASE_GRID (all
    photographs at once), refined by a parabola through the best grid point and its neighbours."""
    P = np.broadcast_to(psi, M.shape)
    r = ((M[:, None, :] - a - b * np.cos(P[:, None, :] + PHASE_GRID[None, :, None])) ** 2).sum(axis=2)
    i = np.argmin(r, axis=1)                            # best grid phase of each photograph
    k = np.arange(len(M))
    y0, y1, y2 = r[k, (i - 1) % len(PHASE_GRID)], r[k, i], r[k, (i + 1) % len(PHASE_GRID)]
    den = y0 - 2 * y1 + y2
    off = np.where(den > 0, 0.5 * (y0 - y2) / np.where(den > 0, den, 1.0), 0.0)
    return PHASE_GRID[i] + off * (PHASE_GRID[1] - PHASE_GRID[0])


def offsets_amplitudes(M, X, dark):
    """For every pixel j, offset a_j and amplitude b_j minimising sum_k (M_kj - a_j - b_j X_kj)^2,
    X = cos(psi_j + delta_k), under the physical limits b_j >= 0 and a_j - b_j >= dark (the darkest
    point of the pixel's fringe cannot be darker than the camera's dark level).
    Written with c = a - b >= dark, the model is M = c + b (1 + X): a two-parameter least squares
    with bounds on both parameters. Its optimum is the unconstrained solution if that is allowed,
    otherwise the best point on one of the two limits; all candidates are compared, per pixel."""
    Z = 1 + X                                           # photographs x pixels
    n = len(M)
    Sz, Szz, Sy, Szy = Z.sum(0), (Z * Z).sum(0), M.sum(0), (Z * M).sum(0)
    det = n * Szz - Sz ** 2
    b_u = np.where(det > 0, (n * Szy - Sz * Sy) / np.where(det > 0, det, 1.0), 0.0)
    c_u = (Sy - b_u * Sz) / n
    cands = [(c_u, b_u, (c_u >= dark) & (b_u >= 0)),                                  # interior
             (np.full_like(Sy, dark), np.maximum(0.0, (Szy - dark * Sz) / Szz), True),   # on c = dark
             (np.maximum(dark, Sy / n), np.zeros_like(Sy), True)]                     # on b = 0
    best_c, best_b, best_sse = None, None, None
    for c, b, ok in cands:
        sse = ((M - c - b * Z) ** 2).sum(0)
        sse = np.where(ok, sse, np.inf)
        if best_sse is None:
            best_c, best_b, best_sse = c, b, sse
        else:
            take = sse < best_sse
            best_c, best_b, best_sse = np.where(take, c, best_c), np.where(take, b, best_b), np.minimum(sse, best_sse)
    return best_c + best_b, best_b


def point_fit(M, psi, dark):
    """The random-points phase. M: intensities, photographs x pixels; psi: each pixel's local
    fringe phase 2 pi s_j / P (from the rest geometry, or per photograph). Model:
        M_kj = a_j + b_j cos(psi_j + delta_k)
    with one unknown phase delta_k per photograph (the same convention as fringes.ring_phase: the
    phase at the ring centre if psi is measured from it) and an unknown offset a_j and amplitude
    b_j per pixel. Alternate until no delta moves by more than POINTS_TOL:
      (i)  a, b known -> delta for every photograph (best_phases);
      (ii) delta known -> a, b for every pixel (offsets_amplitudes, with the physical limits).
    Start: a = the pixel's mean, b = half its range. Returns delta (wrapped), a, b, the rms residual
    (counts) and the number of rounds (POINTS_MAX_ITER means it did not converge)."""
    a, b = M.mean(axis=0), 0.5 * np.ptp(M, axis=0)
    delta = best_phases(M, a, b, psi)
    for it in range(1, POINTS_MAX_ITER + 1):
        a, b = offsets_amplitudes(M, np.cos(np.broadcast_to(psi, M.shape) + delta[:, None]), dark)
        new = best_phases(M, a, b, psi)
        done = np.max(np.abs(wrap(new - delta))) < POINTS_TOL
        delta = new
        if done:
            break
    resid = M - (a + b * np.cos(np.broadcast_to(psi, M.shape) + delta[:, None]))
    return wrap(delta), a, b, float(np.sqrt(np.mean(resid ** 2))), it


def point_phases(fr, M, psi, dark):
    """point_fit, then the phases made continuous and relative to rest (relative_to_rest). Unlike a
    model-free fit, the sign needs no fixing: psi ties the phase to the ring geometry."""
    delta, a, b, rms, it = point_fit(M, psi, dark)
    return relative_to_rest(fr, delta), delta, a, b, rms, it


def relative_to_rest(fr, phase):
    """Make wrapped frame phases continuous and relative to rest. The frames are put in order of
    S (then time) and np.unwrap removes the 2 pi jumps between neighbours (it assumes successive
    frames differ by less than pi); then the mean of the rest frames is subtracted."""
    order = np.lexsort((fr.time_s, fr.S))              # sort by S, ties by time
    du = np.empty_like(phase); du[order] = np.unwrap(phase[order])
    return du - du[(fr.S == 0).to_numpy()].mean()


# ---------------------------------------------------------------- robustness scans

def slope(fr, phi, model):
    """Sagnac slope, its fit error and sigma_x for per-frame phases phi, as sagnac.py fits them
    (same per-setting averaging and the same weighted line). Used for the robustness scans."""
    tab = setting_means(fr.assign(x=phi), 'x')
    S, y, sem = tab.S.to_numpy(), tab.phi.to_numpy(), tab['sem'].to_numpy()
    omega = calibration_part(S, y, None, model)[0]
    a, _, cov, _, sx, *_ = fit_line(omega, y, sem, model['eps'])
    return a, float(np.sqrt(cov[0, 0])), sx
def illustration(fr, s, M, T, D, pts, V, phi_p, delta, a, b, psi, rest_img, centre, path):
    """figure*: (a) cosine fit and (b) template match of one top-speed photograph, (c) the ten fixed
    pixels on the rest pattern, (d) their intensities against the fitted phase, with the model."""
    import matplotlib.pyplot as plt
    import report_style as rs
    rs.apply()
    tex = plt.rcParams['text.usetex']
    w = rs.TEXT_WIDTH
    fig, axs = plt.subplots(2, 2, figsize=(w, 0.66 * w), gridspec_kw=dict(hspace=0.36, wspace=0.24))
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
    # (c) the rest pattern with the ten fixed pixels, numbered
    ax = axs[1, 0]
    cx, cy = centre
    half = 540                                               # the region reaches r = 510 px
    x0, y0 = max(int(cx) - half, 0), max(int(cy) - half, 0)
    crop = rest_img[y0:int(cy) + half, x0:int(cx) + half]
    ax.set_facecolor('black')                                # beyond the image edge
    ax.imshow(crop, cmap='gray', vmin=0, vmax=220, extent=(x0, x0 + crop.shape[1], y0 + crop.shape[0], y0))
    py, px = pts
    cols = plt.get_cmap('viridis')(np.linspace(0.2, 1.0, len(px)))
    for j in range(len(px)):
        ax.plot(px[j], py[j], 'o', ms=4, mfc='none', mec=cols[j], mew=1.0)
        ax.text(px[j] + 14, py[j] - 10, str(j + 1), color=cols[j], fontsize=6.5)
    ax.plot(cx, cy, '+', color=rs.COLORS[4], ms=4, mew=0.6)
    ax.set_xlim(cx - half, cx + half); ax.set_ylim(cy + half, cy - half)
    ax.set_xticks([]); ax.set_yticks([]); ax.minorticks_off()

    # (d) intensity at the ten pixels against the fitted phase, and the fitted cosines
    ax = axs[1, 1]
    rest = (fr.S == 0).to_numpy()
    ref = np.angle(np.mean(np.exp(1j * delta[rest])))       # circular mean of the rest phases
    rel = phi_p                                              # unwrapped, relative to rest
    xx = np.linspace(rel.min() - 0.1, rel.max() + 0.1, 200)  # the model only where there are data
    for j in range(V.shape[1]):
        ax.plot(rel, V[:, j], 'o', color=cols[j], ms=1.6, mew=0, alpha=0.8)
        ax.plot(xx, a[j] + b[j] * np.cos(psi[j] + ref + xx), '-', color=cols[j], lw=0.7)
    ax.set_xlabel('fitted phase, relative to rest (rad)'); ax.set_ylabel('intensity (counts)')
    ax.set_ylim(0, None)
    for ax, lab in zip(axs.flat, 'abcd'):
        ax.text(-0.1, 1.02, f'({lab})', transform=ax.transAxes, va='bottom', ha='left')
    rs.save(fig, path)
    plt.close(fig)


def main():
    # 1. Inputs: frame table and r^2 profiles from fringes.py (same order of photographs), the
    #    rest frames (S = 0) and the calibration model (for the slopes in the scans).
    fr = pd.read_csv(HERE / 'fringes_frames.csv')
    z = np.load(HERE / 'fringes_profiles.npz')
    s, M = z['s'], z['profiles']                       # s: bin centres (px^2); M: photographs x bins
    assert list(z['files']) == list(fr.file), 'profiles and frame table differ: rerun fringes.py'
    rest = (fr.S == 0).to_numpy()
    model = json.loads((HERE / 'calibration.json').read_text())['conversion']

    # 2. Template matching on the r^2 profiles.
    T, P_T, D, phi_t = template_phases(s, M, rest)

    # 3. Random fixed points. The pixels are drawn once, inside the fringe region around the rest
    #    ring centre, and read at the same image positions in every photograph. Their local phase
    #    psi_j = 2 pi s_j / P0 comes from the rest geometry (centre and period of the rest frames).
    #    All draws (the main one and the scan) are read in one pass over the photographs, from the
    #    same smoothed images the other methods use.
    cx0, cy0, P0 = fr.cx[rest].median(), fr.cy[rest].median(), fr.P_px2[rest].median()
    rest_img = np.mean([ndimage.gaussian_filter(load_red(HERE / f), SMOOTH_PX) for f in fr.file[rest]], axis=0)
    dark = dark_level(rest_img)
    ys, xs = fringe_region(rest_img.shape, cx0, cy0)
    seeds = [POINTS_SEED] + list(POINTS_SCAN_SEEDS)
    draws = [draw_points(ys, xs, sd) for sd in seeds]
    Y, X = np.concatenate([d[0] for d in draws]), np.concatenate([d[1] for d in draws])
    V_all = np.array([ndimage.gaussian_filter(load_red(HERE / f), SMOOTH_PX)[Y, X] for f in fr.file])
    n = POINTS_N
    psi_of = lambda d: 2 * np.pi * ((d[1] - cx0) ** 2 + (d[0] - cy0) ** 2) / P0      # rest geometry
    CX, CY, PP = (fr[c].to_numpy()[:, None] for c in ('cx', 'cy', 'P_px2'))
    psi_frame = lambda d: 2 * np.pi * ((d[1][None] - CX) ** 2 + (d[0][None] - CY) ** 2) / PP   # diagnostic
    V = V_all[:, :n]
    phi_p, delta_p, a_p, b_p, rms_p, it_p = point_phases(fr, V, psi_of(draws[0]), dark)

    # 4. Per-photograph phases for sagnac.py, and the illustration figure.
    pd.DataFrame(dict(file=fr.file, phi_template_rad=np.round(phi_t, 5),
                      shift_template_px2=np.round(D, 1), phi_points_rad=np.round(phi_p, 5))
                 ).to_csv(HERE / 'methods_frames.csv', index=False)
    illustration(fr, s, M, T, D, draws[0], V, phi_p, delta_p, a_p, b_p, psi_of(draws[0]), rest_img, (cx0, cy0),
                 HERE / 'figures' / 'methods_illustration')

    # 5. Robustness. Template: other comparison ranges. Random points: the same fit for other
    #    random draws of ten pixels (how much the answer depends on which pixels), and, as a
    #    diagnostic, the same pixels with psi from each photograph's own fitted centre and period
    #    (how much of the result comes from period changes and centre motion, which a fixed pixel
    #    cannot tell apart from a phase shift).
    t_scan = [(m, *slope(fr, template_phases(s, M, rest, m)[3], model)) for m in MATCH_SCAN]
    p_scan, g_scan = [], []
    for i, (sd, d) in enumerate(zip(seeds, draws)):
        Vi = V_all[:, i * n:(i + 1) * n]
        ph, _, _, _, rms, it = point_phases(fr, Vi, psi_of(d), dark)
        p_scan.append((sd, *slope(fr, ph, model), rms, it))
        ph, _, _, _, rms, it = point_phases(fr, Vi, psi_frame(d), dark)
        g_scan.append((sd, *slope(fr, ph, model), rms, it))
    main_draw = dict(seed=POINTS_SEED, n=n, dark=dark, centre=[cx0, cy0], P0=P0,
                     pixels_xy=[[int(a), int(b)] for a, b in zip(draws[0][1], draws[0][0])],
                     offsets=a_p.tolist(), amplitudes=b_p.tolist(), rms=rms_p, iterations=it_p)
    scan_row = lambda r: dict(seed=r[0], slope=r[1], u_fit=r[2], sigma_x=r[3], rms=r[4], iterations=r[5])
    (HERE / 'methods.json').write_text(json.dumps(dict(
        calibration_model=model.get('name', 'polynomial'), calibration_degree=model['degree'],
        template_scan=[dict(range=m, slope=k, u_fit=u, sigma_x=x) for m, k, u, x in t_scan],
        points=main_draw,
        points_scan=[scan_row(r) for r in p_scan],
        points_frame_geometry_scan=[scan_row(r) for r in g_scan]),
        indent=1, default=float))

    # 6. Text report (methods.txt).
    q = lambda v, p: float(np.percentile(v, p))
    ps, gs = np.array([r[1] for r in p_scan]), np.array([r[1] for r in g_scan])
    L = ['Two more phase methods (methods.py)\n',
         f"Calibration: {model.get('name', 'polynomial')}, degree {model['degree']}.\n",
         f'Template matching: template = mean of the {rest.sum()} rest profiles; its own period '
         f'P_T = {P_T / 1e3:.2f}k px^2 (minima median {fr.P_px2.median() / 1e3:.2f}k); shifts D = '
         f'{D.min() / 1e3:.1f} to {D.max() / 1e3:.1f}k px^2 over {MATCH_S[0]:.0f}-{MATCH_S[1]:.0f} px^2.',
         '  slope against the comparison range (fit error, sigma_x):']
    L += [f'    {m[0]:7.0f}-{m[1]:6.0f} px^2: {k:.3f} +- {u:.3f} rad/(rad/s), sigma_x = {x:.3f} rad'
          for m, k, u, x in t_scan]
    L += [f'Random points: {n} pixels drawn with seed {POINTS_SEED} inside the fringe sector '
          f'({THETA_DEG[0]:.0f} to {THETA_DEG[1]:.0f} deg, r^2 < {S2_FIT_MAX:.0f} px^2) about the rest centre '
          f'({cx0:.1f}, {cy0:.1f}) px, read at the same positions in all {len(fr)} photographs (smoothed, '
          f'sigma = {SMOOTH_PX:.0f} px). Local phase from the rest period P0 = {P0 / 1e3:.2f}k px^2; '
          f'dark level {dark:.1f} counts.',
          '  pixels (x, y): ' + ', '.join(f'({p[0]}, {p[1]})' for p in main_draw['pixels_xy']),
          f'  fit: {it_p} rounds{" (not converged)" if it_p >= POINTS_MAX_ITER else ""}, rms residual '
          f'{rms_p:.2f} counts; amplitudes {b_p.min():.1f}-{b_p.max():.1f} counts; phases span '
          f'{np.ptp(phi_p):.2f} rad.',
          f'  slope of this draw: {p_scan[0][1]:.3f} +- {p_scan[0][2]:.3f} rad/(rad/s) (fit), sigma_x = '
          f'{p_scan[0][3]:.3f} rad.',
          f'  {len(seeds)} draws (this one and {len(seeds) - 1} more, seeds {seeds[1]}-{seeds[-1]}): slope median '
          f'{np.median(ps):.3f}, 16-84 % {q(ps, 16):.3f}-{q(ps, 84):.3f}, full range {ps.min():.3f}-{ps.max():.3f}; '
          f'not converged: {sum(r[5] >= POINTS_MAX_ITER for r in p_scan)}.',
          f'  diagnostic, the same pixels with the local phase from each photograph\'s own fitted centre and '
          f'period: median {np.median(gs):.3f}, 16-84 % {q(gs, 16):.3f}-{q(gs, 84):.3f}.',
          '  A fixed pixel reads the local phase at its own radius. A change of the fringe period (the '
          'wavefront curvature) or a move of the ring centre changes that local phase as well, and the '
          'pixel cannot tell either from a Sagnac shift; the diagnostic removes both. With ten pixels and a '
          'sweep of half a fringe, the answer also depends on which pixels are drawn.',
          'Per-frame phases in methods_frames.csv; statistics and the comparison with the other '
          'methods in sagnac.txt.']
    (HERE / 'methods.txt').write_text('\n'.join(L) + '\n', encoding='utf-8')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
