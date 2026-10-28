"""Fixed pixels: a fourth way to read the fringe phase (preparation question 5), as a check of the
ring readout of fringes.py and of methods.py.

Ten pixels at fixed positions in the image are read in every photograph. As the rings move, each
pixel brightens and dims:
    I_kj = a_j + b_j cos(psi_j + delta_k),         pixel j, photograph k,
where psi_j = 2 pi s_j / P0 is the pixel's fringe phase at rest (s_j its r^2 from the rest ring
centre, P0 the rest spacing), delta_k is the photograph's phase, the unknown, and a_j, b_j are the
pixel's offset and amplitude. An alternating least-squares fit finds them all (fit_phases) under
two physical limits: b_j >= 0 and a_j - b_j >= the camera's dark level (a pixel's darkest fringe
cannot be darker than the camera's black). Without them, ten pixels that see only half a fringe
cannot tell a large phase change at low contrast from a small one at high contrast.

A fixed pixel reads the local phase at its own radius: a move of the ring centre or a change of
the spacing between photographs reads as a phase change, and over this run no pixel sweeps a
whole fringe. The answer therefore depends on which pixels are drawn, so the fit is repeated for
many random draws and the spread of the slope between draws is the method's uncertainty.

Reads fringes_frames.csv (fringes.py), calibration.json and the photographs. Writes
pixels_frames.csv (the phase of every photograph for the main draw, which sagnac.py fits like the
other readouts) and pixels.json (the main draw, every draw's slope, and the spread per number of
pixels). Run from this folder after fringes.py:   python pixels.py [draws]
(1000 draws of ten pixels by default, plus 60 draws each of 20, 40 and 80 pixels; a few minutes).
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage, stats

from calibration import omega_of
from fringes import S_MAX, SECTOR_DEG, SMOOTH_PX, load_red

HERE = Path(__file__).resolve().parent
N_PIXELS = 10
MAIN_SEED = 20261011                              # the measurement date, fixed in advance of any result
N_DRAWS = 1000                                    # draws of N_PIXELS pixels (seeds 1..N_DRAWS)
MORE_PIXELS = [(20, 60), (40, 60), (80, 60)]      # (pixels per draw, draws): does the spread fall with more pixels?
PHASE_GRID = np.linspace(-np.pi, np.pi, 361)      # 1 deg grid for the phase search, refined by a parabola
TOL, MAX_ROUNDS = 1e-6, 500                       # rad; rounds of the alternating fit


def wrap(x):
    """Angles wrapped into (-pi, pi]."""
    return np.angle(np.exp(1j * np.asarray(x)))


def fringe_region(shape, cx, cy):
    """(y, x) of every pixel inside the fringe sector SECTOR_DEG with r^2 < S_MAX about (cx, cy)."""
    y, x = np.indices(shape)
    s = (x - cx) ** 2 + (y - cy) ** 2
    angle = np.degrees(np.arctan2(-(y - cy), x - cx))          # minus: image y grows downwards
    inside = (angle > SECTOR_DEG[0]) & (angle < SECTOR_DEG[1]) & (s < S_MAX)
    return y[inside], x[inside]


def draw_pixels(ys, xs, seed, n):
    """n different pixels of the region, every pixel equally likely; the same seed gives the same pixels."""
    i = np.random.default_rng(seed).choice(len(ys), n, replace=False)
    return ys[i], xs[i]


def dark_level(img):
    """The camera's dark level (counts): the median of the four 100 x 100 px corners, far outside the beam."""
    corners = [img[:100, :100], img[:100, -100:], img[-100:, :100], img[-100:, -100:]]
    return float(np.median(np.concatenate([c.ravel() for c in corners])))


def best_phases(I, a, b, psi):
    """The phase delta_k of every photograph with the pixels' a_j, b_j known: the grid point of
    PHASE_GRID with the least sum of squares, refined by a parabola through it and its neighbours."""
    r = ((I[:, None, :] - a - b * np.cos(psi + PHASE_GRID[None, :, None])) ** 2).sum(axis=2)   # photographs x grid
    i = r.argmin(axis=1)
    k = np.arange(len(I))
    y0, y1, y2 = r[k, (i - 1) % len(PHASE_GRID)], r[k, i], r[k, (i + 1) % len(PHASE_GRID)]
    curvature = y0 - 2 * y1 + y2
    offset = np.where(curvature > 0, 0.5 * (y0 - y2) / np.where(curvature > 0, curvature, 1.0), 0.0)
    return PHASE_GRID[i] + offset * (PHASE_GRID[1] - PHASE_GRID[0])


def offsets_amplitudes(I, X, dark):
    """a_j and b_j of every pixel with the phases known (X_kj = cos(psi_j + delta_k)), by least squares
    under b >= 0 and a - b >= dark. With c = a - b the model is I = c + b (1 + X): the best allowed
    point is the free solution if it is allowed, otherwise the best point on one of the two limits."""
    Z = 1 + X
    n = len(I)
    Sz, Szz, Sy, Szy = Z.sum(0), (Z * Z).sum(0), I.sum(0), (Z * I).sum(0)
    det = n * Szz - Sz ** 2
    b_free = np.where(det > 0, (n * Szy - Sz * Sy) / np.where(det > 0, det, 1.0), 0.0)
    c_free = (Sy - b_free * Sz) / n
    candidates = [(c_free, b_free, (c_free >= dark) & (b_free >= 0)),                 # free solution
                  (np.full_like(Sy, dark), np.maximum(0.0, (Szy - dark * Sz) / Szz), True),  # on a - b = dark
                  (np.maximum(dark, Sy / n), np.zeros_like(Sy), True)]                    # on b = 0
    c_best = b_best = sse_best = None
    for c, b, allowed in candidates:
        sse = np.where(allowed, ((I - c - b * Z) ** 2).sum(0), np.inf)
        if sse_best is None:
            c_best, b_best, sse_best = c, b, sse
        else:
            take = sse < sse_best
            c_best, b_best, sse_best = np.where(take, c, c_best), np.where(take, b, b_best), np.minimum(sse, sse_best)
    return c_best + b_best, b_best


def fit_phases(I, psi, dark):
    """Alternating least squares for I_kj = a_j + b_j cos(psi_j + delta_k): phases with a, b held,
    then a, b with the phases held, until no phase moves by more than TOL. Returns the wrapped phases,
    a, b, the rms residual (counts) and the number of rounds (MAX_ROUNDS means it did not converge)."""
    a, b = I.mean(axis=0), 0.5 * np.ptp(I, axis=0)
    delta = best_phases(I, a, b, psi)
    for rounds in range(1, MAX_ROUNDS + 1):
        a, b = offsets_amplitudes(I, np.cos(psi + delta[:, None]), dark)
        new = best_phases(I, a, b, psi)
        converged = np.max(np.abs(wrap(new - delta))) < TOL
        delta = new
        if converged:
            break
    residual = I - (a + b * np.cos(psi + delta[:, None]))
    return wrap(delta), a, b, float(np.sqrt(np.mean(residual ** 2))), rounds


def relative_to_rest(S, time_s, delta):
    """Wrapped photograph phases made continuous (unwrapped in order of S, then time) and relative to
    the mean of the rest photographs."""
    order = np.lexsort((time_s, S))
    continuous = np.empty_like(delta)
    continuous[order] = np.unwrap(delta[order])
    return continuous - continuous[S == 0].mean()


def slope_of(S, omega, phases):
    """The setting means of the phases against omega by an equal-weight line, as sagnac.py fits them."""
    table = pd.DataFrame(dict(S=S, omega=omega, phi=phases)).groupby('S', sort=True).mean()
    fit = stats.linregress(table.omega, table.phi - table.phi.iloc[0])
    return float(fit.slope), float(fit.stderr)


def one_draw(task):
    """Fit one draw and regress it: (slope, fit error, rms, rounds, phases relative to rest)."""
    I, psi, dark, S, time_s, omega = task
    delta, _, _, rms, rounds = fit_phases(I, psi, dark)
    phases = relative_to_rest(S, time_s, delta)
    return (*slope_of(S, omega, phases), rms, rounds, phases)


def main():
    t0 = time.time()
    n_draws = int(sys.argv[1]) if len(sys.argv) > 1 else N_DRAWS
    frames = pd.read_csv(HERE / 'fringes_frames.csv')
    S, time_s = frames.S.to_numpy(), frames.time_s.to_numpy()
    rest = S == 0
    line = json.loads((HERE / 'calibration.json').read_text())['line']
    omega = np.where(rest, 0.0, omega_of(S, line))
    cx0, cy0, P0 = frames.cx[rest].median(), frames.cy[rest].median(), frames.P_px2[rest].median()

    # The draws: the main one, then n_draws of ten pixels, then draws with more pixels.
    shape = load_red(HERE / frames.file[0]).shape
    ys, xs = fringe_region(shape, cx0, cy0)
    plan = [(MAIN_SEED, N_PIXELS)] + [(seed, N_PIXELS) for seed in range(1, n_draws + 1)]
    for n, k in MORE_PIXELS:
        plan += [(100000 * n + seed, n) for seed in range(1, k + 1)]
    draws = [draw_pixels(ys, xs, seed, n) for seed, n in plan]
    Y, X = np.concatenate([d[0] for d in draws]), np.concatenate([d[1] for d in draws])
    edges = np.r_[0, np.cumsum([n for _, n in plan])]

    # One pass over the photographs: every draw's pixel values, and the mean rest image for the dark level.
    values = np.empty((len(frames), len(Y)))
    rest_sum = np.zeros(shape)
    for k, f in enumerate(frames.file):
        img = ndimage.gaussian_filter(load_red(HERE / f), SMOOTH_PX)
        values[k] = img[Y, X]
        if rest[k]:
            rest_sum += img
    dark = dark_level(rest_sum / rest.sum())
    print(f'read {len(frames)} photographs for {len(plan)} draws ({len(Y)} pixels); dark level {dark:.1f} counts; '
          f'{time.time() - t0:.0f} s', flush=True)

    tasks = [(values[:, edges[i]:edges[i + 1]], 2 * np.pi * ((d[1] - cx0) ** 2 + (d[0] - cy0) ** 2) / P0,
              dark, S, time_s, omega) for i, d in enumerate(draws)]
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as pool:
        results = list(pool.map(one_draw, tasks, chunksize=8))
    rows = [dict(seed=int(seed), n=int(n), slope=a, u_fit=u, rms=rms, rounds=rounds)
            for (seed, n), (a, u, rms, rounds, _) in zip(plan, results)]

    # Per number of pixels: the spread of the slope between draws (fits with a negative slope failed).
    summary = []
    for n in sorted({n for _, n in plan}):
        slopes = np.array([r['slope'] for r in rows if r['n'] == n])
        good = slopes[slopes > 0]
        summary.append(dict(n=int(n), draws=len(slopes), failed=int(np.sum(slopes <= 0)), mean=float(good.mean()),
                            sd=float(good.std(ddof=1)), median=float(np.median(good)),
                            median_u_fit=float(np.median([r['u_fit'] for r in rows if r['n'] == n]))))

    main_slope, main_u, main_rms, main_rounds, main_phases = results[0]
    pd.DataFrame(dict(file=frames.file, phi_pixels_rad=np.round(main_phases, 5))).to_csv(
        HERE / 'pixels_frames.csv', index=False)
    (HERE / 'pixels.json').write_text(json.dumps(dict(
        main=dict(seed=MAIN_SEED, n=N_PIXELS, pixels_xy=[[int(x), int(y)] for y, x in zip(*draws[0])], dark=dark,
                  centre=[float(cx0), float(cy0)], P0=float(P0), slope=main_slope, u_fit=main_u, rms=main_rms,
                  rounds=main_rounds),
        summary=summary, draws=rows), indent=1))
    print(f'main draw (seed {MAIN_SEED}, {N_PIXELS} pixels): slope {main_slope:.4f} +- {main_u:.4f} (fit), rms '
          f'{main_rms:.1f} counts, {main_rounds} rounds -> pixels_frames.csv')
    for s_ in summary:
        print(f'{s_["n"]:3d} pixels, {s_["draws"]:4d} draws ({s_["failed"]} failed): slope mean {s_["mean"]:.3f}, '
              f'SD between draws {s_["sd"]:.3f} ({100 * s_["sd"] / s_["mean"]:.1f} %), median fit error {s_["median_u_fit"]:.3f}')
    print(f'-> pixels.json; {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
