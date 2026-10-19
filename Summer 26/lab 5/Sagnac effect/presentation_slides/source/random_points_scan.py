"""Random-points method repeated for many random draws (for the random-pixels slide).

Uses the functions of sagnac_update/methods.py unchanged (fringe region, pixel draw, dark level,
alternating least-squares fit) and the same averaging and weighted line as sagnac.py (methods.slope).
The only differences from methods.py are the number of draws and, as an extra check, draws with
more than ten pixels.

    python random_points_scan.py [draws]        (default 1000 draws of 10 pixels)

Writes random_points_scan.json next to this script. Runs in a few minutes.
"""
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage

HERE = Path(__file__).resolve().parent
DATA = Path(os.environ.get('SAGNAC_DATA', HERE.parents[1]))          # the 'Sagnac effect' folder
sys.path[:0] = [str(DATA / 'sagnac_update'), str(DATA)]
import methods as M                                                  # noqa: E402  (sagnac_update/methods.py)
from fringes import SMOOTH_PX, load_red                              # noqa: E402

N_DRAWS = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
EXTRA_N = [(20, 60), (40, 60), (80, 60)]                          # (pixels per draw, draws)
WORKERS = max(1, os.cpu_count() or 1)
G = {}                                                               # shared with the worker processes (fork)


def one_draw(i):
    """Fit draw i with the local phase from the rest geometry ('fixed', the method as used) and, as a
    diagnostic, from each photograph's own fitted centre and period ('frame')."""
    (sd, n), d = G['plan'][i], G['draws'][i]
    V = G['V_all'][:, G['offsets'][i]:G['offsets'][i + 1]].astype(float)
    psi_rest = 2 * np.pi * ((d[1] - G['cx0']) ** 2 + (d[0] - G['cy0']) ** 2) / G['P0']
    psi_frame = 2 * np.pi * ((d[1][None] - G['CX']) ** 2 + (d[0][None] - G['CY']) ** 2) / G['PP']
    row = dict(seed=int(sd), n=int(n))
    for key, psi in (('fixed', psi_rest), ('frame', psi_frame)):
        ph, _, _, _, rms, it = M.point_phases(G['fr'], V, psi, G['dark'])
        a, u, sx = M.slope(G['fr'], ph, G['model'])
        row[key] = dict(slope=float(a), u_fit=float(u), sigma_x=float(sx), rms=float(rms), iterations=int(it))
    if i % 50 == 0:
        print(f'{i}/{len(G["plan"])}  seed {sd} n {n}: {row["fixed"]["slope"]:.3f}  ({time.time() - G["t0"]:.0f} s)',
              flush=True)
    return row


def linear_model():
    """The linear calibration used by sagnac_update (calibration.json, fit with eps, degree 1)."""
    c = json.loads((DATA / 'calibration.json').read_text())
    lin = c['fits_with_eps']['1']
    return dict(degree=1, coef=lin['coef'], cov=lin['cov'], eps=c['conversion']['eps'])


def main():
    t0 = time.time()
    fr = pd.read_csv(DATA / 'fringes_frames.csv')
    fr['file'] = [str(DATA / f) for f in fr.file]
    rest = (fr.S == 0).to_numpy()
    model = linear_model()
    cx0, cy0, P0 = fr.cx[rest].median(), fr.cy[rest].median(), fr.P_px2[rest].median()
    rest_img = np.mean([ndimage.gaussian_filter(load_red(f), SMOOTH_PX) for f in fr.file[rest]], axis=0)
    dark = M.dark_level(rest_img)
    ys, xs = M.fringe_region(rest_img.shape, cx0, cy0)

    # all draws: the original one, seeds 1..N_DRAWS (seeds 1-30 are the scan in methods.json),
    # then draws with more pixels (seeds offset so they do not repeat the ten-pixel draws)
    plan = [(M.POINTS_SEED, 10)] + [(sd, 10) for sd in range(1, N_DRAWS + 1)]
    for n, k in EXTRA_N:
        plan += [(100000 * n + sd, n) for sd in range(1, k + 1)]
    draws = [M.draw_points(ys, xs, sd, n) for sd, n in plan]
    Y = np.concatenate([d[0] for d in draws]); X = np.concatenate([d[1] for d in draws])
    V_all = np.empty((len(fr), len(Y)), np.float32)
    for k, f in enumerate(fr.file):                                   # one pass over the photographs
        V_all[k] = ndimage.gaussian_filter(load_red(f), SMOOTH_PX)[Y, X]
    print(f'read {len(fr)} photographs, {len(Y)} pixels, {time.time() - t0:.0f} s', flush=True)

    CX, CY, PP = (fr[c].to_numpy()[:, None] for c in ('cx', 'cy', 'P_px2'))
    offsets = np.r_[0, np.cumsum([n for _, n in plan])]
    G.update(fr=fr, model=model, dark=dark, cx0=cx0, cy0=cy0, P0=P0, CX=CX, CY=CY, PP=PP, V_all=V_all,
             plan=plan, draws=draws, offsets=offsets, t0=t0)
    rows = []
    with mp.get_context('fork').Pool(WORKERS) as pool:
        for row in pool.imap(one_draw, range(len(plan)), chunksize=4):
            rows.append(row)
            if len(rows) % 100 == 0:                                  # partial results, in case of interruption
                (HERE / 'random_points_scan.partial.json').write_text(json.dumps(rows))

    out = dict(note='Random-points method of sagnac_update/methods.py repeated for many draws. '
                    "'fixed': local phase from the rest geometry (the method as used). "
                    "'frame': diagnostic, local phase from each photograph's own fitted centre and period.",
               dark=dark, centre=[cx0, cy0], P0=P0, max_iter=M.POINTS_MAX_ITER, draws=rows)
    (HERE / 'random_points_scan.json').write_text(json.dumps(out, indent=0))
    print(f'done, {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
