"""
Largest Lyapunov exponent straight from the time series (Rosenstein 1993).

This is the "proper" calculation the experiment plan says to compare the
return-map estimate against. It needs no map, no maxima and no assumption
that the flow has collapsed to a one-dimensional function: it watches how
fast nearby points of the reconstructed attractor separate.

  1. Reconstruct the state. The measured channels are V1 and V2; the third
     state variable (the inductor current) is not measured, so it is stood
     in for by a delayed copy of V1, a quarter period back. The three
     coordinates are scaled to unit spread so V2, which swings a tenth of
     V1, counts equally.
  2. For every reference point find its nearest neighbour that is at least
     one period away in time. The temporal exclusion (the Theiler window)
     matters: without it the nearest neighbour is simply the next sample on
     the same orbit, which never diverges.
  3. Follow each pair forward k samples and average ln of their distance.
     The mean log-distance rises linearly while the divergence is
     exponential, then bends over when the pairs are as far apart as the
     attractor is wide. The slope of the straight part, per second, is
     lambda.

Quantisation. The scope resolves both channels in 43 mV steps. On the raw
codes the nearest neighbour of almost every point sits on the same code,
the measured distance is quantisation noise for the first periods, and the
slope reads low: on a simulated single scroll with this step the estimate
comes out at half the true exponent. The channels are therefore smoothed
first with the same Savitzky-Golay filter the maxima are found on (window
a twentieth of a period, --smooth), which averages the step down over
~17 samples without touching the dynamics. With it, simulated double
scrolls quantised like this bench's records (benchmark_lyapunov.py) come
out at 0.95 to 1.15 times the true exponent at the 43 mV step and about
1.2 times it at the coarse 201 mV step; a single scroll, whose exponent is
small, reads about 0.85. Quote it to about 20 %.

The fit window is a fixed range in periods (--fit, default 0.5 to 2.5).
Starting after half a period skips the initial stretch where the nearest
pairs are atypically close and relax outward faster than lambda; stopping
before the curve bends avoids saturation. The plot shows both the curve
and the window, so the choice can be checked by eye for every record. The
record is split into four blocks in time whose separate slopes give the
quoted uncertainty.

Usage:
    python rosenstein.py "sweep forward/trace90.csv"
    python rosenstein.py "chaos - 716.9 ohm.csv" --fit 0.5 2
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.spatial import cKDTree


def load_channels(path):
    """(t, V1, V2) from a scope CSV, finite rows only."""
    d = pd.read_csv(path, usecols=[0, 1, 2], dtype=np.float32)
    t = d.iloc[:, 0].values.astype(np.float64)
    v1 = d.iloc[:, 1].values.astype(np.float64)
    v2 = d.iloc[:, 2].values.astype(np.float64)
    ok = np.isfinite(t) & np.isfinite(v1) & np.isfinite(v2)
    return t[ok], v1[ok], v2[ok]


def embed(v1, v2, tau, smooth_window=0):
    """
    3-D state (V1, V2, V1 delayed by tau samples), each scaled to unit std.

    With `smooth_window` (odd, >= 5) both channels are Savitzky-Golay
    smoothed first, which is what takes the quantisation out.
    """
    if smooth_window and smooth_window >= 5:
        w = int(smooth_window) | 1
        v1 = savgol_filter(v1, w, 3)
        v2 = savgol_filter(v2, w, 3)
    n = len(v1) - tau
    X = np.column_stack([v1[tau:], v2[tau:], v1[:n]])
    X -= X.mean(0)
    s = X.std(0)
    s[s == 0] = 1.0
    return X / s


def nearest_recurrent(X, last, theiler, stride, tree_stride, k_query):
    """
    Index pairs (i, j): each reference i and its nearest neighbour j with
    |i - j| > theiler, both below `last` so they can be followed forward.

    The tree holds every `tree_stride`-th sample; nothing needs the neighbour
    to be sample-exact, and it keeps the query manageable. Samples inside the
    Theiler window are nearly always the closest, so many neighbours are
    asked for at once and the first admissible one is taken; a reference
    whose first batch is all excluded is retried with a much longer list.
    """
    tree_idx = np.arange(0, last, tree_stride)
    tree = cKDTree(X[tree_idx])
    refs = np.arange(0, last, stride)
    partner = np.full(len(refs), -1)

    def query(rows, k):
        k = min(k, len(tree_idx))
        _, jj = tree.query(X[refs[rows]], k=k)
        jj = np.atleast_2d(jj)
        j = tree_idx[jj]
        ok = np.abs(j - refs[rows][:, None]) > theiler
        found = ok.any(1)
        first = ok.argmax(1)
        partner[rows[found]] = j[np.arange(len(rows)), first][found]
        return rows[~found]

    for c in range(0, len(refs), 20000):
        rows = np.arange(c, min(c + 20000, len(refs)))
        missing = query(rows, k_query)
        if len(missing):
            query(missing, 8 * k_query)
    keep = partner >= 0
    return refs[keep], partner[keep]


def divergence(X, theiler, k_max, stride=4, tree_stride=4, k_query=256,
               blocks=4, k_step=1):
    """
    <ln d(k)> for k = 0, k_step, ... up to k_max, overall and per time block.

    Returns (ks, curve, block_curves, n_pairs). A zero distance (identical
    smoothed states; rare) is left out of that step's average.
    """
    n = len(X)
    last = n - k_max - 1
    if last <= 2 * theiler:
        raise ValueError('record too short for this Theiler window')
    i, j = nearest_recurrent(X, last, theiler, stride, tree_stride, k_query)
    if len(i) < 100:
        raise ValueError('too few usable pairs')
    block = np.minimum(i * blocks // last, blocks - 1)
    ks = np.arange(0, k_max + 1, max(int(k_step), 1))
    sums = np.zeros((blocks, len(ks)))
    counts = np.zeros((blocks, len(ks)))
    for c, k in enumerate(ks):
        d = np.linalg.norm(X[i + k] - X[j + k], axis=1)
        good = d > 0
        sums[:, c] = np.bincount(block[good], weights=np.log(d[good]),
                                 minlength=blocks)
        counts[:, c] = np.bincount(block[good], minlength=blocks)
    with np.errstate(invalid='ignore', divide='ignore'):
        block_curves = sums / counts
        curve = sums.sum(0) / counts.sum(0)
    return ks, curve, block_curves, len(i)


def rosenstein(v1, v2, dt, period, tau=None, theiler=None, fit=(0.5, 2.5),
               follow=3.5, stride=None, smooth=20):
    """
    Largest Lyapunov exponent of the record, in 1/s.

    `period` is the winding period in samples (lorenz_map.period_samples).
    The channels are smoothed over period/`smooth` samples (0 to skip), the
    pairs followed for `follow` periods and the slope fitted over the `fit`
    window, in periods. `stride` (reference points every n-th sample)
    defaults to whatever gives about 50 000 references. Returns a dict with
    'lam', 'lam_err', the divergence curve and everything needed to draw it.
    """
    period = int(period)
    tau = int(tau) if tau else max(period // 4, 1)
    theiler = int(theiler) if theiler else period
    stride = int(stride) if stride else max(1, round(len(v1) / 50000))
    window = max(5, (period // smooth) | 1) if smooth else 0
    k_max = int(follow * period)
    X = embed(v1, v2, tau, window)
    ks, curve, blocks, n_pairs = divergence(X, theiler, k_max, stride,
                                            k_step=max(1, period // 200))
    t = ks * dt
    sel = (ks >= fit[0] * period) & (ks <= fit[1] * period)
    if sel.sum() < 10:
        raise ValueError('fit window too short')
    slope, icpt = np.polyfit(t[sel], curve[sel], 1)
    block_slopes = []
    for b in blocks:
        ok = sel & np.isfinite(b)
        if ok.sum() > 10:
            block_slopes.append(float(np.polyfit(t[ok], b[ok], 1)[0]))
    err = (float(np.std(block_slopes, ddof=1)) / np.sqrt(len(block_slopes))
           if len(block_slopes) > 1 else np.nan)
    return {'lam': float(slope), 'lam_err': err, 'intercept': float(icpt),
            't': t, 'curve': curve, 'fit_sel': sel, 'n_pairs': n_pairs,
            'tau': tau, 'theiler': theiler, 'stride': stride,
            'smooth_window': window, 'block_slopes': block_slopes}


def draw_curve(ax, res, dt_period, title=None):
    """The divergence curve with the fitted window, time in periods."""
    tp = res['t'] / dt_period
    ax.plot(tp, res['curve'], color='0.3', lw=1)
    sel = res['fit_sel']
    ax.plot(tp[sel], res['intercept'] + res['lam'] * res['t'][sel],
            color='C3', lw=1.5,
            label=f"$\\lambda$ = {res['lam']:.0f} $\\pm$ {res['lam_err']:.0f} s$^{{-1}}$")
    ax.set_xlabel('time  (periods)')
    ax.set_ylabel('$\\langle \\ln d \\rangle$')
    ax.legend(fontsize=8, loc='lower right')
    ax.set_title(title or 'direct (Rosenstein): mean log divergence of neighbours',
                 fontsize=9)


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from lorenz_map import period_samples

    p = argparse.ArgumentParser()
    p.add_argument('files', nargs='+', help='scope CSV records')
    p.add_argument('--fit', nargs=2, type=float, default=(0.5, 2.5),
                   metavar=('LO', 'HI'), help='fit window, in periods')
    p.add_argument('--follow', type=float, default=3.5,
                   help='how long to follow each pair, in periods')
    p.add_argument('--smooth', type=int, default=20,
                   help='smoothing window is period / this (0: none)')
    p.add_argument('--stride', type=int,
                   help='reference point every n-th sample (default ~50k refs)')
    p.add_argument('--tau', type=int, help='embedding delay, samples')
    p.add_argument('--theiler', type=int, help='temporal exclusion, samples')
    a = p.parse_args()

    for path in a.files:
        t, v1, v2 = load_channels(path)
        dt = float(np.median(np.diff(t)))
        per = period_samples(v1)
        if per < 20:
            print(f'{path}: no periodicity found')
            continue
        res = rosenstein(v1, v2, dt, per, a.tau, a.theiler, tuple(a.fit),
                         a.follow, a.stride, a.smooth)
        print(f'{os.path.basename(path)}: period {per} samples '
              f'({per * dt * 1e6:.0f} us), {res["n_pairs"]} pairs, '
              f'lambda = {res["lam"]:.0f} +/- {res["lam_err"]:.0f} /s '
              f'(blocks: {", ".join(f"{b:.0f}" for b in res["block_slopes"])})')
        fig, ax = plt.subplots(figsize=(6, 4))
        draw_curve(ax, res, per * dt, os.path.basename(path))
        fig.tight_layout()
        out = os.path.splitext(path)[0] + '_rosenstein.png'
        fig.savefig(out, dpi=140)
        plt.close(fig)
        print(f'  -> {out}')


if __name__ == '__main__':
    main()
