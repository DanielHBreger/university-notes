"""
Build Lorenz (first-return) maps from a folder of scope CSVs.

For each record:

  1. Smooth V1 lightly with a Savitzky-Golay filter, window ~ 1/20 of a period.
  2. Find the local maxima M_n with a prominence threshold of 2 % of the range.
  3. Plot M_{n+1} against M_n, with the diagonal drawn.

A period-1 limit cycle collapses to a single point on the diagonal, period-2
gives two points, and chaos fills out the characteristic one-humped curve.

The period that sets the smoothing window is measured per record from the
first peak of the signal's autocorrelation, since the sweep folders change
time base partway through (400 ns/sample early, 1 us/sample later) and a
fixed window would be wrong for half the records. Override it with
--period-samples if you want one fixed window throughout.

Records with no detectable period - the circuit not yet oscillating, where the
trace is pure noise - are skipped and listed in the summary; a prominence
threshold means nothing when there is no signal to take a fraction of.

If a "<folder>_rpot.csv" from batch_rpot.py sits next to the folder, the Rpot
values are picked up automatically and used to label the panels.

This module also owns `maxima` and `collect`, which bifurcation.py and
lyapunov.py both build on, so all three agree by construction about what a
record's maxima are and which records count.

Usage:
    python lorenz_map.py                          # pick the folder
    python lorenz_map.py "sweep forward"
    python lorenz_map.py "sweep forward" --each    # + one PNG per record
    python lorenz_map.py "sweep forward" --pooled  # one map, coloured by Rpot
    python lorenz_map.py "sweep forward" --ch CH2 --prominence 0.05
"""
import argparse
import csv
import os
from collections import namedtuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.fft import next_fast_len
from scipy.signal import find_peaks, savgol_filter

from sweeplib import (list_csvs, load_rpot, resolve_folders, sibling,
                      sweep_colors)

# One admitted record: its maxima and everything the callers filter on.
Record = namedtuple('Record', 'name rpot residual M info')

SUMMARY_FIELDS = ['filename', 'rpot_ohm', 'n_maxima', 'f0_hz',
                  'period_samples', 'sg_window', 'm_min', 'm_max',
                  'signal_range', 'clip_pct', 'status']


def period_samples(x, limit=200000, max_lag=50000):
    """
    Oscillation period in samples, from the first autocorrelation peak.

    Done through the FFT: np.correlate(y, y, 'full') is a direct O(N^2)
    convolution, which costs ~7 s per record here and turns a folder into a
    quarter-hour crawl. This is the same autocorrelation to the last sample.
    """
    y = x[:limit] - x[:limit].mean()
    s = y.std()
    if s == 0:
        return 0
    y = y / s
    n = next_fast_len(2 * len(y))
    F = np.fft.rfft(y, n)
    ac = np.fft.irfft(F * np.conj(F), n)[:max_lag]
    if ac[0] <= 0:
        return 0
    ac /= ac[0]
    p, _ = find_peaks(ac, height=0.1)
    return int(p[0]) if len(p) else 0


def _blank(dt=0.0, per=0, status='ok'):
    """The info record, with every field present whichever way maxima exits."""
    return {'dt': dt, 'period_samples': per, 'window': 0,
            'f0_hz': 1.0 / (per * dt) if per and dt else 0.0,
            'range': 0.0, 'clip_pct': 0.0, 'status': status,
            't_peaks': np.empty(0)}


def maxima(path, ch='CH1', prominence=0.02, period=None):
    """
    Return (M, info): the successive local maxima of `ch`.

    Smooth at period/20, then take the maxima above `prominence` * range.
    """
    # Only the time column and the one channel are needed; parsing the other
    # two costs about a fifth of the per-record time on a 28 MB record.
    head = pd.read_csv(path, nrows=0)
    names = [c.split('(')[0].strip() for c in head.columns]
    if ch not in names:
        raise ValueError(f'no column {ch!r} in {names}')
    k = names.index(ch)
    if k == 0:
        raise ValueError(f'{ch!r} is the time column')
    d = pd.read_csv(path, usecols=[0, k], dtype=np.float32)
    t = d.iloc[:, 0].values.astype(np.float64)
    x = d.iloc[:, 1].values.astype(np.float64)
    ok = np.isfinite(t) & np.isfinite(x)
    t, x = t[ok], x[ok]
    if len(x) < 100:
        raise ValueError('too few samples')

    dt = float(np.median(np.diff(t)))
    per = int(period) if period else period_samples(x)
    info = _blank(dt, per)
    # Fraction of samples pinned at either extreme: on this bench the low-R
    # records run off the scope's input range and clip hard, and a clipped
    # peak is a flat plateau at the rail rather than a real maximum. It is a
    # property of the raw trace, so it is worth reporting even for the
    # records that never get as far as having maxima.
    info['clip_pct'] = 100.0 * max(
        np.count_nonzero(x >= x.max() - 1e-3),
        np.count_nonzero(x <= x.min() + 1e-3)) / len(x)

    # A lag this short is the noise floor, not an orbit.
    if per < 20:
        info['period_samples'] = 0
        info['f0_hz'] = 0.0
        info['status'] = 'no periodicity found (not oscillating?)'
        return np.empty(0), info

    # Savitzky-Golay, window about a twentieth of a period, forced odd and >= 5.
    w = max(5, (per // 20) | 1)
    if w >= len(x):
        w = (len(x) - 1) | 1
    xs = savgol_filter(x, w, 3)
    info['window'] = w

    rng = float(xs.max() - xs.min())
    info['range'] = rng
    if rng <= 0:
        info['status'] = 'flat record'
        return np.empty(0), info

    p, _ = find_peaks(xs, prominence=prominence * rng)
    M = xs[p]
    # When each maximum happened, so return times can be measured from these
    # same peaks rather than re-detected somewhere else.
    info['t_peaks'] = t[p]
    if len(M) < 2:
        info['status'] = 'fewer than two maxima'
    return M, info


def collect(folder, ch='CH1', prominence=0.02, period=None,
            max_residual=10.0, max_clip=2.0, max_files=None):
    """
    The records of a sweep that count, with their maxima and their Rpot.

    Returns (records, dropped), where `dropped` pairs each rejected record
    with the reason. This is the single definition of "a record that counts",
    so the bifurcation diagram and the lambda plot are always built from the
    same set rather than from two copies of the same filter chain.
    """
    rpot = load_rpot(folder, required=True)
    files = list_csvs(folder)
    if max_files:
        files = files[:max_files]

    records, dropped = [], []
    for i, name in enumerate(files, 1):
        if name not in rpot:
            dropped.append((name, 'no Rpot entry'))
            continue
        r, resid = rpot[name]
        if resid > max_residual:
            dropped.append((name, f'Rpot fit residual {resid:.1f} %'))
            continue
        print(f'[{i}/{len(files)}] {name}', end=' ', flush=True)
        try:
            M, info = maxima(os.path.join(folder, name), ch, prominence, period)
        except Exception as e:
            print(f'-> ERROR: {e}')
            dropped.append((name, f'error: {e}'))
            continue
        if len(M) < 2:
            print(f'-> skipped ({info["status"]})')
            dropped.append((name, info['status']))
            continue
        if info['clip_pct'] > max_clip:
            print(f'-> clipped ({info["clip_pct"]:.0f} % of samples at the rail)')
            dropped.append((name,
                            f'clipped: {info["clip_pct"]:.0f} % of samples at the rail'))
            continue
        print(f'-> R={r:7.1f} ohm  {len(M):5d} maxima')
        records.append(Record(name, r, resid, M, info))
    return records, dropped


def frame(ax, lim):
    """The square return-map frame: identity diagonal and equal axes."""
    ax.plot(lim, lim, ls='--', lw=0.8, color='0.5', zorder=0)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect('equal', adjustable='box')


def draw(ax, M, title, lim, size=3, alpha=0.5):
    """One Lorenz map panel: M_{n+1} vs M_n, with the diagonal."""
    ax.scatter(M[:-1], M[1:], s=size, alpha=alpha, lw=0, color='C0')
    frame(ax, lim)
    ax.set_title(title, fontsize=8)


def write_csvs(out_dir, results, rpot):
    """Peak pairs and the per-record summary."""
    pairs_csv = os.path.join(out_dir, 'lorenz_pairs.csv')
    with open(pairs_csv, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['filename', 'rpot_ohm', 'n', 'm_n', 'm_next'])
        for name, M, _ in results:
            r = rpot.get(name, '')
            w.writerows((name, r, n, round(float(M[n]), 6),
                         round(float(M[n + 1]), 6)) for n in range(len(M) - 1))

    summary_csv = os.path.join(out_dir, 'lorenz_summary.csv')
    with open(summary_csv, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(SUMMARY_FIELDS)
        for name, M, info in results:
            w.writerow([name, rpot.get(name, ''), len(M),
                        round(info['f0_hz'], 1), info['period_samples'],
                        info['window'],
                        round(float(M.min()), 4) if len(M) else '',
                        round(float(M.max()), 4) if len(M) else '',
                        round(info['range'], 4), round(info['clip_pct'], 2),
                        info['status']])
    return pairs_csv, summary_csv


def plot_grid(out_dir, usable, rpot, lim, ch, folder, rows, cols):
    """Pages of small panels, one per record, on shared axes."""
    per_page = cols * rows
    pages = (len(usable) + per_page - 1) // per_page
    for pg in range(pages):
        chunk = usable[pg * per_page:(pg + 1) * per_page]
        fig, axes = plt.subplots(rows, cols, figsize=(2.3 * cols, 2.5 * rows))
        for ax, (name, M, _) in zip(np.ravel(axes), chunk):
            r = rpot.get(name, '')
            lab = f'{name}\n{r:.0f} $\\Omega$' if r != '' else name
            draw(ax, M, lab, lim)
        for ax in np.ravel(axes)[len(chunk):]:
            ax.axis('off')
        fig.suptitle(f'Lorenz map, {ch} maxima - {os.path.basename(folder)}'
                     + (f' (page {pg + 1}/{pages})' if pages > 1 else ''))
        fig.supxlabel('$M_n$  (V)')
        fig.supylabel('$M_{n+1}$  (V)')
        fig.tight_layout()
        out = os.path.join(out_dir, f'lorenz_grid_{pg + 1:02d}.png')
        fig.savefig(out, dpi=140)
        plt.close(fig)
        print(f'wrote {out}')


def plot_each(out_dir, usable, rpot, lim):
    """One full-size map per record."""
    each_dir = os.path.join(out_dir, 'each')
    os.makedirs(each_dir, exist_ok=True)
    for name, M, _ in usable:
        r = rpot.get(name, '')
        fig, ax = plt.subplots(figsize=(4.5, 4.5))
        draw(ax, M, name + (f'   Rpot = {r:.1f} $\\Omega$' if r != '' else ''),
             lim, size=6, alpha=0.6)
        ax.set_xlabel('$M_n$  (V)')
        ax.set_ylabel('$M_{n+1}$  (V)')
        fig.tight_layout()
        fig.savefig(os.path.join(each_dir, os.path.splitext(name)[0] + '.png'),
                    dpi=150)
        plt.close(fig)
    print(f'wrote {len(usable)} PNGs to {each_dir}')


def plot_pooled(out_dir, usable, rpot, lim, folder):
    """Every record on one map, coloured by Rpot."""
    fig, ax = plt.subplots(figsize=(6, 5.4))
    vals = [rpot[n] for n, _, _ in usable if n in rpot]
    if vals:
        norm = plt.Normalize(min(vals), max(vals))
        cmap = plt.get_cmap('viridis')
        for name, M, _ in usable:
            if name in rpot:
                ax.scatter(M[:-1], M[1:], s=2, alpha=0.4, lw=0,
                           color=cmap(norm(rpot[name])))
        fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                     label='$R_{pot}$  ($\\Omega$)')
    else:
        for _, M, _ in usable:
            ax.scatter(M[:-1], M[1:], s=2, alpha=0.4, lw=0)
    frame(ax, lim)
    ax.set_xlabel('$M_n$  (V)')
    ax.set_ylabel('$M_{n+1}$  (V)')
    ax.set_title(f'Lorenz map, all records - {os.path.basename(folder)}')
    fig.tight_layout()
    out = os.path.join(out_dir, 'lorenz_pooled.png')
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f'wrote {out}')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('folder', nargs='?', help='folder of scope CSVs')
    p.add_argument('--out-dir', help='output dir (default: <folder>_lorenz)')
    p.add_argument('--ch', default='CH1', help='channel to take maxima of')
    p.add_argument('--prominence', type=float, default=0.02,
                   help='peak prominence, fraction of the range')
    p.add_argument('--period-samples', type=int,
                   help='fixed period in samples (default: per-record autocorr)')
    p.add_argument('--each', action='store_true',
                   help='also write one full-size PNG per record')
    p.add_argument('--pooled', action='store_true',
                   help='also write one combined map coloured by Rpot')
    p.add_argument('--cols', type=int, default=6, help='panels per grid row')
    p.add_argument('--rows', type=int, default=4, help='panel rows per page')
    p.add_argument('--max-files', type=int, help='stop after this many records')
    a = p.parse_args()

    folder = resolve_folders([a.folder] if a.folder else [], multi=False)[0]
    files = list_csvs(folder)
    if a.max_files:
        files = files[:a.max_files]
    if not files:
        raise SystemExit(f'no CSV files in {folder}')

    out_dir = a.out_dir or sibling(folder, '_lorenz')
    os.makedirs(out_dir, exist_ok=True)
    rpot = {k: v[0] for k, v in load_rpot(folder).items()}
    if rpot:
        print(f'labelling panels with Rpot from '
              f'{os.path.basename(folder)}_rpot.csv')

    results = []
    for i, name in enumerate(files, 1):
        print(f'[{i}/{len(files)}] {name}', end=' ', flush=True)
        try:
            M, info = maxima(os.path.join(folder, name), a.ch, a.prominence,
                             a.period_samples)
        except Exception as e:
            print(f'-> ERROR: {e}')
            results.append((name, np.empty(0), _blank(status=f'error: {e}')))
            continue
        print(f'-> {len(M):5d} maxima  win={info["window"]:3d}  '
              f'f0={info["f0_hz"]:7.0f} Hz  {info["status"]}')
        results.append((name, M, info))

    pairs_csv, summary_csv = write_csvs(out_dir, results, rpot)

    usable = [r for r in results if len(r[1]) > 2]
    if not usable:
        raise SystemExit('no record produced a usable set of maxima')

    # Common axis limits so panels are comparable across the sweep.
    allM = np.concatenate([M for _, M, _ in usable])
    lo, hi = float(allM.min()), float(allM.max())
    pad = 0.05 * (hi - lo)
    lim = (lo - pad, hi + pad)

    plot_grid(out_dir, usable, rpot, lim, a.ch, folder, a.rows, a.cols)
    if a.each:
        plot_each(out_dir, usable, rpot, lim)
    if a.pooled:
        plot_pooled(out_dir, usable, rpot, lim, folder)

    skipped = len(results) - len(usable)
    print(f'\n{len(usable)}/{len(results)} records mapped'
          + (f'   ({skipped} skipped)' if skipped else '')
          + f'\npairs   -> {pairs_csv}\nsummary -> {summary_csv}')


if __name__ == '__main__':
    main()
