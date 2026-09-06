"""
Lyapunov exponent estimated from the first-return map (M5), with the direct
time-series calculation alongside it.

    lambda ~= < ln |f'(M_n)| > / <T>

f'(M_n) is the local slope of the return map at each visited maximum, from a
linear fit over a window sliding along M_n, and <T> is the mean return time
between successive maxima. Both come from the same peaks, found by
lorenz_map.maxima, so this plot, the Lorenz maps and the bifurcation diagram
agree by construction about what a record's maxima are.

The two traps the plan names, and how they are handled:

  <T> is the MEAN, not the median. A lobe switch takes longer than an
  ordinary winding, so the return-time distribution has a tail and the
  median sits below the mean. Both are reported per record and the gap is
  summarised by regime, so the choice can be defended with the data.

  ln |f'| is averaged over the points ACTUALLY VISITED, including the
  neighbourhood of the turning point where the slope falls toward zero.
  That region pulls the average down and must not be left out - which is
  the trap of this estimator on noisy data, because a slope that small is
  not distinguishable from zero. A slope smaller than its own standard
  error is therefore floored at se/e: for a slope passing uniformly
  through zero within +-se the exact mean of ln|f'| is ln(se) - 1, so the
  floor reproduces that expectation without letting one fluke dominate.
  The share of such points is reported (unresolved_pct).

What it takes to have a slope at all. The formula assumes M_(n+1) is a
function of M_n. A periodic orbit is a few tight clusters of noise, and
the fit is judged WITHIN each branch: against the pooled variance a local
fit that explains nothing inside either cluster still scores R2 ~ 0.99,
because the pooled variance is the distance between the clusters. So R2 is
taken against each branch's own variance (--min-r2), a branch narrower than
--min-spread is a cluster and is not fitted (a period-2 orbit with slow
amplitude drift smears each cluster into a streak that can pass R2 on its
own), and a record keeps its lambda only if at least --min-curve-frac of
its points lie on fitted branches. Periodic windows carry no lambda by
design; it is not that lambda is negative there, it is that this method
cannot measure it.

The window is a fixed WIDTH in M_n (--width, a fraction of the branch
span, floored at --min-width and at three quantisation steps). A fixed
number of nearest points squeezes the x-range of the window to a sliver
while M_(n+1) keeps its full noise, and the slope of pure noise then comes
out above 1.

Quantisation. The scope resolves V1 in 43 mV steps, so a maximum is known
to about that. A map spanning a few volts (the double scroll) is fine; one
spanning a few tenths of a volt (just past the period-doubling cascade, or
a single scroll) is a handful of steps wide, and on a simulated single
scroll with this step the local slopes are fiction: the estimate came out
at 1.7 times the true exponent with clean maxima giving 0.6. On synthetic
one-hump maps quantised the same way the estimate is within 10 % from
about 35 steps up and 20-30 % low at 23 steps, so a map narrower than
--min-steps (30) quantisation steps is refused.

Uncertainty on lambda: the standard error of the mean of ln|f'| (successive
maxima of a chaotic orbit decorrelate within a few windings), half the
change when the window is doubled, and the standard error of <T>, combined
in quadrature. It is a precision, not an accuracy.

Accuracy, and the direct calculation. The return map of maxima of V1 is
only approximately a one-dimensional function of M_n: the lobe switch and
the fold of the band give it a near-vertical stretch whose slope is a
property of the projection, not of the flow, and such a map reads high. On
simulated Chua circuits with this bench's parameters and quantisation
(benchmark_lyapunov.py) the map estimate lands between 1.1 and 1.7 times
the true exponent on the double scroll, while the Rosenstein calculation
on the (smoothed) time series lands within about 20 % on either side. So
--rosenstein computes the direct estimate for every record that has a map
lambda (rosenstein.py) and plots it alongside. Quote the map value as what
it is - an estimate from the return map, an upper figure - and the direct
value as the measurement.

Records: any number of sweep folders (chosen by lorenz_map.collect, the
same gate the bifurcation diagram uses) and/or single CSV records named
with their resistance, like "doublescroll 2 - 656.7ohm.csv".

Usage:
    python lyapunov.py                              # pick folders, repeatedly
    python lyapunov.py "sweep forward" "sweep back"
    python lyapunov.py "sweep forward" --rosenstein --each
    python lyapunov.py "chaos - 716.9 ohm.csv" "doublescroll 2 - 656.7ohm.csv" --rosenstein --each
"""
import argparse
import csv
import os
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from bifurcation import add_record_args
from lorenz_map import Record, collect, frame, maxima
from rosenstein import draw_curve, load_channels, rosenstein
from sweeplib import resolve_folders, sibling, sweep_colors

CSV_FIELDS = [
    ('sweep', 'sweep', None), ('filename', 'filename', None),
    ('rpot_ohm', 'rpot', 2), ('n_maxima', 'n_maxima', None),
    ('n_branches', 'n_branches', None), ('n_slopes', 'n_used', None),
    ('unresolved_pct', 'unresolved_pct', 1), ('cluster_pct', 'cluster_pct', 1),
    ('mean_T_us', 'mean_T_us', 4), ('median_T_us', 'median_T_us', 4),
    ('mean_vs_median_pct', 'T_diff_pct', 2),
    ('mean_ln_abs_fprime', 'mean_ln_slope', 5),
    ('lambda_map_per_s', 'lam', 2), ('lambda_map_err_per_s', 'lam_err', 2),
    ('lambda_map_err_stat', 'lam_stat', 2),
    ('lambda_map_err_width', 'lam_sys', 2),
    ('map_r2_within_branch', 'r2', 4),
    ('lambda_rosenstein_per_s', 'lam_ros', 2),
    ('lambda_rosenstein_err_per_s', 'lam_ros_err', 2),
    ('status', 'status', None),
]


def split_branches(m_n, gap):
    """
    Indices of `m_n` sorted by value, cut into branches at large gaps.

    A double-scroll map lives on two separated clusters; a fit that crosses
    the void between them measures the void. Each branch comes back already
    in ascending m_n order, so callers need not re-sort.
    """
    order = np.argsort(m_n, kind='stable')
    v = m_n[order]
    span = float(v[-1] - v[0])
    if span <= 0:
        return [order]
    cuts = np.flatnonzero(np.diff(v) > gap * span) + 1
    return np.split(order, cuts) if len(cuts) else [order]


def local_fit(x, y, width, min_points):
    """
    Local linear fit at every point over the window |x' - x| <= width / 2.

    `x` must be ascending. Returns (slope, fitted value, slope standard
    error) per point; all three are NaN where the window holds fewer than
    `min_points` points or has no spread in x. Prefix sums make each window
    O(1), so the cost is O(n) whatever the width.
    """
    n = len(x)
    if n == 0:
        return np.empty(0), np.empty(0), np.empty(0)
    h = width / 2.0
    lo = np.searchsorted(x, x - h, 'left')
    hi = np.searchsorted(x, x + h, 'right')
    z = np.zeros(1)
    cx = np.concatenate([z, np.cumsum(x)])
    cy = np.concatenate([z, np.cumsum(y)])
    cxx = np.concatenate([z, np.cumsum(x * x)])
    cxy = np.concatenate([z, np.cumsum(x * y)])
    cyy = np.concatenate([z, np.cumsum(y * y)])
    m = (hi - lo).astype(float)
    Sx, Sy = cx[hi] - cx[lo], cy[hi] - cy[lo]
    Sxx, Sxy, Syy = cxx[hi] - cxx[lo], cxy[hi] - cxy[lo], cyy[hi] - cyy[lo]
    vx = Sxx - Sx * Sx / m                      # sum of squares about x-mean
    ok = (m >= max(min_points, 3)) & (vx > 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        slope = (Sxy - Sx * Sy / m) / vx
        pred = Sy / m + slope * (x - Sx / m)
        ss_res = np.maximum(Syy - Sy * Sy / m - slope * slope * vx, 0.0)
        se = np.sqrt(ss_res / (m - 2) / vx)
    slope[~ok] = np.nan
    pred[~ok] = np.nan
    se[~ok] = np.nan
    return slope, pred, se


def fit_map(m_n, m_next, width, min_width, gap, min_spread, min_points):
    """
    Local slopes of the return map at every visited point.

    Returns a dict with per-point 'slopes', 'fit' and 'se' (NaN off the
    fitted branches), the within-branch 'r2', the number of fitted
    branches and the number of points that lie on them.
    """
    n = len(m_n)
    slopes = np.full(n, np.nan)
    fit = np.full(n, np.nan)
    se = np.full(n, np.nan)
    ss_res = ss_tot = 0.0
    n_branches = n_on_curve = 0
    for idx in split_branches(m_n, gap):
        if len(idx) < min_points:
            continue
        x, y = m_n[idx], m_next[idx]
        if float(np.ptp(x)) < min_spread:
            continue                    # a cluster, not a curve
        w = max(width * float(np.ptp(x)), min_width)
        s, pred, e = local_fit(x, y, w, min_points)
        good = np.isfinite(pred)
        if not good.any():
            continue
        n_branches += 1
        n_on_curve += len(idx)
        slopes[idx], fit[idx], se[idx] = s, pred, e
        ss_res += float(np.sum((y[good] - pred[good]) ** 2))
        ss_tot += float(np.sum((y[good] - y[good].mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return {'slopes': slopes, 'fit': fit, 'se': se, 'r2': r2,
            'n_branches': n_branches, 'n_on_curve': n_on_curve}


def mean_ln_slope(slopes, se):
    """
    <ln|f'|> over the fitted points, its standard error, and the share of
    points whose slope was not resolved from noise (|f'| <= se), floored
    at se/e.
    """
    ok = np.isfinite(slopes) & np.isfinite(se)
    a = np.abs(slopes[ok])
    floor = se[ok] / np.e
    unresolved = a <= se[ok]
    a = np.maximum(a, floor)
    a = a[a > 0]
    if len(a) < 2:
        return np.nan, np.nan, np.nan, 0
    ln = np.log(a)
    return (float(ln.mean()), float(ln.std(ddof=1) / np.sqrt(len(ln))),
            100.0 * float(unresolved.mean()), len(ln))


def lyapunov(M, t_peaks, width=0.1, min_width=0.02, gap=0.05,
             min_spread=0.1, min_r2=0.8, min_points=10, min_curve_frac=0.5,
             quantum=0.0, min_steps=30):
    """
    Return-time statistics, <ln|f'|> and lambda for one record.

    M is the sequence of maxima and t_peaks the time of each; `quantum` is
    the scope's step on the channel. The per-point slopes, fitted map and
    standard errors come back under 'slopes', 'fit' and 'se', indexed like
    M[:-1], so the diagnostic plot shows exactly what was measured.
    """
    min_width = max(min_width, 3.0 * quantum)
    n_pairs = max(len(M) - 1, 0)
    out = {'n_maxima': len(M), 'status': 'ok', 'n_branches': 0,
           'mean_T': np.nan, 'median_T': np.nan, 'T_diff_pct': np.nan,
           'mean_ln_slope': np.nan, 'lam': np.nan, 'lam_err': np.nan,
           'lam_stat': np.nan, 'lam_sys': np.nan, 'n_used': 0,
           'unresolved_pct': np.nan, 'cluster_pct': np.nan, 'spread': 0.0,
           'r2': np.nan, 'lam_ros': np.nan, 'lam_ros_err': np.nan,
           'slopes': np.full(n_pairs, np.nan), 'fit': np.full(n_pairs, np.nan),
           'se': np.full(n_pairs, np.nan)}
    if len(M) < 2 * min_points + 2 or len(t_peaks) != len(M):
        out['status'] = 'too few maxima for a local fit'
        return out

    # Return times, from the same peaks the map is built from.
    T = np.diff(t_peaks)
    T = T[np.isfinite(T) & (T > 0)]
    if len(T) < 2:
        out['status'] = 'no usable return times'
        return out
    mean_T, med_T = float(T.mean()), float(np.median(T))
    out['mean_T'], out['median_T'] = mean_T, med_T
    out['T_diff_pct'] = 100.0 * (mean_T - med_T) / mean_T if mean_T else np.nan
    rel_T = float(T.std(ddof=1) / np.sqrt(len(T)) / mean_T)

    m_n, m_next = M[:-1], M[1:]
    spread = float(m_n.max() - m_n.min())
    out['spread'] = spread
    if spread < min_spread:
        # A period-1 orbit collapses the map onto one point: nothing to fit.
        out['status'] = f'map spread {spread:.3f} V below --min-spread'
        return out
    if quantum > 0 and spread < min_steps * quantum:
        out['status'] = (f'map spans only {spread / quantum:.0f} quantisation '
                         f'steps of {quantum * 1e3:.0f} mV (under --min-steps)')
        return out

    fm = fit_map(m_n, m_next, width, min_width, gap, min_spread, min_points)
    out.update(slopes=fm['slopes'], fit=fm['fit'], se=fm['se'], r2=fm['r2'],
               n_branches=fm['n_branches'],
               cluster_pct=100.0 * (1.0 - fm['n_on_curve'] / n_pairs))
    if not fm['n_branches']:
        out['status'] = ('map is clusters, not a curve (no branch spans '
                         f'{min_spread:.2f} V with {min_points} points)')
        return out
    if fm['n_on_curve'] < min_curve_frac * n_pairs:
        out['status'] = (f'only {100 * fm["n_on_curve"] / n_pairs:.0f} % of '
                         f'points lie on a branch wide enough to fit')
        return out
    if not np.isfinite(fm['r2']) or fm['r2'] < min_r2:
        out['status'] = (f'return map is not a curve within its branches '
                         f'(local fit R2 = {fm["r2"]:.2f})')
        return out

    mean_ln, sem, unresolved_pct, n_used = mean_ln_slope(fm['slopes'], fm['se'])
    if n_used < min_points:
        out['status'] = 'too few usable slopes'
        return out
    # Window sensitivity: the same average with the window doubled; half
    # the change is the systematic on <ln|f'|>. (Halving is not used: on
    # quantised maxima a half window holds one or two distinct M_n values.)
    f2 = fit_map(m_n, m_next, 2.0 * width, min_width, gap, min_spread,
                 min_points)
    alt = mean_ln_slope(f2['slopes'], f2['se'])[0]
    sys_ln = 0.5 * abs(alt - mean_ln) if np.isfinite(alt) else 0.0

    out['mean_ln_slope'] = mean_ln
    out['unresolved_pct'] = unresolved_pct
    out['n_used'] = n_used
    out['lam'] = mean_ln / mean_T
    out['lam_stat'] = sem / mean_T
    out['lam_sys'] = sys_ln / mean_T
    out['lam_err'] = float(np.hypot(np.hypot(out['lam_stat'], out['lam_sys']),
                                    abs(out['lam']) * rel_T))
    return out


def diagnostic(path, M, t_peaks, res, title, gap, ros=None, period_s=None):
    """
    Return map with its local fit and slopes, the return-time histogram,
    and the Rosenstein divergence curve when it was computed.
    """
    m_n, m_next = M[:-1], M[1:]
    slopes, fit, se = res['slopes'], res['fit'], res['se']
    n_panels = 3 if ros else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(5.5 * n_panels, 4.4))
    a1, a2 = axes[0], axes[1]

    fitted = np.isfinite(slopes)
    unres = fitted & (np.abs(slopes) <= se)
    if unres.any():
        a1.scatter(m_n[unres], m_next[unres], s=5, color='0.6', lw=0,
                   label=f'|f\'| below noise ({int(unres.sum())})')
    shown = fitted & ~unres
    sc = a1.scatter(m_n[shown], m_next[shown], c=np.log(np.abs(slopes[shown])),
                    s=4, cmap='coolwarm', lw=0)
    fig.colorbar(sc, ax=a1, label="$\\ln|f'(M_n)|$")
    for idx in split_branches(m_n, gap):
        ok = np.isfinite(fit[idx])
        if ok.sum() > 1:
            a1.plot(m_n[idx][ok], fit[idx][ok], color='k', lw=0.8, alpha=0.7)
    lo = float(min(m_n.min(), m_next.min()))
    hi = float(max(m_n.max(), m_next.max()))
    pad = 0.06 * (hi - lo)
    frame(a1, (lo - pad, hi + pad))
    a1.set_xlabel('$M_n$  (V)')
    a1.set_ylabel('$M_{n+1}$  (V)')
    a1.set_title(f'return map, local fit (black) and slope   '
                 f'R$^2$ = {res["r2"]:.2f}', fontsize=9)
    if unres.any():
        a1.legend(fontsize=7, loc='upper left')
    if not fitted.any():
        a1.scatter(m_n, m_next, s=4, color='C0', lw=0, alpha=0.6)
        a1.set_title(f'return map: {res["status"][:60]}', fontsize=8)

    T = np.diff(t_peaks) * 1e6
    if len(T) > 1:
        a2.hist(T, bins=80, color='0.7', edgecolor='none')
        a2.axvline(res['mean_T'] * 1e6, color='C3', lw=1.5,
                   label=f"mean {res['mean_T'] * 1e6:.1f} $\\mu$s")
        a2.axvline(res['median_T'] * 1e6, color='C0', lw=1.5, ls='--',
                   label=f"median {res['median_T'] * 1e6:.1f} $\\mu$s")
        a2.legend(fontsize=8)
    a2.set_xlabel('return time  ($\\mu$s)')
    a2.set_ylabel('count')
    a2.set_title('return times: the mean is the one to use', fontsize=9)

    if ros:
        draw_curve(axes[2], ros, period_s)

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def num(v, nd):
    """CSV cell: rounded when finite, blank when not."""
    if nd is None:
        return v
    return round(v, nd) if isinstance(v, float) and np.isfinite(v) else (
        v if not isinstance(v, float) else '')


def regime_report(label, rows):
    """
    The mean/median return-time gap, split by regime.

    It is a lobe-switch effect, so it belongs split this way: one branch means
    a single scroll with no lobe switches to skew the distribution, two or
    more means a double scroll where the long switching returns build the tail.
    """
    mine = [r for r in rows if r['sweep'] == label and r['status'] == 'ok']
    for tag, sel in (('single scroll (1 branch) ',
                      [r for r in mine if r['n_branches'] == 1]),
                     ('double scroll (2+ branch)',
                      [r for r in mine if r['n_branches'] >= 2])):
        d = [abs(r['T_diff_pct']) for r in sel if np.isfinite(r['T_diff_pct'])]
        if d:
            print(f'  {tag}: n={len(d):3d}  |mean-median| T '
                  f'avg {np.mean(d):5.2f} %, max {max(d):5.2f} %, '
                  f'{100 * np.mean(np.array(d) > 10):.0f} % of records over 10 %')


def rpot_from_name(name):
    """Resistance from a file name like 'doublescroll 2 - 656.7ohm.csv'."""
    m = re.search(r'(\d+(?:\.\d+)?)\s*ohm', name, re.IGNORECASE)
    return float(m.group(1)) if m else np.nan


def single_records(paths, ch, prominence, period):
    """Named CSV records as a pseudo-sweep: (label, records, dropped)."""
    records, dropped = [], []
    for path in paths:
        name = os.path.basename(path)
        print(f'{name}', end=' ', flush=True)
        try:
            M, info = maxima(path, ch, prominence, period)
        except Exception as e:
            print(f'-> ERROR: {e}')
            dropped.append((name, f'error: {e}'))
            continue
        if len(M) < 2:
            print(f'-> skipped ({info["status"]})')
            dropped.append((name, info['status']))
            continue
        r = rpot_from_name(name)
        print(f'-> R={r:7.1f} ohm  {len(M):5d} maxima')
        records.append(Record(name, r, 0.0, M, info))
    return records, dropped


def main():
    p = argparse.ArgumentParser()
    p.add_argument('inputs', nargs='*',
                   help='sweep folders and/or single CSV records named with '
                        'their resistance (omit to pick folders)')
    p.add_argument('-o', '--out',
                   help='output PNG (default: <first input>_lyapunov.png)')
    add_record_args(p)
    p.add_argument('--width', type=float, default=0.1,
                   help='slope-fit window as a fraction of the branch span in M_n')
    p.add_argument('--min-width', type=float, default=0.02,
                   help='floor on the slope-fit window, V')
    p.add_argument('--min-points', type=int, default=10,
                   help='fewest points a fit window (or a branch) may hold')
    p.add_argument('--gap', type=float, default=0.05,
                   help='branch break: gap in M_n as a fraction of its span')
    p.add_argument('--min-spread', type=float, default=0.1,
                   help='skip records whose map spans less than this, V; a '
                        'branch narrower than this is a cluster and is not fitted')
    p.add_argument('--min-curve-frac', type=float, default=0.5,
                   help='refuse a record unless this fraction of its points '
                        'lie on branches wide enough to fit')
    p.add_argument('--min-r2', type=float, default=0.8,
                   help='skip records whose local fits explain less than this'
                        ' fraction of the WITHIN-branch variance of M_(n+1)')
    p.add_argument('--min-steps', type=float, default=30,
                   help='refuse a map narrower than this many quantisation '
                        'steps of the channel')
    p.add_argument('--rosenstein', action='store_true',
                   help='also compute the direct time-series estimate for '
                        'every record that has a map lambda')
    p.add_argument('--fit', nargs=2, type=float, default=(0.5, 2.5),
                   metavar=('LO', 'HI'),
                   help='Rosenstein fit window, in periods')
    p.add_argument('--each', action='store_true',
                   help='also write a per-record diagnostic PNG')
    a = p.parse_args()

    files = [f for f in a.inputs if f.lower().endswith('.csv')]
    folder_args = [f for f in a.inputs if f not in files]
    for f in files:
        if not os.path.isfile(f):
            raise SystemExit(f'not a file: {f}')
    folders = resolve_folders(folder_args) if (folder_args or not files) else []
    first = folders[0] if folders else os.path.splitext(files[0])[0]
    out = a.out or sibling(first, '_lyapunov.png')
    csv_out = os.path.splitext(out)[0] + '.csv'
    each_dir = os.path.splitext(out)[0] + '_each'
    if a.each:
        os.makedirs(each_dir, exist_ok=True)

    groups = [(os.path.basename(f), 'folder', f) for f in folders]
    if files:
        groups.append(('records', 'files', files))
    colors = sweep_colors(len(groups))

    fig, ax = plt.subplots(figsize=(10, 5.5))
    rows, plotted = [], False

    for k, (label, kind, src) in enumerate(groups):
        print(f'\n--- {label} ---')
        if kind == 'folder':
            records, dropped = collect(src, a.ch, a.prominence,
                                       a.period_samples, a.max_residual,
                                       a.max_clip, a.max_files)
        else:
            records, dropped = single_records(src, a.ch, a.prominence,
                                              a.period_samples)
        if dropped:
            print(f'  dropped {len(dropped)} before fitting')

        xs, ys, es, xr, yr, er = [], [], [], [], [], []
        for rec in records:
            res = lyapunov(rec.M, rec.info['t_peaks'], a.width, a.min_width,
                           a.gap, a.min_spread, a.min_r2, a.min_points,
                           a.min_curve_frac, rec.info['quantum'], a.min_steps)
            ros = None
            if a.rosenstein:
                # The direct estimate needs no map, so it runs for every
                # admitted record, including the ones the map refuses.
                path = (os.path.join(src, rec.name) if kind == 'folder'
                        else next(f for f in src
                                  if os.path.basename(f) == rec.name))
                try:
                    t, v1, v2 = load_channels(path)
                    ros = rosenstein(v1, v2, rec.info['dt'],
                                     rec.info['period_samples'],
                                     fit=tuple(a.fit))
                    res['lam_ros'], res['lam_ros_err'] = ros['lam'], ros['lam_err']
                except Exception as e:
                    print(f'  {rec.name}: rosenstein failed: {e}')
                if ros and np.isfinite(rec.rpot):
                    xr.append(rec.rpot)
                    yr.append(ros['lam'])
                    er.append(ros['lam_err'])
            if res['status'] != 'ok':
                print(f'  {rec.name}: R={rec.rpot:7.1f}  {res["status"]}'
                      + (f'  direct {res["lam_ros"]:+7.0f} +/- '
                         f'{res["lam_ros_err"]:4.0f} /s' if ros else ''))
            else:
                print(f'  {rec.name}: R={rec.rpot:7.1f}  '
                      f'<T>={res["mean_T"] * 1e6:7.2f} us  '
                      f'(median {res["median_T"] * 1e6:7.2f}, '
                      f'{res["T_diff_pct"]:+5.1f} %)  '
                      f'<ln|f\'|>={res["mean_ln_slope"]:+6.3f}  '
                      f'R2={res["r2"]:.3f}  '
                      f'below noise {res["unresolved_pct"]:3.0f} %  '
                      f'map {res["lam"]:+7.0f} +/- {res["lam_err"]:4.0f} /s'
                      + (f'  direct {res["lam_ros"]:+7.0f} +/- '
                         f'{res["lam_ros_err"]:4.0f} /s' if ros else ''))
                if np.isfinite(rec.rpot):
                    xs.append(rec.rpot)
                    ys.append(res['lam'])
                    es.append(res['lam_err'])
            if a.each and (res['n_branches'] or ros):
                diagnostic(
                    os.path.join(each_dir,
                                 os.path.splitext(rec.name)[0] + '.png'),
                    rec.M, rec.info['t_peaks'], res,
                    f'{label} / {rec.name}   '
                    f'$R_{{pot}}$ = {rec.rpot:.1f} $\\Omega$   '
                    + (f'map $\\lambda$ = {res["lam"]:.0f} $\\pm$ '
                       f'{res["lam_err"]:.0f} s$^{{-1}}$'
                       if res['status'] == 'ok' else 'no map $\\lambda$'),
                    a.gap, ros,
                    rec.info['period_samples'] * rec.info['dt'])
            rows.append(dict(res, sweep=label, filename=rec.name,
                             rpot=rec.rpot,
                             mean_T_us=res['mean_T'] * 1e6,
                             median_T_us=res['median_T'] * 1e6))

        if xs:
            o = np.argsort(xs)
            xs, ys, es = (np.asarray(v)[o] for v in (xs, ys, es))
            ax.errorbar(xs, ys, yerr=es, fmt='.-', ms=5, lw=0.8,
                        elinewidth=0.6, capsize=2, color=colors[k],
                        label=f'{label}: return map  ({len(xs)} records)')
            plotted = True
            n_rec = sum(1 for r in rows if r['sweep'] == label)
            print(f'{label}: map lambda for {len(xs)}/{n_rec} records, '
                  f'{ys.min():.0f} to {ys.max():.0f} /s; the rest have no '
                  f'measurable slope (periodic or not a curve)')
            regime_report(label, rows)
        if xr:
            o = np.argsort(xr)
            xr, yr, er = (np.asarray(v)[o] for v in (xr, yr, er))
            ax.errorbar(xr, yr, yerr=er, fmt='o--', ms=4, lw=0.8,
                        elinewidth=0.6, capsize=2, color=colors[k],
                        mfc='none', label=f'{label}: direct (Rosenstein, '
                                          f'{len(xr)} records)')
            plotted = True
            print(f'{label}: direct lambda for {len(xr)} records, '
                  f'{yr.min():.0f} to {yr.max():.0f} /s')

    with open(csv_out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow([head for head, _, _ in CSV_FIELDS])
        w.writerows([[num(r[key], nd) for _, key, nd in CSV_FIELDS]
                     for r in rows])

    if not plotted:
        raise SystemExit(f'no record yielded a lambda\nsummary -> {csv_out}')

    ax.axhline(0, color='0.4', lw=0.8, ls='--')
    ax.set_xlabel('$R_{pot}$  ($\\Omega$)')
    ax.set_ylabel('$\\lambda$  (s$^{-1}$)')
    ax.set_title('Largest Lyapunov exponent\n'
                 "return map: $\\lambda \\simeq \\langle \\ln|f'(M_n)| \\rangle"
                 " / \\langle T \\rangle$   (records whose map is not a curve"
                 ' are left out)', fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.15, lw=0.5)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)

    print(f'\nplot    -> {out}\nsummary -> {csv_out}')
    print('the map value is an estimate from the return map and reads high '
          'where the map\nis not single-valued (see benchmark_lyapunov.py); '
          'periodic windows carry no lambda by design.')


if __name__ == '__main__':
    main()
