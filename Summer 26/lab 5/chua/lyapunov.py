"""
Lyapunov exponent estimated from the first-return map.

    lambda ~= < ln |f'(M_n)| > / <T>

f'(M_n) is the local slope of the return map, taken by a sliding linear fit
along M_n, and <T> is the mean return time between successive maxima.

The two traps this script is written around:

  <T> is the MEAN, not the median. The return-time distribution is skewed
  because a lobe switch takes about twice as long as an ordinary winding, so
  the median sits below the mean. Both are reported per record, with the
  difference split by regime, so the number can be quoted and the choice
  defended: a single scroll has no lobe switches and shows no skew, while the
  double-scroll records do.

  ln |f'| is averaged over the points ACTUALLY VISITED, not from one
  straight-line fit per branch. The slope is re-fitted in a window sliding
  along M_n, so the neighbourhood of the turning point - where the slope falls
  toward zero and ln |f'| dives - is weighted by how often the orbit really
  goes there. A single fit per branch misses this and reads high.

Branches are separated before fitting (--gap): with a double scroll the map
falls into two clusters, and a window straddling the gap between them would
measure the empty space rather than a slope.

The formula assumes M_(n+1) really is a function of M_n. On the transitional
records here it is not - the maxima form a blob, the local fits explain almost
none of it, and the "slope" that comes out is fitted noise. So the fit quality
is measured (R2, reported per record) and a record below --min-r2 is left
without a lambda rather than given a meaningless one.

This is an ESTIMATE FROM THE RETURN MAP, not a direct measurement. A proper
Rosenstein or Wolf calculation on the time series typically comes out 10 % to
20 % lower, and the number should be quoted that way.

Records are chosen by lorenz_map.collect, the same gate the bifurcation
diagram uses, so the two plots cover the same sweep.

Usage:
    python lyapunov.py                             # pick folders, repeatedly
    python lyapunov.py "sweep forward"
    python lyapunov.py "sweep forward" "sweep back"
    python lyapunov.py "sweep forward" --each      # + per-record diagnostics
    python lyapunov.py "sweep forward" --window 31
"""
import argparse
import csv
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy.lib.stride_tricks import sliding_window_view

from bifurcation import add_record_args
from lorenz_map import collect, frame
from sweeplib import resolve_folders, sibling, sweep_colors

CSV_FIELDS = [
    ('sweep', 'sweep', None), ('filename', 'filename', None),
    ('rpot_ohm', 'rpot', 2), ('n_maxima', 'n_maxima', None),
    ('n_branches', 'n_branches', None), ('n_slopes_fitted', 'n_fitted', None),
    ('n_slopes_dropped', 'n_dropped', None), ('mean_T_us', 'mean_T_us', 4),
    ('median_T_us', 'median_T_us', 4),
    ('mean_vs_median_pct', 'T_diff_pct', 2),
    ('mean_ln_abs_fprime', 'mean_ln_slope', 5),
    ('lambda_per_s', 'lam', 2), ('map_r2', 'r2', 4),
    ('status', 'status', None),
]


def split_branches(m_n, gap):
    """
    Indices of `m_n` sorted by value, cut into branches at large gaps.

    A double-scroll map lives on two separated clusters; a sliding fit that
    crosses the void between them measures the void. Each branch comes back
    already in ascending m_n order, so callers need not re-sort.
    """
    order = np.argsort(m_n, kind='stable')
    v = m_n[order]
    span = float(v[-1] - v[0])
    if span <= 0:
        return [order]
    cuts = np.flatnonzero(np.diff(v) > gap * span) + 1
    return np.split(order, cuts) if len(cuts) else [order]


def sliding_fit(x, y, k):
    """
    Local slope and fitted value at every point, from a length-k OLS window.

    Closed-form least squares over sliding windows: for the interior each
    point takes the window centred on it, and near the ends the nearest full
    window, so every visited point gets a slope without shrinking the fit at
    the edges. The fitted value comes back too, so how well the map behaves
    as a function can be measured from the same fits.
    """
    n = len(x)
    if n < k:
        return np.full(n, np.nan), np.full(n, np.nan)
    X = sliding_window_view(x, k)
    Y = sliding_window_view(y, k)
    Sx = X.sum(1)
    Sy = Y.sum(1)
    den = k * (X * X).sum(1) - Sx * Sx
    num = k * (X * Y).sum(1) - Sx * Sy
    with np.errstate(divide='ignore', invalid='ignore'):
        slope = np.where(np.abs(den) > 0, num / den, np.nan)
    # Point i -> the window whose centre is nearest to i.
    idx = np.clip(np.arange(n) - k // 2, 0, n - k)
    s = slope[idx]
    pred = Sy[idx] / k + s * (x - Sx[idx] / k)
    return s, pred


def lyapunov(M, t_peaks, window=21, gap=0.05, min_spread=0.1, min_r2=0.8):
    """
    Return-time statistics, <ln|f'|> and lambda for one record.

    M is the sequence of maxima and t_peaks the time of each. The per-point
    slopes are returned under 'slopes' (indexed like M[:-1]) so the diagnostic
    plot plots exactly what was measured instead of re-deriving it.
    """
    out = {'n_maxima': len(M), 'status': 'ok', 'n_branches': 0,
           'mean_T': np.nan, 'median_T': np.nan, 'T_diff_pct': np.nan,
           'mean_ln_slope': np.nan, 'lam': np.nan, 'n_fitted': 0,
           'n_dropped': 0, 'spread': 0.0, 'r2': np.nan,
           'slopes': np.full(max(len(M) - 1, 0), np.nan)}
    if len(M) < window + 2 or len(t_peaks) != len(M):
        out['status'] = 'too few maxima for a sliding fit'
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

    m_n, m_next = M[:-1], M[1:]
    spread = float(m_n.max() - m_n.min())
    out['spread'] = spread
    if spread < min_spread:
        # A periodic orbit collapses the map onto a few points; there is no
        # slope to measure, and any number here would be fitting noise.
        out['status'] = f'map spread {spread:.3f} V below --min-spread'
        return out

    slopes = out['slopes']
    fitted = []
    for idx in split_branches(m_n, gap):
        if len(idx) < window:
            continue
        s, pred = sliding_fit(m_n[idx], m_next[idx], window)
        slopes[idx] = s
        fitted.append((m_next[idx], pred))
    out['n_branches'] = len(fitted)
    if not fitted:
        out['status'] = f'no branch has {window} points'
        return out

    # How much of M_{n+1} the local fits actually explain. lambda from a
    # return map assumes M_{n+1} IS a function of M_n; when the points are a
    # blob instead of a curve the "slope" is fitted noise, and on the
    # transitional records here that reads as a wildly large exponent.
    obs = np.concatenate([o for o, _ in fitted])
    fit = np.concatenate([f for _, f in fitted])
    ok_fit = np.isfinite(obs) & np.isfinite(fit)
    ss_tot = float(np.sum((obs[ok_fit] - obs[ok_fit].mean()) ** 2))
    ss_res = float(np.sum((obs[ok_fit] - fit[ok_fit]) ** 2))
    out['r2'] = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    if np.isfinite(out['r2']) and out['r2'] < min_r2:
        out['status'] = (f'return map is not a function '
                         f'(local fit R2 = {out["r2"]:.2f})')
        return out

    good = np.isfinite(slopes) & (np.abs(slopes) > 0)
    out['n_fitted'] = int(good.sum())
    out['n_dropped'] = int((~good).sum())
    if out['n_fitted'] < window:
        out['status'] = 'too few usable slopes'
        return out

    mean_ln = float(np.mean(np.log(np.abs(slopes[good]))))
    out['mean_ln_slope'] = mean_ln
    out['lam'] = mean_ln / mean_T
    return out


def diagnostic(path, M, t_peaks, res, title):
    """Return map with its local slopes, and the return-time histogram."""
    m_n, m_next = M[:-1], M[1:]
    slopes = res['slopes']
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4))

    good = np.isfinite(slopes)
    sc = a1.scatter(m_n[good], m_next[good], c=np.log(np.abs(slopes[good])),
                    s=4, cmap='coolwarm', lw=0)
    fig.colorbar(sc, ax=a1, label="$\\ln|f'(M_n)|$")
    lo = float(min(m_n.min(), m_next.min()))
    hi = float(max(m_n.max(), m_next.max()))
    pad = 0.06 * (hi - lo)
    frame(a1, (lo - pad, hi + pad))
    a1.set_xlabel('$M_n$  (V)')
    a1.set_ylabel('$M_{n+1}$  (V)')
    a1.set_title('return map, coloured by local slope', fontsize=9)

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


def main():
    p = argparse.ArgumentParser()
    p.add_argument('folders', nargs='*',
                   help='one or more sweep folders (omit to pick them)')
    p.add_argument('-o', '--out',
                   help='output PNG (default: <first folder>_lyapunov.png)')
    add_record_args(p)
    p.add_argument('--window', type=int, default=21,
                   help='points per sliding linear fit')
    p.add_argument('--gap', type=float, default=0.05,
                   help='branch break: gap in M_n as a fraction of its span')
    p.add_argument('--min-spread', type=float, default=0.1,
                   help='skip records whose map spans less than this, V')
    p.add_argument('--min-r2', type=float, default=0.8,
                   help='skip records whose local fits explain less than this'
                        ' fraction of M_(n+1); below it the map is not a curve')
    p.add_argument('--each', action='store_true',
                   help='also write a per-record diagnostic PNG')
    a = p.parse_args()

    folders = resolve_folders(a.folders)
    out = a.out or sibling(folders[0], '_lyapunov.png')
    csv_out = os.path.splitext(out)[0] + '.csv'
    colors = sweep_colors(len(folders))
    each_dir = os.path.splitext(out)[0] + '_each'
    if a.each:
        os.makedirs(each_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    rows, plotted = [], False

    for k, folder in enumerate(folders):
        label = os.path.basename(folder)
        print(f'\n--- {label} ---')
        records, dropped = collect(folder, a.ch, a.prominence, a.period_samples,
                                   a.max_residual, a.max_clip, a.max_files)
        if dropped:
            print(f'  dropped {len(dropped)} before fitting')

        xs, ys = [], []
        for rec in records:
            res = lyapunov(rec.M, rec.info['t_peaks'], a.window, a.gap,
                           a.min_spread, a.min_r2)
            if res['status'] != 'ok':
                print(f'  {rec.name}: {res["status"]}')
            else:
                print(f'  {rec.name}: R={rec.rpot:7.1f}  '
                      f'<T>={res["mean_T"] * 1e6:7.2f} us  '
                      f'(median {res["median_T"] * 1e6:7.2f}, '
                      f'{res["T_diff_pct"]:+5.1f} %)  '
                      f'<ln|f\'|>={res["mean_ln_slope"]:+6.3f}  '
                      f'R2={res["r2"]:.3f}  lambda={res["lam"]:+8.1f} /s')
                xs.append(rec.rpot)
                ys.append(res['lam'])
                if a.each:
                    diagnostic(
                        os.path.join(each_dir,
                                     os.path.splitext(rec.name)[0] + '.png'),
                        rec.M, rec.info['t_peaks'], res,
                        f'{label} / {rec.name}   '
                        f'$R_{{pot}}$ = {rec.rpot:.1f} $\\Omega$   '
                        f'$\\lambda$ = {res["lam"]:.0f} s$^{{-1}}$')
            rows.append(dict(res, sweep=label, filename=rec.name,
                             rpot=rec.rpot,
                             mean_T_us=res['mean_T'] * 1e6,
                             median_T_us=res['median_T'] * 1e6))

        if xs:
            o = np.argsort(xs)
            xs = np.asarray(xs)[o]
            ys = np.asarray(ys)[o]
            ax.plot(xs, ys, '.-', ms=5, lw=0.8, color=colors[k],
                    label=f'{label}  ({len(xs)} records)')
            plotted = True
            print(f'{label}: lambda from {ys.min():.0f} to {ys.max():.0f} /s, '
                  f'{int(np.sum(ys > 0))}/{len(ys)} positive')
            regime_report(label, rows)

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
    ax.set_title('Lyapunov exponent estimated from the return map\n'
                 "$\\lambda \\simeq \\langle \\ln|f'(M_n)| \\rangle"
                 " / \\langle T \\rangle$", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.15, lw=0.5)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)

    print(f'\nplot    -> {out}\nsummary -> {csv_out}')
    print('quote as an estimate from the return map: a Rosenstein or Wolf\n'
          'calculation on the time series typically lands 10-20 % lower.')


if __name__ == '__main__':
    main()
