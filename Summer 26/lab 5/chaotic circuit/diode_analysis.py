"""
M1/M2: the five-segment fit of the nonlinear element's V-I record.

In the V-I record CH1 = VV (the source side of the shunt) and CH2 = VI (across
the element); the current is i = (VV - VI)/Rs and it is fitted against VI, as
the experiment plan specifies. (The oscillator records use other channels.)

Outputs, in --out-dir (default: this folder):
    diode_fit.json            every fitted number, read by chua/simulate.py
    diode_segments.csv        the five segments
    diode_fit_stats.txt       fit uncertainties and chi-square tests
    current_vs_voltage.pdf/.png   report figure
    optimized_breakpoints.png     diagnostic: samples, segments and residuals
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, find_peaks
from scipy.stats import chi2 as chi2_dist

from chua.scope_data import R0
import report_style

HERE = Path(__file__).resolve().parent
RS = 216.0                       # ohm, the shunt
SOURCE_RANGE = (-9.33, 8.28)     # V, the source swing without its end-of-drive saturation
LABELS = ['Gc left', 'Gb left', 'Ga inner', 'Gb right', 'Gc right']


def fit_segments(x, y, n_segments=5, min_points=50, min_width=0.3):
    """
    Segmented least squares with the global minimum of the total SSE.

    Every sample belongs to exactly one segment, and segments end between
    observed voltage codes (the breakpoint is the midpoint of the two codes).
    Dynamic programming over the codes finds the optimum without any
    stochastic search. Returns (one dict per segment, combined R^2).
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.ndim != 1 or y.shape != x.shape or not np.isfinite([x, y]).all():
        raise ValueError('x and y must be matching finite vectors')
    order = np.argsort(x, kind='stable')
    x, y = x[order], y[order]
    codes, starts = np.unique(x, return_index=True)
    cuts = np.r_[starts, len(x)]
    m = len(codes)
    if m > 2500:
        raise ValueError('more than 2500 voltage codes; bin higher-resolution data first')
    if m < 2 * n_segments:
        raise ValueError('too few distinct voltages for the requested segments')

    # costs[a, b]: SSE of one line through the codes a..b-1, from prefix sums
    prefixes = [np.r_[0., np.cumsum(v)] for v in (x, y, x * x, x * y, y * y)]
    costs = np.full((m + 1, m + 1), np.inf)
    for a in range(m):
        b = np.arange(a + 2, m + 1)
        n = cuts[b] - cuts[a]
        sx, sy, sxx, sxy, syy = [v[cuts[b]] - v[cuts[a]] for v in prefixes]
        vx = sxx - sx * sx / n
        valid = (n >= min_points) & (codes[b - 1] - codes[a] >= min_width) & (vx > 0)
        sse = syy - sy * sy / n - (sxy - sx * sy / n) ** 2 / np.maximum(vx, 1e-300)
        costs[a, b[valid]] = np.maximum(sse[valid], 0)

    # best[k, b]: least SSE of k segments covering codes 0..b-1
    best = np.full((n_segments + 1, m + 1), np.inf)
    prev = np.full(best.shape, -1, int)
    best[0, 0] = 0
    for k in range(1, n_segments + 1):
        for b in range(1, m + 1):
            vals = best[k - 1, :b] + costs[:b, b]
            a = int(np.argmin(vals))
            best[k, b], prev[k, b] = vals[a], a
    if not np.isfinite(best[-1, -1]):
        raise ValueError('no partition satisfies the minimum width/sample count')
    ends, b = [m], m
    for k in range(n_segments, 0, -1):
        b = prev[k, b]
        ends.append(b)
    ends = list(reversed(ends))

    bounds = [float(x[0])] + [float((codes[e - 1] + codes[e]) / 2) for e in ends[1:-1]] + [float(x[-1])]
    fits = []
    for k, (a, b) in enumerate(zip(ends[:-1], ends[1:])):
        xx, yy = x[cuts[a]:cuts[b]], y[cuts[a]:cuts[b]]
        (slope, intercept), cov = np.polyfit(xx, yy, 1, cov=True)
        residual = yy - (slope * xx + intercept)
        sst = np.sum((yy - yy.mean()) ** 2)
        fits.append(dict(label=LABELS[k] if n_segments == 5 else str(k + 1),
                         v_lo=bounds[k], v_hi=bounds[k + 1], n=len(xx),
                         slope_S=float(slope), intercept_A=float(intercept),
                         slope_se_S=float(np.sqrt(cov[0, 0])),
                         intercept_se_A=float(np.sqrt(cov[1, 1])),
                         r2=float(1 - np.sum(residual ** 2) / sst) if sst else None))
    sst = float(np.sum((y - y.mean()) ** 2))
    return fits, float(1 - best[-1, -1] / sst) if sst else None


def slope_ranges(fits, r0=R0):
    """
    M2: the resistance range over which each shoulder can oscillate,
    -1/Ga < Rt < -1/Gb, total and as a potentiometer setting (slopes only).
    """
    ga = fits[2]['slope_S']
    out = []
    for side, i in [('left', 1), ('right', 3)]:
        gb = fits[i]['slope_S']
        valid = ga < gb < 0
        low, high = (-1 / ga, -1 / gb) if valid else (None, None)
        out.append(dict(side=side, valid=valid, total_lo_ohm=low, total_hi_ohm=high,
                        pot_lo_ohm=low - r0 if valid else None,
                        pot_hi_ohm=high - r0 if valid else None))
    return out


def code_means(vi, current, keep, segments, u_v):
    """
    One measurement per voltage code: (codes, mean current, its uncertainty, segment index).

    sigma^2 = (std/sqrt(n))^2 + (G u_v)^2: the standard error of the code's mean
    current plus the code's own voltage quantisation projected through the local
    slope. The per-sample current quantisation (~0.24 mA) is the resolution of a
    single reading, not the uncertainty of a mean over 100-700 samples: the source
    channel steps through its own codes within one element code and dithers it.
    """
    table = pd.DataFrame({'v': vi[keep], 'i': current[keep]}).groupby('v').i.agg(['mean', 'std', 'count'])
    x = table.index.to_numpy()
    seg = np.full(len(x), -1)
    for k, f in enumerate(segments):
        last = k == len(segments) - 1
        seg[(x >= f['v_lo']) & ((x <= f['v_hi']) if last else (x < f['v_hi']))] = k
    slope = np.array([segments[k]['slope_S'] if k >= 0 else 0.0 for k in seg])
    se = (table['std'].fillna(0) / np.sqrt(table['count'])).to_numpy()
    sigma = np.sqrt(se ** 2 + (slope * u_v) ** 2)
    return x, table['mean'].to_numpy(), sigma, seg


def chi2_block(c2, dof):
    """chi2, dof, chi2/dof and p, as stored in the JSON."""
    return dict(chi2=float(c2), dof=int(dof), chi2_red=float(c2 / dof) if dof > 0 else np.nan,
                p_value=float(chi2_dist.sf(c2, dof)) if dof > 0 else np.nan)


def weighted_line(x, y, sigma):
    """(parameters, covariance) of the weighted least-squares line y = slope x + intercept."""
    X = np.c_[x, np.ones(len(x))]
    w = 1 / sigma ** 2
    cov = np.linalg.inv(X.T @ (X * w[:, None]))
    return cov @ (X.T @ (w * y)), cov


def code_stats(vi, current, keep, segments, u_v, overlap):
    """
    Fit uncertainties and chi-square tests on the per-code measurements.

    For each segment: the a-priori parameter covariance (X^T W X)^-1 with the
    per-code sigma, chi2 of the reported (OLS) line against the code means, the
    same chi2 without the code at the drive turnaround (first code of the first
    segment, last code of the last: the source clips there), and the number of
    sign runs in the residuals against its random expectation n/2 + 1 (a run
    count far below that means systematic curvature, not noise).
    """
    x_all, y_all, s_all, seg = code_means(vi, current, keep, segments, u_v)
    out = []
    c2_all = c2_core = 0.0
    dof_all = dof_core = 0
    for k, f in enumerate(segments):
        m = seg == k
        x, y, s = x_all[m], y_all[m], s_all[m]
        z = (y - (f['slope_S'] * x + f['intercept_A'])) / s
        full = chi2_block(np.sum(z ** 2), len(z) - 2)
        core_mask = np.ones(len(z), bool)
        if k == 0:
            core_mask[0] = False
        if k == len(segments) - 1:
            core_mask[-1] = False
        core = chi2_block(np.sum(z[core_mask] ** 2), core_mask.sum() - 2)
        X = np.c_[x[core_mask], np.ones(core_mask.sum())]     # same points as the core chi2
        cov = np.linalg.inv(X.T @ (X / s[core_mask, None] ** 2))
        signs = np.sign(z)
        runs = int(1 + np.sum(signs[1:] != signs[:-1]))
        # local slope over thirds of the segment: how much the "straight" segment actually bends.
        # The statistical u is not inflated by sqrt(chi2/dof): the excess is curvature, not noise.
        xc, yc, sc = x[core_mask], y[core_mask], s[core_mask]
        thirds = []
        for part in np.array_split(np.arange(len(xc)), 3):
            p, c = weighted_line(xc[part], yc[part], sc[part])
            thirds.append(dict(v_lo=float(xc[part][0]), v_hi=float(xc[part][-1]),
                               slope_S=float(p[0]), u_slope_S=float(np.sqrt(c[0, 0]))))
        c2_all += full['chi2']
        dof_all += full['dof']
        c2_core += core['chi2']
        dof_core += core['dof']
        out.append(dict(label=f['label'], n_codes=int(m.sum()), sigma_median_A=float(np.median(s)),
                        u_slope_S=float(np.sqrt(cov[0, 0])), u_intercept_A=float(np.sqrt(cov[1, 1])),
                        chi2_all=full, chi2_core=core, sign_runs=runs,
                        sign_runs_expected=float(len(z) / 2 + 1), max_abs_z=float(np.max(np.abs(z))),
                        local_slopes=thirds))
    direction = None
    if len(overlap):
        d, s = overlap[:, 1] - overlap[:, 2], overlap[:, 3]
        direction = chi2_block(np.sum((d / s) ** 2), len(d))
    return dict(segments=out, chi2_all=chi2_block(c2_all, dof_all), chi2_core=chi2_block(c2_core, dof_core),
                direction_check=direction, u_v_V=float(u_v))


def direction_overlap(vi, current, keep, rising, falling, n_bins=80):
    """
    Rows (bin centre, mean current rising, mean current falling, SE of the
    difference) over element-voltage bins where both sweep directions have
    at least ten samples.
    """
    overlap = []
    edges = np.linspace(vi[keep].min(), vi[keep].max(), n_bins + 1)
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = keep & (vi >= lo) & (vi < hi)
        up, down = current[mask & rising], current[mask & falling]
        if len(up) >= 10 and len(down) >= 10:
            overlap.append(((lo + hi) / 2, up.mean(), down.mean(),
                            np.sqrt(up.var(ddof=1) / len(up) + down.var(ddof=1) / len(down))))
    return np.asarray(overlap).reshape(-1, 4)


def analyze(path, source='CH1(V)', element='CH2(V)', rs=RS, source_range=SOURCE_RANGE):
    """Fit the record; returns (result dict, arrays for the figures)."""
    if not np.isfinite(rs) or rs <= 0:
        raise ValueError('shunt resistance must be positive')
    data = pd.read_csv(path)
    t, vv, vi = (data[c].to_numpy(float) for c in ['Time(s)', source, element])
    if source == element or not np.isfinite([t, vv, vi]).all() or np.any(np.diff(t) <= 0):
        raise ValueError('invalid channels or timestamps')
    current = (vv - vi) / rs
    # drop the drive-end saturation of the source: a source-voltage limit, not a chosen segment boundary
    keep = (vv >= source_range[0]) & (vv <= source_range[1])
    fits, r2 = fit_segments(vi[keep], current[keep])

    # sweep direction from the smoothed source voltage, for the hysteresis check
    window = min(501, (len(vv) - 1) // 2 * 2 + 1)
    smoothed = savgol_filter(vv, window, 3)
    slope = np.gradient(smoothed, t)
    threshold = 0.1 * np.median(np.abs(slope))
    rising, falling = slope > threshold, slope < -threshold
    overlap = direction_overlap(vi, current, keep, rising, falling)

    peaks, _ = find_peaks(smoothed, prominence=0.5 * np.ptp(smoothed))
    troughs, _ = find_peaks(-smoothed, prominence=0.5 * np.ptp(smoothed))
    turns = np.sort(np.r_[peaks, troughs])
    freq = float(1 / (2 * np.median(np.diff(t[turns])))) if len(turns) > 1 else None
    # uniform-within-code quantisation, q/sqrt(12) per channel, both channels independent
    u_v = [float(np.median(np.diff(np.unique(ch))) / np.sqrt(12)) for ch in (vv, vi)]
    u_current = float(np.hypot(*u_v) / rs)
    directional_rms = float(np.sqrt(np.mean((overlap[:, 1] - overlap[:, 2]) ** 2))) if len(overlap) else None
    result = dict(schema_version=1, input=str(Path(path).resolve()),
                  voltage_column=element, source_column=source, voltage_basis='element',
                  current_definition=f'({source} - {element}) / {rs:g} ohm',
                  shunt_ohm=rs, r0_ohm=R0, source_fit_range_V=list(source_range),
                  n_input=len(t), n_fit=int(keep.sum()), combined_r2=r2,
                  segments=fits, ranges=slope_ranges(fits), drive_frequency_hz=freq,
                  directional_rms_difference_A=directional_rms,
                  u_voltage_quantization_V=u_v[1], u_current_quantization_A=u_current,
                  code_stats=code_stats(vi, current, keep, fits, u_v[1], overlap),
                  uncertainty_note='Segment slope_se/intercept_se are OLS standard errors over samples; '
                                   'code_stats carries the per-code error model (SE of the mean plus projected '
                                   'voltage quantisation) with chi-square and p-values. Channel gains, resistor '
                                   'calibration and drift are not included in either.')
    return result, dict(vi=vi, current=current, keep=keep)


def plot_report_figure(result, arrays, out):
    """current_vs_voltage.pdf/.png: one point per voltage code with its sigma, and the five lines."""
    vi, current, keep = arrays['vi'], arrays['current'], arrays['keep']
    with plt.rc_context():
        report_style.apply()
        fig, ax = plt.subplots(figsize=report_style.figsize(aspect=0.8), layout='constrained')
        c = report_style.COLORS
        # both sweep directions pooled (they agree to directional_rms_difference_A, far below one
        # current code); the bars are visible only on the steep outer segments
        x, y, s, _ = code_means(vi, current, keep, result['segments'], result['u_voltage_quantization_V'])
        ax.errorbar(x, y * 1e3, yerr=s * 1e3, ls='none', marker='o', ms=1.8, mew=0, color=c[0], ecolor=c[0],
                    elinewidth=0.5, capsize=0, label='Measured', zorder=2)
        for k, f in enumerate(result['segments']):
            xx = np.array([f['v_lo'], f['v_hi']])
            ax.plot(xx, (f['slope_S'] * xx + f['intercept_A']) * 1e3, color='black', lw=0.9, zorder=3,
                    label='Piecewise-linear fit' if k == 0 else None)
        ax.set(xlabel='Element voltage (V)', ylabel='Element current (mA)', xlim=(-9.2, 8.6))
        ax.set_xticks(np.arange(-8, 9, 4))
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='upper right', bbox_to_anchor=(0.9, 1))
        report_style.save(fig, out / 'current_vs_voltage.png')
        plt.close(fig)


def plot_diagnostic(result, arrays, out):
    """optimized_breakpoints.png: every sample, the segments with their breakpoints, and the residuals."""
    vi, current, keep = arrays['vi'], arrays['current'], arrays['keep']
    fig, (ax, resax) = plt.subplots(2, 1, figsize=(10, 7.8), height_ratios=[3, 1], layout='constrained',
                                    sharex=True)
    ax.scatter(vi[keep], current[keep] * 1e3, s=2, alpha=.12, color='#0072B2', rasterized=True,
               label='Measured samples')
    for k, f in enumerate(result['segments']):
        xx = np.array([f['v_lo'], f['v_hi']])
        color = plt.get_cmap('tab10')(k)
        ax.plot(xx, (f['slope_S'] * xx + f['intercept_A']) * 1e3, color=color, lw=2,
                label=f"{f['label']}: {f['slope_S'] * 1e3:.3f} mS")
        mask = keep & (vi >= xx[0]) & ((vi < xx[1]) if k < 4 else (vi <= xx[1]))
        resax.scatter(vi[mask], (current[mask] - f['slope_S'] * vi[mask] - f['intercept_A']) * 1e3,
                      s=2, alpha=.12, color=color, rasterized=True)
        if k < 4:
            ax.axvline(xx[1], ls=':', color='0.4', lw=.8)
            ax.text(xx[1], .98, f'{xx[1]:.2f} V', rotation=90, va='top', ha='right',
                    transform=ax.get_xaxis_transform(), fontsize=8)
    ax.set(ylabel='Current (mA)',
           title=f"M1 · Optimized element-voltage breakpoints (R² = {result['combined_r2']:.4f})")
    ax.legend(fontsize=8, loc='lower right')
    ax.grid(alpha=.2)
    resax.axhline(0, color='0.4', lw=.8)
    resax.set(xlabel=f"Voltage across nonlinear element, {result['voltage_column'].split('(')[0]} (V)",
              ylabel='Fit residual\n(mA)')
    resax.grid(alpha=.2)
    fig.savefig(out / 'optimized_breakpoints.png', dpi=200)
    plt.close(fig)


def fit_stats_table(result, bootstrap_path=HERE / 'uncertainty' / 'diode_uncertainty.json'):
    """Per-segment fit statistics as plain text; the block-bootstrap u is added when available."""
    boot = {}
    if bootstrap_path.exists():
        boot = {s['label']: s for s in json.loads(bootstrap_path.read_text(encoding='utf-8'))['segments']}
    cs = result['code_stats']
    per_code = {s['label']: s for s in cs['segments']}
    lines = ['1. Segment fits (one measurement per voltage code; sigma = SE of the code mean (+) slope x u(V))',
             f"{'segment':10} {'V range (V)':>17} {'codes':>5} {'sigma':>6} {'slope (mS)':>11} {'u':>7}"
             f" {'intercept (mA)':>15} {'u':>7} {'R2':>6} | {'local slope, thirds (mS)':>26}"]
    lines.append('-' * len(lines[-1]))
    for f in result['segments']:
        s = per_code[f['label']]
        thirds = '  '.join(f"{1e3 * t['slope_S']:.3f}" for t in s['local_slopes'])
        lines.append(f"{f['label']:10} {f['v_lo']:8.3f}..{f['v_hi']:7.3f} {s['n_codes']:5d} {1e3 * s['sigma_median_A']:6.3f}"
                     f" {1e3 * f['slope_S']:11.4f} {1e3 * s['u_slope_S']:7.4f}"
                     f" {1e3 * f['intercept_A']:15.4f} {1e3 * s['u_intercept_A']:7.4f}"
                     f" {f['r2']:6.3f} | {thirds:>26}")
    lines += ['sigma: median per-code uncertainty (mA). u: statistical fit uncertainty from those sigmas',
              '(outer segments: without the turnaround code); NOT inflated by sqrt(chi2/dof), because the',
              'excess scatter on the shoulders is systematic curvature, shown by the local slopes over',
              'thirds of each segment (from the low-voltage to the high-voltage end, each +- ~0.01 mS on',
              'the shoulders). Outer-segment intercepts are the line extrapolated 7-8 V to V = 0 and are',
              'fully correlated with the slope; quote breakpoints instead.',
              '',
              'Breakpoints (midpoint between adjacent codes; true corner within one code, u = q/sqrt(12)):',
              '  ' + ', '.join(f"{f['v_hi']:+.3f}" for f in result['segments'][:-1])
              + f" V, each +- {1e3 * result['u_voltage_quantization_V']:.0f} mV",
              '',
              '2. Chi-square of the fitted lines against the code means',
              f"{'segment':10} {'chi2':>7} {'dof':>4} {'chi2/dof':>9} {'p':>8} | "
              f"{'core chi2':>9} {'dof':>4} {'chi2/dof':>9} {'p':>8}"
              f" | {'runs':>4} {'expect':>6} {'max|z|':>6}"]
    lines.append('-' * len(lines[-1]))
    for s in cs['segments']:
        a, c = s['chi2_all'], s['chi2_core']
        lines.append(f"{s['label']:10} {a['chi2']:7.1f} {a['dof']:4d} {a['chi2_red']:9.2f} {a['p_value']:8.3f} |"
                     f" {c['chi2']:9.1f} {c['dof']:4d} {c['chi2_red']:9.2f} {c['p_value']:8.3f} |"
                     f" {s['sign_runs']:4d} {s['sign_runs_expected']:6.1f} {s['max_abs_z']:6.1f}")
    a, c = cs['chi2_all'], cs['chi2_core']
    lines += [f"{'all five':10} {a['chi2']:7.1f} {a['dof']:4d} {a['chi2_red']:9.2f} {a['p_value']:8.3f} |"
              f" {c['chi2']:9.1f} {c['dof']:4d} {c['chi2_red']:9.2f} {c['p_value']:8.3f} |",
              'core: without the first and last code of the trace (drive turnaround, source clipped).',
              'runs: sign runs of the residuals; expect = n/2 + 1 for random scatter. Far fewer runs = curvature.',
              '',
              '3. Other slope uncertainty estimates (mS)',
              f"{'segment':10} {'OLS se (samples)':>17} {'block bootstrap':>16} {'N samples':>10}"]
    for f in result['segments']:
        b = boot.get(f['label'], {})
        bs = f"{1e3 * b['u_slope_block_S']:16.4f}" if 'u_slope_block_S' in b else f"{'-':>16}"
        lines.append(f"{f['label']:10} {1e3 * f['slope_se_S']:17.4f} {bs} {f['n']:10d}")
    d = cs['direction_check']
    direction = (f"chi2 = {d['chi2']:.1f}, dof = {d['dof']} bins, chi2/dof = {d['chi2_red']:.2f}, p = {d['p_value']:.3f}"
                 if d else 'no overlap bins')
    lines += ['',
              '4. Record',
              f"Points fitted: {result['n_fit']} of {result['n_input']}; drive frequency {result['drive_frequency_hz']:.2f} Hz",
              f"Single-reading resolution: u(V) = {1e3 * result['u_voltage_quantization_V']:.1f} mV, "
              f"u(I) = {1e3 * result['u_current_quantization_A']:.3f} mA (one code / sqrt(12); not the per-point bar)",
              f"Sweep-direction check: rising - falling = {1e3 * result['directional_rms_difference_A']:.4f} mA rms; "
              + direction,
              'Not included anywhere above: channel gain accuracy, shunt resistor calibration, drift.'
              if boot else 'block bootstrap: not available (run uncertainty/calculate.py).']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('csv', nargs='?', type=Path, default=HERE / 'trace1.csv')
    parser.add_argument('--source', default='CH1(V)')
    parser.add_argument('--element', default='CH2(V)')
    parser.add_argument('--rs', type=float, default=RS)
    parser.add_argument('--source-range', type=float, nargs=2, default=SOURCE_RANGE)
    parser.add_argument('--out-dir', type=Path, default=HERE)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    result, arrays = analyze(args.csv, args.source, args.element, args.rs, args.source_range)
    (args.out_dir / 'diode_fit.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    with (args.out_dir / 'diode_segments.csv').open('w', newline='') as fh:
        writer = csv.DictWriter(fh, list(result['segments'][0]))
        writer.writeheader()
        writer.writerows(result['segments'])
    plot_report_figure(result, arrays, args.out_dir)
    plot_diagnostic(result, arrays, args.out_dir)
    stats = fit_stats_table(result)
    (args.out_dir / 'diode_fit_stats.txt').write_text(stats + '\n', encoding='utf-8')

    print(f"M1: voltage={args.element}, source={args.source}\n")
    print(stats + '\n')
    for r in result['ranges']:
        print(f"M2 {r['side']}: {r}")
    print('M2 is a slope-only necessary estimate. Check intersections with the measured offset retained.')
    print(f"Sweep frequency: {result['drive_frequency_hz']} Hz; "
          f"rising/falling RMS difference: {result['directional_rms_difference_A']} A")


if __name__ == '__main__':
    main()
