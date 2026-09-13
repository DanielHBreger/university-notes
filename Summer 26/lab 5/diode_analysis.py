"""M1/M2 analysis. In the VI record CH1=VV (source), CH2=VI (element).

Fit i=(VV-VI)/Rs against VI, as specified in the experiment plan.
Oscillator-mode channel assignments are different from this VI recording.
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

HERE = Path(__file__).resolve().parent
from chua.scope_data import R0
RS = 216.0
LABELS = ['Gc left', 'Gb left', 'Ga inner', 'Gb right', 'Gc right']


def block_average(arr, window):
    arr = np.asarray(arr)
    if window < 1 or int(window) != window:
        raise ValueError('window must be a positive integer')
    n = len(arr) // window * window
    return arr[:n].reshape(-1, window).mean(axis=1)


def fit_segments(x, y, n_segments=5, min_points=50, min_width=0.3):
    """Global minimum of segmented OLS SSE over observed voltage codes.

    Every sample belongs to exactly one segment. Breakpoints are midpoints
    between adjacent codes. Dynamic programming avoids stochastic searches.
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
    prefixes = [np.r_[0., np.cumsum(v)] for v in (x, y, x*x, x*y, y*y)]
    costs = np.full((m + 1, m + 1), np.inf)
    for a in range(m):
        b = np.arange(a + 2, m + 1)
        n = cuts[b] - cuts[a]
        sx, sy, sxx, sxy, syy = [v[cuts[b]] - v[cuts[a]] for v in prefixes]
        vx = sxx - sx*sx/n
        valid = (n >= min_points) & (codes[b-1] - codes[a] >= min_width) & (vx > 0)
        sse = syy - sy*sy/n - (sxy - sx*sy/n)**2 / np.maximum(vx, 1e-300)
        costs[a, b[valid]] = np.maximum(sse[valid], 0)
    dp = np.full((n_segments+1, m+1), np.inf)
    prev = np.full(dp.shape, -1, int)
    dp[0, 0] = 0
    for k in range(1, n_segments+1):
        for b in range(1, m+1):
            vals = dp[k-1, :b] + costs[:b, b]
            a = int(np.argmin(vals))
            dp[k, b], prev[k, b] = vals[a], a
    if not np.isfinite(dp[-1, -1]):
        raise ValueError('no partition satisfies the minimum width/sample count')
    ends, b = [m], m
    for k in range(n_segments, 0, -1):
        b = prev[k, b]
        ends.append(b)
    ends = list(reversed(ends))
    bounds = [float(x[0])] + [float((codes[e-1]+codes[e])/2) for e in ends[1:-1]] + [float(x[-1])]
    fits = []
    for k, (a, b) in enumerate(zip(ends[:-1], ends[1:])):
        xx, yy = x[cuts[a]:cuts[b]], y[cuts[a]:cuts[b]]
        (slope, intercept), cov = np.polyfit(xx, yy, 1, cov=True)
        residual = yy - (slope*xx + intercept)
        sst = np.sum((yy - yy.mean())**2)
        fits.append(dict(label=LABELS[k] if n_segments == 5 else str(k+1),
                         v_lo=bounds[k], v_hi=bounds[k+1], n=len(xx),
                         slope_S=float(slope), intercept_A=float(intercept),
                         slope_se_S=float(np.sqrt(cov[0,0])),
                         intercept_se_A=float(np.sqrt(cov[1,1])),
                         r2=float(1-np.sum(residual**2)/sst) if sst else None))
    sst = float(np.sum((y-y.mean())**2))
    return fits, float(1-dp[-1,-1]/sst) if sst else None


def slope_ranges(fits, r0=R0):
    ga = fits[2]['slope_S']
    out = []
    for side, i in [('left', 1), ('right', 3)]:
        gb = fits[i]['slope_S']
        valid = ga < gb < 0
        low, high = (-1/ga, -1/gb) if valid else (None, None)
        out.append(dict(side=side, valid=valid, total_lo_ohm=low, total_hi_ohm=high,
                        pot_lo_ohm=low-r0 if valid else None,
                        pot_hi_ohm=high-r0 if valid else None))
    return out


def analyze(path, source='CH1(V)', element='CH2(V)', rs=RS, source_range=(-9.33, 8.28)):
    if not np.isfinite(rs) or rs <= 0:
        raise ValueError('shunt resistance must be positive')
    data = pd.read_csv(path)
    t, vv, vi = (data[c].to_numpy(float) for c in ['Time(s)', source, element])
    if source == element or not np.isfinite([t,vv,vi]).all() or np.any(np.diff(t) <= 0):
        raise ValueError('invalid channels or timestamps')
    current = (vv-vi)/rs
    # Preserve the existing exclusion of drive-end saturation artifacts.
    # These are source-voltage limits, not manually chosen segment boundaries.
    keep = (vv >= source_range[0]) & (vv <= source_range[1])
    fits, r2 = fit_segments(vi[keep], current[keep])
    window = min(501, (len(vv)-1)//2*2+1)
    smoothed = savgol_filter(vv, window, 3)
    derivative = np.gradient(smoothed, t)
    threshold = 0.1 * np.median(np.abs(derivative))
    rising, falling = derivative > threshold, derivative < -threshold
    overlap = []
    edges = np.linspace(vi[keep].min(),vi[keep].max(),81)
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = keep & (vi >= lo) & (vi < hi)
        up, down = current[mask & rising], current[mask & falling]
        if len(up) >= 10 and len(down) >= 10:
            overlap.append(((lo+hi)/2,up.mean(),down.mean()))
    overlap = np.asarray(overlap).reshape(-1,3)
    peaks, _ = find_peaks(smoothed,prominence=0.5*np.ptp(smoothed))
    troughs, _ = find_peaks(-smoothed, prominence=0.5*np.ptp(smoothed))
    turns = np.sort(np.r_[peaks,troughs])
    freq = float(1/(2*np.median(np.diff(t[turns])))) if len(turns)>1 else None
    result = dict(schema_version=1, input=str(Path(path).resolve()),
                  voltage_column=element, source_column=source, voltage_basis='element',
                  current_definition=f'({source} - {element}) / {rs:g} ohm',
                  shunt_ohm=rs, r0_ohm=R0, source_fit_range_V=list(source_range),
                  n_input=len(t), n_fit=int(keep.sum()), combined_r2=r2,
                  segments=fits, ranges=slope_ranges(fits), drive_frequency_hz=freq,
                  directional_rms_difference_A=float(np.sqrt(np.mean((overlap[:,1]-overlap[:,2])**2))) if len(overlap) else None,
                  uncertainty_note='OLS standard errors only; quantisation, channel gains, correlated errors and drift are not included.')
    return result, (t,vv,vi,current,keep,rising,falling,overlap)


def plot_analysis(result, arrays, out):
    t,vv,vi,current,keep,rising,falling,overlap = arrays
    fig, (ax, diff) = plt.subplots(2,1,figsize=(10,7.8),height_ratios=[3,1],layout='constrained')
    for mask,label,color in [(rising,'Rising source sweep','#0072B2'),(falling,'Falling source sweep','#D55E00')]:
        table = pd.DataFrame({'v':vi[mask & keep],'i':current[mask & keep]}).groupby('v').i.mean()
        ax.plot(table.index,table.values*1e3,'.-',ms=3,lw=0.8,color=color,label=label)
    for k,f in enumerate(result['segments']):
        xx=np.array([f['v_lo'],f['v_hi']])
        ax.plot(xx,(f['slope_S']*xx+f['intercept_A'])*1e3,color='black',lw=1.6,
                label='Five fitted segments' if k==0 else None)
    ax.set(xlabel=f"Voltage across nonlinear element, {result['voltage_column'].split('(')[0]} (V)",ylabel='Current into nonlinear element (mA)',
           title='M1 · Diode characteristic and sweep-direction check')
    ax.legend(fontsize=9); ax.grid(alpha=.2)
    if len(overlap):
        diff.plot(overlap[:,0],(overlap[:,1]-overlap[:,2])*1e3,color='#0072B2',lw=1)
    diff.axhline(0,color='0.4',lw=.8)
    diff.set(xlabel='Voltage across nonlinear element (V)',ylabel='Rising − falling\n(mA)'); diff.grid(alpha=.2)
    fig.savefig(out/'current_vs_voltage.png',dpi=200); plt.close(fig)
    fig, (ax,resax)=plt.subplots(2,1,figsize=(10,7.8),height_ratios=[3,1],layout='constrained',sharex=True)
    ax.scatter(vi[keep],current[keep]*1e3,s=2,alpha=.12,color='#0072B2',rasterized=True,label='Measured samples')
    for k,f in enumerate(result['segments']):
        xx=np.array([f['v_lo'],f['v_hi']]); color=plt.get_cmap('tab10')(k)
        ax.plot(xx,(f['slope_S']*xx+f['intercept_A'])*1e3,color=color,lw=2,
                label=f"{f['label']}: {f['slope_S']*1e3:.3f} mS")
        mask=keep & (vi>=xx[0]) & ((vi<xx[1]) if k<4 else (vi<=xx[1]))
        resax.scatter(vi[mask],(current[mask]-f['slope_S']*vi[mask]-f['intercept_A'])*1e3,s=2,alpha=.12,color=color,rasterized=True)
        if k<4:
            ax.axvline(xx[1],ls=':',color='0.4',lw=.8)
            ax.text(xx[1],.98,f'{xx[1]:.2f} V',rotation=90,va='top',ha='right',transform=ax.get_xaxis_transform(),fontsize=8)
    ax.set(ylabel='Current (mA)',title=f"M1 · Optimized element-voltage breakpoints (R² = {result['combined_r2']:.4f})")
    ax.legend(fontsize=8,loc='lower right'); ax.grid(alpha=.2)
    resax.axhline(0,color='0.4',lw=.8)
    resax.set(xlabel=f"Voltage across nonlinear element, {result['voltage_column'].split('(')[0]} (V)",ylabel='Fit residual\n(mA)'); resax.grid(alpha=.2)
    fig.savefig(out/'optimized_breakpoints.png',dpi=200); plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv',nargs='?',type=Path,default=HERE/'trace1.csv')
    parser.add_argument('--source',default='CH1(V)'); parser.add_argument('--element',default='CH2(V)')
    parser.add_argument('--rs',type=float,default=RS)
    parser.add_argument('--source-range',type=float,nargs=2,default=(-9.33,8.28))
    parser.add_argument('--out-dir',type=Path,default=HERE)
    args=parser.parse_args(); args.out_dir.mkdir(parents=True,exist_ok=True)
    result,arrays=analyze(args.csv,args.source,args.element,args.rs,args.source_range)
    (args.out_dir/'diode_fit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    with (args.out_dir/'diode_segments.csv').open('w',newline='') as fh:
        w=csv.DictWriter(fh,list(result['segments'][0])); w.writeheader(); w.writerows(result['segments'])
    plot_analysis(result,arrays,args.out_dir)
    print(f"M1: voltage={args.element}, source={args.source}; fit R2={result['combined_r2']:.6f}")
    for f in result['segments']:
        print(f"{f['label']:10}: {f['v_lo']:.3f} .. {f['v_hi']:.3f} V, G={1e3*f['slope_S']:.6f} mS, intercept={1e3*f['intercept_A']:.5f} mA")
    for r in result['ranges']:
        print(f"M2 {r['side']}: {r}")
    print('M2 is a slope-only necessary estimate. Check intersections with the measured offset retained.')
    print(f"Sweep frequency: {result['drive_frequency_hz']} Hz; rising/falling RMS difference: {result['directional_rms_difference_A']} A")


if __name__ == '__main__':
    main()
