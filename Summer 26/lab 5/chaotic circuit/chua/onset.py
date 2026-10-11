"""
N1 from the bench: where the oscillation starts, on the resistance scale of rpot.py.

rpot.py's divider fit (the *_rpot.csv sidecars) reads small oscillations low:
noise on V1 and V2 attenuates the fitted coefficients (errors in variables), B
more than A because V2 swings less, so R0*B/A is biased down, by up to 13 ohm on
the accepted records nearest the onset. The records rpot.py flags (the top of
the dial, where the circuit rests at its equilibrium or barely oscillates) are
excluded. This script leaves rpot.py and the sidecars untouched and, for the
accepted records at the top of each sweep:

1. Reads R from the mean voltages. The divider holds for them too,
       m3 = a (R0 m1 + R m2) / (R0 + R) + C,
   with a = A + B of the divider fit (the midpoint channel's relative gain) and C
   the combined channel offset. C is calibrated so that this reproduces the
   sidecar R on records of the same V/div block whose V1 half-range is at least
   CAL_AMP_V, where the fit's bias is negligible; a block without three such
   records borrows the nearest block that has them.
2. Measures the amplitude A1 of V1's fundamental by a least-squares sine fit at
   the dominant frequency between 2.5 and 3.6 kHz. Averaging over the 200 ms
   record resolves amplitudes far below one ADC step.
3. Locates the onset by extrapolating A1^2 = k (R_H - R) (supercritical Hopf)
   over the oscillating records with A1 < FIT_MAX_V, with chi-square and p. Each
   point carries the uncertainty of A1^2 and k times the calibration rms in R.
   (With the flagged records included, the direct bracket between the last
   record without a limit cycle and the first with one agrees: 907.3-907.8 ohm,
   R_H = 907.8 +- 0.5 ohm.)

Outputs: onset.txt, onset_records.csv, onset.pdf/.png (report figure).

Usage:
    python onset.py
"""
import csv
import os
import sys

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2

from scope_data import R0, read_scope
from sweeplib import list_csvs

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import report_style  # noqa: E402

CAL_AMP_V = 0.6      # V1 half-range above which rpot.py's divider fit is unbiased to < 1 ohm
MIN_CAL = 3          # calibration records needed in a block
OSC_MIN_V = 10e-3    # limit cycle: above the noise-driven precursor (<= 2.4 mV), below the smallest cycle (29 mV)
FIT_MAX_V = 0.2      # A1^2 is linear in R only near the onset
BAND_HZ = (2500, 3600)


def sine_amplitude(t, x):
    """(amplitude, frequency, standard error of the amplitude) of the dominant line."""
    x = x - x.mean()
    n, dt = len(x), float(np.median(np.diff(t)))
    F = np.abs(np.fft.rfft(x * np.hanning(n)))
    f = np.fft.rfftfreq(n, dt)
    band = (f > BAND_HZ[0]) & (f < BAND_HZ[1])
    f0 = f[np.argmax(np.where(band, F, 0))]
    best = None
    for ff in np.linspace(f0 - 5, f0 + 5, 81):            # refine within +- one FFT bin
        M = np.column_stack([np.cos(2 * np.pi * ff * t), np.sin(2 * np.pi * ff * t)])
        c = np.linalg.lstsq(M, x, rcond=None)[0]
        amp = float(np.hypot(*c))
        if best is None or amp > best[0]:
            best = (amp, float(ff), float(np.std(x - M @ c) * np.sqrt(2 / n)))
    return best


def read_record(folder, name, sidecar):
    """Mean voltages, divider gain, V1 amplitude and the sidecar's verdict of one record."""
    t, d = read_scope(os.path.join(folder, name), ('CH1', 'CH2', 'CH3'), min_samples=10)
    X = np.column_stack([d[:, 0] - d[:, 0].mean(), d[:, 1] - d[:, 1].mean()])
    A, B = np.linalg.lstsq(X, d[:, 2] - d[:, 2].mean(), rcond=None)[0]
    a1, f, u = sine_amplitude(t, d[:, 0])
    row = sidecar.get(name, {})
    r_side = float(row['rpot_ohm']) if row.get('rpot_ohm') else np.nan
    codes = np.unique(d[:, 0])
    return dict(name=name, mean=d.mean(0), a_plus_b=float(A + B), half_range=float(0.5 * np.ptp(d[:, 0])),
                q=float(np.median(np.diff(codes))), A1=a1, f=f, uA1=u, r_side=r_side,
                side_ok=row.get('status') == 'ok')


def from_means(mean, a, c):
    """R from the mean voltages through the divider m3 = a (R0 m1 + R m2)/(R0 + R) + c."""
    m1, m2, m3 = mean
    return float(R0 * (a * m1 - (m3 - c)) / ((m3 - c) - a * m2))


def usable_for_calibration(r):
    """Accepted by rpot.py, large enough for its fit to be unbiased, and away from the origin."""
    return r['side_ok'] and r['half_range'] >= CAL_AMP_V and abs(r['mean'][0]) > 1.0


def calibrate(cal):
    """(a, C, rms against rpot.py) of the mean-voltage divider on the calibration records `cal`."""
    a = float(np.mean([x['a_plus_b'] for x in cal]))

    def cost(c):
        return sum((from_means(x['mean'], a, c) - x['r_side']) ** 2 for x in cal)
    c = float(minimize_scalar(cost, bounds=(-0.3, 0.3), method='bounded').x)
    return a, c, float(np.sqrt(cost(c) / len(cal)))


def top_of_dial(sweep):
    """Records from the top-of-dial end of the sweep, in acquisition order, each placed by its mean voltages."""
    folder = os.path.join(HERE, sweep)
    with open(os.path.join(HERE, f'{sweep}_rpot.csv'), newline='') as fh:
        sidecar = {r['filename']: r for r in csv.DictReader(fh)}
    # only the records rpot.py accepts (status ok); the flagged ones at the top of the dial are excluded
    names = [n for n in list_csvs(folder) if sidecar.get(n, {}).get('status') == 'ok']
    # forward starts at the top of the dial, back ends there
    order = names if sweep == 'forward' else names[::-1]

    # Read from the top of the dial down, one V/div block (one quantisation step) at a time,
    # until a complete block holds enough calibration records.
    blocks, current = [], []
    for name in order:
        r = read_record(folder, name, sidecar)
        if current and abs(r['q'] - current[-1]['q']) > 1e-4:
            blocks.append(current)
            if sum(usable_for_calibration(x) for x in current) >= MIN_CAL:
                break
            current = []
        current.append(r)
    else:
        blocks.append(current)

    # calibrate each block, borrowing the nearest block with enough records
    enough = [j for j in range(len(blocks)) if sum(usable_for_calibration(x) for x in blocks[j]) >= MIN_CAL]
    for i, blk in enumerate(blocks):
        j = min(enough, key=lambda j: abs(j - i))
        cal = [x for x in blocks[j] if usable_for_calibration(x)]
        a, c, rms = calibrate(cal)
        span = f"{cal[0]['name']}..{cal[-1]['name']}" if sweep == 'forward' else f"{cal[-1]['name']}..{cal[0]['name']}"
        for x in blk:
            x.update(R=from_means(x['mean'], a, c), a=a, C=c, cal_rms=rms, cal=span)
    recs = [x for blk in blocks for x in blk]
    return recs if sweep == 'forward' else recs[::-1]


def onset(recs):
    """The first limit cycle, the last record without one, and the A1^2 extrapolation."""
    osc = [r for r in recs if r['A1'] >= OSC_MIN_V]
    quiet = [r for r in recs if r['A1'] < OSC_MIN_V]    # empty when only accepted records are used
    first = max(osc, key=lambda r: r['R'])
    last_quiet = min((r for r in quiet if r['R'] > first['R']), key=lambda r: r['R'], default=None)
    out = dict(first_osc=first, last_quiet=last_quiet)
    pts = [r for r in osc if r['A1'] < FIT_MAX_V]
    if len(pts) >= 3:
        out.update(hopf_fit(pts))
        # the same fit over the n smallest cycles, to show how the result depends on FIT_MAX_V
        ordered = sorted(osc, key=lambda r: r['A1'])
        out['by_n'] = [hopf_fit(ordered[:n]) for n in range(3, min(len(ordered), 7) + 1)]
    return out


def hopf_fit(pts):
    """A1^2 = k (R_H - R) by effective variance: sigma_y^2 + (k sigma_R)^2."""
    R = np.array([r['R'] for r in pts])
    y = np.array([r['A1'] ** 2 for r in pts])
    uy_a = np.array([2 * r['A1'] * r['uA1'] for r in pts])
    sR = pts[0]['cal_rms']
    k, b = np.polyfit(R, y, 1)
    for _ in range(3):
        s = np.sqrt(uy_a ** 2 + (k * sR) ** 2)
        W = 1 / s ** 2
        X = np.column_stack([R, np.ones_like(R)])
        cov = np.linalg.inv(X.T @ (X * W[:, None]))
        k, b = cov @ (X.T @ (W * y))
    c2 = float(np.sum(((y - (k * R + b)) / s) ** 2))
    dof = len(R) - 2
    RH = -b / k
    g = np.array([b / k ** 2, -1 / k])                # dRH/dslope, dRH/db
    # the fitted slope is -k of the Hopf law A1^2 = k (R_H - R)
    return dict(RH=float(RH), uRH=float(np.sqrt(g @ cov @ g)), k=float(-k), uk=float(np.sqrt(cov[0, 0])),
                chi2=c2, dof=dof, p=float(chi2.sf(c2, dof)), n_fit=len(R),
                R_fit=(float(R.min()), float(R.max())), A_max=float(max(r['A1'] for r in pts)))


def csv_rows(sweep, recs):
    """The rows of onset_records.csv for one sweep."""
    return [dict(sweep=sweep, filename=r['name'], rpot_means_ohm=round(r['R'], 2),
                 rpot_sidecar_ohm='' if not np.isfinite(r['r_side']) else round(r['r_side'], 2),
                 sidecar_status='ok' if r['side_ok'] else 'flagged',
                 v1_half_range_V=round(r['half_range'], 4), A1_mV=round(1e3 * r['A1'], 4),
                 u_A1_mV=round(1e3 * r['uA1'], 4), f_Hz=round(r['f'], 1), calibration=r['cal'])
            for r in recs]


def sweep_report(sweep, recs, res):
    """The lines of onset.txt for one sweep."""
    f, q = res['first_osc'], res['last_quiet']
    lines = [f'== {sweep} sweep ({len(recs)} records at the top of the dial)',
             f"   mean-voltage divider calibrated on {recs[0]['cal']}: a = {recs[0]['a']:.4f}, "
             f"C = {1e3 * recs[0]['C']:+.1f} mV, rms against rpot.py {recs[0]['cal_rms']:.2f} ohm"
             + ('' if sweep == 'forward' else
                f"; top block calibrated on {recs[-1]['cal']}, rms {recs[-1]['cal_rms']:.2f} ohm"),
             f"   first limit cycle: {f['name']} at {f['R']:.1f} ohm, A1 = {1e3 * f['A1']:.1f} mV, "
             f"f = {f['f']:.1f} Hz (T = {1e6 / f['f']:.1f} us)",
             (f"   last record without: {q['name']} at {q['R']:.1f} ohm, A1 = {1e3 * q['A1']:.2f} mV"
              if q else '   no record without a limit cycle above it')]
    if 'RH' in res:
        lines.append(f"   A1^2 = k (R_H - R) over {res['n_fit']} records, {res['R_fit'][0]:.1f}..{res['R_fit'][1]:.1f} ohm: "
                     f"R_H = {res['RH']:.1f} +- {res['uRH']:.1f} ohm, "
                     f"k = {1e6 * res['k']:.0f} +- {1e6 * res['uk']:.0f} mV^2/ohm, "
                     f"chi2 = {res['chi2']:.1f}, dof = {res['dof']}, chi2/dof = {res['chi2'] / res['dof']:.2f}, "
                     f"p = {res['p']:.3f}")
        lines.append('   the same fit over the n smallest cycles:')
        lines += [f"     n = {v['n_fit']} (A1 up to {1e3 * v['A_max']:.0f} mV, down to {v['R_fit'][0]:.1f} ohm): "
                  f"R_H = {v['RH']:.1f} +- {v['uRH']:.1f} ohm, chi2 = {v['chi2']:.2f}, dof = {v['dof']}, p = {v['p']:.3f}"
                  for v in res['by_n']]
    else:
        lines.append(f"   {sum(r['A1'] < FIT_MAX_V for r in recs)} record(s) below FIT_MAX_V: too few for the A1^2 extrapolation")
    lines.append('')
    return lines


def main():
    lines, rows, results = [], [], {}
    for sweep in ('forward', 'back'):
        recs = top_of_dial(sweep)
        res = onset(recs)
        results[sweep] = (recs, res)
        rows += csv_rows(sweep, recs)
        lines += sweep_report(sweep, recs, res)
    lines += ['Uncertainties are statistical on the rpot.py scale. The calibration choice (which records,',
              'fixed or fitted gain) moves the onset by about +-5 ohm; rpot.py\'s own scale carries the',
              'relative channel gain, known to a few percent.']
    text = '\n'.join(lines)
    print(text)
    with open(os.path.join(HERE, 'onset.txt'), 'w') as fh:
        fh.write(text + '\n')
    with open(os.path.join(HERE, 'onset_records.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    plot(results)


def plot(results):
    """onset.pdf/.png: (a) the amplitude of V1 near the onset with the Hopf law, (b) A1^2 and its fit."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    report_style.apply()
    c = report_style.COLORS
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=report_style.figsize(aspect=1.15), layout='constrained')
    RH = [res['RH'] for _, res in results.values() if 'RH' in res]
    for (sweep, (recs, res)), col, mk, mfc in zip(results.items(), c, 'os', (None, 'none')):
        R = np.array([r['R'] for r in recs])
        A = np.array([r['A1'] for r in recs]) * 1e3
        uR = np.array([r['cal_rms'] for r in recs])      # the per-point R uncertainty the fit uses
        eb = dict(ls='none', marker=mk, color=col, mfc=mfc or col, elinewidth=0.5, capsize=0)
        ax.errorbar(R, A, xerr=uR, ms=2.5, label=f'{sweep.capitalize()} sweep', **eb)
        if 'RH' in res:
            m = (A < 1e3 * FIT_MAX_V) & (A >= 1e3 * OSC_MIN_V)
            ax2.errorbar(R[m], A[m] ** 2 / 1e3, xerr=uR[m], ms=3, **eb)
            rr = np.linspace(R[m].min() - 1, res['RH'], 50)
            ax2.plot(rr, res['k'] * (res['RH'] - rr) * 1e3, color='black', lw=0.8, zorder=0)   # V^2 -> 10^3 mV^2
            ax2.axvline(res['RH'], color='0.5', lw=0.5, ls=':')
            # the Hopf law beyond its fit range, to show where the growth departs from it
            rr = np.linspace(860, res['RH'], 200)
            ax.plot(rr, 1e3 * np.sqrt(res['k'] * (res['RH'] - rr)), color='black', lw=0.7, zorder=0,
                    label=r'$A_1^2 = k(R_\mathrm{H} - R_\mathrm{pot})$')
            ax.axvline(res['RH'], color='0.5', lw=0.5, ls=':')
    xlim = (860, max(RH) + 4 if RH else 910)
    shown = [1e3 * r['A1'] for recs, _ in results.values() for r in recs if xlim[0] <= r['R'] <= xlim[1]]
    ax.set(xlabel=r'$R_\mathrm{pot}$ ($\Omega$)', ylabel=r'Amplitude of $V_1$ (mV)',
           xlim=xlim, ylim=(0, 1.15 * max(shown)))
    handles, labels = ax.get_legend_handles_labels()
    order = sorted(range(len(labels)), key=lambda i: ('sweep' not in labels[i], i))   # sweeps first, then the law
    ax.legend([handles[i] for i in order], [labels[i] for i in order], loc='upper right', bbox_to_anchor=(0.9, 1.0))
    ax2.set(xlabel=r'$R_\mathrm{pot}$ ($\Omega$)', ylabel=r'$A_1^2$ ($10^3$ mV$^2$)')
    ax.text(0.02, 0.05, '(a)', transform=ax.transAxes)
    ax2.text(0.02, 0.05, '(b)', transform=ax2.transAxes)
    report_style.save(fig, os.path.join(HERE, 'onset.png'))
    plt.close(fig)


if __name__ == '__main__':
    main()
