"""
Figures for sections 2.2-2.5 of the report, in the report's style (report_style.py).

Reads only existing outputs and records; changes nothing in the pipeline:
    <sweep>_bifurcation_points.csv   maxima of V1 per record (bifurcation.py)
    <sweep>_rpot.csv                 rpot.py's resistance per record
    onset_records.csv                resistances at the top of the dial (onset.py)
    ../uncertainty/<sweep>_uncertainty.csv   per-record statistical u(Rpot)

Writes to ../report_figures/:
    bifurcation.pdf/.png   maxima of V1 against Rpot, both sweeps, with the
                           equilibrium branch at the top of the dial (figure*)
    gallery.pdf/.png       V2 against V1 for the five regimes (figure*)
    cascade.pdf/.png       maxima of V1 through the period-doubling cascade, forward
                           sweep, with R_1..R_3 from delta_bench.json (column)
    return_maps.pdf/.png   M_(n+1) against M_n, one record per regime along the sweep (figure*)
    lyapunov.pdf/.png      Lyapunov exponent of every forward record (direct) and the three
                           return-map estimates (column)
    simulation.pdf/.png    maxima of V1 of the plan's model and of the identified circuit
                           (simulate.py --tag _nominal / --model bench --datasheet-lag
                           --tag _bench_datasheet) over the measured points, both sweep
                           directions (figure*)

Only the records rpot.py accepts are used (the flagged ones at the top of the dial
are excluded). For accepted records whose V1 half-range is below
onset.CAL_AMP_V, rpot.py's value is biased low, so they are placed at onset.py's
resistance.

Usage:
    python report_figures.py
"""
import csv
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import report_style  # noqa: E402
from lorenz_map import maxima  # noqa: E402
from onset import CAL_AMP_V, OSC_MIN_V  # noqa: E402
from scope_data import read_scope  # noqa: E402

OUT = os.path.join(ROOT, 'report_figures')
SMOOTH_S = 16e-6      # Savitzky-Golay window, T/20 of the ~330 us period (as lorenz_map.maxima)
GALLERY = [('period 1', 'trace54.csv'), ('period 2', 'trace82.csv'), ('single scroll', 'trace164.csv'),
           ('double scroll', 'trace270.csv'), ('large cycle', 'trace411.csv')]


def onset_placed(sweep):
    """Rpot by filename of the accepted small cycles that onset.py places by their mean voltages (rpot.py reads them low)."""
    top = pd.read_csv(os.path.join(HERE, 'onset_records.csv'))
    small = top[(top.sweep == sweep) & (top.sidecar_status == 'ok')
                & (top.v1_half_range_V < CAL_AMP_V) & (top.A1_mV >= 1e3 * OSC_MIN_V)]
    return small.set_index('filename').rpot_means_ohm


def forward_rpot():
    """({filename: Rpot}, u(Rpot) by filename) of the forward records, for the panel titles."""
    with open(os.path.join(HERE, 'forward_rpot.csv'), newline='') as fh:
        side = {r['filename']: float(r['rpot_ohm']) for r in csv.DictReader(fh) if r['rpot_ohm']}
    unc = pd.read_csv(os.path.join(ROOT, 'uncertainty', 'forward_uncertainty.csv')).set_index('filename').u_r_block_ohm
    return side, unc


def bifurcation_data(sweep):
    """(R, maxima) arrays for one sweep's accepted records."""
    pts = pd.read_csv(os.path.join(HERE, f'{sweep}_bifurcation_points.csv'))
    placed = onset_placed(sweep)
    # rpot.py's maxima, except the small cycles, which move to onset.py's resistance
    pts = pts[~pts.filename.isin(placed.index)]
    R, M = [pts.rpot_ohm.to_numpy()], [pts.max_v.to_numpy()]
    for name, r in placed.items():
        m, _ = maxima(os.path.join(HERE, sweep, name), 'CH1', 0.02, None)
        R.append(np.full(len(m), r))
        M.append(m)
    return np.concatenate(R), np.concatenate(M)


def break_axis(top, ax, ratio, d=0.006):
    """Join a short upper panel and the main panel below it (heights 1 : ratio) as one broken vertical axis."""
    top.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    top.tick_params(bottom=False, labelbottom=False, which='both')
    ax.tick_params(top=False, which='both')
    for a, y, h in ((top, 0, 4 * d), (ax, 1, 4 * d / ratio)):     # same slant on both panels
        kw = dict(transform=a.transAxes, color='k', lw=0.6, clip_on=False)
        a.plot((-d, d), (y - h, y + h), **kw)
        a.plot((1 - d, 1 + d), (y - h, y + h), **kw)


def plot_bifurcation():
    """Maxima of V1 against Rpot, both sweeps, with the large cycle on a broken axis (figure*)."""
    import matplotlib.pyplot as plt
    c = report_style.COLORS
    fig, (top, ax) = plt.subplots(2, 1, sharex=True, figsize=report_style.figsize(True, aspect=0.5),
                                  height_ratios=[1, 5], layout='constrained')
    fig.get_layout_engine().set(hspace=0.01, h_pad=0.01)
    style = {'forward': dict(color=c[0], label='Forward (decreasing $R_\\mathrm{pot}$)', z=2),
             'back': dict(color=c[1], label='Back (increasing $R_\\mathrm{pot}$)', z=3)}
    for sweep, st in style.items():
        R, M = bifurcation_data(sweep)
        for a in (top, ax):
            a.scatter(R, M, s=0.15, lw=0, color=st['color'], rasterized=True, zorder=st['z'])
        ax.plot([], [], ls='none', marker='o', ms=3, color=st['color'], label=st['label'])
    top.set_ylim(6.0, 6.7)
    ax.set_ylim(-3.6, 2.9)
    ax.set_xlim(300, 915)
    break_axis(top, ax, ratio=5)
    ax.set_xlabel(r'$R_\mathrm{pot}$ ($\Omega$)')
    fig.supylabel(r'Local maxima of $V_1$ (V)', fontsize=9)
    labels = [(315, 6.52, 'large cycle', top), (500, 2.35, 'double scroll', ax), (705, 0.6, 'single scroll', ax),
              (840, -0.2, 'period 1', ax)]
    for x, y, s, a in labels:
        a.text(x, y, s, fontsize=8, ha='left')
    ax.legend(loc='lower left', markerscale=1.2)
    os.makedirs(OUT, exist_ok=True)
    report_style.save(fig, os.path.join(OUT, 'bifurcation.png'))
    plt.close(fig)


# The start of chaos lies between the last record with a steady period 8 and the first chaotic
# one. period8_check.py shows that of the five records cascade_periods.py labels period 8, only
# trace92, 93 and 95 are steady; trace96 (between these two) reverses phase and is neither.
CHAOS_BRACKET = ('trace93.csv', 'trace97.csv')


def cascade_marks(sweep='forward'):
    """{label: (R, half-bracket)}: R_1..R_3 from delta_bench.py, and the start of chaos from CHAOS_BRACKET."""
    with open(os.path.join(HERE, 'delta_bench.json')) as fh:
        res = json.load(fh)
    marks = {f'$R_{n}$': (res[f'R{n}']['value'], res[f'R{n}']['u']) for n in (1, 2, 3)}
    with open(os.path.join(HERE, f'cascade_periods_{sweep}.txt')) as fh:
        txt = fh.read()
    R = {m[1]: float(m[2]) for m in re.finditer(r'(trace\d+\.csv)\s+([\d.]+)\s+\d+', txt)}
    hi, lo = (R[name] for name in CHAOS_BRACKET)
    marks['chaos'] = (0.5 * (hi + lo), 0.5 * (hi - lo))
    return marks


def plot_cascade():
    """Maxima of V1 through the period-doubling cascade, forward sweep: (a) overview, (b) zoom."""
    import matplotlib.pyplot as plt
    from matplotlib.transforms import blended_transform_factory
    R, M = bifurcation_data('forward')
    marks = cascade_marks()
    fig, axes = plt.subplots(2, 1, figsize=report_style.figsize(aspect=1.1), layout='constrained')
    views = [((741, 791), ['$R_1$', '$R_2$'], [(783, 'period 1'), (765, 'period 2')]),
             ((745.2, 757.8), ['$R_2$', '$R_3$', 'chaos'], [(753.0, 'period 4'), (749.3, '8')])]
    for ax, (lim, labelled, regions), tag in zip(axes, views, '(a) (b)'.split()):
        sel = (R > lim[0]) & (R < lim[1])
        ax.scatter(R[sel], M[sel], s=0.4, lw=0, color=report_style.COLORS[0], rasterized=True, zorder=2)
        tr = blended_transform_factory(ax.transData, ax.transAxes)
        for label, (r, h) in marks.items():
            if not lim[0] < r < lim[1]:
                continue
            ax.axvspan(r - h, r + h, color='0.88', lw=0, zorder=0)
            ax.axvline(r, color='0.35', lw=0.5, ls='--', zorder=1)
            if label in labelled:
                ax.text(r, 1.02, label, transform=tr, ha='center', va='bottom', fontsize=8)
        for x, s in regions:
            ax.text(x, 0.95, s, transform=tr, ha='center', va='top', fontsize=8)
        ax.set_xlim(*lim)
        ax.set_ylim(-1.0, 0.0)
        ax.set_ylabel(r'Local maxima of $V_1$ (V)')
        ax.text(0.97, 0.04, tag, transform=ax.transAxes, ha='right')
    axes[1].set_xlabel(r'$R_\mathrm{pot}$ ($\Omega$)')
    report_style.save(fig, os.path.join(OUT, 'cascade.png'))
    plt.close(fig)


def value_unc(x, u):
    """'804.4(1)': value to the digit of the uncertainty's first significant figure."""
    digits = max(0, -int(np.floor(np.log10(u))))
    return f'{x:.{digits}f}({round(u * 10 ** digits):d})'


def plot_gallery(window_s=0.04):
    """V2 against V1 of the five regimes, 40 ms each (figure*)."""
    import matplotlib.pyplot as plt
    side, unc = forward_rpot()
    fig, axes = plt.subplots(1, 5, figsize=report_style.figsize(True, aspect=0.25), layout='constrained')
    fig.get_layout_engine().set(wspace=0.02, w_pad=0.02)
    for k, (ax, (label, name)) in enumerate(zip(axes, GALLERY)):
        t, d = read_scope(os.path.join(HERE, 'forward', name), ('CH1', 'CH2'), min_samples=10)
        dt = float(np.median(np.diff(t)))
        n = int(window_s / dt)
        # V2 spans only ~40 ADC levels; a Savitzky-Golay window of T/20 removes the steps
        # without changing the orbit's shape
        w = max(5, int(round(SMOOTH_S / dt)) | 1)
        v1, v2 = (savgol_filter(d[:n, j], w, 3) for j in (0, 1))
        ax.plot(v1, v2, lw=0.3, color=report_style.COLORS[0], rasterized=True)
        ax.set_title(f'({"abcde"[k]}) {label}\n{value_unc(side[name], unc[name])} $\\Omega$', fontsize=8, pad=3)
        ax.set_xlabel(r'$V_1$ (V)', labelpad=1)
        ax.tick_params(labelsize=7)
        ax.locator_params(nbins=4)
    axes[0].set_ylabel(r'$V_2$ (V)')
    report_style.save(fig, os.path.join(OUT, 'gallery.png'))
    plt.close(fig)


# one record per regime along the forward sweep (map quality from forward_lyapunov.csv); (f)-(h) are
# the three long double-scroll records whose maps give the Lyapunov exponents (plan M3-M5)
RETURN_MAPS = [('period 4', 'trace90.csv'), ('two-band chaos', 'trace100.csv'), ('single scroll', 'trace118.csv'),
               ('period-3 window', 'trace141.csv'), ('single scroll', 'trace164.csv'),
               ('double scroll', 'trace229.csv'), ('double scroll', 'trace270.csv'), ('double scroll', 'trace342.csv')]


def plot_return_maps():
    """M_(n+1) against M_n, one record per regime along the forward sweep, with the diagonal (figure*)."""
    import matplotlib.pyplot as plt
    side, unc = forward_rpot()
    fig, axes = plt.subplots(2, 4, figsize=report_style.figsize(True, aspect=0.56), layout='constrained')
    fig.get_layout_engine().set(wspace=0.03, hspace=0.04)
    for k, (ax, (label, name)) in enumerate(zip(axes.flat, RETURN_MAPS)):
        m, _ = maxima(os.path.join(HERE, 'forward', name))
        pad = 0.08 * np.ptp(m)
        lo, hi = m.min() - pad, m.max() + pad
        ax.plot([lo, hi], [lo, hi], color='0.5', lw=0.5, ls='--', zorder=0)
        ax.scatter(m[:-1], m[1:], s=0.8, lw=0, color=report_style.COLORS[0], rasterized=True)
        ax.set(xlim=(lo, hi), ylim=(lo, hi), aspect='equal')
        ax.set_title(f'({"abcdefgh"[k]}) {label}\n{value_unc(side[name], unc[name])} $\\Omega$', fontsize=8, pad=3)
        ax.tick_params(labelsize=7)
        ax.locator_params(nbins=4)
        if k >= 4:
            ax.set_xlabel(r'$M_n$ (V)', labelpad=1)
        if k % 4 == 0:
            ax.set_ylabel(r'$M_{n+1}$ (V)')
    report_style.save(fig, os.path.join(OUT, 'return_maps.png'))
    plt.close(fig)


LYAPUNOV_MAP_RECORDS = ['trace229.csv', 'trace270.csv', 'trace342.csv']   # the three long records of plan M3-M5


def plot_lyapunov():
    """Largest Lyapunov exponent of every forward record, measured directly (Rosenstein), with the
    return-map estimates of the three long double-scroll records (column)."""
    import matplotlib.pyplot as plt
    d = pd.read_csv(os.path.join(HERE, 'forward_lyapunov.csv'))
    d = d[d.rosenstein_status == 'ok'].set_index('filename')
    # the small cycles near the onset sit at their mean-voltage resistance, as in the bifurcation diagram
    R = d.rpot_ohm.copy()
    for name, r in onset_placed('forward').items():
        if name in R.index:
            R[name] = r
    c = report_style.COLORS
    fig, ax = plt.subplots(figsize=report_style.figsize(aspect=0.72), layout='constrained')
    ax.axhline(0, color='0.5', lw=0.5)
    ax.errorbar(R, d.lambda_rosenstein_per_s / 1e3, yerr=d.lambda_rosenstein_err_per_s / 1e3, ls='none',
                marker='o', ms=1.5, elinewidth=0.4, capsize=0, color=c[0], label='Direct (Rosenstein)')
    m = d.loc[LYAPUNOV_MAP_RECORDS]
    ax.errorbar(m.rpot_ohm, m.lambda_map_per_s / 1e3, yerr=m.lambda_map_err_per_s / 1e3, ls='none',
                marker='s', ms=3.5, mfc='white', mew=0.8, elinewidth=0.8, capsize=1.5, color=c[1],
                label='Return map', zorder=3)
    ax.set(xlabel=r'$R_\mathrm{pot}$ ($\Omega$)', ylabel=r'$\lambda$ ($10^3$ s$^{-1}$)', xlim=(0, 915), ylim=(-0.5, 3.4))
    ax.legend(loc='upper left', handletextpad=0.3)
    report_style.save(fig, os.path.join(OUT, 'lyapunov.png'))
    plt.close(fig)


# (b) is the identified circuit with node 1's capacitor law and the op-amps at their datasheet
# speed (--datasheet-lag), the model the report describes; the measured stage-B lag stays in RESULTS.md
SIM_MODELS = [('_nominal', '(a) Ideal model: nominal components, element of Table 2'),
              ('_bench_datasheet', '(b) Identified model: components measured from the records')]


def simulated_maxima(tag):
    """{'down': (R, M), 'up': (R, M)} of simulate.py's continuation sweeps, without the resistances where the
    model rests at its equilibrium (the decaying transient there has maxima but no oscillation)."""
    z = np.load(os.path.join(HERE, f'simulated_sweep{tag}.npz'))
    with open(os.path.join(HERE, f'simulated_transitions{tag}.json')) as fh:
        runs = json.load(fh)
    out = {}
    for d in ('down', 'up'):
        rest = [(min(a, b), max(a, b)) for lab, a, b in runs[f'regimes_{d}'] if lab == 'rest']
        R, M = z[f'R_{d}'], z[f'M_{d}']
        keep = np.ones(len(R), bool)
        for lo, hi in rest:
            keep &= ~((R >= lo) & (R <= hi))
        out[d] = (R[keep], M[keep])
    return out


def plot_simulation():
    """Maxima of V1 of the two models over the measured points, both sweep directions (figure*)."""
    import matplotlib.pyplot as plt
    c = report_style.COLORS
    meas = [bifurcation_data(s) for s in ('forward', 'back')]
    Rm, Mm = np.concatenate([m[0] for m in meas]), np.concatenate([m[1] for m in meas])
    fig = plt.figure(figsize=report_style.figsize(True, aspect=0.72), layout='constrained')
    subs = fig.subfigures(2, 1, hspace=0.04)
    for k, (sub, (tag, title)) in enumerate(zip(subs, SIM_MODELS)):
        top, ax = sub.subplots(2, 1, sharex=True, height_ratios=[1, 4])
        sim = simulated_maxima(tag)
        for a in (top, ax):
            a.scatter(Rm, Mm, s=0.15, lw=0, color='0.72', rasterized=True, zorder=1)
            a.scatter(*sim['down'], s=0.15, lw=0, color=c[0], rasterized=True, zorder=2)
            a.scatter(*sim['up'], s=0.15, lw=0, color=c[1], rasterized=True, zorder=3)
        top.set_ylim(5.0, 7.3)
        ax.set_ylim(-4.3, 4.6)
        ax.set_xlim(300, 1000)
        break_axis(top, ax, ratio=4)
        top.set_title(title, fontsize=9, loc='left', pad=2)
        if k == 1:
            ax.set_xlabel(r'$R_\mathrm{pot}$ ($\Omega$)')
        else:
            ax.tick_params(labelbottom=False)
            for col, lab in (('0.72', 'Measured, both sweeps'), (c[0], r'Model, decreasing $R_\mathrm{pot}$'),
                             (c[1], r'Model, increasing $R_\mathrm{pot}$')):
                ax.plot([], [], ls='none', marker='o', ms=3, color=col, label=lab)
            ax.legend(loc='lower left', markerscale=1.2)
        sub.supylabel(r'Local maxima of $V_1$ (V)', fontsize=9)
    report_style.save(fig, os.path.join(OUT, 'simulation.png'))
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use('Agg')
    report_style.apply()
    plot_bifurcation()
    plot_gallery()
    plot_cascade()
    plot_return_maps()
    plot_lyapunov()
    plot_simulation()
    print(f'wrote bifurcation, gallery, cascade, return_maps, lyapunov and simulation to {OUT}')


if __name__ == '__main__':
    main()
