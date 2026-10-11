"""
Are the records that cascade_periods.py labels period 8 really period 8?

cascade_periods.py calls a record period P when its lag-P distance of the V1
maxima is within a factor of two of the record's smallest lag distance. That
says the maxima repeat every P cycles *on average*; it does not say the
pattern is steady through the record. This script tests both, for the
forward records around R_3:

1. Period 8 against period 4. Group the maxima in eights; if the orbit repeats
   every 4, maxima j and j+4 of a group agree to within noise. chi2 over the
   four differences (mean over groups, standard error from their scatter),
   4 dof, with p. The same test at 16 against 8 (8 dof) checks that nothing
   repeats only every 16.
2. Steadiness. The four differences in the first and second half of the
   record, and the lag-1 autocorrelation of each level's deviations. A steady
   orbit keeps the same differences; noise gives an autocorrelation near 0.

The chi2 tests assume a steady orbit, so they are read together with (2).
Figure: the four differences m_j - m_(j+4) group by group through the record,
for a period-4 record, a steady period-8 record, record 94 (the pattern
collapses and re-forms with the opposite phase) and record 96 (repeated phase
reversals, fluctuations far above noise).

Outputs: period8_check.txt, period8_check.pdf/.png (diagnostic, not a report figure).

Usage:
    python period8_check.py
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2

from lorenz_map import maxima
from scope_data import read_scope

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import report_style  # noqa: E402

RECORDS = range(84, 101)
FIGURE = [(90, 'period 4'), (92, 'period 8'), (94, 'period 8, phase reverses'), (96, 'not steady')]


def fold(m, P):
    """The maxima in rows of P, the incomplete last row dropped."""
    k = (len(m) // P) * P
    return m[:k].reshape(-1, P)


def split_test(m, P):
    """Period P against period P/2: chi2 of the P/2 differences between levels j and j+P/2."""
    c = fold(m, P)
    h = P // 2
    d = c[:, :h] - c[:, h:]
    mean, se = d.mean(0), d.std(0, ddof=1) / np.sqrt(len(d))
    c2 = float(np.sum((mean / se) ** 2))
    return mean, c2, h, float(chi2.sf(c2, h))


def steadiness(m, P=8):
    """(newest splittings in the first half, in the second half, mean lag-1 autocorrelation of the levels)."""
    c = fold(m, P)
    h, half = P // 2, len(c) // 2
    d = c[:, :h] - c[:, h:]
    dev = c - c.mean(0)
    ac = float(np.mean([np.corrcoef(dev[:-1, j], dev[1:, j])[0, 1] for j in range(P)]))
    return d[:half].mean(0), d[half:].mean(0), ac


def main():
    r = pd.read_csv(os.path.join(HERE, 'forward_rpot.csv'))
    R = dict(zip(r.filename, r.rpot_ohm))
    lines = [f'{"record":>13} {"R (ohm)":>8} | {"8 vs 4: chi2/dof":>16} {"p":>8} | {"16 vs 8: chi2/dof":>17} {"p":>6} | '
             f'differences m_j - m_(j+4), first half / second half (mV)            | autocorr']
    for n in RECORDS:
        name = f'trace{n}.csv'
        M, _ = maxima(os.path.join(HERE, 'forward', name))
        m = np.asarray(M, float)[len(M) // 4:] * 1e3        # as cascade_periods.py: settle for a quarter
        _, c8, h8, p8 = split_test(m, 8)
        _, c16, h16, p16 = split_test(m, 16)
        s1, s2, ac = steadiness(m)
        lines.append(f'{name:>13} {R[name]:8.2f} | {c8 / h8:16.1f} {p8:8.1e} | {c16 / h16:17.1f} {p16:6.2f} | '
                     + '  '.join(f'{a:+6.1f}/{b:+6.1f}' for a, b in zip(s1, s2)) + f' | {ac:+.2f}')
    lines += ['', 'Period-4 records give p > 0.05 on the 8-vs-4 test (no false period 8). A steady period 8',
              'rejects period 4, shows nothing at 16, keeps its differences from one half to the other,',
              'and has an autocorrelation near 0.']
    text = '\n'.join(lines)
    print(text)
    with open(os.path.join(HERE, 'period8_check.txt'), 'w') as fh:
        fh.write(text + '\n')
    plot(R)


def plot(R):
    """The four differences m_j - m_(j+4), group by group, for the records in FIGURE."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    report_style.apply()
    unc = pd.read_csv(os.path.join(os.path.dirname(HERE), 'uncertainty', 'forward_uncertainty.csv'))
    unc = unc.set_index('filename').u_r_block_ohm
    fig, axes = plt.subplots(len(FIGURE), 1, sharex=True, figsize=report_style.figsize(aspect=1.25),
                             layout='constrained')
    for ax, (n, label) in zip(axes, FIGURE):
        name = f'trace{n}.csv'
        path = os.path.join(HERE, 'forward', name)
        M, _ = maxima(path)
        t, _ = read_scope(path, ('CH1',), min_samples=100)
        dur = (t[-1] - t[0]) * 1e3
        c = fold(np.asarray(M, float) * 1e3, 8)
        tt = (np.arange(len(c)) + 0.5) * 8 * dur / len(M)      # ms, centre of each group of 8
        for j in range(4):
            ax.plot(tt, c[:, j] - c[:, j + 4], lw=0.6, marker='.', ms=1.5, color=report_style.COLORS[j])
        ax.axhline(0, color='0.5', lw=0.4)
        ax.axvline(dur / 4, color='0.6', lw=0.4, ls=':')
        ax.set_ylim(-110, 110)
        ax.text(0.02, 0.93, rf'record {n}, {R[name]:.2f}({round(unc[name] * 100):d}) $\Omega$: {label}',
                transform=ax.transAxes, va='top', fontsize=7.5)
    axes[-1].set_xlabel('Time in record (ms)')
    fig.supylabel(r'$m_j - m_{j+4}$ within each group of 8 maxima (mV)', fontsize=8)
    report_style.save(fig, os.path.join(HERE, 'period8_check.png'))
    plt.close(fig)


if __name__ == '__main__':
    main()
