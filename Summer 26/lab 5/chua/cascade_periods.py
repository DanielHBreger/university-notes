"""
Period of each swept record from lag distances, and the bifurcation points
R_1, R_2, R_3 that follow -- including the 4 -> 8 doubling, which IS in the
bench data.

Why this replaces clustering
----------------------------
Grouping a record's maxima into "levels" with a tolerance cannot see the
deep cascade: by period 8 the new splittings are 10-45 mV while a tolerance
loose enough to survive the scatter is ~30 mV, so levels merge and the count
comes out as 3, or 4, or anything. Nothing is wrong with the DATA -- the
splittings are there.

Instead, use the fact that a settled period-P orbit has a maxima sequence
that repeats exactly every P:

    D(p) = < |m_i - m_(i+p)| >

collapses to the record's noise floor at p = P and every multiple, and is
large at every other lag. The period is then the smallest lag at which D
reaches the floor -- no tolerance anywhere. D(P/2) is the splitting created
by the doubling that produced this orbit, and it obeys the same sqrt law
used in feigenbaum.py, so R_n follows by extrapolating D(P/2)^2 to zero.

D(64)/floor is a stationarity check: on a clean record it is ~1, and it grows
when the pot was still moving during acquisition.

The sweep tag names a "<tag>_bifurcation_points.csv" sidecar written by
bifurcation.py beside this script (forward, back). Results for the current
sweeps are kept in cascade_periods_forward.txt and cascade_periods_back.txt;
the model's values come from feigenbaum.py.

Usage:
    python cascade_periods.py forward
    python cascade_periods.py back --range 745 785
    python cascade_periods.py --record forward/trace81.csv   # one raw record
"""
import argparse
import csv
import os
import sys
from collections import OrderedDict

import numpy as np
from uncertainty import delta_covariance

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

LAGS = (1, 2, 4, 8, 16, 32, 64)


def lag_distances(M, drop=0.25):
    """D(p) for each lag, over the last (1-drop) of the maxima sequence."""
    t = np.asarray(M, float)[int(len(M) * drop):]
    return {p: float(np.mean(np.abs(t[:-p] - t[p:]))) for p in LAGS if len(t) > p}


def period_of(D, factor=2.0):
    """
    (period, splitting, floor, drift). The period is the smallest lag whose
    distance has fallen within `factor` of the record's own noise floor; the
    splitting is D at half that lag -- the doubling that created the orbit.
    """
    floor = min(D.values())
    P = next((p for p in LAGS if p in D and D[p] < factor * floor), 0)
    return P, (D[P // 2] if P >= 2 else np.nan), floor, D[max(D)] / floor


def load_sidecar(tag, base='.'):
    """{filename: (rpot, maxima in acquisition order)} from bifurcation.py."""
    d = OrderedDict()
    with open(os.path.join(base, f'{tag}_bifurcation_points.csv')) as fh:
        for r in csv.DictReader(fh):
            try:
                d.setdefault(r['filename'],
                             [float(r['rpot_ohm']), []])[1].append(float(r['max_v']))
            except (ValueError, KeyError):
                pass
    return {k: (v[0], np.array(v[1])) for k, v in d.items()}


def bifurcation_points(rows, nfit=3):
    """
    R_n from the nearest `nfit` records inside each 2^n window, by the sqrt
    law, with the bracket set by the period labels. Returns
    {n: (R, half-bracket, bracket, used_fit, points)}.
    """
    out = {}
    for n, P in ((1, 2), (2, 4), (3, 8)):
        pts = sorted([(r, s) for r, pp, s in rows if pp == P], key=lambda x: -x[0])
        if len(pts) < 2:
            continue
        lo = pts[0][0]
        above = [r for r, pp, s in rows if pp == P // 2 and r > lo]
        hi = min(above) if above else np.nan
        use = pts[:nfit]
        a = np.polyfit([p[0] for p in use], [p[1] ** 2 for p in use], 1)
        Rc = -a[1] / a[0]
        ok = (lo <= Rc <= hi) if np.isfinite(hi) else Rc >= lo
        out[n] = (Rc if ok else 0.5 * (lo + hi), 0.5 * (hi - lo) if np.isfinite(hi)
                  else np.nan, (lo, hi), ok, use)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('sweep', nargs='?', help='sweep tag, e.g. forward')
    p.add_argument('--record', help='analyse one raw scope record instead')
    p.add_argument('--range', type=float, nargs=2, default=[744, 790],
                   metavar=('RMIN', 'RMAX'))
    p.add_argument('--drift-max', type=float, default=2.5,
                   help='reject records whose D(64)/floor exceeds this')
    p.add_argument('--nfit', type=int, default=3)
    a = p.parse_args()

    if a.record:
        import lorenz_map as lm
        M, info = lm.maxima(a.record)
        D = lag_distances(M)
        P, sp, fl, dr = period_of(D)
        print(f'{a.record}: {len(M)} maxima, quantum {info["quantum"]*1e3:.1f} mV')
        print('  ' + '  '.join(f'D({q})={D[q]*1e3:6.1f}' for q in LAGS if q in D))
        print(f'  -> period {P}, newest splitting {sp*1e3:.1f} mV, '
              f'noise floor {fl*1e3:.1f} mV')
        if P >= 8:
            t = np.asarray(M, float)[len(M)//4:]
            k = (len(t)//P)*P
            cyc = t[:k].reshape(-1, P)
            mu, sd = cyc.mean(0)*1e3, cyc.std(0)*1e3
            print(f'  the {P} levels (mV, mean +- sd over {cyc.shape[0]} repeats):')
            for i in np.argsort(mu):
                print(f'     {mu[i]:8.1f} +- {sd[i]:4.1f}')
        return

    recs = load_sidecar(a.sweep)
    rows = []
    print(f'{"file":>13} {"Rpot":>8} {"P":>3} {"floor":>7} {"split":>8} {"drift":>6}')
    for name, (rp, M) in sorted(recs.items(), key=lambda kv: -kv[1][0]):
        if not (a.range[0] <= rp <= a.range[1]):
            continue
        P, sp, fl, dr = period_of(lag_distances(M))
        flag = '' if dr < a.drift_max else '  (pot moving)'
        print(f'{name:>13} {rp:8.2f} {P:3d} {fl*1e3:7.1f} {sp*1e3:8.1f} {dr:6.1f}{flag}')
        if dr < a.drift_max:
            rows.append((rp, P, sp))

    out = bifurcation_points(rows, a.nfit)
    print()
    for n in sorted(out):
        R, e, br, ok, use = out[n]
        print(f'  R_{n} = {R:7.2f}  half-bracket {e:4.2f} ohm  bracket [{br[0]:.2f},{br[1]:.2f}]'
              f'  {"sqrt fit" if ok else "bracket midpoint"}')
    if {1, 2, 3} <= set(out):
        g1, g2 = out[1][0] - out[2][0], out[2][0] - out[3][0]
        e1 = np.hypot(out[1][1], out[2][1]); e2 = np.hypot(out[2][1], out[3][1])
        d, sensitivity = delta_covariance([out[k][0] for k in (1, 2, 3)],
                                          np.diag([out[k][1]**2 for k in (1, 2, 3)]))
        print(f'\n  g_1 = {g1:5.2f}    g_2 = {g2:5.2f}')
        print(f'  delta_1 = {d:.2f}; propagated half-bracket sensitivity {sensitivity:.2f}'
              f'    (universal 4.669; model: feigenbaum.py)')
        print('  Half-brackets describe resistance resolution, not 1-sigma errors.')
        print('  Shared R2 covariance is included; calibration and drift terms are additional.')


if __name__ == '__main__':
    main()
