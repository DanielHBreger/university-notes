"""
Feigenbaum's first ratio from the forward sweep, delta_1 = (R1 - R2)/(R2 - R3).

Each doubling R_n lies between the last record of the old period and the first
record of the new one. R_n is the midpoint of that bracket and its uncertainty
half the bracket, as for every other transition in the report; delta_1 follows
by linear propagation with R2 shared between the two intervals.

A record counts only if its pattern is steady: the period from the lag
distances of the V1 maxima (cascade_periods.py), and each newest splitting
(maxima j against j + P/2) keeping its sign in all four quarters of the record.
Two records flip the phase of their splitting (trace67, trace94; see
period8_check.py); they sit too close to a doubling to say on which side and
are left out. From trace96 down the maxima no longer repeat to within
FLOOR_MAX_MV: chaotic.

The whole cascade (trace61-100) is on one oscilloscope range, so the
range-to-range scale cancels in delta_1.

Outputs: delta_bench.txt, delta_bench.json (R_n for the cascade figure).

Usage:
    python delta_bench.py
"""
import json
import os
from collections import namedtuple

import numpy as np
import pandas as pd

from cascade_periods import lag_distances, period_of
from lorenz_map import maxima
from uncertainty import delta_covariance

HERE = os.path.dirname(os.path.abspath(__file__))
RECORDS = range(58, 101)
FLOOR_MAX_MV = 10.0      # periodic: maxima repeat at lag P to within this (chaotic records: 15 mV and up)
DELTA_UNIVERSAL = 4.6692

# one record of the cascade; tuples sort by resistance first
Record = namedtuple('Record', 'R period steady name chaotic')


def classify(name):
    """(period, steady, chaotic) of one record."""
    M, _ = maxima(os.path.join(HERE, 'forward', name))
    M = np.asarray(M, float)
    P, _, floor, _ = period_of(lag_distances(M))
    if floor * 1e3 > FLOOR_MAX_MV:
        return P, False, True
    if P == 1:
        return P, True, False
    m = M[len(M) // 4:]
    c = m[:(len(m) // P) * P].reshape(-1, P)
    d = c[:, :P // 2] - c[:, P // 2:]                   # newest splittings, one column per pair
    quarters = np.array([q.mean(0) for q in np.array_split(d, 4)])
    return P, bool(np.all(np.sign(quarters) == np.sign(d.mean(0)))), False


def describe(rec):
    """One line of the record table."""
    if rec.chaotic:
        state = 'chaotic'
    else:
        state = f'period {rec.period}' + ('' if rec.steady else ', phase flips: left out')
    return f'  {rec.name:>13} {rec.R:7.2f} ohm  ' + state


def main():
    R = pd.read_csv(os.path.join(HERE, 'forward_rpot.csv')).set_index('filename').rpot_ohm
    recs = []
    for n in RECORDS:
        name = f'trace{n}.csv'
        P, steady, chaotic = classify(name)
        recs.append(Record(float(R[name]), P, steady, name, chaotic))
    recs.sort(reverse=True)

    lines = ['Forward records 58-100 (26.5 mV-step range), by resistance:']
    lines += [describe(rec) for rec in recs]
    lines.append('')
    res = {}
    for n in (1, 2, 3):
        P = 2 ** n
        new = max((x for x in recs if x.steady and x.period == P), key=lambda x: x.R)       # first of the new period
        old = min((x for x in recs if x.steady and x.period == P // 2 and x.R > new.R), key=lambda x: x.R)
        mid, half = 0.5 * (old.R + new.R), 0.5 * (old.R - new.R)
        res[f'R{n}'] = dict(value=mid, u=half)
        lines.append(f'R{n} (period {P // 2} -> {P}) between {old.name} ({old.R:.2f}) and {new.name} ({new.R:.2f}): '
                     f'{mid:.2f} +- {half:.2f} ohm')
    Rs = [res[f'R{n}']['value'] for n in (1, 2, 3)]
    us = [res[f'R{n}']['u'] for n in (1, 2, 3)]
    d, ud = delta_covariance(Rs, np.diag(np.square(us)))
    res['delta_1'] = dict(value=float(d), u=float(ud))
    # accumulation point if the intervals keep shrinking by the universal ratio
    k = 1 / (DELTA_UNIVERSAL - 1)
    rinf = Rs[2] - k * (Rs[1] - Rs[2])
    urinf = float(np.hypot((1 + k) * us[2], k * us[1]))
    res['R_inf_universal'] = dict(value=float(rinf), u=urinf)
    lines += ['', f'R1 - R2 = {Rs[0] - Rs[1]:.2f} +- {np.hypot(us[0], us[1]):.2f} ohm, '
                  f'R2 - R3 = {Rs[1] - Rs[2]:.2f} +- {np.hypot(us[1], us[2]):.2f} ohm',
              f'delta_1 = {d:.2f} +- {ud:.2f} (half-brackets propagated, R2 shared); universal 4.669 is '
              f'{(DELTA_UNIVERSAL - d) / ud:.1f} of these units away',
              f'accumulation point with the universal ratio: {rinf:.2f} +- {urinf:.2f} ohm']
    text = '\n'.join(lines)
    print(text)
    with open(os.path.join(HERE, 'delta_bench.txt'), 'w') as fh:
        fh.write(text + '\n')
    with open(os.path.join(HERE, 'delta_bench.json'), 'w') as fh:
        json.dump(res, fh, indent=1)


if __name__ == '__main__':
    main()
