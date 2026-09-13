"""
Feigenbaum delta from the BENCH sweeps, with R_3 replaced by the point where
the period-4 orbit turns chaotic.

Why the substitution is needed: the 4 -> 8 splitting is 30-50 mV, against
~10 mV of scatter in the interpolated maxima, so the bench never resolves
period 8. The last structure it sees is period 4, and then broadband chaos.

What the substitution costs: the observable "period-4 goes chaotic" point is
the ACCUMULATION point R_inf of the whole cascade, not R_3. With geometric
scaling R_2 - R_inf = g_2 + g_3 + ... = g_2 * d/(d-1), whereas R_2 - R_3 is
g_2 alone, so

    d_proxy = (R_1 - R_2)/(R_2 - R_inf)  ->  d - 1 = 3.669,

i.e. the proxy converges to delta MINUS ONE. Two corrections are reported:
"+1" (the asymptotic relation above) and a rescaling by the bias that the
identical estimator shows on the simulation, where R_3 and R_inf are both
known. The simulated bias is the larger of the two because the early gaps
have not yet reached the asymptotic ratio (g_1/g_2 = 4.286, not 4.669).

R_1 and R_2 are each taken from the sqrt law s^2 ~ (R_c - R) on the three
points nearest the transition, but CLIPPED to the bracket set by the level
counts (the lowest R still showing 2^(n-1) levels and the highest showing
2^n); where the fit falls outside that bracket -- set 2's R_2, whose three
nearest points are too flat to extrapolate from -- the bracket midpoint is
used instead. The half-bracket width is the uncertainty in both cases.

Chaos onset: set 1 from the Rosenstein exponent changing sign (-10.6 +- 2.2
at 742.3 ohm, +52.3 +- 4.5 at 741.6); set 2 has no Rosenstein column in its
lyapunov file, so from the level structure (clean period 4 down to 750.0,
merged by 749.5, broadband at 749.0).
"""
import csv
import numpy as np
from collections import defaultdict

BASE = "/mnt/user-data/uploads/lab 5/chua"
SETS = {'set 1': 'sweep forward', 'set 2': 'set 2/sweep forward'}
R_CHAOS = {'set 1': (742.0, 0.4), 'set 2': (749.8, 0.5)}
SIM = dict(R1=798.7581, R2=772.8372, R3=766.7894, Rinf=765.115)
D_UNIV = 4.669202


def levels_of(m, tol=0.030, minfrac=0.03):
    m = np.sort(np.asarray(m, float))
    gs = np.split(m, np.flatnonzero(np.diff(m) > tol) + 1)
    return np.array([g.mean() for g in gs if len(g) > minfrac * len(m)])


def newest(L):
    L = np.sort(L)
    return float(np.mean(L[1::2] - L[0::2]))


def load(tag):
    g = defaultdict(list)
    with open(f'{BASE}/{tag}_bifurcation_points.csv') as fh:
        for row in csv.DictReader(fh):
            try:
                g[float(row['rpot_ohm'])].append(float(row['max_v']))
            except (ValueError, KeyError):
                pass
    return {r: levels_of(np.array(v)) for r, v in g.items() if 690 <= r <= 800}


def transition(lev, k, nfit=3):
    below = sorted([r for r, L in lev.items() if len(L) == k], reverse=True)
    above = sorted([r for r, L in lev.items() if len(L) == k // 2])
    lo = below[0]
    hi = min(r for r in above if r > lo)
    pts = [(r, newest(lev[r])) for r in below[:nfit]]
    p = np.polyfit([q[0] for q in pts], [q[1] ** 2 for q in pts], 1)
    fit = -p[1] / p[0]
    ok = lo <= fit <= hi
    return (fit if ok else 0.5 * (lo + hi)), 0.5 * (hi - lo), (lo, hi), ok


print(__doc__)
print(f'{"":7} {"R_1 (ohm)":>26} {"R_2 (ohm)":>26} {"R_chaos":>13}')
res = {}
for name, tag in SETS.items():
    lev = load(tag)
    R1, e1, b1, ok1 = transition(lev, 2)
    R2, e2, b2, ok2 = transition(lev, 4)
    rc, erc = R_CHAOS[name]
    res[name] = (R1, e1, R2, e2, rc, erc)
    f1 = 'fit ' if ok1 else 'mid '
    f2 = 'fit ' if ok2 else 'mid '
    print(f'{name:7} {f1}{R1:7.2f}+-{e1:.2f} [{b1[0]:6.1f},{b1[1]:6.1f}] '
          f'{f2}{R2:7.2f}+-{e2:.2f} [{b2[0]:6.1f},{b2[1]:6.1f}] '
          f'{rc:8.1f}+-{erc:.1f}')

print(f'\n{"":7} {"g_1 = R1-R2":>14} {"R2-R_chaos":>14} {"d_proxy":>16}')
vals, errs = [], []
for name, (R1, e1, R2, e2, rc, erc) in res.items():
    g1, gp = R1 - R2, R2 - rc
    eg1, egp = np.hypot(e1, e2), np.hypot(e2, erc)
    d = g1 / gp
    ed = d * np.hypot(eg1 / g1, egp / gp)
    vals.append(d); errs.append(ed)
    print(f'{name:7} {g1:9.2f}+-{eg1:.2f} {gp:9.2f}+-{egp:.2f} {d:10.2f}+-{ed:.2f}')

v, e = np.array(vals), np.array(errs)
w = 1 / e ** 2
mean = float((v * w).sum() / w.sum())
err = float(1 / np.sqrt(w.sum()))
chi2 = float((((v - mean) / e) ** 2).sum())
if chi2 > len(v) - 1:                      # inflate for the set-to-set spread
    err *= np.sqrt(chi2 / max(len(v) - 1, 1))
print(f'\n   pooled d_proxy = {mean:.2f} +- {err:.2f}   '
      f'(chi2 = {chi2:.1f} for {len(v)-1} dof: the two sweeps disagree by more '
      f'than their\n   internal errors, so the pot calibration differs between '
      f'them; the error is inflated to match)')

dps = (SIM['R1'] - SIM['R2']) / (SIM['R2'] - SIM['Rinf'])
d1s = (SIM['R1'] - SIM['R2']) / (SIM['R2'] - SIM['R3'])
print(f'\nTHE SAME ESTIMATOR ON THE SIMULATION')
print(f'   true       d_1 = g_1/g_2       = {d1s:.3f}')
print(f'   proxy  (R1-R2)/(R2-R_inf)      = {dps:.3f}      <- what the bench '
      f'estimator should return')
print(f'   universal  d                   = {D_UNIV:.3f}')

print(f'\nRESULT')
print(f'   bench  d_proxy                 = {mean:.2f} +- {err:.2f}'
      f'   (model says {dps:.2f})')
print(f'   corrected, "+1"                = {mean + 1:.2f} +- {err:.2f}')
print(f'   corrected, x model bias {D_UNIV/dps:.3f} = '
      f'{mean * D_UNIV / dps:.2f} +- {err * D_UNIV / dps:.2f}')
print(f'   universal delta                = {D_UNIV:.3f}')
