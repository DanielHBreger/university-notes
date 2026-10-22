"""Sagnac phase against angular velocity (task 3), the speed of light, and the uncertainties (task 5).

fringes.py gives the phase at the ring centre of every photograph (and a second reading by a
cosine fit); methods.py adds template matching; calibration.py gives omega(S). This script
averages the photographs of each setting, fits the line Delta phi = a omega + b, compares a with
the manual's Eq. (5), solves Eq. (5) for the speed of light and writes every number with its
uncertainty.   python sagnac.py   (from this folder, after calibration.py, fringes.py, methods.py)

Steps (main(), top to bottom)
  1. setting_means  per setting: n photos, mean phase, SEM = SD/sqrt(n), median ring spacing P.
  2. omega          from the calibration line, exactly 0 at rest, with the prediction uncertainty
                    (calibration.py: the line's uncertainty and the scatter of the timed settings
                    about the line, in quadrature).
  3. regress        Delta phi = a omega + b by ordinary least squares through the ten setting
                    means (Schwartz 2017, Eq. 7). Every point's uncertainty is s, the rms distance
                    of the points from the line, and u(a) = s / sqrt(sum (omega - mean)^2).
  4. checks         the same line weighted by the photos' SEM (chi2 >> dof: the settings differ by
                    more than the photo noise, which is why the regression's own s is used); an
                    added omega^2 term; the line forced through the origin.
  5. methods        the cosine fit and template matching through the same regression; half the
                    range of the three slopes is the method uncertainty of a.
  6. calibration    every moving omega comes from the same line and the same table, so the relative
                    prediction uncertainty, sqrt((u(k)/k)^2 + rel_rms^2), carries straight into a.
  7. result         a with the three parts in quadrature; a_th = 8 pi^2 N R^2 / (lambda c), the
                    comparison, N implied by the slope, the largest shift and its travel time.
  8. c              Eq. (5) solved for c; the relative parts of a, N and R in quadrature.
  9. drift          Delta phi = a omega + D t + b: can a drift in time be told from the slope?
 10. outputs        sagnac_speeds.csv (per setting), sagnac.txt (words and report form), sagnac.json.

Phases are angles, known modulo 2 pi: within a setting each photo is brought within pi of the
setting's circular mean before averaging, and each setting mean within pi of the previous one.
Here neither step changes a number (the whole run spans 3.3 rad); they are safeguards.
"""
import json
import time
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from calibration import omega_of, u_omega_pred

HERE = Path(__file__).resolve().parent

# Fibre loop, centred on the rotation axis (given in the manual; +- taken as standard uncertainties)
N_TURNS, U_N_TURNS = 13, 1
RADIUS, U_RADIUS = 0.570, 0.005       # m
WAVELENGTH = 640e-9                   # m
U_WAVELENGTH = None                   # not given by the lab: left out of the budget, stated in the output
C_LIGHT = 299_792_458.0               # m/s, exact
OTHER_METHODS = [('cosine fit', 'phi_cos_rad'), ('template matching', 'phi_template_rad')]
SCAN = HERE / 'presentation_slides' / 'source' / 'random_points_scan.json'   # optional, prep Q5
SCAN_ORIGINAL_SEED = 20261011


def wrap(x):
    """Angle(s) wrapped into (-pi, pi]."""
    return np.angle(np.exp(1j * np.asarray(x)))


def round1(x, u):
    """'x +- u' with u to one significant figure and x to the same digit, rounding half up on the
    decimal digits as printed (float round() sends ties to even: 2.25 -> 2.2). Decimal works on
    the printed digits, so there are no binary rounding surprises."""
    if not np.isfinite(u) or u <= 0:
        return f'{x:.4g} +- {u}'
    U = Decimal(str(float(u)))
    q = Decimal(1).scaleb(U.adjusted())                  # 10^(position of u's first digit)
    Ur = U.quantize(q, ROUND_HALF_UP)
    if Ur.adjusted() > U.adjusted():                     # 0.096 -> 0.1 moves the digit
        q = q.scaleb(1)
        Ur = U.quantize(q, ROUND_HALF_UP)
    X = Decimal(str(float(x))).quantize(q, ROUND_HALF_UP)
    return f'{abs(X) if X == 0 else X:f} +- {Ur:f}'      # abs(): no "-0.00"


def theory():
    """The predicted slope a_th = 8 pi^2 N R^2 / (lambda c) (rad per rad/s) and its uncertainty
    components: a_th is proportional to N and to R^2, so the relative parts are u_N/N and 2 u_R/R."""
    a = 8 * np.pi ** 2 * N_TURNS * RADIUS ** 2 / (WAVELENGTH * C_LIGHT)
    comps = {'N': a * U_N_TURNS / N_TURNS, 'R': a * 2 * U_RADIUS / RADIUS}
    if U_WAVELENGTH is not None:
        comps['lambda'] = a * U_WAVELENGTH / WAVELENGTH
    return a, comps


def speed_of_light(a, u_a_parts):
    """Eq. (5) solved for c with the measured slope a: c = 8 pi^2 N R^2 / (lambda a). u_a_parts: the
    slope's uncertainty parts by name. Returns c (m/s) and its components (m/s): the slope's parts
    keep their relative size (c ~ 1/a), then N and R."""
    c = 8 * np.pi ** 2 * N_TURNS * RADIUS ** 2 / (WAVELENGTH * a)
    comps = {k: c * u / a for k, u in u_a_parts.items()}
    comps['N'] = c * U_N_TURNS / N_TURNS
    comps['R'] = c * 2 * U_RADIUS / RADIUS
    if U_WAVELENGTH is not None:
        comps['lambda'] = c * U_WAVELENGTH / WAVELENGTH
    return c, comps


def regress(x, y):
    """Step 3. Ordinary least squares y = a x + b, equal weights. The uncertainty of every point is
    taken from the points themselves: s = sqrt(sum r^2 / (n - 2)) is the rms distance from the
    line, and the covariance of (a, b) is s^2 (A^T A)^-1. Returns a dict with a, b, u_a, u_b, s,
    R^2 and the residuals."""
    A = np.c_[x, np.ones_like(x)]
    q, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ q
    s2 = r @ r / (len(y) - 2)
    cov = s2 * np.linalg.inv(A.T @ A)
    return dict(a=float(q[0]), b=float(q[1]), u_a=float(np.sqrt(cov[0, 0])), u_b=float(np.sqrt(cov[1, 1])),
                s=float(np.sqrt(s2)), r2=float(1 - r @ r / np.sum((y - y.mean()) ** 2)), resid=r)


def setting_means(fr, col):
    """Step 1. One row per setting from the per-photo phases in column col: n, mean phase phi, its
    SEM = SD/sqrt(n), the time of the first photo, the median ring spacing P, and dphi = phi minus
    the rest mean. See the module docstring for the 2 pi safeguards."""
    rows = []
    for S, g in fr.groupby('S', sort=True):
        ph = g[col].to_numpy()
        c = np.angle(np.mean(np.exp(1j * ph)))           # circular mean: angle of the mean unit vector
        ph = c + wrap(ph - c)                            # each photo within pi of it
        rows.append(dict(S=S, n=len(ph), phi=ph.mean(), sem=ph.std(ddof=1) / np.sqrt(len(ph)),
                         t_first=g.time_s.min(), P=g.P_px2.median()))
    tab = pd.DataFrame(rows)
    phi = tab.phi.to_numpy().copy()
    for i in range(1, len(phi)):                         # each setting within pi of the previous one
        phi[i] += 2 * np.pi * np.round((phi[i - 1] - phi[i]) / (2 * np.pi))
    tab['phi'] = phi
    tab['dphi'] = phi - phi[0]
    return tab


def hhmmss(t):
    return time.strftime('%H:%M:%S', time.gmtime(t))


def scan_summary(path=SCAN):
    """Optional (prep Q5): the random-pixel method of the presentation stage, repeated for many
    draws of n fixed pixels (random_points_scan.json). Per n: the number of draws, the mean and SD
    of their slopes (failed fits with a negative slope left out) and the median fit error of one
    draw. The slopes in that file were fitted by the presentation's weighted line, not by regress()."""
    if not path.exists():
        return None
    draws = [d for d in json.loads(path.read_text())['draws'] if d['seed'] != SCAN_ORIGINAL_SEED]
    out = {}
    for n in sorted({d['n'] for d in draws}):
        sl = np.array([d['fixed']['slope'] for d in draws if d['n'] == n])
        ok = sl[sl > 0]
        out[int(n)] = dict(draws=len(sl), failed=int(np.sum(sl <= 0)), mean=float(ok.mean()),
                           sd=float(ok.std(ddof=1)), u_fit=float(np.median([d['fixed']['u_fit'] for d in draws if d['n'] == n])))
    return out


def main():
    cal = json.loads((HERE / 'calibration.json').read_text())
    fr = pd.read_csv(HERE / 'fringes_frames.csv')
    mfile = HERE / 'methods_frames.csv'
    if mfile.exists():
        fr = fr.merge(pd.read_csv(mfile)[['file', 'phi_template_rad']], on='file', how='left', validate='one_to_one')

    # 1-2. Setting means; omega from the calibration line (0 at rest), its uncertainty the line's scatter.
    tab = setting_means(fr, 'phi_rad')
    S = tab.S.to_numpy()
    moving = S > 0
    omega = np.where(moving, omega_of(S, cal['line']['coef']), 0.0)
    u_omega = np.where(moving, u_omega_pred(S, cal['line']), 0.0)
    y, sem = tab.dphi.to_numpy(), tab['sem'].to_numpy()      # tab.sem would be pandas' sem() method
    n, dof = len(y), len(y) - 2

    # 3. The regression.
    reg = regress(omega, y)
    a, b, u_stat, s = reg['a'], reg['b'], reg['u_a'], reg['s']

    # 4. Checks: SEM-weighted line, an added omega^2 term, the line through the origin, the intercept.
    A = np.c_[omega, np.ones(n)]
    W = 1 / sem ** 2
    qw = np.linalg.inv(A.T @ (A * W[:, None])) @ (A.T @ (W * y))
    chi_sem = float(np.sum(((y - A @ qw) / sem) ** 2))
    p_sem = float(stats.chi2.sf(chi_sem, dof))
    A2 = np.c_[omega ** 2, omega, np.ones(n)]
    q2, *_ = np.linalg.lstsq(A2, y, rcond=None)
    r2 = y - A2 @ q2
    C2 = (r2 @ r2 / (n - 3)) * np.linalg.inv(A2.T @ A2)
    t_q = float(q2[0] / np.sqrt(C2[0, 0]))
    p_q = float(2 * stats.t.sf(abs(t_q), n - 3))
    a0 = float((omega @ y) / (omega @ omega))
    r0 = y - a0 * omega
    u_a0 = float(np.sqrt(r0 @ r0 / (n - 1) / (omega @ omega)))
    p_b = float(2 * stats.t.sf(abs(b / reg['u_b']), dof))

    # 5. The other phase-reading methods through the same regression; half the range = method part.
    methods = {}
    for label, col in OTHER_METHODS:
        if col in fr:
            t = setting_means(fr, col)
            methods[label] = dict(tab=t, reg=regress(omega, t.dphi.to_numpy()), col=col)
    slopes = {'ring minima': a, **{l: m['reg']['a'] for l, m in methods.items()}}
    u_meth = (max(slopes.values()) - min(slopes.values())) / 2

    # 6. The calibration part: every moving omega is predicted by the same line from the same table,
    #    so the relative prediction uncertainty (k's relative uncertainty and the scatter of the
    #    timed settings about the line) is a common scale error of omega, hence of a. (omega0
    #    shifts the moving points together and is absorbed by the intercept, up to the rest point.)
    u_k_rel = cal['line']['u_k'] / cal['line']['k']
    rel_cal = float(np.hypot(u_k_rel, cal['line']['rel_rms']))
    u_cal = a * rel_cal

    # 7. The result and the comparison with Eq. (5).
    u_a = float(np.sqrt(u_stat ** 2 + u_meth ** 2 + u_cal ** 2))
    a_th, comps = theory()
    u_th = float(np.sqrt(sum(v ** 2 for v in comps.values())))
    z = (a - a_th) / np.hypot(u_a, u_th)
    p_z = float(2 * stats.norm.sf(abs(z)))
    z_stat = (a - a_th) / np.hypot(u_stat, u_th)
    ratio, u_ratio = a / a_th, a / a_th * np.hypot(u_a / a, u_th / a_th)
    N_impl = N_TURNS * ratio
    u_N_impl = N_impl * np.hypot(u_a / a, 2 * U_RADIUS / RADIUS)
    dphi_max, u_max = float(y[-1]), float(np.sqrt(2) * s)          # top setting minus rest, each known to s
    dt = dphi_max * WAVELENGTH / (2 * np.pi * C_LIGHT)
    dt_pred = a_th * omega[-1] * WAVELENGTH / (2 * np.pi * C_LIGHT)

    # 8. The speed of light.
    c_parts = {'regression': u_stat, 'method': u_meth, 'calibration': u_cal}
    c_m, c_comps = speed_of_light(a, c_parts)
    u_c = float(np.sqrt(sum(v ** 2 for v in c_comps.values())))
    z_c = (c_m - C_LIGHT) / u_c
    p_c = float(2 * stats.norm.sf(abs(z_c)))

    # 9. Drift: a phase linear in time fitted together with the slope (settings taken in increasing order).
    tm = (fr.groupby('S').time_s.mean() - fr.time_s.min()).to_numpy()
    Ad = np.c_[omega, tm, np.ones(n)]
    qd, *_ = np.linalg.lstsq(Ad, y, rcond=None)
    rd = y - Ad @ qd
    Cd = (rd @ rd / (n - 3)) * np.linalg.inv(Ad.T @ Ad)
    rho = float(Cd[0, 1] / np.sqrt(Cd[0, 0] * Cd[1, 1]))
    t_D = float(qd[1] / np.sqrt(Cd[1, 1]))
    p_D = float(2 * stats.t.sf(abs(t_D), n - 3))
    scan = scan_summary()

    # 10. Outputs. The per-setting table first.
    tab['omega'] = omega; tab['u_omega'] = u_omega; tab['resid'] = reg['resid']
    for label, m in methods.items():
        key = m['col'].removeprefix('phi_').removesuffix('_rad')
        tab[f'dphi_{key}'] = m['tab'].dphi; tab[f'sem_{key}'] = m['tab']['sem']
    tab.to_csv(HERE / 'sagnac_speeds.csv', index=False, float_format='%.6g')

    E8 = 1e8
    L = ['Sagnac phase against angular velocity (sagnac.py), fibre loop centred on the rotation axis\n',
         'PER SETTING: phase at the ring centre (minima), mean over the photos, SEM = SD/sqrt(n); omega from the',
         'calibration line +- its prediction uncertainty (common to all settings); dphi = phase minus the rest mean.',
         f"{'S':>6} {'n':>3} {'first photo':>11} {'omega (rad/s)':>18} {'dphi (rad)':>18} {'P (10^5 px^2)':>13} {'resid':>7}"]
    for r in tab.itertuples():
        L.append(f'{r.S:6.2f} {r.n:3d} {hhmmss(r.t_first):>11} {r.omega:8.4f} +- {r.u_omega:.4f} '
                 f'{r.dphi:8.3f} +- {r.sem:.3f} {r.P / 1e5:13.3f} {r.resid:+7.3f}')
    L.append(f'(the rest mean, SEM {sem[0]:.3f} rad, is subtracted from every row; resid = dphi minus the line)\n')
    L.append('THE LINE dphi = a omega + b, ordinary least squares through the ten setting means (Schwartz Eq. 7):')
    L.append(f'  a = {a:.4f} +- {u_stat:.4f} rad/(rad/s) (regression), b = {b:.3f} +- {reg["u_b"]:.3f} rad '
             f'(p = {p_b:.2f} for b = 0), R^2 = {reg["r2"]:.4f}')
    L.append(f'  s = {s:.4f} rad: the rms distance of the points from the line, used as every point\'s uncertainty.')
    L.append(f'  Weighted by the photos\' SEM instead: chi2 = {chi_sem:.1f}, dof = {dof}, p = {p_sem:.1g} (slope {qw[0]:.4f}): '
             'the settings differ by')
    L.append(f'  {np.sqrt(chi_sem / dof):.1f} times more than the photo noise allows (the phase wanders between settings), '
             'so the regression\'s s is used.')
    L.append(f'  Added omega^2 term: {q2[0]:+.3f} +- {np.sqrt(C2[0, 0]):.3f} rad/(rad/s)^2, t = {t_q:.2f}, p = {p_q:.2f}: no bend.')
    L.append(f'  Line through the origin (check): a = {a0:.4f} +- {u_a0:.4f}.\n')
    L.append('OTHER PHASE-READING METHODS, same averaging and regression (cosine fit: fringes.py, prep Q5; template')
    L.append('matching: methods.py):')
    for label, m in methods.items():
        rg = m['reg']
        L.append(f'  {label:18s} a = {rg["a"]:.4f} +- {rg["u_a"]:.4f}, s = {rg["s"]:.4f} rad, R^2 = {rg["r2"]:.4f}; '
                 'minus the minima per setting: ' + ', '.join(f'{d:+.2f}' for d in m['tab'].dphi - y) + ' rad')
    L.append(f'  Method part of u(a): half the range of {", ".join(f"{v:.4f}" for v in slopes.values())} = {u_meth:.4f} '
             f'({100 * u_meth / a:.1f} %).')
    L.append(f'  Calibration part of u(a): k is known to {100 * u_k_rel:.2f} % and the timed settings scatter about the '
             f'line by {100 * cal["line"]["rel_rms"]:.2f} % (calibration.txt);')
    L.append(f'  both are common to every omega, so u(a)/a = sqrt({100 * u_k_rel:.2f}^2 + {100 * cal["line"]["rel_rms"]:.2f}^2) % '
             f'= {100 * rel_cal:.2f} %, u = {u_cal:.4f}.\n')
    L.append(f'RESULT  a = {a:.4f} +- {u_stat:.4f} (regression) +- {u_meth:.4f} (method) +- {u_cal:.4f} (calibration)')
    L.append(f'          = {a:.4f} +- {u_a:.4f} rad/(rad/s) ({100 * u_a / a:.1f} %, the three parts in quadrature)')
    L.append(f"Theory, Eq. (5): N = {N_TURNS} +- {U_N_TURNS}, R = {100 * RADIUS:.1f} +- {100 * U_RADIUS:.1f} cm "
             f"(A = pi R^2 = {np.pi * RADIUS ** 2:.4f} m^2), lambda = {WAVELENGTH * 1e9:.0f} nm, c exact.")
    L.append(f'  a_th = 8 pi^2 N R^2 / (lambda c) = {a_th:.4f} +- {u_th:.4f} rad/(rad/s) '
             f'[N: {comps["N"]:.4f} ({100 * U_N_TURNS / N_TURNS:.1f} %), R: {comps["R"]:.4f} ({200 * U_RADIUS / RADIUS:.1f} %)]')
    L.append('  lambda: no uncertainty given, NOT included (unknown, not zero); a_th and c move by '
             f'{100e-9 / WAVELENGTH:.2f} % per nm.')
    L.append(f'  a / a_th = {ratio:.3f} +- {u_ratio:.3f}; a - a_th = {a - a_th:+.4f} = {z:+.2f} combined sigma, p = {p_z:.2f} '
             f'(regression part only: {z_stat:+.2f} sigma, p = {2 * stats.norm.sf(abs(z_stat)):.2f}).')
    L.append(f'  Read the other way (R, lambda, c known): N = {N_impl:.2f} +- {u_N_impl:.2f} turns.')
    L.append(f'  Largest shift (S = 90 minus rest): {dphi_max:.3f} +- {u_max:.3f} rad = {dphi_max / (2 * np.pi):.3f} +- '
             f'{u_max / (2 * np.pi):.3f} fringe at omega = {omega[-1]:.3f} rad/s; predicted {a_th * omega[-1]:.2f} rad.')
    L.append(f'  As a travel-time difference, dt = dphi lambda / (2 pi c) = {dt * 1e15:.2f} +- {dt * 1e15 * u_max / dphi_max:.2f} fs '
             f'(predicted {dt_pred * 1e15:.2f} fs; one optical period is {WAVELENGTH / C_LIGHT * 1e15:.2f} fs).\n')
    L.append('SPEED OF LIGHT, Eq. (5) solved for c = 8 pi^2 N R^2 / (lambda a), in 10^8 m/s:')
    L.append(f'  c = {c_m / E8:.4f}' + ''.join(f' +- {v / E8:.4f} ({k}, {100 * v / c_m:.1f} %)' for k, v in c_comps.items()))
    L.append(f'    = {c_m / E8:.4f} +- {u_c / E8:.4f} ({100 * u_c / c_m:.1f} %); c / c_def = {c_m / C_LIGHT:.3f} +- {u_c / C_LIGHT:.3f}, '
             f'{z_c:+.2f} sigma, p = {p_c:.2f}.')
    L.append('  This is a / a_th read the other way (c / c_def = a_th / a), not an independent test.\n')
    L.append('DRIFT CHECK. The settings were photographed in increasing order with rest photos only at the start, so')
    L.append('omega and time rise together. dphi = a omega + D t + b (t: mean photo time of the setting):')
    L.append(f'  a = {qd[0]:.2f} +- {np.sqrt(Cd[0, 0]):.2f} rad/(rad/s), D = {1e3 * qd[1]:+.2f} +- {1e3 * np.sqrt(Cd[1, 1]):.2f} mrad/s '
             f'(p = {p_D:.2f}), correlation(a, D) = {rho:.3f}: a steady drift cannot be told from the slope.')
    L.append('  Not in the quoted uncertainty. Also not checked: the sign of the shift (rotation direction not recorded).\n')
    if scan:
        L.append('RANDOM PIXELS (prep Q5, optional; presentation stage, random_points_scan.json):')
        for k, v in scan.items():
            L.append(f'  {k:3d} pixels per draw: {v["draws"]} draws, slope mean {v["mean"]:.3f}, SD between draws {v["sd"]:.3f}, '
                     f'fit error of one draw {v["u_fit"]:.3f}' + (f', {v["failed"]} failed' if v['failed'] else ''))
        L.append('')
    L.append('REPORT FORM (uncertainty to one significant figure, value to the same digit; comparisons above use')
    L.append('unrounded values):')
    L.append(f'  k = {round1(cal["line"]["k"], cal["line"]["u_k"])} rad/s per unit S, omega0 = '
             f'{round1(cal["line"]["omega0"], cal["line"]["u_omega0"])} rad/s')
    L.append(f'  a = {round1(a, u_a)} rad/(rad/s) (all parts); {round1(a, u_stat)} (regression only)')
    L.append(f'  a_th = {round1(a_th, u_th)} rad/(rad/s); a / a_th = {round1(ratio, u_ratio)}; N implied = {round1(N_impl, u_N_impl)}')
    L.append('  ' + '; '.join(f'{l}: a = {round1(m["reg"]["a"], m["reg"]["u_a"])}' for l, m in methods.items()))
    L.append(f'  largest shift = {round1(dphi_max, u_max)} rad = {round1(dphi_max / (2 * np.pi), u_max / (2 * np.pi))} fringe, '
             f'dt = {round1(dt * 1e15, dt * 1e15 * u_max / dphi_max)} fs')
    L.append(f'  c = ({round1(c_m / E8, u_c / E8)}) x 10^8 m/s; c / c_def = {round1(c_m / C_LIGHT, u_c / C_LIGHT)}')
    (HERE / 'sagnac.txt').write_text('\n'.join(L) + '\n', encoding='utf-8')
    print('\n'.join(L))

    out = dict(
        slope=dict(a=a, u_stat=u_stat, u_method=u_meth, u_calibration=u_cal, u=u_a, b=b, u_b=reg['u_b'], p_b=p_b,
                   s=s, r2=reg['r2'], dof=dof, chi2_sem_weighted=chi_sem, p_sem_weighted=p_sem,
                   omega2_term=float(q2[0]), u_omega2_term=float(np.sqrt(C2[0, 0])), t_omega2=t_q, p_omega2=p_q,
                   a_through_origin=a0, u_a_through_origin=u_a0, u_k_rel=u_k_rel, rel_calibration=rel_cal),
        methods={l: dict(a=m['reg']['a'], u_a=m['reg']['u_a'], s=m['reg']['s'], r2=m['reg']['r2']) for l, m in methods.items()},
        theory=dict(a_th=a_th, u_a_th=u_th, components=comps, u_wavelength='not given',
                    per_nm=1e-9 / WAVELENGTH),
        comparison=dict(ratio=ratio, u_ratio=u_ratio, z=float(z), p=p_z, z_regression_only=float(z_stat),
                        N_implied=N_impl, u_N_implied=u_N_impl),
        largest_shift=dict(dphi=dphi_max, u=u_max, fringes=dphi_max / (2 * np.pi), dt_fs=dt * 1e15,
                           u_dt_fs=dt * 1e15 * u_max / dphi_max, dt_predicted_fs=dt_pred * 1e15,
                           omega=float(omega[-1]), dphi_predicted=a_th * float(omega[-1])),
        speed_of_light=dict(c=c_m, u_c=u_c, components=c_comps, rel=100 * u_c / c_m, ratio_to_defined=c_m / C_LIGHT,
                            u_ratio_to_defined=u_c / C_LIGHT, z=float(z_c), p=p_c),
        drift=dict(a=float(qd[0]), u_a=float(np.sqrt(Cd[0, 0])), D=float(qd[1]), u_D=float(np.sqrt(Cd[1, 1])),
                   p_D=p_D, correlation=rho),
        random_pixels=scan,
        report_form=L[L.index('unrounded values):') + 1:])
    (HERE / 'sagnac.json').write_text(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
