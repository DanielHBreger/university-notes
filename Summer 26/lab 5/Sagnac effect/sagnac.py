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
  3. linregress     Delta phi = a omega + b through the ten setting means (Schwartz 2017, Eq. 7).
                    The slope's uncertainty comes from how far the points lie from the line.
  4. checks         the line against the photos' SEM (chi2 >> dof: the settings differ by more than
                    the photo noise, which is why the regression's own scatter is used); an added
                    omega^2 term.
  5. methods        the cosine fit and template matching through the same regression; half the
                    range of the three slopes is the method uncertainty of a.
  6. calibration    every moving omega comes from the same line and the same table, so the relative
                    prediction uncertainty, sqrt((u(k)/k)^2 + rel_rms^2), carries straight into a.
  7. result         a with the three parts in quadrature; a_th = 8 pi^2 N R^2 / (lambda c), the
                    comparison, N implied by the slope, the largest shift and its travel time.
  8. c              Eq. (5) solved for c; the relative parts of a, N and R in quadrature.
  9. drift          Delta phi = a omega + D t + b: can a drift in time be told from the slope?
 10. outputs        sagnac_speeds.csv (per setting), sagnac.txt (words and report form), sagnac.json.
"""
import json
import time
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

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
SCAN_ORIGINAL_SEED = 20261011         # the presentation's main draw, left out of the statistics


def round1(x, u):
    """'x +- u' with u to one significant figure and x to the same digit, rounding half up on the
    printed digits (float round() sends ties to even: 2.25 -> 2.2; Decimal works on the printed
    digits, so there are no binary rounding surprises)."""
    if not np.isfinite(u) or u <= 0:
        return f'{x:.4g} +- {u}'
    u_dec = Decimal(str(float(u)))
    digit = Decimal(1).scaleb(u_dec.adjusted())                 # 10^(position of u's first digit)
    u_rounded = u_dec.quantize(digit, ROUND_HALF_UP)
    if u_rounded.adjusted() > u_dec.adjusted():                 # 0.096 -> 0.1 moves the digit
        digit = digit.scaleb(1)
        u_rounded = u_dec.quantize(digit, ROUND_HALF_UP)
    x_rounded = Decimal(str(float(x))).quantize(digit, ROUND_HALF_UP)
    return f'{abs(x_rounded) if x_rounded == 0 else x_rounded:f} +- {u_rounded:f}'      # abs(): no "-0.00"


def theory():
    """The predicted slope a_th = 8 pi^2 N R^2 / (lambda c) (rad per rad/s) and its uncertainty
    parts: a_th is proportional to N and to R^2, so the relative parts are u_N/N and 2 u_R/R."""
    a_th = 8 * np.pi ** 2 * N_TURNS * RADIUS ** 2 / (WAVELENGTH * C_LIGHT)
    parts = {'N': a_th * U_N_TURNS / N_TURNS, 'R': a_th * 2 * U_RADIUS / RADIUS}
    if U_WAVELENGTH is not None:
        parts['lambda'] = a_th * U_WAVELENGTH / WAVELENGTH
    return a_th, parts


def speed_of_light(a, u_a_parts):
    """Eq. (5) solved for c with the measured slope a: c = 8 pi^2 N R^2 / (lambda a). u_a_parts: the
    slope's uncertainty parts by name. Returns c (m/s) and its parts (m/s): the slope's parts keep
    their relative size (c ~ 1/a), then N and R."""
    c = 8 * np.pi ** 2 * N_TURNS * RADIUS ** 2 / (WAVELENGTH * a)
    parts = {name: c * u / a for name, u in u_a_parts.items()}
    parts['N'] = c * U_N_TURNS / N_TURNS
    parts['R'] = c * 2 * U_RADIUS / RADIUS
    if U_WAVELENGTH is not None:
        parts['lambda'] = c * U_WAVELENGTH / WAVELENGTH
    return c, parts


def setting_means(photos, phase_col):
    """Step 1. One row per setting from the per-photo phases in phase_col: n, mean phase phi, its
    SEM = SD/sqrt(n), the time of the first photo, the median ring spacing P, and dphi = phi minus
    the rest mean. Phases are angles: the setting means are unwrapped in order of S (here no 2 pi
    jump occurs; it is a safeguard)."""
    groups = photos.groupby('S', sort=True)
    means = pd.DataFrame(dict(n=groups[phase_col].size(), phi=groups[phase_col].mean(), sem=groups[phase_col].sem(),
                              t_first=groups.time_s.min(), P=groups.P_px2.median())).reset_index()
    means['phi'] = np.unwrap(means.phi.to_numpy())
    means['dphi'] = means.phi - means.phi[0]
    return means


def rms_from_line(x, y, fit):
    """The rms distance of the points from a linregress line, with n - 2 degrees of freedom."""
    residual = y - (fit.slope * x + fit.intercept)
    return float(np.sqrt(residual @ residual / (len(y) - 2)))


def drift_model(omega_and_t, a, D, b):
    omega, t = omega_and_t
    return a * omega + D * t + b


def hhmmss(seconds):
    return time.strftime('%H:%M:%S', time.gmtime(seconds))


def random_pixel_summary(path=SCAN):
    """Optional (prep Q5): the random-pixel method of the presentation stage, repeated for many
    draws of n fixed pixels (random_points_scan.json). Per n: the number of draws, the mean and SD
    of their slopes (failed fits with a negative slope left out) and the median fit error of one
    draw. The slopes in that file were fitted by the presentation's weighted line, not by linregress."""
    if not path.exists():
        return None
    draws = [d for d in json.loads(path.read_text())['draws'] if d['seed'] != SCAN_ORIGINAL_SEED]
    summary = {}
    for n in sorted({d['n'] for d in draws}):
        slopes = np.array([d['fixed']['slope'] for d in draws if d['n'] == n])
        good = slopes[slopes > 0]
        summary[int(n)] = dict(draws=len(slopes), failed=int(np.sum(slopes <= 0)), mean=float(good.mean()),
                               sd=float(good.std(ddof=1)),
                               u_fit=float(np.median([d['fixed']['u_fit'] for d in draws if d['n'] == n])))
    return summary


def main():
    cal_line = json.loads((HERE / 'calibration.json').read_text())['line']
    photos = pd.read_csv(HERE / 'fringes_frames.csv')
    template_file = HERE / 'methods_frames.csv'
    if template_file.exists():
        photos = photos.merge(pd.read_csv(template_file)[['file', 'phi_template_rad']], on='file', how='left',
                              validate='one_to_one')

    # 1-2. Setting means; omega from the calibration line (0 at rest) with its prediction uncertainty.
    settings = setting_means(photos, 'phi_rad')
    S = settings.S.to_numpy()
    moving = S > 0
    omega = np.where(moving, omega_of(S, cal_line), 0.0)
    u_omega = np.where(moving, u_omega_pred(S, cal_line), 0.0)
    dphi, sem = settings.dphi.to_numpy(), settings['sem'].to_numpy()      # settings.sem would be pandas' sem() method
    n_settings, dof = len(dphi), len(dphi) - 2

    # 3. The line through the ten setting means.
    fit = stats.linregress(omega, dphi)
    a, b, u_stat = fit.slope, fit.intercept, fit.stderr
    rms = rms_from_line(omega, dphi, fit)
    residual = dphi - (a * omega + b)
    r2 = fit.rvalue ** 2
    p_b = float(2 * stats.t.sf(abs(b / fit.intercept_stderr), dof))

    # 4. Checks: the line against the photos' SEM, and an added omega^2 term.
    chi2_sem = float(np.sum((residual / sem) ** 2))
    p_sem = float(stats.chi2.sf(chi2_sem, dof))
    quad, quad_cov = np.polyfit(omega, dphi, 2, cov=True)
    t_quad = float(quad[0] / np.sqrt(quad_cov[0, 0]))
    p_quad = float(2 * stats.t.sf(abs(t_quad), n_settings - 3))

    # 5. The other phase-reading methods through the same regression; half the range = method part.
    methods = {}
    for name, col in OTHER_METHODS:
        if col in photos:
            means = setting_means(photos, col)
            method_fit = stats.linregress(omega, means.dphi.to_numpy())
            methods[name] = dict(means=means, a=method_fit.slope, u_a=method_fit.stderr,
                                 rms=rms_from_line(omega, means.dphi.to_numpy(), method_fit),
                                 r2=method_fit.rvalue ** 2, col=col)
    slopes = {'ring minima': a, **{name: m['a'] for name, m in methods.items()}}
    u_method = (max(slopes.values()) - min(slopes.values())) / 2

    # 6. The calibration part: every moving omega is predicted by the same line from the same table, so
    #    the relative prediction uncertainty (k's relative uncertainty and the scatter of the timed
    #    settings about the line) is a common scale error of omega, hence of a. (omega0 shifts the
    #    moving points together and is absorbed by the intercept, up to the rest point.)
    u_k_rel = cal_line['u_k'] / cal_line['k']
    rel_cal = float(np.hypot(u_k_rel, cal_line['rel_rms']))
    u_cal = a * rel_cal

    # 7. The result and the comparison with Eq. (5).
    u_a = float(np.sqrt(u_stat ** 2 + u_method ** 2 + u_cal ** 2))
    a_th, th_parts = theory()
    u_th = float(np.sqrt(sum(v ** 2 for v in th_parts.values())))
    z = (a - a_th) / np.hypot(u_a, u_th)
    p_z = float(2 * stats.norm.sf(abs(z)))
    z_stat = (a - a_th) / np.hypot(u_stat, u_th)
    ratio, u_ratio = a / a_th, a / a_th * np.hypot(u_a / a, u_th / a_th)
    N_implied = N_TURNS * ratio
    u_N_implied = N_implied * np.hypot(u_a / a, 2 * U_RADIUS / RADIUS)
    dphi_max, u_dphi_max = float(dphi[-1]), float(np.sqrt(2) * rms)      # top setting minus rest, each known to rms
    dt = dphi_max * WAVELENGTH / (2 * np.pi * C_LIGHT)
    dt_pred = a_th * omega[-1] * WAVELENGTH / (2 * np.pi * C_LIGHT)

    # 8. The speed of light.
    a_parts = {'regression': u_stat, 'method': u_method, 'calibration': u_cal}
    c_meas, c_parts = speed_of_light(a, a_parts)
    u_c = float(np.sqrt(sum(v ** 2 for v in c_parts.values())))
    z_c = (c_meas - C_LIGHT) / u_c
    p_c = float(2 * stats.norm.sf(abs(z_c)))

    # 9. Drift: a phase linear in time fitted together with the slope (settings taken in increasing order).
    t_mean = (photos.groupby('S').time_s.mean() - photos.time_s.min()).to_numpy()
    (a_drift, D, b_drift), drift_cov = optimize.curve_fit(drift_model, np.vstack([omega, t_mean]), dphi)
    u_a_drift, u_D = np.sqrt(drift_cov[0, 0]), np.sqrt(drift_cov[1, 1])
    corr_aD = float(drift_cov[0, 1] / (u_a_drift * u_D))
    p_D = float(2 * stats.t.sf(abs(D / u_D), n_settings - 3))
    random_pixels = random_pixel_summary()

    # 10. Outputs. The per-setting table first.
    settings['omega'] = omega
    settings['u_omega'] = u_omega
    settings['resid'] = residual
    for name, m in methods.items():
        key = m['col'].removeprefix('phi_').removesuffix('_rad')
        settings[f'dphi_{key}'] = m['means'].dphi
        settings[f'sem_{key}'] = m['means']['sem']
    settings.to_csv(HERE / 'sagnac_speeds.csv', index=False, float_format='%.6g')

    E8 = 1e8
    lines = ['Sagnac phase against angular velocity (sagnac.py), fibre loop centred on the rotation axis\n',
             'PER SETTING: phase at the ring centre (minima), mean over the photos, SEM = SD/sqrt(n); omega from the',
             'calibration line +- its prediction uncertainty (common to all settings); dphi = phase minus the rest mean.',
             f"{'S':>6} {'n':>3} {'first photo':>11} {'omega (rad/s)':>18} {'dphi (rad)':>18} {'P (10^5 px^2)':>13} {'resid':>7}"]
    for row in settings.itertuples():
        lines.append(f'{row.S:6.2f} {row.n:3d} {hhmmss(row.t_first):>11} {row.omega:8.4f} +- {row.u_omega:.4f} '
                     f'{row.dphi:8.3f} +- {row.sem:.3f} {row.P / 1e5:13.3f} {row.resid:+7.3f}')
    lines.append(f'(the rest mean, SEM {sem[0]:.3f} rad, is subtracted from every row; resid = dphi minus the line)\n')
    lines.append('THE LINE dphi = a omega + b, linregress through the ten setting means (Schwartz Eq. 7):')
    lines.append(f'  a = {a:.4f} +- {u_stat:.4f} rad/(rad/s), b = {b:.3f} +- {fit.intercept_stderr:.3f} rad '
                 f'(p = {p_b:.2f} for b = 0), R^2 = {r2:.4f}')
    lines.append(f'  s = {rms:.4f} rad: the rms distance of the points from the line, which sets the slope\'s uncertainty.')
    lines.append(f'  The line against the photos\' SEM: chi2 = {chi2_sem:.1f}, dof = {dof}, p = {p_sem:.1g}: the settings lie '
                 f'{np.sqrt(chi2_sem / dof):.1f} times farther')
    lines.append('  from the line than the photo noise allows (the phase wanders between settings), so s is used, not the SEM.')
    lines.append(f'  Added omega^2 term: {quad[0]:+.3f} +- {np.sqrt(quad_cov[0, 0]):.3f} rad/(rad/s)^2, t = {t_quad:.2f}, '
                 f'p = {p_quad:.2f}: no bend.\n')
    lines.append('OTHER PHASE-READING METHODS, same averaging and regression (cosine fit: fringes.py, prep Q5; template')
    lines.append('matching: methods.py):')
    for name, m in methods.items():
        lines.append(f'  {name:18s} a = {m["a"]:.4f} +- {m["u_a"]:.4f}, s = {m["rms"]:.4f} rad, R^2 = {m["r2"]:.4f}; '
                     'minus the minima per setting: ' + ', '.join(f'{d:+.2f}' for d in m['means'].dphi - dphi) + ' rad')
    lines.append(f'  Method part of u(a): half the range of {", ".join(f"{v:.4f}" for v in slopes.values())} = {u_method:.4f} '
                 f'({100 * u_method / a:.1f} %).')
    lines.append(f'  Calibration part of u(a): k is known to {100 * u_k_rel:.2f} % and the timed settings scatter about the '
                 f'line by {100 * cal_line["rel_rms"]:.2f} % (calibration.txt);')
    lines.append(f'  both are common to every omega, so u(a)/a = sqrt({100 * u_k_rel:.2f}^2 + {100 * cal_line["rel_rms"]:.2f}^2) % '
                 f'= {100 * rel_cal:.2f} %, u = {u_cal:.4f}.\n')
    lines.append(f'RESULT  a = {a:.4f} +- {u_stat:.4f} (regression) +- {u_method:.4f} (method) +- {u_cal:.4f} (calibration)')
    lines.append(f'          = {a:.4f} +- {u_a:.4f} rad/(rad/s) ({100 * u_a / a:.1f} %, the three parts in quadrature)')
    lines.append(f"Theory, Eq. (5): N = {N_TURNS} +- {U_N_TURNS}, R = {100 * RADIUS:.1f} +- {100 * U_RADIUS:.1f} cm "
                 f"(A = pi R^2 = {np.pi * RADIUS ** 2:.4f} m^2), lambda = {WAVELENGTH * 1e9:.0f} nm, c exact.")
    lines.append(f'  a_th = 8 pi^2 N R^2 / (lambda c) = {a_th:.4f} +- {u_th:.4f} rad/(rad/s) '
                 f'[N: {th_parts["N"]:.4f} ({100 * U_N_TURNS / N_TURNS:.1f} %), R: {th_parts["R"]:.4f} ({200 * U_RADIUS / RADIUS:.1f} %)]')
    lines.append('  lambda: no uncertainty given, NOT included (unknown, not zero); a_th and c move by '
                 f'{100e-9 / WAVELENGTH:.2f} % per nm.')
    lines.append(f'  a / a_th = {ratio:.3f} +- {u_ratio:.3f}; a - a_th = {a - a_th:+.4f} = {z:+.2f} combined sigma, p = {p_z:.2f} '
                 f'(regression part only: {z_stat:+.2f} sigma, p = {2 * stats.norm.sf(abs(z_stat)):.2f}).')
    lines.append(f'  Read the other way (R, lambda, c known): N = {N_implied:.2f} +- {u_N_implied:.2f} turns.')
    lines.append(f'  Largest shift (S = 90 minus rest): {dphi_max:.3f} +- {u_dphi_max:.3f} rad = {dphi_max / (2 * np.pi):.3f} +- '
                 f'{u_dphi_max / (2 * np.pi):.3f} fringe at omega = {omega[-1]:.3f} rad/s; predicted {a_th * omega[-1]:.2f} rad.')
    lines.append(f'  As a travel-time difference, dt = dphi lambda / (2 pi c) = {dt * 1e15:.2f} +- '
                 f'{dt * 1e15 * u_dphi_max / dphi_max:.2f} fs (predicted {dt_pred * 1e15:.2f} fs; one optical period is '
                 f'{WAVELENGTH / C_LIGHT * 1e15:.2f} fs).\n')
    lines.append('SPEED OF LIGHT, Eq. (5) solved for c = 8 pi^2 N R^2 / (lambda a), in 10^8 m/s:')
    lines.append(f'  c = {c_meas / E8:.4f}' + ''.join(f' +- {v / E8:.4f} ({name}, {100 * v / c_meas:.1f} %)'
                                                      for name, v in c_parts.items()))
    lines.append(f'    = {c_meas / E8:.4f} +- {u_c / E8:.4f} ({100 * u_c / c_meas:.1f} %); c / c_def = {c_meas / C_LIGHT:.3f} +- '
                 f'{u_c / C_LIGHT:.3f}, {z_c:+.2f} sigma, p = {p_c:.2f}.')
    lines.append('  This is a / a_th read the other way (c / c_def = a_th / a), not an independent test.\n')
    lines.append('DRIFT CHECK. The settings were photographed in increasing order with rest photos only at the start, so')
    lines.append('omega and time rise together. dphi = a omega + D t + b (t: mean photo time of the setting):')
    lines.append(f'  a = {a_drift:.2f} +- {u_a_drift:.2f} rad/(rad/s), D = {1e3 * D:+.2f} +- {1e3 * u_D:.2f} mrad/s '
                 f'(p = {p_D:.2f}), correlation(a, D) = {corr_aD:.3f}: a steady drift cannot be told from the slope.')
    lines.append('  Not in the quoted uncertainty. Also not checked: the sign of the shift (rotation direction not recorded).\n')
    if random_pixels:
        lines.append('RANDOM PIXELS (prep Q5, optional; presentation stage, random_points_scan.json):')
        for n_pixels, v in random_pixels.items():
            lines.append(f'  {n_pixels:3d} pixels per draw: {v["draws"]} draws, slope mean {v["mean"]:.3f}, SD between draws '
                         f'{v["sd"]:.3f}, fit error of one draw {v["u_fit"]:.3f}'
                         + (f', {v["failed"]} failed' if v['failed'] else ''))
        lines.append('')
    lines.append('REPORT FORM (uncertainty to one significant figure, value to the same digit; comparisons above use')
    lines.append('unrounded values):')
    lines.append(f'  k = {round1(cal_line["k"], cal_line["u_k"])} rad/s per unit S, '
                 f'omega0 = {round1(cal_line["omega0"], cal_line["u_omega0"])} rad/s')
    lines.append(f'  a = {round1(a, u_a)} rad/(rad/s) (all parts); {round1(a, u_stat)} (regression only)')
    lines.append(f'  a_th = {round1(a_th, u_th)} rad/(rad/s); a / a_th = {round1(ratio, u_ratio)}; '
                 f'N implied = {round1(N_implied, u_N_implied)}')
    lines.append('  ' + '; '.join(f'{name}: a = {round1(m["a"], m["u_a"])}' for name, m in methods.items()))
    lines.append(f'  largest shift = {round1(dphi_max, u_dphi_max)} rad = '
                 f'{round1(dphi_max / (2 * np.pi), u_dphi_max / (2 * np.pi))} fringe, '
                 f'dt = {round1(dt * 1e15, dt * 1e15 * u_dphi_max / dphi_max)} fs')
    lines.append(f'  c = ({round1(c_meas / E8, u_c / E8)}) x 10^8 m/s; c / c_def = {round1(c_meas / C_LIGHT, u_c / C_LIGHT)}')
    (HERE / 'sagnac.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print('\n'.join(lines))

    results = dict(
        slope=dict(a=float(a), u_stat=float(u_stat), u_method=float(u_method), u_calibration=float(u_cal), u=u_a,
                   b=float(b), u_b=float(fit.intercept_stderr), p_b=p_b, s=rms, r2=float(r2), dof=dof,
                   chi2_vs_sem=chi2_sem, p_vs_sem=p_sem, omega2_term=float(quad[0]), u_omega2_term=float(np.sqrt(quad_cov[0, 0])),
                   t_omega2=t_quad, p_omega2=p_quad, u_k_rel=float(u_k_rel), rel_calibration=rel_cal),
        methods={name: dict(a=float(m['a']), u_a=float(m['u_a']), s=m['rms'], r2=float(m['r2'])) for name, m in methods.items()},
        theory=dict(a_th=a_th, u_a_th=u_th, components=th_parts, u_wavelength='not given', per_nm=1e-9 / WAVELENGTH),
        comparison=dict(ratio=ratio, u_ratio=u_ratio, z=float(z), p=p_z, z_regression_only=float(z_stat),
                        N_implied=N_implied, u_N_implied=u_N_implied),
        largest_shift=dict(dphi=dphi_max, u=u_dphi_max, fringes=dphi_max / (2 * np.pi), dt_fs=dt * 1e15,
                           u_dt_fs=dt * 1e15 * u_dphi_max / dphi_max, dt_predicted_fs=dt_pred * 1e15,
                           omega=float(omega[-1]), dphi_predicted=a_th * float(omega[-1])),
        speed_of_light=dict(c=c_meas, u_c=u_c, components=c_parts, rel=100 * u_c / c_meas, ratio_to_defined=c_meas / C_LIGHT,
                            u_ratio_to_defined=u_c / C_LIGHT, z=float(z_c), p=p_c),
        drift=dict(a=float(a_drift), u_a=float(u_a_drift), D=float(D), u_D=float(u_D), p_D=p_D, correlation=corr_aD),
        random_pixels=random_pixels,
        report_form=lines[lines.index('unrounded values):') + 1:])
    (HERE / 'sagnac.json').write_text(json.dumps(results, indent=1, default=float))


if __name__ == '__main__':
    main()
