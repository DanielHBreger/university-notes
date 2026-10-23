"""Rotation calibration (task 1): the table's angular velocity omega at each controller setting S.

The Sagnac phase is proportional to omega, but during the photographs only the controller
setting S is known. This script turns the stopwatch data in "Rotation calibration.csv" into the
straight line omega = k S + omega0 that sagnac.py uses to convert the photo settings, with the
uncertainty of omega at every photo setting.

Steps (main(), top to bottom)
  1. load_turns     per setting S, the times of 5 or 10 consecutive single turns.
  2. per_setting    mean turn time T with its uncertainty (lap scatter and reaction time),
                    omega = 2 pi / T.
  3. curve_fit      the line, weighted by the timing uncertainties.
  4. save           calibration.json (read by sagnac.py and report_figures.py), calibration.txt.

Uncertainty of a timed setting. The n turns are consecutive laps, so the mean turn time is
(t_n - t_0)/n and only the first and the last button press enter it. With a human reaction time
sigma_r per press that part is sqrt(2) sigma_r / n. The scatter of the single laps gives the
statistical part SD/sqrt(n). The two are added in quadrature:
    u(T)^2 = (SD/sqrt(n))^2 + (sqrt(2) sigma_r / n)^2,     u(omega)/omega = u(T)/T.
sigma_r = REACTION_TIME_S (0.2 s, a conventional value; the lap scatter itself implies a smaller
press-to-press jitter, so this is on the safe side).

The line. scipy's curve_fit with the timing uncertainties as sigma gives k and omega0 with
uncertainties propagated from those (absolute_sigma=True, no rescaling). The line does not
describe the points to within their uncertainties (chi2/dof about 4: the table turns slower than
the line at high settings). That misfit is carried as the rms scatter of the timed settings about
the line: omega read from the line at a setting that was not timed gets, besides the line's own
uncertainty, this scatter (u_omega_pred), and sagnac.py takes the calibration part of the Sagnac
slope's uncertainty as the same relative size.

Run from this folder:  python calibration.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

HERE = Path(__file__).resolve().parent
CSV = HERE / 'Rotation calibration.csv'
REACTION_TIME_S = 0.2          # human reaction time per stopwatch press (s)
PHOTO_SETTINGS = [12.59, 20.74, 31.11, 44.07, 55.83, 64.00, 73.65, 82.18, 90.0]


def load_turns(path=CSV):
    """Step 1. Column 'Speed' is the setting S, columns 'turn 1'...'turn 10' are the times (s) of
    single turns; empty cells are skipped. Returns S (array) and a list of turn-time arrays."""
    table = pd.read_csv(path, encoding='utf-8-sig')         # utf-8-sig: the file starts with a BOM
    S = table['Speed'].to_numpy(float)
    turns = [row[~np.isnan(row)] for row in table.drop(columns='Speed').to_numpy(float)]
    return S, turns


def per_setting(S, turns, sigma_r=REACTION_TIME_S):
    """Step 2. One row per setting: n turns, mean turn time T, SD of the turns, the statistical
    part SD/sqrt(n), the reaction-time part sqrt(2) sigma_r / n, their quadrature sum u_T,
    omega = 2 pi / T and u_omega = omega u_T / T."""
    rows = []
    for setting, times in zip(S, turns):
        n = len(times)
        T, sd = times.mean(), times.std(ddof=1)
        u_T_stat, u_T_reaction = sd / np.sqrt(n), np.sqrt(2) * sigma_r / n
        u_T = np.hypot(u_T_stat, u_T_reaction)
        omega = 2 * np.pi / T
        rows.append(dict(S=setting, n=n, T=T, sd_T=sd, u_T_stat=u_T_stat, u_T_reaction=u_T_reaction, u_T=u_T,
                         omega=omega, u_omega=omega * u_T / T))
    return pd.DataFrame(rows)


def line(S, omega0, k):
    return omega0 + k * S


def omega_of(S, fit):
    """The line at the settings S. fit: the 'line' entry of calibration.json."""
    return line(np.asarray(S, float), fit['omega0'], fit['k'])


def u_omega_pred(S, fit):
    """Uncertainty of omega predicted by the line at a setting that was not timed: the line's own
    uncertainty (from the covariance of omega0 and k) and the rms scatter of the timed settings
    about the line, in quadrature."""
    S = np.asarray(S, float)
    cov = np.asarray(fit['cov'])
    u_line = np.sqrt(cov[0, 0] + 2 * S * cov[0, 1] + S ** 2 * cov[1, 1])
    return np.hypot(u_line, fit['rel_rms'] * omega_of(S, fit))


def main():
    S, turns = load_turns()
    settings = per_setting(S, turns)
    omega, u_omega = settings.omega.to_numpy(), settings.u_omega.to_numpy()

    (omega0, k), cov = optimize.curve_fit(line, S, omega, sigma=u_omega, absolute_sigma=True)
    residual = omega - line(S, omega0, k)
    chi2, dof = float(np.sum((residual / u_omega) ** 2)), len(S) - 2
    fit = dict(k=float(k), u_k=float(np.sqrt(cov[1, 1])), omega0=float(omega0), u_omega0=float(np.sqrt(cov[0, 0])),
               cov=cov.tolist(), chi2=chi2, dof=dof, p=float(stats.chi2.sf(chi2, dof)),
               rel_rms=float(np.sqrt(np.mean((residual / omega) ** 2))), resid_rel=(residual / omega).tolist())
    fit['omega90'], fit['u_omega90'] = float(omega_of(90, fit)), float(u_omega_pred(90, fit))
    results = dict(reaction_time_s=REACTION_TIME_S, settings=settings.to_dict('records'), line=fit,
                   omega90_timed=float(omega[-1]),
                   note='omega in rad/s, k in rad/s per unit S; cov is the covariance of (omega0, k) from curve_fit '
                        'with the timing uncertainties; rel_rms = rms of (omega - line)/omega over the timed settings')
    (HERE / 'calibration.json').write_text(json.dumps(results, indent=1))

    lines = ['Rotation calibration (calibration.py): omega = 2 pi / T at each controller setting S\n',
             f'Uncertainty of the mean turn time: u(T)^2 = (SD/sqrt(n))^2 + (sqrt(2) sigma_r / n)^2, with the human '
             f'reaction time sigma_r = {REACTION_TIME_S} s per press (two presses set the mean of n consecutive laps).\n',
             f"{'S':>6} {'n':>2} {'T (s)':>8} {'SD (s)':>7} {'SD/sqrt(n)':>10} {'reaction':>8} {'u_T (s)':>8} "
             f"{'omega (rad/s)':>13} {'u_omega':>8} {'rel':>6}"]
    for row in settings.itertuples():
        lines.append(f'{row.S:6.2f} {row.n:2d} {row.T:8.3f} {row.sd_T:7.3f} {row.u_T_stat:10.4f} {row.u_T_reaction:8.4f} '
                     f'{row.u_T:8.4f} {row.omega:13.5f} {row.u_omega:8.5f} {100 * row.u_omega / row.omega:5.2f}%')
    lines.append('\nLine omega = k S + omega0, curve_fit weighted by the timing uncertainties:')
    lines.append(f"  k = {k:.5f} +- {fit['u_k']:.5f} rad/s per unit S ({100 * fit['u_k'] / k:.2f} %), "
                 f"omega0 = {omega0:.4f} +- {fit['u_omega0']:.4f} rad/s")
    lines.append(f"  chi2 = {chi2:.1f}, dof = {dof}, chi2/dof = {chi2 / dof:.2f}, p = {fit['p']:.1g}: the points lie "
                 'farther from the line than their uncertainties allow.')
    lines.append('  (omega - line)/omega: ' + ', '.join(f'{100 * v:+.1f}' for v in residual / omega)
                 + f" %  (rms {100 * fit['rel_rms']:.2f} %); the table turns slower than the line at high settings.")
    lines.append(f"  At S = 90: timed {omega[-1]:.4f} rad/s, line {fit['omega90']:.4f} +- {fit['u_omega90']:.4f} "
                 f"({100 * (fit['omega90'] / omega[-1] - 1):+.1f} %).\n")
    lines.append('omega at the photo settings, predicted by the line: u^2 = (line uncertainty)^2 + (rms scatter)^2, '
                 f"the scatter being {100 * fit['rel_rms']:.2f} % of omega; both parts are common to all settings:")
    for setting in PHOTO_SETTINGS:
        omega_line, u_pred = float(omega_of(setting, fit)), float(u_omega_pred(setting, fit))
        lines.append(f'  S = {setting:6.2f}: omega = {omega_line:.4f} +- {u_pred:.4f} rad/s ({100 * u_pred / omega_line:.2f} %)')
    (HERE / 'calibration.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
