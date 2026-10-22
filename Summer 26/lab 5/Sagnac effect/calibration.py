"""Rotation calibration (task 1): the table's angular velocity omega at each controller setting S.

The Sagnac phase is proportional to omega, but during the photographs only the controller
setting S is known. This script turns the stopwatch data in "Rotation calibration.csv" into the
straight line omega = k S + omega0 that sagnac.py uses to convert the photo settings, with the
uncertainty of the line and of omega at every photo setting.

Steps (main(), top to bottom)
  1. load_turns     per setting S, the times of 5 or 8 consecutive single turns.
  2. per_setting    mean turn time T with its uncertainty (lap scatter and reaction time),
                    omega = 2 pi / T.
  3. fit            the line by weighted least squares, weights 1/u(omega)^2; its coefficient
                    covariance is enlarged by chi2/dof when the points scatter more than their
                    uncertainties allow.
  4. save           calibration.json (read by sagnac.py and report_figures.py), calibration.txt.

Uncertainty of a timed setting. The n turns are consecutive laps, so the mean turn time is
(t_n - t_0)/n and only the first and the last button press enter it. With a human reaction time
sigma_r per press, that part is sqrt(2) sigma_r / n. The scatter of the single laps gives the
statistical part SD/sqrt(n). The two are added in quadrature:
    u(T)^2 = (SD/sqrt(n))^2 + (sqrt(2) sigma_r / n)^2,     u(omega)/omega = u(T)/T.
sigma_r = REACTION_TIME_S (0.2 s, a conventional value; the lap scatter itself implies a smaller
press-to-press jitter, so this is on the safe side).

The line. Each point weighted by 1/u(omega)^2. The line does not describe the points to within
these uncertainties (chi2/dof > 1: the table is not exactly linear in S, it turns slower than the
line at high settings), so the coefficient covariance is multiplied by chi2/dof, the usual
remedy, and u(k) then reflects the actual scatter.

omega at a photo setting. A photo setting was not timed: the line predicts its omega. The
prediction carries the line's own uncertainty, (1, S) Cov (1, S)^T, and the scatter of the timed
settings about the line, rel_rms * omega (the table does not repeat a setting better than that):
    u(omega)^2 = u_line^2 + (rel_rms omega)^2        (u_omega_pred)
Both parts are common to all photo settings (the same line, the same table), so sagnac.py takes
the calibration part of the Sagnac slope's uncertainty as the same relative size,
sqrt((u_k/k)^2 + rel_rms^2).

Run from this folder:  python calibration.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
CSV = HERE / 'Rotation calibration.csv'
REACTION_TIME_S = 0.2          # human reaction time per stopwatch press (s)
PHOTO_SETTINGS = [12.59, 20.74, 31.11, 44.07, 55.83, 64.00, 73.65, 82.18, 90.0]


def load_turns(path=CSV):
    """Step 1. Column 'Speed' is the setting S, columns 'turn 1'...'turn 10' are the times (s) of
    single turns; empty cells are skipped. Returns S (array) and a list of turn-time arrays."""
    df = pd.read_csv(path, encoding='utf-8-sig')       # utf-8-sig: the file starts with a BOM
    S = df['Speed'].to_numpy(float)
    turns = [row[~np.isnan(row)] for row in df.drop(columns='Speed').to_numpy(float)]
    return S, turns


def per_setting(S, turns, sigma_r=REACTION_TIME_S):
    """Step 2. One row per setting: n turns, mean turn time T, SD of the turns, the statistical
    part SD/sqrt(n), the reaction-time part sqrt(2) sigma_r / n, their quadrature sum u_T,
    omega = 2 pi / T and u_omega = omega u_T / T."""
    rows = []
    for s, t in zip(S, turns):
        n = len(t)
        T, sd = t.mean(), t.std(ddof=1)
        u_stat, u_react = sd / np.sqrt(n), np.sqrt(2) * sigma_r / n
        uT = np.hypot(u_stat, u_react)
        w = 2 * np.pi / T
        rows.append(dict(S=s, n=n, T=T, sd_T=sd, u_T_stat=u_stat, u_T_reaction=u_react, u_T=uT,
                         omega=w, u_omega=w * uT / T))
    return pd.DataFrame(rows)


def fit_line(x, y, u):
    """Step 3. Weighted least squares y = c0 + c1 x with weights 1/u^2 (the textbook normal
    equations with the design matrix A = [1, x]). The covariance (A^T W A)^-1 is multiplied by
    chi2/dof when chi2/dof > 1, so that the coefficient uncertainties reflect the scatter of the
    points about the line and not only the error bars. Returns the coefficients (constant first),
    their covariance, chi2, dof, p and the scale factor applied."""
    A = np.c_[np.ones_like(x), x]
    W = 1 / u ** 2
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    c = cov @ (A.T @ (W * y))
    chi2 = float(np.sum(((y - A @ c) / u) ** 2))
    dof = len(x) - 2
    scale = max(1.0, chi2 / dof)
    return c, cov * scale, chi2, dof, float(stats.chi2.sf(chi2, dof)), scale


def omega_of(S, coef):
    """The fitted polynomial (coefficients constant first) at the settings S."""
    return np.vander(np.atleast_1d(np.asarray(S, float)), len(coef), increasing=True) @ np.asarray(coef)


def u_omega_of(S, cov):
    """Uncertainty of the line at the settings S: sqrt((1, S) Cov (1, S)^T), correlated between
    settings because they share the same line."""
    A = np.vander(np.atleast_1d(np.asarray(S, float)), 2, increasing=True)
    return np.sqrt(np.einsum('ij,jk,ik->i', A, np.asarray(cov), A))


def u_omega_pred(S, line):
    """Uncertainty of omega predicted by the line at a setting that was not timed: the line's
    uncertainty and the rms scatter of the timed settings about the line, in quadrature.
    line: the 'line' entry of calibration.json (coef, cov, rel_rms)."""
    w = omega_of(S, line['coef'])
    return np.hypot(u_omega_of(S, line['cov']), line['rel_rms'] * w)


def main():
    S, turns = load_turns()
    tab = per_setting(S, turns)
    x, y, uy = tab.S.to_numpy(), tab.omega.to_numpy(), tab.u_omega.to_numpy()

    c, cov, chi2, dof, p, scale = fit_line(x, y, uy)
    k, u_k = c[1], float(np.sqrt(cov[1, 1]))
    w0, u_w0 = c[0], float(np.sqrt(cov[0, 0]))
    resid = y - omega_of(x, c)
    rms = float(np.sqrt(np.mean((resid / y) ** 2)))

    out = dict(
        reaction_time_s=REACTION_TIME_S,
        settings=tab.to_dict('records'),
        line=dict(coef=c.tolist(), cov=cov.tolist(), k=k, u_k=u_k, omega0=w0, u_omega0=u_w0,
                  chi2=chi2, dof=dof, p=p, cov_scale=scale, rel_rms=rms, resid_rel=(resid / y).tolist(),
                  omega90=float(omega_of(90, c)[0]), u_omega90=float(u_omega_of(90, cov)[0])),
        omega90_timed=float(y[-1]),
        note='omega in rad/s; coef constant first (c0 in rad/s, c1 in rad/s per unit S); the fit is weighted by '
             '1/u_omega^2 and cov is multiplied by cov_scale = max(1, chi2/dof); rel_rms = rms of (omega - line)/omega')
    (HERE / 'calibration.json').write_text(json.dumps(out, indent=1))

    L = ['Rotation calibration (calibration.py): omega = 2 pi / T at each controller setting S\n',
         f'Uncertainty of the mean turn time: u(T)^2 = (SD/sqrt(n))^2 + (sqrt(2) sigma_r / n)^2, with the human '
         f'reaction time sigma_r = {REACTION_TIME_S} s per press (two presses set the mean of n consecutive laps).\n',
         f"{'S':>6} {'n':>2} {'T (s)':>8} {'SD (s)':>7} {'SD/sqrt(n)':>10} {'reaction':>8} {'u_T (s)':>8} "
         f"{'omega (rad/s)':>13} {'u_omega':>8} {'rel':>6}"]
    for r in tab.itertuples():
        L.append(f'{r.S:6.2f} {r.n:2d} {r.T:8.3f} {r.sd_T:7.3f} {r.u_T_stat:10.4f} {r.u_T_reaction:8.4f} {r.u_T:8.4f} '
                 f'{r.omega:13.5f} {r.u_omega:8.5f} {100 * r.u_omega / r.omega:5.2f}%')
    L.append('\nLine omega = k S + omega0, weighted least squares (weights 1/u_omega^2):')
    L.append(f'  k = {k:.5f} +- {u_k:.5f} rad/s per unit S ({100 * u_k / k:.2f} %), omega0 = {w0:.4f} +- {u_w0:.4f} rad/s')
    L.append(f'  chi2 = {chi2:.1f}, dof = {dof}, chi2/dof = {chi2 / dof:.2f}, p = {p:.1g}'
             + (f': the points scatter more than their uncertainties allow, so the coefficient covariance is '
                f'multiplied by {scale:.2f}.' if scale > 1 else '.'))
    L.append('  residuals (omega - line)/omega: ' + ', '.join(f'{100 * v:+.1f}' for v in resid / y)
             + f' %  (rms {100 * rms:.2f} %); the table turns slower than the line at high settings.')
    L.append(f'  At S = 90: timed {y[-1]:.4f} rad/s, line {out["line"]["omega90"]:.4f} +- {out["line"]["u_omega90"]:.4f} '
             f'({100 * (out["line"]["omega90"] / y[-1] - 1):+.1f} %).\n')
    L.append('omega at the photo settings, predicted by the line: u^2 = (line uncertainty)^2 + (rms scatter)^2, '
             f'the scatter being {100 * rms:.2f} % of omega; both parts are common to all settings:')
    for s in PHOTO_SETTINGS:
        wl, ul, up = omega_of(s, c)[0], u_omega_of(s, cov)[0], u_omega_pred(s, out['line'])[0]
        L.append(f'  S = {s:6.2f}: omega = {wl:.4f} +- {up:.4f} rad/s ({100 * up / wl:.2f} %; line alone {ul:.4f})')
    (HERE / 'calibration.txt').write_text('\n'.join(L) + '\n', encoding='utf-8')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
