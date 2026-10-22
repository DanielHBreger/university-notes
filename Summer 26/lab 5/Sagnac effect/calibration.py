"""Rotation calibration (task 1): the table's angular velocity omega at each controller setting S.

The Sagnac phase is proportional to omega, but during the photographs only the controller
setting S is known. This script turns the stopwatch data in "Rotation calibration.csv" into the
straight line omega = k S + omega0 that sagnac.py uses to convert the photo settings, and
reports how well that line describes the timed settings.

Steps (main(), top to bottom)
  1. load_turns     per setting S, the times of 5 or 8 consecutive single turns.
  2. per_setting    mean turn time T, its standard error SD/sqrt(n), omega = 2 pi / T.
  3. fit            the line (used, task 1) and a quadratic (check), both by least squares with
                    each point weighted by 1/omega^2, i.e. equal relative weight per setting.
  4. save           calibration.json (read by sagnac.py and report_figures.py), calibration.txt.

Why relative weights. The timed points scatter about any smooth curve by about 1 % of omega,
far more than the stopwatch uncertainty (0.1-0.5 % of omega), and that scatter is roughly a
constant fraction of omega. An unweighted fit would be dominated by the fast settings and miss
the slowest one by 13 %. Minimising the relative residuals treats every setting alike. The
coefficient uncertainties are taken from the residual scatter (the usual regression estimate),
so they describe how well the points fix the line, not the stopwatch.

What the result means. The line is rejected by the stopwatch errors (chi2 >> dof): the table
does not run exactly linearly with S, it turns slower than the line at high settings (2.3 %
below it at S = 90). The photos use the line anyway, as task 1 asks; the quadratic, which fits
to 0.5 %, is the check that sets the calibration uncertainty of the Sagnac slope (sagnac.py).

Run from this folder:  python calibration.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
CSV = HERE / 'Rotation calibration.csv'
PHOTO_SETTINGS = [12.59, 20.74, 31.11, 44.07, 55.83, 64.00, 73.65, 82.18, 90.0]


def load_turns(path=CSV):
    """Step 1. Column 'Speed' is the setting S, columns 'turn 1'...'turn 10' are the times (s) of
    single turns; empty cells are skipped. Returns S (array) and a list of turn-time arrays."""
    df = pd.read_csv(path, encoding='utf-8-sig')       # utf-8-sig: the file starts with a BOM
    S = df['Speed'].to_numpy(float)
    turns = [row[~np.isnan(row)] for row in df.drop(columns='Speed').to_numpy(float)]
    return S, turns


def per_setting(S, turns):
    """Step 2. One row per setting: n turns, mean turn time T, SD of the turns, the standard
    error u_T = SD/sqrt(n), omega = 2 pi / T and u_omega = omega u_T / T."""
    rows = []
    for s, t in zip(S, turns):
        n = len(t)
        T, sd = t.mean(), t.std(ddof=1)
        w = 2 * np.pi / T
        rows.append(dict(S=s, n=n, T=T, sd_T=sd, u_T=sd / np.sqrt(n), omega=w, u_omega=w * sd / np.sqrt(n) / T))
    return pd.DataFrame(rows)


def fit_poly(x, y, deg, weights):
    """Weighted least squares of y = c0 + c1 x + ... + c_deg x^deg with weights W (one per point).
    The coefficient covariance is scaled by the residual scatter, chi2_w/dof, so it reflects how
    far the points lie from the curve. Returns the coefficients (constant first), their covariance,
    the residuals and the rms of the relative residuals (y - curve)/y."""
    A = np.vander(x, deg + 1, increasing=True)          # columns 1, x, x^2, ...
    cov0 = np.linalg.inv(A.T @ (A * weights[:, None]))
    c = cov0 @ (A.T @ (weights * y))
    r = y - A @ c
    dof = len(x) - deg - 1
    cov = cov0 * float(np.sum(weights * r ** 2)) / dof
    return c, cov, r, float(np.sqrt(np.mean((r / y) ** 2)))


def omega_of(S, coef):
    """The fitted polynomial (coefficients constant first) at the settings S."""
    return np.vander(np.atleast_1d(np.asarray(S, float)), len(coef), increasing=True) @ np.asarray(coef)


def main():
    S, turns = load_turns()
    tab = per_setting(S, turns)
    x, y, uy = tab.S.to_numpy(), tab.omega.to_numpy(), tab.u_omega.to_numpy()
    weights = 1 / y ** 2                                 # equal relative weight per setting

    # 3. The line (task 1) and the quadratic (check).
    c1, cov1, r1, rms1 = fit_poly(x, y, 1, weights)
    c2, cov2, r2, rms2 = fit_poly(x, y, 2, weights)
    chi2_sw = float(np.sum((r1 / uy) ** 2))              # the line against the stopwatch errors alone
    dof1 = len(x) - 2
    k, u_k = c1[1], float(np.sqrt(cov1[1, 1]))
    w0, u_w0 = c1[0], float(np.sqrt(cov1[0, 0]))

    # 4. Save.
    out = dict(
        settings=tab.to_dict('records'),
        line=dict(coef=c1.tolist(), cov=cov1.tolist(), k=k, u_k=u_k, omega0=w0, u_omega0=u_w0,
                  rel_rms=rms1, resid_rel=(r1 / y).tolist(), chi2_stopwatch=chi2_sw, dof=dof1,
                  p_stopwatch=float(stats.chi2.sf(chi2_sw, dof1)), omega90=float(omega_of(90, c1)[0])),
        quadratic=dict(coef=c2.tolist(), rel_rms=rms2, resid_rel=(r2 / y).tolist(),
                       omega90=float(omega_of(90, c2)[0])),
        omega90_timed=float(y[-1]),
        note='omega in rad/s; coef constant first (c0 in rad/s, c1 in rad/s per unit S, ...); both fits '
             'weighted by 1/omega^2; cov scaled by the residual scatter; rel_rms = rms of (omega - fit)/omega')
    (HERE / 'calibration.json').write_text(json.dumps(out, indent=1))

    L = ['Rotation calibration (calibration.py): omega = 2 pi / T at each controller setting S\n',
         f"{'S':>6} {'n':>2} {'T (s)':>8} {'SD (s)':>7} {'u_T (s)':>8} {'omega (rad/s)':>13} {'u_omega':>8} {'rel':>6}"]
    for r in tab.itertuples():
        L.append(f'{r.S:6.2f} {r.n:2d} {r.T:8.3f} {r.sd_T:7.3f} {r.u_T:8.4f} {r.omega:13.5f} {r.u_omega:8.5f} '
                 f'{100 * r.u_omega / r.omega:5.2f}%')
    L.append('(u_T = SD/sqrt(n), the standard error of the mean turn time; u_omega = omega u_T / T)\n')
    L.append('Fits of omega(S), least squares with weights 1/omega^2 (equal relative weight per setting):')
    L.append(f'  line      omega = k S + omega0:  k = {k:.5f} +- {u_k:.5f} rad/s per unit S ({100 * u_k / k:.2f} %), '
             f'omega0 = {w0:.4f} +- {u_w0:.4f} rad/s')
    L.append('            residuals (omega - line)/omega: ' + ', '.join(f'{100 * v:+.1f}' for v in r1 / y)
             + f' %  (rms {100 * rms1:.2f} %)')
    L.append(f'            against the stopwatch errors alone: chi2 = {chi2_sw:.0f}, dof = {dof1}, '
             f'p = {stats.chi2.sf(chi2_sw, dof1):.1g}: the table does not follow the line to the stopwatch precision.')
    L.append(f'  quadratic c0 = {c2[0]:.4f}, c1 = {c2[1]:.5f}, c2 = {c2[2]:.3g}: residuals '
             + ', '.join(f'{100 * v:+.1f}' for v in r2 / y) + f' %  (rms {100 * rms2:.2f} %)')
    L.append(f'\nAt S = 90: timed {y[-1]:.4f} rad/s, line {out["line"]["omega90"]:.4f} '
             f'({100 * (out["line"]["omega90"] / y[-1] - 1):+.1f} %), quadratic {out["quadratic"]["omega90"]:.4f} '
             f'({100 * (out["quadratic"]["omega90"] / y[-1] - 1):+.1f} %).')
    L.append('The photos take omega from the line (task 1). Its uncertainty at a setting is taken as the rms '
             f'scatter of the timed points about the line, {100 * rms1:.1f} % of omega. The quadratic is the check '
             'that sets the calibration uncertainty of the Sagnac slope (sagnac.py).\n')
    L.append('omega at the photo settings (line, +- rms scatter; quadratic for comparison):')
    for s in PHOTO_SETTINGS:
        wl, wq = omega_of(s, c1)[0], omega_of(s, c2)[0]
        L.append(f'  S = {s:6.2f}: line {wl:.4f} +- {rms1 * wl:.4f} rad/s, quadratic {wq:.4f} ({100 * (wq / wl - 1):+.1f} %)')
    (HERE / 'calibration.txt').write_text('\n'.join(L) + '\n', encoding='utf-8')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
