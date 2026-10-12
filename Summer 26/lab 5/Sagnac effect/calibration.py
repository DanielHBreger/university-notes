"""Rotation calibration: table angular velocity against the motor speed setting.

Minimum task 1 of the manual ("the (linear) ratio between the voltage on the motor and its
angular velocity"). Input: Rotation calibration.csv, one row per speed setting S with the
times of 5 or 8 consecutive turns.

The turns were timed as stopwatch splits, so neighbouring laps share an endpoint:
lap_i = t_i - t_(i-1), and the mean lap is (t_n - t_0)/n. With an independent timing error
sigma_r on every press, Var(mean lap) = 2 sigma_r^2 / n^2 and E[SD^2 of the laps] =
2 sigma_r^2 (n+1)/n, so the stopwatch uncertainty of the mean lap is
    u_sw(T) = SD / sqrt(n (n+1))
(not SD/sqrt(n), which treats the laps as independent and is 2.4-3x larger here).
omega = 2 pi / T, u_sw(omega) = omega u_sw(T) / T.

Weighted least squares of omega on S (S is the controller reading, taken as exact):
  linear     omega = k S + omega0          (the manual's model, reported as asked)
  quadratic  omega = q S^2 + k S + omega0  (used for the conversion)
  cubic      (model check only)
The settings scatter about any smooth curve by more than the stopwatch allows. That scatter is
modelled as a reproducibility of the table speed at a fixed setting, a fraction eps of omega
added in quadrature (u^2 = u_sw^2 + (eps omega)^2), with eps set by chi2 = dof for the
quadratic. The linear model is then tested with these errors. For the measurement settings,
omega has two parts: the curve's covariance (correlated across settings) and eps omega
(independent at each setting, since the speed during the photographs deviates from the curve
the same way). Every fit is reported with chi2, dof and p.

Outputs: calibration.json (conversion model with covariance and eps), calibration.txt,
figures/calibration.pdf/.png.

Run from this folder:  python calibration.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

HERE = Path(__file__).resolve().parent
CSV = HERE / 'Rotation calibration.csv'
CONVERSION_DEGREE = 2          # see calibration.txt for the model comparison behind this
N_SIM = 20000                  # draws for the lag-1 correlation check
SEED = 1


def load_turns(path=CSV):
    """(S, list of turn-time arrays) from the calibration CSV; empty cells are skipped."""
    df = pd.read_csv(path, encoding='utf-8-sig')
    S = df['Speed'].to_numpy(float)
    turns = [row[~np.isnan(row)] for row in df.drop(columns='Speed').to_numpy(float)]
    return S, turns


def per_setting(S, turns):
    rows = []
    for s, t in zip(S, turns):
        n = len(t)
        T, sd = t.mean(), t.std(ddof=1)
        uT = sd / np.sqrt(n * (n + 1))                 # shared stopwatch endpoints
        w = 2 * np.pi / T
        rows.append(dict(S=s, n=n, T=T, sd_T=sd, sigma_press=sd * np.sqrt(n / (2 * (n + 1))),
                         u_T=uT, omega=w, u_sw=w * uT / T))
    return pd.DataFrame(rows)


def lag1(laps):
    """Pooled lag-1 correlation of laps about their own means; laps: list of 1-D arrays."""
    num = sum(np.sum((t[1:] - t.mean()) * (t[:-1] - t.mean())) for t in laps)
    den = sum(np.sum((t - t.mean()) ** 2) for t in laps)
    return num / den


def lag1_simulated(ns, split, rng):
    """Distribution of the pooled lag-1 correlation for laps of lengths ns (vectorised)."""
    num = np.zeros(N_SIM); den = np.zeros(N_SIM)
    for n in ns:
        laps = np.diff(rng.normal(size=(N_SIM, n + 1)), axis=1) if split \
            else rng.normal(size=(N_SIM, n))
        d = laps - laps.mean(axis=1, keepdims=True)
        num += np.sum(d[:, 1:] * d[:, :-1], axis=1); den += np.sum(d ** 2, axis=1)
    return num / den


def polyfit_weighted(x, y, u, deg):
    """Weighted LSQ y = sum c_j x^j (c[0] constant). Returns c, cov, chi2, dof, p."""
    A = np.vander(x, deg + 1, increasing=True)
    W = 1 / u ** 2
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    c = cov @ (A.T @ (W * y))
    chi2 = float(np.sum(((y - A @ c) / u) ** 2))
    dof = len(x) - deg - 1
    return c, cov, chi2, dof, float(stats.chi2.sf(chi2, dof))


def omega_of(S, c):
    """omega at settings S for increasing-order coefficients c."""
    return np.vander(np.atleast_1d(np.asarray(S, float)), len(c), increasing=True) @ c


def omega_with_cov(S, model, reproducibility=False):
    """omega(S) and its covariance from a calibration.json model dict: the curve part, plus
    (eps omega)^2 on the diagonal when reproducibility=True."""
    c = np.array(model['coef']); C = np.array(model['cov'])
    A = np.vander(np.atleast_1d(np.asarray(S, float)), len(c), increasing=True)
    w = A @ c
    cov = A @ C @ A.T
    if reproducibility:
        cov = cov + np.diag((model['eps'] * w) ** 2)
    return w, cov


def make_figure(tab, fits_tot, u_tot, path):
    import matplotlib.pyplot as plt
    import report_style as rs
    rs.apply()
    tex = plt.rcParams['text.usetex']
    unit = r'rad\,s$^{-1}$' if tex else r'rad s$^{-1}$'
    fig, (ax, axr) = plt.subplots(2, 1, figsize=rs.figsize(aspect=1.0), sharex=True,
                                  gridspec_kw=dict(height_ratios=[2.2, 1], hspace=0.06))
    S, w = tab.S.to_numpy(), tab.omega.to_numpy()
    xs = np.linspace(0, 95, 200)
    ax.errorbar(S, w, u_tot, fmt='o', color='k', ms=3, capsize=0, elinewidth=0.6,
                label='measured', zorder=3)
    ax.plot(xs, omega_of(xs, fits_tot[1]['coef']), '-', color=rs.COLORS[0], label='linear fit')
    ax.plot(xs, omega_of(xs, fits_tot[2]['coef']), '--', color=rs.COLORS[1], label='quadratic fit')
    ax.set_ylabel(rf'$\omega$ ({unit})')
    ax.set_xlim(0, 95); ax.set_ylim(0, None)
    ax.legend(loc='upper left')
    for deg, mk, col, dx, lab in ((1, 'o', rs.COLORS[0], -0.6, 'linear'),
                                  (2, 's', rs.COLORS[1], 0.6, 'quadratic')):
        r = (w - omega_of(S, fits_tot[deg]['coef'])) * 1e3
        axr.errorbar(S + dx, r, u_tot * 1e3, fmt=mk, color=col, ms=2.5, capsize=0,
                     elinewidth=0.6, mfc='white' if deg == 2 else col, label=lab)
    axr.axhline(0, color='0.5', lw=0.5)
    axr.legend(loc='lower left', ncol=2)
    axr.set_xlabel('speed setting $S$')
    axr.set_ylabel(r'res.\ ($10^{-3}$ ' + unit + ')' if tex else r'res. ($10^{-3}$ ' + unit + ')')
    fig.align_ylabels((ax, axr))
    rs.save(fig, path)
    plt.close(fig)


def main():
    S, turns = load_turns()
    tab = per_setting(S, turns)
    x, y, usw = tab.S.to_numpy(), tab.omega.to_numpy(), tab.u_sw.to_numpy()

    # Check of the split model: neighbouring laps anticorrelated
    rng = np.random.default_rng(SEED)
    ns = [len(t) for t in turns]
    r_obs = lag1(turns)
    r_split = lag1_simulated(ns, True, rng)
    r_indep = lag1_simulated(ns, False, rng)
    p_indep = float(np.mean(r_indep <= r_obs))

    # Stopwatch errors only
    fits_sw = {}
    for deg in (1, 2, 3):
        c, cov, chi2, dof, p = polyfit_weighted(x, y, usw, deg)
        fits_sw[deg] = dict(coef=c.tolist(), cov=cov.tolist(), chi2=chi2, dof=dof, p=p)

    # Reproducibility eps from the conversion model
    def excess(e):
        f = polyfit_weighted(x, y, np.hypot(usw, e * y), CONVERSION_DEGREE)
        return f[2] - f[3]
    eps = 0.0 if excess(0.0) <= 0 else optimize.brentq(excess, 0.0, 0.5)
    u_tot = np.hypot(usw, eps * y)
    fits_tot = {}
    for deg in (1, 2, 3):
        c, cov, chi2, dof, p = polyfit_weighted(x, y, u_tot, deg)
        fits_tot[deg] = dict(coef=c.tolist(), cov=cov.tolist(), chi2=chi2, dof=dof, p=p)

    conv = fits_tot[CONVERSION_DEGREE]
    model = dict(degree=CONVERSION_DEGREE, coef=conv['coef'], cov=conv['cov'], eps=eps,
                 S_range=[float(x.min()), float(x.max())],
                 note='omega in rad/s; coef in increasing order (constant first); cov is the '
                      'curve covariance (correlated across settings); eps*omega is the '
                      'independent speed reproducibility at a setting')
    lin = fits_tot[1]
    out = dict(conversion=model,
               linear=dict(k=lin['coef'][1], u_k=float(np.sqrt(lin['cov'][1][1])),
                           omega0=lin['coef'][0], u_omega0=float(np.sqrt(lin['cov'][0][0])),
                           chi2=lin['chi2'], dof=lin['dof'], p=lin['p'],
                           chi2_stopwatch_only=fits_sw[1]['chi2'], p_stopwatch_only=fits_sw[1]['p']),
               fits_stopwatch_only={str(k): v for k, v in fits_sw.items()},
               fits_with_eps={str(k): v for k, v in fits_tot.items()},
               lag1=dict(observed=r_obs, split_mean=float(r_split.mean()),
                         split_sd=float(r_split.std()), p_independent=p_indep))
    (HERE / 'calibration.json').write_text(json.dumps(out, indent=1))

    L = ['Rotation calibration (calibration.py)\n']
    L.append('Stopwatch splits: u_sw(T) = SD/sqrt(n(n+1)) (shared endpoints; statistical). '
             'sigma_press = SD sqrt(n/(2(n+1))) is the implied timing error of one press.\n')
    L.append(f"{'S':>6} {'n':>2} {'T (s)':>8} {'SD (s)':>7} {'s_press':>7} {'u_sw(T)':>8} "
             f"{'omega (rad/s)':>13} {'u_sw(omega)':>11} {'rel':>6}")
    for r in tab.itertuples():
        L.append(f'{r.S:6.2f} {r.n:2d} {r.T:8.3f} {r.sd_T:7.3f} {r.sigma_press:7.3f} {r.u_T:8.4f} '
                 f'{r.omega:13.5f} {r.u_sw:11.5f} {100 * r.u_sw / r.omega:5.2f}%')
    L.append(f'\nLap correlation check: pooled lag-1 correlation {r_obs:+.2f}; split model '
             f'(simulated, same n) {r_split.mean():+.2f} +- {r_split.std():.2f}; independent laps '
             f'would give {r_obs:+.2f} or lower with p = {p_indep:.2f}. Consistent with splits.\n')
    names = {1: 'linear   ', 2: 'quadratic', 3: 'cubic    '}
    for title, fits in (('stopwatch errors only', fits_sw),
                        (f'with reproducibility eps = {100 * eps:.2f} % of omega', fits_tot)):
        L.append(f'Weighted fits omega(S), {title}:')
        for deg, f in fits.items():
            c, e = np.array(f['coef']), np.sqrt(np.diag(f['cov']))
            terms = ', '.join(f'c{j} = {c[j]:.6g} +- {e[j]:.2g}' for j in range(len(c)))
            L.append(f"  {names[deg]} chi2 = {f['chi2']:.1f}, dof = {f['dof']}, "
                     f"chi2/dof = {f['chi2'] / f['dof']:.2f}, p = {f['p']:.2g}   [{terms}]")
        L.append('')
    L.append('(c0 in rad/s, c1 in rad/s per unit S, c2 in rad/s per unit S^2, ...)')
    L.append(f'eps = {100 * eps:.2f} % is set by chi2 = dof for the quadratic: the settings scatter '
             'about the smooth curve by this fraction of omega, far more than the stopwatch allows '
             f'({100 * np.min(usw / y):.2f}-{100 * np.max(usw / y):.2f} %). The cubic does not '
             'remove it, so it is scatter, not a missing term.\n')
    z = (y - omega_of(x, lin['coef'])) / u_tot
    L.append(f"Linear ratio (the manual's task 1): k = d omega / dS = {lin['coef'][1]:.5f} +- "
             f"{np.sqrt(lin['cov'][1][1]):.5f} rad/s per unit S (fit, with eps), "
             f"omega0 = {lin['coef'][0]:.4f} +- {np.sqrt(lin['cov'][0][0]):.4f} rad/s.")
    L.append(f"The linear model is rejected even with eps (chi2 = {lin['chi2']:.1f}, dof = {lin['dof']}, "
             f"p = {lin['p']:.1g}): normalised residuals {', '.join(f'{v:+.1f}' for v in z)} "
             '(in order of S); the table turns slower than linear at high settings.')
    r_lin = y - omega_of(x, lin['coef'])
    L.append(f'Largest linear-model residual: {r_lin[np.argmax(np.abs(r_lin))]:+.4f} rad/s '
             f'({100 * np.max(np.abs(r_lin) / y):.1f} % of omega).')
    L.append('\nConversion used by sagnac.py (quadratic): omega +- curve (correlated) '
             '+- reproducibility eps omega (independent):')
    for s in (12.59, 20.74, 31.11, 44.07, 55.83, 64.00, 73.65, 82.18, 90.0):
        w, C = omega_with_cov([s], model)
        L.append(f'  S = {s:6.2f}: omega = {w[0]:.4f} +- {np.sqrt(C[0, 0]):.4f} (curve) '
                 f'+- {eps * w[0]:.4f} (reproducibility) rad/s')
    (HERE / 'calibration.txt').write_text('\n'.join(L) + '\n', encoding='utf-8')
    print('\n'.join(L))
    make_figure(tab, fits_tot, u_tot, HERE / 'figures' / 'calibration')


if __name__ == '__main__':
    main()
