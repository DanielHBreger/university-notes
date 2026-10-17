"""Rotation calibration: how fast the table turns at each controller setting S.

WHAT THIS SCRIPT IS FOR
The Sagnac phase is proportional to the table's angular velocity omega, but during the
photographs only the controller setting S is known. This script turns the stopwatch data in
"Rotation calibration.csv" into a function omega(S), with uncertainties, which sagnac.py uses
later. It also answers the manual's task 1 (the linear ratio between setting and omega).

HOW IT RUNS (main(), top to bottom)
  1. load_turns     read the CSV: for each setting S, the times of 5 or 8 consecutive turns.
  2. per_setting    mean turn time T, its stopwatch uncertainty, and omega = 2 pi / T.
  3. lag1 check     test that the turn times really are stopwatch "splits" (see below).
  4. fits           fit omega(S) with a straight line, a parabola and a cubic, first with the
                    stopwatch errors only.
  5. eps            the points scatter about any smooth curve more than the stopwatch explains;
                    find the run-to-run speed variation eps (a fraction of omega) that makes the
                    parabola fit statistically acceptable (chi2 = dof).
  6. refit          refit all three curves with stopwatch + eps errors and compare them.
  7. save           calibration.json (the parabola, its covariance and eps, read by sagnac.py),
                    calibration.txt (tables) and figures/calibration.

KEY IDEAS
Stopwatch splits. The turns were timed by pressing "lap" at each pass, so turn i lasts
t_i - t_(i-1) and neighbouring turns share a button press. Their mean is (t_n - t_0)/n: only
the first and last presses matter. If every press has an independent error sigma_r, then
Var(mean) = 2 sigma_r^2 / n^2, while the spread of the individual turns obeys
E[SD^2] = 2 sigma_r^2 (n+1)/n. Eliminating sigma_r:
    u_sw(T) = SD / sqrt(n (n+1)),
which is 2.4-3x smaller than the usual SD/sqrt(n) (that formula assumes independent turns).
A testable consequence: neighbouring split laps are anticorrelated (a late press makes one lap
long and the next one short). Step 3 compares the measured correlation with simulations of
both cases.

Model choice (steps 4-6). "Model" means the curve used for omega(S):
  linear     omega = omega0 + k S            (what the manual asks for)
  quadratic  omega = omega0 + k S + q S^2    (used for the conversion)
  cubic      one more term                   (only to check whether the quadratic is enough)
Each is fitted by weighted least squares (points with smaller errors count more). Its chi2 is
the sum of (residual / error)^2. A curve that describes the data, with correct errors, gives
chi2 close to dof = (number of points - number of parameters); p is the probability of a chi2
at least this large if the curve were right. With stopwatch errors alone every curve fails
(p << 0.05): the table does not reproduce its speed that precisely. That extra scatter, eps,
is added in quadrature, u^2 = u_sw^2 + (eps omega)^2, and chosen so that the quadratic has
chi2 = dof. With these errors the straight line is still rejected (the table runs slower than
linear at high S) and the cubic gains nothing, so the quadratic is used.

Two kinds of omega uncertainty at the photographed settings (omega_with_cov):
  curve            from the fitted coefficients; it moves all settings together (correlated),
                   so sagnac.py propagates it through the slope separately;
  reproducibility  eps * omega; independent at each setting, because the table's actual speed
                   during the photographs differs from the curve by a random amount.

Run from this folder:  python calibration.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

HERE = Path(__file__).resolve().parent
CSV = HERE / 'Rotation calibration.csv'
CONVERSION_DEGREE = 2          # the quadratic; see "Model choice" above and calibration.txt
N_SIM = 20000                  # simulated data sets for the lap-correlation check (step 3)
SEED = 1                       # fixed random seed, so the check gives the same numbers every run


def load_turns(path=CSV):
    """Step 1. Read the CSV: column 'Speed' is the setting S, columns 'turn 1'...'turn 10' are
    the times (s) of single turns. Settings have 5 or 8 turns; the empty cells (NaN) are dropped.
    Returns S (array) and a list with one array of turn times per setting."""
    df = pd.read_csv(path, encoding='utf-8-sig')       # utf-8-sig: the file starts with a BOM
    S = df['Speed'].to_numpy(float)
    turns = [row[~np.isnan(row)] for row in df.drop(columns='Speed').to_numpy(float)]
    return S, turns


def per_setting(S, turns):
    """Step 2. One row per setting: mean turn time T, the spread SD of single turns, the
    stopwatch uncertainty u_T of the mean (split formula, see the module docstring), the
    implied error of one button press, omega = 2 pi / T and its stopwatch uncertainty u_sw."""
    rows = []
    for s, t in zip(S, turns):
        n = len(t)
        T, sd = t.mean(), t.std(ddof=1)                # ddof=1: sample SD (divides by n-1)
        uT = sd / np.sqrt(n * (n + 1))                 # shared stopwatch endpoints
        w = 2 * np.pi / T
        # sigma_press: the per-press timing error implied by SD (from E[SD^2] = 2 sigma^2 (n+1)/n)
        # u_sw: error propagation of omega = 2 pi / T, i.e. |d omega / dT| u_T = omega u_T / T
        rows.append(dict(S=s, n=n, T=T, sd_T=sd, sigma_press=sd * np.sqrt(n / (2 * (n + 1))),
                         u_T=uT, omega=w, u_sw=w * uT / T))
    return pd.DataFrame(rows)


def lag1(laps):
    """Step 3 (measured). How strongly each turn time correlates with the next one, pooled over
    all settings: sum of (lap_i - mean)(lap_(i+1) - mean) over sum of (lap_i - mean)^2.
    Split timing predicts a clearly negative value; independent laps give a value near 0
    (slightly negative for short series)."""
    num = sum(np.sum((t[1:] - t.mean()) * (t[:-1] - t.mean())) for t in laps)
    den = sum(np.sum((t - t.mean()) ** 2) for t in laps)
    return num / den


def lag1_simulated(ns, split, rng):
    """Step 3 (expected). The same pooled correlation for N_SIM simulated data sets with the
    same numbers of turns ns. split=True: laps are differences of independent press times (the
    split model); split=False: independent laps. Vectorised: each row of the arrays is one
    simulated data set. Returns N_SIM correlation values (their distribution)."""
    num = np.zeros(N_SIM); den = np.zeros(N_SIM)
    for n in ns:
        laps = np.diff(rng.normal(size=(N_SIM, n + 1)), axis=1) if split \
            else rng.normal(size=(N_SIM, n))
        d = laps - laps.mean(axis=1, keepdims=True)
        num += np.sum(d[:, 1:] * d[:, :-1], axis=1); den += np.sum(d ** 2, axis=1)
    return num / den


def polyfit_weighted(x, y, u, deg):
    """Weighted least squares of a polynomial y = c0 + c1 x + ... + c_deg x^deg, each point
    weighted by 1/u^2. Returns the coefficients c (constant first), their covariance matrix,
    chi2, the degrees of freedom and p (probability of a chi2 this large if the model is right).

    The solution is the textbook one, written with the design matrix A whose columns are
    1, x, x^2, ...:  cov = (A^T W A)^-1  and  c = cov A^T W y,  with W = diag(1/u^2)."""
    A = np.vander(x, deg + 1, increasing=True)         # columns: x^0, x^1, ..., x^deg
    W = 1 / u ** 2
    cov = np.linalg.inv(A.T @ (A * W[:, None]))         # A * W[:, None] = W A without building diag(W)
    c = cov @ (A.T @ (W * y))
    chi2 = float(np.sum(((y - A @ c) / u) ** 2))
    dof = len(x) - deg - 1
    return c, cov, chi2, dof, float(stats.chi2.sf(chi2, dof))   # sf = 1 - CDF = P(chi2 >= value)


def omega_of(S, c):
    """The polynomial with coefficients c (constant first) evaluated at the settings S."""
    return np.vander(np.atleast_1d(np.asarray(S, float)), len(c), increasing=True) @ c


def omega_with_cov(S, model, reproducibility=False):
    """omega at the settings S from the saved conversion model (calibration.json, 'conversion'),
    with its covariance matrix. Used by sagnac.py for the photographed settings.

    omega = A c, where A is the design matrix at S (columns 1, S, S^2). The coefficients' own
    covariance C propagates as cov(omega) = A C A^T: a full matrix, because all settings share
    the same coefficients (if the curve is a bit high at one setting, it is high at its
    neighbours too). With reproducibility=True the independent part (eps omega)^2 is added on
    the diagonal only."""
    c = np.array(model['coef']); C = np.array(model['cov'])
    A = np.vander(np.atleast_1d(np.asarray(S, float)), len(c), increasing=True)
    w = A @ c
    cov = A @ C @ A.T
    if reproducibility:
        cov = cov + np.diag((model['eps'] * w) ** 2)
    return w, cov


def make_figure(tab, fits_tot, u_tot, path):
    """figures/calibration: measured omega with the linear and quadratic fits (top), and each
    fit's residuals in 10^-3 rad/s (bottom). Error bars include eps."""
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
    # 1-2. Read the turn times; per setting: T, its stopwatch error, omega and its error.
    S, turns = load_turns()
    tab = per_setting(S, turns)
    x, y, usw = tab.S.to_numpy(), tab.omega.to_numpy(), tab.u_sw.to_numpy()

    # 3. Check the split model: neighbouring laps should be anticorrelated. Compare the measured
    #    correlation with simulated split laps and with simulated independent laps; p_indep is
    #    the fraction of independent-lap simulations at least as negative as the measurement.
    rng = np.random.default_rng(SEED)
    ns = [len(t) for t in turns]
    r_obs = lag1(turns)
    r_split = lag1_simulated(ns, True, rng)
    r_indep = lag1_simulated(ns, False, rng)
    p_indep = float(np.mean(r_indep <= r_obs))

    # 4. Fit line, parabola and cubic with the stopwatch errors only. All three come out with
    #    chi2 far above dof: the points scatter more than the stopwatch explains.
    fits_sw = {}
    for deg in (1, 2, 3):
        c, cov, chi2, dof, p = polyfit_weighted(x, y, usw, deg)
        fits_sw[deg] = dict(coef=c.tolist(), cov=cov.tolist(), chi2=chi2, dof=dof, p=p)

    # 5. The speed reproducibility eps: excess(e) = chi2 - dof of the parabola when every point's
    #    error is sqrt(u_sw^2 + (e omega)^2). brentq finds the e where it crosses zero, i.e. the
    #    extra scatter that makes the parabola's chi2 equal to its dof.
    def excess(e):
        f = polyfit_weighted(x, y, np.hypot(usw, e * y), CONVERSION_DEGREE)
        return f[2] - f[3]                              # f = (c, cov, chi2, dof, p)
    eps = 0.0 if excess(0.0) <= 0 else optimize.brentq(excess, 0.0, 0.5)
    u_tot = np.hypot(usw, eps * y)                      # hypot(a, b) = sqrt(a^2 + b^2)

    # 6. Refit all three curves with these errors. Now the parabola has chi2 = dof by
    #    construction, the line is still rejected, and the cubic does not improve on the parabola.
    fits_tot = {}
    for deg in (1, 2, 3):
        c, cov, chi2, dof, p = polyfit_weighted(x, y, u_tot, deg)
        fits_tot[deg] = dict(coef=c.tolist(), cov=cov.tolist(), chi2=chi2, dof=dof, p=p)

    # 7. Save. 'conversion' is what sagnac.py reads: the parabola's coefficients, their
    #    covariance and eps. 'linear' is the manual's task-1 answer. Then the text report.
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

    # Text report: per-setting table, the lap-correlation check, every fit with chi2/dof/p,
    # the linear ratio, and omega at the photographed settings with its two error parts.
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
    z = (y - omega_of(x, lin['coef'])) / u_tot          # line residuals in units of their error
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
