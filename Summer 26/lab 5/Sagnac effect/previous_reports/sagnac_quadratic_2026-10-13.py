"""Sagnac phase against angular velocity: the slope, the comparison with Eq. (5), and c.

WHAT THIS SCRIPT IS FOR
fringes.py (and methods.py) give a phase for every photograph; calibration.py gives omega for
every setting. This script averages the photographs of each setting, fits a straight line
phase = a omega + b, compares the slope a with the manual's Eq. (5), solves Eq. (5) for the
speed of light, and writes all numbers with their uncertainties.
Inputs: fringes_frames.csv (fringes.py), methods_frames.csv (methods.py, if present) and
calibration.json (calibration.py).   python sagnac.py

HOW IT RUNS (main(), top to bottom)
  1. setting_means    per setting: mean phase of its photographs and its standard error.
  2. calibration_part omega for every setting from the calibration curve.
  3. fit_line         weighted straight line through (omega, mean phase); finds sigma_x.
  4. calibration_part again, now for the slope's uncertainty from the calibration curve.
  5. checks           a quadratic term (is the dependence linear?) and ring_by_ring.
  6. theory           a_th from Eq. (5); ratio, difference and its significance.
  7. other methods    the same averaging and fit for the cosine fit, template matching,
                      2D Fourier and AIA; their spread gives the method systematic u_method.
  8. drift test       phase = a omega + D t + b, to see if a drift in time could mimic a.
  9. outputs          sagnac_speeds.csv, then the text report sagnac.txt; c (speed_of_light,
                      from the measured slope with every uncertainty part) is computed while
                      the report is written; then sagnac.json and the figures phase_vs_omega
                      and methods_comparison.

KEY IDEAS
Phase averaging. Phases are angles, so they are only known modulo 2 pi. Within a setting the
photographs are first brought next to their circular mean (the angle of the average of
e^(i phi)), then averaged. Between settings, each mean is shifted by a multiple of 2 pi to lie
within pi of the previous setting's ("unwrapping"), which works because neighbouring settings
differ by much less than pi. omega(S = 0) is exactly 0 (the table at rest).

The weighted line and its errors. Each setting's point gets the error
    u^2 = SEM^2 + (a eps omega)^2 + sigma_x^2
  SEM          scatter of its photographs, SD / sqrt(n);
  a eps omega  the table's speed reproducibility (calibration.py) turned into phase: an error
               delta omega in x shifts the line's y by a delta omega ("effective variance");
               since it depends on a, the fit is iterated until a stops changing;
  sigma_x      extra scatter between settings, the same for all. Without it the line is
               rejected (chi2 = 62 for 8 dof): the phase wanders on minutes-long time scales,
               longer than the 20-60 s spent at each setting. sigma_x is set so that
               chi2 = dof, which makes the slope's error honest.
The slope's uncertainty has two parts:
  u_stat  from the fit (phase scatter and speed reproducibility),
  u_cal   from the calibration curve: the same curve error moves every omega together, so it is
          propagated separately (the derivative of the slope with respect to each calibration
          coefficient, combined with the coefficients' covariance).
The other methods go through the same averaging and fit, and are also compared frame by frame
with the minima (differences on the same photographs cancel their shared noise). Half the range
of the reliable methods' slopes is the method systematic u_method.

Theory, Eq. (5) with A = pi R^2:  a_th = 8 pi N A / (lambda c) = 8 pi^2 N R^2 / (lambda c),
independent of where the rotation axis lies (Schwartz 2017).
Speed of light: Eq. (5) solved for c with the measured slope, c = 8 pi^2 N R^2 / (lambda a), the
slope's parts (fit, calibration curve, method) and N, R carried as separate components.

Outputs: sagnac.txt, sagnac.json, sagnac_speeds.csv, figures/phase_vs_omega.pdf/.png,
figures/methods_comparison.pdf/.png. Run from this folder after calibration.py, fringes.py and
methods.py.
"""
import json
import time
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

from calibration import omega_with_cov

HERE = Path(__file__).resolve().parent

# Fibre loop, centred on the rotation axis (values given by the experimenters; +- taken as
# standard uncertainties)
N_TURNS, U_N_TURNS = 13, 1
RADIUS, U_RADIUS = 0.570, 0.005       # m
# Other phase methods (label, frame column, counted in the method spread). AIA is not counted:
# its slope depends on its window and variant here (methods.txt).
OTHER_METHODS = [('cosine fit', 'phi_cos_rad', True), ('template matching', 'phi_template_rad', True),
                 ('2D Fourier', 'phi_fourier_rad', True), ('AIA', 'phi_aia_rad', False)]
WAVELENGTH = 640e-9                   # m
U_WAVELENGTH = None                   # not given by the lab: left out of the budget, stated in the output
C_LIGHT = 299_792_458.0               # m/s, exact


def wrap(x):
    """Angle(s) wrapped into (-pi, pi]."""
    return np.angle(np.exp(1j * np.asarray(x)))


def round1(x, u):
    """'x +- u' with u to one significant figure and x to the same digit, rounding half up on the
    decimal digits as printed (float round() sends ties to even: 2.25 -> 2.2, and 1.005 -> 1.0).
    Decimal works on the printed digits, so there are no binary rounding surprises."""
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
    """Step 6. The predicted slope a_th = 8 pi^2 N R^2 / (lambda c) (rad per rad/s), and its
    uncertainty components from N and R. a_th is proportional to N and to R^2, so the relative
    errors are u_N/N and 2 u_R/R."""
    a = 8 * np.pi ** 2 * N_TURNS * RADIUS ** 2 / (WAVELENGTH * C_LIGHT)
    comps = {'N': a * U_N_TURNS / N_TURNS, 'R': a * 2 * U_RADIUS / RADIUS}
    if U_WAVELENGTH is not None:
        comps['lambda'] = a * U_WAVELENGTH / WAVELENGTH
    return a, comps


def speed_of_light(a, u_a_parts):
    """Step 8. Eq. (5) solved for c with the measured slope a: c = 8 pi^2 N R^2 / (lambda a).
    u_a_parts: the slope's uncertainty parts by kind, rad/(rad/s). Returns c (m/s) and its uncertainty
    components (m/s): the slope's parts (c ~ 1/a, so the same relative size), then N and R."""
    c = 8 * np.pi ** 2 * N_TURNS * RADIUS ** 2 / (WAVELENGTH * a)
    comps = {k: c * u / a for k, u in u_a_parts.items()}
    comps['N'] = c * U_N_TURNS / N_TURNS
    comps['R'] = c * 2 * U_RADIUS / RADIUS
    if U_WAVELENGTH is not None:
        comps['lambda'] = c * U_WAVELENGTH / WAVELENGTH
    return c, comps


def setting_means(fr, col='phi_rad'):
    """Step 1. One row per setting from the per-photograph phases in column `col`: number of
    photographs n, mean phase, SD, SEM = SD/sqrt(n), start time and time span, the phase drift
    within the setting (slope of phase against time, with its error and p), and the median
    fringe period P. The setting means are then unwrapped in order of S (see KEY IDEAS)."""
    rows = []
    for S, g in fr.groupby('S', sort=True):            # g: the photographs of setting S
        ph = g[col].to_numpy()
        c = np.angle(np.mean(np.exp(1j * ph)))          # circular mean: angle of the mean unit vector
        ph = c + wrap(ph - c)                           # each phase within pi of it, then plain averages work
        t = g.time_s.to_numpy() - g.time_s.min()
        n = len(ph)
        lr = stats.linregress(t, ph)                    # drift within the setting (rad/s)
        rows.append(dict(S=S, n=n, phi=ph.mean(), sd=ph.std(ddof=1), sem=ph.std(ddof=1) / np.sqrt(n),
                         t_start=g.time_s.min(), span_s=t.max(), drift=lr.slope,
                         u_drift=lr.stderr, p_drift=lr.pvalue, P=g.P_px2.median()))
    tab = pd.DataFrame(rows)
    phi = tab.phi.to_numpy().copy()
    for i in range(1, len(phi)):                       # unwrap in order of S
        # add the multiple of 2 pi that brings this setting within pi of the previous one
        phi[i] += 2 * np.pi * np.round((phi[i - 1] - phi[i]) / (2 * np.pi))
    tab['phi'] = phi
    return tab


def wls(A, y, u):
    """Weighted linear least squares y = A p (design matrix A, errors u): the same normal
    equations as calibration.polyfit_weighted. Returns p, its covariance, and chi2."""
    W = 1 / u ** 2
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    p = cov @ (A.T @ (W * y))
    return p, cov, float(np.sum(((y - A @ p) / u) ** 2))


def wls_line(x, y, u):
    """Weighted y = a x + b. Returns (a, b), covariance, chi2."""
    return wls(np.c_[x, np.ones_like(x)], y, u)        # columns: x, 1


def extra_scatter(x, y, u_base):
    """The extra scatter sigma_x >= 0 that, added in quadrature to every point's error u_base,
    makes the line's chi2 equal to its dof (= points - 2). 0 if chi2 <= dof already. Found as
    the root of chi2 - dof with brentq (chi2 falls steadily as sigma_x grows)."""
    dof = len(x) - 2
    f = lambda s: wls_line(x, y, np.sqrt(u_base ** 2 + s ** 2))[2] - dof
    return 0.0 if f(0.0) <= 0 else optimize.brentq(f, 0.0, 10.0)


def fit_line(omega, y, sem, eps):
    """Step 3. The weighted line phase = a omega + b with per-point error
    u_fit^2 = SEM^2 + (a eps omega)^2 + sigma_x^2 (see KEY IDEAS). The middle term needs a, and
    sigma_x depends on all errors, so: start from a fit with the SEMs alone, then repeat
    (errors from the current a -> sigma_x -> refit -> new a) until a no longer changes.
    Returns a, b, their covariance, chi2 (= dof when sigma_x > 0), sigma_x, the errors without
    and with sigma_x, and the chi2 the line would have without sigma_x."""
    a = wls_line(omega, y, sem)[0][0]
    for _ in range(50):
        u_base = np.sqrt(sem ** 2 + (a * eps * omega) ** 2)
        sx = extra_scatter(omega, y, u_base)
        u_fit = np.sqrt(u_base ** 2 + sx ** 2)
        (a_new, b), cov, chi = wls_line(omega, y, u_fit)
        done = abs(a_new - a) < 1e-12
        a = a_new
        if done:
            break
    chi_base = wls_line(omega, y, u_base)[2]
    return a, b, cov, chi, sx, u_base, u_fit, chi_base


def calibration_part(S, y, u_fit, model):
    """Steps 2 and 4. omega for every setting (0 at rest) from the calibration model, with its
    curve error u_curve and reproducibility error u_rep (= eps omega).

    With u_fit given, also the slope's uncertainty from the calibration curve, u_cal. The slope
    is a function of the three calibration coefficients c (through omega); its derivatives J are
    taken numerically (central differences, the fit's weights held fixed), and the coefficients'
    covariance Cc propagates as u_cal^2 = J Cc J^T. This is done once for the slope instead of
    per point because the curve error is shared by all settings."""
    moving = S > 0
    omega = np.zeros(len(S)); u_curve = np.zeros(len(S))
    w, Cw = omega_with_cov(S[moving], model)
    omega[moving] = w
    u_curve[moving] = np.sqrt(np.diag(Cw))
    u_rep = model['eps'] * omega
    coef = np.array(model['coef']); Cc = np.array(model['cov'])

    def slope(c):                    # fit weights held fixed
        om = np.zeros(len(S))
        om[moving] = np.vander(S[moving], len(c), increasing=True) @ c
        return wls_line(om, y, u_fit)[0][0]

    if u_fit is None:
        return omega, u_curve, u_rep, None
    # one small step h per coefficient (relative 1e-6), slope difference over 2h
    J = np.array([(slope(coef + h) - slope(coef - h)) / (2 * h[j])
                  for j, h in enumerate(np.diag(np.abs(coef) * 1e-6 + 1e-12))])
    return omega, u_curve, u_rep, float(np.sqrt(J @ Cc @ J))


def ring_by_ring(fr, omega):
    """Step 5 check: manual Eq. (7) applied to each dark ring separately, with one period P for
    all photographs (the median), instead of the centre phase: dphi = 2 pi (s_rest - s) / P for
    ring k, s = mean R_k^2 per setting. If the period changes between settings, the rings
    disagree with each other; the centre phase is not affected (see fringes.py)."""
    P = float(fr.P_px2.median())
    out = {}
    for k in (1, 2, 3):
        col = f'R{k}'
        s = (pd.to_numeric(fr[col], errors='coerce') ** 2).groupby(fr.S).mean()   # NaN-skipping
        if not np.isfinite(s.iloc[0]):          # ring not seen at rest: no reference
            continue
        d = 2 * np.pi * (s.iloc[0] - s.to_numpy()) / P
        ok = np.isfinite(d)
        out[col] = dict(dphi=d.tolist(), slope=float(np.polyfit(omega[ok], d[ok], 1)[0]))
    return P, out


def make_figure(tab, omega, u_omega, u_fit, a, b, a_th, u_th, path, cos):
    """figures/phase_vs_omega: phase shift relative to rest against omega for the minima (filled)
    and the cosine fit (open), each with its fitted line, and Eq. (5) with its 1-sigma band;
    residuals about each line below. cos = (phases, errors, a, b) of the cosine fit."""
    import matplotlib.pyplot as plt
    import report_style as rs
    rs.apply()
    tex = plt.rcParams['text.usetex']
    unit = r'rad\,s$^{-1}$' if tex else r'rad s$^{-1}$'
    fig, (ax, axr) = plt.subplots(2, 1, figsize=rs.figsize(aspect=1.0), sharex=True,
                                  gridspec_kw=dict(height_ratios=[2.2, 1], hspace=0.06))
    ref = tab.phi.iloc[0]
    d = tab.phi.to_numpy() - ref
    xs = np.linspace(0, 1.05 * omega.max(), 100)
    band = ax.fill_between(xs, (a_th - u_th) * xs, (a_th + u_th) * xs, color=rs.COLORS[1],
                           alpha=0.18, lw=0)
    th_line, = ax.plot(xs, a_th * xs, '--', color=rs.COLORS[1])
    fit_line_, = ax.plot(xs, a * xs + b - ref, '-', color=rs.COLORS[0])
    meas = ax.errorbar(omega, d, u_fit, u_omega, fmt='o', color=rs.COLORS[0], ms=3, capsize=0,
                       elinewidth=0.6, zorder=3)
    y_c, u_fit_c, a_c, b_c = cos
    cos_line, = ax.plot(xs, a_c * xs + b_c - y_c[0], '-.', color=rs.COLORS[2])
    cos_pts = ax.errorbar(omega, y_c - y_c[0], u_fit_c, fmt='s', color=rs.COLORS[2], mfc='white',
                          ms=3, capsize=0, elinewidth=0.6, zorder=3)
    ax.set_ylabel(r'$\Delta\phi$ (rad)')
    ax.set_xlim(-0.03 * xs[-1], xs[-1])
    ax.legend([(meas, fit_line_), (cos_pts, cos_line), (band, th_line)],
              ['minima', 'cosine fit', r'Eq.~(5), $\pm1\sigma$' if tex else r'Eq. (5), $\pm1\sigma$'],
              loc='upper left')
    dx = 0.01 * xs[-1]                                    # side by side in the residual panel only
    axr.errorbar(omega - dx, tab.phi.to_numpy() - (a * omega + b), u_fit, u_omega, fmt='o',
                 color=rs.COLORS[0], ms=3, capsize=0, elinewidth=0.6)
    axr.errorbar(omega + dx, y_c - (a_c * omega + b_c), u_fit_c, fmt='s', color=rs.COLORS[2],
                 mfc='white', ms=3, capsize=0, elinewidth=0.6)
    axr.axhline(0, color='0.5', lw=0.5)
    axr.set_xlabel(rf'$\omega$ ({unit})')
    axr.set_ylabel('res. (rad)')
    fig.align_ylabels((ax, axr))
    rs.save(fig, path)
    plt.close(fig)


def make_comparison(omega, others, a, u_a, a_th, u_th, scan, path):
    """(a) each method minus the minima per setting, relative to rest; (b) the slopes of all methods
    (fit and calibration errors), the range of each robustness scan (grey) and Eq. (5)."""
    import matplotlib.pyplot as plt
    import report_style as rs
    rs.apply()
    tex = plt.rcParams['text.usetex']
    unit = r'rad\,s$^{-1}$' if tex else r'rad s$^{-1}$'
    style = {'cosine fit': (rs.COLORS[2], 's'), 'template matching': (rs.COLORS[3], 'D'),
             '2D Fourier': (rs.COLORS[4], 'v'), 'AIA': (rs.COLORS[5], '^')}
    fig, (ax, axb) = plt.subplots(2, 1, figsize=rs.figsize(aspect=1.15),
                                  gridspec_kw=dict(height_ratios=[1.5, 1], hspace=0.42))
    ax.axhline(0, color=rs.COLORS[0], lw=0.8, label='minima')
    dx = 0.012 * omega.max()
    for i, (label, m) in enumerate(others.items()):
        col, mk = style[label]
        ax.errorbar(omega + (i - (len(others) - 1) / 2) * dx, m['d_setting'], m['u_d_setting'], fmt=mk, color=col,
                    ms=3, capsize=0, elinewidth=0.6, label=label)
    ax.set_xlabel(rf'$\omega$ ({unit})')
    ax.set_ylabel(r'$\Delta\phi - \Delta\phi_{\rm minima}$ (rad)')
    ax.legend(loc='lower left', ncol=2)

    rows = [('minima', a, u_a, rs.COLORS[0], 'o', None)]
    ranges = {'template matching': [r['slope'] for r in scan.get('template_scan', [])],
              '2D Fourier': [r['slope'] for r in scan.get('fourier_scan', [])],
              'AIA': [r['slope'] for r in scan.get('aia_scan', [])]}
    rows += [(l, m['a'], m['u_a'], *style[l], ranges.get(l)) for l, m in others.items()]
    axb.axvspan(a_th - u_th, a_th + u_th, color=rs.COLORS[1], alpha=0.18, lw=0)
    axb.axvline(a_th, color=rs.COLORS[1], ls='--', lw=0.8)
    for yk, (label, val, err, col, mk, rng) in enumerate(rows):
        if rng:
            axb.plot([min(rng), max(rng)], [yk, yk], '-', color='0.7', lw=2.5, solid_capstyle='butt', zorder=1)
        axb.errorbar(val, yk, xerr=err, fmt=mk, color=col, ms=3.5, capsize=0, elinewidth=0.8, zorder=3)
    axb.set_yticks(range(len(rows)), [r[0] for r in rows])
    axb.tick_params(axis='y', which='both', length=0)
    axb.minorticks_off(); axb.xaxis.set_minor_locator(plt.matplotlib.ticker.AutoMinorLocator())
    axb.set_ylim(len(rows) - 0.4, -0.6)
    axb.set_xlabel(rf'slope $a$ (rad per {unit})')
    axb.text(a_th, len(rows) - 0.45, r' Eq.~(5)' if tex else ' Eq. (5)', color=rs.COLORS[1], va='bottom', fontsize=7)
    for axx, lab in ((ax, '(a)'), (axb, '(b)')):
        axx.text(0.02, 0.96, lab, transform=axx.transAxes, va='top', ha='left')
    rs.save(fig, path)
    plt.close(fig)


def hhmm(t):
    """Seconds since midnight as 'hh:mm'."""
    return time.strftime('%H:%M', time.gmtime(t))


def main():
    src = HERE / 'fringes_frames.csv'
    if not src.exists():
        raise SystemExit(f'{src.name} not found: run fringes.py first (see README.md).')

    # 1. Per-photograph phases (minima method) -> one mean phase and SEM per setting.
    fr = pd.read_csv(src)
    model = json.loads((HERE / 'calibration.json').read_text())['conversion']
    eps = model['eps']
    tab = setting_means(fr)
    S = tab.S.to_numpy()
    sem = tab['sem'].to_numpy()
    y = tab.phi.to_numpy()
    dof = len(y) - 2                                   # 10 settings - 2 line parameters

    # 2. omega of every setting (and its two error parts) from the calibration.
    omega, u_curve, u_rep, _ = calibration_part(S, y, None, model)
    # 3. The weighted line. chi_base: chi2 without sigma_x (shows why sigma_x is needed);
    #    p_sem_only: the same line with only the frame SEMs as errors, for comparison.
    a, b, cov, chi, sx, u_base, u_fit, chi_base = fit_line(omega, y, sem, eps)
    p_base = stats.chi2.sf(chi_base, dof)
    p_sem_only = wls_line(omega, y, sem)
    u_stat = float(np.sqrt(cov[0, 0]))
    # 4. The slope's error from the calibration curve, then the total (in quadrature).
    _, _, _, u_cal = calibration_part(S, y, u_fit, model)
    u_a = float(np.hypot(u_stat, u_cal))
    u_omega = np.hypot(u_curve, u_rep)

    # 5. Curvature check: a quadratic term with the same errors; zq = how many sigma from zero.
    q, covq, _ = wls(np.c_[omega ** 2, omega, np.ones_like(omega)], y, u_fit)
    zq = q[0] / np.sqrt(covq[0, 0])

    # 6. Theory and the comparison: ratio, difference in combined sigma (z), its two-sided p,
    #    and the number of turns the measured slope implies if R and lambda are right.
    a_th, comps = theory()
    u_th = float(np.sqrt(sum(v ** 2 for v in comps.values())))
    ratio = a / a_th
    u_ratio = ratio * np.hypot(u_a / a, u_th / a_th)
    diff = a - a_th
    u_diff = np.hypot(u_a, u_th)
    z = diff / u_diff
    p_z = 2 * stats.norm.sf(abs(z))
    N_impl = ratio * N_TURNS
    u_N_impl = N_impl * np.hypot(u_a / a, 2 * U_RADIUS / RADIUS)

    Pglob, rings = ring_by_ring(fr, omega)

    # 7. The other phase methods: same averaging and fit, and frame by frame against the minima
    #    (they read the same photographs, so the shared noise cancels in the difference).
    #    Per method: its slope with fit and calibration errors (am, usm, ucm), and the slope of
    #    (method - minima) per photograph against omega (ad, u_d), which tests whether the two
    #    methods really disagree.
    mfile = HERE / 'methods_frames.csv'
    if mfile.exists():
        fr = fr.merge(pd.read_csv(mfile), on='file', how='left', validate='one_to_one')
    others = {}
    for label, col, _ in OTHER_METHODS:
        if col not in fr:
            continue
        t = setting_means(fr, col)
        ym, semm = t.phi.to_numpy(), t['sem'].to_numpy()
        am, bm, covm, _, sxm, _, ufm, chibm = fit_line(omega, ym, semm, eps)
        usm = float(np.sqrt(covm[0, 0]))
        ucm = calibration_part(S, ym, ufm, model)[3]
        uam = float(np.hypot(usm, ucm))
        td = setting_means(fr.assign(d=wrap(fr[col] - fr.phi_rad)), 'd')
        ad, _, covd, _, sxd, _, _, chibd = fit_line(omega, td.phi.to_numpy(), td['sem'].to_numpy(), 0.0)
        u_dset = np.hypot(td['sem'].to_numpy(), td['sem'].iloc[0]); u_dset[0] = 0.0
        others[label] = dict(col=col, y=ym, sem=semm, u_fit=ufm, a=am, b=bm, u_stat=usm, u_cal=ucm,
                             u_a=uam, sigma_x=sxm, chi_base=chibm, ratio=am / a_th,
                             u_ratio=am / a_th * np.hypot(uam / am, u_th / a_th),
                             z=(am - a_th) / np.hypot(uam, u_th), d_slope=ad,
                             u_d=float(np.sqrt(covd[0, 0])), chi_d=chibd, sx_d=sxd,
                             d_setting=td.phi.to_numpy() - td.phi.iloc[0], u_d_setting=u_dset)
    spread = [a] + [others[l]['a'] for l, _, use in OTHER_METHODS if use and l in others]
    u_method = (max(spread) - min(spread)) / 2              # half the range of the reliable methods

    # 8. Drift test: a phase linear in time fitted together with the slope. Because the settings
    #    were taken in increasing order, omega and time rise together; rho_t near -1 means the
    #    data cannot tell a slope in omega from a drift in time.
    tm = (fr.groupby('S').time_s.mean() - fr.time_s.min()).to_numpy()
    pt, covt, _ = wls(np.c_[omega, tm, np.ones_like(omega)], y, u_fit)
    rho_t = covt[0, 1] / np.sqrt(covt[0, 0] * covt[1, 1])
    increasing = np.all(np.diff(tab.t_start.to_numpy()) > 0)

    # 9. Outputs: the per-setting table, then the text report (sagnac.txt) line by line.
    tab['omega'] = omega; tab['u_omega_curve'] = u_curve; tab['u_omega_repro'] = u_rep
    tab['dphi'] = y - y[0]; tab['u_fit'] = u_fit
    tab['residual'] = y - (a * omega + b)
    for m in others.values():
        key = m['col'].removeprefix('phi_').removesuffix('_rad')
        tab[f'phi_{key}'] = m['y']; tab[f'sem_{key}'] = m['sem']; tab[f'u_fit_{key}'] = m['u_fit']
    tab.to_csv(HERE / 'sagnac_speeds.csv', index=False, float_format='%.6g')

    L = ['Sagnac phase against angular velocity (sagnac.py), fibre loop centred on the rotation axis\n']
    L.append('Per setting: phase at the ring centre, mean over the frames, u = SD/sqrt(n) '
             '(statistical, frame-to-frame); omega from calibration.json +- curve (correlated) '
             '+- speed reproducibility (independent).')
    L.append(f"{'S':>6} {'n':>3} {'omega (rad/s)':>26} {'phi (rad)':>9} {'SD':>6} {'u=SEM':>6} "
             f"{'dphi (rad)':>10} {'u_fit':>6} {'P (k px^2)':>10} {'drift (mrad/s)':>15} {'p':>5}")
    for r, w, uc, ur, uf in zip(tab.itertuples(), omega, u_curve, u_rep, u_fit):
        L.append(f'{r.S:6.2f} {r.n:3d} {w:8.4f} +- {uc:.4f} +- {ur:.4f} {r.phi:9.3f} {r.sd:6.3f} '
                 f'{r.sem:6.3f} {r.phi - y[0]:10.3f} {uf:6.3f} {r.P / 1e3:10.1f} '
                 f'{1e3 * r.drift:7.2f} +- {1e3 * r.u_drift:4.2f} {r.p_drift:5.2f}')
    L.append('(u_fit = sqrt(SEM^2 + (a eps omega)^2 + sigma_x^2), the error used in the fit; drift: '
             'slope of the frame phases against time within the setting, with its standard error and '
             'p for zero slope)\n')
    L.append(f'Line with the frame SEMs and the speed reproducibility (eps = {100 * eps:.2f} %) only: '
             f'chi2 = {chi_base:.0f}, dof = {dof}, chi2/dof = {chi_base / dof:.1f}, p = {p_base:.1g} '
             f'(frame SEMs alone: chi2 = {p_sem_only[2]:.0f}).')
    if sx > 0:
        L.append('  -> rejected: the settings scatter about the line more than their frames and the '
                 f'speed reproducibility explain. Extra scatter sigma_x = {sx:.3f} rad per setting '
                 'added in quadrature (chi2 = dof).')
    L.append(f'  (frame-to-frame SEMs {sem.min():.3f}-{sem.max():.3f} rad, speed reproducibility '
             f'{np.min((a * u_rep)[S > 0]):.3f}-{np.max(a * u_rep):.3f} rad; settings with a significant '
             f'drift within their 20-60 s: {int(np.sum(tab.p_drift < 0.05))} of {len(tab)}.)')
    L.append(f'Slope: a = {a:.4f} +- {u_stat:.4f} (fit) +- {u_cal:.4f} (calibration curve) '
             f'= {a:.4f} +- {u_a:.4f} rad/(rad/s)')
    L.append(f'  intercept b = {b:.3f} +- {np.sqrt(cov[1, 1]):.3f} rad (phase at rest, arbitrary); '
             f'chi2 = {chi:.1f}, dof = {dof}' + (' by construction.' if sx > 0
                                                 else f', p = {stats.chi2.sf(chi, dof):.2f}.'))
    L.append(f'  with the frame SEMs alone the slope would be {p_sem_only[0][0]:.4f} +- '
             f'{np.sqrt(p_sem_only[1][0, 0]):.4f} (fit).')
    L.append(f'Curvature check (phi = q omega^2 + a omega + b, same errors): q = {q[0]:.3f} +- '
             f'{np.sqrt(covq[0, 0]):.3f} rad/(rad/s)^2, {abs(zq):.1f} sigma from zero.\n')
    L.append(f"Theory, Eq. (5): N = {N_TURNS} +- {U_N_TURNS}, R = {100 * RADIUS:.1f} +- "
             f"{100 * U_RADIUS:.1f} cm (A = pi R^2 = {np.pi * RADIUS ** 2:.4f} m^2), "
             f'lambda = {WAVELENGTH * 1e9:.0f} nm, c exact.')
    L.append(f'  a_th = 8 pi^2 N R^2 / (lambda c) = {a_th:.4f} +- {u_th:.4f} rad/(rad/s)  '
             f"[N: {comps['N']:.4f}, R: {comps['R']:.4f}"
             + (f", lambda: {comps['lambda']:.4f}]" if 'lambda' in comps else ']'))
    if U_WAVELENGTH is None:
        L.append('  The lab gave no wavelength uncertainty; it is NOT included (unknown, not zero). '
                 f'Every 1 nm of error in lambda moves a_th by {a_th * 1e-9 / WAVELENGTH:.4f} '
                 f'({100e-9 / WAVELENGTH:.2f} %).')
    L.append(f'\nMeasured / theory: a / a_th = {ratio:.3f} +- {u_ratio:.3f}')
    L.append(f'a - a_th = {diff:+.4f} +- {u_diff:.4f} rad/(rad/s) = {z:+.2f} combined sigma, '
             f'p = {p_z:.2f} (two-sided, normal).')
    L.append(f'Equivalently, the measured slope corresponds to N = {N_impl:.2f} +- {u_N_impl:.2f} turns '
             f"(given R and lambda), against N = {N_TURNS} +- {U_N_TURNS}.")
    u_dmax = np.hypot(u_fit[-1], u_fit[0])
    L.append(f'Largest shift: dphi = {y[-1] - y[0]:.3f} +- {u_dmax:.3f} rad (fit errors of both '
             f'settings) = {(y[-1] - y[0]) / (2 * np.pi):.3f} +- {u_dmax / (2 * np.pi):.3f} fringe, '
             f'at omega = {omega[-1]:.4f} +- {u_omega[-1]:.4f} rad/s.\n')
    L.append('Report form (one significant figure, comparisons above from unrounded values):')
    L.append(f'  a = {round1(a, u_a)} rad/(rad/s) (fit and calibration combined); '
             f'a_th = {round1(a_th, u_th)} rad/(rad/s); a/a_th = {round1(ratio, u_ratio)}\n')
    L.append('Other phase methods, same averaging and fit (cosine fit: fringes.py, cosine_fit, prep Q5; '
             'template matching, 2D Fourier and AIA: methods.py, robustness scans in methods.txt):')
    for label, m in others.items():
        L.append(f'{label}: dphi per setting ' + ', '.join(f'{d:.3f}' for d in m['y'] - m['y'][0])
                 + f' rad; SEM {m["sem"].min():.3f}-{m["sem"].max():.3f} rad.')
        if label == 'cosine fit':
            rms_c = fr.groupby('S').rms_cos.median().to_numpy()
            L.append(f'  fit rms (median per setting): {", ".join(f"{v:.1f}" for v in rms_c)} counts '
                     '(no chi2: neighbouring pixels share speckle, so the bin errors are not independent).')
        L.append(f'  line with SEM and speed reproducibility: chi2 = {m["chi_base"]:.0f}, dof = {dof}, '
                 f'p = {stats.chi2.sf(m["chi_base"], dof):.1g}; sigma_x = {m["sigma_x"]:.3f} rad '
                 f'(minima: {sx:.3f}).')
        L.append(f'  a = {m["a"]:.4f} +- {m["u_stat"]:.4f} (fit) +- {m["u_cal"]:.4f} (calibration curve) '
                 f'= {m["a"]:.4f} +- {m["u_a"]:.4f} rad/(rad/s); a / a_th = {m["ratio"]:.3f} +- '
                 f'{m["u_ratio"]:.3f}, {m["z"]:+.2f} combined sigma from theory, '
                 f'p = {2 * stats.norm.sf(abs(m["z"])):.2f}.')
        L.append(f'  minus minima, frame by frame: slope {m["d_slope"]:+.4f} +- {m["u_d"]:.4f} rad/(rad/s) '
                 f'({abs(m["d_slope"]) / m["u_d"]:.1f} sigma; chi2 without extra scatter {m["chi_d"]:.0f}, '
                 f'dof = {dof}, p = {stats.chi2.sf(m["chi_d"], dof):.1g}, sigma_x = {m["sx_d"]:.3f} rad). '
                 'Per setting, relative to rest: ' + ', '.join(f'{d:+.2f}' for d in m['d_setting']) + ' rad.')
        L.append(f'  Report form: a = {round1(m["a"], m["u_a"])} rad/(rad/s); '
                 f'a/a_th = {round1(m["ratio"], m["u_ratio"])}')
    if len(spread) > 1:
        used = ['minima'] + [l for l, _, use in OTHER_METHODS if use and l in others]
        z_m = diff / np.sqrt(u_a ** 2 + u_method ** 2 + u_th ** 2)
        L.append(f'Spread of the slopes of {", ".join(used)}: {min(spread):.3f}-{max(spread):.3f}; half the '
                 f'range, {u_method:.3f}, is the phase-extraction (method) systematic. AIA is left out: '
                 'its slope depends on its window and variant (methods.txt).')
        L.append(f'  a = {a:.4f} +- {u_a:.4f} (fit and calibration) +- {u_method:.4f} (method) rad/(rad/s); '
                 f'with it, a - a_th = {z_m:+.2f} combined sigma, p = {2 * stats.norm.sf(abs(z_m)):.2f}. '
                 f'Report form: a = {round1(a, np.hypot(u_a, u_method))} rad/(rad/s) (all three parts).\n')
    # Speed of light: Eq. (5) solved for c with the measured slope (minima). The slope's three
    # error parts and N, R are kept as separate components, then combined in quadrature.
    c_parts = {'fit': u_stat, 'calibration curve': u_cal}
    if len(spread) > 1:
        c_parts['method'] = u_method
    c_m, c_comps = speed_of_light(a, c_parts)
    u_c = float(np.sqrt(sum(v ** 2 for v in c_comps.values())))
    u_c_meas = float(np.sqrt(sum(c_comps[k] ** 2 for k in c_parts)))
    u_c_noN = float(np.sqrt(u_c ** 2 - c_comps['N'] ** 2))
    z_c = (c_m - C_LIGHT) / u_c
    p_c = 2 * stats.norm.sf(abs(z_c))
    c_methods = {}
    for label, _, use in OTHER_METHODS:
        if use and label in others:
            m = others[label]
            cm, cc = speed_of_light(m['a'], {'fit': m['u_stat'], 'calibration curve': m['u_cal']})
            c_methods[label] = dict(c=cm, u=float(np.sqrt(sum(v ** 2 for v in cc.values()))))
    E8 = 1e8
    L.append('Speed of light, Eq. (5) solved for c: c = 8 pi^2 N R^2 / (lambda a), with the measured slope a '
             '(minima) and the given N, R, lambda. In units of 10^8 m/s:')
    L.append(f'  c = {c_m / E8:.4f}' + ''.join(f' +- {v / E8:.4f} ({k})' for k, v in c_comps.items())
             + f'\n    = {c_m / E8:.4f} +- {u_c / E8:.4f} ({100 * u_c / c_m:.1f} %, all parts)')
    L.append(f'  The measurement ({", ".join(c_parts)}) contributes {u_c_meas / E8:.4f} '
             f'({100 * u_c_meas / c_m:.1f} %), N = {N_TURNS} +- {U_N_TURNS} contributes '
             f'{c_comps["N"] / E8:.4f} ({100 * c_comps["N"] / c_m:.1f} %); without the N term the total '
             f'is {u_c_noN / E8:.4f} ({100 * u_c_noN / c_m:.1f} %).')
    if U_WAVELENGTH is None:
        L.append('  lambda: no uncertainty given, NOT included (unknown, not zero). Every 1 nm of error in '
                 f'lambda moves c by {c_m * 1e-9 / WAVELENGTH / E8:.4f} ({100e-9 / WAVELENGTH:.2f} %).')
    L.append(f'  Against the defined c = {C_LIGHT:.0f} m/s: c / c_def = {c_m / C_LIGHT:.3f} +- '
             f'{u_c / C_LIGHT:.3f}; c - c_def = {(c_m - C_LIGHT) / E8:+.4f} +- {u_c / E8:.4f} = '
             f'{z_c:+.2f} combined sigma, p = {p_c:.2f} (two-sided, normal).')
    L.append('  This is the comparison a / a_th above read the other way (c / c_def = a_th / a), here with '
             'the method systematic included.')
    if c_methods:
        L.append('  Other methods (fit, calibration curve, N and R; no method term): '
                 + ', '.join(f'{l} {v["c"] / E8:.3f} +- {v["u"] / E8:.3f}' for l, v in c_methods.items())
                 + '.')
    L.append(f'  Report form: c = ({round1(c_m / E8, u_c / E8)}) x 10^8 m/s (all parts).')
    L.append('  The phase fixes the product lambda c (= c^2 / nu). With the given wavelength c follows as above; '
             'if the laser frequency nu were known instead, c = sqrt(8 pi^2 N R^2 nu / a), with half the '
             'relative uncertainty.\n')
    L.append(f'Cross-check, Eq. (7) ring by ring with one period P = {Pglob / 1e3:.1f}k px^2 for all '
             f'frames (no extrapolation to the centre):')
    for k, v in rings.items():
        L.append(f"  {k}: slope {v['slope']:.3f} rad/(rad/s); dphi = "
                 + ', '.join(f'{d:.2f}' for d in v['dphi']))
    if len(rings) >= 2:
        d1, d2 = (np.array(v['dphi']) for v in list(rings.values())[:2])
        Pn = tab.P.to_numpy()
        L.append(f'  The two rings differ by up to {np.nanmax(np.abs(d1 - d2)):.2f} rad at single '
                 'settings, more than their frame scatter; the fringe period (curvature difference '
                 f'of the two beams) ranges over {Pn.min() / 1e3:.1f}-{Pn.max() / 1e3:.1f}k px^2 '
                 f'({100 * (Pn.max() - Pn.min()) / 2 / np.median(Pn):.1f} % half-range) between '
                 'settings, which the per-frame centre phase accounts for.\n')
    L.append('Systematics not in the numbers above:')
    t0, t1 = fr.time_s.min(), fr.time_s.max()
    if increasing:
        L.append(f'  * The settings were taken in increasing order ({hhmm(t0)}-{hhmm(t1)}) with the '
                 'only rest frames at the start, so a phase drift linear in time is indistinguishable '
                 'from the Sagnac slope. A rest frame at the end (or the reverse rotation) would test it.')
    L.append(f'    Fitting phi = a omega + D t + b: a = {pt[0]:.2f} +- {np.sqrt(covt[0, 0]):.2f} '
             f'rad/(rad/s), D = {1e3 * pt[1]:.1f} +- {1e3 * np.sqrt(covt[1, 1]):.1f} mrad/s, '
             f'correlation {rho_t:+.3f}. The within-setting drifts range over '
             f'{1e3 * tab.drift.min():+.1f} to {1e3 * tab.drift.max():+.1f} mrad/s and the residual '
             f'scatter about the line is sigma_x = {sx:.3f} rad, so the phase wanders but does not '
             'random-walk; a steady drift is not excluded by this.')
    L.append('  * omega is the calibration at the same setting, not measured during the photographs '
             '(the reproducibility eps covers the scatter seen in the calibration).')
    L.append('  * The sign of the shift (rings shrink) is not compared: the rotation direction is '
             'not recorded here.')
    (HERE / 'sagnac.txt').write_text('\n'.join(L) + '\n', encoding='utf-8')
    print('\n'.join(L))

    # Machine-readable results, then the two figures.
    out = dict(slope=a, u_stat=u_stat, u_cal=u_cal, u_slope=u_a, intercept=b,
               sigma_x=sx, eps=eps, chi2_without_sigma_x=chi_base, p_without_sigma_x=p_base,
               dof=dof, curvature_q=q[0], u_q=float(np.sqrt(covq[0, 0])),
               theory=a_th, u_theory=u_th, theory_components=comps,
               u_wavelength='not given', ratio=ratio, u_ratio=u_ratio, z=z, p_z=p_z,
               N_implied=N_impl, u_N_implied=u_N_impl,
               ring_by_ring={k: v['slope'] for k, v in rings.items()}, P_global=Pglob,
               u_method=u_method,
               methods={l: dict(slope=m['a'], u_stat=m['u_stat'], u_cal=m['u_cal'], u_slope=m['u_a'],
                                sigma_x=m['sigma_x'], ratio=m['ratio'], u_ratio=m['u_ratio'], z=m['z'],
                                minus_minima_slope=m['d_slope'], u_minus_minima=m['u_d'])
                        for l, m in others.items()})
    out['speed_of_light'] = dict(c=c_m, u_c=u_c, components=c_comps, u_measurement=u_c_meas,
                                 u_without_N=u_c_noN, ratio_to_defined=c_m / C_LIGHT,
                                 u_ratio_to_defined=u_c / C_LIGHT, z=z_c, p_z=p_c,
                                 u_wavelength='not given', methods=c_methods)
    (HERE / 'sagnac.json').write_text(json.dumps(out, indent=1, default=float))
    c = others['cosine fit']
    make_figure(tab, omega, u_omega, u_fit, a, b, a_th, u_th, HERE / 'figures' / 'phase_vs_omega',
                (c['y'], c['u_fit'], c['a'], c['b']))
    scan = HERE / 'methods.json'
    make_comparison(omega, others, a, u_a, a_th, u_th, json.loads(scan.read_text()) if scan.exists() else {},
                    HERE / 'figures' / 'methods_comparison')


if __name__ == '__main__':
    main()
