"""Sagnac phase against angular velocity, compared with the manual's Eq. (5).

Inputs: the frame table from fringes.py and calibration.json (calibration.py).
  python sagnac.py                 large loop, centred on the axis (fringes_frames.csv)
  python sagnac.py --loop small    small loop, off the axis (fringes_frames_small.csv), task 4

Per speed setting S: the frame phases phi_c are unwrapped about their circular mean; the
setting's phase is their mean, with u = SD/sqrt(n) (statistical, frame-to-frame scatter).
The setting means are unwrapped in order of S. omega(S) comes from the calibration
(S = 0 is the table at rest, omega = 0 exactly). Its uncertainty has two parts: the
calibration curve (correlated across settings) and the speed reproducibility eps omega at a
setting (independent), see calibration.py.

Fit phi(omega) = a omega + b (weighted least squares; b absorbs the arbitrary phase at rest).
The per-setting error is u^2 = SEM^2 + (a eps omega)^2 + sigma_x^2: the reproducibility enters
as an error in omega (effective variance, iterated on a). If chi2 without sigma_x rejects the
line, the speed-to-speed scatter exceeds what the frames and the speed reproducibility explain
(the phase wanders on time scales longer than the 20-60 s taken at each setting); sigma_x, the
same for every setting, is then set by chi2 = dof. The slope's uncertainty has two parts:
  u_stat  from the fit (phase scatter and speed reproducibility),
  u_cal   from the calibration curve (its covariance, propagated through the fit).
The second method of the prep (Q5), a cosine fit of each frame's profile (fringes.py), goes
through the same averaging and fit; the two are also compared frame by frame.
Theory, Eq. (5) with A = pi R^2:  a_th = 8 pi N A / (lambda c) = 8 pi^2 N R^2 / (lambda c),
independent of where the rotation axis lies (Schwartz 2017).

Outputs (suffix _small for the small loop): sagnac.txt, sagnac.json, sagnac_speeds.csv,
figures/phase_vs_omega.pdf/.png. Run from this folder after calibration.py and fringes.py.
"""
import argparse
import json
import time
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

from calibration import omega_with_cov

HERE = Path(__file__).resolve().parent

# Fibre loops (values given by the experimenters; +- taken as standard uncertainties)
LOOPS = {
    'big': dict(label='large loop, centred on the axis', N=13, u_N=1, R=0.570, u_R=0.005,
                frames='fringes_frames.csv', suffix=''),
    'small': dict(label='small loop, off the axis', N=46, u_N=1, R=0.155, u_R=0.001,
                  frames='fringes_frames_small.csv', suffix='_small'),
}
# Other phase methods (label, frame column, counted in the method spread). AIA is not counted:
# its slope depends on its window and variant here (methods.txt).
OTHER_METHODS = [('cosine fit', 'phi_cos_rad', True), ('template matching', 'phi_template_rad', True),
                 ('2D Fourier', 'phi_fourier_rad', True), ('AIA', 'phi_aia_rad', False)]
WAVELENGTH = 640e-9                   # m
U_WAVELENGTH = None                   # not given by the lab: left out of the budget, stated in the output
C_LIGHT = 299_792_458.0               # m/s, exact


def wrap(x):
    return np.angle(np.exp(1j * np.asarray(x)))


def round1(x, u):
    """'x +- u' with u to one significant figure and x to the same digit, rounding half up on the
    decimal digits as printed (float round() sends ties to even: 2.25 -> 2.2, and 1.005 -> 1.0)."""
    if not np.isfinite(u) or u <= 0:
        return f'{x:.4g} +- {u}'
    U = Decimal(str(float(u)))
    q = Decimal(1).scaleb(U.adjusted())
    Ur = U.quantize(q, ROUND_HALF_UP)
    if Ur.adjusted() > U.adjusted():                     # 0.096 -> 0.1 moves the digit
        q = q.scaleb(1)
        Ur = U.quantize(q, ROUND_HALF_UP)
    X = Decimal(str(float(x))).quantize(q, ROUND_HALF_UP)
    return f'{abs(X) if X == 0 else X:f} +- {Ur:f}'


def theory(loop):
    """a_th (rad/(rad/s)) and its uncertainty components."""
    a = 8 * np.pi ** 2 * loop['N'] * loop['R'] ** 2 / (WAVELENGTH * C_LIGHT)
    comps = {'N': a * loop['u_N'] / loop['N'], 'R': a * 2 * loop['u_R'] / loop['R']}
    if U_WAVELENGTH is not None:
        comps['lambda'] = a * U_WAVELENGTH / WAVELENGTH
    return a, comps


def setting_means(fr, col='phi_rad'):
    """One row per setting: n, mean phase, SD, SEM, time span and within-setting drift."""
    rows = []
    for S, g in fr.groupby('S', sort=True):
        ph = g[col].to_numpy()
        c = np.angle(np.mean(np.exp(1j * ph)))
        ph = c + wrap(ph - c)
        t = g.time_s.to_numpy() - g.time_s.min()
        n = len(ph)
        lr = stats.linregress(t, ph)
        rows.append(dict(S=S, n=n, phi=ph.mean(), sd=ph.std(ddof=1), sem=ph.std(ddof=1) / np.sqrt(n),
                         t_start=g.time_s.min(), span_s=t.max(), drift=lr.slope,
                         u_drift=lr.stderr, p_drift=lr.pvalue, P=g.P_px2.median()))
    tab = pd.DataFrame(rows)
    phi = tab.phi.to_numpy().copy()
    for i in range(1, len(phi)):                       # unwrap in order of S
        phi[i] += 2 * np.pi * np.round((phi[i - 1] - phi[i]) / (2 * np.pi))
    tab['phi'] = phi
    return tab


def wls(A, y, u):
    """Weighted linear least squares y = A p. Returns p, covariance, chi2."""
    W = 1 / u ** 2
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    p = cov @ (A.T @ (W * y))
    return p, cov, float(np.sum(((y - A @ p) / u) ** 2))


def wls_line(x, y, u):
    """Weighted y = a x + b. Returns (a, b), covariance, chi2."""
    return wls(np.c_[x, np.ones_like(x)], y, u)


def extra_scatter(x, y, u_base):
    """sigma_x >= 0 with chi2(sqrt(u_base^2 + sigma_x^2)) = dof (0 if chi2 <= dof already)."""
    dof = len(x) - 2
    f = lambda s: wls_line(x, y, np.sqrt(u_base ** 2 + s ** 2))[2] - dof
    return 0.0 if f(0.0) <= 0 else optimize.brentq(f, 0.0, 10.0)


def fit_line(omega, y, sem, eps):
    """Line with SEM, the speed reproducibility (a eps omega) and sigma_x; iterated on a."""
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
    """omega per setting, its curve and reproducibility parts, and the slope's u_cal."""
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
    J = np.array([(slope(coef + h) - slope(coef - h)) / (2 * h[j])
                  for j, h in enumerate(np.diag(np.abs(coef) * 1e-6 + 1e-12))])
    return omega, u_curve, u_rep, float(np.sqrt(J @ Cc @ J))


def ring_by_ring(fr, omega):
    """Cross-check: manual Eq. (7) per dark ring with one global period (no centre extrapolation)."""
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
    return time.strftime('%H:%M', time.gmtime(t))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--loop', choices=LOOPS, default='big')
    args = ap.parse_args()
    loop = LOOPS[args.loop]
    sfx = loop['suffix']
    src = HERE / loop['frames']
    if not src.exists():
        raise SystemExit(f'{src.name} not found: run fringes.py --data <folder> --out {src.name} '
                         f'on the photographs of the {loop["label"]} first (see README.md).')

    fr = pd.read_csv(src)
    model = json.loads((HERE / 'calibration.json').read_text())['conversion']
    eps = model['eps']
    tab = setting_means(fr)
    S = tab.S.to_numpy()
    sem = tab['sem'].to_numpy()
    y = tab.phi.to_numpy()
    dof = len(y) - 2

    omega, u_curve, u_rep, _ = calibration_part(S, y, None, model)
    a, b, cov, chi, sx, u_base, u_fit, chi_base = fit_line(omega, y, sem, eps)
    p_base = stats.chi2.sf(chi_base, dof)
    p_sem_only = wls_line(omega, y, sem)
    u_stat = float(np.sqrt(cov[0, 0]))
    _, _, _, u_cal = calibration_part(S, y, u_fit, model)
    u_a = float(np.hypot(u_stat, u_cal))
    u_omega = np.hypot(u_curve, u_rep)

    # curvature check: a quadratic term with the same errors
    q, covq, _ = wls(np.c_[omega ** 2, omega, np.ones_like(omega)], y, u_fit)
    zq = q[0] / np.sqrt(covq[0, 0])

    a_th, comps = theory(loop)
    u_th = float(np.sqrt(sum(v ** 2 for v in comps.values())))
    ratio = a / a_th
    u_ratio = ratio * np.hypot(u_a / a, u_th / a_th)
    diff = a - a_th
    u_diff = np.hypot(u_a, u_th)
    z = diff / u_diff
    p_z = 2 * stats.norm.sf(abs(z))
    N_impl = ratio * loop['N']
    u_N_impl = N_impl * np.hypot(u_a / a, 2 * loop['u_R'] / loop['R'])

    Pglob, rings = ring_by_ring(fr, omega)

    # the other phase methods: same averaging and fit, and frame by frame against the minima
    # (they read the same photographs, so the shared noise cancels in the difference)
    mfile = HERE / f'methods_frames{sfx}.csv'
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

    # drift test: a phase linear in time fitted together with the slope
    tm = (fr.groupby('S').time_s.mean() - fr.time_s.min()).to_numpy()
    pt, covt, _ = wls(np.c_[omega, tm, np.ones_like(omega)], y, u_fit)
    rho_t = covt[0, 1] / np.sqrt(covt[0, 0] * covt[1, 1])
    increasing = np.all(np.diff(tab.t_start.to_numpy()) > 0)

    tab['omega'] = omega; tab['u_omega_curve'] = u_curve; tab['u_omega_repro'] = u_rep
    tab['dphi'] = y - y[0]; tab['u_fit'] = u_fit
    tab['residual'] = y - (a * omega + b)
    for m in others.values():
        key = m['col'].removeprefix('phi_').removesuffix('_rad')
        tab[f'phi_{key}'] = m['y']; tab[f'sem_{key}'] = m['sem']; tab[f'u_fit_{key}'] = m['u_fit']
    tab.to_csv(HERE / f'sagnac_speeds{sfx}.csv', index=False, float_format='%.6g')

    L = [f'Sagnac phase against angular velocity (sagnac.py), {loop["label"]}\n']
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
    L.append(f"Theory, Eq. (5): N = {loop['N']} +- {loop['u_N']}, R = {100 * loop['R']:.1f} +- "
             f"{100 * loop['u_R']:.1f} cm (A = pi R^2 = {np.pi * loop['R'] ** 2:.4f} m^2), "
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
             f"(given R and lambda), against N = {loop['N']} +- {loop['u_N']}.")
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
    (HERE / f'sagnac{sfx}.txt').write_text('\n'.join(L) + '\n', encoding='utf-8')
    print('\n'.join(L))

    out = dict(loop=args.loop, slope=a, u_stat=u_stat, u_cal=u_cal, u_slope=u_a, intercept=b,
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
    (HERE / f'sagnac{sfx}.json').write_text(json.dumps(out, indent=1, default=float))
    c = others['cosine fit']
    make_figure(tab, omega, u_omega, u_fit, a, b, a_th, u_th, HERE / 'figures' / f'phase_vs_omega{sfx}',
                (c['y'], c['u_fit'], c['a'], c['b']))
    scan = HERE / f'methods{sfx}.json'
    make_comparison(omega, others, a, u_a, a_th, u_th, json.loads(scan.read_text()) if scan.exists() else {},
                    HERE / 'figures' / f'methods_comparison{sfx}')


if __name__ == '__main__':
    main()
