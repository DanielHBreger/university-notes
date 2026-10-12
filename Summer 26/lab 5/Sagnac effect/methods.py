"""Two more ways to read the fringe phase, to compare with fringes.py (minima, cosine fit).

Input: fringes_frames.csv and fringes_profiles.npz (fringes.py): for every photograph, the
intensity averaged over the fringe sector in bins of s = r^2 about its ring centre.

  Template matching. The rest frames (S = 0) are averaged into a template T(s). Each profile
  is matched to g T(s + D) + c by least squares over the shift D, with gain g and offset c free.
  In s, an exponential envelope shifted by D only changes by a constant factor, which g absorbs.
  Rings that move inwards give D > 0, and the phase relative to rest is 2 pi D / P_T, where P_T
  is the shift that maps the template onto itself one fringe further out. No fringe model.

  Phase shifting with unknown steps (AIA, Wang & Han, Opt. Lett. 29, 1671 (2004)). Bin j of
  frame k is I_kj = a_j + b_j cos(psi_j + delta_k): every bin brightens and dims with the frame's
  phase delta_k. Alternating least squares: given delta, each bin's a, b cos psi, b sin psi;
  given those, each frame's cos delta, sin delta. Started from the first two principal
  components of the profiles. No fringe period, envelope or cosine shape in r is assumed;
  psi(s) comes out of the fit and should be linear in s (rings equally spaced in r^2). delta is
  defined up to an offset and a sign; the sign is taken so that delta grows with S.

Both give an average phase over the window rather than the phase at the ring centre, so a change
of the fringe period between frames (fringes.py) enters them, but not the minima or the cosine fit.
Each is also run with other comparison ranges (template) or windows and a per-frame brightness
factor (AIA), to show how much its slope depends on those choices.

Outputs: methods_frames.csv (phi_template_rad, phi_aia_rad per frame), methods.txt, methods.json
(the robustness scans), figures/methods_illustration.pdf/.png (figure*).
Run from this folder after calibration.py and fringes.py, before sagnac.py.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import interpolate, ndimage, optimize

from fringes import P_GUESS, SMOOTH_PX, cos_model, cosine_fit, load_red
from sagnac import calibration_part, fit_line, setting_means

HERE = Path(__file__).resolve().parent
MATCH_S = (2e4, 1.8e5)          # px^2, template-matching comparison range
SHIFT_RANGE = (-2e4, 7.8e4)     # px^2, allowed shifts: s + D stays inside the profile
SELF_S = (1e3, 1.1e5)           # px^2, comparison range for the template's own period
AIA_TOL = 1e-10                 # rad, convergence of the frame phases
MATCH_SCAN = [(2e4, 1.8e5), (1e4, 1.2e5), (2e4, 1.4e5), (5e4, 1.8e5)]   # px^2, template robustness
AIA_SCAN = [2.6e5, 1.8e5, 1.3e5, 0.9e5]                                # px^2, largest s used by AIA


def wrap(x):
    return np.angle(np.exp(1j * np.asarray(x)))


def match_residual(D, s, prof, T):
    """Least-squares residual of prof(s) against g T(s + D) + c (g, c solved), and (g, c)."""
    A = np.c_[T(s + D), np.ones_like(s)]
    gc, *_ = np.linalg.lstsq(A, prof, rcond=None)
    return float(np.sum((prof - A @ gc) ** 2)), gc


def best_shift(s, prof, T, lo, hi):
    """Shift D in [lo, hi] minimising match_residual: grid of 1000 px^2, then Brent."""
    grid = np.arange(lo, hi + 1, 1000.0)
    i = int(np.argmin([match_residual(D, s, prof, T)[0] for D in grid]))
    a, b = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    return optimize.minimize_scalar(lambda D: match_residual(D, s, prof, T)[0],
                                    bounds=(a, b), method='bounded').x


def template_phases(s, M, rest, match=MATCH_S):
    """Template T, its period P_T, and the shift D and phase 2 pi D / P_T of every profile,
    compared over match (px^2)."""
    T = interpolate.CubicSpline(s, M[rest].mean(axis=0))
    sel = (s >= SELF_S[0]) & (s <= SELF_S[1])
    P_T = best_shift(s[sel], T(s[sel]), T, 0.8 * P_GUESS, 1.25 * P_GUESS)
    sel = (s >= match[0]) & (s <= match[1])
    D = np.array([best_shift(s[sel], p[sel], T, *SHIFT_RANGE) for p in M])
    return T, P_T, D, 2 * np.pi * D / P_T


def aia(M, gain=False):
    """Frame phases delta, per-bin a, c = b cos psi, d = b sin psi, and the iteration count.
    gain=True also fits a brightness factor per frame: I_kj = g_k (a_j + c_j cos delta_k + ...)."""
    U, S, _ = np.linalg.svd(M - M.mean(axis=0), full_matrices=False)
    delta = np.arctan2(U[:, 1] * S[1], U[:, 0] * S[0])                 # PCA start
    g = np.ones(len(M))
    for it in range(1, 3001):
        D = g[:, None] * np.c_[np.ones_like(delta), np.cos(delta), np.sin(delta)]
        a, c, d = np.linalg.lstsq(D, M, rcond=None)[0]                  # per bin
        if gain:                                                        # per frame: g, g cos, g sin
            gk, u, v = np.linalg.lstsq(np.c_[a, c, d], M.T, rcond=None)[0]
            g = gk / gk.mean()
        else:                                                           # per frame: cos, sin
            u, v = np.linalg.lstsq(np.c_[c, d], (M - a).T, rcond=None)[0]
        new = np.arctan2(v, u)
        done = np.max(np.abs(wrap(new - delta))) < AIA_TOL
        delta = new
        if done:
            break
    return delta, a, c, d, it


def aia_phases(fr, M, gain=False):
    """AIA frame phases relative to the rest frames, unwrapped along increasing S, sign so that
    they grow with S; also delta, a, c, d (same sign convention) and the iteration count."""
    delta, a, c, d, n_it = aia(M, gain)
    order = np.lexsort((fr.time_s, fr.S))
    du = np.empty_like(delta); du[order] = np.unwrap(delta[order])
    if np.polyfit(fr.S, du, 1)[0] < 0:
        delta, d, du = -delta, -d, -du
    rest = (fr.S == 0).to_numpy()
    return du - du[rest].mean(), delta, a, c, d, n_it


def slope(fr, phi, model):
    """Sagnac slope, its fit error and sigma_x for per-frame phases phi, as sagnac.py fits them."""
    tab = setting_means(fr.assign(x=phi), 'x')
    S, y, sem = tab.S.to_numpy(), tab.phi.to_numpy(), tab['sem'].to_numpy()
    omega = calibration_part(S, y, None, model)[0]
    a, _, cov, _, sx, *_ = fit_line(omega, y, sem, model['eps'])
    return a, float(np.sqrt(cov[0, 0])), sx


def illustration(fr, s, M, T, D, phi_aia, delta, a, c, d, path):
    """figure*: (a) cosine fit, (b) template match, (c) AIA bins against the frame phase, (d) psi(s)."""
    import matplotlib.pyplot as plt
    import report_style as rs
    rs.apply()
    tex = plt.rcParams['text.usetex']
    w = rs.TEXT_WIDTH
    fig, axs = plt.subplots(2, 2, figsize=(w, 0.62 * w), gridspec_kw=dict(hspace=0.32, wspace=0.22))
    top = np.flatnonzero(fr.S == fr.S.max())
    k = top[np.argsort(D[top])[len(top) // 2]]             # median-shift frame at the top speed
    x = s / 1e5
    xlab = r'$r^2$ ($10^5$ px$^2$)'

    # (a) cosine fit of that frame
    ax = axs[0, 0]
    row = fr.iloc[k]
    img = ndimage.gaussian_filter(load_red(HERE / row.file), SMOOTH_PX)
    popt, _, prof = cosine_fit(img, row.cx, row.cy)
    ax.plot(x, prof, 'o', color='k', ms=1.6, mew=0, label='profile')
    ax.plot(x, cos_model(s, *popt), '-', color=rs.COLORS[2], label='cosine fit')
    ax.set_xlabel(xlab); ax.set_ylabel('intensity (counts)'); ax.legend(loc='upper right')

    # (b) template matching: rest template, the frame, the shifted template
    ax = axs[0, 1]
    sel = (s >= MATCH_S[0]) & (s <= MATCH_S[1])
    _, (g, c0) = match_residual(D[k], s[sel], M[k][sel], T)
    ax.plot(x, T(s), '--', color='0.55', label='rest template')
    ax.plot(x, M[k], 'o', color='k', ms=1.6, mew=0, label='profile')
    ax.plot(x[sel], g * T(s[sel] + D[k]) + c0, '-', color=rs.COLORS[3],
            label=rf'template shifted by {D[k] / 1e3:.0f}k px$^2$' if tex
            else f'template shifted by {D[k] / 1e3:.0f}k px²')
    ax.axvspan(*(np.array(MATCH_S) / 1e5), color='0.92', zorder=0, lw=0)
    ax.set_ylim(0, 1.2 * T(s).max())
    ax.set_xlabel(xlab); ax.legend(loc='upper right')

    # (c) AIA: three bins oscillate with the frame phase (relative to rest)
    ax = axs[1, 0]
    off = np.angle(np.mean(np.exp(1j * (delta - phi_aia))))
    xx = np.linspace(phi_aia.min() - 0.2, phi_aia.max() + 0.2, 200)
    js = np.searchsorted(s, [3e4, 9e4, 1.5e5])
    for j, col, mk in zip(js, (rs.COLORS[0], rs.COLORS[4], rs.COLORS[5]), ('o', 's', '^')):
        ax.plot(phi_aia, M[:, j], mk, color=col, ms=2, mew=0, alpha=0.8,
                label=rf'$r^2 = {s[j] / 1e5:.1f}\times10^5$ px$^2$' if tex else f'r² = {s[j] / 1e5:.1f}e5 px²')
        ax.plot(xx, a[j] + c[j] * np.cos(xx + off) + d[j] * np.sin(xx + off), '-', color=col, lw=0.7)
    ax.set_xlabel('AIA frame phase, relative to rest (rad)'); ax.set_ylabel('intensity (counts)')
    ax.set_ylim(0, 1.1 * M[:, js].max())
    ax.legend(loc='upper right', fontsize=7)

    # (d) AIA: the fitted phase of every bin against r^2
    ax = axs[1, 1]
    psi = np.unwrap(np.arctan2(d, c))
    b = np.hypot(c, d)
    good = b > 0.2 * b.max()
    k1, k0 = np.polyfit(s[good], psi[good], 1, w=b[good])
    ax.plot(x[good], psi[good], 'o', color=rs.COLORS[3], ms=2, mew=0,
            label=r'$\psi$ from AIA' if tex else 'psi from AIA')
    ax.plot(x, k1 * s + k0, '-', color='0.4', lw=0.7,
            label=rf'linear, $P = {2 * np.pi / abs(k1) / 1e3:.1f}$k px$^2$' if tex
            else f'linear, P = {2 * np.pi / abs(k1) / 1e3:.1f}k px²')
    ax.set_xlabel(xlab); ax.set_ylabel(r'bin phase $\psi$ (rad)'); ax.legend(loc='upper right')
    for ax, lab in zip(axs.flat, 'abcd'):
        ax.text(-0.1, 1.02, f'({lab})', transform=ax.transAxes, va='bottom', ha='left')
    rs.save(fig, path)
    plt.close(fig)
    return 2 * np.pi / abs(k1), np.max(np.abs(psi[good] - (k1 * s[good] + k0)))


def main():
    fr = pd.read_csv(HERE / 'fringes_frames.csv')
    z = np.load(HERE / 'fringes_profiles.npz')
    s, M = z['s'], z['profiles']
    assert list(z['files']) == list(fr.file), 'profiles and frame table differ: rerun fringes.py'
    rest = (fr.S == 0).to_numpy()
    model = json.loads((HERE / 'calibration.json').read_text())['conversion']

    T, P_T, D, phi_t = template_phases(s, M, rest)
    phi_aia, delta, a, c, d, n_it = aia_phases(fr, M)
    pd.DataFrame(dict(file=fr.file, phi_template_rad=np.round(phi_t, 5),
                      shift_template_px2=np.round(D, 1), phi_aia_rad=np.round(phi_aia, 5))
                 ).to_csv(HERE / 'methods_frames.csv', index=False)
    P_aia, dev = illustration(fr, s, M, T, D, phi_aia, delta, a, c, d,
                              HERE / 'figures' / 'methods_illustration')

    # robustness: the same slope from other comparison ranges / windows / an AIA brightness factor
    t_scan = [(m, *slope(fr, template_phases(s, M, rest, m)[3], model)) for m in MATCH_SCAN]
    a_scan = [(smax, gain, *slope(fr, aia_phases(fr, M[:, s < smax], gain)[0], model))
              for gain in (False, True) for smax in AIA_SCAN]
    (HERE / 'methods.json').write_text(json.dumps(dict(
        template_scan=[dict(range=m, slope=k, u_fit=u, sigma_x=x) for m, k, u, x in t_scan],
        aia_scan=[dict(s_max=m, gain=g, slope=k, u_fit=u, sigma_x=x) for m, g, k, u, x in a_scan]),
        indent=1, default=float))

    resid = M - (a + np.outer(np.cos(delta), c) + np.outer(np.sin(delta), d))
    L = ['Two more phase methods (methods.py)\n',
         f'Template matching: template = mean of the {rest.sum()} rest profiles; its own period '
         f'P_T = {P_T / 1e3:.2f}k px^2 (minima median {fr.P_px2.median() / 1e3:.2f}k); shifts D = '
         f'{D.min() / 1e3:.1f} to {D.max() / 1e3:.1f}k px^2 over {MATCH_S[0]:.0f}-{MATCH_S[1]:.0f} px^2.',
         '  slope against the comparison range (fit error, sigma_x):']
    L += [f'    {m[0]:7.0f}-{m[1]:6.0f} px^2: {k:.3f} +- {u:.3f} rad/(rad/s), sigma_x = {x:.3f} rad'
          for m, k, u, x in t_scan]
    L += [f'AIA: {len(fr)} frames x {len(s)} bins, converged in {n_it} iterations; rms residual '
          f'{np.sqrt(np.mean(resid ** 2)):.2f} counts; frame phases span {np.ptp(phi_aia):.2f} rad '
          f'({np.degrees(np.ptp(phi_aia)):.0f} deg).',
          f'  bin phases psi(s): linear fit P = {P_aia / 1e3:.1f}k px^2, largest deviation {dev:.2f} rad '
          '(bins with > 20 % of the peak modulation).',
          '  slope against the largest s used, without and with a brightness factor per frame:']
    L += [f'    s < {m:7.0f} px^2, {"gain " if g else "plain"}: {k:.3f} +- {u:.3f} rad/(rad/s), '
          f'sigma_x = {x:.3f} rad' for m, g, k, u, x in a_scan]
    L += ['  AIA assumes that only the phase changes between frames. On synthetic profiles '
          '(tests/test_methods.py) it recovers a 3.3 rad sweep to 0.01 rad; a 5 % brightness drift '
          'leaves the slope right but puts single frames off by up to 0.1 rad, and a 1.7 % '
          'frame-to-frame change of the fringe period, as measured, gives errors of order 0.2 rad. '
          'Its spread over windows and variants here measures how far the photographs break that '
          'assumption.',
          'Per-frame phases in methods_frames.csv; statistics and the comparison with the other '
          'methods in sagnac.txt.']
    (HERE / 'methods.txt').write_text('\n'.join(L) + '\n', encoding='utf-8')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
