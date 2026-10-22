"""
Uncertainty audit of the Sagnac analysis (measurements of 11 October 2026).

Recomputes every uncertainty component of the final result from the inputs and the existing
outputs, states the model behind each one, and checks the components that rest on an assumption.
It changes nothing in the pipeline: the final numbers stay those of RESULTS.md (computed by
presentation_slides/source/slides.py), and every value here that RESULTS.md also quotes is checked
against presentation_slides/source/numbers.json at the end (section 12).

Reads:  Rotation calibration.csv, calibration.json, fringes_frames.csv,
        sagnac_update/sagnac_speeds.csv, uncertainty/instrument_inputs.json,
        presentation_slides/source/numbers.json, and the final-analysis functions of
        presentation_slides/source/slides.py (fit8, final_result, cal_model_slopes, wander_sigma)
Writes: uncertainty/budget.json     every number below, unrounded
        uncertainty/calculate.txt   the same as text, with report-form values

Sections (the same numbering in calculate.txt and uncertainty/REPORT.md):
     1. stopwatch: u(T) from lap splits, the timing error of one press, the display resolution
     2. speed reproducibility eps: re-derived, and why it is taken from the quadratic
     3. calibration line: coefficients, their correlation, omega at the photo settings
     4. phase of one photo: precision from the ring radii against the photo-to-photo spread
     5. setting means: SEM, the shared rest reference, wander between the halves of a setting
     6. regression: slope, its parts, comparison with the prediction
     7. final slope: the three parts and other ways to combine them
     8. prediction a_th: sensitivity to N, R, lambda; N as a standard uncertainty or as a bound
     9. speed of light: first-order budget and a Monte Carlo check
    10. largest shift and the travel-time difference
    11. inputs that are missing (instrument_inputs.json nulls)
    12. agreement with numbers.json

Usage (a few seconds):   python uncertainty/calculate.py
"""
import json
import math
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

UNC = Path(__file__).resolve().parent
HERE = UNC.parent                                      # the 'Sagnac effect' folder
SRC = HERE / 'presentation_slides' / 'source'
os.environ.setdefault('SAGNAC_DATA', str(HERE))
os.environ.setdefault('SAGNAC_OUT', tempfile.mkdtemp())   # the slide code is only read from
sys.path[:0] = [str(HERE), str(SRC)]

import slides as S  # noqa: E402  (imports deck.py: the data of the final analysis)
import deck as D  # noqa: E402

INP = json.loads((UNC / 'instrument_inputs.json').read_text())
NUM = json.loads((SRC / 'numbers.json').read_text())
P_GUESS = 107.5e3          # px^2, as fringes.P_GUESS: numbers the rings across gaps
N_MC = 400_000             # Monte Carlo draws for section 9
RNG = np.random.default_rng(20261011)

out = {}
txt = []


def line(s=''):
    txt.append(s)


def pm(v, u):
    """Report form: the uncertainty to one significant figure (half up), the value to the same digit."""
    if not (u > 0 and np.isfinite(u)):
        return f'{v:g}'
    e = math.floor(math.log10(u))
    d = round(u / 10 ** e + 1e-12)            # first significant digit, half up
    if d == 10:
        d, e = 1, e + 1
    dec = max(0, -e)
    q = 10 ** e
    vv = math.floor(abs(v) / q + 0.5 + 1e-12) * q * (1 if v >= 0 else -1)
    return f'{vv:.{dec}f} ± {d * q:.{dec}f}'


def p_two_sided(z):
    return float(2 * stats.norm.sf(abs(z)))


# ---- 1. stopwatch -----------------------------------------------------------------------------
S_cal, laps = D.load_turns()
rows = []
for s, lap in zip(S_cal, laps):
    n = len(lap)
    T = lap.mean()
    sd = lap.std(ddof=1)
    u_T = sd / math.sqrt(n * (n + 1))                  # lap splits share their presses
    s_press = sd * math.sqrt(n / (2 * (n + 1)))        # timing error of one press implied by the lap spread
    u_res = INP['stopwatch_display_resolution_s'] / math.sqrt(12) / math.sqrt(n)   # rounding of n displayed laps
    w = 2 * math.pi / T
    rows.append(dict(S=float(s), n=n, T_s=T, SD_s=sd, u_T_s=u_T, sigma_press_s=s_press, u_resolution_s=u_res,
                     omega=w, u_omega_stopwatch=w * u_T / T, rel_pct=100 * u_T / T))
sw = pd.DataFrame(rows)
lag = D.CAL['lag1']
out['1_stopwatch'] = dict(per_setting=sw.to_dict('records'), lag1=lag,
                          rel_pct_range=[float(sw.rel_pct.min()), float(sw.rel_pct.max())])
line('1. STOPWATCH (lap splits)')
line('   T = (t_n - t_0)/n. Each press has an independent error sigma; the laps share presses, so')
line('   E[SD^2] = 2 sigma^2 (n+1)/n and Var(T) = 2 sigma^2/n^2  ->  u(T) = SD/sqrt(n(n+1)).')
line('       S   n    <T> (s)   SD (s)  u(T) (s)  sigma_press (s)  u_display (s)  u(omega)/omega')
for r in rows:
    line(f'   {r["S"]:6.2f} {r["n"]:2d}  {r["T_s"]:8.3f}  {r["SD_s"]:7.3f}  {r["u_T_s"]:7.4f}  {r["sigma_press_s"]:12.3f}'
         f'  {r["u_resolution_s"]:12.4f}  {r["rel_pct"]:9.2f} %')
line(f'   Lag-1 correlation of neighboring laps {lag["observed"]:.2f}; splits predict {lag["split_mean"]:.2f} '
     f'+- {lag["split_sd"]:.2f} (simulated); independent laps give a value this negative with p = '
     f'{lag["p_independent"]:.2f}.')
line(f'   The 0.01 s display resolution: rounding of the n displayed laps adds {sw.u_resolution_s.min():.4f}-'
     f'{sw.u_resolution_s.max():.4f} s to <T>, below every u(T),')
line('   and it is already part of the lap spread that u(T) is computed from (statistical).')
line()

# ---- 2. speed reproducibility eps -------------------------------------------------------------
tab = D.calibration_table()                       # S, n, T, omega, u_sw, u (= sqrt(u_sw^2 + (eps omega)^2))


def chi2_poly(deg, eps):
    """chi2 of the weighted polynomial fit omega(S) with errors sqrt(u_sw^2 + (eps omega)^2)."""
    Sx, w, usw = tab.S.to_numpy(float), tab.omega.to_numpy(float), tab.u_sw.to_numpy(float)
    u = np.hypot(usw, eps * w)
    A = np.vander(Sx, deg + 1, increasing=True)
    W = 1 / u ** 2
    C = np.linalg.inv(A.T @ (A * W[:, None]))
    q = C @ (A.T @ (W * w))
    return float(np.sum(((w - A @ q) / u) ** 2)), q, C


dof = {1: len(tab) - 2, 2: len(tab) - 3, 3: len(tab) - 4}
eps_quad = optimize.brentq(lambda e: chi2_poly(2, e)[0] - dof[2], 1e-6, 0.05)
eps_lin = optimize.brentq(lambda e: chi2_poly(1, e)[0] - dof[1], 1e-6, 0.2)
chi_tab = {deg: dict(stopwatch_only=chi2_poly(deg, 0)[0], with_eps=chi2_poly(deg, eps_quad)[0], dof=dof[deg])
           for deg in (1, 2, 3)}
for deg in chi_tab:
    for k in ('stopwatch_only', 'with_eps'):
        chi_tab[deg]['p_' + k] = float(stats.chi2.sf(chi_tab[deg][k], dof[deg]))
out['2_eps'] = dict(eps_from_quadratic=eps_quad, eps_calibration_json=D.EPS, eps_if_set_by_the_line=eps_lin,
                    fits=chi_tab)
line('2. SPEED REPRODUCIBILITY eps (the table speed at a fixed setting)')
line(f'   eps solves chi2/dof = 1 for the quadratic: {100 * eps_quad:.3f} % of omega (calibration.json: '
     f'{100 * D.EPS:.3f} %). Statistical, independent between settings, added in quadrature to u_stopwatch.')
for deg, name in ((1, 'linear   '), (2, 'quadratic'), (3, 'cubic    ')):
    c = chi_tab[deg]
    line(f'   {name}  stopwatch only chi2/dof = {c["stopwatch_only"]:7.1f}/{c["dof"]}  p = {c["p_stopwatch_only"]:.1e};'
         f'  with eps chi2/dof = {c["with_eps"]:5.1f}/{c["dof"]}  p = {c["p_with_eps"]:.2g}')
line(f'   Set by the line instead, eps would be {100 * eps_lin:.2f} % of omega: the line\'s systematic misfit at high S')
line('   would then be counted as random speed noise. Taking eps from the quadratic keeps the two apart; the')
line('   misfit of the line enters as the calibration-model term (section 7).')
line()

# ---- 3. calibration line ------------------------------------------------------------------------
lin = D.CAL['fits_with_eps']['1']
c0, c1 = lin['coef']
Cl = np.array(lin['cov'])
rho = Cl[0, 1] / math.sqrt(Cl[0, 0] * Cl[1, 1])
sp = D.SPEEDS
rows = []
for s in sp.S[1:]:
    w = c0 + c1 * s
    J = np.array([1.0, s])
    u_curve = float(math.sqrt(J @ Cl @ J))
    rows.append(dict(S=float(s), omega=w, u_curve=u_curve, u_repro=D.EPS * w, u=math.hypot(u_curve, D.EPS * w)))
om = pd.DataFrame(rows)
chk = np.allclose(om.u_curve, sp.u_omega_curve[1:], rtol=1e-3) and np.allclose(om.u_repro, sp.u_omega_repro[1:], rtol=1e-3)
out['3_calibration_line'] = dict(k=c1, u_k=math.sqrt(Cl[1, 1]), omega0=c0, u_omega0=math.sqrt(Cl[0, 0]),
                                 correlation=rho, photo_settings=om.to_dict('records'),
                                 matches_sagnac_speeds_csv=bool(chk))
line('3. CALIBRATION LINE omega = k S + omega0 (used for every photo)')
line(f'   k = {c1:.6f} +- {math.sqrt(Cl[1, 1]):.6f} rad/s per unit S, omega0 = {c0:.5f} +- {math.sqrt(Cl[0, 0]):.5f} rad/s,'
     f' correlation {rho:.2f}')
line('   u_curve = sqrt((1, S) Cov (1, S)^T): calibration, the same line for every setting, so correlated between')
line('   settings; it enters the slope through the redrawn lines (section 6). u_repro = eps omega: independent.')
line('       S   omega (rad/s)  u_curve   u_repro   u (combined)')
for r in rows:
    line(f'   {r["S"]:6.2f}   {r["omega"]:.4f}      {r["u_curve"]:.5f}   {r["u_repro"]:.5f}   {r["u"]:.5f}  '
         f'({100 * r["u"] / r["omega"]:.2f} %)')
line(f'   Same values as sagnac_update/sagnac_speeds.csv: {chk}.')
line()

# ---- 4. phase of one photo ---------------------------------------------------------------------
FR = D.FRAMES


def phi_c(R):
    """Phase at the ring center from the dark-ring radii, exactly as fringes.ring_phase (unwrapped)."""
    s = np.sort(np.asarray(R, float)) ** 2
    m = np.r_[0, np.cumsum(np.rint(np.diff(s) / P_GUESS))]
    P, s0 = np.polyfit(m, s, 1)
    return math.pi - 2 * math.pi * s0 / P


u_photo = []
for _, r in FR.iterrows():
    R = np.array([r.R1, r.R2, r.R3])
    u = np.array([r.uR1, r.uR2, r.uR3])
    ok = np.isfinite(R)
    R, u = R[ok], u[ok]
    o = np.argsort(R)
    R, u = R[o], u[o]
    base = phi_c(R)
    J = np.array([(phi_c(R + 1e-3 * np.eye(len(R))[k]) - base) / 1e-3 for k in range(len(R))])
    u_photo.append(float(math.sqrt(np.sum((J * u) ** 2))))
FR = FR.assign(u_phi_rings=u_photo)
per = []
for s, g in FR.groupby('S'):
    row = sp.set_index('S').loc[s]
    per.append(dict(S=float(s), photos=len(g), median_u_from_rings=float(g.u_phi_rings.median()),
                    SD_of_photos=float(g.phi_rad.std(ddof=1)), drift_mrad_per_s=1e3 * float(row.drift),
                    p_drift=float(row.p_drift)))
per = pd.DataFrame(per)
per['ratio'] = per.SD_of_photos / per.median_u_from_rings
out['4_photo_phase'] = dict(per_setting=per.to_dict('records'),
                            u_from_rings_range=[float(FR.u_phi_rings.min()), float(FR.u_phi_rings.max())],
                            u_from_rings_median=float(FR.u_phi_rings.median()))
line('4. PHASE OF ONE PHOTO: precision from the ring radii against the spread of the photos')
line('   u(phi_c) of a photo: the circle-fit uncertainties of its ring radii (fringes_frames.csv uR1..uR3)')
line('   propagated through R_m^2 = P (m + eta), phi_c = pi - 2 pi eta, as the radii were independent.')
line(f'   Over all 161 photos: {FR.u_phi_rings.min():.3f} to {FR.u_phi_rings.max():.3f} rad, median '
     f'{FR.u_phi_rings.median():.3f} rad.')
line('       S   photos  median u(phi_c)  SD of the photos  SD/u   drift within setting (mrad/s, p)')
for _, r in per.iterrows():
    line(f'   {r.S:6.2f}  {int(r.photos):3d}      {r.median_u_from_rings:.3f}           {r.SD_of_photos:.3f}'
         f'        {r.ratio:4.1f}   {r.drift_mrad_per_s:+6.1f}  p = {r.p_drift:.3f}')
hi = per[per.S >= 55.83]
lo = per[(per.S > 0) & (per.S < 55.83)]
line(f'   From S = 55.83 up the photos of a setting differ by about what the ring fit allows (SD/u '
     f'{hi.ratio.min():.1f}-{hi.ratio.max():.1f});')
line(f'   from 12.59 to 44.07 they differ {lo.ratio.min():.1f}-{lo.ratio.max():.1f} times more, and there the phase '
     f'drifts within the setting (p = {lo.p_drift.min():.0e} to {lo.p_drift.max():.3f}).')
out['4_photo_phase'].update(ratio_high_settings=[float(hi.ratio.min()), float(hi.ratio.max())],
                            ratio_low_settings=[float(lo.ratio.min()), float(lo.ratio.max())])
line('   Caveat: neighboring rays share speckle, so the ray minima are not fully independent and the radius')
line('   uncertainties (and u(phi_c)) are probably somewhat too small. Not used in the result: the setting')
line('   means take their uncertainty from the photos themselves (SEM), and the line from the regression.')
line()

# ---- 5. setting means -------------------------------------------------------------------------
wd = S.wander_sigma()
out['5_setting_means'] = dict(SEM=dict(zip(map(float, sp.S), map(float, sp['sem']))),
                              rest_SEM=float(sp['sem'].iloc[0]),
                              wander_sigma=float(wd['sigma']), u_wander_sigma=float(wd['u_sigma']),
                              halves_difference=dict(zip(map(float, wd['S']), map(float, wd['d']))),
                              halves_u=dict(zip(map(float, wd['S']), map(float, wd['ud']))))
zh = np.abs(np.asarray(wd['d'])) / np.asarray(wd['ud'])
nsig = int(np.sum(zh > 2))
line('5. SETTING MEANS')
line(f'   Mean of the photos, u = SEM = SD/sqrt(n): {sp["sem"].min():.3f} to {sp["sem"].max():.3f} rad (statistical).')
line(f'   The rest mean (SEM {sp["sem"].iloc[0]:.3f} rad) is subtracted from every setting: a common offset, which a line')
line('   with a free intercept absorbs. It enters the slope only as the rest point\'s own deviation.')
line(f'   Second half minus first half of each setting: beyond 2 of its photo-noise errors in {nsig} of 10 settings '
     f'({zh[zh > 2].min():.1f}-{zh.max():.1f} errors).')
line(f'   Implied wander of a setting mean: {wd["sigma"]:.3f} +- {wd["u_sigma"]:.3f} rad (bootstrap over settings),')
line('   so the SEM, which assumes independent photos, understates how well a setting mean is known.')
line()

# ---- 6. regression ----------------------------------------------------------------------------
f8 = S.fit8('regression')
TH, UTH = D.RES['theory'], D.RES['u_theory']
z_reg = (f8['a'] - TH) / math.hypot(f8['u_a'], UTH)
fp = S.fit8('photo')
out['6_regression'] = dict(a=f8['a'], u_stat=f8['u_stat'], u_cal=f8['u_cal'], u_a=f8['u_a'], s=f8['s'],
                           r2=f8['r2'], b=f8['b'] - float(sp.phi.iloc[0]), z_vs_theory=z_reg,
                           p_vs_theory=p_two_sided(z_reg), photo_weighted=dict(a=fp['a'], chi2=fp['chi'],
                                                                               dof=fp['dof'], p=fp['p']))
line('6. REGRESSION Delta phi = a omega + b (ten setting means, equal weights; Schwartz Eq. 7)')
line(f'   s = sqrt(sum(Delta phi - line)^2 / 8) = {f8["s"]:.4f} rad: the uncertainty of every point, taken from the points.')
line(f'   u_stat = s / sqrt(sum (omega - mean)^2) = {f8["u_stat"]:.4f}; calibration line (4000 redrawn lines) '
     f'{f8["u_cal"]:.4f}; together {f8["u_a"]:.4f}.')
line(f'   a = {f8["a"]:.4f} -> {pm(f8["a"], f8["u_a"])} rad/(rad/s); R^2 = {f8["r2"]:.4f}.')
line(f'   With the photos\' SEM and eps as weights instead: chi2/dof = {fp["chi"]:.1f}/{fp["dof"]}, p = {fp["p"]:.0e}:')
line('   those error bars are too small for the differences between settings (section 5), hence the regression.')
line(f'   Against a_th = {TH:.4f} +- {UTH:.4f}: z = {z_reg:.2f}, p = {p_two_sided(z_reg):.2f} (regression error only).')
# What moves the setting means (quoted in RESULTS.md, not in numbers.json): the points' distances from the line
# against each setting's fringe period P. P_mean is the plain mean of the photos' P; sagnac_speeds.csv column P
# (Table 2) is the median.
gP = FR.groupby('S')
P_mean, P_med, Pcos_mean = (gP[c].agg(f).to_numpy() for c, f in (('P_px2', 'mean'), ('P_px2', 'median'),
                                                                  ('P_cos_px2', 'mean')))
assert np.allclose(np.array(sorted(gP.groups)), sp.S.to_numpy())      # same settings, same order as the residuals
res = np.asarray(f8['resid'])
r_mean, p_mean = stats.pearsonr(res, P_mean)
r_med, p_med = stats.pearsonr(res, P_med)
r_cos, p_cos = stats.pearsonr(P_mean, Pcos_mean)
out['6_regression']['distance_vs_P'] = dict(r_mean_P=r_mean, p_mean_P=p_mean, r_median_P=r_med, p_median_P=p_med,
                                            r_ring_vs_cosine_period=r_cos, p_ring_vs_cosine_period=p_cos,
                                            P_mean_range=[float(P_mean.min()), float(P_mean.max())],
                                            P_median_range=[float(P_med.min()), float(P_med.max())])
line(f'   Distance of each point from the line against the setting\'s mean fringe period P: r = {r_mean:+.2f}, '
     f'p = {p_mean:.3f}')
line(f'   (with the median P of Table 2: r = {r_med:+.2f}, p = {p_med:.3f}); ring-fit and cosine-fit mean periods '
     f'agree, r = {r_cos:.2f}.')
line(f'   P per setting, 10^5 px^2: medians {P_med.min() / 1e5:.3f}-{P_med.max() / 1e5:.3f}, '
     f'means {P_mean.min() / 1e5:.3f}-{P_mean.max() / 1e5:.3f}.')
line()

# ---- 7. final slope -----------------------------------------------------------------------------
fin = S.final_result()
a, u_reg, u_meth, u_model, u = fin['a'], fin['u_reg'], fin['u_meth'], fin['u_model'], fin['u']
slopes = np.array(fin['slopes'])
# The line's coefficient covariance (calibration.py) is not scaled by its chi2/dof (47.2/7 = 6.7). Scaled, the
# calibration-line part of the slope grows by sqrt(chi2/dof): the redrawn-line spread is proportional to the
# coefficients' standard deviations to first order (k is known to 0.2 %, 0.5 % scaled). The final analysis covers the
# line's misfit by the calibration-model term instead; using both would count it twice. Shown only as an alternative.
scale_cal = math.sqrt(lin['chi2'] / lin['dof'])
u_reg_scaled = math.hypot(f8['u_stat'], f8['u_cal'] * scale_cal)
alts = {'quadrature (used)': u,
        'straight sum (worst case)': u_reg + u_meth + u_model,
        'method as range/sqrt(12)': math.sqrt(u_reg ** 2 + (np.ptp(slopes) / math.sqrt(12)) ** 2 + u_model ** 2),
        'method as the full range': math.sqrt(u_reg ** 2 + np.ptp(slopes) ** 2 + u_model ** 2),
        'line covariance x chi2/dof': math.sqrt(u_reg_scaled ** 2 + u_meth ** 2 + u_model ** 2)}
alt_z = {k: (a - TH) / math.hypot(v, UTH) for k, v in alts.items()}
out['7_final_slope'] = dict(a=a, parts=dict(regression_and_line=u_reg, method=u_meth, calibration_model=u_model),
                            u=u, rel_pct=100 * u / a, alternatives={k: dict(u=v, z=alt_z[k], p=p_two_sided(alt_z[k]))
                                                                     for k, v in alts.items()})
line('7. FINAL SLOPE')
line(f'   regression and calibration line {u_reg:.4f} ({100 * u_reg / a:.1f} %), statistical')
line(f'   phase-reading method {u_meth:.4f} ({100 * u_meth / a:.1f} %): half the range of {", ".join(f"{x:.4f}" for x in slopes)}')
line(f'   calibration model {u_model:.4f} ({100 * u_model / a:.1f} %): |a(quadratic) - a(linear)|')
line(f'   a = {a:.4f} +- {u:.4f} ({100 * u / a:.1f} %) -> {pm(a, u)} rad/(rad/s)')
for k, v in alts.items():
    line(f'   {k:28s} u = {v:.3f}   z = {alt_z[k]:.2f}   p = {p_two_sided(alt_z[k]):.2f}')
line(f'   The line\'s covariance is not scaled by its chi2/dof ({lin["chi2"] / lin["dof"]:.2f}); scaled, the calibration-line')
line(f'   part would be {f8["u_cal"] * scale_cal:.4f} instead of {f8["u_cal"]:.4f}. The line\'s misfit is already in the')
line('   calibration-model term, so the result keeps the unscaled value; using both would count it twice.')
out['7_final_slope']['u_cal_if_line_covariance_scaled'] = f8['u_cal'] * scale_cal
line()

# ---- 8. prediction a_th -----------------------------------------------------------------------
N, uN = INP['fiber_turns_N'], INP['fiber_turns_u']
R, uR = INP['loop_radius_m'], INP['loop_radius_u_m']
lam = INP['wavelength_m']
c0_def = INP['speed_of_light_defined_m_per_s']
a_th = 8 * math.pi ** 2 * N * R ** 2 / (lam * c0_def)
rel_N, rel_R = uN / N, 2 * uR / R
u_th = a_th * math.hypot(rel_N, rel_R)
per_nm = 1e-9 / lam                                    # relative change of a_th (and c) per nm of lambda
rel_N_bound = uN / math.sqrt(3) / N                    # +-1 turn read as a bound: flat distribution
u_th_bound = a_th * math.hypot(rel_N_bound, rel_R)
z_bound = (a - a_th) / math.hypot(u, u_th_bound)
own = math.sqrt(u_reg ** 2 + u_meth ** 2 + u_model ** 2) / a
out['8_prediction'] = dict(a_th=a_th, u_a_th=u_th, rel_N=rel_N, rel_R=rel_R, rel_per_nm_lambda=per_nm,
                           lambda_nm_matching_own_uncertainty=own / per_nm,
                           N_as_bound=dict(u_a_th=u_th_bound, z=z_bound, p=p_two_sided(z_bound)))
line('8. PREDICTION a_th = 8 pi^2 N R^2 / (lambda c)')
line(f'   a_th = {a_th:.4f} s; N: {100 * rel_N:.2f} %, R (enters squared): {100 * rel_R:.2f} %, together '
     f'{100 * math.hypot(rel_N, rel_R):.2f} % -> {pm(a_th, u_th)} rad/(rad/s). c is exact (defined).')
line(f'   lambda: no uncertainty given (missing, not zero). a_th and c change by {100 * per_nm:.3f} % per nm; lambda would')
line(f'   have to be off by {own / per_nm:.0f} nm to matter as much as all our own measurement uncertainties together.')
line(f'   If 13 +- 1 is a bound (flat distribution, u = 0.58 turn): u(a_th) = {u_th_bound:.3f}, the comparison '
     f'z = {z_bound:.2f}, p = {p_two_sided(z_bound):.2f}.')
line()

# ---- 9. speed of light ------------------------------------------------------------------------
c_val = fin['c']
rel = fin['rel']
rel_tot = fin['rel_tot']
Ns = RNG.normal(N, uN, N_MC)
Rs = RNG.normal(R, uR, N_MC)
As = RNG.normal(a, u, N_MC)
cs = 8 * math.pi ** 2 * Ns * Rs ** 2 / (lam * As)
q16, q50, q84 = np.percentile(cs, [16, 50, 84])
p_mc = 2 * min(np.mean(cs < c0_def), np.mean(cs > c0_def))
rel_bound = math.sqrt(rel_N_bound ** 2 + rel_R ** 2 + own ** 2)
out['9_speed_of_light'] = dict(c=c_val, u_c=fin['u_c'], rel=rel, rel_total=rel_tot, z=fin['z'], p=fin['p_z'],
                               monte_carlo=dict(draws=N_MC, mean=float(cs.mean()), median=float(q50),
                                                sd=float(cs.std()), q16=float(q16), q84=float(q84), p_two_sided=float(p_mc)),
                               rel_if_N_is_a_bound=rel_bound, N_implied=NUM['N_implied'])
line('9. SPEED OF LIGHT c = 8 pi^2 N R^2 / (lambda a)')
line(f'   c = {c_val / 1e8:.4f} x 10^8 m/s; relative parts (first order, c ~ N R^2 / a):')
for k, v in sorted(rel.items(), key=lambda kv: -kv[1]):
    line(f'      {k.replace("$", ""):38s} {100 * v:5.2f} %')
line(f'      total {100 * rel_tot:.2f} % -> ({pm(c_val / 1e8, fin["u_c"] / 1e8)}) x 10^8 m/s, z = {fin["z"]:.2f}, p = {fin["p_z"]:.2f}')
line(f'   Monte Carlo ({N_MC} draws of N, R, a, normal): median {q50 / 1e8:.3f}, 16-84 % range {q16 / 1e8:.3f}-{q84 / 1e8:.3f},')
line(f'   SD {cs.std() / 1e8:.3f} x 10^8 m/s; two-sided p of c0 = {p_mc:.2f}. Below the median {(q50 - q16) / 1e8:.3f}, above '
     f'it {(q84 - q50) / 1e8:.3f} x 10^8 m/s:')
line('   c ~ R^2/a makes the distribution slightly lopsided, too little to change the first-order result.')
line(f'   With N as a bound: {100 * rel_bound:.1f} %. Read the other way, N = {NUM["N_implied"][0]:.1f} +- '
     f'{NUM["N_implied"][1]:.1f} turns (+- {NUM["N_implied"][2]:.1f} with R).')
line()

# ---- 10. largest shift ----------------------------------------------------------------------
top = float(sp.dphi.iloc[-1])
u_top = math.sqrt(2) * f8['s']
dt = top * lam / (2 * math.pi * c0_def)
out['10_largest_shift'] = dict(dphi=top, u=u_top, fringes=top / (2 * math.pi), u_fringes=u_top / (2 * math.pi),
                               dt_fs=dt * 1e15, u_dt_fs=dt * 1e15 * u_top / top,
                               predicted_dphi=TH * float(sp.omega.iloc[-1]))
line('10. LARGEST SHIFT (top setting minus rest)')
line(f'   {top:.3f} +- {u_top:.3f} rad (u = sqrt(2) s: two setting means) = {top / (2 * math.pi):.3f} +- '
     f'{u_top / (2 * math.pi):.3f} fringe; Delta t = Delta phi lambda/(2 pi c) = {dt * 1e15:.3f} +- '
     f'{dt * 1e15 * u_top / top:.3f} fs (lambda taken as exact). Predicted {TH * sp.omega.iloc[-1]:.2f} rad.')
line()

# ---- 11. missing inputs -----------------------------------------------------------------------
missing = [k for k, v in INP.items() if v is None]
out['11_missing'] = missing
line('11. INPUTS NOT KNOWN (instrument_inputs.json: null; never treated as zero)')
for k in missing:
    line(f'   {k}')
line('   Not tested by the data: the sign of the shift (rotation direction), a steady drift in time')
line('   (settings in increasing order), omega during the photographs (only the calibration at other times).')
line()

# ---- 12. agreement with numbers.json ----------------------------------------------------------
checks = {'regression slope': (f8['a'], NUM['reg']['a']), 'regression u': (f8['u_a'], NUM['reg']['u_a']),
          's': (f8['s'], NUM['reg']['s']), 'final u': (u, NUM['final']['u']),
          'c': (c_val, NUM['final']['c']), 'u_c': (fin['u_c'], NUM['final']['u_c']),
          'a_th': (a_th, NUM['theory'][0]), 'u_a_th': (u_th, NUM['theory'][1]),
          'eps': (D.EPS, NUM['eps']), 'wander': (wd['sigma'], NUM['wander'][0])}
bad = {k: v for k, v in checks.items() if not math.isclose(v[0], v[1], rel_tol=1e-6)}
out['12_numbers_json_agreement'] = dict(checked=list(checks), differing=bad)
line('12. AGREEMENT WITH presentation_slides/source/numbers.json')
line(f'   {len(checks)} values checked, {len(bad)} differ' + (f': {bad}' if bad else '.'))

(UNC / 'budget.json').write_text(json.dumps(out, indent=1, default=float))
(UNC / 'calculate.txt').write_text('\n'.join(txt) + '\n', encoding='utf-8')
print('\n'.join(txt))
print('wrote', UNC / 'budget.json', 'and', UNC / 'calculate.txt')
