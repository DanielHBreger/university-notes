"""Every number quoted in RESULTS.md, SUMMARY.md and presenter_guide.md, computed with the presentation's own code
(slides.py, deck.py in this folder). Run from anywhere with Python 3.11+:  python collect_numbers.py
It writes numbers.json next to this file. Set SAGNAC_DATA to the 'Sagnac effect' folder if this folder is not in
'Sagnac effect/presentation_slides/source/'."""
import json, os, sys, tempfile
from pathlib import Path
import numpy as np
from scipy import stats
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
os.environ.setdefault('SAGNAC_OUT', tempfile.mkdtemp())        # slides.py would write images here, none are made
import slides as S
from deck import *

out = {}
# calibration
out['cal_keys'] = list(CAL.keys())
out['cal_fits_keys'] = {k: list(v.keys()) for k, v in CAL['fits_with_eps'].items()}
lin = CAL['linear']; out['lin'] = {k: lin[k] for k in lin if not isinstance(lin[k], (list, dict))}
for d in ('1', '2', '3'):
    f = CAL['fits_with_eps'][d]
    out[f'fit_eps_{d}'] = {k: f[k] for k in f if k in ('chi2', 'dof', 'p', 'coef')}
out['eps'] = EPS
tab = calibration_table()
out['cal_table'] = tab.round(5).to_dict('records')
out['stopwatch_rel_pct'] = [float((tab.u_sw / tab.omega).min() * 100), float((tab.u_sw / tab.omega).max() * 100)]
qc = CAL['fits_with_eps']['2']['coef']
out['omega90_lin'] = float(C0 + C1 * 90); out['omega90_quad'] = float(qc[0] + qc[1] * 90 + qc[2] * 8100)
Sall, laps = load_turns(); j = int(np.argmax(Sall)); lap = laps[j]
out['S90'] = dict(n=len(lap), T=float(lap.mean()), sd=float(lap.std(ddof=1)), omega=float(2 * np.pi / lap.mean()))
# per setting
sp = SPEEDS
first = FRAMES.groupby('S').time_s.min()
rows = []
for i, r in sp.iterrows():
    t = first[r.S]
    rows.append(dict(S=float(r.S), n=int(r.n), first=f'{int(t // 3600):02d}:{int(t % 3600 // 60):02d}:{int(t % 60):02d}',
                     omega=float(r.omega), u_omega=float(np.hypot(r.u_omega_curve, r.u_omega_repro)),
                     dphi=float(r.dphi), sem=float(r['sem']), P=float(r.P)))
out['settings'] = rows
out['spcols'] = list(sp.columns)
# theory
out['theory'] = [RES['theory'], RES['u_theory']]
# slide 8 regression
f = S.fit8('regression')
out['reg'] = {k: (float(v) if np.ndim(v) == 0 and v is not None else None) for k, v in f.items() if np.ndim(v) == 0}
out['reg_resid'] = [float(x) for x in f['resid']]
w = f['w']; y = sp.phi.to_numpy() - sp.phi.to_numpy()[0]
A = np.c_[w, np.ones(len(w))]; q, *_ = np.linalg.lstsq(A, y, rcond=None); r = y - A @ q
C = (r @ r / 8) * np.linalg.inv(A.T @ A)
out['intercept'] = [float(q[1]), float(np.sqrt(C[1, 1])), float(2 * stats.t.sf(abs(q[1] / np.sqrt(C[1, 1])), 8))]
a0 = (w @ y) / (w @ w); r0 = y - a0 * w
out['origin_fit'] = [float(a0), float(np.sqrt(r0 @ r0 / 9 / (w @ w)))]
fp = S.fit8('photo')
out['photo_weighted'] = dict(a=float(fp['a']), chi=float(fp['chi']), dof=fp['dof'], p=float(fp['p']))
# all 161 photos regression
fr = FRAMES
Af = np.c_[fr.omega, np.ones(len(fr))]; qf, *_ = np.linalg.lstsq(Af, fr.dphi, rcond=None); rf = fr.dphi - Af @ qf
Cf = (rf @ rf / (len(fr) - 2)) * np.linalg.inv(Af.T @ Af)
out['all_photos'] = [float(qf[0]), float(np.sqrt(Cf[0, 0]))]
# wander
wd = S.wander_sigma(); out['wander'] = [float(wd['sigma']), float(wd['u_sigma'])]
# final
fr_ = S.final_result()
out['final'] = {k: (float(v) if np.ndim(v) == 0 else [float(x) for x in v]) for k, v in fr_.items() if k not in ('fits', 'rel')}
out['final_rel'] = {k: float(v) for k, v in fr_['rel'].items()}
out['methods'] = [dict(a=float(m['a']), u_a=float(m['u_a']), s=float(m['s']), r2=float(m['r2']), z=float(m['z']), p_z=float(m['p_z'])) for m in fr_['fits']]
mchi = []
for col, sc in (('phi', 'sem'), ('phi_cos', 'sem_cos'), ('phi_template', 'sem_template')):   # each method's own SEM
    y_ = sp[col].to_numpy(); u_ = np.sqrt(sp[sc].to_numpy() ** 2 + (RES['slope'] * EPS * w) ** 2)
    A_ = np.c_[w, np.ones(len(w))]; W_ = 1 / u_ ** 2; C_ = np.linalg.inv(A_.T @ (A_ * W_[:, None])); q_ = C_ @ (A_.T @ (W_ * y_))
    ch = float(np.sum(((y_ - A_ @ q_) / u_) ** 2)); mchi.append([ch, 8, float(stats.chi2.sf(ch, 8))])
out['methods_photo_chi'] = mchi
th, uth = RES['theory'], RES['u_theory']
out['z_full'] = float((fr_['a'] - th) / np.hypot(fr_['u'], uth)); out['p_full'] = float(2 * stats.norm.sf(abs(out['z_full'])))
out['ratio_full'] = [float(fr_['a'] / th), float(fr_['a'] / th * np.hypot(fr_['u'] / fr_['a'], uth / th))]
out['c_over_c0'] = [float(fr_['c'] / C_DEF), float(fr_['u_c'] / C_DEF)]
rel_meas = np.sqrt(fr_['rel']['Phase-reading method'] ** 2 + fr_['rel']['Calibration model'] ** 2 + fr_['rel']['Regression and calibration line'] ** 2)
out['rel_meas'] = float(rel_meas)
out['N_implied'] = [float(N_TURNS * fr_['a'] / th), float(N_TURNS * fr_['a'] / th * rel_meas),
                    float(N_TURNS * fr_['a'] / th * np.hypot(rel_meas, 2 * .005 / R_LOOP))]
# largest shift
s_ = out['reg']['s']
out['max_shift'] = [float(sp.dphi.iloc[-1]), float(np.sqrt(2) * s_), float(sp.dphi.iloc[-1] / (2 * np.pi))]
out['dt_meas_fs'] = float(sp.dphi.iloc[-1] * LAM / (2 * np.pi * C_DEF) * 1e15)
out['dt_pred_fs'] = float(4 * N_TURNS * AREA * sp.omega.iloc[-1] / C_DEF ** 2 * 1e15)
out['period_fs'] = float(LAM / C_DEF * 1e15)
out['phi_pred_top'] = float(th * sp.omega.iloc[-1])
# scan
st = S.scan_stats()
def clean(d):
    return {k: (float(v) if np.ndim(v) == 0 else None) for k, v in d.items() if np.ndim(v) == 0}
out['scan'] = dict(orig=float(st['orig']), n_draws=st['n_draws'], fixed=clean(st['fixed']), frame=clean(st['frame']),
                   by_n={str(k): clean(v) for k, v in st['by_n'].items()})
x = st['fixed']['s']; x = x[x > 0]
out['scan_robust'] = float((np.percentile(x, 84) - np.percentile(x, 16)) / 2)
out['scan_pct'] = [float(np.percentile(x, 16)), float(np.percentile(x, 84))]
# drift
p, cov, pD = S.drift_fit()
out['drift'] = dict(a=float(p[0]), u_a=float(np.sqrt(cov[0, 0])), D=float(p[1] * 1e3), u_D=float(np.sqrt(cov[1, 1]) * 1e3),
                    rho=float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])), p_D=float(pD), r2=float(S.drift_fit.r2))
out['drift_within'] = dict(min=float(sp.drift.min() * 1e3), max=float(sp.drift.max() * 1e3), n_sig=int((sp.p_drift < .05).sum()),
                           mean=float(sp.drift.mean() * 1e3), u_mean=float(sp.drift.std(ddof=1) / np.sqrt(len(sp)) * 1e3))
# rings
d = FRAMES
u = d[['uR1', 'uR2', 'uR3']].values.ravel(); u = u[np.isfinite(u)]
out['rings'] = dict(rms_min=float(d.rms_px.min()), rms_max=float(d.rms_px.max()), rms_med=float(d.rms_px.median()),
                    uR_min=float(u.min()), uR_max=float(u.max()), uR_med=float(np.median(u)),
                    n_min=float(np.nanmin(d[['n1', 'n2', 'n3']].values)), n_med=float(np.nanmedian(d[['n1', 'n2', 'n3']].values)),
                    ucx_med=float(d.u_cx.median()), ucy_med=float(d.u_cy.median()),
                    P_med=float(d.P_px2.median()), P_min_set=float(sp.P.min()), P_max_set=float(sp.P.max()),
                    centre_sd=[float(d.cx.std()), float(d.cy.std())])
out['parabola'] = [float(f['t_q']), float(f['p_q'])]
json.dump(out, open(HERE / 'numbers.json', 'w'), indent=1, default=float)
print('wrote', HERE / 'numbers.json')
