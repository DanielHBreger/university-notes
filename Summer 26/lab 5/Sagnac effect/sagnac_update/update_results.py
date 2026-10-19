"""Refresh the readable results from the selected calibration and saved fit outputs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
from sagnac import round1

HERE = Path(__file__).resolve().parent


def main():
    r = json.loads((HERE/'sagnac.json').read_text())
    c = json.loads((HERE/'calibration.json').read_text())
    optional = json.loads((HERE/'optional_quadratic/sagnac.json').read_text())
    scans = json.loads((HERE/'methods.json').read_text())
    draws = [d['slope'] for d in scans['points_scan']]
    speeds = pd.read_csv(HERE/'sagnac_speeds.csv')
    assert r['calibration_degree'] == c['conversion']['degree'] == scans['calibration_degree'] == 1
    value = lambda x, u: round1(x, u).replace('+-', '±')
    total_a = np.hypot(r['u_slope'], r['u_method'])
    z = (r['slope']-r['theory'])/np.hypot(total_a, r['u_theory'])
    maximum = speeds.dphi.iloc[-1]
    u_max = np.hypot(speeds.u_fit.iloc[-1], speeds.u_fit.iloc[0])
    sc = r['speed_of_light']
    delta = 100*(optional['slope']/r['slope']-1)
    common = [
        '# Sagnac effect: current linear-calibration results',
        '', 'Measurements: 11 October 2026. Standard uncertainties, rounded to one significant figure in report values. '
        'Comparisons use unrounded values. Full numerical output: `calibration.txt`, `sagnac.txt`, `sagnac_speeds.csv`.',
        '', '## Setup and analysis', '',
        'N = 13 ± 1 turns, R = 57.0 ± 0.5 cm, λ = 640 nm. The large loop is centred on the rotation axis. '
        'The data contain 161 photographs at rest and nine moving settings, plus rotation timing at nine settings. '
        'The off-centre-loop task was not performed.', '',
        'Dark minima along rays over −75° to 180° define concentric circles. A fit of squared ring radius '
        'against ring order gives the fringe period P and centre phase φc. Differences from rest give Δφ. '
        'The rings contract with increasing rotation and are equally spaced in r², with P ≈ 1.074 × 10⁵ px².',
        '', '## Default rotation calibration', '',
        f'ω(S) = ω₀ + kS, with k = {value(c["linear"]["k"], c["linear"]["u_k"])} rad/s per unit S '
        f'and ω₀ = {value(c["linear"]["omega0"], c["linear"]["u_omega0"])} rad/s. '
        'This line supplies every moving-setting angular velocity below. The rest setting has ω = 0 exactly.', '',
        'Stopwatch splits give u(T) = SD/√[n(n+1)]. The original speed reproducibility '
        f'ε = {100*c["conversion"]["eps"]:.2f}% is retained, estimated from the quadratic residuals. '
        'The selected linear coefficients and their covariance are unchanged from the saved straight-line fit.', '',
        '## Phase shift and rotation', '',
        '| S | Frames | ω (rad/s) | Δφ (rad) |', '|---|---:|---:|---:|',
    ]
    for row in speeds.itertuples():
        uw = np.hypot(row.u_omega_curve, row.u_omega_repro)
        omega = '0' if row.S == 0 else value(row.omega, uw)
        common.append(f'| {row.S:g} | {row.n} | {omega} | {value(row.dphi, row.u_fit)} |')
    common += [
        '', 'Phase entries carry the per-setting fit weight. They are not independent errors on differences '
        'from the shared rest reference. The fit uses absolute centre phases with a free intercept.', '',
        '| Quantity | Result |', '|---|---|',
        f'| Slope, fit + calibration | {value(r["slope"], r["u_slope"])} rad/(rad/s) |',
        f'| Slope, including method contribution | {value(r["slope"], total_a)} rad/(rad/s) |',
        f'| Theory 8π²NR²/(λc) | {value(r["theory"], r["u_theory"])} rad/(rad/s) |',
        f'| Difference from theory, including method | {z:.2f} combined σ, p = {2*stats.norm.sf(abs(z)):.2f} |',
        f'| Maximum measured shift | {value(maximum, u_max)} rad = {value(maximum/(2*np.pi), u_max/(2*np.pi))} fringe |',
        f'| Quadratic phase term | {value(r["curvature_q"], r["u_q"])} rad/(rad/s)² |',
        '', '## Phase-extraction methods', '',
        '| Method | Slope a, rad/(rad/s) | Phase observable |', '|---|---:|---|',
        f'| Minima (main) | {value(r["slope"], r["u_slope"])} | Ring centre |',
    ]
    for name, observable in [('cosine fit', 'Ring centre'), ('template matching', 'Middle of comparison window'),
                              ('random points', 'Local phase at ten fixed pixels; excluded from method spread')]:
        m = r['methods'][name]
        common.append(f'| {name} | {value(m["slope"], m["u_slope"])} | {observable} |')
    common += [
        '', f'Errors in this table combine phase-fit and calibration covariance. Half the range of the three retained '
        f'slopes gives a method contribution of {r["u_method"]:.5f} rad/(rad/s), or {100*r["u_method"]/r["slope"]:.1f}%. '
        'The procedures use the same photographs and are not independent measurements. The minima method '
        'does not assume an intensity envelope. The random-points method reads ten pixels at fixed image positions; '
        f'it is excluded because its slope depends on which pixels are drawn (central 68 % of {len(draws)} draws: '
        f'{np.percentile(draws, 16):.2f}–{np.percentile(draws, 84):.2f} rad/(rad/s)). '
        'The refreshed robustness scans are in `methods.txt` and `methods.json`.',
        '', '## Inferred speed of light and uncertainty', '',
        'c = 8π²NR²/(λa). This is the same theory comparison expressed as an inferred c, not an independent test.', '',
        f'**c = ({value(sc["c"]/1e8, sc["u_c"]/1e8)}) × 10⁸ m/s**, '
        f'c/c₀ = {value(sc["ratio_to_defined"], sc["u_ratio_to_defined"])}. '
        f'Difference from the defined value: {sc["z"]:.2f} σ, p = {sc["p_z"]:.2f}.', '',
        '| Contribution | Relative standard uncertainty in c |', '|---|---:|',
    ]
    common += [f'| {key} | {100*amount/sc["c"]:.2f}% |' for key, amount in sc['components'].items()]
    common += [f'| Combined | {100*sc["u_c"]/sc["c"]:.2f}% |', '',
        f'The measurement alone contributes {100*sc["u_measurement"]/sc["c"]:.2f}%. '
        f'The per-setting phase error combines SEM, speed reproducibility and σx = {r["sigma_x"]:.4f} rad. '
        f'Without σx, χ² = {r["chi2_without_sigma_x"]:.1f} for {r["dof"]} dof '
        f'(frame SEMs alone: {r["chi2_sem_only"]:.1f}). σx is set by χ² = dof. '
        'Calibration-curve covariance propagates as a correlated contribution to the slope.', '',
        'Steady time drift and wavelength uncertainty are outside the quoted budget. Settings were recorded '
        'in increasing order, with rest frames only at the start. The joint rotation/time fit has parameter '
        f'correlation {r["drift_check"]["correlation"]:.3f}, so it cannot separate a rotation slope from a steady drift. '
        'The rotation direction was not recorded.', '',
        '## Optional quadratic calibration check', '',
        f'The default line retains residual curvature: χ²/dof = {c["linear"]["chi2"]:.1f}/{c["linear"]["dof"]}, '
        f'p = {c["linear"]["p"]:.2g}. The quadratic has χ²/dof = '
        f'{c["fits_with_eps"]["2"]["chi2"]:.1f}/6 by construction because its residuals set ε. '
        'The line is the requested baseline approximation, not the statistically preferred calibration model. '
        'Residual model curvature is separate from the propagated coefficient covariance.', '',
        '| Downstream quantity | Linear baseline | Quadratic check |', '|---|---:|---:|',
        f'| a, rad/(rad/s) | {r["slope"]:.6f} | {optional["slope"]:.6f} |',
        f'| c, 10⁸ m/s | {sc["c"]/1e8:.6f} | {optional["speed_of_light"]["c"]/1e8:.6f} |',
        '', f'The quadratic changes the slope by {delta:+.2f}%. Its complete output is isolated in '
        '`optional_quadratic/`, and does not feed the main figures or numbers. Calibration-model sensitivity '
        'is reported separately from the quoted standard uncertainty.', '',
        '## Presentation assets', '',
        '`presentation_figures/` contains eight main PNGs and one optional model-comparison PNG. '
        'The source is `presentation_source/make_figures.py`. No animation or presentation file is delivered.', '',
        'Reproduce the complete pipeline using `README.md`. Refresh this document and `SUMMARY.md` with '
        '`python update_results.py` after both calibration runs.', '',
    ]
    (HERE/'RESULTS.md').write_text('\n'.join(common), encoding='utf-8')
    summary = [
        '# Sagnac effect — linear baseline summary', '',
        '11 October 2026 data. The linear calibration supplies all default angular velocities. '
        'Full tables, uncertainty details and interpretation are in `RESULTS.md` and `sagnac.txt`.', '',
        '## Setup', '',
        'Centred fibre loop, N = 13 ± 1 turns, R = 57.0 ± 0.5 cm, λ = 640 nm. '
        '161 photographs at 10 settings. Dark-ring radii give the centre phase through a fit in r².', '',
        '## Main results', '',
        f'- Linear k = {value(c["linear"]["k"], c["linear"]["u_k"])} rad/s per unit S, '
        f'ω₀ = {value(c["linear"]["omega0"], c["linear"]["u_omega0"])} rad/s.',
        f'- Phase slope a = {value(r["slope"], r["u_slope"])} (fit + calibration), '
        f'or {value(r["slope"], total_a)} including method spread, rad/(rad/s).',
        f'- Theory slope: {value(r["theory"], r["u_theory"])} rad/(rad/s), '
        f'{abs(z):.2f} combined σ from the measurement.',
        f'- Maximum shift: {value(maximum, u_max)} rad = '
        f'{value(maximum/(2*np.pi), u_max/(2*np.pi))} fringe.',
        f'- Inferred c = ({value(sc["c"]/1e8, sc["u_c"]/1e8)}) × 10⁸ m/s. '
        f'c/c₀ = {value(sc["ratio_to_defined"], sc["u_ratio_to_defined"])}.', '',
        '## Uncertainty and optional comparison', '',
        f'The total relative c uncertainty is {100*sc["u_c"]/sc["c"]:.1f}%, dominated by turn count '
        f'(7.7%). Phase extraction contributes {100*r["u_method"]/r["slope"]:.1f}%. '
        'Steady drift, unknown wavelength error and residual calibration-model curvature are outside the quoted budget.', '',
        f'The optional quadratic calibration changes a from {r["slope"]:.3f} to '
        f'{optional["slope"]:.3f} rad/(rad/s), a {delta:+.1f}% change. It improves the calibration '
        'residuals statistically but does not supply the default conversion. Its output is isolated in '
        '`optional_quadratic/`.', '',
        '## Presentation', '',
        'Eight main static images plus one optional calibration comparison, in `presentation_figures/`. '
        'No animation. Source: `presentation_source/make_figures.py`.', '',
    ]
    (HERE/'SUMMARY.md').write_text('\n'.join(summary), encoding='utf-8')


if __name__ == '__main__':
    main()
