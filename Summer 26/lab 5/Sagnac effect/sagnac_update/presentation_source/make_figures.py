"""Static Sagnac presentation images, using the Chaos presentation's visual style.

The saved linear calibration and analysis are read without modification.
"""
from pathlib import Path
import json
import sys
import zipfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from calibration import load_turns, per_setting, omega_of
from sagnac import round1

OUT = ROOT / 'presentation_figures'
SOURCE = Path(__file__).resolve().parent
ASSETS = SOURCE / 'assets'
INK, TEAL, ORANGE, PURPLE, GREY = '#243444', '#147da3', '#d87530', '#8257a6', '#586673'
plt.rcParams.update({
    'font.family': ['Arial', 'DejaVu Sans'], 'font.size': 18,
    'mathtext.fontset': 'dejavusans', 'text.usetex': False,
    'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.edgecolor': '#85919b', 'axes.linewidth': .8,
    'axes.grid': True, 'grid.color': '#dce2e7', 'grid.alpha': .45,
    'legend.frameon': False, 'legend.fontsize': 15,
    'lines.linewidth': 2.4, 'savefig.facecolor': 'white',
})


def canvas(title):
    fig = plt.figure(figsize=(16, 9), dpi=120, facecolor='white')
    fig.text(.075, .93, title, fontsize=27, weight='bold')
    return fig


def photograph(fig, file, rect, sample=None, rings=False, color=TEAL):
    """Camera photograph, optionally with quantitative dark-ring fit overlays."""
    ax = fig.add_axes(rect)
    ax.imshow(Image.open(ROOT / file).convert('RGB'))
    if sample is not None:
        cx, cy = sample.cx, sample.cy
        ax.set_xlim(cx-520, cx+520)
        ax.set_ylim(cy+520, cy-520)
        if rings:
            theta = np.linspace(0, 2*np.pi, 900)
            for k in (1, 2, 3):
                radius = sample[f'R{k}']
                if np.isfinite(radius):
                    ax.plot(cx+radius*np.cos(theta), cy+radius*np.sin(theta),
                            color=color, linewidth=1.6, alpha=.95)
                    ax.text(cx+.91*radius, cy-.415*radius, f'm = {k-1}',
                            color='white', fontsize=11, ha='center', va='center',
                            bbox={'facecolor':'black', 'edgecolor':'none', 'alpha':.6, 'pad':1.5})
            ax.plot(cx, cy, '+', color='white', markersize=10, markeredgewidth=1.2)
    ax.set_axis_off()
    return ax


def save(fig, name):
    filename = OUT / f'{name}.png'
    fig.savefig(filename, dpi=120)
    plt.close(fig)
    return filename


def report(x, u):
    return round1(x, u).replace('+-', '±')


def representative(row, frames):
    group = frames[frames.S == row.S]
    return group.iloc[np.argmin(np.abs(group.phi_rad.to_numpy() - group.phi_rad.median()))]


def title_image(speeds, frames, results, calibration, profiles):
    fig = canvas('')
    fig.text(.075, .64, 'Measuring rotation\nwith the Sagnac effect', fontsize=42,
             weight='bold', color='#000000', linespacing=1.23)
    photograph(fig, 'Measurements/0/20261011162149.png', [.565, .28, .395, .5])
    fig.text(.075, .155, 'Daniel Haim Breger & Hila Hemo', fontsize=24, color='#000000')
    fig.text(.075, .085, '11 October 2026', fontsize=17, color=GREY)
    return save(fig, '01_title')


def setup_image(speeds, frames, results, calibration, profiles):
    fig = canvas('The Sagnac interferometer')
    ax = fig.add_axes([.075, .205, .435, .655])
    ax.imshow(Image.open(ASSETS / 'optical_setup_manual.png'))
    ax.set_xlim(93, 492)
    ax.set_ylim(355, 0)
    ax.set_axis_off()
    # The two equations show how the optical phase follows from the travel-time difference.
    fig.text(.57, .755, 'Travel time difference', fontsize=18, color=GREY)
    fig.text(.57, .66, r'$\Delta t=\frac{4NA}{c^2}\,\omega$', fontsize=31)
    fig.text(.57, .50, 'Interference phase', fontsize=18, color=GREY)
    fig.text(.57, .39, r'$\Delta\phi=\frac{2\pi c}{\lambda}\,\Delta t$', fontsize=31)
    fig.text(.57, .275, r'$\Delta\phi=\frac{8\pi NA}{\lambda c}\,\omega$',
             fontsize=31, color=TEAL)
    for x, text in [(.11, 'N = 13 ± 1'), (.385, 'R = 57.0 ± 0.5 cm'), (.745, 'λ = 640 nm')]:
        fig.text(x, .105, text, fontsize=22)
    return save(fig, '02_interferometer_and_theory')


def calibration_image(speeds, frames, results, calibration, profiles):
    fig = canvas('Rotation calibration')
    settings, turns = load_turns()
    tab = per_setting(settings, turns)
    model = calibration['conversion']
    assert model['degree'] == 1
    u = np.hypot(tab.u_sw, model['eps'] * tab.omega)
    ax = fig.add_axes([.105, .145, .82, .68])
    xs = np.linspace(0, 94, 200)
    ax.errorbar(tab.S, tab.omega, u, fmt='o', color=TEAL, ms=8, capsize=4, zorder=3,
                label='Timed rotations')
    ax.plot(xs, omega_of(xs, model['coef']), color=ORANGE, label='Linear calibration')
    ax.set(xlim=(0, 95), ylim=(-.04, 2.03), xlabel='Motor setting S', ylabel='Angular velocity ω (rad/s)')
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_major_locator(MultipleLocator(.5))
    ax.legend(loc='lower right', bbox_to_anchor=(.98, .055))
    # A measured timing series motivates the calibration, rather than a prose explanation.
    j = 4
    ins = fig.add_axes([.175, .545, .235, .235], facecolor='white')
    lap = turns[j]
    ins.plot(np.arange(1, len(lap)+1), lap, 'o', color=TEAL, ms=5)
    ins.axhline(lap.mean(), color=ORANGE, linewidth=1.7)
    ins.set(xticks=np.arange(1, len(lap)+1), ylim=(6.84, 7.005),
            xlabel='Turn', ylabel='Time (s)')
    ins.tick_params(labelsize=12)
    ins.xaxis.label.set_size(13)
    ins.yaxis.label.set_size(13)
    ins.set_title(f'S = {settings[j]:.2f}', fontsize=15, pad=9)
    ins.text(.98, .85, f'⟨T⟩ = {lap.mean():.2f} s', transform=ins.transAxes,
             fontsize=13, ha='right', color=ORANGE)
    ax.text(.105, .41, r'$\omega=2\pi/\langle T\rangle$', transform=ax.transAxes, fontsize=23)
    linear = calibration['linear']
    ax.text(.49, .38, r'$\omega(S)=kS+\omega_0$', transform=ax.transAxes, fontsize=26, color=ORANGE)
    ax.text(.49, .285, f'k = {report(linear["k"], linear["u_k"])} rad/s',
            transform=ax.transAxes, fontsize=20)
    ax.text(.49, .21, f'ω₀ = {report(linear["omega0"], linear["u_omega0"])} rad/s',
            transform=ax.transAxes, fontsize=20)
    return save(fig, '03_linear_rotation_calibration')


def patterns_image(speeds, frames, results, calibration, profiles):
    fig = canvas('Rotation shifts the interference rings')
    # The large upper-centre rectangle is intentionally blank for the user's animation.
    # Three equally sized stills remain as a chronological reference along the bottom.
    for x, idx in zip([.075, .365, .655], [0, 4, 9]):
        row = speeds.iloc[idx]
        sample = representative(row, frames)
        photograph(fig, sample.file, [x, .09, .27, .27])
        fig.text(x+.135, .385, f'ω = {row.omega:.2f} rad/s', fontsize=19, ha='center')
    return save(fig, '04_measured_fringe_shift')


def extraction_image(speeds, frames, results, calibration, profiles):
    fig = canvas('The ring displacement gives the phase shift')
    samples = [representative(speeds.iloc[idx], frames) for idx in (0, 9)]
    for sample, color, y, text in zip(samples, [TEAL, ORANGE], [.515, .135], ['At rest', 'ω = 1.81 rad/s']):
        photograph(fig, sample.file, [.105, y, .19, .3378], sample=sample, rings=True, color=color)
        fig.text(.20, y+.352, text, fontsize=18, color=color, ha='center')
    ax = fig.add_axes([.455, .22, .465, .57])
    etas = []
    for sample, color, name in zip(samples, [TEAL, ORANGE], ['At rest', 'Rotating']):
        radii = np.array([sample[f'R{k}'] for k in (1, 2, 3)], float)
        squares = np.sort(radii[np.isfinite(radii)])**2
        order = np.r_[0, np.cumsum(np.rint(np.diff(squares)/results['P_global']))]
        period, s0 = np.polyfit(order, squares, 1)
        eta = s0/period
        etas.append(eta)
        xs = np.linspace(0, 2.1, 100)
        ax.plot(xs, xs+eta, color=color, label=name)
        ax.scatter(order, squares/period, s=80, color=color, zorder=3)
    ax.set(xlim=(-.2, 2.2), ylim=(-.08, 2.85), xlabel='Dark-ring order m', ylabel=r'$R_m^2/P$')
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(.5))
    ax.set_title(r'$R_m^2/P=m+\eta$', fontsize=25, pad=18)
    ax.legend(loc='upper left', fontsize=15)
    # The intercept displacement is the measured fraction of a fringe.
    ax.annotate('', xy=(.13, etas[0]), xytext=(.13, etas[1]),
                arrowprops={'arrowstyle':'|-|', 'color':INK, 'linewidth':1.5})
    frac = etas[0]-etas[1]
    ax.text(.24, .5*(etas[0]+etas[1]), f'{frac:.2f} fringe', fontsize=17,
            color=INK, va='center', bbox={'facecolor':'white', 'edgecolor':'none', 'pad':2})
    fig.text(.455, .095, r'$\Delta\phi=2\pi\times 0.52\;\simeq\;3.3\ \mathrm{rad}$',
             fontsize=27, color=TEAL)
    return save(fig, '05_phase_extraction')


def phase_image(speeds, frames, results, calibration, profiles):
    fig = canvas('Phase shift versus angular velocity')
    ax = fig.add_axes([.105, .145, .82, .685])
    xs = np.linspace(0, 1.89, 200)
    ath, uth = results['theory'], results['u_theory']
    ax.fill_between(xs, (ath-uth)*xs, (ath+uth)*xs, color=GREY, alpha=.13, linewidth=0)
    ax.plot(xs, ath*xs, '--', color=GREY, label='Sagnac prediction ±1σ')
    ax.plot(xs, results['slope']*xs+results['intercept']-speeds.phi.iloc[0],
            color=ORANGE, label='Linear fit')
    ax.errorbar(speeds.omega, speeds.dphi, speeds.u_fit,
                np.hypot(speeds.u_omega_curve, speeds.u_omega_repro),
                fmt='o', color=TEAL, capsize=4, ms=7, zorder=3, label='Measured phase')
    ax.set(xlim=(-.05, 1.99), ylim=(-.2, 3.94), xlabel='Angular velocity ω (rad/s)', ylabel='Phase shift Δφ (rad)')
    ax.xaxis.set_major_locator(MultipleLocator(.4))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.text(.045, .87, r'$\Delta\phi=a\omega+b$', transform=ax.transAxes, fontsize=27)
    ax.text(.045, .765, f'a = {report(results["slope"], results["u_slope"])} rad/(rad/s)',
            transform=ax.transAxes, fontsize=24, color=ORANGE)
    ax.text(.045, .625, r'$a_{\mathrm{th}}=\frac{8\pi NA}{\lambda c}$',
            transform=ax.transAxes, fontsize=23, color=GREY)
    ax.text(.045, .515, f'     = {ath:.2f} ± {uth:.2f} rad/(rad/s)',
            transform=ax.transAxes, fontsize=21, color=GREY)
    maximum = speeds.dphi.iloc[-1]
    umax = np.hypot(speeds.u_fit.iloc[-1], speeds.u_fit.iloc[0])
    ax.annotate(f'{report(maximum, umax)} rad', xy=(speeds.omega.iloc[-1], maximum),
                xytext=(1.47, 3.72), fontsize=20, color=TEAL,
                arrowprops={'arrowstyle':'-', 'color':TEAL, 'linewidth':1.2})
    ax.legend(loc='lower right', fontsize=15)
    return save(fig, '06_phase_vs_rotation')


def methods_image(speeds, frames, results, calibration, profiles):
    fig = canvas('Sensitivity to the phase-extraction method')
    columns = ['phi', 'phi_cos', 'phi_template', 'phi_points']
    labels = ['Ring minima', 'Cosine fit', 'Template', 'Random points']
    colors = [TEAL, ORANGE, PURPLE, GREY]
    entries = [results]+[results['methods'][k] for k in ['cosine fit', 'template matching', 'random points']]
    accepted = [True, True, True, False]        # random points: shown, not in the method range
    ax = fig.add_axes([.105, .255, .42, .535])
    main = speeds.phi-speeds.phi.iloc[0]
    ax.axhline(0, color=TEAL, linewidth=1.8, label='Ring minima')
    for col, name, color in zip(columns[1:], labels[1:], colors[1:]):
        relative = speeds[col]-speeds[col].iloc[0]
        ax.plot(speeds.omega, relative-main, 'o-', ms=5, color=color, label=name, linewidth=1.8)
    ax.set(xlim=(-.04, 1.88), ylim=(-.22, .33), xlabel='Angular velocity ω (rad/s)',
           ylabel='Phase difference from ring minima (rad)')
    ax.xaxis.set_major_locator(MultipleLocator(.4))
    ax.yaxis.set_major_locator(MultipleLocator(.1))
    ax.legend(loc='upper left', ncol=2, fontsize=13, handlelength=1.8, columnspacing=1.1)
    ax.set_title('Same photographs, different phase estimates', fontsize=18, pad=22)
    forest = fig.add_axes([.685, .255, .235, .535])
    ys = np.arange(4)[::-1]
    values = [entry['slope'] for entry, ok in zip(entries, accepted) if ok]
    forest.axvspan(min(values), max(values), color=TEAL, alpha=.07)
    for y, entry, color, ok in zip(ys, entries, colors, accepted):
        forest.errorbar(entry['slope'], y, xerr=entry['u_slope'], fmt='o', color=color,
                        mfc=color if ok else 'white', capsize=5, ms=8, linewidth=2.4)
        forest.text(entry['slope'], y+.28, report(entry['slope'], entry['u_slope']),
                    ha='center', color=color, fontsize=16)
    forest.set(yticks=ys, yticklabels=labels, ylim=(-.55, 3.7), xlim=(1.7, 1.94),
               xlabel='Slope a  [rad/(rad/s)]')
    forest.set_title('Fitted Sagnac slope', fontsize=18, pad=22)
    forest.spines['left'].set_visible(False)
    forest.tick_params(axis='y', length=0, pad=10, labelsize=15)
    forest.xaxis.set_major_locator(MultipleLocator(.1))
    forest.grid(axis='y', visible=False)
    fig.text(.685, .168, 'Open marker: ten fixed pixels, not in the method range', fontsize=13, color=GREY)
    fig.text(.105, .115, r'$u_{\mathrm{method}}=\frac{a_{\max}-a_{\min}}{2}$'
             +f' = {results["u_method"]:.3f} rad/(rad/s)', fontsize=24)
    total = np.hypot(results['u_slope'], results['u_method'])
    fig.text(.685, .115, f'a = {report(results["slope"], total)}', fontsize=25, color=TEAL, weight='bold')
    return save(fig, '07_phase_method_comparison')


def result_image(speeds, frames, results, calibration, profiles):
    fig = canvas('The speed of light from the Sagnac slope')
    c = results['speed_of_light']
    value = report(c['c']/1e8, c['u_c']/1e8).replace('±', r'\pm')
    fig.text(.105, .80, r'$c=\frac{8\pi NA}{\lambda a}=(' + value +
             r')\times 10^{8}\ \mathrm{m/s}$', fontsize=30, color=TEAL)
    # Each row adds an independent standard-uncertainty component in quadrature.
    keys = ['fit', 'calibration curve', 'method', 'R', 'N']
    labels = ['Phase fit', '+ rotation calibration', '+ phase method', '+ loop radius R', '+ turn count N']
    cumulative = np.sqrt(np.cumsum([c['components'][k]**2 for k in keys]))
    ax = fig.add_axes([.285, .16, .56, .50])
    ys = np.arange(5)[::-1]
    ax.axvline(299792458/1e8, color=GREY, linestyle='--', linewidth=1.6)
    for i, (y, u) in enumerate(zip(ys, cumulative)):
        ax.errorbar(c['c']/1e8, y, xerr=u/1e8, fmt='o', color=TEAL,
                    capsize=8, ms=10 if i==4 else 7, linewidth=4 if i==4 else 2.2,
                    markeredgecolor='white', markeredgewidth=.7, zorder=3)
        ax.text(3.265, y, f'{100*u/c["c"]:.1f}%', fontsize=19, va='center',
                weight='bold' if i==4 else 'normal', color=TEAL if i==4 else INK)
    ax.set(xlim=(2.51, 3.18), ylim=(-.55, 4.65), yticks=ys, yticklabels=labels,
           xlabel='Speed of light c  (10⁸ m/s)')
    ax.xaxis.set_major_locator(MultipleLocator(.1))
    ax.tick_params(axis='y', length=0, pad=14, labelsize=18)
    ax.spines['left'].set_visible(False)
    ax.grid(axis='y', visible=False)
    ax.text(299792458/1e8+.008, 4.38, 'Defined c₀', fontsize=17, color=GREY, va='bottom')
    fig.text(.105, .71, 'Uncertainty included', fontsize=17, color=GREY)
    fig.text(.915, .71, 'u(c) / c', fontsize=17, color=GREY, ha='center')
    return save(fig, '08_speed_of_light_and_uncertainty')


def optional_calibration_image(speeds, frames, results, calibration, profiles):
    fig = canvas('Optional: calibration-model sensitivity')
    settings, turns = load_turns()
    tab = per_setting(settings, turns)
    u = np.hypot(tab.u_sw, calibration['conversion']['eps']*tab.omega)
    ax = fig.add_axes([.105, .295, .82, .515])
    for degree, color, marker, offset, name in [(1, ORANGE, 'o', -.65, 'Linear baseline'),
                                               (2, TEAL, 's', .65, 'Quadratic check')]:
        model = calibration['fits_with_eps'][str(degree)]
        residual = 1e3*(tab.omega-omega_of(tab.S, model['coef']))
        ax.errorbar(tab.S+offset, residual, 1e3*u, fmt=marker, color=color, capsize=4,
                    ms=7, label=f'{name}  (χ²/dof = {model["chi2"]:.1f}/{model["dof"]})')
    ax.axhline(0, color=GREY, linewidth=1)
    ax.set(xlim=(0, 95), xlabel='Motor setting S', ylabel='Calibration residual (10⁻³ rad/s)')
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.legend(loc='lower left', fontsize=16)
    alternative = json.loads((ROOT/'optional_quadratic/sagnac.json').read_text())
    for x, text, color in [(.105, 'Linear baseline', ORANGE), (.505, 'Quadratic check', TEAL)]:
        fig.text(x, .155, text, fontsize=18, color=color, weight='bold')
    fig.text(.105, .09, f'a = {results["slope"]:.3f}     c = {results["speed_of_light"]["c"]/1e8:.3f} × 10⁸ m/s', fontsize=20)
    fig.text(.505, .09, f'a = {alternative["slope"]:.3f}     c = {alternative["speed_of_light"]["c"]/1e8:.3f} × 10⁸ m/s', fontsize=20)
    return save(fig, 'optional_A1_quadratic_calibration')


FIGURES = [title_image, setup_image, calibration_image, patterns_image, extraction_image,
           phase_image, methods_image, result_image, optional_calibration_image]


def main():
    OUT.mkdir(exist_ok=True)
    speeds = pd.read_csv(ROOT/'sagnac_speeds.csv')
    frames = pd.read_csv(ROOT/'fringes_frames.csv')
    results = json.loads((ROOT/'sagnac.json').read_text())
    calibration = json.loads((ROOT/'calibration.json').read_text())
    assert results['calibration_degree'] == calibration['conversion']['degree'] == 1
    z = np.load(ROOT/'fringes_profiles.npz')
    profiles = z['s'], z['profiles'], list(z['files'])
    names = [function(speeds, frames, results, calibration, profiles) for function in FIGURES]
    metadata = dict(main_images=8, optional_images=1, image_size=[1920,1080],
                    calibration='linear', calibration_degree=1, animation=False,
                    user_animation_space={'image':'04_measured_fringe_shift.png',
                                          'rectangle_px':[442,119,1037,435],
                                          'background':'white'},
                    files=[p.name for p in names], sources=['calibration.json','sagnac.json',
                    'sagnac_speeds.csv','fringes_frames.csv','fringes_profiles.npz'])
    (OUT/'manifest.json').write_text(json.dumps(metadata, indent=2)+'\n')
    (OUT/'README.md').write_text(
        '# Sagnac slide images\n\n'
        'Eight main PNG images in numbered order, plus one optional calibration comparison. '
        'Each is 1920 × 1080 (16:9). Paste each image onto a blank slide. '
        'Explanations and transitions are in presenter_notes.md. No presentation file is included.\n\n'
        'The sequence follows the measurement: apparatus and phase prediction, timed rotation calibration, '
        'observed ring motion, ring displacement to phase, phase versus rotation, method sensitivity, '
        'and speed of light with the uncertainty build-up.\n\n'
        'Image 04 has a white central area for an animation you add yourself. '
        'Use a rectangle beginning at x = 442 px, y = 119 px, width 1037 px, height 435 px '
        '(top-left image coordinates), or fit a smaller animation inside it. '
        'Three original stills below it show the same progression. No animation is supplied.\n\n'
        'Every main result uses the linear calibration. The quadratic calibration appears only '
        'in optional_A1. The phase fit image quotes fit and calibration uncertainty. '
        'Image 07 adds the phase-method spread. Image 08 then adds the loop geometry. '
        'All displayed uncertainties are standard uncertainties.\n\n'
        'Sources and the single consolidated discussion of limitations are in presenter_notes.md. '
        'Source data and rebuild code remain in the experiment folder.\n', encoding='utf-8')
    notes = SOURCE/'presenter_notes.md'
    (OUT/'presenter_notes.md').write_text(notes.read_text(), encoding='utf-8')
    package = ROOT/'Sagnac slide images.zip'
    with zipfile.ZipFile(package, 'w', zipfile.ZIP_DEFLATED) as archive:
        for p in names+[OUT/'README.md', OUT/'presenter_notes.md', OUT/'manifest.json']:
            archive.write(p, p.name)
    print(f'Created {len(names)} PNG images and {package}')


if __name__ == '__main__':
    main()
