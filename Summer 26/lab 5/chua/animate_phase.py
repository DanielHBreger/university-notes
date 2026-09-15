"""Create matched measurement/simulation phase-space GIFs.

Run: python chua/animate_phase.py --r 550 --window-ms 8
Uses the existing ideal model and a clean forward record nearest --r.
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

from gallery import nearest
from scope_data import read_scope
import simulate

HERE = Path(__file__).resolve().parent


def render(t, xy, title, detail, color, limits, output, seconds=12, fps=25):
    fig, ax = plt.subplots(figsize=(7.2, 6), dpi=100)
    fig.subplots_adjust(left=.13, right=.96, bottom=.16, top=.81)
    fig.suptitle(title, fontsize=18, y=.96)
    ax.set_title(detail, fontsize=10, pad=13)
    ax.set(xlabel='Voltage across C1, V1 (V)', ylabel='Voltage across C2, V2 (V)',
           xlim=limits[0], ylim=limits[1])
    ax.grid(alpha=.2)
    history, = ax.plot([], [], color=color, alpha=.28, lw=.8)
    trail, = ax.plot([], [], color=color, lw=1.6)
    point, = ax.plot([], [], 'o', color=color, ms=6)
    clock = fig.text(.13, .065, '', fontsize=11)
    fig.text(.96, .065, 'Faint: elapsed path  |  Bright: last 0.25 ms',
             ha='right', fontsize=9, color='#555555')
    frame_times = np.linspace(0, t[-1], round(seconds * fps))

    def update(frame):
        now = frame_times[frame]
        stop = min(np.searchsorted(t, now, side='right'), len(t))
        start = np.searchsorted(t, max(0, now - .00025))
        history.set_data(xy[:stop, 0], xy[:stop, 1])
        trail.set_data(xy[start:stop, 0], xy[start:stop, 1])
        point.set_data([xy[stop-1, 0]], [xy[stop-1, 1]])
        clock.set_text(f't = {now*1000:5.2f} ms')
        return history, trail, point, clock

    animation = FuncAnimation(fig, update, frames=len(frame_times), blit=False)
    animation.save(output, writer=PillowWriter(fps=fps))
    update(len(frame_times)//2)
    fig.savefig(output.with_suffix('.png'))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r', type=float, default=550)
    parser.add_argument('--window-ms', type=float, default=8)
    args = parser.parse_args()
    if args.window_ms <= 0:
        parser.error('--window-ms must be positive')
    name, resistance = nearest(str(HERE / 'forward'), args.r)
    mt, measured = read_scope(HERE / 'forward' / name)
    mt = mt - mt[0]
    if args.window_ms * .001 > mt[-1]:
        parser.error('requested window exceeds measured record')
    mask = mt <= args.window_ms * .001
    mt, measured = mt[mask], measured[mask]
    g = simulate.make_g()
    initial, _ = simulate.sweep_starts(g, [resistance], 0.)
    settle = .1
    st, states = simulate.integrate([resistance], None, settle + mt[-1],
                                   .5e-6, 0., g, t_skip=settle, y0=initial)
    st = st - st[0]
    simulated = states[:, :2, 0]
    assert np.isfinite(simulated).all()
    limits = []
    for axis in range(2):
        low = min(measured[:, axis].min(), simulated[:, axis].min())
        high = max(measured[:, axis].max(), simulated[:, axis].max())
        pad = max((high-low)*.08, .05)
        limits.append((low-pad, high+pad))
    out = HERE / 'animations'
    out.mkdir(exist_ok=True)
    render(mt, measured, 'Measurements · C1–C2 phase space',
           f'forward/{name}  |  Rpot = {resistance:.2f} Ω', '#167aab', limits,
           out / 'measurements_c1_c2.gif')
    render(st, simulated, 'Simulation · C1–C2 phase space',
           f'Ideal model  |  Rpot = {resistance:.2f} Ω  |  after 100 ms settling',
           '#cc6530', limits, out / 'simulation_c1_c2.gif')
    metadata = dict(measurement=f'forward/{name}', rpot_ohm=resistance,
                    window_ms=mt[-1]*1000, playback_seconds=12,
                    model='ideal', c1_nF=11.5, c2_nF=100, l_mH=18,
                    r0_ohm=simulate.R0, rL_ohm=0, simulation_dt_s=.5e-6,
                    settling_s=settle,
                    note='Time is relative to each displayed window; trajectories are not phase aligned.')
    (out / 'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
