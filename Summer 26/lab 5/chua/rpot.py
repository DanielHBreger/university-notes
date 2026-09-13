"""
Measure the potentiometer from the midpoint voltage, without opening the circuit.

R0 and the potentiometer are in series between the v_C1 and v_C2 nodes, with
CH3 at the junction. No current flows into the scope probe, so the midpoint is
a plain voltage divider:

    v_m = A * v_C1 + B * v_C2 + C

Fit A and B over the whole record, then

    Rpot = R0 * (B / A) * (g2 / g1)

This form uses both coefficients, so the midpoint channel's gain cancels out
entirely; only the relative gain of the two end channels survives, and that is
close to 1 when they sit on the same V/div. The fitted offset C absorbs any
DC offset in the channels.

Defaults use the retained sweep calibration: R0 = 992 ohm, CH1 = v_C1,
CH2 = v_C2, CH3 = midpoint.

Usage:
    python rpot.py record.csv
    python rpot.py record.csv --g21 0.9760
    python rpot.py *.csv                              # whole sweep at once
"""
import argparse
import glob
import numpy as np
from scope_data import R0, RPOT_MAX, read_scope
from sweeplib import CLEAN_DIVIDER_MAX_PCT


def rpot(path, r0, g21, c1='CH1', c2='CH2', mid='CH3'):
    """Return (Rpot, residual_fraction) for one record."""
    if not np.isfinite([r0, g21]).all() or r0 <= 0 or g21 <= 0:
        raise ValueError('r0 and g21 must be positive and finite')
    _, d = read_scope(path, (c1, c2, mid), min_samples=10)
    return fit_divider(d, r0, g21)


def fit_divider(d, r0=R0, g21=1.0):
    """Fit already validated simultaneous CH1, CH2, CH3 samples."""
    if not np.isfinite([r0,g21]).all() or min(r0,g21)<=0:
        raise ValueError('r0 and g21 must be positive and finite')
    if d.ndim != 2 or d.shape[1] != 3 or len(d) < 10 or not np.isfinite(d).all():
        raise ValueError('need at least ten finite three-channel samples')
    v1, v2, vm = d.T
    M = np.column_stack([v1 - v1.mean(), v2 - v2.mean()])
    scales = M.std(axis=0)
    if np.any(scales <= 1e-12) or np.ptp(vm) <= 1e-12:
        raise ValueError('flat channels cannot identify the divider')
    Z = M / scales
    coef, _, rank, singular = np.linalg.lstsq(Z, vm - vm.mean(), rcond=None)
    if rank < 2 or singular[-1] / singular[0] < 1e-4:
        raise ValueError('end channels are collinear; divider ratio is unidentifiable')
    A, B = coef / scales
    if A <= 1e-10 or B < 0:
        raise ValueError('non-physical divider coefficients; check channels/clipping')
    resid = float(np.sqrt(np.mean((vm - vm.mean() - Z @ coef) ** 2)))
    span = float(vm.max() - vm.min())
    return r0 * (B / A) * g21, resid / max(span, 1e-12)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('csv', nargs='+')
    p.add_argument('--r0', type=float, default=R0, help='fixed resistor, ohm')
    p.add_argument('--g21', type=float, default=1.0,
                   help='gain ratio of the v_C2 channel to the v_C1 channel')
    p.add_argument('--c1', default='CH1', help='column with v_C1')
    p.add_argument('--c2', default='CH2', help='column with v_C2')
    p.add_argument('--mid', default='CH3', help='column with the midpoint')
    a = p.parse_args()

    files = []
    for pat in a.csv:
        files.extend(sorted(glob.glob(pat)) or [pat])

    for f in files:
        try:
            r, res = rpot(f, a.r0, a.g21, a.c1, a.c2, a.mid)
        except Exception as e:
            print(f"{f}: {e}")
            continue
        warn = ("   CHECK: not a clean divider"
                if 100 * res > CLEAN_DIVIDER_MAX_PCT else "")
        if not 0 <= r <= RPOT_MAX:
            warn += "   CHECK: outside 0-1000 ohm potentiometer range"
        print(f"{f}: Rpot = {r:7.1f} ohm   (fit residual {100*res:.2f} %){warn}")


if __name__ == '__main__':
    main()
