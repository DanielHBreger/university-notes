"""
Measurement uncertainty components. A missing calibration is reported as
unknown (None), never as zero.

Sources: scope specification Keysight 5990-6619EN (2014), pp. 18-19, 22;
power supply ISO-TECH IPS X303 manual, pp. 49-50. Specification limits and
standard uncertainties are kept separate throughout.
"""
import numpy as np


def voltage_quantum(v):
    """
    (quantum, consistency): the typical gap between adjacent occupied codes of
    a channel, robust to CSV rounding, and the fraction of gaps that are a
    whole multiple of it.

    Does not infer V/div, ADC span, probe factor or analogue noise, and assumes
    the adjacent codes are well occupied; check the consistency.
    """
    codes = np.unique(np.asarray(v, float))
    gaps = np.diff(codes)
    gaps = gaps[gaps > 0]
    if len(gaps) < 3:
        return np.nan, np.nan
    q = float(np.median(gaps))
    consistency = float(np.mean(np.abs(gaps / q - np.rint(gaps / q)) < .002))
    return q, consistency


def scope_vertical_limit(scale, offset, probe=1.):
    """
    Typical accuracy envelope of a single DC reading, in probe-tip volts.

    `scale` and `offset` are tip-referred; the 2 mV instrument-input term is
    multiplied by the attenuation. Magnified 1-2 mV/div ranges use the 4 mV/div
    input-equivalent scale. Full scale is 8 divisions; gain 2 % of full scale
    plus 0.25 % for the cursor. This is a DC/cursor envelope, not the standard
    deviation of waveform noise, and excludes probe accuracy and drift.
    """
    if any(x is None for x in (scale, offset, probe)):
        return None
    if not np.isfinite([scale, offset, probe]).all() or scale <= 0 or probe <= 0:
        raise ValueError('invalid scope settings')
    effective_scale = max(scale, .004 * probe)
    return .0225 * 8 * effective_scale + .1 * effective_scale + .002 * probe + .01 * abs(offset)


def timebase_limit_fraction(years):
    """Timebase limit: 25 ppm plus 5 ppm per year of ageing; None if the age is unknown."""
    if years is None:
        return None
    if not np.isfinite(years) or years < 0:
        raise ValueError('invalid calibration age')
    return (25 + 5 * years) * 1e-6


def supply_voltage_limit(v):
    """Programming OR readback limit of a CH1/CH2 setting (not their sum), in volts."""
    if not np.isfinite(v) or abs(v) > 30:
        raise ValueError('outside 0-30 V channel range')
    return .0003 * abs(v) + .010


def current_covariance(u_vv, u_vi, rs, current=0., u_rs=0., cov_vv_vi=0.):
    """
    Covariance matrix of (VI, I) with I = (VV - VI)/Rs, for an Rs independent
    of the voltages. The VI term shared by both is kept rather than treating
    the voltage and current errors as independent. Inputs are standard
    uncertainties, not limits.
    """
    if min(u_vv, u_vi, u_rs) < 0 or rs <= 0:
        raise ValueError('invalid uncertainties')
    var_i = (u_vv * u_vv + u_vi * u_vi - 2 * cov_vv_vi + (current * u_rs) ** 2) / rs ** 2
    cov = (cov_vv_vi - u_vi * u_vi) / rs
    return np.array([[u_vi * u_vi, cov], [cov, var_i]])


def divider_standard_uncertainty(A, B, cov_ab, r0=992., g21=1., u_r0=None, u_g21=None):
    """
    (u from the fit alone, combined u) of Rpot = r0 g21 B/A, propagating the
    covariance of A and B; the calibration terms are independent of the fit.
    The combined u is None while either calibration term is unknown. u_g21 is
    the absolute standard uncertainty of g2/g1 (usually near 1).
    """
    if A <= 0 or B < 0 or r0 <= 0 or g21 <= 0:
        raise ValueError('invalid divider')
    grad = np.array([-r0 * g21 * B / A ** 2, r0 * g21 / A])
    var = float(grad @ np.asarray(cov_ab) @ grad)
    u_fit = np.sqrt(max(0, var))
    if u_r0 is None or u_g21 is None:
        return float(u_fit), None
    total = np.sqrt(var + (g21 * B / A * u_r0) ** 2 + (r0 * B / A * u_g21) ** 2)
    return float(u_fit), float(total)


def delta_covariance(r, cov):
    """(delta, u) of delta = (R1 - R2)/(R2 - R3), including the covariance R2 shares with both gaps."""
    r = np.asarray(r, float)
    gap = r[1] - r[2]
    if gap <= 0:
        raise ValueError('R2 must exceed R3')
    delta = (r[0] - r[1]) / gap
    grad = np.array([1, -(1 + delta), delta]) / gap
    return float(delta), float(np.sqrt(max(0, grad @ np.asarray(cov) @ grad)))


def divider_block_bootstrap(d, dt, block_s=.005, draws=400, seed=0, r0=992.):
    """
    Paired, non-overlapping time-block bootstrap of the divider fit
    V3 = A V1 + B V2 + C on the channels d = (V1, V2, V3).

    Simultaneous blocks of all three channels move together. This measures how
    stable the fit is within a record, given the observed voltage scale; it
    does not identify absolute calibration, errors-in-variables bias or drift
    between records. A tail shorter than one block is left out.
    """
    d = np.asarray(d, float)
    length = max(10, round(block_s / dt))
    n = len(d) // length
    if n < 8:
        raise ValueError('need at least 8 time blocks')
    z = d[:n * length, :2] - d[:n * length, :2].mean(0)
    scale = z.std(0)
    if np.any(scale < 1e-12):
        raise ValueError('flat channel')
    x = np.column_stack((z / scale, np.ones(len(z))))
    y = d[:n * length, 2] - d[:n * length, 2].mean()

    # normal equations per block, then one weighted sum per bootstrap draw
    xb = x.reshape(n, length, 3)
    yb = y.reshape(n, length)
    xx = np.einsum('nki,nkj->nij', xb, xb)
    xy = np.einsum('nki,nk->ni', xb, yb)
    rng = np.random.default_rng(seed)
    weights = rng.multinomial(n, np.full(n, 1 / n), size=draws)
    coef = np.linalg.solve(np.einsum('bn,nij->bij', weights, xx),
                           np.einsum('bn,ni->bi', weights, xy)[..., None])[..., 0]

    ab = coef[:, :2] / scale
    valid = (ab[:, 0] > 0) & (ab[:, 1] >= 0) & np.isfinite(ab).all(1)
    ab = ab[valid]
    if len(ab) < .95 * draws:
        raise ValueError('bootstrap divider not identifiable')
    r = r0 * ab[:, 1] / ab[:, 0]
    return dict(u_r_ohm=float(np.std(r, ddof=1)), ci95_ohm=np.quantile(r, [.025, .975]).tolist(),
                cov_ab=np.cov(ab.T).tolist(), ab_mean=ab.mean(0).tolist(), n_blocks=n,
                valid_draws=len(r), block_s=block_s)
