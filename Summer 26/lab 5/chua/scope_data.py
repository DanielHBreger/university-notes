"""Strict readers for the scope's time/channel CSV format.

Removing invalid rows would join unrelated samples and corrupt return times.
Reject them instead. Keep timestamps in float64 even when voltages are coarse.
"""
import numpy as np
import pandas as pd

R0 = 990.0  # ohm, measured value in the Summer 2026 experiment plan
RPOT_MAX = 1000.0


def read_scope(path, channels=('CH1', 'CH2'), min_samples=3):
    head = pd.read_csv(path, nrows=0)
    names = [c.split('(')[0].strip() for c in head.columns]
    if not names or names[0].lower() != 'time':
        raise ValueError('first column must be Time(s)')
    if len(set(names)) != len(names) or len(set(channels)) != len(channels):
        raise ValueError('duplicate channel names')
    if any(c not in names or c == names[0] for c in channels):
        raise ValueError(f'required channels {channels}; found {names}')
    cols = [head.columns[0]] + [head.columns[names.index(c)] for c in channels]
    data = pd.read_csv(path, usecols=cols, dtype=np.float64)[cols].to_numpy()
    if len(data) < min_samples:
        raise ValueError(f'too few samples (need {min_samples})')
    if not np.isfinite(data).all():
        raise ValueError('non-finite sample; cannot bridge gaps in a time series')
    validate_time(data[:, 0])
    return data[:, 0], data[:, 1:]


def validate_time(t):
    steps = np.diff(t)
    if not len(steps) or not np.isfinite(t).all() or np.any(steps <= 0):
        raise ValueError('timestamps must be finite and strictly increasing')
    dt = float(np.median(steps))
    # %.6e scope CSV timestamps have finite decimal precision.
    tolerance = max(0.005 * dt, 2e-6 * float(np.max(np.abs(t))), 1e-15)
    if np.max(np.abs(steps - dt)) > tolerance:
        raise ValueError('timestamps are not uniformly sampled')
    return dt
