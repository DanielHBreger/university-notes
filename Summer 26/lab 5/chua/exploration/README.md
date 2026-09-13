# Exploration (superseded)

Scratch scripts and their outputs from the 11 September 2026 session that searched for the
cause of the model/measurement mismatch: per-record fits of C1 and of the tank from the
records (`audit_dynamics.py`, `calibrate_model.py`, `fit_trajectories.py`), spline and
rounded-knee diodes (`fit_curve.py`, `fit_rounded.py`, `fit_curvature.py`,
`conditional_curve.py`), a first-order diode response (`fit_response.py`,
`response_integration.py`, `check_response.py`), an amplitude-dependent inductor
(`fit_inductor.py`), saturation and onset scans, and a revalidation of every record
(`audit_measurements.py`). `measured_delta.py` is an earlier Feigenbaum estimate that
referred to data sets that are no longer in the tree.

Everything these found is now done properly by `../identify.py` (cycle-averaged
identification of C1, C2, the Rayleigh inductor and the two-op-amp diode) and used by
`../simulate.py --model bench`. Nothing here is run by the pipeline; the scripts are kept for
the record and may not run against the current `simulate.py`. The `.npz` caches are not
tracked.
