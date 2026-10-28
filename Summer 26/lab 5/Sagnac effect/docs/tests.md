# `tests/` and `sagnac_update/tests/`

Five tests in four files at the top level, all passing on 24 October 2026 with the report pipeline (Python 3.13, pytest 9.1). They check template matching and the fixed-pixel fit on synthetic data with a known phase sweep, the report-form rounding, and Eq. 5 solved for c. The statistics behind the final numbers (`linregress`, `setting_means`, the three-part uncertainty) have no test of their own; `sagnac.txt` prints every intermediate value instead.

## Running them

```bash
# from the 'Sagnac effect' folder: 5 tests
python -m pytest -q -p no:cacheprovider tests
```

`sagnac_update/tests/test_methods.py` (4 tests, presentation stage) imports `sagnac_update/methods.py` and `sagnac_update/sagnac.py`, which need the earlier `calibration.py` (`previous_reports/calibration_eps_2026-10-13.py`) one level up; with the report pipeline's `calibration.py` their import fails. They are not part of the report's checks.

## `tests/test_methods.py` (`methods.py`)

`make_profiles(phases)` builds ring profiles in s = r² on the 130 bins of `fringes.py` (`np.arange(1000, 260000, 2000.0)` px²): 80 e^(−s/2.6×10⁵) + 85 e^(−s/1.4×10⁵) cos(2πs/P + φ) + 7, with P = 107500 px² and one phase per frame.

- `test_template_recovers_shift`: 10 rest frames and 11 frames from 0.3 to 3.3 rad, noise σ = 1.5 counts. `methods.template_phases` finds the template's own period within 2 % of P and every phase within 0.1 rad.

## `tests/test_pixels.py` (`pixels.py`)

- `test_fit_recovers_phases`: ten pixels at random fringe phases with random offsets (60–120 counts) and amplitudes (20–50 counts), 10 rest frames and 11 frames from 0.3 to 3.3 rad (a sweep of half a fringe, as in the run), noise σ = 1 count. `pixels.fit_phases` converges and, after `relative_to_rest`, every phase is within 0.05 rad of the truth and every offset and amplitude within 2 counts.

## `tests/test_rounding.py`

`round1(x, u)` of `sagnac.py`, the report form of `sagnac.txt`: u to one significant figure, x to the same digit, half up on the printed digits.

| call | result | what it pins down |
|---|---|---|
| `round1(1.8713, 0.0275)` | `1.87 +- 0.03` | the basic rule |
| `round1(3.79, 0.97)` | `4 +- 1` | a leading 1 gets no second digit (the course convention) |
| `round1(2.25, 0.1)` | `2.3 +- 0.1` | a tie goes up, not to the even digit |
| `round1(1.005, 0.002)`, `round1(1.0049, 0.005)` | `1.005 +- 0.002`, `1.005 +- 0.005` | no binary rounding surprise |
| `round1(2.345, 0.096)` | `2.3 +- 0.1` | 0.096 becomes 0.1, so the value keeps one digit fewer |
| `round1(123.0, 14.0)` | `120 +- 10` | rounding left of the decimal point |
| `round1(-0.004, 0.03)` | `0.00 +- 0.03` | no negative zero |

## `tests/test_speed_of_light.py`

- `test_theory_slope_gives_defined_c`: the predicted slope a_th, put into `speed_of_light`, returns the defined c to a relative 10⁻¹². A slope part of 1 % gives a c part of 1 % (c ∝ 1/a), the N part is c u_N/N, and the R part is 2c u_R/R (A = πR²).
- `test_c_inverse_in_slope`: c for a = 1 is twice c for a = 2.

## Not covered

- `rms_from_line`, `setting_means`, the drift fit and the bend check in `sagnac.py`; the calibration fit; the fixed-pixel fit on real images (only the synthetic check above).
- `fringes.py` (ring finder, ring-minima phase, cosine fit). The method that gives the result is not checked against a known phase.
