"""Report-form rounding: one significant figure, value to the same digit, half up."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sagnac import round1  # noqa: E402


def test_round1():
    assert round1(1.8713, 0.0275) == '1.87 +- 0.03'
    assert round1(3.79, 0.97) == '4 +- 1'          # the course convention, also for a leading 1
    assert round1(2.25, 0.1) == '2.3 +- 0.1'       # a tie goes up, not to the even digit
    assert round1(1.005, 0.002) == '1.005 +- 0.002'
    assert round1(1.0049, 0.005) == '1.005 +- 0.005'
    assert round1(2.345, 0.096) == '2.3 +- 0.1'    # 0.096 rounds to 0.1 and moves the digit
    assert round1(123.0, 14.0) == '120 +- 10'
    assert round1(-0.004, 0.03) == '0.00 +- 0.03'  # no negative zero
