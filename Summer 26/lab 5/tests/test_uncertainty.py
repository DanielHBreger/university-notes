"""Independent propagation identities and correlated resampling checks."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'chua'))
from uncertainty import (scope_vertical_limit,timebase_limit_fraction,supply_voltage_limit,
 current_covariance,divider_standard_uncertainty,delta_covariance,divider_block_bootstrap,voltage_quantum)


def test_full_scale_not_reading_and_probe_conversion():
    assert scope_vertical_limit(1.,0.,1.)==pytest.approx(.282)
    assert scope_vertical_limit(10.,0.,10.)==pytest.approx(2.82)
    assert scope_vertical_limit(None,0.,1.) is None
    assert scope_vertical_limit(.001,0.,1.)==pytest.approx(.00312)


def test_supply_accuracy_and_unknown_clock_age():
    assert supply_voltage_limit(9)==pytest.approx(.0127)
    assert timebase_limit_fraction(None) is None
    assert timebase_limit_fraction(2)==pytest.approx(35e-6)


def test_current_covariance_monte_carlo():
    rng=np.random.default_rng(17); vv=rng.normal(0,.1,300000); vi=rng.normal(0,.2,300000)
    empirical=np.cov(vi,(vv-vi)/216)
    expected=current_covariance(.1,.2,216)
    np.testing.assert_allclose(empirical,expected,rtol=.015)
    assert expected[0,1]<0


def test_common_voltage_error_cancels_current():
    result=current_covariance(.2,.2,216,cov_vv_vi=.04)
    assert result[1,1]==pytest.approx(0,abs=1e-18)


def test_divider_covariance_and_missing_calibration():
    # Perfectly common fractional errors in A,B cancel in B/A.
    ab=np.array([.6,.4]);cov=np.outer(ab,ab)*.01**2
    fit,total=divider_standard_uncertainty(*ab,cov)
    assert fit<1e-6 and total is None
    fit,total=divider_standard_uncertainty(*ab,np.zeros((2,2)),u_r0=2.,u_g21=.01)
    assert total==pytest.approx(np.hypot(2*.4/.6,992*.4/.6*.01))


def test_delta_shared_middle_transition_and_common_scale():
    r=np.array([775.,756.,750.]); cov=np.diag([1.5**2,.75**2,.4**2])
    value,u=delta_covariance(r,cov)
    jac=[]
    for k in range(3):
        eps=np.zeros(3);eps[k]=1e-4
        f=lambda x:(x[0]-x[1])/(x[1]-x[2])
        jac.append((f(r+eps)-f(r-eps))/2e-4)
    assert u==pytest.approx(np.sqrt(np.array(jac)@cov@jac),rel=1e-7)
    _,common=delta_covariance(r,.01**2*np.outer(r,r))
    assert common<1e-6


def test_block_bootstrap_exact_divider_and_drift():
    t=np.arange(40000)*1e-5; a=np.sin(700*t);b=np.cos(923*t)
    vm=.6*a+.4*b+.1
    exact=divider_block_bootstrap(np.column_stack([a,b,vm]),1e-5,draws=100)
    assert exact['u_r_ohm']<1e-7
    vm+=.01*np.sin(13*t)*a
    drift=divider_block_bootstrap(np.column_stack([a,b,vm]),1e-5,draws=100)
    assert drift['u_r_ohm']>.05


def test_code_spacing_handles_csv_rounding():
    x=np.round(np.arange(200)*.0241206,6)
    q,quality=voltage_quantum(x)
    assert q==pytest.approx(.0241206,abs=1e-6)
    assert quality==1
